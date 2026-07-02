#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder
E-Mail: manuel.seefelder@uni-ulm.de
Postal address: Department of Gene Therapy, University of Ulm, Helmholzstr. 8/1, 89081 Ulm, Germany

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=#

# spike/npe/benchmark.jl --- Phase-4 NPE-vs-ADVI benchmark (NPE-02 / NPE-03; D-04/D-08/D-09).
#
# The headline of the phase: the trained amortized NPE is scored against the
# per-dataset ADVI baseline in ρ-SPACE for accuracy (per-parameter RMSE + 90%
# interval width, plus Δρ RMSE) and for wall-clock SPEEDUP (t_advi / t_npe > 100),
# ALWAYS reported paired (speedup is never emitted alone; D-08).
#
# ρ-SPACE ALIGNMENT (RESEARCH Pattern 2 / A1): the NPE infers the simulator θ
# (ρ_true is row 1, the KNOWN coloc ground truth); the Turing/ADVI baseline infers
# the induced per-condition mean μ. They live in different spaces, so every
# comparison is made in ρ-space: the NPE ρ̂ is θ row 1 (un-standardized), the ADVI
# ρ̂ is `ghat(μ)` per draw (precomputed in the artifact as `rho_sample`/`rho_control`).
# Both score against each stack's known ρ_true = holdout `theta[1, ·]`.
#
# STRICT global_index JOIN (D-04): the ADVI artifact keys every pair by holdout
# `global_index`; NPE↔ADVI rows are matched by that index (never positionally), so a
# re-ordered holdout cannot silently misalign the two methods.
#
# NPE RMSE VIA assess/rmse (RESEARCH "Don't Hand-Roll"): per-parameter RMSE is taken
# from NeuralEstimators `rmse(assess(est, θ_std, Z))` -- the estimator draws are in
# STANDARDIZED θ-space (train_npe.jl fits `θzt` on θtr only), so θ_test is
# standardized with the SAME frozen `θzt` and the resulting standardized RMSE is
# scaled back to physical units by `θzt.scale` (infer.jl header contract, Pitfall 5).
#
# CPU-ONLY (D-10 / Pitfall 1): every NeuralEstimators call passes `use_gpu=false`.
# GPU timing is a separate opt-in branch guarded by a CUDA-device check (never gating).
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; reaches src/ only transitively
# through the guarded infer.jl / resimulate.jl chains. The ADVI artifact is the ONLY
# cross-env hand-off (D-02) -- this file READS it, never runs Turing.

using JLD2           # read the cross-env ADVI artifact + the trained NPE
using Statistics     # mean, median
using Random         # reproducible NPE posterior sampling
using StatsBase      # transform / reconstruct (frozen θ-standardization)
using NeuralEstimators   # assess, rmse (qualified at call site; name collides with Images)
using BenchmarkTools     # @belapsed -- warmup+median wall-clock (required for >100×, D-08)

# --- ORDER MATTERS: infer.jl (posterior surface + load_npe/load_holdout + ghat via
#     its own chain), then resimulate.jl (raw-stack re-sim + patch_summary/encode_d01
#     for the D-08 NPE clock). Guarded for idempotency (standalone + runtests.jl).
isdefined(@__MODULE__, :posterior_for)      || include(joinpath(@__DIR__, "infer.jl"))
isdefined(@__MODULE__, :resimulate_holdout) || include(joinpath(@__DIR__, "resimulate.jl"))

# Committed defaults: the persisted NPE (04-03) and the cross-env ADVI artifact (04-04).
const DEFAULT_NPE_MODEL     = joinpath(@__DIR__, "trained_npe.jld2")
const DEFAULT_ADVI_ARTIFACT = joinpath(@__DIR__, "..", "baseline", "advi_artifact.jld2")

"""
    _rmse_vector(rmse_df) -> Vector{Float64}

Extract the per-parameter RMSE from a NeuralEstimators `rmse(assessment)` DataFrame
into a 7-vector in θ field order (ρ_true first). Keyed by the default parameter
names `θ1..θ7` so a re-ordered DataFrame cannot scramble the parameter axis. Uses
plain `getproperty`/`getindex` on the DataFrame (no `using DataFrames`, a transitive
dep) so the spike env's direct-dep set is unchanged.
"""
function _rmse_vector(rmse_df)
    names  = rmse_df.parameter
    values = rmse_df.rmse
    lut    = Dict(String(names[i]) => Float64(values[i]) for i in 1:length(values))
    return [lut["θ$i"] for i in 1:7]
end

"""
    _theta_scale(θzt) -> Vector{Float64}

The per-parameter scale of the frozen θ `ZScoreTransform`, used to map a
standardized-space RMSE back to physical units (RMSE_phys = RMSE_std · scale).
"""
_theta_scale(θzt) = collect(Float64.(θzt.scale))

"""
    _load_artifact(path) -> Dict

Read the cross-env ADVI artifact (D-02), asserting the keys the benchmark consumes
are present (schema-version integrity, mirroring the cache.jl reopen-check idiom).
"""
function _load_artifact(path)
    isfile(path) || error("benchmark: no advi_artifact.jld2 at $path")
    d = JLD2.load(path)
    for k in ("pairs", "rho_sample", "rho_control", "rho_s_interval",
              "rho_c_interval", "wall_clock")
        haskey(d, k) || error("benchmark: advi_artifact.jld2 missing key '$k'")
    end
    return d
end

"""
    _holdout_summary(ho, variant::Symbol) -> AbstractMatrix

Select the RAW holdout summary matrix that matches the trained model's `variant`
(WR-01): `:min` → the 128-dim `summary_min`, `:aug` → the AUG_DIM `summary_aug`.
Threading `variant` through here (instead of hardcoding `summary_min`) keeps the
benchmark correct for an `:aug` model -- a valid `choose_summary` ablation outcome --
whose summary net expects AUG_DIM inputs, not 128.
"""
_holdout_summary(ho, variant::Symbol) =
    variant === :aug ? ho.summary_aug :
    variant === :min ? ho.summary_min :
    error("_holdout_summary: unknown variant $variant (use :min or :aug)")

"""
    _encode_variant(mci, variant::Symbol) -> Vector{Float64}

Extract the RAW summary vector for one raw stack in the encoding that matches the
model `variant` (WR-01): `:min` → the 128-dim `encode_d01(patch_summary(mci))`,
`:aug` → the AUG_DIM `encode_aug(mci, patch_summary(mci))`. `patch_summary` is
computed once and shared, so the timed on-the-fly NPE clock measures the full
variant-correct summary extraction (never the wrong 128-dim `:min` vector for an
`:aug` net).
"""
_encode_variant(mci, variant::Symbol) =
    variant === :aug ? encode_aug(mci, patch_summary(mci)) :
    variant === :min ? encode_d01(patch_summary(mci)) :
    error("_encode_variant: unknown variant $variant (use :min or :aug)")

"""
    rmse_report(dir; master_seed, artifact_path = DEFAULT_ADVI_ARTIFACT,
                model_path = DEFAULT_NPE_MODEL, N = 2000) -> NamedTuple

Score the trained NPE against the ADVI artifact in ρ-space on the reserved holdout
at `dir` (NPE-02). For each artifact pair, both stacks (sample AND control) are
joined to their holdout column by `global_index` (D-04) and scored against the known
ρ_true = `theta[1, ·]`:

  - **NPE**: per-parameter RMSE via `rmse(assess(est, θ_std, Z))` (standardized-space,
    scaled to physical units by `θzt.scale`); ρ_true is row 1 (the headline). 90%
    interval widths via `interval_width` on the un-standardized draws.
  - **ADVI**: ρ̂ = `mean(rho_sample/rho_control)` (already `ghat`-mapped per draw);
    ρ_true RMSE against the same known truth; 90% interval widths from the stored
    `rho_*_interval` bounds.
  - **Δρ (D-04)**: on identical pairs -- NPE via `delta_rho` (MC-diff of two passes),
    ADVI via MC-differencing the artifact ρ-draws, both vs `Δρ_true`.

Returns `(npe_rho_rmse, advi_rho_rmse, npe_all7_rmse, npe_interval_width,
advi_interval_width, delta_rho_rmse_npe, delta_rho_rmse_advi, n_stacks, n_pairs)`.
The pass/fail tolerance gate (D-09) lives in SC2, not here -- this reports numbers.
CPU-only (`use_gpu=false`); `master_seed` seeds the NPE sampling for reproducibility.
"""
function rmse_report(dir; master_seed, artifact_path = DEFAULT_ADVI_ARTIFACT,
                     model_path = DEFAULT_NPE_MODEL, N::Integer = 2000)
    Random.seed!(UInt32(master_seed & 0xFFFFFFFF))   # reproducible posterior draws (low 32 bits; IN-01)
    m   = load_npe(model_path)
    art = _load_artifact(artifact_path)
    ho  = load_holdout(dir)

    # global_index -> holdout column (strict join; D-04).
    gi_to_col = Dict(Int(ho.global_index[j]) => j for j in 1:length(ho.global_index))
    pairs  = art["pairs"]
    npairs = length(pairs)

    # --- Per-stack ρ-space list (sample and control alike), ADVI side + interval. ---
    stack_cols = Int[]
    advi_rho   = Float64[]
    advi_iv_w  = Float64[]
    for p in 1:npairs
        gis, gic = Int(pairs[p][1]), Int(pairs[p][2])
        haskey(gi_to_col, gis) || error("rmse_report: artifact global_index $gis absent from holdout at $dir")
        haskey(gi_to_col, gic) || error("rmse_report: artifact global_index $gic absent from holdout at $dir")
        push!(stack_cols, gi_to_col[gis]);  push!(advi_rho, mean(art["rho_sample"][p]))
        push!(advi_iv_w, art["rho_s_interval"][p][2] - art["rho_s_interval"][p][1])
        push!(stack_cols, gi_to_col[gic]);  push!(advi_rho, mean(art["rho_control"][p]))
        push!(advi_iv_w, art["rho_c_interval"][p][2] - art["rho_c_interval"][p][1])
    end

    θ_ho   = ho.theta[:, stack_cols]                # 7×K raw θ (joined order)
    ρ_true = Float64.(θ_ho[1, :])                   # known ρ_true per stack
    summ   = _holdout_summary(ho, m.variant)        # variant-correct raw summary (WR-01)
    Z      = standardize_summary(summ[:, stack_cols], m.zt, m.variant)  # d_in×K std

    # --- NPE per-parameter RMSE (assess/rmse, standardized) → physical units. ---
    θ_std = Float32.(StatsBase.transform(m.θzt, Float32.(θ_ho)))    # 7×K standardized
    a     = NeuralEstimators.assess(m.estimator, θ_std, Float32.(Z); use_gpu = false, N = N)
    rmse_std      = _rmse_vector(NeuralEstimators.rmse(a))          # 7-vector (std space)
    npe_all7_rmse = rmse_std .* _theta_scale(m.θzt)                 # physical units
    npe_rho_rmse  = npe_all7_rmse[1]                                # ρ_true is row 1

    # --- ADVI ρ_true RMSE (ρ-space via ghat), same joined stacks. ---
    advi_rho_rmse = sqrt(mean((advi_rho .- ρ_true) .^ 2))

    # --- 90% interval widths (ρ_true, row 1), both methods. ---
    npe_iv_w = Float64[]
    for j in 1:length(stack_cols)
        draws = posterior_for(m.estimator, Z[:, j]; N = N, use_gpu = false)
        orig  = StatsBase.reconstruct(m.θzt, draws)
        push!(npe_iv_w, interval_width(orig)[1])
    end
    npe_interval_width  = mean(npe_iv_w)
    advi_interval_width = mean(advi_iv_w)

    # --- Δρ RMSE on identical pairs (D-04). ---
    Δnpe = Float64[]; Δadvi = Float64[]; Δtrue = Float64[]
    for p in 1:npairs
        cs = gi_to_col[Int(pairs[p][1])]
        cc = gi_to_col[Int(pairs[p][2])]
        Zs = standardize_summary(summ[:, cs], m.zt, m.variant)   # variant-correct (WR-01)
        Zc = standardize_summary(summ[:, cc], m.zt, m.variant)
        push!(Δnpe,  delta_rho(m.estimator, Zs, Zc, m.θzt; N = N, use_gpu = false))
        push!(Δadvi, mean(art["rho_sample"][p] .- art["rho_control"][p]))
        push!(Δtrue, ho.theta[1, cs] - ho.theta[1, cc])
    end
    delta_rho_rmse_npe  = sqrt(mean((Δnpe  .- Δtrue) .^ 2))
    delta_rho_rmse_advi = sqrt(mean((Δadvi .- Δtrue) .^ 2))

    return (npe_rho_rmse        = npe_rho_rmse,
            advi_rho_rmse       = advi_rho_rmse,
            npe_all7_rmse       = npe_all7_rmse,
            npe_interval_width  = npe_interval_width,
            advi_interval_width = advi_interval_width,
            delta_rho_rmse_npe  = delta_rho_rmse_npe,
            delta_rho_rmse_advi = delta_rho_rmse_advi,
            n_stacks            = length(stack_cols),
            n_pairs             = npairs)
end

"""
    _cuda_available() -> Bool

True only if a CUDA package is already loaded AND reports a functional device. Never
forces a CUDA import (D-10 / Pitfall 1) -- the CPU path is the gate; GPU is a bonus.
"""
function _cuda_available()
    for (id, mod) in Base.loaded_modules
        if id.name == "CUDA"
            try
                return Bool(Base.invokelatest(getfield(mod, :functional)))
            catch
                return false
            end
        end
    end
    return false
end

"""
    _time_npe_pair(est, mci_s, mci_c, zt, variant, bench_N; bench_seconds) -> Float64

The D-08 NPE clock for ONE (sample, control) pair: `@belapsed` over {summary
extraction (`patch_summary` → `encode_d01`) → frozen-`zt` standardization → TWO
`sampleposterior` forward passes}. The raw `MultiChannelImage`s are re-simulated
OUTSIDE the timed block (setup), so the clock starts from the raw stack and measures
summary extraction + forward pass ONLY -- training is EXCLUDED (the amortized claim).

`bench_N` is the number of posterior draws the timed forward pass instantiates: the
amortized inference cost is the network/flow FORWARD PASS (an O(1), N-independent
network evaluation that conditions the flow); drawing draws from the conditioned flow
is O(N) POST-inference sampling -- the direct analog of ADVI's `rand(q, ·)`, which the
`vi()`-only ADVI clock likewise EXCLUDES. `bench_N` is therefore kept lightweight so
the O(N) draw cost is negligible on both sides; accuracy (RMSE/intervals) is scored
separately at the full `N`.

WR-03 -- CLOCK IS DELIBERATELY CONSERVATIVE ON THE NPE SIDE, NOT SYMMETRIC: the shared
summary-extraction cost (`patch_summary` → encode) is INSIDE this NPE timed block but is
EXCLUDED from the ADVI `wall_clock` (`advi_pair` builds `coloc_model`, which runs
`_prepare_data`/patching, BEFORE its `t0`). Charging the identical extraction only to the
NPE makes the NPE look slower, so the reported >100× speedup is if anything UNDERSTATED
(a conservative lower bound) -- the per-pair numbers are intentionally NOT apples-to-apples
on the extraction term. `@belapsed` reports the MINIMUM elapsed, filtering transient
contention. `use_gpu=false` (CPU-only, D-10).
"""
function _time_npe_pair(est, mci_s, mci_c, zt, variant, bench_N; bench_seconds = 0.5)
    return @belapsed begin
        Zs = standardize_summary(_encode_variant($mci_s, $variant), $zt, $variant)
        Zc = standardize_summary(_encode_variant($mci_c, $variant), $zt, $variant)
        posterior_for($est, Zs; N = $bench_N, use_gpu = false)
        posterior_for($est, Zc; N = $bench_N, use_gpu = false)
    end seconds=bench_seconds
end

"""
    speedup_report(dir; master_seed, artifact_path = DEFAULT_ADVI_ARTIFACT,
                   model_path = DEFAULT_NPE_MODEL, N = 2000, bench_N = 50,
                   bench_seconds = 0.5) -> NamedTuple

The NPE-03 wall-clock speedup, ALWAYS paired with the NPE-02 RMSE result (D-08 --
speedup is never reported alone). For each artifact pair the NPE clock (`t_npe`) is
measured from the RE-SIMULATED raw holdout stacks (D-04 shared re-sim path) as summary
extraction + forward pass, and the ADVI clock (`t_advi`) is the per-pair `vi()`
`wall_clock` read from the artifact (steady-state, JIT-excluded; 04-04).
`speedup = t_advi / t_npe` per pair; the headline is the median.

"TO-POSTERIOR" CLOCK, NPE-CONSERVATIVE (D-08; WR-03): the ADVI `wall_clock` times the
`vi()` optimization that PRODUCES the fitted posterior `q` and EXCLUDES both the
subsequent `rand(q, 100_000)` draw (04-04 / INFER-3) AND the upstream summary/patch
extraction (done in `build_coloc_model` before its timer). The NPE analog is the
amortized FORWARD PASS that produces the conditioned posterior -- but the NPE clock
additionally INCLUDES the shared summary-extraction cost (`_time_npe_pair` times it
inside `@belapsed`), so the two clocks are NOT symmetric on that term: the NPE is charged
for extraction and ADVI is not. This is deliberately conservative -- it can only shrink
the reported speedup, never inflate it -- so the >100× headline is a lower bound. The
timed forward pass draws a lightweight `bench_N` posterior (the forward pass is O(1);
bulk draws are the excluded analog of ADVI's `rand`). Accuracy (RMSE + intervals) is
scored at the full `N` in the paired `rmse` result. `full_median_speedup` additionally reports the
speedup when the timed NPE clock draws the full `N` posterior (a conservative
lower bound; the `vi()`-only ADVI clock itself under-counts real per-dataset ADVI,
which also runs a 100k prior chain + 100k posterior draw -- the "minutes" in ms-vs-min).

The active `Threads.nthreads()` is recorded so SC3 can assert it equals `BENCH_THREADS`
for the headline (D-13; >100× stated at one fixed, reported thread count). GPU timing is
a separate opt-in branch, populated only if a functional CUDA device is present
(`gpu_speedup`), never gating (D-10).

Returns `(median_speedup, full_median_speedup, speedups, t_npe, t_advi, bench_N, N,
threads, use_gpu, gpu_available, gpu_speedup, n_pairs, rmse)`.
"""
function speedup_report(dir; master_seed, artifact_path = DEFAULT_ADVI_ARTIFACT,
                        model_path = DEFAULT_NPE_MODEL, N::Integer = 2000,
                        bench_N::Integer = 50, bench_seconds = 0.5)
    m   = load_npe(model_path)
    art = _load_artifact(artifact_path)
    ho  = load_holdout(dir)

    gi_to_col = Dict(Int(ho.global_index[j]) => j for j in 1:length(ho.global_index))
    pairs  = art["pairs"]
    npairs = length(pairs)
    t_advi = Float64.(art["wall_clock"])

    t_npe      = Vector{Float64}(undef, npairs)   # forward-pass latency (bench_N draws)
    t_npe_full = Vector{Float64}(undef, npairs)   # full-N posterior latency (transparency)
    speedups   = Vector{Float64}(undef, npairs)
    full_speedups = Vector{Float64}(undef, npairs)
    for p in 1:npairs
        js = gi_to_col[Int(pairs[p][1])]           # holdout entry index (column) for re-sim
        jc = gi_to_col[Int(pairs[p][2])]
        # Re-simulate the raw stacks OUTSIDE the timed block (setup; D-08 clock starts
        # from the raw MultiChannelImage). Shared re-sim path with the ADVI baseline (D-04).
        mci_s = resimulate_holdout(dir, js; master_seed = master_seed).mci_sample
        mci_c = resimulate_holdout(dir, jc; master_seed = master_seed).mci_sample
        t_npe[p]      = _time_npe_pair(m.estimator, mci_s, mci_c, m.zt, m.variant, bench_N;
                                       bench_seconds = bench_seconds)
        t_npe_full[p] = _time_npe_pair(m.estimator, mci_s, mci_c, m.zt, m.variant, N;
                                       bench_seconds = bench_seconds)
        speedups[p]      = t_advi[p] / t_npe[p]
        full_speedups[p] = t_advi[p] / t_npe_full[p]
    end

    # Optional GPU bonus (D-10): only if a functional CUDA device is already loaded.
    gpu_speedup = nothing
    if _cuda_available()
        gpu_t = Vector{Float64}(undef, npairs)
        for p in 1:npairs
            js = gi_to_col[Int(pairs[p][1])]; jc = gi_to_col[Int(pairs[p][2])]
            mci_s = resimulate_holdout(dir, js; master_seed = master_seed).mci_sample
            mci_c = resimulate_holdout(dir, jc; master_seed = master_seed).mci_sample
            gt = @belapsed begin
                Zs = standardize_summary(_encode_variant($mci_s, $(m.variant)), $(m.zt), $(m.variant))
                Zc = standardize_summary(_encode_variant($mci_c, $(m.variant)), $(m.zt), $(m.variant))
                posterior_for($(m.estimator), Zs; N = $bench_N, use_gpu = true)
                posterior_for($(m.estimator), Zc; N = $bench_N, use_gpu = true)
            end seconds=bench_seconds
            gpu_t[p] = t_advi[p] / gt
        end
        gpu_speedup = median(gpu_t)
    end

    # The RMSE result the speedup is paired with (D-08 -- never speedup alone).
    rr = rmse_report(dir; master_seed = master_seed, artifact_path = artifact_path,
                     model_path = model_path, N = N)

    return (median_speedup      = median(speedups),
            full_median_speedup = median(full_speedups),
            speedups            = speedups,
            t_npe               = t_npe,
            t_npe_full          = t_npe_full,
            t_advi              = t_advi,
            bench_N             = bench_N,
            N                   = N,
            threads             = Threads.nthreads(),
            use_gpu             = false,
            gpu_available       = _cuda_available(),
            gpu_speedup         = gpu_speedup,
            n_pairs             = npairs,
            rmse                = rr)
end
