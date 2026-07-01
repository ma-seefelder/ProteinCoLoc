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

using JLD2         # read the cross-env ADVI artifact + the trained NPE
using Statistics   # mean, median
using Random       # reproducible NPE posterior sampling
using StatsBase    # transform / reconstruct (frozen θ-standardization)
using NeuralEstimators   # assess, rmse

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
    Random.seed!(UInt32(master_seed % 0xFFFFFFFF))   # reproducible posterior draws
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
    Z      = standardize_summary(ho.summary_min[:, stack_cols], m.zt, m.variant)  # 128×K std

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
        Zs = standardize_summary(ho.summary_min[:, cs], m.zt, m.variant)
        Zc = standardize_summary(ho.summary_min[:, cc], m.zt, m.variant)
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
