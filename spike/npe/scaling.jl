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

# spike/npe/scaling.jl --- Phase-4 amortization SCALING characterization (NPE-03 / D-12).
#
# The amortized thesis IS a scaling claim: the >100× is presented as a CURVE over
# (a) the number of datasets N and (b) the input size (imsize / pixel count), for
# BOTH NPE inference and the ADVI baseline, with empirical scaling exponents fit by
# a log-log least-squares slope. This is CHARACTERIZATION (reported), not a pass/fail
# gate -- the honest, defensible replacement for a single cherry-picked speedup point.
#
# THE AMORTIZATION CURVE OVER N (D-12, the textbook framing). Total wall-clock to go
# from a COLD start to N calibrated posteriors:
#   npe_time_N(N)  = C_train + N · t_fwd_per_dataset     (one-time training amortized
#                                                          over every future dataset)
#   advi_time_N(N) = N · t_advi_per_dataset              (NO amortization -- every
#                                                          dataset pays a full vi() run)
# In log-log over N the ADVI curve is a pure line (exponent ≈ 1.0); the NPE curve is
# FLATTENED by the fixed C_train offset (exponent < 1.0 over any finite grid whenever
# C_train > 0). That gap -- NPE materially flatter than ADVI over N -- IS the amortized
# O(1)/dataset-post-training claim made quantitative. The two curves cross at
# N* = C_train / (t_advi_per_dataset − t_fwd_per_dataset); beyond it the amortized NPE
# wins by the >100× per-dataset ratio the 04-05 headline reports. Reporting the whole
# curve (not just the asymptote) keeps the claim honest -- at very small N the ADVI
# baseline is realistically fast (~0.5 s, 04-05 caveat), so the win is asymptotic.
#
# THE SCALING CURVE OVER INPUT SIZE (D-12). Both methods must extract the SAME fixed
# 8×8 patch-correlation summary from the raw stack; that summary extraction is the ONLY
# imsize-DEPENDENT cost (O(pixels)). The imsize-INDEPENDENT cost is the amortized NPE
# forward pass (tiny) vs the ADVI vi() optimization (large, operates on the fixed 64-dim
# summary, not raw pixels). So:
#   npe_time_sz(sz)  = t_summary(sz) + t_fwd            (t_fwd imsize-independent)
#   advi_time_sz(sz) = t_summary(sz) + t_vi_fixed       (t_vi_fixed from the artifact)
# where t_summary(sz) is MEASURED at each size and t_vi_fixed is the real per-stack
# vi() cost read from advi_artifact.jld2 (04-04). This decomposes ADVI's cost into its
# imsize-dependent (shared summary) and imsize-independent (VI) parts using the real
# artifact for the latter -- an honest MODEL (ADVI itself cannot run in the pinned spike
# env, D-01), clearly labelled. The interesting D-12 finding is on the N axis; the imsize
# axis honestly shows ADVI's pixel-scaling is masked by its large fixed VI cost.
#
# REUSE, DON'T DUPLICATE (plan): the @belapsed timing primitives come from benchmark.jl
# (`_time_npe_pair`, the D-08 to-posterior clock convention: forward pass at a lightweight
# bench_N, bulk draws excluded on both sides). scaling.jl composes them across grid points.
#
# CPU-ONLY (D-10 / Pitfall 1): `use_gpu = false` on every NeuralEstimators call.
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; reaches src/ only transitively
# through the guarded benchmark.jl → infer.jl / resimulate.jl chains.

using JLD2           # atomic curve persistence (cache.jl .tmp+reopen+mv idiom)
using Statistics     # median, mean
using BenchmarkTools # @belapsed (via the benchmark.jl harness)

# benchmark.jl brings the whole chain into scope: _time_npe_pair / load_npe /
# load_holdout / _load_artifact / resimulate_holdout / standardize_summary /
# posterior_for / patch_summary / encode_d01 / sample_prior / simulate_pair /
# build_mci / sample_rng. Guarded for idempotency (standalone + runtests.jl).
isdefined(@__MODULE__, :speedup_report) || include(joinpath(@__DIR__, "benchmark.jl"))

# --- PRE-REGISTERED scaling grids (declared BEFORE the reported run; D-12/A2-style) ---
# Claude's discretion (plan), but LOCKED here as named consts so the reported curve
# cannot be tuned to a desired exponent. Guarded (mirrors ablation.jl) so a standalone
# include still gets the locked values while a harness that declared them first is respected.
if !isdefined(@__MODULE__, :SCALING_N_GRID)
    const SCALING_N_GRID        = [1, 5, 10, 20, 40]          # dataset-count axis (D-12a)
end
if !isdefined(@__MODULE__, :SCALING_IMSIZE_GRID)
    const SCALING_IMSIZE_GRID   = [(128, 128), (256, 256),    # input-size axis (D-12b);
                                   (512, 512), (1024, 1024)]  # every dim ≥ 64 (forward.jl)
end
if !isdefined(@__MODULE__, :SCALING_BENCH_N)
    const SCALING_BENCH_N       = 50                          # forward-pass draws (04-05 bench_N)
end
if !isdefined(@__MODULE__, :SCALING_BENCH_SECONDS)
    const SCALING_BENCH_SECONDS = 0.4                         # @belapsed budget per point
end
if !isdefined(@__MODULE__, :SCALING_TRAIN_EPOCHS)
    const SCALING_TRAIN_EPOCHS  = 200                         # reported C_train measurement epochs
end

const DEFAULT_SCALING_N_PATH      = joinpath(@__DIR__, "scaling_N.jld2")
const DEFAULT_SCALING_IMSIZE_PATH = joinpath(@__DIR__, "scaling_imsize.jld2")
const DEFAULT_SCALING_PATH        = joinpath(@__DIR__, "scaling_curves.jld2")

"""
    _loglog_slope(xs, ys) -> Float64

Empirical scaling exponent: the ordinary least-squares slope of `log(ys)` on
`log(xs)` (the exponent `p` in `y ∝ x^p`). Requires ≥2 points with distinct `xs`
and strictly positive `xs`/`ys`. Hand-rolled closed-form LSQ (no extra dependency);
this is the standard log-log complexity-fit, not a bespoke error metric.
"""
function _loglog_slope(xs, ys)
    lx = log.(Float64.(collect(xs)))
    ly = log.(Float64.(collect(ys)))
    n  = length(lx)
    n >= 2 || return NaN
    mx = sum(lx) / n
    my = sum(ly) / n
    sxx = sum((lx .- mx) .^ 2)
    sxx > 0 || return NaN                # xs all equal → exponent undefined
    return sum((lx .- mx) .* (ly .- my)) / sxx
end

"""
    _atomic_save_curve(path; kwargs...) -> String

Persist scaling curves to `path` via the cache.jl atomic idiom: `jldsave` every
keyword field to a `.tmp`, REOPEN to integrity-check a sentinel key, then
`mv(...; force = true)` (a filesystem-atomic rename). A crash before the `mv`
leaves only a discardable `.tmp`, never a half-written curve file.
"""
function _atomic_save_curve(path; kwargs...)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f              # integrity check before commit
        @assert haskey(f, "schema_version") "scaling: curve integrity check failed: $tmp"
    end
    mv(tmp, path; force = true)              # atomic commit
    return path
end

"""
    _simulate_at_imsize(master_seed, key, sz) -> MultiChannelImage

Simulate one raw 2-channel stack at a FORCED image size `sz` (bypassing the cost-aware
`sample_imsize` categorical so the size is exactly the grid point), keyed reproducibly by
`(master_seed, key)` through the UNCHANGED Phase-2 simulator (`sample_prior` →
`simulate_pair(rng, θ; imsize = sz)` → `build_mci`). Used only to time the imsize-dependent
summary-extraction cost -- θ content is irrelevant to the timing.
"""
function _simulate_at_imsize(master_seed::Integer, key::Integer, sz::Tuple{Int,Int})
    rng = sample_rng(master_seed, key)
    θ   = sample_prior(rng)
    return build_mci(simulate_pair(rng, θ; imsize = sz))
end

"""
    scaling_over_N(dir; master_seed, Ns = SCALING_N_GRID, train_cost = nothing,
                   artifact_path = DEFAULT_ADVI_ARTIFACT, model_path = DEFAULT_NPE_MODEL,
                   bench_N = SCALING_BENCH_N, bench_seconds = SCALING_BENCH_SECONDS,
                   train_epochs = SCALING_TRAIN_EPOCHS, save_path = nothing) -> NamedTuple

The amortization curve over the dataset count N (D-12a). Measures the per-dataset NPE
forward-pass latency once (via the benchmark.jl `_time_npe_pair` D-08 clock, halved to
per-stack) and the per-dataset ADVI cost from the artifact `wall_clock` (median per pair,
halved to per-stack), then builds the total-wall-clock-from-cold curves

    npe_time_N  = C_train + Ns .* t_fwd_per_dataset      (one-time training amortized)
    advi_time_N = Ns .* t_advi_per_dataset               (linear -- no amortization)

and fits a log-log slope of each. WR-04 -- HONESTY NOTE: `npe_exponent_N`/`advi_exponent_N`
are the slopes of these two CLOSED-FORM ANALYTIC curves, NOT empirical exponents measured
from repeated end-to-end runs. Because `advi_time_N` is exactly linear and `npe_time_N` is
affine with a positive `C_train` intercept, `npe_exponent_N < advi_exponent_N` holds BY
CONSTRUCTION for any `C_train > 0` -- it is an analytic amortization MODEL, not a
measurement. The MEASURED inputs are `t_fwd_per_dataset`, `t_advi_per_dataset` and
(optionally) `C_train`; those, and the `crossover_N` they imply, are the honest quantities
to interrogate (see the SC3 test, which asserts on them rather than on the tautological
exponent ordering). `C_train` (the amortized NPE training investment) is either supplied via
`train_cost` or MEASURED by timing a single `train_fold(dir, 1; epochs = train_epochs)` when
`train_cost === nothing` (the reported run; requires a full CV cache at `dir`). Returns
`(N_grid, npe_time_N, advi_time_N, npe_exponent_N, advi_exponent_N, train_cost,
t_fwd_per_dataset, t_advi_per_dataset, crossover_N, threads)`. `use_gpu = false`
throughout. Characterization only -- no pass/fail here.
"""
function scaling_over_N(dir; master_seed, Ns = SCALING_N_GRID, train_cost = nothing,
                        artifact_path = DEFAULT_ADVI_ARTIFACT,
                        model_path = DEFAULT_NPE_MODEL,
                        bench_N::Integer = SCALING_BENCH_N,
                        bench_seconds = SCALING_BENCH_SECONDS,
                        train_epochs::Integer = SCALING_TRAIN_EPOCHS,
                        save_path = nothing)
    m   = load_npe(model_path)
    art = _load_artifact(artifact_path)

    # Per-dataset ADVI cost: the artifact vi() wall_clock is per (sample,control) PAIR,
    # so a single stack is half a pair's cost (the paired joint vi() amortizes trivially
    # over its two stacks). Median over pairs (steady-state, JIT-excluded; 04-04).
    t_advi_per_dataset = median(Float64.(art["wall_clock"])) / 2

    # Per-dataset NPE forward-pass latency via the SHARED benchmark.jl D-08 clock
    # (_time_npe_pair times TWO summary+forward passes; halve for one dataset).
    mci_s = resimulate_holdout(dir, 1; master_seed = master_seed).mci_sample
    mci_c = resimulate_holdout(dir, 2; master_seed = master_seed).mci_sample
    t_pair = _time_npe_pair(m.estimator, mci_s, mci_c, m.zt, m.variant, bench_N;
                            bench_seconds = bench_seconds)
    t_fwd_per_dataset = t_pair / 2

    # One-time amortized training investment C_train: supplied, or measured by timing a
    # single fold train (reported run; needs a full CV cache at `dir`).
    C_train = train_cost === nothing ?
        (@elapsed train_fold(dir, 1; master_seed = master_seed, variant = m.variant,
                             epochs = train_epochs, verbose = false)) :
        Float64(train_cost)

    Ngrid       = collect(Ns)
    npe_time_N  = C_train .+ Ngrid .* t_fwd_per_dataset       # amortized (train + N passes)
    advi_time_N = Ngrid .* t_advi_per_dataset                 # linear (per-dataset vi())

    # Slopes of the ANALYTIC model curves above (WR-04): npe < advi holds by construction
    # for C_train > 0 -- an amortization MODEL, not a measured empirical exponent.
    npe_exponent_N  = _loglog_slope(Ngrid, npe_time_N)        # < 1 (flattened by C_train)
    advi_exponent_N = _loglog_slope(Ngrid, advi_time_N)       # ≈ 1 (pure linear)

    # Amortization crossover: N beyond which the trained NPE's total cost undercuts ADVI.
    denom      = t_advi_per_dataset - t_fwd_per_dataset
    crossover_N = denom > 0 ? C_train / denom : Inf

    result = (N_grid            = Ngrid,
              npe_time_N         = npe_time_N,
              advi_time_N        = advi_time_N,
              npe_exponent_N     = npe_exponent_N,
              advi_exponent_N    = advi_exponent_N,
              train_cost         = C_train,
              t_fwd_per_dataset  = t_fwd_per_dataset,
              t_advi_per_dataset = t_advi_per_dataset,
              crossover_N        = crossover_N,
              threads            = Threads.nthreads())

    save_path === nothing || _atomic_save_curve(save_path;
        schema_version = 1, axis = "N", N_grid = Ngrid,
        npe_time_N = npe_time_N, advi_time_N = advi_time_N,
        npe_exponent_N = npe_exponent_N, advi_exponent_N = advi_exponent_N,
        train_cost = C_train, t_fwd_per_dataset = t_fwd_per_dataset,
        t_advi_per_dataset = t_advi_per_dataset, crossover_N = crossover_N,
        threads = Threads.nthreads(), master_seed = UInt64(master_seed))
    return result
end

"""
    scaling_over_imsize(dir; master_seed, sizes = SCALING_IMSIZE_GRID,
                        artifact_path = DEFAULT_ADVI_ARTIFACT, model_path = DEFAULT_NPE_MODEL,
                        bench_N = SCALING_BENCH_N, bench_seconds = SCALING_BENCH_SECONDS,
                        save_path = nothing) -> NamedTuple

The scaling curve over input size (D-12b). For each `sz` in `sizes` a raw stack is
simulated at that FORCED imsize and the imsize-DEPENDENT summary-extraction cost
`t_summary(sz) = @belapsed encode_d01(patch_summary(mci))` is measured. The imsize-
INDEPENDENT costs are measured/read once: the NPE forward-pass latency `t_fwd` (operates
on the fixed 128-dim summary) and the ADVI fixed vi() cost `t_vi_fixed` (per-stack median
from the artifact). The two curves are then

    npe_time_sz  = t_summary(sz) .+ t_fwd        (measured summary + tiny forward)
    advi_time_sz = t_summary(sz) .+ t_vi_fixed   (measured summary + real vi() fixed cost)

with pixel count `W·H` as the size axis for the exponent fit. `advi_time_sz` is an HONEST
MODEL: the shared summary term is measured at each size, the vi() term is the real
artifact cost (ADVI cannot run in the pinned spike env, D-01). Returns
`(imsize_grid, imsize_pixels, npe_time_sz, advi_time_sz, npe_exponent_sz,
advi_exponent_sz, t_summary, t_fwd, t_vi_fixed, threads)`. `use_gpu = false` throughout.
"""
function scaling_over_imsize(dir; master_seed, sizes = SCALING_IMSIZE_GRID,
                             artifact_path = DEFAULT_ADVI_ARTIFACT,
                             model_path = DEFAULT_NPE_MODEL,
                             bench_N::Integer = SCALING_BENCH_N,
                             bench_seconds = SCALING_BENCH_SECONDS,
                             save_path = nothing)
    m   = load_npe(model_path)
    art = _load_artifact(artifact_path)
    t_vi_fixed = median(Float64.(art["wall_clock"])) / 2      # imsize-independent vi() cost

    szv    = collect(sizes)
    pixels = [Int(sz[1]) * Int(sz[2]) for sz in szv]

    # imsize-DEPENDENT summary-extraction cost, measured at each grid point.
    t_summary = Vector{Float64}(undef, length(szv))
    for (i, sz) in enumerate(szv)
        mci = _simulate_at_imsize(master_seed, i, sz)
        t_summary[i] = @belapsed _encode_variant($mci, $(m.variant)) seconds = bench_seconds
    end

    # imsize-INDEPENDENT NPE forward-pass latency, measured once on a fixed summary vector.
    mci_ref = _simulate_at_imsize(master_seed, length(szv) + 1, first(szv))
    Zref    = standardize_summary(_encode_variant(mci_ref, m.variant), m.zt, m.variant)  # variant-correct (WR-01)
    t_fwd   = @belapsed posterior_for($(m.estimator), $Zref; N = $bench_N, use_gpu = false) seconds = bench_seconds

    npe_time_sz  = t_summary .+ t_fwd
    advi_time_sz = t_summary .+ t_vi_fixed

    npe_exponent_sz  = _loglog_slope(pixels, npe_time_sz)
    advi_exponent_sz = _loglog_slope(pixels, advi_time_sz)

    result = (imsize_grid      = szv,
              imsize_pixels     = pixels,
              npe_time_sz       = npe_time_sz,
              advi_time_sz      = advi_time_sz,
              npe_exponent_sz   = npe_exponent_sz,
              advi_exponent_sz  = advi_exponent_sz,
              t_summary         = t_summary,
              t_fwd             = t_fwd,
              t_vi_fixed        = t_vi_fixed,
              threads           = Threads.nthreads())

    save_path === nothing || _atomic_save_curve(save_path;
        schema_version = 1, axis = "imsize",
        imsize_grid = szv, imsize_pixels = pixels,
        npe_time_sz = npe_time_sz, advi_time_sz = advi_time_sz,
        npe_exponent_sz = npe_exponent_sz, advi_exponent_sz = advi_exponent_sz,
        t_summary = t_summary, t_fwd = t_fwd, t_vi_fixed = t_vi_fixed,
        threads = Threads.nthreads(), master_seed = UInt64(master_seed))
    return result
end

"""
    scaling_curves(dir; master_seed, save_path = DEFAULT_SCALING_PATH, kwargs...) -> NamedTuple

Run BOTH scaling axes (D-12) and merge into the single reported NamedTuple
`(N_grid, npe_time_N, advi_time_N, npe_exponent_N, advi_exponent_N, imsize_grid,
npe_time_sz, advi_time_sz, npe_exponent_sz, advi_exponent_sz)` (plus the auxiliary
timing fields), persisting the combined curves to `save_path` via the atomic idiom.
Any `Ns`/`sizes`/`train_cost`/timing kwargs pass through to the two sub-runs.
"""
function scaling_curves(dir; master_seed, save_path = DEFAULT_SCALING_PATH,
                        Ns = SCALING_N_GRID, sizes = SCALING_IMSIZE_GRID,
                        train_cost = nothing,
                        artifact_path = DEFAULT_ADVI_ARTIFACT,
                        model_path = DEFAULT_NPE_MODEL,
                        bench_N::Integer = SCALING_BENCH_N,
                        bench_seconds = SCALING_BENCH_SECONDS)
    rN = scaling_over_N(dir; master_seed = master_seed, Ns = Ns, train_cost = train_cost,
                        artifact_path = artifact_path, model_path = model_path,
                        bench_N = bench_N, bench_seconds = bench_seconds)
    rS = scaling_over_imsize(dir; master_seed = master_seed, sizes = sizes,
                             artifact_path = artifact_path, model_path = model_path,
                             bench_N = bench_N, bench_seconds = bench_seconds)

    merged = (N_grid           = rN.N_grid,
              npe_time_N        = rN.npe_time_N,
              advi_time_N       = rN.advi_time_N,
              npe_exponent_N    = rN.npe_exponent_N,
              advi_exponent_N   = rN.advi_exponent_N,
              imsize_grid       = rS.imsize_grid,
              npe_time_sz       = rS.npe_time_sz,
              advi_time_sz      = rS.advi_time_sz,
              npe_exponent_sz   = rS.npe_exponent_sz,
              advi_exponent_sz  = rS.advi_exponent_sz,
              train_cost        = rN.train_cost,
              crossover_N       = rN.crossover_N,
              imsize_pixels     = rS.imsize_pixels,
              threads           = Threads.nthreads())

    save_path === nothing || _atomic_save_curve(save_path;
        schema_version = 1, axis = "both",
        N_grid = rN.N_grid, npe_time_N = rN.npe_time_N, advi_time_N = rN.advi_time_N,
        npe_exponent_N = rN.npe_exponent_N, advi_exponent_N = rN.advi_exponent_N,
        imsize_grid = rS.imsize_grid, imsize_pixels = rS.imsize_pixels,
        npe_time_sz = rS.npe_time_sz, advi_time_sz = rS.advi_time_sz,
        npe_exponent_sz = rS.npe_exponent_sz, advi_exponent_sz = rS.advi_exponent_sz,
        train_cost = rN.train_cost, crossover_N = rN.crossover_N,
        threads = Threads.nthreads(), master_seed = UInt64(master_seed))
    return merged
end
