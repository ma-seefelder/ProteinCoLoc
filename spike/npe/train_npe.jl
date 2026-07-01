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

# spike/npe/train_npe.jl --- Phase-4 CPU-only fixed-data NPE training (NPE-01, D-05).
#
# Trains the `build_estimator` PosteriorEstimator on the Phase-3 leak-free loader
# tensors using the FIXED-DATA `train` form (RESEARCH Pattern 1) -- one dataset =
# one column of the standardized 128/142-dim summary. The online sampler/simulator
# `train` form is deliberately NOT used: re-simulating online would bypass the
# leak-free loader and the reserved-holdout exclusion (D-10).
#
# CPU-ONLY (D-10 / Pitfall 1): `use_gpu` DEFAULTS TO TRUE in NeuralEstimators, so
# EVERY `train` call here passes `use_gpu = false`. CUDA is never imported; the
# runtests CPU-only gate (and this file's own verify) asserts CUDA stays out of
# `Base.loaded_modules`.
#
# θ-STANDARDIZATION (RESEARCH Pitfall 5, A6): the loader standardizes only Z, not θ.
# θ components span wildly different ranges (ρ∈[-1,1], spillover∈[0,0.2], …), which
# hurts the flow's standard-Gaussian base. We fit a `ZScoreTransform` on θtr ONLY
# (leak-free, mirroring loader.jl L182 `fit(...; dims=2)`), train in standardized
# θ-space, and FREEZE that transform alongside the model so infer.jl can
# un-standardize draws before any ρ read (`StatsBase.reconstruct`). The θ-transform
# is the θ-space analog of the loader's frozen `zt` for Z.
#
# PERSISTENCE (cache.jl idiom): `save_npe` writes {estimator, θ-transform, zt, meta}
# to a `.tmp`, reopens to integrity-check, then `mv(...; force=true)` (filesystem-
# atomic). `load_npe` is the paired reader. The trained model + BOTH frozen
# transforms are the reproducibility contract the benchmark (04-05) / ablation
# (04-06) consume.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; touches no src/. Guarded
# includes keep the file loadable standalone AND inside runtests.jl.

using NeuralEstimators   # train
using Flux               # Flux.Optimisers.AdamW
using StatsBase          # ZScoreTransform, fit, transform, reconstruct
using JLD2               # atomic model persistence

# --- ORDER MATTERS: architecture first (build_estimator + the NPE_* consts), then
#     the loader (load_fold, the SOLE Z standardization path). Guarded for idempotency.
isdefined(@__MODULE__, :build_estimator) || include(joinpath(@__DIR__, "architecture.jl"))
isdefined(@__MODULE__, :load_fold)       || include(joinpath(@__DIR__, "..", "data", "loader.jl"))

# Phase-4 master seed (mirrors the pre-registered constant in test_npe.jl). Guarded
# so a redefinition to the SAME value under runtests.jl is a silent no-op.
if !isdefined(@__MODULE__, :NPE_MASTER_SEED)
    const NPE_MASTER_SEED = 0xC0FFEE
end

# Persistence schema for trained_npe.jld2 (bump on any breaking layout change).
const NPE_MODEL_SCHEMA = 1

"""
    fit_theta_transform(θtr::AbstractMatrix) -> ZScoreTransform

Fit a leak-free per-parameter `ZScoreTransform` on the TRAIN θ columns ONLY
(`dims = 2`, mirroring loader.jl's Z fit). All 7 θ rows are continuous, so every
row is standardized (no mask bypass, unlike the Z summaries). The returned
transform is FROZEN and persisted so posterior draws can be un-standardized with
`StatsBase.reconstruct` before any ρ read (Pitfall 5).
"""
fit_theta_transform(θtr::AbstractMatrix) = fit(ZScoreTransform, θtr; dims = 2)

"""
    train_fold(dir, fold::Integer; master_seed = NPE_MASTER_SEED, variant = :min,
               use_gpu::Bool = false, epochs::Integer = 200, batchsize::Integer = 64,
               dstar = NPE_DSTAR, depth = NPE_DEPTH, width = NPE_WIDTH,
               savepath::Union{Nothing,AbstractString} = nothing,
               verbose::Bool = false) -> NamedTuple

Train the Phase-4 NPE on fold `fold` of the cache at `dir`, CPU-only.

Loads the leak-free fold tensors via `load_fold` (Z standardized by the loader's
train-only `zt`, mask rows bypassed), fits a leak-free θ `ZScoreTransform` on the
train θ columns, standardizes θtr/θva, builds `build_estimator(size(Ztr,1), ...)`,
and calls the FIXED-DATA `train(est, θtr_std, θva_std, Ztr, Zva; use_gpu=false, ...)`
with AdamW weight decay (regularisation for fixed data; RESEARCH Pattern 1).

Returns `(estimator, θzt, zt, variant, d_in)` where `θzt` is the frozen θ-transform
and `zt` is the loader's frozen Z-transform. If `savepath` is given, the result is
persisted via `save_npe` (the atomic cache.jl idiom). `use_gpu=false` is passed on
every NeuralEstimators call (D-10 / Pitfall 1).
"""
function train_fold(dir, fold::Integer; master_seed = NPE_MASTER_SEED,
                    variant::Symbol = :min, use_gpu::Bool = false,
                    epochs::Integer = 200, batchsize::Integer = 64,
                    dstar::Integer = NPE_DSTAR, depth::Integer = NPE_DEPTH,
                    width::Integer = NPE_WIDTH,
                    savepath::Union{Nothing,AbstractString} = nothing,
                    verbose::Bool = false)
    use_gpu && throw(ArgumentError("train_fold: use_gpu=true is out of scope (D-10 CPU-only gate)"))

    fold_data = load_fold(dir, fold; master_seed = master_seed, variant = variant)

    # Leak-free θ-standardization (Pitfall 5): fit on train θ only, apply to both.
    θzt      = fit_theta_transform(fold_data.θtr)
    θtr_std  = Float32.(StatsBase.transform(θzt, fold_data.θtr))
    θva_std  = Float32.(StatsBase.transform(θzt, fold_data.θva))

    d_in = size(fold_data.Ztr, 1)
    est  = build_estimator(d_in; dstar = dstar, depth = depth, width = width)

    # FIXED-DATA train form on the leak-free cache. use_gpu=false on EVERY call.
    # NB: AdamW args are Float64 to match NeuralEstimators' Float64 CosAnneal
    # lr_schedule -- a Float32-parameterised optimiser state trips `Optimisers.adjust!`
    # when the schedule feeds a Float64 eta.
    est = train(est, θtr_std, θva_std, fold_data.Ztr, fold_data.Zva;
                epochs = epochs, batchsize = batchsize, use_gpu = false,
                optimiser = Flux.Optimisers.AdamW(5e-4, (0.9, 0.999), 1e-4),
                stopping_epochs = 10, verbose = verbose)

    result = (estimator = est, θzt = θzt, zt = fold_data.zt,
              variant = variant, d_in = d_in)

    savepath === nothing || save_npe(savepath, result;
                                     master_seed = master_seed, fold = fold)
    return result
end

"""
    save_npe(path, result::NamedTuple; master_seed, fold) -> String

Atomically persist a trained NPE `result` (as returned by `train_fold`) to `path`:
`jldsave` {estimator, θzt, zt, variant, d_in, meta} to a `.tmp`, reopen to
integrity-check the estimator key, then `mv(...; force=true)` (the filesystem-atomic
cache.jl idiom). A crash before the `mv` leaves only a discardable `.tmp`, never a
half-written model. Returns `path`.
"""
function save_npe(path, result::NamedTuple; master_seed, fold)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp;
        schema_version = NPE_MODEL_SCHEMA,
        estimator = result.estimator,
        theta_transform = result.θzt,
        zt = result.zt,
        variant = result.variant,
        d_in = result.d_in,
        meta = (master_seed = master_seed, fold = fold,
                dstar = NPE_DSTAR, depth = NPE_DEPTH, width = NPE_WIDTH,
                generated = string(Dates_now())))
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "estimator") "save_npe: integrity check failed, $tmp missing estimator"
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end

# Tiny timestamp helper (avoids a Dates dependency in the header `using` list while
# still tagging the artifact with a generation time).
Dates_now() = Base.Libc.strftime("%Y-%m-%dT%H:%M:%S", time())

"""
    load_npe(path) -> NamedTuple

Reload a trained NPE persisted by `save_npe`. Returns
`(estimator, θzt, zt, variant, d_in, meta)`. Integrity-checks the estimator key on
open (the cross-artifact reproducibility contract the benchmark/ablation consume).
"""
function load_npe(path)
    isfile(path) || error("load_npe: no model file at $path")
    d = JLD2.load(path)
    haskey(d, "estimator") || error("load_npe: $path missing estimator key")
    return (estimator = d["estimator"], θzt = d["theta_transform"], zt = d["zt"],
            variant = d["variant"], d_in = d["d_in"], meta = get(d, "meta", nothing))
end
