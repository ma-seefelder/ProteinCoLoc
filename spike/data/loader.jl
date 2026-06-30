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

# spike/data/loader.jl --- DATA-03 leak-free k-fold loader (D-07 / D-08 / D-09 / D-10).
#
# The SOLE standardization path of the pipeline. The cache is RAW (D-07): shards
# store the un-standardized (θ, summary) columns. ALL standardization happens here
# and ONLY here, fit on the TRAIN columns of a fold and applied to the held-out
# fold -- so cross-fold leakage is impossible by construction, NOT by convention
# (D-08). There is deliberately NO public `standardize_all` / global-stats symbol:
# the module exports only `load_main_pool`, `load_fold`, `load_holdout`,
# `all_folds`, so there is no API a caller could use to z-score the whole pool and
# leak val statistics into train. SC-3a asserts `:standardize_all ∉ names(Loader)`.
#
# K-FOLD (D-09): `fold_rng(master_seed ⊻ FOLD_SALT, 0)` (seeding.jl) deterministically
# permutes `1:N`; fold f = `perm[f:K:end]`, train = the complement. Two calls with
# the same `master_seed` yield IDENTICAL membership -- folds are reproducible and
# independent of the generation key stream.
#
# MASK BYPASS (D-01 / Anti-Pattern): the binary present/absent mask rows (65:128 of
# summary_min) are NOT z-scored -- z-scoring a 0/1 mask re-couples folds through the
# mask mean. Only the continuous value rows (1:64 of D-01 + the continuous moment
# rows of the aug variant) are standardized; the mask passes through UNCHANGED.
#
# HOLDOUT (D-10): `load_main_pool` reads ONLY `shard_*.jld2`, NEVER `holdout.jld2`,
# so the reserved ≥20-stack ADVI holdout is excluded from every fold STRUCTURALLY
# (the loader simply never sees it). `load_holdout` is the sole reader of that file.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local read path; touches no src/.
# `loader.jl` is wrapped in `module Loader` so `names(Loader)` is a controlled public
# surface (the D-08 guarantee); the trailing `using .Loader` re-exports the API so a
# plain `include("loader.jl")` makes `load_fold`/`load_main_pool` directly callable.

module Loader

using JLD2          # jldopen / load (column-major shard reader)
using StatsBase     # ZScoreTransform, fit, transform (the ONLY standardization)
using Random        # randperm (deterministic fold permutation)

# seeding.jl supplies FOLD_SALT + the keyed `fold_rng`. Included INSIDE the module
# so the fold key namespace lives in the Loader surface, not leaked to the caller.
include(joinpath(@__DIR__, "seeding.jl"))

export load_main_pool, load_fold, load_holdout, all_folds

"""
    _shard_files(dir) -> Vector{String}

The sorted absolute paths of the main-pool shards (`shard_*.jld2`) in `dir`.
`holdout.jld2` is EXCLUDED by the `shard_` prefix filter -- the structural basis of
the D-10 holdout exclusion (`load_main_pool` can never read the reserved set).
"""
function _shard_files(dir)
    files = filter(readdir(dir; join = true)) do f
        startswith(basename(f), "shard_") && endswith(f, ".jld2")
    end
    return sort(files)
end

"""
    load_main_pool(dir) -> NamedTuple

Load the RAW (un-standardized) main pool by `hcat`-ing every `shard_*.jld2` in
`dir` in lexical (== numeric, zero-padded) order. Returns column-major matrices
`(theta=7×N, summary_min=128×N, summary_aug=AUG_DIM×N, global_index=Vector{Int} N)`
where one sample = one column. `holdout.jld2` is NEVER read (D-10), so the column
count is exactly N, the size of the main pool.
"""
function load_main_pool(dir)
    files = _shard_files(dir)
    isempty(files) && error("load_main_pool: no shard_*.jld2 found in $dir")
    theta        = reduce(hcat, [JLD2.load(f, "theta")        for f in files])
    summary_min  = reduce(hcat, [JLD2.load(f, "summary_min")  for f in files])
    summary_aug  = reduce(hcat, [JLD2.load(f, "summary_aug")  for f in files])
    global_index = reduce(vcat, [JLD2.load(f, "global_index") for f in files])
    return (theta = theta, summary_min = summary_min,
            summary_aug = summary_aug, global_index = global_index)
end

"""
    load_holdout(dir) -> NamedTuple

Load the reserved ADVI holdout from `holdout.jld2` -- the SOLE reader of that file.
Returns the same column-major NamedTuple shape as `load_main_pool`. Kept strictly
separate from the main pool so the holdout is excluded from every CV fold by
construction (D-10).
"""
function load_holdout(dir)
    path = joinpath(dir, "holdout.jld2")
    isfile(path) || error("load_holdout: no holdout.jld2 found in $dir")
    d = JLD2.load(path)
    return (theta = d["theta"], summary_min = d["summary_min"],
            summary_aug = d["summary_aug"], global_index = d["global_index"])
end

"""
    _row_partition(variant::Symbol, nrows::Int) -> (cont_rows, mask_rows)

Split the feature rows of a summary variant into the CONTINUOUS rows (z-scored) and
the binary MASK rows (passed through unchanged):

  - `:min` (128 rows): cont = 1:64 (imputed correlations), mask = 65:128 (D-01 mask).
  - `:aug` (AUG_DIM rows): cont = [1:64; 129:end] (D-01 values + the D-02 continuous
    moments), mask = 65:128. The aug moments (rows 129:end) are all continuous.

z-scoring a 0/1 mask would re-couple folds through the mask mean (Anti-Pattern), so
the mask rows are the SAME index set for both variants and never standardized.
"""
function _row_partition(variant::Symbol, nrows::Int)
    if variant === :min
        nrows == 128 || error("_row_partition: :min expects 128 rows, got $nrows")
        return (collect(1:64), collect(65:128))
    elseif variant === :aug
        nrows >= 128 || error("_row_partition: :aug expects ≥128 rows, got $nrows")
        return (vcat(collect(1:64), collect(129:nrows)), collect(65:128))
    else
        error("_row_partition: unknown variant $variant (use :min or :aug)")
    end
end

"""
    all_folds(N::Integer; K::Int=5, master_seed) -> Vector{Vector{Int}}

The K disjoint validation index vectors of the deterministic k-fold split (D-09):
`perm = randperm(fold_rng(master_seed), N)`; fold f = `perm[f:K:end]`. Their union
is exactly `1:N` and every pair is disjoint (a single permutation, strided K ways).
Reproducible: same `master_seed` ⇒ same folds.
"""
function all_folds(N::Integer; K::Int = 5, master_seed)
    perm = randperm(fold_rng(master_seed), N)
    return [perm[f:K:end] for f in 1:K]
end

"""
    load_fold(dir, fold::Int; K::Int=5, master_seed, variant::Symbol=:min) -> NamedTuple

The leak-free per-fold loader -- the ONLY standardization entry point (D-07/D-08).
Loads the RAW main pool, assigns the deterministic k-fold split (D-09), fits a
`StatsBase.ZScoreTransform` on the TRAIN columns ONLY (continuous rows), applies it
to both train and the held-out fold, passes the binary mask rows (65:128) through
UNCHANGED, and returns Float32 `d×K` tensors plus the fit object:

  `(Ztr, θtr, Zva, θva, zt)`

`zt` is returned so Phase 5 can FREEZE the train-only preprocessing. Because the fit
never sees the val columns, no val statistic can leak into train -- by construction.
`variant` selects the summary matrix (`:min` 128-dim D-01, `:aug` AUG_DIM D-02).
"""
function load_fold(dir, fold::Int; K::Int = 5, master_seed, variant::Symbol = :min)
    pool = load_main_pool(dir)
    Z    = variant === :min ? pool.summary_min :
           variant === :aug ? pool.summary_aug :
           error("load_fold: unknown variant $variant (use :min or :aug)")
    θ    = pool.theta
    N    = size(Z, 2)
    1 <= fold <= K || error("load_fold: fold $fold out of 1:$K")

    perm      = randperm(fold_rng(master_seed), N)   # deterministic split (D-09)
    val_idx   = perm[fold:K:end]
    train_idx = setdiff(1:N, val_idx)

    cont_rows, mask_rows = _row_partition(variant, size(Z, 1))

    # ◄ FIT ON TRAIN ONLY (continuous rows). Never touches the val columns (D-07).
    zt = fit(ZScoreTransform, Z[cont_rows, train_idx]; dims = 2)

    Ztr = Matrix{Float64}(undef, size(Z, 1), length(train_idx))
    Zva = Matrix{Float64}(undef, size(Z, 1), length(val_idx))
    # continuous rows: standardized with the train-fit transform
    Ztr[cont_rows, :] = StatsBase.transform(zt, Z[cont_rows, train_idx])
    Zva[cont_rows, :] = StatsBase.transform(zt, Z[cont_rows, val_idx])
    # binary mask rows: passed through UNCHANGED (no z-scoring; D-01 / Anti-Pattern)
    Ztr[mask_rows, :] = Z[mask_rows, train_idx]
    Zva[mask_rows, :] = Z[mask_rows, val_idx]

    return (Ztr = Float32.(Ztr), θtr = Float32.(θ[:, train_idx]),
            Zva = Float32.(Zva), θva = Float32.(θ[:, val_idx]), zt = zt)
end

end # module Loader

# Re-export the public API so a plain `include("loader.jl")` makes the loader
# directly callable (`load_fold(...)`) while `names(Loader)` stays the controlled
# D-08 surface. NO `standardize_all` symbol exists anywhere -- leakage is impossible
# by construction, not by convention.
using .Loader
