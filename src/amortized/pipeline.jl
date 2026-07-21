#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder

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

#############################################################################################
# src/amortized/pipeline.jl --- the single reusable per-grid training pipeline (PROD-01/02).
#
# `_train_grid_pipeline(grid; ...)` factors the WHOLE per-grid sequence into ONE function that
# `registry.train_and_register(grid)` calls, so the user-definable-grid on-ramp (D-04) is REAL
# code, not per-grid prose:
#
#   datagen → leak-free zt fit → train_npe → train_ratio → fit_ood_nulls
#           → persist all three CPU-resident → EstimatorBundle
#
# IDEMPOTENCY (training-layer critical):
#   • SKIP-IF-DONE: if the three artifacts (npe_G / ratio_G / ood_nulls_G) already exist AND pass
#     their reopen-integrity check (`_estimator_ok`/`_ratio_ok`/`_ood_nulls_ok`, persist.jl), the
#     pipeline LOADS and returns them instead of retraining.
#   • ATOMIC writes: every artifact is written `.tmp` → reopen-integrity → `mv(...; force=true)`
#     (persist.jl), so an interruption mid-write leaves only a discardable `.tmp`.
#
# REPRODUCIBILITY SPLIT (D-06 / 07-RESEARCH §Reproducibility): training MAY use GPU
# (`use_gpu = has_cuda_device()`), but the persisted nets are CPU-resident (`Flux.state`,
# persist.jl) and the ship-gate + shipped inference run CPU-only, so what is pre-registered and
# reproduced is device-independent. Training-time non-determinism is acceptable.
#
# THREADING (-t auto): the `generate_samples` datagen step parallelizes with `Threads.@threads`
# ONLY when `Threads.nthreads() > 1`. Launch Julia with `-t auto` (e.g.
# `julia -t auto --project -e '...'`) or the datagen SILENTLY falls back to serial. Because every
# sample is a pure function of `(master_seed, grid, index)` (Philox-per-index, D-11/D-12), the
# parallel result is BYTE-IDENTICAL to the serial `-t 1` result — `-t auto` only changes wall
# time, never the trained pool.
#############################################################################################

# --- Per-grid defaults (small tables; a novel grid gets a sane bias) ------------------------

"""
    default_artifacts_root() -> String

The package-root `artifacts/` directory holding the per-grid stores (`grid_G/`). A later Phase-7
plan swaps this for an Artifacts.toml lazy-download root; here it is a plain on-disk directory.
"""
default_artifacts_root() = normpath(joinpath(@__DIR__, "..", "..", "artifacts"))

"""
    default_imsize_for(grid::Integer) -> NTuple{N,Tuple{Int,Int}}

Per-grid default image-size set biased so finer grids draw LARGER images (fine patches need more
pixels to clear the ≥15-survivor floor): the `IMSIZE_SET` entries whose SMALLER dimension is
≥ `clamp(32·grid, 256, 1024)`. So `grid ≤ 8` keeps the full 256²-heavy set, `default_imsize_for(16)`
biases ≥512², and `default_imsize_for(32)` biases ≥1024². A novel grid interpolates on the same
rule; an over-fine grid falls back to the largest available size.
"""
function default_imsize_for(grid::Integer)
    min_dim = clamp(32 * Int(grid), 256, 1024)
    set = Tuple{Int,Int}[sz for sz in IMSIZE_SET if min(sz[1], sz[2]) >= min_dim]
    isempty(set) && (set = Tuple{Int,Int}[IMSIZE_SET[end]])
    return Tuple(set)
end

"""
    default_npairs_for(grid::Integer) -> Int

Per-grid default training-pool size (number of simulated `(θ, summary)` samples). Finer grids
have a longer summary vector, so they get a larger default pool. Overridable via the `n_pairs`
keyword of `_train_grid_pipeline`.
"""
function default_npairs_for(grid::Integer)
    grid <= 4  && return 30_000
    grid <= 8  && return 50_000
    grid <= 16 && return 80_000
    return 120_000
end

# Placeholder (not-yet-gated) calibration report — the per-grid D-05 ship-gate plan fills the
# reliability curve + ECE/MCE and flips `passed`. Empty curve + NaN ECE/MCE signal "un-gated".
_placeholder_calibration(grid::Integer) =
    CalibrationMeta(Float64[], Float64[], Float64[], Int[], NaN, NaN, Int(grid),
                    (; gated = false, passed = false))

# Artifact paths for a grid's store.
_grid_dir(root, grid)       = joinpath(root, "grid_$(grid)")
_npe_path(root, grid)       = joinpath(_grid_dir(root, grid), "npe_$(grid).jld2")
_ratio_path(root, grid)     = joinpath(_grid_dir(root, grid), "ratio_$(grid).jld2")
_ood_nulls_path(root, grid) = joinpath(_grid_dir(root, grid), "ood_nulls_$(grid).jld2")

# Rebuild an EstimatorBundle from the three persisted artifacts (the SKIP-IF-DONE / lazy-load
# path). The persisted ood_nulls carries the frozen density (+noise) fits; the PP `:model` handle
# is re-composed by the inference layer from the loaded npe.
function _bundle_from_artifacts(grid, npe_path, ratio_path, ood_path)
    npe   = load_estimator(npe_path)
    ratio = load_ratio(ratio_path)
    nulls = load_ood_nulls(ood_path)
    return EstimatorBundle(Int(grid), npe.estimator, ratio, nulls, npe.zt, npe.θzt,
                           _placeholder_calibration(grid))
end

"""
    _train_grid_pipeline(grid::Integer; imsize_set = default_imsize_for(grid),
                         imsize_weights = nothing, use_gpu = has_cuda_device(),
                         n_pairs = default_npairs_for(grid), ratio_n = n_pairs,
                         master_seed = DEFAULT_MASTER_SEED,
                         artifacts_root = default_artifacts_root(), val_frac = 0.15,
                         npe_epochs = 300, ratio_epochs = 300, npe_kwargs = (;),
                         ratio_kwargs = (;), skip_if_done = true, datagen = nothing,
                         verbose = false) -> EstimatorBundle

Run the FULL per-grid pipeline ONCE for any `grid` and return a populated `EstimatorBundle`:

  1. **SKIP-IF-DONE:** if `skip_if_done` and the three artifacts (`npe_G`/`ratio_G`/`ood_nulls_G`)
     exist AND pass their reopen-integrity check, LOAD and return them (no retraining).
  2. **datagen:** `generate_samples(n_pairs; grid, master_seed, imsize_set, imsize_weights)` — a
     TRAIN pool of `(θ, summary)` pairs. Injectable via `datagen` (a zero-arg function returning
     `(; theta, summary_min)`) for testing without the forward simulator.
  3. **leak-free `zt`:** fit the summary `ZScoreTransform` on the TRAIN continuous rows only.
  4. **train_npe / train_ratio** (both `use_gpu`-plumbed), on the frozen-`zt` standardized pool.
  5. **fit_ood_nulls** on the TRAIN-ONLY ID pool (density channel).
  6. **persist** all three CPU-resident (`save_estimator`/`save_ratio`/`save_ood_nulls`) under
     `artifacts_root/grid_G/`, each carrying the TRAINING IMAGE-SIZE PROVENANCE
     (`imsize_set`/`imsize_weights`/`imsize_source`) in its `meta` — read it back with
     `training_imsize_provenance`. When `datagen` is INJECTED the keywords do not describe the
     pool, so the provenance is recorded as `:unknown` (never guessed).

The returned bundle's `calibration` is a NOT-YET-GATED placeholder the per-grid D-05 ship-gate
fills. **-t auto (peer requirement):** the datagen step parallelizes only when
`Threads.nthreads() > 1`, so launch Julia with `-t auto` to parallelize; the result is
BYTE-IDENTICAL to the serial `-t 1` run (Philox-per-index keying, D-11/D-12), so `-t auto` only
changes wall time, never the trained pool.
"""
function _train_grid_pipeline(grid::Integer;
                              imsize_set = default_imsize_for(grid),
                              imsize_weights = nothing,
                              use_gpu::Bool = has_cuda_device(),
                              n_pairs::Integer = default_npairs_for(grid),
                              ratio_n::Integer = n_pairs,
                              master_seed::Integer = DEFAULT_MASTER_SEED,
                              artifacts_root = default_artifacts_root(),
                              val_frac::Real = 0.15,
                              npe_epochs::Integer = 300, ratio_epochs::Integer = 300,
                              npe_kwargs = (;), ratio_kwargs = (;),
                              skip_if_done::Bool = true, datagen = nothing,
                              verbose::Bool = false)
    grid >= 1 || throw(ArgumentError("_train_grid_pipeline: grid must be ≥ 1, got $grid"))
    npe_path   = _npe_path(artifacts_root, grid)
    ratio_path = _ratio_path(artifacts_root, grid)
    ood_path   = _ood_nulls_path(artifacts_root, grid)

    # (1) SKIP-IF-DONE — all three artifacts present AND integrity-valid ⇒ load, do NOT retrain.
    if skip_if_done && _estimator_ok(npe_path) && _ratio_ok(ratio_path) && _ood_nulls_ok(ood_path)
        verbose && @info "_train_grid_pipeline: artifacts present, loading (skip-if-done)" grid
        return _bundle_from_artifacts(grid, npe_path, ratio_path, ood_path)
    end

    # (2) datagen (real path launches the Philox-per-index simulator; `-t auto` to parallelize).
    #
    # PROVENANCE (07-CALIBRATION-FINDINGS F6): the training image-size distribution is the single
    # most consequential un-recorded training knob — a net trained at one imsize is NOT known to be
    # calibrated at another (F5, covariate shift). We therefore capture the ACTUAL distribution the
    # pool was drawn from and persist it in every artifact's `meta` below. When the caller INJECTS
    # `datagen`, the `imsize_set`/`imsize_weights` keywords do NOT describe the returned pool, so we
    # record them as `:unknown` rather than a guess.
    injected_datagen = datagen !== nothing
    weights = imsize_weights === nothing ?
        Tuple(fill(1.0 / length(imsize_set), length(imsize_set))) : imsize_weights
    if datagen === nothing
        datagen = () -> generate_samples(n_pairs; grid = grid, master_seed = master_seed,
                                         imsize_set = imsize_set, imsize_weights = weights)
    end
    imsize_prov = injected_datagen ?
        (imsize_set = :unknown, imsize_weights = :unknown, imsize_source = :injected_datagen) :
        (imsize_set = Tuple(imsize_set), imsize_weights = Tuple(weights),
         imsize_source = :generate_samples)
    data = datagen()
    θ    = data.theta                                   # 7×N
    Zraw = data.summary_min                             # 2·G²×N (raw, pre-standardization)
    N    = size(Zraw, 2)
    N >= 4 || error("_train_grid_pipeline: need ≥4 training samples, got $N")

    # (3) leak-free train/val split + summary zt (fit on TRAIN continuous rows only).
    nval = clamp(round(Int, val_frac * N), 1, N - 1)
    ntr  = N - nval
    θtr, θva       = θ[:, 1:ntr], θ[:, ntr+1:end]
    Zraw_tr, Zraw_va = Zraw[:, 1:ntr], Zraw[:, ntr+1:end]
    zt      = fit_summary_transform(Zraw_tr; variant = :min)
    Ztr_std = standardize_summary(Zraw_tr, zt, :min)
    Zva_std = standardize_summary(Zraw_va, zt, :min)

    # (4) train NPE + NRE (both use_gpu-plumbed) on the frozen-zt standardized TRAIN pool.
    npe = train_npe(Ztr_std, Zva_std, θtr, θva;
                    use_gpu = use_gpu, epochs = npe_epochs, verbose = verbose, npe_kwargs...)
    ρtr = vec(θtr[1, :])                                 # ρ_true row (D-07 contrast)
    ratio = train_ratio(Ztr_std, ρtr;
                        n = ratio_n, use_gpu = use_gpu, epochs = ratio_epochs,
                        verbose = verbose, ratio_kwargs...)

    # (5) OOD nulls on the TRAIN-ONLY ID pool (density channel; noise/PP composed downstream).
    ood_nulls = (; density = fit_ood_nulls(Ztr_std; variant = :min))

    # (6) persist all three CPU-resident under artifacts_root/grid_G/ (atomic; Pitfall 4).
    save_estimator(npe_path, npe.estimator, npe.θzt, zt, npe.arch;
                   meta = (; grid = Int(grid), n_pairs = N, use_gpu = use_gpu,
                           master_seed = master_seed, imsize_prov...))
    save_ratio(ratio_path, ratio; meta = (; grid = Int(grid), ratio_n = ratio_n, imsize_prov...))
    save_ood_nulls(ood_path, ood_nulls;
                   meta = (; grid = Int(grid), n_train = ntr, imsize_prov...))

    return EstimatorBundle(Int(grid), npe.estimator, ratio, ood_nulls, zt, npe.θzt,
                           _placeholder_calibration(grid))
end
