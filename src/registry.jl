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
# src/registry.jl --- grid-keyed estimator registry (PROD-02, D-04).
#
# The public `num_patches` value is the registry KEY: each shipped grid has its own calibrated
# `EstimatorBundle` (NPE + NRE + OOD nulls + frozen standardizers + its passed D-05 gate report).
# `estimator_for(grid)` returns a registered bundle, lazily loads a shipped grid's artifact, or —
# for any other grid — throws a clear `ArgumentError` that names the shipped grids and points at
# `train_and_register` (T-7-05: validate the grid, NEVER silently fall back to a default grid).
#
# D-04 caps the shipped family to (4, 8, 16, 32) — 64 is DROPPED (finer grids lengthen the input
# but the NPE output stays a single GLOBAL ρ; true per-region maps are Phase 12). `_train_grid_
# pipeline` and `_lazy_load_from_artifact!` are HOOKS filled by later Phase-7 plans (per-grid
# training + Artifacts wiring); here they are honest stubs that error until populated.
#############################################################################################

"""
    EstimatorBundle

Everything the amortized read path needs for ONE patch grid, keyed by `grid` in `_REGISTRY`.

# Fields
- `grid::Int`: patch grid `G` (the registry key).
- `npe`: trained NPE `PosteriorEstimator` (posterior draws).
- `ratio`: trained NRE `RatioEstimator` (amortized log Bayes factor).
- `ood_nulls`: fitted OOD null statistics (density + noise channels).
- `zt`: frozen summary `ZScoreTransform` (applied at inference; mask rows bypassed).
- `θzt`: frozen parameter `ZScoreTransform` (un-standardizes posterior draws before any ρ read).
- `calibration::CalibrationMeta`: the grid's PASSED D-05 ship-gate report.
"""
struct EstimatorBundle
    grid        :: Int
    npe
    ratio
    ood_nulls
    zt
    θzt
    calibration :: CalibrationMeta
end

"""
    _REGISTRY :: Dict{Int, EstimatorBundle}

Process-local map from patch grid `G` to its registered `EstimatorBundle`. Populated by
`register!` (either from `train_and_register` or a lazy artifact load).
"""
const _REGISTRY = Dict{Int, EstimatorBundle}()

"""
    _SHIPPED_GRIDS :: NTuple{4, Int}

The CAPPED shipped grid family `(4, 8, 16, 32)` (D-04). 64 is DROPPED. Only grids that PASS their
independent D-05 ship-gate are populated into `_REGISTRY`; a shipped-but-not-yet-loaded grid lazy-
loads from its artifact, an unshipped grid errors.
"""
const _SHIPPED_GRIDS = (4, 8, 16, 32)

"""
    _lazy_load_from_artifact!(grid::Integer) -> EstimatorBundle

HOOK: load a shipped grid's frozen `EstimatorBundle` from its (Artifacts-backed) store and
`register!` it. Wired by a later Phase-7 plan (Artifacts + `Flux.state`/`loadmodel!` persistence);
until then it errors so a shipped-but-unpopulated grid fails loudly rather than silently.
"""
function _lazy_load_from_artifact!(grid::Integer)
    error("No artifact populated for a $(grid)×$(grid) estimator yet. Shipped grids " *
          "$(_SHIPPED_GRIDS) gain artifacts in a later Phase-7 plan; until then build one " *
          "in-process with `train_and_register($grid)`.")
end

"""
    estimator_for(grid::Integer) -> EstimatorBundle

Return the `EstimatorBundle` for a `grid×grid` patch grid: a registered bundle if present, else a
lazily loaded shipped-grid artifact, else a clear `ArgumentError` naming the shipped grids and
`train_and_register` (T-7-05 — never a silent default-grid fallback).
"""
function estimator_for(grid::Integer)
    haskey(_REGISTRY, grid) && return _REGISTRY[grid]
    grid in _SHIPPED_GRIDS && return _lazy_load_from_artifact!(grid)
    throw(ArgumentError(
        "No estimator registered for a $(grid)×$(grid) patch grid. " *
        "Shipped grids: $(_SHIPPED_GRIDS). " *
        "Train and register one with `train_and_register($grid)`."))
end

"""
    register!(b::EstimatorBundle) -> EstimatorBundle

Register `b` under its `grid`, so subsequent `estimator_for(b.grid)` returns it. Returns `b`.
"""
function register!(b::EstimatorBundle)
    _REGISTRY[b.grid] = b
    return b
end

"""
    has_cuda_device() -> Bool

Weakdep-safe GPU probe: `true` iff a functional CUDA device is available. CUDA is an OPTIONAL
weakdep (D-06), so this never `import`s CUDA in core — it inspects the already-loaded modules and
calls `CUDA.functional()` only if CUDA has been loaded. On a CPU-only process (no CUDA loaded) it
returns `false` without erroring, so the shipped/reproducible path degrades gracefully to CPU.
"""
function has_cuda_device()
    for (id, m) in Base.loaded_modules
        if id.name == "CUDA"
            try
                return Bool(Base.invokelatest(getproperty(m, :functional)))
            catch
                return false
            end
        end
    end
    return false
end

"""
    _train_grid_pipeline(grid::Integer; use_gpu, kwargs...) -> EstimatorBundle

HOOK: run the full per-grid data-gen → NPE + NRE training → OOD fit → D-05 gate pipeline and
return a calibrated `EstimatorBundle`. Wired by the per-grid Phase-7 plans (07-03+); until then it
errors so `train_and_register` fails loudly rather than returning an unpopulated bundle.
"""
function _train_grid_pipeline(grid::Integer; use_gpu = has_cuda_device(), kwargs...)
    error("Per-grid training pipeline for a $(grid)×$(grid) grid is not implemented yet " *
          "(added by a later Phase-7 plan). use_gpu=$(use_gpu).")
end

"""
    train_and_register(grid::Integer; use_gpu = has_cuda_device(), kwargs...) -> EstimatorBundle

Train an estimator for `grid` (GPU-accelerated when available, D-06; graceful CPU fallback) via
`_train_grid_pipeline`, `register!` it, and return the bundle. The documented on-ramp for any grid
outside the shipped family (PROD-02, D-04).
"""
function train_and_register(grid::Integer; use_gpu = has_cuda_device(), kwargs...)
    b = _train_grid_pipeline(grid; use_gpu = use_gpu, kwargs...)
    return register!(b)
end
