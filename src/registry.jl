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
# SHIPPED FAMILY = {8} ONLY (07-10, Go/No-Go Option A). The v2.0 GO
# (`07-GO-NO-GO-UPDATE.md`, 2026-07-24 — "GO with named limits") rests entirely on the 8×8
# reference grid: its coloc targets ρ_true/Δρ are SBC-calibrated with randomized-rank atom
# handling, its Bayes factor is validated by simulation-based discrimination (spike 014, AUC
# 0.994), and its OOD flag passes (amended gate, pooled AUC 1.0). D-04 originally scoped a
# (4, 8, 16, 32) family, but 32×32 was CAPPED (never trained — `gate-32x32.md`), 16×16 FAILED its
# gate with no non-method cause (`gate-16x16.md`), and 4×4 was never post-hoc re-analysed with the
# corrected statistics the GO relies on. So only 8×8 ships. 64×64 stays dropped.
#
# `_lazy_load_from_artifact!` (07-10) resolves the shipped 8×8 bundle from a content-hashed,
# tree-sha1-verified `Artifacts.toml` entry (GitHub-Release-hosted, lazy). `_train_grid_pipeline`
# (07-03, src/amortized/pipeline.jl) is the user-definable-grid on-ramp for any OTHER grid.
#############################################################################################

import Pkg
import Artifacts

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
    _SHIPPED_GRIDS :: NTuple{1, Int}

The shipped grid family — **`(8,)` only** (07-10, Go/No-Go Option A, `07-GO-NO-GO-UPDATE.md`).
The v2.0 GO ("GO with named limits") rests on the 8×8 reference grid alone. 4×4 (never post-hoc
re-analysed with the corrected SBC/BF statistics), 16×16 (gate FAILED with no non-method cause,
`gate-16x16.md`), 32×32 (CAPPED — never trained, `gate-32x32.md`) and 64×64 (dropped, D-04) are all
EXCLUDED. A shipped-but-not-yet-loaded grid lazy-loads from its tree-sha1-verified `Artifacts.toml`
entry (`_lazy_load_from_artifact!`); any other grid errors, pointing at `train_and_register`.
"""
const _SHIPPED_GRIDS = (8,)

# Path to the package-root Artifacts.toml (content-hashed, lazy, tree-sha1-verified shipped bundles).
_artifacts_toml() = normpath(joinpath(@__DIR__, "..", "Artifacts.toml"))

# In-repo dev source for a shipped grid's bundle (gitignored; the GitHub-Release asset is built from
# this exact tree). Used as an integrity-verified LOCAL fallback when the artifact is not installed
# in the Julia artifact store and no network download is available (dev/CI). The shipped 8×8 bundle
# is the retrained bounded-θ / realistic-imsize net the GO rests on (`artifacts/amended_v2/grid_8/`).
_dev_bundle_dir(grid::Integer) =
    normpath(joinpath(@__DIR__, "..", "artifacts", "amended_v2", "grid_$(Int(grid))"))

# Verify a resolved bundle directory's git-tree-sha1 against the value pinned in Artifacts.toml
# (T-7-01: never load a tampered/mismatched model). Errors loudly on mismatch.
function _verify_tree_sha1(dir, expected::Base.SHA1, grid::Integer)
    actual = Base.SHA1(Pkg.GitTools.tree_hash(dir))
    actual == expected || error(
        "Artifact integrity check FAILED for the $(grid)×$(grid) estimator bundle: expected " *
        "git-tree-sha1 $(expected) (Artifacts.toml), got $(actual) for $dir. Refusing to load a " *
        "content-mismatched model (T-7-01).")
    return dir
end

"""
    _lazy_load_from_artifact!(grid::Integer) -> EstimatorBundle

Resolve a shipped grid's frozen `EstimatorBundle` from the content-hashed `Artifacts.toml` entry
`grid_G`, VERIFY its git-tree-sha1 (T-7-01), rebuild the bundle (`_bundle_from_artifacts`), wire in
the grid's recorded pre-registered OOD operating point (`_activate_recorded_ood_threshold`), and
`register!` it. Only `8` ships (`_SHIPPED_GRIDS`); any other grid errors.

Resolution order (all three paths integrity-checked before load):
 1. **artifact store** — if the content-addressed artifact is already installed
    (`Artifacts.artifact_exists`), use `Artifacts.artifact_path` (verified at install time);
 2. **in-repo dev fallback** — else if the in-repo bundle dir exists, VERIFY its tree-sha1 against
    Artifacts.toml and load from it (the dev/CI path — no network);
 3. **lazy download** — else `Pkg.Artifacts.ensure_artifact_installed` fetches the GitHub-Release
    asset (production), then resolve + verify via the store.
"""
function _lazy_load_from_artifact!(grid::Integer)
    g = Int(grid)
    g in _SHIPPED_GRIDS || error(
        "_lazy_load_from_artifact!: $(g)×$(g) is not a shipped grid. Shipped: $(_SHIPPED_GRIDS). " *
        "Build one in-process with `train_and_register($g)`.")
    toml = _artifacts_toml()
    name = "grid_$(g)"
    expected = Artifacts.artifact_hash(name, toml)
    expected === nothing && error(
        "_lazy_load_from_artifact!: Artifacts.toml has no `$name` entry at $toml.")

    dir = if Artifacts.artifact_exists(expected)
        Artifacts.artifact_path(expected)                  # content-addressed store (hash-verified)
    else
        devdir = _dev_bundle_dir(g)
        if isdir(devdir)
            _verify_tree_sha1(devdir, expected, g)         # dev/CI fallback, explicit verification
        else
            Pkg.Artifacts.ensure_artifact_installed(name, toml)  # production: fetch the Release asset
            Artifacts.artifact_path(expected)
        end
    end

    npe_path   = joinpath(dir, "npe_$(g).jld2")
    ratio_path = joinpath(dir, "ratio_$(g).jld2")
    ood_path   = joinpath(dir, "ood_nulls_$(g).jld2")
    for p in (npe_path, ratio_path, ood_path)
        isfile(p) || error("_lazy_load_from_artifact!: resolved artifact $dir is missing $p.")
    end
    b = _bundle_from_artifacts(g, npe_path, ratio_path, ood_path)
    # Wire the grid's RECORDED pre-registered OOD operating point (from gate_report_G.jld2 in the
    # same artifact tree) so the shipped OOD flag is operative rather than vacuously false. The
    # persisted ood_nulls carries only the :density channel, so the recorded density id_threshold is
    # scale-comparable to the fused statistic (scale-guard inside `_activate_recorded_ood_threshold`).
    return register!(_activate_recorded_ood_threshold(b, dir))
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

# `_train_grid_pipeline(grid; ...)` — the full per-grid data-gen → NPE + NRE training → OOD fit →
# persist pipeline that returns a populated `EstimatorBundle` — is defined in
# `src/amortized/pipeline.jl` (07-03), included AFTER this file. `train_and_register` references it
# only at call time, so the forward reference resolves once the module finishes loading.

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
