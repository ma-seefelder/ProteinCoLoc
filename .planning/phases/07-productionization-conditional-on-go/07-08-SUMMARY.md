---
phase: 07-productionization-conditional-on-go
plan: 08
subsystem: amortized-inference
tags: [windowed-inference, local-map, registry, estimator-bundle, sub-tile, ood, phase-12-extension-point]

# Dependency graph
requires:
  - phase: 07-05
    provides: the gated 8×8 bundle on disk (artifacts/grid_8/{npe,ratio,ood_nulls}_8.jld2)
  - phase: 07-03
    provides: persist.jl load_estimator/load_ratio/load_ood_nulls + pipeline.jl _bundle_from_artifacts/_grid_dir
  - phase: 07-02
    provides: the amortized read surface (standardize_summary, delta_rho, ood_verdict)
provides:
  - "_ensure_grid_registered(grid; artifact_dir) — in-process registration bridge from local artifacts (wave-6, no Artifacts.toml)"
  - "local_coloc_map(img, control, channels; grid, tiles, N) -> LocalColocMap — coarse windowed sub-tile Δρ + OOD map"
  - "test/test_local_map.jl (wired into runtests.jl), incl. a conditional real-8×8 smoke"
affects: [07-10, phase-12-spatial-map]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "In-process registry population from a local artifact dir, reusing _bundle_from_artifacts (no restatement of the load path)"
    - "Degenerate-tile finite sentinel + mandatory OOD flag instead of an exception (T-7-04)"
    - "Conditional test block for gitignored artifacts: real bundle exercised where present, honest @info skip otherwise"

key-files:
  created:
    - src/amortized/local_map.jl
    - test/test_local_map.jl
  modified:
    - src/ProteinCoLoc.jl
    - test/runtests.jl

key-decisions:
  - "local_coloc_map takes `grid` as a keyword defaulting to 8 (the gated reference grid) instead of hardcoding 8 — the map inherits whichever grid's gate the caller selects, and the tests can drive a tiny grid-4 bundle"
  - "Sub-tiles are cut with the package's OWN patch(img, nx, ny) tiler, so a sub-tile is exactly the block the summary path would see (same trimming rule)"
  - "Tile MCIs carry the PARENT image's Otsu thresholds — a tile-local Otsu would re-threshold each window and make tiles incomparable"
  - "LOCAL_MAP_SENTINEL = 0.0 with a MANDATORY ood_flag=true, so a sentinel is never indistinguishable from a measured Δρ of zero"
  - "LocalColocMap is deliberately NOT an AbstractColocResult subtype — it has no posterior draws, no Bayes factor and no per-region uncertainty, so implementing that interface would overstate it"

patterns-established:
  - "Wave-ordering bridge: when a later plan owns the public load path, the earlier consumer loads in-process from in-repo outputs and says so in the docstring, rather than pre-empting the later design"

requirements-completed: [PROD-02]

# Metrics
duration: 55min
completed: 2026-07-21
---

# Phase 7 Plan 08: Windowed Sub-Tile Local Colocalization Map

**The gated 8×8 estimator now runs over image sub-tiles to produce a coarse per-tile Δρ + OOD map, registering the frozen bundle in-process from `artifacts/grid_8/` (the Artifacts lazy path is 07-10's job) — no retraining, no new gate, and the Phase-12 per-region extension point left untouched.**

## Performance

- **Duration:** ~55 min
- **Tasks:** 2 (both executed)
- **Files:** 2 created, 2 modified
- **Verification:** `julia --project -e 'using Pkg; Pkg.test()'` → **`Testing ProteinCoLoc tests passed`** (full suite, no failures, no errors)

## Accomplishments

- **`_ensure_grid_registered(grid; artifact_dir)`** — loads `npe_G.jld2` / `ratio_G.jld2` / `ood_nulls_G.jld2` from a local artifact directory via the existing `_bundle_from_artifacts` (which already composes `load_estimator`/`load_ratio`/`load_ood_nulls`) and `register!`s the result. Verified on the real 8×8 bundle: `julia --project -e 'using ProteinCoLoc; ProteinCoLoc._ensure_grid_registered(8); @assert haskey(ProteinCoLoc._REGISTRY, 8)'` → `OK grid=8 zt=StatsBase.ZScoreTransform{Float64, Vector{Float64}}`. It references neither `_lazy_load_from_artifact!` nor Artifacts in executable code (machine-asserted after stripping docstrings/comments).
- **`local_coloc_map(img, control, channels; grid = 8, tiles = (4,4), N = 2000)`** — cuts both images into `r×c` sub-tiles with the package's own `patch(img, nx, ny)`, builds a 2-channel tile `MultiChannelImage` carrying the parent Otsu thresholds, and runs the frozen bundle per tile: `patch_summary(tile, G)` → `encode_d01` → `standardize_summary(…, b.zt)` → `delta_rho(b.npe, Zs, Zc, b.θzt)` → `ood_verdict`. Returns a lightweight `LocalColocMap(grid, tiles, delta_rho::Matrix, ood_flag::Matrix, meta)`.
- **Three finite-guard paths (T-7-04)** — a sub-tile smaller than the patch grid, a tile where no patch survived the ≥15-survivor floor, and a non-finite posterior read each write `LOCAL_MAP_SENTINEL = 0.0` **and** set `ood_flag = true`. The map never throws on a degenerate layout; `meta.n_sentinel` counts them.
- **Honest scope framing in the docstring** — states plainly that each tile is an independent global-ρ read on a smaller window, with no spatial prior, no borrowing of strength between neighbours, and no per-region uncertainty, and that the calibrated per-region GP/CAR map (`SpatialColocResult` / `delta_rho_map` / `uncertainty_map`) is Phase 12 and not implemented here.
- **Tests** — `test/test_local_map.jl` (30 assertions, wired into `runtests.jl`) drives a tiny grid-4 bundle written into a temp dir, plus a **conditional real-8×8 smoke** (5/5 on this checkout) that runs 2×2 windowed inference over the actual gated bundle with zero sentinels.

## Task Commits

1. **Task 1 + Task 2 (source):** `25037e6` — `feat(07-08): windowed sub-tile local coloc map + in-process 8x8 registration bridge` (`src/amortized/local_map.jl`, `src/ProteinCoLoc.jl`)
2. **Task 2 (tests):** `90f612e` — `test(07-08): tile-grid + per-tile finiteness tests for the windowed local map` (`test/test_local_map.jl`, `test/runtests.jl`)

## Files Created/Modified

- `src/amortized/local_map.jl` — **created**: `LOCAL_MAP_SENTINEL`, `LocalColocMap`, `_ensure_grid_registered`, `_sub_tiles`, `_tile_mci`, `_tile_summary`, `local_coloc_map`.
- `src/ProteinCoLoc.jl` — **modified**: `include("amortized/local_map.jl")` after `pipeline.jl` (it composes the read surface *and* the pipeline's artifact loaders), and `export local_coloc_map, LocalColocMap`.
- `test/test_local_map.jl` — **created**: bridge + map contract + Phase-12-cleanliness assertions, and the conditional real-8×8 smoke.
- `test/runtests.jl` — **modified**: `include(joinpath(@__DIR__, "test_local_map.jl"))` after the GPU smoke.

## Decisions Made

- **`grid` is a keyword (default 8), not a hardcode.** The plan's behavior spec names the 8×8 bundle; making the grid a defaulted keyword keeps that as the shipped behaviour while letting the CI test drive a 2-second grid-4 bundle instead of depending on a gitignored 5 MB artifact. Default call `local_coloc_map(img, ctrl, [1,2])` is exactly the 8×8 path.
- **Reuse `_bundle_from_artifacts` rather than restating the three loads.** It already performs `load_estimator`/`load_ratio`/`load_ood_nulls` and assembles the `EstimatorBundle` with the placeholder calibration; duplicating that in the bridge would create a second, drift-prone load path.
- **Parent Otsu thresholds on tile MCIs.** Recomputing Otsu per tile would silently re-threshold every window and make tiles incomparable across the map.
- **Sentinel = `0.0` + forced flag.** A finite sentinel keeps `all(isfinite, delta_rho)` a real invariant (no `NaN` leaking into downstream plotting/aggregation), and the forced OOD flag prevents "unscorable" from reading as "measured no difference".
- **`LocalColocMap` is not an `AbstractColocResult`.** Implementing `delta_rho`/`bayes_factor`/`is_ood`/`posterior_draws` would claim a calibrated single-result contract the windowed readout does not have.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Repo reality] Task 1 and Task 2 share one new file, committed source-then-tests**
- **Found during:** commit step
- **Issue:** The plan assigns `src/amortized/local_map.jl` to BOTH tasks. The bridge and the map function were authored in one pass, so a strictly per-task commit split of that single file was not meaningful.
- **Fix:** Two atomic commits along a real seam instead: `25037e6` = the complete source (+ module wiring), `90f612e` = the tests (+ `runtests.jl` wiring). Task 1's own acceptance criterion was verified independently before the first commit via the plan's `<automated>` command.
- **Verification:** `_ensure_grid_registered(8)` populates `_REGISTRY[8]` from the real artifacts (command output recorded above); full `Pkg.test()` green after the second commit.
- **Committed in:** `25037e6`, `90f612e`

**2. [Rule 3 - Blocking] Test assertions had to strip docstrings before grepping for Phase-12 symbols**
- **Found during:** Task 2 verification (first `Pkg.test()` run, 4 failures)
- **Issue:** The plan's acceptance criterion is `grep -v '^#'` finding no `SpatialColocResult` / `delta_rho_map` **definition**. A line-comment filter alone does not strip Julia `"""…"""` docstrings, and the docstrings deliberately NAME those Phase-12 symbols to explain what the function is not — so the naive check failed on its own scope-honesty prose.
- **Fix:** Added `_code_only(path)` to the test, which drops `"""…"""` docstrings, `#=…=#` blocks and `#` lines, and asserts against the executable code only (plus definition-shaped regexes `delta_rho_map\s*\(`). The source was NOT weakened to satisfy the literal grep — the honest docstring stayed.
- **Files modified:** `test/test_local_map.jl`
- **Committed in:** `90f612e`

**3. [Rule 3 - Blocking] Both test blocks must clear the `_REGISTRY` slot they exercise**
- **Found during:** Task 2 verification (2 further failures)
- **Issue:** Earlier testsets in `runtests.jl` register SYNTHETIC grid-4 and grid-8 bundles. `_ensure_grid_registered` is (correctly) a no-op on an already-registered grid, so the grid-4 pre-condition assertion failed and the real-8×8 smoke silently ran against a stale synthetic bundle with `zt === nothing` (`MethodError: no method matching transform(::Nothing, ::Matrix{Float64})`).
- **Fix:** Both blocks now `delete!` the slot first and restore the prior occupant afterwards. This surfaced a genuine property worth stating: **`_ensure_grid_registered` never overwrites an existing registration** — whatever registered the grid first wins, which is the intended no-op semantics but does mean a caller wanting the artifact bundle specifically must clear the slot.
- **Files modified:** `test/test_local_map.jl`
- **Committed in:** `90f612e`

---

**Total deviations:** 3 auto-fixed (all Rule 3). No pre-registered gate constants, no `test/gate/gate_consts_*.jl`, and nothing under `spike/` was touched; no estimator was retrained and no gate was re-run.
**Impact on plan:** None on scope. The commit seam differs from the literal per-task split, and the Phase-12-cleanliness check is enforced on executable code rather than raw file text.

## Issues Encountered

- **`_ensure_grid_registered` is no-op-on-registered by design, which can mask a stale bundle.** Discovered via the test failure above. It is the right semantics for a registry (an explicitly registered bundle must win over a disk artifact), but it means "ensure" does not mean "ensure *this* artifact". Documented in the docstring; **07-10 should be aware** when it wires the Artifacts path — if a lazily loaded artifact must supersede a registration, that needs an explicit invalidation, not a second `_ensure_*` call.
- **`meta.ood_score` is `NaN` for sentinel tiles and the flags are almost always driven by the sentinel path, not by a fired detector.** The persisted `ood_nulls` from `_train_grid_pipeline` carry only the `:density` channel — no `:noise`, no `:model`, and **no `:thr`** — so `ood_verdict` computes a Mahalanobis score but has no threshold to compare it against and returns `flag = false`. Per-tile OOD flagging is therefore currently *structural only* (below-floor / degenerate tiles). This is inherited from the artifact contract, not introduced here; a plan that wants live per-tile OOD firing must persist the fitted `thr`/`zref` alongside the nulls.

## Next Phase Readiness

- `local_coloc_map` is callable today against the on-disk 8×8 bundle and inherits the 07-05 gate verbatim — **including that gate's honest FAIL residuals**; the map is only ever as calibrated as the grid it rides.
- **07-10 (Artifacts wiring, wave 8)** can replace the bridge's file loads with the content-hashed path; the call site (`_ensure_grid_registered(G)` → `estimator_for(G)`) does not change. Note the no-op-on-registered semantics above.
- **Phase 12** starts from a clean slate: `SpatialColocResult` / `delta_rho_map` / `uncertainty_map` remain sketch-only comments in `src/results.jl`, machine-asserted absent from `local_map.jl`'s executable code.

## Self-Check: PASSED

- `src/amortized/local_map.jl` contains `function local_coloc_map` and `function _ensure_grid_registered` — FOUND, both tracked.
- `test/test_local_map.jl` references `local_coloc_map` and is included from `test/runtests.jl` — FOUND, both tracked.
- Commits `25037e6`, `90f612e` — FOUND.
- `Pkg.test()` → `Testing ProteinCoLoc tests passed` (windowed map 30/30, real-8×8 smoke 5/5).
- `spike/` untouched; no `test/gate/gate_consts_*.jl` modified; no training or gate run.

---
*Phase: 07-productionization-conditional-on-go*
*Completed: 2026-07-21*
