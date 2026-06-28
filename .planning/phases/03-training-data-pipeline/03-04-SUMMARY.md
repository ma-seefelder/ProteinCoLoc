---
phase: 03-training-data-pipeline
plan: 04
subsystem: data-loader
tags: [zscore, kfold, leak-free, standardization, holdout, random123, statsbase]

# Dependency graph
requires:
  - phase: 03-03
    provides: column-major shard_*.jld2 (theta, summary_min[128], summary_aug[AUG_DIM]); separate holdout.jld2 (negative global_index); generate_cache driver; meta.jld2
  - phase: 03-02
    provides: seeding salts (FOLD_SALT, HOLDOUT_SALT) + fold_rng; encode_d01 mask layout (rows 65:128)
provides:
  - "loader.jl: module Loader with load_main_pool / load_fold / load_holdout / all_folds — the SOLE standardization path (D-07/D-08)"
  - "load_fold: ZScoreTransform fit on TRAIN columns only, mask rows (65:128) bypass, Float32 d×K tensors + returned fit object for Phase-5 freeze"
  - "all_folds: deterministic k=5 disjoint/complete fold index vectors from the master seed (D-09)"
  - "SC-3 + SC-4 gates green (no @test_skip placeholders remain)"
affects: [phase-04-training, phase-05-sbc]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Footgun-removal by module boundary: Loader exports only the 4 read APIs; no standardize_all symbol exists anywhere, so cross-fold leakage is impossible by construction (D-08)"
    - "Fit-on-train-only standardization: fit(ZScoreTransform, Z[cont_rows, train_idx]; dims=2), applied to both folds — val statistics never enter the fit (D-07)"
    - "Mask bypass: continuous rows z-scored, binary mask rows (65:128) copied through unchanged so a 0/1 mask never re-couples folds via its mean"
    - "Structural holdout exclusion: load_main_pool globs shard_*.jld2 only — never holdout.jld2 — so the reserved set is outside every fold by construction, not by index arithmetic (D-10)"
    - "module + trailing `using .Loader`: names(Loader) is a controlled D-08 surface while load_fold stays directly callable after a plain include"

key-files:
  created:
    - spike/data/loader.jl
  modified:
    - spike/test/test_data_pipeline.jl

key-decisions:
  - "loader.jl is wrapped in `module Loader` (the ONLY spike module to date) so `names(Loader)` is a meaningful, controlled public surface for the D-08 `:standardize_all ∉ names(Loader)` assertion; a trailing `using .Loader` re-exports the API so the plan's plain-include verify command keeps working"
  - "Row partition derived from variant + nrows: :min ⇒ cont 1:64 / mask 65:128; :aug ⇒ cont [1:64; 129:end] / mask 65:128 (the 14 D-02 moments are all continuous)"
  - "SC-4b holdout exclusion is STRUCTURAL (count==N + hiding holdout.jld2 leaves the pool unchanged + HOLDOUT_SALT≠0), NOT an integer-set intersection of holdout-vs-fold indices (both are 1:N-style ranges that overlap by value and would false-positive)"
  - "A single tiny cache (N=40, master_seed=3, holdout=20) is generated once and shared by SC-3 and SC-4 to keep the quick gate fast"

patterns-established:
  - "Pattern: load_fold returns (Ztr, θtr, Zva, θva, zt) with zt the train-only fit object so Phase 5 freezes the exact preprocessing it cross-validated under"
  - "Pattern: guarded loader include (isdefined(:load_fold) || include) keeps the test idempotent under standalone + runtests loading"

requirements-completed: [DATA-03]

# Metrics
duration: 18min
completed: 2026-06-28
---

# Phase 3 Plan 4: DATA-03 Leak-Free K-Fold Loader Summary

**The sole standardization path: a `module Loader` whose `load_fold` fits a `ZScoreTransform` on TRAIN columns only, bypasses the binary mask, and runs deterministic k=5 CV — with no global-standardize symbol in existence, so cross-fold leakage and holdout contamination are impossible by construction (D-07/D-08/D-09/D-10).**

## Performance

- **Duration:** ~18 min
- **Started:** 2026-06-28
- **Completed:** 2026-06-28
- **Tasks:** 2 (Task 1 implement + automated verify; Task 2 SC-3/SC-4 gates)
- **Files modified:** 2 (1 created, 1 modified)

## Accomplishments
- `loader.jl` — `module Loader` exposing exactly `load_main_pool`, `load_fold`, `load_holdout`, `all_folds`. `load_main_pool` `hcat`s the RAW `shard_*.jld2` columns (never `holdout.jld2`); `load_fold` fits `StatsBase.ZScoreTransform` on the TRAIN columns of a deterministic k=5 fold, applies it to the held-out fold, passes the binary mask rows (65:128) through UNCHANGED, and returns Float32 d×K tensors plus the fit object (so Phase 5 freezes train-only preprocessing).
- D-08 footgun removal is STRUCTURAL: there is no `standardize_all` / global-stats symbol anywhere, and `names(Loader)` is the controlled public surface the SC-3a assertion introspects.
- SC-3 (fit-on-train, mask bypass, no global symbol) and SC-4 (k-fold disjoint/complete + reproducible + structurally-excluded ≥20 holdout) gates are green; the SC-3/SC-4 `@test_skip` placeholders are gone.
- Full suite green (`runtests.jl` exit 0, 68 pipeline tests + simulator tests + resolve-risk gate); decoupling baseline `f581d95 -- src Project.toml Manifest.toml` clean.

## Task Commits

Each task was committed atomically:

1. **Task 1: leak-free k-fold loader (loader.jl)** - `d530dca` (feat)
2. **Task 2: SC-3 + SC-4 leak-free / k-fold / holdout-exclusion gates** - `7c94a38` (test)

**Plan metadata:** (this commit) (docs: complete plan)

## Files Created/Modified
- `spike/data/loader.jl` (created) - `module Loader`: `load_main_pool`/`load_fold`/`load_holdout`/`all_folds` + internal `_shard_files`/`_row_partition`; `using .Loader` re-export
- `spike/test/test_data_pipeline.jl` (modified) - filled SC-3 + SC-4 testsets, added `using Statistics`, guarded loader include, shared `DP_LOAD_DIR` fixture + pre-declared tolerances

## Decisions Made
- **module Loader for a real D-08 surface:** the only spike module so far, justified because `names(Loader)` must be a *controlled* public surface for `:standardize_all ∉ names(Loader)` to be a meaningful structural guarantee rather than introspecting all of `Main`. The trailing `using .Loader` preserves the plain-include callability the plan's verify command relies on.
- **Variant-aware row partition:** `_row_partition(:min, 128)` → cont 1:64 / mask 65:128; `_row_partition(:aug, AUG_DIM)` → cont `[1:64; 129:end]` / mask 65:128. The mask index set is identical for both variants (the 14 augmented moments are all continuous), so the mask is never z-scored under either variant.
- **Structural (not integer-intersection) holdout exclusion:** SC-4b asserts `size(load_main_pool(d).theta,2)==N`, that hiding `holdout.jld2` leaves the pool column count unchanged (proving `load_main_pool` never reads it), and `HOLDOUT_SALT != 0` — deliberately avoiding the value-overlap false-positive an integer intersection of two 1:N-style ranges would produce.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] loader.jl wrapped in `module Loader` (plan implied a "LoaderModule" but the existing spike files are plain includes)**
- **Found during:** Task 1
- **Issue:** The D-08 success criterion and the plan's `!(:standardize_all in names(<LoaderModule>))` require a real module surface to introspect, yet the plan's Task-1 verify command calls `load_fold`/`load_main_pool` *unqualified* after a plain `include` (no `using`). A plain-include-only file gives no controlled `names()` surface; a bare module breaks the unqualified verify call.
- **Fix:** Defined `module Loader` (controlled export surface) and appended `using .Loader` at file end so a plain `include("loader.jl")` both defines `Loader` (for `names(Loader)`) and brings the API into the includer's scope (for unqualified `load_fold`). seeding.jl is included *inside* the module so `FOLD_SALT`/`fold_rng` live on the Loader surface, not leaked.
- **Files modified:** spike/data/loader.jl
- **Verification:** Task-1 verify command (unqualified `load_fold`) green; SC-3a `!(:standardize_all in names(Loader))` green.
- **Committed in:** d530dca

**2. [Rule 3 - Blocking] Shared single cache + Statistics import for the loader gates**
- **Found during:** Task 2
- **Issue:** `using Statistics` (mean/std) cannot live inside an `@testset begin … end` block (it is wrapped in a function); and generating a fresh ≥60-sample cache inside *each* of SC-3 and SC-4 would roughly double the quick-gate generation cost.
- **Fix:** Added `using Statistics` at file top and computed one `const DP_LOAD_DIR` cache shared by both testsets (SC-4 mutates it only by hiding/restoring `holdout.jld2`).
- **Files modified:** spike/test/test_data_pipeline.jl
- **Verification:** Full `test_data_pipeline.jl` green in ~1m (SC-3 ~3.4s, SC-4 ~0.2s); `runtests.jl` exit 0.
- **Committed in:** 7c94a38

---

**Total deviations:** 2 auto-fixed (2 blocking)
**Impact on plan:** Both are structural enablers, not scope changes; all files remain under `spike/`. The TDD RED state for Task 1 was the file's absence (the automated verify command errored pre-implementation) → GREEN after `loader.jl`.

## Threat Model Coverage
- **T-03-10 (cross-fold standardization leakage):** mitigated — `fit(ZScoreTransform, Z[cont_rows, train_idx]; dims=2)` on TRAIN only; no exported global-standardize symbol (SC-3a).
- **T-03-11 (z-scoring the binary mask):** mitigated — mask rows 65:128 bypass standardization, asserted byte-identical to RAW (SC-3b).
- **T-03-12 (holdout contaminating a CV fold):** mitigated — `load_main_pool` never reads `holdout.jld2`; structural exclusion asserted (SC-4b).
- **T-03-13 (non-reproducible folds):** mitigated — `fold_rng(master_seed ⊻ FOLD_SALT, 0)`; two calls give identical membership (SC-4a).

## Known Stubs
None — `load_fold`/`load_main_pool`/`load_holdout`/`all_folds` are fully wired and exercised by green SC-3/SC-4 gates.

## Threat Flags
None — read-only local loader; no new network/auth/file-access surface beyond the planned cache read path already in `<threat_model>`.

## Next Phase Readiness
- DATA-03 complete: Phase 4 (NPE training) can call `load_fold(dir, f; K=5, master_seed)` for honestly cross-validated metrics, and `load_holdout(dir)` for the ADVI benchmark set — with `zt` returned so Phase 5 freezes the exact train-only preprocessing.
- Decoupling holds: `git diff --quiet f581d95 -- src Project.toml Manifest.toml` is clean.

## Self-Check: PASSED
- Files verified present: spike/data/loader.jl, spike/test/test_data_pipeline.jl
- Commits verified: d530dca, 7c94a38
- Gates: `julia --project=spike spike/test/test_data_pipeline.jl` green (SC-3 + SC-4 no longer skipped, 68/68); `runtests.jl` exit 0; decoupling baseline clean

---
*Phase: 03-training-data-pipeline*
*Completed: 2026-06-28*
