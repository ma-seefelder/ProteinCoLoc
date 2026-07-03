---
phase: 07-productionization-conditional-on-go
plan: 01
subsystem: amortized grid-coupling + estimator registry
status: COMPLETE
tags: [grid-parametric, summary-encoder, dimension-helpers, data-generation, content-hash-cache, estimator-registry, PROD-02, D-04]
requirements: [PROD-02]
dependency_graph:
  requires:
    - 07-00 (AbstractColocResult/CalibrationMeta hierarchy; green root Pkg.resolve; NeuralEstimators 0.2.1 pinned)
  provides:
    - grid-parametric patch_summary(mci,G) + encode_d01 + _summary_row_partition(:min,2G^2)
    - single-source dimension helpers summary_dim(G)=2G^2 / cont_rows(G)=G^2 / ratio_input_dim(G)=5G^2
    - grid-parametric data generation + atomic content-hashed cache (src/amortized/datagen.jl)
    - grid-keyed estimator registry (EstimatorBundle, _REGISTRY, _SHIPPED_GRIDS=(4,8,16,32), estimator_for, register!, train_and_register, has_cuda_device)
  affects:
    - all per-grid Phase-7 plans (instantiate the SAME parametrized encoder per grid; register into _REGISTRY)
    - the amortized simulator/training/inference/OOD promotions (derive every dimension from the helpers here)
tech_stack:
  added: []
  patterns: [single-source grid coupling, Base.hash content-addressed cache (no new stdlib dep), atomic .tmp+reopen-integrity+mv persistence, weakdep-safe GPU probe, hook-stub registry skeleton]
key_files:
  created:
    - src/amortized/summary.jl
    - src/amortized/datagen.jl
    - src/registry.jl
    - .planning/phases/07-productionization-conditional-on-go/07-01-SUMMARY.md
  modified:
    - src/ProteinCoLoc.jl
    - test/runtests.jl
decisions:
  - "D-04 grid coupling centralized ONCE: summary_dim/cont_rows/ratio_input_dim are the single source of truth; every downstream module derives dims from them, no per-grid re-hardcoding"
  - "_SHIPPED_GRIDS=(4,8,16,32) — 64 DROPPED per D-04 feasibility verdict"
  - ":min-only shipped path — the spike :aug moment-superset (encode_aug/AUG_DIM) dropped from the productionized datagen (ABL-01 convenience, not shipped)"
  - "content-hash cache uses Base.hash (not the SHA stdlib) to avoid re-resolving the fragile Wave-0 co-resolution Manifest; grid enters the hash via summary_min_dim=summary_dim(G) so grids auto-separate"
  - "registry validates grid against _REGISTRY/_SHIPPED_GRIDS and throws a train_and_register-pointing ArgumentError — never a silent default-grid fallback (T-7-05)"
metrics:
  tasks_completed: 3
  tasks_total: 3
  files_created: 4
  files_modified: 2
  completed_date: 2026-07-03
---

# Phase 7 Plan 01: Grid Coupling Generalization + Estimator Registry Skeleton Summary

**One-liner:** Generalized every grid-coupled site from the hardcoded 8×8/128/64 spike
constants to an arbitrary grid `G` in one place — a grid-parametric summary/encoder with
single-source dimension helpers (`summary_dim`/`cont_rows`/`ratio_input_dim`), a grid-parametric
data-gen + content-hashed cache, and a grid-keyed estimator registry (PROD-02) that returns a
bundle for a shipped grid `(4,8,16,32)` or raises a clear `train_and_register` error for anything
else.

## Status: COMPLETE

All three tasks complete, each committed atomically (two under a TDD RED→GREEN gate). `Pkg.test()`
is fully green (co-resolution 4/4, D-02 10/10, LoadImages 32/32, patch 3/3, Colocalization 29/29,
**amortized summary 20/20**, **estimator registry 7/7**). `spike/` provably byte-untouched.

## What Was Built

### Task 1 — Grid-parametric summary/encoder + dimension helpers (`src/amortized/summary.jl`)
Commits: `5d5da79` (test/RED), `54bad6f` (feat/GREEN)
- `patch_summary(mci, G)` — the `G×G` per-patch Pearson correlation via the UNCHANGED
  `patch()`/`correlation()` (already grid-general; the ≥15-survivor floor is preserved).
- `encode_d01(M)` — body unchanged (column-major `vec` handles any `G`); length `2·G²`.
- `_summary_row_partition(:min, nrows)` — replaces the spike's hardcoded 128-row guard with an
  `iseven(nrows)` check; returns `(1:G², G²+1:2G²)` for every `G`.
- Single-source dimension helpers: `summary_dim(G)=2G²`, `cont_rows(G)=G²`,
  `ratio_input_dim(G)=5G²`. `G=8` reproduces the proven spike constants (128 / 64 / 320).

### Task 2 — Grid-parametric data generation + atomic cache (`src/amortized/datagen.jl`)
Commit: `3433f4a` (feat)
- Promoted the spike `seeding.jl` + `generate.jl` + `cache.jl` into one grid-threaded module:
  Philox4x keyed RNGs (`sample_rng`/`holdout_rng`/`fold_rng` with disjoint XOR-salts), cost-aware
  `sample_imsize`, the in-memory column-major generation core, and the sharded cache driver.
- **Grid delta:** every summary buffer is `Matrix(summary_dim(grid), N)` (not literal 128), and
  `generating_config.summary_min_dim = summary_dim(grid)`, so a 4×4 (dim 32) and a 16×16 (dim 512)
  cache resolve to DIFFERENT content-hash dirs (verified: hashes differ).
- Exposed `imsize_set`/`imsize_weights` keywords so per-grid plans can bias image size upward
  (fine grids need bigger images to clear the ≥15-survivor floor — Pitfall 2 / T-7-03).
- Kept the atomic `.tmp` → reopen-integrity `@assert` → `mv(...; force=true)` idiom verbatim in
  `write_shard`/`write_meta`/`_write_holdout`.

### Task 3 — Grid-keyed estimator registry skeleton (`src/registry.jl`)
Commits: `bf18583` (test/RED), `694a1de` (feat/GREEN)
- `struct EstimatorBundle` (`grid`, `npe`, `ratio`, `ood_nulls`, `zt`, `θzt`,
  `calibration::CalibrationMeta`); `const _REGISTRY`; `const _SHIPPED_GRIDS = (4, 8, 16, 32)`
  (64 dropped, D-04).
- `estimator_for(grid)` — registered bundle → lazy shipped-artifact load → clear `ArgumentError`
  that names the shipped grids and points at `train_and_register` (T-7-05: no silent default).
- `register!`, `train_and_register(grid; use_gpu = has_cuda_device(), kwargs...)`, and a
  weakdep-safe `has_cuda_device()` (inspects `Base.loaded_modules`; returns `false` with no CUDA
  loaded, never errors — D-06 graceful CPU fallback).
- `_lazy_load_from_artifact!` and `_train_grid_pipeline` are honest error-stubs (wired by later
  Phase-7 plans: Artifacts persistence + per-grid training).
- Exported `estimator_for`, `register!`, `train_and_register`.

## Deviations from Plan

### [Scope adaptation] `:aug` moment-superset dropped from the promoted datagen
- **Found during:** Task 2.
- **Issue:** The spike `generate.jl`/`cache.jl` also carry the `:aug` encoding
  (`encode_aug`/`summary_aug`/`AUG_DIM`/`N_AUG_MOMENTS`), an ABL-01 feature-slicing convenience.
  Task 1's `summary.jl` was scoped to promote `encode_d01` only (not `encode_aug`), and the
  shipped D-04 estimators consume the `:min` encoding exclusively.
- **Fix:** Promoted the `:min`-only path (`summary_min` buffers, `_summary_row_partition(:min,…)`);
  dropped the aug superset. Documented in the datagen header. No shipped-path capability lost.

### [Rule 3 - Blocking] Content hash uses `Base.hash` instead of the `SHA` stdlib
- **Found during:** Task 2.
- **Issue:** The spike `hashguard.jl` uses `using SHA` (a stdlib). `SHA` is NOT in the root
  `Project.toml` `[deps]`; `using SHA` would fail precompile, and adding it would force a re-resolve
  of the fragile Wave-0 co-resolution Manifest (07-00's hard-won green gate).
- **Fix:** Implemented a self-contained `cache_hash`/`subhashes` over `Base.hash` (folds file bytes
  + `_canonical(config)`), zero new dependencies. The grid still separates caches via
  `summary_min_dim = summary_dim(G)` in the hashed config (verified: 4×4 vs 16×16 hashes differ).
  The atomic-write wrapper is unchanged.

### [Design] Simulator chain referenced but not yet promoted
- `datagen.jl`'s `generate_sample`/`_write_holdout` reference `sample_prior`/`simulate_pair`/
  `build_mci` — the Phase-2 forward model, promoted into `src/` by a LATER Phase-7 plan. They are
  referenced only inside function bodies, so the module LOADS today; the datagen becomes runnable
  once the simulator is promoted. Documented in the datagen header. Task-2 verification
  (`summary_dim(16)==512` + spike-clean) does not exercise generation, so this is non-blocking.

## Authentication Gates
None.

## Known Stubs
By design (registry hooks filled by later Phase-7 plans, D-04):
- `_lazy_load_from_artifact!(grid)` — errors "No artifact populated…"; Artifacts persistence is a
  later plan. Shipped-but-unloaded grids fail loudly rather than silently.
- `_train_grid_pipeline(grid; …)` — errors "…not implemented yet"; the per-grid data-gen → NPE +
  NRE → OOD → D-05 gate pipeline is added by 07-03+.
- `datagen.jl` simulator references (`sample_prior`/`simulate_pair`/`build_mci`) are undefined
  until the simulator is promoted (see Deviations). These are the plan's intended on-ramp seams,
  not accidental gaps — the grid-coupling and registry deliverables of THIS plan are complete.

## Threat Flags
None new. The plan's `<threat_model>` mitigations are implemented: T-7-05 (validate grid, no
silent fallback) via `estimator_for`'s `ArgumentError`; T-7-03 (fine-grid missingness) via per-grid
content-hash cache separation + the exposed `imsize_set` bias knob.

## Deferred Items
- Runnable end-to-end data generation awaits the simulator promotion (later Phase-7 plan).
- Artifacts-backed `Flux.state`/`loadmodel!` estimator persistence (the `_lazy_load_from_artifact!`
  body) is a later plan.

## Self-Check

Created files:
- FOUND: src/amortized/summary.jl
- FOUND: src/amortized/datagen.jl
- FOUND: src/registry.jl
- FOUND: .planning/phases/07-productionization-conditional-on-go/07-01-SUMMARY.md

Commits:
- FOUND: 5d5da79 (Task 1 test/RED)
- FOUND: 54bad6f (Task 1 feat/GREEN)
- FOUND: 3433f4a (Task 2 feat)
- FOUND: bf18583 (Task 3 test/RED)
- FOUND: 694a1de (Task 3 feat/GREEN)

Verifications:
- `Pkg.test()`: fully green (amortized summary 20/20, estimator registry 7/7, all prior suites unchanged).
- `summary_dim(16)==512`; per-grid `generating_config` dims 4×4→32, 16×16→512; content hashes differ.
- `estimator_for(7)` throws `ArgumentError` naming `train_and_register`; `has_cuda_device()==false`.
- spike/Project.toml + spike/Manifest.toml + all spike/ sources: provably UNTOUCHED (`git diff --quiet spike/`).

## Self-Check: PASSED
