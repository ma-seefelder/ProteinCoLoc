---
phase: 07-productionization-conditional-on-go
plan: 03
subsystem: amortized training layer (NPE/NRE training, CPU-resident persistence, per-grid pipeline)
status: COMPLETE
tags: [amortized-training, npe-architecture, nre-ratio, gpu-plumbed-d06, flux-state-persistence, cpu-reproducible, train-and-register, idempotent-pipeline, grid-general]
requirements: [PROD-01, PROD-02]
dependency_graph:
  requires:
    - 07-00 (co-resolution GREEN: NeuralEstimators 0.2.1 / Flux 0.16.10; D-02 hierarchy: EstimatorBundle host types CalibrationMeta/OODVerdict)
    - 07-01 (grid-parametric summary.jl helpers summary_dim/cont_rows/ratio_input_dim + _summary_row_partition; datagen.jl generate_samples; registry.jl EstimatorBundle/_REGISTRY/_SHIPPED_GRIDS/has_cuda_device/train_and_register)
    - 07-02 (read surfaces: standardize_summary/posterior_for (infer.jl), amortized_log_bf/pair_encode (bf.jl), fit_ood_nulls/maha_score/noise_features (ood.jl))
  provides:
    - "src/amortized/architecture.jl: build_estimator (input-width-agnostic; q a NormalisingFlow INSTANCE positional) + NPE_* Phase-5 constants"
    - "src/amortized/train_npe.jl: train_npe (GPU-plumbed, no throw guard), fit_theta_transform, fit_summary_transform (leak-free zt)"
    - "src/amortized/train_ratio.jl: build_ratio_estimator (input width ratio_input_dim(G)=5G²), assemble_ratio_pairs, measure_log_prior_odds, train_ratio (GPU-plumbed, no custom loss)"
    - "src/amortized/persist.jl: save/load_estimator + save/load_ratio (CPU-resident Flux.state + loadmodel!, Pitfall 4/T-7-07) + save/load_ood_nulls, all atomic .tmp→integrity→mv; _estimator_ok/_ratio_ok/_ood_nulls_ok integrity predicates"
    - "src/amortized/pipeline.jl: _train_grid_pipeline(grid; ...) → EstimatorBundle (SKIP-IF-DONE + atomic persist); default_imsize_for/default_npairs_for; train_and_register(grid) wired to run it"
  affects:
    - "the per-grid D-05 ship-gate (07-05+) trains each shipped grid THROUGH _train_grid_pipeline and fills the placeholder calibration"
    - "the user-definable-grid on-ramp (D-04): train_and_register(6) now runs REAL training code, not a stub"
    - "colocalization_amortized (later plan) consumes a populated EstimatorBundle (npe + ratio handle + ood_nulls density + frozen zt/θzt)"
tech_stack:
  added: []
  patterns: [D-06 use_gpu split (train GPU-plumbed / persist+infer CPU-resident), Flux.state+loadmodel! CPU-resident persistence (Pitfall 4/T-7-07), atomic .tmp→reopen-integrity→mv wrapper reused verbatim, SKIP-IF-DONE idempotency via integrity predicates, injectable datagen seam (simulator-independent testability), grid coupling via single-source ratio_input_dim/_summary_row_partition, AdamW-Float64 CosAnneal gotcha preserved]
key_files:
  created:
    - src/amortized/architecture.jl
    - src/amortized/train_npe.jl
    - src/amortized/train_ratio.jl
    - src/amortized/persist.jl
    - src/amortized/pipeline.jl
    - .planning/phases/07-productionization-conditional-on-go/07-03-SUMMARY.md
  modified:
    - src/ProteinCoLoc.jl
    - src/registry.jl
    - test/runtests.jl
decisions:
  - "train_npe/train_ratio default use_gpu=has_cuda_device() (D-06); the spike's use_gpu && throw guards are REMOVED on the training path; LR/weight_decay stay Float64 (AdamW-CosAnneal gotcha)"
  - "Shipped models persist as Flux.state(cpu(est)) + arch metadata (NOT whole PosteriorEstimator); load rebuilds via build_estimator/build_ratio_estimator then Flux.loadmodel! — CPU-resident, device-independent (Pitfall 4 / T-7-07 / T-7-01)"
  - "The EstimatorBundle.ratio field carries the ratio HANDLE NamedTuple (estimator, log_prior_odds, num_summaries, summary_width, input_dim) — the registry struct has no log_prior_odds field, so the handle keeps the measured prior-odds with its estimator for amortized_log_bf"
  - "EstimatorBundle.ood_nulls persists the frozen density null only (; density=fit_ood_nulls(...)); the PP :model handle + :noise null are composed by the inference/gate layer, not baked into the persisted artifact"
  - "_train_grid_pipeline gained an injectable `datagen` keyword (defaults to generate_samples) so the pipeline is testable end-to-end NOW without the not-yet-promoted forward simulator; the real datagen path is unchanged"
  - "calibration is a not-yet-gated CalibrationMeta placeholder (empty curve, NaN ECE/MCE, gate=(; gated=false, passed=false)) the per-grid D-05 ship-gate fills"
  - "train_ratio/ood/npe train on the frozen-zt standardized TRAIN split only (leak-free); ratio pairs are assembled over the train pool via the grid-general pair_encode"
metrics:
  tasks_completed: 4
  tasks_total: 4
  files_created: 6
  files_modified: 3
  completed_date: 2026-07-03
---

# Phase 7 Plan 03: Amortized Training Layer (NPE/NRE training + CPU-resident persistence + per-grid pipeline) Summary

**One-liner:** Promoted the NPE architecture + NPE/NRE training into `src/amortized/` with the
D-06 GPU-train / CPU-reproducible split plumbed (default `use_gpu=has_cuda_device()`, throw-guards
removed), migrated shipped-model persistence for BOTH estimators and OOD nulls to CPU-resident
`Flux.state`+`loadmodel!` inside the atomic `.tmp`→integrity→`mv` wrapper (Pitfall 4 / T-7-07),
and factored the whole per-grid sequence into a single idempotent `_train_grid_pipeline(grid; ...)`
that `train_and_register` calls — making the user-definable-grid on-ramp (D-04) real code.

## Status: COMPLETE

All 4 tasks complete and committed atomically (Task 3 under a TDD RED→GREEN gate). `Pkg.test()` is
fully green under `-t auto`: co-resolution 4/4 (NeuralEstimators 0.2.1 / Flux 0.16.10 pins intact),
D-02 10/10, all prior suites unchanged, plus the four new suites — **train_npe 11/11**,
**train_ratio 12/12**, **persist 17/17**, **pipeline 19/19**. `spike/` provably byte-untouched
(`git diff --quiet spike/`).

## What Was Built

### Task 1 — NPE architecture + training, GPU plumbed (`architecture.jl`, `train_npe.jl`) — commit d13e351
- `build_estimator(d_in, D=NPE_D; ...)` — the MLP summary net → `NormalisingFlow` INSTANCE passed
  POSITIONALLY (v0.2.1 gotcha preserved); INPUT-WIDTH-AGNOSTIC (reads `d_in = 2·G²`), so a new
  grid needs no change here. Carries the Phase-5 constants (dstar=64, depth=3/width=256, 10
  coupling layers).
- `train_npe(Ztr_std, Zva_std, θtr, θva; use_gpu=has_cuda_device(), ...)` — fits the leak-free θ
  `ZScoreTransform`, standardizes θ, builds the estimator, runs the FIXED-DATA `train`. **DELTA
  (D-06):** the spike's `use_gpu && throw` guard is REMOVED; LR/weight_decay stay Float64. Returns
  `(estimator, θzt, d_in, variant, arch)` — `arch` is the metadata persist rebuilds from.
- `fit_summary_transform` — the leak-free summary `zt`, fit on TRAIN continuous rows only.

### Task 2 — NRE ratio training, grid-general + GPU plumbed (`train_ratio.jl`) — commit 4baacef
- `build_ratio_estimator(input_dim; num_summaries, width)` — the summary conditioner Chain wrapped
  as `RatioEstimator(net, 1; num_summaries)`; `input_dim` derives from `ratio_input_dim(G)=5·G²`
  (grid-general), NOT the literal 320.
- `assemble_ratio_pairs(Zstd, ρ, n)` — pairs distinct columns of a frozen-`zt` standardized pool
  via the grid-general `pair_encode`; labels are the D-07 ρ_true-contrast split (balanced by
  construction). `measure_log_prior_odds` ported verbatim (measured, never assumed 0).
- `train_ratio(Zstd, ρ; use_gpu=has_cuda_device(), ...)` — **DELTA (D-06):** throw-guard removed;
  **no custom loss** passed (v0.2.1 hard-codes `logitbinarycrossentropy`; a passed loss is silently
  ignored). Returns the ratio HANDLE `(estimator, log_prior_odds, num_summaries, summary_width,
  input_dim)`.

### Task 3 — CPU-resident Flux.state persistence (`persist.jl`) — commits bf58127 (test/RED), 519db96 (feat/GREEN)
- `save_estimator`/`load_estimator` — persist `Flux.state(cpu(est))` + `arch` metadata (NOT the
  whole `PosteriorEstimator`), rebuild via `build_estimator` then `Flux.loadmodel!` on load →
  CPU-resident, device-independent (Pitfall 4 / T-7-07; narrower deserialization surface, T-7-01).
- `save_ratio`/`load_ratio` — the ratio variant (state + input_dim/num_summaries/summary_width +
  measured log_prior_odds).
- `save_ood_nulls`/`load_ood_nulls` — the frozen density (+noise) null fits persist DIRECTLY (plain
  arrays / a `Cholesky`, not a Flux model) through the SAME atomic wrapper.
- All three keep the `.tmp`→reopen-integrity-`@assert`→`mv(...; force=true)` idiom verbatim, add a
  `schema_version` key, and expose `_estimator_ok`/`_ratio_ok`/`_ood_nulls_ok` SKIP-IF-DONE
  predicates (a torn file reports not-done).
- **Round-trip proven:** a tiny trained estimator, saved and reloaded, reproduces its pre-save CPU
  inference **bitwise** under a seeded `posterior_for` (17/17 green).

### Task 4 — Reusable per-grid pipeline + wiring (`pipeline.jl`, `registry.jl`) — commit cb8ac48
- `_train_grid_pipeline(grid; ...)` runs the full sequence ONCE: datagen → leak-free `zt` →
  `train_npe` → `train_ratio` → `fit_ood_nulls` → persist all three CPU-resident under
  `artifacts/grid_G/` → return a populated `EstimatorBundle`. **SKIP-IF-DONE:** valid artifacts are
  loaded instead of retraining. `default_imsize_for` biases finer grids to larger images (16→≥512²,
  32→≥1024²); `default_npairs_for` sizes the pool per grid.
- `registry.train_and_register(grid)` = `_train_grid_pipeline(grid)` then `register!` — the D-04
  on-ramp is now REAL code (the erroring stub is gone; the forward reference resolves at call time).
- The docstring documents the **`-t auto`** datagen requirement and its thread-count-independent
  (Philox-per-index, D-11/D-12) byte-identical reproducibility.

## GPU path (D-06)

The GPU **training** path is plumbed (`use_gpu` defaults to `has_cuda_device()`, throw-guards
removed) and degrades gracefully to CPU. In this execution environment **no CUDA device was
available** (`has_cuda_device() == false` — CUDA is an optional weakdep and was not loaded), so the
GPU path was **NOT exercised at runtime**; every test ran on the CPU-reproducible path. This plan
builds the training machinery only — it does NOT train a shipped grid (07-05+ do). The persisted
nets, ship-gate, and shipped inference are CPU-resident by construction, so training on GPU later
cannot leak device non-determinism into the frozen, pre-registered object (T-7-07).

## Deviations from Plan

### [Rule 2 - Testability seam] Injectable `datagen` keyword on `_train_grid_pipeline`
- **Found during:** Task 4.
- **Issue:** The pipeline's real datagen step calls `generate_samples`, which invokes the Phase-2
  forward simulator (`sample_prior`/`simulate_pair`/`build_mci`) — NOT yet promoted into `src/`
  (an established 07-01/07-02 on-ramp seam). So the plan's acceptance "a tiny fixture run produces
  `artifacts/grid_4/...` and a loadable bundle" could not run end-to-end with the real datagen.
- **Fix:** Added a `datagen` keyword (defaults to the real `generate_samples`) that a caller/test
  injects with a synthetic `(; theta, summary_min)` source. This makes the whole pipeline
  (train→persist→bundle→SKIP-IF-DONE) testable NOW without the simulator, and leaves the real
  datagen path unchanged. The `-t auto` byte-identical-reproducibility property of `generate_samples`
  itself is already validated in 03-02/07-02 (Philox-per-index) and documented in the docstring.
- **Files modified:** src/amortized/pipeline.jl.
- **Commit:** cb8ac48.

### [Design] Ratio handle stored in `EstimatorBundle.ratio`; density-only OOD nulls persisted
- **Found during:** Task 4.
- **Issue:** The 07-01 `EstimatorBundle` struct has no field for the ratio's measured
  `log_prior_odds` (needed by `amortized_log_bf`), and the OOD `:model`/`:noise` channels aren't
  a persistable fit.
- **Fix:** `EstimatorBundle.ratio` (an untyped field) carries the ratio HANDLE NamedTuple so the
  measured prior-odds travels with its estimator; the persisted `ood_nulls` carries the frozen
  density null only. The struct is unchanged (no schema churn). The inference/gate layer re-composes
  the PP `:model` handle from the loaded npe and adds `:noise`/`:zref`/`:thr`.
- **Files modified:** src/amortized/pipeline.jl.
- **Commit:** cb8ac48.

### [Design] TDD RED confirmed statically (persist.jl not yet included), not via a failing full-suite run
- **Found during:** Task 3.
- **Issue:** `Pkg.test()` re-precompiles the heavy Flux/NeuralEstimators/GLMakie stack; a full RED
  run in addition to GREEN is costly (mirrors the 07-02 rationale).
- **Fix:** Committed the RED test FIRST (bf58127) while `persist.jl` was on disk but NOT included in
  the module — the referenced `save_estimator`/… symbols are provably undefined until the GREEN
  include (519db96). RED→GREEN commit order is preserved in history; only the intermediate failing
  full-suite run was skipped for cost.

## Authentication Gates
None.

## Known Stubs
- `_train_grid_pipeline`'s DEFAULT datagen path (`generate_samples`) references the Phase-2
  simulator (`sample_prior`/`simulate_pair`/`build_mci`), promoted by a later Phase-7 plan (the same
  07-01/07-02 seam). The pipeline LOADS and is fully exercised today via the injectable `datagen`
  seam; the real simulator-backed run becomes end-to-end once the simulator is promoted.
- `EstimatorBundle.calibration` is a not-yet-gated `CalibrationMeta` placeholder
  (`gate=(; gated=false, passed=false)`); the per-grid D-05 ship-gate (07-05+) fills the reliability
  curve + ECE/MCE and flips `passed`. Intended hand-off seam, not an accidental gap.

## Threat Flags
None new. The plan's `<threat_model>` mitigations are implemented: T-7-07 (GPU-train reproducibility)
via CPU-resident `Flux.state(cpu(est))` persistence + CPU-default gate/inference; T-7-01 (`.jld2`
deserialization surface) via persisting parameter-array state (not an arbitrary object graph) inside
the integrity-checked atomic write.

## Deferred Items
- Training a SHIPPED grid end-to-end (needs the promoted simulator) — 07-05+ per-grid ship-gate plans.
- Filling `EstimatorBundle.calibration` with a real gate report — 07-05+ (D-05).
- `_lazy_load_from_artifact!` (registry.jl) Artifacts-backed lazy download — a later Phase-7 plan;
  `_train_grid_pipeline` already provides the in-process build + persist path it will front.

## Self-Check

Created files:
- FOUND: src/amortized/architecture.jl
- FOUND: src/amortized/train_npe.jl
- FOUND: src/amortized/train_ratio.jl
- FOUND: src/amortized/persist.jl
- FOUND: src/amortized/pipeline.jl
- FOUND: .planning/phases/07-productionization-conditional-on-go/07-03-SUMMARY.md

Commits:
- FOUND: d13e351 (Task 1 feat — NPE architecture + training)
- FOUND: 4baacef (Task 2 feat — NRE ratio training)
- FOUND: bf58127 (Task 3 test/RED — persistence round-trip)
- FOUND: 519db96 (Task 3 feat/GREEN — CPU-resident persistence)
- FOUND: cb8ac48 (Task 4 feat — _train_grid_pipeline + train_and_register wiring)

Verifications:
- `Pkg.test()` fully green under `-t auto` (train_npe 11/11, train_ratio 12/12, persist 17/17,
  pipeline 19/19; all prior suites unchanged; co-resolution gate 4/4 — NeuralEstimators 0.2.1 /
  Flux 0.16.10 pins intact, no new external deps).
- Task 4 verify command: `_train_grid_pipeline` isdefined + `train_and_register` has an Int method
  (nthreads=32, has_cuda_device=false).
- spike/Project.toml + spike/Manifest.toml + all spike/ sources: provably UNTOUCHED
  (`git diff --quiet spike/`).

## Self-Check: PASSED
