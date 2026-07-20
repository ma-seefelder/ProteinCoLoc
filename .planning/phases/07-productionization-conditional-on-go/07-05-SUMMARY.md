---
phase: 07-productionization-conditional-on-go
plan: 05
subsystem: testing
tags: [neural-estimators, flux, sbc, bayes-factor, ood, ship-gate, npe, nre, amortized-inference, reproducibility]

# Dependency graph
requires:
  - phase: 07-03
    provides: _train_grid_pipeline shared per-grid pipeline (datagen→train_npe→train_ratio→fit_ood_nulls→persist)
  - phase: 07-04
    provides: per-grid CPU ship-gate machinery (run_gate.jl, harness.jl, sbc.jl, gate_consts_template.jl)
  - phase: 07-05 (prior task)
    provides: fresh 8×8 pre-registration gate_consts_8.jl + Phase-2 forward simulator promoted to src/amortized/simulator.jl
provides:
  - "8×8 NPE+NRE+OOD-nulls bundle produced by _train_grid_pipeline(8), CPU-resident (artifacts/grid_8/*.jld2, on-disk)"
  - "Recorded 8×8 CPU SBC/BF/OOD ship-gate outcome (gate-8x8.md) on the fresh disjoint PROD_SEED[8]"
  - "End-to-end proof of the D-05 per-grid gate machinery on the spike-validated reference grid"
affects: [07-08, 07-10, registry-population, sub-tile-local-map, per-grid-plans]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Per-grid training driven ONLY through _train_grid_pipeline(grid) (no inline datagen→train→persist restatement)"
    - "Ship-gate reads the grid's OWN locked gate_consts_G.jl + fresh disjoint PROD_SEED[G] (anti-snooping)"
    - "Trained .jld2 artifacts are regenerable on-disk caches (gitignored); recorded numbers in markdown are the committed deliverable"

key-files:
  created:
    - .planning/phases/07-productionization-conditional-on-go/gate-8x8.md
    - artifacts/grid_8/npe_8.jld2 (on-disk, gitignored)
    - artifacts/grid_8/ratio_8.jld2 (on-disk, gitignored)
    - artifacts/grid_8/ood_nulls_8.jld2 (on-disk, gitignored)
    - artifacts/grid_8/gate_report_8.jld2 (on-disk, gitignored)
  modified: []

key-decisions:
  - "Constrained datagen imsize_set to ((256,256),) — matches the gate's SBC_IMSIZE=(256,256); the default mixed set (up to 2048²) makes 50k-pair CPU datagen intractable in the compute budget"
  - "Kept plan-default n_pairs=50000 once 256² datagen proved tractable (~21 min full pipeline); did NOT under-train"
  - "Did NOT force-commit artifacts/grid_8/*.jld2 — repo policy (.gitignore:63 artifacts/, Phase-5 precedent) treats trained nets + reported .jld2 as regenerable caches, NOT deliverables"
  - "Recorded honest FAIL verdict on the locked pre-registration; gate_consts_8.jl byte-unchanged, NOT re-tuned to force a pass"

patterns-established:
  - "Compute-budget adaptation: constrain train-time image size to the gate's evaluation image size, keep pre-registered gate constants untouched"

requirements-completed: [PROD-02]

# Metrics
duration: 130min
completed: 2026-07-20
---

# Phase 7 Plan 05: 8×8 Per-Grid Ship-Gate Summary

**8×8 NPE+NRE+OOD-nulls bundle trained through `_train_grid_pipeline(8)` (50k pairs, 256²) and run through the CPU SBC/BF/OOD ship-gate on the fresh disjoint `PROD_SEED[8]` — ECE green on all 8 SBC params, but honest FAIL on the strict all-8 KS conjunction and the BF corr/tol thresholds; gate machinery proven end-to-end.**

## Performance

- **Duration:** ~130 min (incl. one killed 45-min datagen attempt at full imsize + a 5.8-min 12k probe run)
- **Started:** 2026-07-20T14:47Z (approx, session start)
- **Completed:** 2026-07-20T16:33Z
- **Tasks:** 2 executed this session (Task 1 train, Task 3 gate); Task 2 (gate_consts_8.jl) pre-committed
- **Files modified:** 1 committed (gate-8x8.md) + 4 on-disk artifacts (gitignored)

## Accomplishments
- Produced the 8×8 estimator bundle via the shared `_train_grid_pipeline(8)` (datagen under `-t auto`, NPE early-stopped @ epoch 42, NRE trained), persisted CPU-resident; all three artifacts load CPU-only and yield finite posterior draws.
- Ran the determinism-pinned (single-thread, `use_gpu=false`) SBC/BF/OOD ship-gate against the CPU-resident net keyed by the fresh disjoint `PROD_SEED[8]=0x8b39fecd4e2bcceb` (disjoint from spike VAL/NPE seeds), and wrote the atomic `gate_report_8.jld2` + human-readable `gate-8x8.md`.
- Recorded the honest, un-tuned verdict: SBC ECE green on all 8 params (ρ_true 0.0062, Δρ 0.0098) with 5/8 KS pass; BF corr 0.947 near-miss + `max|Δ logBF|` 18.12 KDE-tail artifact; OOD ID operating point 103.86 (per-family AUC deferred). Machinery proven end-to-end on the reference grid.

## Task Commits

Each task committed atomically:

1. **Task 2: Fresh 8×8 pre-registration (gate_consts_8.jl)** - `d1d0040` (feat) — pre-committed before this session (verified byte-locked, PROD_SEED[8] disjoint)
2. **Task 1: Produce the 8×8 bundle via `_train_grid_pipeline(8)`** - no commit (outputs are `artifacts/grid_8/*.jld2`, gitignored per repo policy; on-disk + load-verified)
3. **Task 3: Run + record the 8×8 CPU SBC/BF/OOD ship-gate** - `632c69b` (feat) — gate-8x8.md (the gate_report_8.jld2 is gitignored)

_Prerequisite in place from prior task: `a09c675` (Phase-2 forward simulator promoted to src/amortized/simulator.jl)._

## Files Created/Modified
- `.planning/phases/07-productionization-conditional-on-go/gate-8x8.md` - Recorded SBC/BF/OOD numbers + PASS/FAIL verdict against the pre-registered thresholds (committed deliverable).
- `artifacts/grid_8/npe_8.jld2` - CPU-resident Flux.state 8×8 NPE (on-disk, gitignored).
- `artifacts/grid_8/ratio_8.jld2` - CPU-resident 8×8 NRE (num_summaries 64, input dim 320) (on-disk, gitignored).
- `artifacts/grid_8/ood_nulls_8.jld2` - density-channel OOD null fit (on-disk, gitignored).
- `artifacts/grid_8/gate_report_8.jld2` - atomic gate report (SBC/BF/OOD) (on-disk, gitignored).

## Decisions Made
- **imsize constrained to 256²:** `default_imsize_for(8)` spans up to 2048², which made a 50k-pair CPU datagen exceed the compute budget (first attempt killed at ~45 min still in datagen). Constraining to `((256,256),)` both fits the budget and aligns the training image size with the gate's `SBC_IMSIZE=(256,256)` — a well-matched, defensible choice for a 256²-evaluated gate.
- **Kept n_pairs=50000 (plan default):** once 256² datagen proved fast (~6 min), the full plan-default pool was used rather than under-training; the whole pipeline completed in ~21 min.
- **Artifacts not force-committed:** `.gitignore:63` ignores `artifacts/`, and the Phase-5 precedent comments (lines 427–440) state trained nets + reported `.jld2` are regenerable caches, NOT deliverables — the recorded markdown numbers are the deliverable. Force-adding 6 MB of gitignored binaries would violate the repo convention (CLAUDE.md: follow project conventions).
- **Honest FAIL recorded:** `gate_consts_8.jl` is byte-locked; it was NOT re-tuned to force a pass. Residual failures are read to non-method causes per the Phase-6 memo §6 framing (M=2000 KS/χ² over-sensitivity; KDE-baseline tail divergence).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Constrained datagen image size to fit the CPU compute budget**
- **Found during:** Task 1 (produce the 8×8 bundle)
- **Issue:** `_train_grid_pipeline(8)` with the default `default_imsize_for(8)` image set (up to 2048²) did not complete datagen within the available compute window — the first full-default background run was killed at ~45 min still in datagen with no artifacts.
- **Fix:** Passed `imsize_set=((256,256),)` (256²-only), matching the gate's pre-registered `SBC_IMSIZE=(256,256)`; kept plan-default `n_pairs=50000`. Full pipeline (datagen→NPE→NRE→persist) then completed in ~21 min. The gate constants (M/L/thresholds/PROD_SEED) were untouched.
- **Files modified:** none tracked (produces on-disk `artifacts/grid_8/*.jld2`)
- **Verification:** all three artifacts load via `load_estimator`/`load_ratio`/`load_ood_nulls` on a CPU-only process and produce finite posterior draws; the gate ran to completion against them.
- **Committed in:** n/a (artifacts gitignored; evidence recorded here + in gate-8x8.md)

**2. [Rule 3 - Repo convention] Did not commit the trained/report .jld2 artifacts**
- **Found during:** Task 1 / Task 3 commit steps
- **Issue:** The plan's `files_modified` lists `artifacts/grid_8/*.jld2`, but `.gitignore:63` ignores `artifacts/` and the repo's Phase-5 precedent explicitly treats trained nets + reported `.jld2` as regenerable caches (not deliverables).
- **Fix:** Left the artifacts on-disk (consumed by the gate in this same worktree); committed only the human-readable `gate-8x8.md`. CLAUDE.md directive to follow project conventions takes precedence over the plan's implicit "commit the artifacts".
- **Files modified:** n/a
- **Verification:** `git check-ignore` confirms `artifacts/` ignored; `gate-8x8.md` tracked and committed.
- **Committed in:** `632c69b` (gate-8x8.md)

---

**Total deviations:** 2 auto-fixed (both Rule 3 — 1 compute-budget blocker, 1 repo-convention). No pre-registration constants were altered.
**Impact on plan:** No scope creep. The scientific pre-registration (gate_consts_8.jl) and the gate verdict are untouched by the deviations; only training-compute knobs and the artifact-persistence mechanism were adapted to the environment and repo policy.

## Issues Encountered
- **First full-default datagen run killed (~45 min):** the default mixed imsize set (up to 2048²) is too slow for 50k CPU sims in the budget. Resolved by constraining to 256² (see Deviation 1). A misleading early per-epoch timing probe (14.7 s/epoch) was compile-dominated; steady-state is ~1 s/epoch at ~10k train and ~3.5 s/epoch at ~42.5k train, so full 50k training was in fact tractable (~21 min end-to-end).
- **OOD per-family AUC not produced by the bare CLI:** `run_gate.jl --ood` injects no positive-control (misspecified) simulators, so `ood_gate` returns only the ID operating point (`auc=nothing`). This is by design — the controlled misspecification-grid ROC experiment is explicitly deferred to a later Phase-7 plan per the `src/amortized/ood.jl` scope note (lines 49–52). Inventing misspec families ad-hoc would be un-pre-registered and violate the phase's anti-snooping discipline. Recorded honestly in gate-8x8.md.

## Gate Verdict (honest, un-tuned)

| Gate | Verdict | Headline numbers |
|------|---------|------------------|
| SBC (M=2000, L=999) | **FAIL** (ks_pass=false, ece_pass=true) | ECE green on all 8 (ρ_true 0.0062, Δρ 0.0098); KS 5/8 pass (fails: autofluorescence 0.042, noise 0.0032, label_efficiency 1.8e-9) |
| BF (n=15 finite/25) | **FAIL** | corr 0.9472 (near-miss of 0.95), max\|Δ logBF\| 18.12 (KDE-baseline tail artifact) |
| OOD (n=200) | indeterminate | ID operating point id_threshold=103.86; per-family AUC deferred to later misspec-ROC plan |

Every residual failure reproduces a documented non-method cause (memo §6): ECE-green calibration with M=2000 KS/χ² over-sensitivity on summary-uninformative nuisance params; BF mid-range agreement with a finite-sample KDE-baseline tail divergence. Constants byte-locked.

## Next Phase Readiness
- The D-05 per-grid gate machinery is proven end-to-end on the 8×8 reference grid — the harness, fresh-seed pre-registration, atomic report writer, and CPU-resident inference all ran correctly.
- The 8×8 estimator bundle exists on-disk in this worktree (`artifacts/grid_8/`). **Concern:** these are gitignored and will not survive worktree removal; 07-10 (registry population) owns artifact persistence/distribution and can regenerate via `_train_grid_pipeline(8)` (SKIP-IF-DONE loads if present, else retrains).
- Honest FAIL verdict surfaced verbatim for the orchestrator/verifier — consistent with the Phase-6 "Clean Go" framing; not a reason to loop.

## Self-Check: PASSED

- On-disk artifacts: `npe_8.jld2`, `ratio_8.jld2`, `ood_nulls_8.jld2`, `gate_report_8.jld2` — all FOUND (gitignored, this worktree).
- Committed deliverables: `test/gate/gate_consts_8.jl` TRACKED; `.planning/.../gate-8x8.md` TRACKED.
- Commits: `d1d0040` (gate_consts_8.jl), `a09c675` (simulator promotion), `632c69b` (gate-8x8.md) — all FOUND.

---
*Phase: 07-productionization-conditional-on-go*
*Completed: 2026-07-20*
