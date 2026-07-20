---
phase: 07-productionization-conditional-on-go
plan: 06
subsystem: testing
tags: [neural-estimators, flux, sbc, bayes-factor, ood, ship-gate, npe, nre, amortized-inference, reproducibility, user-definable-grid]

# Dependency graph
requires:
  - phase: 07-03
    provides: train_and_register + _train_grid_pipeline (the public user-definable-grid path)
  - phase: 07-04
    provides: per-grid CPU ship-gate machinery (run_gate.jl, harness.jl, sbc.jl, gate_consts_template.jl)
  - phase: 07-05
    provides: proven 8×8 pipeline shape + reference gate-8x8.md verdict to interpret against
provides:
  - "4×4 NPE+NRE+OOD-nulls bundle produced by the PUBLIC train_and_register(4) path, CPU-resident (artifacts/grid_4/*.jld2, on-disk)"
  - "Recorded 4×4 CPU SBC/BF/OOD ship-gate outcome (gate-4x4.md) on the fresh disjoint PROD_SEED[4]=0x8c0ad97b99bd6031"
  - "End-to-end proof of the PROD-02/D-04 user-definable-grid happy path (train_and_register + register! + in-process estimator_for)"
  - "Hardened non-clamped KDE BF baseline (bf.jl) against QuadGK domain-overshoot crash"
affects: [07-10, registry-population, per-grid-plans]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Per-grid training driven through the PUBLIC train_and_register(grid) helper (not inline datagen→train→persist)"
    - "Ship-gate reads the grid's OWN locked gate_consts_G.jl + fresh disjoint PROD_SEED[G] (anti-snooping)"
    - "Compute-budget adaptation: constrain train-time imsize to the gate's SBC_IMSIZE, keep pre-registered constants untouched"

key-files:
  created:
    - test/gate/gate_consts_4.jl
    - .planning/phases/07-productionization-conditional-on-go/gate-4x4.md
    - artifacts/grid_4/npe_4.jld2 (on-disk, gitignored)
    - artifacts/grid_4/ratio_4.jld2 (on-disk, gitignored)
    - artifacts/grid_4/ood_nulls_4.jld2 (on-disk, gitignored)
    - artifacts/grid_4/gate_report_4.jld2 (on-disk, gitignored)
  modified:
    - src/amortized/bf.jl
    - test/runtests.jl

key-decisions:
  - "Drove the 4×4 bundle through the PUBLIC train_and_register(4) path (not _train_grid_pipeline directly) — this plan is the designated end-to-end exercise of the PROD-02/D-04 user-definable-grid happy path; estimator_for(4) confirmed the in-process registration"
  - "Constrained datagen imsize_set to ((256,256),) — matches the gate's SBC_IMSIZE=(256,256); 4×4 is robust there (64×64 px/patch ≥4096 px); the default mixed set (up to 2048²) makes CPU datagen intractable"
  - "Kept plan-default n_pairs=30000 (default_npairs_for(4)); did NOT under-train"
  - "Did NOT force-commit artifacts/grid_4/*.jld2 — repo policy (.gitignore:63 artifacts/, 07-05 precedent) treats trained nets + reported .jld2 as regenerable caches, NOT deliverables"
  - "Recorded honest FAIL verdict on the byte-locked gate_consts_4.jl; NOT re-tuned to force a pass"
  - "Hardened kde_log_bf_unclamped with a [0,1] domain clamp (numerical-safety Rule-1 fix) — NOT the forbidden 1e-8 floor; a saturated tail yields the honest ±Inf (dropped), not a finite ±18.4"

patterns-established:
  - "The public train_and_register(grid) helper is the shipping on-ramp exercised end-to-end per grid; the gate rides the grid's own fresh disjoint seed"

requirements-completed: [PROD-02]

# Metrics
duration: 40min
completed: 2026-07-20
---

# Phase 7 Plan 06: 4×4 Per-Grid Ship-Gate Summary

**4×4 NPE+NRE+OOD-nulls bundle produced through the PUBLIC `train_and_register(4)` path (30k pairs,
256², datagen under `-t auto`) — proving the PROD-02/D-04 user-definable-grid happy path end-to-end
— then run through the CPU SBC/BF/OOD ship-gate on the fresh disjoint `PROD_SEED[4]`: ECE green on
all 8 SBC params with 7/8 KS pass (stronger than 8×8's 5/8), honest FAIL on the strict all-8 KS
conjunction and the BF corr/tol thresholds; a QuadGK domain-overshoot crash in the non-clamped KDE
baseline was fixed inline.**

## Performance

- **Duration:** ~40 min (train ~2.7 min compute + gate ~6 min compute + one gate crash-and-fix cycle)
- **Completed:** 2026-07-20
- **Tasks:** 3 executed (Task 1 train via public path, Task 2 gate_consts_4.jl, Task 3 gate)
- **Files committed:** gate_consts_4.jl, gate-4x4.md, bf.jl, runtests.jl (4 on-disk artifacts gitignored)

## Accomplishments
- Produced + registered the 4×4 estimator bundle **through the public `train_and_register(4)`**
  helper (`use_gpu=has_cuda_device()`==false → CPU), datagen launched under `julia -t auto`; NPE
  early-stopped @ epoch 61, NRE trained (num_summaries 64, `ratio_input_dim(4)`=80). `estimator_for(4)`
  returned the same registered bundle in-process — the user-definable-grid path works end-to-end.
- All three artifacts (`npe_4`/`ratio_4`/`ood_nulls_4`) load CPU-only and yield finite draws; the
  4×4 store dir is distinct from grid_8 (`summary_dim(4)=32 ≠ summary_dim(8)=128`).
- Committed the fresh 4×4 pre-registration `gate_consts_4.jl` (disjoint `PROD_SEED[4]=0x8c0ad97b99bd6031`)
  BEFORE the reported gate run.
- Ran the determinism-pinned (single-thread, `use_gpu=false`) SBC/BF/OOD ship-gate keyed by
  `PROD_SEED[4]`, wrote the atomic `gate_report_4.jld2` + human-readable `gate-4x4.md`, and recorded
  the honest un-tuned verdict.

## Task Commits

Each task committed atomically:

1. **Task 2: Fresh 4×4 pre-registration (gate_consts_4.jl)** — `278ed60` (feat) — committed before the gate run; `PROD_SEED[4]` disjoint from the two forbidden spike seeds, `DEFAULT_MASTER_SEED`, and `PROD_SEED[8]`.
2. **Task 1: Produce the 4×4 bundle via `train_and_register(4)`** — no commit (outputs are `artifacts/grid_4/*.jld2`, gitignored per repo policy; on-disk + load-verified).
3. **Rule-1 fix (enabling Task 3): harden `kde_log_bf_unclamped`** — `4ed533b` (fix) — bf.jl domain clamp + runtests regression test.
4. **Task 3: Run + record the 4×4 CPU SBC/BF/OOD ship-gate** — `fc87400` (feat) — gate-4x4.md (the `gate_report_4.jld2` is gitignored).

## Decisions Made
- **Public path (not the internal pipeline):** the plan designates 07-06 as the end-to-end exercise
  of the PROD-02/D-04 user-definable-grid path, so the bundle was produced via
  `ProteinCoLoc.train_and_register(4; use_gpu=..., imsize_set=((256,256),))` and the in-process
  `estimator_for(4)` return was verified — not by calling `_train_grid_pipeline` directly.
- **imsize constrained to 256²:** `default_imsize_for(4)` still spans up to 2048² (the 256²-heavy
  set), which makes CPU datagen intractable in the budget. Constraining to `((256,256),)` both fits
  the budget and aligns training image size with the gate's `SBC_IMSIZE=(256,256)`; 4×4 is robust at
  256² (each 64×64 patch carries ≥4096 px, far above the ≥15-survivor floor). Gate constants untouched.
- **Kept n_pairs=30000 (plan default):** `default_npairs_for(4)`; full pipeline completed in ~2.7 min
  compute, so the default pool was used rather than under-training.
- **Artifacts not force-committed:** `.gitignore:63` ignores `artifacts/`; the recorded markdown
  numbers are the deliverable (07-05 precedent).
- **Honest FAIL recorded:** `gate_consts_4.jl` byte-locked, NOT re-tuned. Residual failures read to
  non-method causes per the Phase-6 memo §6 framing.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Hardened the non-clamped KDE BF baseline against a QuadGK domain-overshoot crash**
- **Found during:** Task 3 (first gate run)
- **Issue:** `run_gate.jl --grid 4 --bf` crashed with a `DomainError` in `kde_log_bf_unclamped`
  (bf.jl:138): `log` was called with a tiny negative argument (`-4.79e-9`). `_p_gt_threshold_unclamped`
  returns `1 - p_le`, where `p_le` is a QuadGK adaptive integral of the KDE that can overshoot a few
  ulp past 1.0 — making the returned probability a tiny negative, which then makes `posterior_odds`
  negative and crashes `log(·)`. (The 8×8 sweep happened to land on exactly-zero → `log(0)=-Inf`,
  already dropped by the gate's `isfinite` filter, so 8×8 never hit this.)
- **Fix:** Clamp the returned tail probability to its mathematically valid `[0,1]` domain
  (`clamp(1 - p_le, 0.0, 1.0)`). This maps out-of-domain roundoff to an exact boundary (0 or 1) →
  `±Inf` logBF, which the gate drops as non-finite — the honest, un-floored one-sided value the
  docstring already promises. It is **NOT** the forbidden `1e-8` clamp (which would introduce a finite
  ±18.4 floor, distorting `max|Δ logBF|`); in-range finite values are byte-identical. Added a
  regression test (fully one-sided negative posterior is domain-safe).
- **Files modified:** `src/amortized/bf.jl`, `test/runtests.jl`
- **Verification:** standalone repro — one-sided positive posterior still exceeds the clamp ceiling,
  one-sided negative posterior no longer throws and yields a large-negative/-Inf value; the full gate
  then ran to completion and wrote `gate_report_4.jld2`.
- **Committed in:** `4ed533b`

**2. [Rule 3 - Blocking] Constrained datagen image size to fit the CPU compute budget**
- **Found during:** Task 1 (produce the 4×4 bundle)
- **Issue:** `train_and_register(4)`'s default `default_imsize_for(4)` image set (256²-heavy, up to
  2048²) would not complete 30k-pair CPU datagen within the compute window (the same intractability
  documented at 8×8).
- **Fix:** Passed `imsize_set=((256,256),)` (256²-only), matching the gate's pre-registered
  `SBC_IMSIZE=(256,256)`; kept plan-default `n_pairs=30000`. The gate constants
  (M/L/thresholds/PROD_SEED) were untouched.
- **Files modified:** none tracked (produces on-disk `artifacts/grid_4/*.jld2`)
- **Verification:** all three artifacts load via `load_estimator`/`load_ratio`/`load_ood_nulls` on a
  CPU-only process and produce finite posterior draws; the gate ran to completion against them.
- **Committed in:** n/a (artifacts gitignored)

**3. [Rule 3 - Repo convention] Did not commit the trained/report .jld2 artifacts**
- **Found during:** Task 1 / Task 3 commit steps
- **Issue:** The plan's `files_modified` lists `artifacts/grid_4/*.jld2`, but `.gitignore:63` ignores
  `artifacts/` and the 07-05 precedent treats trained nets + reported `.jld2` as regenerable caches.
- **Fix:** Left the artifacts on-disk; committed only the human-readable `gate-4x4.md` + `gate_consts_4.jl`.
- **Files modified:** n/a
- **Verification:** `git check-ignore artifacts/grid_4/gate_report_4.jld2` confirms ignored.
- **Committed in:** n/a

---

**Total deviations:** 3 auto-fixed (1 Rule-1 numerical bug in the shipped BF baseline; 2 Rule-3 —
compute-budget blocker + repo convention). No pre-registration constants were altered.
**Impact on plan:** No scope creep. The scientific pre-registration (`gate_consts_4.jl`) and the gate
verdict are untouched by the deviations; the bf.jl fix is a genuine numerical-safety bug fix that
leaves all in-range values byte-identical.

## Gate Verdict (honest, un-tuned)

| Gate | Verdict | Headline numbers |
|------|---------|------------------|
| SBC (M=2000, L=999) | **FAIL** (ks_pass=false, ece_pass=true) | ECE green on all 8 (ρ_true 0.0120, Δρ 0.0074); **KS 7/8 pass** (only shift_dx 0.0326) — stronger than 8×8's 5/8 |
| BF (n=20 finite/25) | **FAIL** | corr 0.9410 (near-miss of 0.95), max\|Δ logBF\| 12.40 (KDE-tail artifact) |
| OOD (n=200) | indeterminate | ID operating point id_threshold=31.79; per-family AUC deferred to later misspec-ROC plan |

Every residual failure reproduces a documented non-method cause (memo §6): ECE-green calibration with
M=2000 KS over-sensitivity on the summary-uninformative `shift_dx` nuisance parameter; BF mid-range
agreement with a finite-sample KDE-baseline tail divergence. Constants byte-locked.

## Next Phase Readiness
- The PROD-02/D-04 **user-definable-grid happy path is proven end-to-end**: `train_and_register(4)`
  ran the full pipeline + `register!`, and `estimator_for(4)` returned the bundle in-process.
- The 4×4 estimator bundle exists on-disk (`artifacts/grid_4/`, gitignored). 07-10 (registry
  population) owns artifact persistence/distribution and can regenerate via `train_and_register(4)`
  (SKIP-IF-DONE loads if present, else retrains).
- Honest FAIL verdict surfaced verbatim — consistent with the Phase-6 "Clean Go" framing; not a
  reason to loop. The 4×4 SBC axis is notably stronger than the 8×8 reference (7/8 vs 5/8 KS pass).

## Known Stubs
None — the trained bundle is real (finite CPU draws), the gate report is real, and the recorded
numbers are the honest gate output. The OOD per-family AUC being `nothing` is by-design deferral (the
bare CLI injects no misspecification families), not a stub.

## Self-Check: PASSED

- On-disk artifacts: `npe_4.jld2`, `ratio_4.jld2`, `ood_nulls_4.jld2`, `gate_report_4.jld2` — all FOUND (gitignored).
- Committed deliverables: `test/gate/gate_consts_4.jl` TRACKED; `.planning/.../gate-4x4.md` TRACKED; `src/amortized/bf.jl` + `test/runtests.jl` modified & committed.
- Commits: `278ed60` (gate_consts_4.jl), `4ed533b` (bf.jl fix), `fc87400` (gate-4x4.md) — all FOUND.

---
*Phase: 07-productionization-conditional-on-go*
*Completed: 2026-07-20*
