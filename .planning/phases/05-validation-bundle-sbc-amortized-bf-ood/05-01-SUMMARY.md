---
phase: 05-validation-bundle-sbc-amortized-bf-ood
plan: 01
subsystem: testing
tags: [sbc, calibration, neural-estimators, random123, hypothesistests, cairomakie, jld2, amortized-inference]

# Dependency graph
requires:
  - phase: 04-npe-training-advi-benchmark-ablation
    provides: frozen trained_npe.jld2 (PosteriorEstimator + zt/θzt transforms), infer.jl amortized read surface (posterior_for/rho_draws/delta_rho/standardize_summary), load_npe
  - phase: 03-training-data-pipeline
    provides: Loader row-partition, encode_d01, Random123 keyed-seeding idiom (seeding.jl)
  - phase: 02-simulator
    provides: sample_prior, simulate_pair, contract.jl build_mci/patch_summary (read-only src/), CairoMakie headless-save idiom
provides:
  - "spike/validation/consts.jl — single pre-registered SBC/BF/OOD constants source + VAL_MASTER_SEED (anti-snooping contract)"
  - "spike/validation/harness.jl — shared θ*~π→simulate→infer chain + load_frozen_model + disjoint val_rng + paired Δρ path"
  - "spike/validation/sbc.jl — 7θ+Δρ rank table, KS/χ² uniformity, coverage curve, ECE/MCE traffic-light, SBC_CAPTION"
  - "spike/validation/figures.jl — shared CairoMakie surface (rank hist/coverage/reliability + ROC/BF stubs for Plans 02/03)"
  - "test_sbc.jl SBC fast fixture gate + test_bf.jl/test_ood.jl scaffolds wired into single runtests.jl gate with a Phase-5 resolve-risk clause"
affects: [05-02-amortized-bf, 05-03-ood, 05-04-reported-runs, phase-06-go-no-go-memo]

# Tech tracking
tech-stack:
  added: []  # no new dependencies — all runtime deps already present from Phases 1-4 (shared-env isolation honored)
  patterns:
    - "Pre-registered constants in one committed consts.jl BEFORE any reported run (anti-snooping)"
    - "Fresh disjoint Random123 stream (VAL_MASTER_SEED ⊻ VAL_SALT) reserved for reported runs; separate VAL_FIX_SEED for fast gates"
    - "Shared harness composes the frozen amortized read surface — never re-fits standardization"
    - "Fast fixture gate (small M/L, 64×64) vs reported run at locked M/L (Wave 3)"

key-files:
  created:
    - spike/validation/consts.jl
    - spike/validation/harness.jl
    - spike/validation/sbc.jl
    - spike/validation/figures.jl
    - spike/test/test_sbc.jl
    - spike/test/test_bf.jl
    - spike/test/test_ood.jl
  modified:
    - spike/test/runtests.jl

key-decisions:
  - "VAL_MASTER_SEED=0x5BC0FFEE reserved for the reported Wave-3 run; disjoint from NPE_MASTER_SEED=0xC0FFEE (D-02); VAL_FIX_SEED=0xF1F7ED drives fast gates so they never consume the reported stream"
  - "SBC_L=999 with SBC_BINS=50 (L+1=1000 divisible → even rank bins, Pitfall 4)"
  - "SBC calibration routed through the ported _bin_calibration via central-credible-interval coverage; ECE scored by pre-registered traffic-light cutoffs (0.05/0.10)"
  - "figures.jl owned by 05-01 with ROC/BF stub signatures so Wave-2 plans call, never edit it"
  - "test_sbc.jl scaffold created in Task 1 (not just Task 3) so the single runtests.jl gate stays green after every task commit"

patterns-established:
  - "Anti-snooping pre-registration: all M/L/bins/thresholds + seeds locked in consts.jl, committed before any reported run"
  - "PIT rank-uniformity via HypothesisTests ExactOneSampleKSTest + one-way ChisqTest(counts) — never hand-rolled"
  - "Bidirectional test: uniformity asserted finite on the fixture AND p≈0 on a deliberately biased (all-zeros) rank vector"

requirements-completed: [SBC-01, SBC-02, SBC-03, SBC-04]

# Metrics
duration: 22min
completed: 2026-07-02
---

# Phase 5 Plan 01: SBC Calibration Harness + Diagnostics Summary

**Shared θ*~π→simulate→infer harness over the frozen Phase-4 net delivering the SBC calibration proof: 7θ + paired-Δρ rank histograms, KS/χ² uniformity, coverage curve, and ECE/MCE traffic-light — with all thresholds and a disjoint reported seed pre-registered in a committed consts.jl before any reported run.**

## Performance

- **Duration:** ~22 min
- **Started:** 2026-07-02T10:00:00Z
- **Completed:** 2026-07-02T10:22:00Z
- **Tasks:** 3
- **Files modified:** 8 (7 created, 1 modified)

## Accomplishments
- Pre-registered `consts.jl` as the single anti-snooping source for the whole phase (SBC/BF/OOD constants + `VAL_MASTER_SEED` disjoint from `NPE_MASTER_SEED`), committed BEFORE any reported run.
- Shared `harness.jl` that loads the frozen net once, threads a provably-disjoint Random123 stream, and exposes the standard and paired (Δρ) draw→simulate→infer chains reused by SBC/BF/OOD — never re-fitting standardization.
- Full SBC diagnostics in `sbc.jl`: M×8 rank table (7 θ + dedicated paired Δρ, D-01), KS + χ² rank uniformity via HypothesisTests, nominal-vs-empirical coverage curve, and a ported `_bin_calibration`/`CalibrationResult` ECE/MCE with a green/yellow/red traffic-light, plus the `SBC_CAPTION` "calibrated under the simulator; pair with the OOD result" (SBC-04).
- Shared `figures.jl` CairoMakie surface (rank histogram, coverage, reliability/ECE) plus ROC and BF-agreement stub signatures owned here so Wave-2 plans call them without editing the module.
- Extended the single `runtests.jl` gate with three Phase-5 test files and a Phase-5 resolve-risk clause (NeuralEstimators pinned 0.2.1; no ROC/flow/Turing package may enter the env). Suite green.

## Task Commits

Each task was committed atomically:

1. **Task 1: Pre-registration scaffold (consts.jl + test scaffolds + runtests wiring + resolve-risk gate)** - `e9c91d3` (feat)
2. **Task 2: Shared harness.jl (draw→simulate→infer + frozen model + disjoint Random123 + paired Δρ)** - `feee8d5` (feat)
3. **Task 3: sbc.jl + figures.jl + SBC gate (7θ+Δρ ranks, KS/χ² uniformity, coverage, ECE traffic-light)** - `0910b9a` (feat, TDD)

**Plan metadata:** _(final docs commit — this SUMMARY + STATE + ROADMAP)_

## Files Created/Modified
- `spike/validation/consts.jl` - Single pre-registered SBC/BF/OOD constants + VAL_MASTER_SEED/VAL_FIX_SEED (anti-snooping contract; Plans 02/03 read from here)
- `spike/validation/harness.jl` - load_frozen_model, val_rng (disjoint stream), draw_simulate_infer, draw_simulate_infer_paired
- `spike/validation/sbc.jl` - sbc_ranks, sbc_uniformity, sbc_coverage, _bin_calibration/CalibrationResult, sbc_traffic_light, sbc_calibration, SBC_CAPTION
- `spike/validation/figures.jl` - plot_rank_histogram/plot_coverage_curve/plot_reliability_ece + plot_roc/plot_bf_agreement stubs (headless PNG save)
- `spike/test/test_sbc.jl` - SBC fast fixture gate (SC1..SC4 + biased-input p≈0 + anti-snooping seed assertion)
- `spike/test/test_bf.jl` - Wave-2 BF gate scaffold (skipped placeholder)
- `spike/test/test_ood.jl` - Wave-2 OOD gate scaffold (skipped placeholder)
- `spike/test/runtests.jl` - Wired the three Phase-5 test files + Phase-5 resolve-risk clause

## Decisions Made
- **Reserved reported seed vs fixture seed.** `VAL_MASTER_SEED=0x5BC0FFEE` (disjoint from `NPE_MASTER_SEED=0xC0FFEE` and from the training/holdout/fold salts via a new `VAL_SALT`) is consumed only by the reported Wave-3 run; the fast per-task gates draw from `VAL_FIX_SEED=0xF1F7ED` so they can never pre-observe or tune against the reported stream.
- **SBC calibration ECE via coverage routing.** The ported `_bin_calibration` is fed the central-credible-interval coverage (predicted = nominal level, positive = in-interval indicator), so the traffic-light ECE measures coverage miscalibration — the natural SBC calibration metric — while still being the verbatim-ported reusable binning function.
- **Coverage curve is rank/PIT-based** (`|u−0.5| ≤ α/2`), guaranteeing monotone non-decreasing empirical coverage by construction (asserted in the gate).
- **test_sbc.jl scaffold created in Task 1.** Since Task 1 wires `test_sbc.jl` into `runtests.jl`, a skipped scaffold was created in the same commit so the single gate stays green after every task; Task 3 then replaced it with the real gate.

## Deviations from Plan

**None affecting behavior.** One structural clarification of the plan's task boundaries:

**1. [Rule 3 - Blocking] Created test_sbc.jl scaffold in Task 1 (plan lists it only under Task 3 files)**
- **Found during:** Task 1 (runtests.jl wiring)
- **Issue:** Task 1's acceptance criterion requires `runtests.jl` to include `test_sbc.jl` AND exit 0 with the Phase-5 testsets enumerated, but `test_sbc.jl` is only listed among Task 3's files. Including a non-existent file would break the gate at the Task 1 commit.
- **Fix:** Created `test_sbc.jl` as a minimal skipped scaffold (matching test_bf.jl/test_ood.jl) in the Task 1 commit; Task 3 replaced its contents with the real fixture gate.
- **Files modified:** spike/test/test_sbc.jl
- **Verification:** `julia --project=spike spike/test/runtests.jl` exits 0 after Task 1 (SBC/BF/OOD testsets skipped) and after Task 3 (SBC real, 28/28).
- **Committed in:** e9c91d3 (Task 1), then filled in 0910b9a (Task 3)

---

**Total deviations:** 1 (structural task-boundary clarification; no scope creep, no threshold changes).
**Impact on plan:** None on scientific content. Every pre-registered value is exactly as specified in the plan; no threshold was tuned.

## TDD Gate Compliance

Task 3 is `tdd="true"`. Because this executor runs sequentially on `main` with commit hooks and must keep the full `runtests.jl` gate green at every commit, the RED (failing test) and GREEN (implementation) states were NOT split into two separate commits — a test-only commit referencing not-yet-existing `sbc.jl` would leave `main` with a failing gate. Instead, `sbc.jl` + `figures.jl` + the real `test_sbc.jl` were committed together in the GREEN state (`0910b9a`) after verifying the fixture gate passes 28/28. The behavior contract from the plan's `<behavior>` block is fully exercised, including the bidirectional uniformity test (finite p on the fixture, p≈0 on a biased all-zeros rank vector).

## Issues Encountered
- Verified the installed `HypothesisTests` `ChisqTest(counts::Vector{Int})` one-way goodness-of-fit form (plan note A3) before wiring — confirmed it returns high p on uniform counts and p≈0 on degenerate counts.
- Ran under concurrent Phase-9/10 agents mutating the shared `spike/` env; honored isolation strictly (no `Pkg.add/update/rm`, stayed within `spike/validation/` and `spike/test/`, used the on-disk frozen net and caches in place). The Phase-5 resolve-risk clause re-asserts NeuralEstimators 0.2.1 after the validation code loads.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Plans 05-02 (amortized BF) and 05-03 (OOD) can now read pre-registered constants from `consts.jl`, ride `harness.jl` (standard + paired chains), and call the `figures.jl` ROC/BF-agreement stubs without editing shared files.
- The reported SBC run (Wave 3, 05-04) has the locked M=2000/L=999 and the reserved `VAL_MASTER_SEED` stream ready; it must report honestly whether uniformity/coverage/ECE pass — thresholds must not be adjusted post-hoc.
- Structural bound to carry forward (STATE blocker): OOD detection power is limited by the fixed patch-correlation summary; summary-orthogonal misspecifications must be named, not hidden (relevant to 05-03 and the SBC-04 caption pairing).

## Self-Check: PASSED

- All 7 created files present on disk (consts.jl, harness.jl, sbc.jl, figures.jl, test_sbc.jl, test_bf.jl, test_ood.jl).
- All 3 task commits verified in git history (e9c91d3, feee8d5, 0910b9a).
- Full `runtests.jl` gate exits 0 with the SBC testset real (28/28) and BF/OOD scaffolds skipped.

---
*Phase: 05-validation-bundle-sbc-amortized-bf-ood*
*Completed: 2026-07-02*
