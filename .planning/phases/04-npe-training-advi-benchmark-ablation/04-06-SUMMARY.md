---
phase: 04-npe-training-advi-benchmark-ablation
plan: 06
subsystem: npe
tags: [ablation, summary-statistic, k-fold-cv, rmse, assess, decision-rule, ood, min-vs-aug, cpu-only, zero-resim]

# Dependency graph
requires:
  - phase: 04-npe-training-advi-benchmark-ablation
    plan: 03
    provides: "build_estimator (arch), train_fold (fixed-data CPU-only training + leak-free theta ZScoreTransform), infer surface"
  - phase: 04-npe-training-advi-benchmark-ablation
    plan: 05
    provides: "assess/rmse per-parameter RMSE pattern (theta_zt.scale physical-unit scaling), rmse_report scoring analog"
  - phase: 04-npe-training-advi-benchmark-ablation
    plan: 01
    provides: "test_npe.jl scaffold + pre-registered consts ABL_REL_MARGIN=0.05, ABL_FOLD_CONSISTENCY=4, NPE_MASTER_SEED"
  - phase: 03-training-data-pipeline
    provides: "Loader.load_fold(...; variant=:min|:aug) leak-free k=5 CV over both cached summaries (zero re-simulation)"
provides:
  - "spike/npe/ablation.jl — ablate(dir) :min-vs-:aug k=5 CV per-parameter RMSE (K x 2 x 7 table + rho_rmse_min/aug + fold_wins_aug) and choose_summary(result) implementing the D-07 rule"
  - "spike/npe/summary_choice.md — ABL-02 written justification: chosen :min, rule inputs, all-7 theta RMSE, explicit OOD-detectability coupling to Phase 5"
  - "SC4 (ABL-01) + SC5 (ABL-02) green in spike/test/test_npe.jl"
  - "spike/data/loader.jl zero-variance z-score guard (constant continuous feature -> centered 0 instead of NaN); benefits every :aug consumer"
affects: [phase-5-sbc-bf-ood]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "arch-controlled ablation: train_fold with variant=:min|:aug on the SAME build_estimator/master_seed/hyperparams — only the cached summary differs (D-06, zero re-simulation)"
    - "per-fold per-variant per-parameter RMSE via NeuralEstimators.assess/rmse in standardized theta-space, scaled to physical units by theta_zt.scale (Pitfall 5); names qualified (Images collision)"
    - "two-part pre-registered decision rule (relative margin AND fold-consistency) — margin alone would pick :aug, consistency guard vetoes fold-unstable gains"
    - "shared ABL_FIX_RESULT const: the fixture ablation runs once, consumed by both SC4 and SC5"

key-files:
  created:
    - spike/npe/ablation.jl
    - spike/npe/summary_choice.md
  modified:
    - spike/test/test_npe.jl
    - spike/data/loader.jl

key-decisions:
  - "assess/rmse qualified as NeuralEstimators.assess / NeuralEstimators.rmse (Images/ImageQualityIndexes also export them under the contract.jl chain) — mirrors the 04-05 decision"
  - "RMSE helpers (_abl_rmse_vector/_abl_theta_scale) re-implemented locally in ablation.jl rather than including benchmark.jl, so the ablation does not drag in BenchmarkTools/resimulate"
  - "ABL_REL_MARGIN/ABL_FOLD_CONSISTENCY guarded-const in ablation.jl (no-op under runtests.jl where test_npe.jl declares them first) so a standalone include still gets the locked values"
  - "reported ablation uses N=200/K=5/epochs=200; SC4/SC5 fixture uses N=24/epochs=30 for a fast deterministic gate"

requirements-completed: [ABL-01, ABL-02]

# Metrics
duration: ~55min
completed: 2026-07-01
---

# Phase 4 Plan 06: Summary-Statistic Ablation Summary

**`ablation.jl` trains the SAME NPE architecture on the minimal (:min, 128-dim) and augmented (:aug, 142-dim) cached summaries under the loader's leak-free k=5 CV with zero re-simulation, scores per-parameter RMSE (all 7 theta, rho_true headline) per fold via assess/rmse, and applies the pre-registered two-part decision rule (D-07) — turning SC4 and SC5 green and choosing `:min`: on the reported N=200 run :aug improves rho_true mean RMSE ~18% (0.099 vs 0.122, passing the 5% relative margin) but wins only 3/5 folds (failing the >=4/5 consistency gate), so the parsimony + OOD-detectability default holds, with the structural OOD blind spot named for Phase 5.**

## Performance

- **Duration:** ~55 min
- **Tasks:** 3
- **Files:** 4 (2 created, 2 modified)

## Accomplishments

- **`ablate(dir; master_seed, K, N, epochs, batchsize)` (Task 1)** — for each fold f in 1:K and each variant in (:min, :aug), trains `build_estimator` via `train_fold` (arch/seed/hyperparams identical; only the cached `variant` differs, D-06) and scores the held-out fold's per-parameter RMSE via `NeuralEstimators.rmse(assess(est, theta_va_std, Zva))` scaled to physical units by `theta_zt.scale` (Pitfall 5). Returns a `K x 2 x 7` `rmse_table` (fold x variant x parameter) plus `rho_rmse_min`/`rho_rmse_aug` (per-variant mean rho_true RMSE), `fold_wins_aug` (folds where :aug beats :min on rho_true), and `rho_per_fold`. CPU-only (`use_gpu=false`).
- **`choose_summary(result)` + `spike/npe/summary_choice.md` (Task 2)** — the D-07 rule returns `:aug` iff `rho_rmse_aug <= (1 - ABL_REL_MARGIN)*rho_rmse_min` AND `fold_wins_aug >= ABL_FOLD_CONSISTENCY`, else `:min`. The docstring notes a SBC-pass-but-high-RMSE outcome is an insufficiency signal, not a pass. `summary_choice.md` records the RMSE table, the rule inputs (margin 0.05, >=4/5 folds), the chosen variant, and an explicit OOD-detectability note naming the structural blind spot (summary-orthogonal misspecifications are undetectable by a fixed patch-correlation summary) — the ABL-02 deliverable Phase 5 consumes.
- **SC4 (ABL-01) + SC5 (ABL-02) green (Task 3)** — SC4 asserts a finite `K x 2 x 7` RMSE table for BOTH variants across k=5 folds with rho_true present per variant and `fold_wins_aug` consistent; SC5 asserts `choose_summary` in {:min,:aug} matching the pre-registered rule applied to the fixture result AND that `summary_choice.md` exists with an OOD note. A shared `ABL_FIX_RESULT` const runs the fixture ablation once. SC1/SC2/SC3 untouched.

## Measured Results (reported run: N=200 main pool, K=5, epochs=200, NPE_MASTER_SEED, CPU-only)

| Fold | :min rho_true RMSE | :aug rho_true RMSE | aug wins |
|------|-------------------|-------------------|----------|
| 1 | 0.1240 | 0.0889 | yes |
| 2 | 0.1059 | 0.1148 | no |
| 3 | 0.1148 | 0.0894 | yes |
| 4 | 0.1444 | 0.0830 | yes |
| 5 | 0.1201 | 0.1207 | no |
| **mean** | **0.1218** | **0.0993** | **3/5** |

Rule: margin threshold (1-0.05)*0.1218 = 0.1158; aug 0.0993 <= 0.1158 (margin PASS) but fold_wins_aug 3 < 4 (consistency FAIL) => **chosen = :min**.

## Task Commits

1. **Task 1: :min vs :aug per-parameter RMSE under k=5 CV** — `9dbbc88` (feat)
2. **Task 2: decision rule + summary-choice justification** — `7570af7` (feat)
3. **Task 3: turn SC4 (ABL-01) + SC5 (ABL-02) green** — `4def573` (test)

## Files Created/Modified

- `spike/npe/ablation.jl` (created) — `ablate`, `choose_summary`, `fold_rmse`, local `_abl_rmse_vector`/`_abl_theta_scale`, guarded `ABL_REL_MARGIN`/`ABL_FOLD_CONSISTENCY` consts.
- `spike/npe/summary_choice.md` (created) — ABL-02 written justification + OOD-detectability note.
- `spike/test/test_npe.jl` (modified) — SC4/SC5 `@test_skip` replaced with live gates; shared `ABL_FIX_RESULT` const + guarded include of ablation.jl; SC1/SC2/SC3 + Wave-0 untouched.
- `spike/data/loader.jl` (modified) — zero-variance z-score guard (Rule-3 deviation, below).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Loader zero-variance z-score produced NaN on the :aug arm**
- **Found during:** Task 1 (first `:aug` training run — `rho_rmse_aug` was NaN).
- **Issue:** `load_fold(...; variant=:aug)` standardizes the 14 augmented moment rows with a `ZScoreTransform` fit on the fold's train columns. The `frac_missing` moment (row 138) is **constant** across the train columns whenever no train stack has a missing patch (std = 0), so the z-score `(x - mean)/scale` divided by zero and yielded NaN, poisoning the whole `:aug` tensor and its RMSE. The `:min` arm was unaffected (its continuous correlation rows vary). A latent loader limitation that only surfaced now that `:aug` is trained for the first time.
- **Fix:** After the `fit(ZScoreTransform, ...)` in `load_fold`, guard the scale vector — any zero or non-finite scale is set to 1, mapping a constant feature to its centered value (0). Leak-free (the fit still sees only train columns); harmless for `:min`; benefits every future `:aug` consumer (Phase 5).
- **Files modified:** spike/data/loader.jl
- **Commit:** 9dbbc88
- **Scope note:** `loader.jl` is a Phase-3 file outside the plan's declared `files_modified`; the change is a minimal, behavior-preserving blocking-issue fix. Full suite (Phase 1-3 incl. SC-3 leak-free standardization) stays green after it.

**Total deviations:** 1 auto-fixed (Rule 3 blocking). No architectural (Rule 4) changes; scope otherwise unchanged.

## Verification Evidence

- **Task 1:** `julia --project=spike -e 'include("spike/npe/ablation.jl"); @assert isdefined(Main,:ablate)'` => "ablate defined OK". Functional smoke on a fixture cache returns a finite `K x 2 x 7` table; CUDA not loaded.
- **Task 2:** `julia --project=spike -e 'include(ablation.jl); @assert isdefined(Main,:choose_summary); @assert isfile("spike/npe/summary_choice.md")'` => OK. Reported ablation (N=200) => choose_summary => :min.
- **Task 3:** `julia --project=spike -t 1 spike/test/test_npe.jl` => 121/121 pass (SC1 68, SC2 15, SC3 9, SC4 12, SC5 6, Wave-0 11); EXIT 0.
- **Full suite:** `julia --project=spike -t 1 spike/test/runtests.jl` => EXIT 0; CPU smoke 11/11, SIM-01/02/03/04 + D-15 green, Phase 3 69/69 (SC-3 leak-free standardization intact), Phase 4 121/121. No CUDA loaded.
- **Decoupling:** `git diff --name-only <base> HEAD` = only spike/data/loader.jl, spike/npe/ablation.jl, spike/npe/summary_choice.md, spike/test/test_npe.jl. 0 src/ edits; 0 spike Project.toml/Manifest.toml changes.

## Known Stubs

None. SC1-SC5 are all live; the scaling characterization (04-07) is the remaining Phase-4 wave, not a stub in this plan's code.

## Next Phase Readiness

- Phase 5 (SBC / amortized BF / OOD) builds on the **`:min`** summary and `trained_npe.jld2`; `summary_choice.md` records the choice and the inherited OOD blind spot (patch-correlation-orthogonal misspecifications undetectable).
- `ablate`/`choose_summary` are re-runnable at production scale; the zero-variance loader guard makes the `:aug` path robust should Phase 5 revisit it.

## Self-Check: PASSED
- Files verified present: spike/npe/ablation.jl, spike/npe/summary_choice.md, spike/test/test_npe.jl, spike/data/loader.jl.
- Commits verified present: 9dbbc88, 7570af7, 4def573.
