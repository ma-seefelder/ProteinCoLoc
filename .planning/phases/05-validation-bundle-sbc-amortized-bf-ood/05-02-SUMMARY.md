---
phase: 05-validation-bundle-sbc-amortized-bf-ood
plan: 02
subsystem: testing
tags: [neuralestimators, ratioestimator, nre, bayes-factor, kde, quadgk, jld2, flux, sbi]

# Dependency graph
requires:
  - phase: 05-01
    provides: shared harness.jl (load_frozen_model, val_rng, draw_simulate_infer_paired) + pre-registered consts.jl (BF_CORR_MIN/BF_LOGBF_TOL/BF_SWEEP_*)
  - phase: 04
    provides: trained_npe.jld2 (frozen NPE + zt/θzt) and the infer.jl read surface (rho_draws, standardize_summary)
  - phase: 04
    provides: isolated spike/baseline/ env pattern (advi_artifact.jld2, run_advi.jl "port not include" precedent)
provides:
  - Amortized one-pass log Bayes factor via a NeuralEstimators RatioEstimator (NRE-as-model-comparison over a binary model index m∈{0,1})
  - measure_log_prior_odds from the training label balance (Pitfall 5 — measured, not assumed 0)
  - A KDE-baseline reproduction pipeline (bf.jl generates Δρ draws → spike/baseline/run_bf_baseline.jl scores the compute_BayesFactor KDE odds-ratio → bf_reproduction compares)
  - test_bf.jl re-runnable BF gate (one-pass finite log-BF, no-kde/quadgk, D-08 locked consts)
affects: [05-04 (reported BF_SWEEP_N=25 run), 05-03 (OOD sibling), 06 (Go/No-Go memo), 07 (productionization)]

# Tech tracking
tech-stack:
  added: [KernelDensity (isolated spike/baseline/ env only), QuadGK (baseline env direct dep)]
  patterns: [NRE model-comparison log-BF = logratio(m=1)-logratio(m=0)-log_prior_odds, cross-env baseline reproduction via persisted Δρ draws]

key-files:
  created:
    - spike/validation/train_ratio.jl
    - spike/validation/bf.jl
    - spike/baseline/run_bf_baseline.jl
  modified:
    - spike/test/test_bf.jl
    - spike/baseline/Project.toml
    - spike/baseline/Manifest.toml
    - .gitignore

key-decisions:
  - "l-POP honestly re-labeled: v0.2.1 RatioEstimator loss is hard-coded logitbinarycrossentropy — used the standard Hermans-2020 NRE model-comparison reduction, which equals the baseline BF by construction"
  - "src/bayes.jl is NOT include()-able in the lean spike env (needs Turing/GLMakie/KernelDensity via the ProteinCoLoc module's imports) — the 3-statement compute_BayesFactor KDE math is PORTED verbatim into the isolated baseline env (run_advi.jl precedent), src/ untouched"
  - "coloc/null split defined on the ρ_true CONTRAST at threshold 0 (≡ baseline {Δμ>0}, ghat strictly monotone) — ghat(0)≈-0.0375 is irrelevant for a two-stack contrast"
  - "log_prior_odds MEASURED = -0.0294 (≈0 by symmetric i.i.d. model-index prior) — measured, never assumed 0"
  - "path (b): KernelDensity+QuadGK added to the ISOLATED spike/baseline/ env only; the shared spike/ env is provably untouched (NeuralEstimators stays 0.2.1)"

patterns-established:
  - "Amortized log-BF in one forward pass: logratio(Z; grid=[0 1]) difference minus measured prior-odds, no quadgk/KDE/shuffle on the amortized side"
  - "Cross-env apples-to-apples reproduction: the spike env persists Δρ posterior+prior draws; the isolated baseline env scores the KDE odds-ratio on those identical draws"

requirements-completed: [BF-01, BF-02]

# Metrics
duration: 19min
completed: 2026-07-02
---

# Phase 5 Plan 02: Amortized Bayes Factor Summary

**Amortized one-pass log Bayes factor via a NeuralEstimators RatioEstimator (NRE model-comparison over a binary Δρ-sign model index), reproducing the src/bayes.jl KDE odds-ratio over a held-out Δρ sweep — correlation 0.961 at fixture scale (magnitude calibration deferred to the reported run).**

## Performance

- **Duration:** 19 min
- **Started:** 2026-07-02T10:37:06Z
- **Completed:** 2026-07-02T10:56:11Z
- **Tasks:** 2 (both TDD auto)
- **Files modified:** 7 (3 created, 4 modified)

## Accomplishments
- `train_ratio.jl`: a `RatioEstimator` (num_parameters=1 = the binary model index) trains CPU-only on 256-dim paired frozen-`m.zt` summaries; `measure_log_prior_odds` reads the coloc prior-odds from the label balance (measured = −0.0294 ≈ 0); atomic `.tmp`+integrity-check+`mv` persistence round-trips the net + zt + measured constant.
- `bf.jl`: `amortized_log_bf` = `logratio(Z;m=1) − logratio(Z;m=0) − log_prior_odds` in ONE forward pass with NO kde/quadgk/shuffle (grep-asserted); `generate_bf_sweep`/`bf_reproduction` drive the D-08 Δρ sweep.
- `run_bf_baseline.jl` (isolated baseline env): the verbatim `compute_BayesFactor` KDE odds-ratio (`kde`→`quadgk`→posterior-odds/prior-odds) on the SAME persisted Δρ draws — the BF-02 reproduction target, reached without editing src/.
- `test_bf.jl`: green re-runnable gate — BF-01 finite one-pass scalar, BF-02 no-kde/quadgk grep, D-08 pre-registered constants locked, and the fixture reproduction (corr ≥ relaxed floor).
- Empirical reproduction at fixture scale: **corr = 0.961 (≥ BF_CORR_MIN 0.95 — D-08a PASS)**; **max|Δ log-BF| = 3.18 (> BF_LOGBF_TOL 0.5 — D-08b NOT met at fixture scale)**.

## Task Commits

1. **Task 1: train_ratio.jl — RatioEstimator model-index net + measured log_prior_odds** - `30cf570` (feat)
2. **Task 2: bf.jl amortized log-BF + KDE-baseline reproduction + test_bf.jl + run_bf_baseline.jl** - `e5a4a5c` (feat)

_Both tasks are `tdd="true"`; each was implemented against its verify snippet / gate before commit._

## Files Created/Modified
- `spike/validation/train_ratio.jl` - RatioEstimator model-index training, measured log_prior_odds, atomic save/load_ratio
- `spike/validation/bf.jl` - amortized_log_bf (one pass), build_bf_pair/generate_bf_sweep/bf_reproduction, save_bf_figure hook into figures.jl
- `spike/baseline/run_bf_baseline.jl` - verbatim compute_BayesFactor KDE odds-ratio in the isolated baseline env (BF-02 target)
- `spike/test/test_bf.jl` - filled BF gate (replaces the Wave-1 @test_skip scaffold)
- `spike/baseline/Project.toml` / `Manifest.toml` - added KernelDensity + QuadGK to the ISOLATED baseline env (path b)
- `.gitignore` - ignore the regenerable BF caches (trained_ratio / bf_sweep_draws / bf_baseline_artifact `.jld2`)

## Decisions Made
- **NRE model-comparison, not literal l-POP** (honesty note carried from research): v0.2.1's `RatioEstimator` loss is hard-coded `logitbinarycrossentropy`; a custom loss is silently ignored. The difference `logratio(m=1)−logratio(m=0)` provably equals the baseline posterior-odds/prior-odds Bayes factor, so it is the correct amortized quantity. Literal l-POP remains a documented fallback only.
- **Split on the ρ_true contrast at threshold 0.** Verified `ghat(0) ≈ −0.0375 ≠ 0`, but the model index is a two-stack CONTRAST and `ghat` is strictly monotone, so `sign(ρ_s − ρ_c) == sign(μ*_s − μ*_c)` — the split reproduces the baseline `{Δμ>0}` event exactly, independent of `ghat(0)`.
- **log_prior_odds measured (−0.0294), not assumed 0** (Pitfall 5). By the symmetric i.i.d. model-index prior it is ≈0, but it is measured from the label balance and subtracted, as D-07 requires.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] src/bayes.jl not include()-able in the lean spike env — ported the KDE math instead**
- **Found during:** Task 2 (bf.jl / run_bf_baseline.jl)
- **Issue:** The plan's `<action>`/acceptance wanted bf.jl to reach `compute_BayesFactor` via a read-only `include(src/bayes.jl)`. But `bayes.jl` has no top-of-file `using`; it relies on the enclosing `ProteinCoLoc` module's `import KernelDensity: kde` / `import QuadGK: quadgk` and `DynamicPPL`/Turing — none of which are in the lean `spike/` env (this is exactly why `contract.jl`'s `induced_mu` re-implements bayes.jl's mean WITHOUT including it). Including it is infeasible in both the main spike env and (cleanly) the baseline env.
- **Fix:** Reproduced the three `compute_BayesFactor` statements (src/bayes.jl:110-135) VERBATIM inside `spike/baseline/run_bf_baseline.jl` (isolated env), mirroring the `run_advi.jl` precedent (which PORTS the src/bayes.jl `@model` rather than including it). Provenance is cited in a header comment. `src/` is never edited.
- **Files modified:** spike/baseline/run_bf_baseline.jl, spike/validation/bf.jl (loads the baseline artifact instead of including bayes.jl)
- **Verification:** `git status --porcelain src/` clean; baseline produced `bf_baseline_artifact.jld2`; reproduction corr computed.
- **Committed in:** e5a4a5c

**2. [Rule 3 - Blocking] KernelDensity absent from the baseline env — added to the ISOLATED baseline Project.toml (path b)**
- **Found during:** Task 2 (run_bf_baseline.jl)
- **Issue:** `KernelDensity` was NOT a top-level node in `spike/baseline/Manifest.toml` (only a FlexiChains weakdep), so `using KernelDensity` failed. QuadGK was present.
- **Fix:** `Pkg.add(["KernelDensity","QuadGK"])` against **`--project=spike/baseline` only** — the sanctioned isolated env (spawner explicitly permits its own Project/Manifest). The SHARED `spike/` env is provably untouched: NeuralEstimators stays 0.2.1 and KernelDensity/QuadGK are NOT in `spike/Project.toml` (verified). This is the default path (b); path (a) — promoting KDE into the shared env — was NOT taken, so the co-resolve landmine is avoided entirely.
- **Files modified:** spike/baseline/Project.toml, spike/baseline/Manifest.toml
- **Verification:** shared-env resolve-risk gate (h) green; `KernelDensity in spike deps? false`, `NeuralEstimators = 0.2.1`.
- **Committed in:** e5a4a5c

---

**Total deviations:** 2 auto-fixed (both Rule 3 blocking). **Impact:** Necessary to reach the BF-02 reproduction target without violating the read-only-src/ and shared-env-isolation hard constraints. No scope creep; no shared-env dependency change.

## Issues Encountered

**Honest reproduction finding (pre-registration discipline — NOT tuned).** At the fixture-scale ratio net (n=1200 pairs, 50 epochs) over a 10-point sweep, the amortized log-BF reproduces the KDE baseline's ORDERING well (**Pearson corr = 0.961 ≥ BF_CORR_MIN = 0.95, D-08a met**) but is MAGNITUDE-INFLATED (amortized ±5.3 vs KDE ±2.6 → **max|Δ log-BF| = 3.18 > BF_LOGBF_TOL = 0.5, D-08b NOT met at fixture scale**). This is the classic under-trained-NRE overconfidence that D-08b's bounded-error criterion exists to catch. The pre-registered constants were NOT changed. The FULL reported reproduction (BF_SWEEP_N=25 at a reported-scale ratio net) is owned by Wave-3 (05-04 run_bf.jl); magnitude calibration is expected to tighten with more training data / a longer sweep, but that is a Wave-3 empirical question, reported honestly here rather than tuned away. The test gate uses a relaxed fixture correlation floor (0.75) and asserts the locked constants separately — it does not assert D-08b at fixture scale.

## User Setup Required
None - CPU-only, offline. The BF caches (`trained_ratio.jld2`, `bf_sweep_draws.jld2`, `bf_baseline_artifact.jld2`) are gitignored and regenerable via `train_ratio.jl` / `bf.jl` + `julia --project=spike/baseline spike/baseline/run_bf_baseline.jl`.

## Known Stubs
None — every function is wired to real data (the frozen NPE, the trained ratio net, and the KDE baseline). The reproduction artifacts are intentionally untracked caches, not stubs.

## Next Phase Readiness
- BF-01/02 machinery complete and green under `runtests.jl`. Ready for Wave-3 (05-04) to run the REPORTED BF_SWEEP_N=25 reproduction at the locked BF_CORR_MIN/BF_LOGBF_TOL and record the headline numbers.
- OPEN for 05-04: D-08b (bounded magnitude) is not met at fixture scale; the reported-scale net must demonstrate whether corr AND magnitude both clear the pre-registered gate. If magnitude remains inflated at reported scale, that is a genuine calibration finding for the Go/No-Go memo (Phase 6), not a code defect.
- `src/` provably untouched; shared `spike/` env unchanged (NeuralEstimators 0.2.1); KDE/QuadGK quarantined in the isolated baseline env.

## Self-Check: PASSED

- Files verified on disk: train_ratio.jl, bf.jl, run_bf_baseline.jl, test_bf.jl, 05-02-SUMMARY.md — all FOUND.
- Commits verified in history: `30cf570` (Task 1), `e5a4a5c` (Task 2) — all FOUND.
- Full `spike/test/runtests.jl` green; shared-env resolve-risk gate green (NeuralEstimators 0.2.1); `src/` clean.

---
*Phase: 05-validation-bundle-sbc-amortized-bf-ood*
*Completed: 2026-07-02*
