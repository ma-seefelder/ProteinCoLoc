---
phase: 05-validation-bundle-sbc-amortized-bf-ood
plan: 04
subsystem: testing
tags: [sbc, calibration, bayes-factor, nre, ood, misspecification, roc, reported-run, pre-registration, honest-negative]

# Dependency graph
requires:
  - phase: 05-01
    provides: harness.jl (load_frozen_model/val_rng/draw_simulate_infer), sbc.jl, figures.jl, LOCKED consts.jl
  - phase: 05-02
    provides: bf.jl (amortized_log_bf/bf_reproduction/generate_bf_sweep), train_ratio.jl (load_ratio), isolated spike/baseline/ KDE reproduction
  - phase: 05-03
    provides: ood.jl (fit_ood_nulls/ood_roc_over_grid/misspec_* / negctrl_* / verify_summary_invariance)
  - phase: 04
    provides: frozen trained_npe.jld2 (estimator + zt + θzt), infer.jl amortized read surface
provides:
  - "REPORTED SBC run at locked M=2000/L=999/bins=50 on VAL_MASTER_SEED → sbc_report.jld2 + per-parameter figures; hard-gated (exits nonzero on fail)"
  - "REPORTED amortized-BF reproduction over the full BF_SWEEP_N=25 Δρ sweep → bf_report.jld2 + agreement figure; hard-gated on corr/error"
  - "REPORTED OOD ROC over the full 4-family × OOD_GRID_LEVELS grid → ood_report.jld2 + ROC figure; hard-gated on AUC + negative-control KS/quiet"
  - "Reported-scale ratio net (n=3000/epochs=100) regenerated to trained_ratio.jld2"
  - "THREE honest-negative reported findings (all pre-registered gates FAIL) feeding the Phase-6 Go/No-Go memo"
affects: [phase-06-go-no-go-memo, phase-07-productionization]

# Tech tracking
tech-stack:
  added: []  # no new dependency; shared spike/ env untouched (NeuralEstimators 0.2.1)
  patterns:
    - "Reported run = script OWNS its artifact + figures + a REAL pre-registered gate (exits nonzero on fail); numbers recorded verbatim, consts never tuned"
    - "Cross-env BF reproduction driven from the reported script: generate sweep (spike) → shell out to isolated baseline env → score reproduction"

key-files:
  created:
    - spike/validation/run_sbc.jl
    - spike/validation/run_bf.jl
    - spike/validation/run_ood.jl
  modified:
    - .gitignore

key-decisions:
  - "ALL THREE reported gates FAIL at locked constants on the fresh VAL_MASTER_SEED stream — recorded verbatim as honest negatives; consts.jl provably UNCHANGED (no tune-until-pass)"
  - "Reported-scale ratio net (n=3000/epochs=100, train_ratio.jl's declared Wave-3 training) regenerated so the BF headline is not confounded by the 05-02 fixture net (n=1200); BF still fails, robust to training scale"
  - "PP (posterior-predictive) OOD channel is REPORTED NON-VIABLE on this frozen net: non-finite posterior θ̂ on strong misspec crashes pp re-simulation (same NPE overconfidence SBC exposed); excluded from OR-fusion (density channel is the reported detector), non-viability evidenced by a θ̂-finiteness probe. Root-cause fix belongs in ood.jl (out-of-scope per file-isolation) — surfaced to Phase-6"
  - "OOD nulls fit on FRESH in-distribution draws through the frozen chain (id_summary, m.zt space) — train-distribution only, never eval/misspec (SC5/D-06), consistent with the ROC's own ID negatives"

requirements-completed: [SBC-01, SBC-02, SBC-03, SBC-04, BF-01, BF-02, OOD-01, OOD-02]

# Metrics
duration: 60min
completed: 2026-07-02
---

# Phase 5 Plan 04: Reported Trifecta (SBC / amortized-BF / OOD) Summary

**The three pre-registered proofs were executed at the LOCKED constants on the fresh, never-tuned-against VAL_MASTER_SEED stream — and ALL THREE gates HONESTLY FAIL: the frozen Phase-4 NPE is overconfident/miscalibrated (SBC KS/χ² p≈0 for every parameter), the amortized log-BF is magnitude-inflated (corr 0.886 < 0.95, max|Δ| 3.60 > 0.5), and the density OOD detector clears only 3 of 4 families (noise blind spot → pooled AUC 0.730 < 0.80). consts.jl is provably untouched; these are negative findings for the Phase-6 Go/No-Go memo, not tuned away.**

## HONESTY CONTRACT COMPLIANCE

- **consts.jl UNCHANGED** by all three run tasks (`git diff` clean on `spike/validation/consts.jl`) — no threshold weakened, no VAL_MASTER_SEED altered, no seed cherry-picked.
- **Reserved stream verified disjoint**: `VAL_MASTER_SEED = 0x5BC0FFEE ≠ NPE_MASTER_SEED = 0xC0FFEE` (training) and `≠ VAL_FIX_SEED = 0xF1F7ED` (fixtures); asserted at the top of every reported script.
- **Every script is a REAL gate**: each exits NONZERO when its pre-registered thresholds fail (verified: all three exited 1).
- **No retry-to-pass**: gates were run ONCE at the locked constants; the only re-run was BF after regenerating the *reported-scale* ratio net (train_ratio.jl's own declared Wave-3 training), and it still failed — reported transparently below.

## REPORTED GATE RESULTS (actual numbers, verbatim)

### SBC (SBC-01..04) — FAIL

Reported run: M=2000 draws × L=999 posterior samples, bins=50, imsize (256,256); rank table computed in 5.07 min on VAL_MASTER_SEED. Artifact: `sbc_report.jld2`; figures: per-parameter rank-hist / coverage / reliability under `spike/validation/figures/`.

| parameter | KS p | χ² p | ECE | verdict | gate (ks>0.05 ∧ χ²>0.05 ∧ ECE<0.05) |
|-----------|------|------|-----|---------|------|
| ρ_true | 2.81e-40 | 5.27e-146 | 0.189 | red | **FAIL** |
| spillover | 3.30e-20 | 1.54e-50 | 0.0459 | green | **FAIL** (KS/χ²) |
| autofluorescence | 2.42e-15 | 3.51e-38 | 0.0469 | green | **FAIL** (KS/χ²) |
| label_efficiency | 1.35e-6 | 7.54e-36 | 0.0302 | green | **FAIL** (KS/χ²) |
| shift_dx | 1.03e-9 | 2.77e-33 | 0.0499 | green | **FAIL** (KS/χ²) |
| shift_dy | 3.80e-12 | 1.89e-49 | 0.068 | yellow | **FAIL** |
| noise | 2.81e-4 | ~0 (2.96e-370) | 0.0408 | green | **FAIL** (KS/χ²) |
| Δρ | 9.49e-40 | 1.45e-147 | 0.192 | red | **FAIL** |

**Verdict: FAIL — every one of the 7 θ parameters AND the paired Δρ fails rank-uniformity** (KS and χ² p-values ≪ 0.05). ρ_true and Δρ additionally show red ECE (~0.19). The rank distributions are strongly non-uniform → the amortized posterior is **overconfident/miscalibrated** (too-narrow credible intervals pile PIT mass at the extremes). Caption stands: *"calibrated under the simulator; pair with the OOD result"* — here the honest reading is **NOT calibrated under the simulator at reported scale**.
*(Note: the reliability-diagram MCE reads 0.99 uniformly — a binning artifact of 19 discrete predicted coverage levels into 50 bins; ECE is count-weighted and unaffected, and the KS/χ² rank tests are the decisive, artifact-free signal.)*

### Amortized Bayes Factor (BF-01, BF-02) — FAIL

Reported sweep: BF_SWEEP_N=25 Δρ points in [−0.6, 0.8] on VAL_MASTER_SEED at imsize (256,256); amortized side `logratio`-only (no kde/quadgk), KDE baseline computed in the isolated `spike/baseline/` env on identical Δρ draws (apples-to-apples). Measured `log_prior_odds = −0.0294`. Artifact: `bf_report.jld2`; figure: `bf_agreement.png`.

| ratio net | corr(amortized, KDE) | max\|Δ log-BF\| | D-08a (corr≥0.95) | D-08b (max\|Δ\|≤0.5) |
|-----------|----------------------|-----------------|-------------------|----------------------|
| reported-scale (n=3000, ep=100) | **0.8861** | **3.598** | **FAIL** | **FAIL** |
| fixture-scale (n=1200, ep=50) — first run | 0.8947 | 4.028 | FAIL | FAIL |

**Verdict: FAIL — both pre-registered criteria fail, robust to training scale.** The amortized log-BF is order-*roughly*-correct but **magnitude-inflated** (amortized range ±5.5/6.5 vs KDE ±3.9/3.2) and correlation dips below 0.95. Consistent with the SBC overconfidence finding and 05-02's fixture-scale prediction. The reported-scale training (train_ratio.jl's declared Wave-3 training) did NOT rescue it.

### OOD ROC (OOD-01, OOD-02) — FAIL

Reported grid: 4 families × OOD_GRID_LEVELS=4 levels, n_id=100, n_pos=40 on VAL_MASTER_SEED at imsize (256,256); nulls fit train-only (fresh ID draws, frozen m.zt space); pre-registered ID-quantile operating point. Artifact: `ood_report.jld2`; figure: `ood_roc_combined.png`.

Summary-density (Mahalanobis) AUC per family (best over levels):

| family | L1 | L2 | L3 | L4 | best | gate (≥0.80) |
|--------|----|----|----|----|------|------|
| texture | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | **PASS** |
| optics | 0.886 | 0.908 | 0.946 | 0.983 | 0.983 | **PASS** |
| background | 0.886 | 0.96 | 0.0 | 0.938 | 0.96 | **PASS** |
| noise | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | **FAIL** (named blind spot) |

- **Density pooled strongest-level AUC = 0.730 < OOD_AUC_MIN 0.80 → FAIL** (dragged down by the noise family = 0.0).
- **Held-out ID fire-rate = 0.05 = 1 − OOD_ID_QUANTILE exactly** → the pre-registered operating point is honest (out-of-sample ~5% FPR). Post-hoc Youden-J j=0.725 (reference only, never the gate).

Summary-orthogonal negative control (D-04 named blind spot):

| transform | max KS | KS<0.05 | density-flag fire-rate | quiet (≤0.05) |
|-----------|--------|---------|------------------------|---------------|
| affine | 0.0156 | yes | 0.0333 | **yes** (provably blind) |
| rotate | 0.0156 | yes | 0.933 | **NO** |
| block | 0.0 | yes | 0.933 | **NO** |

**Verdict: FAIL.** The gate fails on (a) noise-family AUC = 0 (the documented detector-noise blind spot — corruption decorrelates patches toward a central summary), (b) pooled AUC 0.730 < 0.80, and (c) the rotate/block negative controls: all three are **KS-invariant on the summary value distribution** (< 0.05, as designed), but the **position-aware Mahalanobis density flag fires ~93%** on rotate/block — the exact "residual spatial sensitivity" 05-03 predicted (the density null encodes which value sits at which patch position; a spatial permutation moves them). **affine — the only transform that preserves the summary VECTOR bit-identically — is the provable blind spot and correctly stays quiet.** Blind spot MEASURED and NAMED, not hidden.

### PP (posterior-predictive) OOD channel — REPORTED NON-VIABLE

On strongly-OOD inputs the frozen NPE returns a **non-finite posterior mean θ̂** (the flow extrapolates to NaN out of distribution — the same overconfidence SBC exposed), which crashes `pp_mismatch_score` inside `simulate_pair`. The PP channel is therefore **excluded from the reported OR-fusion** (`with_pp=false`); the summary-density channel is the reported detector. Non-viability is evidenced by a θ̂-finiteness probe persisted in `ood_report.jld2`. The root-cause fix (guard non-finite θ̂ in `ood.jl`) is out-of-scope here per the shared-file isolation contract — surfaced to Phase-6 / gap-closure.

## Task Commits

1. **Task 2: run_bf.jl — reported amortized-BF reproduction gate** — `db38e5b` (feat) *(committed first; SBC ran in background)*
2. **Task 1: run_sbc.jl — reported SBC calibration gate at locked M/L** — `b4bbf48` (feat)
3. **Task 3: run_ood.jl — reported OOD ROC gate over the full grid** — `f5bc855` (feat)

## Files Created/Modified

- `spike/validation/run_sbc.jl` (created) — reported SBC at locked M/L, per-parameter KS/χ²/coverage/ECE, atomic `sbc_report.jld2`, figures, hard gate.
- `spike/validation/run_bf.jl` (created) — reported full-sweep amortized-BF reproduction; shells out to the isolated baseline env; atomic `bf_report.jld2`, figure, hard gate.
- `spike/validation/run_ood.jl` (created) — reported density-channel ROC over the full grid + negative control; PP-viability probe; atomic `ood_report.jld2`, figure, hard gate.
- `.gitignore` (modified) — ignore the regenerable reported `*_report.jld2` (mirrors the existing 05-02 artifact pattern; committed deliverables are the scripts).
- `spike/validation/trained_ratio.jld2` (regenerated, gitignored) — reported-scale ratio net (n=3000/epochs=100).

## Deviations from Plan

### Auto-fixed / constraint-driven adjustments

**1. [Rule 1 - Bug, worked around under file-isolation] PP OOD channel crashes on non-finite posterior θ̂**
- **Found during:** Task 3 (first full-grid OOD run, with_pp=true).
- **Issue:** On strong misspec inputs the frozen NPE posterior draws contain NaN/Inf → `θ̂ = mean(...)` is NaN → `_theta_tuple` NaN → `simulate_pair` throws "all θ fields must be finite". A genuine Rule-1 bug in `ood.jl`'s `pp_mismatch_score` (no non-finite-θ̂ guard).
- **Fix:** The shared-file isolation contract forbids editing `ood.jl`, so the fix is worked around in `run_ood.jl`: run the density channel over the full grid (`with_pp=false`, crash-free), exclude PP from the reported OR-fusion, and EVIDENCE the non-viability with a θ̂-finiteness probe. Root-cause fix (guard non-finite θ̂ in `ood.jl`) surfaced to Phase-6.
- **Files modified:** spike/validation/run_ood.jl only.
- **Committed in:** f5bc855.

**2. [Rule 3 - Blocking, honest completeness] Regenerated the ratio net at reported scale before the reported BF run**
- **Found during:** Task 2 (BF).
- **Issue:** The persisted `trained_ratio.jld2` was the 05-02 *fixture-scale* net (n=1200/epochs=50). Reporting the headline BF number on a toy net would understate the method; train_ratio.jl's own banner reserves the reported-scale training (n=3000/epochs=100) for Wave-3 (05-04).
- **Fix:** Ran `julia --project=spike spike/validation/train_ratio.jl` ONCE (the pre-declared reported hyperparameters — NOT a tuning sweep) to produce the reported-scale net, then re-ran run_bf.jl. Result still FAILS (corr 0.886, max|Δ| 3.60). Both numbers reported transparently. consts.jl untouched.
- **Files modified:** spike/validation/trained_ratio.jld2 (gitignored, regenerable).
- **Committed in:** n/a (gitignored artifact); the script is db38e5b.

**Total deviations:** 2 (1 Rule-1 bug worked around under isolation, 1 Rule-3 completeness). **No pre-registered constant in consts.jl was changed. No threshold weakened. No seed cherry-picked.**

## Honest Findings (the whole point of this plan)

All three legs of the trifecta produce **HONEST NEGATIVE results** at the pre-registered constants on the fresh reserved stream. They are mutually consistent and point to a single root cause: **the frozen Phase-4 amortized NPE/NRE is overconfident / miscalibrated at reported (256×256, M=2000) scale.**

1. **SBC:** every parameter's rank distribution is far from uniform (KS/χ² p ≈ 1e-6 … 1e-146); ρ_true & Δρ ECE red (~0.19). The posterior is systematically too narrow.
2. **BF:** the amortized log-BF is magnitude-inflated (max|Δ| ~3.6 vs a 0.5 tolerance) and correlation 0.886 < 0.95 — the ratio net is overconfident in the same direction, and reported-scale training does not fix it.
3. **OOD:** the density detector works on correlation-altering families (texture 1.0, optics 0.983, background 0.96) but is BLIND to detector-noise (AUC 0), and the PP channel is non-viable (NaN θ̂ out-of-distribution — a direct consequence of the same overconfidence). The named summary-orthogonal blind spot is confirmed: affine (vector-preserving) is provably blind & quiet; rotate/block preserve only the value multiset, so the position-aware density flag correctly fires.

These feed the **Phase-6 Go/No-Go memo** as a clear signal: the spike's inference machinery is complete and honest, but the frozen net does not yet meet the pre-registered calibration/reproduction/detection bars — a No-Go-or-iterate input, exactly what the pre-registration discipline exists to surface.

## Known Stubs

None. Every reported script is wired to the real frozen net + real trained ratio net + real KDE baseline. The `*_report.jld2` and PNGs are intentionally gitignored regenerable outputs (not stubs); the numbers are recorded verbatim above.

## Threat Flags

None new. All surface is offline/CPU-only/local-artifact; the pre-registration anti-snooping boundary (T-05-11/12) was honored — consts.jl unchanged, reserved stream disjoint, gates real.

## Next Phase Readiness

- **Phase 6 (Go/No-Go memo):** consumes the three honest-negative reported findings above. The recommended framing is machinery-complete + net-not-yet-calibrated; candidate gap-closure levers (all Phase-6/7, NOT tuned here): retrain/regularize the NPE for calibration (temperature/flow capacity/more data), guard non-finite θ̂ in `ood.jl` to restore the PP channel, and revisit the fixed 8×8 summary's noise blind spot.
- **consts.jl** remains the locked pre-registration; any future run re-uses it unchanged.
- `src/` provably untouched; shared `spike/` env unchanged (NeuralEstimators 0.2.1); KDE/QuadGK remain quarantined in the isolated baseline env.

## Self-Check: PASSED

- Files verified on disk: run_sbc.jl, run_bf.jl, run_ood.jl (FOUND); sbc_report.jld2, bf_report.jld2, ood_report.jld2 (FOUND); figures/sbc_*_*.png, bf_agreement.png, ood_roc_combined.png (FOUND).
- Commits verified in history: db38e5b (run_bf.jl), b4bbf48 (run_sbc.jl), f5bc855 (run_ood.jl) — all FOUND.
- consts.jl `git diff` clean (thresholds never tuned); all three reported scripts exit nonzero (real gates); VAL_MASTER_SEED disjoint from NPE_MASTER_SEED and VAL_FIX_SEED.

---
*Phase: 05-validation-bundle-sbc-amortized-bf-ood*
*Completed: 2026-07-02*
