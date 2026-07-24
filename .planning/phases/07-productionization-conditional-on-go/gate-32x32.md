# 32×32 — CAP DECISION (not trained, not shipped)

**Decision (2026-07-24, user): CAP the shipped family. 32×32 is NOT trained and NOT shipped.**
This closes plan 07-09 without executing its train-then-gate flow. It supersedes the 07-09-PLAN's
`ship-with-caveat vs cap` checkpoint: the Go/No-Go outcome determined the disposition up front, so
the conditional 32×32 training run (`_train_grid_pipeline(32)`) was correctly **skipped**. No
`artifacts/grid_32/*`, no `test/gate/gate_consts_32.jl`, and no `gate_report_32.jld2` were created.

This is a decision note, not a gate outcome. There is no 32×32 SBC / BF / OOD verdict to record,
because no 32×32 net was trained. What follows is the rationale.

## Why 32×32 is capped

1. **The GO rests on 8×8, not on the fine grids.** Per
   `07-GO-NO-GO-UPDATE.md` (Option A — GO with named limits), v2.0 ships as a calibrated
   colocalization tool on the **8×8** reference grid: TARGETS ρ_true/Δρ SBC-calibrated with
   randomized-rank atom handling, the Bayes factor validated by simulation-based discrimination
   (spike 014, AUC 0.994), and the OOD flag passing (amended gate, pooled AUC 1.0). None of that
   evidence depends on a 32×32 estimator. A 32×32 grid would add shipped surface with no bearing on
   the scientific claim the GO certifies.

2. **A 32×32 gate would run the apparatus now known to be defective.** The per-grid ship-gate
   machinery (`test/gate/run_gate.jl`) scores BF against the **KDE baseline** and scores SBC without
   the **prior-atom / randomized-rank** correction and without the **nuisance-equivalence** rule.
   Spikes 006–014 established that:
   - the KDE Bayes-factor reference is **invalid past |log-BF| ≈ log(L) ≈ 6.9** (BF attrition is
     100 % baseline-side; the NRE never goes non-finite) — so a fresh 32×32 BF arm would fail on the
     same baseline artifact that failed 8×8 (amended, corr `:invalid`, 42 % attrition) and 16×16
     (corr 0.9150, max|Δ| 12.72);
   - the ρ_true SBC KS test rejects on **prior atoms** unless randomized ranks are used, and the
     16×16 gate already showed the headline ρ_true rejected (KS p = 7.0e-4) at a finer grid — a
     genuine grid-resolution effect that gets **worse**, not better, at 32×32 (noisier per-patch
     correlations, longer summary);
   - the fixed patch-correlation summary leaves the nuisances **non-identified** (5/8 columns vacuous
     at 8×8), and a finer grid does not add identifiability.

   A 32×32 gate run through that apparatus would therefore reproduce the 16×16 outcome — an honest
   FAIL that contributes no new information — at substantial CPU cost (32×32 needs ~100k pairs at
   ≥1024² images; data-gen is the CPU-bound bottleneck, per 07-RESEARCH §Compute cost reality).

3. **User decision: cap the family.** Rather than ship a calibrated-but-uninformative or
   gate-failing estimator, the user elected to cap the shipped family. 32×32 (and 16×16) are
   excluded; only 8×8 ships (`_SHIPPED_GRIDS = (8,)`). This is the "cap-the-family" branch the
   07-09-PLAN itself names as an expected, non-failure outcome (Task 3, `cap-the-family`), and it is
   consistent with the 16×16 disposition already recorded (`gate-16x16.md`: NOT eligible for default
   registry population).

## Scope confirmations

- **64×64 remains dropped** (D-04), as it was before this plan.
- **16×16** is excluded from `_SHIPPED_GRIDS` (its gate FAILED with no non-method cause for the
  ρ_true rejection; `gate-16x16.md`).
- **4×4** is excluded from the shipped registry: it was never post-hoc re-analysed with the corrected
  statistics (randomized ranks / nuisance equivalence / simulation-based BF) on which the GO rests,
  so its literal-FAIL gate (`gate-4x4.md`) is not superseded by any calibrated re-reading. The GO
  certifies 8×8 only.
- **Only 8×8 ships.** `_SHIPPED_GRIDS = (8,)` in `src/registry.jl` (07-10), with a comment citing
  this decision and `07-GO-NO-GO-UPDATE.md`.

## Reproducibility / anti-snooping

No pre-registered seed was consumed. `PROD_SEED[32]` / `PROD_SEED_V2[32]` were never drawn against
a trained net. No gate was run. The existing frozen artifacts (`amended_v2`, `grid_4/8/16`,
`spike01x`) are byte-untouched.

**Reference:** `07-GO-NO-GO-UPDATE.md` (Option A), `gate-16x16.md`, `gate-8x8-amended.md`,
`Skill("spike-findings-proteincoloc")`.
