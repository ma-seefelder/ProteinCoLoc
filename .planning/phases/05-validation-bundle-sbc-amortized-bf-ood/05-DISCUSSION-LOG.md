# Phase 5: Validation Bundle (SBC + Amortized BF + OOD) - Discussion Log

> **Audit trail only.** Not consumed by downstream agents (researcher/planner/executor).
> Decisions captured in `05-CONTEXT.md` — this log preserves the reasoning.

**Date:** 2026-07-02
**Phase:** 05-validation-bundle-sbc-amortized-bf-ood
**Mode:** discuss (interactive)
**Areas discussed:** SBC coverage scope, Misspecification grid, OOD score & threshold, BF null + agreement

## Areas Selected

User selected all four presented gray areas to lock (rest taken as research/discretion under the
Phase-4 pre-registered discipline).

## Area 1 — SBC coverage scope

- **Question:** Which quantities get a calibration proof?
- **Options:** (a) 7 θ + Δρ contrast [recommended]; (b) 7 per-stack θ only; (c) ρ_true + Δρ only.
- **Selected:** (a) 7 θ + Δρ contrast.
- **Rationale:** Δρ is the quantity the BF and coloc verdict rest on; a dedicated paired-draw Δρ SBC
  demonstrates its calibration directly rather than inheriting it by argument from ρ_true.
  → **D-01**.

## Area 2 — Misspecification grid

- **2a — Positive controls (multiSelect):** Selected ALL four — texture-model mismatch,
  noise-model mismatch, optics/PSF mismatch, background/illumination. → **D-03**.
- **2b — Summary-orthogonal negative control:** Options (a) correlation-preserving transform
  [recommended]; (b) phase-scramble in orthogonal band; (c) empirically-selected null set.
  **Selected (a)** — affine intensity / rotation-flip / distribution-preserving rearrangement, with
  empirical KS verification that the 8×8 summary is unchanged; blind spot demonstrated by
  construction. → **D-04**.

## Area 3 — OOD score & threshold

- **3a — Channel fusion:** Options (a) report both + OR-fire [recommended]; (b) single fused score;
  (c) summary-density primary, PP secondary. **Selected (a)**. → **D-05**.
- **3b — Threshold discipline:** User answered "Perform 1-3" → combine all three: threshold-free
  ROC/AUC headline (opt 2) + pre-registered in-distribution quantile as the reported operating point
  (opt 1) + Youden-J shown as a descriptive post-hoc reference only, never gating (opt 3). Nulls fit
  on TRAIN only; ROC on fully external misspecified test images (SC5). → **D-06**.

## Area 4 — BF null + agreement

- **4a — Null definition:** Options (a) match existing baseline exactly [recommended]; (b)
  interval/ROPE null; (c) point null Δρ=0. **Selected (a)** — read `compute_BayesFactor()`
  (`src/bayes.jl:109`, `ρ_threshold=0.0`) and mirror its Δρ null + prior exactly so BF-02 is a true
  reproduction. → **D-07**.
- **4b — Agreement criterion:** Options (a) corr + bounded log-BF error [recommended]; (b)
  decision-agreement rate; (c) correlation only. **Selected (a)** — pre-register both high
  correlation AND bounded absolute log-BF error over the decision-relevant Δρ range; analogous to
  Phase-4's pre-registered RMSE tolerance. → **D-08**.

## Carried Forward (not re-decided)

- M ≈ 2000, pre-registered M + thresholds, fresh independently-seeded held-out SBC run; SBC reported
  "under the simulator" paired with OOD (SBC-03/04) → **D-02**, following the Phase-4
  pre-registered-consts + keyed holdout-repro-gate discipline.
- `:min` summary (Phase-4 D-07); one NPE, Δρ by differencing two passes (Phase-4 D-03); all
  preprocessing frozen from the loader train split (Phase-3 D-07/D-08; SC5).

## Deferred Ideas

- Demo + Go/No-Go memo (Phase 6); productionization / `src/` edits (Phase 7); external physical
  corpus (Phase 8); cross-method comparator (Phase 9); three-hypothesis BF (Phase 13); adversarial
  nuisance sweep + CI gate (Phase 15).

## Claude's Discretion

- Exact M/L; KS/χ² and ECE/MCE cutoffs; density-vs-flow choice + PP statistic; per-family
  perturbation magnitudes + transform ε; exact ID quantile, BF corr threshold, log-BF tolerance, Δρ
  range (all pre-registered before the reported run); shared-harness module layout; CairoMakie figure
  set.
