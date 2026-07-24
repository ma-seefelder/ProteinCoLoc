# Go/No-Go — Updated Verdict (supersedes the Phase-6 "Clean Go")

**Date:** 2026-07-24. This updates `06-GO-NO-GO-MEMO.md` after the amended grid-8 ship-gate ran and
FAILED, and after the calibration investigation (spikes 006-014) diagnosed every residual. The
Phase-6 memo's "Clean Go" rested on arguments this investigation has overturned; this document
states the honest current position and the decision it forces.

---

## 1. What the Phase-6 "Clean Go" assumed — and what fell

The Phase-6 memo justified a Go by reading each residual gate failure to a non-model cause
(chi²-over-power, clamped-KDE tail artifact, data-scale gap), backed by a promise: an independent
fresh-seed re-pre-registered confirmation ship-gate inside Phase 7. That gate ran. Here is what
changed:

| Phase-6 claim | Status now |
|---|---|
| "A fresh pre-registered gate will confirm the Go" | **The amended gate FAILED** (one-shot, PROD_SEED_V2): SBC ks_pass=false, BF corr_verdict=:invalid, OOD pass. |
| "SBC failures are chi²-over-power at M=2000" | **Partly true, but the multiplicity correction did NOT rescue SBC** — Holm-corrected KS still rejects (Spike 011). Over-power is real for nuisances only. |
| "max\|Δ logBF\| is a clamped-KDE tail artifact" | **INCOMPLETE** — the clamp MASKED a genuine order-of-magnitude NRE-vs-KDE tail disagreement (Spike 008). The KDE baseline was the wrong reference, not a mere artifact. |
| "The model is calibrated" (with ρ_true/Δρ underconfident) | **The coloc targets ARE calibrated; the "overconfidence" was a prior-ATOM artifact** (Spike 010 retraction). But no grid ever passed a gate literally. |

Two Phase-6 pillars — the multiplicity defence and the clamp-artifact defence — do not hold. The
Go can no longer be "clean." The question is whether it can be **honest**.

## 2. The honest current position, with the CORRECT methods

The investigation did not just diagnose; it established the methods the Phase-6 gate lacked. With
those methods applied:

**SBC (post-hoc re-analysis, `posthoc_reanalysis.jl`, independent DEV seed, M=2000):**
- **Targets calibrated.** ρ_true: raw KS 5.9e-4 → **randomized-rank KS 0.72** (the atoms are a
  mixed-distribution artifact; randomized ranks are the textbook fix). Δρ: KS 0.20. Holm over the
  targets: **PASS**.
- **Nuisances: a named limit.** Under a pre-registered δ=0.10-SD equivalence test, 4/6 pass; 2/6
  (autofluorescence, label_efficiency) fail — NOT because the drift is large (point drift ~0.08 SD
  < 0.10) but because at M=2000 the 90% CI edges just past 0.10. No model lever removes it (Spike
  012 capacity redistributes; Spike 013 truncation fixes only the atoms at a real cost). This is a
  marginal-reproduction limit on parameters the summary cannot constrain, not a failure of the
  coloc inference.

**BF (Spike 014, the replaced method):** the KDE baseline is abandoned. Validated the amortized
NRE log-BF SBC-style — by discrimination + monotonicity + decision-calibration, no per-pair
reference:
- **AUC(coloc vs null) = 0.994**, 0 attrition (vs the KDE baseline's 35%).
- Monotone in evidence (Spearman 0.93), decision calibration ECE 0.019, FPR 0.036 at logBF>0.
- This is a genuine methodological upgrade: the BF axis moves from "defective reference" to
  "excellent, calibrated discrimination." (Concept-proof, single seed — full integration is a §6
  step; see §4.)

**OOD (amended gate):** pooled AUC 1.0, ID fire rate 0.05, all four misspecification families
detected. Passes. (Known blind spots documented: correlation-preserving affine transforms.)

**Net:** on the axes that carry the scientific claim — calibrated coloc posteriors, a
discriminating Bayes factor, an honest OOD flag — v2.0 delivers, once the correct statistics are
used. What remains is a bounded, named limit on nuisance-parameter marginal accuracy.

## 3. The §6.4 dead-end (why this is a decision, not more iteration)

Every residual is **test-side fixable** (randomized ranks, nuisance equivalence, simulation-based
BF). But test-side fixes are not a model change, and the first amendment's §6.4 forbids a further
pre-registered gate run without one. The only model change that fits — μ-prior truncation — was
tested (Spike 013) and **rejected by the project** because it breaks ADVI comparability and
degrades |ρ|≥0.95 inference, while randomized ranks achieve the same SBC validity for free.

So there is no path to a fresh "pre-registered gate passed" claim on grid 8. The honest options are
bounded, and this is what makes it a Go/No-Go decision rather than another spike:

## 4. The decision

**Option A — GO with named limits (publish the honest post-hoc picture).**
Claim: v2.0 delivers calibrated colocalization inference (ρ_true, Δρ) verified by SBC with correct
atom handling, a Bayes factor that discriminates coloc from null at AUC 0.994, and an honest OOD
flag. Named limits, stated in the manuscript: (i) nuisance-parameter marginals carry a ~0.08-SD
residual drift not certifiable as negligible at M=2000; (ii) the ρ_true SBC uses randomized ranks
for the prior atoms; (iii) the BF is validated by simulation-based discrimination, not a per-pair
reference; (iv) the gate was amended twice — pre-registration credibility on grid 8 is spent, and
the memo says so. No further training. Fastest honest ship.

**Option B — MORE WORK (raise the calibration bar before ship).**
Pursue the two open roots: (i) integrate the simulation-based BF into the gate as a §6
pre-registration step (the concept is proven; this is the third gate change and further spends
§6.4 credibility); (ii) attack the nuisance drift at its source — the fixed 8×8 patch-correlation
summary cannot constrain the nuisances, so this likely means a **summary-design change** (a richer
or learned summary), which is a substantial re-architecture with its own retrain + revalidation.
Higher calibration ceiling, materially more effort, and — critically — no evidence yet that the
nuisance drift is removable at all (both model levers failed).

**Recommendation:** Option A. The coloc claim is sound; the residual is a bounded nuisance limit
that no tested lever removes, so Option B's cost is high and its payoff on the binding limit is
unproven. The BF concept-proof (Spike 014) already lets the manuscript present a strong, honest BF
story without a fresh gate. Option B is justified only if a reviewer/venue requires all-parameter
SBC — a bar the summary design may not physically support without becoming a different method.

## 5. What either option must carry into the manuscript (non-negotiable honesty)

- The gate was amended twice; the final grid-8 verdict is post-hoc, not a clean pre-registration.
- ρ_true calibration depends on randomized-rank atom handling; the underlying prior places mass on
  non-realizable μ (a documented simulator/prior inconsistency).
- The BF is validated by discrimination/calibration, not magnitude agreement with an independent
  per-pair method (the KDE reference was shown invalid past |logBF|≈log(L)).
- The one-shot gate carries ±1.5-z per-parameter seed variance; single-seed verdicts are fragile.
- Nuisance marginals are a named calibration limit.

Findings and scripts: `Skill("spike-findings-proteincoloc")`, `posthoc_reanalysis.jl`,
`.planning/spikes/006-014`.
