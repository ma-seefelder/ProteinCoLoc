---
status: specification-note
kind: bar-derivation
phase: 11-registration-and-chromatic-uncertainty-as-latent
written: 2026-07-27
applies_to: SC1g / P11_LAMBDA_ABLATION_FACTOR
code_changes: none
likely_consumer: phase 15 (Calibration Operating Envelope and CI Gate)
evidence: 11-DIAGNOSIS.md, 11-CLOSURE.md
---

# Corrected derivation of the SC1g bar — specification note

## This note deliberately changes no code

`spike/validation/p11_consts.jl` is **byte-unchanged** and stays that way.
`P11_LAMBDA_ABLATION_FACTOR = 2.501971634054976` **keeps its value permanently**, as the historical
record of what was pre-registered.

That is not squeamishness. Phase 11 closes with the record *"not one threshold was edited,
`P11_ITERATION_ALLOWANCE` unspent, `p11_consts.jl` byte-unchanged"*, and that record is what lets
the manuscript say the gate was **diagnosed rather than adjusted**. Editing the very constant that
failed — even to correct a real defect — would compromise exactly that claim and invite the obvious
question. The correction belongs in prose, where a future phase can adopt it knowingly.

## 1. The defect, precisely

`lambda_ratio = 5.003943268109952` is the ratio, across the λ ladder, of the **registration-induced
component** of Δρ — the D-06 probe's `drho_eq` curve evaluated at λ_max over λ_min. The bar was set
at half of it (`SC2_SPEARMAN_ATTENUATION = 0.5`), giving 2.502.

The tripwire then applies that bar to the **total Δρ posterior standard deviation**, which contains
the registration component *plus every other source of uncertainty*. **A bar derived from a
component ratio is being applied to a total.** Those are not the same quantity, and the mismatch is
not a matter of degree — it makes the bar unreachable.

**The arithmetic that shows it was unreachable** (measured, `p11_hybrid_report.jld2`):

| quantity | λ_min = 0.25 | λ_max = 3.0 |
|---|---|---|
| net Δρ posterior SD | 0.09088 | 0.09282 |
| σ_reg (marginalised over the shift prior) | 0.00915 | 0.04595 |
| total in quadrature | 0.09134 | 0.10357 |

Total ratio = **1.134**. To instead reach 2.502, the λ_max total would have to be
2.502 × 0.09134 = 0.22853, requiring

> σ_reg = √(0.22853² − 0.09282²) = **0.2088**

against a measured **0.0460** — i.e. **4.5× larger than the forward model says it is**. There is no
estimator, and no summary, that makes a subdominant variance component drive a 2.5× change in a
total it does not dominate. **The gate was unreachable by construction, independent of the network.**

*(An earlier RMS point-estimate of σ_reg gave 0.042 and a total ratio of 1.115, requiring ~5×. The
marginalised figures above supersede it; the conclusion is identical.)*

## 2. The two correct constructions

### (a) Gate the registration component directly

Statistic and bar both live on σ_reg. `lambda_ratio` is then the *right* source for the bar, because
it is a ratio of exactly the quantity being gated.

- **For:** the bar is already derived correctly; no re-derivation needed. It isolates the physical
  effect from everything else, so the gate answers "does registration sensitivity behave as the
  forward model predicts" — a clean, well-posed question.
- **Against:** σ_reg is not directly observable from a trained posterior. It has to be *constructed*
  (as this phase did, by propagating the probe's `drho_eq` curve), so the gate would test the
  construction as much as the network. It also does not answer the question the phase actually
  cares about, which is about the posterior the user receives.

### (b) Gate the total width, with the bar derived from the measured σ_reg

Keep the statistic (total Δρ posterior SD) and re-derive the bar from the measured σ_reg against the
measured total, rather than from `lambda_ratio`.

- **For:** the statistic is the quantity that matters — the width the user actually gets — and it is
  directly observable from the posterior with no intermediate construction.
- **Against:** the bar becomes small and therefore hard to distinguish from noise. Under this phase's
  measurements the honest pre-registered expectation would have been **≈1.1**, not 2.5.

### Which I would choose: **(b)**, with a caveat that changes its form

Construction (b) gates the quantity that is actually reported to a user, and it needs no
intermediate modelling step that could itself be wrong. Construction (a) is cleaner in isolation but
tests a derived object; a gate that can fail because the propagation was mis-built, rather than
because the model misbehaved, has the same category of defect this note exists to correct.

The caveat is in §5 and it is substantial: at the measured magnitudes, (b)'s correct bar sits so
close to 1 that a "must exceed" threshold cannot do the job.

## 3. The transferable rule

> **A bar and the statistic it is applied to must be the same quantity.**
>
> **And when a bar is derived from a component ratio, check that the component actually dominates
> the total before the bar can be reached at all.**

The second clause is the part that would have caught this one. A ratio of 5 on a component that
contributes ~20 % of a total's variance cannot move that total by 5, or by 2.5, or by much at all.
That check is one line of arithmetic and it is worth making mandatory whenever a threshold is
derived from a measurement of something narrower than what it will gate.

A useful generalisation: **derive the bar and the statistic in the same units, on the same object,
in the same expression.** If the derivation and the application appear in different files — as they
did here, `p11_consts.jl` versus `test_lambda_ablation.jl` — the mismatch is invisible at both ends.

## 4. Likely consumer: Phase 15

**Phase 15 (Calibration Operating Envelope and CI Gate)** is the likely consumer of this machinery.

**Any reuse should start from this note, not from `p11_consts.jl`.** The constant in that file is
correct as a historical record and wrong as a specification, and nothing in the file says so — by
design, since it is byte-frozen. This note is the pointer that prevents the value being lifted
verbatim into a new gate.

## 5. Design point for whoever writes Phase 15 — flagged, not resolved

The coverage result (`p11_coverage_report.jld2`) measures the **true value of the gated quantity**:
the RMSE of the Δρ posterior mean against simulated truth is flat across the ladder, ratio
**1.0003**. So the correct widening is ~1.00, and under construction (b) the correct bar would sit
**near 1**.

That breaks the gate's *form*, not just its number:

> A "must exceed 1.0" threshold cannot distinguish **"correctly flat"** from **"broken"**. Both
> produce a ratio near 1. A gate at that level has to be an **equivalence-style test** — two
> one-sided tests against a tolerance band around the predicted value — rather than a one-sided
> "must exceed" threshold.

This phase already owns the machinery for that: `spike/validation/p11_stats.jl` carries a TOST
implementation (`test_p11_tost.jl`), and plan 11-08 was scoped to use per-rung TOST equivalence with
an inverted Holm direction. That is the right shape for a re-derived SC1g.

**Deliberately not resolved here.** Choosing the tolerance band is a pre-registration decision, and
setting it after seeing this phase's measurements would be precisely the "amended twice, credibility
spent" pattern that Phase 7 was burned by. It belongs to whoever writes Phase 15, before they look
at their own data.

## Cross-references

- `11-DIAGNOSIS.md` — the measurements: σ_reg magnitudes, the coverage/RMSE result, the bottleneck
  and redesign evidence.
- `11-CLOSURE.md` — the phase decision, limits, and manuscript honesty items.
- `spike/validation/p11_probe_report.jld2` — `lambda_ratio`, `drho_eq` curves, the source of the
  original derivation.
- `spike/validation/p11_hybrid_report.jld2` — σ_reg per rung and the quadrature totals used above.
- `spike/validation/p11_coverage_report.jld2` — the RMSE-flat result that fixes the true value at ~1.
