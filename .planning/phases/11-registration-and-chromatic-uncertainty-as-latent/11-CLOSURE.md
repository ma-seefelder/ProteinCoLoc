---
status: closed
kind: closure
phase: 11-registration-and-chromatic-uncertainty-as-latent
closed: 2026-07-27
outcome: goal answered negatively; SC1g gate MIS-SPECIFIED; plans 11-08..11-11 SUPERSEDED
evidence: 11-DIAGNOSIS.md
gates: none
---

# Phase 11 — Closure

## What the phase set out to do, and what happened

**Goal.** Promote sub-pixel registration (and an optional chromatic warp) from a fixed simulator
nuisance to an inferred latent, so the coloc posterior widens *honestly* under registration
uncertainty instead of reporting false confidence.

**Outcome.** The goal is **answered, negatively, with evidence** — and the answer is more useful
than the one the phase expected.

1. **Registration at ≤3 px is not inferable** from an 8×8 patch-correlation summary. Not at finer
   resolution (SNR flat: 0.403 / 0.338 / 0.351 at 8×8 / 16×16 / 32×32), not with a self-referential
   probe-shift feature, and not at any colocalization strength (ridge/prior 0.993–1.001 across five
   |ρ| bins of n ≥ 1322).
2. **It also does not need to be.** The width the data actually demands is **flat in λ**: RMSE of
   the Δρ posterior mean against simulated truth is 0.13334 at λ_min and 0.13338 at λ_max, ratio
   **1.0003**. Registration uncertainty at ≤3 px does not measurably degrade Δρ estimation.
3. So the posterior does **not** report false confidence. It reports approximately the right width
   (1.021) for a quantity whose correct width ratio is ~1.00.

## The SC1g gate: MIS-SPECIFIED, not failed

**This is the finding that must not be recorded as "the net failed".**

`P11_LAMBDA_ABLATION_FACTOR = 2.502` was derived as half of `lambda_ratio = 5.0039` — the ratio of
the **registration-induced component alone** — but the tripwire applies it to the **total** Δρ
posterior SD, which contains every other source of uncertainty. Measured σ_reg is subdominant by
~2× at every rung (≤0.046 against a net width of ~0.091), so even exact analytic accounting reaches
only 1.134, and the calibration evidence puts the truth at ~1.00.

> **The net's measured ratios were 1.109 / 0.987 / 0.936. The correct answer is ~1.00. The net was
> right, and the gate failed it anyway** — because the gate demanded 2.502 of a quantity whose true
> value is 1.

**No threshold was edited.** Re-deriving the bar is a pre-registration decision and is left open;
this document records the evidence, not a new constant. Any future phase reusing SC1g must
re-derive it to compare like with like — either gate the registration *component*, or set a
total-width bar from the measured σ_reg rather than from `lambda_ratio`.

## Disposition of plans 11-08 through 11-11: SUPERSEDED

They are **not** silently skipped. The argument is stated here so it can stand or fall in the open.

| Plan | Disposition | Reasoning |
|---|---|---|
| **11-08** — SC2 ladder | **Superseded** | Its entry gate is mis-specified, and its question — does the Δρ posterior widen with λ — is already answered by two independent routes: analytically (σ_reg propagation, 1.134 as an upper bound) and empirically (coverage against simulated truth, RMSE ratio 1.0003). A 7-rung × 500-draw ladder would measure the same quantity at higher precision. It cannot change the conclusion, because the coverage test establishes that the truth itself is flat. |
| **11-09** — SC3 breakdown + D-02 attenuation | **Superseded in substance** | The attenuation question *is* the flat-RMSE finding. The beyond-prior breakdown at λ = 3.0 / 1.0 measures degradation of a quantity shown not to degrade in-prior; its in-prior baseline is the flat curve already measured. |
| **11-10** — D-17 real-image transfer | **Not applicable** | Its premise is that inferred registration transfers to real images. Nothing is inferred to transfer. Note this plan already carried a correction (`ea4a7d3`) for measuring on the wrong channel pair; that correction stands in the record. |
| **11-11** — figures, report, docs | **Superseded by `11-DIAGNOSIS.md`** | The deliverable was a results report plus a D-16 interpretation note. The diagnosis is that report, with stronger evidence than the planned ladder would have produced. |

**The counter-argument, stated fairly:** running 11-08 would produce the pre-registered SC2 numbers
at the pre-registered scale, which is what a reader expecting the registered analysis would look
for. If the manuscript needs the registered ladder *as registered*, 11-08 should be run in reduced
form (fewer rungs, fewer draws) purely for completeness — not because it can change the answer. That
is a manuscript-presentation decision, not a scientific one, and it is the user's.

## Limits and unresolved items

1. **Double-counting bracket.** The hybrid quadrature adds σ_reg to a net posterior that may already
   contain part of it. The truth is bracketed by **[1.021, 1.134]** — net contains all of it, versus
   none of it. **Every conclusion holds across the whole bracket**, since both ends are far below
   2.502. The coverage test resolves it to ~1.00, but the bracket is what makes the conclusion
   robust rather than dependent on a modelling choice.
2. **x/y anisotropy in the paired feature: a false alarm.** 37 % at n = 80 (2.8 SE) did not
   replicate and **reversed sign** at n = 200 (log-ratio +0.317 → −0.159). Scatter, not a warp
   asymmetry. The y-axis SNR of 1.227 was the favourable half of a coin flip.
3. **SNR legs measured in an unfavourable ρ regime.** They reuse the D-06 probe's five frozen theta
   bases, of which |ρ| = 0.18, 0.03, 0.19, 0.99, 0.47 — three of five in the least identifiable bin.
   This does not affect the stratified 128-row result (flat over 50,000 samples), only the SNR
   numbers.
4. **Research-net under-coverage.** Δρ intervals cover 0.72–0.77 against a nominal 0.90 (z-sd ≈ 1.4)
   at *every* λ. It is **flat in λ**, so it cancels in every ratio used here. It is consistent with
   the known `ρ_true` prior-atom artifact (no atom randomisation in this test) rather than a new
   defect, and it is a property of the **research net with its four declared deviations, not the
   shipped bundle**.
5. **Ridge legs are linear probes.** The magnitude legs are estimator-independent and the trained
   non-linear net showed the same flat response, but the strict scope of the ridge nulls is linear
   recoverability.
6. **Sign-aware paired variant is not viable as a summary.** It gains a 19 % shift-RMSE reduction at
   high |ρ| but degrades ρ_true recovery from 0.178 to 0.755, because flipping by sign destroys ρ's
   sign. Recorded so nobody adopts it on the strength of the shift gradient alone.

## Flag for the user — Phase 12's premise needs revisiting (NOT acted on)

Phase 12 depends on Phase 11 and was scoped to train on **"registration-aware θ"**. If registration
cannot be inferred — and this phase establishes it cannot, at ≤3 px from patch correlations — then
that dependency rationale no longer holds as stated. Phase 12's premise should be revisited before
it is planned. **This is noted, not acted on.** No roadmap edit to Phase 12 was made.

## Manuscript honesty items

In the style of the Phase-7 GO memo. Every item below is measured, not inferred.

1. **Registration is not inferred; it must be calibrated externally.** The 8×8 patch-correlation
   summary carries no recoverable information about sub-pixel shifts ≤3 px, at any colocalization
   level. State this positively: real microscopy calibrates channel registration with **fiducial
   beads**, so the method's scope is consistent with standard practice.
2. **State that registration uncertainty is a minor contributor**, with the number: at 8×8, exact
   analytic propagation of the forward model's own measured sensitivity inflates the Δρ posterior by
   ≤13 %, and the calibration evidence puts the true effect at ~0 % (RMSE ratio 1.0003).
3. **Do not report SC1g as a failed gate.** Report it as mis-specified, with the derivation error
   named: a component ratio applied to a total. The net's response was approximately correct.
4. **Scope every number to the research net.** The `d_in = 129`, `D = 8` net carries four declared
   deviations, is not the shipped bundle, and is not drop-in comparable to the shipped 128-row read
   surface. The shipped `amended_v2/grid_8` bundle and the Phase-7 GO (Option A, 2026-07-24) are
   untouched by this phase.
5. **Report the under-coverage with its caveats** if any calibration claim is made from this net —
   flat in λ, consistent with the known prior-atom artifact, research net only.
6. **The finer-grid family is ruled out with evidence**, not opinion: signal and noise floor rise
   together (8.9× vs 10.3× from 8×8 to 32×32).
7. **Report the negative redesign results.** Both candidates were built and tested and both failed;
   publishing the negative saves the field the same two attempts.
8. **Instrument check belongs in the methods:** the 8×8 SNR was measured as 0.403 and 0.406 by two
   independent implementations on separate code paths.

## Verification at closure

- `spike/validation/p11_consts.jl` — **byte-unchanged**. `LAMBDA_MAX`, `LAMBDA_MIN`,
  `P11_LAMBDA_ABLATION_FACTOR`, `SC2_SPEARMAN_FLOOR` all untouched.
- `src/`, `spike/Project.toml`, `spike/Manifest.toml`, `spike/contract.jl`, `spike/simulator/` —
  **byte-unchanged** since `78b170e`. No package installed.
- `P11_ITERATION_ALLOWANCE = 1` — **unspent (0 of 1)**.
- No retraining, no re-seed, no reship. `amended_v2/grid_8` untouched; ship gate not reopened.
- No SC2 ladder run. No wave 7–9 artifact exists.
- `11-07-SUMMARY.md` remains `status: blocked` — it is the honest record of where execution stopped.
