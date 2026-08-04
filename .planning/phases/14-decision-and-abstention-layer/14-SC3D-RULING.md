# Phase 14 — SC3-d Disposition Ruling

**Ruled:** 2026-08-04 by the user (Manuel Seefelder)
**Scope:** the third of `14-14-PLAN.md` Task 2's honesty items — *how the failed SC3-d row is
recorded*. The `resume-signal` asks for exactly this: *"If any SC row FAILED, say how you want it
recorded."*

**This ruling covers SC3-d ONLY.** The other two checkpoint items — the Assumption A2 derivation's
referee-readiness (§4) and the framing of the five judgement-call bars (§6) — are **NOT ruled here**
and remain open. `14-14` must still obtain and record those two answers verbatim before the phase
closes.

---

## RULING

**Phase 14 closes NEGATIVE-BUT-USEFUL. SC3-d is recorded as a named limit SCOPED TO THE DENSITY-ONLY
WIRING — not as a limit of the shipped OOD flag.**

No rebuild. No fresh bar. No re-measurement. `P14_ITERATION_ALLOWANCE` is **not** spent — its single
pre-declared trigger (conformal coverage below the SC1-d band) did not fire, and SC1-d passed at
0.9115 against a band of 0.8799.

**Rejected, with reasons on the record:**

- **Routing SC3-d to Phase 15 as a pre-registered rebuild with a fresh bar.** Phase 15's OOD arm runs
  `gate_ood_roc`, i.e. the **fused density∨noise** detector, so a fresh bar there would measure a
  differently-equipped detector and answer a different question. It would also collide with Phase 15's
  **D-14**, which lists "add a detector channel for the failing mechanism" under Deferred and forbids
  it as an in-phase response.
- **Re-wiring `(:pp, :noise)` into the Phase-14 decision layer and re-measuring SC3-d.** It would very
  likely clear the bar (see the evidence below) — and that is precisely why it is refused. Improving
  the instrument after seeing the result, in order to pass a frozen bar, is the pattern this project
  has already had to correct four times (`.planning/STATE.md`, the four-amendments observation).
- **Recording SC3-d as a general limit of the OOD flag.** Factually too pessimistic, and falsifiable
  by opening the shipped gate report — see below.

---

## THE EVIDENCE THAT FORCED THE SCOPING

**What SC3-d measured.** Abstain rates at the fitted operating point: `texture` 0.999, `background`
0.989, `optics` 0.387, **`misspec_noise` 0.000 (0 of 1000)**, against an in-distribution rate of
0.129. Margin **−0.129** against the frozen floor `0.50`.

**What the Phase-14 decision layer actually wired.** `spike/p14/decide.jl:583-584` and `:979-980`,
and `spike/p14/pools.jl:487,506,563-564`, all record:

    ood_channels_wired     = (:density,)
    ood_channels_not_wired = (:pp, :noise)

and `decide.jl:286,293` switches the posterior-predictive channel off (`with_pp = false`). **One of
three channels.** This was declared in `pools.jl` before any arm was built — the ruling converts it
from a stated limit to a measured one, which is the useful part of the result.

**What the shipped detector does on the same families.** Read directly from
`artifacts/amended_v2/grid_8/gate_report_8.jld2`, key `report[:ood]`, four misspecification levels:

| family | `density_auc` (what P14 wired) | `noise_auc` | `fused_auc` (**shipped**) |
|---|---|---|---|
| `noise` | **0.0 / 0.0 / 0.0 / 0.0** | 1.0 / 1.0 / 1.0 / 1.0 | **1.0 / 1.0 / 1.0 / 1.0** |
| `optics` | 0.58 / 0.676 / 0.657 / 0.67 | 1.0 / 1.0 / 1.0 / 1.0 | 1.0 / 1.0 / 1.0 / 1.0 |
| `texture` | 1.0 / 1.0 / 1.0 / 0.995 | 1.0 / 1.0 / 1.0 / 1.0 | 1.0 / 1.0 / 1.0 / 1.0 |
| `background` | 0.559 / 0.959 / 0.0 / 0.96 | 1.0 / 1.0 / 1.0 / 1.0 | 1.0 / 1.0 / **0.0** / 1.0 |

The density channel is **exactly blind** to the `noise` family — AUC 0.0 at every level, i.e. worse
than chance. That is the *by-design* blind spot the noise channel exists to close:
`test/gate/misspec.jl:32-33` states the reported detector is the OR-fused density∨noise pair because
*"the image-noise channel closes the detector-noise blind spot of a correlation-only summary."*

**The scope is DENSITY-ONLY, not NOISE-FAMILY-ONLY.** `optics` shows the same effect more mildly —
density 0.58-0.68 against a fused 1.0 — and its 0.387 abstain rate is consistent with that. The
limit is a property of the wiring, not of one family.

**What this ruling does NOT claim.** The fused detector's SC3-d margin is **not measured** — Phase 14
never ran the fused detector through `decide_coloc`. It is therefore *unknown* whether SC3-d would
pass with all channels wired, and `14-REPORT.md` must not assert that it would. What is known is that
the channel SC3-d needed was present in the shipped artifact and absent from the decision layer.

---

## HOW `14-REPORT.md` AND `14-14-SUMMARY.md` MUST STATE IT

1. SC3-d is reported **FAILED**, with the measured margin −0.129 and the 0/1000 `misspec_noise`
   count. The failure is not softened, and the abstain-rate table stays.
2. Immediately beside it, the scope: the measurement was taken on a **density-only** detector
   (`ood_channels_wired = (:density,)`), and the shipped fused detector attains AUC 1.0 on the same
   family in `gate_report_8.jld2`. Cite the file and the key.
3. The honest one-line form: **"the Phase-14 decision layer's abstention is blind to detector-noise
   misspecification because it wires one of the shipped flag's three channels — not because the
   shipped flag is blind to it."**
4. State explicitly that the fused detector's SC3-d margin is unmeasured, and that **Phase 15's
   `noise` axis tests the fused detector independently** (Phase 15 D-02/D-06/D-07). Phase 15's
   three-valued verdict for that axis is the number that settles it.
5. Record that `P14_ITERATION_ALLOWANCE` is **unspent** and why (its trigger did not fire).
6. Record the two rejected alternatives above with their reasons, so the refusal to re-wire is
   visible as a choice rather than an omission.

---

## FOR PHASE 15 — a concrete prediction, recorded BEFORE that phase runs

Phase 15 sweeps a `noise` axis on the fused detector and assigns it PROTECTIVE / LATE /
SILENT-BUT-SAFE (D-06). On the AUC evidence above, **PROTECTIVE is expected**. Recording the
expectation now means a `LATE` outcome there would be a genuine surprise with evidence against it,
rather than a result read after the fact.

**A second observation to carry into Phase 15, flagged rather than explained.** `fused_auc[:background]`
drops to **0.0 at level 3** while `noise_auc[:background]` is 1.0 at that level. The z-score OR-fusion
is not monotone in its channels, so this is possible — but a clean zero exactly there is worth a look
before the `background` and new pure-offset autofluorescence axes (D-02a) are read. Not diagnosed
here; noted so it is not discovered in the reported run.

---

*Ruling recorded before `14-14` was resumed. The remaining two honesty items are open.*
