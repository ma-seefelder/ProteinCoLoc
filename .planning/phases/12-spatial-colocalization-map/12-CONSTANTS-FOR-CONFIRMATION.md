# Phase 12 — Constants requiring USER CONFIRMATION before Wave 1 executes

**Written:** 2026-07-28 · **Status:** OPEN — none of these has been confirmed.

> Collected in one place deliberately. Hunting these out of twenty plans would guarantee one gets missed,
> and in this repository a threshold quietly picked to make a plan complete is the same class of defect as
> the Phase-11 bar that was derived from a component ratio and applied to a total. **No value below was
> chosen to make a number pass; each is argued from design.** Where 12-01 already pre-registers a constant
> that covers a need, that is referenced rather than duplicated.

## 0. What is NOT here, and why

The great majority of Phase 12's Tier-1 constants were already pre-registered by `12-01-PLAN.md` before
the seven late plans were written, and the late plans consume them unchanged. Nothing in that set needs
confirmation: seeds and salts, the reserved counters, `P12_G`/`P12_K_DEV`/`P12_D_MINISPIKE`, the r₁
bracket and ladder, `P12_ABLATION_R1`, `P12_SIM02_*`, all the D-12 Stage-1 constants, all the Stage-2
constants (`P12_SBC_M/L/BINS`, `P12_SBC_TARGET_KS_ALPHA`, `P12_SBC_TOST_*`, `P12_COVERAGE_NOMINAL`,
`P12_STAGE2_COVERAGE_TOST_DELTA`, `P12_STAGE2_N_MIN`, `P12_STAGE2_LOGSCORE_MIN`), `P12_MASK_K_SET`, the
reporting-only bars, the budget constants, and `P12_ITERATION_ALLOWANCE`.

Also not here: the four **Tier-2** constants (`P12_CHOSEN_PRIOR`, `P12_K_PROD`, `P12_D_PROD`,
`P12_FISHERZ_NEFF`). Those are *measured* and appended with provenance, by design — confirming a value for
them in advance would be the data-snooping the two-tier structure exists to prevent.

## 1. `P12_R1_PRIOR` — the lattice prior over induced lag-1 correlation

| | |
|---|---|
| **What it gates** | Nothing directly. It **is** π(θ) for the spatial field: 12-07 draws from it, and 12-15's entire CAR-vs-GP comparison rests on both arms being drawn from the *same* prior. |
| **Proposed** | `Uniform(P12_R1_MIN, P12_R1_MAX)` = `Uniform(0.05, 0.95)` |
| **Derivation** | Both bracket endpoints are already pre-registered in 12-01; this only binds the distribution *object* so two plans cannot each construct their own. R-8 is the reason it is uniform on **r₁** rather than on α or ℓ: measured lag-1 correlation runs 0.136 → 0.947 across α ∈ [0.5, 0.999], so a uniform prior on α puts almost all mass on "no spatial structure" and would make D-08's correlation length unidentifiable *by parametrization rather than by physics* — the S-1 prior-echo trap. |
| **Too wide** | r₁ → 1 approaches an intrinsic (improper) CAR limit where the field is nearly constant; the deviation field carries almost no energy and the phase's premise becomes untestable at that rung. |
| **Too narrow** | Excluding the low-r₁ end removes the regime where per-region borrowing is *hardest*, which is precisely where D-12 Stage 1 needs discriminating power. Excluding the high end removes the regime where borrowing should be easiest, i.e. the positive control on the premise. |
| **Confirm** | That `Uniform` (not a Beta or a log-uniform) is the intended shape. A uniform on r₁ is the honest default given no prior belief about biological smoothness — which is also D-08's stated reason for inferring rather than conditioning. |

## 2. `P12_GOLDEN_IMSIZE` and `P12_GOLDEN_KEYS` — the D-06 golden-fixture scope

| | |
|---|---|
| **What they gate** | 12-08's exact-`==` byte-for-byte regression claim. They determine **which bytes are frozen**; a fixture captured at a different size or over different keys is a different fixture, and the regression claim is only meaningful if both were locked before the capture ran. |
| **Proposed** | Mirror Phase 11's `P11_GOLDEN_IMSIZE` / `P11_GOLDEN_KEYS` values exactly, **read from `p11_consts.jl` rather than retyped**, with keys drawn off `P12_FIXTURE_COUNTER = 99`. |
| **Derivation** | Not argued from any Phase-12 measurement — deliberately inherited, so the Phase-12 golden fixture is directly comparable to Phase 11's and any byte difference is attributable to the stage-1 field edit rather than to a changed capture configuration. |
| **Too large** | A large `imsize` makes the fixture slow and bulky in git for no added regression power; the stage-1 mixing path is size-independent. |
| **Too small** | Below roughly 128² the patch grid starts producing masked regions from the ≥15-surviving-pixel floor, so the fixture would freeze a degenerate case and the regression would not exercise the normal path. |
| **Confirm** | Whether inheriting Phase 11's values is right, or whether the fixture should be captured at the F5 mixture's modal size instead. Inheriting is my recommendation — comparability beats representativeness for a byte-equality tripwire. |

## 3. `n_low` — the identified/vacuous boundary in 12-18's SBC row classes

| | |
|---|---|
| **What it gates** | Which deviation coefficients get a **strict point-null KS** test versus a **nuisance-appropriate equivalence** test. Pitfall 8 is explicit that applying a point null to a non-identified row produces a loud, uninformative failure — this project has already paid for that twice. |
| **Proposed** | Declared as a named constant with its derivation recorded, and reported in the artifact as `n_low` alongside the full `row_classes` partition. **I am not proposing a number.** |
| **Derivation** | Should follow the DCT smoothness ordering: the first `n_low` modes are those a reader should *expect* to be identified from a 64-region summary. The principled route is to read it off 12-15's measured truncation curve (the rank past which truncation error falls below posterior noise) rather than to guess — but that curve does not exist until Wave 8, and this is a Tier-1 decision. |
| **Too high** | Rows that are genuinely vacuous get a point null and fail loudly for a reason that is not a calibration defect — exactly the Phase-7 outcome. |
| **Too low** | Rows that *are* identified get only an equivalence test, so a real calibration miss in a mode the map depends on could pass unnoticed. This is the more dangerous direction. |
| **Confirm** | Either (i) a Tier-1 value now, with the reasoning, or (ii) an explicit ruling that `n_low` is a **Tier-2 append** derived from 12-15's truncation curve — in which case 12-18 must record `n_low_source = :tier2_from_minispike` and 12-01's Tier-2 sentinel list needs a third entry. **(ii) is my recommendation**, because it makes the boundary a measurement rather than a guess, and because the plans already have the machinery for a provenanced append. It does mean 12-18 cannot run before 12-15, which the wave order already enforces. |

## 4. Guard-4 pass fraction — 12-20's radial-orthogonalization reading

| | |
|---|---|
| **What it gates** | Nothing mechanically (all four guards are reporting-only, and `P12_RADIAL_ENERGY_CEILING` / `P12_OFFSET_GRID_TOL` are members of `P12_REPORTING_ONLY_CONSTANTS`). But it is the **fraction of `headline_logscore_delta` that `radial_orth_logscore_delta` must retain** for the guard to read as "the advantage survives the confound being removed", and 12-20 requires it to be stated in the header *before* the run. |
| **Proposed** | No value proposed. |
| **Derivation** | Cannot be argued from any existing measurement: no repo artifact reports how much of a per-region log-score advantage is radial, which is the whole reason the guard exists. Any number I picked would be invention dressed as design. |
| **Too high** | The guard reads FAIL whenever any radial component contributes, including a genuinely radial *biological* field — and note the projection removes real radial signal too, so a strict fraction penalises truthful models. |
| **Too low** | A model whose advantage is mostly the chromatic artifact passes, which defeats the guard's entire purpose and would let a confounded win reach the manuscript. |
| **Confirm** | A fraction, with the reasoning. If you prefer, the honest alternative is to declare it **purely descriptive for v2.0** — report the retained fraction with no pass/fail reading at all, and let the Stage-2 verdict's written attribution argument (12-20's burden-shifting rule) carry the judgement. That avoids inventing a bar and is consistent with the guards being reporting-only. **This alternative is my recommendation.** |

## 5. `P12_GUARD_METRICS` placement — Tier 1 or runner-local?

| | |
|---|---|
| **What it gates** | Nothing. It is the frozen metric-symbol tuple that stops the reported guard quantities drifting between runs, mirroring `P11_PROBE_METRIC` (`p11_consts.jl:260`). |
| **Current** | Declared at the top of `run_p12_guards.jl` (12-20), not in Tier 1. |
| **Trade-off** | Phase 11's precedent puts the frozen metric symbol in the consts file, which makes drift a pre-registration breach. Runner-local is defensible for a reporting-only guard and keeps Tier 1 smaller. |
| **Confirm** | Low stakes; my recommendation is to leave it runner-local and say so, because promoting it to Tier 1 implies a gating status the guards deliberately do not have. Flagged only so the deviation from the Phase-11 precedent is a choice rather than an oversight. |

---

## Sequencing

Items **1, 2, 3 and 5** are Tier-1 and must be settled **before 12-01 executes**, because Tier 1 is
append-only afterwards and adding a constant later is itself a pre-registration breach. Item **4** is
needed before 12-20 runs (Wave 12), so it has slack — but it is cheapest to answer now, alongside the
others.

Separately and more urgently: the **Δρ semantics blocker** in `12-02-PLAN.md` is not a constant and is not
listed here. It blocks Wave 1 outright and is recorded in that plan and in `.planning/STATE.md`.
