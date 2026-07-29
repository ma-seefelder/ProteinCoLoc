# Phase 12 — Constants requiring USER CONFIRMATION before Wave 1 executes

**Written:** 2026-07-28 · **CONFIRMED AND CLOSED: 2026-07-29.**

> **STATUS: CLOSED.** All four open items were confirmed by the user on 2026-07-29, each as recommended.
> Item 2 was withdrawn before confirmation and stays withdrawn with its false citation on the record.
> **The Tier-1 blocker on Wave 1 is lifted by this document.** Nothing below is open; the per-item
> **Confirm** rows are retained verbatim as the question that was asked, with the answer recorded beneath.
>
> | # | Item | Outcome |
> |---|---|---|
> | 1 | `P12_R1_PRIOR` | **CONFIRMED** `Uniform(0.05, 0.95)` — uniform on r₁, not a Beta or log-uniform |
> | 2 | `P12_GOLDEN_*` | **WITHDRAWN** (raised in error; false citation recorded, not deleted) |
> | 3 | `n_low` | **CONFIRMED option (ii)** — Tier-2 append derived from 12-15's truncation curve |
> | 4 | Guard-4 pass fraction | **CONFIRMED descriptive** — reported with NO pass/fail reading for v2.0 |
> | 5 | `P12_GUARD_METRICS` | **CONFIRMED runner-local** — deliberately not Tier 1 |

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

Also not here: the **Tier-2** constants — `P12_CHOSEN_PRIOR`, `P12_K_PROD`, `P12_D_PROD`, `P12_N_LOW` and
`P12_FISHERZ_NEFF`, five of them. Those are *measured* and appended with provenance, by design — confirming
a value for them in advance would be the data-snooping the two-tier structure exists to prevent.
**`P12_N_LOW` joined this list as a RESULT of item 3 below**: the question §3 asked was precisely whether
`n_low` belongs in Tier 1 or Tier 2, and the 2026-07-29 answer moved it here. So §3 appearing in this
document is not a contradiction of this paragraph — it is where the paragraph's fifth entry came from.

## 1. `P12_R1_PRIOR` — the lattice prior over induced lag-1 correlation

| | |
|---|---|
| **What it gates** | Nothing directly. It **is** π(θ) for the spatial field: 12-07 draws from it, and 12-15's entire CAR-vs-GP comparison rests on both arms being drawn from the *same* prior. |
| **Proposed** | `Uniform(P12_R1_MIN, P12_R1_MAX)` = `Uniform(0.05, 0.95)` |
| **Derivation** | Both bracket endpoints are already pre-registered in 12-01; this only binds the distribution *object* so two plans cannot each construct their own. R-8 is the reason it is uniform on **r₁** rather than on α or ℓ: measured lag-1 correlation runs 0.136 → 0.947 across α ∈ [0.5, 0.999], so a uniform prior on α puts almost all mass on "no spatial structure" and would make D-08's correlation length unidentifiable *by parametrization rather than by physics* — the S-1 prior-echo trap. |
| **Too wide** | r₁ → 1 approaches an intrinsic (improper) CAR limit where the field is nearly constant; the deviation field carries almost no energy and the phase's premise becomes untestable at that rung. |
| **Too narrow** | Excluding the low-r₁ end removes the regime where per-region borrowing is *hardest*, which is precisely where D-12 Stage 1 needs discriminating power. Excluding the high end removes the regime where borrowing should be easiest, i.e. the positive control on the premise. |
| **Confirm** | That `Uniform` (not a Beta or a log-uniform) is the intended shape. A uniform on r₁ is the honest default given no prior belief about biological smoothness — which is also D-08's stated reason for inferring rather than conditioning. |
| **✅ CONFIRMED 2026-07-29** | `P12_R1_PRIOR = Uniform(0.05, 0.95)`, uniform on **r₁**. The R-8 argument is the reason the shape is *right* and not merely permitted, so it is recorded here rather than left in the Derivation row: measured induced lag-1 correlation runs 0.136 at α = 0.5 to 0.947 at α = 0.999, so a uniform prior on **α** would put roughly 90 % of its mass on "no spatial structure" and would make D-08's correlation length unidentifiable **by parametrization rather than by physics** — the S-1 prior-echo trap, which this milestone has already met once in Phase 11. Parametrizing both arms by induced r₁ is what makes the CAR-vs-GP comparison a question about the kernels rather than about their coordinates. No plan change follows: 12-01's `THE r₁ PRIOR ITSELF` block already pre-registers exactly this binding in Tier 1, and 12-07 reads it rather than re-declaring it. |

## 2. ~~`P12_GOLDEN_IMSIZE` / `P12_GOLDEN_KEYS`~~ — WITHDRAWN, no decision needed

This item was raised in error and is withdrawn rather than deleted, so the record shows why.

It claimed Tier 1 should carry the golden-fixture scope "mirroring `P11_GOLDEN_IMSIZE` /
`P11_GOLDEN_KEYS` in `p11_consts.jl`". **That citation was false.** `grep -c GOLDEN
spike/validation/p11_consts.jl` returns **0**; the analog defines both constants LOCALLY in
`spike/test/capture_p11_golden.jl:51-52` as `(256, 256)` and `1:4`. 12-08 had already followed that
precedent correctly, defining the Phase-12 equivalents locally in `capture_p12_golden.jl`.

Putting them in Tier 1 would also have broken the capture: `capture_p12_golden.jl` guard-includes
`p12_consts.jl` for `P12_FIXTURE_COUNTER`, so a Tier-1 value differing from the local one is an
invalid-redefinition-of-constant error — and since the capture must precede the `forward.jl` edit, that
would have stalled wave 4. Removed from 12-01; 12-08 owns them.

## 3. `n_low` — the identified/vacuous boundary in 12-18's SBC row classes

| | |
|---|---|
| **What it gates** | Which deviation coefficients get a **strict point-null KS** test versus a **nuisance-appropriate equivalence** test. Pitfall 8 is explicit that applying a point null to a non-identified row produces a loud, uninformative failure — this project has already paid for that twice. |
| **Proposed** | Declared as a named constant with its derivation recorded, and reported in the artifact as `n_low` alongside the full `row_classes` partition. **I am not proposing a number.** |
| **Derivation** | Should follow the DCT smoothness ordering: the first `n_low` modes are those a reader should *expect* to be identified from a 64-region summary. The principled route is to read it off 12-15's measured truncation curve (the rank past which truncation error falls below posterior noise) rather than to guess — but that curve does not exist until Wave 8, and this is a Tier-1 decision. |
| **Too high** | Rows that are genuinely vacuous get a point null and fail loudly for a reason that is not a calibration defect — exactly the Phase-7 outcome. |
| **Too low** | Rows that *are* identified get only an equivalence test, so a real calibration miss in a mode the map depends on could pass unnoticed. This is the more dangerous direction. |
| **Confirm** | Either (i) a Tier-1 value now, with the reasoning, or (ii) an explicit ruling that `n_low` is a **Tier-2 append** derived from 12-15's truncation curve — in which case 12-18 must record `n_low_source = :tier2_from_minispike` and 12-01's Tier-2 sentinel list needs a third entry. **(ii) is my recommendation**, because it makes the boundary a measurement rather than a guess, and because the plans already have the machinery for a provenanced append. It does mean 12-18 cannot run before 12-15, which the wave order already enforces. |
| **✅ CONFIRMED 2026-07-29 — OPTION (ii)** | `n_low` is a **Tier-2 append**, not a Tier-1 number: `P12_N_LOW`, derived from 12-15's measured truncation curve and appended by 12-15 under a **third** Tier-2 sentinel `:P12_N_LOW` (12-01). 12-18 records `n_low_source = :tier2_from_minispike` as a mandatory artifact key. **This also RESOLVES the blocker escalated in `4d4c28a`**, and resolves it structurally rather than by agreement: if `n_low` had been fixed in Tier 1 and the measured `P12_K_PROD` had come back smaller, 12-18's `dev_low`/`dev_high` split would have been ill-defined and its strict row partition unbuildable — with the only escapes being a post-hoc clamp (moving a pre-registered boundary after the measurement) or a stalled wave. Under (ii) both numbers come off the **same** truncation curve in the **same** commit, so `P12_N_LOW <= P12_K_PROD` holds **by construction**, and 12-15 asserts it at append time so the construction is checked rather than trusted. Everything else in row 3 stands: Δρ is **not** a θ row, `targets = [p12_row_c0]`, Δρ is carried as a separate appended rank column, and the Δρ atom rate — higher than the 6.79 % single-map figure, because a difference of two `ghat`-clamped fields inherits atoms from **both** sides — is **measured and reported, never assumed**. |
| **Interaction with the Δρ decision** | Now settled and no longer open — but **not** in the form an earlier draft of this row stated. That draft said the ruling "carries a Δρ row, so `targets` = `c₀` **plus Δρ**"; it is corrected here rather than deleted, because it is the row a confirming reader would otherwise act on. **Δρ is NOT a θ row.** θ is 72 rows (`c₀`, 63 deviations, 7 nuisances, `r₁`) and Δρ is a *derived* quantity from a paired read, so `targets = [p12_row_c0]` and 12-18's strict-partition assertion over `1:D` would fail on a phantom Δρ index. Δρ is carried instead as a **separate appended rank column** outside the θ-row table, in the shape `sbc.jl:169-171` already uses. Consequence for `n_low`: the target class is **not** widened, the `dev_low` / `dev_high` split is unchanged, and `n_low` is untouched either way. What does survive from the old row is the atom point: the Δρ column needs randomized ranks in ρ space, because a difference of two `ghat`-clamped fields inherits atoms from BOTH sides and its atom rate is therefore HIGHER than the 6.79 % measured for a single map. 12-18 requires that rate to be measured and reported rather than assumed. |

## 4. Guard-4 pass fraction — 12-20's radial-orthogonalization reading

| | |
|---|---|
| **What it gates** | Nothing mechanically (all four guards are reporting-only, and `P12_RADIAL_ENERGY_CEILING` / `P12_OFFSET_GRID_TOL` are members of `P12_REPORTING_ONLY_CONSTANTS`). It *was* the **fraction of `headline_logscore_delta` that `radial_orth_logscore_delta` must retain** for the guard to read as "the advantage survives the confound being removed", and 12-20 *used to* require it stated in the header before the run. **That requirement was removed by the 2026-07-29 confirmation below** — recorded in the past tense here because this row is the question as asked, and the answer changed the plan. |
| **Proposed** | No value proposed. |
| **Derivation** | Cannot be argued from any existing measurement: no repo artifact reports how much of a per-region log-score advantage is radial, which is the whole reason the guard exists. Any number I picked would be invention dressed as design. |
| **Too high** | The guard reads FAIL whenever any radial component contributes, including a genuinely radial *biological* field — and note the projection removes real radial signal too, so a strict fraction penalises truthful models. |
| **Too low** | A model whose advantage is mostly the chromatic artifact passes, which defeats the guard's entire purpose and would let a confounded win reach the manuscript. |
| **Confirm** | A fraction, with the reasoning. If you prefer, the honest alternative is to declare it **purely descriptive for v2.0** — report the retained fraction with no pass/fail reading at all, and let the Stage-2 verdict's written attribution argument (12-20's burden-shifting rule) carry the judgement. That avoids inventing a bar and is consistent with the guards being reporting-only. **This alternative is my recommendation.** |
| **✅ CONFIRMED 2026-07-29 — PURELY DESCRIPTIVE for v2.0** | Guard 4 reports the retained fraction of `headline_logscore_delta` with **no pass/fail reading**. There is therefore **no fraction to state in advance**, and 12-20's "state the fraction in the header before running" requirement is **removed** rather than left as a requirement with nothing to satisfy it. The judgement is carried by 12-20's burden-shifting rule instead: crossing a bar obliges a **written attribution argument** in the Stage-2 verdict, and absent that argument the verdict is ITERATE or DESCOPE. The reason for refusing to propose a number is preserved deliberately: **no repo artifact reports how much of a per-region log-score advantage is radial** — which is precisely why the guard exists — so any value would have been invention dressed as design. This project has twice paid for a bar that could not do its job (the Phase-7 over-powered χ², the twice-amended Phase-11 gate); declining to invent a third is the correct call, not a deferral. Note the asymmetry that makes a descriptive reading honest here: the orthogonal projection removes *real* radial signal along with the artifact, so a strict fraction penalises a truthful model with genuinely radial biology. |

## 5. `P12_GUARD_METRICS` placement — Tier 1 or runner-local?

| | |
|---|---|
| **What it gates** | Nothing. It is the frozen metric-symbol tuple that stops the reported guard quantities drifting between runs, mirroring `P11_PROBE_METRIC` (`p11_consts.jl:260`). |
| **Current** | Declared at the top of `run_p12_guards.jl` (12-20), not in Tier 1. |
| **Trade-off** | Phase 11's precedent puts the frozen metric symbol in the consts file, which makes drift a pre-registration breach. Runner-local is defensible for a reporting-only guard and keeps Tier 1 smaller. |
| **Confirm** | Low stakes; my recommendation is to leave it runner-local and say so, because promoting it to Tier 1 implies a gating status the guards deliberately do not have. Flagged only so the deviation from the Phase-11 precedent is a choice rather than an oversight. |
| **✅ CONFIRMED 2026-07-29 — RUNNER-LOCAL, as a CHOICE** | `P12_GUARD_METRICS` stays declared at the top of `run_p12_guards.jl` and is **deliberately not** promoted to Tier 1. The reason is semantic, not convenience: Tier 1 is this phase's pre-registration of things that **gate**, and `P12_GUARD_METRICS` gates nothing — all four guards are reporting-only, `P12_RADIAL_ENERGY_CEILING` and `P12_OFFSET_GRID_TOL` are members of `P12_REPORTING_ONLY_CONSTANTS`, and 12-01's testset 6 asserts that set disjoint from `P12_GATING_CONSTANTS`. Putting the tuple in Tier 1 would imply a gating status the guards do not have and blur the one distinction this phase went to trouble to make machine-checkable. **This is a deliberate deviation from the Phase-11 precedent** (`p11_consts.jl:260` puts `P11_PROBE_METRIC` in the consts file), recorded as such so a later auditor reads it as a decision rather than an oversight. The anti-drift protection is unchanged and is not weakened by the placement: 12-20's verify asserts the artifact's metric keys are **exactly** `Set(P12_GUARD_METRICS)` (T-12-83), so a run reporting a different quantity under the same name still fails. |

---

## Sequencing

Items **1, 2, 3 and 5** are Tier-1 and must be settled **before 12-01 executes**, because Tier 1 is
append-only afterwards and adding a constant later is itself a pre-registration breach. Item **4** is
needed before 12-20 runs (Wave 12), so it has slack — but it is cheapest to answer now, alongside the
others.

**The Δρ semantics question is RESOLVED (2026-07-28)** and no longer blocks Wave 1: ship all three maps,
Δρ primary. It is recorded in `12-02-PLAN.md` as the §5 the frozen amendment must carry, and propagated to
12-09, 12-12, 12-16, 12-18 and 12-19.

**One finding from the corpus audit belongs on the user's desk even though it is not a constant.** The
Phase-8 external corpus supplies **zero** images today: `corpus/data/` is empty; all 32 rows have
`bytes = 0`; the two `physical-primary` rows carry `sha256 = PENDING-FETCH` and `split = sealed_holdout`
(reserved for Phase 16); the 30 `simulated-secondary` rows carry an empty `sha256` and cannot support a
real-image claim at any n.

So the real-image evidence base is the six committed TIFFs = **two specimens**, and two consequences follow
that the user should see:

- **SC3 is unaffected and runs at n = 2.** D-09's leave-region-out construction is single-stack end to end
  (it masks one region of ONE image and scores that image's own observed entry), so predictive coverage
  needs no control and runs on both specimens.
- **Δρ is NOT scored on real data**, and this supersedes an earlier note here that called
  `positive/`-against-`negative/` a legitimate pairing at n = 1. It is not legitimate: the simulated Δρ is
  a difference of two **exchangeable** draws from one prior, whereas the two specimens are deliberately
  non-exchangeable (`truth = coloc` against `truth = segregated`, different specimen types), so their
  contrast is between-population where the calibrated quantity is within-population. 12-19 records
  `real_delta_rho_computed = false`; the phase carries it as a named limit.

Neither is a threshold to set.

**THE CORPUS FETCH — RULED 2026-07-29: not a Phase-12 decision, and it is recorded here so PHASE 16
inherits the framing rather than rediscovering it.**

- **Phase 12 needs nothing from the corpus.** Its only two `physical-primary` rows are
  `split = sealed_holdout`, reserved for the Phase-16 blind evaluation and untouchable here, so fetching
  today would not give this phase a single usable image. The other 30 rows are `simulated-secondary` and
  cannot support a real-image claim at any n. Phase 12's real-image arm is the six committed TIFFs =
  two specimens, and that is settled.
- **It IS a Phase-16 prerequisite.** The `PENDING-FETCH` sentinels on the two physical rows must be
  replaced with real digests **before** Phase 16 opens the sealed holdout — an unfetched sentinel cannot
  verify the bytes a blind evaluation is supposed to be blind to. The positive anchor is a **~6.3 GB**
  archive whose download was deliberately left as an explicit human decision rather than automated, and
  that decision is still open and still belongs to the user.
- **Fetching would not create a matched pair**, so it does not reopen the Δρ question at any point: the
  sealed anchors are single specimens too, and the manifest's `role` column is positive/negative/benchmark
  — experimental controls on different specimens, not sample/control pairs of one specimen.

So: no action in Phase 12, and Phase 16 must treat the fetch as a gate on its own opening.
