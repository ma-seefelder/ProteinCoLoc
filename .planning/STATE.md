---
gsd_state_version: 1.0
milestone: v2.0
milestone_name: milestone
status: in-progress
stopped_at: Phase 15 context gathered
last_updated: "2026-08-04T12:15:34.523Z"
last_activity: 2026-08-04
progress:
  total_phases: 16
  completed_phases: 11
  total_plans: 124
  completed_plans: 104
  percent: 69
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-06-26)

**Core value:** A trained NPE/NRE produces calibrated, amortized colocalization inference (posterior + Bayes factor) in a single forward pass, >100x faster than per-dataset ADVI, with a demonstrated SBC/coverage proof and an honest OOD flag.
**Current focus (2026-07-24):** Phase 07 productionization is **COMPLETE (11/11 plans)**, built on the
GO decision (Option A — GO with named limits, `07-GO-NO-GO-UPDATE.md`). The v2.0 amortized-only
public API ships on the **8×8 reference grid only**. The prior calibration investigation (spikes
006-014) stands as the honest basis: coloc TARGETS ρ_true/Δρ are SBC-calibrated (randomized-rank atom
handling); BF is simulation-validated (spike 014, AUC 0.994, KDE baseline replaced); OOD passes;
nuisance marginal drift (~0.08 SD), the twice-amended gate, and atom handling are documented as named
limits (`docs/amortized.md`). No further training/gate iteration; Option B (BF §6 gate integration +
summary redesign) deferred. Findings in `Skill("spike-findings-proteincoloc")`.

## Current Position

Three phases are active concurrently. All lines are authoritative — do not overwrite one with another.

Phase: 14 (decision-and-abstention-layer) — EXECUTING
Plan: 11 of 14 complete (14-01…14-10; wave 6 in progress — SC1-b MET at all four
pre-registered levels, SC1-c reported and prior-sensitive)
executing** (started 2026-07-29 from plan HEAD `79c66d0`; 12 of 20 plans complete, waves 1-6 run)~~
**SUPERSEDED 2026-08-03: EXECUTION IS OVER — 16 of 20 plans complete, and the remaining four are
FORECLOSED by rulings on the record, not pending.** `12-VERIFICATION.md` exists (`gaps_found`); the
phase **did not achieve its stated goal**, and what remains is a bookkeeping scope-closure ruling, not
work. Struck through rather than deleted. See the dated section below.

Phase: 11 (registration-and-chromatic-uncertainty-as-latent) — **CLOSED 2026-07-27,
negative-but-useful. NOT active; nothing pending.** Restored 2026-08-03 (this line had been lost to
the documented single-slot clobber, `deferred-items.md`). Goal answered negatively with evidence:
registration ≤3 px is not inferable from the 8×8 patch summary at any coloc level, and does not need
to be (Δρ RMSE flat in λ, ratio 1.0003). SC1g **mis-specified, not failed**. Plans 11-08…11-11
superseded and now marked as such in their own frontmatter — **do not pick them up as live work; do
not retrain a registration-aware net.** No retraining/re-seed/reship occurred; `p11_consts.jl` and
`amended_v2/grid_8` byte-unchanged; iteration allowance unspent (0 of 1). Authority: `11-CLOSURE.md`.

## ⚠️ FOR PHASES 14/15/16 — a pool-overlap hazard that fires OUTSIDE Phase 12 (DEF-12-05, C-04)

**Rule:** `.planning/CONVENTIONS.md` **C-04**, added 2026-07-31 and verified at all three defining
sites — Phase 11, 12 and 13 **all** key their per-sample RNG on the global sample index, so this is
project-wide, not a Phase-12 quirk.

**Index-keyed generation has two faces.** The property that makes a lost pool regenerate
byte-identically (why Phase 11's 54 MB rebuild cost no iteration allowance) is the SAME property
that makes two pools at the same configuration share samples `1:min(n,m)` byte-identically.
**A "separate scoring pool" is a SUPERSET of the training pool, never a held-out set.**

**The live instance:** indices `1:10000` of 12-17's 50 000-sample production pool ARE 12-15's
10 000-sample mini-spike pool, byte-for-byte, **including the block 12-15 held out for scoring**.

**Harmless today** — nothing scores a 12-17 net on a 12-15 index range. **It BECOMES a leak if any
one of these happens:** (1) anything scores a 12-17-trained net on index range `<= 10000` of the
production pool; (2) anything reuses 12-15's held-out block as a held-out set for the *production*
net; (3) **a later phase generates "a fresh evaluation pool" at `arm = P12_CHOSEN_PRIOR`** — under
C-04 its first `min(n, 50000)` samples ARE the production training set. Any coverage, SBC or BF
number so obtained is measured on training data and is not a calibration result.

**Correct constructions:** draw fresh through `harness.jl` on the **validation** stream (structurally
safe — different salt, asserted disjoint at load time), or carve from the ONE pool via
`train_p12_npe`'s `pool_indices` and assert against the bundle's recorded `train_indices` **and**
`val_indices` — both, because the val block drives early stopping and is contaminated for model
selection *while passing every index check*. Full entry: Phase 12 `deferred-items.md` DEF-12-05.

---

## ✅ RESOLVED 2026-07-31 — the Stage-1 blocker below was ruled by the user (`ee78cfe`)

**`VERDICT: PROCEED`.** `.planning/phases/12-spatial-colocalization-map/12-STAGE1-VERDICT.md` now
exists, so `p12_stage1_verdict()` returns `:proceed` and `p12_require_proceed` admits 12-17, 12-18,
12-19 and 12-20. **Nothing is descoped** — SPAT-07 stays with 12-18 at full scope, SPAT-08 with 12-20.

**The ruling, stated without softening: the positive control did NOT pass as written.** `control_live`
is still `false` in the artifact and `stage1_pass` is still `false`; the ruling adjudicates those
values rather than changing them. The user ruled the frozen rule's *premise* false — a control above
the ceiling does not here mean a dead harness — because the own-row ridge attains the analytic
information limit `sqrt(1−corr²)` to within 5·10⁻⁴ at every rung, a bound computed with no ridge,
split or standardizer, so no defect in the audited machinery could manufacture agreement with it.

**One claim in the blocker below is CORRECTED by that commit, and the blocker text is left standing
as the historical record rather than edited.** The condensed retelling — including the "0.766 / 0.712
/ 0.637 at r₁ ≤ 0.50" sentence in the paragraph below — says the information limit exceeds the 0.5
ceiling at the **three** shortest rungs. It exceeds it at **four** (0.05 / 0.25 / 0.50 / 0.75), and at
the fourth the control **cleared the ceiling anyway** (0.49659), because the own-row bound does not
bound the full-128 control, which borrows from the other 63 regions.
`12-STAGE1-ADJUDICATION-BLOCKED.md` §4 stated this correctly with the qualifier *"from that region's
own row"*; condensation dropped it. Found by asserting the argument executably; both tuples are now
recorded in Tier-2 and re-derived by assertion, so the "three rungs" phrasing fails loudly if
reintroduced. **It sharpens rather than weakens the ruling**: the only thing that carries the control
under 0.5 anywhere is borrowing — the effect under test — so the control passes where that effect is
strongest and fails where it is weakest, and cannot certify the instrument independently of it.

**Tier 1 untouched, and the diff proves it** (171 added, 0 deleted). `P12_STAGE1_CONTROL_CEILING`
keeps its value of **0.5** as the historical record, as `SPEEDUP_GATE = 100.0` and
`P11_LAMBDA_ABLATION_FACTOR = 2.502` kept theirs. `P12_ITERATION_ALLOWANCE` stands at **1, UNSPENT**.
The adjudication is Tier-2 block 4 (sentinel `:P12_STAGE1_CONTROL_ADJUDICATION`) — every constant in
it is a measurement or a record, none is a bar, none enters `P12_GATING_CONSTANTS`.

**FOR PHASE 16, RECORDED SO IT IS NOT ASSEMBLED BY A REFEREE INSTEAD:** this is the **fourth** time
this milestone has concluded *"the bar was wrong, not the model"* — SC1g's component-vs-total, the
n=2-against-271 real arm, the Wald-labelled-Wilson sizing (DEF-12-04), and now a gating ceiling with
no derivation. The honest framing is that pre-registration is what SURFACES these and all four were
caught before corrupting a result — but the count is now high enough that "we pre-registered it" no
longer settles an argument alone. The common licensing standard, and the thing a referee should hold
each amendment to, is that every one was shown **from evidence independent of the machinery under
audit** to measure something other than what it named. Belongs in the manuscript's methods as a
stated observation. `12-STAGE1-VERDICT.md` §7 carries the table.

---

## ✅ RESOLVED 2026-07-31 (wave 8) — the NONE-BEATS-ABLATION blocker was ruled by the user: **DESCOPE ONTO THE D-13 ABLATION**

**USER RULING: accept the negative; the full spatial map is given up and the D-13 ablation becomes the
deliverable.** The blocker recorded at `### Blockers/Concerns` is resolved by this entry and is left
standing there as the record of the halt.

**THE RESULT HELD IN THREE RUNS.** `select_prior` returned `NONE-BEATS-ABLATION` in the original
unseeded 18-epoch run, in a seeded 18-epoch control, and in a seeded 100-epoch treatment. **The
qualifier "at mini-spike scale" stays on it.** Artifacts:
`spike/validation/p12_minispike_report.rerun_control_18ep.jld2` and `…rerun_treatment_100ep.jld2`;
verdict `12-15-RERUN-VERDICT.md`; reading pre-declared in `12-15-RERUN-PREDECLARATION.md` (`ea1d2b7`,
15:37:40+02:00, **committed before both artifacts**). Commit `cb52271`.

**THE UNDER-TRAINING HYPOTHESIS WAS TESTED AND ELIMINATED — this is what the ruling rests on.** The
18-epoch run ended on its epoch budget, not on convergence, so "under-trained" was the strongest
available objection to the negative. At a 100-epoch budget **early stopping fired at 56, 99 and 36
epochs** and validation risk fell **11–27 %**. Trained to convergence, **no arm became admissible and
no arm reached positive skill.** The best counter-explanation is gone, which is what turned a
provisional negative into one worth acting on.

**THE CAVEAT, AND IT MUST TRAVEL WITH THE NEGATIVE — THE DISQUALIFICATION WAS A CALIBRATION
DISQUALIFICATION, NOT AN ACCURACY ONE.** At 100 epochs **both spatial arms BEAT the ablation on
RMSE** — `:car` 0.99811 and `:gp` 1.10068 against `:none` 1.14663, i.e. by **12.95 %** and **4.01 %**.
They were disqualified **solely on coverage**. And coverage **crossed the admissible band without
landing in it**: all three arms **over**-cover at 18 epochs (0.99734 / 0.99998 / 0.97544) and all
three **under**-cover at 100 (0.79816 / 0.80028 / 0.82403), against the band [0.87, 0.93]. **On this
evidence admissibility tracks TRAINING DURATION, not the prior.** The honest statement is therefore:
***no spatial arm was ever both calibrated and better than the ablation — as gated, by a gate whose
behaviour on this evidence follows training duration.***

**NAMED LIMIT, and a candidate for the `12-STAGE1-VERDICT.md` §7 families.** The Stage-2 coverage
criterion is, on this evidence, **the same shape as the Stage-1 control ceiling: a criterion whose
pass/fail tracks something other than what it was meant to measure.** We are **not** pursuing it, and
not pursuing it is a choice this record shows rather than hides. The reason for refusing to hunt the
epoch count that lands coverage inside the band, verbatim: **that would be tuning a model until it
passes its own calibration gate — a gate that has then stopped measuring anything.**

**A CLAIM IN `12-MINISPIKE-VERDICT.md` IS WITHDRAWN, AND THE WITHDRAWAL IS PART OF THE RULING'S
BASIS.** That document's §1 called it "the stronger finding" that **no arm learned the field at all**
— every arm worse than a constant-zero predictor. **THAT DOES NOT SURVIVE A RE-SEED.** At the *same*
18 epochs with only the weight init changed, `:car` reaches skill **+0.08629**. "Every arm has
negative skill" is a property of **one initialisation**, not of the method. **The robust claim that
replaces it:** *no spatial arm is ever both calibrated and better than the ablation* — across the
three runs `admissible` was `[:gp]`, `[]`, `[]` and `beats_ablation` was `false` every time. Recorded
additively as `12-MINISPIKE-VERDICT.md` §7, with §1–§6 byte-unchanged.

**THE INIT-VARIANCE RESULT, ON ITS OWN LINE.** Same pool, same 18 epochs, only the init different:
`:car` RMSE **−36.68 %**, `:gp` **+51.05 %**, and **the two spatial arms swap which is worse**. That
is **two orders of magnitude larger than the 0.33 % margin** that decided the earlier rule-3 branch.
`12-MINISPIKE-VERDICT.md` §2.1 declined to read a ranking out of one run and argued from the
mechanism; **that refusal is now observation rather than argument.** LIMIT CARRIED AS DECLARED: this
is **n = 1 versus n = 1, it bounds nothing**, and nothing rests on it.

**CARRIED VERBATIM, UNCHANGED BY THE RULING:** ***a comparison between two failed arms is not
evidence about priors — the CAR-vs-GP question is UNRESOLVED, not answered.*** In the two new runs
**neither** arm was admissible, which is *further* from a prior comparison than the original, not
closer.

**NOTHING WAS APPENDED OR RELAXED.** `p12_consts.jl` byte-unchanged; `P12_CHOSEN_PRIOR`,
`P12_K_PROD`, `P12_D_PROD`, `P12_N_LOW` all absent, **all three sentinels intact**; `select_prior`
untouched; `src/`, `spike/Project.toml`, `spike/Manifest.toml`, `artifacts/`, `spike/data/cache/p11/`
clean. `K = 35` reported, **not** appended. Budgets, each against its own constant: `elapsed_min`
**11.29** and **29.29** against `P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150`; `datagen_min` **per arm**
against `P12_DATAGEN_WALLCLOCK_CEILING_MIN = 150` (0.09 / 0.00 / 0.00 in both runs) — their sum is
not a quantity anything compares to 150. `test_p12_train.jl` 107 pass, 0 fail.

### ⛔ NEW BLOCKER OPENED BY THIS RULING — THE STAGE-2 DESCOPE HAS NO EXECUTABLE ROUTE. USER DECISION REQUIRED.

**The plans do not cover the state the phase is actually in, and this is reported rather than
improvised.** Two descopes are defined; neither is ours:

| route | trigger | producer of the D-13 bundle | consumer |
|---|---|---|---|
| **Stage-1 descope** | `12-STAGE1-VERDICT.md` carries `VERDICT: DESCOPE` | **12-14 Task 3** fires | 12-16 → `:descope_ablation_only` |
| **Stage-2 descope** | 12-16's `lro_pass` returns `:fail`/`:disqualified` **after a full two-arm training** | — | spend `P12_ITERATION_ALLOWANCE` or descope (`12-VALIDATION.md:140`) |
| **ACTUAL STATE** | Stage-1 **PROCEED** + mini-spike selects **no prior** | **NONE — no plan fires** | 12-16 is specified to **THROW** |

- **`p12_stage1_verdict()` returns `:proceed`** (`12-STAGE1-VERDICT.md:3` = `VERDICT: PROCEED`), so
  `p12_require_proceed` still **admits** 12-17/12-18/12-19/12-20. Nothing is mechanically blocked —
  but **12-17 is impossible as written** (`12-17:137` keys its 50,000-pair pool on `P12_CHOSEN_PRIOR`,
  which does not exist).

- **12-14 Task 3 is the only executable producer of the D-13 deliverable, and it is keyed on
  `:descope`.** On `:proceed` it does nothing — it already ran as a no-op (`12-14-SUMMARY.md:31`), and
  re-running it changes nothing. **It has no branch for this state.**

- **12-16 resolves its bundle by FILE EXISTENCE, not by the verdict** (`12-16-PLAN.md:77-81`):
  `p12_train_full_report.jld2` → `:full_two_arm`; else the ablation bundle recorded in
  `p12_ablation_report.jld2` → `:descope_ablation_only`; **else THROW**, because "a third, unhandled
  case is not permitted to fall through silently". **Neither file exists.** So the descope *consumer*
  is intact and runnable; only the *producer trigger* is mis-keyed.

- **`12-CONTEXT.md:310-312` defines D-12 Stage 1 as the CAR-vs-GP mini-spike itself** — *"If neither
  does, descope before the full training run"* — which is verbatim what happened. **The implemented
  Stage-1 gate is a different measurement** (12-11's ridge borrowing probe) that already returned
  PROCEED. That divergence between the specified and the implemented gate is the root of the gap.

- **`12-15-PLAN.md:137-139`** requires this outcome be *"carried into the Stage-2 verdict, not silently
  overwritten"*. The Stage-2 verdict lives in 12-16, which cannot produce one without a bundle.

**WHAT WOULD BE IMPROVISING, AND WAS NOT DONE:** writing `VERDICT: DESCOPE` into
`12-STAGE1-VERDICT.md` to make 12-14 Task 3 fire — that would be **falsifying a gate record to obtain
a route**, and the Stage-1 adjudication is on the record; running Task 3's `:descope` body manually
under a `:proceed` verdict; or inventing a `stage2_gate` value 12-16 does not define.

**THE DECISION, FOR THE USER.** The D-13 recipe **is** fully defined executably — 12-14 Task 3's
`:descope` body specifies pool generation (`generate_p12_pool(P12_MINISPIKE_N; arm = :none)`, F5
mixture), training, durable save under `P12_PRIMARY_CHECKOUT` with `sha256`, and the
`p12_ablation_report.jld2` schema; 12-16 then scores it with `spat07_scope = :reduced_descope` and
`stage2_gate = :not_applicable_descope`. **What is missing is only the authorisation to spend that
compute (~11 min datagen + training) while the Stage-1 verdict reads `:proceed`.** Note
`train_p12_npe`'s `_assert_spend_allowed` **permits** `arm = :none` here, so it is not blocked — it is
simply not instructed.

**Two consequences the user must also rule on:** (a) **SPAT-08** (12-20's radial/offset/ε guards) has
**no descope coverage anywhere** and `12-11-PLAN.md:257` says to defer it **to v2.1** — but that
instruction is written for a *Stage-1* descope and does not literally apply here; **SPAT-07** would be
carried in reduced scope by 12-16. (b) **D-13 shipped from this route is the SPIKE-SCALE fallback**
(10,000 pairs), not the 50,000-pair version — `12-CONTEXT.md:323-324` warns the Stage-1-descope
ablation "exists only at spike scale", and the same reduction applies here.

**`P12_ITERATION_ALLOWANCE` — STATUS IS ITSELF AN OPEN QUESTION, flagged rather than ruled.** The
constant says the single allowance "is cheapest spent BEFORE the mini-spike — on the D-12 Stage-1
arm". The 100-epoch re-run **was** a documented iteration on the mini-spike arm with a pre-declared
reading, and the honest reading is that **it should be counted as SPENT**; it was not a threshold
relaxation. **This is the user's call, not an executor's**, and it does not block the ruling either
way.

---

## ⛔ 2026-08-03 — PHASE 12 HAS NO DISPATCHABLE WORK LEFT. A SCOPE-CLOSURE RULING IS REQUIRED.

**ADDITIVE ENTRY. Nothing above or below is edited.** An `execute-phase 12` run was dispatched on
2026-08-03 and **dispatched no executors**, deliberately. This records why, so the same run is not
attempted again in the belief that it was merely interrupted.

**`12-VERIFICATION.md` was produced** (status `gaps_found`, 8/9 must-haves; 6 verified, 2 accepted by
the recorded v2.1 deferral rulings, 1 failed). It is the first phase-level verification this phase has
had.

**THE FOUR REMAINING PLANS ARE FORECLOSED — BY RULINGS ALREADY ON THE RECORD, NOT BY ANY NEW
DECISION.** `init.execute-phase` reports 5 incomplete plans. None is dispatchable:

| plan | requirement | why not dispatchable |
|---|---|---|
| **12-11** | SPAT-02 | **Complete in fact.** All three artifacts on disk and `12-STAGE1-VERDICT.md` reads `VERDICT: PROCEED` (`ee78cfe`). Only its SUMMARY frontmatter still says `status: blocked` — **stale, and left byte-unchanged** rather than edited, per this phase's append-never-overwrite discipline. It is what makes the SDK count 5 rather than 4. |
| **12-17** | SPAT-03, SPAT-04 | **Impossible as written, verified on disk:** it keys its 50,000-pair pool on `P12_CHOSEN_PRIOR`, which does not exist and is *asserted absent* at `spike/validation/p12_consts.jl:798`. Dispatching it throws. Superseded by the 2026-07-31 descope. |
| **12-18** | SPAT-07 | **DEFERRED TO v2.1** by user ruling (`12-D13-AUTHORISATION.md` §9). |
| **12-19** | SPAT-06 real arm | **Hard-blocked:** `12-19-PLAN.md:202` resolves and hash-verifies **both bundles** from `p12_train_full_report.jld2`, which only 12-17 produces and which does not exist. |
| **12-20** | SPAT-08, SPAT-09 | SPAT-08 **DEFERRED TO v2.1** by user ruling (§7). |

**THE HEADLINE, AND IT IS A NEGATIVE: THE PHASE DID NOT ACHIEVE ITS STATED GOAL.** Two of the three
things the goal sentence promises are false on the record, and both were already known — verification
assembled them in one place rather than discovering them:

- **"spatial"** — no CAR or GP prior ever beat the neutralized ablation (`NONE-BEATS-ABLATION`, three
  runs). **The shipped deliverable is the ablation** — the same network with spatial borrowing off —
  carrying the `SpatialColocResult` type. No trained spatial arm exists at any scale.

- **"calibrated per-region uncertainty"** — SBC calibration (SPAT-07) is deferred and never ran, and
  the one in-scope calibration-adjacent measurement **failed its own pre-registered band**: pooled
  leave-region-out coverage **0.9707** against [0.87, 0.93] (`12-16-SUMMARY.md:29`, *"SPAT-06 IS NOT
  MET"*). Read this beside `12-D13-AUTHORISATION.md` §10.1 — **no arm met that band**, and the
  selected arm was never tested against it.

**WHAT WAS DELIVERED AND STANDS:** SPAT-01, SPAT-02, SPAT-05 and SPAT-09, plus `p12_coloc_map` — a
per-region Δρ map **with** per-region uncertainty, which is strictly more than `LocalColocMap` offers
today. **At spike scale (10,000 pairs), not the 50,000-pair version**, per §4.

**HARD CONSTRAINTS RE-VERIFIED INDEPENDENTLY AND ALL HOLD:** `src/` untouched since `ca02b0e`
(2026-07-25, *before* Phase 12 began); `spike/Project.toml` / `Manifest.toml` unchanged since Phase 9;
all three sentinels (`P12_CHOSEN_PRIOR`, `P12_N_LOW`, `P12_K_PROD`) still absent.

**THE RULING NEEDED IS BOOKKEEPING, NOT SCIENCE — and it is NOT taken by an executor.** Every
scientific question here is closed; §10.5 records "document only" and declined even the no-compute
option. What is open is only how the phase is **recorded**: SPAT-03, SPAT-04 and SPAT-06 are **not
delivered** and, unlike SPAT-07/08, **carry no v2.1 deferral ruling**. Phase 12 must not be marked
complete against 9 requirements while 5 are undelivered. **ROADMAP.md still shows all of
SPAT-01..SPAT-09 as `Pending` and Phase 12's plan list as `TBD`** — neither was touched by this run.

**NOTHING WAS RUN, RE-RUN, TUNED OR RELAXED.** No compute; no threshold moved; `p12_consts.jl`,
`src/`, every `*-VERDICT.md` and `12-D13-AUTHORISATION.md` byte-unchanged. The `for`-loop admissibility
bound at `run_p12_minispike.jl:249` is **left wrong and documented as wrong**, per §10.1.

---

## ⛔ BLOCKER (HISTORICAL — RESOLVED ABOVE 2026-07-31; retained as the record of the halt)

**12-11 measured the Stage-1 gate and deliberately wrote NO verdict.** `12-STAGE1-VERDICT.md` does
not exist, so `p12_stage1_verdict()` returns `:absent` and `p12_require_proceed` correctly refuses
12-17, 12-18, 12-19 and 12-20. Full evidence: `12-STAGE1-ADJUDICATION-BLOCKED.md`. Nothing was
tuned; `p12_consts.jl` is byte-unchanged and `P12_ITERATION_ALLOWANCE` stands at **1, UNSPENT**.

**The question, in one sentence:**
> Does a positive control that is provably **at the information limit of the data**, but above
> `P12_STAGE1_CONTROL_CEILING = 0.5` at r₁ ≤ 0.50, mean the harness is not live — or does it mean
> the ceiling is unreachable in that unit at those rungs?

**Why it cannot be answered by an agent.** The gate has two halves and they disagreed.
`borrowing_ok = TRUE` at **5 of 5** rungs (`ratio_mean` 0.926 / 0.847 / 0.719 / 0.562 / 0.336, all
under the 0.95 ceiling) — the half the phase is actually about passed by a margin. But
`control_live = FALSE` (max `control_ratio_mean` 0.732 vs a ceiling of 0.5). The frozen rule says a
failed control means the run is uninformative — record NEITHER verdict — and spend the single
iteration allowance "fixing the harness". That instruction rests on a premise the rule states rather
than tests, and 12-11 tested it: all three causes the rule names (broken standardizer,
target/predictor misalignment, wrong-`r1` pool) are ruled out **executably**, and the ridge attains
the analytic information limit `sqrt(1−corr²)` to within **5·10⁻⁴ at every rung** — a bound computed
with no ridge, split or standardizer, so no defect in the machinery could produce it. That limit is
itself **0.766 / 0.712 / 0.637** at r₁ ≤ 0.50, i.e. **above the 0.5 bar**, so at those rungs *no
estimator of any kind* can clear it. The Phase-11 benchmark also reproduces (global-level control
0.19687 vs Phase 11's 0.157). So all three contents the plan permits — PROCEED, DESCOPE, "harness
not live" — are untrue, which is why no verdict was written.

**The two defensible readings** (both are in the record; this is why it is not an agent's call):

- **A — the rule stands.** The ceiling is Tier 1 and in `P12_GATING_CONSTANTS`; reinterpreting a
  gate constant after seeing the number it failed is the "amend after the fact" pattern this project
  has already paid for twice. ⇒ run uninformative, allowance spent on the Stage-1 arm.

- **B — the ceiling is mis-scaled for the unit.** `P12_STAGE1_RATIO_CEILING` carries a derivation
  *and* an explicit recorded warning that the unit it is applied to is not the unit it was derived
  for. **`P12_STAGE1_CONTROL_CEILING` carries no derivation anywhere in the pre-registration** — its
  only stated rationale is 12-RESEARCH's qualitative "must be far below 1.0", which the measured
  0.324-0.732 (a 27-68 % error reduction) meets at every rung; and its plausible implicit anchor is
  Phase 11's *global* 0.157, which 12-CONTEXT S-2 explicitly warns must not be read as a per-region
  expectation. ⇒ control is live, gate reduces to `borrowing_ok`, TRUE at 5/5.

If the ruling is B, the bar would need expressing per rung or against the measured information
limit, and that is a **Tier-2 append** naming this artifact — never an edit to append-only Tier 1.

**The finding stands either way, and it is the deliverable Pitfall 2 named in advance.** Neighbours
DO inform a held-out region, monotonically in r₁. The crossing where the other 63 regions beat a
region's OWN measurement sits essentially exactly at **r₁ = 0.75** (parity, Δ = 0.0063; at r₁ = 0.95
neighbours win 0.336 vs 0.472). Carry into 12-20 regardless: `radial_r2 = 0.749` at r₁ = 0.95 is an
S-4 warning, with a competing benign reading (CAR corner/edge variance geometry, which 12-07
measured) that this run does not separate.

**Forbidden until the ruling exists:** running 12-17/18/19/20; editing any `p12_consts.jl` constant;
re-running the ladder on a different seed/size/arm to obtain a nicer control; spending the iteration
allowance; and reading `stage1_pass = false` as a DESCOPE (it is an AND over one TRUE and one FALSE
component, and the frozen rule forbids reading the FALSE half as a descope).

**Not blocked by the gate:** 12-13, 12-14, 12-15, 12-16. The orchestrator halted them anyway rather
than assume the unresolved question has no blast radius on work built from this pool — resuming them
early is available if the user wants progress while deciding.

---

Phase: 12 (spatial-colocalization-map) — execution record
Plans: 20 of 20 written and committed (`927cf7a`); FOUR plan-checker gates run 2026-07-28/29.
**EXECUTION MODE — NO WORKTREES, sequential on the MAIN working tree, deliberately.** Same ruling as
13-11 and for the same reason, plus one Phase-12-specific reason that is stronger: (a) 12-09's pool
lands in gitignored `spike/data/cache/p12/` under a 150-min append-only ceiling and is consumed by
12-11/13/14/17 in *later waves*, so a worktree cleanup between waves would destroy it exactly as it
destroyed Phase 11's 54 MB pool; (b) `spike/data/cache/p11` — the 54 MB READ-ONLY Phase-11 pool that
12-15 and 12-17 depend on — is gitignored and therefore **does not exist inside a fresh worktree at
all**, so those plans could not read it there. Executors are told not to write STATE.md/ROADMAP.md;
the orchestrator owns those writes, additively, because a phase-13 executor is live on this branch.
Wave order follows the plans' DECLARED `wave:` frontmatter (12 waves), not the SDK's DAG-derived
regrouping: the declared order deliberately staggers 12-15 (wave 8) and 12-16 (wave 10) so their two
Tier-2 appends to `p12_consts.jl` serialize; the DAG collapse would have put both in wave 8 and raced
them on the same append-only file.
**STATUS: Planning complete and UNBLOCKED. All user decisions received and propagated: the two of
2026-07-28, and all four remaining constants CONFIRMED 2026-07-29 (12-CONSTANTS-FOR-CONFIRMATION.md is
now CLOSED). No Tier-1 constant awaits confirmation. Wave 1 has no remaining blocker.**
Gate 1 found 7 blockers (5 defects, fixed; 2 escalated as the decisions below). Gate 2 found the
prose-vs-executable Δρ gap: the three-map ruling had been propagated as PROSE into six plans while the
executable specs stayed single-map. Gate 3 found 6 blockers and 11 concerns, all closed. Gate 4 (2026-07-29,
commits 50a83bf a6863c3 b7bfe2b ed92505 + this one) found 8 more, of which the largest was NOT fix-induced
and had survived all three earlier passes: **the production head was unbuildable** — `train_p12_npe` had no
width keyword, so a `P12_K_PROD < 63` bundle would have been written at 72 rows and then refused by its own
loader for every downstream consumer. Gate 4 also found the frozen amendment carrying a two-line
self-contradiction, six verifies that could not fail, and Δρ — the primary deliverable — getting no SBC
test, key or verdict. **Gate 4's own fixes were then audited and had introduced 3 blockers and 11 concerns
of their own, all fixed in the same session; that trend has still not broken.** The two decisions are now
answered:

  1. **Δρ semantics — RESOLVED: ship all three maps.** `region_delta_rho` stays the primary named
     deliverable per `src/results.jl:181`, computed as the per-region MC difference mirroring
     `infer.jl:108`; `region_rho_sample` and `region_rho_control` are additionally exposed so a reader
     can see WHERE a difference comes from. Costs zero extra forward passes — both single-stack reads
     already exist inside a paired Δρ. 12-12's `[-2,2]` range guard is REPLACED (not tightened) by a
     structural identity `delta ≈ sample − control`, because with a true Δρ the range genuinely IS
     [-2,2] and a range can never distinguish the two quantities. Propagated to 12-02 §5, 12-12,
     12-16, 12-18, 12-19. **The pool stays SINGLE-STACK** — an earlier version of this line said
     "12-09 (paired control draw)" and that was wrong: the Δρ\* truth is built at SCORING time from two
     independent prior draws (`harness.jl:124-136`), and a paired pool would have pushed datagen from a
     measured 65-80 min/50k to ~130-160 min against an append-only ceiling of 150.

  2. **Real-image SC3 arm — REPORTED, not gated, confirmed by a counted corpus audit.**
     `corpus/data/` is EMPTY: all 32 manifest rows have `bytes = 0`; the two `physical-primary` rows
     carry `sha256 = PENDING-FETCH` and are sealed for Phase 16, while the 30 CBS `simulated-secondary`
     rows carry an EMPTY `sha256`. (An earlier version of this line said all 32 were `PENDING-FETCH`;
     only the 2 physical rows are.) So the Phase-8 corpus supplies ZERO images today. The corpus `role`
     column is positive/negative/benchmark — experimental controls, i.e. DIFFERENT SPECIMENS, not
     sample/control pairs. The real data stays six FILES = **2 specimens**, and they are **not** a
     legitimate (sample, control) pair — an earlier version of this line called them one, and
     12-CONSTANTS-FOR-CONFIRMATION.md now says plainly that it is not legitimate. The SC3 criterion is
     SINGLE-STACK (D-09 masks one region of ONE image and scores that image's own observed entry), so
     predictive coverage IS computable on both specimens: n = 2, ~64 regions each.
     Δρ is NOT computed on real data, and the reason is **EXCHANGEABILITY, not shared nuisances.** An
     earlier version of this line claimed "12-09's simulated pair shares all seven nuisances by
     construction"; that is FALSE — `harness.jl:124-136` draws two independent priors, so nuisances
     differ between sample and control in the simulated Δρ too. The correct argument: the simulated pair
     is two *exchangeable* draws from one prior, so Δρ\* is a WITHIN-population difference, which is what
     the net is calibrated against; the two anchors are deliberately non-exchangeable (truth=coloc vs
     truth=segregated, different specimen types), so their contrast is BETWEEN-population. Scoring a
     between-population contrast against a within-population gate would be two arms measuring different
     quantities under one gate. That makes it a named limit, not an omission.
     A ±0.03 tolerance at n=2 is still noise, so the Stage-2 gate rests on 12-16's simulated arm
     (N ≥ 271). Recorded as clause (h) plus §5 of the frozen amendment.
     The effective independent n is **2**, not 128: the 128 region-draws are pseudo-replicates.
  **CORPUS FETCH — RULED 2026-07-29: not a Phase-12 decision.** Phase 12 needs nothing from the corpus:
  both `physical-primary` rows are `sealed_holdout`, reserved for Phase 16 and untouchable here, so
  fetching would not yield this phase a single usable image. It IS a **Phase-16 prerequisite** — the
  `PENDING-FETCH` sentinels must be replaced with real digests before Phase 16 opens the sealed holdout —
  and the positive anchor is a ~6.3 GB archive whose download was deliberately left as an explicit human
  decision, still open and still the user's. Fetching would not create a matched pair either: the sealed
  anchors are single specimens too. Framing recorded in 12-CONSTANTS-FOR-CONFIRMATION.md so Phase 16
  inherits it.
  **CLOSED 2026-07-29 — the `n_low` / `P12_K_PROD` blocker raised by gate 4 is RESOLVED structurally.**
  The user confirmed option (ii): `n_low` is a Tier-2 append (`P12_N_LOW`), appended by 12-15 under a third
  reserved sentinel `:P12_N_LOW` and read off the SAME truncation curve that yields `P12_K_PROD`. So
  `P12_N_LOW <= P12_K_PROD` holds BY CONSTRUCTION rather than by luck, and 12-15 asserts it at append time.
  12-18 reads the binding, declares no local `n_low`, records `n_low_source = :tier2_from_minispike`, and
  stops rather than defaulting if the constant is absent. No clamp was needed, so no pre-registered
  boundary moved after a measurement.
No code has been written, no seed consumed, no compute spent. `spike/` and `src/` untouched by Phase 12.

Phase: 13 (three-hypothesis-amortized-bayes-factor) — **ALL 16 PLANS COMPLETE**
Plan: 16 of 16 complete (13-01…13-16); the phase report 13-14 is written
**STATUS: Complete — waves 1-9 done (the gate arm, the α-series arm, the real-image arm and the
report). τ MEASURED (0.15), Phase-11 basis BOUND, the three-way evidence net is TRAINED (13-11), THE
AMENDED GATE HAS RUN AND PASSED (13-12, one run), the REPORTED α-ladder has run (13-13), the REPORTED
real-image arm has run (13-16), and `13-REPORT.md` is written (13-14). D-04 allowance still UNSPENT.**

**APPENDED 2026-07-29 by plan 13-17 (ADDITIVE — nothing above this block was altered).
PHASE 13 IS NOW 17 OF 17 PLANS. THE PRE-REGISTRATION WAS AMENDED AND THE REAL ARM RE-RUN.**

`13-D15-AMENDMENT.md` was applied to `spike/p13/consts.jl` in **ONE edit carrying TWO separately
justified changes**, per the user's ruling of 2026-07-29:

> **M1 — fix at the source, inside the same authorised amendment. M2 (the isolated-module read) is
> the tool for files that CANNOT be edited because no amendment authorises opening them.**

The reasoning, recorded because a mechanism chosen for a reason is auditable and one chosen by
default is not: *the byte-lock breaks either way, and `h.consts_sha == p13_consts_sha()` has to be
re-derived in the gate, alpha and realimage runners regardless — so the include-guard fix rides
along at ZERO additional pre-registration cost. That coupling is precisely why the answer is M1 and
not M2.*

- **CHANGE A — the channel pair.** The constant named the **wrong physical object**: `c1` is the
  DAPI/Hoechst nuclear counterstain (`test/runtests.jl:105`), so both pre-registered pairs measured
  counterstain-versus-protein overlap. `P13_REAL_CHANNEL_PAIR` `(1,2)` → **`(2,3)`**;
  `P13_REAL_REDUNDANCY_PAIR` **DROPPED** with no replacement; `P13_REAL_ANCHOR_MBAR` →
  **`0.4603 / 0.3815`, UNMASKED**. The masked `+0.8238 / −0.0338` are **forbidden** and appear
  nowhere in the arm, though they separate ten times better and no prior source forbids them.

- **CHANGE B — the include-guard sentinel (DEF-12-03).** `consts.jl:100` guarded the whole Tier-1
  body on `:P13_DEV_SEED`, a name `spike/validation/p12_consts.jl:109` **must legitimately mirror**
  to assert seed disjointness. Re-pointed to **`:P13_DECLARED_DEVIATIONS`** at the body wrapper and
  **all eight callers** — a caller-only fix would have been provably useless. **CHANGE B changed no
  value.** `spike/validation/p12_consts.jl` is **BYTE-UNCHANGED** and `P13_DEV_SEED` is
  byte-unchanged. **The 28 OTHER poisoned guards found by the same sweep are DOCUMENTED with a
  remedy each in `.planning/CONVENTIONS.md` C-01 and were deliberately NOT FIXED.**

**The `consts_sha` guard was re-derived, not weakened.** `spike/p13/net.jl` now defines the named,
dated `P13_CONSTS_SHA = (pre_amendment = "100e97a3…", post_amendment = "d6a6e63b…")`, asserted on
**BOTH sides at all three runner sites**. Strictly stronger than the single equality it replaces,
which only required the two sides to agree with each other. The widened `||` form is rejected in a
code comment. No guard deleted; `gate_report.jld2` / `alpha_report.jld2` byte-unchanged; no net
retrained; `P13_ITERATION_ALLOWANCE` **1 of 1, UNSPENT**.

**Suite signatures.** BEFORE: aborted at `runtests.jl:174` in the **Phase-12** block — *not* where
the amendment predicted — so Phase 13 was never reached; a concurrent Phase-12 fix (`8880c29`)
changed that between runs, so **the raw before/after comparison is confounded and is not claimed**.
The clean evidence for CHANGE B is a **targeted reproduction**: mirror loaded, then
`test_p13_consts.jl` — **22 pass / 1 fail / 114 `UndefVarError` before, 137 pass / 0 error after**.
AFTER: the Phase-13 block **executes** (`consts` 137/137, `labels` 164/164, `alpha` 129/129). The
§5B.5 `UInt32`/`UInt64` redefinition knock-on was exercised and is a **non-event**; **no seed was
edited** and the STOP RULE did not trigger.

**The re-run (one run, no RNG stream, no counter collision).** **The conclusion is UNCHANGED** —
qualitative, n = 2, unlabelled, OOD-bound. The verdict moved **RANDOM → COLOC** reading the positive
fixture as sample (`log BF(C:R)` −0.4658 → **+5.70311**); the OOD headline moved **417.2974 /
167.5446 / 2.491× → 703.2995 / 167.5446 / 4.198×**, i.e. the binding limit got **WORSE**. Reported
as a correction, never as a rescue.

**BLOCKERS CLOSED by this plan:** the two channel-pair blockers of 2026-07-25 (`ea4a7d3`), and the
`P13_DEV_SEED` include-guard blocker / **DEF-12-03** — closed **by reference**, without editing
Phase 12's `deferred-items.md`, because a live Phase-12 executor owns that file.

**⚠ NEW BLOCKER OPENED — A PRE-REGISTRATION RULING IS NEEDED, AND 13-17 DID NOT RESOLVE IT.**
On the corrected pair **two claims this phase made do not survive**, and the three
`spike/test/test_p13_real.jl` assertions encoding them were **left FAILING rather than rewritten**,
because rewriting a substrate expectation to match what was measured *after* the measurement is the
act the amendment exists not to be:

1. **The D-16 α-ladder no longer crosses zero** (`alpha_star_real` `nothing` / `nothing`, was
   0.875 / 0.625; `m-bar` stays positive on both fixtures). So *"negative induced μ is
   CONSTRUCTIBLE from real microscopy pixels via the D-16 mask-based reassignment"* — a sentence
   written for the manuscript — **is RETRACTED**; it was an artefact of segregating a target
   channel against a nuclear counterstain.

2. **The D-05 coherence check that previously PASSED now DISAGREES**: contrast **0.09819** is inside
   `tau = 0.15` so the label rule assigns RANDOM, while the net's argmax says COLOC.

3. **Consequence:** the abort moved to `runtests.jl:223`, masking **seven** Phase-13 files. All
   seven were run individually and pass (`test_p13_correction.jl` retains its 2 known deliberate
   misses).
**DECISION NEEDED:** accept the three misses as named limits, amend the D-16 real-substrate
expectation, or re-order the harness so they stop masking siblings. **No seed, bar, floor, band or
allowance may be moved to resolve this**, and `P13_ITERATION_ALLOWANCE` cannot be spent on a
real-image observation (`13-SC2-AMENDMENT.md` §7).

**✅ BLOCKER CLOSED — RULED 2026-07-29. APPENDED by plan 13-17 (ADDITIVE — nothing above this line
was altered). THE RED STAYS, AND THE MASKING IS FIXED.**

**The ruling: accept the three misses as NAMED LIMITS, leave them RED, and do NOT amend the D-16
real-substrate expectation.** Recorded verbatim, because the reasoning is the point:

> *Phase 13's own two-sentence verdict already says the three-way Bayes factor works on simulated
> data and is NOT SHOWN to work on real microscopy. A test that asserts the real-substrate
> expectation and FAILS is therefore TELLING THE TRUTH. Amending it to expect the new measurement
> would produce a GREEN TEST STANDING NEXT TO A CONCLUSION THAT SAYS "NOT SHOWN" — a test that
> passes while the science says otherwise is worse than a red one, and this project's credibility
> rests on exactly that not happening. **The red is not a defect to be cleared; it is the finding,
> encoded where a future reader will trip over it.***

- **THE RETRACTION IS WRITTEN, PLAINLY.** The manuscript-bound sentence *"negative induced μ is
  CONSTRUCTIBLE from real microscopy pixels via the D-16 mask-based reassignment"* is **RETRACTED**.
  It was an **artefact of segregating a target channel against a nuclear counterstain**:
  `alpha_star_real` went from 0.875 / 0.625 to **`nothing` / `nothing`**, and `m-bar` stays
  **positive at every rung on both fixtures** (+0.4603 → +0.1381, +0.3815 → +0.2041). **The
  correction COST the claim rather than revealing a new one** — the same register as the amendment
  taking the worse number. The converse is **not** asserted either; n = 2 cannot settle whether the
  negative regime is physically reachable, only that it is **not demonstrable on this substrate**.

- **CORROBORATED BY TWO INDEPENDENT ROUTES**, both verified against the *Quick Tasks Completed*
  table in this file and against `spike/simulator/ghat.jl:84-91` at HEAD: **260725-vl8** (`151ad79`)
  **WITHDREW — did not invert** — the unqualified *"negative tail not physically reachable /
  PRIOR-ONLY"* claim; **260725-wb7** (`78dc37f`) found `neg_reachable` flipping to `true` on a masked
  reading of **−0.0338 ≈ 0** to be a **PREDICATE ARTEFACT rather than evidence**. Three routes, one
  direction — and that agreement is what makes this robust rather than one surprising measurement.

- **The D-05 coherence disagreement stays RED on the same principle, and is recorded as an
  INFORMATIVE red.** Contrast **0.09819** inside `tau = 0.15` assigns RANDOM by the phase's own
  label rule while the net's argmax says COLOC — on fixtures flagged OOD at **4.198×** their own
  threshold, i.e. above the maximum of all 96,000 simulated acquisitions the null was fit on. A
  disagreement there is **close to what one should expect**; continued agreement would have been the
  surprising outcome.

- **THE MASKING IS FIXED STRUCTURALLY, and re-ordering was rejected as a fix.** The documented house
  pattern (`runtests.jl`, *"THE CORRECTION ARM IS LAST, AND DELIBERATELY SO"*) **supports exactly ONE
  throwing file**; there are now two, so whichever ran first masked the other and re-ordering would
  only choose which red is hidden. Each known-red include is now wrapped and a **ledger** runs after
  the last one, with four verified properties: every Phase-13 sibling **runs and reports**; the
  named-limit failures **still surface**; the suite **still exits non-zero**; and a known-red file
  that ever **starts passing fails loudly**, naming itself. Anything that is not a
  `Test.TestSetException` is **rethrown immediately** and never recorded as an expected red.

- **MEASURED, one full-suite run each side, 2026-07-29.** BEFORE: aborts at `runtests.jl:223`,
  **4** Phase-13 file-level testsets report (703 pass / 3 fail), **7 files masked**, exit 1.
  AFTER: runs to completion, **11** Phase-13 file-level testsets report (**1,901 pass / 5 fail**)
  plus the ledger's 3 pass, **0 masked**, **exit 1**, on **exactly 5 `Test Failed at` lines** and no
  non-Phase-13 failure. `test_p12_suite.jl` keeps FIRST position and `test_p12_consts.jl` testset 9
  still passes **44 / 44**. **A full-suite green was neither achieved nor claimed — the named limits
  keep it red BY DESIGN.** *Concurrency disclosed:* a live Phase-12 executor landed five commits
  (`986d518` … `406015e`) between the two runs, none touching any file under `spike/test/` or
  `spike/p13/`.

- **NOTHING WAS RELAXED.** The five failing assertions (`test_p13_real.jl:226` ×1 and `:292` ×2;
  `test_p13_correction.jl:245` and `:246`) are **byte-unchanged** — not rewritten, not
  `@test_skip`-ed, not `@test_broken`-ed. No seed, bar, floor, band or allowance moved.
  `P13_ITERATION_ALLOWANCE` **1 of 1, UNSPENT**. `spike/p13/consts.jl` byte-unchanged at git blob
  `d82caea6`. `spike/validation/p12_consts.jl`, `src/`, `spike/Project.toml` and
  `spike/Manifest.toml` byte-unchanged; the Phase-16 seal stays SHUT. Written up in `13-REPORT.md`
  §8c (the retraction), §9 **Limit E**, **§9a** (the five red assertions and the harness fix) and
  §14.9 (the ruling). **`13-17-SUMMARY.md` is now `status: complete`.**

- **2026-07-29 — 13-14 COMPLETE (`0c1a0bd` report + `48fb808` ROADMAP tick-off). THE PHASE IS
  CLOSED.** `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-REPORT.md` (1,086 lines):
  what was CLAIMED, what was MEASURED, what was DECLINED, what is OWED.
  **THE TWO-SENTENCE VERDICT, both halves up front and neither buried:** the three-way amortized
  Bayes factor **works on simulated data** — the amended gate cleared all six pre-registered criteria
  on ONE run, on thresholds frozen before the net existed, without spending the iteration allowance —
  and it is **not shown to work on real microscopy**, because both real fixtures are flagged
  out-of-distribution at ~2.5× their own ID thresholds and the registration-aware basis did not fix
  that. The real arm demonstrates COHERENCE, not CORRECTNESS, on n = 2.
  **The `## Scope of evidence` block is a HEADER element**, positionally asserted to precede section 1,
  carrying the three-arm standing table (exactly ONE gating arm) and the greppable sentence *No
  pass/fail threshold is defined for any real-image quantity.*
  **FOUR NAMED LIMITS**, in the liftable `docs/amortized.md` shape (`docs/` NOT edited): **A** the
  physical segregation anchor is n = 0 in practice (sealed for Phase 16); **B** exclusion is negative
  intensity correlation, not disjoint localization, and α is a CONSTRUCTION parameter so α\*
  CALIBRATES the gap rather than MEASURING it; **C** research-lane, inherits no gate lineage, and the
  shipped artifact is an **epoch-4 checkpoint of an untuned early-overfitting run**; **D** the
  exclusion hypothesis has no labelled real-data validation. A conditional FIFTH entry (the ~2.5 : 1
  prior-atom asymmetry) was **evaluated and NOT added** — its condition was exclusion-end degradation
  and the gate showed the opposite (deep-tail AUC 0.997257 vs 0.979025 nearest τ).
  **WHAT WAS DECLINED is written as first-class content** (ten items): the KDE+quadgk reference,
  continuity as a criterion, MCE as the gate statistic, the confusion matrix as a decision rule, any
  real-data confusion matrix, the unspent iteration, **a third amendment**, the `src/` rename, the
  corpus anchor and CBS, and a third channel pair.
  **EVERY NUMBER WAS RE-READ FROM ITS ARTIFACT AT REPORT TIME**, and that caught things a copy would
  not have. §14 records five discrepancies: (1) **`13-12-SUMMARY.md` §5 is wrong** that the deep-tail
  band is "the best of the four" — measured `band_auc` is `[0.979025, 0.995709, **0.998144**,
  0.997257]`, so the deep tail is the SECOND best; nothing substantive moves and the trigger still did
  not fire; (2) `consts.jl` has two stable fingerprints and they are NOT in conflict — `70fe66df…` is
  the git blob (git normalizes CRLF→LF before hashing) and `5a4ea222…` is `git hash-object
  --no-filters`, the raw working-tree bytes the runners record, both describing one file whose sha256
  is `100e97a3…`; (3) `deferred-items.md`'s `runtests.jl:169` has drifted to **line 189** because the
  concurrent Phase-12 agent wired its aggregator in first — the abort (`test_npe.jl:230`,
  `SPEEDUP_GATE`) is unchanged; (4) the advisory `consts.jl` COMMENT claiming a summary-level crossing
  at α ≈ 0.49 on the positive fixture does not reproduce (measured ≈ 0.674; the negative 0.46 does) —
  a comment, not a bar, and `consts.jl` was NOT edited; (5) `13-06-SUMMARY.md` and `13-07-SUMMARY.md`
  carry **no `status:` frontmatter field**.
  **OWED, with owners:** the `log_bf_simplex` → `log_bf_vs_random` rename at `src/results.jl:172,178`
  (comment-only, zero provenance cost, deliberately NOT done — D-01 forbids a `src/` edit, and until
  it lands the sketch must not be quoted as authoritative); the physical anchor → Phase 16; reshipping
  the net → a productionization phase; corpus extension with real graded segregation anchors → a
  roadmap change; the **F5 confident-tail precision shortfall → Phase 14**, whose abstention layer
  consumes exactly those magnitudes; and the **Phase-4 `SPEEDUP_GATE`** red (measured 68.35 vs bar
  100.0, drifting 92.50 → 83.97 → 68.35, with a 50.402 reading recorded by 13-08) → a **user
  decision**, still logged as `deferred-items.md` D-13-A.
  **ALL EIGHT `13-RESEARCH` OPEN QUESTIONS ARE CLOSED** in §12.
  **INTEGRITY:** `spike/p13/consts.jl` byte-unchanged with only TWO commits in its whole history
  (`c42cc8e` Tier-1, `bfba6ac` Tier-2 τ) — **no commit that produced a Phase-13 result touched it**;
  `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` exit 0; `test/`, `docs/`,
  `corpus/` byte-unchanged; root `Pkg.test()` re-run at report time and **green** (exit 0,
  `Testing ProteinCoLoc tests passed`); the seal stays SHUT; `P13_ITERATION_ALLOWANCE` **1 of 1,
  UNSPENT**. The full spike suite was NOT re-run — it aborts at the pre-existing Phase-4 gate, so the
  plan's own "runtests.jl exits 0" verification is unsatisfiable for reasons predating the plan.
  **STATE.md was written ADDITIVELY BY HAND and no mutating state handler was run**, because
  `## Current Position` carries three concurrent phases and `state.advance-plan` /
  `state.update-progress` would have flattened the Phase-11 and Phase-12 lines while a Phase-12 agent
  is live on this branch. Both commits used explicit paths; the ROADMAP was re-read immediately before
  editing and no race occurred.

- **2026-07-29 — 13-16 COMPLETE (`2112bed` runner + reported artifact).** The amended D-15
  QUALITATIVE real-image arm ran ONCE on the six committed `test/test_images/` TIFFs, on the MAIN
  working tree, no worktree, same ruling as 13-11/13-12/13-13. **This arm consumes NO RNG stream at
  all** — six committed files, a deterministic transform and a deterministic forward pass — so no
  reserved counter was touched and none could be pre-observed.
  **HEADLINE: on BOTH unmodified real pairs, at EVERY pre-registered λ rung, in BOTH read
  directions, the descriptive verdict is RANDOM — and every one of those reads is OOD-FLAGGED.**
  At the headline rung `λ = 3.0` (`:widest_rung`, asserted equal to Phase-11's `LAMBDA_MAX` and to
  the frozen `P13_REAL_LAMBDA_HEADLINE_EXPECTED`): positive-as-sample `log BF(C:R) = −0.4658`,
  `log BF(E:R) = −10.2727`; negative-as-sample `−3.9432` / `−9.0534`. Both log-BFs negative ⇒ the
  `argmax(logBF_C, 0, logBF_E)` rule lands on the RANDOM reference for both specimens.
  **AND THAT IS THE COHERENT ANSWER BY THE PHASE'S OWN LABEL RULE.** Through the frozen `ghat` the
  positive fixture sits at ρ = **+0.3304** and the negative at **+0.2148**, so the D-05 contrast is
  **±0.1156** — INSIDE the measured `P13_TAU = 0.15` dead zone. Both specimens are individually
  correlated but to SIMILAR degrees, and "less colocalized than the control" is not segregation.
  A consistency check that passed — **NOT** a correctness check: there is no ground-truth label.
  **THE OOD FINDING IS THE ARM'S HEADLINE HONESTY ITEM, AND BOTH DETECTORS AGREE.** Phase-13 net:
  density **417.2974** vs its own ID threshold **167.5446** (**2.491×**). Shipped reference:
  **433.6884** vs **179.1368** (**2.421×**), reproducing the frozen `P13_REAL_OOD_SHIPPED_DENSITY`
  / `..._THRESHOLD` (433.69 / 179.14) EXACTLY. **A12 outcome: `agreement`.** So the
  registration-aware Phase-11 basis did NOT move real microscopy back inside the training
  distribution. The real fixtures score above the **maximum** of 96,000 simulated acquisitions
  (ID q50 40.17, q95 167.54, max 290.69). Not suppressed, not softened, not "fixed" — every number
  in every table is printed and persisted on the SAME ROW as its density, threshold and verdict.
  **λ RESPONSE IS ESSENTIALLY FLAT** — 0.0117 / 0.0585 nats on the exclusion head across the whole
  ladder, 0.2595 / 0.1750 on the coloc head — which Phase 11's own flat-in-λ result had already made
  a PRE-AUTHORIZED expected outcome, not a defect.
  **THE REAL α-LADDER, its own arm, NEVER averaged with 13-13's:** `alpha_star_real` = **0.875**
  (positive) / **0.625** (negative), each a RUNG of the frozen ladder, never interpolated. On the
  negative fixture that is EXACTLY the first rung whose induced ρ clears −τ (−0.2193; the rung below
  sits at −0.1450, inside the dead zone by 0.005); on the positive fixture it is ONE RUNG beyond the
  boundary. With n = 2 that is a two-point observation, not a bias estimate. The gap from 13-13's
  `α* = 0.25` is a difference of STARTING POINT (α = 0 is moderately colocalized here, random there),
  not of net quality. All EIGHT transform invariants hold on BOTH conditions, verified BEFORE any
  curve was read (worst intensity residual 2.86e-16, zero absent patches at every rung).
  **REDUNDANCY arm (`P13_REAL_REDUNDANCY_PAIR = (1,3)`):** same qualitative outcome, both log-BFs
  negative at every rung, at a HIGHER density (607.5728, 3.63× over). Corroborates; adds nothing
  independent.
  **NOT A GATE, AND NONE WAS ADDED.** `P13_REAL_IS_GATED = false` / `P13_REAL_QUALITATIVE_ONLY =
  true` were frozen first; the runner has zero `@test` lines (source-grep asserted) and no threshold
  exists for any real-image quantity. `consts.jl` byte-unchanged (git blob `5a4ea222…`, identical to
  13-11/13-12/13-13); `src/`, `test/` and both spike manifests byte-unchanged, `src/` asserted again
  at RUN TIME. **`P13_ITERATION_ALLOWANCE` remains 1 of 1 UNSPENT.**
  **THE SEAL STAYS SHUT.** The sealed-holdout accessor is never called, no path under the sealed
  provenance tree is constructed, the runner contains the string "corpus" **zero** times (not even
  in a comment), and it asserts that about its own source text at load. `test/test_images/` is
  proven byte-identical before and after: digest
  `eaee22f9185460fd910788026e7a33511f6de373f313459b96bb43f3fc1ef276` both times, matching 13-15's.
  **DEVIATION worth carrying:** the plan named `artifacts/grid_8/` as the frozen OOD comparison
  reference; that bundle records `id_threshold = 99.22` and scores these fixtures at 85.04. The
  frozen 433.69 / 179.14 belong to `artifacts/amended_v2/grid_8/`. The runner does not hard-code
  either — it selects the candidate whose OWN recorded threshold reproduces the frozen constant, so
  **the pre-registration is the selector**. Note `artifacts/` is UNTRACKED, so `consts.jl` is the
  durable record and the recomputation is corroboration.
  **NAMED LIMITS for 13-14 to lift verbatim (13-16-SUMMARY §8 and §9):** n = **2 specimens** (six
  files = 3 channels each), so no rate and no coverage claim; no colocalization ground-truth label,
  so BEHAVIOUR and never CORRECTNESS; both specimens OOD-flagged by both detectors; **both**
  pre-registered channel pairs include the DAPI nuclear counterstain — a limit of the frozen
  pre-registration, NAMED rather than repaired by choosing a new pair after the fixtures had been
  measured; and the one physically-segregated anchor stays `sealed_holdout`, deferred to Phase 16.
  Open Questions 6, 7 and 8 are CLOSED (13-16-SUMMARY §10). The manuscript sentence: **negative
  induced μ is CONSTRUCTIBLE from real microscopy pixels via the D-16 construction, but has NOT been
  observed to occur naturally in the images this project holds.**
  Artifact `spike/p13/realimage_report.jld2` (81,649 bytes, 57 keys) is COMMITTED.
  `spike/figures/p13_realimage.png` is **gitignored** under the house `*.png` rule; regenerating it
  costs a ~1.0–1.3 min re-run that is BYTE-IDENTICAL (deterministic substrate), so a loss is a
  compute cost, not a re-seed. **No bulk cache directory was created:** the 48,000-pair training pool
  was opened READ-ONLY shard by shard with the cache layer's `open_or_invalidate` deliberately NOT
  called, and `spike/data/cache/p11/` (54 MB) was neither read nor written.
  Test files run directly (the full suite still aborts at the Phase-4 `SPEEDUP_GATE`, pre-existing):
  `test_p13_real.jl`, `test_p13_alpha.jl` and `test_p13_calibration.jl` all exit 0 — and
  `test_p13_real.jl`'s seal testset is now **34 pass / 0 broken**, because the `@test_skip` 13-15 left
  for this plan's runner is LIVE and passing.
  **13-11's OPEN QUESTION IS NOW CLOSED FOR EVERY ARM:** the epoch-4 best-validation checkpoint
  behaves coherently here too (strictly monotone ladders, no degeneracy, no non-finite value). The
  binding constraint on this arm is **not** the checkpoint — it is the OOD flag.

- **2026-07-29 — 13-13 COMPLETE (`0943a02` runner + reported artifact).** The REPORTED, NON-GATED
  α-graded segregation series ran ONCE on 64 simulated items drawn at `ρ_true = 0` on the reserved
  stream at `P13_ALPHA_COUNTER = 4`, disjoint from both the training pool's counter 2 and the gate's
  counter 3. Executed on the MAIN working tree, no worktree, same ruling as 13-11/13-12.
  **HEADLINE: `α* = 0.25`** — the smallest rung of the frozen `P13_ALPHA_LADDER` at which the mean
  `log BF(exclusion : random)` first exceeds 0. Reported as a RUNG, never interpolated, so the true
  crossing lies somewhere in `(0.125, 0.25]` and this design cannot say where inside it.
  **THE CROSSING IS NOT ARBITRARY — IT LANDS EXACTLY ON THE PRE-REGISTERED HYPOTHESIS BOUNDARY.**
  Mapping the persisted `m̄` through the frozen `ghat`: α = 0.125 induces ρ = **−0.1443**, INSIDE the
  measured `P13_TAU = 0.15` dead zone (by 0.006); α = 0.25 induces ρ = **−0.2499**, the FIRST rung to
  clear it. The evidence head declines below the boundary — correctly, by D-05's own label rule — and
  calls exclusion at the first rung the rule permits, on a spatially-constructed ladder it never saw
  in training. The summary's own sign flip happens ONE RUNG EARLIER (m̄ crosses between α = 0 and
  0.125), and that one rung of lag IS the τ dead zone, not an unexplained gap.
  **THE LADDER** (mean log BF(E:R) / mean log BF(C:R) / m̄, over 9 rungs): −4.074/−3.239/+0.0678 →
  −1.137/−5.694/−0.0182 → **+1.505**/−7.310/−0.1043 → +3.678/−8.068/−0.1860 → +5.446/−8.403/−0.2598
  → +6.843/−8.634/−0.3233 → +7.913/−8.852/−0.3764 → +8.725/−9.061/−0.4197 → +9.345/−9.250/−0.4546.
  Per-image positive rate for the exclusion head rises **0% → 27% → 86% → 98% → 100%** and the
  per-rung spread NARROWS 1.88 → 0.82 nats. The coloc head is negative on 64/64 items at every rung
  above α = 0.
  **THE LADDER IS PROVABLY NOT SHAPED BY A MECHANICAL ARTIFACT.** All EIGHT transform invariants hold
  on ALL 64 REPORTED images (not merely on a fixture), verified BEFORE any curve was read: bitwise
  identity at α = 0, no NEW zero, intensity conserved (worst residual 9.63e-16 against an rtol bar of
  1e-10), mask α-invariant, mask fraction 0.284–0.343 inside the (0.01, 0.40) band. Pitfall-2 signs:
  **ZERO absent patches at every rung** (0 of 4,096 at both ends), the count did not rise, and the
  summary's MASK ROWS are **byte-identical across every rung on every image**.
  **NO GATE, AND NONE WAS ADDED.** `P13_ALPHA_GATED = false` was frozen before any Phase-13 result;
  the runner carries no assertion on any curve (asserted by source grep, 0 occurrences).
  `consts.jl` byte-unchanged (git blob `5a4ea222…`, identical to 13-11's and 13-12's records);
  `src/`, `corpus/` and both spike manifests byte-unchanged, `src/` asserted again at RUN TIME.
  **`P13_ITERATION_ALLOWANCE` remains 1 of 1 UNSPENT.**
  **NAMED LIMIT — the physical anchor stays deferred and the seal stays shut.** The one real
  segregated anchor (`neg-lightmycells-01`) is triply unavailable — `sha256 = "PENDING-FETCH"`, no
  bytes on disk (`bytes = 0`), and `split = sealed_holdout` behind the anti-snooping accessor. Phase
  13 did NOT open it: that control exists for Phase 16's BLIND evaluation. No download-size claim
  attaches (assumption A6 withdrawn — Phase 13 fetches nothing). Nothing under `corpus/` was checked,
  fetched or referenced; the runner contains the string "corpus" **zero** times, not even in a
  comment. 13-13-SUMMARY §6 carries the deferral as a paragraph 13-14 can lift verbatim.
  **SCOPE: this is the `:simulated` arm ONLY.** The `:real` arm is 13-16's, over the six committed
  `test/test_images/` TIFFs, and the two are **NEVER averaged** — α = 0 means "random by
  construction" here (both acquisitions at ρ_true = 0) and "moderately colocalized" there (m̄ +0.3292
  / +0.2481). The CBS `cbs-RG-000` arm is **DROPPED**, not deferred.
  **DECLARED DEVIATION worth carrying:** the simulated arm declares an **UNBOUNDED** dynamic range
  rather than the frozen `P13_ALPHA_MAX_VALUE_BOUND = 1.0`, which is a Gray-TIFF property of the REAL
  substrate; simulator intensities measured **4.08–17.96** here. `max_value` is a caller-supplied
  substrate property by `alpha_series.jl`'s own design, so nothing is clamped and no new bar was
  invented — but the consequence is that the `max_value_ok` invariant is **vacuously true on this
  arm** and is never quoted as evidence. Realized maxima are recorded per image instead.
  **13-11's OPEN QUESTION IS NOW ANSWERED FOR 13-13 TOO:** the epoch-4 best-validation checkpoint is
  sufficient here as well — monotone means on both heads, 100% positive rate from α ≥ 0.5, narrowing
  spread. It remains open only for 13-16.
  **SCOPE LIMITS on α\*:** α is a CONSTRUCTION parameter, not a physical quantity, so α\* CALIBRATES
  the correlation-versus-localization gap rather than MEASURING it; the ladder spans induced ρ
  [−0.705, −0.026] and therefore does **not** probe the −0.99 atom (that tail was cleared separately
  by 13-12's gate at AUC 0.997257); and the result is specific to one mask rule and one
  redistribution rule. The exclusion hypothesis still has NO observed real-data instance anywhere in
  this project.
  Artifact `spike/p13/alpha_report.jld2` (53,391 bytes, 65 keys) is COMMITTED and carries the three
  curves, the per-image scores, the invariant residuals, α\*, the seed/counter and the consts
  fingerprint. `spike/figures/p13_alpha_ladder.png` is **gitignored** under the house `*.png` rule;
  regenerating it costs a ~2.5 min re-run that is BYTE-IDENTICAL (counter-based Philox per index), so
  a loss is a compute cost, not a re-seed. **No bulk cache directory was created** and
  `spike/data/cache/p11/` (54 MB) was neither read nor written.
  Test files run directly (the full suite still aborts at the Phase-4 `SPEEDUP_GATE`, pre-existing):
  `test_p13_alpha.jl`, `test_p13_real.jl` and `test_p13_calibration.jl` all exit 0 — the latter two
  matter because both discover `spike/p13/` with `readdir` and source-grep every `.jl` in it,
  including the new runner.

- **2026-07-29 — 13-12 COMPLETE (`e5f9569` runner + reported artifact).** The amended D-12/D-13
  gate ran ONCE on a FRESH 4,000-pair evaluation set at `P13_GATE_COUNTER = 3`, disjoint from the
  training pool's counter 2, and **PASSED all six pre-registered assertions** (exit 0, 6/6).
  Executed on the MAIN working tree, no worktree, same ruling as 13-11.
  **MEASURED, against thresholds byte-locked before the net existed:**
  per-class one-vs-random AUC **0.990169** (coloc) and **0.988262** (exclusion) vs floor **0.9**;
  per-head ECE **0.0119808** and **0.012886** vs green band **0.05**, both `:green`, **0 of 10
  empty bins on both heads**, MCE 0.0816278 / 0.116682 (reported, never gated); **neither head
  vacuous** (both AUCs far above `P13_VACUOUS_AUC_FLOOR = 0.6`). Realized classes E 1334 / R 1333 /
  C 1333 from 6,322 draws (overhead 1.58×); restricted per-head sets 2,666 / 2,667, both far above
  `P13_MIN_EVAL_PER_HEAD = 1000`. Generation 10.86 min.
  **The DESCRIPTIVE argmax confusion matrix** (rows true, cols predicted, order E/R/C):
  `[1265 69 0; 101 1148 84; 0 66 1267]` — ZERO exclusion-called-coloc and ZERO coloc-called-
  exclusion; every error is an adjacent confusion with the shared `random` reference. It carries
  no assertion: decisions and abstention are Phase 14's.
  **D-04 ITERATION TRIGGER: DID NOT FIRE, evaluated mechanically.** The exclusion AUC restricted to
  the deep tail |ρ| > 0.9 is **0.997257** (n = 239) — not below its floor, and in fact the BEST of
  the four |ρ| bands (0.979025 / 0.995709 / 0.998144 / 0.997257 over (0.15,0.5] / (0.5,0.75] /
  (0.75,0.9] / (0.9,1.0]). Discrimination IMPROVES toward the −0.99 atom, which is the opposite of
  the failure the allowance was reserved for. **`P13_ITERATION_ALLOWANCE` remains 1 of 1 UNSPENT
  and no retraining is authorised.** `consts.jl` byte-unchanged (git blob `5a4ea222…`, identical to
  13-11's record); `src/` and both spike manifests byte-unchanged, asserted at RUN TIME by step 0.
  **λ RESPONSE (D-03, reported not gated) — FLAT *AND* PERFECTLY MONOTONE, and the distinction
  matters.** Across Phase 11's own frozen `SC2_RUNGS` ladder (0.25→3.0), `mean|logBF|` moves
  7.36726→7.38413 (coloc head) and 8.26806→8.33566 (exclusion head): spans of **0.017** and
  **0.068 nats**, i.e. **0.2% and 0.8%**, with **Spearman 1.0 on both heads** (all seven rungs in
  order, in opposite directions). So "the λ input has no effect" would be FALSE — it is wired, read,
  and monotone; the 8×8 summary simply does not carry enough registration information for it to
  move the evidence materially. This is Phase 11's own negative result one level up, pre-authorised
  as a reportable finding, and a future summary redesign has something to move rather than a dead
  input.
  **BINARY-NRE CONTINUITY (reported, NEVER gated, `P13_CONTINUITY_GATED = false`):** over the 2,666
  coloc/random pairs at the reference λ = 3.0, corr **0.62216**, max|Δ| **22.02**, mean Δ **−4.21**.
  Read against the SPIKE-LANE frozen binary NRE (`spike/validation/trained_ratio.jld2`), not
  `src/amortized/bf.jl` (which imports KernelDensity/QuadGK, absent from the lean spike env) and not
  `artifacts/grid_8/` (reserved by `consts.jl` §I2 for the OOD comparison only). The summaries were
  round-tripped EXACTLY between the two frozen z-score bases, so both nets saw the same
  acquisitions; the offset is what D-02 and the two-different-logits argument predict.
  **13-11's OPEN QUESTION IS ANSWERED FOR THESE TWO CRITERIA:** the persisted epoch-4
  best-validation checkpoint of an early-overfitting run IS sufficient for D-12 discrimination and
  D-13 calibration — comfortably, not marginally. It remains open for 13-13 and 13-16.
  **DISCLOSED so the ordering is checkable:** a load smoke at `m = 90, min_per_head = 30` was run on
  the gate stream before the reported run to prove the runner executes. It writes to a `*_smoke`
  path, prints `SMOKE MODE -- NOT A REPORTED RUN` and does NOT run the gate assertions; both smoke
  outputs were deleted. **Not one byte of the runner changed afterwards.** Its ECE read yellow
  (0.0738 / 0.0789 at 60 items/head vs 0.0120 / 0.0129 at 2,666) — recorded as evidence that
  `P13_MIN_EVAL_PER_HEAD` is load-bearing, with the caveat that the observed direction is the
  OPPOSITE of the downward small-n bias the pre-registration warns about (per-bin noise dominates
  when ten bins share sixty samples).
  **WHAT THE PASS DOES NOT ESTABLISH:** simulator ground truth only; well-specified regime
  (train-joint == eval-joint by F5); not a decision rule; and the exclusion hypothesis still has NO
  observed real-data instance anywhere in this project. SC3's other two arms (13-13 α-ladder,
  13-16 real images) are explicitly not gates.

- **2026-07-29 — 13-11 COMPLETE (`133558d` tests → `41a24f4` trainer + artifact; Task 1 was
  pre-committed at `4839c66` by a predecessor that died on a session limit).** Executed on the
  MAIN working tree with NO worktree, deliberately: the pool lives in gitignored
  `spike/data/cache/p13/` and a worktree cleanup is exactly what destroyed Phase 11's 54 MB pool.
  **MEASURED:** realized class frequencies EXACTLY equal thirds (whole-class accept/reject to an
  integer quota of 16,000 each; acceptance overhead 1.545×, 48,000 accepted from 74,149 draws).
  Frozen-zt check **PASS / verdict `inherited`** — max|per-row mean| 0.057353463, max|per-row sd-1|
  0.035997152 over 96,000 columns × 64 continuous rows, with the gate's re-fit negative control
  confirming a genuine re-fit lands at machine zero (<1e-8) while the inherited residual is >1e-4.
  Per-head corrections MEASURED on the training subset only: **coloc +0.0029420, exclusion
  +0.0023543**, against π-level −0.7765 / −0.4090 (π masses E 0.31272 / R 0.47073 / C 0.21655,
  matching 13-RESEARCH E1 to ~3 dp) — stratification removes the prior imbalance as designed, and
  the values are measured rather than assumed because the contiguous split does not divide a
  stratified pool into exact thirds.
  **HONEST FINDING, NOT REPAIRED — the run overfits early.** Best validation risk **0.13625 at
  epoch 4**; validation then rose monotonically to 0.38632 while training risk fell to 0.03655
  (~10× gap), early stopping at epoch 45 of 300. The persisted artifact IS the epoch-4 best
  checkpoint (NeuralEstimators `train` returns `trainstate_best` on CPU, `train.jl:643`), so the
  overfitting cost usable budget, not shipped weights. The verbatim recipe was NOT tuned;
  `consts.jl` byte-unchanged, `P13_ITERATION_ALLOWANCE` still 1 of 1 **UNSPENT**. Whether a
  best-at-epoch-4 checkpoint suffices is left to the downstream gates (13-12/13/14/16), because
  answering it by retraining is what the allowance forbids absent its pre-declared trigger.
  **ASSUMPTION DRIFT worth carrying:** training took **2.37 min**, pool generation ~120 min — the
  spend is ~98% datagen, so the marginal cost of a second training run is minutes. The iteration
  allowance is therefore a matter of pre-registration integrity, NOT compute budget; it is cheap,
  and that is not the reason it is forbidden.
  Artifact `spike/p13/three_way_net.jld2` (953,027 bytes) carries head_log_odds, τ + its Tier-2
  provenance shas, the Phase-11 path + sha256, realized/target class freq, the consts sha256 AND
  git blob sha, the full recipe, seed/salt/counter and pool dir. Gate
  `spike/test/test_p13_datagen.jl` 706/706 with ZERO skips; root suite green; `src/` and both spike
  manifests byte-unchanged.
  **The 50 MB pool at `spike/data/cache/p13/7c65a1a9…/` is GITIGNORED and unbacked.** Regenerating
  costs ~2 h CPU but is BYTE-IDENTICAL (counter-based Philox per global index), so a loss is a
  compute cost, not a re-seed. Resume-by-skip was exercised for real twice (5 predecessor shards
  skipped; a mid-run process kill at 10 shards lost nothing).

- **2026-07-28 — 13-10 closed out from its committed artifact (`85a75c2`).** The prior executor
  landed all three 13-10 commits (`581602a` artifact → `bfba6ac` constant → `45f5a44` tests) in a
  worktree and died on a session limit before writing the SUMMARY. Resolved via the
  `safe_resume_gate` "close out manually" path rather than re-dispatch, because re-running the
  probe would be a SECOND reported measurement of a pre-registered quantity and would have to be
  declared as spending the D-04 iteration allowance reserved for the stratification switch.
  **MEASURED `P13_TAU = 0.15`** at A(τ) = 0.916356 against the pre-registered bar 0.9, reference
  λ = 3.0 (widest rung). Grid NOT extended, bar NOT relaxed, D-04 allowance UNSPENT. τ now flows
  into 13-11/12/13/14/16 and must always be quoted with its reference λ.

- **2026-07-28 — 13-09 ORCHESTRATOR RULING: Phase 13 binds the Phase-11 research NPE even though
  Phase 11 closed NEGATIVE.** 13-09's checkpoint asks whether Phase 11 "landed". Taken literally
  it FAILS: Phase 11 is CLOSED not COMPLETE, `11-07-SUMMARY` is `status: blocked`, and 11-08…11-11
  are SUPERSEDED. Proceeded anyway, on these grounds, recorded here and in
  `spike/p13/preconditions.jl` so the decision is auditable rather than implicit:
  (1) Phase 11's own closure records SC1g as **MIS-SPECIFIED, not failed** — the bar 2.502 was
  derived from a component ratio but applied to a total, and the true width ratio is ~1.00, so
  "the net was right and the gate failed it anyway"; (2) the λ conditioning is measurably ALIVE
  (shift marginals track λ at SD ratios 4.84–9.05); (3) all six items of 13-09's own artifact
  contract were empirically verified present before dispatch; (4) 13-12 already pre-authorises a
  FLAT λ response as a reportable finding, so proceeding does not launder a negative result;
  (5) D-02 claims Phase 13 trains on Phase 11's registration-aware INPUT SURFACE, not that
  registration is inferable — so Phase 11's negative result leaves D-02 intact.
  **This is distinct from the Phase-12 premise flagged below, which remains NOT acted on.**
  Bound: `spike/npe/p11_research_npe.jld2`, frozen `zt` (64 continuous rows), frozen
  `P11BoundedThetaTransform` (arity 8, chromatic ε 8th at ±0.02), `n_cond = 1` DERIVED via
  `length(encode_lambda(λ))`, λ ∈ (0.25, 3.0), F5 imsize provenance. No grid-8 fallback exists.

- **2026-07-28 — two defects found and fixed while landing 13-09, both worth carrying forward.**
  (a) `spike/p13/consts.jl` reserves the name `P11_DEV_SEED`, which is exactly the guard sentinel
  `spike/validation/p11_consts.jl` keys its whole Tier-1 block on — so loading the Phase-13
  pre-registration first made `p11_consts.jl` skip its own block and die on its out-of-guard
  self-check. Repaired by reading it through a private `module _P11C` (the existing `module _GC`
  precedent); neither frozen consts file was edited. **General lesson: "guarded include ⇒
  order-free" is FALSE whenever two frozen consts files reserve the same name.**
  (b) `test_p13_preconditions.jl` drew its `assert_frozen_zt` fixture pools from an UNSEEDED
  `randn()`, while the verdict is decided from the pool's measured moments — an intermittent
  suite-reddening flake (observed 72/74 once, then 74/74 four times). Reseeded from
  `p13_fix_rng(P13_FIXTURE_COUNTER)` (`b9b688e`); 74/74 on five consecutive runs.

- **KNOWN RED, pre-existing, NOT a regression:** `spike/test/runtests.jl` exits non-zero at
  `test_p13_correction.jl` on exactly 2 assertions — `maxabs` 0.5498 / 0.3783 against
  `P13_F5_MAXABS_TOL = 0.25`. Plan 13-07 committed these deliberately as an honest finding
  ("Recorded as an honest finding, not repaired"; its own verification table logs "exits 1: 88
  pass, 2 fail"). Correlations pass comfortably (0.99751 / 0.9988 vs bar 0.99). Because this
  abort sits at `runtests.jl:213`, it masks everything after it — a full-suite green is therefore
  NOT available as a gate signal for the rest of Phase 13, and per-file runs must be used instead.
  **Re-deriving that bar is a pre-registration decision and is left for the user.**

- **CORRECTION to the note above, measured 2026-07-29 during 13-11: the suite never REACHES line
  213.** It aborts EARLIER, at `runtests.jl:169` → `test_npe.jl:230`, on the Phase-4 `SPEEDUP_GATE`
  (median speedup **68.35** against the pre-registered bar 100.0 — the blocker already recorded
  below from the Phase-11 wave-3 post-merge gate, where it measured 92.50× and 83.97×, so the
  number has drifted further down). Consequence: the masked region is the **ENTIRE Phase-13 include
  block** (all eleven `test_p13_*` files), not merely what follows the correction arm. The
  conclusion is unchanged and now has two independent causes — per-file runs are the gate signal
  for Phase 13. NOT fixed by 13-11: `SPEEDUP_GATE` is a pre-registered threshold (user decision)
  and `test_npe.jl` is Phase-4 code that plan does not own. Detail:
  `.planning/phases/13-three-hypothesis-amortized-bayes-factor/deferred-items.md`.

Phase: 11 (registration-and-chromatic-uncertainty-as-latent) — CLOSED 2026-07-27
Plan: 7 of 11 executed; 11-08 through 11-11 SUPERSEDED by the diagnosis (reasoning per plan in
`11-CLOSURE.md`), not silently skipped.
**STATUS: Closed on a negative-but-useful result. Registration at ≤3 px is NOT inferable from the
8×8 patch-correlation summary — not at finer resolution, not with a self-referential probe feature,
and not at any colocalization strength (ridge/prior 0.993–1.001 across five |ρ| bins, n ≥ 1322).
It also does NOT need to be: the width the data demands is flat in λ (RMSE ratio 1.0003), so the
posterior does not report false confidence. The SC1g gate is recorded as MIS-SPECIFIED, not failed
— it was derived from a component ratio (`lambda_ratio`) but applied to a total, and the net's
measured 1.109/0.987/0.936 sits against a correct answer of ~1.00. No threshold edited,
`P11_ITERATION_ALLOWANCE` still 1 of 1 UNSPENT, shipped `amended_v2/grid_8` untouched.
Evidence: `11-DIAGNOSIS.md`. Decisions and limits: `11-CLOSURE.md`.
FLAGGED FOR THE USER, NOT ACTED ON: Phase 12 was scoped to train on "registration-aware θ"; that
premise no longer holds as stated and should be revisited before Phase 12 is planned.**

- **2026-07-25 — Phase 11 probe verdict (plan 11-06): ABOVE RESOLUTION, DECISION `proceed`.** The D-06 abort criterion did not fire: `S_probe = 1.0` against `P11_PROBE_S_FLOOR = 0.9`, and a Δρ_eq ladder span of 0.046885 against `P11_PROBE_SPAN_FLOOR = 0.02` (F5 mixture arm; both independently re-derived from `p11_probe_report.jld2`, agreeing to 1e-16). Plans 11-07 through 11-11 are authorised as planned. `P11_ITERATION_ALLOWANCE = 1` is UNSPENT (0 of 1); no threshold was altered. The ~1.5 h datagen-plus-training spend is authorised, and `P11_DATAGEN_WALLCLOCK_CEILING_MIN = 150` is a BLOCKER threshold, not a downgrade trigger. Checkpoint resolved autonomously by the orchestrator under the pre-registered rule (background session, no user present), conditional on the re-derivation agreeing — which it did. See `11-PROBE-VERDICT.md`.

**Phase 11 (Registration and Chromatic Uncertainty as Latent) — IN PROGRESS.**
Research-lane only (D-01): a Phase-11 net trains in `spike/` on a fresh DEV seed; the shipped
`amended_v2/grid_8` artifact and the Phase-7 GO (Option A, 2026-07-24) stay untouched — no reship,
no new artifact, no ship gate reopened, no public API change. ROADMAP SC1 is already nominally
satisfied (`shift_dx`/`shift_dy` are existing latents), so the phase owns three genuine gaps: the
missing chromatic term, the never-run SC2 monotone-widening sweep, and SC3.
Deliverable (D-16) = a results report under `.planning/` + an interpretation note in
`docs/amortized.md`.

Three ordering constraints are load-bearing and encoded in the wave graph:
(1) the D-06 simulator-only pre-flight probe runs strictly **before** any datagen/training and
supplies the SC2 threshold; (2) the D-11 `ε` mirror into `src/amortized/simulator.jl` and the D-12
named limit #8 in `docs/amortized.md` land in **one commit**; (3) the pre-edit golden fixture and
the pre-`ε` commit sha are captured **before** the first simulator edit.
Plan 11-06 is a blocking decision checkpoint (probe verdict) gating ~1.5 h of CPU spend and the
single declared D-04 iteration allowance.

Eight researcher open questions were resolved as orchestrator rulings at plan time (recorded in each
implementing plan's `<context>`): D-17 retargets from `corpus/` (sealed-holdout, bytes absent) to
`test/test_images/`; D-11 mirrors `ε` only (src `SHIFT_PRIOR` stays `Uniform(-1.0, 1.0)`); one λ
conditioning input, so SC2's widening claim is registration-only; `BoundedThetaTransform` is ported
into the spike trainer as a named deviation; breakdown curve primary at λ = 3.0, secondary at λ = 1.0;
no new REQ-IDs; train and evaluate on the F5 image-size mixture (Tier-1 locked); and the narrower
decoupling claim replaces ROADMAP.md:126's literal `git status`-clean-on-`src/` wording.

**Phase 07 (productionization) — COMPLETE (11/11 plans).**

**07-09 (2026-07-24) — CAP DECISION (no training):** 32×32 is NOT trained and NOT shipped
(`gate-32x32.md`, `07-09-SUMMARY.md`). The GO rests on 8×8, and a 32×32 gate would run the KDE-BF /
non-atom-corrected apparatus spikes 006-014 proved defective (would FAIL like 16×16, contribute
nothing). User decision: cap-the-family. No `_train_grid_pipeline(32)` run, no `grid_32` artifacts, no
`gate_consts_32.jl`, no `PROD_SEED[32]` consumed. 64×64 stays dropped.

**07-10 (2026-07-24) — productionize the amortized-only public API:** shipped
`colocalization_amortized(img, control, channels; num_patches=8, N=2000)` → `AmortizedColocResult`
(`src/amortized/api.jl`), composing the frozen read surfaces (registry → summary → NPE posterior → Δρ
cloud → NRE log-BF → density-channel OOD), CPU-default (D-06), with the D-02 accessors dispatching on
the result. **`_SHIPPED_GRIDS = (8,)`** (4 never post-hoc re-analysed, 16 gate-FAILED, 32 CAPPED, 64
dropped). The shipped 8×8 bundle (`artifacts/amended_v2/grid_8/`, the retrained bounded-θ /
realistic-imsize net the GO rests on) loads via a content-hashed **`Artifacts.toml`** entry
(`git-tree-sha1 = 90e6b63a…`, `lazy = true`, Release-hosted; tarball sha256 `17162904…`), and
`_lazy_load_from_artifact!(8)` resolves it store→dev-fallback→download, **tree-sha1 verified before
load** (T-7-01), then wires the recorded density-channel OOD operating point. `Artifacts`+`Pkg`
stdlibs added to `[deps]` (resolve clean; NeuralEstimators 0.2.1 / Flux 0.16.10 unchanged). Turing
path stays internal (not exported). `test/test_integration.jl` exercises the full 8×8 public path +
error paths; `docs/amortized.md` carries the API and the §5 named limits verbatim.

**Two carry-forward notes:** (1) the shipped OOD flag is the **density Mahalanobis channel** only —
the fused noise/PP channels that reached gate AUC 1.0 are gate-time constructs not in the frozen
bundle (documented). (2) `AmortizedColocResult.calibration` is the `_bundle_from_artifacts`
placeholder (`gated=false`); gate provenance + named limits live in `docs/amortized.md` and the memo.

**07-08 (2026-07-21):** the windowed sub-tile local colocalization map shipped
(`src/amortized/local_map.jl`, commits 25037e6 / 90f612e). `local_coloc_map(img, control,
channels; grid=8, tiles=(r,c))` cuts both images into r×c sub-tiles with the package's own
`patch()` tiler and scores each tile against its control tile through the FROZEN, already-gated
bundle — no training, no new gate, so the map inherits the 8×8 gate verdict verbatim (residual
FAILs included). Because 07-10 has not wired Artifacts yet, `_ensure_grid_registered(grid;
artifact_dir)` bridges in-process from the local `artifacts/grid_G/*.jld2` via
`_bundle_from_artifacts` + `register!` (touches neither `_lazy_load_from_artifact!` nor
Artifacts.toml). Degenerate tiles (below the grid size, zero survivors of the ≥15 floor, or a
non-finite read) take `LOCAL_MAP_SENTINEL = 0.0` AND a forced OOD flag — never a crash (T-7-04).
Phase-12 `SpatialColocResult`/`delta_rho_map`/`uncertainty_map` remain sketch-only in
`src/results.jl` (machine-asserted absent from the new module's executable code). `Pkg.test()`
green: windowed-map 30/30 plus a conditional real-8×8 smoke 5/5 (zero sentinels at 2×2 on 256²).

**Two carry-forward notes for 07-10:** (1) `_ensure_grid_registered` is a NO-OP when the grid is
already registered — an explicit registration always wins over a disk artifact, so superseding one
needs explicit invalidation, not a second `_ensure_*` call. (2) `_train_grid_pipeline` persists
ood_nulls with the `:density` channel ONLY (no `:noise`, no `:model`, **no `:thr`**), so
`ood_verdict` currently computes a score but has no threshold and returns `flag = false`; per-tile
OOD flagging is structural (degenerate tiles) only until the fitted `thr`/`zref` are persisted.
The 16×16 ship-gate RAN to completion (~16.5 min, `gate_report_16.jld2`, `status = :ran`) against the
byte-locked `gate_consts_16.jl` (7e2318b) on the fresh disjoint `PROD_SEED[16] = 0xb906f369f6cacf91`.
Recorded verdict is an **honest FAIL** (`gate-16x16.md`, commits 7f093a1 / e8f6a9d): SBC
`ks_pass=false` / `ece_pass=true` with **only 4/8 KS pass and the headline ρ_true itself rejected at
p = 7.0e-4**; BF corr 0.9150 < 0.95 AND max|Δ logBF| 12.72 > 0.5 (n=15 finite of 25); OOD
**NOT RUN / INCONCLUSIVE** (`auc = nothing`, `passed = nothing` — no `pos_sim` injected; ID operating
point 378.44 only, explicitly NOT a pass). No pre-registered constant was tuned.
**07-10 consequence:** grid 16 is recommended **NOT eligible** for default registry population — the
ρ_true rejection has no documented non-method cause, unlike the 8×8/4×4 residuals.

**Phase 08 (external-physical-ground-truth-corpus) — COMPLETE (2026-07-21).** 5 of 5 plans done.
Wave 5 (08-05) is unblocked and executed: both physical anchors human-verified, pinned, and sealed.
Offline gate `julia --project=. corpus/test/runtests.jl` → 217/217, ZERO network; `src/` provably
untouched; `git ls-files corpus/data` empty. NOTE: this Current Position was previously clobbered by a
parallel Phase-8 run — Phase 7, not Phase 8, is the active phase.

Last activity: 2026-08-04

Progress: [█████████░] 90%

## Resolved (2026-07-03): 07-00 CO-RESOLUTION GATE — GREEN

The Phase-7 Wave-0 co-resolution HARD GATE (07-00 Task 3) initially failed on a SECONDARY
GLMakie/Makie conflict (NOT the Finding-1 NeuralEstimators downgrade, which the Turing→weakdep
surgery had already cleared). Root cause: `GLMakie = "0.10.5"` forced Makie 0.21.18, while
NeuralEstimators 0.2.1 requires Makie 0.24.x.

**User decision (Option 2):** bump GLMakie compat (`0.10.5 → 0.13`, Makie 0.24.x) and keep
GLMakie as a normal CORE dependency (NOT an extension). Applied; the root `Pkg.resolve()` is now
GREEN and `Pkg.test()` passes fully.

Resolved versions (regenerated root Manifest.toml): NeuralEstimators **0.2.1** (not downgraded),
Flux **0.16.10**, Makie **0.24.12**, GLMakie **0.13.12**. Turing/CUDA remain weakdeps; the
`ext/ProteinCoLocTuringExt.jl` extension is unchanged. No `src/plot.jl` API changes were needed
at precompile/load (Makie 0.24 runtime plotting correctness deferred with the plotting tests).
spike/Project.toml + spike/Manifest.toml provably UNTOUCHED throughout.

## Performance Metrics

**Velocity:**

- Total plans completed: 22
- Average duration: - min
- Total execution time: 0.0 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01 | 4 | - | - |
| 02 | 4 | - | - |
| 03 | 5 | - | - |
| 04 | 7 | - | - |
| 06 | 2 | - | - |

**Recent Trend:**

- Last 5 plans: -
- Trend: -

*Updated after each plan completion*
| Phase 01 P02 | 22min | 2 tasks | 4 files |
| Phase 01 P03 | 9 | 2 tasks | 2 files |
| Phase 01 P04 | 13min | 2 tasks | 1 files |
| Phase 02 P01 | 15min | 2 tasks | 5 files |
| Phase 02 P02 | 11min | 2 tasks | 4 files |
| Phase 02 P03 | 32min | 2 tasks | 5 files |
| Phase 02 P04 | 22min | 2 tasks | 5 files |
| Phase 03 P01 | 8min | 2 tasks | 4 files |
| Phase 03 P02 | 12min | 3 tasks | 5 files |
| Phase 03 P03 | 26min | 3 tasks | 6 files |
| Phase 03 PP04 | 18min | 2 tasks tasks | 2 files files |
| Phase 05 P01 | 22min | 3 tasks | 8 files |
| Phase 05 P05-02 | 19min | 2 tasks | 7 files |
| Phase 05 P05-03 | 70min | 2 tasks | 2 files |
| Phase 05 P05-04 | 60min | 3 tasks | 4 files |
| Phase 06 P01 | 40 | 3 tasks | 1 files |
| Phase 06 P02 | 20min | 2 tasks | 2 files |
| Phase 07 P07-02 | 45min | 3 tasks | 8 files |
| Phase 07 P07-03 | 50min | 4 tasks | 9 files |
| Phase 07 P07-04 | 55min | 3 tasks | 7 files |
| Phase 7 P6 | 40min | 3 tasks | 4 files |
| Phase 07 P07-08 | 55min | 2 tasks | 4 files |
| Phase 14 P09 | 35m | 2 tasks | 3 files |
| Phase 14 P10 | 25m | 2 tasks | 1 files |
| Phase 14 P11 | 35 min | 2 tasks | 1 files |

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Roadmap]: 7-layer horizontal build order; ENV smoke gate (Phase 1) is the hard gate before all downstream investment (pre-1.0 NeuralEstimators API is the top unknown)
- [Roadmap]: SBC/BF/OOD grouped into one Phase 5 — sibling plans sharing one θ*~π→simulate→infer harness over the trained nets
- [Roadmap]: Phase 7 (productionization) is conditional on a Go decision in the Phase 6 memo and is the ONLY phase that edits `src/`; spike phases (1–6) keep `src/` provably untouched
- [Phase 01-02]: Spike env is minimal+isolated (NeuralEstimators 0.2.1, Flux 0.16.10, Distributions 0.25.128, no CUDA); CPU-only NPE smoke recovers theta within tol=0.3 under a re-runnable stdlib Test gate. CUDA-absence asserted via Pkg.dependencies (installed set) not a Manifest regex, since Flux/NNlib/Zygote/NeuralEstimators declare CUDA as inert weakdep extensions.
- [Phase ?]: [Phase 01-03]: Spike env frozen as reproducibility artifact — Julia 1.12.6 pinned via juliaup directory override (enforcing) + spike/.julia-version (doc-only); Manifest pins NeuralEstimators 0.2.1 + Flux 0.16.10, no top-level CUDA. NeuralEstimators default / BayesFlow fallback recorded in spike/NOTES.md (ENV-03, ENV-04).
- [Phase ?]: [Phase 01-04]: include() fallback chosen over Pkg.develop for ENV-01 — the parent's GLMakie/Makie 0.21 + Turing 0.44 + GraphNeuralNetworks 1.1.0 tree caps NeuralEstimators <0.2.1, silently downgrading the pinned 0.2.1 to 0.1.4 and breaking the v0.2.1-API smoke; later phases reach src/ read-only via include('../src/<file>.jl'). Package boundary deferred to Phase 4 (D-02). Decoupling proven byte-identical to baseline f581d95.
- [Phase ?]: SIM-03 summary contract proven on a synthetic strictly-positive image via UNCHANGED src/ patch()/correlation() at fixed 8x8 (D-10) before any physics exists
- [Phase ?]: CairoMakie/Images/ImageFiltering/HypothesisTests/StatsBase co-resolve cleanly with NeuralEstimators 0.2.1 (no downgrade, no CUDA); resolve-risk gate automated in runtests.jl
- [Phase ?]: Phase 2 simulator: shared-latent standardized-smooth-field generator (D-15) with sign(rho) flip drives monotone, sign-correct induced patch-correlation (-0.43 to +0.60)
- [Phase ?]: SIM-02 (02-03): frozen monotone ĝ (isotonic PAVA + clamped piecewise-linear inverse) maps μ↦ρ_true; induced μ matches Turing μ-prior W1=0.052<0.10 over realized range [-0.68,0.847]; negative μ tail prior-only (real anchor pos=0.33/neg=0.25)
- [Phase ?]: SIM-04 Spearman monotonicity threshold fixed at 0.95 (achieved 1.0); perturbation effects asserted as paired shared-seed differences with OneSampleTTest p<0.05 (T-02-GATE)
- [Phase ?]: CairoMakie plausibility figures are spike-local headless (D-13); .gitignore scoped negation tracks spike/figures/plausibility.png as the evidence deliverable
- [Phase ?]: [Phase 03-01]: JLD2 (cache backend, D-04) + Random123 (counter-based seeding, D-11) added to isolated spike env; JLD2 already transitive (promoted to direct), only Random123 v1.7.1 + RandomNumbers v1.6.0 newly installed; NeuralEstimators stays pinned v0.2.1 (no co-resolve downgrade); resolve-risk gate extended to cover both new deps; five-SC MISSING scaffold (test_data_pipeline.jl) wired into single runtests.jl gate
- [Phase ?]: [Phase 03-02]: DATA-01 generation core — encode_d01 (128-dim: 64 imputed corr + 64 binary mask, fully-missing kept per D-13), encode_aug (AUG_DIM=142, 14 moments, D-02), Philox4x per-sample keyed RNG keyed by (master_seed,idx) with disjoint HOLDOUT/FOLD salts (D-10/D-11), cost-aware imsize sampler (>=1024^2 capped 10%, E[cost]~4.68x, D-03); generate_samples parallel==serial byte-identical across -t 1 and -t 4 (D-12)
- [Phase ?]: [Phase 03-03]: DATA-02 sharded JLD2 cache — SHA-256 content hash over data-defining source bytes + canonical(config) names the cache dir and lives in meta.jld2 (D-05); atomic .tmp+integrity-check+mv shard writes with shard_done resume-by-skip (D-04/D-06); reserved >=20 ADVI holdout from disjoint XOR-salted holdout_rng into separate holdout.jld2 with negative global_index, structurally disjoint from main pool (D-10)
- [Phase ?]: [Phase 03-04]: DATA-03 leak-free loader — module Loader is the SOLE standardization path; ZScoreTransform fit on TRAIN columns only (D-07), no standardize_all symbol exists so leakage is impossible by construction (D-08), deterministic k=5 folds from fold_rng(master_seed XOR FOLD_SALT) (D-09), mask rows 65:128 bypass, holdout excluded structurally (D-10)
- [Phase 05-01]: Phase-5 anti-snooping contract locked in committed spike/validation/consts.jl BEFORE any reported run — M=2000/L=999/bins=50, all SBC/BF/OOD thresholds, VAL_MASTER_SEED=0x5BC0FFEE disjoint from NPE_MASTER_SEED=0xC0FFEE (via VAL_SALT); fixture gates use a separate VAL_FIX_SEED so they never consume the reported stream (D-02/SBC-03)
- [Phase 05-01]: One shared harness.jl (draw_simulate_infer + paired Δρ path) composes the frozen Phase-4 read surface over the loaded-once net; sbc.jl delivers 7θ+Δρ ranks, KS/χ² uniformity, coverage curve, and ported _bin_calibration ECE/MCE traffic-light; figures.jl owns ROC/BF stubs so Wave-2 plans call, never edit it (SBC-01..04)
- [Phase 05-02]: amortized log-BF = logratio(m=1)-logratio(m=0)-measured_log_prior_odds via NRE model-comparison (v0.2.1 has no l-POP; honest relabel)
- [Phase 05-02]: KDE compute_BayesFactor math ported verbatim to isolated spike/baseline/ env; shared spike/ env untouched (NeuralEstimators 0.2.1)
- [Phase 05-02]: fixture reproduction corr=0.961 (D-08a PASS) but max|dlogBF|=3.18 (D-08b NOT met at fixture scale) - reported honestly, consts not tuned; reported gate is 05-04
- [Phase 05]: [Phase 05-03]: two-channel OOD flag (train-only Mahalanobis density + posterior-predictive per-feature discrepancy, OR-fused); ID-quantile operating point committed on a HELD-OUT ID pool (out-of-sample honest ~5% FPR); hand-rolled ROC/AUC, no new package
- [Phase 05]: [Phase 05-03]: fixed-summary OOD power structurally bounded — TWO measured named blind spots: D-04 correlation-preserving transforms (affine exact, block-permute residual) AND detector-noise mismatch (density+PP AUC~0); optics/PSF robustly detected; full four-family ROC is 05-04 reported job
- [Phase ?]: [Phase 06-01]: spike/demo.jl two-tier seeded runner — fast tier re-exercises NPE+BF/NRE+OOD at fixture scale on VAL_FIX_SEED with per-layer twin-run bit-reproducibility asserts (global RNG seeded alongside Philox since sampleposterior threads no rng); step-0 git-status assertion certifies src/ untouched; --full spawns reported gates as independent subprocesses; locked files byte-unchanged
- [Phase ?]: [Phase 06-02]: Go/No-Go memo renders a Clean Go — literal pre-registered gates did NOT pass as written, Go justified by reading each residual failure to a non-method cause (chi2-over-power at M=2000, residual data-scale gap, clamped-KDE tail artifact); credibility backed by the D-05 independent fresh-seed re-pre-registered confirmation ship-gate inside Phase 7
- [Phase ?]: [Phase 06-02]: BF max|d logBF|=10.9535 quoted after runtime-confirming against demo.jl printed bf_report max_abs_err (A1 resolved); memo reports BOTH number sets — Set 1 pre-registered all-FAIL (0.8861/0.730/2.81e-40) verbatim from 05-04-SUMMARY.md, Set 2 post-hoc (0.0164/0.9358/AUC 1.0) with data-snooping exposure + 3 mitigations (consts.jl vs e9c91d3, DEV seed 0xDE7C0DE, single confirmatory VAL run)
- [Phase ?]: [Phase 06-02]: demo.jl machine-gates DEMO-03 — assert isfile(MEMO_PATH) + occursin content asserts (Clean Go, falsification, 0.8861, 0.9358, ship-gate, Phase 11, Phase 13); SC3/DEMO-03 table rows now report real PASS
- [Phase 07-01]: Grid coupling centralized ONCE in src/amortized/summary.jl — summary_dim(G)=2G^2 / cont_rows(G)=G^2 / ratio_input_dim(G)=5G^2 are the single source of truth (8x8 reproduces 128/64/320); patch_summary(mci,G) + encode_d01 + _summary_row_partition(:min,2G^2) grid-general for G in {4,8,16,32}. Every downstream amortized module derives dims from these helpers (no per-grid re-hardcoding)
- [Phase 07-01]: Estimator registry (PROD-02) keyed by patch grid — _SHIPPED_GRIDS=(4,8,16,32), 64 DROPPED (D-04); estimator_for validates grid and throws a train_and_register-pointing ArgumentError (T-7-05, no silent default-grid fallback); has_cuda_device() weakdep-safe (D-06 graceful CPU fallback); _lazy_load_from_artifact!/_train_grid_pipeline are honest hook-stubs for later plans
- [Phase 07-01]: Grid-parametric datagen (src/amortized/datagen.jl) — summary buffers Matrix(summary_dim(G),N), generating_config.summary_min_dim=summary_dim(G) so grids auto-separate into distinct content-hash cache dirs; imsize_set exposed for per-grid image-size bias (>=15-survivor floor, T-7-03); :min-only (spike :aug superset dropped); content hash uses Base.hash to avoid re-resolving the fragile Wave-0 Manifest; simulator chain referenced-but-promoted-later; spike/ byte-untouched
- [Phase 07-02]: Amortized READ surfaces promoted to src/amortized/{infer,bf,ood}.jl (PROD-01), grid-general + CPU-default (use_gpu=false everywhere, D-06). infer.jl: standardize_summary (frozen zt, mask bypass) + posterior_for/rho_draws/delta_rho (StatsBase.reconstruct BEFORE ρ read, Pitfall 5). bf.jl: amortized_log_bf (one NRE pass, measured log_prior_odds subtracted) + pair_encode (nc=G² derived ⇒ 5G²) + kde_log_bf_unclamped (Memo §5/T-7-06: compute_BayesFactor KDE math WITHOUT the 1e-8 _clampp floor so max|Δ logBF| is artifact-free). ood.jl: density (fit_ood_nulls continuous-rows+ridge) + noise (10 invariant features) + re-enabled posterior-predictive channel (ood_verdict defaults with_pp=true, Memo §5/T-7-04, kept crash-free by the iter1 _finite_or/_theta_tuple finite-guard) OR-fused by ood_verdict → OODVerdict. OOD ship-gate experiment (misspec families/ood_roc_over_grid, need simulator+ImageFiltering) deferred to 07-03+.
- [Phase 07-02]: Declared LinearAlgebra as a direct stdlib dep (cholesky/Symmetric/I for the OOD Mahalanobis) — a stdlib already in the Manifest with no version to resolve, so the co-resolution gate stays 4/4 green (NeuralEstimators 0.2.1 / Flux 0.16.10 pins intact); Pkg.test fully green (infer 11/11, bf 18/18, ood 31/31); spike/ byte-untouched.
- [Phase 07-04]: Per-grid CPU-reproducible ship-gate machinery (D-05) delivered to test/gate/{harness,sbc,gate_consts_template,run_gate}.jl + test/gpu_smoke.jl (PROD-02) — no grid trained yet; this is the gate the per-grid plans invoke. harness.jl draw_simulate_infer(m,rng;G,...) + paired-Δρ path is grid-parametrized (patch_summary(mci,G)), CPU-only (use_gpu=false on EVERY NeuralEstimators call), and loads a per-grid net via load_estimator (not the fixed spike trained_npe.jld2); frozen zt/θzt applied never re-fit (Pitfall 5). sbc.jl: M×8 rank table (7 θ + dedicated paired-draw Δρ column), KS+χ² uniformity via HypothesisTests (never hand-rolled), coverage curve, ported ECE/MCE CalibrationResult traffic-light, aggregated sbc_gate verdict. gate_consts_template.jl is the FRESH-per-grid pre-registration (SBC/BF/OOD consts) carrying a disjoint PROD_SEED[G] via a Philox stream salted (PROD_SALT) off the FORBIDDEN VAL_MASTER_SEED=0x5BC0FFEE and NPE_MASTER_SEED=0xC0FFEE — unit-asserted PROD_SEED[G]∉{those} for G∈{4,8,16,32} + 4 distinct seeds (anti-snooping Pitfall 3/T-7-08); BF gate uses the non-clamped kde_log_bf_unclamped baseline (Memo §5/T-7-06), OOD gate with_pp=true. run_gate.jl per-grid CLI --grid G [--sbc --bf --ood] loads the grid net + selected artifacts, runs the gates (use_gpu=false), writes an atomic .tmp→integrity→mv gate report; invocable before any grid is trained (missing NPE → :not_trained, no crash). bf_gate is self-contained (amortized_log_bf vs non-clamped KDE, no Turing); ood_gate computes the pre-registered id_threshold now + control-separability AUC when the per-grid plan injects pos_sim. gpu_smoke.jl: train_npe(...;use_gpu=true) runs on GPU when present and DEGRADES CLEANLY to CPU when CUDA absent (no error, CLAUDE.md graceful-fallback), asserts CPU-resident Flux.state persistence reloads/infers CPU-side; CUDA-present branch guarded by has_cuda_device() (==false here → fallback path exercised). Forward simulator NOT yet promoted (referenced-only in datagen/ood bodies) → harness takes an injectable `sim` seam (default_simulator resolves the promoted chain, errors clearly until then), mirroring 07-03's injectable datagen; the fixture SBC smoke injects a lightweight fake simulator. Pkg.test green (gate harness+SBC 29/29, GPU smoke 8/8, co-resolution 4/4 intact, no new external deps); spike/ byte-untouched.
- [Phase 07-07]: 16x16 (SHIP-WITH-CAVEAT fine grid) trained via _train_grid_pipeline(16) (80k pairs, 68k train, NPE d_in=512/dstar=64/10 coupling, NRE input_dim=1280) and CPU-gated on fresh disjoint PROD_SEED[16]=0xb906f369f6cacf91 with SBC_IMSIZE raised to 512^2. HONEST FAIL, constants byte-locked: SBC ks_pass=false/ece_pass=true with only 4/8 KS pass and the HEADLINE rho_true itself rejected (KS p=7.0e-4, ECE 0.0240) plus spillover 0.0027 -- both were comfortably uniform at 8x8 (0.567) and 4x4 (0.427), so the documented 'M=2000 over-sensitivity on summary-uninformative nuisance parameters' reading does NOT cover this; the paired Delta-rho column survives (KS 0.161, chi2 0.897, ECE 0.0060). BF corr 0.9150 (LOWEST of the three grids: 8x8 0.9472, 4x4 0.9410) AND max|d logBF| 12.72 (documented KDE-tail artifact), n=15 finite of 25. OOD NOT RUN / INCONCLUSIVE -- auc=nothing, passed=nothing (no pos_sim positive control injected); id_threshold=378.44 only, explicitly NOT a pass. CONSEQUENCE for 07-10: grid 16 recommended NOT eligible for default registry population (leave out of _SHIPPED_GRIDS or ship only behind explicit sign-off), and any entry must carry the >=512^2 minimum-image-size caveat (256 px/patch at 256^2 is marginal after background exclusion). PROVENANCE GAP: _train_grid_pipeline does not persist imsize_set, so the realized training image distribution is not recoverable from artifacts/grid_16/* alone.
- [Phase 07-08]: Local localisation ships as WINDOWED sub-tile inference over the frozen 8x8 bundle (src/amortized/local_map.jl), NOT a new fine-grid model — local_coloc_map returns a lightweight LocalColocMap (grid, tiles, delta_rho::Matrix, ood_flag::Matrix, meta), deliberately NOT an AbstractColocResult (no draws / no BF / no per-region uncertainty to back that interface). grid is a defaulted KEYWORD (8) so the map inherits whichever grid's gate the caller picks and CI can drive a tiny grid-4 bundle. Sub-tiles come from the package's OWN patch(img,nx,ny) (same trimming as the summary path) and tile MCIs carry the PARENT Otsu thresholds (tile-local Otsu would make tiles incomparable). Degenerate tiles take LOCAL_MAP_SENTINEL=0.0 AND a forced ood_flag=true so 'unscorable' can never read as 'measured no difference' (T-7-04). Wave-6 bridge _ensure_grid_registered(grid; artifact_dir) reuses _bundle_from_artifacts + register! rather than restating the three loads; it is a NO-OP on an already-registered grid.
- [Phase 07-03]: Amortized TRAINING layer promoted to src/amortized/{architecture,train_npe,train_ratio,persist,pipeline}.jl (PROD-01/02). build_estimator input-width-agnostic (q a NormalisingFlow INSTANCE positional; NPE_* Phase-5 consts); train_npe/train_ratio default use_gpu=has_cuda_device() with the spike use_gpu&&throw guards REMOVED (D-06) and LR/decay kept Float64 (AdamW-CosAnneal gotcha); ratio conditioner width from ratio_input_dim(G)=5G² (not 320), no custom loss (v0.2.1 hard-codes logit-BCE). SHIPPED PERSISTENCE migrated to CPU-resident Flux.state(cpu(est))+arch metadata → build_estimator/build_ratio_estimator + Flux.loadmodel! on load (Pitfall 4/T-7-07: device-independent, narrower deserialization surface T-7-01); OOD nulls persist directly; all through the atomic .tmp→reopen-integrity-@assert→mv(force=true) wrapper + schema_version + _estimator_ok/_ratio_ok/_ood_nulls_ok SKIP-IF-DONE predicates; seeded save→load→CPU inference bitwise-equal. _train_grid_pipeline(grid;...) factors datagen→zt→train_npe→train_ratio→fit_ood_nulls→persist→EstimatorBundle (SKIP-IF-DONE loads valid artifacts; injectable datagen seam makes it testable without the not-yet-promoted simulator; default_imsize_for 16→≥512²/32→≥1024²); train_and_register(grid) now runs it (D-04 on-ramp is real code, stub removed). Docstring documents -t auto datagen + Philox-per-index thread-count-independent byte-identical repro. GPU path plumbed but NOT exercised (has_cuda_device()==false, no CUDA loaded); CPU path fully green. Pkg.test green under -t auto (train_npe 11/11, train_ratio 12/12, persist 17/17, pipeline 19/19; co-resolution 4/4, no new external deps); spike/ byte-untouched.
- [Phase 07-06]: 4x4 bundle produced via the PUBLIC train_and_register(4) path (PROD-02/D-04 user-definable-grid happy path proven end-to-end; estimator_for(4) confirmed in-process), imsize_set constrained to ((256,256),) for the CPU budget; gate on fresh disjoint PROD_SEED[4]=0x8c0ad97b99bd6031 recorded honest FAIL (SBC ECE-green all 8, 7/8 KS pass vs 8x8 5/8, only shift_dx 0.0326; BF corr 0.941 near-miss + max|dlogBF| 12.40 KDE-tail artifact; OOD ID-op 31.79 only), gate_consts_4.jl byte-locked
- [Phase 07-09]: 32×32 CAPPED (user cap-the-family) — NOT trained, NOT shipped (gate-32x32.md). The GO rests on 8×8; a 32×32 gate would run the KDE-BF/non-atom-corrected apparatus spikes 006-014 proved defective. No training run, no grid_32 artifacts, no gate_consts_32.jl, no PROD_SEED[32] consumed. 64×64 stays dropped.
- [Phase 07-10]: Amortized-only public API shipped — colocalization_amortized(img,control,channels;num_patches=8,N=2000)→AmortizedColocResult (src/amortized/api.jl), composing registry→frozen summary→NPE posterior→Δρ cloud→NRE log-BF→density-channel OOD, CPU-default (D-06), D-02 accessors dispatching. _SHIPPED_GRIDS=(8,) ONLY (Go/No-Go Option A; 4 never post-hoc re-analysed, 16 gate-FAILED, 32 CAPPED, 64 dropped). Content-hashed lazy Artifacts.toml [grid_8] (git-tree-sha1=90e6b63a…, lazy, Release-hosted tarball sha256=17162904…) → _lazy_load_from_artifact!(8) resolves store→dev-fallback→download, tree-sha1 VERIFIED before load (T-7-01), then wires the recorded density-channel OOD op-point 179.14. Shipping bundle = amended_v2/grid_8 (retrained bounded-θ/realistic-imsize net the GO rests on). Artifacts+Pkg stdlibs added to [deps] (resolve clean, NeuralEstimators 0.2.1/Flux 0.16.10 intact). Turing path internal-only. test/test_integration.jl + docs/amortized.md (named limits §5 verbatim). Named limit: shipped OOD flag = density Mahalanobis channel only (fused noise/PP are gate-time constructs).
- [Phase 07-06]: fixed a Rule-1 numerical crash in the shipped non-clamped KDE BF baseline (bf.jl kde_log_bf_unclamped) — QuadGK adaptive integral overshoots the tail probability a few ulp outside [0,1], crashing log(); clamped to the valid [0,1] domain so a saturated tail yields the honest ±Inf (dropped) NOT the forbidden 1e-8 finite floor; in-range values byte-identical
- [Phase 08-05]: BOTH physical ground-truth anchors human-verified and pinned (SC1). POSITIVE = TetraSpeck 100 nm multicolor fiducial beads (RegiSTORM v1.0.0 sample data, Zenodo DOI 10.5281/zenodo.5509861, CC-BY-4.0, Karlsson et al. 2023 BMC Bioinformatics 10.1186/s12859-023-05320-1): one physical bead emits in EVERY colour channel ⇒ coloc BY CONSTRUCTION and state-independent (multi-channel STORM .tif frame stacks + ThunderSTORM .csv; a mean/max projection is owed in Phase 16). NEGATIVE = "Light My Cells" (France-BioImaging / ISBI 2024, BioImage Archive S-BIAD1047, CC-BY-4.0, OME-TIFF): nucleus (DNA) vs mitochondria are disjoint BY BIOLOGY from the acquisition design (fields must be FILTERED in Phase 16 to those carrying BOTH channels). Neither label is a computed coloc score (Pitfall 1) and neither is an environment-quenched tandem (Pitfall 2).
- [Phase 08-05]: D-01 DEVIATION, human-accepted and recorded (never papered over) — the literal decision specified a CELLULAR tandem-fluorophore construct, but no open-licensed, NON-environment-quenched tandem-FP dataset could be verified in any public archive; the deposited tandem-FP sets are quenched mRFP/mCherry-EGFP-LC3 autophagy reporters, which the acceptance bar disqualifies. Substituted a multicolor-bead dataset: the SAME physical particle in both channels is a strictly STRONGER physical positive (no pH/quenching failure mode at all), but it IS a deviation from the literal wording.
- [Phase 08-05]: D-02 DEVIATION — the matched-same-study preference is NOT satisfiable; no qualifying single-study deposit carrying both a same-particle positive and a segregated negative could be verified. The anchors are cross-study / cross-archive (`ANCHORS_MATCHED = false`). KNOWN LIMITATION: imaging-condition confounds (microscope, objective, exposure, detector, prep) between positive and negative are NOT controlled; Phase 16 must report this.
- [Phase 08-05]: Hash discipline (T-08-16) — `sha256` uses a deliberately NON-hex `"PENDING-FETCH"` sentinel rather than an empty string (an empty hash is indistinguishable from an unfilled CBS row and could be silently accepted downstream), with `is_real_sha256`/`is_pending_hash` making "real 64-hex digest OR sentinel, never in between" a testable invariant. The bootstrap fetch was DELIBERATELY not executed — the positive anchor is a ~6.3 GB archive, so the bulk download stays an explicit human decision. `bootstrap_anchor_hashes()` is implemented and must be run online+authorized BEFORE Phase 16. No digest was fabricated.
- [Phase 08-05]: `finalized_manifest()` (corpus/anchor_rows.jl) SUPERSEDES the 08-04 placeholder `committed_manifest()`; `anchor_rows()` is parameterized on sha256/bytes so the pending path and the post-bootstrap path are the SAME code path and both are covered by the offline gate. All anchors remain `physical-primary` + `sealed_holdout`, reachable ONLY via `open_sealed_holdout(; reason)` (D-09) — asserted.
- [Phase ?]: 14-09: SC1-d PASSES -- realized split-conformal coverage 0.9115 (1823/2000) vs the derived band 0.8799; q-hat 0.4498; iteration_trigger_fired=false so P14_ITERATION_ALLOWANCE stays unspent. Qualifier: ambiguous rate is exactly 0.0, so coverage comes from classifier accuracy not set-valued caution; the guarantee is simulator-derived (D-03a: the intended real-data bound is ABSENT, not merely loose).
- [Phase ?]: 14-10: SC1-b MET — realized FDP inside the binomial 95% band of the rule's own predicted value at all four levels of P14_ALPHA_FDR_GRID, decided fraction 0.888 (1776/2000) at every level. SC1-c REPORTED, not gated: the guarantee is genuinely prior-sensitive — at pi_coloc = 0.70 the realized FDP is 0.347 against a predicted 0.200, so an FDR number may never be quoted without the prior it assumed. Accepted null mass is ~all random, ~0 exclusion.
- [Phase ?]: 14-11: SC3-a/b/c all MET on the shared 2000-item pool (skill 0.9017 vs floor 0.60, spearman 0.9702 vs 0.95, AUC_hard 0.9117 vs 0.80); all four bars recorded in-artifact as judgement calls with no derivation; coverage floor frozen in a6c8258 before the runner existed; every number is a well-specified-regime number

### Roadmap Evolution

- v2.0 feature expansion appended to the current milestone (Phases 8–16), derived from `next_project_analysis/11_plan_proteincoloc_v2.md` + scope decision log. Consolidated as ONE v2.0 paper (not split); multiplex/copula deferred to v2.1; expert-concordance dropped.
- Phase 8 added: External Physical Ground-Truth Corpus (Depends on: nothing — parallelizable now)
- Phase 9 added: Cross-Method Comparator Harness (Depends on: nothing — parallelizable now)
- Phase 10 added: Manuscript Skeleton + Related-Work Positioning (Depends on: nothing — parallelizable now)
- Phase 11 added: Registration + Chromatic Uncertainty as Latent (Depends on: Phase 7)
- Phase 12 added: Spatial Colocalization Map (GP/CAR) (Depends on: Phases 7, 11 — serialized behind 11 due to shared PosteriorEstimator/training code; descope-to-v2.1 candidate)
- Phase 13 added: Three-Hypothesis Amortized Bayes Factor (Depends on: Phase 7)
- Phase 14 added: Decision + Abstention Layer (Depends on: Phases 11, 12, 13)
- Phase 15 added: Calibration Operating Envelope + CI Gate (Depends on: Phases 11, 12)
- Phase 16 added: External Validation + Manuscript Assembly (Depends on: Phases 8, 9, 10, 14, 15)
- Dependency analysis (analyze-dependencies): sole correction vs. draft was serializing Phase 12 behind Phase 11 (file overlap on the estimator); Waves — A {8,9,10} now, B {11 ∥ 13, then 12}, C {14,15}, D {16}.

### Pending Todos

[From .planning/todos/pending/ — ideas captured during sessions]

None yet.

### Blockers/Concerns

[Issues that affect future work]

- **[SPIKE-WIDE — INCLUDE-GUARD POISONING: MECHANICAL SWEEP DONE, 29 HITS, AND THE FIX IS COUPLED TO THE BLOCKED M1/M2 RULING, 2026-07-29]** Ran the sweep asked for after `fb76b84` (`:SBC_M`) and the `:P13_DEV_SEED` abort: for every `isdefined(..., :NAME) || include(...)` and `if !isdefined(..., :NAME)` under `spike/`, is `:NAME` declared as a `const` by more than one file? **220 guards over 118 files; 29 POISONED, in eight families.** This is a repo-wide latent defect, not a third instance.
  | Sentinel | Guards | Declared by |
  |---|---|---|
  | `:P13_DEV_SEED` | **9** | `p13/consts.jl`, `validation/p12_consts.jl` |
  | `:LAMBDA_MIN` | **9** | `p13/preconditions.jl`, `validation/p11_consts.jl` |
  | `:P11_DEV_SEED` | 4 | `p13/consts.jl`, `p11_consts.jl`, `p12_consts.jl` |
  | `:NPE_MASTER_SEED` | 1 | **SIX** files (`run_advi.jl`, `train_npe.jl`, `p13/consts.jl`, `test_npe.jl`, `p11_consts.jl`, `p12_consts.jl`) |
  | `:ABL_REL_MARGIN`, `:ABL_FOLD_CONSISTENCY`, `:BENCH_THREADS` | 3 | each also declared in `test_npe.jl` |
  | `:MU_PRIOR` | 1 | `simulator/calibration.jl`, `simulator/prior.jl` |
  | `:P11_IMSIZE_SET` | 1 | `p13/preconditions.jl`, `p11_consts.jl` |
  | `:P13_REPO_ROOT` | 1 | `p13/real_images.jl`, `test/test_p13_result.jl` |
  **SEEDS DOMINATE, AND THAT IS STRUCTURAL:** 23 of the 29 guard on a seed or a prior bound — precisely the names every other phase MUST mirror to assert stream disjointness. A seed is the worst possible guard sentinel in this project, and four separate phases independently reached for one. (`p13/consts.jl:93` is a false positive of the scan — a commented illustration of the idiom, not a live guard. The live owner-side wrapper is `:100`.)
  **MECHANISM VERIFIED EMPIRICALLY, BOTH DIRECTIONS:** defining `P13_DEV_SEED` first exactly as `p12_consts.jl:109` does, then running the guarded include, reproduces the abort — guard short-circuits, `P13_GATE_M` never defined. Reading the same frozen file through an isolated module under the SAME poisoning recovers it — `_P13C.P13_GATE_M = 4000`, 100+ constants restored, file untouched.
  **WHY "GUARD ON THE NAME YOU OWN" (`fb76b84`) CANNOT SIMPLY BE APPLIED HERE, AND WHY THAT COUPLES THIS TO THE CHANNEL-PAIR DECISION.** `fb76b84` only had to fix CALLERS, because `validation/consts.jl` already guarded its own body on `:SBC_M`, the name it owns. Here BOTH ends use the poisoned name, and the owner-side end is `spike/p13/consts.jl:100` — **inside the byte-locked pre-registration**. A caller-only fix is provably useless: the caller would correctly call `include`, and the body wrapper would still see `P13_DEV_SEED` defined by Phase 12 and skip the entire body. So applying the precedent REQUIRES editing `consts.jl`, which changes its sha256 and trips `h.consts_sha == p13_consts_sha()` in `run_three_way_gate.jl:353` (**the gate**), `run_alpha_series.jl:451` and `run_p13_realimage.jl:634` — **the exact M1/M2 question already open above.** The two tasks are therefore one decision, not two.
  **RECOMMENDATION, NOT APPLIED (the ruling is the user's):** under **M2** (consts.jl stays byte-unchanged) the fix is the ISOLATED-MODULE READ, which is not a departure from the precedent but the only form of it available when the owner file is frozen — and it is already this repo's house pattern for this exact collision, used three times: `module _P11C` (`p13/preconditions.jl:160`, whose comment states *"NO GATE IS RUN AND THE FILE IS NOT MODIFIED -- this is a read"*), `module _GC` (`p13/consts.jl`) and `module GateV2` (`p11_consts.jl`). Under **M1** the wrapper at `:100` is re-pointed to a Tier-1 name `p13/consts.jl` exclusively owns — **100 such constants exist**, e.g. `P13_STRATIFICATION`, `P13_GATE_STATISTIC`, `P13_ITERATION_ALLOWANCE` — in the SAME edit that lands the channel-pair amendment, and the isolated module is then unnecessary.
  **COST OF THE M2 FORM, STATED SO IT IS NOT A SURPRISE:** `spike/test/test_p13_consts.jl` asserts ~114 constants UNQUALIFIED, so an isolated read needs them re-bound into the test module. That is mechanical but not free, and it is why this was reported rather than hand-applied mid-flight.
  **NOT FIXED, DELIBERATELY:** `p12_consts.jl` is append-only Tier-1 and its mirror is CORRECT — asserting seed disjointness requires naming the seeds. The defect is on the guard side every time.

- **[Phase 13 — M1 RULED; the channel-pair amendment and the include-guard fix now ride in ONE edit; PLANNED, NOT YET EXECUTED, 2026-07-29]** **The user has ruled on the blocked `13-17-PLAN.md` Task-1 mechanism question: M1 — fix at the source, inside the same authorised amendment. M2 (the isolated-module read) is the tool for files that CANNOT be edited because no amendment authorises opening them.** The reasoning, recorded because a mechanism chosen for a reason is auditable and one chosen by default is not: **the byte-lock breaks either way** — `consts.jl`'s sha256 changes for the channel-pair correction regardless, and `h.consts_sha == p13_consts_sha()` has to be re-derived in `run_three_way_gate.jl:353` (**the gate**), `run_alpha_series.jl:451` and `run_p13_realimage.jl:634` regardless — **so the include-guard fix rides along at ZERO additional pre-registration cost. That coupling is precisely why the answer is M1 and not M2.**
  **`13-D15-AMENDMENT.md` IS NOW FROZEN AND CARRIES TWO SEPARATELY DISCLOSED CHANGES (§0.4), with two distinct justifications, so that "we opened the file for X and also did Y" is VISIBLE ON THE PAGE rather than inferable from a diff.** **CHANGE A** (§5) corrects a constant that named the **WRONG PHYSICAL OBJECT**: `test/runtests.jl:105` records the fixture channels as `["blue","green","red"]`, so `c1` is the DAPI/Hoechst counterstain, not a target protein. **CHANGE B** (§5B) corrects a sentinel that **another phase must legitimately mirror**: `p12_consts.jl:109` declares `P13_DEV_SEED` to assert seed disjointness, that mirror is CORRECT, and `p12_consts.jl` is NOT touched. **Neither is a threshold moved after proving inconvenient. No seed, bar, floor, band, ladder, tolerance or allowance changes in either**, and `P13_DEV_SEED = 0x0000_0000_0B13_DE71` stays byte-unchanged at `consts.jl:188`.
  **SENTINEL CHOSEN AND ARGUED: `P13_DECLARED_DEVIATIONS` (`consts.jl:662`).** Verified by fresh repo-wide grep to be declared by `spike/p13/consts.jl` alone (the only other file that even mentions it is `test_p13_consts.jl`); Tier-1 and declared UNCONDITIONALLY inside the guarded body `100-731`; **not a seed and not a prior bound**; and about the FILE'S IDENTITY rather than a value another phase might quote — a Phase-14 deviation register would be `P14_DECLARED_DEVIATIONS`. **Rejected on the page:** `P13_TAU` (Tier-2 and legitimately SOMETIMES ABSENT — `consts.jl:745` tolerates its absence, and a sometimes-absent sentinel is worse than the bug), `P13_GATE_STATISTIC` (a gate knob), `P13_CUT_VARIANT` (persisted into the net artifact, so a later phase has a motive to mirror it), `P13_LAMBDA_PLACEMENT` and `P13_INPUT_WIDTH_RULE` (architecture rules on the cross-phase compatibility surface), `P13_STRATIFICATION` (a design knob, weaker on the identity criterion).
  **BOTH ENDS PLUS SIX SIBLINGS, IN ONE EDIT.** A caller-only fix is provably useless here — the caller would correctly call `include` and the wrapper at `:100` would still see the Phase-12 mirror and skip the whole body — and this project's signature failure is a fix applied at one end of a pair. The nine sites are `consts.jl:93` (commented illustration), `consts.jl:100` (the wrapper), `test_p13_consts.jl:42` (the aborting caller) and `labels.jl:73`, `net.jl:89`, `preconditions.jl:191`, `result.jl:79`, `tau_probe.jl:92`, `toy_gaussian.jl:79`.
  **`consts_sha` RE-DERIVATION IS FIRST-CLASS PLAN WORK, NOT A SIDE EFFECT, and the guard ends up STRONGER.** A widened `h.consts_sha == p13_consts_sha() || h.consts_sha == PRE` is **REJECTED**: it would accept ANY future drift silently. Instead `spike/p13/net.jl` gains one named, dated `const P13_CONSTS_SHA = (pre_amendment = ..., post_amendment = ...)` (both literals copied from command output, never retyped, never copied out of a document) and each of the three sites asserts **BOTH sides** against it. The original only required the two sides to agree with each other — a retrain plus an edit would have passed it. The replacement fails on a net trained under any other pre-registration AND on `consts.jl` at any third sha, including a future undisclosed edit. **Retraining is NOT authorised; `gate_report.jld2` and `alpha_report.jld2` stay byte-unchanged; the allowance stays 1 of 1 UNSPENT.**
  **ORDERING: the guard fix and the sha re-derivation LAND AND COMMIT BEFORE the 13-16 re-run.** Until the guard is fixed nobody can use the suite to verify anything, including the corrected 13-16 — so fixing it first is what makes the re-run verifiable instead of another per-file-only result. The plan captures the BEFORE suite signature (Task 1) and the AFTER signature (Task 3) side by side. **The suite is still expected to exit 1** (paused Phase-4 gate + `test_p13_correction.jl`'s two committed MEASURED MISSES); what must change is that the Phase-13 block EXECUTES instead of being masked.
  **NEW FINDING, RECORDED BECAUSE THE GUARD FIX MAKES IT LIVE (no user decision needed; verified with a stop rule):** while the guard was skipping the body, ~104 `const` declarations never executed in a full-suite process. Repairing it makes them execute into `Main` alongside earlier declarations of the same names. **Ten names collide; three disagree in TYPE**, because in Julia a hex literal's width is its type: `NPE_MASTER_SEED`, `VAL_MASTER_SEED` and `VAL_FIX_SEED` are `UInt32` in `npe/train_npe.jl:65` / `test/test_npe.jl:81` / `baseline/run_advi.jl:74` / `validation/consts.jl:76,79` and `UInt64` in the three pre-registrations. The VALUES agree; only the widths differ, and all three are FORBIDDEN foreign seeds Phase 13 reserves rather than uses. Julia 1.12.6 permits constant redefinition so this is expected to be a non-event — **but it is an expectation that has never been executed**, so `13-17-PLAN.md` Task 3 verifies it. **STOP RULE: no seed literal may be edited to resolve it**; if `invalid redefinition of constant` fires, record it here and STOP.
  **THE 28 OTHER POISONED GUARDS: DOCUMENTED, NOT FIXED.** Fixing another phase's frozen file on a planner's own initiative is exactly the act this project's amendment discipline exists to prevent. **DEF-12-03** (`.planning/phases/12-spatial-colocalization-map/deferred-items.md`) is the right ledger for this collision and is accumulating correctly alongside `13-09`, DEF-12-01 (FIXED, `fb76b84`) and DEF-12-02; `13-17` closes it **BY REFERENCE** and does **not** edit Phase 12's file, because a live Phase-12 executor owns it and a cross-phase write for a one-line status change is not worth the merge.
  **THE STANDING RULE NOW HAS A DURABLE HOME A FUTURE PHASE WILL ACTUALLY READ: `.planning/CONVENTIONS.md`, registered from `.planning/PROJECT.md`'s `## Constraints` block** — chosen because every GSD plan's `<context>` references `@.planning/PROJECT.md`, so a Phase-14 planner reaches it without being told to. STATE.md was rejected as the home (a 900-line running log whose blocker entries get closed — a rule buried in it is a record, not a convention), as were Phase 12's `deferred-items.md` (Phase 14 will not open it) and any Phase-13 plan. `CLAUDE.md`'s `## Conventions` section is the natural second home and is RECOMMENDED TO THE USER rather than written by the planner. The rule, stated once for every future phase: **never guard an include on a seed or a prior bound, because those are exactly the names other phases are required to mirror** — 23 of the 29 hits do exactly that, and that one sentence would have prevented all 29. `.planning/CONVENTIONS.md` also records **C-02** (a hex literal's width is its type — the finding above) and **C-03** (`git diff` context lines are not changes: gate on `-U0` + `^[+-]`, or a negative criterion fails on CORRECT work and trains an executor to override the one check that proves nothing moved).
  **STATUS: planned only.** `13-17-PLAN.md` is rewritten to six tasks and validated; NO `.jl` file has been edited, NO experiment has been run, and the re-run has NOT happened. The re-run stays ONE run, reported whatever it says, and is required to answer IN WORDS whether the conclusion changes — **"corrected, conclusion unchanged" is pre-authorised as a GOOD result and rescue framing is forbidden.**

- [Phase 1]: NeuralEstimators v0.2.x is pre-1.0 and the Julia ML ecosystem is mid Flux→Lux/Reactant migration — the CPU-only Flux path must be confirmed exercised at the smoke gate; pinned Manifest is the mitigation
- [Phase 5]: OOD detection power is structurally bounded by the fixed patch-correlation summary — summary-orthogonal misspecifications are provably undetectable and must be named, not hidden
- [Phase 5]: "Tune until calibrated" is a data-snooping hazard — M and SBC threshold must be pre-registered before Phase 5 planning runs
- [Phase 5-04]: ALL THREE reported pre-registered gates FAIL at locked consts on fresh VAL_MASTER_SEED (honest negatives, consts.jl untouched): SBC every-parameter KS/chi2 p~0 (overconfident/miscalibrated NPE); BF corr=0.886<0.95 AND max|dlogBF|=3.60>0.5; OOD density pooled AUC=0.730<0.80 (noise family AUC=0 blind spot). Feeds Phase-6 Go/No-Go as No-Go-or-iterate.
- [Phase 5-iter1 (2026-07-02)]: BOUNDED calibration iteration retrained the NPE on the SAME 50k cache with a higher-capacity flow (dstar 32→64, summary 2×128→3×256, coupling 6→10) + stabler recipe (LR 5e-4→2.5e-4, batch 64→128, 300 epochs/patience 40). Selected on a DISJOINT dev seed (0xDE7C0DE); consts.jl BYTE-UNCHANGED. Calibration improved dramatically — DEV SBC ρ_true KS D 0.158→0.038, Δρ 0.155→0.017; the frozen retrain overwrites trained_npe.jld2. CONFIRMATORY (VAL_MASTER_SEED, run once): SBC still FAIL overall but ALL 8 params now GREEN on ECE (ρ_true ECE=0.016) and 3/8 (spillover/autofluorescence/noise) pass KS — the strict all-8 KS+50bin-χ² conjunction rejects (label_efficiency/shift genuinely non-uniform, χ² hyper-sensitive at M=2000). BF FAIL corr=0.918 (↑ from 0.886) but max|dlogBF|=14.4 (↑ worse — sharper NPE Δρ tails inflate the KDE-baseline gap). OOD FAIL unchanged density pooled AUC=0.730 (Mahalanobis channel is NPE-independent; noise-family blind spot structural). Verdict: PARTIAL. Next levers: scale data toward 200k (SBC), retrain NRE (BF).
- [Phase 5-iter1]: OOD PP channel now VIABLE — finite-guard added in ood.jl (_theta_tuple maps non-finite posterior θ̂ to in-range fallbacks) so strong-misspec re-simulation no longer crashes; PP θ̂ finite on all 4 families. Still excluded from the reported OR-fusion (with_pp=false unchanged in run_ood.jl); available to re-enable in Phase-6/7.
- [Phase 5-iter2 (2026-07-02)]: TASK 1 — retrained the NRE on the SAME 50k (θ,summary) cache the NPE used (assemble_ratio_data_cache; 16× the prior 3k fresh-sim pairs, leak-free: cache gen-seed 0x134d8f3 disjoint from VAL/NPE seeds), higher capacity (num_summaries 32→64, conditioner 64→3×256), stabler recipe (LR 2.5e-4, batch 128, 300ep/pat40), + research-A7 difference encoding (pair_encode = concat + 64-dim continuous-row correlation CONTRAST, input dim 256→320; shared with bf.jl). Develop DEV 0xDE7C0DE (corr 0.862→0.903); CONFIRMATORY VAL_MASTER_SEED once (run_bf.jl unchanged): corr 0.918→0.936 (↑), max|Δ logBF| 14.4→10.95 (↓). BF gate still FAILS. VERDICT on max|Δ|: it is a KDE-BASELINE TAIL/CLAMP ARTIFACT, not an NRE deficiency — every large-|Δ| sweep point (|Δρ|≳0.4) is where the sharpened NPE Δρ posterior is fully one-sided, forcing the KDE baseline's P(Δρ>0) to its _clampp=1e-8 floor → logBF=±18.41, while the amortized NRE (correctly) stays bounded ~±8; mid-range (|Δρ|<0.4) agrees within ~1-2. D-08b (max|Δ|≤0.5) is structurally unclearable against a clamped-KDE baseline at the sweep tails — a pre-registration nuance, consts.jl untouched. BF fast test 17/17.
- [Phase 5-iter2]: TASK 2 — added an AUXILIARY image-noise OOD channel (ood.jl Channel 3: noise_features/fit_noise_null/noise_score) OR-fused with the density channel, closing the detector-noise blind spot. 10 scale/rotation/permutation-INVARIANT features (HF energy ratio, robust HF scale ratio, outlier fraction, HF excess kurtosis, lag-1 autocorr; ∇² via finite difference, NO FFT/new dep) designed from the SMOOTH training distribution (SC5, not the held-out test set). Null + robust-z fit on a fresh TRAIN-ONLY ID pool; reported detector = per-sample max of density/noise robust-z (continuous OR, D-05). maha_auc/combined_auc now report the FUSED detector; density maha_thr unchanged so run_ood's density-only neg-control flag is byte-identical. Develop DEV (noise-family AUC 0.0→1.0); CONFIRMATORY VAL once (run_ood.jl unchanged): OOD gate PASSES — all 4 families best-AUC=1.0 (noise 1.0 at every level), combined pooled AUC=1.0, ID fire-rate 0.05, neg-controls KS-invariant + density-quiet. OOD fast test 27/27. Updated blind-spot framing: affine+rotate stay quiet on ALL channels (true blind spots); block-permute stays density-quiet (frozen gate passes) but the noise channel correctly FIRES on its tile-seam HF artifacts (~0.93) — so block is no longer summary-orthogonal to the fused detector. Remaining blind spot is NARROWER: orthogonal to BOTH the 8×8 correlation summary AND the image-noise features (e.g. a pure positive-affine rescale, still exactly invariant). consts.jl BYTE-UNCHANGED (verified vs e9c91d3).
- [Phase 7 — F5 COVARIATE SHIFT (HIGHEST SEVERITY), 2026-07-21]: **The SBC calibration proof holds at 256² and is NOT demonstrated to transfer to the real 1376×1028 paper regime that Phase-2 decision D-08 anchors the design on.** Training runs 07-05/07-06 constrained `imsize_set` to `((256,256),)` as a documented compute-budget deviation, and the gates evaluate at `SBC_IMSIZE` = (256,256) for grids 4/8 and (512,512) for grid 16 — train and gate are self-consistent, so the recorded verdicts are internally valid, but they are valid *at 256²*. MEASURED on the grid-8 net: evaluating on the mixed imsize set instead of 256² moves `label_efficiency` mean(u) from 0.546 → **0.807** (z ≈ 16.8 at M=250), bias −0.130. At FIXED θ, changing image size shifts the mean per-patch correlation by roughly the ENTIRE `label_efficiency` prior range (G=8: 0.513 @2048² vs 0.470 @256² at le=0.95, Δ≈0.043, while le 0.65→0.95 moves it only 0.045) — the confound is as large as the signal. **Consequence: the chain "SBC passes ⇒ posteriors on real images are calibrated" does not currently hold, and the Phase-8 physical anchors and the Phase-16 blind evaluation both depend on exactly that chain.** Closing it requires training/gating at the real image dimensions, or demonstrating and gating summary invariance to image size — both need retraining and are OUT of the current scope. Full evidence: `.planning/phases/07-productionization-conditional-on-go/07-CALIBRATION-FINDINGS.md` (F5). Related F6: `_train_grid_pipeline` never persisted `imsize_set`, so grid 16's actual training image-size distribution is UNRECORDED; a zt forensic mildly favours 512²-only (0.4249 vs 0.4288 for a 512²-only pool, 0.405 for the mixture) but is explicitly NOT CONCLUSIVE at n=100 and the value is recorded as `unknown`, never guessed. Provenance is now persisted going forward (`training_imsize_provenance`); existing bundles honestly report `recorded=false`.
- [Phase 7 — F3 VACUOUS SBC PASS, 2026-07-21]: **grid 4's `label_efficiency` SBC pass is VACUOUS and must NOT be read as calibration evidence anywhere** (memo, manuscript, or any per-grid comparison table). MEASURED: post_sd/prior_sd = **0.994**, shrinkage centre = **0.8034** = the prior mean (0.80), tercile mean(u) = 0.18/0.49/0.82 — the textbook "posterior = prior" signature. Grid 4 passes SBC on this parameter (KS p = 0.567) **because it learns nothing about it**: a posterior that reproduces the prior yields uniform ranks BY CONSTRUCTION. Rank uniformity is NECESSARY BUT NOT SUFFICIENT for calibration, and a vacuous pass and an informative pass must never be tabulated as the same result. Related F1: where the SBC *does* fail on `label_efficiency` (grid 8 z=+7.4 p=1.84e-9; grid 16 z=+8.6 p=7.15e-13) the failure is a LOCATION BIAS, not overconfidence — tail mass 0.199/0.231/0.215 vs expected 0.20, histograms flat not U-shaped — so the standard overconfidence levers (widening, tempering, more capacity) are the WRONG remedy; the NPE systematically UNDERESTIMATES the parameter. F2 mechanism (INFERRED): the unbounded `NormalisingFlow` cannot represent the `Uniform(0.6,1.0)` box and the forward-KL objective constrains nothing about the implied marginal, so information gained, centre displacement and rank bias all grow monotonically 4→8→16 — **the bias grows with learning, not with ignorance**. MITIGATION SHIPPED (reporting only, no threshold or verdict changed): `sbc_gate` now reports per-parameter `shrinkage = post_sd/prior_sd` + a `vacuous` flag (cutoff `SBC_VACUOUS_SHRINKAGE = 0.95` in `test/gate/sbc.jl`, deliberately NOT in the frozen `gate_consts_<G>.jl`).
- [Phase 8-05 — RESOLVED 2026-07-21]: The human-verify blocker is CLEARED. The human confirmed both anchors and 08-05 Task 3 executed (commit 8eb9147): `corpus/anchor_rows.jl` + `corpus/test/test_anchors.jl` created, `corpus/manifest.csv` finalized (30 CBS + 2 sealed physical anchors), offline gate 217/217, `git ls-files corpus/data` empty, `src/` untouched. **Phase 8 is COMPLETE.** THREE residual items carried to Phase 16, all recorded honestly in `ANCHOR_PROVENANCE_NOTES`: (1) **sha256 = "PENDING-FETCH" on BOTH anchors** — the bootstrap fetch was DELIBERATELY not run because the positive anchor is a ~6.3 GB archive and a bulk download must stay an explicit human decision; no digest was fabricated (T-08-16). `bootstrap_anchor_hashes()` is implemented and must be run when online + authorized, and the sentinel replaced BEFORE Phase 16 opens the sealed holdout. (2) **D-01 substitution (human-accepted)** — no open-licensed, non-environment-quenched tandem-fluorophore dataset exists in any public archive (every deposited tandem-FP set is a quenched mRFP/mCherry-EGFP-LC3 autophagy reporter, disqualified by Pitfall 2). POSITIVE = TetraSpeck 100 nm multicolor beads (RegiSTORM sample data, Zenodo 10.5281/zenodo.5509861, CC-BY-4.0): the same physical particle emits in both channels, so coloc is by construction AND state-independent. (3) **D-02 preference unmet** — NEGATIVE = "Light My Cells" (BioImage Archive S-BIAD1047, CC-BY-4.0, nucleus vs mitochondria). Cross-study / cross-archive (`ANCHORS_MATCHED=false`); imaging-condition confounds between the two anchors are NOT controlled and Phase 16 MUST report this limitation. Phase-16 conversion work still owed: mean/max frame projection for the positive STORM stacks, and field-filtering the negative to fields carrying BOTH the nucleus and mitochondria channels.

- [Phase 11 wave-3 post-merge gate, 2026-07-25]: **`SPEEDUP_GATE` (NPE-03) FAILS — median speedup 92.50× and 83.97× on two consecutive serial runs against the pre-registered `> 100.0×` bar** (`spike/test/test_npe.jl:230`, gate at `:74`). Reproducible, not noise. It is the ONLY remaining red in `spike/test/runtests.jl`; it was invisible until commit `b298c00` because a stale θ-arity assertion in `test_simulator.jl` aborted the suite four includes earlier (Phase 3 is now 70/70, Phase 5 SBC/BF/OOD and Phase 9 comparator all green). **Likely cause, NOT yet confirmed by measurement:** Phase-11 D-09 appended `chromatic_eps`, so the spike NPE now trains an **8-marginal** flow instead of 7 and the wider flow is a slower forward pass; confirming needs an A/B of the benchmark at `D_flow` 7 vs 8, which has not been run. **Scope: this does NOT invalidate the shipped >100× claim** — the shipped `amended_v2/grid_8` bundle still carries a 7-marginal flow (`NPE_D = 7`, unmoved per D-16) and its `Artifacts.toml` pin is byte-unchanged; the failing number comes from a spike-side net that is no longer the shipped one. **`SPEEDUP_GATE` was deliberately NOT lowered** — relaxing a pre-registered threshold after seeing it fail is the "amended twice, credibility spent" pattern of `07-GATE-AMENDMENT.md` §6.4. Two honest options, both user decisions, neither taken during execution: (a) confirm the 7-vs-8 hypothesis and record it as a named limit scoped to the research net, or (b) retire the legacy Phase-4 speedup gate explicitly as not meaningful against an 8-column prior. Must be carried into the Phase-11 report (plan 11-11). Full detail: `.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/deferred-items.md`.

- [Phase 4 `SPEEDUP_GATE` — 7-vs-8 A/B MEASURED, 2026-07-29, commit `ca994b9`]: **The "8-marginal flow" hypothesis recorded in the blocker above is FALSIFIED, on two independent grounds. No threshold, test or gate was touched.** **(1) The premise does not hold.** `speedup_report` loads `DEFAULT_NPE_MODEL = spike/npe/trained_npe.jld2`, whose persisted flow is **`q.d = 7`, `d_in = 128`**. That file last changed in `200e971` (the Phase-5 retrain); Phase 11 never touched it. The 8-marginal net (`q.d = 8`, `d_in = 129`) exists only in `p11_research_npe.jld2`, which the benchmark never loads. The failing gate has been timing a **7-marginal** flow all along. **(2) The width effect was measured anyway, and is ~16x too small to matter.** Interleaved, order-rotated, 9 reps at the gate's own `bench_N = 50` / `bench_seconds = 0.4`, CPU-only at 1 thread: `real7` (the actual gate model) median **97.04** [82.77 .. 110.09]; `fresh7` = `build_estimator(128, 7)` median **96.05** [87.18 .. 111.19]; `fresh8` = `build_estimator(128, 8)` median **94.27** [82.42 .. 109.13]. The 8th marginal costs **1.67%** of forward-pass latency (5.208 -> 5.295 ms), i.e. **1.78 speedup points**, against a **28.8-point** run-to-run spread. The arms' ranges OVERLAP and their per-rep ordering swaps repeatedly, so the width effect is not resolvable above machine noise and cannot move a 100x bar. **What the data DOES show:** the benchmark sits **ON** the bar rather than far below it — every arm, D=8 included, clears 100 in its best repetition, and the gate model needs 4.962 ms to clear 100 while measuring 5.178 ms (**1.04x**). Pass/fail is being decided by concurrent machine load (3 julia processes at start and end, observed varying 2-4 mid-run as peer agents started and finished; the sub-88 reps are the 4-process ones). That also explains the recorded drift 92.50 -> 83.97 -> 84.44 -> 89.77 -> 68.35 on this **unchanged** model, and answers the third question left open in `.planning/phases/13-three-hypothesis-amortized-bayes-factor/deferred-items.md`. **SCOPE CORRECTION THE USER SHOULD SEE:** the existing scope note states the failing number comes from "a spike-side net that is no longer the shipped one". On the width axis that is **not** the case — `src/amortized/architecture.jl:169-175` (`NPE_D = 7`, dstar 64, depth 3, width 256, coupling 10, flow 2/128) plus `summary_dim(8) = 2*8^2 = 128` make the shipped `grid_8` **topology identical** to the benchmarked one, and `trained_npe.jld2`'s stored metadata matches it field for field. The shortfall is therefore **not** insulated from the shipped >100x claim by the width argument. (Weights were NOT compared: the shipped `grid_8` artifact is a lazy remote pin and was deliberately not fetched.) **THE DECISION REMAINS THE USER'S and is NOT taken here.** Options (a) and (b) in the blocker above both rested on the now-falsified 8-column premise, so a third option is on the table: treat >100x as a claim whose measurement needs a **quiet-machine, pre-registered protocol** — the bar was derived on an unloaded machine and is being scored on a shared one. Reproduce with `julia --project=spike --threads=1 spike/npe/flow_width_ab.jl`; artifact `spike/npe/flow_width_ab.jld2` (carries per-rep speedups, latencies, and the julia-process counts).

- [Phase 4 `SPEEDUP_GATE` — **PAUSED, NOT RETIRED** (user ruling), 2026-07-29, commit `a494247`]: **The user has ruled: pause the gate's execution, do not retire it.** Rationale: a wall-clock benchmark only makes sense once development is finished, and must then be measured properly; until then the measurement is deferred. **PAUSED vs RETIRED is the entire point and is deliberate** — retiring spends pre-registration credibility this project cannot afford to spend a third time (`07-GATE-AMENDMENT.md` §6.4, "amended twice, credibility spent"); pausing spends none, because the threshold stands as the record and only the measurement is deferred. **What changed:** in `spike/test/test_npe.jl` SC3, the single assertion `sr.median_speedup > SPEEDUP_GATE` is now `@test_skip` rather than `@test`, tagged with the greppable marker **`DEFERRED-NPE-03-WALLCLOCK`** (grep that marker rather than a line number — the explanatory block shifts line numbers below it). **`SPEEDUP_GATE = 100.0` is BYTE-UNCHANGED at `:74`**, the assertion expression is unweakened and still present in source, and nothing was commented out or deleted — re-arming is a one-word edit (`@test_skip` -> `@test`). The RMSE half of the D-08 joint condition is deliberately **left live** (it is an accuracy claim, not wall-clock dependent, so a real accuracy regression still fails the suite). **The Phase-11 "8-marginal flow" hypothesis in the blocker above is FALSIFIED AT THE ROOT** and the entry is retained unedited because how long it was believed is itself part of the record: the benchmark never loaded that net at all (`trained_npe.jld2` is `q.d = 7`, `d_in = 128`, last changed in `200e971`, untouched by Phase 11; the 8-marginal net lives only in `p11_research_npe.jld2`). **NEW FINDING — only one side of the ratio is benchmarked.** `BenchmarkTools` IS used, but only on the denominator (`_time_npe_pair` -> `@belapsed`, minimum of many samples). The numerator `t_advi` is a **single un-replicated `time_ns()`** per pair (`spike/baseline/run_advi.jl:136-137`) frozen into `advi_artifact.jld2` on 2026-07-01 under that machine's then-current load, and never re-measured. A ratio of a benchmarked denominator to a one-shot numerator cannot be tightened by re-running the testset — which is a second, independent reason the number is unreliable, on top of concurrent-load sensitivity. **DEFERRED OBLIGATION (must be re-run, not dropped):** at the END of v2.0 development, on a QUIET machine (no concurrent julia processes), with proper benchmarking tooling on **both** sides of the ratio (re-measure `t_advi` under `BenchmarkTools` instead of reusing the frozen one-shot), re-arm the assertion against the **unchanged** bar of 100.0 and report the outcome. A failure under those conditions would be a real result about the >100x claim and must be reported as one. **SUITE EFFECT (the reason this was worth doing) — CONFIRMED BY RUN:** with the gate paused, `julia --project=spike --threads=1 spike/test/runtests.jl` no longer aborts at Phase 4, and **the previously-masked phases now EXECUTE and are GREEN**: Phase 4 **149 pass / 1 broken** (the broken is exactly the one `@test_skip`; SC3 itself is 8 pass / 1 broken, so the live joint-RMSE half still passes), Phase 5 SBC **28/28**, Phase 5 BF **17/17**, Phase 5 OOD **27/27**, Phase 9 comparator **45/45**, and the Phase-12 block runs in full (Tier-1 pre-registration 143/143, CNN estimator surface 1141/1141, decoupling 24/24, all 10 `P12-SUITE-RAN` markers printed). **The suite still exits 1**, but for a DIFFERENT and newly-unmasked reason recorded in its own entry below (`test_p13_consts.jl`, P13_DEV_SEED include-guard collision) — that abort is at `runtests.jl:220` and is NOT caused by this pause, NOT caused by Phase-4 code, and was simply invisible behind the Phase-4 abort until now.

- **[Phase 13 — `test_p13_consts.jl` ABORTS THE SUITE: `P13_DEV_SEED` include-guard collision, NEWLY UNMASKED, OPEN, 2026-07-29]**: **A pre-existing defect that was invisible behind the Phase-4 `SPEEDUP_GATE` abort and surfaced the moment that gate was paused (`a494247`). It is NOT caused by the pause, and NOT a Phase-4 problem.** In a full-suite run `test_p13_consts.jl` reports **22 pass / 1 fail / 114 error**, the errors all `UndefVarError` on Phase-13 constants (`P13_DECLARED_DEVIATIONS`, `P13_TAU`, …), and the thrown testset aborts `runtests.jl:220` — which masks the remaining **eight** Phase-13 includes (`test_p13_labels` … `test_p13_correction`, lines 221-239). **MECHANISM, fully traced:** (1) `runtests.jl:174` includes `test_p12_suite.jl` FIRST — deliberately, so Phase 12 cannot be masked by anyone else's abort; (2) that chain reaches `test_p12_consts.jl:42`, which **unconditionally** does `include(validation/p12_consts.jl)`; (3) `spike/validation/p12_consts.jl:109` defines `const P13_DEV_SEED = 0x0000_0000_0B13_DE71` — a deliberate **MIRROR** of Phase 13's seed, held so Phase 12 can assert stream disjointness; (4) at `runtests.jl:220`, `test_p13_consts.jl:42`'s guard `isdefined(@__MODULE__, :P13_DEV_SEED) || include(".../p13/consts.jl")` sees that mirrored symbol already present in the shared `Main` and therefore **SKIPS the include entirely**; (5) `spike/p13/consts.jl` wraps its whole body in `if !isdefined(@__MODULE__, :P13_DEV_SEED)` (lines 100-731), so **none** of the ~114 Tier-1 constants are ever defined. The guard symbol is thus a name Phase 12 legitimately owns a copy of, which makes it useless as a "have I loaded Phase 13's consts?" sentinel. **THIS EXACT BUG CLASS IS ALREADY DOCUMENTED IN-REPO, IN THE MIRROR DIRECTION:** `spike/p13/preconditions.jl:121-135` records that `p13/consts.jl` reserves the name `P11_DEV_SEED`, that this breaks a naive guard "in the suite while the standalone run passes", and repairs it locally in that one file. The same hazard has now recurred with `P13_DEV_SEED`, so the local repair did not generalise. **WHY IT WAS INVISIBLE:** per-file runs pass (nothing has defined the mirror yet, so the include fires normally), and the practice STATE.md prescribes for Phase 13 is exactly per-file runs — so both the per-file signal and the previously-aborted suite hid it. **NOT FIXED HERE, DELIBERATELY:** it sits in Phase-12 and Phase-13 files that a live Phase-13 executor owns, the fix is a design choice between at least three options (guard on a symbol Phase 12 does NOT mirror, e.g. `P13_DECLARED_DEVIATIONS`; load Phase-13 consts into their own module instead of `Main`; or have `p12_consts.jl` stop shadowing foreign seed names), and choosing among them is not this measurement task's call. **Consequence for verification:** a full-suite green is still NOT available as a Phase-13 signal, now for this reason rather than the Phase-4 one; per-file runs remain the working gate for Phase 13.

- **[Phase 13 — WRONG CHANNEL PAIR, OPEN, 2026-07-25]: `13-01-PLAN.md` and `13-15-PLAN.md` pre-register the real-image arm on the DAPI counterstain and MUST be corrected before Phase 13 executes.** `P13_REAL_CHANNEL_PAIR = (1, 2)` and `P13_REAL_REDUNDANCY_PAIR = (1, 3)`, justified as "the exact pair the frozen `ghat` calibration anchor was measured on" — but that anchor was itself measured on the wrong pair and has since been corrected (quick tasks 260725-vl8, 260725-wb7). `test/runtests.jl:105` records these fixtures as `["blue", "green", "red"]`, so **`c1` is the DAPI/Hoechst nuclear counterstain, not a target protein**; the colocalization pair is `c2`/`c3`. MEASURED induced μ — masked: `(2,3)` positive **+0.8238** (n=22) / negative **−0.0338** (n=20) vs `(1,2)` positive −0.0729 / negative +0.1149; unmasked: `(2,3)` **0.4603 / 0.3815** vs `(1,2)` 0.3292 / 0.2481. **Phase 13 needs the UNMASKED `(2,3)` values (0.4603 / 0.3815)** because it feeds the net, which was trained on unmasked summaries (`patch_summary` applies no Otsu mask; `src/amortized/summary.jl:54` matches) — masking the real images would put them out of distribution. Affected sites: `13-01` lines ~321, 364, 397, 442, 558, 580, 587 and `13-15` lines ~17, 119, 139, including a `must_have` and several verify criteria that make reproducing the WRONG numbers a pass condition. `P13_REAL_REDUNDANCY_PAIR = (1, 3)` must be DROPPED (DAPI-vs-red, −0.0925 masked on the positive fixture — noise, not a redundancy arm); with `c1` excluded and only three channels, no second pair exists. **Scope note:** hand-patching was started and deliberately reverted — this spans `must_haves` and verify criteria and needs a planner revision with the plan-checker, not spot edits. **Phase 11's `11-10-PLAN.md` had the same defect and IS fixed** (this commit). Honest consequence to carry into the Phase-13 report: unmasked, the two fixtures barely separate (0.4603 vs 0.3815), so the real-image arm is weak evidence — consistent with D-15 already scoping it QUALITATIVE ONLY.

- **[Phase 13 — WRONG CHANNEL PAIR: THE PRECONDITION WAS MISSED. 13-16 HAS ALREADY RUN. BLOCKED ON A NEW USER DECISION, 2026-07-29]** The blocker above says the correction "MUST be corrected before Phase 13 executes". **It was not, and Phase 13 executed to completion.** `spike/p13/consts.jl:552-553` still read `P13_REAL_CHANNEL_PAIR = (1, 2)` / `P13_REAL_REDUNDANCY_PAIR = (1, 3)` at the moment `2112bed` ran, and `spike/p13/run_p13_realimage.jl` consumes `P13_REAL_CHANNEL_PAIR` at ten sites (`:724`, `:734-735`, `:797`, `:830-852`). **So the real-image arm as shipped scores the DAPI/Hoechst nuclear counterstain against green, not the two target proteins.** Cause: the orchestrator dispatched 13-16 without reading this Blockers section, and the user's authorisation to correct arrived after 13-16 and 13-14 had both completed. Not an executor defect — **13-15 and 13-16 BOTH detected the problem independently and correctly declined to fix it**: `13-16-SUMMARY.md:276-279` names `c1` as the counterstain and records that "choosing a new channel pair" after the fixtures had been measured would be the wrong act, so it named the limit rather than editing a frozen constant. That was the right call under the pre-registration as it stood.
  **THE USER'S AUTHORISATION (received 2026-07-29, NOT YET APPLIED — nothing has been hand-patched):** correct to `P13_REAL_CHANNEL_PAIR = (2, 3)`; DROP `P13_REAL_REDUNDANCY_PAIR` entirely (with `c1` excluded and only three channels there is no second pair); use the **UNMASKED** values **0.4603** (positive) / **0.3815** (negative), because the net was trained on unmasked summaries and masking would push the real images out of distribution — the masked `+0.8238 / −0.0338` figures are the trap and must not be used however much better they separate; and document the whole thing openly as a **pre-registration amendment**, stating explicitly that correcting a factual mis-designation of the physical object is a *different act* from relaxing a bar that proved inconvenient. The honest consequence stands and must not be softened: unmasked, the fixtures barely separate (0.4603 vs 0.3815), so the corrected arm is **weak** evidence. Correcting the pair makes it honest, not strong.
  **WHY THIS IS STOPPED RATHER THAN DONE — three decisions the authorisation does not cover, because it was written on the assumption the correction would land BEFORE 13-16 ran:**

  1. **Does 13-16 RE-RUN, or are its numbers superseded in place?** The authorisation asks this only about 13-15. 13-16 is the plan that actually consumed the wrong constant.
  2. **Does Phase 13 REOPEN?** It is marked ALL 16 PLANS COMPLETE, `13-VERIFICATION.md` is written (`gaps_found`, no blocking gaps), and ROADMAP is ticked across all sixteen plans (`48fb808`, `a0efd7e`).
  3. **The real arm's HEADLINE NUMBER CHANGES.** The reported binding limit is the OOD flag at density 417.2974 vs ID threshold 167.5446 (**2.491×**) — measured on `(1,2)`. `13-16-SUMMARY.md:143` already records that a *different* channel pair lands at 607.57 (**3.63×**). So `13-REPORT.md`'s central real-image claim is pair-dependent and would need revision, not just a constant swap.
  **WHAT IS NOT AFFECTED:** 13-12's amended gate (6/6, AUC 0.990169 / 0.988262, ECE 0.0119808 / 0.012886) and 13-13's α-ladder are **SIMULATED** arms that never touch `P13_REAL_CHANNEL_PAIR`. The phase's gating verdict is untouched by this. Only the qualitative real-image arm (13-15, 13-16) and the sections of `13-REPORT.md` that quote it are in scope.
  **CONSTRAINT ON WHOEVER PICKS THIS UP:** STATE.md is explicit that a hand-patch was started and deliberately reverted — the change spans `must_have`s and verify criteria, several of which currently make **reproducing the wrong numbers a pass condition**. It needs a planner revision with the plan-checker, not spot edits. `consts.jl` remains BYTE-UNCHANGED (blob `70fe66df`, two commits in its whole history, `c42cc8e` and `bfba6ac`, neither postdating a result) and `P13_ITERATION_ALLOWANCE` remains 1 of 1 UNSPENT.

- **[Phase 13 — CHANNEL-PAIR AMENDMENT WRITTEN AND PLANNED; ONE DECISION ESCALATED, 2026-07-29]**
  The two channel-pair blockers above are now answered by a planning artifact rather than a
  hand-patch. **Nothing executable was touched: no `.jl` file was edited, nothing was run, and
  `spike/p13/consts.jl` remains BYTE-UNCHANGED.**
  **WRITTEN:** `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-D15-AMENDMENT.md` —
  the pre-registration amendment, in the `13-SC2-AMENDMENT.md` house shape. It opens by conceding
  the ordering is against it (results exist; the SC2 amendment could say they did not) and then
  carries four verified facts: (1) `c1` is the DAPI counterstain (`test/runtests.jl:105`), so the
  arm measured the wrong physical object; (2) **the strongest fact** — `spike/simulator/ghat.jl:58-68`
  and `spike/simulator/calibration.jl:351,433` had ALREADY marked `0.3292 / 0.2481` SUPERSEDED and
  already recorded the unmasked c2/c3 anchors `0.4603 / 0.3815` in `78dc37f` (2026-07-25 23:28),
  and `consts.jl` was locked in `c42cc8e` (2026-07-27 14:43) **38 h 47 min later** still quoting the
  superseded figures, so this is drift from an already-corrected upstream, not a result being
  chased; (3) **no gating threshold moves** — every untouched bar is enumerated by name and the
  gating arm (13-12) is grep-proven to contain ZERO references to either pair constant; (4) the
  consequence is **worse-or-equal** — unmasked separation 0.0788 against 0.0811 before, and the
  better-separating masked `+0.8238 / −0.0338` figures are named as THE TRAP and forbidden, because
  `patch_summary` applies no Otsu mask so masking would push already-OOD fixtures further out.
  **THE AUTHORISATION TIMELINE IS RECORDED IN BOTH DIRECTIONS AND NEITHER IS SOFTENED:** the
  *factual* correction was on the record from 2026-07-25 23:56 (`ea4a7d3`, which fixed
  `11-10-PLAN.md` in the same commit) — before `consts.jl` existed and 3 d 17 h before 13-16 ran;
  the *specific disposition* (drop the redundancy arm, use unmasked, handle as an amendment) arrived
  2026-07-29 17:48 (`f411ee7`), **41 minutes after 13-16 finished at 17:07** (`2112bed`).
  **REVISED (annotated, never rewritten — the executed record stays byte-auditable):**
  `13-01-PLAN.md` (6 sites) and `13-15-PLAN.md` (6 sites), including the `13-15` `must_have` and the
  `13-01`/`13-15` verify criteria that made **reproducing the WRONG numbers a pass condition**.
  **NEW PLAN:** `13-17-PLAN.md` (wave 10, `autonomous: false`) — amend the constants, correct the
  runner, DROP the redundancy arm everywhere (code, banner, artifact field, return value, report),
  re-run 13-16, and revise `13-REPORT.md` plus superseding notes in `13-15-SUMMARY.md` /
  `13-16-SUMMARY.md`. **The re-run consumes NO RNG stream** (six committed files, deterministic
  transform, deterministic forward pass), so no reserved counter can collide — τ-probe 1, training
  pool 2, gate 3, α-series 4, continuity 5 are all untouched, and 13-16 used none of them either.
  **The plan REQUIRES the report to answer in words: does the CONCLUSION change?** If the corrected
  pair lands in the same place the report must say *"we corrected it and the conclusion is
  unchanged"* — a GOOD result, because it shows the limitation is STRUCTURAL rather than an artefact
  of the wrong channels. Overselling the correction as a rescue is forbidden. The headline OOD number
  is EXPECTED to move (417.2974 vs 167.5446, 2.491× on `(1,2)`; a different pair already landed at
  607.57, 3.63×) and must be reported as a CHANGE, never restated.
  **13-15 does NOT re-execute, and the evidence supports the ruling:** every ingestion entry point in
  `spike/p13/real_images.jl` takes `channels` as a KEYWORD and `real_tif(cond, i)` builds both path
  and channel name from the index, so the ingestion is pair-agnostic BY CONSTRUCTION. Only its
  recorded lineage numbers are superseded, and those are re-measured inside the 13-16 re-run.
  **NOT AFFECTED, CONFIRMED BY GREP:** `run_three_way_gate.jl`, `train_three_way.jl`,
  `run_alpha_series.jl`, `labels.jl`, `datagen.jl`, `net.jl` and `alpha_series.jl` contain **zero**
  occurrences of `P13_REAL_CHANNEL_PAIR` or `P13_REAL_REDUNDANCY_PAIR`. The gate verdict (6/6, AUC
  0.990169 / 0.988262, ECE 0.0119808 / 0.012886) is untouched, `P13_ITERATION_ALLOWANCE` remains
  1 of 1 UNSPENT, and the Phase-16 seal stays SHUT.
  **ONE DECISION IS ESCALATED AND WAS DELIBERATELY NOT RESOLVED BY THE PLANNER — it is Task 1 of
  13-17, a blocking `checkpoint:decision`: WHERE do the amended constants live?**
  `run_p13_realimage.jl:634` asserts `h.consts_sha == p13_consts_sha()` where `h` is 13-11's TRAINED
  NET and `consts_sha` was recorded at train time; the identical assertion sits at
  `run_three_way_gate.jl:353` (**the gate**) and `run_alpha_series.jl:451`. **So editing `consts.jl`
  in place makes the 13-16 re-run ABORT AT LOAD and makes the gate and α runners un-reproducible
  against their own net** — and retraining is not available (not authorised; it would invalidate
  13-12). (`test_p13_net.jl:380` is NOT affected — it round-trips a net saved inside the test.)
  Two options, both defensible, both fully specified in the plan and in `13-D15-AMENDMENT.md` §7:
  **M1** edit `consts.jl` and widen the three sha guards (justified because training reads no
  `P13_REAL_*` constant — but it touches an integrity guard to produce a result, and retires the
  "two commits, neither postdating a result" claim); **M2** leave `consts.jl` byte-unchanged and add
  a separately named `spike/p13/consts_d15_amendment.jl` (every guard keeps holding, the gating arm's
  reproducibility is untouched, the superseded pre-registration stays byte-present and citable per
  `13-SC2-AMENDMENT.md` §6, and `test_p13_consts.jl:235-242` keeps asserting the frozen values
  TRUTHFULLY with the amended ones in a new testset). **The planner's reading is that M2 dominates on
  every stated constraint, and the planner did not adopt it: this is a provenance-narrative decision
  about a frozen pre-registration and it is the user's.** Under M2 the new module must be guarded on a
  symbol no other phase mirrors — see the `P13_DEV_SEED` include-guard blocker below, which is exactly
  this trap already live in the suite.

- **2026-07-29 — D-15 amendment planning passed an independent plan-check; five defects fixed; the
  M1/M2 mechanism question is STILL OPEN and was NOT resolved.** The check's verdict on the central
  question was that `13-D15-AMENDMENT.md` is **not** laundering: the "no gating threshold moves"
  claim was independently grep-verified TRUE (all seven files exactly zero for both pair constants),
  the 39-hour upstream-drift fact VERIFIED, `rng_stream = "NONE"` VERIFIED, the mask-fraction
  reportable-finding path VERIFIED REAL, the 13-15 no-re-execute ruling VERIFIED CORRECT, and the
  13-01 / 13-15 revisions confirmed surgical and byte-auditable. Fixed in this planning pass, no
  `.jl` file touched and `spike/p13/consts.jl` still at blob `70fe66df`:
  **(D1)** `spike/p13/real_images.jl` was in NO task's file list while its header blocks (b2)
  `:39-49` and (e) `:66-71` and the `real_mbar` docstring `:277-283` carry stale prose — the same
  drift class the amendment exists to repair, one level down. Now assigned to `13-17` Task 3 as a
  COMMENT-ONLY edit with explicit M1/M2-conditional wording (`13-D15-AMENDMENT.md` §6.4a).
  **(D2)** `.planning/ROADMAP.md` — the phase's most-read surface — quoted the superseded-pair
  result unmarked at four sites and was absent from the inventory. Now §6.6. The planner marked all
  four SUPERSEDED and corrected `16 plans (9 waves)` -> `17 plans (10 waves)` with a 13-17 bullet,
  as one pathspec-scoped commit after re-reading the file (a Phase-12 executor is live on it); Task
  5 adds the amended figures beside them after the re-run. **The gate sentences were not touched.**
  **(D3)** Task 5's "four honesty items survive" criterion was a single alternated `grep -c ... >= 4`
  that counts lines matching ANY branch — four lines saying only "Phase 16" would have satisfied it
  while the other three were deleted. T-13-72, the plan's central anti-rescue guarantee, rested on
  it. Replaced with four separate greps, each `>= 1`.
  **(D4)** Under M2 the not-a-gate source gate would have gone blind: `test_p13_real.jl:336` builds
  its `P13_REAL_*` name set from `consts.jl` alone, although `:86-88` already discovers every `.jl`
  in `spike/p13/` by `readdir` precisely "so a file added later cannot slip past the gate" (`:241`).
  Fixed by sourcing `:336` from the amendment module too, via `get(P13_CODE, ..., "")` at `:90` — one
  line, mechanism-independent — rather than by the weaker hand-written testset originally planned.
  **(D5)** Task 3 never named the `include` of the amendment module. Under M2 the runner must
  include it, and the guard-symbol hazard Task 2 flags applies to that include site, which lives in
  Task 3's file (`run_p13_realimage.jl:132-136`). Now named, guarded on
  `P13_REAL_CHANNEL_PAIR_AMENDED`, placed above the `real_images.jl` include.
  Also folded in: the check's strongest observation, now `13-D15-AMENDMENT.md` **§2.1** — **the
  amendment picks the WORSE number.** `ghat.jl:69-71` says the masked read of the corrected pair
  "separates the controls far more sharply" (0.8576) and `:79-83` RETRACTS the earlier "MASKED is
  biased NEGATIVE" framing as overstated, so upstream neither forbids masked nor prefers unmasked.
  The amendment forbids it anyway on a code property (`patch_summary` applies no Otsu mask; the net
  was trained unmasked), taking 0.0788 — narrower even than the 0.0811 it replaces — and then
  deleting an arm. Selecting against its own arm's interest when the better number was available
  and unforbidden is the check that distinguishes correction from relaxation.
  Corrected as smaller findings: §3's occurrence counts for the two CONSUMING files were low
  (`run_p13_realimage.jl` is 16 + 8, not "10 + 4"; `real_images.jl` is 11 and NOT "all as keyword
  defaults" — 5 defaults, 5 docstring echoes, 1 header line) — **the ZERO rows, which carry the
  gating-arm-is-blind claim, were exact and are re-verified exact**; five superseded-value sites
  annotated (`13-01:353`, `13-15:278, 294-295, 488, 547`); `13-15`'s frontmatter `must_have` at `:17`
  restored to ANNOTATED form (it had been REWRITTEN, contradicting that plan's own line-70 claim);
  and the timeline made exact — **38 h 46 min 18 s** (`ea4a7d3` 23:56:54 -> `c42cc8e` 14:43:12), with
  `13-01:62`'s figure re-attributed to `78dc37f` at **39 h 15 min 06 s**.
  **STILL OPEN, and the reason this entry exists: the M1/M2 ruling at `13-17` Task 1.** The plan's
  claim that "Tasks 2-5 name their M1/M2 deltas explicitly; nothing else in this plan changes with
  the answer" was **FALSE and is retracted on the page**. M1 additionally requires: **(B1)** an
  M1-conditional rewrite of `13-REPORT.md` §13, whose `git log` "two commits only … no commit that
  produced a Phase-13 result touched it", `rev-parse` `70fe66df…`, `hash-object` `5a4ea222…` and
  `sha256` `100e97a3…` rows all become verifiably false — i.e. the section that exists to prove
  nothing was tampered with must itself be rewritten — while Task 5 currently orders §13 left
  byte-unchanged; **(B2)** an explicit `channels = (1, 2)` pin on `test_p13_real.jl:144-146` against
  literal `0.3292 / 0.2481`, because it currently rides the default, so under M1 both sides move
  together, it STILL PASSES, and the c1/c2 lineage proof is silently destroyed while the executor
  believes it was kept as instructed; **(B3)** a `git diff -U0 | grep '^[+-]'` form of the
  no-threshold-moved criterion, because `P13_REAL_ANCHOR_TOL` (`consts.jl:579`) sits three lines
  from `P13_REAL_ANCHOR_MBAR` (`:576`) and appears as diff CONTEXT, failing that criterion on
  correct work. **These three are recorded, NOT resolved and NOT pre-fixed** — under M2 all three
  are no-ops, so pre-fixing would pre-decide the question. They materially raise M1's cost and the
  user must see them WITH the Task-1 choice.

- Phase 11 plan 11-07: the SC1g lambda-ablation tripwire FAILED and the prescribed diagnosis came back GREEN. Datagen (50,000 pairs, F5 mixture) took 56.37 min vs P11_DATAGEN_WALLCLOCK_CEILING_MIN = 150 (ceiling NOT exceeded, image-size arm NOT downgraded). Training completed CPU-only in 16.08 min at d_in = 129, D = 8; validation risk 15.6452 -> 6.0490 (best epoch 22), early stopping epoch 63 of 300. spike/test/test_lambda_ablation.jl then failed (exit 1, 9 pass / 8 fail): per-dataset delta-rho posterior-SD ratios 1.109 / 0.987 / 0.936 and 90% HDI ratios 1.122 / 0.984 / 0.940 against the pre-registered P11_LAMBDA_ABLATION_FACTOR = 2.502. The plan's prescribed diagnosis is green on all three legs: (a) on the REALIZED pool 0 of 50,000 samples violate |shift| <= lambda, cor(lambda,|shift|) = 0.604/0.611 vs a > 0.3 bar, cor(lambda,|chromatic_eps|) = -0.007; (b) training inputs were 129 rows with row 129 varying over 39,968 of 40,000 columns; (c) augment_input ran AFTER standardization (row 129 == encode_lambda(lambda) exactly, rows 1:128 == standardized summary exactly). A fourth diagnostic settles the ambiguity: the shift marginals, whose lambda response is analytically known, widen by 4.84/6.99/8.05 (dx) and 3.98/6.52/9.05 (dy) against an ideal prior-SD ratio of 12.0. THE LAMBDA CONDITIONING IS ALIVE AND USED; it is the rho_true width that does not track lambda (1.20 / 1.06 / 0.98). Consequence: wave 7 (plan 11-08, the SC2 ladder) stays GATED and no ladder was run; waves 8 and 9 (11-09, 11-10, 11-11) are gated behind it. Nothing was relaxed: P11_LAMBDA_ABLATION_FACTOR unchanged, LAMBDA_MAX unchanged (the R11 reduce-LAMBDA_MAX remedy was NOT executed - premise falsified, LAMBDA_MAX is Tier-1 pre-registration, and a narrower span makes the bar harder), no capacity raised, no re-seed, no re-run. p11_consts.jl and src/ byte-unchanged. P11_ITERATION_ALLOWANCE remains 1 of 1 UNSPENT. DECISION REQUIRED (user/orchestrator, not an executor call): SC1g's statistic is delta-rho width, the same quantity SC2 measures, so as authored it cannot separate "conditioning dead" from "effect null" - the separation its own header claims. Options: (i) re-scope SC1g onto the shift marginals (diagnostic d) as a corrected liveness tripwire, keeping delta-rho as an SC2 question; (ii) read the flat rho response as the SC2 result and proceed to the ladder to measure it at pre-registered scale; (iii) spend P11_ITERATION_ALLOWANCE. Reusable at zero further compute for any option keeping LAMBDA_MIN/LAMBDA_MAX and the F5 mixture: the 54 MB pool (spike/data/cache/p11/, untracked, content-hash-keyed on lambda_max) and spike/npe/p11_research_npe.jld2 (4,933,362 bytes, committed, generated 2026-07-27T13:49:47.473Z). Full detail: .planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-07-SUMMARY.md

- **Phase 12 plan 12-15 (wave 8): `SELECTED PRIOR: NONE-BEATS-ABLATION`. BLOCKS WAVE 9 ONWARD — USER RULING REQUIRED.** At matched induced lag-1 correlation through a full-rank `K = P12_K_DEV = 63` head, neither the CAR nor the GP lattice prior beats a neutralized prior at mini-spike scale. Rule 3 of the pre-registered `select_prior` fired: the best admissible arm (`:gp`, pooled per-region RMSE **1.06141**) does not beat the neutralized ablation (**1.05793**). Per-arm — car RMSE 1.43483 / coverage 0.99981 / **inadmissible**; gp 1.06141 / 0.89666 / admissible; none 1.05793 / 0.95872. **The load-bearing finding is not the ordering: ALL THREE ARMS HAVE NEGATIVE SKILL** against a constant-zero predictor (car −0.44305, gp −0.06856, none −0.06181; trivial RMSE ≈ 0.9943 = the field's own marginal sd). The best net is 6.2 % *worse* than predicting nothing. **NO TIER-2 CONSTANT WAS APPENDED AND NO SENTINEL REMOVED** — `P12_CHOSEN_PRIOR`, `P12_K_PROD`, `P12_D_PROD`, `P12_N_LOW` do not exist; all three negative assertions stand; `p12_consts.jl` byte-unchanged; `src/`, `spike/Project.toml`, `spike/Manifest.toml`, `corpus/`, `spike/data/cache/p11/` untouched. Budgets, separately and each against its own constant (both read 150, which is what hides a conflation): `elapsed_min` = **10.667** vs `P12_MINISPIKE_WALLCLOCK_CEILING_MIN` = 150, in budget, `BLOCKER` guard did not fire; `datagen_min` PER ARM, each judged independently inside `generate_p12_pool` vs `P12_DATAGEN_WALLCLOCK_CEILING_MIN` = 150 — car 0.084, gp 0.000, none 0.000 (pools already complete, resume-by-skip) — their sum is not a quantity anything compares to 150. No contention inflation, measured not assumed (train+score 10.574 vs the first run's 10.156, ratio 1.041). **THE PREDECESSOR'S RUN (2026-07-31T09:51Z) WAS HARVEST-READY AND SCIENTIFICALLY VOID, AND WAS DISCARDED UNHARVESTED.** It scored STANDARDIZED posterior draws against RAW truth — `run_p12_minispike.jl` never applied `bundle.theta_zt`, violating the contract written verbatim at `spike/npe/infer.jl:35-38`. It passed every structural check (all 40 keys, index disjointness, per-region pooling, budget reconciling to 0.006 min) **because none of those look at scale**. It selected **CAR** with `beats_ablation = true`, `K_PROD = 35`, `N_LOW = 3`; the corrected run reverses all of it and the two spatial arms swap admissibility (car 0.92730 → 0.99981 out of band; gp 0.99983 → 0.89666 in band). Had it been harvested, `P12_CHOSEN_PRIOR = :car` and `P12_N_LOW = 3` would have entered an append-only pre-registration permanently — and `n_low = 3` was not "the field is smooth" but "the bar equals the field's own marginal sd because the posterior carries no information in the space it was scored in", a number that would have reached the manuscript as a physical finding about biology. The truncation curve is byte-identical across both runs (arithmetic on the drawn field ensemble, touching no net) — the one quantity the defect could not reach, which makes the diagnosis a control rather than a coincidence. **SECOND FINDING, RECORDED NOT FIXED: training is NOT seeded.** `build_p12_estimator` sets no RNG and `NeuralEstimators.train` is called without one (`train_p12_npe.jl:465-469`), so Flux weight init draws from the unseeded global RNG; only the masking streams and the deterministic tail-block split are reproducible. The margin deciding rule 3 is **0.33 %** and is therefore not a stable quantity — the verdict does not rest on it (NONE-BEATS-ABLATION holds under any reordering of `:gp`/`:none`, because the negative skill carries it) but **no ranking among the three arms is established by this run**. This deviates from the project's fixed-seed constraint; it is 12-14's already-committed surface. **DECISION REQUIRED (user, not an executor or orchestrator call): 12-17 must not be started — it is not merely inadvisable but impossible as written, since `12-17:137` generates the 50,000-pair pool with `arm = P12_CHOSEN_PRIOR`, a constant that does not exist on this branch; guessing it burns ~73 min into a directory keyed on a guess. Waves 10-12 (12-16, 12-18/12-19, 12-20) consume a trained spatial bundle this run does not license. The D-13 ablation is now the likely deliverable.** Options: (i) accept the mini-spike result and re-scope the phase onto the D-13 ablation deliverable; (ii) treat the negative skill as a capacity/epoch artifact (18 epochs, 10,000 samples, 72-row full-rank head) and authorize a larger mini-spike before re-deciding — this is a compute spend and touches `P12_ITERATION_ALLOWANCE` (1, UNSPENT); (iii) seed training first so the 0.33 % margin becomes a measurable quantity, then re-run. Nothing was relaxed, no threshold moved, no arm dropped, no pool shrunk, `select_prior` byte-identical. Full detail: `.planning/phases/12-spatial-colocalization-map/12-MINISPIKE-VERDICT.md` and `12-15-SUMMARY.md`. Commits `f076da6` (fix + mechanical contract), `6176849` (artifact + verdict + summary), `605062e` (plan/briefing edits).

- **↑ RESOLVED 2026-07-31 by the user ruling recorded in `## ✅ RESOLVED 2026-07-31 (wave 8)` above: DESCOPE ONTO THE D-13 ABLATION.** Option (i) of the three listed above was taken; option (ii) was executed first as one pre-declared run and came back negative, which is what licensed (i). The entry above supersedes two claims made in this bullet: **(1) "ALL THREE ARMS HAVE NEGATIVE SKILL" does NOT survive a re-seed** — at the same 18 epochs with only the init changed, `:car` reaches skill **+0.08629**; the robust claim is instead *no spatial arm is ever both calibrated and better than the ablation*. **(2) The disqualification at the converged budget is a CALIBRATION disqualification, not an accuracy one** — at 100 epochs both spatial arms BEAT the ablation on RMSE (by 12.95 % and 4.01 %) and were rejected solely on coverage. This bullet is retained unedited as the record of the halt. **A NEW BLOCKER IS OPEN**: the Stage-2 descope has no executable route (12-14 Task 3 is keyed on `:descope`, the Stage-1 verdict is `:proceed`, and 12-16 is specified to throw with neither bundle present) — see the section above; it needs a user decision.

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 260725-vl8 | Real-anchor diagnostic measures MASKED as well as unmasked; the unqualified "negative tail not physically reachable / PRIOR-ONLY" claim in `ghat.jl` is WITHDRAWN (not inverted) | 2026-07-25 | 151ad79 | [260725-vl8-fix-the-real-anchor-diagnostic-in-spike-](./quick/260725-vl8-fix-the-real-anchor-diagnostic-in-spike-/) |
| 260725-wb7 | Real-anchor measured on the COLOC pair c2/c3 (green/red), not c1/c2 — `c1` is the DAPI counterstain (`test/runtests.jl:105`). Pair is now a named parameter. Masked: positive **+0.8238**, negative **−0.0338** (superseding 0.3292 / 0.2481). Withdrawal of the PRIOR-ONLY claim STANDS — `neg_reachable` flips to `true` on a value ≈ 0, which is a predicate artifact, not evidence | 2026-07-25 | 78dc37f | [260725-wb7-correct-the-real-anchor-diagnostic-chann](./quick/260725-wb7-correct-the-real-anchor-diagnostic-chann/) |
| 260803-jm5 | ROADMAP Phase 12 amended to **CLOSED NEGATIVE** to match `12-VERIFICATION.md` (`status: gaps_found`, 8/9, 1 FAILED). Doc-only, `.planning/ROADMAP.md` the sole file touched. Records: goal NOT achieved; `NONE-BEATS-ABLATION` ×3 so the shipped deliverable IS the ablation and no trained spatial arm exists; SPAT-06 measured **0.9707** outside **[0.87, 0.93]**; SPAT-07/SPAT-08 deferred to v2.1 (2026-07-31 rulings); real 20-plan list with 12-17/12-19 FORECLOSED and 12-18/12-20 DEFERRED. 12-11's stale `status: blocked` frontmatter documented in the roadmap and left **byte-unchanged** (why the SDK counts 5 not 4). §10.5 carried verbatim as a named limit; the three refused re-tune moves named; phase stays `partial` **by design** and re-running execute-phase 12 is NOT indicated | 2026-08-03 | d860461 | [260803-jm5-amend-roadmap-md-phase-12-to-closed-nega](./quick/260803-jm5-amend-roadmap-md-phase-12-to-closed-nega/) |

## Deferred Items

Items acknowledged and carried forward from previous milestone close:

| Category | Item | Status | Deferred At |
|----------|------|--------|-------------|
| Backend (v2) | BACK-01 Turing→RxInfer.jl feasibility evaluation | Deferred (post-Go, never a spike dependency) | 2026-06-26 |
| Summary net (v2) | BACK-02 DeepSet permutation-invariant summary upgrade | Deferred (upgrade if MLP summary proves insufficient) | 2026-06-26 |

## Session Continuity

Last session: 2026-08-04T10:02:46.294Z
Stopped at: Phase 15 context gathered
Resume file: .planning/phases/15-calibration-operating-envelope-and-ci-gate/15-CONTEXT.md
Resume action: continue the ACTIVE phase — Phase 7, plan 07-08, wave 6 of 8. Phase 8 remains COMPLETE.
NOTE: plan 07-09 (wave 7, the 32x32 grid) is `autonomous: false` and carries a long training run —
it MUST NOT be started without explicit human authorization. Separately, before Phase 16: run `bootstrap_anchor_hashes()` (corpus/anchor_rows.jl)
when online and authorized to replace the `PENDING-FETCH` sentinel on both physical anchors with real SHA-256
digests (the positive anchor is a ~6.3 GB download — an explicit human decision).
