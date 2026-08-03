# 12-D13 — STANDALONE AUTHORISATION FOR THE D-13 ABLATION

**The user authorised production of the D-13 ablation deliverable on 2026-07-31.** This document is
the authorisation record. It is a **NEW record, not an edited one**: no gate file was modified to
create a route.

---

## 1. Why no plan step fires — the gate divergence

**The descope trigger fired at a decision point the plans' executable routing does not model.**

`12-CONTEXT.md:310-312` defines **D-12 Stage 1** as the mini-spike itself:

> **Stage 1 (cheap, before the training spend):** at the CAR-vs-GP mini-spike, require the better
> prior to show a stated coverage improvement over the neutralized-prior ablation on simulated data.
> **If neither does, descope before the full training run.**

**That is verbatim what happened** — 12-15 returned `NONE-BEATS-ABLATION` in three runs. But the
**implemented** Stage-1 gate is a *different measurement*: 12-11's ridge borrowing probe, whose
verdict `12-STAGE1-VERDICT.md:3` reads `VERDICT: PROCEED`, and which is the only thing
`p12_stage1_verdict()` reads.

**The specified gate and the implemented gate diverged, and the trigger fired on the specified one.**
One name — "D-12 Stage 1" — meaning two different measurements in two places. That is why:

| | |
|---|---|
| **12-14 Task 3** — the only executable producer of the D-13 bundle | keyed on `p12_stage1_verdict() === :descope`; on `:proceed` it **does nothing**, and it already ran as a no-op (`12-14-SUMMARY.md:31`) |
| **12-17** | impossible as written — `12-17:137` keys its 50,000-pair pool on `P12_CHOSEN_PRIOR`, which does not exist |
| **12-16** | resolves by **file existence**; with neither bundle present it is specified to **throw** (`12-16-PLAN.md:81`) |

**Recorded as a §7 candidate.** This is the same family as the entries already in
`12-STAGE1-VERDICT.md` §7: a name that denotes two different things in two places, caught only
because the discrepancy was checked against something outside itself.

## 2. What was NOT done, and why

**`12-STAGE1-VERDICT.md` is NOT touched and remains `:proceed`, byte-unchanged.**

Writing `VERDICT: DESCOPE` into it would have made 12-14 Task 3 fire and required no new code. **It
would also have been falsifying a gate record to obtain a route** — the Stage-1 adjudication is on
the record, was ruled by the user on independent evidence, and did not change. The route is therefore
a **new record**, which is the same discipline this phase has applied to every withdrawn claim: add,
never overwrite.

Also not done: running Task 3's `:descope` body under a `:proceed` verdict as though the branch had
fired, and inventing a `stage2_gate` value that 12-16 does not define.

## 3. Why 12-16 will still find the bundle without any gate value being invented

**12-16 resolves its bundle by FILE EXISTENCE, not by the verdict** (`12-16-PLAN.md:77-81`): score
`p12_train_full_report.jld2` if present (`:full_two_arm`); **else** the bundle recorded in
`p12_ablation_report.jld2` (`:descope_ablation_only`); else throw. So writing that artifact under its
contracted name, with `bundle_path` carrying its contracted meaning, is **sufficient** — the consumer
side of the descope route was always intact, and only the producer trigger was mis-keyed.

On that branch 12-16 records `stage2_gate = :not_applicable_descope` (there is nothing to compare, and
a "fail" would be read as evidence against the spatial prior, which would be wrong) and
`spat07_scope = :reduced_descope`.

**Status note:** **12-16 has not been implemented.** No `p12_coverage.jl`, `lro_pass`, runner, or
report exists. Scoring this bundle requires executing the whole 12-16 plan, which is a separate scope
and a separate decision.

## 4. THE SCALE LIMIT, NAMED HERE RATHER THAN LEFT TO BE DISCOVERED

**This D-13 is the SPIKE-SCALE fallback: `P12_MINISPIKE_N` = 10,000 pairs, NOT the 50,000-pair
version 12-17 would have built.** `12-CONTEXT.md:323-324` already warns:

> **Caveat:** if the Stage-1 trigger fires, the ablation exists only at spike scale — this fallback is
> fully available only at the Stage-2 gate.

The same reduction applies on this route. The deliverable is still strictly more than `LocalColocMap`
offers today — a per-region Δρ map **with** per-region uncertainty, where today there is none — but it
is the reduced release, and any claim made from it inherits that limit.

## 5. Deviations from 12-14 Task 3's schema, declared

1. **`stage1_verdict = :proceed`, not `:descope`.** Written from `p12_stage1_verdict()` at run time.
   Task 3 hard-codes `:descope` because it can only run on that branch; recording it here would be
   false. `authorisation` and `authorisation_doc` name the real route instead.

2. **Task 3's `realized_r1_quantiles` acceptance criterion cannot be met by Task 3's own
   instruction — a plan defect, not an executor choice.** The criterion demands the quantiles be
   *"degenerate at `P12_ABLATION_R1`"* (0.0). The instruction,
   `generate_p12_pool(P12_MINISPIKE_N; arm = :none)`, leaves `r1 = nothing` — **drawn** from
   `Uniform(P12_R1_MIN, P12_R1_MAX)`. **Verified by loading 12-15's existing `:none` pool: r1 =
   0.078, 0.120, 0.908.** Neutralization here is **structural, not a pin**: `p12_lattice.jl:311`
   returns `Symmetric(Matrix(1.0I))` for `arm === :none` *whatever r1 is*.

   **This is a Family-B instance** — a check built for a different mechanism than the one in use. The
   criterion's *purpose* (evidence the prior really was neutralized) is met by the arm itself; its
   *literal test* is a proxy this generator never touches.

   **Resolved by following the instruction and recording the truth**, not by pinning r1 to make the
   check pass — that would be moving the mechanism to satisfy a criterion. The drawn-r1 configuration
   is also **the one 12-15 measured**, so the deliverable matches the ablation whose numbers were
   reported. The artifact carries honest quantiles, a
   `realized_r1_degenerate_at_ablation_r1 = false` flag, and a `neutralization_mechanism` field.

3. **DECLARED DEVIATION FROM THE PLAN'S DEFAULT — EPOCHS = 100, NOT 12-15's 18. A LATER READER MUST
   NOT ASSUME PLAN DEFAULTS WERE USED.** The re-run established that at an 18-epoch budget every arm
   ends **on its budget** rather than on convergence, while at 100 the `:none` arm converges and
   early-stops. **A shipped deliverable should be trained to convergence, not halted by its budget.**

   **REALIZED, FROM THE PERSISTED `risk_trace`: the budget was 100 and training stopped at epoch
   27.** (`risk_trace` is 28 × 2 — row 1 is the initial validation risk, then 27 epochs.) **The real
   cost is 27 epochs, not 100**; the budget is a ceiling that was not reached, which is the whole
   point of setting it above the expected stopping point. Final train 73.549 / val 90.382,
   **val/train = 1.2289**.

   Note the stopping epoch differs from the 36 observed for `:none` in the 12-15 treatment. That run
   trained on 9,000 samples of a **pinned 512²** pool; this one trains on the full 10,000 of the
   **F5 image-size mixture**. Different data, so a different stopping point — expected, and recorded
   so the two are not read as inconsistent.

## 6. `P12_ITERATION_ALLOWANCE` — COUNTED AS **SPENT**

**The user ruled it SPENT.** Recorded here so nobody later reads an unspent allowance as headroom
that was never there.

`p12_consts.jl:489-494` authorises **one** documented iteration and says it *"is cheapest spent
BEFORE the mini-spike — on the D-12 Stage-1 arm — never on relaxing a threshold after training."*

**The 100-epoch re-run was a documented iteration on precisely that arm.** It re-ran the mini-spike
with one axis changed under a reading pre-declared before the numbers existed. **That it was also a
defect fix (unseeded init, unpersisted risk trace) does not make it not an iteration.** The
convenient reading — "auditability fixes don't count" — was available and was not taken.

**No threshold was relaxed, no bar moved, `select_prior` is byte-identical, and `p12_consts.jl` is
byte-unchanged.** The allowance was spent on *evidence*, which is what it was reserved for. **A second
iteration is NOT authorised**, and `p12_consts.jl:490-491` is explicit that it cannot be authorised by
amending that file.

## 7. SPAT-08 — **DEFERRED TO v2.1**, as a decision taken on 2026-07-31

`12-11-PLAN.md:255-259` directs that SPAT-08 be recorded as **DEFERRED TO v2.1** *"on a Stage-1
descope"*. **That wording does not literally apply here** — this is not a Stage-1 descope. **So the
deferral is EXTENDED on the record rather than stretched in silence:** the user's ruling of
2026-07-31 defers SPAT-08 to v2.1 on the Stage-2 descope path as well.

**What is thereby NOT delivered, stated plainly:**

- **SPAT-08's radial-energy, offset-grid and ε = 0 guards do not run.** They are claimed only by
  12-20, which this route does not reach, and `12-11-PLAN.md:256` records that they **cannot be
  carried anywhere else** — 12-06's ε ridge covers only the identifiability half.
- **Therefore the phase makes NO radial-confound claim.** The S-4 warning carried forward from
  Stage 1 (`12-STAGE1-VERDICT.md` §6, `radial_r2`) is **left standing and unaddressed**, and nothing
  in the D-13 release should be read as having ruled a radial confound in or out.
- **The D-13 ablation ships without the S-4 guard suite**, which `12-11-PLAN.md:258` calls *"a real
  reduction in the honesty apparatus"*. It belongs in the memo, not in a later discovery.

**SPAT-07** — *this paragraph as first written said "is not deferred: 12-16 carries it in reduced
scope". **THAT IS SUPERSEDED BY §9**, written after a later user ruling the same day. The original
sentence is corrected here rather than left standing, because unlike the epochs figure it is not a
stale estimate but a statement that now contradicts the decision on record.*

## 8. What travels with this deliverable, unchanged

- **The qualifier "at mini-spike scale" stays on the 12-15 result.**
- > **A comparison between two failed arms is not evidence about priors — the CAR-vs-GP question is
  > UNRESOLVED, not answered.**
- **The disqualification was a CALIBRATION disqualification, not an accuracy one.** At 100 epochs
  **both spatial arms beat the ablation on RMSE** (12.95 % and 4.01 %) and were rejected **solely on
  coverage**. Coverage crossed the band [0.87, 0.93] without landing in it — all three arms
  over-cover at 18 epochs, all three under-cover at 100. **On this evidence admissibility tracks
  training duration, not the prior.** The honest claim is *never both calibrated and better — **as
  gated***, by a gate whose behaviour here follows training duration. **Named limit, and a §7
  candidate: the same shape as the Stage-1 control ceiling — a criterion whose pass/fail tracks
  something other than what it was meant to measure.** Not pursuing it is a choice this record shows
  rather than hides, and the reason is: **hunting the epoch count that lands coverage inside the band
  would be tuning a model until it passes its own calibration gate — a gate that has then stopped
  measuring anything.**
- **The withdrawn claim.** `12-MINISPIKE-VERDICT.md` §1's "no arm learned the field at all / all arms
  negative skill" **does not survive a re-seed** — `:car` reached **+0.08629** at the same 18 epochs.
  It is a property of one initialisation. The robust claim is *no spatial arm is ever both calibrated
  and better than the ablation*.
- **Init variance.** Same pool, same epochs, only the init different: `:car` **−36.68 %**, `:gp`
  **+51.05 %**, and the two spatial arms **swap which is worse** — two orders of magnitude larger than
  the 0.33 % margin that decided the earlier rule-3 branch. **LIMIT: n = 1 versus n = 1. It bounds
  nothing.**

---

## 9. SPAT-07 — **DEFERRED TO v2.1**, and the evidence that this is genuinely cheap

**User ruling, 2026-07-31 (later the same day than §7, and superseding the SPAT-07 sentence there).**
The authorised 12-16 scope is **SPAT-05** (`p12_coloc_map`) and **SPAT-06** (leave-region-out
predictive coverage in Fisher-z, with `P12_FISHERZ_NEFF` appended). **SPAT-07 is deferred to v2.1
alongside SPAT-08.**

### 9.1 What is thereby NOT delivered, stated plainly

- **No Gaussian-space SBC, no randomized-rank ρ-space SBC, no nuisance-appropriate equivalence
  testing.**
- **THE PHASE MAKES NO PER-REGION CALIBRATION CLAIM BEYOND PREDICTIVE COVERAGE.** What SPAT-06
  establishes is that the leave-region-out predictive intervals cover what was actually measured. It
  does **not** establish that the posterior is calibrated in the SBC sense, and the two must not be
  conflated in the memo.
- `spat07_scope` is recorded as **deferred**, not `:reduced_descope`.

### 9.2 THE EVIDENCE THAT SPAT-07 EXTENDS RATHER THAN REBUILDS

Recorded because *"we deferred it and it will be cheap to add"* is a claim a future reader is entitled
to see backing for, rather than a reassurance.

1. **It is a different STATISTIC over the SAME machinery.** SPAT-07's reduced pass is
   `sbc_uniformity` from the **already-existing** `spike/validation/sbc.jl`, computed over the target
   rows in Gaussian space. It consumes the same posterior-read path (`lro_arm` / the sample pass) that
   SPAT-05 and SPAT-06 require and that 12-16 Task 1 builds regardless. **No new inference machinery
   is implied by adding it later.**
2. **12-16 was designed for it to stand alone.** `12-16-PLAN.md:91-93` explicitly forbids branching on
   whether 12-18's `p12_sbc.jl` exists — *"on the descope route 12-18 is skipped, so it provably will
   not, and a dead branch invites an executor to assume the full path is available"* — and directs the
   pass be computed directly from `sbc.jl`. **The SPAT-07 half was written to have no dependency on
   the half that is skipped.**
3. **The artifact absorbs the addition backward-compatibly.** `spat07_scope` is already a required key
   whose value merely changes; the SBC arrays are **added** keys. This is the same pattern as adding
   `risk_trace` and `train_seed` to the training bundle in `faa1871`, where `schema_version` was
   deliberately **not** bumped because the loader checks for the *presence* of required keys, so an
   added key is compatible while a bumped version would invalidate every existing artifact.

**Conclusion: LOW LOCK-IN. Building only SPAT-05 + SPAT-06 now does not force a rebuild when v2.1
picks up SPAT-07 — it forces an extension.** The cheaper path is genuinely cheaper in total, not just
sooner.

### 9.3 The one piece of genuinely unexercised work, named as such

**`lro_pass` is built but CANNOT FIRE on this route.** It compares a spatial arm against the ablation;
there is no spatial arm, so `stage2_gate = :not_applicable_descope` **by construction**.

It is built anyway, deliberately: small next to `p12_coloc_map`, it keeps the artifact contract intact
for 12-20 Guard 3 (`headline_logscore_delta`), and it is the first thing a v2.1 spatial arm needs.
**Its tests exercise the BRANCH LOGIC, not the gate, and must say so in their names or comments — a
green test must never be read as the gate having run.**

---

## 10. APPENDIX, added 2026-08-03: WHERE THE ABLATION'S EXEMPTION FROM THE CALIBRATION GATE ACTUALLY LIVES

**APPENDED, NOT MERGED. §8 ABOVE IS UNTOUCHED AND MUST STAY UNTOUCHED**, including its sentence that
*"all three arms over-cover at 18 epochs, all three under-cover at 100"*. §8 treating the three arms
alike on coverage is **the evidence for what follows**; editing it to agree with the code would erase
the finding instead of recording it.

**The question.** §8 disqualifies both spatial arms on coverage against
`P12_COVERAGE_NOMINAL ± P12_STAGE2_COVERAGE_TOST_DELTA` = **[0.87, 0.93]**. The ablation's own coverage
on the repaired control is **0.98694** — also outside that band. Is the ablation exempt, and where is
the exemption written?

**The answer: nowhere in this document, and nowhere in any document.** It exists as the range of a
`for` loop at `spike/validation/run_p12_minispike.jl:249`, which tests admissibility over `(:car, :gp)`
and never over `:none` — while the pre-registered prose at `12-15-PLAN.md:127` and `select_prior`'s own
docstring both say *"An arm is admissible only if…"*, unqualified.

**It did not manufacture the verdict.** Held to the same gate the ablation fails too, giving *"no arm
is admissible"* — which selects no spatial prior either. `NONE-BEATS-ABLATION` survives.

**Full record, including why no change to the loop is proposed:** `12-STAGE1-VERDICT.md` §7.6, filed
as an instance of the §7.2 family.

### §10.1 USER RULING, 2026-08-03: **CORRECT THE RECORDED REASON. LEAVE THE VERDICT.**

**Appended, not merged. §10 above stands, and §8 remains byte-unchanged.** This subsection exists
because §10 as first written carried only *half* of what the repaired numbers support, and half of it
reads as exoneration.

**`NONE-BEATS-ABLATION` STANDS AS THE OUTCOME AND IS NOT REOPENED.** What is corrected is the *reason
on record*, which as written implies the spatial arms failed something the ablation passed. **On the
repaired basis that reading is false.**

**The repaired treatment (100 epochs, `eb3ced3`), on its face:**

| arm | coverage | distance from nominal 0.90 | rmse | skill |
|---|---|---|---|---|
| `car` | 0.94920 | **0.0492** | 0.63941 | +0.35693 |
| `gp` | 0.96503 | 0.0650 | 0.65540 | +0.34018 |
| `none` *(selected)* | **0.95089** | **0.0509** | 0.79471 | +0.20237 |

- **All three arms are outside [0.87, 0.93].** Not two.
- **`car` is CLOSER to nominal than the selected `none`** — 0.0492 against 0.0509.
- **Both spatial arms beat the selected arm on RMSE**, by 19.5 % and 17.5 %.

> **THE HONEST STATEMENT IS THEREFORE NOT "the spatial priors are miscalibrated." IT IS: _NO ARM IS
> CALIBRATED AT 100 EPOCHS, AND THE FALLBACK WAS RETAINED WITHOUT BEING TESTED._**

**BOTH HALVES, BECAUSE EITHER ALONE MISLEADS.**

1. **The loop bound did NOT change WHICH arm was selected.** Hold the ablation to the same gate and it
   fails too; the outcome becomes *"no arm is admissible"*, which selects no spatial prior either.
   The verdict is not an artifact of the exemption, and nobody should read this as a discovered error
   in the result.
2. **The loop bound DID change WHAT THE SELECTION MEANS.** The arm that won **won by default**, and it
   **fails the same gate it was used to disqualify the others with** (0.95089). "The ablation beat the
   spatial priors" and "nothing was calibrated and the ablation was never asked" are different claims,
   and only the second is supported.

**NO CODE CHANGE. THE LOOP AT `run_p12_minispike.jl:249` IS NOT TOUCHED.** Widening it to
`(:car, :gp, :none)` *after seeing that `:none` fails it* would be changing the procedure because of
the result — the move this phase has refused at every previous opportunity. **The loop stays wrong and
documented as wrong.** Any change is a v2.1 proposal, written as a proposal, not applied here.

**NO RE-DERIVATION OF THE TOLERANCE.** Re-examining the [0.87, 0.93] band itself was offered to the
user and **declined**, on the ground that it changes a rule after the numbers are known. Same
discipline as §7.4(b).

### §10.2 THE PRE-DECLARATION'S PREMISE WAS REMOVED BY `f039729` — RECORDED HERE, NOT THERE

**`12-16-PREDECLARATION.md` IS UNTOUCHED AND STAYS UNTOUCHED.** It was written before the fact and
does not get edited after it. The supersession is recorded here instead.

**§1 of that document reads:**

> *"that second branch is **live**, not hypothetical: the D-13 ablation's measured skill on
> latent-field recovery is **−0.15084**, i.e. worse than predicting a constant zero."*

**On the repaired basis that same quantity is `+0.20237`** — better than predicting a constant zero,
not worse. **The branch §1 called "live, not hypothetical" had its factual premise removed** by
`f039729`, the repair of the smoothness-permutation defect. The sign flip is not marginal and is not
confined to the ablation: `car` −0.00383 → +0.35693, `gp` −0.10810 → +0.34018, `none` −0.15084 →
+0.20237. **Every arm crossed zero.**

**CORROBORATED INDEPENDENTLY, BY A DIFFERENT ROUTE.** 12-16 measured the ablation beating
`prior_only_floor` by **+0.19195 nats/region** (`12-16-SUMMARY.md:48`), against the pre-registered
`P12_STAGE2_LOGSCORE_MIN` = 0.02. That sits right beside the **+0.20237** skill figure. **Two
different statistics, computed by two different runners on two different quantities, agree that the
ablation is informative** — which is worth more than either number alone, because the pre-declaration's
whole concern was that a model could be *calibrated because uninformative*.

**What §1 was built to guard against therefore did not occur.** The document is not thereby wrong: it
pre-committed a discriminator against a branch that turned out not to be taken, which is what
pre-declaration is *for*. **It is superseded in its premise, not in its judgement.**

### §10.3 **MANUSCRIPT-LEVEL NAMED LIMITATION — not a phase-internal note**

**FLAGGED EXPLICITLY FOR PICKUP WITH THE v2.0 NAMED LIMITS AND HONESTY ITEMS.** This does not belong
buried in a phase verdict, and it is not a housekeeping remark. The item, in the form it should travel:

> **The spatial-prior comparison selected its ablation arm by default rather than on merit. At 100
> epochs no arm — neither spatial prior nor the neutralized ablation — met the pre-registered coverage
> band of 0.90 ± 0.03; both spatial arms were more accurate than the selected arm (RMSE lower by
> 19.5 % and 17.5 %) and one was closer to nominal coverage. The selection rule tested admissibility
> only over the two spatial arms, so the ablation was never held to the criterion that disqualified
> them. This did not change which arm was selected — held to the same gate, no arm qualifies — but it
> means the result must be read as "no arm was calibrated at this scale", NOT as "the ablation
> outperformed the spatial priors".**

Full technical account: `12-STAGE1-VERDICT.md` §7.6.

### §10.4 THE COVERAGE **DIRECTION** IN §8 IS A WRONG-BASIS ARTIFACT — and it is the THIRD finding this one repair inverted

**Appended 2026-08-03. §8 REMAINS BYTE-UNCHANGED, INCLUDING THE SENTENCE THIS SUBSECTION FALSIFIES.**
That is deliberate and it is the same rule as §10: the false sentence is the evidence, and a record
that quietly acquires the right direction teaches nothing about how it acquired the wrong one.

**§8:167-169 reads:**

> *"Coverage crossed the band [0.87, 0.93] without landing in it — all three arms over-cover at 18
> epochs, all three **under**-cover at 100. **On this evidence admissibility tracks training duration,
> not the prior.**"*

**ON THE REPAIRED BASIS, NOTHING UNDER-COVERS AND NOTHING CROSSES.**

| | `car` | `gp` | `none` | all vs. band |
|---|---|---|---|---|
| 18 epochs (repaired control, `75e3b6d`) | 0.99855 | 1.00000 | 0.98694 | **all above 0.93** |
| 100 epochs (repaired treatment, `eb3ced3`) | **0.94920** | **0.96503** | **0.95089** | **all above 0.93** |

**All three arms over-cover at BOTH epoch counts. The band is approached FROM ABOVE and never
crossed.** The intervals are too *wide* at 18 epochs and still too *wide* at 100 — just less so.

**WHAT THIS COSTS, NAMED PRECISELY RATHER THAN LEFT TO BE INFERRED.**

- **SURVIVES:** *"never both calibrated and better — **as gated**"*. That claim rests on the
  admissibility outcome, which is unchanged.
- **DOES NOT SURVIVE IN THE FORM WRITTEN:** *"admissibility tracks training duration, not the
  prior."* **Its evidence was the crossing.** An overshoot — too wide, then too narrow — is what makes
  a gate look like it is tracking something other than calibration, and is what makes *"more training
  is not the answer"* follow. **There is no overshoot.** A monotone approach from above that stops
  short is instead the ordinary signature of a model that has **not trained long enough**.

**SO THIS REMOVES ONE OF THE TWO SUPPORTS FOR "UNDER-TRAINING IS ELIMINATED". THE OTHER IS UNTOUCHED,
AND BOTH ARE STATED HERE WITHOUT BEING RECONCILED.**

- **Removed:** the coverage-direction argument above.
- **Standing:** early stopping fired at **56 / 99 / 36** epochs against a 100-epoch budget
  (`12-MINISPIKE-VERDICT.md:238`, `12-15-RERUN-VERDICT.md:253`) — training halted on its own before
  the budget in all three arms.
- **AND THAT SUPPORT IS NOW CONFIRMED BASIS-INDEPENDENT BY EXECUTION.** The repaired 100-epoch run
  reproduced the stopping epochs **exactly — 56, 99 and 36** (`treatment_100ep.log`, arms in order).
  This is the direct confirmation of §7.4(a)'s claim that the defect never reached the trained
  bundles: **training consumes θ as stored and never reconstructs a field**, so the permutation could
  not touch it. The early-stopping evidence is therefore genuinely untainted, unlike the coverage
  evidence beside it.
- **NOT ADJUDICATED HERE:** `gp` stopping at **99 of 100** is not obviously convergence, and a
  patience of 5 on a still-descending curve is a weak stopping criterion. Whether that leaves
  "under-training eliminated" standing is **not decided in this subsection** and is a user question.

**NO COMPUTE WAS RUN FOR THIS SUBSECTION AND NO REMEDY IS PROPOSED.** In particular this is **not** a
recommendation to train longer. The epoch count is the one knob whose effect on the gate has already
been observed, so choosing it now would be selecting a hyperparameter to pass a calibration gate —
the hunt §8 itself refused, and the refusal binds harder with better evidence, not less.

**THE COUNT, PUT IN ONE PLACE SO ITS SIZE IS VISIBLE.** `f039729` — a single unreversed permutation in
one scoring function — inverted **three separate stated findings** of this phase:

| # | the finding as stated | on the repaired basis |
|---|---|---|
| 1 | every arm has **negative** skill; "no arm learned the field at all" | **all three positive**; −0.004/−0.108/−0.151 → +0.357/+0.340/+0.202 |
| 2 | the spatial arms were disqualified on coverage **that the ablation met** | **no arm meets it**; the ablation was never tested (§10.1) |
| 3 | all three **under**-cover at 100 epochs; the band was crossed | **all three over-cover**; the band is never crossed |

**They are three faces of one defect, not three unrelated corrections**, and a reader who meets them
separately will never see the size of it. Each was invisible to structural verification for the same
reason: **a permutation is orthogonal**, so cardinality, pooling identities, disjointness, budget
reconciliation and even Parseval energy all survive it intact (§7.1, §7.4(a)).

**MANUSCRIPT-LEVEL, alongside §10.3.** The item as it should travel:

> **A single scoring-basis defect inverted three of this phase's stated findings — the sign of every
> skill measurement, which arms failed the calibration gate, and the direction of the coverage miss.
> All three passed every structural check the phase could make, because the defect was an orthogonal
> permutation and structural checks are invariant to it. The corrected reading is that no arm was
> calibrated at this scale, that coverage was converging toward nominal from above rather than
> overshooting, and that one of the two arguments for having ruled out under-training does not
> survive.**
