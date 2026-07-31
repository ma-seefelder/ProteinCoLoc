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

**SPAT-07** is not deferred: 12-16 carries it in reduced scope (`spat07_scope = :reduced_descope`) —
**if and when 12-16 is implemented.**

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
