---
phase: 12-spatial-colocalization-map
plan: 15
status: complete-blocked-downstream
wave: 8
date: 2026-07-31
---

# 12-15 SUMMARY — the mini-spike ran, and the answer is NONE-BEATS-ABLATION

**`status:` is `complete-blocked-downstream` deliberately.** This plan completed every task it owns
and reached a pre-registered outcome. What it blocks is *wave 9 onward*, not itself.

---

## 1. Outcome

`SELECTED PRIOR: NONE-BEATS-ABLATION`. Neither the CAR nor the GP lattice prior beats a neutralized
prior at mini-spike scale. **No Tier-2 constant was appended and no sentinel was removed** — which is
what the branch requires and what the two separate sentinels were reserved to record.

| arm | RMSE | trivial (predict 0) | **skill** | coverage | admissible |
|---|---|---|---|---|---|
| `:car`  | 1.43483 | 0.99430 | **−0.44305** | 0.99981 | no |
| `:gp`   | 1.06141 | 0.99331 | **−0.06856** | 0.89666 | yes |
| `:none` | 1.05793 | 0.99634 | **−0.06181** | 0.95872 | — |

Rule 3 fired: the best admissible arm (`:gp`, 1.06141) does not beat the ablation (1.05793).

**The load-bearing finding is not the ordering — it is that all three arms have NEGATIVE skill.**
Every net is a worse predictor of the drawn field than a constant zero; the best is 6.2 % worse than
predicting nothing. Full numbers, both budgets, and the branch trace: `12-MINISPIKE-VERDICT.md`.

---

## 2. THE DEFECT THAT MADE THE FIRST RUN VOID

The predecessor's run (`2026-07-31T09:51Z`) was complete, internally consistent, and **scientifically
void**. It was never harvested.

`train_p12_npe` trains on standardized θ (`:458`) and carries `theta_zt` in the bundle (`:476`) for
inversion. `run_p12_minispike.jl` never applied it — `StatsBase` was not imported and
`reconstruct` never called; the sole grep hit was the English word "reconstructing" in a docstring.
**It scored standardized posterior draws against raw truth.** The contract it violated is written
verbatim at `spike/npe/infer.jl:35-38`: *"reading a standardized draw as ρ would be silently wrong."*

**It passed every structural check** — all 40 keys, index disjointness (0 overlap, three arms), 64
per-region values pooling exactly to the reported scalars, budget reconciling to 0.006 min — because
**none of those look at scale**.

Three independent lines converged, and no one of them alone would have been conclusive:

1. the missing call;
2. GP scoring **1.5894** against a trivial-zero RMSE of **0.9943** — 60 % *worse* than predicting
   nothing, which a trained net cannot be if scored in the right space;
3. `r1_shrinkage = 4.298 = 1.117 / 0.2598` — a standardized posterior's sd is ≈ 1, so that ratio *is*
   the signature of the missing inverse.

### What would have been written into an append-only pre-registration

`P12_CHOSEN_PRIOR = :car`, `P12_K_PROD = 35`, `P12_N_LOW = 3`, `beats_ablation = true`. The corrected
run **reverses the substance of all of them**, and the two spatial arms swap admissibility (CAR
0.92730 → 0.99981, out of band; GP 0.99983 → 0.89666, in band).

**`n_low = 3` was not "the field is smooth". It was "the bar equals the field's own marginal sd,
because the posterior carries no information in the space it was scored in."** That number would have
entered the manuscript as a physical finding about biology. It is a property of an untrained net.

The truncation curve is **byte-identical across both runs** — arithmetic on the drawn field ensemble,
touching no net. It is the one quantity the defect could not reach, and it is what makes the
diagnosis a control rather than a coincidence.

---

## 3. The two evidential cases — why the machinery is worth anything

A pre-declared tie-break that only fires when it agrees with the numbers is not a tie-break; it is a
preference wearing a rule's name. Two cases are the only evidence this rule would act **against the
more attractive number**:

- **a tie resolving to CAR when GP's RMSE is nominally lower** (0.4975 vs 0.5000, inside the 1 % band);
- **a miscalibrated arm with the LOWER RMSE being DISQUALIFIED** (GP at 0.30 against CAR's 0.70, and
  CAR still selected, because coverage 0.99 is outside the band).

**Neither arose in the real run** — GP was both miscalibrated *and* worse in the invalid run, and in
the corrected run rule 3 returned before the tie-break was reached. So the synthetic tables are the
*only* place that property is established.

**Commit `4287862` claimed all five branches were exercised, but nothing in `spike/test/` referenced
`select_prior`.** They had been demonstrated ad hoc in a shell command and lost. A demonstration that
is not persisted is not a test. Now written, with 19 assertions covering: both-admissible either way,
the two evidential cases, both-inadmissible, ablation-wins, the `n_test ≥ P12_STAGE2_N_MIN` guard, a
missing arm, and purity.

**Deviation:** these went in `spike/test/test_p12_train.jl`, not the plan's declared
`test_p12_consts.jl`. The runner transitively loads NeuralEstimators, and `test_p12_consts.jl` is the
Tier-1 pre-registration gate and the **first** file the aggregator runs — making it depend on the
neural stack would mean a NeuralEstimators problem breaks the check that the frozen constants are
intact. `test_p12_train.jl` already loads that stack.

---

## 4. The contract is now mechanical, and it is scoped

Prose did not hold: the contract existed, was honoured by `infer.jl:110,122` and `benchmark.jl:188`,
and was violated silently by the one consumer that did not read it. A sweep over
`spike/validation/run_p12_*.jl` now asserts every runner that samples a posterior also calls
`reconstruct(`, and — **only where the runner's own source reports an rmse** — reports `trivial_rmse`.

**The conditional half is deliberate and is documented at the assertion.** An unconditional demand
would go RED on runners that are *correct*: 12-16 reports coverage and a log score, 12-18 reports
ranks. A rule that fires on correct code teaches people to suppress it. (I wrote the unconditional
version first and had to correct it — the same error one layer down.)

Scope was checked against the alternative: a blanket "every file calling `sampleposterior` must call
`reconstruct`" would be red on `test_p12_architecture.jl:264` (asserts only `size(draws)`) and on
Phase 11's `test_lambda_ablation.jl` (reports *ratios* of two posterior sds of the same row, which an
affine z-transform leaves invariant).

**`run_p12_minispike.jl` is the only `run_p12_*` runner that samples a posterior at all** — the other
four are ridge estimators. That is the answer to "how did this go unnoticed": there was no sibling to
be inconsistent with.

The sweep is guarded against vacuity (non-empty glob, ≥1 file actually sampling, minispike named).
**Both new guards were falsified**: removing the `reconstruct` line turned them red (2 failed, 0
passed), and the runner was reverted byte-exactly — `git hash-object` `59f98780…` before and after.

---

## 5. A SECOND finding: training is not seeded

`build_p12_estimator` sets no RNG and `NeuralEstimators.train` is called without one
(`train_p12_npe.jl:465-469`), so Flux weight init draws from the unseeded global RNG. Only the masking
streams and the deterministic tail-block split are reproducible.

**Consequence for this verdict:** the margin that decided rule 3 is **0.33 %** (gp 1.06141 vs none
1.05793) and is *not* a stable quantity. The verdict does not rest on it — `NONE-BEATS-ABLATION`
holds under any reordering of `:gp` and `:none`, because what carries it is that every arm has
negative skill. But **no ranking among the three arms is established by this run**, and the summary
says so rather than quoting the ordering as a result.

This deviates from the project's fixed-seed reproducibility constraint. It is 12-14's surface, already
committed. **Recorded, not fixed** — fixing it would change a committed training surface mid-phase.

---

## 6. Budgets — reported separately, each against its own constant

- **`elapsed_min` = 10.667** vs `P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150`. Within budget; the
  `BLOCKER` guard did not fire.
- **`datagen_min` PER ARM**, each judged independently inside `generate_p12_pool` vs
  `P12_DATAGEN_WALLCLOCK_CEILING_MIN = 150`: car 0.084, gp 0.000, none 0.000 (pools already complete;
  resume-by-skip). **Their sum is not a quantity anything compares to 150.**

Both constants read 150, which is exactly the condition that hides a conflation. They are never added
and never crossed.

**No contention inflation, measured not assumed:** train+score 10.574 min this run vs 10.156 min in
the first — ratio 1.041. I earlier speculated the run was inflated by my own concurrent probes;
process inspection and this ratio show it was not, and I am not recording a cause for 4 %.

---

## 7. Pre-registration audit trail

```
spike/validation/p12_consts.jl        UNCHANGED (git diff --quiet HEAD: clean)
sentinels present:                    3 of 3  (:P12_CHOSEN_PRIOR, :P12_N_LOW, :P12_FISHERZ_NEFF)
src spike/Project.toml spike/Manifest.toml corpus   UNCHANGED
spike/data/cache/p11                  untouched
```

The plan's own verify prints `NONE-BEATS-ABLATION no-append`. `test_p12_train.jl`: **107 pass, 0
fail** (31 new).

---

## 8. What this blocks — WAVE 9 ONWARD, FOR THE USER TO RULE ON

- **12-17 must not be started.** Not merely inadvisable — *impossible as written*: `12-17:137`
  generates the 50,000-pair pool with `arm = P12_CHOSEN_PRIOR`, a constant that does not exist on this
  branch. Guessing it burns ~73 min into a directory keyed on a guess.
- **12-16, 12-18/12-19, 12-20** consume a trained spatial bundle this run does not license.
- **The D-13 ablation is now the likely deliverable.**

This is a **pre-registration outcome, not an executor call**. Written to STATE.md as a blocker.

Plan documents edited under the orchestrator's ruling (+2/−0 each, zero deletions): `12-16`, `12-18`,
`12-20` carry the un-standardization requirement; `12-20` additionally carries `trivial_rmse` (it is
the only one of the four with a literal per-region RMSE). `12-19` was deliberately **not** edited — it
adds no runner, only the `--real` branch of 12-16's file.

**"Pitfall 5" is overloaded** and is now recorded in `12-EXECUTOR-BRIEFING.md`: spike-wide it is the
θ-un-standardization contract; `12-RESEARCH.md:867` numbers a different one, which is what
`12-16-PLAN.md` already means by it. Cite `spike/npe/infer.jl:35-38`, never the number.

---

## 9. Commits

| # | scope |
|---|---|
| 1 | `run_p12_minispike.jl` — the un-standardization fix + `trivial_rmse`/`skill`; `test_p12_train.jl` — `select_prior` branch tests, the Pitfall-5 guard, the sweep |
| 2 | `p12_minispike_report.jld2`, `12-MINISPIKE-VERDICT.md`, `12-15-SUMMARY.md` |
| 3 | plan/briefing doc edits |
