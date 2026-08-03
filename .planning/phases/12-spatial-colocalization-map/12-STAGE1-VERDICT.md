# Phase 12 — D-12 Stage 1: VERDICT

VERDICT: PROCEED

**This verdict was ruled by the user, not by the executing agent, and it overrides one half of a
frozen interpretation rule.** The reasoning is written out in full below rather than referenced,
because the honest summary of this document is uncomfortable in one specific way: *the run's positive
control did not clear the ceiling written for it, and this document rules nonetheless that the phase
proceeds.* A reader who is handed only the verdict line would be entitled to suspect a threshold was
moved after the number was seen. It was not — `spike/validation/p12_consts.jl` Tier 1 is
byte-unchanged, `P12_STAGE1_CONTROL_CEILING` still reads `0.5`, and `P12_ITERATION_ALLOWANCE` still
reads `1`, unspent. §4 is the argument; §7 is the part a referee should read.

---

## §0 Provenance

| Item | Value |
|---|---|
| Measurement date | 2026-07-29 |
| Ruling date | 2026-07-31 |
| Ruled by | the user (a pre-registration question; the executing agent has no authority over it) |
| `git rev-parse HEAD` at measurement | `e88b97f435daf1113899c5c2bd229da6c3c6720d` |
| `git rev-parse HEAD` at ruling | `d2f0124523ea59d9e6f99dfae0d0c64825b71623` |
| Artifact | `spike/validation/p12_stage1_report.jld2` |
| Artifact `generated` | `2026-07-29T19:14:32.274Z` |
| Runner | `spike/validation/run_p12_stage1_ridge.jl` (committed at `e88b97f`, **before** any pool was generated) |
| Pre-registration | `spike/validation/p12_consts.jl` — Tier 1, frozen before any Phase-12 field existed; **byte-unchanged** |
| Arm / size | `arm = :car`, `imsize = (512, 512)` pinned, `imsize_tag = :pinned_512x512` |
| z-scoring arm | `:per_row` (R-3, `P12_ZSCORE_ARM`) |
| Iterations spent | **0 of `P12_ITERATION_ALLOWANCE` = 1** — carried forward UNSPENT |
| Full blocker evidence | `12-STAGE1-ADJUDICATION-BLOCKED.md` (the halt that produced this ruling; superseded as a *decision* document, retained as the *evidence* record) |

Every number below was re-derived by loading the `.jld2`, not copied from a console log or from the
blocked document. Five pools, one per rung, `r1` **PINNED** and reaching the content hash, in five
distinct content-hash directories under `spike/data/cache/p12/` (`pool_dirs` in the artifact):

```
0.05  d5562341054cb89d254b3f4afb269856f94df2aad5aae84bce1acd016bdd5120
0.25  026051006ff8e32694d82d0a35fce26c6df59fc8a7a3e618f31fb01a3de42c47
0.50  6e1aa27e7b0fbeee439b9c02585965f5d9e2fbc2316d7622d07a8608863cc3f8
0.75  89cdd7400eea4e3a637f760b60605ae6c11bd92316e8e6ebb1936d531630684c
0.95  18a985ecbdbe1eed165df96ea0e4344281562db546da610e72b49fd7371e8d08
```

## §1 The pre-registered thresholds, quoted from the artifact's own keys

| Artifact key | Value | Tier-1 name |
|---|---|---|
| `ratio_ceiling` | **0.95** | `P12_STAGE1_RATIO_CEILING` |
| `control_ceiling` | **0.5** | `P12_STAGE1_CONTROL_CEILING` |
| `min_rungs` | **1** | `P12_STAGE1_MIN_RUNGS` |

All three were fixed in Tier 1 **before any Phase-12 field existed** — before the lattice kernels,
before the simulator edit, before this runner — and all three are members of `P12_GATING_CONSTANTS`.

`P12_STAGE1_RATIO_CEILING` carries its derivation verbatim in the file (`p12_consts.jl:358`;
n_test = 1250 per rung, sampling sd ≈ 1/√(2·1250) = 0.0200, 2.5 sd of separation from 1.0 ⇒ 0.950),
and that derivation is *re-executed* by an assertion at `p12_consts.jl:763-764`.
**`P12_STAGE1_CONTROL_CEILING` carries no derivation anywhere in the pre-registration.** Its only
in-file guard, `p12_consts.jl:765`, is the ordering constraint `0.0 < CONTROL < RATIO < 1.0` — which
0.4, 0.5 and 0.6 satisfy identically. That is a fact about the file, checked by reading it, and it is
load-bearing for §4.

## §2 The ladder, as measured

`n_per_rung = 5000`, realized `n_test = 1250` at every rung — equal to
`P12_STAGE1_N_TEST_PER_RUNG_DERIVED` and asserted at every rung before the ridge ran, so the gate has
the power its ceiling was derived for.

| r₁ | `ratio_mean` | `ratio_min` | `ratio_max` | `control_ratio_mean` | `rho_arm_ratio_mean` | `radial_r2` | rung vs 0.95 |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.92557 | 0.84632 | 1.01020 | 0.73188 | 0.92439 | 0.34522 | **PASS** |
| 0.25 | 0.84660 | 0.75726 | 0.97470 | 0.67740 | 0.84641 | 0.37843 | **PASS** |
| 0.50 | 0.71870 | 0.61900 | 0.85285 | 0.59827 | 0.71985 | 0.25225 | **PASS** |
| 0.75 | 0.56184 | 0.48333 | 0.70133 | 0.49659 | 0.56276 | 0.41096 | **PASS** |
| 0.95 | 0.33569 | 0.29381 | 0.42932 | 0.32369 | 0.32311 | 0.74890 | **PASS** |

The gated target is the **Gaussian-space** lattice value `z_field[r]` (R-2, atom-free). The
`rho_arm` column is **REPORTED ONLY**; its ratio is flattered by the elementwise `ghat` clamp's
measured 6.79 % per-region atom mass, and in the event it tracks the Gaussian arm to within 0.002 at
four of five rungs, so nothing in this reading depends on which is quoted.

Supporting measurements, same artifact:

| r₁ | `baseline_rmse_mean` | `within_image_sd_rows` | `within_image_sd_rho` | `within_image_sd_zfield` |
|---|---|---|---|---|
| 0.05 | 0.99798 | 0.24666 | 0.46954 | 0.99276 |
| 0.25 | 0.99568 | 0.26154 | 0.46229 | 0.97781 |
| 0.50 | 0.99008 | 0.26783 | 0.43346 | 0.91534 |
| 0.75 | 0.97952 | 0.23320 | 0.33716 | 0.70830 |
| 0.95 | 0.99772 | 0.17013 | 0.15571 | 0.32801 |

## §3 The two gate components, as the artifact records them

```
borrowing_ok = count(ratio_mean .<= 0.95) = 5  >=  min_rungs = 1     ->  TRUE
control_live = maximum(control_ratio_mean) = 0.73188  <=  0.5        ->  FALSE
stage1_pass  = borrowing_ok && control_live                          ->  FALSE
```

**Stated plainly, and this sentence is not to be softened anywhere downstream: the positive control
did NOT pass as written.** `control_live` is `false` in the artifact and stays `false`; nothing in
this ruling changes a stored value. The half the phase is about — `borrowing_ok` — passed at 5 of 5
rungs by a margin.

## §4 The adjudication: why a failed control produces a PROCEED

The frozen rule (`12-11-PLAN.md:231-235`) reads, verbatim:

> **`control_live` false** (control ratio above `P12_STAGE1_CONTROL_CEILING`) ⇒ the harness is not
> live and the whole run is **uninformative**. Do not record a PROCEED or a DESCOPE. Record
> "harness not live", name the most likely causes in order (a broken standardizer, a
> target/predictor misalignment, a pool with the wrong `r1`), and spend the single
> `P12_ITERATION_ALLOWANCE` on fixing the harness — not on relaxing a threshold.

The rule contains a premise and a consequence. The premise — *a control above the ceiling means the
harness is not live* — is **stated by the rule, not tested by it**. The ruling is that the premise is
false in this instance, on five pieces of evidence, four of which are in the artifact and one of
which is in the file.

### (1) The control is provably at the information limit of the data

`own_row_bound = sqrt(1 − corr(row_r, z_field[r])²)` is the RMSE-over-prior ratio of the **best
possible linear predictor** of a region's lattice value from that region's own summary row. It is a
property of the DATA. It is computed with **no ridge, no penalty, no train/test split and no
standardizer** — that is, without any part of the machinery the positive control exists to audit.
No defect in that machinery could manufacture agreement with it.

| r₁ | own-row ridge (`ownrow_ratio_mean`) | information limit (`own_row_bound_mean`) | difference | full-128 control | global-level control |
|---|---|---|---|---|---|
| 0.05 | 0.76588 | 0.76564 | 0.00024 | 0.73188 | 0.77173 |
| 0.25 | 0.71236 | 0.71201 | 0.00035 | 0.67740 | 0.61400 |
| 0.50 | 0.63716 | 0.63667 | 0.00049 | 0.59827 | 0.38141 |
| 0.75 | 0.55590 | 0.55551 | 0.00039 | 0.49659 | 0.23234 |
| 0.95 | 0.47199 | 0.47183 | 0.00016 | 0.32369 | 0.19687 |

The estimator attains the analytic limit to within **5·10⁻⁴ at every rung**. It is not merely
working; it is optimal against a bound computed without it. This is **proof of liveness**, not an
absence of evidence of death.

### (2) The control, as specified, is not a liveness test at all

**A correction to how this was stated when the ruling was framed, found by asserting it
executably.** The condensed retelling of the evidence — carried in `STATE.md`, in the blocker
briefing, and in my own first draft of this document — said *"at r₁ ≤ 0.50 the limit is 0.766 /
0.712 / 0.637, above the 0.5 bar"*, i.e. three rungs. The full ladder is:

| r₁ | own-row information limit | above the 0.5 ceiling? | full-128 control | clears 0.5? |
|---|---|---|---|---|
| 0.05 | 0.76564 | **yes** | 0.73188 | no |
| 0.25 | 0.71201 | **yes** | 0.67740 | no |
| 0.50 | 0.63667 | **yes** | 0.59827 | no |
| 0.75 | 0.55551 | **yes** | 0.49659 | **yes** |
| 0.95 | 0.47183 | no | 0.32369 | **yes** |

The own-row limit is above the ceiling at **four** of five rungs, not three — and at the fourth
(r₁ = 0.75) the control **cleared the ceiling anyway**, at 0.49659. The two facts are not in
tension, and keeping them apart is the whole point: **the own-row bound is not an upper bound on the
full-128 control**, which sees all 128 rows and legitimately beats it at every rung by borrowing
from the other 63 regions. `12-STAGE1-ADJUDICATION-BLOCKED.md` §4 also names three rungs, but it is
NOT wrong: it carries the qualifier *"from that region's own row"*, never claims those three are the
only ones, and says outright that the control clears the ceiling at r₁ = 0.75 by borrowing. What was
lost was the QUALIFIER, in condensing that section into a briefing — and losing it turned a
true-but-partial statement into an overclaim about what no estimator can do. Both tuples are now
recorded in the Tier-2 append
(`P12_STAGE1_OWNROW_LIMIT_ABOVE_CEILING_RUNGS`, `P12_STAGE1_CONTROL_CLEARS_CEILING_RUNGS`) with an
assertion that re-derives each from the measured vectors, so the file refuses to load if the
"three rungs" phrasing is ever reintroduced.

**The correction sharpens the diagnosis rather than weakening it.** Read the table by column:

- No estimator of any kind can reach 0.5 from a region's **own row** at four of the five rungs.
- The **only** thing that carries the control under 0.5 anywhere is **borrowing from the other 63
  regions** — and borrowing is precisely the effect `borrowing_ok` exists to measure.
- So the control passes at exactly the two rungs where the effect under test is strongest, and
  fails at the three where it is weakest. Its pass/fail pattern is a **monotone function of the
  quantity being gated**.

A positive control is supposed to certify that the instrument works *independently* of whether the
effect is present — that is what makes it a control. This one cannot: it is a second, harder
borrowing test wearing a positive control's name. Its FALSE therefore carries no information about
instrument health, which is the exact inference the frozen rule draws from it.

That is a stronger objection than "the bar is unreachable". Unreachability would make the control
uninformative at three rungs; **structural dependence on the effect under test makes it uninformative
as a control at every rung**, including the two where it passed.

### (3) The rationale that is actually on the record IS met, at every rung

`P12_STAGE1_CONTROL_CEILING` has no derivation. The only rationale recorded anywhere for it is
12-RESEARCH's qualitative sketch — the control *"must be far below 1.0, else the harness is dead"*.
The measured 0.324–0.732 is a **27–68 % error reduction** against the empirical prior-mean baseline,
at every rung. On the stated rationale, the control passes everywhere.

The reading that fails rests on an implicit anchor: Phase 11's **global** ρ_true control at 0.157.
`12-CONTEXT.md` S-2 warns in advance, in its own words, that that number *"does not prove per-region
recovery"* and must not be read as a per-region expectation. Applying it as one is the error the
context document pre-emptively named.

### (4) All three causes the rule names are ruled out executably

| Named cause | Evidence against it |
|---|---|
| a broken standardizer | Fitted on the FIT split only, at every rung and for every region, with the analog's zero-variance guard. Verified directly: after masking, region r's continuous and mask rows are **exactly 0.0** raw *and* **exactly 0.0** after standardization, while the unmasked row is not. The control and masked designs differ **only** in those two rows. |
| a target/predictor misalignment | `z_field` and `rho_field` are reshaped `G²×n` column-major — the identical `p12_idx` ordering `encode_d01` produces for the summary rows — so summary row r and lattice row r are the same region with no transpose. Independently confirmed by the own-row correlation being high and rising with r₁ (0.635 → 0.881); a misaligned target would give a correlation near zero at every rung. |
| a pool with the wrong `r1` | Every realized `r1` asserted `==` its pin at every rung; the five rungs resolved to five distinct content-hash directories, asserted **before** generation. A collapse to one directory was demonstrated to be detectable by the same assertion. |

And the Phase-11 benchmark reproduces where the two are comparable: `global_control_ratio` at
r₁ = 0.95 measures **0.19687** against Phase 11's 0.157, at the noisiest single image size in the F5
set. It degrades to 0.77173 at r₁ = 0.05 for the reason the physics predicts — with a rough field,
the image-mean of 64 near-independent regions is a small target against unchanged measurement noise.

### (5) The file itself distinguishes the two ceilings

Verified independently of the executing agent, by reading `p12_consts.jl`:

- `P12_STAGE1_RATIO_CEILING` (`:353`) is followed at `:358` by its derivation recorded verbatim, and
  at `:763-764` by an assertion that **re-derives** it from `P12_STAGE1_N_TEST_PER_RUNG_DERIVED`. If
  that ceiling were mistyped, the file would refuse to load.
- `P12_STAGE1_CONTROL_CEILING` (`:354`, the adjacent line) has **no derivation block**, and its only
  guard (`:765`) is the ordering constraint `0.0 < CONTROL < RATIO < 1.0` — satisfied identically by
  0.4, 0.5 and 0.6. Nothing in the pre-registration would have caught a wrong value for it.

Both nonetheless sit in `P12_GATING_CONSTANTS`. That is the defect: a constant with the authority of
a gate and the provenance of a guess.

### The ruling

> The harness is provably live, certified by evidence the gate never asked for: the own-row ridge
> attains the analytic information limit at every rung. The gate's own control cannot certify that
> independently, because the only thing that brings it under its ceiling is the very borrowing
> effect under test — it clears at exactly the two rungs where that effect is strongest, and the
> own-row limit forbids reaching the ceiling at four of five rungs without it. The criterion as
> written is **mis-scaled for the unit it is applied to** and does not test what it was meant to
> test. Its failure is not evidence that the harness is dead.

The gate accordingly reduces to `borrowing_ok`, which is **TRUE at 5 of 5 rungs** against a
`P12_STAGE1_MIN_RUNGS` of 1. Under the second interpretation rule — *`borrowing_ok` true ⇒
`VERDICT: PROCEED`* — the verdict is **PROCEED**.

Note on the shape, as that rule requires it be stated: clearing at the short-correlation end only
would have been the *expected* shape and sufficient under `min_rungs = 1`, because the gate asks
whether spatial borrowing exists **anywhere** in the prior's support. In the event the ladder cleared
at **all five** rungs, so the gate is not resting on its minimum.

### What this ruling does NOT do

- It does **not** edit Tier 1. `P12_STAGE1_CONTROL_CEILING` **keeps its value of 0.5** as the
  historical record — the same treatment `SPEEDUP_GATE = 100.0` and
  `P11_LAMBDA_ABLATION_FACTOR = 2.502` received when their bars were found mis-specified.
- It does **not** spend `P12_ITERATION_ALLOWANCE`. The ruling is that the harness needs no repair, so
  there is nothing to spend it on. It stands at **1, UNSPENT**.
- It does **not** re-run the ladder on a different seed, size or arm. No number in §2 changed.
- It does **not** rewrite `control_live` in the artifact. It remains `false`, and `stage1_pass`
  remains `false`. This document is the adjudication of those values, not a replacement for them.
- The measured own-row limits, the control ratios and this adjudication's provenance are recorded as
  a **Tier-2 APPEND** (`P12_STAGE1_CONTROL_ADJUDICATION` and companions) naming
  `spike/validation/p12_stage1_report.jld2`. Every constant in that append is a MEASUREMENT or a
  RECORD; **none of them is a new pass/fail bar**, and none enters `P12_GATING_CONSTANTS`. The
  Stage-1 gate is adjudicated once, here, by a human — it is not re-applied under a friendlier
  number.

## §5 The scientific deliverable, untouched by any of the above

**Neighbouring regions do carry information about a held-out region, and the amount rises
monotonically with the correlation length.** `ratio_mean` runs 0.926 → 0.847 → 0.719 → 0.562 → 0.336
across r₁ = 0.05 → 0.95, clearing the pre-registered 0.95 ceiling at **all five** rungs, with
`ratio_max` (the worst single region of 64) clearing it at four of five and missing only at the
shortest rung (1.010 at r₁ = 0.05).

**The crossing point, in the gate's own units.** Comparing the masked arm (the other 63 regions, no
own row) against the own-row information limit — i.e. asking *when do the neighbours become more
informative about a region than that region's own measurement*:

| r₁ | neighbours only (`ratio_mean`) | own row only (limit) | neighbours better? |
|---|---|---|---|
| 0.05 | 0.92557 | 0.76564 | no |
| 0.25 | 0.84660 | 0.71201 | no |
| 0.50 | 0.71870 | 0.63667 | no |
| 0.75 | 0.56184 | 0.55551 | **parity** (Δ = 0.0063) |
| 0.95 | 0.33569 | 0.47183 | **yes**, by a wide margin |

The crossing sits at essentially exactly **r₁ = 0.75**. 12-RESEARCH Pitfall 2 named this crossing "a
genuine scientific deliverable" and predicted its shape — improving monotonically with r₁ and
crossing "somewhere in the middle" — **before it was measured**. The crossing is quoted against the
own-row information limit, a measured quantity in this artifact, rather than against the research's
0.0662 per-draw noise sd: that figure was measured at 1376×1028 while this ladder is pinned to 512²,
and Pitfall 2 itself records the noise sd as image-size dependent, so the two are not comparable and
are deliberately not compared.

Pitfall 2 also predicted the control's shape in advance, in these words: *"under the D-05 copula
every region's marginal is `MU_PRIOR` regardless of the correlation length, so the within-image
contrast is controlled entirely by the correlation length."* The measured `within_image_sd_zfield`
falls 0.993 → 0.328 across the ladder while the control ratio falls 0.732 → 0.324, in step.

## §6 Carried forward into 12-20 — the S-4 warning, UNSEPARATED

`radial_r2` — the fraction of the 64-region recoverability spread explained by radius alone — is
0.345 / 0.378 / 0.252 / 0.411 / **0.749**. **At r₁ = 0.95 three quarters of it is radial.** This is
an early S-4 warning and must be carried into 12-20's guard interpretation.

It has a competing **benign** reading, and this run does **not** separate the two: at long
correlation lengths the per-region spread is small and dominated by lattice geometry, and 12-07
measured the CAR marginal sd to be largest at low-degree (corner and edge) cells — corner (1,1)
0.9922 vs interior (4,4) 0.6575 at α = 0.95 — which is itself a radial pattern with **no chromatic
content**. 12-20's ε = 0 and offset-grid arms are what separate them.

The confound is **LIVE**, not closed, and 12-06 is why: it measured ε to be genuinely *identified*
(`eps_ratio = 0.96307`, 5.8 sampling-sd below 1.0 at n_test = 12500 — the first nuisance in this
project distinguishable from its prior) but **pre-registered-vacuous** (`ridge_residual_shrinkage =
0.96870` against the frozen floor of 0.90, so `vacuous = true`; 96.9 % of the ε prior spread
survives). Under rules frozen before that number was seen, S-4's radial chromatic confound stays
LIVE and **the 12-20 guards carry the SC3 claim**. Named limit #4 is *sharpened* by this, not
revised.

## §7 The pattern a referee will notice, recorded so Phase 16 inherits it

**This is the fourth time in this milestone that the conclusion has been "the bar was wrong, not the
model".** Written down here deliberately, because the pattern is more visible from outside the
project than from inside it, and a reader who assembles it themselves will read it less charitably
than it deserves:

| # | Where | What was wrong with the bar |
|---|---|---|
| 1 | Phase 11, SC1g | The λ-ablation tripwire compared a **component** against a **total** — mis-specified, not failed. `P11_LAMBDA_ABLATION_FACTOR = 2.502` kept its value; the finding closed negative-but-useful. |
| 2 | Phase 12, real-image arm | A gate sized for n = 271 applied to an arm with `effective_independent_n = 2`. Resolved by REPORTING that arm, not gating it. |
| 3 | Phase 12/13, Z_TWO_SIDED_90 | A **Wald**-labelled quantity used to size a **Wilson** interval (DEF-12-04, a silent function default). |
| 4 | Phase 12, Stage 1 (**this document**) | A gating ceiling with **no derivation**, on a "positive control" whose only route to passing is the very effect it was meant to control *for*. |

### The four are TWO recurring design errors, not four unrelated slips

This is the part that is actually worth a referee's attention, and it is not visible from the list.
The four entries fall into **two families with two distinct root causes**, and naming the families
says something a list of four cannot.

**FAMILY A — a control or tripwire whose statistic IS the quantity under test.** (Entries 1 and 4.)

Such a check cannot do the one job it exists for, because its failure is indistinguishable from the
outcome it is supposed to be independent of.

| | The check | Its statistic | What it therefore could not separate |
|---|---|---|---|
| 1 | Phase 11's SC1g λ-ablation tripwire | Δρ posterior width — **the same quantity SC2 measures** | "conditioning is dead" from "the effect is null" — exactly the separation its own header claimed |
| 4 | This document's Stage-1 positive control | per-region borrowing performance — **the same quantity `borrowing_ok` measures** | "the instrument is dead" from "borrowing is weak" |

In both cases the check was, structurally, a second and harder version of the measurement it was
meant to certify, wearing a control's name. **The generalisable rule, stated so no later phase has
to rediscover it: a control must be measurable when the effect under test is ABSENT. If it is not,
it is not a control.** Applied to entry 4: the Stage-1 control could only come under its ceiling by
borrowing, so at a rung where borrowing is weak it must fail whether the instrument works or not.

**FAMILY B — a quantity compared against a bar or baseline constructed for a different purpose.**
(Entries 2 and 3, plus a third instance found after this table was first written.) The bar was
sound; it was applied to something it was not built for — an `effective_independent_n` of 2 against a
bar sized for N ≥ 271, and a Wald-derived quantity used to size a Wilson interval.
`P12_STAGE1_RATIO_CEILING` carries a recorded warning about exactly this hazard in its own
derivation block, which is why it was *not* a further instance.

**The third Family-B instance: 12-13's per-r₁-bin table.** Its per-bin ratio divides each bin's
ridge RMSE by the RMSE of the **global** prior mean *within that bin* — so the denominator is
essentially the bin's distance from that mean. Measured:
`corr(|bin midpoint − prior mean|, reported prior RMSE) = 0.99682`, and the central bin's midpoint
sits **0.0006** from the prior mean, making its baseline near-perfect by construction so that any
estimator must lose there. The 1.17 in that cell is binning, not physics.

**Family B's failure condition is checkable IN ADVANCE, which is what makes the family actionable
rather than merely descriptive:** the analogous `run_p12_eps_ridge.jl` escapes this because it bins
by `|ρ|` while *targeting* ε — two different quantities — whereas in 12-13 **the binning variable IS
the target**. Ask of any binned report: *is the variable I am binning on the same as the one I am
predicting, and is the baseline computed globally?* If both, the central bins are uninterpretable.

**12-13 did NOT re-specify that metric after seeing its output.** The plan specified the global
baseline, the headline result does not depend on the per-bin table, and the artifact records the
metric as run. The trap is documented, not tuned away — which is the whole difference between a
documented artifact and a fitted one.

**A Family-B NEAR-MISS, recorded as evidence that the family is a LIVE hazard rather than a
historical list.** In wave 8 the executor described ~11 min of datagen as leaving the net training
*"as the dominant cost against the 150-min ceiling"* — conflating **two separate budgets that
happen to carry the same number**:

```
P12_DATAGEN_WALLCLOCK_CEILING_MIN   = 150   # p12_consts.jl:476
P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150   # p12_consts.jl:477
```

They govern different things and are enforced differently: the datagen ceiling is checked **inside
`generate_p12_pool`, per call, by projection** (`p12_generate.jl:603, 665`), so each arm's pool is
judged on its own; the mini-spike ceiling governs 12-15's **three-arm training run** and has **no
enforcer in the codebase at all** — 12-15's runner is required to implement it. Caught in review
before any number was reported against the wrong bar, so it is **not** a fourth instance. It is
recorded because **identical values are one of the conditions that hide a Family-B error**: a
conflation that would be obvious at 150 vs 90 is invisible at 150 vs 150, and survives review
precisely because both numbers check out.

Two families, three root-cause instances in Family B (plus this near-miss) and two in Family A — and
every one was caught by the same mechanism: writing the claim down in advance made it checkable
against something outside itself.

**A FOURTH FAMILY-B INSTANCE, added 2026-07-31 — and it UPDATES THE COUNT in the paragraph above to
four in Family B.** Recorded here rather than by editing that sentence, so the record stays additive.

**12-14 Task 3's `realized_r1_quantiles` acceptance criterion cannot be satisfied by Task 3's own
instruction.** The criterion (`12-14-PLAN.md:410-411`) requires the quantiles be *"degenerate at
`P12_ABLATION_R1`"* (= 0.0), *"proving the ablation really was trained on the neutralized prior"*. The
instruction (`12-14-PLAN.md:377`) is `generate_p12_pool(P12_MINISPIKE_N; arm = :none)`, which leaves
`r1 = nothing` — **drawn** from `Uniform(P12_R1_MIN, P12_R1_MAX)`.

**The bar is sound; it was built for a mechanism this generator does not use.** Neutralization here is
**structural**: `p12_lattice.jl:311` returns `Symmetric(Matrix(1.0I, G*G, G*G))` for `arm === :none`
**whatever r₁ is**. Pinning r₁ is simply not how the `:none` arm is neutralized, so a test on r₁'s
distribution measures nothing about neutralization.

**Verified by execution, not by reading** — twice. First against 12-15's existing `:none` pool
(r₁ = 0.078, 0.120, 0.908, plainly not 0.0), and then confirmed by the D-13 production run itself,
whose realized quantiles are **[0.05017, 0.27345, 0.49885, 0.72227, 0.94992]** — a clean uniform,
exactly as the instruction implies and exactly what the criterion forbids.

**THIS INSTANCE WAS PREDICTED BEFORE THE RUN AND THEN OBSERVED, WHICH IS EPISTEMICALLY DIFFERENT FROM
THE OTHER FOUR AND STRONGER.** Entries 1-4 of the families above were all diagnosed *after* a number
looked wrong. Here the defect was derived from the mechanism (`p12_lattice.jl:311` neutralizes on the
arm, so r₁'s distribution cannot evidence neutralization), the uniform quantiles were **stated in
advance as the expected outcome**, the run was executed unchanged, and the prediction held.

**A taxonomy that only ever explains failures after the fact is a narrative; one that predicts the
next instance is a finding.** Phase 16 should report this instance as the family's first *predictive*
confirmation, and should say plainly that the prediction was recorded in the runner's header and in
`12-D13-AUTHORISATION.md` §5.2 **before** the production run wrote its artifact.

**The disposition, and why it is the Family-B-correct one:** the criterion's **purpose** — evidence
that the prior really was neutralized — **is met by the arm itself**, and better than r₁ could ever
show it. Its **literal test** is a proxy the generator never touches. **The instruction was followed
and the truth recorded** — honest quantiles, an explicit
`realized_r1_degenerate_at_ablation_r1 = false` flag, and a `neutralization_mechanism` field naming
`p12_lattice.jl:311`. **Pinning r₁ to make the check pass would have been moving the mechanism to
satisfy the measurement** — the precise inversion of what a check is for, and a thing this phase has
refused three times now under other names.

Two further things must be said about the list together, and neither excuses the other.

**The honest framing.** Every one of the four was caught **before it corrupted a result**, and it was
caught *because* the bar had been written down in advance and could therefore be checked against the
physics. A project that set its thresholds after seeing its numbers would have had none of these four
findings — not because the mistakes would not have happened, but because nothing would have surfaced
them. The pre-registration discipline is what makes these visible; the count is a property of the
audit, not of the science.

**The uncomfortable part, stated without hedging.** Four is enough that "we pre-registered it" no
longer settles an argument on its own, because the pre-registration itself has now been wrong four
times. The distinguishing feature in all four cases was **not** that the number was inconvenient — it
was that the bar could be shown, from evidence independent of the machinery under audit, to be
measuring something other than what it named. In this case that independent evidence is the analytic
information limit; in Phase 11 it was the component/total decomposition. **That standard, not the
inconvenience of the number, is what licensed each amendment**, and it is the standard a referee
should hold the fourth one to.

**Duty on Phase 16:** this belongs in the manuscript's methods as a stated observation about the
pre-registration process, not left for a reader to assemble from four separate phase records. Report
it as **two recurring design errors with four instances**, not as a list of four — the families are
the finding, and they turn "they keep getting bars wrong" into "the pre-registration discipline
surfaced two recurring design errors, both caught before either corrupted a result." Report the
mechanism of each family, the common licensing standard above, and Family A's rule verbatim: *a
control must be measurable when the effect under test is absent, or it is not a control.*

### §7.1 A THIRD family, added 2026-07-31 by 12-15: STRUCTURAL VERIFICATION CANNOT SEE A WRONG-SPACE ERROR

Families A and B are both about **bars**. This one is about **verification itself**, and it is the
generalisable lesson of the 12-15 episode.

12-15's first run produced an artifact that was **complete, internally consistent, and scientifically
void**. It scored standardized posterior draws against raw truth, having never applied
`bundle.theta_zt` (`spike/npe/infer.jl:35-38`). It nonetheless passed **every** check the phase knew
how to make:

- all 40 expected keys present, `schema_version` correct;
- held-out indices disjoint from `train_indices` **and** `val_indices`, all three arms, zero overlap;
- 64 per-region values pooling **exactly** to the reported scalar, every arm;
- the wall-clock budget reconciling to **0.006 min** against its own parts.

**It passed all of them because none of them look at SCALE.** That is the family, stated generally:

> **Every structural property of a result — cardinality, disjointness, aggregation identity,
> conservation, internal reconciliation — is PRESERVED BY AN AFFINE TRANSFORM OF THE VALUES. So no
> amount of structural verification can detect that the numbers are in the wrong space. A structural
> check answers "is this table well-formed?", never "is this table about the right quantity?"**

The failure mode this creates is worse than a loud error: **it produces an artifact that looks
harvestable.** Had it been harvested, `P12_CHOSEN_PRIOR = :car` and `P12_N_LOW = 3` would have entered
an **append-only** pre-registration permanently, and `n_low = 3` — which meant only "the bar equals
the field's own marginal sd because the posterior carries no information in the space it was scored
in" — would have reached the manuscript as a physical finding about biology.

**The cheap fix is a check that DOES look at scale: report the null predictor's score beside your
own.** `trivial_rmse` costs one line. An RMSE reported alone is uninterpretable; reported beside what
a *constant* predictor achieves it is self-checking, because "60 % worse than predicting zero" is not
a number a trained model can produce in the right space. It is now mandatory for every
`spike/validation/run_p12_*.jl` runner that reports an RMSE, enforced by a sweep in
`test_p12_train.jl` rather than by prose — **because prose is exactly what failed here.** The contract
existed, was written down, was honoured by `infer.jl:110,122` and `benchmark.jl:188`, and was violated
silently by the one consumer that did not read it.

**Duty on Phase 16:** report this as the third family. Families A and B are about writing bars down
correctly; this one is about the limits of checking. A referee who is told "we verified the artifact"
should be told *which* class of error that verification could and could not have caught.

### §7.2 A FOURTH family, added 2026-07-31 by the 12-15 descope: ONE NAME, TWO MEASUREMENTS, IN TWO PLACES

Families A and B are about **bars**; §7.1 is about **verification**. This one is about **naming**, and
it is the only entry so far that reached into *routing* rather than into a number.

**`12-CONTEXT.md:310-312` defines D-12 Stage 1 as the CAR-vs-GP mini-spike itself:**

> *"**Stage 1 (cheap, before the training spend):** at the CAR-vs-GP mini-spike, require the better
> prior to show a stated coverage improvement over the neutralized-prior ablation on simulated data.
> **If neither does, descope before the full training run.**"*

**The implemented Stage-1 gate is a different measurement entirely** — 12-11's ridge borrowing probe,
whose `VERDICT: PROCEED` is the only thing `p12_stage1_verdict()` reads. Both are called "D-12
Stage 1". They measure different quantities, at different times, on different data.

**The consequence was not cosmetic.** On 2026-07-31 the mini-spike returned `NONE-BEATS-ABLATION` in
three runs — *precisely the condition `12-CONTEXT` names as the Stage-1 descope trigger* — while the
implemented gate had already returned PROCEED. So the phase entered a state **no plan step routes**:
12-14 Task 3, the only executable producer of the D-13 deliverable, is keyed on `:descope` and did
nothing; 12-17 is impossible as written; and 12-16 is specified to throw. **The descope had to be
authorised by a standalone user ruling** (`12-D13-AUTHORISATION.md`) because no branch existed to
carry it.

> **The family, stated generally: when a pre-registered decision point is given a NAME, and that name
> is later attached to a different measurement, the pre-registration silently acquires two readings —
> and nothing fails loudly, because each reading is individually coherent. The divergence is visible
> only where the two are put side by side.**

**Why it survived four review passes:** each document is internally consistent. `12-CONTEXT` describes
a coherent Stage 1; 12-11 implements a coherent Stage 1; `p12_stage1_verdict` reads a coherent
verdict. **Nothing is wrong within any one file.** This is the same structural blindness §7.1
identifies for scale — a property that no single-document check can see, because it is a property of
the *relation between* documents.

**The near-miss that makes it concrete:** the available shortcut was to write `VERDICT: DESCOPE` into
`12-STAGE1-VERDICT.md`, which would have fired Task 3 and required no new code. **That would have
falsified a gate record — a ruling the user made on independent evidence, which did not change — in
order to obtain a route.** It was refused, and the descope was carried by a new record instead. Recorded
because the cheapness of that shortcut is exactly what makes this family dangerous: the wrong fix is
one line, and it looks like routing maintenance.

**Distinct from Family B**, which is a quantity compared against a bar built for a different purpose.
Here **no bar was misapplied at all** — the correct bar was applied to the correct quantity, and the
*name* of the gate was what carried two meanings.

**Duty on Phase 16:** report as the fourth family, and report the count honestly — this makes **five
root-cause instances across four families**, all surfaced by pre-registration and all caught before
corrupting a result. The mitigation is cheap and generalisable: **when a decision point is renamed,
re-sited, or re-implemented, the document that originally defined it must be amended in the same
commit, or the two definitions must be explicitly cross-referenced.**

### §7.3 THE BUDGET CONSTANTS ARE A FINDING IN THEIR OWN RIGHT, added 2026-07-31

Four separate entries in this record are about **wall-clock ceilings**, and taken singly each reads as
an isolated slip. **Taken together they are one finding about how this project's budget constants were
written**, and Phase 16 should report them as one rather than let a reader assemble four accidents:

| # | what | where |
|---|---|---|
| 1 | **A gating ceiling with NO DERIVATION** — `P12_STAGE1_RATIO_CEILING = 0.5`, on a control whose only route to passing was the effect it was meant to control for | §4 of this document |
| 2 | **TWO ceilings with IDENTICAL VALUES guarding DIFFERENT budgets** — `P12_DATAGEN_WALLCLOCK_CEILING_MIN` and `P12_MINISPIKE_WALLCLOCK_CEILING_MIN`, both 150. A conflation obvious at 150-vs-90 is invisible at 150-vs-150 | §7 Family-B near-miss |
| 3 | **A ceiling with NO ENFORCER** — the mini-spike ceiling had none in the codebase at all until 12-15's runner was required to implement it | §7 Family-B near-miss |
| 4 | **A REAL SPEND WITH NO CEILING AT ALL** — 12-16 simulates and scores at length, and **no Tier-1 constant governs it.** `P12_DATAGEN_WALLCLOCK_CEILING_MIN` fires only inside `generate_p12_pool`, which 12-16 never calls (it extends `harness.jl`'s `draw_simulate_infer`); `P12_MINISPIKE_WALLCLOCK_CEILING_MIN` is explicitly 12-15's three-arm run. **Yet 12-16's artifact schema records `elapsed_min` — a measured cost with nothing to compare it to.** | found 2026-07-31 |

**The pattern: the ceilings were written where a spend was ANTICIPATED, not where spends actually
occur, and their values were chosen before the thing they govern existed.** Entry 4 is the sharpest
form — the schema author knew a cost was worth recording and recorded it, and no one noticed there was
no bar for it to be recorded against.

**Consequence for entry 4, stated so it is not mistaken for compliance:** a 12-16 run **cannot** be
reported as "within budget" or "no breach", because there is no budget. It can only be reported as a
**measurement**. Saying "no breach" would imply a bar was cleared that does not exist.

**Duty on Phase 16:** report as one finding about budget-constant construction, not four slips. The
generalisable rule: **a constant that records a cost must name the bar that cost is judged against, or
it is telemetry rather than a budget — and telemetry must not be reported as compliance.**

### §7.4 A FOURTH FAMILY-B INSTANCE and a SECOND §7.1 instance, added 2026-07-31 by the 12-15 scorer audit

Two entries, from one audit of `_p12ms_score_arm`. **Both are ADDITIONS to existing families, not new
families** — which is itself the point worth reporting: the families are now predicting what turns up.

---

**(a) A SECOND §7.1 WRONG-SPACE DEFECT, IN THE SAME FUNCTION, FOUND BY LOOKING FOR IT.**

§7.1 was written after `_p12ms_score_arm` scored **standardized** draws against **raw** truth. The
audit that followed found a **second, independent** wrong-space error in the same function: the θ
coefficient rows are stored **permuted into `p12_dct_order` smoothness order** (`p12_generate.jl:193`),
and the scorer applied `p12_idct_vec` **directly to the permuted vector**, never inverting the
permutation. Every arm was therefore scored against a **scrambled reconstruction of the truth**.

**The magnitude, measured rather than argued** (real `sample_p12_prior` draw, idx 7, `:car`):

| | max error vs. the true field | RMS error / field's own RMS |
|---|---|---|
| correct inverse (`p12_region_field`) | **3.11e-15** | — |
| the scorer's naive inverse | **3.72** | **1.302** |

> **A PERFECT POSTERIOR WOULD STILL HAVE SCORED ~1.30× THE TRIVIAL PREDICTOR. The three arms were
> ranked on how well each net's coefficients happened to survive a scramble.**

**Why §7.1's own mitigation could not catch it, and this is the sharpening §7.1 needs.** §7.1
prescribed `trivial_rmse` — report the null predictor beside your own. **It fired.** Skill was
negative for every arm in all three runs. What a null baseline establishes is *that* the numbers are
in a wrong space; **it cannot say WHICH wrong space**, so the first defect's fix left a second one
wearing the first one's symptoms. And **a permutation is orthogonal**, so on top of §7.1's affine
argument even *Parseval* holds: total energy is preserved exactly.

> **Sharpened rule for §7.1: a null baseline detects a wrong space; it does not localise one. When a
> wrong-space defect is found, audit every space the read path crosses, not only the one that
> failed.** Applied here, that audit compared generation against scoring on fifteen agreement points —
> the standardization inverse, the permutation inverse, the DCT direction, the θ row ordering, ρ vs
> Fisher-z vs Gaussian space, region indexing, batch shape, and the rest — and executed all fifteen
> rather than reading them.

**What the defect did and did not reach, established by execution.** *Affected:* `rmse`, `skill`,
`rmse_per_region`, `coverage`, `coverage_per_region`, hence **admissibility**, hence the
`NONE-BEATS-ABLATION` verdict. *Unaffected:* `c0_rmse` / `c0_coverage` (because `p12_dct_order(G)[1]
== 1`, verified bit-identical on both sides), `r1_shrinkage`, the truncation curve (flat index space
throughout), and **the trained bundles** (training consumes θ as stored and never reconstructs a
field). Recording which numbers survive is what lets a write-up be precise instead of gesturing.

**A third instance was hunted and NOT found.** The batched `sampleposterior` fallback
(`smp isa AbstractVector ? smp : [smp]`) would have scored 1 dataset per chunk of 100 while dividing
by 1000, had the library returned a single matrix for a batch. It returns a `Vector`, verified by
execution. **"No third defect" is a finding here because it was looked for, not because nothing was
noticed.**

Repaired at `f039729`, routed through the pre-existing `p12_region_field` (`spike/p12/result.jl:562`)
rather than a second inverse — the repository containing **both a right and a wrong reconstruction of
the same object** is the proximate cause and is not reproduced. Guarded by a test that is numeric as
well as textual, because the file legitimately contains `p12_idct_vec` for the truncation curve and
source text alone cannot separate the correct use from the wrong one.

---

**(b) A FOURTH FAMILY-B INSTANCE: `P12_STAGE2_COVERAGE_TOST_DELTA` IS APPLIED TO A PSEUDO-REPLICATED
STATISTIC. NAMED, NOT CORRECTED.**

`P12_STAGE2_COVERAGE_TOST_DELTA = 0.03` is derived in `p12_consts.jl` from a binomial standard error
at **N ≥ 271 DATASETS** (`P12_STAGE2_N_MIN`), and `select_prior` **asserts** that `n_test ≥ 271`
before applying it. But the quantity it is applied to — `coverage` — is the mean over **64 regions ×
1000 datasets**, and the 64 within-dataset trials are **strongly correlated**: they are 64 readings of
one posterior over one field. The effective independent n is therefore neither 1000 nor 64000, and the
bar is sized for neither.

This is Family B exactly: **a sound bar applied to something it was not built for.** It is the same
shape as entry 2 (`effective_independent_n` of 2 against a bar sized for N ≥ 271) and as the SC1g
component-vs-total error — **the third independent recurrence of pseudo-replication in this
milestone**, which is the finding, more than any single instance is.

**IT IS NAMED AND NOT MOVED, AND THE REASON IS THE ONLY ONE THAT LICENSES LEAVING IT.** It is Tier-1
pre-registered and append-only. Correcting a bar *after* seeing that the runs it gated came back
negative — and negative for an unrelated, now-repaired reason — is indistinguishable from tuning a
gate until it passes, which is exactly what this phase refused to do when it declined to hunt the
epoch count that lands coverage inside the band. **A bar that is wrong and pre-registered is reported
as a named limit; it is not quietly repaired by the party it inconveniences.**

**Duty on Phase 16:** report (a) as a second §7.1 instance **with the sharpened rule**, and (b) as the
fourth Family-B instance, carrying the note that it was left standing deliberately. The honest
aggregate sentence: **the families are no longer only a retrospective classification — §7.1 predicted
where to look for (a), and Family B predicted the shape of (b).** That is a stronger claim about the
audit than the count of instances is a weakness in the science.

### §7.5 A FIFTH family, added 2026-08-03: AN ASSERTION WHOSE SUBJECT IS THE **ABSENCE** OF STATE

Families A and B are about **bars**; §7.1 is about **verification**; §7.2 is about **naming**. This one
is about **preconditions**, and it is the first family in this record that can be triggered by a commit
that never touches the test or the code the test covers.

**The instance.** `test_p12_train.jl`'s Stage-1-guard testset distinguishes *permitted by the guard*
from *blocked by the guard* by the message of the error that comes **next**:

```julia
# `arm = :none` is permitted -- it gets PAST the guard and fails later, on the absent
# :none pool, which is a DIFFERENT error.
derr = try train_p12_npe(arm = :none, epochs = 18, n = 10_000, verdict_dir = td); "" catch e ... end
@test !occursin("may not be trained", derr)      # NOT the guard's refusal
@test occursin("no COMPLETE pool", derr)         # the pool check, i.e. past the guard
```

That reasoning is sound **only while the `:none` pool at n = 10 000 does not exist**, and
**nothing in this repository maintains the absence of a pool.** `88ddfe2` (2026-07-31, the D-13
ablation) built exactly that pool. From that commit onward the call sailed past the pool check,
**trained a real net for 2m15.7s** (measured; the estimate carried in the handoff was ~3.5 min, and
the smaller figure is the one that was executed), and then failed on an empty `derr`. Measured on the current
tree, `p12_pool_complete` returns **`true`** for `:none` at n = 10 000 and **`false`** for `:car` — and
the `:car` sibling three lines below was still green for that reason alone, which is what makes the
diagnosis executable rather than argued.

> **The family, stated generally: an assertion whose subject is the ABSENCE of an artifact has a
> precondition that no file, test or convention owns. Every other precondition a test relies on lives
> somewhere that a reviewer can be pointed at — a fixture it builds, a constant it reads, a function it
> calls. Absence lives nowhere. Any commit, in any plan step, can falsify it by CREATING something,
> with no edit to the test and no edit to the code under test — so no diff review, on either side, has
> anything to look at.**

**It is distinct from every family above.** Families A and B misapply a bar; §7.1 measures the right
thing in the wrong space; §7.2 attaches one name to two measurements. Here the bar, the space and the
name are all correct — **the test simply stopped being about its subject**, silently, from a distance.

**What makes it worse than a stale test:** the failure mode is not "passes vacuously". It is
**consumes the phase's compute and then fails for a reason unrelated to anything anyone changed.**
The next person to see it red reads it as a regression in the guard, or in the trainer, or in whatever
commit happened to be in flight.

**Repaired 2026-08-03** by pointing both assertions at a `mktempdir()` path that is never created, so
they test the **guard** rather than the **environment**. `isdir` alone would not have sufficed:
`p12_pool_dir` resolves through `open_or_invalidate`, which **creates** the directory it resolves
(`train_p12_npe.jl:356-359`), so the real `:car` cache path already exists as an empty hash directory
and only `p12_pool_complete` sees through it. The `:car` sibling was repaired in the same move rather
than left standing as a demonstration: **12-17's production pool is specified to build n = 10 000
`:car`**, so it was one planned commit from the identical failure.

**THE COUNTER-EXAMPLE, AND IT IS THE LOAD-BEARING HALF OF THIS ENTRY.** The obvious lesson —
*"avoid absence-dependent assertions, decouple every call from the cache"* — **licenses exactly the
wrong fix**, and the same testset contains the case that proves it.

Two assertions in it are guard-**refusals**: `@test_throws Exception train_p12_npe(arm = :car,
epochs = 18, n = 10_000, verdict_dir = td)`, under `:absent` and under `:descope`. They carry **no
message assertion**, and they currently throw **from the guard**. Handing *them* the absent pool
would mean the `no COMPLETE pool` error satisfies `@test_throws Exception` just as happily as the
guard's refusal does — so **a guard regression would be MASKED rather than caught.** Applying this
family's "fix" uniformly would have converted a live test into a vacuous one **in the same commit that
repaired a vacuous one**, and the suite would have gone greener while testing less.

> **The rule is therefore NOT "decouple every call from the state of the working tree." It is: LET AN
> ASSERTION DEPEND ON STATE ONLY WHERE THAT DEPENDENCE IS WHAT PROVES THE POINT.** For the two "past
> the guard" assertions the pool's absence is *incidental* — they are about the guard, and the pool is
> scenery, so the scenery must be owned by the testset. For the two refusal assertions the real
> cache is *load-bearing* — it is what forces the throw to come from the guard and nothing else, so
> touching it would destroy the discrimination the assertion exists to make.

Which is why the two were treated **oppositely in one commit**: the discriminating question is never
"does this depend on state?" but **"if this dependence were removed, would the assertion still be able
to fail for the reason it names?"**

**THE PROCESS FINDING, which is the more general half.** This phase's gate signal is the **per-file
test run** — the file a plan step touched is the file that gets run. `88ddfe2` had no reason to run
`test_p12_train.jl`, and did not. **Fifteen commits** separate `88ddfe2` from the discovery, and only
one of them (`f039729`) touched `test_p12_train.jl` at all. So a unit test could spend **2m15.7s
training a neural network** on every full-suite invocation, and there were no full-suite invocations
to notice. **The whole file runs in 1m57.1s once repaired** — the defect was costing more than the
entire rest of the file put together.

> **A per-file gate cannot see a cross-file precondition break, because the file that breaks it is
> never the file that is run.** The cheap mitigation is not "always run the full suite" — that is what
> the per-file gate exists to avoid — but: **a test that spends real compute must be unable to spend it
> by accident.** Every assertion in this repository that a call *fails fast* should be pointed at state
> the testset itself owns, so the fast path is guaranteed by construction rather than by the state of
> the working tree.

**A sweep for siblings was run across the whole p12 test suite** (`test_p12_{architecture, consts,
coverage, datagen, decoupling, lattice, prior, result, sbc, suite, train}.jl` and
`capture_p12_golden.jl`). **The two repaired here are the only assertions whose PASS depends on
something not existing.** The suite does contain eight *presence-conditional* branches — the mirror
shape, which skips while a file is absent and arms itself when it lands — and **all eight are currently
armed and live**, every file they condition on having been built. That is a clean sweep, and it is
reported as a result rather than as an absence of one.

### §7.6 A §7.2 INSTANCE, added 2026-08-03: THE ABLATION'S EXEMPTION FROM THE CALIBRATION GATE EXISTS ONLY AS THE RANGE OF A `for`

**This is an ADDITION to §7.2's family (one rule, two readings), not a new family, and it is recorded
ADDITIVELY. Nothing above it is edited.** In particular `12-D13-AUTHORISATION.md` §8 is left exactly
as written — **§8 treating all three arms alike on coverage is the evidence**, and making it agree
with the code retroactively would destroy the finding rather than resolve it.

**THE QUESTION THAT SURFACED IT.** The `NONE-BEATS-ABLATION` verdict disqualifies both spatial arms on
coverage against a band of `P12_COVERAGE_NOMINAL ± P12_STAGE2_COVERAGE_TOST_DELTA` = **[0.87, 0.93]**.
On the repaired control the ablation's own coverage is **0.98694** — also outside that band. So: is
the ablation exempt from the gate, and if so, where is that written?

**THE PROSE SAYS ALL ARMS.** `12-15-PLAN.md:127`, the pre-registration:

> *"1. An arm is **admissible** only if its Gaussian-space per-region coverage is within
> `P12_STAGE2_COVERAGE_TOST_DELTA` of `P12_COVERAGE_NOMINAL`."*

**"An arm"**, unqualified. `select_prior`'s own docstring (`run_p12_minispike.jl:222`) repeats it
unqualified. `12-D13-AUTHORISATION.md` §8 likewise treats the three alike — *"all three arms over-cover
at 18 epochs, all three under-cover at 100"*.

**THE CODE SAYS TWO ARMS.** `run_p12_minispike.jl:249`:

```julia
adm = Symbol[]
for a in (:car, :gp)
    abs(scores[a].coverage - P12_COVERAGE_NOMINAL) <= P12_STAGE2_COVERAGE_TOST_DELTA && push!(adm, a)
end
```

`:none` is **never entered into the admissibility test at all.** It appears in `select_prior` only as
the RMSE bar of rule 3. The `n_test >= P12_STAGE2_N_MIN` assertion immediately above *does* loop over
all three — so the file is not uniformly two-armed; the narrowing is specific to admissibility.

**THE DIRECTION OF THE EFFECT, STATED PLAINLY BECAUSE IT CUTS THE REASSURING WAY. THE LOOP BOUND DID
NOT MANUFACTURE THE VERDICT.** Hold the ablation to the same gate and it fails too (0.98694 against
[0.87, 0.93]). The outcome is then *"no arm is admissible"* — which selects **no spatial prior**
either. **`NONE-BEATS-ABLATION` survives holding the ablation to the gate**, so nobody should read
this as a discovered error in the result. The arithmetic is not in question and is not re-derived
here.

**WHAT IT ACTUALLY IS: A PROVENANCE DEFECT, NOT AN ARITHMETIC ONE.** The ablation's role as the
**unconditional fallback** — the thing you land on when nothing else qualifies, rather than a
candidate that must itself qualify — is a real and defensible design choice. It is simply **implemented
as a loop bound and asserted nowhere.** A load-bearing verdict in this phase rests on a step that
exists only as the range of a `for`.

**AND IT IS UNTESTED, WHICH IS CHECKABLE RATHER THAN ASSERTED.** The eight synthetic score tables in
`test_p12_train.jl`'s `select_prior` branch testset pass `none = (0.90, …)` — coverage **exactly at
nominal — in all eight**. No test ever hands the ablation an out-of-band coverage, so **every one of
those tests would pass identically whether or not `:none` were in the admissibility loop.** The
branch coverage that reads as thorough does not touch this decision at any point.

> **The §7.2 shape, in its cheapest possible form: the divergence is not between two documents here but
> between the pre-registered prose and the three lines implementing it — and it is invisible from
> inside either one. The prose is coherent. The code is coherent. Only the pair disagrees, and no test
> puts them side by side.**

**No change to the loop is proposed and none is made.** Consistent with §7.4(b), a pre-registered rule
is not quietly rewritten by the party who found it inconvenient; it is **named**. What this entry
establishes is that the ablation's exemption is a **choice that was never recorded as one** — and
therefore that any write-up saying the ablation "was not required to be calibrated" is describing a
loop bound, not a decision anyone can be pointed at.

## §8 Routing consequences of PROCEED

- `p12_stage1_verdict()` now returns `:proceed`, so `p12_require_proceed` admits **12-17, 12-18,
  12-19 and 12-20**. All four run.
- 12-13, 12-14, 12-15 and 12-16 were never excluded by this gate and run as written.
- Nothing is descoped. **SPAT-07** stays with 12-18 at full scope (12-16 does **not** record
  `spat07_scope = :reduced_descope`), and **SPAT-08** stays with 12-20 — it is **not** deferred to
  v2.1. §7 of `12-STAGE1-ADJUDICATION-BLOCKED.md` enumerates what a DESCOPE would have deferred; on
  this route none of it applies, and that section is retained only as the counterfactual it was
  written as.

---

*Phase: 12-spatial-colocalization-map*
*Plan: 12-11 (measurement) — adjudicated 2026-07-31 by user ruling*
