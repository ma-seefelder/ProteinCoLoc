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

Two families, three root-cause instances in Family B and two in Family A — and every one was caught
by the same mechanism: writing the claim down in advance made it checkable against something outside
itself.

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
