# Phase 12 — D-12 Stage 1: measurement complete, ADJUDICATION BLOCKED

**No verdict is recorded by this document, and that is deliberate.** The pre-registered rule
resolves to neither `PROCEED` nor `DESCOPE` on the numbers measured, and its third branch cannot be
written truthfully either. Which of the two remaining readings is correct is a **pre-registration
question**, which the executing agent has no authority to settle.

`12-STAGE1-VERDICT.md` is therefore **NOT created**. That is the load-bearing consequence, not an
omission: `p12_stage1_verdict()` (`p12_consts.jl:637`) returns `:absent` when the file is missing,
and its own docstring defines `:absent` as *"the Stage-1 gate has not been adjudicated yet"* — which
is exactly and literally the state of this phase. `p12_require_proceed` therefore refuses 12-17,
12-18, 12-19 and 12-20, as it should until a human rules.

---

## §0 Provenance

| Item | Value |
|---|---|
| Date | 2026-07-29 |
| `git rev-parse HEAD` at measurement | `e88b97f435daf1113899c5c2bd229da6c3c6720d` |
| Artifact | `spike/validation/p12_stage1_report.jld2` |
| Artifact `generated` | `2026-07-29T19:14:32.274Z` |
| Runner | `spike/validation/run_p12_stage1_ridge.jl` (committed at `e88b97f`, **before** any pool was generated) |
| Pre-registration | `spike/validation/p12_consts.jl` — Tier 1, frozen before any Phase-12 field existed; **byte-unchanged by this plan** |
| Reserved stream | pools ride `P12_DEV_SEED ⊻ P12_DATAGEN_SALT` keyed per global index; the ridge consumes **no** random draw (deterministic tail-block split). `P12_STAGE1_COUNTER = 2` is recorded as the audit tag, not consumed |
| Arm / size | `arm = :car`, image size PINNED to `P12_STAGE1_IMSIZE = (512, 512)`, `imsize_tag = :pinned_512x512` |
| z-scoring arm | `:per_row` (R-3, `P12_ZSCORE_ARM`) |
| Iterations spent | **0 of `P12_ITERATION_ALLOWANCE` = 1** — the allowance is carried forward UNSPENT |

Five pools, one per rung, `r1` **PINNED** and reaching the content hash, in five distinct
content-hash directories under `spike/data/cache/p12/`:

```
0.05  d5562341054cb89d254b3f4afb269856f94df2aad5aae84bce1acd016bdd5120
0.25  026051006ff8e32694d82d0a35fce26c6df59fc8a7a3e618f31fb01a3de42c47
0.50  6e1aa27e7b0fbeee439b9c02585965f5d9e2fbc2316d7622d07a8608863cc3f8
0.75  89cdd7400eea4e3a637f760b60605ae6c11bd92316e8e6ebb1936d531630684c
0.95  18a985ecbdbe1eed165df96ea0e4344281562db546da610e72b49fd7371e8d08
```

Every realized `r1` in every pool was asserted equal to its pin, and every realized image size
equal to `(512, 512)`, before any number below was computed.

## §1 The pre-registered thresholds, quoted from the artifact's own keys

| Artifact key | Value | Tier-1 name |
|---|---|---|
| `ratio_ceiling` | **0.95** | `P12_STAGE1_RATIO_CEILING` |
| `control_ceiling` | **0.5** | `P12_STAGE1_CONTROL_CEILING` |
| `min_rungs` | **1** | `P12_STAGE1_MIN_RUNGS` |

All three were fixed in Tier 1 **before any Phase-12 field existed** — before the lattice kernels,
before the simulator edit, before this runner. All three are members of `P12_GATING_CONSTANTS`.
`p12_consts.jl` is byte-unchanged by this plan and nothing here proposes editing it.

`P12_STAGE1_RATIO_CEILING` carries its derivation verbatim in the file (n_test = 1250 per rung,
sampling sd ≈ 1/√(2·1250) = 0.0200, 2.5 sd of separation from 1.0 ⇒ 0.950), together with its own
recorded caveat that the bar is applied to a 64-region MEAN while the arithmetic is for a
SINGLE-REGION ratio, making it conservative. **`P12_STAGE1_CONTROL_CEILING` carries no derivation
anywhere in the pre-registration.** That is a fact about the file, checked by reading it, and it
matters to the ruling below.

## §2 The ladder, as measured

`n_per_rung = 5000`, realized `n_test = 1250` at every rung — equal to
`P12_STAGE1_N_TEST_PER_RUNG_DERIVED`, asserted at every rung before the ridge ran, so the gate has
the power its ceiling was derived for.

| r₁ | `ratio_mean` | `ratio_min` | `ratio_max` | `control_ratio_mean` | `control_ratio_max` | `rho_arm_ratio_mean` | `radial_r2` | rung vs 0.95 |
|---|---|---|---|---|---|---|---|---|
| 0.05 | 0.92557 | 0.84632 | 1.01020 | 0.73188 | 0.87510 | 0.92439 | 0.34522 | **PASS** |
| 0.25 | 0.84660 | 0.75726 | 0.97470 | 0.67740 | 0.81498 | 0.84641 | 0.37843 | **PASS** |
| 0.50 | 0.71870 | 0.61900 | 0.85285 | 0.59827 | 0.72389 | 0.71985 | 0.25225 | **PASS** |
| 0.75 | 0.56184 | 0.48333 | 0.70133 | 0.49659 | 0.61621 | 0.56276 | 0.41096 | **PASS** |
| 0.95 | 0.33569 | 0.29381 | 0.42932 | 0.32369 | 0.40648 | 0.32311 | 0.74890 | **PASS** |

The gated target is the **Gaussian-space** lattice value `z_field[r]` (R-2, atom-free). The
`rho_arm` column is **REPORTED ONLY**; its ratio is flattered by the elementwise `ghat` clamp's
measured 6.79 % per-region atom mass, and in the event it tracks the Gaussian arm to within 0.002
at four of five rungs, so nothing in the reading depends on which is quoted.

Supporting measurements, same artifact:

| r₁ | `baseline_rmse_mean` | `within_image_sd_rows` | `within_image_sd_rho` | `within_image_sd_zfield` |
|---|---|---|---|---|
| 0.05 | 0.99798 | 0.24666 | 0.46954 | 0.99276 |
| 0.25 | 0.99568 | 0.26154 | 0.46229 | 0.97781 |
| 0.50 | 0.99008 | 0.26783 | 0.43346 | 0.91534 |
| 0.75 | 0.97952 | 0.23320 | 0.33716 | 0.70830 |
| 0.95 | 0.99772 | 0.17013 | 0.15571 | 0.32801 |

## §3 The two gate components

```
borrowing_ok = count(ratio_mean .<= 0.95) = 5  >=  min_rungs = 1     ->  TRUE
control_live = maximum(control_ratio_mean) = 0.73188  <=  0.5        ->  FALSE
stage1_pass  = borrowing_ok && control_live                          ->  FALSE
```

**The half the phase is about passed, at every rung, by a margin.** The half that exists to prove
the instrument works did not clear its ceiling at the three shortest correlation lengths (it does
clear it at r₁ = 0.75 and 0.95).

## §4 Why no verdict line is written

The interpretation rules were fixed in the plan **before the numbers were seen**. The first one
governs here, quoted verbatim from `12-11-PLAN.md`:

> **`control_live` false** (control ratio above `P12_STAGE1_CONTROL_CEILING`) ⇒ the harness is not
> live and the whole run is **uninformative**. Do not record a PROCEED or a DESCOPE. Record
> "harness not live", name the most likely causes in order (a broken standardizer, a
> target/predictor misalignment, a pool with the wrong `r1`), and spend the single
> `P12_ITERATION_ALLOWANCE` on fixing the harness — not on relaxing a threshold.

So the rule forbids both verdicts. It then directs the phase's single, unrepeatable iteration
allowance at **fixing the harness** — and that instruction rests on a premise the rule states
rather than tests. The premise is testable, the three causes it names are individually checkable,
and all three are ruled out below. Recording "harness not live" would therefore put a claim into
this phase's decision record that its own artifact contradicts.

That leaves three writable contents the plan specifies — `PROCEED`, `DESCOPE`, `harness not live` —
and **none of them is true**. That, and only that, is why this document exists under a different
name and carries no verdict.

### The three named causes, each ruled out executably

| Named cause | Evidence against it |
|---|---|
| a broken standardizer | Fitted on the FIT split only, at every rung and for every region, with the analog's zero-variance guard. Verified directly: after masking, region r's continuous and mask rows are **exactly 0.0** raw *and* **exactly 0.0** after standardization, while the unmasked row is not. The control and masked designs differ **only** in those two rows. |
| a target/predictor misalignment | `z_field` and `rho_field` are reshaped `G²×n` column-major, which is the identical `p12_idx` ordering `encode_d01` produces for the summary rows, so summary row r and lattice row r are the same region with no transpose. Independently confirmed by the own-row correlation being high and rising with r₁ (0.635 → 0.881); a misaligned target would produce a correlation near zero at every rung. |
| a pool with the wrong `r1` | Every realized `r1` asserted `== ` its pin at every rung; the five rungs resolved to five distinct content-hash directories, asserted **before** generation. A collapse to one directory was demonstrated to be detectable by the same assertion. |

### The evidence the rule does not ask for, and which decides the question

`own_row_bound = sqrt(1 − corr(row_r, z_field[r])²)` is the RMSE-over-prior ratio of the **best
possible linear predictor** of a region's lattice value from that region's **own** summary row. It
is a property of the DATA — no ridge, no penalty, no split, no standardizer — so it cannot be
produced or distorted by any defect in the machinery it is used to audit.

| r₁ | own-row ridge | its information limit | difference | full-128 control | global-level control |
|---|---|---|---|---|---|
| 0.05 | 0.76588 | 0.76564 | 0.00024 | 0.73188 | 0.77173 |
| 0.25 | 0.71236 | 0.71201 | 0.00035 | 0.67740 | 0.61400 |
| 0.50 | 0.63716 | 0.63667 | 0.00049 | 0.59827 | 0.38141 |
| 0.75 | 0.55590 | 0.55551 | 0.00039 | 0.49659 | 0.23234 |
| 0.95 | 0.47199 | 0.47183 | 0.00016 | 0.32369 | 0.19687 |

Three readings, none of which requires a judgement call:

1. **The ridge attains the analytic information limit to within 5·10⁻⁴ at every rung.** The
   estimator is not merely working, it is optimal against a bound computed without it.
2. **The information limit itself lies above the ceiling at the three shortest rungs**
   (0.766, 0.712, 0.637 against a bar of 0.5). At r₁ ≤ 0.50, **no estimator of any kind** — correct,
   broken, or neural — can bring a per-region ratio to 0.5 from that region's own row, because the
   information is not in the data. The full-128 control beats the own-row limit at every rung by
   borrowing from the other regions, and that is what carries it under the ceiling at r₁ = 0.75.
3. **The Phase-11 benchmark reproduces.** 12-CONTEXT S-2 records global ρ_true recovered from this
   same 128-row summary at ratio **0.157** on a harness known to work. The global-level control here
   measures **0.19687** at r₁ = 0.95 — the rung where the two are comparable, since the global level
   of a near-constant field is the analogue of a constant ρ — at the noisiest single image size in
   the F5 set. It degrades to 0.772 at r₁ = 0.05 for the reason the same physics predicts: with a
   rough field, the image-mean of 64 near-independent regions is a small target against unchanged
   measurement noise.

12-RESEARCH Pitfall 2 predicted this in advance, in these words: *"Per-region recovery is therefore
possible only when the injected field's per-region contrast is comparable to or larger than ~0.066
in correlation units"*, and *"under the D-05 copula every region's marginal is `MU_PRIOR` regardless
of the correlation length, so the within-image contrast is controlled entirely by the correlation
length"*. The measured `within_image_sd_zfield` falls 0.993 → 0.328 across the ladder while the
control ratio falls 0.732 → 0.324, exactly in step.

## §5 The finding this run delivers regardless of the ruling

**Neighbouring regions do carry information about a held-out region, and the amount rises
monotonically with the correlation length.** `ratio_mean` runs 0.926 → 0.847 → 0.719 → 0.562 →
0.336 across r₁ = 0.05 → 0.95, clearing the pre-registered 0.95 ceiling at **all five** rungs, with
`ratio_max` (the worst single region of 64) clearing it at four of five and missing at only the
shortest rung (1.010 at r₁ = 0.05).

**The crossing point, stated in the gate's own units.** Comparing the masked arm (the other 63
regions, no own row) against the own-row information limit — i.e. asking *when do the neighbours
become more informative about a region than that region's own measurement*:

| r₁ | neighbours only (`ratio_mean`) | own row only (limit) | neighbours better? |
|---|---|---|---|
| 0.05 | 0.92557 | 0.76564 | no |
| 0.25 | 0.84660 | 0.71201 | no |
| 0.50 | 0.71870 | 0.63667 | no |
| 0.75 | 0.56184 | 0.55551 | **parity** (Δ = 0.0063) |
| 0.95 | 0.33569 | 0.47183 | **yes**, by a wide margin |

The crossing lies **between r₁ = 0.75 and r₁ = 0.95**, at essentially exactly r₁ = 0.75. Pitfall 2
named this crossing "a genuine scientific deliverable" and predicted its shape — improving
monotonically with r₁ and crossing "somewhere in the middle" — before it was measured. Note the
crossing is quoted here against the **own-row information limit**, a measured quantity in the
artifact, rather than against the research's 0.0662 per-draw noise sd: that figure was measured at
1376×1028 while this ladder is pinned to 512², and Pitfall 2 itself records that the noise sd is
image-size dependent, so the two are not directly comparable and are deliberately not compared.

**The S-4 early warning, recorded as the plan requires.** `radial_r2` — the fraction of the
64-region recoverability spread explained by radius alone — is 0.345 / 0.378 / 0.252 / 0.411 /
**0.749**. At r₁ = 0.95 three quarters of it is radial. This must be carried into 12-20's guard
interpretation **whatever the ruling**, together with its competing benign explanation: at long
correlation lengths the per-region spread is small and dominated by lattice geometry, and 12-07
measured the CAR marginal sd to be largest at low-degree (corner and edge) cells — which is itself
a radial pattern with no chromatic content. The two readings are not separated by this run, and
12-20's ε = 0 and offset-grid arms are what separate them.

## §6 What the user must rule on

The measurement is complete and is not in dispute. The open question is one sentence:

> **Does a positive control that is provably at the information limit of the data, but above
> `P12_STAGE1_CONTROL_CEILING = 0.5` at r₁ ≤ 0.50, mean the harness is not live — or does it mean
> the ceiling is unreachable in that unit at those rungs?**

Both readings are defensible from the record, which is precisely why this is not the executing
agent's call:

**Reading A — the rule stands as written.** The control ceiling is Tier 1, it is in
`P12_GATING_CONSTANTS`, and reinterpreting a gate constant after seeing the number it failed is the
"amend after the fact" pattern this project has already paid for twice. Consequence: the run is
uninformative, and `P12_ITERATION_ALLOWANCE` is spent on the Stage-1 arm.

**Reading B — the ceiling is mis-scaled for the unit it is applied to.** `P12_STAGE1_RATIO_CEILING`
carries a derivation and an explicit recorded warning that *"THE UNIT THE BAR IS APPLIED TO IS NOT
THE UNIT IT WAS DERIVED FOR"*. `P12_STAGE1_CONTROL_CEILING` carries **no derivation at all**; its
only stated rationale, in 12-RESEARCH's own code sketch, is qualitative — *"must be far below 1.0,
else the harness is dead"* — a bar the measured 0.324–0.732 (a 27–68 % error reduction) meets at
every rung. Its plausible implicit anchor is Phase 11's global ρ_true control at 0.157, which
12-CONTEXT **S-2 explicitly warns must not be read as a per-region expectation**: *"It does not
prove per-region recovery."* Consequence: the control is live, and the gate reduces to
`borrowing_ok`, which is TRUE at 5 of 5 rungs.

**Which is which is a pre-registration question, and the executing agent will not answer it.** No
threshold has been altered, no constant appended, no iteration spent, and no verdict written.

If the ruling is Reading B, note that a bar which the pre-registration cannot reach at short
correlation lengths would need to be expressed per rung or against the measured information limit,
and any such change is an append to Tier 2 with provenance naming this artifact — never an edit to
Tier 1, which is append-only by construction.

## §7 If the eventual ruling is DESCOPE, this is what it defers

Recorded here so a later reader does not have to derive it, per the plan. **Nothing below is a
verdict; it is the consequence table for one of the possible rulings.**

- **SPAT-07** (per-region calibration: Gaussian-space + randomized-rank SBC with equivalence
  testing) is claimed by **12-18**, which the descope routing skips. It is carried in **REDUCED**
  scope by **12-16** instead, which records `spat07_scope = :reduced_descope` in its artifact.
- **SPAT-08** (the radial-energy, offset-grid and ε = 0 guards) is claimed by **12-20**, which the
  descope routing skips, and **cannot be carried anywhere else** — 12-06's ε ridge covers only the
  identifiability half, and it measured `ridge_residual_shrinkage = 0.96870` against a floor of
  0.90, i.e. `vacuous = true`, leaving S-4's radial chromatic confound **LIVE**. On a Stage-1
  descope, **SPAT-08 is DEFERRED TO v2.1**, in those words. The D-13 ablation would then ship
  without the S-4 guard suite, which is a real reduction in the honesty apparatus and belongs in the
  memo rather than being discovered later. §5's `radial_r2 = 0.749` at r₁ = 0.95 makes that
  deferral more consequential, not less.
- On the descope route the phase continues with **12-12, 12-14 and 12-16** to build and calibrate
  the D-13 ablation model — a per-region Δρ map **with** per-region uncertainty, deferring only the
  spatial borrowing to v2.1 — which is strictly more than `LocalColocMap` offers today, since that
  carries no uncertainty at all. D-13's own caveat applies honestly: at a **Stage-1** descope the
  ablation exists only at spike scale, so the deliverable is smaller than the Stage-2 version of the
  same fallback.

## §8 What is forbidden until a ruling exists

- **Running 12-17, 12-18, 12-19 or 12-20.** `p12_require_proceed` refuses all four while
  `12-STAGE1-VERDICT.md` is absent, and that refusal is correct — the gate is unadjudicated, not
  passed.
- **Editing any constant in `p12_consts.jl`.** Tier 1 is append-only; a diff that MODIFIES a line is
  by construction a pre-registration breach.
- **Re-running the ladder on a different seed, size or arm to obtain a nicer control.** The reserved
  stream is fixed, the pools are content-hash addressed, and they regenerate byte-identically.
- **Spending `P12_ITERATION_ALLOWANCE`** before the ruling. It stands at **1**, unspent.
- **Reading the banner or `stage1_pass` as a DESCOPE.** `stage1_pass = false` here is an AND over
  one TRUE and one FALSE component, and the frozen rule forbids reading the FALSE half as a descope.

---

*Phase: 12-spatial-colocalization-map*
*Plan: 12-11*
