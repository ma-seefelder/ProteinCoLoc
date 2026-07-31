---
phase: 12-spatial-colocalization-map
plan: 16
status: complete
requirements_delivered: [SPAT-05, SPAT-06]
requirements_deferred: [SPAT-07, SPAT-08]
commits:
  - 914dd78  feat(12-16) p12_coverage.jl — the D-09 construction and p12_coloc_map
  - 4ee0bb5  feat(12-16) the LRO run at N=271, and the P12_FISHERZ_NEFF Tier-2 append
  - 1c80c5b  test(12-16) test_p12_coverage.jl
---

# 12-16 — SUMMARY

**Scope built: SPAT-05 (`p12_coloc_map`) and SPAT-06 (leave-region-out predictive coverage in
Fisher-z), with `P12_FISHERZ_NEFF` appended. SPAT-07 and SPAT-08 are DEFERRED TO v2.1 by the user
ruling of 2026-07-31.** One run at `P12_STAGE2_N_MIN` = 271, 19.99 min, exit 0.

---

## 1. THE HEADLINE — the pre-declared branch that fired

**`predeclared_branch = :miscalibrated`. SPAT-06 IS NOT MET.**

The reading was fixed in `12-16-PREDECLARATION.md`, committed in its own commit (`df9ef91`,
22:56:50) **before this run existed**, and the runner **computes** the branch rather than leaving a
reader to apply it afterwards. Branch 3.3, verbatim from the artifact:

> SPAT-06 IS NOT MET: pooled leave-region-out coverage 0.9707 lies outside [0.87, 0.93].
> **REPORTED, NOT TUNED** — moving epochs, N, thresholds or the noise model to bring it into band
> would be tuning a model until it passes its own calibration gate, and a gate that has been tuned
> to has stopped measuring anything.

**No epoch, N, threshold or noise-model change was made or attempted in response to this number.**
`P12_ITERATION_ALLOWANCE` was already SPENT (`12-D13-AUTHORISATION.md` §6). This is not a boundary
call: the interval on the independent unit is `[0.9675, 0.9738]`, clear of the band by ≈ 0.037.

---

## 2. THE PER-ARM TABLE — width beside coverage, floor beside the log score, always

| arm | coverage | **mean width** | log score (nats/region) | CRPS |
|---|---|---|---|---|
| `ablation` — *no spatial borrowing, full nuisance and global borrowing* | **0.9707** | **1.3673** | **−0.6455** | 0.1669 |
| `prior_only_floor` — genuinely prior-only, no conditioning on the image at all | 0.9999 | 1.7412 | −0.8375 | 0.1590 |
| `spatial` | **ABSENT BY DESCOPE**, not missing by error | — | — | — |

- **Discriminator (ablation − floor): `+0.19195` nats/region**, against the pre-registered
  `P12_STAGE2_LOGSCORE_MIN` = 0.02.
- `clamped_region_count` = **0**, `masked_region_count` = **0**. No region was degenerate; every one
  of the 17,344 region-draws was a measurement.
- `logscore_by_imsize`: −0.7350 (512²) · −0.6133 (1024²) · −0.5852 (1376×1028) · −0.4767 (2048²).

**The log score is a density in Fisher-z units.** The change-of-variables Jacobian is omitted
deliberately and visibly: it depends only on the observation, so it is identical across arms at a
given held-out region and cancels exactly in every difference. The 0.02 bar reads a difference, so
it is unaffected — but the absolute figure must not be quoted as a correlation-space density. CRPS
is reported alongside precisely because it **is** in correlation units.

---

## 3. THE COMPARISON THIS PLAN EXISTS TO MAKE — predictive versus parameter coverage, in words

| | value |
|---|---|
| **predictive** coverage (observed correlation entry, noise model included) | **0.9707** |
| **parameter** coverage (drawn `z_field`, NO observation noise, atom-free Gaussian space) | **0.8203** |

**The two miscalibrations point in OPPOSITE directions, and only simulation could show it.**

- The posterior is **over-confident about the latent field**: parameter coverage 0.820 against a
  nominal 0.90, i.e. its Gaussian-space intervals are too narrow. Per-region range 0.697 – 0.897,
  median 0.832 — **no region reaches nominal**, so this is systematic, not a tail effect.
- The observation-noise convolution then **over-widens**: predictive coverage 0.971 against 0.90.
- The noise term more than compensates the posterior's narrowness, and the sum lands **above**
  nominal rather than at it.

The plan anticipated the case where these two cancel *to* nominal — "predictive coverage near
nominal while parameter coverage is far from nominal means the observation-noise model is absorbing
a miscalibrated posterior, which is precisely the failure D-09 cannot detect on a real image."
**What happened is the same mechanism, over-corrected rather than exactly cancelling.** The
substantive point is unchanged and is now measured rather than hypothesised: **on a real image only
the predictive number is visible, and here it moves away from nominal for a reason that has nothing
to do with the posterior being too wide.** A reader given only `coverage = 0.971` would conclude
"intervals too conservative"; the truth is "posterior too narrow, noise model too generous". That
cross-check is the deliverable of this plan, and it earned its place.

---

## 4. WHAT MAY NOT BE CONCLUDED

- **The map IS informative relative to the prior, and that does NOT rescue it.** The ablation beats
  `prior_only_floor` by 0.192 nats/region, twenty times the 0.02 bar — so the pre-declaration's
  §3.2 *"calibrated because uninformative"* branch is **not** what happened, and the phrase must not
  be used for this result. But **calibration is the QUALIFYING condition**: a miscalibrated arm is
  DISQUALIFIED, not "better", and branch 3.3 fires on coverage alone. Reporting the informativeness
  as a mitigation would be exactly the softening §3.2's wording was fixed to prevent, applied to a
  different branch.
- **§6 of the pre-declaration stands.** Predictive calibration and latent-field skill are different
  quantities in different spaces. **This result neither confirms nor contradicts the 12-15
  mini-spike in either direction**, and nothing here rescues or worsens the −0.15084.
- **`stage2_gate = :not_applicable_descope`, and the gate is UNEXERCISED.** `lro_pass` compares a
  spatial arm against the ablation; there is no spatial arm. A "fail" would be read as evidence
  against the spatial prior, which would be wrong. **`lro_pass` is built and unit-tested on
  synthetic inputs — a green testset is NOT the gate having run**, and the testset's name says so.
- **`headline_logscore_delta = nothing`** by construction: it is the SPATIAL-minus-ablation quantity
  that 12-20's Guard 3 reads, and 12-20 does not run on this route.
  `discriminator_logscore_delta` is a **different quantity** and must not be substituted for it.
- **Spike scale.** The bundle is the D-13 fallback at `P12_MINISPIKE_N` = 10,000 pairs, not the
  50,000-pair version 12-17 would have built. Every claim inherits that limit.

---

## 5. THE INDEPENDENT UNIT — 271, never 17,344

`effective_independent_n = 271`, `n_region_draws = 17344`. **The 64 regions of one image share every
nuisance and the global term, so the region-draws are PSEUDO-REPLICATES.**

| interval | value | status |
|---|---|---|
| **cluster, n = 271 datasets** (between-dataset sd) | **[0.9675, 0.9738]**, half-width 0.0032 | **THE ONE THAT MAY BE QUOTED** |
| Wilson at n = 271 | [0.9484, 0.9833] | reported |
| Wilson at n = 17,344 | [0.9685, 0.9727] | **INVALID — recorded only to show the tightening it manufactures** |

All three land outside `[0.87, 0.93]`, so the verdict does not depend on the choice — but the choice
was made before the numbers, and the invalid one is labelled as such in the artifact itself.

---

## 6. THE NOISE MODEL — `P12_FISHERZ_NEFF`, and why the pooled value was NOT applied

Calibrated by holding a parameter FIXED and re-observing it: `n_theta = 5` prior draws per image
size, each observed `n_obs = 20` times through independent simulations. Redrawing the field each
time would have measured the prior convolved with the noise and fitted an `n_eff` far too small.

| image size | `n_eff` | realized `Var(z)` | between-θ sd |
|---|---|---|---|
| (512, 512) | **36.84** | 0.029547 | 0.003113 |
| (1024, 1024) | **128.31** | 0.007980 | 0.000408 |
| (1376, 1028) | **166.43** | 0.006119 | 0.000781 |
| (2048, 2048) | **454.93** | 0.002213 | 0.000225 |
| **pooled** | **90.22** | 0.011465 | — |

**`varz_fit_residual` = 0.018082 — LARGER than the pooled variance itself (0.011465).** One constant
does not fit these four sizes. The run therefore applied the **per-size** value and records
`n_eff_applied = :per_imsize`; the pooled scalar is appended and reported but **not applied**.

The image size is **observed at read time, not latent**, so conditioning on it is free and is
available on a real image too. Applying the pooled constant would have made every 512² interval too
narrow and every 2048² interval too wide, and the pooled coverage could then have sat at nominal as
the **average of two opposite miscalibrations while neither size was calibrated**. Per-size
calibration is `12-16-PLAN.md:168`'s own specification, for Pitfall 2's reason; the measurement
shows how much it mattered.

**`P12_FISHERZ_NEFF` is a DECLARED MODELLING ASSUMPTION, not a threshold.** It parametrizes the
noise model whose JOINT with the posterior D-09 validates. It gates nothing by itself, no branch
compares anything against it, and it joins neither `P12_GATING_CONSTANTS` nor
`P12_REPORTING_ONLY_CONSTANTS`.

### The Tier-2 append, and the paired removal

`spike/validation/p12_consts.jl` diffstat: **123 additions, 1 deletion** — and that single deletion
**is** the reserved `@assert !isdefined(@__MODULE__, :P12_FISHERZ_NEFF)`. Its `test_p12_consts.jl`
mirror was removed **in the same commit** (`4ee0bb5`) and replaced with literal positive assertions;
dropping only one of the three is the documented way this goes wrong.

**`:P12_CHOSEN_PRIOR` and `:P12_N_LOW` STAY CLOSED.** 12-15's NONE-BEATS-ABLATION branch appended
neither, so their absence is a positive record of which branch fired, not an unfinished task.

---

## 7. SPAT-07 AND SPAT-08 — DEFERRED TO v2.1, and what is thereby NOT delivered

Deliberate decisions of **2026-07-31** (`12-D13-AUTHORISATION.md` §7 and §9), recorded in the
artifact as `spat07_scope = spat08_scope = :deferred_to_v2_1` and in prose here:

- **No Gaussian-space SBC, no randomized-rank ρ-space SBC, no nuisance-appropriate equivalence
  testing.** **THE PHASE MAKES NO PER-REGION CALIBRATION CLAIM BEYOND PREDICTIVE COVERAGE.** What
  SPAT-06 addresses is whether the leave-region-out predictive intervals cover what was actually
  measured; it does **not** establish SBC calibration, and the two must not be conflated in the
  memo. *(On this run the predictive claim was not met either.)*
- **No radial-energy, offset-grid or ε = 0 guard.** The phase makes **NO radial-confound claim**;
  the S-4 warning carried from Stage 1 is left standing and unaddressed, and the D-13 ablation ships
  without the S-4 guard suite — `12-11-PLAN.md:258` calls that *"a real reduction in the honesty
  apparatus"*.

---

## 8. NAMED LIMITS OF THIS RESULT

1. **The evaluation joint is the model's own.** Datasets were drawn at `arm = :none`, the same joint
   the D-13 bundle was trained on (F5: train-joint == eval-joint). These numbers are a calibration
   statement **under the model's own generative assumptions**; they say nothing about behaviour under
   a spatially correlated truth, which is a misspecification question and is not asked here.
   Recorded in the artifact as `eval_arm` + `eval_arm_note`.
2. **The log score is in Fisher-z units** (§2).
3. **`suspect_mask_ood = true`, in its SINGLE-ARM form only.** Pitfall 4's check is that *both* arms
   far off nominal implicates the MCAR masking augmentation rather than the prior. There is only one
   model arm here, so that form **cannot be evaluated**. The rates are recorded beside the flag:
   training `realized_mask_rate` = 0.0627 against D-09's held-out configuration of exactly one
   masked region in 64 = 0.0156. **A reader should look at that pair before the prior** — the net saw
   ~4× the masking rate D-09 presents it with, and the augmentation is a live candidate explanation
   for the over-wide predictive that no arm on this route can rule in or out.
4. **Spike scale**, 10,000 pairs (§4).
5. **`p12_coloc_map` was exercised on synthetic images, not on real specimens.** The `--real` path is
   wired and deliberately UNEXERCISED; it is 12-19's, and 12-19 does not run on this route.

---

## 9. DEFECTS FOUND AND HOW THEY WERE HANDLED

**A wrong-basis reconstruction in `run_p12_minispike.jl:326`, reported not fixed.** θ rows are
`p12_dct_vec(z)[p12_dct_order(G)]` — the DCT coefficients **permuted into smoothness order** — and
`p12_dct_order(8)` is not the identity (`[1, 2, 9, 10, 3, 17, 11, 18, …]`). The 12-15 scorer applies
`p12_idct_vec` **directly to a θ column, without inverting that permutation**. Verified executably:
the round-trip is exact to `3.1e-15` through `p12_region_field` and wrong by `3.17` without the
inverse. Because a permutation is orthogonal this is not a small perturbation — it scores the
estimate against a scrambled truth. It touches the 12-15 latent-field skill numbers (−0.15084,
+0.08629) and the Gaussian-space per-region coverage; it does **not** touch `c0_rmse`/`c0_coverage`
(order[1] == 1), nor the trained bundle (training never reconstructs), so `p12_ablation_none.jld2`
and its sha256 are unaffected and this plan's scope was safe. **Reported immediately and NOT acted
on**: it is outside 12-16's scope, and §6 of the pre-declaration already rules that the −0.15084
predicts neither of this plan's outcomes. It does mean the pre-declaration's §1 sentence *"the
uninformative branch is live because latent-field skill is −0.15084"* rests on a number computed in
the wrong basis — a matter for the orchestrator and the user, not for this plan. **This file's own
reconstruction uses `p12_region_field` throughout.**

**A plan defect, reported not silently improved.** `12-16-PLAN.md:201` requires `p12_coloc_map` to
take `N`, `nominal` **and** `n_eff`, but the function reports posterior sd per region and builds no
predictive interval, so `nominal` and `n_eff` have nothing to consume. The signature is kept (it
stops a map and its coverage numbers being built under different settings) and both are validated,
but the docstring says plainly that no predictive interval is stored in the returned object.

**Three defects in this plan's own first draft of the tests**, found by running them: a fixture
sitting exactly on a floating-point boundary (`0.90 - 0.03` = `0.8699999999999999`, so an arm meant
to be *calibrated at the edge* was silently miscalibrated and the testset reported `:disqualified`
where it meant `:fail`); **a line-wrap that broke an acceptance literal** — 12-02's recorded defect,
met again; and a test asserting something the construction never produces (two independent draw sets
from one posterior do not difference to exactly zero). All three fixed at the mechanism, none by
loosening a check. Full detail in `1c80c5b`.

---

## 10. VERIFICATION

| check | result |
|---|---|
| Task 1 plan command | exit 0; Fisher-z round-trips to 1e-12; all three `isdefined` true |
| Task 2 plan command | exit 0 — `not_applicable_descope descope_ablation_only 271` |
| `test_p12_coverage.jl` | **175 pass / 0 fail / 0 error**, 29–35 s (budget 120 s) |
| `test_p12_consts.jl` | **208 pass / 0 fail**, exit 0 |
| `grep -c P12_PENDING_SCAFFOLD` | **0** |
| `p12_consts.jl` diff | **+123 / −1**, the one deletion being the reserved assertion |
| mutation: `lro_pass` margin 0.02 → 0.0 | testset 2 goes **3 RED**; restored byte-exact (sha256 verified) |
| mutation: drop the clamped-±1 degeneracy test | testset 4 goes **3 RED**; restored byte-exact (sha256 verified) |
| `git diff --quiet HEAD -- src spike/Project.toml spike/Manifest.toml corpus` | **clean** |
| `test/test_images`, `spike/data/cache/p11/`, `artifacts/amended_v2` | untouched |
| bundle `sha256` | `4bdf5b03…f017` **re-verified before scoring**; a mismatch refuses to score |

**A positive property worth claiming, because a referee will ask:** these coverage numbers **cannot**
have been computed on training data, structurally rather than incidentally. This plan never opens a
cached pool at all — it draws fresh through `p12_rng(P12_COVERAGE_COUNTER)` on the **validation**
salt, while pools ride the **datagen** salt, and the two are asserted disjoint at load time in both
`p12_consts.jl` and `p12_generate.jl`. The answer is not "we were careful"; it is **"that path does
not exist."**
