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

### 2.1 Per image size — INFORMATIVE, NOT GATED, and the tolerance does not apply at these n

| image size | datasets | coverage | log score | **is [0.87, 0.93] applicable?** |
|---|---|---|---|---|
| (512, 512) | **113** | 0.9638 | −0.7350 | **no** |
| (1024, 1024) | **76** | 0.9714 | −0.6133 | **no** |
| (1376, 1028) | **57** | 0.9742 | −0.5852 | **no** |
| (2048, 2048) | **25** | 0.9913 | −0.4767 | **no** |
| **pooled** | **271** | **0.9707** | **−0.6455** | **yes — this is the gated quantity** |

**`P12_STAGE2_COVERAGE_TOST_DELTA` = 0.03 was derived for N ≥ 271 independent datasets. It does NOT
apply to any row above the last one.** At n = 25 the Wald half-width at p = 0.9 is ≈ 0.099 — more
than three times the tolerance — so a per-size figure printed beside a [0.87, 0.93] band would be
read against a bar it is nowhere near powered for. **This is the pseudo-replication family in a
fourth guise: a bar derived for N applied to a fraction of N.** The pre-declaration's branches apply
to the **pooled** number over 271 and to nothing else.

*(Correcting an assumption in the dispatch: the sizes are **not** ~68 each. `P12_IMSIZE_WEIGHTS` is
an uneven mixture and the realized split was 113 / 76 / 57 / 25 — so the smallest cell is 25
datasets, not 68, and the mis-powering is worse than a quarter-of-N estimate suggests.)*

The per-size numbers are nonetheless worth reading for **direction**: coverage rises monotonically
with image size, 0.9638 → 0.9913, while the log score improves. That is consistent with the
over-widening being strongest where the noise model has the least to do — but with n = 25 in the top
cell it is a pattern, not a measurement, and it is recorded as one.

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

| image size | `n_eff` | realized `Var(z)` | **between-θ sd** | **rel.** | **genuine θ component** |
|---|---|---|---|---|---|
| (512, 512) | **36.84** | 0.029547 | 0.003113 | 10.5 % | **9.7 %** |
| (1024, 1024) | **128.31** | 0.007980 | 0.000408 | 5.1 % | 3.1 % *(barely above the floor)* |
| (1376, 1028) | **166.43** | 0.006119 | 0.000781 | 12.8 % | **12.1 %** |
| (2048, 2048) | **454.93** | 0.002213 | 0.000225 | 10.1 % | **9.3 %** |
| **pooled** | **90.22** | 0.011465 | — | — | — |

**The between-θ spread is reported beside `P12_FISHERZ_NEFF` wherever that constant appears — in the
artifact, in the Tier-2 block, and here — and never alone.** Same rule as width-beside-coverage and
floor-beside-log-score.

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

### 6.1 `n_eff` IS A SUMMARY OF A DISTRIBUTION, NOT A CONSTANT — and why θ-replication was added

The θ-replication was **not** added for robustness. It was added because **the between-θ spread is a
direct test of the assumption that justifies the transform.** Fisher-z is adopted precisely because
`Var(z)` should be approximately independent of the underlying correlation; measuring that spread
measures the premise. **The premise does not hold exactly** — genuine θ-dependence of ≈ 9–12 % at
three of four sizes. That is a measurement of an assumption, not a variance reduction, and the very
first two blocks already showed the premise is not exact.

**The consequence that must be written down: a predictive interval built from a single `n_eff` is
too wide for some θ and too narrow for others, so POOLED coverage can sit inside [0.87, 0.93] while
θ-CONDITIONAL coverage does not — and this run does not measure the latter.** That is the same shape
as the trap the pre-declaration's §1 names: a number that looks calibrated for a reason other than
the one claimed. *(It did not arise here — pooled coverage is outside the band anyway — but it is a
live hazard for any future run that lands inside it.)*

**Two sources of variation, and they are different kinds of thing. We conditioned on what is
observable and named what is not:**

| source | status | what was done |
|---|---|---|
| **image size** | **OBSERVED** at read time, on simulated and real data alike | **conditioned away** — `applied = :per_imsize` |
| **θ** | **LATENT** — cannot be conditioned on, on any data | **IRREDUCIBLE limit**, named, not fixable |

**What `n_θ = 5` bounds, so five numbers are not read as a distribution.** Five blocks characterise a
spread **loosely**: the relative uncertainty on an sd estimate from 5 draws is `1/√(2·4)` = **35 %**.
These figures establish that θ-dependence **exists** and give its **order of magnitude**; they do not
pin it, and no interval should be built from them.

**And part of the observed spread is sampling noise, not θ.** Each block's pooled `Var(z)` is itself
estimated from `n_obs` = 20 observations over 64 regions, a relative sampling sd of
`√(2/19)/√64` = **4.06 %**. The "genuine θ component" column subtracts that floor in quadrature. **At
(1024, 1024) the observed 5.1 % is barely above the 4.06 % floor, so θ-dependence is only weakly
evidenced there** — recorded rather than smoothed over, because quoting the raw 5.1 % as if it were
all θ would *overstate* the very limit this discloses.

**A figure in circulation that should not be quoted: the ~30 % spread** cited during the build came
from the **throwaway smoke run** at `n_θ = 1, n_obs = 4`, where the variance estimate's own sampling
error dominates completely (floor 10.2 % per block, two blocks). **It is not a measurement.** The
reported run gives ≈ 10 %, and that is the number.

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
   training `realized_mask_rate` = 0.0627 against D-09's held-out configuration of exactly one masked
   region in 64 = 0.0156. **But see §8.1 — the mask mismatch predicts the WRONG SIGN, and my earlier
   framing of it as the lead candidate was wrong.**
4. **Spike scale**, 10,000 pairs (§4).
5. **`p12_coloc_map` was exercised on synthetic images, not on real specimens.** The `--real` path is
   wired and deliberately UNEXERCISED; it is 12-19's, and 12-19 does not run on this route.
6. **`P12_FISHERZ_NEFF` IS A CENTRAL ESTIMATE OVER θ, NOT A CONSTANT, AND THE RESIDUAL θ-DEPENDENCE
   IS IRREDUCIBLE.** Genuine between-θ variation of ≈ 9–12 % in `Var(z)` at three of four image
   sizes. **θ is latent, so it cannot be conditioned on** — unlike image size, which is observed and
   *was* conditioned away. **Every coverage figure in this report is therefore a θ-MARGINAL
   statement; θ-conditional coverage is not measured.** Characterised loosely at `n_θ = 5` (±35 % on
   each spread). Full treatment in §6.1.
7. **Per-size coverage is informative, not gated** (§2.1). The pre-registered 0.03 tolerance was
   derived for N ≥ 271 and does not apply at n = 113 / 76 / 57 / 25.

### 8.1 WHICH ERROR DOES THE MASK MISMATCH ACTUALLY PREDICT? — it explains NEITHER, and the sign is why

**A reader meeting `suspect_mask_ood = true` beside two opposite-signed errors will assume it covers
both. It does not cover either, and this section exists so that assumption cannot be made.**

The two errors require **opposite corrections to the same posterior**:

| to fix | the posterior must | effect on the other |
|---|---|---|
| parameter coverage 0.820 → 0.90 | **WIDEN** | pushes predictive further ABOVE 0.90 |
| predictive coverage 0.971 → 0.90 | total must **NARROW** | narrowing the posterior pushes parameter coverage further BELOW 0.90 |

**So no single-signed effect on the posterior can explain both.** They are reconciled only by *the
posterior being too narrow AND the observation-noise term being too wide, with the noise error the
larger of the two.*

**The mask mismatch acts on the posterior, and the sign it predicts is the wrong one.** Fewer masked
regions means **more** information, so the correct posterior at k = 1 is **narrower** than at the
training mean of k ≈ 4. A net that under-conditions on the mask would therefore carry the wider,
average-case posterior into a k = 1 read and be **too WIDE**. The posterior measured here is **too
NARROW**. The hypothesis predicts the opposite of the observation.

**And the mismatch is weaker than "4×" makes it sound.** `P12_MASK_K_SET` is `0:8`, so **k = 1 is
in-distribution, not out of it** — it is one of nine equiprobable values, seen in roughly 1/9 of
training samples, ≈ 1,100 of the 10,000. Under-represented relative to the mean, **not unseen**.
Pitfall 4's concern is covariate shift into a configuration the net has essentially never met; that
is not this configuration. The 0.0627-vs-0.0156 comparison is a comparison of a **mean** against a
**point**, and the point lies inside the training support.

**I over-stated this in my own earlier reporting** — "look at that pair before the prior" implied the
augmentation is the lead candidate. **On the sign analysis it is not**, and a partial explanation
named as partial is worth more than one that sounds complete.

**What the arithmetic does point at, with the right sign: the observation-noise model.** Since the
posterior contributes too *little* width, the excess must come from the noise term, i.e. `n_eff`
fitted too small. There is a structural reason it would be, and it is checkable:

> **`calibrate_neff` measures the MARGINAL observation noise at fixed θ; the LRO predictive needs the
> noise CONDITIONAL on having observed the other 63 regions of the same image.** Any component of
> observation variability that is **shared across regions of one image** — the per-image Otsu
> threshold, the background floor, any global illumination or gain realization — is already pinned
> down by those 63 observed regions, and is then **added a second time** at the held-out region.
> That is a double-count, and it over-widens the predictive **by construction**.

**This is a hypothesis with a predicted sign, not a finding.** It is consistent with both
observations (it inflates the predictive without touching the posterior), but it is **not measured
here**, and it does not by itself explain the over-confident posterior either — that remains
**unaccounted for**.

**What would settle it, stated so it is not re-derived later:** decompose the fixed-θ variance of
`fisherz(observed)` into a **shared** (across-region, within-simulation) component and an
**independent** one, and refit `n_eff` on the independent part alone. If the shared component is
material, the corrected `n_eff` rises, the predictive narrows, and the parameter-coverage shortfall
is left standing as a separate, genuine posterior defect. **That is a v2.1 measurement — it is new
compute on a spent iteration allowance, and it would change a modelling constant after seeing a
result, so it must not be done inside this phase.**

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
*(Subsequently repaired outside this plan, in `f039729`, and recorded as a §7.4 family instance in
`914d6e1`. Both landed after this plan's own work and neither is 12-16's; noted here only so a later
reader does not take the defect to be still live. **The 12-15 numbers it affects have not been
re-run**, so anything quoting −0.15084 or +0.08629 still carries the wrong-basis caveat.)*

**A plan defect, reported not silently improved.** `12-16-PLAN.md:201` requires `p12_coloc_map` to
take `N`, `nominal` **and** `n_eff`, but the function reports posterior sd per region and builds no
predictive interval, so `nominal` and `n_eff` have nothing to consume. The signature is kept (it
stops a map and its coverage numbers being built under different settings) and both are validated,
but the docstring says plainly that no predictive interval is stored in the returned object.

**A criterion met by a construction that dominates it, recorded so it is not "restored" to satisfy a
grep.** `12-16-PLAN.md`'s acceptance criterion asks that a source read show **two `sampleposterior`
calls** in `p12_coloc_map`. The file does contain exactly two — but in `p12_coloc_map` they arrive as
**one `_p12cov_single_stack` helper called twice**, not as two inlined copies. That is strictly
stronger for what the criterion protects: as one helper called twice the sample and control passes
**cannot** differ, whereas two inlined copies could drift apart. Testset 8 asserts the stronger
property directly (exactly two `_p12cov_single_stack(` call sites in `p12_coloc_map`, exactly one
`sampleposterior(` inside the helper). **Inlining two copies to make a literal grep pass would be a
regression, not a fix.**

**Three defects in this plan's own first draft of the tests**, found by **running** them rather than
reading them. Full detail in `1c80c5b`; two of the three generalise:

1. **A fixture sat exactly on a floating-point boundary.** `0.90 - 0.03` = `0.8699999999999999`, so
   `|coverage − nominal|` = `0.030000000000000027 > 0.03` and an arm meant to be *calibrated, at the
   edge* was silently **miscalibrated** — the testset reported `:disqualified` where it meant
   `:fail`, and stayed green because the surrounding assertions were consistent with the wrong
   branch. **The general rule: a fixture for a BRANCH must sit clearly inside it, never on its
   boundary — a test that lands on a threshold is testing floating-point representation, not
   behaviour.** Use half the tolerance. Boundary behaviour may be tested, but then it is a *boundary*
   test, named as one, with its expected value derived the way the implementation derives it. This is
   the phase's recurring "a gate that measures something other than what it names" failure, at test
   scale and harder to see because the suite stays green. *(Recorded here rather than in
   `CONVENTIONS.md`: the orchestrator's authorisation for a conventions entry was scoped to the
   raw/stripped rule below and to nothing else, and I kept it to that.)*
2. **A line-wrap broke an acceptance literal** — 12-02's recorded defect, reproduced independently.
   Now **`CONVENTIONS.md` C-05**, authorised by the orchestrator and citing both instances: **wording
   rules are checked on the RAW source, structural bans on the COMMENT-STRIPPED source**, because a
   rule that lives in prose is deleted by stripping while a ban a comment can satisfy is no ban at
   all. Corollary: a required literal must never be line-wrapped.
3. **A test asserting something the construction never produces.** Two independent draw sets from one
   posterior do not difference to exactly zero, so scoring an image against itself does not yield a
   measured Δρ of 0.0. Replaced by a deterministic construction through 12-12's own
   `p12_region_maps`, which is sharper: region 7 and region 8 carry the **same value** in
   `region_delta_rho` and are unambiguously different states. **The distinguisher is the flag, not
   the value** — which is the property a downstream consumer actually depends on.

All three were fixed at the mechanism; none by loosening a check.

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

---

## 11. THE PATH-AGREEMENT CHECKLIST — every point where a READ path must match the GENERATION path

Requested by the orchestrator after this plan found the **second** wrong-space defect in the same
scoring function (the missing `reconstruct` in the afternoon, the missing inverse permutation at
night). Two independent instances of one family in one function make a third likelier than not. This
is the checklist **this plan's own read path was built against**, and it is the one the repair agent
can audit the old scorer with. **The last column is the point:** several of these are *implemented*
correctly in some files without being *stated* anywhere as a contract — and an ambiguity that two
files resolved differently is exactly how both defects happened.

| # | pairing | generation side | read side | documented as a CONTRACT? |
|---|---|---|---|---|
| 1 | **θ standardization** | pool stores RAW θ; trainer fits `theta_zt`, trains on `transform(theta_zt, θ)` | `reconstruct(bundle.theta_zt, ·)` **before any row is touched** | **YES, loudly** — `infer.jl:35-38`, and mechanically swept on every `run_p12_*` runner |
| 2 | **summary standardization** | pool stores RAW `summary_min` (128 rows, never z-scored at rest) | `standardize_p12(Zraw, bundle.zt)` — z-scores rows 1:64, passes 65:128 through | **YES** — `train_p12_npe.jl:250-261`, `summary.jl:85-87` |
| 3 | **mask → standardize ORDER** | augment on RAW, then standardize (`train_p12_npe.jl:428-445`) | `mask_regions` on RAW, then `standardize_p12` | **YES, loudly** — `train_p12_npe.jl:30-42`; asserted at the producer (12-14 T3) and now at the consumer (12-16 T6) |
| 4 | **the SMOOTHNESS PERMUTATION** | `p12_dct_vec(vec(z_field))[p12_dct_order(G)]` — θ rows are **permuted** | must **scatter back through `p12_dct_order` BEFORE** `p12_idct_vec` — i.e. `p12_region_field` | **NO — AND THIS IS THE DEFECT'S ROOT.** See below. |
| 5 | **DCT direction** | `p12_dct` = `C F C'` | `p12_idct` = `C' Ĉ C`, orthonormal, exactly invertible | **YES** — asserted, `p12_idct∘p12_dct == id` |
| 6 | **flat index convention** | `encode_d01` `vec`s column-major; lattice uses `p12_idx(i,j) = (j-1)G+i` | `reshape(rows[1:G²], G, G)` with **NO transpose**; same `p12_idx` | **YES** — `p12_architecture.jl:31-40`, asserted by `test_p12_architecture` T1 |
| 7 | **Gaussian vs ρ space** | θ carries the **Gaussian** field `z`; ρ is `ghat(quantile(MU_PRIOR, Φ(z)))` | truth scored in **Gaussian** space (atom-free, R-2); observed entries are **correlations** | **PARTLY** — the *distinction* is documented; the *composition* was not named. See below. |
| 8 | **θ row layout** | one assembler, `p12_theta_column` | `P12_THETA_ROWS` / `p12_theta_index` / `p12_row_*`, never inline arithmetic | **YES** — single assembler by design |
| 9 | **head width / truncation** | `p12_truncate_theta`, provably a no-op at `K_dev = 63` | `load_p12_npe` **refuses** a bundle whose `theta_rows` disagrees with its `D` | **YES** — refusal, not repair |
| 10 | **image-size joint (F5)** | pool drawn from `P12_IMSIZE_SET`/`WEIGHTS` | evaluation drawn from the **same** mixture | **YES** — the binding F5 invariant |
| 11 | **prior arm** | pool at `arm = :none` | evaluation at `arm = :none` (§8 limit 1) | **YES** — but only as prose; nothing refuses a mismatch |
| 12 | **held-out region encoding** | naturally unusable patch ⇒ value 0, mask 0 | `mask_regions` produces **byte-identically** that | **YES** — `p12_architecture.jl:129-141` |

### The two that were ambiguous rather than merely unimplemented — say so, as asked

**(4) THE PERMUTATION HAS NO STATED CONTRACT, AND ONE DOCSTRING ACTIVELY MISLEADS.**
`p12_idct_vec`'s own docstring reads *"Inverse of `p12_dct_vec`, on the same column-major flat
layout."* That is **true of `p12_dct_vec`'s output and false of a θ column** — and a θ column *looks
exactly like* a flat coefficient vector: same length `G²`, same element type, same plausible
magnitudes. Nothing at either site says **"θ rows are permuted; do not feed them to `p12_idct_vec`
directly."** `p12_theta_column` documents that it permutes, and `p12_region_field` documents that it
inverts, but **neither names the other**, so a reader holding only one of them reaches the wrong
conclusion with no warning. That is the whole mechanism of the defect. **The cheapest durable fix is
one sentence in `p12_idct_vec`'s docstring** — "θ columns are in `p12_dct_order`; use
`p12_region_field`, not this function" — and it is the repair agent's to make, not mine.

**(7) THE READ-TIME COPULA WAS IMPLEMENTED THREE TIMES INLINE AND NAMED NOWHERE.**
`ghat(quantile(MU_PRIOR, cdf(Normal(), ·)))` appears in `rho_field` (elementwise, over a field), in
`theta_scalar_view` (over the field **mean**), and in `p12_architecture.jl`'s comment for the derived
scalar (over `c0/G`). Three call sites, three different arguments, one composition, **no shared
function** — so nothing enforces that a fourth consumer composes it the same way, and the difference
between "apply to each region" and "apply to the mean" is invisible at the call site. I added
**`p12_z_to_rho`** in `p12_coverage.jl` for exactly this reason and used it everywhere in this
plan's read path; it is documented there as *not a new map* but as the existing composition, named.

**(11) is a real gap with no defect attached yet.** The evaluation arm must match the training arm
for any calibration claim to mean anything, and **nothing refuses a mismatch** — a bundle trained at
`:none` scored against `:car` draws would run happily and produce a misspecification measurement
labelled as a calibration one. This run records `eval_arm` in the artifact and states the limit, but
that is prose. **A cheap executable guard would be to compare `bundle.arm` against the draw arm and
refuse.** Named here rather than built, because adding a refusal to a shared surface is outside
12-16's scope.

**What this plan's own code does about all twelve:** every row is honoured; #1, #3, #4, #6 and #7 are
additionally asserted or exercised in `test_p12_coverage.jl`; #4 and #7 are honoured through the
single named functions `p12_region_field` and `p12_z_to_rho` rather than by inline re-derivation, so
a future reader of this file cannot resolve either ambiguity the wrong way by accident.

---

**A positive property worth claiming, because a referee will ask:** these coverage numbers **cannot**
have been computed on training data, structurally rather than incidentally. This plan never opens a
cached pool at all — it draws fresh through `p12_rng(P12_COVERAGE_COUNTER)` on the **validation**
salt, while pools ride the **datagen** salt, and the two are asserted disjoint at load time in both
`p12_consts.jl` and `p12_generate.jl`. The answer is not "we were careful"; it is **"that path does
not exist."**
