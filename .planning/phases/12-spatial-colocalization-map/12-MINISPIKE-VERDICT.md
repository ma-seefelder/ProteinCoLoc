# 12-MINISPIKE-VERDICT — the D-03 prior choice and the D-04 production rank

Every number below is **re-derived by loading `spike/validation/p12_minispike_report.jld2`**, never
copied from stdout. Generated `2026-07-31T12:57Z`, `schema_version = 1`, Julia 1.12.6.

---

SELECTED PRIOR: NONE-BEATS-ABLATION

PRODUCTION RANK: K = 24

**Both lines are reported, and NEITHER is appended to Tier 2.** On this branch the plan appends no
constant and removes no sentinel. `P12_CHOSEN_PRIOR`, `P12_K_PROD`, `P12_D_PROD` and `P12_N_LOW` do
not exist. The `PRODUCTION RANK` line records what the truncation arithmetic *recommends*; it is not
a constant, because the arm it would be scaled against did not qualify.

---

## 1. The result in one paragraph

At matched induced lag-1 correlation, through a head that cannot under-fit either arm
(`K = P12_K_DEV = 63`, full rank on the deviation space), **neither the CAR nor the GP lattice prior
beats a neutralized prior at mini-spike scale.** The best *admissible* arm — GP, pooled per-region
RMSE **1.06141** — does not beat the neutralized ablation at **1.05793**. Rule 3 therefore returns
`NONE-BEATS-ABLATION`, and under the rule fixed in code before any number existed, the spatial
premise is **unsupported at mini-spike scale**.

**The stronger finding is not the ordering — it is that no arm learned the field at all.** All three
arms score *worse* than predicting a constant zero field:

| arm | pooled per-region RMSE | trivial (predict 0) | **skill** | coverage | admissible |
|---|---|---|---|---|---|
| `:car`  | 1.43483 | 0.99430 | **−0.44305** | 0.99981 | **no** |
| `:gp`   | 1.06141 | 0.99331 | **−0.06856** | 0.89666 | yes |
| `:none` | 1.05793 | 0.99634 | **−0.06181** | 0.95872 | — (ablation) |

`skill = 1 − rmse / trivial_rmse`. Negative means the posterior mean is a *worse* predictor of the
drawn field than a constant zero. The best net is **6.2 % worse than predicting nothing.**

---

## 2. Which branch of the selection rule fired

The rule, fixed in `select_prior` before the run and persisted into the artifact as `selection_rule`:

1. **Admissible** only if Gaussian-space per-region coverage is within
   `P12_STAGE2_COVERAGE_TOST_DELTA = 0.03` of `P12_COVERAGE_NOMINAL = 0.90` — band **[0.87, 0.93]** —
   asserted at `n_test = 1000 ≥ P12_STAGE2_N_MIN = 271`. A miscalibrated arm is DISQUALIFIED, not
   "better".
   - `:car` coverage **0.99981**, |Δ| = 0.09981 → **inadmissible**
   - `:gp` coverage **0.89666**, |Δ| = 0.00334 → **admissible**
   - `admissible = [:gp]`
2. Among admissible arms, the lower pooled per-region RMSE → only one candidate, `:gp`.
3. **The winner must beat the neutralized `:none` arm by any margin.** `1.06141 < 1.05793` is
   **false**. → `NONE-BEATS-ABLATION`. **THIS IS THE BRANCH THAT FIRED.**
4. Tie-break (within 1 %, resolves to CAR): never reached — rule 3 returned first.

Recorded `selection_reason`: *"the best admissible arm :gp (RMSE 1.06141) does NOT beat the
neutralized ablation (RMSE 1.05793). The spatial premise is UNSUPPORTED at mini-spike scale
(rule 3)."*

### 2.1 The margin that decided rule 3 is 0.33 %, and it is NOT a stable quantity

`gp − none = 0.00348`, i.e. **0.329 %** of the ablation's RMSE. That margin must not be read as a
measurement of prior quality, for a reason found while harvesting this run:

**TRAINING IS NOT SEEDED.** `build_p12_estimator` sets no RNG and `NeuralEstimators.train` is called
without one (`train_p12_npe.jl:465-469`), so Flux's weight initialisation draws from the unseeded
global RNG. Only the masking streams (`p12_rng(P12_MINISPIKE_COUNTER)`) and the deterministic
tail-block train/val split are reproducible. A re-run would produce different nets, and a 0.33 % gap
is well inside that variation.

**The verdict does not rest on that margin.** `NONE-BEATS-ABLATION` would hold under any reordering
of `:gp` and `:none`, because the finding that carries it is that **every arm has negative skill**.
A spatial arm that edged the ablation by 0.3 % would still be worse than predicting a constant. What
is *not* established by this run is any ranking *among* the three arms.

This reproducibility gap is a deviation from the project's stated constraint that everything be
reproducible from a fixed Random123 seed. It is `train_p12_npe`'s surface (12-14, already committed),
not 12-15's, and is **recorded here, not fixed here**.

---

### 2.2 The training traces — the run hit its EPOCH BUDGET, not convergence

Recorded here because **the bundle does not persist them.** `train_p12_npe.jl:472-500` stores the
estimator, both transforms, `theta_rows`, `D`, `epochs`, `batchsize`, the index sets and the mask
rate — **no per-epoch risk** — and `NeuralEstimators.train` returns only `est`. These numbers survive
only because this run's stdout was redirected to a log. A training run whose convergence cannot be
audited afterwards cannot support a claim about convergence, so persisting the trace is a recorded
gap in the same class as the unseeded init.

| arm | initial val | ep 1 train / val | ep 9 train / val | ep 18 train / val | last Δval | **val/train @18** |
|---|---|---|---|---|---|---|
| `:car` | 852,488.1 | 481157.53 / 156.202 | 123.628 / 123.079 | 116.786 / 117.628 | −0.019 (falling) | **1.007** |
| `:gp` | 1,306,758.8 | 255.372 / 134.495 | 99.528 / 100.727 | 91.054 / 94.340 | **+0.029 (rose)** | **1.036** |
| `:none` | 97,665.2 | 312.342 / 133.862 | 109.367 / 111.160 | 101.903 / 105.751 | −0.076 (falling) | **1.038** |

**Early stopping never fired. All three arms ran exactly 18 of 18 — the run ended on its epoch
budget, not on convergence.**

**THE RATIO IS THE HARDER CLAIM, AND IT IS THE ONE THAT MATTERS.** The slope says only "still
falling". The **val/train ratio of 1.007–1.038 says these nets are nowhere near capacity** — they are
not fitting their own training data well enough to separate from validation at all. "Cannot learn this
from the summary" has a characteristic and *different* shape: training risk falling well below
validation while validation stalls, the model memorising what it cannot generalise. **None of that is
present.**

Combined with the learning rate annealing **33× to 1.51E-05 across exactly the 18-epoch budget**, this
is a textbook under-trained run rather than an information ceiling.

**But the experiment as designed CONFOUNDS "converged" with "annealed to a stop"**, because the LR
schedule is tied to the epoch count. The traces alone therefore cannot separate *under-trained* from
*cannot-learn*; the no-overfitting fact is what tips it, and that is where this run's evidence ends.

## 3. Global term, shrinkage, and the truncation curve

| arm | `c0_rmse` | `c0_coverage` | `r1_shrinkage` | `vacuous` |
|---|---|---|---|---|
| `:car`  | 5.20592 | 1.000 | 2.96625 | true |
| `:gp`   | 1.23340 | 0.995 | 1.48239 | true |
| `:none` | 0.79834 | 0.986 | 1.38125 | true |

The `c0` split shows the arm difference is **not** confined to the spatial part: CAR is worst on the
*global* term too (5.206 against the ablation's 0.798), and its coverage of 1.000 with an RMSE of
5.206 is the signature of an interval so wide it always contains the truth.

`r1_shrinkage = post_sd / prior_sd` against `P12_VACUOUS_SHRINKAGE_FLOOR = 0.90`: **all three arms
vacuous**, and all three **greater than 1** — the r₁ posterior is *wider* than its prior. Reporting-only
(F3). **This is NOT 12-13's `ridge_residual_shrinkage`**: that is a linear point predictor's residual
sd, this is a posterior WIDTH. Read as a pair under `12-13:131-152`; never averaged.

**Truncation curve** (a property of the field ensemble, computed with no net — `p12_dct_order` on the
drawn `z_field`):

| K | 3 | 8 | 15 | 24 | 35 | 48 | 63 |
|---|---|---|---|---|---|---|---|
| truncation RMSE | 0.76443 | 0.68626 | 0.60974 | 0.52662 | 0.43031 | 0.29864 | 6.92e-16 |

Read against `winner_rmse = 1.06141` (the minimum of the two spatial arms, since none qualified):
`k_prod` = smallest K with truncation RMSE ≤ 0.5 × 1.06141 = 0.53071 → **K = 24**.
`n_low` = smallest K with truncation RMSE ≤ 1.0 × 1.06141 → **K = 3**.

**These are recommendations only and were NOT appended.** The bar they are read against is a
posterior RMSE from a net with negative skill, so `n_low = 3` does not mean "the field is smooth" —
it means the bar is larger than the field's own marginal sd (0.994), which the very first grid point
clears trivially. Appending `n_low = 3` as an identified/vacuous boundary would have written a
property of an untrained net into the pre-registration as though it were a property of the biology.

---

## 4. Budgets — TWO SEPARATE BUDGETS THAT HAPPEN TO CARRY THE SAME NUMBER

Reported separately, each against its own constant. Never summed, never crossed.

- **`elapsed_min` = 10.667** — THIS RUN'S TOTAL, against
  `P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150`. **Within budget**; the `BLOCKER` guard did not fire.
- **`datagen_min` PER ARM** — each pool judged **independently** inside `generate_p12_pool` against
  `P12_DATAGEN_WALLCLOCK_CEILING_MIN = 150`: `:car` 0.084, `:gp` 0.000, `:none` 0.000. All three
  near zero because the pools were already complete and resume-by-skip cost nothing. **Their sum is
  not a quantity anything compares to 150.**

**No contention inflation detected, measured rather than assumed.** Train+score this run = 10.574 min
against the first run's 10.156 min — a ratio of **1.041**. The 4 % is not attributed to a cause; it is
simply too small to need one.

---

## 5. What this branch means for the phase

Stated plainly, as the plan requires on this branch:

- **The spatial premise is unsupported at mini-spike scale.** Neither lattice prior beats a
  neutralized one, and no arm beats a constant predictor.
- **The D-13 ablation is now the likely deliverable.**
- **12-17 MUST NOT BE STARTED until the user has seen this.** It is not merely inadvisable — it is
  *impossible as written*: `12-17:137` generates the 50,000-pair pool with `arm = P12_CHOSEN_PRIOR`,
  and that constant does not exist on this branch. Guessing it would burn ~73 min into a directory
  keyed on a guess.
- Waves 10-12 (12-16, 12-18/12-19, 12-20) consume a trained spatial bundle that this run does not
  license.

This is a **reportable outcome, not a failure to resolve** — recorded under a rule written before the
numbers existed, exactly so that it could not be resolved by picking a winner anyway.

---

## 6. Provenance — why this artifact and not the first one

**The 2026-07-31T09:51Z run of this same file was harvest-ready and scientifically void, and was
discarded unharvested.** It scored STANDARDIZED posterior draws against RAW truth, because it never
applied `bundle.theta_zt` — the contract at `spike/npe/infer.jl:35-38`. It passed every structural
check (all 40 keys, index disjointness, per-region pooling, budget reconciliation to 0.006 min)
because **none of those look at scale**.

It selected **CAR** and reported `beats_ablation = true`, `K_PROD = 35`, `N_LOW = 3`.

**The corrected run reverses the substance of every one of those.** The winner does not merely
change — the branch changes, and the two spatial arms swap admissibility:

| | invalid (09:51Z) | corrected (12:57Z) |
|---|---|---|
| selected | **CAR** | **NONE-BEATS-ABLATION** |
| `beats_ablation` | true | **false** |
| admissible | `[:car]` | `[:gp]` |
| `:car` coverage | 0.92730 (in band) | 0.99981 (**out of band**) |
| `:gp` coverage | 0.99983 (out of band) | 0.89666 (**in band**) |
| `k_prod` | 35 | 24 |

Had the first artifact been harvested, `P12_CHOSEN_PRIOR = :car` and `P12_N_LOW = 3` would have been
written into an **append-only** pre-registration, the phase would have proceeded to spend ~73 minutes
generating a 50,000-pair CAR pool, and the phase's conclusion would have been the opposite of what its
own data supports.

The truncation curve is **byte-identical between the two runs** — it is arithmetic on the drawn field
ensemble and touches no net. That is the one quantity the defect could not reach, and it is the
control that confirms the diagnosis rather than a coincidence.

**No partial output from the corrected run was read before the artifact was written.** After a first
run that looked complete and was void, that is part of the provenance and not merely good manners.
