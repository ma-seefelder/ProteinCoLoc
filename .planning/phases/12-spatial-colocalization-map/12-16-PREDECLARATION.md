# 12-16 — PRE-DECLARED INTERPRETATION OF THE LRO SCORING RUN

**COMMITTED BEFORE THE SCORING RUN, IN ITS OWN COMMIT. The git timestamp of this commit precedes the
artifact, and that ordering is the point: a reading fixed after the numbers exist is not a
pre-declaration.**

Scope authorised by the user: **SPAT-05** (`p12_coloc_map`) and **SPAT-06** (leave-region-out
predictive coverage in Fisher-z), with `P12_FISHERZ_NEFF` appended. **SPAT-07 is deferred to v2.1.**

---

## 1. THE TRAP THIS DOCUMENT EXISTS TO CLOSE

**Coverage is a property of interval WIDTH as much as of centre.** A predictive interval that is wide
and centred near the marginal mean will cover ≈ 90 % of observations **while carrying no information
about which region is which.**

So **coverage alone cannot separate _calibrated and informative_ from _calibrated because
uninformative_** — and that second branch is **live**, not hypothetical: the D-13 ablation's measured
skill on latent-field recovery is **−0.15084**, i.e. worse than predicting a constant zero.

**Predictive coverage and latent-field skill are DIFFERENT QUANTITIES in different spaces.** A model
can be honest about its own uncertainty while its point map is useless. Deciding what that means
*after* seeing it would be indefensible, so it is decided here.

## 2. THE DISCRIMINATOR — already in 12-16's contract, not invented here

**`prior_only_floor`** (`12-16-PLAN.md:240-242`): predict the held-out region **from its marginal
prior alone** — draws from `MU_PRIOR` pushed through `ghat`, with **no conditioning on the image at
all**. The plan's own words: *"the floor that lets a reader separate nuisance/global borrowing from
spatial borrowing."*

**That is exactly the separation this pre-declaration needs, and it already exists.** No new statistic
is minted.

**The margin is `P12_STAGE2_LOGSCORE_MIN = 0.02` nats per held-out region — the EXISTING
pre-registered constant, reused rather than invented**, the same discipline as reusing
`select_prior`'s 1 % tie-break in the 12-15 re-run. **It is not a new bar.**

## 3. THE BRANCHES, FIXED BEFORE ANY NUMBER EXISTS

### 3.1 SUPPORTS SPAT-05 / SPAT-06 — **both must hold**
1. Pooled LRO coverage within **`P12_COVERAGE_NOMINAL ± P12_STAGE2_COVERAGE_TOST_DELTA` = [0.87, 0.93]**; **and**
2. mean LRO log score **beats `prior_only_floor` by ≥ `P12_STAGE2_LOGSCORE_MIN` = 0.02 nats/region**.

Only then does the phase claim a **calibrated and informative** per-region map.

### 3.2 CALIBRATED BUT UNINFORMATIVE — coverage in band, floor NOT cleared by 0.02

**Reported, verbatim, in these words:**

> **The uncertainty quantification is honest; the map carries no information beyond the prior.**

**THIS IS NOT A PASS, AND IT MAY NOT BE WRITTEN UP AS "PROMISING."** It is the direct analogue of the
12-15 re-run's §2.3(a) — *the method works and the spatial prior adds nothing* — one layer out. It is
the most confusing outcome available and the one most open to being softened, which is precisely why
its wording is fixed here.

### 3.3 MISCALIBRATED — coverage outside [0.87, 0.93]
**SPAT-06 is not met.** Reported as such. **No tuning of epochs, N, thresholds or the noise model to
move it into band** — that would be tuning a model until it passes its own calibration gate, and a
gate that has been tuned to has stopped measuring anything.

## 4. THE EFFECTIVE INDEPENDENT n — DECLARED BEFORE THE NUMBERS, NOT AFTER

**`N = P12_STAGE2_N_MIN = 271` datasets. The independent unit is the DATASET.**

At 64 regions per dataset the run produces **271 × 64 = 17,344 region-draws — but only 271 independent
units.** The 64 regions of one image **share every nuisance and the global term**, so they are
**PSEUDO-REPLICATES**, exactly as `12-19-PLAN.md:148-150` states for the real-image arm (*"regions
within an image share every nuisance and the global term"*).

**MANDATORY, AND THE REASON IS THAT THIS MILESTONE HAS ALREADY MADE THIS EXACT ERROR ONCE** — the
real-image arm, where a gate sized for N ≥ 271 was applied to an arm whose `effective_independent_n`
is **2** (§7 Family B, entry 2):

- **The coverage interval is computed from n = 271, NEVER from 17,344.** An interval built on 17,344
  would be ≈ 8× too tight and would manufacture significance out of pseudo-replication.
- **`effective_independent_n = 271` is a mandatory artifact key**, alongside
  `n_region_draws = 17344`. **One key for the independent unit**, per `12-19-PLAN.md:292`.
- **The 17,344 figure is never reported alone.** Wherever it appears it carries the sentence that the
  region-draws are pseudo-replicates.

**CITATION CORRECTION, recorded rather than propagated:** `12-02-PLAN.md:130` attributes the
independent-unit registration to `12-01:172-178`. **That location is the r₁ prior block and does not
contain it.** The rule is stated at `12-19-PLAN.md:148-150` and `12-02-PLAN.md:235`, which is what this
document cites. The requirement is unaffected; only the pointer was stale. (`12-02-SUMMARY` already
records "three stale line citations" as a defect found in that plan — this is a fourth.)

## 5. MANDATORY REPORTING RULES — §7.1 APPLIED TO ITSELF

**A number whose null baseline is absent is uninterpretable.** This phase has been bitten by that
twice (the wrong-space RMSE that `trivial_rmse` would have caught; the coverage that needed a width).
Therefore, always and without exception:

- **Mean interval WIDTH is reported beside coverage.** Coverage without width cannot be read.
- **`prior_only_floor` is reported beside the log score.** A log score without its floor cannot be read.
- **`crps` is reported alongside** as the second, bounded scoring rule.
- **Clamped regions are counted and reported** (`clamped_region_count`) — a silently clamped ±1
  correlation is a degenerate region, not a measurement.

## 6. WHAT MAY NOT BE CONCLUDED, IN EITHER DIRECTION

**The ablation's latent-field skill of −0.15084 PREDICTS NEITHER OUTCOME ABOVE.** Latent-field
recovery (Gaussian-space DCT coefficients against the drawn field) and LRO predictive calibration
(observed per-region correlation in Fisher-z) are different quantities measured in different spaces.

- A **good** coverage result does **not** overturn, soften, or "rescue" the 12-15 negative.
- A **bad** coverage result is **not** additional evidence for it.

**Neither result may be read as confirming or contradicting the mini-spike.**

## 7. `lro_pass` IS BUILT BUT CANNOT FIRE ON THIS ROUTE

`lro_pass` compares a spatial arm against the ablation. **There is no spatial arm**, so
`stage2_gate = :not_applicable_descope` **by construction**, and the gate is **UNEXERCISED**.

It is built anyway — it is small next to `p12_coloc_map`, it keeps the artifact contract intact for
12-20 Guard 3 (`headline_logscore_delta`), and it is the first thing a v2.1 spatial arm would need.
**Its tests exercise the BRANCH LOGIC, not the gate.** A green test must not be read as the gate having
run, and the test names say so.

## 8. WHAT DOES NOT DRIFT

- **The qualifier "at mini-spike scale" stays on the 12-15 result.**
- > **A comparison between two failed arms is not evidence about priors — the CAR-vs-GP question is
  > UNRESOLVED, not answered.**
- **The disqualification at 100 epochs was a CALIBRATION disqualification, not an accuracy one** —
  both spatial arms beat the ablation on RMSE (12.95 %, 4.01 %) and were rejected solely on coverage.
- **The withdrawn claim stays withdrawn:** "no arm learned the field at all" does not survive a
  re-seed (`:car` reached **+0.08629** at the same 18 epochs).
- **No Tier-2 append except `P12_FISHERZ_NEFF`.** `P12_CHOSEN_PRIOR` and `P12_N_LOW` stay closed.
- **The D-13 deliverable is SPIKE SCALE** (10,000 pairs), and any claim from it inherits that limit.
- **ONE run at N = 271.** Not a search. `P12_ITERATION_ALLOWANCE` is **SPENT**.
