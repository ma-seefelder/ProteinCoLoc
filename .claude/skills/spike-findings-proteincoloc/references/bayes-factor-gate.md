# Reference: The Bayes-Factor Ship-Gate (spikes 006, 007, 008)

**Question chain:** The amended grid-8 confirmatory ship-gate's BF arm returned only **58 of 100**
finite pairs (below the pre-registered `BF_GATE_N_MIN = 69`) → `corr_verdict = :invalid`. These
three spikes diagnose *why*, whether the survivors are trustworthy, and whether any drop-in baseline
fixes it. **Bottom line: the weak link is the KDE baseline, NOT the amortized NRE under test.**

## The diagnosis (what these spikes established)

### 006 — attrition is 100% baseline-side and boundary-triggered
- The gate keeps a pair only if **both** log-BF values are finite (`isfinite(a) && isfinite(b)`,
  `test/gate/run_gate.jl`). On a DEV seed (n=100): **40 dropped, 40/40 baseline-side (KDE), 0
  NRE-side, 0 both.** The amortized NRE **never** produced a non-finite log-BF.
- Mechanism: `_p_gt_threshold_unclamped` (`src/amortized/bf.jl:113-124`) returns
  `clamp(1 - p_le, 0, 1)` with **no `1e-8` floor** (the clamp was removed deliberately, T-7-06 /
  Memo §5, so `max|Δ logBF|` would be free of the saturation artifact). QuadGK roundoff on a
  one-sided posterior lands `p_post` on **exactly** 0.0 or 1.0 → `log(p/(1-p)) = ±Inf` → the pair is
  filtered. Of 40 drops: 37 at `p_post == 1.0` exactly, 3 at `0.0`, **0 strictly inside (0,1)**. The
  boundary is a floating-point event, not a clean statistical one (one-sidedness is *necessary but
  not sufficient*: 31 of the 60 *kept* pairs are also fully one-sided).
- **The drops are the high-signal pairs, not random.** Median `|Δρ_true|`: dropped **0.691** vs kept
  **0.255**. The gate discards precisely where the Bayes factor is easiest/most interesting and
  computes its correlation on the ambiguous middle.
- Attrition **rises with image size** (512²→2048²: 33%→71%): more pixels/patch → sharper posteriors
  → more fully one-sided draws → more boundary hits. The realistic image-size mixture (finding F5)
  *aggravated* an attrition the 256²-only regime had masked.

### 007 — the survivors-only filter biases r̂ down, and the clamped all-pairs number is also unsound
- Every clamped-baseline variant is a **closed-form function of the `p_post` values Spike 006
  already stored** → 007 answered a new question with **zero additional inference**. (This is why the
  convention "store the intermediate quantities, not just the verdict" pays off.)
- The survivor filter depresses r̂ by **~0.03–0.05** across clamp ε ∈ {1e-8…1e-15}:
  survivors r̂ ≈ 0.924–0.930 vs all-pairs r̂ ≈ 0.954–0.978. **Decision-relevant against a 0.95
  threshold** — survivors below it, all-pairs above.
- But the all-pairs number is *not* trustworthy either: the clamp maps every one-sided posterior
  onto ±log((1−ε)/ε) **regardless of how one-sided it is**. The 40 dropped pairs collapse to just
  **two values** (e.g. {−18.43, +18.41} at ε=1e-8; {−27.64, +27.62} at ε=1e-12). **Sign survives,
  magnitude is destroyed.** r̂(all) moves from 0.9775→0.9539 purely by changing ε — a quantity with
  no scientific content. Spearman is stable (~0.94) but ties the saturated pairs, also discarding
  the tail. **Neither treatment is a clean measurement.**
- (Note a caught self-error: the first script called the clamped dropped-pairs "a constant / zero
  variance" while printing `sd=9.83` — they are *two-valued* (±split), not constant. Corrected in
  the committed script. Report refuted sub-claims, don't quietly rewrite them.)

### 008 — the clamp MASKED a genuine order-of-magnitude NRE-vs-KDE tail disagreement (VERDICT: PARTIAL)
- **The hypothesis was WRONG and the refutation is the finding.** Proposed remedy C computes the
  tail in **log space** (see recipe below): 100% finite, 0 saturation, and — validated against the
  unclamped baseline where both are defined — the **same estimator** (median |log-space − unclamped|
  = **0.002** over the 65 finite pairs). Yet the correlation *collapsed* to **0.4448**.
- Why: on the dropped pairs the log-space baseline spans **[−1045, +9524]** while the NRE stays
  bounded near **±8**. Clamping both into ±18.42 made them look concordant; dropping them made the
  disagreement invisible; computing them honestly makes it unmissable. **So the Memo §5 "clamp
  artifact" framing is INCOMPLETE — the clamp is not manufacturing the discrepancy, it is *masking*
  one.**
- **But the KDE baseline can't be trusted at those magnitudes either.** With `L = 999` posterior
  draws the smallest resolvable tail probability is 1/999, so **|log-BF| beyond ~log(999) ≈ 6.9 is
  extrapolation of the Gaussian kernel's assumed tails, not information in the draws.** A log-BF of
  9523 from 999 draws is not a measurement. **In the tail regime the KDE Bayes factor is not a valid
  reference at all** — clamped it saturates, unclamped it diverges, log-space it extrapolates; all
  three are the same fact: the baseline is asked for tail probabilities its estimator cannot resolve
  from 999 draws.

## How to build it (the log-space tail baseline — the reusable artifact)

From `sources/008-bf-baseline-remedies/logspace_baseline.jl`. Both linear-space failure modes
(clamp-saturation and ±Inf) share one root cause: `1 − p` underflows to 0 once `p` is within ~1e-16
of 1. Compute the tail in log space instead — `logcdf(Normal(), z)` is stable far into the tail:

```julia
# Stable log-sum-exp (avoids a LogExpFunctions dep in the fragile root project)
function logsumexp(v)
    m = maximum(v); isfinite(m) || return m
    return m + log(sum(exp.(v .- m)))
end

# Gaussian-KDE tail in LOG space. P(X>t) = (1/n) Σ_i Φ((x_i - t)/h).
function log_tail(draws; threshold::Real = 0.0)
    x = collect(float.(draws)); h = _bw(x)            # _bw = Silverman rule
    z = (x .- threshold) ./ h
    lP   = logsumexp(logcdf.(Normal(),  z)) - log(length(x))
    l1mP = logsumexp(logcdf.(Normal(), -z)) - log(length(x))
    return lP, l1mP
end
logspace_logbf(draws, log_prior_odds; threshold=0.0) =
    (first(log_tail(draws; threshold)) - last(log_tail(draws; threshold))) - log_prior_odds
```

This is diagnostically correct (finite + unsaturated). **It does not rescue the gate** — it exposes
that no draw-based KDE reference can adjudicate the tail. Use it to *measure agreement where the
baseline is resolvable*, not to compare magnitudes past |log-BF| ≈ log(L).

## Build guidance for the gate (PROPOSALS — none applied)

1. **Restrict the BF comparison to the regime where the baseline is resolvable** (pairs whose tail
   probability is resolvable at the given `L`), and report the excluded fraction explicitly rather
   than silently dropping it.
2. **Raise `L`** if deeper tails must be compared — but the ceiling scales as log(L): resolving
   |log-BF| ≈ 18 needs L ≈ 6.6e7 draws (impractical). This bounds what any draw-based KDE reference
   can do.
3. **Validate magnitude agreement where measurable; compare sign/rank beyond that.** Demanding
   magnitude agreement in an unresolvable regime is testing the reference, not the estimator.
4. **Reconsider gating `max|Δ logBF|` (BF_LOGBF_TOL = 0.5) on a tail the reference cannot resolve.**

## What to avoid

- **Don't read attrition as an NRE failure.** It is 100% baseline-side; the NRE never went
  non-finite. The estimator under test is not the thing breaking.
- **Don't trust the clamped all-pairs r̂** — it is driven by the arbitrary ε position, which has no
  scientific content.
- **Don't compare BF magnitudes past |log-BF| ≈ log(L) ≈ 6.9** (L=999) — that is Gaussian-tail
  extrapolation, not data.
- **Don't treat the clamp as the root cause** — it masks a real NRE-vs-KDE disagreement; removing it
  starves n and biases survivors, both worse.
- **Don't apply any gate change here.** Changing the gate is a **pre-registration matter** under
  `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md` §6 — not a spike
  decision.

## Constraints / honesty

- **Verdicts:** 006 VALIDATED, 007 VALIDATED, **008 PARTIAL** (mechanism answered decisively; the
  proposed remedy *failed* — and that failure is the finding).
- **Single DEV seed** for each (006 `0x00B7A771`, 008 `0x00B7A772`), all asserted disjoint from
  `PROD_SEED_V2`/`VAL_MASTER_SEED`/`NPE_MASTER_SEED` and prior dev seeds. n=100. Different seeds mean
  006's 40% and the gate's 42% are the same *mechanism/magnitude*, not the same draws.
- Bandwidth used (Silverman) approximates `KernelDensity.jl`'s default — immaterial in the normal
  regime (0.002 agreement), not tuned for the tail.
- These spikes do **not** establish the NRE is *correct* in the tail — only that the KDE reference
  cannot adjudicate there.

## Origin

`sources/006-bf-attrition-mechanism/` (README + `measure_attrition.jl` + `figures.jl` +
`attrition.png`), `sources/007-bf-survivorship-bias/` (README + `quantify_bias.jl`),
`sources/008-bf-baseline-remedies/` (README + `logspace_baseline.jl` + `figures.jl` +
`remedies.png`).
