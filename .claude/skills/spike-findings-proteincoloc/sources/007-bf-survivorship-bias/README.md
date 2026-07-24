---
spike: 007
name: bf-survivorship-bias
type: standard
validates: "Given 40% baseline-side attrition, when r̂ is computed on survivors only, then it differs systematically from r̂ over all pairs under a finite-preserving baseline"
verdict: VALIDATED
related: [006, 008]
tags: [bayes-factor, ship-gate, survivorship-bias, statistics]
---

# Spike 007: Survivorship Bias in the BF Correlation

## What This Validates

**Given** the 40% baseline-side attrition established in Spike 006, **when** the gate computes
`corr(amortized, baseline)` on survivors only, **then** that estimate differs systematically from
the correlation over all pairs.

## Research

None needed — pure post-processing. Every clamped-baseline variant is a closed-form function of
the `p_post` values already stored by Spike 006, so **no new inference run was required**.

## How to Run

```bash
julia --project=. .planning/spikes/007-bf-survivorship-bias/quantify_bias.jl
```

Requires `.planning/spikes/006-bf-attrition-mechanism/attrition_rows.jld2`.

## Investigation Trail

1. **Reconstructed clamped baselines** from stored `p_post` at ε ∈ {1e-8, 1e-10, 1e-12, 1e-15}.
   A finite-preserving baseline is what makes "r̂ over ALL pairs" computable at all.
2. **Compared survivors-only vs all-pairs** at each ε.
3. **Wrote a wrong claim and caught it.** The first version of the script asserted the clamped
   dropped-pairs form "a constant" with "zero variance" — while printing `sd = 9.83` beside it.
   They are **two-valued** (+18.41 for the 37 pairs at `p_post = 1`, −18.43 for the 3 at
   `p_post = 0`); the sd reflects that ± split, not within-group spread. Corrected in the script
   and below. The substantive point survives, but weaker than first stated: sign is preserved,
   *magnitude* is destroyed.
4. **Added Spearman** as a tail-mapping-invariant cross-check.

## Results

**VERDICT: VALIDATED.** n = 100 (60 survivors, 40 dropped), `log_prior_odds = 0.00885`.

### The survivor filter depresses r̂ by roughly 0.03–0.05

| clamp ε | r̂ survivors | r̂ all pairs | difference |
|---|---|---|---|
| 1e-8 | 0.9296 | 0.9775 | **+0.0479** |
| 1e-10 | 0.9243 | 0.9728 | +0.0485 |
| 1e-12 | 0.9243 | 0.9661 | +0.0417 |
| 1e-15 | 0.9243 | 0.9539 | +0.0295 |

The gate's survivors-only estimate is **lower** than the all-pairs estimate at every clamp
position. Against a 0.95 threshold this is decision-relevant: the survivors-only value sits below
it, the all-pairs value above.

### What the survivor subsample is

| group | median \|Δρ_true\| | IQR | n |
|---|---|---|---|
| survivors | 0.255 | [0.130, 0.521] | 60 |
| dropped | 0.691 | [0.478, 0.888] | 40 |

The gate validates the amortized estimator **only where the signal is weak**.

### But the all-pairs number is not trustworthy either

The clamp maps every one-sided posterior onto ±log((1−ε)/ε) *regardless of how one-sided it is*:

| ε | distinct values taken by the 40 dropped pairs |
|---|---|
| 1e-8 | {−18.4295, +18.4118} |
| 1e-12 | {−27.6399, +27.6222} |

Sign survives, magnitude does not. Within that group the baseline has no resolution left, so the
"all pairs" correlation is driven by where the arbitrary clamp sits — visible in the table above,
where r̂(all) moves from 0.9775 to 0.9539 purely by changing ε, a quantity with no scientific
content. **Neither number is a sound estimate.**

Spearman (invariant to the monotone tail mapping) is stable at ρ_s ≈ 0.935–0.957 across all ε, but
the saturated pairs are tied, so it too discards the tail information rather than using it.

### Conclusion

The BF gate's correlation is **not** a clean measurement under either treatment:
- **survivors-only** (status quo): biased low, computed on the ambiguous middle, and starved of n;
- **clamped all-pairs**: inflated/deflated by an arbitrary ε with no resolution in the tail.

The tail regime needs a baseline that is finite *and* unsaturated — pursued in Spike 008.

### Figure

Panels A and B of `../008-bf-baseline-remedies/remedies.png`.
