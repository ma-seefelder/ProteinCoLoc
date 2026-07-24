---
spike: 008
name: bf-baseline-remedies
type: comparison
validates: "Given the clamp-vs-Inf tension, when clamped / unclamped / log-space baselines are compared, then one preserves both finiteness and artifact-freedom"
verdict: PARTIAL
related: [006, 007]
tags: [bayes-factor, kde, numerics, baseline, remedy]
---

# Spike 008: Baseline Remedies — clamp vs Inf vs log-space

## What This Validates

**Given** the clamp-vs-Inf tension (Spikes 006/007), **when** three baselines are compared on the
same pairs, **then** one of them preserves both finiteness and artifact-freedom.

**The hypothesis was WRONG, and the refutation is the finding.**

## Research

No external research. The remedy is standard numerics: both failure modes share one root cause —
the tail probability is computed in **linear** space, where `1 − p` underflows to 0 once `p` is
within ~1e-16 of 1.

Remedy C computes the tail in log space. For a Gaussian KDE with bandwidth `h` on draws `xᵢ`:

```
log P    = logsumexp_i logΦ( (xᵢ − t)/h ) − log n
log(1−P) = logsumexp_i logΦ(−(xᵢ − t)/h ) − log n
logBF    = (log P − log(1−P)) − log_prior_odds
```

`logcdf(Normal(), z)` is stable far into the tail, so no clamp and no Inf. `logsumexp` is
implemented locally (3 lines) rather than adding `LogExpFunctions` to the root project, whose
co-resolution is fragile.

## How to Run

```bash
SPIKE008_N=100 julia --project=. -t auto .planning/spikes/008-bf-baseline-remedies/logspace_baseline.jl
julia --project=.planning/spikes/figenv .planning/spikes/008-bf-baseline-remedies/figures.jl
```

## Investigation Trail

1. Implemented all three baselines and ran them on the same 100 pairs (DEV seed `0x00B7A772`).
2. **Validated the new baseline against the old one where both are well-defined** — this was the
   check that made the result interpretable: median |log-space − unclamped| over the 65 finite
   pairs = **0.002**. It is the same estimator, not a different one.
3. **Expected** the log-space baseline to restore the correlation. It did the opposite:
   r̂ collapsed to **0.4448**.
4. **Chased the collapse** to its cause: on the dropped pairs the log-space baseline spans
   **[−1045, +9524]** while the amortized NRE stays bounded near ±8.
5. **Checked whether those magnitudes are meaningful.** With `L = 999` posterior draws the
   smallest empirically resolvable tail probability is 1/999, so |log-BF| beyond ~log(999) ≈ **6.9**
   is extrapolation of the Gaussian kernel's tails, not information in the draws.

## Results

**VERDICT: PARTIAL** — the mechanism question is answered decisively; the proposed remedy failed.

| baseline | finite | r̂ (all) | r̂ (survivors) | saturated |
|---|---|---|---|---|
| unclamped (status quo) | 65/100 | n/a | 0.9566 | → ±Inf |
| clamped 1e-8 | 100/100 | 0.9873 | 0.9601 | 61 at ±18.42 |
| **log-space (remedy C)** | **100/100** | **0.4448** | 0.3237 | **0** |

Agreement check: median |log-space − unclamped| = **0.002** over the 65 finite pairs (max 3831,
which is entirely the near-degenerate tail — exactly the regime under investigation).

### The finding: the clamp was hiding a real, large disagreement

Remedy C delivers what it promised — 100% finite, zero saturation, same estimator in the normal
regime — and the correlation still collapses. That is not a defect of the remedy. It reveals that
in the tail the amortized NRE (~±8) and the KDE baseline (up to ~9500) **genuinely disagree by
orders of magnitude**. Clamping both into ±18.42 made them look concordant; dropping them made the
disagreement invisible. Computing them honestly makes it unmissable.

The Phase-5 framing that `max|Δ logBF|` is "a clamped-KDE tail artifact" is therefore **incomplete**.
The clamp is not manufacturing the discrepancy; it is *masking* one.

### But the KDE baseline cannot be trusted at those magnitudes either

A log-BF of 9523 from 999 draws is not a measurement. Past |log-BF| ≈ 6.9 the value is determined
by the Gaussian kernel's assumed tail shape, not by the data. So the honest reading is:

> **In the tail regime, the KDE Bayes factor is not a valid reference at all.** Clamped it
> saturates, unclamped it diverges, log-space it extrapolates. All three failures are the same
> underlying fact: the baseline is being asked for tail probabilities its estimator cannot resolve
> from 999 draws.

The weak link in the BF gate is the **baseline**, not the amortized estimator under test.

### Implications for the gate (proposals, not decisions)

1. **Restrict the BF comparison to the regime where the baseline is valid** — e.g. pairs whose
   posterior tail probability is resolvable at the given `L` — and report the excluded fraction
   explicitly rather than silently dropping it.
2. **Raise `L`** if deeper tails must be compared: the ceiling scales as log(L), so resolving
   |log-BF| ≈ 18 needs L ≈ 6.6e7 draws — impractical. This bounds what any draw-based KDE
   reference can do.
3. **Validate agreement where it is measurable, and sign/rank agreement beyond** rather than
   demanding magnitude agreement in an unresolvable regime.
4. Reconsider whether `BF_LOGBF_TOL = 0.5` on a tail statistic is meaningful when the reference
   itself has no resolution there.

**None of these were applied.** Changing the gate is a pre-registration matter (see
`07-GATE-AMENDMENT.md` §6), not a spike decision.

### Figure

`remedies.png` — (A) survivorship bias in r̂ vs clamp ε, (B) survivors are the ambiguous middle,
(C) the three baselines against the amortized estimator (symlog), (D) |log-BF| against the
L = 999 resolution ceiling.

### What this does NOT establish

- It does not show the amortized NRE is *correct* in the tail — only that the KDE reference cannot
  adjudicate there.
- The bandwidth rule used (Silverman) approximates `KernelDensity.jl`'s default; the 0.002 median
  agreement shows this is immaterial in the normal regime but it was not tuned for the tail.
- n = 100 on one DEV seed.
