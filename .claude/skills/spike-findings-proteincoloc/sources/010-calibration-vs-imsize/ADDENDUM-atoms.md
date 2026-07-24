# Spike 010 ADDENDUM — the ρ_true "overconfidence" was an ATOM ARTIFACT (retraction)

**This addendum RETRACTS the overconfidence conclusion in the parent README.** The U-shaped
ρ_true rank histogram is not caused by a too-narrow posterior. It is caused by atoms in the prior.

## The mechanism, confirmed to 3 decimal places

`sample_prior` draws `μ* ~ Truncated(Cauchy(0, 0.3), −1, 1)` (the Turing μ-prior, a hard project
constraint) and sets `ρ_true = ghat(μ*)`. `ghat` (`src/amortized/simulator.jl`, frozen SIM-02) is a
**clamped** piecewise-linear inverse over the physically realized sweep range
`μ ∈ [−0.67976, 0.847149]`. Cauchy tails extend well past that, so all such μ map to exactly
±0.99.

Predicted from the prior alone vs observed at M = 2000:

| | predicted | observed |
|---|---|---|
| P(ρ_true = −0.99) | 0.0485 | 0.0450 (90/2000) |
| P(ρ_true = +0.99) | 0.0191 | 0.0195 (39/2000) |
| total atom mass | **6.76%** | **6.45%** |

## Why this invalidates the SBC test for ρ_true

SBC assumes a continuous prior. When θ* sits exactly on an atom, essentially every posterior draw
lies on one side of it, so the rank is forced to an extreme. Measured:

| subset | n | tails | middle | ratio |
|---|---|---|---|---|
| at the lower atom | 93 | 0.968 | 0.000 | degenerate |
| at the upper atom | 40 | 0.975 | 0.000 | degenerate |
| **interior** | **1867** | **0.196** | **0.187** | **1.05 (flat)** |

KS p for ρ_true (DEV seed, M = 2000):

| treatment | n | KS p | verdict |
|---|---|---|---|
| all draws (what the gate does) | 2000 | 5.88e-04 | FAIL |
| **atoms excluded** | 1871 | **0.667** | **PASS** |
| **randomized ranks** (standard correction for mixed distributions) | 2000 | **0.460** | **PASS** |

Counter-check against overconfidence: if the posterior were too narrow, tail cases would have
visibly narrower posteriors. Median post_sd is 0.0654 on tail cases vs 0.0678 on middle cases —
no difference. **ρ_true is calibrated.**

ρ_true is the ONLY parameter carrying atoms.

## The deeper issue: the prior places mass outside the forward model's achievable support

ρ is a correlation, bounded by ±1; the ρ-sweep to ±0.99 realized only μ ∈ [−0.68, 0.847]. So 6.8%
of the μ-prior asks for a mean patch correlation **no ρ_true can induce**. Clamping is the correct
local response to an impossible request, but the prior itself is inconsistent with the simulator.

Consequences beyond the test:
1. ~6.5% of every training set has ρ_true pinned at exactly ±0.99.
2. The trained net has **learned** this: 7.6% of posterior means sit at |ρ̂| > 0.95 (112 of them on
   atom draws), and on atom draws the median ρ̂ is −0.9701 against a true −0.99 — reproduced almost
   exactly. The net learned a degenerate prior faithfully.
3. Real data will essentially never have ρ_true exactly ±0.99.
4. The asymmetry (4.85% low vs 1.91% high) puts the heavier clamp on the negative μ tail, which
   STATE.md records for SIM-02 as **"prior-only"** — i.e. the range never anchored by real data.

## Remedies

- **Test-side (no retraining):** randomized ranks, or exclude atoms and report the excluded
  fraction. Mathematically required whenever atoms exist. Does NOT remove the degenerate prior.
- **Simulator-side (requires retraining):** truncate the μ-prior to the achievable support
  `[GHAT_MU_MIN, GHAT_MU_MAX]`. ρ_true then becomes continuous, SBC is valid by construction, and
  no boundary spike enters training. **Note:** this deviates from the stated project constraint
  that the simulator prior match the Turing `@model` μ-prior — a deviation that needs an explicit
  decision, though the current setup arguably already violates that constraint's intent by
  implying μ values the forward model cannot produce.

Neither was applied. Both are pre-registration matters under `07-GATE-AMENDMENT.md` §6.
