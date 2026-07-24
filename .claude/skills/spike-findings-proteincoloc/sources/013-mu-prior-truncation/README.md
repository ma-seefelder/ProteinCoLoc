---
spike: 013
name: mu-prior-truncation
type: standard
validates: "Given the ρ_true atoms from the clamped ghat, when the μ-prior is truncated to the achievable support, then the atoms vanish and ρ_true SBC is valid without randomized ranks — at what cost?"
verdict: VALIDATED
related: [010, 012]
tags: [simulator, prior, atoms, sbc, calibration, tradeoff]
---

# Spike 013: μ-prior truncation — does it fix the atoms at the root, and what does it cost?

## What This Validates

**Given** the ρ_true atoms at exactly ±0.99 (~6.5% of prior draws, from the clamped `ghat` acting on
Truncated-Cauchy μ tails, Spike 010 addendum), **when** the μ-prior is truncated to the achievable
support `[GHAT_MU_MIN, GHAT_MU_MAX] = [-0.67976, 0.847149]`, **then** ρ_true becomes continuous over
[-0.99, 0.99] and SBC is valid without the randomized-rank correction — **at what cost?**

This CHANGES the prior, so the truncated net is a different statistical model; the amended_v2
comparison on ρ_true is not apples-to-apples (different priors). It also breaks the CLAUDE.md
ADVI-comparability constraint (a known, separately-decided cost). The spike measures the technical
outcome and the high-correlation inference cost.

## How to Run

```bash
julia --project=. -t auto .../train_trunc.jl   # ~80 min: truncated-prior mixture datagen + train
julia --project=. -t auto .../eval_trunc_vs_amended.jl
```

Custom `sample_prior_trunc` + `trunc_datagen` injected via `_train_grid_pipeline`'s `datagen` seam;
src/ NOT edited. New root `artifacts/spike013_trunc/grid_8/`. Architecture = amended_v2 defaults
(10/2/128) so ONLY the μ-prior differs — isolating the truncation effect from Spike 012's capacity
effect. Baseline artifacts verified byte-identical (amended_v2 npe `198bb078`).

## Results

**VERDICT: VALIDATED — truncation removes the atoms and makes ρ_true SBC natively clean, but
degrades inference near |ρ|≈0.95–0.99 and does not touch the nuisance drift. Randomized ranks
achieve the same SBC validity at none of these costs.**

### §A — Atoms: eliminated

| net | atom fraction (ρ_true = ±0.99 exactly) | range |
|---|---|---|
| truncated | **0 / 1000 = 0.00%** | [-0.9876, 0.9809] |
| amended_v2 | 64 / 1000 = 6.40% | [-0.99, 0.99] |

### §B — ρ_true SBC WITHOUT any test correction (the headline)

| net | KS raw (all draws) | atoms-excluded | randomized-rank |
|---|---|---|---|
| **truncated** | **0.342 (clean)** | 0.342 | 0.342 |
| amended_v2 | 0.0096 (**reject**) | 0.086 | 0.428 (clean) |

The truncated net passes ρ_true SBC **natively** — all three columns identical because there are no
atoms to correct. amended_v2 needs randomized ranks to pass (raw 0.0096 → randomized 0.428). Both Δρ
calibrated (trunc 0.117, amended 0.304). This is the truncation's real win: no test-side correction
needed for ρ_true.

### §C — THE COST: high-|ρ| inference (paired common test set, 60 reps/level)

Δ(mean|bias|) = truncated − amended_v2 (positive ⇒ truncation WORSE):

| ρ* | −0.99 | −0.95 | −0.90 | −0.85 | −0.80 | +0.80 | +0.85 | +0.90 | +0.95 | +0.99 |
|---|---|---|---|---|---|---|---|---|---|---|
| Δ\|bias\| | **+0.039** | +0.013 | −0.016 | −0.014 | −0.008 | −0.003 | −0.015 | −0.003 | +0.012 | **+0.041** |

**The trade-off is clear and localized:** truncation is WORSE at the extremes (|ρ|≥0.95 — bias
roughly DOUBLES at ±0.99: 0.020→0.059 and 0.017→0.058), because it removed the training mass that
taught the net perfect correlation. At |ρ|≤0.90 truncation is actually equal or slightly BETTER —
it no longer spends capacity on the ±0.99 spike and concentrates on the achievable range.
(Coverage=0.000 at exactly ±0.99 for BOTH nets is a boundary artifact: θ* sits on the support edge,
so no symmetric interval covers it — not a truncation effect.)

### §D — Nuisance drift: unchanged (as predicted)

Truncation changes the μ-prior, not the nuisance priors, so the F2 nuisance drift persists at the
same 0.04–0.09 SD magnitude, merely reshuffled (trunc rejects label_efficiency 0.085 / shift_dy
0.088; amended rejects spillover 0.076 / shift_dy 0.062). **Truncation does NOT fix the residual
nuisance failures** — those remain a test-side (equivalence) matter regardless.

### Training

80 min, early-stopped epoch 45, val risk drifted 0.55→0.61 (mild overfit, contained) — same profile
as amended_v2. Architecture identical to amended_v2 defaults.

## Decision input: truncation vs randomized ranks

Both make ρ_true SBC valid. The comparison:

| | randomized ranks (test-side) | μ-prior truncation (model-side) |
|---|---|---|
| ρ_true SBC | clean (0.43) | clean (0.34) |
| cost — high-|ρ| inference | none (net unchanged) | **worse at \|ρ\|≥0.95** |
| cost — ADVI comparability | none (prior unchanged) | **broken** (Phase-4 benchmark) |
| nuisance drift | unchanged | unchanged |
| statistical standing | textbook correction for a mixed (discrete+continuous) prior | changes the model to simplify the test |

**Randomized ranks solve the same problem at neither cost.** The atoms are a legitimate feature of
the Turing prior (μ beyond the realizable range ⇒ maximal correlation); the standard treatment of a
mixed distribution is the randomized rank, not a prior change. Truncation's only compensating merit
is a manuscript that needs no "we randomized atom ranks" sentence — a cosmetic gain against a real
inference cost at the extremes plus the ADVI break.

## Implication for the §6.4 model change

Spike 012 (capacity) brought no reliable benefit; Spike 013 (truncation) fixes only the atoms —
already solved test-side — at a real cost. **Neither lever earns the model change that §6.4 requires
to legitimize a new confirmatory run.** That sharpens the strategic fork: either take truncation as
the §6.4 change anyway (accepting its costs for a correction-free ρ_true), or accept that no model
change is warranted and route to Go/No-Go with the existing FAIL, publishing with named limits.

## What this does NOT establish

- Single DEV seed for the SBC comparison (the ±1.5-z seed variance applies).
- The high-|ρ| cost is measured at 60 reps/level; the ±0.99 doubling is clear but the mid-range
  "slightly better" differences are within noise.
- It does not weigh the ADVI-comparability break — that is a project decision, not a spike result.
