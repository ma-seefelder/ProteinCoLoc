---
spike: 011
name: sbc-power-nuisance
type: standard
validates: "Given the residual nuisance SBC failures, when the KS test's power vs marginal drift and M is computed, then it is shown whether the requirement is achievable or is test over-power"
verdict: VALIDATED
related: [009, 010]
tags: [sbc, power-analysis, over-power, nuisance, calibration]
---

# Spike 011: Is the M=2000 SBC test over-powered for non-identified nuisances?

## What This Validates

**Given** the residual SBC rejections on `autofluorescence` and `label_efficiency` after the atom
correction (Spike 010), **when** the KS test's rejection power is computed as a function of
marginal drift and M, **then** it is settled whether the requirement is achievable by a better
model or is test over-power on parameters the data cannot constrain.

No training, no model, no seed consumed — pure Monte-Carlo power analysis plus a decomposition of
the stored ranks.

## How to Run

```bash
julia --project=. .../power_curve.jl        # power table + persists power_data.jld2
julia --project=.planning/spikes/figenv .../figure.jl   # power.png
```

## Results

**VERDICT: VALIDATED — the residual is a small marginal LOCATION misplacement on non-identified
nuisances, sitting exactly at the M=2000 detection edge. It is over-power on nuisances, not a
pathology on the targets.**

### 1. What the residual failure IS (decomposition of the stored ranks)

Removing the mean shift and re-testing the shape:

| parameter | KS p (raw) | KS p (mean-centred) | reading |
|---|---|---|---|
| spillover | 2.51e-2 | 1.88e-1 | **pure location** (shape fine) |
| autofluorescence | 4.08e-4 | 1.52e-1 | **pure location** (shape fine) |
| label_efficiency | 1.00e-3 | 4.79e-2 | mostly location, a little shape |

The residual is overwhelmingly a **location misplacement of the learned marginal** — the flow
centres a non-identified parameter's marginal ~0.05–0.07 posterior-SD off the prior mean. Not
overconfidence, not image size.

### 2. What M=2000 demands (power curve)

Rejection power vs marginal shift s = mean(u)−0.5:

| s | M=250 | M=500 | M=1000 | M=2000 |
|---|---|---|---|---|
| 0.00 | 0.04 | 0.06 | 0.06 | 0.04 (≈α ✓) |
| 0.02 | 0.10 | 0.17 | 0.34 | 0.64 |
| 0.05 | 0.53 | 0.87 | 1.00 | 1.00 |
| 0.07 | 0.87 | 1.00 | 1.00 | 1.00 |

Smallest drift each M resolves at 50% power:

| M | 50%-power drift |
|---|---|
| 250 | 0.168 SD |
| 500 | 0.118 SD |
| 1000 | 0.084 SD |
| **2000** | **0.059 SD** |

### 3. The verdict, with a self-correction

The observed flow marginal error is **0.04–0.07 SD** — straddling the M=2000 50%-power line of
**0.059 SD**. At M=500 the same error would sit well below the 0.118-SD detection floor and would
NOT reject.

**Self-correction:** earlier in this session I stated the test "requires <0.02 SD marginal
accuracy". That was the 1/√M rule-of-thumb and was too pessimistic — the actual 50%-power
threshold is **0.059 SD**. The corrected reading matters: the requirement is sharp but **not
impossible**, so a capacity/training remedy is not chasing an unreachable bar.

### 4. Two legitimate remedies — both now defensible

- **Model-side (capacity):** halving the marginal error from ~0.07 to ~0.03 SD would drop
  rejection power from ~50% to ~15%; under Holm across parameters that could yield a pass. Requires
  retraining; outcome uncertain — the flow must centre a non-identified marginal to ~3% of its SD.
- **Test-side (nuisance-appropriate power):** M=2000 is correct for the TARGETS (ρ_true, Δρ) where
  high power is wanted, but demanding 0.059-SD marginal fidelity on **non-identified** parameters
  is arguably testing the wrong null. A parameter the data cannot inform has posterior = prior by
  design; whether a trained flow reproduces that prior to 3% of an SD is a property of the
  approximator, not of calibration in the sense that matters scientifically.

  **Note the pre-registration history:** `07-GATE-AMENDMENT.md` §5 *deliberately refused* a
  practical-equivalence band, to avoid a p-hacking appearance. This spike supplies what was missing
  then — an **outcome-independent** justification (the power analysis does not depend on the
  results) for treating non-identified nuisances differently. Any such change is still a
  pre-registration matter under §6.

### Figure

`power.png` — (A) rejection power vs drift for M ∈ {250, 500, 1000, 2000} with the observed
nuisance-drift band; (B) the smallest drift each M resolves, against the observed flow error band.

### What this does NOT establish

- It does not prove capacity WILL close the gap — only that the gap is not unreachable in
  principle.
- The shift model is leading-order; `label_efficiency` carries a small extra shape component the
  model omits.
- It does not decide which remedy to take — that is a project/publication call.
