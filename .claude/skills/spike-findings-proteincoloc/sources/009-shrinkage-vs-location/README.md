---
spike: 009
name: shrinkage-vs-location
type: standard
validates: "Given spillover rejects SBC (KS p=6.8e-6) at shrinkage 1.01, when the rank distribution is decomposed into location vs shape, then the rejection is explained without invoking a bug in the rank machinery"
verdict: VALIDATED
related: [010, 011]
tags: [sbc, calibration, diagnostics, vacuity]
---

# Spike 009: Shrinkage measures width, not location

## What This Validates

**Given** `spillover` rejects SBC uniformity (KS p = 6.8e-6) while its shrinkage is 1.010 — i.e.
the posterior is as wide as the prior and the parameter is flagged *vacuous* — **when** the rank
distribution is decomposed into a location component and a shape component, **then** the rejection
is explained without any defect in the rank machinery.

This began as a suspected contradiction: "if the posterior *is* the prior, ranks must be uniform by
construction, so a rejection implies a bug." That framing was wrong, and finding out why took one
command rather than a spike — recorded here because the conclusion changes how every other SBC
table in this project must be read.

## How to Run

```bash
julia --project=. .planning/spikes/009-shrinkage-vs-location/location_vs_shape.jl
```

Reads the committed confirmatory report `artifacts/amended_v2/grid_8/gate_report_8.jld2`. No
inference, no training, no seed consumed.

## Investigation Trail

1. Suspected a bug in rank computation or a marginal-vs-prior inconsistency.
2. Computed the normalized-rank mean `u = (rank+0.5)/(L+1)` per parameter and its z-score against
   the uniform null (`SE = sqrt((1/12)/M)`, M = 2000 → SE = 0.00645).
3. The suspicion dissolved immediately: `spillover` has z = +4.3. It is a **location shift**.
4. Realised the premise was the error: **shrinkage = post_sd / prior_sd is a ratio of WIDTHS.**
   A posterior can match the prior's width exactly and still be displaced. "Vacuous" (as the gate
   defines it) means *uninformative about spread*, not *equal to the prior*.

## Results

**VERDICT: VALIDATED — no bug. The rejection is a location shift.**

| parameter | mean(u) | z (location) | shrinkage | KS p | location explains it? |
|---|---|---|---|---|---|
| **ρ_true** | 0.4666 | **−5.2** | 0.130 | 1.12e-07 | **yes** |
| **spillover** | 0.5278 | **+4.3** | 1.010 | 6.82e-06 | **yes** |
| autofluorescence | 0.5165 | +2.6 | 0.997 | 7.61e-02 | (not rejected) |
| **label_efficiency** | 0.5119 | +1.8 | 0.923 | 1.92e-04 | **no → shape** |
| shift_dx | 0.4946 | −0.8 | 0.976 | 8.83e-01 | (not rejected) |
| **shift_dy** | 0.4857 | −2.2 | 0.987 | 1.67e-03 | **no → shape** |
| **noise** | 0.4880 | −1.9 | 1.025 | 1.92e-02 | **no → shape** |
| Δρ | 0.5081 | +1.3 | 0.136 | 1.96e-01 | (not rejected) |

### Two distinct failure modes coexist

- **Location shifts** — ρ_true (−5.2) and spillover (+4.3). This is the finding-F2 marginal-drift
  mechanism (NPE's forward-KL constrains nothing about ∫q(θ|Z)p(Z)dZ = π(θ)). The bounded-θ fix
  reduced it but did **not** eliminate it. Crucially it acts on *location*, so it reaches even
  parameters whose *width* the summary never constrains — which is exactly `spillover`.
- **Shape failures** — label_efficiency, shift_dy, noise reject with |z| < 3. Their rank
  histograms are misshapen in a way a mean cannot capture. Mechanism unknown; **not** explained by
  F2.

### The finding that matters most

**ρ_true has acquired a location shift of z = −5.2 under the realistic image-size mixture.** The
same net's ρ_true was z = −1.25 at 256²-only in the bounded-θ DEV diagnostic. This is a direct lead
into Spike 010 (why calibration degrades with image size) and reframes it: the question is not only
"why does it get worse" but "why does ρ_true drift *low*".

### Consequence for reading every SBC table in this project

`vacuous = shrinkage ≥ 0.95` (introduced in commit `73f32b2`) flags parameters that are
uninformative **about spread**. It does **not** mean "posterior = prior", and a vacuous parameter
can still legitimately reject on location. Conversely a clean KS p on a vacuous column remains
weak evidence, as originally intended. Both readings must be stated whenever the flag is used.

### What this does NOT establish

- It does not explain the shape failures (label_efficiency, shift_dy, noise).
- It does not explain *why* ρ_true drifts low under the mixture — that is Spike 010.
- z uses the uniformity-null SE, which is anti-conservative when ranks are non-uniform; read the
  z values as effect sizes, not as tests.
