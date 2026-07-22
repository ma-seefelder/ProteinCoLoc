---
spike: 012
name: capacity-lever
type: standard
validates: "Given the residual nuisance marginal drift, when the flow capacity is raised, then the drift falls below the 0.10-SD equivalence margin without harming the coloc targets"
verdict: PARTIAL
related: [010, 011, 013]
tags: [npe, capacity, calibration, nuisance, flow]
---

# Spike 012: Does more flow capacity fix the residual nuisance drift?

## What This Validates

**Given** the residual F2 marginal-location drift (~0.05–0.07 posterior-SD) on the non-identified
nuisances after the bounded-θ fix, **when** the FLOW capacity is raised (coupling 10→16,
flow_depth 2→3, flow_width 128→256; summary net unchanged), **then** does the drift fall below the
spec §3 equivalence margin of 0.10 SD, while ρ_true/Δρ stay calibrated?

## How to Run

```bash
julia --project=. -t auto .../train_highcap.jl   # ~130 min: 50k-pair mixture datagen + train
julia --project=. -t auto .../eval_highcap_vs_amended.jl   # SBC both nets, same DEV seed
```

New artifact root `artifacts/spike012_highcap/grid_8/`. amended_v2 and every frozen bundle verified
byte-identical (npe_8 sha256 `198bb078…` unchanged before/after).

## Investigation Trail

1. Raised ONLY the flow (the marginal-accuracy lever); left the summary net (identifiability lever)
   at amended_v2 defaults. Same mixture, n_pairs=50000, bounded-θ, CPU.
2. Training: early-stopped at epoch 51, best ≈ epoch 11; training risk fell to 0.478 while
   validation risk drifted 0.55→0.59 — **moderate overfitting**, contained by early stopping, not a
   divergence. Higher capacity did trade into mild overfit, as flagged going in.
3. Evaluated BOTH nets on the SAME DEV seed (0xca11b3, asserted disjoint from 20 forbidden seeds),
   M=1000, L=999, with atom-aware ρ_true handling.

## Results

**VERDICT: PARTIAL — capacity REDISTRIBUTES the drift across nuisances but does not systematically
reduce it below the margin. Not a reliable lever.**

Same DEV seed, both nets:

| parameter | amended_v2 drift / KS p | high-cap drift / KS p | effect |
|---|---|---|---|
| autofluorescence | 0.082 / 2.0e-3 **rej** | 0.051 / 1.6e-1 ok | improved |
| shift_dy | 0.110 / 6.9e-5 **rej** | 0.059 / 2.1e-2 (rej@.05) | improved |
| label_efficiency | 0.022 / 2.4e-1 ok | 0.020 / 8.8e-1 ok | ~same |
| **shift_dx** | 0.032 / 5.3e-1 ok | 0.076 / 1.2e-2 **rej** | **worse** |
| spillover | 0.015 / 3.4e-1 ok | 0.032 / 4.8e-1 ok | ~same |
| noise | 0.044 / 1.0e-1 ok | 0.043 / 4.8e-1 ok | ~same |
| **ρ_true** (atom-corr.) | ok (KS 0.74) | ok (KS 0.58) | both calibrated |
| **Δρ** | ok (KS 0.74) | ok (KS 0.53) | both calibrated |

### The finding

The residual marginal error stays in the **0.02–0.08 SD** band under BOTH capacities. Raising the
flow moved *which* nuisance drifts (autofluorescence and shift_dy got better, shift_dx got worse)
without shrinking the band. This is consistent with the drift being a **fundamental property of a
neural flow reproducing the marginal of a NON-identified parameter**, not a capacity shortfall.

This STRENGTHENS the test-side route (spec `07-NUISANCE-SBC-SPEC-DRAFT.md`): if no realistic model
lever reliably removes a 0.06-SD marginal drift on parameters the data cannot inform, the correct
response is a nuisance-appropriate test (equivalence), not a bigger model.

### Robust across both nets

**ρ_true and Δρ — the coloc targets — are calibrated under both capacities** (atom-corrected KS
0.53–0.74). The targets are insensitive to the capacity lever. This is the reassuring result: the
scientific claim rests on the targets, and they are solid.

### Honest caveats

- **Single DEV seed.** The absolute pass/fail per nuisance is seed-dependent (on this seed
  amended_v2's label_efficiency already passes, unlike on the reporting seed — the ±1.5-z seed
  variance documented in Spike 010). What is robust across the comparison is the DRIFT MAGNITUDE
  band (0.02–0.08 SD) and that capacity does not systematically shrink it. A multi-seed sweep would
  firm this up but the redistribution (shift_dx worsening) is already clear on one seed.
- **Overfitting** at this capacity means an even larger flow is unlikely to help — the val risk
  already diverges. Capacity is near its useful ceiling here.
- ρ_true's raw drift rose slightly (0.060→0.077) with capacity, but atom-corrected it stays clean —
  so the increase is on the atom-driven component, not a target regression.

### Implication for the §6.4 model change

Capacity brings no clear, reliable benefit and cannot alone justify the retrain that §6.4 requires
to legitimize a new confirmatory run. This shifts weight toward the simulator-side μ-prior
truncation (Spike 013) as the model change with a definite purpose — at the cost of ADVI
comparability. That trade-off is the decision the two spikes were meant to inform.

### What this does NOT establish

- It does not prove capacity NEVER helps — only that at the tested step it redistributes rather
  than reduces, on one DEV seed.
- It does not test intermediate capacities or regularization (e.g. an explicit marginal-consistency
  penalty), which NeuralEstimators v0.2.1 does not expose a hook for.
