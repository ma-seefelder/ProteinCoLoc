---
spike: 010
name: calibration-vs-imsize
type: standard
validates: "Given rho_true's SBC collapse under the image-size mixture, when image size is swept with everything else held fixed, then the mechanism is identified"
verdict: VALIDATED
related: [009, 011]
tags: [sbc, calibration, overconfidence, imsize, covariate-shift]
---

# Spike 010: What actually breaks ρ_true's calibration

## What This Validates

**Given** the confirmatory amended gate's ρ_true failure (KS p = 1.12e-7, z = −5.2) under the
realistic image-size mixture, **when** image size is swept with everything else held fixed,
**then** the mechanism is identified.

**Answer: the ρ_true posterior is OVERCONFIDENT (too narrow). Image size is not the driver.**

Three hypotheses were entered and two were refuted. The eventual answer came from a statistic none
of the first three probes measured.

## How to Run

```bash
SPIKE010_M=200   julia --project=. -t auto .../calib_vs_imsize.jl      # per-size ranks + H1 probe
SPIKE010B_M=200  julia --project=. -t auto .../location_error.jl        # low-variance location bias
SPIKE010C_M=2000 julia --project=. -t auto .../mixture_replicate.jl     # gate replication, DEV seed
                 julia --project=. .../shape_analysis.jl                # rank-histogram shape
```

All runs use DEV seeds asserted disjoint from `PROD_SEED_V2`, every v1 `PROD_SEED[G]`,
`VAL_MASTER_SEED` and `NPE_MASTER_SEED`. Nothing retrained, no gate report written.

## Investigation Trail

1. **H1 — summary not size-invariant.** At fixed θ the mean patch correlation *does* move with
   image size. But over four θ the drift **changes sign** (θ1 positive, θ2–θ4 negative; mean
   −0.010 to −0.017). It is a θ×size interaction, not an additive offset.
   *Self-correction:* an earlier single-θ probe gave −0.032 and was reported as "comparable to the
   whole label_efficiency prior range". That was an overstatement from n=1.
2. **Per-size SBC ranks (M=200).** Underpowered by construction: the SE of z is 1.0 at M=200 and
   the effect sought is ≈1.6 in those units. No conclusion possible. *Design error on my part* —
   and the script stored only ranks and posterior SDs, not θ and posterior means, violating the
   "store the intermediates" convention this repo learned in Spike 006.
3. **Location error (M=200), the low-variance route.** SE ≈ 0.005 on the raw bias — an order of
   magnitude sharper at the same cost. Result: ρ_true's per-size bias is small and **non-monotone**
   (+0.021 / +0.001 / +0.006 / −0.006 / +0.011 across 256²→2048²). Mixture-weighted standardized
   bias ≈ **+0.014**, whereas the gate implies ≈ **−0.085**. Wrong sign, wrong magnitude.
   **H1 refuted as the explanation.** H2 (sharpening) refuted too: ρ_true's width falls only ~20%
   and label_efficiency's *increases*.
4. **Replication of the gate's SBC arm on a DEV seed (M=2000, same mixture/weights/L).**
   ρ_true: mean(u) 0.4666 → **0.4891**, z −5.2 → **−1.69**, standardized location bias
   **+0.0013 ± 0.0226 — null**. Yet KS still rejects (5.88e-4).
   *A distribution cannot be both centred and rejected for uniformity unless its SHAPE is wrong.*
5. **Shape analysis** (free — the ranks were stored this time).

## Results

**VERDICT: VALIDATED — ρ_true's posterior is overconfident.**

Rank-histogram shape on the DEV replication (M = 2000; uniform expects tails = middle = 0.20;
binomial SE on each fraction = 0.009):

| parameter | tails | middle | ratio | reading |
|---|---|---|---|---|
| **ρ_true** | **0.247** | **0.174** | **1.42** | **U-shaped → posterior TOO NARROW (overconfident)** |
| spillover | 0.183 | 0.208 | 0.88 | flat |
| autofluorescence | 0.213 | 0.210 | 1.01 | flat |
| label_efficiency | 0.218 | 0.196 | 1.12 | flat |
| shift_dx | 0.207 | 0.194 | 1.07 | flat |
| shift_dy | 0.204 | 0.209 | 0.98 | flat |
| noise | 0.186 | 0.202 | 0.92 | flat |

ρ_true's tail mass is **5.2σ** above expectation. This is a **dispersion** failure, not location
and not image size. Every earlier probe in this spike measured the *first* moment while the defect
lives in the *second*.

Consistently, ρ_true has by far the narrowest posterior (shrinkage 0.130 — 7.7× narrower than its
prior): the net compresses hardest exactly the parameter the summary constrains best, and
overshoots.

### This is a recurrence, not a new defect

`.planning/STATE.md` records for Phase 5-04: *"SBC every-parameter KS/chi2 p~0
(**overconfident**/miscalibrated NPE)"*. The iter1 capacity increase reduced it and the bounded-θ
fix (F2) removed the *location* component from label_efficiency — but ρ_true's overconfidence has
persisted since Phase 5 and was masked by the location effects that dominated the diagnostics.

### Secondary result: image size IS a real covariate, via generalization not drift

`label_efficiency` at 256² carries a standardized bias of **+1.100 ± 0.12** (raw +0.083),
an order of magnitude larger than at any other size — and 256² is **outside** the amended net's
training mixture (512²–2048²). Large bias off the training support is expected behaviour, but it
confirms the net does **not** generalize across image size. Finding F5's concern stands; its
mechanism is out-of-support generalization failure, not within-range drift.

### Methodological result: the one-shot gate has real seed variance

Gate vs DEV replication, per-parameter z: ρ_true −5.2 → −1.69, spillover 4.3 → 2.13,
autofluorescence 2.6 → 3.77, label_efficiency 1.8 → 3.25, shift_dx −0.8 → −1.85,
shift_dy −2.2 → −0.80, noise −1.9 → −1.61. Mean |difference| ≈ 1.6 z-units, broadly **consistent
with sampling noise** (two independent M=2000 runs are expected to differ by ~1.1–1.4).

So the gate's numbers are **not** an artifact — but z = −5.2 was a high draw, and any single
pre-registered run reports per-parameter statistics with roughly ±1.5 z-units of seed noise.
For a one-shot protocol whose verdict is binding, that is decision-relevant and should be stated
wherever those per-parameter numbers are quoted.

### What this does NOT establish

- It does not explain the *other* parameters' KS rejections (all flat in shape and near-null in
  location, yet autofluorescence rejects at 4.1e-4). Unresolved.
- It does not show that fixing the overconfidence would make the gate pass.
- Shape was assessed with a tail/middle ratio, not a formal dispersion test.
- Single DEV seed for the replication.

### Follow-ups (not actioned)

Overconfidence is a training-side property (flow capacity, epochs, the KL objective's tendency to
under-disperse). Candidate remedies: posterior-width regularization, an ensemble, or explicit
calibration of the ρ_true marginal. All require retraining and would need a fresh pre-registration
per amendment §6.
