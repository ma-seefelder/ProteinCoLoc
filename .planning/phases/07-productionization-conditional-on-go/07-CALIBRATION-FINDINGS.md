# 07 — Calibration Findings: what the SBC failures actually are

**Date:** 2026-07-21
**Scope:** read-only forensic investigation of the SBC arm of the frozen grid-4 / grid-8 / grid-16
ship-gate bundles. No retraining, no re-run of any reported gate arm, no pre-registered constant
touched.
**Evidence base:** `artifacts/grid_{4,8,16}/gate_report_{4,8,16}.jld2`, field `report.sbc.ranks`
(M = 2000 SBC draws, L = 999 posterior draws per draw), plus targeted read-only re-inference
against the same frozen nets for the shrinkage and covariate-shift probes.

**Epistemic labels used throughout — read them literally:**

| Label | Meaning |
|---|---|
| **MEASURED** | A number computed directly from the frozen artifacts or from a read-only re-run against the frozen nets. Reproducible. |
| **INFERRED** | A mechanism proposed to explain measured numbers. Consistent with the evidence and, where stated, supported by sign/monotonicity checks — but not itself directly measured. |
| **NOT CONCLUSIVE** | Evidence was gathered and it does **not** settle the question. Explicitly not a finding. |

---

## Severity ordering

| # | Finding | Severity |
|---|---|---|
| **F5** | **Covariate-shift risk: the calibration proof holds at 256² and is not demonstrated to transfer to the 1376×1028 paper regime** | **HIGHEST — blocks the Phase-8 anchor / Phase-16 blind-evaluation reading of the SBC proof** |
| F3 | grid 4's `label_efficiency` SBC pass is vacuous (posterior = prior) | High — a pass that is not calibration evidence |
| F1 | `label_efficiency` SBC failure is a location bias, not overconfidence | Medium — reframes the failure, changes the remedy |
| F2 | Mechanism: learned-conditional location displacement | Medium — explains F1/F3 |
| F6 | Provenance gap: `imsize_set` was never persisted | Medium — blocks verification of F5 for grid 16 |
| F4 | Ruled-out hypotheses, with evidence | Informational |

---

## F5 — COVARIATE-SHIFT RISK (highest severity)

**The claim:** the SBC calibration proof is established at a training/evaluation image size of
256×256 (512×512 for grid 16). It is **not demonstrated** to transfer to the real paper regime of
**1376×1028**, which is precisely the regime Phase-2 decision **D-08** anchors the design on and
which the Phase-8 anchors and the Phase-16 blind evaluation depend upon.

**MEASURED:**

- Training runs 07-05 / 07-06 constrained `imsize_set` to `((256,256),)` as a *documented
  compute-budget deviation*.
- The gates evaluate at `SBC_IMSIZE` = (256,256) for grids 4 and 8, and (512,512) for grid 16.
  Train and gate are therefore **self-consistent** — the recorded verdicts are internally valid.
- Re-evaluating the **grid-8** net on the mixed `imsize` set instead of 256²-only moves
  `label_efficiency` mean(u) from **0.546 → 0.807** (z ≈ 16.8 at M = 250), bias **−0.130**.
- At **fixed θ**, changing image size shifts the mean per-patch correlation by roughly the
  **entire `label_efficiency` prior range**. At G = 8: mean per-patch correlation is **0.513 at
  2048²** vs **0.470 at 256²** at `label_efficiency` = 0.95 (Δ ≈ 0.043), while sweeping
  `label_efficiency` across 0.65 → 0.95 moves it only **0.045**.

**INFERRED:** image size acts on the summary statistic with a magnitude comparable to the physical
parameter the summary is supposed to identify. A net trained at one image size therefore has no
demonstrated calibration at another, and the confound is strong enough to swamp the signal.

**Why this is the highest-severity item:** D-08 anchors the design on the real paper regime
(`test/test_images/{positive,negative}/*_c{1,2}.tif`, 16-bit, 1376×1028) *precisely so the summary
is comparable to real inputs*. The measured shift means the chain "SBC passes ⇒ posteriors on real
images are calibrated" **does not currently hold**. Every downstream use of the calibration proof
on real-image dimensions — the Phase-8 physical anchors, the Phase-16 blind evaluation — inherits
this gap.

**What would close it (not done here, requires retraining — out of scope):** train and gate at (or
across) the real image dimensions, or demonstrate summary invariance to image size and gate that
invariance explicitly.

---

## F3 — grid 4's `label_efficiency` SBC pass is VACUOUS

**MEASURED (grid 4, `label_efficiency`):**

| Quantity | Value | Reference |
|---|---|---|
| post_sd / prior_sd | **0.994** | 1.0 = posterior is the prior |
| shrinkage centre | **0.8034** | prior mean = 0.80 |
| tercile mean(u) | 0.18 / 0.49 / 0.82 | the "posterior = prior" signature |
| KS p | 0.567 | "passes" |

**INFERRED:** grid 4 passes SBC on `label_efficiency` **because it learns nothing about it**. A
posterior that reproduces the prior yields uniform ranks *by construction* — rank uniformity is a
**necessary but not sufficient** condition for calibration.

> **This must NOT be read as calibration evidence anywhere** — not in the Go/No-Go memo, not in the
> manuscript, not in any per-grid comparison table. A vacuous pass and an informative pass are not
> the same result and must never be tabulated as if they were.

**Mitigation shipped with this document (reporting only):** `sbc_gate` now reports a per-parameter
`shrinkage = post_sd / prior_sd`, its `post_sd` / `prior_sd` components, and a boolean `vacuous`
flag (cutoff `SBC_VACUOUS_SHRINKAGE = 0.95`, defined in `test/gate/sbc.jl`, **not** in any frozen
`gate_consts_<G>.jl`). No pass/fail threshold and no recorded verdict changed.

---

## F1 — the `label_efficiency` SBC failure is a LOCATION BIAS, not overconfidence

Normalized rank `u = (rank + 0.5) / (L + 1)`; under uniformity E[u] = 0.5 and
SE = sqrt((1/12)/M) = sqrt((1/12)/2000) = **0.00645**.

**MEASURED — `label_efficiency`:**

| Grid | mean(u) | z | KS p |
|---|---|---|---|
| 4 | 0.496 | −0.6 | 0.567 |
| 8 | **0.548** | **+7.4** | **1.84e-9** |
| 16 | **0.556** | **+8.6** | **7.15e-13** |

**MEASURED — dispersion (the overconfidence discriminator):**

Tail mass = 0.199 / 0.231 / 0.215 (grids 4 / 8 / 16) against an expected **0.20**. The rank
histograms are **flat, not U-shaped**.

**INFERRED:** overconfidence produces a **U-shaped** rank histogram (excess mass in both tails).
That is not what is present. The dispersion is normal; the *location* is displaced. A high mean
rank means θ* falls above most posterior draws, i.e. the NPE **systematically underestimates**
`label_efficiency`.

**Secondary cases with the same signature (MEASURED):** `noise` at grid 8 (z = −3.5), `ρ_true` at
grid 16 (z ≈ +2.3).

**Consequence for the remedy:** this is *not* fixed by widening posteriors, tempering, or more
capacity — the standard overconfidence levers. It is a bias/location problem.

---

## F2 — Mechanism: learned-conditional location displacement

**INFERRED mechanism.** The `label_efficiency` prior is `Uniform(0.6, 1.0)`. The
`NormalisingFlow` used by the NPE is **unbounded** and cannot represent the truncated box. The
NPE's forward-KL training objective constrains nothing about the implied marginal
∫ q(θ|Z) p(Z) dZ = π(θ) — nothing in the loss forces the average posterior to reproduce the prior.
So the learned conditional is free to sit off-centre, and does.

**MEASURED — shrinkage and centre for `label_efficiency` (prior mean = 0.80):**

| Grid | post_sd / prior_sd | shrinkage centre | E[posterior] − θ* |
|---|---|---|---|
| 4 | 0.99 (= prior) | 0.8034 | — |
| 8 | 0.83 | 0.7877 | **−0.018** |
| 16 | 0.74 | 0.7651 | **−0.023** |

**MEASURED — out-of-support leakage (sign-matched to each parameter's bias):**

| Grid | mass < 0.6 | mass > 1.0 |
|---|---|---|
| 4 | symmetric | symmetric |
| 8 | 0.0143 | 0.0059 |
| 16 | 0.0217 | 0.0061 |

**The key observation:** information gained, centre displacement, and rank bias all move together
**monotonically 4 → 8 → 16**. The bias **grows with learning, not with ignorance** — the better the
net identifies the parameter, the further off-centre the learned conditional sits. This is the
opposite of an under-training signature and is the single strongest argument that the mechanism is
the unbounded-flow / unconstrained-marginal one above.

---

## F4 — Ruled out, with evidence

Each of these was a live hypothesis and each is **MEASURED-ruled-out**, not assumed away.

1. **Empty / degenerate patch masking — RULED OUT.**
   Measured missing-patch fraction = **0.0000** and both-channel survivor fraction ≈ **1.0000**
   for *every* (G, imsize) combination. The smallest patch is 32×32 = **1024 px**, never
   approaching the ≥15-survivor floor. **Consequence:** the mask rows `G²+1 : 2G²` are constant
   1.0 and therefore carry **no information** at all.

2. **Mask-row / standardization mishandling at G ≠ 8 — RULED OUT.**
   `_summary_row_partition(:min, nrows)` derives the split as `nrows ÷ 2` — there is no
   8-specific literal anywhere in the partition logic. And per (1) the mask half is inert
   regardless, so it cannot be the source of a grid-dependent bias.

3. **Identifiability collapse as the CAUSE — RULED OUT as a cause, RETAINED as the F3
   explanation.** A prior-equal posterior yields uniform ranks by construction, so collapse cannot
   *cause* a rank failure — but it *is* exactly the explanation for grid 4's vacuous pass (F3).

4. **Prior-boundary leakage as a sufficient cause — RULED OUT.** The leakage is real and
   sign-consistent (F2 table) but amounts to only ~2% of posterior mass — far too small to produce
   the measured −0.023 shift on its own.

5. **Under-training at fine grids — RULED OUT as posed, NOT fully excluded as a contributor.**
   As posed ("fine grids are under-trained") it is contradicted by the evidence: the pool grows
   30k / 50k / 80k across grids 4 / 8 / 16 and the net demonstrably learns **more** at G = 16
   (shrinkage 0.99 → 0.83 → 0.74). It cannot be *fully* excluded as a contributing factor and is
   not claimed to be.

---

## F6 — Provenance gap (and an explicitly INCONCLUSIVE forensic)

**MEASURED:** `_train_grid_pipeline` accepts an `imsize_set` keyword but **did not persist it**.
Confirmed against the frozen artifacts: `meta` on `npe_{4,8,16}.jld2` is
`(grid = G, n_pairs = N, use_gpu = false)` — no image-size key. **Grid 16's actual training
image-size distribution is therefore unrecorded.**

**NOT CONCLUSIVE — the zt forensic.** A comparison of the persisted summary `ZScoreTransform`
row-scales mildly favours a 512²-only training pool:

| Hypothesis | mean row-scale |
|---|---|
| persisted grid-16 `zt` | 0.4249 |
| simulated 512²-only pool | 0.4288 |
| simulated mixture | 0.405 |

**This is NOT conclusive at n = 100.** The gap between the persisted value and the 512²-only
reference (0.0039) is not large relative to the sampling noise at that n. **The grid-16 training
image-size distribution must be treated as `unknown`, and is recorded as such.** It has **not**
been back-filled with a guess, and `default_imsize_for(16)` — what a run *would* use today — is
explicitly **not** evidence of what the existing artifact was trained on.

**Mitigation shipped with this document:** `_train_grid_pipeline` now persists
`imsize_set` / `imsize_weights` / `imsize_source` into the `meta` of all three artifacts
(`npe_G` / `ratio_G` / `ood_nulls_G`), read back via
`ProteinCoLoc.training_imsize_provenance(...)`. Pre-provenance artifacts — including all three
frozen bundles — report `recorded = false` and `imsize_set = :unknown`. When a caller injects a
custom `datagen`, the keywords do not describe the pool and the provenance is likewise recorded as
`:unknown` rather than guessed.

---

## What was deliberately NOT done

- **No retraining.** Every finding above is read-only against the frozen bundles.
- **No re-run of the SBC or BF gate arms.** The recorded verdicts in
  `gate_report_{4,8,16}.jld2` are byte-unchanged.
- **No pre-registered constant touched.** `test/gate/gate_consts_{4,8,16}.jl` are frozen
  pre-registration and remain byte-identical. The two mitigations shipped alongside this document
  are a metadata addition and a reporting addition; neither changes a threshold or a verdict.
