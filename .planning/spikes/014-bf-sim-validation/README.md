---
spike: 014
name: bf-sim-validation
type: standard
validates: "Given the frozen amended_v2 grid-8 NRE, when the amortized log-BF is scored on labelled simulated coloc/null pairs (the way SBC validates a posterior — no per-pair reference), then it discriminates H1 vs H0 (AUC) and is monotone in evidence (Δρ), so the KDE-tail baseline can be replaced by a simulation-based BF check"
verdict: VALIDATED
related: [006, 007, 008]
tags: [bayes-factor, nre, simulation-based-validation, kde, calibration, proof-of-concept]
---

# Spike 014: Simulation-Based Bayes-Factor Validation

**Proof of concept — NOT a gate change.** A gate change would be a §6 pre-registration step; this
spike only demonstrates the concept on a DEV seed and reports the real numbers.

## What This Validates

**Given** the frozen amended_v2 grid-8 NRE (`artifacts/amended_v2/grid_8/ratio_8.jld2`), **when**
the amortized log-BF is scored on labelled simulated coloc (H1) / null (H0) pairs — the way SBC
validates a posterior, WITHOUT a per-pair reference — **then** it (1) discriminates H1 from H0
(ROC/AUC) and (2) is monotone in evidence strength (signed Δρ). This is the replacement for the
KDE-tail baseline, which attrits and saturates (spikes 006–008) and is semi-circular.

## The coloc/null definition (extracted — the crux)

The RatioEstimator's model index was defined at training time, `src/amortized/train_ratio.jl`:

- `:55` — `const RATIO_SPLIT_THRESHOLD = 0.0` (D-07 split on the ρ_true contrast).
- `:130-131` — for each training column, pair two **distinct** pool indices `(i1, i2)` and set
  `Z_pair = pair_encode(Zstd[:,i1], Zstd[:,i2])`, `model_index = Float32((ρ[i1] − ρ[i2]) > 0.0)`.

So a coloc/null "pair" is **two independent prior draws**, labelled by the **sign of their ρ_true
contrast**:

> **H1 (coloc) ⟺ ρ_true_sample > ρ_true_control ;  H0 (null) ⟺ ρ_true_sample ≤ ρ_true_control.**
> Evidence strength = signed **Δρ = ρ_sample − ρ_control**.

Our test pairs are constructed **identically** — two i.i.d. `sample_prior` draws per pair,
sign-of-contrast label — so the validation is exactly on the trained definition, not a guess. The
deployed `log_prior_odds = −0.0102` (measured, ≈0 by i.i.d. symmetry) is used in `amortized_log_bf`
exactly as shipped.

## Image sizes (F5 training-joint consistency)

Simulated at the **F5 mixture read from the frozen `test/gate/gate_consts_8_v2.jl`**
(`SBC_IMSIZE_SET = ((512,512),(1024,1024),(1376,1028),(2048,2048))`, weights `(0.4,0.25,0.25,0.1)`).
This was **verified equal** to the training-joint imsize provenance recorded in the amended_v2
artifacts' `meta` (`training_imsize_provenance(npe_8)`). NRE calibration holds only under the
training joint (F5, covariate shift), so simulating at the training mixture is the valid choice —
not a convenience.

## How to Run

```bash
# analysis (simulate + score); artifacts/ is gitignored
SPIKE014_M=400 SPIKE014_NPOST=1000 julia --project=. -t auto \
  .planning/spikes/014-bf-sim-validation/run_sim_bf.jl
# figure (isolated figenv; CairoMakie is not a root dep)
julia --project=.planning/spikes/figenv .planning/spikes/014-bf-sim-validation/figures.jl
```

## Seed discipline

DEV seed **`0x00BF5014`** ("BF-Sim-014"), `@assert`-ed disjoint from the full forbidden set:
`DEFAULT_MASTER_SEED` (0x1), `NPE_MASTER_SEED` (0xC0FFEE), `VAL_MASTER_SEED` (0x5BC0FFEE),
`RATIO_PAIR_SEED` (0x4A7107), every enumerated prior DEV seed, and **every v1 `PROD_SEED[G]` and
v2 `PROD_SEED_V2[G]`** (read out of the frozen `gate_consts_8_v2.jl` loaded into an isolated
module — no gate run, file unmodified). PROD_SEED_V2 is never consumed. The NRE log-BF and the
AUC/monotonicity result are **deterministic** functions of this seed (the ratio net's `logratio`
is a forward pass, no sampling); the global RNG is pinned only for the NPE draws feeding the KDE
contrast.

## Investigation Trail

1. **Extracted the coloc/null construction from the trainer** (`train_ratio.jl:107,131,55`), the
   stated crux — the label is `sign(Δρ_true)` on two i.i.d. prior draws, threshold 0.
2. **Verified the image-size joint**: read the F5 mixture from the frozen amendment and asserted it
   equals the net's recorded training provenance (both nets carry
   `imsize_set=((512,512)…(2048,2048))`, `w=(0.4,0.25,0.25,0.1)`).
3. **Simulated 400 pairs** (800 i.i.d. members, Philox-per-index ⇒ threaded & byte-reproducible)
   through the deployed `sample_prior → sample_imsize → simulate_pair → patch_summary(·,8) →
   encode_d01 → frozen-zt standardize` path; label = `sign(ρ_s − ρ_c)`.
4. **Scored each pair three ways on the SAME inputs**: amortized NRE log-BF (`amortized_log_bf` on
   `pair_encode`, one forward pass), the src non-clamped KDE baseline (`kde_log_bf_unclamped`), and
   a clamped KDE variant (1e-8 floor → ±18.42 rail).
5. **Measured** AUC(H1 vs H0), Spearman(logBF, Δρ) + binned-median curve, decision error rates at
   the logBF>0 threshold, reliability of p̂(coloc), and the KDE baseline's non-finite/saturation
   fractions on the identical pairs.

## Results

**VERDICT: VALIDATED.** M = 400 pairs, DEV seed `0x00BF5014`, N_post = 1000, F5 imsize mixture.

### 1 — Discrimination (no per-pair reference)

| statistic | value |
|---|---|
| **AUC(NRE log-BF, H1 vs H0)** | **0.9939** |
| AUC(KDE baseline, finite pairs only) | 0.9822 on **261/400** pairs (survivorship-biased subset) |
| NRE non-finite fraction | **0.0** — keeps **all** 400 pairs |

The NRE separates coloc from null almost perfectly and scores **every** pair. The KDE baseline can
only be evaluated on the 261 finite pairs — and those are the *easy* (low-|Δρ|) ones (see §3), so
its 0.982 is measured on a biased, easier subsample (spike 007 survivorship effect).

### 2 — Monotonicity in evidence Δρ

| statistic | value |
|---|---|
| **Spearman(NRE log-BF, Δρ)** | **0.9256** |
| binned-median log-BF vs Δρ | monotone across the meaningful range; crosses 0 at Δρ≈0 |

Binned medians (9 equal-count Δρ bins), Δρ-center → median log-BF:

```
−1.02→−7.38  −0.70→−7.35  −0.40→−7.33  −0.18→−3.56  +0.01→+0.26
+0.20→+4.87  +0.42→+8.73  +0.72→+8.79  +1.13→+8.76
```

The curve rises monotonically from a −7.4-nat null plateau through **≈0 exactly at Δρ≈0** to a
+8.8-nat coloc plateau. The only non-monotone step is the last bin (8.79 → 8.76, a **0.024-nat**
wiggle on the saturated high-Δρ plateau) — negligible tail noise where the bounded ratio net has
already saturated, not a calibration defect. The decision boundary `logBF>0` therefore coincides
with the true evidence boundary `Δρ>0`.

### 3 — Decision calibration

| statistic | value |
|---|---|
| accuracy @ logBF>0 | 0.960 |
| FPR (declare H1 \| true H0) | 0.0355 |
| FNR (declare H0 \| true H1) | 0.0443 |
| reliability ECE of p̂(coloc) | **0.0188** |

At the natural `logBF>0` threshold the error rates are low and near-symmetric, and predicted
P(coloc)=σ(log-BF) tracks the empirical coloc fraction (ECE 0.019). The residual FPR/FNR are
concentrated at |Δρ|≈0 pairs, where H0 and H1 are genuinely indistinguishable — the correct place
for errors.

### 4 — KDE-baseline pathology on the SAME pairs (what the sim route sidesteps)

| quantity | value |
|---|---|
| **non-finite (attrited by the gate's isfinite filter)** | **0.3475 (139/400)** |
| saturated at the ±18.42 rail (clamped variant) | 0.325 (130/400) |
| median \|Δρ\| of **dropped** pairs | **0.736** |
| median \|Δρ\| of **kept** pairs | **0.283** |

The non-clamped KDE baseline is **non-finite on 35% of the pairs** and, with the old 1e-8 clamp,
**saturates on 33%** — and precisely on the **high-|Δρ| (high-signal) pairs** (dropped median
0.74 vs kept 0.28), reproducing spike 006's mechanism and spike 007's survivorship bias directly.
The simulation-based NRE check has **no per-pair KDE**, so it keeps all 400 pairs and scores the
high-signal cases the baseline is forced to discard.

### Figure

`sim_bf_validation.png` — (A) NRE ROC (AUC 0.994); (B) NRE log-BF vs signed Δρ with binned-median
overlay and the Δρ=0 / logBF=0 decision lines; (C) reliability of p̂(coloc); (D) |Δρ| histogram of
KDE-dropped vs KDE-kept pairs (the baseline discards the high-signal tail).

## What This Establishes / Does Not Establish

- **Establishes** that the amortized NRE log-BF can be validated *without* a per-pair reference:
  it discriminates coloc/null (AUC 0.994) and is monotone in evidence (Spearman 0.926) with good
  decision calibration (ECE 0.019), on the exact training-consistent coloc/null definition and the
  training-joint image-size mixture. The KDE baseline's 35% attrition / 33% saturation on the same
  pairs — concentrated on the high-signal tail — is exactly what the simulation-based route avoids.
- **Does not** re-register the gate: replacing `BF_CORR_MIN` / `max|Δ logBF|` (KDE-referenced) with
  an AUC + monotonicity + decision-calibration criterion is a §6 pre-registration change, out of
  scope here. This spike is the outcome-independent *evidence* such a proposal would cite.
- **Single-seed caveat** (CONVENTIONS 009–013): one DEV seed; per-cell numbers carry seed noise.
  The reportable results are the **magnitudes and directions** (near-1 AUC, strong positive
  Spearman, ~⅓ KDE attrition on the high-signal tail), not any third-decimal value.
