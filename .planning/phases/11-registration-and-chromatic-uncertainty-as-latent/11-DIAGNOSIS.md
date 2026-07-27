---
status: complete
kind: diagnosis
phase: 11-registration-and-chromatic-uncertainty-as-latent
source_blocker: 11-07-SUMMARY.md (SC1g tripwire failed, prescribed diagnosis green)
authorised_by: user decision relayed 2026-07-27 (two work items, research-lane only)
gates: none
artifacts:
  - spike/validation/p11_bottleneck_report.jld2
  - spike/validation/p11_recovery_report.jld2
  - spike/validation/p11_hybrid_report.jld2
scripts:
  - spike/validation/run_p11_bottleneck.jl
  - spike/validation/run_p11_recovery.jl
  - spike/validation/run_p11_hybrid.jl
---

# Phase 11 — Why the registration ladder is flat

## The answer, in plain language

**Is the 8×8 summary the bottleneck? YES. Confidence: high.**

The method does not look at the image directly. It chops each image into an 8×8 grid of patches
and records a single correlation number per patch — 64 numbers, plus 64 more saying which patches
had usable signal. That compressed picture is all the network ever sees.

On the real-data image size (1376×1028) each of those patches is about **172 × 128 pixels**. The
misalignment this phase is trying to detect is at most **3 pixels** — under **2 % of a patch
width**. Sliding one channel by that much barely changes any patch's correlation. By the time the
image has been reduced to those 128 numbers, the evidence about a 3-pixel shift has essentially
been averaged away.

So when the network's colocalization uncertainty fails to widen as registration gets worse, it is
not being overconfident and it is not miswired. **It is being asked to report an uncertainty it has
no information about.** No amount of retraining, extra capacity, or architecture change can recover
something that is not present in the input.

**There is also a second, independent problem**, found while producing the comparison value, and it
matters just as much: **even a perfect estimator could not have passed the test as written.** See
"The second finding" below. Fixing the summary alone would not have made the tripwire pass.

---

## What was measured

Three scripts, all research-lane, all gating nothing. Two of them contain **no trained network at
all**, so a null result from them cannot be blamed on training, capacity or optimisation.

### 1. Is the effect bigger than the noise? No — it is a fraction of it.

The D-06 pre-flight probe could detect a 3 px shift only because it uses a **paired** design: it
simulates the same scene twice from the same random key and changes only the geometry, so its
reference has *exactly zero* noise (the probe asserts this: the zero rung comes back `0.0`).

**The network never gets that paired reference.** It is handed one image per dataset and must judge
from that single draw. The noise it actually faces is the spread between two *independent* draws of
the same scene. That floor existed nowhere in any artifact — only as a claim in a source comment —
so it was measured from scratch, on the probe's own metric (`paired_l2_rows_1_64`) and its own five
frozen base thetas:

> **Independent-key noise floor: mean ‖Δs‖₂ = 0.9765 (sd 0.575)**

Against that floor, on the F5 image-size mixture:

| injected shift (px) | mean ‖Δs‖₂ | effect ÷ noise floor |
|---|---|---|
| 0.25 | 0.0938 | 0.096 |
| 1.0 | 0.1156 | 0.118 |
| 2.0 | 0.2441 | 0.250 |
| **3.0 (= `LAMBDA_MAX`)** | **0.3966** | **0.406** |
| 5.0 | 0.7063 | 0.723 |
| 8.0 | 1.2837 | 1.315 |

And the **corrected** λ test — theta fixed *excluding* the shift, shift marginalised over
`Uniform(-λ, λ)`, which is the version that is actually informative:

| λ | mean spread | spread ÷ noise floor |
|---|---|---|
| 0.25 | 0.0920 | 0.094 |
| 1.0 | 0.1067 | 0.109 |
| 2.0 | 0.1718 | 0.176 |
| 3.0 | 0.3114 | 0.319 |

**Across the entire pre-registered ladder the signal is 9 %–32 % of the per-draw noise.** It only
overtakes the noise at 8 px, far outside anything this phase claims to measure.

Two bookkeeping notes, recorded rather than smoothed over:
- The literal "does the summary move with λ at fixed theta" test is **null by construction** — λ has
  no causal path to the image except through the shift, and theta already contains the shift. It was
  run as a leak check and returned exactly `0.0`, confirming no λ leak into the forward model. It is
  not evidence about the bottleneck either way.
- The source comment in `run_p11_probe.jl` quotes an independent-key floor of 2.653 against a paired
  3 px effect of 1.281. **Those numbers did not reproduce** and are not used here; they appear to
  come from a different configuration (they match the 8 px / 256²-arm scale, not the F5 3 px scale).
  The freshly measured values above are the ones this report stands behind. At 3 px the measured
  floor-to-effect ratio is **2.46×**, which is consistent in magnitude with the comment's "~2.1×"
  claim even though the individual numbers differ.

### 2. Can any estimator extract it? No — and the control proves the test works.

Ridge regression of the shift on the summary. **Row 129 — the encoded λ — is excluded from the
predictors, and that exclusion is the entire point.** Because the training joint draws
`shift ~ Uniform(-λ, λ)`, anything that sees λ can reproduce the conditional prior spread perfectly
without extracting one bit from the image: `Uniform(-λ, λ)` has SD `λ/√3`, so the prior-SD ratio
across the ladder is `LAMBDA_MAX / LAMBDA_MIN = 3.0 / 0.25 = 12.0` **by construction**. A diagnostic
that admits row 129 measures *prior echo*, not learning. (The pool stores the raw 128-row summary
and λ in separate fields, so this exclusion is structural rather than a filter that could slip.)

30,000 training samples, penalty chosen on a held-out validation slice, evaluated on 12,500 held-out
samples, against the baseline of simply predicting the conditional prior mean (zero):

| λ bin | n | ridge RMSE | prior RMSE | ratio |
|---|---|---|---|---|
| [0.25, 0.5) | 1172 | 0.22599 | 0.21812 | **1.036** |
| [0.5, 1.0) | 2306 | 0.44997 | 0.44393 | **1.014** |
| [1.0, 1.5) | 2260 | 0.73289 | 0.73173 | **1.002** |
| [1.5, 2.0) | 2299 | 1.02480 | 1.02745 | **0.997** |
| [2.0, 2.5) | 2241 | 1.30824 | 1.31155 | **0.997** |
| [2.5, 3.0) | 2222 | 1.56502 | 1.56929 | **0.997** |

Pooled: `shift_dx` 0.998, `shift_dy` 0.998. **Ridge does not beat the prior at any rung.** In the two
smallest-λ bins it is measurably *worse* than predicting zero — the signature of fitting pure noise.

**Positive control** — the same ridge, same predictors, same split, predicting `ρ_true` instead:

> ridge RMSE **0.07412** vs prior RMSE **0.47083** → ratio **0.157**, an **84 % reduction in error**.

This is what makes the null trustworthy. The harness works, and the summary is *richly* informative
about ρ_true — which is exactly the premise of the whole method. It is specifically and only the
sub-pixel shift that is absent.

### 3. The narrow claim that survives from the earlier diagnostic

The 11-07 summary reported that the shift marginals widen 4.84–9.05 "against an ideal 12.0" and read
this as the λ conditioning being "alive and used". **12.0 is the prior-echo ceiling, not an ideal.**
It equals `LAMBDA_MAX/LAMBDA_MIN` exactly, so a network that ignores the image and echoes the λ it
was handed reproduces it precisely.

The claim that survives, and the only one that should ever be restated, is the narrow one:

> **The wiring works.** Row 129 reaches the output heads — a network ignoring it would show a
> near-constant width ratio of ~1.0, not 4.84–9.05. That is a plumbing fact and nothing more. It is
> **not** evidence of learning, and it must not be re-inflated into one.

---

## The second finding: the bar could not have been met

This came out of the analytic-propagation work and is **independent of the summary bottleneck**.

`P11_LAMBDA_ABLATION_FACTOR = 2.502` was derived as half the forward model's registration
sensitivity ratio (`lambda_ratio = 5.0039 × SC2_SPEARMAN_ATTENUATION = 0.5`). But `lambda_ratio` is
the ratio of the **registration-induced component alone**. The tripwire applies the resulting bar to
the **total Δρ posterior SD**, which also contains every other source of uncertainty.

Measured magnitudes (12 datasets, 2000 draws per arm):

| λ | net-only Δρ SD | σ_reg (forward model) | hybrid SD | hybrid ratio | fwd-model ratio |
|---|---|---|---|---|---|
| 0.25 | 0.09088 | 0.00956 | 0.09138 | 1.000 | 1.000 |
| 1.0 | 0.09131 | 0.01405 | 0.09238 | 1.011 | 1.289 |
| 2.0 | 0.09178 | 0.03111 | 0.09691 | 1.061 | 2.988 |
| 3.0 | 0.09282 | 0.04196 | 0.10186 | **1.115** | **5.004** |

**σ_reg is subdominant by roughly a factor of two at every rung.** Even accounting for the
registration term *exactly* and adding it in quadrature, the total width grows by **11 %**, not
150 %. To reach 2.502× the registration term would need to be ≈ 0.21 at λ_max — **five times larger
than the forward model says it is**.

So the tripwire demanded that a minor variance component drive a 2.5× change in a total it does not
dominate. **Even a perfect estimator with a perfect summary would have failed it.** This is a
specification defect in how the bar was derived, not a property of the trained net.

*(The hybrid column is an* **analytic correction, not amortized inference** *— no coverage
guarantee, not produced in one forward pass, not an SC2 result, not a ship path. The quadrature step
assumes independence between the net's residual uncertainty and the registration perturbation; that
assumption is stated in the artifact rather than buried. The net-only column is the only amortized
quantity in the table.)*

---

## What the summary would need

Reasoning from the measured curve in §1, not from a generic menu.

The effect crosses the single-draw noise floor at roughly **6–8 px** and sits at 41 % of it at 3 px.
To make a 3 px shift as visible as an 8 px shift is now, the shift-to-noise ratio must improve by
about **3×**. Options, in the order the measurement supports them:

1. **A finer patch grid (highest leverage, and the direct implication of the arithmetic).** The
   signal is destroyed by averaging over 172×128 px patches. At a 32×32 grid the patches are
   ~43×32 px, so 3 px becomes ~7 % of a patch instead of 1.74 % — a ~4× improvement in the ratio
   that matters, which is the right order. Cost: the summary grows from 64 to 1024 correlation rows,
   which changes `d_in` and therefore invalidates the frozen `zt`, the trained net, and the SBC
   calibration.
2. **A multi-scale summary.** Keep the 8×8 grid for ρ_true (where it demonstrably works — 84 % error
   reduction) and *add* a fine-grid block carrying the registration signal. Preserves the existing
   ρ_true performance instead of trading it away, at the cost of a wider input.
3. **An explicit correlation-vs-probe-shift curve as a summary channel.** Compute the patch
   correlation at several small probe offsets and append the curve. This is the most direct encoding
   — it measures exactly the quantity §1 shows is informative under pairing — and it is cheap in
   dimensions. It is essentially the paired probe design promoted into the summary, which is
   attractive precisely *because* the paired design is the one that can see the effect.

**Option 3 deserves the first look**: §1 shows the effect *is* detectable when a reference is
available, and a probe-shift curve manufactures that reference inside the summary instead of hoping
a single draw carries it.

**But note the second finding still applies.** Any of these fixes the *input*; none of them changes
the fact that σ_reg is subdominant in the total Δρ posterior. A redesign should be paired with a
re-derived bar that compares like with like — either gate the registration *component*, or set a
total-width bar from the measured σ_reg, not from `lambda_ratio`.

### Rough compute estimate for the deferred Option B redesign

Anchored on measured numbers from this phase: datagen **56.37 min / 50k pairs**, training
**16.08 min**, both CPU-only on 32 threads.

| Stage | Estimate | Assumption |
|---|---|---|
| Datagen (regenerate) | **~56 min**, unchanged | Simulation dominates; the summary step is a small tail. A finer grid adds correlation work but does not re-simulate. |
| Datagen if imaging cost grows | up to ~75 min | Only if a probe-shift curve requires *k* extra correlation passes; scales with *k*, not with the simulator. |
| Training | **~20–40 min** | Scales with input width. 8×8→32×32 takes `d_in` 129→1025; the first layer grows ~8×, the rest is unchanged. A multi-scale or probe-curve summary lands nearer the low end. |
| SBC / calibration re-run | **~30–60 min** | The dominant unknown. Must be re-run from scratch — a changed `d_in` invalidates every existing calibration claim. |
| **Total per iteration** | **~2–3 h** | Excludes analysis and reporting. |

Two costs that are **not** compute and should dominate the decision:
- **Every calibration claim in the milestone is scoped to the current summary.** Changing `d_in`
  invalidates the frozen `zt`, the trained net, the SBC coverage proof and the OOD flag. This is a
  re-validation project, not a retrain.
- **The shipped `amended_v2/grid_8` bundle and the Phase-7 GO are untouched by any of this** and must
  stay that way. This is research-lane only.

---

## What was NOT done

- **No ladder run.** Plans 11-08 through 11-11 remain gated; waves 7–9 were not dispatched.
- **No threshold edited.** `spike/validation/p11_consts.jl` is byte-unchanged, including
  `LAMBDA_MAX`, `LAMBDA_MIN` and `P11_LAMBDA_ABLATION_FACTOR`. The SC1g bar was **not** relaxed even
  though this report argues it is mis-specified — re-deriving it is a pre-registration decision.
- **No retraining, no re-seed, no capacity change.** `P11_ITERATION_ALLOWANCE` remains **1 of 1
  unspent**.
- **Option (i) (re-scoping SC1g onto the shift marginals) was not implemented** — rejected by the
  user as near-tautological, and §2 independently confirms the shift channel carries no image
  information, so it would have measured prior echo over an empty channel.
- **`src/`, `spike/Project.toml`, `spike/Manifest.toml` byte-unchanged.** No package installed.
- **Not verified:** that a *non-linear* estimator also fails to recover the shift. §2 is a linear
  probe. The magnitude argument in §1 is estimator-independent and the trained non-linear net showed
  the same flat response, which is why confidence is high — but the strict scope of §2 is linear
  recoverability.
- **Not verified:** any coverage or calibration property of the hybrid column. It is a comparison
  value only.

## Provenance note

The 50,000-pair pool was destroyed by worktree cleanup between plan 11-07 and this diagnosis (see
`deferred-items.md`). It was **restored, not re-drawn**: `p11_generating_config` is a pure function
of byte-unchanged constants, and sampling is counter-based Philox keyed per sample index. The
restoration was verified before use — all four realized image-size counts reproduce exactly
(512²=20069, 1024²=12604, 1376×1028=12289, 2048²=5038), 50,000 samples, `global_index` complete over
1:50000. Four exact count matches cannot arise from a different draw.

`11-07-SUMMARY.md` records the pool hash as `c12f8ab2…`; the canonical hash for this config is
`bff8550959d5…`. Since every config input is byte-unchanged, the summary's prose hash was a
mis-transcription. No constant was altered to reconcile it.
