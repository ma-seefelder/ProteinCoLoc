# Reference: SBC Calibration of the NPE (spikes 009, 010, 011, 012, 013)

**Question chain:** The amended grid-8 confirmatory gate FAILED its SBC arm — every-parameter
KS/χ² rejections. These five spikes decompose *why*, and whether any model lever fixes it. **Bottom
line: the coloc TARGETS (ρ_true, Δρ) ARE calibrated. The failures are (a) a prior-atom artifact on
ρ_true, fixable test-side, and (b) a ~0.06-SD marginal location drift on NON-identified nuisances,
sitting exactly at the M=2000 over-power edge, which no realistic model lever reliably removes.**

## The recommended calibration approach (build this)

For the SBC arm of the productionized gate:

1. **Coloc TARGETS (ρ_true, Δρ): strict point-null SBC at M=2000.** High power is *wanted* here — the
   scientific claim rests on the targets. On every DEV/replication seed tested (009–013) they are
   calibrated (atom-corrected KS 0.34–0.77 for ρ_true; Δρ 0.12–0.30).
2. **ρ_true carries prior atoms → use randomized ranks** (the standard correction for a
   discrete+continuous mixed prior). ρ_true is the **only** atom-carrying parameter. Do NOT truncate
   the prior to remove the atoms — that costs high-|ρ| inference accuracy and ADVI comparability for
   no test-side benefit (spike 013).
3. **NON-identified nuisances (spillover, autofluorescence, label_efficiency, shift_dx/dy, noise):
   nuisance-appropriate equivalence test, not a strict point null.** M=2000 point-null SBC is
   *over-powered* for parameters the data cannot inform. Spec draft:
   `.planning/phases/07-productionization-conditional-on-go/07-NUISANCE-SBC-SPEC-DRAFT.md`.

### The randomized-rank correction (recipe)

From `sources/010-calibration-vs-imsize/atom_ks_correction.jl`. `u = (rank+0.5)/(L+1)`; θ* on an
atom is forced to a rank extreme, so replace atom-case ranks with fresh uniforms:

```julia
using HypothesisTests, Distributions, Random
u = (R[:,p] .+ 0.5) ./ (L + 1)               # normalized SBC ranks for parameter p
x = theta[:,p]; lo, hi = minimum(x), maximum(x)
atom = (x .== lo) .| (x .== hi)              # θ* sitting exactly on a prior atom
# Standard correction for a mixed (discrete+continuous) prior:
Random.seed!(1234)
ur = copy(u); ur[atom] = rand(count(atom))
pvalue(ExactOneSampleKSTest(ur, Uniform(0,1)))   # ρ_true: 0.0096 (raw) → 0.46 (randomized) PASS
# equivalently, exclude atoms and report the excluded fraction:
pvalue(ExactOneSampleKSTest(u[.!atom], Uniform(0,1)))  # 0.667 PASS
```

## What each spike established

### 009 — shrinkage measures WIDTH, not LOCATION (two failure modes coexist)
- Premise "if the posterior *is* the prior (shrinkage≈1), ranks must be uniform, so a rejection is a
  bug" is **wrong**. `spillover` has shrinkage 1.010 (flagged *vacuous*) yet rejects SBC — because it
  has a **location shift** (z = +4.3). **Shrinkage = post_sd/prior_sd is a ratio of WIDTHS; a vacuous
  parameter can still legitimately reject on location.**
- Two distinct mechanisms: **location shifts** (ρ_true z=−5.2, spillover z=+4.3 — the finding-F2
  marginal-drift mechanism: NPE's forward-KL constrains nothing about ∫q(θ|Z)p(Z)dZ = π(θ)) and
  **shape failures** (label_efficiency, shift_dy, noise, |z|<3 — misshapen histograms a mean can't
  capture; mechanism was then unknown).
- **Consequence for reading every SBC table in this project:** `vacuous = shrinkage ≥ 0.95` means
  *uninformative about spread*, NOT *posterior = prior*. Both readings must be stated wherever the
  flag is used. (z uses the uniformity-null SE, anti-conservative under non-uniformity — read z as
  effect size, not a test.)

### 010 — ρ_true "overconfidence" was a PRIOR-ATOM artifact (parent README RETRACTED by its addendum)
- **Read `sources/010-calibration-vs-imsize/ADDENDUM-atoms.md` — it retracts the parent README's
  "overconfident posterior" conclusion.** The U-shaped ρ_true rank histogram is NOT a too-narrow
  posterior. `sample_prior` draws `μ* ~ Truncated(Cauchy(0,0.3), −1, 1)` and sets `ρ_true = ghat(μ*)`;
  `ghat` (frozen SIM-02) is a **clamped** piecewise-linear inverse over the realized sweep range
  `μ ∈ [−0.67976, 0.847149]`. Cauchy tails extend past that, so ~6.5% of μ draws map to **exactly
  ±0.99** — atoms. Predicted atom mass **6.76%** vs observed **6.45%** (confirmed to 3 dp).
- SBC assumes a continuous prior; θ* on an atom forces the rank to an extreme. **Excluding atoms →
  KS 0.667; randomized ranks → KS 0.460. Both PASS. ρ_true is CALIBRATED.** Counter-check vs
  overconfidence: median post_sd is 0.0654 (tail cases) vs 0.0678 (middle) — no narrowing.
- **Image size is a real covariate via out-of-support GENERALIZATION, not within-range drift.**
  `label_efficiency` at 256² carries standardized bias +1.100 — but 256² is *outside* the amended
  net's training mixture (512²–2048²). Large bias off-support is expected; it confirms the net does
  not generalize across image size. F5's concern stands; its mechanism is out-of-support failure.
- **Methodological result: the one-shot gate has ~±1.5-z seed variance** (ρ_true z −5.2→−1.69 on a
  replication). The gate numbers are not an artifact, but z=−5.2 was a high draw — state the seed
  noise wherever per-parameter numbers are quoted in a binding one-shot protocol.

### 011 — residual nuisance failure = 0.06-SD location drift at the M=2000 over-power edge
- Decomposing stored ranks (mean-centre then re-test shape): spillover, autofluorescence,
  label_efficiency are **pure location** (KS-raw rejects, mean-centred KS passes). The residual is a
  **location misplacement of the learned marginal — the flow centres a non-identified parameter's
  marginal ~0.05–0.07 posterior-SD off the prior mean.** Not overconfidence, not image size.
- **Power curve (pure Monte-Carlo, no training):** smallest drift each M resolves at 50% power —
  M=250: 0.168 SD, M=500: 0.118 SD, M=1000: 0.084 SD, **M=2000: 0.059 SD**. The observed flow error
  (0.04–0.07 SD) **straddles the M=2000 line** — it would NOT reject at M=500. **Self-correction:** an
  earlier "<0.02 SD required" claim was the 1/√M rule-of-thumb, too pessimistic; the real threshold
  is 0.059 SD, so the bar is sharp but **not impossible**.
- Two defensible remedies: capacity (halving 0.07→0.03 SD drops power ~50%→~15%; retrain, uncertain)
  or a **nuisance-appropriate equivalence test** (M=2000 is right for TARGETS; demanding 0.059-SD
  marginal fidelity on NON-identified nuisances is arguably the wrong null). The power analysis is
  **outcome-independent** — it supplies the justification §5 of the amendment deliberately refused a
  results-dependent band for.

### 012 — raising flow capacity REDISTRIBUTES the drift, does not reduce it (VERDICT: PARTIAL)
- Raised ONLY the flow (coupling 10→16, depth 2→3, width 128→256; summary net unchanged) — the
  marginal-accuracy lever, not the identifiability lever. New artifact root
  `artifacts/spike012_highcap/`; amended_v2 and frozen bundles verified **byte-identical**
  (npe_8 sha256 `198bb078…` unchanged before/after).
- **The residual drift stays in the 0.02–0.08 SD band under BOTH capacities.** Capacity moved *which*
  nuisance drifts (autofluorescence, shift_dy improved; **shift_dx got worse**) without shrinking the
  band. Consistent with the drift being a **fundamental property of a neural flow reproducing the
  marginal of a NON-identified parameter**, not a capacity shortfall. This **strengthens the
  test-side (equivalence) route.**
- **ρ_true and Δρ — the targets — stay calibrated under both capacities** (atom-corrected KS
  0.53–0.74). Targets are insensitive to the capacity lever. Overfitting appeared at this capacity
  (val risk drifted 0.55→0.59, contained by early stopping) → an even larger flow is unlikely to
  help; capacity is near its useful ceiling. **Single DEV seed** (`0xca11b3`); what's robust is the
  drift-magnitude band and the redistribution, not per-nuisance pass/fail.

### 013 — μ-prior truncation removes the atoms but costs high-|ρ| inference + ADVI comparability
- Truncating the μ-prior to `[GHAT_MU_MIN, GHAT_MU_MAX] = [−0.67976, 0.847149]` makes ρ_true
  continuous. Custom `sample_prior_trunc`/`trunc_datagen` injected via `_train_grid_pipeline`'s
  `datagen` seam — **src/ NOT edited**; new root `artifacts/spike013_trunc/`; architecture =
  amended_v2 defaults (10/2/128) so ONLY the μ-prior differs. Baselines byte-identical (`198bb078`).
- **§A Atoms: eliminated** (0/1000 vs 6.40%). **§B ρ_true SBC natively clean: KS 0.342 with NO
  correction** (all three columns identical — no atoms to correct), vs amended_v2 needing randomized
  ranks (0.0096→0.428). **§C THE COST:** high-|ρ| inference **degrades — bias roughly DOUBLES at
  ±0.99** (0.020→0.059, 0.017→0.058) because truncation removed the training mass that taught perfect
  correlation; at |ρ|≤0.90 it's equal-or-slightly-better. **§D nuisance drift unchanged** (same
  0.04–0.09 SD, merely reshuffled). Also **breaks ADVI comparability** (the CLAUDE.md constraint).
- **Decision: randomized ranks achieve the same SBC validity at NEITHER cost.** The atoms are a
  legitimate feature of the Turing prior; the textbook treatment of a mixed distribution is the
  randomized rank, not a model change. Truncation's only merit is a manuscript with no "we
  randomized atom ranks" sentence — cosmetic, against a real inference cost + the ADVI break.

## What to avoid

- **Don't call ρ_true "overconfident."** That reading was RETRACTED (010 addendum). It is calibrated;
  the U-shape is the prior atoms.
- **Don't run strict M=2000 point-null SBC on the non-identified nuisances** — it is over-powered
  (resolves 0.059 SD; the flow's irreducible marginal drift is ~0.06 SD). Use an equivalence test.
- **Don't reach for flow capacity to fix the nuisance drift** — it redistributes, doesn't reduce
  (012 PARTIAL), and overfits near its ceiling.
- **Don't truncate the μ-prior to "fix" ρ_true SBC** — it doubles high-|ρ| bias and breaks ADVI
  comparability; randomized ranks solve the same problem for free.
- **Don't confuse shrinkage with location** — a vacuous (wide) parameter can still reject on a
  location shift.
- **Don't quote one-shot per-parameter z/KS without the ~±1.5-z seed-variance caveat.**

## Constraints / honesty

- **§6.4 blocks iterating the current amended_v2 gate.** Any model change requires a *fresh
  pre-registration on a retrained model* to legitimize a new confirmatory run; all fixes must bundle
  into **one final amendment**, not incremental tweaks. Neither capacity (012) nor truncation (013)
  earns that model change. The strategic fork: take truncation as the §6.4 change anyway (accepting
  its costs), or accept that no model change is warranted and route to Go/No-Go with the existing
  FAIL, publishing with named limits.
- **Verdicts:** 009 VALIDATED, 010 VALIDATED (parent README's overconfidence claim retracted by its
  ADDENDUM — read the addendum as the conclusion), 011 VALIDATED, **012 PARTIAL**, 013 VALIDATED.
- **Single DEV seed** per SBC comparison; the ±1.5-z seed variance applies to all per-parameter
  pass/fail. Robust results are the *drift-magnitude band* and *target calibration*, not individual
  nuisance verdicts. No production seed consumed; no gate report written; src/ untouched.
- Superseded readings: the earlier `phase5-gates-failed-npe-overconfident` memory
  ("overconfident"/"clamp-artifact") was CORRECTED here — see `phase7-amended-gate-diagnosis`.

## Origin

`sources/009-shrinkage-vs-location/`, `sources/010-calibration-vs-imsize/` (README **+
ADDENDUM-atoms.md** + `atom_ks_correction.jl` + `mixture_replicate.jl` + `location_error.jl` +
`shape_analysis.jl`), `sources/011-sbc-power-nuisance/` (README + `power_curve.jl` + `power.png`),
`sources/012-capacity-lever/` (README + `train_highcap.jl` + `eval_highcap_vs_amended.jl` +
`sha256_before.txt`), `sources/013-mu-prior-truncation/` (README + `train_trunc.jl` +
`trunc_datagen.jl` + `eval_trunc_vs_amended.jl`). Spec draft:
`.planning/phases/07-productionization-conditional-on-go/07-NUISANCE-SBC-SPEC-DRAFT.md`.
