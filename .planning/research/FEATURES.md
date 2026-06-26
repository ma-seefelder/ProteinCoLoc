# Feature Research

**Domain:** Amortized, calibrated simulation-based inference (SBI) for fluorescence-microscopy colocalization (ProteinCoLoc v2.0 "AmortizedColoc" spike)
**Researched:** 2026-06-26
**Confidence:** HIGH for SBI table-stakes (multiple recent peer-reviewed sources + method docs); MEDIUM for the colocalization-specific novelty claim (grounded by absence of prior art in searches, not by a positive citation of someone failing to do it)

> "Users" here are two audiences: (1) the **author** as solo operator who must trust the spike enough to make a Go/No-Go call, and (2) **reviewers / the SBI+microscopy community** who decide whether the calibration story is credible and novel enough to publish. Table stakes = what *any* amortized SBI tool must demonstrate to be taken seriously. Differentiators = the specific combination that makes this publishable in the colocalization niche.

---

## Feature Landscape

### Table Stakes (Credibility — Reviewers Will Reject Without These)

Features the SBI community now treats as mandatory. Missing these = "you trained a network and trusted it blindly," which post-2021 reviewers explicitly distrust (the Hermans et al. "trust crisis").

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Amortized posterior (NPE) for ρ_true and Δρ** | The headline capability; the whole point is ms-per-dataset vs ADVI-minutes. Single forward pass after training is the defining property of amortized SBI. | MEDIUM | NeuralEstimators `PosteriorEstimator` + normalizing-flow conditioner. Inferring **both** ρ_true (for SBC) and Δρ (for BF comparability) is a design choice already locked in PROJECT.md. |
| **Quantitative speedup benchmark vs ADVI** | A "faster" claim without wall-clock numbers + matched RMSE/interval-width is not publishable. Amortization only pays off if accuracy holds. | LOW | Run real `colocalization()` ADVI on 20–30 stacks; report RMSE, interval width, wall-clock ratio. Target >100× (DoD). Cheap once NPE exists. |
| **Simulation-Based Calibration (SBC) rank histograms** | Talts et al. SBC is *the* standard self-consistency check for Bayesian computation. Uniform ranks ⟺ correct average coverage. Expected in every credible SBI paper. | MEDIUM | M≈2000 draws θ*~π → simulate → L posterior draws → rank of θ*. Shape diagnoses bias direction: U-shape = overconfident (too narrow), ∩-shape = underconfident, skew = biased. |
| **Uniformity hypothesis test on ranks (KS / χ²)** | Eyeballing histograms is not enough; reviewers want a p-value or gray-band test. | LOW | HypothesisTests.jl KS/χ². Per-parameter. Note the multiple-comparison caveat across params. |
| **Coverage curve (nominal vs empirical credible-interval coverage)** | The operational meaning of calibration: does the x% interval contain truth x% of the time? More interpretable than ranks to a biology audience. | MEDIUM | Falls out of the same SBC simulation loop. Overconfidence = empirical below diagonal. |
| **ECE / MCE + traffic-light verdict** | Turns the calibration diagnostics into a single honest go/amber/stop signal — the project's "Seefelder signature." Reusable from BayesInteractomics `_bin_calibration`/`CalibrationResult`. | LOW | Reuse existing template; mostly wiring. |
| **Amortized Bayes factor (NRE / Evidence Network)** | The existing tool's BF is KDE+quadgk per dataset; an amortized log-BF in one forward pass is the parallel headline to NPE. Evidence Networks (Jeffrey & Wandelt 2023/24, l-POP-Exp loss) are the current standard for fast neural model comparison. | MEDIUM-HIGH | `RatioEstimator` or dedicated evidence net. Must be **validated against** `compute_BayesFactor()` in the well-specified regime — agreement is the credibility proof. |
| **OOD / misspecification flag** | Post-2022 SBI consensus: amortized posteriors are silently wrong off-distribution. A flag that fires on misspecified inputs and stays quiet in-distribution is now expected, not optional (Schmitt et al. 2024). | MEDIUM | Summary-space density score (Mahalanobis / flow-likelihood) + posterior-predictive mismatch. Detectability depends on the summary statistic (documented caveat). |
| **Reproducibility from a single seeded script** | `demo.jl` with fixed Random123 seed. Non-reproducible SBI is uncitable. | LOW | Already a DoD item. |
| **Held-out evaluation set** | Train/test split (5k held-out) so accuracy claims aren't in-sample. Basic ML hygiene reviewers assume. | LOW | Part of AP3 data generator. |
| **Prior consistency with the existing Turing model** | The simulator prior π(θ) must match the published `@model` ranges (μ/ν/σ/τ) or the ADVI comparison is apples-to-oranges and the whole baseline collapses. | MEDIUM | Cross-read `colocalization.jl` priors into `spike/NOTES.md`. A correctness constraint, not a feature, but failure here invalidates everything. |

### Differentiators (Novelty / Publishability)

The publishable claim is **the combination**, not any single piece. Each ingredient exists in the SBI literature (astro/cosmology/MRI), but searches surface **no prior application of amortized SBI to fluorescence colocalization**, and none combining amortized BF + SBC/coverage + OOD flag in this domain.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **First amortized SBI for colocalization** | SBI is established in cosmology, X-ray spectra, diffusion MRI, exoplanets — but not colocalization microscopy (no hits in search). Bringing it to a high-volume bioimaging task with an established Bayesian baseline is the core novelty. | HIGH (the whole spike) | Novelty rests on the *domain transfer + validated baseline*, not on inventing SBI. Frame honestly: "first to apply," not "first to invent." |
| **Amortized Bayes factor that reproduces a published per-dataset BF** | Most SBI papers do posteriors *or* model comparison; few validate an amortized BF against an independent, already-published KDE-based BF on the same statistic. The agreement-in-well-specified-regime check is a strong credibility lever. | MEDIUM-HIGH | Differentiator = the *cross-validation against `compute_BayesFactor()`*, not the evidence net per se. |
| **SBC/coverage proof bundled with the amortized BF and OOD flag** | The trifecta — calibrated posterior + amortized evidence + honest misspecification flag — as one shipped story is rare even in mature SBI fields and unseen in colocalization. | MEDIUM | Integration cost, not new algorithms. This bundle is the publishable headline. |
| **Summary-statistic ablation (efficacy, not assumption)** | Measuring whether the *existing* patch-correlation vector is a sufficient statistic — vs + Manders/moments, vs DeepSet permutation-invariant pooling — turns "we reused the old summary" into a principled finding. Detectability of misspecification *depends on the summary* (literature), so ablation also strengthens the OOD claim. | MEDIUM | DeepSets (Zaheer 2017) for permutation-invariant patch pooling; compare RMSE/coverage/BF-agreement across summary sets. Start minimal (identical to ADVI pipeline) for comparability. |
| **Honest negative-result calibration reporting** | Pre-registering a threshold and reporting miscalibration as a *finding* (with traffic-light) rather than hiding it is itself differentiating and aligns with the field's reliability turn (balanced NRE, conservative posteriors). | LOW | Define pass/fail threshold before running SBC (open decision in plan). |
| **Misspecification gap made measurable, not hidden** | The simulator-realism gap is the obvious attack vector; the OOD flag converts it from a weakness into a reported, quantified property. Tested with deliberately misspecified images (double spillover, foreign PSF). | MEDIUM | This is the spike's defensive-credibility move. |

### Anti-Features (Deliberately NOT in the Spike)

These look like obvious adds but would blow the 4-week spike budget or dilute the single falsifiable claim. Most are already named in the plan's Scope-Guard.

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| **Full 4-level hierarchy (Pixel→Segment→Celltype→Neighborhood)** | The eventual scientific model; "do it properly." | Multiplies simulator + flow complexity; SBC across a deep hierarchy is a research project of its own. Kills the spike timeline. | Prove 2-channel 2D first; defer hierarchy to conditional full build-out. |
| **3D / Z-stacks** | Real data is volumetric. | 3D simulator (PSF, registration, noise) and memory/compute cost explode; not needed to prove the principle. | 2D spike; flag 3D in Go/No-Go full-build plan. |
| **GUI** | Usability. | Zero relevance to the inference/calibration claim; pure sink. | Seeded `demo.jl` script. |
| **Multi→2-channel generalization** | Flexibility. | Variable channel count breaks the fixed summary dimension the NPE needs. | Fix 2 channels in spike; generalize in production only. |
| **User-definable `num_patches` in the spike** | Configurability. | NPE needs a *fixed* summary dimension; variable grid = retrain-per-config. | Fix 8×8 in spike; expose in productionized API (post-Go). |
| **Editing `src/bayes.jl` / `colocalization.jl` / main pipeline** | "Just integrate it now." | Breaks the strict-decoupling constraint; risks the two live manuscript pipelines. | Everything in `spike/` with own Project.toml; clean cut documented for later merge. |
| **GPU/CUDA as a requirement** | "Train faster." | Windows Flux/CUDA friction; reproducibility risk; 2D is small enough for CPU. | CPU-only baseline; CUDA.jl optional accelerator that degrades gracefully. |
| **BayesFlow / Python (PythonCall) in the core path** | Mature, battle-tested SBI library. | Adds Conda/Python coupling, undermines Julia-native + Windows-tauglich goal. | NeuralEstimators.jl default; BayesFlow a *documented fallback* only if Julia flow won't converge. |
| **Wet-lab data, AlphaFold3/Boltz/Chai deps** | Realism / "validate on real images." | Out of scope for a simulation-only proof-of-principle; introduces external dependencies and confounds the calibration story. | Simulation-only; real-data validation belongs to full build-out. |
| **Sequential / focused SBI rounds (SNPE/SNRE)** | Better sample efficiency. | Breaks amortization (per-observation retraining) — directly contradicts the ms-inference claim. | Stay **amortized** (single-round NPE/NRE); sample efficiency is not the bottleneck for a 2D toy. |
| **TARP / full joint-coverage machinery** | Strongest coverage diagnostic (joint, not just marginal). | Nice-to-have; marginal SBC + coverage curve already clear the table-stakes bar for a spike. Adds scope. | Marginal SBC + coverage now; note TARP as a v2+ hardening option. |

---

## Feature Dependencies

```
[AP2 Physics simulator g(θ)→MultiChannelImage]
    └──requires──> consistent prior π(θ) matched to Turing @model ranges
    └──enables──> [AP3 Training-data generator (θ~π → summary vector, JLD2 cache)]
                       └──requires──> existing summary stat (correlation/patch/_prepare_data)
                       └──enables──> [AP4 NPE training]
                                          ├──enables──> [ADVI speedup benchmark]
                                          ├──enables──> [AP5 SBC / coverage / ECE-MCE / traffic-light]
                                          └──enables──> [AP6 NRE amortized BF] ──validated-against──> compute_BayesFactor()
                                                             └──enables──> [AP6 OOD / misspecification flag]

[Summary-statistic ablation] ──enhances──> [AP4 NPE] and ──strengthens──> [OOD flag]
                                 (detectability of misspecification depends on the summary)

[Amortized / single-round design] ──conflicts──> [Sequential SBI (SNPE/SNRE)]
[Fixed summary dimension (8×8)]   ──conflicts──> [user-definable num_patches] and [multi-channel generalization]
```

### Dependency Notes

- **Everything downstream requires the simulator (AP2) + prior consistency.** If π(θ) drifts from the Turing priors, the ADVI baseline comparison and the BF validation both become meaningless. This is the highest-leverage correctness gate.
- **SBC, BF, and OOD all consume the trained NPE/summary network.** They can share the same simulation loop (θ*~π → simulate → infer), so build the simulate→infer harness once.
- **OOD flag depends on the summary statistic choice.** Literature is explicit: misspecification is only detectable if the summary retains the discrepancy. This couples the ablation (differentiator) to the OOD flag (table stakes) — run them aware of each other.
- **Amortization conflicts with sequential SBI.** Choosing single-round NPE/NRE is what *buys* the ms-inference claim; do not drift into SNPE for accuracy without re-examining the core value proposition.
- **Fixed summary dimension conflicts with configurability.** The NPE's fixed-input requirement is why `num_patches` and multi-channel are deferred — not arbitrary scope-cutting.

---

## MVP Definition

### Launch With (the spike = falsifiable proof-of-principle)

Minimum to support a defensible Go/No-Go and a publishable claim.

- [ ] **Physics simulator → MultiChannelImage, prior-consistent** — without it nothing else is valid
- [ ] **NPE for ρ_true + Δρ, benchmarked vs ADVI (RMSE + interval width + >100× wall-clock)** — the headline capability
- [ ] **SBC rank histograms + KS/χ² + coverage curve + ECE/MCE + traffic-light** — the calibration credibility, the project's differentiator anchor
- [ ] **Amortized log-BF (NRE/Evidence Net) validated against `compute_BayesFactor()`** — the second headline
- [ ] **OOD/misspecification flag (fires on misspecified, quiet in-distribution)** — defensive credibility against the realism-gap critique
- [ ] **Seeded `demo.jl` chaining all of it; main repo provably untouched** — reproducibility + decoupling

### Add After Validation (only on a Go decision — productionization)

- [ ] **Summary-statistic ablation written up as a finding** — strengthens the paper; can be partial in spike, full at build-out — trigger: NPE works and time remains
- [ ] **Integrate NPE/NRE into `src/` as a shipped feature** — trigger: Go memo
- [ ] **User-definable `num_patches`** — trigger: productionization (after fixed-dim spike succeeds)
- [ ] **Turing→RxInfer.jl backend evaluation** — trigger: separate feasibility; decide complementary-vs-redundant with amortized SBI

### Future Consideration (full build-out, v2+)

- [ ] **4-level hierarchy, 3D/Z-stacks, multi-channel** — defer: each is a research effort; spike must prove 2D-2-channel first
- [ ] **TARP / joint-coverage diagnostics** — defer: marginal SBC suffices for the spike's calibration claim
- [ ] **Balanced/conservative NRE (Delaunoy 2022)** — defer: adopt only if SBC reveals overconfidence that simple flow-tuning can't fix
- [ ] **Real wet-lab validation** — defer: simulation-only is the deliberate spike boundary

---

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Prior-consistent physics simulator | HIGH | MEDIUM-HIGH | P1 |
| NPE (ρ_true + Δρ) | HIGH | MEDIUM | P1 |
| ADVI speedup + accuracy benchmark | HIGH | LOW | P1 |
| SBC ranks + uniformity test | HIGH | MEDIUM | P1 |
| Coverage curve + ECE/MCE + traffic-light | HIGH | LOW | P1 |
| Amortized BF (NRE/Evidence Net) vs KDE-BF | HIGH | MEDIUM-HIGH | P1 |
| OOD / misspecification flag | HIGH | MEDIUM | P1 |
| Seeded reproducible demo | HIGH | LOW | P1 |
| Summary-statistic ablation (+DeepSet) | MEDIUM | MEDIUM | P2 |
| Honest negative-result threshold (pre-registered) | MEDIUM | LOW | P2 |
| Productionize into `src/` | HIGH | HIGH | P2 (post-Go) |
| User-definable num_patches | MEDIUM | MEDIUM | P3 |
| RxInfer backend evaluation | MEDIUM | MEDIUM | P3 |
| TARP / hierarchy / 3D / real data | MEDIUM-HIGH | HIGH | P3 |

**Priority key:** P1 = must-have for the spike; P2 = strengthens publishability / post-Go; P3 = future build-out.

---

## Competitor / Prior-Art Feature Analysis

| Capability | Classical coloc tools (Coloc2, JACoP) | Existing ProteinCoLoc v1 (Turing ADVI + KDE-BF) | SBI in other domains (astro, MRI, exoplanets) | AmortizedColoc v2.0 (this) |
|------------|----------------------------------------|--------------------------------------------------|-----------------------------------------------|-----------------------------|
| Colocalization metric | Pearson, Manders M1/M2, Costes | Hierarchical Bayesian ρ posterior | n/a (different domains) | Amortized ρ_true / Δρ posterior |
| Significance test | Costes Monte-Carlo block-resampling (per image) | KDE Bayes factor on Δρ + quadgk (per dataset) | Often posterior-only, no BF | Amortized log-BF (one forward pass), validated vs v1 KDE-BF |
| Per-dataset speed | seconds–minutes (resampling) | ADVI minutes per dataset | ms after training | **ms after training** |
| Calibration proof | none (frequentist p-value only) | none reported | SBC/coverage increasingly standard | **SBC + coverage + ECE/MCE + traffic-light** |
| Misspecification handling | none | none | emerging (Schmitt 2024, Mahalanobis/MMD) | **explicit OOD flag, validated on misspecified images** |
| Amortization | no | no | yes | yes |

**Reading:** classical tools give a metric + a frequentist significance test but no calibrated uncertainty; ProteinCoLoc v1 adds Bayesian uncertainty + a BF but pays per-dataset ADVI cost and reports no calibration. SBI elsewhere has amortization and (increasingly) calibration but has not touched colocalization. v2.0's white space = **amortized + calibrated + misspecification-aware, in colocalization, validated against an existing published Bayesian baseline.**

---

## Sources

- Talts et al., *Validating Bayesian Inference Algorithms with Simulation-Based Calibration* (SBC foundation) — via sbi docs: https://sbi-dev.github.io/sbi/0.22/tutorial/13_diagnostics_simulation_based_calibration/ and https://sbi.readthedocs.io/en/stable/how_to_guide/16_sbc.html (HIGH)
- Modrák et al., *SBC Checking: The Choice of Test Quantities Shapes Sensitivity* — https://www.researchgate.net/publication/375884312 (MEDIUM)
- Jeffrey & Wandelt, *Evidence Networks: simple losses for fast, amortized, neural Bayesian model comparison* (l-POP-Exp loss), Mach. Learn. Sci. Technol. 5 015008 (2024) — https://arxiv.org/abs/2305.11241 ; https://iopscience.iop.org/article/10.1088/2632-2153/ad1a4d (HIGH)
- Schmitt et al., *Detecting Model Misspecification in Amortized Bayesian Inference with Neural Networks* (summary-space MMD/Mahalanobis) — https://www.researchgate.net/publication/381190138 (HIGH)
- *Model Misspecification in SBI — Recent Advances and Open Challenges*, ICLR Blogposts 2026 — https://iclr-blogposts.github.io/2026/blog/2026/model-misspecification-in-sbi/ (MEDIUM)
- Hermans et al., *A Trust Crisis in Simulation-Based Inference? Your Posterior Approximations Can Be Unfaithful* — https://arxiv.org/abs/2110.06581 (HIGH)
- Delaunoy et al., *Towards Reliable SBI with Balanced Neural Ratio Estimation* (BNRE, conservative posteriors) — https://arxiv.org/abs/2208.13624 (HIGH)
- TARP / expected coverage — *Calibrating Neural SBI* https://arxiv.org/pdf/2310.13402 ; *SBI: A Practical Guide* https://arxiv.org/pdf/2508.12939 (MEDIUM)
- Zaheer et al., *Deep Sets* (permutation-invariant summary networks) — https://www.inference.vc/deepsets-modeling-permutation-invariance/ (HIGH)
- Colocalization domain (Pearson/Manders/Costes, Coloc2, JACoP) — https://imagej.net/imaging/colocalization-analysis ; *Statistical tests for measures of colocalization* https://pubmed.ncbi.nlm.nih.gov/24117417/ (HIGH)
- Absence of prior art: SBI/NPE applications found in X-ray, diffusion MRI, exoplanets, cosmology — none in fluorescence colocalization (https://arxiv.org/abs/2401.06061 and related) → basis for the "first to apply" framing (MEDIUM — negative evidence)

---
*Feature research for: amortized, calibrated SBI for fluorescence colocalization (ProteinCoLoc v2.0 spike)*
*Researched: 2026-06-26*
