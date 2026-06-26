# Pitfalls Research

**Domain:** Amortized simulation-based inference (NPE/NRE) on a 2-channel 2D fluorescence-colocalization physics simulator, with SBC calibration proof and OOD/misspecification flag (ProteinCoLoc v2.0 "AmortizedColoc" spike — solo, CPU-first, Windows, Julia/NeuralEstimators.jl)
**Researched:** 2026-06-26
**Confidence:** HIGH for SBI/SBC/misspecification methodology (Context7-equivalent: arXiv primary sources + NeuralEstimators docs); MEDIUM for NeuralEstimators.jl API specifics (v0.2.0, pre-1.0, actively churning); MEDIUM for Windows/Flux specifics (general ecosystem signals, not a project-specific repro)

> Scope note: these are pitfalls *specific to this spike*, mapped to the AP1–AP7 / M0–M5 structure in `.planning/10_plan_amortizedcoloc.md`. Generic ML advice (use a validation set, normalize inputs) is included only where this project has a concrete, non-obvious way to get it wrong.

---

## Critical Pitfalls

### Pitfall 1: The OOD/misspecification flag silently fails on the very misspecifications it should catch

**What goes wrong:**
The plan computes the OOD score from the **fixed, pre-existing patch-correlation summary** (Mahalanobis distance / flow-likelihood on that vector — AP6). A misspecification is only detectable if it *moves the summary distribution*. The headline tests (double spillover, foreign PSF) may be partially or wholly **orthogonal** to a patch-correlation vector: spillover and a wider PSF both inflate inter-channel correlation in ways the simulator's own `spillover`/`PSF` nuisance parameters already span, so the misspecified summary can land squarely inside the training summary cloud. Result: a green "in-distribution" flag on data that is genuinely out-of-model, and a confidently wrong posterior — the textbook amortized-SBI silent failure.

**Why it happens:**
Schmitt et al. (2024) show a hard **detection-vs-inference trade-off**: detection power exists only along the directions the summary statistic encodes. Misspecifications that are "summary-sufficient" (leave the summary distribution unchanged) are *provably undetectable* by any summary-space distance. The published misspecification detector (BayesFlow) gets around this by *learning* the summary network with an MMD-to-N(0,I) regularizer so the summary space is structured for detection. This spike deliberately uses a **fixed, hand-engineered summary not trained for detection**, so it inherits the failure mode without the mitigation.

**How to avoid:**
- Treat OOD validation as an *experiment with a positive and a negative control*, not a checkbox. Build misspecifications along **multiple distinct axes**: (a) summary-shifting (foreign/heavy-tailed noise, gross intensity scaling) — should fire; (b) **summary-orthogonal** (a perturbation chosen to leave patch-corr moments ~unchanged) — measure whether it (correctly) does *not* fire, and report this as a known blind spot rather than hiding it.
- Add a **posterior-predictive-check (PPC) channel** to the flag (the plan already lists "Posterior-Predictive-Mismatch"): simulate from the inferred posterior, recompute the summary, compare to the observed summary. PPC catches some misspecifications the marginal summary-distance misses.
- Consider augmenting the summary with moments/Manders (the AP4 ablation) *partly to widen detection coverage*, not only for inference accuracy — note this couples Pitfall 1 and Pitfall 3.
- Quantify detection: report ROC/AUC of the flag over a grid of misspecification strengths, not a single binary "it fired."

**Warning signs:**
Flag fires on trivial cases (huge intensity offset) but the AUC near realistic misspecification strengths is ~0.5; OOD score distribution for "foreign PSF" overlaps the in-distribution score histogram; posterior stays narrow and confident on misspecified inputs.

**Phase to address:** AP6 (M4). Design the OOD validation grid *before* tuning the detector; pre-register which axes are expected to be blind spots.

---

### Pitfall 2: "Tune until calibrated" becomes data-snooping — SBC overfit to the test that is supposed to certify it

**What goes wrong:**
The stated policy (PROJECT.md, Key Decisions; AP5) is **"SBC miscalibration → tune until calibrated."** If you repeatedly (a) run SBC on the *same* simulation draws, (b) read the rank histogram / KS p-value, (c) adjust flow depth, training size, or learning rate, and (d) re-test on the same draws, you are optimizing the model *against the calibration test itself*. KS p>0.05 then certifies nothing — you have selected the one configuration that happens to pass on that particular finite sample. The published claim ("demonstrated SBC/coverage proof") becomes circular.

**Why it happens:**
SBC is a hypothesis test; running it as a tuning signal turns its p-value into a selection statistic. With enough tuning iterations on a fixed M, some configuration passes by chance (multiple-comparisons-over-time). This is the calibration analogue of training on the test set.

**How to avoid:**
- **Pre-register the acceptance threshold and M** before looking at any ranks (the plan's open decision "Schwelle vorab definieren" — make this binding). Decide M (e.g. 2000), bins, the uniformity statistic, and the pass criterion up front.
- **Separate the tuning SBC from the certifying SBC.** Tune against an SBC run on one seed/draw set; report the *final, untouched* number from a **fresh, independently-seeded SBC run with new simulations** that was never used for any decision. Only the fresh run goes in the memo.
- Cap the number of tuning iterations and log every one (config + result) so the effective multiple-testing is visible and the final p-value can be honestly discounted.
- Embrace the plan's own honest framing: **miscalibration is a publishable finding, not a failure.** A clean "calibrated after N pre-registered tuning steps, confirmed on held-out SBC" is far stronger than a tuned-to-pass p-value.

**Warning signs:**
KS p-value hovers just above 0.05; small changes flip pass/fail; you have re-run SBC >5 times on the same draws; the config that passes has no mechanistic reason to be better, it was just the one that passed.

**Phase to address:** AP5 (M3). This is the #1 honesty risk of the whole spike and the quality gate calls it out explicitly.

---

### Pitfall 3: Summary-statistic insufficiency — patch-correlation throws away information the posterior needs

**What goes wrong:**
NPE/NRE are only as good as the summary they consume. The patch-correlation vector was designed for the *existing* ADVI pipeline, not as a sufficient statistic for ρ_true under the full nuisance parameterization (spillover, autofluorescence, label efficiency, sub-pixel shift, noise). If patch-corr is insufficient for these nuisances, the NPE posterior for ρ_true will be **biased or inflated**, the speedup claim becomes "fast but wrong," and — subtly — **SBC can still pass** because SBC certifies calibration of the *posterior given the summary*, not recovery of θ from the full image. You get a calibrated posterior over a statistic that has discarded the signal.

**Why it happens:**
Sufficiency is assumed rather than measured. The patch-corr vector marginalizes over spatial structure; spillover and PSF effects can be partially confounded with ρ_true in correlation-space, so the data needed to disentangle them never reaches the network.

**How to avoid:**
- The AP4 ablation is the right instrument — but make it **diagnostic, not decorative**: compare minimal patch-corr vs. + Manders/moments vs. (optionally) a DeepSet over patches, scoring each on RMSE *and* interval width vs. the ADVI baseline, per parameter. A summary is "insufficient" if adding statistics measurably tightens/de-biases the posterior.
- Use a **permutation-invariant DeepSet** over patches rather than a fixed-order MLP so adding patch-level features doesn't blow up dimensionality or impose spurious ordering (NeuralEstimators supports DeepSet summary networks).
- Cross-check sufficiency *independently of SBC*: posterior-predictive coverage of held-out **summaries** AND posterior RMSE against known θ on simulated held-out data. SBC-pass + high-RMSE = insufficiency smoking gun.

**Warning signs:**
NPE RMSE for ρ_true plateaus well above ADVI's even with more training data; adding moments/Manders sharply improves RMSE (proves the minimal summary was lossy); posterior intervals are calibrated (SBC ok) yet systematically wide or off-center vs. ADVI.

**Phase to address:** AP4 (M2), with the ablation; revisited at AP6 because summary choice also drives OOD detection power (Pitfall 1).

---

### Pitfall 4: NeuralEstimators.jl / Flux API churn and the Flux→Lux/Reactant ecosystem shift breaks the spike mid-flight

**What goes wrong:**
NeuralEstimators.jl is **v0.2.0 (pre-1.0)** — semver gives no API-stability guarantee, and minor releases can break `PosteriorEstimator`/`RatioEstimator`/`train`/`sampleposterior` signatures. The broader Julia ML ecosystem is mid-migration from **Flux.jl to Lux.jl + Reactant.jl** (SciML has fully dropped Flux in docs; NeuralEstimators' own examples now show `device = reactant_device()`). The plan pins the spike to the **Flux** backend. Risk: the Flux path becomes second-class or undocumented, examples assume Reactant, and a `Pkg.update` silently moves you to an incompatible combination — costing days of a 4-week budget.

**Why it happens:**
Solo dev on a fast-moving pre-1.0 package; copy-pasting docs examples that target a newer backend than the pinned one; not freezing the Manifest.

**How to avoid:**
- **AP1 smoke test is the gate — keep it that way.** Do not write AP2+ until `00_smoke.jl` (1-param Gaussian, `PosteriorEstimator` + NormalisingFlow, `train`+`sampleposterior`) is green on *this* machine with *this* Manifest.
- **Commit the `spike/Manifest.toml`** with pinned versions (`NeuralEstimators`, `Flux`, `Reactant` if pulled transitively) and never `Pkg.update` without re-running the smoke test. Reproducibility constraint demands this anyway.
- Decide the backend explicitly: if NeuralEstimators v0.2 examples lean on Reactant, evaluate whether the **CPU-only + Flux** path you need is actually exercised by the package's own tests; if not, that is a real risk to surface in the stack-decision note.
- Keep the documented fallbacks warm but *out of the critical path*: NormalizingFlows.jl (TuringLang, Bijectors-compatible) for a pure-Julia flow, and BayesFlow-via-PythonCall/CondaPkg for the nuclear option. Write a 10-line note on the switch cost for each *now*, not when blocked.

**Warning signs:**
Smoke test passes only after pinning to an older NeuralEstimators; docs examples won't run without `reactant_device()`; `train` errors mentioning Lux/Reactant types while you pass Flux models; transitive `Reactant_jll`/XLA artifacts fail to build on Windows.

**Phase to address:** AP1 (M0). This pitfall is precisely why the smoke-gate exists; enforce Manifest pinning as part of M0's definition of done.

---

### Pitfall 5: Windows + CUDA/Reactant/XLA toolchain breakage burns time the spike doesn't have

**What goes wrong:**
Windows is the most fragile Julia GPU/AI-compiler target. CUDA.jl, and especially the newer Reactant.jl/XLA artifact path, have a history of Windows build/runtime friction. If GPU acceleration is reached for opportunistically (the constraint says GPU is an *optional* accelerator), a half-configured CUDA or a Windows-incompatible XLA artifact can break the *whole* run rather than degrading to CPU.

**Why it happens:**
GPU code paths assume a working CUDA/XLA stack; "optional accelerator" wiring often hard-fails instead of falling back; Windows artifact coverage lags Linux.

**How to avoid:**
- Make **CPU-only the literal default and the reproducible baseline** (already the constraint). The 2D spike is small (50–200k pairs, "Training in Minuten" per the plan) — CPU is sufficient; do not block on GPU.
- Wire GPU as a **guarded, fall-through option**: `try` the GPU device, `catch` → log + continue on CPU. Never let GPU selection be load-bearing.
- `spike/demo.jl` (the reproducibility artifact) must run **CPU-only, no CUDA import required**. Verify the Manifest doesn't *force* CUDA/Reactant artifacts onto the critical path.
- If you do test GPU, do it in a *separate* script after the CPU path is proven, so a Windows GPU failure can't invalidate the demo.

**Warning signs:**
`CUDA.functional()` false but code proceeds into GPU calls; XLA/Reactant artifact build errors on Windows; demo.jl import-fails on a machine without a GPU; runtime asks for a CUDA driver during what should be a CPU run.

**Phase to address:** AP1 (M0) for the CPU-default wiring; revisit only at full build-out (plan defers GPU explicitly).

---

### Pitfall 6: SBC certifies the prior-predictive (well-specified) world, then gets quoted as if it certifies real data

**What goes wrong:**
SBC draws θ*~π, simulates, infers, and checks rank uniformity — *entirely inside the simulator*. A perfect SBC pass proves the NPE is self-consistent **under the model**, and says **nothing** about performance on real microscopy images where the simulator is wrong (the simulator-realism gap). Quoting "calibrated (SBC proven)" as a guarantee for real data conflates the well-specified and misspecified regimes — exactly the conflation flagged in the research question, and the same point as "coverage is not enough" (Anau Montel et al.).

**Why it happens:**
SBC and the OOD flag answer different questions; it's tempting to let the strong SBC result imply real-world trustworthiness. The spike is simulation-only (no wet-lab), so there is no real-data check to puncture the illusion.

**How to avoid:**
- State the scope explicitly in the memo: **SBC = calibration under π and g; the OOD flag (Pitfall 1) is the *only* line of defense for the realism gap.** Pair every SBC claim with the OOD result; never report SBC alone as fitness-for-real-data.
- Run SBC under at least one **deliberately perturbed simulator** to show what miscalibration *looks like* (sanity that the diagnostic can fail), strengthening the honest framing.
- Make the realism gap a measured deliverable, not a caveat: this is the plan's stated spike value ("OOD/Misspez-Meldung macht den Gap messbar").

**Warning signs:**
Memo language like "calibrated posteriors for colocalization" without "under the simulator"; reviewers/collaborators inferring real-data guarantees; no perturbed-simulator SBC to contextualize the pass.

**Phase to address:** AP5 + AP7 (M3, M5). It's a framing/reporting pitfall — guard it in the Go/No-Go memo.

---

### Pitfall 7: Amortized Bayes-factor (NRE) miscalibration and prior sensitivity vs. the KDE/quadgk baseline

**What goes wrong:**
The amortized log-BF (AP6, RatioEstimator / evidence network) can be **systematically over/under-confident**: standard NRE trained with binary cross-entropy is known to produce overconfident ratios, and **sigmoid saturation** makes BCE numerically unreliable when the BF spans many orders of magnitude — exactly the regime of a colocalization-vs-null test. Separately, a Bayes factor is **intrinsically prior-dependent** (it integrates the likelihood over the prior); if the simulator prior π(θ) differs in width from the implicit prior of the KDE/quadgk baseline (`compute_BayesFactor()` on Δρ), the two BFs will disagree *for correct reasons*, and the disagreement can be misread as an NRE bug (or vice-versa, a real bug excused as "prior difference").

**Why it happens:**
NRE overconfidence is a documented default-loss failure; BF prior-sensitivity is fundamental Bayesian model comparison; the validation target (KDE-BF) and the NRE-BF must share the *same* Δρ definition and *same* effective prior to be comparable, which is easy to violate silently.

**How to avoid:**
- Use the **l-POP / exponential-style loss (Jeffrey & Wandelt 2023)** the plan already specifies — it maps directly to log-BF and avoids sigmoid saturation. Avoid vanilla BCE for the ratio.
- Consider **Balanced NRE (BNRE)** if ratios look overconfident: it yields conservative ratios with the same Bayes-optimal solution and is a small loss-side change.
- **Pin the comparison:** validate NRE-BF against `compute_BayesFactor()` *on the same Δρ and the same prior ranges* (PROJECT.md constraint already requires π(θ) consistent with the Turing priors — extend this to the BF baseline). Report agreement only in the well-specified regime (DoD already scopes this).
- **Calibrate the ratio itself**: SBC/coverage the NRE posterior, or check that the classifier is calibrated (reliability diagram of the discriminator), before trusting the log-BF magnitude.

**Warning signs:**
NRE-BF and KDE-BF agree in sign but diverge by orders of magnitude; log-BF saturates / clips at extreme values; ratio reliability diagram is off-diagonal; changing prior width changes NRE-BF but you assumed it shouldn't.

**Phase to address:** AP6 (M4).

---

### Pitfall 8: Standardization / normalization leakage between training and held-out (and SBC) sets

**What goes wrong:**
AP3 standardizes the summary vectors and holds out 5k for evaluation; AP5 generates *fresh* SBC draws. If the standardization statistics (mean/std, or the Mahalanobis covariance used for the OOD score) are computed over the **full pool including held-out / SBC data**, evaluation is contaminated: held-out RMSE looks too good, and — critically — the **OOD covariance learns about the very out-of-distribution data it's meant to flag**, deflating the flag.

**Why it happens:**
Convenience: fit the scaler once on everything. The OOD case is especially insidious because the covariance is itself a learned quantity that leaks.

**How to avoid:**
- Fit **all** preprocessing (summary standardization, OOD covariance/flow) on the **training split only**; apply frozen transforms to held-out, SBC, and OOD inputs.
- For the OOD detector, the in-distribution reference statistics must come from a **held-out in-distribution validation set never seen by the detector fit**, and misspecified test images must be *fully external*.
- Encode the split discipline in the data-loader (AP3) so it's structural, not a per-script convention.

**Warning signs:**
Held-out RMSE ≈ training RMSE with suspiciously small gap; OOD score on known-misspecified inputs is lower than expected; changing what's "in the pool" changes evaluation numbers.

**Phase to address:** AP3 (data generator) sets the discipline; verified in AP4/AP5/AP6.

---

## Technical Debt Patterns

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Start with the minimal patch-corr summary only (skip the ablation) | Fastest first NPE; identical to ADVI input | May ship a biased/insufficient posterior that SBC won't catch (Pitfall 3) and a blind OOD flag (Pitfall 1) | Acceptable as the *first* loop; the AP4 ablation is mandatory before the memo |
| Single SBC run, M≈2000, read p-value once | Fast calibration verdict | Underpowered for subtle miscalibration; tempts re-runs → snooping (Pitfall 2) | Acceptable only with pre-registered M + ECDF + a held-out confirmatory run |
| `Pkg.add` latest, no Manifest pin | Newest features/bugfixes | Reproducibility breaks; API churn (Pitfall 4) | Never for this spike — reproducibility is a hard constraint |
| GPU as default for "speed" | Faster training | Windows/XLA fragility breaks the whole run (Pitfall 5) | Never as default; only as guarded opt-in after CPU path proven |
| One global scaler over all data | One line, simple | Train/eval/OOD leakage (Pitfall 8) | Never |
| Report SBC pass as "calibrated for colocalization" | Strong-sounding claim | Conflates well-specified vs real (Pitfall 6); reputational risk on a published claim | Never — always pair with OOD + "under the simulator" |

## Integration Gotchas

| Integration | Common Mistake | Correct Approach |
|-------------|----------------|------------------|
| Existing `correlation()`/`patch()`/`_prepare_data()` as NPE input | Re-implementing the summary in the spike, drifting from the ADVI pipeline → invalidates the "same input" comparability claim | Call the *existing* functions unchanged (strict decoupling); the whole comparability argument rests on identical summaries |
| `MultiChannelImage` as simulator output | Producing arrays the spike massages, so `correlation()`/`patch()` silently take a different code path | `simulate_pair(θ)` must return a real `MultiChannelImage` that `correlation`/`patch`/`otsu_thresholds` consume *unmodified* (plan AP2) |
| `compute_BayesFactor()` (KDE+quadgk) as BF baseline | Comparing NRE-BF to a baseline with a different Δρ definition or different prior | Lock identical Δρ and prior; compare only in well-specified regime (Pitfall 7) |
| Turing `@model` priors → simulator π(θ) | Simulator priors drift from μ/ν/σ/τ ranges → ADVI vs NPE no longer comparable | Copy prior ranges into `spike/NOTES.md` from the `@model` and assert them in `sample_prior()` (plan Tag-1 step 4) |
| BayesInteractomics `_bin_calibration`/`CalibrationResult` as ECE/MCE template | Copying code that edits or imports from a sibling repo → couples projects | Use as a *template* re-implemented inside `spike/`; no cross-repo import (decoupling constraint) |

## Performance Traps

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| Regenerating training data every run | Slow iteration; non-reproducible | JLD2 cache keyed by seed + π version (plan AP3) | As soon as you iterate on the network |
| Under-training the flow then blaming the summary | NPE RMSE high, looks like insufficiency | Confirm convergence (val loss plateau) *before* concluding Pitfall 3 | Conflated diagnosis wastes the AP4 budget |
| Too-small training N for a deeper flow | Calibration fails / overfits | Scale 50k→200k as M2 accuracy demands (plan); tie flow capacity to N | When you grow flow depth to chase SBC pass |
| SBC with M too small | Noisy rank histogram, false pass/fail | Pre-register M (≥~2000); use ECDF, not just histogram bars | Subtle miscalibration regime |

## Security Mistakes

Not a meaningful axis for a solo, simulation-only, offline research spike (no untrusted input, no network service, no PII). The nearest analogue is **scientific integrity**, covered under Pitfalls 2 and 6 (calibration honesty) — the real "security" risk here is publishing an overstated calibration claim.

## UX Pitfalls

The "user" is the author + future readers of the Go/No-Go memo. Domain-relevant pitfalls:

| Pitfall | User Impact | Better Approach |
|---------|-------------|-----------------|
| Reporting a single KS p-value as the calibration result | Reader can't judge robustness; hides snooping | Show ECDF + rank histogram + coverage curve + traffic-light per parameter |
| Binary "OOD flag fired: yes" | Hides the blind-spot structure (Pitfall 1) | Report OOD score distributions + ROC over misspecification strengths, name the blind spots |
| Speedup headline without accuracy caveat | "100× faster" read as "100× better" | Always pair wall-clock speedup with RMSE/interval-width vs ADVI |

## "Looks Done But Isn't" Checklist

- [ ] **NPE trained:** Often missing the *sufficiency* check — verify RMSE vs ADVI per parameter AND the AP4 ablation, not just "it samples."
- [ ] **SBC passed:** Often missing the held-out confirmatory run — verify the reported p-value comes from a fresh, never-tuned-against SBC (Pitfall 2).
- [ ] **OOD flag works:** Often missing the negative control — verify it *stays quiet in-distribution* AND has a measured ROC, AND a documented blind-spot axis (Pitfall 1).
- [ ] **Amortized BF validated:** Often missing prior/Δρ alignment — verify NRE-BF and KDE-BF share identical Δρ and prior before comparing (Pitfall 7).
- [ ] **Reproducible demo:** Often missing seed coverage — verify Random123 seeds *every* stochastic step (sim, training init, SBC draws) and that `demo.jl` runs CPU-only with a pinned Manifest.
- [ ] **Main repo untouched:** Often missing proof — verify `git status` on `src/`, `bayes.jl`, `colocalization.jl` is clean (hard decoupling constraint).

## Recovery Strategies

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| OOD flag blind to a misspecification (P1) | MEDIUM | Add PPC channel; augment summary with moments/Manders to span the missed axis; re-measure ROC; document residual blind spots honestly |
| SBC snooping discovered (P2) | LOW–MEDIUM | Freeze current config; run one fresh independently-seeded SBC; report *that* number only; log prior tuning as multiple-testing context |
| Summary insufficiency (P3) | MEDIUM | Switch summary net to DeepSet + moments; retrain; re-run ablation + SBC + OOD |
| NeuralEstimators API break (P4) | LOW (if pinned) / HIGH (if not) | Restore pinned Manifest; if package is unusable, switch to NormalizingFlows.jl, else BayesFlow-via-PythonCall (warm fallback note written in AP1) |
| Windows GPU breakage (P5) | LOW | Drop to CPU-only (already the baseline); GPU is optional |
| NRE-BF miscalibrated (P7) | LOW–MEDIUM | Swap BCE→l-POP/exponential loss; try BNRE; recalibrate discriminator; re-validate vs KDE in well-specified regime |

## RxInfer migration pitfalls (only if the backend-evaluation track is attempted)

Flagged separately because it is an *open consideration*, not spike-critical (PROJECT.md "Backend evaluation").

- **Non-conjugate Student-t hierarchy:** RxInfer's reactive message passing is fastest with conjugate/exponential-family structure. The hierarchical **Student-t** likelihood (heavy tails, ν parameter) is non-conjugate; naive translation forces approximate/iterative messages or hybrid CVI nodes, which can be slow or fail to converge — potentially *negating* the speedup that motivated the migration. Warning sign: messages requiring custom approximations for the Student-t / ν node. Prevention: prototype the Student-t node in isolation before committing; benchmark convergence, not just per-iteration speed.
- **Message-passing non-convergence:** loopy factor graphs from the hierarchy may not converge or may converge to a different fixed point than ADVI. Validate RxInfer posteriors against the *existing ADVI* on identical data before trusting them — and note the irony that this whole spike may make the RxInfer track **redundant** (amortized SBI replaces per-dataset inference entirely). Prevention: decide complementary-vs-redundant *before* sinking time into RxInfer.
- **Phase to address:** out-of-spike; a separate feasibility task. Do not let it consume the 4-week spike budget.

## Pitfall-to-Phase Mapping

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| P1 OOD silent failure | AP6 / M4 | ROC over misspecification grid; positive + summary-orthogonal negative controls |
| P2 SBC snooping | AP5 / M3 | Pre-registered M+threshold; fresh held-out SBC for the reported number; tuning log |
| P3 Summary insufficiency | AP4 / M2 | Ablation table (RMSE + interval width per param); SBC-pass-but-high-RMSE check |
| P4 NeuralEstimators/Flux churn | AP1 / M0 | Green smoke test + committed pinned `spike/Manifest.toml` |
| P5 Windows/GPU breakage | AP1 / M0 | `demo.jl` runs CPU-only, no CUDA import; guarded GPU fallback |
| P6 SBC↔real conflation | AP5+AP7 / M3,M5 | Memo states "calibrated under the simulator"; SBC always paired with OOD |
| P7 NRE-BF miscalibration / prior sensitivity | AP6 / M4 | l-POP loss; identical Δρ+prior vs KDE-BF; discriminator reliability diagram |
| P8 Standardization leakage | AP3 / M1–M2 | Scaler/covariance fit on train split only; external misspecified test images |

## Sources

- Schmitt, Bürkner, Köthe, Radev (2024). *Detecting Model Misspecification in Amortized Bayesian Inference with Neural Networks: An Extended Investigation* — summary-space MMD detector, detection-vs-inference trade-off, summary-sufficient (undetectable) misspecifications, low-N power loss. https://arxiv.org/html/2406.03154v2 [HIGH]
- Model Misspecification in SBI — Recent Advances and Open Challenges (ICLR Blogposts 2026) — silent failures, out-of-simulation, overconfident posteriors. https://iclr-blogposts.github.io/2026/blog/2026/model-misspecification-in-sbi/ [HIGH]
- Talts, Betancourt, Simpson, Vehtari, Gelman (2018). *Validating Bayesian Inference Algorithms with Simulation-Based Calibration* — rank histograms, uniformity. https://sites.stat.columbia.edu/gelman/research/unpublished/sbc.pdf [HIGH]
- Säilynoja, Bürkner, Vehtari (2022). *Graphical Test for Discrete Uniformity…* and *SBC Checking: The Choice of Test Quantities Shapes Sensitivity* — ECDF preferred over histogram; SBC insensitive to some mismatches; test-quantity choice. https://pmc.ncbi.nlm.nih.gov/articles/PMC12490788/ [HIGH]
- Delaunoy, Hermans, Rozet, Wehenkel, Louppe (2022). *Towards Reliable SBI with Balanced Neural Ratio Estimation (BNRE)* — NRE overconfidence; conservative balanced loss. https://arxiv.org/abs/2208.13624 [HIGH]
- *Calibrating Bayesian Tension Statistics using Neural Ratio Estimation* — exponential loss → direct log-BF, avoids sigmoid saturation when BF spans orders of magnitude. https://arxiv.org/pdf/2407.15478 [MEDIUM]
- Anau Montel et al., *Coverage is not enough…* — coverage/SBC under prior ≠ guarantee under real (misspecified) data. https://arxiv.org/html/2605.00980 [MEDIUM]
- NeuralEstimators.jl docs (v0.2.0) — PosteriorEstimator/RatioEstimator/NormalisingFlow/train/sampleposterior/assess; Reactant.jl device; Flux/Lux/SimpleChains backends. https://msainsburydale.github.io/NeuralEstimators.jl/dev/ and /methodology/ [MEDIUM — pre-1.0, API churning]
- State of the SciML Ecosystem 2025 — Flux→Lux migration signal across the Julia ML ecosystem. https://sciml.ai/news/2025/06/26/state_of_sciml/ [MEDIUM]
- Project inputs: `.planning/PROJECT.md`, `.planning/10_plan_amortizedcoloc.md` [HIGH — authoritative for scope/decisions]

---
*Pitfalls research for: amortized SBI (NPE/NRE) on a fluorescence-colocalization physics simulator with SBC + OOD flag (ProteinCoLoc v2.0 spike)*
*Researched: 2026-06-26*
