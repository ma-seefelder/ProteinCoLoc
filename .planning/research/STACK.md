# Stack Research

**Domain:** Amortized simulation-based inference (NPE + NRE) on a patch-correlation summary statistic, added to a Julia fluorescence-colocalization package (Windows, CPU-first, optional CUDA GPU)
**Researched:** 2026-06-26
**Confidence:** HIGH for the core SBI stack (NeuralEstimators/Flux/Lux, current versions verified against official docs + GitHub releases); MEDIUM for RxInfer Student-t feasibility (verified mechanism, but no published colocalization-shaped example).

## Executive Verdict (read first)

1. **NeuralEstimators.jl is the correct primary engine.** v0.2.1 (April 2026) ships `PosteriorEstimator` (NPE via a built-in `NormalisingFlow` approximate distribution), `RatioEstimator` (NRE / amortized likelihood-to-evidence ratio = your amortized Bayes factor), plus `train`, `sampleposterior`, `estimate`, and `assess`. One package covers NPE **and** NRE on your existing summary vector. Julia-native, Windows-friendly, CPU-sufficient for the 2D spike. This single decision de-risks the whole spike.
2. **Flux.jl is the default NN backend; Lux.jl is the upgrade path, not a competing choice.** NeuralEstimators is backend-agnostic (Flux, Lux, SimpleChains, and Reactant/XLA). Start on Flux (most examples, smallest mental overhead, fine on CPU). Keep Lux+Reactant in reserve for GPU/throughput if AP4 training becomes the bottleneck.
3. **Do NOT pull in TuringLang/NormalizingFlows.jl or InvertibleNetworks.jl for the spike.** NeuralEstimators already owns the flow. The other two solve problems you do not have yet (custom VI flows / GPU-scalable image-sized flows for the 3D full build-out).
4. **RxInfer.jl as a Turing replacement: feasible only after a model reformulation, and largely *redundant* with the amortized SBI goal.** A hierarchical Student-t is non-conjugate as written; RxInfer needs either CVI/ExponentialFamilyProjection (slower, stochastic, accuracy cost) or a Gaussian scale-mixture rewrite (conjugate, fast, but real modelling work). Given that the entire point of v2.0 is to replace per-dataset ADVI with millisecond amortized inference, investing in a faster *per-dataset* MCMC/VMP baseline is low-value. Treat RxInfer as an **optional, deferred baseline experiment**, not a spike dependency. Verdict: **MAYBE-feasible, recommend DEFER.**

## Recommended Stack

### Core Technologies

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| **NeuralEstimators.jl** | **v0.2.1** (Apr 2026) | NPE (`PosteriorEstimator`) + NRE (`RatioEstimator`) + `train`/`sampleposterior`/`estimate`/`assess` | The only Julia-native package that delivers both amortized posterior **and** amortized ratio/Bayes-factor estimation from one API, designed exactly for "simulate → summary → train → infer on many datasets." Built-in `NormalisingFlow` distribution removes the need for a separate flow library. Backend-agnostic. Has a companion CRAN R package, signalling maintained, documented, multi-language maturity. |
| **Flux.jl** | **v0.16.x** | Neural-network backend for the summary net (MLP/DeepSet) and flow conditioner | Pure-Julia, Windows-tauglich, the most-documented backend for NeuralEstimators examples. Since v0.14, CUDA is an optional, separately-loaded extension — so CPU-only installs are clean and GPU degrades gracefully (exactly your constraint). |
| **Distributions.jl** | current | Prior π(θ) sampling, must match the Turing `@model` ranges (μ/ν/σ/τ) | Already the lingua franca; lets the simulator prior stay provably consistent with the ADVI baseline for comparability. |
| **JLD2.jl** | current | Cache the 50k–200k simulated (θ, summary) training pairs and trained models | Pure-Julia, no HDF5 binary dependency (cleaner on Windows), the de-facto Julia serialization choice and already the pattern in the sibling `BayesInteractomics` simulation engine. |
| **Random123.jl** | current | Counter-based, reproducible seeding of the simulator and data generator | Stateless, splittable, reproducible across processes — the right tool for "everything reproducible from `demo.jl` with a fixed seed," more robust than `Random.seed!` for parallel simulation. |
| **HypothesisTests.jl** | current | KS / one-sample uniformity and χ² goodness-of-fit on SBC rank statistics | Standard Julia stats package; `ExactOneSampleKSTest` / `ChisqTest` give the p-values for the SBC rank-uniformity claim (the calibration proof). |

### Supporting Libraries

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| **Images.jl** | current | Image container/ops to assemble `MultiChannelImage`-compatible arrays | AP2 simulator output, so existing `correlation()`/`patch()` apply unchanged. |
| **ImageFiltering.jl** | current | Gaussian PSF convolution (`imfilter`, `Kernel.gaussian`), sub-pixel shift | AP2: PSF blur + registration shift in the forward model. |
| **StatsBase.jl** | current | Standardization of summary vectors, histograms, rank computation | AP3 (z-scoring summaries) and AP5 (rank histograms, coverage curves). |
| **KernelDensity.jl** | current | KDE for the existing `compute_BayesFactor()` baseline | AP6 only, to reproduce the baseline BF you validate the NRE against (already a dependency). |
| **QuadGK.jl** | current | Numerical integration in the existing KDE Bayes-factor baseline | AP6 baseline reproduction (already a dependency). |
| **CUDA.jl** | current (optional) | Optional GPU acceleration of `train` | Only if AP4 training time hurts. Loaded as a Flux extension; spike must run without it. |
| **Lux.jl** | current (v1.x) | Alternative NN backend for GPU/Reactant throughput | Reserve. Switch only if Flux GPU path or compile/perf becomes a blocker; NeuralEstimators supports it natively. |
| **MLUtils.jl** | current | Batching/`DataLoader` for training pairs | If you need custom minibatching beyond what `train` provides. |
| **CairoMakie.jl** *(or Plots.jl)* | current | SBC rank histograms, coverage curves, reliability/traffic-light plots | AP5/AP6 diagnostics and the demo figures. |

### Development Tools

| Tool | Purpose | Notes |
|------|---------|-------|
| Isolated `spike/Project.toml` + `Manifest.toml` | Provable decoupling from main package | `Pkg.activate("spike")`; never touch the root `Project.toml`/`Manifest.toml`. This is a hard constraint, not a nicety. |
| `BenchmarkTools.jl` | Wall-clock NPE-ms vs ADVI-minutes claim (AP4) | Needed to substantiate the ">100× speedup" Definition-of-Done. |
| Context7 / `npx ctx7` | Live API lookup for NeuralEstimators during AP1 smoke | API surface is young (v0.x); verify signatures against `dev` docs, not memory. |

## Installation

```julia
# In ProteinCoLoc/spike/ — own environment, do NOT touch the root project
using Pkg
Pkg.activate(".")            # spike/

Pkg.add([
    "NeuralEstimators",      # NPE + NRE engine (v0.2.1)
    "Flux",                  # NN backend (v0.16.x)
    "Distributions",         # priors matching the Turing @model
    "JLD2",                  # cache training pairs + models
    "Random123",             # seeded, reproducible simulation
    "Images", "ImageFiltering",  # PSF convolution, shift, noise
    "StatsBase",             # standardization, ranks, histograms
    "KernelDensity", "QuadGK",   # baseline KDE Bayes factor (AP6)
    "HypothesisTests",       # SBC uniformity (KS / chi-square)
    "BenchmarkTools",        # speedup measurement
    "CairoMakie",            # diagnostic plots
])

# Optional GPU accelerator (only if training is slow; spike must run without it)
# Pkg.add("CUDA")           # loaded as a Flux extension, separate from CPU path
```

## Key API Notes (AP1 smoke target)

- **NPE path:** build a `DeepSet` summary network (permutation-invariant over patches) -> wrap a built-in `NormalisingFlow` as the approximate distribution -> `PosteriorEstimator(network, q)` -> `train(estimator, sampler, simulator; ...)` -> `sampleposterior(estimator, Z, N)`. Infer **both** ρ_true (for SBC) and Δρ (for the BF comparison).
- **NRE / amortized Bayes factor:** `RatioEstimator` approximates the likelihood-to-evidence ratio; the log-ratio at the observed summary **is** the amortized log Bayes factor (coloc vs null) in one forward pass — no `quadgk`/KDE/shuffle. Validate against `compute_BayesFactor()` in the well-specified regime (AP6).
- **Calibration:** `assess()` gives bias/RMSE; SBC ranks are computed manually (draw θ*~π, simulate, draw L posterior samples, rank θ* among them per parameter) then tested with `HypothesisTests` KS/χ². Reuse `_bin_calibration`/`CalibrationResult` (ECE/MCE + traffic-light) from `BayesInteractomics`.
- **Naming gotcha:** NeuralEstimators' internal flow type is `NormalisingFlow` (British spelling) and is **distinct** from the standalone `NormalizingFlows.jl` package — do not conflate them or add the latter (see "What NOT to Use").

## RxInfer.jl Feasibility Assessment (explicit, per requirements)

**Question:** Can the hierarchical Student-t colocalization model migrate Turing -> RxInfer.jl (v1.17.0, June 2026) for a faster Bayesian baseline, and is it complementary to or redundant with amortized SBI?

**Verdict: MAYBE-feasible with reformulation; recommend DEFER — largely redundant with the v2.0 amortized goal.** Confidence: MEDIUM.

**Why feasible-but-hard:**
- RxInfer's speed advantage comes from exploiting *local conjugate* relationships on the factor graph (BP/VMP). A **Student-t likelihood is non-conjugate to a Gaussian prior**, so it is not natively fast.
- Two escape hatches exist, both with costs:
  1. **CVI / `ExponentialFamilyProjection`** — re-projects non-conjugate factors back into the exponential family. Official docs explicitly warn this uses stochastic gradient approximations, increases runtime "significantly," and costs accuracy. You also must hand-write `@constraints`, `@initialization`, and projection targets.
  2. **Gaussian scale-mixture reparameterization** (the principled route): a Student-t is exactly a Gaussian whose precision is Gamma-distributed (`t_ν = ∫ N(x|μ, σ²/τ)·Gamma(τ|ν/2, ν/2) dτ`). With ν **fixed**, the Gamma scale node is conjugate and VMP runs fast and cleanly — this is the standard message-passing robust-regression trick. With ν **also inferred** (your model has ν as a parameter), the ν update is non-conjugate and pushes you back to CVI/projection or fixing/gridding ν.
- Net: it is achievable, but it is a genuine model-rebuild + diagnostics-revalidation effort, not a drop-in backend swap. RxInfer programs are factor-graph `@model`s with different semantics from Turing's, so the existing `@model` does not port verbatim.

**Why redundant (the decisive point):**
- The entire thesis of v2.0 is to *replace* per-dataset ADVI-minutes with amortized millisecond inference. A faster per-dataset VMP/CVI baseline optimizes the very axis amortized SBI makes obsolete. The ADVI you already have is sufficient as the accuracy/RMSE ground-truth baseline for AP4.
- Where RxInfer *could* add value is purely as an **independent cross-check** of the ADVI posterior (a second non-amortized reference), strengthening the calibration story. That is a "nice-to-have for the paper," not a spike gate.

**Recommendation:** Do not add RxInfer to the spike environment or critical path. If raised at full build-out, scope it as a separate, time-boxed experiment using the **fixed-ν Gaussian-scale-mixture** formulation first (the only version likely to actually be faster than Turing). Document it as complementary-baseline, not replacement.

## Alternatives Considered

| Recommended | Alternative | When to Use Alternative |
|-------------|-------------|-------------------------|
| NeuralEstimators built-in `NormalisingFlow` | **TuringLang/NormalizingFlows.jl** | Only if you need a fully custom flow / variational objective Bijectors-integration that NeuralEstimators cannot express. Not needed for the spike; adds a second flow stack. |
| NeuralEstimators flow | **InvertibleNetworks.jl** | Memory-scalable, GPU-first flows for *image-sized* densities (it beats PyTorch at >1024×1024 on an A100). Reserve for the 3D / full-hierarchy build-out, not the 2D summary-vector spike. |
| Flux.jl backend | **Lux.jl (+ Reactant/XLA)** | If GPU throughput or AOT-compiled performance becomes the AP4 bottleneck; NeuralEstimators supports it natively, so switching is a backend flag, not a rewrite. |
| NeuralEstimators (Julia-native) | **BayesFlow (Python via PythonCall/CondaPkg)** | Documented **fallback only**: if the Julia flow fails to converge or calibrate after tuning. Integration cost is real — CondaPkg env, PythonCall data marshalling, two-language reproducibility/seed management, Windows Conda friction. A CondaPkg pattern exists in `BayesInteractomics`, but keep it off the core path. |
| Turing ADVI baseline (keep) | **RxInfer.jl VMP/CVI** | Only as a deferred independent baseline cross-check, and only via fixed-ν Gaussian scale-mixture (see feasibility section). Not for the spike. |
| JLD2 | BSON / HDF5.jl | If you need cross-language/HDF5 interop; otherwise JLD2 is cleaner on Windows. |

## What NOT to Use

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| **NormalizingFlows.jl (TuringLang)** *for the spike* | Separate flow stack; duplicates NeuralEstimators' built-in `NormalisingFlow`; noted GPU-compat caveats (needs Lux/Reactant work). Adds dependency surface for zero spike benefit. | NeuralEstimators' built-in `NormalisingFlow`. |
| **InvertibleNetworks.jl** *for the spike* | Solves high-dimensional image-sized flow scalability you don't have on an 8×8 summary vector; pure overhead now. | NeuralEstimators flow; revisit for 3D/hierarchy build-out only. |
| **RxInfer.jl as a Turing replacement** *in the spike* | Non-conjugate Student-t needs reformulation; optimizes the per-dataset axis amortized SBI is designed to retire. | Keep Turing ADVI as the accuracy baseline; defer RxInfer. |
| **BayesFlow on the core path** | Two-language complexity, Conda/Windows friction, reproducibility burden. | Julia-native NeuralEstimators; BayesFlow strictly as last-resort fallback. |
| **GPU as a requirement** | Windows Flux/CUDA friction is the stated risk; 2D summary-vector training is small. | CPU-only baseline; CUDA.jl as optional Flux extension that degrades gracefully. |
| Global env install of spike deps | Violates the hard decoupling constraint (main `Project.toml`/`Manifest.toml` must stay untouched). | `Pkg.activate("spike")` with its own manifest. |

## Stack Patterns by Variant

**If AP1 smoke converges and CPU training is fast enough (expected for 2D):**
- Flux backend, CPU-only, built-in `NormalisingFlow`, JLD2 cache. Ship the whole spike on this. Simplest reproducible path.

**If AP4 training is too slow on CPU:**
- Add `CUDA.jl` as a Flux extension (no code rewrite, graceful CPU fallback retained), or switch backend to **Lux.jl + Reactant** via NeuralEstimators' backend support. Try CUDA-on-Flux first (smaller change).

**If the Julia flow won't calibrate after tuning (AP5 SBC fails persistently, not just reports honest miscalibration):**
- Escalate to **BayesFlow via PythonCall/CondaPkg** as the documented fallback. Budget for Conda env setup + PythonCall marshalling + dual-language seed/reproducibility handling.

**If, at full build-out (post-Go), you move to 3D / hierarchical / image-sized densities:**
- Re-evaluate **InvertibleNetworks.jl** (memory-scalable GPU flows) and Lux/Reactant. Not before.

## Version Compatibility

| Package | Compatible With | Notes |
|---------|-----------------|-------|
| NeuralEstimators v0.2.1 | Flux v0.16.x, Lux v1.x, SimpleChains, Reactant | Backend-agnostic; pick one backend per environment. Verify exact compat bounds in the spike Manifest after resolve. |
| Flux v0.14+ | CUDA.jl (as extension) | Since v0.14 CUDA loads as a package extension — CPU installs stay lean, GPU is opt-in and graceful. |
| RxInfer v1.17.0 | Julia 1.12.x (docs built on 1.12.6) | Non-conjugate path needs ExponentialFamilyProjection; expect extra setup. |
| JLD2 / Random123 | pure Julia | No native binary deps — Windows-clean. |
| HypothesisTests | Distributions, StatsBase | KS/χ² for SBC ranks; standard, stable. |

> Pin exact versions in `spike/Manifest.toml` after the first `Pkg.resolve()` and commit it — reproducibility is a Definition-of-Done item.

## Sources

- https://msainsburydale.github.io/NeuralEstimators.jl/dev/ — official docs: PosteriorEstimator/RatioEstimator/train/sampleposterior/assess, Flux+Lux+SimpleChains+Reactant backend support, GPU/TPU via Reactant. HIGH.
- https://github.com/msainsburydale/NeuralEstimators.jl — README + releases: v0.2.1 (released 2026-04-07), supports both Flux and Lux; NPE/NRE/NBE methods; active (605 commits, 7 releases). HIGH.
- https://msainsburydale.github.io/NeuralEstimators.jl/dev/methodology/ + /API/core/ — NPE via KL-divergence min, NRE via likelihood-to-evidence ratio, `sampleposterior` signatures. HIGH.
- https://cran.r-project.org/package=NeuralEstimators — companion CRAN R package, indicates maintained/documented project. MEDIUM (maturity signal).
- http://fluxml.ai/Flux.jl/dev/guide/gpu/ + fluxml.ai — Flux v0.16.5 stable; CUDA separate since v0.14, Windows NVIDIA path via CUDA+cuDNN. HIGH.
- https://lux.csail.mit.edu/dev/manual/gpu_management — Lux GPU/device management, conditional GPU use. HIGH.
- https://github.com/TuringLang/NormalizingFlows.jl + https://arxiv.org/abs/2312.13480 (InvertibleNetworks.jl JOSS) — roles/tradeoffs: NormalizingFlows.jl = modular Bijectors-compatible VI flows; InvertibleNetworks.jl = memory-scalable GPU image flows (>1024² vs PyTorch 480² on A100). HIGH.
- https://docs.rxinfer.com/stable/ — RxInfer v1.17.0 (generated 2026-06-18, Julia 1.12.6); conjugate-exponential-family core, BP/VMP. HIGH.
- https://github.com/ReactiveBayes/RxInfer.jl/blob/main/docs/src/manuals/inference/nonconjugate.md — non-conjugate handled via ExponentialFamilyProjection with explicit runtime/accuracy warnings and required @constraints/@initialization. HIGH.
- Gaussian-scale-mixture representation of Student-t for message passing — standard robust-inference result (Student-t = Gaussian with Gamma-distributed precision); applied judgment, MEDIUM confidence as the recommended RxInfer reformulation.
- https://arxiv.org/pdf/1804.06788 (Talts et al., SBC) + sbi-dev SBC tutorial — SBC rank-histogram + KS-uniformity methodology mapped to HypothesisTests.jl. HIGH (method), MEDIUM (Julia-specific wiring is hand-rolled).

---
*Stack research for: amortized SBI (NPE+NRE) on a colocalization summary statistic, Julia, Windows, CPU-first*
*Researched: 2026-06-26*
