<!-- GSD:project-start source:PROJECT.md -->
## Project

**ProteinCoLoc v2.0 — AmortizedColoc**

ProteinCoLoc is a Julia package for Bayesian colocalization analysis of multi-channel
fluorescence microscopy images (basis for the author's *Scientific Reports* publication).
**v2.0 ("AmortizedColoc")** extends it with **amortized simulation-based inference (SBI)**:
a 2-channel 2D physics forward-simulator trains a neural posterior estimator (NPE) and a
neural ratio estimator (NRE) on the package's *existing* patch-correlation summary statistic,
delivering calibrated colocalization posteriors and Bayes factors in **milliseconds per
dataset** instead of per-dataset ADVI minutes — with a Simulation-Based-Calibration (SBC)
coverage proof and an out-of-distribution (OOD) / misspecification flag.

The work begins as a **strictly decoupled spike** (`spike/`, own `Project.toml`, no edits to
`src/`) that proves the core principle and ends in a Go/No-Go memo. If the spike succeeds,
v2.0 **productionizes** the amortized inference into the main `src/` package.

**Core Value:** A trained NPE/NRE produces **calibrated, amortized colocalization inference** (posterior +
Bayes factor) in a single forward pass, **>100× faster than the existing per-dataset ADVI**,
with a **demonstrated SBC/coverage calibration proof** and an honest OOD flag — the first
colocalization tool to combine amortized inference, calibration proof, and misspecification
detection.

### Constraints

- **Decoupling**: Spike lives entirely in `ProteinCoLoc/spike/` with its own `Project.toml`/`Manifest.toml`; the main package, `src/`, and both manuscript pipelines must remain **provably untouched** during the spike.
- **Compute**: CPU-only is the portable, reproducible baseline (2D is small enough). A **GPU is available** — CUDA.jl is an *optional accelerator* for training, not a requirement; the spike must run without it.
- **Reproducibility**: Everything reproducible from `spike/demo.jl` with a fixed (Random123) seed.
- **Platform**: Windows-tauglich — Julia-native stack chosen partly to avoid Windows CUDA/Flux friction; GPU use must degrade gracefully to CPU.
- **Summary dimension**: Fixed (8×8 patch grid) during the spike; user-definable only after productionization.
- **Priors**: Simulator prior π(θ) must stay consistent with the existing Turing `@model` prior ranges (μ/ν/σ/τ) for ADVI comparability.
<!-- GSD:project-end -->

<!-- GSD:stack-start source:research/STACK.md -->
## Technology Stack

## Executive Verdict (read first)
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
# In ProteinCoLoc/spike/ — own environment, do NOT touch the root project
# Optional GPU accelerator (only if training is slow; spike must run without it)
# Pkg.add("CUDA")           # loaded as a Flux extension, separate from CPU path
## Key API Notes (AP1 smoke target)
- **NPE path:** build a `DeepSet` summary network (permutation-invariant over patches) -> wrap a built-in `NormalisingFlow` as the approximate distribution -> `PosteriorEstimator(network, q)` -> `train(estimator, sampler, simulator; ...)` -> `sampleposterior(estimator, Z, N)`. Infer **both** ρ_true (for SBC) and Δρ (for the BF comparison).
- **NRE / amortized Bayes factor:** `RatioEstimator` approximates the likelihood-to-evidence ratio; the log-ratio at the observed summary **is** the amortized log Bayes factor (coloc vs null) in one forward pass — no `quadgk`/KDE/shuffle. Validate against `compute_BayesFactor()` in the well-specified regime (AP6).
- **Calibration:** `assess()` gives bias/RMSE; SBC ranks are computed manually (draw θ*~π, simulate, draw L posterior samples, rank θ* among them per parameter) then tested with `HypothesisTests` KS/χ². Reuse `_bin_calibration`/`CalibrationResult` (ECE/MCE + traffic-light) from `BayesInteractomics`.
- **Naming gotcha:** NeuralEstimators' internal flow type is `NormalisingFlow` (British spelling) and is **distinct** from the standalone `NormalizingFlows.jl` package — do not conflate them or add the latter (see "What NOT to Use").
## RxInfer.jl Feasibility Assessment (explicit, per requirements)
- RxInfer's speed advantage comes from exploiting *local conjugate* relationships on the factor graph (BP/VMP). A **Student-t likelihood is non-conjugate to a Gaussian prior**, so it is not natively fast.
- Two escape hatches exist, both with costs:
- Net: it is achievable, but it is a genuine model-rebuild + diagnostics-revalidation effort, not a drop-in backend swap. RxInfer programs are factor-graph `@model`s with different semantics from Turing's, so the existing `@model` does not port verbatim.
- The entire thesis of v2.0 is to *replace* per-dataset ADVI-minutes with amortized millisecond inference. A faster per-dataset VMP/CVI baseline optimizes the very axis amortized SBI makes obsolete. The ADVI you already have is sufficient as the accuracy/RMSE ground-truth baseline for AP4.
- Where RxInfer *could* add value is purely as an **independent cross-check** of the ADVI posterior (a second non-amortized reference), strengthening the calibration story. That is a "nice-to-have for the paper," not a spike gate.
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
- Flux backend, CPU-only, built-in `NormalisingFlow`, JLD2 cache. Ship the whole spike on this. Simplest reproducible path.
- Add `CUDA.jl` as a Flux extension (no code rewrite, graceful CPU fallback retained), or switch backend to **Lux.jl + Reactant** via NeuralEstimators' backend support. Try CUDA-on-Flux first (smaller change).
- Escalate to **BayesFlow via PythonCall/CondaPkg** as the documented fallback. Budget for Conda env setup + PythonCall marshalling + dual-language seed/reproducibility handling.
- Re-evaluate **InvertibleNetworks.jl** (memory-scalable GPU flows) and Lux/Reactant. Not before.
## Version Compatibility
| Package | Compatible With | Notes |
|---------|-----------------|-------|
| NeuralEstimators v0.2.1 | Flux v0.16.x, Lux v1.x, SimpleChains, Reactant | Backend-agnostic; pick one backend per environment. Verify exact compat bounds in the spike Manifest after resolve. |
| Flux v0.14+ | CUDA.jl (as extension) | Since v0.14 CUDA loads as a package extension — CPU installs stay lean, GPU is opt-in and graceful. |
| RxInfer v1.17.0 | Julia 1.12.x (docs built on 1.12.6) | Non-conjugate path needs ExponentialFamilyProjection; expect extra setup. |
| JLD2 / Random123 | pure Julia | No native binary deps — Windows-clean. |
| HypothesisTests | Distributions, StatsBase | KS/χ² for SBC ranks; standard, stable. |
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
<!-- GSD:stack-end -->

<!-- GSD:conventions-start source:CONVENTIONS.md -->
## Conventions

Conventions not yet established. Will populate as patterns emerge during development.
<!-- GSD:conventions-end -->

<!-- GSD:architecture-start source:ARCHITECTURE.md -->
## Architecture

Architecture not yet mapped. Follow existing patterns found in the codebase.
<!-- GSD:architecture-end -->

<!-- GSD:skills-start source:skills/ -->
## Project Skills

No project skills found. Add skills to any of: `.claude/skills/`, `.agents/skills/`, `.cursor/skills/`, `.github/skills/`, or `.codex/skills/` with a `SKILL.md` index file.
<!-- GSD:skills-end -->

<!-- GSD:workflow-start source:GSD defaults -->
## GSD Workflow Enforcement

Before using Edit, Write, or other file-changing tools, start work through a GSD command so planning artifacts and execution context stay in sync.

Use these entry points:
- `/gsd-quick` for small fixes, doc updates, and ad-hoc tasks
- `/gsd-debug` for investigation and bug fixing
- `/gsd-execute-phase` for planned phase work

Do not make direct repo edits outside a GSD workflow unless the user explicitly asks to bypass it.
<!-- GSD:workflow-end -->



<!-- GSD:profile-start -->
## Developer Profile

> Profile not yet configured. Run `/gsd-profile-user` to generate your developer profile.
> This section is managed by `generate-claude-profile` -- do not edit manually.
<!-- GSD:profile-end -->
