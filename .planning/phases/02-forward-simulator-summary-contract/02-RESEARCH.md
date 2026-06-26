# Phase 2: Forward Simulator + Summary Contract - Research

**Researched:** 2026-06-26
**Domain:** 2D physics forward-simulation for SBI; microscopy image synthesis; prior-consistency calibration; Julia image/stats stack
**Confidence:** HIGH (summary contract verified by source read; stack verified against CLAUDE.md + NeuralEstimators registry; calibration method is a defensible MEDIUM-confidence design)

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

- **D-01:** **Full induced-μ calibration.** Fit/transform the `ρ_true` prior so the simulator's **induced** distribution of the per-patch correlation mean (μ — the quantity the Turing model infers) matches the Turing μ-prior `Truncated(Cauchy(0,0.3),-1,1)` as closely as possible — not merely adopting the range by assertion. One generative prior must serve both the ADVI benchmark (infers μ) and the NPE/BF validation (infers ρ_true / Δρ).
- **D-02:** The calibration is a concrete step: sweep `ρ_true` → measure induced per-patch-correlation μ via the **real** `patch()`/`correlation()` summary on simulated pairs → fit the transform/prior so induced μ matches the Turing μ-prior. The fitted mapping AND the calibration evidence (ρ_true→μ curve, induced-vs-target μ overlap) are documented in `spike/NOTES.md`.
- **D-03:** Noise/scatter nuisances of θ map to Turing scale hyperpriors: per-patch correlation spread ↔ `σ ~ Truncated(Cauchy(0.1,0.3),1e-4,1)` / `τ`; tail/heaviness ↔ `ν ~ Exponential()`. Keep consistent with Turing ranges **where the simulator induces them**; exact mapping form is research/planning detail.
- **D-04:** θ is the **7-parameter** vector: `θ = (ρ_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise)`.
- **D-05:** **PSF is a FIXED Gaussian nuisance** — NOT an inferred θ-component. Fixed realistic width tied to pixel scale (researcher's discretion).
- **D-06:** Pipeline order (ROADMAP SC-1): correlated bivariate densities → **Bernoulli thinning** (label efficiency = keep-probability) → **Gaussian PSF** (`ImageFiltering.imfilter`, `Kernel.gaussian`) → **2×2 directional spillover** (bleed-through, directional ≠ symmetric) → **autofluorescence** offset → **sub-pixel shift** (dx/dy) → **Poisson (shot) + Gaussian (read) noise**.
- **D-07:** Noise is **combined Poisson + Gaussian** (shot + read), not Gaussian-only.
- **D-08:** **Anchor on the real paper size: 1376 × 1028, 16-bit** (`test/test_images/{positive,negative}/*_c{1,2,3}.tif`).
- **D-09:** **Also sweep other sensible microscopy sizes** — image dimensions are a **configurable** parameter; evaluate a range of common confocal/widefield dimensions (exact set is a research item).
- **D-10:** **The 8×8 patch grid stays FIXED** for the spike, feeding the unchanged `patch()`/`correlation()`.
- **D-11:** `MultiChannelImage.pixel_size` carries the **(width, height) pixel-dimension tuple** matching `src/LoadImages.jl:444` (`pixel_size = size(data[1])`) — NOT micrometers despite the docstring.
- **D-12:** **Quantitative plausibility gate.** SIM-04 passes only if: (a) ρ_true sweep shows **monotone** ρ_true↑ → mean patch-correlation↑ (Spearman ≥ threshold; threshold = discretion), AND (b) spillover and sub-pixel shift **measurably** perturb the pair, AND (c) output verifies as a valid `MultiChannelImage` (SIM-03).
- **D-13:** **Plotting backend = CairoMakie** (static, headless, vector). **Spike-local only** (under `spike/`), zero `src/` risk. Does NOT resolve the main-package `src/plot.jl` GLMakie question.
- **D-14:** Simulator accepts an **explicit `rng::AbstractRNG`** argument (threaded through `sample_prior`/`simulate_pair`). Phase 2 uses **stdlib `Random`** (Random123 deferred to Phase 3).

### Claude's Discretion

- Exact fixed PSF Gaussian width (tie to pixel scale).
- Exact Poisson gain + Gaussian read-noise constants (plausible microscopy ranges).
- The Spearman threshold for the SIM-04 monotonicity gate.
- The exact calibration-fit method for D-01/D-02 (regression / quantile-matching of ρ_true→μ).
- The exact set of "sensible microscopy sizes" to sweep (D-09).
- Internal simulator file/module layout within `spike/`.

### Deferred Ideas (OUT OF SCOPE)

- **Random123 counter-based seeding** — Phase 3 (D-14).
- **`src/plot.jl` GLMakie → CairoMakie/Plots migration** for the MAIN package — separate post-spike decision.
- **Training-data generation at 50k–200k scale, NPE/NRE training, ADVI baseline, Bayes-factor reproduction** — Phases 3, 4, 6.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| **SIM-01** | `simulate_pair(θ) → MultiChannelImage` generates a 2-channel pair from the 7-param θ via the D-06 pipeline. | §Forward-Simulation Pipeline gives a verified Julia stage-by-stage recipe (correlated fields → thinning → PSF → spillover → autofluorescence → sub-pixel shift → Poisson+Gaussian noise) with exact `Images`/`ImageFiltering`/`Distributions` calls. |
| **SIM-02** | π(θ) / `sample_prior()` consistent with Turing μ/ν/σ/τ ranges, documented in `spike/NOTES.md`. | §Induced-μ Calibration gives the concrete sweep→fit→evidence method (monotone-inverse / quantile matching) and the σ/τ/ν consistency checks. |
| **SIM-03** | Output verifies as a valid `MultiChannelImage` so existing `correlation()`/`patch()`/`_prepare_data()` apply unchanged. | §Summary Contract documents the exact verified data flow, array shapes, element types, and constructor requirements — the simulator must emit `Vector{Matrix{Float64}}` of 2 channels. |
| **SIM-04** | Plausibility plots: ρ_true↑→patch-corr↑; spillover & shift visibly affect the pair. | §Plausibility Gate gives the quantitative bar (Spearman monotonicity + paired perturbation effect) and the CairoMakie figure recipe; §Validation Architecture gives the re-runnable harness. |
</phase_requirements>

## Summary

Phase 2 builds a pure-Julia, CPU-only 2D forward simulator that emits real `MultiChannelImage` pairs which the **unmodified** existing summary stack (`patch()` → `correlation()` → `_prepare_data()`) ingests. The summary contract is now fully verified by source read (see §Summary Contract): the simulator must produce a 2-element `Vector{Matrix{Float64}}` whose two channel matrices, after `patch.(_, 8)` and `correlation(...)`, yield an 8×8 grid of per-patch Pearson correlations; the **mean of those 64 values is the induced μ** that the Turing `@model` infers as `μ ~ Truncated(Cauchy(0,0.3),-1,1)`. That identity is what makes SIM-02 a *mapping* problem (physics ρ_true → induced summary-stat μ), not a copy of ranges.

The entire stage stack is covered by mature, Windows-clean Julia packages already named in CLAUDE.md: `Distributions` (already in the spike env) for correlated bivariate densities, Bernoulli thinning, and Poisson+Gaussian noise; `ImageFiltering` (`imfilter`, `Kernel.gaussian`) for the fixed Gaussian PSF; `ImageTransformations`/`Interpolations` (via the `Images` umbrella) for sub-pixel shift; `StatsBase`/`Statistics` for the correlation/standardization the `src/` code calls; `CairoMakie` for the static plausibility figures; `HypothesisTests` for the SBC-style uniformity tests reused later. **Critical de-risking finding:** the flagged CairoMakie resolve risk is *low* — NeuralEstimators 0.2.1 declares `Makie = "0.24"` as a weak-dependency compat bound `[VERIFIED: github.com/msainsburydale/NeuralEstimators.jl Project.toml]`, so a *current* CairoMakie (Makie 0.24 series) **aligns** with the pinned estimator rather than conflicting. The Phase-1 downgrade was caused by the parent's *old* GLMakie 0.10 / Makie 0.21 tree, which the lean spike env never pulls in.

**Primary recommendation:** Implement `simulate_pair(rng, θ; imsize)` and `sample_prior(rng)` in `spike/` as pure functions over `Vector{Matrix{Float64}}`, reaching `src/` summary functions read-only via the proven Phase-1 `include()` coupling; build the induced-μ calibration as a deterministic offline sweep that fits a monotone ρ_true→μ map and is frozen into `sample_prior`; gate the phase on a quantitative Spearman-monotonicity + paired-perturbation test in an extended `spike/test/` harness, with CairoMakie evidence figures. Add `Images`, `ImageFiltering`, `StatsBase`, `Statistics`, `CairoMakie`, `HypothesisTests` to the spike env, re-resolve, **re-run the Phase-1 smoke to confirm NeuralEstimators stays at 0.2.1 with no CUDA**, and re-freeze the Manifest.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| θ prior sampling `sample_prior(rng)` | Simulator (spike) | — | New code; must be statistically consistent with the `src/bayes.jl` Turing prior but lives entirely in `spike/`. |
| Physics image synthesis `simulate_pair` | Simulator (spike) | — | New forward model; pure function producing `Vector{Matrix{Float64}}`. |
| Image container `MultiChannelImage` | `src/` (read-only) | Simulator constructs instances | Constructor is the contract boundary; reached via `include()`. |
| Summary statistic `patch`/`correlation`/`_prepare_data` | `src/` (read-only, UNCHANGED) | Simulator calls it | SIM-03 mandate: no `src/` edits; the simulator is a *client* of the summary. |
| Induced-μ calibration | Simulator (spike) | uses `src/` summary as the measurement instrument | Calibration *measures* via the real summary, then fits a transform stored in the spike prior. |
| Plausibility plots / gate | Spike test + CairoMakie | — | Spike-local validation; never touches `src/plot.jl`. |
| RNG / reproducibility | Simulator (spike, stdlib `Random`) | — | Explicit `rng` threaded (D-14); Random123 deferred to Phase 3. |

## Summary Contract (SIM-03) — VERIFIED data flow

Read directly from `src/colocalization.jl` and `src/bayes.jl` this session. This is the load-bearing contract; the simulator output MUST satisfy it with **zero `src/` edits**.

### The exact call sequence

```
MultiChannelImage.data :: Vector{Matrix{Float64}}   # 2 channels, e.g. data[1], data[2]
        │
        │  channels = [1, 2]   (integer indices into data)
        ▼
x = image.data[1];  y = image.data[2]                # each :: Matrix{Float64}, size (W, H)
        │
        ▼  x, y = patch.([x, y], 8)                   # broadcast patch(img, num_patches=8)
patch(::Matrix, 8) ⇒ Array{Union{Float64,Missing},4} # shape: (8, 8, W÷8, H÷8)
        │                                             # dims 1–2 = patch grid; dims 3–4 = pixels-in-patch
        ▼  correlation(x, y; method=:pearson)
correlation(::Array{T,4}, ::Array{T,4})              # T = Union{Float64,Missing}
        │   per patch (i,j): a = vec(x[i,j,:,:]); b = vec(y[i,j,:,:])
        │   _exclude_zero(a,b)  → drops 0.0 and NaN pairs
        │   length(a) ≤ 15 ? missing : cor(a, b)
        ▼
ρ :: Matrix{Union{Float64,Missing}}  size (8, 8)     # 64 per-patch Pearson correlations
        │
        ▼  (_prepare_data path, on a MultiChannelImageStack of 1 image)
reshape → (1, 64), filter !ismissing & !isnan
        ▼
sample_data :: Vector{Vector{Float64}}               # one vector of ≤64 per-patch ρ per image

INDUCED μ  ==  mean(skipmissing(ρ))                   # the Turing model's μ for this image
```

### Hard requirements the simulator must meet

| Requirement | Source evidence | Implication for simulator |
|-------------|-----------------|---------------------------|
| `data :: Vector{Matrix{Float64}}`, length 2 | `MultiChannelImage` struct `data::Vector{Matrix{T}}`, T<:Union{Missing,Float64}; constructor asserts `length(data)==length(channels)==length(otsu_threshold)` (`LoadImages.jl:118-130`) | Emit exactly 2 `Matrix{Float64}` channels + 2 channel names + 2 otsu thresholds. `Float64 <: Union{Float64,Missing}` so plain `Matrix{Float64}` is accepted. |
| `pixel_size = size(data[1])` as `(W,H)::Tuple{Int,Int}` | `LoadImages.jl:444` convenience constructor; D-11 | Set `pixel_size = size(data[1])`, NOT micrometers. |
| `otsu_threshold :: Vector{<:AbstractFloat}`, length 2 | constructor assert | Compute via `Images.otsu_threshold.(data)` OR supply plausible floats. **Not** used by `_prepare_data`/`colocalization` (only by `_calculate_mask`, which the Turing path does not call), so any valid float passes the contract. |
| Each 8×8 patch must have **> 15** non-zero, non-NaN pixel pairs | `correlation` sets `missing` when `length(a) ≤ 15` after `_exclude_zero` | Trivially satisfied at all swept sizes (smallest is 256² → patch 32×32 = 1024 px). Simulator must emit mostly **non-zero** intensities (zeros are dropped by `_exclude_zero`); positive-valued intensity fields satisfy this. |
| No masking is applied before the summary | `colocalization()` calls `_prepare_data` directly on raw `image.data`; `apply_mask!` is NOT invoked on the Turing path | The simulator need not pre-mask; positive intensities flow straight through. |
| `patch(img, 8)` requires `size(img,1)÷8` and `size(img,2)÷8 ≥ 1` and trims remainder | `patch` lines 37-60 (`rows ÷ num_patches`, `@view img[1:8*patch_size_x, ...]`) | Any size ≥ 8 in each dim works; non-multiples of 8 are silently trimmed (e.g. 1028÷8=128 → uses 1024 rows). |

### How to exercise the full contract in the test (proves SIM-03 end-to-end)

```julia
mci   = MultiChannelImage(data, ["ch1","ch2"], "sim", ["",""], size(data[1]),
                          Images.otsu_threshold.(data))
stack = MultiChannelImageStack([mci], "sim_stack")
ρvec  = _prepare_data(stack, [1, 2], 8)          # exercises the real summary end-to-end
μ_ind = Statistics.mean(reduce(vcat, ρvec))      # induced μ for this pair
```
Also call `patch`/`correlation` directly to assert the 8×8 shape independently.

## Standard Stack

### Core (add to `spike/Project.toml`)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `Distributions` | 0.25.x (already pinned 0.25.128) | Correlated bivariate densities (`MvNormal`), Bernoulli thinning (`Bernoulli`), Poisson shot noise (`Poisson`), Gaussian read noise (`Normal`); truncated Cauchy/Exponential priors | Already in the spike env; the JuliaStats lingua franca; keeps π(θ) provably consistent with the Turing `@model` `[CITED: CLAUDE.md]`. |
| `ImageFiltering` | current (Images-ecosystem) | Gaussian PSF via `imfilter(img, Kernel.gaussian(σ))` | Named in CLAUDE.md; the standard Julia Gaussian-convolution path; pure-Julia, Windows-clean `[CITED: CLAUDE.md]`. |
| `Images` (umbrella) | current | Re-exports `imfilter`, `Kernel`, `warp`, `imresize`, `Translation`, `otsu_threshold`; one dependency instead of many | Single mature umbrella; matches the `src/` code's own `Images.otsu_threshold` usage `[ASSUMED re-export set — verify at install]`. |
| `StatsBase` | current | `corspearman`/`corkendall` referenced by `correlation`; summary standardization (Phase 3); Spearman gate | **Required by `src/colocalization.jl`** before `include()` (per `spike/NOTES.md` §4) `[VERIFIED: spike/NOTES.md]`. |
| `Statistics` (stdlib) | bundled | `mean`/`median`/`quantile` referenced by `src/colocalization.jl`; induced-μ computation | stdlib, no version `[VERIFIED: spike/NOTES.md]`. |
| `CairoMakie` | Makie 0.24-series (e.g. 0.15.x) | Static/vector plausibility figures (D-13) | Headless, no OpenGL → Windows-clean; **aligns** with NeuralEstimators' `Makie="0.24"` weakdep compat `[VERIFIED: NeuralEstimators Project.toml]`. |
| `Random` (stdlib) | bundled | Explicit `rng::AbstractRNG` threading (D-14) | stdlib; consistent with Phase-1 `Random.seed!` pattern. |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `HypothesisTests` | current | KS / χ² uniformity for SBC ranks (Phase 5) and optional KS distance for induced-μ match | Named in CLAUDE.md; optional in Phase 2 (a Wasserstein/KS distance is a clean SIM-02 pass metric) `[CITED: CLAUDE.md]`. |
| `ImageTransformations` | current (via `Images`) | `warp` + `Translation` for sub-pixel shift; `imresize` | If sub-pixel shift uses the `warp` idiom (recommended) rather than a hand-rolled Fourier shift. |
| `Interpolations` | current (via `Images`) | Interpolation backend for `warp` (`BSpline(Linear())`/`Cubic`) | Implicit dependency of the `warp` sub-pixel path. |
| `CoordinateTransformations` | current (via `Images`) | `Translation(dx, dy)` for the warp transform | Provides the affine shift; re-exported by `Images` `[ASSUMED — verify re-export]`. |
| `TiffImages` | current | OPTIONAL: load the real anchor TIFFs as sanity reference | Only if the plausibility step overlays real `test/test_images/*` patch-correlation; pure-Julia TIFF reader, Windows-clean. Prefer over ImageMagick. |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `warp` + `Translation` sub-pixel shift | Fourier-shift theorem (FFTW phase ramp) | Exact band-limited shift, but adds FFTW dep and complexity; `warp` with linear/cubic interpolation is the idiomatic Images.jl path and adequate for a registration-error nuisance. |
| `MvNormal` per-pixel correlated latent | Correlated Gaussian *random fields* (smooth, e.g. low-pass-filtered noise) | A shared smoothed latent gives biologically-plausible spatial structure AND controllable correlation; pure per-pixel MvNormal gives white-noise images (still valid for the contract but less realistic). Recommend correlated-then-smoothed fields. |
| `CairoMakie` | `Plots.jl` (GR backend) | Plots/GR is lighter to resolve but lower-quality vector output; CONTEXT locked CairoMakie (D-13). |
| `Images` umbrella | individual `ImageFiltering`+`ImageTransformations`+`Interpolations` | Umbrella is heavier but avoids version-skew between the small packages; either is fine — planner's call. |

**Installation (spike env only — never touch root):**
```bash
julia --project=spike -e 'using Pkg; Pkg.add(["Images","ImageFiltering","StatsBase","CairoMakie","HypothesisTests"])'
# Statistics + Random are stdlib (no add). Distributions already present.
# THEN: re-run spike/test/runtests.jl to confirm NeuralEstimators==0.2.1 & no CUDA, then commit the re-frozen Manifest.
```

**Version verification:** Julia packages resolve from the General registry; pin exact versions in `spike/Manifest.toml` after the green re-resolve (the Phase-1 reproducibility-artifact pattern, ENV-03). The binding constraint is `NeuralEstimators 0.2.1` → Makie 0.24 series; let the resolver pick the matching CairoMakie and **assert** post-resolve that NeuralEstimators did not move off 0.2.1.

## Package Legitimacy Audit

> Julia General registry; `slopcheck`/`npm` do not cover Julia. All packages are long-established, in-registry, and either already in the spike env or named in CLAUDE.md's verified stack. No hallucination vector — these are not discovered via web search but from the existing verified stack and direct `src/` import requirements.

| Package | Registry | Maturity | Source Repo | Disposition |
|---------|----------|----------|-------------|-------------|
| Distributions | Julia General | mature (JuliaStats), already pinned 0.25.128 | github.com/JuliaStats/Distributions.jl | Approved (in env) |
| Images / ImageFiltering / ImageTransformations / Interpolations | Julia General | mature (JuliaImages) | github.com/JuliaImages | Approved |
| StatsBase / Statistics | Julia General / stdlib | mature; **required** by `src/colocalization.jl` | github.com/JuliaStats/StatsBase.jl | Approved (mandatory) |
| CairoMakie | Julia General | mature (Makie org); compat-aligned with NeuralEstimators weakdep | github.com/MakieOrg/Makie.jl | Approved — planner adds post-resolve assertion |
| HypothesisTests | Julia General | mature (JuliaStats) | github.com/JuliaStats/HypothesisTests.jl | Approved (optional) |
| TiffImages | Julia General | mature (JuliaIO) | github.com/JuliaIO/TiffImages.jl | Approved (optional) |

**Packages removed due to slopcheck [SLOP]:** none (slopcheck N/A for Julia).
**Packages flagged [SUS]:** none.
**Resolve-risk gate (planner MUST add):** a `checkpoint`-style task that, after `Pkg.add`, asserts `Pkg.dependencies()` still shows `NeuralEstimators v0.2.1` and no installed `CUDA` node, mirroring `spike/test/runtests.jl` checks (b)/(c). If CairoMakie forces a NeuralEstimators downgrade, isolate plotting into a **separate** Julia env/process that reads simulator outputs from disk (JLD2/PNG) rather than co-resolving — the documented fallback.

## Architecture Patterns

### System Architecture Diagram

```
                         sample_prior(rng)                 [SIM-02]
                                │  draws μ* ~ Truncated(Cauchy(0,0.3),-1,1)
                                │  ρ_true = ĝ(μ*)   (frozen induced-μ inverse map)
                                ▼
   θ = (ρ_true, spillover, autofluor, label_eff, shift_dx, shift_dy, noise)   [D-04]
                                │
                                ▼   simulate_pair(rng, θ; imsize)             [SIM-01]
   ┌──────────────────────────────────────────────────────────────────────┐
   │ (1) correlated bivariate intensity fields  (MvNormal corr ρ_true,     │
   │      then shared spatial smoothing → nonneg transform)                 │
   │ (2) Bernoulli thinning        (label_eff = keep-probability)           │
   │ (3) Gaussian PSF              imfilter(·, Kernel.gaussian(σ_psf))  [D-05 fixed] │
   │ (4) 2×2 directional spillover  ch1' = ch1 + s·ch2 (directional ≠ symm) │
   │ (5) autofluorescence           additive offset per channel            │
   │ (6) sub-pixel shift            warp(ch2, Translation(dx,dy))           │
   │ (7) Poisson(shot) + Gaussian(read) noise   rand(Poisson(g·I)) + N(0,r) │
   └──────────────────────────────────────────────────────────────────────┘
                                │  data :: Vector{Matrix{Float64}} (2 ch)
                                ▼
            MultiChannelImage(data, names, ..., size(data[1]), otsu)   [SIM-03]
                                │  (UNCHANGED src/ summary — include() coupling)
                                ▼
        patch.(_,8) → correlation(_) → 8×8 per-patch ρ → mean = induced μ
                                │
              ┌─────────────────┴─────────────────┐
              ▼                                     ▼
     Plausibility gate (SIM-04)            Induced-μ calibration (SIM-02)
     • Spearman(ρ_true, mean ρ) ≥ thr      • sweep ρ_true → E[μ] curve, fit ĝ
     • spillover/shift paired effect       • KS/Wasserstein(induced μ, target)
     • valid MultiChannelImage             • freeze ĝ into sample_prior
              │                                     │
              ▼                                     ▼
     CairoMakie figures (D-13)            spike/NOTES.md §3 evidence
```

### Recommended Project Structure (`spike/`, layout = discretion)
```
spike/
├── simulator/
│   ├── prior.jl          # sample_prior(rng), the θ NamedTuple, ranges, frozen ĝ map
│   ├── forward.jl        # simulate_pair(rng, θ; imsize) — the 7-stage pipeline
│   └── calibration.jl    # offline ρ_true→μ sweep + fit ĝ + KS/Wasserstein evidence
├── contract.jl           # include()s src/colocalization.jl + src/LoadImages.jl; summary helpers
├── 02_simulator_demo.jl  # runnable demo: sample → simulate → summary → CairoMakie figs
├── figures/              # CairoMakie .png/.pdf outputs (gitignored or committed evidence)
├── NOTES.md              # §3 calibration mapping + evidence (D-02)
└── test/
    ├── runtests.jl       # extend: include test_simulator.jl
    └── test_simulator.jl # SIM-01..04 assertions, CPU-only, seeded
```

### Pattern 1: Controllable-correlation intensity fields
**What:** Generate two channels with a target pixel-level correlation ρ_true while retaining spatial structure.
**When to use:** Stage (1) of `simulate_pair`.
**Example:**
```julia
# Source: Distributions.jl MvNormal + ImageFiltering smoothing (idiom)  [ASSUMED]
using Distributions, ImageFiltering, Random
function correlated_fields(rng, W, H, ρ; struct_σ=4.0)
    Σ = [1.0 ρ; ρ 1.0]                       # 2×2 pixelwise correlation
    L = cholesky(Σ).L
    z1 = randn(rng, W, H); z2 = randn(rng, W, H)
    f1 = z1
    f2 = L[2,1].*z1 .+ L[2,2].*z2            # corr(f1,f2) = ρ exactly (Gaussian)
    # shared spatial smoothing preserves correlation, adds structure:
    k  = Kernel.gaussian(struct_σ)
    s1 = imfilter(f1, k); s2 = imfilter(f2, k)
    return softplus.(s1), softplus.(s2)      # nonneg intensities (distorts ρ → calibration)
end
```
Note: the nonneg transform + PSF + noise make the *induced* per-patch μ a monotone-nonlinear function of ρ_true — exactly why SIM-02 calibration is needed.

### Pattern 2: Fixed Gaussian PSF (D-05)
**What:** Microscopy point-spread blur as a fixed nuisance.
**Example:**
```julia
# Source: ImageFiltering.jl  [CITED: CLAUDE.md]
psf = Kernel.gaussian(σ_psf)        # σ_psf ≈ 1.0–1.5 px (see PSF sizing below)
ch  = imfilter(ch, psf)             # replicate-pad default; Windows-clean
```
**PSF sizing:** diffraction-limited FWHM ≈ 0.51·λ_em / NA. For λ_em≈520 nm, NA≈1.4 → FWHM≈190 nm. At Nyquist sampling (~65 nm/px, 100× oil) that is ≈3 px → σ = FWHM/2.355 ≈ **1.3 px**. For larger CCD pixels (the 1376×1028 anchor, coarser sampling) σ≈0.5–1 px. Recommend a fixed **σ_psf = 1.0–1.5 px** for the spike `[ASSUMED — discretion D-05]`.

### Pattern 3: Sub-pixel shift via `warp` (stage 6)
**What:** Channel-registration error of dx/dy fractional pixels.
**Example:**
```julia
# Source: ImageTransformations.warp + CoordinateTransformations.Translation  [ASSUMED — verify re-export]
using ImageTransformations, CoordinateTransformations, Interpolations
shifted = warp(ch2, Translation(dx, dy), axes(ch2);
               method=BSpline(Linear()), fillvalue=0.0)
```
Apply to **one** channel only (directional registration error). `fillvalue=0.0` keeps the contract (zeros are dropped by `_exclude_zero`, not NaN).

### Pattern 4: Poisson + Gaussian noise (D-07, stage 7)
```julia
# Source: Distributions.jl  [ASSUMED]
noisy = [rand(rng, Poisson(gain * max(I, 0.0))) / gain for I in ch] .+
        rand(rng, Normal(0.0, read_noise), size(ch))   # shot + read
```
Keep `gain` and `read_noise` in plausible microscopy ranges (discretion). Clamp to ≥0 if downstream wants non-negative intensities.

### Anti-Patterns to Avoid
- **Editing `src/` to make the contract fit.** Forbidden (SIM-03, decoupling). Adapt the simulator output to the contract, never the reverse.
- **Asserting ρ_true == induced μ.** They differ after the nonneg transform + PSF + noise; that gap is the whole point of D-01/D-02 calibration.
- **Pre-masking the simulated image.** The Turing path does not mask; masking would silently drop patches and break comparability.
- **Emitting zeros as "background."** `_exclude_zero` drops 0.0 pairs — large zero regions shrink the per-patch sample below 15 → spurious `missing`. Use small positive autofluorescence/offset for background, not exact zeros.
- **Resolving CairoMakie blindly.** Always re-assert NeuralEstimators 0.2.1 + no-CUDA after adding it; isolate plotting if it forces a downgrade.
- **Using `Random123` now.** Deferred to Phase 3 (D-14); use stdlib `Random` with an explicit `rng`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Per-patch Pearson correlation | A custom patch+correlation loop | `src/` `patch()` + `correlation()` via `include()` | SIM-03 *requires* the unmodified functions; re-implementing breaks comparability and the contract proof. |
| Gaussian convolution / PSF | Manual kernel + nested loops | `ImageFiltering.imfilter` + `Kernel.gaussian` | Edge handling, separability, performance all solved. |
| Sub-pixel translation | Manual bilinear resampling | `ImageTransformations.warp` + `Translation` | Interpolation, bounds, fill all handled. |
| Correlated bivariate draws | Hand-rolled correlation injection | `MvNormal` / Cholesky of the 2×2 Σ | Numerically exact target correlation. |
| Poisson / truncated-Cauchy / Exponential sampling | Inverse-CDF by hand | `Distributions.rand(rng, dist)` | Matches the Turing `@model` distributions exactly → prior consistency. |
| Monotone ρ_true→μ fit | Bespoke optimizer | Isotonic regression / PCHIP monotone spline (Interpolations / a small monotone fit) | Guarantees monotonicity, defensible, cheap. |
| Otsu threshold (if needed) | Histogram thresholding | `Images.otsu_threshold` | The `src/` constructor's own choice; keeps the object identical to real images. |

**Key insight:** Phase 2's value is *fidelity of the contract and the calibration*, not novel image-processing code. Every pixel op has a mature JuliaImages/JuliaStats primitive; the only genuinely new code is the physics composition, the prior, and the induced-μ fit.

## Induced-μ Calibration (SIM-02 / D-01, D-02) — the load-bearing method

**Goal:** make the *induced* distribution of per-patch-correlation mean μ (measured through the real `patch()`/`correlation()`) match the Turing μ-prior `Truncated(Cauchy(0,0.3),-1,1)`.

**Recommended method — monotone-inverse / transform sampling (defensible, deterministic):**

1. **Sweep.** Choose a grid of `ρ_true ∈ [ρ_min, ρ_max]` (≥ ~25 points). At each, draw `N` (e.g. 200–500) simulated pairs with the *other 6 nuisances drawn from their own prior* (so the induced μ reflects the real generative spread, not a single nuisance setting). Compute induced μ per pair via the real summary; record `E[μ | ρ_true]` and its conditional CDF.
2. **Fit ĝ.** Fit a **monotone** map between ρ_true and induced μ — isotonic regression or a PCHIP/monotone cubic on `(ρ_true, E[μ|ρ_true])`. Invert to get `ĝ : μ ↦ ρ_true`.
3. **Define `sample_prior`.** Draw `μ* ~ Truncated(Cauchy(0,0.3),-1,1)`, set `ρ_true = ĝ(μ*)`, draw the 6 nuisances from their priors → θ. By transform-sampling, the induced μ distribution matches the target by construction (up to the conditional spread, which the σ/τ mapping accounts for).
4. **Evidence artifacts (freeze in `spike/NOTES.md` §3 and CairoMakie figs):**
   - the `ρ_true → E[μ]` sweep curve with a CI band (monotonicity visible);
   - overlay histogram / QQ-plot of induced μ (after applying ĝ) vs the `Truncated(Cauchy(0,0.3),-1,1)` density;
   - a **quantitative match metric**: one-sample KS or Wasserstein-1 distance between the induced-μ sample and the target distribution, with a pre-declared tolerance as the SIM-02 pass bar.

**Alternative (simpler, also defensible):** direct quantile/histogram matching — simulate a large induced-μ pool under a base ρ_true prior, then remap via the empirical-CDF transform so the induced-μ CDF equals the target CDF. The monotone-inverse approach is cleaner because ĝ is reusable inside `sample_prior` and Phase 3's large-scale generator.

**σ / τ / ν consistency (D-03 — "where induced," not directly set):**
- **σ / τ** ↔ the spread of the 64 per-patch correlations within/across images. Finite per-patch pixel count and the noise level set this: the sampling SD of a patch Pearson r ≈ `(1−r²)/√n_eff`. Choose noise/label-efficiency ranges so the induced per-patch-correlation SD lands inside `Truncated(Cauchy(0.1,0.3),1e-4,1)`. **Verify** by measuring induced σ under the prior and checking range overlap; report in NOTES.
- **ν** (Student-t tail heaviness, `Exponential()`) ↔ heavy tails / outlier patches in the per-patch-correlation distribution (few-effective-pixel patches, spillover-driven outliers). Document that ν is *induced and checked*, not set: confirm the induced per-patch-correlation distribution is plausibly heavier-tailed than Gaussian. This is a consistency check, not a fit — the SIM-02 gate is on μ.

## Image-Size Sweep Set (D-09)

Survey of common scientific microscopy sensor formats (sCMOS/CCD/EMCCD, confocal + widefield) `[VERIFIED: WebSearch — 2048² sCMOS standard; MEDIUM for the smaller formats]`:
- **2048 × 2048** — modern sCMOS standard (Hamamatsu ORCA-Flash, Andor Zyla; 6.5 µm pixels), widefield + spinning-disk confocal.
- **1024 × 1024** — common laser-scanning-confocal default; EMCCD iXon 1024.
- **512 × 512** — fast/cropped confocal & EMCCD (iXon 512) — high frame-rate regime.
- **1376 × 1028** — the **paper anchor** (D-08), a CCD/CMOS sensor format (~Sony ICX285 class, 1376×1040), 16-bit.
- **256 × 256** — fast resonant-scanner confocal / heavy crop (smallest sensible; patch = 32×32 = 1024 px, still ≫ 15-pixel floor).

**Recommended sweep set:** `{256², 512², 1024², 1376×1028 (anchor), 2048²}`.

**Compute implications (8×8 grid fixed, D-10):** the `patch` array is `Array{Union{Float64,Missing},4}` ≈ `8·8·(W÷8)·(H÷8)` elements ≈ W·H elements per channel. At 2048²: ~4.2 M elements × ~16 bytes (Union boxing) ≈ ~67 MB per channel patch array, transient — fine on CPU. `correlation` is O(pixels) per patch (Pearson is single-pass). Per-pair summary at 2048² is sub-second; the sweep/calibration (thousands of sims) is the cost driver — keep calibration sims at a **smaller** size (e.g. 512²) and validate size-invariance separately. Non-multiples of 8 (1028→1024, 1376→1376) are silently trimmed by `patch`.

## Plausibility Gate (SIM-04 / D-12)

Three measurable sub-gates, all in the re-runnable `spike/test/` harness:

**(a) Monotonicity.** Sweep `ρ_true` over ≥15 points (N sims each, nuisances at prior medians to isolate ρ), compute `mean_patch_corr` per point, assert `corspearman(ρ_true_grid, mean_patch_corr) ≥ threshold`. **Recommended threshold: 0.95** (the construction is strongly monotone; a clean sim should easily clear it — leave headroom so noise doesn't flake the gate) `[ASSUMED — discretion]`. Use `StatsBase.corspearman`.

**(b) Perturbation visibility (spillover & shift).** With ρ_true fixed:
- *Spillover:* increasing the 2×2 directional spillover should **measurably raise** cross-channel mean patch correlation (bleed-through adds shared signal). Assert a paired difference `mean_corr(spill>0) − mean_corr(spill=0) > δ_spill` beyond Monte-Carlo noise (paired t / sign test over seeds).
- *Sub-pixel shift:* increasing `|shift|` should **measurably lower** mean patch correlation (decorrelates fine structure). Assert monotone decrease / `mean_corr(|shift|=0) − mean_corr(|shift|=k) > δ_shift`.
Quantify with an **effect size exceeding a tolerance**, not a visual check (D-12 "measurably").

**(c) Valid `MultiChannelImage` (SIM-03).** Constructor succeeds (2 channels, 2 otsu, `pixel_size == size(data[1])`), data are finite `Float64`, and `_prepare_data(stack,[1,2],8)` returns a non-empty `Vector{Vector{Float64}}` with ≤64 finite per-patch ρ.

**CairoMakie figure recipe (D-13):**
```julia
# Source: CairoMakie / Makie  [ASSUMED — verify save() formats on Windows]
using CairoMakie
CairoMakie.activate!()                       # headless, no OpenGL
fig = Figure(size=(1100, 800))
heatmap(fig[1,1], ch1_pos); heatmap(fig[1,2], ch2_pos)   # positive (high ρ) pair
heatmap(fig[1,3], ch1_neg); heatmap(fig[1,4], ch2_neg)   # negative (low ρ) pair
scatterlines(fig[2,1], ρ_grid, mean_patch_corr)          # (a) monotonicity sweep
hist(fig[2,2], induced_μ_samples; normalization=:pdf)    # (c) induced μ vs target overlay
scatterlines(fig[2,3], spill_grid, mean_corr_vs_spill)   # (b) spillover effect
scatterlines(fig[2,4], shift_grid, mean_corr_vs_shift)   # (b) shift effect
save(joinpath(@__DIR__,"figures","plausibility.png"), fig)  # raster
save(joinpath(@__DIR__,"figures","plausibility.pdf"), fig)  # vector
```
CairoMakie renders to PNG/PDF/SVG with **no display/OpenGL** → Windows-tauglich. `save()` infers format from extension.

## Common Pitfalls

### Pitfall 1: Zeros silently destroy patch correlations
**What goes wrong:** Large exact-zero "background" regions get dropped by `_exclude_zero`; patches fall below the 15-pixel floor → `missing` → induced μ computed from too few patches, unstable.
**Why:** `correlation` excludes 0.0 and NaN before `cor`.
**Avoid:** Give background a small positive autofluorescence/offset (stage 5) so pixels are non-zero; never emit hard zeros as background. Verify each patch keeps ≫15 pairs.
**Warning signs:** Many `missing` entries in the 8×8 ρ grid; induced μ jumps between sims.

### Pitfall 2: Confusing ρ_true with induced μ
**What goes wrong:** Asserting the gate against ρ_true directly, or skipping calibration.
**Why:** Nonneg transform + PSF + noise make induced μ a monotone-nonlinear, attenuated function of ρ_true.
**Avoid:** Always measure μ through the real summary; calibrate (D-01/D-02). Gate monotonicity on the *measured* μ, not ρ_true.

### Pitfall 3: CairoMakie / Makie resolve downgrading NeuralEstimators
**What goes wrong:** Adding a plotting stack drags an incompatible Makie/StructArrays in and silently downgrades the pinned NeuralEstimators 0.2.1 (the exact Phase-1 failure mode with the parent's GLMakie 0.10).
**Why:** Co-resolution picks the latest mutually-compatible versions; an old Makie series caps StructArrays < 0.7 which NeuralEstimators 0.2.1 needs.
**Avoid:** Use a **current** CairoMakie (Makie 0.24 series, which NeuralEstimators 0.2.1 *itself* declares as its weakdep compat). After `Pkg.add`, assert NeuralEstimators stays 0.2.1 and no CUDA, then re-freeze the Manifest. If it downgrades, isolate plotting into a separate env/process reading sim outputs from disk.
**Warning signs:** `Pkg.status` shows `NeuralEstimators v0.1.x`; the Phase-1 smoke errors on `NormalisingFlow(::Int; num_summaries=...)`.

### Pitfall 4: `pixel_size` set to micrometers
**What goes wrong:** Following the docstring (which says micrometers) instead of the code.
**Why:** `LoadImages.jl:444` sets `pixel_size = size(data[1])` — the (W,H) pixel tuple. D-11 makes this explicit.
**Avoid:** `pixel_size = size(data[1])` always.

### Pitfall 5: Sub-pixel shift introducing NaN/Inf at borders
**What goes wrong:** `warp` extrapolation produces NaN at out-of-bounds pixels → `_exclude_zero` keeps NaN? (it drops NaN, but NaN in `cor` of a too-small patch still risks instability).
**Why:** Default extrapolation can yield NaN outside the source axes.
**Avoid:** Pass `fillvalue=0.0` (dropped cleanly by `_exclude_zero`) and crop to the common valid region; never leave NaN borders large enough to thin a patch below 15.

## Code Examples

### End-to-end contract check (the SIM-03 proof, runnable in the test)
```julia
# include() coupling (Phase-1 proven path) — read-only, no src/ edits
include(joinpath(@__DIR__, "..", "src", "LoadImages.jl"))     # MultiChannelImage(Stack)
include(joinpath(@__DIR__, "..", "src", "colocalization.jl")) # patch, correlation, _prepare_data

data  = [ch1, ch2]                       # Vector{Matrix{Float64}}, 2 channels
mci   = MultiChannelImage(data, ["ch1","ch2"], "sim", ["",""],
                          size(data[1]), Images.otsu_threshold.(data))
stack = MultiChannelImageStack([mci], "sim")
ρvec  = _prepare_data(stack, [1, 2], 8)              # real summary, unchanged
@assert !isempty(ρvec) && length(ρvec[1]) ≤ 64
μ_ind = Statistics.mean(reduce(vcat, ρvec))          # induced μ
```

### Prior consistent with the Turing model (SIM-02 skeleton)
```julia
# Source: src/bayes.jl priors  [VERIFIED: src/bayes.jl:274-291]
using Distributions, Random
const MU_PRIOR = Truncated(Cauchy(0.0, 0.3), -1.0, 1.0)   # Turing μ-prior
# ĝ : μ ↦ ρ_true  — the frozen monotone inverse from the calibration sweep
function sample_prior(rng::AbstractRNG, ĝ)
    μ_star = rand(rng, MU_PRIOR)
    ρ_true = ĝ(μ_star)
    return (ρ_true        = ρ_true,
            spillover     = rand(rng, Uniform(0.0, 0.3)),   # ranges = discretion
            autofluor     = rand(rng, Uniform(0.0, 0.1)),
            label_eff     = rand(rng, Uniform(0.5, 1.0)),
            shift_dx      = rand(rng, Uniform(-1.0, 1.0)),
            shift_dy      = rand(rng, Uniform(-1.0, 1.0)),
            noise         = rand(rng, Uniform(0.0, 1.0)))
end
```

## State of the Art

| Topic | Current Approach | Note |
|-------|------------------|------|
| Microscopy forward models for SBI | Physics-based simulators (PSF + shot/read noise + spillover) feeding summary stats; calibrated to match a reference model's prior | Standard in amortized-SBI image pipelines; the induced-prior matching here is the rigorous form of "consistent prior." |
| Julia GPU-less imaging | Pure-Julia JuliaImages stack (no ImageMagick needed for synthetic arrays) | Windows-clean; TIFF reading (if used) via pure-Julia TiffImages, not ImageMagick. |
| Makie backends | Makie 0.24 series; CairoMakie for static/vector (headless) | Aligns with NeuralEstimators 0.2.1 weakdep; the old GLMakie 0.10/Makie 0.21 ecosystem is what broke Phase-1 resolve. |

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | A current CairoMakie (Makie 0.24) co-resolves with NeuralEstimators 0.2.1 without downgrade | Stack / Pitfall 3 | If it still conflicts, plotting must be isolated to a separate env/process — adds a task but not a blocker (mitigation documented). |
| A2 | `Images` umbrella re-exports `warp`, `Translation`, `Kernel`, `otsu_threshold` (so one dep suffices) | Stack | If not all re-exported, add `ImageTransformations`+`CoordinateTransformations`+`Interpolations` explicitly — trivial. |
| A3 | σ_psf ≈ 1.0–1.5 px is a realistic fixed PSF for the swept sizes | Pattern 2 | Too-wide PSF over-smooths fine structure → reduced dynamic range of induced μ; tune during calibration. Discretion (D-05). |
| A4 | Spearman monotonicity threshold 0.95 is achievable and not flaky | Plausibility Gate | If the clean sim cannot clear 0.95 under nuisance spread, lower to ~0.9; discretion (D-12). |
| A5 | Poisson+Gaussian constants (gain, read_noise) in chosen ranges induce σ/τ within the Turing scale-prior ranges | Calibration D-03 | If induced σ falls outside `Truncated(Cauchy(0.1,0.3),1e-4,1)`, retune noise/label-efficiency ranges; consistency check, not a hard gate. |
| A6 | The monotone-inverse calibration (transform sampling) makes induced μ match the target up to conditional spread | Calibration | If the ρ_true→μ relation is non-monotone in some regime, isotonic regression still yields a usable map; falls back to quantile matching. |
| A7 | `warp` with `fillvalue=0.0` keeps the contract clean (zeros dropped, no NaN) | Pattern 3 / Pitfall 5 | If borders introduce NaN, crop to valid region; trivial mitigation. |

## Open Questions (RESOLVED)

1. **Exact ρ_true sweep range for calibration.**
   - **RESOLVED:** sweep ρ_true ∈ [−0.99, 0.99]; calibration.jl extends if induced μ saturates before the ±0.9 tails (documented in NOTES §3).
   - What we know: induced μ must cover `Truncated(Cauchy(0,0.3),-1,1)` (heavy mass near 0, tails to ±1).
   - What's unclear: the ρ_true span needed to induce μ out to ±0.9 after attenuation by PSF/noise.
   - **Recommendation:** sweep ρ_true ∈ [−0.99, 0.99]; extend if induced μ saturates before reaching the target tails. Determine empirically in the first calibration run; document the chosen span in NOTES §3.

2. **Whether to calibrate at the anchor size (1376×1028) or a smaller size.**
   - **RESOLVED:** calibrate at 512², then verify induced-μ size-invariance across the sweep set {256²,512²,1024²,1376×1028,2048²}.
   - What we know: calibration needs thousands of sims; large sizes are slower.
   - What's unclear: whether induced μ is size-invariant (it should be — Pearson is sample-size-robust above the 15-px floor).
   - **Recommendation:** calibrate at 512², then **verify size-invariance** of induced μ across the sweep set (a cheap check); if invariant, the 512² map transfers to all sizes.

3. **Does the negative test image actually exhibit low patch-correlation (sanity anchor)?**
   - **RESOLVED:** OPTIONAL overlay of test/test_images/{positive,negative} induced-μ; sanity-only, not a gate.
   - What we know: `test/test_images/{positive,negative}` are the paper's ground-truth coloc/non-coloc pairs.
   - What's unclear: their exact mean patch-correlation values.
   - **Recommendation:** OPTIONAL — load them via TiffImages, run the real summary, and use their (μ_pos, μ_neg) as plausibility anchors for the simulator's μ range. Not a gate; a nice sanity overlay for the SIM-04 figure.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Julia | everything | ✓ (Phase 1) | 1.12.6 (juliaup override on `spike/`) | — |
| spike env (`spike/Project.toml`) | isolation | ✓ | NeuralEstimators 0.2.1, Flux 0.16.10, Distributions 0.25.128 | — |
| Distributions.jl | priors, noise, correlated densities | ✓ | 0.25.128 (pinned) | — |
| Images / ImageFiltering / ImageTransformations | PSF, sub-pixel shift, otsu | ✗ (to add) | Julia General current | manual conv/shift (not recommended) |
| StatsBase + Statistics | `src/colocalization.jl` deps, Spearman | ✗ StatsBase (to add) / ✓ Statistics (stdlib) | current / bundled | none — StatsBase is **mandatory** for the `include()` |
| CairoMakie | SIM-04 figures (D-13) | ✗ (to add) | Makie 0.24 series | Plots/GR, or isolate to separate env if resolve conflicts |
| HypothesisTests | optional KS for induced-μ match | ✗ (to add) | current | StatsBase/manual KS statistic |
| `src/` files via include() | summary contract | ✓ (Phase-1 proven) | baseline f581d95 | — |
| TiffImages | OPTIONAL real-image sanity anchor | ✗ | current | skip the real-image overlay |

**Missing dependencies with no fallback:** `StatsBase` (required by the `src/colocalization.jl` `include()` — without it the summary contract cannot load).
**Missing dependencies with fallback:** `CairoMakie` (isolate to separate env if it conflicts), `HypothesisTests` (manual KS), `Images` family (manual ops, discouraged).

**Decoupling guard (carry-forward):** all adds touch **only** `spike/Project.toml` + `spike/Manifest.toml`. After resolve, re-assert `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` is clean (the Phase-1 / DEMO-02 proof command).

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | stdlib `Test` (`@testset`), extending Phase-1 `spike/test/runtests.jl` |
| Config file | none — single runner (`spike/test/runtests.jl`), Phase-1 pattern |
| Quick run command | `julia --project=spike spike/test/runtests.jl` (subset: small size, few sweep points) |
| Full suite command | `julia --project=spike spike/test/runtests.jl` (full sweep + calibration check) |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| SIM-03 | Sim output builds a valid `MultiChannelImage`; `_prepare_data` returns ≤64 finite per-patch ρ | unit | `julia --project=spike spike/test/runtests.jl` | ❌ Wave 0 (`test_simulator.jl`) |
| SIM-01 | `simulate_pair(rng,θ)` runs all 7 stages, returns 2× `Matrix{Float64}`; θ is the 7-tuple | unit | same | ❌ Wave 0 |
| SIM-02 | Induced-μ KS/Wasserstein distance to `Truncated(Cauchy(0,0.3),-1,1)` < tol; ĝ monotone | integration | same (full) | ❌ Wave 0 |
| SIM-04a | Spearman(ρ_true, mean patch-corr) ≥ threshold over the sweep | integration | same (full) | ❌ Wave 0 |
| SIM-04b | Spillover raises & sub-pixel shift lowers mean patch-corr beyond MC noise (paired) | integration | same (full) | ❌ Wave 0 |
| (guard) | NeuralEstimators stays 0.2.1, no CUDA after env change | unit | reuse Phase-1 checks (b)/(c) | ✅ exists in `runtests.jl` |

### Sampling Rate
- **Per task commit:** quick subset — 1 sim at 256², a 5-point monotonicity check, contract assertion (seconds).
- **Per wave merge:** full sweep (≥15 points) + induced-μ KS at calibration size + perturbation tests.
- **Phase gate:** full suite green + CairoMakie figures regenerated, before `/gsd:verify-work`.

### Wave 0 Gaps
- [ ] `spike/test/test_simulator.jl` — SIM-01..04 assertions (included from `runtests.jl`).
- [ ] Add `StatsBase`, `Images`/`ImageFiltering`, `CairoMakie` (+ optional `HypothesisTests`) to `spike/Project.toml`; re-resolve; re-freeze Manifest; re-run Phase-1 smoke as the regression guard.
- [ ] `spike/simulator/{prior,forward,calibration}.jl` + `spike/contract.jl` — the units under test.
- [ ] Seed strategy: explicit `rng = Random.Xoshiro(seed)` threaded through (D-14), so the gate is deterministic & CPU-only.

*(Existing infra reused: the Phase-1 `runtests.jl` testset harness and the CUDA-absence / NeuralEstimators-pin checks.)*

## Security Domain

> `security_enforcement` not set in config → treated as enabled. This is a **local, offline scientific simulation spike** — no auth, no network, no untrusted input, no persistence of secrets. Most ASVS categories are N/A.

### Applicable ASVS Categories
| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | No auth surface (local Julia script). |
| V3 Session Management | no | No sessions. |
| V4 Access Control | no | No multi-user surface. |
| V5 Input Validation | partial | Validate θ ranges / `imsize ≥ 8` in `simulate_pair`; assert finite, non-negative intensities before constructing `MultiChannelImage` (also a correctness gate). |
| V6 Cryptography | no | No crypto; RNG is for reproducibility (stdlib `Random`), not security — do not treat the seed as a secret. |

### Known Threat Patterns for this stack
| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Out-of-range θ (e.g. ρ outside [−1,1], negative size) producing NaN/garbage images | Tampering (data integrity) | Argument validation + `@assert` finite/in-range at simulator entry; covered by the SIM-01/03 tests. |
| Accidental `src/` mutation breaking the decoupling guarantee | Tampering | `include()` is read-only; re-run the `git diff --quiet f581d95 -- src/` proof after the phase. |
| Supply-chain (new spike deps) | Tampering | All adds are mature JuliaStats/JuliaImages/Makie packages from the General registry; pin exact versions in the committed Manifest (ENV-03 pattern). |

## Sources

### Primary (HIGH confidence)
- `src/colocalization.jl` (read this session) — `patch` (37/78), `correlation` (221), `_prepare_data` (bayes.jl 164) signatures, 4-D array shapes, `_exclude_zero` 15-pixel floor, element types.
- `src/bayes.jl` (read this session, lines 130–329) — Turing `@model` priors `μ~Truncated(Cauchy(0,0.3),-1,1)`, `ν~Exponential()`, `σ/τ~Truncated(Cauchy(0.1,0.3),1e-4,1)`; `_prepare_data` flow.
- `src/LoadImages.jl` (read this session) — `MultiChannelImage`/`MultiChannelImageStack` structs + constructor asserts; `pixel_size = size(data[1])` (444); `otsu_threshold` usage.
- `spike/NOTES.md` (read this session) — include() coupling path, StatsBase+Statistics requirement, Turing prior seed, decoupling proof command (f581d95).
- `spike/00_smoke.jl`, `spike/test/runtests.jl` — Phase-1 test harness pattern, CUDA/NeuralEstimators pin checks, AGPL header, stdlib `Random` seeding.
- NeuralEstimators.jl `Project.toml` (github.com/msainsburydale/NeuralEstimators.jl) — weakdeps Makie/Plots/ColorSchemes; `compat Makie="0.24"` → CairoMakie resolve-alignment finding.

### Secondary (MEDIUM confidence)
- WebSearch (microscopy sensor formats) — 2048² sCMOS standard; corroborates the image-size sweep set; smaller formats from domain knowledge.
- CLAUDE.md Technology Stack — Images.jl / ImageFiltering.jl / Distributions.jl / CairoMakie / HypothesisTests as the named, prior-researched stack.

### Tertiary (LOW confidence)
- Distributions.jl / ImageFiltering.jl / ImageTransformations `warp` API details — from training knowledge (stable, mature APIs); marked `[ASSUMED]` and to be confirmed at implementation against installed docstrings.

## Metadata

**Confidence breakdown:**
- Summary contract (SIM-03): **HIGH** — verified by direct source read of `patch`/`correlation`/`_prepare_data`/constructor.
- Standard stack: **HIGH** — all packages named in CLAUDE.md or required by the `include()`; CairoMakie resolve-alignment verified against NeuralEstimators Project.toml.
- Forward pipeline (SIM-01): **MEDIUM-HIGH** — stage decomposition locked by D-06; exact API calls from mature-but-training-sourced libs (`[ASSUMED]`, verify at install).
- Induced-μ calibration (SIM-02): **MEDIUM** — method is defensible and concrete; the exact fit and ρ_true span are empirical (resolved in the first calibration run).
- Plausibility gate (SIM-04): **HIGH** — quantitative bar fully specified; CairoMakie save path standard.

**Research date:** 2026-06-26
**Valid until:** ~2026-07-26 (stable Julia ecosystem; re-verify CairoMakie/NeuralEstimators co-resolve if either is bumped).
