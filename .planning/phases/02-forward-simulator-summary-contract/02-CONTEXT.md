# Phase 2: Forward Simulator + Summary Contract - Context

**Gathered:** 2026-06-26
**Status:** Ready for planning

<domain>
## Phase Boundary

Build a 2-channel 2D **physics forward-simulator** that emits real `MultiChannelImage` pairs from a physics parameter vector θ, such that the **unmodified** existing summary functions (`correlation()` / `patch()` / `_prepare_data()`) ingest them directly, with a prior π(θ) **provably consistent** with the existing Turing `@model`, validated by plausibility plots.

In scope: `simulate_pair(θ)`, `sample_prior()` / π(θ), the prior-consistency calibration, and the plausibility validation. **Out of scope:** large-scale training-data generation (Phase 3), NPE/NRE training and inference (Phase 4+), Bayes-factor reproduction (Phase 6). No `src/` edits — the simulator *uses* the existing summary functions via the Phase-1 `include()` coupling.

Requirements: SIM-01, SIM-02, SIM-03, SIM-04.
</domain>

<decisions>
## Implementation Decisions

### Prior consistency — the ρ_true ↔ μ mapping (SIM-02)
- **D-01:** **Full induced-μ calibration.** Fit/transform the `ρ_true` prior so the simulator's **induced** distribution of the per-patch correlation mean (μ, the quantity the Turing model infers) matches the Turing μ-prior `Truncated(Cauchy(0, 0.3), -1, 1)` as closely as possible — not merely adopting the range by assertion. This is the rigorous "provably consistent" path required so the ADVI benchmark (infers μ) and the NPE/BF validation (infers ρ_true / Δρ) share **one** generative prior.
- **D-02:** The calibration is a concrete step: sweep `ρ_true` → measure the induced per-patch-correlation μ (via the real `patch()`/`correlation()` summary on simulated pairs) → fit the transform/prior so the induced μ distribution matches the Turing μ-prior. The fitted mapping AND the calibration evidence (the ρ_true→μ curve, the induced-vs-target μ overlap) are documented in `spike/NOTES.md` (the file Phase 1 seeded with the μ/ν/σ/τ ranges for exactly this purpose).
- **D-03:** The noise/scatter nuisances of θ map to the Turing scale hyperpriors: per-patch correlation spread ↔ `σ ~ Truncated(Cauchy(0.1,0.3),1e-4,1)` / `τ`, and tail/heaviness ↔ `ν ~ Exponential()`. Keep these consistent with the Turing ranges where the simulator induces them; exact mapping form is research/planning detail.

### Physics structure (SIM-01)
- **D-04:** θ is the **7-parameter** vector exactly as the roadmap lists: `θ = (ρ_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise)`. (`shift` counts as dx/dy.)
- **D-05:** **PSF is a fixed Gaussian nuisance** — NOT an inferred θ-component. A realistic Gaussian PSF width is fixed for the spike (exact width tied to pixel scale = Claude's discretion / researcher).
- **D-06:** Pipeline order (per ROADMAP SC-1): correlated bivariate densities → **Bernoulli thinning** (label efficiency = Bernoulli keep-probability) → **Gaussian PSF** convolution (`ImageFiltering.imfilter`, `Kernel.gaussian`) → **2×2 directional spillover** mixing (bleed-through, directional not symmetric) → **autofluorescence** offset → **sub-pixel shift** (dx/dy registration error) → **Poisson (shot) + Gaussian (read) noise**.
- **D-07:** Noise is the **combined Poisson + Gaussian** model (shot + read), not Gaussian-only — chosen for fidelity over per-sim speed.

### Image geometry
- **D-08:** **Anchor on the real paper data size: 1376 × 1028, 16-bit** (the actual *Scientific Reports* images at `test/test_images/{positive,negative}/*_c{1,2,3}.tif`). The simulated pairs should match this native regime so the 8×8 summary is comparable to real inputs.
- **D-09:** **Also sweep/investigate other sensible microscopy sizes** — image dimensions are a **configurable** parameter, and Phase 2 evaluates a range of common confocal/widefield dimensions (the exact set — e.g. 512², 1024², 1376×1028, 2048² — is a research item: survey typical microscopy image sizes). Goal: confirm the simulator + summary behave across the realistic size range, not just one size.
- **D-10:** **The 8×8 patch grid stays FIXED** for the spike (CLAUDE.md constraint), feeding the unchanged `patch()` / `correlation()`.
- **D-11:** The `MultiChannelImage.pixel_size` field carries the **(width, height) pixel-dimension tuple**, matching the existing `src/LoadImages.jl:444` convention (`pixel_size = size(data[1])`) — it is NOT micrometers despite the docstring. The simulator must follow this so output verifies as a valid `MultiChannelImage`.

### Validation + plotting (SIM-04)
- **D-12:** **Quantitative plausibility gate**, not visual-only. SIM-04 passes only if: (a) a `ρ_true` sweep shows **monotone** ρ_true↑ → mean patch-correlation↑ (Spearman rank correlation ≥ a threshold; threshold = Claude's discretion), AND (b) **spillover and sub-pixel shift measurably/visibly perturb** the pair, AND (c) output verifies as a valid `MultiChannelImage` (SIM-03). Mirrors Phase 1's green+correctness gate philosophy.
- **D-13:** **Plotting backend for the spike's plausibility figures = CairoMakie** — lightweight, static, headless, vector output; near drop-in from the package's existing GLMakie idiom. **Spike-local only** (under `spike/`), zero `src/` risk. This does **NOT** resolve the main-package `src/plot.jl` GLMakie→? backend question, which remains a separate, deferred, post-spike decision.

### Reproducibility / seeding
- **D-14:** The simulator accepts an **explicit `rng::AbstractRNG`** argument (RNG threaded through `sample_prior` / `simulate_pair`) so the seeding strategy is pluggable. Phase 2 uses stdlib `Random` (consistent with Phase 1); **Random123 counter-based seeding is adopted in Phase 3** for reproducible large-scale/parallel simulation (Phase 1 deferred Random123; also addresses the Phase-1 code-review thread-reproducibility note WR-02).

### Claude's Discretion
- Exact fixed PSF Gaussian width (tie to pixel scale).
- Exact Poisson gain + Gaussian read-noise constants (plausible microscopy ranges).
- The Spearman threshold for the SIM-04 monotonicity gate.
- The exact calibration-fit method for D-01/D-02 (e.g. regression or quantile-matching of ρ_true→μ).
- The exact set of "sensible microscopy sizes" to sweep (D-09) — researcher surveys typical confocal/widefield dimensions.
- Internal simulator file/module layout within `spike/`.
</decisions>

<specifics>
## Specific Ideas

- Real anchor data is in-repo: `test/test_images/{positive,negative}/*_c{1,2,3}.tif`, native **1376×1028, 16-bit, 3-channel**. The positive/negative pair are the package's own colocalized / non-colocalized ground-truth examples — useful sanity references for "what realistic patch-correlation looks like".
- The Turing model does NOT model images — it models the **per-patch correlation summary** as a Student-t sample: `control[idx]/sample[idx] ~ TDist(ν)·σ + μ` (`src/bayes.jl:296,301`). This is WHY SIM-02 is a mapping problem (physics θ → induced summary-stat μ), not a copy of ranges.
- Carry-forward (Phase 1): spike isolation is a hard guarantee — decoupling baseline is commit `f581d95`; nothing under `src/` or the root manifests may change. The simulator reaches `src/` functions via the D-01 `include()` coupling proven in Phase 1.
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & requirements
- `.planning/ROADMAP.md` §"Phase 2: Forward Simulator + Summary Contract" — goal + 4 success criteria (the simulator pipeline is spelled out in SC-1).
- `.planning/REQUIREMENTS.md` — SIM-01, SIM-02, SIM-03, SIM-04.

### The consistency target (SIM-02)
- `src/bayes.jl` §`@model function model` (≈ lines 251–316) — the Turing model: per-patch-correlation Student-t likelihood (`TDist(ν)·σ + μ`) and the μ/ν/σ/τ hyperpriors (`μ ~ Truncated(Cauchy(0,0.3),-1,1)`, `ν ~ Exponential()`, `σ/τ ~ Truncated(Cauchy(0.1,0.3),1e-4,1)`). This is the prior π(θ) must be made consistent with.
- `spike/NOTES.md` — Phase-1 seed of the Turing μ/ν/σ/τ ranges + the `include()`-coupling path to reach `src/` functions; this is where the SIM-02 calibration mapping + evidence get documented.

### The summary contract (SIM-03)
- `src/colocalization.jl` — `patch()` (lines 37, 78), `correlation(x::Array{T,4}, y::Array{T,4}; method)` (line 221), `_prepare_data()`. Applied UNCHANGED to simulator output; the 8×8 patch grid feeds these.
- `src/LoadImages.jl` — `MultiChannelImage` struct + constructor (lines 118–130); note `pixel_size = size(data[1])` convention (line 444). Defines what a valid simulated image must satisfy.

### Constraints
- `CLAUDE.md` §Constraints + §Technology Stack — priors-consistency mandate (π(θ) ↔ Turing μ/ν/σ/τ), fixed 8×8 summary during the spike, reproducibility-from-seed, the Images.jl / ImageFiltering.jl / CairoMakie supporting libs, and the "What NOT to Use" list.
- `.planning/phases/01-environment-smoke-gate/01-CONTEXT.md` + `01-BASELINE.md` — Phase-1 decisions: spike isolation, decoupling baseline `f581d95`, the D-01 `include()` coupling the simulator relies on.
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets (used UNCHANGED — no src/ edits)
- `correlation()` / `patch()` / `_prepare_data()` (`src/colocalization.jl`) — the summary-statistic the simulator output must feed; `correlation` takes 4-D arrays, `patch` builds the patch grid.
- `MultiChannelImage` constructor (`src/LoadImages.jl:126`) — signature `(data::Vector{Matrix{Union{Missing,Float64}}}, channels::Vector{String}, name::String, path::Vector{String}, pixel_size::Tuple{Int,Int}, otsu_threshold::Vector{<:AbstractFloat})`; constructor asserts `length(data)==length(channels)==length(otsu_threshold)`.
- Turing `@model` + μ/ν/σ/τ priors (`src/bayes.jl`) — the SIM-02 consistency target.
- `spike/NOTES.md` — already holds the prior-range seed from Phase 1.

### Established Patterns
- Spike reaches `src/` via `include()` (Phase-1 D-01 fallback — `Pkg.develop` downgrades NeuralEstimators, so include() is the proven path).
- Spike-local figures via CairoMakie; reproducible-from-seed with an explicitly-threaded RNG.
- Validation = quantitative gate (Phase-1 green+correctness philosophy), re-runnable under the spike `Test` harness.

### Integration Points
- Simulator output must construct a valid `MultiChannelImage` (2 channels): `data` = `Vector{Matrix{Float64}}`, `channels` = 2 names, `pixel_size` = the (W,H) pixel-dim tuple, `otsu_threshold` = per-channel (compute via existing `otsu_thresholds` or supply). Real anchor regime: 1376×1028, 16-bit (from `test/test_images/`).
- ImageFiltering.jl (`imfilter`, `Kernel.gaussian`) for the PSF + sub-pixel shift; Distributions.jl for the bivariate correlated densities + Poisson/Gaussian noise.
</code_context>

<deferred>
## Deferred Ideas

- **Random123 counter-based seeding** — adopted in Phase 3 (training-data scale / parallel sims), not Phase 2 (D-14).
- **The `src/plot.jl` GLMakie → CairoMakie/Plots backend migration** for the MAIN package — a separate post-spike decision; out of scope here (the spike's CairoMakie choice D-13 is independent and touches no `src/`).
- **Training-data generation at 50k–200k scale, NPE/NRE training, ADVI baseline, Bayes-factor reproduction** — Phases 3, 4, 6 respectively.

### Reviewed Todos (not folded)
None — no pending todos matched this phase.
</deferred>

---
*Phase: 02-forward-simulator-summary-contract*
*Context gathered: 2026-06-26*
