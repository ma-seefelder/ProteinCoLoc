# Phase 2: Forward Simulator + Summary Contract - Pattern Map

**Mapped:** 2026-06-26
**Files analyzed:** 9 (7 new, 2–3 modified)
**Analogs found:** 7 / 9 (2 stages have NO in-repo analog — flagged explicitly)

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `spike/simulator/prior.jl` | model (prior π(θ)) | transform / sampling | `src/bayes.jl` Turing priors (μ/ν/σ/τ) + `spike/NOTES.md §3` | role-match (the prior is the consistency *target*, not a copyable sampler) |
| `spike/simulator/forward.jl` | service (physics simulator) | transform (θ → image) | **NO in-repo analog** — see RESEARCH §Architecture Patterns 1–4 | none (new physics code) |
| `spike/simulator/calibration.jl` | service (offline fit) | batch / sweep | `spike/00_smoke.jl` (offline-script + correctness-gate shape) | partial (structure only) |
| `spike/contract.jl` | utility (include() coupling + summary helpers) | request-response (calls src/) | `spike/NOTES.md §4` include() recipe + `src/colocalization.jl` / `src/LoadImages.jl` call signatures | exact (coupling pattern proven in Phase 1) |
| `spike/02_simulator_demo.jl` | config / runnable demo | batch | `spike/00_smoke.jl` (AGPL header, `Random.seed!`, CPU-only top-level script) | exact (script shape) |
| `spike/test/test_simulator.jl` | test | request-response | `spike/test/runtests.jl` (`@testset verbose=true`, assert-tolerance gate) | exact (harness shape) |
| `spike/test/runtests.jl` (MODIFY) | test (harness) | — | itself (Phase 1) — append `include("test_simulator.jl")`, PRESERVE CUDA/NeuralEstimators guards | exact |
| `spike/Project.toml` / `spike/Manifest.toml` (MODIFY) | config | — | `spike/Project.toml` (minimal-deps `[deps]` block) | exact |
| `spike/NOTES.md` (MODIFY) | docs | — | `spike/NOTES.md §3` (pre-seeded Turing ranges) — append §3 calibration evidence | exact |
| `spike/figures/*.png/.pdf` (NEW output) | artifact | file-I/O | **NO in-repo analog** — see RESEARCH §CairoMakie recipe | none (CairoMakie new to repo) |

## Pattern Assignments

### `spike/02_simulator_demo.jl` + `spike/simulator/*.jl` headers (offline script shape)

**Analog:** `spike/00_smoke.jl`

**AGPL header + CPU-only seeded-script pattern** — every new `spike/*.jl` file MUST open with the verbatim `#= … =#` AGPL block (lines 1-19 of `spike/00_smoke.jl`), then a `# spike/<file>.jl --- <purpose>` comment, then `using …, Random`, then `Random.seed!(<seed>)`. The smoke uses:
```julia
using NeuralEstimators, Flux, Distributions, Random
Random.seed!(2026)   # determinism: stdlib Random (Random123 is a later phase)
```
For Phase 2 mirror this but thread an **explicit `rng`** (D-14) rather than only the global seed:
```julia
using Distributions, Images, ImageFiltering, StatsBase, Statistics, Random
rng = Random.Xoshiro(2026)   # explicit rng threaded into sample_prior / simulate_pair
```
Keep the same `@assert <condition> "msg"` + `println("... OK: ...")` correctness-gate idiom (lines 77-81) so the demo self-checks.

---

### `spike/test/test_simulator.jl` (new) + `spike/test/runtests.jl` (modify)

**Analog:** `spike/test/runtests.jl`

**Harness shape to mirror** (lines 40-61): `using Test`, top-level `@testset "<name>" verbose = true begin … end`, nested `@testset` per requirement, `include(joinpath(@__DIR__, "..", "<script>.jl"))` to pull script-defined bindings into scope:
```julia
using Test
@testset "Phase 2 simulator" verbose = true begin
    include(joinpath(@__DIR__, "..", "02_simulator_demo.jl"))  # defines bindings under test
    @testset "SIM-03 valid MultiChannelImage" begin
        @test !isempty(ρvec) && length(ρvec[1]) ≤ 64
    end
    @testset "SIM-04a monotonicity" begin
        @test corspearman(ρ_grid, mean_patch_corr) ≥ threshold
    end
end
```

**CRITICAL — PRESERVE the Phase-1 guard tests** (lines 52-59). The CUDA-absence + NeuralEstimators-pin `@testset` is the resolve-risk gate RESEARCH §Resolve-risk demands after adding Images/CairoMakie. Do NOT delete or weaken it; `runtests.jl` should keep running it AND newly `include("test_simulator.jl")`:
```julia
@testset "D-04 CPU-only (no CUDA dependency, none loaded)" begin
    @test !haskey(Pkg.project().dependencies, "CUDA")
    @test !any(p -> occursin("CUDA", p.name), values(Pkg.dependencies()))
    @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
end
```
RESEARCH adds the explicit post-resolve assertion that NeuralEstimators stays `v0.2.1` — add as a sibling `@test` in the same block (e.g. `@test Pkg.dependencies()[NE_uuid].version == v"0.2.1"`).

---

### `spike/contract.jl` (new — the include() coupling boundary)

**Analog:** `spike/NOTES.md §4` (proven Phase-1 fallback) + `src/colocalization.jl` / `src/LoadImages.jl` call signatures.

**include() coupling — read-only, NEVER edit src/** (verbatim recipe from `spike/NOTES.md` lines 142-145 and RESEARCH §Code Examples):
```julia
include(joinpath(@__DIR__, "..", "src", "LoadImages.jl"))     # MultiChannelImage(Stack)
include(joinpath(@__DIR__, "..", "src", "colocalization.jl")) # patch, correlation, _prepare_data
```
StatsBase + Statistics MUST be loaded before the `include` (`src/colocalization.jl` references `corspearman`/`corkendall`/`mean`/`median`/`quantile`) — see `spike/NOTES.md` lines 147-149.

**MultiChannelImage construction — mirror the convenience-constructor field tuple** (`src/LoadImages.jl:444-454`). The struct constructor (lines 126-130) asserts `length(data)==length(channels)==length(otsu_threshold)`:
```julia
# data :: Vector{Matrix{Float64}}, exactly 2 channels (Float64 <: Union{Missing,Float64})
mci = MultiChannelImage(
    data, ["ch1","ch2"], "sim", ["",""],
    size(data[1]),                       # pixel_size = (W,H) tuple — D-11, NOT micrometers
    Images.otsu_threshold.(data))        # per-channel Float; mirrors src/LoadImages.jl:445
```
Note constructor's `I <: Int` constraint (line 126): `pixel_size` must be `Tuple{Int,Int}` — `size(data[1])` already returns this.

**Summary call sequence — invoke UNCHANGED, do NOT reimplement** (`src/colocalization.jl`):
- `patch(img::AbstractMatrix{T}, 8)` (line 37) → `Array{Union{Float64,Missing},4}` shape `(8,8,W÷8,H÷8)`; non-multiples of 8 silently trimmed (lines 38-43).
- `correlation(x::Array{T,4}, y::Array{T,4}; method=:pearson)` (line 221) → `Matrix{Union{Float64,Missing}}` (8×8). Per patch: `_exclude_zero(a,b)` then `length(a) <= 15 ? missing : cor(a,b)` (lines 234-235).
- Induced μ = `Statistics.mean(skipmissing(ρ))` / `mean(reduce(vcat, ρvec))`.

End-to-end contract exercise (RESEARCH §Code Examples, the SIM-03 proof):
```julia
stack = MultiChannelImageStack([mci], "sim")
ρvec  = _prepare_data(stack, [1, 2], 8)     # real summary, unchanged
μ_ind = Statistics.mean(reduce(vcat, ρvec))
```

---

### `spike/simulator/prior.jl` (new)

**Analog:** `src/bayes.jl` Turing priors (the consistency *target*) + `spike/NOTES.md §3` (pre-seeded ranges).

**The Turing μ/ν/σ/τ ranges to stay consistent with** — already copied into `spike/NOTES.md` lines 188-193, sourced from `src/bayes.jl` ~274-291:
```julia
mu    ~ Truncated(Cauchy(0.0, 0.3), -1.0, 1.0)
nu    ~ Exponential()
sigma ~ Truncated(Cauchy(0.1, 0.3), 1e-4, 1.0)
tau   ~ Truncated(Cauchy(0.1, 0.3), 1e-4, 1.0)
```

**Prior skeleton — explicit rng, frozen ĝ map** (D-01/D-02/D-14; RESEARCH §Code Examples). The 7-param θ (D-04) is sampled by drawing μ* from the Turing μ-prior then mapping ρ_true = ĝ(μ*):
```julia
const MU_PRIOR = Truncated(Cauchy(0.0, 0.3), -1.0, 1.0)   # = src/bayes.jl μ-prior
function sample_prior(rng::AbstractRNG, ĝ)
    μ_star = rand(rng, MU_PRIOR)
    ρ_true = ĝ(μ_star)                     # frozen monotone inverse from calibration sweep
    return (ρ_true=ρ_true, spillover=…, autofluor=…, label_eff=…,
            shift_dx=…, shift_dy=…, noise=…)   # nuisance ranges = discretion
end
```
**Trap (do NOT assert ρ_true == μ):** the nonneg transform + PSF + noise make induced μ a monotone-nonlinear, attenuated function of ρ_true — that gap is the whole reason calibration exists (RESEARCH Pitfall 2 / Anti-Patterns).

---

### `spike/simulator/forward.jl` (new) — NO in-repo analog

**Analog:** **NONE.** The 7-stage physics pipeline (correlated bivariate fields → Bernoulli thinning → Gaussian PSF → 2×2 directional spillover → autofluorescence → sub-pixel shift → Poisson+Gaussian noise) has no existing analog in `src/` or `spike/`. The executor MUST follow **RESEARCH §Architecture Patterns 1–4** (the verified `Distributions`/`ImageFiltering`/`ImageTransformations` recipes), NOT invent a local analog. Key stage primitives (RESEARCH "Don't Hand-Roll"): `MvNormal`/Cholesky of 2×2 Σ for correlation; `imfilter(·, Kernel.gaussian(σ_psf))` for PSF; `warp(ch2, Translation(dx,dy); fillvalue=0.0)` for sub-pixel shift; `rand(rng, Poisson(g·I))/g .+ rand(rng, Normal(0,r))` for noise.

**Contract traps the forward model must respect (RESEARCH §Common Pitfalls):**
1. **Background = small positive, NOT 0.0.** `_exclude_zero` (`src/colocalization.jl:234`, called inside `correlation`) drops 0.0 pairs; large zero regions thin a patch below the 15-px floor → spurious `missing`. Emit positive autofluorescence offset as background.
2. **≥15-px patch floor.** `correlation` sets `missing` when `length(a) <= 15` after `_exclude_zero` (line 235). Trivially met at all swept sizes (smallest 256² → 32×32 = 1024 px/patch) provided intensities stay mostly non-zero.
3. **`warp` fillvalue=0.0** (dropped cleanly), never leave NaN borders (RESEARCH Pitfall 5).

---

### `spike/simulator/calibration.jl` (new)

**Analog:** `spike/00_smoke.jl` (offline-script + correctness-assert shape) for *structure only*; the calibration *method* is **RESEARCH §Induced-μ Calibration** (sweep ρ_true → measure induced μ via real summary → fit monotone ĝ via isotonic/PCHIP → KS/Wasserstein evidence). No in-repo analog for the fit itself. Evidence (ρ_true→μ curve, induced-vs-target overlap, match metric) is appended to `spike/NOTES.md §3`.

---

### `spike/figures/*` (new CairoMakie outputs) — NO in-repo analog

**Analog:** **NONE.** The repo's existing plotting is GLMakie in `src/plot.jl` (interactive, NOT applicable). Use the **RESEARCH §CairoMakie figure recipe** verbatim (`CairoMakie.activate!()`, `Figure`, `heatmap`/`scatterlines`/`hist`, `save(joinpath(@__DIR__,"figures","plausibility.png"), fig)`). D-13: spike-local only; this does NOT resolve the `src/plot.jl` backend question.

## Shared Patterns

### AGPL header (every new `spike/*.jl`)
**Source:** `spike/00_smoke.jl` lines 1-19 (identical block in `spike/test/runtests.jl`, `src/*.jl`).
**Apply to:** `prior.jl`, `forward.jl`, `calibration.jl`, `contract.jl`, `02_simulator_demo.jl`, `test_simulator.jl`. Copy the `#= … =#` block verbatim.

### Explicit rng threading (D-14)
**Source:** evolves `spike/00_smoke.jl:41` (`Random.seed!(2026)`) into an explicit `rng::AbstractRNG` argument.
**Apply to:** `sample_prior(rng, …)`, `simulate_pair(rng, θ; imsize)`, and every stochastic stage in `forward.jl`. Use `Random.Xoshiro(seed)`; Random123 is deferred (Phase 3).

### include() read-only coupling
**Source:** `spike/NOTES.md §4` (lines 142-149) — the Phase-1 proven fallback.
**Apply to:** `contract.jl` (and any test/demo needing the summary). Load StatsBase+Statistics first; never edit `src/`.

### Decoupling proof (carry-forward gate)
**Source:** `spike/NOTES.md §5` (lines 162-180), baseline `f581d95`.
**Apply to:** post-phase verification — re-run `git diff --quiet f581d95 -- Project.toml Manifest.toml src/`. All new deps touch ONLY `spike/Project.toml` + `spike/Manifest.toml`.

### CUDA-absence / NeuralEstimators-pin guard
**Source:** `spike/test/runtests.jl` lines 52-59.
**Apply to:** PRESERVE in `runtests.jl` after adding Images/CairoMakie; extend with a `NeuralEstimators == v0.2.1` assertion (RESEARCH §Resolve-risk gate).

## No Analog Found

| File / Stage | Role | Data Flow | Reason | Use Instead |
|------|------|-----------|--------|-------------|
| `spike/simulator/forward.jl` (7-stage physics) | service | transform | No microscopy forward-simulator exists in repo | RESEARCH §Architecture Patterns 1–4 (Distributions / ImageFiltering / ImageTransformations recipes) |
| `spike/figures/*` (CairoMakie) | artifact | file-I/O | Repo plotting is GLMakie (`src/plot.jl`), interactive, inapplicable; CairoMakie is new | RESEARCH §CairoMakie figure recipe (D-13, spike-local) |
| `spike/simulator/calibration.jl` (ĝ fit) | service | batch | No isotonic/PCHIP induced-prior fit exists in repo | RESEARCH §Induced-μ Calibration (monotone-inverse / quantile matching) |

## Metadata

**Analog search scope:** `spike/` (Phase-1 files), `src/colocalization.jl`, `src/LoadImages.jl`, `src/bayes.jl` (priors via NOTES seed), `spike/NOTES.md`.
**Files scanned:** 7 (00_smoke.jl, test/runtests.jl, Project.toml, NOTES.md, colocalization.jl, LoadImages.jl, + RESEARCH/CONTEXT).
**Pattern extraction date:** 2026-06-26
