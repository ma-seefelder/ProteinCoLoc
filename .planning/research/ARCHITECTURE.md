# Architecture Research

**Domain:** Amortized simulation-based inference (SBI) for fluorescence-microscopy colocalization — a decoupled research spike (`spike/`) followed by a conditional productionization into `src/`
**Researched:** 2026-06-26
**Confidence:** MEDIUM-HIGH (HIGH on component structure / data flow / decoupling — anchored in existing code + two concrete reference engines in `BayesInteractomics`; MEDIUM on exact `NeuralEstimators.jl` API surface, which is intentionally gated behind the AP1 smoke test)

---

## Executive Orientation

This is **not** a CRUD/web architecture. It is a **scientific-computing pipeline** with two strictly separated halves:

1. **The spike (`spike/`)** — a horizontal-layers research pipeline, each layer a standalone `.jl` script with one job, chained by a final `demo.jl`. Build bottom-up; every layer is independently runnable and its output is a checkpoint (plot, `.jld2`, or printed verdict). The decoupling constraint (no edits to `src/`) is an **architectural invariant**, not a preference.
2. **Productionization (conditional on Go)** — fold the validated estimators into `src/` behind a clean function boundary that **coexists** with the existing Turing/ADVI `colocalization()` path rather than replacing it.

The single most important architectural fact: **the existing summary statistic is the contract between the simulator and the neural network.** `_prepare_data` → `patch` → `correlation` already produce the exact `num_patches × num_patches` Pearson-correlation grid that the NPE/NRE consume as input. The simulator's only job is to emit a `MultiChannelImage` that those *unmodified* functions can ingest. This is what makes the spike comparable to the ADVI baseline and what keeps `src/` untouched.

---

## Standard Architecture

### System Overview — Spike (horizontal layers, bottom-up build)

```
┌──────────────────────────────────────────────────────────────────────┐
│  L6  DEMO / REPORT          spike/demo.jl  (fixed Random123 seed)      │
│      chains every layer end-to-end → Go/No-Go memo                     │
├──────────────────────────────────────────────────────────────────────┤
│  L5  VALIDATION & DIAGNOSTICS                                          │
│   ┌────────────────┐  ┌────────────────┐  ┌────────────────────────┐  │
│   │ sbc.jl         │  │ evidence_net.jl│  │ ood.jl                 │  │
│   │ SBC ranks /    │  │ amortized      │  │ misspecification flag  │  │
│   │ coverage / ECE │  │ log-BF (NRE)   │  │ (summary density)      │  │
│   └───────┬────────┘  └───────┬────────┘  └──────────┬─────────────┘  │
│           │                   │                       │                │
├───────────┴───────────────────┴───────────────────────┴───────────────┤
│  L4  NEURAL INFERENCE        spike/npe_model.jl  (+ nre_model.jl)      │
│   summary network (MLP / DeepSet over patches)                        │
│      → NormalizingFlow conditioner  =  PosteriorEstimator (NPE)        │
│      → ratio classifier            =  RatioEstimator    (NRE)          │
│   artifacts: npe_model.jld2, nre_model.jld2                            │
├──────────────────────────────────────────────────────────────────────┤
│  L3  TRAINING-DATA PIPELINE  spike/training_data.jl                    │
│   θ~π → simulate_pair → [_prepare_data/correlation] → 8×8 vector       │
│      → standardize → JLD2 cache  (train 50–200k + held-out 5k)         │
├──────────────────────────────────────────────────────────────────────┤
│  L2  SUMMARY STATISTIC  ===  REUSED FROM src/, UNMODIFIED  ===         │
│   patch() · correlation() · _prepare_data()   ← the simulator/NN      │
│                                                  contract              │
├──────────────────────────────────────────────────────────────────────┤
│  L1  FORWARD SIMULATOR       spike/simulator.jl                        │
│   sample_prior() → simulate_pair(θ)::MultiChannelImage                 │
│   π(θ) ranges mirror the Turing @model priors (μ/ν/σ/τ)               │
├──────────────────────────────────────────────────────────────────────┤
│  L0  ENVIRONMENT + SMOKE     spike/Project.toml · 00_smoke.jl          │
│   isolated env; NeuralEstimators+Flux toy NPE train+sample = green     │
└──────────────────────────────────────────────────────────────────────┘
                              ▲   READ-ONLY   ▲
        ┌─────────────────────┴───────────────┴────────────────────┐
        │  src/ (UNTOUCHED during spike) — provides L2 functions    │
        │  colocalization.jl · bayes.jl · LoadImages.jl             │
        │  also: real colocalization()/ADVI = the L4 benchmark      │
        └───────────────────────────────────────────────────────────┘
```

### Component Responsibilities — Spike

| Component | Responsibility (owns) | Does NOT own | Typical Implementation |
|-----------|----------------------|--------------|------------------------|
| `simulator.jl` | Generative physics `g(θ)→MultiChannelImage`; `sample_prior()`; π(θ) consistency with Turing priors | Summary stats, NN, validation | Distributions.jl + ImageFiltering.jl `imfilter`; returns existing `MultiChannelImage` struct |
| **`src/` summary fns (reused)** | Map image → fixed `num_patches²` correlation vector | Anything new; **never edited** | Called, not rewritten: `patch`, `correlation`, `_prepare_data` |
| `training_data.jl` | θ→sim→summary→standardize→JLD2; train/held-out split; seeding | Simulation internals, NN | Loop + `jldsave(compress=true)`; Random123 `Philox`/`Threefry` |
| `npe_model.jl` | Summary net + flow conditioner; `train`/`sampleposterior`; ADVI benchmark | BF, SBC, OOD | `PosteriorEstimator(DeepSet(ψ,ϕ), q=NormalizingFlow)` |
| `nre_model.jl` / `evidence_net.jl` | Amortized log-BF via likelihood-ratio classifier; validate vs `compute_BayesFactor()` | Posterior, OOD | `RatioEstimator`; l-POP / BCE loss |
| `sbc.jl` | Rank histograms, KS/χ² uniformity, coverage curve, ECE/MCE traffic-light | Training, simulation | Reuse `_bin_calibration`/`CalibrationResult` template |
| `ood.jl` | Misspecification flag from summary-space density / posterior-predictive mismatch | Posterior accuracy | Mahalanobis or flow log-density on summary vector |
| `demo.jl` | Deterministic end-to-end chaining; emits metrics table + memo | Any computation logic (only orchestrates) | Single script, fixed seed |

### System Overview — Productionization (post-Go integration surface)

```
┌──────────────────────────────────────────────────────────────────────┐
│  PUBLIC API  (ProteinCoLoc exported)                                  │
│   colocalization(img, control, channels, num_patches; ...)  ← EXISTING│
│   colocalization_amortized(img, control, channels, num_patches;       │
│                            estimator=load_default(), backend=:npe)    │
├──────────────────────────────────────────────────────────────────────┤
│  INFERENCE BACKENDS (dispatch boundary — shared input contract)       │
│   ┌──────────────┐  ┌──────────────┐  ┌───────────────────────────┐   │
│   │ :turing/advi │  │ :amortized   │  │ :rxinfer (eval / future)  │   │
│   │ (existing)   │  │ (NPE/NRE)    │  │                           │   │
│   └──────┬───────┘  └──────┬───────┘  └────────────┬──────────────┘   │
│          └─────────────────┴──── shared ───────────┘                  │
│                  _prepare_data → summary vector                       │
├──────────────────────────────────────────────────────────────────────┤
│  SHIPPED ARTIFACTS  (src/inference/amortized/)                        │
│   trained estimator(s) per num_patches grid · loader · calibration    │
└──────────────────────────────────────────────────────────────────────┘
```

### Component Responsibilities — Productionization

| Component | Responsibility | Boundary note |
|-----------|---------------|---------------|
| `colocalization_amortized()` | New public entry; same input contract (`img, control, channels, num_patches`) as `colocalization()`; returns a `CoLocResult`-compatible object | Must return the SAME result type family so `compute_BayesFactor` keeps working unchanged |
| Estimator registry/loader | Resolve a trained estimator for the requested `num_patches`; lazy-load `.jld2`; error if no estimator for that grid | The `num_patches` flexibility lives HERE, not in the NN |
| Backend dispatch | Select `:amortized` / `:turing` (/ `:rxinfer`); the two coexist | No shared mutable state; pure functions over the summary vector |
| Promoted simulator/SBC | Migrated from spike with provenance; used for re-training & shipped calibration proof | Only promoted after Go; spike code is the reference, not a live dependency |

---

## Recommended Project Structure

### Spike (decoupled — own environment)

```
ProteinCoLoc/
├── src/                       # ← UNTOUCHED during spike (read-only contract source)
├── Project.toml               # ← main env, NOT modified by spike
└── spike/                     # entire spike lives here
    ├── Project.toml           # OWN env: NeuralEstimators, Flux, Distributions,
    │                          #   JLD2, Images, ImageFiltering, KernelDensity,
    │                          #   QuadGK, StatsBase, Random123, HypothesisTests
    ├── Manifest.toml          # own, isolated
    ├── NOTES.md               # copied signatures + Turing prior ranges (Day-1)
    ├── 00_smoke.jl            # L0  — API/Flux gate (toy 1-param Gaussian NPE)
    ├── simulator.jl           # L1  — sample_prior() + simulate_pair(θ)
    ├── training_data.jl       # L3  — generator + JLD2 cache + loader
    ├── npe_model.jl           # L4  — PosteriorEstimator (NPE)
    ├── nre_model.jl           # L4  — RatioEstimator (NRE) ┐ may merge into
    ├── evidence_net.jl        # L5  — amortized log-BF      ┘ one file
    ├── sbc.jl                 # L5  — calibration/coverage
    ├── ood.jl                 # L5  — misspecification flag
    ├── demo.jl                # L6  — deterministic end-to-end + memo
    ├── artifacts/             # generated (gitignore large blobs)
    │   ├── training_data.jld2
    │   ├── npe_model.jld2
    │   └── nre_model.jld2
    └── figures/               # plausibility plots, SBC histograms, reliability
```

**To reuse `src/` functions without coupling:** in `spike/` scripts, `include("../src/colocalization.jl")` and `include("../src/LoadImages.jl")` (read-only), OR `Pkg.develop(path="..")` the main package into the spike env and call qualified names. **Recommendation: `Pkg.develop` the parent package** — it gives clean `ProteinCoLoc.correlation(...)` calls, guarantees you exercise the *real* shipped functions, and makes the eventual promotion a no-op import change. `include`-ing raw files risks double-definition and load-order fragility.

### Productionization (post-Go, inside `src/`)

```
src/
├── ProteinCoLoc.jl            # add exports: colocalization_amortized
├── colocalization.jl          # EXISTING — untouched logic; summary fns now shared
├── bayes.jl                   # EXISTING — compute_BayesFactor reused as-is
├── LoadImages.jl              # EXISTING
└── inference/                 # NEW subtree
    ├── amortized.jl           # colocalization_amortized() entry + dispatch
    ├── simulator.jl           # promoted from spike (provenance preserved)
    ├── estimator_io.jl        # registry/loader keyed by num_patches
    ├── calibration.jl         # promoted SBC + traffic-light
    └── artifacts/             # shipped trained estimators (per grid)
```

### Structure Rationale

- **`spike/` flat script layout:** research code optimizes for *runnability and inspection*, not module hygiene. One script = one layer = one checkpoint. A reader (or future-you, post-pause) can run any layer top-to-bottom and see a plot or printed verdict. Premature modularization here is an anti-pattern (see below).
- **Own `Project.toml`:** the decoupling invariant. A heavy, possibly-Flux/CUDA-fragile dependency tree must never leak into the published package's manifest until Go. The two manuscript pipelines stay reproducible.
- **`src/inference/` subtree in production:** isolates the new dependency-heavy code path so the *core* package (image loading, correlation, ADVI) can still be loaded without pulling NeuralEstimators if a user only wants the classic path (consider an extension/weakdep — see Integration Points).

---

## Architectural Patterns

### Pattern 1: Summary statistic as the simulator↔network contract

**What:** The simulator's output type and the network's input type are decoupled by a *fixed, pre-existing* transform: `MultiChannelImage → _prepare_data → standardized 8×8 (=64-dim) vector`. Neither the simulator nor the network knows about the other; both know only the contract.

**When to use:** Always, here. It is the reason `src/` stays untouched and the reason NPE output is comparable to the ADVI baseline (same input).

**Trade-offs:** (+) Maximum decoupling, baseline comparability, zero `src/` edits. (−) The summary is *fixed* and possibly not sufficient — mitigated by the planned ablation (patch-corr only vs + Manders/moments). The fixed `num_patches=8×8` is a hard constraint on the trained net; flexibility is deferred to production via the estimator registry.

```julia
# spike/training_data.jl — the contract in code (illustrative)
img = simulate_pair(θ)                                  # L1 output: MultiChannelImage
stack = MultiChannelImageStack([img], "sim")
s = _prepare_data(stack, [1,2], 8)                      # L2 reused, UNMODIFIED
z = standardize(reshape_to_fixed(s, 8*8))               # L3: fixed 64-dim summary
```

### Pattern 2: DeepSet permutation-invariant aggregation over patches/replicates

**What:** NeuralEstimators' core architecture is the **DeepSet**: an inner network ψ applied per-element, a permutation-invariant pool (mean/sum), then an outer network ϕ. Patches (and independent images in a stack) are *exchangeable* — order carries no information — so a DeepSet summary network is the principled choice over a flat MLP on a fixed vector.

**When to use:** Use a plain MLP first (simplest, matches the fixed 64-vector, gets you to M2 fastest). Upgrade to DeepSet when (a) the patch-corr-only summary proves insufficient, or (b) you want a single net to handle a *variable number of images per stack* — which is also exactly what unlocks the production `num_patches`/variable-stack flexibility.

**Trade-offs:** (+) Permutation invariance baked in; handles variable-size inputs; the natural path to the productionization flexibility. (−) More moving parts than an MLP; harder to debug if the smoke layer hasn't proven the API. Start MLP, keep DeepSet as the documented upgrade.

```julia
# L4 — conceptual NeuralEstimators construction (verify exact API at AP1 smoke)
ψ = Chain(Dense(1=>32, relu), Dense(32=>32, relu))      # per-patch
ϕ = Chain(Dense(32=>64, relu), Dense(64=>d_θ))          # after pooling
summary_net = DeepSet(ψ, ϕ)                             # permutation-invariant
npe = PosteriorEstimator(summary_net; q = NormalizingFlow(...))
train(npe, sample_prior, simulate_summary)              # sampler + simulator
θ_draws = sampleposterior(npe, z_obs)                   # ms-scale forward pass
```

### Pattern 3: Horizontal layers with checkpoint artifacts (research-pipeline shape)

**What:** Each layer writes a durable artifact (`.jld2`, figure, or printed metric) consumed by the next. `demo.jl` re-chains them deterministically. This is the BayesInteractomics simulation-engine shape: compute → `jldsave(compress=true, …, param_hash, timestamp)` → version-checked `load`.

**When to use:** Any expensive, resumable, pausable research pipeline (explicitly a project constraint — "lowest priority, pausable").

**Trade-offs:** (+) Resumable across pauses; each layer independently inspectable; expensive training cached. (−) Cache-invalidation discipline required — adopt the reference engine's `cache_version` + `param_hash` guard so a changed simulator/prior auto-invalidates stale data.

```julia
# L3 cache guard — mirror BayesInteractomics save/load_simulation_cache
jldsave(path; compress=true, cache_version=SPIKE_CACHE_VERSION,
        param_hash=hash(prior_spec), timestamp=string(now()),
        theta=Θ, summaries=Z, held_out=HO)
# loader returns nothing on version/hash mismatch → forces regeneration
```

### Pattern 4: Coexisting inference backends behind one input contract (production)

**What:** `colocalization()` (Turing/ADVI) and `colocalization_amortized()` (NPE/NRE) are sibling functions sharing the `_prepare_data` summary and returning the same `CoLocResult` family, selected by argument/dispatch. A future `:rxinfer` backend slots in identically.

**When to use:** Productionization only. Keeps the published, peer-reviewed ADVI path as the trusted default while shipping amortized inference as an opt-in accelerator.

**Trade-offs:** (+) No regression risk to the published method; users choose speed vs familiarity; clean A/B for the paper. (−) Two code paths to maintain; the amortized path carries trained-artifact + calibration baggage. Mitigate by returning a result type `compute_BayesFactor` already understands.

---

## Data Flow

### Training-time flow (build/offline)

```
sample_prior()            π(θ) mirrors Turing μ/ν/σ/τ ranges
     ↓  θ
simulate_pair(θ)          bivariate corr densities → Bernoulli thinning →
     ↓  MultiChannelImage   PSF convolution → 2×2 spillover → autofluor →
     │                      sub-pixel shift → Poisson/Gaussian noise
[_prepare_data / correlation]   ← REUSED src/, UNMODIFIED
     ↓  num_patches² vector
standardize + cache       JLD2 (train 50–200k, held-out 5k), Random123 seed
     ↓
train(estimator, …)       summary net → flow (NPE) / ratio (NRE)
     ↓
npe_model.jld2 / nre_model.jld2
```

### Inference-time flow (amortized, ms-scale)

```
observed MultiChannelImage(s)
     ↓ [_prepare_data]            ← same contract as training
standardized summary vector z
     ↓
┌─ sampleposterior(npe, z) → posterior draws (ρ_true, Δρ)   → SBC-checked
└─ logratio(nre, z)        → amortized log Bayes factor      → vs compute_BayesFactor()
     ↓
ood(z): summary-density score → in-dist? quiet : raise flag
```

### Validation flow (SBC / coverage)

```
for m in 1:~2000:
    θ* ~ π  →  simulate  →  z  →  L posterior draws  →  rank(θ*) among draws
ranks → histogram per parameter → KS/χ² uniformity
     → coverage curve (nominal vs empirical)
     → _bin_calibration → CalibrationResult(ECE, MCE) → traffic light
```

### Key Data Flows (named)

1. **The contract flow:** `MultiChannelImage → _prepare_data → fixed-dim summary` — identical at train and inference time; the linchpin of decoupling and baseline comparability.
2. **The benchmark flow:** the SAME held-out images go through (a) `npe` forward pass and (b) real `src/colocalization()` ADVI; metrics = RMSE, interval width, **wall-clock speedup**.
3. **The calibration loop:** simulate→infer→rank, repeated ~2000×, feeding the traffic-light verdict — the project's headline differentiator.

---

## Build Order & Dependencies

Strict bottom-up; each layer gated by the one below. Maps directly to AP1–AP7 / M0–M5.

| Order | Layer / File | Depends on | Gate (must pass before next) | AP / M |
|-------|--------------|-----------|------------------------------|--------|
| 1 | `Project.toml` + `00_smoke.jl` | nothing (isolated env) | Toy NPE trains+samples; Flux backend green | AP1 / M0 |
| 2 | `simulator.jl` | smoke (proves stack), `src` types (read) | Monotone ρ_true→corr; spillover/shift visible in plots | AP2 / M1 |
| 3 | `training_data.jl` | simulator + reused `_prepare_data` | Cache writes/loads; held-out split; deterministic under seed | AP3 |
| 4 | `npe_model.jl` | training cache | NPE RMSE ≈ ADVI at **>100× speedup** | AP4 / M2 |
| 5a | `sbc.jl` | trained NPE + simulator | Ranks uniform (KS p>0.05) *or* miscalibration cleanly reported | AP5 / M3 |
| 5b | `nre_model.jl` / `evidence_net.jl` | training cache (+NPE) | log-BF matches `compute_BayesFactor()` in well-spec regime | AP6 / M4 |
| 5c | `ood.jl` | summaries (+ trained nets) | Fires on misspecified images, quiet in-distribution | AP6 / M4 |
| 6 | `demo.jl` + memo | all above | Reproducible from fixed seed; Go/No-Go decision | AP7 / M5 |
| 7 | **(conditional)** `src/inference/` | Go decision + promoted spike code | Coexists with ADVI; `num_patches` user-definable | Production |

**Dependency notes:**
- L0 smoke is a **hard gate before any investment** — it de-risks the single biggest unknown (NeuralEstimators flow convergence + Flux-on-Windows). Do not build the simulator before it is green.
- 5a/5b/5c are siblings — all depend on L4 artifacts and the simulator, but not on each other; can be built in any order or parallel.
- Production (7) depends ONLY on a Go memo; the spike remains the reference implementation and is *promoted* (copied with provenance), never imported live.

---

## Scaling Considerations

"Scale" here = simulation budget, training-set size, and inference throughput — not concurrent users.

| Scale | Architecture adjustments |
|-------|--------------------------|
| Spike (2D, 50k pairs) | CPU-only; single-threaded sim fine; in-memory + one JLD2 cache. The portable baseline. |
| Spike accuracy push (→200k) | Multithread `simulate_pair` loop (`Threads.@threads`); optional CUDA.jl *only* to accelerate `train` — must degrade gracefully to CPU (Windows constraint). |
| Production (variable `num_patches`, many stacks) | Estimator registry per grid; DeepSet to handle variable stack size; pre-trained artifacts shipped so users pay only the ms inference cost. |
| Full build-out (3D / hierarchy — out of scope) | InvertibleNetworks.jl (GPU flows), possibly BayesFlow via PythonCall; re-architect summary net for hierarchy. Explicitly deferred. |

### Scaling Priorities

1. **First bottleneck:** simulation throughput at 200k pairs. Fix: thread the generator loop; cache aggressively with version/hash guard. Training is secondary (minutes on CPU for 2D).
2. **Second bottleneck:** fixed `num_patches` blocks reuse across image sizes. Fix (production): registry keyed by grid + DeepSet for variable input — the explicit production deliverable.

---

## Anti-Patterns

### Anti-Pattern 1: Editing or "lightly refactoring" `src/` to make the simulator fit

**What people do:** Tweak `correlation`/`_prepare_data` signatures "just a little" so simulated images flow more conveniently.
**Why it's wrong:** Violates the decoupling invariant; contaminates the published manuscript pipelines; destroys baseline comparability (the NPE would no longer consume the *same* statistic as ADVI).
**Do this instead:** Make the simulator emit a real `MultiChannelImage` that the *unmodified* functions accept. If a helper is genuinely needed, write it in `spike/`. `Pkg.develop` the parent package and call qualified, read-only.

### Anti-Pattern 2: Premature modularization of the spike

**What people do:** Build a polished `SpikeSBI` module with abstract types, interfaces, and a package skeleton before the principle is proven.
**Why it's wrong:** The spike's job is a Go/No-Go answer in ~4 weeks; module hygiene is wasted effort if the answer is No, and obstructs the layer-by-layer inspectability that makes research debuggable.
**Do this instead:** Flat scripts, one per layer, each runnable top-to-bottom with a visible checkpoint. Modularize *only* at promotion into `src/inference/`.

### Anti-Pattern 3: Prior drift between simulator and Turing model

**What people do:** Pick "reasonable" prior ranges for `sample_prior()` independently of the existing `@model`.
**Why it's wrong:** The whole benchmark (NPE vs ADVI) and the BF validation assume a shared generative prior. Divergent priors make every comparison meaningless and can silently bias SBC.
**Do this instead:** Copy the exact truncations/ranges from the Turing `@model` (μ ~ Truncated(Cauchy(0,0.3),−1,1); ν ~ Exponential; σ,τ ~ Truncated Cauchy on (0,1)) into `spike/NOTES.md` on Day 1; have `sample_prior()` reference those constants. Treat any change as a cache-invalidating event.

### Anti-Pattern 4: Building before the smoke test is green

**What people do:** Write the simulator and data pipeline first because they feel productive, deferring the NeuralEstimators integration.
**Why it's wrong:** If the flow API is unstable or Flux misbehaves on Windows, weeks of simulator work sit on a broken foundation. The plan explicitly makes AP1 the gate.
**Do this instead:** 20-line toy NPE (`00_smoke.jl`) MUST pass before L1. If it fails, trigger the documented fallback (NormalizingFlows.jl, then BayesFlow via PythonCall) *before* investing.

### Anti-Pattern 5: Hiding the simulator-realism gap

**What people do:** Tune the simulator until in-distribution metrics look great and quietly hope real images match.
**Why it's wrong:** The sim-to-real gap is the central scientific risk; concealing it produces an overconfident, non-generalizing tool.
**Do this instead:** Make the gap *measurable* — the OOD/misspecification flag (L5) is designed precisely to surface it. A miscalibration or OOD finding is a publishable result, not a failure.

---

## Integration Points

### External Tools / Libraries

| Tool | Integration pattern | Notes / gotchas |
|------|--------------------|-----------------|
| NeuralEstimators.jl | Primary SBI engine: `PosteriorEstimator` (NPE) + `RatioEstimator` (NRE) + `train`/`sampleposterior`/`assess` | Verify exact constructor/flow API at AP1 smoke — training data is staleness-prone. DeepSet for exchangeable inputs is its core idiom. |
| Flux.jl | NN backend for summary net + flow | Windows/CUDA friction is the named risk → spike is CPU-only; GPU optional accelerator that must degrade gracefully. |
| ImageFiltering.jl / Images.jl | PSF convolution (`imfilter`), sub-pixel shift, noise | Already a project dependency; reuse to keep image semantics consistent with `load_tiff`. |
| Distributions.jl | π(θ) sampling, bivariate correlated densities, Poisson/Gaussian noise | Ranges must mirror Turing priors. |
| JLD2.jl | Cache for training data + trained estimators | Adopt BayesInteractomics `jldsave(compress=true, …, cache_version, param_hash, timestamp)` + version-checked loader. |
| HypothesisTests.jl | KS/χ² uniformity of SBC ranks | Feeds traffic-light. |
| Random123.jl | Counter-based RNG for full reproducibility | Seed once in `demo.jl`; thread-safe streams for parallel simulation. |
| NormalizingFlows.jl / BayesFlow (PythonCall) | Documented fallbacks | Kept off the core path; only if NeuralEstimators flow won't converge. |

### Reusable reference templates (from sibling `BayesInteractomics`)

| Template | Path | What to copy |
|----------|------|--------------|
| `_bin_calibration` + `CalibrationResult` | `BayesInteractomics/src/diagnostics/calibration.jl`, `types.jl` | ECE/MCE binning + struct + `show` traffic-light → drop into `sbc.jl` nearly verbatim (replace "empirically positive" proxy with SBC coverage indicator). |
| Parametric simulation engine | `BayesInteractomics/src/simulation/simulation.jl` | Sweep-over-scenario grid + replicates + confidence bands + `save/load_*_cache` with `cache_version`/`param_hash` guard → the shape for `training_data.jl` and the SBC sweep. |

### Internal Boundaries

| Boundary | Communication | Notes |
|----------|---------------|-------|
| `spike/` ↔ `src/` | **Read-only**, via `Pkg.develop` qualified calls | The decoupling invariant. No writes, no edits, ever, during spike. |
| simulator ↔ network | The fixed summary vector (the contract) | Neither side knows the other's internals. |
| L4 ↔ L5 (validation) | `.jld2` artifacts (trained estimators) + simulator | SBC/NRE/OOD all consume the same trained nets + sampler; mutually independent. |
| (production) `colocalization` ↔ `colocalization_amortized` | Shared `_prepare_data` input + shared `CoLocResult` output type | Coexistence, not replacement; pure functions, no shared mutable state. `compute_BayesFactor` consumes either unchanged. |
| (production) core pkg ↔ `inference/amortized` | Consider a package extension / weakdep | Lets users load the classic ADVI path without pulling the NeuralEstimators/Flux tree. |

---

## Open Questions for the Roadmap

- **MLP vs DeepSet for L4 first cut:** recommendation is MLP-first (fastest to M2), DeepSet as the documented upgrade and the production variable-size enabler. Roadmap should flag L4 as the layer most likely to need a second iteration.
- **NRE/NPE file split:** `nre_model.jl` + `evidence_net.jl` may collapse into one file; low risk, decide during AP6.
- **Production weakdep/extension:** whether the amortized path ships as a Julia package extension (keeps core lightweight) is a productionization-phase decision, not a spike concern.
- **RxInfer backend:** purely a coexisting-backend question; architecturally it slots into the same dispatch boundary as `:amortized`. No bearing on spike structure.

## Sources

- Existing code (HIGH): `src/colocalization.jl`, `src/bayes.jl`, `src/LoadImages.jl` — read directly; summary-statistic contract, Turing prior ranges, `MultiChannelImage`/`CoLocResult` types.
- Project planning (HIGH): `.planning/PROJECT.md`, `.planning/10_plan_amortizedcoloc.md` — AP/M structure, constraints, decoupling invariant.
- Reference engines (HIGH): `BayesInteractomics/src/diagnostics/calibration.jl` + `types.jl` (`_bin_calibration`/`CalibrationResult`), `BayesInteractomics/src/simulation/simulation.jl` (sweep + JLD2 cache version/hash guard) — read directly.
- NeuralEstimators.jl (MEDIUM): https://msainsburydale.github.io/NeuralEstimators.jl/dev/ and repo README — confirms `PosteriorEstimator`/`RatioEstimator`/`PointEstimator`, `train(estimator, sampler, simulator)`, `sampleposterior`, `assess`, DeepSet idiom and normalizing-flow approximate distributions. Exact constructor signatures intentionally re-verified at the AP1 smoke gate (training-data staleness).

---
*Architecture research for: amortized-SBI colocalization (decoupled spike → conditional productionization)*
*Researched: 2026-06-26*
