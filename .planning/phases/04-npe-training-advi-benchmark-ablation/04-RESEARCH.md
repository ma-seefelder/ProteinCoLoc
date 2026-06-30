# Phase 4: NPE Training + ADVI Benchmark + Ablation - Research

**Researched:** 2026-06-30
**Domain:** Amortized simulation-based inference (NeuralEstimators.jl NPE) + variational baseline (Turing/AdvancedVI ADVI) + benchmarking + summary-statistic ablation, all inside the decoupled `spike/`.
**Confidence:** HIGH for the NeuralEstimators v0.2.1 API and the modern Turing/AdvancedVI `vi` API (both verified against installed source and current docs); MEDIUM for the cross-method target-space alignment and benchmark unit-of-comparison (methodological design, flagged); MEDIUM for the pre-registered tolerances/margins (Claude's discretion, concrete defensible values proposed).

## Summary

This phase trains one `PosteriorEstimator` (MLP summary net → `NormalisingFlow`) on the Phase-3 leak-free standardized loader output to infer the 7-dim simulator θ (ρ_true is row 1), derives Δρ by Monte-Carlo differencing two single-stack passes (D-03), benchmarks it against the existing `colocalization()` ADVI for per-parameter RMSE / interval width and a >100× speedup at comparable accuracy (NPE-01/02/03), and ablates the minimal (128-dim) vs augmented (142-dim) cached summaries under k=5 CV scoring per-parameter RMSE as a gating sufficiency diagnostic (ABL-01/02).

The NeuralEstimators v0.2.1 API needed is fully confirmed from the installed package source: `train(estimator, θ_train, θ_val, Z_train, Z_val; ...)` trains directly on fixed pre-simulated `d×K` / `inputdim×K` column-major Float32 matrices — exactly the shape `load_fold` already returns (`Ztr`/`θtr`). The ADVI baseline runs in its own isolated `spike/baseline/` env (D-01) because Turing cannot co-resolve with the pinned NeuralEstimators 0.2.1 (Phase-1 landmine) and the parent env's Turing is broken; the removed `vi(m, ADVI(n,iter))`/`syms` calls port to the modern `vi(model, q_meanfield_gaussian, max_iter; adtype=...)` → `VIResult` API (Turing v0.45 / AdvancedVI v0.7).

**Primary recommendation:** Train a CPU-only (`use_gpu=false`) `PosteriorEstimator(Chain(Dense(128|142→…→dstar)), NormalisingFlow(7; num_summaries=dstar))` on the loader's per-fold `(Ztr, θtr)` tensors; benchmark in **ρ-space** — push ADVI's `μ_sample`/`μ_control` posteriors through the frozen `ghat` (μ→ρ) so both methods land on the known ρ_true ground truth; pre-register **NPE ρ_true RMSE ≤ 1.2× ADVI** (D-09) and **aug-wins margin δ=0.05 relative + ≥4/5-fold consistency** (D-07); land everything as re-runnable `spike/test/` gates.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**ADVI ground-truth baseline (NPE-02/03)**
- **D-01:** Dedicated isolated baseline env, ported API — NOT the spike env, NOT `src/`. The ADVI baseline runs in its own isolated Julia environment (e.g. `spike/baseline/` with its own `Project.toml`) carrying a modern Turing/AdvancedVI; the hierarchical `@model` is copied/`include()`d read-only and the removed `vi(m, ADVI(int,int))` / `syms` calls are ported to the current API there. It is NOT added to the spike env (Turing co-resolve would cap NeuralEstimators below the pinned 0.2.1 — the Phase-1 landmine) and NOT fixed in `src/` (spike decoupling holds until Phase 7).
- **D-02:** ADVI runs offline → JLD2 artifact the spike consumes. The baseline runs ADVI on the reserved holdout pairs and serializes `{per-parameter posterior summaries, interval widths, wall-clock}` to a JLD2 artifact; the spike's NPE benchmark reads that artifact for the RMSE and >100× comparison. The parent project's own env is also broken (Manifest pins uninstalled Turing; `src` uses the removed API) — the baseline must NOT depend on repairing it.

**NPE targets — ρ_true and Δρ (NPE-01)**
- **D-03:** One NPE over per-stack θ; Δρ by differencing two single-stack passes. A single `PosteriorEstimator` infers the per-stack θ (including ρ_true). The Δρ posterior = Monte-Carlo difference of two independent single-stack posterior passes (sample − control). Reuses the Phase-3 single-stack cache unchanged, stays fully amortized (2 forward passes), and assumes between-stack independence (matches the sample-vs-control experimental design).
- **D-04:** Δρ benchmark = paired reserved-holdout stacks, symmetric on both methods. Sample/control pairs are formed from the reserved holdout; `Δρ_true = ρ_true(sample) − ρ_true(control)`. Both the NPE (difference of two passes) and ADVI (joint on the pair) run on identical pairs → apples-to-apples Δρ RMSE. The per-stack ρ_true benchmark stays per-stack.

**Summary network**
- **D-05:** MLP summary net (DeepSet deferred). An MLP conditioner over the fixed standardized loader output (128-dim D-01 / 142-dim D-02) → NeuralEstimators `NormalisingFlow`. DeepSet (BACK-02) is deferred until MLP sufficiency is measured, not assumed. The D-01 mask rows (65:128) are part of the input vector and are bypassed by loader standardization (Phase-3 D-08). Depth/width/activation = research/discretion.

**Summary ablation (ABL-01/02)**
- **D-06:** Ablate both cached variants under leak-free k=5 CV, zero re-simulation. Train the NPE on the minimal (D-01) and augmented (D-02) summaries — both already cached (Phase-3 D-02) — and compare per-parameter RMSE under the leak-free 5-fold CV the loader enforces (DATA-03).
- **D-07:** ρ_true/Δρ-led, advisory-with-rule gate. The headline decision is on ρ_true (and Δρ) RMSE, with all 7 θ reported. Choose augmented only if it materially beats minimal on ρ_true RMSE by a pre-set margin; otherwise keep minimal for parsimony + OOD detectability (couples to Phase 5). A SBC-pass-but-high-RMSE outcome is treated as an insufficiency signal, not a pass. Research sets the exact margin/tie-break.

**Speedup benchmark protocol (NPE-03)**
- **D-08:** Raw-stack → posterior on identical inputs; summary time counted; training excluded. Time both methods from the same raw `MultiChannelImage` stack to posterior. NPE clock = summary extraction + forward pass (training excluded — the amortized claim); ADVI clock = the full `vi()` run. CPU-only, `BenchmarkTools`, over the ≥20 reserved holdout stacks. Speedup is always reported paired with RMSE, never alone.
- **D-09:** "Comparable RMSE" = pre-registered tolerance. NPE per-parameter RMSE must fall within a pre-set, pre-registered tolerance of ADVI's (e.g. ≤1.2×; research fixes the exact value before running).
- **D-10:** CPU gates; GPU is an optional bonus. All pass/fail criteria (>100×, RMSE) are decided on CPU. GPU timings reported only if a CUDA device is present, skipping gracefully otherwise (Flux CUDA extension, no forced import).
- **D-11:** KernelAbstractions.jl only if non-perturbing, summary/simulator kernels only. KA evaluated for cpu/gpu auto-dispatch of the summary/simulator kernels only (the network's dispatch is already handled by Flux), and adopted only if it does not perturb the pinned spike env / NeuralEstimators 0.2.1. Otherwise CPU/GPU dispatch stays Flux-native.
- **D-12:** Complexity scaling — the >100× as a curve, both methods. Characterize wall-clock vs (a) number of datasets N (NPE ~flat O(1)/dataset post-training, ADVI ~linear in N) and (b) input size (imsize / patch grid), for both NPE inference and ADVI; report empirical scaling exponents.
- **D-13:** CPU core-scaling — reported characterization, not a gate. Sweep thread counts (e.g. 1/2/4/8) reporting NPE (summary+forward) and ADVI wall-clock vs `Threads.nthreads()` on the reserved holdout; identify parallel speedup / Amdahl saturation. The >100× headline is stated at a fixed, reported thread count.

### Claude's Discretion
- MLP depth/width/activation; flow type/coupling/transform stack; optimizer/LR/epochs/early-stopping/batch size.
- Exact comparable-RMSE tolerance (D-09), ablation margin/tie-break (D-07), thread-count sweep (D-13), and the N / imsize ranges for the scaling study (D-12) — all pre-registered before the reported run.
- Baseline env layout/location, the ADVI `iter`/`num_samples` defaults used for the baseline, and the JLD2 artifact schema (D-01/D-02).
- The exact KA evaluation depth (D-11) and how the optional GPU path is wired (D-10).
- Whether other θ parameters join ρ_true/Δρ in the headline table (all 7 are reported regardless).

### Deferred Ideas (OUT OF SCOPE)
- SBC, amortized Bayes factor (`RatioEstimator`), OOD flag — Phase 5. Δρ correctness vs `compute_BayesFactor()` is reproduced in BF-02 (Phase 5), not here.
- DeepSet permutation-invariant summary (BACK-02) — upgrade only if the MLP (D-05) proves insufficient by the ablation.
- RxInfer.jl independent cross-check baseline (BACK-01) — post-Go, never a spike dependency.
- Productionization / user-definable `num_patches` — Phase 7.
- Making GPU a required deliverable / putting CUDA in the reproducibility path — out of scope; CPU is the gating baseline.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| NPE-01 | `PosteriorEstimator` (summary net → NormalisingFlow) trained to infer both ρ_true and Δρ | `PosteriorEstimator(network, NormalisingFlow(7; num_summaries=dstar))` over d=7 θ; ρ_true is row 1 of `θtr` (NamedTuple field order). Δρ by MC-differencing two `sampleposterior` passes (D-03). Both API constructors verified from installed source. |
| NPE-02 | For ≥20 held-out stacks, single-pass posteriors benchmarked vs real `colocalization()` ADVI on 20–30 stacks for RMSE + interval width, CV-reported | `load_holdout` returns the reserved ≥20-stack set; `sampleposterior`/`posteriormean`/`interval`/`assess`+`rmse` give per-parameter RMSE and interval width. ADVI via the isolated baseline env → JLD2 artifact. RMSE in ρ-space via `ghat`. |
| NPE-03 | NPE wall-clock >100× faster than per-dataset ADVI at comparable RMSE (ms vs min) | `BenchmarkTools.@belapsed` on raw-stack→posterior (D-08); ADVI = full `vi()` run; NPE = summary extraction + forward pass. Speedup curve over N and input size (D-12). Paired with RMSE always (D-09: ≤1.2×). |
| ABL-01 | Per-parameter RMSE compared between minimal and augmented summaries as a gating sufficiency diagnostic | `load_fold(dir, f; variant=:min)` vs `variant=:aug`, same arch, k=5 CV, `assess`+`rmse`. Pre-registered aug-wins margin (proposed δ=0.05 relative + ≥4/5 folds). |
| ABL-02 | Chosen summary justified by ablation; OOD-detectability interaction noted | Decision rule outputs `:min` or `:aug`; note couples to Phase-5 OOD (minimal = more parsimonious/orthogonal-discrepancy-detectable). |
</phase_requirements>

## Architectural Responsibility Map

This phase is a pipeline of stages, not multi-tier app code. The "tier" is the owning subsystem/environment.

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Leak-free standardized training tensors | `spike/data/loader.jl` (`Loader`) | — | Sole standardization path (DATA-03/D-07/D-08); already returns `d×K` Float32 |
| NPE architecture + training | NeuralEstimators v0.2.1 + Flux v0.16.10 (spike env) | — | `PosteriorEstimator`/`NormalisingFlow`/`train` (fixed-data form) |
| Posterior inference (ρ_true, all 7 θ) | NeuralEstimators `sampleposterior`/`posteriormean`/`interval` | — | Single forward pass, amortized |
| Δρ posterior | Spike-local MC differencing of two single-stack passes (D-03) | NeuralEstimators sampling | Reuses single-stack cache; assumes between-stack independence |
| ADVI ground-truth baseline | **Isolated `spike/baseline/` env** (modern Turing/AdvancedVI) | read-only `include()` of `src/bayes.jl` `@model` | Turing cannot co-resolve with pinned NeuralEstimators 0.2.1 (Phase-1 landmine); parent env broken |
| Cross-env hand-off | JLD2 artifact (D-02) | — | The ONLY coupling between the broken-Turing baseline and the pinned spike |
| Target-space alignment (μ↔ρ) | `spike/simulator/ghat.jl` `ghat(μ)→ρ` (frozen) | — | Puts ADVI's `μ_sample` on the NPE's known ρ_true axis |
| Wall-clock benchmark | `BenchmarkTools.jl` (NEW spike dep) | `Threads.nthreads()` sweep | D-08/D-12/D-13 |
| Ablation | `Loader` `variant=:min|:aug` + `assess`/`rmse` | — | Zero re-simulation (both variants cached) |
| Re-runnable gates | stdlib `Test` under `spike/test/` | — | Mirrors Phase-3 `test_data_pipeline.jl` pattern |

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| NeuralEstimators.jl | **0.2.1** (pinned) | `PosteriorEstimator`, `NormalisingFlow`, `train`, `sampleposterior`, `posteriormean/median/quantile`, `interval`, `assess`, `rmse` | The locked stack (CLAUDE.md); fixed-data `train` form matches loader output exactly. `[VERIFIED: installed source ~/.julia/packages/NeuralEstimators/gFxuZ]` |
| Flux.jl | **0.16.10** (pinned) | NN backend for the MLP summary net + flow conditioner; `NormalisingFlow` is **Flux-only** | CPU-clean, CUDA as opt-in extension (D-10). `[VERIFIED: spike/Manifest.toml]` |
| Turing.jl | ~**0.45** (baseline env only) | The hierarchical `@model` + modern `vi()` ADVI | The benchmarked ground truth; isolated per D-01. `[CITED: turinglang.org/docs/tutorials/variational-inference]` |
| AdvancedVI.jl | ~**0.7.0** (transitive of Turing) | `q_meanfield_gaussian`, `KLMinRepGradDescent`, the VI engine | Modern replacement for removed `ADVI(n,iter)`. `[CITED: github.com/TuringLang/AdvancedVI.jl]` |
| BenchmarkTools.jl | current (**NEW spike dep**) | `@belapsed`/`@benchmark` wall-clock for the >100× claim | CLAUDE.md names it for NPE-03. **Not yet in spike Manifest** — Wave-0 add + resolve-risk gate. `[VERIFIED: absent — grep spike/Manifest.toml = 0]` |
| JLD2.jl | current (pinned) | Cross-env ADVI artifact (D-02); trained-model + benchmark-result persistence | Pure-Julia, already the cache backend. `[VERIFIED: spike/Manifest.toml]` |
| StatsBase.jl | current (pinned) | `ZScoreTransform` (loader already uses it; reuse `zt` to freeze train-only preprocessing) | Already a dep. `[VERIFIED: spike/Manifest.toml]` |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| DataFrames.jl | transitive (NeuralEstimators dep) | `assess` returns an `Assessment` wrapping a DataFrame of estimates | Reading per-parameter RMSE; available without a direct add |
| ForwardDiff.jl / ReverseDiff.jl | current (baseline env) | AD backend for `vi(...; adtype=AutoForwardDiff())` | ADVI port; ForwardDiff is the safe default for this small model |
| Random123.jl | current (pinned) | Counter-based reproducibility of benchmark/scaling/holdout re-simulation | Already a dep; keep all runs seed-reproducible |
| CUDA.jl | current (**optional, NOT added**) | GPU timings only if a device is present (D-10) | Loaded as a Flux extension; never a forced import; never in the reproducibility path |
| KernelAbstractions.jl | 0.9.41 (transitive only) | (D-11) cpu/gpu kernel auto-dispatch — **evaluation only** | Already transitively present; adding as a DIRECT dep is a resolve-risk — recommend NOT adopting for the spike (see Pitfall 6) |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Online `train(est, sampler, simulator)` | **Fixed-data `train(est, θ_train, θ_val, Z_train, Z_val)`** | Phase-3 already cached the pairs and enforces leak-free folds; the fixed-data form is the correct one (re-simulating online would bypass the leak-free loader and the holdout exclusion). |
| `q_meanfield_gaussian` ADVI | `q_fullrank_gaussian` | Full-rank captures posterior correlations (better interval-width fidelity) at higher cost; mean-field matches the original `ADVI` intent and is faster. Recommend mean-field for the baseline default; note full-rank as a sensitivity check. |
| Benchmark in ρ-space (via `ghat`) | Benchmark in μ-space | ρ_true is the known simulator ground truth per stack; μ-space has no exact ground truth (only the induced mean). ρ-space is the defensible common axis (see Architecture Pattern 2). |

**Installation (spike env — Wave 0):**
```julia
# In spike/ env ONLY. Then RE-FREEZE Manifest and run the resolve-risk gate.
using Pkg; Pkg.activate("spike")
Pkg.add("BenchmarkTools")
# Assert pin held: NeuralEstimators UUID 38f6df31-… still v0.2.1 (runtests.jl gate (d)/(e)).
```

**Installation (baseline env — separate):**
```julia
# In spike/baseline/ — its OWN Project.toml. Turing must NEVER enter the spike env.
using Pkg; Pkg.activate("spike/baseline")
Pkg.add(["Turing", "DataFrames", "JLD2", "ForwardDiff", "StatsBase"])
# Pin + commit spike/baseline/Manifest.toml as a reproducibility artifact.
```

**Version verification:** NeuralEstimators 0.2.1 and Flux 0.16.10 confirmed from `spike/Manifest.toml` and the installed source tree. BenchmarkTools/Turing/AdvancedVI versions resolve fresh on `Pkg.add`; pin them in the respective Manifests after resolve. Turing v0.45 / AdvancedVI v0.7.0 are the current docs versions (June 2026) — confirm exact bounds in `spike/baseline/Manifest.toml` after resolve (the `vi` return type differs across 0.43→0.45; see Pitfall 3).

## Package Legitimacy Audit

> slopcheck is pip-based and does not cover the Julia (General registry) ecosystem; it was unavailable in this session (`SLOPCHECK_MISSING`). Julia packages are verified instead by (a) presence in the already-resolved/pinned `spike/Manifest.toml`, (b) official docs, and (c) being flagship JuliaStats/FluxML packages. No package here is novel or low-trust.

| Package | Registry | Status | Source Repo | Disposition |
|---------|----------|--------|-------------|-------------|
| NeuralEstimators 0.2.1 | Julia General | pinned, installed | github.com/msainsburydale/NeuralEstimators.jl | Approved (locked stack) |
| Flux 0.16.10 | Julia General | pinned, installed | github.com/FluxML/Flux.jl | Approved |
| BenchmarkTools | Julia General | flagship, NEW direct add | github.com/JuliaCI/BenchmarkTools.jl | Approved — pin after resolve |
| Turing (~0.45) | Julia General | flagship (baseline env) | github.com/TuringLang/Turing.jl | Approved — isolated env only |
| AdvancedVI (~0.7) | Julia General | flagship (Turing dep) | github.com/TuringLang/AdvancedVI.jl | Approved |
| ForwardDiff | Julia General | flagship | github.com/JuliaDiff/ForwardDiff.jl | Approved |
| JLD2 / DataFrames / StatsBase / Random123 | Julia General | pinned/flagship | JuliaIO / JuliaData / JuliaStats | Approved |
| CUDA (optional) | Julia General | flagship, NOT added | github.com/JuliaGPU/CUDA.jl | Optional — never in reproducibility path |

**Packages removed due to slopcheck [SLOP] verdict:** none.
**Packages flagged [SUS]:** none. All are established Julia ecosystem packages; the only NEW direct add is `BenchmarkTools` (CLAUDE.md-mandated), which is a JuliaCI flagship.

## Architecture Patterns

### System Architecture Diagram

```
                          .planning + cache dir (Phase-3 outputs)
                                      │
        ┌─────────────────────────────┴──────────────────────────────┐
        │                                                             │
   shard_*.jld2 (CV pool)                                    holdout.jld2 (≥20 reserved)
        │                                                             │
   Loader.load_fold(dir,f;variant=:min|:aug)              Loader.load_holdout(dir)
   → (Ztr[d×Ktr], θtr[7×Ktr],                             → (theta, summary_min,
      Zva[d×Kva], θva[7×Kva], zt)                            summary_aug, global_index)
        │  (Float32, column-major, mask rows bypass z-score)         │
        ▼                                                            │
   PosteriorEstimator(                                              │ θ + global_index + master key
     Chain(Dense(d→…→dstar)),                                       ▼
     NormalisingFlow(7; num_summaries=dstar))          RE-SIMULATE raw MultiChannelImage
        │ train(est, θtr,θva, Ztr,Zva;                  via UNCHANGED simulate_pair (Random123)
        │       epochs, batchsize, use_gpu=false)        │
        ▼                                                ├──────────────┬─────────────────┐
   trained NPE  ──────────────────────────┐             ▼              ▼                 ▼
        │                                  │      ┌─ SPIKE ENV ─┐  ┌── BASELINE ENV ──┐  (timing
        │  sampleposterior(est, Z; N)      │      │ NPE path:   │  │ colocalization() │   harness)
        │  → 7×N draws per stack           │      │ patch_summary│ │ → vi(m,q_mf,iter)│
        ▼                                  │      │ +encode     │  │ → μ_sample,μ_ctrl │
   ρ_true posterior (row 1)                │      │ +zt.transform│ │  → VIResult.q     │
   Δρ = MC diff(sample,control) (D-03)     │      │ +forward    │  │  rand(q,N)        │
   per-param RMSE via assess/rmse          │      │  (D-08 clock)│ │  (D-08 clock)     │
        │                                  │      └──────┬──────┘  └────────┬──────────┘
        │                                  │             │      JLD2 artifact (D-02)
        ▼                                  │             │      {μ summaries, intervals,
   ABLATION :min vs :aug (k=5 CV)          │             │       wall-clock} ───────────┐
   decision rule (D-07) → :min | :aug      │             ▼                              ▼
                                           └──────►  BENCHMARK: speedup=ADVI/NPE,  ghat(μ)→ρ̂_ADVI
                                                     RMSE(ρ̂ vs known ρ_true),  Δρ RMSE,
                                                     scaling curves (N, imsize), thread sweep
```

### Recommended Project Structure
```
spike/
├── npe/                      # NEW — Phase-4 spike-side code
│   ├── architecture.jl       # build_estimator(d, dstar; depth,width,...) → PosteriorEstimator
│   ├── train_npe.jl          # train over folds (fixed-data form), θ-standardization, save model
│   ├── infer.jl              # sampleposterior wrappers; Δρ MC-diff (D-03); ρ-space mapping (ghat)
│   ├── benchmark.jl          # BenchmarkTools harness (D-08), scaling (D-12), thread sweep (D-13)
│   └── ablation.jl           # :min vs :aug k=5 CV, per-param RMSE, decision rule (D-07)
├── baseline/                 # NEW — ISOLATED env (D-01), own Project.toml/Manifest.toml
│   ├── Project.toml          # Turing, DataFrames, JLD2, ForwardDiff, StatsBase
│   ├── model.jl              # read-only include of src/bayes.jl @model OR a copied lift
│   ├── run_advi.jl           # modern vi() port; re-simulate holdout; write JLD2 artifact (D-02)
│   └── advi_artifact.jld2    # the cross-env hand-off (committed or regenerable fixture)
└── test/
    └── test_npe.jl           # NEW — SC1..SC5 re-runnable gates, included by runtests.jl
```

### Pattern 1: Train the NPE on the fixed leak-free cache
**What:** Use the fixed-data `train` method on the loader's column-major Float32 tensors. One dataset = one column; single replicate (the 128/142 summary vector IS the data for that stack), exactly like the smoke test's `inputdim × K` form.
**When to use:** Always — re-simulating online would bypass the leak-free loader and holdout exclusion.
**Example:**
```julia
# Source: VERIFIED installed NeuralEstimators v0.2.1 src/train.jl + Estimators/PosteriorEstimator.jl
using NeuralEstimators, Flux
include(joinpath(@__DIR__, "..", "data", "loader.jl"))   # Loader (re-exported)

const D      = 7        # θ dimension (ρ_true is ROW 1; NamedTuple field order)
const DSTAR  = 32       # learned-summary dim (discretion; examples use ~3d → 21–64)

fold = load_fold(CACHE_DIR, 1; K=5, master_seed=MASTER_SEED, variant=:min)  # :min → d=128
d_in = size(fold.Ztr, 1)                                                    # 128 (or 142 :aug)

network = Chain(Dense(d_in, 128, gelu), Dense(128, 128, gelu), Dense(128, DSTAR))
q       = NormalisingFlow(D; num_summaries = DSTAR)        # Flux-only; d=7, dstar=DSTAR
est     = PosteriorEstimator(network, q)                   # q passed as an INSTANCE, positionally

# Fixed parameters AND fixed data — the correct form for a pre-simulated cache.
est = train(est, fold.θtr, fold.θva, fold.Ztr, fold.Zva;
            epochs = 200, batchsize = 64, use_gpu = false,   # use_gpu DEFAULTS true → MUST set false
            optimiser = Flux.Optimisers.AdamW(5e-4, (0.9,0.999), 1e-4),  # weight decay vs overfit
            stopping_epochs = 10, verbose = false)
```
Key facts (all `[VERIFIED: installed source]`):
- `train(estimator, θ_train, θ_val, Z_train, Z_val; ...)` — fixed-params-and-data method.
- `θ_train`/`θ_val` are `d×K`; `Z_train`/`Z_val` are `inputdim×K` (matches `θtr`/`Ztr`).
- Keyword defaults: `epochs=100`, `stopping_epochs=5`, `batchsize=32`, `optimiser=Adam(5e-4)`, `lr_schedule=CosAnneal` (auto), `use_gpu=true`, `adtype=AutoZygote()`, `savepath=tempdir()`, `verbose=true`. Docstring explicitly: *"When the training data or parameters are fixed, one may wish to use regularisation to help prevent overfitting."*
- `PosteriorEstimator(summary_network, q::ApproximateDistribution)` or `(summary_network, num_parameters; num_summaries, q=nothing)`.
- `NormalisingFlow(d, dstar; num_coupling_layers=6, use_act_norm=true)` — note positional `(d, dstar)`; the smoke test used the `(d; num_summaries=dstar)` keyword convenience form (both work).

### Pattern 2: Cross-method target-space alignment (the key methodological decision)
**What:** The NPE infers the **simulator** θ (ρ_true, the known ground-truth coloc knob). The Turing/ADVI model infers the **induced** per-condition means `μ_sample`/`μ_control`. They live in different spaces. Benchmark in **ρ-space**: push ADVI's μ-posterior through the frozen `ghat` (μ→ρ) so both methods have a per-parameter RMSE against the **known ρ_true** of each holdout stack.
**When to use:** For every NPE-vs-ADVI RMSE and the Δρ comparison (D-04 "identical pairs, apples-to-apples").
**Example:**
```julia
# Source: spike/simulator/prior.jl (ρ_true = ghat(μ*)) + spike/simulator/ghat.jl (ghat(μ::Real)→Float64)
# ADVI posterior samples of μ_sample → implied ρ̂ posterior:
ρ̂_advi_samples = ghat.(μ_sample_samples)         # monotone map, applied per draw
ρ̂_advi         = mean(ρ̂_advi_samples)
# Δρ (both methods on identical pairs, D-04):
Δρ_true = θ_true_sample[1] - θ_true_control[1]    # known
Δρ_npe  = mean(ρ_npe_sample_draws .- ρ_npe_control_draws)        # MC diff (D-03)
Δρ_advi = mean(ghat.(μ_s_draws) .- ghat.(μ_c_draws))
```
**Rationale:** ρ_true is exactly known per holdout stack (it is the simulator input); μ-space has no exact ground truth. `ghat` is already frozen and validated (SIM-02). This is MEDIUM-confidence as a *design choice* — flag for the planner/discuss as the central benchmark-axis decision (`[ASSUMED]`, see Assumptions Log A1).

### Pattern 3: The unit of comparison — ADVI is inherently paired, NPE is per-stack
**What:** `colocalization(img, control, channels, num_patches; iter, posterior_samples)` jointly models a **sample stack + control stack** and emits `μ_sample` and `μ_control` from one `vi()` run. The NPE's atomic unit is a **single** stack (one forward pass). So:
- **Δρ benchmark (D-04):** atomic unit = one (sample,control) PAIR. ADVI = one `vi()` call; NPE = two `sampleposterior` passes. This is the clean apples-to-apples timing for the >100× headline (full `vi()` vs 2 forward passes).
- **Per-stack ρ_true benchmark (NPE-02):** read both `μ_sample` and `μ_control` out of the paired `vi()` run, map each through `ghat`, and score every stack (sample and control alike) against its known ρ_true. NPE scores each stack from its single pass.
**When to use:** Forming the holdout into pairs (D-04); reporting per-pair wall-clock for the speedup, per-stack RMSE for accuracy.
**Note:** `colocalization()` also runs `sample(m, Prior(), posterior_samples)` (a 100k prior chain) which is **not** part of the ADVI cost — the baseline should time **only** the `vi()` call (D-08: "ADVI clock = the full `vi()` run"), and may skip the prior chain entirely (INFER-3 in spike 003 README notes it is wasteful/dead).

### Pattern 4: ADVI port to the modern Turing/AdvancedVI API (baseline env)
**What:** The removed pre-0.35 calls port to Turing v0.45 / AdvancedVI v0.7.
**Example:**
```julia
# Source: CITED turinglang.org/docs/tutorials/variational-inference (v0.45) + AdvancedVI.jl README (v0.7)
#         cross-checked with .planning/spikes/003-bayes-model-modularity/README.md (INFER-1)
using Turing
m = coloc_model(ctrl_data, sample_data)        # the lifted/included @model (ports VERBATIM — body
                                               # captures no enclosing var, verified by spike 003)
# REMOVED: q = vi(m, ADVI(num_latent, iter))   # ADVI(n,iter) gone; vi(m,alg) 2-arg gone
result = vi(m, q_meanfield_gaussian, ITER;     # family passed by reference; ITER = max iterations
            adtype = AutoForwardDiff(),
            show_progress = false)
q = result.q                                   # VIResult struct: .q (approx) + .info (diagnostics)
draws = rand(q, POSTERIOR_SAMPLES)             # d_latent × N
```
Critical facts:
- `ADVI(n, iter)` → `const ADVI = KLMinRepGradDescent`; first positional arg is an **AD backend**, not two ints. The default algorithm is `KLMinRepGradProxDescent`. `[CITED: AdvancedVI.jl]`
- `vi(model, family, max_iter; adtype, algorithm, show_progress, ...)` returns a **`VIResult`** with `.q`, `.info`, `.state`, `.ldf`. `[CITED: turinglang.org]`
- `DynamicPPL.syms(VarInfo(m))` is **removed** — replace with `keys(DynamicPPL.VarInfo(m))` or extract the 8 globals (`μ_sample`, `μ_control`, …) by known declaration order.
- **Verify constrained vs unconstrained space** (Pitfall 3): the model has `Truncated` priors; confirm `rand(q, N)` yields **constrained** μ ∈ [−1,1] (apply the model's inverse bijector if not). This is the single biggest correctness risk in the port.
- The `@model` itself ports **verbatim** (spike 003 verified it captures no enclosing variable). Only the inference call + name extraction change.

### Anti-Patterns to Avoid
- **Adding Turing to the spike env.** Co-resolve caps NeuralEstimators below the pinned 0.2.1 (Phase-1 GLMakie-class regression). Keep it in `spike/baseline/` (D-01). The `runtests.jl` resolve-risk gate (assert UUID `38f6df31-…` == v0.2.1) must extend to cover the BenchmarkTools add.
- **Benchmarking in μ-space or comparing ρ_true (NPE) directly against μ_sample (ADVI).** Different quantities — always map through `ghat` (Pattern 2).
- **Timing `train` as part of the speedup.** D-08: training is excluded (the amortized claim). NPE clock = summary extraction + forward pass only.
- **Leaving `use_gpu` at its default.** It defaults to **true** — every `train`/inference call in the spike must pass `use_gpu=false` (or `device=cpu_device()`). CPU is the hard gate (D-10).
- **Feeding un-standardized or Float64 θ/Z to the flow.** Use the loader's Float32 output; standardize θ leak-free if adopted (Pitfall 5).
- **Re-z-scoring the 0/1 mask rows (65:128).** The loader already bypasses them (Phase-3 D-08); never standardize them downstream.
- **Running ADVI on the cached *summaries*.** `colocalization()` needs raw `MultiChannelImage`s — re-simulate the holdout from stored θ + index under the recorded Random123 key (Open Question 1).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Normalising-flow posterior | A custom coupling-flow stack | `NormalisingFlow(d, dstar; num_coupling_layers, use_act_norm)` | Built-in, Flux-backed, ActNorm+permutation handled; the locked stack |
| Posterior sampling / point estimates | Manual flow inversion + sampling loop | `sampleposterior`, `posteriormean`, `posteriormedian`, `posteriorquantile` | Verified API; handles single-dataset vs Vector-of-datasets dispatch |
| Per-parameter RMSE / bias | Hand-rolled error loops over draws | `assess(est, θ_test, Z_test)` → `rmse(assessment)` / `bias` / `risk` | Per-parameter by default (`average_over_parameters` kwarg); returns a tidy DataFrame |
| Credible intervals / interval width | Manual quantile bookkeeping | `interval(θ; probs=[0.05,0.95])` or `posteriorquantile(est, Z, probs)` | d×2 matrix; named rows; interval width = col2−col1 |
| Variational inference | A custom ELBO optimiser | `vi(model, q_meanfield_gaussian, iter; adtype=...)` | The whole point of the ADVI baseline; modern AdvancedVI |
| Wall-clock timing | `@time` / manual `time_ns()` deltas | `BenchmarkTools.@belapsed` / `@benchmark` | Warmup, sampling, GC control, median+IQR — required to substantiate >100× |
| Leak-free standardization | Re-standardizing in the NPE code | `load_fold(...).zt` (frozen `ZScoreTransform`) | The loader is the SOLE standardization path (DATA-03); reuse `zt` |
| μ→ρ mapping | A new induced-mean inverse | `ghat(μ)` (frozen, SIM-02-validated) | Already calibrated; do not recompute |

**Key insight:** Almost everything this phase needs is a one-liner against NeuralEstimators v0.2.1 or AdvancedVI; the genuinely novel work is the *wiring* (loader→train, ρ-space alignment, paired benchmark, isolated baseline env, JLD2 hand-off) and the *pre-registered decision rules* — not algorithms.

## Common Pitfalls

### Pitfall 1: `use_gpu` defaults to true
**What goes wrong:** `train`/inference silently tries the GPU path; on a machine without a working CUDA stack it errors or (worse) forces a CUDA load into the pinned spike env.
**Why it happens:** NeuralEstimators `train(...; use_gpu::Bool = true)` default.
**How to avoid:** Pass `use_gpu=false` (or `device=cpu_device()`) on every `train`, `sampleposterior`, `posteriormean`, `assess` call. The GPU path (D-10) is a separate opt-in branch guarded by a CUDA-device check.
**Warning signs:** CUDA appears in `Base.loaded_modules`; the `runtests.jl` CPU-only testset (b)/(c) goes red.

### Pitfall 2: NeuralEstimators 0.2.1 downgraded by the BenchmarkTools add
**What goes wrong:** `Pkg.add("BenchmarkTools")` co-resolves the pinned NeuralEstimators below 0.2.1 (the Phase-1 GLMakie landmine class).
**Why it happens:** Pre-1.0 packages have tight compat bounds.
**How to avoid:** After the add, re-assert the v0.2.1 pin (extend `runtests.jl` gate (e)); re-freeze and commit `spike/Manifest.toml`. If it downgrades, pin NeuralEstimators explicitly or add BenchmarkTools only to a benchmark-scoped env.
**Warning signs:** `Pkg.dependencies()[UUID("38f6df31-…")].version != v"0.2.1"`.

### Pitfall 3: ADVI samples come back in unconstrained space
**What goes wrong:** `rand(result.q, N)` yields μ values outside [−1,1] (the model's `Truncated` support), corrupting RMSE and `ghat` (which expects μ in the calibrated range).
**Why it happens:** `vi` optimizes in unconstrained space (`unconstrained=true` by default); whether `result.q` auto-applies the inverse bijector differs across Turing 0.43→0.45.
**How to avoid:** After the port, sanity-check `extrema(μ_sample_draws) ⊂ [−1,1]`; if not, apply the model's `Bijectors` inverse transform before reading marginals. Add an assertion to the baseline run.
**Warning signs:** μ draws with |value| > 1; `ghat` clamping or NaNs; nonsensical RMSE.

### Pitfall 4: ADVI runs on cached summaries instead of raw images
**What goes wrong:** The baseline tries to feed `holdout.jld2` summaries into `colocalization()`, which needs raw `MultiChannelImage`s.
**Why it happens:** The cache stores summaries (and θ, index, imsize), not raw images.
**How to avoid:** Re-simulate each holdout stack from its stored θ + `global_index` under the recorded Random123 master key via the UNCHANGED `simulate_pair`/`build_mci`, reproducing byte-identical images, then run ADVI. The NPE timing path (D-08) must start from the SAME re-simulated raw stack (summary extraction included).
**Warning signs:** Type errors in `_prepare_data`; mismatched stacks between methods.

### Pitfall 5: Heterogeneous θ scales hurt flow training
**What goes wrong:** The flow's standard-Gaussian base struggles when θ components span wildly different ranges (ρ∈[−1,1], spillover∈[0,0.2], autofl∈[0,0.1], label_eff∈[0.6,1], shift∈[−1,1], noise∈[0,1]); posterior mass leaks outside physical ranges.
**Why it happens:** The loader returns **raw** θ (Float32, un-standardized); only the summaries `Z` are standardized.
**How to avoid:** Standardize θ **leak-free** (fit a `ZScoreTransform` on `θtr` only, transform `θtr`/`θva`), train in standardized space, and **un-standardize posterior draws** before RMSE / `ghat`. ActNorm helps but does not replace this. Document the θ-transform as part of the frozen preprocessing (parallels `zt` for Z; Phase 5 will need it frozen too).
**Warning signs:** ρ̂ draws outside [−1,1]; poor RMSE on small-range params; unstable training loss.

### Pitfall 6: KernelAbstractions adopted as a direct dep
**What goes wrong:** Adding KA to `[deps]` to "auto-dispatch" the summary/simulator kernels risks perturbing the pinned spike env for marginal benefit on CPU-cheap kernels.
**Why it happens:** KA is already transitively present (v0.9.41), so it *looks* free — but a direct dep changes the resolve graph.
**How to avoid:** Per D-11, keep CPU/GPU dispatch **Flux-native** (the network already dispatches via Flux). Do NOT add KA directly unless a resolve check proves it non-perturbing AND a profiled bottleneck justifies it. The summary/simulator kernels are small (8×8 grid); KA is unlikely to pay off in the spike.
**Warning signs:** Manifest churn; NeuralEstimators pin moves after `Pkg.add("KernelAbstractions")`.

### Pitfall 7: `Threads.nthreads()` is fixed per process
**What goes wrong:** The D-13 thread sweep tries to change thread count inside one Julia session.
**Why it happens:** `nthreads()` is set at launch (`julia -t N`) and immutable at runtime.
**How to avoid:** Run the sweep as separate processes (`julia -t 1`, `-t 2`, `-t 4`, `-t 8`), each writing its timings to JLD2/CSV; aggregate afterward. State the >100× headline at one fixed, reported thread count.
**Warning signs:** Identical timings across "different" thread counts.

## Runtime State Inventory

This phase is greenfield spike code plus a new isolated env; it renames nothing in existing systems. The only stateful concerns are environment isolation and the cross-env artifact.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | Phase-3 cache (`shard_*.jld2`, `holdout.jld2`, `meta.jld2`) is READ-ONLY input; content-hash-named dir. New writes: trained NPE model(s), `advi_artifact.jld2`, benchmark/ablation result files — all NEW under `spike/npe/`,`spike/baseline/`. | None to existing data; add new artifacts only |
| Live service config | None — no external services. | None — verified (CPU-only, offline) |
| OS-registered state | None — no scheduled tasks/daemons. The D-13 thread sweep spawns transient `julia -t N` processes only. | None |
| Secrets/env vars | None. | None |
| Build artifacts / installed packages | spike env gains BenchmarkTools (re-freeze `spike/Manifest.toml`); NEW `spike/baseline/` env with its own `Project.toml`/`Manifest.toml`. Parent env's broken Turing pins are NOT touched (D-02). | Re-freeze + commit both Manifests; resolve-risk gate |

**Nothing found in categories** Live service / OS-registered / Secrets — verified: the phase is CPU-only, offline, simulation-driven.

## Code Examples

### Inference + per-parameter RMSE + interval width on the holdout
```julia
# Source: VERIFIED installed NeuralEstimators v0.2.1 src/inference.jl + assess.jl
ho = load_holdout(CACHE_DIR)                          # reserved ≥20-stack set
Z  = Float32.(standardize_holdout(ho.summary_min, fold.zt))   # apply FROZEN train-fit zt
# Single-pass posteriors (Vector{Matrix}, one 7×N per stack):
draws = sampleposterior(est, Z; N = 2000)             # use_gpu handled by est/device
ρ̂      = [posteriormean(d)[1] for d in draws]          # row 1 = ρ_true (un-standardize first!)
ints  = [interval(d; probs=[0.05,0.95]) for d in draws]
width = [I[1,2] - I[1,1] for I in ints]               # ρ_true 90% interval width

# Tidy per-parameter RMSE via assess (θ_test = ho.theta as 7×K):
assessment = assess(est, Float32.(ho.theta), Z)       # Assessment(DataFrame)
per_param_rmse = rmse(assessment)                      # per-parameter by default
```

### BenchmarkTools harness (D-08), per pair
```julia
# Source: julia-development skill + CLAUDE.md (BenchmarkTools for the >100× claim)
using BenchmarkTools
# NPE clock = summary extraction + forward pass, from the RE-SIMULATED raw stack:
t_npe = @belapsed begin
    Zs = $(standardize)(encode_d01(patch_summary($sample_mci)), $fold.zt)
    Zc = $(standardize)(encode_d01(patch_summary($control_mci)), $fold.zt)
    sampleposterior($est, Zs; N=$N); sampleposterior($est, Zc; N=$N)
end
# ADVI clock = full vi() run on the pair (baseline env, separate process → JLD2):
t_advi = advi_artifact["wall_clock"][pair_id]
speedup = t_advi / t_npe                               # ALWAYS reported with RMSE (D-08/D-09)
```

### Proposed JLD2 ADVI artifact schema (D-02, Claude's discretion)
```julia
# spike/baseline/run_advi.jl writes:
jldsave("advi_artifact.jld2";
  schema_version = 1,
  meta = (turing_version = "0.45.x", adtype = "AutoForwardDiff", iter = ITER,
          posterior_samples = N, master_seed = MASTER_SEED, bayes_git_ref = "<sha>",
          generated = "2026-…"),
  pairs = pair_ids,                                    # align to holdout global_index
  mu_sample = mu_s,  mu_control = mu_c,                # per-pair posterior SAMPLE vectors
  rho_sample = ghat.(mu_s), rho_control = ghat.(mu_c), # ρ-space (precomputed convenience)
  rho_s_interval = ints_s, rho_c_interval = ints_c,    # 90% interval bounds
  wall_clock = t_vi_seconds)                           # per-pair vi() time only
```
Store raw μ **sample vectors** (not just means) so the spike can recompute Δρ via MC differencing and any interval the benchmark needs.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `vi(m, ADVI(num_latent, iter))` | `vi(model, q_meanfield_gaussian, max_iter; adtype=…)` → `VIResult` | Turing ≥0.35 (ADVI≡KLMinRepGradDescent) | The baseline port; `result.q` + `rand(q,N)` |
| `DynamicPPL.syms(VarInfo(m))` | `keys(DynamicPPL.VarInfo(m))` / explicit global selection | DynamicPPL 0.40/0.41 | Parameter-name extraction in the baseline |
| Per-dataset ADVI minutes | Amortized NPE single forward pass (ms) | This phase's thesis | The >100× claim (NPE-03) |
| `train(est, sampler, simulator)` online | `train(est, θ_train, θ_val, Z_train, Z_val)` on the leak-free cache | Phase-3 cache exists | Leak-free, holdout-excluded training |

**Deprecated/outdated:**
- `ADVI(n, iter)` two-int constructor — removed; `ADVI` is now an alias for `KLMinRepGradDescent`.
- 2-arg `vi(m, alg)` — removed; needs `vi(model, family, max_iter; …)`.
- Parent project's pinned Turing 0.42.8 / DynamicPPL 0.39.13 — **not installed**; baseline resolves fresh (D-02), never repairs the parent.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Benchmark in **ρ-space** via `ghat(μ_sample)` is the right common axis for NPE-vs-ADVI RMSE | Pattern 2 | If μ-space is preferred, all RMSE numbers and the comparable-RMSE gate change basis; flag to discuss-phase |
| A2 | Pre-registered **comparable-RMSE tolerance = 1.2×** ADVI ρ_true RMSE (D-09) | Validation Architecture | A looser/tighter bar changes the pass/fail verdict; must be locked BEFORE the run |
| A3 | Pre-registered **aug-wins margin: δ=0.05 relative + improvement in ≥4/5 folds** (D-07) | Validation Architecture | A different margin flips the summary choice (:min vs :aug), which couples to Phase-5 OOD |
| A4 | **MLP arch**: 2 hidden layers × 128, gelu, `dstar=32`, flow `num_coupling_layers=6`, AdamW(5e-4, wd=1e-4), epochs≈200, batch 64, `stopping_epochs=10` | Pattern 1 | Discretion (D-05/Claude's); may need tuning if RMSE poor — not a correctness risk |
| A5 | ADVI baseline defaults: `q_meanfield_gaussian`, `ITER≈1000`, `AutoForwardDiff`, `posterior_samples≈10000` | Pattern 4 | Affects baseline accuracy/runtime; full-rank or more iters may tighten intervals |
| A6 | θ should be **leak-free standardized** before flow training, draws un-standardized after | Pitfall 5 | If skipped and the flow handles raw θ fine via ActNorm, harmless; if needed and skipped, poor small-range-param RMSE |
| A7 | Holdout raw images are **reproducible** from stored θ + `global_index` + master key | Open Question 1 | If not byte-reproducible, must persist raw holdout images (extra cache/storage) |
| A8 | Δρ atomic benchmark unit = a (sample,control) **pair**; >100× headline uses full `vi()` vs 2 NPE passes | Pattern 3 | If a per-single-stack ADVI is expected instead, the timing protocol changes |
| A9 | BenchmarkTools/Turing/AdvancedVI resolve to current versions without perturbing the pin | Standard Stack | Resolve-risk; mitigated by the gate, but a downgrade would block |

## Open Questions (RESOLVED)

1. **Holdout raw-image reproducibility for the ADVI baseline.**
   - What we know: `generate.jl` uses Random123 keyed-per-`global_index` seeding; `simulate_pair`/`build_mci` are deterministic under a fixed rng; the cache stores θ, `global_index`, and `imsize`.
   - What's unclear: whether the exact per-stack rng key is recoverable from `(global_index, master_seed)` alone to reproduce byte-identical holdout images.
   - Recommendation: in Wave 0, prove round-trip reproducibility (re-simulate a holdout stack, assert summary == cached summary). If it fails, persist raw holdout images alongside `holdout.jld2`.
   - **RESOLVED:** Planned as 04-01 Task 2 — a Wave-0 round-trip reproducibility gate re-simulates a holdout stack and asserts `summary == cached summary`, with raw-image persistence as the documented fallback if it fails.

2. **ADVI per-stack ρ without a control.**
   - What we know: the Turing model is inherently sample+control; `μ_sample` is the per-condition mean of interest.
   - What's unclear: the cleanest way to score per-stack ρ_true for a stack used as `control` vs `sample`.
   - Recommendation: form pairs (D-04), run one `vi()` per pair, read BOTH `μ_sample` and `μ_control`, map each via `ghat`, and score every stack against its own known ρ_true (Pattern 3). No single-stack ADVI needed.
   - **RESOLVED:** Planned as 04-04 (Pattern 3) — pairs are formed per D-04, one `vi()` runs per pair, both `μ_sample`/`μ_control` are read and mapped via `ghat`, and each stack is scored against its own known ρ_true. No single-stack ADVI is required.

3. **Constrained vs unconstrained ADVI draws** (Pitfall 3) — resolve empirically in the baseline port with an assertion on μ-support.
   - **RESOLVED:** Planned as 04-04 Task 1 — the baseline port carries an explicit bijector/μ-support assertion (`μ ∈ [−1,1]`) so the draw space is verified empirically rather than assumed.

4. **Interval-width comparability.** NPE 90% interval (flow quantiles) vs ADVI 90% interval (Gaussian VI) measure different posterior shapes; report both honestly, note mean-field VI tends to under-disperse (a known VI artifact) — consider `q_fullrank_gaussian` as a sensitivity check (A5).
   - **RESOLVED:** Planned as 04-05 Task 1 — both methods' interval widths are reported honestly side-by-side with the mean-field under-dispersion caveat noted; `q_fullrank_gaussian` remains an available sensitivity check (A5).

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Julia | everything | ✓ | 1.12.6 | — |
| NeuralEstimators | NPE-01/02 | ✓ | 0.2.1 (pinned) | — |
| Flux | summary net + flow | ✓ | 0.16.10 (pinned) | — |
| JLD2 | artifact + cache | ✓ | pinned | — |
| StatsBase | standardization | ✓ | pinned | — |
| KernelAbstractions | (D-11 eval only) | ✓ (transitive) | 0.9.41 | Flux-native dispatch (recommended) |
| BenchmarkTools | NPE-03 timing | ✗ | — | **Add to spike env (Wave 0)** + resolve-risk gate |
| Turing / AdvancedVI | ADVI baseline | ✗ | — | **Create `spike/baseline/` env** (D-01) |
| ForwardDiff | ADVI adtype | ✗ | — | Add to baseline env |
| CUDA | GPU bonus (D-10) | ✗ (assumed) | — | Skip gracefully; CPU is the gate |

**Missing dependencies with no fallback:** none block the CPU-gated deliverables.
**Missing dependencies with fallback:** BenchmarkTools (add), Turing stack (isolated baseline env), CUDA (optional, graceful skip).

## Validation Architecture

> nyquist_validation is enabled (config.json `workflow.nyquist_validation: true`).

### Test Framework
| Property | Value |
|----------|-------|
| Framework | stdlib `Test` (+ `BenchmarkTools` for timings) — mirrors `spike/test/test_data_pipeline.jl` |
| Config file | none — `spike/test/runtests.jl` includes child testset files |
| Quick run command | `julia --project=spike spike/test/test_npe.jl` (a single testset during dev) |
| Full suite command | `julia --project=spike spike/test/runtests.jl` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| NPE-01 | Trained `PosteriorEstimator` returns a 7×N posterior in one pass for ≥20 holdout stacks; Δρ = MC-diff of two passes | integration | `julia --project=spike spike/test/test_npe.jl` (SC1 testset) | ❌ Wave 0 |
| NPE-02 | Per-parameter RMSE + interval width vs ADVI artifact, CV-reported | integration | SC2 testset (reads `advi_artifact.jld2`) | ❌ Wave 0 |
| NPE-03 | `t_advi/t_npe > 100` at comparable RMSE, fixed thread count | benchmark | SC3 testset (`@belapsed`) | ❌ Wave 0 |
| NPE-03 | Scaling curves over N and imsize; empirical exponents (NPE ~flat, ADVI ~linear) | characterization | SC3-scaling testset (D-12) | ❌ Wave 0 |
| ABL-01 | `:min` vs `:aug` per-parameter RMSE under k=5 CV | integration | SC4 testset | ❌ Wave 0 |
| ABL-02 | Decision rule emits `:min`/`:aug` per the pre-registered margin; OOD note recorded | unit | SC5 testset | ❌ Wave 0 |

### Pre-registered constants (declared BEFORE any reported run — mirror Phase-3 `const` fixtures)
```julia
const NPE_RMSE_TOLERANCE   = 1.2     # D-09 (A2): NPE ρ_true RMSE ≤ 1.2× ADVI ρ_true RMSE
const SPEEDUP_GATE         = 100.0   # NPE-03: t_advi/t_npe > 100 at the reported thread count
const ABL_REL_MARGIN       = 0.05    # D-07 (A3): aug wins iff RMSE_aug ≤ (1−δ)·RMSE_min on ρ_true …
const ABL_FOLD_CONSISTENCY = 4       # … AND aug improves in ≥4 of 5 folds (else keep :min, parsimony)
const BENCH_THREADS        = 1       # D-13: headline thread count (sweep {1,2,4,8} reported separately)
const NPE_MASTER_SEED      = 0xC0FFEE  # Random123 reproducibility
```

### Sampling Rate
- **Per task commit:** `julia --project=spike spike/test/test_npe.jl` (the touched SC testset; trained model loaded from cache, not retrained).
- **Per wave merge:** `julia --project=spike spike/test/runtests.jl` (full suite incl. Phase 1–3 gates + resolve-risk).
- **Phase gate:** full suite green; the `advi_artifact.jld2` regenerated-or-verified; all five SC testsets green at the pre-registered constants.

### Wave 0 Gaps
- [ ] `spike/test/test_npe.jl` — SC1..SC5 testsets as named `@test_skip` placeholders (mirror `test_data_pipeline.jl`), included by `runtests.jl`
- [ ] `spike/Project.toml` — add `BenchmarkTools`; re-freeze `spike/Manifest.toml`; extend resolve-risk gate to re-assert NeuralEstimators v0.2.1 with BenchmarkTools present
- [ ] `spike/baseline/Project.toml` + `Manifest.toml` — new isolated env (Turing, DataFrames, JLD2, ForwardDiff, StatsBase), pinned and committed
- [ ] Holdout raw-image reproducibility check (Open Question 1) — Wave-0 round-trip assertion
- [ ] Pre-registered constants committed before any reported run (A2/A3)

## Security Domain

> `security_enforcement` is not set in config.json. This phase is an offline, CPU-only, simulation-driven research spike with no network, auth, untrusted input, secrets, or persistence of user data. Standard application-security (ASVS) categories do not meaningfully apply.

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V5 Input Validation | partial | `simulate_pair` ArgumentErrors propagate (never swallowed) — Phase-2 convention; JLD2 artifact schema-versioned + integrity-checked (Phase-3 cache pattern) |
| V6 Cryptography | no | — |
| V2/V3/V4 (auth/session/access) | no | No users, sessions, or access control in a local spike |

**Reproducibility-as-integrity (the relevant analog):** the content-hash cache guard (Phase-3), the pinned Manifests (both envs), and the Random123 seeding are the integrity controls here. The cross-env JLD2 artifact (D-02) should carry a schema version + meta (Turing version, seed, git ref) and be integrity-checked on read, mirroring `cache.jl`'s reopen-check pattern.

## Sources

### Primary (HIGH confidence)
- `~/.julia/packages/NeuralEstimators/gFxuZ/src/{train.jl, Estimators/PosteriorEstimator.jl, ApproximateDistributions/NormalisingFlow.jl, inference.jl, assess.jl}` — installed v0.2.1 source; exact signatures for `train` (fixed-data form), `PosteriorEstimator`, `NormalisingFlow`, `sampleposterior`, `posteriormean/median/quantile`, `interval`, `assess`, `rmse`/`bias`/`risk`. `[VERIFIED]`
- `spike/Manifest.toml` — NeuralEstimators 0.2.1, Flux 0.16.10, KernelAbstractions 0.9.41 (transitive), BenchmarkTools absent. `[VERIFIED]`
- `spike/00_smoke.jl`, `spike/test/runtests.jl` — confirmed NPE API call shape (`inputdim×K` Z, `use_gpu=false`), CPU-only gate, resolve-risk gate pattern. `[VERIFIED]`
- `spike/data/{loader.jl, encode.jl, cache.jl, generate.jl}`, `spike/contract.jl`, `spike/simulator/{prior.jl, ghat.jl, forward.jl}` — loader output shape `(Ztr 128|142×K, θtr 7×K, …, zt)`, θ field order (ρ_true row 1), `ghat(μ)→ρ`, holdout/re-simulation path. `[VERIFIED]`
- `src/bayes.jl` — `colocalization()` (line 234), removed `vi(m, ADVI(num_latent, iter))` (line 325), `CoLocResult` (line 75), `_prepare_data` (line 164), `compute_BayesFactor` (line 109), the `@model` (lines 267–303). `[VERIFIED]`
- `.planning/spikes/003-bayes-model-modularity/README.md` — INFER-1 ADVI port analysis (the 0.43→0.45 `vi` return-type caveat), MODEL-1 verbatim-port verification, INFER-3 prior-chain waste. `[VERIFIED: repo]`

### Secondary (MEDIUM confidence)
- turinglang.org/docs/tutorials/variational-inference (v0.45) — `q_meanfield_gaussian`, `vi(model, family, max_iter; adtype, algorithm, show_progress)` → `VIResult.q`, `rand(q,n)`. `[CITED]`
- github.com/TuringLang/AdvancedVI.jl (v0.7.0, 2026-06-07) — `KLMinRepGradDescent`, `q_fullrank_gaussian`/`q_meanfield_gaussian`, ADVI alias. `[CITED]`
- msainsburydale.github.io/NeuralEstimators.jl/dev (navigation/index) — confirms API page structure; detailed pages fetched via installed source instead (dev pages 404'd on fetch). `[CITED]`

### Tertiary (LOW confidence)
- WebSearch on flow θ-standardization (general SBI practice) — supports Pitfall 5 directionally; the leak-free θ-transform recommendation is applied judgment. `[ASSUMED]`

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — NeuralEstimators API verified from installed v0.2.1 source; modern Turing `vi` API verified from current official docs.
- Architecture: MEDIUM-HIGH — wiring follows verified loader/simulator contracts; the ρ-space alignment (A1) and paired-unit benchmark (A8) are sound but are design choices to confirm.
- Pitfalls: HIGH — each traces to verified source behavior (use_gpu default, resolve-risk gate, removed API, raw-vs-summary, θ scales).
- Pre-registered tolerances: MEDIUM — concrete defensible values proposed (1.2×, δ=0.05/4-of-5); locked by Claude's discretion before the run.

**Research date:** 2026-06-30
**Valid until:** ~2026-07-30 for the pinned spike stack (Manifest frozen); the baseline Turing/AdvancedVI API is fast-moving — re-verify `vi` signature/return type against the resolved `spike/baseline/Manifest.toml` (7-day currency on that env).
