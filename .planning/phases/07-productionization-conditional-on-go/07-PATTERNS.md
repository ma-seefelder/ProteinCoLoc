# Phase 7: Productionization (conditional on Go) - Pattern Map

**Mapped:** 2026-07-03
**Files analyzed:** 13 new/modified src+test targets (+ Project.toml surgery)
**Analogs found:** 13 / 13 (every new src/ file has a direct spike/ or existing src/ analog)

> **Read first (planner):** Phase 7 is a **copy-and-adapt** promotion, not a green-field build. Almost
> every new `src/` file is a lightly-generalized copy of a proven `spike/` file. The core adaptation
> deltas repeated across files are: (1) **hardcoded `8`/`128`/`64` → grid `G` / `2G²` / `G²`**, (2)
> **`use_gpu=false` forced everywhere → `use_gpu` a plumbed keyword** (train GPU, gate/infer CPU),
> (3) **whole-object `jldsave` → `Flux.state`+`loadmodel!`** for shipped artifacts, (4) **flat spike
> top-level functions → module-scoped `src/` API behind the D-02 type hierarchy**, and (5) **the single
> `VAL_MASTER_SEED`/`consts.jl` → per-grid fresh `PROD_SEED[G]`/`gate_consts_G.jl`**. Copy the spike
> body; apply these five deltas.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `src/results.jl` | model (types) | transform | `src/bayes.jl::CoLocResult` (L75-82) | role-match (struct→hierarchy) |
| `src/registry.jl` | store/config | CRUD (lookup+register) | `spike/npe/train_npe.jl` save/load (L152-201) + `spike/data/cache.jl` atomic idiom | role-match |
| `src/amortized/infer.jl` | service | request-response (single-pass read) | `spike/npe/infer.jl` (whole file) | exact |
| `src/amortized/train_npe.jl` + architecture | service | batch (training) | `spike/npe/train_npe.jl` + `spike/npe/architecture.jl` | exact |
| `src/amortized/bf.jl` | service | request-response | `spike/validation/bf.jl` (L82-100) | exact |
| `src/amortized/train_ratio.jl` | service | batch | `spike/validation/train_ratio.jl` | exact |
| `src/amortized/ood.jl` | service | request-response + batch | `spike/validation/ood.jl` (whole file) | exact |
| `patch_summary`/`encode_d01` grid-parametrization | utility | transform | `spike/contract.jl:85-90` + `spike/data/encode.jl:71-75` | exact (add `G` arg) |
| `colocalization_amortized()` public entry | controller | request-response | `spike/npe/infer.jl` surface + `src/bayes.jl::colocalization()` (L234-242 signature) | role-match |
| `ext/ProteinCoLocTuringExt.jl` | provider (weakdep ext) | request-response | existing `src/bayes.jl` body (L234-336) + `src/ProteinCoLoc.jl` imports | role-match (move, not rewrite) |
| `Project.toml` (deps+weakdeps+ext+version) | config | — | current root `Project.toml` + `spike/Project.toml` (proven resolve set) | role-match |
| `test/gate/run_gate.jl` + harness | test | batch | `spike/validation/{harness,sbc}.jl` | exact |
| `test/gate/gate_consts_G.jl` (per grid) | config (test) | — | `spike/validation/consts.jl` | exact (fresh values) |

---

## Pattern Assignments

### `src/results.jl` (model/types, transform) — the D-02 governing deliverable

**Analog:** `src/bayes.jl::CoLocResult` (L75-82) — becomes the internal `AdviColocResult` subtype.

**Existing struct to demote** (`src/bayes.jl:75-82`):
```julia
struct CoLocResult
    img::MultiChannelImageStack
    control::MultiChannelImageStack
    channels::Vector{Int64}
    num_patches::Int64
    posterior::DataFrame
    advi_result
end
```
Rename to `AdviColocResult <: AbstractColocResult`, keep fields verbatim (it is the Turing path's result), move it beside the Turing extension so it is NOT exported. The `Δρ` accessor maps onto the existing `μ_sample`/`μ_control` columns the current `compute_BayesFactor` already uses (`src/bayes.jl:110`):
```julia
# existing field access to preserve in the AdviColocResult accessor
Δρ_post = posterior.posterior.μ_sample .- posterior.posterior.μ_control   # src/bayes.jl:110
```

**Target hierarchy** (from RESEARCH §Pattern 1, D-02 — this is the sketch to implement):
```julia
abstract type AbstractColocResult end
delta_rho(r::AbstractColocResult)       = _iface_error(r, :delta_rho)
bayes_factor(r::AbstractColocResult)    = _iface_error(r, :bayes_factor)
is_ood(r::AbstractColocResult)          = _iface_error(r, :is_ood)
posterior_draws(r::AbstractColocResult) = _iface_error(r, :posterior_draws)

struct AmortizedColocResult <: AbstractColocResult
    grid            :: Int
    posterior       :: Matrix{Float64}   # 7×N physical-θ draws (row 1 = ρ_true)
    delta_rho_draws :: Vector{Float64}
    log_bayes_factor:: Float64
    ood             :: OODVerdict
    calibration     :: CalibrationMeta   # the grid's PASSED D-05 gate report
    meta            :: NamedTuple
end
```
Extension-point stubs (Phase 12 `SpatialColocResult`, Phase 13 `ThreeHypothesisColocResult`) are sketched but NOT implemented (RESEARCH §Pattern 1, D-02).

**Reusable metadata struct to port for `calibration::CalibrationMeta`:** the `CalibrationResult` struct already ported into the spike from BayesInteractomics (`spike/validation/sbc.jl:65-72`):
```julia
struct CalibrationResult
    bin_midpoints::Vector{Float64}
    predicted_rate::Vector{Float64}
    observed_rate::Vector{Float64}
    bin_counts::Vector{Int}
    ece::Float64
    mce::Float64
end
```

---

### `src/registry.jl` (store, CRUD) — PROD-02, D-04

**Analog (target shape):** RESEARCH §Pattern 2. **Analog (persistence mechanics):** `spike/npe/train_npe.jl` `save_npe`/`load_npe` (L152-201) and `spike/data/cache.jl` atomic idiom.

**Registry skeleton to build** (RESEARCH §Pattern 2 — copy directly):
```julia
struct EstimatorBundle
    grid       :: Int
    npe;  ratio;  ood_nulls;  zt;  θzt
    calibration :: CalibrationMeta
end
const _REGISTRY = Dict{Int, EstimatorBundle}()
const _SHIPPED_GRIDS = (4, 8, 16, 32)   # 64 DROPPED (D-04); only gate-PASSED grids populate

function estimator_for(grid::Integer)
    haskey(_REGISTRY, grid) && return _REGISTRY[grid]
    grid in _SHIPPED_GRIDS && (return _lazy_load_from_artifact!(grid))
    throw(ArgumentError(
        "No estimator registered for a $(grid)×$(grid) patch grid. " *
        "Shipped grids: $(_SHIPPED_GRIDS). Train and register with `train_and_register($grid)`."))
end
register!(b::EstimatorBundle) = (_REGISTRY[b.grid] = b)
train_and_register(grid::Integer; kwargs...) = (b = _train_grid_pipeline(grid; kwargs...); register!(b); b)
```

**Atomic-persistence pattern to REUSE for the artifact writer** — copy the `.tmp`+integrity-check+`mv` idiom verbatim (`spike/npe/train_npe.jl:160-180`):
```julia
function save_npe(path, result::NamedTuple; master_seed, fold)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp; schema_version = NPE_MODEL_SCHEMA, estimator = result.estimator,
        theta_transform = result.θzt, zt = result.zt, variant = result.variant,
        d_in = result.d_in, meta = (...))
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "estimator") "save_npe: integrity check failed, $tmp missing estimator"
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end
```
Same idiom appears 4× in the spike (`save_npe`, `save_ratio` at `train_ratio.jl:331-347`, `write_shard` at `cache.jl:97-108`, `write_meta` at `cache.jl:118-134`) — it is the established persistence convention.

**CRITICAL DELTA (Pitfall 4, RESEARCH Standard Stack):** the shipped registry must migrate the model bytes from **whole-object `estimator = result.estimator`** to **`Flux.state(cpu(result.estimator))`** + architecture metadata + `Flux.loadmodel!` on load, moving parameters to CPU before persistence (so artifacts load device-independently). Keep the atomic `.tmp`+integrity+`mv` wrapper unchanged.

**Content-hash cache-dir pattern** (if the registry versions caches per grid) — `spike/data/cache.jl:147-164` `open_or_invalidate` recomputes the hash and ERRORs on stored-vs-recomputed mismatch. Note the cache key already includes `summary_min_dim` (`spike/data/generate.jl:160`), so distinct grids auto-separate into distinct dirs (RESEARCH §Grid Coupling).

---

### `src/amortized/infer.jl` (service, request-response) — the amortized read surface

**Analog:** `spike/npe/infer.jl` — copy near-verbatim. Key functions to promote:

**`standardize_summary`** (`spike/npe/infer.jl:76-84`) — frozen-`zt` application, mask rows bypassed:
```julia
function standardize_summary(S::AbstractMatrix, zt, variant::Symbol = :min)
    cont, mask = _summary_row_partition(variant, size(S, 1))
    out = Matrix{Float32}(undef, size(S))
    out[cont, :] = Float32.(StatsBase.transform(zt, S[cont, :]))
    out[mask, :] = Float32.(S[mask, :])     # binary mask NEVER re-z-scored
    return out
end
```

**`posterior_for` / `rho_draws` / `delta_rho`** (`spike/npe/infer.jl:96-138`) — single-pass posterior, un-standardize with frozen `θzt` before any ρ read (Pitfall 5):
```julia
function posterior_for(est, Z; N::Integer = 2000, use_gpu::Bool = false)
    Zc = Z isa AbstractVector ? reshape(Z, :, 1) : Z   # coerce to d_in×1 = ONE dataset
    return sampleposterior(est, Zc; N = N, use_gpu = use_gpu)
end
function rho_draws(est, Z, θzt; N::Integer = 2000, use_gpu::Bool = false)
    draws = posterior_for(est, Z; N = N, use_gpu = use_gpu)
    orig  = StatsBase.reconstruct(θzt, draws)          # un-standardize (Pitfall 5)
    return vec(orig[1, :])                              # row 1 = ρ_true
end
function delta_rho(est, Zsample, Zcontrol, θzt; N = 2000, use_gpu = false)
    ρs = rho_draws(est, Zsample, θzt; N, use_gpu)
    ρc = rho_draws(est, Zcontrol, θzt; N, use_gpu)
    return mean(ρs .- ρc)                              # D-03 MC difference
end
```

**CRITICAL DELTA:** `use_gpu` is already a plumbed keyword here (defaults `false`) — keep the CPU default for the shipped/reproducible inference path (RESEARCH §Reproducibility split). The spike's `use_gpu=true` HARD-THROW guards (present in `train_npe.jl:114`, `train_ratio.jl:286`, `bf.jl:95`) must be REMOVED for the training path but the gate/inference default stays `false`.

**Also note the row-partition delegation** (`spike/npe/infer.jl:65`): `_summary_row_partition` delegates to the authoritative `Loader._row_partition` — grid-parametrization of that split (see below) automatically flows here.

---

### `src/amortized/train_npe.jl` + architecture (service, batch)

**Analog:** `spike/npe/architecture.jl` (`build_estimator`) + `spike/npe/train_npe.jl` (`train_fold`).

**`build_estimator`** (`spike/npe/architecture.jl:87-115`) — **already input-width-agnostic** (reads `d_in`), so NO change for grid coupling (RESEARCH §Grid Coupling "no change ✓"). The v0.2.1 API gotcha to preserve: `q` is a constructed `NormalisingFlow` INSTANCE passed POSITIONALLY:
```julia
layers = Dense[Dense(d_in, width, gelu)]
for _ in 2:depth; push!(layers, Dense(width, width, gelu)); end
push!(layers, Dense(width, dstar))
network = Chain(layers...)
q = NormalisingFlow(D; num_summaries = dstar, num_coupling_layers = num_coupling_layers,
                    depth = flow_depth, width = flow_width)
return PosteriorEstimator(network, q)   # INSTANCE positional, NOT q=Type
```
Architecture consts to carry (`architecture.jl:60-66`): `NPE_D=7`, `NPE_DSTAR=64`, `NPE_DEPTH=3`, `NPE_WIDTH=256`, `NPE_COUPLING=10`, `NPE_FLOW_DEPTH=2`, `NPE_FLOW_WIDTH=128` (the Phase-5-calibrated higher-capacity defaults).

**`train_fold`** (`spike/npe/train_npe.jl:102-149`) — fixed-data train form, leak-free θ standardization, `build_estimator(size(Ztr,1), ...)`:
```julia
θzt      = fit_theta_transform(fold_data.θtr)          # fit(ZScoreTransform, θtr; dims=2)
θtr_std  = Float32.(StatsBase.transform(θzt, fold_data.θtr))
d_in = size(fold_data.Ztr, 1)                          # data-driven; grid-agnostic
est  = build_estimator(d_in; dstar, depth, width, num_coupling_layers, flow_depth, flow_width)
est = train(est, θtr_std, θva_std, fold_data.Ztr, fold_data.Zva;
            epochs, batchsize, use_gpu = false,        # ← DELTA: plumb use_gpu for GPU training
            optimiser = Flux.Optimisers.AdamW(learning_rate, (0.9, 0.999), weight_decay),
            stopping_epochs, verbose)
```
**DELTA (D-06):** default `use_gpu = has_cuda_device()` for training; remove the `use_gpu && throw(...)` guard at `train_npe.jl:114`. Note the AdamW-Float64 gotcha comment (`train_npe.jl:129-131`): keep LR/decay Float64 to match NeuralEstimators' Float64 CosAnneal schedule.

---

### `src/amortized/bf.jl` + `train_ratio.jl` (service) — amortized Bayes factor (NRE)

**Analog:** `spike/validation/bf.jl` + `spike/validation/train_ratio.jl`.

**`amortized_log_bf`** — one forward pass, no quadgk/KDE (`spike/validation/bf.jl:94-100`):
```julia
function amortized_log_bf(est_bf, Z_pair, log_prior_odds; use_gpu::Bool = false)
    Zc = Z_pair isa AbstractVector ? reshape(Z_pair, :, 1) : Z_pair
    grid = reshape(Float32[0.0 1.0], 1, 2)                       # null (m=0), coloc (m=1)
    lr = logratio(est_bf, Zc; grid = grid, use_gpu = false)      # 1×2
    return (lr[1, 2] - lr[1, 1]) - log_prior_odds               # measured log_prior_odds subtracted
end
```

**`pair_encode`** (`spike/validation/train_ratio.jl:185-188`) — the A7 difference encoding; the grid-coupled site (`RATIO_CONT_ROWS`, `RATIO_INPUT_DIM`):
```julia
function pair_encode(Zs::AbstractVector, Zc::AbstractVector)
    d = @view(Zs[1:RATIO_CONT_ROWS]) .- @view(Zc[1:RATIO_CONT_ROWS])
    return Float32.(vcat(Zs, Zc, d))            # 128+128+64 = 320 at 8×8
end
```
**Grid-coupling DELTA (RESEARCH §Grid Coupling):** `RATIO_CONT_ROWS = G²`, `RATIO_INPUT_DIM = 2·(2G²)+G² = 5G²` (was 320 at G=8). The ratio net loss is hard-coded `logitbinarycrossentropy` — a custom loss passed to `train` is SILENTLY IGNORED (`train_ratio.jl:29-34`); do not attempt an l-POP loss.

**Ratio training** (`spike/validation/train_ratio.jl:278-321`) — same fixed-data `train`, cache-paired data (`assemble_ratio_data_cache`), measured `log_prior_odds` (`measure_log_prior_odds`, L240-244). Same `save_ratio`/`load_ratio` atomic idiom (L331-364).

**Reproduction-gate note:** the BF gate compares amortized log-BF vs the KDE baseline; the baseline is the demoted `src/bayes.jl::compute_BayesFactor` (L109-136) run through the Turing extension, on IDENTICAL Δρ draws (`bf.jl:208-237`). RESEARCH memo §5 ride-along: add a **non-clamped BF baseline** so `max|Δ logBF|` is testable without the KDE `1e-8` floor artifact.

---

### `src/amortized/ood.jl` (service) — OOD / misspecification flag

**Analog:** `spike/validation/ood.jl` — copy whole file. Zero new dependencies (hand-rolled Mahalanobis + ROC via LinearAlgebra/Statistics).

**Continuous-rows-only Mahalanobis** (`spike/validation/ood.jl:112-132`) — the Σ-singularity guard is grid-coupled through `_row_partition`:
```julia
function fit_ood_nulls(Ztrain::AbstractMatrix; variant::Symbol = :min)
    cont, _ = Loader._row_partition(variant, size(Ztrain, 1))   # continuous rows only
    Zc = Float64.(Ztrain[cont, :])
    μS = vec(mean(Zc; dims = 2)); ΣS = cov(Zc; dims = 2)
    C  = cholesky(Symmetric(ΣS + 1e-6 * I))                     # ridge (Pitfall 2)
    return (cont = cont, μS = μS, C = C)
end
```

**OR-fusion flag** (`spike/validation/ood.jl:375-380`), **hand-rolled ROC/AUC** (Mann-Whitney U, `L319-335`), **pre-registered `id_threshold`** (L344-345, the gate — NOT Youden-J). The **noise channel** (L251-303) computes features DIRECTLY from images (not the frozen summary) to narrow the detector-noise blind spot.

**Finite-guard to preserve** (`spike/validation/ood.jl:145-159`) — `_finite_or`/`_theta_tuple` keep the PP channel computable on misspecified inputs (memo §5 hardening item: re-enable the PP channel in the reported OR-fusion, currently `with_pp=false` at `ood.jl:590`).

---

### Grid-parametrized summary/encoder (utility, transform) — the ONE shared edit

**Analog:** `spike/contract.jl:85-90` (`patch_summary`) + `spike/data/encode.jl:71-75` (`encode_d01`).

**Current (8 hardcoded):**
```julia
# spike/contract.jl:85-90
function patch_summary(mci::MultiChannelImage)
    x = mci.data[1]; y = mci.data[2]
    xp, yp = patch.([x, y], 8)                        # ← hardcoded 8
    return correlation(xp, yp; method = :pearson)
end
# spike/data/encode.jl:71-75
function encode_d01(M::AbstractMatrix)
    vals = vec(coalesce.(M, 0.0))                     # 64-dim at 8×8
    mask = vec(Float64.(.!ismissing.(M)))             # 64-dim
    return vcat(vals, mask)                            # 128-dim
end
```

**Target (grid `G` a parameter, RESEARCH §Code Examples):**
```julia
function patch_summary(mci::MultiChannelImage, G::Integer)
    x, y = mci.data[1], mci.data[2]
    xp, yp = patch.([x, y], G)                        # G-parametric; patch() already accepts any G
    return correlation(xp, yp; method = :pearson)     # G×G Union{Float64,Missing}
end
# encode_d01 unchanged in body (vec order handles any G); dim becomes 2·G²
```

**Underlying src is already grid-general:** `src/colocalization.jl:37` `patch(img, num_patches)` accepts any grid; the ≥15-survivor floor lives in `correlation`/`_exclude_zero` (RESEARCH §Feasibility cites `colocalization.jl:235`). The existing ADVI path `_prepare_data(img, channels, num_patches)` (`src/bayes.jl:164`) already threads `num_patches` through `patch.([x,y], num_patches)` (L171) — the amortized path mirrors this input contract.

**Row-partition coupling** (`spike/data/loader.jl:126-136`) — the `:min` split is hardcoded `1:64`/`65:128`:
```julia
if variant === :min
    nrows == 128 || error("_row_partition: :min expects 128 rows, got $nrows")
    return (collect(1:64), collect(65:128))           # ← DELTA: (1:G², G²+1:2G²)
```
**Grid-coupling site table (RESEARCH §Grid Coupling), all to generalize once in Wave 0:**
| Site | 8×8 hardcoded | Generalize to `G` |
|------|---------------|-------------------|
| `patch_summary` | `patch.(…, 8)` | `patch.(…, G)` |
| `encode_d01` | 128 = 64+64 | `2G²` (vals 1:G², mask G²+1:2G²) |
| `generate.jl` buffers (`generate.jl:108-110,172-174`) | `Matrix(128, N)` / `Matrix(AUG_DIM, N)` | `Matrix(2G², N)` |
| `_row_partition` | `1:64`, `65:128` | `1:G²`, `G²+1:2G²` |
| `pair_encode`/`RATIO_INPUT_DIM` | 320 = 128+128+64 | `5G²` |
| `build_estimator` | reads `size(Ztr,1)` | no change ✓ |
| cache `generating_config` (`generate.jl:160`) | `summary_min_dim = 128` | `2G²` (auto-separates caches) |

---

### `colocalization_amortized()` public entry (controller, request-response)

**Analog (behavior):** the `spike/npe/infer.jl` read surface composed end-to-end.
**Analog (signature/entry-point shape):** existing `src/bayes.jl::colocalization()` (L234-242):
```julia
function colocalization(
    img::MultiChannelImageStack, control::MultiChannelImageStack,
    channels::Vector{T}, num_patches::T = 1;
    iter::T = 1000, posterior_samples::T = 100_000,
    cor_method::Symbol = :pearson) where T <: Int
```
**Target (RESEARCH §Code Examples — name is discretion, `num_patches` is the registry key):**
```julia
function colocalization_amortized(img, control, channels; num_patches::Integer = 8, N = 2000)
    b   = estimator_for(num_patches)                                   # registry (PROD-02)
    Zs  = standardize_summary(encode_d01(patch_summary(img,     num_patches)), b.zt)
    Zc  = standardize_summary(encode_d01(patch_summary(control, num_patches)), b.zt)
    post = posterior_for(b.npe, Zs; N = N)
    Δρ   = rho_draws(b.npe, Zs, b.θzt; N) .- rho_draws(b.npe, Zc, b.θzt; N)
    logbf = amortized_log_bf(b.ratio, pair_encode(Zs, Zc), b.ratio_log_prior_odds)
    ood   = ood_verdict(b.ood_nulls, Zs, Zc)
    return AmortizedColocResult(num_patches, reconstruct(b.θzt, post), Δρ, logbf, ood,
                                b.calibration, (; N))
end
```
Reuses the exact frozen read chain the harness rides (`spike/validation/harness.jl:106-113`).

---

### `ext/ProteinCoLocTuringExt.jl` (weakdep extension) — Finding-1 co-resolution fix

**Analog:** the entire existing Turing body of `src/bayes.jl` (the `@model`, `vi(m, ADVI(...))` at L325, `convert_posterior_samples` L42-54, `compute_BayesFactor` L109-136) + the `using Turing` imports in `src/ProteinCoLoc.jl:35-36`.

**Current module imports to split** (`src/ProteinCoLoc.jl:22-52`):
```julia
import KernelDensity: kde       # L30  → move to ext (compute_BayesFactor)
import QuadGK: quadgk           # L31  → move to ext
using Turing                    # L35  → weakdep
using Turing: Variational       # L36  → weakdep
export … colocalization … compute_BayesFactor … CoLocResult   # L49-51 → UNEXPORT (D-01)
```
**Target (RESEARCH §Pattern 3):** move the ADVI `colocalization()`, `compute_BayesFactor()`, `convert_posterior_samples()`, `_prepare_data()` and the `AdviColocResult` (renamed `CoLocResult`) into `ext/ProteinCoLocTuringExt.jl`, gated by:
```toml
# root Project.toml
[weakdeps]
Turing = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"
[extensions]
ProteinCoLocTuringExt = "Turing"
```
**Anti-pattern (RESEARCH):** do NOT `export` the Turing path or `AdviColocResult` — D-01 makes it internal reference only.

---

### `Project.toml` (config) — dependency surgery + version bump

**Analog:** current root `Project.toml` (13 deps, `version = "1.0.1"`, `GLMakie = "0.10.5"` compat pin) + the proven spike resolve set.

**Deltas (D-03, RESEARCH §Standard Stack + Runtime State Inventory):**
- Add hard deps: `NeuralEstimators` (pin **0.2.1**), `Flux` (pin **0.16.10**), `JLD2`, `HypothesisTests` (Random123/StatsBase/Distributions/Images already present).
- Move `Turing` → `[weakdeps]` + `[extensions]` (Finding 1 — Turing+GLMakie caps NeuralEstimators to 0.1.4). Consider `GLMakie` → extension too if resolve still fails.
- Add `CUDA` as an **optional weakdep** (GPU training only, graceful CPU fallback — D-06).
- Bump `version = "1.0.1"` → `"2.0.0"` (breaking, D-01).
- **HARD GATE (Wave 0):** `Pkg.resolve` must pin NeuralEstimators 0.2.1 (assert in `test/runtests.jl`). This is the single most plan-reshaping risk (RESEARCH A3, Pitfall 1) — de-risk BEFORE any src inference code.

---

### `test/gate/run_gate.jl` + harness (test, batch) — the per-grid ship-gate (D-05)

**Analog:** `spike/validation/harness.jl` + `spike/validation/sbc.jl`.

**Harness core** (`spike/validation/harness.jl:106-136`) — the θ*~π→simulate→infer chain and paired-Δρ path, all frozen-stats + CPU:
```julia
function draw_simulate_infer(m, rng; imsize = SBC_IMSIZE, N = SBC_L)
    θ   = sample_prior(rng)
    mci = build_mci(simulate_pair(rng, θ; imsize = imsize))
    Z   = standardize_summary(encode_d01(patch_summary(mci)), m.zt, :min)   # ← DELTA: patch_summary(mci, G)
    draws_std = posterior_for(m.estimator, Z; N = N, use_gpu = false)       # CPU gate (mandatory)
    draws     = StatsBase.reconstruct(m.θzt, draws_std)
    return (θ = θ, Z = Z, draws = draws)
end
```

**SBC rank table** (`spike/validation/sbc.jl:156-172`) — M×8 (7 θ + Δρ column), rank = `count(<(θ*), draws)`:
```julia
for p in 1:7; rank_table[i, p] = count(<(t.θ[p]), @view t.draws[p, :]); end
pr = draw_simulate_infer_paired(m, rng; imsize, N = L)
rank_table[i, 8] = count(<(pr.θs.ρ_true - pr.θc.ρ_true), pr.ρs .- pr.ρc)   # dedicated Δρ SBC
```

**KS/χ² uniformity** (`spike/validation/sbc.jl:184-196`) — HypothesisTests, never hand-rolled:
```julia
u = (ranks .+ 0.5) ./ (L + 1)
ks_p = pvalue(ExactOneSampleKSTest(u, Uniform(0.0, 1.0)))
# … χ² on equal-width rank bins → pvalue(ChisqTest(counts))
```
**CRITICAL DELTAS:** (1) every gate call stays **`use_gpu = false`** against a CPU-resident frozen net (RESEARCH §CPU-reproducibility — the pre-registered numbers must be deterministic regardless of training hardware). (2) `patch_summary(mci)` → `patch_summary(mci, G)`. (3) load the per-grid net, not the fixed `spike/npe/trained_npe.jld2`.

---

### `test/gate/gate_consts_G.jl` (config, per grid) — fresh pre-registration (D-05)

**Analog:** `spike/validation/consts.jl` (the values), `spike/validation/harness.jl:60-94` (the disjoint-seed salt idiom).

**Template to re-express FRESH per grid** (`spike/validation/consts.jl:42-80`): `SBC_M`, `SBC_L`, `SBC_BINS`, `SBC_KS_ALPHA`, `SBC_ECE_GREEN/YELLOW`, `SBC_IMSIZE`, `BF_CORR_MIN`, `BF_LOGBF_TOL`, `BF_SWEEP_*`, `OOD_ID_QUANTILE`, `OOD_AUC_MIN`, `OOD_GRID_LEVELS`, `OOD_PP_REPS`.

**CRITICAL DELTA (Pitfall 3 — anti-snooping):** each grid gets a **new disjoint `PROD_SEED[G]`**, NOT the spike's `VAL_MASTER_SEED = 0x5BC0FFEE` (`consts.jl:76`) and NOT `NPE_MASTER_SEED = 0xC0FFEE`. The net was iterated against `VAL_MASTER_SEED`, so reusing it is data-snooping. Also raise `SBC_IMSIZE` per grid so fine grids clear the ≥15-survivor floor (Pitfall 2). The Random123 keyed-stream idiom to reuse (`harness.jl:93-94`):
```julia
val_rng(master_seed = VAL_MASTER_SEED) =
    Philox4x(UInt64, (UInt64(master_seed) ⊻ VAL_SALT, UInt64(0)))   # provably disjoint stream
```

---

## Shared Patterns

### Atomic JLD2 persistence (`.tmp` + reopen-integrity-check + `mv force=true`)
**Source:** `spike/data/cache.jl:97-108` (`write_shard`), `spike/npe/train_npe.jl:160-180` (`save_npe`).
**Apply to:** the registry artifact writer, any per-grid cache, the gate-report writer.
```julia
tmp = path * ".tmp"
jldsave(tmp; <keys...>)
JLD2.jldopen(tmp, "r") do f; @assert haskey(f, "<key>") "integrity check failed"; end
mv(tmp, path; force = true)   # filesystem-atomic commit; a crash leaves only a discardable .tmp
```
**DELTA for shipped models:** persist `Flux.state(cpu(est))` + arch metadata (not the whole object); reload via `Flux.loadmodel!` (Pitfall 4).

### Frozen-stats discipline (never re-fit standardization at inference)
**Source:** `spike/npe/train_npe.jl:80` (fit `θzt` on train θ only), `spike/npe/infer.jl:76-84` (`standardize_summary` applies the frozen `zt`, mask rows bypassed), `spike/validation/harness.jl:29-32` (invariant banner).
**Apply to:** `colocalization_amortized`, every gate call, ratio training, OOD scoring. The `EstimatorBundle` carries the frozen `zt` + `θzt`; downstream ALWAYS `StatsBase.transform`/`reconstruct` with them, never `fit(...)`.

### Disjoint Random123 keyed streams (reproducible seeding)
**Source:** `spike/validation/harness.jl:60-94`, `spike/validation/consts.jl:72-79`.
**Apply to:** every per-grid gate (`PROD_SEED[G]`), data-gen, ratio pairing — each provably disjoint (XOR-salt) from training and from each other. Never `Random.seed!` (RESEARCH §Don't Hand-Roll).

### Hand-rolled OOD ROC/AUC + continuous-rows Mahalanobis (no new dependency)
**Source:** `spike/validation/ood.jl:112-132, 319-335`.
**Apply to:** the promoted `src/amortized/ood.jl` — keeps LinearAlgebra/Statistics-only footprint, no ROC package enters the resolve.

### `use_gpu` reproducibility split (train GPU, gate/infer CPU)
**Source:** the `use_gpu` keyword is already plumbed through `sampleposterior`/`logratio`/`train` (`spike/npe/infer.jl:96`, `spike/validation/bf.jl:98`, `spike/npe/train_npe.jl:138`); the spike forces `false` everywhere.
**Apply to:** training paths default `use_gpu = has_cuda_device()` (D-06); gate + shipped inference default `use_gpu = false` (deterministic pre-registered numbers). Remove the spike's `use_gpu && throw(...)` guards (`train_npe.jl:114`, `train_ratio.jl:286`, `bf.jl:95`) on the training path only.

---

## No Analog Found

None. Every Phase-7 file maps to a proven spike/ or existing src/ analog — this phase is a copy-and-adapt promotion, not new-algorithm work (RESEARCH §Don't Hand-Roll "Key insight"). The genuinely NEW *design* work (not new code patterns) is:
- The `AbstractColocResult` supertype + accessor interface (D-02) — no prior abstract hierarchy exists in src/; `CoLocResult` is a bare struct. The interface shape is fully specified in RESEARCH §Pattern 1.
- `Flux.state`/`loadmodel!` persistence — a mechanism SWAP inside the existing atomic-write wrapper (the wrapper is the analog; the serialization call is the delta).
- The `ext/` weakdep extension file — a MOVE of existing `src/bayes.jl` code into Julia's package-extension structure (the code is the analog; the packaging is new).

---

## Conventions

Convention derivation skipped (`gsd-tools verify conventions --derive` returned `skipped: no-readable-files` — the deterministic module targets JS/TS sources and does not parse Julia `.jl` files). The repo's conventions are instead read directly from the codebase and stated below as observed (not machine-derived).

| Axis | Dominant (observed) | Share | Status |
|------|---------------------|-------|--------|
| File-name casing | lowercase / snake-ish (`train_npe.jl`, `bayes.jl`); PascalCase for module-named files (`LoadImages.jl`, `ProteinCoLoc.jl`) | ~mixed by role | named contract (per-role) |
| Identifier casing | `snake_case` functions (`patch_summary`, `build_estimator`, `amortized_log_bf`); `PascalCase` types (`CoLocResult`, `EstimatorBundle`, `MultiChannelImage`); `SCREAMING_SNAKE` consts (`NPE_DSTAR`, `VAL_MASTER_SEED`) | dominant | named contract |
| Export style | explicit `export` list at module bottom (`src/ProteinCoLoc.jl:49-52`); flat top-level functions in spike files, no per-file module wrapper | dominant | named contract |
| Include style | `include(joinpath(@__DIR__, ...))` with `isdefined(@__MODULE__, :sym) || include(...)` idempotency guards (spike); plain ordered `include("...")` in `src/ProteinCoLoc.jl` | two contexts | contested hotspot (see below) |

**Contested hotspots (author's choice).** The **spike/ vs src/ include discipline is the prototype intentional-contested split**, directly analogous to the CJS↔SDK dual-resolver pattern: `spike/**` files are **standalone-loadable, guarded, ordered `isdefined || include`** modules-of-functions (so each file loads under both `demo.jl` and `runtests.jl` without redefinition warnings); `src/**` is a **single `module ProteinCoLoc` with plain top-level `include`s and one export list**. Each half is internally consistent per-directory, contested only repo-wide. **Planner/reviewer guidance:** promoted `src/amortized/*.jl` files must adopt the **src/ module convention** (no per-file `isdefined` guards, rely on the module's ordered includes and `export`), NOT the spike's standalone-guard idiom — match the destination directory's local style, not the source file's. Likewise: `use_gpu` guards, `VAL_MASTER_SEED`, and flat-function layout are spike-local conventions that must be TRANSLATED, not copied, on promotion.

---

## Metadata

**Analog search scope:** `spike/npe/`, `spike/validation/`, `spike/data/`, `spike/contract.jl`, `src/bayes.jl`, `src/colocalization.jl`, `src/ProteinCoLoc.jl`, root + spike `Project.toml`.
**Files scanned:** 18 (13 read in full or targeted, 5 grepped for structure).
**Pattern extraction date:** 2026-07-03
