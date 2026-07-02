# Phase 4: NPE Training + ADVI Benchmark + Ablation - Pattern Map

**Mapped:** 2026-06-30
**Files analyzed:** 12 (10 new, 2 modified)
**Analogs found:** 11 / 12 (1 new-dep harness has only a partial analog)

Every new file has a strong in-repo analog: the spike already contains the NPE
construction call (`00_smoke.jl`), the leak-free loader the trainer consumes
(`loader.jl`), the re-simulate-from-keyed-RNG path the ADVI baseline needs
(`generate.jl` `_write_holdout`), the atomic JLD2 + integrity-check write
(`cache.jl`), the read-only `src/` `include()` boundary (`contract.jl`), and the
`@testset` + `const`-fixture gate harness (`test_data_pipeline.jl`). The genuinely
new surface is wiring, not algorithms.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `spike/npe/architecture.jl` | component (builder) | transform | `spike/00_smoke.jl` L60-63 | exact |
| `spike/npe/train_npe.jl` | service (training) | batch / transform | `spike/00_smoke.jl` L66-67 + `spike/data/loader.jl` | exact |
| `spike/npe/infer.jl` | service (inference) | transform / request-response | `spike/00_smoke.jl` L72-74 + `spike/simulator/ghat.jl` | exact |
| `spike/npe/benchmark.jl` | service (benchmark) | batch | `spike/data/generate.jl` (threading) + RESEARCH Code Examples | partial (BenchmarkTools NEW) |
| `spike/npe/ablation.jl` | service (CV scoring) | batch | `spike/data/loader.jl` `variant=` + `test_data_pipeline.jl` fold loop | role-match |
| `spike/baseline/Project.toml` | config | — | `spike/Project.toml` | exact |
| `spike/baseline/model.jl` | model (read-only lift) | — | `spike/contract.jl` L45-47 | exact |
| `spike/baseline/run_advi.jl` | service (port + migration) | file-I/O | `src/bayes.jl` `colocalization()` L234-339 + `generate.jl` `_write_holdout` L171-198 | exact (port target) |
| `spike/baseline/advi_artifact.jld2` | data artifact | file-I/O | `holdout.jld2` schema (`generate.jl` L191) | role-match |
| `spike/test/test_npe.jl` | test | — | `spike/test/test_data_pipeline.jl` | exact |
| `spike/test/runtests.jl` (MOD) | test (harness) | — | `spike/test/runtests.jl` L52-86 | exact (self) |
| `spike/Project.toml` (MOD) | config | — | `spike/Project.toml` | exact (self) |

## Pattern Assignments

### `spike/npe/architecture.jl` (component, transform)

**Analog:** `spike/00_smoke.jl` L60-63 — the only verified-working `PosteriorEstimator` construction in the repo.

**Core pattern** (`00_smoke.jl:60-63`):
```julia
# summary network: maps an m-vector data set to num_summaries summaries.
network = Chain(Dense(m, 64, gelu), Dense(64, 64, gelu), Dense(64, num_summaries))
q       = NormalisingFlow(d; num_summaries = num_summaries)
est     = PosteriorEstimator(network, q)
```
Copy directly. For Phase 4: `m → d_in` (128 `:min` / 142 `:aug`, from `size(fold.Ztr,1)`), `d → D=7` (the 7-field θ; ρ_true is row 1), `num_summaries → DSTAR`. The smoke file's docstring (L31-37) records the two API gotchas to preserve: **`q` is a constructed `NormalisingFlow` INSTANCE passed POSITIONALLY** (the `q=` keyword convenience form expects a TYPE), and `NormalisingFlow(d; num_summaries=dstar)`.

**Const-fixture convention** (`00_smoke.jl:45-48`): module-level `const` for every architecture dimension declared before use:
```julia
const d             = 1
const num_summaries = 8
```
→ Phase 4 declares `const D = 7`, `const DSTAR = 32` (RESEARCH Pattern 1).

---

### `spike/npe/train_npe.jl` (service, batch/transform)

**Analog (train call):** `spike/00_smoke.jl` L66-67. **Analog (fold input):** `spike/data/loader.jl` `load_fold`.

**The CPU-only train call** (`00_smoke.jl:65-67`):
```julia
# CPU-only training (D-04: use_gpu = false is mandatory).
est = train(est, prior_sampler, gaussian_simulator;
            K = 8000, epochs = 60, use_gpu = false, verbose = false)
```
**Critical change for Phase 4:** the smoke used the *online* `train(est, sampler, simulator; ...)` form. Phase 4 MUST use the *fixed-data* form on the leak-free cache (RESEARCH Pattern 1, "Alternatives Considered"):
```julia
est = train(est, fold.θtr, fold.θva, fold.Ztr, fold.Zva;
            epochs=200, batchsize=64, use_gpu=false, verbose=false,
            optimiser=Flux.Optimisers.AdamW(5e-4, (0.9,0.999), 1e-4),
            stopping_epochs=10)
```
`use_gpu` **defaults to true** — every `train`/inference call must pass `use_gpu=false` (RESEARCH Pitfall 1, anti-pattern). The smoke comment at L37-38 documents this exact gotcha.

**Loader input shape** (`loader.jl:152-195`): `load_fold(dir, fold; K=5, master_seed, variant=:min|:aug)` returns `(Ztr, θtr, Zva, θva, zt)` — column-major Float32 `d×K` tensors that drop straight into `train`. **Reuse `fold.zt`** (the frozen train-only `ZScoreTransform`) for any later holdout standardization; never re-standardize (`loader.jl:160-163` — `zt` is returned precisely so downstream phases freeze preprocessing).

**θ-standardization (NEW, RESEARCH Pitfall 5):** the loader standardizes only `Z`, not `θ` (`loader.jl:193` returns raw `Float32.(θ[...])`). If θ leak-free standardization is adopted, mirror the loader's fit-on-train discipline exactly (`loader.jl:182` `fit(ZScoreTransform, ...; dims=2)` on train columns only) and un-standardize draws before RMSE/`ghat`.

**Module/guarded-include structure** for any new multi-include spike file — copy `generate.jl:42-51`:
```julia
isdefined(@__MODULE__, :load_fold) || include(joinpath(@__DIR__, "..", "data", "loader.jl"))
```
Guarded `isdefined`/`isfile` includes keep the file loadable both standalone and inside `runtests.jl`.

---

### `spike/npe/infer.jl` (service, transform / request-response)

**Analog (sampling):** `spike/00_smoke.jl` L72-74. **Analog (μ→ρ map):** `spike/simulator/ghat.jl`.

**Single-pass posterior** (`00_smoke.jl:72-74`):
```julia
draws  = sampleposterior(est, Z_obs; N = 1000)   # d x N matrix; N is a KEYWORD
mu_hat = posteriormean(draws)[1]                 # mean over the N draws
```
For Phase 4: `ρ̂ = posteriormean(draws)[1]` (row 1 = ρ_true). Δρ = MC-diff of two single-stack passes (D-03, RESEARCH Pattern 2):
```julia
Δρ_npe = mean(ρ_npe_sample_draws .- ρ_npe_control_draws)
```

**ρ-space alignment via the frozen ghat** (`ghat.jl:59-67`) — the cross-method common axis (D-04 / RESEARCH Pattern 2). `ghat` is FROZEN and auto-generated (`ghat.jl:22-23` "DO NOT EDIT BY HAND"); call it, never recompute:
```julia
ρ̂_advi = mean(ghat.(μ_sample_samples))   # monotone map, applied per draw
```
`ghat` clamps μ outside `[GHAT_MU_MIN, GHAT_MU_MAX]` (`ghat.jl:48-49`); ADVI draws must land in μ∈[−1,1] first (RESEARCH Pitfall 3).

**Don't hand-roll** (RESEARCH "Don't Hand-Roll"): use `interval(d; probs=[0.05,0.95])` for interval width and `assess(est, θ_test, Z)` → `rmse(assessment)` for per-parameter RMSE — never manual quantile/error loops.

---

### `spike/npe/benchmark.jl` (service, batch) — PARTIAL ANALOG

**Analog (threading + re-sim inputs):** `spike/data/generate.jl`. **No in-repo BenchmarkTools usage exists** — BenchmarkTools is a NEW spike dep (RESEARCH Standard Stack). Use the RESEARCH "Code Examples → BenchmarkTools harness (D-08)" block as the template:
```julia
using BenchmarkTools
t_npe = @belapsed begin
    Zs = $(standardize)(encode_d01(patch_summary($sample_mci)), $fold.zt)
    sampleposterior($est, Zs; N=$N)
end
speedup = t_advi / t_npe    # ALWAYS reported with RMSE (D-08/D-09)
```

**Threading / scaling analog** (`generate.jl:103-137`): the `parallel`/serial byte-identical pattern and per-global-index keying are the template for the D-12 scaling sweep and D-13 thread sweep. Note `Threads.nthreads()` is fixed per process (RESEARCH Pitfall 7) — the thread sweep runs as separate `julia -t N` processes writing JLD2/CSV, NOT inside one session.

**Re-simulate holdout raw stacks** (the D-08 raw-stack→posterior input) by copying the keyed re-sim loop from `generate.jl:177-187` (see `run_advi.jl` below) so NPE and ADVI time the SAME raw `MultiChannelImage`.

**Pre-registered `const` gates** (mirror `test_data_pipeline.jl:60-67` fixture style; RESEARCH Validation Architecture):
```julia
const NPE_RMSE_TOLERANCE = 1.2; const SPEEDUP_GATE = 100.0
const ABL_REL_MARGIN = 0.05; const ABL_FOLD_CONSISTENCY = 4
const BENCH_THREADS = 1; const NPE_MASTER_SEED = 0xC0FFEE
```

---

### `spike/npe/ablation.jl` (service, batch CV scoring)

**Analog:** `spike/data/loader.jl` `variant=` switch + the `all_folds`/`load_fold` k=5 loop, scored with NeuralEstimators `assess`/`rmse`.

**Variant switch — zero re-simulation** (`loader.jl:166-170`): both summaries are already cached; ablation just toggles `variant`:
```julia
fold_min = load_fold(dir, f; K=5, master_seed=NPE_MASTER_SEED, variant=:min)
fold_aug = load_fold(dir, f; K=5, master_seed=NPE_MASTER_SEED, variant=:aug)
```
Same architecture, same k=5 CV, per-parameter `rmse(assess(est, θva, Zva))`.

**Decision rule (NEW, pre-registered):** `:aug` wins iff `RMSE_aug ≤ (1−ABL_REL_MARGIN)·RMSE_min` on ρ_true AND improves in `≥ABL_FOLD_CONSISTENCY` of 5 folds, else keep `:min` (parsimony + OOD; D-07, RESEARCH A3). Couples to Phase 5 (ABL-02 OOD note).

---

### `spike/baseline/Project.toml` (config)

**Analog:** `spike/Project.toml` — same `[deps]` UUID-list format. The baseline gets its OWN isolated env (D-01): `Turing`, `DataFrames`, `JLD2`, `ForwardDiff`, `StatsBase`. **Turing MUST NEVER enter `spike/Project.toml`** — co-resolve caps NeuralEstimators below the pinned 0.2.1 (RESEARCH anti-pattern; Phase-1 landmine). Pin + commit `spike/baseline/Manifest.toml` as a reproducibility artifact, exactly as `spike/Manifest.toml` is pinned.

---

### `spike/baseline/model.jl` (model, read-only lift)

**Analog:** `spike/contract.jl` L45-47 — the established read-only `src/` `include()` boundary.

**The read-only coupling discipline** (`contract.jl:45-47`):
```julia
# --- include() coupling boundary -- READ-ONLY, never edit src/ -----------------
include(joinpath(@__DIR__, "..", "src", "LoadImages.jl"))     # MultiChannelImage(Stack)
include(joinpath(@__DIR__, "..", "src", "colocalization.jl")) # patch, correlation, _exclude_zero
```
The baseline reaches `src/bayes.jl`'s `@model` only by read-only copy/`include()` — never edits `src/` (CLAUDE.md hard constraint; spike decoupling holds until Phase 7). The `@model` ports VERBATIM (spike 003 verified it captures no enclosing variable); only the inference call changes (see `run_advi.jl`).

**The @model to lift** (`src/bayes.jl:267-303`): the hierarchical Turing model with `μ_control`/`μ_sample ~ Truncated(Cauchy(0,0.3),-1,1)` globals (L274/L279) + per-image locals + TDist likelihood (L294-302). `μ_sample`/`μ_control` are the per-condition means the baseline reads and maps through `ghat`.

---

### `spike/baseline/run_advi.jl` (service, port + migration, file-I/O)

**Analog (port target):** `src/bayes.jl` `colocalization()` L234-339. **Analog (re-sim holdout):** `generate.jl` `_write_holdout` L171-198. **Analog (atomic write):** `cache.jl` `write_shard`/`write_meta`.

**The removed API to port** (`src/bayes.jl:325`):
```julia
q = vi(m, ADVI(num_latent, iter))     # REMOVED pre-0.35 — must be ported
```
→ modern AdvancedVI (RESEARCH Pattern 4):
```julia
result = vi(m, q_meanfield_gaussian, ITER; adtype=AutoForwardDiff(), show_progress=false)
q      = result.q                      # VIResult: .q + .info
draws  = rand(q, POSTERIOR_SAMPLES)    # ASSERT extrema(μ draws) ⊂ [−1,1] (Pitfall 3)
```
Also: `num_latent` via `DataFrames.DataFrame(prior_chain)` (`bayes.jl:323`) and the `syms`/name extraction are removed → use `keys(DynamicPPL.VarInfo(m))` or the 8 known globals (`bayes.jl:313-316`). The baseline should time ONLY the `vi()` call (D-08) and may SKIP the 100k `sample(m, Prior(), posterior_samples)` prior chain (`bayes.jl:308`) — it is dead/wasteful cost (RESEARCH Pattern 3).

**Re-simulate the holdout raw stacks** — copy the keyed re-sim loop from `generate.jl:177-187` so ADVI runs on byte-identical raw images (the cache stores summaries, not raw `MultiChannelImage`s; RESEARCH Pitfall 4):
```julia
rng = holdout_rng(master_seed, j)               # disjoint key namespace (D-10)
θ   = sample_prior(rng); isz = sample_imsize(rng)
mci = build_mci(simulate_pair(rng, θ; imsize=isz))
```
Wave-0 must prove round-trip reproducibility (re-sim → summary == cached summary; RESEARCH Open Question 1).

**Atomic + integrity-checked JLD2 artifact** — copy `cache.jl write_shard` (the `.tmp` → reopen-assert → `mv(...; force=true)` idiom) and the schema/namespace tagging from `generate.jl:191-196`:
```julia
jldsave(tmp; schema_version=1, meta=(turing_version=…, adtype=…, iter=ITER,
        master_seed=MASTER_SEED, bayes_git_ref="<sha>"),
        pairs=pair_ids, mu_sample=mu_s, mu_control=mu_c,
        rho_sample=ghat.(mu_s), rho_control=ghat.(mu_c), wall_clock=t_vi)
JLD2.jldopen(tmp,"r") do f; @assert haskey(f,"mu_sample"); end
mv(tmp, path; force=true)
```
Store raw μ **sample vectors** (not just means) so the spike can recompute Δρ by MC differencing (RESEARCH Code Examples).

---

### `spike/test/test_npe.jl` (test)

**Analog:** `spike/test/test_data_pipeline.jl` — copy its structure exactly.

**File skeleton** (`test_data_pipeline.jl:1-70`): GNU-AGPL license header `#= … =#` → `using Test` (+ `Random`, `Statistics`, `BenchmarkTools`) → guarded `include`s of the units under test → pre-declared `const` fixtures → one outer `@testset "Phase 4 — …" verbose=true` wrapping named child testsets.

**Const-fixture-before-gate discipline** (`test_data_pipeline.jl:55-67`):
```julia
const DP_LOAD_SEED = 3; const DP_K = 5
const DP_MEAN_TOL  = 1e-6; const DP_STD_TOL = 1e-3   # tolerances NOT tuned to pass
const DP_LOAD_DIR  = generate_cache(mktempdir(); N=DP_N_LOAD, master_seed=DP_LOAD_SEED, n_holdout=20)
```
→ Phase 4 declares `NPE_RMSE_TOLERANCE`, `SPEEDUP_GATE`, `ABL_REL_MARGIN`, `ABL_FOLD_CONSISTENCY`, `BENCH_THREADS`, `NPE_MASTER_SEED` BEFORE any reported run (RESEARCH A2/A3).

**Wave-0 skipped-placeholder pattern** (`test_data_pipeline.jl:70-91`): each SC testset ships as a named `@testset` holding `@test_skip true` tagged with the wave that fills it, so the single `runtests.jl` gate enumerates SC1..SC5 from Wave 0. Map: SC1=NPE-01, SC2=NPE-02, SC3=NPE-03 (+ scaling), SC4=ABL-01, SC5=ABL-02 (RESEARCH Test Map).

**Surface-don't-swallow assertion** (`test_data_pipeline.jl:87-90`): the `@test_throws ArgumentError simulate_pair(...)` pattern — a bad input must raise, never corrupt a column (ASVS V5). Preserve for any new validating path.

---

### `spike/test/runtests.jl` (MODIFIED, test harness)

**Analog:** itself, L52-86 — extend the existing resolve-risk gate + child-include chain.

**Add the BenchmarkTools resolve-risk assertion** alongside the existing JLD2/Random123 gate (`runtests.jl:65-74`):
```julia
@test haskey(_deps_by_name, "BenchmarkTools")    # NEW Wave-0 dep present
@test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
```
The NeuralEstimators v0.2.1 pin (UUID-keyed, `runtests.jl:64`) must be re-asserted with BenchmarkTools present (RESEARCH Pitfall 2). **Add the include** at the bottom, mirroring L81-86:
```julia
include(joinpath(@__DIR__, "test_npe.jl"))
```

---

### `spike/Project.toml` (MODIFIED, config)

**Analog:** itself. Add one `[deps]` line — `BenchmarkTools = "<uuid>"` — alongside the existing UUID list (note: Turing/AdvancedVI go in `spike/baseline/Project.toml`, NEVER here). Re-freeze + commit `spike/Manifest.toml` after `Pkg.add`.

## Shared Patterns

### Read-only `src/` coupling (decoupling, CLAUDE.md hard constraint)
**Source:** `spike/contract.jl:45-47`
**Apply to:** `spike/baseline/model.jl` (and any baseline access to `src/bayes.jl`)
```julia
include(joinpath(@__DIR__, "..", "src", "colocalization.jl"))  # READ-ONLY, never edit src/
```
The spike and baseline reach `src/` only by read-only copy/`include()`; no `src/` edits until Phase 7.

### Guarded idempotent includes
**Source:** `spike/data/generate.jl:42-53`, `test_data_pipeline.jl:46-53`
**Apply to:** every new multi-include file (`train_npe.jl`, `infer.jl`, `benchmark.jl`, `ablation.jl`, `test_npe.jl`)
```julia
isdefined(@__MODULE__, :load_fold) || include(joinpath(@__DIR__, "..", "data", "loader.jl"))
```
Keeps files loadable both standalone and inside `runtests.jl`.

### Random123 keyed-per-index reproducibility
**Source:** `spike/data/seeding.jl:75-96`
**Apply to:** `run_advi.jl` (holdout re-sim), `benchmark.jl` (scaling/thread sweeps)
```julia
sample_rng(ms, i)  = Philox4x(UInt64, (UInt64(ms), UInt64(i)))
holdout_rng(ms, i) = Philox4x(UInt64, (UInt64(ms) ⊻ HOLDOUT_SALT, UInt64(i)))   # disjoint stream
```
`Philox4x` is an `AbstractRNG` → drops straight into `sample_prior`/`simulate_pair`. Use `holdout_rng` to re-simulate the reserved set; the stream is a pure function of `(master_seed, idx)`, thread/order-independent.

### Atomic + integrity-checked JLD2 writes
**Source:** `spike/data/cache.jl:97-106` (`write_shard`), `generate.jl:191-196` (holdout)
**Apply to:** `run_advi.jl` (`advi_artifact.jld2`), any trained-model/benchmark-result persistence
```julia
jldsave(tmp; …, schema_version = SCHEMA_VERSION)
JLD2.jldopen(tmp, "r") do f; @assert haskey(f, "<key>") "integrity check failed: $tmp"; end
mv(tmp, path; force = true)    # filesystem-atomic commit
```
Schema-version + meta tag every artifact; integrity-check on reopen (the cross-env JLD2 hand-off is the ONLY coupling between the broken-Turing baseline and the pinned spike — D-02).

### CPU-only / `use_gpu=false` everywhere
**Source:** `spike/00_smoke.jl:65-67`, `runtests.jl:52-58`
**Apply to:** every `train`/`sampleposterior`/`assess` call + the runtests CUDA gate
`use_gpu` defaults to **true**; pass `use_gpu=false` on every call. The `runtests.jl` CPU-only testset (no CUDA direct dep / resolved node / loaded module) must stay green. GPU (D-10) is a separate opt-in branch guarded by a device check, never a forced import.

## No Analog Found

| File | Role | Data Flow | Reason |
|------|------|-----------|--------|
| (none fully analog-less) | — | — | `benchmark.jl` is the weakest match: no existing BenchmarkTools usage in-repo, so its `@belapsed` harness follows RESEARCH "Code Examples" rather than a repo file; its threading/re-sim halves DO have analogs (`generate.jl`). |

## Conventions

> Deterministic convention-derivation tool skipped (`reason: no-readable-files`) — `gsd-tools verify conventions` is JS/TS-oriented and does not parse Julia (`.jl`) sources. Conventions below are observed manually across the 9 spike `.jl` files read.

| Axis | Dominant | Share | Entropy | Status |
|------|----------|-------|---------|--------|
| File-name casing | `snake_case.jl` (`loader.jl`, `test_data_pipeline.jl`, `00_smoke.jl`) | 100% | low | named contract |
| Identifier casing | `snake_case` funcs (`load_fold`, `sample_prior`); `SCREAMING_SNAKE` consts (`AUG_DIM`, `FOLD_SALT`); Unicode math (`θ`, `μ`, `ρ_true`, `Ztr`) | ~95% | low | named contract |
| Module/export style | single `module Xxx` + `export …` + trailing `using .Xxx` re-export (`loader.jl:51-203`); flat files use guarded top-level `include` | mixed by file role | medium | named contract per role |
| Include style | `include(joinpath(@__DIR__, "..", …))` always; guarded by `isdefined`/`isfile` in shared files | 100% | low | named contract |

**Contested hotspots (author's choice):** none material in the spike subtree — `spike/**` is uniformly snake_case + `@__DIR__`-relative guarded includes + AGPL license header + `const`-fixtures-before-gates. The one **intentional cross-directory split** to respect is the **spike ↔ baseline ↔ src env boundary**: `spike/**` is the pinned-NeuralEstimators env (snake_case spike modules, read-only into `src/`), `spike/baseline/**` is a separate isolated Turing env (own Project/Manifest), and `src/**` is the untouched parent package (PascalCase types like `MultiChannelImage`, `CoLocResult`; `colocalization`/`compute_BayesFactor` lower-case funcs). Each half is internally consistent; match the directory's local style and NEVER edit `src/` from the spike (read-only `include()` only). Every new `.jl` file MUST carry the GNU-AGPL `#= … =#` header block (all 9 spike files have it).

## Metadata

**Analog search scope:** `spike/` (data, simulator, test, root), `src/bayes.jl` (read-only target)
**Files scanned:** 12 read in full/targeted (`00_smoke.jl`, `runtests.jl`, `test_data_pipeline.jl`, `loader.jl`, `encode.jl`, `ghat.jl`, `prior.jl`, `seeding.jl`, `generate.jl`, `Project.toml`, `src/bayes.jl` §§, plus `cache.jl`/`forward.jl`/`contract.jl` signatures via grep)
**Pattern extraction date:** 2026-06-30
