# Phase 3: Training-Data Pipeline - Research

**Researched:** 2026-06-27
**Domain:** Reproducible large-scale simulation-based training-data generation in Julia (sharded JLD2 cache, counter-based RNG, leak-free k-fold loader) for NeuralEstimators NPE
**Confidence:** HIGH (architecture/patterns), MEDIUM (NeuralEstimators tensor layout — Phase-4 owns final reshape)

## Summary

Phase 3 builds a **generator → sharded cache → leak-free loader** pipeline entirely under
`spike/` (hard decoupling: `src/` and root manifests stay byte-identical to baseline `f581d95`;
`src/` is reached only via the read-only `include()` in `spike/contract.jl`). The generator drives
the already-built `sample_prior(rng)` → `simulate_pair(rng, θ; imsize)` → UNCHANGED
`patch_summary` chain, encodes the 8×8 correlation matrix into the locked **128-dim** D-01 vector
(64 imputed values + 64 mask) plus an augmented-moment variant (D-02), and writes RAW (θ, summaries)
to sharded JLD2 chunks under a content-hash-named directory. A separate loader is the *only* path to
standardized tensors; it fits `StatsBase.ZScoreTransform` on TRAIN folds only — making cross-fold
leakage impossible by construction (D-07/D-08).

The four young-stack primitives were verified live: **Random123** `Philox4x`/`Threefry4x` are
`AbstractRNG`s with per-stream keying + `set_counter!` (slots into the existing `rng::AbstractRNG`
simulator interface with zero simulator changes); **NeuralEstimators** expects parameters as a
**d×K matrix** (parameters in columns, datasets in columns for Z) — so the cache should store
column-major (one sample = one column); **JLD2** is pure-Julia (no HDF5 binary dep, Windows-clean);
**SHA** (stdlib) supplies `sha256` for the D-05 content hash. **JLD2 and Random123 are NOT yet in
the spike environment** — they must be added (extends `spike/Project.toml` + re-freeze Manifest,
exactly the Phase-2 Wave-0 pattern).

**Primary recommendation:** One generation pass writes RAW column-major (θ, summary_min[128],
summary_aug) into atomically-written `shard_NNNN.jld2` files inside `spike/data/cache/<sha256-hash>/`;
resume = skip shards whose final file exists and passes an integrity check; a SHA-256 hash over the
simulator/prior/ghat/contract **source bytes + canonical config** names the directory and
auto-invalidates stale data; per-sample RNG is `Philox4x(UInt64, (master_seed, global_index))`
(keyed by index → identical regardless of thread count); the loader is the sole standardization path.

## User Constraints (from CONTEXT.md)

### Locked Decisions

- **D-01:** 128-dim summary = 64 correlation values (missing→0 imputed) + 64-dim binary present/absent
  mask, distinguishing genuine-zero from missing patches. 8×8 patch grid stays FIXED. Mask convention
  must be documented so Phase-4's summary-net input dimension is unambiguous.
- **D-02:** Cache BOTH summary variants in ONE pass — the minimal D-01 summary AND the augmented-moment
  summary (+Manders/median/IQR moments) for Phase-4's ABL-01 ablation, so Phase 4 runs zero
  re-simulation. Exact augmented-moment set = researcher/planner recommendation; the commitment to
  cache both is locked.
- **D-03/D-09:** `imsize` is a per-sample sampled nuisance over a range of microscopy sizes
  (e.g. 256²–2048²), not fixed. Trains a size-robust net. Cost implication is real (large sizes
  dominate sim time) — budget/parallelism must absorb it. Exact size set / sampling distribution =
  Claude's discretion.
- **D-04:** Sharded JLD2 chunks (~10k samples/shard, BayesInteractomics pattern). Resume = skip
  completed shards, restart only the in-flight one. Exact shard size = Claude's discretion.
- **D-05:** Content hash over simulator/prior/ghat SOURCE + generating CONFIG (imsize-range, summary
  spec incl. D-01 mask, θ-dim, master seed). Any code OR config change auto-invalidates stale data.
  Config-only / hand-bumped-integer guard rejected as too weak. Exact hash function + which source
  files feed it = Claude's discretion.
- **D-06:** Default budget 50k; `N` is a single config knob; adding shards extends the dataset without
  regenerating existing ones; sharded scale-up to 200k.
- **D-07:** Cache stores RAW summaries + θ (standardization-free). Loader fits z-scoring mean/std on
  TRAIN folds only, applies to held-out fold. Storing globally-standardized summaries or a global-stats
  sidecar is REJECTED (bakes global stats into every fold = the prohibited leakage).
- **D-08:** The loader owns the split; there is NO public global-standardize path. Single loader API
  returns per-fold standardized tensors; cross-fold leakage impossible by construction. Convention +
  regression-test approach rejected (leaves the footgun in place).
- **D-09 (fold):** k=5 fold CV; fold assignment deterministic from master seed; k configurable.
- **D-10:** Reserved ADVI-benchmark holdout ≥20 stacks from its OWN disjoint Random123 counter range,
  excluded from every CV fold by construction, sampled i.i.d. from π(θ). Not carved post-hoc; not a
  stratified grid.
- **D-11:** Random123 counter-based per-sample seeding `f(master_seed, global_sample_index)` —
  independent and reproducible regardless of execution order or thread count; reserved set occupies a
  disjoint counter range. Simulator already accepts explicit `rng::AbstractRNG` (Phase-2 D-14).
- **D-12:** `Threads.@threads` multithreaded generation, identical results regardless of thread count,
  serial fallback on a single thread. Distributed.jl reserved.
- **D-13:** θ sampled i.i.d. from calibrated prior π(θ). Negative-correlation tail sparsity (Phase-2
  D-16) documented as a known caveat, NOT engineered away by stratification.

### Claude's Discretion

- Exact shard size (D-04) and imsize set / sampling distribution (D-03/D-09).
- Exact augmented-moment set for the second summary variant (D-02).
- Internal module/file layout within `spike/` (e.g. `spike/data/`) and the precise JLD2 key schema.
- Exact hash function and which source files feed the D-05 content hash.
- Mask dtype / 128-dim layout in the tensor handed to NeuralEstimators.

### Deferred Ideas (OUT OF SCOPE)

- NPE/NRE training, the ADVI benchmark, the summary-statistic ablation → Phase 4 (Phase 3 only
  *provisions* the augmented summary D-02 and reserved holdout D-10).
- SBC / amortized Bayes factor / OOD flag → Phase 5.
- Distributed.jl multiprocess generation → reserved (D-12) unless single-node threading is too slow.
- Fixing the negative-correlation tail sparsity (stratification / importance weighting) → explicitly
  NOT done (D-13).

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| DATA-01 | Generator maps θ~π → `simulate_pair` → existing summary functions → fixed-dim standardized summary vector | §The 128-Dim D-01 Encoding + §Architecture Patterns (generator); standardization deferred to the loader (D-07) |
| DATA-02 | 50k–200k pairs, version-guarded JLD2 cache (BayesInteractomics pattern), Random123-seeded loader; changed simulator/prior auto-invalidates | §Sharded JLD2 Cache + Content-Hash Guard + §Random123 Counter-Based Seeding |
| DATA-03 | Leak-free k-fold CV (fit-on-train-only standardization in the loader) + reserved ≥20-stack ADVI holdout | §Leak-Free k-Fold Loader; §Reserved Holdout via Disjoint Key Namespace |

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| θ sampling from π(θ) | Existing `spike/simulator/prior.jl` (`sample_prior`) | — | Already built & calibrated (Phase 2); generator only drives it |
| Per-pair physics + summary | Existing `simulate_pair` + `contract.jl` (UNCHANGED) | `src/` via read-only `include()` | Phase-2 outputs; Phase 3 must not edit them |
| 128-dim D-01 / augmented encoding | NEW `spike/data/` encode layer | reads `patch_summary` matrix | Pure transform of the 8×8 matrix; no src edit |
| RNG / reproducibility | NEW seeding layer (Random123) | simulator's `rng::AbstractRNG` slot | Counter-based keys make threads order-independent |
| Cache persistence + resume + hash guard | NEW `spike/data/` cache layer (JLD2 + SHA) | filesystem | Sharded, atomic, content-hash-invalidated |
| Split / standardization | NEW `spike/data/` loader (StatsBase) | — | Sole standardization path; leak-free by construction |
| Validation gates | NEW `spike/test/` testset | stdlib `Test` | Re-runnable quantitative gates (Nyquist) |

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| JLD2 | latest (pure-Julia, no HDF5 C dep) | Sharded cache of (θ, summaries); model/meta serialization | CLAUDE.md-mandated; Windows-clean (no binary HDF5); de-facto Julia serialization; the cited BayesInteractomics pattern `[CITED: CLAUDE.md]` |
| Random123 | latest | Counter-based per-sample reproducible RNG (Philox4x/Threefry4x) | CLAUDE.md-mandated; stateless/splittable; `AbstractRNG` → drops into `simulate_pair(rng,…)` `[VERIFIED: RandomNumbers.jl docs]` |
| StatsBase | already in spike env | `ZScoreTransform` fit-on-train standardization; moments for D-02 | CLAUDE.md-mandated; `fit(ZScoreTransform, X; dims=2)` / `transform` is the standard z-scoring API `[CITED: StatsBase docs]` |
| SHA | stdlib (no add) | `sha256` over source+config for the D-05 content hash | Stdlib, deterministic, available (verified `using SHA` works) `[VERIFIED: julia -e]` |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Statistics | stdlib | `mean`/`median`/`quantile`/`std` for moments | D-01 imputation, D-02 augmented moments |
| Base.Threads | stdlib | `Threads.@threads` over shards (D-12) | Generation parallelism; serial when `nthreads()==1` |
| NeuralEstimators | 0.2.1 (pinned) | Consumer of the loader's tensors (Phase 4) | Phase 3 only matches its **d×K column-major** layout; no training here |

**Packages to ADD to the spike env (NOT currently present — verified via `Pkg.project().dependencies`):**
`JLD2`, `Random123`. (`StatsBase`, `NeuralEstimators` already present.) This is a Wave-0 task:
extend `spike/Project.toml`, `Pkg.resolve()`, re-run the Phase-2 resolve-risk gate
(`runtests.jl` asserts NeuralEstimators stays pinned at v0.2.1 — verify JLD2/Random123 do not
downgrade it), and re-freeze `spike/Manifest.toml`.

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| JLD2 | BSON / HDF5.jl | HDF5 adds a binary C dep (Windows friction); BSON is slower/less typed. JLD2 is the locked choice. |
| Random123 key-per-sample | `set_counter!` stride per sample | Keying by `(master_seed, index)` gives a fresh independent stream per sample with zero collision risk; `set_counter!` stride works but needs a guaranteed-non-overlapping stride budget. Recommend keying. |
| `Threads.@threads` | Distributed.jl | Multiprocess reserved (D-12) until single-node threading proves too slow. |

**Installation:**
```bash
# In ProteinCoLoc/spike/ — own environment, NEVER touch root Project.toml/Manifest.toml
julia --project=spike -e 'using Pkg; Pkg.add(["JLD2","Random123"]); Pkg.resolve()'
# Then re-freeze the Manifest and re-run the resolve-risk gate (NeuralEstimators must stay v0.2.1)
```

**Version verification:** JLD2 + Random123 are both in the General registry (registry update succeeded).
Pin exact resolved versions into `spike/Manifest.toml` after resolve, as Phases 1–2 did. Confirm the
resolve does NOT downgrade NeuralEstimators 0.2.1 (the Phase-1 GLMakie-style co-resolve regression risk;
`runtests.jl` already gates this by UUID).

## Package Legitimacy Audit

| Package | Registry | Age | Source Repo | slopcheck | Disposition |
|---------|----------|-----|-------------|-----------|-------------|
| JLD2 | Julia General | ~8 yrs | github.com/JuliaIO/JLD2.jl | n/a (Julia) | Approved — CLAUDE.md-mandated, pure-Julia |
| Random123 | Julia General | ~8 yrs | github.com/JuliaRandom/Random123.jl | n/a (Julia) | Approved — CLAUDE.md-mandated |
| SHA | stdlib | bundled | JuliaLang stdlib | n/a | Approved — stdlib |
| StatsBase | Julia General | already resolved | github.com/JuliaStats/StatsBase.jl | n/a | Already in spike env |

slopcheck targets the npm/PyPI hallucination vector and does not apply to the Julia General registry
(curated, no arbitrary publish). All four packages are long-established, CLAUDE.md-named or stdlib, with
known source repos — no slop risk. Legitimacy established by CLAUDE.md mandate + registry resolve.

## Architecture Patterns

### System Architecture Diagram

```
                       master_seed, N, k, imsize-spec, summary-spec  (CONFIG)
                                          │
                                          ▼
            ┌──────────────── D-05 CONTENT HASH ────────────────┐
            │ sha256( source_bytes[forward,prior,ghat,contract] │
            │         + canonical(CONFIG) )  → <hash>           │
            └───────────────────────┬───────────────────────────┘
                                    │ names the cache dir
                                    ▼
   global_index ─► Philox4x(UInt64,(master_seed, idx))  ──► rng_i (independent stream)
                                    │
        ┌─────────────── Threads.@threads over SHARDS ───────────────┐
        │  for each sample idx in shard:                              │
        │    θ      = sample_prior(rng_i)        (UNCHANGED prior.jl) │
        │    imsize = sample_imsize(rng_i)       (D-03 nuisance)      │
        │    imgs   = simulate_pair(rng_i, θ; imsize)  (UNCHANGED)    │
        │    M8x8   = patch_summary(build_mci(imgs))   (UNCHANGED)    │
        │    s_min  = encode_d01(M8x8)   → 128-dim (64 vals + 64 mask)│
        │    s_aug  = encode_aug(imgs,M8x8) → +moments (Manders/...)  │
        └──────────────┬─────────────────────────────────────────────┘
                       ▼  column-major: one sample = one column
        jldsave(shard_NNNN.jld2.tmp; theta, summary_min, summary_aug,
                global_index, imsize) ──atomic mv──► shard_NNNN.jld2
                       │
   RESERVED HOLDOUT ───┤  (disjoint key namespace: Philox4x(UInt64,(master_seed ⊻ SALT, idx)),
   ≥20 i.i.d. π(θ)     │   written to holdout_*.jld2 — NEVER in the main pool)
                       ▼
   meta.jld2 { hash, per-file subhashes, config, N, shard_size, schema_version }
                       │
                       ▼
        ┌──────────── LOADER (sole standardization path) ────────────┐
        │ load_fold(cache, k; K=5, master_seed):                     │
        │   perm   = randperm(seeded_rng, N)   (deterministic folds) │
        │   train,val = split(perm, k)                               │
        │   zt   = fit(ZScoreTransform, Z[:,train]; dims=2)  ◄ TRAIN  │
        │   Ztr  = transform(zt, Z[:,train]) ; Zva = transform(zt, Z[:,val]) │
        │   return (Ztr, θtr, Zva, θva, zt)   NO global-standardize API│
        └────────────────────────────────────────────────────────────┘
                       │
                       ▼  d×K matrices (Float32) → Phase-4 NeuralEstimators train()
```

### Recommended Project Structure
```
spike/
├── data/
│   ├── generate.jl     # generator: prior→simulate→encode→shard write; @threads; resume
│   ├── encode.jl       # encode_d01 (128-dim) + encode_aug (moments) — pure transforms
│   ├── seeding.jl      # Random123 keyed RNG per (master_seed, global_index); holdout namespace
│   ├── cache.jl        # JLD2 shard read/write, atomic write, resume-by-skip, meta.jld2
│   ├── hashguard.jl    # D-05 sha256 over source files + canonical config
│   └── loader.jl       # leak-free k-fold loader (StatsBase fit-on-train); NO global-standardize
├── data/cache/<hash>/  # generated artifacts (gitignored except a tiny smoke fixture)
│   ├── shard_0001.jld2 ...
│   ├── holdout.jld2
│   └── meta.jld2
└── test/
    └── test_data_pipeline.jl   # Validation Architecture gates (included by runtests.jl)
```

### Pattern 1: Column-major cache matching NeuralEstimators d×K layout
**What:** Store θ as a `d×n` matrix and each summary as a `feature×n` matrix per shard (one sample =
one column). **Why:** NeuralEstimators expects parameters as a **d×K matrix with parameters in columns**
`[VERIFIED: NeuralEstimators docs — "a d×B matrix θ … parameters in columns"]`, and `ZScoreTransform`
with `dims=2` standardizes each feature-row across sample-columns — both align with column-major. The
loader then `hcat`s the needed shards and slices columns by fold with no transpose.
```julia
# Source: pattern derived from NeuralEstimators d×K convention + StatsBase dims=2
# per shard, in-memory then atomic write:
theta        = Matrix{Float64}(undef, 7, n)          # θ columns (D-04 7-vector)
summary_min  = Matrix{Float64}(undef, 128, n)        # D-01: rows 1:64 vals, 65:128 mask
summary_aug  = Matrix{Float64}(undef, A, n)          # D-02 augmented (A = 128 + #moments)
global_index = Vector{Int}(undef, n)
imsize       = Vector{Tuple{Int,Int}}(undef, n)
```

### Pattern 2: Atomic shard write + resume-by-skip (D-04)
**What:** Build a whole shard in memory, `jldsave` to `shard_NNNN.jld2.tmp`, then `mv` to the final
name. A shard is "done" iff its final file exists and loads (integrity check). Resume skips done shards
and regenerates only the in-flight one (its `.tmp` is discarded). **Why:** A partially-written final
file would silently corrupt the dataset; the tmp+rename makes completion atomic on the filesystem.
```julia
# Source: standard atomic-write idiom; JLD2 jldsave one-shot per shard
tmp = path * ".tmp"
jldsave(tmp; theta, summary_min, summary_aug, global_index, imsize, schema_version=1)
# integrity check by reopening before commit
JLD2.jldopen(tmp, "r") do f; @assert haskey(f, "theta"); end
mv(tmp, path; force=true)            # atomic commit
```
**Resume check:** `isfile(path) && _loads_ok(path)` → skip. One-shot `jldsave` per shard is recommended
over in-file append: a 10k-sample shard is ~10–20 MB in memory (128+A+7 floats × 10k), trivially fits,
and avoids JLD2 append/group complexity.

### Pattern 3: Content-hash version guard (D-05)
**What:** `sha256` over the concatenated **source bytes** of the files that define the data semantics,
plus a **canonical serialization of the config**. **Which files feed it (recommendation):**
`spike/simulator/forward.jl`, `spike/simulator/prior.jl`, `spike/simulator/ghat.jl`,
`spike/contract.jl` (the encode logic file `spike/data/encode.jl` too, once it exists). Optionally the
two frozen `src/` files (`colocalization.jl`, `LoadImages.jl`) as a defensive check — they are
baseline-frozen but *do* define the summary, so hashing them catches an accidental decoupling breach.
```julia
# Source: SHA stdlib (verified available)
using SHA
function cache_hash(config::NamedTuple; src_files)
    h = SHA.SHA256_CTX()
    for f in src_files                      # FIXED order
        SHA.update!(h, read(f)); SHA.update!(h, UInt8['\0'])
    end
    SHA.update!(h, Vector{UInt8}(canonical(config)))  # sorted-key deterministic string
    return bytes2hex(SHA.digest!(h))
end
```
Store the full hash + **per-file sub-hashes** + the raw config in `meta.jld2`. On load: recompute,
compare; mismatch → refuse to reuse stale data (error with a diff of which file/config key changed, or
regenerate). The per-file sub-hashes make "which source changed" diagnosable — strictly better than a
single opaque hash. **Config fields hashed:** N, shard_size, imsize-spec (set + sampling weights),
summary-spec (D-01 mask convention + augmented-moment list), θ-dim (7), master_seed, k.

### Pattern 4: Random123 per-sample keyed RNG (D-11)
**What:** Each sample constructs its own RNG **keyed by its global index**, so the draw stream is a pure
function of `(master_seed, global_index)` — independent of thread count and execution order.
```julia
# Source: RandomNumbers.jl Random123 docs (VERIFIED constructors + AbstractRNG)
using Random123
# independent stream per sample, keyed by index (no counter-stride collision risk):
rng_i = Philox4x(UInt64, (UInt64(master_seed), UInt64(global_index)))
θ = sample_prior(rng_i)                      # AbstractRNG → drops into existing simulator
imgs = simulate_pair(rng_i, θ; imsize=isz)   # same rng threaded through (Phase-2 D-14 interface)
```
`Philox4x`/`Threefry4x` are `AbstractRNG`s usable with `rand(rng, …)` `[VERIFIED]`; the simulator's
`rng::AbstractRNG` slot (forward.jl L110, prior.jl L71) accepts them with NO simulator change.
**Why key-per-sample beats one shared RNG + `set_counter!`:** under `@threads`, a shared RNG is mutable
shared state (data race); a fresh per-sample keyed RNG is race-free and order-independent by
construction. `set_counter!(rng, idx)` is the documented alternative if a single keyed generator with a
reserved stride is preferred, but keying is simpler and collision-free.

### Pattern 5: Reserved holdout via disjoint key namespace (D-10)
**What:** Draw the ≥20-stack ADVI holdout from a **provably disjoint** key space so it can never appear
in the main pool or any CV fold.
```julia
const HOLDOUT_SALT = 0xA11D1F0_holdout_const  # any fixed nonzero UInt64
# main pool:   Philox4x(UInt64, (master_seed,           idx))   idx ∈ 1:N
# holdout:     Philox4x(UInt64, (master_seed ⊻ HOLDOUT_SALT, idx))  idx ∈ 1:H,  H ≥ 20
```
Because the holdout uses a different *key* (first key word XOR-salted) AND lives in separate
`holdout.jld2` files never indexed into the main pool, exclusion from every fold is structural — the
loader only ever iterates main-pool indices `1:N`. Holdout θ is i.i.d. from π(θ) (same `sample_prior`),
satisfying D-10's "prior-representative, not a stratified grid."

### Pattern 6: Leak-free k-fold loader (D-07/D-08) — the structural footgun-removal
**What:** A single `load_fold` is the **only** function that returns standardized tensors. There is NO
exported `standardize_all` / global-stats function, so no caller *can* leak.
```julia
# Source: StatsBase ZScoreTransform fit/transform (dims=2 = per-feature across samples)
using StatsBase
function load_fold(cache, fold::Int; K::Int=5, master_seed)
    Z, θ = load_main_pool(cache)                 # RAW d×N / feature×N, columns = samples
    perm = randperm(fold_rng(master_seed), size(Z,2))   # deterministic from master seed
    val_idx   = perm[fold:K:end]                 # or contiguous block of perm
    train_idx = setdiff(1:size(Z,2), val_idx)
    zt = fit(ZScoreTransform, Z[:, train_idx]; dims=2)   # ◄ FIT ON TRAIN ONLY
    Ztr = StatsBase.transform(zt, Z[:, train_idx])
    Zva = StatsBase.transform(zt, Z[:, val_idx])         # applied to val, never fit on it
    return (Ztr=Float32.(Ztr), θtr=θ[:,train_idx],
            Zva=Float32.(Zva), θva=θ[:,val_idx], zt)
end
```
`fold_rng` = a Random123 generator keyed by `(master_seed ⊻ FOLD_SALT, 0)` so folds are reproducible
and independent of the generation key stream. The leak-impossibility argument: the cache is RAW (no
global stats stored, D-07); the only standardization entry point requires choosing a fold and fits
inside on that fold's train columns; `zt` is returned so Phase-5 can FREEZE train-only preprocessing
(ROADMAP Phase-5 SC-5). **The mask channel (rows 65:128) should NOT be z-scored** — it is binary
indicator data; recommend standardizing only the value rows (1:64 of D-01, and the continuous moment
rows of D-02) and passing the mask through unchanged. Document this in the loader (z-scoring a 0/1 mask
distorts it and couples folds through the mask mean). This is a concrete planning detail.

### Anti-Patterns to Avoid
- **Storing globally-standardized summaries or a global-stats sidecar** (D-07 REJECT): bakes global
  statistics into every fold = the exact prohibited leakage.
- **Exposing any `standardize_all`/global-mean-std function** (D-08 REJECT): leaves the footgun in
  place; convention+test is weaker than structural impossibility.
- **Carving the holdout post-hoc from the main pool** (D-10 REJECT): risks fold overlap; use the
  disjoint key namespace.
- **Sharing one RNG across `@threads` iterations**: data race + thread-count-dependent results; use
  per-sample keyed RNGs.
- **In-file JLD2 append per sample**: needless complexity + partial-write risk; build a shard in
  memory and write once atomically.
- **z-scoring the 64-dim binary mask**: distorts the present/absent signal and re-introduces a
  cross-fold coupling through the mask mean.
- **Stratifying θ to fix the negative-tail sparsity** (D-13 REJECT): distorts π(θ), breaks SBC.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| z-score standardization | manual mean/std loops | `StatsBase.fit(ZScoreTransform; dims=2)` + `transform` | Correct per-feature dims handling; the fit object is reusable + freezable for Phase 5 |
| Reproducible per-stream RNG | hashing seeds into `MersenneTwister` | `Random123.Philox4x`/`Threefry4x` keyed by index | Counter-based generators are designed for exactly "independent stream per index", high quality, splittable `[VERIFIED]` |
| Content hashing | custom checksum | `SHA.sha256` (stdlib) | Cryptographic, deterministic, collision-resistant; zero deps |
| Typed binary cache | CSV / custom binary | JLD2 `jldsave`/`jldopen` | Pure-Julia, preserves array types/shapes, Windows-clean |
| Atomic file commit | write-in-place | write `.tmp` + `mv` | Filesystem rename is atomic; prevents half-written shards on crash/resume |
| Deterministic fold split | ad-hoc index math | `randperm(seeded_rng, N)` + slice | One seeded permutation gives reproducible, balanced, leak-controlled folds |

**Key insight:** Every primitive this phase needs (z-scoring, counter RNG, hashing, typed serialization,
atomic write) has a standard, battle-tested implementation. The phase's *value* is the **structural
composition** (leak-impossible loader, hash-invalidated cache, order-independent parallelism), not any
custom algorithm.

## Runtime State Inventory

> Rename/refactor categories — Phase 3 is greenfield *generation* but writes persistent artifacts, so
> the cache-state questions apply.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | The JLD2 cache under `spike/data/cache/<hash>/` is NEW persistent state keyed by content hash. | Gitignore the bulk cache (large); the hash-named dir auto-separates incompatible versions — no manual migration. Keep a tiny committed smoke fixture for tests. |
| Live service config | None — pure local file generation, no external service. | None. |
| OS-registered state | None — no scheduled tasks/daemons; generation is a script run under `julia -t`. | Document `JULIA_NUM_THREADS`/`-t` requirement (default `nthreads()==1` verified). |
| Secrets/env vars | `JULIA_NUM_THREADS` (or `-t auto`) controls D-12 parallelism; `master_seed` is a config value, not a secret. | Document in NOTES; thread count must NOT change results (gate it). |
| Build artifacts | `spike/Manifest.toml` changes when JLD2/Random123 are added; `spike/data/cache/` artifacts. | Re-freeze Manifest (Wave-0); gitignore cache bulk. Verify NeuralEstimators stays 0.2.1. |

**Decoupling verified:** the generator touches ONLY `spike/`; `src/` reached read-only via
`contract.jl`'s `include()`. No `src/`/root-manifest edits — baseline `f581d95` stays byte-identical
(the standing DEMO-02 gate).

## The 128-Dim D-01 Encoding (DATA-01)

`patch_summary(mci)` returns an `8×8 Matrix{Union{Float64,Missing}}` (a patch is `missing` when
`_exclude_zero` leaves ≤15 surviving pixels — verified `src/colocalization.jl:235`). Encode:

```julia
# Source: D-01 spec + Julia column-major vec()
function encode_d01(M::AbstractMatrix)   # M is 8×8 Union{Float64,Missing}
    vals = vec(coalesce.(M, 0.0))                    # 64-dim, missing→0  (column-major order)
    mask = vec(Float64.(.!ismissing.(M)))            # 64-dim, 1=present / 0=missing
    return vcat(vals, mask)                           # 128-dim; rows 1:64 vals, 65:128 mask
end
```
- **Layout (document for Phase 4):** rows 1:64 = imputed correlations in **column-major** vec order of
  the 8×8 grid; rows 65:128 = the parallel binary mask, same ordering. Phase-4's summary net input dim
  is unambiguously **128**.
- **Mask dtype:** store RAW as `Float64` in the cache; convert the whole tensor to `Float32` at the
  loader boundary (Flux/NeuralEstimators default eltype is Float32 for speed). The mask rows must be
  passed through standardization UNCHANGED (see Pattern 6).
- **Edge case:** a fully-missing 8×8 (all patches below the 15-px floor) → vals all 0, mask all 0. This
  is rare given the `BG_FLOOR` background (forward.jl L66 keeps pixels >0) but possible at the smallest
  imsize; the mask makes it self-describing to the net. `induced_mu` already returns `NaN` here
  (contract.jl L104) — the generator should still cache the sample (the mask encodes the degeneracy);
  do NOT drop it silently (that would distort the π(θ)-faithful distribution, D-13).

## The Augmented-Moment Variant (D-02) — Recommended Set

Cache `summary_aug` as a **superset**: the 128 D-01 dims **plus** per-pair scalar moments, so Phase-4's
ABL-01 slices columns with zero re-simulation. All computable from the two channel matrices + the 8×8
matrix WITHOUT editing `src/` (spike-local helpers in `encode.jl`).

| Moment | Definition | Source |
|--------|-----------|--------|
| Manders M1, M2 | M1 = Σ ch1ᵢ·[ch2ᵢ>t₂] / Σ ch1ᵢ ; M2 symmetric. Thresholds t₁,t₂ = `mci.otsu_threshold` (already on the MCI). | `[CITED: Manders et al. 1993; thresholds from MultiChannelImage]` |
| Whole-image Pearson | single `cor(vec(ch1), vec(ch2))` over non-excluded pixels | overall coloc scalar |
| Patch-grid median | `median(skipmissing(M))` | distribution location |
| Patch-grid IQR | `quantile(v,0.75) - quantile(v,0.25)` of `skipmissing(M)` | distribution spread |
| Patch-grid mean (=induced_μ) | `mean(skipmissing(M))` | matches the Turing μ summary |
| Patch-grid std | `std(skipmissing(M))` | ↔ σ/τ scale (D-03) |
| Patch-grid skewness, excess kurtosis | `StatsBase.skewness/kurtosis(skipmissing(M))` | ↔ tail/ν |
| Fraction-missing | `count(ismissing,M)/64` | sparsity / signal-coverage |
| Per-channel intensity median + IQR (×2) | on `vec(ch1)`, `vec(ch2)` | intensity-regime nuisance descriptors |

Recommended `A = 128 + ~12` (round to a documented constant). Exact composition is Claude's discretion
(D-02); the **commitment to cache both variants** is locked. **Tag:** the specific moment list is an
`[ASSUMED]`/recommendation — the planner/Phase-4 may trim it; only "cache D-01 + an augmented superset"
is locked. Manders requires the otsu thresholds already stored on `MultiChannelImage` (no src edit).

## imsize Distribution (D-03/D-09) — Recommendation + Cost

**Recommended discrete set:** `{256², 512², 1024², 1376×1028 (real anchor, Phase-2 D-08), 2048²}` — all
≥64 (forward.jl validation) and 8-divisible-friendly (1376÷8=172 exact, 1028÷8=128 trims 4 px,
verified). **Cost model:** `simulate_pair` cost is dominated by per-pixel work (3 smooth-field
`imfilter`s + Bernoulli + PSF + Poisson over W·H) → **cost ∝ W·H**. Relative per-sample cost vs 256²:

| imsize | pixels | rel. cost |
|--------|--------|-----------|
| 256² | 65k | 1× |
| 512² | 262k | 4× |
| 1024² | 1.05M | 16× |
| 1376×1028 | 1.41M | ~21.6× |
| 2048² | 4.19M | 64× |

A **uniform** draw over these sizes is dominated by the 2048² tail (mean ≈ 21× the 256² cost).
**Recommendation:** cost-aware categorical weights favoring smaller sizes (e.g. weights ∝ 1/pixels,
giving roughly 0.62/0.15/0.04/… normalized — heavily 256²/512²) OR cap the large-size fraction (e.g.
≤10% at ≥1024²) so the 50k-default budget stays tractable on CPU. The planner should pick concrete
weights and record the **expected per-sample cost** = Σ wᵢ·costᵢ to size the generation budget. This is
the main lever balancing D-03 size-robustness against D-06's 50k→200k budget.

## Threads.@threads Generation (D-12) — Memory & Reproducibility

- **Granularity:** `Threads.@threads` over **shards** (coarse). Each thread builds one shard into
  thread-local arrays, then atomically writes it — no shared mutable state, no locks. Sample-level
  `@threads` within a shard also works but shard-level bounds memory more cleanly.
- **Order-independence (D-12 core claim):** because each sample's RNG is keyed by its *global index*
  (Pattern 4), the `(θ, imgs, summary)` for index `i` is identical no matter which thread computes it
  or in what order → byte-identical dataset for `nthreads ∈ {1, 2, …}`. Serial fallback = the same loop
  with `nthreads()==1` (verified default on this machine).
- **Memory bound (the real driver):** peak ≈ `T × (largest in-flight simulate_pair working set)`. A
  2048²×2-channel pair plus the 3 smooth fields ≈ 4.19M × 8 B × ~5 arrays ≈ **~170 MB per in-flight
  sample**; with `T` threads each possibly on a 2048² sample, peak ≈ `T × ~170 MB`. The shard buffer
  itself is small (~10–20 MB). **Recommendation:** the imsize weighting (above) is also the memory
  governor — capping large-size frequency bounds simultaneous large allocations. Document the
  `JULIA_NUM_THREADS`/`-t` launch requirement (default is 1; generation must be launched with `-t auto`
  or an explicit count).

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|------------------|--------|
| `Random.seed!` global state for parallel sims | Random123 counter/key-based per-stream RNG | Thread-count-independent reproducibility (D-11/D-12) — the standard SBI-at-scale pattern |
| Monolithic single-file dataset | Sharded JLD2 + resume-by-skip | Bounded memory, parallel writers, incremental scale-up (D-04/D-06) |
| Integer/version-bump cache guard | Content hash over source + config | Auto-invalidation on any code/config change — no silent stale reuse (D-05) |
| Standardize-then-cache | Cache RAW, standardize per-fold in loader | Leak-free CV by construction (D-07/D-08) — the modern ML-hygiene default |

## Code Examples

### Reproducible per-sample generation (verified primitives)
```julia
# Source: Random123 docs (VERIFIED) + existing simulator interface
using Random123
function generate_sample(master_seed::UInt64, idx::Int)
    rng = Philox4x(UInt64, (master_seed, UInt64(idx)))     # independent keyed stream
    θ   = sample_prior(rng)                                  # UNCHANGED prior.jl
    isz = sample_imsize(rng)                                 # D-03 nuisance
    mci = build_mci(simulate_pair(rng, θ; imsize=isz))       # UNCHANGED simulator + contract
    M   = patch_summary(mci)                                 # UNCHANGED 8×8 summary
    return (θ=collect(values(θ)), s_min=encode_d01(M), s_aug=encode_aug(mci, M), isz=isz)
end
```

### Hash-guarded cache load
```julia
# Source: SHA stdlib (verified) + JLD2
function open_or_invalidate(cache_root, config; src_files)
    h   = cache_hash(config; src_files)
    dir = joinpath(cache_root, h)
    if isdir(dir)
        meta = JLD2.load(joinpath(dir, "meta.jld2"))
        meta["hash"] == h || error("stale cache: hash mismatch at $dir")  # auto-invalidate
        return dir                                                         # reuse
    end
    return dir   # fresh: generator will create + write meta.jld2
end
```

## Validation Architecture

> `workflow.nyquist_validation: true` (verified in `.planning/config.json`) → section REQUIRED.
> All gates are re-runnable under the spike `Test` harness (`spike/test/test_data_pipeline.jl`,
> included by `spike/test/runtests.jl` exactly as Phase-2's `test_simulator.jl` is).

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Julia stdlib `Test` (`@testset`/`@test`) — the established spike harness |
| Config file | none — `spike/test/runtests.jl` is the entry (Phase-1/2 pattern) |
| Quick run command | `julia --project=spike spike/test/test_data_pipeline.jl` (tiny N≈64 fixture) |
| Full suite command | `julia --project=spike spike/test/runtests.jl` (smoke + simulator + data) |
| Thread-repro run | `julia -t1 … ` vs `julia -t4 …` then assert identical arrays |

### Phase Requirements / Success Criteria → Test Map
| SC / Req | Behavior | Test Type | Automated Gate | File |
|----------|----------|-----------|----------------|------|
| **SC-1 / DATA-01** | generator maps θ~π→simulate→summary→fixed-dim vector | unit | Generate N=64; assert `size(summary_min)==(128,64)`; rows 65:128 ∈ {0,1}; mask=0 ⟺ value row==0; all `isfinite`; θ is the 7-vector | ❌ Wave 0 |
| **SC-2 / DATA-02** | write+reload from hash-guarded JLD2; resume; auto-invalidate | integration | (a) round-trip: write shards, reload, `@test theta == theta_reloaded` & summaries equal; (b) resume: delete 1 shard, re-run, only it regenerates, final dataset identical; (c) invalidation: perturb a tracked source byte (or a config field) → `cache_hash` differs → loader refuses stale | ❌ Wave 0 |
| **SC-3 / DATA-03** | fit-on-train-only standardization, no global path | unit | (a) `@test !( :standardize_all in names(DataLoaderModule) )` (no global-standardize symbol exported); (b) `Ztr` columns have ~0 mean/unit std per row; (c) `Zva` mean/std are NOT exactly 0/1 (proves val not used to fit); (d) mask rows pass through unchanged | ❌ Wave 0 |
| **SC-4 / DATA-03** | k-fold disjoint + reserved ≥20 holdout excluded by construction | unit | (a) `union(folds)==1:N` & pairwise `∩==∅`; (b) `holdout_idx ∩ every fold == ∅`; (c) `length(holdout) ≥ 20`; (d) holdout key namespace ≠ main (XOR-salt) so no index collision; (e) fold assignment reproducible across two loader calls with same master_seed | ❌ Wave 0 |
| **D-11/D-12** | thread-count-independent reproducibility | integration | Generate the same small dataset with `nthreads()==1` and `>1`; `@test` byte-identical θ + summaries (key-per-index ⇒ order-independent) | ❌ Wave 0 |
| **Decoupling (standing)** | `src/` + root manifests untouched | integration | `git diff --quiet f581d95 -- src ':!spike'` style check stays clean (DEMO-02 gate) | reuse existing |

### Sampling Rate
- **Per task commit:** quick gate on a tiny N≈64 fixture (`test_data_pipeline.jl`) — sub-30s.
- **Per wave merge:** full `runtests.jl` (smoke + simulator + data) including resolve-risk gate
  (NeuralEstimators stays v0.2.1 after adding JLD2/Random123).
- **Phase gate:** full suite green + a real 50k generation run reloads + a thread-repro check before
  `/gsd:verify-work`.

### Wave 0 Gaps
- [ ] `spike/test/test_data_pipeline.jl` — covers SC-1..SC-4 + thread-repro (new).
- [ ] Wire it into `spike/test/runtests.jl` (one `include`, mirroring `test_simulator.jl`).
- [ ] Add `JLD2`, `Random123` to `spike/Project.toml`; `Pkg.resolve()`; re-freeze Manifest; confirm
      NeuralEstimators stays 0.2.1 (extend the existing resolve-risk gate to cover the new deps).
- [ ] A tiny committed cache fixture (a few samples) for the round-trip/invalidation tests, with the
      bulk cache gitignored.

## Security Domain

> `security_enforcement` key absent in config → treated as enabled. This is a local data-generation
> pipeline (no network, no auth, no user input surface) — most ASVS categories N/A.

### Applicable ASVS Categories
| Category | Applies | Standard Control |
|----------|---------|------------------|
| V5 Input Validation | yes (already done) | `simulate_pair` validates all θ fields + imsize at entry (forward.jl L116–134) — the generator must surface, not swallow, these `ArgumentError`s |
| V6 Cryptography | partial | SHA-256 used for *integrity/versioning* (D-05), not secrecy — appropriate; do not repurpose as a security boundary |
| V2/V3/V4/V7+ | no | No auth/session/access-control/network surface in a local generator |

### Known Threat Patterns
| Pattern | STRIDE | Mitigation |
|---------|--------|-----------|
| Stale-cache silent reuse (wrong-data-in) | Tampering | D-05 content hash auto-invalidates; per-file sub-hashes diagnose which source changed |
| Untrusted `.jld2` deserialization (JLD2 reconstructs arbitrary types) | Tampering/EoP | Cache is self-generated & hash-checked; never load third-party `.jld2`. Note JLD2 type reconstruction risk if that ever changes |
| Decoupling breach (accidental `src/` edit) | Tampering | Standing `git diff` vs `f581d95` gate (DEMO-02); optionally hash the two frozen `src/` files into D-05 |
| Half-written shard on crash | Availability/Integrity | Atomic `.tmp`+`mv` write + reload-integrity check on resume |

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Exact augmented-moment set (Manders/median/IQR/skew/kurtosis/fraction-missing/intensity moments) | D-02 variant | LOW — D-02 locks only "cache both"; Phase-4 ablation can trim columns. Manders thresholds assumed = MCI otsu_threshold |
| A2 | imsize set `{256²,512²,1024²,1376×1028,2048²}` + cost-aware weighting | imsize Distribution | MEDIUM — affects budget/runtime; planner picks concrete weights; cost∝W·H is the verified physics |
| A3 | NeuralEstimators consumes a `feature×K` matrix for vector summaries (single replicate/dataset) | Pattern 1 | MEDIUM — θ d×K verified; Z-as-matrix for scalar-summary datasets is the documented quick-start pattern but Phase-4 owns the final reshape (could be a vector-of-arrays). Phase 3 stores column-major RAW either way |
| A4 | BayesInteractomics JLD2 sharded-cache idiom matches the atomic-tmp+mv / one-shot-jldsave pattern | Pattern 2 | LOW — repo not accessible this session; pattern is standard and CLAUDE.md-cited. `[ASSUMED]` |
| A5 | Hashing the 4 spike source files (+ optionally 2 frozen src files) is sufficient for D-05 | Pattern 3 | LOW — covers all data-defining code; planner confirms the file list |
| A6 | Mask rows (65:128) should bypass z-scoring | Pattern 6 / D-01 encoding | LOW-MEDIUM — z-scoring a 0/1 mask is standard-incorrect; but Phase-4 net design may prefer otherwise. Flag for the planner |

## Open Questions (RESOLVED)

1. **Δρ target representation (NPE-01 needs ρ_true AND Δρ).**
   - What we know: cache stores the full 7-vector θ (incl. ρ_true); the loader returns θ columns.
   - What's unclear: Δρ is a *contrast* (difference of correlations / coloc-vs-null) whose definition
     lives in Phase 4. Phase 3 must store enough to derive it.
   - Recommendation: store full θ + summaries (done); leave Δρ derivation to the Phase-4 loader
     consumer. Confirm with the planner that ρ_true + the summary suffice to construct Δρ downstream
     (likely yes — Δρ is computed from paired/contrastive summaries at train time).
   - **RESOLVED:** D-07 (store RAW full 7-vector θ + both summary variants) covers this. The
     generator persists ρ_true within θ; Δρ is a Phase-4 *consumer* concern — Phase 4 constructs
     contrastive coloc-vs-null pairs from the prior and derives Δρ at train time from the stored
     raw summaries. **No Phase-3 schema change needed** (confirmed by the planner: plans 03-02/03-03
     store the full θ-vector and raw summaries, which is sufficient).

2. **Exact NeuralEstimators tensor shape for a fixed-length summary (matrix vs vector-of-arrays).**
   - What we know: θ is d×K (verified); the simple-summary quick-start treats Z column-wise.
   - What's unclear: whether a single-replicate scalar-vector dataset is passed as one big `feature×K`
     matrix or a `Vector` of `feature×1` arrays to `train`/`assess`.
   - Recommendation: store column-major RAW (works for both); Phase-4's Wave-0 smoke confirms the exact
     `train` signature against the pinned v0.2.1 (the CLAUDE.md "verify API against dev docs, not memory"
     note). No Phase-3 blocker.
   - **RESOLVED:** Store column-major RAW (one sample = one column) — this layout feeds either a
     `feature×K` matrix or a `Vector` of `feature×1` arrays without re-shaping the cache. The exact
     `train()`/`assess()` Z container is a **Phase-4 Wave-0** smoke-test detail against pinned
     NeuralEstimators v0.2.1, **explicitly out of Phase-3 scope** and not a Phase-3 blocker.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Julia | everything | ✓ | 1.12.6 (pinned) | — |
| SHA (stdlib) | D-05 content hash | ✓ | bundled | — |
| Base.Threads | D-12 parallelism | ✓ | bundled (default `nthreads`==1) | serial (D-12) |
| StatsBase | standardization/moments | ✓ | in spike env | — |
| NeuralEstimators | tensor-layout target (Phase 4) | ✓ | 0.2.1 pinned | — |
| JLD2 | sharded cache | ✗ (must add) | latest (registry) | none — must add (Wave 0) |
| Random123 | counter-based RNG | ✗ (must add) | latest (registry) | stdlib `Random` (worse repro; D-11 wants Random123) |

**Missing with no fallback:** `JLD2` (cache is the phase deliverable). **Missing with fallback:**
`Random123` (stdlib `Random` is a degraded fallback, but D-11 mandates Random123). Both are a Wave-0
add + re-freeze, with the resolve-risk gate confirming NeuralEstimators 0.2.1 is not downgraded.

## Sources

### Primary (HIGH confidence)
- RandomNumbers.jl / Random123 docs (`juliarandom.github.io/RandomNumbers.jl/stable/man/random123/`) —
  `Philox4x`/`Threefry4x` constructors, `set_counter!`, AbstractRNG + `rand` compatibility. VERIFIED.
- NeuralEstimators.jl dev docs (`API/core`, framework) — θ as **d×K matrix (parameters in columns)**;
  `train(estimator, sampler, simulator)`; `sampleposterior(estimator, Z, N)`; column-wise Z. HIGH.
- Local verification (`julia -e`): SHA stdlib works; `nthreads()` default 1; JLD2/Random123 absent from
  spike env; registry reachable; 1376×1028 8-divisibility. VERIFIED this session.
- Existing spike source (read this session): `forward.jl`, `prior.jl`, `calibration.jl`, `ghat.jl`,
  `contract.jl`, `src/colocalization.jl` (patch/correlation/_exclude_zero), `src/LoadImages.jl`
  (MultiChannelImage), `runtests.jl`. HIGH.

### Secondary (MEDIUM confidence)
- JLD2.jl docs (`juliaio.github.io/JLD2.jl/stable/`) — pure-Julia, HDF5-format, `@save`/`jldsave`;
  exact `jldopen` modes not fully quoted on the home page (CITED, standard idiom).
- CLAUDE.md §Technology Stack — JLD2/Random123/StatsBase mandate, BayesInteractomics cache pattern,
  "What NOT to Use". CITED.

### Tertiary (LOW confidence)
- BayesInteractomics sharded-cache idiom — referenced by CLAUDE.md/CONTEXT but repo not accessible this
  session; pattern reconstructed from standard atomic-write/JLD2 practice. ASSUMED (A4).
- Augmented-moment composition and imsize weights — recommendations, not locked (A1, A2).

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — JLD2/Random123/StatsBase/SHA are CLAUDE.md-mandated or stdlib; APIs verified.
- Architecture (sharded cache, hash guard, keyed RNG, leak-free loader): HIGH — composed from verified
  primitives + locked decisions.
- NeuralEstimators tensor layout: MEDIUM — d×K verified; exact Z container is a Phase-4 detail (A3).
- Pitfalls/validation: HIGH — gates map 1:1 to the 4 success criteria and reuse the established harness.

**Research date:** 2026-06-27
**Valid until:** ~2026-07-27 (stable Julia stack; re-verify NeuralEstimators v0.2.1 layout at Phase-4
Wave-0 against the pinned Manifest, per the CLAUDE.md "verify young API, not memory" note).

## RESEARCH COMPLETE
