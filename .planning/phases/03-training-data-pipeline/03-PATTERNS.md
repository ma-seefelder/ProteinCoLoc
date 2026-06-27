# Phase 3: Training-Data Pipeline - Pattern Map

**Mapped:** 2026-06-27
**Files analyzed:** 11 (6 new modules, 1 new test, 4 modified)
**Analogs found:** 11 / 11 (all in-repo; BayesInteractomics external pattern noted from RESEARCH.md)

> **Hard constraint (every file below):** ALL new code lives under `spike/`. NEVER edit `src/`
> or the root `Project.toml`/`Manifest.toml` (decoupling baseline `f581d95`, CLAUDE.md). New
> generator modules reach `src/` ONLY transitively through `spike/contract.jl`'s read-only
> `include()`. The standing DEMO-02 gate (`git diff` vs `f581d95 -- src`) must stay clean.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `spike/data/seeding.jl` (new) | utility | transform | `spike/simulator/prior.jl` (const block + `rng::AbstractRNG`) | role-match |
| `spike/data/encode.jl` (new) | utility | transform | `spike/contract.jl` (`patch_summary`/`induced_mu` pure transforms of the 8×8 matrix) | exact |
| `spike/data/hashguard.jl` (new) | utility / config | transform | `spike/simulator/calibration.jl` `emit_ghat` (read source bytes, deterministic serialization) | role-match |
| `spike/data/cache.jl` (new) | persistence | file-I/O | `spike/simulator/calibration.jl` `emit_ghat` (`write(path, content)` atomic-ish artifact write) | partial (new JLD2 API) |
| `spike/data/generate.jl` (new) | generator / script | batch + event-driven (resume) | `spike/simulator/forward.jl` (validation-at-entry) + `prior.jl` (`sample_prior` drive loop) | role-match |
| `spike/data/loader.jl` (new) | loader | transform / CRUD-read | `spike/contract.jl` (module of pure functions over loaded data) | role-match |
| `spike/test/test_data_pipeline.jl` (new) | test | request-response | `spike/test/test_simulator.jl` | exact |
| `spike/test/runtests.jl` (modify) | test harness | request-response | itself (existing `include` wiring at L71) | exact |
| `spike/Project.toml` (modify) | config | — | itself (existing `[deps]` block) | exact |
| `spike/Manifest.toml` (modify) | config | — | Phase-2 Wave-0 resolve+refreeze | exact |
| `.gitignore` (modify) | config | — | root `.gitignore` | exact |

**Shared header (ALL new `.jl` files):** every spike `.jl` file opens with the 19-line AGPL
`#= … =#` license block (verbatim from `spike/contract.jl:1-19`), then a `# spike/<path> --- <one-line role>`
comment, then a short design-rationale paragraph citing the governing D-xx decision. Copy this exactly.

## Pattern Assignments

### `spike/data/seeding.jl` (utility, transform — D-11/D-10)

**Analog:** `spike/simulator/prior.jl` — for the const-config-block + `rng::AbstractRNG` docstring style.

**Header + imports pattern** (`prior.jl:34-37`): license block, then narrow `using`:
```julia
using Random                 # AbstractRNG
include(joinpath(@__DIR__, "ghat.jl"))   # sibling-include via @__DIR__
```
New file uses `using Random123` instead. Per-sample keyed RNG (RESEARCH Pattern 4):
```julia
# main pool:   Philox4x(UInt64, (UInt64(master_seed),                 UInt64(idx)))
# holdout:     Philox4x(UInt64, (UInt64(master_seed) ⊻ HOLDOUT_SALT,  UInt64(idx)))
const HOLDOUT_SALT = 0x... # fixed nonzero UInt64 (RESEARCH Pattern 5)
const FOLD_SALT    = 0x... # disjoint key stream for loader fold permutation
```
**Why this analog:** `prior.jl`'s `sample_prior(rng::AbstractRNG)` (L71) is the exact interface
the keyed RNG feeds — `Philox4x` is an `AbstractRNG`, so it drops into `sample_prior(rng_i)` and
`simulate_pair(rng_i, θ; imsize)` with ZERO simulator change (RESEARCH §D-11, VERIFIED).

**Const-naming convention to copy** (`prior.jl:40-57`): `SCREAMING_SNAKE_CASE` module consts with
an inline `#` rationale comment per const. Apply to `HOLDOUT_SALT`, `FOLD_SALT`, `IMSIZE_SET`, weights.

---

### `spike/data/encode.jl` (utility, pure transform — D-01/D-02)

**Analog:** `spike/contract.jl` — `patch_summary` (L85-90) and `induced_mu` (L100-107) are the exact
pure-transform-of-the-8×8-matrix pattern `encode_d01`/`encode_aug` extend.

**Core transform pattern** (`contract.jl:100-107`, the missing-handling idiom to mirror):
```julia
function induced_mu(mci::MultiChannelImage)
    vals = collect(skipmissing(patch_summary(mci)))
    isempty(vals) && return NaN          # CR-02: fully-missing degeneracy handled, NOT dropped
    return Statistics.mean(vals)
end
```
**D-01 encode (RESEARCH §The 128-Dim D-01 Encoding):**
```julia
function encode_d01(M::AbstractMatrix)   # M is 8×8 Union{Float64,Missing}
    vals = vec(coalesce.(M, 0.0))                 # 64-dim, missing→0, column-major
    mask = vec(Float64.(.!ismissing.(M)))         # 64-dim, 1=present / 0=missing
    return vcat(vals, mask)                        # 128-dim; rows 1:64 vals, 65:128 mask
end
```
**Critical conventions inherited from the contract:**
- The 8×8 grid stays FIXED (D-10 / CLAUDE.md); `encode_*` wraps `patch_summary(mci)`'s matrix, never
  re-computes or edits the `src/` summary.
- `missing` is a first-class value here: `patch_summary` returns `Matrix{Union{Float64,Missing}}`
  (a patch is `missing` when `_exclude_zero` leaves ≤15 px, `src/colocalization.jl:235`). The mask
  channel makes that self-describing. A fully-missing 8×8 → vals all 0, mask all 0 — cache it, do
  NOT drop (mirrors `induced_mu`'s `NaN`-not-throw choice; preserves the π(θ)-faithful distribution, D-13).
- Augmented moments (D-02) use `mci.otsu_threshold` (already on the struct, `contract.jl:68`) for
  Manders, and `StatsBase`/`Statistics` (`median`/`quantile`/`std`/`skewness`/`kurtosis`) over
  `skipmissing(M)` — spike-local helpers, NO `src/` edit. Store as a SUPERSET column (128 + ~12).

**Docstring style to copy:** `contract.jl` triple-quoted docstrings with a signature line, prose, and
the explicit src line references (e.g. `src/colocalization.jl:221`). Document the mask layout here so
Phase-4's summary-net input dim is unambiguously 128.

---

### `spike/data/hashguard.jl` (utility/config, transform — D-05)

**Analog:** `spike/simulator/calibration.jl` `emit_ghat` (L250-322) — the only existing code that
reads/serializes source-defining content deterministically and writes a derived artifact.

**Content-hash pattern (RESEARCH Pattern 3):**
```julia
using SHA   # stdlib, no add
function cache_hash(config::NamedTuple; src_files)
    h = SHA.SHA256_CTX()
    for f in src_files                              # FIXED order = deterministic
        SHA.update!(h, read(f)); SHA.update!(h, UInt8['\0'])
    end
    SHA.update!(h, Vector{UInt8}(canonical(config)))   # sorted-key deterministic string
    return bytes2hex(SHA.digest!(h))
end
```
**Recommended `src_files`** (RESEARCH Pattern 3, A5): `spike/simulator/forward.jl`,
`spike/simulator/prior.jl`, `spike/simulator/ghat.jl`, `spike/contract.jl`, plus `spike/data/encode.jl`
(once it exists). Optionally the two frozen `src/` files (`colocalization.jl`, `LoadImages.jl`) as a
defensive decoupling-breach check. Use `@__DIR__`-relative joins exactly as `prior.jl:37` /
`contract.jl:46-47` do. Store full hash + per-file sub-hashes + raw config in `meta.jld2`.

**Config fields hashed** (D-05): N, shard_size, imsize-spec (set + weights), summary-spec (D-01 mask
convention + augmented-moment list), θ-dim (7), master_seed, k.

---

### `spike/data/cache.jl` (persistence, file-I/O — D-04/D-06)

**Analog:** `calibration.jl` `emit_ghat` L320 (`write(path, content)`) for the "build artifact, write
to path, return path" shape. The JLD2 + atomic-rename API itself is new to the repo (no existing JLD2
call site — this is the Wave-0 dependency add).

**Atomic shard write + resume-by-skip (RESEARCH Pattern 2):**
```julia
tmp = path * ".tmp"
jldsave(tmp; theta, summary_min, summary_aug, global_index, imsize, schema_version=1)
JLD2.jldopen(tmp, "r") do f; @assert haskey(f, "theta"); end   # integrity check
mv(tmp, path; force=true)                                       # atomic commit
# resume: isfile(path) && _loads_ok(path)  → skip this shard
```
**Column-major layout (RESEARCH Pattern 1 — matches NeuralEstimators d×K):** one sample = one column.
```julia
theta        = Matrix{Float64}(undef, 7, n)      # θ columns (D-04 7-vector)
summary_min  = Matrix{Float64}(undef, 128, n)    # D-01: rows 1:64 vals, 65:128 mask
summary_aug  = Matrix{Float64}(undef, A, n)      # D-02 augmented superset (A ≈ 128 + 12)
```
**Cache dir:** `spike/data/cache/<sha256>/` with `shard_NNNN.jld2`, `holdout.jld2`, `meta.jld2`
(RESEARCH §Recommended Project Structure). One-shot `jldsave` per shard (a 10k-sample shard ≈ 10–20 MB,
fits trivially) — NOT in-file append (RESEARCH Anti-Pattern). BayesInteractomics is the cited reference
for this sharded pattern (CLAUDE.md; repo not accessible this session, so pattern is from RESEARCH §A4).

---

### `spike/data/generate.jl` (generator/script, batch + resume — D-12/D-03)

**Analog (validation-at-entry):** `spike/simulator/forward.jl:110-134` — surface, don't swallow, the
simulator's `ArgumentError`s (ASVS V5, RESEARCH §Security). **Analog (drive loop):** `prior.jl`
`sample_prior` is the per-sample entry the generator calls.

**Per-sample chain (RESEARCH §Code Examples, drives UNCHANGED Phase-2 code):**
```julia
rng = Philox4x(UInt64, (master_seed, UInt64(idx)))   # seeding.jl
θ   = sample_prior(rng)                                # UNCHANGED prior.jl
isz = sample_imsize(rng)                               # D-03 nuisance (seeding.jl / encode config)
mci = build_mci(simulate_pair(rng, θ; imsize=isz))     # UNCHANGED forward.jl + contract.jl
M   = patch_summary(mci)                               # UNCHANGED 8×8 summary
s_min = encode_d01(M); s_aug = encode_aug(mci, M)      # encode.jl
```
**Threaded generation (D-12, RESEARCH §Threads.@threads):** `Threads.@threads` over SHARDS (coarse;
each thread builds one shard into thread-local arrays, then atomically writes — no shared mutable
state, no locks). Key-per-global-index ⇒ byte-identical dataset for any `nthreads`; serial when
`nthreads()==1` (the verified machine default). Document the `-t auto` / `JULIA_NUM_THREADS` launch
requirement. Reserved holdout (D-10) drawn from the XOR-salted disjoint key namespace into a separate
`holdout.jld2`, NEVER indexed into the main pool.

**Sibling-include convention** (copy from `test_simulator.jl:43,92,176`): `include(joinpath(@__DIR__, "...jl"))`
to pull in `contract.jl` (→ `src/` read-only), `simulator/prior.jl`, `simulator/forward.jl`, and the
new `data/*.jl` modules. Order matters: `contract.jl` first (brings `build_mci`/`patch_summary` into scope).

---

### `spike/data/loader.jl` (loader, transform — D-07/D-08)

**Analog:** `spike/contract.jl` as a module of pure functions over loaded data; the leak-free
fit-on-train standardization is a new structural composition (RESEARCH Pattern 6).

**The footgun-removal pattern (RESEARCH Pattern 6 — the ONLY standardization entry point):**
```julia
using StatsBase
function load_fold(cache, fold::Int; K::Int=5, master_seed)
    Z, θ = load_main_pool(cache)                           # RAW feature×N, columns = samples
    perm = randperm(fold_rng(master_seed), size(Z,2))      # deterministic from master seed
    val_idx   = perm[fold:K:end]
    train_idx = setdiff(1:size(Z,2), val_idx)
    zt  = fit(ZScoreTransform, Z[:, train_idx]; dims=2)    # ◄ FIT ON TRAIN ONLY
    Ztr = StatsBase.transform(zt, Z[:, train_idx])
    Zva = StatsBase.transform(zt, Z[:, val_idx])           # applied to val, never fit on it
    return (Ztr=Float32.(Ztr), θtr=θ[:,train_idx],
            Zva=Float32.(Zva), θva=θ[:,val_idx], zt)        # return zt → Phase-5 freezes train-only prep
end
```
**Structural constraints (locked):**
- NO public `standardize_all` / global-stats function exists — cross-fold leakage is impossible by
  construction, not by convention (D-08; RESEARCH Anti-Patterns). The SC-3 test asserts the symbol
  is NOT exported.
- The cache is RAW (D-07); the loader is the SOLE standardization path.
- **Mask rows (65:128) are NOT z-scored** — they are binary indicators; standardize only the value
  rows (1:64 of D-01 + continuous moment rows of D-02), pass the mask through unchanged (RESEARCH
  Pattern 6 / Anti-Pattern; z-scoring a 0/1 mask re-couples folds through the mask mean).
- `fold_rng` keyed by `(master_seed ⊻ FOLD_SALT, 0)` (seeding.jl) — folds reproducible and independent
  of the generation key stream. Holdout (`holdout.jld2`) is never in `load_main_pool`'s `1:N` index range.

---

### `spike/test/test_data_pipeline.jl` (test — SC-1..SC-4 + thread-repro)

**Analog:** `spike/test/test_simulator.jl` — copy its EXACT structure.

**Structure to mirror** (`test_simulator.jl:37-43`): license header, `using Test` / `using Random`,
then `include(joinpath(@__DIR__, "..", "contract.jl"))` and the new `data/*.jl` modules, then nested
`@testset "..." verbose = true begin … end` blocks. Determinism via fixed seeds (`test_simulator.jl:121-124`):
```julia
@testset "determinism (...)" begin
    @test simulate_pair(Random.Xoshiro(7), _θ(0.3)) == simulate_pair(Random.Xoshiro(7), _θ(0.3))
end
```
**Pre-declared tolerance consts** (`test_simulator.jl:240-241` style): declare gate thresholds as
`const` BEFORE the testset, NOT tuned to pass.

**Gate map (RESEARCH §Validation Architecture, copy the test idioms):**
- SC-1: `size(summary_min)==(128,64)`; rows 65:128 ∈ {0,1}; mask=0 ⟺ value row==0; all `isfinite`; θ is the 7-vector.
- SC-2: (a) write→reload round-trip array equality; (b) delete 1 shard, re-run, only it regenerates,
  dataset identical; (c) perturb a tracked source byte/config field → `cache_hash` differs → loader refuses stale.
- SC-3: `@test !(:standardize_all in names(LoaderModule))`; `Ztr` rows ≈0 mean/unit std; `Zva` NOT exactly 0/1; mask rows unchanged.
- SC-4: `union(folds)==1:N` & pairwise `∩==∅`; `holdout ∩ every fold == ∅`; `length(holdout) ≥ 20`; reproducible across two loader calls.
- D-11/D-12: generate with `nthreads()==1` and `>1` → byte-identical θ + summaries.

Use a tiny N≈64 fixture (sub-30s quick gate). The fully-finite, strictly-positive synthetic image
idiom (`test_simulator.jl:50-54`) is the reference for any hand-built fixture image.

---

### `spike/test/runtests.jl` (modify — wire in the new testset)

**Analog:** itself, line 71 — the existing one-line include that wires `test_simulator.jl` into the harness:
```julia
include(joinpath(@__DIR__, "test_simulator.jl"))
```
Add the parallel line for `test_data_pipeline.jl` so a single `julia --project=spike spike/test/runtests.jl`
remains the gate. The resolve-risk gate at `runtests.jl:64` (asserts NeuralEstimators stays pinned at
v0.2.1 by UUID) MUST be re-verified after adding JLD2/Random123 — extend it to confirm the new deps do
NOT downgrade NeuralEstimators (the Phase-1 GLMakie co-resolve regression class).

---

### `spike/Project.toml` (modify — Wave-0 dependency add)

**Analog:** itself — the existing `[deps]` block. Add two UUID lines (alphabetical), matching the exact format:
```toml
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
Random123 = "74087812-796a-5b5d-8853-05524746bad3"
```
(StatsBase, NeuralEstimators already present.) Then `Pkg.resolve()`, re-freeze `spike/Manifest.toml`,
re-run the resolve-risk gate. NEVER touch the root `Project.toml`/`Manifest.toml` (CLAUDE.md hard constraint).

---

### `.gitignore` (modify — gitignore the bulk cache)

**Analog:** root `.gitignore`. Add a rule to ignore `spike/data/cache/` bulk artifacts while keeping a
tiny committed smoke fixture for the round-trip/invalidation tests (RESEARCH §Runtime State Inventory).

## Shared Patterns

### License + module header (ALL new `.jl` files)
**Source:** `spike/contract.jl:1-19` (the 19-line AGPL `#= … =#` block), then `spike/contract.jl:21-47`
(role one-liner + design-rationale paragraph + narrow `using` + `@__DIR__`-relative `include`s).
**Apply to:** `seeding.jl`, `encode.jl`, `hashguard.jl`, `cache.jl`, `generate.jl`, `loader.jl`,
`test_data_pipeline.jl`. Verbatim license block; module-specific one-liner cites the governing D-xx.

### Read-only `src/` coupling
**Source:** `spike/contract.jl:45-47` — the ONLY place new code may reach `src/`, via `include()`.
```julia
include(joinpath(@__DIR__, "..", "src", "LoadImages.jl"))      # MultiChannelImage
include(joinpath(@__DIR__, "..", "src", "colocalization.jl"))  # patch, correlation, _exclude_zero
```
**Apply to:** All generator/encode/loader files reach `build_mci`/`patch_summary`/`correlation`
transitively by including `contract.jl` — they NEVER include `src/` directly and NEVER edit it.

### Explicit `rng::AbstractRNG` threading
**Source:** `spike/simulator/prior.jl:71` (`sample_prior(rng::AbstractRNG)`), `forward.jl:110`
(`simulate_pair(rng::AbstractRNG, θ; imsize)`).
**Apply to:** `seeding.jl` (produces `Philox4x` RNGs), `generate.jl` (threads one keyed rng per sample),
`loader.jl` (`fold_rng`). Random123 generators ARE `AbstractRNG`s → drop in with zero callee changes.

### Const-config block with per-const rationale
**Source:** `spike/simulator/prior.jl:40-57` and `forward.jl:55-66` — `SCREAMING_SNAKE_CASE` consts,
each with an inline `#` rationale citing the D-xx it satisfies.
**Apply to:** `seeding.jl` (salts, imsize set/weights), `cache.jl` (shard size, schema version),
`hashguard.jl` (src-file list).

### Validation-at-entry, surface-don't-swallow
**Source:** `spike/simulator/forward.jl:110-134` — `throw(ArgumentError(...))` on every out-of-range
input at function entry.
**Apply to:** `generate.jl` must let `simulate_pair`'s `ArgumentError`s propagate (ASVS V5); `cache.jl`
integrity check `@assert haskey(f, "theta")` before atomic commit; `loader.jl` `error()` on hash mismatch.

### Test harness (stdlib Test, deterministic-seed gates)
**Source:** `spike/test/test_simulator.jl` (whole file) + `runtests.jl:71` (include wiring).
**Apply to:** `test_data_pipeline.jl` (mirror structure) and the one-line `runtests.jl` addition.
Pre-declared `const` thresholds; `verbose = true` testsets; fixed seeds for reproducible (non-flaky) gates.

## No Analog Found

No file lacks an analog, but two patterns are NEW to the repo (no existing call site) — the planner
should source them from RESEARCH.md, not a repo file:

| Pattern | Role | Data Flow | Reason | Source |
|---------|------|-----------|--------|--------|
| JLD2 `jldsave`/`jldopen` cache I/O | persistence | file-I/O | No existing JLD2 call site (Wave-0 dep add); `calibration.jl`'s `write(path, content)` is the closest "write artifact" shape but not JLD2 | RESEARCH Patterns 1-2; JLD2 docs |
| Random123 `Philox4x` keyed RNG | utility | transform | No existing Random123 call site (Wave-0 dep add); existing code uses `Xoshiro`/`MersenneTwister` | RESEARCH Pattern 4; RandomNumbers.jl docs (VERIFIED) |

Both are CLAUDE.md-mandated, registry-legitimate, and VERIFIED in RESEARCH §Package Legitimacy Audit.
The interface they plug into (`rng::AbstractRNG`, `write`-then-`mv` atomicity) IS analog-backed above.

## Metadata

**Analog search scope:** `spike/` (all `.jl`), `spike/test/`, `spike/simulator/`, `spike/Project.toml`,
root `.gitignore`. `src/` inspected read-only via existing `contract.jl` references (not edited).
**Files scanned:** 9 spike `.jl` files + Project.toml + .gitignore (full repo `.jl` inventory under spike).
**External reference noted (not accessible this session):** BayesInteractomics JLD2 sharded-cache idiom
(CLAUDE.md-cited; reconstructed from RESEARCH §A4 standard atomic-write practice).
**Pattern extraction date:** 2026-06-27
