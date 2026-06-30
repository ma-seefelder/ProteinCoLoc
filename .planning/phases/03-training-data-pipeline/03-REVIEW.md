---
phase: 03-training-data-pipeline
reviewed: 2026-06-30T00:00:00Z
depth: deep
files_reviewed: 7
files_reviewed_list:
  - spike/data/encode.jl
  - spike/data/seeding.jl
  - spike/data/generate.jl
  - spike/data/cache.jl
  - spike/data/hashguard.jl
  - spike/data/loader.jl
  - spike/test/test_data_pipeline.jl
findings:
  critical: 1
  warning: 4
  info: 1
  total: 6
status: issues_found
---

# Phase 03: Training-Data Pipeline — Code Review Report

**Reviewed:** 2026-06-30
**Depth:** deep (cross-file call chain + module-boundary analysis)
**Files Reviewed:** 7
**Status:** issues_found

## Summary

The Phase-3 pipeline is structurally sound. The atomic shard write mechanism, the
content-hash version guard, and the leak-free k-fold loader's design are all
correct. The decoupling constraint is fully satisfied: `git diff f581d95 -- src/
Project.toml Manifest.toml` shows zero changes to `src/` or the root manifests.
The `Threads.@threads`-over-shards concurrency model is race-free by construction
(each thread owns a distinct shard path, no shared mutable RNG, per-global-index
Philox4x keying). The `module Loader` namespace encapsulation and the absence of
any `standardize_all` symbol make val-into-train leakage structurally impossible.
The SHA-256 content hash (NUL-separated file bytes + canonical config string) is
correct and all seven hashed source files exist.

Five issues were found: one BLOCKER (bare catch swallows InterruptException in the
resume predicate), three WARNINGs (theta type inconsistency in load_fold, unguarded
includes that risk const-redefinition warnings, write_meta missing pre-commit
integrity check), one WARNING in the test suite (tautological reproducibility
assertion), and one CONVENTION (undocumented non-negative precondition on
master_seed).

---

## Critical Issues

### CR-01: Bare `catch` in `_loads_ok` swallows `InterruptException`

**File:** `spike/data/cache.jl:68-78`

**Issue:** The resume-by-skip predicate uses a bare `catch` with no binding
variable. In Julia, bare `catch` catches every exception type including
`InterruptException` (Ctrl-C), `OutOfMemoryError`, and `TaskFailedException`. When
a user presses Ctrl-C during the shard integrity check inside a long 50k-sample
`generate_cache` run, the interrupt is silently consumed and the function returns
`false`, causing the "done" shard to be regenerated instead of aborting the run.
The program continues rather than stops.

Current code:
```julia
function _loads_ok(path)
    isfile(path) || return false
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "theta") && haskey(f, "summary_min") &&
                haskey(f, "summary_aug") && haskey(f, "schema_version")
        end
    catch          # <-- swallows ALL exceptions
        return false
    end
end
```

**Fix:** Rebind the exception and rethrow interrupt-class exceptions:
```julia
function _loads_ok(path)
    isfile(path) || return false
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "theta") && haskey(f, "summary_min") &&
                haskey(f, "summary_aug") && haskey(f, "schema_version")
        end
    catch e
        e isa InterruptException && rethrow()
        return false
    end
end
```

---

## Warnings

### WR-01: `load_fold` returns `theta` as `Float64` — type inconsistency with docstring and Phase-4 expectations

**File:** `spike/data/loader.jl:193-194`

**Issue:** `Ztr` and `Zva` are explicitly converted to `Float32` via `Float32.(Ztr)`
and `Float32.(Zva)`, but `θtr = θ[:, train_idx]` and `θva = θ[:, val_idx]` remain
`Matrix{Float64}`. The docstring states "returns Float32 d×K tensors", which a
Phase-4 caller reading the contract would interpret as all returned matrices being
Float32. If NeuralEstimators' `train(estimator, θ, Z)` expects Float32 for both
inputs (the common default for GPU-targeting code, and NeuralEstimators v0.2.1
uses Float32 internally), the Float64 theta will cause an automatic scalar promotion
per sample, or a type error, depending on the batch dispatch path. The mismatch is
silent at this call site but surfaces as a runtime type error in Phase 4.

```julia
# current — inconsistent dtypes
return (Ztr = Float32.(Ztr), θtr = θ[:, train_idx],   # θtr is Float64
        Zva = Float32.(Zva), θva = θ[:, val_idx], zt = zt)  # θva is Float64
```

**Fix:** Either convert theta consistently:
```julia
return (Ztr = Float32.(Ztr), θtr = Float32.(θ[:, train_idx]),
        Zva = Float32.(Zva), θva = Float32.(θ[:, val_idx]), zt = zt)
```
Or document explicitly that theta stays Float64 (if NeuralEstimators accepts it):
update the docstring to "Ztr/Zva are Float32; θtr/θva are Float64."

---

### WR-02: `seeding.jl` and `encode.jl` included unconditionally in `generate.jl`

**File:** `spike/data/generate.jl:47-48`

**Issue:** The two includes are not guarded by `isdefined`, unlike the four
neighboring includes (contract.jl, prior.jl, forward.jl, cache.jl). Every call to
`include("generate.jl")` unconditionally re-includes both files, triggering Julia's
"WARNING: replacing module" / "WARNING: redefinition of constant" for
`HOLDOUT_SALT`, `FOLD_SALT`, `IMSIZE_SET`, `IMSIZE_WEIGHTS`, `IMSIZE_EXPECTED_COST`,
`N_AUG_MOMENTS`, `AUG_DIM`. In interactive development sessions (Revise.jl or
repeated `include` calls), this produces noisy warnings. With `--depwarn=error`
(tightened CI flag), const redefinition is an error.

```julia
# guarded (correct pattern):
isdefined(@__MODULE__, :build_mci)     || include(...)
isdefined(@__MODULE__, :sample_prior)  || include(...)
isdefined(@__MODULE__, :simulate_pair) || include(...)
# NOT guarded (inconsistency):
include(joinpath(@__DIR__, "seeding.jl"))   # line 47 — always re-includes
include(joinpath(@__DIR__, "encode.jl"))    # line 48 — always re-includes
isdefined(@__MODULE__, :write_shard) || include(...)
```

**Fix:** Apply the same isdefined guard pattern as the neighbours:
```julia
isdefined(@__MODULE__, :HOLDOUT_SALT)   || include(joinpath(@__DIR__, "seeding.jl"))
isdefined(@__MODULE__, :N_AUG_MOMENTS)  || include(joinpath(@__DIR__, "encode.jl"))
```

---

### WR-03: `write_meta` commits without a pre-integrity-check, unlike `write_shard`

**File:** `spike/data/cache.jl:117-130`

**Issue:** `write_shard` follows a defensive pattern: write `.tmp`, reopen to
verify at least one key is present, then atomically rename. `write_meta` writes
`.tmp` and immediately renames without reopening:

```julia
function write_meta(dir, config)
    path = joinpath(dir, "meta.jld2")
    tmp  = path * ".tmp"
    jldsave(tmp; hash=..., subhashes=..., ...)   # write
    mv(tmp, path; force = true)                   # commit — no integrity check
    return path
end
```

If `jldsave` flushes the file descriptor but the OS writes a partial page before
the `mv` (e.g., disk-full mid-write where `jldsave` didn't throw), the resulting
`meta.jld2` will be corrupt. On the next run, `JLD2.load(meta_path)` in
`open_or_invalidate` will throw a JLD2 parse error, making the entire cache
directory inaccessible without manual cleanup. The impact is low on well-behaved
storage but the inconsistency with `write_shard` is a latent robustness gap.

**Fix:** Mirror the `write_shard` pattern:
```julia
function write_meta(dir, config)
    path = joinpath(dir, "meta.jld2")
    tmp  = path * ".tmp"
    jldsave(tmp; hash = cache_hash(config), subhashes = subhashes(),
            config = config, N = config.N, shard_size = SHARD_SIZE,
            schema_version = SCHEMA_VERSION)
    JLD2.jldopen(tmp, "r") do f          # integrity check before commit
        @assert haskey(f, "hash") "meta integrity check failed: $tmp"
    end
    mv(tmp, path; force = true)
    return path
end
```

---

### WR-04: SC-4 `load_fold` reproducibility assertion is tautological

**File:** `spike/test/test_data_pipeline.jl:220-221`

**Issue:** The test asserting that `load_fold` is reproducible reads:

```julia
@test load_fold(d, 2; K = DP_K, master_seed = DP_LOAD_SEED).θva ==
      load_fold(d, 2; K = DP_K, master_seed = DP_LOAD_SEED).θva
```

Both sides are the same expression (identical arguments, same deterministic
function). For any purely deterministic function `f`, `f(x) == f(x)` is always
true — this is not evidence of the D-09 property "same master_seed ⇒ same
membership". The only scenario where this would fail is if `load_fold` internally
consumed a shared global RNG (which it does not, correctly), making the test a
vacuous check against the one bug pattern that the design already eliminates by
construction. An accidental bug in fold membership logic (e.g., an off-by-one in
the stride) would pass this test silently because both sides would produce the
same wrong output.

Note: the preceding line `@test all_folds(DP_N_LOAD; ...) == folds` does
provide a genuine reproducibility check for `all_folds`. The gap is specifically
the `load_fold` path.

**Fix:** Capture two independent calls:
```julia
fa = load_fold(d, 2; K = DP_K, master_seed = DP_LOAD_SEED)
fb = load_fold(d, 2; K = DP_K, master_seed = DP_LOAD_SEED)
@test fa.θva == fb.θva   # same args → identical fold membership
```
For a stricter gate, also assert that a different `master_seed` produces a
different split:
```julia
fc = load_fold(d, 2; K = DP_K, master_seed = DP_LOAD_SEED + 1)
@test fa.θva != fc.θva   # different seed → different fold (probabilistically)
```

---

## Info

### IN-01: `sample_rng`/`holdout_rng`/`fold_rng` silently throw for negative `master_seed`

**File:** `spike/data/seeding.jl:75, 87, 96`

**Issue:** All three RNG constructors call `UInt64(master_seed)`, which throws
`InexactError: convert(UInt64, -1)` for any negative `Integer` argument. The type
signature `master_seed::Integer` implies acceptance of any integer, but the
functions only operate correctly for non-negative values. All current call sites
use non-negative literals, so this is not a current bug but an undocumented
precondition that will surface unexpectedly if the API is extended.

**Fix:** Document the precondition in the docstrings or add a runtime guard:
```julia
sample_rng(master_seed::Integer, idx::Integer) =
    (master_seed >= 0 || error("master_seed must be non-negative, got $master_seed");
     Philox4x(UInt64, (UInt64(master_seed), UInt64(idx))))
```
Alternatively, restrict the signature to `UInt64` directly, matching the Philox4x
key type and making the constraint self-documenting.

---

## Decoupling Constraint Verification (PASS)

`git diff f581d95 -- src/ Project.toml Manifest.toml` returns empty: zero edits to
`src/` or the root manifests since the Phase-3 baseline. The `hashguard.jl`
inclusion of `src/colocalization.jl` and `src/LoadImages.jl` in `HASH_SRC_FILES`
is read-only (SHA-256 of file bytes); these files are not modified by any Phase-3
code. The spike's own `Manifest.toml` (new file, `spike/Manifest.toml`) is
appropriately scoped. Decoupling constraint: SATISFIED.

---

## Concurrency Model Verification (PASS)

`generate_cache` launches `Threads.@threads for s in 1:nshards`. Each thread
works on a unique shard index `s`, producing a unique temp path `shard_NNNN.jld2.tmp`
and a unique final path `shard_NNNN.jld2`. No shared mutable state between threads.
The generation core (`generate_samples`) creates a fresh `Philox4x` RNG per column
keyed by `(master_seed, global_idx)`, making the output a pure function of the
global index with no RNG state shared across columns or threads. The `write_meta`
call is sequential after the `@threads` block (Julia's `@threads` is synchronous).
Race-freedom claim: VERIFIED.

---

## Leak-Free Loader Verification (PASS)

The `load_fold` function:
1. fits `ZScoreTransform` on `Z[cont_rows, train_idx]` only — never sees val columns
2. applies the transform to both train and val columns (val is standardized with
   train statistics, which is the correct and intended behavior)
3. passes binary mask rows (65:128) through unchanged — z-scoring a 0/1 mask would
   re-couple folds through the mask mean; the bypass is correct
4. `load_main_pool` globs only `shard_*.jld2`, never `holdout.jld2` (the D-10
   structural exclusion)

No path exists through the public API (`load_main_pool`, `load_fold`, `load_holdout`,
`all_folds`) that would allow val statistics to contaminate train standardization.
The `module Loader` namespace prevents accidental export of any internal symbol.
`names(Loader)` surface check in SC-3a is correct. Leak-freedom: VERIFIED by construction.

---

_Reviewed: 2026-06-30_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: deep_
