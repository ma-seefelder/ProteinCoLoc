---
phase: 03-training-data-pipeline
verified: 2026-06-30T00:00:00Z
status: passed
score: 14/14 must-haves verified
has_blocking_gaps: false
re_verification: null
gaps: []
human_verification: []
---

# Phase 3: Training-Data Pipeline Verification Report

**Phase Goal:** A reproducible, resumable generator turns the prior π(θ) + Phase-2 simulator + reused
(UNCHANGED) summary functions into cached standardized training vectors, with leak-free split discipline
baked into the loader structurally.
**Verified:** 2026-06-30
**Status:** PASSED
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

All 14 truths derived from the five plan `must_haves` blocks (DATA-01 via Plans 01–02, DATA-02 via Plans
01, 03, 05, DATA-03 via Plan 04) plus the hard decoupling gate.

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `using JLD2, Random123` succeeds in the spike env; NeuralEstimators stays pinned at v0.2.1 | VERIFIED | `spike/Project.toml` [deps] contains both UUIDs; runtests.jl resolve-risk gate (lines 71–74) re-asserts v0.2.1 with those deps present; full suite green |
| 2 | `test_data_pipeline.jl` is wired into the single `runtests.jl` gate | VERIFIED | `runtests.jl:86` — `include(joinpath(@__DIR__, "test_data_pipeline.jl"))` |
| 3 | `encode_d01(M)` turns an 8×8 patch-correlation matrix into a 128-dim vector (64 imputed values + 64 binary mask) | VERIFIED | `encode.jl:71–74`: `vcat(coalesce.(M,0.0)\|vec, Float64.(.!ismissing.(M))\|vec)` = 128-dim; spot-check confirmed (length=128, rows 1:64==vec(M), rows 65:128 all 1.0, missing position → val=0/mask=0, fully-missing → 128-zeros, not dropped) |
| 4 | `encode_aug` produces AUG_DIM=142 superset = 128 D-01 dims + 14 moments in one pass | VERIFIED | `encode.jl:96–137`: `vcat(encode_d01(M), moments)`, `const N_AUG_MOMENTS=14`, `const AUG_DIM=142`; SC-1 asserts `size(summary_aug)==(142,64)` (passes) |
| 5 | Each sample's RNG is a pure function of `(master_seed, global_index)` — order- and thread-independent | VERIFIED | `seeding.jl:75–76`: `Philox4x(UInt64, (UInt64(master_seed), UInt64(idx)))`; D-11/D-12 testset (9 assertions): parallel==serial byte-identical AND shuffled-index re-sorted == in-order — all pass |
| 6 | `generate_samples(N)` maps θ~π → simulate_pair → patch_summary → 128-dim column-major vectors | VERIFIED | `generate.jl:72–85` (`generate_sample`) chains the UNCHANGED prior/simulator/contract; `generate.jl:103–137` (`generate_samples`) pre-allocates column-major matrices, fills per-column via keyed sample, dispatches serial/parallel |
| 7 | 50k–200k samples write to sharded JLD2 chunks and reload byte-identically; resume skips completed shards | VERIFIED | `cache.jl:96–107` atomic `.tmp`→integrity-check→`mv`; `shard_done` predicate; SC-2a/b tests pass; at-scale: 52-file cache dir on disk (50 shards + holdout.jld2 + meta.jld2, NOTES.md §6 records 18.7 min / reload exactly 50000 cols / resume regenerates only deleted shard) |
| 8 | A changed source byte OR config field auto-invalidates the stale cache | VERIFIED | `hashguard.jl:76–84`: SHA-256 over FIXED-order source bytes + `canonical(config)`; `open_or_invalidate` errors on mismatch; SC-2c `@test_throws ErrorException` passes; hashguard testset (5 assertions) passes |
| 9 | The reserved ≥20-stack ADVI holdout is in a separate file from a disjoint key namespace, never in the main pool | VERIFIED | `generate.jl:171–198` (`_write_holdout`) uses `holdout_rng` (XOR-salted, `seeding.jl:86–87`); stores `global_index` as negative; `cache.jl` never reads `holdout.jld2`; SC-2 asserts `intersect(pool.global_index, holdout.global_index)==∅` (passes) |
| 10 | A single loader API returns per-fold standardized tensors with z-scoring fit on TRAIN folds only | VERIFIED | `loader.jl:166–195` (`load_fold`): `zt = fit(ZScoreTransform, Z[cont_rows, train_idx]; dims=2)` — val columns never enter the fit; SC-3b asserts `Ztr` rows 1:64 ≈0/unit-std AND `Zva` rows not exactly 0/1 (passes) |
| 11 | No public global-standardize function exists — cross-fold leakage impossible by construction | VERIFIED | `loader.jl:61`: `export load_main_pool, load_fold, load_holdout, all_folds` (4 functions only); `module Loader` boundary; spot-check confirms `names(Loader)` = {Loader, all_folds, load_fold, load_holdout, load_main_pool}; SC-3a `!(:standardize_all in names(Loader))` passes |
| 12 | The binary mask rows (65:128) pass through standardization unchanged | VERIFIED | `loader.jl:126–136` (`_row_partition`): cont=1:64, mask=65:128 for `:min`; mask rows copied raw: `Ztr[mask_rows,:] = Z[mask_rows, train_idx]` (line 190); SC-3b mask-bypass assertions pass |
| 13 | k=5 folds disjoint, complete, reproducible from master seed; reserved ≥20 holdout excluded from every fold by construction | VERIFIED | `loader.jl:146–149` (`all_folds`): single `randperm(fold_rng(master_seed), N)` strided K ways; `_shard_files` globs `shard_*` only (never `holdout.jld2`); SC-4a (disjoint+complete+reproducible) and SC-4b (structural holdout exclusion: count==N, hide holdout leaves pool unchanged, HOLDOUT_SALT≠0) — all 19 assertions pass |
| 14 | HARD GATE: `git diff --quiet f581d95 -- src Project.toml Manifest.toml` is clean (spike decoupling) | VERIFIED | `git diff --quiet f581d95 -- src Project.toml Manifest.toml` → EXIT 0 / output: CLEAN (verified during this session) |

**Score: 14/14 truths verified**

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `spike/data/encode.jl` | 128-dim D-01 encoding + AUG_DIM=142 D-02 superset | VERIFIED | 138 lines; `encode_d01` (line 71) + `encode_aug` (line 96); `const AUG_DIM=142`; AGPL header present |
| `spike/data/seeding.jl` | Philox4x keyed RNG, HOLDOUT_SALT, FOLD_SALT, sample_imsize | VERIFIED | 115 lines; `sample_rng`/`holdout_rng`/`fold_rng` all present; `HOLDOUT_SALT=0x9E3779B97F4A7C15` (non-zero); `IMSIZE_SET` + `IMSIZE_WEIGHTS` cost-aware categorical |
| `spike/data/generate.jl` | generate_sample/generate_samples + generate_cache driver | VERIFIED | 243 lines; `generate_sample` (line 72), `generate_samples` (line 103), `generating_config` (line 150), `_write_holdout` (line 171), `generate_cache` (line 216) |
| `spike/data/hashguard.jl` | SHA-256 content hash over source files + canonical config | VERIFIED | 96 lines; `using SHA`; `cache_hash` (line 76) folds `SHA256_CTX` over FIXED-order `HASH_SRC_FILES` + `canonical(config)`; `subhashes` (line 93) for diagnosability |
| `spike/data/cache.jl` | Atomic shard write, resume-by-skip, meta.jld2, open_or_invalidate | VERIFIED | 161 lines; `write_shard` (line 96) uses `.tmp`→integrity-check→`mv`; `shard_done` (line 86); `write_meta` (line 117); `open_or_invalidate` (line 143) errors on hash mismatch |
| `spike/data/loader.jl` | module Loader: sole standardization path (D-07/D-08/D-09/D-10) | VERIFIED | 204 lines; `module Loader` (line 51); exports 4 functions only; `load_fold` fits `ZScoreTransform` on train only; `_row_partition` bypasses mask; `all_folds` uses `fold_rng` |
| `spike/test/test_data_pipeline.jl` | SC-1..SC-4 + D-11/D-12 + hashguard + fixture testsets | VERIFIED | 294 lines; 7 testsets, 0 `@test_skip` placeholders (all filled); 68 assertions total |
| `spike/test/runtests.jl` | wired include of test_data_pipeline.jl + extended resolve-risk gate | VERIFIED | Line 86: `include(joinpath(@__DIR__, "test_data_pipeline.jl"))`; lines 71–74: JLD2/Random123 presence + NE v0.2.1 re-assert |
| `spike/Project.toml` | JLD2 + Random123 in [deps] | VERIFIED | Line 11: `JLD2 = "033835bb..."`, Line 12: `Random123 = "74087812..."` |
| `spike/data/cache/fixture/shard_0001.jld2` | Committed schema-stability fixture | VERIFIED | File exists; fixture testset (6 assertions) confirms theta(7×?), summary_min(128×?), summary_aug(AUG_DIM×?), schema_version=1, mask rows ∈ {0,1} |
| `spike/NOTES.md` §6 | At-scale 50k run evidence: N, wall-clock, shard count, thread-repro | VERIFIED | §6 records N=50000, 18.7 min / 32 threads, 50 shards, reload exactly 50000 cols, resume (only deleted shard regenerated), cross-process -t1 vs -t auto byte-identical; at-scale cache dir exists on disk (52 files) |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `spike/test/runtests.jl` | `test_data_pipeline.jl` | `include(joinpath(@__DIR__, "test_data_pipeline.jl"))` | WIRED | Line 86 |
| `spike/test/runtests.jl` | NeuralEstimators v0.2.1 | UUID-keyed `Pkg.dependencies()` assertion | WIRED | Lines 64 + 74 (both assert v0.2.1 after adding JLD2/Random123) |
| `spike/data/generate.jl` | `seeding.jl` + `encode.jl` | `include` of siblings (lines 47–48) | WIRED | Guarded `isdefined` includes keep file idempotent under standalone + runtests loading |
| `spike/data/encode.jl` | `patch_summary` 8×8 matrix | `coalesce` + `ismissing` over the matrix | WIRED | Lines 72–73 |
| `spike/data/cache.jl` | `cache/<hash>/meta.jld2` | `cache_hash` names the dir; `write_meta` stores full hash | WIRED | `open_or_invalidate` (line 143) + `write_meta` (line 117) |
| `spike/data/generate.jl` | `holdout.jld2` | `holdout_rng` disjoint namespace via `_write_holdout` | WIRED | Lines 171–198; XOR-salted RNG key, separate file, negative `global_index` |
| `spike/data/loader.jl` | `fit(ZScoreTransform, Z[cont_rows, train_idx]; dims=2)` | fit-on-train-only standardization | WIRED | Line 182 |
| `spike/data/loader.jl` | `fold_rng(master_seed)` | `randperm(fold_rng(...), N)` deterministic split | WIRED | Lines 147 + 175 |

---

### Data-Flow Trace (Level 4)

All dynamic-data artifacts (generate.jl, cache.jl, loader.jl) were verified at Level 4 via behavioral spot-checks and the test suite. No hollow wiring found.

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|--------------------|--------|
| `generate.jl generate_sample` | `s_min` (128-dim) | `patch_summary(build_mci(simulate_pair(rng,θ)))` chain | Yes — UNCHANGED Phase-2 simulator + contract functions | FLOWING |
| `cache.jl write_shard` | `theta/summary_min/summary_aug` | `generate_samples(...).{theta,summary_min,summary_aug}` | Yes — column-major matrices from the keyed generator | FLOWING |
| `loader.jl load_fold` | `Ztr/Zva` | `load_main_pool(dir)` → `hcat` of all `shard_*.jld2` | Yes — reads real shard files; ZScoreTransform on real train data | FLOWING |

---

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| Full Phase-3 test suite (68 assertions) | `julia --project=spike spike/test/runtests.jl` | EXIT 0; 68/68 Phase-3 pass; 8 smoke + 62 simulator tests also pass | PASS |
| `encode_d01` 128-dim output + mask semantics | inline Julia one-liner | length=128; rows 1:64==vec(M); rows 65:128 all 1.0; missing→val=0/mask=0; fully-missing→128-zeros (not dropped) | PASS |
| D-08: Loader exports do not include `standardize_all` | `julia --project=spike -e 'include("spike/data/loader.jl"); println(names(Loader))'` | `[Loader, all_folds, load_fold, load_holdout, load_main_pool]`; `:standardize_all in names(Loader) = false` | PASS |
| HARD GATE: decoupling baseline | `git diff --quiet f581d95 -- src Project.toml Manifest.toml` | EXIT 0 (CLEAN) | PASS |

---

### Probe Execution

No conventional probe scripts (`scripts/*/tests/probe-*.sh`) exist for this phase. The equivalent gate is `julia --project=spike spike/test/runtests.jl` (run above, EXIT 0).

---

### Requirements Coverage

| Requirement | Source Plans | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| DATA-01 | Plans 02 | θ~π → simulate_pair → existing summary functions → fixed 128-dim standardized vector | SATISFIED | `encode.jl` (128-dim D-01 + AUG_DIM=142 D-02), `seeding.jl` (Philox4x D-11), `generate.jl` (thread-safe D-12); SC-1 + D-11/D-12 testsets pass (17 assertions) |
| DATA-02 | Plans 01, 03, 05 | 50k–200k pairs in a version-guarded JLD2 cache with Random123 loader | SATISFIED | `hashguard.jl` (D-05 SHA-256), `cache.jl` (D-04 atomic shards, resume), `generate_cache` driver (D-06 N-parameterized, D-10 holdout); SC-2 + hashguard + fixture testsets pass (26 assertions); NOTES.md §6 records real 50k run; at-scale cache dir exists (52 files) |
| DATA-03 | Plan 04 | Leak-free k-fold CV: fit-on-train-only ZScoreTransform, no global-standardize, ≥20-stack holdout excluded | SATISFIED | `loader.jl` (`module Loader`, `load_fold` D-07, no `standardize_all` D-08, `fold_rng` D-09, `_shard_files` glob D-10); SC-3 + SC-4 testsets pass (25 assertions) |

All three v1 requirements mapped to Phase 3 in REQUIREMENTS.md are SATISFIED. No orphaned requirements found (traceability table has 0 unmapped Phase-3 IDs).

---

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `spike/test/test_data_pipeline.jl` | 24 | "skipped placeholders" in comment | Info | Dead-code documentation — describes the Wave-0 design; all 7 testsets are now fully implemented (0 `@test_skip` in the executed paths) |
| `spike/test/test_data_pipeline.jl` | 271 | `@test false  # RED: spike/data/hashguard.jl not yet implemented` | Info | Inside `if !(@isdefined cache_hash)` branch — this branch is NEVER reached at runtime because `hashguard.jl` IS defined via `generate.jl` → `cache.jl` → `hashguard.jl`; 5/5 hashguard assertions pass |

No TBD / FIXME / XXX debt markers in any Phase-3 source file. No stubs, no placeholder returns, no hardcoded empty data flowing to rendered output.

---

### Human Verification Required

None. All must-haves are verifiable from the codebase and the one checkpoint:human-verify task (Plan 05) is resolved by:
(a) the at-scale cache directory on disk (`spike/data/cache/eca37f54.../`, 52 files confirming the 50-shard + holdout + meta layout), and
(b) NOTES.md §6 recording specific, detailed gate results (N=50000, 18.7 min, reload=50000 cols, resume=only deleted shard, cross-process -t1 vs -t auto byte-identical, EXIT 0 68/68, git diff clean).

---

### Gaps Summary

None. All 14 must-have truths are VERIFIED. The phase goal is fully achieved.

- DATA-01 (encoding + generator): `encode.jl` / `seeding.jl` / `generate.jl` implement all decision nodes (D-01 through D-03, D-11 through D-13); tested by SC-1 + D-11/D-12 (17 green assertions).
- DATA-02 (cache + persistence): `hashguard.jl` / `cache.jl` / `generate_cache` driver implement the full resumable, version-guarded, sharded pipeline (D-04 through D-06, D-10); tested by SC-2 + hashguard + fixture (26 green assertions) + confirmed at 50k scale on disk.
- DATA-03 (loader + split discipline): `loader.jl` (`module Loader`) implements the sole standardization path with structural cross-fold-leakage prevention (D-07 through D-10); tested by SC-3 + SC-4 (25 green assertions).
- HARD GATE (decoupling): `git diff --quiet f581d95 -- src Project.toml Manifest.toml` → clean.

---

_Verified: 2026-06-30_
_Verifier: Claude (gsd-verifier)_
