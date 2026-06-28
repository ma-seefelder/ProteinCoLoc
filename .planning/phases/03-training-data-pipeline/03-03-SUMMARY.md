---
phase: 03-training-data-pipeline
plan: 03
subsystem: database
tags: [jld2, sha256, sharded-cache, content-hash, random123, holdout, atomic-write]

# Dependency graph
requires:
  - phase: 03-02
    provides: generate_sample/generate_samples column-major (theta, summary_min[128], summary_aug[AUG_DIM]); seeding salts (HOLDOUT_SALT) + holdout_rng; encode_d01/encode_aug
  - phase: 03-01
    provides: JLD2 + Random123 added to isolated spike env; resolve-risk gate; MISSING test_data_pipeline.jl scaffold
provides:
  - "hashguard.jl: SHA-256 content hash over data-defining source bytes + canonical(config), with per-file sub-hashes (D-05)"
  - "cache.jl: atomic sharded JLD2 writes, resume-by-skip, content-hash-named dir + meta.jld2, open_or_invalidate (D-04)"
  - "generate_cache driver: Threads.@threads over shards, N-parameterized scale-up, reserved disjoint-namespace holdout (D-06/D-10)"
  - "SC-2 round-trip/resume/invalidation gate green; committed cache fixture pinning the on-disk schema"
affects: [03-04-loader, phase-04-training]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Content-hash version guard: sha256(source bytes + canonical config) names the cache dir and is re-checked on reopen"
    - "Atomic shard write: jldsave .tmp -> reopen-integrity-check -> mv(force=true)"
    - "Resume-by-skip: shard_done (isfile && reload-ok) -> skip completed shards, regenerate only the in-flight one"
    - "Disjoint-namespace holdout: XOR-salted RNG + separate file + negative global_index (structural exclusion)"

key-files:
  created:
    - spike/data/hashguard.jl
    - spike/data/cache.jl
    - spike/data/cache/fixture/shard_0001.jld2
  modified:
    - spike/data/generate.jl
    - spike/test/test_data_pipeline.jl
    - .gitignore

key-decisions:
  - "Cache dir NAMED by D-05 content hash; open_or_invalidate's mismatch-error is the defensive tripwire for a corrupted/tampered manifest (versions otherwise auto-separate by dir name)"
  - "Holdout stored global_index is NEGATIVE (-1,-2,...) so it shares no value with the positive 1:N main pool — structural disjointness, satisfying Task-2 'shares NO global_index' without contradicting the loader's never-read-holdout exclusion"
  - "Holdout is also resume-skipped (shard_done on holdout.jld2) so re-runs don't re-simulate the >=20 reserved stacks"
  - "Frozen src/colocalization.jl + src/LoadImages.jl hashed into the D-05 digest as a decoupling-breach tripwire"
  - "Committed fixture written via direct write_shard (no driver/holdout) to keep it tiny and at a stable non-hash-named path for direct-reload schema-stability assertion"

patterns-established:
  - "Pattern: generating_config NamedTuple centralizes every D-05-hashed field (N, shard_size, master_seed, k, theta_dim, imsize-spec, summary-spec, schema_version); canonical() sorts by name so field order is irrelevant"
  - "Pattern: guarded sibling includes (isdefined || include) keep generate.jl->cache.jl->hashguard.jl idempotent under standalone + runtests loading"

requirements-completed: [DATA-02]

# Metrics
duration: 26min
completed: 2026-06-28
---

# Phase 3 Plan 3: DATA-02 Sharded Version-Guarded JLD2 Cache Summary

**Resumable, content-hash-versioned sharded JLD2 cache: atomic .tmp+mv shard writes with resume-by-skip, a SHA-256 source+config guard that auto-invalidates stale data, and a reserved disjoint-namespace ADVI holdout.**

## Performance

- **Duration:** ~26 min
- **Started:** 2026-06-28
- **Completed:** 2026-06-28
- **Tasks:** 3 (Task 1 TDD: 2 commits)
- **Files modified:** 6 (3 created, 3 modified)

## Accomplishments
- `hashguard.jl` — deterministic SHA-256 over FIXED-order data-defining source bytes + `canonical(config)`, with per-file `subhashes` so a mismatch names exactly which source changed (D-05)
- `cache.jl` — atomic sharded JLD2 writes (`.tmp` → reopen-integrity-check → `mv force=true`), `shard_done` resume predicate, content-hash-named dir + `meta.jld2`, and `open_or_invalidate` that errors on a stored-vs-recomputed hash mismatch (D-04 / T-03-06 / T-03-07)
- `generate_cache` driver in `generate.jl` — `Threads.@threads` over shards with per-global-index keying (byte-identical for any thread count), N-parameterized scale-up that appends shards without regeneration (D-06), and a reserved ≥20-stack holdout drawn from the disjoint `holdout_rng` namespace into a separate `holdout.jld2` (D-10)
- SC-2 round-trip / resume-by-skip / invalidation gate is green and no longer skipped; a tiny committed fixture pins the column-major on-disk schema; bulk cache gitignored

## Task Commits

Each task was committed atomically:

1. **Task 1 (TDD RED): failing hashguard test** - `c497c93` (test)
2. **Task 1 (TDD GREEN): hashguard content-hash guard** - `053e43c` (feat)
3. **Task 2: cache.jl + sharded driver + holdout + .gitignore** - `58c3cb6` (feat)
4. **Task 3: SC-2 tests + committed fixture** - `44c9c88` (test)

**Plan metadata:** (this commit) (docs: complete plan)

## Files Created/Modified
- `spike/data/hashguard.jl` (created) - `cache_hash`/`subhashes`/`canonical` + `HASH_SRC_FILES`; SHA-256 D-05 content hash
- `spike/data/cache.jl` (created) - `shard_path`/`write_shard`/`shard_done`/`write_meta`/`open_or_invalidate`; `SHARD_SIZE`/`SCHEMA_VERSION`
- `spike/data/cache/fixture/shard_0001.jld2` (created) - tiny committed schema fixture (4 samples, seed 20240628)
- `spike/data/generate.jl` (modified) - `generating_config`, `_write_holdout`, `generate_cache` driver; guarded `include` of cache.jl
- `spike/test/test_data_pipeline.jl` (modified) - hashguard testset (Task 1) + SC-2 round-trip/resume/invalidation + fixture testset (Task 3)
- `.gitignore` (modified) - ignore `spike/data/cache/*` with `!spike/data/cache/fixture/` negation

## Decisions Made
- Cache dir is content-hash-named, so changing config/source resolves to a *different* dir (auto-separation); `open_or_invalidate`'s error path is the defensive guard against a corrupted/tampered `meta.jld2` — SC-2c exercises it by rewriting the stored hash.
- Holdout `global_index` stored as negative to make disjointness from the 1:N main pool structural and value-level (satisfies Task-2 acceptance "shares NO global_index").
- Holdout draws are resume-skipped on an existing `holdout.jld2` to avoid re-simulating the reserved stacks on every driver call.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Hardened resume + holdout disjointness beyond the literal plan text**
- **Found during:** Task 2
- **Issue:** (a) The plan's holdout uses the natural 1:H index, but Task-2 acceptance requires "shares NO global_index with any main-pool shard"; (b) resume-by-skip only covered shards, so each driver call would re-simulate the ≥20 holdout stacks.
- **Fix:** Store holdout `global_index` as negative (`-j`) for value-level disjointness; gate holdout generation on `_loads_ok(holdout.jld2)` so it is also resume-skipped.
- **Files modified:** spike/data/generate.jl
- **Verification:** SC-2 asserts `intersect(pool.global_index, holdout.global_index) == ∅`; resume re-run regenerates only the deleted shard.
- **Committed in:** 58c3cb6 (Task 2 commit)

**2. [Rule 3 - Blocking] Added hashguard testset in Task 1 (touches test_data_pipeline.jl earlier than the plan's Task-3 file list)**
- **Found during:** Task 1 (TDD requires a test home)
- **Issue:** Task 1 is `tdd="true"` but its only listed file is `hashguard.jl`; a RED/GREEN cycle needs a committed test.
- **Fix:** Added a dedicated `hashguard content-hash version guard (D-05)` testset (covering per-file sub-hash diagnosability, which SC-2 does not) with a guarded include so the suite runs RED then GREEN. The SC-2 placeholder is still filled in Task 3 as planned.
- **Files modified:** spike/test/test_data_pipeline.jl
- **Verification:** Testset RED (file absent) → GREEN (after hashguard.jl); full suite green.
- **Committed in:** c497c93 (RED) / 053e43c (GREEN)

---

**Total deviations:** 2 auto-fixed (2 blocking)
**Impact on plan:** Both strengthen correctness/acceptance coverage. No scope creep; all files remain under `spike/` (+ the planned `.gitignore`).

## Issues Encountered
- `gsd-sdk query state.record-metric` / `add-decision` take named flags (`--phase`/`--summary`), not positional args — adjusted invocation.

## Known Stubs
None — the cache, driver, and holdout are fully wired and exercised by green gates. SC-3/SC-4 (loader) remain `@test_skip` placeholders owned by Plan 03-04 (Wave 4), not this plan.

## Threat Flags
None — no new trust boundary beyond the planned cache/holdout/content-hash surface already in `<threat_model>` (T-03-06/07/08 mitigated; T-03-09 accepted).

## Next Phase Readiness
- DATA-02 persistence is complete: a 50k→200k generation is now resumable and version-guarded. Plan 03-04 (loader, DATA-03) can consume the column-major shards + `holdout.jld2` and implement the leak-free k-fold standardization (SC-3/SC-4).
- Decoupling holds: `git diff --quiet f581d95 -- src Project.toml Manifest.toml` is clean.

## Self-Check: PASSED
- Files verified present: spike/data/hashguard.jl, spike/data/cache.jl, spike/data/generate.jl, spike/test/test_data_pipeline.jl, spike/data/cache/fixture/shard_0001.jld2, .gitignore
- Commits verified: c497c93, 053e43c, 58c3cb6, 44c9c88
- Gates: `julia --project=spike spike/test/test_data_pipeline.jl` green (SC-2 no longer skipped); `julia --project=spike spike/test/runtests.jl` exit 0; decoupling baseline clean

---
*Phase: 03-training-data-pipeline*
*Completed: 2026-06-28*
