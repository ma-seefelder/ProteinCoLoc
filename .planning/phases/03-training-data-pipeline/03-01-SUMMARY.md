---
phase: 03-training-data-pipeline
plan: 01
subsystem: infra
tags: [julia, jld2, random123, pkg-resolve, testing, spike, decoupling]

# Dependency graph
requires:
  - phase: 01-env-smoke
    provides: isolated spike Project.toml/Manifest with NeuralEstimators v0.2.1 pinned + UUID-keyed resolve-risk gate in runtests.jl
  - phase: 02-simulator
    provides: contract.jl include-coupling boundary, test_simulator.jl harness structure, Wave-0 resolve+refreeze pattern
provides:
  - JLD2 (sharded-cache backend, D-04) declared in spike [deps]
  - Random123 (counter-based per-sample seeding, D-11) declared in spike [deps]
  - re-frozen spike/Manifest.toml with NeuralEstimators still pinned v0.2.1
  - resolve-risk gate extended to cover the two new deps (JLD2/Random123 present + NE v0.2.1 re-assert)
  - test_data_pipeline.jl Wave-0 scaffold (SC-1..SC-4 + D-11/D-12 placeholders) wired into the single runtests.jl gate
affects: [03-02 generator, 03-03 cache, 03-04 loader, 03-05 holdout, training-data-pipeline]

# Tech tracking
tech-stack:
  added: [JLD2 (already transitive, now direct), Random123 v1.7.1 (+ RandomNumbers v1.6.0)]
  patterns: [Wave-0 dependency-add + Manifest refreeze, MISSING test scaffold of named skipped SC testsets, UUID-keyed resolve-risk gate extension]

key-files:
  created:
    - spike/test/test_data_pipeline.jl
  modified:
    - spike/Project.toml
    - spike/Manifest.toml
    - spike/test/runtests.jl

key-decisions:
  - "JLD2 was already a transitive dependency of the spike env; the add only promoted it to a direct [deps] entry — no new resolve churn. Random123 v1.7.1 + RandomNumbers v1.6.0 were the only newly-installed nodes."
  - "Resolve-risk gate extended in place (within the existing D-04 testset, branch (e)) rather than a new testset, keeping a single coherent supply-chain/co-resolve guard keyed on the NeuralEstimators UUID."
  - "Wave-0 scaffold uses @test_skip true placeholders (reported as Broken, EXIT=0) so every SC is enumerated in the gate before any pipeline code exists; each placeholder is tagged with the wave (W2/W3) that fills it."

patterns-established:
  - "MISSING test scaffold: named child testsets with @test_skip true + wave tag, replaced by real gates in later waves"
  - "Dependency-add Wave-0: edit spike Project.toml [deps] alphabetically -> Pkg.resolve+precompile -> re-assert NE v0.2.1 -> extend resolve-risk gate before any code targets the new deps"

requirements-completed: [DATA-02]

# Metrics
duration: 8min
completed: 2026-06-28
---

# Phase 3 Plan 01: Pipeline Dependency + Test Scaffold Foundation Summary

**JLD2 + Random123 added to the isolated spike env (NeuralEstimators still pinned v0.2.1), resolve-risk gate extended to cover them, and a five-SC MISSING test scaffold wired into the single runtests.jl gate — all spike-local, src/ + root manifests byte-identical to baseline f581d95.**

## Performance

- **Duration:** ~8 min
- **Started:** 2026-06-28T09:24:55Z
- **Completed:** 2026-06-28T09:32:34Z
- **Tasks:** 2
- **Files modified:** 4 (1 created, 3 modified)

## Accomplishments
- Added `JLD2` (033835bb-…) and `Random123` (74087812-…) to `spike/Project.toml` `[deps]` in alphabetical position; re-froze `spike/Manifest.toml` via `Pkg.resolve()` + `Pkg.precompile()` (only Random123 v1.7.1 + RandomNumbers v1.6.0 newly installed; JLD2 was already transitive).
- Confirmed the Phase-1 GLMakie-class co-resolve regression did NOT recur: NeuralEstimators stays pinned at v0.2.1 after the new deps resolve.
- Extended the D-04 resolve-risk testset with branch (e): asserts JLD2 + Random123 are present in the resolved set and re-asserts NeuralEstimators v0.2.1 with them resolved (now 7 asserts).
- Created `spike/test/test_data_pipeline.jl` — a Wave-0 MISSING scaffold with five named child testsets (SC-1 generator→128-vector, SC-2 cache round-trip/resume/invalidation, SC-3 leak-free standardization, SC-4 k-fold+holdout, D-11/D-12 thread independence), each a `@test_skip true` placeholder tagged with the wave that fills it.
- Wired `include(test_data_pipeline.jl)` into `runtests.jl` so a single `julia --project=spike spike/test/runtests.jl` remains the gate; full harness green (EXIT=0).

## Task Commits

Each task was committed atomically:

1. **Task 1: Add JLD2 + Random123 to the spike env and re-freeze the Manifest** - `e06414e` (chore)
2. **Task 2: Extend the resolve-risk gate and scaffold test_data_pipeline.jl wired into runtests.jl** - `be26006` (test)

**Plan metadata:** committed separately (docs: complete plan)

## Files Created/Modified
- `spike/Project.toml` - added JLD2 + Random123 UUID lines to `[deps]` (alphabetical)
- `spike/Manifest.toml` - re-frozen: +Random123 v1.7.1, +RandomNumbers v1.6.0; NeuralEstimators unchanged at v0.2.1
- `spike/test/runtests.jl` - extended D-04 resolve-risk testset (branch (e), JLD2/Random123 present + NE re-assert) and added `include(test_data_pipeline.jl)`
- `spike/test/test_data_pipeline.jl` - new Wave-0 MISSING scaffold: five named SC testsets with skipped placeholders

## Decisions Made
- Edited the resolve-risk gate in place inside the existing D-04 testset (branch (e)) rather than creating a parallel testset — keeps one coherent UUID-keyed supply-chain guard.
- Kept JLD2's promotion-to-direct as the only Project.toml change for it (it was already resolved transitively), avoiding any Manifest churn beyond Random123/RandomNumbers.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None. The resolve completed cleanly; the only newly-installed packages were Random123 v1.7.1 and its dependency RandomNumbers v1.6.0 (JLD2 was already present transitively in the spike Manifest).

## User Setup Required
None - no external service configuration required. (Pkg resolve/precompile ran offline against the already-populated Julia depot.)

## Next Phase Readiness
- JLD2 (cache backend, D-04) and Random123 (counter-based seeding, D-11) are now loadable in the spike env — Wave 2 (`generate.jl`, `seeding.jl`, `encode.jl`) and Wave 3 (`cache.jl`, `hashguard.jl`, `loader.jl`) can build against them.
- The five SC placeholders are enumerated in the gate; later waves replace each `@test_skip true` with its real gate.
- Decoupling preserved: `git diff --quiet f581d95 -- src Project.toml Manifest.toml` is clean.

## Self-Check: PASSED

- FOUND: spike/test/test_data_pipeline.jl
- FOUND: .planning/phases/03-training-data-pipeline/03-01-SUMMARY.md
- FOUND commit: e06414e (Task 1)
- FOUND commit: be26006 (Task 2)

---
*Phase: 03-training-data-pipeline*
*Completed: 2026-06-28*
