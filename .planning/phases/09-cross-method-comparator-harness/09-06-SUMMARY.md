---
phase: 09-cross-method-comparator-harness
plan: 06
subsystem: comparator
tags: [seeded-entry-point, reproducibility, content-addressed-artifact, tapqir-graceful-skip, fixture-oracle, traffic-light, determinism-gate]

# Dependency graph
requires:
  - phase: 09-01
    provides: Wave-0 @test_skip scaffold + pre-declared thresholds (config.jl) + runtests.jl wiring
  - phase: 09-02
    provides: classical estimator battery (manders, pearson_whole, spearman_whole, costes_p)
  - phase: 09-03
    provides: seeded regime-labelled shared-input builder (build_shared_inputs)
  - phase: 09-04
    provides: optional Tapqir anchor (tapqir_anchor) + isolated tapqir_env
  - phase: 09-05
    provides: build_table / write_table / comparator_audit + CMP_SRC_FILES content hash
provides:
  - Single seeded entry point (run_comparator) composing inputs+estimators+table+audit+optional Tapqir into one bit-reproducible artifact
  - Filled D-13 acceptance gate proving SC1 (finite per-method table), SC2 (Tapqir anchor-or-skip), SC3 (seeded/reproducible)
affects: [phase-10 related-work delta, phase-16 blind external eval]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Guarded-include composition (contract.jl FIRST) mirroring spike/data/generate.jl — no second src/ include"
    - "Optional external anchor stored as run meta, self-caught + defence-in-depth try/catch so a skip never fails the run"
    - "Determinism gate uses run_tapqir=false to stay subprocess-free + fully seeded; full Costes N on a tiny 3-input set for speed"
    - "Independent hand-built fixture oracles (encode_aug M1/M2 slice; known divergence→color; known regime-vs-verdict mismatches) — never re-derive the function's own formula"

key-files:
  created:
    - spike/comparator/run_comparator.jl
  modified:
    - spike/test/test_comparator.jl

key-decisions:
  - "run_comparator seeds build_table off its own master_seed kwarg (concrete even for the D-03 passthrough, whose shared.master_seed is missing) so Costes_p stays reproducible either way"
  - "Content-hash config = (master_seed, n_inputs, tapqir_status) — all deterministic — so two same-seed runs resolve to the identical artifact dir; the timestamped audit.md is written OUTSIDE the hash inputs and never perturbs determinism"
  - "Determinism gate compares the reloaded JLD2 table with isequal (bit-exact incl. missing NPE columns) plus == on the pure-Float Costes_p column, rather than raw file bytes"

patterns-established:
  - "One seeded entry point returns a NamedTuple (table, artifact_dir, audit, tapqir) — the callable deliverable"
  - "Reduced Costes n kwarg in the direct-estimator gates; tiny seeded inputs elsewhere keep the full suite fast"

requirements-completed: [CMP-08, CMP-09]

# Metrics
duration: ~15min
completed: 2026-07-02
---

# Phase 9 Plan 06: Seeded run_comparator Entry Point + Filled D-13 Gate

**A single seeded `run_comparator` entry point composes the frozen chain + the new comparator modules into one bit-reproducible comparison artifact + audit with a non-blocking optional Tapqir anchor, and the Wave-0 `@test_skip` scaffold is replaced with the real D-13 gate (45 assertions) proving SC1/SC2/SC3 green in the single spike test entry point.**

## Performance

- **Duration:** ~15 min
- **Completed:** 2026-07-02
- **Tasks:** 2
- **Files modified:** 2 (1 created, 1 modified)

## Accomplishments

- `run_comparator(; master_seed, n_per_regime, imsize, inputs, outdir, run_tapqir)` composes `build_shared_inputs` → `build_table` → `tapqir_anchor` → `write_table` → `comparator_audit` behind one call, returning `(table, artifact_dir, audit, tapqir)`.
- Bit-reproducibility proven: two same-seed runs resolve to the identical content-hash `artifact_dir` and a byte-identical JLD2 table payload (the Task-1 verify printed `OK dir=39ae7bfb… tapqir=skipped`; the CMP-05/08 gate asserts `isequal` on the reloaded table + `==` on `Costes_p`).
- The optional Tapqir anchor is genuinely non-blocking: `tapqir_anchor` self-catches, `run_comparator` adds a defence-in-depth try/catch, and `run_tapqir=false` skips the out-of-process probe entirely. On this machine the documented clean-skip path returns `:skipped` (no materialized Python), and the run still ships a full green classical table.
- The D-13 gate replaces all seven `@test_skip` placeholders with 45 real assertions using **independent hand-built oracles**: CMP-02 checks `manders(mci) == encode_aug(mci, patch_summary(mci))[129:130]` with `encode.jl` untouched; CMP-04 pins traffic-light colors at the exact pre-declared band boundaries and the divergence semantics against known regime-vs-verdict mismatches (random+high-Pearson→red, true coloc→green, missed coloc→red, `:unknown`→`missing`).
- `run_comparator.jl` reaches the frozen `src/` math only transitively through `contract.jl` — it contains **no** `include(` of `src/`.

## Task Commits

Each task was committed atomically:

1. **Task 1: Seeded `run_comparator` entry point (D-12, CMP-08)** — `78f934f` (feat)
2. **Task 2: Fill the D-13 comparator test gate (CMP-09)** — `cc2114d` (test)

## Files Created/Modified

- `spike/comparator/run_comparator.jl` (created) — guarded-include composition (contract.jl first) + `run_comparator` single seeded entry point returning `(table, artifact_dir, audit, tapqir)`.
- `spike/test/test_comparator.jl` (modified) — replaced the Wave-0 `@test_skip` scaffold with the filled CMP-01..08 gates (shared seeded fixtures + hand-built oracles).

## Verification

- **Task 1 verify:** `run_comparator(n_per_regime=2, imsize=(96,96))` twice → identical `artifact_dir` + `Costes_p`, `tapqir.status=:skipped` → printed `OK dir=… tapqir=skipped`.
- **Task 2 acceptance gate:** `julia --project=spike spike/test/runtests.jl` exits **0**; the full suite is green including the Phase-9 resolve-risk gate (g) and the concurrent Phase-5 testsets. Phase-9 comparator testset: **45/45 pass** (CMP-01 10, CMP-02 5, CMP-03 7, CMP-04 11, CMP-05/08 4, CMP-06 2, CMP-07 6).
- **Decoupling:** `git diff spike/data/encode.jl`, `git diff -- src/`, and `git diff spike/test/runtests.jl` are all empty (byte-unchanged).

## Decisions Made

- `run_comparator` passes its own `master_seed` kwarg to `build_table` so the Costes null is seeded even for the D-03 passthrough (whose `shared.master_seed` is `missing`) — otherwise `costes_rng(missing, …)` would throw on `UInt64(missing)`.
- The content-hash config is `(master_seed, n_inputs, tapqir_status)` (all deterministic), matching the plan's "master_seed + input count + Tapqir status"; the timestamped `audit.md` is written to the artifact dir *outside* the hash inputs, so it never perturbs artifact determinism.
- The determinism gate compares the reloaded JLD2 `table` via `isequal` (bit-exact, handles the `missing` NPE columns) plus `==` on the pure-Float `Costes_p` column, rather than raw JLD2 file bytes — a robust "byte-identical payload" check.
- The CMP-05/08 and CMP-07 gates use tiny seeded inputs (`n_per_regime=1`, `imsize=(64,64)`) and `run_tapqir=false` to keep the full suite fast (comparator testset ran in ~3.7s); the single real Tapqir probe lives in CMP-06.

## Deviations from Plan

None — plan executed exactly as written. Both tasks' verify commands passed as authored (Task 1 `OK dir=… tapqir=skipped`; Task 2 full suite exit 0).

## Issues Encountered

- The isolated `spike/comparator/tapqir_env/` Manifest is fully resolved (PythonCall/CondaPkg present), so `tapqir_anchor` spawns an out-of-process probe rather than short-circuiting on a missing Manifest. Handled by (a) exercising the real probe only once (CMP-06) and (b) using `run_tapqir=false` in the determinism gate — the probe warms once, returns `:skipped`, and the comparator testset stays at ~1.4s for CMP-06.

## Shared-Environment Isolation

- Staged only `spike/comparator/run_comparator.jl` and `spike/test/test_comparator.jl` by explicit path (never `git add -A`); the concurrent Phase-5 agent's uncommitted work in `spike/baseline/`, `spike/validation/`, and `.planning/` was left untouched.
- No package added/removed/resolved in the shared `spike/` env; `spike/data/encode.jl`, `src/`, and `spike/test/runtests.jl` are byte-unchanged.
- STATE.md / ROADMAP.md intentionally NOT updated (orchestrator-owned).

## Next Phase Readiness

- The callable deliverable (`run_comparator`) + the green D-13 gate close Phase 9's SC1/SC2/SC3. Phase 10 (related-work delta) and Phase 16 (blind external eval) can consume the seeded, content-addressed comparison artifact + audit without redesign.

## Self-Check: PASSED

- FOUND: spike/comparator/run_comparator.jl
- FOUND: spike/test/test_comparator.jl (no `@test_skip` remaining)
- FOUND commit: 78f934f (Task 1)
- FOUND commit: cc2114d (Task 2)
- VERIFIED: `julia --project=spike spike/test/runtests.jl` exits 0 (Phase-9 45/45)
- VERIFIED: spike/data/encode.jl, src/, spike/test/runtests.jl byte-unchanged (git diff empty)

---
*Phase: 09-cross-method-comparator-harness*
*Completed: 2026-07-02*
