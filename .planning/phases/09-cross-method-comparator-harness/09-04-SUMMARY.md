---
phase: 09-cross-method-comparator-harness
plan: 04
subsystem: testing
tags: [tapqir, pythoncall, condapkg, cosmos, sbi, graceful-degradation, subprocess-isolation, decoupling]

# Dependency graph
requires:
  - phase: 09-01
    provides: "spike/comparator/config.jl with TAPQIR_PUBLISHED_VALUE / TAPQIR_TOL placeholders and the runtests.jl gate (g) resolve-risk assertions"
  - phase: 01
    provides: "CPU/CUDA graceful-degradation idiom (optional heavy stack never hard-fails the core path) mirrored here for the Python/Tapqir edge"
provides:
  - "Isolated Tapqir sub-env (spike/comparator/tapqir_env/) pinning PythonCall + CondaPkg in its OWN Project/Manifest — never in the main spike env"
  - "spike/comparator/tapqir_bridge.jl: out-of-process, subprocess-only Tapqir anchor with graceful skip-with-flag (tapqir_anchor / _tapqir_env_available)"
  - "config.jl Tapqir anchor resolved as a documented NaN clean-skip (no Conda backend on this machine)"
affects: [09-06, 09-05, phase-10, phase-16]

# Tech tracking
tech-stack:
  added: ["PythonCall 0.9 (isolated sub-env only)", "CondaPkg 0.2 (isolated sub-env only)"]
  patterns:
    - "Isolated Python sub-env mirroring spike/baseline/ (own Project.toml/Manifest.toml)"
    - "Subprocess-only external bridge: no in-process PythonCall on the classical path (OQ2)"
    - "Hang-proof subprocess (poll+kill timeout) + CondaPkg offline probe (no network install)"

key-files:
  created:
    - "spike/comparator/tapqir_env/Project.toml"
    - "spike/comparator/tapqir_env/Manifest.toml"
    - "spike/comparator/tapqir_bridge.jl"
  modified:
    - "spike/comparator/config.jl"

key-decisions:
  - "Take the documented clean-skip (TAPQIR_PUBLISHED_VALUE = NaN): no Conda backend on this Windows machine, so per D-08 the optional anchor skips-with-flag rather than blocking the phase"
  - "Subprocess-only Tapqir invocation under the isolated env (OQ2) — PythonCall is NEVER imported in the main harness process"
  - "Corrected the CondaPkg UUID (the plan's value was unregistered; registry canonical UUID used instead)"

patterns-established:
  - "Optional external bridge with skip-with-flag (RESEARCH Pattern 4): probe availability first, map ANY failure to status=:skipped"
  - "CondaPkg offline + Null-backend probe so an unmaterialized env reports unavailable instead of triggering a hanging download"

requirements-completed: [CMP-06]

# Metrics
duration: 45min
completed: 2026-07-02
---

# Phase 9 Plan 04: Tapqir Sanity Anchor (Isolated, Graceful Skip-with-Flag) Summary

**An out-of-process, strictly isolated Tapqir bridge (PythonCall/CondaPkg in their own sub-env) that degrades to a clean `:skipped` anchor with zero Python present, plus a documented NaN clean-skip for the published anchor value on this Conda-less Windows machine.**

## Performance

- **Duration:** ~45 min
- **Started:** 2026-07-02T09:57:46Z
- **Completed:** 2026-07-02T~10:40Z
- **Tasks:** 3 (2 auto-executed + 1 checkpoint auto-approved)
- **Files modified:** 4

## Accomplishments
- Stood up `spike/comparator/tapqir_env/` — an isolated sub-env (own `Project.toml` + resolved/committed `Manifest.toml`) pinning `PythonCall` + `CondaPkg`, mirroring the `spike/baseline/` decoupling precedent. These deps are provably ABSENT from the main `spike/Project.toml` (runtests.jl gate (g) already asserts this; Task 1 verify re-confirmed via `Pkg.project().dependencies`).
- Implemented `spike/comparator/tapqir_bridge.jl`: `_tapqir_env_available()` probes the sub-env + a subprocess PythonCall/tapqir import and returns `false` on ANY error (never throws); `tapqir_anchor(; tol)` returns `(status, value, reason)` and maps every failure to a clean `:skipped`. Tapqir runs strictly OUT-OF-PROCESS in a subprocess activated against the isolated env — no in-process PythonCall on the classical path (OQ2 RESOLVED). A poll+kill timeout and CondaPkg-offline probe make the bridge hang-proof (no network/Conda install can stall the core path).
- Resolved the Tapqir anchor value as a documented clean-skip: no Conda/mamba/micromamba backend is available on this machine (system Python is 3.14, unsupported by tapqir 1.1.19), so `TAPQIR_PUBLISHED_VALUE` stays `NaN` with a comment recording the eLife-2022 Part II `cosmos` dataset + how to pin a real value later. `TAPQIR_TOL = 0.05` retained.
- Verified: `tapqir_anchor()` returns `status=:skipped` (reason "Tapqir/Python env unavailable") without throwing; `PythonCall`/`CondaPkg` absent from the main env; anchor value `NaN` accepted by both verify one-liners.

## Task Commits

1. **Task 1: Isolated Tapqir sub-env + graceful skip-with-flag bridge** - `7e8ce9c` (feat)
2. **Task 2: Document Tapqir anchor clean-skip** - `45a3a6d` (docs)
3. **Task 3: Verify anchor value / documented clean-skip (checkpoint:human-verify)** - AUTO-APPROVED (no separate commit; see below)

**Plan metadata:** committed with this SUMMARY (docs)

## Files Created/Modified
- `spike/comparator/tapqir_env/Project.toml` - Isolated PythonCall + CondaPkg sub-env spec ([compat]-pinned)
- `spike/comparator/tapqir_env/Manifest.toml` - Resolved+precompiled reproducibility artifact for the sub-env
- `spike/comparator/tapqir_bridge.jl` - Out-of-process, subprocess-only Tapqir anchor with skip-with-flag
- `spike/comparator/config.jl` - Tapqir anchor resolved as a documented NaN clean-skip (comment updated; consts unchanged in value)

## Decisions Made
- **Documented clean-skip over a materialization attempt.** Per the autonomous-execution directive and locked decision D-08, and given no Conda backend + an incompatible system Python, materializing tapqir 1.1.19 would have been a long/likely-failing network install. Left `TAPQIR_PUBLISHED_VALUE = NaN` (documented, not a guessed number — D-14 discipline) so the 09-06 D-13 gate accepts `status=:skipped`.
- **Subprocess-only bridge (OQ2).** Chose the RESEARCH-recommended subprocess invocation under the isolated env over in-process PythonCall, so no Python ever loads on the classical battery path.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Corrected the CondaPkg package UUID**
- **Found during:** Task 1 (resolving the isolated sub-env)
- **Issue:** The plan specified CondaPkg UUID `992eb4ea-22a4-4c86-a8db-e01b23e3f0c6`, which is not registered in the General registry — `Pkg.instantiate()` failed with a `check_registered` error that reported the canonical UUID.
- **Fix:** Used the registry-canonical CondaPkg UUID `992eb4ea-22a4-4c89-a5bb-47a3300528ab` in `tapqir_env/Project.toml`. (Does not affect the main-env absence gate, which is keyed by NAME.)
- **Files modified:** spike/comparator/tapqir_env/Project.toml
- **Verification:** `Pkg.instantiate()` resolved + precompiled cleanly (PythonCall/CondaPkg + deps); Manifest committed.
- **Committed in:** `7e8ce9c` (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** The UUID correction was required for the isolated env to resolve at all; no scope creep. All other work followed the plan.

## Issues Encountered
- The full `julia --project=spike spike/test/runtests.jl` gate (checkpoint how-to-verify step 1) runs the NPE smoke train and is slow; it was started but abandoned per instruction in favor of the direct, quick invariant checks. The two auto-approval invariants were confirmed independently: (a) Tasks 1-2 verify one-liners pass, and (b) `PythonCall`/`CondaPkg` are absent from the main `spike/Project.toml` (asserted in the Task 1 verify) and `tapqir_anchor()` returns a clean `:skipped` without throwing. Gate (g) that lives inside runtests.jl already encodes the same PythonCall/CondaPkg-absence assertion.

## Task 3 Checkpoint (auto-approved)
Task 3 is a `checkpoint:human-verify` gate. This is an autonomous background run with no human available, so the checkpoint was **auto-approved** per the execution directive. The auto-approval preconditions were met: Tasks 1-2 passed their automated `<verify>` commands, PythonCall/CondaPkg are absent from the main spike env, and `tapqir_anchor()` returns a clean status without throwing. **Path taken: the documented clean-skip** (no Conda/Tapqir env on this Windows machine) — `TAPQIR_PUBLISHED_VALUE = NaN`, `tapqir_anchor() -> status=:skipped`. To capture a real anchor value later, materialize `spike/comparator/tapqir_env/` on a Conda-capable machine and re-run Task 2.

## User Setup Required
Optional only. To report the Tapqir anchor in the manuscript, install a Conda backend (set `JULIA_CONDAPKG_BACKEND`), materialize `spike/comparator/tapqir_env/` (installs tapqir 1.1.19), run the eLife-2022 Part II `cosmos` tutorial once to capture the recovered scalar, and pin it as `TAPQIR_PUBLISHED_VALUE`/`TAPQIR_TOL` in `config.jl`. The classical battery runs green without any of this.

## Next Phase Readiness
- CMP-06 satisfied: the Tapqir bridge is fully isolated and gracefully degrades; the classical battery is unaffected by the optional Python edge.
- 09-06 D-13 gate can assert `status ∈ {:passed, :skipped}` — on this machine it will observe `:skipped`.
- No blockers. STATE.md / ROADMAP.md intentionally NOT modified (orchestrator owns those writes in worktree mode).

## Self-Check: PASSED
- All 5 created/modified files present on disk.
- All 3 commits present in git history (`7e8ce9c`, `45a3a6d`, `60d1734`).
- PythonCall/CondaPkg confirmed ABSENT from the committed main `spike/Project.toml`.

---
*Phase: 09-cross-method-comparator-harness*
*Completed: 2026-07-02*
