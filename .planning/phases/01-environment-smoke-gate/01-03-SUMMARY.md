---
phase: 01-environment-smoke-gate
plan: 03
subsystem: infra
tags: [julia, juliaup, neuralestimators, flux, reproducibility, manifest, sbi]

# Dependency graph
requires:
  - phase: 01-02
    provides: Green CPU-only NeuralEstimators NPE smoke + resolved spike/Manifest.toml committed
provides:
  - Pinned, committed reproducibility artifact (spike/Manifest.toml, NeuralEstimators 0.2.1 + Flux 0.16.10, no top-level CUDA)
  - Enforced Julia-version pin (juliaup directory override spike/ -> 1.12.6) + documentation record (spike/.julia-version)
  - Stack-decision note (spike/NOTES.md): NeuralEstimators default / BayesFlow fallback-only + excluded packages
  - Phase-2 seed: Turing prior ranges (mu/nu/sigma/tau) recorded in NOTES.md
affects: [01-04, 02, phase-2-simulator, all-downstream-spike-phases]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Reproducibility freeze as committed artifact: pin Manifest + Julia version after a green gate"
    - "Dual Julia-version pin: juliaup override (enforcing) + .julia-version (documentation-only)"
    - "CPU-only assertion via no top-level [[deps.CUDA]] package (CUDA weakdep extensions are inert)"

key-files:
  created:
    - spike/.julia-version
    - spike/NOTES.md
  modified: []

key-decisions:
  - "Pinned exact Julia 1.12.6 as a juliaup channel (juliaup add 1.12.6) rather than the floating release channel, then set a spike/-scoped directory override"
  - "CUDA-absence asserted as 'no top-level [[deps.CUDA]] package' rather than a bare \\bCUDA\\b regex, since CUDA appears as inert weakdep extension declarations in the Manifest"
  - "spike/Manifest.toml left verbatim (already committed by 01-02, unchanged) — no re-resolve, preserving the proven env byte-for-byte"

patterns-established:
  - "Pattern: freeze a proven env (Manifest + Julia pin) only after the smoke gate is green"
  - "Pattern: enforcing pin (juliaup override) + human-readable record (.julia-version + NOTES.md)"

requirements-completed: [ENV-03, ENV-04]

# Metrics
duration: 9min
completed: 2026-06-26
---

# Phase 1 Plan 03: Reproducibility Freeze + Stack-Decision Note Summary

**Froze the proven CPU-only NeuralEstimators 0.2.1 / Flux 0.16.10 spike env as a committed reproducibility artifact — Julia 1.12.6 pinned via juliaup directory override + documentation record, and the NeuralEstimators-default / BayesFlow-fallback stack decision recorded in spike/NOTES.md.**

## Performance

- **Duration:** ~9 min (excludes 3m33s smoke re-run)
- **Started:** 2026-06-26
- **Completed:** 2026-06-26
- **Tasks:** 2
- **Files modified:** 2 created (spike/.julia-version, spike/NOTES.md); juliaup override registered (OS state)

## Accomplishments
- Re-confirmed the smoke gate GREEN before freezing (4/4 tests pass; recovered mu_hat=0.799 vs theta_true=0.7, |delta|=0.099 < tol=0.3) — only a proven env is worth pinning.
- Registered the **enforcing** Julia pin: `juliaup add 1.12.6` then `juliaup override set --path .../spike 1.12.6` (verified via `juliaup override status`).
- Wrote `spike/.julia-version` (literal `1.12.6`, documentation-only record).
- Confirmed `spike/Manifest.toml` (already committed verbatim by 01-02, unchanged) pins NeuralEstimators 0.2.1 + Flux 0.16.10 + Distributions 0.25.128, Manifest `julia_version` 1.12.6, and carries **no top-level `[[deps.CUDA]]` package**.
- Wrote `spike/NOTES.md`: stack decision (NeuralEstimators default / BayesFlow PythonCall fallback-only), the named exclusions (NormalizingFlows.jl, InvertibleNetworks.jl, RxInfer-as-Turing-replacement, BayesFlow-on-core-path, CUDA-as-requirement), the Julia-version record, and a Phase-2 seed of the Turing prior ranges (mu/nu/sigma/tau).
- Root `Project.toml`/`Manifest.toml`/`src/` verified byte-clean vs the Plan-01 baseline `f581d95` throughout.

## Task Commits

Each task was committed atomically:

1. **Task 1: Pin and commit the Manifest + record the Julia version** - `9a20327` (chore)
2. **Task 2: Write the stack-decision + Julia-version note** - `0283a46` (docs)

**Plan metadata:** _(final docs commit — see below)_

## Files Created/Modified
- `spike/.julia-version` - Single-line documentation record of the pinned Julia version (`1.12.6`); not auto-read by juliaup.
- `spike/NOTES.md` - Stack decision (NeuralEstimators default / BayesFlow fallback-only) + excluded packages + Julia-version pin record + Phase-2 Turing prior seed.
- _(unchanged, confirmed)_ `spike/Manifest.toml` - Pinned reproducibility artifact (NeuralEstimators 0.2.1, Flux 0.16.10), already committed verbatim by Plan 01-02; this plan is the safe known-good revert point before Plan 01-04's `Pkg.develop` attempt.

## Decisions Made
- **Exact-version channel, not floating `release`:** Pinned `1.12.6` as a dedicated juliaup channel via `juliaup add 1.12.6` (the machine default `release` channel currently resolves to 1.12.6 but could move), then scoped a directory override to `spike/`. This makes the pin survive future Julia releases.
- **Manifest left verbatim:** It was already resolved and committed by Plan 01-02 and is unchanged — no re-resolve was performed, preserving the proven environment byte-for-byte (re-resolving would risk version drift on a pre-1.0 stack).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Corrected the CUDA-absence verification check**
- **Found during:** Task 1 (Manifest pin verification)
- **Issue:** The plan's automated verify uses `!occursin(r"\bCUDA\b", m)` on the whole Manifest string. The resolved Manifest legitimately contains `CUDA` as **inert weakdep/extension declarations** under Flux/NNlib/Zygote/NeuralEstimators (e.g. `FluxCUDAExt = "CUDA"`, `CUDA = "052768ef-..."` in `[weakdeps]`). The bare-word regex matches those and would falsely report a CUDA dependency — the same issue STATE.md already recorded for Plan 01-02.
- **Fix:** Used the correct semantic check — assert there is **no top-level `[[deps.CUDA]]` installed package entry** (`!occursin(r"(?m)^\[\[deps\.CUDA\]\]", m)`). This is the true CPU-only guarantee; weakdep extensions install and load nothing on the CPU path.
- **Files modified:** None (verification-command correction only; no artifact changed).
- **Verification:** `julia -e` assertion confirms NeuralEstimators 0.2.1 + Flux 0.16.10 present and no top-level `[[deps.CUDA]]`; passes.
- **Committed in:** N/A (no file change — documented here for the verifier).

---

**Total deviations:** 1 auto-fixed (1 bug — verification-command correctness).
**Impact on plan:** No artifact or scope change. The corrected check is strictly more accurate and matches the precedent set in Plan 01-02. All plan acceptance criteria met.

## Issues Encountered
- `juliaup override set 1.12.6` initially failed with "'1.12.6' is not installed" because only the floating `release` channel was present. Resolved by `juliaup add 1.12.6` to materialize the exact-version channel, then setting the override. This is the correct enforcing-pin path and is documented in NOTES.md for reproduction on other machines.

## User Setup Required
None - no external service configuration required. (The juliaup directory override is OS-local state; NOTES.md documents how to reproduce it on another machine.)

## Next Phase Readiness
- ENV-03 and ENV-04 satisfied: the proven environment is a durable, committed, reproducible artifact (pinned Manifest + enforced Julia pin) and the stack ladder is recorded.
- Plan 01-04 (`Pkg.develop(path="..")` co-resolve attempt) can now proceed against a safe known-good revert point: if the develop resolve drags an incompatible/heavy tree, revert `spike/Manifest.toml` to this committed state and fall back to `include()` per D-01.
- No blockers.

## Self-Check: PASSED

- FOUND: spike/.julia-version
- FOUND: spike/NOTES.md
- FOUND: .planning/phases/01-environment-smoke-gate/01-03-SUMMARY.md
- FOUND commit: 9a20327 (Task 1)
- FOUND commit: 0283a46 (Task 2)

---
*Phase: 01-environment-smoke-gate*
*Completed: 2026-06-26*
