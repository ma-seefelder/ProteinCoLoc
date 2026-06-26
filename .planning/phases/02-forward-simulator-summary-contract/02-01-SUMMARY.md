---
phase: 02-forward-simulator-summary-contract
plan: 01
subsystem: testing
tags: [julia, neuralestimators, statsbase, images, imagefiltering, cairomakie, hypothesistests, summary-contract, sbi, decoupling]

# Dependency graph
requires:
  - phase: 01-environment-smoke-gate
    provides: "Pinned spike env (NeuralEstimators 0.2.1 / Flux 0.16.10, CPU-only, julia 1.12.6), include() coupling fallback (D-01), baseline f581d95, re-runnable runtests.jl smoke gate"
provides:
  - "Spike env extended with the simulator/summary stack: StatsBase, Images, ImageFiltering, CairoMakie, HypothesisTests — NeuralEstimators 0.2.1 intact, no CUDA"
  - "spike/contract.jl: read-only include() boundary to src/LoadImages.jl + src/colocalization.jl with build_mci/summary/induced_mu helpers"
  - "Frozen SIM-03 summary contract (fixed 8x8 per-patch Pearson rho) proven on a synthetic image — the known-good target Wave-2 simulator builds against"
  - "Automated resolve-risk gate (NeuralEstimators==v0.2.1) added to runtests.jl"
affects: [forward-simulator, prior, calibration, sbc, bayes-factor, advi-baseline]

# Tech tracking
tech-stack:
  added: [StatsBase, Images, ImageFiltering, CairoMakie, HypothesisTests]
  patterns: ["include() read-only coupling to frozen src/", "StatsBase-before-include ordering", "background-not-zero + >=15-px-floor contract traps as test assertions", "resolve-risk version-pin gate in test harness"]

key-files:
  created:
    - spike/contract.jl
    - spike/test/test_simulator.jl
  modified:
    - spike/Project.toml
    - spike/Manifest.toml
    - spike/test/runtests.jl

key-decisions:
  - "summary() extends Base.summary dispatched on MultiChannelImage (avoids shadowing Base's generic while keeping the plan's summary(mci) call site)"
  - "induced_mu = mean(skipmissing(summary(mci))) replicates src/bayes.jl::_prepare_data WITHOUT including bayes.jl (Turing/GLMakie absent from the lean env) — decoupling-faithful"
  - "CairoMakie co-resolved cleanly with NeuralEstimators 0.2.1 (Makie 0.24, NE-aligned weakdep) — no fallback to isolated plotting process needed"

patterns-established:
  - "include() read-only coupling: src/ reached via include(), never edited; root stays byte-identical to f581d95"
  - "StatsBase + Statistics + Images loaded BEFORE the bare-file includes (correlation() builds its method Dict at call time)"
  - "Contract traps encoded as test assertions: background-not-zero (positive image -> ~0 missing) + negative control (hard-zero background -> materially more missing)"

requirements-completed: [SIM-03]

# Metrics
duration: 15min
completed: 2026-06-26
---

# Phase 2 Plan 01: Forward-Simulator Summary Contract Bootstrap Summary

**Spike env extended with the imaging/plotting stack (no NeuralEstimators 0.2.1 downgrade, no CUDA), and the fixed 8x8 SIM-03 summary contract proven end-to-end on a synthetic image through the UNCHANGED src/ patch()/correlation() via a read-only include() boundary.**

## Performance

- **Duration:** ~15 min
- **Started:** 2026-06-26T18:46:40Z
- **Completed:** 2026-06-26T19:01:48Z
- **Tasks:** 2
- **Files modified:** 5 (2 created, 3 modified)

## Accomplishments
- Extended the isolated spike env with StatsBase, Images, ImageFiltering, CairoMakie, HypothesisTests — resolver kept NeuralEstimators at exactly v0.2.1 and Flux at 0.16.10, pulled NO CUDA install node; spike/Manifest.toml re-frozen.
- Established `spike/contract.jl`: the read-only `include()` coupling to `src/LoadImages.jl` + `src/colocalization.jl` with `build_mci`, `summary`, `induced_mu` helpers calling the source functions UNCHANGED.
- Proved the SIM-03 summary contract on a synthetic strictly-positive correlated 2-channel image: a valid `MultiChannelImage` (pixel_size = pixel-dim tuple, D-11) yields a fixed 8x8 per-patch Pearson rho (D-10) with all-finite entries and a finite induced mu; a hard-zero-background negative control proves the background-not-zero trap matters.
- Hardened `runtests.jl` with an automated resolve-risk gate (NeuralEstimators==v0.2.1 by UUID) and wired the new SIM-03 testset into the same harness. Full `julia --project=spike spike/test/runtests.jl` is green (19 tests: Phase-1 smoke + CUDA-absence + pin + SIM-03), exit 0.
- Decoupling invariant holds: root Project.toml / Manifest.toml / src/ byte-identical to baseline f581d95.

## Task Commits

Each task committed atomically (Task 2 followed RED → GREEN TDD):

1. **Task 1: Extend spike env + resolve-risk/CUDA/pin guards** - `caf33e1` (feat)
2. **Task 2 (RED): Failing SIM-03 summary-contract testset** - `f68cc4d` (test)
3. **Task 2 (GREEN): include() contract boundary, SIM-03 proven** - `329d93e` (feat)

No REFACTOR commit — the GREEN implementation was already clean.

## Files Created/Modified
- `spike/contract.jl` - Read-only include() coupling boundary to src/; `build_mci` (D-11 pixel-dim tuple), `summary` (8x8 patch/correlation, extends Base.summary on MultiChannelImage), `induced_mu` (mean of skipmissing rho).
- `spike/test/test_simulator.jl` - Wave-0 SIM-03 contract testset: synthetic positive image, D-11 construction, 8x8 finite-rho + background-not-zero trap, hard-zero negative control.
- `spike/Project.toml` - Added StatsBase, Images, ImageFiltering, CairoMakie, HypothesisTests.
- `spike/Manifest.toml` - Re-frozen reproducibility artifact (NeuralEstimators 0.2.1, Flux 0.16.10, no top-level CUDA, julia 1.12.6).
- `spike/test/runtests.jl` - Added NeuralEstimators==v0.2.1 resolve-risk gate to the D-04 testset; included test_simulator.jl.

## Decisions Made
- `summary(mci)` extends `Base.summary` dispatched on `MultiChannelImage` — keeps the plan's `summary(mci)` call site without shadowing or clobbering Base's generic (a plain top-level `summary(...) = ...` would error: "function Base.summary must be explicitly imported to be extended").
- `induced_mu = mean(skipmissing(summary(mci)))` replicates `src/bayes.jl::_prepare_data`'s reshape/!ismissing-filter + mean WITHOUT including bayes.jl (Turing/GLMakie not in the lean env). Decoupling-faithful per the PLAN <interfaces> note.
- Used the resolver-default current CairoMakie (Makie 0.24 series, NeuralEstimators 0.2.1's declared weakdep compat). It co-resolved and precompiled cleanly alongside the pinned stack — the RESEARCH "isolate plotting in a separate process" fallback was NOT needed.

## Deviations from Plan

None - plan executed exactly as written. All acceptance criteria met without auto-fixes:
- env extended, NeuralEstimators stayed v0.2.1, no CUDA, Manifest re-frozen;
- contract.jl includes both src/ files read-only; build_mci/summary/induced_mu work;
- SIM-03 proven (<=64 finite per-patch rho at 8x8, background-not-zero + >=15-px traps respected);
- runtests.jl green (exit 0); root byte-identical to f581d95; no reimplementation of patch/correlation in spike/.

## Issues Encountered
- The `summary` name collides with `Base.summary`; resolved by extending `Base.summary` with a `MultiChannelImage`-typed method (clean dispatch, no ambiguity).

## User Setup Required
None - no external service configuration required. All work is local Julia env + spike/ files.

## Next Phase Readiness
- The SIM-03 summary contract is frozen and proven: Wave-2's forward simulator (prior.jl / forward.jl) can target `build_mci -> summary -> induced_mu` as a known-good interface rather than a moving target.
- The imaging stack (Images/ImageFiltering for PSF/shift) and plotting (CairoMakie for SBC/diagnostic figures) are installed and co-resolved with the pinned SBI stack — Wave-2/3 need no further env changes for these.
- No blockers. Decoupling guarantee (root == f581d95) intact for the rest of the phase.

## Self-Check: PASSED
- FOUND: spike/contract.jl
- FOUND: spike/test/test_simulator.jl
- FOUND commit: caf33e1
- FOUND commit: f68cc4d
- FOUND commit: 329d93e

---
*Phase: 02-forward-simulator-summary-contract*
*Completed: 2026-06-26*
