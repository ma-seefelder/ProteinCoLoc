---
phase: 09-cross-method-comparator-harness
plan: 03
subsystem: testing
tags: [comparator, shared-inputs, seeded, random123, philox, multichannelimage, regime-labels, reproducibility]

# Dependency graph
requires:
  - phase: 09-01
    provides: "spike/comparator/config.jl (MASTER_SEED, pre-declared thresholds)"
  - phase: 02-simulator
    provides: "sample_prior, simulate_pair, build_mci frozen chain (read-only via contract.jl)"
  - phase: 03-data
    provides: "sample_rng Philox keyed-stream seeding primitive"
provides:
  - "build_shared_inputs: seeded θ-grid → Vector{MultiChannelImage} with ground-truth regime labels"
  - "regime_of + pre-declared regime band consts (:coloc/:random/:exclusion)"
  - "D-03 external-vector passthrough (labels :unknown, no Phase-8 dependency)"
affects: [09-05 table assembly, 09-06 reproducibility gate, 16 blind external eval]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Reject-and-redraw regime routing keyed purely by (master_seed, global_index) — order/thread independent"
    - "Guarded-include composition mirroring generate.jl (contract.jl FIRST) to reach the frozen chain read-only"

key-files:
  created:
    - spike/comparator/inputs.jl
  modified: []

key-decisions:
  - "Regime bands fixed as consts (REGIME_COLOC_BAND=0.30, REGIME_EXCLUSION_BAND=-0.30) — pre-declared, anti-data-snooping (D-14)"
  - "Reject-and-redraw over increasing global indices to span all three regimes while staying byte-reproducible"
  - "Assembly reconstructs a fresh keyed rng per index so sample_prior→simulate_pair matches generate_sample byte-for-byte"
  - "External passthrough labels :unknown / rho_true missing — never fabricate ground truth (T-09-07)"

patterns-established:
  - "Pattern: seeded regime-labelled shared-input source with input-source-agnostic passthrough"

requirements-completed: [CMP-03]

# Metrics
duration: 12min
completed: 2026-07-02
---

# Phase 9 Plan 03: Seeded Regime-Labelled Shared-Input Builder Summary

**`build_shared_inputs` emits a reproducible `Vector{MultiChannelImage}` spanning :coloc/:random/:exclusion via a reject-and-redraw θ grid on the frozen simulator chain, plus a D-03 external-vector passthrough that labels untrusted inputs `:unknown`.**

## Performance

- **Duration:** ~12 min
- **Started:** 2026-07-02T09:42:00Z
- **Completed:** 2026-07-02T09:54:26Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Seeded θ-grid builder that spans all three colocalization regimes (verified: 3/3/3 coverage at n_per_regime=3) by reject-and-redraw routing keyed purely by `(master_seed, global_index)`.
- Two builds under the same master seed produce byte-identical images AND labels (D-12 / T-09-08), reusing the FROZEN `sample_prior → simulate_pair → build_mci` chain read-only via contract.jl.
- Pre-declared fixed regime band consts (D-14 anti-data-snooping); `regime_of` labels consistent with ρ_true bands verified across a full grid.
- D-03 input-source-agnostic passthrough: any `Vector{MultiChannelImage}` returned unchanged with `:unknown` labels and `missing` ρ_true — no Phase-8 dependency taken.

## Task Commits

Each task was committed atomically:

1. **Task 1: Seeded regime-labelled shared-input builder with external-vector passthrough** - `2f4f417` (feat)

**Plan metadata:** committed with this SUMMARY.

_Note: single-file plan (files_modified restricted to `spike/comparator/inputs.jl`); TDD gate handled via the plan's inline julia smoke (see TDD Gate Compliance)._

## Files Created/Modified
- `spike/comparator/inputs.jl` - `build_shared_inputs` (keyword + passthrough methods), `regime_of`, pre-declared regime band consts; guarded-include of the frozen chain (contract.jl → prior.jl → forward.jl → seeding.jl) plus sibling config.jl for MASTER_SEED.

## Decisions Made
- **Regime bands** fixed at ±0.30 on ρ_true (which spans [-0.99, 0.99] via GHAT_RHO_KNOTS): coloc ≥ +0.30, exclusion ≤ -0.30, random between. Declared as `const`s before any table exists (D-14).
- **Reject-and-redraw** by scanning increasing global indices and bucketing each draw into its target regime keeps the builder fully seeded and order/thread-independent (acceptance depends only on the draw's own `(master_seed, idx)`), so it stays byte-reproducible.
- **Assembly reconstructs a fresh keyed rng per accepted index** so the `sample_prior → simulate_pair` continuation is byte-identical to `generate_sample`; the routing draw uses a separate discarded rng object.
- **Passthrough master_seed = missing** to keep the NamedTuple shape stable while signalling "no seed / no ground truth" for external corpora.

## Deviations from Plan

None - plan executed exactly as written. Band boundaries and `n_per_regime` default (8) were within Claude's discretion per the plan (`n_per_regime::Int=…`, bands "declare as consts here").

## TDD Gate Compliance

The task carried `tdd="true"`, but the plan restricts `files_modified` to `spike/comparator/inputs.jl` only and the isolation constraint forbids touching `spike/test/`. Consequently there is no separate `test(...)` commit — the RED/GREEN cycle was exercised via the plan's `<verify><automated>` inline julia smoke:
- **RED:** smoke failed (file absent) before implementation.
- **GREEN:** smoke printed `OK n=6` after implementation; extended checks confirmed full regime coverage, label↔band consistency, and passthrough `:unknown`.

The formal reproducibility/determinism test is wired into the phase gate in plan **09-06** (per this plan's `<verification>` note); no test file is created here by design.

## Issues Encountered
None. (Git emitted a benign LF→CRLF warning on Windows; no action needed.)

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- The comparator now has a reproducible, regime-labelled shared-input source that is also input-source-agnostic for later external corpora (Phase 8 / Phase 16).
- Ready for 09-05 table assembly (consumes `images` + `regime` + `rho_true`) and the 09-06 reproducibility gate.

## Self-Check: PASSED

- FOUND: spike/comparator/inputs.jl
- FOUND: commit 2f4f417

---
*Phase: 09-cross-method-comparator-harness*
*Completed: 2026-07-02*
