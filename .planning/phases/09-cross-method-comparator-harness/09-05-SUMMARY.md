---
phase: 09-cross-method-comparator-harness
plan: 05
subsystem: comparator
tags: [dataframes, csv, jld2, sha256, traffic-light, costes, ks-uniformity, content-addressed-cache]

# Dependency graph
requires:
  - phase: 09-01
    provides: pre-declared thresholds (DIVERGENCE_WARN/FAIL, COSTES_ALPHA, MASTER_SEED)
  - phase: 09-02
    provides: classical estimator battery (manders, pearson_whole, spearman_whole, costes_p)
  - phase: 09-03
    provides: seeded regime-labelled shared-input builder (build_shared_inputs)
provides:
  - Tidy per-input comparison DataFrame (build_table) with classical + ground-truth + optional NPE columns
  - Deterministic divergence + green/amber/red traffic-light positioning column (traffic_light/_regime_score/classical_verdict/divergence)
  - Atomic content-addressed CSV+JLD2 artifact writer (write_table) over the comparator's OWN SHA-256 src_files list
  - BayesInteractomics-style audit report (comparator_audit) with band counts + Costes-p KS-uniformity
affects: [09-06 determinism/oracle gates, phase-10 related-work delta, phase-16 blind external eval]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Comparator-own SHA-256 content hash over a SEPARATE src_files list — never mutates Phase-3 HASH_SRC_FILES"
    - "Pre-declared config bands read into traffic_light (no inline threshold literals, D-14)"
    - "Atomic .tmp -> integrity-check -> mv(force=true) artifact writes mirroring spike/data/cache.jl"

key-files:
  created:
    - spike/comparator/table.jl
    - spike/comparator/audit.jl
  modified: []

key-decisions:
  - "write_table uses a local _cmp_canonical + comparator_hash mirroring hashguard.jl, so no Phase-3 hash code is called or mutated (T-09-13)"
  - "KS-uniformity reimplemented locally (_ks_uniform) rather than imported from the read-only BayesInteractomics template"
  - "NPE columns modelled as an optional npe NamedTuple (rho_hat, ood_flag); missing when absent (D-09 optional, OQ3)"

patterns-established:
  - "Guarded includes (isdefined || include) so table.jl loads standalone or after siblings"
  - "Threshold header emitted as BOTH CSV comment lines and a JLD2 meta entry (D-14 pre-registration record)"

requirements-completed: [CMP-01, CMP-04, CMP-05, CMP-07]

# Metrics
duration: ~12min
completed: 2026-07-02
---

# Phase 9 Plan 05: Comparison Table + Divergence/Traffic-Light Summary

**A tidy per-input comparison DataFrame with a deterministic divergence + green/amber/red traffic-light positioning column, persisted atomically to a content-addressed CSV+JLD2 artifact over the comparator's own SHA-256 hash inputs, plus a BayesInteractomics-style band-count + Costes-p KS-uniformity audit.**

## Performance

- **Duration:** ~12 min
- **Completed:** 2026-07-02
- **Tasks:** 2
- **Files modified:** 2 (both created)

## Accomplishments
- `build_table` assembles one tidy `DataFrame` row per shared input with the full classical estimator battery (Costes_p, M1, M2, Pearson, Spearman), the simulator ground-truth regime + rho_true, optional NPE columns, and the divergence/light positioning columns.
- Divergence semantics pinned deterministically (`_regime_score`/`classical_verdict`/`divergence`): a `:random` regime with a spuriously significant Pearson flags red (false-positive coloc); a `:coloc` regime with a non-significant Pearson flags red (missed real coloc); a `:coloc` regime with a significant high Pearson stays green.
- `traffic_light` reads the pre-declared `DIVERGENCE_WARN`/`DIVERGENCE_FAIL` bands from config — no inline literals (D-14) — and is monotone at the exact band boundaries.
- `write_table` persists both `table.csv` (with a threshold-header comment block) and `table.jld2` (exact reload + meta) atomically under a content-addressed dir named by a comparator-OWN SHA-256 over `[config.jl, classical.jl, inputs.jl, table.jl]` + canonical config — Phase-3 `hashguard.jl` left byte-unchanged.
- `comparator_audit` renders a Markdown report with traffic-light band counts and a Costes-p KS-uniformity statistic in the `generate_diagnostics_report` style.

## Task Commits

Each task was committed atomically:

1. **Task 1: DataFrame row assembly + divergence + traffic-light** - `299e5eb` (feat)
2. **Task 2: Atomic content-addressed CSV/JLD2 writer + audit** - `4594665` (feat)

## Files Created/Modified
- `spike/comparator/table.jl` - traffic_light/_regime_score/classical_verdict/divergence/build_table (Task 1) + comparator_hash/_cmp_canonical/_threshold_header/write_table (Task 2)
- `spike/comparator/audit.jl` - comparator_audit (band counts + Costes-p KS-uniformity) + local _ks_uniform helper

## Decisions Made
- Reimplemented the content-hash machinery locally (`_cmp_canonical`, `comparator_hash`, `CMP_SRC_FILES`) instead of calling `hashguard.jl:cache_hash`/`canonical`, so the comparator never touches Phase-3 hash code and `HASH_SRC_FILES` stays byte-unchanged (T-09-13).
- Reimplemented KS-uniformity (`_ks_uniform`) locally, mirroring BayesInteractomics `_ks_test_uniform`, since the sibling package is a read-only template.
- Modelled the optional NPE columns as an `npe` NamedTuple with `rho_hat`/`ood_flag` vectors; both columns are `missing` when no NPE artifact is supplied.

## Deviations from Plan

None - plan executed exactly as written. Both tasks' fixture-oracle verify commands passed as authored (`OK rows=6`, `OK dir=78d95e3a...`).

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- The comparison table + divergence/traffic-light artifact and audit are ready for the formal determinism + traffic-light oracle gates wired in plan 09-06.
- Phase-3 hash state verified untouched (`git diff spike/data/hashguard.jl` empty).

## Self-Check: PASSED

- FOUND: spike/comparator/table.jl
- FOUND: spike/comparator/audit.jl
- FOUND commit: 299e5eb (Task 1)
- FOUND commit: 4594665 (Task 2)
- VERIFIED: spike/data/hashguard.jl byte-unchanged (git diff empty)

---
*Phase: 09-cross-method-comparator-harness*
*Completed: 2026-07-02*
