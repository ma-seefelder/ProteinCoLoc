---
phase: 09-cross-method-comparator-harness
plan: 02
subsystem: comparator
tags: [julia, colocalization, costes, manders, pearson, spearman, block-scramble, philox, nan-safe]

# Dependency graph
requires:
  - phase: 09-cross-method-comparator-harness
    plan: 01
    provides: comparator/config.jl pre-declared consts (COSTES_N_SCRAMBLE/COSTES_BLOCK_PX/COSTES_SALT) + spike deps
  - phase: 02-summary-contract
    provides: spike/contract.jl read-only include of frozen src/ (patch, correlation, MultiChannelImage, build_mci)
  - phase: 03-training-data-pipeline
    provides: encode_aug Manders/Pearson formula (reproduced, not edited) + Philox disjoint-salt seeding convention
provides:
  - manders / pearson_whole / spearman_whole / patch_correlation NaN-safe classical estimator callables
  - costes_p seeded block-scramble randomization p-value (the one new Phase-9 algorithm, D-05)
  - costes_rng / _block_permute reproducible Costes null primitives on a disjoint COSTES_SALT Philox stream
affects: [09-cross-method-comparator-harness later waves (table assembly 09-05, D-13 gates 09-06), comparator-table-artifact]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Reproduce-not-refactor: hash-guarded encode.jl left byte-identical; encode_aug formula re-expressed as a new callable, consistency enforced by test"
    - "Costes block-scramble null (grid-cell position permutation, not pixel scramble) with Davison-Hinkley +1/+1 correction"
    - "Disjoint COSTES_SALT Philox stream so the significance null never collides with data/CV RNG streams"

key-files:
  created:
    - spike/comparator/classical.jl
  modified: []

key-decisions:
  - "manders/pearson_whole reproduce encode.jl:102-109 byte-for-byte in a NEW file (D-04 intent by equality, not by editing HASH_SRC_FILES entry 5) — verified aug[129/130/131] == manders/pearson_whole"
  - "_block_permute permutes BLOCK POSITIONS on a cld(W,block)xcld(H,block) grid; partial edge blocks copied with a min-overlap clamp seeded from copy(ch2) so no pixel is unassigned (Pitfall 2: never a pixel scramble)"
  - "Costes statistic = whole-image Pearson to match pearson_whole exactly; block side 7px pre-declared in config.jl (~resolution element for sigma_psf=1.3)"
  - "classical.jl brings its own using Random123/Random so it loads under the plan's minimal include chain (contract.jl + config.jl + classical.jl) without seeding.jl"

patterns-established:
  - "Every estimator is total: _safe(f,v) NaN-not-throw guard mirrored from encode.jl:80; all-zero/all-missing input returns NaN, never an exception"
  - "Frozen src/ reached ONLY via contract.jl's existing include; no second include of src/colocalization.jl (would double-define frozen functions)"

requirements-completed: [CMP-01, CMP-02]

# Metrics
duration: 9min
completed: 2026-07-02
---

# Phase 9 Plan 02: Classical Estimator Battery Summary

**A spike-local `classical.jl` exposing five NaN-safe colocalization estimators — Manders M1/M2 and whole-image Pearson reproduced byte-for-byte from the hash-guarded `encode_aug` (encode.jl untouched), whole-image Spearman and frozen 8x8 patch-grid correlation reached read-only through `contract.jl`, plus the one genuinely new algorithm: a seeded, bit-reproducible Costes block-scramble randomization p-value on a disjoint `COSTES_SALT` Philox stream.**

## Performance
- **Duration:** ~9 min
- **Completed:** 2026-07-02
- **Tasks:** 2
- **Files modified:** 1 (1 created, 0 modified)

## Accomplishments
- Created `spike/comparator/classical.jl` with `manders(mci) -> (M1, M2)` using the SAME `mci.otsu_threshold[1]/[2]` source as `encode_aug`, reproducing the encode.jl:102-106 formula exactly — verified equal to the augmented-encoder moments at positions 129/130 (and 131 for Pearson) with `encode.jl` left byte-identical (`git diff` empty).
- Implemented `pearson_whole` (= `encode_aug`'s `pearson_whole`), `spearman_whole` (StatsBase `corspearman`), and `patch_correlation(mci; method=:spearman)` reducing the FROZEN `patch(...,8)` + `correlation(...; method)` 8x8 grid via `_safe(mean, skipmissing)` — reached read-only through `contract.jl`, no second `include` of `src/`.
- Implemented the NEW Costes randomization p-value `costes_p(master_seed, idx, mci)` (D-05): a `_block_permute` that permutes block POSITIONS on the grid (canonical Costes null, NOT a pixel scramble), the Davison-Hinkley `(ge+1)/(n+1)` correction (p in `[1/(n+1), 1]`, never 0), a `NaN` guard on non-finite observed correlation, and a `costes_rng` keyed by `(master_seed XOR COSTES_SALT, idx)` disjoint from `sample_rng`/`holdout_rng`/`fold_rng`.
- Confirmed block-scramble calibration end-to-end: `costes_p` returns ~0.49 on a random-regime image (Pearson approx 0) and bottoms at the `1/201` floor on a Pearson approx 1.0 colocalized image — the calibrated behavior a pixel scramble would have destroyed (p approx 0 everywhere).

## Task Commits
Each task was committed atomically:

1. **Task 1: Reproduce the NaN-safe classical estimators (Manders, Pearson, Spearman, patch)** - `5154d24` (feat)
2. **Task 2: Implement the seeded Costes block-scramble p-value (D-05, CMP-02)** - `1fff884` (feat)

## Files Created/Modified
- `spike/comparator/classical.jl` - AGPL header + decoupling note; `_safe`, `manders`, `pearson_whole`, `spearman_whole`, `patch_correlation` (Task 1) and `costes_rng`, `_block_permute`, `costes_p` (Task 2).

## Decisions Made
- **Reproduce-not-refactor (Pitfall 1):** `encode.jl` is entry 5 of `hashguard.jl:HASH_SRC_FILES`; any byte change would invalidate the completed Phase-3 >=50k cache + Phase-4 NPE. The Manders/Pearson formula is therefore re-expressed as a new callable and D-04's consistency intent is enforced by the equality assertion (`aug[129/130/131] == manders/pearson_whole`), not by editing the hash-guarded file. The 09-06 gate will formalize this.
- **Block vs pixel scramble (Pitfall 2):** `_block_permute` permutes whole-block grid positions, preserving within-block spatial autocorrelation. Partial edge blocks are handled by a min-overlap clamp over an output seeded from `copy(ch2)`, so every pixel is assigned and no bias toward significance is introduced.
- **Self-sufficient includes:** `classical.jl` declares its own `using Random123` / `using Random` so it loads under the plan's minimal three-file include chain (the verify commands do not include `seeding.jl`). `COSTES_SALT`/`COSTES_N_SCRAMBLE`/`COSTES_BLOCK_PX` are referenced from `config.jl` (included first), never redeclared.

## Deviations from Plan

None - plan executed exactly as written.

## TDD Gate Compliance
Both tasks are marked `tdd="true"`, but the plan's declared `files_modified` is `spike/comparator/classical.jl` ONLY — the test scaffold `spike/test/test_comparator.jl` is out of scope and its CMP-01/CMP-02 gates are filled in plan 09-06 (per that file's own header and the 09-RESEARCH validation map). The RED/GREEN cycle was therefore run against the plan's inline `<verify>` smokes as the executable gates: RED confirmed before each implementation (functions absent at the prior commit — a genuine git-observable RED state), GREEN confirmed by the `julia --project=spike` verify one-liners printing `OK`. No separate `test(...)` commit was produced because creating one would require editing an out-of-scope file (isolation constraint); this is the intended split — the formal per-requirement test gates land in 09-06.

## Issues Encountered
None. Git reported the usual LF->CRLF warnings on the Windows checkout (cosmetic, no content impact).

## Threat Model Compliance
- **T-09-04 (Tampering — Phase-3/4 content-hash cache):** mitigated — the formula is reproduced in a NEW file; `git diff --quiet -- spike/data/encode.jl src/` asserted clean before every commit (encode.jl byte-identical).
- **T-09-05 (DoS — estimator on degenerate input):** mitigated — `_safe` NaN-not-throw on every reduction; verified all-zero input yields `NaN` for both `manders` and `costes_p` without throwing.
- **T-09-06 (Info Disclosure — biased Costes p via pixel scramble):** mitigated — `_block_permute` is a block-position permutation; the ~0.49 p on a random image (vs the p approx 0 a pixel scramble would produce) is direct evidence the autocorrelation-preserving null is in force.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- The five estimator callables are ready to be wired as table columns by 09-05 (`run_comparator`/table assembly) and asserted by the 09-06 D-13 finiteness/Manders-equality/Costes-reproducibility gates that replace the CMP-01/CMP-02 skipped scaffolds.
- No blockers.

## Self-Check: PASSED

---
*Phase: 09-cross-method-comparator-harness*
*Completed: 2026-07-02*
