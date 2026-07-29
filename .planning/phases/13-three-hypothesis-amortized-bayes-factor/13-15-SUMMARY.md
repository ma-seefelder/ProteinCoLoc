---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 15
subsystem: real-image-arm
tags: [d-15-amended, d-16, real-images, ingestion, alpha-ladder, read-only, anti-snooping]
requires:
  - spike/p13/consts.jl (13-01, Tier-1 pre-registration)
  - spike/p13/alpha_series.jl (13-04, the shared alpha transform)
  - spike/contract.jl (frozen summary contract, read-only src/ include chain)
  - spike/data/encode.jl (the 128-row D-01 encoding)
provides:
  - P13_REPO_ROOT
  - real_tif, load_real, real_pair
  - real_mbar_of, real_mbar
  - real_grid_provenance
  - real_provenance
  - real_alpha_ladder
  - verify_real_readonly_digest
affects:
  - 13-16 (the Phase-11-dependent three-way read runs on this ingestion)
tech-stack:
  added: []
  patterns:
    - "@__DIR__-anchored repo root (spike/simulator/calibration.jl:59), never cd()"
    - "include-time self-check on read(@__FILE__, String) (11-10-PLAN.md:191-195)"
    - "comment-stripped source-grep gate (test_p13_alpha.jl:123-124)"
    - "forbidden literals assembled from fragments so a gate cannot trip on its own text"
key-files:
  created:
    - spike/p13/real_images.jl
    - spike/test/test_p13_real.jl
  modified: []
decisions:
  - "The frozen c1/c2 ghat anchor regression is a READ-CHAIN IDENTITY check, not a colocalization claim -- ghat.jl now records the c1/c2 figures as superseded by c2/c3 for colocalization purposes, so the caveat is written into the file header rather than the pre-registration being reinterpreted."
  - "real_mbar takes G as a keyword and REJECTS anything but the frozen 8x8 grid, because patch_summary (spike/contract.jl:85) is grid-8 and takes no grid argument; Phase 13 introduces no second grid."
  - "The not-a-gate source assertion matches threshold words per underscore SEGMENT, not as a substring, because a substring match reports P13_REAL_NAMING_CORRECTION as a bar (naMINg) and would force the honesty string onto the exemption list."
metrics:
  duration: ~35 min
  completed: 2026-07-27
  tasks: 2
  commits: 2
  files: 2
---

# Phase 13 Plan 15: Real-Image Arm (Ingestion + Ladder) Summary

> ## ⚠ SUPERSEDING NOTE — PLAN 13-17 AMENDED THE CHANNEL PAIR; THIS DELIVERABLE NEEDED NO RE-AUTHORING
>
> **The measured tables below are RETAINED IN FULL and nothing in them is rewritten.** The measured
> `0.329163 / 0.248052` are a **READ-CHAIN IDENTITY CHECK ON A SUPERSEDED PAIR**: they prove this
> ingestion is bit-for-bit the same load path the frozen `ghat` calibration used. **That claim is
> unaffected by the amendment** — the identity holds on whichever pair it is exercised on. What is
> superseded is only their status as *the colocalization reference*, because `c1` is the DAPI/Hoechst
> nuclear counterstain rather than a target protein.
>
> **This deliverable was pair-agnostic by construction and needed no re-authoring.** Every entry point
> takes `channels` as a keyword, so `13-D15-AMENDMENT.md` (2026-07-29) was a re-point rather than a
> rewrite: plan 13-17 changed **not one executable line** in `spike/p13/real_images.jl`, and the five
> `channels = P13_REAL_CHANNEL_PAIR` keyword defaults are byte-identical. Only the header prose blocks
> (b2), (e) and the `real_mbar` docstring were updated, to state what is **operative** rather than what
> is **frozen**.
>
> **The amended colocalization anchors were re-measured through this same code path by 13-17** and
> reproduce `P13_REAL_ANCHOR_MBAR = (positive = 0.4603, negative = 0.3815)` on the operative `(2, 3)`
> pair within `P13_REAL_ANCHOR_TOL = 1e-3`. `spike/test/test_p13_real.jl` now exercises **both**: the
> c1/c2 identity check with the pair **and** the expected values pinned to literals on both sides — so
> it cannot follow the constants and quietly change meaning — and the operative c2/c3 regression.
>
> See `13-D15-AMENDMENT.md` §5, `13-REPORT.md` §8c, and `13-17-SUMMARY.md`.

The six committed microscopy TIFFs now load into the amortized summary path through the frozen
`MultiChannelImage(name, paths, channels)` constructor from a `spike/p13/` depth without touching
the process CWD, reproduce the frozen `ghat` lineage anchors to 5e-5, feed the fixed 8x8 grid with
the truncation derived by `patch()`'s own integer division, and drive the D-16 alpha ladder through
the one shared `alpha_segregate` code path -- with the sealed holdout unreachable by an include-time
self-check and `test/test_images/` proven byte-unchanged by a git-index digest.

## What Was Built

**`spike/p13/real_images.jl`** (474 lines) --- `P13_REPO_ROOT`, `real_tif`, `load_real`,
`real_pair`, `real_mbar_of`, `real_mbar`, `real_grid_provenance`, `real_provenance`,
`real_alpha_ladder`, `verify_real_readonly_digest`. Flat top-level functions, no module wrapper,
guarded includes in the documented order (`consts.jl` -> `contract.jl` -> `data/encode.jl` ->
`alpha_series.jl`).

**`spike/test/test_p13_real.jl`** (378 lines) --- the five D-15 testsets named by the research Test
Map, the A14 transposition falsifier, a no-RNG-stream testset and the CPU-only postscript.

## Measured Results (all observed by running the code, not inferred)

### Ingestion vs the frozen `ghat` anchors

| Condition | frozen `P13_REAL_ANCHOR_MBAR` | measured `real_mbar` | abs difference | tol `P13_REAL_ANCHOR_TOL` |
|---|---|---|---|---|
| positive | +0.3292 | +0.329163 | 3.7e-5 | 1e-3 |
| negative | +0.2481 | +0.248052 | 4.8e-5 | 1e-3 |

Frame size `(1028, 1376)`, element type `Matrix{Float64}`, values inside `[0, 1]` on every channel
(first channel measures min 6.1036e-5, max 0.19794). Zero patches `missing` on either fixture.

### Derived grid truncation (never a literal)

| quantity | derived value |
|---|---|
| grid | 8 x 8 |
| patch size | 128 x 172 |
| rows dropped | 4 (0.39%), matches `P13_REAL_GRID_TRUNCATION_ROWS` |
| columns dropped | 0 |
| pixels per patch | 22016 |
| missing patches | 0 |

No crop, resize or pad step exists anywhere in the arm.

### Real alpha ladders (both conditions, `P13_REAL_ALPHA_GRID` = `P13_ALPHA_LADDER`)

Positive: Otsu ch1 = 0.0286607, mask fraction 13.493%, background floor b = 0.0062409,
**source zero count = 2** (of 1 414 528).

| alpha | m-bar | n_missing | realized max | sum ratio | n_zero |
|---|---|---|---|---|---|
| 0.000 | +0.32916 | 0 | 0.65669 | 1.0 | 2 |
| 0.125 | +0.27312 | 0 | 0.57538 | 1.0000000000000002 | 2 |
| 0.250 | +0.21231 | 0 | 0.57779 | 1.0000000000000002 | 2 |
| 0.375 | +0.14958 | 0 | 0.61037 | 1.0000000000000002 | 2 |
| 0.500 | +0.08690 | 0 | 0.64296 | 1.0000000000000002 | 2 |
| 0.625 | +0.02465 | 0 | 0.67554 | 1.0000000000000002 | 2 |
| 0.750 | -0.03851 | 0 | 0.70812 | 1.0000000000000002 | 2 |
| 0.875 | -0.10423 | 0 | 0.74070 | 1.0000000000000002 | 2 |
| 1.000 | -0.16787 | 0 | 0.77329 | 1.0000000000000002 | 2 |

Negative: Otsu ch1 = 0.0357099, mask fraction 22.919%, background floor b = 0.0059052,
**source zero count = 3** (of 1 414 528).

| alpha | m-bar | n_missing | realized max | sum ratio | n_zero |
|---|---|---|---|---|---|
| 0.000 | +0.24805 | 0 | 0.58686 | 1.0 | 3 |
| 0.125 | +0.18357 | 0 | 0.51424 | 1.0000000000000002 | 3 |
| 0.250 | +0.11407 | 0 | 0.44162 | 1.0000000000000002 | 3 |
| 0.375 | +0.04554 | 0 | 0.41901 | 1.0 | 3 |
| 0.500 | -0.01866 | 0 | 0.44459 | 1.0000000000000002 | 3 |
| 0.625 | -0.07918 | 0 | 0.47018 | 1.0000000000000002 | 3 |
| 0.750 | -0.13912 | 0 | 0.49576 | 1.0000000000000002 | 3 |
| 0.875 | -0.20059 | 0 | 0.52134 | 1.0000000000000002 | 3 |
| 1.000 | -0.25198 | 0 | 0.54693 | 1.0000000000000002 | 3 |

Both mask fractions sit inside `P13_ALPHA_MASK_FRACTION_BOUNDS = (0.01, 0.40)`; every realized
maximum is below `P13_ALPHA_MAX_VALUE_BOUND = 1.0` and is RECORDED, never clamped. m-bar is
strictly monotone decreasing on both and crosses zero strictly inside the grid (between alpha 0.625
and 0.750 on positive, between 0.375 and 0.500 on negative). `alpha = 0` reproduces the input
BITWISE (exact `==`) on both fixtures, and `count(iszero, .)` is constant at the source count at
every rung.

**These are behaviour numbers, not correctness numbers.** The fixtures carry no colocalization
ground-truth label, so nothing above is scored against a threshold.

### Read-only proof

`sha256(git ls-files -s test/test_images)` = `eaee22f9185460fd910788026e7a33511f6de373f313459b96bb43f3fc1ef276`

`git status --porcelain test/test_images/` is empty and `git diff --quiet HEAD -- src/ test/ spike/Project.toml spike/Manifest.toml` exits 0, both observed after each commit.

### Suite result

`julia --project=spike spike/test/test_p13_real.jl` -> **exit 0, 245 pass, 0 fail, 0 error, 1 broken**
(the deliberate `@test_skip` for plan 13-16's not-yet-existing `run_p13_realimage.jl`), **25 s wall
clock**, well inside the 60 s budget. `julia --project=spike -e 'include(...)'` also exits 0, so the
file is harness-includable as well as standalone-runnable.

## Verification (each command was RUN; the observed result is what is reported)

| # | Check | Observed |
|---|---|---|
| 1 | `julia --project=spike spike/test/test_p13_real.jl` | exit 0 |
| 2 | `real_mbar("positive")` / `real_mbar("negative")` within 1e-3 of 0.3292 / 0.2481 | 0.329163 / 0.248052, PASS |
| 3 | `grep -v '^\s*#' spike/p13/real_images.jl \| grep -c 'open_sealed_holdout\|_apply_mask!\|apply_mask!\|imresize'` | 0 |
| 4 | `grep -v '^\s*#' spike/p13/real_images.jl \| grep -c '1028\|1376\|22016\|321'` | 0 (every dimension read from consts or derived) |
| 5 | `grep -c 'cd(' spike/p13/real_images.jl` | 0 |
| 6 | `grep -c 'all(.*\.> 0)' spike/test/test_p13_real.jl` | 0 (forbidden invariant never written) |
| 7 | `grep -c 'count(iszero' spike/test/test_p13_real.jl` | 3 |
| 8 | `grep -c 'startswith(strip' spike/test/test_p13_real.jl` | 1 |
| 9 | the five Test-Map testset names present | 5 |
| 10 | `git status --porcelain test/test_images/` | empty |
| 11 | `git diff --quiet HEAD -- src/ test/ spike/Project.toml spike/Manifest.toml` | exit 0 |
| 12 | post-commit deletion scan over both commits | no deletions |

**No dependency was added.** `spike/Project.toml` and `spike/Manifest.toml` are byte-unchanged; the
TIFF read path uses `Images` (already resolved through `spike/contract.jl`) and `SHA`, a Julia
standard library. **The seal was not opened.** `open_sealed_holdout` is never called, no path under
the sealed provenance directory is constructed anywhere, and the property is enforced twice: by an
include-time self-check inside `real_images.jl` on its own text, and by a testset that greps every
`.jl` file discovered under `spike/p13/` with comment lines stripped.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] `patch_summary` takes no grid argument**
- **Found during:** Task 1
- **Issue:** The plan specifies `real_mbar(...; G = 8)` computing
  `encode_d01(patch_summary(build_mci(real_pair(cond)), G))`, but the frozen summary contract
  exposes only `patch_summary(mci)` (`spike/contract.jl:85`), which hardcodes the 8x8 D-10 grid.
  There is no two-argument method anywhere in `spike/` or `src/`.
- **Fix:** `real_mbar` keeps `G` as a keyword so the grid is stated at the call site, validates
  `G == 8` with a named `ArgumentError` explaining that Phase 13 changes no `num_patches` and
  introduces no second grid, and calls the frozen one-argument method. Adding a grid-parameterised
  `patch_summary` would have been a summary-contract change, which is out of scope and forbidden by
  the CLAUDE.md fixed-8x8 constraint.
- **Files modified:** `spike/p13/real_images.jl`
- **Commit:** 8f98c0c
- **Extra:** the suite asserts `@test_throws ArgumentError real_mbar(P13_REAL_SAMPLE; G = 4)`, so the
  refusal is a checked property rather than a comment.

**2. [Rule 1 - Bug] The not-a-gate source assertion mis-flagged the honesty string**
- **Found during:** Task 2 (caught by running the suite, which failed on the first run)
- **Issue:** Matching `FLOOR|MIN|MAX|THRESHOLD|TOL` as a SUBSTRING against `P13_REAL_*` names
  reports `P13_REAL_NAMING_CORRECTION` as a threshold constant, because "naMINg" contains `MIN`.
  That would have forced the honesty string onto the exemption list -- growing the very list the
  assertion exists to keep from growing.
- **Fix:** match per underscore-delimited SEGMENT instead (`split(n, '_')`), which is the stronger
  and semantically correct form: a threshold constant names its role as a whole segment. The
  word list was also widened to `BOUND` and `CUTOFF` so a future bar cannot enter under a synonym.
  The exemption list is unchanged at exactly two entries, both named and reasoned in-place.
- **Files modified:** `spike/test/test_p13_real.jl`
- **Commit:** 9785f95

**3. [Rule 2 - Missing critical] The forbidden invariant could not be quoted in the comment**
- **Found during:** Task 2
- **Issue:** The plan asks for an in-test comment "stating in full that `all(y2 .> 0)` is
  deliberately NOT asserted", but the plan's own acceptance criterion requires
  `grep -c 'all(.*\.> 0)'` to output 0. Writing the literal would have failed the gate.
- **Fix:** the comment states the rule in words ("an all-pixels-strictly-positive check") without
  quoting the forbidden form, which is exactly the technique `spike/p13/alpha_series.jl:182-191`
  already uses ("whose name is not even quoted here, so a source-grep gate cannot be defeated by
  prose"). The substance of the explanation -- 2 source zeros in `positive_c2`, 3 in `negative_c2`,
  of 1 414 528, and that the naive form fails a CORRECT algorithm at alpha = 0 -- is written in
  full.
- **Files modified:** `spike/test/test_p13_real.jl`
- **Commit:** 9785f95

## Assumption Drift (advisory)

**The frozen c1/c2 anchor is no longer this project's colocalization reference for these fixtures.**

- **Found during:** Task 1, while reading `spike/simulator/ghat.jl` to verify the lineage claim.
- **Planned:** the plan's `<interfaces>` block quotes `spike/simulator/ghat.jl:39` as
  `real anchor (D-16): faithful LoadImages.jl load_tiff (NOT luminance): positive mu = 0.3292,
  negative mu = 0.2481`, and the pre-registration freezes `P13_REAL_CHANNEL_PAIR = (1, 2)` and
  `P13_REAL_ANCHOR_MBAR = (positive = 0.3292, negative = 0.2481)`.
- **Actual:** `ghat.jl` now carries a SECOND hand-patch entry (2026-07-25, landed as commit ea4a7d3
  from the concurrent Phase-11 work) that moves the measured anchor pair from c1/c2 to c2/c3 and
  states in the file that the c1/c2 figures "were measured against the counterstain and are
  SUPERSEDED, not merely updated". `test/runtests.jl:105` records the fixture channels as
  `["blue", "green", "red"]`, so c1 is the DAPI/Hoechst nuclear counterstain, not a target protein.
  The current unmasked c2/c3 figures are positive mu = 0.4603, negative mu = 0.3815.
- **Why it matters:** the numbers themselves reproduce exactly (measured 0.329163 / 0.248052), so
  nothing failed and nothing was worked around. What drifted is the *interpretation*: the anchor
  regression is a READ-CHAIN IDENTITY check -- it proves this ingestion is bit-for-bit the same load
  path the frozen calibration used -- and it is NOT evidence that the c1/c2 pair measures
  colocalization.
- **Handling:** `P13_REAL_CHANNEL_PAIR` is Tier-1 frozen and was NOT edited. The caveat is written
  into the `real_images.jl` header as paragraph (b2) and into the `real_mbar` docstring and the
  ingestion testset comment, and every entry point (`load_real`, `real_pair`, `real_mbar`,
  `real_alpha_ladder`, `real_provenance`) takes `channels` as a keyword so plan 13-16's reported run
  can read the c2/c3 coloc pair without any change here.
- **Advisory only.** This is recorded so the drift is visible, not gated. Note that
  `P13_REAL_REDUNDANCY_PAIR = (1, 3)` also involves the counterstain channel; plan 13-16 should
  decide, in the open, which pair it headlines.

## Known Stubs

None. Every function returns measured values from the committed fixtures; there is no placeholder,
no hardcoded empty container and no mock data path.

## Threat Flags

None. This plan adds no network endpoint, no auth path and no schema at a trust boundary. It adds
one filesystem read surface (`test/test_images/`), which is already in the plan's threat register as
T-13-50/T-13-51 and is mitigated by the digest, the porcelain-status assertion, the mutating-call
source grep and the frame-size/element-type validation at ingestion.

## Notes for Plan 13-16

- No Phase-11 artifact path is hard-coded anywhere here, and nothing depends on
  `spike/p13/preconditions.jl` (plan 13-09, not yet written).
- `spike/test/test_p13_real.jl` is NOT wired into `spike/test/runtests.jl` -- that is plan 13-08's
  job, and `runtests.jl` was not touched.
- The `@test_skip` in `P13 corpus seal untouched` becomes a real assertion the moment
  `spike/p13/run_p13_realimage.jl` exists: that runner must carry its own
  `read(@__FILE__, String)` self-check.
- The reported runner should record `verify_real_readonly_digest()` before and after the run and
  assert equality.

## Self-Check: PASSED

- `spike/p13/real_images.jl` -- FOUND
- `spike/test/test_p13_real.jl` -- FOUND
- commit 8f98c0c -- FOUND
- commit 9785f95 -- FOUND
- `git status --porcelain test/test_images/` -- empty
- `git diff --quiet HEAD -- src/ test/ spike/Project.toml spike/Manifest.toml` -- exit 0
- `spike/p13/consts.jl`, `spike/test/runtests.jl`, `.planning/STATE.md`, `.planning/ROADMAP.md` --
  not modified (absent from both commits' file lists)
