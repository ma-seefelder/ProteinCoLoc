---
phase: 12-spatial-colocalization-map
plan: 05
status: complete
subsystem: spike-test
tags: [decoupling, pre-registration, sealed-holdout, resolve-risk, git]
requires:
  - "spike/validation/p12_consts.jl (P12_REPO_ROOT, Tier-1)"
  - "spike/test/test_p12_suite.jl (wiring, wave 1)"
provides:
  - "spike/test/test_p12_decoupling.jl — executable src/-untouched, env-frozen, sealed-holdout assertions"
affects:
  - "every later Phase-12 wave: the scan grows with the phase's source surface"
tech-stack:
  added: []
  patterns:
    - "guarded git shell-out (Cmd(...; dir = P12_REPO_ROOT)) with a _HAS_GIT availability predicate"
    - "comment-stripped source scan (_strip_comment_lines) with a declared-as-data allowlist"
    - "run-time-assembled needles so the scan can include the file that states the rule"
key-files:
  created: []
  modified:
    - "spike/test/test_p12_decoupling.jl (scaffold body replaced)"
decisions:
  - "P12_REPO_ROOT is READ from p12_consts.jl §12, not re-derived; only _HAS_GIT stays local"
  - "The sealed-holdout rule is three checks (image-load scan, allowlisted token scan, corpus byte-equality), never one grep"
  - "Testset 5 skips on this checkout because artifacts/ is gitignored — recorded, not faked"
metrics:
  tasks: 1
  commits: 2
  tests_added: 24
  duration_min: 27
  completed: 2026-07-29
requirements: [SPAT-09]
---

# Phase 12 Plan 05: Decoupling, Frozen Deps and the Sealed Holdout — Summary

Replaced the wave-1 pending scaffold in `spike/test/test_p12_decoupling.jl` with seven testsets
(24 assertions) that turn CLAUDE.md's decoupling rule, the frozen 16-package spike environment and
D-11's sealed-holdout exclusion from prose into assertions that fail the suite the moment they stop
being true — and every one of them was demonstrated to actually fail when its subject was violated.

**Commits**
- `d7195cc` — `test(12-05): make decoupling, the frozen dep set and the sealed holdout executable`
- (this SUMMARY, committed separately)

## What was built

`C:\Users\Manuel\Documents\GitHub\ProteinCoLoc\spike\test\test_p12_decoupling.jl`, 322 insertions /
13 deletions, seven testsets under one `verbose = true` parent:

| # | Testset | Assertions |
|---|---------|-----------|
| 1 | `src/ is byte-unchanged (CLAUDE.md hard constraint)` | 6 |
| 2 | `the spike environment is byte-frozen` | 2 |
| 3 | `the dependency name set is exactly the frozen sixteen (resolve-risk clause k)` | 5 |
| 4 | `the Phase-16 sealed holdout is not consumed (D-11 AMENDED)` | 7 |
| 5 | `the shipped bundle is untouched` | 0 (skipped, see below) |
| 6 | `the Phase-11 pool is read-only to Phase 12` | 3 |
| 7 | `decoupling checks ran CPU-only` | 1 |

`P12_REPO_ROOT` is read through the guarded `p12_consts.jl` include, per that file's §12
("12-05 READS these bindings rather than re-deriving its own"). Only `_HAS_GIT` is local — that is
a git-availability question Tier 1 does not answer. All 13 shell-outs use
`Cmd(...; dir = P12_REPO_ROOT)`; zero bare backticks (verified by count: 13 `Cmd(` occurrences,
13 `dir = P12_REPO_ROOT` occurrences, 0 bare-backtick `run/success/readchomp/read` calls).

The sealed-holdout rule is implemented as the three checks the plan specified, not as one grep:
**(a)** every Phase-12 source that exists at run time is comment-stripped and scanned for the sole
reachable opener or for an image-loading call (`load(`, `FileIO.load`, `build_mci`) on a line naming
the corpus; **(b)** non-allowlisted sources may not mention the corpus at all, while the four
allowlisted metadata files get the narrower "never the sealed image directory" property; **(c)**
`corpus/` byte-equality plus a positive `test/test_images` assertion that arms itself when 12-19
lands.

**Why the scan can include the file that states the rule.** Three needles
(`open_` + `sealed_holdout`, `cor` + `pus/data`, `CORPUS_` + `DATA_DIR`) are assembled at run time,
so the contiguous literals never appear in this source. That is what lets testset 4(a) scan this
file rather than exempt it — and Demo C below proves the self-scan is live: an injection into this
very file was caught by this very file.

## Verification — actually run, with observed output

Standalone (the plan's `<verify>` command):

```
julia --project=spike -e 'using Test; include("spike/test/test_p12_decoupling.jl")'
P12 decoupling: src/, deps and the sealed holdout (12-05) | 24  24  1.0s
```
24 passed, 0 failed, 0 errored. Three `@info` skips printed (testsets 4c, 5, 6 — each named below).

Whole Phase-12 aggregator (`include("spike/test/test_p12_suite.jl")`): all ten `P12-SUITE-RAN`
markers printed, every file green, including `test_p12_consts.jl`'s 44-assertion wiring manifest
(so this file is still correctly wired) and `test_p12_architecture.jl` (1141 passes, 51 s).

`grep -c P12_PENDING_SCAFFOLD spike/test/test_p12_decoupling.jl` → `0`.

**The full suite was NOT run.** `spike/test/runtests.jl` aborts at `test_npe.jl`'s pre-existing
Phase-4 `SPEEDUP_GATE` measurement, which is a user-owned pre-registration decision and out of scope
here. The Phase-12 aggregator is included FIRST in `runtests.jl` (line 174), ahead of that abort, so
these testsets do execute under the full suite; per-file and per-aggregator runs were used as the
gate signal.

## Falsification — every assertion, demonstrated or reported

Each mutation was made once, observed, and reverted byte-exactly (`git diff --quiet HEAD -- <path>`
confirmed clean after each). Nothing from the battery was committed.

| Demo | Mutation | Observed |
|---|---|---|
| A | appended a blank line to `src/results.jl` | testset 1 failed 3 assertions, **naming `src/results.jl`** by pathspec, plus the whole-tree and porcelain halves. Reverted; `src/results.jl` diff-clean. |
| B | created untracked `src/_p12_falsification_scratch.jl` | `git diff --quiet HEAD -- src` **still passed** (exactly the blind spot the plan predicted); the `git status --porcelain -- src` half failed. Deleted; `src` clean. |
| C | injected `load(joinpath("corpus", "data", ...))` into non-allowlisted `test_p12_lattice.jl` **and** `open_sealed_holdout(...)` + `"corpus/data"` into this file | all three halves of testset 4 failed, each naming file and offending line: `load_hits` caught **both** files (self-scan live), `token_hits` caught the non-allowlisted one, `image_dir_hits` caught the allowlisted one. Both reverted clean. |
| D | appended a whole-line **comment** containing `corpus/data`, `open_sealed_holdout`, `manifest.csv`, `CORPUS_MASTER_SEED`, `load(joinpath("corpus"))` to `test_p12_lattice.jl` | suite stayed **green** (24/24) — `_strip_comment_lines` runs before every `occursin`, so prose about the holdout is free and code is not. Reverted clean. |
| E | no mutation — predicate discrimination evaluated directly | `Set(keys(deps)) == Set(["Flux"])` → false; `!haskey(deps, "Flux")` → false; `!haskey(deps, "Images")` → false; version pin against `v"0.1.4"` → false; the testset-6 source tie against `"cache", "p99"` → false while `"cache", "p11"` → true; the loaded-module predicate returns true for a module that *is* loaded. Every predicate discriminates. |

### Assertions I could NOT falsify, stated plainly

1. **Testset 5 (`the shipped bundle is untouched`) is inert on this checkout.** `artifacts/` is
   gitignored (`.gitignore:63`) and `git ls-files -- artifacts/amended_v2/grid_8` is empty, so the
   file takes the `@info`-and-skip branch the plan prescribed for exactly this case. The testset
   therefore reports **0 tests**, not a passing test. The bundle is Release-hosted and content-hashed
   (`git-tree-sha1 = 90e6b63a…`), so there is nothing in this tree for git to guard. I verified the
   branch condition is real (the `git ls-files` probe returns empty) but could not exercise the
   asserting branch without tracking the bundle, which is out of scope and would itself be a change
   to the shipped lane. **This is the one plan requirement whose executable half is not live here.**
2. **Testset 6's path-string inequality is a tautology by construction.**
   `normpath(.../cache/p12) != normpath(.../cache/p11)` compares two literals built in this file and
   cannot fail without editing the test. That is exactly what the plan asked for ("by comparing the
   two path strings"), so it is implemented as specified — and then given falsifiable content: see
   the deviation below.
3. **Testset 7 (CPU-only) cannot be driven red here**, because CUDA is not installed in the spike
   environment. Demo E shows the predicate detects a module that *is* loaded, so it is not vacuous;
   it simply has no way to be violated on a CPU-only box. This matches every sibling CPU-only
   testset in the suite.

## Deviations from Plan

### Auto-fixed / strengthened

**1. [Rule 2 — missing critical check] Testset 6 given a falsifiable half**
- **Found during:** Task 1 falsification review
- **Issue:** The plan's testset-6 assertion (compare the two cache-root path strings) is
  unfalsifiable without editing the test itself, which makes it documentation rather than a guard —
  and the thing it protects (Phase 11's 54 MB pool) has already been destroyed once in this
  milestone.
- **Fix:** Implemented the plan's comparison verbatim, then tied each literal to its DECLARING site:
  `spike/data/p11_generate.jl` is asserted to declare `joinpath(@__DIR__, "cache", "p11")` (live
  today — falsified in Demo E), and a symmetric check on `spike/data/p12_generate.jl` (asserting it
  declares `"cache", "p12"` and never `"p11"`) arms itself automatically when 12-09 lands, with an
  `@info` until then. Also added the nesting check `!startswith(p11_root, p12_root * sep)`.
- **Files modified:** `spike/test/test_p12_decoupling.jl`
- **Commit:** `d7195cc`

**2. [Rule 2 — missing critical check] Allowlist-integrity and non-empty-scan assertions**
- **Found during:** Task 1
- **Issue:** Two ways the sealed-holdout scan could pass while checking nothing: an allowlist entry
  that names a file not in the scanned surface (a typo or a rename would silently exempt nothing and
  hide that it exempts nothing), and an empty source enumeration.
- **Fix:** `@test !isempty(sources)` and `@test all(a -> !isfile(a) || a in sources, allowlist)`.
- **Commit:** `d7195cc`

**3. [Rule 2 — missing critical check] `CORPUS_DATA_DIR` added to the narrow allowlist needle**
- **Found during:** Task 1, reading `corpus/config.jl:112`
- **Issue:** The plan's narrow property for the allowlisted four is "`corpus/data` appears nowhere".
  But `corpus/config.jl` binds the sealed image directory to `CORPUS_DATA_DIR =
  normpath(joinpath(@__DIR__, "data"))` — so the literal path string never appears anywhere in this
  repository, and an allowlisted file could reach the images through the constant while satisfying
  the plan's check exactly.
- **Fix:** The narrow check bans both spellings, and both needles are run-time-assembled so this
  file can be scanned by them. Serves the plan's stated property ("no `sealed_holdout` image path")
  rather than only its stated spelling.
- **Commit:** `d7195cc`

### Not a deviation, recorded for the auditor

- The plan's `<action>` for testset 1 says to define `P12_REPO_ROOT` by walking up from `@__DIR__`
  **and** (in bold) to read it from `p12_consts.jl`. These read as contradictory; the bold
  parenthetical is the ruling and the walk description is a restatement of how Tier 1 derives it.
  Implemented as: guarded include, read the binding, no local re-derivation. Consistent with
  `p12_consts.jl:568`.
- `spike/npe/train_p12_npe.jl`, `spike/data/p12_generate.jl`, `spike/p12/`,
  `spike/validation/run_p12_*.jl` and `spike/test/capture_p12_golden.jl` do not exist yet. The
  enumeration is pattern-based over directories that exist, so the scanned surface grows with the
  phase instead of freezing a list that would go stale. Today it covers 13 files.
- `spike/data/cache/p12/` was **not** created. It does not exist yet by design (12-09, wave 5) and
  no assertion fabricates it.

## Assumption Drift (advisory)

None material. The one thing worth flagging is the artifacts skip (testset 5), which is recorded
above as an unfalsifiable assertion rather than drift: the plan anticipated the untracked case and
prescribed the skip, so the outcome matches the plan even though the executable half is inert.

## Read-only invariants, verified after every mutation and after the commit

`git status --porcelain` and `git diff --quiet HEAD` are both clean for: `src`,
`spike/Project.toml`, `spike/Manifest.toml`, `corpus`, `artifacts`, `spike/validation/p12_consts.jl`
(FROZEN Tier-1, unmodified), `spike/test/test_p12_lattice.jl`, `.planning/STATE.md`,
`.planning/ROADMAP.md`. `spike/data/cache/p11` is intact at 54 MB. No Tier-2 block was appended.
The sealed holdout was not fetched and `corpus/data/` remains empty (both `physical-primary` rows
still `split = sealed_holdout`, `sha256 = PENDING-FETCH`, `bytes = 0`).

## Known Stubs

None. Three `@info` skip branches exist (testset 4c pending 12-19, testset 5 untracked bundle,
testset 6 pending 12-09); each is a documented, self-arming branch rather than a stub, and each is
named in "Assertions I could NOT falsify" or in the deviations above.

## Threat Flags

None. This plan adds no network endpoint, no auth path and no schema; it adds read-only `git`
shell-outs and read-only source-text reads inside the repository, all pinned to `P12_REPO_ROOT`.
The threat register's five mitigations (T-12-03, T-12-20, T-12-09, T-12-11, T-12-12) are each
implemented; T-12-12's is inert on this checkout for the reason recorded above.

## Self-Check: PASSED

- `spike/test/test_p12_decoupling.jl` — FOUND, 24/24 green standalone and in the aggregator
- `.planning/phases/12-spatial-colocalization-map/12-05-SUMMARY.md` — FOUND
- commit `d7195cc` — FOUND in `git log`
