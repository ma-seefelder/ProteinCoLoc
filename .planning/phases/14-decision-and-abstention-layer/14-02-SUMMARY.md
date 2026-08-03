---
phase: 14
plan: 02
subsystem: threshold-provenance
status: COMPLETE
tags: [d-07, provenance, tau, sha-pinning, wave-1]
requires:
  - "14-01 — spike/p14/consts.jl frozen at a6c8258 (P14_DEV_SEED, P14_SALT, the five ruled bars)"
  - "spike/p13/three_way_net.jld2, tau_probe_report.jld2, three_way_gate_report.jld2 (read-only)"
provides:
  - "spike/p14/provenance.jl — p14_load_tau, the only sanctioned route to the inherited three-way cut"
  - "p14_provenance_record — the block every Phase-14 .jld2 embeds; caller extras merge UNDER required keys"
  - "p14_assert_lane_clean — the D-01 runner step 0 (T-14-11)"
  - "p14_consts_sha / p14_blob_sha — the Phase-14 pre-registration hashes"
  - "P14_P13_DIR, P14_REPO_ROOT, P14_CONSTS_PATH, P14_AMENDMENT_NOTICE, P14_PROVENANCE_LOADED"
affects:
  - "every later Phase-14 plan: the threshold may be obtained ONLY by calling p14_load_tau()"
  - "every Phase-14 runner: step 0 is p14_assert_lane_clean()"
tech-stack:
  added: []          # no dependency added; spike/Project.toml and Manifest.toml byte-unchanged
  patterns:
    - "the net artifact is the SOURCE OF TRUTH for a threshold; the constant is a cross-check"
    - "one assertion per equality, so a four-way disagreement names WHICH site moved"
    - "sha pinned on BOTH sides to named, dated literals, never merely to each other"
    - "assert the DIVERGENCE of the probe's third sha; the widened disjunction is recorded as rejected"
    - "caller extras merged UNDER required keys (the p13_calibration_meta idiom)"
    - "comment-stripped source grep, so the file's own prose cannot satisfy or defeat its guard"
key-files:
  created:
    - spike/p14/provenance.jl            # 349 lines
    - spike/test/test_p14_provenance.jl  # 269 lines, 59 assertions
  modified: []
decisions: [D-07, D-01]
metrics:
  duration: ~25 min
  completed: 2026-08-03
  tasks_completed: 2
  tasks_total: 2
---

# Phase 14 Plan 02: Load τ With Its Provenance Asserted — Summary

**τ can now be obtained in Phase 14 only by loading it.** `p14_load_tau()` reads the Phase-13 net
artifact, asserts that all four on-disk sites agree, pins seven sha values on both sides, and
asserts that the tau probe's *third* `consts_sha256` **diverges** from both pinned literals — the
trap the research flagged. No Phase-14 source contains the threshold as a literal, and a running
test enforces that lane-wide.

## Commits, in order

| # | Task | Commit | Files |
|---|---|---|---|
| 1 | The loader | `473aac3` | `spike/p14/provenance.jl` (349 lines) |
| 2 | The D-07 table as a running test | `3d94438` | `spike/test/test_p14_provenance.jl` (269 lines) |

Both are descendants of the Wave-0 freeze `a6c8258` — verified with
`git merge-base --is-ancestor a6c8258 HEAD` → **YES**. Neither commit contains a measured number:
this plan produces no result, only the audited route to an inherited one.

## Verification — actually run, with observed output

| Check | Command | Observed |
|---|---|---|
| The loader returns the inherited threshold and its agreement flag | `julia --project=spike -e 'include("spike/p14/provenance.jl"); p = p14_load_tau(); …'` | printed `tau loaded: 0.15`, **exit 0** |
| The D-07 table as a test | `julia --project=spike spike/test/test_p14_provenance.jl` | **59 pass, 0 fail, 0 error**, 2.3 s test time, **exit 0** |
| No hardcoded literal in the loader's executable code | `grep -v '^\s*#' spike/p14/provenance.jl \| grep -c '0\.15'` | `0` |
| The divergence assertion exists, exactly once, in the right direction | `grep -n 'values(P13_CONSTS_SHA)' spike/p14/provenance.jl` | one line, operator `∉` |
| A missing artifact throws rather than defaulting | `… try p14_load_tau(net_path="nope.jld2"); exit(1) catch; exit(0) end` | **exit 0** (it threw) |
| Phase-13 artifacts untouched | `git status --porcelain -- spike/p13` | empty |
| `src`, `spike/Project.toml`, `spike/Manifest.toml`, `corpus`, `spike/p14/consts.jl` byte-unchanged | `git diff --exit-code HEAD -- …` | **exit 0** |
| Wave-0 suite still green (a new file landed in `spike/p14/`) | `julia --project=spike spike/test/test_p14_consts.jl` | **273 pass, 0 fail**, exit 0 |
| CPU-only testset present | `grep -c 'CUDA' spike/test/test_p14_provenance.jl` | `1` |
| The trap literal is asserted by value | `grep -c '640ec90e…aca4cf' spike/test/test_p14_provenance.jl` | `1` |

### Falsification checks — run, observed, and reverted

A test that has never been seen to fail is not evidence. Two were run and both bit:

1. **Flipping the divergence assertion's operator.** `probe_consts_sha ∉ values(P13_CONSTS_SHA)`
   → `∈`, patched in place. `test_p14_provenance.jl` **exit 1**. Reverted with `git checkout --`;
   `git diff --exit-code` confirmed clean. (The first attempt at this patch used `perl -i -pe`
   with `\x{2209}`, which silently substituted nothing on the UTF-8 file — `grep -c` returned `0`
   and the run never happened. Redone through Julia's own string replace with an
   `@assert s2 != s`, which is why the patch is provably real rather than assumed.)
2. **Planting a hardcoded threshold in the lane.** A scratch
   `spike/p14/_falsify_scratch.jl` containing `const P14_SCRATCH_TAU = 0.15` made the
   "tau is loaded, not written here" testset report **1 failure at line 98**, whole run **exit 1**.
   Scratch file deleted; `git status --porcelain -- spike/p14` empty afterwards.

The two checks are independent: (1) exercises the loader's own integrity assertion, (2) exercises
the lane-wide source scan, which is the half that catches a *future* file rather than this one.

## What was built

### `spike/p14/provenance.jl`

**`p14_load_tau(; net_path, probe_path, gate_path)`** — the only sanctioned route to the
inherited threshold. In order: every artifact must exist (`error`, never a default, never a
literal); the **net is the source of truth** because the cut travels with the net trained under
it; the probe and gate reports are hard preconditions (unlike `_recorded_ood_threshold`'s
deliberately inert `nothing` path, and the docstring says why the two differ); then

- **four-way agreement**, one assertion per equality — net vs probe, net vs gate, net vs
  `P13_TAU`, net vs `p13_tau()` — so a failure names *which* site moved rather than reporting
  that "something" drifted;
- **the sha table**, seven assertions: `h.consts_sha == P13_CONSTS_SHA.pre_amendment`,
  `p13_consts_sha() == .post_amendment`, `gate["net_consts_sha"] == h.consts_sha`,
  `probe["repo_sha"] == P13_TAU_PROBE_SHA`, `probe["simulator_sha"] == P13_TAU_SIMULATOR_SHA`,
  `probe["tau_bar"] == P13_TAU_AUC`, `P13_TAU_MEASURED_AUC >= P13_TAU_AUC`;
- **the divergence, asserted rather than tolerated** — `probe["consts_sha256"] ∉
  values(P13_CONSTS_SHA)`, with the reason (the probe ran *before* the Tier-2 append, which is
  the ordering evidence `spike/p13/consts.jl:842-846` cites) and the explicit statement that
  widening it to a disjunction destroys the integrity argument and is forbidden. The rejected
  form is named at the assertion, mirroring `spike/p13/net.jl:549-554`.

`probe_repo_dirty_at_run` is returned **honestly**. It is `true` on disk — the Phase-13 probe ran
on a dirty tree — and the docstring states that Phase 14 surfaces it and does not assert it false.

**`p14_provenance_record(prov; extra)`** — merges caller extras **under** a required block
(`merge(extra, required)`, the `p13_calibration_meta` idiom), so a runner that supplied its own
`tau` or its own `amendment` cannot overwrite the inherited ones. The required block adds
`p14_consts_sha`, `p14_consts_git_blob_sha`, `master_seed`, `salt`, `generated`, and the
SC1/SC2 `amendment` citation **in code**, so it travels with every artifact rather than living
only in a document.

**`p14_assert_lane_clean()`** — the D-01 step-0 guard (T-14-11), promoted from
`run_three_way_gate.jl:270-280`: `src` and the frozen spike environment asserted byte-unchanged
*while the run happens*, degrading to `@info` only when `git` is genuinely unavailable.

### `spike/test/test_p14_provenance.jl`

Nine testsets, 59 assertions, 2.3 s. Testset 2 **loads all four τ sites independently** rather
than reading them back off the loader's return value, so it is a check *on* the loader and not a
restatement of it. Testset 4 asserts the trap in both directions: the value is not in
`values(P13_CONSTS_SHA)`, it equals the recorded third sha, and the **loader's comment-stripped
source** is grepped to confirm the operator is `∉`, that no ` in values(...)` or `∈ values(...)`
form appears, that there is exactly one such line, and that its condition carries no `||`.

## Deviations from Plan

### Structural additions the plan implies but does not name

None is a choice about a value; each exists only so an assertion the plan *does* demand can be
expressed, or so a surface the plan *does* require is not shipped untested.

| Name | Why it exists |
|---|---|
| `_p14_load_report(path, what)` / `_p14_report_key(d, key, what, path)` | the plan asks for a house-`catch` load that errors rather than returning `nothing`, and for messages naming the failing artifact; expressed once instead of six times |
| `p14_consts_sha()` | the required record key `p14_consts_sha` needs a named accessor; mirrors `p13_consts_sha()` (`net.jl:524`) so the two are recognisably the same device |
| `P14_CONSTS_PATH`, `P14_AMENDMENT_NOTICE` | the record's `p14_consts_git_blob_sha` and `amendment` keys, hoisted to named constants so the amendment text is greppable rather than buried in a function body |
| `P14_PROVENANCE_LOADED` | the include-guard sentinel, declared unconditionally at the foot of the const block. **Not** a seed and not a path constant — CONVENTIONS C-01: a sentinel a later phase might legitimately mirror is how a whole guard body silently skips |
| testset "the provenance record cannot be overwritten by a caller" | the plan's Task-2 action lists six testsets, none covering `p14_provenance_record`. Its merge-under semantics are the *only* thing preventing a runner from overwriting the inherited threshold in its own artifact — shipping that untested would leave the plan's own `must_haves` half-enforced |
| testset "a missing artifact is a block, never a default" | the plan states this as a Task-1 *acceptance criterion* (a shell one-liner) rather than a testset; promoted into the suite so it keeps running after this plan closes |

### Not done, deliberately

- **`spike/test/runtests.jl` not touched.** It exits 1 early at the Phase-4 `SPEEDUP_GATE`, so
  wiring Phase-14 tests into it would make them unrunnable. Per-file runs are the signal
  (14-PATTERNS §4.8). The file mentions `runtests.jl` once, in a comment saying it is *not* a
  gate.
- **`STATE.md` / `ROADMAP.md` not written** — the orchestrator owns those.
- **`14-VALIDATION.md` not touched** — its D-07 row describes what this plan implements; it was
  not amended, and nothing here required amending it.

### Auto-fixed issues

None. No bug, no missing critical functionality, no blocker was encountered. The one wrong turn —
the `perl -i -pe` substitution that silently did nothing — was in a *falsification check*, caught
by its own `grep -c` returning `0`, and is recorded above rather than quietly redone.

## Assumption Drift (advisory)

**1. Execution environment — planned and actual now agree**

- **Found during:** startup
- **Planned:** main working tree, no worktrees (`execution_environment: main-working-tree-no-worktrees`)
- **Actual:** main working tree, no worktrees. `git worktree` was never invoked.
- **Why it matters:** 14-01's original BLOCKED record flagged a worktree dispatch and warned it
  would break later Phase-14 plans needing the gitignored `spike/data/cache/p13/` pool. That
  advisory still stands for the plans that *read the pool*; this plan reads only tracked sources
  and three committed `.jld2` artifacts, all present.

**2. The plan's file list vs `14-PATTERNS.md` §1**

- **Planned:** the plan itself flags that `provenance.jl` is an addition to the seven module files
  enumerated in `14-PATTERNS.md` §1.
- **Actual:** confirmed — `spike/p14/` now holds `consts.jl` and `provenance.jl`, and the test's
  lane scan globs the directory rather than enumerating it, so the six later modules
  (`result.jl`, `fuse.jl`, `fdr.jl`, `conformal.jl`, `posterior.jl`, `decide.jl`) are picked up
  automatically when they land. Six `@info` lines announce this on every run.
- **Why it matters:** a hand-maintained file list is how a new file escapes a lane-wide guard.
  Recorded so the choice is visible rather than inferred.

Nothing else drifted materially.

## Constraint Verification (all hold, all checked on disk)

| Constraint | Check | Status |
|---|---|---|
| `src/` untouched (D-01) | `git diff --exit-code HEAD -- src` | exit 0 |
| `spike/Project.toml` + `Manifest.toml` byte-unchanged (`test_p12_decoupling.jl:176-178`) | `git diff --exit-code HEAD -- …` | exit 0 |
| No dependency added or changed | `JLD2` (frozen sixteen), `SHA`, `Dates`, `Test` (stdlib) | ok |
| Phase-16 sealed holdout not consumed (`test_p12_decoupling.jl:210`) | `corpus/` never read; `git diff --exit-code HEAD -- corpus` | exit 0 |
| `spike/p14/consts.jl` not modified (frozen at `a6c8258`) | `git diff --exit-code HEAD -- spike/p14/consts.jl` | exit 0 |
| Phase-13 lane read-only | `git status --porcelain -- spike/p13` | empty |
| Every result-producing commit descends from the freeze | `git merge-base --is-ancestor a6c8258 HEAD` | YES |
| Corrections APPENDED, never overwritten | both files are first writes; this summary is a first write | ok |
| `runtests.jl` not wired as a gate | mentioned once, in a comment saying it is *not* a gate | ok |
| `STATE.md` / `ROADMAP.md` not modified by this executor | not staged, not written | ok |
| Files modified outside the allowed three | none | ok |

## Known Stubs

None. Neither file contains a placeholder, a TODO, or a value awaiting a later run. The six
`@info "… does not exist yet"` lines are **not** stubs: they are the arms-itself-later idiom the
plan's Task-2 action asks for, and the guard they announce is already live for every file that
does exist.

## Threat Flags

None — no network surface, no file write at runtime, no new trust boundary. `p14_load_tau` and
its helpers only `read`. The plan's four register entries are discharged:

- **T-14-14** (tampering with τ) — the loader errors on a missing artifact instead of defaulting,
  the four-way agreement is asserted one equality at a time, and the lane-wide grep test
  (falsified above) forbids a literal in Phase-14 executable code.
- **T-14-15** (tampering with the divergence assertion) — asserted in BOTH the loader and the
  test; the loader's source is grepped for a widened form; the rejection is recorded in the error
  message so a later reader knows it was considered rather than missed. Falsified above.
- **T-14-10 / T-14-11** (`src/`, `spike/p13/` writes; clean-start-dirty-finish) —
  `p14_assert_lane_clean()` exists and is exercised by the test, plus the git testset.
- **T-14-SC** (package installs) — nothing fetched, nothing installed.

## Self-Check

- `spike/p14/provenance.jl` — **FOUND** (349 lines)
- `spike/test/test_p14_provenance.jl` — **FOUND** (269 lines)
- commit `473aac3` — **FOUND**
- commit `3d94438` — **FOUND**
- `git status --porcelain` limited to the permitted paths — **clean**. The only other entries are
  `.planning/HANDOFF.json`, `.planning/config.json` and three untracked planning paths, none
  written or staged by this executor.

## Self-Check: PASSED

Plan 14-02's success criteria are met and observed, not inferred:

- **No Phase-14 code path can obtain τ except by loading it with provenance asserted** — enforced
  by the lane-wide literal scan, which was falsified and bit.
- **A future divergence in any of the four τ sites, or in any pinned sha, fails a test in under
  10 seconds** — the suite runs in **2.3 s** wall-clock test time.
