---
phase: 14
plan: 03
subsystem: executable-constraints
status: COMPLETE
tags: [decoupling, d-01, d-02, d-03a, d-06, d-07, source-grep, wave-1]
requires:
  - "14-01 — spike/p14/consts.jl frozen at a6c8258"
  - "spike/test/test_p12_decoupling.jl (read-only analog; NOT edited)"
provides:
  - "spike/test/test_p14_decoupling.jl — the executable hard-constraint guard for Phase 14"
  - "_p14_sources() — the GLOBBED Phase-14 source surface that grows with the phase"
  - "the SC2-c / SC2-d / withdrawn-figure / forbidden-costes-seed greps, armed before the code they guard"
affects:
  - "every later Phase-14 plan: a violation now fails on the commit that introduces it, not at review"
  - "Wave 0+1 is complete — result-producing waves may run"
tech-stack:
  added: []          # no dependency added; spike/Project.toml and Manifest.toml byte-unchanged
  patterns:
    - "glob the lane, never enumerate it — an enumerated list is how a new file escapes a guard"
    - "collect offending `path:line: code` strings, never booleans"
    - "run-time-assembled needle so a source scan can scan itself"
    - "escaped-regex needles are self-non-matching by construction"
    - "arms-itself-later: `if isfile(...) … else @info … end`, derived from the FILE not a list"
    - "every negative paired with a POSITIVE twin, so no guard can pass vacuously"
    - "positional exemption located from the source, and failing loudly when it cannot be located"
key-files:
  created:
    - spike/test/test_p14_decoupling.jl   # 638 lines, 42 assertions, 14 testsets
  modified: []
decisions: [D-01, D-02, D-03a, D-06, D-07]
validates: ["D-01/env", "SC2-c", "SC2-d", "SC1-f"]
metrics:
  duration: ~40 min
  completed: 2026-08-03
  tasks_completed: 2
  tasks_total: 2
  per_file_pass_count: 42
---

# Phase 14 Plan 03: Every Hard Constraint, Made Executable — Summary

**Phase 14's constraints are no longer prose.** `spike/test/test_p14_decoupling.jl` turns
14-CONTEXT.md's "Constraints that bind this phase", CLAUDE.md's decoupling rule and D-01 / D-02 /
D-03a / D-06 / D-07 into **42 running assertions across 14 testsets, in 1.1 s**. It runs in Wave 1,
before six of the seven modules it guards exist, so a violation fails on the commit that
**introduces** it rather than at review.

## Commits, in order

| # | Task | Commit | Files |
|---|---|---|---|
| 1 | Environment / `src/` / corpus half | `312329f` | `spike/test/test_p14_decoupling.jl` (+392) |
| 2 | The Phase-14-specific greps | `d802650` | `spike/test/test_p14_decoupling.jl` (+246) |

Both descend from the Wave-0 freeze — `git merge-base --is-ancestor a6c8258 HEAD` → **YES**.
Neither commit contains a measured number: this plan produces no result, only the guard that later
results are produced under.

## Per-file pass count (Phase 13's precedent — the aggregate suite cannot be used)

```
P14 decoupling: src/, deps, the sealed holdout and the phase's own greps (14-03) |  42  42  1.1s
  src/ is byte-unchanged (CLAUDE.md hard constraint, D-01)                       |   6
  the spike environment is byte-frozen (D-02)                                    |   2
  the dependency name set is exactly the frozen sixteen (D-02)                   |   8
  the sealed holdout is never consumed (D-03a)                                   |  13
  the shipped bundle and the Phase-13 artifacts are read-only                    |   2
  Phase-14 caches use a NEW root, never Phase 13's                               |   3
  SC2-c: no `false` OOD flag is read as in-distribution (D-06)                   |   1
  SC2-d: every classical statistic compared to tau passes through ghat           |   1
  the WITHDRAWN corpus-derived alpha floor is not an executable default (D-03a)  |   1
  costes_p never consumes the comparator MASTER_SEED (Pitfall 7)                 |   1
  shared helpers are REUSED, never re-implemented (D-02)                         |   1
  the two alphas are never derived from one another (D-07)                       |   1
  the aggregate suite is NOT a Phase-14 gate                                     |   1
  decoupling checks ran CPU-only                                                 |   1
```

**35 after Task 1 → 42 after Task 2**, strictly greater, as the acceptance criterion requires.

The six single-assertion greps are **not** thin: each one scans every comment-stripped line of every
Phase-14 source, and each carries a positive twin that arms itself when the module it concerns
lands. Six `@info` lines announce the not-yet-armed halves on every run, so what is *not* being
checked yet is visible rather than assumed.

## Verification — actually run, with observed output

| Check | Command | Observed |
|---|---|---|
| The guard itself | `julia --project=spike spike/test/test_p14_decoupling.jl` | **42 pass, 0 fail, 0 error**, 1.1 s, **exit 0** |
| Wave 0 still green | `julia --project=spike spike/test/test_p14_consts.jl` | **273 pass, 0 fail**, exit 0 |
| Wave 1 sibling still green | `julia --project=spike spike/test/test_p14_provenance.jl` | **59 pass, 0 fail**, exit 0 |
| The Phase-12 guard was not edited | `git diff --exit-code HEAD -- spike/test/test_p12_decoupling.jl` | **exit 0** |
| `src`, both env files, `corpus`, `spike/p13`, `spike/p14` byte-unchanged | `git diff --exit-code HEAD -- …` | **exit 0** |
| Required literals present | `grep -cF` per string | `p14_ood_state` 9, `not_checked` 2, `ghat(patch_correlation(` 1, `costes_p(MASTER_SEED` 1, `0.032` 3, `runtests.jl` 3 |
| The sealed needles are concatenated | `grep -c '"open_" \* "sealed_holdout"'` | `1` |
| `ConformalPrediction` is inside a `!haskey` | `grep -n '!haskey(…"ConformalPrediction")'` | one line, `:222` |
| No count-based gate on unfiltered input | `grep -v '^\s*#' … \| grep -c '== 0$'` | `0` — every gate is `== String[]` on comment-stripped input |
| Both `key_links` from the plan frontmatter | `grep -cF` | `("spike/p14", r"^.*\.jl$")` → 1; `git diff --quiet HEAD` → 8 |

### Falsification checks — every one RUN, the mutation VERIFIED ON DISK first, then reverted

A test that has never been seen to fail is not evidence, and a falsification that never executed is
worse than none. 14-02 hit exactly that (`perl -i -pe` silently substituted nothing behind an `&&`),
so **each mutation below was confirmed present with a `grep -c` / `test -f` before the run**, and
each revert was confirmed with `git diff --exit-code` or `git status --porcelain` after.

| # | Mutation | Landed? | Result | Reverted? |
|---|---|---|---|---|
| 1 | scratch `src/zzz_probe.jl` (untracked) | `?? src/zzz_probe.jl` | **exit 1**, failure at `:179` — the untracked-case line, not the diff lines | file removed, `git status -- src` empty |
| 2 | `# 0.032` appended to `spike/p14/consts.jl` **as a comment** | `grep -c` → 1 | **exit 0** — comment-stripping is real, not claimed | `git checkout --`, diff clean |
| 3 | `x = verdict.flag == false` appended to `consts.jl` **as code** | `grep -c` → 1 | **exit 1**, SC2-c reported `spike/p14/consts.jl:173: x = verdict.flag == false` | `git checkout --`, diff clean |
| 4 | scratch `spike/p14/_falsify_scratch.jl` with one violation per remaining grep | 8 lines on disk | **exit 1**, and **six separate testsets** each named their own offending line (see below) | deleted, `git status -- spike/p14` empty |
| 5a | scratch `spike/p14/fuse.jl` with `.flag == false` **inside** `p14_ood_state` | `grep -c` → 1 | **exit 0**, and the SC2-c testset grew 1 → **4** as the positive `:not_checked` / long-form assertions armed | — |
| 5b | plus `leaked = other.flag == false` **outside** the function | `grep -c` → 2 | **exit 1**, reported `spike/p14/fuse.jl:5: leaked = other.flag == false` | deleted, `git status -- spike/p14` empty |

Check 4's six reported hits, verbatim:

```
["spike/p14/_falsify_scratch.jl:1: z = patch_correlation(mci) > tau"]
["spike/p14/_falsify_scratch.jl:2: w = 0.032", "spike/p14/_falsify_scratch.jl:3: v = 1/31"]
["spike/p14/_falsify_scratch.jl:4: q = costes_p(MASTER_SEED, idx, mci)"]
["spike/p14/_falsify_scratch.jl:5: function roc_auc(x)"]
["spike/p14/_falsify_scratch.jl:7: alpha_fdr = P14_ALPHA_CONFORMAL"]
["spike/p14/_falsify_scratch.jl:8: include(\"runtests.jl\")"]
```

**5a and 5b together are the important pair.** A positional exemption that is never exercised is
indistinguishable from a whole-file hole. 5a proves the exemption *protects* the one function
allowed to inspect a raw flag; 5b proves it protects **nothing else in the same file**. And 5a
independently proves the positive half is not decorative: the assertion count moved the moment
`fuse.jl` appeared.

## Every "Constraint that binds this phase" → the assertion that runs it

14-CONTEXT.md's `<canonical_refs>` section names four. This is the mapping, and it is the plan's
headline success criterion:

| 14-CONTEXT constraint | Now enforced by |
|---|---|
| `test_p12_decoupling.jl:176-178` — `spike/Project.toml` + `Manifest.toml` byte-unchanged (why D-02 forbids `ConformalPrediction.jl`) | testset "the spike environment is byte-frozen (D-02)" (diff **and** untracked-clean) **plus** testset 3, which bans `ConformalPrediction` / `MLJ` / `ROCAnalysis` / `MultipleTesting` / `Distances` / `CUDA` **by name**, one `@test` each |
| `test_p12_decoupling.jl:210` — the Phase-16 sealed holdout is not consumed (bounds D-03) | testset "the sealed holdout is never consumed (D-03a)" — all four checks, plus the allowlist self-consistency check and the six-TIFF positive assertion |
| `src/amortized/local_map.jl:70-72, :93` — the NOT-CHECKED semantics D-06 exists to handle | testset "SC2-c" (source grep, positional exemption, positive `:not_checked` twin). **Note:** the *behavioural* half (SC2-b — `p14_ood_state` actually returns `:not_checked` on absent/`Inf`/`NaN`) is `test_p14_fuse.jl`, a wave-2 file. This plan enforces that no source can *read around* the three-valued state; it does not and cannot yet assert the state is computed correctly. |
| `corpus/manifest.csv` / `anchor_rows.jl` — "the n=30 that caps α in D-03" | **superseded by D-03a**, and the supersession is itself executable: testset "the WITHDRAWN corpus-derived alpha floor is not an executable default" bans `0.032` / `1/31` / `1 / 31` / `1/(30+1)` from Phase-14 code, with the general `1/(n+1)` expression deliberately exempt |

Plus D-01 (`src/` byte-unchanged **including the untracked case**), D-02 (exactly sixteen deps), D-07
(the two alphas never assigned from one another), Pitfall 7 (`costes_p` and the forbidden seed),
Pitfall 5 (SC2-d's `ghat` unit guard), the no-re-implementation ban, the Phase-13-read-only /
cache-root separation, the aggregate-suite ban, and CPU-only.

## Deviations from Plan

### Auto-fixed issues

**1. [Rule 3 — Blocking] `_P14_CORPUS_METADATA_ALLOWED` needed four entries, not the two the plan named**

- **Found during:** Task 1, before the first run, from a comment-stripped grep of the two existing
  Phase-14 files.
- **Issue:** The plan specifies the allowlist as exactly
  `("spike/p14/run_p14_real_images.jl", "spike/test/test_p14_decoupling.jl")`. Two **already
  committed, frozen** files mention corpus tokens in *executable* text and would have been reported
  as violations on the very first run:
  - `spike/p14/consts.jl:330` — the D-03a substitution is recorded as a **declared deviation
    string**, and that record has to name what was substituted away *from*. The file is frozen at
    `a6c8258` and may not be edited.
  - `spike/test/test_p14_consts.jl:130,140` — `@test UInt64(P14_DEV_SEED) != UInt64(CORPUS_MASTER_SEED)`,
    the seed-disjointness assertion. Banning the token would forbid the pre-registration's own test
    from naming the seed it exists to forbid.
- **Fix:** both added to the allowlist, with the reason for **each** written at the allowlist and
  again in testset 4's header. This is not an invention: `test_p12_decoupling.jl:214-218` documents
  these two exact categories as the legitimate exemptions, and its own allowlist carries four
  entries for the same reasons. The allowlist keeps the "adding an entry is a DECISION, not a fix"
  note, and both additions remain bound by the **narrower** property (check (c)): metadata yes, the
  sealed image directory / symbol / opener never — verified passing.
- **Files modified:** `spike/test/test_p14_decoupling.jl` only. Neither frozen file was touched.
- **Commit:** `312329f`

**2. [Rule 1 — Bug] The `runtests.jl` testset flagged itself**

- **Found during:** Task 2, first run — **exit 1**, one failure, the lane's only offender being
  `spike/test/test_p14_decoupling.jl:391: @testset "runtests.jl is NOT a Phase-14 gate" begin`.
- **Issue:** a `@testset` title is a **string in code**, not a comment. Comment-stripping does not
  strip strings, so naming the runner in the title put the banned literal on an executable line.
- **Fix:** the testset is renamed to "the aggregate suite is NOT a Phase-14 gate", and the literal
  is kept in comments (where the acceptance criterion's `grep` still finds it, 3 occurrences). The
  incident is recorded **in the file, at the testset**, rather than quietly renamed away — it is the
  cheapest available demonstration that the scan is live and that comment-stripping is not the same
  thing as string-stripping.
- **Commit:** `d802650`

### Structural additions the plan implies but does not name

None is a choice about a value; each exists only so an assertion the plan *does* demand can be
expressed, or so a guard the plan *does* require cannot pass vacuously.

| Name | Why it exists |
|---|---|
| `_ARTIFACT_WRITERS`, `_CACHE_WRITERS`, `_P13_CACHE_LITERAL` | the plan asks for "no write under `artifacts/`" and "no write on `"cache", "p13"`" while `h.meta.pool_dir` stays **readable**; both are therefore conjunctions (path token AND write call), and the needles need names |
| `_p14_line_hits(pred)` | the five simple greps are the same loop five times; expressed once, and it is what makes every failure a `path:index: line` string rather than a boolean |
| `_P14_TABLE` / `_p14_entry(path)` | the `test_p13_tau.jl` `SRC`/`LINES`/`CODE` triple generalized **per file**, because Phase 14's greps are lane-wide rather than about one probe; the index is what makes the SC2-c positional exemption expressible |
| `_p14_ood_state_range` raising on an unmatched `end` | an exemption that cannot be located must fail loudly. Returning `nothing` would silently protect *nothing* while looking like it protected something — the exact failure mode the allowlist self-consistency check exists for |
| the six-TIFF `isfile` assertions | the plan's positive `test/test_images` check is gated behind `run_p14_real_images.jl`, which does not exist yet. Asserting the six files are **present on disk** runs *now* and is what the later check will depend on |
| `@test occursin(":ghat_rho_true_scale", decide.code)` | 14-RESEARCH §D.5 requires the basis to be **recorded on the result**; a `ghat` call with no recorded basis makes the rescaling unauditable, which is half the point of SC2-d |
| `@test occursin("function p14_ood_state", fuse.code)` | the long form is what makes the exemption delimitable; asserting it turns an implicit trap into a stated requirement |
| `git status --porcelain -- spike/p13` | the plan asks for it in prose under "the shipped bundle … are read-only"; kept as its own assertion so the failure names the Phase-13 lane rather than the bundle |

### Not done, deliberately

- **`spike/test/test_p12_decoupling.jl` not edited.** Verified `git diff --exit-code` → 0. Each phase
  states its own claim in its own file; that file's scan cannot see `spike/p14/**` or `test_p14_*.jl`
  (`:78-88`), which is precisely the gap this plan closes.
- **`spike/test/runtests.jl` not touched and not invoked.** It exits 1 early at the Phase-4
  `SPEEDUP_GATE`, masking every later include block, so a green aggregate run would be no evidence.
  Every number above is from a per-file run. The new file also *asserts* that no Phase-14 source
  reaches for it.
- **`spike/p14/consts.jl` and `provenance.jl` not modified.** Both frozen (`a6c8258`, `473aac3`);
  `git diff --exit-code HEAD -- spike/p14` → 0. The two scratch appends in falsification checks 2
  and 3 were reverted with `git checkout --` and confirmed clean.
- **`STATE.md` / `ROADMAP.md` / `14-VALIDATION.md` not written** — the orchestrator owns the first
  two, and nothing here required amending the third.
- **`corpus/` never read or fetched.** `git diff --exit-code HEAD -- corpus` → 0.

## Assumption Drift (advisory)

**1. The allowlist's size — planned two, actual four**

- **Found during:** Task 1, pre-first-run.
- **Planned:** the plan's `<action>` names a two-entry `_P14_CORPUS_METADATA_ALLOWED` and calls
  adding an entry "a DECISION, not a fix".
- **Actual:** four entries. The two additions are *already-frozen* files whose corpus mentions are
  required by other running assertions, not new files asking for an exemption.
- **Why it matters:** the plan's framing implies the exemption surface is a forward-looking budget.
  It is partly a backward-looking fact: two entries were already spent before this file was written.
  Recorded so a reader does not read four entries as two decisions having been taken here.

**2. The plan assumed the file could name what it forbids in a testset title**

- **Found during:** Task 2, first run.
- **Planned:** the plan's acceptance criterion asks that the file *contain* the literal
  `runtests.jl`, alongside a check banning that literal from Phase-14 code.
- **Actual:** those two are compatible only if the literal lives in a **comment**. A `@testset`
  title is code. Same shape as 14-01's τ-absence contradiction, resolved the same way (the
  house run-time-assembly / comment-only idiom).
- **Why it matters:** it is the second time in this phase that "the guard must name the thing it
  forbids" produced a literal contradiction in a plan's own acceptance criteria. Worth knowing
  before the next guard is planned.

Nothing else drifted materially. No architectural question arose; no checkpoint was hit.

## Constraint Verification (all hold, all checked on disk)

| Constraint | Check | Status |
|---|---|---|
| `src/` untouched (D-01) | `git diff --exit-code HEAD -- src` + `git status --porcelain -- src` | exit 0, empty |
| `spike/Project.toml` + `Manifest.toml` byte-unchanged | `git diff --exit-code HEAD -- …` | exit 0 |
| No dependency added or changed | only `Test` + `Pkg` (both stdlib) | ok |
| Phase-16 sealed holdout not consumed | `corpus/` never read; `git diff --exit-code HEAD -- corpus` | exit 0 |
| `spike/p14/consts.jl` frozen at `a6c8258`, `provenance.jl` at `473aac3` | `git diff --exit-code HEAD -- spike/p14` | exit 0 |
| Phase-13 lane read-only | `git status --porcelain -- spike/p13` | empty |
| Corrections APPENDED, never overwritten | this file is a first write; `consts.jl` untouched | ok |
| `runtests.jl` not a gate | never invoked; asserted absent from Phase-14 code | ok |
| `STATE.md` / `ROADMAP.md` not modified by this executor | `git status --porcelain -- .planning/STATE.md .planning/ROADMAP.md` | **empty** |
| Files modified outside the allowed two | none — only `spike/test/test_p14_decoupling.jl` and this summary | ok |
| Every result-producing commit descends from the freeze | `git merge-base --is-ancestor a6c8258 HEAD` | **YES** |
| Main working tree, no worktrees | `git worktree` never invoked | ok |

## Known Stubs

None. The six `@info "… does not exist yet"` lines are **not** stubs — they are the house
arms-itself-later idiom, derived from the FILE rather than from a fixed list, and each announces a
*positive* half whose *negative* half is already live for every file that exists. Falsification
check 5a confirms one of them arms correctly the moment its file appears.

## Threat Flags

None — no network surface, no file write at runtime, no new trust boundary. The file shells out to
`git` read-only, reads Phase-14 source **text**, and queries the active `Pkg` environment. It opens
no image and loads no artifact.

The plan's nine register entries are discharged:

- **T-14-05** (env tampering) — diff + untracked-clean + exactly-sixteen + six by-name `!haskey`.
- **T-14-10** (`src/` productionization by stealth) — diff **and** `git status --porcelain`;
  falsification 1 dropped a scratch file and the untracked line is the one that failed.
- **T-14-04** (sealed holdout) — the four checks over a globbed surface, plus the allowlist
  self-consistency check so an exemption cannot silently protect nothing.
- **T-14-07** (the withdrawn `1/31 ≈ 0.032`) — banned from code; the general `1/(n+1)` exemption is
  deliberate and documented at the needle. Falsified (check 4).
- **T-14-08** (`verdict.flag == false` as in-distribution) — grep with a **positional** exemption,
  positive `:not_checked` twin. Falsified in **both** directions (5a/5b).
- **T-14-16** (unit mismatch vs τ) — per-line `ghat(` requirement plus the positive
  `ghat(patch_correlation(` and `:ghat_rho_true_scale` assertions. Falsified (check 4).
- **T-14-03** (`costes_p` consuming `NPE_MASTER_SEED`) — banned; positive check requires
  `P14_DEV_SEED` at every call site once `decide.jl` lands. Falsified (check 4).
- **T-14-17** (helper re-implementation) — five names banned, positive twin for the two that must
  actually be used. Falsified (check 4).
- **T-14-SC** (package installs) — nothing fetched, nothing installed; this plan **is** the
  enforcement.

## Note on in-file phase references

The file's comments cite plan numbers, wave numbers and D-/SC-/T- identifiers (e.g. "runs in WAVE 1",
"14-03", "D-03a", "Pitfall 7"). These are **research pre-registration and audit-trail references**,
not process bookkeeping: they are how this repository's frozen constants, amendments and threat
register are cited, they match `test_p12_decoupling.jl` / `test_p13_tau.jl` / `spike/p14/consts.jl`
verbatim in style, and the plan's own `must_haves` require several of them by name. They are
deliberate and consistent with the surrounding code.

## Self-Check

- `spike/test/test_p14_decoupling.jl` — **FOUND** (638 lines, ≥ 260 required; contains `_p14_sources`)
- commit `312329f` — **FOUND** (1 file, +392)
- commit `d802650` — **FOUND** (1 file, +246)
- both `key_links` from the plan frontmatter — **FOUND** (`("spike/p14", r"^.*\.jl$")` ×1,
  `git diff --quiet HEAD` ×8)
- `git status --porcelain` limited to the permitted paths — **clean**. The only other entries are
  `.planning/HANDOFF.json`, `.planning/config.json` and three untracked planning paths, none written
  or staged by this executor.

## Self-Check: PASSED

Plan 14-03's success criteria are met and **observed**, not inferred:

- **Every hard constraint in 14-CONTEXT.md "Constraints that bind this phase" is a running assertion
  in a Phase-14 file** — mapped one by one in the table above, with the one honest qualification
  named there (SC2-b's behavioural half is a wave-2 file).
- **The guard exists before the code it guards** — six of the seven modules do not exist yet, and
  falsification 5a shows the guard arming itself the moment one appears.

**Wave 0 + Wave 1 are complete.** `test_p14_consts.jl` (273), `test_p14_provenance.jl` (59) and
`test_p14_decoupling.jl` (42) all exit 0. Result-producing plans may now run.
