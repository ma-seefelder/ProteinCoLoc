---
phase: 11-registration-and-chromatic-uncertainty-as-latent
plan: 03
subsystem: testing
tags: [julia, regression-test, exact-equality, affine-warp, provenance, artifacts, decoupling, jld2]

# Dependency graph
requires:
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 01
    provides: "spike/test/fixtures/p11_stage6_golden.jld2 (the PRE-EDIT golden bytes), PHASE11_BASE_SHA = 17ebd1e, and p11_consts.jl's P11_STAGE6_EXACT / P11_STAGE6_TOLERANCE_FALLBACK / P11_FIXTURE_COUNTER / p11_rng"
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 02
    provides: "the composed single-AffineMap stage 6, chromatic_eps as the 8th theta field in both simulators, the defensive hasproperty read, _theta_tuple's length(v) >= 8 guard, and docs/amortized.md named limit #8 pinning 17ebd1e"
  - phase: 07-productionization-conditional-on-go
    provides: "the shipped amended_v2/grid_8 artifact, Artifacts.toml's git-tree-sha1 pin, registry.jl's three-branch lazy load, and the colocalization_amortized public entry point"
provides:
  - "spike/test/test_stage6_regression.jl — the standing, runnable SC1c/SC1d/SC1e exact-equality regression (51 assertions, no tolerance anywhere)"
  - "11-PROVENANCE.md — the recorded D-12 sha verification, the five F20c decoupling command outputs, the shipped-bundle regression, and the adopted narrower Q8 decoupling claim"
  - "Measured proof that the 8th theta column does not reach the shipped read surface (persisted D == 7, 7-element frozen BoundedThetaTransform)"
  - "A characterised, pre-existing Artifacts.toml tree-sha1 / installed-store mode-bit discrepancy (blobs byte-identical)"
affects: [11-05-probe, 11-07-research-trainer, 11-ladder-and-breakdown, 11-report, 11-11-interpretation-note]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Source-level assertion of a call-count invariant (`count(\"warp(\", stripped_source) == 1`) when the property being defended lives in the source, not in a value"
    - "A golden fixture that refuses itself: the regression test asserts the fixture LACKS the field the edit added, so a post-edit regeneration fails loudly instead of passing vacuously"
    - "Recording which resolution branch an integrity-checked loader actually took, rather than claiming the verification the plan assumed would run"

key-files:
  created:
    - spike/test/test_stage6_regression.jl
    - .planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-PROVENANCE.md
  modified: []

key-decisions:
  - "Exact `==` held everywhere; P11_STAGE6_TOLERANCE_FALLBACK was NOT used, so no tolerance deviation is carried into the Phase-11 report"
  - "estimator_for(8) resolved via branch 1 (installed artifact store), where _verify_tree_sha1 does NOT run; integrity was therefore verified independently and the branch recorded rather than assumed"
  - "The Q8 narrower decoupling wording is corrected from three edited src/amortized files to FOUR — train_npe.jl is in the measured diff (plan 11-02 deviation 1)"
  - "The Artifacts.toml tree-sha1 finding is RECORDED and DEFERRED, not fixed: repairing it would move Artifacts.toml and break the D-16 no-artifact-change boundary this plan audits"
  - "deferred-items.md was deliberately NOT edited — plan 11-04 runs in parallel and this plan declares only two files"

requirements-completed: [D-01, D-10, D-11, D-12, D-16]

# Metrics
duration: ~35min
completed: 2026-07-25
---

# Phase 11 Plan 03: Stage-6 Exact-Equality Regression and the D-12 Provenance Record Summary

**The D-10 regression is now a standing, runnable test rather than a one-off measurement: 51 assertions at exact `==` (no tolerance anywhere in the file) prove the composed stage 6 is byte-identical to the pre-edit golden at `chromatic_eps = 0` for both an 8-field and a legacy 7-field θ — and the recorded provenance document proves the pinned sha is the ε commit's first-parent lineage, that the four decoupling diffs are empty, and that the shipped `grid_8` bundle still loads with a persisted `D == 7` and still answers `colocalization_amortized(...)` on real images.**

## Performance

- **Duration:** ~35 min
- **Tasks:** 2
- **Files created:** 2 (no file modified, no file deleted)
- **Worktree base:** `2aae21f` → 2 commits

## Task Commits

1. **Task 1: stage-6 exact-equality regression test** — `8d7e45c` (test)
2. **Task 2: sha verification, decoupling commands, shipped-bundle regression** — `153e2b9` (docs)

## Accomplishments

- **Exact equality held, everywhere, first try.** All four golden Philox keys reproduce bit-for-bit
  through the full `simulate_pair` pipeline, for both `merge(θ7, (chromatic_eps = 0.0,))` and the
  stored 7-field θ passed unmodified; the composed warp matches the legacy `Translation` warp
  exactly across all four shift pairs. `P11_STAGE6_TOLERANCE_FALLBACK` was **not** used, so **no
  tolerance deviation is carried into the Phase-11 report**.
- **The fixture cannot silently invalidate itself.** The test asserts the golden θ has **no**
  `:chromatic_eps` key, has exactly 7 fields, ends at `:noise`, carries a 40-character sha, and rode
  `P11_FIXTURE_COUNTER` off `P11_DEV_SEED`/`P11_SALT`. A golden regenerated after the edit — the one
  failure mode that would make the whole regression pass vacuously — fails the first three of those.
- **The one-interpolation-pass requirement is asserted where it lives.** It is a property of the
  *source* (how many times `warp` is invoked), so it is asserted at source level:
  comment-lines stripped, `count("warp(", code) == 1`. The composed map is also pinned to the
  collapsed `AffineMap` form and asserted **not** to be a `ComposedTransformation`.
- **The shipped read surface is now measured, not argued.** RESEARCH §A5 argued from code that the
  8th θ column cannot reach the shipped bundle. It now has a regression: `arch.D == 7` read off
  disk while `theta_prior_bounds()` returns 8, `length(θzt.lo) == 7` with the **un-widened**
  `shift` bounds `[-1.0, 1.0]` frozen in, and `colocalization_amortized(...)` returning 2000 finite
  Δρ draws on the committed real TIFFs.
- **The Q8 decoupling claim is written down beside the ROADMAP text it replaces** — quoted verbatim,
  declared unavailable, and replaced with a narrower claim every clause of which has a recorded
  command behind it.

## Verification — observed results

Every command below was run and its **real output** observed.

| Check | Command | Observed |
|---|---|---|
| Stage-6 regression | `julia --project=spike spike/test/test_stage6_regression.jl` | **51 pass / 0 fail**, exit 0 (7.7 s) |
| — SC1c | nested testset | 3 pass |
| — SC1d(a) composed == legacy warp | nested testset | 4 pass |
| — SC1d(b) pipeline == golden | nested testset | 18 pass |
| — SC1e legacy 7-field θ | nested testset | 16 pass |
| — hostile `chromatic_eps` guards | nested testset | 9 pass |
| No tolerance anywhere | `grep -c 'isapprox\|atol\|rtol' spike/test/test_stage6_regression.jl` | `0` |
| Centre derived, not typed | `grep -c '128\.5' …` | `1` (only the SC1c type assertion) |
| Fixture consumed | `grep -c 'p11_stage6_golden.jld2' …` | `1` |
| AffineMap asserted | `grep -c 'AffineMap' …` | `4` |
| ε referenced | `grep -c 'chromatic_eps' …` | `14` |
| Root suite | `julia --project -t auto -e 'using Pkg; Pkg.test()'` | **passed**, exit 0, zero failures |
| First parent of the ε commit | `git rev-parse ca02b0e…^` | `a2ced26d93c7409a465aa50320e188868b3ab462` |
| Pinned sha predates ε | `git show 17ebd1e…:src/amortized/simulator.jl \| grep -c chromatic_eps` | `0` |
| Parent predates ε | `git show a2ced26…:src/amortized/simulator.jl \| grep -c chromatic_eps` | `0` |
| `simulator.jl` identical across the span | `git diff --stat 17ebd1e… a2ced26… -- src/amortized/simulator.jl` | empty |
| Pin quoted in the limit | `grep -c 17ebd1e… docs/amortized.md` | `2` |
| Parent quoted in the limit | `grep -c a2ced26… docs/amortized.md` | `1` |
| F20c-1 manuscript pipelines | `git diff --stat 17ebd1e…..HEAD -- src/bayes.jl src/colocalization.jl src/LoadImages.jl src/plot.jl` | empty |
| F20c-2 nothing outside `src/amortized/` | `git diff --name-only 17ebd1e…..HEAD -- src/ \| grep -v '^src/amortized/'` | empty |
| F20c-3 image bytes | `git diff --stat 17ebd1e…..HEAD -- test/test_images/` | empty |
| F20c-3 image digest | `git ls-files -s test/test_images/ \| sha256sum` | `eaee22f9185460fd910788026e7a33511f6de373f313459b96bb43f3fc1ef276` |
| F20c-4 artifact pin | `git diff --stat 17ebd1e…..HEAD -- Artifacts.toml` | empty |
| Shipped load | `ProteinCoLoc.estimator_for(8)` | loaded, `EstimatorBundle`, `grid = 8` (16.8 s) |
| Resolution branch | `Artifacts.artifact_exists(expected)` / `isdir(devdir)` | `true` / `false` → **installed artifact store** |
| Persisted flow marginals | `JLD2.load(npe_8.jld2)["arch"].D` | `7` (while `theta_prior_bounds()` is `8`) |
| Frozen θ transform | `length(bundle.θzt.lo)` / `.hi` | `7` / `7`; `lo[5:6] = -1.0`, `hi[5:6] = 1.0` |
| Real-image run | `colocalization_amortized(pos, neg, [1,2]; num_patches = 8)` | `AmortizedColocResult`, 2000 draws, **all finite**, 8.95 s |
| Δρ point estimate | mean of draws | `0.11021147163771093` (median `0.10957…`, 90% CI `[-0.0773, 0.3013]`) |
| log Bayes factor | `bayes_factor(res)` | `1.987931489944458` |
| OOD flag | `is_ood(res)` | `true` (density statistic `433.69`) |
| R3 / T-11-08 hinge | `_theta_tuple(v7)` → `simulate_pair` | accepted, 2×(128,128), all finite; `out7 == out8` at ε = 0 → `true` |
| Blobs vs pinned tree | SHA-256 of all four bundle files, store vs in-repo dev | **all four identical** |
| Store tree-sha1 | `Pkg.GitTools.tree_hash(store)` | `7d72b47e…` ≠ pin `90e6b63a…` (mode-bit only — see below) |
| Dev tree-sha1 | `Pkg.GitTools.tree_hash(devdir)` | `90e6b63a8a234d067b407fefd7914f2ae4845448` = pin |
| No deletions | `git diff --diff-filter=D --name-only HEAD~1 HEAD` (both commits) | empty |
| Clean tree | `git status --short` after each commit | empty |

**Not verified here (out of this plan's scope):** nothing in this plan trains a net, runs the D-06
probe, runs an SC2 ladder or evaluates coverage — so **no calibration, width, monotonicity,
vacuity or identifiability claim is made or checked**. The §3.5 real-image numbers are a **smoke
regression fingerprint on an OOD-flagged image pair**, explicitly not a colocalization result. No
download path was exercised (the installed store branch was taken), so nothing here is evidence
about the GitHub-Release fetch.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing critical] `estimator_for(8)` did not exercise `_verify_tree_sha1`, so integrity was verified independently**
- **Found during:** Task 2, Group 3
- **Issue:** The plan's acceptance criterion states that calling `ProteinCoLoc.estimator_for(8)`
  "exercises `_lazy_load_from_artifact!(8)` **including** `_verify_tree_sha1`". Measured: it does
  not. `_verify_tree_sha1` is called only on branch 2 (in-repo dev fallback,
  `registry.jl:141-142`). On this machine `Artifacts.artifact_exists(expected)` is `true`, so
  branch 1 (installed content-addressed store) is taken and no verification runs at load time.
  Reporting a "tree-sha1-verified load" would have been a claim about code that did not execute —
  precisely the failure the plan's own "never claim a download-verified load that did not happen"
  instruction forbids.
- **Fix:** Recorded the branch taken explicitly in `11-PROVENANCE.md` §3.1, stated plainly that
  `_verify_tree_sha1` did not run there, and supplied independent integrity evidence instead:
  SHA-256 of all four bundle blobs compared store-vs-pinned-tree (all identical), plus direct
  invocation of `_verify_tree_sha1` on both directories.
- **Files modified:** `11-PROVENANCE.md` (§3.1, §3.6, §6)
- **Verification:** the branch-detection output and both `_verify_tree_sha1` invocations are
  transcribed verbatim in the document.
- **Committed in:** `153e2b9`

**2. [Rule 2 - Missing critical] The Q8 narrower decoupling wording names three `src/amortized/` files; the measured diff has four**
- **Found during:** Task 2, Group 2 (command 2)
- **Issue:** The plan's prescribed wording is *"only `src/amortized/` is edited, and within it only
  `simulator.jl`, plus a defensive read in `ood.jl` and the theta-arity constants in
  `datagen.jl`"*. `git diff --name-only 17ebd1e…..HEAD -- src/` returns **four** paths:
  `datagen.jl`, `ood.jl`, `simulator.jl` and **`train_npe.jl`**. Adopting the three-file wording
  verbatim would have made the phase's central honesty claim inaccurate on its face — checkable in
  one command by any reviewer.
- **Fix:** Adopted the narrower claim with the enumeration corrected to four, and wrote §4.2
  explaining *why* the fourth file is there (plan 11-02's recorded deviation 1: after the prior
  gained an 8th row the src trainer accepted no θ arity at all, and the repair derives the flow's
  marginal count from `size(θtr, 1)`), together with the measured evidence that `NPE_D`,
  `build_estimator`'s default `D` and `architecture.jl` are untouched — so the shipped read surface
  is unaffected.
- **Files modified:** `11-PROVENANCE.md` (§4.2)
- **Verification:** the unfiltered `git diff --name-only` output is transcribed in §2.2.
- **Committed in:** `153e2b9`

**3. [Rule 3 - Blocking, minor] The golden's `phase11_base_sha` is a `SubString{String}`, not a `String`**
- **Found during:** Task 1
- **Issue:** The plan directs the test to assert the fixture's `phase11_base_sha` "is a 40-character
  string". `JLD2` reads it back as `SubString{String}` (the capture script stored the output of
  `readchomp`), so `isa String` is `false` and a literal reading of the instruction would have
  produced a test that fails on a perfectly valid fixture.
- **Fix:** Asserted `sha isa AbstractString` together with `length(sha) == 40` — the substantive
  property the plan is after, without pinning an incidental storage type.
- **Files modified:** `spike/test/test_stage6_regression.jl`
- **Verification:** SC1d(b) testset passes 18/18.
- **Committed in:** `8d7e45c`

### Discovered, recorded, NOT fixed (out of scope)

**4. [Pre-existing] The `Artifacts.toml` `git-tree-sha1` pin does not match the installed store copy**

Not a deviation from the plan and not caused by Phase 11 — surfaced only because deviation 1 forced
an independent integrity check. Recorded in full in `11-PROVENANCE.md` §3.6.

- **What:** `Pkg.GitTools.tree_hash` over the installed store dir gives `7d72b47e…`; the
  `Artifacts.toml` pin is `90e6b63a…`. The in-repo dev tree hashes to exactly the pin.
- **Root cause, measured:** a git tree hash covers the file **mode**, and on Windows
  `Pkg.GitTools.gitmode` derives it from `Sys.isexecutable`, which tracks the write bit. The
  writable in-repo copies hash as `100755`; Julia's artifact install makes store files read-only, so
  they hash as `100644`. **All four blobs are byte-identical (SHA-256 verified).**
- **Why it is latent:** branch 1 is taken on this machine and never calls `_verify_tree_sha1`. If
  the store copy were removed and the dev dir present, branch 2 would run the check — and it
  **passes** there. The configuration that would actually bite is a non-Windows machine
  reproducing the pin from a writable tree.
- **Why not fixed here:** any repair moves `Artifacts.toml`, which would break the D-16
  "no artifact change" boundary and invalidate the §2.4 byte-unchanged claim this plan exists to
  record. Deliberately left alone.
- **Handover:** this finding belongs in the phase's `deferred-items.md`. It was **not** written
  there by this plan because plan 11-04 runs in parallel in a sibling worktree and this plan
  declares only two files; the orchestrator should carry §3.6 across at merge time.

---

**Total deviations:** 3 auto-fixed (1 blocking-minor, 2 missing-critical) + 1 pre-existing finding
recorded and deferred. **No Rule 4 (architectural) situation arose; no checkpoint was hit.**
**Impact on scope:** none. Deviations 1 and 2 both make the recorded claims *more* accurate rather
than weaker; deviation 3 is a type technicality. No `src/` file, no dependency, no artifact and no
public API was touched by this plan.

## Assumption Drift (advisory)

**1. "Exercising the shipped load path exercises its integrity check."**
- **Found during:** Task 2, Group 3
- **Planned:** the plan (and the T-11-11 threat-register row) treats `estimator_for(8)` as running
  `_verify_tree_sha1`, i.e. as an end-to-end integrity exercise.
- **Actual:** the shipped load has three branches and only one of them verifies. The common
  developer configuration — artifact already in the store — verifies **nothing** at load time.
- **Why it matters:** any future statement of the form "the shipped bundle is tree-sha1-verified on
  every load" is false as written. It is verified at *install* time, and on the dev-fallback branch.
  The Phase-11 report should say which.

**2. "Byte-identical contents imply an identical git tree hash."**
- **Found during:** Task 2, §3.6
- **Planned:** implicitly, that a tree-sha1 mismatch would mean tampered or mismatched bytes (the
  T-11-11 framing).
- **Actual:** a git tree hash also covers the file mode, and on Windows that mode is synthesised
  from the write bit. Identical bytes in a read-only store hash differently from the same bytes in a
  writable checkout.
- **Why it matters:** a future reader hitting the §3.6 error message would reasonably conclude the
  artifact was corrupted. It is not. The distinction between a content failure and a mode failure
  needs to survive into whatever release-engineering pass picks this up.

**3. "The `src/amortized/` blast radius is the three files the research enumerated."**
- **Found during:** Task 2, Group 2
- **Planned:** `simulator.jl` + `ood.jl` + `datagen.jl`.
- **Actual:** four files; `train_npe.jl` joined via 11-02's deviation 1.
- **Why it matters:** 11-02's summary already flagged this, but the plan text for 11-03 had not
  absorbed it. The claim wording is the one place where a stale enumeration becomes a false public
  statement, so it is corrected at the source (§4.2) rather than only in a summary.

## Issues Encountered

- **`_verify_tree_sha1` never ran on the branch taken** — diagnosed by reading
  `Artifacts.artifact_exists(expected)` directly rather than assuming, then confirmed by invoking
  the guard manually on both candidate directories. Resolved as deviation 1.
- **The manual guard invocation on the store directory ERRORED** with a genuine-looking integrity
  failure. Fully diagnosed to a file-mode difference with byte-identical contents before anything
  was written down; recorded as finding 4 rather than reported as tampering.
- **No exact-equality assertion failed at any point.** The tolerance escape hatch was never
  approached, so §A3's measured result is confirmed rather than merely re-asserted.

## Known Stubs

None. Every assertion in `test_stage6_regression.jl` is a real comparison against real data — there
is no placeholder value, no mock, no empty container, no TODO, and no skipped testset. The entry-guard
testset deliberately exercises **both** directions (three rejections and two acceptances), so a guard
stuck in the always-throw state would be caught.

## Threat Flags

None — this plan adds no network endpoint, no auth path, no schema at a trust boundary, and its only
new file access is a read of the already-committed golden fixture and the already-committed TIFFs.
The registered threats were handled as planned:

| Threat ID | Handling | Evidence |
|---|---|---|
| T-11-09 (repudiation — wrong pinned sha) | `git rev-parse "<EPS>^"` and both `git show <SHA>:… \| grep -c chromatic_eps` recorded verbatim | §1.1-1.3; all three agree |
| T-11-11 (shipped artifact integrity) | `Artifacts.toml` asserted byte-unchanged; the resolution branch **recorded, not assumed**; blob-level SHA-256 comparison supplied where `_verify_tree_sha1` did not run | §2.4, §3.1, §3.6 |
| T-11-08 (7-row vector into an 8-field `simulate_pair`) | explicit regression: `_theta_tuple(v7)` → `simulate_pair` accepted, and `out7 == out8` at ε = 0 | §3.7 |
| T-11-13 (exact equality silently relaxed) | `grep -c 'isapprox\|atol\|rtol'` = `0`; the tolerance constant was never invoked | verification table |
| T-11-14 (reading outside the committed image set) | only `test/test_images/` read; `corpus/` and its `sealed_holdout` split untouched | §3.5 |
| T-11-SC (package installs) | none attempted; no dependency added to either environment | `git status --short` empty after each commit |

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- **Ready:** the D-10 regression is now a standing test any later wave can re-run in ~8 s. Plans
  11-05 onward can edit the probe/ladder/trainer surfaces knowing a stage-6 behavioural change will
  be caught immediately rather than at report time.
- **Ready:** the Phase-11 report has its §F20c command outputs, its pinned sha, its `test_images`
  digest and its `skip_if_done` caveat sentence already recorded and quotable from
  `11-PROVENANCE.md`.
- **Constraint for the report:** the §3.5 real-image numbers must be quoted **only** as a smoke
  fingerprint, never as a colocalization result — the shipped OOD flag fired on that pair.
- **Handover for the orchestrator:** carry `11-PROVENANCE.md` §3.6 into the phase
  `deferred-items.md` at merge time (not written here, to respect the parallel-worktree file
  boundary with plan 11-04).
- **Unspent:** `P11_ITERATION_ALLOWANCE = 1` is still unspent; nothing in this plan trained,
  evaluated or gated anything.

## Self-Check: PASSED

- `spike/test/test_stage6_regression.jl` exists on disk (197 lines).
- `.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-PROVENANCE.md` exists on
  disk (486 lines).
- Both task commits exist in `git log --oneline`: `8d7e45c`, `153e2b9`.
- Neither commit deleted a tracked file (`git diff --diff-filter=D --name-only HEAD~1 HEAD` empty
  for both).
- `git status --short` is empty after each commit; no untracked residue.

---
*Phase: 11-registration-and-chromatic-uncertainty-as-latent*
*Completed: 2026-07-25*
