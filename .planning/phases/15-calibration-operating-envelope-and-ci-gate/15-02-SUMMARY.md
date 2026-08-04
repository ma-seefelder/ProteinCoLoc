---
phase: 15-calibration-operating-envelope-and-ci-gate
plan: 02
status: blocked
subsystem: infra
tags: [github-actions, ci, julia, glmakie, xvfb, headless, pkg-instantiate, artifacts]

# Dependency graph
requires:
  - phase: 07-productionization-conditional-on-go
    provides: "the co-resolution hard gate in test/runtests.jl and the committed Manifest.toml pinning NeuralEstimators 0.2.1 / Flux 0.16.10"
  - phase: 07-productionization-conditional-on-go
    provides: "ProteinCoLoc.estimator_for(8) and its hash-verified lazy-artifact resolution path in src/registry.jl"
provides:
  - "The repository's first .github/ directory and its first CI workflow"
  - "A PROVEN headless environment recipe (apt list + Julia pin + xvfb prefix) for the two-tier calibration CI"
  - "A measured minute budget: ~1 min warm / ~15 min cold fixed overhead per job"
  - "Confirmation of research assumptions A4 and A5 against real runner logs"
  - "The finding that the shipped grid_8 artifact is not fetchable by anyone — a hard blocker for 15-06 and 15-08"
affects: [15-06, 15-08, v2.0-release, reproducibility]

# Tech tracking
tech-stack:
  added: [GitHub Actions, actions/checkout@v4, julia-actions/setup-julia@v3, julia-actions/cache@v3, xvfb]
  patterns:
    - "De-risk the highest-risk CI unknown in an isolated throwaway-shaped job before building tiers on top of it"
    - "Measure environment requirements with an explicit control rather than inferring them from a passing run"

key-files:
  created:
    - .github/workflows/headless-smoke.yml
    - .planning/phases/15-calibration-operating-envelope-and-ci-gate/15-CI-BASELINE.md
  modified: []

key-decisions:
  - "Left the failing artifact step RED rather than softening it — the red run is an accurate report of a genuinely broken shipping path"
  - "Ran an explicit no-xvfb control to prove the prefix is load-bearing instead of inheriting an untested prefix into two more workflows"
  - "Did NOT publish a v2.0.0 release to make the run green — that is a public, irreversible act outside this plan's scope"
  - "Kept all ten apt packages despite necessity being untested; sufficiency is proven and narrowing on one observation would trade a known-good recipe for an untested one"

patterns-established:
  - "Workflow paths filter includes the workflow's own path, so committing it triggers its own proof run at zero cost on ordinary source pushes"
  - "Least-privilege CI: workflow-level permissions contents:read, push/workflow_dispatch only, no github.event.* interpolation into shell"

requirements-completed: [D-08, D-08a]

# Metrics
duration: 62min
completed: 2026-08-04
---

# Phase 15 Plan 02: Headless CI Smoke Baseline Summary

**The package is PROVEN to load on a headless ubuntu-24.04 runner under `xvfb-run` with the manifest-pinned NeuralEstimators 0.2.1 / Flux 0.16.10 — and the same runs prove the shipped `grid_8` net cannot be downloaded by anyone, because `Artifacts.toml` points at a `v2.0.0` release that does not exist.**

## Performance

- **Duration:** ~62 min (dominated by three real GitHub Actions runs)
- **Started:** 2026-08-04T17:28Z
- **Completed:** 2026-08-04T18:30Z
- **Tasks:** 2 of 2
- **Files modified:** 2 created, 0 modified

## Accomplishments

- Created the repository's first `.github/` directory with one pinned, least-privilege workflow.
- **Research assumption A5 CONFIRMED** — the named highest-risk CI unknown. `ProteinCoLoc loaded`
  printed on a headless runner. The whole CI tier is unblocked.
- **Research assumption A4 CONFIRMED** — `Pkg.instantiate()` on the committed manifest delivered
  `NeuralEstimators = 0.2.1` and `Flux = 0.16.10`, exactly what the co-resolution hard gate demands
  as literal strings. `[compat]` was not touched (D-08a).
- **Proved `xvfb-run` is load-bearing** with an explicit control run rather than assuming it: the
  bare load dies in GLFW's `__init__` with `"X11: The DISPLAY environment variable is missing"`.
- Measured a real minute budget for 15-06 / 15-08: **~1 min warm, ~15 min cold** fixed overhead,
  with `Pkg.instantiate()` at 11 m 41 s cold vs 3 s warm (a 537 MB depot cache is worth ~230× there).
- **Discovered the artifact-fetch blocker** that gates the rest of the phase's CI work.

## Task Commits

1. **Task 1: The isolated headless-load workflow** — `c5ff0f3` (ci)
2. **Task 2: Observe the run and record the measured baseline** — four commits, because the
   observation itself iterated:
   - `02689e6` (docs) — record the measured status in the workflow header after run 1
   - `738684b` (test) — add the temporary bare-load control
   - `9fbd861` (test) — retire the control, keeping its measured result
   - `1f729c3` (docs) — `15-CI-BASELINE.md`, the measured baseline

## Files Created/Modified

- `.github/workflows/headless-smoke.yml` — the isolated headless-load proof job and the reusable
  apt/xvfb/Julia-pin recipe the two later workflows copy from.
- `.planning/phases/15-calibration-operating-envelope-and-ci-gate/15-CI-BASELINE.md` — real run URLs,
  measured cold/warm durations, printed versions, the working apt list, per-step xvfb necessity, the
  A4/A5 verdicts, and the artifact finding.

## Verification Performed

Both automated `<verify>` commands were RUN and observed to pass:

- Task 1 YAML assertion (permissions, action pins, Julia 1.12.6, instantiate/xvfb/estimator_for
  present, no `pull_request_target`) → `ok`. Additionally checked and confirmed absent:
  `Pkg.test(`, `Pkg.resolve(`, `Pkg.add(`, `@v2`, `version: '1'`, `${{ github.event.`.
- Task 2 file/`A5`/`cold` assertion → `ok`.
- `git status --porcelain Project.toml Manifest.toml` → empty. `[compat]` untouched (D-08a).

Plan `<verification>` item 2 — "a run of *headless load smoke* is green" — is **NOT satisfied**.
All three runs are RED. The plan's stated alternative applies: the failure and its cause are
recorded. Steps 1-7 of the job are green in every run; only step 8 fails.

## Decisions Made

- **The failing artifact step was left failing.** Making it `continue-on-error`, or deleting it,
  would convert a true statement about a broken shipping path into a green badge. This project's
  discipline is measured-not-chosen and corrections-appended-not-hidden; a red run that accurately
  reports a broken artifact is the correct output.
- **No `v2.0.0` release was published.** Building and uploading the tarball would have turned the run
  green, but publishing a public release is irreversible, needs the author's decision, and is well
  outside a plan whose scope is one workflow plus a baseline document.
- **The xvfb control was worth an extra run.** The plan asks "whether `xvfb-run` was needed on each
  step". Without a control that question can only be answered by correlation, and two more workflows
  were about to inherit the answer.
- **All ten apt packages retained.** Sufficiency proven, necessity untested; explicitly stated as
  such in the baseline rather than implied.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Invalid YAML: a plain scalar `run:` containing `": "`**
- **Found during:** Task 1
- **Issue:** The `estimator_for(8)` step used an unquoted `run:` scalar whose Julia code contains
  `println("grid_8 bundle resolved: ", ...)`. YAML rejects `: ` inside a plain scalar, so the file
  did not parse at all.
- **Fix:** Converted that step to a `run: |` block scalar.
- **Verification:** The Task-1 assertion command parses the file and exits 0.
- **Committed in:** `c5ff0f3`

**2. [Rule 3 - Blocking] The header comment tripped the plan's own security assertion**
- **Found during:** Task 1
- **Issue:** The security-posture comment spelled out `pull_request_target` while explaining its
  absence, which failed the acceptance criterion "the string does not appear in the file". The
  criterion is right — a string assertion cannot distinguish prose from a trigger.
- **Fix:** Rewrote the comment to describe the privileged fork-trigger variant without the literal
  token.
- **Verification:** Assertion command exits 0.
- **Committed in:** `c5ff0f3`

**3. [Rule 2 - Missing Critical] Added an explicit no-xvfb control, then removed it**
- **Found during:** Task 2
- **Issue:** Runs 1 and 2 prove the load works *with* `xvfb-run`; they do not establish the prefix is
  needed. The plan's how-to-verify asks for exactly that, and 15-06/15-08 would have inherited an
  untested prefix.
- **Fix:** Added a temporary `continue-on-error` step running the identical load bare (run
  30937534659). It failed with `GLFW.GLFWError(65550, "X11: The DISPLAY environment variable is
  missing")`. Control removed; result kept as a comment beside the steps it constrains.
- **Verification:** Read the control step's log directly — its API `conclusion` reads `success`
  because of `continue-on-error`, which would have inverted the answer if trusted.
- **Committed in:** `738684b`, `9fbd861`

---

**Total deviations:** 3 auto-fixed (2 blocking, 1 missing-critical verification).
**Impact on plan:** No scope creep. Two were mechanical fixes to make the plan's own assertions pass;
the third strengthened the recipe the plan exists to produce.

## Assumption Drift (advisory)

**1. The shipped net was assumed fetchable; it is not.**
- **Found during:** Task 2
- **Planned:** 15-RESEARCH obstacle (3) states the release tarball "contains the amended_v2 net" and
  concludes CI "gets the shipped net with hash verification and works both locally (dev dir) and on a
  runner (download)". The plan's Task-1 step 8 is written on that basis.
- **Actual:** The release does not exist. `gh release list` returns only `v1.0.1` and
  `v1.0.0.compiled`; both the asset URL and the package-server mirror return HTTP 404, confirmed
  independently of the runner. The research verification compared a *local* tarball against the
  local dev dir, which cannot detect a missing publication.
- **Why it matters:** It is not only a CI problem. No third party can obtain the v2.0 network today,
  so the milestone's headline claim is currently unreproducible outside the author's checkout.

Advisory only — recorded, not gating. The blocker it implies is tracked below.

## Issues Encountered

- **BLOCKER (open): `grid_8` is not downloadable.** `Artifacts.toml` declares
  `git-tree-sha1 = 90e6b63a…`, `sha256 = 17162904…`, pointing at
  `.../releases/download/v2.0.0/grid_8.tar.gz`, which 404s. The files exist only in the git-ignored
  `artifacts/amended_v2/grid_8/`. **Resolution requires publishing a `v2.0.0` release whose tarball
  matches both recorded hashes** — a mismatched upload will still be rejected by
  `ensure_artifact_installed`, so this cannot be fixed carelessly.
- Run 2's `checkout` step took 1 m 31 s against 4 s in runs 1 and 3 — runner-side network variance,
  not a workflow property. Run 3 is quoted as the clean warm measurement, and the discrepancy is
  disclosed in the baseline rather than averaged away.
- `workflow_dispatch` is declared but not usable from the Actions UI/API yet, because GitHub only
  exposes dispatch for workflows present on the default branch. Runs were triggered by pushing to
  `gsd/p15-ci-smoke`, which the `paths` filter catches.

## User Setup Required

**Two items need the author, and one of them gates the rest of this phase's CI work.**

1. **Publish a `v2.0.0` release with `grid_8.tar.gz`** (BLOCKING for 15-06 and 15-08). Build it from
   `artifacts/amended_v2/grid_8/` and verify the upload's sha256 is
   `17162904c69b30e2a99a6dd934434553d08ebaadbfd55dea664b36efcf92b6d5` and its unpacked tree-sha1 is
   `90e6b63a8a234d067b407fefd7914f2ae4845448`.
2. **Delete the throwaway CI branch** `gsd/p15-ci-smoke` once its runs are no longer needed as
   evidence. It was pushed solely to trigger the proof runs; nothing was pushed to `main` or
   `gsd/v2.0-milestone`.

## Next Phase Readiness

**Ready:**
- The environment recipe is proven, not hoped: Ubuntu 24.04.4, the ten apt packages, Julia pinned to
  `1.12.6`, `julia-actions/cache@v3` at defaults, `Pkg.instantiate()`, and `xvfb-run -a` on every
  step that loads the package. 15-06 and 15-08 can copy it verbatim.
- The minute budget is real: ~1 min warm / ~15 min cold fixed overhead, on top of whatever the golden
  script costs.

**Blockers and concerns:**
- **15-06 and 15-08 cannot go green** until the `v2.0.0` release exists, because the D-09 golden is
  designed around `estimator_for(8)`. They can be written now; they cannot be proven now.
- A full sweep does not fit GitHub-hosted CI: `ubuntu-latest` is 4-core against the ~8.5 h budget's
  32-thread assumption, and the job cap is 6 h. The slow tier needs either a reduced scale or a
  self-hosted runner — decide this in 15-08, do not discover it there.
- **SHA-pinning the three third-party actions is recorded as a pre-v2.0-release hardening item**
  (threat T-15-06). Major tags are mutable refs.

## Self-Check: PASSED

- All three claimed files exist on disk: `.github/workflows/headless-smoke.yml`,
  `15-CI-BASELINE.md`, `15-02-SUMMARY.md`.
- All six claimed commits are present in `git log`: `c5ff0f3`, `02689e6`, `738684b`, `9fbd861`,
  `1f729c3`, `b2648b0`.
- STATE.md and ROADMAP.md were NOT modified — the orchestrator owns those writes.

---
*Phase: 15-calibration-operating-envelope-and-ci-gate*
*Completed: 2026-08-04*
