---
phase: 06-reproducible-demo-go-no-go-memo
plan: 01
subsystem: testing
tags: [demo, reproducibility, decoupling, spike, sbc, bayes-factor, ood, julia, jld2, random123]

# Dependency graph
requires:
  - phase: 04-npe-training-benchmark
    provides: frozen NPE (spike/npe/trained_npe.jld2) + harness draw_simulate_infer
  - phase: 05-validation-bundle-sbc-amortized-bf-ood
    provides: bf.jl/ood.jl validation surfaces, trained_ratio.jld2, *_report.jld2 headline numbers, consts.jl pre-registration
provides:
  - "spike/demo.jl — two-tier seeded end-to-end demo runner (DEMO-01)"
  - "in-script git-status decoupling assertion certifying src/ untouched (DEMO-02)"
  - "fast-tier NPE+BF/NRE+OOD chain proof, each with a twin-run bit-reproducibility assert (SC1)"
  - "--full subprocess dispatch to the reported-scale gates (run_sbc/run_bf/run_ood)"
affects: [06-02-go-no-go-memo, 07-productionization]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Twin-run reproducibility proof: seed both the Philox stream AND the global RNG, run the chain twice, assert bit-identical output"
    - "Fast/reported split reused: fixture-scale chain proof (VAL_FIX_SEED) vs --full reported gates (subprocess, VAL_MASTER_SEED)"
    - "isfile-guarded report loading with @warn fallback for gitignored artifacts on a fresh clone"

key-files:
  created:
    - spike/demo.jl
    - .planning/phases/06-reproducible-demo-go-no-go-memo/deferred-items.md
  modified: []

key-decisions:
  - "Seed the global RNG (Random.seed!) alongside the Philox stream so the stochastic NPE posterior-sampling layer is deterministic — required because posterior_for/sampleposterior thread no rng (infer.jl)"
  - "OOD fixture runs with_pp=false (mirrors run_ood.jl), so its twin-run proof is deterministic without touching the global RNG"
  - "Hard-assert only the committed frozen NPE at close; ratio net + reports are gitignored and handled by isfile guards / tiny-retrain fallback"
  - "Memo presence is reported but NOT hard-gated here — plan 06-02 authors and gates it"

patterns-established:
  - "Two-tier entry script (fast default + --full) mirroring 02_simulator_demo.jl idiom with a git-status decoupling assertion as step 0"

requirements-completed: [DEMO-01, DEMO-02]

# Metrics
duration: 40min
completed: 2026-07-03
---

# Phase 6 Plan 01: Reproducible Demo Runner Summary

**`spike/demo.jl` — a two-tier seeded runner that re-exercises the whole NPE+NRE/BF+OOD stack at fixture scale on VAL_FIX_SEED with per-layer twin-run bit-reproducibility asserts, certifies `src/` was never mutated via an in-script `git status` assertion, loads the POST-ITERATION headline numbers, and dispatches the reported-scale gates under `--full`.**

## Performance

- **Duration:** ~40 min
- **Started:** 2026-07-03T10:45Z
- **Completed:** 2026-07-03T11:21Z
- **Tasks:** 3
- **Files modified:** 1 created (spike/demo.jl) + 1 planning artifact (deferred-items.md)

## Accomplishments
- `spike/demo.jl` (312 lines) chains prior→simulator→summary→NPE→NRE/BF→OOD at fixture scale (64×64) from a fixed Random123 seed, each inference layer carrying a twin-run reproducibility `@assert` proving bit-identical output on the same seed (SC1).
- Step-0 decoupling assertion (`git status --porcelain -- src/ src/bayes.jl src/colocalization.jl`, argv form, no interpolation) proves `src/` untouched and fails fast before any compute (SC2/DEMO-02).
- Fast tier loads the three `*_report.jld2` headline numbers under `isfile` guards, annotated POST-ITERATION / Set 2, kept visually DISTINCT from the this-run chain-proof rows in the success-criteria table. Prints `bf max_abs_err = 10.9535` — confirming the memo's 10.95 figure (RESEARCH Assumption A1) before plan 06-02 quotes it.
- `--full` spawns run_sbc/run_bf/run_ood as independent subprocesses (try/catch keeps exit codes independent so one gate's throw cannot abort the others).
- Verified `runtests.jl` stays green with demo.jl present (OOD 27/27); locked files byte-unchanged.

## Task Commits

Each task was committed atomically:

1. **Task 1: demo.jl skeleton — header, includes, --full, decoupling assert, NPE chain proof** - `ab2fba3` (feat)
2. **Task 2: fast-tier BF/NRE + OOD fixture chains with twin-run asserts** - `a10a5aa` (feat)
3. **Task 3: headline-number load, success-criteria table, --full dispatch, self-assert close** - `c460d86` (feat)

## Files Created/Modified
- `spike/demo.jl` - Two-tier seeded end-to-end demo runner (DEMO-01) + in-script decoupling assertion (DEMO-02); composes the frozen harness/bf/ood surfaces, re-implements nothing.
- `.planning/phases/06-reproducible-demo-go-no-go-memo/deferred-items.md` - Logs D6-DEFER-01 (pre-existing flaky test_ood.jl SC5, out of scope).

## Decisions Made
- **Global-RNG seeding for the NPE twin-run.** `posterior_for` → `sampleposterior` threads no rng and draws the flow's base samples from Julia's global default RNG, so two runs on the same Philox seed still differed. Seeding `Random.seed!(seed)` alongside `val_rng(seed)` makes the whole stochastic chain deterministic — the honest meaning of "reproducible from a fixed seed" (SC1). (Rule 1 fix, see Deviations.)
- **OOD `with_pp=false`.** Mirrors run_ood.jl and avoids the non-finite-θ̂ PP crash (RESEARCH Pitfall 3); as a bonus the OOD twin-run needs no global reseed (the density/ROC path is Philox-only).
- **Memo not hard-gated here.** Presence is reported (`pending-06-02`); plan 06-02 authors the memo and wires the presence + content assertion.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] NPE twin-run reproducibility assert failed — global RNG not seeded**
- **Found during:** Task 1 (NPE fast-tier chain proof)
- **Issue:** The plan's twin-run `@assert r1.draws == r2.draws` failed: `draw_simulate_infer` → `posterior_for` → `sampleposterior` draws the normalizing-flow base samples from Julia's global default RNG (infer.jl threads no rng), so two calls on the same Philox seed produced different posterior draws.
- **Fix:** Seed the global RNG deterministically (`Random.seed!(seed)`) at the top of `run_npe_fixture` alongside the Philox `val_rng(seed)`. This makes the entire stochastic NPE chain reproducible from the fixed seed — the SC1 claim — without altering any frozen artifact or locked constant.
- **Files modified:** spike/demo.jl
- **Verification:** `julia --project=spike spike/demo.jl` → "NPE twin-run reproducible: true"; exits 0 CPU-only.
- **Committed in:** ab2fba3 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** The fix is necessary for the SC1 reproducibility proof to hold and stays within the plan's stated intent ("twin-run bit-reproducibility assert"). No scope creep; no locked file touched.

## Issues Encountered
- **Pre-existing flaky aggregate test (out of scope).** `runtests.jl` intermittently reports 1 failure in `test_ood.jl` SC5 (block-permute OR-flag). Root cause: the OOD fixture's `pp_mismatch_score` draws from the global RNG, which the aggregate suite never seeds, so it varies run-to-run. Proven independent of demo.jl (runtests.jl never loads demo.jl; the aggregate passed EXIT=0 with demo.jl both absent and present). Logged as D6-DEFER-01; suggested Phase-7 fix is to seed the global RNG in `_build_ood_fixture`. Not fixed (out of scope per SCOPE BOUNDARY).

## Verification Evidence
- `julia --project=spike spike/demo.jl` → exit 0; prints SC1/SC2/SC3 + DEMO-01/02/03 table; all three twin-run asserts PASS.
- Headline numbers loaded (Set 2): SBC ρ_true ECE=0.0164 (green); BF corr=0.9358, max|Δ logBF|=10.9535; OOD combined AUC=1.0, ID fire-rate=0.05.
- `git status --porcelain -- src/ src/bayes.jl src/colocalization.jl` → empty (asserted in-script).
- `git diff --exit-code spike/validation/consts.jl spike/Project.toml spike/Manifest.toml` → clean.
- `julia --project=spike spike/test/runtests.jl` → green with demo.jl present (OOD 27/27, EXIT 0).

## Next Phase Readiness
- **Ready for plan 06-02 (Go/No-Go memo):** demo.jl prints the confirmed `bf max_abs_err = 10.9535` (A1) and the full Set-2 headline numbers the memo must transcribe; the success-criteria table gives the SC1/SC2 evidence. The memo is authored + hard-gated in 06-02 (demo.jl already surfaces `MEMO_PATH` presence as `pending-06-02`).
- **No blockers.** src/ decoupling clean; locked files byte-unchanged.

## Self-Check: PASSED
- FOUND: spike/demo.jl
- FOUND: .planning/phases/06-reproducible-demo-go-no-go-memo/06-01-SUMMARY.md
- FOUND: .planning/phases/06-reproducible-demo-go-no-go-memo/deferred-items.md
- FOUND commits: ab2fba3, a10a5aa, c460d86

---
*Phase: 06-reproducible-demo-go-no-go-memo*
*Completed: 2026-07-03*
