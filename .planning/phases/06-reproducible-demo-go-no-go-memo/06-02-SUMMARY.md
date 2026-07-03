---
phase: 06-reproducible-demo-go-no-go-memo
plan: 02
subsystem: documentation
tags: [memo, go-no-go, decision, spike, demo, pre-registration, data-snooping, ship-gate, julia]

# Dependency graph
requires:
  - phase: 06-01
    provides: spike/demo.jl two-tier runner with pending-06-02 memo presence row + confirmed bf max_abs_err=10.9535
  - phase: 05-validation-bundle-sbc-amortized-bf-ood
    provides: 05-04-SUMMARY.md Set-1 pre-registered FAIL numbers; on-disk *_report.jld2 Set-2 post-hoc numbers
provides:
  - "06-GO-NO-GO-MEMO.md — 2-3 page Clean Go decision memo (DEMO-03)"
  - "demo.jl DEMO-03 memo presence + content-keyword hard gate (occursin on Clean Go / falsification / 0.8861 / 0.9358 / ship-gate / Phase 11 / Phase 13)"
affects: [07-productionization]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Machine-gated decision memo: demo.jl reads the memo read-only and @assert occursin each required keyword, so the spike cannot report closed without the verdict/falsification/both-number-set anchors present"
    - "Runtime-confirmed headline figure: BF max|Δ logBF|=10.9535 quoted in the memo only after confirming it against demo.jl's printed bf_report[\"max_abs_err\"] (A1 resolved)"

key-files:
  created:
    - .planning/phases/06-reproducible-demo-go-no-go-memo/06-GO-NO-GO-MEMO.md
  modified:
    - spike/demo.jl

key-decisions:
  - "Quoted BF max|Δ logBF|=10.9535 (demo.jl's exact printed value) rather than STATE.md's rounded 10.95 prose — the BF FIGURE VERIFICATION constraint; both are noted (10.9535 ≈ 10.95)"
  - "Memo content gate uses explicit occursin @asserts for the 7 required keywords plus a MEMO_OK conjunction, so a memo that exists but omits required content still FAILS DEMO-03"
  - "Added lowercase 'falsification condition' phrasing so the case-sensitive keyword gate and the plan's grep both match (heading is 'Falsification condition')"

requirements-completed: [DEMO-03]

# Metrics
duration: 20min
completed: 2026-07-03
---

# Phase 6 Plan 02: Go/No-Go Memo + DEMO-03 Gate Summary

**Authored the 2–3-page Go/No-Go memo (`06-GO-NO-GO-MEMO.md`) rendering a defensible Clean Go — reporting BOTH number sets honestly (Set 1 pre-registered all-FAIL verbatim from 05-04-SUMMARY.md; Set 2 post-hoc with the named data-snooping exposure + three mitigations), an explicit falsification condition, the D-05 independent confirmation ship-gate, and the Phase 8–16 DAG with Phase 11 ∥ 13 lead — then wired a machine-checked presence + content gate into `spike/demo.jl` so DEMO-03 now reports a real PASS.**

## Performance

- **Duration:** ~20 min
- **Started:** 2026-07-03T11:25Z
- **Completed:** 2026-07-03T11:45Z
- **Tasks:** 2
- **Files:** 1 created (memo) + 1 modified (demo.jl)

## Accomplishments

- **`06-GO-NO-GO-MEMO.md` (243 lines):** Clean Go verdict that states explicitly the literal pre-registered gates did NOT pass as written, justifying the Go by reading each residual failure to a characterized non-method cause (χ²-over-power at M=2000, residual data-scale gap, clamped-KDE baseline tail artifact).
  - **Set 1 (primary, pre-registered, all FAIL)** transcribed verbatim from 05-04-SUMMARY.md: SBC ρ_true KS p=2.81e-40 / χ² 5.27e-146 / ECE 0.189 red, Δρ ECE 0.192 red, 6/8 ECE-green; BF corr 0.8861 / max|Δ logBF| 3.598 / log_prior_odds −0.0294; OOD pooled AUC 0.730 FAIL, noise family 0.0 blind spot, ID fire-rate 0.05. Labeled NOT reproducible from current artifacts.
  - **Set 2 (post-hoc)** from on-disk reports + STATE.md: SBC ECE green all 8 (ρ_true 0.0164), BF corr 0.9358 + clamped-KDE diagnosis, OOD combined pooled AUC 1.0. Data-snooping exposure named plainly (frozen net retrained twice after the pre-registered result was seen) + three mitigations (consts.jl byte-unchanged vs `e9c91d3`, disjoint DEV seed `0xDE7C0DE`, single confirmatory VAL run each).
  - **Falsification condition (D-02):** SBC ECE staying red OR BF mid-range disagreement OR an OOD family undetectable by any summary-orthogonal channel — none fired; each becomes a Phase-7 pre-registered pass criterion.
  - **Ship-gate (D-05):** commits Phase 7 to an independent fresh-seed re-pre-registered SBC/BF/OOD confirmation run before amortized inference ships into `src/`.
  - **SBC framing:** "calibrated under the simulator" paired with OOD in the same breath (carried-forward hard requirement).
  - **Thesis (NPE-03):** >100× faster than ADVI (median 325× forward-pass, ~16× full workload) at comparable RMSE, ρ recovery corr 0.983, ADVI ~0.5 s/pair caveat.
  - **Build-out DAG (D-07/D-08):** Wave A {8,9,10}, Wave B {11 ∥ 13, then 12}, Wave C {14,15}, Wave D {16}; Phase 11 ∥ Phase 13 lead, Phase 12 follows 11; hierarchy/3D/multi-channel longer horizon.
- **`spike/demo.jl` DEMO-03 gate:** `@assert isfile(MEMO_PATH)` plus explicit `occursin` content asserts for "Clean Go", "falsification", "0.8861" (Set 1), "0.9358" (Set 2), "ship-gate", "Phase 11", "Phase 13"; SC3 + DEMO-03 table rows now report real PASS/FAIL (was `pending-06-02`). demo.jl re-run exits 0 with DEMO-03 PASS.

## BF Figure Verification (critical constraint)

Before quoting the BF figure, ran `julia --project=spike spike/demo.jl` and confirmed it prints `max_abs_err(max|Δ logBF|)=10.9535`. The memo quotes **10.9535** (noted ≈10.95), matching demo.jl's runtime output exactly — resolving RESEARCH Assumption A1. corr=0.9358 was already read directly from bf_report.jld2.

## Task Commits

1. **Task 1: Author the Go/No-Go memo** — `ef30dfc` (docs)
2. **Task 2: Wire memo presence + content gate into demo.jl** — `530b287` (feat)

## Files Created/Modified

- `.planning/phases/06-reproducible-demo-go-no-go-memo/06-GO-NO-GO-MEMO.md` (created) — the 2–3-page Clean Go decision memo.
- `spike/demo.jl` (modified) — MEMO_KEYWORDS/MEMO_CONTENT/MEMO_OK/DEMO03_PASS block, updated SC3 + DEMO-03 table rows, hard `@assert isfile` + `occursin` content gate at close.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Julia parse error — escaped quotes inside string interpolation**
- **Found during:** Task 2 (first demo.jl re-run after adding the memo gate)
- **Issue:** `$(join(MEMO_MISSING, \", \"))` inside an `@assert` message string is a `ParseError` in Julia — `\"` is not valid inside a `$(...)` interpolation nested in a string literal.
- **Fix:** Simplified the message to `$(MEMO_MISSING)` (interpolates the vector directly). The seven explicit per-keyword `occursin` asserts already give a precise failure message; the conjunction assert only needs to report the residual missing set.
- **Files modified:** spike/demo.jl
- **Verification:** `julia --project=spike spike/demo.jl` → exit 0, DEMO-03 PASS.
- **Committed in:** 530b287 (Task 2 commit)

**2. [Rule 1 - Bug] Lowercase "falsification" keyword absent (case-sensitive gate would fail)**
- **Found during:** Task 1 (memo content self-check)
- **Issue:** The memo used "Falsification condition" (capital F), "falsifiable", "falsifiers" — but the demo.jl content gate and the plan's verify grep both match the case-sensitive literal "falsification". The gate would have FAILED DEMO-03.
- **Fix:** Added an explicit lowercase "falsification condition" phrase in §6 prose.
- **Files modified:** 06-GO-NO-GO-MEMO.md
- **Verification:** `grep -F "falsification"` → present; demo.jl DEMO-03 PASS.
- **Committed in:** ef30dfc (Task 1 commit)

**Total deviations:** 2 auto-fixed (1 blocking parse error, 1 keyword-match bug). No locked file touched; no threshold or number changed.

## Verification Evidence

- `julia --project=spike spike/demo.jl` → exit 0; SC3 and DEMO-03 rows report PASS; prints `bf_report: corr=0.9358 max_abs_err(max|Δ logBF|)=10.9535`.
- Memo content: all required anchors present (Clean Go, 0.8861, 3.598, 0.730, 2.81e-40, 0.189, 0.0164, 0.9358, combined pooled AUC = 1.0, falsification, e9c91d3, 0xDE7C0DE, single confirmatory VAL run, ship-gate, Phase 11, Phase 12, Phase 13, calibrated under the simulator, 10.9535).
- `git diff --exit-code spike/validation/consts.jl spike/Project.toml spike/Manifest.toml` → clean.
- `git status --porcelain -- src/ src/bayes.jl src/colocalization.jl` → empty.
- `julia --project=spike spike/test/runtests.jl` → green (EXIT 0; OOD 27/27, Phase-9 45/45).

## Known Stubs

None. The memo is complete prose; demo.jl reads real on-disk artifacts and the committed memo. No placeholder values, no unwired data paths.

## Threat Flags

None new. All surface is offline/CPU-only/local-file; demo.jl's memo read is read-only (T-6-02 mitigated); Set 1 transcribed verbatim + labeled primary, Set 2 explicitly labeled post-hoc (T-6-04 mitigated); no package installs (T-6-03/T-6-SC).

## Next Phase Readiness

- **Phase 7 (productionization, conditional on Go):** the memo renders a Clean Go and commits Phase 7 to the D-05 independent confirmation ship-gate (fresh-seed, re-pre-registered SBC/BF/OOD) as the hard precondition before amortized inference integrates into `src/`. Two hardening items ride inside it: re-enable the OOD PP channel (with_pp=false today), and an optional non-clamped BF baseline so D-08b is testable.
- **Phase 6 complete:** both plans done (06-01 demo runner, 06-02 memo + gate); DEMO-01/02/03 all PASS; src/ provably untouched; locked files byte-unchanged.

## Self-Check: PASSED

- FOUND: .planning/phases/06-reproducible-demo-go-no-go-memo/06-GO-NO-GO-MEMO.md
- FOUND: spike/demo.jl (modified)
- FOUND commits: ef30dfc, 530b287

---
*Phase: 06-reproducible-demo-go-no-go-memo*
*Completed: 2026-07-03*
