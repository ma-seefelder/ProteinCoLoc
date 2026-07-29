---
phase: 13-three-hypothesis-amortized-bayes-factor
verified: 2026-07-29T18:30:00Z
status: gaps_found
score: 5.5/6 must-haves verified
has_blocking_gaps: false
overrides_applied: 0
gaps:
  - truth: "Known issues are honestly and completely recorded"
    status: partial
    severity: minor
    reason: >
      13-REPORT.md §14 point 5 claims only "Two Phase-13 SUMMARY files carry no `status:`
      frontmatter field (13-06-SUMMARY.md, 13-07-SUMMARY.md)". Independent verification found
      NINE files missing the field: 13-01 through 13-08 and 13-15-SUMMARY.md. The underlying
      fact pattern (a latent-risk, non-blocking gap) IS disclosed and IS recorded, but the
      report's own accuracy-check section undercounts it by 7 files. This does not affect any
      technical claim in the phase (all 16 plans are independently confirmed complete via
      ROADMAP.md [x] marks, git history, and artifact/test evidence), but it is a factual error
      inside the very section whose stated purpose is catching factual errors.
    artifacts:
      - path: ".planning/phases/13-three-hypothesis-amortized-bayes-factor/13-REPORT.md"
        issue: "§14 point 5 says 'Two' SUMMARYs lack status: frontmatter; nine actually do (13-01..13-08, 13-15)."
    missing:
      - "Correct §14 point 5 to name all nine files (13-01 through 13-08, 13-15), or replace the count with an accurate one."
---

# Phase 13: Three-Hypothesis Amortized Bayes Factor Verification Report

**Phase Goal:** Extend the amortized evidence network from two- to three-way model comparison —
colocalized / random / mutually-exclusive — so segregation becomes a first-class testable
hypothesis, replacing the fragile KDE+quadgk Bayes factor.

**Verified:** 2026-07-29
**Status:** gaps_found (1 minor, non-blocking gap; phase goal achieved)
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|---|---|---|
| 1 | A three-way evidence net exists, trains, and emits `log BF(coloc:random)` / `log BF(exclusion:random)` in one forward pass | ✓ VERIFIED | `spike/p13/net.jl:157-164` defines `ThreeWayEvidenceNet <: NeuralEstimator` with one shared trunk + one `Dense(64,2)` head, `(m::ThreeWayEvidenceNet)(Z) = m.head(m.trunk(Z))`. `spike/p13/three_way_net.jld2` (953,027 bytes) is present. `spike/p13/result.jl` defines `ThreeHypothesisColocResult <: AbstractColocResult` with inner constructor enforcing `random === 0.0` (structural zero, D-08). |
| 2 | The KDE+quadgk Bayes-factor reference is not depended on by the new path, and the report states honestly what was and wasn't compared | ✓ VERIFIED | `spike/Project.toml` (direct deps) contains no `QuadGK`/`KernelDensity` entry; both appear only transitively in `spike/Manifest.toml` (via `Distributions`/`CairoMakie`). `grep` for `KernelDensity\|QuadGK\|quadgk\|kde(` in `spike/p13/*.jl` returns only comments explaining why they are *not* used. `spike/test/runtests.jl:108-150` asserts (executable, ran green) `!haskey(Pkg.project().dependencies, "QuadGK"/"KernelDensity")`. Report §1 and "What was DECLINED" §1 state the KDE reference's `log(999)≈6.9` ceiling and its 34.75%/32.5% attrition, sourced from a prior spike (014), and explicitly declare it was declined as a validation target rather than silently dropped. |
| 3 | Pre-registration integrity: `consts.jl` has exactly two commits, both predating all Phase-13 results; `P13_ITERATION_ALLOWANCE` is 1 and unspent; no threshold edited after a result was seen | ✓ VERIFIED | `git log --oneline -- spike/p13/consts.jl` → exactly two commits: `c42cc8e` (2026-07-27 14:43:12) and `bfba6ac` (2026-07-27 21:46:26). Both predate all four result commits: `41a24f4` (2026-07-29 15:17), `e5f9569` (2026-07-29 16:13), `0943a02` (2026-07-29 16:36), `2112bed` (2026-07-29 17:07). `spike/p13/consts.jl:348` reads `const P13_ITERATION_ALLOWANCE = 1`; the file's own architecture (Tier-1/Tier-2 append-only guard blocks, `@assert P13_ITERATION_ALLOWANCE == 1`) makes in-place edits self-detecting. Report §10 shows the one pre-declared trigger evaluated mechanically and NOT fired (`iteration_trigger_fired = false`, deep-tail exclusion AUC 0.997257 is far above the 0.90 floor, so the trigger's own condition does not hold), and the allowance remains 1 of 1, unspent. |
| 4 | Decoupling: `src/`, `spike/Project.toml`, `spike/Manifest.toml` untouched by Phase 13 | ✓ VERIFIED | `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` → exit 0. `git log --oneline --all --since="2026-07-27" -- src/` → empty (no commit touched `src/` in the Phase-13 window). |
| 5 | The report puts BOTH the 6/6 gate pass and the real-image OOD flag in its opening, does not let the simulated pass imply real-microscopy validation, and discloses the epoch-4 checkpoint | ✓ VERIFIED | `13-REPORT.md` lines 4-6 (the `**Status:**` line, first content after the title) state both "cleared all six pre-registered criteria on one run (`6 pass / 6 total`, exit 0)" and "both real microscopy fixtures are simultaneously flagged out-of-distribution by both amortized detectors at ~2.5× their own in-distribution thresholds" in the same sentence — not in a trailing limitations section. The very next block, "## Scope of evidence" (lines 22-49), states explicitly: "The real-data arm can show the three-way verdicts behave sensibly on real microscopy; **it cannot show they are correct**." §2 discloses "The shipped artifact is an epoch-4 checkpoint of an early-overfitting run" with the measured best-validation risk (0.13625102 at epoch 4), the subsequent rise to 0.38632473, and the ~10× train/validation gap. Independently re-derived via `julia --project=spike spike/test/test_p13_correction.jl`, which reproduced the report's exact F5 numbers (coloc maxabs 0.5498217821206905, exclusion maxabs 0.37831267936810864, both against the 0.25 bar — 88 pass / 2 fail, matching §5's claimed "MISSED on both heads, and that is recorded, not repaired"). |
| 6a | Known issue: the Phase-4 `SPEEDUP_GATE` masks the Phase-13 include block in a full-suite run | ✓ VERIFIED | `.planning/phases/13-three-hypothesis-amortized-bayes-factor/deferred-items.md` §D-13-A records exactly this, with the measured speedup drift (92.50 → 83.97 → 68.35) and the consequence that per-file runs are the only available signal. Report §13 restates it and records that per-file Phase-13 tests were run directly. Independently re-ran `julia --project=spike spike/test/test_p13_labels.jl` (164/164 pass) and `julia --project=spike spike/test/test_p13_consts.jl` (137/137 pass) as spot checks; both green. |
| 6b | Known issue: SUMMARYs missing `status:` frontmatter are recorded, not silently dropped | ⚠ PARTIAL | Nine files genuinely lack a `status:` key (`13-01` through `13-08`, `13-15`; confirmed by direct grep and `head` inspection of each file's frontmatter). `13-REPORT.md` §14 point 5 DOES raise this issue, but undercounts it: it names only two files (`13-06`, `13-07`) and says "Two". See gap above. |

**Score:** 5.5/6 truths fully verified; 1 truth (6b) partially verified — the underlying fact is disclosed, but the report's own count is wrong by 7 files. This is a documentation-accuracy slip, not a defect in the delivered artifact, architecture, or gate result.

### Required Artifacts

| Artifact | Expected | Status | Details |
|---|---|---|---|
| `spike/p13/net.jl` | `ThreeWayEvidenceNet` architecture, two-head trunk | ✓ VERIFIED | Present, 30,715 bytes; defines the struct, `build_three_way_net`, `train_three_way`, `three_way_log_bf`, `save_three_way`. |
| `spike/p13/result.jl` | Result type carrying the three-way evidence triple | ✓ VERIFIED | Present, 20,804 bytes; `ThreeHypothesisColocResult <: AbstractColocResult` with `bayes_factor`, `log_bf_vs_random`, `is_ood`, `delta_rho` accessors; inner constructor asserts `random === 0.0`. |
| `spike/p13/three_way_net.jld2` | Trained net checkpoint | ✓ VERIFIED | Present, 953,027 bytes, dated 2026-07-29 (commit `41a24f4`). |
| `spike/p13/three_way_gate_report.jld2` | Persisted gate result artifact | ✓ VERIFIED | Present, 425,996 bytes (commit `e5f9569`). |
| `spike/p13/alpha_report.jld2` | Alpha-ladder result artifact | ✓ VERIFIED | Present, 53,391 bytes (commit `0943a02`). |
| `spike/p13/realimage_report.jld2` | Real-image qualitative-check artifact | ✓ VERIFIED | Present, 81,649 bytes (commit `2112bed`). |
| `spike/p13/tau_probe_report.jld2` | Measured tau curve | ✓ VERIFIED | Present, 15,343 bytes. |
| `spike/p13/consts.jl` | Two-tier pre-registration, append-only | ✓ VERIFIED | Present, 55,956 bytes; two-commit history confirmed independently. |
| `13-REPORT.md` | Phase results report | ✓ VERIFIED | Present, 1087 lines, substantive (not a stub); every headline number cross-checked against a persisted artifact or re-derived test run. |
| `13-SC2-AMENDMENT.md` | Pre-registered amendment, frozen before results | ✓ VERIFIED | Present; dated 2026-07-27, sha `df16e769...` recorded, and independently confirmed to predate all four Phase-13 result commits (2026-07-29). |

### Key Link Verification

| From | To | Via | Status | Details |
|---|---|---|---|---|
| `spike/p13/consts.jl` (Tier-1, `c42cc8e`) | Phase-13 results (`41a24f4`, `e5f9569`, `0943a02`, `2112bed`) | commit ordering | ✓ WIRED | All four result commits postdate both `consts.jl` commits by ~2 days; independently confirmed via `git log --format="%H %ai %s"`. |
| `spike/p13/net.jl` (`ThreeWayEvidenceNet`) | `spike/p13/run_three_way_gate.jl` | `build_three_way_net` / `train_three_way` calls | ✓ WIRED | The gate runner trains and evaluates the net type defined in `net.jl`; the persisted `three_way_net.jld2` metadata (`net_consts_sha`) is asserted equal to the live `consts.jl` sha at every reported run. |
| `spike/p13/consts.jl` | `spike/test/test_p13_*.jl` | literal-value assertions | ✓ WIRED | `test_p13_consts.jl` (137/137 pass, re-run) and `test_p13_labels.jl` (164/164 pass, re-run) both exercise constants and label logic defined in `consts.jl`/`labels.jl`. |
| `three_way_gate_report.jld2` | `13-REPORT.md` §6 | numbers re-read at report time | ✓ WIRED | Report states "every number...re-read from a persisted artifact...at report time rather than retyped"; §14 documents a discrepancy the report itself corrected against the artifact (deep-tail band ranking), which is affirmative evidence the artifact — not a SUMMARY — was actually re-read. |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|---|---|---|---|
| D-05 label trap regression holds (the counter-example the phase's central hypothesis-boundary decision rests on) | `julia --project=spike spike/test/test_p13_labels.jl` | 164/164 pass, incl. `three_way_label(0.8, 0.9) === RANDOM` (not EXCLUSION) and its mirror | ✓ PASS |
| Tier-1/Tier-2 pre-registration literal-value assertions hold | `julia --project=spike spike/test/test_p13_consts.jl` | 137/137 pass | ✓ PASS |
| F5 correction verification reproduces the report's claimed miss exactly | `julia --project=spike spike/test/test_p13_correction.jl` | 88 pass, 2 fail — coloc maxabs 0.549822 vs bar 0.25, exclusion maxabs 0.378313 vs bar 0.25, bit-for-bit matching §5's table | ✓ PASS (confirms report's honesty about the F5 miss, not a phase defect) |
| `src/`/`spike/Project.toml`/`spike/Manifest.toml` are byte-unchanged by Phase 13 | `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` | exit 0 | ✓ PASS |
| No commit touched `src/` in the Phase-13 window | `git log --oneline --all --since="2026-07-27" -- src/` | empty | ✓ PASS |

### Probe Execution

No `scripts/*/tests/probe-*.sh` convention or PLAN/SUMMARY-declared probes exist for this phase; Phase 13 uses Julia `@testset`/`@test` files as its verification mechanism, which were spot-checked above instead. Step 7c: SKIPPED (no probe-shaped entry points; this is a Julia test-file project, not a probe-script project).

### Requirements Coverage

ROADMAP.md records `**Requirements**: TBD` for Phase 13 — no formal `REQ-*` IDs are mapped to this phase in `.planning/REQUIREMENTS.md`. Coverage is instead expressed as the D-01 through D-16 decision log, all of which are cited against specific plans in the ROADMAP checklist and cross-referenced in the report. No orphaned requirement IDs found.

### Anti-Patterns Found

None. `grep` for `TBD|FIXME|XXX` and `TODO|HACK|PLACEHOLDER` across `spike/p13/*.jl` and `spike/test/test_p13_*.jl` returns zero debt markers; the only "placeholder" hits are legitimate prose describing a JLD2-reconstructed-placeholder failure mode the code guards against (`preconditions.jl`), not stub implementations.

### Human Verification Required

None. This phase's deliverable is a research artifact + report, verifiable entirely from git history, persisted `.jld2` artifacts, and re-runnable Julia test files — all of which were checked directly rather than deferred to a human.

### Gaps Summary

One minor, non-blocking gap: `13-REPORT.md` §14 point 5 undercounts the number of Phase-13 SUMMARY files missing a `status:` frontmatter field (says "Two", actual count is nine: `13-01` through `13-08` and `13-15`). The underlying issue — that some plan completions are attested only by the ROADMAP's `[x]` mark and the file's own prose rather than a machine-checkable `status:` field — IS disclosed in the report, and every one of those nine plans is independently corroborated by other evidence (commit history, artifact existence, or a later plan's dependency chain reading from it). This does not affect the phase's technical deliverable, the pre-registration integrity chain, the decoupling guarantee, or the gate result, all of which were independently re-verified above. Recommended fix: correct §14 point 5's count and file list, or delete the specific count and replace with "nine SUMMARY files (13-01 through 13-08, 13-15)".

**This looks intentional but incomplete rather than fabricated** — the report's author caught the pattern (a missing `status:` field being a latent risk) but scoped the grep too narrowly when writing the count. To accept this as-is without a follow-up correction, add to this file's frontmatter:

```yaml
overrides:
  - must_have: "Known issues are honestly and completely recorded"
    reason: "Undercount is a documentation slip in a non-gating section; all nine files are independently verified complete and the underlying risk pattern is disclosed"
    accepted_by: "<name>"
    accepted_at: "<ISO timestamp>"
```

---

*Verified: 2026-07-29T18:30:00Z*
*Verifier: Claude (gsd-verifier)*
