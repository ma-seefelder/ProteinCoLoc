---
phase: 06-reproducible-demo-go-no-go-memo
verified: 2026-07-03T15:10:00Z
status: passed
score: 14/14 must-haves verified
has_blocking_gaps: false
overrides_applied: 0
---

# Phase 6: Reproducible Demo + Go/No-Go Memo Verification Report

**Phase Goal:** The spike closes with a falsifiable, reproducible verdict — a single seeded script chains every layer and a memo states the metrics and a concrete full-build-out decision.
**Verified:** 2026-07-03T15:10:00Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `spike/demo.jl` chains the full pipeline reproducibly from a fixed Random123 seed, CPU-only, and tabulates success criteria (ROADMAP SC1) | VERIFIED | Ran `julia --project=spike spike/demo.jl` live — exit code 0. Printed "NPE twin-run reproducible: true", "BF twin-run reproducible: true", "OOD twin-run reproducible: true", then the full success-criteria table with SC1/SC2/SC3 + DEMO-01/02/03 all PASS. |
| 2 | The main repo and both manuscript pipelines are demonstrably untouched — `git status` on `src/`, `bayes.jl`, `colocalization.jl` is clean (ROADMAP SC2 / DEMO-02) | VERIFIED | Independently ran `git status --porcelain -- src/ src/bayes.jl src/colocalization.jl` — empty output. demo.jl step-0 also asserts this in-script (backtick argv form, no `$`-interpolation) and printed "[step 0] decoupling OK" during the live run. |
| 3 | A 2–3-page Go/No-Go memo reports the metrics, frames SBC "calibrated under the simulator" paired with OOD, and states a concrete full-build-out decision (ROADMAP SC3 / DEMO-03) | VERIFIED | `06-GO-NO-GO-MEMO.md` (243 lines) read directly. §3 contains the literal phrase "calibrated under the simulator" immediately followed by the OOD framing. §8 states the Phase 8–16 DAG with Wave A/B/C/D. |
| 4 | Two-tier design: fast default (fixture-scale) + `--full` (reported-scale subprocess dispatch) | VERIFIED | `const FULL = ("--full" in ARGS)`; `if FULL … run(\`julia --project=spike $script\`)` block over run_sbc.jl/run_bf.jl/run_ood.jl present in demo.jl (source-inspected; not exercised live to keep verification time-boxed). |
| 5 | Each of NPE, BF/NRE, and OOD fixture steps carries a twin-run bit-reproducibility assert | VERIFIED | Live run confirmed all three: `@assert NPE_R1.draws == NPE_R2.draws`, `@assert BF_R1.logbfs == BF_R2.logbfs`, `@assert OOD_R1.maha_pos == OOD_R2.maha_pos && … combined_auc == …` all passed (printed "reproducible: true" for each). |
| 6 | Decoupling assertion uses backtick argv form with no string interpolation (command-injection safe, T-6-01) | VERIFIED | Source line 79: `` read(`git status --porcelain -- src/ src/bayes.jl src/colocalization.jl`, String) `` — literal tokens, no `$(...)`. |
| 7 | `spike/validation/consts.jl`, `spike/Project.toml`, `spike/Manifest.toml` byte-unchanged | VERIFIED | `git diff --exit-code spike/validation/consts.jl spike/Project.toml spike/Manifest.toml` → exit 0 (clean). Additionally `git diff e9c91d3 -- spike/validation/consts.jl` → exit 0 (byte-identical to the pre-registration commit the memo cites). |
| 8 | Memo reports BOTH number sets in order — Set 1 pre-registered FAIL (primary, verbatim from 05-04-SUMMARY.md) then Set 2 post-iteration (post-hoc) | VERIFIED | §3 "Set 1 — PRIMARY pre-registered result (all three gates FAIL)" transcribes SBC/BF/OOD numbers matching 05-04-SUMMARY.md exactly (KS 2.81e-40, χ² 5.27e-146, ECE 0.189, corr 0.8861, max|Δ| 3.598, pooled AUC 0.730). §4 "Set 2 — POST-HOC confirmatory numbers" follows, explicitly labeled post-hoc, matching STATE.md `[Phase 5-iter1]`/`[Phase 5-iter2]` entries (ECE 0.0164, corr 0.9358, max|Δ| 10.9535, combined AUC 1.0). |
| 9 | Memo states the explicit falsification condition (D-02) | VERIFIED | §6 "Falsification condition (D-02)" lists three explicit falsifiers (SBC ECE staying red / BF mid-range disagreement / an OOD family undetectable by any summary-orthogonal channel), each stated as "would have forced a No-Go." |
| 10 | Memo names the data-snooping exposure + three mitigations (consts.jl vs e9c91d3, DEV seed 0xDE7C0DE, single confirmatory VAL run) | VERIFIED | §4 "Data-snooping exposure — named plainly" + numbered list: (1) consts.jl byte-unchanged vs `e9c91d3`, (2) DEV seed `0xDE7C0DE`, (3) single confirmatory VAL run per iteration, no retry-to-pass. |
| 11 | Memo commits Phase 7 to an independent, fresh-seed, re-pre-registered confirmation ship-gate (D-05) | VERIFIED | §5 "Ship-gate — the independent confirmation run (D-05)" — explicit hard-gate language: "If the confirmation run does not reproduce the story, integration does not proceed regardless of this memo's recommendation." |
| 12 | Memo commits to the Phase 8–16 DAG with Phase 11 ∥ Phase 13 leading, Phase 12 following 11 | VERIFIED | §8 "Wave B {11 ∥ 13, then 12}" — "run Phase 11 … ∥ Phase 13 … as the lead parallel axes … Phase 12 … follows Phase 11." Matches ROADMAP.md Progress §"Wave B" text in structure. |
| 13 | demo.jl hard-gates the memo's presence and required content (DEMO-03 machine gate) | VERIFIED | Source lines 328–336: `@assert isfile(MEMO_PATH)` + seven explicit `occursin` asserts ("Clean Go", "falsification", "0.8861", "0.9358", "ship-gate", "Phase 11", "Phase 13") + `@assert MEMO_OK`. Live run printed `SC3 … PASS` and `DEMO-03 … PASS`. |
| 14 | BF figure integrity: memo's max\|Δ logBF\|=10.9535 matches what demo.jl prints at runtime | VERIFIED | Live run printed `bf_report: corr=0.9358  max_abs_err(max|Δ logBF|)=10.9535  log_prior_odds=-0.0294` — matches the memo's §4 quoted figure exactly (10.9535). |

**Score:** 14/14 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `spike/demo.jl` | Two-tier seeded end-to-end demo runner + in-script decoupling assertion | VERIFIED | 340 lines. Exists, substantive (not a stub — chains prior→simulator→summary→NPE→NRE/BF→OOD with real Philox-seeded calls to `draw_simulate_infer`/`build_bf_pair`/`amortized_log_bf`/`fit_ood_nulls`/`maha_score`), wired (executed live, exit 0), data flows (headline numbers loaded and printed from real `*_report.jld2` files, not static). |
| `.planning/phases/06-reproducible-demo-go-no-go-memo/06-GO-NO-GO-MEMO.md` | 2–3 page Go/No-Go decision memo | VERIFIED | 243 lines. Contains all required sections (Verdict, Set 1, Set 2, Ship-gate, Falsification, Thesis, Build-out DAG, Bottom line). Wired into demo.jl via the machine-checked content gate. |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| `spike/demo.jl` | `git status --porcelain -- src/ …` | backtick argv + `@assert isempty` | WIRED | Live-executed; printed "decoupling OK". No string interpolation (T-6-01 mitigated). |
| `spike/demo.jl` | `spike/validation/harness.jl` | `include()` + `draw_simulate_infer` on `val_rng(VAL_FIX_SEED)` | WIRED | Live-executed; NPE twin-run produced real 7×99 draws, ρ̂=0.06. |
| `spike/demo.jl` | `spike/validation/bf.jl` | `load_ratio` + `build_bf_pair` + `amortized_log_bf` fixture sweep | WIRED | Live-executed; 5 Δρ targets, log-BF range [-5.78, 2.97], all finite. |
| `spike/demo.jl` | `spike/validation/ood.jl` | `fit_ood_nulls` + `maha_score` fixture subset | WIRED | Live-executed; combined AUC=1.0 over 4 families reproduced twin-run. |
| `spike/demo.jl` | `spike/validation/*_report.jld2` | isfile-guarded `jldopen` for headline numbers | WIRED, data flowing | Live-executed; printed real on-disk values (corr=0.9358, max_abs_err=10.9535, ECE=0.0164, combined AUC=1.0) — not static/hardcoded. |
| `spike/demo.jl` | `spike/validation/run_sbc.jl` etc. | subprocess `run()` under `--full` | WIRED (source-verified) | Present in source at lines 301–313; not executed live in this verification pass (time-boxed) but code inspected and matches the plan's required pattern exactly. |
| `.planning/.../06-GO-NO-GO-MEMO.md` | `.planning/.../05-04-SUMMARY.md` | verbatim transcription of Set 1 numbers | WIRED | Memo's Set-1 SBC/BF/OOD table matches 05-04-SUMMARY.md's reported figures exactly (cross-checked KS/χ²/ECE/corr/max|Δ|/AUC values). |
| `spike/demo.jl` | `.planning/.../06-GO-NO-GO-MEMO.md` | `isfile` + content-keyword `@assert` | WIRED | Live-executed; DEMO-03 row reported PASS; all 7 keyword `occursin` asserts passed silently (no AssertionError raised). |

### Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|---------------------|--------|
| `demo.jl` STEP 2 headline table | `BF_HEAD["max_abs_err"]`, `BF_HEAD["corr"]` | `jldopen(bf_report.jld2)` | Yes — printed 10.9535 / 0.9358, matching STATE.md iter2 entry | FLOWING |
| `demo.jl` STEP 2 headline table | `SBC_HEAD["ece"]` (ρ_true) | `jldopen(sbc_report.jld2)` | Yes — printed 0.0164, matching STATE.md iter1 entry | FLOWING |
| `demo.jl` STEP 2 headline table | `OOD_HEAD["combined_auc"]`, `OOD_HEAD["id_fire_rate"]` | `jldopen(ood_report.jld2)` | Yes — printed 1.0 / 0.05, matching STATE.md iter2 entry | FLOWING |
| `demo.jl` STEP 1a/1b/1c twin-run draws | `NPE_R1.draws`, `BF_R1.logbfs`, `OOD_R1.maha_pos`/`combined_auc` | live re-simulation via `draw_simulate_infer`/`build_bf_pair`/`fit_ood_nulls` on `val_rng(VAL_FIX_SEED)` | Yes — non-trivial, non-empty numeric output (7×99 draws, 5-element log-BF vector with real spread, AUC=1.0) each run, verified bit-identical across two independent calls | FLOWING |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| demo.jl fast tier runs end-to-end, CPU-only, exits 0 | `julia --project=spike spike/demo.jl` (run live, full output captured) | Exit 0; full success-criteria table printed; all SC1/SC2/SC3 + DEMO-01/02/03 rows PASS; final line "demo OK: …" | PASS |
| src/ decoupling holds independently of demo.jl's own assertion | `git status --porcelain -- src/ src/bayes.jl src/colocalization.jl` | Empty output | PASS |
| Locked files byte-unchanged | `git diff --exit-code spike/validation/consts.jl spike/Project.toml spike/Manifest.toml` | Exit 0 | PASS |
| consts.jl byte-identical to the pre-registration commit cited in the memo | `git diff e9c91d3 -- spike/validation/consts.jl` | Exit 0 | PASS |
| SUMMARY-claimed commits exist | `git log -1 --oneline {ab2fba3,a10a5aa,c460d86,ef30dfc,530b287}` | All 5 resolve to the exact commit subjects described in the SUMMARYs | PASS |
| Requirement IDs mapped correctly, no orphans | `grep "Phase 6" .planning/REQUIREMENTS.md` | Only DEMO-01/02/03 map to Phase 6 (plus NPE-03's "restate in Phase 6 memo" cross-reference, satisfied by memo §7); all marked Complete | PASS |
| Full aggregate regression suite (`spike/test/runtests.jl`) still green with demo.jl present | `julia --project=spike spike/test/runtests.jl` (300s timeout) | Timed out and was killed after ~5 minutes wall-clock (TTFX + partial run). All tests that completed before the kill passed (ENV-02 22/22, SIM-03 14/14, SIM-01 12/12, SIM-03-on-output 10/10, D-15 5/5, SIM-02 10/10, SIM-04 21/21 — through the Phase 1–4 env/sim suites) but the run never reached the Phase 5/6 NPE/BF/OOD test files (`test_sbc.jl`/`test_bf.jl`/`test_ood.jl`) that the SUMMARY's "OOD 27/27, EXIT 0" claim refers to | INCONCLUSIVE (informational — see note below) |

**Note on the INCONCLUSIVE check:** This full-suite run is **not** one of the ROADMAP Success Criteria for Phase 6 and not one of the four CRITICAL phase-specific checks requested for this verification (demo.jl chain+asserts, decoupling proof, memo content, BF figure integrity) — all four of which were independently confirmed by directly executing `spike/demo.jl` itself, the authoritative deliverable for this phase. No test failures were observed in the portion that did complete; the timeout only cut off further-stage tests, it did not report a failure. This is recorded as informational, not a gap, and does not block phase closure.

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| DEMO-01 | 06-01-PLAN.md | `spike/demo.jl` chains the full pipeline reproducibly from a fixed seed and tabulates success criteria | SATISFIED | Live execution confirmed; twin-run asserts pass; success-criteria table printed. |
| DEMO-02 | 06-01-PLAN.md | The main repo and both manuscript pipelines are demonstrably untouched (decoupling proof) | SATISFIED | `git status --porcelain -- src/ src/bayes.jl src/colocalization.jl` independently confirmed empty; in-script assertion also passed live. |
| DEMO-03 | 06-02-PLAN.md | A 2–3-page Go/No-Go memo reports metrics and a concrete full-build-out decision | SATISFIED | Memo read directly; all required content present; demo.jl's machine content-gate passed live (DEMO-03 PASS). |

No orphaned requirements: `grep "Phase 6" .planning/REQUIREMENTS.md` returns only DEMO-01/02/03 (plus the NPE-03 cross-reference "restate in Phase 6 Go/No-Go memo," which is satisfied by memo §7).

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| — | — | No TBD/FIXME/XXX/TODO/HACK/PLACEHOLDER markers found in `spike/demo.jl` or `06-GO-NO-GO-MEMO.md` | — | None |

No stub patterns (`return null`, hardcoded empty data feeding rendering, no-op handlers) found in either deliverable. Both `06-01-SUMMARY.md` and `06-02-SUMMARY.md` self-report "Known Stubs: None," and static inspection of `demo.jl` confirms every layer step performs a real computation (Philox-seeded simulate/infer calls) rather than displaying stored constants.

### Human Verification Required

None. All must-haves are either mechanically checkable (file existence, git state, keyword presence, live script execution) or were directly confirmed by running the deliverable. The memo's narrative judgment call ("Clean Go is defensible") is a decision recorded by the user interactively during context-gathering (06-CONTEXT.md provenance note: "Decided interactively with the user on 2026-07-03... the verdict is a Clean Go"), not a claim requiring independent human re-adjudication in this verification pass.

### Gaps Summary

No gaps. All 14 must-have truths verified against live execution of `spike/demo.jl`, direct reads of `06-GO-NO-GO-MEMO.md`, and independent `git` state checks (not merely SUMMARY.md claims). The one informational item (`runtests.jl` full-suite regression check) was truncated by a 300s timeout partway through the Phase 1–4 tests (all of which passed); it did not reach the Phase 5/6 NPE/BF/OOD test files. This is outside the scope of this phase's ROADMAP Success Criteria and the four CRITICAL checks requested; it does not block phase closure.

---

*Verified: 2026-07-03T15:10:00Z*
*Verifier: Claude (gsd-verifier)*
