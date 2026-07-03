---
phase: 6
slug: reproducible-demo-go-no-go-memo
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-07-03
---

# Phase 6 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (`@testset`/`@test`) + `@assert` self-checks (`02_simulator_demo.jl` idiom) |
| **Config file** | none — `spike/test/runtests.jl` is the aggregate gate; `demo.jl` self-asserts |
| **Quick run command** | `julia --project=spike spike/demo.jl` (fast tier — must run clean, seconds-to-minutes) |
| **Full suite command** | `julia --project=spike spike/demo.jl --full` (re-runs reported-scale gates) |
| **Estimated runtime** | ~30–120 s fast tier; minutes for `--full` |

---

## Sampling Rate

- **After every task commit:** Run `julia --project=spike spike/demo.jl` (fast tier green + decoupling assert passes)
- **After every plan wave:** Run `julia --project=spike spike/test/runtests.jl` (existing aggregate gate stays green — `demo.jl` must not break it) + `git diff --exit-code spike/Project.toml spike/Manifest.toml spike/validation/consts.jl` (locked files unchanged)
- **Before `/gsd:verify-work`:** `julia --project=spike spike/demo.jl` clean AND `git status --porcelain -- src/` empty AND memo present with both number sets
- **Max feedback latency:** ~120 seconds (fast tier)

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 6-01-01 | 01 | 1 | DEMO-02 | T-6-01 (cmd injection) | git shell-out uses backtick argv form, no string interpolation | in-script `@assert` | `julia --project=spike spike/demo.jl` (exit 0) | ❌ W0 | ⬜ pending |
| 6-01-02 | 01 | 1 | DEMO-01 | — | reads `src/` read-only via `include()` only | smoke + self-assert | `julia --project=spike spike/demo.jl` (exit 0) | ❌ W0 | ⬜ pending |
| 6-01-03 | 01 | 1 | DEMO-01 | — | fast-tier chain bit-reproducible on `VAL_FIX_SEED` (twin-draw equality `@assert`) | unit (in-script) | `julia --project=spike spike/demo.jl` | ❌ W0 | ⬜ pending |
| 6-01-04 | 01 | 1 | DEMO-01 (full) | T-6-03 (dep tamper) | reported gates re-run at locked consts; no installs | integration (hard gates) | `julia --project=spike spike/demo.jl --full` | ❌ W0 (run_*.jl exist ✅) | ⬜ pending |
| 6-02-01 | 02 | 2 | DEMO-03 | — | memo reports BOTH number sets, verdict, falsification condition, build-out DAG | doc review + `@assert isfile(memo)` + keyword grep | manual review + `julia --project=spike spike/demo.jl` (memo-present assert) | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `spike/demo.jl` — the deliverable (DEMO-01/02); does not yet exist
- [ ] The Go/No-Go memo (DEMO-03) — does not yet exist
- [ ] A guarded `isfile()` check in `demo.jl` for the gitignored reports (fresh-clone robustness)
- [ ] (Optional) an in-demo memo-content assertion (grep for both number sets, "falsification", "Clean Go")

*Framework install: none — `Test` is stdlib, everything else is pinned in `spike/Manifest.toml`.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Memo is 2–3 pages, prose is accurate and honest (both number sets correctly labeled, falsification condition sound, DAG matches ROADMAP) | DEMO-03 | Prose quality/accuracy is not machine-checkable beyond keyword presence | Read the memo; confirm Set 1 (pre-registered FAIL) numbers match `05-04-SUMMARY.md`, Set 2 (post-iteration) numbers match on-disk reports/STATE.md, verdict = Clean Go, D-05 ship-gate named, D-07/D-08 DAG present |

*All script behaviors (DEMO-01/02) have automated verification via `demo.jl` self-asserts.*

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 120s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
