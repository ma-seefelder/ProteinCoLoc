---
phase: 1
slug: environment-smoke-gate
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-06-26
---

# Phase 1 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (ships with Julia 1.12.6) |
| **Config file** | none — `spike/test/runtests.jl` (created Phase 1, D-06) |
| **Quick run command** | `julia --project=spike spike/test/runtests.jl` |
| **Full suite command** | `julia --project=spike spike/test/runtests.jl` (single gate in Phase 1) |
| **Estimated runtime** | ~seconds to low minutes on CPU (1-param Gaussian NPE) |

---

## Sampling Rate

- **After every task commit:** Run `julia --project=spike spike/test/runtests.jl`
- **After every plan wave:** Run `julia --project=spike spike/test/runtests.jl` (same single gate)
- **Before `/gsd:verify-work`:** Harness green + clean `git status` on root files vs agreed baseline + committed `spike/Manifest.toml`
- **Max feedback latency:** ~120 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 1-01-01 | 01 | 1 | ENV-01 | — | Root `Project.toml`/`Manifest.toml`/`src/` untouched vs agreed baseline | integration (git) | `git status --porcelain src/ Project.toml Manifest.toml` (clean) | ❌ W0 | ⬜ pending |
| 1-02-01 | 02 | 2 | ENV-02 | — | NPE trains + samples + recovers θ within tolerance | smoke/unit | `julia --project=spike spike/test/runtests.jl` | ❌ W0 | ⬜ pending |
| 1-02-02 | 02 | 2 | ENV-02 / D-04 | — | No forced CUDA on CPU path (`use_gpu=false`, CUDA not loaded) | unit (guard `@test`) | `julia --project=spike spike/test/runtests.jl` | ❌ W0 | ⬜ pending |
| 1-03-01 | 03 | 2 | ENV-03 | — | `spike/Manifest.toml` pinned + committed; Julia version recorded | artifact/manual | inspect committed `spike/Manifest.toml` + `spike/.julia-version` | ❌ W0 | ⬜ pending |
| 1-04-01 | 04 | 2 | ENV-04 | — | Stack-decision note records NeuralEstimators default / BayesFlow fallback | artifact | inspect `spike/NOTES.md` (or `STACK-DECISION.md`) | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `spike/00_smoke.jl` — readable NPE smoke script (ENV-02)
- [ ] `spike/test/runtests.jl` — re-runnable hard gate wrapping the smoke + CUDA-absence guard (ENV-02 + D-04/D-06)
- [ ] `spike/Project.toml` + resolved/committed `spike/Manifest.toml` — pinned env (ENV-03)
- [ ] `spike/.julia-version` + `juliaup override set 1.12.6` + NOTES.md record — Julia pin (ENV-03)
- [ ] `spike/NOTES.md` (or `STACK-DECISION.md`) — stack decision (ENV-04)
- [ ] Root-baseline reconciliation step — precondition for ENV-01's untouched proof

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| `spike/Manifest.toml` pins exact NeuralEstimators/Flux + transitive versions | ENV-03 | Reproducibility artifact — content review, not a runtime assertion | After first `Pkg.resolve()`, open `spike/Manifest.toml`, confirm `NeuralEstimators` + `Flux` entries have exact pinned versions; confirm it is git-committed |
| Julia version recorded | ENV-03 | juliaup override is environment state, not test-assertable in-process | Confirm `juliaup override set 1.12.6` applied for `spike/` dir + `spike/.julia-version` documents `1.12.6` |
| Stack-decision note content | ENV-04 | Prose decision record | Confirm `spike/NOTES.md` states NeuralEstimators.jl default, BayesFlow (PythonCall) fallback-only |

---

## Validation Sign-Off

- [ ] All tasks have automated verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 120s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
