---
phase: 5
slug: validation-bundle-sbc-amortized-bf-ood
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-07-02
---

# Phase 5 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> SBC / amortized-BF / OOD gates land as re-runnable tests under `spike/test/`
> with **pre-registered consts committed before the reported run** (Phase-4 pattern).

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (`@testset`), run via `Pkg.activate("spike")` |
| **Config file** | `spike/Project.toml` (own environment — decoupling constraint) |
| **Quick run command** | `julia --project=spike -e 'using Pkg; Pkg.test()'` (single testset via `spike/test/runtests.jl` filter) |
| **Full suite command** | `julia --project=spike spike/test/runtests.jl` |
| **Estimated runtime** | ~minutes (SBC M≈2000 + ratio-net eval; CPU-only) |

---

## Sampling Rate

- **After every task commit:** Run the relevant single `@testset` (SBC / BF / OOD harness unit).
- **After every plan wave:** Run `spike/test/runtests.jl` (full spike suite — must stay green).
- **Before `/gsd:verify-work`:** Full suite green AND pre-registration file committed before the reported run.
- **Max feedback latency:** harness-unit tests < 60s; full SBC/OOD reported run is a longer gated job.

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 5-01-xx | 01 | 1 | SBC-01..04 | — | N/A | unit | `julia --project=spike spike/test/test_sbc.jl` | ❌ W0 | ⬜ pending |
| 5-02-xx | 02 | 2 | BF-01, BF-02 | — | N/A | unit | `julia --project=spike spike/test/test_bf.jl` | ❌ W0 | ⬜ pending |
| 5-03-xx | 03 | 2 | OOD-01, OOD-02 | — | N/A | unit | `julia --project=spike spike/test/test_ood.jl` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*
*Exact task IDs assigned by planner; this map is the validation contract skeleton.*

---

## Wave 0 Requirements

- [ ] `spike/validation/harness.jl` — shared θ*~π→simulate→infer harness (drives SBC/BF/OOD)
- [ ] `spike/validation/preregistered.jl` (or `.toml`) — pinned consts (M, L, KS/χ² thresholds, ECE/MCE cutoffs, ID quantile, BF correlation + log-BF error tolerance, Δρ sweep range) committed BEFORE the reported run
- [ ] `spike/test/test_sbc.jl`, `spike/test/test_bf.jl`, `spike/test/test_ood.jl` — re-runnable gates wired into `spike/test/runtests.jl`

*Existing infrastructure (`spike/test/test_npe.jl`, `spike/npe/infer.jl`) is the reuse template.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Visual traffic-light / rank-histogram / ROC figures render correctly | SBC-02, OOD-01 | Figure aesthetics are human-judged | Open CairoMakie PNGs under `spike/validation/figures/`; confirm rank histograms ~uniform, ROC monotone, traffic-light color matches computed verdict |

*All numeric gates (KS/χ² p-values, ECE/MCE, AUC, BF correlation/error) have automated assertions.*

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (harness + pre-registration + test files)
- [ ] No watch-mode flags
- [ ] Feedback latency < 60s for harness units
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
