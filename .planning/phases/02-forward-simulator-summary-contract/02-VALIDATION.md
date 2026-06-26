---
phase: 2
slug: forward-simulator-summary-contract
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-06-26
---

# Phase 2 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (`@testset`), extending Phase-1 `spike/test/runtests.jl` |
| **Config file** | none — single runner (`spike/test/runtests.jl`) |
| **Quick run command** | `julia --project=spike spike/test/runtests.jl` (quick subset: small size, few sweep points) |
| **Full suite command** | `julia --project=spike spike/test/runtests.jl` (full sweep + calibration check) |
| **Estimated runtime** | quick ~seconds; full sweep + induced-μ calibration ~minutes (CPU) |

---

## Sampling Rate

- **After every task commit:** quick subset — 1 sim at 256², a ~5-point monotonicity check, contract assertion (seconds)
- **After every plan wave:** full sweep (≥15 ρ_true points) + induced-μ KS at calibration size + perturbation tests
- **Before `/gsd:verify-work`:** full suite green + CairoMakie figures regenerated
- **Max feedback latency:** quick subset < ~30s; full < ~5 min

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------------|-----------|-------------------|-------------|--------|
| 2-01-01 | 01 | 1 | (env) | Spike env extended (StatsBase/Images/ImageFiltering/CairoMakie) WITHOUT downgrading NeuralEstimators 0.2.1 / adding CUDA; root untouched vs `f581d95` | unit (guard) | reuse Phase-1 CUDA-absence + NE-pin checks in `runtests.jl` | ✅ partial | ⬜ pending |
| 2-02-01 | 02 | 2 | SIM-01 | `simulate_pair(rng,θ)` runs all 7 stages, returns 2× `Matrix{Float64}` from the 7-tuple θ | unit | `julia --project=spike spike/test/runtests.jl` | ❌ W0 | ⬜ pending |
| 2-02-02 | 02 | 2 | SIM-03 | Sim output builds a valid `MultiChannelImage`; `patch`/`correlation` yield ≤64 finite per-patch ρ (background small-positive, not 0.0; ≥15-px patch floor) | unit | same | ❌ W0 | ⬜ pending |
| 2-03-01 | 03 | 3 | SIM-02 | Induced-μ KS/Wasserstein distance to `Truncated(Cauchy(0,0.3),-1,1)` < tol; fitted ĝ monotone | integration | same (full) | ❌ W0 | ⬜ pending |
| 2-04-01 | 04 | 4 | SIM-04 | Spearman(ρ_true, mean patch-corr) ≥ threshold over the sweep | integration | same (full) | ❌ W0 | ⬜ pending |
| 2-04-02 | 04 | 4 | SIM-04 | Spillover & sub-pixel shift perturb mean patch-corr beyond MC noise (paired) | integration | same (full) | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] Extend `spike/Project.toml`: add `StatsBase` (REQUIRED before `include(src/colocalization.jl)` — calls `corspearman`/`corkendall`), `Images`/`ImageFiltering`, `CairoMakie` (+ optional `HypothesisTests`); re-resolve; re-freeze `spike/Manifest.toml`; re-run Phase-1 smoke as the regression guard
- [ ] `spike/test/test_simulator.jl` — SIM-01..04 assertions, `include`d from `runtests.jl`
- [ ] `spike/simulator/{prior,forward,calibration}.jl` + `spike/contract.jl` — the units under test
- [ ] Seed strategy: explicit `rng = Random.Xoshiro(seed)` threaded through (D-14) so the gate is deterministic and CPU-only

*Existing infra reused: the Phase-1 `runtests.jl` testset harness + the CUDA-absence / NeuralEstimators-pin guards.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Plausibility figures look right | SIM-04 | Visual confirmation complements the quantitative gate | Open the CairoMakie-generated `spike/figures/` plots: image panels at low/high ρ_true; the ρ_true→mean-patch-corr sweep curve; spillover/shift perturbation panels |
| Induced-μ overlay vs real data | SIM-02 (optional) | Sanity-only, not a gate | Optionally overlay induced-μ on the `test/test_images/{positive,negative}` patch-correlation μ |

---

## Validation Sign-Off

- [ ] All tasks have automated verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Quick feedback latency < ~30s
- [ ] Root untouched vs baseline `f581d95` asserted (decoupling preserved)
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
