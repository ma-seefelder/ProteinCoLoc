---
phase: 7
slug: productionization-conditional-on-go
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-07-03
---

# Phase 7 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> **The per-grid ship-gate IS the validation strategy (D-05).** Training may run on GPU;
> the ship-gate and shipped default inference are **CPU-reproducible** (D-06).

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (+ HypothesisTests, already a dep) — no new install |
| **Config file** | root `test/runtests.jl` (new); per-grid gate under `test/gate/` |
| **Quick run command** | `julia --project -e 'using Pkg; Pkg.test()'` (fixture-scale: types, registry, resolve-risk, tiny-M SBC smoke) |
| **Full suite command** | per-grid reported gate: `julia --project test/gate/run_gate.jl --grid G` (M=2000, full sweeps, `use_gpu=false`) |
| **Estimated runtime** | quick ~tens of s; per-grid reported gate minutes–hours (grid-dependent) |

---

## Sampling Rate

- **After every task commit:** Run `Pkg.test()` (fixture-scale — types, registry, resolve, tiny-M SBC smoke).
- **After every plan wave:** Run the merged grid's full `run_gate.jl --grid G` (reported scale).
- **Before `/gsd:verify-work`:** Every **shipped** grid green on its fresh-seed CPU gate; unregistered/failed grids excluded from `_SHIPPED_GRIDS`.
- **Max feedback latency:** quick suite < ~60 s.

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 7-00 | 00 | 0 | PROD-01 | T-7-02 (silent downgrade) | Root resolves NeuralEstimators 0.2.1, no downgrade to 0.1.4 | integration | `Pkg.resolve` assertion in `runtests.jl` | ❌ W0 | ⬜ pending |
| 7-00 | 00 | 0 | PROD-01 | V5 input-validation | `colocalization_amortized` → `AmortizedColocResult`; accessors dispatch on `AbstractColocResult` | unit | `Pkg.test()` (type + accessor tests) | ❌ W0 | ⬜ pending |
| 7-00 | 00 | 0 | PROD-02 | V5 | Registry returns bundle for shipped grid; clear error + `train_and_register` hint for unregistered grid | unit | `Pkg.test()` (registry tests) | ❌ W0 | ⬜ pending |
| 7-00 | 00 | 0 | PROD-01 | A7 (GPU) | `train(...; use_gpu=true)` works with CUDA present AND degrades to CPU when absent; frozen net persists CPU-resident | integration | GPU-train smoke (mirror ENV-02) | ❌ W0 | ⬜ pending |
| per-grid | 01+ | 1..n | PROD-02 | numerical robustness | Fresh re-pre-registered **SBC** rank-uniformity (KS/χ²/ECE) at `PROD_SEED[G]`, `use_gpu=false` | reported gate | `run_gate.jl --grid G --sbc` | ❌ per-grid | ⬜ pending |
| per-grid | 01+ | 1..n | PROD-02 | BF artifact | Fresh **BF** reproduction vs **non-clamped** KDE reference (no 1e-8 floor artifact) | reported gate | `run_gate.jl --grid G --bf` | ❌ per-grid | ⬜ pending |
| per-grid | 01+ | 1..n | PROD-02 | OOD | Fresh **OOD** ROC incl. re-enabled posterior-predictive channel | reported gate | `run_gate.jl --grid G --ood` | ❌ per-grid | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

**Shipped-grid scope (finalized D-04/D-05):** 4×4, 8×8, 16×16 gate to ship; 32×32 gates **conditional** (large-image caveat + user sign-off); **64×64 is dropped — not trained, not gated, not shipped.** The 8×8-sub-tile windowed local map inherits the 8×8 gate (no separate gate).

---

## Wave 0 Requirements

- [ ] `test/runtests.jl` — root test entry (resolve-risk assertion + type/registry/accessor units)
- [ ] `test/gate/run_gate.jl` + `gate_consts_G.jl` template — per-grid fresh **CPU** gate harness (`use_gpu=false`, keyed by fresh disjoint `PROD_SEED[G]`, Random123)
- [ ] Promote `spike/validation/{harness,sbc,bf,ood}.jl` → `src/validation` or `test/gate` (grid-parametrized; `use_gpu=false` on gate + shipped inference)
- [ ] GPU-train smoke (mirror ENV-02 CPU smoke): assert `train(...; use_gpu=true)` works with CUDA present AND degrades cleanly to CPU when absent; assert frozen net persists CPU-resident
- [ ] Re-enable OOD posterior-predictive channel (iter1 finite-θ̂ guard) + add non-clamped BF baseline (memo §5 hardening)
- [ ] Framework: no new install — stdlib `Test` + HypothesisTests (existing dep); `CUDA.jl` as optional weakdep for training only

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| 32×32 conditional ship / cap-the-family decision | PROD-02 | Requires user sign-off on the large-image feasibility caveat after the gate reports | Present 32×32 gate result + feasibility caveat; user decides ship-with-caveat vs exclude from `_SHIPPED_GRIDS` |
| Root co-resolution outcome escalation | PROD-01 | Wave-0 `Pkg.resolve` is a hard gate; if it fails after Turing→extension, escalate to user (candidate: GLMakie→extension too) | Run resolve spike; if downgrade/conflict persists, stop and escalate |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (resolve gate, gate harness, GPU smoke)
- [ ] No watch-mode flags
- [ ] Ship-gate is CPU-reproducible (`use_gpu=false`, fresh `PROD_SEED[G]`) regardless of training hardware
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
