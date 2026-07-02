---
phase: 9
slug: cross-method-comparator-harness
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-07-02
---

# Phase 9 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (`@testset`/`@test`) |
| **Config file** | none — single entry `spike/test/runtests.jl` |
| **Quick run command** | `julia --project=spike spike/test/runtests.jl` |
| **Full suite command** | `julia --project=spike spike/test/runtests.jl` (all phase testsets) |
| **Estimated runtime** | ~60-180 seconds (small Costes N in tests) |

New testset `spike/test/test_comparator.jl` is `include`d from `runtests.jl` (append after `test_npe.jl`), mirroring the Phase-2/3/4 scaffold-then-fill idiom. Wave 0 lands named `@test_skip` placeholders for every SC; later waves replace them with real gates.

---

## Sampling Rate

- **After every task commit:** Run `julia --project=spike spike/test/runtests.jl` (fast: small fixture, small N for Costes in tests).
- **After every plan wave:** Full suite must be green.
- **Before `/gsd:verify-work`:** Full suite green + a real `run_comparator` produces the CSV/JLD2 artifact with a stable content hash.
- **Max feedback latency:** ~180 seconds

---

## Per-Task Verification Map

| Req ID | Behavior | Test Type | Automated Command | File Exists | Status |
|--------|----------|-----------|-------------------|-------------|--------|
| CMP-01 | Costes_p, M1, M2, Pearson, Spearman all **finite** on a non-degenerate fixture (D-13b) | unit | `julia --project=spike spike/test/runtests.jl` | ❌ W0 | ⬜ pending |
| CMP-02 | Costes-p uses block scramble; p bit-reproducible under fixed seed; p∈[1/(n+1),1] | unit | same | ❌ W0 | ⬜ pending |
| CMP-02 | `manders(mci)` **equals** M1/M2 that `encode_aug` computes for same mci (D-04 consistency, no encode.jl edit) | unit/consistency | same | ❌ W0 | ⬜ pending |
| CMP-03 | Seeded θ grid → `Vector{MultiChannelImage}` reproducible; regime labels attached; harness accepts external vector | unit | same | ❌ W0 | ⬜ pending |
| CMP-04 | Traffic-light bands use pre-declared consts; `traffic_light` monotone; red only above `DIVERGENCE_FAIL` | unit | same | ❌ W0 | ⬜ pending |
| CMP-05/08 | **Determinism:** two seeded `run_comparator` runs → byte-identical table (JLD2) + identical content-hash dir (D-13a); thread-count independence | integration | same | ❌ W0 | ⬜ pending |
| CMP-06 | Tapqir anchor: reproduces published value ± tol **OR** `status=:skipped` cleanly when env absent (D-13c) | integration | same | ❌ W0 | ⬜ pending |
| CMP-07 | Audit summary emits band counts + KS-uniformity in BayesInteractomics report style | unit | same | ❌ W0 | ⬜ pending |
| — | **Resolve-risk gate:** NeuralEstimators stays v0.2.1 after DataFrames/CSV promotion; Turing/PythonCall absent from main env | env | same (extend existing block) | ⚠️ extend `runtests.jl` | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `spike/test/test_comparator.jl` — SC1/SC2/SC3 skipped scaffolds + pre-registered `const`s (mirror `test_npe.jl` Wave 0).
- [ ] `spike/comparator/config.jl` — pre-declared thresholds (`COSTES_N_SCRAMBLE`, `COSTES_BLOCK_PX`, `DIVERGENCE_WARN/FAIL`, `MASTER_SEED`, `TAPQIR_TOL`) committed BEFORE any reported run (D-14).
- [ ] Extend `runtests.jl` resolve-risk block to re-assert NeuralEstimators v0.2.1 after DataFrames/CSV promotion and to assert PythonCall/CondaPkg absent from the main env.
- [ ] Env: promote DataFrames + CSV to direct deps; re-freeze `spike/Manifest.toml`.
- [ ] `.gitignore` the CondaPkg materialized env dir; commit isolated Tapqir sub-env `Project.toml`+`Manifest.toml`.

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Tapqir published anchor value/tolerance capture | CMP-06 | Requires running the published Tapqir tutorial once to record the ground-truth anchor + tolerance (Open Question 1) | Run the eLife-2022 Tapqir tutorial in the isolated sub-env; record recovered quantity and set `TAPQIR_TOL` before wiring the automated anchor test |

---

## Validation Sign-Off

- [ ] All tasks have automated verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 180s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
