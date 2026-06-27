---
phase: 03
slug: training-data-pipeline
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-06-27
---

# Phase 03 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Derived from `03-RESEARCH.md` §"Validation Architecture". All gates are
> re-runnable under the established spike `Test` harness.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (`@testset`/`@test`) — the established spike harness (Phase-1/2 pattern) |
| **Config file** | none — `spike/test/runtests.jl` is the entry point |
| **Quick run command** | `julia --project=spike spike/test/test_data_pipeline.jl` (tiny N≈64 fixture, sub-30s) |
| **Full suite command** | `julia --project=spike spike/test/runtests.jl` (smoke + simulator + data + resolve-risk gate) |
| **Estimated runtime** | ~30 seconds (quick) / ~minutes (full, excludes the real 50k generation run) |

---

## Sampling Rate

- **After every task commit:** Run `julia --project=spike spike/test/test_data_pipeline.jl` (tiny N≈64 fixture)
- **After every plan wave:** Run `julia --project=spike spike/test/runtests.jl` (full suite incl. resolve-risk gate — NeuralEstimators stays v0.2.1 after adding JLD2/Random123)
- **Before `/gsd:verify-work`:** Full suite green + a real 50k generation run reloads + a thread-count-independence check
- **Max feedback latency:** 30 seconds (quick gate)

---

## Per-Task Verification Map

> Task IDs are assigned by the planner; rows below bind each success criterion to its
> automated gate and target wave. The planner maps these onto concrete `{NN}-{plan}-{task}` IDs.

| Criterion | Plan/Wave | Requirement | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|-----------|-----------|-------------|-----------------|-----------|-------------------|-------------|--------|
| SC-1: generator → fixed-dim 128-vector | W1 | DATA-01 | Generate N=64; `size(summary_min)==(128,64)`; rows 65:128 ∈ {0,1}; mask=0 ⟺ value row==0; all `isfinite`; θ is the 7-vector | unit | `julia --project=spike spike/test/test_data_pipeline.jl` | ❌ W0 | ⬜ pending |
| SC-2a: cache round-trip | W2 | DATA-02 | Write shards, reload, `@test theta == theta_reloaded` & summaries equal | integration | `…test_data_pipeline.jl` | ❌ W0 | ⬜ pending |
| SC-2b: resume-by-skip | W2 | DATA-02 | Delete 1 shard, re-run, only it regenerates, final dataset identical | integration | `…test_data_pipeline.jl` | ❌ W0 | ⬜ pending |
| SC-2c: hash auto-invalidation | W2 | DATA-02 | Perturb a tracked source byte (or config field) → `cache_hash` differs → loader refuses stale | integration | `…test_data_pipeline.jl` | ❌ W0 | ⬜ pending |
| SC-3a: no global-standardize path | W3 | DATA-03 | `@test !( :standardize_all in names(DataLoaderModule) )` — no global symbol exported | unit | `…test_data_pipeline.jl` | ❌ W0 | ⬜ pending |
| SC-3b: fit-on-train-only | W3 | DATA-03 | `Ztr` columns ~0 mean/unit std per row; `Zva` mean/std NOT exactly 0/1 (val not used to fit); mask rows pass through unchanged | unit | `…test_data_pipeline.jl` | ❌ W0 | ⬜ pending |
| SC-4a: k-fold disjoint+complete | W3 | DATA-03 | `union(folds)==1:N` & pairwise `∩==∅`; reproducible across two loader calls (same master_seed) | unit | `…test_data_pipeline.jl` | ❌ W0 | ⬜ pending |
| SC-4b: reserved holdout excluded (STRUCTURAL) | W3 | DATA-03 | `size(load_main_pool(d).theta,2)==N` (main pool never counts the holdout); `load_main_pool` never reads `holdout.jld2`; `length(holdout) ≥ 20`; `HOLDOUT_SALT != 0` (disjoint XOR-salted key namespace). Do NOT integer-intersect holdout vs fold indices — both are 1:N-style ranges that overlap by value. | unit | `…test_data_pipeline.jl` | ❌ W0 | ⬜ pending |
| Thread-count independence | W1 | DATA-02/D-11/D-12 | Generate same small dataset with `nthreads()==1` vs `>1`; `@test` byte-identical θ + summaries | integration | `julia -t1 … ` vs `julia -t4 …` then assert identical | ❌ W0 | ⬜ pending |
| Spike decoupling (standing) | all | DEMO-02 | `src/` + root manifests byte-identical to baseline `f581d95` | integration | `git diff --quiet f581d95 -- src` (+ manifest check) | ✅ reuse | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `spike/test/test_data_pipeline.jl` — new; covers SC-1..SC-4 + thread-repro
- [ ] Wire `test_data_pipeline.jl` into `spike/test/runtests.jl` (one `include`, mirroring `test_simulator.jl`)
- [ ] Add `JLD2`, `Random123` to `spike/Project.toml`; `Pkg.resolve()`; re-freeze `spike/Manifest.toml`; extend the existing resolve-risk gate to confirm NeuralEstimators stays pinned at v0.2.1
- [ ] A tiny committed cache fixture (a few samples) for round-trip/invalidation tests; bulk cache gitignored

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Full 50k generation run | DATA-02 | Wall-clock + disk cost too large for CI fixture | Run generator at N=50k on a dev machine; confirm shards write, reload equals source, no OOM |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
