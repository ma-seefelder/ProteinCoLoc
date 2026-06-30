---
phase: 4
slug: npe-training-advi-benchmark-ablation
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-06-30
---

# Phase 4 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Derived from 04-RESEARCH.md §Validation Architecture.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (+ `BenchmarkTools` for timings) — mirrors `spike/test/test_data_pipeline.jl` |
| **Config file** | none — `spike/test/runtests.jl` includes child testset files |
| **Quick run command** | `julia --project=spike spike/test/test_npe.jl` |
| **Full suite command** | `julia --project=spike spike/test/runtests.jl` |
| **Estimated runtime** | quick ~tens of seconds (model loaded from cache, not retrained); full suite minutes (incl. resolve-risk + benchmark) |

---

## Sampling Rate

- **After every task commit:** Run `julia --project=spike spike/test/test_npe.jl` (the touched SC testset; trained model loaded from cache, not retrained)
- **After every plan wave:** Run `julia --project=spike spike/test/runtests.jl` (full suite incl. Phase 1–3 gates + resolve-risk)
- **Before `/gsd:verify-work`:** Full suite green; `advi_artifact.jld2` regenerated-or-verified; all five SC testsets green at the pre-registered constants
- **Max feedback latency:** quick run under ~60s (cached model)

---

## Per-Task Verification Map

| Req ID | Behavior | Test Type | Automated Command | File Exists |
|--------|----------|-----------|-------------------|-------------|
| NPE-01 | Trained `PosteriorEstimator` returns a 7×N posterior in one pass for ≥20 holdout stacks; Δρ = MC-diff of two passes | integration | `julia --project=spike spike/test/test_npe.jl` (SC1 testset) | ❌ W0 |
| NPE-02 | Per-parameter RMSE + interval width vs ADVI artifact, CV-reported | integration | SC2 testset (reads `advi_artifact.jld2`) | ❌ W0 |
| NPE-03 | `t_advi/t_npe > 100` at comparable RMSE, fixed thread count | benchmark | SC3 testset (`@belapsed`) | ❌ W0 |
| NPE-03 | Scaling curves over N and imsize; empirical exponents (NPE ~flat, ADVI ~linear) | characterization | SC3-scaling testset (D-12) | ❌ W0 |
| ABL-01 | `:min` vs `:aug` per-parameter RMSE under k=5 CV | integration | SC4 testset | ❌ W0 |
| ABL-02 | Decision rule emits `:min`/`:aug` per the pre-registered margin; OOD note recorded | unit | SC5 testset | ❌ W0 |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

### Pre-registered constants (declared BEFORE any reported run — mirror Phase-3 `const` fixtures)

```julia
const NPE_RMSE_TOLERANCE   = 1.2     # D-09 (A2): NPE ρ_true RMSE ≤ 1.2× ADVI ρ_true RMSE
const SPEEDUP_GATE         = 100.0   # NPE-03: t_advi/t_npe > 100 at the reported thread count
const ABL_REL_MARGIN       = 0.05    # D-07 (A3): aug wins iff RMSE_aug ≤ (1−δ)·RMSE_min on ρ_true …
const ABL_FOLD_CONSISTENCY = 4       # … AND aug improves in ≥4 of 5 folds (else keep :min, parsimony)
const BENCH_THREADS        = 1       # D-13: headline thread count (sweep {1,2,4,8} reported separately)
const NPE_MASTER_SEED      = 0xC0FFEE  # Random123 reproducibility
```

---

## Wave 0 Requirements

- [ ] `spike/test/test_npe.jl` — SC1..SC5 testsets as named `@test_skip` placeholders (mirror `test_data_pipeline.jl`), included by `runtests.jl`
- [ ] `spike/Project.toml` — add `BenchmarkTools`; re-freeze `spike/Manifest.toml`; extend resolve-risk gate to re-assert NeuralEstimators v0.2.1 with BenchmarkTools present
- [ ] `spike/baseline/Project.toml` + `Manifest.toml` — new isolated env (Turing, DataFrames, JLD2, ForwardDiff, StatsBase), pinned and committed
- [ ] Holdout raw-image reproducibility check (Open Question 1) — Wave-0 round-trip assertion (else persist raw images)
- [ ] Pre-registered constants committed before any reported run (A2/A3)

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Optional GPU timing add-on | NPE-03 (D-10) | Requires a CUDA device; skips gracefully on CPU-only CI | On a CUDA host, run the SC3 testset with GPU path enabled; confirm reported separately from the CPU gate |

*All pass/fail gates are CPU-automated; GPU is an optional, manually-run bonus that never gates.*

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency acceptable (cached-model quick run)
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
