# Phase 4: NPE Training + ADVI Benchmark + Ablation - Discussion Log

> **Audit trail only.** Not consumed by downstream agents (researcher, planner, executor).
> Decisions captured in `04-CONTEXT.md` — this log preserves the discussion.

**Date:** 2026-06-30
**Phase:** 04-npe-training-advi-benchmark-ablation
**Mode:** discuss `--all` (all gray areas auto-selected, discussed interactively)
**Areas:** ADVI baseline, Δρ design, Summary net, Ablation gate, Speedup protocol (→ GPU / scaling / core-scaling)

## Round 1

### ADVI baseline (→ D-01/D-02)
- **Options:** Separate env/process → JLD2 *(chosen)* · Spike-local ADVI (add Turing) · Port src/bayes.jl in place
- **Selected:** Separate env/process → JLD2.
- **Why it mattered:** `src/bayes.jl:325` uses the removed `vi(m, ADVI(int,int))` API and the spike env excludes Turing (adding it caps NeuralEstimators <0.2.1 — the Phase-1 landmine); editing `src/` violates spike decoupling. An isolated ported baseline → JLD2 hand-off is the only path satisfying all three constraints.

### Δρ design (→ D-03)
- **Options:** Difference of two single-stack passes *(chosen)* · Paired-summary joint estimator · You decide
- **Selected:** Difference of two single-stack posterior passes (sample − control), reusing the single-stack cache unchanged.

### Summary net (→ D-05)
- **Options:** MLP summary net *(chosen)* · DeepSet now · You decide
- **Selected:** MLP (DeepSet deferred per BACK-02 until MLP sufficiency is measured).

### Ablation gate (→ D-06/D-07)
- **Options:** ρ_true/Δρ-led advisory-with-rule *(chosen)* · Strict all-parameter gate · You decide
- **Selected:** ρ_true/Δρ-led, advisory-with-rule (pre-set margin; parsimony + OOD tie-break; report all 7 θ; SBC-pass-but-high-RMSE = insufficiency).

## Round 2

### Speedup protocol (→ D-08/D-09)
- **Options:** Raw-stack→posterior, summary time counted *(chosen + additions)* · Forward-pass-only NPE clock · You decide
- **Selected:** Raw-stack→posterior, NPE clock = summary + forward (training excluded), comparable-RMSE = pre-registered tolerance.
- **User additions (drove Round 3):** also benchmark on GPU; determine Big-O scaling; determine CPU-core dependence; consider KernelAbstractions.jl for cpu/gpu auto-dispatch.

### Baseline env (→ D-01/D-02)
- **Options:** Dedicated baseline Project.toml *(chosen)* · Repair + use parent env · You decide
- **Selected:** Dedicated isolated baseline env with ported API; the parent env is itself broken so the baseline must not depend on repairing it.

### Δρ benchmark (→ D-04)
- **Options:** Pair reserved holdout stacks *(chosen)* · ρ_true headline, Δρ → Phase 5 · You decide
- **Selected:** Pair reserved holdout stacks; both NPE and ADVI on identical pairs → symmetric RMSE.

## Round 3 (scoping the speedup additions vs CLAUDE.md CPU-only constraint)

### GPU scope (→ D-10/D-11)
- **Options:** CPU gates; GPU optional add-on *(chosen)* · GPU required deliverable · You decide
- **Selected:** CPU gates all pass/fail; GPU reported only if a CUDA device is present (graceful skip). KernelAbstractions evaluated only for summary/simulator kernels and only if it doesn't perturb the pinned env.

### Scaling study (→ D-12)
- **Options:** Amortization (N datasets) + input-size, both methods *(chosen)* · NPE inference only · You decide
- **Selected:** Wall-clock vs N datasets (amortization) and input size, for both NPE and ADVI; report empirical exponents — >100× as a curve.

### Core scaling (→ D-13)
- **Options:** Sweep threads, report-only *(chosen)* · Single representative thread count · You decide
- **Selected:** Sweep 1/2/4/8 threads, report-only characterization (not gating); >100× headline stated at a fixed, reported thread count.

## Deferred Ideas Raised
- DeepSet (BACK-02), RxInfer cross-check (BACK-01), SBC/BF/OOD (Phase 5, incl. Δρ-vs-`compute_BayesFactor` in BF-02), productionization/`num_patches` (Phase 7), GPU-as-required-deliverable (out — CPU-only constraint).
