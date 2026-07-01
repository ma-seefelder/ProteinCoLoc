---
phase: 04-npe-training-advi-benchmark-ablation
plan: 04
subsystem: baseline
tags: [julia, turing, advancedvi, variational-inference, advi, jld2, cross-env-handoff, decoupling, benchmark-baseline]

# Dependency graph
requires:
  - phase: 04-npe-training-advi-benchmark-ablation
    plan: 02
    provides: "isolated spike/baseline env (Turing 0.45.0 / AdvancedVI 0.6.2 / DynamicPPL 0.41.8), read-only @model lift build_coloc_model, guarded resimulate hook"
  - phase: 04-npe-training-advi-benchmark-ablation
    plan: 01
    provides: "resimulate_holdout keyed raw re-simulation path, NPE_MASTER_SEED=0xC0FFEE pre-registered constant"
provides:
  - "spike/baseline/run_advi.jl: modern AdvancedVI vi() port (vi(m, q_meanfield_gaussian, ITER; adtype=AutoForwardDiff()) -> VIResult), constrained-support guard, holdout materialization, paired-holdout sweep, atomic integrity-checked artifact writer"
  - "spike/baseline/advi_artifact.jld2: the ONLY cross-env hand-off (D-02) the spike benchmark reads — per-pair mu SAMPLE VECTORS, ghat ρ-space mappings, 90% intervals, vi()-only steady-state wall_clock, keyed to holdout global_index"
  - "spike/baseline/holdout/holdout.jld2: reserved 20-stack holdout fixture (10 pairs) so 04-05 joins by global_index against a stable holdout"
affects: [04-05-benchmark, 05-sbc-bf-ood]

# Tech tracking
tech-stack:
  added: [Distributions, ImageTransformations, CoordinateTransformations, Interpolations (promoted to DIRECT deps of the ISOLATED spike/baseline env only)]
  patterns: [modern-AdvancedVI-vi-port, VarNamedTuple constrained-draw extraction via @varname, warm-up-before-timed-loop for steady-state wall-clock, atomic reopen-assert-mv artifact write, self-simulating baseline (resimulate raw holdout in the baseline env)]

key-files:
  created:
    - spike/baseline/run_advi.jl
    - spike/baseline/advi_artifact.jld2
    - spike/baseline/holdout/holdout.jld2
  modified:
    - spike/baseline/Project.toml
    - spike/baseline/Manifest.toml

key-decisions:
  - "Turing 0.45 rand(VIResult, N) returns CONSTRAINED VarNamedTuple draws (applies the inverse link via result.ldf) — so no manual Bijectors inverse transform is needed; the [-1,1] assertion becomes an invariant check, satisfying Pitfall 3 / T-04-CONSTR by construction"
  - "Added the 4 simulator deps to the isolated baseline env so 04-01 resimulate_holdout loads there and ADVI runs on RAW re-simulated holdout (a hard acceptance criterion); Turing 0.45.0 / AdvancedVI 0.6.2 pins provably unchanged (only project_hash moved in the Manifest)"
  - "Warm up vi() once before the timed per-pair loop so the stored wall_clock is steady-state (excludes one-time JIT compile) — the >100x headline (04-05, D-08) must pair a compiled vi() against NPE passes, not a compile-contaminated first shot"
  - "Committed the 20-stack holdout fixture (regenerable, 52KB) so 04-05 can join the artifact by a stable global_index without re-deriving the holdout"

patterns-established:
  - "Baseline is self-simulating: run_advi.jl re-simulates each holdout stack byte-for-byte from its stored key (Pitfall 4) rather than reading cached summaries; the same shared resimulate path 04-05 uses guarantees identical raw stacks across envs"
  - "Artifact provenance meta (turing_version, adtype, family, iter, N, master_seed, bayes_git_ref, timestamp) + reopen-assert-then-mv atomic write mirrors the Phase-3 cache.jl integrity idiom (T-04-ARTIFACT)"

requirements-completed: [NPE-02, NPE-03]

# Metrics
duration: ~40min
completed: 2026-07-01
---

# Phase 4 Plan 04: Modern ADVI Baseline Port + Cross-Env Artifact Summary

**The removed pre-0.35 `vi(m, ADVI(n,iter))` ADVI call is ported to the modern `vi(m, q_meanfield_gaussian, ITER; adtype=AutoForwardDiff())` -> `VIResult` API inside the isolated baseline env, run on RAW re-simulated paired holdout stacks with a hard constrained-support guard, and serialized to an integrity-checked JLD2 artifact (per-pair μ sample vectors, ρ-space mappings, 90% intervals, steady-state vi()-only wall-clock) — the single cross-env hand-off the spike benchmark consumes for the NPE-02 RMSE and NPE-03 >100× comparison.**

## Performance
- **Duration:** ~40 min
- **Completed:** 2026-07-01
- **Tasks:** 2
- **Files:** 3 created, 2 modified

## Accomplishments
- **Ported the ADVI call to the resolved baseline API (Turing 0.45.0 / AdvancedVI 0.6.2).** Verified the actual `vi` signature from installed source — `vi(model, family, max_iter; adtype, algorithm, unconstrained, ...)` — and that `q_meanfield_gaussian` / `q_fullrank_gaussian` / `AutoForwardDiff` exist while `ADVI(n,iter)` is gone. The 100k `sample(m, Prior(), …)` prior chain is skipped entirely (INFER-3 dead cost; D-08 times `vi()` only).
- **Resolved the constrained-space question empirically (Pitfall 3 / T-04-CONSTR).** In Turing 0.45, `rand(result::VIResult, N)` returns `DynamicPPL.VarNamedTuple`s of RAW (constrained) parameter values — the inverse link is applied internally via `result.ldf`. So μ draws arrive in `[-1,1]` and no manual `Bijectors` transform is needed; `_assert_constrained!` hard-checks the invariant per pair. Confirmed across all 20 stacks: global μ extrema ⊂ [-1,1].
- **Ran ADVI on RAW re-simulated holdout pairs (Pitfall 4, D-04).** `run_all_pairs` forms consecutive (sample,control) pairs, re-simulates each stack byte-for-byte via 04-01's `resimulate_holdout`, and runs one `vi()` per pair (Pattern 3 — ADVI is inherently paired). μ is read via `@varname(μ_sample)` / `@varname(μ_control)` (the removed `DynamicPPL.syms` replacement).
- **Wrote the integrity-checked cross-env artifact (D-02, T-04-ARTIFACT).** 10 pairs, all 11 required keys, full length-10000 μ sample vectors (not scalars), ρ-space mappings via the frozen `ghat`, 90% intervals, steady-state per-pair `wall_clock` (~0.5 s), and provenance `meta`. Written via the atomic `jldsave(tmp) → reopen-assert(mu_sample,wall_clock) → mv(force=true)` idiom.

## Task Commits
1. **Task 1: Port ADVI to modern AdvancedVI + constrained-support guard** — `512473e` (feat)
2. **Task 2: Paired-holdout sweep + integrity-checked JLD2 artifact** — `86c81a1` (feat)

## Files Created/Modified
- `spike/baseline/run_advi.jl` (created) — consts (ITER=1000, N=10000, q_meanfield_gaussian, AutoForwardDiff, NPE_MASTER_SEED); `_assert_constrained!`; `advi_pair`; `ensure_holdout`; `run_all_pairs`; `write_advi_artifact`; script-mode artifact regeneration guard.
- `spike/baseline/advi_artifact.jld2` (created) — the D-02 hand-off: `schema_version, meta, pairs, mu_sample, mu_control, rho_sample, rho_control, rho_s_interval, rho_c_interval, wall_clock, delta_rho_true`.
- `spike/baseline/holdout/holdout.jld2` (created) — reserved 20-stack holdout fixture (disjoint `holdout_rng` namespace, global_index -1..-20).
- `spike/baseline/Project.toml` / `Manifest.toml` (modified) — 4 simulator deps promoted to direct (isolated env only; only `project_hash` moved — no version churn).

## Decisions Made
- **rand(VIResult) is already constrained (0.45 API):** the biggest port-correctness risk (Pitfall 3) is resolved by using the constraining `rand(result, N)` path (VarNamedTuple reconstruction applies the inverse link), not by a manual bijector; the assertion enforces the invariant.
- **Simulator deps into the isolated baseline env:** required so ADVI runs on RAW re-simulated holdout images (acceptance criterion); legitimate because the env is isolated (not the pinned spike env, not `src/`, not root) and the packages are flagship JuliaImages/JuliaStats already used by the spike.
- **Warm-up before timing:** the first `vi()` JIT-compiles the whole AdvancedVI/ForwardDiff specialization (~5 s single-shot); a throwaway warm-up on the smoke pair makes every stored `wall_clock` steady-state for a fair >100× headline.
- **Commit the holdout fixture:** stabilizes the `global_index` join for 04-05 (regenerable, tiny).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Simulator deps absent from the isolated baseline env**
- **Found during:** Task 1 (activating the 04-01 `resimulate_holdout` hook 04-02 left guarded).
- **Issue:** `resimulate_holdout` transitively loads `spike/simulator/forward.jl` + `prior.jl`, which `using` Distributions, ImageTransformations, CoordinateTransformations, Interpolations — none of which were direct deps of the minimal baseline env (04-02 D-11 deliberately excluded them, so the hook stayed `RESIMULATE_AVAILABLE[]=false`). Without them ADVI could not run on RAW re-simulated holdout images (a hard acceptance criterion).
- **Fix:** Promoted the 4 packages to direct deps of the ISOLATED `spike/baseline` env only. They were already present transitively (Manifest reported no packages added/removed — only `project_hash` changed); Turing 0.45.0 / AdvancedVI 0.6.2 / Bijectors 0.15.24 / DynamicPPL 0.41.8 verified unchanged after the add.
- **Files modified:** spike/baseline/Project.toml, spike/baseline/Manifest.toml
- **Committed in:** 512473e

**2. [Rule 2 - Missing critical functionality] No holdout cache existed to re-simulate from**
- **Found during:** Task 1 (needed a holdout to feed `resimulate_holdout`).
- **Issue:** No Phase-3 cache / `holdout.jld2` is materialized in the repo, so `load_holdout` had nothing to read; the baseline had no pairs to run ADVI over.
- **Fix:** Added `ensure_holdout` — idempotently generates the reserved holdout from the disjoint `holdout_rng` namespace (D-10) via the UNCHANGED simulator (`_write_holdout`). Committed the resulting 20-stack fixture so the artifact is regenerable/verifiable and 04-05 joins by a stable `global_index`.
- **Files modified:** spike/baseline/run_advi.jl (+ spike/baseline/holdout/holdout.jld2)
- **Committed in:** 512473e (helper), 86c81a1 (fixture)

**3. [Rule 1 - Bug] First-pair wall-clock contaminated by JIT compile**
- **Found during:** Task 2 (first artifact showed wall_clock[1]≈5.1 s vs ~0.3 s for the rest).
- **Issue:** The one-time `vi()` compile inflated the first pair's single-shot `time_ns`, which would corrupt the >100× speedup headline (04-05 reads `wall_clock`).
- **Fix:** A throwaway warm-up `vi()` on the dependency-light smoke pair before the timed loop; every stored `wall_clock` is now steady-state (~0.5 s).
- **Files modified:** spike/baseline/run_advi.jl
- **Committed in:** 86c81a1

**Total deviations:** 3 auto-fixed (2× Rule 3/2 blocking-or-missing, 1× Rule 1 measurement bug). No architectural (Rule 4) changes; scope unchanged.

## Verification Evidence
- **Task 1 verify:** `include("spike/baseline/run_advi.jl")` loads; `advi_pair` defined; a single re-simulated holdout pair runs (μ ⊂ [-1,1], wall-clock measured, ρ recovered).
- **Task 2 verify (plan command):** artifact reopens with `schema_version, mu_sample, rho_sample, wall_clock, pairs` — passed. Deep check: all 11 keys present; 10 pairs keyed to holdout global_index (-1..-20); `mu_sample[1]` length 10000 (sample vectors, not scalars); global μ extrema ⊂ [-1,1]; `meta` complete (turing_version=0.45.0, adtype, family, iter=1000, N=10000, master_seed, bayes_git_ref, timestamp); steady-state `wall_clock` ~0.5 s/pair.
- **Decoupling:** `git diff --name-only <base> HEAD` touches only `spike/baseline/`; `src/` and the spike-env `Project.toml`/`Manifest.toml` pin are byte-clean; working tree clean.

## Issues Encountered
- **AdvancedVI is 0.6.2, not the RESEARCH-anticipated ~0.7** (as 04-02 flagged). The `vi(model, family, max_iter; adtype, ...)` → `VIResult` shape holds, but the return of `rand(VIResult, N)` is a `Vector{VarNamedTuple}` (constrained), not a `d×N` matrix — params are read by `@varname` indexing, not positional `syms`. This is the version-accurate port and is recorded for any future re-resolve (baseline env currency ~7 days per RESEARCH).
- Mean-field ADVI is approximate (a smoke pair recovered ρ mean −0.41 vs true −0.64); expected for VI and orthogonal to this plan — accuracy is scored in 04-05 with the pre-registered NPE_RMSE_TOLERANCE. `q_fullrank_gaussian` remains the documented A5 sensitivity option.

## Threat Surface
- **T-04-ARTIFACT (mitigate):** artifact carries `schema_version` + provenance `meta` (turing version, seed, git ref) and is written via reopen-assert-then-`mv(force=true)`; consumed by 04-05 keyed on `global_index`. Satisfied.
- **T-04-CONSTR (mitigate):** hard `_assert_constrained!` on every μ draw; verified ⊂ [-1,1] across all pairs. Satisfied.
- No new threat surface introduced (offline, CPU-only, file I/O within `spike/baseline/`).

## Self-Check: PASSED
- Files verified present: spike/baseline/run_advi.jl, spike/baseline/advi_artifact.jld2, spike/baseline/holdout/holdout.jld2, 04-04-SUMMARY.md
- Commits verified in branch history: 512473e (Task 1), 86c81a1 (Task 2)
- Decoupling: src/ and spike-env pin byte-clean; only spike/baseline/ changed

---
*Phase: 04-npe-training-advi-benchmark-ablation*
*Completed: 2026-07-01*
