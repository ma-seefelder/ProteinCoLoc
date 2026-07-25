---
phase: 11
slug: registration-and-chromatic-uncertainty-as-latent
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-07-25
---

# Phase 11 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Derived from `11-RESEARCH.md` § Validation Architecture (lines 1948-2019).

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (`@testset` / `@test`) |
| **Config file** | `Project.toml:42-46` (`[extras] Test`, `[targets] test = ["Test"]`). The spike has no separate test target — files are `include`d / run directly. |
| **Quick run command (spike)** | `julia --project=spike spike/test/<file>.jl` |
| **Quick run command (src)** | `julia --project -e 'using Pkg; Pkg.test()'` |
| **Full suite command** | `julia --project -t auto -e 'using Pkg; Pkg.test()'` |
| **Estimated runtime** | spike unit files < 30 s each; root suite: minutes |
| **Existing gate machinery reused** | `test/gate/sbc.jl`, `test/gate/harness.jl`, `spike/validation/{harness,sbc}.jl`, `test/gate/gate_consts_8_v2.jl` |

---

## Sampling Rate

- **After every task commit:** `julia --project=spike spike/test/<touched>.jl` for spike-only tasks;
  `julia --project -e 'using Pkg; Pkg.test()'` for any task touching `src/`.
- **After every plan wave:** full root suite **plus** every `spike/test/` file.
- **Before `/gsd:verify-work`:** full root suite green **AND** the five §F20c provenance/decoupling
  commands producing their expected (empty) output **AND** every Wave-2/3 experiment script having
  written its report artifact.
- **Max feedback latency:** < 30 s for spike unit files; root suite on `src/`-touching commits only.

**Experiments are NOT tests.** The ladder / breakdown / attenuation / real-image runs are *reported
measurements* with pre-registered pass criteria, executed once each against frozen consts. They must
NOT be wired into the test suite (a stochastic multi-minute experiment is a flaky gate). Their
*harness functions* must have fast fixture-scale unit tests, following the existing `SBC_FIX_*`
pattern (`spike/validation/consts.jl:53-56`, `test/gate/gate_consts_8_v2.jl:169-171`).

---

## Per-Task Verification Map

Task IDs are assigned by the planner; the SC → test mapping below is binding regardless of how
tasks are partitioned across plans. Every row must be claimed by at least one task.

| SC / Decision | Behaviour verified | Wave | Test type | Automated command | File exists |
|---|---|---|---|---|---|
| **SC1a** | `shift_dx`/`shift_dy` are prior-drawn, applied in stage 6, read off the NPE (already satisfied — regression only) | — | regression | `julia --project -e 'using Pkg; Pkg.test()'` (θ round-trip, `runtests.jl:541-563`) | ✅ |
| **SC1b** | `sample_prior` returns 8 fields; `theta_prior_bounds()` returns 8; `ε` reaches stage 6 | 1 | unit | root suite → new `@testset "chromatic ε (D-09)"` | ❌ W0 |
| **SC1c** | the stage-6 `tform` is a single `AffineMap` (not a `ComposedTransformation`) and `warp` is called exactly once | 1 | unit | root suite; `@test A isa CoordinateTransformations.AffineMap` | ❌ W0 |
| **SC1d** | ε = 0 regression: composed warp `==` legacy `Translation`, and full `simulate_pair` `==` a **pre-edit golden fixture** | 1 | regression, exact `==` | `julia --project=spike spike/test/test_stage6_regression.jl` | ❌ W0 (golden captured **pre-edit**) |
| **SC1e** | backward compat: `simulate_pair(rng, θ_7field)` still runs, bit-identical | 1 | unit | same | ❌ W0 |
| **SC1f** | `d_in = 129`, `D = 8` estimator trains on 200 synthetic samples without error | 1 | smoke | `julia --project=spike spike/test/test_p11_smoke.jl` | ❌ W0 |
| **SC1g** | λ conditioning is live: same `Z₁₂₈` at λ = 0.25 vs λ = 3.0 gives materially different Δρ posterior SD | 3 | **behavioural — BLOCKS the ladder** | `julia --project=spike spike/test/test_lambda_ablation.jl` | ❌ W0 |
| **SC2a** | monotone widening: permutation Spearman(λ, Δρ-width) passes and `S ≥ SC2_SPEARMAN_FLOOR` | 4 | experiment (reported) | `julia --project=spike -t auto spike/validation/run_p11_ladder.jl` | ❌ |
| **SC2b** | honest coverage: per-rung TOST equivalence to 0.90 within ±3 pp, Holm across rungs, direction `all(<=(fwer), …)` | 4 | experiment (reported) | same run | ❌ |
| **SC3a** | in-prior coverage — identical statistic to SC2b, **not a separate run** | 4 | experiment | same run | ❌ |
| **SC3b** | breakdown curve: first beyond-prior rung with `wilson_ci(...).hi < 0.90` | 4 | experiment (reported, **not gated**) | `julia --project=spike -t auto spike/validation/run_p11_breakdown.jl` | ❌ |
| **D-02** | ρ_true attenuation: paired frozen-vs-research RMSE / post_sd / coverage on one shared eval set | 4 | experiment (reported) | `julia --project=spike spike/validation/run_p11_attenuation.jl` | ❌ |
| **D-05 / F3** | per-rung shrinkage + vacuity flag for `shift_dx`, `shift_dy`, `ε`; `prior_sd` computed **within each rung** | 4 | experiment (reported) | ladder run artifact | ❌ |
| **D-17** | real-image transfer: Δρ width monotone in injected shift on the positive arm, steeper than negative | 5 | qualitative (reported) | `julia --project=spike spike/validation/run_p11_realimage.jl` | ❌ |
| **D-12** | pinned sha is the ε commit's first parent and its `simulator.jl` contains no `chromatic_eps` | 1 | shell assertion | `git rev-parse "<EPS>^"` + `git show "<SHA>":src/amortized/simulator.jl \| grep -c chromatic_eps` | ❌ |
| **D-01 / D-16** | no shipped change: `Artifacts.toml` unchanged; shipped bundle still loads; `colocalization_amortized(...)` still runs; v1.0 pipeline files byte-unchanged | 1, 6 | shell + suite | the five §F20c commands + root suite | ❌ |
| **TOST direction** | `p11_tost_pass` uses `<=` (NOT the inverted `sbc_holm_pass` semantics), verified on a hand-computed vector | 0 | unit | `julia --project=spike spike/test/test_p11_tost.jl` | ❌ W0 |
| **θ row alignment** | `merge(sample_prior(rng), forced)` preserves field order so `collect(values(θ))` stays row-aligned | 0 | unit | `julia --project=spike spike/test/test_p11_forced_theta.jl` | ❌ W0 |
| **Pre-registration** | `P11_DEV_SEED ∉ _p11_forbidden()`; salt distinctness; ladder monotonicity; `N_PER_RUNG ≥ 271` | 0 | unit | `julia --project=spike spike/test/test_p11_consts.jl` | ❌ W0 |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `spike/validation/p11_consts.jl` — **Tier-1 pre-registration** (DEV seed, rung ladders, N per
      rung, ±3 pp band, α, Holm FWER, the D-04 one-iteration allowance, the locked F5 image-size
      arm) + executable seed-disjointness `@assert`s (RESEARCH §G22)
- [ ] `spike/test/test_p11_consts.jl` — asserts `P11_DEV_SEED ∉ _p11_forbidden()`, salt distinctness,
      ladder monotonicity, `N_PER_RUNG ≥ 271`, `(SBC_L + 1) % bins == 0`
- [ ] **Pre-edit golden fixture** `spike/test/fixtures/p11_stage6_golden.jld2` — captured **before**
      any stage-6 change (SC1d). **Ordering is load-bearing: capture it or SC1d is unverifiable.**
- [ ] **Recorded pre-ε commit sha** (`git rev-parse HEAD` before any `src/`/simulator edit) — D-12
- [ ] `spike/test/test_stage6_regression.jl` — SC1c / SC1d / SC1e, exact `==`
- [ ] Root-suite `@testset "chromatic ε (D-09)"` in `test/runtests.jl` — SC1b
- [ ] **Update the ~17 literal-`7` θ assertions in `test/runtests.jl`** to
      `length(ProteinCoLoc.theta_prior_bounds())` (RESEARCH §A4)
- [ ] `spike/test/test_p11_smoke.jl` — `d_in = 129, D = 8` trains on 200 synthetic samples (SC1f);
      **must run before the ~1.2 h datagen**
- [ ] `spike/test/test_lambda_ablation.jl` — SC1g, the λ-independence tripwire
- [ ] `spike/test/test_p11_tost.jl` — Wilson CI, TOST p-value, and the Holm **direction**
- [ ] `spike/test/test_p11_forced_theta.jl` — forced-θ field-order preservation
- [ ] Holm adjuster available in `spike/` — copy `sbc_holm_adjusted` (`test/gate/sbc.jl:352-362`)
      with attribution, or `include` the gate file. **No new dependency.**

---

## Manual-Only Verifications

| Behavior | Decision | Why Manual | Test Instructions |
|---|---|---|---|
| D-17 real-image width response | D-17 | No ground-truth ρ exists on real data, so coverage cannot be tested. Only the *direction and shape* of the width response can be compared against the simulated sweep. | Run `run_p11_realimage.jl`; inspect the width-vs-injected-shift curve on the positive arm against the simulated ladder; confirm monotone and steeper than the negative arm. Reported qualitatively. |
| Probe "below resolution" verdict | D-06 | The abort branch requires human-legible judgement of the probe curves against the pre-registered floor. | Inspect probe output; if `S_probe < 0.9` or ladder span `< 0.02 Δρ_eq`, take the RESEARCH §C10(c) branch (one re-parameterised probe, then stop or proceed) and record the decision. |
| SC2/SC3 tension statement | D-13 | A prose honesty requirement, not a computable assertion. | The report must state in plain language that widening the shift prior deliberately pulls mis-registered summaries into the training distribution, eroding the density-channel detector. |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or a Wave 0 dependency
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all ❌ references above
- [ ] No watch-mode flags
- [ ] Feedback latency < 30 s for spike unit files
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
