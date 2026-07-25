---
phase: 13
slug: three-hypothesis-amortized-bayes-factor
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-07-25
---

# Phase 13 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Derived from `13-RESEARCH.md` §Validation Architecture. Planner fills the per-task rows.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (`@testset` / `@test`) — no external test package |
| **Config file** | none (Julia convention); spike gate entry point is `spike/test/runtests.jl` |
| **Quick run command** | `julia --project=spike spike/test/runtests.jl` |
| **Full suite command** | `julia --project=spike spike/test/runtests.jl` **and** `julia --project=. -e 'using Pkg; Pkg.test()'` (root must stay green — D-01 requires `src/` untouched) |
| **Estimated runtime** | ~60–120 s (fixture-scale; no training, no reported streams) |

**Wiring convention:** new test files are `include`d at the bottom of `spike/test/runtests.jl` so one
command stays the gate (precedent: `test_simulator.jl`, `test_npe.jl`, `test_sbc.jl`, `test_bf.jl`,
`test_ood.jl`).

**Reported-run convention:** long/reported runs are **scripts** (`spike/validation/run_*.jl` pattern),
never tests. The suite uses fast fixtures on fixture seeds that never consume the reported stream.

---

## Sampling Rate

- **After every task commit:** `julia --project=spike spike/test/runtests.jl`
- **After every plan wave:** spike suite **plus** root `Pkg.test()` (proves `src/` untouched and the
  root Manifest still resolves)
- **Before `/gsd:verify-work`:** full suite green; reported D-12/D-13 runs executed once on the frozen
  `spike/p13/consts.jl`
- **Max feedback latency:** 120 seconds

---

## Per-Task Verification Map

*Planner populates Task IDs. Criteria below are keyed to the amended success criteria and the locked
decisions, since no REQ-IDs are mapped to Phase 13 (ROADMAP Requirements: TBD).*

| Task ID | Plan | Wave | Criterion | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-----------|------------|-----------------|-----------|-------------------|-------------|--------|
| TBD | TBD | 1 | D-01/D-04 seeds: `P13_DEV_SEED ∉ _p13_forbidden()` incl. recomputed `PROD_SEED`/`PROD_SEED_V2`; consts re-include is a no-op | T-13-seed | Reserved-seed-stream consumption blocked by executable assert | unit | `julia --project=spike spike/test/runtests.jl` | ❌ W0 | ⬜ pending |
| TBD | TBD | 1 | D-12 SC2 amendment document exists **before** any result | — | Repudiation control (pre-registration) | source assertion | file presence check | ❌ W0 | ⬜ pending |
| TBD | TBD | 1 | D-05 labels: `three_way_label(0.8, 0.9) != EXCLUSION`; sign symmetry; τ dead-zone → `RANDOM` | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | TBD | 1 | D-06 τ probe deterministic under fixed key; unpaired streams disjoint | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | TBD | 1 | D-16 α-series: `alpha_segregate(x,y,M,0.0;b) == y` **bitwise**; `all(y' .> 0)`; `sum(y') ≈ sum(y)`; `m̄(α)` monotone non-increasing | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | TBD | 1/2 | D-02/D-03 preconditions: script errors with the Phase-11 pointer when the net is absent; asserts 8-row θ, λ range, imsize provenance | T-13-input | V5 input validation — explicit `ArgumentError` at boundary | unit | same | ❌ W0 | ⬜ pending |
| TBD | TBD | 2 | SC1: `ThreeWayEvidenceNet` returns a 2×n logit matrix in one call; input width `== ratio_input_dim(G) + n_cond` | T-13-input | Shape asserts at boundary | unit | same | ❌ W0 | ⬜ pending |
| TBD | TBD | 2 | D-11 masked loss: `head_targets` gives `(1,1,0,0)/(0,0,1,1)/(0,1,0,1)`; perturbing an inactive logit leaves the loss bit-identical | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | TBD | 2 | D-09 train smoke: NeuralEstimators `train` accepts the custom estimator + custom loss; 2 epochs on a 200-sample toy reduce validation risk | T-13-deser | Narrow deserialization surface (`Flux.state` + arch metadata) | smoke | same | ❌ W0 | ⬜ pending |
| TBD | TBD | 2 | **D-07** corrected logits recover the analytic log BF on the closed-form Gaussian toy (`cor ≥ 0.99`, `max\|Δ\| ≤ 0.25`) **and** the uncorrected logit FAILS the same bars | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | TBD | 2 | D-07 §F2: scalar-only correction **fails** under a within-class reshaped proposal (negative control) | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | TBD | 2 | D-13 `_bin_calibration` green on synthetic calibrated input, red on degenerate; empty-bin MCE trap asserted | — | N/A | unit | same | ⚠ partial (`spike/test/test_sbc.jl:77-86`) | ⬜ pending |
| TBD | TBD | 2 | D-01 decoupling: `git diff --quiet HEAD -- src/`; `spike/Project.toml`/`Manifest.toml` unchanged; `NeuralEstimators == v"0.2.1"` by UUID; no new deps | T-13-dep | Supply-chain / scope containment | integration | same (gate clause (i)) | ⚠ extend existing | ⬜ pending |
| TBD | TBD | N | **D-12 gate**: per-class AUC over {E,R,C} + full confusion matrix at pre-registered thresholds on frozen consts | — | N/A | reported script | `julia --project=spike -t auto spike/p13/run_three_way_gate.jl` | ❌ Wave N | ⬜ pending |
| TBD | TBD | N | D-12 reported: binary-NRE continuity at fixed λ (reported, NOT gated) | — | N/A | reported script | same runner, separate section | ❌ Wave N | ⬜ pending |
| TBD | TBD | N | **D-13 gate**: per-head reliability ECE ≤ `P13_ECE_GREEN`, per-head AUC reported beside it (vacuous-pass guard) | — | N/A | reported script | same runner | ❌ Wave N | ⬜ pending |
| TBD | TBD | N | D-15/D-16 reported: α-ladder log-BF curves + crossing point `α*` | — | N/A | reported script | `julia --project=spike spike/p13/run_alpha_series.jl` | ❌ Wave N | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `spike/p13/consts.jl` — D-04 pre-registration (Tier-1 constants + `_p13_forbidden()` + executable
      `@assert` self-checks), committed **before** anything runs
- [ ] `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-SC2-AMENDMENT.md` — **before any
      result exists** (D-12; follow the `07-GATE-AMENDMENT.md` precedent)
- [ ] `spike/p13/preconditions.jl` — the Phase-11 hard block (RESEARCH §D2)
- [ ] `spike/p13/labels.jl` — `three_way_label`, `head_targets`, `measure_head_log_odds`
- [ ] `spike/p13/tau_probe.jl` — the D-06 probe
- [ ] `spike/p13/alpha_series.jl` — `alpha_segregate` + invariants
- [ ] `spike/p13/net.jl` — `ThreeWayEvidenceNet`, `build_three_way_net`, `masked_two_head_bce`
- [ ] `spike/p13/result.jl` — `ThreeHypothesisColocResult <: AbstractColocResult` (RESEARCH §G3)
- [ ] `spike/test/test_p13.jl` — every unit testset above, `include`d from `spike/test/runtests.jl`
- [ ] Extend `spike/test/runtests.jl` resolve-risk gate with Phase-13 clause (i)

*No framework install needed — Julia stdlib `Test` is already the harness.*

---

## Manual-Only Verifications

| Behavior | Criterion | Why Manual | Test Instructions |
|----------|-----------|------------|-------------------|
| Physical segregation anchor (n=1) qualitative check | D-15 | `sha256 = "PENDING-FETCH"`, `corpus/data/` empty, and `split = sealed_holdout` behind the anti-snooping seal — **not readable in Phase 13** | Do NOT bypass the seal. Record as a declared deferral (RESEARCH §J2 recommends Phase 16) and add to named limits. |
| α* crossing-point interpretation (correlation-vs-localization gap) | D-16 / Specifics | The scientific reading of where intensity-correlation starts tracking disjoint localization is a judgement call, not a threshold | Inspect the α-ladder log-BF curves from `run_alpha_series.jl`; report `α*` descriptively in the phase summary. |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 120s
- [ ] Reported gates ran once on byte-locked `spike/p13/consts.jl`, sha quoted in the report
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
