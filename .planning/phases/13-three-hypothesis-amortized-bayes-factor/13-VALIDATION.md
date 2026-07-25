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

*Criteria are keyed to the amended success criteria and the locked decisions, since no REQ-IDs are
mapped to Phase 13 (ROADMAP Requirements: TBD).*

| Task ID | Plan | Wave | Criterion | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-----------|------------|-----------------|-----------|-------------------|-------------|--------|
| TBD | 13-01 | 1 | D-01/D-04 seeds: `P13_DEV_SEED ∉ _p13_forbidden()` incl. recomputed `PROD_SEED`/`PROD_SEED_V2`; consts re-include is a no-op | T-13-01 | Reserved-seed-stream consumption blocked by executable assert | unit | `julia --project=spike spike/test/runtests.jl` | ❌ W0 | ⬜ pending |
| TBD | 13-01 | 1 | **D-15 not-a-gate sentinels**: `P13_REAL_IS_GATED == false`, `P13_REAL_QUALITATIVE_ONLY == true`, and a reflective assertion that no `P13_REAL_*` binding is a pass/fail bar (two named, reasoned exemptions only) | T-13-54 | Repudiation control — a qualitative arm cannot silently become a gate | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-02 | 1 | D-12 SC2 amendment document exists **before** any result, and its §5 records SC3's three arms with only the simulator arm gating | T-13-09 | Repudiation control (pre-registration) | source assertion | file presence + `occursin` checks | ❌ W0 | ⬜ pending |
| TBD | 13-03 | 2 | D-05 labels: `three_way_label(0.8, 0.9) != EXCLUSION`; sign symmetry; τ dead-zone → `RANDOM` | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-05 | 2 | D-06 τ probe deterministic under fixed key; unpaired streams disjoint | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-04 | 2 | D-16 α-series: `alpha_segregate(x,y,M,0.0;b) == y` **bitwise**; **`count(iszero, y′) == count(iszero, y)`** (the CORRECTED invariant — `all(y′ .> 0)` is FALSE on the committed real fixtures, which carry 2–3 source zeros before any α, and must never be written); `sum(y′) ≈ sum(y)` at `rtol = 1e-10`; `maximum(y′) ≤ P13_ALPHA_MAX_VALUE_BOUND` recorded not clamped; mask α-invariant and fraction-guarded; `m̄(α)` monotone non-increasing; missing count α-invariant | T-13-52 | Semantic-tampering control — the zero SET drives `_exclude_zero`, so a changed zero set silently changes the summary | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-04 | 2 | **D-15 substrate-agnosticism**: `alpha_segregate` is typed on `AbstractMatrix{Float64}` with no `substrate ==` branch, no `load_tiff` and no `test_images` reference, so ONE code path serves the simulated and real arms and their rungs are comparable | T-13-56 | Spoofing control — a provenance branch would make the two arms different experiments | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-09 | 6 | D-02/D-03 preconditions: script errors with the Phase-11 pointer when the net is absent; asserts 8-row θ, λ range, imsize provenance | T-13-input | V5 input validation — explicit `ArgumentError` at boundary | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-06 | 3 | SC1: `ThreeWayEvidenceNet` returns a 2×n logit matrix in one call; input width `== ratio_input_dim(G) + n_cond` | T-13-input | Shape asserts at boundary | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-06 | 3 | D-11 masked loss: `head_targets` gives `(1,1,0,0)/(0,0,1,1)/(0,1,0,1)`; perturbing an inactive logit leaves the loss bit-identical | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-06 | 3 | D-09 train smoke: NeuralEstimators `train` accepts the custom estimator + custom loss; 2 epochs on a 200-sample toy reduce validation risk | T-13-deser | Narrow deserialization surface (`Flux.state` + arch metadata) | smoke | same | ❌ W0 | ⬜ pending |
| TBD | 13-07 | 4 | **D-07** corrected logits recover the analytic log BF on the closed-form Gaussian toy (`cor ≥ 0.99`, `max\|Δ\| ≤ 0.25`) **and** the uncorrected logit FAILS the same bars | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-07 | 4 | D-07 §F2: scalar-only correction **fails** under a within-class reshaped proposal (negative control) | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-08 | 5 | D-13 `_bin_calibration` green on synthetic calibrated input, red on degenerate; empty-bin MCE trap asserted | — | N/A | unit | same | ⚠ partial (`spike/test/test_sbc.jl:77-86`) | ⬜ pending |
| TBD | 13-08 | 5 | D-01 decoupling: `git diff --quiet HEAD -- src/ test/`; `spike/Project.toml`/`Manifest.toml` unchanged; `NeuralEstimators == v"0.2.1"` by UUID; no new deps | T-13-dep | Supply-chain / scope containment | integration | same (gate clause (i)) | ⚠ extend existing | ⬜ pending |
| TBD | 13-12 | 8 | **D-12 gate**: per-class AUC over {E,R,C} + full confusion matrix at pre-registered thresholds on frozen consts | — | N/A | reported script | `julia --project=spike -t auto spike/p13/run_three_way_gate.jl` | ❌ Wave 8 | ⬜ pending |
| TBD | 13-12 | 8 | D-12 reported: binary-NRE continuity at fixed λ (reported, NOT gated) | — | N/A | reported script | same runner, separate section | ❌ Wave 8 | ⬜ pending |
| TBD | 13-12 | 8 | **D-13 gate**: per-head reliability ECE ≤ `P13_ECE_GREEN`, per-head AUC reported beside it (vacuous-pass guard) | — | N/A | reported script | same runner | ❌ Wave 8 | ⬜ pending |
| TBD | 13-13 | 8 | D-15/D-16 reported: simulated α-ladder log-BF curves + crossing point `α*` | T-13-40 | Reported-not-gated: the runner contains zero `@test` lines | reported script | `julia --project=spike spike/p13/run_alpha_series.jl` | ❌ Wave 8 | ⬜ pending |
| TBD | 13-15 | 3 | **D-15 real-image ingestion** (`test/test_images/positive\|negative_c{1,2,3}.tif`, 1028×1376): `P13_REPO_ROOT`, `real_tif`, `load_real`; grid-truncation math; frozen `ghat` anchor regression (`m̄ = +0.3292 / +0.2481`); no new dependency | T-13-51 | V5 shape/dimension validation on the real fixtures | unit | `julia --project=spike spike/test/runtests.jl` | ❌ W0 | ⬜ pending |
| TBD | 13-15 | 3 | **D-15 real α-ladder**: substrate-agnostic transform; `count(iszero, y′) == count(iszero, y)`; mask-fraction + max-value guards; monotone `m̄(α)` crossing zero strictly inside the ladder | — | N/A | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-15 | 3 | **D-15 anti-snooping**: every `.jl` under `spike/p13/` (via `readdir`, not an enumeration) is free of `open_sealed_holdout` and the corpus dir name; runner self-asserts on `read(@__FILE__, String)` | T-13-05 | **V4 access control — the `sealed_holdout` seal must not be bypassed; `test/test_images/` is the sanctioned substitute** | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-15 | 3 | **D-15 read-only discipline**: `verify_real_readonly_digest()` over `git ls-files -s test/test_images/` equal before and after; no write/mv/rm touching `test_images` anywhere in `spike/p13/` | T-13-50 | Integrity of committed fixtures | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-15 | 3 | **D-15 not-a-gate**: the real-image arm asserts it is reported, never gated | T-13-54 | Repudiation control | unit | same | ❌ W0 | ⬜ pending |
| TBD | 13-16 | 8 | **D-15 reported real-image run**: λ sweep with every log-BF printed beside its OOD verdict; Phase-13-vs-shipped OOD comparison (`P13_REAL_OOD_SHIPPED_DENSITY = 433.69` vs `_THRESHOLD = 179.14`); real α-ladder through the net; naming correction | T-13-05 | Shipped bundle is comparison-reference-only, never an evidence basis (preserves D-02) | reported script | `julia --project=spike spike/p13/run_p13_realimage.jl` | ❌ Wave 8 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

**Threshold note (D-15):** every row above tagged `13-15` or `13-16` is a **mechanics** check —
ingestion shape, invariant preservation, discipline, absence of a gate. **No row asserts that any
real-image evidence value falls on one side of a bar, because no such bar exists.** That is the
intended state, not an omission.

---

## Wave 0 Requirements

- [ ] `spike/p13/consts.jl` — D-04 pre-registration (Tier-1 constants + `_p13_forbidden()` + executable
      `@assert` self-checks), committed **before** anything runs
- [ ] `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-SC2-AMENDMENT.md` — **before any
      result exists** (D-12; follow the `07-GATE-AMENDMENT.md` precedent)
- [ ] `spike/p13/preconditions.jl` — the Phase-11 hard block (RESEARCH §D2)
- [ ] `spike/p13/labels.jl` — `three_way_label`, `head_targets`, `measure_head_log_odds`
- [ ] `spike/p13/tau_probe.jl` — the D-06 probe
- [ ] `spike/p13/alpha_series.jl` — `alpha_segregate` + invariants (substrate-agnostic: `P13_ALPHA_SUBSTRATE = (:simulated, :real)`, `P13_ALPHA_ZERO_INVARIANT = :count_iszero_preserved`)
- [ ] `spike/p13/real_images.jl` — D-15 real-image ingestion (`P13_REPO_ROOT`, `real_tif`, `load_real`)
- [ ] `spike/p13/run_p13_realimage.jl` — the reported-not-gated real-image runner
- [ ] `spike/test/test_p13_real.jl` — the five D-15 testsets, `include`d from `spike/test/runtests.jl`
- [ ] `spike/p13/net.jl` — `ThreeWayEvidenceNet`, `build_three_way_net`, `masked_two_head_bce`
- [ ] `spike/p13/result.jl` — `ThreeHypothesisColocResult <: AbstractColocResult` (RESEARCH §G3)
- [ ] `spike/test/test_p13.jl` — every unit testset above, `include`d from `spike/test/runtests.jl`
- [ ] Extend `spike/test/runtests.jl` resolve-risk gate with Phase-13 clause (i)

*No framework install needed — Julia stdlib `Test` is already the harness.*

---

## Manual-Only Verifications

| Behavior | Criterion | Why Manual | Test Instructions |
|----------|-----------|------------|-------------------|
| Corpus physical anchor — **NOT used in Phase 13** | D-15 (AMENDED 2026-07-25) | Both `physical-primary` rows are `sha256 = "PENDING-FETCH"`, `bytes = 0`, `split = sealed_holdout` — reserved for the **Phase-16 blind evaluation**. Consuming the segregated anchor here would irreversibly burn Phase 16 on the very hypothesis it evaluates. | Do NOT bypass the seal; `corpus/manifest.csv` is read for **provenance only**. Record as a declared deferral to Phase 16 and add to named limits. |
| α* crossing-point interpretation (correlation-vs-localization gap) | D-16 / Specifics | The scientific reading of where intensity-correlation starts tracking disjoint localization is a judgement call, not a threshold | Inspect the α-ladder log-BF curves from `run_alpha_series.jl` (simulated) and `run_p13_realimage.jl` (real); report `α*` descriptively in the phase summary. **Do not average the two ladders** — their α = 0 semantics differ (simulated is random by construction; real measures `m̄ = +0.329 / +0.248`). |
| `test/test_images/` real-data qualitative check | D-15 (AMENDED) | **Qualitative only** — the six committed TIFFs carry **no colocalization ground-truth label**, so the check can show exclusion verdicts behave sensibly on real microscopy, not that they are correct. Inherently non-thresholded and partly interpretive (RESEARCH §J6.1) — **reported, never gated**. | Read `run_p13_realimage.jl`'s artifact: λ-sweep log-BFs beside their OOD verdicts. Report in the phase report (13-14 §8c + Limit D): the OOD-flagging finding (both fixtures ~2.4× over threshold), the positive/negative naming correction (the "negative" pair is in fact positively correlated, m̄ = +0.2481), the qualitative-only statement, and the Phase-16 deferral. |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 120s
- [ ] Reported gates ran once on byte-locked `spike/p13/consts.jl`, sha quoted in the report
- [ ] The falsified `all(y′ .> 0)` invariant appears in no plan, no test and no validation row
- [ ] No pass/fail threshold is defined for any real-image quantity
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
