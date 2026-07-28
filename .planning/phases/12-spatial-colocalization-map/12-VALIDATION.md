---
phase: 12
slug: spatial-colocalization-map
status: planned
nyquist_compliant: true
wave_0_complete: false
created: 2026-07-27
amended: 2026-07-28
---

# Phase 12 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Derived from `12-RESEARCH.md` § Validation Architecture. Amended premises are recorded in
> `12-CONTEXT.md` § AMENDMENT 2026-07-27 — read that block before planning against this file.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (`@testset`), invoked directly — no `Pkg` test target for the spike |
| **Config file** | none — `spike/test/runtests.jl` is the single entry point |
| **Quick run command** | `julia --project=spike spike/test/runtests.jl` |
| **Full suite command** | `julia --project=spike -t auto spike/test/runtests.jl` |
| **Decoupling proof** | `julia --project=. -e 'using Pkg; Pkg.test()'` — must stay green AND show `src/` untouched |
| **Reported-run commands** | `julia --project=spike -t auto spike/validation/run_p12_*.jl` (each writes one `.jld2`) |
| **Estimated runtime** | unit + smoke tiers ~60-120 s; reported runners are minutes-to-hours and are **never** invoked from the suite |

### CRITICAL ordering constraint

`spike/test/runtests.jl` ends with `include("test_p13_correction.jl")`, whose outer `@testset`
**throws** on two committed measured misses. That file's own comment records the consequence: *"a
thrown testset aborts the remaining includes, so any sibling placed after it would silently never
run."*

**Every Phase-12 `include` must be inserted BEFORE that line.** A Phase-12 test placed after it would
report green by never executing — a silent-pass failure mode, not a cosmetic ordering nit.

---

## Sampling Rate

**AMENDED 2026-07-28 after the plan-checker gate.** The original per-commit tier was
`julia --project=spike spike/test/runtests.jl`, and that is wrong twice over:

1. **It cannot signal a pass.** `runtests.jl` ends with `test_p13_correction.jl`, whose outer
   `@testset` **throws** on two committed measured misses, so the suite already exits non-zero by
   design. Any per-commit check whose signal is that process's exit status is unusable, which is why
   eight plan verifies had been written as `runtests.jl | grep <name>` — and a pipeline's status is
   `grep`'s, while `grep` matches a *failing* testset's name just as happily as a passing one.
2. **It cannot meet the stated latency.** Summing the plans' own budgets — prior ≤60 s, architecture
   ≤60 s with a 1-epoch flow build, datagen ≤90 s, train ≤120 s, coverage ≤120 s, SBC ≤120 s, plus
   40,000 Cholesky-based lattice draws, five `git` shell-outs in the decoupling tests, and the whole
   pre-existing Phase-1-5/11/13 suite ahead of them — the realistic figure is 8-10 minutes.

Amended tiers:

- **After every task commit:** run the ONE test file the task touched, directly, so a failing testset
  throws and the process exits non-zero:
  `julia --project=spike -e 'using Test; include("spike/test/test_p12_<X>.jl")'`. This is falsifiable
  and lands in 5-30 s.
- **After every plan wave:** the full spike suite `julia --project=spike -t auto spike/test/runtests.jl`
  **plus** `julia --project=. -e 'using Pkg; Pkg.test()'` (the decoupling proof — every wave, not once
  at the end). At this tier the operator reads the Phase-12 testset results directly; the suite's own
  exit status is expected non-zero because of the Phase-13 correction arm.
- **Before `/gsd:verify-work`:** full spike suite with every Phase-12 testset green + every reported
  `.jld2` artifact present and re-derivable + the D-12 Stage-2 verdict recorded.
- **Max feedback latency:** ~30 s for the per-commit tier; ~10 min for the per-wave tier.

---

## Per-Task Verification Map

No REQ-IDs are mapped to Phase 12 (`phase_req_ids` is null). The map is keyed on the **amended**
success criteria and the CONTEXT decisions instead. Task IDs are assigned by the planner; the
Behaviour column is the stable key.

| Behaviour | Source | Wave | Test Type | Automated Command | File Exists | Status |
|---|---|---|---|---|---|---|
| Tier-1 constants literal; seeds provably disjoint (incl. **Phase-13** seeds) | D-12, carried-forward | 0 | unit | `test_p12_consts.jl` | ❌ W0 | ⬜ pending |
| CAR/GP kernels SPD; per-cell marginal sd = 1 after the rescale step | D-03, D-05 | 0 | unit | `test_p12_lattice.jl` | ❌ W0 | ⬜ pending |
| DCT-II basis orthonormal; `idct(dct(F)) == F` to float tolerance | D-04, D-07 | 0 | unit | `test_p12_lattice.jl` | ❌ W0 | ⬜ pending |
| Prior parametrized by induced lag-1 correlation; r₁ monotone in the parameter | Research Pattern 1 | 0 | unit | `test_p12_lattice.jl` | ❌ W0 | ⬜ pending |
| **SIM-02 per region**: max over 64 regions of W1(induced μ, `MU_PRIOR`) ≤ pre-declared tol | D-05 | 0 | integration | `run_p12_sim02.jl` → asserted by `test_p12_prior.jl` | ❌ W0 | ⬜ pending |
| Per-region `ghat` atom mass measured and recorded (**reported, not gated**) | D-05, S-3 | 0 | integration | same artifact | ❌ W0 | ⬜ pending |
| `reshape_summary` round-trips `encode_d01` exactly (no transpose) | D-02 | 0 | unit | `test_p12_architecture.jl` | ❌ W0 | ⬜ pending |
| Mask rows never z-scored; mask is channel 2 | D-02 | 0 | unit | `test_p12_architecture.jl` | ❌ W0 | ⬜ pending |
| CNN summary net + wide `NormalisingFlow` builds, trains 1 epoch, `sampleposterior` returns D×N | D-01, D-02 | 0 | smoke | `test_p12_architecture.jl` | ❌ W0 | ⬜ pending |
| Field at constant ρ reproduces the current simulator **bit-for-bit** | D-06 | 0 | unit | `test_p12_prior.jl` (golden fixture captured **before** the first simulator edit) | ❌ W0 | ⬜ pending |
| Separable upsample equals a reference interpolant to tolerance | D-06 | 0 | unit | `test_p12_prior.jl` | ❌ W0 | ⬜ pending |
| Decoupling: `src/`, `spike/Project.toml`, `spike/Manifest.toml` byte-unchanged; `corpus/` untouched | CLAUDE.md, D-11 | 0 | integration | `test_p12_decoupling.jl` | ❌ W0 | ⬜ pending |
| **D-12 Stage-1**: leave-region-out ridge beats the per-region prior baseline, with positive control | D-12 | 1 | reported | `run_p12_stage1_ridge.jl` | ❌ W1 | ⬜ pending |
| Correlation-length identifiability ridge (raw rows + Moran's I) vs prior-mean baseline | D-08, S-1 | 1 | reported | `run_p12_ell_ridge.jl` | ❌ W1 | ⬜ pending |
| `chromatic_eps` identifiability ridge on the Phase-11 pool | Research Pitfall 1, OQ6 | 1 | reported | `run_p12_eps_ridge.jl` | ❌ W1 | ⬜ pending |
| SBC: global term + deviation-energy summary; **randomized ranks**; equivalence test on nuisances | D-01, D-07 | 3 | reported | `run_p12_sbc.jl` | ❌ W3 | ⬜ pending |
| **Leave-region-out predictive coverage** on simulation, spatial vs matched ablation — coverage **and** a proper scoring rule | D-09, D-10 | 3 | reported | `run_p12_coverage.jl` | ❌ W3 | ⬜ pending |
| Same on the 6 real TIFFs; per-region n reported; sealed holdout untouched | D-11 AMENDED | 3 | reported | `run_p12_coverage.jl --real` | ❌ W3 | ⬜ pending |
| **Radial-energy guard**: fraction of deviation-field energy that is radial | S-4 | 3 | reported | `run_p12_guards.jl` | ❌ W3 | ⬜ pending |
| **ε = 0 ablation**: spatial advantage survives with chromatic held at 0 | S-4 | 3 | reported | `run_p12_guards.jl` | ❌ W3 | ⬜ pending |
| **Offset-grid guard**: no degradation when the field lattice is offset half a cell | D-06 | 3 | reported | `run_p12_guards.jl` | ❌ W3 | ⬜ pending |
| Root suite green and `src/` provably untouched | CLAUDE.md | every wave | integration | `julia --project=. -e 'using Pkg; Pkg.test()'` | ✅ exists | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

### Reported vs gated — a distinction this phase must keep

`reported` runners write a `.jld2` and print numbers. They are **not** run by `runtests.jl` and they
do not turn the suite red. Only the pre-registered D-12 Stage-1 and Stage-2 thresholds gate. Every
other reported number is evidence, and several (atom mass, ε identifiability, vacuous-column
shrinkage) are explicitly **reporting-only** — recording them as gates would manufacture the
over-powered-χ² failure mode this project already hit once in Phase 7.

---

## Wave 0 Requirements

- [ ] `spike/validation/p12_consts.jl` — Tier-1 pre-registration: fresh DEV seed, reserved counters,
      Stage-1/Stage-2 thresholds, `P12_ITERATION_ALLOWANCE = 1`, per-region W1 tolerance, jitter,
      radial-energy ceiling, offset-grid tolerance, wall-clock ceiling
- [ ] `spike/test/test_p12_consts.jl` — literal re-assertion of every Tier-1 constant plus an
      **executable** seed-disjointness check that includes Phase 13's seeds
- [ ] `spike/simulator/p12_lattice.jl` + `spike/test/test_p12_lattice.jl`
- [ ] `spike/simulator/p12_prior.jl` + `spike/test/test_p12_prior.jl` — including the constant-ρ
      bit-for-bit golden fixture, **captured before the first simulator edit**
- [ ] `spike/npe/p12_architecture.jl` + `spike/test/test_p12_architecture.jl`
- [ ] `spike/test/test_p12_decoupling.jl`
- [ ] Wire every Phase-12 include into `spike/test/runtests.jl` **before** `test_p13_correction.jl`

*No framework install needed — Julia stdlib `Test` is already the harness.*

---

## Manual-Only Verifications

| Behaviour | Source | Why Manual | Instructions |
|---|---|---|---|
| D-12 **Stage-1** descope verdict | D-12 | A pre-registered judgement call on a reported number, taken **before** the training spend | Read `run_p12_stage1_ridge.jl`'s artifact, compare against the Tier-1 threshold, record PROCEED/DESCOPE in a verdict doc with the numbers quoted. Follow Phase 11's `11-PROBE-VERDICT.md` pattern. |
| D-12 **Stage-2** calibration verdict | D-12, D-13 | Same, after training; one documented iteration allowed | Compare per-region calibration against the Tier-1 gate. On failure, either spend the single `P12_ITERATION_ALLOWANCE` or descope to the D-13 ablation model. |
| SC3 amendment pre-declaration | D-09 | Must be frozen **before** any result exists, so it cannot be verified by a test that runs after | Commit the amendment doc in Wave 0/1, before `run_p12_coverage.jl` is ever executed. Mirrors Phase 13's `13-SC2-AMENDMENT.md`. |

---

## Validation Sign-Off

Signed off 2026-07-28, after the plan-checker gate and the fixes it forced.

- [x] All tasks have an `<automated>` verify or a declared Wave 0 dependency
- [x] **Every verify is falsifiable.** Eight `runtests.jl | grep` verifies were replaced with direct
      per-file includes (12-03, 12-04, 12-05, 12-07, 12-08 x2, 12-09, 12-12) and a literal `|| true`
      was removed from 12-18. 12-01's `grep -q` is retained deliberately: its criterion genuinely is
      "did this testset execute at all", which presence proves and which is the whole point of the
      include-ordering trap it guards.
- [x] Sampling continuity: no 3 consecutive tasks without a falsifiable automated verify. Before the
      fix, waves 2, 3 and 5 had runs of 4, 4 and 3.
- [x] Wave 0 covers every MISSING reference above
- [x] No watch-mode flags
- [x] Every Phase-12 include sits before `test_p13_correction.jl` in `runtests.jl` (12-01 Task 3, and
      the insertion point is located by CONTENT rather than line number, because Phase 13 is executing
      concurrently and may have added includes)
- [x] Decoupling proof runs every wave, not once
- [x] Feedback latency: ~30 s per-commit tier (see the amended Sampling Rate)
- [x] `nyquist_compliant: true` set in frontmatter

Every behaviour row above is owned by a plan: Wave-1 rows by 12-06/12-11/12-13, and the Wave-3 rows
(SBC, leave-region-out coverage on simulation and on the real TIFFs, and the three S-4/D-06 guards) by
12-18, 12-16, 12-19 and 12-20 respectively — the seven plans that did not exist when this document was
written.

**Approval:** signed off with one execution blocker outstanding — the Δρ semantics question recorded in
`12-02-PLAN.md` and `.planning/STATE.md`. That is a scope decision for the user, not a validation gap;
no wave may execute until it is answered, because 12-02 (Wave 1) freezes the amendment write-once.
