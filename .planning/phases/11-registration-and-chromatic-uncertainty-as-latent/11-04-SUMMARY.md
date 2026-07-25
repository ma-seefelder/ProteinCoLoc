---
phase: 11-registration-and-chromatic-uncertainty-as-latent
plan: 04
subsystem: npe
tags: [julia, neuralestimators, normalising-flow, logit-transform, conditioning-input, tripwire, smoke-test]

# Dependency graph
requires:
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 01
    provides: "spike/validation/p11_consts.jl — LAMBDA_MIN/LAMBDA_MAX, P11_FIXTURE_COUNTER, p11_rng, P11_RESEARCH_NET_DEVIATIONS (the four names this file's header must match), P11_PROBE_COMPARABILITY_IMSIZE"
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 02
    provides: "chromatic_eps as the 8th spike prior field and the widened SHIFT_PRIOR — the two deviations p11_theta_prior_bounds() DERIVES rows 5/6 and 8 from"
  - phase: 07-productionization-conditional-on-go
    provides: "src/amortized/architecture.jl BoundedThetaTransform (the F2 remedy ruling Q4 authorises porting) and src/amortized/train_npe.jl's fit-side construction"
provides:
  - "spike/npe/p11_architecture.jl — P11BoundedThetaTransform (+ p11_theta_to_unbounded / p11_theta_from_unbounded / fit_p11_theta_transform / p11_theta_prior_bounds), encode_lambda / decode_lambda, an order-enforcing augment_input, and build_p11_estimator"
  - "MEASURED: a d_in = 129, D = 8 estimator trains and samples through the real NeuralEstimators v0.2.1 API on CPU (R10 retired)"
  - "spike/test/test_lambda_ablation.jl — the SC1g R4 tripwire, authored and refusing to pass vacuously"
  - "spike/test/runtests.jl clause (i) — Phase 11 adds no dependency (D-01), asserted as an explicit name set"
affects: [11-05-probe-tier2-append, 11-07-research-trainer, 11-08-ladder, 11-report]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Copy-with-attribution of a src model component the spike environment cannot depend on, carrying the whole rationale block plus a note naming which part of that rationale does NOT carry over"
    - "Runtime @assert as the enforcement of a documented ordering invariant (standardize, then vcat), paired with a unit test asserting the wrong order throws"
    - "A behavioural tripwire that exits non-zero when its input artifact is absent, so it can never report success vacuously"
    - "Tier-2 constant read defensively with an in-file placeholder that names the plan obliged to supersede it"

key-files:
  created:
    - spike/npe/p11_architecture.jl
    - spike/test/test_p11_architecture.jl
    - spike/test/test_p11_smoke.jl
    - spike/test/test_lambda_ablation.jl
  modified:
    - spike/test/runtests.jl
    - .planning/phases/11-registration-and-chromatic-uncertainty-as-latent/deferred-items.md

key-decisions:
  - "The src cross-check is a REAL test, not a skip: src/amortized/architecture.jl was MEASURED to load standalone under --project=spike (its only imports are Flux, NeuralEstimators, StatsBase), so agreement is asserted against the actual shipped implementation to 1e-12"
  - "No bare single-letter epsilon identifier appears anywhere in p11_architecture.jl (grep count 0); the clamp constant is P11_THETA_LOGIT_EPS and the atom-handling prose is explicitly scoped to row 1, with chromatic_eps named as atom-free"
  - "The F2 containment guarantee is asserted CLOSED-interval on extreme flow-space values and STRICT-interior on a deterministic grid — the plan's literal 'strictly within' is false at 50-sigma for a floating-point reason, not a model reason"
  - "P11_LAMBDA_ABLATION_FACTOR placeholder = 1.15, with the Tier-2 read already wired so plan 11-05's append needs no edit to the tripwire"
  - "spike/test/test_simulator.jl was left red: it is outside this plan's declared file set and plan 11-03 was executing in parallel; logged to deferred-items.md instead"

requirements-completed: [D-01, D-02, D-03, D-04]

# Metrics
duration: ~55min
completed: 2026-07-25
---

# Phase 11 Plan 04: Research-Net Scaffold, the F2 Port and the Two Tripwires Summary

**The Phase-11 research net's model surface exists and is proven: the F2 bounded theta transform is
ported into the spike and agrees with the real `src` implementation to 1e-12, the lambda encoder
and a runtime-enforced 129-row augmenter implement D-03, a `d_in = 129, D = 8` estimator was
MEASURED to train and sample through NeuralEstimators v0.2.1 on CPU in ~55 s (retiring R10 before
any 1.2 h datagen), and the R4 lambda-ablation tripwire is authored and exits non-zero rather than
passing vacuously.**

## Headline numbers requested by the plan's `<output>`

| Quantity | Value |
|---|---|
| **SC1f smoke-train wall time** | **54.45 s** (last run; 55.67 / 54.98 / 54.45 s across three runs) — 200 train + 50 val synthetic samples, 3 epochs, batch 32, CPU |
| **SC1f end-to-end wall time** | ~74 s including Julia start and compile (well inside the 5-minute bar) |
| **`P11_LAMBDA_ABLATION_FACTOR` placeholder** | **1.15**, held in `P11_LAMBDA_ABLATION_FACTOR_FALLBACK`; **plan 11-05 must supersede it** with the Tier-2 value derived from the D-06 probe's measured effect scale. The Tier-2 read is already wired (`isdefined(@__MODULE__, :P11_LAMBDA_ABLATION_FACTOR) ? … : …`), so the append needs no edit to the tripwire. |

## Task Commits

1. **Task 1: port the F2 bounded theta transform + define the lambda encoder** — `851d64f` (feat)
2. **Task 2: SC1f smoke (d_in = 129, D = 8 trains) + runtests clause (i)** — `b1d6228` (test)
3. **Task 3: the SC1g lambda-ablation tripwire** — `f353fda` (test)

## Accomplishments

- **Ruling Q4 is implemented and the port is provably the same transform.** `src/amortized/architecture.jl`
  was measured to load standalone under `--project=spike`, so the agreement testset compares against
  the REAL shipped implementation rather than a restatement: the logit forward map, the fitted
  `ZScoreTransform` mean/scale, `StatsBase.transform` and `StatsBase.reconstruct` all agree to
  `1e-12`, and `P11_THETA_LOGIT_EPS == SrcArch.THETA_LOGIT_EPS`. No dependency was added to make
  that work.
- **The F2 property is asserted, not assumed.** Arbitrary flow-space values (`50·randn`) map back
  inside the prior box for every one of the 8 rows; a deterministic `-4:0.5:4` grid maps back
  strictly interior; the map is strictly monotone per row (so SBC ranks are invariant to the
  change of space); `Float32` draws keep their element type.
- **The PATTERNS §11 epsilon trap is closed by construction.** `grep -cE '(^|[^_[:alnum:]])(eps|ε)([^_[:alnum:]]|$)' spike/npe/p11_architecture.jl`
  is **0**. The copied atom-handling rationale is explicitly re-scoped: it applies to ROW 1
  (`rho_true = ghat(mu*)`, a clamped map with genuine endpoint atoms) and is stated NOT to apply
  to `chromatic_eps`, which is drawn straight from `Uniform(-0.02, 0.02)` and has no atoms.
- **The standardize-then-vcat order is enforced at runtime, not merely documented.** `augment_input`
  carries `@assert size(Z128, 1) == 128`; the test asserts that 129-, 142- and 64-row inputs all
  throw, so a double-append or a "fixed" `Loader._row_partition` is caught rather than silently
  z-scoring the conditioning variable against the training pool.
- **R10 is retired for ~74 s instead of ~1.2 h.** The smoke builds the 129-row input PER SAMPLE
  (the lambda row genuinely varies, `length(unique(Z[129, :])) > 1`), draws the shift rows
  CONDITIONAL on each sample's lambda (the D-13 shape, asserted `|shift| <= lam` for all 200
  samples), fits the ported transform train-only, trains through the fixed-data `train` form, and
  reads 8 marginals back inside the prior box. A 128-row and a 130-row input both throw.
- **The R4 tripwire cannot report a false pass.** With no trained net it prints
  `research net not trained yet — this tripwire cannot pass vacuously` and exits 1. There is no
  skip path in the file (`@test_skip` count 0). It asserts a MATERIAL inequality against a
  pre-registered factor on TWO spread statistics (`std` and the sanctioned `interval_width`) across
  3 independent datasets, all of which must clear the bar, and it proves the statistic can fail via
  a deliberately dead same-lambda conditioning.
- **D-01's no-dependency claim is now executable.** `runtests.jl` clause (i) asserts the project's
  dependency NAME SET equals an explicit 16-name literal (so an addition names itself in the
  failure) and re-asserts the NeuralEstimators v0.2.1 UUID pin with the Phase-11 code on disk.

## Verification — observed results

Every command below was run and its real output observed.

| Check | Command | Observed |
|---|---|---|
| Task 1 unit gate | `julia --project=spike spike/test/test_p11_architecture.jl` | **66 pass / 0 fail**, exit **0** |
| Task 2 smoke (SC1f) | `julia --project=spike spike/test/test_p11_smoke.jl` | **26 pass / 0 fail**, exit **0**, train 54.45 s, total ~74 s |
| Task 3 authoring gate | the plan's `occursin` verify one-liner | prints `authored ok` |
| Task 3 non-vacuous refusal | `julia --project=spike spike/test/test_lambda_ablation.jl` | exit **1**, message names `spike/npe/p11_research_npe.jld2` and plan 11-07 |
| Aggregate spike suite | `julia --project=spike spike/test/runtests.jl` | exit **1** — see "Deferred Issues"; the failure is PRE-EXISTING and unrelated |
| Clause (i) itself | the `D-04 CPU-only (no CUDA dependency, none loaded)` testset in that run | **23 pass / 23**, up from 21 assertions |
| Attribution present | `grep -c 'src/amortized/architecture.jl' spike/npe/p11_architecture.jl` | `7` (bar: ≥1) |
| Port named distinctly | `grep -c 'P11BoundedThetaTransform' spike/npe/p11_architecture.jl` | `8` (bar: ≥3) |
| No bare epsilon | `grep -cE '(^\|[^_[:alnum:]])(eps\|ε)([^_[:alnum:]]\|$)' spike/npe/p11_architecture.jl` | `0` |
| File length | `wc -l spike/npe/p11_architecture.jl` | `373` (bar: ≥90) |
| 129 present in smoke | `grep -c '129' spike/test/test_p11_smoke.jl` | `13` (bar: ≥2) |
| CPU-only on every call | `grep -c 'use_gpu *= *false' spike/test/test_p11_smoke.jl` | `6` (bar: ≥1) |
| Clause (i) names its decision | `grep -c 'D-01' spike/test/runtests.jl` | `2` (bar: ≥1) |
| Phase-11 files NOT wired into the runner | `grep -c 'test_p11_smoke\|test_lambda_ablation\|test_stage6_regression' spike/test/runtests.jl` | `0` |
| Tripwire has no skip path | `grep -c '@test_skip' spike/test/test_lambda_ablation.jl` | `0` |
| Global-RNG pin present | `grep -c 'Random.seed!' spike/test/test_lambda_ablation.jl` | `1` (bar: ≥1) |
| Tier-2 read + placeholder | `grep -c 'P11_LAMBDA_ABLATION_FACTOR' spike/test/test_lambda_ablation.jl` | `6` (bar: ≥2) |
| Second spread statistic | `grep -c 'interval_width' spike/test/test_lambda_ablation.jl` | `4` (bar: ≥1) |
| Material inequality, not `!=` | `grep -c '!=' spike/test/test_lambda_ablation.jl` | `0` |
| No dependency change | `git diff --stat 2aae21f..HEAD -- spike/Project.toml spike/Manifest.toml` | empty |
| No `src/` change | `git diff --stat 2aae21f..HEAD -- src/` | empty |
| Files touched by this plan | `git diff --name-only 2aae21f..HEAD` | the 4 new spike files + `spike/test/runtests.jl` + `deferred-items.md` |
| No deletions in any commit | `git diff --diff-filter=D --name-only HEAD~1 HEAD` (×3) | empty each time |

**NOT verified here, and deliberately not claimed.** This plan trains nothing that means anything.
The 3-epoch, 200-sample smoke runs on `randn` noise standing in for a summary; it establishes API
and SHAPE conformance ONLY. **No calibration, coverage, width, monotonicity, attenuation or
identifiability claim is made or supported by anything in this plan.** In particular:

- The lambda-ablation tripwire has **never been executed against a trained net** — it has only been
  shown to refuse to run without one. Whether the conditioning is actually alive is **unknown** and
  is measured for the first time in plan 11-07.
- The `1.15` ablation factor is a **placeholder**, not a measured or pre-registered threshold.
- The `d_in = 129, D = 8` estimator is proven to TRAIN, not to CONVERGE, and not to be calibrated.
- No ship gate was run, no artifact was produced, no `src/` byte moved, and the shipped
  `amended_v2/grid_8` bundle and the Phase-7 GO are untouched (D-16).

## Decisions Made

1. **The `src` cross-check is a real assertion, not a documented skip.** The plan allowed falling
   back to a property-only test if `src/amortized/architecture.jl` could not load standalone from
   `--project=spike`. It was measured directly and it DOES load (imports: Flux, NeuralEstimators,
   StatsBase — all already spike deps), so the stronger check was used. This matters because the
   port's entire purpose is to be the SAME transform the shipped net has; agreeing with a local
   restatement would have proven nothing. An in-file comment records that if it ever stops loading,
   the correct response is to drop the cross-check, never to add a dependency (D-01).
2. **`p11_theta_prior_bounds()` is a spike-local derivation, not a call into `src`.** Its rows 5/6
   deliberately disagree with `src`'s `theta_prior_bounds()` because D-02 widened `SHIFT_PRIOR` in
   the spike only. The test asserts both the derivation (`== (minimum(SHIFT_PRIOR), maximum(SHIFT_PRIOR))`)
   and the resulting literal `(-3.0, 3.0)`, so the divergence is pinned rather than incidental.
3. **The ablation tripwire re-pins the global RNG to the SAME seed at the top of each lambda arm**,
   so both arms consume an identical flow base-sample stream and the only difference between them
   is the 129th row. This also gives the negative-direction clause a byte-exact form
   (`dead_a == dead_b`), which doubles as the reproducibility proof for the pin itself.
4. **`spike/validation/consts.jl` and `harness.jl` are deliberately NOT included by the tripwire.**
   They bind `VAL_MASTER_SEED`/`VAL_FIX_SEED` at a different integer width than `p11_consts.jl`
   does, and a `const` redeclaration at a different type is a hard error. The three
   simulator/contract includes (`forward.jl`, `contract.jl`, `encode.jl`) are taken directly, which
   is exactly what `harness.jl:55-57` does.
5. **The tripwire uses `P11_PROBE_COMPARABILITY_IMSIZE` rather than a fresh literal image size**,
   with an in-file note that it therefore makes no F5 train-joint/eval-joint claim. It is a
   behavioural check, not a reported number.
6. **Every stochastic draw in this plan rides `P11_FIXTURE_COUNTER` (99).** `p11_consts.jl` reserves
   it as FIXTURES ONLY and asserts it disjoint from all six reported counters, which is the correct
   home for a unit fixture, a smoke and a tripwire. No reported stream was consumed; `PROD_SEED_V2`,
   `VAL_MASTER_SEED`, `NPE_MASTER_SEED`, `RATIO_PAIR_SEED`, `VAL_FIX_SEED` and `CORPUS_MASTER_SEED`
   were not touched.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] The plan's "strictly within `[lo, hi]`" containment assertion is false at extreme flow-space values, for a floating-point reason**
- **Found during:** Task 1 (`bounded transform keeps mass inside the prior box`)
- **Issue:** The plan directs asserting that for `50.0 .* randn(8, 200)` every reconstructed theta
  lies **strictly** within `[lo, hi]`. Measured: it does not. `theta = lo + width·logistic(y)`; at
  `|y|` of order 40+ `logistic` underflows and the product is smaller than `eps(lo)`, so the sum
  rounds to exactly `lo`. The `chromatic_eps` row makes this easy to hit — its width is 0.04. This
  is float rounding at the boundary, not a support leak.
- **Fix:** Asserted the F2 guarantee in its true CLOSED-interval form (`lo <= x <= hi`), which is
  exactly what the `src` twin asserts for the same reason (`test/runtests.jl:563`), and added a
  SEPARATE strict-interior assertion on a deterministic `-4:0.5:4` grid — the range the flow
  actually works in — where strictness is exactly reproducible rather than tail-dependent. An
  in-file comment records the mechanism so nobody "restores" the strict form later.
- **Files modified:** `spike/test/test_p11_architecture.jl`
- **Verification:** 10/10 in that testset; `julia --project=spike spike/test/test_p11_architecture.jl` exit 0.
- **Committed in:** `851d64f`

**2. [Rule 1 - Bug] `sampleposterior` returns a VECTOR of draw matrices for K > 1 datasets, so the plan's shape assertion read K instead of D**
- **Found during:** Task 2 (`sampleposterior returns 8 marginals`)
- **Issue:** The multi-dataset call `posterior_for(est, Z[:, 1:3]; N = 8)` returned an object whose
  `size(·, 1)` is **3**, not 8 — measured, `3 == 8` failed. NeuralEstimators v0.2.1 returns a
  K-element vector of `D×N` matrices for K datasets, not a single matrix.
- **Fix:** Asserted the shape the API actually returns (`many isa AbstractVector`,
  `length(many) == 3`, `all(size(m, 1) == 8 …)`, `all(size(m, 2) == 8 …)`), with a comment naming
  the measured API shape so a later reader does not "simplify" it back to a `size(·, 1)` check.
  The single-dataset assertion (the one that matters for the D-03 read path) was already correct
  and is unchanged.
- **Files modified:** `spike/test/test_p11_smoke.jl`
- **Verification:** 26/26, exit 0.
- **Committed in:** `b1d6228`

**3. [Rule 3 - Blocking, then deliberately NOT fixed] `spike/test/runtests.jl` cannot exit 0 because of a pre-existing plan-11-02 defect**
- **Found during:** Task 2 (plan acceptance criterion "runtests.jl exits 0")
- **Issue:** `spike/test/test_simulator.jl:200` still asserts the PRE-11-02 seven-field theta
  NamedTuple. `spike/simulator/prior.jl` gained `chromatic_eps` in `ca02b0e` (plan 11-02) and this
  spike test was never updated the way its sibling `test_p11_forced_theta.jl` was (11-02 deviation
  3). The file THROWS, so `runtests.jl` aborts there and the six later includes never run at all.
- **Established as pre-existing, not caused by this plan:** `test_simulator.jl` was last modified
  in `02c3bee` (Phase 2); `git diff --name-only HEAD -- spike/test/test_simulator.jl` is empty for
  this plan; and the only edit this plan made to the suite (clause (i)) lives in the FIRST testset,
  which passes 23/23.
- **Decision:** NOT fixed. `spike/test/test_simulator.jl` is outside this plan's declared
  `files_modified` set, the orchestrator's parallel-execution instruction was explicit
  ("Touch ONLY the files your plan declares"), and plan 11-03 was executing concurrently against
  the same phase — an undeclared edit risked a merge collision on a file neither plan owns.
  Logged to `deferred-items.md` with the exact failing expression, the blast radius (six test
  files never run), and the mechanical fix.
- **Consequence for this plan's acceptance criteria:** the criterion
  "`julia --project=spike spike/test/runtests.jl` exits 0" is **NOT MET**, and is reported as such
  rather than worked around. Every other acceptance criterion in all three tasks is met.
- **Files modified:** `.planning/phases/11-…/deferred-items.md`
- **Committed in:** `b1d6228`

**4. [Rule 2 - Missing critical] The declared-deviations header is pinned by a test, not just written**
- **Found during:** Task 1
- **Issue:** The plan requires the file's ALL-CAPS `DECLARED DEVIATIONS` paragraph to enumerate
  exactly the four `P11_RESEARCH_NET_DEVIATIONS` entries (T-11-16, the D-04 discipline), but
  nothing would have caught that paragraph drifting from the Tier-1 list later.
- **Fix:** Added a testset that asserts `length(P11_RESEARCH_NET_DEVIATIONS) == 4` and that the
  source text contains `DECLARED DEVIATIONS`, the `src/amortized/architecture.jl` attribution, the
  `NOT DROP-IN` / `SHIPPED READ SURFACE` consequence, and NO bare epsilon identifier. Booleans are
  computed before `@test` so a failure prints `false` rather than dumping the whole file.
- **Files modified:** `spike/test/test_p11_architecture.jl`
- **Verification:** 6/6 in that testset.
- **Committed in:** `851d64f`

---

**Total deviations:** 4 (2 bugs in plan-specified assertions, 1 blocking-and-deferred, 1
missing-critical). **No Rule 4 (architectural) situation arose; no checkpoint was hit; no package
was installed.**
**Impact on scope:** Deviations 1 and 2 correct assertions whose literal form was wrong about the
world (floating-point boundary behaviour; the v0.2.1 return shape) — neither weakens what is
proven. Deviation 4 strengthens the pre-registration. Deviation 3 is the only unmet acceptance
criterion and is a pre-existing defect in a file this plan is forbidden to touch.

## Assumption Drift (advisory)

**1. "`build_estimator`'s `D` was the only thing standing between the spike and an 8-marginal flow."**
- **Found during:** Task 1
- **Planned:** the plan (and RESEARCH §D12) frame the D = 8 move as "the positional `D` already
  exists, so only an argument changes".
- **Actual:** true for the ARCHITECTURE, but the spike's `fit_theta_transform`
  (`train_npe.jl:80`) is a plain arity-agnostic z-score with no prior box, so nothing in the spike
  previously coupled theta arity to a bounds tuple. The port INTRODUCES that coupling:
  `fit_p11_theta_transform` hard-requires `size(θ, 1) == length(bounds)`. This is the same coupling
  that broke `src`'s trainer in 11-02 deviation 1.
- **Why it matters:** plan 11-07's datagen must emit exactly 8 theta rows in `sample_prior` field
  order or the transform throws a `DimensionMismatch` — a loud failure, but one worth expecting
  rather than debugging.

**2. "The lambda-ablation tripwire is a cheap unit test."**
- **Found during:** Task 3
- **Planned:** the plan treats it as a test file alongside the others.
- **Actual:** it is the heaviest include chain in the phase — `p11_consts` (which itself loads the
  frozen gate constants through `module GateV2`), `p11_architecture`, `infer.jl` (→ `train_npe.jl`
  → `loader.jl`), `forward.jl`, `contract.jl`, `encode.jl` — and once a net exists it will simulate
  6 image pairs and take 12 posterior passes at 2000 draws. It is a minutes-scale gate, not a
  seconds-scale one.
- **Why it matters:** plan 11-07's "run the tripwire immediately after training" step should be
  budgeted as a real run, and the tripwire's include chain is now known-good (it was exercised in
  full by the exit-1 path, which loads everything before checking for the net).

**3. "The four declared deviations are four independent changes."**
- **Found during:** Task 1
- **Planned:** `P11_RESEARCH_NET_DEVIATIONS` reads as a flat list of four.
- **Actual:** they are not independent for the purposes of the report. Deviations (1) and (2) are
  PRIOR changes that this file only consumes (through `p11_theta_prior_bounds()`); deviation (4)
  is a change whose whole purpose is to REMOVE a confound from measuring deviation (2)'s effect;
  and deviation (3) is the only one that changes the read surface's shape.
- **Why it matters:** the report should not present the four as a uniform list of "things that
  differ" — (4) exists to make (2) measurable, and only (1) and (3) are what make the research net
  non-drop-in.

## Issues Encountered

- **`spike/test/runtests.jl` exits 1** on a pre-existing `test_simulator.jl` theta-arity assertion
  (deviation 3). Diagnosed by direct execution, confirmed pre-existing by `git log` on the two
  files, and logged rather than fixed.
- **Two plan-specified assertions were wrong about the world** (deviations 1 and 2). Both were
  found by running the tests, not by inspection, and both were corrected to the measured truth with
  an in-file comment naming the mechanism.
- **A naive `@test occursin(...)` on a source file dumps the entire file into the failure output.**
  Rewritten to compute the booleans first.

## Known Stubs

**One, declared and pre-registered:** `P11_LAMBDA_ABLATION_FACTOR_FALLBACK = 1.15` in
`spike/test/test_lambda_ablation.jl`. It is an **intentional placeholder**, not an oversight:
`P11_LAMBDA_ABLATION_FACTOR` is a Tier-2 constant that MUST NOT exist until the D-06 probe has run
(11-01's test asserts Tier 2 has not been pre-empted), so the tripwire cannot read a real value
today. The file states in ALL-CAPS prose that **plan 11-05 must supersede it**, and the Tier-2 read
is already wired (`isdefined(@__MODULE__, :P11_LAMBDA_ABLATION_FACTOR) ? … : …`) so the append will
be picked up with no edit here. This does not prevent this plan's goal — the tripwire's job in this
plan is to EXIST and to refuse to pass vacuously, both of which are verified.

No other placeholder value, empty container, mock data path, TODO or unwired component was
introduced. `augment_input`, `encode_lambda`/`decode_lambda`, the ported transform and
`build_p11_estimator` are all fully implemented and exercised.

## Threat Flags

None. This plan adds no network endpoint, no auth path, no schema at a trust boundary, and no file
write outside a test's own read of source text. The registered threats were handled as planned:

| Threat ID | Handling | Evidence |
|---|---|---|
| T-11-15 (`augment_input` before standardization) | Runtime `@assert size(Z128, 1) == 128` inside `augment_input`; `Loader._row_partition`'s 128-row error deliberately NOT "fixed" | `@test_throws AssertionError` on 129-, 142- and 64-row inputs; `spike/data/loader.jl` byte-unchanged |
| T-11-16 (unnamed model change contaminating D-02) | The file's ALL-CAPS `DECLARED DEVIATIONS` header enumerates exactly the four pre-registered entries and is PINNED by a test | `length(P11_RESEARCH_NET_DEVIATIONS) == 4` plus four `occursin` assertions on the header text |
| T-11-17 (a tripwire that passes vacuously) | `test_lambda_ablation.jl` prints a named message and `exit(1)`; no skip path exists in the file | measured exit code **1**; `grep -c '@test_skip'` = `0` |
| T-11-18 (non-reproducible posterior draws) | Deterministic `Random.seed!(ABLATION_GLOBAL_SEED)` re-pinned before every arm, derived from `p11_rng(P11_FIXTURE_COUNTER)` | `grep -c 'Random.seed!'` = `1` at the single read helper both arms go through; the `dead_a == dead_b` byte-equality clause proves the pin works |
| T-11-SC (package-manager installs) | None attempted | `git diff --stat 2aae21f..HEAD -- spike/Project.toml spike/Manifest.toml` empty; `runtests.jl` clause (i) asserts the 16-name dependency set as a literal and re-asserts the NeuralEstimators v0.2.1 UUID pin |

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- **Ready — plan 11-05 (probe + Tier-2 append):** it must append `P11_LAMBDA_ABLATION_FACTOR` to the
  Tier-2 block of `p11_consts.jl`, derived from the probe's measured effect scale. The consuming
  read is already written; no edit to `test_lambda_ablation.jl` is needed.
- **Ready — plan 11-07 (research trainer):** the model surface it needs exists and is proven to
  train. It must (a) emit exactly **8** theta rows in `sample_prior` field order or
  `fit_p11_theta_transform` throws, (b) draw **lambda FIRST** and then `shift ~ Uniform(-lambda, lambda)`
  (R4/Pitfall 1 — the tripwire exists precisely because this is easy to get wrong), (c) standardize
  the 128 rows with the frozen `zt` and only THEN call `augment_input`, and (d) save the trained net
  to `spike/npe/p11_research_npe.jld2`, which is the path the tripwire hard-codes.
- **Blocked until 11-07 — plan 11-08 (SC2 ladder):** the tripwire must be RUN and must PASS before
  any ladder result is meaningful. It has never been run against a net.
- **For the report:** the four declared deviations should be presented with their relationships
  (see Assumption Drift 3), and every number must name which net it came from. Nothing in this plan
  produces a reportable measurement.
- **Left red for someone else:** `spike/test/runtests.jl` exits 1 on the pre-existing
  `test_simulator.jl` theta-arity assertion. See `deferred-items.md`.

## Self-Check: PASSED

- All 4 created artifact paths exist on disk: `spike/npe/p11_architecture.jl`,
  `spike/test/test_p11_architecture.jl`, `spike/test/test_p11_smoke.jl`,
  `spike/test/test_lambda_ablation.jl`.
- All 3 task commits exist in `git log --oneline --all`: `851d64f`, `b1d6228`, `f353fda`.
- `git diff --name-only 2aae21f..HEAD` lists exactly the 4 new files plus `spike/test/runtests.jl`
  and `deferred-items.md` — no `src/` file, no manifest, no orchestrator-owned artifact
  (`STATE.md` / `ROADMAP.md` untouched), and no file owned by the parallel plan 11-03.

---
*Phase: 11-registration-and-chromatic-uncertainty-as-latent*
*Completed: 2026-07-25*
