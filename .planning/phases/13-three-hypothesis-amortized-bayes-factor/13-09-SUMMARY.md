---
status: complete
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 09
subsystem: Phase-11 artifact binding (the single precondition site)
tags: [d-02, d-03, d-01, d-04, no-fallback, frozen-zt, derived-input-width, silent-failure-tripwire, order-independence]

# Dependency graph
requires:
  - phase: 13-three-hypothesis-amortized-bayes-factor
    plan: 01
    provides: "spike/p13/consts.jl Tier-1 pre-registration — P13_IMSIZE_SET/P13_IMSIZE_WEIGHTS (the F5 mixture, READ not retyped), P13_TAU_REFERENCE_LAMBDA_RULE = :widest_rung, P13_LAMBDA_PLACEMENT, P13_USE_GPU"
  - phase: 13-three-hypothesis-amortized-bayes-factor
    plan: 08
    provides: "spike/p13/net.jl — p13_ratio_input_dim, three_way_input_dim, p13_pair_encode with its iseven guard, encode_conditioned_pair; and the Phase-13 include block in spike/test/runtests.jl"
  - phase: 13-three-hypothesis-amortized-bayes-factor
    plan: 10
    provides: "the frozen Tier-2 P13_TAU_REFERENCE_LAMBDA = 3.0 that :widest_rung resolves to"
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 07
    provides: "spike/npe/p11_research_npe.jld2 — the trained registration-aware research NPE, its FROZEN zt and FROZEN P11BoundedThetaTransform; spike/npe/p11_architecture.jl — P11BoundedThetaTransform, encode_lambda, p11_theta_prior_bounds; spike/validation/p11_consts.jl — LAMBDA_MIN/LAMBDA_MAX, P11_IMSIZE_SET/WEIGHTS, SC2_RUNGS"
provides:
  - "spike/p13/preconditions.jl — p13_require_phase11, load_p13_basis, assert_frozen_zt, p13_conditioning_length, p13_encode_lambda, p13_encode_pair, and the guarded include preamble every later spike/p13/*.jl uses"
  - "the SINGLE site binding Phase-11 names: P13_PHASE11_NET, P13_PHASE11_LAMBDA_RANGE, P13_PHASE11_REFERENCE_LAMBDA"
  - "MEASURED: n_cond = 1, derived as length(encode_lambda(3.0)); the Phase-13 input width at G = 8 is three_way_input_dim(8, 1) = 5*64 + 1, never written as a literal"
  - "spike/test/test_p13_preconditions.jl — 74/74 pass with the artifact present, 46 pass + 2 explicit skips with it absent"
  - "the D-02 tripwire: a missing Phase-11 net is a hard error naming Phase 11 and forbidding the grid-8 fallback"
affects: [13-11-datagen, 13-12-gate, 13-13-alpha-series, 13-16-realimage, 13-14-report]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Assert the CONCRETE persisted type, not its field names: JLD2 answers a missing type with a WARNING and a ReconstructedMutable placeholder whose fields match, so a hasproperty check passes while the frozen basis is silently wrong"
    - "Close the input arithmetic instead of asserting the width: d_in == 2*length(zt.mean) + p13_conditioning_length() proves 129 = 128 + 1 without typing either number"
    - "Bind provenance from a frozen pre-registration rather than from net metadata — a pre-registered consts file cannot drift with a retrain and sits inside the audit trail"
    - "Read a shadowed guarded-include file through a private module: when two pre-registrations collide on a guard sentinel and neither may be edited, an isolated read restores order-independence locally"
    - "Split an inherited-versus-re-fit check into a HARD object-identity assertion and a SOFT reported moment heuristic, so the unambiguous half errors and the ambiguous half is reported rather than silently accepted"

key-files:
  created:
    - spike/p13/preconditions.jl
    - spike/test/test_p13_preconditions.jl
  modified:
    - spike/test/runtests.jl

key-decisions:
  - "Task 1's checkpoint was auto-approved on the orchestrator's five measured bindings; every literal in the file was read off the artifact and the frozen consts on disk, none was guessed"
  - "Phase 11 is CLOSED-negative rather than COMPLETE, and plan 13-09 proceeds anyway on a five-point argument written into preconditions.jl as a header comment — nothing in the file asserts Phase 11 is complete"
  - "The lambda range and imsize provenance are bound from spike/validation/p11_consts.jl, not from the handle, because those two properties DO NOT EXIST on the handle and the pre-registration is the stronger source"
  - "assert_frozen_zt requires the handle to have come from load_p13_basis(): a second p13_require_phase11 call is a second read and a different object, and it fails on purpose — the discipline is one load per session, not merely equal numbers"
  - "n_cond is 1, and the chromatic epsilon is NOT a second conditioning input (13-RESEARCH Open Question 4 resolves NO): it is an inferred theta field, row 8 of the frozen box"

requirements-completed: [D-02, D-03, D-01, D-04]

# Metrics
duration: ~35min
completed: 2026-07-28
---

# Phase 13 Plan 09: The Phase-11 Precondition Binding Summary

**`spike/p13/preconditions.jl` now binds the Phase-11 registration-aware research NPE in exactly
one place, asserts all six contract items against the artifact and the frozen pre-registration
rather than trusting them, derives the conditioning width from Phase 11's own encoder
(`n_cond = 1`, closing `d_in = 129 = 2*64 + 1`), and turns a missing net into a hard error that
names Phase 11 and forbids the grid-8 fallback — and along the way it caught two ways the load
was silently wrong.**

## Performance

- **Duration:** ~35 min
- **Completed:** 2026-07-28T10:09:01+02:00
- **Tasks:** 3 (Task 1 checkpoint auto-approved, Tasks 2-3 auto)
- **Files:** 2 created, 1 modified (+831 lines, 0 deletions to pre-existing content)

## Phase-11 bindings

Task 1's deliverable, recorded verbatim. Every item was MEASURED on 2026-07-28 at repo commit
`85a75c2`, by loading the artifact and reading the frozen consts — not by reading a summary
document.

**1. Artifact path:** `spike/npe/p11_research_npe.jld2`

Confirmed present (4 933 362 bytes, git-tracked). It is NOT `spike/npe/trained_npe.jld2` (the
pre-Phase-11 net, which sits in the same directory at 4 900 171 bytes) and NOT under
`artifacts/`. Trained CPU-only per 11-07-SUMMARY at `d_in = 129`, `D = 8`.

**2. Persisted handle field names:** `propertynames(h)` is exactly
`(:estimator, :θzt, :zt, :variant, :d_in, :meta)`.

| Field | Measured |
|---|---|
| `h.estimator` | `PosteriorEstimator` (NeuralEstimators) |
| `h.zt` | `ZScoreTransform{Float64, Vector{Float64}}`, `length(h.zt.mean) == 64` — the FROZEN summary transform. **64, not 128**: the mask rows are bypassed, so only the 64 continuous rows are standardized. |
| `h.θzt` | `P11BoundedThetaTransform{Float64, ZScoreTransform{Float64, Vector{Float64}}}`, `propertynames` `(:lo, :hi, :zt)` — the FROZEN theta transform |
| `h.variant` | `:min` |
| `h.d_in` | `129` |
| `h.meta` keys | `(:pool_dir, :n_pairs, :n_train, :n_val, :dstar, :depth, :width, :num_coupling_layers, :flow_depth, :flow_width, :research_lane, :shipped, :generated)`; `n_pairs = 50000`, `research_lane = true`, `shipped = false` |

**3. Theta arity and field order:** `length(h.θzt.lo) == 8`. Bounds, in `sample_prior` field order
(`ρ_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise,
chromatic_eps`):

```
lo = [-0.99, 0.0, 0.0, 0.6, -3.0, -3.0, 0.0, -0.02]
hi = [ 0.99, 0.2, 0.1, 1.0,  3.0,  3.0, 1.0,  0.02]
```

The 8th field is the chromatic epsilon, bounded ±0.02. Fields 5 and 6 are the shift components at
the D-02 widened ±3.0 half-width.

**4. Conditioning encoder and `n_cond`:** `encode_lambda` is defined at
`spike/npe/p11_architecture.jl:297` as `encode_lambda(lam) = (lam - LAMBDA_MIN) / (LAMBDA_MAX -
LAMBDA_MIN)`. It returns a **scalar `Float64`** (`encode_lambda(3.0) == 1.0`), so

> **`n_cond = 1`, derived as `length(encode_lambda(P13_PHASE11_REFERENCE_LAMBDA))`.**

`129 == 128 + 1` closes the arithmetic. **The chromatic epsilon is NOT a second conditioning
input** — 13-RESEARCH Open Question 4 resolves NO: it is an inferred theta field (item 3), not a
conditioning row. There is ONE conditioning input, consistent with the plan's "one lambda, not
two" ruling. The realized width is never written down: neither `321` nor `129` appears as a
literal in executable code (both asserted absent).

**5. Reference lambda for `:widest_rung`:** Phase-11's `LAMBDA_MAX = 3.0`
(`spike/validation/p11_consts.jl:171`, where the comment records that it equals the widened
`SHIFT_PRIOR` half-width). `LAMBDA_MIN = 0.25` (`p11_consts.jl:166` — deliberately never 0,
because the 0 → 0.25 step is interpolation onset rather than misalignment sensitivity). This
agrees with the already-frozen `P13_TAU_REFERENCE_LAMBDA = 3.0` from plan 13-10, and the
agreement is asserted at runtime rather than assumed.

**Note on `h.meta.pool_dir`:** it records an absolute path inside a since-deleted worktree
(`...\.claude\worktrees\agent-adcdf94a3723622aa\spike\data\cache\p11\c12f8ab2...`) whose hash key
differs from the pool now on disk (`spike/data/cache/p11/bff8550959...`), because the original
pool was destroyed by a worktree cleanup and regenerated. It is treated as a historical provenance
string only; `isdir` is deliberately NOT asserted on it.

## Phase-11 status ruling

Recorded here and, more importantly, as a prominent header comment in `preconditions.jl` itself,
so a future reader who runs `git log`, finds Phase 11 closed-negative, and wonders why Phase 13
bound its net anyway can find the answer in the file rather than in a planning archive.

**Taken literally, Task 1's first verify item FAILS, and that is not papered over:**

- Phase 11 is **CLOSED (2026-07-27), not COMPLETE** (`11-CLOSURE.md`).
- `11-07-SUMMARY.md` carries `status: blocked` — it ends at its Task 3 gate because the
  pre-registered SC1g lambda-ablation tripwire failed (measured Δρ posterior-SD ratios
  1.109 / 0.987 / 0.936 against a required `P11_LAMBDA_ABLATION_FACTOR = 2.502`).
- Plans 11-08 through 11-11 are **SUPERSEDED**, so 11-11 is not and will not be checked off.

**Why this does not block plan 13-09 — the five-point argument, stated so it can stand or fall in
the open:**

1. **The gate that failed was MIS-SPECIFIED, not failed.** `P11_LAMBDA_ABLATION_FACTOR = 2.502`
   was derived as half a ratio of the *registration-induced component alone*, but the tripwire
   applied it to the *total* Δρ posterior SD. Phase 11's independent coverage evidence puts the
   true width ratio at ~1.00 (RMSE 0.13334 at λ_min vs 0.13338 at λ_max, ratio 1.0003).
   11-CLOSURE's own words: "The net's measured ratios were 1.109 / 0.987 / 0.936. The correct
   answer is ~1.00. The net was right, and the gate failed it anyway." **No threshold was edited
   to reach that conclusion.**
2. **The conditioning input is ALIVE and MEASURED.** 11-07-SUMMARY records shift marginals whose
   SDs track lambda with ratios 4.84–9.05 against an ideal prior ratio of 12.0, plus a fourth
   decisive check that the 129th input row is read and used. What does NOT track lambda is the
   ρ_true posterior *width* — which is the finding, not a defect.
3. **What 13-09 actually needs is the six-item artifact contract**, and all six are verified
   present and asserted at runtime. The checkpoint existed to stop artifact names being GUESSED;
   nothing was guessed.
4. **Phase 13 is already designed to absorb this honestly.** Plan 13-12 carries the must_have
   "(D-03) The log-BF pair's response to the registration-uncertainty level is measured and
   reported honestly, **including a flat response as a finding**." Phase 11 having established
   that the width is genuinely flat in λ makes a flat Phase-13 λ response a *pre-authorized,
   expected* outcome, not something to explain away. Proceeding does not launder a negative
   result.
5. **Note precisely what D-02 claims.** It claims Phase 13 trains on Phase 11's registration-aware
   INPUT SURFACE rather than the shipped grid-8 basis, so the comparison is like-for-like. It does
   NOT claim registration is inferable. Phase 11's negative result leaves D-02 intact.

Nothing in the code asserts "Phase 11 is complete", and the no-fallback rule is not weakened: a
missing net is still a hard `error`, never a skip, and no grid-8 path exists.

## Task Commits

| Task | Name | Commit | Files |
|---|---|---|---|
| 1 | Confirm Phase 11 has landed and bind its realized artifact names | *(checkpoint, auto-approved — no commit)* | — |
| 2 | Write the Phase-11 precondition binding and its loud failure | `8a22a29` | `spike/p13/preconditions.jl` (528 lines) |
| 3 | Write the precondition test and wire it into the gate | `cac93af` | `spike/test/test_p13_preconditions.jl` (298 lines), `spike/test/runtests.jl` (+5, −0) |
| — | Rule-1/3 fix discovered by Task 3's suite-order run | `071b9c3` | `spike/p13/preconditions.jl` (+62, −13) |

## Files Created/Modified

- **`spike/p13/preconditions.jl`** (577 lines) — AGPL header; a capitalized banner stating that
  Phase 13 trains on Phase 11's registration-aware FROZEN `zt` and not the shipped grid-8 one,
  that there is no fallback, that a missing net is a BLOCK, and that this is the single Phase-11
  binding site; the one-lambda-not-two modelling statement; the DECOUPLING line; the Phase-11
  status ruling; the three deviations. Then the isolated `module _P11C` read of Phase 11's Tier-1
  pre-registration, the guarded order-dependent include preamble, and:
  - `P13_PHASE11_NET`, `P13_PHASE11_LAMBDA_RANGE`, `P13_PHASE11_REFERENCE_LAMBDA`
  - `p13_require_phase11(path)` — the block, then the six contract assertions plus the row
    arithmetic, then an assert-then-print provenance line; returns the handle
  - `load_p13_basis(; path)` — the memoized single-load accessor (the `harness.jl:70-80`
    `load_frozen_model` analog, pointed at the Phase-11 artifact)
  - `p13_conditioning_length()`, `p13_encode_lambda(lambda)`, `p13_encode_pair(Zs, Zc, lambda)`
  - `assert_frozen_zt(basis, Zstd)` — hard `zt` object identity plus the soft Pitfall-5 moment
    warning sign, with vector and matrix methods
- **`spike/test/test_p13_preconditions.jl`** (298 lines) — eight testsets, 74 assertions.
- **`spike/test/runtests.jl`** — one additive include plus a four-line rationale comment, placed
  after `test_p13_calibration.jl` and before the throwing correction arm. `git diff --numstat`
  reports **5 insertions, 0 deletions**.

## Verification (observed output, not paraphrase)

Every command below was RUN and its real output read.

| Check | Observed |
|---|---|
| `julia --project=spike spike/test/test_p13_preconditions.jl` | **74 pass / 74 total, 0 fail**, exit 0 |
| Same file with the artifact-guarded branch forced absent | **46 pass, 2 Broken (explicit `@test_skip`), 0 fail**, exit 0 — and `absence is a loud, instructional block` passes 6/6 with no skip |
| `p13_require_phase11()` | prints `Phase-11 basis bound OK: p11_research_npe.jld2 \| research_lane=true shipped=false \| d_in=129 = 2*64 + 1 \| G=8 \| theta arity=8 \| lambda in (0.25, 3.0) \| imsize provenance ((512,512),(1024,1024),(1376,1028),(2048,2048)) w=(0.4,0.25,0.25,0.1)` |
| `p13_require_phase11("does_not_exist.jld2")` | throws; message contains `BLOCKED on Phase 11`, `13-CONTEXT D-02`, `grid-8` → prints `block ok` |
| `p13_conditioning_length()` | `1` (matches the Task-1 binding) |
| `grep -v '^\s*#' … \| grep -c 'trained_npe.jld2\|artifacts/\|grid_8'` | `0` |
| `grep -v '^\s*#' … \| grep -c 'try'` | `0` |
| `grep -v '^\s*#' … \| grep -c '321'` | `0` |
| `grep -c 'test_p13_preconditions.jl' spike/test/runtests.jl` | `1` |
| `grep -c 'include(joinpath(@__DIR__, "test_p13' spike/test/runtests.jl` | `10` (nine from 13-08 plus this one) |
| `grep -c 'test_p13' spike/test/runtests.jl` | `10` (≥ 10 as required) |
| `git diff --numstat -- spike/test/runtests.jl` | `5  0` — zero deletions |
| `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` | exit `0` — byte-unchanged |
| Phase-13 include block in suite order through this gate | **1077 pass / 0 fail** (1 pre-existing Broken in `test_p13_real.jl`), exit 0 |

### What was NOT verified, and why

**`julia --project=spike spike/test/runtests.jl` does NOT exit 0 from this worktree**, and the
cause is upstream of this plan. The suite aborts at `spike/test/test_bf.jl:62` with
`resolve_cache_dir: expected exactly one main-pool cache dir, found 0`. `test_bf.jl` loads a
persisted `trained_ratio.jld2` when one exists and otherwise falls into a tiny-retrain path that
needs the main-pool cache. Both of those are gitignored bulk artifacts, and `git worktree add`
does not copy gitignored files — so this worktree has only `spike/data/cache/fixture`, while the
primary checkout has `spike/data/cache/eca37f54…`, `spike/data/cache/p11` and
`spike/validation/trained_ratio.jld2`. In the primary checkout `test_bf.jl` takes the persisted
path and never calls `train_ratio`.

Evidence that it is pre-existing rather than caused here: `git diff --name-only 85a75c2 HEAD --
spike/test/test_bf.jl spike/validation/train_ratio.jl` is **empty**, and the abort is at
`runtests.jl:175` — before the Phase-13 block this plan appends to. As a substitute, the entire
Phase-13 include block was run in suite order through this gate (1077 pass / 0 fail, above), which
covers the composition risk the full-suite run would have covered. **The full suite should be
re-run in the primary checkout after merge**; it will still turn red at the pre-existing
`test_p13_correction.jl` measured misses, which plan 13-07 committed deliberately.

## Deviations from Plan

### Auto-fixed and forced deviations

**1. [Task 2 — forced by the artifact] `p11_architecture.jl` must be in scope BEFORE `load_npe`**

- **Found during:** Task 2, before writing the file (flagged by the orchestrator, then reproduced)
- **Issue:** `P11BoundedThetaTransform` is defined in `spike/npe/p11_architecture.jl`. Loading the
  artifact without it in scope does **not** fail. JLD2 emits `Warning: type
  Main.P11BoundedThetaTransform{...} does not exist in workspace; reconstructing` and returns
  `JLD2.ReconstructedMutable{..., (:lo, :hi, :zt), ...}` — a placeholder whose field names match,
  so `hasproperty(h, :θzt)` passes and the frozen theta basis is silently not the real one.
  Reproduced verbatim in this worktree.
- **Fix:** the include preamble orders `p11_architecture.jl` first, and `p13_require_phase11`
  asserts the **concrete type** (`h.θzt isa P11BoundedThetaTransform`) with an error message that
  names the cause and the remedy. The test asserts the same thing, so a future regression is loud.
- **Impact:** strictly stronger than the plan's `hasproperty` check. This is exactly the class of
  silent-wrong-basis failure the plan exists to prevent.

**2. [Task 2 — forced by the artifact] `h.meta` carries neither `lambda_range` nor
`training_imsize_provenance`**

- **Found during:** Task 2
- **Issue:** the plan's `<interfaces>` skeleton asserts `h.lambda_range isa Tuple` and
  `h.training_imsize_provenance.recorded == true`. **Neither property exists** on this handle;
  both asserts would throw. The skeleton was written before Phase 11 executed.
- **Fix:** both are bound from Phase 11's Tier-1 pre-registration instead —
  `LAMBDA_MIN`/`LAMBDA_MAX` for the range, `P11_IMSIZE_SET`/`P11_IMSIZE_WEIGHTS` for the
  provenance — with the argument stated in a comment: **a frozen, pre-registered consts file is a
  STRONGER provenance source than net metadata**, because it cannot drift with a retrain and it
  sits inside the pre-registration audit trail, whereas a metadata field is whatever the trainer
  happened to write. What *can* be read off the handle is asserted off the handle: theta arity 8,
  the positional prior box, the ±3.0 shift half-width, `variant === :min`, the row arithmetic and
  the `research_lane`/`shipped` flags — so the binding is not merely documentary.
- **Impact:** the F5 covariate-shift check and the `:widest_rung` resolution are both still
  asserted, from a source that is harder to drift.

**3. [Rule 1 + Rule 3 — a real bug in the first commit] the two pre-registrations collide on a
guard sentinel, so the naive include ordering passed standalone and failed in the suite**

- **Found during:** Task 3, running the file after `test_p13_consts.jl` as `runtests.jl` does
- **Issue:** `spike/p13/consts.jl` RESERVES the name `P11_DEV_SEED` in its forbidden-seed block
  (so two research lanes provably cannot share a Philox key) — and that name is **precisely** the
  guard sentinel `spike/validation/p11_consts.jl` keys its whole Tier-1 block on. Once the
  Phase-13 pre-registration is loaded into a module, `p11_consts.jl` believes it has already
  installed itself there, skips its Tier-1 block, and then dies on its own out-of-guard
  self-check with `UndefVarError: SC2_SPEARMAN_ATTENUATION`. Commit `8a22a29` therefore passed
  standalone and would have failed inside `runtests.jl` — the worst shape of failure to ship.
- **Fix (commit `071b9c3`):** neither file may be edited (both are frozen pre-registrations, and
  the `P11_DEV_SEED` reservation is there for a good reason), so the repair is entirely local to
  `preconditions.jl`: `p11_consts.jl` is read through a private `module _P11C` — the same
  precedent as `module _GC` in `p13/consts.jl` and `module GateV2` in `p11_consts.jl` itself — and
  `LAMBDA_MIN`, `LAMBDA_MAX`, `SC2_RUNGS`, `P11_IMSIZE_SET`, `P11_IMSIZE_WEIGHTS` are bound from
  there, guarded. With the two lambda bounds already in scope, `p11_architecture.jl`'s own guard
  makes it skip `p11_consts.jl` and the shadowed-sentinel path is never taken. Measured first:
  `p11_architecture.jl`'s executable code references exactly `LAMBDA_MIN` and `LAMBDA_MAX` from
  that file, and nothing else in the include chain references it at all.
- **Impact:** the file is now **order-independent** — verified passing 74/74 with
  `p13/consts.jl` loaded first, with `p11_consts.jl` loaded first, and standalone. That is the
  property the guarded-include idiom is supposed to buy everywhere else in the spike, so the fix
  restores an invariant rather than adding a constraint. It also rejected the tempting alternative
  of making `runtests.jl` order-dependent.

**4. [Naming] `P13_TAU_REFERENCE_LAMBDA_EXPECTED` does not exist; `P13_TAU_REFERENCE_LAMBDA` was
used**

- **Found during:** Task 3
- **Issue:** the plan's Task-3 action names a constant `P13_TAU_REFERENCE_LAMBDA_EXPECTED` for the
  lambda-range upper-end assertion. No such constant exists; plan 13-10 froze the value as
  `P13_TAU_REFERENCE_LAMBDA = 3.0`.
- **Fix:** the assertion compares against `P13_TAU_REFERENCE_LAMBDA`, re-exposed locally as
  `P13_PHASE11_REFERENCE_LAMBDA` and asserted equal to Phase-11's live `LAMBDA_MAX`.
- **Impact:** none on substance — the criterion the plan asked for is asserted, under the name the
  repository actually uses. τ is read through `p13_tau()` conventions and `0.15` is never inlined
  by this plan (it is never referenced at all).

**5. [Task 2 — minor] `p13_conditioning_length` derives from `encode_lambda`, not from
`p13_encode_lambda`**

- **Found during:** Task 2
- **Issue:** the plan's item 6 specifies `p13_conditioning_length() =
  length(p13_encode_lambda(<lambda>))`, while its item 5 specifies that `p13_encode_lambda`
  asserts its return length **against** `p13_conditioning_length()`. Taken together those are
  mutually recursive and would stack-overflow.
- **Fix:** `p13_conditioning_length` calls Phase 11's `encode_lambda` directly (which is what
  "derived from the encoder's return length" means anyway) and memoizes; `p13_encode_lambda` keeps
  its assertion. Documented inline.
- **Impact:** none. The success criterion "`n_cond` is DERIVED via `length(encode_lambda(...))`"
  is met literally.

---

**Total deviations:** 5 — two forced by the artifact as built (1, 2), one genuine bug caught and
fixed (3), two naming/wiring corrections (4, 5).
**Impact on plan:** none on substance. Every must_have and every success criterion is met, and
deviation 3 makes the deliverable strictly more robust than the plan specified.

## Assumption Drift (advisory)

**1. The plan assumed the Phase-11 handle would self-describe its training joint.**

- **Found during:** Task 2
- **Planned:** the persisted handle would carry `lambda_range` and
  `training_imsize_provenance.recorded`, so the precondition could interrogate the artifact alone.
- **Actual:** the handle carries only `(:estimator, :θzt, :zt, :variant, :d_in, :meta)`, and
  `meta` is training-recipe bookkeeping (pool dir, pair counts, layer widths) rather than
  scientific provenance. The training joint had to be bound from the frozen consts file instead.
- **Why it matters to a reader:** the resulting guarantee has a different *shape* from the one the
  plan described. It is not "the artifact declares its own joint and we check the declaration"; it
  is "the pre-registration declares the joint, and the artifact is checked for the properties that
  *are* observable on it." The former would catch a net trained against an edited consts file; the
  latter would not. That gap is real, and it is bounded by the fact that `consts.jl` is
  fingerprinted by sha256 in the τ probe artifact and by git history.

**2. The plan assumed guarded includes made ordering an optimization.**

- **Found during:** Task 3
- **Planned:** copy `harness.jl`'s guarded preamble; ordering documents dependencies but is not
  load-bearing (`runtests.jl`'s own comment says exactly this).
- **Actual:** for this particular pair of files the ordering is load-bearing, because the two
  pre-registrations collide on a guard sentinel and one file's out-of-guard tail depends on its
  own skipped block. The idiom's promise fails precisely where two frozen files reserve the same
  name.
- **Why it matters to a reader:** it means "guarded include ⇒ order-free" cannot be assumed
  anywhere else in the spike either, whenever two `consts.jl` files share a reserved name.
  Deviation 3 restored order-independence for this file specifically; it did not fix the general
  property.

## Issues Encountered

- **The full suite cannot be run to completion from a worktree.** See "What was NOT verified"
  above: gitignored bulk caches are absent, so `test_bf.jl`'s fallback retrain path aborts the
  include chain before the Phase-13 block. Not caused by this plan; the abort site is
  byte-untouched.
- **`spike/data/cache/p11` is absent from this worktree** (gitignored, not copied by
  `git worktree add`). It remains intact in the primary checkout. Nothing in this plan reads or
  writes it.
- **The Phase-11 artifact was moved to the scratchpad and restored** while attempting to
  demonstrate the artifact-absent branch empirically; the follow-up run was blocked by the
  permission layer, so it was restored immediately (verified: same 4 933 362 bytes, `git status`
  reports no modification) and the branch was instead exercised in memory by forcing the `isfile`
  guard false. No bulk artifact was regenerated or destroyed.

## Known Stubs

None. Every function in `preconditions.jl` is fully implemented and exercised by the test; the
only `@test_skip` in the suite is the deliberate artifact-absent branch, which is a documented
skip rather than a stub and does not fire in the current repository state.

## Threat Flags

None. This plan adds no network endpoint, no auth path, no schema and no new file-access pattern.
It reads one existing project-produced `.jld2` from a `@__DIR__`-derived path (T-13-04 mitigation:
the handle's expected properties are asserted immediately after load) and one frozen consts file.
The register's four `mitigate` dispositions for this plan are all implemented: T-13-08 (hard
`error` plus the no-fallback source assertions), T-13-30 (`assert_frozen_zt`), T-13-31 (imsize
provenance asserted equal), T-13-07 (`p13_conditioning_length` derives the width; the `iseven`
guard is retained and tested). T-13-SC holds: no package was installed and the manifests are
byte-unchanged.

## Self-Check: PASSED

- `spike/p13/preconditions.jl` — FOUND (577 lines)
- `spike/test/test_p13_preconditions.jl` — FOUND (298 lines)
- `spike/test/runtests.jl` — FOUND, modified additively (+5, −0)
- commit `8a22a29` — FOUND
- commit `cac93af` — FOUND
- commit `071b9c3` — FOUND
- `src/`, `spike/Project.toml`, `spike/Manifest.toml` — byte-unchanged (`git diff --quiet` exit 0)

## Next Phase Readiness

- **Every later Phase-13 file now has one door to Phase 11.** Plans 13-11 onward include
  `spike/p13/preconditions.jl` and call `load_p13_basis()`; none may hard-code an artifact path,
  and `_P11C.NAME` is the documented route to any other Phase-11 Tier-1 constant.
- **The input width is settled and must stay derived.** `three_way_input_dim(8,
  p13_conditioning_length())`. The datagen and training plans must not write the realized number.
- **`assert_frozen_zt` expects the memoized handle.** Pass `load_p13_basis()`, not a fresh
  `p13_require_phase11()` result — a second read is a different object and fails the identity
  check by design.
- **A flat λ response in 13-12 is pre-authorized, not a failure.** Phase 11 established the width
  is genuinely flat in λ (ratio 1.0003); 13-12's must_have already admits a flat response as a
  finding. Report it as one.
- **Re-run the full suite in the primary checkout after merge.** It will still be red at the
  pre-existing `test_p13_correction.jl` measured misses (0.5498 and 0.3783 against a 0.25
  ceiling), which plan 13-07 committed deliberately and which this plan does not touch.

---
*Phase: 13-three-hypothesis-amortized-bayes-factor*
*Completed: 2026-07-28*
