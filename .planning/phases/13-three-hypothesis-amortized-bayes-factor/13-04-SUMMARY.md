---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 04
subsystem: spike-p13-alpha-series
tags: [d-15, d-16, semi-synthetic, image-transform, invariants, pitfall-2]
requires:
  - spike/p13/consts.jl          # P13_ALPHA_* pre-registration (Tier-1, frozen)
  - spike/contract.jl            # build_mci, patch_summary, the read-only src/ include chain
  - spike/data/encode.jl         # encode_d01, the 128-row vals/mask layout
  - spike/simulator/forward.jl   # simulate_pair (FIXTURE ONLY, not a dependency of the transform)
provides:
  - alpha_background_floor
  - ch1_object_mask
  - alpha_segregate
  - alpha_ladder
  - alpha_ladder_summaries
  - verify_alpha_invariants
affects:
  - 13-13   # run_alpha_series.jl, the SIMULATED arm (reported)
  - 13-16   # run_p13_realimage.jl, the REAL arm (qualitative, non-gating)
  - 13-08   # runtests.jl wiring of spike/test/test_p13_*.jl
tech-stack:
  added: []            # NO package installed; spike/Project.toml + Manifest.toml byte-unchanged
  patterns:
    - "graded perturbation family with a continuous knob (spike/validation/ood.jl:477-502)"
    - "zero-set-preserving transform (negctrl_affine, spike/validation/ood.jl:513-528)"
    - "reporting invariance verifier (verify_summary_invariance, ood.jl:567-581)"
    - "frozen already-shipped Otsu mask rule (_calculate_mask, src/LoadImages.jl:235-241)"
    - "comment-stripped source-grep gate (test_p13_consts.jl:49-50)"
key-files:
  created:
    - spike/p13/alpha_series.jl        # 433 lines
    - spike/test/test_p13_alpha.jl     # 344 lines, 129 assertions
  modified: []
decisions:
  - "The declared dynamic range the value bound is stated in is a FOURTH substrate-dependent knob, so P13_ALPHA_MAX_VALUE_BOUND became the DEFAULT of a keyword rather than an inlined absolute post-condition"
  - "Every guard throws ArgumentError, never @assert, because @assert is elided under -O3 and the plan's own acceptance criteria demand ArgumentError"
  - "verify_alpha_invariants REPORTS mask_fraction_ok / max_value_ok instead of throwing; enforcement lives in alpha_segregate"
  - "The test fixture is 128x128, not 64x64: at 64x64 the Otsu mask covers 0.66% of the frame, below the pre-registered lower bound"
metrics:
  duration: ~45 min
  completed: 2026-07-27
  tasks: 2
  commits: 2
  files: 2
---

# Phase 13 Plan 04: Alpha-Graded Random-to-Exclusion Series Summary

An intensity-conserving, zero-set-preserving, **spatially** constructed alpha ladder that moves
ch2 signal out of the ch1 Otsu object mask into its complement — bitwise the input at
`alpha = 0` with no special case, monotone in `m-bar`, and proved by 129 executable assertions
on a simulated fixture in 11 s CPU-only.

## What Was Built

`spike/p13/alpha_series.jl` (433 lines) defines six flat top-level functions, all with contract
docstrings citing D-15/D-16 and Pitfall 2:

| Function | Role |
|---|---|
| `alpha_background_floor(y; q = P13_ALPHA_BG_QUANTILE)` | `quantile(vec(y), 0.05)`, `ArgumentError` naming `_exclude_zero` if not strictly positive |
| `ch1_object_mask(pair)` | `_calculate_mask(build_mci(pair))[1]` — the frozen shipped Otsu rule, computed once from the unmodified ch1 |
| `alpha_segregate(x, y, M, alpha; b, max_value)` | the transform, typed on `AbstractMatrix{Float64}`, no substrate branch |
| `alpha_ladder(pair; ladder, mask, floor, max_value)` | `(pairs, mask, floor, ladder)`; mask and floor are ladder-level constants |
| `alpha_ladder_summaries(pair; ladder, max_value)` | per rung `(alpha, mbar, n_missing, max_value, sum_ratio, n_zero)` |
| `verify_alpha_invariants(pair; ladder, rtol, max_value)` | 8 booleans + 4 measured residuals, reports rather than enforces |

`m-bar` is the **mask-weighted** mean of the present continuous summary rows (rows 1:64 of
`encode_d01` weighted by rows 65:128). Unweighted would let `encode_d01`'s missing-to-0
imputation manufacture a decaying ladder out of missingness.

No plotting function was defined (figures belong to 13-13 / 13-16). No `corpus/` path and no
`open_sealed_holdout` call exists anywhere in the file — asserted by a comment-stripped
source grep in the suite.

## Measured Results

### The `m-bar(alpha)` curve on the fixture

Fixture: `simulate_pair(p13_fix_rng(P13_FIXTURE_COUNTER), theta; imsize = (128,128))` at
`rho_true = 0.0`, `spillover = 0.05`, `autofluorescence = 0.02`, `label_efficiency = 0.9`,
`shift_dx = shift_dy = 0.0`, `noise = 0.3`, `chromatic_eps = 0.0`.
Ladder-level constants: mask fraction **0.28638** (in the frozen band `(0.01, 0.40)`),
floor `b = 0.181661`, `maximum(y) = 4.84488`, source zero count **1** of 16 384.

| alpha | `m-bar` | `n_missing` | `max(y_alpha)` | `sum_ratio` | `n_zero` |
|---|---|---|---|---|---|
| 0.000 | **+0.103299** | **0** | 4.844881 | 1.0 | 1 |
| 0.125 | +0.029803 | 0 | 4.261978 | 1.0 | 1 |
| 0.250 | −0.043034 | 0 | 4.401067 | 1.0000000000000002 | 1 |
| 0.375 | −0.112201 | 0 | 4.573875 | 1.0 | 1 |
| 0.500 | −0.174874 | 0 | 4.746682 | 1.0 | 1 |
| 0.625 | −0.228604 | 0 | 4.919490 | 1.0000000000000002 | 1 |
| 0.750 | −0.272831 | 0 | 5.092297 | 1.0 | 1 |
| 0.875 | −0.309004 | 0 | 5.265105 | 1.0 | 1 |
| 1.000 | **−0.336522** | **0** | 5.437912 | 0.9999999999999997 | 1 |

**Missing-patch count: 0/64 at the first rung and 0/64 at the last rung** — the Pitfall-2
warning sign never fires. The curve is strictly monotone decreasing and **crosses zero between
alpha = 0.125 and alpha = 0.250** (linear interpolation gives `alpha* ≈ 0.176`), so the sign
transition is well inside the ladder rather than pinned at an endpoint. `max_rel_intensity_err`
= **2.92e-16**, eleven orders of magnitude under the frozen `P13_ALPHA_INVARIANT_RTOL = 1e-10`.

`verify_alpha_invariants` reports all eight booleans `true`.

### The Gray-range companion arm

The same pair rescaled by one positive scalar to `maximum = 0.5` (mirroring the real fixtures'
measured regime) reproduces the `m-bar` curve to **~1e-16** at every rung — direct evidence that
the per-patch Pearson summary is exactly affine-invariant, so the value bound is a substrate
*declaration* and not part of the mechanism. On that arm the **default** (pre-registered) bound
applies and every rung stays inside it: `realized_max_value = 0.561202 <= 1.0`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] `P13_ALPHA_MAX_VALUE_BOUND` is a real-substrate range, so it became a keyword default rather than an inlined absolute post-condition**

- **Found during:** Task 1, confirmed by measurement before a line of the transform was written.
- **Issue:** The plan and 13-VALIDATION.md require the post-condition
  `maximum(y_out) <= P13_ALPHA_MAX_VALUE_BOUND` (= 1.0), and Task 2 requires that post-condition
  to hold **on a simulated fixture**. It cannot. The constant's measured basis is the committed
  real Gray TIFFs, whose values live in `[0,1]` and whose post-boost maximum stayed at or below
  0.75. The simulator's intensities are unnormalized softplus outputs and measure **4.84 at
  128x128 and 5.34 at 256x256 before any alpha is applied**, so an absolute 1.0 bound rejects
  every simulated input. The two requirements are mutually unsatisfiable as literally written.
- **Fix:** `max_value::Real = P13_ALPHA_MAX_VALUE_BOUND` is a keyword on `alpha_segregate`,
  threaded through `alpha_ladder`, `alpha_ladder_summaries` and `verify_alpha_invariants`. The
  frozen constant remains the **default**, so the real arm (13-16) gets the pre-registered bar
  with no call-site knob and the default code path is the one under test on the Gray-range
  companion fixture. The simulated arm declares no range explicitly at the call site, and the
  realized maximum is **still recorded** in `alpha_ladder_summaries(...).max_value` and
  `verify_alpha_invariants(...).realized_max_value`. Nothing is clamped, and no auto-derived
  bound was introduced — a bound derived from `(maximum(y), boost)` would be mathematically
  vacuous, since `maximum(y_out) <= boost * maximum(y)` holds identically.
  `spike/p13/consts.jl` is **byte-unchanged**.
- **Why this is not a pre-registration breach:** the pre-registered value is neither edited nor
  relaxed; it is the default and it is enforced on the substrate it was measured on. What moved
  into the caller is the *declaration of which dynamic range applies* — which the plan's own
  D-15 ruling already places in caller-supplied provenance ("the three substrate-dependent knobs
  ... belong in a provenance record supplied by the caller, never in the algorithm"). This is a
  fourth such knob that measurement discovered.
- **Files modified:** `spike/p13/alpha_series.jl`, `spike/test/test_p13_alpha.jl`
- **Commits:** `9422620`, `472b463`

**2. [Rule 1 - Bug] Every guard throws `ArgumentError`; no `@assert` survives from the research snippet**

- **Found during:** Task 1.
- **Issue:** 13-RESEARCH §J3's snippet uses `@assert` for `b > 0` and `S_out > 0`. The plan's own
  Task-2 testset 6 requires `@test_throws ArgumentError alpha_segregate(x, y, M, 0.5; b = 0.0)`,
  which an `AssertionError` does not satisfy. Independently, `@assert` is elided at `-O3`, so the
  single guard standing between this phase and a silently flat ladder would vanish in exactly the
  configuration a long reported run would use.
- **Fix:** all five entry guards plus the post-condition throw `ArgumentError` with the reasoning
  in the message (`_exclude_zero` named for the floor, "nowhere to redistribute to" for the
  complement, the measured boost table for the fraction band, the realized value for the range).
- **Commit:** `9422620`

**3. [Rule 2 - Missing validation] Three input-validation guards the plan did not enumerate**

- **Found during:** Task 1.
- **Added:** frame-shape agreement across `x`, `y` and `M`; `q` in `(0,1)` in
  `alpha_background_floor`; `length(pair) == 2` in `alpha_ladder`. A mismatched frame would
  otherwise broadcast into a silently wrong mask application (T-13-14 is a semantic-tampering
  threat, so a shape mismatch that does not throw is a mitigation gap, ASVS V5).
- **Commit:** `9422620`

**4. [Rule 2 - Correctness] `verify_alpha_invariants` reports an out-of-band mask instead of throwing**

- **Found during:** Task 1.
- **Issue:** the plan requires a `mask_fraction_ok` **boolean field**. If the verifier built the
  ladder through the enforcing path, an out-of-band mask would raise inside `alpha_segregate` and
  the field could never be `false` — the reported invariant would be unobservable.
- **Fix:** the verifier evaluates the fraction itself and short-circuits to a report with
  `mask_fraction_ok = false` and the realized fraction; it builds the ladder with the range check
  disabled and then compares `realized_max_value` against `max_value` itself. Enforcement stays
  in `alpha_segregate`; the verifier's job is to say what was measured. Asserted in testset 10.
- **Commit:** `9422620`

**5. [Rule 3 - Blocking] The test fixture is 128x128; 64x64 trips a pre-registered guard**

- **Found during:** Task 2 (measured in Task 1's pre-write probe).
- **Issue:** the plan offers "64 by 64 or 128 by 128". At 64x64 the Otsu mask of this theta covers
  **0.66%** of the frame, below `first(P13_ALPHA_MASK_FRACTION_BOUNDS) = 0.01`, so the guard
  correctly refuses the fixture.
- **Fix:** 128x128 (mask 28.6%, comfortably in band). A fixture that trips a pre-registered guard
  is a bad fixture, not a reason to move the bound; the in-file comment says so explicitly.
- **Commit:** `472b463`

**6. [Additive] A Gray-range companion fixture so the pre-registered bound's own code path is tested**

- **Found during:** Task 2.
- **What:** the simulated pair rescaled by one positive scalar to `maximum = 0.5`. A positive
  affine rescale leaves every per-patch Pearson correlation exactly unchanged (verified: the two
  `m-bar` curves agree to ~1e-16), so this is the same ladder on a substrate that *has* a
  declared range. It lets the **default** `max_value = 1.0` path be exercised for real, rather
  than only through a hand-built toy — which is what keeps deviation 1 honest. The plan's real
  substrate (`test/test_images/`) is owned by plans 13-15/13-16 and was deliberately not read.
- **Commit:** `472b463`

## Assumption Drift (advisory)

**1. The naive all-positive invariant is false on the SIMULATED substrate too, not only on real TIFFs**

- **Found during:** Task 1 pre-write measurement.
- **Planned:** 13-RESEARCH §J5.3 finding 4 and Pitfall 2b both state the all-pixels-strictly-positive
  property "holds for simulator output (`BG_FLOOR = 0.02` guarantees positivity)" and fails only on
  the committed real TIFFs.
- **Actual:** the simulated fixture carries **1 exact-zero ch2 pixel of 16 384** at 128x128 and
  **8 of 65 536** at 256x256, so `minimum(y) == 0.0` on simulated input as well. `BG_FLOOR` is
  applied at stage 5 and the stage-7 Poisson/Gaussian noise clamp can push a pixel back to exactly
  zero afterwards.
- **Why it matters:** it *strengthens* the locked rule rather than weakening it. Writing the naive
  form would have failed a correct algorithm on **both** substrates, not just one, so
  `P13_ALPHA_ZERO_INVARIANT = :count_iszero_preserved` is even more load-bearing than the
  pre-registration argued. Recorded in the suite's own comments so the next reader does not
  re-derive it.

**2. `P13_ALPHA_MAX_VALUE_BOUND` is not substrate-agnostic**

- **Planned:** the plan and 13-VALIDATION.md treat the bound as an absolute post-condition of the
  algorithm, alongside three named substrate-dependent knobs.
- **Actual:** it is a fourth substrate-dependent knob. See deviation 1 for the resolution.

**3. The alpha ladder's reach on simulated substrate is narrower than the real ladder's**

- **Planned:** 13-RESEARCH §J5.3 finding 3 records the real ladder spanning
  `m-bar` ∈ [−0.25, +0.33] and warns the series does not probe the deep-exclusion tail.
- **Actual:** the simulated fixture spans `m-bar` ∈ [−0.337, +0.103] — the exclusion end reaches
  slightly *further* than the real arm, but the alpha = 0 end starts near random by construction
  rather than at +0.33. The scope limit ("this series tests the sign transition and the
  near-exclusion regime, not `rho ≈ −0.99`") therefore holds on both arms and must still be
  stated in the report. The two arms remain separate experiments and must never be averaged.

## Verification Performed

Every command below was RUN and its real output observed.

| Check | Command | Observed |
|---|---|---|
| suite green | `julia --project=spike spike/test/test_p13_alpha.jl` | **129 pass / 0 fail / 0 error**, exit 0 |
| wall clock | `time` on the above | **11.5 s** (budget 60 s) |
| includable from a harness | `julia --project=spike -e 'include("spike/test/test_p13_alpha.jl")'` | exit 0 |
| transform loads | `julia --project=spike -e 'include("spike/p13/alpha_series.jl"); println("loaded")'` | `loaded` |
| bitwise identity, synthetic | one-liner, quarter-field mask | `bitwise=true` |
| alpha = 1 synthetic invariants | same one-liner | zeros preserved `true`, `sum` rtol 1e-10 `true`, `max = 0.7269 <= 1.0` `true` |
| mask-fraction guard, empty mask | `alpha_segregate(zeros(10,10), ones(10,10), falses(10,10), 0.5; b=0.1)` | `fraction guard ok` |
| full-field mask guard | `alpha_segregate(zeros(4,4), ones(4,4), trues(4,4), 0.5; b=0.1)` | `guard ok` |
| forbidden invariant absent | `grep -c 'all(.*\.> 0)' spike/p13/alpha_series.jl` | `0` |
| corrected field present | `grep -c 'no_new_zeros' ...` / `grep -c 'all_positive' ...` | `4` / `0` |
| no sealed-holdout reach | `grep -v '^\s*#' ... \| grep -c 'open_sealed_holdout\|corpus/'` | `0` |
| no zero-alpha special case | `grep -v '^\s*#' ... \| grep -c 'alpha == 0\|iszero(alpha)'` | `0` |
| comment-strip idiom in suite | `grep -c 'startswith(strip' spike/test/test_p13_alpha.jl` | `1` |
| exact `==` at the alpha-0 rung | `grep -A3 'BITWISE' spike/test/test_p13_alpha.jl` | shows `== ALPHA_FIX_Y` |
| decoupling | `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` | exit 0 |
| scope | `git diff --name-only cc90d79 HEAD` | exactly the two planned files |

`spike/p13/consts.jl`, `spike/test/runtests.jl`, `.planning/STATE.md` and `.planning/ROADMAP.md`
were not touched. No package was installed.

## Known Stubs

None. Every function is fully wired; the only deliberately deferred surface is the reporting
runner (`m-bar` / `log BF` curves and the `alpha*` crossing figure), which the plan explicitly
assigns to 13-13 (simulated arm) and 13-16 (real arm).

## Threat Flags

None. No new network endpoint, auth path, file-access pattern or trust-boundary schema was
introduced. The two boundaries in the plan's register are both mitigated as specified:
`unmodified image -> transformed ladder` (T-13-14/15/16/53) by the strictly-positive floor, the
`count(iszero, ...)` comparison against the unmodified channel, ladder-level mask/floor
constants, intensity conservation at rtol 1e-10 and the two named range/band guards; and
`Phase-13 code -> corpus sealed holdout` (T-13-05) by the comment-stripped source-grep testset.
Every path is `joinpath(@__DIR__, ...)` with no interpolated strings (T-13-06).

## Self-Check: PASSED

- `spike/p13/alpha_series.jl` — FOUND (433 lines, `min_lines: 170` satisfied)
- `spike/test/test_p13_alpha.jl` — FOUND (344 lines)
- commit `9422620` — FOUND in `git log --all`
- commit `472b463` — FOUND in `git log --all`
- both `key_links` patterns present: `contract\.jl|_calculate_mask` in `alpha_series.jl`,
  `alpha_series\.jl` in `test_p13_alpha.jl`
