---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 05
subsystem: spike/p13 (D-06 summary-resolution probe)
tags: [tau, resolution-probe, unpaired, fit-free, pre-registration, D-06, D-04, D-01]
requires:
  - spike/p13/consts.jl (Tier-1 pre-registration; P13_TAU_* knobs, P13_IMSIZE_*, p13_rng/p13_fix_rng)
  - spike/simulator/prior.jl (sample_prior, MU_PRIOR, SHIFT_PRIOR, ghat knots)
  - spike/simulator/forward.jl (simulate_pair)
  - spike/contract.jl (build_mci, patch_summary - read-only src/)
  - spike/data/encode.jl (encode_d01, the 128-row layout)
  - spike/data/seeding.jl (HOLDOUT_SALT, the salt-XOR arm idiom)
  - spike/validation/ood.jl (the hand-rolled tie-aware roc_auc, REUSED)
provides:
  - mbar_from_summary (the mask-weighted D-06 statistic)
  - sample_imsize(rng, set, weights) (additive 3-arg method over an explicit mixture)
  - tau_probe_arm_keys (the two disjoint unpaired arm keys)
  - summary_draws (one arm, rho pinned, nuisances from their full priors)
  - tau_auc_at (A(delta) in both directions plus the reported magnitude and bootstrap median)
  - tau_curve (the WHOLE frozen grid, always)
  - tau_from_curve (the abort criterion; :measured / :below_resolution)
  - ghat_knot_spacing_near_zero (the sub-knot caution threshold, derived from the knots)
  - simulator_provenance_guard (reports the Phase-11 surgery markers; never throws)
affects:
  - 13-10 (runs the reported probe and appends the Tier-2 P13_TAU)
  - 13-08 (wires spike/test/test_p13_* into spike/test/runtests.jl)
tech-stack:
  added: []
  patterns:
    - guarded-include preamble (harness.jl idiom)
    - counter-based per-index Philox keys (seeding.jl idiom)
    - salt-XOR-into-the-first-key-word to carve a disjoint second arm
    - comment-stripped source-grep gate (test_p13_alpha.jl:123-124 idiom)
    - guarded `abspath(PROGRAM_FILE) == @__FILE__` script entry point
key-files:
  created:
    - spike/p13/tau_probe.jl
    - spike/test/test_p13_tau.jl
  modified: []
decisions:
  - "The two probe arms ride `reference = P13_DEV_SEED xor P13_SALT` and `contrast = reference xor HOLDOUT_SALT`. HOLDOUT_SALT was reused rather than a fresh salt invented, because carving a provably disjoint second stream is exactly what that constant already means in this repository, and inventing a stream-defining value in a non-pre-registration file would put it outside the frozen record."
  - "`sample_imsize` is extended with a THREE-POSITIONAL-ARGUMENT method rather than a keyword one, because `spike/data/seeding.jl:106` already owns `sample_imsize(rng::AbstractRNG)` and keyword arguments do not participate in Julia dispatch - the keyword form would have silently overwritten the generation sampler (or been overwritten by it, depending on include order)."
  - "`tau_curve` draws the rho = 0 reference arm ONCE and shares it across the grid (common random numbers). Each individual comparison stays unpaired; the sharing makes the curve's shape a statement about delta rather than about reference-arm noise, and costs one arm instead of seven."
  - "Each direction is reported as `max(A, 1 - A)`. Verified necessary: m-bar decreases with rho on the negative side, so the raw AUC of the rho = -delta arm ran to 0.0 at delta = 0.5 - perfect separation that a bar on the raw value would have scored as worse than chance."
  - "Non-finite m-bar draws (all-absent patch grid) are filtered and COUNTED in `n_degenerate` rather than passed to the AUC, because NaN compares false against everything and would silently contribute a zero to the U statistic."
metrics:
  duration_min: 34
  tasks_completed: 2
  files_created: 2
  files_modified: 0
  tests_added: 98
  completed: 2026-07-27
---

# Phase 13 Plan 05: D-06 Summary-Resolution Probe Summary

A fit-free, UNPAIRED, mask-weighted resolution probe that measures the smallest rho
difference the fixed 8x8 patch-correlation summary can reliably distinguish, needing no `zt`
and no trained network, reusing the in-repo tie-aware AUC, always reporting the whole frozen
delta grid, and returning `:below_resolution` rather than extending its own grid.

## What was built

**`spike/p13/tau_probe.jl` (657 lines)** — the probe machinery, flat top-level, no module
wrapper:

| Function | Role |
|---|---|
| `mbar_from_summary(Z128)` | the mask-weighted mean of the PRESENT continuous rows; `NaN` on an all-absent grid |
| `sample_imsize(rng, set, weights)` | inverse-CDF draw from an explicit categorical, one `rand` consumed |
| `tau_probe_arm_keys(; seed)` | the two disjoint Philox first-key-words that make the design unpaired |
| `summary_draws(rho, R, key; counter_base, imsize_set, imsize_weights)` | one arm: `rho_true` pinned, all 7 nuisances from their full priors |
| `tau_auc_at(delta; ...)` | `(neg, pos, auc = max, auc_median_boot, n_reference, n_neg, n_pos, n_degenerate)` |
| `tau_curve(; grid, R, bootstrap, ...)` | the WHOLE frozen grid, one row per delta |
| `tau_from_curve(curve; bar)` | `(tau, status, delta_index, bar, sub_knot)` — the abort criterion |
| `ghat_knot_spacing_near_zero()` | 0.0825, derived from `GHAT_RHO_KNOTS`, never retyped |
| `simulator_provenance_guard()` | `(post_p11, reason, has_chromatic_eps, theta_arity, shift_halfwidth, expected_halfwidth)` |

The banner states, in the house voice, what tau means; why the design is unpaired (with the
explicit contrast against Phase 11's deliberately *paired* displacement probe); why the
statistic is `m-bar` and not a learned readout; that the probe stops before the frozen summary
standardization so it needs no `zt`; and the DECOUPLING line. The script entry point is guarded
and prints only a usage note pointing at `spike/p13/run_tau_probe.jl` (plan 13-10) — including
the file simulates nothing.

**`spike/test/test_p13_tau.jl` (355 lines, 98 assertions)** — the fixture-scale gate: 10
testsets covering the fixture-stream rule, mask weighting, bitwise determinism under a fixed
key, the unpaired property, the discriminability magnitude, a labelled fixture-power sanity
guard, the abort criterion on synthetic curves plus a source grep against grid extension, the
frozen-knob references, the provenance guard's report-never-throw contract, and the CPU-only
postscript.

## Observed verification output (actually run, not inferred)

Every `<verify>` command was executed and its real output observed.

| Check | Observed |
|---|---|
| `julia --project=spike spike/test/test_p13_tau.jl` | **exit 0**, `98 Pass / 98 Total`, 0 fail, 0 error, **12.3 s** wall (budget 120 s) |
| `include("spike/p13/tau_probe.jl"); println("loaded")` | printed `loaded`, 8.4 s, no simulation |
| `mbar_from_summary(vcat([0.5,0.1,0.9,0.0],[1.0,1.0,1.0,0.0]))` | printed `mbar ok` (= 0.5; the absent 4th patch excluded, not averaged in as 0) |
| `tau_from_curve` synthetic-curve pair | printed `abort ok`; the sub-knot NOTE and the ABORT CRITERION FIRED block both printed as designed |
| `grep -v '^\s*#' tau_probe.jl \| grep -c 'standardize_summary'` | `0` |
| `grep -v '^\s*#' tau_probe.jl \| grep -c 'ROCAnalysis\|MLJ\|Pkg.add'` | `0` |
| `grep -c 'roc_auc' tau_probe.jl` / `grep -c 'function roc_auc'` (exec) | `3` / `0` — reused, never re-implemented |
| `grep -c 'P13_FIX_SEED' test_p13_tau.jl` | `4` |
| `grep -c 'not the reported measurement\|fixture power' test_p13_tau.jl` | `2` |
| `include("spike/test/test_p13_tau.jl"); @assert !isdefined(Main, :P13_TAU)` | printed `tier2 still absent`, exit 0 |
| `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` | **exit 0** |
| `git diff --name-only 318826a HEAD` | exactly `spike/p13/tau_probe.jl`, `spike/test/test_p13_tau.jl` |

`spike/p13/consts.jl`, `spike/test/runtests.jl`, `.planning/STATE.md`, `.planning/ROADMAP.md`,
`src/` and both spike manifests are byte-unchanged.

## Fixture-scale AUC values observed — NOT THE REPORTED MEASUREMENT

Measured off `P13_FIX_SEED` at `R = 16` per arm on a single 128x128 frame, `bootstrap = 64`:

| delta | `auc_neg` (magnitude) | `auc_pos` (magnitude) | reported `auc` = max | `auc_median_boot` |
|---|---|---|---|---|
| 0.02 | 0.789 | 0.715 | **0.789** | 0.803 |
| 0.50 | 1.000 | 0.996 | **1.000** | 1.000 |
| 0.60 | — | — | **1.000** | — |

`n_degenerate = 0` on every arm.

**These are NOT the reported measurement and must not be quoted as the summary's resolution.**
They were taken at `R = 16` (AUC standard error roughly 0.09, versus roughly 0.03 at the
frozen `P13_TAU_R = 400`), on ONE small 128x128 frame instead of the frozen F5 mixture (which
starts at 512x512 and has a 2048x2048 tail), and off the FIXTURE stream. Their only job is to
prove the code path computes something with the right monotone sense at fixture power. The
reported curve is plan 13-10's, at `P13_TAU_R` over `P13_IMSIZE_SET`, off `P13_DEV_SEED` at
`P13_TAU_COUNTER`.

**`P13_TAU` remains UNDEFINED.** No Tier-2 value was measured, computed, guessed or written.
`p13_tau()` still fails loudly, `spike/p13/consts.jl` is byte-unchanged, and the suite asserts
`!isdefined(:P13_TAU)`. The raw `A_neg = 0.0` observed at delta = 0.5 before the
`max(A, 1 - A)` magnitude step is recorded here because it is the empirical justification for
that step, not a result.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] The plan's `merge(..., (; rho_true = rho))` would have silently disabled the probe**

- **Found during:** Task 1
- **Issue:** `sample_prior` names the colocalization knob `ρ_true` (Unicode rho), not
  `rho_true` (`spike/simulator/prior.jl:96`). `merge` on a NamedTuple replaces an existing key
  in place but APPENDS an unknown one, so `merge(sample_prior(rng), (; rho_true = rho))` would
  have produced a NINE-field theta with `ρ_true` still at its prior draw. Both arms would then
  have been identically distributed, `A(delta)` would have come back at chance for every
  delta, and the probe would have reported `:below_resolution` for a reason with nothing to do
  with the summary's resolution — a silent, plausible-looking, phase-killing failure.
- **Fix:** pinned via `merge(sample_prior(rng), (; ρ_true = rho))`, and added two runtime
  assertions inside `summary_draws` (`length(theta) == 8` and `theta.ρ_true == rho`) so this
  class of bug fails loudly instead of flattening a curve. Asserted in the suite too.
- **Files modified:** `spike/p13/tau_probe.jl`, `spike/test/test_p13_tau.jl`
- **Commit:** `06f2f85`, `e6ca9f5`

**2. [Rule 3 - Blocking] `sample_imsize` keyword form would have clobbered the generation sampler**

- **Found during:** Task 1
- **Issue:** the plan specifies `sample_imsize(rng; set = P13_IMSIZE_SET, weights = P13_IMSIZE_WEIGHTS)`.
  `spike/data/seeding.jl:106` already defines `sample_imsize(rng::AbstractRNG)` over the
  generation-time `IMSIZE_SET`, and that file IS in scope here (the validation include chain
  pulls it). Keyword arguments do not participate in Julia dispatch, so the two definitions
  share one positional signature: whichever loaded second would silently overwrite the other.
  Either the probe would have drawn from the 256-squared-inclusive generation mixture, or
  `spike/data/generate.jl` would have started drawing from the F5 mixture — both silent.
- **Fix:** extended the generic with an ADDITIVE three-positional-argument method
  `sample_imsize(rng, set, weights)`; the F5 mixture remains the *default* at the
  `summary_draws` / `tau_auc_at` / `tau_curve` keyword level, so the frozen knobs are still the
  defaults and a fixture can substitute a small frame. Rationale recorded in the docstring.
- **Files modified:** `spike/p13/tau_probe.jl`
- **Commit:** `06f2f85`

**3. [Rule 2 - Missing critical functionality] NaN draws were not specified to be filtered**

- **Found during:** Task 1
- **Issue:** `mbar_from_summary` returns `NaN` on an all-absent grid (correct, per the plan).
  Handed straight to `roc_auc`, a `NaN` compares false against everything and contributes
  exactly 0 to the Mann-Whitney U — a silent downward bias on one arm with no trace in the
  output.
- **Fix:** `tau_auc_at` filters non-finite draws and REPORTS the count as `n_degenerate`;
  `summary_draws` still returns raw `NaN`s so the caller can count them rather than receive a
  silently shortened arm.
- **Files modified:** `spike/p13/tau_probe.jl`
- **Commit:** `06f2f85`

**4. [Rule 2 - Missing critical functionality] the two directional contrast arms shared a counter block**

- **Found during:** Task 1
- **Issue:** the plan carves the arms by KEY only. With one key and one counter range the
  `-delta` and `+delta` arms would have shared every nuisance draw, making the two directional
  AUCs — whose MAX is the reported statistic — needlessly correlated.
- **Fix:** the negative arm takes counters `1:R` and the positive arm `R+1:2R` on
  `arm_keys.contrast`, so no two arms anywhere in one comparison share a draw. Asserted
  (different counter block on the same key gives a different vector).
- **Files modified:** `spike/p13/tau_probe.jl`, `spike/test/test_p13_tau.jl`
- **Commit:** `06f2f85`, `e6ca9f5`

### Adaptations to the plan's literal text (no behaviour reopened)

- **`keys` keyword renamed `arm_keys`.** `keys` shadows `Base.keys`, which the CPU-only
  postscript itself calls (`keys(Base.loaded_modules)`). Same semantics.
- **Test calls pass fixture keywords explicitly.** The plan's testset text writes
  `tau_auc_at(0.5).auc`; at defaults that is `R = 400` over the frozen F5 mixture (512x512 to
  2048x2048) — minutes, and off the *reported* stream. Every simulating call in the suite goes
  through the `TAU_FIX` bundle instead, which is what "fixture scale is mandatory" requires.
  `tau_auc_at(0.0)` / `tau_auc_at(-0.1)` ARE called at defaults, because the `ArgumentError`
  fires before any simulation.
- **`tau_auc_at` returns a SUPERSET NamedTuple.** The plan's `(neg, pos, auc)` plus
  `auc_median_boot`, `n_reference`, `n_neg`, `n_pos`, `n_degenerate`. `tau_from_curve` reads
  only `delta` and `auc`, so the synthetic-curve acceptance case works unchanged.
- **`tau_auc_at` gained a `reference` keyword** so `tau_curve` draws the rho = 0 arm once
  instead of seven times. At the reported scale that is 2 400 fewer full simulate-and-summarize
  calls at F5 image sizes.

## Assumption Drift (advisory)

**The Phase-11 simulator surgery has ALREADY MERGED on this base.** Non-blocking; recorded
because it changes how two artefacts should be read.

- **Planned:** the plan (Task 2, item 8) says the suite must assert the guard "does not throw
  on today's *pre*-Phase-11 simulator", and 13-RESEARCH E4 frames the probe as blocking on a
  future simulator merge.
- **Actual:** at base `318826a`, `spike/simulator/prior.jl` already returns an 8-field theta
  including `chromatic_eps` (D-09) and `SHIFT_PRIOR = Uniform(-3.0, 3.0)` (D-02 widened). So
  `simulator_provenance_guard()` returns
  `post_p11 = true, reason = "POST-Phase-11 simulator confirmed: theta arity 8 with chromatic_eps present, shift-prior half-width 3.0 == expected 3.0."`
- **Why it matters:** the guard is written to be side-agnostic and the suite asserts only the
  report-never-throw contract plus, conditionally, the post-merge markers — so it stays green
  either way. But plan 13-10's Task 1 checkpoint should expect a `true` immediately rather than
  a block, and the fixture AUC values above were measured against the *widened* nuisance joint,
  not the easier pre-merge one.
- **Not a deviation:** nothing in this plan's scope changed, and no threshold moved.

## Known Stubs

None. Every function is fully implemented; nothing returns a placeholder, and no data source
is left unwired. `P13_TAU` is absent by *design* (Tier-2, plan 13-10) rather than stubbed —
`p13_tau()` fails loudly with the reason, and no code in this plan reads it.

## Threat Flags

None. The probe adds no network endpoint, no auth path, no schema change and no new file read
beyond the guarded includes the plan specifies. Threat-register dispositions all `mitigate`,
all implemented: T-13-17 (abort criterion + source grep + literal assertion), T-13-18
(`simulator_provenance_guard`), T-13-01 (fixture stream, `P13_FIX_SEED != P13_DEV_SEED`
asserted), T-13-07 (`ArgumentError` on non-positive delta / bad arity / bad rho; `NaN` on an
all-absent summary), T-13-19 (`roc_auc` reused, package greps clean), T-13-SC (no package
installed, manifests byte-unchanged).

## Handoff to plan 13-10

The seam is clean and unconsumed:

1. `simulator_provenance_guard()` returns `post_p11 = true` today — turn a `false` into the
   hard stop, quoting `.reason`.
2. `tau_curve()` at defaults IS the reported spec (`P13_TAU_DELTA_GRID`, `P13_TAU_R`,
   `P13_TAU_BOOTSTRAP_B`, `P13_TAU_BOTH_DIRECTIONS`, `P13_IMSIZE_SET`), off
   `P13_DEV_SEED xor P13_SALT` with the bootstrap keyed by `p13_rng(P13_TAU_COUNTER)`. No
   argument needs to be supplied to get the reported run.
3. **Budget warning, measured:** the reported run is 7 deltas x 2 directions x `R = 400` =
   5 600 contrast draws plus 400 reference draws, each a full simulate-and-summarize at an F5
   image size (512x512 to 2048x2048). The 128x128 fixture rate observed here was roughly
   0.3 s per draw; F5 frames are 16x to 256x the pixel count. Confirm the wall-clock estimate
   at the plan-13-10 Task-1 checkpoint before starting.
4. `tau_from_curve(curve)` gives the verdict; `.sub_knot` is the Pitfall-8 caution flag and
   `ghat_knot_spacing_near_zero()` (= 0.0825) is the threshold it uses.
5. Persist the whole curve BEFORE evaluating the verdict; every row already carries
   `auc_neg`, `auc_pos`, `auc`, `auc_median_boot`, `n_reference`, `n_neg`, `n_pos` and
   `n_degenerate`.

`spike/test/test_p13_tau.jl` is NOT yet wired into `spike/test/runtests.jl` — plan 13-08 owns
that, and `runtests.jl` was not touched here.

## Self-Check: PASSED

- `spike/p13/tau_probe.jl` — FOUND (657 lines)
- `spike/test/test_p13_tau.jl` — FOUND (355 lines)
- commit `06f2f85` — FOUND (`feat(13-05): add the fit-free unpaired D-06 summary-resolution probe`)
- commit `e6ca9f5` — FOUND (`test(13-05): gate the D-06 tau probe at fixture scale`)
- `git diff --name-only 318826a HEAD` lists exactly the two planned files — CONFIRMED
- `spike/p13/consts.jl`, `spike/test/runtests.jl`, `src/`, `spike/Project.toml`,
  `spike/Manifest.toml`, `.planning/STATE.md`, `.planning/ROADMAP.md` — byte-unchanged
