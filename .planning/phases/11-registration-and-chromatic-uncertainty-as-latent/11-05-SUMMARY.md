---
phase: 11-registration-and-chromatic-uncertainty-as-latent
plan: 05
subsystem: validation
tags: [julia, random123, philox, paired-design, ols-calibration, spearman, jld2, pre-registration, threading]

# Dependency graph
requires:
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 01
    provides: "spike/validation/p11_consts.jl Tier-1 — P11_PROBE_* spec, p11_rng/P11_PROBE_COUNTER, SC2_RUNGS, P11_IMSIZE_SET/WEIGHTS, P11_PROBE_S_FLOOR/SPAN_FLOOR, SC2_SPEARMAN_ATTENUATION, and the assertion that Tier 2 was not pre-empted"
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 02
    provides: "chromatic_eps as the 8th theta field and the widened SHIFT_PRIOR — the axes this probe sweeps"
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 03
    provides: "the single composed affine stage 6 (D-10/D-15), so misalignment is simulator-injected in ONE interpolation pass"
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 04
    provides: "spike/test/test_lambda_ablation.jl and its P11_LAMBDA_ABLATION_FACTOR placeholder wiring, which this plan supersedes"
provides:
  - "spike/validation/run_p11_probe.jl — the D-06 paired pre-flight probe (shift / chromatic / dRho axes, F5 mixture + 256-squared arm, produces thresholds, gates nothing)"
  - "spike/validation/p11_probe_report.jld2 — the single reported probe run, 2026-07-25T21:26:36.787Z, 7.109 min"
  - "MEASURED: S_probe = 1.0, dRho_eq ladder span = 0.046885 on the F5 mixture arm — the Tier-1 abort criterion did NOT fire"
  - "spike/validation/p11_consts.jl Tier-2 block — SC2_SPEARMAN_FLOOR, the dRho_eq calibration coefficients for both arms, P11_LADDER_IMSIZE_ARM, P11_LAMBDA_ABLATION_FACTOR, P11_PROBE_S_MEASURED, P11_PROBE_SPAN_MEASURED"
affects: [11-06-probe-verdict, 11-07-research-trainer, 11-08-ladder, 11-breakdown, 11-report]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Paired-key probing: a fresh RNG object re-derived from ONE key per replicate, so every rung shares stages 1-5 bit-for-bit and only the swept geometry differs"
    - "A produces-thresholds experiment script: house run_sbc.jl skeleton minus the test-set ending, plus a printed copy-paste-ready const block and no write path to the pre-registration"
    - "Named single-field theta setters instead of inline merge, because a missing trailing comma silently reaches merge's keyword method"
    - "Two-tier append-only pre-registration realised: a second guard block on a distinct sentinel, with the Tier-1 literals re-asserted in the same test run"

key-files:
  created:
    - spike/validation/run_p11_probe.jl
    - spike/test/test_p11_probe.jl
    - spike/validation/p11_probe_report.jld2
  modified:
    - spike/validation/p11_consts.jl
    - spike/test/test_p11_consts.jl
    - spike/test/test_lambda_ablation.jl

key-decisions:
  - "The shift arm is the UNION of P11_PROBE_SHIFT_RUNGS, P11_PROBE_SHIFT_EXT and SC2_RUNGS — S_probe is defined across SC2_RUNGS, so those magnitudes must be evaluated for the statistic to exist"
  - "The dRho calibration step is taken TOWARD the interior of [-1, 1] per theta base, so no rung clamps and flattens the curve (theta base 4 has rho_true = -0.99)"
  - "The bare main() call is guarded on P11_PROBE_LOAD_ONLY so the unit test can reach the helpers without re-consuming the reserved stream (T-11-21)"
  - "The probe threads over (theta, replicate) units; determinism is unaffected because every unit writes a fixed index and derives its own RNG from a key"
  - "P11_LAMBDA_ABLATION_FACTOR = 2.502 — materially stricter than the 1.15 placeholder, derived exactly as the plan specifies and recorded as such in advance"

patterns-established:
  - "A probe that supplies a threshold must prove it cannot write the file that threshold lands in — asserted by a source scan in the unit test, not only by convention"
  - "Pairing determinism is measured, not assumed: the zero rung is re-simulated rather than reused, so a non-bit-identical re-derivation shows up as a non-zero displacement"

requirements-completed: [D-04, D-06, D-07, D-08, D-09]

# Metrics
duration: ~50min
completed: 2026-07-25
---

# Phase 11 Plan 05: The D-06 Pre-Flight Probe and the Tier-2 Append Summary

**The cheapest risk control in the phase ran and came back clean: on the F5 mixture the 128-dim
summary moves strictly monotonically across all seven SC2 rungs (`S_probe = 1.0`) and the ladder is
worth `0.0469` in equivalent-`dRho` units against a pre-registered floor of `0.02`, so the "SC2 fails
flat / below resolution" branch is not the one this phase is on — and the SC2 threshold is now
DERIVED from that measurement in an append-only Tier-2 block rather than guessed.**

## ABORT-CRITERION: **NOT FIRED**

The Tier-1 criterion (frozen in `d336699`, before the probe existed) is:

> Declare "SC2 below resolution" if `S_probe < P11_PROBE_S_FLOOR` **OR**
> `dRho_eq(lambda_max) - dRho_eq(lambda_min) < P11_PROBE_SPAN_FLOOR`.

| Quantity | Measured (F5 mixture arm) | Tier-1 floor | Comparison |
|---|---|---|---|
| `S_probe` = `corspearman(SC2_RUNGS, mean dRho_eq)` | **1.0** | `P11_PROBE_S_FLOOR` = 0.9 | at or above |
| dRho_eq ladder span across `SC2_RUNGS` | **0.046884564903982365** | `P11_PROBE_SPAN_FLOOR` = 0.02 | at or above |

**Neither leg fired. Recorded, not acted on** — the branch decision belongs to plan 11-06, and the
`P11_ITERATION_ALLOWANCE = 1` remains unspent.

## Performance

- **Duration:** ~50 min wall (including the crashed first attempt and its fix)
- **Reported probe run:** 7.109 min, 32 threads, `generated = 2026-07-25T21:26:36.787Z`
- **Tasks:** 3 (plus one auto-fix commit)
- **Files created:** 3; **modified:** 3; **deleted:** 0

## The reported run, verbatim

Command: `julia --project=spike -t auto spike/validation/run_p11_probe.jl`
Stream: `p11_rng(P11_PROBE_COUNTER = 1)` off `P11_DEV_SEED = 0x000000000b11de71`.
Arms: F5 mixture `((512,512),(1024,1024),(1376,1028),(2048,2048))` w = `(0.40,0.25,0.25,0.10)`,
plus the `(256,256)` comparability arm. `theta bases = 5`, `R = 32`, 160 paired units of 24
evaluations per arm, `metric = paired_l2_rows_1_64`.

**Frozen theta base points** (geometry zeroed; nuisances held fixed for the whole run):

| # | rho_true | spillover | autofluorescence | label_efficiency | noise |
|---|---|---|---|---|---|
| 1 |  0.1797 | 0.0999 | 0.0772 | 0.8000 | 0.4546 |
| 2 | -0.0270 | 0.0815 | 0.0486 | 0.7159 | 0.4384 |
| 3 | -0.1926 | 0.1854 | 0.0504 | 0.9160 | 0.3464 |
| 4 | -0.9900 | 0.1423 | 0.0618 | 0.7760 | 0.5900 |
| 5 |  0.4656 | 0.1884 | 0.0678 | 0.6799 | 0.9174 |

### F5 mixture arm — SHIFT sweep (diagonal, dx = dy = |s|/sqrt(2))

| rung (px) | mean ‖ds‖ | median | sd | dRho_eq |
|---|---|---|---|---|
| 0.00 | 0.0 | 0.0 | 0.0 | -0.002807 |
| 0.25 | 0.09375 | 0.07003 | 0.0672 | 0.01171 |
| 0.50 | 0.09724 | 0.07463 | 0.0679 | 0.01225 |
| 1.00 | 0.1156 | 0.09156 | 0.0739 | 0.01509 |
| 1.50 | 0.2099 | 0.1705 | 0.133 | 0.02970 |
| 2.00 | 0.2441 | 0.1984 | 0.138 | 0.03499 |
| 2.50 | 0.2942 | 0.2472 | 0.145 | 0.04275 |
| 3.00 | 0.3966 | 0.3266 | 0.193 | 0.05859 |
| 4.00 | 0.5248 | 0.4921 | 0.221 | 0.07845 |
| 5.00 | 0.7063 | 0.7670 | 0.295 | 0.1066 |
| 6.00 | 0.9008 | 1.048 | 0.386 | 0.1367 |
| 8.00 | 1.284 | 1.267 | 0.605 | 0.1960 |

### F5 mixture arm — CHROMATIC sweep (shift = 0; SEPARATE AXIS)

| rung (chromatic_eps) | mean ‖ds‖ | median | sd | dRho_eq |
|---|---|---|---|---|
| 0.0000 | 0.0 | 0.0 | 0.0 | -0.002807 |
| 0.0025 | 0.09395 | 0.0923 | 0.0288 | 0.01174 |
| 0.0050 | 0.2172 | 0.1864 | 0.110 | 0.03082 |
| 0.0100 | 0.5513 | 0.3957 | 0.378 | 0.08256 |
| 0.0200 | 1.242 | 0.8061 | 0.877 | 0.1895 |
| 0.0300 | 1.744 | 1.375 | 1.17 | 0.2672 |
| 0.0500 | 2.294 | 1.870 | 1.41 | 0.3524 |

### F5 mixture arm — dRho CALIBRATION sweep (shift = 0, chromatic_eps = 0)

| dRho | mean ‖ds‖ | median | sd | dRho_eq |
|---|---|---|---|---|
| 0.02 | 0.1423 | 0.1378 | 0.0237 | 0.01923 |
| 0.05 | 0.3459 | 0.3339 | 0.0660 | 0.05074 |
| 0.10 | 0.6657 | 0.6630 | 0.0896 | 0.1003 |
| 0.20 | 1.308 | 1.304 | 0.142 | 0.1998 |

**Fit:** `‖ds‖_2 = 6.45851 * dRho + 0.0181265`, **R² = 0.999931**.

### 256-squared comparability arm — SHIFT sweep

| rung (px) | mean ‖ds‖ | median | sd | dRho_eq |
|---|---|---|---|---|
| 0.00 | 0.0 | 0.0 | 0.0 | -0.01367 |
| 0.25 | 0.3050 | 0.3002 | 0.0962 | 0.02800 |
| 0.50 | 0.3165 | 0.3167 | 0.0942 | 0.02957 |
| 1.00 | 0.3683 | 0.3853 | 0.103 | 0.03666 |
| 1.50 | 0.6291 | 0.6350 | 0.181 | 0.07230 |
| 2.00 | 0.7015 | 0.7212 | 0.186 | 0.08218 |
| 2.50 | 0.7940 | 0.8318 | 0.192 | 0.09482 |
| 3.00 | 1.031 | 1.046 | 0.241 | 0.1273 |
| 4.00 | 1.246 | 1.273 | 0.233 | 0.1566 |
| 5.00 | 1.581 | 1.551 | 0.253 | 0.2024 |
| 6.00 | 1.907 | 1.860 | 0.278 | 0.2469 |
| 8.00 | 2.471 | 2.446 | 0.376 | 0.3240 |

### 256-squared comparability arm — CHROMATIC sweep

| rung (chromatic_eps) | mean ‖ds‖ | median | sd | dRho_eq |
|---|---|---|---|---|
| 0.0000 | 0.0 | 0.0 | 0.0 | -0.01367 |
| 0.0025 | 0.1306 | 0.1352 | 0.0247 | 0.004174 |
| 0.0050 | 0.1742 | 0.1850 | 0.0398 | 0.01014 |
| 0.0100 | 0.2938 | 0.3188 | 0.0796 | 0.02648 |
| 0.0200 | 0.5773 | 0.6236 | 0.146 | 0.06521 |
| 0.0300 | 0.8884 | 0.9183 | 0.172 | 0.1077 |
| 0.0500 | 1.529 | 1.526 | 0.195 | 0.1953 |

### 256-squared comparability arm — dRho CALIBRATION sweep

| dRho | mean ‖ds‖ | median | sd | dRho_eq |
|---|---|---|---|---|
| 0.02 | 0.2374 | 0.2250 | 0.0514 | 0.01876 |
| 0.05 | 0.4781 | 0.4107 | 0.171 | 0.05166 |
| 0.10 | 0.8299 | 0.7754 | 0.218 | 0.09974 |
| 0.20 | 1.563 | 1.627 | 0.300 | 0.1998 |

**Fit:** `‖ds‖_2 = 7.31855 * dRho + 0.100025`, **R² = 0.999766**.

### Derived quantities

```
dRho_eq across SC2_RUNGS (F5 arm) : [0.01171, 0.01225, 0.01509, 0.0297, 0.03499, 0.04275, 0.05859]
dRho_eq across SC2_RUNGS (256 arm): [0.028,   0.02957, 0.03666, 0.0723, 0.08218, 0.09482, 0.1273]
lambda_max / lambda_min dRho_eq ratio (F5): 5.003943268109952
PAIRING DETERMINISM: max ||ds|| at the zero rungs = 0.0   (must be exactly 0.0)
elapsed_min = 7.11
```

### The printed copy-paste-ready const block, verbatim

```julia
    const P11_PROBE_S_MEASURED     = 1.0
    const P11_PROBE_SPAN_MEASURED  = 0.046884564903982365
    const SC2_SPEARMAN_FLOOR       = SC2_SPEARMAN_ATTENUATION * P11_PROBE_S_MEASURED
    const P11_DRHO_EQ_SLOPE        = 6.458505483772493
    const P11_DRHO_EQ_INTERCEPT    = 0.018126542853519556
    const P11_DRHO_EQ_R2           = 0.9999313512418424
    const P11_DRHO_EQ_SLOPE_256    = 7.318546621771334
    const P11_DRHO_EQ_INTERCEPT_256 = 0.10002486401477406
    const P11_DRHO_EQ_R2_256       = 0.9997658248160909
    const P11_LADDER_IMSIZE_ARM    = :f5_mixture
    const P11_LAMBDA_ABLATION_FACTOR = 2.501971634054976
    #   raw ratio dRho_eq(LAMBDA_MAX)/dRho_eq(LAMBDA_MIN) = 5.003943268109952
    #   attenuated: SC2_SPEARMAN_ATTENUATION * ratio      = 2.501971634054976
    #   floored at 1.05                                   = 2.501971634054976
```

## Task Commits

1. **Task 1: the paired pre-flight probe + its fixture-scale test** — `69ac6ae` (feat)
1a. **Auto-fix: named theta setters (Rule 1)** — `4c9e435` (fix)
2. **Task 2: the reported probe run and its artifact** — `66a48db` (test)
3. **Task 3: the Tier-2 append** — `2780a06` (chore)

## Accomplishments

- **The probe is genuinely paired, and the pairing is PROVEN rather than assumed.** Every rung of a
  replicate is simulated from a freshly re-derived RNG object on ONE key, so stages 1-5 are
  bit-identical and only the stage-6 geometry differs. The zero rung is deliberately **re-simulated**
  rather than reused, which turns the pairing into a measured invariant: `max ‖ds‖ at the zero rungs
  = 0.0`, exactly. The unit test additionally shows a DIFFERENT key on the same theta does **not**
  reproduce, so that zero is a property of the pairing and not of a degenerate simulator.
- **The threshold is derived, in the units D-06 asked for.** The calibration fit is near-perfect on
  both arms (R² = 0.99993 / 0.99977), so `dRho_eq` is a well-conditioned inversion rather than a
  strained one. "3 px of misalignment looks like a dRho of 0.059 on the realistic image mixture" is
  a directly quotable, reviewer-legible number.
- **Both axes are comfortably above resolution on the realistic image sizes**, which is the risk
  D-06 exists to retire. Registration: 0.0117 -> 0.0586 dRho_eq across the SC2 ladder (5.0x).
  Chromatic at its prior edge (0.02): 0.190 dRho_eq on the F5 mixture — roughly 3.2x what the
  ladder's widest rung is worth, and about 2.9x the same rung's response at 256².
- **The F5 mixture is measurably LESS sensitive than 256², as physics predicts.** At grid 8 a patch
  is 32 px at 256² but 256 px at 2048², so the same pixel shift is a far smaller relative
  displacement: dRho_eq(3 px) is 0.0586 on the mixture against 0.1273 at 256², a 2.2x attenuation.
  The indicative research probe was 256²-only, so it was optimistic by roughly that factor — this is
  precisely why ruling Q7 put the probe on the same mixture the net will train on.
- **The Tier-2 append is structurally honest.** Two guard blocks on distinct sentinels; additions
  only (`git diff HEAD~1` shows zero deletion lines in the consts file); a provenance comment per
  constant naming the artifact field it came from; and `SC2_SPEARMAN_FLOOR` written as a PRODUCT so
  the Tier-1 judgment allowance and the Tier-2 measurement can never be confused for one another.
- **The SC1g tripwire is no longer a placeholder.** `P11_LAMBDA_ABLATION_FACTOR` moves from the
  authored `1.15` to a probe-derived `2.502`, picked up with no edit to the tripwire's logic.

## Verification — observed results

Every command below was run and its real output observed.

| Check | Command | Observed |
|---|---|---|
| Task 1 unit gate | `julia --project=spike spike/test/test_p11_probe.jl` | **69 pass / 0 fail**, exit **0**, **9.16 s** wall |
| Task 2 artifact readback | the plan's `JLD2.load` verify one-liner | `S_probe=1.0  span=0.046884564903982365`; `generated=2026-07-25T21:26:36.787Z`, `elapsed_min=7.1087` |
| Task 2 governing consts recorded | the plan's 5-key `haskey` one-liner | prints `consts recorded` |
| Task 2 artifact is tracked | `git ls-files --error-unmatch spike/validation/p11_probe_report.jld2` | exit **0** |
| Task 2 consts byte-unchanged | `git diff --stat -- spike/validation/p11_consts.jl` (end of task) | empty |
| Task 3 consts gate | `julia --project=spike spike/test/test_p11_consts.jl` | **84 Tier-1 + 32 Tier-2 = 116 pass / 0 fail**, exit **0** |
| Task 3 downstream consumer | `julia --project=spike spike/test/test_p11_tost.jl` | exit **0** (36 pass, 1 deliberate skip) |
| Task 3 inline assertion | the plan's `include(...); @assert ...` one-liner | prints `tier2 ok` |
| Two guard blocks | `grep -c 'if !isdefined' spike/validation/p11_consts.jl` | `2` |
| Additions only | `git diff HEAD~1 -- spike/validation/p11_consts.jl \| grep -c '^-[^-]'` | `0` |
| Per-constant provenance | `grep -c 'p11_probe_report' spike/validation/p11_consts.jl` | `1` |
| Append notice present | `grep -c 'APPENDED, NEVER EDITED' spike/validation/p11_consts.jl` | `1` |
| Probe is not a gate | `grep -c '@testset' spike/validation/run_p11_probe.jl` | `0` |
| One prior draw only | `grep -c 'sample_prior' spike/validation/run_p11_probe.jl` | `1` (line 366, the frozen theta-base draw) |
| No write path to the consts file | `grep -cE '(open\|write\|jldsave)\(.*p11_consts' spike/validation/run_p11_probe.jl` | `0` |
| Consts include present | `grep -c 'p11_consts' spike/validation/run_p11_probe.jl` | `1` |
| Row restriction present | `grep -c '1:64' spike/validation/run_p11_probe.jl` | `2` |
| Knobs read, not retyped | `grep -c 'P11_PROBE_R\|P11_PROBE_SHIFT_RUNGS\|P11_PROBE_EPS_RUNGS\|P11_PROBE_DRHO_RUNGS' spike/validation/run_p11_probe.jl` | `16` |
| Atomic write | `grep -c 'mv(' spike/validation/run_p11_probe.jl` | `1` |
| Script length | `wc -l spike/validation/run_p11_probe.jl` | `559` (bar: >= 150) |
| Tier-2 read reaches the tripwire | direct eval of the tripwire's `isdefined` branch | `ABLATION_FACTOR=2.501971634054976  IS_TIER2=true` |
| Tripwire still refuses vacuously | `julia --project=spike spike/test/test_lambda_ablation.jl` | exit **1**, names `p11_research_npe.jld2` and plan 11-07 |
| No `src/` change | `git diff --stat 7cae358..HEAD -- src/` | empty |
| No dependency change | `git diff --stat 7cae358..HEAD -- spike/Project.toml spike/Manifest.toml` | empty |
| Files this plan touched | `git diff --name-only 7cae358..HEAD` | exactly the 6 declared files |
| No deletions in any commit | `git diff --diff-filter=D --name-only HEAD~1 HEAD` per commit | empty each time |

**NOT verified here, and deliberately not claimed.** Nothing in this plan trains, evaluates or
reads a net. **No calibration, coverage, posterior-width, monotonicity or identifiability claim is
made or supported by anything above.** In particular:

- `S_probe = 1.0` is a property of the **forward model's summary displacement**, not of any
  posterior. Whether posterior width actually tracks lambda is unknown and is measured for the
  first time in plan 11-07/11-08.
- `SC2_SPEARMAN_FLOOR = 0.5` is the bar the ladder must clear; nothing has been scored against it.
- The `dRho_eq` axis is a **reporting** transform. No pass/fail reads its coefficients.
- `P11_LAMBDA_ABLATION_FACTOR = 2.502` has **never been evaluated against a trained net**.
- The shift and chromatic theta columns remain **vacuous by design** (D-05). The probe measures
  DISPLACEMENT of the summary, never recoverability of `dx`/`dy`/`chromatic_eps`.
- `spike/test/runtests.jl` was not run to green — the pre-existing `SPEEDUP_GATE` blocker stands
  (see "Deferred Issues"), and it is out of this plan's scope.

## Decisions Made

1. **The shift arm is the UNION of three pre-registered rung sets.** `S_probe` is defined as
   `corspearman(SC2_RUNGS, mean dRho_eq per rung)`, so the shift arm must be evaluated at exactly
   the SC2 magnitudes or the statistic does not exist. `SC2_RUNGS` contributes 2.5 px, the probe
   ladder contributes the 0.0 px reference; the union adds no invented rung and the unit test
   asserts its length equals the union's cardinality.
2. **The dRho calibration step is taken toward the interior of `[-1, 1]`, once per theta base.**
   Theta base 4 drew `rho_true = -0.99`; a blind `+dRho` would have clamped at the 0.20 rung and
   silently flattened the calibration curve. The direction is fixed **before** any evaluation and
   is recorded in the artifact as `drho_direction`. The norm is a magnitude, so the sign of the
   step does not bias it — a clamped step would.
3. **`main()` is guarded on `P11_PROBE_LOAD_ONLY`, diverging from the house bare-`main()` shape.**
   The plan puts the probe's helpers in the script and its unit test in a separate file, so the
   test must `include` the script to reach them. A bare call would make the unit test re-run the
   reported probe on the reserved stream (T-11-21) and cost ~7 minutes per test invocation. The
   guard name exists nowhere else in the repository, so the reported invocation is otherwise
   unconditional.
4. **The probe threads over `(theta, replicate)` units.** Single-threaded the F5 arm was measured
   at ~0.86 s per simulate+summary call (512² 0.401 s, 1024² 0.643 s, 1376x1028 1.079 s, 2048²
   2.714 s), i.e. ~53 min for 3840 calls — past the plan's 30-minute anomaly bar. Threading is
   determinism-neutral here: each unit derives its own RNG from a key and writes a fixed index, so
   the result does not depend on scheduling. Measured: 6.74 min for the F5 arm.
5. **The imsize draw is keyed on the REPLICATE, on a second key family disjoint from the
   simulation keys.** One replicate is one image size held fixed across every rung, so a
   rung-to-rung difference can never be an image-size change. Disjointness is asserted at run time
   and in the unit test.
6. **Rows 65:128 are excluded and no z-scoring is applied.** The mask rows are near-constant and
   are exactly the rows `_summary_row_partition` excludes; all 64 retained rows are the same
   physical quantity on the same scale, and no fitted pool exists before training anyway.
7. **`P11_LAMBDA_ABLATION_FACTOR` is implemented exactly as the plan specifies** —
   `max(1.05, SC2_SPEARMAN_ATTENUATION * ratio)` — with the full arithmetic recorded in the
   constant's comment. See "Consequence to carry forward" below.

## Consequence to carry forward (stated in advance, not after seeing a result)

`P11_LAMBDA_ABLATION_FACTOR` moves from the authored placeholder `1.15` to a derived `2.502`. That
is a **materially harder bar**: plan 11-07's tripwire now demands that the posterior spread at
`lambda = 3.0` exceed the spread at `lambda = 0.25` by a factor of 2.5 on **two** spread statistics
across **three** datasets. This is what the plan directs (derive from the measured ratio, attenuate,
floor at 1.05) and what D-04 requires (fix it before any result exists). Two honest caveats belong
in the report:

- The bar is a **summary-displacement** ratio halved by an attenuation allowance. The quantity the
  tripwire measures is a **posterior-width** ratio, which passes through a learned, lossy stage the
  probe cannot see. The 0.5 allowance is the only thing standing between the two.
- The measured ratio 5.004 is if anything **conservative**: `dRho_eq(lambda = 0.25)` sits on the
  interpolation-onset plateau (0.0117, 0.0123, 0.0151 for 0.25, 0.5, 1.0 px), which inflates the
  denominator and shrinks the ratio.

If the tripwire fails at 2.502 in plan 11-07, that is a **result**, not a licence to relax the
constant — relaxing it would be exactly the Phase-7 amend-after-seeing pattern D-04 exists to
prevent.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] A one-field `merge` without a trailing comma reached `merge`'s keyword method and threw inside a worker thread**
- **Found during:** Task 2 (the first probe invocation)
- **Issue:** `merge(theta0, (chromatic_eps = e))` — no trailing comma — is not a NamedTuple merge.
  Julia parses `(chromatic_eps = e)` as an assignment expression, so the call reached
  `merge(::NamedTuple; kwargs...)` and raised
  `MethodError: no method matching merge(::@NamedTuple{...}; chromatic_eps::Float64)`. The failure
  surfaced only at run time, inside `Threads.@threads`, after the F5 arm had already started. The
  shift and dRho mutations were unaffected (two fields, and an explicit trailing comma).
- **Fix:** Introduced `_set_shift` / `_set_eps` / `_set_rho` as the **only** sanctioned theta
  mutations, each carrying the comma once, with an in-file comment naming the trap. Added a
  16-assertion testset proving each setter changes exactly the field it names, preserves the
  eight-field arity, and returns something `simulate_pair` accepts.
- **Files modified:** `spike/validation/run_p11_probe.jl`, `spike/test/test_p11_probe.jl`
- **Verification:** unit test 53 -> **69 pass**, exit 0; the reported run then completed cleanly.
- **Committed in:** `4c9e435`
- **Not a snooped run:** the crashed invocation produced **no numbers and no artifact**; it aborted
  before any displacement was computed. The single reported run is `66a48db`'s artifact, and the
  `P11_ITERATION_ALLOWANCE = 1` was not touched.

**2. [Rule 3 - Blocking] The house bare-`main()` ending had to be guarded**
- **Found during:** Task 1
- **Issue:** PATTERNS prescribes `main()` as a bare last line, and the plan puts the probe's
  helpers in the script while its unit test lives in a separate file. Including the script to reach
  the helpers would therefore execute the reported run — a ~7-minute test that also re-consumes the
  reserved stream (T-11-21).
- **Fix:** `isdefined(@__MODULE__, :P11_PROBE_LOAD_ONLY) || main()`, with an in-file paragraph
  naming the divergence, the reason, and the fact that no other file in the repository defines that
  name.
- **Files modified:** `spike/validation/run_p11_probe.jl`
- **Verification:** the unit test loads the helpers in 9.16 s without running the probe; the
  reported invocation ran unconditionally from the command line.
- **Committed in:** `69ac6ae`

**3. [Rule 2 - Missing critical] The probe's own no-write-path property is asserted, not merely conventional**
- **Found during:** Task 1
- **Issue:** T-11-19 (the probe writing its own thresholds) was covered by an acceptance-criterion
  `grep`, which is a one-time check that nothing re-runs.
- **Fix:** Added a testset that scans the probe's source and asserts (a) no `open`/`write`/`jldsave`
  call targets the consts file, (b) the file contains no test-set block, and (c) the prior sampler
  occurs at most once. Booleans are computed before `@test` so a failure prints `false` instead of
  dumping the whole source file (the 11-04 lesson).
- **Files modified:** `spike/test/test_p11_probe.jl`
- **Verification:** 4/4 in that testset.
- **Committed in:** `69ac6ae`, strengthened in `4c9e435`

**4. [Rule 2 - Missing critical] Threading, because the single-threaded F5 arm exceeded the plan's own anomaly bar**
- **Found during:** Task 1 (pre-run budgeting)
- **Issue:** Measured single-thread cost over the F5 mixture is ~0.86 s per simulate+summary call,
  so the pre-registered 3840-call arm would have taken ~53 min — past the plan's "stop it and record
  an anomaly" 30-minute threshold. RESEARCH §C11 budgets under 10 min **at 32 threads**, and the
  plan's own run line carries `-t auto`, so the threaded form is the intended one.
- **Fix:** `Threads.@threads` over `(theta, replicate)` units, each deriving its own RNG from a key
  and writing a fixed index, so the result is scheduling-independent. Progress logging via an
  atomic counter.
- **Files modified:** `spike/validation/run_p11_probe.jl`
- **Verification:** F5 arm 6.74 min, 256² arm 0.37 min, total 7.11 min — inside the RESEARCH budget.
- **Committed in:** `69ac6ae`

**5. [Rule 3 - Blocking] `LinearAlgebra` dropped in favour of an explicit norm**
- **Found during:** Task 1
- **Issue:** `spike/Project.toml` lists neither `LinearAlgebra` nor `Statistics`; relying on stdlib
  resolution for a **new** import was an avoidable resolve risk on a locked manifest (D-01 forbids
  adding a dependency).
- **Fix:** `sqrt(sum(abs2, ...))` instead of `norm(...)`; the import was removed. `Statistics` is
  retained because `run_sbc.jl` already depends on it resolving.
- **Files modified:** `spike/validation/run_p11_probe.jl`
- **Verification:** both the unit test and the reported run load cleanly; `spike/Project.toml` and
  `spike/Manifest.toml` are byte-unchanged.
- **Committed in:** `69ac6ae`

---

**Total deviations:** 5 auto-fixed (1 bug, 2 blocking, 2 missing-critical). **No Rule 4
(architectural) situation arose; no checkpoint was hit; no package was installed.**
**Impact on scope:** none. Deviation 1 is a real defect found by running the code. Deviations 2, 4
and 5 are execution mechanics with no effect on what is measured. Deviation 3 strengthens the
anti-snooping guarantee. No locked value was changed, no dependency added, no `src/` byte moved.

## Assumption Drift (advisory)

**1. "The indicative 256² probe numbers are representative of what the pre-registered probe will find."**
- **Found during:** Task 2
- **Planned:** RESEARCH §C9b's indicative table (dRho_eq 0.052 at 0.25 px, 0.191 at 3 px, a ~3.7x
  ladder headroom) framed the expected scale, and Tier-1's own band derivation quotes "a factor of
  order 3.7x".
- **Actual:** on the **F5 mixture** the ladder is 0.0117 -> 0.0586 (a 5.0x ratio but a *much smaller
  absolute* span, 0.047 rather than 0.139), because a fixed pixel shift is a far smaller relative
  displacement when a patch is 128-256 px instead of 32. The 256² arm reproduced the indicative
  regime as expected (0.028 -> 0.127).
- **Why it matters:** the ratio went **up** while the absolute span went **down** by ~3x. Anything
  in the phase that reasons about "how much dRho_eq the ladder is worth" must quote the mixture
  number (0.047), not the 256² one — including the D-08 band's "1/30 of the claimed effect"
  sentence, which was written against the 3.7x/256² framing.

**2. "The dRho_eq of a zero perturbation is zero."**
- **Found during:** Task 2
- **Planned:** implicit in reading `dRho_eq` as "the equivalent dRho".
- **Actual:** `dRho_eq` is an **affine** inversion with a fitted positive intercept, so a zero
  displacement maps to `-intercept/slope` = **-0.0028** (F5) and **-0.0137** (256²), not 0.
- **Why it matters:** a reader scanning a ladder table will see a small negative headline number in
  the zero row and could mistake it for a negative effect. It is an artifact of the inversion. The
  note is recorded in the Tier-2 block's comment so it travels with the coefficients.

**3. "The chromatic axis is comparable in magnitude to the registration axis, as D-09 designed."**
- **Found during:** Task 2
- **Planned:** D-09 chose the ±0.02 chromatic range explicitly to "match the ±3 px translation
  headroom so both misalignment axes are comparable in magnitude".
- **Actual:** on the F5 mixture they are **not** comparable: `chromatic_eps = 0.02` is worth 0.190
  dRho_eq while 3 px of shift is worth 0.059 — the chromatic term at its prior edge is ~3.2x the
  registration term at its prior edge. At 256² the two are much closer (0.065 vs 0.127, the other
  way round). The mechanism is exactly the one D-09 names: corner displacement scales with the
  half-diagonal, so the chromatic term grows with image size while a fixed pixel shift does not.
- **Why it matters:** the report should not repeat "both axes are comparable in magnitude" as a
  property of the realistic regime. It is true at 256² and false on the F5 mixture, and the
  direction of the asymmetry reverses.

## Issues Encountered

- **The first probe invocation crashed after ~1 minute** on the missing-comma `merge` (deviation 1).
  It produced no numbers and no artifact; diagnosed from the thrown `MethodError` and fixed before
  the reported run.
- **`stdout` is block-buffered when piped**, so the backgrounded run showed no interim progress
  until completion. Cosmetic; the progress lines are all present in the completed log.

## Deferred Issues

- **`spike/test/runtests.jl` still exits 1** on the pre-existing `SPEEDUP_GATE` (NPE-03) blocker
  recorded in `deferred-items.md` (92.5x / 84.0x against a `> 100x` bar). It is not this plan's to
  fix and the gate was **not** lowered. This plan added no file to the aggregate runner.

## Known Stubs

**None.** Every function in `run_p11_probe.jl` is fully implemented and exercised, and every Tier-2
constant carries a measured value from the committed artifact. The `1.15` fallback in
`spike/test/test_lambda_ablation.jl` is **no longer a placeholder** — it is retained only as the
standalone-load path (the `isdefined` guard was kept per the plan) and the Tier-2 value governs
whenever the pre-registration is in scope, which is the normal path. The file's prose was updated to
say exactly that.

## Threat Flags

None. This plan adds no network endpoint, no auth path, no schema at a trust boundary, and no file
write outside the probe's own report artifact. The registered threats were handled as planned:

| Threat ID | Handling | Evidence |
|---|---|---|
| T-11-19 (probe writing its own thresholds) | No write path exists; the values are PRINTED only | `grep -cE '(open\|write\|jldsave)\(.*p11_consts'` = `0`, plus a source-scan testset asserting the same |
| T-11-20 (Tier-1 silently edited during the append) | Additions-only append + Tier-1 literals re-asserted in the Tier-2 test run | `git diff HEAD~1 -- p11_consts.jl \| grep -c '^-[^-]'` = `0`; 11 Tier-1 literals re-asserted, all green |
| T-11-21 (probe re-run until convenient) | ONE reported run on the reserved counter; the crashed attempt produced nothing; the unit test is structurally prevented from re-running it | `elapsed_min`, `seed`, `salt`, `P11_PROBE_COUNTER` and a UTC stamp persisted in the artifact; `P11_PROBE_LOAD_ONLY` guard |
| T-11-22 (torn report artifact) | `.tmp` -> `jldsave` -> reopen + two `@assert haskey` -> `mv(...; force = true)`; persisted BEFORE any comparison | `grep -c 'mv('` = `1`; the artifact loads and every asserted key is present |
| T-11-23 (reading data outside the simulator) | Accepted by construction — no corpus, no real image, no network | the probe's only inputs are the pre-registration and the simulator |
| T-11-SC (package-manager installs) | None attempted | `git diff --stat 7cae358..HEAD -- spike/Project.toml spike/Manifest.toml` empty |

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- **Ready — plan 11-06 (the probe verdict).** Both raw numbers are frozen in the pre-registration
  (`P11_PROBE_S_MEASURED`, `P11_PROBE_SPAN_MEASURED`) so the verdict is reproducible from the file
  rather than from this document. The comparison is recorded above and, per the plan, **not acted
  on here**. `P11_ITERATION_ALLOWANCE = 1` is unspent.
- **Ready — plan 11-07 (research trainer).** Training is unblocked pending 11-06's verdict. Note the
  harder tripwire bar (2.502) described under "Consequence to carry forward"; the tripwire needs no
  edit, it already reads the Tier-2 constant.
- **For plan 11-08 (the SC2 ladder).** `SC2_SPEARMAN_FLOOR = 0.5` and the `dRho_eq` coefficients are
  frozen; the ladder's dRho_eq axis is the **F5 mixture** fit
  (`6.458505 * dRho + 0.018127`), not the 256² one.
- **For the report.** Assumption drifts 1-3 are all report-facing: the mixture-vs-256² span
  correction, the negative-zero-row artifact of the affine inversion, and the reversal of the
  "both axes comparable" claim on realistic image sizes.
- **Left red for someone else:** the pre-existing `SPEEDUP_GATE` blocker in
  `spike/test/runtests.jl`.

## Self-Check: PASSED

- All 3 created artifact paths exist on disk: `spike/validation/run_p11_probe.jl`,
  `spike/test/test_p11_probe.jl`, `spike/validation/p11_probe_report.jld2`.
- All 4 commits exist in `git log --oneline --all`: `69ac6ae`, `4c9e435`, `66a48db`, `2780a06`.
- `git diff --name-only 7cae358..HEAD` lists exactly the 6 files this plan declares — no `src/`
  file, no manifest, and no orchestrator-owned artifact (`STATE.md` / `ROADMAP.md` untouched).

---
*Phase: 11-registration-and-chromatic-uncertainty-as-latent*
*Completed: 2026-07-25*
