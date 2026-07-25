# Phase 11 — D-06 Pre-Flight Probe, Verdict

**Verdict date:** 2026-07-25
**Probe artifact:** `spike/validation/p11_probe_report.jld2`, `generated = 2026-07-25T21:26:36.787Z`, `elapsed_min = 7.1087`
**Pre-registration:** `spike/validation/p11_consts.jl` — Tier 1 frozen at commit `d336699` (before the probe existed); Tier 2 appended at `2780a06` (measurements only, additions only)
**Reserved stream:** `p11_rng(P11_PROBE_COUNTER = 1)` off `P11_DEV_SEED = 0x000000000b11de71`, salt `P11_SALT = 0xA24BAED4663EE121`
**Arms:** F5 mixture `((512,512),(1024,1024),(1376,1028),(2048,2048))`, `w = (0.40, 0.25, 0.25, 0.10)` — the arm the ladder trains and evaluates on — plus the `(256,256)` comparability arm
**Metric:** `P11_PROBE_METRIC = :paired_l2_rows_1_64` (fixed in Tier 1, not chosen after the fact)
**Protocol:** one probe run, criterion evaluated against constants frozen before the run. No re-run, no re-seed, no threshold change.

---

## Verdict: **ABOVE RESOLUTION**

The Tier-1 abort criterion **did NOT fire**. Both legs clear their pre-registered floor.

| Quantity | Measured (F5 mixture arm) | Pre-registered threshold | Leg fires? |
|---|---|---|---|
| `S_probe` = `corspearman(SC2_RUNGS, mean Δρ_eq)` | **1.0** | needs >= 0.9 (`P11_PROBE_S_FLOOR`) | no |
| ladder span (`Δρ_eq`) across `SC2_RUNGS` | **0.046884564903982** | needs >= 0.02 (`P11_PROBE_SPAN_FLOOR`) | no |

The criterion, as frozen in Tier 1 §6:

> Declare "SC2 below resolution" — do NOT train — if either holds:
> (i) `S_probe < P11_PROBE_S_FLOOR`, or (ii) `Δρ_eq(λ_max) − Δρ_eq(λ_min) < P11_PROBE_SPAN_FLOOR`.

Evaluated mechanically from the frozen file rather than transcribed:

```
$ julia --project=spike -e 'include("spike/validation/p11_consts.jl");
    fired = (P11_PROBE_S_MEASURED < P11_PROBE_S_FLOOR) ||
            (P11_PROBE_SPAN_MEASURED < P11_PROBE_SPAN_FLOOR);
    println("abort_criterion_fired=", fired)'
abort_criterion_fired=false
```

### Independent re-derivation from the raw artifact

The two Tier-2 numbers are not taken on trust from the 11-05 SUMMARY or from the Tier-2 block.
Both were **recomputed from the raw per-`(θ, replicate, rung)` tables** stored in
`p11_probe_report.jld2` (`shift_norm_f5`, dims `(5, 32, 12)`; `drho_norm_f5`, dims `(5, 32, 4)`),
using per-rung means, an OLS fit, a `Δρ_eq` inversion, midrank Spearman and a span, all written
independently of the probe script's own helpers. No simulation, no RNG, no stream consumed.

| Quantity | Re-derived from raw tables | Tier-2 frozen constant | Artifact headline field |
|---|---|---|---|
| `S_probe` | 1.0 | `P11_PROBE_S_MEASURED` = 1.0 | `S_probe` = 1.0 |
| ladder span | 0.046884564903982295 | `P11_PROBE_SPAN_MEASURED` = 0.046884564903982365 | `ladder_span` = 0.046884564903982365 |
| `Δρ_eq` slope (F5) | 6.4585054837724964 | `P11_DRHO_EQ_SLOPE` = 6.458505483772493 | — |
| `Δρ_eq` intercept (F5) | 0.018126542853519334 | `P11_DRHO_EQ_INTERCEPT` = 0.018126542853519556 | — |
| `λ_max/λ_min` ratio | 5.003943268109935 | (Tier 2 records 5.003943268109952) | `lambda_ratio` |

The residual differences are at the `1e-16` level and are float summation-order artifacts of a
differently written reduction; they are five orders of magnitude below the `1e-12` agreement
tolerance and cannot move either leg of the criterion. **Re-derived verdict = frozen verdict =
NOT FIRED.**

---

## What the probe measured

All numbers below are reproduced from `11-05-SUMMARY.md` and the committed artifact. `Δρ_eq` is a
**reporting** transform — the calibrated conversion of a summary displacement into "the Δρ that
would have moved the summary this far". No pass/fail reads its coefficients.

### F5 mixture arm — SHIFT sweep (diagonal, `dx = dy = |s|/√2`)

| rung (px) | mean ‖Δs‖ | median | sd | Δρ_eq |
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

**The 0 → 0.25 px step is interpolation onset, not misalignment sensitivity.** At exactly integer
offsets the backward map samples on-grid and `BSpline(Linear())` is an identity lookup; any
sub-pixel offset engages the interpolation kernel. The 0 → 0.25 jump therefore mixes a
qualitatively different (no-interpolation) regime with the misalignment response, and it is
precisely why the SC2 ladder starts at `LAMBDA_MIN = 0.25` and **carries no λ = 0 rung**. It is
also why the Tier-1 criterion is defined on the **ladder span** rather than on every adjacent
pair: the 0.25 → 0.50 → 1.00 stretch (0.01171, 0.01225, 0.01509) is a known local plateau on that
same onset, and an adjacent-pair criterion would have failed the design on a simulator artifact.

A second reading note, recorded so it travels with the tables: because the calibration fit has a
positive intercept, a **zero** displacement inverts to a slightly **negative** `Δρ_eq`
(`−intercept/slope` = −0.0028 on F5, −0.0137 at 256²). That is an artifact of the affine
inversion, not a negative effect.

### F5 mixture arm — CHROMATIC sweep (shift = 0; SEPARATE AXIS)

| rung (`chromatic_eps`) | mean ‖Δs‖ | median | sd | Δρ_eq |
|---|---|---|---|---|
| 0.0000 | 0.0 | 0.0 | 0.0 | -0.002807 |
| 0.0025 | 0.09395 | 0.0923 | 0.0288 | 0.01174 |
| 0.0050 | 0.2172 | 0.1864 | 0.110 | 0.03082 |
| 0.0100 | 0.5513 | 0.3957 | 0.378 | 0.08256 |
| 0.0200 | 1.242 | 0.8061 | 0.877 | 0.1895 |
| 0.0300 | 1.744 | 1.375 | 1.17 | 0.2672 |
| 0.0500 | 2.294 | 1.870 | 1.41 | 0.3524 |

### F5 mixture arm — Δρ CALIBRATION sweep (shift = 0, `chromatic_eps` = 0)

| Δρ | mean ‖Δs‖ | median | sd | Δρ_eq |
|---|---|---|---|---|
| 0.02 | 0.1423 | 0.1378 | 0.0237 | 0.01923 |
| 0.05 | 0.3459 | 0.3339 | 0.0660 | 0.05074 |
| 0.10 | 0.6657 | 0.6630 | 0.0896 | 0.1003 |
| 0.20 | 1.308 | 1.304 | 0.142 | 0.1998 |

**Fit:** `‖Δs‖₂ = 6.458505 · Δρ + 0.018127`, **R² = 0.999931**.

### 256² comparability arm — SHIFT sweep

| rung (px) | mean ‖Δs‖ | median | sd | Δρ_eq |
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

### 256² comparability arm — CHROMATIC sweep

| rung (`chromatic_eps`) | mean ‖Δs‖ | median | sd | Δρ_eq |
|---|---|---|---|---|
| 0.0000 | 0.0 | 0.0 | 0.0 | -0.01367 |
| 0.0025 | 0.1306 | 0.1352 | 0.0247 | 0.004174 |
| 0.0050 | 0.1742 | 0.1850 | 0.0398 | 0.01014 |
| 0.0100 | 0.2938 | 0.3188 | 0.0796 | 0.02648 |
| 0.0200 | 0.5773 | 0.6236 | 0.146 | 0.06521 |
| 0.0300 | 0.8884 | 0.9183 | 0.172 | 0.1077 |
| 0.0500 | 1.529 | 1.526 | 0.195 | 0.1953 |

### 256² comparability arm — Δρ CALIBRATION sweep

| Δρ | mean ‖Δs‖ | median | sd | Δρ_eq |
|---|---|---|---|---|
| 0.02 | 0.2374 | 0.2250 | 0.0514 | 0.01876 |
| 0.05 | 0.4781 | 0.4107 | 0.171 | 0.05166 |
| 0.10 | 0.8299 | 0.7754 | 0.218 | 0.09974 |
| 0.20 | 1.563 | 1.627 | 0.300 | 0.1998 |

**Fit:** `‖Δs‖₂ = 7.318547 · Δρ + 0.100025`, **R² = 0.999766**.

### Three properties of the measurement worth carrying forward

1. **The pairing is proven, not assumed.** `pairing_zero_max = 0.0` exactly: the zero rung was
   re-simulated from a freshly re-derived RNG object rather than reused, so a non-bit-identical
   re-derivation would have shown as a non-zero displacement. Every rung of a replicate shares
   stages 1-5 bit-for-bit; only the stage-6 geometry differs.
2. **The mixture is measurably LESS sensitive than 256², as physics predicts.** At grid 8 a patch
   is 32 px at 256² but 256 px at 2048², so the same pixel shift is a far smaller relative
   displacement: `Δρ_eq(3 px)` is 0.0586 on the mixture against 0.1273 at 256², a 2.2× attenuation.
   Anything downstream that quotes "how much Δρ_eq the ladder is worth" must quote the **mixture**
   number (0.047 span), not the 256² one.
3. **The two misalignment axes are NOT comparable in magnitude on the realistic arm.** D-09 chose
   the ±0.02 chromatic range to match the ±3 px translation headroom, and at 256² they are close
   (0.065 vs 0.127). On the F5 mixture the relation reverses and widens: `chromatic_eps = 0.02` is
   worth 0.1895 Δρ_eq against 0.0586 for 3 px of shift, roughly 3.2×. Corner displacement scales
   with the half-diagonal, so the chromatic term grows with image size while a fixed pixel shift
   does not.

---

## Scope of the SC2 claim (ruling Q3)

**The ladder is a registration ladder only.** The three probe axes were swept separately on
purpose, and the SC2 monotone-widening claim covers the **registration** axis and nothing else:

- λ is defined as the half-width of the **shift** prior, and it is λ — not `chromatic_eps` — that
  is appended to the network input as the 129th row (D-03). The chromatic term is **not λ-scaled**
  and **carries no conditioning input**, so there is no chromatic ladder to be monotone in.
- The chromatic dose-response therefore appears only on its own fixed-injection curve, in the
  SC3 / D-14 breakdown arm, never inside the SC2 ladder.
- A reader must **not** take SC2 as covering chromatic misalignment. "Posterior width increases
  monotonically with declared uncertainty" is a statement about declared **registration**
  uncertainty.

Two further scope limits, restated so the verdict cannot be over-read:

- `S_probe = 1.0` is a property of the **forward model's summary displacement**, not of any
  posterior. Whether posterior width actually tracks λ is unknown and is measured for the first
  time in plans 11-07/11-08. Passing the probe buys the right to try; it does not predict SC2.
- The shift and `chromatic_eps` θ columns remain **vacuous by design** (D-05). The probe measures
  displacement of the summary, never recoverability of `dx`/`dy`/`chromatic_eps`. Nothing here is
  an identifiability claim.

---

## The three pre-registered branches

Reproduced from RESEARCH §C10(c) / the plan context. There is no fourth branch.

1. **PROCEED.** The criterion did not fire. Training is authorised as planned.
2. **SPEND THE ITERATION.** The criterion fired. Do **not** train. Re-parameterise the *ladder*,
   not the *criterion*: raise the maximum toward the D-14 extension range (up to 8 px) and/or move
   the ladder to the F5 realistic image sizes where the chromatic response is much larger. Re-run
   the probe **once**. This consumes `P11_ITERATION_ALLOWANCE`.
3. **STOP SC2.** The criterion fired again after the single re-run, or the developer chooses not
   to spend the allowance. Report SC2 as **NOT DEMONSTRABLE at the fixed 8x8 summary resolution**,
   with the probe curve as the evidence; route the remaining budget to SC3 (a coverage claim,
   which does not require a resolvable width gradient) plus the D-16 docs note; and link the
   finding explicitly to the deferred "summary extension for registration identifiability" idea,
   since a below-resolution probe is direct evidence *for* that deferred work.

### Concrete consequence of each branch

| Branch | Plans that run | Plans skipped or cancelled | What the report claims |
|---|---|---|---|
| 1 — PROCEED | 11-07 … 11-11 as planned | none | The full SC2 ladder, the SC3 breakdown curve, the D-02 attenuation A/B, the D-17 real-image check, and the D-16 docs note |
| 2 — SPEND THE ITERATION | 11-05 Task 2 re-opened once with a re-parameterised ladder; return to this checkpoint | 11-07 BLOCKED until the re-run verdict | Nothing yet; the re-measured Tier-2 values are appended as a **THIRD** guard block, never editing the second |
| 3 — STOP SC2 | 11-11 only | 11-07, 11-08, the SC2 portion of 11-11; and 11-09's SC3 arm and 11-10's real-image arm, which cannot run without a research net | SC2 NOT DEMONSTRABLE at the fixed 8x8 summary resolution, with the probe curve as evidence; SC3's coverage claim not reachable without the research net; the finding linked to the deferred summary-extension idea |

Note that branch 3 is not free of collateral: SC3 and the real-image check both need the research
net that branch 3 declines to train, so choosing it concedes more than the SC2 headline alone.

---

## Iteration ledger

| Item | Value |
|---|---|
| `P11_ITERATION_ALLOWANCE` (declared in Tier 1, before any result existed) | **1** |
| Iterations spent so far | **0** |
| Iterations remaining | **1** |

A **second** iteration is not authorised by the pre-registration file, and — per Tier 1 §10 — it
cannot be authorised by amending that file. The single allowance is cheapest spent on a
re-parameterised probe *before* training, never on relaxing a threshold *after* it.

For the record, the one crashed probe invocation during plan 11-05 (a missing-comma `merge` that
threw inside a worker thread) produced **no numbers and no artifact** — it aborted before any
displacement was computed — so it is not an iteration and did not touch the allowance. The single
reported run is the one stamped `2026-07-25T21:26:36.787Z`.

DECISION: pending

---

## What is forbidden under every branch

- **Widening the ladder a second time.** One re-parameterisation is authorised; a second is the
  Phase-7 amend-twice pattern D-04 exists to avoid.
- **Relaxing `SC2_SPEARMAN_FLOOR`** (or `P11_PROBE_S_FLOOR`, `P11_PROBE_SPAN_FLOOR`, or
  `P11_LAMBDA_ABLATION_FACTOR`) after seeing a result it failed. `spike/test/test_p11_consts.jl`
  re-asserts the Tier-1 literals, so a floor edit breaks the suite loudly.
- **Re-running the probe on a different seed** to obtain a nicer verdict. The reserved stream is
  `p11_rng(P11_PROBE_COUNTER = 1)`; the counter is pre-registered and the burned research keys are
  enumerated as forbidden.
- **Editing any constant in either tier.** Constants are only ever appended, so the git history of
  `p11_consts.jl` is itself the audit trail; a diff that *modifies* a line is by construction a
  pre-registration breach.

---

*Phase: 11-registration-and-chromatic-uncertainty-as-latent*
*Plan: 11-06*
