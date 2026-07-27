---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 07
subsystem: amortized-evidence / verification
tags: [D-07, D-11, D-01, D-04, correction, negative-control, closed-form, honest-finding]
requires:
  - spike/p13/consts.jl        # frozen Tier-1 pre-registration (P13_F5_* bars, fixture stream)
  - spike/p13/labels.jl        # ThreeWayClass, target_matrix, head_log_odds
  - spike/p13/net.jl           # ThreeWayEvidenceNet, masked_two_head_bce, train_three_way
provides:
  - toy_log_evidence
  - toy_prior_class_masses
  - toy_analytic_log_bf
  - toy_class_of
  - toy_draw_dataset
  - toy_reshaped_dataset
  - toy_ideal_logits
  - toy_fit_two_head
  - toy_head_grid
  - toy_compare
affects:
  - spike/test/runtests.jl     # NOT edited here; plan 13-08 wires this file in
tech-stack:
  added: []                    # no package installed, no dependency added
  patterns:
    - closed-form conjugate-Gaussian ground truth instead of quadrature or KDE
    - corrected/uncorrected as ONE parameter of ONE comparison function
    - Bayes-risk reference logits to separate "wrong algebra" from "under-trained"
key-files:
  created:
    - spike/p13/toy_gaussian.jl          # 545 lines
    - spike/test/test_p13_correction.jl  # 454 lines
  modified: []
decisions:
  - "The D-07 per-head correction is CONFIRMED against analytic truth on the offset it exists to remove (mean deviation 0.069 / 0.067 nat corrected, versus 0.568 / 0.977 uncorrected)."
  - "The pre-registered P13_F5_MAXABS_TOL bar is MISSED on both heads (0.5498 / 0.3783 vs 0.25). Recorded as an honest finding, not repaired."
  - "The miss is NOT under-training: the trained net sits 0.0018 nat above the Bayes risk of the same problem."
  - "Each head is scored on its OWN restricted class pair (D-11), matching the convention the frozen consts already commit to for the reported gate."
metrics:
  duration: ~75 min
  completed: 2026-07-27
  tasks: 2
  tests: 90 assertions (88 pass, 2 fail -- both are the pre-registered max-deviation bar)
  runtime: 40.8 s wall for the whole file (bar: under 120 s)
---

# Phase 13 Plan 07: D-07 Three-Way Evidence-Scale Correction Verification Summary

A closed-form conjugate-Gaussian toy gives the three-way log Bayes factors analytically, and the
D-07 per-head correction is verified against them: the correction removes the evidence-scale
offset to within 0.07 nat and both mandatory negative controls behave exactly as the derivation
predicts, but the pre-registered `max |Δ log BF| ≤ 0.25` bar is **missed on both heads** because
logit-BCE is nearly flat in the confident tail — a shortfall recorded honestly rather than tuned
away.

## VERDICT: PARTIAL — one pre-registered bar PASSED, one MISSED, both negative controls PASSED

`julia --project=spike spike/test/test_p13_correction.jl` exits **1**: 88 assertions pass, 2 fail.
The two failures are exactly the two `maxabs <= P13_F5_MAXABS_TOL` assertions. Nothing else in the
file fails.

### The required four-row table (observed, D-11 restricted per-head grids)

Bars, read from the frozen `spike/p13/consts.jl`: `P13_F5_CORR_MIN = 0.99`,
`P13_F5_MAXABS_TOL = 0.25`, over the central `P13_F5_CENTRAL_FRAC = 0.90` of the evaluation grid.

| Arm | Head | `corr` (bar ≥ 0.99) | `maxabs` (bar ≤ 0.25) | `meandev` | `maxabs_demeaned` | Verdict |
|---|---|---|---|---|---|---|
| CORRECTED | coloc | **0.99751** ✅ | **0.5498** ❌ | −0.0693 | 0.4805 | corr PASS, max MISS |
| CORRECTED | exclusion | **0.99880** ✅ | **0.3783** ❌ | +0.0668 | 0.4452 | corr PASS, max MISS |
| UNCORRECTED (control 1) | coloc | 0.99751 | **1.0482** (must fail) ✅ | **−0.5676** | — | control PASSED |
| UNCORRECTED (control 1) | exclusion | 0.99880 | **1.1726** (must fail) ✅ | **+0.9773** | — | control PASSED |
| RESHAPED (control 2) | exclusion | 0.99831 | **1.0289** (must fail) ✅ | −0.7873 | **0.6745** | control PASSED |
| RESHAPED (control 2) | coloc *(classes untouched)* | 0.99873 | 0.3607 | −0.0535 | 0.3071 | unaffected, as predicted |

Measured per-head corrections on the training subset:
`head_log_odds = (coloc = −0.49833, exclusion = +0.91050)`.
Evaluation bands: coloc `Z ∈ [−0.756, 1.765]` (n = 1440), exclusion `Z ∈ [−2.188, 0.504]` (n = 3060).

## What was actually established

**1. The correction is right.** This is the substantive result and it is unambiguous. Without the
correction the two heads are off by a systematic **−0.568** and **+0.977** nat. With it they are off
by **−0.069** and **+0.067** nat. The residual offsets sit inside `P13_F5_MAXABS_TOL`, and the test
asserts that separately from the max (`abs(meandev) <= P13_F5_MAXABS_TOL`, both heads, PASS).

**2. Negative control 1 fires, and fires with the specific expected magnitude.** The uncorrected
signed mean deviation lands on the head's own measured log-odds to within 0.07 nat on both heads
(asserted at `atol = 0.15`), and the identity `uncorrected − corrected == log(q_pos/q_neg)` holds to
1e-9. This is what rules out the three failure modes of 13-RESEARCH Pitfall 1 — wrong sign, wrong
denominator, double application — each of which would also "fail" the bars but would not land on
that number.

The control also produced a finding worth carrying forward: **correlation cannot detect a missing
correction at all.** A missing correction is a constant shift and correlation is shift-invariant, so
the uncorrected `corr` is bit-comparable to the corrected one (0.99751 both). A verification that
gated on `P13_F5_CORR_MIN` alone would have passed a completely uncorrected net. This is asserted in
the file so no future reader adds a correlation-only gate.

**3. Negative control 2 fires, and fires irreparably.** Under a within-class reshaped proposal at
*identical class frequencies*, the exclusion head's deviation grows to 1.0289 while the coloc head —
whose own two classes were untouched — stays at 0.3607, i.e. essentially its un-reshaped value. And
after removing **the best constant in hindsight**, the exclusion head's residual is still **0.6745**,
far above the tolerance. That is 13-RESEARCH F2's theorem as a number: the discrepancy is a function
of `Z`, and no scalar — not the measured one, not a better-measured one, not an oracle one — could
have repaired it. `P13_STRATIFICATION = :class_frequency` is therefore evidence-backed rather than
asserted, and `P13_STRATIFICATION_FALLBACK = :importance_weighted` is the right shape of fallback
(a per-sample weight, not a different constant).

**4. The correction uses the head's own denominator.** `head_log_odds.coloc == log(n_C/n_R)` to
1e-12, and is asserted *not* to equal `log(n_C/(n_R+n_E))` (the one-vs-rest form that would silently
break D-08) nor `log(n_C/(n_C+n_R+n_E))` (the whole-set form the shipped binary read surface uses).
It is measured on the training subset only, and the test asserts the training and full-dataset
counts actually differ so that assertion is not vacuous.

## The honest finding: why the max-deviation bar is missed

**The miss is not under-training, and this is measured, not argued.** The test scores the analytic
Bayes-optimal logits (`log BF + head_log_odds`, the exact inverse of the D-07 correction) under the
same masked loss on the same held-out data:

```
held-out masked BCE:  net = 0.408983   Bayes-optimal = 0.407191   excess = 0.001792
```

The trained toy net is **0.0018 nat (0.4 %) above the Bayes risk of its own classification problem**.
It has extracted essentially all the information the objective can express.

**The mechanism is confident-tail flatness of logit-BCE.** Both heads' worst deviations sit at the
outer edge of the retained central band, where the true log Bayes factor is large:

| Head | worst-deviation location | true log BF | predicted | deviation |
|---|---|---|---|---|
| coloc | `Z = 1.765` (band edge, `hi = 1.765`) | +4.564 | +4.014 | −0.550 |
| exclusion | `Z = −2.188` (band edge, `lo = −2.188`) | +6.953 | +6.574 | −0.378 |

Where the head is already right with probability near one, a half-nat logit error costs almost
nothing in loss — so the optimizer has no gradient pressure to place it correctly, and a net at the
Bayes risk is still visibly off in nats out there. The bulk is excellent; only the extreme tail of
the band misses. Deviation quantiles measured in a scratch diagnostic under the *identical*
configuration (not printed by the committed test, stated here as a characterization):

| Head | \|dev\| q50 | q90 | q99 | max |
|---|---|---|---|---|
| coloc | 0.103 | 0.243 | 0.397 | 0.550 |
| exclusion | 0.050 | 0.248 | 0.262 | 0.378 |

So roughly **90 % of the band already meets the 0.25 bar** on both heads; the bar is a maximum, and
the top decile carries it over.

**The miss is robust to the choices I was free to make**, which is why it should be read as a real
property and not as a configuration artifact:

| Variation | coloc `maxabs` | exclusion `maxabs` |
|---|---|---|
| shipped config (n = 6000, 200 epochs, per-head grids) | 0.5498 | 0.3783 |
| pooled three-class evaluation grid instead of per-head | 0.4115 | 0.4340 |
| n = 12000, 300 epochs | 0.4767 | 0.1188 |
| n = 20000, 400 epochs | 1.7966 | 0.5774 |
| n = 48000, 200 epochs | 0.509 | 0.159 |
| n = 20000, 400 epochs, width 64 / 16 summaries | 1.3388 | 0.4338 |

The exclusion head — the data-rich one under `P13_F5_SKEW_FREQ` — does clear the bar at n = 12000 and
n = 48000. The coloc head, deliberately starved to 15 % class mass so the correction would be large,
**never** clears it: 0.55 at n = 6000 and 0.51 at n = 48000, an 8× data increase for no improvement.
More data is not the lever.

**This is not a new class of problem for the project.** It is the same phenomenon as named limit #3
— Bayes-factor magnitudes are least determined exactly where the evidence is most confident. Here it
is measured in a setting where the truth is known exactly, which is strictly more information than
the project had before.

**What was NOT done, deliberately.** The toy, the bars, the widths and the epochs were not adjusted
after seeing the miss (threat T-13-25). The failing assertions are left in the committed file, not
softened or commented out. Per D-04 the single pre-declared iteration allowance belongs to the
stratification switch and not to this test, so it is not spent here.

## Recommendation to the phase (not a decision taken here)

Three readings are available and the choice is not this plan's to make:

1. **Accept the shortfall as characterized** and record it as a named limit: the corrected log Bayes
   factors are unbiased in the bulk (mean deviation < 0.07 nat) and lose sub-nat precision in the
   confident tail, where the sign of the verdict is not in doubt anyway.
2. **Amend the F5 statistic** from a maximum to a high quantile (q90 already passes on both heads).
   This would be a *third* amendment in this project's history and the credibility cost is exactly
   why `P13_ITERATION_ALLOWANCE` exists; it should not be done casually and must be declared before
   any further result.
3. **Treat it as a genuine capability limit of BCE heads** and note it in Phase 14, whose abstention
   layer consumes these magnitudes. Phase 14 thresholds sitting in the confident tail would be
   resting on the least-determined part of the estimate.

Option 1 costs nothing and is already fully documented by the committed test's own output.

## Deviations from Plan

### Auto-fixed issues

**1. [Rule 1 - Bug] The plan's symmetry acceptance criterion has a sign error**
- **Found during:** Task 1
- **Issue:** The plan's acceptance criterion asserts
  `isapprox(f(2.0).coloc, -f(-2.0).exclusion; atol = 1e-8)`. That is unsatisfiable for any non-zero
  value: the toy's exact identity is `coloc(Z) == exclusion(−Z)` (reflecting `Z` *swaps the two
  heads*), not a negation. Measured: `f(2.0).coloc = 5.827570`, `f(-2.0).exclusion = 5.827570`,
  `f(2.0).exclusion = -5.610378`. The plan's form would demand `5.8276 == -5.8276`.
- **Fix:** Implemented and asserted the exact form `f(Z).coloc == f(-Z).exclusion` (and its mirror)
  over four `Z` values, with the docstring on `toy_analytic_log_bf` spelling out why the negation
  form is wrong so the error is not reintroduced.
- **Files:** `spike/p13/toy_gaussian.jl`, `spike/test/test_p13_correction.jl`
- **Commits:** 4e7c1fc, a24a913

**2. [Rule 2 - Missing critical] Per-head evaluation grids (`toy_head_grid`)**
- **Found during:** Task 2
- **Issue:** The plan passes ONE `Zgrid` to `toy_compare`. Scoring the coloc head on the pooled
  three-class sample evaluates it deep in the exclusion region, where — under D-11's per-head class
  restriction — it saw *no training data at all* and is extrapolating. That measures a network's
  extrapolation behaviour, not the correctness of a scalar correction.
- **Fix:** Added `toy_head_grid(dataset, positive)` returning a head's own class pair. This is not a
  new convention: the frozen `spike/p13/consts.jl` already commits to it for the reported gate
  ("each head is scored only on its own class pair, per D-11", the sizing rationale for `P13_GATE_M`).
  Both conventions are reported in the table above; **the verdict is unchanged either way.**
- **Files:** `spike/p13/toy_gaussian.jl`

**3. [Rule 2 - Missing critical] `maxabs_demeaned`, so control 2 states the theorem**
- **Found during:** Task 2
- **Issue:** Asserting only that the reshaped head's `maxabs` exceeds the tolerance leaves open the
  reply "then measure a better constant" — which is precisely the claim F2's theorem denies.
- **Fix:** `toy_compare` now also returns the max deviation after removing the best constant in
  hindsight, and control 2 asserts *that* exceeds the tolerance (0.6745). No correction of the D-07
  form could have saved the reshaped head.
- **Files:** `spike/p13/toy_gaussian.jl`, `spike/test/test_p13_correction.jl`

**4. [Rule 2 - Missing critical] Bayes-risk reference (`toy_ideal_logits`) and its testset**
- **Found during:** Task 2
- **Issue:** The plan instructs that a missed bar be recorded as "a phase-level finding about the
  derivation". Without a convergence reference that attribution is unsupported — an under-trained
  toy would produce an identical-looking miss and would say nothing about the derivation.
- **Fix:** Added `toy_ideal_logits` (the analytic Bayes-optimal logits) and the testset
  "a maxabs shortfall is not a training shortfall", which scores them under the same masked loss on
  the same held-out data and asserts the trained net is within 0.01 nat of the Bayes risk. Measured
  excess: **0.001792**. This is what makes the honest finding above attributable.
- **Files:** `spike/p13/toy_gaussian.jl`, `spike/test/test_p13_correction.jl`

**5. [Rule 3 - Blocking] `toy_fit_two_head` takes `seed`, not `rng`**
- **Found during:** Task 1
- **Issue:** The plan sketches `toy_fit_two_head(...; rng)`. A threaded `AbstractRNG` cannot reach
  either stochastic part of a fit: Flux's `glorot_uniform` initializer draws from the **global** RNG
  when the layers are constructed, and NeuralEstimators builds its `DataLoader` with
  `shuffle = true`, which also draws from the global RNG. A threaded `rng` would have been a dead
  parameter and the fit would have been irreproducible.
- **Fix:** The function takes `seed` (defaulting to `P13_FIX_SEED`) and applies layer 2 of
  `P13_GLOBAL_RNG_DISCIPLINE` twice — `Random.seed!` immediately before the layers are built, and
  again inside `train_three_way` immediately before `train`. The determinism testset proves it works
  by deliberately perturbing the global RNG between two fits and asserting they still agree to 1e-6.
  All *data* draws remain threaded through `p13_fix_rng(P13_FIXTURE_COUNTER)` as the plan requires.
- **Files:** `spike/p13/toy_gaussian.jl`

**6. [Rule 1 - Bug] Two of my own first-draft assertion bugs**
- **Found during:** Task 2 (first run)
- **Issue:** (a) `count(===(COLOC), cl)` — Julia does not support partial application of `===`,
  which errored the whole testset. (b) I asserted the reshaped arm's `head_log_odds` equals the
  primary arm's to 1e-9; the two arms shuffle the *same class counts* into the train/validation
  split differently, so the training-subset counts differ slightly (0.8797 vs 0.9105).
- **Fix:** (a) switched to `==(...)`, matching `labels.jl`'s own idiom. (b) the assertion now checks
  both arms against the frozen `P13_F5_SKEW_FREQ` ratio at `atol = 0.1` — which is the claim actually
  being made ("the frequencies did not move, so the failure is attributable to the shape"), and is
  two orders of magnitude tighter than the 0.79-nat deviation the control then exhibits.
- **Files:** `spike/test/test_p13_correction.jl`
- **Commit:** a24a913 (fixed before the file was committed)

### Assumption Drift (advisory)

**A. The pre-registered `maxabs` bar was assumed reachable; it is not, on the starved head.**
- **Planned:** 13-RESEARCH F5 and the plan both present `cor ≥ 0.99` **and** `max|Δ| ≤ 0.25` as
  bars a correct correction clears, with the negative control as the only expected failure.
- **Actual:** the correction is verifiably correct (the offset it removes is removed to 0.07 nat)
  and the max bar is still missed, because a maximum over a band is a confident-tail statistic and
  logit-BCE does not constrain the confident tail.
- **Why it matters:** a reader of the pre-registration would take a `maxabs` miss as evidence
  against the derivation. It is not. The Bayes-risk diagnostic added under deviation 4 is what
  separates the two readings, and without it this SUMMARY would have had to report an ambiguous
  result.

**B. A single evaluation grid was assumed sufficient; D-11 makes it per-head.**
- **Planned:** `toy_compare(net, head_log_odds, Zgrid; ...)`, one grid for both heads.
- **Actual:** under D-11's per-head class restriction the two heads have different support, so a
  shared grid scores at least one of them outside its training distribution.
- **Why it matters:** the fix is a convention the frozen consts already commit to, so nothing was
  invented — but the plan's single-grid signature implied a symmetry the architecture does not have.
  The verdict is identical under both conventions (table above), so nothing rests on it.

## Threat register outcomes

| Threat | Disposition | Outcome |
|---|---|---|
| T-13-22 (wrong sign / denominator / double application) | mitigated | Control 1 asserts the uncorrected offset lands on the head's own log-odds (±0.15) and the exact subtraction identity to 1e-9; a separate testset asserts the denominator is the head's pair and explicitly *not* the one-vs-rest or whole-set form. |
| T-13-23 (vacuous pass) | mitigated | `P13_F5_SKEW_FREQ` asserted lopsided enough that both per-head log-odds exceed the tolerance; control 1 fires on both heads. Correction measured on the training subset only, with the train/full counts asserted to differ. |
| T-13-24 (asserted rather than evidenced stratification) | mitigated | Control 2 fires, and still fires after the best constant is removed (0.6745), which is the theorem rather than an anecdote. |
| T-13-25 (tuning the toy or the bars) | mitigated | Bars asserted as literals against the frozen consts; the two failing assertions are committed as-is; the sensitivity sweep above is reported rather than used to select a passing configuration. |
| T-13-01 (fixture stream reuse) | mitigated | All data draws ride `p13_fix_rng(P13_FIXTURE_COUNTER)`; `Random.seed!` derives from `P13_FIX_SEED` and is called immediately before each build and each `train`. No reported Phase-13 counter is consumed. |
| T-13-SC (package-manager installs) | mitigated | Nothing installed. `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` exits 0. The closed form avoids `QuadGK` and `KernelDensity`, neither of which is in the environment. |

## Verification performed (commands actually run, results actually observed)

| Check | Result |
|---|---|
| `julia --project=spike -e 'include("spike/p13/toy_gaussian.jl"); a = toy_analytic_log_bf(1.5; tau=0.3, s=1.0, sigma=0.5); @assert a.random == 0.0; @assert a.coloc > a.exclusion; println("toy ok")'` | prints `toy ok` ✅ |
| Symmetry sanity (corrected form, see deviation 1) | prints `symmetry ok` ✅ |
| `grep -v '^\s*#' spike/p13/toy_gaussian.jl \| grep -c 'quadgk\|kde('` | `0` ✅ |
| `grep -c 'ThreeWayEvidenceNet\|train_three_way' spike/p13/toy_gaussian.jl` | `6` (≥ 1) ✅ |
| `grep -v '^\s*#' spike/p13/toy_gaussian.jl \| grep -c '0.99\|0.25'` | `0` ✅ |
| `grep -c 'NEGATIVE CONTROL' spike/test/test_p13_correction.jl` | `3` (≥ 2) ✅ |
| `grep -c 'P13_F5_CORR_MIN' spike/test/test_p13_correction.jl` | `4` (≥ 1) ✅ |
| `grep -v '^\s*#' spike/test/test_p13_correction.jl \| grep -c '0\.99'` | `1` (≤ 1) ✅ |
| Both heads' `corr` and `maxabs` printed by the run | ✅ (see table) |
| Runtime under 120 s | **40.8 s** ✅ |
| `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml` | exit 0 ✅ |
| `julia --project=spike spike/test/test_p13_correction.jl` exits 0 with zero failures | ❌ **exits 1: 88 pass, 2 fail** — the two pre-registered `maxabs` bars. This is the honest finding, not an unrun check. |

Not modified, verified by `git status`: `.planning/STATE.md`, `.planning/ROADMAP.md`,
`spike/p13/consts.jl`, `spike/test/runtests.jl`, `src/`, `spike/Project.toml`,
`spike/Manifest.toml`. `spike/test/test_p13_correction.jl` is **not** wired into
`spike/test/runtests.jl` — that is plan 13-08's job, and 13-08 should know that wiring it in as-is
turns the suite red on the two bars above.

## Known Stubs

None. Both files are complete and executable; nothing is placeholdered.

## Threat Flags

None. This plan adds no network endpoint, no auth path, no file access pattern beyond reading
`spike/p13/consts.jl` for its provenance hash (pre-existing behaviour in `net.jl`), and no schema
change at a trust boundary.

## Self-Check: PASSED

- `spike/p13/toy_gaussian.jl` — FOUND (545 lines, ≥ 130 required; contains `toy_analytic_log_bf`)
- `spike/test/test_p13_correction.jl` — FOUND (454 lines; contains `P13_F5_CORR_MIN`)
- commit `4e7c1fc` — FOUND
- commit `a24a913` — FOUND
