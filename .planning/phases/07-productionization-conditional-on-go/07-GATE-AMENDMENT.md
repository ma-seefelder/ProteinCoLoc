# 07 — PRE-REGISTRATION AMENDMENT to the grid-8 ship-gate

**Date:** 2026-07-21
**Status:** **FROZEN SPECIFICATION. Nothing has been run against it.**
**Implements:** `test/gate/gate_consts_8_v2.jl` (new, frozen)
**Supplements — does NOT replace:** `test/gate/gate_consts_8.jl` (and `_4` / `_16`)
**Depends on:** `07-CALIBRATION-FINDINGS.md` (F1–F6), Phase-2 decisions **D-08**, **D-09**

---

## 0. Provenance disclosure — read this first

**This amendment was written AFTER the original grid-8 ship-gate results were seen.** The recorded
verdicts (SBC `ks_pass = false`; BF `corr = 0.9471` and `max|Δ logBF|` both failing) were known to
the author at the time of writing. That is exactly why every change below is argued **only** from
properties of the test design — family-wise error rate, estimator precision, extreme-value
instability, and the validity conditions of Simulation-Based Calibration — and **never** from what
the observed numbers were or from what any net would or would not pass.

The test of legitimacy applied throughout is:

> *Would this change have been made, with this justification, by someone who had seen the gate's
> code but none of its outputs?*

Each item below states its answer to that question explicitly. Where the answer is weaker than
"unambiguously yes", §5 says so.

**The original pre-registration is not erased.** `test/gate/gate_consts_4.jl`,
`gate_consts_8.jl` and `gate_consts_16.jl` remain **byte-unchanged** and their recorded FAIL
verdicts in `artifacts/grid_{4,8,16}/gate_report_{4,8,16}.jld2` remain citable as-is. The amended
gate is a **second, separately named** pre-registration (`gate_consts_8_v2.jl`) scored on a
**fresh, disjoint seed**. Any report that cites the amended result must cite the original FAIL
alongside it. Reporting only the amended run would be exactly the thing this document exists to
prevent.

---

## 1. A1 — the all-8 KS conjunction is uncorrected for multiplicity

### The defect

`sbc_gate` (`test/gate/sbc.jl:328`) computes

```julia
ks_pass = all(x -> x.ks_p > SBC_KS_ALPHA, per)   # per has 8 entries, SBC_KS_ALPHA = 0.05
```

Eight independent-ish hypothesis tests are each run at α = 0.05 and the gate requires **all** to
survive. The probability that a **perfectly calibrated** net fails this conjunction by chance is

```
FWER = 1 − (1 − 0.05)^8 = 1 − 0.95^8 = 1 − 0.66342 = 0.33658 ≈ 33.7 %
```

A correct model fails this gate roughly **one run in three**. The nominal α of the gate as a whole
is not 0.05; it is 0.34.

### Outcome-independent justification

This is arithmetic on the gate's own structure. It requires the number of tests (8, fixed in
`SBC_PARAM_LABELS`) and the per-test α (0.05, frozen in `gate_consts_8.jl:37`) — nothing else. It
is computable from the source file with the artifacts directory deleted. Multiplicity correction
for a conjunction of simultaneous tests is textbook practice that should have been in the original
pre-registration and was not; its absence is a straightforward omission, not a post-hoc discovery.

*Would a blind reviewer have made this change?* **Yes, unambiguously.** An 8-fold uncorrected
conjunction is a recognised design error on sight.

### The new rule

**Holm–Bonferroni step-down, FWER target 0.05, applied separately to the 8 KS p-values and to the
8 χ² p-values.**

Let p₍₁₎ ≤ … ≤ p₍₈₎ be the sorted p-values, m = 8. The Holm adjusted p-values are

```
p̃₍ᵢ₎ = min( 1, max_{j ≤ i} (m − j + 1) · p₍ⱼ₎ )        (running maximum enforces monotonicity)
```

and the arm **passes** iff `all(p̃ > 0.05)` — i.e. Holm rejects no uniformity null at FWER 0.05.
Equivalently in step-down form: reject H₍ᵢ₎ if p₍ᵢ₎ ≤ 0.05 / (8 − i + 1), stopping at the first
non-rejection. The critical values are

| i | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| 0.05/(9−i) | 0.00625 | 0.007143 | 0.008333 | 0.01000 | 0.01250 | 0.016667 | 0.02500 | 0.05000 |

Constants: `SBC_KS_FWER = 0.05`, `SBC_CHI2_FWER = 0.05`, `SBC_N_TESTS = 8`,
`SBC_MULTIPLICITY = :holm`. Implementation: `holm_adjusted` / `holm_pass` in
`gate_consts_8_v2.jl`.

### Why Holm and not something else

| Candidate | Verdict | Reason |
|---|---|---|
| **Holm–Bonferroni** | **CHOSEN** | Controls FWER ≤ 0.05 under **arbitrary dependence** — required, because the 8 columns share the same M simulated datasets and the same posterior draws, so they are dependent by construction and the dependence structure is unknown. Uniformly at least as powerful as plain Bonferroni (step-down dominance: it never rejects less). |
| Plain Bonferroni | rejected | Strictly dominated by Holm at identical FWER. No reason to prefer it. |
| Šidák | rejected | Requires independence of the tests. Not satisfied here (shared draws). |
| Benjamini–Hochberg (FDR) | rejected | Controls the *expected proportion* of false rejections. The gate's claim is a **conjunction** ("all 8 parameters are calibrated"), for which the relevant error concept is family-wise, not false-discovery. FDR would also be *more* permissive here, which is the wrong direction for an amendment written after seeing a failure. |

### Sensitivity of the change

Holm is **more permissive than the original rule** — that is the whole point (the original had a
34% false-failure rate). To bound how much more permissive: under the amended rule the smallest of
eight p-values must fall below **0.00625** before anything is rejected, versus **0.05** before.
This is a ~8× relaxation on the single most extreme column and no relaxation at all on the
conjunction as a *family-wise* statement — the family-wise false-failure rate goes from 33.7% to
≤ 5%, i.e. the amended gate is calibrated to the α it always claimed.

### What is explicitly NOT changed

- `SBC_KS_ALPHA` / `SBC_CHI2_ALPHA` are **retired as gate inputs**, not re-tuned. There is no new,
  looser per-test α; there is a per-*family* α of 0.05, which is the number the original
  pre-registration already claimed to be enforcing.
- `SBC_M = 2000` and `SBC_L = 999` are **unchanged**. Power is not altered in either direction by
  changing the sample size.
- The **ECE conjunction is unchanged**: `ece_pass = all(ece ≤ SBC_ECE_GREEN)` with
  `SBC_ECE_GREEN = 0.05`. ECE is a thresholded point estimate, not a level-α rejection, so no
  family-wise rejection rate accrues and Holm does not apply. `SBC_ECE_YELLOW = 0.10` unchanged.
- The set of 8 tested columns (7 θ + Δρ) is unchanged. No column is dropped.
- The `shrinkage` / `vacuous` diagnostics (F3) remain **reporting-only**, exactly as shipped.

---

## 2. A1b — the χ² bin count (50 bins at M = 2000)

### The defect

`SBC_BINS = 50` is used for **two different things** at once in `sbc_gate`: the χ² rank-bin count
in `sbc_uniformity`, and the reliability-diagram bin count `n_bins` in `sbc_calibration` (which
feeds ECE/MCE). These are unrelated statistics with unrelated binning requirements, and coupling
them means neither can be set on its own merits.

For the **χ²** use specifically, 50 bins at M = 2000 spends **49 degrees of freedom** and yields an
expected count of 40 per bin.

### Outcome-independent justification

Two design arguments, both computable without any result:

1. **Degrees-of-freedom allocation.** The alternatives SBC is designed to detect are
   *low-frequency*: a location shift produces a monotone ramp across the rank histogram; an
   overconfident posterior produces a U; an underconfident one produces an inverted U. None of
   them produces high-frequency bin-to-bin structure. A 50-bin χ² spreads its 49 df across 49
   contrasts, of which only the first two or three carry any calibration meaning; the remaining
   ~46 df are pure noise absorption that **dilutes** power against the alternatives of interest
   while adding power against structure that has no calibration interpretation. Reducing to 20
   bins (19 df) reallocates the same M toward the alternatives SBC exists to detect. This is a
   **power reallocation, not a uniform power reduction**, and its direction is determined entirely
   by which alternatives matter — a question answerable before any data.
2. **Asymptotic-reference accuracy.** The χ² reference distribution is asymptotic in the per-bin
   expected count. At 20 bins the expected count is 100, versus 40 at 50 bins — a materially
   better approximation, and comfortably above every conventional floor.

Divisibility is preserved: `L + 1 = 1000` is divisible by 20 (exactly 50 ranks per bin), so the
rank binning stays exactly even, which the original pre-registration explicitly required
(`gate_consts_8.jl:36`).

*Would a blind reviewer have made this change?* **Probably yes for the decoupling; "defensible
judgement" for 20 vs 50.** 20 is a standard SBC bin count (Talts et al. 2018 use bin counts of
this order); it is not an unusual or self-serving choice. But the amendment does not pretend that
20 is uniquely derivable — see §5.

### The new rule

- `SBC_CHI2_BINS = 20` — the χ² rank-bin count. Used **only** by `sbc_uniformity`.
- `SBC_BINS = 50` — **unchanged**, and now used **only** by `sbc_calibration` (ECE/MCE reliability
  binning). The ECE arm is therefore numerically identical to the original pre-registration.

**Required code change (specified, NOT executed):** `sbc_gate` in `test/gate/sbc.jl` currently
passes the same `bins` to both helpers (`sbc.jl:320-321`). It must be changed to pass
`bins = SBC_CHI2_BINS` to `sbc_uniformity` and `n_bins = SBC_BINS` to `sbc_calibration`. Both
call sites already take the value as a keyword, so this is a call-site change only.

### On "χ² over-power" more honestly

The deeper issue is that a **point null of exact uniformity is known false a priori** for any
finite-capacity normalising flow: at M → ∞ every approximate posterior fails an exact-uniformity
test. Bin count does not fix that; only a practical-equivalence formulation would.

**The amendment deliberately does NOT introduce an equivalence band.** Adding a
"statistically-significant-but-practically-negligible" escape hatch *after* seeing a significant
failure is precisely the move this document must not make, however defensible it would have been
before. `SBC_M = 2000` is therefore held at its original value and the null stays a point null.
The bin-count change is limited to the df-allocation and asymptotic-accuracy arguments above,
which do not depend on the alternative being real or negligible. This limitation is recorded
here so a future pre-registration — written before its own results — can address the point-null
problem properly.

### What is explicitly NOT changed

- `SBC_M`, `SBC_L`, the ECE binning, `SBC_ECE_GREEN`/`SBC_ECE_YELLOW`, the fixture constants
  `SBC_FIX_*`.
- The χ² test itself (`HypothesisTests.ChisqTest` on equal-width rank-bin counts; never
  hand-rolled).

---

## 3. A2 — the BF correlation gate is underpowered by construction

### The defect

`bf_gate` (`test/gate/run_gate.jl:122`) defaults to `n = BF_SWEEP_N = 25` paired draws and then
compares a **point estimate** `cor(amort, base)` to a fixed threshold `BF_CORR_MIN = 0.95`
(`run_gate.jl:146`). After non-finite pairs are dropped the recorded n was 15–20.

At those n, the Fisher-z 95% CI for a correlation near 0.95 is roughly ±0.07–0.10 wide. The CIs of
all three recorded grids **contain the 0.95 threshold**:

| grid | r̂ | Fisher-z 95% CI | contains 0.95? |
|---|---|---|---|
| 4 | 0.941 | [0.854, 0.977] | yes |
| 8 | 0.9471 | [0.845, 0.983] | yes |
| 16 | 0.915 | [0.758, 0.972] | yes |

The gate as designed **cannot discriminate ρ = 0.94 from ρ = 0.95**. Its verdict on that
comparison is a coin flip with the sampling noise as the coin.

A second, independent defect in the same arm: `max_abs_delta = maximum(abs.(amort .- base))` is an
**extreme-value statistic**. Its expectation grows monotonically with n. A fixed tolerance applied
to a maximum therefore silently becomes a *stricter* gate whenever n is raised, for reasons that
have nothing to do with the agreement between the two log-BFs. Any amendment that raises n (as
this one must, to fix the precision defect) would corrupt that component unless it is restated on
an n-stable estimand.

### Outcome-independent justification

The precision defect is derivable from `n` and the Fisher-z variance formula
Var(atanh r̂) ≈ 1/(n − 3) alone. It needs the *sample size* (25, hard-coded in
`gate_consts_8.jl:54`) and the *threshold* (0.95, line 50) — both frozen pre-registration values,
both known before any run. It does not need r̂. The table above uses the observed r̂ only to
*illustrate* a width that is essentially constant over the relevant range of r̂; the finding is
"at n = 15–25 the CI half-width near r = 0.95 is ~0.04–0.10, which exceeds the 0.01 gap the gate
is being asked to resolve" — a statement about the design.

The extreme-value defect is a property of `maximum` as an estimator: E[max] is non-decreasing in
n for any non-degenerate distribution. No result is involved.

*Would a blind reviewer have made these changes?* **Yes.** "You are testing a threshold your
sample size cannot resolve" and "you are gating on a maximum whose distribution depends on n" are
both visible from the source.

### The new rule — sample size

**Precision target (pre-registered):** the two-sided 95% Fisher-z confidence interval for ρ, at
ρ = `BF_CORR_MIN` = 0.95, must have a **half-width ≤ 0.03** on the correlation scale — i.e. the
gate must be able to resolve the threshold to within ±0.03, three times finer than the 0.10-wide
region [0.90, 1.00] the decision lives in, and comfortably finer than the 0.01 discrimination the
original gate attempted at n = 25.

**Derivation** (the binding side is the *lower* limb, because atanh compresses upward):

```
require   r0 − tanh( atanh(r0) − z₀.₉₇₅ / √(n−3) )  ≤  h        with r0 = 0.95, h = 0.03

atanh(0.95) = ½·ln(1.95/0.05) = ½·ln 39      = 1.831781
atanh(0.92) = ½·ln(1.92/0.08) = ½·ln 24      = 1.589027
Δz          = 1.831781 − 1.589027            = 0.242754

√(n−3) ≥ z₀.₉₇₅ / Δz = 1.959964 / 0.242754  = 8.07386
n − 3  ≥ 65.19  ⇒  n ≥ 68.19  ⇒  n_min = 69
```

**Set `BF_GATE_N = 100`** — the derived minimum rounded up to the next round number, which buys
margin without a further doubling of cost. Verification at n = 100:

```
1.959964/√97 = 0.199001
z ∈ [1.831781 ∓ 0.199001] = [1.632780, 2.030782]
r ∈ [tanh 1.632780, tanh 2.030782] = [0.92646, 0.96613]
half-widths: 0.02354 (lower), 0.01613 (upper)   →  both ≤ 0.03  ✓
```

The precision target is stated **at the threshold** ρ = 0.95, which is where discrimination is
required. It is not claimed to hold uniformly: at ρ = 0.90 the lower half-width at n = 100 is
0.045. That is acceptable — a gate needs precision where its decision boundary is, not everywhere.

`BF_GATE_N` is the number of paired draws **attempted**; the recorded `n` after dropping
non-finite pairs must also satisfy `n ≥ BF_GATE_N_MIN = 69` (the derived minimum) or the arm is
recorded `:invalid`, never `:pass`. This closes the failure mode where attrition silently returns
the gate to the underpowered regime.

**Cost note.** n = 25 → 100 quadruples the BF arm's paired draw→simulate→infer work (2 simulations
+ 2 posterior sample calls per draw). Under the F5 image-size mixture (§4) this is ≈ 200 forward
simulations at ≈ 0.083 s each ≈ 17 s of simulation, plus the flow passes. Negligible against the
retraining budget.

### The new rule — decision rule

**CI-based, one-sided, with a three-way verdict:**

```
LB₉₅ = tanh( atanh(r̂) − 1.644854 / √(n−3) )        # one-sided 95% lower confidence bound
UB₉₅ = tanh( atanh(r̂) + 1.644854 / √(n−3) )

corr_verdict = :pass          if LB₉₅ ≥ BF_CORR_MIN            (= 0.95)
             = :fail          if UB₉₅ <  BF_CORR_MIN
             = :inconclusive  otherwise
```

`corr_pass` is true **only** for `:pass`.

**Justification.** The ship claim is "the amortized log-BF agrees with the KDE baseline at
correlation ≥ 0.95". A ship gate should require *evidence for* that claim, not merely *failure to
refute* it; the burden of proof belongs on the thing being shipped. A lower-confidence-bound rule
is the standard formalisation of that burden.

Note the direction: **this makes the BF correlation gate strictly harder than the original.** At
n = 100 the rule requires a point estimate

```
atanh(r̂) ≥ atanh(0.95) + 1.644854/√97 = 1.831781 + 0.167016 = 1.998797
r̂ ≥ tanh(1.998797) = 0.96394
```

versus r̂ ≥ 0.95 before. An amendment written to manufacture a pass would not have done this. The
`:inconclusive` verdict is added because "we cannot tell" is a real and reportable outcome that
the original two-way rule silently converted into a `:fail`; distinguishing it is an honesty gain,
not a leniency gain — `passed` is false in both cases.

**Rejected alternative — keep the point-estimate rule at larger n.** At n = 100 a point rule would
be *more permissive* than the original in expectation (it accepts r̂ = 0.951 as readily as
r̂ = 0.999) while inheriting the original's silence about estimation error. The LB rule was
preferred precisely because it is the conservative option.

### The new rule — `max|Δ logBF|` (kept, restated, not dropped)

The `max|Δ logBF|` component is **retained**. It is restated on an n-stable estimand and the
maximum is retained as a reported quantity:

- **GATING:** `q95_abs_delta = quantile(|Δ logBF|, 0.95) ≤ BF_LOGBF_TOL` with
  **`BF_LOGBF_TOL = 0.5` unchanged** and `BF_LOGBF_Q = 0.95`.
- **REPORTED, NOT GATING:** `max_abs_delta = maximum(|Δ logBF|)`, recorded in the report exactly
  as before, together with `q50` and `q95`, so the tail is fully visible and any claim about it is
  auditable.

**Justification.** The tolerance *value* (0.5 nats — roughly the width of one Jeffreys evidence
category) is **not touched**; only the estimand it applies to changes, from an n-dependent extreme
to an n-stable quantile. Since n rises 25 → 100 in this amendment, leaving the gate on `maximum`
would have tightened it by an amount driven purely by sample size. At n = 100 the 95th percentile
is the 95th-of-100 order statistic — well inside the sample, stable, and still a genuine tail
statistic (it is not a mean, and it does not let a handful of large disagreements hide).

**On the documented clamped-KDE tail artifact:** the gate already scores against the **non-clamped**
baseline `kde_log_bf_unclamped` (`run_gate.jl:136`), which was introduced specifically so
`max|Δ logBF|` would be free of the 1e-8 clamp-tail artifact. The residual artifact is that a KDE
log-density evaluated in the far tail of a finite sample is numerically unstable *regardless* of
clamping — the baseline itself has large variance there, so a single tail draw can dominate a
maximum without indicating any disagreement in the regime that matters. That is an argument about
the **baseline estimator's** tail variance, available from the estimator's definition alone, and
it is the second independent reason the gated statistic is a quantile rather than a maximum. The
maximum is reported so the artifact remains inspectable rather than being defined away.

### What is explicitly NOT changed

- `BF_CORR_MIN = 0.95` — **not re-tuned**. The threshold value is identical.
- `BF_LOGBF_TOL = 0.5` — **not re-tuned**. Identical.
- The baseline: still the **non-clamped** KDE log-BF, still no Turing, still `prior_n = 4000`.
- `BF_SWEEP_LO / BF_SWEEP_HI / BF_SWEEP_N` retained verbatim for the reported Δρ **sweep**
  (a separate reporting product); `BF_GATE_N` is a **new, separate** constant governing the
  paired-draw correlation sample. `BF_SWEEP_N` is no longer overloaded as the gate's n — that
  overloading (`run_gate.jl:122`, `n::Integer = BF_SWEEP_N`) was itself the mechanism by which the
  sample size was never chosen on statistical grounds.

---

## 4. F5 — the gate's simulate distribution must equal the training joint

### The defect

SBC is valid **only** when the ranks are computed under the same joint distribution
p(θ, y) = π(θ)·p(y|θ) the estimator was trained against. The gate's `SBC_IMSIZE = (256,256)`
(`gate_consts_8.jl:41`) is a **scalar** image size, inherited from the *documented compute-budget
deviation* in run 07-05 that constrained training to `imsize_set = ((256,256),)`. It was never a
scientific choice; it was a budget choice that propagated into the pre-registration.

Phase-2 **D-08** anchors the design on the real paper regime (1376×1028). Phase-2 **D-09** states
image dimensions are a configurable parameter to be validated across a realistic range. Training
is moving to a realistic mixture. **If the gate's simulate distribution does not move with it, the
SBC arm is not a calibration proof of the shipped net** — it is a calibration proof of that net on
a distribution it will never see.

### Outcome-independent justification

This one is the strongest of the three: **it is a validity condition of the method, not a
threshold choice.** "SBC ranks are uniform under the training joint" is the theorem (Talts et al.
2018); evaluating ranks under a *different* joint tests a proposition SBC never asserted. The
requirement `gate joint == training joint` is derivable from the definition of SBC with no data
whatsoever, and the decision to move training to a mixture is a **user decision already taken**
on D-09 grounds, independent of this amendment.

*Would a blind reviewer have made this change?* **Yes — a reviewer would have flagged the original
as invalid the moment training moved off 256².** Note the direction: matching the gate to a
harder, more heterogeneous training joint is if anything *harder* to pass than a single fixed
size, and F5 already measured (read-only, on the old net) that the mixture moves rank means
substantially. This amendment knowingly adopts the more demanding regime.

### The new rule

`SBC_IMSIZE` is replaced by a **pre-registered categorical mixture**, drawn per SBC draw:

| size | weight | rationale | measured cost (32 threads, warmed) |
|---|---|---|---|
| 512 × 512 | **0.40** | cheap lower bracket; keeps the mixture affordable and covers small-field acquisitions | 0.0274 s/pair |
| 1024 × 1024 | **0.25** | intermediate bracket, directly below the anchor | 0.0766 s/pair |
| **1376 × 1028** | **0.25** | **the D-08 real-data anchor** — must carry substantial mass or the proof does not speak to the paper regime | 0.0974 s/pair |
| 2048 × 2048 | **0.10** | upper tail; provides coverage above the anchor so the shipped range is bracketed on both sides rather than extrapolated | ≈ 0.29 s/pair (area-extrapolated from the anchor; ×2.965) |

Σ weights = 1.00. **256 × 256 is deliberately excluded** — it is below the realistic acquisition
range and its inclusion in the original was the budget artifact being corrected.

Constants: `SBC_IMSIZE_SET`, `SBC_IMSIZE_WEIGHTS`. The scalar `SBC_IMSIZE` is set to the sentinel
`:mixture` so that any code path still treating it as a scalar tuple **fails loudly** instead of
silently reverting to a single size.

**Weight rationale.** Mass is concentrated on the anchor and its two neighbours (0.50 combined on
1024²/1376×1028) because those are the regimes the manuscript claims cover; 512² carries the
largest single weight because it bounds the mixture's cost and represents a genuinely common
acquisition size; 2048² is capped at 0.10 because it is ~3× the anchor's cost per pair and its
role is bracketing, not representation. Expected cost:

```
E[cost] = 0.40·0.0274 + 0.25·0.0766 + 0.25·0.0974 + 0.10·0.29
        = 0.01096 + 0.01915 + 0.02435 + 0.02900  =  0.08346 s/pair
50 000 pairs ≈ 4173 s ≈ 1.16 h        (vs ≈ 1.35 h for 50k at the anchor alone)
```

The mixture is therefore **cheaper** than an anchor-only pool while covering four sizes — the
weights were chosen to keep the realistic-range requirement affordable, which is a budget argument
about *feasibility*, not about results.

### The binding invariant

**`SBC_IMSIZE_SET` / `SBC_IMSIZE_WEIGHTS` MUST be byte-identical to the `imsize_set` /
`imsize_weights` used to train the net being gated.** The gate must **assert** this against
`ProteinCoLoc.training_imsize_provenance(...)` (persisted per F6) and refuse to run — status
`:provenance_mismatch` — on any mismatch or on `recorded = false`. A net whose training image-size
distribution is `:unknown` (which includes all three currently frozen bundles, per F6) **cannot be
SBC-gated under this amendment at all.** This is the mechanical enforcement of the F5 validity
condition; without it the condition is a comment.

### Required code changes (specified, NOT executed)

The gate currently threads `imsize` as a **scalar** to `sim.simulate_pair(rng, θ; imsize = imsize)`
in four places. Required changes:

1. **`test/gate/harness.jl:126` `draw_simulate_infer`** — replace the scalar `imsize` keyword with
   a per-call draw
   `isz = ProteinCoLoc.sample_imsize(rng; imsize_set = SBC_IMSIZE_SET, imsize_weights = SBC_IMSIZE_WEIGHTS)`
   taken from the **gate rng** (`prod_rng`), so the size draw is part of the reproducible stream.
   Return `isz` in the NamedTuple so it is recordable.
2. **`test/gate/harness.jl:148` `draw_simulate_infer_paired`** — draw the size **once** and use it
   for **both** members of the pair. A real sample/control pair comes from one acquisition
   configuration; drawing two independent sizes would inject a confound the forward model does not
   contain.
3. **`test/gate/sbc.jl` `sbc_ranks_and_spread`** — thread the mixture rather than a scalar; record
   the realised size per draw alongside the rank table so the realised mixture is auditable
   post-hoc (the F6 lesson: record what was actually drawn, never reconstruct it).
4. **`test/gate/run_gate.jl` `bf_gate` / `ood_gate` / `run_gate`** — same substitution; the
   `imsize = SBC_IMSIZE` defaults become `imsize_set = SBC_IMSIZE_SET,
   imsize_weights = SBC_IMSIZE_WEIGHTS`. Add the provenance assertion of the previous subsection to
   `run_gate` before any arm executes.
5. **New report fields:** `imsize_set`, `imsize_weights`, `realised_imsize_counts`.

`ProteinCoLoc.sample_imsize` (`src/amortized/datagen.jl:108`) already implements exactly the
required inverse-CDF categorical draw and is reused verbatim — the gate must not hand-roll a
second sampler, or the gate and training joints could drift apart in implementation even when the
constants agree.

### What is explicitly NOT changed

- The forward model, the prior π(θ), the summary statistic, `patch_summary`, the 8×8 grid, the
  `encode_d01` / `standardize_summary` chain.
- `SBC_M`, `SBC_L`, the OOD constants (`OOD_ID_QUANTILE`, `OOD_KS_EPS`, `OOD_AUC_MIN`,
  `OOD_GRID_LEVELS`, `OOD_PP_REPS`) — all carried over verbatim. The OOD arm PASSED and nothing
  about it is being amended; it is re-run only because the image-size mixture changes its input
  distribution.

---

## 5. Why this is not p-hacking — and the residual risk, stated honestly

### The argument

1. **Every change is justified from a property of the test, not a property of a result.** A1 from
   the family-wise error rate of an 8-fold conjunction (arithmetic on the source). A1b from
   degrees-of-freedom allocation and χ² asymptotics. A2 from the Fisher-z variance of a
   correlation at n = 25 and from E[max] growing with n. F5 from the validity condition of SBC and
   from decisions D-08/D-09 that predate every gate run. Delete the `artifacts/` directory and
   every argument above still stands verbatim.
2. **Two of the three changes make the gate HARDER.** The BF correlation rule moves from
   "r̂ ≥ 0.95" to "one-sided 95% LB ≥ 0.95", which at n = 100 demands r̂ ≥ 0.9639. The F5 mixture
   replaces a single small image size with a heterogeneous realistic mixture that F5 already
   measured (read-only) to be substantially more demanding for this net. A gate amended to
   manufacture a pass does not do either of these.
3. **No threshold value was re-tuned.** `BF_CORR_MIN` stays 0.95. `BF_LOGBF_TOL` stays 0.5.
   `SBC_ECE_GREEN` stays 0.05. `SBC_M` stays 2000, `SBC_L` stays 999. The family-wise α stays
   0.05 — it is now actually 0.05 rather than nominally 0.05 and really 0.34.
4. **No component was dropped.** `max|Δ logBF|` is retained (as a reported quantity with a
   quantile-based gating twin). The χ² arm is retained. All 8 columns are retained. The ECE
   conjunction is retained unchanged.
5. **The original record is intact and must be co-cited.** §0.
6. **The amended gate runs ONCE, on a fresh seed.** §6.

### The residual risk — not argued away

- **The amendment is post-hoc. That is unfixable.** No amount of outcome-independent argument
  makes an amendment written after seeing a failure equivalent to one written before. The correct
  epistemic weight on the amended result is **lower** than on a genuinely prior registration, and
  it should be reported that way.
- **A1 is the change with real leniency in it.** It reduces the family-wise false-failure rate from
  33.7% to ≤ 5%, which necessarily means some outcomes that would have failed the old rule pass
  the new one. The justification (the old rate was wrong) is sound, but the *direction* of the
  benefit is not neutral and the author knew which direction he needed.
- **"20 bins" is a defensible choice, not a derived one.** The decoupling of χ² binning from ECE
  binning is forced; the specific value 20 is judgement, informed by convention. A different
  reviewer might have chosen 10 or 25. Nothing in the argument uniquely selects 20.
- **The mixture weights are judgement.** 0.40/0.25/0.25/0.10 is reasoned (anchor mass, bracketing,
  cost) but not uniquely determined. A different defensible weighting exists. What is *not*
  judgement is that the gate mixture must equal the training mixture — that part is forced.
- **This document does not fix the point-null problem** (§2). The SBC arm still tests a null known
  to be false asymptotically. Fixing it after seeing a failure would have been indefensible;
  leaving it means the SBC arm remains harsher than it ideally should be, in a way that is *not*
  in the author's favour.
- **The bounded-θ retrain is a confound.** The net that will face this gate is not the net that
  failed it: the F2 remedy (`BoundedThetaTransform`) has been implemented since. The amended run
  therefore changes **two** things at once — the gate rules and the model. It cannot attribute a
  future pass to either alone. This is accepted deliberately (re-running the old net under the new
  rules would burn the fresh seed on a net that is not the ship candidate), but it must be stated
  in any report of the outcome.

---

## 6. Required confirmation protocol (binding)

1. **Retrain first.** Grid 8 is retrained with the bounded-θ parameterisation on the F5 image-size
   mixture (`SBC_IMSIZE_SET` / `SBC_IMSIZE_WEIGHTS`), with `imsize_set` / `imsize_weights` /
   `imsize_source` **persisted** into the artifact `meta` (F6 mitigation).
2. **Freeze before running.** `test/gate/gate_consts_8_v2.jl` and this document are committed
   **before** the amended gate is invoked. The commit hash of that freeze is cited in the result.
3. **ONE run.** The amended gate is run **exactly once** on the **fresh** `PROD_SEED_V2[8]`,
   unit-asserted disjoint from every existing `PROD_SEED[G]` (G ∈ 4, 8, 16, 32), from
   `VAL_MASTER_SEED = 0x5BC0FFEE`, from `NPE_MASTER_SEED = 0xC0FFEE`, from
   `DEFAULT_MASTER_SEED = 0x1`, and from the dev seeds `0xDE7C0DE` / `0xDE7C0DE2` used by the F2
   remedy diagnostic.
4. **No iteration.** If the amended gate fails, the failure is recorded and reported. It is **not**
   followed by a second amendment, a second seed, or a re-tune. A further amendment to this gate
   would exhaust whatever credibility the pre-registration mechanism still carries and the
   Go/No-Go memo must say so.
5. **Co-citation.** Any report of the amended verdict cites the original grid-8 FAIL, this
   amendment, and the confound in §5 (rules and model both changed).
6. **Diagnostic separation.** All development-time evaluation of the retrained net continues on
   the dev seeds. `PROD_SEED_V2[8]` is touched exactly once, by the confirmation run.

---

## 7. Scope of this amendment

- **Amends:** the grid-8 ship-gate only, via the new frozen `test/gate/gate_consts_8_v2.jl`.
- **Does NOT amend:** `gate_consts_4.jl`, `gate_consts_8.jl`, `gate_consts_16.jl`,
  `gate_consts_template.jl` — all byte-unchanged. Grids 4 and 16 are **not** re-gated by this
  document; if they are ever re-gated the same amendment must be instantiated for them
  explicitly, with their own fresh seeds.
- **Does NOT touch:** `spike/` (byte-untouched), `src/` model code, any recorded artifact.
- **Nothing was run against these rules.** No gate arm was invoked, no net was trained, no net was
  evaluated against any threshold in this document. The numbers in §1–§4 are all derived from
  closed-form design formulas (FWER arithmetic, Fisher-z variance, atanh/tanh, the measured
  datagen cost table) and none of them is an outcome.
