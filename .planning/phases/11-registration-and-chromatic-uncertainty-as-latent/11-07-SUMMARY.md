---
phase: 11-registration-and-chromatic-uncertainty-as-latent
plan: 07
subsystem: research-net
tags: [datagen, training, npe, lambda-conditioning, tripwire, blocker, d-01, d-02, d-03, d-09]

# Dependency graph
requires:
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 06
    provides: "DECISION proceed — the D-06 abort criterion did not fire, authorising the ~1.5 h datagen-plus-training spend"
  - phase: 11-registration-and-chromatic-uncertainty-as-latent
    plan: 04
    provides: "build_p11_estimator / fit_p11_theta_transform / encode_lambda / augment_input, and the SC1g tripwire spike/test/test_lambda_ablation.jl"
provides:
  - "spike/data/cache/p11/c12f8ab2…/ — the 50 000-pair lambda-hierarchical training pool on the F5 mixture (gitignored bulk cache, resumable)"
  - "spike/npe/p11_research_npe.jld2 — the trained research net (d_in = 129, D = 8) with realized training_imsize_provenance and the four declared deviations"
  - "MEASURED: the lambda conditioning input is ALIVE — the shift marginals track lambda (SD ratios 4.84–9.05 against an ideal prior ratio of 12.0)"
  - "MEASURED: the rho_true posterior width does NOT track lambda (SD ratios 1.20 / 1.06 / 0.98)"
  - "BLOCKER: the pre-registered SC1g tripwire FAILED at P11_LAMBDA_ABLATION_FACTOR = 2.502; the plan's prescribed (a)/(b)/(c) diagnosis came back GREEN on all three legs"
affects: [11-08-ladder, 11-09-breakdown, 11-10-real-image, 11-11-report]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "A behavioural tripwire whose statistic is the SAME quantity the phase's headline criterion measures cannot separate 'the mechanism is dead' from 'the effect is null' — the separation needs a statistic on a channel whose response is analytically known"
    - "A stall-signature diagnostic keyed on the early-stopping epoch fires on benign overfitting, so its sanctioned remedy must be gated on the mechanism actually having failed"

key-files:
  created:
    - spike/npe/p11_research_npe.jld2
  modified:
    - spike/test/test_lambda_ablation.jl

key-decisions:
  - "The tripwire failure is NOT reported as an implementation bug: all three prescribed diagnostic legs came back green and a fourth, decisive check proves the 129th input row is read and used"
  - "The R11 stall remedy (reduce LAMBDA_MAX) was NOT executed — its premise is falsified by the diagnosis, LAMBDA_MAX is Tier-1 pre-registration that may never be edited, and it would have cost another ~72 min of compute"
  - "P11_LAMBDA_ABLATION_FACTOR was NOT relaxed, and no threshold, seed or pre-registered constant was moved"

requirements-completed: []

# Metrics
duration: ~80min
completed: 2026-07-27
---

# Phase 11 Plan 07: Lambda-Hierarchical Pool and Research-Net Training Summary

**The 50 000-pair lambda-hierarchical pool generated in 56.37 min against a pre-registered
150-min ceiling and the research net trained CPU-only in 16.08 min at `d_in = 129, D = 8`, but the
pre-registered SC1g lambda-ablation tripwire FAILED (Δρ posterior-SD ratios 1.109 / 0.987 / 0.936
against a required 2.502) — and the plan's prescribed (a)/(b)/(c) diagnosis came back GREEN on all
three legs, while a fourth check shows the 129th input row IS read and used (shift marginals widen
by 4.84–9.05× against an ideal prior ratio of 12.0). The conditioning is alive; it is the ρ_true
width that does not track λ. This plan therefore ends BLOCKED, not complete.**

## Status: **BLOCKED at Task 3's gate**

Tasks 1 and 2 were completed by a prior executor and are re-verified green here. Task 3's three
steps all executed and produced real measurements, but its `<verify>` command
(`julia --project=spike spike/test/test_lambda_ablation.jl`) exits **1**. The plan's own
instruction on that branch is to diagnose in a prescribed order and record a BLOCKER, which is what
this summary does. **No ladder was run and no threshold was relaxed.**

## Performance — measured, not estimated

| Stage | Measured | Pre-registered ceiling | Verdict |
|---|---|---|---|
| Datagen (pool generation) | **56.37 min** | `P11_DATAGEN_WALLCLOCK_CEILING_MIN = 150` | **well under** — no blocker |
| Datagen (process wall, incl. startup/compile) | 56.78 min (12:33:20Z → 13:30:07Z) | — | — |
| Training | **16.08 min** (process wall 16.5 min, 13:33:20Z → 13:49:50Z) | none pre-registered | — |
| Tripwire | 3.0 s of assertions after ~40 s of setup | — | FAIL (exit 1) |
| **Total compute** | **~72.5 min** | ~1.5 h budgeted (§D14) | on budget |

**The wall-clock ceiling was NOT exceeded and the image-size arm was NOT downgraded.** Datagen ran
at 32 threads, 5 shards × 10 000, serial-over-shards with parallel fill:

| Shard | Samples | Shard min | Cumulative min |
|---|---|---|---|
| 1/5 | 10 000 | 11.21 | 11.21 |
| 2/5 | 10 000 | 10.90 | 22.11 |
| 3/5 | 10 000 | 11.39 | 33.51 |
| 4/5 | 10 000 | 11.54 | 45.04 |
| 5/5 | 10 000 | 11.31 | 56.35 |

Measured rate 0.0676 s/pair against the frozen cost model's 0.08346 s/pair — i.e. the run came in
**19 % cheaper** than the model that set the ceiling. The resume-by-skip path was never exercised
(`5 generated, 0 resumed`); the run was not interrupted.

## Realized image-size distribution (F5 / T-11-28)

Persisted as `training_imsize_provenance` with **realized counts**, not intended weights:

| Image size | Realized count | Share | Intended weight |
|---|---|---|---|
| (512, 512) | 20 069 | 0.4014 | 0.40 |
| (1024, 1024) | 12 604 | 0.2521 | 0.25 |
| **(1376, 1028)** | **12 289** | 0.2458 | 0.25 |
| (2048, 2048) | 5 038 | 0.1008 | 0.10 |
| total | 50 000 | 1.0000 | — |

The **1376×1028 real-data anchor is present** at 12 289 samples, so the D-17 real-image check
retains its basis for transfer. The persisted provenance shows a **mixture, not a single size** —
no 256²-only fallback was used anywhere.

## Training — validation-risk trajectory and early stopping

- **Initial validation risk:** 15.645162
- **Best validation risk:** **6.0490127**, at **epoch 22**
- **Epochs run:** **63** of a 300 budget
- **Early stopping fired at epoch 63** ("validation loss has not improved in 40 epochs",
  `stopping_epochs = 40` → best at 22, patience exhausted at 62/63)
- Recipe unchanged and CPU-only: `epochs=300 batchsize=128 lr=2.5e-4 wd=1e-4 patience=40`,
  `use_gpu = false` on every NeuralEstimators call, `[ Info: Running on CPUDevice`
- Split: 40 000 train / 10 000 val (K = 5, fold 1 held out)
- Epoch time ~13.3–15.3 s after a 56.3 s first epoch (compile)

Trajectory, abridged (train, validation):

| Epoch | 0 | 1 | 5 | 10 | 15 | **22** | 30 | 40 | 50 | 63 |
|---|---|---|---|---|---|---|---|---|---|---|
| train | 15.645 | 8.580 | 6.615 | 6.121 | 5.735 | **5.264** | 4.949 | 4.512 | 4.053 | 3.465 |
| val | 15.645 | 7.660 | 6.914 | 6.492 | 6.325 | **6.049** | 6.394 | 6.821 | 8.233 | 10.892 |

The shape is **benign overfitting after epoch 22**, not a failure to learn: validation risk fell
by 61 % from its initial value, then diverged while training risk kept falling.

### The R11 stall signature fired, and was deliberately NOT acted on

`_stall_diagnostic` printed the R11 banner because early stopping fired at epoch 63 < 80. Its
first leg did **not** fire (`validation risk improved: yes by epoch 50`). The sanctioned response
in the plan is "reduce `LAMBDA_MAX` and report the reduced ladder". **It was not executed**, for
three reasons recorded here as a declared deviation (see Deviations, item 2):

1. **`LAMBDA_MAX` is Tier-1 pre-registration.** `spike/validation/p11_consts.jl:36-37` states that
   no constant in either guard block is ever *edited* and that "a diff that MODIFIES a line below
   is, by construction, a pre-registration breach". Reducing it in place is exactly that.
2. **The remedy's premise is falsified by the diagnosis below.** R11 exists for a net that cannot
   learn the conditioning. This net learned it: the shift marginals recover 40–75 % of the ideal
   λ-response ratio. There is no learning failure for a narrower ladder to remedy.
3. **It would move the measurement in the wrong direction.** The tripwire's bar is a ratio measured
   *across* the λ span. Shrinking `LAMBDA_MAX` shrinks that span, making a 2.502× ratio strictly
   harder to clear, and would cost a fresh 56-min pool (the content hash keys on `lambda_max`) plus
   a fresh training run.

## The SC1g tripwire — FAILED

`julia --project=spike spike/test/test_lambda_ablation.jl` → **exit 1**, 9 passed / 8 failed.

Required factor **2.501971634054976** (Tier-2, probe-derived, pre-registered). 2000 draws/arm,
3 datasets, imsize (256, 256), global seed `0x71040a563a8dabea`, arms differing **only** in the
129th row (identical summaries, re-pinned global RNG).

| dataset | sd(Δρ\|λ=0.25) | sd(Δρ\|λ=3.0) | **sd ratio** | hdi(0.25) | hdi(3.0) | **hdi ratio** |
|---|---|---|---|---|---|---|
| 1 | 0.08601 | 0.09536 | **1.10872** | 0.27969 | 0.31385 | **1.12210** |
| 2 | 0.08958 | 0.08840 | **0.98674** | 0.29702 | 0.29228 | **0.98405** |
| 3 | 0.10430 | 0.09762 | **0.93589** | 0.34282 | 0.32213 | **0.93964** |

All three datasets fail on both spread statistics; two of the three actually **narrow** slightly.
The three testsets that passed are the ones that matter for trusting the number: the
both-directions control ("the statistic can FAIL") passed, including the byte-identical-draws
reproducibility proof for the global-RNG pin; the factor-provenance testset confirmed the Tier-2
value was in force, not the 1.15 fallback; and the CPU-only clause passed.

## The prescribed diagnosis — all three legs GREEN

The plan directs a specific diagnostic order on failure. Every check below was run and its real
output observed.

### (a) Is λ drawn FIRST, on the **realized pool** (not just on fresh draws)?  **GREEN**

Loaded all 50 000 generated columns from the content-hash pool directory:

| Check | Observed | Bar |
|---|---|---|
| `abs(shift_dx) <= lambda`, all samples | **true**, 0 violations of 50 000 | all |
| `abs(shift_dy) <= lambda`, all samples | **true**, 0 violations of 50 000 | all |
| `cor(lambda, abs(shift_dx))` | **0.6042** | > 0.3 |
| `cor(lambda, abs(shift_dy))` | **0.6114** | > 0.3 |
| `cor(lambda, abs(chromatic_eps))` | **−0.00655** | ≈ 0 (D-09) |
| λ range | (0.2500, 3.0000) | `[LAMBDA_MIN, LAMBDA_MAX]` |
| theta rows / summary rows | 8 / 128 | 8 / 128 |

The hierarchical joint is intact in the data that was actually trained on, and `chromatic_eps` is
confirmed **not** λ-scaled.

### (b) Were the training inputs 129 rows, with a **varying** 129th row?  **GREEN**

| Check | Observed |
|---|---|
| training input rows | **129** |
| distinct values in row 129 | **39 968** of 40 000 columns |
| row-129 range | (2.1663e-5, 0.99996) |
| row-129 std | 0.28868 — i.e. exactly `1/sqrt(12)`, the SD of `Uniform(0,1)` |

Independently corroborated by the training log's own runtime assertion:
`[3/5] 129-row inputs assembled; lambda row varies over 39968 distinct values`.

### (c) Was `augment_input` applied **AFTER** standardization?  **GREEN**

| Check | Observed |
|---|---|
| row 129 `==` `encode_lambda(lambda)` **exactly** | **true** (un-z-scored, data-independent affine map) |
| rows 1:128 `==` the standardized summary **exactly** | **true** |
| mean / std of row 129 | 0.49897 / 0.28868 — **not** ≈ 0 / ≈ 1, so it was never z-scored |
| `encode_lambda(LAMBDA_MIN)`, `encode_lambda(LAMBDA_MAX)` | 0.0, 1.0 |

Had the append happened before standardization, row 129 would carry a ≈0-mean/≈1-SD signature and
would not reproduce `encode_lambda` exactly. It does, to the bit.

### (d) The added, decisive check: **does the net read row 129 at all?**  **IT DOES**

The three prescribed legs cannot distinguish "the input is dead" from "the input is alive but ρ
width does not respond", because both leave (a), (b) and (c) green. One channel settles it. The
**shift** marginals have an analytically known λ response: their conditional prior *is*
`Uniform(-λ, +λ)` and the fixed 8×8 correlation summary is vacuous about them (D-05), so a net
that reads row 129 **must** widen them by roughly `LAMBDA_MAX/LAMBDA_MIN = 12.0`, while a net that
ignores row 129 must leave them flat.

Posterior SD per marginal, λ = 0.25 vs λ = 3.0, same pinned RNG, 2000 draws:

| marginal | ds1 ratio | ds2 ratio | ds3 ratio | reading |
|---|---|---|---|---|
| **shift_dx** | **4.841** | **6.988** | **8.045** | **steep — row 129 is read** |
| **shift_dy** | **3.984** | **6.525** | **9.047** | **steep — row 129 is read** |
| rho_true | 1.203 | 1.062 | 0.983 | flat |
| label_eff | 1.424 | 0.943 | 1.143 | flat |
| spillover | 1.109 | 0.883 | 1.005 | flat |
| autofl | 1.051 | 0.925 | 0.925 | flat |
| noise | 1.044 | 1.037 | 1.046 | flat |
| chrom_eps | 1.065 | 0.997 | 0.936 | flat |
| *(ideal, from the conditional prior)* | *12.0* | *12.0* | *12.0* | — |

Raw shift SDs move from 0.171–0.463 at λ = 0.25 to 1.212–1.933 at λ = 3.0, against conditional
prior SDs of 0.1443 and 1.7321 — i.e. the net's shift posterior tracks the declared uncertainty
into the right absolute neighbourhood, not merely in the right direction.

**Conclusion of the diagnosis.** The λ conditioning is **alive and demonstrably used**. There is no
implementation bug in the hierarchical sampler, the 129-row assembly, or the standardize-then-append
ordering. The tripwire fails because it measures **Δρ** width, and ρ_true width does not track λ on
this net.

## Why this is a BLOCKER and not a finding this executor may resolve

The plan's Task 3 pre-commits to a reading: "a failure here is an implementation bug, not a
scientific result, and must not be reported as one." The diagnosis contradicts that reading on its
own prescribed evidence. Resolving the contradiction is a scientific and pre-registration decision,
not execution mechanics, so it is escalated rather than decided here. Specifically, **all** of the
following are outside this executor's authority and were **not** done:

- `P11_LAMBDA_ABLATION_FACTOR` (2.502) was **not** relaxed, and the 11-06 verdict's condition that
  "a failure there in plan 11-07 is a **result**, not licence to relax the constant" was honoured.
- `LAMBDA_MAX`, `LAMBDA_MIN`, `SC2_SPEARMAN_FLOOR` and every other Tier-1/Tier-2 constant are
  **byte-unchanged**; `spike/validation/p11_consts.jl` is not in this plan's diff.
- No re-seed, no re-run of the tripwire at a different seed, no second training run.
- **No ladder was run.** Wave 7 remains gated, as designed.
- The failure is **not** recorded as "SC2 is below the fixed summary's resolution". That branch has
  a pre-registered Tier-1 owner (the D-06 abort criterion, which did not fire) and it is not this
  plan's to declare.

The substantive open question for the orchestrator/user is that **the tripwire's statistic is the
same quantity SC2 itself measures**. As authored, SC1g can only pass if SC2's headline effect
already holds, so it cannot perform the separation its own header claims ("the ONLY thing in the
phase that tells them apart"). The 11-PATTERNS.md §4 note anticipated the risk of picking the
"right spirit, wrong statistic" and this is that case. Diagnostic (d) is a candidate statistic that
*does* separate the two, on a channel with an analytically known answer.

## Task Commits

| Task | Name | Commit | Files |
|---|---|---|---|
| 1 | λ-hierarchical datagen module + fixture test | `efb30f3` (prior executor) | `spike/data/p11_generate.jl`, `spike/test/test_p11_generate.jl` |
| 2 | spike-local research trainer at d_in = 129, D = 8 | `8fcc9b8` (prior executor) | `spike/npe/train_p11_npe.jl` |
| — | Rule 1 fix: tripwire column formatter | `59d3ed8` | `spike/test/test_lambda_ablation.jl` |
| 3 | datagen + training run (gate FAILED) | `fdf32ed` | `spike/npe/p11_research_npe.jld2` |

## The persisted artifact

`spike/npe/p11_research_npe.jld2` — **committed** (4 933 362 bytes). Repository policy tracks spike
net artifacts under `spike/npe/`: `git check-ignore -v` returns exit 1 for both this file and its
sibling `spike/npe/trained_npe.jld2`, and `git ls-files` confirms `trained_npe.jld2` is tracked.
The plan's parenthetical ("as it does for `spike/npe/trained_npe.jld2`") assumed an ignore rule that
does not exist; the actual policy was checked rather than assumed. **No ignore rule was broadened
or narrowed.**

- `generated` = `2026-07-27T13:49:47.473Z`
- `d_in` = 129, `D` = 8, `length(P11_RESEARCH_NET_DEVIATIONS)` = 4, `elapsed_min` = 16.0800500
- `training_imsize_provenance` = `Dict((2048,2048)=>5038, (1024,1024)=>12604, (512,512)=>20069, (1376,1028)=>12289)`
- The 54 MB pool cache under `spike/data/cache/p11/` stays **untracked** — `.gitignore:spike/data/cache/*`
  already excludes it (a pre-existing rule, untouched).

## Verification — observed results

Every command below was run and its real output observed. Nothing is inferred from code reading.

| # | Check | Command | Observed |
|---|---|---|---|
| 1 | Task 1 re-verify | `julia --project=spike spike/test/test_p11_generate.jl` | **47 pass / 0 fail**, `EXIT=0`, 14.5 s |
| 2 | Task 2 re-verify | `julia --project=spike -e 'include(".../train_p11_npe.jl"); @assert isdefined(…, :load_p11_npe); …'` | `trainer loads` |
| 3 | Datagen within ceiling | the run's own per-shard guard | 56.37 min vs 150; guard never threw |
| 4 | Net payload | the plan's `d_in`/`D`/provenance/4-deviations assertion | `net ok, elapsed_min=16.080050003528594` |
| 5 | Provenance non-empty + anchor | the plan's `sum(values(p)) > 0` assertion | prints the 4-size mixture incl. (1376,1028) |
| 6 | **Tripwire (the gate)** | `julia --project=spike spike/test/test_lambda_ablation.jl` | **exit 1 — 9 pass / 8 FAIL** |
| 7 | Diagnosis (a)/(b)/(c)/(d) | ad-hoc read-only script on the realized pool + net | all three legs green; (d) shows the input alive |
| 8 | Decoupling | `git diff --stat -- src/ spike/Project.toml spike/Manifest.toml` | **empty** |
| 9 | No deletions | `git diff --diff-filter=D --name-only HEAD~1 HEAD` per commit | empty each time |
| 10 | Clean tree | `git status --short` | empty after both commits |

**NOT verified, and deliberately not claimed:**

- **No SC2 claim of any kind.** The ladder was not run. The flat ρ_true response above is a
  3-dataset, single-imsize, fixture-counter diagnostic at (256, 256) — it is **not** the SC2
  measurement, which is 500 paired draws per rung across 7 rungs on the F5 mixture.
- No calibration, coverage, SBC, monotonicity, TOST or Bayes-factor claim.
- The shift and `chromatic_eps` θ columns remain **vacuous by design** (D-05); diagnostic (d) uses
  them precisely *because* the summary is uninformative about them, and their λ response is
  evidence about the **network's plumbing**, not about identifiability.
- `spike/test/runtests.jl` was **not** run to green — the pre-existing `SPEEDUP_GATE` blocker
  stands and was not lowered.
- The pool's byte-identity-under-thread-count property is asserted by `test_p11_generate.jl` at
  fixture scale; it was **not** re-verified at 50 000 scale (that would cost a second 56-min run).

## Decisions Made

1. **Report the tripwire failure as a BLOCKER with a diagnosis that contradicts the plan's
   pre-committed reading, rather than forcing the evidence to fit it.** The plan says a failure is
   an implementation bug. Its own three diagnostics say otherwise. Writing "implementation bug"
   anyway would have been the transcription error the phase's whole two-tier pre-registration
   discipline exists to prevent.
2. **Add a fourth diagnostic rather than stopping at three.** (a)+(b)+(c) green is consistent with
   *both* readings, so the prescribed set is not decisive on its own. The shift channel has an
   analytically known λ response, which makes it a genuine discriminator.
3. **Do not execute the R11 `LAMBDA_MAX` reduction.** Premise falsified, pre-registration would be
   breached, and the change moves the measured ratio the wrong way (three independent reasons, any
   one sufficient).
4. **Do not relax, re-seed or re-run.** One run, no retry — the protocol the phase inherited from
   the 07 gate amendment.
5. **Commit the trained net despite the gate failing.** 72.5 min of authorised compute produced a
   real, atomically written, provenance-complete artifact. Discarding it would force a re-spend to
   ask any follow-up question, and the artifact is the evidence the diagnosis rests on.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 — Bug] The tripwire's column formatter crashed on its own header labels**
- **Found during:** Task 3, Step 3 (first-ever execution of the tripwire with a net present)
- **Issue:** `_f(x) = rpad(string(round(x; digits = 5)), 12)` at
  `spike/test/test_lambda_ablation.jl:179` was applied to the `String` header labels at line 191,
  throwing `MethodError: no method matching round(::String, ::RoundingMode{:Nearest}; digits)`.
  The throw landed **after** `ABLATION_ROWS` had already been computed, so a completed ablation
  measurement was discarded at the print step — the file was unrunnable end-to-end and could never
  have reported a verdict.
- **Fix:** added an `_f(x::AbstractString)` method that pads a label verbatim; numeric cells
  unchanged. No assertion, threshold or statistic was touched.
- **Files modified:** `spike/test/test_lambda_ablation.jl`
- **Commit:** `59d3ed8`

### Declared deviations from the plan's instructions

**2. The R11 stall remedy was NOT executed**
- **Found during:** Task 3, Step 2
- **Plan instruction:** on early stopping before epoch 80, "reduce `LAMBDA_MAX` and report the
  reduced ladder".
- **Actual:** early stopping fired at epoch 63, so the signature did fire; `LAMBDA_MAX` was left at
  3.0 and no reduced ladder was produced.
- **Why:** (i) `LAMBDA_MAX` is Tier-1 pre-registration whose own header forbids in-place edits;
  (ii) the remedy targets a learning failure that diagnostic (d) shows did not occur; (iii) a
  narrower span makes the 2.502× ratio harder, not easier, and costs a fresh 56-min pool plus a
  fresh training run. Capacity was **not** raised either (spike 012, R11).
- **Files modified:** none

**3. `.planning/STATE.md` not written (orchestrator-owned)**
- **Found during:** Task 3, Step 3
- **Issue:** Task 3 directs a BLOCKER entry in `.planning/STATE.md`; the parallel-execution
  contract reserves that file to the orchestrator.
- **Resolution:** left untouched; the exact proposed BLOCKER text is reproduced below for the
  orchestrator. `.planning/ROADMAP.md` and `.planning/REQUIREMENTS.md` likewise untouched — and
  **no requirement is marked complete**, because the plan's gate did not pass.
- **Files modified:** none

**4. This plan was RESUMED from Task 3**
- Tasks 1 and 2 were executed and committed by a prior executor (`efb30f3`, `8fcc9b8`) whose
  worktree was merged. They were **not** re-authored; both were re-verified green (Verification
  rows 1–2) before any compute was spent.
- The worktree spawned at `0d9b87c`, an **ancestor** of the declared base. The startup base
  assertion caught it and `git reset --hard 8fcc9b8` corrected it on a clean tree before any work.
  This mattered: executing from the stale base would have found neither `p11_generate.jl` nor
  `train_p11_npe.jl`.

## Assumption Drift (advisory)

**1. "A failing λ-ablation tripwire means the λ conditioning is dead."**
- **Found during:** Task 3, Step 3
- **Planned:** the plan `<action>`, the tripwire's own header and the datagen module's header all
  assume the tripwire is a *liveness* check on the 129th input row — "the ONLY thing in the phase
  that tells [a bug and a null] apart".
- **Actual:** the tripwire is a **Δρ-width** check. Its pass condition is the SC2 effect itself, so
  a live-but-unresponsive-in-ρ net fails it identically to a dead-input net. The liveness question
  is answered on a different channel (the shift marginals), where the answer is unambiguously
  "alive".
- **Why it matters:** it changes what a reader may conclude from the red gate. Under the planned
  assumption the correct response is "find the bug"; under the measured one there is no bug to
  find, and the open question becomes whether SC1g is measuring the right thing. It also means the
  wave-7 gate, as currently wired, cannot be cleared by any net whose ρ width does not track λ —
  which is the very hypothesis SC2 was built to test.

**2. "Early stopping before epoch 80 indicates a training stall."**
- **Found during:** Task 3, Step 2
- **Planned:** R11 treats an early-stopping epoch below 80 as a stall signature warranting a ladder
  reduction.
- **Actual:** validation risk fell 15.645 → 6.049 (−61 %) by epoch 22 and then diverged while
  training risk kept falling. That is textbook overfitting on a 40 000-sample pool, not a stall;
  the signature's epoch-count leg fires on it regardless.
- **Why it matters:** taken literally the signature would have triggered a pre-registration edit
  and ~72 min of re-spend in response to a healthy fit.

## Issues Encountered

- Julia block-buffers stdout when redirected to a file, so neither long run could be monitored via
  its log. Progress was tracked instead by polling shard files on disk and, during training, the
  `loss_per_epoch.csv` NeuralEstimators writes into its `mktempdir`. No behaviour was changed.
- Reading the net through a bare `julia -e` (outside the trainer's include chain) emits JLD2
  "type … does not exist in workspace; reconstructing" warnings for `ZScoreTransform` and
  `P11BoundedThetaTransform`. Cosmetic and expected — the types are spike-local. Loading through
  `load_p11_npe` / `load_npe` inside the proper include chain is clean, which is how the tripwire
  and the diagnostic both read it.

## Deferred Issues

- **`spike/test/runtests.jl` still exits 1** on the pre-existing `SPEEDUP_GATE` (NPE-03) blocker
  (92.5× / 84.0× against a `> 100×` bar). Not this plan's to fix; not lowered. The tripwire was
  not added to the aggregate runner by this plan.

## Known Stubs

**None.** Both modules are complete and both ran end-to-end at full pre-registered scale. The one
incomplete thing is not a stub but a red gate, reported above as a blocker.

## Threat Flags

None. No network endpoint, auth path or schema at a trust boundary was added. Registered threat
handling:

| Threat ID | Handling | Evidence |
|---|---|---|
| T-11-27 (λ drawn independently of the shift) | Mitigated **and independently confirmed on the realized 50 000-column pool**, not only on fresh draws | 0/50 000 invariant violations; `cor(λ,\|dx\|)` = 0.6042, `cor(λ,\|dy\|)` = 0.6114; row 129 varies over 39 968 values |
| T-11-28 (unrecorded realized training image-size distribution) | Mitigated | `training_imsize_provenance` persisted with realized counts incl. the 1376×1028 anchor; asserted non-empty |
| T-11-29 (silent downgrade of the image-size arm) | Not triggered | 56.37 min vs a 150-min ceiling; guard never fired; persisted provenance is a 4-size mixture |
| T-11-30 (torn artifact) | Mitigated | `.tmp` → reopen-and-assert → `mv(force=true)` for both the 5 shards and the net; reopen check passed |
| T-11-31 (whole-object persistence) | Accepted as pre-registered | header and this summary both state the net is not shipped and not distributable |
| T-11-32 (research net mistaken for the shipped net) | Mitigated | payload carries the 4 `P11_RESEARCH_NET_DEVIATIONS`, `research_lane = true`, `shipped = false`; lives under `spike/npe/`, not `artifacts/` |
| T-11-SC (package-manager installs) | None attempted | `git diff --stat -- spike/Project.toml spike/Manifest.toml` empty |

## STATE.md — proposed BLOCKER text for the orchestrator (NOT written by this executor)

> **BLOCKER — 2026-07-27 — Phase 11 plan 11-07: the SC1g λ-ablation tripwire FAILED, and the
> prescribed diagnosis came back GREEN.** Datagen (50 000 pairs, F5 mixture) completed in
> **56.37 min** against `P11_DATAGEN_WALLCLOCK_CEILING_MIN = 150` — the ceiling was **not**
> exceeded and the image-size arm was **not** downgraded. Training completed CPU-only in
> **16.08 min** at `d_in = 129, D = 8`; validation risk 15.6452 → **6.0490** (best epoch 22), early
> stopping at **epoch 63** of 300. `spike/test/test_lambda_ablation.jl` then **failed** (exit 1,
> 9 pass / 8 fail): per-dataset Δρ posterior-SD ratios **1.109 / 0.987 / 0.936** and 90 % HDI ratios
> **1.122 / 0.984 / 0.940**, against the pre-registered `P11_LAMBDA_ABLATION_FACTOR = 2.502`.
>
> The plan's prescribed diagnosis was run in order and **all three legs are green**: (a) on the
> **realized** pool, 0 of 50 000 samples violate `|shift| <= λ` and `cor(λ,|shift|)` = 0.604/0.611
> (bar > 0.3), with `cor(λ,|chromatic_eps|)` = −0.007; (b) the training inputs were **129 rows**
> with the 129th varying over **39 968** of 40 000 columns; (c) `augment_input` ran **after**
> standardization — row 129 equals `encode_lambda(λ)` exactly and rows 1:128 equal the standardized
> summary exactly. A fourth check settles the ambiguity those three leave: the **shift** marginals,
> whose λ response is analytically known, widen by **4.84 / 6.99 / 8.05** (dx) and
> **3.98 / 6.52 / 9.05** (dy) against an ideal prior-SD ratio of 12.0. **The λ conditioning is
> ALIVE and used; it is the ρ_true width that does not track λ (ratios 1.20 / 1.06 / 0.98).**
>
> Consequence: **wave 7 (plan 11-08, the SC2 ladder) stays GATED** and no ladder was run. Nothing
> was relaxed: `P11_LAMBDA_ABLATION_FACTOR` is unchanged, `LAMBDA_MAX` is unchanged (the R11
> "reduce `LAMBDA_MAX`" remedy was **not** executed — its premise is falsified, `LAMBDA_MAX` is
> Tier-1 pre-registration that may not be edited, and a narrower span makes the bar harder), no
> capacity was raised, no seed re-drawn, no re-run performed. `spike/validation/p11_consts.jl` and
> `src/` are byte-unchanged.
>
> **Decision required (not this executor's to take):** SC1g's statistic is Δρ width, i.e. the same
> quantity SC2 measures, so as authored it cannot separate "the conditioning is dead" from "the
> effect is null" — the separation its header claims to provide. Options are (i) re-scope SC1g onto
> a channel with a known response (the shift marginals, per diagnostic (d)) as a corrected liveness
> tripwire, keeping Δρ as an SC2 question; (ii) read the flat ρ response as the SC2 result and
> proceed to the ladder to measure it at pre-registered scale; or (iii) spend
> `P11_ITERATION_ALLOWANCE` (still **1 of 1 unspent**). Artifact:
> `spike/npe/p11_research_npe.jld2` (`generated = 2026-07-27T13:49:47.473Z`, 4 933 362 bytes,
> committed). See `11-07-SUMMARY.md`.

## User Setup Required

None.

## Next Phase Readiness

- **BLOCKED — plan 11-08 (the SC2 ladder).** Its entry gate (SC1g) is red. Do not run it until the
  decision above is taken.
- **Available for reuse, at no further compute cost:** the 50 000-pair pool
  (`spike/data/cache/p11/c12f8ab2…`, 54 MB, resumable, content-hash-keyed on the generating config
  including `lambda_max`) and the trained net. Any option that keeps `LAMBDA_MIN`/`LAMBDA_MAX` and
  the F5 mixture can reuse both **without regenerating**; changing `lambda_max` invalidates the
  pool by content hash and costs a fresh ~56 min.
- **Carried forward UNSPENT:** `P11_ITERATION_ALLOWANCE = 1` (0 of 1 spent). This plan spent no
  iteration — it relaxed nothing and re-ran nothing.
- **For the report (plan 11-11), three items are now measured and quotable:** the λ conditioning is
  demonstrably alive on the shift channel with named ratios; the ρ_true width response is flat at
  fixture scale; and the realized F5 training mixture is recorded with the real-data anchor at
  12 289 of 50 000.
- **Left red for someone else:** the pre-existing `SPEEDUP_GATE` blocker in `spike/test/runtests.jl`.

## Self-Check: PASSED

- Created file exists on disk: `spike/npe/p11_research_npe.jld2` (4 933 362 bytes) — **FOUND**.
- Modified file exists: `spike/test/test_lambda_ablation.jl` — **FOUND**.
- Commits exist in `git log`: `59d3ed8` — **FOUND**; `fdf32ed` — **FOUND**. Prior-executor commits
  `efb30f3`, `8fcc9b8` — **FOUND** (reachable from HEAD).
- `git diff --stat -- src/ spike/Project.toml spike/Manifest.toml` — **empty**, as required.
- `.planning/STATE.md` and `.planning/ROADMAP.md` untouched, per the parallel-execution contract.
- **Plan status recorded honestly as BLOCKED**, not complete: Task 3's `<verify>` exits 1 and this
  summary says so in its first line, its Status section and its verification table.

---
*Phase: 11-registration-and-chromatic-uncertainty-as-latent*
*Completed: 2026-07-27 (BLOCKED at the SC1g gate)*
