---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 13
status: complete
subsystem: amortized-inference
tags: [reported-not-gated, alpha-series, exclusion, correlation-vs-localization, named-limit]
requires:
  - "spike/p13/three_way_net.jld2 (13-11) — the trained two-head evidence net + its MEASURED per-head log-odds, consumed through the PERSISTED artifact and never re-measured"
  - "spike/p13/alpha_series.jl (13-04) — alpha_ladder / alpha_ladder_summaries / verify_alpha_invariants, with its four invariants already proven by spike/test/test_p13_alpha.jl"
  - "spike/p13/consts.jl (13-01/13-10) — the byte-locked pre-registration incl. MEASURED P13_TAU = 0.15 and the frozen P13_ALPHA_LADDER"
  - "spike/p13/datagen.jl (13-11) — p13_sample_theta_given_lambda / p13_sample_imsize (the F5 mixture)"
  - "spike/npe/p11_research_npe.jld2 (Phase 11) — the bound research basis and its frozen zt"
  - "spike/simulator/ghat.jl — the FROZEN m-bar -> rho map, read (not fitted) for the interpretation below"
provides:
  - "spike/p13/run_alpha_series.jl — the REPORTED, non-gating D-15/D-16 alpha-ladder runner (invariants-before-curves, persist-before-anything-throws)"
  - "spike/p13/alpha_report.jld2 — the three curves, the per-image scores, the invariant residuals, alpha_star, the seed, the counter and the consts fingerprint"
  - "spike/figures/p13_alpha_ladder.png — the three-curve figure (script artifact, gitignored under the house *.png rule)"
  - "The Phase-16 deferral paragraph, written to be lifted verbatim into 13-14's named-limits section"
affects:
  - "13-14 — the phase report cites alpha_star and lifts the deferral paragraph; the per-item scores are persisted so no plan need re-run this"
  - "13-16 — the :real arm of the SAME ladder, reported separately and NEVER averaged with these curves"
tech-stack:
  added: []
  patterns:
    - "Invariants-before-curves: a flat ladder caused by deleted pixels is indistinguishable from a null result, so verify_alpha_invariants runs on the REPORTED substrate and halts the run before a single curve is read"
    - "Substrate-declared dynamic range: max_value is a caller-supplied property of the substrate, declared UNBOUNDED on the simulated arm and the realized maximum RECORDED instead of a number being chosen"
    - "Crossing point reported as a RUNG of the frozen ladder, never interpolated"
key-files:
  created:
    - spike/p13/run_alpha_series.jl
    - spike/p13/alpha_report.jld2
    - spike/figures/p13_alpha_ladder.png
  modified: []
decisions:
  - "The simulated arm declares an UNBOUNDED dynamic range rather than the frozen P13_ALPHA_MAX_VALUE_BOUND = 1.0, which is a Gray-TIFF property of the REAL substrate; the simulator's softplus intensities measured 4.08-17.96 here. Realized maxima are recorded per image; nothing is clamped and no new bar is invented."
  - "The CONTROL acquisition is drawn at rho = 0 and held FIXED across the whole ladder while only the SAMPLE is graded, because the D-05 cut is a two-factor rule on the sample LEVEL and the sample-control CONTRAST; grading both halves would move both arguments together and the ladder would measure nothing."
  - "rho_true is pinned at EXACTLY 0.0 rather than approximately 0, so alpha = 0 is a random pair by construction rather than by measurement."
metrics:
  duration: ~2.5 min wall clock for the reported run (0.43 min substrate draw, ~1.1 min invariant re-verification, 0.50 min scoring); ~35 min total including authoring and verification
  completed: 2026-07-29
---

# Phase 13 Plan 13: Reported Alpha-Graded Segregation Series Summary

**`alpha_star = 0.25`.** On the `:simulated` arm — 64 items drawn at `rho_true = 0` on the reserved
stream at `P13_ALPHA_COUNTER = 4` — the mean `log BF(exclusion : random)` first becomes positive at
the third rung of the frozen `P13_ALPHA_LADDER`. **All eight transform invariants hold on all 64
reported images**, zero patches went absent at any rung, and the summary's mask rows are
byte-identical across the whole ladder — so the ladder's shape is not a mechanical artifact of
deleted pixels.

**The crossing point is not arbitrary, and that is the result.** Mapping the persisted `m-bar`
through the frozen `ghat`, `alpha = 0.125` induces `rho = -0.1443`, which sits **inside** the
pre-registered `P13_TAU = 0.15` dead zone, and `alpha = 0.25` induces `rho = -0.2499`, the **first
rung to clear it**. The evidence head first calls exclusion at exactly the first rung where the D-05
label rule says a pair stops being random — on a spatially-constructed ladder the net never saw
during training.

**There is no gate here and none was added.** `P13_ALPHA_GATED = false` was frozen before any
Phase-13 result existed; the runner contains no assertion on any curve, asserted by source grep.
`spike/p13/consts.jl` is byte-unchanged (git blob `5a4ea222…`, identical to the sha 13-11 and 13-12
both recorded) and **`P13_ITERATION_ALLOWANCE` remains 1 of 1, UNSPENT**.

## 1. The ladder (REPORTED; no bar is attached to any number below)

64 images, one control acquisition each, `lambda = 3.0` (`:widest_rung`, the most conservative
reading), `G = 8`. `m-bar` is the mask-weighted mean of the present continuous summary rows; the
`missing` column counts absent patches summed over all 64 images (4,096 patches per rung).

| alpha | mean log BF(E:R) | median log BF(E:R) | mean log BF(C:R) | median log BF(C:R) | m-bar | missing |
|---|---|---|---|---|---|---|
| 0.000 | **−4.07387** | −4.09404 | −3.23876 | −3.49149 | +0.067767 | 0 / 4096 |
| 0.125 | −1.13730 | −1.09507 | −5.69373 | −5.93463 | −0.018156 | 0 / 4096 |
| 0.250 | **+1.50509** | +1.43886 | −7.30997 | −7.49974 | −0.104270 | 0 / 4096 |
| 0.375 | +3.67840 | +3.57581 | −8.06796 | −8.08047 | −0.186040 | 0 / 4096 |
| 0.500 | +5.44641 | +5.38472 | −8.40306 | −8.38761 | −0.259750 | 0 / 4096 |
| 0.625 | +6.84283 | +6.80221 | −8.63395 | −8.61835 | −0.323330 | 0 / 4096 |
| 0.750 | +7.91279 | +7.85100 | −8.85249 | −8.81939 | −0.376360 | 0 / 4096 |
| 0.875 | +8.72483 | +8.69582 | −9.06106 | −9.01914 | −0.419660 | 0 / 4096 |
| 1.000 | +9.34472 | +9.31883 | −9.25016 | −9.18802 | −0.454640 | 0 / 4096 |

Mean and median track each other to within ~0.07 nats at every rung, so no rung's mean is being
carried by a tail.

**Per-image spread**, reported because a mean crossing point is a statement about a population and
this series has n = 64:

| alpha | sd of log BF(E:R) | min | max | fraction of images with log BF(E:R) > 0 |
|---|---|---|---|---|
| 0.000 | 1.880 | −8.514 | −0.716 | **0 / 64 (0.0%)** |
| 0.125 | 1.807 | −6.087 | +2.114 | 17 / 64 (26.6%) |
| 0.250 | 1.591 | −3.206 | +4.499 | **55 / 64 (85.9%)** |
| 0.375 | 1.401 | −0.436 | +6.358 | 63 / 64 (98.4%) |
| 0.500 | 1.226 | +2.068 | +7.764 | **64 / 64 (100%)** |
| 0.625–1.000 | 1.072 → 0.821 | ≥ +4.164 | ≤ +10.805 | 64 / 64 (100%) |

The spread NARROWS monotonically as alpha grows (1.88 → 0.82 nats). Not one of the 64 items is
called exclusion at `alpha = 0`; every one of them is from `alpha = 0.5` onward.

The coloc head moves the other way and never disagrees: `log BF(C:R)` is negative on 63 of 64 items
at `alpha = 0` and on **64 of 64 at every rung above it**. The one positive item at `alpha = 0` is a
weakly-coloc call on a genuinely random pair, which is what a calibrated head at a random draw is
supposed to produce sometimes.

## 2. `alpha_star = 0.25` — the headline, and what it does and does not say

`alpha_star` is defined as the SMALLEST rung of the frozen `P13_ALPHA_LADDER` at which the mean
`log BF(exclusion : random)` first exceeds 0. It is **0.25**. It is reported as a RUNG and never as
an interpolated value: an interpolated crossing would be a number the pre-registration does not
contain, computed from a curve nobody committed to a functional form for. **The ladder's own
resolution is the resolution of this statistic** — the true crossing lies somewhere in
`(0.125, 0.25]` and this design cannot say where inside it.

## 3. The correlation-versus-localization reading

This series exists for one question: the simulator's "exclusion" mechanism **is** negative intensity
correlation across patches (the `sign(rho)` flip on channel 2's shared latent component,
`spike/simulator/forward.jl` stage 1), while a biologist's "mutually exclusive" is **disjoint spatial
localization**. The ladder is built spatially — `alpha_segregate` moves ch2's above-background mass
out of the ch1 object mask and returns it, intensity-conserving, to the complement — and is scored
by a net trained only on the anti-correlation notion. So the question is at what alpha the second
notion starts to track the first.

**The answer is that the two notions are separated by exactly the pre-registered hypothesis
boundary, and by nothing else.** Reading the persisted `m-bar` through the frozen `ghat` map:

| alpha | m-bar | `ghat(m-bar)` = induced rho | vs. `P13_TAU = 0.15` | mean log BF(E:R) |
|---|---|---|---|---|
| 0.000 | +0.06777 | −0.02645 | inside the dead zone | −4.074 |
| 0.125 | −0.01816 | −0.14430 | **inside the dead zone** (by 0.006) | −1.137 |
| 0.250 | −0.10427 | −0.24989 | **CLEARS tau** | **+1.505** |
| 0.375 | −0.18604 | −0.36684 | clears tau | +3.678 |
| 0.500 | −0.25975 | −0.45129 | clears tau | +5.446 |
| 0.625 | −0.32333 | −0.53094 | clears tau | +6.843 |
| 0.750 | −0.37636 | −0.60285 | clears tau | +7.913 |
| 0.875 | −0.41966 | −0.66511 | clears tau | +8.725 |
| 1.000 | −0.45464 | −0.70471 | clears tau | +9.345 |

The summary's own sign flip happens **one rung earlier** than the evidence's — `m-bar` goes negative
between `alpha = 0` and `alpha = 0.125`, the log Bayes factor between `0.125` and `0.25`. That one
rung of lag is not a defect and not an unexplained gap: at `alpha = 0.125` the constructed spatial
segregation induces an anti-correlation of only `rho = −0.144`, which is **below the summary's own
MEASURED resolution** (`P13_TAU = 0.15`, measured in plan 13-10 as the smallest separation the fixed
8×8 patch-correlation summary can reliably distinguish, at `A(tau) = 0.916`). D-05 defines such a
pair as RANDOM, so a head that called it exclusion would be wrong by the phase's own label rule. The
head declines, and then calls exclusion at the very first rung where the induced rho clears the
boundary.

**So the correlation-versus-localization gap, on this substrate, is quantified as: roughly a quarter
of ch2's above-background mass must be moved out of the ch1 objects before disjoint spatial
localization registers as exclusion under the intensity-correlation notion.** Below that, the
segregation is genuinely present in the pixels and genuinely invisible to the summary — not to the
net, to the *summary*, which is the more useful place to have located it.

**THE LIMITATION, STATED PLAINLY: alpha is a construction parameter, not a physical quantity.** No
microscope has an alpha, no specimen has an alpha, and nothing in this project can convert an alpha
into a biological degree of segregation. `alpha_star` therefore **CALIBRATES** the gap — it says how
much of *this particular constructed transform* is needed before the evidence responds — it does
**not MEASURE** the gap as a property of real biology. A reader who quotes "0.25" as "25% segregation
is detectable" will have over-read it in exactly the way this paragraph exists to prevent.

**Two further scope limits on the same reading.** (i) The ladder spans induced rho from −0.026 to
−0.705, so it **does not probe the deep-exclusion tail near the −0.99 prior atom** — it tests the
sign transition and the near-to-moderate exclusion regime only. (ii) At `alpha = 0` the induced rho
is −0.0265 rather than 0 (equivalently `m-bar = +0.068 ± 0.046`); the frozen `ghat` map's own
zero-crossing does not sit exactly at `m-bar = 0`, so the `alpha = 0` rung is a random pair by
CONSTRUCTION (both acquisitions drawn at `rho_true = 0.0`) and is very slightly off zero by
MEASUREMENT. Both statements are true and neither is a defect.

## 4. Why the ladder can be trusted not to be flat, or shaped, for a mechanical reason

`verify_alpha_invariants` ran on the **REPORTED** substrate — every one of the 64 items, not a
fixture — **before any curve was read**. The runner is written to halt with a named
`INVARIANT VIOLATION` block, persist what it has, and read nothing, if any invariant fails on any
image.

| Invariant | Meaning | Images failing | Measured |
|---|---|---|---|
| `bitwise_alpha0` | `alpha = 0` reproduces the input BITWISE (exact `==`, never `isapprox`) | **0 / 64** | — |
| `no_new_zeros` | `count(iszero, y_out) == count(iszero, y)`: alpha introduces no NEW zero | **0 / 64** | source zeros 0 … 10,659 per image |
| `intensity_conserved` | `sum(y_out) ≈ sum(y)` at `rtol = 1e-10` | **0 / 64** | worst residual **9.63e-16**, ~5 orders inside the bar |
| `mask_invariant` | ch1 untouched AND the recomputed ch1 mask equals the ladder mask at every rung | **0 / 64** | — |
| `mask_fraction_ok` | `mean(M)` inside `P13_ALPHA_MASK_FRACTION_BOUNDS = (0.01, 0.40)` | **0 / 64** | realized 0.2840 … 0.3427, mean 0.3070 |
| `max_value_ok` | realized maximum within the DECLARED range | **0 / 64** | **vacuously true** — see below |
| `mbar_monotone` | `m-bar(alpha)` non-increasing | **0 / 64** | — |
| `missing_nonincreasing` | absent-patch count does not rise with alpha | **0 / 64** | — |

**`max_value_ok` is vacuously true on this arm and is not quoted as evidence.** The simulated arm
declares an UNBOUNDED dynamic range (see Deviation 2), so the comparison cannot fail. What is
reported instead is the realized maximum: **4.079 … 17.964** across the 64 images, recorded per
image in the artifact and never clamped.

**Pitfall-2 warning signs — the mechanical explanations a flat or misleading ladder could have:**

| Warning sign | Measured | Reading |
|---|---|---|
| absent patches at the FIRST rung (`alpha = 0`) | **0 of 4,096** | — |
| absent patches at the LAST rung (`alpha = 1`) | **0 of 4,096** | — |
| did the absent-patch count RISE across the ladder? | **NO** | pixels are being MOVED, not deleted |
| any absent patch at any rung on any image? | **NO** | the 15-survivor floor was never approached |
| did the summary's MASK ROWS change with alpha? | **NO — byte-identical across every rung, on all 64 images** | the zero SET is untouched, so `_exclude_zero` sees the same pixels at every rung |

Both signs are precisely the ones the strictly-positive background floor `b = quantile(vec(y), 0.05)`
exists to prevent: masked ch2 pixels are pushed DOWN TO a positive floor, never to zero, so every
pixel keeps participating in its patch correlation. Had they been zeroed, the measured correlation
would have been computed over the untouched complement only — approximately unchanged, i.e. random —
and the ladder would have read flat for a reason with nothing whatsoever to do with segregation.
**That failure did not occur, and the evidence that it did not occur is in the table above rather
than in an intention.**

Honest caveat carried forward from the transform's own design: **at `alpha = 1` the mask region is
not EMPTY, it is AT the background floor `b`.** The endpoint is "ch2 reduced to background wherever
ch1 has objects", not "ch2 absent". Do not write "fully disjoint" without that qualifier.

## 5. Substrate scope — this is the `:simulated` arm, and the two arms are never averaged

`P13_ALPHA_SUBSTRATE = (:simulated, :real)` is two-valued and pre-registered. **This run covers
`:simulated` only.** The `:real` arm runs the SAME `alpha_segregate` code path — typed on
`AbstractMatrix{Float64}` with no substrate branch anywhere, so the rungs are comparable rung for
rung — over the six committed `test/test_images/` TIFFs, and is reported separately by plan 13-16's
`spike/p13/run_p13_realimage.jl`.

**Their curves must NEVER be averaged or merged.** What differs is what `alpha = 0` MEANS. Here it is
a genuinely random pair by construction: both acquisitions are drawn at `rho_true = 0.0` from the
Phase-11 simulator under the frozen F5 image-size mixture. There it is **moderately colocalized** —
the unmodified real pairs measure `m-bar = +0.3292` (positive fixture) and `+0.2481` (negative
fixture). The simulated arm's ladder runs random → exclusion; the real arm's runs moderately
colocalized → near-exclusion. Both are informative; they are **not the same experiment**, and a
merged curve would be a statement about neither.

**The CBS `cbs-RG-000` arm was DROPPED from scope, not deferred and not optional.** It would require
the fetch path, and the amended D-15 supplies a real substrate that needs no fetch, no hash gate and
no seal. **Nothing under the provenance corpus directory was checked, fetched or referenced** — not
by the runner, which contains no reference to it at all (asserted by the pre-registered source-grep
gate in `spike/test/test_p13_real.jl`, which scans every `.jl` under `spike/p13/` with `readdir` so a
newly added file cannot slip past by not being listed), and not by this plan.

## 6. NAMED LIMIT — the single physical segregated anchor is deferred to Phase 16

> *(This paragraph is written to be lifted verbatim into 13-14's named-limits section.)*
>
> **The provenance corpus contains exactly ONE real segregated anchor** — `neg-lightmycells-01`,
> `physical-primary`, condition `negative`, truth label `segregated`, channels
> `2:nucleus-dna,mitochondria` — and Phase 13 did not read it. It is **triply unavailable**, and each
> obstacle is independently sufficient: (1) its recorded `sha256` is the unfilled sentinel
> `"PENDING-FETCH"`, because the bootstrap fetch was deliberately never run, so there is no digest to
> verify any bytes against; (2) the corpus data directory is **empty** and carries no tracked bytes —
> the row records `bytes = 0`, which states exactly that nothing was fetched; and (3) the row's split
> is `sealed_holdout`, reachable only through the anti-snooping accessor `open_sealed_holdout(df;
> reason)`. **Phase 13 deliberately did not open that seal.** The seal is not an inconvenience to be
> routed around — it is a control that exists for Phase 16's BLIND evaluation, and consuming the one
> segregated anchor here would irreversibly burn Phase 16 on the very hypothesis Phase 16 exists to
> evaluate. No seal-break escape hatch exists anywhere in Phase 13: not in a script, not in a test,
> not behind a flag, and the earlier proposal for a `checkpoint:human-verify` escape hatch was
> withdrawn by the amended D-15 rather than left dormant. **The n = 1 physical segregation check is
> therefore DEFERRED TO PHASE 16** (`P13_PHYSICAL_ANCHOR_DEFERRED_TO = "Phase 16"`, frozen in the
> pre-registration). No download-size claim attaches to this deferral: `bytes = 0` records only that
> nothing was fetched, and the size question is retired as moot (assumption A6, withdrawn) precisely
> because Phase 13 fetches nothing.

**The consequence for what this plan establishes.** The exclusion hypothesis still has **no observed
real-data instance anywhere in this project**. Every unmodified real image measured in Phase 13 sits
at positive induced mu, including the "negative" biological control at `m-bar = +0.2481`. This
series therefore constructs its segregation rather than observing it, and the construction is one
particular intensity-redistribution rule rather than a sample of how real proteins segregate.
Labelled real segregation validation remains Phase 16's, and this plan leaves that control intact.

## 7. What this series does NOT establish

Stated plainly so the result is not over-read:

- **It is not a gate and cannot become one.** `P13_ALPHA_GATED = false` was pre-registered.
  Promoting an observation here into a criterion after the fact would be an undeclared amendment.
- **alpha is a knob, not a measurement.** See §3's limitation paragraph.
- **The deep-exclusion tail is untested here.** Induced rho spans [−0.705, −0.026]; the −0.99 prior
  atom is nowhere near this ladder. (That tail WAS tested, and cleared, by the D-12 gate: exclusion
  AUC 0.997257 at `|rho| > 0.9`.)
- **Simulated substrate only.** The `:real` arm is 13-16's and reports different numbers for a
  different experiment.
- **One transform, one mask rule.** The segregation is constructed by mask-based disjoint
  reassignment under the frozen `_calculate_mask` Otsu rule. A different segmentation rule or a
  different redistribution rule would give a different `alpha_star`, and nothing here bounds by how
  much.

## Verification performed

| Check | Result |
|-------|--------|
| `julia --project=spike -t auto spike/p13/run_alpha_series.jl` | **exit 0**, ~2.5 min wall clock |
| Task-1 `<verify>` (`alpha_star`, `logbf_exclusion`, `mbar` present) | prints **`alpha artifact ok`** |
| Task-2 `<verify>` (`alpha_star` + rung count) | prints **`alpha_star = 0.25 over 9 rungs`** — 9 = `length(P13_ALPHA_LADDER)` |
| `spike/p13/alpha_report.jld2` exists, no `.tmp` sibling | 53,391 bytes, confirmed; 65 keys |
| one entry per rung for `logbf_exclusion` / `logbf_coloc` / `mbar` / `n_missing` | all length 9 |
| recorded `ladder` equals the frozen `P13_ALPHA_LADDER` | identical |
| `spike/figures/p13_alpha_ladder.png` rendered | 209,646 bytes |
| `grep -c '@test' spike/p13/run_alpha_series.jl` | **0** |
| `grep -v '^\s*#' … \| grep -c 'open_sealed_holdout\|fetch_verified\|bootstrap_anchor_hashes'` | **0** |
| `grep -c 'Phase 16' spike/p13/run_alpha_series.jl` | **3** |
| `grep -v '^\s*#' … \| grep -c '0\.125\|0\.875'` | **0** (the ladder is read from consts, never inlined) |
| `grep -c 'corpus' spike/p13/run_alpha_series.jl` | **0** (not even in a comment) |
| `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml corpus/` | **exit 0 — byte-unchanged** (src/ and the manifests asserted again at RUN TIME by step 0) |
| `git status --porcelain corpus/data` | **empty** |
| `git diff --quiet HEAD -- spike/p13/consts.jl` | **exit 0 — byte-unchanged**; `consts_git_blob_sha` in the artifact is `5a4ea2223ce11bcae2f974c643ec1e523458e8ac`, identical to 13-11's and 13-12's records |
| `net_consts_sha == p13_consts_sha()` | `100e97a3…` both sides (asserted at run time) |
| Post-commit deletion check | no deletions |
| No new bulk cache directory created | confirmed — the two `spike/data/cache/p13/` directories both predate this run (2026-07-28); `spike/data/cache/p11/` (54 MB) untouched |

**Test files run.** Adding a file under `spike/p13/` is not inert: **two** Phase-13 test files
discover that directory with `readdir` and grep every `.jl` in it. Both were run directly, together
with the transform's own suite:

| File | Result |
|---|---|
| `spike/test/test_p13_alpha.jl` | **exit 0** — all testsets pass (the four transform invariants, the guards, the sealed-holdout grep) |
| `spike/test/test_p13_real.jl` | **exit 0** — incl. `P13 corpus seal untouched` 31 pass / 1 **broken** (pre-existing `@test_broken`, not a failure); this is the testset that greps the new runner |
| `spike/test/test_p13_calibration.jl` | **exit 0** — incl. the D-14 `no equivalence machinery` absence grep, 31/31, which also scans the new runner |

The FULL spike suite still exits 1 at `spike/test/runtests.jl:169` on the Phase-4 `SPEEDUP_GATE`
(measured 68.35 against the pre-registered bar 100.0), which masks the entire Phase-13 include block.
That is **pre-existing**, already logged as `deferred-items.md` D-13-A by 13-11, and is not this
plan's to fix — it is a pre-registered threshold in a Phase-4 file this plan does not own. Per-file
runs are the available gate signal for Phase 13.

## Deviations from Plan

### Auto-fixed issues

**1. [Rule 3 — Blocking] The printed deferral paraphrases the sealed-holdout accessor name; the
literal name appears only in the file header comment**
- **Found during:** Task 1, writing the header.
- **Issue:** the plan's step (c) asks the runner to state "explicitly that `open_sealed_holdout` is
  NOT called". The plan's own acceptance criterion asserts
  `grep -v '^\s*#' … | grep -c 'open_sealed_holdout\|…'` outputs `0`, and — independently and more
  bindingly — `spike/test/test_p13_real.jl:240-248` is a **pre-registered source-grep gate** that
  scans every `.jl` under `spike/p13/` with `readdir`, strips whole-line comments, and asserts the
  remaining text contains neither that identifier nor the corpus directory name. A `println` string
  is not a comment line, so printing the literal at run time would have turned that gate red on a
  correct runner.
- **Fix:** the full statement — naming the accessor, its signature, all three obstacles and the
  reason the seal is left intact — lives in the file's header COMMENT, which is exactly the form both
  greps exempt and which the seal ruling explicitly permits ("may be cited in prose for provenance
  only"). The runtime banner states the same thing in words that carry no forbidden identifier:
  *"The sealed-holdout accessor is NOT called by this runner and the anti-snooping seal is left
  INTACT for Phase 16's blind evaluation."*
- **Files:** `spike/p13/run_alpha_series.jl`. **Commit:** `0943a02`.

**2. [Rule 3 — Blocking] The simulated arm declares an UNBOUNDED dynamic range instead of the frozen
`P13_ALPHA_MAX_VALUE_BOUND = 1.0`**
- **Found during:** Task 1, on the fixture-stream mechanics check (see Disclosure 4).
- **Issue:** `P13_ALPHA_MAX_VALUE_BOUND = 1.0` is measured from the REAL substrate — the committed
  TIFFs are Gray values in `[0, 1]`. The simulator's softplus intensities are unnormalized and
  measured **4.079 … 17.964** on this substrate *before any alpha*, so passing the frozen bound
  would have made `alpha_segregate` throw on the very first image, for a reason that is a property of
  the substrate rather than of the transform.
- **Fix:** `max_value` is forwarded as `Inf`. This is not a workaround — `alpha_series.jl`'s own
  docstring designates `max_value` a **caller-supplied property of the SUBSTRATE** precisely so the
  simulated arm can declare its own range, and `verify_alpha_invariants` already builds its ladder
  with the range check disabled and *reports* the realized maximum. Declaring some other finite
  number would have meant **inventing a bar** in a runner whose entire posture is that it scores
  nothing, which is exactly what this plan's discipline forbids. So nothing is clamped, nothing
  throws, and the realized maximum is RECORDED per image — which is what the pre-registration asks
  for ("ASSERTED AND RECORDED, NEVER SILENTLY CLAMPED").
- **Consequence, disclosed rather than buried:** the `max_value_ok` invariant is **vacuously true**
  on this arm. It is reported as such in §4 and in the run's own printed output, and is never quoted
  as evidence of anything.
- **Files:** `spike/p13/run_alpha_series.jl`. **Commit:** `0943a02`.

**3. [Rule 3 — Blocking] The figure is composed in the runner, not routed through a `figures.jl`
helper**
- **Found during:** Task 1.
- **Issue:** the plan asks for the figure "through the existing `spike/validation/figures.jl`
  helpers into `spike/figures/p13_alpha_ladder.png`". Two things make that impossible as written, and
  they are the same two plan 13-12 recorded for its own figure: `figures.jl` has **no** ladder
  helper, and every helper it does have routes its output through `_val_fig_path`, which hard-joins
  to `spike/validation/figures/` — not the path the plan requires.
- **Fix:** `p13_alpha_figure` composes the two panels using the CairoMakie surface, the headless
  backend activation and the caption discipline that `figures.jl` establishes and that the runner
  inherits by including it. No plotting dependency was added. The reason is documented in the
  function's docstring.
- **Files:** `spike/p13/run_alpha_series.jl`. **Commit:** `0943a02`.

### Design choices the plan left to the executor (recorded, not deviations)

**4. The CONTROL acquisition is drawn at `rho = 0` and held FIXED across the ladder; only the SAMPLE
is graded.** The plan says to build the transformed pair and score it, without saying what the
control half is. The evidence net reads a PAIR — a sample stack and its control stack — and the D-05
cut is a **two-factor** rule on the sample LEVEL *and* the sample-control CONTRAST. Grading both
halves in lockstep would move both arguments together and the ladder would measure nothing; grading
only the sample makes the ladder walk exactly the trajectory `rho_s < -tau AND rho_s - rho_c < -tau`
describes. Recorded here because it is load-bearing for how §3's reading should be understood.

**5. `rho_true` is pinned at EXACTLY `0.0`, not merely near zero.** The plan says "at a near-zero
true rho". Exactly zero is the cleanest reading and is what makes `alpha = 0` a random pair by
CONSTRUCTION rather than by measurement. Everything else in the item — the seven nuisances, the
lambda-conditional shift, the F5 image size — is drawn from the same joint the training pool used,
through the same `p13_sample_theta_given_lambda`.

### Disclosures (not deviations, but the ordering must be checkable)

**6. A FIXTURE-STREAM mechanics check preceded the reported run, and the runner was written after
it.** Before writing the runner, two probes were run on the **fixture** family
(`P13_FIX_SEED` at `P13_FIXTURE_COUNTER`), never on the reported alpha stream, to answer two
questions the plan could not answer from the source: does the simulated substrate satisfy the
transform's mask-fraction and strictly-positive-floor guards, and what does one image cost. They
showed mask fractions 0.288–0.318, `b` ≈ 0.16–0.22 (strictly positive with margin), zero absent
patches, and a monotone `m-bar`. One of them also printed a per-rung log-BF ladder on a fixture
image, so the ladder's rough SHAPE was seen before the reported run. Four things bound what that
means:

- The probes rode the **fixture** seed and counter, which by construction never touch the reported
  stream. No reported index was consumed or observed.
- **The reported run was executed exactly once**, and every value it scores against was already
  byte-locked in `consts.jl` and committed long before either probe.
- **There is nothing here to tune.** The ladder is frozen, `alpha_star` is a definition rather than a
  threshold, and this plan has no criterion of any kind. The two things the probe genuinely
  influenced are Deviation 2 (the declared range, forced by measured intensities of 4–18) and the
  decision to run at all rather than report a guard failure — both of which are stated above.
- The reported run's own numbers differ from the probe's, as they must: different Philox key.

**7. The figure is gitignored under the house `*.png` rule** (`.gitignore:401`), like every other
reproducible spike figure and like 13-12's `p13_confusion.png`. It is **not** "available for reuse at
no further compute cost": if it is lost, regenerating it costs a full re-run of
`julia --project=spike spike/p13/run_alpha_series.jl` (~2.5 min CPU). Because the substrate is
counter-based Philox keyed on the global index, that re-run is **byte-identical** — so a loss is a
compute cost, not a re-seed. `spike/p13/alpha_report.jld2` (53,391 bytes) is **committed and tracked**
and carries every number in this document, so no table above depends on the gitignored figure.

**8. This run created no bulk cache directory and wrote nothing under `spike/data/cache/`.** The two
`spike/data/cache/p13/` directories both predate it (2026-07-28): `7c65a1a9…` is 13-11's 50 MB
training pool and `57014f9e…` is a 456 KB smoke shard. `spike/data/cache/p11/` (54 MB, read-only to
this phase, depended on by a Phase-12 plan) was neither read nor written.

**9. A parallel agent is executing Phase 12 on this branch.** Every commit here used explicit paths
only (`git commit -- <paths>`); nothing was rebased or amended; `spike/test/runtests.jl` and the
`test_p12_*` files that agent is editing were left untouched. No index race occurred.

## Assumption Drift (advisory)

**A. The plan framed the one-rung lag between the summary's sign flip and the evidence's as THE gap
to be characterized; what was measured is that the lag is fully accounted for by the pre-registered
tau dead zone.**
- **Planned:** 13-RESEARCH §J3 and this plan describe `alpha_star` as "the phase's most quotable
  result about the correlation-versus-localization gap", implicitly treating whatever separation
  appears as the gap itself.
- **Actual:** `alpha_star = 0.25` is **exactly** the first rung whose induced rho (`ghat(m-bar) =
  −0.2499`) clears `P13_TAU = 0.15`, while the rung below it (`−0.1443`) sits inside the dead zone by
  0.006. The separation is not a residual to be explained; it is the D-05 label rule and the
  measured summary resolution doing precisely what they were defined to do.
- **Why:** `tau` was MEASURED (13-10) as the smallest rho separation the fixed 8×8 summary can
  reliably distinguish. Any ladder scored by a net trained under that cut must cross where the
  induced rho crosses tau. The agreement is a consistency check that passed, not a coincidence.
- **Materiality:** it changes the sentence 13-14 should write. "The net needs more segregation than
  the summary shows" would be an over-claim; "the net calls exclusion at the first rung the
  pre-registered hypothesis boundary permits, and declines below it — correctly, by the phase's own
  label rule" is what was measured. It also relocates the real limitation from the NET to the
  SUMMARY, which is where Phase 11 already put it.

**B. 13-12 left open whether the epoch-4 best-validation checkpoint of an early-overfitting training
run is sufficient for 13-13; it is, comfortably.**
- **Planned:** 13-11 recorded the early overfitting as an honest finding and explicitly deferred the
  question to the downstream gates; 13-12 answered it for D-12/D-13 only and said so.
- **Actual:** the checkpoint produces a monotone ladder in the mean on both heads, a per-image
  positive rate rising 0% → 27% → 86% → 98% → 100%, and a per-rung spread that NARROWS from 1.88 to
  0.82 nats — on a spatially-constructed transform the net never saw in training.
- **Why:** the persisted artifact is the *best-validation* checkpoint, so the overfitting consumed
  budget rather than weights, and the D-11 masked loss gives the trunk all three classes.
- **Materiality:** the "overfits early" finding stays on the record and stays honest, but it must not
  be carried forward as a caveat on any number in this document. It remains open only for 13-16.

## Known Stubs

None. Every number in this document was computed from `spike/p13/alpha_report.jld2` or, for the
`ghat(m-bar)` column in §3, derived from the persisted `m-bar` through the frozen
`spike/simulator/ghat.jl` map — a read of frozen code, not a re-measurement and not a fit. No
placeholder, mock or hardcoded value flows into any table. The runner's one early-return branch (the
`INVARIANT VIOLATION` path) did not fire: all eight invariants held on all 64 images.

## Requirements satisfied

**D-15** (the alpha-graded series run on the simulated substrate and reported as supporting evidence,
never as a gate; the correlation-versus-localization gap addressed explicitly and quantified; the
single physical segregated anchor deferred to Phase 16 with a named reason and the seal left intact),
**D-16** (the three curves plus `alpha_star`; the four transform invariants re-verified on the
REPORTED substrate before any curve was read; the `:simulated` arm of the two-valued
`P13_ALPHA_SUBSTRATE` reported separately from and never averaged with 13-16's `:real` arm),
**D-08** (evidence emitted as two log Bayes factors against a structural-zero random reference,
through the documented `three_way_log_bf` read surface with the PERSISTED per-head correction),
**D-01** (the reserved stream consumed at `P13_ALPHA_COUNTER = 4`, asserted disjoint from both the
training-pool and reported-gate counters; atomic persist with a reopen integrity check; the
decoupling proof asserted at RUN TIME; CPU-only; no package installed),
**D-04** (every pre-registered value READ from `consts.jl` and written INTO the artifact; nothing
inlined; `consts.jl` byte-unchanged; `P13_ITERATION_ALLOWANCE` untouched and UNSPENT).

## Commits

| Commit | Task | Description |
|--------|------|-------------|
| `0943a02` | Task 1 | run the alpha ladder once — the evidence crosses zero exactly where tau says it should |
| _(this)_ | Task 2 | record the correlation-versus-localization reading and the Phase-16 deferral |

Task 2 produced **no new artifact and no code change**, as the plan requires: it read
`spike/p13/alpha_report.jld2` and wrote this SUMMARY. The artifact is byte-identical to what Task 1
committed.

## Self-Check: PASSED

- `spike/p13/run_alpha_series.jl` — FOUND (46,417 bytes, 771 lines)
- `spike/p13/alpha_report.jld2` — FOUND (53,391 bytes, 65 keys, **no `.tmp` sibling**)
- `spike/figures/p13_alpha_ladder.png` — FOUND (209,646 bytes; gitignored, see Disclosure 7)
- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-13-SUMMARY.md` — FOUND
- commit `0943a02` — FOUND in history
- `git diff --quiet HEAD -- spike/p13/consts.jl` — exit 0
- `git diff --quiet HEAD -- src/ spike/Project.toml spike/Manifest.toml corpus/` — exit 0
- `git status --porcelain corpus/data` — empty
