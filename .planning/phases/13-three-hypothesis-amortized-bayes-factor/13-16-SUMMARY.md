---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 16
status: complete
subsystem: amortized-inference
tags: [d-15-amended, d-16, real-images, reported-not-gated, qualitative-only, ood, named-limit, phase-16-deferral]
requires:
  - "spike/p13/three_way_net.jld2 (13-11) -- the trained two-head evidence net; its MEASURED per-head log-odds were consumed through the PERSISTED artifact and never re-measured"
  - "spike/p13/real_images.jl (13-15) -- real_pair / real_alpha_ladder / real_provenance / real_grid_provenance / verify_real_readonly_digest"
  - "spike/p13/alpha_series.jl (13-04) -- alpha_ladder / verify_alpha_invariants, the ONE shared transform code path"
  - "spike/p13/preconditions.jl (13-09) -- p13_require_phase11 / load_p13_basis / p13_encode_pair, and the Phase-11 SC2_RUNGS lambda ladder"
  - "spike/p13/consts.jl (13-01/13-10) -- the byte-locked pre-registration incl. MEASURED P13_TAU = 0.15"
  - "spike/data/cache/p13/7c65a1a9... (13-11) -- the 48,000-pair training pool, READ-ONLY, as the Phase-13 in-distribution OOD reference"
  - "artifacts/amended_v2/grid_8/ -- the frozen shipped bundle, READ-ONLY, as the OOD COMPARISON REFERENCE ONLY (never an evidence basis, D-02)"
provides:
  - "spike/p13/run_p13_realimage.jl -- the reported, non-gating real-image runner (no @test, no bar, seal self-check, digest before/after)"
  - "spike/p13/realimage_report.jld2 -- the lambda sweep with paired OOD verdicts, both detectors side by side, the real alpha ladders, alpha_star_real, both fixtures' provenance and the digest pair"
  - "spike/figures/p13_realimage.png -- the two-panel figure (script artifact, gitignored under the house *.png rule)"
  - "The named-limit entry and the report-header block, written to be lifted VERBATIM by plan 13-14"
affects:
  - "13-14 -- the phase report lifts the named limit, the report header and the four honesty paragraphs from here"
  - "Phase 16 -- the sealed holdout is left INTACT; the labelled real segregation check is deferred there"
tech-stack:
  added: []
  patterns:
    - "The evidence read and the OOD read are computed together and emitted on ONE row, so no printed or persisted quantity can be separated from its misspecification verdict"
    - "The shipped comparison bundle is SELECTED BY THE PRE-REGISTRATION: the chosen directory is the one whose own recorded id_threshold reproduces the frozen P13_REAL_OOD_SHIPPED_THRESHOLD"
    - "Persist BEFORE the read-only digest comparison can throw; the figure is wrapped so a backend failure cannot cost the artifact"
key-files:
  created:
    - spike/p13/run_p13_realimage.jl
    - spike/p13/realimage_report.jld2
    - spike/figures/p13_realimage.png
  modified: []
decisions:
  - "The frozen shipped OOD reference lives in artifacts/amended_v2/grid_8/, NOT artifacts/grid_8/ as the plan's key_links stated. Which bundle is the reference is DERIVED, not chosen: the runner selects the candidate whose recorded report.ood.id_threshold reproduces the frozen P13_REAL_OOD_SHIPPED_THRESHOLD, and it reproduced 433.68840 / 179.13678 exactly."
  - "The Phase-13 OOD channel is built by the IDENTICAL shipped recipe (continuous-rows-only Mahalanobis, ridge 1e-6, operating point at the copied OOD_ID_QUANTILE = 0.95) on the net's OWN training pool, read from the trained artifact's meta.pool_dir. Two detectors built by different recipes would compare recipes, not nets."
  - "No third channel pair was added. Both pre-registered pairs include the nuclear counterstain; that is a limit of the frozen pre-registration and is NAMED rather than repaired by choosing a new pair after the fixtures had been measured."
metrics:
  duration: ~1.0-1.3 min wall clock per reported run; ~70 min total including authoring and verification
  completed: 2026-07-29
  tasks: 2
  commits: 2
---

# Phase 13 Plan 16: Real-Image Arm (D-15 AMENDED) Summary

**On both unmodified real pairs, at every pre-registered lambda rung, in both read directions, the
three-way net's descriptive verdict is RANDOM -- and every one of those reads is OOD-flagged.** Both
log Bayes factors are negative everywhere: `log BF(C:R)` runs -0.73 to -0.47 (positive as sample) and
-4.12 to -3.94 (negative as sample), `log BF(E:R)` sits near -10.3 and -9.0. The Phase-13 OOD channel
scores these fixtures at density **417.30 against its own in-distribution threshold 167.54 (2.49x
over)**; the shipped detector, read as a frozen comparison reference, scores **433.69 against 179.14
(2.42x over)**, reproducing the pre-registered constants exactly. **A12 outcome: agreement** -- both
detectors flag both fixtures, in both directions.

**Nothing above is a correctness result, and it cannot become one.** `P13_REAL_IS_GATED = false` and
`P13_REAL_QUALITATIVE_ONLY = true` were frozen before any Phase-13 number existed; the runner
contains zero `@test` lines (asserted by source grep) and no threshold is defined for any real-image
quantity. `spike/p13/consts.jl` is byte-unchanged (git blob `5a4ea222...`, identical to the sha
13-11, 13-12 and 13-13 all recorded) and **`P13_ITERATION_ALLOWANCE` remains 1 of 1, UNSPENT**. The
sealed holdout was not opened.

---

## 1. Section 1 -- the unmodified read, every rung, every verdict beside its OOD flag

Primary arm, `P13_REAL_CHANNEL_PAIR = (1, 2)`. The sweep is `P13_REAL_LAMBDA_READS_RULE =
:full_phase11_ladder`, realized as Phase 11's `SC2_RUNGS = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)`; the
headline is `P13_REAL_LAMBDA_HEADLINE_RULE = :widest_rung`, realized at **3.0**, asserted equal to
Phase 11's `LAMBDA_MAX` and to the frozen `P13_REAL_LAMBDA_HEADLINE_EXPECTED`.

| direction | lambda | log BF(C:R) | log BF(E:R) | OOD density | OOD threshold | is_ood |
|---|---|---|---|---|---|---|
| positive as sample | 0.25 | -0.725239 | -10.261014 | 417.2974 | 167.5446 | **true** |
| positive as sample | 0.5 | -0.701851 | -10.262200 | 417.2974 | 167.5446 | **true** |
| positive as sample | 1.0 | -0.654947 | -10.264497 | 417.2974 | 167.5446 | **true** |
| positive as sample | 1.5 | -0.607885 | -10.266683 | 417.2974 | 167.5446 | **true** |
| positive as sample | 2.0 | -0.560664 | -10.268779 | 417.2974 | 167.5446 | **true** |
| positive as sample | 2.5 | -0.513290 | -10.270785 | 417.2974 | 167.5446 | **true** |
| **positive as sample** | **3.0 (HEADLINE)** | **-0.465771** | **-10.272712** | **417.2974** | **167.5446** | **true** |
| negative as sample | 0.25 | -4.118168 | -8.994916 | 417.2974 | 167.5446 | **true** |
| negative as sample | 0.5 | -4.102732 | -9.000272 | 417.2974 | 167.5446 | **true** |
| negative as sample | 1.0 | -4.071577 | -9.010961 | 417.2974 | 167.5446 | **true** |
| negative as sample | 1.5 | -4.040039 | -9.021613 | 417.2974 | 167.5446 | **true** |
| negative as sample | 2.0 | -4.008124 | -9.032229 | 417.2974 | 167.5446 | **true** |
| negative as sample | 2.5 | -3.975831 | -9.042814 | 417.2974 | 167.5446 | **true** |
| **negative as sample** | **3.0 (HEADLINE)** | **-3.943164** | **-9.053374** | **417.2974** | **167.5446** | **true** |

**The headline verdict in one sentence.** At the widest, most conservative rung `lambda = 3.0`, both
log Bayes factors are negative in both read directions, so the descriptive `argmax(logBF_C, 0,
logBF_E)` rule lands on the RANDOM reference for both specimens -- read out of a net whose
misspecification channel is simultaneously reporting that these inputs sit 2.49x outside its own
training distribution.

**The density is constant across the sweep, and that is a property of the channel, not an error.**
The OOD density reads the standardized summaries; the lambda conditioning row is appended *after* the
pair encoding (`P13_LAMBDA_PLACEMENT = :append_after_pair_encode`), so it never enters the summary the
density is computed on. It is printed on every row regardless, because the standing D-15 requirement
is that no real-image number is ever quotable without its verdict beside it. The per-acquisition
scores behind the pair max are **positive 417.2974** and **negative 321.1505**.

**The lambda response is essentially flat, and that is a PRE-AUTHORIZED expected outcome.** Across
the whole ladder the exclusion evidence moves by **0.0117 nats** (positive as sample) and **0.0585
nats** (negative as sample); the coloc evidence moves by **0.2595** and **0.1750** nats, and in the
direction a conservative reading predicts -- widening the stated registration uncertainty makes the
coloc evidence *less* negative, i.e. weaker against the random reference. Phase 11 already
established that the `rho_true` posterior width is genuinely flat in lambda (RMSE ratio 1.0003), so a
flat Phase-13 lambda response was recorded as an expected finding before this run
(`preconditions.jl` ruling 4, and 13-12's must_have "INCLUDING A FLAT RESPONSE AS A FINDING"). It is
not evidence that the conditioning row is dead: Phase 11 measured shift marginals whose SDs track
lambda with ratios 4.84-9.05.

### 1b. Why RANDOM is the coherent answer here, by the phase's own label rule

Mapping the measured `m-bar` through the FROZEN `spike/simulator/ghat.jl` map (a read of frozen code,
not a fit): the positive fixture sits at `rho = +0.33042`, the negative at `rho = +0.21480`. The D-05
two-factor cut is `COLOC iff rho_s > tau AND rho_s - rho_c > tau` at the MEASURED `P13_TAU = 0.15`:

| read | rho_s | rho_c | contrast | D-05 label | net's descriptive verdict |
|---|---|---|---|---|---|
| positive as sample | +0.33042 | +0.21480 | **+0.11562** (inside the tau dead zone) | RANDOM | RANDOM |
| negative as sample | +0.21480 | +0.33042 | **-0.11562** (inside the tau dead zone) | RANDOM | RANDOM |

Both specimens are individually correlated (`rho_s > tau` on both), but they are correlated to
*similar degrees*, and "less colocalized than the control" is not segregation -- which is precisely
the counter-example the two-factor cut exists to prevent. So the pre-registered label rule assigns
RANDOM to both reads, and the net agrees on both. **This is a consistency check that passed, not a
correctness check, and it must not be quoted as one:** there is no ground-truth label here, so what
was checked is that the net agrees with a rule derived from the same summary it reads.

### 1c. Section 5 -- the REDUNDANCY arm (`P13_REAL_REDUNDANCY_PAIR = (1, 3)`)

Labelled REDUNDANCY, not a separate claim: it re-reads the SAME two specimens through the second
pre-registered channel pair, so it adds no independent specimen.

| direction | lambda (min -> max) | log BF(C:R) | log BF(E:R) | OOD density | OOD threshold | is_ood |
|---|---|---|---|---|---|---|
| positive as sample | 0.25 -> **3.0** | -2.405902 -> **-2.178877** | -9.887478 -> **-9.886481** | 607.5728 | 167.5446 | **true** |
| negative as sample | 0.25 -> **3.0** | -1.676072 -> **-1.450446** | -9.913825 -> **-10.010441** | 607.5728 | 167.5446 | **true** |

Same qualitative outcome -- both log Bayes factors negative at every rung, verdict RANDOM, both reads
OOD-flagged -- through a different channel pair, at a **higher** OOD density (607.57, 3.63x over
threshold; per-acquisition 410.24 positive / 607.57 negative). It corroborates the primary arm and
adds nothing independent.

---

## 2. Section 3 -- the real alpha ladder through the net, invariants first

Both conditions, at the headline `lambda = 3.0`, on the frozen `P13_REAL_ALPHA_GRID` (= the simulated
arm's `P13_ALPHA_LADDER`, so the rungs are comparable rung for rung). `verify_alpha_invariants` ran on
the REPORTED substrate BEFORE any curve was read; all eight invariants hold on both conditions, so
the runner's `INVARIANT VIOLATION` early-return branch did not fire.

| condition | all 8 invariants | source zeros | mask fraction | max rel. intensity err | realized max |
|---|---|---|---|---|---|
| positive | hold | 2 | 0.134933 | 2.86e-16 | 0.773286 |
| negative | hold | 3 | 0.229189 | 1.72e-16 | 0.586862 |

### positive (control held fixed at the unmodified `negative` fixture)

| alpha | log BF(E:R) | log BF(C:R) | m-bar | ghat(m-bar) | n_missing | max_value | is_ood |
|---|---|---|---|---|---|---|---|
| 0.000 | -10.2727 | -0.465771 | +0.32916 | +0.33042 | 0 | 0.65669 | **true** |
| 0.125 | -9.56408 | -1.14050 | +0.27312 | +0.24648 | 0 | 0.57538 | **true** |
| 0.250 | -8.65155 | -1.88244 | +0.21231 | +0.16963 | 0 | 0.57779 | **true** |
| 0.375 | -7.49962 | -2.61873 | +0.14958 | +0.07102 | 0 | 0.61037 | **true** |
| 0.500 | -6.02779 | -3.32545 | +0.08690 | +0.00046 | 0 | 0.64296 | **true** |
| 0.625 | -4.16637 | -4.06013 | +0.02465 | -0.08711 | 0 | 0.67554 | **true** |
| 0.750 | -1.90658 | -4.93801 | -0.03851 | -0.17080 | 0 | 0.70812 | **true** |
| **0.875** | **+0.70625** | -5.94974 | -0.10423 | -0.24982 | 0 | 0.74070 | **true** |
| 1.000 | +3.40055 | -6.83244 | -0.16787 | -0.34701 | 0 | 0.77329 | **true** |

**`alpha_star_real("positive") = 0.875`** -- a RUNG of the frozen ladder, never interpolated.

### negative (control held fixed at the unmodified `positive` fixture)

| alpha | log BF(E:R) | log BF(C:R) | m-bar | ghat(m-bar) | n_missing | max_value | is_ood |
|---|---|---|---|---|---|---|---|
| 0.000 | -9.05337 | -3.94316 | +0.24805 | +0.21480 | 0 | 0.58686 | **true** |
| 0.125 | -7.55256 | -5.17191 | +0.18357 | +0.12265 | 0 | 0.51424 | **true** |
| 0.250 | -5.55352 | -6.38048 | +0.11407 | +0.03105 | 0 | 0.44162 | **true** |
| 0.375 | -3.19712 | -7.33669 | +0.04554 | -0.05786 | 0 | 0.41901 | **true** |
| 0.500 | -0.79561 | -8.01521 | -0.01866 | -0.14497 | 0 | 0.44459 | **true** |
| **0.625** | **+1.54869** | -8.53514 | -0.07918 | -0.21929 | 0 | 0.47018 | **true** |
| 0.750 | +3.92881 | -9.02240 | -0.13912 | -0.30803 | 0 | 0.49576 | **true** |
| 0.875 | +6.38335 | -9.57032 | -0.20059 | -0.38272 | 0 | 0.52134 | **true** |
| 1.000 | +8.24502 | -10.10820 | -0.25198 | -0.44184 | 0 | 0.54693 | **true** |

**`alpha_star_real("negative") = 0.625`** -- a RUNG of the frozen ladder, never interpolated.

### Does the net's verdict flip where the summary's sign flips? No -- it flips where TAU says it should

The summary's own `m-bar` sign flip and the evidence's zero crossing are at different rungs, and the
gap is accounted for by the measured `P13_TAU = 0.15` dead zone rather than left unexplained:

| condition | m-bar zero crossing (linear, between rungs) | first rung whose `ghat` clears -tau | `alpha_star_real` (net) | lag |
|---|---|---|---|---|
| positive | ~0.674 | 0.750 (`ghat` = -0.17080) | **0.875** | one rung beyond the boundary |
| negative | ~0.464 | 0.625 (`ghat` = -0.21929) | **0.625** | none -- exactly at the boundary |

On the negative fixture the net calls exclusion at *precisely* the first rung where the induced rho
clears the pre-registered hypothesis boundary; the rung below it (`ghat` = -0.14497) sits inside the
dead zone by 0.005. On the positive fixture the net is one rung more conservative than the boundary.
**With n = 2 specimens that is a two-point observation, not an estimate of a bias**, and it is
reported as such. It matches the pattern 13-13 measured on simulated substrate, where the crossing
also landed at the first rung clearing tau.

### The two arms are NEVER averaged

`P13_ALPHA_SUBSTRATE = (:simulated, :real)` is two-valued and pre-registered, and what differs is
what `alpha = 0` MEANS. On the simulated substrate (plan 13-13) both acquisitions are drawn at
`rho_true = 0.0`, so the ladder runs **random -> exclusion by construction** and crosses at
`alpha_star = 0.25`. Here `alpha = 0` is a **moderately colocalized** pair (`m-bar` +0.3292 / +0.2481)
whose control is the *other* colocalized specimen, so the ladder runs **moderately colocalized ->
near-exclusion** and crosses much later (0.875 / 0.625). Those are different experiments with
different starting points; a merged curve would be a statement about neither, and the difference
between 0.25 and 0.875 is a difference of starting point, not of net quality.

---

## 3. The OOD paragraph -- the arm's headline honesty item

**Both fixtures are flagged out-of-distribution by both detectors, in both read directions, on both
pre-registered channel pairs.**

| detector | density | threshold | ratio to threshold | is_ood |
|---|---|---|---|---|
| **Phase-13** (this net, Phase-11 registration-aware basis) | 417.2974 | 167.5446 | **2.491x** | true |
| **shipped** (frozen comparison reference, `artifacts/amended_v2/grid_8/`) | 433.6884 | 179.1368 | **2.421x** | true |
| frozen pre-registration constants | `P13_REAL_OOD_SHIPPED_DENSITY` = 433.69 | `P13_REAL_OOD_SHIPPED_THRESHOLD` = 179.14 | 2.42x | true |

The shipped read **reproduced the frozen constants exactly** (433.68840 against a recorded operating
point of 179.13678), which is what makes the comparison trustworthy rather than merely adjacent.
**The observed A12 outcome is `agreement`.**

**This is not a bug and it is not a reason to skip the check.** It is the net's own misspecification
channel reporting, correctly, that real microscopy at this frame size and these statistics sits
outside the simulator's training distribution. Suppressing it, softening it or "fixing" it would be
the single most damaging thing this arm could do, which is why every number in every table above
carries its density, its threshold and its boolean verdict on the same row, in the printed output and
in the persisted artifact alike.

**What the agreement buys, stated narrowly.** The registration-aware Phase-11 basis did NOT move real
microscopy back inside the training distribution: 2.49x over its own threshold against the shipped
net's 2.42x over its own. The two numbers are not directly comparable in magnitude -- different
standardizer, different training pool, different threshold -- but the *verdict* is the same, and that
agreement corroborates the named limit rather than weakening it. Had the two diverged, that would
itself have been the finding.

**How the Phase-13 detector was built, so nothing here is a new invention.** It is the shipped recipe
applied to the Phase-13 net's own inputs: continuous-rows-only Mahalanobis with the same 1e-6 ridge,
fit on the 96,000 standardized acquisitions of the net's own 48,000-pair training pool (read from the
trained artifact's `meta.pool_dir`, opened read-only, never regenerated), with the operating point at
the copied shipped `OOD_ID_QUANTILE = 0.95`. Realized ID score quantiles: q50 = 40.17, q90 = 149.57,
q95 = **167.54**, q99 = 199.47, max = 290.69 -- so the real fixtures score above the *maximum* of
96,000 simulated acquisitions. **Disclosed rather than buried:** that pool is class-frequency
stratified (D-07), so its rho marginal is not the prior's; the two detectors are each referenced to
their own net's training distribution, which is the comparison this arm is actually making.

---

## 4. The naming-correction paragraph (Open Question 8, answered YES)

**`positive`/`negative` are the original package's BIOLOGICAL test conditions, not colocalization
labels.** The "negative" pair measures mean patch correlation **+0.2481** (`rho_true` ~ +0.215) -- a
positively correlated pair, not an anti-correlated one. Nothing in this phase treats `negative/` as an
exclusion example, and **a reader who assumes `negative` means "not colocalized" will misread every
real-image figure in the phase**: they would read the RANDOM verdict on the "negative" specimen as a
miss, when the specimen is in fact moderately colocalized and RANDOM is what the pre-registered label
rule assigns it. `P13_REAL_NAMING_CORRECTION` is printed verbatim in the runner's banner before any
code runs and is persisted verbatim into the artifact.

**A second naming trap, named here for the first time in an executed plan.** Both pre-registered
channel pairs -- `P13_REAL_CHANNEL_PAIR = (1, 2)` and `P13_REAL_REDUNDANCY_PAIR = (1, 3)` -- include
channel 1, which `test/runtests.jl:105` records as the DAPI/Hoechst **nuclear counterstain**, not a
target protein; `spike/simulator/ghat.jl` records the c1/c2 figures as SUPERSEDED for colocalization
purposes by the c2/c3 green/red pair. This plan did NOT add a c2/c3 arm. Choosing a new channel pair
after the fixtures had been measured is exactly the move the D-04 anti-snooping contract forbids, so
the limit is named rather than repaired, and it belongs in the manuscript alongside the folder-name
correction.

---

## 5. The qualitative-only paragraph

**The six committed TIFFs carry no colocalization ground-truth label.** There is no recorded rho, no
recorded Delta-rho, no segregation degree and no provenance metadata asserting one. So this arm can
show that the three-way verdicts behave sensibly on real microscopy -- that ingestion works, that the
evidence responds monotonically to a constructed ladder, that the verdict tracks the pre-registered
hypothesis boundary, that nothing degenerates -- and **it cannot show the verdicts are correct**,
because correctness requires a label the data does not have. **No pass/fail threshold is defined for
any real-image quantity**, so this arm cannot accidentally become a gate and no Phase-13 verdict turns
on it.

**State the n plainly.** This is **two specimens** -- one `positive` and one `negative` biological
condition, six TIFF files in total because each specimen carries three channels. The two
pre-registered channel pairs give four reads, but they are four reads of the *same two specimens*, so
the effective independent n is **2**, not 4 and not 6. Two specimens cannot support a rate, a
coverage claim, a confidence interval or a comparison of conditions. What they can support is exactly
what is claimed above: that the read chain works end to end on real files, that the numbers it
produces are coherent with the phase's own label rule, and that the misspecification channel fires.
**Labelled real segregation validation remains deferred to Phase 16's blind evaluation, which is what
the sealed holdout exists for and precisely why Phase 13 did not open it.**

---

## 6. The Phase-02 lineage paragraph

`spike/simulator/ghat.jl:40` and `spike/NOTES.md:238-247` already record, as a frozen Phase-02
finding, that *"negative tail : physically reachable in real fluorescence = false => negative mu-prior
tail documented PRIOR-ONLY; consistency scoped to realized range"* -- real anti-correlation was never
observed, and even the "negative" biological control sits at `+0.25`. **The exclusion hypothesis,
which Phase 13 promotes to a first-class hypothesis, therefore has no observed real-data instance
anywhere in this project.** Every unmodified real image measured here confirms it: `m-bar` +0.32916
and +0.24805, both positive.

The D-16 construction sharpens that verdict rather than refuting it. At `alpha = 1` it produces
`m-bar` = -0.16787 / -0.25198 (`rho` ~ -0.347 / -0.442) **from real microscopy pixels** -- the first
negative induced mu this project has obtained from real data. The sentence intended for the
manuscript, verbatim:

> **Negative induced mu is CONSTRUCTIBLE from real microscopy pixels via the D-16 mask-based
> reassignment, but it has NOT been observed to occur naturally in the images this project holds.**

---

## 7. The scope-limit note

The real ladder spans `ghat(m-bar)` from **+0.33042 down to -0.44184**, i.e. `rho_true` in roughly
`[-0.44, +0.33]` -- a densely covered interior region of `GHAT_RHO_KNOTS`, comfortably away from the
+/-0.99 clamp atoms where the prior-atom asymmetry (~4.9% of prior mass at the negative clamp against
~1.9% at the positive one, about 2.5:1 against the exclusion end) lands hardest. **It probes the sign
transition and the near-exclusion regime; it does NOT probe the deep-exclusion tail.** A reader who
assumes otherwise will over-read the result. (That tail was tested, and cleared, by the D-12 gate on
simulator ground truth: exclusion AUC 0.997257 at `|rho| > 0.9`.)

**Honest caveat on the ladder endpoint**, carried forward from the transform's own design: at
`alpha = 1` the mask region is not EMPTY, it is AT the background floor (0.006241 / 0.005905). The
endpoint is "ch2 reduced to background wherever ch1 has objects", not "ch2 absent" -- do not write
"fully disjoint" without that qualifier.

---

## 8. NAMED LIMIT -- liftable verbatim into 13-14 and `docs/amortized.md`

> **Named limit (Phase 13): the exclusion hypothesis has no labelled real-data validation.** The
> three-way evidence network is gated on simulator ground truth. The real-image arm runs on six
> committed microscopy TIFFs (`test/test_images/`) that carry **no colocalization ground-truth
> label**; it is a qualitative behaviour check, not a correctness check, and it is explicitly not a
> gate -- no pass/fail threshold is defined for any quantity in it. It covers **two specimens**, so it
> can support no rate and no coverage claim. Both specimens are additionally flagged
> **out-of-distribution** by both amortized detectors, in both read directions and on both
> pre-registered channel pairs: the Phase-13 net scores density **417.30 against its own
> in-distribution threshold 167.54 (2.49x)** and the shipped 8x8 bundle **433.69 against 179.14
> (2.42x)** -- agreement, not divergence, so the registration-aware basis did not move real microscopy
> back inside the training distribution. Both pre-registered channel pairs include the nuclear
> counterstain, which is a further limit of the frozen pre-registration and was named rather than
> repaired. The provenance manifest holds exactly **one** physically-segregated anchor; it is
> `split = sealed_holdout`, `sha256 = "PENDING-FETCH"`, `bytes = 0`, and reserved for the **Phase-16**
> blind evaluation, so it was deliberately **not** consumed here -- consuming it would irreversibly
> burn Phase 16 on the very hypothesis Phase 16 exists to evaluate. **Labelled real segregation
> validation is deferred to Phase 16.**

---

## 9. REPORT-HEADER TEXT -- liftable verbatim into 13-14

> **Scope of evidence.** The Phase-13 gate is simulator ground truth (`rho_true < -tau`, exact labels,
> well-powered, `P13_GATE_M = 4000`). Two further arms are reported and **neither is a gate**: the
> D-16 semi-synthetic alpha-graded random-to-exclusion series (`P13_ALPHA_GATED = false`), and a
> qualitative check on six committed real microscopy TIFFs (`P13_REAL_IS_GATED = false`,
> `P13_REAL_QUALITATIVE_ONLY = true`). The TIFFs carry no colocalization label -- the folder names
> `positive`/`negative` are the original package's biological test conditions, and the "negative" pair
> in fact has mean patch correlation **+0.25**. The real-data arm can show the three-way verdicts
> behave sensibly on real microscopy; **it cannot show they are correct.** It covers two specimens,
> and both of them are flagged out-of-distribution by both the Phase-13 and the shipped detector, so
> every real-image number in this report is printed beside its OOD verdict and its lambda. A real
> image has no known registration uncertainty, so the whole pre-registered lambda ladder is reported
> and the widest, most conservative rung (`lambda = 3.0`) is the headline; lambda was never estimated
> from the images. The sealed-holdout rows of the provenance manifest were **not** opened.

---

## 10. Open Questions 6, 7 and 8 -- closed

**OQ 6 (at which lambda is the real check read, and which rung is the headline?)** -- CLOSED. The
pre-registered rules resolved without ambiguity: `:full_phase11_ladder` realized as Phase-11's
`SC2_RUNGS = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)`, and `:widest_rung` realized at **3.0**, asserted
at run time to equal both Phase-11's `LAMBDA_MAX` and the frozen
`P13_REAL_LAMBDA_HEADLINE_EXPECTED = 3.0`. Lambda was never estimated from the images. The response
across the ladder is essentially flat (0.0117-0.0585 nats on the exclusion head, 0.1750-0.2595 on the
coloc head), which Phase 11's own flat-in-lambda result had already made an expected outcome.

**OQ 7 (should the Phase-13 net's OOD verdict be compared to the shipped net's?)** -- CLOSED, and the
answer was YES. Both were measured and reported side by side: Phase-13 **417.30 / 167.54 (2.49x)**
against shipped **433.69 / 179.14 (2.42x)**, with the shipped read reproducing the frozen constants
exactly. **A12 outcome: agreement.** The shipped bundle was opened read-only as a comparison reference
only and was never an evidence basis; every log Bayes factor in this document came from the Phase-13
net on Phase 11's registration-aware basis, with no fallback path.

**OQ 8 (does the `positive`/`negative` naming need an explicit correction?)** -- CLOSED, and the
answer was YES. `P13_REAL_NAMING_CORRECTION` is printed verbatim in the runner's ALL-CAPS banner
before any code runs, persisted verbatim into `realimage_report.jld2` as `naming_correction`, and
restated in section 4 above with its measured number (+0.2481) for 13-14 and the manuscript to lift.

---

## 11. The read-only proof

`sha256(git ls-files -s test/test_images)`, recorded by the runner before and after the run:

| | digest |
|---|---|
| before | `eaee22f9185460fd910788026e7a33511f6de373f313459b96bb43f3fc1ef276` |
| after | `eaee22f9185460fd910788026e7a33511f6de373f313459b96bb43f3fc1ef276` |

**They match**, and they match the digest 13-15 recorded, so the committed fixtures are byte-identical
across both plans. `git status --porcelain test/test_images/` is empty and
`git diff --quiet HEAD -- test/ src/ spike/Project.toml spike/Manifest.toml` exits 0, both observed
after the commit. The comparison is deliberately made AFTER the artifact is persisted: the digest pair
belongs in the record whatever it says, and an artifact lost to an assertion would be the worse
failure.

---

## Verification performed (every command was RUN; the observed result is what is reported)

| Check | Observed |
|-------|----------|
| `julia --project=spike spike/p13/run_p13_realimage.jl` | **exit 0**, 1.01 min wall clock |
| Task-1 `<verify>` (`is_gated == false`, four keys present, digests equal) | prints **`realimage artifact ok`** |
| Task-1 pairing criterion (`:ood_density`, `:ood_is_ood`, `:lambda` on every sweep row) | prints **`paired ok (14 rows)`** |
| Task-2 `<verify>` (headline lambda, both OOD records, per-condition `alpha_star_real`) | prints headline `3.0`, both records, `(positive = 0.875, negative = 0.625)` |
| `spike/p13/realimage_report.jld2` exists, no `.tmp` sibling | 81,649 bytes, 57 keys, confirmed |
| `spike/figures/p13_realimage.png` rendered | 294,821 bytes |
| `grep -c '@test' spike/p13/run_p13_realimage.jl` | **0** |
| `grep -v '^\s*#' ... \| grep -c 'open_sealed_holdout\|fetch_verified\|bootstrap_anchor_hashes'` | **0** |
| the sealed provenance directory name in executable text | **0 occurrences anywhere in the file**, comments included |
| `grep -c '@__FILE__' ...` | **2** (the source self-check is present and live) |
| `grep -c 'Phase 16' ...` | **5** |
| `grep -v '^\s*#' ... \| grep -c '\b321\b\|1028\|1376'` | **0** (widths and dimensions derived or read from consts) |
| `git status --porcelain test/test_images/` | empty |
| `git diff --quiet HEAD -- test/ src/ spike/Project.toml spike/Manifest.toml` | **exit 0** (src/ and the manifests also asserted at RUN TIME by step 0) |
| `git diff --quiet HEAD -- spike/p13/consts.jl` | **exit 0 -- byte-unchanged**; artifact `consts_git_blob_sha` = `5a4ea2223ce11bcae2f974c643ec1e523458e8ac`, identical to 13-11/13-12/13-13 |
| `net_consts_sha == p13_consts_sha()` | `100e97a3...` both sides (asserted at run time) |
| derived input width vs. the trained net's | `three_way_input_dim(8, 1)` matched `h.input_dim`, asserted, never a literal |
| Post-commit deletion scan | no deletions |
| Training pool / Phase-11 cache | `spike/data/cache/p13/7c65a1a9...` opened READ-ONLY (12 shards, 96,000 acquisitions); nothing written, nothing regenerated; `spike/data/cache/p11/` neither read nor written |

**Test files run.** Adding a file under `spike/p13/` is not inert: three Phase-13 test files discover
that directory with `readdir` and grep every `.jl` in it. All three were run directly:

| File | Result |
|---|---|
| `spike/test/test_p13_real.jl` | **exit 0** -- and its seal testset is now **34 pass / 0 broken**: the `@test_skip` that 13-15 left for this plan's runner is LIVE and passing |
| `spike/test/test_p13_alpha.jl` | **exit 0** -- all testsets, including the sealed-holdout grep over the new file |
| `spike/test/test_p13_calibration.jl` | **exit 0** -- including the D-14 absence grep, which also scans the new file |

The FULL spike suite still exits 1 at `spike/test/runtests.jl:169` on the Phase-4 `SPEEDUP_GATE`
(measured 68.35 against the pre-registered bar 100.0), which masks the entire Phase-13 include block.
That is **pre-existing**, already logged as `deferred-items.md` D-13-A by 13-11, and is not this
plan's to fix -- it is a pre-registered threshold in a Phase-4 file this plan does not own. Per-file
runs are the available signal for Phase 13.

---

## Deviations from Plan

### Auto-fixed issues

**1. [Rule 1 - Bug] The plan named the wrong shipped bundle for the OOD comparison reference**
- **Found during:** Task 1, on a read-only mechanics probe before the runner was written.
- **Issue:** the plan's `key_links` and step 7 name `artifacts/grid_8/ood_nulls_8.jld2` as the frozen
  comparison reference. That bundle does NOT carry the frozen numbers: its recorded
  `report.ood.id_threshold` is **99.2207**, and the real fixtures score **85.04** against it. The
  frozen `P13_REAL_OOD_SHIPPED_DENSITY = 433.69` / `P13_REAL_OOD_SHIPPED_THRESHOLD = 179.14` belong to
  `artifacts/amended_v2/grid_8/` -- the AMENDED v2 ship-gate bundle, which reproduces
  **433.68840347** against **179.13677813**. Reading the wrong bundle would have printed a comparison
  against a net the pre-registration is not quoting, and the discrepancy would have looked like a
  Phase-13 finding.
- **Fix:** the runner does not hard-code either path. It searches a documented candidate list and
  selects the bundle whose OWN recorded operating point reproduces the frozen constant -- so **the
  pre-registration is the selector**, not the executor. `reproduced_frozen = true` is persisted, and
  when no bundle on disk reproduces the constant the run says so and quotes the frozen constants
  alone rather than silently comparing against something else.
- **Files:** `spike/p13/run_p13_realimage.jl`. **Commit:** `2112bed`.
- **Note for 13-14:** `artifacts/` is UNTRACKED in this repository, so the comparison reference is a
  local artifact. The frozen constants in `consts.jl` are the durable record; the recomputation is
  corroboration, not the source.

**2. [Rule 3 - Blocking] There is no Phase-13 OOD null bundle, so one had to be constructed**
- **Found during:** Task 1.
- **Issue:** the plan requires the Phase-13 net's OOD density AND threshold beside the shipped
  reference, but Phase 13 ships no OOD null artifact and `src/amortized/ood.jl` is not loadable
  standalone (its signatures mention `MultiChannelImage`), while the spike twin
  `spike/validation/ood.jl` drags the whole Phase-5 harness (a frozen Phase-4 model load) into a file
  that needs one covariance.
- **Fix:** `p13_real_fit_density` / `p13_real_maha` / `p13_real_cont_rows` are ATTRIBUTED COPIES of
  the shipped `fit_ood_nulls` / `maha_score` / `_summary_row_partition`, the same technique
  `spike/p13/net.jl` already uses for `pair_encode`. They are fit on the net's OWN training pool,
  located from the trained artifact's `meta.pool_dir` (never guessed, never routed through the cache
  layer, so no pool can be created, invalidated or regenerated), at the COPIED shipped
  `OOD_ID_QUANTILE = 0.95`. **No number was chosen:** the recipe, the ridge and the quantile are all
  copied from shipped, pre-registered code, and the resulting flag scores nothing -- it is a reported
  diagnostic that travels beside every quantity.
- **Files:** `spike/p13/run_p13_realimage.jl`. **Commit:** `2112bed`.

**3. [Rule 1 - Bug] The figure's per-point alpha vector failed the pinned Makie at draw time**
- **Found during:** Task 1, on the second run (the first run's figure used a scalar alpha and lost the
  coloc-head markers entirely).
- **Issue:** marking OOD-flagged points hollow by passing a per-point `alpha` vector to `scatter!`
  throws `Failed to resolve cairo_attributes` in the pinned CairoMakie/Makie; the runner's `try/catch`
  around the figure did its job -- the artifact was already persisted and the run still exited 0 --
  but the PNG was not produced.
- **Fix:** the flagged and unflagged points are drawn as two INDEX-SPLIT `scatter!` calls (hollow with
  a stroke, filled) instead of one call with a per-point alpha. All 28 markers on the sweep panel are
  hollow, which is itself the finding.
- **Files:** `spike/p13/run_p13_realimage.jl`. **Commit:** `2112bed` (the runner was committed after
  the corrected render).

**4. [Rule 3 - Blocking] The figure is composed in the runner, not routed through a `figures.jl`
helper**
- **Found during:** Task 1.
- **Issue:** the plan asks for the figure "through the existing `spike/validation/figures.jl`
  helpers". Same two obstacles plans 13-12 and 13-13 already recorded: that file has no sweep or
  ladder helper, and its `_val_fig_path` hard-joins every output into `spike/validation/figures/`,
  which is not where this deliverable lives.
- **Fix:** `p13_real_figure` composes the two panels on the CairoMakie surface, the headless backend
  activation and the caption discipline that `figures.jl` establishes and that the runner inherits by
  including it. **No plotting dependency was added.**
- **Files:** `spike/p13/run_p13_realimage.jl`. **Commit:** `2112bed`.

### Design choices the plan left to the executor (recorded, not deviations)

**5. The alpha ladder's CONTROL half is the OTHER unmodified real fixture, held fixed across the
ladder.** The plan says to score each rung but does not say what the control half is. The evidence net
reads a PAIR and the D-05 cut is a two-factor rule on the sample LEVEL and the sample-control
CONTRAST, so grading both halves would move both arguments together and the ladder would measure
nothing. Grading only the sample is the same choice 13-13 made on simulated substrate, which is what
keeps the two arms rung-comparable. The consequence, stated because it is load-bearing for section 2:
the real control is itself moderately colocalized (`rho_c` = +0.2148 / +0.3304), not random, so the
contrast factor is *already* satisfied well before the level factor is -- the crossing point is
driven by the sample level clearing `-tau`.

**6. No third channel pair was added.** See section 4. Both pre-registered pairs include the nuclear
counterstain; adding a c2/c3 arm after the fixtures had been measured would be an un-pre-registered
scope widening, so the limit is named instead.

### Disclosures (not deviations, but the ordering must be checkable)

**7. A READ-ONLY mechanics probe preceded the runner, and it is what caught deviation 1.** Before
writing the runner, the shipped-bundle read path was probed against the committed fixtures to answer
one question the plan could not answer from source: can the frozen 433.69 / 179.14 reference be
reproduced without reconstructing an estimator object? Four things bound what that means: the probe
consumed **no RNG stream** (this whole arm consumes none); it read only **frozen shipped artifacts and
committed fixtures**; the numbers it produced (433.688 / 179.137) were **already published in
13-RESEARCH §J4.6 and frozen into `consts.jl`**, so nothing was pre-observed that was not already
committed; and the one thing it influenced is deviation 1, which is stated above rather than folded
in silently. No Phase-13 evidence number was seen before the runner was written.

**8. The reported run was executed three times, and the numbers are byte-identical across all three.**
Runs 2 and 3 exist only because of the figure bug (deviation 3); the arm is deterministic -- six
committed files, a deterministic transform, a deterministic forward pass -- so every log Bayes factor,
every OOD density and both digests are identical in all three runs. The only field that differs is the
UTC `generated` stamp. Nothing was re-run in response to a result.

**9. The figure is gitignored under the house `*.png` rule** (`.gitignore:401`), like every other
reproducible spike figure and like 13-12's and 13-13's. It is **NOT** "available for reuse at no
further compute cost": if it is lost, regenerating it costs a full re-run of
`julia --project=spike spike/p13/run_p13_realimage.jl` (~1.0-1.3 min CPU). Because the substrate is
six committed files and the whole path is deterministic, that re-run is **byte-identical** -- so a loss
is a compute cost, not a re-seed. `spike/p13/realimage_report.jld2` (81,649 bytes) is **committed and
tracked** and carries every number in this document, so no table above depends on the gitignored
figure. The Phase-13 in-distribution reference is recomputed from the gitignored 50 MB training pool
on every run and is **not** persisted as a separate artifact; if that pool were lost, rebuilding it
would cost a full re-generation (13-11 reports ~10.9 min per 4,000 pairs), and because the pool is
counter-based Philox keyed on the global index that regeneration would also be byte-identical -- again
a compute cost, not a re-seed.

**10. This run created no bulk cache directory and wrote nothing under `spike/data/cache/`.** The
Phase-13 training pool was opened read-only through `load_pool` on each of its 12 shards; the cache
layer's `open_or_invalidate` was deliberately NOT called, so no directory could be created or
invalidated. `spike/data/cache/p11/` (54 MB, read-only to this phase, depended on by a Phase-12 plan)
was neither read nor written.

**11. A parallel agent is executing Phase 12 on this branch.** Every commit here used explicit paths
only (`git commit -- <paths>`); nothing was rebased or amended; `spike/test/runtests.jl` and the
`test_p12_*` files were left untouched. No index race occurred.

## Assumption Drift (advisory)

**A. The pre-registration's parenthetical summary-crossing figure for the `positive` fixture does not
match the measured ladder; the one for `negative` does.**
- **Found during:** Task 2, while building the section-2 crossing table.
- **Planned:** `spike/p13/consts.jl` (the `P13_REAL_ALPHA_GRID` comment block) states *"The measured
  summary-level sign crossing sits at alpha ~ 0.49 (positive) and ~ 0.46 (negative) -- strictly
  INTERIOR to this grid"*.
- **Actual:** measured on the committed fixtures through the frozen ingestion, `m-bar` crosses zero at
  **alpha ~ 0.674** on the positive fixture (between rungs 0.625 and 0.750) and **~0.464** on the
  negative one (between 0.375 and 0.500). The negative figure reproduces; the positive one does not.
  13-15 measured the same ladder independently and reported the same rungs.
- **Why it matters, and how far:** it is a **comment**, not a constant, and nothing is scored against
  it -- the executable pre-registration (`P13_REAL_ALPHA_GRID`, the invariant bounds, the substrate
  tuple) is unaffected, and the claim the comment was making ("strictly interior to this grid rather
  than pinned at an endpoint") remains TRUE on both fixtures. Most likely it was measured during
  research under a different construction. `consts.jl` was NOT edited: the pre-registration is
  append-only and a comment is not a bar.
- **Handling:** recorded here so a reader comparing the comment with the measured table sees the
  discrepancy stated rather than discovering it. Advisory only.

**B. 13-13 left open whether the epoch-4 best-validation checkpoint is sufficient for THIS arm; on
real substrate the answer is "it behaves coherently, but the OOD flag is the more binding
constraint".**
- **Planned:** 13-11 recorded the early overfitting (best validation at epoch 4) as an honest finding
  and deferred the question downstream; 13-12 answered it for the gate, 13-13 for the simulated
  ladder, and both explicitly left 13-16 open.
- **Actual:** the checkpoint produces a strictly monotone ladder on both heads and both conditions, a
  crossing point that lands exactly at the pre-registered hypothesis boundary on one fixture and one
  rung beyond it on the other, and a RANDOM verdict on both unmodified pairs that agrees with the
  D-05 label rule. No degeneracy, no saturation, no non-finite value.
- **Why it matters:** the checkpoint is not what limits this arm. **The OOD flag is** -- every read
  above is 2.49x outside the net's own training distribution, which bounds what any checkpoint of any
  quality could claim here. The "overfits early" finding stays on the record and stays honest, but it
  should not be carried forward as a caveat on any number in this document; the OOD caveat should be.

## Known Stubs

None. Every number in this document was read from `spike/p13/realimage_report.jld2` or, for the
`ghat(m-bar)` columns, derived from the persisted `m-bar` through the frozen `spike/simulator/ghat.jl`
map -- a read of frozen code, not a re-measurement and not a fit. No placeholder, mock or hardcoded
value flows into any table. The runner's two graceful-degradation branches (the `INVARIANT VIOLATION`
early return and the unavailable-reference paths) did not fire: both references were built and all
eight invariants held on both conditions.

## Threat Flags

None. This plan adds no network endpoint, no auth path and no schema at a trust boundary. It adds two
read surfaces, both already in the plan's threat register.

- The committed fixture directory under `test/` (T-13-50). Handled by the before/after git-index
  digest, the empty porcelain status, the mutating-call source grep, and the frame-size and
  element-type validation at ingestion.
- The frozen shipped bundle (T-13-57). Handled by opening it read-only for a comparison only, by
  printing and persisting that distinction, and by `p13_require_phase11()` running first with no
  fallback net.

## Requirements satisfied

**D-15 AMENDED** (the real-image arm run as a QUALITATIVE, REPORTED-NOT-GATED check; every log Bayes
factor recorded on the same row as its OOD verdict and its lambda; no pass/fail threshold defined for
any real-image quantity; the naming correction, the target-substitution record and the
qualitative-only statement carried in the banner, in the artifact metadata and in this SUMMARY; the
OOD finding measured and reported rather than suppressed, with both detectors side by side; the
sealed holdout left intact and the labelled check deferred to Phase 16),
**D-16** (the real alpha ladder read through the trained net as its OWN arm, with its invariants
re-verified on the reported substrate before any curve was read, `alpha_star_real` computed per
condition, the never-average ruling and the scope limit stated),
**D-03** (the registration-uncertainty response measured across the full pre-registered ladder and
reported honestly, INCLUDING its flatness as a finding),
**D-08** (evidence emitted as two log Bayes factors against a structural-zero random reference,
through the documented `three_way_log_bf` read surface with the PERSISTED per-head correction),
**D-01** (src/ decoupling asserted at RUN TIME; atomic persist with a reopen integrity check; the
read-only digest recorded before and after; CPU-only; no package installed; **no RNG stream consumed
at all**),
**D-04** (every pre-registered value READ from `consts.jl` and written INTO the artifact; nothing
inlined; `consts.jl` byte-unchanged; `P13_ITERATION_ALLOWANCE` untouched and UNSPENT).

## Commits

| Commit | Task | Description |
|--------|------|-------------|
| `2112bed` | Task 1 | read the real pairs -- the net says RANDOM on both, and both reads are OOD-flagged |
| _(this)_ | Task 2 | record the qualitative reading, the four honesty paragraphs and the liftable named limit |

Task 2 produced **no new artifact and no code change**, as the plan requires: it read
`spike/p13/realimage_report.jld2` and wrote this SUMMARY. The artifact is byte-identical to what
Task 1 committed.

## Self-Check: PASSED

- `spike/p13/run_p13_realimage.jl` -- FOUND (1,081 lines)
- `spike/p13/realimage_report.jld2` -- FOUND (81,649 bytes, 57 keys, **no `.tmp` sibling**)
- `spike/figures/p13_realimage.png` -- FOUND (294,821 bytes; gitignored, see Disclosure 9)
- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-16-SUMMARY.md` -- FOUND
- commit `2112bed` -- FOUND in history
- `git diff --quiet HEAD -- spike/p13/consts.jl` -- exit 0
- `git diff --quiet HEAD -- test/ src/ spike/Project.toml spike/Manifest.toml` -- exit 0
- `git status --porcelain test/test_images/` -- empty
- `spike/p13/consts.jl`, `spike/test/runtests.jl`, `spike/data/cache/p11/` -- not modified (absent
  from the commit's file list; the p11 cache was neither read nor written)
