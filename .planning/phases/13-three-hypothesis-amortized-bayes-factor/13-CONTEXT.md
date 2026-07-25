# Phase 13: Three-Hypothesis Amortized Bayes Factor - Context

**Gathered:** 2026-07-25 (discuss --analyze --all)
**Status:** Ready for planning

<domain>
## Phase Boundary

Phase 13 extends the amortized evidence network from two- to three-way model comparison —
**colocalized / random / mutually-exclusive** — so segregation becomes a first-class testable
hypothesis.

**Four scouting findings reshape this phase. Read them before planning.**

1. **Exclusion is already in the training distribution.** `GHAT_RHO_KNOTS` spans `[-0.99, +0.99]`
   and `RATIO_SPLIT_THRESHOLD = 0.0` (`src/amortized/train_ratio.jl:55`) splits `m=1 ⟺ Δρ > 0` vs
   `m=0 ⟺ Δρ ≤ 0`. The existing "null" class **already lumps random together with exclusion**. This
   phase is largely a relabelling plus a head swap — far cheaper than Phase 12, and no new physics.

2. **NeuralEstimators blocks a native 3-class head.** `RatioEstimator(net, 1; num_summaries)` takes
   `num_parameters = 1` (a *binary* model index), and per `train_ratio.jl:38-39` its loss is
   hard-coded `logitbinarycrossentropy` — **a passed loss is SILENTLY IGNORED in v0.2.1**. A 3-way
   softmax is not reachable by passing a loss. This is an architectural fork, not a tuning knob.

3. **The prior-atom problem lands hardest on the exclusion end.** `ghat` clamps asymmetrically
   (`GHAT_MU_MIN = -0.67976 → ρ = -0.99`; `GHAT_MU_MAX = 0.847149 → ρ = +0.99`). Under
   `Truncated(Cauchy(0, 0.3), -1, 1)` that is **≈4.9% of prior mass at the negative clamp versus
   ≈1.9% at the positive one** — about 2.5×. Named limit #2's atom artifact therefore concentrates on
   precisely the hypothesis this phase promotes to first-class.

4. **ROADMAP SC2 cites a retired reference.** It asks for reproduction of `compute_BayesFactor()`
   (`src/bayes.jl:109`), but named limit #3 records that the KDE baseline is **invalid past
   |logBF| ≈ log(L) ≈ 6.9**, that attrition is 100% baseline-side, and that the BF is now validated
   by simulation-based discrimination (spike 014, AUC 0.994). **SC2 is amended by D-12.**

**Out of scope:** the decision/abstention rule that turns evidence into a chosen hypothesis
(Phase 14), the spatial per-region map (Phase 12), and any corpus extension (Phase 8, COMPLETE).

</domain>

<decisions>
## Implementation Decisions

### Retraining posture and lane

- **D-01: Research lane.** Train the 3-way evidence net in `spike/` on a **fresh DEV seed** asserted
  disjoint from `PROD_SEED_V2`, `VAL_MASTER_SEED` (`0x5BC0FFEE`), `NPE_MASTER_SEED` (`0xC0FFEE`),
  `RATIO_PAIR_SEED` (`0x00000000004A7107`), `VAL_FIX_SEED` (`0xF1F7ED`), `CORPUS_MASTER_SEED`
  (`0x0000000000c05eed`) and all prior dev seeds. The shipped artifact, the shipped `bayes_factor`
  accessor, and the Phase-7 GO stay untouched.
  *Noted during discussion:* a reship would be materially cheaper here than in Phase 11/12 — Phase 13
  touches neither the NPE, the summary, nor any `DATAGEN_HASH_SRC_FILES` entry — but research lane
  was chosen to keep the GO intact.

- **D-02: Train on Phase-11's registration-aware `zt`,** not the shipped grid-8 `zt`.
  `train_ratio.jl:144` trains on pairs assembled from the FROZEN `zt` of an already-trained NPE, so
  the basis is inherited rather than created.
  **⚠ CONSEQUENCE FOR PLANNING: this adds a Phase 11 dependency that ROADMAP does not declare**
  (ROADMAP lists Phase 7 only). Phase 13 *planning* may proceed now; **execution blocks on Phase 11's
  research net existing.** It also weakens the SC2 overlap check (different input surface) — handled
  by D-12.

- **D-03: Condition the evidence net on the registration-uncertainty level** (Phase 11 D-03's extra
  input). The log-BF pair weakens as registration uncertainty grows. Rationale: Phase 11 D-07 commits
  to posterior width growing with uncertainty; a Bayes factor that stayed confident while the
  posterior widened would be an internal contradiction inside one codebase.

- **D-04: Pre-register thresholds, one documented iteration.** Carried forward from Phase 11 D-04 and
  **not re-asked**: lock every Phase-13 threshold in a Phase-13 consts file with the fresh DEV seed
  **before any run**, and declare the single-iteration allowance in advance.

### Hypothesis boundaries

- **D-05: Two-factor cut — `ρ_sample` sign × control contrast.** Exclusion is claimed only when the
  sample is genuinely anti-correlated (`ρ_sample < -τ`) **and** stands apart from the control; coloc
  when `ρ_sample > +τ` and above control; random otherwise.
  **Why not the Δρ contrast:** `infer.jl:117-121` defines `Δρ = mean(ρ_sample_draws .- ρ_control_draws)`
  — a sample-vs-control contrast. A sample at `ρ = 0.8` under a control at `ρ = 0.9` has `Δρ < 0`.
  The binary net calls that "null" (harmless); a naive three-way extension on the same quantity would
  publish it as **"mutually exclusive"** while the sample is strongly colocalized. "Less colocalized
  than the control" is not segregation.

- **D-06: τ comes from measured summary resolution.** Probe the smallest ρ difference the fixed
  patch-correlation summary can reliably distinguish and set τ to it, so "random" means
  *indistinguishable from zero at this summary's resolution*. Same measure-don't-guess posture as
  Phase 11 D-06. **τ must be frozen in the Phase-13 consts file before any labelled data is
  generated.**

- **D-07: Stratify training coverage; restore scale with a 3-way prior-odds correction.** Build the
  evidence net's training set with stratified coverage across `ρ ∈ [-0.99, -τ]` rather than sampling
  straight from π(θ), then restore the correct evidence scale at read time via a **3-way
  generalization of `measure_log_prior_odds`** — the mechanism `bf.jl:74`
  (`amortized_log_bf(est_bf, Z_pair, log_prior_odds)`) already uses for the binary case.
  **π(θ) is untouched**, so the CLAUDE.md constraint that π(θ) match the Turing μ/ν/σ/τ ranges holds.
  **The 3-way correction must be DERIVED AND VERIFIED, not assumed to generalize.**
  Rejected: μ-prior truncation (spike 013) — VALIDATED *with cost*, it breaks that same
  prior-consistency constraint and degrades high-|ρ| coverage.

- **D-08: Output = log-BFs against the random reference.** Emit `log BF(coloc : random)` and
  `log BF(exclusion : random)`; the random entry is identically 0. Same semantics as the shipped
  binary `bayes_factor`, and **no model prior over the three hypotheses has to be invented**
  (turning evidence into a chosen hypothesis is Phase 14's job).
  **API note:** the sketched field name `log_bf_simplex` (`src/results.jl:172`) becomes a misnomer —
  it is not a simplex. Rename it when implementing `ThreeHypothesisColocResult`.

### Three-way head architecture

- **D-09: Shared trunk, two BCE heads.** One Flux conditioner feeding two binary-cross-entropy heads
  (coloc-vs-random, exclusion-vs-random). Satisfies SC1's "one forward pass", preserves the logit/BCE
  semantics that make the D-07 prior-odds correction valid, and the shared conditioner makes the two
  evidences structurally consistent. **Cost accepted:** it lives outside `RatioEstimator`, inherits
  no gate lineage from the shipped binary net, and must be validated on its own terms (D-12/D-13).

- **D-10: Reuse the exact existing trunk.** `RATIO_SUMMARY_WIDTH = 256`, `RATIO_NUM_SUMMARIES = 64`,
  topology `Chain(Dense(in,256,gelu), Dense(256,256,gelu), Dense(256,256,gelu), Dense(256,64))`
  verbatim. Only the head and the input width change (`ratio_input_dim(G) = 5·G² = 320` at G=8,
  **plus 1** for the D-03 uncertainty input). Any difference in the 2-way overlap is then attributable
  to the head alone. Backed by spike 012 (**PARTIAL**): capacity *redistributes* error rather than
  removing it, so capacity is not a reliable lever in this codebase.

### Validation target (SC2 amended)

- **D-12: Gate on simulation-based 3-way discrimination; report binary-NRE continuity.**
  **Gate:** per-class AUC plus a full confusion matrix over {exclusion, random, coloc} at known
  ground-truth labels — the method spike 014 validated and the GO adopted.
  **Reported, not gated:** correspondence with the **shipped binary NRE** in the overlapping 2-way
  regime at a fixed uncertainty level, as a continuity check.
  **ROADMAP SC2 is amended** — it must NOT be validated against `compute_BayesFactor()`, which named
  limit #3 established is invalid past |logBF| ≈ 6.9, precisely the confident-call tail where
  three-way verdicts matter most. **Document the amendment up front, before results.**

- **D-13: Also gate evidence calibration.** Reliability of the implied class probabilities
  (ECE/MCE via the existing `_bin_calibration` / `CalibrationResult` traffic-light), populating the
  `CalibrationMeta` field already sketched on `ThreeHypothesisColocResult`.
  **Why this does not contradict named limit #3:** that retreat was forced by a broken *reference*
  (the KDE baseline), not by a finding that NRE magnitudes are wrong — nothing was ever measured.
  Classifier calibration is checked directly against ground-truth labels with **no external
  reference**, so the obstacle does not apply. Phase 14's abstention layer needs validated magnitudes.

### Exclusion ground truth (SC3)

- **D-15 (AMENDED 2026-07-25): Simulator gate + semi-synthetic graded series + a `test/test_images/`
  qualitative check. The corpus physical anchor is NOT used.**

  > **Amendment rationale — a correction to the original audit.** The first version of this decision
  > named the physical mitochondria anchor as the real-data check. Verification during Phase 11
  > planning showed the audit was incomplete: it counted the anchors correctly but did not read their
  > `split` and `sha256` columns. Both `physical-primary` rows are `sha256 = PENDING-FETCH`,
  > `bytes = 0`, **`split = sealed_holdout`** — unfetched, and reserved for the **Phase-16 blind
  > evaluation**. Consuming the segregated anchor here would irreversibly burn Phase 16 on the very
  > hypothesis Phase 16 would be evaluating. Mirrors the ruling Phase 11's planner made for its D-17
  > and the amended Phase 12 D-11.

  **Corpus audit (still valid, and still the reason the semi-synthetic series is load-bearing):** the
  corpus holds exactly one segregated anchor and one coloc anchor, plus 30 `simulated-secondary` CBS
  Red-Green rows whose degrees are colocalization *percentages* — so the CBS `0.0` end is **random, not
  exclusion**. **The corpus has no graded segregation series at all**, and now also no *usable* real
  segregation anchor for this phase.

  Therefore: **gate** on simulator ground truth (`ρ_true < -τ`, exact labels, well-powered); add the
  D-16 semi-synthetic graded random→exclusion series as supporting evidence — it is now the **only**
  graded segregation evidence that will exist; and run the qualitative real-data check on the six
  committed microscopy TIFFs in `test/test_images/` (`positive_c{1,2,3}.tif`, `negative_c{1,2,3}.tif`,
  1028×1376).

  **Honesty consequence to state in the report:** those images carry no colocalization ground-truth
  label, so the real-data check is **qualitative only** — it can show that exclusion verdicts behave
  sensibly on real microscopy, not that they are correct. Real segregation validation is deferred to
  Phase 16's blind evaluation, which is precisely what the sealed holdout exists for.

  **Do NOT touch `corpus/` sealed_holdout rows in this phase.**
  **Conceptual gap this addresses:** the simulator's "exclusion" is *negative intensity correlation
  across patches*; a biologist's "mutually exclusive" is *disjoint spatial localization*. These
  coincide often but not necessarily, and the semi-synthetic series is the only evidence that probes
  it.

- **D-16: Semi-synthetic series = mask-based disjoint reassignment.** Threshold ch1 into an object
  mask and progressively redistribute ch2 signal into its complement, graded by a mixing parameter
  `α ∈ [0, 1]` from random (α=0, unchanged) to fully disjoint (α=1). Segregation is constructed
  **spatially**, so the series tests the simulator's assumption rather than restating it.
  **Rejected: displacement of ch2 objects** — a global displacement is exactly the mis-registration
  nuisance Phase 11 marginalizes over, so the net could read the series as shift; a null result would
  be uninterpretable and a positive result unattributable.
  **Rejected: intensity anti-correlation** — reproduces the simulator's own mechanism, so it tests
  the net against its own training assumption (circular) and leaves the localization gap untouched.
  The segmentation/threshold rule is a pre-registered choice (D-04 applies).

### Claude's Discretion

- **D-11: Per-head class restriction and joint trunk training.** Each head trains **only on its own
  pair of classes** (coloc head sees coloc + random; exclusion head sees exclusion + random). One-vs-rest
  training would make the logit a coloc-vs-*mixture* ratio rather than `log BF(coloc : random)`,
  silently breaking D-08. The **trunk trains jointly** (both head losses summed) so it sees all three
  classes — that is what makes the two evidences structurally consistent per D-09.
  *Decided by Claude as following necessarily from D-08/D-09; user may override.*

- **D-14: No over-power workaround needed for the calibration gate.** `_bin_calibration` /
  `CalibrationResult` score ECE/MCE against a **traffic-light band**, not a p-value, so the D-13 gate
  cannot be over-powered the way Phase 7's M=2000 SBC point-null was. No TOST/equivalence machinery
  is required here (unlike Phase 11 D-08). *Decided by Claude; user may override.*

- **Sequencing.** Phase 13 planning proceeds now; execution blocks on Phase 11's research net
  existing (see D-02). *Decided by Claude as GSD-internal sequencing.*

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### The evidence-net surface this phase extends
- `src/amortized/train_ratio.jl` — `build_ratio_estimator`, `RATIO_SPLIT_THRESHOLD = 0.0` (:55),
  `RATIO_SUMMARY_WIDTH = 256` / `RATIO_NUM_SUMMARIES = 64` (:53-54),
  `RATIO_PAIR_SEED = 0x00000000004A7107` (:60), `measure_log_prior_odds`.
  **:38-39 — the hard-coded-loss landmine: a passed loss is SILENTLY IGNORED in v0.2.1.**
- `src/amortized/bf.jl` — `amortized_log_bf(est_bf, Z_pair, log_prior_odds)` (:62-74), `pair_encode`,
  `kde_log_bf_unclamped` (:128-148, the retired baseline), and the §35 honesty note on the
  `RatioEstimator` loss.
- `src/amortized/infer.jl:94-121` — `rho_draws` and `delta_rho`; **the Δρ definition D-05 turns on**.
- `src/amortized/api.jl:140-145` — how `logbf` is composed into `AmortizedColocResult` today.
- `src/amortized/summary.jl` — `ratio_input_dim(G) = 5·G²`, `_summary_row_partition`.

### Result-type extension point
- `src/results.jl:161-192` — the SKETCH-ONLY extension block. `ThreeHypothesisColocResult` at :169-178
  with `log_bf_simplex :: NTuple{3,Float64}` and `bayes_factor_simplex`. **D-08 renames the field.**
  Phase 7 D-02 requires new variants to slot in as new `AbstractColocResult` subtypes, never as
  bolted-on fields.

### Prior / atom machinery (D-06, D-07)
- `spike/simulator/ghat.jl:42-49` — `GHAT_MU_KNOTS`, `GHAT_RHO_KNOTS`, `GHAT_MU_MIN = -0.67976`,
  `GHAT_MU_MAX = 0.847149`, and the clamping behaviour that creates the asymmetric atoms.
- `spike/simulator/prior.jl:40, 71-82` — `MU_PRIOR = Truncated(Cauchy(0,0.3),-1,1)`,
  `ρ_true = ghat(μ*)`, and the SIM-02 prior-consistency contract.
- `src/bayes.jl` ~274-291 — the Turing μ-prior that π(θ) must stay consistent with (CLAUDE.md
  constraint). `src/bayes.jl:109` — `compute_BayesFactor()`, the reference **D-12 declines to gate on**.

### Honesty contract and prior findings
- `docs/amortized.md` §Named limits — limit #2 (ρ_true atoms / randomized ranks) and limit #3
  (BF validated by discrimination, KDE invalid past |logBF| ≈ 6.9). Both are load-bearing here.
- Skill `spike-findings-proteincoloc` → `references/bayes-factor-gate.md` (spikes 006-008, why the KDE
  reference is invalid) and `references/sbc-calibration.md` (spike 012 PARTIAL — capacity
  redistributes; spike 013 — μ-truncation costs prior consistency).
- `.planning/phases/07-productionization-conditional-on-go/07-CONTEXT.md` — D-02 result hierarchy,
  D-05 per-grid re-pre-registered ship gate, D-06 GPU-optional/CPU-reproducible.

### Upstream dependency introduced by D-02
- `.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-CONTEXT.md` — D-01 research
  lane and seed set, **D-03 the uncertainty conditioning input this phase inherits**, D-06 the probe
  posture D-06 here mirrors, D-07/D-08 the widening-and-coverage commitments D-03 here must not
  contradict.

### Ground-truth data (D-15 AMENDED, D-16)
- `test/test_images/positive/positive_c{1,2,3}.tif`, `test/test_images/negative/negative_c{1,2,3}.tif`
  — **the real images D-15 now uses** for the qualitative check, and the source images D-16's
  mask-based α-graded series is constructed from. Six committed microscopy TIFFs (1028×1376), already
  exercised by `test/runtests.jl`. **No colocalization ground-truth labels.**
- `corpus/manifest.csv` — read for provenance only. **Both `physical-primary` rows are
  `split = sealed_holdout`, `sha256 = PENDING-FETCH`, `bytes = 0`** — reserved for the Phase-16 blind
  evaluation and **must not be consumed by this phase**. `CORPUS_MASTER_SEED = 0x0000000000c05eed`.
- `spike/validation/consts.jl:72-79` — the pre-registered seed constants and the disjointness
  rationale; the file pattern D-04's Phase-13 consts file should mirror.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`measure_log_prior_odds` + `amortized_log_bf`** — the training-balance-to-evidence correction that
  makes D-07's stratification legitimate. Already built and shipped for the binary case; D-07 requires
  its 3-way generalization.
- **The `build_ratio_estimator` trunk topology** — reused verbatim per D-10; only the head changes.
- **`_bin_calibration` / `CalibrationResult`** (ECE/MCE traffic-light) and the `CalibrationMeta` field
  — D-13's calibration gate should populate these rather than invent a new diagnostic surface.
- **`pair_encode`** (`bf.jl`) — grid-general (sample, control) summary pairing, independent of any
  cache/loader; reusable unchanged.
- **`spike/validation/consts.jl`** — the pre-registration file pattern (named seed constants with
  explicit `≠ <other seed>` disjointness assertions) for D-04.

### Established Patterns
- **Evidence semantics are logit-based.** Everything downstream (`amortized_log_bf`, prior-odds
  correction) assumes a binary-cross-entropy logit. D-09's two-BCE-head design preserves this; a
  softmax head would have broken it, which is why D-08 rejected model probabilities.
- **New result variants are new subtypes** (Phase 7 D-02), never new fields on existing structs.
- **Frozen `zt` discipline.** The ratio net always trains on an already-frozen standardized summary
  basis from a trained NPE — it never creates its own.
- **Pre-register before running.** Every gate in this project locked thresholds and seeds first;
  D-04 continues it, D-12's SC2 amendment must be declared *before* results.
- **Named limits are the honesty mechanism.** The registration-robustness gap and the n=1 physical
  segregation anchor both belong on that list.

### Integration Points
- `src/amortized/api.jl:140-145` — where a three-hypothesis result would eventually compose, if a
  later phase productionizes this. **Not touched in Phase 13** (D-01 research lane).
- `src/results.jl:169-178` — the `ThreeHypothesisColocResult` slot, with the D-08 field rename.
- Phase 11's research-net output (summaries + uncertainty conditioning input) is this phase's input
  contract (D-02, D-03).
- Phase 14 consumes D-13's calibrated magnitudes as its abstention thresholds.

</code_context>

<specifics>
## Specific Ideas

- The correctness trap in D-05 should be stated explicitly in the manuscript, not just the code: a
  three-way test cut on the sample-vs-control contrast would label "less colocalized than control" as
  "mutually exclusive". Naming the trap is itself a contribution, since a reader may well assume the
  naive extension.
- D-16's α-graded series is wanted specifically because it probes the **correlation-vs-localization**
  gap. The interesting output is not just "does exclusion get detected" but "at what α does the
  intensity-correlation notion start to track disjoint localization".
- The `log_bf_simplex` rename (D-08) is a small thing but should happen at implementation time — a
  field named "simplex" that is not a simplex will be quoted as one.

</specifics>

<deferred>
## Deferred Ideas

- **Reshipping the 3-way evidence net.** Noted as materially cheaper here than in Phase 11/12 (no NPE
  retrain, no summary change, no `DATAGEN_HASH_SRC_FILES` edit), but deferred to keep the Phase-7 GO
  intact. A future productionization phase should revisit it.
- **Extending the corpus with real graded segregation anchors.** The honest fix for the n=1 problem in
  D-15, but Phase 8 is COMPLETE — reopening it is a roadmap change, not a Phase-13 decision.
- **μ-prior truncation to remove the atoms at source** (spike 013). VALIDATED with a real cost to
  π(θ)/Turing consistency and high-|ρ| coverage; a separately-decided change.
- **Full pairwise log-BF matrix** including a direct coloc-vs-exclusion contrast. Rejected for D-08
  (only two of three are independent), but the coloc-vs-exclusion comparison may be worth deriving and
  reporting in the manuscript.
- **Posterior model probabilities over the three hypotheses.** Requires model priors that do not exist;
  properly belongs to Phase 14's decision layer if anywhere.
- **Training the 3-way net on the shipped grid-8 basis as a second, input-comparable net.** Would give
  an exact 2-way overlap comparison; rejected as a second training run, but it is the clean way to
  strengthen D-12's continuity check if it proves weak.

</deferred>

---

*Phase: 13-three-hypothesis-amortized-bayes-factor*
*Context gathered: 2026-07-25*
