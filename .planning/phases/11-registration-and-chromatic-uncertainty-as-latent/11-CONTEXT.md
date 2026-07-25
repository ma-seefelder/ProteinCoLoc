# Phase 11: Registration and Chromatic Uncertainty as Latent - Context

**Gathered:** 2026-07-25
**Status:** Ready for planning

<domain>
## Phase Boundary

Phase 11 delivers a **research-lane** demonstration that the amortized coloc posterior widens
*honestly* under registration and chromatic uncertainty, and that this widening preserves
credible-interval coverage up to a measured breakdown point.

**Scope correction — read before planning.** ROADMAP SC1 ("dx/dy is added to θ … and the NPE is
retrained") is **already nominally satisfied**. `shift_dx`/`shift_dy` are sampled from
`SHIFT_PRIOR = Uniform(-1, 1)` in `spike/simulator/prior.jl`, applied in `spike/simulator/forward.jl`
stage 6, mirrored in `src/amortized/simulator.jl`, and read back off the NPE at
`src/amortized/ood.jl:129-130`. Planning must NOT re-implement this.

The three genuine gaps Phase 11 owns:

1. **No chromatic term exists.** `forward.jl` stage 6 is translation-only.
2. **The SC2 monotone-widening sweep was never run.**
3. **SC3 mis-registration handling is unaddressed.**

What Phase 11 explicitly does NOT do: no retraining of a shipped model, no new artifact, no
re-opened ship gate, no public API change, and no claim that `shift_dx`/`shift_dy`/`ε` are
*identified*.

</domain>

<decisions>
## Implementation Decisions

### Retraining posture

- **D-01: Research lane only.** Train a Phase-11 research net in `spike/` on a **fresh DEV seed**
  asserted disjoint from `PROD_SEED_V2`, `VAL_MASTER_SEED` (`0x5BC0FFEE`), `NPE_MASTER_SEED`
  (`0xC0FFEE`), `VAL_FIX_SEED` (`0xF1F7ED`), `CORPUS_MASTER_SEED` (`0x0000000000c05eed`) and all
  prior dev seeds. The shipped `amended_v2/grid_8` artifact and the Phase-7 GO (Option A,
  2026-07-24) stay untouched — no reship, no gate reopened. Science advances without spending
  pre-registration credibility.

- **D-02: Widen `SHIFT_PRIOR` to `Uniform(-3, 3)` px** for the research net (from `Uniform(-1, 1)`).
  Covers realistic multi-channel error — chromatic aberration plus stage/filter-cube drift — and
  gives the SC2 ladder real headroom (0.25 → 3 px). **A wider nuisance prior adds marginal spread
  that may attenuate `ρ_true` inference; this must be MEASURED against the frozen net, not assumed.**

- **D-03: Amortize over the uncertainty level.** Append the registration-uncertainty level (the
  shift prior's half-width) to the network input so one training run yields a continuous, controlled
  sweep — feed the same image at level 0.25 → 3 px and read posterior width directly. **Consequence:
  the input dimension changes, so this net is NOT drop-in comparable to the shipped read surface.**

- **D-04: Pre-register, allow one iteration.** Lock the SC2 monotonicity and SC3 no-overconfidence
  thresholds in a Phase-11 consts file (pattern: `spike/validation/consts.jl`) with the fresh DEV
  seed **BEFORE any run**, and state up front that this is a research-lane result permitting **one
  documented iteration**. Naming the allowance in advance avoids repeating the Phase-7
  "amended twice, credibility spent" pattern.

### Identifiability strategy

- **D-05: Marginalize only — no identification claim.** The patch-correlation summary stays
  unchanged. `gate-8x8-amended.md:42` lists
  `vacuous_params = [spillover, autofluorescence, shift_dx, shift_dy, noise]` at shrinkage 0.92–1.03
  — the posterior IS the prior. Widening the shift prior propagates registration uncertainty into
  `ρ_true`/`Δρ` as honest marginal spread; the shift columns (and `ε`) remain vacuous and are
  **reported as such**. Per `gate-8x8-amended.md:46`, a clean KS p on a vacuous column is
  meaningless and must not be quoted as calibration evidence.
  **SC2 is a claim about the COLOC posterior widening, not about recovering dx/dy/ε.**

- **D-06: Simulator-only pre-flight probe BEFORE training.** Risk being controlled: at grid 8 on
  256² images each patch is ~32 px, so ±3 px is ~9% displacement; after PSF blur the summary may
  barely move and SC2 would fail *flat* — below resolution, not wrong design. Probe = sweep
  shift `0 → 3` px and `ε` (see D-09) through `forward.jl` at fixed θ and measure how far the
  128-dim summary actually moves, relative to the movement induced by a known `Δρ`. Minutes, no
  training. **The probe supplies the SC2 threshold; the threshold is not guessed.**

- **D-07: SC2 pass criterion is two-part.** (a) Posterior SD (or 90% HDI width) of `Δρ` is
  **monotone non-decreasing** across the uncertainty ladder, AND (b) empirical coverage of the 90%
  credible interval stays within nominal at **every** rung. Part (b) is what makes the widening
  *honest* rather than merely larger — it rules out a net that inflates variance without tracking
  truth, and it is what the ROADMAP goal ("widens honestly instead of reporting false confidence")
  actually asserts.

- **D-08: Per-rung coverage tested by equivalence, not by a point null.** TOST-style test for
  equivalence to nominal with a **pre-registered tolerance band** (e.g. ±3 percentage points on the
  90% CI), **Holm-corrected across rungs**. This is the exact remedy the Phase-7 spikes recommended
  after the M=2000 over-power failure (a ~0.06-SD drift rejected by a strict point-null). Passing
  then means "coverage is provably close to nominal", not "we failed to reject".
  **The tolerance band must be justified from the D-06 probe, not chosen after seeing results.**

### Chromatic warp scope

- **D-09: Add a 1-parameter radial chromatic scale `ε`,** prior `Uniform(-0.02, 0.02)`, symmetric.
  Rationale: lateral chromatic aberration is predominantly a **radial magnification difference** —
  zero at the image centre, maximal at the corners — so a global `dx`/`dy` translation *cannot*
  represent it. It is a distinct failure mode, not a reparameterisation. On a 256² field the
  half-diagonal is ~181 px, so corner displacement ≈ `ε × 181`: 0.005 → ~0.9 px (well-corrected
  apochromat), 0.01 → ~1.8 px (typical achromat), 0.02 → ~3.6 px (poorly corrected). The ±0.02 range
  matches the ±3 px translation headroom so both misalignment axes are comparable in magnitude.
  **Symmetric because channel ordering is the user's choice, not the model's** — a one-sided prior
  would silently assume ch2 is always the longer wavelength, which is a correctness trap the API
  does not enforce. `ε` joins θ with its own entry in `prior.jl`.
  Note: `ε` is expected to be another **vacuous** column (D-05 applies to it too).

- **D-10: Single composed affine warp.** Build one `AffineMap = Translation(dy, dx) ∘ LinearMap(
  scale about image centre)` and pass it to a **single** `warp` call, replacing the current
  `Translation`-only call at `forward.jl:168-174`. **Rationale is a correctness issue, not style:**
  every resampling pass applies interpolation smoothing, and smoothing decorrelates — a second
  `warp` would widen the posterior for a reason unrelated to misalignment, contaminating exactly
  the quantity SC2 measures. One interpolation pass also keeps the `ε = 0, shift = 0` case
  byte-comparable with today's simulator, which gives a **regression check that the refactor changed
  nothing**. The `WR-06` axis-order comment must be **re-derived** for the composed map, since the
  `(dy, dx)` argument convention now interacts with the scale centre.

- **D-11: Mirror the `ε` extension into `src/amortized/simulator.jl` now** (not spike-only), so the
  two simulators do not diverge.

- **D-12: Record the provenance consequence of D-11 as a named limit.** `src/amortized/simulator.jl`
  is listed in `DATAGEN_HASH_SRC_FILES` (`src/amortized/datagen.jl:208-213`); its source **bytes**
  feed the content-hash digest that names the training-data cache directory.
  **Accurate scope of the cost:** the shipped artifact is pinned by `git-tree-sha1` in
  `Artifacts.toml` (`90e6b63a…`), **not** by the datagen digest — so editing `simulator.jl` does
  **not** break the artifact download or `colocalization_amortized(...)`. What it breaks is
  narrower: it orphans existing training-data caches and severs the claim that current `src/`
  reproduces the training distribution the shipped net came from.
  **Repair:** add a named limit to `docs/amortized.md` stating that the shipped `grid_8` bundle was
  trained by `src/amortized/simulator.jl` **as of commit `<sha>` (pre-`ε`)**, recording that sha
  explicitly so the training distribution stays exactly recoverable via git. Provenance is repaired
  **by reference rather than by bytes**, consistent with how named limits #1–#4 are already handled.

### Mis-registration handling

- **D-13: SC3 is a coverage claim, not a flagging claim.** "Handled without silent overconfidence"
  means: for injected misalignment, the `Δρ` credible interval must still cover the truth. **No OOD
  flag requirement.** A covered interval *is* the absence of silent overconfidence.
  **Structural tension this resolves (must be stated in the report):** the shipped OOD density
  channel fires when a summary looks unlike the *training* summaries; widening the shift prior
  (D-02) and adding `ε` (D-09) deliberately pull mis-registered summaries *into* the training
  distribution. **The wider the prior, the less mis-registration is detectable as OOD** — SC2's
  mechanism directly erodes any density-channel detector. SC2 and SC3 are not independent criteria.
  **Accepted consequence:** badly mis-registered data returns a wide-but-unflagged answer; the user
  is not told *why* the interval is wide. The failure is honest but silent.

- **D-14: Report the breakdown curve.** Test coverage well past the training prior (0 → 8 px,
  `|ε| → 0.05`) and report the point where 90% coverage first drops below nominal. **SC3 passes if
  coverage holds throughout the training prior; the breakdown point beyond it is REPORTED, not
  gated.** This yields a directly quotable limit — "coverage holds to X px, degrades beyond" — which
  pre-empts the obvious reviewer question rather than inviting it, and costs extra evaluation rungs
  only, no extra training.

- **D-16: Deliverable = report + a `docs/amortized.md` interpretation note.** Phase 11 ships a
  results report (probe, SC2 ladder, SC3 breakdown curve, figures) under `.planning/`, **plus** a
  short interpretation section in `docs/amortized.md` telling users how to read posterior width
  under registration uncertainty and quoting the breakdown point. **No API change, no artifact
  change.** The wording must be scrupulous about which model each number came from, since the
  section describes a **research** net that is not the shipped one.

- **D-17: Validation is simulation + a qualitative real-image check.** Primary evidence is simulated
  data with known ground-truth misalignment. **Additionally**, on the existing real images already
  used by the manuscript pipelines (`corpus/`, pre-registered manifest D-08,
  `CORPUS_MASTER_SEED = 0x0000000000c05eed`), deliberately shift one channel by a known amount and
  confirm the posterior widens as the simulated sweep predicts. There is no ground-truth `ρ` on real
  data, so **this cannot be a coverage test** — but a matching width response on real images is
  strong evidence the simulator's degradation mechanism transfers. **Read-only access to those
  images; the manuscript pipelines must be provably unmodified, and the report must say so.**

### Claude's Discretion

- **D-15: Misalignment for SC3 is simulator-injected, not post-hoc.** Force `shift`/`ε` to the test
  value and run the single composed affine warp (D-10); do **not** post-hoc warp already-simulated
  images. Post-hoc injection would add the second interpolation pass D-10 exists to eliminate, and
  its smoothing would decorrelate independently of misalignment, contaminating the coverage curve.
  *Decided by Claude for consistency with D-10; user may override.*

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Simulator and prior (the files this phase edits)
- `spike/simulator/prior.jl` — `SHIFT_PRIOR = Uniform(-1, 1)`, `sample_prior(rng)` returning the
  7-field θ NamedTuple. D-02 widens `SHIFT_PRIOR`; D-09 adds `ε`.
- `spike/simulator/forward.jl` §stage 6 (lines 168-174) — the current `Translation(θ.shift_dy,
  θ.shift_dx)` warp and the `WR-06` axis-order comment that D-10 must re-derive.
- `src/amortized/simulator.jl` — the productionized mirror (`SHIFT_PRIOR` :94, sampling :116-117,
  bounds tuple :145-146, warp :234-238). D-11 mirrors `ε` here; D-12 covers the consequence.

### Provenance and artifact machinery (why D-11/D-12 are constrained)
- `src/amortized/datagen.jl:199-224` — `DATAGEN_HASH_SRC_FILES` content-hash version guard.
  `simulator.jl` is a hashed data-defining source; any byte change re-digests the cache directory.
- `src/registry.jl:145-161` — artifact resolution and `_activate_recorded_ood_threshold`.
  **Line 158: "The persisted `ood_nulls` carries only the `:density` channel."**
- `Artifacts.toml` — the shipped `amended_v2/grid_8` pin (`git-tree-sha1 = 90e6b63a…`, `lazy = true`).

### OOD read surface (background for D-13)
- `src/amortized/ood.jl` — three implemented channels (`:density` Mahalanobis on continuous summary
  rows, `:pp` posterior-predictive, `:noise` image-noise), `id_threshold`, `OOD_ID_QUANTILE = 0.95`,
  fused `ood_verdict`. Lines 129-130 read `shift_dx`/`shift_dy` off the NPE output.
  **Only `:density` is fitted in the shipped bundle.**

### Calibration and gate history (why D-05/D-08 are what they are)
- `.planning/phases/07-productionization-conditional-on-go/gate-8x8-amended.md` — §line 42
  `vacuous_params` list and shrinkage values; §line 46 why `shift_dx`'s clean KS p is meaningless;
  §line 125 "5/8 parameters are vacuous … any future SBC pass on those columns must be read as
  vacuous."
- `.planning/phases/07-productionization-conditional-on-go/07-CALIBRATION-FINDINGS.md` — finding F3
  (vacuous pass ≠ calibration evidence), `shrinkage = post_sd / prior_sd` diagnostic definition.
- `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md` — §6.4 blocks
  iterating the current amended_v2 gate; §121 shrinkage/vacuous diagnostics are reporting-only.
- `.planning/phases/07-productionization-conditional-on-go/07-NUISANCE-SBC-SPEC-DRAFT.md` — the
  nuisance-appropriate equivalence-test design D-08 follows.
- Skill `spike-findings-proteincoloc` → `references/sbc-calibration.md` — randomized-rank recipe,
  M=2000 over-power diagnosis, why capacity/truncation levers do not fix nuisance drift
  (spikes 012, 013 are PARTIAL).

### Phase-7 posture that D-01 must not violate
- `.planning/phases/07-productionization-conditional-on-go/07-CONTEXT.md` — D-01 breaking
  amortized-only public API; D-02 `AbstractColocResult` hierarchy; D-05 per-grid re-pre-registered
  ship gate; D-06 GPU-optional / CPU-reproducible.
- `docs/amortized.md` §Named limits (lines ~106-135) — limit #4 (nuisance-marginal drift on
  parameters the fixed summary cannot constrain) is the direct antecedent of D-05, and the section
  D-12 extends.
- User memory `v2-go-decision-named-limits` — the 2026-07-24 GO (Option A): no further training,
  no gate iteration.

### Seed discipline (D-01)
- `spike/validation/consts.jl:72-79` — `VAL_MASTER_SEED = 0x5BC0FFEE`, `VAL_FIX_SEED = 0xF1F7ED`,
  and the D-02 disjointness rationale. Also the file pattern D-04's Phase-11 consts file follows.
- `spike/npe/train_npe.jl:65` — `NPE_MASTER_SEED = 0xC0FFEE`.
- `corpus/manifest.csv` header — `CORPUS_MASTER_SEED = 0x0000000000c05eed` (D-17 target data).
- `spike/data/seeding.jl` — Philox-per-index counter-based seeding helper.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`spike/simulator/forward.jl` stage pipeline** — `ε` slots into the existing stage 6 rather than
  adding a stage; the surrounding PSF/spillover/autofluorescence/noise stages are untouched.
- **`spike/validation/consts.jl`** — the established pre-registration pattern (named seed constants
  + rationale comments in-file). D-04's Phase-11 consts file should mirror its structure, including
  the explicit `≠ <other seed>` disjointness assertions.
- **`spike/data/seeding.jl`** — counter-based Philox seeding, already reproducible per index.
- **`src/amortized/ood.jl`** `roc_auc` / `id_threshold` / `youden_j` helpers — available if the
  breakdown curve (D-14) wants a separability statistic, though D-13 removes the flag requirement.
- **Existing SBC/coverage harness** (`spike/validation/harness.jl`) — the per-rung coverage runs of
  D-07/D-08 should extend it rather than build a parallel harness.

### Established Patterns
- **Vacuous-column reporting (F3).** `shrinkage = post_sd / prior_sd` with a boolean `vacuous` flag
  is already a reporting-only diagnostic. `ε` and the widened `shift_*` columns must be reported
  through it, and any uniform rank result on them explicitly labelled non-evidence.
- **Named limits as the honesty mechanism.** `docs/amortized.md` already carries limits #1–#4 in a
  consistent voice; D-12 and D-16 both extend that list rather than inventing a new surface.
- **Content-hash provenance.** Any edit under `DATAGEN_HASH_SRC_FILES` is a provenance event, not a
  routine change — planning must treat D-11 as such and sequence D-12 in the same commit.
- **Pre-registration before execution.** Every prior gate in this project locked thresholds and
  seeds in a consts file before running. D-04 continues it; the one-iteration allowance is the only
  deviation and it is declared in advance.

### Integration Points
- `forward.jl` stage 6 ↔ `src/amortized/simulator.jl:234-238` — the two warps must stay
  behaviourally identical after D-10/D-11.
- `prior.jl` θ NamedTuple ↔ `src/amortized/ood.jl:129-130` θ read-back ↔ `theta_prior_bounds()`
  (consumed by `src/amortized/architecture.jl:65`) — adding `ε` changes θ arity, so the bounds tuple
  and the flow's `D` (default 7, `architecture.jl:187`) both move.
- **D-03's input-dimension change is confined to the research net.** It must not propagate into the
  shipped `EstimatorBundle` read surface.
- `corpus/load.jl` / `corpus/fetch.jl` — the read-only path for D-17's real-image check.

</code_context>

<specifics>
## Specific Ideas

- The chromatic term is deliberately chosen as the **only spatially-varying misalignment** in the
  model — corner patches decorrelate while centre patches do not. This makes it the natural bridge
  to **Phase 12's calibrated per-region map**, which `docs/amortized.md` already flags as the owner
  of the fine-grained local map (`local_coloc_map(...)` is documented as coarse-only today).
- The `ε = 0, shift = 0` byte-comparability check (D-10) is wanted explicitly as a **regression
  test**, not just a nice property: it proves the stage-6 refactor changed nothing for existing
  behaviour.
- The report should state the SC2/SC3 tension (D-13) in plain language rather than burying it — that
  widening the prior to demonstrate honest uncertainty is *precisely* what makes mis-registration
  undetectable by the density channel.

</specifics>

<deferred>
## Deferred Ideas

- **Summary extension for registration identifiability** — adding lag / cross-correlation-peak
  features so `dx`/`dy`/`ε` become estimable rather than marginalized. Rejected for Phase 11 (D-05)
  because it breaks summary reuse: the 128-dim read surface changes, invalidating the frozen net,
  the OOD Mahalanobis reference, `local_coloc_map`, and the artifact bundle, and it spends the
  fixed-summary-dimension constraint. Revisit only in a phase that already accepts a full re-cut.
- **A dedicated registration OOD channel** (cross-correlation-peak displacement features with their
  own Mahalanobis null and pre-registered threshold) — would break the SC2/SC3 tension outright
  instead of accepting it. Deferred because it needs fitted nulls in the shipped bundle, which
  contradicts the no-reship posture.
- **Staged marginalize-then-identify comparison** — run marginalize-only first, then a lag-feature
  net, and report how much identification buys over marginalization. Scientifically the strongest
  answer to "is the summary the bottleneck?"; deferred as roughly double the Phase-11 effort.
- **Off-axis optical centre** for the chromatic warp (currently pinned to image centre). A second
  chromatic parameter; out of scope for a "1-param warp".
- **Dedicated bead / registration-target acquisition** to calibrate `ε` against a real optical path.
  Strongest possible evidence but instrument-access-blocked; would turn Phase 11 into a wet-lab
  dependency.
- **API-level width advisory** (surfacing an interpretation hint when the `Δρ` posterior is unusually
  wide). Rejected for Phase 11 (D-16) because it touches the frozen public API and would need its
  own release.
- **Productionizing the Phase-11 research net** — reconciling the research input surface (D-03's
  extra conditioning input) with the shipped bundle, and re-cutting an artifact. Explicitly a
  future phase.

</deferred>

---

*Phase: 11-registration-and-chromatic-uncertainty-as-latent*
*Context gathered: 2026-07-25*
