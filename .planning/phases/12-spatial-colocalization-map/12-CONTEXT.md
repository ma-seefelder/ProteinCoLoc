# Phase 12: Spatial Colocalization Map (GP/CAR) - Context

**Gathered:** 2026-07-25
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace exchangeable patch pooling with a spatial lattice prior over the `correlation()` grid,
producing an **amortized per-region Δρ map with calibrated per-region uncertainty**. Fills the
`SpatialColocResult` extension point already sketched at `src/results.jl:181-191`.

**The constraint that defines this phase** — `local_map.jl:27-32`, citing `07-RESEARCH` Finding 2:

> "a finer patch grid lengthens the summary vector but **the NPE output stays a single GLOBAL ρ**,
> and windowing does not add a spatial prior, borrow strength between neighbours, or propagate
> per-region uncertainty."

Phase 12 is therefore *not* reachable by making the grid finer. It needs a different output head
**and** a simulator that produces spatially-varying ρ. ROADMAP names it the highest effort-risk phase
and the natural descope-to-v2.1 candidate — hence the pre-registered descope trigger (D-12/D-13).

**Two ROADMAP SC clauses are corrected by this discussion:**
- SC1 offers a "CNN/**DeepSet**" summary, but a DeepSet is *permutation-invariant over patches* —
  precisely the exchangeable pooling this phase's own goal says it replaces. Only the CNN half is
  coherent (D-02).
- SC3's "beats independent pooling in **coverage** on ≥1 real image" is undefined as written: there is
  no ground-truth ρ field on a real image, and the independent-pooling baseline has no uncertainty to
  compare (D-09/D-10).

**Out of scope:** the three-hypothesis evidence net (Phase 13), registration/chromatic latents
(Phase 11), the decision/abstention layer (Phase 14).

</domain>

<decisions>
## Implementation Decisions

### Spatial model + amortization

- **D-01: Hierarchical head — global ρ + low-rank per-region deviation field.** The global term stays
  exactly the quantity the shipped estimator reports, so continuity is preserved; calibration stays
  tractable (SBC on the global term plus a low-dimensional deviation summary, **not** 64 independent
  marginals); and shrinkage of deviations toward zero **is** the borrow-strength mechanism the CAR/GP
  prior is meant to supply.
  Rejected: a joint flow over all 64 lattice dims (the SBC apparatus already strained at 8 parameters
  with 5 vacuous); a per-region conditional over region index (yields no joint posterior — that *is*
  the independent-pooling baseline SC3 requires the model to beat).

- **D-02: Reshape the summary to G×G×2 and read it with a CNN. Summary CONTENT is unchanged.**
  **Correction recorded during discussion:** an earlier framing claimed Phase 12 must spend the
  fixed-summary constraint. It must not. The summary is `2·G² = 128` dims at G=8 (64 continuous +
  64 mask rows, **one entry per patch**), so per-region information is **already present** — flattening
  discards the spatial *inductive bias*, not the spatial *information*. SC1's CNN is an architecture
  choice, not an information necessity.
  Consequence: `patch_summary`, the OOD density Mahalanobis reference, `local_coloc_map`, the shipped
  artifact and `DATAGEN_HASH_SRC_FILES` all stay valid. Only the network's *view* of the same numbers
  changes. Comparability claims against the shipped MLP conditioner still need care, since the input
  contract differs.

- **D-03: Run the full CAR-vs-GP mini-spike, SC1 verbatim.** Chosen over a reduced spike specifically
  to keep a pre-registered criterion intact given Phase 7's amendment history.

- **D-04: The mini-spike runs with a deliberately HIGH-RANK deviation head; production head rank is
  chosen afterwards.** Rationale: a low-rank head bounds how *rough* a field it can represent, so a
  full-rank CAR field (sparse precision, lattice-native, capable of sharp neighbour jumps) would be
  systematically under-fit by it — a CAR loss would then be a **head artifact reported as a prior
  verdict**, the same misreading class that cost the Phase-7 SBC arm its credibility. Running the
  comparison through an expressive head makes it measure the prior. Cost accepted: the spike's head is
  not the production head.

### Spatially-varying simulator

- **D-05: Marginal-preserving copula construction.** Draw a correlated standard-Gaussian field from
  the CAR/GP lattice prior → push each cell through Φ → then the `MU_PRIOR` quantile function → apply
  `ghat` **elementwise** to obtain the ρ-field.
  **Why this and not a field on ρ directly:** `prior.jl:71-82` draws `μ* ~ Truncated(Cauchy(0,0.3),-1,1)`
  (verbatim the Turing μ-prior) and sets `ρ_true = ghat(μ*)`, so the induced per-patch correlation mean
  matches that prior *by construction* (SIM-02). Under the copula, **every region's marginal μ is
  exactly `MU_PRIOR`**, so SIM-02 holds per region and the CLAUDE.md constraint that π(θ) match the
  Turing μ/ν/σ/τ ranges is untouched — the lattice prior contributes **only spatial dependence**.
  Rejected: a field directly on ρ (breaks SIM-02 and ADVI comparability outright); global μ* plus a
  zero-mean deviation field (marginals become a convolution, so SIM-02 holds only approximately).

- **D-06: Apply the field by smooth interpolation to pixel resolution** (bilinear/bicubic upsample of
  the G×G field), not piecewise-constant per cell.
  **Rationale is an artifact risk, not aesthetics:** piecewise-constant cells create hard
  discontinuities aligned **exactly** with the patch grid the summary reads. The network could learn
  grid-aligned block-edge cues that do not exist in real images — passing SC2 in simulation and failing
  on the real images SC3 requires, with a failure mode that is hard to diagnose afterwards.

- **D-07: Per-region ground truth = the drawn G×G lattice values.** Interpolation (D-06) is a rendering
  detail with no bearing on scoring. SBC ranks the sampled parameter among posterior draws; making the
  scored quantity anything else inserts a deterministic transform between what was sampled and what is
  scored — precisely the structure that produced the ρ_true atom artifacts of named limit #2.
  Rejected: the cell mean of the interpolated field (a derived quantity); the induced per-cell patch
  correlation (a *measured, noisy* quantity, not a parameter — SBC ranks parameters).

- **D-08: The spatial correlation length is INFERRED as part of θ,** not fixed and not conditioned on.
  **Why not fixed:** if smoothness is fixed wrongly, the per-region uncertainty is wrong in a way that
  matched-simulation calibration **cannot** reveal — the simulations would be drawn at the same assumed
  length, making the test self-confirming. "Calibrated per-region uncertainty" *is* this phase's
  deliverable, so that assumption cannot be baked in.
  **Why not conditioned on (unlike Phase 11 D-03 / Phase 13 D-03):** registration uncertainty is
  something a user genuinely knows about their instrument; the spatial correlation length of biological
  colocalization is not — conditioning would ask for a number the caller has no basis to supply.
  Expect this column may prove unidentifiable from the patch-correlation summary (5/8 vacuous history);
  report it as such if so.

### Real-image coverage claim (SC3 amended)

- **D-09: SC3 becomes leave-region-out PREDICTIVE coverage.** Hold out lattice regions, predict them
  from their neighbours, and check the predictive interval against the held-out region's observed
  summary. **No ground truth is required anywhere**, and "beats independent pooling" becomes
  well-defined because independent pooling has no mechanism to predict a held-out region beyond the
  prior. It measures the borrowing-strength mechanism the phase claims to add.
  **SC3's wording is amended — pre-declare the amendment before results.** Predictive coverage is not
  parameter coverage; the claim is "predicts held-out regions honestly".

- **D-10: The "independent pooling" baseline is a MATCHED ABLATION** — the same Phase-12 network, same
  training, same summary, with the spatial prior neutralized (correlation length driven to zero).
  A coverage win is then attributable to the spatial prior and nothing else.
  Rejected as the *gate*: per-tile independent NPE reads — confounded by a different net reading
  smaller windows, so a win could stem from architecture, training or window size.
  **Note:** `LocalColocMap` (`local_map.jl:56-79`) cannot serve as a coverage baseline as-is — its own
  docstring states it "carries no posterior draws, no Bayes factor and no per-region uncertainty".

- **D-11: Score on BOTH `physical-primary` anchors** — the `coloc` anchor and the `segregated`
  (mitochondria) anchor. At G=8 that is ~128 scored regions across two genuinely different spatial
  regimes, and a spatial prior that improves coverage on a colocalized specimen but not a segregated
  one is an important finding, not a footnote.
  **The 30 CBS Red-Green rows do NOT qualify** — they are labelled `simulated-secondary`
  (Zinchuk & Grossenbacher-Zinchuk benchmark images are computer-generated), so they cannot support a
  "real image" claim. n=2 images is small, but coverage is scored **per region**, so the effective
  sample is regions, not images. State both facts.

### Descope trigger and v2.1 fallback

- **D-12: Two-stage pre-registered trigger.**
  **Stage 1 (cheap, before the training spend):** at the CAR-vs-GP mini-spike, require the better prior
  to show a stated coverage improvement over the neutralized-prior ablation on simulated data. If
  neither does, descope **before** the full training run.
  **Stage 2 (after training):** the per-region calibration gate must pass within one documented
  iteration.
  Rationale: a single post-training gate would place the descope decision *after* the phase's largest
  cost is sunk — the moment when an honest no-go is hardest to call. The mini-spike already computes
  the needed numbers, so the early gate is nearly free.

- **D-13: On descope, ship the ABLATION MODEL** — a per-region Δρ map **with** per-region uncertainty,
  deferring only the spatial borrowing to v2.1. It is built regardless as the D-10 gate baseline, so it
  costs nothing extra, and it is strictly more than `LocalColocMap` offers today (which has no
  uncertainty at all). The descope becomes a reduced release rather than a write-off.
  **Caveat:** if the Stage-1 trigger fires, the ablation exists only at spike scale — this fallback is
  fully available only at the Stage-2 gate.

### Carried forward (not re-asked)

- **Research-lane posture.** Not explicitly discussed for Phase 12, but Phase 11 D-01 and Phase 13 D-01
  both chose it, and the 2026-07-24 GO locked "no further training, no gate iteration". **Planners
  should assume research lane** — train in `spike/` on a fresh DEV seed asserted disjoint from
  `PROD_SEED_V2`, `VAL_MASTER_SEED` (`0x5BC0FFEE`), `NPE_MASTER_SEED` (`0xC0FFEE`), `RATIO_PAIR_SEED`
  (`0x00000000004A7107`), `VAL_FIX_SEED` (`0xF1F7ED`), `CORPUS_MASTER_SEED` (`0x0000000000c05eed`) and
  prior dev seeds; shipped artifact and GO untouched. **Flagged as an assumption, not a locked answer**
  — confirm before executing if it matters.
- **Pre-registration discipline** (Phase 11 D-04): lock thresholds in a Phase-12 consts file with the
  fresh DEV seed before any run; one documented iteration allowed. D-12 already assumes this.
- **Phase 7 D-02:** new result variants slot in as new `AbstractColocResult` subtypes, never as
  bolted-on fields.

### Claude's Discretion

- Nothing was delegated to Claude's discretion in this discussion. Every gray area was decided by the
  user.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### The constraint and the extension point
- `src/amortized/local_map.jl:19-40` — the honest Phase-7 answer to local localisation, and **:27-32
  the 07-RESEARCH Finding 2 statement** that a finer grid does not produce a per-region posterior.
  **:56-79 `LocalColocMap`** — explicitly no posterior draws, no BF, **no per-region uncertainty**
  (why D-10's baseline had to be constructed).
- `src/results.jl:180-191` — the `SpatialColocResult` sketch (`region_delta_rho`, `region_sd`, `ood`,
  `calibration`, `meta`) with `delta_rho_map` / `uncertainty_map` accessors. Phase 12 fills this.
  `:161-167` — Phase 7 D-02's rule that variants are new subtypes.

### Summary and architecture (D-01, D-02)
- `src/amortized/summary.jl` — `patch_summary`, `summary_dim(grid)`, `_summary_row_partition`
  (the single grid-coupling site), `encode_d01`. **The 2·G² = 64+64 layout D-02 reshapes.**
- `src/amortized/architecture.jl:65, 187` — `theta_prior_bounds()` provenance and the
  `NormalisingFlow` over the D-dim (default 7) θ. Adding the deviation field and the correlation
  length (D-08) changes θ arity, so both move.

### Simulator and the SIM-02 contract (D-05 … D-08)
- `spike/simulator/prior.jl:21-32, 40, 71-82` — the SIM-02 consistency statement, `MU_PRIOR`,
  `ρ_true = ghat(μ*)`, and the explicit warning **"Do NOT assert ρ_true == μ"** that D-07 turns on.
  `:42-45` — the pre-declared W1 tolerance pattern for a consistency claim.
- `spike/simulator/ghat.jl:42-67` — `GHAT_MU_KNOTS` / `GHAT_RHO_KNOTS`, `GHAT_MU_MIN = -0.67976`,
  `GHAT_MU_MAX = 0.847149`, and the clamping that creates the prior atoms. D-05 applies `ghat`
  elementwise, so the atoms become **per-region** — plan for that.
- `spike/simulator/forward.jl` — the stage pipeline the ρ-field enters. **Phase 11 reworks stage 6
  into a single composed affine warp; D-06's field application stacks on top of that.**
- `src/bayes.jl` ~274-291 — the Turing μ-prior π(θ) must stay consistent with (CLAUDE.md constraint).

### Ground-truth data (D-11)
- `corpus/manifest.csv` — `truth_label` / `tier` columns. **Only two rows are `physical-primary`**
  (`positive|coloc` and `negative|segregated|mitochondria`); the 30 CBS Red-Green rows are
  `simulated-secondary`. `CORPUS_MASTER_SEED = 0x0000000000c05eed` in the header.
- `corpus/fetch.jl`, `corpus/load.jl` — the read path.

### Honesty contract and prior findings
- `docs/amortized.md` §Windowed sub-tile local map (~95-104) — states the coarse map is **not** a
  calibrated per-region posterior and that Phase 12 owns that. §Named limits — limit #2 (ρ_true atoms)
  and limit #4 (vacuous nuisance columns) both bear on D-07 and D-08.
- Skill `spike-findings-proteincoloc` → `references/sbc-calibration.md` — randomized-rank recipe,
  M=2000 over-power diagnosis, spike 012 (capacity redistributes), spike 013 (μ-truncation costs prior
  consistency).
- `.planning/phases/07-productionization-conditional-on-go/07-CONTEXT.md` — D-02 result hierarchy,
  D-04 capped grid family, D-05 per-grid re-pre-registered ship gate.
- `spike/validation/consts.jl:72-79` — the pre-registered seed constants and the file pattern a
  Phase-12 consts file should mirror.

### Upstream and sibling phases
- `.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-CONTEXT.md` — **ROADMAP
  makes Phase 12 depend on 11** ("the spatial map trains on the registration-aware θ"); Phase 11 D-10's
  composed affine warp is the stage-6 baseline D-06 builds on, and D-05 there mirrors the
  summary-unchanged posture D-02 here preserves.
- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-CONTEXT.md` — sibling research-lane
  phase; its D-16 semi-synthetic construction is the pattern the deferred parameter-coverage option
  here would follow.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`patch_summary` / `correlation()` / `patch()`** — untouched by D-02; only the network's view of
  their output changes, so the whole upstream summary path is reused verbatim.
- **`ghat`** — applied elementwise in D-05; no change to the frozen knots or the inverse itself.
- **`LocalColocMap`** — not a coverage baseline (no uncertainty), but its sentinel/OOD discipline
  (`LOCAL_MAP_SENTINEL`, forced OOD flag on unscorable tiles, T-7-04) is the pattern the spatial map's
  degenerate-region handling should follow.
- **`spike/validation/consts.jl`** — the pre-registration file pattern (named seed constants with
  explicit disjointness assertions) for D-12's two thresholds.
- **Existing SBC/coverage harness** (`spike/validation/harness.jl`) — the per-region coverage runs
  should extend it rather than build a parallel harness.

### Established Patterns
- **SIM-02 consistency is proved, not assumed.** `prior.jl:42-45` pre-declares a W1 tolerance for the
  induced-μ claim *before* the calibration run. D-05's per-region marginal claim should be verified the
  same way, per region.
- **Vacuous-column reporting (F3).** `shrinkage = post_sd / prior_sd` with a boolean `vacuous` flag is
  a reporting-only diagnostic. The correlation length (D-08) and the deviation field components must be
  reported through it.
- **Grid coupling lives in one place** — `_summary_row_partition` in `summary.jl`. D-02's reshape
  belongs there or immediately adjacent, not scattered.
- **Pre-register before running.** Every gate in this project locked thresholds and seeds first.
- **Named limits are the honesty mechanism** — the SC1 DeepSet correction, the SC3 amendment, and the
  n=2 real-image base all belong on that list.

### Integration Points
- `src/results.jl:180-191` — the `SpatialColocResult` slot.
- `spike/simulator/forward.jl` stage 6 — Phase 11's composed affine warp; D-06's spatially-varying
  mixing stacks on it. **Sequencing matters: Phase 11 lands first.**
- `theta_prior_bounds()` ↔ `architecture.jl:65,187` — θ arity changes with the deviation field and the
  correlation length.
- The ablation model (D-10) is simultaneously the SC3 gate baseline and the D-13 descope deliverable —
  it should be built as a first-class artifact, not a throwaway comparison arm.

</code_context>

<specifics>
## Specific Ideas

- The grid-alignment artifact (D-06) is the failure mode most worth guarding against explicitly: a map
  that scores well in simulation because it detected block boundaries would pass SC2 and fail SC3, and
  the diagnosis would be expensive after the fact. Worth an explicit check that the trained model's
  performance does not degrade when the field grid is offset from the summary grid.
- The self-confirming-calibration argument behind D-08 should be stated in the manuscript, not just the
  code. "We inferred the correlation length rather than fixing it, because a fixed smoothness makes the
  calibration claim circular" is a methodological point a reviewer will credit.
- D-13 turns the descope from a write-off into a release. Worth framing that way in the Go/No-Go memo
  if the trigger fires — the phase would still have delivered per-region uncertainty, which the package
  has never had.

</specifics>

<deferred>
## Deferred Ideas

- **A new CNN summary over the raw `correlation()` grid**, unconstrained by the current 128-dim
  encoding (SC1's literal reading). Rejected for D-02 because it changes summary *content* and would
  invalidate the shipped read surface, OOD reference, `local_coloc_map` and artifact — the cost Phase 11
  D-05 declined to pay. Revisit only in a phase that already accepts a full re-cut.
- **A DeepSet summary.** Named in SC1 but incoherent with the phase goal (permutation-invariance over
  patches *is* the exchangeable pooling being replaced). Recorded so nobody re-proposes it.
- **Joint flow over all 64 lattice dimensions.** Fully expressive and able to represent sharp
  boundaries; deferred because the calibration burden (64 marginals) exceeds what this project's SBC
  apparatus has demonstrated at 8.
- **Per-region conditional head over region index.** Cheap and grid-agnostic; deferred because it
  produces no joint posterior and would fail SC3 by construction.
- **Per-tile independent NPE reads as a product-level baseline.** Rejected as the *gate* (confounded),
  but reporting it alongside the ablation would answer "is this better than what we ship?" as distinct
  from "does the spatial prior help?".
- **Semi-synthetic imposed-field series on real images** for *parameter* coverage (as opposed to D-09's
  predictive coverage), following Phase 13 D-16's construction pattern. The cleanest way to strengthen
  SC3 if predictive coverage proves a weak signal.
- **CBS Red-Green benchmark as a graded supplementary series.** 30 images spanning coloc degree 0.0-0.9,
  but `simulated-secondary` — usable only as clearly-labelled supplementary evidence, never for the
  "real image" claim.
- **Matching a head to each prior in the mini-spike** (full-rank for CAR, low-rank for GP). Rejected for
  D-04 because it confounds in the opposite direction, but it is the fairer comparison if the high-rank
  spike proves inconclusive.
- **Productionizing the spatial map** — reconciling the CNN input contract with the shipped bundle and
  re-cutting an artifact. Explicitly a future phase.

</deferred>

---

*Phase: 12-spatial-colocalization-map*
*Context gathered: 2026-07-25*
