# Phase 16: External Validation and Manuscript Assembly - Context

**Gathered:** 2026-08-04 (discuss --analyze --all)
**Status:** Ready for planning

<domain>
## Phase Boundary

Close v2.0 with a blind external evaluation against the physical corpus (Phase 8) and the
comparator harness (Phase 9), assemble every manuscript figure into the Phase-10 Typst skeleton,
and ship a one-command seeded reproduction plus a Zenodo/DOI release.

**Seven scouting findings, verified on disk. Five change what can be built.**

1. **The physical anchors have no hashes.** Both rows carry the literal `PENDING-FETCH` sentinel
   (`corpus/anchor_rows.jl:85`, `is_pending_hash`). The positive anchor is a ~6.3 GB archive whose
   bootstrap fetch was **deliberately never run** — an explicit human decision recorded in Phase 8's
   completion note, which also states the fetch must happen **before Phase 16 opens the sealed
   holdout**. The sealed holdout cannot honestly be unsealed against data whose provenance digest is
   a sentinel string. Resolved by D-04.

2. **The positive anchor is 100 nm TetraSpeck beads, not a tandem-FP construct.** Phase 8 searched
   the public archives and found no open-licensed, non-environment-quenched tandem-FP dataset (every
   deposited set is a quenched autophagy reporter, disqualified by its Pitfall 2). Beads are
   sub-diffraction point sources; the v2.0 simulator generates correlated density fields. **The OOD
   flag may legitimately fire on the project's only physical ground truth.** Resolved by D-05.

3. **The two anchors are cross-study and cross-archive** (`ANCHORS_MATCHED=false`). Imaging-condition
   confounds between positive (Zenodo `10.5281/zenodo.5509861`) and negative (BioImage Archive
   `S-BIAD1047`) are **not controlled**; any separation measured could be a batch effect. Phase 8
   already ruled that Phase 16 must report this. Resolved by D-06.

4. **Two specced figure slots have no positive result behind them.** F5 (registration-UQ widening) —
   Phase 11 closed negative-but-useful: registration ≤3 px is not inferable from the 8×8 summary at
   any coloc level, and Δρ RMSE is flat in λ (ratio 1.0003) (`11-CLOSURE.md`). F6 (spatial per-region
   Δρ map) — Phase 12's goal was **not achieved**; after the DCT-permutation correction all three arms
   fail the band and the winning ablation was never tested (`12-VERIFICATION.md`, `gaps_found`).
   Resolved by D-07.

5. **F7's three-way BF and the Phase-14 decision layer live in `spike/`, not `src/`.** The shipped
   v2.0 public API is amortized-only on the `amended_v2/grid_8` reference grid. A figure drawn from
   research-lane code is not a figure of the released tool. Resolved by D-08.

6. **The C-04 pool-overlap hazard fires here by name.** `.planning/STATE.md` flags Phases 14/15/16
   explicitly: index-keyed generation means *"a later phase generates 'a fresh evaluation pool' at
   `arm = P12_CHOSEN_PRIOR` — under C-04 its first `min(n, 50000)` samples ARE the production training
   set. Any coverage, SBC or BF number so obtained is measured on training data and is not a
   calibration result."* Resolved structurally by D-11.

7. **The deposit's licences do not all compose.** The repo is **AGPL-3.0**; `artifacts/` is **82 MB**
   (`npe_8.jld2`, `ratio_8.jld2`, `ood_nulls_8.jld2`, `gate_report_8.jld2`); `git ls-files corpus/data`
   is **empty** — third-party imagery has never been tracked. CBS is **CC-BY-NC-SA-4.0**
   (non-commercial *and* share-alike); both physical anchors are **CC-BY-4.0**. Bundling CBS would
   propagate NC-SA to the collection. Resolved by D-10.

**Out of scope:** retraining or retuning any net; promoting `spike/p13`/`spike/p14` into `src/`;
adding OOD detector channels; re-opening Phase 8's anchor search; the Phase-15 envelope/CI work
itself (Phase 15 owns it — Phase 16 consumes its output as F8).

</domain>

<decisions>
## Implementation Decisions

### Blind-evaluation protocol and verdict rules

- **D-01: Pre-registration + single authorized shot.** The evaluation script **and** a frozen
  expected-outcome table are committed **before** the sealed holdout is unsealed. One authorized run.
  Any re-run must be disclosed in the manuscript with its reason. This mirrors the pre-registration
  discipline already used in this project (`test/gate/gate_consts_8_v2.jl` literal-value assertions;
  Phase 13's Tier-1/Tier-2 const tests; Phase 15's pre-declared consts). A reviewer must be able to
  diff the pre-registration commit against the result commit.

- **D-02: The gate is ORDINAL SEPARATION + HONEST UNCERTAINTY; point accuracy and the decision-layer
  outcome are computed and reported but do NOT gate.**
  - **Gating:** the positive anchor's posterior must sit above the negative's with non-overlapping
    declared intervals; a well-founded ABSTAIN counts as honest rather than as failure.
  - **Non-gating, reported:** point accuracy against the asserted ρ=1.0 / ρ=0.0 truths, and the
    Phase-14 decision-layer outcome (`coloc` / `not` / `ABSTAIN`).
  - **Why ordinal is the gate:** the anchors sit at the *extremes*, which is exactly where the
    ρ_true prior-atom artifact lives (a named limit in `docs/amortized.md`). A point-accuracy miss at
    an extreme may measure the prior rather than the tool. The anchors are additionally
    cross-study confounded (finding 3), so a point bar would not be a fair test. Ordinal separation
    is what this corpus can genuinely falsify.
  - Gating on the Phase-14 layer was explicitly rejected: it is unshipped (`spike/p14`) and its own
    SC3-d came in at **−0.129 against a floor of 0.50**.

- **D-03: On gate failure — one bounded remediation attempt, disclosed, then re-run.**
  A diagnosed cause may be fixed once (harness-side *or* model-side) and the evaluation re-run.
  **The remediation and the re-run are both disclosed in the manuscript.**
  **Consequence that must be stated, not inferred:** a post-remediation re-run **is no longer blind**.
  The manuscript reports the original blind result *and* the post-remediation result, each labelled
  for what it is. Exactly one remediation is authorized; a second failure is published as-is.

### Physical-anchor readiness and the OOD question

- **D-04: Full fetch, real digests, then a pre-registered seeded subset for evaluation.**
  A dedicated gated task runs `bootstrap_anchor_hashes()` to download both anchors (~6.3 GB
  positive), records real SHA-256 over the **full archives**, and replaces the `PENDING-FETCH`
  sentinel in `corpus/manifest.csv`. Evaluation then runs on a subset selected by a seed that is
  committed in the **same pre-registration commit as D-01**, so the subset choice cannot become a
  tuning knob. Full-archive digests mean provenance covers everything, not just what was scored.

- **D-05: Two-track OOD reporting with an explicit anti-vacuity clause.**
  - **Track 1:** the OOD verdict is a **first-class reported result** — both anchors scored against
    the calibrated OOD null. A flag that fires on beads but not on a simulator-matched control is
    *positive evidence the flag works* (claim C4), not a failure.
  - **Track 2:** the coloc posteriors are reported **regardless**, explicitly labelled "computed
    under a raised OOD flag", so the ordinal test in D-02 retains content.
  - **Anti-vacuity clause (pre-registered):** *both anchors flagged with no separation demonstrated*
    is declared **INCONCLUSIVE**, never a pass. Without this, D-02's abstain-is-honest clause would
    silently become a free pass.
  - Suppressing the flag for the anchor run was rejected outright — it would disable the paper's
    named differentiator on the one dataset where it matters.

- **D-06: Both controls for the cross-study confound.**
  - Run the Phase-9 comparator battery (Costes-p, Manders M1/M2, Pearson, Spearman) on the
    **identical** anchor images. If the classics separate the anchors the same way, the separation is
    a property of the images rather than of v2.0 — that is the discriminating evidence. The Phase-9
    harness is input-source-agnostic by construction (`09-03`), so this is a wiring task.
  - Use the **CBS 0–100% series as a condition-matched graded ladder** (single source, matched
    imaging), reported strictly at the `simulated-secondary` tier and **never conflated** with
    `physical-primary` — the manifest's D-06 tier guard already enforces the split.
  - The confound is stated explicitly in the manuscript regardless of outcome.

### Manuscript figure assembly

- **D-07: F5 becomes a rendered NEGATIVE panel; F6 becomes text-only with its spec retained.**
  Evidentiary strength drives the format. F5 has a clean quantified null (Δρ RMSE flat in λ, ratio
  1.0003; registration not inferable at any coloc level) that supports a real methods recommendation
  — *calibrate registration externally with fiducial beads* — which the TetraSpeck positive anchor
  directly exemplifies. F6's result is a corrected sign, three failing arms and an untested
  exemption; a panel would imply a cleaner finding than exists, so it goes in Discussion.
  **Numbering is unchanged and F7 keeps its number** — renumbering would break the `fig:` linkage in
  `manuscript/claims.typ` that the Phase-10 machinery exists to enforce.

- **D-08: Add F8 and F9; every panel carries a provenance lane tag.**
  - **F8 — Phase-15 operating envelope:** coverage vs. the six typed nuisance axes, with in-prior
    extrapolation held visually distinct from out-of-model misspecification, plus the
    "OOD fires before coverage breaks" result.
  - **F9 — Phase-14 decision/abstention layer** at controlled Bayesian FDR, **including the missed
    SC3-d (−0.129 vs a 0.50 floor)**.
  - **Lane tag on every panel** (one field in the figure YAML, rendered in the caption):
    `shipped (src/)` vs `research lane (spike/)`. This makes the shipped-vs-prototype boundary a
    property of each figure rather than a footnote a reader can miss — F7 and F9 are research-lane.
  - The claim table (`manuscript/claims.typ`) status column is updated in the same pass: C3/C4 move
    off `reserved slot`, and any claim whose supporting phase closed negative is restated honestly.

### Reproducibility

- **D-09: `run_v2.jl` ships THREE tiers, with measured wall-clock per tier recorded in the repo.**
  - `--smoke` — minutes, tiny N, wiring assertions only. **Wired into the Phase-15 CI gate**, so the
    entry point is continuously proven rather than asserted once.
  - `--figures` — restore pinned artifacts, re-run inference, render every panel. This is the tier a
    reviewer actually runs.
  - `--full` — regenerate pools and retrain from the fixed seed. Counter-based Philox makes this
    **byte-exact**, not merely statistically similar.
  - Wall-clock for each tier is **measured and recorded**, not estimated — untimed tiers are
    marketing.

- **D-10: Zenodo deposit = code + manifests + the 82 MB derived artifacts + figure sources +
  `run_v2.jl`. All third-party imagery by accession + SHA-256.**
  Nothing third-party is redistributed. This keeps the AGPL-3.0 deposit free of CBS's CC-BY-NC-SA-4.0
  non-commercial/share-alike terms, and matches how the repo already works (`corpus/data` has never
  been tracked). A digest plus `fetch_verified` is **stronger** provenance than a re-upload — a
  re-upload can drift, a digest cannot.

- **D-11: Simulator-side draws come from the shipped nulls first, the validation stream second, and
  NEVER from a fresh pool at the production arm.**
  - Score against the shipped `artifacts/amended_v2/grid_8/ood_nulls_8.jld2` wherever it suffices —
    zero new draws, so C-04 cannot fire, and it measures the *released* tool rather than a bespoke
    re-fit.
  - Where the bundle is insufficient, draw fresh on the **validation stream via `test/gate/harness.jl`**
    (different salt, disjointness asserted at load time) — the correction `.planning/STATE.md`
    itself prescribes.
  - Generating at `arm = P12_CHOSEN_PRIOR` is **prohibited by an executable assertion in the
    evaluation script**, not by a comment. This project has five recorded instances of guards
    satisfied by a comment rather than by code (ROADMAP Backlog 999.1) — the assertion must fail red
    when violated.

### Claude's Discretion

- Exact subset size for D-04's seeded anchor subset (must be pre-registered; size itself is a
  planning call).
- Panel layout, palette and Typst plumbing for F8/F9 — clone the F1 pattern
  (`manuscript/figures/f1_speedup.{yaml,typ}` + `lib/figstyle.typ`) per the Phase-10 D-14/D-16
  convention.
- Whether F8 is one multi-panel figure or splits into F8a/F8b, depending on what Phase 15 returns.
- CLI surface details of `run_v2.jl` beyond the three tier flags.

</decisions>

<specifics>
## Specific Ideas

- The bead/registration link is worth making explicit in the writing: Phase 11's negative concludes
  *"calibrate registration externally with fiducial beads"*, and the positive physical anchor **is**
  100 nm fiducial beads. F5's negative panel and the external-validation section should reference
  each other rather than sit as unrelated results.
- An OOD flag that fires on beads is a **better** paper result than a coloc number would be, provided
  the anti-vacuity clause (D-05) keeps it from being a free pass. Write it that way, not
  apologetically.
- The manuscript must never quote a track-2 posterior (D-05) without its "under a raised OOD flag"
  label attached.

</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope and prior rulings
- `.planning/ROADMAP.md` §"Phase 16: External Validation and Manuscript Assembly" — goal, deps, SC1–SC3
- `.planning/ROADMAP.md` §"Phase 8 … Completion note (2026-07-21)" — the D-01 substitution, the unmet
  D-02 preference (`ANCHORS_MATCHED=false`), and the **hashes-pending** precondition binding this phase
- `.planning/ROADMAP.md` §"Backlog" 999.1 — the parked `spike/p14` OOD-wiring traps and the
  comment-satisfied guard pattern D-11 must not repeat
- `.planning/STATE.md` §"⚠️ FOR PHASES 14/15/16 — a pool-overlap hazard" — the C-04 statement in full
- `.planning/CONVENTIONS.md` **C-04** — index-keyed generation and why a "separate scoring pool" is a
  superset of the training pool
- `.planning/phases/07-productionization-conditional-on-go/07-GO-NO-GO-UPDATE.md` — the GO-with-named-limits
  decision the manuscript's claims must not exceed
- `docs/amortized.md` — the shipped named limits (nuisance marginal drift, twice-amended gate, atom
  handling) that the paper is required to state

### Corpus (Phase 8) — the sealed holdout
- `.planning/phases/08-external-physical-ground-truth-corpus/08-CONTEXT.md` — **D-09** sealed-holdout
  split policy, **D-06** physical/simulated tier guard, **D-03** anchor acceptance predicate
- `corpus/anchor_rows.jl` — `PENDING_FETCH` sentinel, `is_pending_hash`, `anchor_rows(...)`,
  `bootstrap_anchor_hashes()` (D-04 acts here)
- `corpus/manifest.csv` — the pre-registered manifest; anchor rows carry `split=sealed_holdout`
- `corpus/config.jl`, `corpus/fetch.jl`, `corpus/hash.jl`, `corpus/manifest.jl`, `corpus/load.jl`
- `.planning/phases/08-external-physical-ground-truth-corpus/08-VERIFICATION.md`

### Comparators (Phase 9)
- `.planning/phases/09-cross-method-comparator-harness/09-CONTEXT.md` — D-13/D-14 consts, seeded
  Costes block-scramble, input-source-agnostic shared-input builder (D-06 wires the anchors in here)
- `.planning/phases/09-cross-method-comparator-harness/09-05-SUMMARY.md` — comparison table,
  divergence/traffic-light, content-addressed CSV/JLD2 output

### Manuscript (Phase 10)
- `manuscript/main.typ`, `manuscript/claims.typ` — the claim spine; `fig:` linkage D-07 preserves and
  the `status` column D-08 updates
- `manuscript/figures/specs.typ` — F1–F7 specs; F8/F9 are added here
- `manuscript/figures/f1_speedup.{yaml,typ}` + `manuscript/lib/figstyle.typ` — the rendered-figure
  pattern every new panel clones (Phase-10 D-14/D-16)
- `.planning/phases/10-manuscript-skeleton-and-related-work-positioning/10-CONTEXT.md` — D-01..D-12
- `.planning/phases/10-manuscript-skeleton-and-related-work-positioning/10-04-PLAN.md` — figure
  specifications enumerating per-phase panels linked to claim rows

### Upstream results the figures report
- `.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-CLOSURE.md` — F5's negative
  (flat RMSE in λ, ratio 1.0003; SC1g mis-specified, not failed)
- `.planning/phases/12-spatial-colocalization-map/12-VERIFICATION.md` — F6's `gaps_found`; goal not
  achieved; DCT-permutation correction; untested winning ablation
- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-VERIFICATION.md` — F7 (goal achieved;
  research lane, `spike/`)
- `.planning/phases/14-decision-and-abstention-layer/14-REPORT.md` and `14-VERIFICATION.md` — F9,
  including the **SC3-d miss at −0.129 vs a 0.50 floor** D-08 requires to be shown
- `.planning/phases/14-decision-and-abstention-layer/14-REVIEW.md` — CR-01/CR-02/WR-01, the traps
  parked as Backlog 999.1
- `.planning/phases/15-calibration-operating-envelope-and-ci-gate/15-CONTEXT.md` — F8's source; the
  six typed axes (D-02), in-prior vs out-of-model typing (D-03), and the two-tier CI design (D-08)
  that `run_v2.jl --smoke` hooks into

### Shipped artifacts and gate machinery
- `artifacts/amended_v2/grid_8/` — `npe_8.jld2`, `ratio_8.jld2`, `ood_nulls_8.jld2`,
  `gate_report_8.jld2` (82 MB total; D-10 deposits these, D-11 scores against `ood_nulls_8`)
- `test/gate/harness.jl` — the seeded draw→simulate→infer path; D-11's validation-stream route
- `test/gate/misspec.jl` — `OOD_FAMILIES = (texture, noise, optics, background)`
- `test/gate/run_gate.jl` — the `PROD_SEED_V2[8]` single-run binding (already spent; not for this phase)
- `LICENSE` — AGPL-3.0, the constraint behind D-10

### Packaged findings
- `Skill("spike-findings-proteincoloc")` — BF/SBC calibration diagnosis and the src/ refactor roadmap;
  read before touching the BF or calibration narrative

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`corpus/` (Phase 8, 217/217 offline gate green):** `fetch_verified` (atomic, timeout,
  skip-vs-abort asymmetry), SHA-256 content hashing, manifest schema + validation, tier guard,
  sealed-holdout guard, read-only `MultiChannelImage` conversion. D-04 and D-10 are wiring on top of
  this, not new machinery.
- **Phase-9 comparator harness:** already input-source-agnostic (`09-03`) and seeded — D-06 feeds
  anchor images into the existing shared-input builder rather than writing a parallel path.
- **Phase-10 figure pattern:** `f1_speedup.{yaml,typ}` + `lib/figstyle.typ` + `claims.typ`
  data-plus-render split. Every value in a panel flows from the YAML; F5/F8/F9 clone this.
- **`test/gate/harness.jl` / `misspec.jl` / `sbc.jl`:** seeded draw→simulate→infer, OOD families,
  rank/ECE/coverage — D-11's safe draw route and F8's source.
- **Content-addressing discipline throughout** (Phase-3 cache → Phase-8 manifest → Phase-9 CSV/JLD2):
  SHA-256 over data-defining bytes, atomic `.tmp → verify → mv`. `run_v2.jl` artifact restore (D-09
  `--figures`) should reuse it rather than invent a download step.

### Established Patterns
- **Pre-registration by literal-value assertion** (`gate_consts_8_v2.jl`, `test_p13_consts.jl` 137/137,
  `test_p14_decoupling.jl`) — the mechanism D-01 and D-05's anti-vacuity clause must use.
- **Structural, not conventional, separation** — holdouts enforced by disjoint indexing / distinct
  namespaces, never by comments. Binds D-11.
- **Negative-but-useful closure** — Phases 11 and 12 are closed negatives on the record. D-07's
  handling of F5/F6 continues that posture rather than quietly dropping them.
- **`src/` is not edited outside Phase 7** — Phases 13/14 built in `spike/`; Phase 15 builds in
  `test/gate/`. Phase 16 writes to `corpus/`, `manuscript/`, and a top-level `run_v2.jl`; it has no
  mandate to touch `src/`.

### Integration Points
- `corpus/manifest.csv` — D-04 rewrites the two anchor `sha256`/`bytes` fields in place
- `corpus/load.jl` sealed-holdout accessor — the single unsealing point, gated by D-01
- Phase-9 shared-input builder — D-06's anchor wiring
- `manuscript/figures/specs.typ` + `claims.typ` — D-07/D-08 edits
- Phase-15 CI workflow (in flight) — D-09's `--smoke` hook

</code_context>

<deferred>
## Deferred Ideas

- **Bundling the CC-BY anchor imagery into the Zenodo deposit** (~6.3 GB) as insurance against
  upstream link rot — rejected for the v2.0 deposit under D-10; revisit if either archive shows
  instability.
- **A condition-matched physical anchor pair** (positive and negative from one study on one
  microscope) that would remove the D-06 confound rather than merely attribute it — Phase 8
  established none exists publicly; this is a wet-lab acquisition, i.e. its own project.
- **Fixing Backlog 999.1** (`spike/p14` OOD-wiring traps, comment-satisfied guard) — stays parked;
  Phase 16 must not repeat the pattern (D-11) but does not fix the existing instances.
- **Promoting `spike/p13`/`spike/p14` into `src/`** — out of scope; D-08's lane tags exist precisely
  because this has *not* happened.
- **Option B from the Go/No-Go** (BF §6 gate integration + summary redesign) — deferred at the GO
  decision and still deferred.

</deferred>

---

*Phase: 16-external-validation-and-manuscript-assembly*
*Context gathered: 2026-08-04*
