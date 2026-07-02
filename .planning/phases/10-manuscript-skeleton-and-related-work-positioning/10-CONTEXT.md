# Phase 10: Manuscript Skeleton and Related-Work Positioning - Context

**Gathered:** 2026-07-02 (--auto mode; recommended defaults selected autonomously)
**Status:** Ready for planning

<domain>
## Phase Boundary

Deliver a **compiling manuscript skeleton** for the ProteinCoLoc v2.0 paper with related-work
positioning drafted early — the explicit delta versus **Tapqir / Costes / Manders** — so the
downstream experiment phases (5–9, 11–15) are shaped by the claims they must support.

**In scope (from ROADMAP.md Phase 10 success criteria):**
1. A Typst skeleton that **compiles** with section scaffolding and a claim table.
2. Related-work prose drafting the **Tapqir differentiation** along four axes: amortization +
   SBC/calibration + registration-UQ + spatial colocalization map.
3. **Figure specifications** enumerating the panels each later phase must deliver.

**Explicitly NOT in scope for this phase:**
- Writing final Results/Methods prose or real numbers (those land as later phases produce them;
  this phase only scaffolds and reserves slots).
- Generating actual figures (Phase 16 assembles them; this phase specifies them).
- Running any experiment or touching `spike/` model code or `src/`.
- Choosing/locking a final target journal (skeleton stays venue-neutral; candidates noted only).
</domain>

<decisions>
## Implementation Decisions

### Manuscript tooling & format
- **D-01:** Author in **Typst** (already the ROADMAP-mandated format; Typst **0.15.0** is
  installed at `~/AppData/Local/.../typst`). The skeleton must actually compile with the
  installed binary — a broken skeleton is not a deliverable.
- **D-02:** Keep the skeleton **venue-neutral** — a generic single-column Typst `article`
  scaffold, not a specific journal template. Record candidate venues (e.g., Bioinformatics,
  Nature Methods, eLife-style methods paper) as a comment/note only; do not lock one. Rationale:
  target journal is undecided and Phase 16 assembles the final; premature template lock-in adds
  rework risk.

### Manuscript location & decoupling
- **D-03:** Create a **new top-level `manuscript/` directory** at repo root (sibling to `spike/`,
  `src/`, `test/`). Do **not** place manuscript sources under `spike/` or `src/`. This honors the
  project decoupling constraint: the spike keeps `src/`, `bayes.jl`, `colocalization.jl` and both
  manuscript pipelines provably untouched. The v2.0 manuscript is a genuinely new artifact and
  belongs in its own tree.
- **D-04:** Suggested internal layout: `manuscript/main.typ` (entry), `manuscript/sections/*.typ`
  (per-section includes), `manuscript/refs.bib`, `manuscript/claims.typ` (or a data file feeding
  the claim table), `manuscript/figures/` (specs + later PDFs). Planner may refine, but section
  files must be modular so later phases append prose without merge collisions.

### Section scaffolding
- **D-05:** Scaffold the standard methods-paper sections as empty-but-labelled stubs with a
  one-line intent comment each: Abstract, Introduction, Related Work, Methods (forward simulator,
  summary statistic, NPE/NRE, SBC, OOD), Results (placeholder subsections mirroring the claim
  table), Discussion, Data/Code Availability, References. Each stub carries a `// TODO(phase N)`
  marker pointing at the phase that will fill it.

### Claim table structure
- **D-06:** The claim table is the **spine** of the phase — a structured table mapping each
  **scientific claim → differentiator axis → supporting phase(s) → figure/experiment ID →
  status**. It must be legible both as a rendered Typst table and as a plain source block a
  human/agent can diff. Seed it with v2.0's headline claims (amortized ms-scale inference,
  >100× speedup at comparable RMSE, SBC calibration-under-the-simulator, honest OOD flag,
  registration-UQ, spatial map, three-hypothesis BF) each linked to its owning phase.
- **D-07:** Every claim row must name the artifact that will substantiate it (a figure spec ID
  and/or a phase deliverable), so "experiments are shaped by the claims they must support" is
  enforced by construction, not left to prose.

### Related-work delta (Tapqir / Costes / Manders)
- **D-08:** Draft related work as **comparison matrix + prose**. The matrix has methods as rows
  (ProteinCoLoc v2.0, Tapqir, Costes randomization, Manders M1/M2, Pearson/Spearman, original
  ProteinCoLoc v1 ADVI) and capability axes as columns: amortized inference, calibrated
  uncertainty / SBC, Bayes-factor / model comparison, registration uncertainty as latent,
  spatial (per-region) map, misspecification/OOD flag, speed.
- **D-09:** The **Tapqir differentiation** paragraph is mandatory and must explicitly cover all
  four axes from the success criterion: amortization, SBC, registration-UQ, spatial map. Prose
  may be draft-quality but claims must be accurate and defensible (this is the positioning the
  reviewers will scrutinize). Frame honestly — v2.0's calibration is "under the simulator" and
  the fixed-summary OOD blind spot is named, consistent with the project's honesty posture.

### Figure specifications
- **D-10:** Each planned figure gets a **spec entry** (not a rendered figure): figure ID, working
  title, panel list, the phase that produces its data, the claim row(s) it supports, and a status
  placeholder. Enumerate at minimum the figures implied by the claim table (e.g., speedup-vs-RMSE,
  SBC rank histograms + coverage, OOD ROC, cross-method comparison, registration-UQ widening,
  spatial map, three-hypothesis BF). Later phases (esp. 16) fill the actual panels.

### Bibliography
- **D-11:** Use a BibTeX `refs.bib` (Typst-native support). Seed with the load-bearing citations:
  Tapqir, Costes randomization, Manders coefficients, SBC (Talts et al. 2018), NeuralEstimators.jl,
  BayesFlow/amortized SBI, and the existing ProteinCoLoc Scientific Reports / `CITATION.bib` self-cite.
  Missing DOIs/fields may be left as clearly-marked TODO placeholders rather than blocking compile.

### Compile gate (Definition of Done)
- **D-12:** The phase's hard DoD is **`typst compile manuscript/main.typ` succeeds** producing a
  PDF, with the claim table and related-work matrix rendering. This compile check should be
  runnable/scriptable so later phases can regression-check that the skeleton still builds.

### Claude's Discretion
- Exact Typst version-pinning approach (note installed 0.15.0; a simple version comment is fine).
- Precise section ordering nuances and stub wording.
- Whether the claim table lives inline in a `.typ` or in a small data file rendered by a Typst
  function — planner picks based on Typst ergonomics.
- Exact figure count/IDs beyond the claim-table-implied minimum.
</decisions>

<specifics>
## Specific Ideas

- The claim table is the mechanism that makes this phase valuable: it is the contract that
  "shapes experiments by the claims they must support." Treat it as the primary deliverable,
  with sections and figure specs hanging off it.
- Honesty posture must be visible in the skeleton itself: reserve explicit slots for the
  "calibrated under the simulator" caveat and the "summary-orthogonal OOD blind spot" so later
  phases cannot quietly drop them (mirrors STATE.md Phase-5 concerns and the Phase-4 speedup caveat).
- The >100× speedup headline must carry the Phase-4 honest caveat (median 325× forward-pass,
  ~16× on the full posterior-sample workload; realized ADVI baseline ~0.5s/pair) — the claim-table
  row and any speedup figure spec should reference this, not the bare headline.
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & claims
- `.planning/ROADMAP.md` §"Phase 10" — goal, three success criteria (Typst skeleton compiles;
  Tapqir differentiation across amortization+SBC+registration-UQ+spatial-map; figure specs per
  later phase).
- `.planning/ROADMAP.md` §"Phase 8/9/11/12/13/14/15/16" — the downstream phases whose deliverables
  the claim table and figure specs must anchor (comparator harness = Phase 9; physical corpus =
  Phase 8; registration-UQ = Phase 11; spatial map = Phase 12; 3-way BF = Phase 13; decision layer
  = Phase 14; operating envelope = Phase 15; assembly = Phase 16).
- `.planning/REQUIREMENTS.md` — the v1 spike requirements (NPE/SBC/BF/OOD) whose results become
  claim-table rows; note Phase 10 has no dedicated REQ IDs yet (TBD in roadmap) — claims derive
  from existing REQ outcomes.
- `.planning/PROJECT.md` §"Core Value" and §"Key Decisions" — the honesty posture and headline
  claims the skeleton must reserve slots for.

### Existing assets
- `CITATION.bib` (repo root) — existing ProteinCoLoc self-citation to fold into `manuscript/refs.bib`.

### Notes
- No external SPEC/ADR docs exist for this phase. The external planning input referenced in
  STATE.md (`next_project_analysis/11_plan_proteincoloc_v2.md`) is NOT present in the repo — do
  not block on it; ROADMAP.md + REQUIREMENTS.md + PROJECT.md are the authoritative sources.
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **Typst 0.15.0** installed and on PATH — the compile toolchain exists; no new install needed.
- **`CITATION.bib`** at repo root — seed entry for the bibliography.
- **Phase 4/5 metrics** (in STATE.md and phase SUMMARYs) — the honest numbers/caveats that
  populate claim-table status and figure specs (speedup caveat, SBC pre-registration concern).

### Established Patterns
- **Decoupling discipline** (project-wide): new work lives in its own tree and must not edit
  `src/`, `bayes.jl`, `colocalization.jl`, or the existing manuscript pipelines. `manuscript/` as
  a new top-level dir follows the same pattern the spike used for `spike/`.
- **Phase artifacts as contracts**: prior phases used pre-registered constants and gate tests;
  the claim table + compile gate play the analogous "contract + gate" role for the manuscript.

### Integration Points
- The claim table and figure specs are the integration surface: Phases 8, 9, 11–16 append prose,
  numbers, and rendered figures into slots this phase reserves. Modular per-section `.typ` files
  keep those later appends collision-free.
</code_context>

<deferred>
## Deferred Ideas

- **Locking a target journal + applying its Typst template** — deferred; skeleton stays
  venue-neutral until the paper is closer to submission (Phase 16 territory).
- **Writing actual Results/Methods prose and inserting real figures** — belongs to the phases
  that produce the data (5–9, 11–15) and final assembly (Phase 16).
- **CI regression gate that compiles the manuscript on every change** — a good idea, but the
  broader CI/calibration gate is Phase 15's scope; this phase only needs a runnable local compile
  check.

None of these are in scope for Phase 10.
</deferred>

---

*Phase: 10-manuscript-skeleton-and-related-work-positioning*
*Context gathered: 2026-07-02*
