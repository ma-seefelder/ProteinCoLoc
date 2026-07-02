# Phase 9: Cross-Method Comparator Harness - Context

**Gathered:** 2026-07-02 (--auto mode; recommended defaults auto-selected)
**Status:** Ready for planning

<domain>
## Phase Boundary

A reproducible, **spike-local** harness that runs the classical colocalization
estimators (**Costes-p, Manders M1/M2, Pearson, Spearman**) and a **Tapqir
bridge** on **shared inputs**, then emits a **per-method comparison table** — so
v2.0 can be positioned as **"knows when the classics are wrong,"** not merely
"agrees with them."

**In scope:** the estimator battery on one shared input set; a Tapqir bridge that
reproduces a *published* Tapqir example as a sanity anchor; a seeded/reproducible
comparison-table artifact reusing the BayesInteractomics comparator/audit pattern.

**Out of scope (other phases):** building the external physical corpus (Phase 8),
manuscript positioning prose (Phase 10), the three-hypothesis BF (Phase 13), the
decision/abstention layer (Phase 14), and blind external evaluation (Phase 16).
This phase builds the *harness and table*, not the paper's final head-to-head run.

**Hard constraint (project CLAUDE.md):** Phase 9 runs **before a Phase-6 Go**, so
it stays entirely inside `spike/` with the spike's own `Project.toml`/`Manifest`.
`src/` is consumed **read-only** (as the spike already does via `include`/develop);
**no `src/` edits, no root `Project.toml` edits.**
</domain>

<decisions>
## Implementation Decisions

### Shared input source & format
- **D-01:** The shared input container is the existing **`MultiChannelImage`**
  (2-channel), built via `spike/contract.jl:build_mci`, so every estimator AND the
  NPE consume byte-identical inputs. This is the "shared inputs" requirement (SC1).
- **D-02:** Primary shared input set is **simulator-generated** via
  `spike/simulator/forward.jl:simulate_pair` over a **seeded θ grid** spanning the
  colocalization regimes (coloc / random / exclusion, i.e. across the ρ_true range
  already used in the spike). The simulator's θ gives each row a **ground-truth
  regime label** — the anchor that lets the table show *when* a classic is wrong.
- **D-03:** Keep the harness input-source-agnostic: accept any `Vector{MultiChannelImage}`
  so the Phase-8 external corpus / CBS can be fed later **without redesign**, but do
  NOT take a Phase-8 dependency now (Phase 9 is parallelizable-now).

### Classical estimators — reuse vs new
- **D-04:** **Reuse, do not reimplement**, what already exists:
  - **Manders M1/M2** — reuse the exact formula in `spike/data/encode.jl:encode_aug`
    (Otsu thresholds from `mci.otsu_threshold`); factor into a callable rather than
    copy-paste so the comparator and the augmented encoder stay consistent.
  - **Pearson / Spearman** — reuse `src/colocalization.jl:correlation` (method
    `:pearson` / `:spearman`; supports patch-grid) and/or whole-image `cor`; report
    the **whole-image** scalar per method for the table, matching `encode_aug`'s
    `pearson_whole`.
- **D-05:** **Implement new (missing everywhere): Costes-p** — the Costes
  randomization significance p-value (block/pixel-scramble null → fraction of
  scrambles with correlation ≥ observed). New code lives **only in the spike-local
  comparator module** (decoupling). Seed the randomization with the spike RNG.
- **D-06:** All estimator wrappers live in a **new spike module**, e.g.
  `spike/comparator/classical.jl`, each returning a NaN-safe scalar with the same
  degenerate-input guarding style as `encode_aug` (`_safe`, guarded reductions).

### Tapqir bridge scope
- **D-07:** The Tapqir bridge is a **minimal sanity anchor only** (SC2): it
  reproduces a **published Tapqir example/tutorial dataset** and checks the
  recovered quantity against the published value. It is **NOT** a comparator column
  run on our simulator images — Tapqir models CoSMoS single-molecule *spot* data,
  a different data regime from our diffuse 2-channel coloc, so forcing it onto our
  inputs would be invalid. It anchors "we can invoke Tapqir correctly," nothing more.
- **D-08:** Integrate via **PythonCall/CondaPkg**, isolated in the spike (the
  CondaPkg pattern already exists in BayesInteractomics; documented fallback path in
  project stack notes). The bridge must **degrade gracefully**: if the Python/Tapqir
  env is unavailable the harness **skips the Tapqir anchor with a clear flag** and
  the classical battery still runs green (mirrors the CPU/CUDA graceful-degradation
  constraint). Tapqir is not a hard gate on the classical table.

### Comparison table schema & "knows when classics are wrong"
- **D-09:** Emit a **tidy per-input table** (`DataFrame`): one row per shared input,
  columns = `{input_id, regime/θ ground-truth, Costes_p, M1, M2, Pearson, Spearman,
  NPE_ρ̂ (or Δρ), NPE_OOD_flag}`. Persist as a versioned artifact (CSV for humans +
  JLD2 for exact reload), seeded and content-addressed like the Phase-3 cache.
- **D-10:** Add the positioning payload: a **disagreement/divergence indicator** —
  a derived column/flag marking rows where a classical verdict diverges from the
  simulator ground-truth (and, where available, from the NPE), i.e. the concrete
  "classic is wrong here" cases. Present with the **BayesInteractomics traffic-light
  pattern** (green/amber/red) rather than a bare number.
- **D-11:** Reuse the **BayesInteractomics comparator/audit table + diagnostics
  pattern** (`_bin_calibration`/`CalibrationResult`/traffic-light thresholds;
  `_ks_test_uniform`/`model_diagnostics` reporting style) for table assembly and the
  audit summary, per SC3.

### Reproducibility & harness structure
- **D-12:** Single seeded entry point `spike/comparator/run_comparator.jl` that
  builds the shared inputs, runs all estimators + the (optional) Tapqir anchor, and
  writes the table artifact. Seed with **Random123 (Philox)** exactly as the rest of
  the spike, so the table is bit-reproducible across runs and thread counts.
- **D-13:** A spike test (wired into the existing `spike/test` gate) asserts:
  (a) determinism of the table across two seeded runs; (b) Costes-p, M1, M2, Pearson,
  Spearman all emit **finite** values on a non-degenerate shared fixture (SC1);
  (c) the Tapqir bridge either reproduces the published anchor within tolerance OR
  is cleanly skipped-with-flag when the env is absent (SC2).

### Anti-data-snooping guard (carried from roadmap blocker)
- **D-14:** The "knows when classics are wrong" thresholds (divergence cutoff,
  traffic-light bands) must be **fixed/pre-declared** in the harness config, not
  tuned post-hoc against the table — same pre-registration discipline flagged for
  Phase 5. Document the chosen thresholds and their rationale in the artifact header.

### Claude's Discretion
- Exact Costes randomization scheme (whole-image pixel scramble vs block/Van Steensel
  shift) and number of scrambles — pick the standard, cite it, keep it seeded.
- Table column ordering, CSV vs additional Arrow output, figure rendering (if any).
- Module/file layout under `spike/comparator/` and how estimators are registered.
</decisions>

<specifics>
## Specific Ideas

- The framing sentence to satisfy is literally the phase goal: the table must let a
  reader point at rows and say "here the classical estimator disagrees with the
  ground-truth regime, and v2.0's NPE/OOD flag catches it." Design the divergence
  column to make those rows obvious.
- Manders in `encode_aug` and the comparator MUST use the **same** Otsu-threshold
  source (`mci.otsu_threshold`) so the ablation-augmented summary and the comparator
  never silently diverge on M1/M2.
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase definition & requirements
- `.planning/ROADMAP.md` §"Phase 9: Cross-Method Comparator Harness" (lines ~172–181) — goal, the 3 success criteria (SC1 per-method table on shared input; SC2 Tapqir published-example anchor; SC3 reuse BayesInteractomics comparator/audit pattern + seeded/reproducible).
- `.planning/REQUIREMENTS.md` — Phase 9 requirements are **TBD** (planner should propose CMP-* IDs); `ABL-01` documents Manders/median/IQR as sufficiency diagnostics already exercised.
- `.planning/10_plan_amortizedcoloc.md` (lines ~19–27) — BayesInteractomics reuse table (`_bin_calibration`/`CalibrationResult`, simulation/audit engine) and the PythonCall/CondaPkg fallback precedent for the Tapqir bridge.

### Decoupling constraint
- `CLAUDE.md` §Constraints (Decoupling; Platform graceful-degradation) — spike-only, no `src/` edits, CPU/optional-Python must degrade gracefully.

### Reusable estimator/summary code (this repo)
- `spike/data/encode.jl:encode_aug` (lines ~96–137) — canonical Manders M1/M2 + whole-image Pearson formulas and NaN-safe reduction style to factor out and reuse.
- `src/colocalization.jl:correlation` (line ~221) — Pearson/Spearman/Kendall over the patch grid (read-only reuse); `patch`/`unpatch` helpers.
- `spike/contract.jl:build_mci` / `patch_summary` / `induced_mu` (lines ~50–100) — the `MultiChannelImage` shared-input constructor and summary contract.
- `spike/simulator/forward.jl:simulate_pair` (line ~110) — seeded shared-input generator with ground-truth θ.

### Comparator / audit pattern (sibling repo, read-only template)
- `~/Documents/GitHub/BayesInteractomics/src/diagnostics/calibration.jl` — `_bin_calibration`, `_compute_calibration` (ECE/MCE + traffic-light).
- `~/Documents/GitHub/BayesInteractomics/src/diagnostics/types.jl` §`CalibrationResult` (line ~227).
- `~/Documents/GitHub/BayesInteractomics/src/diagnostics/predictive_checks.jl` — `model_diagnostics`, `_ks_test_uniform`, `generate_diagnostics_report` (audit-table/reporting style to mirror).
- `~/Documents/GitHub/BayesInteractomics/src/simulation/simulation.jl` — scenario-grid sweep / replicate / JLD2-cache template for the seeded harness.
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **Manders M1/M2 + whole-image Pearson**: already implemented in
  `spike/data/encode.jl:encode_aug` — factor into a shared callable; do not re-derive.
- **Pearson/Spearman/Kendall**: `src/colocalization.jl:correlation` (method kwarg),
  read-only via the spike's existing `include`/develop of the main package.
- **Shared input plumbing**: `spike/contract.jl` (`build_mci`, `patch_summary`) +
  `spike/simulator/forward.jl` (`simulate_pair`) give seeded `MultiChannelImage`s
  with ground-truth θ — the natural "shared inputs" and regime labels.
- **Audit/table pattern**: BayesInteractomics diagnostics (calibration, traffic-light,
  KS-uniformity, scenario-sweep + JLD2 cache) — direct template for SC3.
- **Seeding**: Random123/Philox pattern already established in Phase 3
  (`spike/data/*`) — reuse for bit-reproducibility.

### Established Patterns
- **Decoupling**: everything new goes under `spike/` (new `spike/comparator/`);
  `src/` stays untouched — same discipline as Phases 1–4.
- **NaN-safe reductions** (`_safe`, guarded `skipmissing`) for degenerate inputs —
  mirror `encode_aug`.
- **Graceful degradation** for optional external stacks (CUDA today; Python/Tapqir
  here) — skip-with-flag, never a hard failure of the core path.
- **Content-addressed, atomic, resume-safe artifacts** (Phase 3 JLD2 cache) —
  reuse for the comparison-table artifact.

### Integration Points
- Consumes `simulate_pair` outputs and `MultiChannelImage` (shared with the NPE path).
- Reuses `src/colocalization.jl` and `spike/data/encode.jl` estimator math.
- New Python edge: PythonCall/CondaPkg env for the Tapqir bridge only, isolated and optional.
- Feeds Phase 10 (related-work delta) and Phase 16 (blind external eval) downstream.
</code_context>

<deferred>
## Deferred Ideas

- **Run classical estimators + NPE on the Phase-8 physical corpus / CBS** for the
  final head-to-head — Phase 16 (external validation); Phase 9 only builds the
  harness and validates it on simulator inputs.
- **Costes automatic threshold (Costes regression) determination** beyond Otsu —
  keep Otsu for M1/M2 now; revisit only if the comparator needs it.
- **Full Tapqir integration as a live comparator column** — explicitly rejected as
  invalid for our data regime (D-07); anchor only.
- **Manuscript prose framing of the "when classics are wrong" delta** — Phase 10.

### Reviewed Todos (not folded)
- None — roadmap "Pending Todos" is empty.
</deferred>

---

*Phase: 09-cross-method-comparator-harness*
*Context gathered: 2026-07-02*
