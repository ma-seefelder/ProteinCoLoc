# Phase 8: External Physical Ground-Truth Corpus - Context

**Gathered:** 2026-07-20
**Status:** Ready for planning

<domain>
## Phase Boundary

Assemble a **non-circular** validation corpus whose truth does NOT come from the v2.0 simulator, so calibration can be checked against external reality rather than self-consistency. Scope is fixed by ROADMAP.md Phase 8:

1. ≥1 physical 100%-coloc anchor (positive) + ≥1 segregated anchor (negative), archived with provenance + content hashes.
2. The Colocalization Benchmark Source (CBS) ingested and explicitly flagged simulated/secondary, never conflated with the physical anchors.
3. A versioned validation data contract + manifest (schema, provenance, split policy) checked in.

This phase acquires/curates and contracts external data. It does NOT train models, run evaluation, or edit `src/` model code — blind evaluation against this corpus is Phase 16.
</domain>

<decisions>
## Implementation Decisions

### Anchor sourcing (non-circularity)
- **D-01:** The physical 100%-coloc anchor is a **named published tandem-fluorophore construct** dataset, pinned by DOI/accession as THE canonical positive anchor. Truth is asserted by biology + citation, never by any computed colocalization score (that is what keeps it non-circular).
- **D-02:** The segregated (negative) anchor is a **distinct-compartment two-fluorophore construct** — two markers targeted to non-overlapping compartments (e.g. nuclear-vs-membrane, or mito-vs-nucleus), pinned by DOI/accession. **Prefer sourcing it from the same archive/study as the positive anchor** when a matched segregated construct exists (controls imaging-condition/instrument confounds); otherwise use a separate named dataset.
- **New lab acquisition was explicitly rejected** — out of scope for a reproducible validation-corpus phase (non-reproducible-from-repo, slow).

### Vendoring policy
- **D-04:** Anchor (and CBS) image **bytes stay OUT of git**. The deliverable is a manifest (accession/DOI/URL) + a **seeded fetch script** + **SHA-256 content hashes** that verify integrity on download. Keeps the repo small and license-safe while staying reproducible.
- **D-05:** Integrity/robustness semantics: a **content-hash MISMATCH is a HARD error** (tamper/corruption — abort). A **source being UNREACHABLE emits a graceful skip-with-flag** so offline/CI runs stay green — mirrors the Phase-9 Tapqir bridge skip-with-flag pattern.

### CBS ingestion boundary (SC2)
- **D-06:** CBS lives in its **own manifest tier/namespace** with a **mandatory provenance field**: `simulated-secondary` for CBS vs `physical-primary` for the anchors. Any code path that mixes tiers must **opt in explicitly** — structural separation, not a comment. This is the mechanism that satisfies SC2 (never conflated).
- **D-07:** Ingest the **FULL CBS benchmark corpus** (spanning its labelled colocalization degrees), via the same manifest+fetch+hash mechanism, flagged `simulated-secondary`.

### Data contract + splits
- **D-08:** The versioned contract **reuses the Phase-3 / Phase-9 content-addressed pattern**: a human-readable manifest (CSV, with DataFrames/CSV already promoted in Phase 9) carrying provenance/license/tier/hash/split, a **SHA-256 content hash over data-defining bytes** for naming + verification, and JLD2 for any cached arrays.
- **D-09:** Split policy — the **physical anchors are a SEALED blind holdout**: never touched during any model development/tuning, opened only at Phase-16 evaluation. CBS may carry a seeded dev/eval split. Enforced by a manifest **split field + a guard**, mirroring the Phase-3 structural-holdout discipline (e.g. structurally-disjoint indexing, not a convention).

### Claude's Discretion
- **D-03:** The concrete **acceptance predicate** for qualifying a candidate dataset as an anchor is left to the planner/researcher (user said "you decide"). Recommended default: require (a) an unambiguous biological coloc/segregation label from the source publication AND (b) an open/redistributable license; **prefer** raw two-channel data that is retrievable and convertible to a `MultiChannelImage` so `src/` summary functions apply unchanged and the anchors are directly usable in Phase 16.
- Manifest serialization details (exact column set, TOML vs CSV for the human layer) within the D-08 pattern.
</decisions>

<specifics>
## Specific Ideas

- Positive anchor concept: a documented tandem-dimer / single-protein-two-channel fusion (biology forces 100% colocalization).
- Negative anchor concept: distinct-compartment fusion; nuclear-vs-membrane is the ROADMAP's stated example of a mutually-exclusive pair.
- Prefer a **matched pair from one study/archive** to minimize cross-dataset imaging-condition confounds between the positive and negative anchors.
- Reproducibility framing throughout mirrors the spike's anti-snooping + content-addressing discipline (SHA-256 naming, seeded selection, sealed holdout).
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

No external specs/ADRs exist for this phase (Requirements = TBD in ROADMAP.md; no `Canonical refs:` listed). Requirements are captured in the decisions above and the ROADMAP success criteria.

### Requirements source
- `.planning/ROADMAP.md` §"Phase 8: External Physical Ground-Truth Corpus" — the three success criteria this phase must make TRUE.

### In-repo patterns to reuse (see code_context)
- `.planning/ROADMAP.md` §"Phase 9: Cross-Method Comparator Harness" — the skip-with-flag / seeded shared-input / content-addressed output patterns this phase mirrors.
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **Phase-3 content-addressed cache** (`src/amortized/datagen.jl` + the Phase-3 loader): SHA-256 content hash over data-defining bytes + canonical(config) for cache naming; atomic `.tmp → integrity-check → mv(force=true)` writes; **structural holdout** via disjoint (negative) indexing so leakage is impossible by construction. Reuse the hashing, atomic-write, and structural-holdout mechanics for the manifest + sealed-anchor holdout (D-08/D-09).
- **Phase-9 comparator harness** (`run_comparator` entry point, `config.jl`, `test_comparator.jl`): the **skip-with-flag graceful-degradation** pattern (Tapqir bridge over an isolated sub-env), seeded reproducibility, and content-addressed CSV/JLD2 outputs + audit. Reuse skip-with-flag for D-05 unreachable-source handling and the content-addressed output convention for D-08.
- **DataFrames.jl / CSV.jl** already promoted to direct deps in Phase 9 — available for the manifest table (D-08).
- **Random123 (Philox)** counter-based seeding pattern — for any seeded CBS subset selection / split assignment.
- **`src/` `patch()` / `correlation()` summary functions** at fixed 8×8 (unchanged since Phase 2) — anchors should convert to `MultiChannelImage` so these apply unchanged (D-03 preference).

### Established Patterns
- **Content-addressing + atomic writes + integrity check** is the repo-wide durability idiom (Phases 3, 5, 7) — the manifest/fetch layer should follow it (hash mismatch = hard fail, D-05).
- **Structural (not conventional) separation** for anti-snooping — holdouts are enforced by construction (disjoint indexing / distinct namespaces), not by comments. Applies to both the sealed anchor holdout (D-09) and the physical/simulated tier split (D-06).
- **Isolated sub-environment + graceful skip** for external/optional integrations (Phase-9 Tapqir) — a model for any network/tooling dependency here.

### Integration Points
- **Phase 16 (External Validation + Manuscript Assembly)** is the primary consumer: it opens the sealed physical-anchor holdout for blind evaluation and uses CBS for calibration coverage.
- **Phase 9 comparator harness** accepts `Vector{MultiChannelImage}` — anchors converted per D-03 can feed it directly.
</code_context>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope. (New wet-lab acquisition was raised as an option and explicitly rejected as out of scope, not deferred.)
</deferred>

---

*Phase: 08-external-physical-ground-truth-corpus*
*Context gathered: 2026-07-20*
