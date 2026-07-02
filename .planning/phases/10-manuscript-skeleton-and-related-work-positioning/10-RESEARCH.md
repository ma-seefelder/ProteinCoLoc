# Phase 10: Manuscript Skeleton and Related-Work Positioning - Research

**Researched:** 2026-07-02
**Domain:** Scientific manuscript authoring (Typst 0.15.0), related-work positioning, bibliography engineering, compile-gate scripting
**Confidence:** HIGH (Typst syntax verified against the installed binary; citations verified against publisher DOIs; related-work deltas partially flagged for verification)

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Author in **Typst** (Typst **0.15.0** installed). The skeleton must actually compile with the installed binary — a broken skeleton is not a deliverable.
- **D-02:** Keep the skeleton **venue-neutral** — a generic single-column Typst `article` scaffold, not a specific journal template. Record candidate venues (Bioinformatics, Nature Methods, eLife-style methods paper) as a comment/note only; do not lock one.
- **D-03:** Create a **new top-level `manuscript/` directory** at repo root (sibling to `spike/`, `src/`, `test/`). Do **not** place manuscript sources under `spike/` or `src/`. Honors decoupling: `src/`, `bayes.jl`, `colocalization.jl` and both manuscript pipelines stay provably untouched.
- **D-04:** Suggested internal layout: `manuscript/main.typ` (entry), `manuscript/sections/*.typ` (per-section includes), `manuscript/refs.bib`, `manuscript/claims.typ` (or a data file feeding the claim table), `manuscript/figures/` (specs + later PDFs). Section files must be modular so later phases append prose without merge collisions.
- **D-05:** Scaffold standard methods-paper sections as empty-but-labelled stubs with a one-line intent comment each: Abstract, Introduction, Related Work, Methods (forward simulator, summary statistic, NPE/NRE, SBC, OOD), Results (placeholder subsections mirroring the claim table), Discussion, Data/Code Availability, References. Each stub carries a `// TODO(phase N)` marker.
- **D-06:** The claim table is the **spine** — a structured table mapping **scientific claim → differentiator axis → supporting phase(s) → figure/experiment ID → status**. Legible both as a rendered Typst table and as a plain source block a human/agent can diff. Seed with v2.0 headline claims each linked to its owning phase.
- **D-07:** Every claim row must name the artifact that will substantiate it (a figure spec ID and/or a phase deliverable).
- **D-08:** Draft related work as **comparison matrix + prose**. Matrix rows: ProteinCoLoc v2.0, Tapqir, Costes randomization, Manders M1/M2, Pearson/Spearman, original ProteinCoLoc v1 ADVI. Columns: amortized inference, calibrated uncertainty/SBC, Bayes-factor/model comparison, registration uncertainty as latent, spatial (per-region) map, misspecification/OOD flag, speed.
- **D-09:** The **Tapqir differentiation** paragraph is mandatory and must explicitly cover all four axes: amortization, SBC, registration-UQ, spatial map. Prose may be draft-quality but claims must be accurate and defensible. Frame honestly — calibration is "under the simulator" and the fixed-summary OOD blind spot is named.
- **D-10:** Each planned figure gets a **spec entry** (not a rendered figure): figure ID, working title, panel list, the phase that produces its data, the claim row(s) it supports, status placeholder. Enumerate at minimum the figures implied by the claim table.
- **D-11:** Use a BibTeX `refs.bib` (Typst-native support). Seed load-bearing citations: Tapqir, Costes randomization, Manders coefficients, SBC (Talts et al. 2018), NeuralEstimators.jl, BayesFlow/amortized SBI, and the existing ProteinCoLoc Scientific Reports / `CITATION.bib` self-cite. Missing DOIs/fields may be clearly-marked TODO placeholders rather than blocking compile.
- **D-12:** Hard DoD is **`typst compile manuscript/main.typ` succeeds** producing a PDF, with the claim table and related-work matrix rendering. The compile check should be runnable/scriptable so later phases can regression-check the build.

### Claude's Discretion
- Exact Typst version-pinning approach (note installed 0.15.0; a simple version comment is fine).
- Precise section ordering nuances and stub wording.
- Whether the claim table lives inline in a `.typ` or in a small data file rendered by a Typst function — pick based on Typst ergonomics.
- Exact figure count/IDs beyond the claim-table-implied minimum.

### Deferred Ideas (OUT OF SCOPE)
- **Locking a target journal + applying its Typst template** — deferred until closer to submission (Phase 16 territory).
- **Writing actual Results/Methods prose and inserting real figures** — belongs to the phases that produce the data (5–9, 11–15) and final assembly (Phase 16).
- **CI regression gate that compiles the manuscript on every change** — the broader CI/calibration gate is Phase 15's scope; this phase needs only a runnable local compile check.
</user_constraints>

<phase_requirements>
## Phase Requirements

**No dedicated REQ IDs exist for Phase 10** (Requirements are TBD in ROADMAP for the v2.0 feature-expansion phases). Claims that populate the claim table derive from **existing requirement outcomes and prior-phase metrics**, not from Phase-10-specific requirements. The mapping below is the substitute traceability the planner should encode into the claim table (D-06/D-07):

| Claim (claim-table row) | Derives From | Owning / Supporting Phase | Status as of research |
|-------------------------|--------------|---------------------------|-----------------------|
| Amortized ms-scale inference (single forward pass) | NPE-01, NPE-03 | Phase 4 (done) / Phase 7 (productionize) | Supported (measured) |
| >100× speedup at comparable RMSE | NPE-02, NPE-03 | Phase 4 (done) | Supported **with caveat** (see Honest Numbers) |
| SBC calibration-under-the-simulator | SBC-01..04 | Phase 5 (pending) | **Reserved slot** (not yet proven) |
| Honest OOD / misspecification flag | OOD-01, OOD-02 | Phase 5 (pending) | Reserved slot; blind-spot named |
| Amortized Bayes factor (2-way → 3-way) | BF-01, BF-02 | Phase 5 (2-way) / Phase 13 (3-way) | Reserved slot |
| Registration uncertainty as latent | (new) | Phase 11 (pending) | Reserved slot |
| Spatial per-region colocalization map | (new) | Phase 12 (pending, descope-to-v2.1 candidate) | Reserved slot |
| Cross-method "knows when the classics are wrong" | (new) | Phase 9 (pending) | Reserved slot |
| Non-circular external validation | (new) | Phase 8 (pending) | Reserved slot |

The planner MUST make the claim table carry a **status** column that distinguishes *supported (measured)* from *reserved slot (pending phase)* so the skeleton does not overclaim. This is the mechanism that "shapes experiments by the claims they must support."
</phase_requirements>

## Summary

This is a **documentation phase**: build a compiling, modular Typst 0.15.0 manuscript skeleton under a new top-level `manuscript/` tree, with related-work positioning drafted early and a claim table as the structural spine. There is no code to run, no external packages to install, and no security surface — the single hard gate is `typst compile manuscript/main.typ` returning exit 0 and producing a PDF that renders the claim table and related-work matrix.

Every load-bearing Typst construct the planner will rely on was **verified against the installed 0.15.0 binary** in this session (modular `#include`/`#import`, native `table()` with `table.header`, `#figure` + `<label>` cross-references, `@citekey` citations, and `#bibliography("refs.bib", style: "ieee")`). Two compile-gate gotchas were discovered empirically and are the most important operational findings: **(1)** an undefined citation key or reference label makes `typst compile` **fail with exit 1** (excellent — turns a dangling `@cite` into a hard build error the regression gate catches for free), and **(2)** piping `typst compile ... | head` **masks the real exit code** — the gate script must read `typst`'s own exit status (`PIPESTATUS[0]` or no pipe at all), or it will silently pass on a broken build.

The related-work content is factually anchored: Tapqir (Ordabayev et al., eLife 2022, DOI 10.7554/eLife.73860), Costes 2004 (Biophys J, DOI 10.1529/biophysj.103.038422), and Manders 1993 (J Microsc, DOI 10.1111/j.1365-2818.1993.tb03313.x) all have verified DOIs seeded below. The **honesty posture is a hard constraint on the prose**: the ">100× speedup" is a forward-pass figure (median 325×) that drops to ~16× on the full posterior-sample workload against a realized ~0.5 s/pair ADVI baseline — the claim-table row and speedup figure spec must carry this caveat, not the bare headline.

**Primary recommendation:** Scaffold `manuscript/{main.typ, sections/*.typ, claims.typ, refs.bib, figures/, build.sh|build.jl}`. Put the claim data in `claims.typ` as a Typst array-of-dicts rendered by a `claim_table()` function (satisfies D-06's "rendered table AND diffable source" with one source of truth). Wire the compile gate as a tiny script that runs `typst compile` and checks its own exit code. Draft the Tapqir delta honestly around the one unambiguous axis (amortization) and flag the three domain-dependent axes (SBC-vs-Tapqir, registration-UQ, spatial-map) for verification against the eLife paper before locking prose.

## Architectural Responsibility Map

This is a document build, not a multi-tier app; tiers below are the manuscript's conceptual layers. The planner uses this to keep each concern in its own modular file so later phases append without collisions.

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Document config / template (page, font, headings, title) | Rendering config (`main.typ` set-rules) | — | Venue-neutral template lives once at the entry; journal template swap later touches only this (D-02) |
| Section prose stubs | Authoring layer (`sections/*.typ`) | — | One file per section → later phases edit disjoint files, no merge collisions (D-04/D-05) |
| Claim spine (data + rendered table) | Data/spine layer (`claims.typ`) | Authoring layer (Results mirrors it) | Single source of truth rendered two ways (D-06); Results subsections mirror rows (D-05) |
| Figure specifications | Data/spine layer (`figures/specs.typ` or in `claims.typ`) | Authoring layer | Specs are structured data enumerated per phase (D-10) |
| Bibliography | Bibliography layer (`refs.bib`) | Rendering (`#bibliography` call in `main.typ`) | BibTeX is the diffable citation store; Typst renders it (D-11) |
| Compile gate | Build/gate layer (`build.sh`/`build.jl` + `typst compile`) | — | The falsifiable DoD; scriptable for regression (D-12) |

## Standard Stack

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| **Typst** | **0.15.0** (installed, `typst --version` → `typst 0.15.0 (3ae52774)`) | Typeset the manuscript; `typst compile` is the DoD gate | ROADMAP-mandated (D-01); already installed; native BibTeX + tables + cross-refs; single self-contained binary (no LaTeX toolchain) [VERIFIED: `typst --version` on this host] |
| **BibTeX `.bib`** | n/a (data format) | `manuscript/refs.bib` citation store | Typst has native `.bib` support via `#bibliography("refs.bib")`; diffable plain text (D-11) [VERIFIED: compiled locally] |

### Supporting
No supporting libraries. This phase installs **nothing** — Typst is preinstalled and the deliverables are `.typ`/`.bib`/`.pdf` files plus a shell/Julia gate script.

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Typst native `table()` | LaTeX-style external table pkg | Not applicable — Typst is mandated (D-01) and has a capable native `table` |
| `claims.typ` array-of-dicts | Inline literal table in a `.typ`, or an external JSON/YAML read via `json()`/`yaml()` | Array-of-dicts keeps data+render in one Typst file (simplest, no extra file, still diffable). External JSON adds a file but is language-neutral for downstream tooling — reserve only if a non-Typst agent must machine-read the claims. Typst *can* `#json("claims.json")` if desired. |
| `style: "ieee"` bibliography | `"apa"`, `"nature"`, `"chicago-author-date"`, or a CSL file | Venue-neutral (D-02) → pick a neutral built-in (`ieee` or `apa`); a CSL swap is a one-line change later (Phase 16) |

**Installation:** None. Verify the toolchain only:
```bash
typst --version    # expect: typst 0.15.0 (...)
```

## Package Legitimacy Audit

**Not applicable — this phase installs no external packages.** Typst is a preinstalled system binary (verified 0.15.0); `refs.bib` is a data file authored in-repo. No npm/PyPI/crates dependency is added. The slopcheck / registry-verification gate has nothing to audit for Phase 10.

## Architecture Patterns

### System Architecture Diagram

```
                         manuscript/build.sh (or build.jl)
                                    │  runs
                                    ▼
   manuscript/main.typ  ──set-rules (page/font/headings/title, venue-neutral)
        │
        ├── #import "claims.typ": claims, claim_table   ── data + render fn (single source of truth)
        │         │                          │
        │         ▼                          ▼
        │   diffable array-of-dicts    rendered #table  ──► Claim-table figure <tab:claims>
        │
        ├── #include "sections/abstract.typ"        ┐
        ├── #include "sections/introduction.typ"    │
        ├── #include "sections/related_work.typ" ───┼─ per-section stubs, // TODO(phase N)
        │        (comparison matrix + Tapqir delta) │   each file edited by exactly one later phase
        ├── #include "sections/methods.typ"         │
        ├── #include "sections/results.typ"  ───────┤  (subsections mirror claim rows)
        ├── #include "sections/discussion.typ"      │
        ├── #include "sections/availability.typ"    ┘
        │
        ├── #include "figures/specs.typ"   ── figure spec entries (ID/title/panels/phase/claim/status)
        │
        └── #bibliography("refs.bib", style: "ieee")   ── @cite keys resolve here
                     │
                     ▼
             typst compile  ──► main.pdf   (exit 0 = DoD met; exit 1 = broken build)
```

### Recommended Project Structure
```
manuscript/
├── main.typ              # entry: set-rules + imports/includes + #bibliography
├── claims.typ            # claim data (array of dicts) + claim_table() render fn (D-06)
├── refs.bib              # BibTeX citation store (D-11)
├── sections/
│   ├── abstract.typ          # // TODO(phase 16): final abstract
│   ├── introduction.typ      # // TODO(phase 16)
│   ├── related_work.typ      # comparison matrix + Tapqir delta (drafted THIS phase, D-08/D-09)
│   ├── methods.typ           # subsections: simulator/summary/NPE-NRE/SBC/OOD  // TODO(phase 7,5)
│   ├── results.typ           # subsections mirror claim rows  // TODO(phase 5,8,9,11-15)
│   ├── discussion.typ        # // TODO(phase 16)
│   └── availability.typ      # data/code availability  // TODO(phase 16)
├── figures/
│   ├── specs.typ             # figure spec entries (D-10)
│   └── (later: F1.pdf, F2.pdf … dropped in by phase 16)
└── build.sh                  # (or build.jl) compile gate (D-12)
```

### Pattern 1: Data-driven claim table (rendered + diffable, one source) — D-06
**What:** Store claim rows as a Typst array of dicts; render with a function. The array *is* the diffable plain-source block; the function *is* the rendered table. Optionally also emit a `#raw(...)` verbatim block so the diffable form appears in the PDF too.
**When to use:** The claim spine and the figure-spec list.
**Example (VERIFIED — compiled to PDF, exit 0):**
```typst
// claims.typ
#let claims = (
  (id: "C1", claim: "Amortized ms-scale inference", axis: "Amortization",
   phase: "4,7", fig: "F1", status: "supported"),
  (id: "C2", claim: ">100x speedup at comparable RMSE (caveated)", axis: "Speed",
   phase: "4", fig: "F1", status: "supported (caveat)"),
)

#let claim_table(claims) = table(
  columns: 6,
  table.header([ID], [Claim], [Axis], [Phase], [Fig], [Status]),
  ..claims.map(c => (c.id, c.claim, c.axis, c.phase, c.fig, c.status)).flatten(),
)
```
```typst
// in main.typ
#import "claims.typ": claims, claim_table
#figure(claim_table(claims), caption: [Claim spine.]) <tab:claims>
See @tab:claims.
```

### Pattern 2: Modular sections with `#include` — D-04/D-05
**What:** `main.typ` composes the document by `#include`-ing one file per section. `#include` splices the file's *content* at that point (headings render in order). `#import` pulls *definitions* (functions/variables) without emitting content.
**When to use:** Section prose (`#include`), shared data/functions (`#import`).
**Example (VERIFIED):**
```typst
#include "sections/introduction.typ"      // emits the section's content here
#import "claims.typ": claims, claim_table // pulls definitions, emits nothing
```
**Why it prevents merge collisions:** each later phase edits a different `sections/*.typ` file; `main.typ` only lists includes and rarely changes.

### Pattern 3: Figure spec as a stub `#figure` (no image yet) — D-10
**What:** Reserve each figure's slot with a captioned placeholder that carries the spec (ID, title, panels, producing phase, supported claim rows, status), so a later phase swaps in the real image without renumbering.
**Example:**
```typst
#figure(
  rect(width: 100%, height: 4cm, stroke: (dash: "dashed"))[
    *SPEC F1 — Speedup vs RMSE* \
    Panels: (a) wall-clock NPE vs ADVI; (b) RMSE parity. \
    Data: Phase 4. Supports: C1, C2. Status: PENDING.
  ],
  caption: [Placeholder — Phase 16 inserts the rendered panel.],
) <fig:speedup>
```

### Pattern 4: Cross-references and citations — VERIFIED
```typst
See @tab:claims and @fig:speedup.     // label refs
Tapqir uses a physics-based model @ordabayev2022.  // bib citation
```
All three resolve at compile; an **unresolved** `@key` fails the build (exit 1) — see Pitfall 2.

### Anti-Patterns to Avoid
- **Monolithic `main.typ`** holding all prose → guarantees merge collisions when phases 5,8,9,11–16 all append. Use `sections/*.typ`.
- **Hand-maintaining the claim table as literal table cells** → the data and its rendering drift. Use the array-of-dicts + function (Pattern 1).
- **Locking a journal template now** → violates D-02 and creates rework. Keep set-rules generic; note candidate venues in a comment only.
- **Real numbers/prose in Results** → out of scope (deferred); Results holds only mirrored subsection stubs.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Bibliography formatting | Manual reference list, hand-numbered `[1]` | `#bibliography("refs.bib", style: "ieee")` + `@key` | Native; auto-numbers, sorts, and only lists cited entries; style swap is one line [VERIFIED] |
| Table rendering from data | Literal cell-by-cell tables | `table()` + `.map().flatten()` over a data array | One source of truth, diffable, renderable (Pattern 1) [VERIFIED] |
| Cross-reference numbering | Manual "Table 1/Figure 2" text | `<label>` + `@label` | Auto-renumbers; dangling refs become build errors (free regression check) [VERIFIED] |
| Compile-success detection | Grepping stdout for "error" | `typst compile`'s **exit code** | Exit 0 = success, exit 1 = failure, verified for clean/undefined-cite/missing-file cases |

**Key insight:** Typst's exit code and its treatment of unresolved references as hard errors give you a *falsifiable* build gate almost for free — the discipline is to (a) never pipe away the exit code and (b) put every claim/figure behind a labelled `<...>` so a broken link fails the build.

## Common Pitfalls

### Pitfall 1: Piping `typst compile` masks the exit code (compile gate lands broken)
**What goes wrong:** `typst compile main.typ | head` (or `| tee log`) makes the shell report `head`'s exit status (0), not `typst`'s. A broken build silently "passes" the gate.
**Why it happens:** In a pipeline `$?` is the last command's status. Observed directly this session: `typst compile does_not_exist.typ 2>&1 | head` then `$?` = 0, but the direct `typst compile does_not_exist.typ` = **exit 1**.
**How to avoid:** In the gate script, run `typst compile` **without a pipe** and test `$?`, or capture `${PIPESTATUS[0]}` immediately. In Julia, use `run(`typst compile ...`)` (throws on nonzero) or `success(cmd)`.
**Warning signs:** A gate that "passes" but no `main.pdf` was produced / updated.

### Pitfall 2: An undefined `@citekey` or `<label>` reference fails the whole build (exit 1)
**What goes wrong:** A stub that cites `@tapqir2022` while `refs.bib` uses key `ordabayev2022`, or a `@fig:foo` with no `<fig:foo>` anywhere, aborts compilation.
**Why it happens (VERIFIED):** `error: label <...> does not exist in the document` → exit 1.
**How to avoid:** Keep BibTeX keys and `@` usages in lockstep. During scaffolding, either (a) don't write a `@cite`/`@ref` until its target exists, or (b) add the `refs.bib` entry / `<label>` in the same edit. This is a *feature* for regression — but it means seeded stubs must not reference not-yet-existing keys.
**Warning signs:** First compile of a fresh stub fails on a citation/label the author "meant to add later."

### Pitfall 3: Windows path handling in the compile gate
**What goes wrong:** Backslash paths or spaces in `C:\Users\Manuel\...` break naive scripts; the shell here is Git Bash (POSIX), not cmd/PowerShell.
**How to avoid:** Use forward slashes and quote paths in the gate script. `typst compile "manuscript/main.typ" "manuscript/main.pdf"`. Run the gate from repo root (relative paths keep `#include`/`#bibliography` resolution simple — Typst resolves includes relative to the *importing file*, and `refs.bib` relative to the file that calls `#bibliography`, so keeping `main.typ`, `refs.bib`, `sections/` under one tree avoids `--root` gymnastics). Absolute paths need `--root <DIR>` (env `TYPST_ROOT`) — avoid by staying repo-relative.
**Warning signs:** Compiles from repo root but not from another CWD → path assumption baked in.

### Pitfall 4: `#import` vs `#include` confusion
**What goes wrong:** `#include "claims.typ"` dumps the file's *content* (any top-level `#table(...)` calls render at the include site); `#import "claims.typ": ...` pulls *definitions* only.
**How to avoid:** Put only `#let` definitions in `claims.typ` (no bare content) and `#import` it; use `#include` for `sections/*.typ` which are meant to emit prose.

### Pitfall 5: Overclaiming in the seeded skeleton (honesty posture)
**What goes wrong:** The claim table or Tapqir prose asserts SBC/OOD/registration/spatial results as done when Phases 5/11/12 are pending, or states the bare ">100×" without the caveat.
**How to avoid:** Status column distinguishes *supported (measured)* vs *reserved slot*; the speedup row carries the Phase-4 caveat; reserve explicit slots for "calibrated under the simulator" and the "summary-orthogonal OOD blind spot" (Blockers in STATE.md). Flag Tapqir-comparison cells that are not yet verified against the eLife paper (see Related Work below).

## Code Examples

All examples below were **compiled with the installed Typst 0.15.0 and produced a PDF at exit 0** in this research session.

### Minimal compiling `main.typ` (venue-neutral scaffold)
```typst
// Typst 0.15.0 — venue-neutral single-column article
// Candidate venues (NOT locked, D-02): Bioinformatics | Nature Methods | eLife methods
#set document(title: "ProteinCoLoc v2.0 — AmortizedColoc", author: "M. Seefelder")
#set page(paper: "a4", numbering: "1")
#set text(size: 10pt)
#set heading(numbering: "1.1")

#import "claims.typ": claims, claim_table

#include "sections/abstract.typ"
#include "sections/introduction.typ"
#include "sections/related_work.typ"

= Claims
#figure(claim_table(claims), caption: [Claim spine.]) <tab:claims>

#include "sections/methods.typ"
#include "sections/results.typ"
#include "sections/discussion.typ"
#include "sections/availability.typ"

#bibliography("refs.bib", style: "ieee")
```

### Diffable-block companion to the rendered table (optional, satisfies D-06 "plain source block")
```typst
// renders the same claim data as a verbatim block a human/agent can diff in the PDF
#raw(
  claims.map(c => c.id + " | " + c.claim + " | " + c.axis + " | " + c.status).join("\n"),
  block: true,
)
```
*(Note: the array-of-dicts in `claims.typ` is itself the primary diffable artifact in source control; this `#raw` block is only needed if the diffable view must also appear in the rendered PDF.)*

### Compile gate (`manuscript/build.sh`) — exit-code-correct (Pitfall 1)
```bash
#!/usr/bin/env bash
set -euo pipefail
# Run from repo root. No pipe on the compile line — preserve typst's exit code.
typst compile "manuscript/main.typ" "manuscript/main.pdf"
echo "OK: manuscript/main.pdf built (exit 0)"
```
Julia variant (fits the repo's Julia tooling; `run` throws on nonzero):
```julia
# manuscript/build.jl  —  julia manuscript/build.jl
run(`typst compile manuscript/main.typ manuscript/main.pdf`)  # throws if exit != 0
@info "manuscript/main.pdf built"
```

## Related Work — Factual Accuracy and Defensible Delta (D-08 / D-09)

### What each method actually is (verified against publisher sources)

- **Tapqir (cosmos model)** — Ordabayev, Friedman, Gelles, Theobald, *eLife* 2022;11:e73860 [CITED: elifesciences.org/articles/73860, DOI 10.7554/eLife.73860]. An **unsupervised ML** method with a **holistic, physics-based causal generative model of CoSMoS** (colocalization single-molecule spectroscopy) image data. Accounts for photon + camera noise, optical non-uniformities, non-specific binding, and spot detection; outputs **per-spot classification probabilities** (not a binary spot/no-spot call). Inference is **per-dataset variational (Pyro/SVI)**. Domain: **single-molecule TIRF spot data (AOIs around candidate spots)**.
- **Costes randomization** — Costes et al., *Biophysical Journal* 2004;86(6):3993–4003 [CITED: DOI 10.1529/biophysj.103.038422]. Automatic threshold selection (below which pixels show no correlation) + a **pixel/block randomization test** producing a **p-value** for colocalization significance (Pearson's R above threshold, compared to scrambled-image null). Frequentist significance, not a posterior.
- **Manders M1/M2** — Manders, Verbeek, Aten, *Journal of Microscopy* 1993;169(3):375–382 [CITED: DOI 10.1111/j.1365-2818.1993.tb03313.x]. Two **co-occurrence** coefficients: fraction of channel-1 signal overlapping channel-2, and vice versa. **Threshold-dependent**, intensity-co-occurrence (not correlation), **no uncertainty quantification** on its own.
- **Pearson / Spearman** — global (rank-)correlation of the two channels' intensities. Single scalar, **no UQ, no significance thresholding**, confounded by background.
- **ProteinCoLoc v1 (ADVI)** — the existing package: hierarchical Student-t Turing model, per-dataset ADVI posterior over Δρ, KDE-based Bayes factor. **Per-dataset** inference (~0.5 s/pair realized `vi()` clock; full pipeline adds a 100k prior chain + `rand(q,100k)`), no formal SBC proof, fixed registration.

### Comparison matrix (seed for `related_work.typ`) — with defensibility flags

| Method | Amortized | Calibrated UQ / SBC | BF / model comparison | Registration UQ as latent | Spatial per-region map | OOD flag | Speed |
|--------|-----------|---------------------|-----------------------|---------------------------|------------------------|----------|-------|
| **ProteinCoLoc v2.0** | ✓ (single forward pass) | ✓ *under the simulator* (SBC, Ph5) | ✓ 2-way (Ph5) → 3-way (Ph13) | ✓ (Ph11) | ✓ (Ph12) | ✓ *summary-bounded* (Ph5) | **ms** (fwd pass) |
| **Tapqir** | ✗ per-dataset SVI | ~ Bayesian per-spot prob.; **formal SBC unclear** ⚠ | ~ probabilistic classification, **not explicit BF** ⚠ | ~ models spot xy; **drift handling** ⚠ | ~ per-AOI/per-spot (single-molecule regime) ⚠ | ✗ **no explicit misspec flag** ⚠ | slow (per-dataset SVI) |
| **Costes randomization** | ✗ | ✗ (p-value only) | ✗ | ✗ | ✗ (global) | ✗ | moderate (randomization) |
| **Manders M1/M2** | ✗ | ✗ | ✗ | ✗ | ✗ (global) | ✗ | fast |
| **Pearson / Spearman** | ✗ | ✗ | ✗ | ✗ | ✗ (global) | ✗ | fast |
| **ProteinCoLoc v1 (ADVI)** | ✗ per-dataset ADVI | ~ Bayesian posterior, **no SBC proof** | ✓ KDE Δρ BF | ✗ (fixed) | ~ patch-level, no spatial prior | ✗ | ~0.5 s–min/dataset |

**⚠ = cell NOT independently verified against the Tapqir paper internals this session — [ASSUMED] from the abstract + domain knowledge. The planner MUST add a task to verify each ⚠ cell against Ordabayev et al. 2022 before locking the prose (D-09 requires accuracy under reviewer scrutiny).**

### The four-axis Tapqir differentiation (D-09) — honest framing

1. **Amortization** — *Strongly defensible.* Tapqir fits a variational posterior **per dataset** (SVI). ProteinCoLoc v2.0 trains once and infers new datasets in a **single amortized forward pass (ms)**. This is the cleanest, least-contestable delta.
2. **SBC / calibration** — *Defensible with care.* v2.0 provides a **formal SBC rank-uniformity + coverage-curve proof (under the simulator)**. Frame the delta as "provides a formal SBC coverage proof" rather than "Tapqir has no calibration" — Tapqir is Bayesian and validates on simulated data; whether it runs **formal SBC** is the ⚠ item to verify. State v2.0's calibration as *under the simulator* and pair it with OOD (honesty posture).
3. **Registration uncertainty as latent** — *Defensible with care.* v2.0 (Phase 11) promotes sub-pixel registration to an **inferred latent so the posterior widens honestly** under registration uncertainty. Tapqir models spot xy positions and typically relies on **external/preprocessing drift correction**; verify (⚠) before claiming it does not propagate registration uncertainty into the coloc posterior.
4. **Spatial (per-region) map** — *Frame as complementary domains, not "better."* v2.0 (Phase 12) yields an **amortized per-region Δρ map with per-region UQ over a dense correlation lattice**. Tapqir operates on **sparse single-molecule AOIs** — a different data regime. The defensible statement is that v2.0 targets **dense/diffuse two-channel fluorescence** and adds a spatial map there, whereas Tapqir is **single-molecule-spot-specialized**; they address overlapping but distinct problems.

**Overarching honesty point for the prose:** Tapqir and ProteinCoLoc target **partly different data regimes** (single-molecule CoSMoS spots vs general patch-correlation microscopy). Position v2.0's contribution as *amortization + formal calibration proof + spatial map + registration-UQ + OOD flag on the general/dense colocalization problem*, not as a blanket superiority claim over Tapqir. This survives reviewer scrutiny; a blanket claim will not.

## Citation Entries to Seed `refs.bib` (D-11)

BibTeX-ready. DOIs verified against publisher pages this session; **author lists and page numbers marked TODO where not fully confirmed** — clearly-marked placeholders are acceptable per D-11.

```bibtex
@article{ordabayev2022,
  title   = {Bayesian machine learning analysis of single-molecule fluorescence colocalization images},
  author  = {Ordabayev, Yerdos A. and Friedman, Larry J. and Gelles, Jeff and Theobald, Douglas L.},
  journal = {eLife},
  volume  = {11},
  pages   = {e73860},
  year    = {2022},
  doi     = {10.7554/eLife.73860}
}

@article{costes2004,
  title   = {Automatic and quantitative measurement of protein-protein colocalization in live cells},
  author  = {Costes, Sylvain V. and Daelemans, Dirk and Cho, Edward H. and Dobbin, Zachary and Pavlakis, George and Lockett, Stephen},
  journal = {Biophysical Journal},
  volume  = {86},
  number  = {6},
  pages   = {3993--4003},
  year    = {2004},
  doi     = {10.1529/biophysj.103.038422}
}

@article{manders1993,
  title   = {Measurement of co-localization of objects in dual-colour confocal images},
  author  = {Manders, E. M. M. and Verbeek, F. J. and Aten, J. A.},
  journal = {Journal of Microscopy},
  volume  = {169},
  number  = {3},
  pages   = {375--382},
  year    = {1993},
  doi     = {10.1111/j.1365-2818.1993.tb03313.x}
}

@article{talts2018,
  title   = {Validating Bayesian inference algorithms with simulation-based calibration},
  author  = {Talts, Sean and Betancourt, Michael and Simpson, Daniel and Vehtari, Aki and Gelman, Andrew},
  journal = {arXiv preprint arXiv:1804.06788},
  year    = {2018},
  doi     = {10.48550/arXiv.1804.06788}
}

@article{sainsburydale2024,
  title   = {Likelihood-Free Parameter Estimation with Neural Bayes Estimators},
  author  = {Sainsbury-Dale, Matthew and Zammit-Mangion, Andrew and Huser, Rapha{\"e}l},
  journal = {The American Statistician},
  volume  = {78},
  number  = {1},
  pages   = {1--14},
  year    = {2024},
  doi     = {10.1080/00031305.2023.2249522}
}

@article{radev2020,
  title   = {BayesFlow: Learning complex stochastic models with invertible neural networks},
  author  = {Radev, Stefan T. and Mertens, Ulf K. and Voss, Andreas and Ardizzone, Lynton and K{\"o}the, Ullrich},
  journal = {IEEE Transactions on Neural Networks and Learning Systems},
  volume  = {33},
  number  = {4},
  pages   = {1452--1466},
  year    = {2022},
  doi     = {10.1109/TNNLS.2020.3042395}
}

@article{radev2023joss,
  title   = {BayesFlow: Amortized Bayesian Workflows With Neural Networks},
  author  = {Radev, Stefan T. and Schmitt, Marvin and Schumacher, Lukas and Elsem{\"u}ller, Lasse and Pratz, Valentin and Sch{\"a}lte, Yannik and K{\"o}the, Ullrich and B{\"u}rkner, Paul-Christian},
  journal = {Journal of Open Source Software},
  year    = {2023},
  doi     = {10.21105/joss.05702}
}

@misc{seefelder2023proteincoloc,
  title        = {{ProteinCoLoc.jl}: A Julia package for co-localization analysis},
  author       = {Seefelder, Manuel},
  year         = {2023},
  howpublished = {\url{https://github.com/ma-seefelder/ProteinCoLoc}}
}

@article{seefelder_scirep,
  title   = {TODO: exact title of the ProteinCoLoc Scientific Reports paper},
  author  = {Seefelder, Manuel and others},
  journal = {Scientific Reports},
  year    = {TODO},
  doi     = {TODO: Scientific Reports DOI not present in repo; fill before submission}
}
```

**Provenance:** DOIs for `ordabayev2022`, `costes2004`, `manders1993`, `sainsburydale2024`, `radev2020`, `radev2023joss` are [CITED] from publisher/DOI pages surfaced in web search. Author lists are [ASSUMED] where not opened one-by-one — mark TODO if a full-author-list verification is desired. The Scientific Reports self-citation DOI is **not in the repo** (`CITATION.bib` only has the GitHub `@misc` entry) — left as an explicit TODO placeholder (allowed by D-11).

## Honest Headline Numbers (must appear caveated in the claim table + figure specs)

From STATE.md and Phase 4 Plan-05 summary [VERIFIED: `04-05-SUMMARY.md`, measured on committed fixtures, `BENCH_THREADS=1`, CPU-only]:

| Number | Honest form | Do NOT state as |
|--------|-------------|-----------------|
| Speedup | **Median ~325× forward-pass** (min pair ~125×); **~16–17× on the full posterior-sample (N=2000) workload**; against a **realized ADVI baseline ~0.5 s/pair** (seconds, not minutes) | bare ">100× vs minutes" |
| RMSE parity | NPE ρ_true RMSE **0.1271 ≈ ADVI 0.1270** (ratio ≈ 1.00 ≤ 1.2×) — speedup is only meaningful *paired* with this | speedup alone |
| Interval width | NPE 90% width **0.884 > ADVI 0.605** (mean-field VI under-disperses) — report side-by-side | "NPE intervals match ADVI" |
| Calibration | **"calibrated under the simulator"** (SBC pending Phase 5); always paired with the OOD result | "calibrated" unqualified |
| OOD | flag is **structurally bounded by the fixed patch-correlation summary**; summary-orthogonal misspecifications are **provably undetectable** and must be **named, not hidden** | "detects all misspecification" |

The claim-table row for speedup and the F1 (speedup-vs-RMSE) figure spec must reference the caveated form. The skeleton must reserve **explicit slots** for the "calibrated under the simulator" caveat and the "summary-orthogonal OOD blind spot" so later phases cannot quietly drop them.

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|------------------|--------|
| LaTeX for scientific typesetting | **Typst** single-binary, native BibTeX + programmable tables + fast incremental compile | No LaTeX toolchain; compile gate is a single exit-code check |
| Per-dataset Bayesian coloc inference (ADVI/SVI, incl. Tapqir, v1) | **Amortized SBI** (train once, ms per dataset) — the v2.0 thesis | The primary defensible related-work delta |
| Colocalization significance via randomization p-value (Costes) or co-occurrence (Manders) | Calibrated posterior + amortized Bayes factor + SBC coverage proof | The positioning axes of the claim table |

**Deprecated/outdated:** none relevant — this is a greenfield document.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Tapqir does **not** run formal SBC | Related Work matrix (⚠) | Overclaimed calibration delta; reviewer rejection — verify vs eLife paper |
| A2 | Tapqir has **no explicit Bayes-factor / model-comparison** output | Related Work matrix (⚠) | Weakens/miscasts a matrix cell — verify |
| A3 | Tapqir relies on **external drift correction**, not registration-UQ as a propagated latent | Related Work matrix + axis 3 (⚠) | Registration-UQ delta may be softer than stated — verify |
| A4 | Author lists / page numbers in the seeded `refs.bib` entries are complete/correct | Citation Entries | Wrong metadata in refs; low risk (TODO-marked, fixable pre-submission) |
| A5 | `ieee` (or `apa`) built-in bib style is acceptable for a venue-neutral draft | Standard Stack | Cosmetic; one-line swap later |
| A6 | Scientific Reports self-cite DOI must be filled later (not in repo) | Citation Entries | Missing citation until sourced; TODO-marked, non-blocking for compile |
| A7 | No existing `.typ`/`.tex` manuscript pipeline sits under a path that `manuscript/` could collide with | Decoupling | Low — `find` at depth 3 found none; `manuscript/` is a fresh top-level tree |

## Open Questions

1. **Does Tapqir perform formal SBC and/or emit an explicit Bayes factor?** (RESOLVED — Plan 10-03 Task 1 verifies the three ⚠ cells against Ordabayev et al. 2022 / DOI 10.7554/eLife.73860 and records a per-cell verdict before locking prose.)
   - What we know: physics-based Bayesian generative model, per-spot probabilities, per-dataset SVI (eLife abstract).
   - What's unclear: formal SBC rank-uniformity? explicit BF? registration uncertainty propagated to the coloc call?
   - Recommendation: planner adds a task to read Ordabayev et al. 2022 (methods) and confirm the ⚠ matrix cells before locking `related_work.typ` prose (D-09).
2. **Claim data location: `claims.typ` array vs external `claims.json`?** (RESOLVED — plan set uses the `claims.typ` array-of-dicts idiom consistently across Plans 10-01/10-02, per Claude's Discretion under D-06; no non-Typst machine-read requirement exists.)
   - What we know: both compile (`#json()` is supported); array-of-dicts is simplest and diffable.
   - Recommendation: use `claims.typ` array-of-dicts unless a non-Typst downstream agent must machine-read claims, in which case `claims.json` + `#json()`. (Claude's Discretion per D-06.)
3. **Scientific Reports self-citation DOI** — not in the repo. (RESOLVED as allowed TODO placeholder per D-11 — Plan 10-01 seeds refs.bib with a clearly-marked TODO DOI; source before Phase 16 submission. Non-blocking for compile since the key is uncited until then.)

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Typst | Compile gate (DoD, D-12) | ✓ | 0.15.0 (`3ae52774`) | none needed |
| Bash / Git Bash | `build.sh` gate | ✓ (this environment) | — | `build.jl` (Julia) variant |
| Julia | `build.jl` gate variant (optional) | ✓ (repo uses Julia 1.12.6) | 1.12.6 | `build.sh` |

**Missing dependencies with no fallback:** none.
**Missing dependencies with fallback:** none — all required tooling is present.

## Validation Architecture

`nyquist_validation` is enabled (config: `true`). For a document phase, the **compile gate is the validation surface** and phase success is fully falsifiable.

### Test Framework
| Property | Value |
|----------|-------|
| Framework | `typst compile` exit-code gate (no unit-test framework; the compiler *is* the test) |
| Config file | `manuscript/main.typ` (entry) + `manuscript/build.sh` (or `build.jl`) |
| Quick run command | `typst compile manuscript/main.typ manuscript/main.pdf` |
| Full suite command | `bash manuscript/build.sh` (compile + optional label/section assertions) |

### Phase Requirements → Test Map
| Success Criterion | Behavior | Test Type | Automated Command | Exists? |
|-------------------|----------|-----------|-------------------|---------|
| SC1: skeleton compiles with sections + claim table | `typst compile` returns exit 0 and produces `main.pdf` rendering `<tab:claims>` | build/smoke | `typst compile manuscript/main.typ manuscript/main.pdf; echo $?` | ❌ Wave 0 |
| SC2: Tapqir differentiation across 4 axes drafted | `related_work.typ` contains the matrix + the four-axis paragraph | content check | `grep -q -i "amortiz\|SBC\|registration\|spatial" manuscript/sections/related_work.typ` (or `typst query` for labels) | ❌ Wave 0 |
| SC3: figure specs enumerate per-phase panels | `figures/specs.typ` (or `claims.typ`) lists the claim-implied figures with owning phase | content check | `grep -c "SPEC F" manuscript/figures/specs.typ` ≥ N | ❌ Wave 0 |
| Claim table diffable + rendered (D-06) | `claims.typ` array present; renders in PDF | build + source check | compile (above) + `test -f manuscript/claims.typ` | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** `typst compile manuscript/main.typ manuscript/main.pdf` (must stay exit 0 after every edit — cheap, sub-second incremental).
- **Per wave merge:** `bash manuscript/build.sh` (compile + content assertions).
- **Phase gate:** `typst compile` green with claim table + related-work matrix rendering in the PDF (D-12).

### Wave 0 Gaps
- [ ] `manuscript/main.typ` — entry with set-rules + includes + `#bibliography` (SC1)
- [ ] `manuscript/claims.typ` — claim array + `claim_table()` (D-06)
- [ ] `manuscript/refs.bib` — seeded citations (D-11)
- [ ] `manuscript/sections/*.typ` — 7–9 stub files with `// TODO(phase N)` markers (D-05)
- [ ] `manuscript/sections/related_work.typ` — comparison matrix + Tapqir 4-axis prose (D-08/D-09)
- [ ] `manuscript/figures/specs.typ` — figure spec entries (D-10)
- [ ] `manuscript/build.sh` (and/or `build.jl`) — exit-code-correct compile gate (D-12)
- [ ] Framework install: none — Typst 0.15.0 already present

## Security Domain

`security_enforcement` is not set in config (treated as enabled). **This phase has no meaningful security surface**: it produces static `.typ`/`.bib`/`.pdf` documents and a local compile script. No authentication, session, access control, cryptography, network I/O, or untrusted input processing occurs.

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | — (no auth surface) |
| V3 Session Management | no | — |
| V4 Access Control | no | — |
| V5 Input Validation | no | — (inputs are author-written `.typ`/`.bib` in-repo) |
| V6 Cryptography | no | — |

**Only residual consideration:** the compile gate script must not pipe away exit codes (Pitfall 1) — a correctness/reliability concern, not a security one. No threat patterns apply to a static document build.

## Sources

### Primary (HIGH confidence)
- Installed **Typst 0.15.0** binary — all syntax/behavior claims (`#include`/`#import`, `table()`/`table.header`, `#figure`+`<label>`, `@cite`, `#bibliography("refs.bib", style:"ieee")`, exit codes for clean/undefined-cite/missing-file/uncited-bib, `typst query` availability) verified by compiling test files this session.
- eLife — Tapqir: https://elifesciences.org/articles/73860 (DOI 10.7554/eLife.73860)
- Biophysical Journal — Costes 2004: DOI 10.1529/biophysj.103.038422 (PubMed 15189895)
- Journal of Microscopy — Manders 1993: https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2818.1993.tb03313.x
- The American Statistician — Neural Bayes Estimators: https://www.tandfonline.com/doi/full/10.1080/00031305.2023.2249522
- IEEE TNNLS / JOSS — BayesFlow: DOI 10.1109/TNNLS.2020.3042395; DOI 10.21105/joss.05702
- Repo artifacts: `.planning/phases/04-.../04-05-SUMMARY.md` (honest speedup/RMSE numbers), `STATE.md`, `ROADMAP.md`, `PROJECT.md`, `CITATION.bib`, `10-CONTEXT.md`

### Secondary (MEDIUM confidence)
- arXiv — Talts et al. SBC (arXiv:1804.06788); BayesFlow arXiv:2306.16015 / 2003.06281

### Tertiary (LOW confidence — flagged for verification)
- Tapqir internals (formal SBC? explicit BF? registration-UQ handling?) — inferred from the eLife abstract + domain knowledge; ⚠ cells in the comparison matrix must be verified against the full paper before locking prose.

## Metadata

**Confidence breakdown:**
- Typst authoring / compile gate: **HIGH** — every construct compiled against the installed 0.15.0 binary; exit codes measured directly.
- Citations / DOIs: **HIGH** for DOIs (publisher-verified), MEDIUM for author-list/page completeness (TODO-marked).
- Related-work deltas: **MEDIUM** — method identities verified; Tapqir per-axis comparison cells flagged ASSUMED pending paper verification.
- Honest numbers: **HIGH** — copied from the measured Phase-4 summary.

**Research date:** 2026-07-02
**Valid until:** ~2026-08-01 (Typst is fast-moving pre-1.0; re-verify syntax if the installed version changes. Citations are stable.)
