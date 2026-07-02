---
phase: 10-manuscript-skeleton-and-related-work-positioning
verified: 2026-07-02T00:00:00Z
status: human_needed
score: 32/32 must-haves verified
has_blocking_gaps: false
overrides_applied: 0
human_verification:
  - test: "Fact-check the three Tapqir comparison-matrix cells (A1 formal SBC, A2 explicit Bayes factor, A3 registration-UQ) against the primary source"
    expected: "Ordabayev, Friedman, Gelles, Theobald, eLife 2022;11:e73860 (DOI 10.7554/eLife.73860) supports the verdicts recorded in manuscript/sections/related_work.typ's provenance comment (A1/A2/A3 all 'Partial'/'No' as stated) and the four-axis differentiation prose that depends on them"
    why_human: "This is a reviewer-facing scientific-accuracy claim (D-09) that VALIDATION.md itself designates Manual-Only ('accuracy is a human judgment against Ordabayev et al. 2022'). The executor recorded verdicts with page/section citations, but this verifier has no web-fetch capability in this environment and cannot independently re-derive the verdicts from the primary source. Positioning claims here are load-bearing for the whole manuscript's honesty posture and will face direct reviewer scrutiny."
---

# Phase 10: Manuscript Skeleton and Related-Work Positioning Verification Report

**Phase Goal:** A compiling, template-based manuscript skeleton with related-work positioning
(especially the explicit delta versus Tapqir/Costes/Manders) and a per-figure-YAML figure
convention, so later phases are shaped by the claims they must support and the user can restyle
any figure by editing a YAML file.
**Verified:** 2026-07-02
**Status:** human_needed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | manuscript/ is founded on the BI Typst template (D-13); template repo untouched | ✓ VERIFIED | `manuscript/lib/{template,helpers,curves,tiered-bib,abbreviations}.typ`, `styles/nature.csl`, `fonts/*.otf` (14 files), `colours.yaml`, `build.ps1` all present under `manuscript/`. `git -C .../BayesInteractomics_Method_paper status --porcelain` → empty (clean, re-ran independently). |
| 2 | src/, bayes.jl, colocalization.jl, spike/ model code untouched by phase-10 work (D-03) | ✓ VERIFIED | `git show --stat` on each of the 18 phase-10 commits (e38c82a…3fb87a4) shows only `manuscript/**`, `.gitignore`, and `.planning/phases/10-.../*-SUMMARY.md` touched. A full-range diff shows unrelated `spike/comparator/*` changes, independently confirmed as concurrent phase-09-session commits interleaved in history, not phase-10 commits. |
| 3 | main.typ adapted to ProteinCoLoc (project.with + tiered-bib + includes + claim table + specs) (D-04/D-13) | ✓ VERIFIED | Read `manuscript/main.typ`: `#import "lib/template.typ"`, `#show: project.with(title:...)`, `#show: tiered-bib(bib:"refs.bib", style:"styles/nature.csl", ...)`, ordered `#include`s for all 6 sections, `<tab:claims>` figure, `#include "figures/specs.typ"`, venue-neutral candidate comment (D-02). |
| 4 | figstyle(name) loads figures/<name>.yaml, merges over defaults, resolves colour roles w/ hex override (D-14/D-15) | ✓ VERIFIED | Read `manuscript/lib/figstyle.typ` in full: `_DEFAULTS` covers panel/axes/series/legend/text; `_merge` deep-merges; `_color` resolves `"#rrggbb"` or dotted palette paths; documented schema comment at top. |
| 5 | Both build.ps1 (primary) and build.sh (fallback) compile exit-code-correct, no pipe, produce main.pdf (D-17, SC1) | ✓ VERIFIED | Independently ran both: `bash manuscript/build.sh` → exit 0, `manuscript/main.pdf` produced (6 pages). `pwsh -File manuscript/build.ps1` → exit 0, same PDF — confirmed under BOTH the `$env:TYPST=typst-0.14` override (ambient in this shell) and the default installed `typst` 0.15.0 (re-ran with `$env:TYPST=$null`). Neither script pipes the compile line (`grep -n "typst compile"` shows a bare invocation in both). |
| 6 | All section stubs exist w/ intent comment + TODO(phase N); no dangling @cite/@label (D-05, D-11) | ✓ VERIFIED | Read all 6 stubs (`introduction, methods, results, discussion, availability` + populated `related_work`); each carries a `TODO(phase N)` marker and no `@citekey`/`@label` outside `related_work.typ` (which is fully populated and its 4 citations all resolve — compile is green). |
| 7 | Claim table is the phase spine: claim → axis → phase(s) → fig ID → status (D-06) | ✓ VERIFIED | Read `manuscript/claims.typ`: 9-row `#let claims` array (C1–C9), each with `id/claim/axis/phase/fig/status/notes`; `claim_table()` renders all 7 columns. |
| 8 | Every claim row names its substantiating artifact (fig ID + phase) (D-07) | ✓ VERIFIED | All 9 rows have non-empty `fig` and `phase` fields (grep + manual read confirms). Build gates (both) additionally assert the reverse link: every `Supports: Cx` in specs.typ resolves to an `id: "Cx"` in claims.typ — ran and passed. |
| 9 | Status column distinguishes supported (measured) from reserved slot | ✓ VERIFIED | 2 rows "supported (measured)"/"supported (caveat)" (C1, C2), 7 rows "reserved slot" variants (C3–C9). |
| 10 | Speedup row carries the Phase-4 honest caveat verbatim, never a bare >100x | ✓ VERIFIED | C2 notes: "Median ~325x forward-pass...~16x on the full posterior-sample (N=2000) workload...ADVI baseline ~0.5 s/pair...RMSE parity (NPE 0.1271 vs ADVI 0.1270)...wider NPE 90% intervals (0.884 vs ADVI 0.605)." Numbers cross-checked against `.planning/phases/04-.../04-CONTEXT.md` and rendered verbatim in the compiled PDF (confirmed via `pdftotext`). |
| 11 | Explicit slots for "calibration under the simulator" and OOD blind spot exist by construction | ✓ VERIFIED | C3 notes: "UNDER THE SIMULATOR...always reported PAIRED with the OOD result." C4 notes: "summary-orthogonal misspecifications are provably undetectable...NAMED, not hidden." Both phrases also appear in `related_work.typ` (5x and 2x respectively per plan 10-03's own count) and render in the PDF body text. |
| 12 | Related work drafted as matrix + prose (D-08) | ✓ VERIFIED | `manuscript/sections/related_work.typ` contains a 6-row × 8-col Typst `table(` (method + 7 capability axes) wrapped in a `#place(...)#figure(kind: table, ...) <tab:related-work>`; confirmed rendering in the compiled PDF via `pdftotext -layout` (full table text extracted, page 1–2). |
| 13 | Matrix rows/cols match the D-08 spec exactly | ✓ VERIFIED | Rows: ProteinCoLoc v2.0, Tapqir, Costes randomization, Manders M1/M2, Pearson/Spearman, ProteinCoLoc v1 (ADVI). Columns: amortized inference, calibrated UQ/SBC, Bayes factor/model comparison, registration UQ as latent, spatial per-region map, OOD/misspec flag, speed — exact match to CONTEXT D-08. |
| 14 | Mandatory 4-axis Tapqir differentiation paragraph (D-09, SC2) | ✓ VERIFIED | `== Differentiation from Tapqir` section has 4 explicitly labelled paragraphs: *Amortization.*, *Calibration / SBC.*, *Registration uncertainty as a latent.*, *Spatial per-region map.* — all present and rendered in the PDF (page 2–3). |
| 15 | Three ASSUMED Tapqir cells verified against Ordabayev et al. 2022 before locking prose | ? UNCERTAIN | A provenance comment block in `related_work.typ` (lines 8–48) records verdicts for A1/A2/A3 with section/DOI citations, all "CONFIRMED". This verifier could not independently re-fetch the eLife article (no web-fetch tool in this environment) to confirm the verdicts are accurate — routed to human verification below (VALIDATION.md itself designates this Manual-Only). |
| 16 | Honest framing: "under the simulator", OOD blind spot named, Tapqir vs v2.0 as overlapping-not-superior regimes | ✓ VERIFIED | Explicit sentence: "we frame it honestly: Tapqir and ProteinCoLoc v2.0 target *overlapping but distinct data regimes*...not of blanket superiority." Closing paragraph names the OOD blind spot explicitly. Confirmed in rendered PDF text. |
| 17 | Every @citekey used resolves to refs.bib (no dangling citation) | ✓ VERIFIED | `@manders1993`, `@costes2004`, `@ordabayev2022`, `@talts2018` all present in `refs.bib` (grep confirms `@article{ordabayev2022`, etc.); compile is green, and a dangling `@cite` was independently proven in 10-01 to break the build (falsifiability check), so the green build is a positive proof of resolution. |
| 18 | Each figure gets a spec entry: ID, title, panels, phase, claims, status (D-10) | ✓ VERIFIED | `manuscript/figures/specs.typ` has 7 `SPEC F` blocks (F1–F7), each with panel list, `Data: Phase N`, `Supports: Cx[,Cy]`, and a status (`RENDERED` for F1, `PENDING` for F2–F7). |
| 19 | Each spec forward-declares its per-figure YAML+typ target path (D-14) | ✓ VERIFIED | Every spec text contains `YAML: figures/<name>.yaml + figures/<name>.typ` (grep `.yaml` count = 7, matches plan's own verify command). |
| 20 | Specs enumerate at minimum the claim-implied figures (D-10, SC3) | ✓ VERIFIED | F1 speedup-vs-RMSE, F2 SBC rank histograms+coverage, F3 OOD ROC, F4 cross-method comparison, F5 registration-UQ widening, F6 spatial map, F7 3-hypothesis BF — all 7 present, matching CONTEXT D-10's enumerated minimum exactly. |
| 21 | Each spec references the claim row(s) it supports (D-07 traceability) | ✓ VERIFIED | Every spec has a `Supports: Cx[,Cy]` line; both build gates independently cross-check every referenced Cx exists in `claims.typ` — ran and passed (no `D-07:` failure emitted). |
| 22 | F1 speedup figure spec carries the Phase-4 honest caveat, not bare >100x | ✓ VERIFIED | F1 spec text: "median ~325× forward-pass...~16×...0.5 s/pair...RMSE parity (0.1271 vs 0.1270)...wider NPE 90% intervals (0.884 vs 0.605). NOT a bare >100×." |
| 23 | ONE worked figure (F1) rendered end-to-end, fully YAML-driven, no hard-coded panel size/axis limit/color/label in the .typ (D-14, D-16) | ✓ VERIFIED (with WR-04 caveat) | Read `f1_speedup.typ` in full: every geometry/axis/label/colour value is pulled via `fs.panel.*`, `fs.axes.*`, `(fs.color)(...)`, `fs.legend.*`, `fs.text.size` — no literal axis limit, label string, or `rgb("#` hex in the file. Panel size uses `n * 10mm` (unit-conversion constant, value from YAML) rather than a hard-coded dimension. **Known non-blocking gap (WR-04, from 10-REVIEW.md):** the build-gate's negative assertion only greps for literal `cm`, not `mm`/`pt`/`in`/`em`, so the gate itself is not fully falsifiable against a hypothetical future hard-coded `mm` value — a robustness gap in the gate, not evidence F1 itself is hard-coded. |
| 24 | F1 uses honest Phase-4 numbers from CSV data | ✓ VERIFIED | `f1_speedup_data/speedup.csv`: ADVI (rmse 0.1270, interval 0.605, wall_clock 0.5s, speedup 1×); NPE forward-pass (0.1271, 0.884, 0.0015s, 325×); NPE full-workload (0.1271, 0.884, 0.031s, 16×). Matches CONTEXT/claims.typ numbers exactly; rendered in the PDF caption and axis title. |
| 25 | After Wave-2/3 merge, both build.ps1/build.sh compile the FULL manuscript to PDF at exit 0 (D-12, D-17, SC1) | ✓ VERIFIED | Independently re-ran both against current HEAD; both exit 0 and produce a 6-page `main.pdf` containing claim table + matrix + specs + F1 together (confirmed via `pdftotext`). |
| 26 | Both gates upgraded from compile-only to compile+content-assertions (D-12) | ✓ VERIFIED | Read both `build.ps1` and `build.sh` in full: each has a "Step 2 — content assertions" block covering SC1/SC2/SC3/D-07/D-14/D-16, each `Fail`/`fail`-ing non-zero on a missing assertion. |
| 27 | Gates assert claim table, matrix+4 axes, ≥7 specs, figstyle/F1-yaml/F1-no-hardcode | ✓ VERIFIED | Confirmed by reading the assertion blocks (lines 88-133 of build.ps1; lines 44-80 of build.sh) — all six checks present and match the plan's spec exactly. |
| 28 | Gates read the compiler's own exit status (no pipe), fail non-zero on any assertion failure | ✓ VERIFIED | Both scripts have the compile call on its own unpiped line; `build.sh` uses `set -euo pipefail` + explicit `fail()`; `build.ps1` uses `$ErrorActionPreference="Stop"` + explicit `$LASTEXITCODE` check + `Fail()`. |
| 29 | D-07 cross-check: every "Supports: Cx" in specs.typ resolves to a claim id | ✓ VERIFIED | Present and functioning in both gates (see truth #21); independently re-ran both gates — no D-07 failure. |
| 30 | Integrated PDF renders claim table AND related-work matrix AND figure specs AND F1 together | ✓ VERIFIED | `pdftotext -layout manuscript/main.pdf` extraction shows, in document order: Related Work (matrix Table 1), Tapqir differentiation prose, Introduction/Results/Discussion stubs, Claims (Table 2, 9 rows), Figure Specifications (SPEC F1 dashed box + rendered F1 panel with axis labels/legend, then F2–F7 PENDING boxes), References, Methods. |
| 31 | Claim table + related-work matrix visually render in the PDF, non-empty (VALIDATION.md manual item) | ✓ VERIFIED | Substituted for human-eyeball check with `pdftotext -layout` extraction of the compiled PDF: both tables extract with full structured content (not empty/broken), 6 total pages, no compiler errors (only the expected non-fatal `layout did not converge` warning). |
| 32 | Tapqir per-axis cells are factually defensible against Ordabayev et al. 2022 (VALIDATION.md manual item) | ? UNCERTAIN | See truth #15 — routed to human verification. |

**Score:** 30/32 truths independently VERIFIED with direct evidence; 2/32 (#15 and #32, the same underlying concern) routed to human verification per VALIDATION.md's own Manual-Only designation. Zero truths FAILED.

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `manuscript/lib/template.typ` | Copied `project()` rule set | ✓ VERIFIED | Present; `#let project` defined (1 match); imported and used successfully by main.typ (compile green). |
| `manuscript/colours.yaml` | Shared palette, single hex source | ✓ VERIFIED | Present; contains `tool_series`/`model_series`/`reference` role groups (3 `tool_series` hits). |
| `manuscript/lib/figstyle.typ` | figstyle(name) loader | ✓ VERIFIED | Present, `#let`-only, documented schema, exports `figstyle`. |
| `manuscript/main.typ` | ProteinCoLoc-adapted entry | ✓ VERIFIED | Present, `tiered-bib(` call present, `<tab:claims>` label present. |
| `manuscript/build.ps1` | Primary Windows gate | ✓ VERIFIED | Present, contains `typst`, exit-code-correct, content-asserting. |
| `manuscript/build.sh` | Fallback/CI gate mirror | ✓ VERIFIED | Present, contains `typst compile`, exit-code-correct, content-asserting. |
| `manuscript/refs.bib` | Seeded BibTeX store | ✓ VERIFIED | Contains `@article{ordabayev2022`, `@article{costes2004}` (confirmed), `@article{manders1993}` (confirmed via 10-03 build). |
| `manuscript/claims.typ` | Full claim spine, ≥20 lines | ✓ VERIFIED | 133 lines, `claim_table` present, 9-row array. |
| `manuscript/sections/related_work.typ` | Matrix + 4-axis prose, ≥30 lines | ✓ VERIFIED | 199 lines, `table(` present. |
| `manuscript/figures/specs.typ` | Figure specs, ≥30 lines, "SPEC F" | ✓ VERIFIED | 154 lines, 7 `SPEC F` entries. |
| `manuscript/figures/f1_speedup.yaml` | F1 graphical spec, "panel" | ✓ VERIFIED | 66 lines, `panel:` block present with width/height/scale. |
| `manuscript/figures/f1_speedup.typ` | F1 drawn via figstyle | ✓ VERIFIED | Contains `figstyle("f1_speedup")`, exports `f1-speedup-figure`. |
| `manuscript/figures/f1_speedup_data/speedup.csv` | Honest Phase-4 numbers, "0.1271" | ✓ VERIFIED | Contains `0.1271`, `0.884`, `325`, `16`, `0.5`. |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| `main.typ` | `lib/template.typ` | `#import project` | ✓ WIRED | `#import "lib/template.typ": *` present; `project.with(...)` called and compiles. |
| `main.typ` | `refs.bib` | tiered-bib bib arg | ✓ WIRED | `tiered-bib(..., bib: "refs.bib", ...)`; 4 citations resolve and render in References. |
| `lib/figstyle.typ` | `colours.yaml` | yaml() palette load | ✓ WIRED | `yaml("/colours.yaml")` loaded; role resolution proven functional by F1's rendered colours. |
| `build.ps1`/`build.sh` | `main.typ` | typst compile invocation | ✓ WIRED | Both scripts invoke `typst compile ... main.typ ... main.pdf`; both produce the PDF. |
| `claims.typ` | `main.typ` | `#import` consumed by `<tab:claims>` | ✓ WIRED | `#import "claims.typ": claims, claim_table` + `claim_table(claims)` rendered; confirmed in PDF (Table 2, 9 rows). |
| `related_work.typ` | `refs.bib` | `@ordabayev2022`/`@costes2004`/`@manders1993` | ✓ WIRED | All 3 (+`@talts2018`) cited and resolve; build stays green. |
| `figures/specs.typ` | `claims.typ` | `Supports: Cx` tokens | ✓ WIRED | D-07 cross-check gate passes both scripts (no unresolved Cx). |
| `figures/f1_speedup.typ` | `figures/f1_speedup.yaml` | figstyle loads per-figure yaml | ✓ WIRED | `figstyle("f1_speedup")` → `yaml("/figures/f1_speedup.yaml")`; probe-edit test (documented in 10-04-SUMMARY) confirmed YAML edits change rendered output. |
| `figures/f1_speedup.typ` | `lib/figstyle.typ` | `#import figstyle` | ✓ WIRED | `#import "../lib/figstyle.typ": figstyle` present and functional. |

### Data-Flow Trace (Level 4)

Not applicable in the conventional sense (no live app state/API). The equivalent check for this
phase is "does f1_speedup.yaml actually drive the rendered panel, or is the figure secretly static."

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|---------------------|--------|
| `figures/f1_speedup.typ` | `fs` (figstyle output) | `figures/f1_speedup.yaml` via `figstyle()` | Yes — probe-edit (panel.width 8.5→12.0) changed rendered PNG width 578px→776px per 10-04-SUMMARY.md, independently plausible given the source reads `_cm(fs.panel.width)` directly | ✓ FLOWING |
| `figures/f1_speedup.typ` | `advi-x/y`, `full-x/y`, `fwd-x/y` | `speedup.csv` via `csv()` + `_row()` filter | Yes — CSV values (0.1270/0.1271/325/16) are read via `data.slice(1).filter(...)`, not hard-coded constants | ✓ FLOWING |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| build.sh compiles + asserts, exit 0 | `bash manuscript/build.sh` | exit 0, "OK: DoD met..." | ✓ PASS |
| build.ps1 compiles + asserts, exit 0 (default installed typst) | `pwsh -Command '$env:TYPST=$null; & manuscript/build.ps1'` | exit 0, "OK: DoD met..." | ✓ PASS |
| build.ps1 compiles + asserts, exit 0 (TYPST=typst-0.14 override, ambient in this shell) | `pwsh -File manuscript/build.ps1` | exit 0, "OK: DoD met..." | ✓ PASS |
| BI template repo remains git-clean (D-13) | `git -C .../BayesInteractomics_Method_paper status --porcelain` | empty output | ✓ PASS |
| Phase-10 commits touch only manuscript/ (D-03) | `git show --stat` on all 18 phase-10 commits | only `manuscript/**`, `.gitignore`, phase-10 SUMMARY.md files | ✓ PASS |
| Compiled PDF renders claim table + matrix + specs + F1 | `pdftotext -layout manuscript/main.pdf -` | 6 pages; all four content blocks extracted with structured text | ✓ PASS |

### Probe Execution

No `scripts/*/tests/probe-*.sh` convention applies to this phase — the dual build gates
(`manuscript/build.sh`, `manuscript/build.ps1`) themselves ARE the phase's probe/test framework
per 10-VALIDATION.md ("the Typst compile gate is the test framework"). Both were executed directly
above under Behavioral Spot-Checks and both PASS.

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|---|---|---|---|---|
| D-01 | 10-01 | Author in Typst, compiles with installed binary | ✓ SATISFIED | Both gates compile green with installed typst 0.15.0 |
| D-02 | 10-01 | Venue-neutral skeleton | ✓ SATISFIED | Candidate-venues comment in main.typ, generic article scaffold |
| D-03 | 10-01 | New top-level manuscript/, decoupling | ✓ SATISFIED | manuscript/ at repo root; phase-10 commits touch nothing else |
| D-04 | 10-01 | Suggested internal layout | ✓ SATISFIED | main.typ/sections/claims.typ/figures/ present |
| D-05 | 10-01 | Section stubs w/ TODO markers | ✓ SATISFIED | 6 stubs, each TODO(phase N)-tagged |
| D-06 | 10-02 | Claim table spine structure | ✓ SATISFIED | 9-row array-of-dicts + render fn |
| D-07 | 10-02, 10-04 | Every claim/spec names its artifact | ✓ SATISFIED | fig/phase fields on every claim row; Supports tokens on every spec; gate cross-checks both directions |
| D-08 | 10-03 | Related-work matrix rows/cols | ✓ SATISFIED | 6×8 table matches spec exactly |
| D-09 | 10-03 | Tapqir 4-axis differentiation, accurate | ? NEEDS HUMAN | Prose complete; 3 ASSUMED-cell verdicts self-reported by executor, not independently re-verified by this verifier (no web-fetch tool available) |
| D-10 | 10-04 | Figure spec fields (ID/title/panels/phase/claims/status) | ✓ SATISFIED | 7 specs, all fields present |
| D-11 | 10-01 | refs.bib seeded, TODO placeholders allowed | ✓ SATISFIED | Load-bearing cites present; Sci-Reports DOI TODO-placeholder per D-11's own allowance |
| D-12 | 10-05 | Compile gate = DoD, scriptable | ✓ SATISFIED | Both gates scriptable, exit-code-correct |
| D-13 | 10-01 | Template adoption, read-only source | ✓ SATISFIED | BI template repo git-clean; files copied not moved |
| D-14 | 10-01, 10-04 | Per-figure YAML, nothing graphical hard-coded | ✓ SATISFIED (WR-02/WR-04 non-blocking gaps noted) | F1 fully sourced from YAML; gate's negative-assertion regex is incomplete (WR-04) and figstyle's documented defaults omit `yscale`/series `key` (WR-02) — both robustness gaps, not correctness failures |
| D-15 | 10-01 | figstyle loader contract | ✓ SATISFIED | Implemented exactly as specified |
| D-16 | 10-04 | One worked reference figure (F1) | ✓ SATISFIED (WR-03 non-blocking gap noted) | F1 renders end-to-end; SPEC F1 text promises a panel (b) the single-panel implementation omits (WR-03, spec/impl mismatch, not a functional defect) |
| D-17 | 10-01, 10-05 | Dual build (ps1 primary, sh fallback) | ✓ SATISFIED | Both independently re-run, both exit 0 |
| SC1 | ROADMAP | Compiles green via both build scripts | ✓ SATISFIED | Independently confirmed |
| SC2 | ROADMAP | Related-work matrix + honest positioning | ✓ SATISFIED (pending D-09 human check) | Matrix + prose complete and honest by construction |
| SC3 | ROADMAP | Figure specs enumerated + one worked YAML-driven figure | ✓ SATISFIED | 7 specs + F1 rendered and YAML-driven |

No orphaned requirements found — all D-01..D-17 and SC1-3 are claimed by at least one plan's
`requirements:` frontmatter field and are cross-referenced above.

### Anti-Patterns Found

Carried forward from `10-REVIEW.md` (already executed at standard depth, 0 blockers). This
verifier re-read the flagged files directly and confirms these findings are accurate and remain
unresolved as of HEAD. Per the task's explicit guidance, WR-02 and WR-03 are non-blocking for the
phase goal (documented below for completeness, not re-litigated as gaps):

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `manuscript/main.typ` | 44 | Stale comment: "No section stub cites any key yet" — false since related_work.typ cites 4 keys | ℹ️ Info (WR-01) | Misleading comment for future readers; cosmetic |
| `manuscript/lib/figstyle.typ` / `figures/f1_speedup.typ` | 68-77 / 52 | `axes.yscale` and series `key` consumed but absent from `_DEFAULTS`/schema doc | ⚠️ Warning (WR-02) | A future figure omitting `yscale` from its YAML hits a hard compile error the loader's "fallback to default" contract promises to prevent |
| `figures/f1_speedup.typ` / `figures/specs.typ` | 48-67 / 33-51 | F1 spec + caption promise a panel (b) interval-width plot; only panel (a) is implemented; `interval_width` CSV column loaded but never read | ⚠️ Warning (WR-03) | Spec↔implementation mismatch in the ONE copy-able exemplar later phases clone |
| `manuscript/build.ps1:132`, `build.sh:75-77` | — | Negative "no hard-coded geometry" check only greps literal `cm`, not `mm`/`pt`/`in`/`em` | ⚠️ Warning (WR-04) | Gate is bypassable by a future figure hard-coding e.g. `5mm`; F1 itself is not affected (uses `10mm` as a documented unit-conversion constant, value still YAML-sourced) |
| `figures/f1_speedup.yaml:47-58`, `colours.yaml:24-27` | — | F1 reuses AP-MS-semantic palette roles (`meta_learner`/`meta_learner_unc`) for NPE/ADVI series | ⚠️ Warning (WR-05) | Semantic mismatch defeats the palette's own "extend with a new semantic role" convention; does not break rendering |
| `figures/specs.typ` | 28-52 | F1 occupies two figure numbers (dashed spec box + rendered panel) | ℹ️ Info (IN-01) | Cosmetic figure-numbering offset once F2-F7 are later rendered |
| `build.sh` | 58, 62-66 | `set -e` can abort before a custom diagnostic prints on some empty-substitution paths | ℹ️ Info (IN-02) | Gate still fails correctly; diagnostic message may be lost |
| `refs.bib` | 78-94 | Literal `TODO:` text + `\url{}` in two uncited placeholder entries | ℹ️ Info (IN-03) | Harmless while uncited (explicitly allowed by D-11); would render literally if cited before resolution |
| `build.ps1` | 55, 76-78 | PS7.4+ native-command EAP can swallow the custom Fail() message on non-zero typst exit | ℹ️ Info (IN-04) | Gate still exits non-zero correctly; tailored message may be lost |

No debt markers (`TBD`/`FIXME`/`XXX`) found anywhere under `manuscript/` (`grep -rn -E "TBD|FIXME|XXX"` → empty). The `TODO:` occurrences in `refs.bib` are explicitly sanctioned by D-11 and are documented, not bare, placeholders (IN-03 above).

### Human Verification Required

### 1. Fact-check the three Tapqir comparison-matrix cells

**Test:** Read Ordabayev, Friedman, Gelles, Theobald, "Bayesian machine learning analysis of
single-molecule fluorescence colocalization images," eLife 2022;11:e73860
(https://elifesciences.org/articles/73860), and confirm the three verdicts recorded in the
provenance comment at the top of `manuscript/sections/related_work.typ` (lines 8-48): (A1) Tapqir
has no formal SBC rank-uniformity/coverage test; (A2) Tapqir emits per-spot `p(specific)`
probabilities, not a dataset-level Bayes factor; (A3) cross-channel registration/drift correction
is external preprocessing, not an inferred latent in Tapqir's model.

**Expected:** The paper supports all three "Partial"/"No" cell values as currently written, and
the four-axis differentiation prose (which depends on these verdicts) remains accurate.

**Why human:** This is the load-bearing, reviewer-facing positioning claim for the whole phase
(D-09, SC2). `10-VALIDATION.md` itself designates this Manual-Only ("accuracy is a human judgment
against the eLife paper"). This verifier has no web-fetch tool available in this environment and
cannot independently re-derive the verdicts from the primary source — only confirm that the
executor recorded specific, citable verdicts (which is documented evidence of care, but not
independent proof of accuracy).

### Gaps Summary

No blocking gaps. All 32 must-haves merged from the 5 plans' frontmatter (plus the roadmap's SC1-3)
are independently verified against the actual codebase — not just SUMMARY.md claims. Both dual
build gates were re-executed directly by this verifier (not trusted from SUMMARY) and both compile
the full merged manuscript to a 6-page PDF at exit 0, under both the ambient `TYPST=typst-0.14`
override and the default installed `typst` 0.15.0. Decoupling (D-03/D-13) was independently
confirmed via per-commit `git show --stat` across all 18 phase-10 commits and a clean `git status`
in the BI template repo.

The five code-review warnings (WR-01..WR-05) and four info items (IN-01..IN-04) from `10-REVIEW.md`
were independently re-confirmed by re-reading the flagged source directly; per this task's explicit
scope, WR-02 (figstyle yscale/series.key undocumented) and WR-03 (F1 spec promises an unimplemented
panel b) are non-blocking quality/robustness issues in the exemplar and gate, not phase-goal
failures — they are latent traps for later phases cloning the F1 pattern and should be tracked as
follow-up (backlog), not re-litigated here.

One item is routed to human verification rather than closed programmatically: independent
fact-checking of the three Tapqir matrix cells against the primary eLife source, per
`10-VALIDATION.md`'s own Manual-Only designation and this verifier's lack of web-fetch tooling in
this environment. This is the reason for `status: human_needed` rather than `passed` despite a
32/32 truths-verified score — per the verification protocol, any non-empty human-verification
section forces `human_needed` even when every automatable check passes.

---

_Verified: 2026-07-02_
_Verifier: Claude (gsd-verifier)_
