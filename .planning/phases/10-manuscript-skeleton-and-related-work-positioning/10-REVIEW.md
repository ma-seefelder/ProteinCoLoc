---
phase: 10-manuscript-skeleton-and-related-work-positioning
reviewed: 2026-07-02T00:00:00Z
depth: standard
files_reviewed: 17
files_reviewed_list:
  - manuscript/main.typ
  - manuscript/claims.typ
  - manuscript/refs.bib
  - manuscript/colours.yaml
  - manuscript/lib/figstyle.typ
  - manuscript/sections/introduction.typ
  - manuscript/sections/related_work.typ
  - manuscript/sections/methods.typ
  - manuscript/sections/results.typ
  - manuscript/sections/discussion.typ
  - manuscript/sections/availability.typ
  - manuscript/figures/specs.typ
  - manuscript/figures/f1_speedup.typ
  - manuscript/figures/f1_speedup.yaml
  - manuscript/figures/f1_speedup_data/speedup.csv
  - manuscript/build.ps1
  - manuscript/build.sh
findings:
  critical: 0
  warning: 5
  info: 4
  total: 9
status: issues_found
---

# Phase 10: Code Review Report

**Reviewed:** 2026-07-02
**Depth:** standard
**Files Reviewed:** 17
**Status:** issues_found

## Summary

This is a compile-safe Typst/YAML/BibTeX manuscript skeleton plus two regression build gates
(`build.ps1` primary Windows, `build.sh` CI fallback). Reviewed against the five focus areas:
build-gate correctness, the `figstyle` per-figure YAML loader contract, honesty of the F1
numbers, dangling Typst `@cite`/`@label` targets, and cross-file consistency between
`claims.typ`, `specs.typ`, and `related_work.typ`.

**No blockers.** The skeleton is structurally sound for its stated DoD:

- **No dangling citations or references.** Exactly four keys are cited (`@manders1993`,
  `@costes2004`, `@ordabayev2022`, `@talts2018`) — all present in `refs.bib`. The one
  cross-reference `@tab:related-work` resolves to its `<tab:related-work>` label. Spec/figure
  labels (`<fig:f1_speedup>`, `<fig:f2_sbc>`, …) are defined but not `@`-referenced (correct
  per Pitfall 2).
- **F1 honesty holds.** No bare `>100×`. The 325×/16×/0.5 s/RMSE-parity (0.1271 vs 0.1270)/
  interval-width (0.884 vs 0.605) numbers are consistent across `claims.typ` (C2),
  `specs.typ`, `f1_speedup.yaml`, and `speedup.csv`, and reconcile to 2 sig figs
  (0.0015 s × 325 ≈ 0.5 s; 0.031 s × 16 ≈ 0.5 s).
- **`figstyle` deep-merge is correct.** Nested dicts merge recursively; the `series` array
  replaces the default wholesale; role → hex resolution and `#rrggbb` override both work; no
  aliasing/mutation leak into `base`.
- **Exit codes are preserved.** Both gates keep the `typst` compile line pipe-free; warnings
  are tolerated (exit-code driven, not stderr-grep driven).

The findings below are quality/robustness issues — latent traps for the later phases that
clone this skeleton, plus a stale comment and a spec-vs-implementation mismatch in the worked
F1 exemplar.

## Warnings

### WR-01: `main.typ` bibliography comment is factually wrong

**File:** `manuscript/main.typ:44`
**Issue:** The comment asserts *"No section stub cites any key yet, so the lists render empty
until later phases add @cites."* This is false: `sections/related_work.typ` already cites four
keys (`@manders1993`, `@costes2004`, `@ordabayev2022`, `@talts2018`), all at the default H
tier (they precede the `#set-tier("M")` in `main.typ:69`). The H-tier `print-references`
call therefore renders a non-empty "References" list. A load-bearing comment about compile/bib
behaviour that no longer matches the code will mislead the Phase-16 author reasoning about the
tiered bibliography.
**Fix:** Update the comment, e.g.:
```typst
// Related Work already cites 4 keys (H tier); the Methods list stays empty until Phase 5/7
// wire @cites under the M tier. Each key cited must exist in refs.bib (a dangling @cite fails).
```

### WR-02: `figstyle` schema/consumer contract gap — `axes.yscale` has no default

**File:** `manuscript/lib/figstyle.typ:68-77` (defaults) and `manuscript/figures/f1_speedup.typ:52`
**Issue:** `f1_speedup.typ` reads `fs.axes.yscale`, but `yscale` is absent from
`_DEFAULTS.axes` **and** from the documented schema comment (`figstyle.typ:34-42`). F1 only
compiles because `f1_speedup.yaml:23` happens to set `yscale: "log"`. Any later figure that
clones the F1 pattern but omits `yscale` from its YAML will hit a hard compile error
(`dictionary does not contain key "yscale"`) — the exact failure the loader's "any key you
omit falls back to the default" contract promises to prevent. The same undocumented-but-required
pattern applies to `series[].key`, which the consumer relies on (`f1_speedup.typ:28`,
`_ser.insert(s.key, s)`) but which the schema (`figstyle.typ:43-51`) never lists.
**Fix:** Add the key to defaults and document both:
```typst
axes: (
  xlim: none, ylim: none, xticks: none, yticks: none,
  xlabel: none, ylabel: none, title: none,
  yscale: "linear", xscale: "linear",   // <- add; consumed by figure .typ
),
```
and note `key:` as a required per-series field in the schema comment.

### WR-03: Worked F1 exemplar renders only panel (a); spec/caption promise two panels

**File:** `manuscript/figures/f1_speedup.typ:48-67`, `manuscript/figures/specs.typ:33-51`,
`manuscript/figures/f1_speedup_data/speedup.csv`
**Issue:** SPEC F1 (`specs.typ:33-34`) and the rendered-figure caption (`specs.typ:46-51`)
both describe *"(b) paired 90% interval width, NPE vs ADVI"* and cite the honesty-critical
interval comparison (0.884 vs 0.605). The figure body builds a **single** `lq.diagram`
(speedup-vs-RMSE only) — there is no panel (b). Consequently the `interval_width` column in
`speedup.csv` is loaded by `csv(...)` but never read by `_row` (which extracts only cols 2/5),
and the "wider NPE intervals" point survives only as caption prose, not as a rendered panel.
For the ONE worked reference exemplar (D-16) that later phases clone, the spec↔implementation
mismatch is worth closing now.
**Fix:** Either render panel (b) (a paired interval-width bar/point using the existing
`interval_width` column), or reconcile SPEC F1 + the caption to a single-panel figure so the
exemplar and its contract agree.

### WR-04: Build-gate "no hard-coded geometry" check is trivially bypassable

**File:** `manuscript/build.ps1:132`, `manuscript/build.sh:75-77`
**Issue:** The D-14 negative assertion greps only for the literal unit `cm`
(`[0-9]+(\.[0-9]+)?cm`). `f1_speedup.typ` deliberately uses `n * 10mm` (`f1_speedup.typ:22`)
to convert the YAML cm value, which passes. But the check does **not** catch hard-coded
`mm`, `pt`, `in`, or `em` sizes — a future figure that hard-codes `width: 5mm` would sail
through the gate while violating the very rule (nothing graphical hard-coded in the `.typ`)
the gate exists to enforce. Weak falsifiability for a regression gate.
**Fix:** Broaden the pattern to catch any numeric-literal length unit outside the sanctioned
conversion constant, e.g. flag `[0-9]+(\.[0-9]+)?(cm|mm|in|pt|em)` and whitelist the single
`* 10mm` / `* 1pt` conversion idioms, or assert positively that geometry expressions reference
`fs.`.

### WR-05: F1 colours NPE/ADVI series with AP-MS-semantic palette roles

**File:** `manuscript/figures/f1_speedup.yaml:47-58`, `manuscript/colours.yaml:24-27`
**Issue:** F1 assigns the NPE regimes to `model_series.meta_learner` and
`model_series.meta_learner_unc`, whose documented meanings in the copied palette are
*"multi-evidence stack (tr_ddi)"* and *"stack + dropout uncertainty (tr_ddi_mc)"* — semantics
from the BayesInteractomics AP-MS project, not colocalization. `colours.yaml:11-14` states the
file's own convention: colours are named by **semantic role**, and one should *"extend this
file with a new semantic role"* rather than reuse an unrelated one. Mapping an NPE-forward-pass
curve onto a "meta-learner" role defeats that contract (an editor changing `meta_learner`
would not expect to recolour an NPE curve) and leaves the palette carrying ~200 lines of roles
irrelevant to this manuscript. Note the tension with D-13 (`colours.yaml` copied read-only):
if the palette truly cannot be extended, F1 should at least use role names whose semantics are
neutral (e.g. `sample_role`/`reference`) rather than an inference-method-mislabelled one.
**Fix (recommend):** Add colocalization-semantic roles
(e.g. `inference_method: { npe_forward, npe_full, advi }`) to `colours.yaml` and point
`f1_speedup.yaml` at them; if D-13 forbids editing the copied palette, document the exception
and pick semantically neutral roles.

## Info

### IN-01: F1 consumes two figure numbers (dashed spec box + rendered panel)

**File:** `manuscript/figures/specs.typ:28-52`
**Issue:** `<fig:f1_speedup_spec>` (the dashed SPEC box) and `<fig:f1_speedup>` (the rendered
panel) are both `#figure(...)` with captions, so F1 occupies two consecutive figure numbers —
the rendered F1 becomes "Figure 2". F2–F7 each add one more spec-box figure. When later phases
swap real panels in for F2–F7, the spec-box duplicates linger unless removed, and the
F-number ↔ figure-number mapping stays off by the F1 spec box.
**Fix (CONVENTION):** Drop the F1 spec box now that its panel is rendered, or render the spec
boxes as non-`#figure` placeholders (plain `rect`) so they do not consume figure numbers.

### IN-02: `build.sh` `set -e` suppresses the custom SC3/D-07 diagnostics on failure

**File:** `manuscript/build.sh:25,58,62-66`
**Issue:** Under `set -euo pipefail`, `SPEC_COUNT=$(… | grep -cE …)` aborts the script when
the count is 0 (`grep -c` exits 1) *before* the line-59 `[ "$SPEC_COUNT" -ge 7 ]` check runs,
so the intended `fail "SC3: only $SPEC_COUNT figure specs found"` message never prints. Same
class of issue for empty command substitutions in the D-07 loop. The gate still fails
(correctly), but with a bare `set -e` abort rather than the diagnostic. `build.ps1` reports the
precise message here.
**Fix (CONVENTION):** Capture with `|| true` and test explicitly, e.g.
`SPEC_COUNT=$(noc … | grep -cE "SPEC F[0-9]" || true)`.

### IN-03: `refs.bib` placeholder entries carry literal `TODO:` / LaTeX `\url{}`

**File:** `manuscript/refs.bib:78-94`
**Issue:** `@seefelder_scirep` has `title = {TODO: …}`, `doi = {TODO: …}`, `note = {TODO: …}`;
`@seefelder2023proteincoloc` uses `howpublished = {\url{…}}`. Both are currently uncited so
neither renders (harmless, and D-11 explicitly allows the TODO placeholder). But if either is
cited before the TODO is resolved, the literal `TODO:` text and a possibly-unparsed `\url{}`
(hayagriva does not interpret LaTeX macros) would render into the reference list.
**Fix (CONVENTION):** Resolve the Scientific Reports DOI/title before Phase 16 wires the
self-citations; replace `\url{}` with a plain `url = {…}` field.

### IN-04: `build.ps1` non-zero-exit path skips the custom `Fail` message

**File:** `manuscript/build.ps1:55,76-78`
**Issue:** With `$ErrorActionPreference = "Stop"` and PowerShell 7.4+'s default
`$PSNativeCommandUseErrorActionPreference = $true`, a non-zero `typst` exit throws a
terminating error at line 76, so the line-78 `Fail "typst compile exited $LASTEXITCODE"` never
runs. The gate still exits non-zero (correct), but the tailored message is lost.
**Fix (CONVENTION):** Wrap the invocation in `try { & $typst … } catch { Fail "typst compile
failed: $_" }`, or set `$PSNativeCommandUseErrorActionPreference = $false` locally so the
explicit `$LASTEXITCODE` check governs.

---

_Reviewed: 2026-07-02_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
