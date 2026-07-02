---
phase: 10-manuscript-skeleton-and-related-work-positioning
fixed_at: 2026-07-02T10:42:36Z
review_path: .planning/phases/10-manuscript-skeleton-and-related-work-positioning/10-REVIEW.md
iteration: 1
findings_in_scope: 5
fixed: 5
skipped: 0
status: all_fixed
---

# Phase 10: Code Review Fix Report

**Fixed at:** 2026-07-02T10:42:36Z
**Source review:** .planning/phases/10-manuscript-skeleton-and-related-work-positioning/10-REVIEW.md
**Iteration:** 1

**Summary:**
- Findings in scope: 5 (0 critical, 5 warning — `fix_scope = critical_warning`)
- Fixed: 5
- Skipped: 0

All work was done inside an isolated git worktree on `main` and each fix was
re-validated with `bash manuscript/build.sh` (compile + DoD content gate) before
commit. The four IN-* Info findings were out of scope and left untouched.

## Fixed Issues

### WR-01: `main.typ` bibliography comment is factually wrong

**Files modified:** `manuscript/main.typ`
**Commit:** c9a7c26
**Applied fix:** Rewrote the load-bearing comment above the `tiered-bib` call. It no
longer claims "no section stub cites any key yet"; it now states that Related Work
already cites 4 keys at the H tier (so the main "References" list renders non-empty),
that the Methods list stays empty until Phase 5/7 wire M-tier `@cites`, and that every
cited key must exist in `refs.bib`. Comment-only change; build stays green.

### WR-02: `figstyle` schema/consumer contract gap — `axes.yscale` has no default

**Files modified:** `manuscript/lib/figstyle.typ`
**Commit:** 93072d9
**Applied fix:** Added `yscale: "linear"` and `xscale: "linear"` to `_DEFAULTS.axes`
(the consumer `f1_speedup.typ:52` reads `fs.axes.yscale`, previously only satisfied
because the YAML happened to set it) and documented both keys in the schema comment.
Also documented `key:` as a REQUIRED per-series field (the consumer indexes series via
`_ser.insert(s.key, s)`), closing the "any key you omit falls back to a default"
contract gap for figures that clone the F1 pattern.

### WR-03: Worked F1 exemplar renders only panel (a); spec/caption promised two panels

**Files modified:** `manuscript/figures/specs.typ`
**Commit:** d5ab187
**Applied fix:** Chose the "reconcile to a single-panel figure" option (the lower-risk
of the two the reviewer offered, and the one that keeps the honesty numbers intact).
The SPEC F1 box now describes a single panel (a) and states explicitly that the
honesty-critical 90% interval-width comparison (NPE 0.884 vs ADVI 0.605, the
`interval_width` column of `speedup.csv`) is reported in the caveat/caption rather
than as a separate rendered panel — so the spec and the one-diagram implementation now
agree exactly. The rendered caption already carried these numbers and promised no "(b)"
panel, so it was left unchanged. `SPEC F1` / `Supports: C1, C2` tokens preserved (gate
assertions still pass).

### WR-04: Build-gate "no hard-coded geometry" check is trivially bypassable

**Files modified:** `manuscript/figures/f1_speedup.typ`, `manuscript/figures/f1_speedup.yaml`,
`manuscript/lib/figstyle.typ`, `manuscript/build.sh`, `manuscript/build.ps1`
**Commit:** d82e105
**Applied fix:** Strengthened the D-14 "nothing graphical hard-coded" convention that
the user explicitly values, in two coordinated moves:

1. **Moved every remaining hard-coded length in `f1_speedup.typ` into YAML.** The legend
   swatch size/radius/baseline, swatch–label gap, legend column/row gutters, the
   legend↔panel gap, and the standalone-preview page margin are now read from the YAML
   (`fs.legend.*`, `fs.panel.preview_margin`) with the unit applied only via a sanctioned
   idiom (`* 1em`, `* 1pt`, or `_cm(...)` → `* 10mm`). New defaults for all of these were
   added to `figstyle`'s `_DEFAULTS` and documented in its schema comment (so the loader's
   fallback contract covers every key F1 consumes), and surfaced in `f1_speedup.yaml`.
2. **Tightened both build gates.** The old check greped only the literal unit `cm`. It now
   strips the sanctioned conversion idioms (`* 10mm` / `* 1pt` / `* 1em`) and flags any
   REMAINING numeric length literal in `cm|mm|in|pt|em` — so a future `width: 5mm` or
   `width: 200pt` (exactly the bypasses the reviewer named) trips the gate. A positive
   assertion was also added: the panel width/height must be YAML-driven via
   `_cm(fs.panel.*)`. Falsifiability was verified (probe `width: 5mm`, `x: 200pt`, `y: 2in`
   all trip; the sanctioned idioms do not), and both `build.sh` and `build.ps1` run green.

### WR-05: F1 colours NPE/ADVI series with AP-MS-semantic palette roles

**Files modified:** `manuscript/colours.yaml`, `manuscript/figures/f1_speedup.yaml`
**Commit:** 4b52342
**Applied fix:** Took the reviewer's recommended option — added a colocalization-native
`inference_method` role group (`npe_forward`, `npe_full`, `advi`) to `colours.yaml` and
pointed `f1_speedup.yaml` at those roles instead of the AP-MS `model_series.meta_learner`
/ `meta_learner_unc` roles (and the ADVI point at `inference_method.advi` instead of
`reference.axis_chrome`). The hex values were preserved, so the figure renders identically;
only the semantic role names change, so an editor recolouring an NPE curve now edits a role
whose name means "NPE". On the D-13 "copied read-only" tension: this was done by EXTENDING
the palette with a new semantic role — which `colours.yaml`'s own header explicitly sanctions
("extend this file with a new semantic role instead") — not by re-meaning any existing role.
The dashed 1× reference line keeps `reference.chance_line`, which is semantically correct for
a baseline/chance line.

## Skipped Issues

None — all in-scope findings were fixed.

Out of scope (Info tier, `critical_warning` scope): IN-01 (F1 spec box consumes a figure
number), IN-02 (`build.sh set -e` suppresses SC3/D-07 diagnostics), IN-03 (`refs.bib`
placeholder `TODO:`/`\url{}`), IN-04 (`build.ps1` non-zero-exit skips custom `Fail`). These
were not addressed.

---

_Fixed: 2026-07-02T10:42:36Z_
_Fixer: Claude (gsd-code-fixer)_
_Iteration: 1_
