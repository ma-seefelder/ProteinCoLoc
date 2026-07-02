---
phase: 10-manuscript-skeleton-and-related-work-positioning
plan: 01
subsystem: manuscript
tags: [typst, bibliography, tiered-bib, figstyle, yaml, colours, libertinus, build-gate]

# Dependency graph
requires:
  - phase: 04-npe-training-advi-benchmark-ablation
    provides: honest Phase-4 speedup/RMSE numbers seeded into the claim/figure honesty caveats
provides:
  - "manuscript/ tree founded on the BayesInteractomics Typst template (lib/, styles/nature.csl, fonts/, colours.yaml, build.ps1) — a green, template-based skeleton"
  - "main.typ adapted to ProteinCoLoc: project.with() + two-tier tiered-bib(refs.bib) + modular section includes + <tab:claims> claim table + figure-specs include"
  - "lib/figstyle.typ per-figure YAML loader (D-14/D-15): figstyle(name) merges figures/<name>.yaml over documented defaults and resolves colour roles against colours.yaml"
  - "claims.typ claim-spine data + claim_table() render fn (single reserved-slot placeholder)"
  - "6 modular section stubs (introduction/related_work/methods/results/discussion/availability) + figures/specs.typ, all compile-safe"
  - "refs.bib seeded with load-bearing citations + folded self-cite + TODO Sci-Reports DOI"
  - "dual exit-code-correct compile gate: build.ps1 (primary) + build.sh (fallback)"
affects: [10-02 claim set, 10-03 related-work prose, 10-04 figure specs + worked F1, 10-05 content assertions, phase-16 assembly]

# Tech tracking
tech-stack:
  added: [Typst 0.15.0 manuscript toolchain, bundled Libertinus OTF fonts, Typst pkgs (abbr/zero/smartaref/wordometer/lilaq/alexandria via template), Nature CSL]
  patterns: [per-figure YAML styling convention (figstyle + colours.yaml roles), data-driven claim table (array-of-dicts + render fn), modular #include sections, two-tier bibliography, exit-code-correct dual build gate]

key-files:
  created:
    - manuscript/main.typ
    - manuscript/lib/figstyle.typ
    - manuscript/claims.typ
    - manuscript/refs.bib
    - manuscript/colours.yaml
    - manuscript/build.ps1
    - manuscript/build.sh
    - manuscript/figures/specs.typ
    - manuscript/sections/*.typ (6 stubs)
    - manuscript/lib/{template,helpers,curves,tiered-bib,abbreviations}.typ
    - manuscript/styles/nature.csl
    - manuscript/fonts/*.otf (14)
  modified:
    - .gitignore

key-decisions:
  - "figstyle colour helper is called as (fs.color)(role) — Typst forbids method-call syntax on dict-stored functions; documented in the module schema"
  - "Geometry values in per-figure YAML are plain numbers with documented units (cm/pt); the figure .typ applies the unit — keeps YAML round-trippable"
  - "Tiered references placed BI-style: H-tier list after Discussion, M-tier list after Methods (#set-tier(\"M\")); harmless with zero citations, later plans refine"
  - "manuscript/main.pdf gitignored as regenerable; per-figure PDFs (D-16 F1) stay trackable"

patterns-established:
  - "Per-figure YAML (D-14/D-15): every graphical parameter lives in figures/<name>.yaml; figstyle merges over defaults + resolves colour roles; nothing graphical hard-coded in .typ"
  - "Claim spine as diffable array-of-dicts rendered by claim_table() (single source of truth)"
  - "Compile gate is the falsifiable DoD: dangling @cite/@label -> exit 1; scripts never pipe the compile line"

requirements-completed: [D-01, D-02, D-03, D-04, D-05, D-11, D-12, D-13, D-15, D-17, SC1]

# Metrics
duration: 20min
completed: 2026-07-02
---

# Phase 10 Plan 01: Manuscript Skeleton (Template-Founded, YAML-Ready) Summary

**A green, BI-template-founded `manuscript/` tree compiling to `main.pdf` via both build.ps1 and build.sh, wired with a two-tier bibliography, modular section stubs, a data-driven claim table, and the new `figstyle` per-figure-YAML styling loader.**

## Performance

- **Duration:** ~20 min
- **Started:** 2026-07-02T11:15:00Z
- **Completed:** 2026-07-02T11:35:00Z
- **Tasks:** 2
- **Files modified:** 35 (34 created, 1 modified)

## Accomplishments
- Copied the BayesInteractomics Typst template foundation into `manuscript/` (D-13) — BI source repo verified untouched (git status clean at HEAD 3d06465).
- Re-authored `main.typ` for ProteinCoLoc using the template wiring: `#show: project.with(...)`, two-tier `tiered-bib(bib: "refs.bib", style: "styles/nature.csl", ...)`, ordered section `#include`s, the `<tab:claims>` claim table, and the figure-specs include; venue-neutral candidate-venues comment (D-02).
- Implemented `lib/figstyle.typ`, the per-figure YAML convention (D-14/D-15): `figstyle(name)` loads `/colours.yaml` + `/figures/<name>.yaml`, deep-merges over documented defaults (panel/axes/series/legend/text), and resolves colour roles against the palette with `#rrggbb` override. Full schema documented at the top of the file. Validated standalone (merge, defaults, role resolution, hex override, legend all assert-pass).
- Seeded `claims.typ` (spine data + `claim_table()`), `figures/specs.typ` stub, six citation-free section stubs, and `refs.bib` (load-bearing cites + folded self-cite + TODO Sci-Reports DOI).
- Proved the green baseline: `build.ps1` and `build.sh` both exit 0 and produce `manuscript/main.pdf`, under Typst 0.15.0 (installed) AND 0.14 (template pin, via the `$env:TYPST`/`TYPST` override). Gate is falsifiable: an injected dangling `@cite` drops the build to exit 1.

## Task Commits

Each task was committed atomically:

1. **Task 1: Copy BI template foundation + adapt main.typ + section stubs + refs.bib** - `e38c82a` (feat)
2. **Task 2: figstyle YAML loader + claims/specs stubs + build.sh + green-baseline proof** - `4a0c2f5` (feat)

## Files Created/Modified
- `manuscript/main.typ` - ProteinCoLoc-adapted entry: project.with + tiered-bib + section includes + claim table + figure specs.
- `manuscript/lib/figstyle.typ` - Per-figure YAML loader; `figstyle(name)` -> merged style dict + `(fs.color)(role_or_hex)` resolver; documented schema.
- `manuscript/lib/{template,helpers,curves,tiered-bib,abbreviations}.typ` - Copied verbatim from the BI template (page/text/heading/figure rules, helpers, curves, two-tier bib, abbreviations).
- `manuscript/claims.typ` - Claim-spine array-of-dicts + `claim_table()` render fn (single reserved-slot row; full set in Plan 10-02).
- `manuscript/refs.bib` - Seeded BibTeX store (ordabayev2022, costes2004, manders1993, talts2018, sainsburydale2024, radev2020, radev2023joss, self-cite, TODO Sci-Reports).
- `manuscript/sections/*.typ` - 6 modular, citation-free stubs with intent comments + TODO(phase N) markers.
- `manuscript/figures/specs.typ` - Compile-safe figure-spec stub (F1–F7 + worked F1 land in Plan 10-04).
- `manuscript/colours.yaml`, `manuscript/styles/nature.csl`, `manuscript/fonts/*.otf` - Copied palette / CSL / bundled Libertinus fonts.
- `manuscript/build.ps1` - Primary Windows gate: default installed `typst` (keep `$env:TYPST` override), `--root`, `--font-path`, exit-code-correct.
- `manuscript/build.sh` - Fallback/CI gate mirror (`set -euo pipefail`, no pipe on compile line).
- `.gitignore` - Ignore regenerable `manuscript/main.pdf`.

## Decisions Made
- **figstyle call syntax:** the colour helper is invoked as `(fs.color)(role)` — Typst forbids calling a dict-stored function with method syntax (`fs.color(x)` errors). Documented explicitly in the module schema; discovered during standalone validation.
- **Geometry as plain numbers:** per-figure YAML holds unit-free numbers (cm for width/height, pt for gutter, unitless scale); the figure `.typ` applies units. Keeps YAML language-neutral and round-trippable.
- **Reference placement:** followed the BI two-tier pattern (H-tier list after Discussion; M-tier after Methods with `#set-tier("M")`). With zero citations in the skeleton this is observationally inert; later plans that add `@cites` refine placement/tiers.
- **PDF artifact:** `manuscript/main.pdf` gitignored (regenerable); per-figure deliverable PDFs remain trackable.

## Deviations from Plan

None - plan executed exactly as written. No Rule 1–4 deviations were required (no pre-existing bugs, missing critical functionality, blocking issues, or architectural changes). Template packages auto-downloaded from the Typst registry and compiled cleanly under both 0.15.0 and 0.14, so no minimal-adjustment TYPST-COMPAT fix was needed.

## Issues Encountered
- **figstyle name-shadowing (self-authored, fixed pre-commit):** the colour closure was initially named `color`, shadowing Typst's built-in `color` type used in its own `type(role) == color` guard. Renamed the closure to `_color` and exposed it as the `color` dict key. Caught before the Task 2 commit via the standalone figstyle self-test.
- **figstyle call syntax:** standalone validation surfaced that `fs.color(...)` errors ("cannot directly call dictionary keys as functions"); corrected the documented usage to `(fs.color)(...)`. No code change needed beyond docs — the helper itself was already correct.

## User Setup Required
None - no external service configuration required. Typst 0.15.0 is preinstalled; template packages auto-download from the Typst registry.

## Next Phase Readiness
- Green, wired baseline is ready: Plan 10-02 seeds the full claim set into `claims.typ`; Plan 10-03 authors the related-work matrix + Tapqir four-axis prose in `sections/related_work.typ` (and must verify the three ⚠ Tapqir cells against Ordabayev et al. 2022 before locking); Plan 10-04 writes figure specs + builds the worked F1 figure via the `figstyle` pattern; Plan 10-05 adds content assertions to the build scripts.
- `main.typ` is fully wired so later plans edit only their own files (section stubs, claims.typ, figures/*), avoiding merge collisions.
- No blockers. TYPST-COMPAT constraint satisfied without degradation (no figure disabled, no from-scratch fallback, no YAML-loader bypass).

---
*Phase: 10-manuscript-skeleton-and-related-work-positioning*
*Completed: 2026-07-02*
