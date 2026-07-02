---
phase: 10-manuscript-skeleton-and-related-work-positioning
plan: 05
subsystem: manuscript-build-gate
tags: [typst, regression-gate, ci, build-script, dod]
requires:
  - "manuscript/main.typ (compiles green)"
  - "manuscript/claims.typ (claim_table)"
  - "manuscript/sections/related_work.typ (matrix + four Tapqir axes)"
  - "manuscript/figures/specs.typ (>=7 SPEC-F entries, Supports tokens)"
  - "manuscript/lib/figstyle.typ + figures/f1_speedup.{yaml,typ}"
provides:
  - "manuscript/build.ps1 — Windows primary content-asserting regression gate"
  - "manuscript/build.sh — CI/fallback content-asserting regression gate (mirror)"
affects:
  - "all downstream manuscript phases (single scriptable regression gate for the skeleton)"
tech-stack:
  added: []
  patterns:
    - "exit-code-correct compile (no pipe) + post-compile content assertions"
    - "comment-aware grep (strip Typst // lines only; # is code in Typst)"
    - "falsifiable gate: injected bad value flips exit non-zero"
key-files:
  created:
    - .planning/phases/10-manuscript-skeleton-and-related-work-positioning/10-05-SUMMARY.md
  modified:
    - manuscript/build.ps1
    - manuscript/build.sh
decisions:
  - "Strip only Typst // comment lines (not #) — in Typst `#let`/`#table(` are CODE, so a `#` filter wrongly deleted the load-bearing claim_table/table( content"
  - "One falsifiability proof (inject 2cm into f1_speedup.typ) is sufficient per plan's OR clause; both gates flipped to exit 1 and reverted clean"
  - "Do NOT grep stderr for 'warning' — the non-fatal tiered-bib convergence warning is expected; success is decided by the compiler's exit code"
metrics:
  duration: ~20m
  completed: 2026-07-02
requirements: [D-12, D-14, D-16, D-17, SC1, SC2, SC3]
---

# Phase 10 Plan 05: Content-Asserting Dual Build Gate Summary

Upgraded both `manuscript/build.ps1` (Windows primary) and `manuscript/build.sh` (CI/fallback)
from compile-only gates to exit-code-correct compile + content-asserting regression gates that
verify the merged skeleton (claim table, related-work matrix + four Tapqir axes, ≥7 figure specs,
the D-07 Supports→claims cross-check, and the D-14/D-16 YAML-driven-F1 machinery), closing the
phase Definition of Done.

## What Was Built

Both gates now run in two steps:

1. **Compile (exit-code-correct, D-17/Pitfall 1):** `typst compile --root manuscript --font-path …`
   on a line with **no pipe**, so the compiler's own exit status is authoritative. The non-fatal
   `layout did not converge` warning from `lib/tiered-bib.typ` is expected and does **not** fail
   the gate — success is decided by exit code, never by grepping stderr for "warning".

2. **Content assertions (each exits non-zero on failure):**
   - **SC1** — `main.pdf` produced AND `claims.typ` contains `claim_table`.
   - **SC2** — `related_work.typ` contains all four Tapqir axes (`amortiz`; `sbc|calibrat`;
     `registration`; `spatial`) AND a `table(` matrix.
   - **SC3** — ≥7 `SPEC F[0-9]` entries in `figures/specs.typ` (found 7).
   - **D-07** — every `Supports: Cx` token in `specs.typ` (C1–C9) resolves to an `id: "Cx"` in
     `claims.typ`.
   - **D-14/D-16** — `lib/figstyle.typ` exists and defines `figstyle`; `figures/f1_speedup.yaml`
     exists; `f1_speedup.typ` loads it via `figstyle("f1_speedup")` and contains **no** hard-coded
     `cm` geometry and **no** inline `rgb("#…")` hex (proving F1 is fully YAML-driven).

Both gates are behaviorally mirrored (same assertion set, platform-idiomatic: bash uses
`set -euo pipefail` + `grep`; PowerShell uses `Get-NonComment`/`Select-String`/`Test-Path` +
`exit`).

## Verification

- `bash manuscript/build.sh` → exit 0, "OK: DoD met — compiles + claim table + matrix + >=7 specs (7) + YAML-driven F1", `manuscript/main.pdf` produced.
- `pwsh -File manuscript/build.ps1` → exit 0, same message and PDF (compiled here under the
  `$env:TYPST=typst-0.14` override hook; bash used `typst` 0.15.0 — both compile clean, confirming
  the D-17 override works).
- **Falsifiability:** appending `#let _injected_bad = 2cm` to `f1_speedup.typ` made **both** gates
  exit 1 (bash reported `GATE FAIL: D-14: f1_speedup.typ contains hard-coded cm geometry`); after
  revert the gate returned to exit 0 with a clean `git diff`.
- The integrated `main.pdf` renders the claim table, related-work matrix, figure specs, and the
  rendered F1 figure together (all four artifacts are compiled into the single PDF the gate asserts).

## Deviations from Plan

None — plan executed exactly as written. One implementation refinement worth noting (not a scope
deviation): the comment-strip guard filters only Typst `//` lines, **not** `#` lines. In Typst `#`
introduces code (`#let claim_table`, `#table(`), so a naive `#`-filter would delete the exact
load-bearing content the assertions check. YAML `#` comments are never grepped for content (the
`.yaml` is only `Test-Path`/`[ -f ]`-checked), so no YAML comment can self-satisfy an assertion.

## Constraint Compliance

- Only `manuscript/build.ps1` and `manuscript/build.sh` modified (`git status --short` clean of
  `src/`, `spike/`, `bayes.jl`, `colocalization.jl`, `STATE.md`, `ROADMAP.md`).
- BI template (`BayesInteractomics_Method_paper`) not written to.
- Compile line has no pipe in either gate; gate treats exit 0 as success and does not fail on the
  expected non-fatal convergence warning.

## Known Stubs

None — the gates assert against the fully-merged, completed Wave-2/3 manuscript content; no
placeholder or empty-value stubs were introduced.
