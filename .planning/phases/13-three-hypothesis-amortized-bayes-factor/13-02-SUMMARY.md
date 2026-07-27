---
phase: 13-three-hypothesis-amortized-bayes-factor
plan: 02
subsystem: pre-registration governance (planning artifacts only — no executable byte)
tags: [amendment, pre-registration, honesty, roadmap, d-12, d-15, d-02, docs]
requires:
  - ".planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md (READ-ONLY precedent: header shape, §0 provenance disclosure, legitimacy test, original-not-erased clause, numbered-defect structure)"
  - "docs/amortized.md §Named limits (READ-ONLY: limit 2 prior atoms / randomized ranks, limit 3 KDE invalid past |logBF| ~ 6.9 — quoted verbatim)"
  - ".planning/spikes/014-bf-sim-validation/README.md (READ-ONLY: 34.75% = 139/400 non-finite, 32.5% = 130/400 saturated, dropped median |Delta rho| 0.736 vs kept 0.283, AUC 0.9939)"
  - "spike/p13/consts.jl from plan 13-01 (READ-ONLY: every gate number is cited by constant name, never re-stated as a fresh literal)"
  - "13-RESEARCH.md sections A3, F1, H1-H3, J1, J2, J4-J6 and 13-PATTERNS.md lines 910-966 + Pattern S10"
provides:
  - ".planning/phases/13-three-hypothesis-amortized-bayes-factor/13-SC2-AMENDMENT.md — the frozen SC1/SC2/SC3 amendment (486 lines, sections 0-7 plus a provenance appendix)"
  - "the greppable declaration `No pass/fail threshold is defined for any real-image quantity` that plans 13-13/13-15/13-16 and the 13-14 report bind to"
  - "the three-arm SC3 standing table (simulator ground truth GATES; alpha series and test/test_images/ do not)"
  - "a ROADMAP Phase-13 block that points at the amendment from all three criteria and names the D-02 Phase-11 dependency"
affects:
  - "plan 13-14's report: it must print the original SC1/SC2/SC3 beside the amended criteria (section 6 makes that an obligation, not a courtesy)"
  - "plans 13-13, 13-15, 13-16: section 5(g) forbids any real-image pass/fail bar; section 5(d) requires every real-image number beside its OOD verdict and its lambda"
  - "plan 13-12: the gated/reported split of section 3 is the criterion it reports against"
  - "/gsd:verify-work: the Phase-13 roadmap block no longer scores the phase against superseded SC2 text"
tech-stack:
  added: []
  patterns:
    - "frozen-specification amendment document (07-GATE-AMENDMENT.md:1-34)"
    - "per-defect outcome-independence argument + the blind-reviewer question (07-GATE-AMENDMENT.md:38-68)"
    - "original-is-not-erased clause (07-GATE-AMENDMENT.md:28-34)"
    - "reported-not-gated real-data arm, five obligations (13-PATTERNS Pattern S10)"
    - "every real-image number ships with its OOD verdict and its lambda (13-PATTERNS Pattern S11)"
    - "gate numbers cited by constant name, never re-typed as literals (spike/p13/consts.jl as the single binding source)"
key-files:
  created:
    - .planning/phases/13-three-hypothesis-amortized-bayes-factor/13-SC2-AMENDMENT.md
  modified:
    - .planning/ROADMAP.md
decisions:
  - "The amendment INVERTS Phase 7's disclosure and says so in section 0: Phase 7 amended a gate it had already missed; Phase 13 declines to adopt a reference a prior phase already retired and published as named limit 3. The document explicitly instructs the reader not to blur the two positions."
  - "Both defects are argued so they are decidable with every artifacts directory deleted: Defect A is arithmetic on L = 999 in src/bayes.jl; Defect B is which arrays enter which function plus the algebra of two training objectives (shuffled-theta NRE identity vs plain BCE-on-class-labels, so the D-09 heads need the -log(q_A/q_B) correction the shipped net does not)."
  - "The section-0.1 sub-disclosure records that the SC3 scope was itself amended (commit 5b4da6d) after an audit correction, rather than quietly adopting the new scope."
  - "Section 5 makes the real-image honesty consequence part of the FROZEN document instead of leaving it to the report, so a qualitative reading cannot later be quoted as a validation of the exclusion hypothesis."
  - "Section 5(g) names the two documented exemptions (P13_REAL_ANCHOR_MBAR/TOL and the P13_ALPHA_* invariant bounds) as ingestion/transform regressions, so 'no real-image threshold' is not contradicted by the loader's own self-check."
  - "Section 7 scopes the single iteration allowance to the GATE and forbids spending it on a real-image observation — a real-image quantity has no bar to miss."
  - "An appendix labels the provenance of every number in the file (shipped net / frozen simulator calibration / prior spike / pre-registered choice), so no reader can mistake a quoted number for a Phase-13 result."
  - "ROADMAP's SC2 keeps its literal `reproduces compute_BayesFactor()` wording; the amendment pointer redirects the reader rather than rewriting the criterion."
metrics:
  duration: ~30 min
  completed: 2026-07-27
  tasks: 2
  files: 2
  commits: 2
  tests: n/a (documentation plan — verified by the plan's own julia assertion harness, both PASS)
---

# Phase 13 Plan 02: SC2/SC3 Pre-Registration Amendment Summary

The Phase-13 success criteria are now amended on the record **before any Phase-13 number exists**,
at commit `df16e76941cefd28147efb4d4c91278597836cbd`, with two outcome-independent defects, the
replacement criterion stated by constant name, and SC3 split into three arms of which only the
simulator ground truth gates.

## Provenance — the ordering claim, checkable

| Field | Value |
|---|---|
| Commit sha recorded in the amendment's section 0 | **`df16e76941cefd28147efb4d4c91278597836cbd`** |
| Date written | 2026-07-27 |
| Phase-13 results in existence at that sha | **none** |

**No Phase-13 result existed at that sha.** Verified directly, not asserted: `ls spike/p13/` returned
exactly `consts.jl` — no trained net, no `*report*.jld2`, no evaluation artifact of any kind. The only
other Phase-13 executable in the tree is `spike/test/test_p13_consts.jl` (plan 13-01's literal-assertion
test). `P13_TAU` is still absent from the Tier-1 block, which is a second independent check that this
document precedes the phase's measurements.

## What was built

### Task 1 — `13-SC2-AMENDMENT.md` (486 lines, commit `3bd88e4`)

Sections 0 through 7 plus a provenance appendix, mirroring `07-GATE-AMENDMENT.md`'s shape:

- **§0 Provenance disclosure** — states the inversion of Phase 7's posture, records the date and the
  full 40-character sha, quotes the legitimacy test verbatim, and answers it with "stronger than yes:
  the justification is not a Phase-13 result at all". **§0.1** additionally discloses that the §5 scope
  was itself amended (commit `5b4da6d`) after an audit correction — the first reading of D-15 counted
  the corpus anchors correctly but did not read their `split` and `sha256` columns.
- **§1 Defect A** — the KDE reference's validity ceiling: L = 999 draws ⇒ smallest resolvable tail
  probability 1/999 ⇒ `|log BF| > log(999) ≈ 6.9` is kernel-tail extrapolation, not measurement. Cites
  `references/bayes-factor-gate.md` §008 and quotes spike 014's attrition table (34.75% = 139/400
  non-finite, 32.5% = 130/400 saturated, dropped median `|Δρ|` 0.736 vs kept 0.283, NRE attrition 0),
  notes attrition is 100% baseline-side **and not random** (it fails preferentially on high-signal
  pairs, so any agreement statistic runs on a survivorship-biased 261/400). Quotes named limit 3
  exactly. Closes on the point that the invalid band **is** the confident-call tail.
- **§2 Defect B** — architectural, needs no results: D-02 puts the net on Phase 11's registration-aware
  frozen `zt` plus a λ conditioning input, an input surface shared with neither
  `compute_BayesFactor()`'s Turing/ADVI path nor the shipped binary NRE. Records the subtler
  consequence — the shipped logit is a shuffled-θ likelihood-to-evidence ratio while the D-09 heads are
  plain BCE classifiers carrying `log(q_A/q_B)`, so the two are **not numerically interchangeable even
  at identical inputs**, which is exactly why the continuity check is reported and not gated. Explicit
  non-change clause for `src/amortized/bf.jl` (frozen by D-01; the `−0.0102` term is below the
  reporting floor).
- **§3 Replacement criterion** — gated: per-class one-vs-random AUC with `P13_AUC_FLOOR_COLOC = 0.90`
  and `P13_AUC_FLOOR_EXCLUSION = 0.90` on `P13_GATE_M = 4000` pairs with `P13_MIN_EVAL_PER_HEAD = 1000`;
  per-head `ECE ≤ P13_ECE_GREEN = 0.05` with MCE reported-beside-its-empty-bin-count and the
  `P13_VACUOUS_AUC_FLOOR = 0.60` guard. Reported-not-gated: the 3×3 confusion matrix
  (`:argmax_descriptive`, decisions are Phase 14), binary-NRE continuity, the λ response, the α series,
  and every real-image quantity. Includes the D-05 label rule with the `ρ = 0.8` under `ρ = 0.9`
  counter-example.
- **§4 SC1's surviving content** — `RatioEstimator` ruled out (fixed 1-unit head, hard-coded
  logit-BCE, a passed loss silently ignored in v0.2.1), "simplex" ruled out (two free numbers plus a
  structural zero). Bindable criterion: one network, one forward pass, `log BF(coloc : random)` and
  `log BF(exclusion : random)` with the random entry identically 0.
- **§5 SC3 scope** — the three-arm table plus (a) why the α series is load-bearing, (b) why the corpus
  anchor is not used and that the seal-break escape hatch is **withdrawn**, (c) the honesty consequence,
  (d) the two strengthening findings, (e) the λ rule, (f) the ladder scope limit, (g) the declaration.
- **§6 The original is not erased** — ROADMAP SC1/SC2/SC3 stay byte-present and citable; any report
  citing an amended criterion must cite the original alongside it, and plan 13-14 is required to print
  both.
- **§7 Iteration allowance** — `P13_ITERATION_ALLOWANCE = 1`, the single pre-declared trigger
  (exclusion AUC below floor *specifically at high `|ρ|`* ⇒ switch to `:importance_weighted`, D-07-ii),
  no second amendment authorised, and the allowance explicitly cannot be spent on a real-image
  observation.
- **Appendix** — a provenance table labelling every quoted number as shipped-net / frozen-simulator /
  prior-spike / pre-registered-choice. **No measured Phase-13 number appears anywhere in the file.**

### Task 2 — ROADMAP Phase-13 annotation (commit `6dbdbc7`, +7 −6)

Three edits, nothing else, all inside the `### Phase 13:` block:

1. `**Depends on**: Phase 7` → `**Depends on**: Phases 7, 11 (D-02: … execution blocks on Phase 11's
   research net existing; the original entry named Phase 7 only)`.
2. An `**AMENDED**` pointer line under `**Success Criteria**`, plus a per-criterion marker
   (`(see amendment section 4)`, `(see amendment sections 1-3)`, `(see amendment section 5 — three
   arms; …)`). The three criteria are otherwise byte-unchanged, including SC2's literal
   `reproduces \`compute_BayesFactor()\`` citation.
3. Two description-only checklist corrections: the `13-01-PLAN.md` line now names **D-15** (matching
   that plan's `requirements:` frontmatter), and the `13-04-PLAN.md` line's falsified
   **strict positivity** invariant is replaced in place with `no NEW zero (count(iszero) preserved)`
   plus `substrate-agnostic (one code path for simulated and real)`.

`**Requirements**: TBD` left as-is (no REQ-IDs minted). No other phase block touched.

## Section 5 standing table (reproduced so the split is visible without opening the amendment)

| Arm | Substrate | Standing |
|---|---|---|
| **Simulator ground truth** | `ρ_true < −τ`, exact labels, well-powered, `P13_GATE_M = 4000` | **GATES** — the D-12/D-13 floors of §3 |
| **D-16 α-graded series** | mask-based disjoint reassignment on simulated **and** real pairs | Supporting evidence — **not a gate** (`P13_ALPHA_GATED = false`) |
| **`test/test_images/` qualitative check** | six committed 1028 × 1376 microscopy TIFFs | Required deliverable — **explicitly not a gate** |

## Verification — observed results, not paraphrase

Every command below was run and its real output observed.

| # | Check | Observed |
|---|---|---|
| T1 | `julia -e '… @assert occursin(…) … println("AMENDMENT_OK")'` (7 assertions) | printed **`AMENDMENT_OK`** |
| T1 | 40-char hex sha present | **`df16e76941cefd28147efb4d4c91278597836cbd`** |
| T1 | `^## ` heading count ≥ 8 | **9** (§0-§7 + Appendix) |
| T1 | `GATES\|not a gate` line count ≥ 3 | **3** |
| T1 | `no colocalization ground-truth label\|cannot show the verdicts are correct` ≥ 1 | **2** |
| T1 | `Phase 16` ≥ 2 | **3** |
| T1 | `biological` ≥ 1 / `433.69\|179.14` ≥ 1 | **1** / **3** |
| T1 | `No pass/fail threshold is defined for any real-image quantity` ≥ 1 | **1** |
| T1 | `open_sealed_holdout` ≥ 1, every occurrence in prose stating it is NOT called | **1**, in §5(b) prose; no code block |
| T1 | `named limit` (key-link to `docs/amortized.md`) | **6** |
| T1 | file length vs `min_lines: 140` | **486 lines** |
| T2 | `julia -e '… println("ROADMAP_OK")'` (11 assertions incl. the 16-line count) | printed **`ROADMAP_OK`** |
| T2 | `13-15-PLAN.md\|13-16-PLAN.md` whole-file count ≥ 2 | **2** |
| T2 | `^\*\*Plans\*\*: 16 plans` in block | **1** |
| T2 | `^- \[ \] 13-` in block | **16** |
| T2 | `real-image ingestion` / `real-image run` on the `**Plans**:` line | **1** / **1** |
| T2 | 13-15 blocked / 13-16 blocked | **0** / **1** (correct) |
| T2 | `**Plans**:` line byte-unchanged (`git diff -U0 \| grep '^-\*\*Plans\*\*'`) | **0 matches** — not rewritten |
| T2 | removed checklist lines (`git diff -U0 \| grep '^-- \[ \] 13-'`) | **2**, and they are the 13-01 and 13-04 lines |
| T2 | `Phases 7, 11` / `^\*\*AMENDED\*\*` / `see amendment` | **1** / **1** / **3** |
| T2 | original SC2 wording retained in block | **1** |
| T2 | `strict positivity` / `count(iszero) preserved` in block | **0** / **1** |
| T2 | 13-01 line names D-15 | **1** |
| T2 | `run /gsd:plan-phase 13 to break down` (stale placeholder) | **0** |
| T2 | `git diff -U0 -- .planning/ROADMAP.md` hunk bodies contain another `### Phase` heading | **no** — all four hunks are inside the Phase-13 block |
| Both | `git diff --quiet HEAD -- src/ spike/ test/ corpus/` | **exit 0** after each task |
| Plan v8 | no `spike/p13/*report*.jld2` at commit time | `ls spike/p13/` → **`consts.jl`** only |
| Both | post-commit deletion check (`git diff --diff-filter=D HEAD~1 HEAD`) | **empty** on both commits |

`spike/Project.toml` and `spike/Manifest.toml` are byte-unchanged; **no package was installed**. This
plan touched no `.jl` file and no `src/`, `test/` or `corpus/` byte.

## Deviations from Plan

None — the plan executed exactly as written. Both tasks' `<action>` blocks were followed literally,
including the explicit prohibition on regenerating the ROADMAP plan inventory.

### Observation recorded, deliberately NOT fixed

The `13-04-PLAN.md` ROADMAP line still opens `— α-graded disjoint-reassignment transform + its **four**
invariants (…)` while the list it introduces now carries **five** entries, because edit 3 mandated both
replacing `strict positivity` **and** appending `substrate-agnostic (one code path for simulated and
real)`. The plan's edit 3 says *"Change nothing else on the line"* twice, and none of its acceptance
criteria touch the word "four", so the count word was left alone rather than silently widening the
edit. Flagged here so plan 13-14's end-of-phase tick-off pass can correct the numeral in the same
sweep if wanted; it is a cosmetic numeral in a checklist description and gates nothing.

## Assumption Drift (advisory)

**Planned:** the plan's `<interfaces>` block asserts `**Depends on**: Phase 7` is *the* Phase-13
dependency line to edit.
**Actual:** that exact string occurs **twice** in `.planning/ROADMAP.md` (line 234 in the Phase-11
block, line 267 in Phase 13), so a bare targeted edit on it would have been ambiguous and could have
landed in another executor's concurrently-edited block.
**Why it matters / how handled:** the edit was anchored on the Phase-13 `**Goal**` line's tail
(`replacing the fragile KDE+quadgk Bayes factor`) together with the criteria block, making the match
unique and provably inside Phase 13. The resulting diff confirms it: all four hunks sit at lines
267-278. No behavioural change to what the plan asked for — only to how it was located.

## Known Stubs

None. This plan produced no code and no placeholder values; both artifacts are complete documents.

## Threat Flags

None. This plan created no network endpoint, no auth path, no file-access pattern and no schema
change. It reduced one threat surface: §5(b) of the amendment records in prose that
`open_sealed_holdout` is not called anywhere in Phase 13 (T-13-05), and §5(g)'s verbatim declaration
plus §7's scope clause close T-13-54.

## Self-Check: PASSED

- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-SC2-AMENDMENT.md` — FOUND (486 lines)
- `.planning/ROADMAP.md` — FOUND, modified, Phase-13 block only
- commit `3bd88e4` (`docs(13-02): freeze the SC2/SC3 amendment before any Phase-13 result`) — FOUND
- commit `6dbdbc7` (`docs(13-02): annotate the Phase-13 roadmap block …`) — FOUND
- recorded provenance sha `df16e76941cefd28147efb4d4c91278597836cbd` — FOUND, is the real parent of
  commit `3bd88e4`
