---
phase: quick-260803-jm5
plan: 01
subsystem: planning-record
tags: [roadmap, phase-12, closure, documentation-only, negative-result]
requires: []
provides:
  - ".planning/ROADMAP.md — Phase-12 closure record matching 12-VERIFICATION.md"
affects: []
tech-stack:
  added: []
  patterns: []
key-files:
  created: []
  modified:
    - ".planning/ROADMAP.md"
decisions:
  - "Cite the artifact's literal `spat07_scope = spat08_scope = :deferred_to_v2_1` rather than the plan's paraphrased `:deferred` — source-of-truth rule."
  - "Cite `12-17-PLAN.md:118` (the `p12_full_arms()` line that actually reads `P12_CHOSEN_PRIOR`) rather than the paraphrased `:137`."
  - "Cite `p12_consts.jl:798` for the `P12_CHOSEN_PRIOR` negative assertion specifically; `:799` is the `P12_N_LOW` assertion."
metrics:
  duration: ~25 min
  completed: 2026-08-03
status: checkpoint-pending
---

# Quick Task 260803-jm5: Amend ROADMAP.md — Phase 12 to CLOSED NEGATIVE — Summary

`.planning/ROADMAP.md` now records Phase 12 as executed-and-closed-negative in agreement with
`12-VERIFICATION.md`, replacing an open `[ ]` bullet and a `**Plans**: TBD` placeholder that read
"unplanned, pending" where the truth is "executed, closed, answered NO."

## What changed (one file, four places)

1. **Summary bullet (Phases list)** — flipped `[ ]` → `[x]`, existing text and the
   `(descope-to-v2.1 candidate)` parenthetical preserved verbatim, then a **CLOSED NEGATIVE
   2026-08-03** annotation appended in Phase 11's em-dash-then-bold shape carrying all four required
   claims: goal NOT achieved; `NONE-BEATS-ABLATION` in three independent mini-spike runs, so the
   shipped deliverable IS the ablation wearing the `SpatialColocResult` type and no trained spatial
   arm exists at any scale; SPAT-06 measured **0.9707** against **[0.87, 0.93]**; SPAT-07 and SPAT-08
   DEFERRED TO v2.1 by recorded user rulings dated 2026-07-31. Ends with a `12-VERIFICATION.md`
   pointer.
2. **Three success criteria** — annotated inline with bold verdicts: (1) MECHANISM BUILT; NO ARM WAS
   EVER SELECTED, (2) DELIVERED with the §SPAT-05 named caveat (non-spatial D-13 ablation, spike
   scale 10,000 pairs), (3) NOT MET, AND UNREACHABLE AS WRITTEN with both independent reasons. An
   `OUTCOME 2026-08-03` header sits above them, mirroring Phase 11's. **The `**AMENDED**` paragraph
   pointing at `12-SC3-AMENDMENT.md` is untouched** — it still precedes the criteria, byte-unchanged.
3. **Plan list** — `**Plans**: TBD` / `- [ ] TBD (run /gsd:plan-phase 12 …)` replaced with all 20
   real plans: 16 `[x]` (12-01 … 12-16), 4 `[S]` FORECLOSED and never run (12-17 … 12-20), each
   description derived from that plan file's own `<objective>` block, each foreclosure given its
   specific reason.
4. **Closure block** — mirrors Phase 11's, cites both source documents with the 8/9 must-have score
   broken out, gives the legend, and carries the two things the record says must survive.

## Verification actually run (not inferred)

| Gate | Command | Observed |
|---|---|---|
| Task 1 | the plan's five `grep` assertions | printed `TASK1_OK` |
| Task 2 | the plan's eleven assertions incl. both isolation checks | printed `TASK2_OK` |
| §10.5 sentence | Python exact-substring compare, source (unwrapped) vs ROADMAP | `in source: True`, `in roadmap: True`, `occurrences: 1` |
| Isolation | `git diff --name-only -- spike src test .planning/phases .planning/STATE.md` | **zero lines** |
| Scope | `git diff --name-only -- .planning/ROADMAP.md` | exactly 1 |
| Deletions | `git diff --diff-filter=D --name-only HEAD~1 HEAD` | none |

Section-level counts confirmed: zero `TBD` in the Phase-12 section, exactly 20 plan bullets
(16 `[x]` + 4 `[S]`), exactly one `**Closure**` block.

## Commits

| Hash | Message |
|---|---|
| `ae24464` | `docs(quick-260803-jm5): Phase 12 summary bullet + SC verdicts -- CLOSED NEGATIVE` |
| `d860461` | `docs(quick-260803-jm5): Phase 12 real 20-plan list + Closure block` |

Base commit for the run was `35848a8`. Both commits touch `.planning/ROADMAP.md` and nothing else.

## Source-vs-plan discrepancies found (source wins, per the plan's absolute rule)

Three. None changes any verdict; all three are the plan's convenience paraphrase being slightly off
the artifact, and the roadmap was written from the artifact in each case.

1. **`spat07_scope` value.** The plan (and `12-VERIFICATION.md`'s frontmatter prose) writes
   `spat07_scope = :deferred`. The **artifact and the code write `:deferred_to_v2_1`** —
   `spike/validation/run_p12_coverage.jl:698` (`"spat07_scope" => :deferred_to_v2_1`) and
   `12-16-SUMMARY.md:245` (`spat07_scope = spat08_scope = :deferred_to_v2_1`). ROADMAP.md carries
   the literal artifact value.
2. **12-17's `P12_CHOSEN_PRIOR` line.** The plan and STATE.md cite `12-17:137`. Line 137 is the pool
   resolution step; the constant is actually read at **`12-17-PLAN.md:118`**
   (`p12_full_arms()` = `(P12_CHOSEN_PRIOR, :none)`), with further uses at `:161` and `:172`.
   ROADMAP.md cites `:118`.
3. **The negative-assertion line.** The plan cites `spike/validation/p12_consts.jl:798-799` for
   `P12_CHOSEN_PRIOR`. Line **798** is `@assert !isdefined(@__MODULE__, :P12_CHOSEN_PRIOR)`;
   **799** is the `P12_N_LOW` assertion. ROADMAP.md cites `:798`.

Everything else re-read at write time and confirmed exact: `0.9707`, `[0.87, 0.93]`, the
`[0.9675, 0.9738]` interval and its ≈0.037 clearance, W1 = `0.00675` at M = 20,000 over 11
arm×rung configurations against `P12_SIM02_W1_TOL_PERREGION = 0.10`, `eps_ratio = 0.96307`,
`ridge_residual_shrinkage = 0.96870`, `+0.19195` nats/region against `P12_STAGE2_LOGSCORE_MIN = 0.02`,
the 12/12 golden match with `maximum(abs, Δ) = 0.0`, `car_sigma`@`p12_lattice.jl:125` /
`gp_sigma`@`:153`, `build_p12_summary_net`@`p12_architecture.jl:253`, `p12_coloc_map`@
`p12_coverage.jl:990`, `12-STAGE1-VERDICT.md:3` = `VERDICT: PROCEED` at `ee78cfe`,
`12-19-PLAN.md:202`, `run_p12_minispike.jl:249` = `for a in (:car, :gp)`, the 0.0492/0.0509 nominal
distances, and `p12_train_full_report.jld2` absent on disk (`ls` → No such file or directory).

## Deviations from Plan

**None** — plan executed as written. No deviation rule fired. No auto-fix was required.

## Assumption Drift (advisory)

None. Every figure written was read from the source at write time; where the plan's paraphrase and
the source disagreed, the source was used and the divergence is recorded above rather than absorbed.

## Constraint compliance

- Only `.planning/ROADMAP.md` modified. `spike/`, `src/`, `test/`, `.planning/phases/**` and
  `.planning/STATE.md` all show zero diff lines.
- **`12-11-SUMMARY.md` left BYTE-UNCHANGED** with its stale `status: blocked` frontmatter; the
  staleness is documented in ROADMAP.md prose, including an explicit "do NOT fix the plan file"
  instruction and the 5-vs-4 SDK count explanation.
- No threshold examined or moved; no compute run; nothing re-trained or re-seeded.
- The task's deliberate override of quick-mode's "do not touch ROADMAP.md" rule was honored, as the
  plan directs.
- STATE.md and the PLAN/SUMMARY docs are **not** committed here — left for the orchestrator.

## Task 3 — PENDING HUMAN VERIFICATION (BLOCKING)

Task 3 is `checkpoint:human-verify` and this run is non-interactive, so it was **not** executed.
Tasks 1 and 2 are complete and committed. The user must confirm, by eye rather than by grep, because
what is at stake is tone and the rulings transcribed are theirs:

1. **The one question that matters: does the Phase-12 section read as a NO?** Compare against
   `12-VERIFICATION.md`'s headline — *"NO — the phase did not achieve its stated goal."* If the
   roadmap reads softer than the verification, name the sentences and they get sharpened.
2. **The §10.5 sentence** — *"We cannot exclude that longer training would have brought the spatial
   arms into the coverage band."* — appears verbatim in the Closure block and is labelled a **NAMED
   LIMIT of v2.0, not a caveat**. Confirm the wording and the label.
3. **The no-retune ruling** names all three refused moves (the `run_p12_minispike.jl:249`
   admissibility loop, the `[0.87, 0.93]` band re-derivation, the early-stopping criterion) and
   states that **the loop stays wrong and documented as wrong**. Confirm nothing was softened.
4. **The 12-11 entry** states its `status: blocked` frontmatter is stale, is deliberately left
   byte-unchanged, and is why the SDK counts 5 incomplete plans rather than 4.
5. **Nothing over-claimed in the other direction** — no added blame, no editorializing about
   execution quality. The `OUTCOME` header carries the verification's own *"not read here as an
   execution failure"* framing.

**Resume signal:** "approved", or name the sentences that are too soft / too strong / wrong.

## Self-Check: PASSED

- `FOUND: .planning/ROADMAP.md`
- `FOUND: ae24464`
- `FOUND: d860461`
