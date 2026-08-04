# Phase 16: External Validation and Manuscript Assembly - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in `16-CONTEXT.md` — this log preserves how they were reached.

**Date:** 2026-08-04
**Phase:** 16-external-validation-and-manuscript-assembly
**Mode:** discuss `--analyze --all` (all gray areas auto-selected; trade-off table before each question)
**Areas discussed:** Blind-eval protocol & verdict rules; Physical-anchor readiness & OOD-on-anchors;
Figure assembly under negative results; `run_v2.jl` reproduction scope; Zenodo packaging & licensing
(C-04 handling folded in)

---

## Area 1 — Blind-eval protocol & verdict rules

### Q1: How should the blind evaluation be locked down before the sealed holdout is opened?

| Option | Selected |
|---|---|
| Pre-reg script + frozen expected-outcome table, single shot, disclosed-rerun rule *(recommended)* | ✅ |
| Pre-registered script, no frozen predictions | |
| Structural guard only (existing Phase-8 `sealed_holdout` manifest guard) | |

**Rationale offered:** the project already has five recorded instances of guards satisfied by a
comment rather than by code (Backlog 999.1); pre-registered literal-value tables are the mechanism
that has actually worked here (Phases 5, 13, 15). → **D-01**

### Q2: What should the pre-declared PASS/FAIL rule be?

| Option | Selected |
|---|---|
| Ordinal separation + honest UQ gates; point accuracy and Phase-14 decision outcome reported non-gating *(recommended)* | ✅ |
| Point accuracy is the gate | |
| No gate — descriptive only | |
| You decide | |

**Rationale offered:** the anchors sit at the extremes, where the ρ_true prior-atom artifact lives —
a point-accuracy miss there may measure the prior, not the tool; and the anchors are cross-study
confounded. Gating on the Phase-14 layer was argued against: unshipped, and its own SC3-d came in at
−0.129 vs a 0.50 floor. → **D-02**

### Q3: What happens if the pre-declared gate fails?

| Option | Selected |
|---|---|
| Publish as bounded-scope result; one disclosed re-run only for a demonstrated *harness* defect *(recommended)* | |
| Failure blocks the manuscript | |
| Bounded remediation (any diagnosed cause), disclosed, then re-run | ✅ |

**User chose broader remediation than recommended** — model-side fixes are permitted, not only
harness defects. Consequence recorded explicitly in D-03 rather than left implicit: **the
post-remediation re-run is no longer blind** and must be labelled as such alongside the original
blind result. Exactly one remediation is authorized. → **D-03**

---

## Area 2 — Physical-anchor readiness & OOD-on-anchors

### Q4: How should the `PENDING-FETCH` anchor hashes be closed?

| Option | Selected |
|---|---|
| Full fetch (~6.3 GB positive) + real SHA-256 over full archives + pre-registered seeded subset for eval *(recommended)* | ✅ |
| Full fetch, evaluate everything | |
| Fetch as a `checkpoint:human-verify` | |

**Rationale offered:** Phase 8's own completion note makes the fetch a precondition for opening the
sealed holdout; full-archive digests cover everything, and putting the subset seed in the D-01
pre-registration commit stops it becoming a tuning knob. → **D-04**

### Q5: How should the evaluation handle the OOD flag firing on the physical anchors?

| Option | Selected |
|---|---|
| Two-track (OOD verdict first-class + posteriors labelled "under a raised OOD flag") with a pre-registered anti-vacuity clause *(recommended)* | ✅ |
| OOD fire voids the coloc numbers | |
| Two-track, no anti-vacuity clause | |

**The vacuous-pass hole was the crux:** under D-02 abstention counts as honest, so both anchors
being flagged would pass the gate having tested nothing. The anti-vacuity clause declares that
outcome INCONCLUSIVE. Suppressing the flag was raised and rejected — it would disable claim C4 on
the one dataset where it matters. → **D-05**

### Q6: How should the uncontrolled cross-study confound be handled?

| Option | Selected |
|---|---|
| Both controls — Phase-9 comparators on identical inputs + CBS as tier-labelled condition-matched ladder *(recommended)* | ✅ |
| Comparators only | |
| CBS ladder only | |
| Report the limitation only | |

**Rationale offered:** the comparators attribute any separation (if the classics separate identically,
it's a property of the images, not of v2.0); CBS restores a graded dose-response axis the physical
anchors cannot provide, provided it stays at `simulated-secondary` tier. → **D-06**

---

## Area 3 — Figure assembly under negative results

### Q7: What happens to F5 (registration) and F6 (spatial map)?

| Option | Selected |
|---|---|
| F5 rendered as a negative panel; F6 text-only with spec retained; numbering unchanged *(recommended)* | ✅ |
| Render both as negative panels | |
| Drop both and renumber | |
| Both text-only | |

**Rationale offered:** evidentiary strength differs — F5 has a clean quantified null (RMSE flat in λ,
ratio 1.0003) supporting a real recommendation that the bead anchor exemplifies; F6 has a corrected
sign, three failing arms and an untested exemption. Renumbering was argued against: it breaks the
`fig:` linkage `claims.typ` exists to enforce. → **D-07**

### Q8: How should the figure set be extended, and how is shipped-vs-research-lane provenance marked?

| Option | Selected |
|---|---|
| Add F8 (Phase-15 envelope) + F9 (Phase-14 decision layer incl. the SC3-d miss) + per-panel lane tags *(recommended)* | ✅ |
| F8 only + lane tags | |
| F8 + F9, no lane tags | |
| You decide | |

**Rationale offered:** the operating envelope is the differentiator no competing tool has, and a
one-field YAML lane tag permanently prevents the paper implying `spike/p13`/`spike/p14` results ship
in `src/`. → **D-08**

---

## Area 4 — `run_v2.jl` reproduction scope

### Q9: What should `run_v2.jl` actually reproduce?

| Option | Selected |
|---|---|
| Three tiers (`--smoke` CI-wired / `--figures` / `--full`) with measured wall-clock per tier *(recommended)* | ✅ |
| Two tiers: figures + full | |
| Artifact-restore only | |

**Rationale offered:** counter-based Philox makes `--full` byte-exact rather than merely statistically
similar, so the strongest claim is available and worth making; `--smoke` in the Phase-15 CI gate keeps
the entry point continuously proven; untimed tiers are marketing. → **D-09**

---

## Area 5 — Zenodo packaging & licensing (+ C-04)

### Q10: What should the Zenodo/DOI deposit contain?

| Option | Selected |
|---|---|
| Code + manifests + 82 MB derived artifacts + figure sources + `run_v2.jl`; all third-party imagery by accession + SHA-256 *(recommended)* | ✅ |
| Also bundle the CC-BY anchors (~6.3 GB) | |
| Code only | |

**Measured facts presented:** repo is AGPL-3.0; `artifacts/` is 82 MB; `git ls-files corpus/data` is
empty; CBS is CC-BY-NC-SA-4.0 (bundling it would propagate NC-SA to the collection); both anchors are
CC-BY-4.0. → **D-10**

### Q11: Where should simulator-side draws come from, given C-04?

| Option | Selected |
|---|---|
| Shipped `ood_nulls_8.jld2` first; validation-stream draws via `harness.jl` as fallback; generation at `P12_CHOSEN_PRIOR` prohibited by an executable assertion *(recommended)* | ✅ |
| Validation-stream draws throughout | |
| Shipped nulls only | |

**Rationale offered:** zero new draws means C-04 cannot fire, and scoring against the shipped null
measures the released tool rather than a bespoke re-fit; the executable-assertion requirement comes
straight from the comment-satisfied-guard pattern in Backlog 999.1. → **D-11**

---

## Claude's Discretion (explicitly left open)

- Subset size for D-04's seeded anchor subset (the seed is pre-registered; the size is a planning call)
- F8/F9 panel layout and Typst plumbing — clone the F1 pattern
- Whether F8 splits into F8a/F8b depending on what Phase 15 returns
- `run_v2.jl` CLI surface beyond the three tier flags

## Deferred Ideas Raised

- Bundling the CC-BY anchor imagery as link-rot insurance (rejected for this deposit)
- A condition-matched physical anchor pair (would need wet-lab acquisition)
- Fixing Backlog 999.1 (stays parked)
- Promoting `spike/p13`/`spike/p14` into `src/`
- Option B from the Go/No-Go
