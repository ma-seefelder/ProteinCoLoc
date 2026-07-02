# Phase 10: Manuscript Skeleton and Related-Work Positioning - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-07-02
**Phase:** 10-manuscript-skeleton-and-related-work-positioning
**Mode:** --auto (recommended defaults selected autonomously; no interactive questions)
**Areas discussed:** Manuscript tooling & format, Manuscript location & decoupling, Claim-table structure, Related-work delta (Tapqir/Costes/Manders), Figure specifications, Bibliography, Compile gate (DoD)

---

## Manuscript tooling & format

| Option | Description | Selected |
|--------|-------------|----------|
| Venue-neutral Typst article | Generic single-column Typst scaffold; compiles with installed Typst 0.15.0 | ✓ |
| Specific journal Typst template | Lock a target-journal template now | |
| LaTeX instead of Typst | Use LaTeX toolchain | |

**Auto choice:** Venue-neutral Typst article.
**Notes:** ROADMAP mandates Typst and success criterion #1 requires a compiling Typst skeleton; Typst 0.15.0 is installed. Journal undecided, so template lock-in deferred to avoid rework.

---

## Manuscript location & decoupling

| Option | Description | Selected |
|--------|-------------|----------|
| New top-level `manuscript/` dir | Sibling to spike/, src/, test/; honors decoupling | ✓ |
| Under `spike/` | Co-locate with spike artifacts | |
| Under `.planning/` | Keep with planning docs | |

**Auto choice:** New top-level `manuscript/` directory.
**Notes:** Project decoupling constraint keeps src/ and manuscript pipelines untouched; the v2.0 manuscript is a new artifact and gets its own tree, mirroring how the spike used `spike/`.

---

## Claim-table structure

| Option | Description | Selected |
|--------|-------------|----------|
| Structured claim → axis → phase → figure/exp → status table | The phase spine; diffable + renderable | ✓ |
| Prose-only claims | Narrative claims without a table | |
| No claim tracking | Section stubs only | |

**Auto choice:** Structured claim table as the phase spine.
**Notes:** This is the mechanism that "shapes experiments by the claims they must support"; every row names its substantiating artifact.

---

## Related-work delta (Tapqir / Costes / Manders)

| Option | Description | Selected |
|--------|-------------|----------|
| Comparison matrix + prose (all 4 Tapqir axes) | Methods × capability-axes matrix plus positioning prose | ✓ |
| Prose only | Narrative differentiation | |
| Matrix only | Table without prose | |

**Auto choice:** Comparison matrix + prose covering amortization + SBC + registration-UQ + spatial map.
**Notes:** Mandatory Tapqir differentiation paragraph; honest framing (calibration under the simulator, named OOD blind spot).

---

## Figure specifications

| Option | Description | Selected |
|--------|-------------|----------|
| One spec entry per planned figure (ID/panels/phase/claim/status) | Specs, not rendered figures | ✓ |
| Placeholder image slots only | Empty figure environments | |
| Defer figure planning | No specs this phase | |

**Auto choice:** One structured spec entry per planned figure, anchored to claim rows and owning phases.

---

## Bibliography

| Option | Description | Selected |
|--------|-------------|----------|
| BibTeX `refs.bib`, seeded with key cites | Typst-native; seed Tapqir/Costes/Manders/SBC/NeuralEstimators/self-cite | ✓ |
| Typst Hayagriva YAML | Native Typst bib format | |
| Inline citations | No structured bib | |

**Auto choice:** BibTeX `refs.bib` extending existing `CITATION.bib`.

---

## Compile gate (DoD)

| Option | Description | Selected |
|--------|-------------|----------|
| `typst compile` must succeed + render tables | Runnable/scriptable compile gate | ✓ |
| Manual visual check | No automated gate | |

**Auto choice:** Scriptable compile gate producing a PDF with claim table + related-work matrix rendering.

---

## Claude's Discretion

- Typst version-pinning approach (note installed 0.15.0).
- Exact section ordering and stub wording.
- Whether the claim table is inline `.typ` or a data file rendered by a Typst function.
- Exact figure count/IDs beyond the claim-table-implied minimum.

## Deferred Ideas

- Locking a target journal + applying its Typst template — deferred to near-submission (Phase 16).
- Writing real Results/Methods prose and inserting rendered figures — later data phases + Phase 16.
- CI regression compile gate — Phase 15 CI scope; this phase needs only a local compile check.
