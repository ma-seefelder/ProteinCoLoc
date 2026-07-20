# Phase 8: External Physical Ground-Truth Corpus - Discussion Log

> **Audit trail only.** Not consumed by planning/research/execution agents.
> Decisions are captured in `08-CONTEXT.md` — this log preserves the discussion.

**Date:** 2026-07-20
**Phase:** 08-external-physical-ground-truth-corpus
**Mode:** discuss (interactive, default)
**Areas discussed:** Anchor sourcing, Vendoring policy, CBS ingestion boundary, Data contract + splits

## Area: Anchor sourcing

**Q1 — Physical 100%-coloc anchor source**
- Options: Public archive reuse / Named tandem-fluorophore / New lab acquisition
- **Selected:** Named tandem-fluorophore (DOI/accession-pinned) → D-01

**Q2 — Segregated (negative) anchor source**
- Options: Nuclear-vs-membrane pair / Same-archive matched pair / Distinct-compartment fusion
- **Selected:** "2+3" — distinct-compartment fusion, sourced from the same archive/study as the positive anchor when a matched pair exists → D-02

**Q3 — Acceptance bar for a candidate anchor**
- Options: Biological label + open license / Also require raw 2-channel access / You decide
- **Selected:** You decide → D-03 (Claude's discretion; recommended default = biological label + open license, raw 2-channel retrievability preferred)

## Area: Vendoring policy

**Q4 — Storage relative to repo**
- Options: Manifest + fetch + hash / Vendor bytes in-repo / Hybrid tiny-sample + fetch
- **Selected:** Manifest + fetch + hash (bytes out of git) → D-04

**Q5 — Fetch/verify failure behavior**
- Options: Hard-fail on mismatch + skip-flag on unreachable / Hard-fail on both / You decide
- **Selected:** Hard-fail on mismatch, skip-flag on unreachable (Phase-9 pattern) → D-05

## Area: CBS ingestion boundary

**Q6 — Structural separation from physical anchors**
- Options: Separate tier + explicit label field / Separate file/directory only / Type-level separation
- **Selected:** Separate tier + explicit provenance label field (opt-in to mix) → D-06

**Q7 — Ingestion scope**
- Options: Representative labelled subset / Full benchmark / Manifest-only stub now
- **Selected:** Full benchmark → D-07

## Area: Data contract + splits

**Q8 — Contract/manifest format**
- Options: Reuse content-addressed JLD2+CSV+SHA-256 / Human-first TOML-JSON manifest / You decide
- **Selected:** Reuse content-addressed JLD2+CSV+SHA-256 (Phase-3/9 pattern) → D-08

**Q9 — Split policy for Phase-16 blind eval**
- Options: Sealed blind holdout / Fixed dev-eval split for all / You decide
- **Selected:** Sealed blind holdout (anchors sealed; CBS may carry dev/eval split; manifest field + guard) → D-09

## Deferred ideas

- New wet-lab acquisition — raised, explicitly rejected as out of scope (not deferred).

## Claude's discretion items

- D-03 acceptance predicate (default recommendation recorded).
- Manifest serialization details within the D-08 content-addressed pattern.
