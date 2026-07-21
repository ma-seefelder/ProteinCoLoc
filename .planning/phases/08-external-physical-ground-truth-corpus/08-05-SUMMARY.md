---
phase: 08-external-physical-ground-truth-corpus
plan: 05
subsystem: corpus
tags: [corpus, physical-anchors, ground-truth, sealed-holdout, provenance, non-circularity, human-verify, pending-hash-sentinel, offline-gate, decoupling]
requires:
  - corpus/config.jl consts (MANIFEST_COLUMNS, TIERS, SPLITS, ROLES, CORPUS_DATA_DIR) from 08-01
  - corpus/manifest.jl (validate_anchor_row, build_manifest, validate_manifest, physical_anchors, corpus_default, cbs_benchmark, open_sealed_holdout) from 08-03
  - corpus/fetch.jl (fetch_verified bootstrap mode) from 08-02
  - corpus/cbs.jl (cbs_rows, assign_cbs_split, write_committed_manifest, COMMITTED_MANIFEST_PATH) from 08-04
  - human-verified anchor accessions (Task 1 + Task 2 checkpoint:human-verify)
provides:
  - corpus/anchor_rows.jl (anchor_rows the two human-verified physical anchors; PENDING_FETCH sentinel + is_pending_hash/is_real_sha256; bootstrap_anchor_hashes; finalized_manifest; anchor_hash_notice; write_finalized_manifest; ANCHOR_PROVENANCE_NOTES)
  - corpus/manifest.csv (FINALIZED contract — 30 CBS simulated-secondary rows + 2 human-verified physical-primary sealed_holdout anchors with full provenance)
  - corpus/test/test_anchors.jl wired into the single offline gate
affects:
  - corpus/test/runtests.jl
  - corpus/cbs.jl (committed_manifest docstring marked superseded)
requirements: [SC1]
tech-stack:
  added: []
  patterns:
    - "human-verified accession/license/citation pinned as named consts, re-checked structurally by validate_anchor_row at row assembly (fail-closed error, not a silent skip)"
    - "explicit non-hex PENDING-FETCH sentinel instead of an empty or fabricated sha256, with is_real_sha256/is_pending_hash making 'digest OR sentinel, never in between' testable"
    - "bootstrap fetch defined but NOT invoked by assembly or the test gate — a ~6.3 GB pull stays an explicit human action"
    - "deviations recorded as a machine-readable ANCHOR_PROVENANCE_NOTES const so caveats travel with the contract, not only with the plan summary"
    - "finalized_manifest supersedes the 08-04 placeholder committed_manifest; the superseded function keeps a pointer docstring rather than being silently left ambiguous"
key-files:
  created:
    - corpus/anchor_rows.jl
    - corpus/test/test_anchors.jl
  modified:
    - corpus/manifest.csv
    - corpus/test/runtests.jl
    - corpus/cbs.jl
decisions:
  - "D-01 DEVIATION (human-accepted substitution): no open-licensed, non-environment-quenched tandem-fluorophore dataset exists in any public archive — every deposited tandem-FP dataset is a quenched mRFP/mCherry-EGFP-LC3 autophagy reporter, disqualified by Pitfall 2. Substituted TetraSpeck 100 nm multicolor beads (RegiSTORM sample data, Zenodo 10.5281/zenodo.5509861, CC-BY-4.0): the SAME physical particle emits in both channels, so colocalization is by construction AND state-independent (no pH/quenching failure mode at all) — a strictly stronger physical positive than the rejected alternative, but a deviation from the literal D-01 wording."
  - "D-02 DEVIATION (preference not satisfiable): no qualifying single-study deposit carrying BOTH a same-particle positive and a segregated negative could be verified, so the negative (Light My Cells, S-BIAD1047, CC-BY-4.0) is cross-study/cross-archive. ANCHORS_MATCHED = false. KNOWN LIMITATION recorded in code + manifest provenance: imaging-condition confounds (microscope, objective, exposure, detector, prep) are NOT controlled between the anchors; Phase 16 must report this."
  - "The bootstrap fetch was DELIBERATELY NOT executed: the positive anchor is a ~6.3 GB archive and a bulk download must remain an explicit human decision, not a side effect of automated plan execution. Both anchors therefore carry the explicit PENDING-FETCH sentinel. No digest was fabricated (T-08-16). bootstrap_anchor_hashes() is implemented and callable when online and authorized."
  - "sha256 uses a deliberately NON-hex sentinel string rather than empty-string: an empty hash is indistinguishable from 'not yet filled' CBS rows and could be silently accepted by a downstream integrity check, whereas 'PENDING-FETCH' can never be mistaken for a 64-hex digest and is asserted as such by is_real_sha256."
  - "anchor_rows() is parameterized on sha256/bytes rather than mutating a global, so the post-bootstrap path (real digests) and the current pending path are the SAME code path and both are covered by the offline gate."
  - "Anchor ids are descriptive (pos-tetraspeck-01 / neg-lightmycells-01) rather than the generic 08-04 placeholders (pos-anchor-01 / neg-anchor-01); the test asserts the placeholder ids and any ',PENDING,' provenance are GONE from the committed manifest."
metrics:
  duration: ~20 min
  tasks: 3 (2 human-verify checkpoints satisfied by the human + 1 auto)
  files: 5
  completed: 2026-07-21
---

# Phase 8 Plan 05: Physical Ground-Truth Anchors Finalized Summary

The ONE non-autonomous, non-circular deliverable of Phase 8. Both physical anchors were
human-verified against the D-03 acceptance predicate and are now pinned in the committed corpus
contract as `physical-primary` / `sealed_holdout` rows (D-09) with full provenance, honest
recorded deviations, and an explicit pending-hash sentinel in place of any fabricated digest.

## Confirmed Anchors

**POSITIVE — TetraSpeck 100 nm multicolor fiducial beads** (RegiSTORM v1.0.0 sample data)
- Accession: Zenodo DOI `10.5281/zenodo.5509861`; URL `https://zenodo.org/api/records/5509861/files/Sample%20Data.zip/content`
- License: CC-BY-4.0 (verified on the Zenodo record `license` field)
- Citation: Karlsson et al. 2023, "RegiSTORM: channel registration for multi-color STORM",
  BMC Bioinformatics, DOI 10.1186/s12859-023-05320-1
- Truth (PHYSICS, not a score): a single physical bead emits in EVERY colour channel ⇒ the two
  channels are colocalized BY CONSTRUCTION and **state-independent**. Not an environment-quenched
  tandem (Pitfall 2 explicitly avoided).
- Format: multi-channel STORM `.tif` frame stacks (~30,000 frames/channel) + ThunderSTORM `.csv`
  reconstructions; a mean/max projection is required to yield the 2-channel intensity image
  (conversion step owned by Phase 16). ~6.3 GB.

**NEGATIVE — "Light My Cells"** (Bright Field to Fluorescence Imaging Challenge, France-BioImaging / ISBI 2024)
- Accession: BioImage Archive `S-BIAD1047`; URL `https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/047/S-BIAD1047`
- License: CC-BY-4.0 (verified in study metadata)
- Truth (BIOLOGY, not a score): nucleus (DNA stain) vs mitochondria (cytoplasmic,
  nucleus-excluded) are spatially disjoint by compartment identity from the acquisition design.
- Format: OME-TIFF, REMBI-compliant metadata. Fields MUST be filtered to those carrying BOTH the
  nucleus and mitochondria channels (the challenge guarantees only "at least one" target/field) —
  filter step owned by Phase 16.
- `matched? = NO` (different study, different archive).

## What Was Built

**Task 1 / Task 2 — `checkpoint:human-verify` (SATISFIED by the human)**
Both anchors confirmed against the D-03 predicate: biological/physical truth label from the
source publication (not a computed score), open redistributable license, resolving accession,
retrievable two-channel-convertible data.

**Task 3 — anchor rows + finalized manifest + tests (commit 8eb9147)**
- `corpus/anchor_rows.jl`:
  - `anchor_rows(; sha256_pos, sha256_neg, bytes_pos, bytes_neg)` — the two rows
    (`pos-tetraspeck-01` positive/coloc/1.0, `neg-lightmycells-01` negative/segregated/0.0), both
    `physical-primary` + `sealed_holdout`, each ASSERTED to pass `validate_anchor_row` (fail-closed
    `error`) before return.
  - `PENDING_FETCH` = `"PENDING-FETCH"` sentinel + `is_pending_hash` / `is_real_sha256` predicates.
  - `bootstrap_anchor_hashes(; datadir)` — the real `fetch_verified(url, dest, nothing)` bootstrap
    per anchor; **not invoked** by assembly or the gate (see Deviations).
  - `finalized_manifest(; kwargs...)` = `assign_cbs_split(cbs_rows())` (30) + `anchor_rows()` (2),
    `build_manifest` + `validate_manifest` ⇒ 32 rows; supersedes 08-04's placeholder
    `committed_manifest()`.
  - `anchor_hash_notice(df)` / `write_finalized_manifest(df)` — writes the repo-tracked
    `corpus/manifest.csv` (D-08 header, atomic `.tmp→mv`) and emits the pending-hash notice.
  - `ANCHOR_PROVENANCE_NOTES` — machine-readable D-01/D-02/hash caveats.
- `corpus/manifest.csv` — FINALIZED: 30 CBS rows unchanged + the two anchors with real
  accession/URL/citation/license and `sha256=PENDING-FETCH`. No `PENDING` provenance remains.
- `corpus/test/test_anchors.jl` (79 asserts) wired into `corpus/test/runtests.jl`.
- `corpus/cbs.jl` — `committed_manifest` docstring marked SUPERSEDED, pointing at
  `finalized_manifest()` / `write_finalized_manifest()`.

## Verification Results

- `julia --project=. corpus/test/runtests.jl` → **217/217 pass** in ~26s
  (fixture 12 + hashing 17 + fetch 20 + manifest 34 + load 9 + CBS 46 + **anchors 79**).
  Precisely: the gate is **offline-green / network-optional**, not literally zero-network — the
  inherited `fetch_cbs_smoke()` (08-04) attempts one real GET and asserts
  `status ∈ (:skipped, :present)`. On this run the source was unreachable, so the observed run was
  network-free and `corpus/data/` was empty afterwards. Nothing added by 08-05 touches the network:
  `bootstrap_anchor_hashes()` is never called by the gate.
- `git ls-files corpus/data` → **empty** (D-04 / Pitfall 3).
- `git status --porcelain src/` → **empty**; no `src/` file touched (hard CLAUDE.md constraint).
- Default accessors verified sealed: `physical_anchors(df)` → 0 rows, `corpus_default(df)` carries
  no `physical-primary`/`sealed_holdout` row, `cbs_benchmark(df)` carries no anchor;
  `open_sealed_holdout(df; reason="…")` → exactly the 2 anchors and throws on an empty reason (D-09).
- Committed CSV round-trip: re-read via `CSV.read(...; comment="#")` re-validates; both physical
  rows carry non-empty accession + license + citation and a hash that is either 64-hex or the
  explicit sentinel.

## Deviations from Plan

1. **D-01 substitution (human-accepted, recorded, not papered over).** The plan/decision specified
   a cellular *tandem-fluorophore construct*. No open-licensed, non-environment-quenched tandem-FP
   dataset could be verified in any public archive — the deposited tandem-FP datasets are quenched
   autophagy reporters (mRFP/mCherry-EGFP-LC3), which the acceptance bar disqualifies (Pitfall 2:
   GFP is quenched in the acidic autolysosome, so the deposited channels are NOT 100% colocalized).
   The human accepted a **multicolor-bead** substitution instead. Recorded in the file header, in
   `ANCHOR_PROVENANCE_NOTES.d01`, and asserted by the test.
2. **D-02 preference unmet.** No qualifying single-study deposit with both a same-particle positive
   and a segregated negative could be verified, so the anchors are **cross-study / cross-archive**
   (`ANCHORS_MATCHED = false`). KNOWN LIMITATION: imaging-condition confounds are NOT controlled
   between positive and negative. Recorded in `ANCHOR_PROVENANCE_NOTES.d02`; Phase 16 must report it.
3. **Bootstrap fetch deliberately skipped.** The plan's Task 3 calls `fetch_verified(...)` in
   bootstrap mode. The positive anchor is a **~6.3 GB** archive, so the bulk download was left as an
   explicit human decision rather than a side effect of automated execution. Both anchors therefore
   record the plan-specified offline outcome: `sha256 = "PENDING-FETCH"`. **No digest was
   fabricated** (T-08-16). This is the plan's documented offline branch, taken by policy rather than
   by connectivity.

## Known Stubs

- **`sha256 = "PENDING-FETCH"` on both anchors.** `bootstrap_anchor_hashes()` is implemented and
  ready; it must be run **when online and authorized**, and the sentinel replaced by the returned
  digests + byte sizes (via `finalized_manifest(; sha256_pos=…, sha256_neg=…, bytes_pos=…,
  bytes_neg=…)` + `write_finalized_manifest`) **BEFORE Phase 16 opens the sealed holdout**.
  `anchor_hash_notice` emits this loudly on every manifest write while a sentinel is present.
- **`bytes = 0` on both anchors** for the same reason (recorded at bootstrap, never estimated).
- **Conversion/filtering steps owned by Phase 16**: the positive needs a frame-stack mean/max
  projection; the negative needs field filtering to those carrying BOTH nucleus and mitochondria.
  Both are documented in the row `format`/`channels` comments; no converter is fabricated here.

## Threat Flags

None outstanding. Plan `<threat_model>` coverage:
- **T-08-14** (circular/wrong positive): human-verify + `validate_anchor_row` fail-closed; truth is
  physical (one bead, both channels), asserted 1.0/0.0, never a Pearson/Manders score.
- **T-08-15** (environment-quenched tandem mis-pinned): explicitly rejected; the substitution is
  state-independent by construction.
- **T-08-16** (fabricated hash): mitigated by the non-hex `PENDING-FETCH` sentinel + `is_real_sha256`
  assertions; bootstrap not run, nothing invented.
- **T-08-07** (sealed anchor readable pre-Phase-16): asserted — no default accessor reaches them.
- **T-08-17** (missing provenance): accession + citation + license non-empty and non-`PENDING`,
  asserted per row and on the re-read CSV.
- **T-08-SC** (supply chain): no new dependency.

## Self-Check: PASSED

- FOUND: corpus/anchor_rows.jl
- FOUND: corpus/test/test_anchors.jl
- FOUND: corpus/manifest.csv (finalized, 30 CBS + 2 anchors)
- FOUND: `include("test_anchors.jl")` wired into corpus/test/runtests.jl
- FOUND commit: 8eb9147 (Task 3)
- `git ls-files corpus/data` empty; `src/` untouched
