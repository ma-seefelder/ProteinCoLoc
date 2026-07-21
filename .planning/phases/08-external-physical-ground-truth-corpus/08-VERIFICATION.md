---
phase: 08-external-physical-ground-truth-corpus
verified: 2026-07-21T09:20:00Z
status: human_needed
score: 11/11 must-haves verified
has_blocking_gaps: false
overrides_applied: 3
overrides:
  - must_have: "SC1 — physical anchors archived with provenance + SHA-256 content hashes"
    reason: "Both anchors carry the explicit non-hex sentinel sha256=\"PENDING-FETCH\" because the bootstrap fetch was deliberately NOT run (the positive anchor is a ~6.3 GB archive; a bulk download must stay an explicit human authorization). No digest was fabricated (T-08-16). bootstrap_anchor_hashes() is implemented, callable, and the sentinel is asserted never to pass as a digest. Sentinel MUST be replaced before Phase 16 opens the sealed holdout."
    accepted_by: "Manuel Seefelder"
    accepted_at: "2026-07-21T09:08:00Z"
  - must_have: "D-01 — positive anchor is a cellular tandem-fluorophore construct"
    reason: "Human-accepted substitution: every deposited tandem-FP dataset is a quenched mRFP/mCherry-EGFP-LC3 autophagy reporter, disqualified by the acceptance bar (GFP quenched in the acidic autolysosome ⇒ channels are NOT 100% colocalized). Substituted TetraSpeck 100 nm multicolor beads (Zenodo 10.5281/zenodo.5509861, CC-BY-4.0): the SAME physical particle emits in both channels ⇒ colocalization by construction AND state-independent."
    accepted_by: "Manuel Seefelder"
    accepted_at: "2026-07-21T09:08:00Z"
  - must_have: "D-02 — prefer a matched segregated negative from the SAME study/archive"
    reason: "Preference not satisfiable: no qualifying single-study deposit carrying both a same-particle positive and a segregated negative could be verified. Anchors are cross-study/cross-archive (ANCHORS_MATCHED=false); imaging-condition confounds are NOT controlled and are recorded as an explicit KNOWN LIMITATION in code + manifest provenance. Phase 16 must report it."
    accepted_by: "Manuel Seefelder"
    accepted_at: "2026-07-21T09:08:00Z"
human_verification:
  - test: "Authorize and run corpus bootstrap: `include(\"corpus/anchor_rows.jl\"); bootstrap_anchor_hashes()` while online, then rewrite the contract via `write_finalized_manifest(finalized_manifest(; sha256_pos=…, sha256_neg=…, bytes_pos=…, bytes_neg=…))`"
    expected: "Both anchors' PENDING-FETCH sentinel replaced by real 64-hex digests + non-zero byte sizes; anchor_hash_notice(df) returns \"\"; `git ls-files corpus/data` stays empty. MUST happen BEFORE Phase 16 opens the sealed holdout."
    why_human: "A ~6.3 GB download from Zenodo plus a BioImage Archive pull is an explicit human authorization decision and cannot be run inside verification; the digest must never be fabricated."
  - test: "Re-confirm at Phase 16 write-up that the D-01 bead substitution and the D-02 cross-study confound are reported in the manuscript's limitations"
    expected: "Both ANCHOR_PROVENANCE_NOTES.d01 / .d02 caveats appear in the Phase-16 external-validation reporting"
    why_human: "Editorial/scientific reporting judgement, not codebase-checkable"
  - test: "CBS bytes: run the CBS fetch path online once the per-file zip filenames are resolved from the downloads page"
    expected: "CBS rows get real sha256 + byte sizes; safe_extract wired to a concrete zip reader; offline gate remains green"
    why_human: "Requires network + scraping the CBS downloads page for exact filenames (deliberately not invented — RESEARCH A6)"
---

# Phase 8: External Physical Ground-Truth Corpus — Verification Report

**Phase Goal:** Assemble a non-circular validation corpus whose truth does not come from the model's own simulator, so v2.0's calibration can be checked against external reality rather than mere self-consistency.
**Verified:** 2026-07-21
**Status:** human_needed (all automated must-haves verified; 3 human actions remain, chiefly the authorized hash bootstrap before Phase 16)
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | **SC1** — ≥1 physical 100%-coloc anchor + ≥1 segregated anchor pinned with provenance + content hashes | ✓ VERIFIED (override: hash sentinel, D-01/D-02 deviations) | `corpus/anchor_rows.jl:186-228` emits exactly `pos-tetraspeck-01` (role=positive, truth_label=coloc, gt=1.0, Zenodo `10.5281/zenodo.5509861`, CC-BY-4.0, RegiSTORM citation) and `neg-lightmycells-01` (role=negative, truth_label=segregated, gt=0.0, BioImage Archive `S-BIAD1047`, CC-BY-4.0). Both rows are present verbatim in the committed `corpus/manifest.csv:35-36`. `sha256` = `PENDING-FETCH` sentinel on both (accepted override). |
| 2 | **SC2** — CBS ingested and explicitly flagged simulated/secondary, never conflated with the physical anchors | ✓ VERIFIED | `corpus/cbs.jl:81-105` tags every row `tier="simulated-secondary"`, `role="benchmark"`, `truth_label="simulated-degree"`; 30 such rows in `manifest.csv:5-34`. Structural non-conflation: `corpus/manifest.jl:113-172` has **no** default both-tier accessor — `cbs_benchmark` returns only simulated rows, `corpus_default` structurally drops sealed rows (with an `@assert`), and the ONLY union path is the explicitly-named `mixed_tier`. Runtime check: `count(physical-primary in cbs_benchmark(df))=0`, `corpus_default(df)` physical count = 0. |
| 3 | **SC3** — versioned validation data contract + manifest checked in | ✓ VERIFIED | `corpus/manifest.csv` is git-tracked (`git ls-files corpus/`), header-stamped with `MANIFEST_SCHEMA_VERSION = 1` + `CORPUS_MASTER_SEED = 0x0000000000c05eed`, 14 fixed columns, 33 data lines (1 header + 30 CBS + 2 anchors). **Regeneration check I ran independently: `write_committed_manifest(finalized_manifest())` output is BYTE-IDENTICAL to the committed file** — the contract is code-derived, not hand-edited. |
| 4 | Single offline command runs the gate green | ✓ VERIFIED | I ran `julia --project=. corpus/test/runtests.jl` → **217 / 217 pass, 26.0 s**, testsets: fixture 7, conversion 5, hashing 17, fetch 20, manifest 34, load 9, CBS 46, anchors 79. |
| 5 | **D-04** — no image bytes in git | ✓ VERIFIED | `git ls-files corpus/data` → **empty**. `.gitignore:444` `corpus/data/*` with scoped negations `!corpus/manifest.csv`, `!corpus/test/fixtures/`. `git check-ignore` confirms `corpus/data/cbs_smoke.zip` is ignored. `corpus/data/` is empty on disk after the full test run. |
| 6 | **D-05** — mismatch = hard abort, unreachable = skip-with-flag | ✓ VERIFIED | `corpus/fetch.jl:87-112`: the ONLY `try/catch` wraps the download and returns `(status=:skipped, reason="unreachable: …")`; the digest `!=` check at L106 sits **outside** any try/catch and `error(...)`s, removing the `.tmp`. `InterruptException` rethrown. Timeout bound via `DOWNLOAD_TIMEOUT_S=120`. Covered by 20 offline asserts (`test_fetch.jl:76,84`). |
| 7 | **D-09** — physical anchors are a structurally sealed holdout | ✓ VERIFIED | `validate_manifest` (`manifest.jl:104-109`) errors on any `physical-primary` row whose split ≠ `sealed_holdout`. `physical_anchors(df)` → 0 rows by construction; `open_sealed_holdout(df; reason)` is the sole path and asserts a non-empty audit reason. Independently re-ran: `nrow(physical_anchors(df))=0`, `corpus_default` physical count `=0`. |
| 8 | **D-03** — anchor acceptance predicate is fail-closed | ✓ VERIFIED | `validate_anchor_row` (`manifest.jl:187-197`) requires a biological label ∈ {coloc, segregated}, non-empty license, non-empty accession, and an ASSERTED 1.0/0.0 ground truth (a computed score is rejected). `anchor_rows` **errors** if either row fails (`anchor_rows.jl:221-226`) — not a silent skip. |
| 9 | Two-channel TIFF converts to `MultiChannelImage` and reduces through the UNCHANGED `src/` summary path | ✓ VERIFIED | `corpus/load.jl:60-93` calls the frozen `MultiChannelImage(name, paths, channels)` ctor and `patch`/`correlation` at 8×8 via read-only `include` of `src/LoadImages.jl` + `src/colocalization.jl`; `has_signal` guards the ≤15-valid-pixel floor. Exercised green by the 9 `load:` asserts + 5 fixture-conversion asserts against the committed `fixture_ch1/ch2.tif`. |
| 10 | `src/` provably untouched (hard CLAUDE.md decoupling constraint) | ✓ VERIFIED | Union of files touched by ALL Phase-8 commits (`08-01`…`08-05` + phase-08 docs) = `.gitignore`, `.planning/**`, `corpus/**` only. **No `src/` file, no root `Project.toml`/`Manifest.toml`.** `git status --porcelain src/` → empty. (The `src/` diff vs `main` is Phase-7 productionization, the phase licensed to edit `src/`.) |
| 11 | **D-07** full CBS enumeration + **D-09/seeded** reproducible split | ✓ VERIFIED | `CBS_PAIRS` × `CBS_DEGREES` = 3 × 10 = 30 rows (`cbs.jl:58-105`), all labelled degrees 0–90 % as archive-labelled fractions (never computed). `assign_cbs_split` keys a `Philox4x` on `CORPUS_MASTER_SEED ⊻ CBS_SPLIT_SALT` and the 1-based row index; I re-ran it twice → **identical assignment** (dev 14 / eval 16, sealed 2, no CBS row ever `sealed_holdout`). Zip-slip rejection (`_safe_target_path`) rejects escaping entries before any write. |

**Score:** 11/11 truths verified (3 accepted overrides applied, see frontmatter)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `corpus/config.jl` | Pre-declared consts before any manifest exists | ✓ VERIFIED | `CORPUS_MASTER_SEED`, `MANIFEST_SCHEMA_VERSION`, `TIERS`, `SPLITS`, `ROLES`, `MANIFEST_COLUMNS`, `CBS_SPLIT_SALT`, `DOWNLOAD_TIMEOUT_S`, `MIN_PATCH_VALID_PIXELS`, `CORPUS_DATA_DIR` — all present, all consumed downstream |
| `corpus/hash.jl` | SHA-256 content hash / manifest hash / canonical config | ✓ VERIFIED | 17 passing asserts; `file_sha256` consumed by `fetch.jl:96`, `manifest_hash` by `manifest.jl:244` |
| `corpus/fetch.jl` | `fetch_verified` with the two distinct D-05 paths + timeout | ✓ VERIFIED | 113 lines, structural asymmetry, network seam isolated in `_DOWNLOAD_HOOK` Ref so tests stub it |
| `corpus/manifest.jl` | Schema, tier accessors, split guard, D-03 predicate, content-addressed atomic writer | ✓ VERIFIED | 276 lines; atomic `.tmp → @assert → mv(force=true)` for both CSV and JLD2 twin |
| `corpus/load.jl` | Two-channel TIFF → `MultiChannelImage`, read-only `src/` reuse | ✓ VERIFIED | read-only guarded includes; no `src/` mutation |
| `corpus/cbs.jl` | CBS enumeration, seeded split, zip-slip-safe extraction, RGB split | ✓ VERIFIED | 30 rows; `committed_manifest` correctly marked SUPERSEDED with a pointer to `finalized_manifest()` |
| `corpus/anchor_rows.jl` | The two human-verified anchors + sentinel + bootstrap + finalized manifest | ✓ VERIFIED | 317 lines; `PENDING_FETCH` / `is_pending_hash` / `is_real_sha256`; `ANCHOR_PROVENANCE_NOTES` machine-readable deviations |
| `corpus/manifest.csv` | The committed versioned contract | ✓ VERIFIED | 32 rows, byte-identical to code regeneration, no `pos-anchor-01` placeholder, no `,PENDING,` provenance |
| `corpus/test/runtests.jl` + 6 sub-suites + 2 fixtures | Single offline gate | ✓ VERIFIED | all six sub-suites wired via `include`; 217 asserts |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `runtests.jl` | all sub-suites | `include(...)` | ✓ WIRED | 6 includes, lines 104-109 |
| `fetch.jl` | `hash.jl` | `file_sha256` post-download | ✓ WIRED | `fetch.jl:96` |
| `manifest.jl` | `hash.jl` | `manifest_hash` content-addressing | ✓ WIRED | `manifest.jl:244` |
| `load.jl` | `src/LoadImages.jl`, `src/colocalization.jl` | read-only `include` | ✓ WIRED | `load.jl:44-45`, exercised by 9 passing asserts |
| `cbs.jl` | `manifest.jl` | `build_manifest` + `validate_manifest` | ✓ WIRED | `cbs.jl:231-234` |
| `cbs.jl` | `fetch.jl` | `fetch_verified` (offline-skip smoke) | ✓ WIRED | `cbs.jl:277-279` |
| `anchor_rows.jl` | `manifest.jl` | `validate_anchor_row` fail-closed at assembly | ✓ WIRED | `anchor_rows.jl:221-226` |
| `anchor_rows.jl` | `corpus/manifest.csv` | `write_finalized_manifest` → `write_committed_manifest` | ✓ WIRED | Regenerated output byte-identical to the committed file |
| `.gitignore` | `corpus/data` | scoped negation | ✓ WIRED | `.gitignore:444-449`, `git check-ignore` confirms |

### Data-Flow Trace (Level 4)

| Artifact | Data variable | Source | Produces real data | Status |
|----------|---------------|--------|--------------------|--------|
| `corpus/manifest.csv` | 32 manifest rows | `finalized_manifest()` = `assign_cbs_split(cbs_rows())` + `anchor_rows()` | Yes — regenerated output is byte-identical | ✓ FLOWING |
| anchor rows | `sha256` | `bootstrap_anchor_hashes()` (deliberately not run) | **No — explicit `PENDING-FETCH` sentinel** | ⚠️ STATIC (accepted override; never fabricated, loudly flagged by `anchor_hash_notice`) |
| CBS rows | `sha256`, `bytes` | first authorized fetch (not run) | No — empty / 0 | ⚠️ STATIC (by D-04/D-05 design; deferred to Phase 16 fetch) |
| `summary_grid(img)` | 8×8 correlation grid | frozen `src/` `patch`/`correlation` on the committed fixture | Yes — finite grid, `has_signal` true in tests | ✓ FLOWING |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| Offline gate green | `julia --project=. corpus/test/runtests.jl` | 217/217 pass, 26.0 s | ✓ PASS |
| No data bytes in git | `git ls-files corpus/data` | empty | ✓ PASS |
| `src/` untouched | `git status --porcelain src/` + per-commit `--name-only` over all 08-* commits | empty / no `src/` path | ✓ PASS |
| Contract is code-derived | regenerate `write_committed_manifest(finalized_manifest())` into a temp path and compare | `regenerated == committed : true` | ✓ PASS |
| Seeded split reproducible | `assign_cbs_split(cbs_rows())` twice | identical (dev 14 / eval 16) | ✓ PASS |
| Sealed holdout unreachable by default | `nrow(physical_anchors(df))`, `corpus_default(df)` tier count | 0 and 0 | ✓ PASS |
| Pending-hash notice emitted | `anchor_hash_notice(df)` | "2 physical anchor(s) carry the 'PENDING-FETCH' sentinel …" | ✓ PASS |

### Probe Execution

No `scripts/*/tests/probe-*.sh` probes are declared or implied for this phase; `corpus/test/runtests.jl` is the phase's single declared gate and was executed above.

### Requirements Coverage

| Requirement | Source | Description | Status | Evidence |
|-------------|--------|-------------|--------|----------|
| SC1 | ROADMAP (local) | Physical coloc + segregated anchors with provenance + content hashes | ✓ SATISFIED (with hash-sentinel override) | Truth 1 |
| SC2 | ROADMAP (local) | CBS ingested, flagged simulated-secondary, never conflated | ✓ SATISFIED | Truth 2 |
| SC3 | ROADMAP (local) | Versioned data contract + manifest checked in | ✓ SATISFIED | Truth 3 |

No orphaned requirements: `REQUIREMENTS.md` marks Phase-8 requirements as local ROADMAP success criteria (TBD in the global file), and all three are covered.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `corpus/cbs.jl` | 93 | `# verified prefix (filename TODO at fetch)` | ⚠️ Warning | Honest marker, not hidden debt: the exact per-zip filenames are deliberately not invented (RESEARCH A6). Consequence: CBS `source_url` is a prefix only, so a per-row CBS fetch is not yet actionable — deferred to the Phase-16 fetch. No TBD/FIXME/XXX blocker markers anywhere in `corpus/`. |
| `corpus/cbs.jl` | 204-235 | `anchor_placeholder_rows` / `committed_manifest` retained | ℹ️ Info | Correctly marked SUPERSEDED with a pointer to `finalized_manifest()`; the committed CSV provably comes from the finalized path (byte-identical regeneration) and the test asserts `pos-anchor-01` and `,PENDING,` are absent. Dead-ish but documented, not misleading. |
| `corpus/manifest.csv` | 5-34 | CBS rows carry empty `sha256`, `bytes=0` | ⚠️ Warning | By design (hashes recorded only on the first authorized fetch), but SC2's "ingested" is contract-level enumeration, not byte-level ingestion. Phase 16 must fetch + fill. |
| `corpus/test/test_cbs.jl` | 132-133 | `fetch_cbs_smoke()` with `@test r.status in (:skipped, :present)` | ⚠️ Warning | The gate is **offline-green / network-optional**, not literally zero-network: it attempts one real GET to the CBS prefix and tolerates either outcome. On this run the source was unreachable (`corpus/data/` is empty afterwards), so the observed run was effectively network-free. The SUMMARY's "ZERO network" wording is slightly stronger than the code. |
| `.planning/ROADMAP.md` | Phase-8 block | `08-05-PLAN.md` checkbox unchecked, Phase 8 listed `4/5 In Progress` | ⚠️ Warning | Stale phase bookkeeping — 08-05 is committed (`8eb9147`) and its SUMMARY exists. Orchestrator should close the checkbox at phase close. No code impact. |

### Human Verification Required

**1. Authorize the anchor hash bootstrap (blocks Phase 16, not Phase 8)**
- **Test:** while online, `include("corpus/anchor_rows.jl"); r = bootstrap_anchor_hashes()`, then `write_finalized_manifest(finalized_manifest(; sha256_pos=r.pos.sha256, sha256_neg=r.neg.sha256, bytes_pos=r.pos.bytes, bytes_neg=r.neg.bytes))`
- **Expected:** both `PENDING-FETCH` sentinels replaced by real 64-hex digests + non-zero byte counts; `anchor_hash_notice(df) == ""`; `git ls-files corpus/data` still empty; the offline gate still 217/217.
- **Why human:** a ~6.3 GB Zenodo pull plus a BioImage Archive pull is an explicit authorization decision; fabricating a digest is forbidden (T-08-16).

**2. Carry the D-01/D-02 caveats into the Phase-16 write-up**
- **Test:** confirm `ANCHOR_PROVENANCE_NOTES.d01` (bead substitution for the tandem-FP construct) and `.d02` (`ANCHORS_MATCHED=false`, uncontrolled cross-study imaging confounds) are reported in the external-validation section.
- **Why human:** scientific reporting judgement.

**3. Resolve CBS per-file zip names and fetch CBS bytes**
- **Test:** scrape the CBS downloads page for the exact zip filenames, wire a concrete zip reader into `safe_extract(...; entries=…)`, fetch, and fill CBS `sha256`/`bytes`.
- **Why human:** requires network and a page-structure judgement the plan deliberately refused to guess.

### Gaps Summary

**No blocking gaps.** All three ROADMAP success criteria are true in the codebase, verified independently of the SUMMARY narrative: I re-ran the gate (217/217, 26 s), re-derived the committed manifest from code (byte-identical), re-ran the sealed-holdout and seeded-split checks in a fresh Julia process, and confirmed the two hard constraints (`git ls-files corpus/data` empty; no Phase-8 commit touches `src/` or the root project files).

Three deviations are recorded as **accepted overrides**, not failures:
1. `sha256 = "PENDING-FETCH"` on both anchors — the bootstrap fetch was deliberately not run; the sentinel is non-hex by construction, asserted never to pass `is_real_sha256`, and loudly surfaced by `anchor_hash_notice` on every manifest write. This is the one item that **must** close before Phase 16 opens the sealed holdout.
2. D-01 satisfied by a human-accepted substitution (TetraSpeck multicolor beads instead of a cellular tandem-FP construct) — physically a stronger, state-independent positive; the deviation is recorded in code, in `ANCHOR_PROVENANCE_NOTES.d01`, and asserted by the test suite.
3. D-02's matched-same-study preference is unmet (`ANCHORS_MATCHED = false`); the cross-study imaging confound is an explicitly recorded limitation.

The residual warnings (CBS rows without digests, prefix-only CBS URLs, the network-optional CBS smoke, the stale ROADMAP checkbox) are non-goal-blocking and are either by design or bookkeeping.

---

_Verified: 2026-07-21_
_Verifier: Claude (gsd-verifier)_
