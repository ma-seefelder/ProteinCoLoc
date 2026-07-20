# Phase 8: External Physical Ground-Truth Corpus - Research

**Researched:** 2026-07-20
**Domain:** Reproducible external-data acquisition + content-addressed data contracts (no model code)
**Confidence:** HIGH on in-repo mechanics (verified in code) / MEDIUM-LOW on specific dataset accessions (see Assumptions Log)

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01:** Physical 100%-coloc anchor = a **named published tandem-fluorophore construct** dataset, pinned by DOI/accession as THE canonical positive anchor. Truth asserted by biology + citation, never by any computed colocalization score (this is what keeps it non-circular).
- **D-02:** Segregated (negative) anchor = a **distinct-compartment two-fluorophore construct** (non-overlapping compartments, e.g. nuclear-vs-membrane, mito-vs-nucleus), pinned by DOI/accession. **Prefer a matched segregated construct from the same archive/study as the positive anchor** when one exists (controls imaging/instrument confounds); otherwise a separate named dataset. New lab acquisition was explicitly rejected.
- **D-04:** Anchor (and CBS) image **bytes stay OUT of git**. Deliverable = manifest (accession/DOI/URL) + **seeded fetch script** + **SHA-256 content hashes** verifying integrity on download.
- **D-05:** content-hash **MISMATCH = HARD error (abort)**; source **UNREACHABLE = graceful skip-with-flag** (offline/CI stays green), mirroring the Phase-9 Tapqir bridge.
- **D-06:** CBS lives in its **own manifest tier/namespace** with a **mandatory provenance field** (`simulated-secondary` for CBS vs `physical-primary` for anchors). Any code path mixing tiers must **opt in explicitly** — structural separation, not a comment. Satisfies SC2.
- **D-07:** Ingest the **FULL CBS benchmark corpus** (all labelled colocalization degrees) via the same manifest+fetch+hash mechanism, flagged `simulated-secondary`.
- **D-08:** Versioned contract **reuses the Phase-3 / Phase-9 content-addressed pattern**: human-readable manifest (CSV via DataFrames/CSV.jl) carrying provenance/license/tier/hash/split, a **SHA-256 content hash over data-defining bytes** for naming+verification, JLD2 for cached arrays.
- **D-09:** Split policy — **physical anchors are a SEALED blind holdout** (opened only at Phase-16); CBS may carry a seeded dev/eval split. Enforced by a manifest **split field + a guard** (structural, not a convention).

### Claude's Discretion
- **D-03:** Concrete **acceptance predicate** for qualifying an anchor is left to the researcher/planner. Recommended default (adopted below, §Acceptance Predicate): require (a) an unambiguous biological coloc/segregation label from the source publication AND (b) an open/redistributable license; **prefer** raw two-channel data retrievable and convertible to a `MultiChannelImage` so `src/` summary functions apply unchanged.
- Manifest serialization details (exact column set, TOML vs CSV for the human layer) within the D-08 pattern.

### Deferred Ideas (OUT OF SCOPE)
- None. New wet-lab acquisition was raised and **explicitly rejected** (not deferred).
</user_constraints>

<phase_requirements>
## Phase Requirements

No REQ-IDs are mapped for Phase 8 (ROADMAP: "Requirements: TBD"). The three ROADMAP success criteria are the requirements this phase must make TRUE:

| ID (local) | Description | Research Support |
|------------|-------------|------------------|
| SC1 | ≥1 physical 100%-coloc anchor (single-protein-two-channel / tandem fluorophore) AND ≥1 segregated anchor (nuclear-vs-membrane) archived with provenance + content hashes | §Candidate Datasets (positive/negative), §Fetch + Integrity Mechanics, §Manifest Schema |
| SC2 | CBS ingested and explicitly flagged simulated/secondary, never conflated with physical anchors | §CBS (verified), §Manifest Schema (tier column + structural guard), Pattern 3 |
| SC3 | Versioned validation data contract + manifest (schema, provenance, split policy) checked in | §Manifest Schema, §Fetch + Integrity Mechanics, D-08/D-09 reuse map |
</phase_requirements>

## Summary

This phase is a **data-acquisition + data-contract** phase, not a modelling phase. The deliverable is machinery — a seeded fetch script, a SHA-256 integrity layer, and a versioned CSV manifest — plus a small set of **externally-truthed** image anchors whose colocalization ground truth comes from biology/citation, never from a computed score. The entire integrity/manifest/skip-with-flag/atomic-write toolchain **already exists in this repo** (Phase-3 `src/amortized/datagen.jl`, Phase-9 `spike/comparator/`) and should be reused near-verbatim. No new external Julia packages are required: `SHA` and `Downloads` are stdlib; `CSV`, `DataFrames`, `JLD2`, `Images`, `Random123` are already direct deps.

The one genuinely open, load-bearing question is **which concrete datasets** to pin. One dataset is fully verified: the **Colocalization Benchmark Source (CBS)** — a real, citable, downloadable simulated benchmark (colocalization-benchmark.com; Zinchuk & Grossenbacher-Zinchuk; CC BY-NC-SA 4.0; lossless TIFF; labelled 0–90% coloc degrees). It maps cleanly to the `simulated-secondary` tier (D-06/D-07). The **physical anchors** could not be pinned to a single verified accession in this session; I give ranked candidate leads with exact search paths and honest confidence, and recommend the planner gate each anchor selection behind a `checkpoint:human-verify` task. This is low-risk because the manifest+fetch+hash architecture is **dataset-agnostic** — the machinery can be built and tested (against CBS + a synthetic fixture) before the physical-anchor URLs/hashes are finalized.

**Primary recommendation:** Build a decoupled `corpus/` module (NOT under `src/` model code) that reuses the Phase-3/Phase-9 content-addressing verbatim: a tiered CSV manifest (`tier ∈ {physical-primary, simulated-secondary}`, `split ∈ {sealed_holdout, dev, eval}`), a `Downloads.download → SHA-256 verify → atomic mv` fetch script (hash-mismatch = hard error, unreachable = `:skipped` flag), and a `MultiChannelImage` loader for the anchors. Pin CBS now; treat each physical-anchor accession as a human-verified checkpoint, defaulting the positive anchor to a **tandem fluorescent-protein construct** (D-01) with **TetraSpeck multicolor beads** as the strongest fallback-or-second positive.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| External data acquisition (download) | Fetch script (stdlib `Downloads`) | — | Network I/O isolated in one seam so it can be stubbed/skipped offline |
| Integrity verification | Content-hash layer (`SHA`) | Fetch script | SHA-256 over fetched bytes; mismatch aborts (D-05) |
| Provenance + tier + split contract | CSV manifest (`DataFrames`/`CSV`) | JLD2 cache | Human-readable + machine-checkable; tier/split are structural guards (D-06/D-09) |
| Physical/simulated separation | Manifest `tier` field + guard | Directory namespace | Structural (not comment) opt-in to mix tiers (D-06, SC2) |
| Sealed-holdout enforcement | Manifest `split` field + guard | Directory namespace | Anchors inaccessible without explicit opt-in until Phase 16 (D-09) |
| Image → summary conversion | `src/LoadImages.jl` loaders (read-only) + `MultiChannelImage` | `src/colocalization.jl` `patch()`/`correlation()` (read-only) | Anchors must reduce through the UNCHANGED summary path (D-03 preference) |
| Reproducible CBS subset/split selection | `Random123` Philox (seeded) | — | Bit-reproducible split assignment, mirrors Phase-3/9 |

**Decoupling note (architectural decision for the planner):** This phase must NOT edit `src/` *model* code and is "parallelizable now, alongside Phases 5–7". The Phase-7 **ship-gate machinery** precedent lives in `test/gate/` (validation infra, not model code); the Phase-9 comparator lives in `spike/comparator/`. Recommend the corpus machinery live in a **new decoupled directory** — `corpus/` at repo root (or `test/corpus/`) — with its own fetch/manifest/loader modules and a test file wired into the existing test runner. It may `include("../src/LoadImages.jl")` / `include("../src/colocalization.jl")` **read-only** for the `MultiChannelImage` conversion (the Phase-1/2 `include()` read-only contract), exactly as spike phases reach `src/` without editing it. Confirm the exact location in planning. `[ASSUMED — architectural placement]`

## Standard Stack

### Core (all already available — NO new deps required)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `SHA` | stdlib | SHA-256 content hash over fetched bytes + manifest data-defining bytes | Already used in `spike/comparator/table.jl` (`SHA.SHA256_CTX`/`update!`/`digest!`/`bytes2hex`); never hand-roll crypto [VERIFIED: repo grep] |
| `Downloads` | stdlib | Fetch remote files over HTTPS (TLS via libcurl) | Julia stdlib; no add, Windows-clean; not yet used in repo but the canonical fetch primitive [VERIFIED: stdlib] |
| `CSV` | current (root+spike dep) | Human-readable manifest table writer/reader | Already a direct dep, used for Phase-9 `table.csv` [VERIFIED: Project.toml] |
| `DataFrames` | current (root+spike dep) | Tidy manifest as a `DataFrame` | Already a direct dep, Phase-9 pattern [VERIFIED: Project.toml] |
| `JLD2` | current (root+spike dep) | Cache decoded arrays / manifest meta (`jldsave` + reopen-integrity) | Repo-wide durable-cache idiom (Phase-3/5/7/9) [VERIFIED: Project.toml] |
| `Images` | current (root+spike dep) | TIFF load + RGB→channel split → `Matrix{Float64}` | Already the loader in `src/LoadImages.jl:load_tiff` (`Images.load` / `Images.Gray`) [VERIFIED: src grep] |
| `Random123` | current (spike dep) | Philox-seeded CBS subset/split selection, bit-reproducible | Repo-wide seeding idiom (Phase-3/9 `MASTER_SEED`) [VERIFIED: spike Project.toml] |

### Supporting (add only if a format demands it)
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `TiffImages.jl` | current | Lazy/BigTIFF/OME-TIFF and robust multi-page TIFF decode | ONLY if `Images.load` chokes on an anchor's OME-TIFF or BigTIFF; verify need on the actual downloaded bytes before adding [ASSUMED] |
| `p7zip`/`ZipFile.jl` or stdlib | current | Unzip CBS `.zip` set archives | CBS ships `.zip` bundles; prefer decompressing to a temp dir then hashing individual TIFFs. Verify whether `Downloads` + a zip reader suffices before adding a dep [ASSUMED] |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| stdlib `Downloads` | `HTTP.jl` | HTTP.jl adds a heavy dep for no benefit here; `Downloads` covers HTTPS GET with redirects. Use HTTP.jl only if an archive needs custom headers/auth (none expected) |
| `Images.load` | `TiffImages.jl` | `Images.load` already works for standard TIFF (repo `load_tiff` proves it); TiffImages only for OME/Big/lazy |
| CSV manifest | TOML manifest | D-08 says "human-readable"; CSV matches the Phase-9 `table.csv` precedent and DataFrames round-trips. TOML is fine per D-03 discretion but CSV is the established pattern — recommend CSV |

**Installation:** None. All required packages are stdlib or existing direct deps. If a supporting dep is later proven necessary, add it to the **corpus module's own environment / the relevant Project.toml only after verifying the format actually requires it** — do not preemptively bloat the dependency surface.

## Package Legitimacy Audit

No new external packages are introduced by the recommended stack (all stdlib or pre-existing verified direct deps). slopcheck not required for this phase.

| Package | Registry | Status | Disposition |
|---------|----------|--------|-------------|
| SHA, Downloads | Julia stdlib | Ships with Julia | Approved |
| CSV, DataFrames, JLD2, Images, Random123 | General registry | Already resolved in root/spike Manifest (Phases 3–9) | Approved (reuse) |
| TiffImages, ZipFile | General registry | Only-if-needed; verify against actual downloaded bytes first | Deferred — gate behind a `checkpoint:human-verify` if adopted |

## Candidate Datasets (load-bearing)

> **Honesty flag:** Only CBS is fully verified with a canonical source in this session. Every *physical-anchor* pick below is `[ASSUMED]` and must be confirmed (exact accession + license + channel layout + a downloadable URL) before the fetch script pins its URL and hash. Do NOT let the planner hard-code an unverified accession — gate each with `checkpoint:human-verify`. The architecture is dataset-agnostic, so machinery can be built and tested against CBS + a synthetic fixture first.

### Acceptance Predicate (D-03 — RECOMMENDED, adopt as the phase's rule)
A candidate qualifies as an anchor iff **ALL** of:
1. **Biological truth label** — the source publication states, from construct/targeting biology (not from a coloc measurement), that the two channels are either fully colocalized (positive) or targeted to non-overlapping compartments (negative).
2. **Open/redistributable license** — CC-BY / CC-BY-NC / CC0 / equivalent, compatible with citing + re-fetching (bytes stay out of git per D-04, so NC is acceptable for validation use; record the exact license in the manifest).
3. **Retrievable + convertible (PREFERRED, near-required)** — raw two-channel image data downloadable via a stable URL/accession and convertible to a `MultiChannelImage` (two `Matrix{Float64}` channels) so `patch()`/`correlation()` apply unchanged.
4. **Provenance pinned** — a DOI or repository accession (IDR `idr####`, BioImage Archive `S-BIAD####`/`S-BSST####`, Zenodo DOI, or CellImageLibrary ID).

A candidate meeting (1)+(2)+(4) but only offering pre-composited/derived images (fails 3) is a **fallback** anchor: usable for a coarse check but flag it can't feed the unchanged summary path.

### CBS — Colocalization Benchmark Source (simulated-secondary) — VERIFIED
- **What:** Computer-simulated fluorescence image sets with **precisely known, pre-defined colocalization degrees** — the exact "labelled-simulated" corpus the ROADMAP names. `[VERIFIED: colocalization-benchmark.com]`
- **Provenance / citation:** Zinchuk V, Wu Y, Grossenbacher-Zinchuk O. Primary methods citation: Zinchuk & Grossenbacher-Zinchuk (2014), *Current Protocols in Cell Biology* 62:4.19.1–4.19.14; foundational: Zinchuk et al. (2013) *Scientific Reports* 3:1365; Zinchuk et al. (2011) *Nature Protocols*. `[CITED: colocalization-benchmark.com/about]`
- **Ground-truth levels:** colocalization degrees in **10% increments: 0,10,20,…,90%** (per the downloads page). `[VERIFIED: colocalization-benchmark.com/downloads]`  *(Note: the exact top level — whether a 100% set exists and the precise numeric definition of each "degree" — should be confirmed on download; some CBS documentation frames degrees via Pearson Rr / overlap coefficient. Record the per-image labelled value verbatim from the archive.)* `[ASSUMED — exact numeric semantics]`
- **Channel pairs / packaging:** 3 sets — **Red-Green, Red-Blue, Green-Blue** — each 10 images; individual images ≈1.7 MB, full set ≈15.7 MB, packaged as `.zip`. `[VERIFIED: downloads page]`
- **Format:** Lossless TIFF. `[VERIFIED]`
- **License:** **Creative Commons Attribution-NonCommercial-ShareAlike 4.0** (CC BY-NC-SA 4.0). `[VERIFIED: downloads page]` — NC is fine for this validation use; bytes stay out of git (D-04). Record `CC-BY-NC-SA-4.0` in the manifest `license` column and honor attribution in any figure/manuscript.
- **Download URLs:** direct links under `colocalization-benchmark.com/download/[filename].zip` (per-set and per-image zips). Recommended citation string per site: "Image set CBS001…–CBS010… from the Colocalization Benchmark Source (www.colocalization-benchmark.com)". Confirm exact zip filenames by scraping the downloads page at fetch-script build time. `[VERIFIED — pattern; exact filenames to confirm on fetch]`
- **D-03 verdict:** Meets (1) simulated-but-labelled, (2) open, (4) citable. **Tier = `simulated-secondary` (D-06/D-07). Ingest the FULL set across all degrees.**

### Positive physical anchor (physical-primary, 100%-coloc) — CANDIDATES
Ranked by fit to D-01 ("named published tandem-fluorophore construct") + the acceptance predicate:

1. **Tandem fluorescent-protein construct** (LITERAL D-01, PREFERRED) — a single polypeptide bearing two spectrally-distinct fluorophores (e.g. a tandem `EGFP-mCherry`/`dTomato` fusion, or a "tandem fluorescent timer" tFT). Because both fluorophores are on the same molecule, the two channels image the same physical species → 100% colocalization by construction, truth from biology+citation. **Caveat to avoid:** environment-quenched tandems (e.g. autophagy `mCherry-GFP-LC3`, where lysosomal GFP quenching deliberately *breaks* coloc) are NOT valid 100%-anchors — pick a cytosolic/non-quenched tandem. `[ASSUMED — needs a pinned accession]`
   - **Search path:** IDR (`idr.openmicroscopy.org`) study search "tandem fluorescent"/"dual-labelled fusion"; BioImage Archive (`ebi.ac.uk/bioimage-archive`) query "tandem fluorescent protein" filtered to two-channel; Zenodo "tandem fluorescent timer imaging". Verify: two raw channels, open license, stable accession.

2. **TetraSpeck multicolor beads** (STRONGEST NON-CIRCULARITY, looser D-01 wording) — 100 nm beads containing **4 fluorophores co-located in one physical particle** → every channel is 100% colocalized by bead chemistry; the field's standard multicolor fiducial/registration/coloc positive control. Truth is even more unambiguous than a molecular tandem. **But** beads are not a "fluorophore *construct*" in the molecular-biology sense of D-01, so confirm in discuss-phase whether they qualify as THE positive anchor or serve as a **second/confirmatory** positive. `[ASSUMED — needs a pinned image dataset]`
   - **Search path:** Zenodo "TetraSpeck beads" (registration/SMLM benchmark datasets frequently deposit bead stacks); IDR registration studies; single-molecule-localization-microscopy benchmark repositories. TetraSpeck = Invitrogen/Thermo T7280.

3. **CBS highest-coloc level as positive** — **REJECT as the physical positive.** It is `simulated-secondary`; using it as the physical anchor would reintroduce exactly the circularity this phase exists to eliminate. It stays in the CBS tier only.

### Negative (segregated) physical anchor (physical-primary) — CANDIDATES
Ranked; **prefer a matched pair from the SAME study/archive as the positive** (D-02):

1. **Matched dual-compartment two-fluorophore construct** (LITERAL D-02, PREFERRED) — the same cells expressing a **nuclear-targeted FP** (e.g. `H2B-GFP` or an NLS-FP) + a **membrane/mito-targeted FP** (e.g. `Lyn-mCherry`, `mito-DsRed`) → physically non-overlapping compartments → segregation truth by targeting biology. If the positive-anchor study also provides such a construct, use it (controls instrument/imaging confounds). `[ASSUMED — needs a pinned accession]`
   - **Search path:** same IDR/BioImage Archive studies as the positive; filter for co-expressed organelle-targeted FP pairs.

2. **Broad Bioimage Benchmark Collection (BBBC)** two-channel cell sets (e.g. DAPI nucleus + F-actin/cytoskeleton channel) — many BBBC sets are openly licensed with stable `BBBC0xx` IDs; nucleus-vs-cytoskeleton is a clean segregation. `[ASSUMED — needs the specific BBBC id + license confirmed]`
   - **Search path:** `bbbc.broadinstitute.org` — browse two-channel fluorescence sets; confirm license (many are CC0/attribution) + channel layout.

3. **FluoCells / BPAE prepared-slide standard images** (DAPI nucleus + MitoTracker mito + phalloidin actin) — widely distributed, nucleus-vs-mito largely segregated, BUT these are often **vendor sample images with restrictive/unclear licensing** and truth is "largely" not "perfectly" segregated. Fallback only; flag the license risk. `[ASSUMED]`

**Segregated-anchor note:** perfect physical segregation is more available than perfect physical colocalization; the harder anchor to source is the positive. If a single matched study yielding BOTH a tandem positive and a segregated negative can be found, that is the ideal (D-02 preference) — make that the top search priority.

## Fetch + Integrity Mechanics (reuse in-repo patterns — VERIFIED)

All primitives exist in the repo; the fetch script composes them. Cite these exact files/functions to the planner:

### Content hash (SHA-256) — reuse `spike/comparator/table.jl:comparator_hash`
```julia
# Source: spike/comparator/table.jl (VERIFIED) — the repo's SHA-256 content-hash idiom.
using SHA
# Hash of a downloaded FILE's bytes (naming + integrity verification, D-04/D-05/D-08):
file_sha256(path) = bytes2hex(SHA.sha256(read(path)))          # or stream in chunks for large files

# Hash over the manifest's data-defining bytes (fixed order, NUL-separated) — mirrors comparator_hash:
function manifest_hash(src_files::Vector{String}, canonical_config::String)
    ctx = SHA.SHA256_CTX()
    for f in src_files
        SHA.update!(ctx, read(f)); SHA.update!(ctx, UInt8[0x00])   # inter-file separator
    end
    SHA.update!(ctx, Vector{UInt8}(canonical_config))
    return bytes2hex(SHA.digest!(ctx))
end
```

### Fetch + verify + atomic commit — compose `Downloads` with the Phase-3/9 atomic-write idiom
```julia
# Sources: src/amortized/datagen.jl:write_shard + spike/comparator/table.jl:write_table (VERIFIED atomic idiom)
#          .tmp -> integrity check -> mv(force=true).  Skip-with-flag mirrors spike/comparator/tapqir_bridge.jl.
using Downloads
function fetch_verified(url::AbstractString, dest::AbstractString, expected_sha256::AbstractString)
    tmp = dest * ".tmp"
    try
        Downloads.download(url, tmp)                              # HTTPS GET (TLS via libcurl)
    catch e
        e isa InterruptException && rethrow()
        isfile(tmp) && rm(tmp; force=true)
        return (status=:skipped, reason="unreachable: $(sprint(showerror, e))")  # D-05 graceful skip → CI stays green
    end
    got = bytes2hex(SHA.sha256(read(tmp)))
    if got != lowercase(expected_sha256)                         # D-05 HARD error on mismatch (tamper/corruption)
        rm(tmp; force=true)
        error("content-hash MISMATCH for $url: expected $(expected_sha256), got $got")   # ABORT
    end
    mv(tmp, dest; force=true)                                    # atomic commit (filesystem rename)
    return (status=:present, sha256=got)
end
```

Key contract points for the planner:
- **Hash MISMATCH ⇒ `error(...)` / abort** (D-05, tamper/corruption is not recoverable). Matches `datagen.jl:open_or_invalidate` which `error`s on stale-cache hash mismatch.
- **Unreachable ⇒ return `(status=:skipped, reason=...)`**, never throw — the offline/CI path stays green. Directly mirrors `spike/comparator/tapqir_bridge.jl:tapqir_anchor` (`:skipped|:passed|:failed` NamedTuple; probe-first, catch-all-to-`false`).
- **Atomic `.tmp → integrity → mv(force=true)`** for every persisted artifact (never a half-written file). Verified in BOTH `datagen.jl:write_shard/write_meta` and `table.jl:write_table`.
- **Bootstrapping the expected hash:** on the FIRST authorized fetch of a newly-pinned dataset (no prior hash), compute-and-record the hash into the manifest (a `checkpoint:human-verify` moment — a human confirms the source is the intended one before its hash becomes the trusted baseline). Thereafter mismatch = abort. Make this two-mode behavior explicit and structural, not accidental.
- **Seeded selection/split:** use `Random123` Philox with a corpus `MASTER_SEED` (Phase-9 `config.jl` pattern) for any CBS subset or dev/eval split assignment, so the split is bit-reproducible.

## Manifest Schema (D-08)

Recommended CSV columns (one row per image OR per channel-file; pick per-image with a `channels` descriptor — confirm granularity in planning). Written as a `DataFrame` → `CSV.write`, mirroring `spike/comparator/table.jl:write_table`, with a pre-registration comment header quoting the corpus `MASTER_SEED` and schema version verbatim.

| Column | Type | Purpose | Notes |
|--------|------|---------|-------|
| `anchor_id` | String | Stable local identifier | e.g. `pos-tandem-01`, `neg-nucmem-01`, `cbs-RG-050` |
| `tier` | String | **MANDATORY provenance (D-06)** | `physical-primary` \| `simulated-secondary`. The structural separation field (SC2) |
| `role` | String | Anchor role | `positive` \| `negative` \| `benchmark` |
| `truth_label` | String | Biological truth | `coloc` \| `segregated` \| `simulated-degree` |
| `coloc_ground_truth` | Float/String | Numeric/known label | CBS: the labelled degree (e.g. `0.50`); anchors: `1.0`/`0.0` asserted by biology (NOT a measured score) |
| `accession` | String | DOI / IDR / S-BIAD / Zenodo / CBS id | Provenance pin (D-01/D-02/D-04) |
| `source_url` | String | Direct fetch URL | Consumed by the seeded fetch script |
| `citation` | String | Bibliographic / DOI citation | For manuscript + non-circularity audit |
| `license` | String | Redistribution terms | e.g. `CC-BY-NC-SA-4.0` (CBS), `CC-BY-4.0`, `CC0` |
| `format` | String | On-disk format | `tiff` \| `ome-tiff` \| `zip` |
| `channels` | String | Channel layout | e.g. `R,G` or `nucleus:GFP;membrane:mCherry` — drives `MultiChannelImage` assembly |
| `sha256` | String | **SHA-256 content hash** (D-04/D-05/D-08) | Over the data-defining bytes; naming + integrity |
| `split` | String | **Split policy (D-09)** | `sealed_holdout` (anchors — opened only Phase 16) \| `dev` \| `eval` (CBS) |
| `bytes` | Int | Expected file size | Cheap secondary integrity signal |
| `fetch_status` | String | Runtime (in a separate fetch report, not the committed manifest) | `present` \| `skipped_unreachable` |

**Structural guards (not comments — D-06/D-09):**
- **Tier guard (SC2, D-06):** the loader API returns anchors and CBS through **separate typed accessors** (e.g. `physical_anchors(manifest)` vs `cbs_benchmark(manifest)`); a function that returns BOTH tiers must be a distinct, explicitly-named `mixed_tier(...)` entry the caller must opt into. Do not expose a default "give me everything" that silently blends tiers. Mirror `datagen.jl`'s "no `standardize_all` symbol exists so leakage is impossible by construction" discipline.
- **Split guard (D-09):** `sealed_holdout` rows are unreachable from any default corpus accessor. Provide an explicit `open_sealed_holdout(; reason)` that is called ONLY by Phase-16 code; a unit test asserts the default corpus iterator never yields a `sealed_holdout` row (structural anti-snooping, mirrors Phase-3 disjoint-index holdout and Phase-5 `VAL_MASTER_SEED` separation).
- **Manifest content hash:** the committed manifest is itself content-addressed (a `manifest_hash` over the manifest-defining sources + a canonical config), so any edit to the contract is detectable — reuse `datagen.jl:cache_hash`/`subhashes` for which-source-changed diagnosis.

## Image Format → MultiChannelImage Conversion (VERIFIED)

The `src/` constructor and loaders are already suited to two-channel anchors — reach them **read-only** (no `src/` edits).

**Constructor (verified `src/LoadImages.jl:126`):**
```julia
MultiChannelImage(
    data::Vector{Matrix{T}},              # one Matrix{Float64} per channel; T <: Union{Missing,Float64}
    channels::Vector{S},                  # e.g. ["ch1","ch2"]   (S <: AbstractString)
    name::S,
    path::Vector{S},                      # source paths, one per channel
    pixel_size::Tuple{I,I},               # I <: Int (must be Tuple{Int,Int})
    otsu_threshold::Vector{F},            # one per channel; F <: AbstractFloat
)
# asserts length(data)==length(channels)==length(otsu_threshold)
# A convenience constructor exists that fills pixel_size = size(data[1]) and
# otsu_threshold = Images.otsu_threshold.(data)  (see 02-01-PLAN interfaces).
```

**TIFF → channels (verified `src/LoadImages.jl:214`):**
```julia
# Source: src/LoadImages.jl (VERIFIED)
load_tiff(path) = Float64.(Images.Gray.(Images.load(path)))     # grayscale single-channel TIFF → Matrix{Float64}
```
- **Two separate single-channel TIFFs** (typical anchor layout): `load_tiff` each → `data = [ch1, ch2]`.
- **One RGB/composite TIFF** (typical CBS layout — R/G/B planes): split planes with `Images.channelview` / `Images.red.(img)`,`green`,`blue`, each `→ Float64` matrix, then pick the two relevant channels per the CBS pair (RG/RB/GB).
- **Summary path is UNCHANGED (D-03, SIM-03 precedent):** once a `MultiChannelImage` exists, `patch(img, 8)` + `correlation(x,y)` from `src/colocalization.jl` apply verbatim. Watch the documented `correlation` semantics: `_exclude_zero` drops `0.0`/`NaN`/`missing`, and a patch with `length ≤ 15` valid pixels → `missing`. Anchor images must have enough non-zero signal per 8×8 patch — verify on the real anchor, not just the fixture.
- **TiffImages.jl** only if `Images.load` fails on an OME-TIFF/BigTIFF anchor — verify on actual bytes before adding.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SHA-256 hashing | A custom digest loop | `SHA` stdlib (`SHA.sha256` / `SHA256_CTX`) as in `table.jl` | Crypto correctness; already the repo idiom |
| HTTPS download | A curl subprocess / socket code | `Downloads.download` (stdlib) | TLS, redirects, Windows-clean, no dep |
| Atomic durable write | A bare `write`/`jldsave` to the final path | `.tmp → integrity check → mv(force=true)` from `datagen.jl:write_shard` / `table.jl:write_table` | Crash never leaves a half-written artifact |
| Skip-with-flag on failure | ad-hoc `try/catch` that returns bool | The `(status=:skipped\|:present\|..., reason=...)` NamedTuple from `tapqir_bridge.jl` | Uniform, testable graceful-degradation contract (D-05) |
| Reproducible split/subset | `rand()` / `Random.seed!` | `Random123` Philox keyed by a corpus `MASTER_SEED` | Bit-reproducible across runs/threads (Phase-3/9) |
| Content-addressed cache dir | Timestamped/ad-hoc dir names | `cache_hash(config)`-named dir from `datagen.jl` | Same `(config, bytes)` ⇒ same dir; stale auto-invalidates |
| TIFF decode | Manual byte parsing | `Images.load` (repo `load_tiff`) | Already proven on this stack |

**Key insight:** every mechanic this phase needs already exists and is battle-tested in Phases 3, 5, 7, and 9. The phase's real work is **curation + contract**, not new infrastructure. The largest risk is picking a wrong/unlicensed/circular dataset — mitigated by the acceptance predicate + human-verify checkpoints, NOT by more code.

## Runtime State Inventory

> This is an acquisition/contract phase, not a rename/refactor. Included for completeness since it touches on-disk artifacts.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | Downloaded image bytes land OUTSIDE git (D-04) in a gitignored corpus data dir; only manifest + fetch script + hashes are committed | Ensure `.gitignore` excludes the data dir but tracks the manifest/script (mirror the Phase-2 `.gitignore` scoped-negation for `spike/figures/plausibility.png`) |
| Live service config | None — external archives (CBS, IDR, BioImage Archive) are read-only public HTTP; no service config stored | None — verified: fetch is one-way GET |
| OS-registered state | None | None |
| Secrets/env vars | None — public archives need no auth; do NOT introduce API keys | None — verified: all candidate sources are anonymous public download |
| Build artifacts | JLD2 decoded-array caches (if used) are content-addressed and regenerable; stale auto-invalidates via `cache_hash` | None beyond honoring the content-hash guard |

## Common Pitfalls

### Pitfall 1: Circular positive anchor
**What goes wrong:** using CBS (or any coloc-score-selected image) as the "physical" positive anchor.
**Why:** CBS is simulated/secondary; its truth is definitional, not physical — using it defeats the phase's entire purpose (non-circularity).
**Avoid:** enforce the tier guard (D-06) structurally; the positive/negative physical anchors MUST be `physical-primary` with truth from construct biology + citation.
**Warning sign:** a `physical-primary` row whose `coloc_ground_truth` was derived from a Pearson/Manders computation rather than asserted from biology.

### Pitfall 2: Environment-quenched tandem construct
**What goes wrong:** picking `mCherry-GFP-LC3` (autophagy tandem) as the 100%-coloc positive.
**Why:** that construct is DESIGNED so GFP quenches in acidic lysosomes → the channels deliberately DE-colocalize; ground truth is not 100%.
**Avoid:** choose a cytosolic/non-quenched tandem where both fluorophores report the same molecule everywhere.
**Warning sign:** the source paper uses the construct to *measure* a process via changing coloc.

### Pitfall 3: License leakage into git
**What goes wrong:** committing CBS/anchor image bytes (CC-BY-NC-SA / vendor-restricted) into the repo.
**Why:** violates D-04 and possibly the license.
**Avoid:** bytes-out-of-git is structural — the fetch dir is gitignored; only manifest+hashes+script are tracked. A test asserts no image bytes are staged.
**Warning sign:** `git status` shows `.tif`/`.zip` under the corpus data dir.

### Pitfall 4: Hash-mismatch swallowed as skip
**What goes wrong:** treating a hash mismatch like an unreachable source (skip-with-flag).
**Why:** D-05 is asymmetric — unreachable = skip (green), but **mismatch = tamper/corruption = HARD abort**. Conflating them hides integrity failures.
**Avoid:** two distinct code paths — `catch` on `Downloads.download` → `:skipped`; explicit `!=` on the digest → `error(...)`.
**Warning sign:** a single `try/catch` wrapping both the download and the hash check.

### Pitfall 5: Sealed holdout accidentally readable
**What goes wrong:** the default corpus iterator yields the physical anchors, so model-dev code can touch them before Phase 16.
**Why:** D-09 anti-snooping requires structural inaccessibility, not a naming convention.
**Avoid:** anchors only reachable via `open_sealed_holdout(; reason)`; a unit test asserts the default accessor never returns a `sealed_holdout` row.
**Warning sign:** any Phase 8–15 code reading a `physical-primary` image.

### Pitfall 6: `correlation()` returns all-missing on an anchor
**What goes wrong:** an anchor image reduces to an all-`missing` 8×8 correlation grid.
**Why:** `_exclude_zero` drops zeros and patches with ≤15 valid pixels → `missing`; a sparse/low-signal or wrongly-scaled anchor can starve every patch.
**Avoid:** verify on the REAL anchor bytes (not just the synthetic fixture) that patches retain enough non-zero signal; check intensity scaling/background.
**Warning sign:** `mean(skipmissing(ρ))` is `NaN` on the anchor.

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|------------------|--------|
| Validate coloc tools only on self-simulated data | Pair simulated benchmark (CBS) with **externally-truthed physical anchors** | Non-circular calibration check — the phase's whole point |
| Vendor bytes vendored into repos | **Manifest + content-hash + on-demand fetch** (D-04) | Small, license-safe, reproducible repos |
| Trust downloaded data implicitly | **SHA-256 verify, hard-fail on mismatch** | Tamper/corruption detection |

**Deprecated/outdated for this phase:** none. All tooling is current and in-repo.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Specific positive tandem-FP anchor accession | Candidate Datasets → Positive | HIGH — the phase's central anchor; must be human-verified before pinning URL/hash. Mitigated: architecture is dataset-agnostic |
| A2 | TetraSpeck bead image dataset exists openly on Zenodo/IDR | Positive #2 | MEDIUM — strong conceptually; specific archived set unverified |
| A3 | Matched segregated anchor available in the same study as the positive | Negative #1 | MEDIUM — D-02 *prefers* matched; separate datasets are an acceptable fallback |
| A4 | BBBC / FluoCells specifics (exact id + license + channels) | Negative #2/#3 | MEDIUM — leads not pinned; FluoCells license may be restrictive |
| A5 | CBS exact numeric degree semantics + whether a 100% set exists | CBS | LOW — CBS verified as a source; only the precise per-image label value needs confirming on download |
| A6 | CBS zip filenames / exact direct URLs | CBS | LOW — URL *pattern* verified; exact filenames confirmable by scraping the downloads page at fetch time |
| A7 | Corpus machinery placement (`corpus/` vs `test/corpus/`) | Architectural Responsibility Map | LOW — planner decides; must stay decoupled from `src/` model code |
| A8 | `Images.load` suffices for all anchor/CBS TIFFs (no TiffImages/zip dep) | Standard Stack | LOW — verify on actual bytes; add supporting dep only if proven necessary |

**These A-items should drive `checkpoint:human-verify` tasks in the plan** — especially A1–A4 (dataset selection) before the fetch script pins any URL or hash.

## Open Questions

1. **Which exact physical-anchor datasets?** (A1–A4)
   - Known: CBS fully verified for the simulated tier; the acceptance predicate is defined.
   - Unclear: no single verified accession for the physical positive/negative in this session.
   - Recommendation: planner adds a `checkpoint:human-verify` task per anchor; the user/discuss-phase confirms accession+license+channel layout; the fetch script then computes+records the trusted baseline hash on first authorized fetch.

2. **Corpus machinery location** (A7) — `corpus/` at root vs `test/corpus/`. Recommendation: decoupled top-level `corpus/` reaching `src/LoadImages.jl`+`src/colocalization.jl` read-only via `include()`; confirm in planning.

3. **CBS 100% level / degree definition** (A5) — confirm on download whether a 100% set exists and record each image's labelled value verbatim.

4. **Does the positive anchor qualify under D-01's literal "construct" wording if beads are chosen?** — discuss-phase decision: molecular tandem as THE positive (beads as confirmatory second) vs beads as THE positive.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| `SHA` (stdlib) | Content hashing | ✓ | ships with Julia | — |
| `Downloads` (stdlib) | Fetch | ✓ | ships with Julia | — |
| `CSV`/`DataFrames`/`JLD2`/`Images`/`Random123` | Manifest, cache, conversion, seeding | ✓ | resolved in root+spike Manifest | — |
| Network access to CBS/IDR/BioImage Archive | Live fetch of real bytes | ✗ at CI/offline | — | **skip-with-flag (D-05)**: fetch returns `:skipped`, tests run against a committed synthetic fixture → CI green |
| `TiffImages.jl` / zip reader | Only if an anchor is OME/BigTIFF or CBS zip needs in-Julia extraction | ✗ (not added) | — | Add only if `Images.load`/stdlib prove insufficient on real bytes |

**Missing with fallback:** network — fully covered by the D-05 skip-with-flag design; the corpus tests MUST NOT require network to pass.

## Validation Architecture

> `nyquist_validation: true` in `.planning/config.json` → this section is required.

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Julia stdlib `Test` (`@testset`/`@test`), the repo's existing pattern (`spike/test/runtests.jl`, `src`/`test` runners) |
| Config file | none — Julia `Test`; new corpus tests wired into a `runtests.jl` (Wave 0) |
| Quick run command | `julia --project=<corpus-env> corpus/test/runtests.jl` (offline: no-network subset) |
| Full suite command | same, with network available (adds one live small-CBS-fetch smoke) |

### Phase Requirements → Test Map
| Req | Behavior | Test Type | Automated Command | File Exists? |
|-----|----------|-----------|-------------------|-------------|
| SC3 | Manifest schema valid (required columns, tier ∈ {physical-primary,simulated-secondary}, split ∈ {sealed_holdout,dev,eval}) | unit | `runtests.jl` manifest testset | ❌ Wave 0 |
| SC2/D-06 | Tier guard: default accessor never blends tiers; mixing requires explicit `mixed_tier` opt-in | unit | tier-separation testset | ❌ Wave 0 |
| D-09 | Split guard: default corpus iterator never yields a `sealed_holdout` row; anchors only via `open_sealed_holdout` | unit | holdout-guard testset | ❌ Wave 0 |
| D-05 | Hash MISMATCH = hard error (inject a flipped byte on a fixture → expect `error`) | unit | integrity testset | ❌ Wave 0 |
| D-05 | Unreachable = skip-with-flag (bad URL → `(status=:skipped)`, green) | unit | skip-with-flag testset | ❌ Wave 0 |
| SC1/D-04 | Content hash deterministic + names artifact (same bytes ⇒ same sha256) | unit | hashing testset | ❌ Wave 0 |
| SC1/D-03 | Fixture two-channel TIFF → `MultiChannelImage` → `patch(·,8)`/`correlation` apply unchanged, finite result | unit | conversion testset (synthetic fixture, no network) | ❌ Wave 0 |
| D-04 | No image bytes staged in git (corpus data dir gitignored) | unit | gitignore/staging assertion | ❌ Wave 0 |
| SC2/D-07 | CBS live fetch smoke (one small set) verifies + tier-tags — SKIPS cleanly offline | integration | full-suite only | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** offline unit subset (manifest, guards, hashing, skip-with-flag, fixture conversion) — no network, always green.
- **Per wave merge:** full offline suite.
- **Phase gate:** offline suite green; live CBS-fetch smoke green when online (skip-with-flag otherwise) before `/gsd:verify-work`.

### Wave 0 Gaps
- [ ] `corpus/test/runtests.jl` — the testsets above
- [ ] A committed **synthetic two-channel TIFF fixture** (generated in-test or tiny tracked file) so conversion + hashing + integrity tests run with ZERO network
- [ ] `.gitignore` scoped so the corpus **data** dir is excluded but the manifest/fetch-script/tests are tracked (mirror the Phase-2 scoped-negation idiom)
- [ ] Confirm the corpus module's environment (reuse spike/root deps — no new adds expected)

## Security Domain

> `security_enforcement` not set in config (absent = enabled). Scope here is minimal (read-only public downloads), but integrity is central.

### Applicable ASVS Categories
| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V5 Input Validation | yes | Validate manifest fields (tier/split/license enums, URL well-formedness) before use; validate archive-extracted paths |
| V6 Cryptography | yes | SHA-256 integrity via stdlib `SHA` — **never hand-roll**; treat mismatch as abort (D-05) |
| V12 File/Resource | yes | Atomic `.tmp→mv`; extract zips only into a scoped temp dir; reject entries escaping it (zip-slip) |
| V2/V3/V4 (auth/session/access) | no | All sources are anonymous public HTTP; introduce no credentials |

### Known Threat Patterns
| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Tampered/corrupted download | Tampering | SHA-256 verify, hard-fail on mismatch (D-05) |
| Zip-slip on CBS `.zip` extraction | Tampering | Validate every extracted path stays within the target dir before writing |
| Supply-chain (new dep) | Tampering | No new external deps (stdlib + existing); any supporting dep gated by human-verify |
| Executing downloaded content | Elevation | Never `eval`/run downloaded bytes — they are image data only (mirrors tapqir_bridge "never eval untrusted content") |
| Hang on stalled download | Denial of Service | Consider a download timeout (mirror `tapqir_bridge:_wait_timeout`) so CI never blocks |

## Sources

### Primary (HIGH confidence)
- Repo code (VERIFIED via grep/Read): `src/LoadImages.jl` (`MultiChannelImage` ctor L126, `load_tiff` L214), `src/colocalization.jl` (`patch`/`correlation` semantics via 02-01-PLAN interfaces), `src/amortized/datagen.jl` (`cache_hash`/`subhashes`/`write_shard`/`write_meta`/`open_or_invalidate` atomic+hash+invalidate), `spike/comparator/table.jl` (`comparator_hash`/`write_table` SHA-256+atomic content-addressed writer), `spike/comparator/tapqir_bridge.jl` (skip-with-flag `(status=:skipped|:passed|:failed)` + subprocess timeout), `spike/comparator/config.jl` (Philox `MASTER_SEED` + pre-declared consts). `Project.toml`/`spike/Project.toml` (CSV/DataFrames/JLD2/Images present; Random123 in spike).
- colocalization-benchmark.com `/downloads` + `/about` — CBS format (lossless TIFF), degrees (0–90% by 10%), 3 channel pairs, sizes, CC BY-NC-SA 4.0 license, citation (Zinchuk & Grossenbacher-Zinchuk 2014 Curr Protoc Cell Biol; Sci Rep 2013 3:1365). HIGH.

### Secondary (MEDIUM confidence)
- Wikipedia "Colocalization Benchmark Source" — corroborates format/version (v2.0, 2012, lossless TIFF), references Pearson Rr / overlap coefficient. MEDIUM.
- TetraSpeck beads as multicolor coloc/registration fiducial (T7280) — multiple method papers (Nat Commun tessellation coloc; multi-color tracking). MEDIUM (concept), LOW (specific archived dataset).
- "Light My Cells" (Nature Sci Data 2026; BioImage Archive) — nucleus/mitochondria/tubulin/actin OME-TIFF channels — candidate negative-anchor SOURCE. MEDIUM (existence), LOW (fit as a co-registered two-channel segregation anchor).

### Tertiary (LOW confidence — needs validation before pinning)
- Specific IDR/BioImage Archive/Zenodo/BBBC accessions for the physical anchors — NOT verified in this session; search paths given. Gate with `checkpoint:human-verify`.

## Metadata

**Confidence breakdown:**
- In-repo mechanics (hash, atomic write, skip-with-flag, manifest, conversion): **HIGH** — read directly from committed code.
- CBS (simulated-secondary tier): **HIGH** — verified against the official source.
- Physical-anchor dataset selection: **LOW–MEDIUM** — strong candidate leads + acceptance predicate + search paths, but no single verified accession; deliberately deferred to human-verify checkpoints.
- Standard stack / no-new-deps claim: **HIGH** — confirmed against both Project.toml files and stdlib.

**Research date:** 2026-07-20
**Valid until:** ~2026-08-20 (stable; CBS/archive URLs and stdlib APIs are slow-moving. Re-confirm exact CBS zip filenames and any pinned anchor accession at fetch-script build time.)
