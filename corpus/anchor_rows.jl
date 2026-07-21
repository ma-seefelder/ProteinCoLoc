#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder
E-Mail: manuel.seefelder@uni-ulm.de
Postal address: Department of Gene Therapy, University of Ulm, Helmholzstr. 8/1, 89081 Ulm, Germany

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=#

# corpus/anchor_rows.jl --- the two HUMAN-VERIFIED physical ground-truth anchors (SC1, D-01,
# D-02, D-03, D-09).
#
# These are the ONE non-autonomous, non-circular deliverable of Phase 8. Both accessions below
# were confirmed by a human against the D-03 acceptance predicate before being pinned here:
#   (1) the truth label comes from BIOLOGY/PHYSICS + a citation, never from a computed
#       Pearson/Manders colocalization score (Pitfall 1);
#   (2) the licence is open/redistributable;
#   (3) the data is retrievable and convertible to two channels;
#   (4) the provenance is pinned by a resolving DOI / archive accession.
# `validate_anchor_row` (corpus/manifest.jl, D-03) re-checks (1),(2),(4) structurally at
# assembly time, and `validate_manifest` enforces the D-09 sealed-holdout invariant.
#
# ---------------------------------------------------------------------------------------------
# RECORDED DEVIATIONS FROM THE LITERAL PHASE DECISIONS (honest provenance, do NOT paper over)
# ---------------------------------------------------------------------------------------------
# D-01 DEVIATION (human-accepted substitution). D-01 literally specified a CELLULAR
#   "tandem-fluorophore construct" as the positive anchor. No open-licensed, NON-environment-
#   quenched tandem-FP dataset could be verified in any public archive: the deposited tandem-FP
#   datasets are the quenched autophagy reporters (mRFP/mCherry-EGFP-LC3), which the acceptance
#   bar DISQUALIFIES (Pitfall 2 — GFP is quenched in the acidic autolysosome, so the two channels
#   are NOT 100% colocalized in the deposited images). The human therefore ACCEPTED a documented
#   substitution: a MULTICOLOR-BEAD dataset. A single physical TetraSpeck 100 nm bead emits in
#   every colour channel, so the two channels are colocalized BY CONSTRUCTION and, unlike a
#   tandem FP, the label is state-INDEPENDENT (no pH/environment quenching failure mode at all).
#   This is a strictly STRONGER physical positive than the rejected alternative, but it is a
#   deviation from the literal D-01 wording and is recorded as such.
#
# D-02 DEVIATION (preference NOT satisfiable). D-02 PREFERRED a matched segregated negative from
#   the SAME study/archive as the positive, to control imaging/instrument confounds. No
#   qualifying single-study deposit carrying BOTH a same-particle positive AND a segregated
#   negative could be verified. The pinned negative is therefore from a DIFFERENT study and a
#   DIFFERENT archive (`matched = false` below). KNOWN LIMITATION: cross-study anchors mean
#   imaging-condition confounds (microscope, objective, exposure, detector, sample prep) are NOT
#   controlled between the positive and the negative. Phase-16 evaluation MUST report this.
#
# ---------------------------------------------------------------------------------------------
# CONTENT HASHES (D-04/D-05, T-08-16)
# ---------------------------------------------------------------------------------------------
# `sha256` is recorded ONLY from a REAL bootstrap fetch (`fetch_verified(url, dest, nothing)`).
# A digest is NEVER fabricated. Until the bootstrap fetch runs, both anchors carry the explicit
# `PENDING_FETCH` sentinel. The bootstrap was DELIBERATELY NOT run in this plan: the positive
# anchor is a ~6.3 GB archive, and a bulk download must remain an explicit human decision rather
# than a side effect of an automated plan execution. Run `bootstrap_anchor_hashes()` when online
# and authorized; the sentinel MUST be replaced by real digests before Phase 16.
#
# DECOUPLING (hard constraint, CLAUDE.md): corpus-local. Touches no src/, no root Project.toml.
# No new dependency.

# Guarded includes so this file loads standalone or after a sibling already pulled the corpus
# config / manifest / fetch / CBS utilities into scope (mirror cbs.jl:51-53).
isdefined(@__MODULE__, :MANIFEST_COLUMNS)    || include(joinpath(@__DIR__, "config.jl"))
isdefined(@__MODULE__, :validate_anchor_row) || include(joinpath(@__DIR__, "manifest.jl"))
isdefined(@__MODULE__, :fetch_verified)      || include(joinpath(@__DIR__, "fetch.jl"))
isdefined(@__MODULE__, :cbs_rows)            || include(joinpath(@__DIR__, "cbs.jl"))

# --- Pending-hash sentinel (T-08-16) --------------------------------------------
"""
    PENDING_FETCH

The explicit sentinel recorded in `sha256` when no authorized bootstrap fetch has produced a
trusted baseline digest yet. It is deliberately NOT a 64-hex string, so it can never be mistaken
for — or silently accepted as — a real SHA-256. Fabricating a digest is forbidden (D-04/D-05).
"""
const PENDING_FETCH = "PENDING-FETCH"

"""
    is_pending_hash(s) -> Bool

`true` iff `s` is the explicit `PENDING_FETCH` sentinel (no trusted baseline digest yet).
"""
is_pending_hash(s) = !ismissing(s) && string(s) == PENDING_FETCH

"""
    is_real_sha256(s) -> Bool

`true` iff `s` is a syntactically valid 64-character lowercase hexadecimal SHA-256 digest.
Used by the anchor tests to assert every physical row carries EITHER a real digest OR the
explicit pending sentinel — never anything in between (T-08-16).
"""
is_real_sha256(s) = !ismissing(s) &&
                    length(string(s)) == 64 &&
                    all(c -> c in "0123456789abcdef", lowercase(string(s)))

# --- The confirmed POSITIVE anchor (D-01, human-verified) -----------------------
# TetraSpeck 100 nm multicolor fiducial beads — RegiSTORM v1.0.0 sample data.
# PHYSICAL truth: one physical bead emits in EVERY colour channel, so the two channels are
# colocalized BY CONSTRUCTION and state-independent. This is NOT a computed coloc score and NOT
# an environment-quenched tandem (Pitfall 2 explicitly avoided).
const POS_ANCHOR_ACCESSION = "10.5281/zenodo.5509861"
const POS_ANCHOR_URL       = "https://zenodo.org/api/records/5509861/files/Sample%20Data.zip/content"
const POS_ANCHOR_LICENSE   = "CC-BY-4.0"
const POS_ANCHOR_CITATION  = "Karlsson et al. 2023, \"RegiSTORM: channel registration for " *
                             "multi-color STORM\", BMC Bioinformatics, DOI 10.1186/s12859-023-05320-1; " *
                             "sample data Zenodo DOI 10.5281/zenodo.5509861"
# Multi-channel STORM .tif frame stacks (~30,000 frames/channel) + ThunderSTORM .csv
# reconstructions. A mean/max projection over frames is REQUIRED to yield the 2-channel
# intensity image the src/ summary path consumes (conversion step owned by Phase 16).
const POS_ANCHOR_FORMAT    = "tiff-stack+csv(zip)"
const POS_ANCHOR_CHANNELS  = "2:storm-ch1,storm-ch2"

# --- The confirmed NEGATIVE anchor (D-02, human-verified) -----------------------
# "Light My Cells" — Bright Field to Fluorescence Imaging Challenge (France-BioImaging / ISBI 2024).
# BIOLOGICAL truth: nucleus (DNA stain) vs mitochondria (cytoplasmic, nucleus-excluded) are
# spatially DISJOINT BY BIOLOGY. Compartment identity comes from the acquisition design, NOT from
# a computed coloc score.
const NEG_ANCHOR_ACCESSION = "S-BIAD1047"
const NEG_ANCHOR_URL       = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/047/S-BIAD1047"
const NEG_ANCHOR_LICENSE   = "CC-BY-4.0"
const NEG_ANCHOR_CITATION  = "Light My Cells — Bright Field to Fluorescence Imaging Challenge, " *
                             "France-BioImaging / ISBI 2024; BioImage Archive S-BIAD1047"
# OME-TIFF, REMBI-compliant metadata. The challenge guarantees only "AT LEAST ONE" fluorescence
# target per field, so fields MUST be FILTERED to those carrying BOTH the nucleus and the
# mitochondria channel before use (filter step owned by Phase 16).
const NEG_ANCHOR_FORMAT    = "ome-tiff"
const NEG_ANCHOR_CHANNELS  = "2:nucleus-dna,mitochondria"

# Whether the negative anchor comes from the SAME study/archive as the positive (D-02 preference).
const ANCHORS_MATCHED = false

"""
    ANCHOR_PROVENANCE_NOTES

Machine-readable, greppable record of the two deviations from the literal D-01/D-02 decisions
(see the file header for the full rationale). Carried alongside the rows so the honest caveats
travel with the contract rather than living only in a plan summary.
"""
const ANCHOR_PROVENANCE_NOTES = (
    d01 = "DEVIATION (human-accepted): D-01 specified a cellular tandem-fluorophore construct. " *
          "No open-licensed, non-environment-quenched tandem-FP dataset could be verified; the " *
          "deposited tandem-FP datasets are the quenched mRFP/mCherry-EGFP-LC3 autophagy " *
          "reporters, disqualified by Pitfall 2. Substituted a multicolor-bead dataset " *
          "(TetraSpeck 100 nm): the SAME physical particle emits in both channels, so " *
          "colocalization is by construction and state-independent.",
    d02 = "DEVIATION (preference not satisfiable): D-02 preferred a matched same-study negative. " *
          "No qualifying single-study deposit with both a same-particle positive and a " *
          "segregated negative could be verified. The negative is from a DIFFERENT study and " *
          "archive (matched=false). KNOWN LIMITATION: imaging-condition confounds between the " *
          "positive and negative anchors are NOT controlled; Phase 16 must report this.",
    hashes = "sha256 is recorded only from a real bootstrap fetch. The bootstrap was NOT run " *
             "in plan 08-05 (the positive anchor is a ~6.3 GB archive; a bulk download stays an " *
             "explicit human decision). Both anchors carry the explicit '$(PENDING_FETCH)' " *
             "sentinel; run bootstrap_anchor_hashes() when online and authorized. A digest is " *
             "NEVER fabricated (T-08-16).",
)

# --- The anchor rows ------------------------------------------------------------
"""
    anchor_rows(; sha256_pos=PENDING_FETCH, sha256_neg=PENDING_FETCH,
                  bytes_pos=0, bytes_neg=0) -> Vector{<:NamedTuple}

The TWO human-verified physical ground-truth anchors as manifest rows (SC1): exactly one
`role="positive"` (`truth_label="coloc"`, `coloc_ground_truth=1.0`) and one `role="negative"`
(`truth_label="segregated"`, `coloc_ground_truth=0.0`). Both are `tier="physical-primary"` and
`split="sealed_holdout"` (D-09), so `validate_manifest` accepts them and no default accessor can
reach them.

The `coloc_ground_truth` values are ASSERTED FROM PHYSICS/BIOLOGY (one bead emits in both
channels; nucleus and mitochondria are disjoint compartments) — they are NEVER a computed
Pearson/Manders score (Pitfall 1, D-03).

`sha256_*`/`bytes_*` default to the `PENDING_FETCH` sentinel / 0 and are supplied by
`bootstrap_anchor_hashes()` after a REAL authorized fetch. Every row is asserted to pass
`validate_anchor_row` (the D-03 gate) before it is returned.
"""
function anchor_rows(; sha256_pos = PENDING_FETCH, sha256_neg = PENDING_FETCH,
                       bytes_pos::Integer = 0, bytes_neg::Integer = 0)
    pos = (;
        anchor_id          = "pos-tetraspeck-01",
        tier               = "physical-primary",
        role               = "positive",
        truth_label        = "coloc",
        coloc_ground_truth = 1.0,              # asserted: one physical bead in BOTH channels
        accession          = POS_ANCHOR_ACCESSION,
        source_url         = POS_ANCHOR_URL,
        citation           = POS_ANCHOR_CITATION,
        license            = POS_ANCHOR_LICENSE,
        format             = POS_ANCHOR_FORMAT,
        channels           = POS_ANCHOR_CHANNELS,
        sha256             = sha256_pos,
        split              = "sealed_holdout",  # D-09
        bytes              = bytes_pos,
    )
    neg = (;
        anchor_id          = "neg-lightmycells-01",
        tier               = "physical-primary",
        role               = "negative",
        truth_label        = "segregated",
        coloc_ground_truth = 0.0,              # asserted: nucleus vs mitochondria are disjoint
        accession          = NEG_ANCHOR_ACCESSION,
        source_url         = NEG_ANCHOR_URL,
        citation           = NEG_ANCHOR_CITATION,
        license            = NEG_ANCHOR_LICENSE,
        format             = NEG_ANCHOR_FORMAT,
        channels           = NEG_ANCHOR_CHANNELS,
        sha256             = sha256_neg,
        split              = "sealed_holdout",  # D-09
        bytes              = bytes_neg,
    )
    rows = [pos, neg]
    for r in rows
        validate_anchor_row(r) ||
            error("D-03 anchor acceptance FAILED for '$(r.anchor_id)': a physical anchor needs a " *
                  "biological {coloc,segregated} truth label, an open license, a resolving " *
                  "accession, and an ASSERTED 1.0/0.0 ground truth (never a computed score).")
    end
    return rows
end

# --- Bootstrap the trusted baseline digests (NOT run automatically) -------------
"""
    bootstrap_anchor_hashes(; datadir=CORPUS_DATA_DIR) -> NamedTuple

Run the ONE authorized bootstrap fetch per anchor — `fetch_verified(url, dest, nothing)` — which
computes and returns the content digest WITHOUT comparing (no trusted baseline exists yet), and
commits the bytes atomically into `datadir` (gitignored, D-04).

THIS IS DELIBERATELY NOT CALLED BY THE TEST GATE OR BY MANIFEST ASSEMBLY. The positive anchor is
a ~6.3 GB archive; pulling it must be an explicit, deliberate human action. Offline / unreachable
⇒ `fetch_verified` returns `(status=:skipped, ...)` and this function reports the `PENDING_FETCH`
sentinel rather than a digest. A digest is NEVER fabricated (T-08-16).

Returns `(pos=..., neg=...)` status NamedTuples. Feed the resulting digests + byte sizes into
`anchor_rows(; sha256_pos=..., sha256_neg=..., bytes_pos=..., bytes_neg=...)` and rewrite the
manifest with `write_finalized_manifest`.
"""
function bootstrap_anchor_hashes(; datadir::AbstractString = CORPUS_DATA_DIR)
    out = NamedTuple[]
    for (id, url, fname) in (("pos-tetraspeck-01", POS_ANCHOR_URL,  "pos-tetraspeck-01.zip"),
                             ("neg-lightmycells-01", NEG_ANCHOR_URL, "neg-lightmycells-01.bin"))
        dest = joinpath(datadir, fname)
        r = fetch_verified(url, dest, nothing)     # bootstrap: expected=nothing ⇒ record, no compare
        if r.status == :present
            push!(out, (; anchor_id = id, status = :present, sha256 = r.sha256,
                          bytes = filesize(dest)))
        else
            @warn "anchor bootstrap fetch SKIPPED (source unreachable) — recording $(PENDING_FETCH)" anchor_id=id reason=r.reason
            push!(out, (; anchor_id = id, status = :skipped, sha256 = PENDING_FETCH, bytes = 0))
        end
    end
    return (; pos = out[1], neg = out[2])
end

# --- Finalized committed manifest (SC1 + SC3) -----------------------------------
"""
    finalized_manifest(; kwargs...) -> DataFrame

The FINALIZED committed contract: the full CBS benchmark with its seeded dev/eval split
(`assign_cbs_split(cbs_rows())`, 30 rows, `simulated-secondary`) PLUS the two human-verified
physical anchors (`anchor_rows()`, 2 rows, `physical-primary` + `sealed_holdout`). 32 rows total.
Supersedes the 08-04 `committed_manifest()` placeholder assembly (which pinned "PENDING"
provenance). Assembled with `build_manifest` and `validate_manifest`; `kwargs` are forwarded to
`anchor_rows` (i.e. bootstrapped digests / byte sizes).
"""
function finalized_manifest(; kwargs...)
    rows = vcat(assign_cbs_split(cbs_rows()), anchor_rows(; kwargs...))
    df = build_manifest(rows)
    validate_manifest(df)
    return df
end

"""
    anchor_hash_notice(df) -> String

The human-facing notice listing every physical row still carrying the `PENDING_FETCH` sentinel,
i.e. every anchor whose trusted baseline digest MUST still be bootstrapped (online) BEFORE
Phase 16 opens the sealed holdout. Returns `""` when all anchors carry a real digest.
"""
function anchor_hash_notice(df)
    pending = String[]
    for i in 1:nrow(df)
        df.tier[i] == "physical-primary" || continue
        is_pending_hash(df.sha256[i]) && push!(pending, string(df.anchor_id[i]))
    end
    isempty(pending) && return ""
    return "NOTICE (D-04/D-05, T-08-16): $(length(pending)) physical anchor(s) carry the " *
           "'$(PENDING_FETCH)' sentinel instead of a trusted SHA-256 baseline: " *
           join(pending, ", ") * ". No digest was fabricated. Run " *
           "`bootstrap_anchor_hashes()` when online and authorized (the positive anchor is a " *
           "~6.3 GB archive — an explicit human decision), then rewrite the manifest with the " *
           "returned digests. The sentinel MUST be replaced before Phase 16 opens the sealed holdout."
end

"""
    write_finalized_manifest(df=finalized_manifest(); path=COMMITTED_MANIFEST_PATH) -> String

Validate and write the finalized manifest to the repo-tracked `corpus/manifest.csv` (D-08
pre-registration header + atomic `.tmp` → check → `mv`), reusing `write_committed_manifest`.
Prints `anchor_hash_notice(df)` if any anchor still carries the pending-hash sentinel.
"""
function write_finalized_manifest(df = finalized_manifest(); path = COMMITTED_MANIFEST_PATH)
    p = write_committed_manifest(df; path = path)
    notice = anchor_hash_notice(df)
    isempty(notice) || @warn notice
    return p
end
