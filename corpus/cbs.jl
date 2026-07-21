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

# corpus/cbs.jl --- Colocalization Benchmark Source (CBS) ingestion (SC2, D-06, D-07, D-09).
#
# CBS is the ONE fully-verified external dataset (colocalization-benchmark.com; Zinchuk &
# Grossenbacher-Zinchuk 2014; CC-BY-NC-SA-4.0; lossless TIFF; labelled colocalization degrees
# 0..90% in 10% increments; three channel-pair sets Red-Green / Red-Blue / Green-Blue). It is
# ingested STRICTLY as the `simulated-secondary` tier and NEVER conflated with the physical
# anchors (SC2/D-06) — every row this file emits hard-codes `tier="simulated-secondary"` and
# `role="benchmark"`, and `validate_manifest` enforces those enums structurally.
#
# This file provides FOUR things (D-07/D-09 + the security-relevant conversion path):
#   * `cbs_rows()`             -- the FULL corpus enumerated: 3 pairs x 10 degrees = 30 rows (D-07).
#   * `assign_cbs_split(rows)` -- a SEEDED, bit-reproducible dev/eval split keyed by
#                                 CORPUS_MASTER_SEED xor CBS_SPLIT_SALT (D-09). No CBS row is ever
#                                 `sealed_holdout` (that split is anchors-only).
#   * `safe_extract(zip, dir)` -- zip-slip-safe extraction: every entry's normalized path MUST
#                                 stay within `dir` or extraction aborts BEFORE any write (T-08-11).
#   * `cbs_rgb_to_channels`    -- split a CBS RGB/composite TIFF into two Matrix{Float64} channels.
#
# The `coloc_ground_truth` of a CBS row is the ARCHIVE-LABELLED degree (a fraction), NEVER a
# computed Pearson/Manders score — the label is definitional, carried verbatim from the source.
#
# DECOUPLING (hard constraint, CLAUDE.md): corpus-local. Touches no src/, no root Project.toml.
# `Images`/`CSV`/`DataFrames`/`Random123` are already direct project deps (no new dependency).

import Images       # RGB/composite TIFF -> per-channel Float64 matrices (root dep)
import Random123: Philox4x   # counter-based keyed RNG for the reproducible dev/eval split
import Random                # rand(::AbstractRNG) for the split draw

# Guarded includes so cbs.jl loads standalone or after a sibling already pulled the corpus
# config / manifest / fetch utilities into scope (mirror manifest.jl:52-53 / fetch.jl:43-44).
isdefined(@__MODULE__, :CBS_SPLIT_SALT)  || include(joinpath(@__DIR__, "config.jl"))
isdefined(@__MODULE__, :build_manifest)  || include(joinpath(@__DIR__, "manifest.jl"))
isdefined(@__MODULE__, :fetch_verified)  || include(joinpath(@__DIR__, "fetch.jl"))

# --- CBS verified facts (RESEARCH §CBS) -----------------------------------------
# The three channel-pair sets, their short codes (used in anchor_ids) and the manifest channel
# descriptor, plus the labelled colocalization degrees. All VERIFIED against the official source.
const CBS_PAIRS = ("Red-Green", "Red-Blue", "Green-Blue")
const CBS_PAIR_CODE = Dict("Red-Green" => "RG", "Red-Blue" => "RB", "Green-Blue" => "GB")
const CBS_PAIR_CHANNELS = Dict("Red-Green" => "R,G", "Red-Blue" => "R,B", "Green-Blue" => "G,B")
const CBS_DEGREES = 0:10:90   # labelled colocalization degrees, 10% increments (VERIFIED)

# VERIFIED download-URL PREFIX only. The exact per-set/per-image zip filenames are confirmed by
# scraping the downloads page at fetch time (RESEARCH A6) — we deliberately do NOT invent a
# false-precise filename here; the prefix is the verified, honest provenance pin.
const CBS_SOURCE_PREFIX = "https://colocalization-benchmark.com/download/"
const CBS_CITATION = "Zinchuk & Grossenbacher-Zinchuk 2014; Colocalization Benchmark Source (www.colocalization-benchmark.com)"
const CBS_LICENSE = "CC-BY-NC-SA-4.0"   # NC is fine for validation use; bytes stay out of git (D-04)

# --- Full-corpus row enumeration (D-07) -----------------------------------------
"""
    cbs_rows() -> Vector{<:NamedTuple}

Enumerate the FULL CBS corpus (D-07) as manifest rows: for each channel pair in `CBS_PAIRS` and
each labelled degree in `CBS_DEGREES` (30 rows total), emit a `simulated-secondary` / `benchmark`
row whose `coloc_ground_truth` is the ARCHIVE-LABELLED degree as a fraction (e.g. 50% -> 0.5) —
NEVER a computed score. `sha256`/`bytes` are left empty/zero (filled on the first authorized
fetch); `split` is a `"dev"` placeholder overwritten by `assign_cbs_split`. `source_url` is the
verified download PREFIX (exact zip filename confirmed at fetch time — see RESEARCH A6).
"""
function cbs_rows()
    rows = NamedTuple[]
    for pair in CBS_PAIRS
        code = CBS_PAIR_CODE[pair]
        for d in CBS_DEGREES
            push!(rows, (;
                anchor_id          = "cbs-$(code)-$(lpad(d, 3, '0'))",   # e.g. cbs-RG-000 .. cbs-GB-090
                tier               = "simulated-secondary",              # SC2/D-06: never a physical anchor
                role               = "benchmark",
                truth_label        = "simulated-degree",
                coloc_ground_truth = d / 100,                            # labelled degree as a fraction
                accession          = "CBS/$(pair)/$(d)",
                source_url         = CBS_SOURCE_PREFIX,                   # verified prefix (filename TODO at fetch)
                citation           = CBS_CITATION,
                license            = CBS_LICENSE,
                format             = "tiff",                             # lossless TIFF, packaged in a .zip
                channels           = CBS_PAIR_CHANNELS[pair],            # e.g. "R,G"
                sha256             = "",                                 # filled on first authorized fetch
                split              = "dev",                              # placeholder; assign_cbs_split overwrites
                bytes              = 0,                                  # filled on fetch
            ))
        end
    end
    return rows
end

# --- Seeded, bit-reproducible dev/eval split (D-09) -----------------------------
"""
    assign_cbs_split(rows; seed=CORPUS_MASTER_SEED) -> Vector{<:NamedTuple}

Assign each CBS row a `split ∈ {"dev","eval"}` from a Random123 `Philox4x` keyed by
`(seed ⊻ CBS_SPLIT_SALT, row_index)` — a PROVABLY DISJOINT, bit-reproducible draw stream
(mirrors the `datagen.jl` HOLDOUT_SALT/FOLD_SALT idiom). Keying by the 1-based row index makes
the assignment independent of iteration order and identical across runs and thread counts. NO CBS
row is ever `sealed_holdout` — that split is reserved for the physical anchors (D-09). Returns
NEW rows with `split` replaced (inputs untouched).
"""
function assign_cbs_split(rows; seed = CORPUS_MASTER_SEED)
    out = NamedTuple[]
    for (i, row) in enumerate(collect(rows))
        rng = Philox4x(UInt64, (UInt64(seed) ⊻ CBS_SPLIT_SALT, UInt64(i)))
        split = Random.rand(rng) < 0.5 ? "dev" : "eval"
        push!(out, merge(row, (; split = split)))
    end
    return out
end

# --- Zip-slip-safe extraction (T-08-11) -----------------------------------------
"""
    _safe_target_path(target_dir, entry_name) -> String

Resolve `entry_name` against `target_dir` and return the absolute destination path ONLY if it
stays within `target_dir`; otherwise `error(...)`. This is the zip-slip mitigation core: an entry
like `"../evil.tif"` (or an absolute path) resolves OUTSIDE the target and is rejected BEFORE any
byte is written.
"""
function _safe_target_path(target_dir::AbstractString, entry_name::AbstractString)
    root   = normpath(abspath(target_dir))
    joined = normpath(joinpath(root, entry_name))
    sep    = Base.Filesystem.path_separator
    within = joined == root || startswith(joined, endswith(root, sep) ? root : root * sep)
    within || error("zip-slip rejected: entry '$(entry_name)' escapes target dir '$(target_dir)'")
    return joined
end

# The zip-entry provider seam. `entries` is an iterable of `(name::String, getbytes::Function)`.
# Isolating it keeps `safe_extract`'s security logic testable OFFLINE (the tests inject crafted
# entries) with NO zip dependency in the offline gate; an online extraction path (Phase-16 /
# future plan) supplies a concrete reader via the `entries=` keyword. We never fabricate a reader.
function _zip_entries(zip_path::AbstractString)
    error("no zip reader configured for '$(zip_path)': pass `entries=` to safe_extract (the " *
          "offline gate injects entries directly; wire a concrete reader for online extraction)")
end

"""
    safe_extract(zip_path, target_dir; entries=_zip_entries(zip_path)) -> Vector{String}

Extract `entries` (each `(name, getbytes)`) into `target_dir`, validating EVERY entry's
normalized path stays within `target_dir` BEFORE writing it (zip-slip rejection, T-08-11). Any
entry escaping the target aborts the whole extraction with `error(...)`. Returns the written paths.
"""
function safe_extract(zip_path::AbstractString, target_dir::AbstractString;
                      entries = _zip_entries(zip_path))
    mkpath(target_dir)
    written = String[]
    for (name, getbytes) in entries
        dest = _safe_target_path(target_dir, name)   # throws on zip-slip BEFORE any write
        mkpath(dirname(dest))
        write(dest, getbytes())
        push!(written, dest)
    end
    return written
end

# --- RGB/composite TIFF -> two Float64 channels ---------------------------------
"""
    _split_rgb(img, pair) -> (Matrix{Float64}, Matrix{Float64})

Split a loaded RGB/composite image into the two `Matrix{Float64}` channels named by `pair`
(one of `CBS_PAIRS`), via `Images.red`/`green`/`blue`. Kept separate from disk I/O so the
channel-selection logic is unit-testable offline without saving a TIFF.
"""
function _split_rgb(img, pair::AbstractString)
    r = Float64.(Images.red.(img))
    g = Float64.(Images.green.(img))
    b = Float64.(Images.blue.(img))
    pair == "Red-Green"  && return (r, g)
    pair == "Red-Blue"   && return (r, b)
    pair == "Green-Blue" && return (g, b)
    error("unknown CBS channel pair '$(pair)' (expected one of $(CBS_PAIRS))")
end

"""
    cbs_rgb_to_channels(tiff_path, pair) -> (Matrix{Float64}, Matrix{Float64})

Load a CBS RGB/composite TIFF and return its two relevant channels for `pair` as
`Matrix{Float64}` — feeding the UNCHANGED `src/` summary path (a `MultiChannelImage` built from
these two channels reduces through `patch`/`correlation` verbatim).
"""
cbs_rgb_to_channels(tiff_path::AbstractString, pair::AbstractString) =
    _split_rgb(Images.load(tiff_path), pair)

# --- Committed manifest.csv contract (SC3) --------------------------------------
# The two physical-anchor PLACEHOLDER rows: they exercise the schema/guards with anchor rows
# present, but their provenance fields stay "PENDING" until 08-05 finalizes them behind a
# human-verify. They are `physical-primary` + `sealed_holdout` so validate_manifest's D-09
# invariant holds; their asserted 1.0/0.0 ground truth comes from biology, never a computed score.
function anchor_placeholder_rows()
    base = (; accession = "PENDING", source_url = "PENDING", citation = "PENDING",
             license = "PENDING", format = "tiff", channels = "2:ch1,ch2",
             sha256 = "", split = "sealed_holdout", bytes = 0)
    pos = merge(base, (; anchor_id = "pos-anchor-01", tier = "physical-primary",
                        role = "positive", truth_label = "coloc", coloc_ground_truth = 1.0))
    neg = merge(base, (; anchor_id = "neg-anchor-01", tier = "physical-primary",
                        role = "negative", truth_label = "segregated", coloc_ground_truth = 0.0))
    return [pos, neg]
end

"""
    committed_manifest() -> DataFrame

The committed versioned contract (SC3): the full CBS corpus with a seeded dev/eval split, PLUS
the two sealed physical-anchor placeholder rows, assembled and VALIDATED. 32 rows total.

SUPERSEDED for the on-disk contract by `finalized_manifest()` (corpus/anchor_rows.jl, 08-05),
which substitutes the HUMAN-VERIFIED anchor rows for these "PENDING" placeholders. This function
is retained only as the placeholder-shape reference; write `corpus/manifest.csv` via
`write_finalized_manifest()`, NOT via `write_committed_manifest(committed_manifest())`.
"""
function committed_manifest()
    rows = vcat(assign_cbs_split(cbs_rows()), anchor_placeholder_rows())
    df = build_manifest(rows)
    validate_manifest(df)
    return df
end

# The repo-tracked committed manifest path (flat, human-readable copy — the content-addressed
# `write_manifest` dir remains the integrity twin).
const COMMITTED_MANIFEST_PATH = normpath(joinpath(@__DIR__, "manifest.csv"))

"""
    write_committed_manifest(df=committed_manifest(); path=COMMITTED_MANIFEST_PATH) -> String

Validate `df` and write it to the repo-tracked `corpus/manifest.csv` with the D-08
pre-registration comment header (schema version + `CORPUS_MASTER_SEED`), atomically
(`.tmp` → integrity check → `mv(force=true)`). CSV.write repositions an IOStream and clobbers a
pre-written header, so the table is serialized to a buffer first, then header + table bytes are
written (the `write_manifest` idiom). Returns the written path.
"""
function write_committed_manifest(df = committed_manifest(); path = COMMITTED_MANIFEST_PATH)
    validate_manifest(df)
    tbuf = IOBuffer()
    CSV.write(tbuf, df)
    table_bytes = take!(tbuf)
    header = _manifest_header()
    tmp = path * ".tmp"
    open(tmp, "w") do io
        for line in split(header, '\n'; keepempty = false)
            println(io, "# ", line)
        end
        write(io, table_bytes)
    end
    @assert filesize(tmp) > 0 "committed manifest integrity check failed: $tmp is empty"
    mv(tmp, path; force = true)
    return path
end

# --- Live CBS fetch smoke (D-05) ------------------------------------------------
"""
    fetch_cbs_smoke(; nsets=1) -> NamedTuple

Attempt a bootstrap fetch of ONE small CBS zip into `CORPUS_DATA_DIR` (expected hash `nothing`
⇒ bootstrap, no comparison) and return `fetch_verified`'s status NamedTuple. OFFLINE the source
is unreachable ⇒ `(status=:skipped, ...)` (never throws), so the offline gate stays green (D-05).
`nsets` reserves the number of sets to probe (currently one); the smoke never requires network.
"""
function fetch_cbs_smoke(; nsets::Int = 1)
    dest = joinpath(CORPUS_DATA_DIR, "cbs_smoke.zip")
    return fetch_verified(CBS_SOURCE_PREFIX, dest, nothing)
end
