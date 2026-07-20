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

# corpus/manifest.jl --- the versioned validation data contract (SC3) + structural guards
# (SC2/D-06, D-09) + the anchor acceptance predicate (D-03) + a content-addressed writer (D-08).
#
# This file is where the corpus tier/split policy is made TRUE BY CONSTRUCTION rather than by
# convention, mirroring the Phase-3 "there is no `standardize_all` symbol, so cross-split
# leakage is impossible" discipline (src/amortized/datagen.jl):
#
#   * D-06 / SC2 (never conflate tiers): there is NO default accessor returning both tiers.
#     `physical_anchors` and `cbs_benchmark` are SEPARATE typed accessors; the ONLY function
#     that unions the physical + simulated tiers is the explicitly-named `mixed_tier`, so a
#     caller can never silently blend physical ground truth with the simulated CBS benchmark.
#   * D-09 (sealed anchor holdout): every physical-primary row is `sealed_holdout`
#     (enforced by `validate_manifest`), the default iterator `corpus_default` structurally
#     filters those rows out, and `open_sealed_holdout(; reason)` is the SOLE reachable path
#     to the anchors -- called ONLY by Phase-16 evaluation code, never during model dev.
#   * D-03 (non-circularity): `validate_anchor_row` rejects any physical anchor whose truth is
#     not a biological {coloc,segregated} label with an open license + accession, and whose
#     `coloc_ground_truth` is not an asserted 1.0/0.0 (a Pearson/Manders-derived score is
#     rejected -- Pitfall 1).
#   * D-08 (detectable contract edits): `write_manifest` content-addresses the committed
#     manifest (source bytes + manifest content digest) so any edit resolves to a new dir.
#
# DECOUPLING (hard constraint, CLAUDE.md): corpus-local. Touches no src/, no root Project.toml.
# DataFrames/CSV/JLD2 are already direct project deps (promoted in Phase 9).

using DataFrames   # manifest table container + typed row accessors
using CSV          # human-readable manifest.csv writer (D-08)
using JLD2         # exact-reload manifest.jld2 twin (D-08)

# Guarded includes so manifest.jl loads standalone or after a sibling pulled the corpus
# config + hash utilities into scope (mirror table.jl:47-50 / hash.jl:37).
isdefined(@__MODULE__, :MANIFEST_COLUMNS) || include(joinpath(@__DIR__, "config.jl"))
isdefined(@__MODULE__, :manifest_hash)    || include(joinpath(@__DIR__, "hash.jl"))

# --- Row assembly ---------------------------------------------------------------
"""
    build_manifest(rows) -> DataFrame

Assemble `rows` (an iterable of per-image rows, each a `NamedTuple` / `DataFrameRow` /
`Dict` carrying the 14 `MANIFEST_COLUMNS` fields) into a `DataFrame` with EXACTLY those
columns, in the declared order — the tidy `DataFrame(col=..., ...)` idiom from
`spike/comparator/table.jl:build_table`. Column order is fixed so the committed CSV is a
stable, reviewable contract artifact. This function does NOT validate; call
`validate_manifest` on the result to enforce the schema.
"""
function build_manifest(rows)
    rowvec = collect(rows)
    pairs = Pair{Symbol,Vector}[]
    for c in MANIFEST_COLUMNS
        push!(pairs, c => [_rowget(r, c) for r in rowvec])
    end
    return DataFrame(pairs)
end

# Column extraction across NamedTuple / DataFrameRow (both support `hasproperty`) and Dict.
_rowget(row, k::Symbol) = hasproperty(row, k) ? getproperty(row, k) : row[k]

# --- Schema validation (SC3, D-09) ----------------------------------------------
"""
    validate_manifest(df) -> df

Throw `error(...)` on any schema violation and otherwise return `df` unchanged:
  * a missing required `MANIFEST_COLUMNS` column,
  * a `tier` not in `TIERS`, a `split` not in `SPLITS`, or a `role` not in `ROLES`,
  * the D-09 invariant: a `physical-primary` row whose `split` != `sealed_holdout`
    (physical anchors are a SEALED blind holdout — a dev/eval physical row is a leak).
The enum checks are STRUCTURAL (validated against the pre-declared `config.jl` vocabulary),
not conventional.
"""
function validate_manifest(df)
    cols = propertynames(df)
    for c in MANIFEST_COLUMNS
        c in cols || error("manifest schema violation: missing required column :$c")
    end
    for (i, t) in enumerate(df.tier)
        t in TIERS || error("manifest schema violation: row $i tier '$t' not in $(TIERS)")
    end
    for (i, s) in enumerate(df.split)
        s in SPLITS || error("manifest schema violation: row $i split '$s' not in $(SPLITS)")
    end
    for (i, r) in enumerate(df.role)
        r in ROLES || error("manifest schema violation: row $i role '$r' not in $(ROLES)")
    end
    for i in 1:nrow(df)
        if df.tier[i] == "physical-primary" && df.split[i] != "sealed_holdout"
            error("D-09 violation: physical-primary row $i has split '$(df.split[i])' " *
                  "(every physical anchor MUST be sealed_holdout)")
        end
    end
    return df
end

# --- Typed tier accessors (D-06 / SC2 structural separation) ---------------------
# NOTE: there is deliberately NO default accessor that returns BOTH tiers. The only both-tier
# path is the explicitly-named `mixed_tier` below — mirroring datagen.jl's "no standardize_all
# symbol exists" note, so physical ground truth can never silently blend with the CBS benchmark.

"""
    physical_anchors(df) -> DataFrame

Rows that are `physical-primary` AND NOT `sealed_holdout`. Because `validate_manifest` forces
every physical anchor to be sealed, on a valid manifest this returns ZERO rows BY CONSTRUCTION —
the physical anchors are unreachable through a default accessor. Access to the anchors is ONLY
via `open_sealed_holdout(; reason)` (D-09 anti-snooping).
"""
physical_anchors(df) = df[(df.tier .== "physical-primary") .& (df.split .!= "sealed_holdout"), :]

"""
    cbs_benchmark(df) -> DataFrame

Rows in the `simulated-secondary` tier (the CBS benchmark). NEVER returns a physical anchor.
"""
cbs_benchmark(df) = df[df.tier .== "simulated-secondary", :]

"""
    corpus_default(df) -> DataFrame

The DEFAULT corpus iterator: every non-`sealed_holdout` row. Structurally filters out the
sealed physical anchors so no default model-development path can read them (D-09). Asserts the
result carries no `sealed_holdout` row.
"""
function corpus_default(df)
    out = df[df.split .!= "sealed_holdout", :]
    @assert !any(==("sealed_holdout"), out.split) "corpus_default leaked a sealed_holdout row"
    return out
end

"""
    open_sealed_holdout(df; reason::AbstractString) -> DataFrame

The SOLE reachable path to the sealed physical anchors (D-09). Requires a non-empty `reason`
(an audit trail for why the blind holdout is being opened) and returns the `physical-primary`
`sealed_holdout` rows. This is called ONLY by Phase-16 blind-evaluation code — never during any
model development or tuning in Phases 8-15.
"""
function open_sealed_holdout(df; reason::AbstractString)
    @assert !isempty(reason) "open_sealed_holdout requires a non-empty reason (D-09 audit)"
    return df[(df.tier .== "physical-primary") .& (df.split .== "sealed_holdout"), :]
end

"""
    mixed_tier(df) -> DataFrame

The ONLY function that returns both tiers together — the physical anchors UNIONED with the
CBS benchmark. Callers must name it explicitly; there is no default "everything" accessor, so
tier blending is always a deliberate, greppable act (SC2 / D-06).
"""
function mixed_tier(df)
    phys = df[df.tier .== "physical-primary", :]
    sim  = df[df.tier .== "simulated-secondary", :]
    return vcat(phys, sim)
end

# --- Anchor acceptance predicate (D-03) -----------------------------------------
"""
    validate_anchor_row(row) -> Bool

Return `true` only if `row` qualifies as a physical anchor under the D-03 acceptance
predicate (non-circularity):
  (1) `truth_label` is a BIOLOGICAL label in {"coloc","segregated"} — a computed score
      (e.g. "0.73") is rejected;
  (2) `license` is a non-empty (redistributable) tag;
  (4) `accession` is a non-empty provenance pin (DOI / archive id);
and the anchor's `coloc_ground_truth` is an ASSERTED 1.0/0.0 rather than a measured score —
a Pearson/Manders-derived value is rejected (Pitfall 1). Any failure returns `false`.
"""
function validate_anchor_row(row)
    truth = _rowget(row, :truth_label)
    lic   = _rowget(row, :license)
    acc   = _rowget(row, :accession)
    gt    = _rowget(row, :coloc_ground_truth)
    (truth in ("coloc", "segregated"))                || return false   # (1) biological label
    (!ismissing(lic) && !isempty(string(lic)))        || return false   # (2) open license
    (!ismissing(acc) && !isempty(string(acc)))        || return false   # (4) accession
    (!ismissing(gt) && (gt == 1.0 || gt == 0.0))      || return false   # asserted, not computed
    return true
end

# --- Content-addressed atomic writer (D-08) -------------------------------------
# FIXED-ORDER contract source list: the corpus contract code whose bytes co-define the manifest
# content address. Any edit to these files (or the manifest content) resolves to a NEW dir, so a
# contract edit is detectable (D-08). @__DIR__-relative (this file lives in corpus/).
const MANIFEST_SRC_FILES = String[
    joinpath(@__DIR__, "config.jl"),
    joinpath(@__DIR__, "hash.jl"),
    joinpath(@__DIR__, "manifest.jl"),
]

# SHA-256 digest of the manifest CONTENT (its CSV serialization), folded into the content-address
# config so two different manifests never collide on the same dir even for identical source code.
function _manifest_content_digest(df)
    io = IOBuffer()
    CSV.write(io, df)
    return bytes2hex(SHA.sha256(take!(io)))
end

# Pre-registration comment header (D-08): quotes the schema version + corpus master seed VERBATIM
# so the contract that governed a committed manifest is recorded alongside it.
function _manifest_header()
    io = IOBuffer()
    println(io, "ProteinCoLoc v2.0 external ground-truth corpus manifest — pre-registration (D-08)")
    println(io, "MANIFEST_SCHEMA_VERSION = $(MANIFEST_SCHEMA_VERSION)   # bump on schema shape/meaning change")
    println(io, "CORPUS_MASTER_SEED      = $(repr(CORPUS_MASTER_SEED))  # Random123 (Philox) corpus master seed")
    return String(take!(io))
end

"""
    write_manifest(df; outdir, src_files=MANIFEST_SRC_FILES) -> String

Validate `df`, then persist it atomically to a CONTENT-ADDRESSED directory
`joinpath(outdir, manifest_hash(src_files, canonical_config(config)))`, writing both:
  * `manifest.csv` — human-readable, prefixed with the D-08 pre-registration comment header
    (schema version + `CORPUS_MASTER_SEED`), and
  * `manifest.jld2` — an exact-reload twin (the DataFrame + header meta + content config).
Each file is written via the `.tmp` → integrity-check → `mv(...; force=true)` idiom
(`spike/comparator/table.jl:write_table`), so a crash never leaves a half-written artifact.
Two writes of the same `(df, sources)` resolve to the IDENTICAL dir. Returns the dir path.
"""
function write_manifest(df; outdir, src_files = MANIFEST_SRC_FILES)
    validate_manifest(df)
    config = (schema_version = MANIFEST_SCHEMA_VERSION,
              master_seed    = CORPUS_MASTER_SEED,
              content        = _manifest_content_digest(df))
    dir = joinpath(outdir, manifest_hash(src_files, canonical_config(config)))
    isdir(dir) || mkpath(dir)
    header = _manifest_header()

    # --- manifest.csv (humans): pre-registration header as CSV comment lines, then the table.
    # CSV.write repositions an IOStream (clobbering a pre-written header), so serialize the table
    # to a dedicated buffer first, then write the comment header + table bytes to the file.
    tbuf = IOBuffer()
    CSV.write(tbuf, df)
    table_bytes = take!(tbuf)
    csv_path = joinpath(dir, "manifest.csv")
    csv_tmp  = csv_path * ".tmp"
    open(csv_tmp, "w") do io
        for line in split(header, '\n'; keepempty = false)
            println(io, "# ", line)
        end
        write(io, table_bytes)
    end
    @assert filesize(csv_tmp) > 0 "manifest csv integrity check failed: $csv_tmp is empty"
    mv(csv_tmp, csv_path; force = true)

    # --- manifest.jld2 (exact reload): the DataFrame + the header meta + the content config.
    jld_path = joinpath(dir, "manifest.jld2")
    jld_tmp  = jld_path * ".tmp"
    JLD2.jldsave(jld_tmp; manifest = df, meta = header, config = config)
    JLD2.jldopen(jld_tmp, "r") do f
        @assert haskey(f, "manifest") "manifest jld2 integrity check failed: $jld_tmp missing manifest"
    end
    mv(jld_tmp, jld_path; force = true)

    return dir
end
