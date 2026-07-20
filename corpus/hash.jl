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

# corpus/hash.jl --- SHA-256 content-hash utilities for the external corpus (D-04, D-08).
#
# The integrity spine of the corpus module. Every download-verification and manifest
# content-address goes through these primitives. They use CRYPTOGRAPHIC SHA-256
# (`SHA.sha256` / `SHA256_CTX`), NEVER Base `hash`: a content-hash MISMATCH must be a
# tamper/corruption signal a caller can HARD-abort on (D-05), which a non-cryptographic
# digest cannot support. The SHA-256 fold + canonical-config serializer mirror the
# VERIFIED spike/comparator/table.jl `comparator_hash`/`_cmp_canonical` idiom.
#
# DECOUPLING (hard constraint, CLAUDE.md): corpus-local. Touches no src/, no root
# Project.toml. `SHA` is a Julia stdlib (ships with Julia; no dependency added).

using SHA   # stdlib: SHA256_CTX / sha256 / update! / digest! for the content hash (D-04)

# Guarded include so hash.jl loads standalone or after a sibling already pulled the
# corpus config into scope (mirror table.jl:47-50).
isdefined(@__MODULE__, :MANIFEST_SCHEMA_VERSION) || include(joinpath(@__DIR__, "config.jl"))

# Streaming block size for file_sha256. A >100 MB CBS archive must never be fully read
# into memory at once, so file_sha256 folds fixed-size blocks through the SHA context.
const _FILE_SHA256_CHUNK = 1024 * 1024   # 1 MiB read blocks

"""
    file_sha256(path::AbstractString) -> String

Return the 64-char lowercase hex SHA-256 digest of the file at `path`, streaming the
file in `_FILE_SHA256_CHUNK`-byte blocks so a large archive is never loaded fully into
memory. The result is byte-for-byte equal to `bytes2hex(SHA.sha256(read(path)))` (SHA-256
is a streaming hash, so chunking does not change the digest) and is deterministic across
repeated calls. This is THE download-verification primitive (D-04/D-05).
"""
function file_sha256(path::AbstractString)
    ctx = SHA.SHA256_CTX()
    open(path, "r") do io
        while !eof(io)
            SHA.update!(ctx, read(io, _FILE_SHA256_CHUNK))
        end
    end
    return bytes2hex(SHA.digest!(ctx))
end

"""
    canonical_config(config::NamedTuple) -> String

Deterministic, sort-by-key serialization of `config` (each value via `repr`, joined
`"k=repr(v)"` by `;`) so field INSERTION ORDER can never perturb a content hash. A local
mirror of `spike/comparator/table.jl:_cmp_canonical` — reimplemented here so the corpus
module never reaches into the comparator's hashing machinery (decoupling).
"""
function canonical_config(config::NamedTuple)::String
    ks = sort(collect(keys(config)))
    return join(("$(k)=$(repr(getproperty(config, k)))" for k in ks), ";")
end

"""
    manifest_hash(src_files::Vector{String}, canonical_config::String) -> String

Fold `SHA.update!` over the raw bytes of each `src_files` entry (in the FIXED order given,
NUL-separated) then over the `canonical_config` string, returning the 64-char hex digest.
Mirrors `spike/comparator/table.jl:comparator_hash`.

FIXED-ORDER CONTRACT: the digest is order-sensitive BY DESIGN — reordering `src_files`
changes the hash. The caller MUST pass `src_files` in a stable, documented order (e.g. the
declared `MANIFEST_COLUMNS`/source list) so identical contract inputs always content-address
to the identical digest, and any edit to a source or the config is detectable.
"""
function manifest_hash(src_files::Vector{String}, canonical_config::AbstractString)
    ctx = SHA.SHA256_CTX()
    for f in src_files                        # FIXED order = deterministic content address
        SHA.update!(ctx, read(f))
        SHA.update!(ctx, UInt8[0x00])         # inter-file separator (no boundary ambiguity)
    end
    SHA.update!(ctx, Vector{UInt8}(String(canonical_config)))
    return bytes2hex(SHA.digest!(ctx))
end

"""
    subhashes(src_files) -> Dict{String,String}

Per-file SHA-256 hex digests keyed by source path, so a `manifest_hash` mismatch is
DIAGNOSABLE — the caller can diff stored vs recomputed sub-hashes to name exactly WHICH
source changed. Re-expressed over cryptographic SHA-256 (`file_sha256`), NOT the Base
`hash` variant from `src/amortized/datagen.jl:subhashes`: corpus download integrity
requires a cryptographic digest (D-04/D-05 caveat, 08-PATTERNS L91-96).
"""
function subhashes(src_files)
    return Dict{String,String}(f => file_sha256(f) for f in src_files)
end
