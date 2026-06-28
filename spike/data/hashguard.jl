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

# spike/data/hashguard.jl --- DATA-02 content-hash version guard (D-05).
#
# The cache-identity boundary: a SHA-256 over the data-DEFINING source bytes PLUS
# a canonical serialization of the generating config. The digest names the cache
# directory and is stored in meta.jld2, so ANY change to the simulator/prior/ghat/
# encode source OR to a hashed config field (N, shard_size, imsize-spec,
# summary-spec, θ-dim, master_seed, k) yields a different digest and auto-
# invalidates stale data (D-05). Per-file sub-hashes make a mismatch DIAGNOSABLE
# -- the loader can report exactly WHICH source file changed (T-03-06 mitigation).
#
# WHY hash the two frozen src/ files too: they are baseline-frozen (f581d95) but
# DO define the summary (patch/correlation), so hashing them is a defensive
# decoupling-breach tripwire -- an accidental src/ edit auto-invalidates the cache.
#
# DECOUPLING (hard constraint, CLAUDE.md): reads source bytes only; touches no
# src/, edits nothing. `using SHA` is the Julia stdlib -- no Pkg.add.

using SHA   # stdlib: sha256 / SHA256_CTX / update! / digest! (no add)

# --- FIXED-order data-defining source files (D-05) ------------------------------
# Order is part of the hash CONTRACT -- never reorder. @__DIR__-relative
# (hashguard.jl lives in spike/data/). The last two are the frozen src/ tripwire.
const HASH_SRC_FILES = String[
    joinpath(@__DIR__, "..", "simulator", "forward.jl"),
    joinpath(@__DIR__, "..", "simulator", "prior.jl"),
    joinpath(@__DIR__, "..", "simulator", "ghat.jl"),
    joinpath(@__DIR__, "..", "contract.jl"),
    joinpath(@__DIR__, "encode.jl"),
    joinpath(@__DIR__, "..", "..", "src", "colocalization.jl"),
    joinpath(@__DIR__, "..", "..", "src", "LoadImages.jl"),
]

"""
    canonical(config::NamedTuple) -> String

Deterministic, sort-by-field-name serialization of `config` so that field ORDER
can never perturb the hash. Each value is rendered with `repr` (stable for the
Int / UInt64 / Tuple / Vector / Symbol fields the generating config carries).
Works on ANY NamedTuple (full or partial), so callers may hash a subset of
fields without committing to the full generating schema.
"""
function canonical(config::NamedTuple)::String
    ks = sort(collect(keys(config)))                       # Symbols, sorted by name
    return join(("$(k)=$(repr(getproperty(config, k)))" for k in ks), ";")
end

"""
    cache_hash(config::NamedTuple; src_files=HASH_SRC_FILES) -> String

The D-05 content hash: fold `SHA.update!` over the raw bytes of each `src_files`
entry (FIXED order, NUL-separated to remove file-boundary ambiguity) then over
`canonical(config)`, returning the hex digest. Identical (config, source bytes)
=> identical digest; perturbing ANY source byte OR ANY hashed config field flips
it. This digest names the cache directory and is stored in meta.jld2.
"""
function cache_hash(config::NamedTuple; src_files = HASH_SRC_FILES)
    ctx = SHA.SHA256_CTX()
    for f in src_files                                     # FIXED order = deterministic
        SHA.update!(ctx, read(f))
        SHA.update!(ctx, UInt8[0x00])                      # inter-file separator
    end
    SHA.update!(ctx, Vector{UInt8}(canonical(config)))
    return bytes2hex(SHA.digest!(ctx))
end

"""
    subhashes(src_files=HASH_SRC_FILES) -> Dict{String,String}

Per-file SHA-256 hex digests keyed by source path, so a `cache_hash` mismatch is
DIAGNOSABLE -- the caller can diff stored vs recomputed sub-hashes to name exactly
WHICH source file changed (D-05 / T-03-06). One entry per `src_files` member.
"""
function subhashes(src_files = HASH_SRC_FILES)
    return Dict{String,String}(f => bytes2hex(SHA.sha256(read(f))) for f in src_files)
end
