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

# corpus/fetch.jl --- verified download primitive for the external corpus (D-04, D-05).
#
# `fetch_verified` composes three VERIFIED in-repo idioms into the integrity-enforcing
# fetch: the atomic `.tmp -> verify -> mv(force=true)` write (datagen.jl:write_shard /
# table.jl:write_table), the skip-with-flag NamedTuple on an unreachable source
# (tapqir_bridge.jl:tapqir_anchor), and the bounded download timeout (DOWNLOAD_TIMEOUT_S,
# mirroring tapqir_bridge.jl:_wait_timeout).
#
# THE D-05 ASYMMETRY IS STRUCTURAL (Pitfall 4): the download lives in the ONLY try/catch and
# an UNREACHABLE source degrades to (status=:skipped) so offline/CI stays green; the digest
# `!=` check lives OUTSIDE any try/catch and a content-hash MISMATCH is a HARD `error(...)`
# ABORT (tamper/corruption is not recoverable). The two paths are DELIBERATELY separate code
# regions and must never share one try/catch.
#
# DECOUPLING (hard constraint, CLAUDE.md): corpus-local. `Downloads` and `SHA` are Julia
# stdlibs (no dependency added). Touches no src/, no root Project.toml. Downloaded bytes are
# read as DATA only — never eval'd/included/run (mirrors tapqir "never eval untrusted").

using Downloads   # stdlib: HTTPS GET with TLS via libcurl; the network fetch primitive

# Guarded includes so fetch.jl loads standalone or after a sibling already pulled the
# corpus config / hash utilities into scope (mirror table.jl:47-50).
isdefined(@__MODULE__, :DOWNLOAD_TIMEOUT_S) || include(joinpath(@__DIR__, "config.jl"))
isdefined(@__MODULE__, :file_sha256)         || include(joinpath(@__DIR__, "hash.jl"))

# The ONE network seam. Held in a `Ref` so the tests can swap in a local-copy stub by
# ASSIGNING the ref (a runtime mutation, NOT a method redefinition), keeping the real
# network call isolated to exactly this default closure and sidestepping Julia world-age
# surprises. The default bounds every download by DOWNLOAD_TIMEOUT_S (T-08-04 DoS guard).
const _DOWNLOAD_HOOK = Ref{Function}(
    (url, dest) -> Downloads.download(url, dest; timeout = DOWNLOAD_TIMEOUT_S),
)

# Invoke the current hook through `invokelatest` so a test-supplied stub (defined in the
# same top-level block that calls fetch_verified) is always visible regardless of world age.
_download(url::AbstractString, dest::AbstractString) =
    Base.invokelatest(_DOWNLOAD_HOOK[], url, dest)

"""
    fetch_verified(url, dest, expected_sha256) -> NamedTuple

Download `url` to `dest`, verifying integrity, committing atomically. Enforces the D-05
asymmetry structurally:

- Reachable source whose SHA-256 == `expected_sha256` ⇒ atomic `mv(.tmp -> dest)`, returns
  `(status=:present, sha256=got)`.
- `expected_sha256 === nothing` or `""` (BOOTSTRAP — the first authorized fetch of a
  newly-pinned dataset, before a trusted baseline hash exists) ⇒ compute-and-return the hash
  WITHOUT comparing, still committing atomically; returns `(status=:present, sha256=got,
  bootstrapped=true)`. A human confirms the source is the intended one before its hash becomes
  the trusted baseline recorded in the manifest.
- Reachable source whose SHA-256 != `expected_sha256` ⇒ HARD `error(...)` ABORT (tamper/
  corruption); the `.tmp` is removed. This is NEVER a skip.
- UNREACHABLE source (the download throws) ⇒ returns `(status=:skipped, value=missing,
  reason="unreachable: ...")` and does NOT throw; leaves no `.tmp` behind (offline/CI green).

`InterruptException` is always rethrown (never swallowed by the skip path). Comparison is
case-insensitive (`lowercase`).
"""
function fetch_verified(url::AbstractString, dest::AbstractString,
                        expected_sha256::Union{AbstractString,Nothing})
    mkpath(dirname(dest))
    tmp = dest * ".tmp"

    # --- Path A (skip-with-flag): the download. The ONLY try/catch, wrapping ONLY the
    # network fetch. An unreachable source degrades gracefully — never throws (D-05).
    try
        _download(url, tmp)
    catch e
        e isa InterruptException && rethrow()
        isfile(tmp) && rm(tmp; force = true)
        return (status = :skipped, value = missing,
                reason = "unreachable: $(sprint(showerror, e))")
    end

    got = file_sha256(tmp)

    # --- Bootstrap: no trusted baseline hash yet ⇒ record without comparing, commit atomically.
    if expected_sha256 === nothing || isempty(expected_sha256)
        mv(tmp, dest; force = true)
        return (status = :present, sha256 = got, bootstrapped = true)
    end

    # --- Path B (hard abort): the integrity check. DELIBERATELY OUTSIDE any try/catch — a
    # content-hash MISMATCH is tamper/corruption and a HARD abort (D-05), never a skip.
    if got != lowercase(expected_sha256)
        rm(tmp; force = true)
        error("content-hash MISMATCH for $url: expected $(expected_sha256), got $got")
    end

    mv(tmp, dest; force = true)   # atomic commit (filesystem rename)
    return (status = :present, sha256 = got)
end
