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

# corpus/test/test_hash.jl --- SHA-256 content-hash testset (SC1, D-04, D-08).
#
# Proves the download-integrity primitives are a DETERMINISTIC cryptographic digest:
# same bytes -> same hash, one flipped byte -> different hash, key order never perturbs
# the canonical config, and subhashes are full 64-char SHA-256 (never a 16-char Base
# hash). Runs fully OFFLINE against the committed fixture.

# Guarded include so this testset runs standalone or wired into runtests.jl.
isdefined(@__MODULE__, :file_sha256) || include(joinpath(@__DIR__, "..", "hash.jl"))

@testset "hashing (SC1/D-04/D-08)" begin
    fx = joinpath(@__DIR__, "fixtures", "fixture_ch1.tif")

    # (a) file_sha256 equals the one-shot reference and is deterministic + lowercase hex.
    ref = bytes2hex(SHA.sha256(read(fx)))
    h1  = file_sha256(fx)
    @test h1 == ref
    @test file_sha256(fx) == h1            # deterministic across repeated calls
    @test length(h1) == 64
    @test h1 == lowercase(h1)
    @test all(ch -> ch in "0123456789abcdef", h1)

    # (b) same bytes hash equal; a single flipped byte changes the hash.
    d = mktempdir()
    a = joinpath(d, "a.bin"); b = joinpath(d, "b.bin"); c = joinpath(d, "c.bin")
    raw = UInt8[0x01, 0x02, 0x03, 0x04, 0x05]
    write(a, raw); write(b, copy(raw))
    flipped = copy(raw); flipped[3] ⊻= 0xff
    write(c, flipped)
    @test file_sha256(a) == file_sha256(b)
    @test file_sha256(a) != file_sha256(c)

    # (c) canonical_config is stable across key insertion order.
    nt1 = (alpha = 1, beta = "x", gamma = 3.0)
    nt2 = (gamma = 3.0, alpha = 1, beta = "x")
    @test canonical_config(nt1) == canonical_config(nt2)
    @test occursin("alpha=", canonical_config(nt1))

    # (d) manifest_hash: fixed-order determinism, order-sensitivity, config-sensitivity.
    m1 = manifest_hash(String[a, c], canonical_config(nt1))
    @test length(m1) == 64
    @test m1 == manifest_hash(String[a, c], canonical_config(nt1))          # deterministic
    @test m1 != manifest_hash(String[c, a], canonical_config(nt1))          # order-sensitive
    @test m1 != manifest_hash(String[a, c], canonical_config((alpha = 2,))) # config-sensitive

    # (e) subhashes are full 64-char SHA-256 hex (NOT a 16-char Base hash) and match file_sha256.
    sh = subhashes(String[a, c])
    @test length(sh) == 2
    @test all(v -> length(v) == 64, values(sh))
    @test sh[a] == file_sha256(a)
    @test sh[c] == file_sha256(c)
end
