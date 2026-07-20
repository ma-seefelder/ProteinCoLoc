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

# corpus/test/test_fetch.jl --- fetch_verified offline testset (D-04, D-05).
#
# Proves the D-05 ASYMMETRY structurally, with ZERO network: the network call is isolated
# in the overridable `_DOWNLOAD_HOOK` seam, which these tests swap for a local-copy stub.
#   * reachable + correct hash  -> (:present) atomic commit, no .tmp left behind
#   * reachable + bootstrap      -> (:present, bootstrapped) records hash without comparing
#   * reachable + wrong hash     -> HARD error (ErrorException) ABORT, .tmp removed
#   * unreachable                -> (:skipped, reason=...) NEVER throws (offline/CI stays green)

# Guarded include so this testset runs standalone or wired into runtests.jl.
isdefined(@__MODULE__, :fetch_verified) || include(joinpath(@__DIR__, "..", "fetch.jl"))

@testset "fetch_verified (D-04/D-05)" begin
    fx   = joinpath(@__DIR__, "fixtures", "fixture_ch1.tif")
    good = file_sha256(fx)

    orig_hook = _DOWNLOAD_HOOK[]
    try
        # Stub the network seam: "download" == copy the committed fixture to dest.
        _DOWNLOAD_HOOK[] = (url, dest) -> cp(fx, dest; force = true)

        # (a) success — correct expected hash commits atomically, no .tmp remains.
        d1   = mktempdir()
        dest = joinpath(d1, "sub", "out.tif")     # nested dir exercises mkpath(dirname(dest))
        r = fetch_verified("stub://ok", dest, good)
        @test r.status == :present
        @test r.sha256 == good
        @test isfile(dest)
        @test !isfile(dest * ".tmp")

        # Uppercase expected hash still matches (case-insensitive compare via lowercase()).
        d1u   = mktempdir()
        destu = joinpath(d1u, "up.tif")
        ru = fetch_verified("stub://ok", destu, uppercase(good))
        @test ru.status == :present
        @test isfile(destu)

        # (b) bootstrap — nothing / empty expected hash records without comparing.
        d2 = mktempdir()
        dest2 = joinpath(d2, "boot.tif")
        r2 = fetch_verified("stub://boot", dest2, nothing)
        @test r2.status == :present
        @test r2.sha256 == good
        @test r2.bootstrapped
        @test isfile(dest2)
        @test !isfile(dest2 * ".tmp")

        r2e = fetch_verified("stub://boot", joinpath(d2, "boot2.tif"), "")
        @test r2e.status == :present
        @test r2e.bootstrapped

        # (c) hard-fail — wrong expected hash ABORTS (ErrorException); .tmp and dest removed.
        d3 = mktempdir()
        dest3 = joinpath(d3, "bad.tif")
        wrong = "0"^64
        @test_throws ErrorException fetch_verified("stub://bad", dest3, wrong)
        @test !isfile(dest3 * ".tmp")
        @test !isfile(dest3)

        # (d) skip-with-flag — an unreachable source returns :skipped and NEVER throws.
        _DOWNLOAD_HOOK[] = (url, dest) -> error("simulated network failure")
        d4 = mktempdir()
        dest4 = joinpath(d4, "skip.tif")
        r4 = fetch_verified("stub://unreachable", dest4, good)
        @test r4.status == :skipped
        @test occursin("unreachable", r4.reason)
        @test !isfile(dest4 * ".tmp")
        @test !isfile(dest4)
    finally
        _DOWNLOAD_HOOK[] = orig_hook   # restore the real network seam
    end
end
