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

# corpus/test/test_load.jl --- two-channel TIFF -> MultiChannelImage conversion testset (SC1, D-03).
#
# Proves the D-03 conversion path on the committed fixture: two single-channel TIFFs load via the
# FROZEN src/ constructor into a MultiChannelImage and reduce through the UNCHANGED patch/correlation
# summary path to a FINITE (non-all-missing) 8x8 grid, with the Pitfall-6 has_signal guard true.
# Runs fully OFFLINE against the committed fixture.

# Guarded include so this testset runs standalone or wired into runtests.jl.
isdefined(@__MODULE__, :load_anchor)  || include(joinpath(@__DIR__, "..", "load.jl"))
isdefined(@__MODULE__, :make_fixture) || include(joinpath(@__DIR__, "make_fixture.jl"))

@testset "load: two-channel TIFF -> MultiChannelImage (SC1/D-03)" begin
    # Ensure the committed fixture exists (regenerate deterministically on a fresh checkout).
    fp1 = joinpath(@__DIR__, "fixtures", "fixture_ch1.tif")
    fp2 = joinpath(@__DIR__, "fixtures", "fixture_ch2.tif")
    if !(isfile(fp1) && isfile(fp2))
        make_fixture(joinpath(@__DIR__, "fixtures"))
    end

    # (a) load_anchor builds a MultiChannelImage from the two single-channel TIFFs.
    img = load_anchor("fixture_anchor", fp1, fp2)
    @test img isa MultiChannelImage
    @test length(img.data) == 2
    @test img.data[1] isa Matrix{Float64}
    @test img.channels == ["ch1", "ch2"]

    # (b) summary_grid reduces through the UNCHANGED src/ path to a finite 8x8 grid.
    grid = summary_grid(img; num_patches = 8)
    @test size(grid) == (8, 8)
    @test count(!ismissing, grid) > 0          # not an all-missing grid (SC1)
    @test isfinite(mean(skipmissing(grid)))

    # (c) has_signal Pitfall-6 guard is true for the committed fixture.
    @test has_signal(grid) === true

    # (d) has_signal returns false on a starved (all-missing) grid.
    @test has_signal(fill(missing, 8, 8)) === false
end
