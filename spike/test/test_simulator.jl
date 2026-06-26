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

# spike/test/test_simulator.jl --- Phase-2 Wave-0 scaffold: the SIM-03 summary
# contract proven on a synthetic image BEFORE any physics simulator exists, so
# Wave-2 builds against a known-good, frozen summary contract.
#
# The contract is the UNCHANGED src/ summary pipeline reached read-only via
# spike/contract.jl (include() coupling, NOTES §4):
#     2-channel image -> build_mci -> patch()/correlation() -> 8x8 rho -> induced mu
#
# Traps encoded as assertions (src/colocalization.jl):
#   * background-not-zero: _exclude_zero drops 0.0/NaN/missing (line 197); a
#     strictly-positive image keeps every pixel, so all 64 patches survive.
#   * >=15-px patch floor: correlation sets missing when <=15 survivors (line 235).
#     At 256^2 / 8x8 each patch is 32x32 = 1024 px, trivially above the floor.
#   * negative control: a hard-zero background region drops those pixels, starving
#     whole patches below the floor -> materially more `missing` (proves the trap).

using Test
using Random

# include() coupling boundary -> src/ summary functions, read-only (Task 2 impl).
include(joinpath(@__DIR__, "..", "contract.jl"))

@testset "SIM-03 summary contract" verbose = true begin

    # --- Synthetic 2-channel, strictly-positive, correlated image -------------
    # Shared smoothed positive field couples the channels; per-channel positive
    # jitter + a +1.0 baseline guarantees every intensity > 0 (no hard zeros).
    rng   = MersenneTwister(2026)
    field = imfilter(abs.(randn(rng, 256, 256)), Kernel.gaussian(3))
    ch1   = Matrix{Float64}(field .+ 0.5 .* abs.(randn(rng, 256, 256)) .+ 1.0)
    ch2   = Matrix{Float64}(field .+ 0.5 .* abs.(randn(rng, 256, 256)) .+ 1.0)
    data  = [ch1, ch2]

    mci = build_mci(data)

    @testset "MultiChannelImage built per D-11 (pixel_size = pixel-dim tuple)" begin
        @test mci isa MultiChannelImage
        @test mci.pixel_size == size(data[1])        # (256,256) pixels, NOT micrometers
        @test mci.pixel_size isa Tuple{Int, Int}
        @test length(mci.channels) == 2
        @test length(mci.otsu_threshold) == 2
        @test mci.channels == ["ch1", "ch2"]
    end

    rho = summary(mci)

    @testset "8x8 contract: <=64 finite rho, background-not-zero trap respected" begin
        @test size(rho) == (8, 8)
        @test eltype(rho) == Union{Float64, Missing}
        @test count(!ismissing, rho) <= 64               # SIM-03: at most 64 per-patch rho
        @test all(isfinite, skipmissing(rho))            # every present rho is finite
        @test count(ismissing, rho) <= 2                 # strictly-positive image: few/zero missing
        @test isfinite(induced_mu(mci))                  # induced mu is a finite Float64
        @test induced_mu(mci) isa Float64
    end

    @testset "negative control: hard-zero background -> more missing (trap matters)" begin
        ch1z = copy(ch1); ch2z = copy(ch2)
        ch1z[1:200, :] .= 0.0                            # zero most rows of both channels
        ch2z[1:200, :] .= 0.0
        mciz = build_mci([ch1z, ch2z]; name = "zero-bg")
        rhoz = summary(mciz)
        @test count(ismissing, rhoz) > count(ismissing, rho)
    end

end
