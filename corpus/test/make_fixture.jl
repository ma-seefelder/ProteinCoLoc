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

# corpus/test/make_fixture.jl --- deterministic synthetic two-channel TIFF fixture (SC1, D-03).
#
# Generates a tiny, network-free, byte-reproducible two-channel fixture that stands in for a
# real physical anchor so the whole corpus test suite runs OFFLINE. The two channels share a
# smooth positive signal field (strong positive correlation) with per-channel jitter, and EVERY
# pixel is strictly > 0 so that every 8x8 patch clears the src/colocalization.jl `_exclude_zero`
# >15-valid-pixel floor (Pitfall 6): the fixture never degenerates to an all-`missing` grid.
#
# DECOUPLING (hard constraint, CLAUDE.md): corpus-local test asset. Uses Images (already a root
# dep) read-only and Random123 for keyed reproducibility; touches no src/, no root Project.toml.

import Images
using Random123

# Load the pre-declared corpus consts (CORPUS_MASTER_SEED) if not already in scope, so this file
# works both standalone and after the runtests.jl include of ../config.jl.
isdefined(@__MODULE__, :CORPUS_MASTER_SEED) || include(joinpath(@__DIR__, "..", "config.jl"))

"""
    make_fixture(dir; seed = CORPUS_MASTER_SEED) -> (path_ch1, path_ch2)

Deterministically write two single-channel 64x64 `Float64` TIFFs (`fixture_ch1.tif`,
`fixture_ch2.tif`) into `dir` and return their paths.

The channels are built from a shared smooth field `s(i,j)` that varies WITHIN every 8x8 patch
(so per-patch variance is non-zero and Pearson correlation is finite) plus small per-channel
jitter drawn from a Random123 Philox stream keyed by `seed`. All intensities lie strictly inside
`(0, 1)`, so no pixel is dropped by `_exclude_zero` and every patch retains all 64 valid pixels.
Two invocations with the same `seed` produce byte-identical TIFFs.
"""
function make_fixture(dir::AbstractString; seed = CORPUS_MASTER_SEED)
    isdir(dir) || mkpath(dir)
    n = 64
    rng = Philox4x(UInt64, (UInt64(seed), UInt64(0)))

    ch1 = Matrix{Float64}(undef, n, n)
    ch2 = Matrix{Float64}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        # Shared smooth field, strictly positive, period-16 so it varies across each 8x8 patch.
        s  = 0.5 + 0.25 * sin(2π * i / 16) * cos(2π * j / 16)   # ∈ (0.25, 0.75)
        e1 = 0.05 * (rand(rng, Float64) - 0.5)                  # per-channel jitter
        e2 = 0.05 * (rand(rng, Float64) - 0.5)
        ch1[i, j] = clamp(s + e1, 0.05, 0.95)
        ch2[i, j] = clamp(0.9 * s + 0.05 + e2, 0.05, 0.95)      # correlated 2nd channel
    end

    path_ch1 = joinpath(dir, "fixture_ch1.tif")
    path_ch2 = joinpath(dir, "fixture_ch2.tif")
    Images.save(path_ch1, Images.Gray.(ch1))
    Images.save(path_ch2, Images.Gray.(ch2))
    return (path_ch1, path_ch2)
end
