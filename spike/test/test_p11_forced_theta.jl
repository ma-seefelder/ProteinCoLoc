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

# spike/test/test_p11_forced_theta.jl --- the θ row-alignment invariant (D-15, D-09).
#
# `collect(values(θ))` is the house θ→vector map, so the NamedTuple's FIELD ORDER is a
# load-bearing contract that nothing in the compiler checks. Two separate pieces of Phase-11
# design rest on it:
#
#   D-15 forces a fixed misalignment with `merge(sample_prior(rng), forced)`, which is only
#        safe because `merge` overwrites an EXISTING key IN PLACE rather than moving it;
#   D-09 appends `chromatic_eps` as the 8th θ column, which is only safe because `merge`
#        appends a NEW key AT THE END, leaving rows 1..7 (and every positional read of them)
#        untouched.
#
# Both facts are asserted here, EXACTLY (no tolerance): a row-order claim that only holds
# approximately is not a row-order claim.

using Test
using Random123

# Units under test. Guarded for idempotency (S2); consts first, then the prior.
isdefined(@__MODULE__, :P11_DEV_SEED) ||
    include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
isdefined(@__MODULE__, :sample_prior) ||
    include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))

# The fixture rides the FIXTURE counter, never a counter a reported number consumes (D-01).
const P11_FT_THETA = sample_prior(p11_rng(P11_FIXTURE_COUNTER))

@testset "P11 forced-θ row alignment (D-15, D-09)" verbose = true begin

    @testset "NamedTuple merge field-order invariant" begin
        # An EXISTING key is overwritten IN PLACE — the field order does not move. This is
        # what makes D-15's forced-θ injection safe.
        @test keys(merge((a = 1, b = 2, c = 3), (b = 9,))) == (:a, :b, :c)
        @test values(merge((a = 1, b = 2, c = 3), (b = 9,))) == (1, 9, 3)
        # A NEW key APPENDS AT THE END — every prior positional index is preserved. This is
        # what makes D-09's append-at-end θ layout safe.
        @test keys(merge((a = 1, b = 2), (z = 3,))) == (:a, :b, :z)
        @test values(merge((a = 1, b = 2), (z = 3,))) == (1, 2, 3)
        # Both at once: overwrite in place AND append, in one call.
        m = merge((a = 1, b = 2), (b = 9, z = 3))
        @test keys(m)   == (:a, :b, :z)
        @test values(m) == (1, 9, 3)
    end

    @testset "sample_prior row alignment under a forced subset" begin
        θ = P11_FT_THETA
        @test keys(θ) == (:ρ_true, :spillover, :autofluorescence, :label_efficiency,
                          :shift_dx, :shift_dy, :noise)

        forced = (shift_dx = 1.25, shift_dy = -0.5)
        θf = merge(θ, forced)
        # Forcing a SUBSET must not disturb the key order at all.
        @test keys(θf) == keys(θ)
        @test length(θf) == length(θ)
        # ... and the forced values land on rows 5 and 6 of the house θ→vector map.
        @test collect(values(θf))[5:6] == [1.25, -0.5]
        # Every other row is bit-identical to the unforced draw.
        @test collect(values(θf))[1:4] == collect(values(θ))[1:4]
        @test collect(values(θf))[7]   == collect(values(θ))[7]
        # The forced fields are readable by name too, so a positional and a named read agree.
        @test θf.shift_dx == 1.25
        @test θf.shift_dy == -0.5
    end

    @testset "θ arity is derived, never literal" begin
        θ = P11_FT_THETA
        @test length(collect(values(θ))) == length(keys(θ))
        # TODAY the last field is `noise`. Asserting that here makes the append-at-end
        # contract for `chromatic_eps` checkable BEFORE and AFTER the θ-arity edit: after it,
        # `last(keys(θ))` becomes `:chromatic_eps` and rows 1..7 must be unchanged.
        @test last(keys(θ)) == :noise
        # An appended 8th field leaves the first seven rows exactly where they were.
        θ8 = merge(θ, (chromatic_eps = 0.0,))
        @test length(keys(θ8)) == length(keys(θ)) + 1
        @test last(keys(θ8)) == :chromatic_eps
        @test collect(values(θ8))[1:length(keys(θ))] == collect(values(θ))
        # The draw is deterministic under the fixed DEV stream (D-01 reproducibility).
        @test collect(values(sample_prior(p11_rng(P11_FIXTURE_COUNTER)))) == collect(values(θ))
    end

    @testset "P11 forced-θ ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
