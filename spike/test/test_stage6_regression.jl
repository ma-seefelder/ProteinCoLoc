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

# spike/test/test_stage6_regression.jl --- D-10 stage-6 regression (SC1c / SC1d / SC1e).
#
# EXACT EQUALITY IS THE PRE-REGISTERED CRITERION AND IT IS NOT NEGOTIABLE. `P11_STAGE6_EXACT`
# IS `true` IN THE FROZEN TIER-1 PRE-REGISTRATION, SO EVERY COMPARISON BELOW USES `==` ON THE
# RAW Float64 BYTES. BIT-FOR-BIT EQUALITY WAS *MEASURED* TO HOLD ACROSS TWELVE INDEPENDENT
# PHILOX KEYS AND FOUR SHIFT PAIRS BEFORE THE EDIT LANDED (11-RESEARCH.md §A3), AND THE REASON
# IS STRUCTURAL, NOT LUCKY: `recenter(LinearMap(I), c)` COMPOSES TO AN EXACT IDENTITY LINEAR
# PART WITH A ZERO TRANSLATION, SO THE COMPOSED MAP PRODUCES THE IDENTICAL SOURCE COORDINATES
# AND HENCE THE IDENTICAL B-SPLINE INTERPOLATION. A FAILURE HERE IS THEREFORE A GENUINE SIGNAL
# -- FOR EXAMPLE A StaticArrays / Interpolations REDUCTION-ORDER CHANGE -- AND MUST BE
# INVESTIGATED, NOT SMOOTHED AWAY. THE *ONLY* SANCTIONED ESCAPE HATCH IS THE PRE-REGISTERED
# CONSTANT `P11_STAGE6_TOLERANCE_FALLBACK` (1e-13, p11_consts.jl §9), AND ANY USE OF IT MUST BE
# RECORDED AS A NAMED DEVIATION IN THE PHASE-11 PROVENANCE DOCUMENT AND THE PHASE-11 REPORT.
# DO NOT SOFTEN A FAILING ASSERTION TO A TOLERANCE COMPARISON IN THIS FILE.
#
# What the three success criteria mean here:
#   SC1c -- the stage-6 transform collapses to ONE `AffineMap`, so `warp` runs exactly ONE
#           interpolation pass (D-10: a second pass would smooth, and smoothing decorrelates,
#           widening the very posterior SC2 measures for a reason unrelated to misalignment).
#   SC1d -- at `chromatic_eps = 0` the composed warp reproduces the legacy translation-only
#           warp, and the full `simulate_pair` pipeline reproduces the PRE-EDIT golden bytes
#           captured at PHASE11_BASE_SHA before a single stage-6 byte moved.
#   SC1e -- a legacy 7-field θ (what `src/amortized/ood.jl` `_theta_tuple` builds from a 7-row
#           posterior mean) is still accepted and still returns those same bytes.
#
# SEED DISCIPLINE (D-01): this file rides `p11_rng(P11_FIXTURE_COUNTER + k)` only -- the
# reserved FIXTURE counter, never a counter a reported number consumes.
#
# CPU-only (D-10), spike-local; reaches `src/` not at all. Run:
#     julia --project=spike spike/test/test_stage6_regression.jl

using Test
using JLD2
using Random
using Random123
using CoordinateTransformations      # Translation, LinearMap, recenter, AffineMap
using ImageTransformations           # warp
using Interpolations                 # BSpline, Linear

# --- ORDER MATTERS: the pre-registered stream and the stage-6 contract constants must exist
#     before the prior and the simulator that ride them are loaded. Guarded for idempotency (S2).
isdefined(@__MODULE__, :P11_DEV_SEED) ||
    include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
isdefined(@__MODULE__, :sample_prior) ||
    include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair) ||
    include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))

const P11_GOLDEN_FIXTURE = joinpath(@__DIR__, "fixtures", "p11_stage6_golden.jld2")
const FORWARD_SRC_PATH   = joinpath(@__DIR__, "..", "simulator", "forward.jl")

# The four shift pairs §A3 measured exact equality across: the null shift, the golden θ's own
# sub-pixel pair, a >1 px pair, and the widened-prior (D-02) endpoint pair.
const SHIFT_PAIRS = ((0.0, 0.0), (0.37, -0.82), (2.5, 1.25), (-3.0, 3.0))

# Source text with whole-line comments stripped, so a `warp(` mentioned in prose cannot be
# mistaken for a call. The one-interpolation-pass requirement is a property of the SOURCE
# (how many times `warp` is invoked), so it is asserted at source level, not by instrumentation.
_strip_comment_lines(src::AbstractString) =
    join(filter(l -> !startswith(strip(l), "#"), split(src, '\n')), '\n')

@testset "D-10 stage-6 regression" verbose = true begin

    @testset "SC1c: the stage-6 transform is a SINGLE AffineMap" begin
        # Exactly the composition stage 6 builds, with the identity linear part the
        # chromatic_eps = 0 case produces and the geometric centre of a 256² field.
        A = Translation(0.37, -0.82) ∘ recenter(LinearMap([1.0 0.0; 0.0 1.0]), (128.5, 128.5))
        @test A isa CoordinateTransformations.AffineMap
        # If the composition did NOT collapse, `warp` would still accept it -- and would still
        # run one pass -- but the map would be a lazily-composed object, which is the shape a
        # future edit could silently turn into two passes. Pin the collapsed form.
        @test !(A isa CoordinateTransformations.ComposedTransformation)

        # ONE resampling pass, asserted on the simulator source itself.
        code = _strip_comment_lines(read(FORWARD_SRC_PATH, String))
        @test count("warp(", code) == 1
    end

    # A fixed field standing in for a stage-5 output: strictly positive, background-floored,
    # explicit-RNG-constructor (never `seed!` global state).
    X = rand(Random.Xoshiro(2026), 256, 256) .+ BG_FLOOR

    @testset "SC1d(a): composed warp at chromatic_eps = 0 == legacy Translation warp" begin
        # `c` is DERIVED from `axes(X)` exactly as stage 6 derives it -- never a literal and
        # never `size(X)/2`, both of which are off by half a pixel for 1-based axes.
        c = map(ax -> (first(ax) + last(ax)) / 2, axes(X))
        for (dy, dx) in SHIFT_PAIRS
            composed = warp(X, Translation(dy, dx) ∘ recenter(LinearMap([1.0 0.0; 0.0 1.0]), c),
                            axes(X); method = BSpline(Linear()), fillvalue = BG_FLOOR)
            legacy   = warp(X, Translation(dy, dx),
                            axes(X); method = BSpline(Linear()), fillvalue = BG_FLOOR)
            @test collect(composed) == collect(legacy)     # exact `==` (P11_STAGE6_EXACT)
        end
    end

    # --- The PRE-EDIT golden bytes -------------------------------------------------------
    golden        = JLD2.load(P11_GOLDEN_FIXTURE)
    golden_theta  = golden["theta"]
    golden_imsize = golden["imsize"]
    golden_chans  = golden["channels"]
    golden_keys   = golden["keys_used"]

    # `merge` OVERWRITES an existing field in place and APPENDS a new one at the end, so
    # `merge(golden_theta, (chromatic_eps = 0.0,))` yields exactly the 8-field θ layout
    # `sample_prior` now produces -- chromatic_eps last -- preserving `collect(values(θ))`
    # row order for rows 1..7. That is the invariant `test_p11_forced_theta.jl` pins.
    eps0_theta = merge(golden_theta, (chromatic_eps = 0.0,))

    @testset "SC1d(b): full simulate_pair == the PRE-EDIT golden fixture" begin
        # A fixture REGENERATED after the ε edit would carry the 8th field and would compare
        # the new code against itself, passing vacuously. Refuse that fixture outright.
        @test !hasproperty(golden_theta, :chromatic_eps)
        @test length(golden_theta) == 7
        @test last(keys(golden_theta)) == :noise

        sha = golden["phase11_base_sha"]
        @test sha isa AbstractString
        @test length(sha) == 40

        # The fixture must have ridden the reserved FIXTURE counter, not a reported one.
        @test golden["fixture_counter"] == P11_FIXTURE_COUNTER
        @test golden["seed"] == UInt64(P11_DEV_SEED)
        @test golden["salt"] == UInt64(P11_SALT)

        @test P11_STAGE6_EXACT              # the contract this testset scores against
        @test length(golden_chans) == length(golden_keys)

        for (i, k) in enumerate(golden_keys)
            out = simulate_pair(p11_rng(P11_FIXTURE_COUNTER + k), eps0_theta;
                                imsize = golden_imsize)
            @test out[1] == golden_chans[i][1]     # exact `==`, both channels
            @test out[2] == golden_chans[i][2]
        end
    end

    @testset "SC1e: backward compatibility with a legacy 7-field theta" begin
        # The stored θ is passed UNMODIFIED: no chromatic_eps field at all. This is what the
        # defensive `hasproperty` read in simulate_pair means -- an absent field is exactly
        # "no chromatic aberration", not an error and not a different answer.
        for (i, k) in enumerate(golden_keys)
            legacy = simulate_pair(p11_rng(P11_FIXTURE_COUNTER + k), golden_theta;
                                   imsize = golden_imsize)
            eps0   = simulate_pair(p11_rng(P11_FIXTURE_COUNTER + k), eps0_theta;
                                   imsize = golden_imsize)
            @test legacy[1] == golden_chans[i][1]
            @test legacy[2] == golden_chans[i][2]
            @test legacy[1] == eps0[1]
            @test legacy[2] == eps0[2]
        end
    end

    @testset "entry guards reject a hostile chromatic_eps" begin
        # BOTH directions (the test_sbc.jl:63-68 discipline): a guard that always threw would
        # pass the rejection clauses and be caught by the acceptance clauses below.
        @test_throws ArgumentError simulate_pair(Random.Xoshiro(1),
            merge(golden_theta, (chromatic_eps = NaN,)); imsize = golden_imsize)
        @test_throws ArgumentError simulate_pair(Random.Xoshiro(1),
            merge(golden_theta, (chromatic_eps = Inf,)); imsize = golden_imsize)
        # -1.0 is the collapse boundary: the backward scale is 1/(1 + chromatic_eps), so at
        # -1 it divides by zero and below -1 the affine map MIRRORS channel 2.
        @test_throws ArgumentError simulate_pair(Random.Xoshiro(1),
            merge(golden_theta, (chromatic_eps = -1.0,)); imsize = golden_imsize)

        for e in (-0.5, 0.5)
            out = simulate_pair(Random.Xoshiro(1),
                merge(golden_theta, (chromatic_eps = e,)); imsize = golden_imsize)
            @test out isa Vector{Matrix{Float64}}
            @test length(out) == 2
            @test all(c -> all(isfinite, c), out)
        end
    end

    @testset "stage-6 regression ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
