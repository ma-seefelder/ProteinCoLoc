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

# spike/test/test_p13_alpha.jl --- the D-16 alpha-ladder invariant gate (D-15, D-16).
#
# WHAT THIS FILE IS FOR. The alpha ladder is the ONLY graded segregation evidence Phase 13 will
# ever have (the corpus holds one segregated anchor, it is sealed for Phase 16, and there is no
# graded series in it at all), so it is load-bearing evidence -- and it has one dominant failure
# mode that makes it SILENTLY useless rather than visibly broken.
#
#   PITFALL 2 -- the trap. `_exclude_zero` (src/colocalization.jl:154-169) DROPS any pixel where
#   EITHER channel is zero. So zeroing ch2 inside the ch1 object mask does NOT anti-correlate
#   those pixels, it DELETES them: the correlation is then measured over the untouched
#   complement, object-dense patches can fall below the 15-survivor floor and go `missing`, and
#   the ladder reads FLAT. A flat ladder looks exactly like the finding "the net cannot see
#   spatial segregation" while the net was never shown any. The strictly-positive background
#   floor and the intensity conservation below are what stand between this phase and that
#   conclusion, so they are asserted here rather than trusted.
#
#   PITFALL 2b -- the trap in the OTHER direction. The naive repair of Pitfall 2 is to assert
#   that every output pixel is strictly positive. That property is FALSE on real microscopy
#   before any alpha is applied at all, and asserting it turns a CORRECT algorithm into a red
#   suite at alpha = 0, where the output IS the input. The pre-registered rule is
#   `P13_ALPHA_ZERO_INVARIANT = :count_iszero_preserved` -- "alpha introduces no NEW zero".
#
# Mirrors test_bf.jl / test_p13_labels.jl: license header -> using Test -> guarded include of the
# unit under test -> fixtures computed ONCE at module level -> one outer @testset. No training,
# no figure call (figures are script artifacts and belong to plans 13-13 / 13-16). CPU-only
# (D-10), spike-local, touches no src/ byte and never reaches the sealed corpus.

using Test
using Statistics

# Unit under test: the alpha transform (pulls p13/consts.jl -> contract.jl -> data/encode.jl).
isdefined(@__MODULE__, :alpha_segregate) ||
    include(joinpath(@__DIR__, "..", "p13", "alpha_series.jl"))
# The fixture generator. Guarded separately because alpha_series.jl deliberately does NOT depend
# on the simulator -- it is substrate-agnostic (D-15), and only this suite needs a pair to run on.
isdefined(@__MODULE__, :simulate_pair) ||
    include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))

# --- Fixtures computed ONCE (the test_bf.jl:49-65 idiom) --------------------------------------
# ALL fixtures ride p13_fix_rng(P13_FIXTURE_COUNTER), so this gate NEVER consumes (and so never
# pre-observes) a REPORTED Phase-13 stream (D-01). P13_FIX_SEED != P13_DEV_SEED is asserted
# inside the suite as the explicit statement of that rule.
#
# theta is written out rather than drawn from sample_prior so the fixture is a FIXED, readable
# pair rather than a prior draw that shifts if the prior is ever touched: rho_true = 0.0 makes
# the alpha = 0 rung a genuinely RANDOM pair BY CONSTRUCTION, which is the simulated substrate's
# distinguishing property (on real substrate alpha = 0 is moderately COLOCALIZED, measured
# m-bar = +0.329 / +0.248 -- the two arms are never averaged). The nuisances are mild and the
# shifts are exactly zero so the ladder is not confounded with mis-registration.
const ALPHA_FIX_THETA = (ρ_true           = 0.0,
                         spillover        = 0.05,
                         autofluorescence = 0.02,
                         label_efficiency = 0.9,
                         shift_dx         = 0.0,
                         shift_dy         = 0.0,
                         noise            = 0.3,
                         chromatic_eps    = 0.0)

# 128x128, not 64x64. At 64x64 the Otsu mask of this theta covers only 0.66% of the frame, BELOW
# the pre-registered lower bound of P13_ALPHA_MASK_FRACTION_BOUNDS -- so the guard correctly
# refuses it, and a fixture that trips a pre-registered guard is a bad fixture, not a reason to
# move the bound. At 128x128 the mask covers 28.6%, comfortably in band, and the whole file runs
# in a few seconds.
const ALPHA_FIX_IMSIZE = (128, 128)
const ALPHA_FIX_PAIR   = simulate_pair(p13_fix_rng(P13_FIXTURE_COUNTER), ALPHA_FIX_THETA;
                                       imsize = ALPHA_FIX_IMSIZE)
const ALPHA_FIX_X      = ALPHA_FIX_PAIR[1]
const ALPHA_FIX_Y      = ALPHA_FIX_PAIR[2]

# The two ladder-level constants, computed ONCE from the UNMODIFIED pair -- exactly as
# alpha_ladder does internally, so the direct alpha_segregate calls below use the same pair of
# constants the ladder does and every assertion is about the same experiment.
const ALPHA_FIX_MASK  = ch1_object_mask(ALPHA_FIX_PAIR)
const ALPHA_FIX_FLOOR = alpha_background_floor(ALPHA_FIX_Y)

# THE DECLARED DYNAMIC RANGE IS A PROPERTY OF THE SUBSTRATE, NOT OF THE ALGORITHM (D-15).
# `P13_ALPHA_MAX_VALUE_BOUND` = 1.0 was pre-registered against the real Gray TIFFs, whose values
# live in [0,1] and whose measured post-boost maximum stayed at or below 0.75. The SIMULATOR's
# intensities are unnormalized softplus outputs and reach 4.84 on this fixture BEFORE any alpha,
# so 1.0 is not this substrate's range and the simulated arm declares none -- the realized
# maximum is still RECORDED (`realized_max_value`), never clamped.
#
# The pre-registered bound is therefore exercised on a SECOND fixture: the same pair rescaled by
# one positive scalar into the Gray range with headroom (maximum 0.5, mirroring the real
# fixtures' regime). A positive affine rescale leaves every per-patch Pearson correlation
# EXACTLY unchanged, so this is the SAME ladder read on a substrate that HAS a declared range --
# which testset 3 checks explicitly rather than asserting. That is what lets the default code
# path (bound = 1.0) be tested for real instead of only through a hand-built toy.
const ALPHA_FIX_SCALE = 0.5 / maximum(maximum.(ALPHA_FIX_PAIR))
const ALPHA_GRAY_PAIR = [Matrix{Float64}(ch .* ALPHA_FIX_SCALE) for ch in ALPHA_FIX_PAIR]

# Verifier reports and per-rung readouts computed ONCE (each rebuilds the whole ladder).
const ALPHA_FIX_INV  = verify_alpha_invariants(ALPHA_FIX_PAIR; max_value = Inf)
const ALPHA_FIX_SUMS = alpha_ladder_summaries(ALPHA_FIX_PAIR; max_value = Inf)
const ALPHA_GRAY_INV = verify_alpha_invariants(ALPHA_GRAY_PAIR)
const ALPHA_GRAY_SUMS = alpha_ladder_summaries(ALPHA_GRAY_PAIR)
const ALPHA_FIX_LADDER = alpha_ladder(ALPHA_FIX_PAIR; mask = ALPHA_FIX_MASK,
                                      floor = ALPHA_FIX_FLOOR, max_value = Inf)

# Transform source with WHOLE-LINE COMMENTS STRIPPED, so the file's own prose -- which must NAME
# the sealed holdout in order to state that it is never opened -- cannot invalidate its own gate
# (the test_stage6_regression.jl:79-80 / test_p13_consts.jl:49-50 idiom).
const ALPHA_SRC  = read(joinpath(@__DIR__, "..", "p13", "alpha_series.jl"), String)
const ALPHA_CODE = join(filter(l -> !startswith(strip(l), "#"), split(ALPHA_SRC, '\n')), '\n')

@testset "P13 alpha-graded segregation series (D-15, D-16)" verbose = true begin

    @testset "fixtures never consume a reported stream (D-01)" begin
        @test P13_FIX_SEED != P13_DEV_SEED
        @test P13_FIXTURE_COUNTER ∉ (P13_TAU_COUNTER, P13_DATAGEN_COUNTER, P13_GATE_COUNTER,
                                     P13_ALPHA_COUNTER, P13_CONTINUITY_COUNTER)
    end

    @testset "alpha = 0 is BITWISE the input" begin
        # Exact `==`, never `isapprox`:
        @test alpha_segregate(ALPHA_FIX_X, ALPHA_FIX_Y, ALPHA_FIX_MASK, 0.0;
                              b = ALPHA_FIX_FLOOR, max_value = Inf) == ALPHA_FIX_Y
        # `isapprox` is DELIBERATELY not used here, mirroring the Phase-11 zero-epsilon
        # regression pattern. An approximate check would pass a transform that perturbs every
        # pixel by a rounding step -- and a perturbed zero set is a CHANGED summary, because
        # _exclude_zero keys on which pixels are exactly zero. The property is also structural,
        # not special-cased: alpha = 0 makes `removed` exactly 0.0 and `boost` exactly 1.0, so no
        # branch on a zero alpha exists in the transform.
        @test ALPHA_FIX_LADDER.pairs[1][2] == ALPHA_FIX_Y
        @test ALPHA_FIX_INV.bitwise_alpha0
        @test ALPHA_GRAY_INV.bitwise_alpha0
    end

    @testset "alpha introduces no NEW zero" begin
        # THE _exclude_zero GUARD (Pitfall 2). A violation means pixels were DELETED rather than
        # MOVED, which is the failure that makes the ladder read flat for a reason unrelated to
        # segregation.
        #
        # THE ASSERTION IS DELIBERATELY NOT AN ALL-PIXELS-STRICTLY-POSITIVE CHECK (Pitfall 2b).
        # That form is FALSE on the committed real TIFFs at alpha = 0, where the output IS the
        # input: positive_c2 carries 2 exact-zero pixels and negative_c2 carries 3, of 1 414 528,
        # in the SOURCE files before any alpha. It is false on this SIMULATED fixture too -- the
        # noise stage leaves 1 exact-zero pixel of 16 384 despite BG_FLOOR, which is recorded
        # here because 13-RESEARCH J5.3 expected the simulator substrate to be strictly positive.
        # Writing the naive form would fail a CORRECT algorithm on BOTH substrates. The
        # pre-registered rule is P13_ALPHA_ZERO_INVARIANT = :count_iszero_preserved, compared
        # against the count on the UNMODIFIED ch2.
        @test P13_ALPHA_ZERO_INVARIANT === :count_iszero_preserved
        src_zeros = count(iszero, ALPHA_FIX_Y)
        for a in P13_ALPHA_LADDER
            z = alpha_segregate(ALPHA_FIX_X, ALPHA_FIX_Y, ALPHA_FIX_MASK, a;
                                b = ALPHA_FIX_FLOOR, max_value = Inf)
            @test count(iszero, z) == src_zeros
        end
        @test ALPHA_FIX_INV.no_new_zeros
        @test ALPHA_FIX_INV.source_zero_count == src_zeros
        @test all(s.n_zero == src_zeros for s in ALPHA_FIX_SUMS)
    end

    @testset "the mask-fraction guard and the value bound fire" begin
        # Below the band (an empty mask) and above it (a mask covering the whole field). Outside
        # the band the complement cannot absorb the redistributed mass, the boost blows up and
        # the ladder is uninterpretable rather than merely noisy -- so it is a NAMED error, not a
        # warning.
        @test_throws ArgumentError alpha_segregate(zeros(16, 16), ones(16, 16), falses(16, 16),
                                                   0.5; b = 0.1)
        @test_throws ArgumentError alpha_segregate(zeros(16, 16), ones(16, 16), trues(16, 16),
                                                   0.5; b = 0.1)

        # On the fixture the realized fraction is IN band, and the verifier reports it rather
        # than only asserting it.
        @test first(P13_ALPHA_MASK_FRACTION_BOUNDS) <= ALPHA_FIX_INV.realized_mask_fraction <=
              last(P13_ALPHA_MASK_FRACTION_BOUNDS)
        @test ALPHA_FIX_INV.mask_fraction_ok
        @test ALPHA_FIX_INV.realized_mask_fraction == mean(ALPHA_FIX_MASK)

        # THE MAXIMUM IS RECORDED, NEVER CLAMPED. On the Gray-range fixture (the substrate the
        # pre-registered bound was measured on) every rung stays under P13_ALPHA_MAX_VALUE_BOUND
        # on the DEFAULT code path.
        @test ALPHA_GRAY_INV.max_value_ok
        @test all(s.max_value <= P13_ALPHA_MAX_VALUE_BOUND for s in ALPHA_GRAY_SUMS)
        @test ALPHA_GRAY_INV.realized_max_value ==
              maximum(s.max_value for s in ALPHA_GRAY_SUMS)

        # And the bound FIRES rather than clamping when the realized maximum exceeds it: the
        # unnormalized simulated fixture reaches 4.84 before any alpha, so the Gray-range default
        # rejects it outright instead of quietly rescaling it into range.
        @test maximum(ALPHA_FIX_Y) > P13_ALPHA_MAX_VALUE_BOUND
        @test_throws ArgumentError alpha_segregate(ALPHA_FIX_X, ALPHA_FIX_Y, ALPHA_FIX_MASK, 1.0;
                                                   b = ALPHA_FIX_FLOOR)
        # The error message RECORDS the realized value, so a reader learns the measurement from
        # the failure instead of only that it failed.
        err = try
            alpha_segregate(ALPHA_FIX_X, ALPHA_FIX_Y, ALPHA_FIX_MASK, 1.0; b = ALPHA_FIX_FLOOR)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("realized maximum", err)
        @test occursin("never clamped", err)
    end

    @testset "total ch2 intensity is conserved" begin
        # So an alpha-response can never be a brightness change wearing a segregation costume.
        for a in P13_ALPHA_LADDER
            z = alpha_segregate(ALPHA_FIX_X, ALPHA_FIX_Y, ALPHA_FIX_MASK, a;
                                b = ALPHA_FIX_FLOOR, max_value = Inf)
            @test isapprox(sum(z), sum(ALPHA_FIX_Y); rtol = 1e-10)
        end
        @test ALPHA_FIX_INV.intensity_conserved
        @test ALPHA_FIX_INV.max_rel_intensity_err <= P13_ALPHA_INVARIANT_RTOL
        @test ALPHA_GRAY_INV.max_rel_intensity_err <= P13_ALPHA_INVARIANT_RTOL
        @test all(isapprox(s.sum_ratio, 1.0; rtol = 1e-10) for s in ALPHA_FIX_SUMS)
    end

    @testset "ch1 is untouched and the mask is alpha-invariant" begin
        # A per-alpha threshold would make the RUNGS incomparable the same way src/local_map.jl:208
        # says a tile-local Otsu makes tiles incomparable -- and comparability of the rungs is the
        # entire content of the series. ch1 is never written, so re-deriving the mask from each
        # transformed pair is a genuine leak check: it fails the moment a future edit touches
        # channel 1.
        @test P13_ALPHA_MASK_RULE === :calculate_mask_ch1_once
        for rung in ALPHA_FIX_LADDER.pairs
            @test rung[1] == ALPHA_FIX_X
            @test ch1_object_mask(rung) == ALPHA_FIX_MASK
        end
        @test ALPHA_FIX_INV.mask_invariant
        @test ALPHA_GRAY_INV.mask_invariant
        @test ALPHA_FIX_LADDER.mask == ALPHA_FIX_MASK
        @test ALPHA_FIX_LADDER.floor == ALPHA_FIX_FLOOR
    end

    @testset "m-bar is monotone non-increasing and missing counts do not rise" begin
        mbars = [s.mbar for s in ALPHA_FIX_SUMS]
        @test length(mbars) == length(P13_ALPHA_LADDER)
        @test all(isfinite, mbars)
        for i in 2:length(mbars)
            @test mbars[i] <= mbars[i - 1] + 1e-9     # float-noise allowance, not a loosened bar
        end
        @test ALPHA_FIX_INV.mbar_monotone
        @test ALPHA_GRAY_INV.mbar_monotone
        # The ladder must actually MOVE, and on simulated substrate it must cross zero: a
        # constant curve would satisfy non-increase while carrying no evidence at all.
        @test first(mbars) > last(mbars)
        @test last(mbars) < 0.0

        # A RISING missing count is the Pitfall-2 warning sign that pixels are being deleted
        # rather than moved. It must not rise across the ladder, and on this fixture no patch is
        # missing at any rung.
        miss = [s.n_missing for s in ALPHA_FIX_SUMS]
        @test last(miss) <= first(miss)
        @test all(m -> m <= first(miss), miss)
        @test ALPHA_FIX_INV.missing_nonincreasing
        @test ALPHA_GRAY_INV.missing_nonincreasing
    end

    @testset "guards fire" begin
        @test_throws ArgumentError alpha_segregate(ALPHA_FIX_X, ALPHA_FIX_Y, ALPHA_FIX_MASK, 1.5;
                                                   b = ALPHA_FIX_FLOOR, max_value = Inf)
        @test_throws ArgumentError alpha_segregate(ALPHA_FIX_X, ALPHA_FIX_Y, ALPHA_FIX_MASK, -0.1;
                                                   b = ALPHA_FIX_FLOOR, max_value = Inf)
        # A zero floor is the Pitfall-2 trap itself, so it is refused at the boundary.
        @test_throws ArgumentError alpha_segregate(ALPHA_FIX_X, ALPHA_FIX_Y, ALPHA_FIX_MASK, 0.5;
                                                   b = 0.0, max_value = Inf)
        @test_throws ArgumentError alpha_segregate(ALPHA_FIX_X, ALPHA_FIX_Y, ALPHA_FIX_MASK, 0.5;
                                                   b = -1.0, max_value = Inf)
        # A full-field mask leaves nowhere to redistribute to.
        @test_throws ArgumentError alpha_segregate(zeros(16, 16), ones(16, 16), trues(16, 16),
                                                   0.5; b = 0.1, max_value = Inf)
        # Mismatched frames.
        @test_throws ArgumentError alpha_segregate(zeros(16, 16), ones(8, 8), falses(16, 16),
                                                   0.5; b = 0.1, max_value = Inf)
        # A single-channel input is not a pair.
        @test_throws ArgumentError alpha_ladder([ALPHA_FIX_X])
    end

    @testset "floor is strictly positive" begin
        @test ALPHA_FIX_FLOOR > 0.0
        @test alpha_background_floor(ALPHA_FIX_Y) == ALPHA_FIX_FLOOR   # deterministic
        # An image with more than 5% exact-zero pixels silently yields b = 0 and re-opens the
        # _exclude_zero trap, which is why this is a runtime error and not a comment.
        @test_throws ArgumentError alpha_background_floor(zeros(8, 8))
        @test_throws ArgumentError alpha_background_floor(ALPHA_FIX_Y; q = 0.0)
        @test_throws ArgumentError alpha_background_floor(ALPHA_FIX_Y; q = 1.0)
        # The pre-registered quantile, not a local choice.
        @test ALPHA_FIX_FLOOR == quantile(vec(ALPHA_FIX_Y), P13_ALPHA_BG_QUANTILE)
    end

    @testset "verify_alpha_invariants reports all green on the fixture" begin
        for f in (:bitwise_alpha0, :no_new_zeros, :intensity_conserved, :mask_invariant,
                  :mask_fraction_ok, :max_value_ok, :mbar_monotone, :missing_nonincreasing)
            @test getfield(ALPHA_FIX_INV, f) === true
            @test getfield(ALPHA_GRAY_INV, f) === true
        end
        # The measured residuals are reported, not merely flagged.
        @test isfinite(ALPHA_FIX_INV.max_rel_intensity_err)
        @test isfinite(ALPHA_FIX_INV.realized_max_value)
        @test ALPHA_FIX_INV.source_zero_count isa Integer
        # The field is `no_new_zeros`; an all-positive field would be invalid on both substrates.
        @test hasproperty(ALPHA_FIX_INV, :no_new_zeros)

        # An out-of-band mask is REPORTED rather than thrown, so a caller learns WHICH invariant
        # failed instead of getting an exception from inside the transform.
        bad = verify_alpha_invariants([zeros(16, 16), ones(16, 16)]; max_value = Inf)
        @test bad.mask_fraction_ok === false
        @test bad.realized_mask_fraction == 0.0
    end

    @testset "the sealed holdout is untouched" begin
        # The one real segregated anchor is split = sealed_holdout, sha256 = PENDING-FETCH,
        # bytes = 0, reserved for the Phase-16 blind evaluation. Consuming it here would burn
        # Phase 16 on the very hypothesis Phase 16 evaluates, so no executable line may reach it.
        @test !occursin("open_sealed_holdout", ALPHA_CODE)
        @test !occursin("corpus/", ALPHA_CODE)
        @test !occursin("corpus\\", ALPHA_CODE)
        # D-15 substrate-agnosticism: ONE code path, no provenance branch. A branch would make
        # the simulated and real arms different experiments and their rungs incomparable.
        @test !occursin("load_tiff", ALPHA_CODE)
        @test !occursin("test_images", ALPHA_CODE)
        @test occursin("AbstractMatrix{Float64}", ALPHA_CODE)
        # No special case on a zero alpha: the bitwise property must hold structurally.
        @test !occursin("alpha == 0", ALPHA_CODE)
        @test !occursin("iszero(alpha)", ALPHA_CODE)
    end

    @testset "P13 alpha ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end
end
