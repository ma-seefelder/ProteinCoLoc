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

# spike/test/test_p13_labels.jl --- the Phase-13 label-surface gate (D-05, D-07, D-11).
#
# WHAT THIS FILE IS FOR. spike/p13/labels.jl holds the two places Phase 13 is most likely to
# be SILENTLY wrong, and in both cases the wrong version looks exactly like the right one:
#
#   D-05 -- labelling on the sample-minus-control contrast alone would publish a strongly
#           COLOCALIZED sample sitting under a MORE colocalized control as "mutually
#           exclusive". The counter-example is one line of arithmetic and zero lines of
#           training, so it is asserted here rather than discovered in a confusion matrix.
#   D-11 -- one-vs-rest head training would make the coloc logit a coloc-vs-MIXTURE ratio
#           instead of log BF(coloc : random), silently breaking D-08's published semantics
#           while every shape, every type and every test that only checks shapes still passes.
#   D-07 -- the per-head correction must be counted over the head's OWN class pair. Copying
#           the binary form gives the wrong denominator, and because the binary term is ~0.01
#           nat the mistake is invisible in the code it was copied from.
#
# It also proves, executably, that labelled three-way data CANNOT be produced before the
# Tier-2 tau measurement exists (D-01/D-06).
#
# Mirrors test_sbc.jl / test_bf.jl / test_p13_consts.jl: license header -> using Test ->
# guarded include of the unit under test -> fixtures computed ONCE at module level -> one
# outer @testset. No training, no image simulation, no figure call -- this gate is under a
# second (the label depends only on theta, which is exactly why it is cheap).

using Test
using Statistics

# Unit under test: the label surface (pulls p13/consts.jl -> simulator/prior.jl).
isdefined(@__MODULE__, :three_way_label) || include(joinpath(@__DIR__, "..", "p13", "labels.jl"))

# THE FIXTURE TAU IS A TEST CONSTANT, NOT THE PRE-REGISTERED P13_TAU. P13_TAU is Tier-2 and
# does not exist yet (D-06 requires it MEASURED, plan 13-10). This suite must stay runnable
# BEFORE that commit, so every call below passes tau EXPLICITLY and nothing here may be read
# as a pre-registered threshold. The value 0.10 is chosen only because 13-RESEARCH E1's
# measured class-mass table quotes that row, which is what testset 8 reproduces.
const TAU_FIX = 0.10

# The FROZEN class-mass table (13-RESEARCH E1, 2e6 i.i.d. draw pairs through
# ghat(Truncated(Cauchy(0,0.3),-1,1)), variant (a)). Quoted, never recomputed-and-blessed.
const LBL_FROZEN_MASSES = (
    (tau = 0.05, exclusion = 0.3850, random = 0.3413, coloc = 0.2738),
    (tau = 0.10, exclusion = 0.3482, random = 0.4108, coloc = 0.2410),
    (tau = 0.15, exclusion = 0.3119, random = 0.4713, coloc = 0.2168),
)

# n = 20_000 gives a binomial standard error of ~0.0034 at a mass of ~0.35, so 0.02 is ~6 SE:
# MONTE-CARLO SLOP AT THIS n, NOT A LOOSENED THRESHOLD. (The frozen numbers come from 2e6
# draws; reproducing them to 2e6 precision in a per-task gate would cost minutes.)
const LBL_FIX_N   = 20_000
const LBL_FIX_TOL = 0.02

# --- Fixtures computed ONCE (the test_bf.jl:49-65 idiom) ------------------------------------
# ALL fixtures ride p13_fix_rng(P13_FIXTURE_COUNTER), so the quick gate NEVER consumes (and so
# never pre-observes) a REPORTED Phase-13 stream (D-01, T-13-01). P13_FIX_SEED != P13_DEV_SEED
# is asserted inside the suite as the explicit statement of that rule.
#
# Each call re-derives the SAME fresh deterministic stream, so the three rows below are the
# SAME 20_000 rho pairs RELABELLED at three taus -- which is what makes the whole frozen table
# reproducible for the cost of one draw set.
const LBL_FIX_MASSES = [
    class_masses(prior_class_draws(LBL_FIX_N; rng = p13_fix_rng(P13_FIXTURE_COUNTER),
                                   tau = row.tau))
    for row in LBL_FROZEN_MASSES
]

# The rho pairs BEHIND that table, drawn from the same fresh fixture stream in the same order,
# so the Pitfall-3 invariant below is checked on exactly the pairs the masses were counted on.
const LBL_FIX_RHO_PAIRS = let rng = p13_fix_rng(P13_FIXTURE_COUNTER)
    [(sample_prior(rng).ρ_true, sample_prior(rng).ρ_true) for _ in 1:LBL_FIX_N]
end

# The sign-mirror of a class: COLOC <-> EXCLUSION, RANDOM -> RANDOM.
_lbl_mirror(c::ThreeWayClass) =
    c === COLOC ? EXCLUSION : (c === EXCLUSION ? COLOC : RANDOM)

@testset "P13 labels and per-head targets (D-05, D-07, D-11)" verbose = true begin

    @testset "fixtures never consume a reported stream (D-01)" begin
        # The explicit statement of the fixture-stream rule: a gate that drew from the
        # REPORTED stream would pre-observe it, and a pre-observed stream can never govern a
        # reported number.
        @test P13_FIX_SEED != P13_DEV_SEED
        @test P13_FIXTURE_COUNTER ∉ (P13_TAU_COUNTER, P13_DATAGEN_COUNTER, P13_GATE_COUNTER,
                                     P13_ALPHA_COUNTER, P13_CONTINUITY_COUNTER)
    end

    @testset "D-05 counter-example (Pitfall 3)" begin
        # THE TRAP: Δρ = 0.8 - 0.9 < 0, so a naive three-way extension of the shipped binary
        # net's contrast (src/amortized/infer.jl:117-121) would publish a strongly COLOCALIZED
        # sample as MUTUALLY EXCLUSIVE.
        @test three_way_label(0.8, 0.9; tau = TAU_FIX) != EXCLUSION
        @test three_way_label(0.8, 0.9; tau = TAU_FIX) === RANDOM

        # The stronger sweep: a strongly colocalized sample is NEVER exclusion, whatever the
        # control does -- including controls far above it.
        for rho_c in -0.9:0.3:0.9
            @test three_way_label(0.8, rho_c; tau = TAU_FIX) != EXCLUSION
        end
    end

    @testset "level and contrast are both required" begin
        # Level passes (0.5 > tau), contrast fails (0.05 < tau) -> RANDOM.
        @test three_way_label(0.5, 0.45; tau = TAU_FIX) === RANDOM
        # Contrast passes (0.95 > tau), level fails (0.05 < tau) -> RANDOM.
        @test three_way_label(0.05, -0.9; tau = TAU_FIX) === RANDOM
        # Both pass, in each sign.
        @test three_way_label(0.9, 0.1; tau = TAU_FIX) === COLOC
        @test three_way_label(-0.9, -0.1; tau = TAU_FIX) === EXCLUSION
    end

    @testset "dead zone and sign symmetry" begin
        # Anything inside the level dead zone is RANDOM regardless of the control: "random"
        # means indistinguishable from zero AT THIS SUMMARY'S RESOLUTION (D-06).
        for rho_s in (-TAU_FIX, -0.05, 0.0, 0.05, TAU_FIX)
            for rho_c in (-0.99, -0.5, 0.0, 0.5, 0.99)
                @test abs(rho_s) <= TAU_FIX
                @test three_way_label(rho_s, rho_c; tau = TAU_FIX) === RANDOM
            end
        end

        # Sign symmetry: negating BOTH rho values mirrors the class. This is a property of
        # variant (a) that variant (b) does not have, and it is why (-0.8, -0.9) is RANDOM --
        # the exact mirror of the D-05 counter-example -- and not EXCLUSION.
        for a in -0.9:0.3:0.9, b in -0.9:0.3:0.9
            @test three_way_label(-a, -b; tau = TAU_FIX) ===
                  _lbl_mirror(three_way_label(a, b; tau = TAU_FIX))
        end
        @test three_way_label(-0.8, -0.9; tau = TAU_FIX) === RANDOM
    end

    @testset "tau must be positive" begin
        # A zero or negative dead zone degenerates variant (a) into the REJECTED variant (b)
        # (strict inequality, no dead zone), so it is refused at the boundary.
        @test_throws ArgumentError three_way_label(0.5, 0.0; tau = 0.0)
        @test_throws ArgumentError three_way_label(0.5, 0.0; tau = -0.1)
    end

    @testset "D-11 participation weights" begin
        # The literal 4-row targets, legend [y_C; w_C; y_E; w_E].
        @test head_targets(COLOC)     == Float32[1, 1, 0, 0]
        @test head_targets(EXCLUSION) == Float32[0, 0, 1, 1]
        @test head_targets(RANDOM)    == Float32[0, 1, 0, 1]

        # RANDOM is the shared NEGATIVE class of BOTH heads: both participate, both target 0.
        # That shared reference is what makes the two logits ratios against the same
        # denominator, i.e. what makes D-08's two published numbers well defined.
        let t = head_targets(RANDOM)
            @test t[2] == 1f0 && t[4] == 1f0      # w_C, w_E both active
            @test t[1] == 0f0 && t[3] == 0f0      # y_C, y_E both negative
        end
        # Each class switches the OTHER head off, which is the per-head class restriction.
        @test head_targets(COLOC)[4]     == 0f0   # w_E off on a coloc sample
        @test head_targets(EXCLUSION)[2] == 0f0   # w_C off on an exclusion sample

        # The batch contract the masked loss indexes positionally.
        @test size(target_matrix([COLOC, RANDOM, EXCLUSION])) == (4, 3)
        @test target_matrix([COLOC, RANDOM, EXCLUSION])[:, 2] == head_targets(RANDOM)
        @test eltype(target_matrix([COLOC, RANDOM, EXCLUSION])) === Float32
        @test_throws ArgumentError target_matrix(ThreeWayClass[])
    end

    @testset "D-07 head log-odds counts only its own pair" begin
        labels = vcat(fill(COLOC, 30), fill(RANDOM, 20), fill(EXCLUSION, 50))
        h = head_log_odds(labels)
        @test isapprox(h.coloc,     log(30 / 20); atol = 1e-12)
        @test isapprox(h.exclusion, log(50 / 20); atol = 1e-12)

        # THE SILENT D-08 BREAK, ASSERTED AGAINST. log(30/70) is the one-vs-rest denominator
        # (the coloc COMPLEMENT: random + exclusion). A head trained that way would emit a
        # coloc-vs-MIXTURE ratio, and correcting it with the mixture log-odds would make the
        # published number "log BF(coloc : random)" a name for a different quantity.
        @test !isapprox(h.coloc, log(30 / 70); atol = 1e-6)
        # ...and the correction is nowhere near 0, so it cannot be waved away as the binary
        # term (~0.01 nat) can.
        @test abs(h.coloc) > 0.3

        # The keyword form agrees with the NamedTuple surface, and each head's reference is
        # RANDOM -- never the other positive class.
        @test measure_head_log_odds(labels; positive = COLOC,     negative = RANDOM) == h.coloc
        @test measure_head_log_odds(labels; positive = EXCLUSION, negative = RANDOM) == h.exclusion
        # A head cannot be its own reference.
        @test_throws ArgumentError measure_head_log_odds(labels; positive = COLOC, negative = COLOC)
    end

    @testset "D-07 degenerate balance errors" begin
        # No RANDOM members: BOTH heads lose their reference class, so both must fail loudly
        # rather than return log(n/0) = Inf and poison every corrected log BF downstream.
        no_random = vcat(fill(COLOC, 10), fill(EXCLUSION, 10))
        @test_throws ErrorException measure_head_log_odds(no_random; positive = COLOC,
                                                          negative = RANDOM)
        @test_throws ErrorException measure_head_log_odds(no_random; positive = EXCLUSION,
                                                          negative = RANDOM)
        @test_throws ErrorException head_log_odds(no_random)
        # A missing POSITIVE class is equally degenerate (log(0/n) = -Inf).
        @test_throws ErrorException measure_head_log_odds(fill(RANDOM, 10); positive = COLOC,
                                                          negative = RANDOM)
        @test_throws ArgumentError class_masses(ThreeWayClass[])
    end

    @testset "class masses reproduce the frozen table" begin
        # The frozen row values are the 2e6-draw numbers of 13-RESEARCH E1. The tolerance is
        # the MONTE-CARLO SLOP at n = 20_000 (~6 binomial SE), not a loosened threshold: the
        # point is that the label rule in this repo IS the rule those masses were measured
        # under, so a silent change of variant, sign or dead zone shows up here.
        for (row, got) in zip(LBL_FROZEN_MASSES, LBL_FIX_MASSES)
            @test isapprox(got.exclusion, row.exclusion; atol = LBL_FIX_TOL)
            @test isapprox(got.random,    row.random;    atol = LBL_FIX_TOL)
            @test isapprox(got.coloc,     row.coloc;     atol = LBL_FIX_TOL)
            @test isapprox(got.exclusion + got.random + got.coloc, 1.0; atol = 1e-12)
        end

        # The measured per-head corrections at the fixture tau land near the pi-level values
        # 13-RESEARCH E1 quotes (-0.5333 coloc, -0.1654 exclusion) -- 15-50x the shipped
        # binary net's -0.0102, which is the whole reason D-07 exists.
        let h = head_log_odds(prior_class_draws(LBL_FIX_N;
                                                rng = p13_fix_rng(P13_FIXTURE_COUNTER),
                                                tau = TAU_FIX))
            @test isapprox(h.coloc,     -0.5333; atol = 0.05)
            @test isapprox(h.exclusion, -0.1654; atol = 0.05)
        end

        # PITFALL 3's WARNING SIGN, re-derived from the stored rho pairs: NO exclusion-class
        # member may come from a positive rho_sample. This is the invariant that would break
        # first if the cut ever drifted back onto the contrast alone. Violations are COUNTED
        # (rather than asserted per pair) so one failure reports how widespread the drift is
        # instead of emitting 20_000 near-identical records.
        let bad_e = 0, bad_c = 0, n_e = 0, n_c = 0
            for (rho_s, rho_c) in LBL_FIX_RHO_PAIRS
                c = three_way_label(rho_s, rho_c; tau = TAU_FIX)
                if c === EXCLUSION
                    n_e += 1
                    rho_s < -TAU_FIX || (bad_e += 1)
                elseif c === COLOC
                    n_c += 1
                    rho_s > TAU_FIX || (bad_c += 1)
                end
            end
            @test n_e > 0 && n_c > 0        # the invariant is not vacuous
            @test bad_e == 0                # no EXCLUSION member with rho_sample >= -tau
            @test bad_c == 0                # no COLOC member with rho_sample <= +tau
        end
    end

    @testset "tau is measured, and the label default resolves to it (Tier 2)" begin
        # FLIPPED BY PLAN 13-10, exactly as the previous version of this testset said it would
        # be. Until 13-10 ran, the THROW was the D-06 ordering guarantee: labelled three-way
        # data could not be generated at a tau that had not been measured. That guarantee is
        # now carried by git history instead -- artifact commit 581602a, which contains no
        # Tier-2 value, strictly precedes the Tier-2 append. What must be asserted HERE is the
        # thing that replaced it: the label boundary's default really is the MEASURED tau and
        # not some local constant that happens to be in scope.
        @test isdefined(@__MODULE__, :P13_TAU)
        @test p13_tau() === P13_TAU
        @test three_way_label(0.5, 0.0) === three_way_label(0.5, 0.0; tau = P13_TAU)
        # ...and the default is genuinely load-bearing: the same pair labels differently under
        # a tau wide enough to swallow it, so a silently wrong default could not pass unnoticed.
        @test three_way_label(0.5, 0.0) !== three_way_label(0.5, 0.0; tau = 0.6)
        # The generator no longer throws. It is called on the FIXTURE stream on purpose: its
        # default rng rides P13_DATAGEN_COUNTER, and a gate that drew from that reported
        # sub-stream would pre-observe the very draws the labelled data is generated from (D-01).
        @test length(prior_class_draws(4; rng = p13_fix_rng(P13_FIXTURE_COUNTER))) == 4
        # The frozen variant this file tests is still the frozen variant.
        @test P13_CUT_VARIANT === :tau_contrast
    end

    @testset "P13 labels ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
