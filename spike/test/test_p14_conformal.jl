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

# spike/test/test_p14_conformal.jl --- the hand-rolled split-conformal hedge (SC1-e).
#
# THE ONE SUBSTITUTION THIS FILE EXISTS TO CATCH IS THE INTERPOLATED QUANTILE (Pitfall 6).
# The conformal coverage guarantee is stated for an ORDER STATISTIC:
#
#     qhat = sort(s)[ceil(Int, (n + 1) * (1 - alpha))]
#
# `Statistics.quantile(s, 1 - alpha)` is a DIFFERENT number. Its default type-7 rule interpolates
# LINEARLY between two neighbouring order statistics, so it returns something strictly smaller for
# any non-degenerate sample, and a smaller qhat admits fewer classes into the set, which breaks the
# >= 1 - alpha finite-sample bound by O(1/n). AT n = 2000 THAT IS A FIVE-IN-TEN-THOUSAND SHIFT IN
# THE REPORTED COVERAGE -- entirely invisible next to the binomial noise of the evaluation set, and
# therefore undetectable without a DIRECTED test. This file is that directed test, and it works by
# choosing n so small that the two rules cannot be confused:
#
#     n = 9, alpha = 0.10  ->  k = ceil(10 * 0.90) = 9 = n   ->  qhat == maximum(s)
#                              Statistics.quantile(s, 0.9) interpolates at index 1 + 8*0.9 = 8.2,
#                              i.e. strictly between s[8] and s[9]  ->  strictly SMALLER than max
#
# BOTH assertions are required and neither is redundant: `== maximum(s)` alone would still pass a
# unit that happened to return the max by accident, and `!= Statistics.quantile(...)` alone would
# pass any wrong-but-different number. Either equality flipping is the failure (14-VALIDATION
# SC1-e).
#
# THE SECOND THING IT EXISTS TO CATCH is the EMPTY conformal set being folded into :ambiguous. An
# empty set says the observed point is more nonconforming than (1 - alpha) of everything the
# calibration set contained -- a distribution-free misspecification signal on a completely
# different channel from the Mahalanobis/PP/noise OOD detector, and the only signal in this phase
# that does not need the simulator to be right about densities, only about exchangeability.
# Collapsing it into "ambiguous" would throw that away while still abstaining, so the loss would be
# silent.
#
# NO DEPENDENCY WAS ADDED (D-02). SC1's named library is DROPPED and conformal is met IN SUBSTANCE:
# the whole construction under test is `sort` + `ceil` + integer indexing. `Statistics` is imported
# HERE, in the test, for exactly one purpose -- to construct the WRONG answer that the
# discrimination is asserted against.
#
# Mirrors spike/test/test_p14_fdr.jl: license header -> using Test -> guarded include of the unit
# under test -> fixtures computed ONCE at module level -> one outer @testset. No training, no
# simulation, no figure call; the only random draws ride the FIXTURE stream.
#
# RUN IT PER FILE. The aggregate spike suite exits 1 early at the Phase-4 speedup gate and masks
# every later include block, so a green aggregate run would be no evidence at all:
#
#     julia --project=spike spike/test/test_p14_conformal.jl

using Test
import Statistics       # ONLY to build the interpolated answer the order statistic is asserted
                        # AGAINST. Qualified rather than `using`, so nothing in this file can call
                        # the interpolating quantile by accident.

# Unit under test (pulls spike/p14/consts.jl and spike/p14/posterior.jl transitively).
isdefined(@__MODULE__, :p14_conformal_quantile) ||
    include(joinpath(@__DIR__, "..", "p14", "conformal.jl"))

# --- Fixtures built ONCE at module level ---------------------------------------------------

# A FIXTURE level. It is deliberately a literal here rather than a read of the pre-registered
# constant, so this file's expectations cannot move when a pre-registered value does, and so that
# neither of the phase's two alphas is ever derived from the other (D-07).
const CONF_FIX_ALPHA = 0.10

# n = 9 is chosen BECAUSE ceil(10 * 0.9) == 9 == n, which pins qhat to the maximum exactly.
const CONF_FIX_N_SMALL = 9
const CONF_FIX_SCORES  = rand(p14_fix_rng(P14_FIXTURE_COUNTER), CONF_FIX_N_SMALL)

# A larger, INTERIOR case: k lands strictly inside the sorted vector, so an off-by-one in the
# index cannot hide behind "the answer is the maximum anyway".
const CONF_FIX_N_MID    = 19
const CONF_FIX_SCORES_MID = rand(p14_fix_rng(P14_FIXTURE_COUNTER), CONF_FIX_N_MID)

# The infeasibility case. ceil(31 * 0.97) = 31 > 30, so the call must THROW rather than clamp --
# the deliberate difference from `ghat`, which is the closest in-repo analog and DOES clamp.
const CONF_FIX_N_INFEASIBLE     = 30
const CONF_FIX_ALPHA_INFEASIBLE = 0.03
const CONF_FIX_SCORES_30        = rand(p14_fix_rng(P14_FIXTURE_COUNTER), CONF_FIX_N_INFEASIBLE)

# The three set-size regimes, hand-built so each one is reached for a stated arithmetic reason.
const CONF_FIX_P_SHARP   = (coloc = 0.95, random = 0.03, exclusion = 0.02)  # -> :singleton
const CONF_FIX_P_SPLIT   = (coloc = 0.50, random = 0.45, exclusion = 0.05)  # -> :ambiguous
const CONF_FIX_P_DIFFUSE = (coloc = 0.40, random = 0.35, exclusion = 0.25)  # -> :empty
const CONF_FIX_QHAT_TIGHT = 0.10   # scores 0.05 / 0.97 / 0.98  -> only coloc admitted
const CONF_FIX_QHAT_LOOSE = 0.60   # scores 0.50 / 0.55 / 0.95  -> coloc AND random admitted
const CONF_FIX_QHAT_EMPTY = 0.50   # scores 0.60 / 0.65 / 0.75  -> nothing admitted; note 0.5 < 2/3

# THE SAME POSTERIOR, DECLARED IN A DIFFERENT KEY ORDER. Five class orderings coexist in this
# codebase and every headline metric is invariant to a CONSISTENT relabelling, so a positional read
# produces numbers that look right. This fixture is what makes "read by name" executable.
const CONF_FIX_P_SHUFFLED = (exclusion = 0.02, coloc = 0.95, random = 0.03)

# A labelled calibration batch for the end-to-end `p14_conformal_calibrate` shape.
const CONF_FIX_CAL_POSTERIORS = [(coloc = 0.90, random = 0.07, exclusion = 0.03),
                                 (coloc = 0.10, random = 0.85, exclusion = 0.05),
                                 (coloc = 0.05, random = 0.15, exclusion = 0.80),
                                 (coloc = 0.60, random = 0.30, exclusion = 0.10),
                                 (coloc = 0.20, random = 0.70, exclusion = 0.10),
                                 (coloc = 0.45, random = 0.45, exclusion = 0.10),
                                 (coloc = 0.80, random = 0.15, exclusion = 0.05),
                                 (coloc = 0.30, random = 0.25, exclusion = 0.45),
                                 (coloc = 0.15, random = 0.75, exclusion = 0.10)]
const CONF_FIX_CAL_CLASSES = [:coloc, :random, :exclusion, :coloc, :random,
                              :coloc, :coloc, :exclusion, :random]

@testset "P14 hand-rolled split conformal (SC1-e, D-02)" verbose = true begin

    @testset "the order statistic, NOT the interpolated quantile (Pitfall 6)" begin
        # The fixture has to be non-degenerate or the two rules coincide and prove nothing.
        @test length(unique(CONF_FIX_SCORES)) == CONF_FIX_N_SMALL

        q = p14_conformal_quantile(CONF_FIX_SCORES, CONF_FIX_ALPHA)
        # (a) the index arithmetic, spelled out so a failure says WHICH half broke
        @test ceil(Int, (CONF_FIX_N_SMALL + 1) * (1 - CONF_FIX_ALPHA)) == CONF_FIX_N_SMALL
        # (b) BOTH assertions are required; either equality flipping is the failure (SC1-e)
        @test q == maximum(CONF_FIX_SCORES)
        @test q != Statistics.quantile(CONF_FIX_SCORES, 1 - CONF_FIX_ALPHA)
        # and the DIRECTION matters too: interpolation returns something strictly SMALLER, which is
        # what silently breaks the >= 1 - alpha bound rather than merely disagreeing with it.
        @test Statistics.quantile(CONF_FIX_SCORES, 1 - CONF_FIX_ALPHA) < q

        # (c) THE INTERIOR CASE. k = ceil(20 * 0.9) = 18 sits strictly inside a 19-long vector, so
        #     an off-by-one cannot hide behind "the answer is the maximum anyway".
        k_mid = ceil(Int, (CONF_FIX_N_MID + 1) * (1 - CONF_FIX_ALPHA))
        @test k_mid == 18
        @test k_mid < CONF_FIX_N_MID
        qm = p14_conformal_quantile(CONF_FIX_SCORES_MID, CONF_FIX_ALPHA)
        @test qm == sort(CONF_FIX_SCORES_MID)[k_mid]
        @test qm != Statistics.quantile(CONF_FIX_SCORES_MID, 1 - CONF_FIX_ALPHA)
        @test qm isa Float64

        # (d) the input is not mutated -- `sort`, never `sort!`, on a caller's vector
        before = copy(CONF_FIX_SCORES_MID)
        p14_conformal_quantile(CONF_FIX_SCORES_MID, CONF_FIX_ALPHA)
        @test CONF_FIX_SCORES_MID == before
    end

    @testset "an infeasible alpha THROWS (it does not clamp)" begin
        @test_throws ArgumentError p14_conformal_quantile(CONF_FIX_SCORES_30,
                                                          CONF_FIX_ALPHA_INFEASIBLE)
        err = try
            p14_conformal_quantile(CONF_FIX_SCORES_30, CONF_FIX_ALPHA_INFEASIBLE)
        catch e
            e
        end
        # The message names the FUNCTION, the observed alpha, the observed n, and the GENERAL
        # feasibility bound -- computed from the argument n, never a hard-coded figure (D-03a).
        @test occursin("p14_conformal_quantile", err.msg)
        @test occursin(string(CONF_FIX_ALPHA_INFEASIBLE), err.msg)
        @test occursin(string(CONF_FIX_N_INFEASIBLE), err.msg)
        @test occursin("1/(n + 1)", err.msg)

        # Just inside the bound it does NOT throw. The bound is written as an expression in n, so
        # this line stays correct if the fixture size ever changes.
        feasible = 1 / (CONF_FIX_N_INFEASIBLE + 1) + 1e-9
        @test p14_conformal_quantile(CONF_FIX_SCORES_30, feasible) isa Float64
        @test p14_conformal_quantile(CONF_FIX_SCORES_30, feasible) == maximum(CONF_FIX_SCORES_30)

        # The other argument failures, each naming the function.
        @test_throws ArgumentError p14_conformal_quantile(Float64[], CONF_FIX_ALPHA)
        @test_throws ArgumentError p14_conformal_quantile(CONF_FIX_SCORES, 0.0)
        @test_throws ArgumentError p14_conformal_quantile(CONF_FIX_SCORES, 1.0)
        @test_throws ArgumentError p14_conformal_quantile([0.1, NaN, 0.3], CONF_FIX_ALPHA)
        @test_throws ArgumentError p14_conformal_quantile([0.1, Inf, 0.3], CONF_FIX_ALPHA)
    end

    @testset "the three set-size regimes, and :empty is NOT :ambiguous" begin
        s1 = p14_conformal_set(CONF_FIX_P_SHARP, CONF_FIX_QHAT_TIGHT)
        @test s1.status === :singleton
        @test s1.set == (:coloc,)
        @test s1.qhat == CONF_FIX_QHAT_TIGHT

        s2 = p14_conformal_set(CONF_FIX_P_SPLIT, CONF_FIX_QHAT_LOOSE)
        @test s2.status === :ambiguous
        @test s2.set == (:coloc, :random)
        @test length(s2.set) >= 2

        # THE EMPTY CASE, constructed so max(p) < 1 - qhat -- the arithmetic asserted, not assumed.
        pmax = max(CONF_FIX_P_DIFFUSE.coloc, CONF_FIX_P_DIFFUSE.random, CONF_FIX_P_DIFFUSE.exclusion)
        @test pmax < 1 - CONF_FIX_QHAT_EMPTY
        @test CONF_FIX_QHAT_EMPTY < 2 / 3        # with 3 classes an empty set needs qhat < 2/3
        s3 = p14_conformal_set(CONF_FIX_P_DIFFUSE, CONF_FIX_QHAT_EMPTY)
        @test s3.status === :empty
        @test s3.set == ()
        @test isempty(s3.set)
        # THE WHOLE POINT: the empty set has its OWN symbol and is not reported as ambiguous.
        @test s3.status !== :ambiguous

        # The vocabulary is closed and its three members are exactly the three regimes above.
        @test P14_CONFORMAL_STATUSES === (:singleton, :ambiguous, :empty)
        @test s1.status in P14_CONFORMAL_STATUSES
        @test s2.status in P14_CONFORMAL_STATUSES
        @test s3.status in P14_CONFORMAL_STATUSES

        @test_throws ArgumentError p14_conformal_set((coloc = 0.5, random = 0.5), CONF_FIX_QHAT_TIGHT)
        @test_throws ArgumentError p14_conformal_set(CONF_FIX_P_SHARP, NaN)
    end

    @testset "scores are read BY NAME, never positionally" begin
        @test p14_lac_score(CONF_FIX_P_SHARP, :coloc) == 1.0 - 0.95
        @test p14_lac_score(CONF_FIX_P_SHARP, :random) == 1.0 - 0.03
        @test p14_lac_score(CONF_FIX_P_SHARP, :exclusion) == 1.0 - 0.02

        # The same posterior declared in a DIFFERENT key order gives the same three answers.
        for y in P14_CLASS_KEYS
            @test p14_lac_score(CONF_FIX_P_SHUFFLED, y) == p14_lac_score(CONF_FIX_P_SHARP, y)
        end
        # ... and the SET is emitted in the frozen P14_CLASS_KEYS order, not the declaration order.
        @test p14_conformal_set(CONF_FIX_P_SHUFFLED, CONF_FIX_QHAT_TIGHT).set == (:coloc,)

        @test_throws ArgumentError p14_lac_score(CONF_FIX_P_SHARP, :colocalisation)
        @test_throws ArgumentError p14_lac_score(CONF_FIX_P_SHARP, :null)
    end

    @testset "the vacuity flag is a LABEL, not a verdict (Pitfall 9)" begin
        # A hedge that never fires looks like a pass and measures nothing. Same design as
        # `vacuous_pass` for ECE (spike/p13/result.jl:284-298): a label beside the number.
        allsingle = [p14_conformal_set(CONF_FIX_P_SHARP, CONF_FIX_QHAT_TIGHT) for _ in 1:50]
        d = p14_hedge_diagnostics(allsingle, CONF_FIX_QHAT_TIGHT)
        @test d.n == 50
        @test d.qhat == CONF_FIX_QHAT_TIGHT
        @test d.singleton_rate == 1.0
        @test d.ambiguous_rate == 0.0
        @test d.empty_rate == 0.0
        @test d.vacuous_hedge === true          # a LABEL...
        @test d isa NamedTuple                  # ...returned beside the rates, never thrown

        # A hedge that DOES fire is not labelled vacuous.
        mixed = vcat(allsingle,
                     [p14_conformal_set(CONF_FIX_P_SPLIT, CONF_FIX_QHAT_LOOSE) for _ in 1:5],
                     [p14_conformal_set(CONF_FIX_P_DIFFUSE, CONF_FIX_QHAT_EMPTY) for _ in 1:5])
        m = p14_hedge_diagnostics(mixed, CONF_FIX_QHAT_TIGHT)
        @test m.n == 60
        @test isapprox(m.ambiguous_rate, 5 / 60; atol = 1e-12)
        @test isapprox(m.empty_rate, 5 / 60; atol = 1e-12)
        @test isapprox(m.singleton_rate, 50 / 60; atol = 1e-12)
        @test isapprox(m.singleton_rate + m.ambiguous_rate + m.empty_rate, 1.0; atol = 1e-12)
        @test m.ambiguous_rate + m.empty_rate >= P14_AMBIGUOUS_RATE_FLOOR
        @test m.vacuous_hedge === false

        # Bare status symbols are accepted as well as full set tuples -- one vocabulary, two shapes.
        @test p14_hedge_diagnostics([:singleton, :ambiguous], CONF_FIX_QHAT_TIGHT).n == 2
        # An EMPTY collection is a bookkeeping error, not a rate: 0 + 0 < the floor would label a
        # hedge that was never exercised as "vacuous", which is the right word for the wrong reason.
        @test_throws ArgumentError p14_hedge_diagnostics([], CONF_FIX_QHAT_TIGHT)
        @test_throws ArgumentError p14_hedge_diagnostics([:singleton, :enormous],
                                                         CONF_FIX_QHAT_TIGHT)
    end

    @testset "calibration is sort + ceil + integer indexing, and nothing else (D-02)" begin
        c = p14_conformal_calibrate(CONF_FIX_CAL_POSTERIORS, CONF_FIX_CAL_CLASSES, CONF_FIX_ALPHA)
        @test c.n_cal == length(CONF_FIX_CAL_POSTERIORS)
        @test c.alpha_conformal == CONF_FIX_ALPHA

        # The qhat is recomputed here from the labelled pairs BY HAND, so the assertion is against
        # the construction rather than against the unit's own answer.
        manual = sort([1.0 - getproperty(p, y)
                       for (p, y) in zip(CONF_FIX_CAL_POSTERIORS, CONF_FIX_CAL_CLASSES)])
        k = ceil(Int, (length(manual) + 1) * (1 - CONF_FIX_ALPHA))
        @test k == length(manual)
        @test c.qhat == manual[k]
        @test c.qhat == maximum(manual)
        @test c.scores_summary.min == minimum(manual)
        @test c.scores_summary.max == maximum(manual)
        @test c.scores_summary.min <= c.scores_summary.median <= c.scores_summary.max

        @test_throws ArgumentError p14_conformal_calibrate(CONF_FIX_CAL_POSTERIORS,
                                                           CONF_FIX_CAL_CLASSES[1:3],
                                                           CONF_FIX_ALPHA)
        @test_throws ArgumentError p14_conformal_calibrate(CONF_FIX_CAL_POSTERIORS,
                                                           CONF_FIX_CAL_CLASSES, 0.0)
    end

    @testset "the hedge scores THE SAME posterior the FDR rule sorts on" begin
        # A conformal set built on a different object than the FDR rule's `v` can contradict it
        # (set = {coloc} while `v` says reject). Asserting they compose at least once is what keeps
        # the two layers coherent by construction rather than by coincidence.
        prior = (coloc = 0.2, random = 0.5, exclusion = 0.3)
        p = p14_class_posterior((coloc = 6.0, random = 0.0, exclusion = -1.0), prior)
        @test isapprox(p14_lac_score(p, :coloc) + p.coloc, 1.0; atol = 1e-12)
        # `p14_confidence` is the max class mass; the LAC score of the arg-max is its complement,
        # so the hedge and the risk-coverage curve order the batch the SAME way.
        @test isapprox(p14_lac_score(p, :coloc), 1.0 - p14_confidence(p); atol = 1e-12)
        @test p14_conformal_set(p, 0.5).status === :singleton
        @test p14_conformal_set(p, 0.5).set == (:coloc,)
        # The composite null the FDR rule sorts on is the same 1 - p.coloc quantity.
        @test isapprox(p14_null_posterior(p), p14_lac_score(p, :coloc); atol = 1e-12)
    end

    @testset "P14 conformal ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
