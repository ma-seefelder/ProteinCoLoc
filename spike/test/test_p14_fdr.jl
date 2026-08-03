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

# spike/test/test_p14_fdr.jl --- the Bayesian-FDR prefix rule (SC1-a).
#
# THE ONE THING THIS FILE EXISTS TO SEPARATE. The controlled quantity is `E[FDP | data]`, the
# posterior EXPECTED false discovery proportion, and an expectation is a MEAN. A running MAXIMUM
# over the same sorted vector is a different and much more conservative rule -- it controls the
# per-comparison posterior error probability instead -- and on most vectors the two agree, which
# is exactly why a badly chosen fixture proves nothing. The fixture below is chosen BECAUSE it
# separates them:
#
#     v = [0.01, 0.02, 0.15], alpha = 0.10
#       running MEAN: 0.010, 0.015, 0.060   -> all three at or below 0.10 -> k* = 3
#       running max:  0.010, 0.020, 0.150   -> the third exceeds 0.10     -> k* = 2
#
# and the NON-discriminating fixture is recorded beside it, so a later reader can see why the
# first one was chosen rather than having to rediscover it:
#
#     v = [0.01, 0.02, 0.30], alpha = 0.10
#       running MEAN: 0.010, 0.015, 0.110   -> k* = 2
#       running max:  0.010, 0.020, 0.300   -> k* = 2      (the two rules AGREE; useless as a test)
#
# The running-maximum rule is written out ONCE in this file, as `_fdr_running_max_k`, purely as
# the reference the discrimination is asserted against. It is DELIBERATELY absent from the unit
# under test, and the plan's acceptance check greps `spike/p14/fdr.jl` to keep it that way.
#
# THE PREFIX PROPERTY IS ASSERTED, NOT ASSUMED. The running mean of an ASCENDING sequence is
# non-decreasing, which is what makes the accepted set a genuine prefix and the rule equivalent to
# a single threshold on v. If it were not, `findlast` would be selecting a k whose predecessors
# are not all admissible and the "rule" would silently be a subset search.
#
# Mirrors spike/test/test_p13_result.jl: license header -> using Test -> guarded include of the
# unit under test -> fixtures computed ONCE at module level -> one outer @testset. No training, no
# simulation, no figure call; the only random draws ride the FIXTURE stream.
#
# RUN IT PER FILE. The aggregate spike suite exits 1 early at the Phase-4 speedup gate and masks
# every later include block, so a green aggregate run would be no evidence at all:
#
#     julia --project=spike spike/test/test_p14_fdr.jl

using Test

# Unit under test (pulls spike/p14/consts.jl and spike/p14/posterior.jl transitively).
isdefined(@__MODULE__, :p14_bayes_fdr) || include(joinpath(@__DIR__, "..", "p14", "fdr.jl"))

# --- Fixtures built ONCE at module level ---------------------------------------------------

# THE DISCRIMINATING FIXTURE. See the banner for both rules' arithmetic.
const FDR_FIX_DISCRIMINATING = [0.01, 0.02, 0.15]
# The non-discriminating one, recorded so the choice above is auditable rather than arbitrary.
const FDR_FIX_AGREEING = [0.01, 0.02, 0.30]
# A FIXTURE level. It is neither of the phase's two pre-registered alphas: it is written as a
# literal here precisely so this file's expectations cannot move when a pre-registered constant
# does, and so that neither alpha is ever derived from the other (D-07).
const FDR_FIX_ALPHA = 0.10

# The prefix-property sweep. 10^4 vectors is enough that a rule which is a prefix only "usually"
# would be caught, and it costs milliseconds because the whole rule is a sort and a cumulative sum.
const FDR_FIX_N_RANDOM = 10_000
const FDR_FIX_MAX_LEN  = 20

"""
The RUNNING-MAXIMUM rule's k*, present ONLY as the reference the discriminating fixture is
asserted against. This is what `p14_bayes_fdr` must NOT be: it controls the per-comparison
posterior error probability rather than the posterior expected false discovery proportion.
"""
function _fdr_running_max_k(v::AbstractVector{<:Real}, alpha::Real)
    isempty(v) && return 0
    sv = sort(v)
    k  = findlast(<=(alpha), accumulate(max, sv))
    return k === nothing ? 0 : k
end

@testset "P14 Bayesian-FDR prefix rule (SC1-a)" verbose = true begin

    @testset "the running MEAN, not the running max" begin
        r = p14_bayes_fdr(FDR_FIX_DISCRIMINATING, FDR_FIX_ALPHA)
        @test r.k_star == 3
        @test isapprox(r.fdp, 0.06; atol = 1e-12)
        @test r.accepted == [1, 2, 3]
        @test r.t_star == 0.15
        @test isapprox(r.cost_ratio, 0.15 / 0.85; atol = 1e-12)
        @test r.n == 3
        # THE DISCRIMINATION, ASSERTED RATHER THAN ASSERTED-IN-A-COMMENT: the two rules really do
        # differ on this fixture, so a silent substitution of one for the other cannot pass.
        @test _fdr_running_max_k(FDR_FIX_DISCRIMINATING, FDR_FIX_ALPHA) == 2
        @test r.k_star != _fdr_running_max_k(FDR_FIX_DISCRIMINATING, FDR_FIX_ALPHA)

        # And the fixture that does NOT discriminate, recorded so the choice above is auditable.
        a = p14_bayes_fdr(FDR_FIX_AGREEING, FDR_FIX_ALPHA)
        @test a.k_star == 2
        @test _fdr_running_max_k(FDR_FIX_AGREEING, FDR_FIX_ALPHA) == 2
        @test isapprox(a.fdp, 0.015; atol = 1e-12)
    end

    @testset "the accepted set is a genuine PREFIX (asserted, not assumed)" begin
        rng = p14_fix_rng(P14_FIXTURE_COUNTER)
        bad = String[]
        for t in 1:FDR_FIX_N_RANDOM
            n  = rand(rng, 2:FDR_FIX_MAX_LEN)
            u  = rand(rng, n)                       # UNSORTED on purpose (see the next testset)
            sv = sort(u)
            run_mean = cumsum(sv) ./ (1:n)
            issorted(run_mean) ||
                push!(bad, "draw $t: the running mean of an ascending vector is not sorted")
            r = p14_bayes_fdr(u, FDR_FIX_ALPHA)
            expect_k = something(findlast(<=(FDR_FIX_ALPHA), run_mean), 0)
            r.k_star == expect_k ||
                push!(bad, "draw $t: k* = $(r.k_star), expected $expect_k")
            # The accepted items ARE the k* smallest values -- compared by VALUE so ties in `u`
            # cannot make a correct answer look wrong.
            sort(u[r.accepted]) == sv[1:r.k_star] ||
                push!(bad, "draw $t: the accepted set is not the k* smallest items")
            if r.k_star > 0
                (r.fdp <= FDR_FIX_ALPHA) ||
                    push!(bad, "draw $t: the realized fdp $(r.fdp) exceeds alpha")
                (r.t_star == sv[r.k_star]) ||
                    push!(bad, "draw $t: t* is not the k*-th order statistic")
            end
        end
        @test bad == String[]
    end

    @testset "accepted holds ORIGINAL indices, not sorted positions" begin
        r = p14_bayes_fdr([0.9, 0.01], FDR_FIX_ALPHA)
        @test r.accepted == [2]
        @test r.k_star == 1
        @test r.t_star == 0.01
        # The same two values in the other order, so the assertion cannot pass by coincidence.
        s = p14_bayes_fdr([0.01, 0.9], FDR_FIX_ALPHA)
        @test s.accepted == [1]
    end

    @testset "reject nothing, and the degenerate input" begin
        z = p14_bayes_fdr([0.5, 0.6], FDR_FIX_ALPHA)
        @test z.k_star == 0
        @test z.accepted == Int[]
        @test z.accepted isa Vector{Int}
        @test isnan(z.t_star)
        @test isnan(z.cost_ratio)
        @test z.fdp == 0.0
        @test z.n == 2
        # An empty batch is answered BEFORE the arithmetic, not by an empty `cumsum`.
        e = p14_bayes_fdr(Float64[], FDR_FIX_ALPHA)
        @test e.k_star == 0
        @test e.accepted == Int[]
        @test e.n == 0
        @test isnan(e.t_star)
    end

    @testset "validation names the function and the observed value" begin
        for bad_alpha in (0.0, 1.0, -0.1, 1.5)
            @test_throws ArgumentError p14_bayes_fdr(FDR_FIX_DISCRIMINATING, bad_alpha)
        end
        err = try
            p14_bayes_fdr(FDR_FIX_DISCRIMINATING, 1.5)
        catch e
            e
        end
        @test occursin("p14_bayes_fdr", err.msg)
        @test occursin("1.5", err.msg)
        # `v` is a vector of POSTERIOR PROBABILITIES. A NaN would sort to the end and turn every
        # running mean after it into a NaN, silently unordering the very sequence the prefix
        # property is asserted on; a value outside [0,1] is not a probability at all.
        @test_throws ArgumentError p14_bayes_fdr([0.01, NaN], FDR_FIX_ALPHA)
        @test_throws ArgumentError p14_bayes_fdr([0.01, 1.5], FDR_FIX_ALPHA)
        @test_throws ArgumentError p14_bayes_fdr([-0.01, 0.5], FDR_FIX_ALPHA)
    end

    @testset "the decided-subset wrapper cannot omit its denominator" begin
        # STRUCTURAL, NOT DOCUMENTARY. `(alpha, n_decided/n_total)` is ONE quantity: a rule that
        # abstains on 95 % of a batch hits any alpha trivially, so an FDR number without its
        # decided fraction is not a claim. The keyword is REQUIRED for the same reason
        # `p13_calibration_meta` makes `auc` required -- a calibration verdict cannot ship
        # without its discrimination number (spike/p13/result.jl:313-315).
        @test_throws UndefKeywordError p14_fdr_over_decided(FDR_FIX_DISCRIMINATING, FDR_FIX_ALPHA)

        w = p14_fdr_over_decided(FDR_FIX_DISCRIMINATING, FDR_FIX_ALPHA; n_total = 12)
        @test w.fdr_scope === :decided_subset_only
        @test w.n_total == 12
        @test w.n_decided == 3
        @test isapprox(w.decided_fraction, 3 / 12; atol = 1e-12)
        # The wrapper does not change the rule; it only refuses to report it naked.
        @test w.k_star == 3
        @test isapprox(w.fdp, 0.06; atol = 1e-12)
        @test w.accepted == [1, 2, 3]
        # A decided subset larger than the batch it came from is a bookkeeping error, not a
        # coverage above 1.
        @test_throws ArgumentError p14_fdr_over_decided(FDR_FIX_DISCRIMINATING, FDR_FIX_ALPHA;
                                                        n_total = 2)
        @test_throws ArgumentError p14_fdr_over_decided(FDR_FIX_DISCRIMINATING, FDR_FIX_ALPHA;
                                                        n_total = 0)
        # An all-abstain batch is a legal, and very informative, answer.
        allout = p14_fdr_over_decided(Float64[], FDR_FIX_ALPHA; n_total = 7)
        @test allout.n_decided == 0
        @test allout.decided_fraction == 0.0
        @test allout.fdr_scope === :decided_subset_only
    end

    @testset "the rule composes with the composite null it is fed" begin
        # The end-to-end shape: null posteriors come from `p14_null_posterior`, never from a
        # hand-built vector, so the two primitives are asserted to fit together at least once.
        prior = (coloc = 0.2, random = 0.5, exclusion = 0.3)
        ps = [p14_class_posterior((coloc = c, random = 0.0, exclusion = -1.0), prior)
              for c in (6.0, 4.0, 1.0, -2.0)]
        v = [p14_null_posterior(p) for p in ps]
        @test all(0 .<= v .<= 1)
        @test issorted(v)                      # more coloc evidence, less null mass
        r = p14_bayes_fdr(v, 0.20)
        @test r.k_star >= 1
        @test r.fdp <= 0.20
        @test r.accepted == collect(1:r.k_star)
    end

    @testset "P14 FDR ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
