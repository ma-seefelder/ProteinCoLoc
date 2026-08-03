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

# spike/test/test_p14_posterior.jl --- the three-class posterior gate.
#
# WHAT THIS FILE EXISTS TO CATCH, AND WHY IT IS A TASK RATHER THAN A NOTE.
# `p14_class_posterior` is where Phase 14 introduces the one piece of information Phase 13
# refused to introduce -- the class prior -- and it is therefore where a SILENT CLASS
# PERMUTATION would enter. FIVE class orderings coexist upstream: `three_way_probs` returns
# (coloc, exclusion); `three_way_log_bf` returns (coloc, exclusion, random); `ThreeWayLogBF`
# stores (coloc, random, exclusion); the Phase-13 gate's confusion order is
# ["exclusion","random","coloc"]; and `meta.pi_class_masses` is keyed (exclusion, random,
# coloc). ECE, AUC and FDR are ALL invariant to a *consistent* relabelling, so a positional
# index anywhere would produce numbers that look right and mean something else. That is the
# exact shape of the Phase-12 permutation artifact: an orthogonal permutation preserves norms
# and spectra, so the bug was O(1) and yet invisible to structural verification -- it silently
# inverted a phase's headline finding.
#
# AN EQUAL-MASS FIXTURE CANNOT CATCH A PERMUTATION. That is the whole reason the fixture below
# is DELIBERATELY UNEQUAL (0.7 / 0.2 / 0.1) and the prior it rides is unequal too (0.2 / 0.5 /
# 0.3): under equal masses every permutation is the identity on the numbers, so the test would
# pass against a broken implementation. The flat fixture is kept as a separate sanity check and
# is labelled as such, so no reader can mistake it for the class-order proof.
#
# THE FIXTURE IS BUILT BY INVERTING THE ARITHMETIC, not by running the implementation and
# recording what it said. Target posteriors are CHOSEN, then the two log Bayes factors that
# produce them are SOLVED FOR in closed form. A fixture recorded from the unit under test would
# be a golden of the bug as readily as of the behaviour.
#
# Mirrors spike/test/test_p13_result.jl: license header -> using Test -> guarded include of the
# unit under test -> fixtures computed ONCE at module level -> one outer @testset. No training,
# no simulation, no figure call, and no REPORTED stream is consumed -- the one place this file
# draws random numbers rides the FIXTURE stream, which is what P14_FIX_SEED exists for.
#
# RUN IT PER FILE. The aggregate spike suite exits 1 early at the Phase-4 speedup gate and masks
# every later include block, so a green aggregate run would be no evidence at all:
#
#     julia --project=spike spike/test/test_p14_posterior.jl

using Test

# Unit under test. Pulls spike/p14/consts.jl, spike/p14/provenance.jl and the Phase-13 net,
# result and label surfaces transitively, all through guarded includes.
isdefined(@__MODULE__, :p14_class_posterior) ||
    include(joinpath(@__DIR__, "..", "p14", "posterior.jl"))

# --- Fixtures built ONCE at module level ---------------------------------------------------
#
# THE UNEQUAL PRIOR. Unequal in its own right, so a permutation of the PRIOR (not only of the
# posterior) is caught as well: with a flat prior, reading `prior.exclusion` where
# `prior.coloc` was meant would be a no-op.
const POST_FIX_PRIOR = (coloc = 0.2, random = 0.5, exclusion = 0.3)

# THE UNEQUAL TARGET. Chosen first; the log Bayes factors are solved for below.
const POST_FIX_TARGET = (coloc = 0.7, random = 0.2, exclusion = 0.1)

# THE INVERSION. p(k) is proportional to prior(k) * exp(logbf(k)) with logbf(random) == 0, so
#     logbf(k) = log(p(k)/p(random)) - log(prior(k)/prior(random))
# reproduces the chosen target exactly. The `random` entry is the D-08 structural zero: the
# reference class scored against itself.
const POST_FIX_LOGBF = (
    coloc     = log(POST_FIX_TARGET.coloc / POST_FIX_TARGET.random) -
                log(POST_FIX_PRIOR.coloc / POST_FIX_PRIOR.random),
    random    = 0.0,
    exclusion = log(POST_FIX_TARGET.exclusion / POST_FIX_TARGET.random) -
                log(POST_FIX_PRIOR.exclusion / POST_FIX_PRIOR.random),
)

# THE SAME VALUES IN A DIFFERENT DECLARATION ORDER. If the implementation read by POSITION
# instead of by NAME, these two would disagree -- which is precisely the assertion.
const POST_FIX_PRIOR_SHUFFLED = (exclusion = POST_FIX_PRIOR.exclusion,
                                 coloc     = POST_FIX_PRIOR.coloc,
                                 random    = POST_FIX_PRIOR.random)
const POST_FIX_LOGBF_SHUFFLED = (exclusion = POST_FIX_LOGBF.exclusion,
                                 random    = POST_FIX_LOGBF.random,
                                 coloc     = POST_FIX_LOGBF.coloc)

# The flat sanity fixture. IT IS NOT THE CLASS-ORDER PROOF -- see the banner.
const POST_FIX_FLAT_PRIOR = (coloc = 1/3, random = 1/3, exclusion = 1/3)
const POST_FIX_FLAT_LOGBF = (coloc = 0.0, random = 0.0, exclusion = 0.0)

# The MEASURED class prior, read off the net's own metadata rather than retyped.
const POST_FIX_NET = load_three_way(joinpath(P14_P13_DIR, "three_way_net.jld2"))
const POST_FIX_PI  = POST_FIX_NET.meta.pi_class_masses

# The re-derivation scale. Large enough that three binomial standard errors is a tight band
# (about 0.009 at the coloc mass), small enough that the whole draw costs a fraction of a
# second: the label depends only on theta, so no image is simulated.
const POST_FIX_N_REDERIVE = 20_000

@testset "P14 three-class posterior (D-04, class order, prior re-derivation)" verbose = true begin

    @testset "the flat fixture is a SANITY check, not the class-order proof" begin
        # Equal prior, zero evidence: every class must come back at exactly one third. This
        # catches a normalization slip and nothing else -- under equal masses a permutation is
        # the identity on the numbers, which is why the unequal fixture below exists.
        p = p14_class_posterior(POST_FIX_FLAT_LOGBF, POST_FIX_FLAT_PRIOR)
        @test isapprox(p.coloc,     1/3; atol = 1e-12)
        @test isapprox(p.random,    1/3; atol = 1e-12)
        @test isapprox(p.exclusion, 1/3; atol = 1e-12)
        @test isapprox(p.coloc + p.random + p.exclusion, 1.0; atol = 1e-12)
        # The key ORDER of the returned NamedTuple is the frozen contract, so a downstream
        # reader who does destructure positionally still gets the documented layout.
        @test keys(p) === P14_CLASS_KEYS
    end

    @testset "the class order is not permuted (the Phase-12 artifact shape)" begin
        p = p14_class_posterior(POST_FIX_LOGBF, POST_FIX_PRIOR)
        # THREE SEPARATE ASSERTIONS, so a failure NAMES the field that landed on the wrong
        # class. A single fused `@test p == target` would report only that "something" moved.
        @test isapprox(p.coloc,     0.7; atol = 1e-12)
        @test isapprox(p.random,    0.2; atol = 1e-12)
        @test isapprox(p.exclusion, 0.1; atol = 1e-12)
        @test isapprox(p.coloc + p.random + p.exclusion, 1.0; atol = 1e-12)

        # READ BY NAME, NOT BY POSITION: same values, different declaration order, same answer.
        q = p14_class_posterior(POST_FIX_LOGBF_SHUFFLED, POST_FIX_PRIOR_SHUFFLED)
        @test q.coloc     === p.coloc
        @test q.random    === p.random
        @test q.exclusion === p.exclusion
        @test keys(q) === P14_CLASS_KEYS
    end

    @testset "the keys are validated, never silently defaulted" begin
        @test_throws ArgumentError p14_class_posterior((coloc = 0.0, random = 0.0),
                                                       POST_FIX_PRIOR)
        @test_throws ArgumentError p14_class_posterior(POST_FIX_LOGBF,
                                                       (coloc = 0.5, random = 0.5))
        @test_throws ArgumentError p14_class_posterior(
            (coloc = 0.0, random = 0.0, exclusion = 0.0, extra = 1.0), POST_FIX_PRIOR)
        # A non-positive prior mass has no logarithm; that must be an ArgumentError naming the
        # function rather than a silent -Inf that quietly zeroes a class.
        @test_throws ArgumentError p14_class_posterior(POST_FIX_FLAT_LOGBF,
                                                       (coloc = 0.5, random = 0.5, exclusion = 0.0))
    end

    @testset "the D-08 structural zero is enforced with === semantics" begin
        # `===` and not `==`: -0.0 == 0.0 is true but they are different bit patterns, and the
        # claim is that the reference class was scored against ITSELF, giving exactly +0.0.
        @test_throws ArgumentError p14_class_posterior(
            (coloc = 0.0, random = 1e-300, exclusion = 0.0), POST_FIX_FLAT_PRIOR)
        @test_throws ArgumentError p14_class_posterior(
            (coloc = 0.0, random = -0.0, exclusion = 0.0), POST_FIX_FLAT_PRIOR)
        # The message names the decision, so a reader who hits it lands on D-08 rather than on
        # a guess about what "structural zero" meant.
        err = try
            p14_class_posterior((coloc = 0.0, random = 1e-300, exclusion = 0.0),
                                POST_FIX_FLAT_PRIOR)
        catch e
            e
        end
        @test occursin("D-08", err.msg)
    end

    @testset "no overflow: log-sum-exp is not optional" begin
        # Raw log Bayes factors span roughly +-22 nats on this net, and `exp` on a raw value of
        # that size is already large; 800 is the deliberately extreme probe. Without the
        # max-subtraction this is Inf/Inf = NaN.
        p = p14_class_posterior((coloc = 800.0, random = 0.0, exclusion = 0.0), POST_FIX_PRIOR)
        @test isfinite(p.coloc) && isfinite(p.random) && isfinite(p.exclusion)
        @test isapprox(p.coloc + p.random + p.exclusion, 1.0; atol = 1e-12)
        @test p.coloc == 1.0                     # saturates AT one, never above it
        # The mirror case underflows instead of overflowing, and must not produce a NaN either.
        m = p14_class_posterior((coloc = -800.0, random = 0.0, exclusion = -800.0),
                                POST_FIX_PRIOR)
        @test isfinite(m.coloc) && isfinite(m.random) && isfinite(m.exclusion)
        @test m.random == 1.0
    end

    @testset "the composite null is {random OR exclusion}, and it is decomposable" begin
        p = p14_class_posterior(POST_FIX_LOGBF, POST_FIX_PRIOR)
        v = p14_null_posterior(p)
        @test v === 1.0 - p.coloc
        @test isapprox(v, p.random + p.exclusion; atol = 1e-12)
        # A batch whose false discoveries are mostly EXCLUSION mass is scientifically different
        # from one whose are mostly RANDOM mass, so the split travels beside the total.
        s = p14_null_split(p)
        @test s.v === v
        @test s.p_random    === p.random
        @test s.p_exclusion === p.exclusion
        @test isapprox(s.p_random + s.p_exclusion, s.v; atol = 1e-12)
    end

    @testset "confidence is the three-term max, taken by name" begin
        p = p14_class_posterior(POST_FIX_LOGBF, POST_FIX_PRIOR)
        @test p14_confidence(p) === p.coloc
        @test isapprox(p14_confidence(p), 0.7; atol = 1e-12)
        # Each class in turn, so the accessor cannot be a constant that happens to be right on
        # one fixture.
        e = p14_class_posterior((coloc = 0.0, random = 0.0, exclusion = 5.0), POST_FIX_PRIOR)
        @test p14_confidence(e) === e.exclusion
        r = p14_class_posterior((coloc = -5.0, random = 0.0, exclusion = -5.0), POST_FIX_PRIOR)
        @test p14_confidence(r) === r.random
        @test p14_confidence(POST_FIX_FLAT_PRIOR) === 1/3
    end

    @testset "the prior is RE-DERIVABLE, not taken on trust" begin
        # Assumptions Log A7: the FDR guarantee rests on `pi_class_masses` having been computed
        # under the inherited tau with the D-05 two-factor cut, and that claim was corroborated
        # only by a table in a docstring. This re-derives it from the simulator prior and the
        # frozen label rule -- cheap, because the D-05 label depends only on theta and needs no
        # forward simulation at all.
        m = p14_rederive_class_masses(n = POST_FIX_N_REDERIVE)
        @test keys(m) === P14_CLASS_KEYS
        @test isapprox(m.coloc + m.random + m.exclusion, 1.0; atol = 1e-12)
        println("    re-derived class masses at n = $(POST_FIX_N_REDERIVE): $m")
        println("    h.meta.pi_class_masses:                      $(POST_FIX_PI)")
        for k in P14_CLASS_KEYS
            p_ref = getproperty(POST_FIX_PI, k)
            p_hat = getproperty(m, k)
            three_se = 3 * sqrt(p_ref * (1 - p_ref) / POST_FIX_N_REDERIVE)
            println("      $k: pi = $p_ref  re-derived = $p_hat  |diff| = " *
                    "$(abs(p_hat - p_ref))  3 binomial SE = $three_se")
            @test abs(p_hat - p_ref) <= three_se
        end
        # The re-derivation rides the FIXTURE stream, never a reported one, so it can never
        # pre-observe draws a reported number will later ride.
        @test P14_FIXTURE_COUNTER ∉ P14_REPORTED_COUNTERS
        # Reproducible: same construction, same masses.
        @test p14_rederive_class_masses(n = 2_000) == p14_rederive_class_masses(n = 2_000)
        @test_throws ArgumentError p14_rederive_class_masses(n = 0)
    end

    @testset "the measured prior composes with the measured evidence" begin
        # The end-to-end shape the decision layer will actually call: the MEASURED prior, and a
        # log Bayes factor triple in the key order `three_way_log_bf` emits it in.
        p = p14_class_posterior((coloc = 2.5, exclusion = -3.25, random = 0.0),
                                (coloc     = POST_FIX_PI.coloc,
                                 random    = POST_FIX_PI.random,
                                 exclusion = POST_FIX_PI.exclusion))
        @test isapprox(p.coloc + p.random + p.exclusion, 1.0; atol = 1e-12)
        @test p.coloc > p.random > p.exclusion
        @test isapprox(p14_null_posterior(p), p.random + p.exclusion; atol = 1e-12)
    end

    @testset "P14 posterior ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
