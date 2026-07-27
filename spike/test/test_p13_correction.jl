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

# spike/test/test_p13_correction.jl --- D-07: the three-way evidence-scale correction, VERIFIED.
#
# THIS IS THE ONLY CHECK IN THE PHASE THAT MEETS AN ANALYTIC GROUND TRUTH. Everywhere else a
# number is compared against another estimate: the gate compares the net against simulator
# labels, the continuity check compares it against a second net, and the retired KDE baseline was
# itself invalid past |logBF| ~ 6.9 (named limit 3). D-07 refuses to let the three-way correction
# be assumed to generalize from the binary formula, and the conjugate-Gaussian toy of
# spike/p13/toy_gaussian.jl is the only construction here whose three-way log Bayes factors are
# known exactly. So this file is where the algebra either survives contact with truth or does
# not.
#
# WHAT IS BEING VERIFIED, IN ONE LINE: for a plain BCE head separating A from B at training
# frequencies q_A, q_B, the log Bayes factor is the logit MINUS log(q_A/q_B) -- with q_A and q_B
# counted over that head's OWN restricted class pair (13-RESEARCH F1), and exactly so only when
# stratification changed class FREQUENCIES and not within-class shape (13-RESEARCH F2). Both
# clauses are tested; the second is tested by making it FAIL.
#
# BOTH NEGATIVE CONTROLS ARE MANDATORY (P13_F5_NEGCTRL_REQUIRED). Control 1 rules out a vacuous
# pass; control 2 rules out an asserted design choice. See spike/p13/toy_gaussian.jl's banner.
#
# THE BARS ARE READ FROM THE FROZEN PRE-REGISTRATION AND ARE NEVER RELAXED. If a bar is missed
# the shortfall is recorded in the plan SUMMARY as an honest finding; it is NOT repaired by
# editing this file, the toy, the widths or the epochs. That is threat T-13-25, and the whole
# point of freezing P13_F5_CORR_MIN / P13_F5_MAXABS_TOL before any run was to make the failure
# mode visible instead of negotiable.
#
# Mirrors test_p13_net.jl / test_sbc.jl: license header -> using Test -> guarded include of the
# unit under test -> fixtures computed ONCE at module level -> one outer @testset. All draws ride
# the FIXTURE stream, so this file never consumes a reported Phase-13 counter.
#
# RUNNABLE BEFORE THE TIER-2 TAU COMMIT: the toy carries its own TOY_TAU and never calls the
# p13_tau() accessor, which errors by design until plan 13-10 measures the real tau.

using Test
using Statistics
using Random

# Unit under test: the closed-form toy (pulls p13/consts.jl -> p13/labels.jl -> p13/net.jl).
isdefined(@__MODULE__, :toy_analytic_log_bf) ||
    include(joinpath(@__DIR__, "..", "p13", "toy_gaussian.jl"))

# --- Test constants, NOT pre-registered thresholds -------------------------------------------
# These describe the TOY, not the gate. They were chosen for numerical conditioning and for a
# fast file, and they are stated here so a reader can see at a glance which numbers in this file
# are frozen pre-registration (every P13_* name, read from spike/p13/consts.jl) and which are
# fixture parameters (every TOY_* name, declared right here).
#
# TOY_TAU / TOY_S / TOY_SIGMA: the cut sits at 0.6 prior standard deviations and the observation
# noise is half a prior standard deviation, so all three classes carry real mass, the log Bayes
# factors span several nats across the sample, and the closed-form CDF differences stay far from
# underflow over the whole evaluated band.
const TOY_TAU   = 0.30
const TOY_S     = 1.0
const TOY_SIGMA = 0.50
# TOY_N is in the low thousands so the whole file stays a unit test. Under P13_F5_SKEW_FREQ the
# coloc head -- the starved one, by design -- still trains on 0.40*TOY_N samples.
const TOY_N     = 6000
const TOY_EPOCHS    = 200
const TOY_BATCHSIZE = 128
# Deliberately smaller than the D-10 recipe's P13_SUMMARY_WIDTH / P13_NUM_SUMMARIES because the
# toy's input is ONE-dimensional. Passed explicitly at every call site, so nothing here can be
# mistaken for a change to the pre-registered architecture.
const TOY_WIDTH = 32
const TOY_NSUM  = 8
const TOY_GRID_N = 4000
# The determinism check needs only that two identically-seeded fits agree, which is a property of
# the seeding and not of the epoch count, so it runs a short fit.
const TOY_DET_EPOCHS = 20

# --- Fixtures, computed ONCE ------------------------------------------------------------------
# EVERY draw rides p13_fix_rng(P13_FIXTURE_COUNTER), the FIXTURE stream, so this gate never
# consumes (and so never pre-observes) a reported Phase-13 stream (D-01). The stream is threaded
# sequentially through the training draw and then the held-out draw, which is what makes the
# evaluation grid genuinely held out rather than a re-draw of the same sub-stream.
const FIX_RNG = p13_fix_rng(P13_FIXTURE_COUNTER)

const TOY_TRAIN = toy_draw_dataset(TOY_N; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA,
                                   rng = FIX_RNG)
const TOY_EVAL  = toy_draw_dataset(TOY_GRID_N; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA,
                                   rng = FIX_RNG)
# The within-class RESHAPED training set of negative control 2, at the SAME class frequencies.
const TOY_RESHAPED = toy_reshaped_dataset(TOY_N; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA,
                                          rng = FIX_RNG)

# Each head is scored on its OWN class pair (D-11); see toy_head_grid's docstring.
const GRID_C = toy_head_grid(TOY_EVAL, COLOC)
const GRID_E = toy_head_grid(TOY_EVAL, EXCLUSION)

_fit(ds; epochs = TOY_EPOCHS) =
    toy_fit_two_head(ds; epochs = epochs, batchsize = TOY_BATCHSIZE,
                     width = TOY_WIDTH, num_summaries = TOY_NSUM, verbose = false)

# `train_three_way` calls Random.seed! from its `seed` immediately before `train`, and
# `toy_fit_two_head` calls it again immediately before the layers are built -- layer 2 of
# P13_GLOBAL_RNG_DISCIPLINE, required because Flux's initializer AND its shuffling DataLoader
# both draw from the GLOBAL RNG no matter how carefully a counter-based stream is threaded.
const FIT_PRIMARY  = _fit(TOY_TRAIN)
const FIT_RESHAPED = _fit(TOY_RESHAPED)

const HLO = FIT_PRIMARY.head_log_odds

_cmp(fit, grid; corrected = true) =
    toy_compare(fit.net, fit.head_log_odds, grid; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA,
                corrected = corrected)

# Numerically stable log-sum-exp, written out here rather than imported: LogExpFunctions is not a
# direct dependency of the spike environment and this file installs nothing.
function _logsumexp(xs)
    m = maximum(xs)
    return m + log(sum(exp.(xs .- m)))
end

@testset "P13 D-07 three-way evidence-scale correction" verbose = true begin

    @testset "analytic evidence is self-consistent" begin
        # THE CLOSED FORM IS THE MARGINAL DECOMPOSITION, not merely a plausible expression. The
        # three class evidences, mixed at the exact prior class masses, must reproduce the
        # marginal N(0, s^2 + sigma^2) density for EVERY Z. A wrong normalizer, a swapped
        # numerator and denominator, or a mishandled open-ended interval all break this identity
        # while still producing smooth, believable-looking curves.
        q = toy_prior_class_masses(; tau = TOY_TAU, s = TOY_S)
        @test isapprox(q.exclusion + q.random + q.coloc, 1.0; atol = 1e-12)
        # Symmetric prior, symmetric cut: the two tails carry equal mass.
        @test isapprox(q.exclusion, q.coloc; atol = 1e-12)

        marg = Normal(0.0, sqrt(TOY_S^2 + TOY_SIGMA^2))
        for Z in (-3.5, -2.0, -0.4, 0.0, 0.4, 1.5, 3.5)
            lE = toy_log_evidence(Z, -Inf, -TOY_TAU; s = TOY_S, sigma = TOY_SIGMA)
            lR = toy_log_evidence(Z, -TOY_TAU, TOY_TAU; s = TOY_S, sigma = TOY_SIGMA)
            lC = toy_log_evidence(Z,  TOY_TAU,  Inf; s = TOY_S, sigma = TOY_SIGMA)
            mix = _logsumexp([log(q.exclusion) + lE, log(q.random) + lR, log(q.coloc) + lC])
            @test isapprox(mix, logpdf(marg, Z); atol = 1e-8)
        end

        # SIGN SYMMETRY, in its EXACT form: reflecting Z swaps the two heads. It is an equality
        # between heads at reflected arguments, NOT a negation of either -- at Z = 2 the coloc
        # entry is about +5.83 and the exclusion entry about -5.61, which are not negatives of
        # one another, so a negation form would be satisfiable only by zero.
        for Z in (0.5, 1.25, 2.0, 3.0)
            f, g = toy_analytic_log_bf(Z; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA),
                   toy_analytic_log_bf(-Z; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA)
            @test isapprox(f.coloc, g.exclusion; atol = 1e-8)
            @test isapprox(f.exclusion, g.coloc; atol = 1e-8)
            # The random reference is identically zero by construction (D-08), not measured.
            @test f.random === 0.0 && g.random === 0.0
        end

        # The evidence is monotone in the direction each class points, which is the weakest
        # sanity property a reader would check by hand.
        @test toy_analytic_log_bf(2.0; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA).coloc >
              toy_analytic_log_bf(0.0; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA).coloc
        @test toy_analytic_log_bf(0.0; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA).coloc >
              toy_analytic_log_bf(-2.0; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA).coloc

        # Degenerate inputs throw rather than silently returning NaN.
        @test_throws ArgumentError toy_log_evidence(0.0, 1.0, 1.0; s = TOY_S, sigma = TOY_SIGMA)
        @test_throws ArgumentError toy_log_evidence(0.0, -1.0, 1.0; s = 0.0, sigma = TOY_SIGMA)
        @test_throws ArgumentError toy_analytic_log_bf(0.0; tau = 0.0, s = TOY_S, sigma = TOY_SIGMA)
    end

    @testset "the bars come from the frozen consts" begin
        # THE RUN IS SCORED AGAINST COMMITTED VALUES AND NEVER TUNES THEM. Asserting the literals
        # here -- and only here -- is what makes an after-the-fact edit of spike/p13/consts.jl
        # show up as a test failure rather than as a quietly better result (threat T-13-25).
        @test P13_F5_CORR_MIN == 0.99
        @test P13_F5_MAXABS_TOL == 0.25
        @test P13_F5_CENTRAL_FRAC == 0.90
        @test P13_F5_NEGCTRL_REQUIRED == true
        # The toy's class frequencies are the frozen skewed ones, and they are lopsided enough
        # that the per-head correction is several times the tolerance it must be measured against.
        @test P13_F5_SKEW_FREQ == (0.60, 0.25, 0.15)
        @test isapprox(sum(P13_F5_SKEW_FREQ), 1.0; atol = 1e-12)
        @test abs(log(P13_F5_SKEW_FREQ[3] / P13_F5_SKEW_FREQ[2])) > P13_F5_MAXABS_TOL
        @test abs(log(P13_F5_SKEW_FREQ[1] / P13_F5_SKEW_FREQ[2])) > P13_F5_MAXABS_TOL
        # The design the toy's primary dataset implements, and its pre-declared fallback.
        @test P13_STRATIFICATION === :class_frequency
        @test P13_STRATIFICATION_FALLBACK === :importance_weighted
    end

    @testset "the toy realizes the frozen skewed class frequencies" begin
        m = class_masses(TOY_TRAIN.classes)
        @test isapprox(m.exclusion, P13_F5_SKEW_FREQ[1]; atol = 1e-3)
        @test isapprox(m.random,    P13_F5_SKEW_FREQ[2]; atol = 1e-3)
        @test isapprox(m.coloc,     P13_F5_SKEW_FREQ[3]; atol = 1e-3)
        # Within-class shape is pi(theta | class) EXACTLY: every theta lies in its own class.
        @test all(j -> toy_class_of(TOY_TRAIN.theta[j]; tau = TOY_TAU) === TOY_TRAIN.classes[j],
                  eachindex(TOY_TRAIN.theta))
        # The reshaped control differs in within-class SHAPE and in nothing else: same
        # frequencies, same class membership, different exclusion-class distribution.
        mr = class_masses(TOY_RESHAPED.classes)
        @test (mr.exclusion, mr.random, mr.coloc) == (m.exclusion, m.random, m.coloc)
        @test all(j -> toy_class_of(TOY_RESHAPED.theta[j]; tau = TOY_TAU) ===
                       TOY_RESHAPED.classes[j], eachindex(TOY_RESHAPED.theta))
        exc_p = TOY_TRAIN.theta[findall(==(EXCLUSION), TOY_TRAIN.classes)]
        exc_r = TOY_RESHAPED.theta[findall(==(EXCLUSION), TOY_RESHAPED.classes)]
        # The truncated normal piles up just below -tau; the uniform proposal does not.
        @test mean(exc_r) < mean(exc_p) - 0.5
        # The classes NOT reshaped are drawn by the identical code path in both datasets.
        @test minimum(TOY_RESHAPED.theta[findall(==(COLOC), TOY_RESHAPED.classes)]) > TOY_TAU
    end

    @testset "CORRECTED heads recover the analytic log BF" begin
        cC = _cmp(FIT_PRIMARY, GRID_C)
        cE = _cmp(FIT_PRIMARY, GRID_E)
        println("  [D-07 CORRECTED]   coloc     corr = ", round(cC.coloc.corr; digits = 5),
                "  maxabs = ", round(cC.coloc.maxabs; digits = 4),
                "  meandev = ", round(cC.coloc.meandev; digits = 4),
                "  maxabs_demeaned = ", round(cC.coloc.maxabs_demeaned; digits = 4),
                "  (n = ", cC.n, ", band = [", round(cC.lo; digits = 3), ", ",
                round(cC.hi; digits = 3), "])")
        println("  [D-07 CORRECTED]   exclusion corr = ", round(cE.exclusion.corr; digits = 5),
                "  maxabs = ", round(cE.exclusion.maxabs; digits = 4),
                "  meandev = ", round(cE.exclusion.meandev; digits = 4),
                "  maxabs_demeaned = ", round(cE.exclusion.maxabs_demeaned; digits = 4),
                "  (n = ", cE.n, ", band = [", round(cE.lo; digits = 3), ", ",
                round(cE.hi; digits = 3), "])")
        println("  [D-07 head_log_odds] ", HLO)

        # Enough evaluation points for the central-band statistics to mean anything.
        @test cC.n > 500 && cE.n > 500

        # THE PRE-REGISTERED BARS, both heads, read from the frozen consts.
        @test cC.coloc.corr       >= P13_F5_CORR_MIN
        @test cE.exclusion.corr   >= P13_F5_CORR_MIN
        @test cC.coloc.maxabs     <= P13_F5_MAXABS_TOL
        @test cE.exclusion.maxabs <= P13_F5_MAXABS_TOL

        # THE SUBSTANTIVE CLAIM, stated separately from the max: after the correction the SIGNED
        # mean deviation sits inside the tolerance, i.e. the systematic offset the correction
        # exists to remove is gone. The max is a tail statistic and can miss for reasons that
        # have nothing to do with the algebra; the mean cannot.
        @test abs(cC.coloc.meandev)     <= P13_F5_MAXABS_TOL
        @test abs(cE.exclusion.meandev) <= P13_F5_MAXABS_TOL
    end

    @testset "NEGATIVE CONTROL 1: the UNCORRECTED logit fails the same bars" begin
        # WITHOUT THIS CONTROL THE SUITE CANNOT DISTINGUISH "the correction works" from "the
        # correction is unnecessary here". P13_F5_SKEW_FREQ exists to make the omission large.
        # Same trained net, same grids, same statistics -- `corrected` is a PARAMETER of
        # toy_compare, so the two arms are provably one computation differing in one subtraction.
        uC = _cmp(FIT_PRIMARY, GRID_C; corrected = false)
        uE = _cmp(FIT_PRIMARY, GRID_E; corrected = false)
        cC = _cmp(FIT_PRIMARY, GRID_C)
        cE = _cmp(FIT_PRIMARY, GRID_E)
        println("  [D-07 UNCORRECTED] coloc     corr = ", round(uC.coloc.corr; digits = 5),
                "  maxabs = ", round(uC.coloc.maxabs; digits = 4),
                "  meandev = ", round(uC.coloc.meandev; digits = 4))
        println("  [D-07 UNCORRECTED] exclusion corr = ", round(uE.exclusion.corr; digits = 5),
                "  maxabs = ", round(uE.exclusion.maxabs; digits = 4),
                "  meandev = ", round(uE.exclusion.meandev; digits = 4))

        @test P13_F5_NEGCTRL_REQUIRED

        # (i) THE CONTROL MUST FAIL THE MAX BAR -- on BOTH heads here, since both frozen
        #     frequencies are lopsided enough to guarantee it.
        @test uC.coloc.maxabs     > P13_F5_MAXABS_TOL
        @test uE.exclusion.maxabs > P13_F5_MAXABS_TOL

        # (ii) AND IT MUST FAIL BY THE SPECIFIC EXPECTED OFFSET, not by some unspecified amount.
        #      The uncorrected signed mean deviation must land on the head's OWN measured
        #      log-odds. This is what rules out a wrong sign, a wrong denominator (all three
        #      classes instead of the head's pair) or a double application -- each of which would
        #      also "fail", and none of which would land here (13-RESEARCH Pitfall 1).
        @test isapprox(uC.coloc.meandev,     HLO.coloc;     atol = 0.15)
        @test isapprox(uE.exclusion.meandev, HLO.exclusion; atol = 0.15)
        # And the offset is exactly the subtraction, to floating point: this is the identity
        # uncorrected = corrected + log(q_pos/q_neg).
        @test isapprox(uC.coloc.meandev - cC.coloc.meandev, HLO.coloc; atol = 1e-9)
        @test isapprox(uE.exclusion.meandev - cE.exclusion.meandev, HLO.exclusion; atol = 1e-9)

        # (iii) WHAT THE CONTROL CANNOT FAIL ON, asserted so a future reader does not add a
        #       correlation-only gate. A missing correction is a CONSTANT shift and correlation
        #       is shift-invariant, so the uncorrected correlation is bit-comparable to the
        #       corrected one. A verification that gated on correlation alone would have passed a
        #       completely uncorrected net.
        @test isapprox(uC.coloc.corr,     cC.coloc.corr;     atol = 1e-9)
        @test isapprox(uE.exclusion.corr, cE.exclusion.corr; atol = 1e-9)
    end

    @testset "NEGATIVE CONTROL 2: the scalar fails under a within-class reshaped proposal" begin
        # THIS IS 13-RESEARCH F2's THEOREM MADE RUNNABLE, and its failure is the RESULT.
        #
        # The theorem is an IFF: the scalar correction -log(q_A/q_B) recovers the pi-scale Bayes
        # factor if and only if q(theta|A) == pi(theta|A) and q(theta|B) == pi(theta|B). The
        # reshaped dataset holds the class FREQUENCIES fixed -- so the measured scalar is
        # unchanged in kind -- and replaces the exclusion class's within-class SHAPE. The change
        # then enters INSIDE the evidence integral p_q(Z|A) = int p(Z|theta) q(theta|A) dtheta,
        # where it is a function of Z. No constant can repair a function of Z, so the exclusion
        # head's corrected output is wrong by a Z-DEPENDENT amount while the coloc head -- whose
        # own two classes were untouched -- is unaffected.
        #
        # WHY THIS TESTSET IS NOT OPTIONAL. Without it, the ruling
        # P13_STRATIFICATION = :class_frequency (and the rejection of D-07's literal within-range
        # rho-stratification, which is a within-class reshaping) would rest on the derivation
        # alone. D-07 says the correction must be DERIVED AND VERIFIED. This is the verification
        # of the second clause, and it is also why P13_STRATIFICATION_FALLBACK is
        # :importance_weighted rather than "a different scalar": a reshaped proposal needs a
        # per-sample weight pi(theta)/g(theta), which is not a scalar at all.
        rE = _cmp(FIT_RESHAPED, GRID_E)
        rC = _cmp(FIT_RESHAPED, GRID_C)
        println("  [D-07 RESHAPED]    coloc     corr = ", round(rC.coloc.corr; digits = 5),
                "  maxabs = ", round(rC.coloc.maxabs; digits = 4),
                "  meandev = ", round(rC.coloc.meandev; digits = 4),
                "  maxabs_demeaned = ", round(rC.coloc.maxabs_demeaned; digits = 4))
        println("  [D-07 RESHAPED]    exclusion corr = ", round(rE.exclusion.corr; digits = 5),
                "  maxabs = ", round(rE.exclusion.maxabs; digits = 4),
                "  meandev = ", round(rE.exclusion.meandev; digits = 4),
                "  maxabs_demeaned = ", round(rE.exclusion.maxabs_demeaned; digits = 4))
        println("  [D-07 RESHAPED head_log_odds] ", FIT_RESHAPED.head_log_odds)

        # THE SCALAR ITSELF DID NOT MATERIALLY MOVE, which is what makes the failure below
        # attributable to the within-class SHAPE and to nothing else. Both arms measure a
        # log-odds consistent with the frozen frequency ratio; they differ only by the sampling
        # noise of two different shuffles of the same class counts into the training split, which
        # is a hundredth of the deviation the control is about to exhibit.
        freq_ratio_E = log(P13_F5_SKEW_FREQ[1] / P13_F5_SKEW_FREQ[2])
        freq_ratio_C = log(P13_F5_SKEW_FREQ[3] / P13_F5_SKEW_FREQ[2])
        @test isapprox(FIT_RESHAPED.head_log_odds.exclusion, freq_ratio_E; atol = 0.1)
        @test isapprox(HLO.exclusion, freq_ratio_E; atol = 0.1)
        @test isapprox(FIT_RESHAPED.head_log_odds.coloc, freq_ratio_C; atol = 0.1)
        @test isapprox(HLO.coloc, freq_ratio_C; atol = 0.1)

        # THE CONTROL: the reshaped head must now FAIL the pre-registered bar.
        @test rE.exclusion.maxabs > P13_F5_MAXABS_TOL

        # AND NO SCALAR COULD HAVE SAVED IT. `maxabs_demeaned` removes the BEST constant in
        # hindsight -- a better correction than any that could be measured -- and the residual
        # still exceeds the tolerance. That is F2's theorem as a number: the discrepancy is a
        # function of Z, and a function of Z is not repairable by a constant. Asserting only
        # `maxabs` would leave open the reply "then measure a better constant".
        @test rE.exclusion.maxabs_demeaned > P13_F5_MAXABS_TOL

        # It is materially worse than the SAME head on the un-reshaped arm, so the control is
        # measuring the reshaping and not generic estimation noise.
        cE = _cmp(FIT_PRIMARY, GRID_E)
        @test rE.exclusion.maxabs > cE.exclusion.maxabs
        @test abs(rE.exclusion.meandev) > abs(cE.exclusion.meandev)
    end

    @testset "the correction is the head's own denominator" begin
        # THE SILENT D-08 BREAK. A one-vs-rest denominator -- n_random + n_exclusion under the
        # coloc head -- produces code that looks identical and publishes a MIXTURE Bayes factor
        # under D-08's name. It is asserted absent by arithmetic, not by reading.
        cl = FIT_PRIMARY.classes_train
        nC = count(==(COLOC), cl)
        nR = count(==(RANDOM), cl)
        nE = count(==(EXCLUSION), cl)
        @test nC > 0 && nR > 0 && nE > 0

        @test isapprox(HLO.coloc,     log(nC / nR); atol = 1e-12)
        @test isapprox(HLO.exclusion, log(nE / nR); atol = 1e-12)
        # EXPLICITLY NOT the one-vs-rest form.
        @test !isapprox(HLO.coloc,     log(nC / (nR + nE)); atol = 1e-6)
        @test !isapprox(HLO.exclusion, log(nE / (nR + nC)); atol = 1e-6)
        # And not the whole-set form the shipped binary read surface measures, either.
        @test !isapprox(HLO.coloc, log(nC / (nC + nR + nE)); atol = 1e-6)

        # MEASURED ON THE TRAINING SUBSET ONLY. Measuring it on the full dataset would fold an
        # evaluation-set artifact into every reported Bayes factor (threat T-13-23), so the two
        # counts are asserted to actually differ -- otherwise this assertion would be vacuous.
        @test FIT_PRIMARY.n_train < TOY_N
        @test length(cl) == FIT_PRIMARY.n_train
        full = head_log_odds(TOY_TRAIN.classes)
        @test HLO.coloc != full.coloc || HLO.exclusion != full.exclusion

        # RANDOM is the shared negative of BOTH heads -- the structural reason the two logits are
        # ratios against the SAME reference (D-08).
        @test head_log_odds(cl) == HLO
        @test measure_head_log_odds(cl; positive = COLOC, negative = RANDOM) == HLO.coloc
    end

    @testset "a maxabs shortfall is not a training shortfall" begin
        # A DIAGNOSTIC THAT MAKES A MISS ATTRIBUTABLE. If a bar above is missed, the very first
        # question is whether the toy net simply had not converged -- in which case the miss says
        # nothing about the correction. This testset answers it with a number: the analytic
        # Bayes-optimal logits (log BF + the head's log-odds, i.e. the exact inverse of the D-07
        # correction) are scored under the SAME masked loss on the SAME held-out data as the
        # trained net. Their loss is the Bayes risk of the classification problem.
        #
        # It also demonstrates why a max-in-nats bar is hard here at all: logit-BCE is nearly
        # FLAT in the confident tail, so a net can sit essentially AT the Bayes risk while still
        # being visibly off in nats where the head is already right with probability ~0.99.
        tg = target_matrix(TOY_EVAL.classes)
        ideal = toy_ideal_logits(TOY_EVAL.Z, HLO; tau = TOY_TAU, s = TOY_S, sigma = TOY_SIGMA)
        keep = [all(isfinite, view(ideal, :, j)) for j in 1:size(ideal, 2)]
        # The closed form underflows only in the deep tail, so essentially everything survives.
        @test count(keep) >= size(ideal, 2) - 5

        loss_net   = Float64(masked_two_head_bce(FIT_PRIMARY.net(TOY_EVAL.Z[:, keep]), tg[:, keep]))
        loss_ideal = Float64(masked_two_head_bce(ideal[:, keep], tg[:, keep]))
        println("  [D-07 BCE] held-out masked loss: net = ", round(loss_net; digits = 6),
                "   Bayes-optimal = ", round(loss_ideal; digits = 6),
                "   excess = ", round(loss_net - loss_ideal; digits = 6))

        @test isfinite(loss_net) && isfinite(loss_ideal)
        # The Bayes-optimal logits really are optimal, up to the sampling noise of one held-out
        # draw -- so this reference is a reference and not an arbitrary comparison.
        @test loss_ideal <= loss_net + 1e-3
        # And the trained net is within a fraction of a percent of it: whatever the max-deviation
        # numbers say, the net is not under-trained.
        @test loss_net - loss_ideal < 0.01
    end

    @testset "determinism" begin
        # Two fits from the SAME seed on the SAME data must agree. Both stochastic parts of a fit
        # -- the Flux initializer and the shuffling DataLoader -- draw from the GLOBAL RNG, so
        # this asserts layer 2 of P13_GLOBAL_RNG_DISCIPLINE is actually in force and not merely
        # documented. Perturbing the global RNG between the two fits is deliberate: if the
        # re-seeding were missing, that perturbation is exactly what would leak in.
        a = _fit(TOY_TRAIN; epochs = TOY_DET_EPOCHS)
        Random.seed!(1234)
        randn(64)
        b = _fit(TOY_TRAIN; epochs = TOY_DET_EPOCHS)

        Zg = reshape(Float32.(GRID_C), 1, :)
        @test maximum(abs.(a.net(Zg) .- b.net(Zg))) <= 1e-6
        @test a.head_log_odds == b.head_log_odds
        # The comparison statistics are therefore reproducible too.
        ca = _cmp(a, GRID_C)
        cb = _cmp(b, GRID_C)
        @test isapprox(ca.coloc.maxabs, cb.coloc.maxabs; atol = 1e-6)
        @test isapprox(ca.coloc.corr,   cb.coloc.corr;   atol = 1e-6)
    end

    @testset "P13 correction ran CPU-only (D-10)" begin
        @test P13_USE_GPU == false
        # CUDA is never imported by this path; the toy is deliberately small enough that CPU is
        # the whole story (CLAUDE.md: CPU-only is the reproducible baseline).
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
        @test_throws ArgumentError three_way_log_bf(FIT_PRIMARY.net, Float32[0.0], HLO;
                                                    use_gpu = true)
    end

end
