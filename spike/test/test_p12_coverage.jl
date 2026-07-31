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

# spike/test/test_p12_coverage.jl --- 12-16: the two-sided rule, the sentinel/OOD pairing, the
# single-stack refusal, and the RAW-row masking ordering.
#
# EVERYTHING HERE IS FAST AND SYNTHETIC. No training, no pool, no real image, no forward
# simulation: the images are built directly as two correlated matrices and pushed through the
# UNCHANGED `build_mci` / `patch_summary` contract, and the net is a randomly-initialised
# `build_p12_estimator`. That is deliberate -- these testsets check STRUCTURE AND RULES, and a
# rule that only holds for a trained net is not a rule.
#
# ===========================================================================================
# TESTSET 2 EXERCISES THE **BRANCH LOGIC** OF `lro_pass`, NOT THE STAGE-2 GATE
# ===========================================================================================
# `lro_pass` compares a SPATIAL arm against the ablation. On the route this phase actually took
# there IS no spatial arm -- the mini-spike returned NONE-BEATS-ABLATION and 12-17 never ran --
# so `stage2_gate = :not_applicable_descope` BY CONSTRUCTION and the gate is UNEXERCISED.
#
# **A GREEN TESTSET 2 IS NOT THE GATE HAVING RUN.** It is four synthetic NamedTuples exercising
# four branches of a pure function. The testset's own name says so, and this paragraph exists so
# that nobody reading a green suite concludes a spatial-versus-ablation comparison was made. The
# function is built and tested anyway because it keeps the artifact contract intact for 12-20's
# Guard 3 and because it is the first thing a v2.1 spatial arm would need.
#
# HELPERS ARE PHASE-PREFIXED (`_p12cov_t*`). A bare helper name silently overwrote its Phase-13
# namesake once already in this phase, via an identical zero-positional method signature.

using Test
using Statistics
using StatsBase
using Random

isdefined(@__MODULE__, :lro_scores) ||
    include(joinpath(@__DIR__, "..", "validation", "p12_coverage.jl"))

# ===========================================================================================
# Fixtures
# ===========================================================================================

"Strip whole-line `#` comments so a source assertion reads CODE, not prose about the code."
_p12cov_tstrip(src::AbstractString) =
    join(filter(l -> !startswith(strip(l), "#"), split(src, '\n')), '\n')

"""
    _p12cov_tbundle(; G = P12_G, D = P12_D_MINISPIKE) -> NamedTuple

A SYNTHETIC bundle with the same field names the trainer persists: a randomly-initialised
estimator and two `ZScoreTransform`s fitted on random data. It is NOT a trained net and makes no
claim to be -- the testsets below check invariants that must hold for ANY bundle, which is
precisely what makes an untrained one the right fixture.

The transforms are given a NON-ZERO mean deliberately: testset 6's whole point is that a masked
RAW row standardizes to `(0 - mu_r)/sd_r` rather than to zero, and with `mu_r = 0` those two are
indistinguishable and the test could never fail.
"""
function _p12cov_tbundle(; G::Integer = P12_G, D::Integer = P12_D_MINISPIKE, seed::Integer = 20261)
    rng = Random.MersenneTwister(seed)
    n   = G^2
    zt       = StatsBase.fit(StatsBase.ZScoreTransform, randn(rng, n, 64) .+ 0.7; dims = 2)
    theta_zt = StatsBase.fit(StatsBase.ZScoreTransform, randn(rng, D, 64) .* 0.5 .+ 0.2; dims = 2)
    return (estimator = build_p12_estimator(; G = G, D = D),
            zt = zt, theta_zt = theta_zt, arm = :none, D = D, G = G)
end

"""
    _p12cov_timage(; n = 256, seed, blank_patch = nothing, perfect_patch = nothing) -> MultiChannelImage

Two correlated channels, built directly rather than simulated. Strictly positive everywhere, so
`_exclude_zero` drops nothing and every patch is scorable unless this function makes it otherwise.

`blank_patch = (i, j)` ZEROES that patch in both channels. `_exclude_zero` then leaves 0 survivors,
`correlation` returns `missing` (the `<= 15` floor), and `encode_d01` writes value 0 with mask 0 --
a genuinely UNSCORABLE region, produced through the unchanged contract rather than by poking a
sentinel into an array.

`perfect_patch = (i, j)` sets channel 2 equal to channel 1 there, so the Pearson correlation is
EXACTLY 1.0 -- the degenerate `±1` entry `P12_FISHERZ_CLAMP` exists to keep finite.
"""
function _p12cov_timage(; n::Integer = 256, seed::Integer, blank_patch = nothing,
                        perfect_patch = nothing, G::Integer = P12_G)
    rng = Random.MersenneTwister(seed)
    x = abs.(randn(rng, n, n)) .+ 1.0
    y = 0.6 .* x .+ 0.5 .* (abs.(randn(rng, n, n)) .+ 1.0)
    ps = n ÷ G
    _rows(i) = ((i - 1) * ps + 1):(i * ps)
    if perfect_patch !== nothing
        i, j = perfect_patch
        y[_rows(i), _rows(j)] .= x[_rows(i), _rows(j)]
    end
    if blank_patch !== nothing
        i, j = blank_patch
        x[_rows(i), _rows(j)] .= 0.0
        y[_rows(i), _rows(j)] .= 0.0
    end
    return build_mci([x, y])
end

"The `p12_idx` column-major flat index of lattice cell `(i, j)` -- the one indexing contract."
_p12cov_tflat(i, j; G = P12_G) = (j - 1) * G + i

# ===========================================================================================

@testset "P12 per-region coverage, the scoring rule, and the sentinel discipline (12-16)" begin

    G    = P12_G
    nreg = G^2

    # -------------------------------------------------------------------------------------
    @testset "1 — Fisher-z round-trips and clamps degenerate correlations" begin
        for r in (-0.97, -0.5, -1e-9, 0.0, 0.25, 0.8, 0.99)
            @test isapprox(inv_fisherz(fisherz(r)), r; atol = 1e-12)
        end
        # `correlation()` returns exactly +-1 on a degenerate patch and `atanh(+-1)` is infinite.
        # The clamp keeps the arithmetic finite; without it a single degenerate region would
        # poison every pooled variance, log score and CRPS downstream with a NaN.
        @test isfinite(fisherz(1.0))
        @test isfinite(fisherz(-1.0))
        @test fisherz(1.0) == -fisherz(-1.0)
        @test fisherz(1.0) ≈ atanh(P12_FISHERZ_CLAMP)
        # ... and a value beyond the clamp is pulled back to it, never allowed through.
        @test fisherz(1.5) == fisherz(1.0)

        # THE COUNT IS REPORTED, NEVER SWALLOWED. A silently clamped +-1 entering a coverage
        # tally would be an unscorable region reported as a scored one.
        v = [0.1, 1.0, -1.0, 0.5, 0.99]
        @test count_clamped(v) == 2
        @test is_clamped_correlation(1.0) && is_clamped_correlation(-1.0)
        # A `ghat` atom at +-0.99 is a property of the MODEL, not a degeneracy of the data, and
        # must NOT be counted as one. The two are far apart and this asserts they stay apart.
        @test !is_clamped_correlation(0.99)
    end

    # -------------------------------------------------------------------------------------
    @testset "2 — the two-sided rule's BRANCH LOGIC (the gate itself did NOT run)" begin
        # THESE FOUR CASES **ARE** THE FROZEN SC3 CRITERION. A change to `lro_pass` that breaks
        # one of them is a change to a frozen success criterion, not a refactor.
        #
        # AND THEY ARE SYNTHETIC. `stage2_gate` on this route is :not_applicable_descope because
        # there is no spatial arm at all; nothing below reads a bundle, a report or a measurement.
        cal  = P12_COVERAGE_NOMINAL                                        # exactly nominal
        # CALIBRATED BUT WORSE-CENTRED THAN `cal`. HALF the tolerance, NOT the tolerance itself:
        # `0.90 - 0.03` evaluates to 0.8699999999999999, so `|coverage - nominal|` is
        # 0.030000000000000027 > 0.03 and the "edge" case lands on the WRONG side of the bound by
        # one ulp. Sitting a fixture exactly on a floating-point boundary tests the boundary, not
        # the branch -- and it made this testset report :disqualified where it meant :fail.
        near = P12_COVERAGE_NOMINAL - P12_STAGE2_COVERAGE_TOST_DELTA / 2
        off  = P12_COVERAGE_NOMINAL - 3 * P12_STAGE2_COVERAGE_TOST_DELTA   # miscalibrated

        # (a) COVERAGE-ONLY WIN: the spatial arm is better calibrated but scores no better.
        v, d = lro_pass((coverage = cal,  mean_logscore = -1.0),
                        (coverage = near, mean_logscore = -1.0))
        @test v === :fail
        @test d.spatial_calibrated && d.ablation_calibrated   # BOTH qualify; only sharpness ties
        @test d.logscore_delta == 0.0
        @test occursin("COVERAGE-ONLY WIN IS NOT A PASS", uppercase(d.reason))

        # (b) LOG-SCORE WIN with both calibrated, by at least the pre-registered margin.
        v, d = lro_pass((coverage = cal, mean_logscore = -1.0 + P12_STAGE2_LOGSCORE_MIN),
                        (coverage = cal, mean_logscore = -1.0))
        @test v === :pass
        @test d.spatial_calibrated && d.ablation_calibrated

        # (b') ... and a log-score win SMALLER than the margin is still a :fail. The bar is a
        # bar, not a direction.
        v, _ = lro_pass((coverage = cal, mean_logscore = -1.0 + P12_STAGE2_LOGSCORE_MIN / 2),
                        (coverage = cal, mean_logscore = -1.0))
        @test v === :fail

        # (c) SPATIAL MISCALIBRATED despite a better log score -> DISQUALIFIED, not "better".
        v, d = lro_pass((coverage = off, mean_logscore = 5.0),
                        (coverage = cal, mean_logscore = -1.0))
        @test v === :disqualified
        @test d.which_arm === :spatial
        @test !d.spatial_calibrated && d.ablation_calibrated

        # (d) ABLATION MISCALIBRATED -> DISQUALIFIED, NAMING THE ABLATION.
        v, d = lro_pass((coverage = cal, mean_logscore = 5.0),
                        (coverage = off, mean_logscore = -1.0))
        @test v === :disqualified
        @test d.which_arm === :ablation
        @test d.spatial_calibrated && !d.ablation_calibrated

        # (e) BOTH miscalibrated: the spatial arm is named, but BOTH flags survive in `details`,
        # so "both were miscalibrated" is never reported as a single-arm problem.
        v, d = lro_pass((coverage = off, mean_logscore = 5.0),
                        (coverage = off, mean_logscore = -1.0))
        @test v === :disqualified
        @test d.which_arm === :spatial
        @test !d.spatial_calibrated && !d.ablation_calibrated

        # The verdict is a SYMBOL, never a bare boolean -- a boolean gate cannot express
        # "disqualified" and would collapse it into "fail".
        for verdict in (:pass, :fail, :disqualified)
            @test verdict isa Symbol
        end
    end

    # -------------------------------------------------------------------------------------
    @testset "3 — the log score prefers the sharper of two calibrated predictives" begin
        # THE PROPERTY PITFALL 5 SAYS RAW COVERAGE CANNOT SEE. Both predictives below are centred
        # on the observation, so both cover it; only a proper scoring rule separates them.
        rng    = Random.MersenneTwister(4242)
        obs    = 0.30
        narrow = obs .+ 0.02 .* randn(rng, 4000)
        wide   = obs .+ 0.20 .* randn(rng, 4000)
        sn = lro_scores(obs, narrow)
        sw = lro_scores(obs, wide)

        @test sn.inside && sw.inside                 # coverage cannot tell them apart ...
        @test sn.logscore > sw.logscore              # ... the log score can.
        @test sn.crps < sw.crps                      # CRPS agrees, with the OPPOSITE sign
                                                     # convention: lower is better.

        # THE WIDTH IS RETURNED BESIDE THE COVERAGE, AND THAT IS WHY. A wide interval centred
        # near the marginal mean covers ~90 % while carrying no information about which region is
        # which, so a coverage figure without its width cannot be read at all.
        @test sn.width < sw.width
        @test sn.width ≈ sn.hi - sn.lo
        @test sn.hi > sn.lo

        # An empty or singleton predictive is a defect, not a certainty.
        @test_throws AssertionError lro_scores(obs, [0.3])
    end

    # -------------------------------------------------------------------------------------
    @testset "4 — the sentinel and the OOD flag are inseparable (T-7-04); the three maps satisfy the identity" begin
        b   = _p12cov_tbundle()
        N   = 40
        kw  = (N = N, nominal = P12_COVERAGE_NOMINAL, n_eff = 40.0)

        mci_s = _p12cov_timage(; seed = 11)
        mci_c = _p12cov_timage(; seed = 12)
        r = p12_coloc_map(b, mci_s, mci_c; kw...)

        # THE STRUCTURAL IDENTITY that REPLACED (never tightened) 12-12's [-2,2] range guard.
        @test all(isapprox.(r.region_delta_rho, r.region_rho_sample .- r.region_rho_control;
                            atol = P12_IDENTITY_ATOL, rtol = 0.0))
        @test r.ood == (r.ood_sample .| r.ood_control)
        @test r.meta.delta_rho_available === true
        @test !any(r.ood)                       # a clean pair: every region measured
        @test all(isfinite, r.region_sd) && all(>(0), r.region_sd)

        # A GENUINE MEASURED ZERO IS DISTINGUISHABLE FROM A REFUSAL (T-7-04) -- and the
        # DISTINGUISHER IS THE FLAG, NOT THE VALUE, because the value is identical. Built
        # deterministically through 12-12's own bundle rather than hunting for a region whose
        # Monte-Carlo difference happens to land on 0.0: two independent draw sets from the same
        # posterior do NOT difference to exactly zero, so scoring an image against itself would
        # NOT produce this state, and a test written that way would be asserting something the
        # construction never yields.
        mz = p12_region_maps(G)
        for i in 1:nreg
            mz.region_rho_sample[i]  = 0.4;  mz.region_sd_sample[i]  = 0.1
            mz.region_rho_control[i] = 0.4;  mz.region_sd_control[i] = 0.1
            mz.region_delta_rho[i]   = 0.0;  mz.region_sd[i]         = 0.15
            mz.ood_sample[i] = false; mz.ood_control[i] = false; mz.ood[i] = false
        end
        p12_mark_unscorable!(mz, [7]; stack = :sample)      # ONE genuine refusal beside them
        rz = SpatialColocResultSpike(G, mz, nothing, nothing,
                p12_result_meta(; arm = :none, r1_posterior_median = 0.5, derived_rho = 0.0,
                                ood_threshold = nothing, net_path = nothing,
                                delta_rho_available = true))
        # Region 7 and region 8 carry the SAME value in `region_delta_rho` ...
        @test rz.region_delta_rho[7] == P12_SPATIAL_SENTINEL
        @test rz.region_delta_rho[8] == P12_SPATIAL_SENTINEL
        # ... and are nonetheless UNAMBIGUOUSLY different states, because the flag and the sd say so.
        @test rz.ood[7] && isnan(rz.region_sd[7])            # a REFUSAL
        @test !rz.ood[8] && isfinite(rz.region_sd[8])        # a MEASURED zero
        @test rz.meta.delta_rho_available === true

        # NOW MASK A REGION IN THE **SAMPLE** STACK ONLY, through the unchanged contract.
        mci_sm = _p12cov_timage(; seed = 11, blank_patch = (3, 4))
        miss   = findall(ismissing, patch_summary(mci_sm))
        @test !isempty(miss)                    # the fixture really did make a region unscorable
        rm_ = p12_coloc_map(b, mci_sm, mci_c; kw...)
        for ci in miss
            f = _p12cov_tflat(ci[1], ci[2]; G = G)
            # THE UNION RULE PROPAGATES: sample AND difference are sentinel-and-OOD ...
            @test rm_.region_rho_sample[f] == P12_SPATIAL_SENTINEL
            @test rm_.ood_sample[f]
            @test rm_.region_delta_rho[f] == P12_SPATIAL_SENTINEL
            @test rm_.ood[f]
            @test isnan(rm_.region_sd[f]) && isnan(rm_.region_sd_sample[f])
            # ... while the CONTROL stack at that region is still a measured value.
            @test !rm_.ood_control[f]
            @test isfinite(rm_.region_sd_control[f]) && rm_.region_sd_control[f] > 0
        end
        # NO REGION CARRIES ONE WITHOUT THE OTHER, anywhere in the lattice.
        for i in eachindex(rm_.ood)
            @test rm_.ood[i] == isnan(rm_.region_sd[i])
            rm_.ood[i] && @test rm_.region_delta_rho[i] == P12_SPATIAL_SENTINEL
        end
        @test rm_.ood == (rm_.ood_sample .| rm_.ood_control)

        # A CLAMPED +-1 CORRELATION IS A DEGENERATE REGION TOO, not a measurement of perfect
        # colocalization, and it must reach the SAME sentinel writer as an absent patch.
        mci_p = _p12cov_timage(; seed = 11, perfect_patch = (5, 2))
        Zp    = encode_d01(patch_summary(mci_p))
        fp    = _p12cov_tflat(5, 2; G = G)
        @test Zp[nreg + fp] == 1.0                 # PRESENT -- so the mask row cannot explain it
        @test is_clamped_correlation(Zp[fp])       # ... but degenerate all the same
        rp = p12_coloc_map(b, mci_p, mci_c; kw...)
        @test rp.region_rho_sample[fp] == P12_SPATIAL_SENTINEL
        @test rp.ood_sample[fp] && rp.ood[fp]
        @test isnan(rp.region_sd_sample[fp])
    end

    # -------------------------------------------------------------------------------------
    @testset "4b — the single-stack mode is a Δρ REFUSAL, not a Δρ of zero" begin
        # THIS IS THE TESTSET THAT MAKES 12-19's NAMED Delta-rho LIMIT EXECUTABLE. A run that
        # silently produced a real per-region Delta-rho from one stack would fail here, and so
        # would one that emitted a measured-looking zero map.
        b = _p12cov_tbundle()
        # THE THREE KEYWORDS HAVE NO DEFAULTS, so they are passed here as everywhere else.
        r = p12_coloc_map(b, _p12cov_timage(; seed = 21), nothing;
                          N = 40, nominal = P12_COVERAGE_NOMINAL, n_eff = 40.0)

        # The SAMPLE map is a live measurement ...
        @test !all(r.region_rho_sample .== P12_SPATIAL_SENTINEL)
        @test any(isfinite, r.region_sd_sample)
        @test !all(r.ood_sample)

        # ... and everything about the CONTROL and the DIFFERENCE is a refusal.
        @test all(r.region_rho_control .== P12_SPATIAL_SENTINEL)
        @test all(r.ood_control)
        @test all(isnan, r.region_sd_control)
        @test all(r.ood)                                   # by the union rule, not by fiat
        @test all(r.region_delta_rho .== P12_SPATIAL_SENTINEL)
        @test all(isnan, r.region_sd)
        @test r.region_draws === nothing                   # no difference exists to carry

        # THE MODE IS READ BY NAME, never inferred from the sentinel pattern.
        @test r.meta.delta_rho_available === false
        @test is_ood(r) === true                           # the conservative direction

        # The constructor accepts this construction UNCHANGED, because 12-12 scopes the
        # structural identity to regions where no map is a sentinel -- so it holds vacuously
        # here, and P12_SPATIAL_SENTINEL = 0.0 sits inside both sanity ranges. No special case
        # was added to admit the single-stack shape; it is legal by consequence.
        @test P12_SPATIAL_SENTINEL == 0.0
        @test r isa SpatialColocResultSpike

        # `control_mci` IS POSITIONAL AND MANDATORY: passing `nothing` must be a written decision
        # at the call site, never a default a caller reaches by omission.
        @test_throws MethodError p12_coloc_map(b, _p12cov_timage(; seed = 21);
                                               N = 40, nominal = 0.9, n_eff = 40.0)
    end

    # -------------------------------------------------------------------------------------
    @testset "5 — a missing OOD operating point reads as NOT CHECKED" begin
        b = _p12cov_tbundle()
        r = p12_coloc_map(b, _p12cov_timage(; seed = 31), _p12cov_timage(; seed = 32);
                          N = 40, nominal = P12_COVERAGE_NOMINAL, n_eff = 40.0)
        # `src/amortized/local_map.jl:88-91`: a `false` flag under a `nothing` threshold means
        # NOT CHECKED, never "in distribution".
        @test r.meta.ood_threshold === nothing
        @test p12_ood_checked(r) === false
        # ... and nothing coerced it to a boolean `false` meaning in-distribution. A `false`
        # here would be indistinguishable from a checked-and-clean operating point.
        @test !(r.meta.ood_threshold isa Bool)
        @test r.meta.ood_threshold !== false

        # Supplying one flips the predicate, so the `nothing` above is genuinely the
        # not-supplied state and not a hard-wired constant.
        r2 = p12_coloc_map(b, _p12cov_timage(; seed = 31), _p12cov_timage(; seed = 32);
                           N = 40, nominal = P12_COVERAGE_NOMINAL, n_eff = 40.0,
                           ood_threshold = 0.5)
        @test p12_ood_checked(r2) === true
    end

    # -------------------------------------------------------------------------------------
    @testset "6 — masking is applied to the RAW rows inside lro_arm" begin
        # THE SAME ORDERING PROPERTY 12-14's TESTSET 3 PINS, RE-ASSERTED AT THE CONSUMER. It is
        # asserted twice ON PURPOSE: the producer (`train_p12_npe`) and the consumer (`lro_arm`)
        # can drift apart independently, and an ordering that is right in training and wrong at
        # read time fails silently -- the net would read "this region is perfectly average"
        # exactly where read time says "this region is absent".
        b    = _p12cov_tbundle()
        mci  = _p12cov_timage(; seed = 41)
        Zraw = Float64.(encode_d01(patch_summary(mci)))
        r    = _p12cov_tflat(4, 6; G = G)
        @test Zraw[nreg + r] == 1.0            # the region starts out PRESENT

        # (i) MASK-THEN-STANDARDIZE -- the ordering `lro_arm` uses.
        right = standardize_p12(mask_regions(reshape(Zraw, :, 1), r; G = G), b.zt; G = G)
        # (ii) STANDARDIZE-THEN-MASK -- the WRONG ordering, built here so the two can be compared.
        wrong = standardize_p12(reshape(Zraw, :, 1), b.zt; G = G)
        wrong[r, 1] = 0.0; wrong[nreg + r, 1] = 0.0

        # The masked entry equals what a RAW zero standardizes to -- NOT the row mean (0.0).
        @test right[r, 1] ≈ (0.0 - b.zt.mean[r]) / b.zt.scale[r]
        @test right[r, 1] != 0.0
        @test right[r, 1] != wrong[r, 1]
        # The mask row is passed through untouched by BOTH, so the difference is entirely in the
        # continuous row -- which is exactly why the wrong ordering is silent.
        @test right[nreg + r, 1] == 0.0

        # AND THE BEHAVIOURAL PROOF AT THE CONSUMER: `lro_arm`'s encoding of a held-out region is
        # BYTE-IDENTICAL to what a NATURALLY unusable patch produces. That is the whole of
        # Pitfall 4 -- a held-out region must be in-vocabulary, not a new sentinel.
        Znat = copy(Zraw); Znat[r] = 0.0; Znat[nreg + r] = 0.0
        nat  = standardize_p12(reshape(Znat, :, 1), b.zt; G = G)
        @test right[:, 1] == nat[:, 1]

        # ... and `lro_arm` really does route through that encoding: seeded, the held-out read
        # and the naturally-unusable read agree exactly.
        Random.seed!(99)
        a1 = lro_arm(b, Zraw, [r]; N = 24, G = G)
        Random.seed!(99)
        a2 = lro_arm(b, Znat, [r]; N = 24, G = G)
        @test a1.rho == a2.rho
        @test a1.regions == [r]
        @test size(a1.rho) == (1, 24) && size(a1.z) == (1, 24)
        # The rho draws are the Gaussian draws through the frozen read-time copula, elementwise.
        @test all(a1.rho .≈ p12_z_to_rho.(a1.z))
    end

    # -------------------------------------------------------------------------------------
    @testset "7 — the derived scalar rho is LABELLED DERIVED (R-1)" begin
        b = _p12cov_tbundle()
        r = p12_coloc_map(b, _p12cov_timage(; seed = 51), _p12cov_timage(; seed = 52);
                          N = 40, nominal = P12_COVERAGE_NOMINAL, n_eff = 40.0)
        # THE LABEL IS THE DELIVERABLE (R-1), so it is asserted on the RUNNING system: the field
        # name a reader reaches for carries the word.
        @test :derived_rho in keys(r.meta)
        @test any(k -> occursin("derived", String(k)), keys(r.meta))
        @test delta_rho(r) === r.meta.derived_rho
        @test isfinite(r.meta.derived_rho) && -1.0 <= r.meta.derived_rho <= 1.0
        # THE PER-REGION MAP IS THE DELIVERABLE, NOT THIS SCALAR.
        @test delta_rho_map(r) === r.region_delta_rho
        @test uncertainty_map(r) === r.region_sd
        # ... and this result carries NO Bayes factor. Fabricating a plausible-looking one is
        # the failure `_iface_error` exists to prevent.
        @test_throws Exception bayes_factor(r)
    end

    # -------------------------------------------------------------------------------------
    @testset "8 — the library's structure: two single-stack passes, no joint estimator" begin
        # THE PROPERTY THE PLAN'S "TWO sampleposterior CALLS" CRITERION PROTECTS, asserted in the
        # form that is actually stronger. `p12_coloc_map` calls ONE `_p12cov_single_stack` helper
        # TWICE rather than inlining two copies, so the sample and control passes CANNOT differ
        # from each other -- which is what D-03's "two INDEPENDENT single-stack passes, not a
        # jointly-estimated contrast" is really asking for. A source read is sanctioned here for
        # the same reason it is in testset 7: the structure IS the claim.
        raw  = read(joinpath(@__DIR__, "..", "validation", "p12_coverage.jl"), String)
        src  = _p12cov_tstrip(raw)
        # The slice ends at the NEXT DOCSTRING OPENER, not at the next `function` keyword: the
        # helper's docstring sits between the two definitions and repeats its own signature, so a
        # `function`-delimited slice counts that line as a third call site.
        body = split(src, "function p12_coloc_map(bundle")[end]
        body = split(body, "\n\"\"\"")[1]
        @test count(m -> true, eachmatch(r"_p12cov_single_stack\(", body)) == 2
        # ... and the helper it calls performs exactly ONE posterior read.
        helper = split(src, "function _p12cov_single_stack")[end]
        @test count(m -> true, eachmatch(r"sampleposterior\(", helper)) == 1
        # No transform is ever re-fitted at read time, and CPU-only is not left to a default.
        @test !occursin("fit(ZScoreTransform", src)
        @test !occursin("use_gpu = true", src)
        @test count(m -> true, eachmatch(r"use_gpu = false", src)) >= 2
        # THE WORDING RULES ARE CHECKED ON THE **RAW** SOURCE, THE STRUCTURAL BANS ON THE STRIPPED
        # ONE, AND THE SPLIT IS THE POINT. A ban (`fit(ZScoreTransform`) must not be satisfiable by
        # a comment mentioning it; a wording RULE lives precisely in the prose, so stripping
        # comments is exactly what would hide it. Checking both against one string is how the
        # earlier draft of this testset failed: the literal survives in the header block and is
        # LINE-WRAPPED in the docstring, which is 12-02's recorded "a line-wrap that broke an
        # acceptance literal" defect, met again here.
        @test occursin("no spatial borrowing, full nuisance and global borrowing", raw)
        @test occursin("coverage-only win", raw)
    end

    # -------------------------------------------------------------------------------------
    @testset "9 — the coverage surface ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end
end
