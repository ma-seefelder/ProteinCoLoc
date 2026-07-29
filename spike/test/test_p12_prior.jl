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

# spike/test/test_p12_prior.jl --- unit gate for `spike/simulator/p12_prior.jl`, the D-05 copula.
# Run:
#     julia --project=spike -e 'using Test; include("spike/test/test_p12_prior.jl")'
#
# WHAT THIS FILE IS FOR. Every failure mode of a marginal-preserving copula is SILENT. A missing
# per-cell rescale leaves the field looking perfectly plausible and only bends the per-region
# prior. A generative order drawn backwards produces a null indistinguishable from "the
# correlation length is below the 8x8 summary's resolution". A per-pixel interpolant and a
# separable matmul produce the same picture and differ by five hours over a 50 000-pair pool.
# None of those throws. So each is asserted, and where an assertion could pass vacuously it is
# paired with a FALSIFIER -- a companion assertion that fails if the property under test is
# removed. That discipline is inherited from `test_p12_lattice.jl:25-41`.
#
# A MEASURED CORRECTION THIS FILE CARRIES, BECAUSE THE TEST IT CHANGES WOULD OTHERWISE BE UNABLE
# TO FAIL. 12-RESEARCH.md Pattern 1 and `p12_lattice.jl:41-48` both record the CAR marginal-sd
# spread as "0.658 at an edge cell against 0.992 in the interior". MEASURED THIS SESSION on the
# committed `car_sigma`, the two are the other way round: at alpha = 0.95 the CORNER cell (1,1)
# has sd 0.9922 and the INTERIOR cell (4,4) has sd 0.6575; the argmax of sqrt.(diag(Sigma)) is
# always a corner and the argmin is always an interior cell. The consequence is not cosmetic --
# the corners are the cells whose UNRESCALED marginal is already closest to 1, so they are exactly
# the regions that do NOT break when the rescale is deleted (measured corner W1 without the
# rescale at r1 = 0.5: 0.0112, against 0.0122 WITH it). A testset that asserted the falsifier only
# on the four corners could therefore never fail. Testset 2 keeps the corner assertions (they
# hold, and the boundary claim is worth stating) and adds the assertions that actually bite: the
# interior, the bottom r1 rung, and the per-cell sd of the field the copula consumes. Every one of
# those was verified BY MUTATION -- see 12-07-SUMMARY.md for the measured before/after table.
#
# SEED DISCIPLINE. All randomness rides `p12_fix_rng(P12_FIXTURE_COUNTER)` -- the FIXTURE stream
# on the FIXTURE counter (p12_consts.jl §2). Nothing here touches `p12_rng`: a test must never
# pre-observe a stream a reported number rides on.
#
# THE MONTE-CARLO SIZES ARE SCREENING SIZES. The pre-registered per-region SIM-02 claim is
# `run_p12_sim02.jl` at `P12_SIM02_M = 20 000`; this file runs M = 2 000 so a broken copula fails
# in seconds rather than in the reported run. The bar itself is NOT re-derived here: it is
# `P12_SIM02_W1_TOL_PERREGION`, read from the frozen Tier-1 pre-registration.
#
# CPU-only. No training, no file write. One `simulate_pair` call at 128^2 (testset 10).

using Test
using Statistics
using LinearAlgebra
using Distributions

# ORDER MATTERS: the Tier-1 pre-registration, then the unit under test, then the LEGACY forward
# model. The forward include is stated EXPLICITLY rather than left to another file having loaded
# it first -- `p12_prior.jl` guard-includes only `prior.jl` and `p12_lattice.jl`, and include
# order inside the aggregated suite is not a contract.
isdefined(@__MODULE__, :P12_G)            || include(joinpath(@__DIR__, "..", "validation", "p12_consts.jl"))
isdefined(@__MODULE__, :sample_p12_prior) || include(joinpath(@__DIR__, "..", "simulator", "p12_prior.jl"))
isdefined(@__MODULE__, :simulate_pair)    || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))

# --- the scoped SIM-02 reference and the W1 estimator ---------------------------------------
# SCOPED TO [GHAT_MU_MIN, GHAT_MU_MAX] EXACTLY AS `prior.jl:42-46` SCOPES ITS OWN INDUCED-mu
# CLAIM (and as `test_simulator.jl:226` implements that scoping). Conditioning a truncated Cauchy
# on a sub-interval is again a truncated Cauchy, so restricting BOTH sides to the realized range
# is exact rather than an approximation -- the comparison stays like-for-like.
const _P12P_REF = Truncated(Cauchy(0.0, 0.3), GHAT_MU_MIN, GHAT_MU_MAX)

# Wasserstein-1 by quantile transport, the `calibration.jl:137-138` / `test_simulator.jl:179-184`
# estimator, with the reference given ANALYTICALLY instead of as a second finite sample: the
# target is known in closed form here, so drawing a reference sample would only add Monte-Carlo
# noise to a screening statistic.
function _p12p_w1(x::AbstractVector{<:Real})
    s = sort([v for v in x if GHAT_MU_MIN <= v <= GHAT_MU_MAX])
    n = length(s)
    n < 5 && return NaN
    return mean(abs(s[k] - quantile(_P12P_REF, (k - 0.5) / n)) for k in 1:n)
end

# M draws of the copula with the sampler HOISTED. `sample_p12_prior` re-solves the r1 bisection on
# every call, which is right for a datagen loop of independent r1 draws and wrong for a fixed-r1
# Monte Carlo; the object under test (`rho_field`) is identical either way. Returns 64xM matrices
# in the p12_idx column-major layout.
function _p12p_fields(M::Integer; arm::Symbol = :car, r1::Real = 0.5)
    rng = p12_fix_rng(P12_FIXTURE_COUNTER)
    f   = field_sampler(lattice_sigma(arm, r1))
    RH  = Matrix{Float64}(undef, P12_G^2, M)
    MU  = similar(RH)
    ZZ  = similar(RH)
    for k in 1:M
        d = rho_field(rng, f)
        RH[:, k] = vec(d.rho); MU[:, k] = vec(d.mu); ZZ[:, k] = vec(d.z)
    end
    return (rho = RH, mu = MU, z = ZZ)
end

# The same, but through the full `sample_p12_prior` entry point (nuisances included).
function _p12p_draws(M::Integer; arm::Symbol = :car, r1 = nothing)
    rng = p12_fix_rng(P12_FIXTURE_COUNTER)
    ZZ  = Matrix{Float64}(undef, P12_G^2, M)
    MU  = similar(ZZ)
    for k in 1:M
        d = sample_p12_prior(rng; arm = arm, r1 = r1)
        ZZ[:, k] = vec(d.z_field); MU[:, k] = vec(d.mu_field)
    end
    return (z = ZZ, mu = MU)
end

_p12p_per_region_w1(MU) = [_p12p_w1(view(MU, r, :)) for r in 1:size(MU, 1)]

const _P12P_M       = 2_000
const _P12P_CORNERS = (p12_idx(1, 1), p12_idx(1, P12_G), p12_idx(P12_G, 1), p12_idx(P12_G, P12_G))
const _P12P_BOUNDARY = Set(p12_idx(i, j) for j in 1:P12_G for i in 1:P12_G
                           if i in (1, P12_G) || j in (1, P12_G))

# Drawn ONCE and shared by testsets 1-3: three passes over the same 2 000 draws, not three
# independent Monte Carlos, so the reported numbers describe one sample.
const _P12P_F05 = _p12p_fields(_P12P_M; arm = :car, r1 = 0.50)
const _P12P_W05 = _p12p_per_region_w1(_P12P_F05.mu)

@testset "the copula preserves the per-region marginal (D-05, screening)" begin
    @test length(_P12P_W05) == P12_G^2
    @test all(isfinite, _P12P_W05)
    # THE MAX OVER THE 64 REGIONS, NEVER THE MEAN. A mean lets one badly-mismatched region hide
    # behind 63 good ones, and the per-region claim is exactly what the spatial map is selling
    # (p12_consts.jl:339-344 states the same rule on the reported runner).
    @test maximum(_P12P_W05) <= P12_SIM02_W1_TOL_PERREGION
    @info "p12 screening W1 (M = $_P12P_M, arm = :car, r1 = 0.5)" max_w1 = maximum(_P12P_W05) mean_w1 = mean(_P12P_W05) argmax_region = argmax(_P12P_W05)
    # AT M = 2 000 THIS IS A SCREEN, NOT THE PRE-REGISTERED CLAIM. The pre-registered claim is
    # `run_p12_sim02.jl` at P12_SIM02_M = 20 000 draws; this testset exists so a broken copula
    # fails in seconds instead of in the reported run.
    @test P12_SIM02_M == 20_000
    @test _P12P_M < P12_SIM02_M
end

@testset "the marginal claim holds at every cell -- and the falsifier is the INTERIOR (D-05)" begin
    # (a) The boundary claim, stated on the four corner regions as the plan asks. It holds.
    for r in _P12P_CORNERS
        @test _P12P_W05[r] <= P12_SIM02_W1_TOL_PERREGION
    end

    # (b) BUT THE CORNERS ARE NOT WHERE A MISSING RESCALE SHOWS, AND THAT IS MEASURED, NOT
    # ASSUMED. `Sigma = (D - alpha*W)^-1` gives LOW-degree cells the LARGEST variance, so on this
    # lattice the corner is the cell closest to unit sd and the interior is the compressed one --
    # the reverse of the "0.658 edge / 0.992 interior" pairing recorded in 12-RESEARCH.md
    # Pattern 1 and repeated in p12_lattice.jl:41-48. Asserted structurally so the record is
    # executable rather than a footnote.
    sd = sqrt.(diag(lattice_sigma(:car, 0.50)))
    @test argmax(sd) in _P12P_CORNERS            # the LARGEST marginal sd is a corner
    @test !(argmin(sd) in _P12P_BOUNDARY)        # the SMALLEST is an interior cell
    @test sd[p12_idx(1, 1)] > sd[p12_idx(1, 4)] > sd[p12_idx(4, 4)]
    @test 0.98 < sd[p12_idx(1, 1)] < 1.00        # measured 0.9863 at r1 = 0.50
    @test 0.64 < sd[p12_idx(4, 4)] < 0.67        # measured 0.6534 at r1 = 0.50

    # (c) THE FALSIFIER FOR THE `./ sd` RESCALE, PART 1: the field the copula consumes must be
    # marginally standard PER CELL. Verified by mutation (rescale deleted): the observed sd range
    # collapses to [0.644, 1.002] and this assertion fails on 44 of the 64 regions.
    zsd = vec(std(_P12P_F05.z; dims = 2))
    @test all(s -> 0.93 <= s <= 1.07, zsd)       # M = 2 000 => sd-of-sd ~ 1.6 %, so ~4 SE
    # ... paired with the assertion that makes it non-vacuous: Sigma itself is NOT in that band,
    # so the band is a property of the SAMPLER's division and not of the covariance.
    @test !(0.93 <= minimum(sd) <= 1.07)

    # (d) THE FALSIFIER, PART 2, ON THE PRE-REGISTERED BAR ITSELF. At r1 = 0.50 a deleted rescale
    # measures a max per-region W1 of 0.0789 -- WRONG, but still inside the 0.10 bar, so testset 1
    # alone cannot catch it. At the BOTTOM rung of P12_R1_LADDER the same deletion measures
    # 0.1264 and breaches the bar. The rung is pre-registered; the bar is pre-registered; only the
    # choice to also screen at the bottom rung is this file's.
    w1_lo = _p12p_per_region_w1(_p12p_fields(_P12P_M; arm = :car, r1 = first(P12_R1_LADDER)).mu)
    @test first(P12_R1_LADDER) == P12_R1_MIN
    @test maximum(w1_lo) <= P12_SIM02_W1_TOL_PERREGION
    @info "p12 screening W1 at the bottom rung" r1 = first(P12_R1_LADDER) max_w1 = maximum(w1_lo)
end

@testset "per-region atom mass is measured and recorded, not gated (S-3)" begin
    atom_rate = count(x -> abs(x) >= 0.99 - 1e-12, _P12P_F05.rho) / length(_P12P_F05.rho)
    # A DELIBERATELY LOOSE SANITY BAND. The number that matters is REPORTED by 12-10; gating it
    # would manufacture exactly the over-powered failure mode Phase 7 already paid for twice.
    @test 0.0 <= atom_rate <= 0.25
    @info "p12 per-region ghat atom mass (rho at +/-0.99)" atom_rate = atom_rate
    # The research measured 0.0679 per region against a GLOBAL 1-D reference of 0.0678: the
    # per-region RATE is unchanged by D-05. What D-05 changes is that the atom now occurs 64x more
    # often PER DRAW, which is why the randomized-rank budget is planned from the start rather
    # than discovered in the SBC run.
    @test 0.03 <= atom_rate <= 0.12
    # And the atoms really are the ghat clamp's endpoints, not some other saturation.
    @test maximum(_P12P_F05.rho) <= GHAT_RHO_KNOTS[end]
    @test minimum(_P12P_F05.rho) >= GHAT_RHO_KNOTS[1]
end

@testset "r1 is drawn FIRST and the field is conditional on it" begin
    # (a) BEHAVIOURAL, NOT DOCUMENTARY: pinning r1 must actually move the field's spatial
    # statistics. A source-order check alone is weak -- it proves where a line sits, not what the
    # draw depends on.
    Zhi = _p12p_draws(200; r1 = 0.90).z
    Zlo = _p12p_draws(200; r1 = 0.05).z
    r1_hi = induced_lag1(Symmetric(cov(Zhi; dims = 2)))
    r1_lo = induced_lag1(Symmetric(cov(Zlo; dims = 2)))
    @test r1_hi > r1_lo                       # measured 0.8939 vs 0.0493
    @test r1_hi > 0.7
    @test r1_lo < 0.2
    @info "p12 empirical lattice lag-1 under a pinned r1" pinned_090 = r1_hi pinned_005 = r1_lo

    # (b) PINNING CONSUMES LESS OF THE STREAM THAN DRAWING DOES. If r1 came from a side channel
    # the two calls would consume identical amounts and the next draw off each stream would be
    # equal. It is not, so r1 rides the SAME rng that the field and the nuisances ride.
    ra = p12_fix_rng(P12_FIXTURE_COUNTER); sample_p12_prior(ra);            a = rand(ra)
    rb = p12_fix_rng(P12_FIXTURE_COUNTER); sample_p12_prior(rb; r1 = 0.5);  b = rand(rb)
    rc = p12_fix_rng(P12_FIXTURE_COUNTER); sample_p12_prior(rc);            c = rand(rc)
    @test a != b
    @test a == c                              # and the whole thing is deterministic (D-14)

    # (c) the drawn r1 lands inside the Tier-1 prior, and the arm is CARRIED so a persisted
    # sample records which prior produced it (T-12-25).
    d = sample_p12_prior(p12_fix_rng(P12_FIXTURE_COUNTER))
    @test P12_R1_MIN <= d.r1 <= P12_R1_MAX
    @test d.arm === :car
    @test sample_p12_prior(p12_fix_rng(P12_FIXTURE_COUNTER); arm = :gp).arm === :gp
    @test sample_p12_prior(p12_fix_rng(P12_FIXTURE_COUNTER); r1 = 0.33).r1 == 0.33
end

@testset "the ablation arm draws an independent field (D-10)" begin
    # THE ABLATION IS A GENUINE DRAW, not a special case bolted on at read time (R-4). It goes
    # through the same sampler, the same copula and the same datagen path as the spatial arms,
    # which is what makes D-13's descope deliverable real.
    #
    # M = 1 000 RATHER THAN THE PLAN'S 200, AND THE REASON IS RECORDED: the :none arm is the
    # identity covariance, so it costs no bisection and 1 000 draws take 0.05 s. At M = 200 the
    # screening W1 measures 0.0685 against a 0.10 bar -- a two-thirds-of-the-bar reading that is
    # Monte-Carlo noise, not signal. At M = 1 000 it is 0.0281. Raising N keeps the honest
    # pre-declared bar intact instead of loosening it (`test_simulator.jl:213-216`, same move).
    d = _p12p_draws(1_000; arm = :none)
    @test abs(induced_lag1(Symmetric(cov(d.z; dims = 2)))) < 0.05
    @test maximum(_p12p_per_region_w1(d.mu)) <= P12_SIM02_W1_TOL_PERREGION
    # The ablation r1 is OUTSIDE the prior by design, and any arm AT it is the same independent
    # field -- so a pool drawn with `arm = :car, r1 = P12_ABLATION_R1` is byte-for-byte the
    # ablation rather than a near-miss.
    @test P12_ABLATION_R1 < P12_R1_MIN
    @test lattice_sigma(:car, P12_ABLATION_R1) == lattice_sigma(:none, P12_ABLATION_R1)
end

@testset "the separable upsample equals a reference interpolant to tolerance (D-06)" begin
    rng = p12_fix_rng(P12_FIXTURE_COUNTER)
    F   = randn(rng, P12_G, P12_G) ./ 4
    # THE REFERENCE IS WRITTEN HERE, IN THE TEST, ON PURPOSE. It is the slow per-pixel
    # implementation the production path must agree with; having it in the test is what licenses
    # the fast path. The output size is deliberately non-square and not a multiple of G, so a
    # separability bug cannot hide behind a symmetric shape.
    function _ref_bilinear(F, m, n; offset = false)
        G1, G2 = size(F)
        xs = offset ? range(1 - 0.5, G1 + 0.5; length = m) : range(1, G1; length = m)
        ys = offset ? range(1 - 0.5, G2 + 0.5; length = n) : range(1, G2; length = n)
        out = Matrix{Float64}(undef, m, n)
        for (cc, y) in enumerate(ys), (rr, x) in enumerate(xs)
            i = clamp(floor(Int, x), 1, G1 - 1); t = x - i
            j = clamp(floor(Int, y), 1, G2 - 1); s = y - j
            out[rr, cc] = (1 - t) * (1 - s) * F[i, j]     + t * (1 - s) * F[i + 1, j] +
                          (1 - t) * s       * F[i, j + 1] + t * s       * F[i + 1, j + 1]
        end
        return out
    end
    fast = p12_upsample(F, (129, 97))
    @test size(fast) == (129, 97)
    @test maximum(abs, fast .- _ref_bilinear(F, 129, 97)) < 1e-12
    # The offset arm must agree with the SAME reference under the same keyword, or the guard and
    # the thing it guards have already diverged.
    @test maximum(abs, p12_upsample(F, (129, 97); offset = true) .-
                       _ref_bilinear(F, 129, 97; offset = true)) < 1e-12
    # Hoisted operators are the same operators (a datagen loop must be able to lift them out).
    A, B = p12_upsample_ops(P12_G, (129, 97))
    @test maximum(abs, (A * F * B') .- fast) == 0.0
end

@testset "the offset grid is a keyword, and it actually moves the lattice (D-06 guard)" begin
    rng = p12_fix_rng(P12_FIXTURE_COUNTER)
    F   = randn(rng, P12_G, P12_G) ./ 4
    @test p12_upsample(F, (256, 256); offset = true) != p12_upsample(F, (256, 256); offset = false)
    # On a CONSTANT field both are the constant: the offset changes only the sampling POSITIONS,
    # never the field's level. (The offset rows carry weights (1.5, -0.5) at the two half-cell
    # margins -- linear extrapolation -- and still sum to 1, which is why this holds.)
    C = fill(0.42, P12_G, P12_G)
    @test maximum(abs, p12_upsample(C, (256, 256); offset = true)  .- 0.42) < 1e-12
    @test maximum(abs, p12_upsample(C, (256, 256); offset = false) .- 0.42) < 1e-12
    @test all(x -> isapprox(x, 1.0; atol = 1e-12), sum(p12_bilinear_op(P12_G, 256; offset = true); dims = 2))
    # ONE KEYWORD, NOT A SECOND CODE PATH: the guard arm is reachable from the production
    # function, so it cannot drift from the thing it guards.
    @test p12_bilinear_op(P12_G, 64; offset = true) != p12_bilinear_op(P12_G, 64; offset = false)
end

@testset "a constant field renders to a constant image, up to weight rounding" begin
    C = fill(0.42, P12_G, P12_G)
    @test maximum(abs, p12_upsample(C, (512, 512)) .- 0.42) < 1e-12
    # WHY THIS IS A TOLERANCE AND NOT `==`. The row weights are `1-t` and `t`, whose
    # floating-point sum is 1 only up to rounding, so `A * F * B'` on a constant field is exactly
    # 0.42 in real arithmetic and within a few ULP in IEEE. This is precisely why the
    # exact-equality golden criterion in 12-08 is stated on the LEGACY SCALAR-theta path, which
    # does not go through the upsample at all -- and why p12_consts.jl:710-714 carries
    # P12_STAGE1_CONSTFIELD_TOL as a CRITERION rather than an excuse.
    @test P12_STAGE1_EXACT === true
    @test P12_STAGE1_CONSTFIELD_TOL == 1e-12
end

@testset "the field validates itself (ASVS V5)" begin
    @test_throws ArgumentError validate_rho_field([2.0 0.0; 0.0 0.0])
    @test_throws ArgumentError validate_rho_field([NaN 0.0; 0.0 0.0])
    @test_throws ArgumentError validate_rho_field(zeros(2, 3))
    # ... and a legitimate field passes through unchanged, so the guard is not simply always-throw.
    ok = [0.5 -0.5; 1.0 -1.0]
    @test validate_rho_field(ok) === ok
    # The guard is called from INSIDE rho_field, which is what makes it unavoidable rather than
    # advisory: `simulate_pair`'s scalar `-1 <= rho_true <= 1` check cannot see inside a field.
    @test all(-1.0 .<= _P12P_F05.rho .<= 1.0)
end

@testset "the derived scalar view is labelled derived and is a valid simulate_pair theta (R-1)" begin
    draw = sample_p12_prior(p12_fix_rng(P12_FIXTURE_COUNTER))
    θs   = theta_scalar_view(draw)
    @test length(θs) == 8
    @test propertynames(θs)[1] === :ρ_true
    @test last(propertynames(θs)) === :chromatic_eps
    @test -1 <= θs.ρ_true <= 1
    @test θs.chromatic_eps == draw.chromatic_eps      # the nuisances are carried, not redrawn
    @test θs.shift_dx == draw.shift_dx
    # The legacy path still accepts it UNCHANGED -- `forward.jl` is byte-untouched by this plan.
    ch = simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER), θs; imsize = (128, 128))
    @test length(ch) == 2
    @test all(c -> size(c) == (128, 128), ch)
    @test all(c -> all(isfinite, c), ch)
    @test all(c -> all(>=(0), c), ch)
end

@testset "prior surface ran CPU-only" begin
    @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
end

# =====================================================================================
# THE D-06 STAGE-1 GOLDEN REGRESSION (12-08) -- testsets 12-16, APPENDED below the eleven
# copula testsets above, which this plan does not touch.
# =====================================================================================
# EXACT EQUALITY IS THE PRE-REGISTERED CRITERION FOR THE NO-FIELD ARM AND IT IS NOT NEGOTIABLE.
# `P12_STAGE1_EXACT` IS `true` IN THE FROZEN TIER-1 PRE-REGISTRATION, SO TESTSET 13 COMPARES RAW
# Float64 BYTES WITH `==`. It holds BY CONSTRUCTION rather than by luck: a theta carrying no
# `rho_field` takes stage 1's SCALAR branch, which never enters the upsample at all, and the three
# `_smooth_field` calls sit BEFORE the branch so RNG consumption is identical on both paths. A
# FAILURE THERE IS A GENUINE SIGNAL -- an upstream numerics change, or a moved line in stage 1 --
# AND MUST BE INVESTIGATED, NOT SMOOTHED AWAY. THE ONLY SANCTIONED ESCAPE HATCH IS THE
# PRE-REGISTERED `P12_STAGE1_GOLDEN_TOLERANCE_FALLBACK` (p12_consts.jl §17), AND ANY USE OF IT MUST
# BE RECORDED AS A NAMED DEVIATION IN THE PLAN SUMMARY. DO NOT SOFTEN A FAILING ASSERTION TO A
# TOLERANCE COMPARISON IN THIS FILE.
#
# TWO MEASURED CORRECTIONS THIS SECTION CARRIES, BOTH BECAUSE THE ASSERTION AS SPECIFIED WOULD
# OTHERWISE BE UNABLE TO PASS OR UNABLE TO FAIL.
#
# (1) THE CONSTANT-FIELD DIFFERENCE IS NOT "A FEW ULP". IT IS EXACTLY ZERO. 12-VALIDATION.md, and
# `p12_consts.jl:710-714` in its wake, reason that a constant-VALUED field must differ from the
# scalar path in the last few ULP, because the separable upsample's row weights `1-t` and `t` sum
# to 1 only up to IEEE rounding. MEASURED THIS SESSION on the committed operator, that rounding
# does not materialize at these operating points: `_p12_upsample_local(fill(0.7, 8, 8), sz)` is
# EXACTLY 0.7 in every pixel at 128^2, 256^2 and 512^2 (max deviation 0.0, all entries `==`), so
# `a`, `b` and `sgn` are bit-identical to the scalar branch's and the whole seven-stage pipeline
# reproduces the golden bytes with max|delta| = 0.0 across all six keys. 12-07's testset 8 measured
# the same thing independently at 0.42 / 512^2. The pre-registered criterion
# `max|delta| <= P12_STAGE1_CONSTFIELD_TOL` therefore HOLDS (0.0 <= 1e-12) and is asserted
# unchanged -- but the falsifier the plan pairs with it, `0 < max|delta|`, is measurably FALSE, and
# its stated reasoning ("a difference of precisely zero proves the branch did NOT run") is false in
# the same breath: the branch DID run and the difference IS precisely zero. Testset 14 therefore
# proves the branch ran the way that actually bites -- a constant field at a value DIFFERENT from
# `theta.rho_true` must OVERRIDE it. If the field branch were skipped that theta would return the
# golden bytes and the difference would collapse to 0; measured, it is 3.04.
#
# (2) THE TWO VALIDATORS DIVERGE ON EXACTLY ONE OF THE FIVE SHARED BAD INPUTS, AND THE DIVERGENCE
# IS RECORDED RATHER THAN ASSERTED AWAY. Measured: on `> 1`, `< -1`, `NaN` and a non-square matrix
# both `validate_rho_field` and a field-carrying `simulate_pair` throw `ArgumentError`. On a 1x1
# matrix `validate_rho_field` does NOT throw (1x1 is square, finite and in range), while
# `simulate_pair` does. That is not a bug in either: the "at least two cells per axis" rule is a
# constraint of the RENDERER, and on the `p12_prior.jl` side it lives in `p12_bilinear_op`'s
# `G >= 2` guard rather than in the validator -- `p12_upsample(fill(0.5,1,1), ...)` throws too.
# `forward.jl`'s inlined `_p12_check_field` merges the two into one entry guard because it has no
# separate upsample to fail in. Testset 16 asserts the four-way agreement AND this one divergence
# explicitly, so a future edit that silently changes either side is caught.
#
# CPU-only. Two `simulate_pair` calls at 128^2 x 6 keys plus four at 256^2 (testset 15).

using JLD2

# `contract.jl` is reached for the SIM-03 patch summary (testset 15) and is GUARDED: inside the
# aggregated suite an earlier file has already loaded it, and its two `src/` includes are READ-ONLY
# (CLAUDE.md decoupling). Measured cost when it is not already loaded: 0.84 s.
isdefined(@__MODULE__, :build_mci) || include(joinpath(@__DIR__, "..", "contract.jl"))

# PHASE-PREFIXED, deliberately: a bare `_res`-style helper in this repository has already silently
# overwritten a Phase-13 namesake through an identical zero-positional signature.
const _P12G_PATH   = joinpath(@__DIR__, "fixtures", "p12_stage1_golden.jld2")
const _P12G        = JLD2.load(_P12G_PATH)
const _P12G_THETA  = _P12G["theta"]
const _P12G_IMSIZE = _P12G["imsize"]
const _P12G_CHANS  = _P12G["channels"]
const _P12G_KEYS   = _P12G["golden_keys"]
# The off-value constant field of testset 14's falsifier: a legal constant field that is NOT
# `theta.rho_true`, so "the field was read" and "the field was ignored" have different answers.
const _P12G_OFFVAL = 0.1

"""
    _p12g_source_between(startmark, endmark) -> String

The lines of THIS file strictly between two sentinel comments, with whole-line comments stripped --
the `test_stage6_regression.jl:79-80` idiom, applied to this file's own source so a claim about how
testset 13 is written is checkable from the running suite instead of by review.
"""
function _p12g_source_between(startmark::AbstractString, endmark::AbstractString)
    lines = split(replace(read(@__FILE__, String), "\r\n" => "\n"), '\n')
    i = findfirst(l -> occursin(startmark, l), lines)
    j = findfirst(l -> occursin(endmark, l), lines)
    (i === nothing || j === nothing || j <= i + 1) && return ""
    return join(filter(l -> !startswith(strip(l), "#"), lines[(i + 1):(j - 1)]), '\n')
end

"Across-region spread of a patch-correlation grid, `missing` patches dropped."
_p12g_spread(M) = std(collect(skipmissing(vec(M))))

@testset "the golden fixture PRECEDES the stage-1 edit (D-06)" begin
    # A fixture REGENERATED after the edit would have been captured through a Phase-12 draw, would
    # carry `rho_field`, and would compare the new code against itself -- passing vacuously and
    # silently. Refuse that fixture outright, exactly as test_stage6_regression.jl:129-134 refuses
    # a post-eps fixture by the presence of the 8th field.
    @test !hasproperty(_P12G_THETA, :rho_field)
    @test length(_P12G_THETA) == 8
    @test last(keys(_P12G_THETA)) === :chromatic_eps
    @test _P12G_THETA.chromatic_eps == 0.0        # the stage-6 term is pinned INERT by the literal

    sha = _P12G["p12_base_sha"]
    @test sha isa AbstractString
    @test length(sha) == 40

    # The fixture rode the FIXTURE stream on the FIXTURE counter, never a reported one.
    @test _P12G["fixture_counter"] == P12_FIXTURE_COUNTER
    @test _P12G["seed"] == UInt64(P12_FIX_SEED)
    @test _P12G["salt"] == UInt64(P12_SALT)
    @test length(_P12G_CHANS) == length(_P12G_KEYS) == 6

    # AND THE ORDERING IS CHECKABLE, NOT MERELY CLAIMED: the capture-time sha must be an ANCESTOR
    # of HEAD. Guarded by the test_p13_result.jl:199-201 predicate so a sandbox without git skips
    # the executable half rather than failing for an unrelated reason.
    _has_git = Sys.which("git") !== nothing &&
               (isdir(joinpath(P12_REPO_ROOT, ".git")) || isfile(joinpath(P12_REPO_ROOT, ".git")))
    if _has_git
        @test success(Cmd(`git merge-base --is-ancestor $sha HEAD`; dir = P12_REPO_ROOT))
    else
        @info "git unavailable — skipping the executable fixture-precedes-edit ancestry assertion"
    end
end

@testset "a theta with no rho_field reproduces the golden bytes EXACTLY (P12_STAGE1_EXACT)" begin
    @test P12_STAGE1_EXACT == true      # the contract this testset scores against, asserted where
                                        # it is enforced rather than only where it is declared
    # --- P12-T13-EXACT-ARM-BEGIN ---
    for (i, k) in enumerate(_P12G_KEYS)
        out = simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER + k), _P12G_THETA;
                            imsize = _P12G_IMSIZE)
        @test collect(out[1]) == collect(_P12G_CHANS[i][1])
        @test collect(out[2]) == collect(_P12G_CHANS[i][2])
    end
    # --- P12-T13-EXACT-ARM-END ---

    # SOURCE-LEVEL TRIPWIRE, mirroring the `count("warp(", code) == 1` assertion
    # test_stage6_regression.jl:94-96 makes on the simulator: the comparison above must stay an
    # exact one. A future edit that softens it to a tolerance changes the SOURCE, so the check is
    # made on the source. The strings below sit OUTSIDE the sentinels on purpose -- inside them
    # they would match themselves and the tripwire would fail on its own text.
    exact_arm = _p12g_source_between("P12-T13-EXACT-ARM-BEGIN", "P12-T13-EXACT-ARM-END")
    @test !isempty(exact_arm)                     # the extraction really found the block ...
    @test occursin("_P12G_CHANS", exact_arm)      # ... and the block really is the comparison
    @test occursin("==", exact_arm)
    @test !occursin("isapprox", exact_arm)
    @test !occursin("≈", exact_arm)
    @test !occursin("atol", exact_arm)
    @test !occursin("rtol", exact_arm)
end

@testset "a CONSTANT field agrees with the scalar path to the pre-registered tolerance" begin
    ρ  = _P12G_THETA.ρ_true
    θf = merge(_P12G_THETA, (rho_field = fill(ρ, 8, 8),))
    @test length(θf) == 9 && hasproperty(θf, :rho_field)

    dmax = 0.0
    for (i, k) in enumerate(_P12G_KEYS)
        out = simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER + k), θf; imsize = _P12G_IMSIZE)
        dmax = max(dmax,
                   maximum(abs, collect(out[1]) .- collect(_P12G_CHANS[i][1])),
                   maximum(abs, collect(out[2]) .- collect(_P12G_CHANS[i][2])))
    end
    # THE PRE-REGISTERED CRITERION, UNCHANGED AND NOT MOVED. It is a TOLERANCE by construction
    # rather than a softened `==`, because the upsample's row weights `1-t` and `t` sum to 1 only up
    # to IEEE rounding -- see p12_consts.jl:710-714.
    @test dmax <= P12_STAGE1_CONSTFIELD_TOL
    # MEASURED THIS SESSION: the rounding does not materialize here and dmax is EXACTLY 0.0 (header
    # correction 1). Recorded with @info so the phase report can quote the number rather than the
    # bound, and asserted only as the bound, which is what was pre-registered.
    @info "p12 stage-1 constant-field agreement with the scalar path" max_abs_delta = dmax tol = P12_STAGE1_CONSTFIELD_TOL exactly_zero = (dmax == 0.0)

    # THE FALSIFIER, AND IT IS NOT `0 < dmax` -- THAT ONE IS MEASURABLY FALSE (header correction 1).
    # A constant field at a value DIFFERENT from theta.rho_true must OVERRIDE the scalar. If stage
    # 1's field branch never ran, this theta would take the rho_true = 0.7 path and return the
    # golden bytes, collapsing `off_vs_golden` to 0.0 and failing the bound below. Measured: 3.04.
    θoff = merge(_P12G_THETA, (rho_field = fill(_P12G_OFFVAL, 8, 8),))
    o_off = simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER + _P12G_KEYS[1]), θoff;
                          imsize = _P12G_IMSIZE)
    off_vs_golden = max(maximum(abs, collect(o_off[1]) .- collect(_P12G_CHANS[1][1])),
                        maximum(abs, collect(o_off[2]) .- collect(_P12G_CHANS[1][2])))
    @test off_vs_golden > 0.5
    @test _P12G_OFFVAL != ρ            # ... and the two values really are different, so the
                                       # falsifier cannot pass by comparing rho_true with itself

    # ... and the same constant field reproduces the SCALAR path at ITS OWN value to the same
    # pre-registered bar, which is the claim "a constant field is the scalar case" in full.
    o_scal = simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER + _P12G_KEYS[1]),
                           merge(_P12G_THETA, (ρ_true = _P12G_OFFVAL,)); imsize = _P12G_IMSIZE)
    off_vs_scalar = max(maximum(abs, collect(o_off[1]) .- collect(o_scal[1])),
                        maximum(abs, collect(o_off[2]) .- collect(o_scal[2])))
    @test off_vs_scalar <= P12_STAGE1_CONSTFIELD_TOL
    @info "p12 stage-1 constant field at an OFF value" off_value = _P12G_OFFVAL vs_golden = off_vs_golden vs_scalar_at_same_value = off_vs_scalar
end

@testset "a NON-constant field changes the image, and changes it spatially" begin
    # A REAL draw from the D-05 copula, not a hand-built field, so what reaches stage 1 is what
    # 12-09's datagen will send it.
    d  = sample_p12_prior(p12_fix_rng(P12_FIXTURE_COUNTER); r1 = 0.5)
    θF = merge(_P12G_THETA, (rho_field = d.rho_field,))
    @test size(θF.rho_field) == (P12_G, P12_G)
    @test !all(==(first(d.rho_field)), d.rho_field)     # it really is non-constant

    # 256^2 so each of the 64 patches is 32x32 px, far above the >15-surviving-pixel floor.
    for k in _P12G_KEYS[1:2]
        ctrl = simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER + k), _P12G_THETA; imsize = (256, 256))
        fld  = simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER + k), θF;          imsize = (256, 256))
        # (a) the field changes the image at all
        @test collect(fld[1]) != collect(ctrl[1])
        @test collect(fld[2]) != collect(ctrl[2])

        # (b) AND IT REACHES THE PATCH SUMMARY rather than being averaged away by the PSF, which is
        # the only part of (a) that is actually about D-06 delivering a SPATIAL map. Directional
        # check, NOT the D-12 gate: no bar here is pre-registered, and none is invented -- the
        # comparison is against the constant-rho control measured on the same key.
        Sc = patch_summary(build_mci(ctrl))
        Sf = patch_summary(build_mci(fld))
        sc, sf = _p12g_spread(Sc), _p12g_spread(Sf)
        @test sf > sc
        # ... and the field arm reaches regions the constant-rho control never does. At
        # rho_true = 0.7 the control's grid is entirely positive; a drawn field spans the prior, so
        # its weakest region sits below the control's. This is the assertion that separates
        # "spatially varying" from "noisier".
        @test minimum(skipmissing(vec(Sf))) < minimum(skipmissing(vec(Sc)))
        @info "p12 stage-1 across-region spread, field vs constant-rho control" key = k control_spread = sc field_spread = sf ratio = sf / sc
    end
    # 12-RESEARCH Pitfall 2 measured a constant-rho across-region spread of 0.0847 against a
    # per-draw noise sd of 0.0662 and warned that only a modest excess is expected. That pairing is
    # NOT this configuration and is not quoted as if it were: measured here at 256^2, rho = 0.7 and
    # this theta's nuisances, the control spread is ~0.155-0.22 and the field arm runs ~2.2x it.
end

@testset "the two field validators agree (the deliberate duplication is bounded)" begin
    _p12g_bad = (("an entry > 1",  [2.0 0.0; 0.0 0.0]),
                 ("an entry < -1", [-2.0 0.0; 0.0 0.0]),
                 ("a NaN",         [NaN 0.0; 0.0 0.0]),
                 ("a non-square",  zeros(2, 3)))
    # THE FOUR-WAY AGREEMENT. This is the test `forward.jl`'s `_p12_check_field` docstring names as
    # the thing keeping the inlined guard honest with `p12_prior.jl`'s `validate_rho_field`: the
    # duplication exists because forward.jl must not include p12_prior.jl (that would invert the
    # dependency and pull the lattice into every consumer of the forward model), so the two are
    # pinned to one another here instead.
    for (name, F) in _p12g_bad
        @test_throws ArgumentError validate_rho_field(F)
        @test_throws ArgumentError simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER),
                                                 merge(_P12G_THETA, (rho_field = F,));
                                                 imsize = _P12G_IMSIZE)
    end

    # THE ONE DIVERGENCE, ASSERTED RATHER THAN ASSUMED AWAY (header correction 2). A 1x1 field is
    # square, finite and in range, so `validate_rho_field` accepts it; `simulate_pair` rejects it
    # because a bilinear upsample interpolates BETWEEN cell pairs and there is no pair. On the
    # p12_prior.jl side that rule lives in the RENDERER (`p12_bilinear_op`'s `G >= 2`), not in the
    # validator -- so nothing is missing there either, and the two sides together reject exactly
    # what forward.jl's merged guard rejects.
    _p12g_tiny = fill(0.5, 1, 1)
    @test validate_rho_field(_p12g_tiny) === _p12g_tiny            # accepts -- by design
    @test_throws ArgumentError p12_upsample(_p12g_tiny, (128, 128))  # the renderer refuses
    @test_throws ArgumentError simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER),
                                             merge(_P12G_THETA, (rho_field = _p12g_tiny,));
                                             imsize = _P12G_IMSIZE)

    # A LEGITIMATE FIELD STILL PASSES, so neither guard is simply always-throw.
    @test simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER),
                        merge(_P12G_THETA, (rho_field = fill(0.3, 8, 8),));
                        imsize = _P12G_IMSIZE) isa Vector{Matrix{Float64}}

    # THE UPSAMPLE HALF OF THE SAME DUPLICATION. `_p12_upsample_local` is a private copy of
    # `p12_upsample`'s `offset = false` path; if the copy ever drifts, every field-carrying
    # simulation silently renders a different image from the one 12-09's pool believes it drew.
    _p12g_F = randn(p12_fix_rng(P12_FIXTURE_COUNTER), P12_G, P12_G) ./ 4
    @test _p12_upsample_local(_p12g_F, (128, 128)) == p12_upsample(_p12g_F, (128, 128))
    @test _p12_upsample_local(_p12g_F, (129, 97))  == p12_upsample(_p12g_F, (129, 97))

    # The golden-regression section ran CPU-only too (the eleven copula testsets assert this above,
    # but this section is the one that loads contract.jl and src/).
    @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
end
