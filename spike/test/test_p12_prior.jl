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
