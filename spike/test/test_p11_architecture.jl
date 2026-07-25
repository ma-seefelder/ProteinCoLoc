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

# spike/test/test_p11_architecture.jl --- unit gate for the Phase-11 research-net model surface.
#
# Covers the two deviations `spike/npe/p11_architecture.jl` implements: the ported F2 bounded
# theta transform (ruling Q4) and the 129th lambda conditioning row (D-03). Run:
#     julia --project=spike spike/test/test_p11_architecture.jl
#
# WHAT THIS FILE IS FOR. The port exists to remove a CONFOUND from the D-02 rho_true attenuation
# measurement, so "we ported it" is not enough -- it has to be shown to be the SAME transform the
# shipped net uses, and to actually deliver the F2 property (no mass outside the prior box). Both
# are asserted here against the real `src` implementation, read-only.
#
# CPU-only (D-10). No training, no simulation, no file write.

using Test
using Statistics
using StatsBase
using Random
using Random123
using Pkg

# The shipped model surface, loaded into an ISOLATED module and used READ-ONLY: it defines the
# same concepts under unprefixed names (`BoundedThetaTransform`, `THETA_LOGIT_EPS`, `NPE_D`) and
# we want ONLY to compare against them. Nothing is trained, written or mutated through it.
#
# THIS SITS OUTSIDE ANY TESTSET ON PURPOSE: Julia rejects a `module` expression that is not at
# top level (the same constraint `spike/validation/p11_consts.jl` documents for `module GateV2`).
# MEASURED: `src/amortized/architecture.jl` DOES load standalone under `--project=spike` -- its
# only imports are Flux, NeuralEstimators and StatsBase, all already spike dependencies. No
# dependency was added to make this work; if it ever stops loading, the correct response is to
# drop the cross-check, not to add a dependency (D-01).
module SrcArch
    include(joinpath(@__DIR__, "..", "..", "src", "amortized", "architecture.jl"))
end

# ORDER MATTERS: the Tier-1 pre-registration (LAMBDA_MIN/LAMBDA_MAX, P11_FIXTURE_COUNTER,
# P11_RESEARCH_NET_DEVIATIONS) before the unit under test. Guarded for idempotency.
isdefined(@__MODULE__, :LAMBDA_MIN)   || include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
isdefined(@__MODULE__, :encode_lambda) || include(joinpath(@__DIR__, "..", "npe", "p11_architecture.jl"))

# A deterministic synthetic theta matrix INSIDE the research prior box, drawn on the FIXTURE
# counter (99) so it can never consume a counter a reported number rides on (p11_consts.jl §2).
const P11A_BOUNDS = p11_theta_prior_bounds()
const P11A_LO     = Float64[b[1] for b in P11A_BOUNDS]
const P11A_HI     = Float64[b[2] for b in P11A_BOUNDS]

function _p11a_theta(n::Integer)
    rng = p11_rng(P11_FIXTURE_COUNTER)
    Θ = Matrix{Float64}(undef, length(P11A_BOUNDS), n)
    @inbounds for j in 1:n, i in eachindex(P11A_LO)
        Θ[i, j] = P11A_LO[i] + (P11A_HI[i] - P11A_LO[i]) * rand(rng)
    end
    return Θ
end

const P11A_THETA = _p11a_theta(200)

@testset "Phase 11 — research-net model surface (D-02, D-03, ruling Q4)" verbose = true begin

    @testset "ported transform agrees with src" begin
        # The clamp constant must be numerically identical, or the two transforms differ at the
        # atoms even when everything else matches.
        @test P11_THETA_LOGIT_EPS == SrcArch.THETA_LOGIT_EPS

        # Construct BOTH transforms from the SAME theta matrix and the SAME bounds, each through
        # its own forward map, then compare end to end.
        U_src = SrcArch.theta_to_unbounded(P11A_THETA, P11A_LO, P11A_HI)
        U_p11 = p11_theta_to_unbounded(P11A_THETA, P11A_LO, P11A_HI)
        @test maximum(abs.(U_src .- U_p11)) < 1e-12

        zt_src = StatsBase.fit(StatsBase.ZScoreTransform, U_src; dims = 2)
        for i in eachindex(zt_src.scale)
            (isfinite(zt_src.scale[i]) && zt_src.scale[i] > 0) || (zt_src.scale[i] = 1.0)
        end
        t_src = SrcArch.BoundedThetaTransform(P11A_LO, P11A_HI, zt_src)
        t_p11 = fit_p11_theta_transform(P11A_THETA; bounds = P11A_BOUNDS)

        @test t_p11 isa P11BoundedThetaTransform
        @test t_p11.lo == P11A_LO && t_p11.hi == P11A_HI
        @test maximum(abs.(t_p11.zt.mean  .- t_src.zt.mean))  < 1e-12
        @test maximum(abs.(t_p11.zt.scale .- t_src.zt.scale)) < 1e-12

        # Forward agreement.
        Y_src = StatsBase.transform(t_src, P11A_THETA)
        Y_p11 = StatsBase.transform(t_p11, P11A_THETA)
        @test maximum(abs.(Y_src .- Y_p11)) < 1e-12

        # `reconstruct` agreement, on flow-space values the flow could plausibly emit.
        W = 6.0 .* randn(Random.Xoshiro(20261125), length(P11A_BOUNDS), 300)
        @test maximum(abs.(StatsBase.reconstruct(t_src, W) .-
                           StatsBase.reconstruct(t_p11, W))) < 1e-12

        # Round-trip: theta -> flow space -> theta.
        @test maximum(abs.(StatsBase.reconstruct(t_p11, Y_p11) .- P11A_THETA)) < 1e-4

        # Vector methods dispatch and agree with the single-column matrix form.
        v = P11A_THETA[:, 1]
        @test StatsBase.transform(t_p11, v) ≈ vec(StatsBase.transform(t_p11, reshape(v, :, 1)))
        @test StatsBase.reconstruct(t_p11, Y_p11[:, 1]) ≈
              vec(StatsBase.reconstruct(t_p11, Y_p11[:, 1:1]))
    end

    @testset "bounded transform keeps mass inside the prior box" begin
        # THE F2 PROPERTY THE PORT EXISTS TO BUY. Asserted, never assumed.
        t = fit_p11_theta_transform(P11A_THETA)
        wild = 50.0 .* randn(Random.Xoshiro(20261126), length(P11A_BOUNDS), 200)
        back = StatsBase.reconstruct(t, wild)

        # Containment is CLOSED-interval by construction: `p11_theta_from_unbounded` ends in a
        # `clamp` into [lo, hi], and at |y| large enough `logistic` underflows to exactly 0 or 1,
        # so an extreme flow-space value legitimately lands ON an endpoint. This mirrors the
        # assertion the src twin makes (test/runtests.jl:563) -- the guarantee is "no draw leaves
        # the support", not "no draw touches the boundary".
        @test all(P11A_LO[p] <= back[p, j] <= P11A_HI[p]
                  for p in axes(back, 1), j in axes(back, 2))

        # STRICT interior on the range the flow actually works in. A DETERMINISTIC grid, not a
        # random draw: at very large |y| `logistic` underflows and `lo + width*logistic(y)`
        # rounds to exactly `lo` in Float64 (the chromatic row's width is only 0.04), so a
        # random-tail strict assertion would be flaky for a purely floating-point reason. The
        # substantive claim -- values in the flow's working range land strictly inside -- is
        # asserted on a fixed grid where it is exactly reproducible.
        mid   = repeat(collect(-4.0:0.5:4.0)', length(P11A_BOUNDS), 1)
        backm = StatsBase.reconstruct(t, mid)
        @test all(P11A_LO[p] < backm[p, j] < P11A_HI[p]
                  for p in axes(backm, 1), j in axes(backm, 2))

        # Float32 draws (what `sampleposterior` returns) keep their element type.
        @test eltype(StatsBase.reconstruct(t, Float32.(wild))) == Float32

        # Strictly monotone per parameter => SBC ranks are invariant to the change of space.
        mono = StatsBase.reconstruct(t, repeat(collect(-6.0:0.5:6.0)', length(P11A_BOUNDS), 1))
        @test all(issorted(mono[p, :]) for p in axes(mono, 1))

        # The box is DERIVED, never retyped: 8 rows, in sample_prior field order, each row the
        # prior object's own support. Row 5/6 DELIBERATELY differ from src (D-02, spike-only).
        @test length(P11A_BOUNDS) == 8
        @test P11A_BOUNDS[1] == (first(GHAT_RHO_KNOTS), last(GHAT_RHO_KNOTS))
        @test P11A_BOUNDS[5] == (minimum(SHIFT_PRIOR), maximum(SHIFT_PRIOR))
        @test P11A_BOUNDS[5] == (-3.0, 3.0)
        @test P11A_BOUNDS[8] == (minimum(CHROMATIC_PRIOR), maximum(CHROMATIC_PRIOR))
        @test all(lo < hi for (lo, hi) in P11A_BOUNDS)
    end

    @testset "encode_lambda" begin
        @test encode_lambda(LAMBDA_MIN) == 0.0
        @test encode_lambda(LAMBDA_MAX) == 1.0
        # Strictly increasing on the pre-registered SC2 rung ladder.
        enc = [encode_lambda(r) for r in SC2_RUNGS]
        @test all(enc[i] < enc[i + 1] for i in 1:(length(enc) - 1))
        @test all(0.0 .<= enc .<= 1.0)
        # Invertible: the map is affine and data-independent, so the inverse is exact.
        for lam in SC2_RUNGS
            @test decode_lambda(encode_lambda(lam)) ≈ lam
        end
        @test decode_lambda(0.0) == LAMBDA_MIN
        @test decode_lambda(1.0) == LAMBDA_MAX
        # NOT clamped outside the trained range: extrapolation must stay visible to the caller.
        @test encode_lambda(LAMBDA_MAX + 1.0) > 1.0
        @test encode_lambda(0.0) < 0.0
    end

    @testset "augment_input ORDER invariant" begin
        Z = zeros(128, 5)
        A = augment_input(Z, 1.0)
        @test size(A) == (129, 5)
        @test eltype(A) == Float32
        # The lambda row is CONSTANT across columns and carries the encoded value.
        @test all(A[129, j] == Float32(encode_lambda(1.0)) for j in axes(A, 2))
        @test length(unique(A[129, :])) == 1
        # The first 128 rows pass through untouched.
        @test A[1:128, :] == Float32.(Z)

        # THE GUARD (T-11-15): a 129-row input throws, so appending twice -- or appending BEFORE
        # standardization and then "fixing" the partition -- is caught at runtime, not silently
        # z-scored against the training pool's lambda distribution (Pitfall 2).
        @test_throws AssertionError augment_input(zeros(129, 5), 1.0)
        @test_throws AssertionError augment_input(zeros(64, 5), 1.0)
        @test_throws AssertionError augment_input(zeros(142, 5), 1.0)

        # The vector form is the single-column matrix form.
        @test size(augment_input(zeros(128), 2.0)) == (129, 1)
        @test augment_input(zeros(128), 2.0)[129, 1] == Float32(encode_lambda(2.0))

        # Distinct lambdas produce distinct conditioning rows (the row carries information).
        @test augment_input(Z, LAMBDA_MIN)[129, 1] != augment_input(Z, LAMBDA_MAX)[129, 1]
    end

    @testset "build_p11_estimator constructs at d_in = 129, D = 8" begin
        est = build_p11_estimator()
        @test est !== nothing
        # The capacity knobs are the unchanged Phase-5-calibrated defaults, so the comparison to
        # the shipped net is capacity-controlled.
        @test NPE_DSTAR == 64 && NPE_DEPTH == 3 && NPE_WIDTH == 256
        @test NPE_COUPLING == 10 && NPE_FLOW_DEPTH == 2 && NPE_FLOW_WIDTH == 128
        @test NPE_DSTAR >= 8                    # the `dstar >= D` guard, at D = 8
        # D-16: the SHIPPED read surface is untouched. `NPE_D` stays 7 in BOTH copies, and the
        # research D = 8 is passed explicitly rather than by moving the default.
        @test NPE_D == 7
        @test SrcArch.NPE_D == 7
        # An estimator at the shipped width/marginal count still builds from the same builder --
        # proof the wrapper added an argument, not a signature change.
        @test build_estimator(128) !== nothing
        # Guards still fire on nonsense.
        @test_throws ArgumentError build_p11_estimator(d_in = 0)
        @test_throws ArgumentError build_p11_estimator(D = 0)
        @test_throws ArgumentError build_p11_estimator(D = 65)   # dstar (64) >= D violated
    end

    @testset "the four declared deviations are the pre-registered four (D-04)" begin
        # The file's DECLARED DEVIATIONS header must not drift from the Tier-1 enumeration.
        @test length(P11_RESEARCH_NET_DEVIATIONS) == 4
        # Booleans are computed FIRST so a failure prints `false`, not the whole source file.
        srctext         = read(joinpath(@__DIR__, "..", "npe", "p11_architecture.jl"), String)
        has_declared    = occursin("DECLARED DEVIATIONS", srctext)
        has_attribution = occursin("src/amortized/architecture.jl", srctext)
        has_not_dropin  = occursin("NOT DROP-IN", srctext)          # D-03/D-16 consequence stated
        has_read_surf   = occursin("SHIPPED READ SURFACE", srctext)
        # No bare single-letter epsilon identifier was introduced (the chromatic-vs-logit name
        # collision the PATTERNS map warns about).
        has_bare_eps    = occursin(r"(^|[^_[:alnum:]])(eps|ε)([^_[:alnum:]]|$)", srctext)
        @test has_declared
        @test has_attribution
        @test has_not_dropin
        @test has_read_surf
        @test !has_bare_eps
    end

    @testset "model surface ran CPU-only (D-10)" begin
        @test !haskey(Pkg.project().dependencies, "CUDA")
        @test !any(p -> occursin("CUDA", p.name), values(Pkg.dependencies()))
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
