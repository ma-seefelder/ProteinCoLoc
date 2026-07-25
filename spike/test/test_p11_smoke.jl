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

# spike/test/test_p11_smoke.jl --- SC1f: d_in = 129, D = 8 trains (retires R10).
#
# THIS FILE MUST PASS BEFORE THE ROUGHLY 1.2 HOUR PHASE-11 DATAGEN IS LAUNCHED.
# Risk R10 is that the NeuralEstimators v0.2.1 conditioning-input API cannot express D-03 --
# that a 129-row input and an 8-marginal flow simply do not go through `train`. Discovering
# that AFTER datagen costs 1.2 h of simulation plus the training attempt; discovering it here
# costs seconds, because the whole exercise runs on SYNTHETIC in-memory data. No simulator is
# called, no cache is read, no artifact is written. The point is that the SHAPES AND THE API
# WORK -- NOT that the net converges. A 3-epoch budget cannot and does not claim otherwise, and
# nothing in this file may be quoted as a calibration, width or accuracy result.
#
# CPU-only (D-10). Run:
#     julia --project=spike spike/test/test_p11_smoke.jl
#
# DELIBERATELY NOT WIRED INTO spike/test/runtests.jl -- the Phase-11 files are run directly per
# the validation contract, and several of the sibling checks need a trained research net that
# does not exist yet.

using Test
using Statistics
using StatsBase
using Random
using Random123
using Pkg
using Flux                # Flux.Optimisers.AdamW (the Float64 optimiser-arg gotcha below)
using NeuralEstimators    # train

# ORDER MATTERS: the Tier-1 pre-registration (LAMBDA_MIN/LAMBDA_MAX, P11_FIXTURE_COUNTER), then
# the research model surface (augment_input / fit_p11_theta_transform / build_p11_estimator),
# then the spike read surface (posterior_for). Guarded for idempotency.
isdefined(@__MODULE__, :LAMBDA_MIN)    || include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
isdefined(@__MODULE__, :encode_lambda) || include(joinpath(@__DIR__, "..", "npe", "p11_architecture.jl"))
isdefined(@__MODULE__, :posterior_for) || include(joinpath(@__DIR__, "..", "npe", "infer.jl"))

# --- Smoke resolution knobs. NOT pre-registered thresholds: these are the smallest numbers that
#     still exercise every shape (a batch boundary, a train/val split, a multi-epoch loop).
const SMOKE_N_TRAIN  = 200
const SMOKE_N_VAL    = 50
const SMOKE_EPOCHS   = 3
const SMOKE_BATCH    = 32
const SMOKE_N_DRAWS  = 64

const SMOKE_BOUNDS = p11_theta_prior_bounds()
const SMOKE_LO     = Float64[b[1] for b in SMOKE_BOUNDS]
const SMOKE_HI     = Float64[b[2] for b in SMOKE_BOUNDS]

"""
    _smoke_split(rng, n) -> NamedTuple

Synthesize `n` (theta, 129-row input) pairs entirely in memory.

The lambda row is built PER SAMPLE and the samples' lambdas DIFFER, which is the whole point: a
single scalar lambda broadcast across the matrix would leave the 129th row constant and would
therefore not exercise the conditioning path at all -- the very failure mode R4 describes. The
shift rows are drawn CONDITIONAL on that sample's lambda (`Uniform(-lam, lam)`, D-13), so the
smoke rehearses the real generative story's shape rather than an independent draw.

`Z` is `randn` noise standing in for a standardized summary; the smoke asserts API and shape
conformance, so the summary's CONTENT is irrelevant here and deliberately not simulated.
"""
function _smoke_split(rng, n::Integer)
    lam = [LAMBDA_MIN + (LAMBDA_MAX - LAMBDA_MIN) * rand(rng) for _ in 1:n]
    Θ   = Matrix{Float64}(undef, length(SMOKE_BOUNDS), n)
    cols = Vector{Matrix{Float32}}(undef, n)
    @inbounds for j in 1:n
        for i in eachindex(SMOKE_LO)
            Θ[i, j] = SMOKE_LO[i] + (SMOKE_HI[i] - SMOKE_LO[i]) * rand(rng)
        end
        # Rows 5/6 are shift_dx/shift_dy: drawn CONDITIONAL on this sample's lambda (D-13).
        Θ[5, j] = -lam[j] + 2 * lam[j] * rand(rng)
        Θ[6, j] = -lam[j] + 2 * lam[j] * rand(rng)
        cols[j] = augment_input(randn(rng, Float32, 128, 1), lam[j])
    end
    return (Θ = Θ, Z = hcat(cols...), lam = lam)
end

const SMOKE_RNG = p11_rng(P11_FIXTURE_COUNTER)   # FIXTURE counter: never a reported stream
const SMOKE_TR  = _smoke_split(SMOKE_RNG, SMOKE_N_TRAIN)
const SMOKE_VA  = _smoke_split(SMOKE_RNG, SMOKE_N_VAL)

# Leak-free theta standardization: the ported transform is fitted on the TRAIN theta ONLY and
# then applied to both splits.
const SMOKE_ΘZT     = fit_p11_theta_transform(SMOKE_TR.Θ)
const SMOKE_ΘTR_STD = Float32.(StatsBase.transform(SMOKE_ΘZT, SMOKE_TR.Θ))
const SMOKE_ΘVA_STD = Float32.(StatsBase.transform(SMOKE_ΘZT, SMOKE_VA.Θ))

println("="^78)
println("SC1f smoke — d_in = 129, D = 8 (retires R10). SYNTHETIC data, $(SMOKE_EPOCHS) epochs.")
println("  train $(size(SMOKE_TR.Z)) / val $(size(SMOKE_VA.Z))  θ $(size(SMOKE_TR.Θ))")
println("  λ ∈ [$(round(minimum(SMOKE_TR.lam); digits = 3)), " *
        "$(round(maximum(SMOKE_TR.lam); digits = 3))]  (LAMBDA_MIN=$LAMBDA_MIN, LAMBDA_MAX=$LAMBDA_MAX)")
println("  NOT a convergence claim, NOT a calibration claim.")
println("="^78)

const SMOKE_T0 = time()
# FIXED-DATA train form, use_gpu = false on EVERY NeuralEstimators call (D-10 / Pitfall 1).
# NB: the AdamW arguments stay Float64 to match NeuralEstimators' Float64 CosAnneal lr_schedule
# -- a Float32-parameterised optimiser state trips `Optimisers.adjust!` when the schedule feeds a
# Float64 eta (spike/npe/train_npe.jl:129-131).
const SMOKE_EST = train(build_p11_estimator(),
                        SMOKE_ΘTR_STD, SMOKE_ΘVA_STD, SMOKE_TR.Z, SMOKE_VA.Z;
                        epochs = SMOKE_EPOCHS, batchsize = SMOKE_BATCH, use_gpu = false,
                        optimiser = Flux.Optimisers.AdamW(2.5e-4, (0.9, 0.999), 1e-4),
                        verbose = false)
const SMOKE_TRAIN_SECONDS = time() - SMOKE_T0
println("smoke train wall time: $(round(SMOKE_TRAIN_SECONDS; digits = 2)) s")

@testset "SC1f — d_in = 129, D = 8 trains (R10)" verbose = true begin

    @testset "the 129-row input surface is what was trained on" begin
        @test size(SMOKE_TR.Z, 1) == 129
        @test size(SMOKE_VA.Z, 1) == 129
        @test size(SMOKE_TR.Z, 2) == SMOKE_N_TRAIN
        @test size(SMOKE_TR.Θ, 1) == 8
        # The conditioning row genuinely VARIES across samples — a constant row would mean the
        # smoke never exercised the conditioning path (R4's failure mode, rehearsed here).
        @test length(unique(SMOKE_TR.Z[129, :])) > 1
        @test all(0.0 .<= SMOKE_TR.Z[129, :] .<= 1.0)
        # The shift rows respect their own lambda (the D-13 conditional draw).
        @test all(abs(SMOKE_TR.Θ[5, j]) <= SMOKE_TR.lam[j] for j in 1:SMOKE_N_TRAIN)
        @test all(abs(SMOKE_TR.Θ[6, j]) <= SMOKE_TR.lam[j] for j in 1:SMOKE_N_TRAIN)
    end

    @testset "train returned an estimator without error" begin
        @test SMOKE_EST !== nothing
        @test SMOKE_TRAIN_SECONDS > 0
        @test isfinite(SMOKE_TRAIN_SECONDS)
    end

    @testset "sampleposterior returns 8 marginals on an augmented column" begin
        draws = posterior_for(SMOKE_EST, SMOKE_TR.Z[:, 1]; N = SMOKE_N_DRAWS, use_gpu = false)
        @test size(draws, 1) == 8
        @test size(draws, 2) == SMOKE_N_DRAWS
        @test all(isfinite, draws)

        # Reading rho_true (row 1) back through the PORTED transform lands inside the prior box
        # -- the F2 guarantee holding end to end through a real (if untrained) flow.
        orig = StatsBase.reconstruct(SMOKE_ΘZT, draws)
        @test size(orig) == size(draws)
        @test all(SMOKE_LO[p] <= orig[p, j] <= SMOKE_HI[p]
                  for p in axes(orig, 1), j in axes(orig, 2))
        ρ = vec(orig[1, :])
        @test all(-1.0 .<= ρ .<= 1.0)          # row 1 is a correlation

        # The matrix form works too (many datasets in one pass). MEASURED API SHAPE: for K > 1
        # datasets NeuralEstimators v0.2.1 returns a K-element VECTOR of D×N draw matrices, not
        # a single D×N matrix — so a naive `size(·, 1)` on the result reads K, not D. Asserted in
        # the form the API actually returns.
        many = posterior_for(SMOKE_EST, SMOKE_VA.Z[:, 1:3]; N = 8, use_gpu = false)
        @test many isa AbstractVector
        @test length(many) == 3
        @test all(size(m, 1) == 8 for m in many)
        @test all(size(m, 2) == 8 for m in many)
    end

    @testset "the dimension contract is real (a 128-row input throws)" begin
        # If this ever stops throwing, the 129th row has silently stopped being required and the
        # whole D-03 conditioning story would be unverifiable.
        @test_throws Exception posterior_for(SMOKE_EST, randn(Float32, 128); N = 4, use_gpu = false)
        @test_throws Exception posterior_for(SMOKE_EST, randn(Float32, 130); N = 4, use_gpu = false)
    end

    @testset "smoke ran CPU-only (D-10)" begin
        @test !haskey(Pkg.project().dependencies, "CUDA")
        @test !any(p -> occursin("CUDA", p.name), values(Pkg.dependencies()))
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
