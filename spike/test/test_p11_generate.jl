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

# spike/test/test_p11_generate.jl --- fixture-scale unit gate for the λ-hierarchical sampler.
#
# WHAT THIS FILE DEFENDS: Pitfall 1 / risk R4. If λ were drawn independently of the shift, the
# 129th input row would carry no information, the net would correctly learn to ignore it, and
# SC2 would fail flat for an implementation reason indistinguishable from "below the fixed 8×8
# summary's resolution". The correlation testset below therefore asserts BOTH directions: the
# hierarchical sampler must show a materially positive λ-to-|shift| correlation, AND a
# deliberately INDEPENDENT sampler built inline must NOT — proving the statistic can tell the two
# designs apart, which is exactly what the R4 bug looks like.
#
# Fixture scale only: tiny images, small n, no training, no reported number. Every draw rides
# P11_FIXTURE_COUNTER / a local `Xoshiro`, never a reported counter (D-01).
#
# CPU-only (D-10). Run:
#     julia --project=spike spike/test/test_p11_generate.jl

using Test
using Statistics
using Random
using Pkg

isdefined(@__MODULE__, :sample_p11_theta) ||
    include(joinpath(@__DIR__, "..", "data", "p11_generate.jl"))

# --- Fixture knobs. NOT pre-registered thresholds — Monte-Carlo resolution only. ------------
const GEN_FIX_IMSIZE   = (64, 64)   # tiny on purpose: this file must stay seconds-scale
const GEN_FIX_N_THETA  = 2000       # θ draws for the correlation clauses (no simulation)
const GEN_FIX_N_IMSIZE = 4000       # imsize draws for the weight clause (no simulation)
const GEN_FIX_N_POOL   = 200        # simulated pool rows for the layout clause
const GEN_FIX_N_REPRO  = 50         # simulated pool rows for the thread-independence clause

# A fixed tiny-image sampler, injected so the pool clauses never pay the F5 mixture's cost.
# It is a FIXTURE override, not a second copy of the mixture.
_fix_imsize(rng) = (rand(rng); GEN_FIX_IMSIZE)

@testset "P11 λ-hierarchical datagen (D-03, D-09)" verbose = true begin

    @testset "lambda is drawn first and the shift is conditional" begin
        rows = [sample_p11_theta(p11_datagen_rng(i)) for i in 1:GEN_FIX_N_THETA]
        λ    = [r[2] for r in rows]
        adx  = [abs(r[1].shift_dx) for r in rows]
        ady  = [abs(r[1].shift_dy) for r in rows]

        # The per-sample invariant, on every draw.
        @test all(i -> adx[i] <= λ[i], eachindex(λ))
        @test all(i -> ady[i] <= λ[i], eachindex(λ))
        # λ stays inside its declared range.
        @test all(l -> LAMBDA_MIN <= l <= LAMBDA_MAX, λ)

        # MATERIALLY positive, not merely nonzero: under shift = λ·U(-1,1) the population value
        # is ≈ 0.60, so a 0.3 bar is comfortably clear of Monte-Carlo noise at n = 2000.
        @test cor(λ, adx) > 0.3
        @test cor(λ, ady) > 0.3

        # --- THE OPPOSITE DIRECTION. A deliberately INDEPENDENT sampler — λ and the shift drawn
        #     separately from their own marginals — is exactly what the R4 bug looks like. The
        #     statistic must NOT be materially positive there, or it is stuck at "pass" and the
        #     clause above proves nothing.
        rng_ind = Random.Xoshiro(20261107)
        λ_ind   = [rand(rng_ind, Uniform(LAMBDA_MIN, LAMBDA_MAX)) for _ in 1:GEN_FIX_N_THETA]
        s_ind   = [abs(rand(rng_ind, SHIFT_PRIOR))                for _ in 1:GEN_FIX_N_THETA]
        @test !(cor(λ_ind, s_ind) > 0.3)
        @test abs(cor(λ_ind, s_ind)) < 0.1
    end

    @testset "chromatic_eps is not lambda-scaled" begin
        rows = [sample_p11_theta(p11_datagen_rng(i)) for i in 1:GEN_FIX_N_THETA]
        λ    = [r[2] for r in rows]
        ce   = [r[1].chromatic_eps for r in rows]
        # D-09 / ruling Q3: its own FIXED prior. No λ-dependence, and the realized range stays
        # inside the declared support.
        @test abs(cor(λ, abs.(ce))) < 0.1
        @test all(x -> minimum(CHROMATIC_PRIOR) <= x <= maximum(CHROMATIC_PRIOR), ce)
        @test minimum(ce) < -0.015 && maximum(ce) > 0.015   # the prior is actually exercised
    end

    @testset "theta field order is preserved" begin
        # `collect(values(θ))` is the house θ→vector map, so NamedTuple field order is a
        # load-bearing contract: `fit_p11_theta_transform` aligns rows to the prior box by it.
        θ_p11, _ = sample_p11_theta(p11_datagen_rng(7))
        θ_ref    = sample_prior(Random.Xoshiro(7))
        @test keys(θ_p11) == keys(θ_ref)
        @test length(θ_p11) == P11_THETA_DIM
        @test keys(θ_p11)[1] === :ρ_true            # ρ_true stays ROW 1
        @test keys(θ_p11)[end] === :chromatic_eps   # D-09 appended LAST (positional θ)
    end

    @testset "imsize sampler respects the F5 weights" begin
        rng  = p11_rng(P11_FIXTURE_COUNTER)
        draw = [sample_p11_imsize(rng) for _ in 1:GEN_FIX_N_IMSIZE]
        @test all(z -> z in P11_IMSIZE_SET, draw)
        for (i, z) in enumerate(P11_IMSIZE_SET)
            freq = count(==(z), draw) / GEN_FIX_N_IMSIZE
            # ±0.04 is ≳ 4 binomial SDs at n = 4000 for every weight in the mixture.
            @test isapprox(freq, P11_IMSIZE_WEIGHTS[i]; atol = 0.04)
        end
    end

    @testset "pool rows are 128-row raw summaries with lambda stored separately" begin
        root = mktempdir()
        dir  = generate_p11_pool(GEN_FIX_N_POOL; cache_root = root, shard_size = 100,
                                 imsize_sampler = _fix_imsize, imsize_tag = :fixture,
                                 verbose = false)
        pool = load_p11_pool(dir)

        @test size(pool.summary_min, 1) == 128        # RAW, un-standardized, no λ row appended
        @test size(pool.summary_min, 2) == GEN_FIX_N_POOL
        @test size(pool.theta, 1) == P11_THETA_DIM
        @test length(pool.lambda) == GEN_FIX_N_POOL
        @test pool.global_index == collect(1:GEN_FIX_N_POOL)

        # The λ append did NOT happen at generation time: no summary row is the λ vector, in
        # either its raw or its [0, 1]-encoded form.
        enc = (pool.lambda .- LAMBDA_MIN) ./ (LAMBDA_MAX - LAMBDA_MIN)
        @test !any(r -> pool.summary_min[r, :] == pool.lambda, 1:128)
        @test !any(r -> pool.summary_min[r, :] == enc, 1:128)
        # Nor was anything z-scored here: rows 65:128 are the binary present/absent mask, so on
        # a non-degenerate fixture they are exactly 1.0 — a z-score would have moved them.
        @test all(x -> x == 0.0 || x == 1.0, pool.summary_min[65:128, :])

        # The hierarchical invariant survives the round trip through the pool.
        @test all(j -> abs(pool.theta[5, j]) <= pool.lambda[j], 1:GEN_FIX_N_POOL)
        @test all(j -> abs(pool.theta[6, j]) <= pool.lambda[j], 1:GEN_FIX_N_POOL)
        @test cor(pool.lambda, abs.(pool.theta[5, :])) > 0.3

        # Realized-imsize provenance is derivable from what was stored (F6).
        counts = realized_imsize_counts(pool.imsize)
        @test sum(values(counts)) == GEN_FIX_N_POOL
    end

    @testset "thread-count independence" begin
        # Per-global-index Philox keying means the parallel and serial fills must be
        # BYTE-IDENTICAL. Two separate cache roots so neither run can resume the other's shards.
        a = generate_p11_pool(GEN_FIX_N_REPRO; cache_root = mktempdir(), shard_size = 25,
                              imsize_sampler = _fix_imsize, imsize_tag = :fixture,
                              parallel = false, verbose = false) |> load_p11_pool
        b = generate_p11_pool(GEN_FIX_N_REPRO; cache_root = mktempdir(), shard_size = 25,
                              imsize_sampler = _fix_imsize, imsize_tag = :fixture,
                              parallel = true, verbose = false) |> load_p11_pool
        @test a.theta       == b.theta
        @test a.lambda      == b.lambda
        @test a.summary_min == b.summary_min
        @test a.imsize      == b.imsize
        @test a.global_index == b.global_index
    end

    @testset "the datagen stream is disjoint from every reserved stream (D-01)" begin
        @test P11_DATAGEN_SALT != P11_SALT
        @test P11_DATAGEN_SALT != HOLDOUT_SALT
        @test P11_DATAGEN_SALT != FOLD_SALT
        @test P11_DATAGEN_SALT != PROD_SALT
        @test P11_DATAGEN_SALT != AMEND_SALT
        # The key WORD, not just the salt, differs from this phase's own p11_rng family.
        @test (UInt64(P11_DEV_SEED) ⊻ P11_DATAGEN_SALT) != (UInt64(P11_DEV_SEED) ⊻ P11_SALT)
        # And the DEV seed itself is still none of the forbidden streams.
        @test !(UInt64(P11_DEV_SEED) in _p11_forbidden())
    end

    @testset "the wall-clock ceiling is a BLOCKER that throws (T-11-29)" begin
        # A ceiling of zero minutes must abort the run rather than shrink the image-size arm.
        @test_throws ErrorException generate_p11_pool(
            50; cache_root = mktempdir(), shard_size = 25, imsize_sampler = _fix_imsize,
            imsize_tag = :fixture, wallclock_ceiling_min = 0.0, verbose = false)
    end

    @testset "datagen ran CPU-only (D-10)" begin
        @test !haskey(Pkg.project().dependencies, "CUDA")
        @test !any(p -> occursin("CUDA", p.name), values(Pkg.dependencies()))
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
