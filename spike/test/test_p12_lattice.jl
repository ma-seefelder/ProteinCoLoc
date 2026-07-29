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

# spike/test/test_p12_lattice.jl --- unit gate for `spike/simulator/p12_lattice.jl`, the only
# genuinely new arithmetic in Phase 12. Run:
#     julia --project=spike -e 'using Test; include("spike/test/test_p12_lattice.jl")'
#
# WHAT THIS FILE IS FOR. Every failure mode of the lattice arithmetic is SILENT. A transposed
# index rotates the lattice relative to the image and passes every symmetric test. A missing
# per-cell rescale breaks SIM-02 only at the boundary. A monotone-but-mis-solved r1 gives two
# arms that are compared at different spatial correlations while both report the same number.
# None of those throws. So each is asserted here, and where the assertion could pass
# vacuously it is paired with a FALSIFIER -- a companion assertion that FAILS if the property
# under test is removed:
#
#   testset 2 also asserts the UNRESCALED marginals are NOT in the pass band, so the test
#     cannot be satisfied by a sampler that never divides;
#   testset 4 also asserts alpha = 0.5 buys almost no spatial structure, which is the measured
#     fact (R-8) that made the r1 reparametrization necessary in the first place;
#   testset 8 asserts the column-major contract on a concrete matrix rather than restating it.
#
# The unit-marginal claim is asserted on DRAWS, never on `Sigma ./ (sd*sd')` -- that expression
# has a unit diagonal for ANY SPD Sigma, so it is green whether or not `field_sampler` applies
# the rescale at all. The rescale lives in the sampler, so only a draw can see it.
#
# SEED DISCIPLINE. All randomness rides `p12_fix_rng(P12_FIXTURE_COUNTER)` -- the FIXTURE
# stream on the FIXTURE counter (p12_consts.jl §2). Nothing here touches `p12_rng`: a test
# must never pre-observe a stream a reported number rides on.
#
# CPU-only. No training, no simulation, no file write.

using Test
using Statistics
using LinearAlgebra
using SparseArrays

# ORDER MATTERS: the Tier-1 pre-registration before the unit under test. Both guarded for
# idempotency (the file is included from `test_p12_suite.jl` after `test_p12_consts.jl`).
isdefined(@__MODULE__, :P12_G)         || include(joinpath(@__DIR__, "..", "validation", "p12_consts.jl"))
isdefined(@__MODULE__, :field_sampler) || include(joinpath(@__DIR__, "..", "simulator", "p12_lattice.jl"))

# Draw `m` fields from `Σ` on the FIXTURE stream, returned as a 64×m matrix (one field per
# column, matching the column-major layout everything else in this phase uses).
function _p12l_draws(Σ, m::Integer)
    rng = p12_fix_rng(P12_FIXTURE_COUNTER)
    f   = field_sampler(Σ)
    X   = Matrix{Float64}(undef, size(Σ, 1), m)
    for k in 1:m
        X[:, k] = f(rng)
    end
    return X
end

@testset "CAR and GP covariances are SPD" begin
    for α in (0.5, 0.9, 0.95, 0.99)
        Σ = car_sigma(P12_G, α)
        @test isposdef(Matrix(Σ))
        @test issymmetric(Σ)
        @test all(isfinite, Matrix(Σ))
    end
    for ℓ in (0.5, 1.0, 2.0, 4.0, 8.0)
        Σ = gp_sigma(P12_G, ℓ)
        @test isposdef(Matrix(Σ))
        @test issymmetric(Σ)
        @test all(isfinite, Matrix(Σ))
    end
    # α = 1 is the improper (intrinsic) limit, excluded by P12_CAR_ALPHA_BRACKET; α < 0 is not
    # a CAR at all. Both are guard errors, not silently-returned matrices.
    @test_throws ArgumentError car_sigma(8, 1.0)
    @test_throws ArgumentError car_sigma(8, -0.1)
end

@testset "the sqrt(diag) rescale gives per-cell marginal sd 1 (R-8)" begin
    Σ = car_sigma(P12_G, 0.95)
    X = _p12l_draws(Σ, 20_000)
    sds = vec(std(X; dims = 2))
    @test length(sds) == P12_G^2
    @test all(s -> 0.97 <= s <= 1.03, sds)
    # THE FALSIFIER. Σ itself does NOT have unit marginals -- measured [0.658, 0.992] at
    # α = 0.95 -- so the band above is a property of the SAMPLER's `./ sd`, not of Σ. If the
    # rescale were removed the draws would inherit this spread and the assertion above would
    # fail (measured max |sd − 1| = 0.363 without the rescale against 0.026 with it).
    @test !(0.97 <= minimum(sqrt.(diag(car_sigma(8, 0.95)))) <= 1.03)
end

@testset "induced lag-1 correlation is monotone in the arm parameter" begin
    car_r1 = [induced_lag1(car_sigma(P12_G, α)) for α in (0.5, 0.9, 0.95, 0.99, 0.999)]
    gp_r1  = [induced_lag1(gp_sigma(P12_G, ℓ))  for ℓ in (0.5, 1.0, 2.0, 4.0, 8.0)]
    @test all(>(0), diff(car_r1))
    @test all(>(0), diff(gp_r1))
    # Monotonicity is what makes the bisection in `car_alpha_for_r1` / `gp_ell_for_r1` a
    # well-posed inverse rather than a search that can land anywhere.
    @test all(r -> 0 < r < 1, car_r1)
    @test all(r -> 0 < r < 1, gp_r1)
end

@testset "the r1 reparametrization solves both arms (R-8)" begin
    for r1 in P12_R1_LADDER
        @test abs(induced_lag1(lattice_sigma(:car, r1)) - r1) <= 10 * P12_R1_SOLVE_TOL
        @test abs(induced_lag1(lattice_sigma(:gp,  r1)) - r1) <= 10 * P12_R1_SOLVE_TOL
    end
    # THE FALSIFIER THAT MOTIVATES THE WHOLE REPARAMETRIZATION: half of α's nominal range buys
    # almost no spatial structure, so α is NOT the knob and a uniform prior on it would put
    # most of its mass on "no spatial field" (R-8). Asserted here so the suite itself records
    # why α was rejected.
    @test induced_lag1(car_sigma(8, 0.5)) < 0.2
    # Out-of-bracket requests are refused rather than clamped to an endpoint.
    @test_throws ArgumentError car_alpha_for_r1(0.99999999)
    @test_throws ArgumentError lattice_sigma(:nonsense, 0.5)
end

@testset "the ablation arm is a genuine independent field (D-10)" begin
    @test lattice_sigma(:none, P12_ABLATION_R1) == Symmetric(Matrix(1.0I, 64, 64))
    # Any arm AT the ablation r1 is the same independent field, so the ablation goes through
    # the identical sampler/copula/datagen path rather than being special-cased at read time.
    @test lattice_sigma(:car, P12_ABLATION_R1) == Symmetric(Matrix(1.0I, 64, 64))
    X = _p12l_draws(lattice_sigma(:none, P12_ABLATION_R1), 20_000)
    @test abs(induced_lag1(Symmetric(cov(X; dims = 2)))) < 0.03
    @test all(s -> 0.97 <= s <= 1.03, vec(std(X; dims = 2)))
end

@testset "DCT-II is an exact orthonormal bijection (D-04/D-07)" begin
    C = p12_dct_matrix(8)
    @test maximum(abs, C'C - I) < 1e-12
    rng = p12_fix_rng(P12_FIXTURE_COUNTER)
    F   = randn(rng, 8, 8)
    @test maximum(abs, p12_idct(p12_dct(F)) - F) < 1e-12
    @test sum(abs2, p12_dct(F)) ≈ sum(abs2, F)
    # The flat (column-major) form the θ rows are actually stored in must agree with the
    # matrix form, and must round-trip just as exactly.
    @test p12_dct_vec(vec(F)) ≈ vec(p12_dct(F))
    @test maximum(abs, p12_idct_vec(p12_dct_vec(vec(F))) - vec(F)) < 1e-12
    # The smoothness ordering is a frozen permutation with the DC coefficient always first,
    # which is what makes "the first K modes" a defined object.
    @test length(p12_dct_order(8)) == 64
    @test first(p12_dct_order(8)) == 1
    @test p12_dct_order(8) |> unique |> length == 64
    @test sort(p12_dct_order(8)) == collect(1:64)
    # D-04's "high rank" is therefore literally "every non-constant mode".
    @test P12_K_DEV == 64 - 1
end

@testset "the DC coefficient is G times the field mean" begin
    rng = p12_fix_rng(P12_FIXTURE_COUNTER)
    F   = randn(rng, 8, 8)
    @test p12_dct(F)[1, 1] ≈ 8 * mean(F)
    # R-1's derived read-time scalar ρ is this identity, so it is asserted rather than assumed:
    # the global term is a rescaling of the field mean, not a separately-defined quantity.
    @test p12_dct(fill(1.0, 8, 8))[1, 1] ≈ 8.0
end

@testset "column-major contract matches encode_d01 (no transpose)" begin
    # `src/amortized/summary.jl` is NOT included here: src/ stays untouched AND unloaded during
    # the spike, so the contract is asserted DIRECTLY against `vec`/`reshape` -- which is the
    # operation `encode_d01` performs -- rather than against the shipped function.
    rng = p12_fix_rng(P12_FIXTURE_COUNTER)
    M   = Matrix{Union{Float64,Missing}}(randn(rng, 8, 8))
    M[3, 6] = missing                     # the Union eltype is exercised, not merely declared
    M[8, 1] = missing
    v = vec(M)
    for j in 1:8, i in 1:8
        @test isequal(v[p12_idx(i, j)], M[i, j])
    end
    @test isequal(reshape(v, 8, 8), M)
    # And the transposed variant is genuinely a DIFFERENT map, so the assertion above has
    # content: an off-by-transpose would rotate the lattice relative to the image.
    @test p12_idx(3, 6) == (6 - 1) * 8 + 3
    @test p12_idx(3, 6) != p12_idx(6, 3)
    # The lattice geometry that rides on the contract: adjacency and the radial design matrix
    # are indexed the same way, so the S-4 guard cannot disagree with the prior about which
    # cell is which.
    W = lattice_adjacency(8)
    @test W[p12_idx(1, 1), p12_idx(2, 1)] == 1.0    # (1,1) and (2,1) are 4-neighbours
    @test W[p12_idx(1, 1), p12_idx(2, 2)] == 0.0    # diagonal neighbours are not
    @test size(radial_basis(8)) == (64, 2)
    @test radial_basis(8)[p12_idx(1, 1), 2] ≈ radial_basis(8)[p12_idx(8, 8), 2]  # corners tie
end

@testset "lattice ran CPU-only" begin
    @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
end
