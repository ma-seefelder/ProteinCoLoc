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

# spike/test/test_p12_architecture.jl --- unit + smoke gate for `spike/npe/p12_architecture.jl`.
# Run:
#     julia --project=spike -e 'using Test; include("spike/test/test_p12_architecture.jl")'
#
# WHAT THIS FILE IS FOR. Both failure modes this file guards are SILENT.
#
#   (1) A TRANSPOSED reshape rotates the lattice relative to the image. It produces an array of
#       the right size, full of the right numbers, and every downstream shape check passes. It
#       is caught here only because testset 1 rebuilds the column-major encoding LOCALLY and
#       demands exact equality against it, on a deliberately non-symmetric fixture (T-12-14).
#   (2) AN UNDERSPECIFIED "a CNN" silently selects the 7x more expensive topology. The measured
#       cost is 12.5 h against 1.8 h at D = 72 over 50 000 pairs, and nothing about the wrong
#       one is wrong -- it just trains for a week of afternoons. Testset 3 pins the layer chain
#       and names the rejected 1 137 760-parameter alternative as a literal (T-12-17).
#
# Testset 3 additionally asserts, on the COMMENT-STRIPPED source, that no permutation-invariant
# pooling re-entered the summary net (T-12-18): exchangeable pooling over patches is the thing
# this phase exists to replace, so it must not reappear as a convenience layer.
#
# SEED DISCIPLINE. Every DATA draw rides `p12_fix_rng(P12_FIXTURE_COUNTER)` -- the FIXTURE seed
# on the FIXTURE counter (p12_consts.jl section 2). Nothing here touches `p12_rng`: a test must
# never pre-observe a stream a reported number rides on. Flux weight initialisation draws from
# the task-local RNG instead, which is seeded below from a LOCALLY CONSTRUCTED `Xoshiro` -- a
# different RNG family entirely, carrying no reported number and unable to collide with the
# Philox key space.
#
# CPU-only: `use_gpu = false` on EVERY NeuralEstimators call in this file, no exceptions.

using Test
using Random
using Flux
using NeuralEstimators

# ORDER MATTERS: the Tier-1 pre-registration before the unit under test. Both guarded for
# idempotency (this file is included from `test_p12_suite.jl` after `test_p12_consts.jl`).
isdefined(@__MODULE__, :P12_G)              || include(joinpath(@__DIR__, "..", "validation", "p12_consts.jl"))
isdefined(@__MODULE__, :build_p12_estimator) || include(joinpath(@__DIR__, "..", "npe", "p12_architecture.jl"))

const P12_ARCH_SRC = joinpath(@__DIR__, "..", "npe", "p12_architecture.jl")

# Strip everything that is DOCUMENTATION rather than CODE: `#= =#` blocks, `\"\"\" \"\"\"`
# docstrings and `#`-prefixed lines. The anti-pattern names ARE deliberately written in the
# prose of the file under test -- that is how the prohibition stays legible -- so an assertion
# on the raw text would forbid documenting the rule it enforces. Stripping first keeps the
# check pointed at the executable topology, which is the only place the pooling could do harm.
function _p12_strip_noncode(path::AbstractString)
    src = read(path, String)
    src = replace(src, r"(?s)#=.*?=#" => "")
    src = replace(src, r"(?s)\"\"\".*?\"\"\"" => "")
    return join(filter(l -> !startswith(strip(l), "#"), split(src, "\n")), "\n")
end

# Rebuild `encode_d01`'s encoding LOCALLY, with the same column-major rule
# (src/amortized/summary.jl:72-76), so the round-trip is checked against an independent
# construction rather than against the code under test calling itself.
_p12_encode_local(M::AbstractMatrix) =
    vcat(vec(coalesce.(M, 0.0)), vec(Float64.(.!ismissing.(M))))

# `verbose = true` so all NINE sub-testsets print with their own pass counts even when
# everything is green. A collapsed summary cannot distinguish "nine testsets passed" from "one
# testset passed and eight were never written".
@testset verbose = true "P12 CNN estimator surface (12-04)" begin

    # =========================================================================================
    @testset "reshape_summary round-trips encode_d01 exactly (no transpose)" begin
        rng = p12_fix_rng(P12_FIXTURE_COUNTER)
        M = Matrix{Union{Float64,Missing}}(randn(rng, P12_G, P12_G))
        # A deliberately NON-SYMMETRIC missingness pattern: on a symmetric fixture a transposed
        # implementation passes every assertion below.
        M[3, 6] = missing
        M[1, 8] = missing
        M[7, 2] = missing
        Z = _p12_encode_local(M)

        R = reshape_summary(Z, P12_G)
        @test size(R) == (P12_G, P12_G, 2, 1)
        @test eltype(R) === Float32
        # THE assertion a transposed implementation fails. It passes visual inspection either
        # way -- the array is the right size and full of the right numbers -- so equality
        # against the locally rebuilt matrix is the only thing that separates them.
        @test R[:, :, 1, 1] == Float32.(coalesce.(M, 0.0))
        @test R[:, :, 2, 1] == Float32.(.!ismissing.(M))
        # And the transposed variant is genuinely different on this fixture, so the two
        # assertions above are not both satisfiable at once.
        @test R[:, :, 1, 1] != permutedims(Float32.(coalesce.(M, 0.0)))
        @test R[:, :, 2, 1] != permutedims(Float32.(.!ismissing.(M)))

        # The column-major contract, stated through the shared `p12_idx` rather than re-derived.
        for (i, j) in ((1, 1), (3, 6), (7, 2), (8, 8), (2, 5))
            @test R[i, j, 1, 1] == Float32(coalesce(M[i, j], 0.0))
            @test R[i, j, 2, 1] == Float32(!ismissing(M[i, j]))
            @test Float32(Z[p12_idx(i, j)]) == R[i, j, 1, 1]
            @test Float32(Z[P12_G^2 + p12_idx(i, j)]) == R[i, j, 2, 1]
        end

        # Multi-column: data sets live along the LAST dimension and do not bleed into each other.
        M2 = Matrix{Union{Float64,Missing}}(randn(rng, P12_G, P12_G))
        M2[4, 4] = missing
        Zmat = hcat(Z, _p12_encode_local(M2))
        R2 = reshape_summary(Zmat, P12_G)
        @test size(R2) == (P12_G, P12_G, 2, 2)
        @test R2[:, :, 1, 1] == Float32.(coalesce.(M, 0.0))
        @test R2[:, :, 1, 2] == Float32.(coalesce.(M2, 0.0))
        @test R2[:, :, 2, 2] == Float32.(.!ismissing.(M2))

        # A wrong row count is a throw, not a silent reinterpretation.
        @test_throws DimensionMismatch reshape_summary(randn(127, 3), 8)
        @test_throws DimensionMismatch reshape_summary(randn(2 * 8^2 + 1, 3), 8)
    end

    # =========================================================================================
    @testset "mask rows are channel 2 and are never z-scored" begin
        rng = p12_fix_rng(P12_FIXTURE_COUNTER)
        n = P12_G^2
        M = Matrix{Union{Float64,Missing}}(randn(rng, P12_G, P12_G))
        M[5, 2] = missing
        Zm = hcat(_p12_encode_local(M), _p12_encode_local(M), _p12_encode_local(M))

        idxs = [p12_idx(2, 3), p12_idx(7, 7)]
        masked = mask_regions(Zm, idxs)
        @test masked !== Zm                       # a COPY, the input is not mutated
        @test Zm == hcat(_p12_encode_local(M), _p12_encode_local(M), _p12_encode_local(M))
        for r in idxs
            @test all(masked[r,     :] .== 0.0)   # the continuous row
            @test all(masked[n + r, :] .== 0.0)   # AND the mask row
        end
        # Every other entry is bit-identical: masking is local, not a rescaling.
        keep = setdiff(1:(2n), vcat(idxs, idxs .+ n))
        @test masked[keep, :] == Zm[keep, :]

        # The mask lands in channel 2 and stays binary -- nothing z-scored it on the way.
        Rm = reshape_summary(masked, P12_G)
        @test all(v -> v === 0.0f0 || v === 1.0f0, Rm[:, :, 2, :])
        @test Rm[2, 3, 2, 1] == 0.0f0 && Rm[7, 7, 2, 1] == 0.0f0
        @test Rm[5, 2, 2, 1] == 0.0f0            # the genuinely-missing patch, same vocabulary
        @test Rm[1, 1, 2, 1] == 1.0f0

        @test_throws ArgumentError mask_regions(Zm, [0])
        @test_throws ArgumentError mask_regions(Zm, [n + 1])

        # --- the MCAR training-time augmentation (Pitfall 4) ---------------------------------
        Mfull = randn(rng, P12_G, P12_G)               # fully PRESENT, so masked count == k
        Zfull = repeat(vcat(vec(Mfull), ones(n)), 1, 16)
        A = copy(Zfull)
        augment_mask!(A, p12_fix_rng(P12_FIXTURE_COUNTER))
        ks = [count(==(0.0), A[(n+1):(2n), c]) for c in axes(A, 2)]
        @test all(k -> k in P12_MASK_K_SET, ks)
        @test any(k -> k > 0, ks)                      # the augmentation is not a no-op
        # A masked region is masked in BOTH rows, never in one of them.
        for c in axes(A, 2), r in 1:n
            @test (A[r, c] == 0.0) == (A[n + r, c] == 0.0) || Mfull[r] == 0.0
        end
        # Reproducible on a fixed rng: same stream in, same augmentation out.
        B = copy(Zfull)
        augment_mask!(B, p12_fix_rng(P12_FIXTURE_COUNTER))
        @test A == B
        @test sample_mask_k(p12_fix_rng(P12_FIXTURE_COUNTER)) ==
              sample_mask_k(p12_fix_rng(P12_FIXTURE_COUNTER))
        @test sample_mask_k(p12_fix_rng(P12_FIXTURE_COUNTER)) in P12_MASK_K_SET
    end

    # =========================================================================================
    @testset "the lean topology is the pre-registered one (R-7)" begin
        net = build_p12_summary_net()
        nparams = sum(length, Flux.trainables(net))
        @test nparams < 100_000
        # THE REJECTED ALTERNATIVE, NAMED: the fat 2 -> 32 -> 64 -> 64 stack with no channel
        # bottleneck carries 1 137 760 summary parameters and measured 12.5 h at D = 72 over
        # 50 000 pairs, against 70 872 parameters and 1.8 h here. This ceiling is what makes a
        # later "small" widening announce itself instead of showing up as a training bill.
        @test sum(length, Flux.trainables(build_p12_summary_net())) < 200_000
        @test nparams == 70_872

        layers = collect(net.layers)
        convs = filter(l -> l isa Conv, layers)
        @test length(convs) == 3
        # Conv weights are (k1, k2, in, out).
        @test [size(c.weight, 4) for c in convs] == [16, 32, 8]
        @test [size(c.weight, 3) for c in convs] == [2, 16, 32]
        @test size(convs[1].weight)[1:2] == (3, 3)
        @test size(convs[2].weight)[1:2] == (3, 3)
        @test size(convs[3].weight)[1:2] == (1, 1)     # the channel bottleneck

        dens = filter(l -> l isa Dense, layers)
        @test length(dens) == 1
        @test size(dens[1].weight) == (128, 8 * P12_G^2)   # flattened width DERIVED from G
        @test size(net(randn(Float32, P12_G, P12_G, 2, 4))) == (128, 4)

        # The head width follows G rather than a hard-coded 512.
        net6 = build_p12_summary_net(; G = 6, dstar = 32)
        @test size(filter(l -> l isa Dense, collect(net6.layers))[1].weight) == (32, 8 * 36)
        @test size(net6(randn(Float32, 6, 6, 2, 3))) == (32, 3)

        # T-12-18: no permutation-invariant pooling re-entered the executable topology.
        code = _p12_strip_noncode(P12_ARCH_SRC)
        @test !occursin("GlobalMeanPool", code)
        @test !occursin("DeepSet", code)
        @test !occursin("MeanPool", code)
        @test !occursin("mean(", code)
        @test !occursin("Unet", code)
        # ...and the expensive head is absent from the file entirely, prose included.
        raw = read(P12_ARCH_SRC, String)
        @test !occursin("Dense(4096", raw)
        @test !occursin("=> 64,", raw)
    end

    # =========================================================================================
    @testset "build_p12_estimator constructs at D = 72 and guards its arguments" begin
        est = build_p12_estimator()
        @test est !== nothing
        @test est isa PosteriorEstimator
        @test P12_D_MINISPIKE == 72

        @test_throws ArgumentError build_p12_estimator(D = 1)
        @test_throws ArgumentError build_p12_estimator(G = 1)
        @test_throws ArgumentError build_p12_estimator(dstar = 0)
        @test_throws ArgumentError build_p12_estimator(num_coupling_layers = 0)
    end

    # =========================================================================================
    @testset "CNN + wide flow trains one epoch and samples D x N (SMOKE)" begin
        # Flux weight init rides the task-local RNG. Seed it from a locally constructed Xoshiro
        # so this smoke is reproducible; it is NOT a Philox stream and carries no reported
        # number. The data below still rides the FIXTURE Philox stream.
        Random.seed!(Random.default_rng(), 20260729)
        rng = p12_fix_rng(P12_FIXTURE_COUNTER)

        Ztrain = randn(rng, Float32, P12_G, P12_G, 2, 64)
        Zval   = randn(rng, Float32, P12_G, P12_G, 2, 16)
        θtrain = randn(rng, Float32, P12_D_MINISPIKE, 64)
        θval   = randn(rng, Float32, P12_D_MINISPIKE, 16)

        est = build_p12_estimator()
        est = train(est, θtrain, θval, Ztrain, Zval;
                    epochs = 1, batchsize = 16, use_gpu = false)
        @test est !== nothing

        # `use_gpu = false` IS LOAD-BEARING HERE, NOT DECORATION. The pinned v0.2.1 signature is
        # `sampleposterior(estimator::PosteriorEstimator, Z; N = 1000, device = nothing,
        # use_gpu::Bool = true, kwargs...)` (PosteriorEstimator.jl:130) -- the LIBRARY DEFAULT IS
        # TRUE. Omitting it would silently opt this call into the GPU path on any machine that
        # has one, which is the D-10 / Pitfall-1 trap and would make the CPU-only claim of the
        # testset below true only by accident of hardware.
        draws = sampleposterior(est, Ztrain[:, :, :, 1:1]; N = 32, use_gpu = false)
        @test size(draws) == (P12_D_MINISPIKE, 32)
        @test all(isfinite, draws)
    end

    # =========================================================================================
    @testset "theta row layout is stable and appended-last" begin
        rows = P12_THETA_ROWS()
        @test length(rows) == 72
        @test length(rows) == P12_D_MINISPIKE
        @test first(rows) === :c0
        @test last(rows) === :r1
        @test rows[65:71] == [:spillover, :autofluorescence, :label_efficiency,
                              :shift_dx, :shift_dy, :noise, :chromatic_eps]
        @test rows[2] === :c1
        @test rows[1 + P12_K_DEV] === Symbol("c", P12_K_DEV)
        @test length(unique(rows)) == length(rows)
        # rho_true is NOT a theta row in the spike lane -- it is replaced by the field, and the
        # shipped-comparable scalar rho is a DERIVED read-time quantity (R-1).
        @test !(:ρ_true in rows) && !(:rho_true in rows)

        @test p12_theta_index(:r1) == 72
        @test p12_theta_index(:c0) == 1
        @test p12_theta_index(:chromatic_eps) == 71
        @test_throws ArgumentError p12_theta_index(:not_a_row)

        @test length(p12_rows_dev) == P12_K_DEV
        @test p12_row_c0 == 1
        @test p12_rows_dev == 2:64
        @test p12_rows_nuisance == 65:71
        @test p12_row_r1 == 72
        # The named positions agree with the names, so neither can drift alone.
        @test [rows[i] for i in p12_rows_nuisance] == rows[65:71]
        @test rows[p12_row_r1] === :r1
        @test rows[p12_row_c0] === :c0

        # The layout is K-generic, not an 8x8 literal.
        @test length(P12_THETA_ROWS(3)) == 1 + 3 + 7 + 1
        @test last(P12_THETA_ROWS(3)) === :r1
    end

    # =========================================================================================
    @testset "the five declared deviations are the pre-registered five" begin
        @test length(P12_RESEARCH_NET_DEVIATIONS) == 5
        # Each entry must carry its decision tag, so the enumeration cannot drift into prose
        # that no longer says which decision it implements.
        @test occursin("D-02",      P12_RESEARCH_NET_DEVIATIONS[1])
        @test occursin("D-01",      P12_RESEARCH_NET_DEVIATIONS[2])
        @test occursin("D-08",      P12_RESEARCH_NET_DEVIATIONS[2])
        @test occursin("D-05",      P12_RESEARCH_NET_DEVIATIONS[3])
        @test occursin("D-06",      P12_RESEARCH_NET_DEVIATIONS[3])
        @test occursin("Pitfall 4", P12_RESEARCH_NET_DEVIATIONS[4])
        @test occursin("mask",      lowercase(P12_RESEARCH_NET_DEVIATIONS[4]))
        @test occursin("R-3",       P12_RESEARCH_NET_DEVIATIONS[5])
        @test all(d -> d isa AbstractString && !isempty(d), P12_RESEARCH_NET_DEVIATIONS)
    end

    # =========================================================================================
    @testset "z-scoring arm is recorded (R-3)" begin
        # T-12-19: an UNRECORDED arm here makes every later comparison unreadable, because the
        # two arms differ in a way no downstream artifact would show.
        @test P12_ZSCORE_ARM === :per_row
        @test P12_ZSCORE_ARM isa Symbol
    end

    # =========================================================================================
    @testset "architecture surface ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
