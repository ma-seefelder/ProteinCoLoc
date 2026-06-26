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

# spike/test/test_simulator.jl --- Phase-2 Wave-0 scaffold: the SIM-03 summary
# contract proven on a synthetic image BEFORE any physics simulator exists, so
# Wave-2 builds against a known-good, frozen summary contract.
#
# The contract is the UNCHANGED src/ summary pipeline reached read-only via
# spike/contract.jl (include() coupling, NOTES §4):
#     2-channel image -> build_mci -> patch()/correlation() -> 8x8 rho -> induced mu
#
# Traps encoded as assertions (src/colocalization.jl):
#   * background-not-zero: _exclude_zero drops 0.0/NaN/missing (line 197); a
#     strictly-positive image keeps every pixel, so all 64 patches survive.
#   * >=15-px patch floor: correlation sets missing when <=15 survivors (line 235).
#     At 256^2 / 8x8 each patch is 32x32 = 1024 px, trivially above the floor.
#   * negative control: a hard-zero background region drops those pixels, starving
#     whole patches below the floor -> materially more `missing` (proves the trap).

using Test
using Random

# include() coupling boundary -> src/ summary functions, read-only (Task 2 impl).
include(joinpath(@__DIR__, "..", "contract.jl"))

@testset "SIM-03 summary contract" verbose = true begin

    # --- Synthetic 2-channel, strictly-positive, correlated image -------------
    # Shared smoothed positive field couples the channels; per-channel positive
    # jitter + a +1.0 baseline guarantees every intensity > 0 (no hard zeros).
    rng   = MersenneTwister(2026)
    field = imfilter(abs.(randn(rng, 256, 256)), Kernel.gaussian(3))
    ch1   = Matrix{Float64}(field .+ 0.5 .* abs.(randn(rng, 256, 256)) .+ 1.0)
    ch2   = Matrix{Float64}(field .+ 0.5 .* abs.(randn(rng, 256, 256)) .+ 1.0)
    data  = [ch1, ch2]

    mci = build_mci(data)

    @testset "MultiChannelImage built per D-11 (pixel_size = pixel-dim tuple)" begin
        @test mci isa MultiChannelImage
        @test mci.pixel_size == size(data[1])        # (256,256) pixels, NOT micrometers
        @test mci.pixel_size isa Tuple{Int, Int}
        @test length(mci.channels) == 2
        @test length(mci.otsu_threshold) == 2
        @test mci.channels == ["ch1", "ch2"]
    end

    rho = summary(mci)

    @testset "8x8 contract: <=64 finite rho, background-not-zero trap respected" begin
        @test size(rho) == (8, 8)
        @test eltype(rho) == Union{Float64, Missing}
        @test count(!ismissing, rho) <= 64               # SIM-03: at most 64 per-patch rho
        @test all(isfinite, skipmissing(rho))            # every present rho is finite
        @test count(ismissing, rho) <= 2                 # strictly-positive image: few/zero missing
        @test isfinite(induced_mu(mci))                  # induced mu is a finite Float64
        @test induced_mu(mci) isa Float64
    end

    @testset "negative control: hard-zero background -> more missing (trap matters)" begin
        ch1z = copy(ch1); ch2z = copy(ch2)
        ch1z[1:200, :] .= 0.0                            # zero most rows of both channels
        ch2z[1:200, :] .= 0.0
        mciz = build_mci([ch1z, ch2z]; name = "zero-bg")
        rhoz = summary(mciz)
        @test count(ismissing, rhoz) > count(ismissing, rho)
    end

end

# --- Wave-2 unit under test: the forward physics simulator (SIM-01) -----------
# include() AFTER contract.jl so build_mci/summary/induced_mu are already in scope.
include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))

# A literal 7-field θ (D-04). Helper rebuilds it with a swapped ρ_true for sweeps.
const θ_BASE = (ρ_true = 0.7, spillover = 0.1, autofluorescence = 0.05,
                label_efficiency = 0.85, shift_dx = 0.3, shift_dy = -0.4,
                noise = 0.5)
_θ(ρ) = merge(θ_BASE, (ρ_true = ρ,))

@testset "SIM-01 forward pipeline" verbose = true begin

    out = simulate_pair(Random.Xoshiro(2026), θ_BASE; imsize = (256, 256))

    @testset "shape / type / finiteness (returns 2× 256×256 Matrix{Float64})" begin
        @test out isa Vector{Matrix{Float64}}
        @test length(out) == 2
        @test all(c -> size(c) == (256, 256), out)
        @test all(c -> all(isfinite, c), out)     # no NaN/Inf (warp fillvalue trap held)
        @test all(c -> all(>=(0.0), c), out)       # non-negative intensities
    end

    @testset "θ / imsize validation rejects bad input (T-02-IV)" begin
        @test_throws ArgumentError simulate_pair(Random.Xoshiro(1), _θ(1.5))         # ρ_true ∉ [-1,1]
        @test_throws ArgumentError simulate_pair(Random.Xoshiro(1), _θ(-1.5))        # ρ_true ∉ [-1,1]
        @test_throws ArgumentError simulate_pair(Random.Xoshiro(1), θ_BASE; imsize = (4, 4))  # dim < 8
        @test_throws ArgumentError simulate_pair(Random.Xoshiro(1), merge(θ_BASE, (noise = NaN,)))  # non-finite field
    end

    @testset "determinism (equal fresh Xoshiro ⇒ identical output, D-14)" begin
        @test simulate_pair(Random.Xoshiro(7), _θ(0.3)) ==
              simulate_pair(Random.Xoshiro(7), _θ(0.3))
    end

end

@testset "SIM-03 on simulator output" verbose = true begin

    out = simulate_pair(Random.Xoshiro(2026), θ_BASE; imsize = (256, 256))
    mci = build_mci(out)

    @test mci isa MultiChannelImage
    @test mci.pixel_size == (256, 256)             # D-11 pixel-dim tuple
    @test length(mci.channels) == 2

    rho = summary(mci)
    @test size(rho) == (8, 8)
    @test eltype(rho) == Union{Float64, Missing}
    @test count(!ismissing, rho) <= 64             # at most 64 per-patch ρ
    @test all(isfinite, skipmissing(rho))          # every present ρ finite (no NaN borders)
    @test count(ismissing, rho) <= 2               # background-not-zero trap held (≥15-px floor)
    @test isfinite(induced_mu(mci))                # induced μ is a finite Float64
    @test induced_mu(mci) isa Float64

end

@testset "D-15 monotone ρ_true → induced-correlation (anti-correlation reachable)" verbose = true begin

    # Sweep ρ_true negative → ~zero → positive with the OTHER nuisances fixed and a
    # fixed fresh seed per point, isolating ρ_true's effect through the REAL summary.
    ρ_grid = [-0.7, 0.0, 0.7]
    μ_ind  = map(ρ_grid) do ρ
        induced_mu(build_mci(simulate_pair(Random.Xoshiro(2026), _θ(ρ); imsize = (256, 256))))
    end

    @test all(isfinite, μ_ind)
    @test issorted(μ_ind)        # monotone increasing induced correlation across the grid
    @test μ_ind[1] < 0.0         # negative ρ_true ⇒ anti-correlation (sign-flip rule, D-15)
    @test μ_ind[end] > 0.0       # positive ρ_true ⇒ positive correlation
    @test μ_ind[end] > μ_ind[1]  # the two tails are well-separated (not saturated near 0)

end

# --- Wave-2 SIM-02: the prior-consistency gate (D-01/D-02/D-16) ----------------
# prior.jl draws μ*~Truncated(Cauchy(0,0.3),-1,1) and sets ρ_true=ghat(μ*) (the
# frozen inverse from calibration.jl) so the induced μ matches the Turing μ-prior.
# Bundles MU_PRIOR + ghat (via ghat.jl) + sample_prior + the pre-declared SIM02_W1_TOL.
include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))

# Wasserstein-1 between two empirical samples (quantile transport), test-local.
_w1_test(a::Vector{Float64}, b::Vector{Float64}) = begin
    n = min(length(a), length(b)); qs = ((1:n) .- 0.5) ./ n
    mean(abs.(quantile(a, qs) .- quantile(b, qs)))
end

@testset "SIM-02 prior consistency" verbose = true begin

    @testset "ĝ monotone over its μ domain (corspearman ≈ 1, clamped tails)" begin
        μgrid = collect(range(GHAT_MU_MIN, GHAT_MU_MAX; length = 50))
        ρ     = ghat.(μgrid)
        @test corspearman(μgrid, ρ) ≥ 0.999     # monotone non-decreasing (T-02-CAL)
        @test issorted(ρ)
        @test ghat(-5.0) == ghat(GHAT_MU_MIN)    # clamp below realized range
        @test ghat(5.0)  == ghat(GHAT_MU_MAX)    # clamp above realized range
        @test all(-1.0 .≤ ρ .≤ 1.0)              # ρ_true stays in simulate_pair's domain
    end

    @testset "sample_prior: 7-field θ, in-range ρ_true, deterministic (D-14)" begin
        θ = sample_prior(Random.Xoshiro(11))
        @test keys(θ) == (:ρ_true, :spillover, :autofluorescence,
                          :label_efficiency, :shift_dx, :shift_dy, :noise)
        @test -1.0 ≤ θ.ρ_true ≤ 1.0
        @test all(isfinite, values(θ))
        @test sample_prior(Random.Xoshiro(11)) == sample_prior(Random.Xoshiro(11))
    end

    @testset "induced μ matches Turing μ-prior within tol (realized range, D-16)" begin
        rng = Random.Xoshiro(404)
        N   = 120
        induced = Float64[]
        for _ in 1:N
            θ = sample_prior(rng)
            m = induced_mu(build_mci(simulate_pair(rng, θ; imsize = (512, 512))))
            isfinite(m) && push!(induced, m)
        end
        # Scope the consistency claim to the physically-realized μ range (D-16):
        # the negative μ-prior tail past GHAT_MU_MIN is prior-only (real anchor > 0).
        target = Truncated(Cauchy(0.0, 0.3), GHAT_MU_MIN, GHAT_MU_MAX)
        ref    = rand(rng, target, length(induced))
        w1     = _w1_test(induced, ref)
        @test w1 < SIM02_W1_TOL                  # pre-declared SIM-02 bar (NOT tuned to pass)
    end

end
