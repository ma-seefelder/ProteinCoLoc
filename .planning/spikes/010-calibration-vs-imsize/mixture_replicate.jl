#########################################################################################
# Spike 010c --- does the confirmatory gate's rho_true miscalibration REPLICATE?
#
# The one-shot amended gate reported, on PROD_SEED_V2 over the image-size MIXTURE:
#     rho_true  mean(u) = 0.4666, z = -5.2, KS p = 1.12e-7
# Spike 010b then showed the PER-SIZE location biases are small, non-monotone, and aggregate
# (mixture-weighted) to a standardized bias of about +0.014 -- the WRONG SIGN and an order of
# magnitude too small to produce that rank shift.
#
# So either the effect has a mechanism the per-size decomposition misses, or the gate value was
# seed-specific. This script settles which, by reproducing the gate's SBC arm EXACTLY -- same
# mixture, same weights, same M, same L -- on a DEV seed.
#
# DIAGNOSTIC ONLY. PROD_SEED_V2 is NOT touched; no gate report is written.
#########################################################################################

using Random, Random123, Statistics, Printf, JLD2
using HypothesisTests, Distributions

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const GATE = joinpath(REPO, "test", "gate")
using ProteinCoLoc

include(joinpath(GATE, "gate_consts_8_v2.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(GATE, "harness.jl"))

const DEV = UInt64(0x00A1CE01)
const FORBIDDEN = UInt64[0xC0FFEE, 0x5BC0FFEE, 0x1, 0xDE7C0DE, 0xDE7C0DE2,
                         0x0DE7C0D3, 0x00B7A771, 0x00B7A772, 0x00CA11B0, 0x00CA11B1]
@assert !(DEV in FORBIDDEN) && DEV != PROD_SEED_V2 "DEV seed collides"

const ROOT = joinpath(REPO, "artifacts", "amended_v2", "grid_8")
m   = ProteinCoLoc.load_estimator(joinpath(ROOT, "npe_8.jld2"))
sim = default_simulator()

const M = parse(Int, get(ENV, "SPIKE010C_M", "2000"))
const L = SBC_L
const LAB = ["ρ_true", "spillover", "autofluorescence", "label_efficiency",
             "shift_dx", "shift_dy", "noise"]

println("Replicating the gate SBC arm on a DEV seed.")
println("  mixture = ", SBC_IMSIZE_SET, "  weights = ", SBC_IMSIZE_WEIGHTS)
println("  M = $M, L = $L   (gate used M=", SBC_M, ")")

rng = Random123.Philox4x(UInt64, (DEV, UInt64(8)))
Random.seed!(DEV)                     # sampleposterior draws come from the GLOBAL stream

ranks = zeros(Int, M, 7)
pmean = zeros(M, 7); psd = zeros(M, 7); θtr = zeros(M, 7)
szs   = Vector{Tuple{Int,Int}}(undef, M)
for i in 1:M
    d = draw_simulate_infer(m, rng; G = 8, imsize = SBC_IMSIZE,
                            imsize_set = SBC_IMSIZE_SET,
                            imsize_weights = SBC_IMSIZE_WEIGHTS, N = L, sim = sim)
    θv = (d.θ.ρ_true, d.θ.spillover, d.θ.autofluorescence, d.θ.label_efficiency,
          d.θ.shift_dx, d.θ.shift_dy, d.θ.noise)
    for p in 1:7
        ranks[i, p] = count(<(θv[p]), view(d.draws, p, :))
        θtr[i, p]   = θv[p]
        pmean[i, p] = mean(view(d.draws, p, :))
        psd[i, p]   = std(view(d.draws, p, :))
    end
    szs[i] = d.imsize
    i % 250 == 0 && @printf("  %4d/%d\n", i, M)
end

se = sqrt((1 / 12) / M)
println("\n============ DEV-seed mixture SBC (M=$M) vs the reported gate ============")
@printf("%-18s %10s %8s %12s %14s %14s\n",
        "parameter", "mean(u)", "z", "KS p", "gate mean(u)", "gate z")
gate_u = [0.4666, 0.5278, 0.5165, 0.5119, 0.4946, 0.4857, 0.4880]
gate_z = [-5.2, 4.3, 2.6, 1.8, -0.8, -2.2, -1.9]
for p in 1:7
    u  = (ranks[:, p] .+ 0.5) ./ (L + 1)
    ks = pvalue(ExactOneSampleKSTest(u, Uniform(0, 1)))
    @printf("%-18s %10.4f %8.2f %12.2e %14.4f %14.1f\n",
            LAB[p], mean(u), (mean(u) - 0.5) / se, ks, gate_u[p], gate_z[p])
end

println("\n-- standardized location bias (comparable to spike 010b) --")
for p in 1:7
    z = (pmean[:, p] .- θtr[:, p]) ./ max.(psd[:, p], eps())
    @printf("  %-18s %+.4f  (SE %.4f)\n", LAB[p], mean(z), std(z) / sqrt(M))
end

println("\n-- realised mixture --")
for isz in SBC_IMSIZE_SET
    @printf("  %-12s %4d  (%.3f, target %.3f)\n", string(isz), count(==(isz), szs),
            count(==(isz), szs) / M, SBC_IMSIZE_WEIGHTS[findfirst(==(isz), SBC_IMSIZE_SET)])
end

JLD2.save(joinpath(@__DIR__, "mixture_replicate.jld2"),
          "ranks", ranks, "pmean", pmean, "psd", psd, "theta", θtr, "sizes", szs,
          "M", M, "L", L, "dev_seed", DEV, "labels", LAB)
println("\nsaved mixture_replicate.jld2")
