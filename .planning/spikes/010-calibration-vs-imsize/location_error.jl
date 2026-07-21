#########################################################################################
# Spike 010b --- the LOW-VARIANCE version of the same question.
#
# calib_vs_imsize.jl measured per-size SBC RANKS at M=200. That statistic discards magnitude,
# so its standard error is 1.0 in z units by construction -- exactly the size of the effect
# under investigation. Underpowered by design.
#
# This script measures the POSTERIOR LOCATION ERROR directly:
#     err_i = E[theta_hat_i] - theta*_i      (per parameter, per draw)
# SE(mean err) = sd(err)/sqrt(M), which for rho_true (post_sd ~0.065) is ~0.005 at M=200 --
# an order of magnitude more sensitive than the rank route at the same cost.
#
# Also reports the STANDARDIZED error err/post_sd, which is the quantity a rank shift reflects,
# so the two analyses can be compared directly.
#
# DIAGNOSTIC ONLY -- DEV seed, nothing retrained, no gate report.
#########################################################################################

using Random, Random123, Statistics, Printf, JLD2

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const GATE = joinpath(REPO, "test", "gate")
using ProteinCoLoc

include(joinpath(GATE, "gate_consts_8_v2.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(GATE, "harness.jl"))

const DEV_SEED = UInt64(0x00CA11B1)
const FORBIDDEN = UInt64[0xC0FFEE, 0x5BC0FFEE, 0x1, 0xDE7C0DE, 0xDE7C0DE2,
                         0x0DE7C0D3, 0x00B7A771, 0x00B7A772, 0x00CA11B0]
@assert !(DEV_SEED in FORBIDDEN) && DEV_SEED != PROD_SEED_V2 "DEV_SEED collides"

const ROOT = joinpath(REPO, "artifacts", "amended_v2", "grid_8")
m   = ProteinCoLoc.load_estimator(joinpath(ROOT, "npe_8.jld2"))
sim = default_simulator()

const M_PER = parse(Int, get(ENV, "SPIKE010B_M", "200"))
const L     = SBC_L
const SIZES = ((256, 256), (512, 512), (1024, 1024), (1376, 1028), (2048, 2048))
const LAB   = ["ρ_true", "spillover", "autofluorescence", "label_efficiency",
               "shift_dx", "shift_dy", "noise"]

store = Dict{Tuple{Int,Int},Any}()

for isz in SIZES
    rng = Random123.Philox4x(UInt64, (DEV_SEED, UInt64(hash(isz) & 0xffff)))
    Random.seed!(DEV_SEED + UInt64(hash(isz) & 0xff))
    θtrue = zeros(M_PER, 7); pmean = zeros(M_PER, 7); psd = zeros(M_PER, 7)
    for i in 1:M_PER
        d = draw_simulate_infer(m, rng; G = 8, imsize = isz, imsize_set = (isz,),
                                imsize_weights = (1.0,), N = L, sim = sim)
        θv = (d.θ.ρ_true, d.θ.spillover, d.θ.autofluorescence, d.θ.label_efficiency,
              d.θ.shift_dx, d.θ.shift_dy, d.θ.noise)
        for p in 1:7
            θtrue[i, p] = θv[p]
            pmean[i, p] = mean(view(d.draws, p, :))
            psd[i, p]   = std(view(d.draws, p, :))
        end
    end
    store[isz] = (θtrue = θtrue, pmean = pmean, psd = psd)
    err = pmean[:, 1] .- θtrue[:, 1]
    @printf("%-12s rho_true bias = %+.5f +- %.5f  (%.1f sigma)   median post_sd = %.4f\n",
            string(isz), mean(err), std(err)/sqrt(M_PER),
            mean(err)/(std(err)/sqrt(M_PER)), median(psd[:, 1]))
end

println("\n=========== RAW LOCATION BIAS  E[theta_hat] - theta*  (+- SE) ===========")
@printf("%-18s", "parameter")
for isz in SIZES; @printf("%18s", string(isz)); end; println()
for p in 1:7
    @printf("%-18s", LAB[p])
    for isz in SIZES
        e = store[isz].pmean[:, p] .- store[isz].θtrue[:, p]
        @printf("  %+.4f(%.4f)", mean(e), std(e)/sqrt(M_PER))
    end
    println()
end

println("\n=========== STANDARDIZED bias  mean( (E[theta_hat]-theta*) / post_sd ) ===========")
println("(this is what a rank shift reflects; |value| >~ 0.1 is a meaningful miscentring)")
@printf("%-18s", "parameter")
for isz in SIZES; @printf("%14s", string(isz)); end; println()
for p in 1:7
    @printf("%-18s", LAB[p])
    for isz in SIZES
        s = store[isz]
        z = (s.pmean[:, p] .- s.θtrue[:, p]) ./ max.(s.psd[:, p], eps())
        @printf("%9.3f(%.2f)", mean(z), std(z)/sqrt(M_PER))
    end
    println()
end

println("\n=========== posterior WIDTH vs size (sharpening check, H2) ===========")
@printf("%-18s", "parameter")
for isz in SIZES; @printf("%12s", string(isz)); end; println()
for p in 1:7
    @printf("%-18s", LAB[p])
    for isz in SIZES; @printf("%12.4f", median(store[isz].psd[:, p])); end
    println()
end

JLD2.save(joinpath(@__DIR__, "location_error.jld2"),
          "store", store, "sizes", SIZES, "M", M_PER, "labels", LAB, "dev_seed", DEV_SEED)
println("\nsaved location_error.jld2")
