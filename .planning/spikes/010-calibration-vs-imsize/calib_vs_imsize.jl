#########################################################################################
# Spike 010 --- why does rho_true calibration degrade with image size?
#
# The amended confirmatory gate (mixture of 512^2 / 1024^2 / 1376x1028 / 2048^2) gave
# rho_true KS p = 1.12e-7 with a location shift of z = -5.2, where the SAME net at 256^2-only
# was clean (z = -1.25, spike F2 diagnostic). This spike holds everything else fixed and sweeps
# ONLY the image size, computing per-size SBC ranks plus the summary statistics that feed them.
#
# Hypotheses under test:
#   H1 SUMMARY NOT SIZE-INVARIANT -- the patch correlation itself drifts with pixel count at
#      fixed theta, so the net sees a covariate it was not told about. (Design defect.)
#   H2 SHARPENING ONLY -- the posterior narrows with more pixels while staying centred; the
#      rank shift would then come from location drift interacting with narrower posteriors.
#   H3 TRAINING-MIX IMBALANCE -- the net is calibrated at the mixture's dominant size and
#      mis-centred at the rare ones.
#
# DIAGNOSTIC ONLY -- DEV seed, no gate report written, nothing retrained.
#########################################################################################

using Random, Random123, Statistics, Printf, JLD2

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const GATE = joinpath(REPO, "test", "gate")
using ProteinCoLoc

include(joinpath(GATE, "gate_consts_8_v2.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(GATE, "harness.jl"))

const DEV_SEED = UInt64(0x00CA11B0)     # "CALIBration"
const FORBIDDEN = UInt64[0xC0FFEE, 0x5BC0FFEE, 0x1, 0xDE7C0DE, 0xDE7C0DE2,
                         0x0DE7C0D3, 0x00B7A771, 0x00B7A772]
@assert !(DEV_SEED in FORBIDDEN) && DEV_SEED != PROD_SEED_V2 "DEV_SEED collides"

const ROOT = joinpath(REPO, "artifacts", "amended_v2", "grid_8")
m   = ProteinCoLoc.load_estimator(joinpath(ROOT, "npe_8.jld2"))
sim = default_simulator()

const M_PER = parse(Int, get(ENV, "SPIKE010_M", "250"))
const L     = SBC_L
# 256^2 is included as the ORIGINAL training regime even though the amendment excluded it --
# it is the comparison point that makes the trend interpretable.
const SIZES = ((256, 256), (512, 512), (1024, 1024), (1376, 1028), (2048, 2048))

se(M) = sqrt((1 / 12) / M)
labels = ["ρ_true", "spillover", "autofluorescence", "label_efficiency",
          "shift_dx", "shift_dy", "noise"]

# ---- H1 probe: at FIXED theta, does the summary drift with image size? ---------------------
# Uses one fixed prior draw, re-simulated at each size with a fresh stream per size.
println("== H1: is the summary size-invariant at fixed theta? ==")
# Robustness: repeat over SEVERAL fixed theta so the drift cannot be theta-specific.
const N_THETA = 4
h1 = Dict{Tuple{Int,Int},Vector{Float64}}()   # per size: one mean per theta
for isz in SIZES; h1[isz] = Float64[]; end
for t in 1:N_THETA
    rng0 = Random123.Philox4x(UInt64, (DEV_SEED, UInt64(90 + t)))
    θfix = sim.sample_prior(rng0)
    @printf("theta %d: rho_true=%+.4f label_eff=%.4f noise=%.4f
",
            t, θfix.ρ_true, θfix.label_efficiency, θfix.noise)
    for isz in SIZES
        vals = Float64[]
        for r in 1:6
            rng = Random123.Philox4x(UInt64, (DEV_SEED, UInt64(1000 + 100t + r)))
            mci = sim.build_mci(sim.simulate_pair(rng, θfix; imsize = isz))
            cr  = ProteinCoLoc.patch_summary(mci, 8)[1:64]
            push!(vals, mean(filter(isfinite, cr)))
        end
        push!(h1[isz], mean(vals))
        @printf("   %-12s mean patch-corr = %+.4f (sd over 6 reps %.4f)
",
                string(isz), mean(vals), std(vals))
    end
end
println("
  DRIFT vs 256^2, per theta (a size-INVARIANT summary would give all zeros):")
@printf("  %-12s", "size")
for t in 1:N_THETA; @printf("%10s", "theta$t"); end; println("      mean")
for isz in SIZES
    @printf("  %-12s", string(isz))
    d = h1[isz] .- h1[(256,256)]
    for t in 1:N_THETA; @printf("%10.4f", d[t]); end
    @printf("%10.4f
", mean(d))
end

# ---- per-size SBC ranks --------------------------------------------------------------------
println("\n== per-size SBC (M=$M_PER per size, L=$L) ==")
results = Dict{Tuple{Int,Int},Any}()
for isz in SIZES
    rng = Random123.Philox4x(UInt64, (DEV_SEED, UInt64(hash(isz) & 0xffff)))
    Random.seed!(DEV_SEED + UInt64(hash(isz) & 0xff))   # sampleposterior uses the global stream
    ranks = zeros(Int, M_PER, 7)
    psd   = zeros(M_PER, 7)
    for i in 1:M_PER
        d = draw_simulate_infer(m, rng; G = 8, imsize = isz, imsize_set = (isz,),
                                imsize_weights = (1.0,), N = L, sim = sim)
        θv = (d.θ.ρ_true, d.θ.spillover, d.θ.autofluorescence, d.θ.label_efficiency,
              d.θ.shift_dx, d.θ.shift_dy, d.θ.noise)
        for p in 1:7                                   # identical to sbc.jl's rank definition
            ranks[i, p] = count(<(θv[p]), view(d.draws, p, :))
            psd[i, p]   = std(view(d.draws, p, :))
        end
    end
    results[isz] = (ranks = ranks, psd = psd)
    u = (ranks[:, 1] .+ 0.5) ./ (L + 1)
    @printf("  %-12s rho_true: mean(u)=%.4f  z=%+6.2f   median post_sd=%.4f\n",
            string(isz), mean(u), (mean(u) - 0.5) / se(M_PER), median(psd[:, 1]))
end

println("\n== per-size z(location) for every parameter ==")
@printf("%-18s", "parameter")
for isz in SIZES; @printf("%12s", string(isz)); end
println()
for j in 1:7
    @printf("%-18s", labels[j])
    for isz in SIZES
        u = (results[isz].ranks[:, j] .+ 0.5) ./ (L + 1)
        @printf("%12.2f", (mean(u) - 0.5) / se(M_PER))
    end
    println()
end

JLD2.save(joinpath(@__DIR__, "calib_vs_imsize.jld2"),
          "results", results, "h1", h1, "sizes", SIZES, "M", M_PER, "L", L,
          "dev_seed", DEV_SEED, "labels", labels)
println("\nsaved calib_vs_imsize.jld2")
