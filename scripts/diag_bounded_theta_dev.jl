#########################################################################################
# scripts/diag_bounded_theta_dev.jl --- DEV-SEED diagnostic for the F2 bounded-θ experiment.
#
# ANTI-SNOOPING (read this before changing a single constant):
#   * This is a DEVELOPMENT diagnostic, NOT a ship-gate arm. It runs on a DEV seed that is
#     ASSERTED disjoint from every pre-registered `PROD_SEED[G]`, from the spike validation seed
#     `VAL_MASTER_SEED = 0x5BC0FFEE`, from the spike training seed `NPE_MASTER_SEED = 0xC0FFEE`
#     and from the datagen seed `DEFAULT_MASTER_SEED = 0x1`.
#   * It does NOT touch `test/gate/gate_consts_*.jl`, does NOT write a gate report, and must
#     never be reported as a gate result. Burning `PROD_SEED[8]` on a development iteration is
#     exactly the data-snooping failure this project has disciplined itself against.
#
# WHAT IT MEASURES, paired (both nets see the SAME θ* draws, the SAME simulated images and the
# SAME per-draw flow base-RNG seed) so the only difference is the net:
#   * per-parameter SBC mean normalized rank u = (rank+0.5)/(L+1) and z = (mean(u)−0.5)/√((1/12)/M)
#   * per-parameter shrinkage = mean(post_sd)/prior_sd  (a "fix" that merely destroys information
#     makes ranks uniform VACUOUSLY — 07-CALIBRATION-FINDINGS F3 — and must be called out, not
#     celebrated)
#   * out-of-support leakage: the fraction of posterior draws outside the prior box
#
# Run:  julia --project -t auto scripts/diag_bounded_theta_dev.jl [M] [L]
#########################################################################################
using ProteinCoLoc
using StatsBase, Statistics, Random, Printf
import Random123: Philox4x

const M = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 500
const L = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 999
const G = 8
const IMSIZE = (256, 256)

# --- DEV seed (clearly marked, asserted disjoint from every reported/reserved stream) --------
# `0xDE7C0DE` is the develop seed this project already used for the Phase-5 iter1/iter2 probes. A
# 3rd CLI arg supplies an ALTERNATIVE dev seed (used to check that a diagnostic finding is not a
# single-dev-seed fluke); the disjointness assertions below apply to whatever is passed.
const DEV_SEED = length(ARGS) >= 3 ? parse(UInt64, ARGS[3]; base = 16) : 0x0000_0000_0DE7_C0DE
include(joinpath(@__DIR__, "..", "test", "gate", "gate_consts_8.jl"))   # PROD_SEED / forbidden seeds
@assert DEV_SEED != VAL_MASTER_SEED
@assert DEV_SEED != NPE_MASTER_SEED
@assert DEV_SEED != UInt64(ProteinCoLoc.DEFAULT_MASTER_SEED)
for (g, s) in PROD_SEED
    @assert DEV_SEED != s "DEV_SEED collides with the pre-registered PROD_SEED[$g]"
end
@info "DEV seed asserted disjoint from PROD_SEED, VAL_MASTER_SEED, NPE_MASTER_SEED, DEFAULT_MASTER_SEED" DEV_SEED

dev_rng() = Philox4x(UInt64, (DEV_SEED, UInt64(0)))
# Per-draw seed for the GLOBAL stream `sampleposterior` draws its flow base samples from
# (it threads no rng — see test/gate/harness.jl). Derived from DEV_SEED, never from PROD_SEED.
dev_global_seed(i::Integer) = rand(Philox4x(UInt64, (DEV_SEED, UInt64(i))), UInt64)

const LABELS = ("ρ_true", "spillover", "autofluor", "label_eff", "shift_dx", "shift_dy", "noise")
const BOUNDS = ProteinCoLoc.theta_prior_bounds()

nets = Dict{String,Any}(
    "old (unbounded θ)" => ProteinCoLoc.load_estimator(
        normpath(joinpath(@__DIR__, "..", "artifacts", "grid_8", "npe_8.jld2"))),
    "new (bounded θ)"   => ProteinCoLoc.load_estimator(
        normpath(joinpath(@__DIR__, "..", "artifacts", "bounded_theta", "grid_8", "npe_8.jld2"))),
)
order = ["old (unbounded θ)", "new (bounded θ)"]

ranks   = Dict(k => Matrix{Int}(undef, M, 7)     for k in order)
post_sd = Dict(k => Matrix{Float64}(undef, M, 7) for k in order)
leak    = Dict(k => zeros(Float64, 7)            for k in order)
prior_draws = Matrix{Float64}(undef, M, 7)

rng = dev_rng()
t0 = time()
for i in 1:M
    θ   = ProteinCoLoc.sample_prior(rng)                       # SAME θ* for both nets
    mci = ProteinCoLoc.build_mci(ProteinCoLoc.simulate_pair(rng, θ; imsize = IMSIZE))
    S   = ProteinCoLoc.encode_d01(ProteinCoLoc.patch_summary(mci, G))
    θv  = collect(values(θ))
    prior_draws[i, :] = θv
    for k in order
        m = nets[k]
        Random.seed!(dev_global_seed(i))                        # SAME flow base stream per net
        Z  = ProteinCoLoc.standardize_summary(S, m.zt, :min)
        dz = ProteinCoLoc.posterior_for(m.estimator, Z; N = L, use_gpu = false)
        dr = StatsBase.reconstruct(m.θzt, dz)                   # 7×L physical θ
        for p in 1:7
            ranks[k][i, p]   = count(<(θv[p]), @view dr[p, :])
            post_sd[k][i, p] = std(@view dr[p, :])
            leak[k][p]      += count(x -> x < BOUNDS[p][1] || x > BOUNDS[p][2], @view dr[p, :]) / L
        end
    end
    i % 50 == 0 && @info "progress" i minutes=round((time() - t0) / 60, digits = 1)
end

se = sqrt((1 / 12) / M)
println("\n", "="^108)
println("DEV-SEED SBC DIAGNOSTIC  (M = $M, L = $L, imsize = $IMSIZE, grid = $G, DEV_SEED = 0x",
        string(DEV_SEED, base = 16), ")")
println("SE of mean(u) under uniformity = ", round(se, digits = 5))
println("="^108)
@printf("%-11s | %-32s | %-32s | %s\n", "", "OLD (unbounded θ)", "NEW (bounded θ)", "")
@printf("%-11s | %7s %7s %8s %6s | %7s %7s %8s %6s |\n",
        "param", "mean(u)", "z", "shrink", "leak%", "mean(u)", "z", "shrink", "leak%")
println("-"^108)
for p in 1:7
    row = Any[LABELS[p]]
    for k in order
        u  = (ranks[k][:, p] .+ 0.5) ./ (L + 1)
        mu = mean(u)
        z  = (mu - 0.5) / se
        sh = mean(post_sd[k][:, p]) / std(prior_draws[:, p])
        push!(row, mu, z, sh, 100 * leak[k][p] / M)
    end
    @printf("%-11s | %7.4f %7.2f %8.3f %6.2f | %7.4f %7.2f %8.3f %6.2f |\n", row...)
end
println("="^108)
println("shrink = mean(post_sd)/prior_sd  (≈1 ⇒ posterior == prior ⇒ a vacuous rank pass, F3)")
println("leak%  = % of posterior draws outside the prior box (0 by construction for bounded θ)")
