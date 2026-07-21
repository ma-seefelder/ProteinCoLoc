#########################################################################################
# Spike 006 --- BF-arm attrition mechanism.
#
# Reproduces the amended BF gate's pair loop (test/gate/run_gate.jl bf_gate) but INSTRUMENTED:
# for every attempted pair it records which side (amortized NRE vs non-clamped KDE baseline) is
# non-finite, the raw tail probability p_post that produced it, the true |Delta rho|, the
# posterior one-sidedness, and the drawn image size.
#
# DIAGNOSTIC ONLY -- runs on a DEV seed, never PROD_SEED_V2. Writes no gate report.
#########################################################################################

using Pkg
using Random, Random123, Statistics, Printf, JLD2

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const GATE = joinpath(REPO, "test", "gate")

using ProteinCoLoc

# Load the amended pre-registration + harness exactly as run_gate does.
include(joinpath(GATE, "gate_consts_8_v2.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(GATE, "harness.jl"))

# ---- DEV seed (asserted disjoint from every reported seed) --------------------------------
const DEV_SEED = UInt64(0x00B7A771)            # "BF ATTRition"
const FORBIDDEN = UInt64[0xC0FFEE, 0x5BC0FFEE, 0x1, 0xDE7C0DE, 0xDE7C0DE2, 0x0DE7C0D3]
@assert !(DEV_SEED in FORBIDDEN) "DEV_SEED collides with a reported/dev seed"
@assert DEV_SEED != PROD_SEED_V2 "DEV_SEED collides with the confirmatory PROD_SEED_V2"
dev_rng() = Random123.Philox4x(UInt64, (DEV_SEED, UInt64(8)))

# ---- load the amended-regime grid-8 bundle ------------------------------------------------
const ROOT = joinpath(REPO, "artifacts", "amended_v2", "grid_8")
m     = ProteinCoLoc.load_estimator(joinpath(ROOT, "npe_8.jld2"))
ratio = ProteinCoLoc.load_ratio(joinpath(ROOT, "ratio_8.jld2"))
sim   = default_simulator()

const N_PAIRS = parse(Int, get(ENV, "SPIKE006_N", "100"))

rng = dev_rng()
Random.seed!(DEV_SEED)   # sampleposterior threads no rng (see harness.jl caveat)

# Prior Delta rho sample for the KDE denominator -- same recipe as bf_gate.
prior_draws = Float64[sim.sample_prior(rng).ρ_true - sim.sample_prior(rng).ρ_true for _ in 1:4000]
p_prior     = ProteinCoLoc._p_gt_threshold_unclamped(prior_draws)

rows = NamedTuple[]
for i in 1:N_PAIRS
    pr = draw_simulate_infer_paired(m, rng; G = 8, imsize = SBC_IMSIZE,
                                    imsize_set = SBC_IMSIZE_SET,
                                    imsize_weights = SBC_IMSIZE_WEIGHTS,
                                    N = SBC_L, sim = sim)
    dr = pr.ρs .- pr.ρc                      # posterior Delta rho draws
    a  = ProteinCoLoc.amortized_log_bf(ratio.estimator,
             ProteinCoLoc.pair_encode(pr.Zs, pr.Zc), ratio.log_prior_odds)
    p_post = ProteinCoLoc._p_gt_threshold_unclamped(dr)
    b  = ProteinCoLoc.kde_log_bf_unclamped([dr], prior_draws)[1]

    push!(rows, (i = i,
                 a = a, b = b,
                 a_finite = isfinite(a), b_finite = isfinite(b),
                 p_post = p_post,
                 frac_pos = mean(dr .> 0),            # empirical one-sidedness of the draws
                 dr_med = median(dr), dr_sd = std(dr),
                 dtrue = pr.θs.ρ_true - pr.θc.ρ_true, # TRUE Delta rho
                 imsize = pr.imsize))
    if i % 10 == 0
        @printf("  %3d/%d  drop=%d\n", i, N_PAIRS, count(r -> !(r.a_finite && r.b_finite), rows))
    end
end

JLD2.save(joinpath(@__DIR__, "attrition_rows.jld2"),
          "rows", rows, "p_prior", p_prior, "dev_seed", DEV_SEED, "n", N_PAIRS)

# ---------------------------------- summary -------------------------------------------------
drop     = filter(r -> !(r.a_finite && r.b_finite), rows)
keep     = filter(r ->   r.a_finite && r.b_finite,  rows)
a_bad    = count(r -> !r.a_finite, rows)
b_bad    = count(r -> !r.b_finite, rows)
both_bad = count(r -> !r.a_finite && !r.b_finite, rows)

println("\n===================== SPIKE 006 -- ATTRITION SOURCE =====================")
@printf("attempted            : %d\n", length(rows))
@printf("dropped              : %d  (%.1f%%)\n", length(drop), 100*length(drop)/length(rows))
@printf("  amortized non-finite: %d\n", a_bad)
@printf("  baseline  non-finite: %d\n", b_bad)
@printf("  BOTH non-finite     : %d\n", both_bad)
@printf("p_prior              : %.6g\n", p_prior)

if !isempty(drop)
    println("\n-- dropped pairs: p_post at the boundary? --")
    @printf("  p_post == 0.0 exactly : %d\n", count(r -> r.p_post == 0.0, drop))
    @printf("  p_post == 1.0 exactly : %d\n", count(r -> r.p_post == 1.0, drop))
    @printf("  p_post strictly inside: %d\n", count(r -> 0.0 < r.p_post < 1.0, drop))
    println("\n-- one-sidedness of the posterior draws --")
    @printf("  dropped : frac_pos in {0,1} exactly = %d/%d\n",
            count(r -> r.frac_pos == 0.0 || r.frac_pos == 1.0, drop), length(drop))
    @printf("  kept    : frac_pos in {0,1} exactly = %d/%d\n",
            count(r -> r.frac_pos == 0.0 || r.frac_pos == 1.0, keep), length(keep))
end

println("\n-- |Delta rho_true|: dropped vs kept --")
@printf("  dropped : median=%.3f  mean=%.3f  n=%d\n",
        median(abs.(getfield.(drop, :dtrue))), mean(abs.(getfield.(drop, :dtrue))), length(drop))
@printf("  kept    : median=%.3f  mean=%.3f  n=%d\n",
        median(abs.(getfield.(keep, :dtrue))), mean(abs.(getfield.(keep, :dtrue))), length(keep))

println("\n-- attrition by image size --")
for isz in SBC_IMSIZE_SET
    sub = filter(r -> r.imsize == isz, rows)
    isempty(sub) && continue
    d = count(r -> !(r.a_finite && r.b_finite), sub)
    @printf("  %-12s n=%3d dropped=%3d (%.1f%%)\n", string(isz), length(sub), d, 100*d/length(sub))
end
println("=========================================================================")
