#########################################################################################
# Spike 008 --- remedies for the clamp-vs-Inf tension in the KDE Bayes-factor baseline.
#
# The tension (spikes 006/007):
#   * CLAMPED at 1e-8  -> every fully one-sided posterior SATURATES at logBF ~ +-18.42.
#                         A constant column carries no information (Phase-5 "tail artifact").
#   * UNCLAMPED        -> the same pairs give exactly p=0 or p=1 -> +-Inf -> the gate DROPS them,
#                         starving n AND biasing the survivors toward small |Delta rho|.
#
# Both failures share ONE root cause: the tail probability is computed in LINEAR space, where
# 1 - p underflows to 0 as soon as p is within ~1e-16 of 1.
#
# REMEDY C: compute the tail in LOG space. For a Gaussian KDE with bandwidth h centred on the
# draws x_i,  P(X > t) = (1/n) sum_i Phi((x_i - t)/h), so
#     log P     = logsumexp_i logPhi( (x_i - t)/h ) - log n
#     log(1-P)  = logsumexp_i logPhi(-(x_i - t)/h ) - log n
#     logBF     = (log P - log(1-P)) - log_prior_odds
# `logcdf(Normal(), z)` is stable into the far tail, so logBF stays FINITE and UNSATURATED
# for one-sided posteriors -- no clamp, no Inf, no dropped pair.
#
# DIAGNOSTIC ONLY -- DEV seed, writes no gate report.
#########################################################################################

using Random, Random123, Statistics, Printf, JLD2
using Distributions, KernelDensity, QuadGK

"Stable log-sum-exp (avoids a LogExpFunctions dependency in the root project)."
function logsumexp(v)
    m = maximum(v)
    isfinite(m) || return m
    return m + log(sum(exp.(v .- m)))
end

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const GATE = joinpath(REPO, "test", "gate")
using ProteinCoLoc

include(joinpath(GATE, "gate_consts_8_v2.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(GATE, "harness.jl"))

# ---------------- the three baselines ------------------------------------------------------

"Silverman bandwidth, matching KernelDensity.jl's default rule closely enough for comparison."
function _bw(x)
    n = length(x); s = std(x)
    iqr = quantile(x, 0.75) - quantile(x, 0.25)
    a = iqr > 0 ? min(s, iqr / 1.349) : s
    a <= 0 && (a = eps())
    return 0.9 * a * n^(-1 / 5)
end

"REMEDY C -- log-space Gaussian-KDE tail. Returns (logP, log1mP)."
function log_tail(draws; threshold::Real = 0.0)
    x = collect(float.(draws)); h = _bw(x)
    z = (x .- threshold) ./ h
    lP   = logsumexp(logcdf.(Normal(), z))   - log(length(x))
    l1mP = logsumexp(logcdf.(Normal(), -z))  - log(length(x))
    return lP, l1mP
end

logspace_logbf(draws, log_prior_odds; threshold = 0.0) = begin
    lP, l1mP = log_tail(draws; threshold = threshold)
    (lP - l1mP) - log_prior_odds
end

"STATUS QUO ANTE -- clamped at eps (the spike baseline's _clampp)."
function clamped_logbf(draws, log_prior_odds; threshold = 0.0, eps_ = 1e-8)
    dist = kde(collect(float.(draws)))
    p_le, _ = quadgk(x -> pdf(dist, x), -Inf, threshold)
    p = clamp(1 - p_le, eps_, 1 - eps_)
    log(p / (1 - p)) - log_prior_odds
end

# ---------------- DEV seed -----------------------------------------------------------------
const DEV_SEED = UInt64(0x00B7A772)
const FORBIDDEN = UInt64[0xC0FFEE, 0x5BC0FFEE, 0x1, 0xDE7C0DE, 0xDE7C0DE2, 0x0DE7C0D3, 0x00B7A771]
@assert !(DEV_SEED in FORBIDDEN) && DEV_SEED != PROD_SEED_V2 "DEV_SEED collides"
dev_rng() = Random123.Philox4x(UInt64, (DEV_SEED, UInt64(8)))

const ROOT = joinpath(REPO, "artifacts", "amended_v2", "grid_8")
m     = ProteinCoLoc.load_estimator(joinpath(ROOT, "npe_8.jld2"))
ratio = ProteinCoLoc.load_ratio(joinpath(ROOT, "ratio_8.jld2"))
sim   = default_simulator()

const N = parse(Int, get(ENV, "SPIKE008_N", "100"))

rng = dev_rng(); Random.seed!(DEV_SEED)
prior_draws = Float64[sim.sample_prior(rng).ρ_true - sim.sample_prior(rng).ρ_true for _ in 1:4000]
p_prior     = ProteinCoLoc._p_gt_threshold_unclamped(prior_draws)
log_prior_odds = log(p_prior / (1 - p_prior))

rows = NamedTuple[]
for i in 1:N
    pr = draw_simulate_infer_paired(m, rng; G = 8, imsize = SBC_IMSIZE,
                                    imsize_set = SBC_IMSIZE_SET,
                                    imsize_weights = SBC_IMSIZE_WEIGHTS,
                                    N = SBC_L, sim = sim)
    dr = pr.ρs .- pr.ρc
    a  = ProteinCoLoc.amortized_log_bf(ratio.estimator,
             ProteinCoLoc.pair_encode(pr.Zs, pr.Zc), ratio.log_prior_odds)
    b_unc  = ProteinCoLoc.kde_log_bf_unclamped([dr], prior_draws)[1]
    b_clp  = clamped_logbf(dr, log_prior_odds)
    b_log  = logspace_logbf(dr, log_prior_odds)
    push!(rows, (i = i, a = a, b_unc = b_unc, b_clp = b_clp, b_log = b_log,
                 dtrue = pr.θs.ρ_true - pr.θc.ρ_true, imsize = pr.imsize,
                 frac_pos = mean(dr .> 0)))
    i % 10 == 0 && @printf("  %3d/%d\n", i, N)
end

JLD2.save(joinpath(@__DIR__, "remedies_rows.jld2"), "rows", rows,
          "log_prior_odds", log_prior_odds, "dev_seed", DEV_SEED)

a     = getfield.(rows, :a)
b_unc = getfield.(rows, :b_unc); b_clp = getfield.(rows, :b_clp); b_log = getfield.(rows, :b_log)
fin   = isfinite.(b_unc)

println("\n=============== SPIKE 008 -- BASELINE REMEDIES ===============")
@printf("pairs = %d\n\n", N)
@printf("%-26s %-10s %-12s %-12s %-12s\n", "baseline", "finite", "r (all)", "r (survivors)", "saturated?")
@printf("%-26s %-10s %-12s %-12s %-12s\n", "unclamped (status quo)",
        string(count(fin), "/", N), "n/a", @sprintf("%.4f", cor(a[fin], b_unc[fin])), "-> +-Inf")
sat_clp = count(x -> abs(abs(x) - 18.42) < 0.2, b_clp)
@printf("%-26s %-10s %-12s %-12s %-12s\n", "clamped 1e-8",
        string(N, "/", N), @sprintf("%.4f", cor(a, b_clp)),
        @sprintf("%.4f", cor(a[fin], b_clp[fin])), string(sat_clp, " at +-18.42"))
sat_log = count(x -> abs(abs(x) - 18.42) < 0.2, b_log)
@printf("%-26s %-10s %-12s %-12s %-12s\n", "LOG-SPACE (remedy C)",
        string(count(isfinite, b_log), "/", N), @sprintf("%.4f", cor(a, b_log)),
        @sprintf("%.4f", cor(a[fin], b_log[fin])), string(sat_log, " at +-18.42"))

println("\n-- agreement check: in the NON-degenerate regime the log-space baseline must")
println("   reproduce the existing unclamped baseline (else it is a different estimator) --")
if any(fin)
    d = abs.(b_log[fin] .- b_unc[fin])
    @printf("   max|log-space - unclamped| over the %d finite pairs = %.4g  (median %.4g)\n",
            count(fin), maximum(d), median(d))
end

println("\n-- range of the log-space baseline where the others fail --")
if any(.!fin)
    @printf("   on the %d pairs the gate DROPS: log-space logBF spans [%.2f, %.2f]\n",
            count(.!fin), minimum(b_log[.!fin]), maximum(b_log[.!fin]))
    @printf("   clamped gives them a constant %.2f (sd %.2e)\n",
            mean(b_clp[.!fin]), std(b_clp[.!fin]))
end
println("==============================================================")
