#########################################################################################
# Spike 011 --- is the M=2000 SBC KS test OVER-POWERED for non-identified nuisances?
#
# Finding (spike 010 residual diagnosis): autofluorescence/label_efficiency reject SBC with a
# marginal misplacement of only ~0.05-0.07 posterior-SD. For a NON-identified parameter the
# posterior is ~ the prior for every dataset, so its calibration is entirely a question of how
# accurately the trained flow reproduces the PRIOR MARGINAL. This spike asks, WITHOUT any
# training: what marginal accuracy does the KS test demand at M=2000, and how does that scale
# with M?
#
# Model: a non-identified parameter with a marginal shifted by s (in units of the rank scale,
# i.e. mean(u) = 0.5 + s). This is the leading-order effect of a misplaced learned marginal.
# We measure KS rejection power vs (s, M).  Pure numerics, no model, no seed consumed.
#########################################################################################

using Random, Statistics, Printf, HypothesisTests, Distributions
Random.seed!(0x5BC0)   # local, diagnostic only

# One SBC replicate for a shifted-uniform rank distribution, return KS p.
function ks_p(M, s; L = 999)
    # u = clamp(v + s) with v ~ U(0,1); discretize to L+1 rank levels as SBC does
    u = clamp.(rand(M) .+ s, 0.0, 1.0)
    pvalue(ExactOneSampleKSTest(u, Uniform(0, 1)))
end

power(M, s; reps = 400, α = 0.05) =
    mean(ks_p(M, s) < α for _ in 1:reps)

println("SBC KS rejection power vs marginal shift s (=mean(u)-0.5) and M")
println("(a perfectly calibrated non-identified param has s=0; power there should be ~alpha=0.05)\n")
Ms = (250, 500, 1000, 2000)
ss = (0.0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10)
@printf("%-8s", "shift")
for M in Ms; @printf("%10d", M); end; println()
for s in ss
    @printf("%-8.3f", s)
    for M in Ms; @printf("%10.2f", power(M, s)); end
    println()
end

# The observed nuisance shifts, translated to s:
println("\nObserved nuisance shifts (from spike 010 mixture replicate):")
obs = [("autofluorescence", 0.0165, "REJECT"), ("label_efficiency", 0.0119, "REJECT"),
       ("spillover", 0.0137, "borderline"), ("noise", -0.0120, "ok"),
       ("shift_dx", -0.0054, "ok")]
for (nm, s, st) in obs
    @printf("  %-18s s=%+.4f  power@M2000=%.2f  power@M500=%.2f   (gate: %s)\n",
            nm, s, power(2000, abs(s)), power(500, abs(s)), st)
end

# What shift is "just detectable" (50% power) at each M?
println("\nMarginal shift at 50% rejection power (the smallest miscalibration each M 'sees'):")
for M in Ms
    lo, hi = 0.0, 0.2
    for _ in 1:40
        mid = (lo+hi)/2
        power(M, mid; reps = 300) < 0.5 ? (lo = mid) : (hi = mid)
    end
    s50 = (lo+hi)/2
    # translate s (rank-scale) back to posterior-SD units for a Uniform prior:
    # s = delta_sd / sqrt(12)  =>  delta_sd = s*sqrt(12)
    @printf("  M=%-5d  s50 = %.4f   (~ %.3f posterior-SD of marginal drift)\n", M, s50, s50*sqrt(12))
end

# --- persist curves for the figure (figenv has no HypothesisTests) ---
using JLD2
ss_fine = collect(0.0:0.005:0.10)
curves = Dict(M => [power(M, s; reps=500) for s in ss_fine] for M in Ms)
s50v = Float64[]
for M in Ms
    lo,hi=0.0,0.2
    for _ in 1:36; mid=(lo+hi)/2; (power(M,mid;reps=300)<0.5) ? (lo=mid) : (hi=mid); end
    push!(s50v, (lo+hi)/2*sqrt(12))
end
JLD2.save(joinpath(@__DIR__,"power_data.jld2"),"ss",ss_fine,"Ms",collect(Ms),
          "curves",curves,"s50",s50v)
println("saved power_data.jld2")
