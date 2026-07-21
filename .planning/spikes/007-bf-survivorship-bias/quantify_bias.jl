#########################################################################################
# Spike 007 --- survivorship bias in the BF gate's correlation estimate.
#
# The gate keeps only pairs where BOTH log-BF values are finite. Spike 006 showed the drops
# are entirely baseline-side and concentrated at large |Delta rho|. This spike asks: does the
# correlation computed on SURVIVORS differ systematically from the correlation over ALL pairs?
#
# Pure post-processing of .planning/spikes/006-bf-attrition-mechanism/attrition_rows.jld2 --
# every clamped baseline variant is a closed-form function of the stored p_post, so no new
# inference run is needed.
#
# DIAGNOSTIC ONLY.
#########################################################################################

using JLD2, Statistics, Printf

const ROWS_PATH = normpath(joinpath(@__DIR__, "..", "006-bf-attrition-mechanism",
                                    "attrition_rows.jld2"))
d       = JLD2.load(ROWS_PATH)
rows    = d["rows"]
p_prior = d["p_prior"]
log_prior_odds = log(p_prior / (1 - p_prior))

"logBF from a tail probability with an explicit clamp floor eps."
clamped_logbf(p, eps) = begin
    pc = clamp(p, eps, 1 - eps)
    log(pc / (1 - pc)) - log_prior_odds
end

a      = getfield.(rows, :a)          # amortized NRE log-BF (always finite, per spike 006)
p_post = getfield.(rows, :p_post)
dtrue  = getfield.(rows, :dtrue)
keepm  = getfield.(rows, :a_finite) .& getfield.(rows, :b_finite)   # the gate's own filter

@assert all(isfinite, a) "amortized side unexpectedly non-finite -- contradicts spike 006"

println("========== SPIKE 007 -- SURVIVORSHIP BIAS ==========")
@printf("pairs: %d total, %d survive the gate filter (%.0f%% dropped)\n",
        length(rows), count(keepm), 100 * (1 - count(keepm) / length(rows)))
@printf("log_prior_odds = %.6g\n\n", log_prior_odds)

# ---- 1. correlation on survivors vs on ALL pairs, under a finite-preserving baseline -------
println("-- corr(amortized, baseline) : survivors-only vs ALL pairs --")
println("   (a finite-preserving CLAMPED baseline lets 'ALL pairs' be computed at all)")
@printf("   %-12s %-18s %-18s %s\n", "clamp eps", "r survivors", "r ALL pairs", "difference")
for eps in (1e-8, 1e-10, 1e-12, 1e-15)
    b_all  = clamped_logbf.(p_post, eps)
    r_surv = cor(a[keepm], b_all[keepm])
    r_all  = cor(a, b_all)
    @printf("   %-12.0e %-18.4f %-18.4f %+.4f\n", eps, r_surv, r_all, r_all - r_surv)
end

# ---- 2. what the survivor subsample actually is --------------------------------------------
println("\n-- the survivor subsample is the ambiguous middle --")
@printf("   |Δρ_true| survivors : median %.3f   IQR [%.3f, %.3f]   n=%d\n",
        median(abs.(dtrue[keepm])), quantile(abs.(dtrue[keepm]), 0.25),
        quantile(abs.(dtrue[keepm]), 0.75), count(keepm))
@printf("   |Δρ_true| dropped   : median %.3f   IQR [%.3f, %.3f]   n=%d\n",
        median(abs.(dtrue[.!keepm])), quantile(abs.(dtrue[.!keepm]), 0.25),
        quantile(abs.(dtrue[.!keepm]), 0.75), count(.!keepm))

# ---- 3. the clamped baseline saturates -> the dropped pairs carry no information ------------
println("\n-- why the clamp is not a fix either: it SATURATES to a TWO-VALUED column --")
for eps in (1e-8, 1e-12)
    b_all = clamped_logbf.(p_post, eps)
    sat   = b_all[.!keepm]
    u     = sort(unique(round.(sat, digits = 6)))
    @printf("   eps=%.0e : the %d dropped pairs take only %d distinct values: %s\n",
            eps, length(sat), length(u), string(u))
    @printf("             (mean %.4f, sd %.4f -- the sd reflects the +/- SPLIT, not within-group spread)\n",
            mean(sat), std(sat))
end
println("   The clamp maps every one-sided posterior onto +-log((1-eps)/eps) regardless of HOW")
println("   one-sided it is. Sign survives; magnitude is destroyed. Within that group the")
println("   baseline has no resolution left, so any correlation computed with it is driven by")
println("   the arbitrary clamp position, not by the estimator being validated.")

# ---- 4. rank correlation as a tail-robust alternative --------------------------------------
spearman(x, y) = cor(invperm(sortperm(x)), invperm(sortperm(y)))
println("\n-- Spearman (rank) correlation is invariant to the tail mapping --")
for eps in (1e-8, 1e-12, 1e-15)
    b_all = clamped_logbf.(p_post, eps)
    @printf("   eps=%.0e : rho_s survivors = %.4f   rho_s ALL = %.4f\n",
            eps, spearman(a[keepm], b_all[keepm]), spearman(a, b_all))
end
println("   (Ranks of the saturated pairs are tied, so Spearman is stable across eps --")
println("    but ties still cost information; this is a diagnostic, not the remedy.)")
println("====================================================")

JLD2.save(joinpath(@__DIR__, "bias_summary.jld2"),
          "a", a, "p_post", p_post, "dtrue", dtrue, "keep", keepm,
          "log_prior_odds", log_prior_odds)
