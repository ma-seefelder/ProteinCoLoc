#=
Spike 014 figures --- rendered from the ISOLATED figenv (CairoMakie is NOT a root dep).

    julia --project=.planning/spikes/figenv .planning/spikes/014-bf-sim-validation/figures.jl

Reads the plain-array results (artifacts/spike014/results.jld2) written by run_sim_bf.jl and
draws a 4-panel figure:
  (A) ROC of the NRE log-BF (H1 vs H0) with AUC.
  (B) NRE log-BF vs signed Δρ (evidence) -- monotonicity, binned-median overlay, decision line.
  (C) Reliability of p̂(coloc) -- decision calibration.
  (D) KDE-baseline pathology on the SAME pairs: |Δρ| of dropped (non-finite) vs kept pairs.
=#
import CairoMakie as MK
import JLD2

const HERE = @__DIR__
d = JLD2.load(joinpath(HERE, "..", "..", "..", "artifacts", "spike014", "results.jld2"))

drho   = d["drho"]; label = Bool.(d["label"])
lognre = d["logbf_nre"]
logunc = d["logbf_kde_unclamped"]
finite = isfinite.(logunc)
auc_nre = d["auc_nre"]; spear = d["spear"]; ece = d["ece"]
frac_nonfinite = d["frac_nonfinite"]; frac_saturated = d["frac_saturated"]

fig = MK.Figure(size = (1080, 820))

# (A) ROC ----------------------------------------------------------------------------------
axA = MK.Axis(fig[1, 1], title = "(A) Discrimination: NRE log-BF ROC",
              xlabel = "false-positive rate", ylabel = "true-positive rate")
MK.lines!(axA, [0, 1], [0, 1], color = (:gray, 0.6), linestyle = :dash)
MK.lines!(axA, d["roc_fpr"], d["roc_tpr"], color = :steelblue, linewidth = 2.5)
MK.text!(axA, 0.42, 0.12; text = "AUC = $(round(auc_nre; digits=3))", fontsize = 18)
MK.xlims!(axA, -0.02, 1.02); MK.ylims!(axA, -0.02, 1.02)

# (B) log-BF vs Δρ -------------------------------------------------------------------------
axB = MK.Axis(fig[1, 2], title = "(B) Monotonicity in evidence  (Spearman = $(round(spear; digits=3)))",
              xlabel = "true signed Δρ = ρ_sample − ρ_control", ylabel = "NRE log Bayes factor")
MK.hlines!(axB, [0.0], color = (:gray, 0.5), linestyle = :dash)
MK.vlines!(axB, [0.0], color = (:gray, 0.5), linestyle = :dash)
MK.scatter!(axB, drho, lognre, color = [l ? :firebrick : :navy for l in label],
            markersize = 5, alpha = 0.5)
MK.lines!(axB, d["bin_centers"], d["bin_med"], color = :black, linewidth = 2.5)
MK.scatter!(axB, d["bin_centers"], d["bin_med"], color = :black, markersize = 9)
MK.band!(axB, d["bin_centers"], d["bin_lo"], d["bin_hi"], color = (:black, 0.12))

# (C) reliability --------------------------------------------------------------------------
axC = MK.Axis(fig[2, 1], title = "(C) Decision calibration: p̂(coloc)  (ECE = $(round(ece; digits=3)))",
              xlabel = "predicted P(coloc)", ylabel = "empirical coloc fraction")
MK.lines!(axC, [0, 1], [0, 1], color = (:gray, 0.6), linestyle = :dash)
MK.scatterlines!(axC, d["rel_p"], d["rel_emp"], color = :seagreen, markersize = 9, linewidth = 2)
MK.xlims!(axC, -0.02, 1.02); MK.ylims!(axC, -0.02, 1.02)

# (D) KDE pathology ------------------------------------------------------------------------
axD = MK.Axis(fig[2, 2],
              title = "(D) KDE baseline pathology (same pairs): $(round(100*frac_nonfinite))% dropped",
              xlabel = "|Δρ|", ylabel = "count")
absd = abs.(drho)
MK.hist!(axD, absd[finite], bins = 20, color = (:steelblue, 0.7), label = "KDE finite (kept)")
MK.hist!(axD, absd[.!finite], bins = 20, color = (:firebrick, 0.7),
         label = "KDE non-finite (dropped)")
MK.axislegend(axD; position = :rt, framevisible = false)

MK.Label(fig[0, :],
         "Spike 014 — simulation-based BF validation: NRE keeps 100% of pairs (AUC $(round(auc_nre; digits=3))); " *
         "KDE drops $(round(100*frac_nonfinite))% / saturates $(round(100*frac_saturated))% on the SAME pairs",
         fontsize = 15, font = :bold)

out = joinpath(HERE, "sim_bf_validation.png")
MK.save(out, fig)
println("saved: ", out)
