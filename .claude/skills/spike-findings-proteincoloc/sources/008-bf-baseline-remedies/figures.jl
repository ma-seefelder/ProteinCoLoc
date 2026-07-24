#########################################################################################
# Spike 007/008 --- figures: survivorship bias and the three baseline behaviours.
# Headless CairoMakie (.planning/spikes/CONVENTIONS.md): scoped with_theme, no global set_theme!.
#########################################################################################

using JLD2, CairoMakie, Statistics, Printf

const HERE  = @__DIR__
const S006  = normpath(joinpath(HERE, "..", "006-bf-attrition-mechanism"))
const S007  = normpath(joinpath(HERE, "..", "007-bf-survivorship-bias"))

rem  = JLD2.load(joinpath(HERE, "remedies_rows.jld2"))["rows"]
b006 = JLD2.load(joinpath(S006, "attrition_rows.jld2"))
r006 = b006["rows"]
lpo  = JLD2.load(joinpath(S006, "attrition_rows.jld2"))["p_prior"]
log_prior_odds = log(lpo / (1 - lpo))

a     = getfield.(rem, :a)
b_unc = getfield.(rem, :b_unc); b_clp = getfield.(rem, :b_clp); b_log = getfield.(rem, :b_log)
dt    = getfield.(rem, :dtrue)
fin   = isfinite.(b_unc)

const C_KEEP = RGBf(0.20, 0.49, 0.72)
const C_DROP = RGBf(0.85, 0.37, 0.01)
const C_ACC  = RGBf(0.30, 0.60, 0.30)

# resolution ceiling from L draws: smallest resolvable tail prob is 1/L
const L_DRAWS = 999
const RES_CEIL = log(L_DRAWS)

clamped_from_p(p, eps) = (pc = clamp(p, eps, 1 - eps); log(pc / (1 - pc)) - log_prior_odds)
p_post = getfield.(r006, :p_post)
keep6  = getfield.(r006, :a_finite) .& getfield.(r006, :b_finite)
a6     = getfield.(r006, :a)

with_theme(theme_minimal()) do
    fig = Figure(size = (1180, 940), fontsize = 14)
    Label(fig[0, 1:2],
          "Spikes 007/008 — the survivor filter biases r̂, and the clamp hides a real tail disagreement",
          fontsize = 17, font = :bold, halign = :left)

    # --- A: survivorship bias in r ------------------------------------------------------------
    axA = Axis(fig[1, 1], title = "A  r̂ on survivors vs on all pairs (clamped baseline)",
               xlabel = "clamp ε", ylabel = "corr(amortized, baseline)",
               xscale = log10, xticks = ([1e-15, 1e-12, 1e-10, 1e-8], ["1e-15","1e-12","1e-10","1e-8"]))
    epss = [1e-15, 1e-12, 1e-10, 1e-8]
    rs = [cor(a6[keep6], clamped_from_p.(p_post, e)[keep6]) for e in epss]
    ra = [cor(a6,        clamped_from_p.(p_post, e))        for e in epss]
    lines!(axA, epss, rs, color = C_KEEP, linewidth = 2.5, label = "survivors only (what the gate does)")
    scatter!(axA, epss, rs, color = C_KEEP, markersize = 11)
    lines!(axA, epss, ra, color = C_DROP, linewidth = 2.5, label = "all pairs")
    scatter!(axA, epss, ra, color = C_DROP, markersize = 11)
    hlines!(axA, [0.95], color = :black, linestyle = :dash, linewidth = 1)
    text!(axA, 1e-13, 0.952, text = "gate threshold 0.95", fontsize = 11, color = :gray30)
    axislegend(axA, position = :lb, framevisible = false, labelsize = 11)

    # --- B: what the survivors are ------------------------------------------------------------
    axB = Axis(fig[1, 2], title = "B  The survivor filter keeps the ambiguous middle",
               xlabel = "|Δρ_true|", ylabel = "pairs")
    ak = abs.(getfield.(r006, :dtrue))[keep6]
    ad = abs.(getfield.(r006, :dtrue))[.!keep6]
    hist!(axB, ak, bins = range(0, 1.2, length = 19), color = (C_KEEP, 0.75), label = "kept")
    hist!(axB, ad, bins = range(0, 1.2, length = 19), color = (C_DROP, 0.75), label = "dropped")
    axislegend(axB, position = :rt, framevisible = false)

    # --- C: the three baselines against the amortized estimator --------------------------------
    axC = Axis(fig[2, 1], title = "C  Baseline vs amortized log-BF (symlog y)",
               xlabel = "amortized log-BF (NRE)", ylabel = "baseline log-BF")
    sgn(y) = sign.(y) .* log10.(1 .+ abs.(y))
    scatter!(axC, a[fin],  sgn(b_unc[fin]), color = (C_KEEP, 0.85), markersize = 9,
             label = "unclamped (survives)")
    scatter!(axC, a[.!fin], sgn(b_clp[.!fin]), color = (C_DROP, 0.85), markersize = 9,
             marker = :rect, label = "clamped (dropped pairs) → ±18.42")
    scatter!(axC, a[.!fin], sgn(b_log[.!fin]), color = (C_ACC, 0.9), markersize = 9,
             marker = :diamond, label = "log-space (dropped pairs)")
    axislegend(axC, position = :lt, framevisible = false, labelsize = 10)
    axC.ylabel = "sign·log₁₀(1+|log-BF|)"

    # --- D: the resolution ceiling -------------------------------------------------------------
    axD = Axis(fig[2, 2], title = "D  |log-BF| vs what L=999 draws can resolve",
               xlabel = "|Δρ_true|", ylabel = "|log-BF| (log₁₀ scale)", yscale = log10)
    scatter!(axD, abs.(dt), max.(abs.(a), 1e-2), color = (C_KEEP, 0.8), markersize = 9,
             label = "amortized (NRE)")
    scatter!(axD, abs.(dt), max.(abs.(b_log), 1e-2), color = (C_ACC, 0.8), markersize = 9,
             marker = :diamond, label = "log-space KDE baseline")
    hlines!(axD, [RES_CEIL], color = :red, linestyle = :dash, linewidth = 2)
    text!(axD, 0.02, RES_CEIL * 1.25, text = "log(L)=6.9 — resolution ceiling for 999 draws",
          fontsize = 11, color = :red)
    axislegend(axD, position = :lt, framevisible = false, labelsize = 10)

    Label(fig[3, 1:2],
          "Anything above the red line is Gaussian-kernel tail EXTRAPOLATION, not information in the 999 draws. " *
          "DIAGNOSTIC — DEV seeds, not the confirmatory gate.",
          fontsize = 12, color = :gray40, halign = :left)

    save(joinpath(HERE, "remedies.png"), fig, px_per_unit = 2)
    println("wrote ", joinpath(HERE, "remedies.png"))
end
