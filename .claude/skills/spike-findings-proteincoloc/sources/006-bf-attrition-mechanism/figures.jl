#########################################################################################
# Spike 006 --- figures for the BF-arm attrition mechanism.
#
# Reads attrition_rows.jld2 (written by measure_attrition.jl) and renders a 4-panel PNG.
# Headless CairoMakie per .planning/spikes/CONVENTIONS.md (scoped with_theme, never a global
# set_theme!).
#########################################################################################

using JLD2, CairoMakie, Statistics, Printf

const HERE = @__DIR__
d      = JLD2.load(joinpath(HERE, "attrition_rows.jld2"))
rows   = d["rows"]
n      = length(rows)

keep = filter(r ->   r.a_finite && r.b_finite, rows)
drop = filter(r -> !(r.a_finite && r.b_finite), rows)

adrop = abs.(getfield.(drop, :dtrue))
akeep = abs.(getfield.(keep, :dtrue))

# a colourblind-safe pair
const C_KEEP = RGBf(0.20, 0.49, 0.72)   # blue
const C_DROP = RGBf(0.85, 0.37, 0.01)   # orange

with_theme(theme_minimal()) do
    fig = Figure(size = (1180, 900), fontsize = 14)

    Label(fig[0, 1:2],
          "Spike 006 — BF-arm attrition: the non-clamped KDE baseline drops the high-signal pairs",
          fontsize = 17, font = :bold, halign = :left)

    # --- Panel A: where does the non-finiteness come from? ---------------------------------
    axA = Axis(fig[1, 1], title = "A  Source of non-finite values",
               ylabel = "pairs", xticks = (1:3, ["amortized\n(NRE)", "baseline\n(KDE)", "both"]))
    cnt = [count(r -> !r.a_finite, rows), count(r -> !r.b_finite, rows),
           count(r -> !r.a_finite && !r.b_finite, rows)]
    barplot!(axA, 1:3, cnt, color = [C_KEEP, C_DROP, :gray70], width = 0.6)
    for (i, c) in enumerate(cnt)
        text!(axA, i, c, text = string(c), align = (:center, :bottom), offset = (0, 4))
    end
    ylims!(axA, 0, max(maximum(cnt), 1) * 1.25)

    # --- Panel B: |Delta rho_true| -- dropped vs kept ---------------------------------------
    axB = Axis(fig[1, 2], title = "B  True |Δρ| of dropped vs kept pairs",
               ylabel = "|Δρ_true|", xticks = (1:2, ["kept", "dropped"]))
    for (i, (v, c)) in enumerate(((akeep, C_KEEP), (adrop, C_DROP)))
        isempty(v) && continue
        boxplot!(axB, fill(i, length(v)), v, color = (c, 0.35), width = 0.5,
                 mediancolor = c, whiskerwidth = 0.4)
        scatter!(axB, fill(i, length(v)) .+ 0.16 .* (rand(length(v)) .- 0.5), v,
                 color = (c, 0.75), markersize = 7)
    end
    !isempty(akeep) && !isempty(adrop) && text!(axB, 1.5, max(maximum(akeep), maximum(adrop)),
        text = @sprintf("median %.2f  vs  %.2f", median(akeep), median(adrop)),
        align = (:center, :top), fontsize = 12, color = :gray30)

    # --- Panel C: p_post lands exactly on the boundary --------------------------------------
    axC = Axis(fig[2, 1], title = "C  KDE tail probability p_post (dropped pairs sit at 1.0)",
               xlabel = "p_post", ylabel = "pairs")
    pk = getfield.(keep, :p_post); pd = getfield.(drop, :p_post)
    !isempty(pk) && hist!(axC, pk, bins = range(0, 1, length = 26), color = (C_KEEP, 0.75),
                          label = "kept")
    !isempty(pd) && hist!(axC, pd, bins = range(0, 1, length = 26), color = (C_DROP, 0.85),
                          label = "dropped")
    vlines!(axC, [0.0, 1.0], color = :black, linestyle = :dash, linewidth = 1)
    axislegend(axC, position = :lt, framevisible = false)

    # --- Panel D: attrition vs |Delta rho| ---------------------------------------------------
    axD = Axis(fig[2, 2], title = "D  Attrition rate rises with true |Δρ|",
               xlabel = "|Δρ_true| bin", ylabel = "fraction dropped")
    edges = 0.0:0.2:1.2
    ctrs, fracs, ns = Float64[], Float64[], Int[]
    for i in 1:(length(edges) - 1)
        lo, hi = edges[i], edges[i + 1]
        sub = filter(r -> lo <= abs(r.dtrue) < hi, rows)
        isempty(sub) && continue
        push!(ctrs, (lo + hi) / 2)
        push!(fracs, count(r -> !(r.a_finite && r.b_finite), sub) / length(sub))
        push!(ns, length(sub))
    end
    barplot!(axD, ctrs, fracs, width = 0.17, color = C_DROP)
    for (c, f, k) in zip(ctrs, fracs, ns)
        text!(axD, c, f, text = "n=$k", align = (:center, :bottom), offset = (0, 3), fontsize = 10)
    end
    ylims!(axD, 0, 1.15)

    Label(fig[3, 1:2],
          @sprintf("n = %d attempted, %d dropped (%.0f%%). DIAGNOSTIC run on a DEV seed — not the confirmatory gate.",
                   n, length(drop), 100 * length(drop) / n),
          fontsize = 12, color = :gray40, halign = :left)

    save(joinpath(HERE, "attrition.png"), fig, px_per_unit = 2)
    println("wrote ", joinpath(HERE, "attrition.png"))
end
