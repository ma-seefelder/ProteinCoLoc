using CairoMakie, JLD2, Statistics, Printf
const HERE = @__DIR__
d = JLD2.load(joinpath(HERE,"power_data.jld2"))
ss, Ms, curves, s50 = d["ss"], d["Ms"], d["curves"], d["s50"]
cols = Dict(250=>RGBf(0.6,0.6,0.6),500=>RGBf(0.30,0.60,0.80),
            1000=>RGBf(0.90,0.60,0.10),2000=>RGBf(0.80,0.25,0.20))
with_theme(theme_minimal()) do
    fig = Figure(size=(1100,640), fontsize=14)
    Label(fig[0,1:2], "Spike 011 — M=2000 SBC KS resolves a ~0.06-SD marginal drift; the residual nuisance drift sits right at that edge",
          fontsize=15, font=:bold, halign=:left)
    ax = Axis(fig[1,1], xlabel="marginal drift  s = mean(u) − 0.5", ylabel="KS rejection power",
              title="A  Rejection power vs drift, by M")
    for M in Ms
        lines!(ax, ss, curves[M], color=cols[M], linewidth=2.6, label="M=$M")
        scatter!(ax, ss, curves[M], color=cols[M], markersize=5)
    end
    hlines!(ax, [0.05], color=:black, linestyle=:dot, linewidth=1)
    hlines!(ax, [0.5], color=:gray, linestyle=:dash, linewidth=1)
    text!(ax, 0.002, 0.53, text="50% power", fontsize=10, color=:gray40)
    vspan!(ax, 0.011, 0.017, color=(RGBf(0.80,0.25,0.20),0.10))
    text!(ax, 0.014, 0.10, text="observed\nnuisance drift", fontsize=10,
          color=RGBf(0.6,0.2,0.15), align=(:center,:bottom))
    axislegend(ax, position=:rb, framevisible=false)

    ax2 = Axis(fig[1,2], xlabel="posterior-SD of marginal drift at 50% power", ylabel="M",
               title="B  Smallest drift each M resolves")
    barplot!(ax2, 1:length(Ms), s50, color=[cols[M] for M in Ms], width=0.6, direction=:x)
    ax2.yticks = (1:length(Ms), string.(Ms))
    for (i,v) in enumerate(s50)
        text!(ax2, v, i, text=@sprintf("%.3f", v), align=(:left,:center), offset=(4,0), fontsize=11)
    end
    vspan!(ax2, 0.04, 0.07, color=(RGBf(0.6,0.2,0.15),0.12))
    text!(ax2, 0.055, length(Ms)+0.4, text="observed flow error", fontsize=10,
          color=RGBf(0.6,0.2,0.15), align=(:center,:bottom))
    xlims!(ax2, 0, 0.20)
    Label(fig[2,1:2],
        "Residual drift (0.04–0.07 SD) is a location misplacement of the learned marginal on NON-identified nuisances. " *
        "M=2000 rejects it at ~50% power; M=500 would not. Over-power on nuisances, not a pathology on the targets ρ_true/Δρ.",
        fontsize=11.5, color=:gray35, halign=:left)
    save(joinpath(HERE,"power.png"), fig, px_per_unit=2)
    println("wrote power.png")
end
