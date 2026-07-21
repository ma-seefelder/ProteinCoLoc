using JLD2, Printf, Statistics
r = JLD2.load("artifacts/amended_v2/grid_8/gate_report_8.jld2")["report"]
R = r.sbc.ranks; L = r.sbc.L; M = size(R,1)
se = sqrt((1/12)/M)
labels = ["ρ_true","spillover","autofluor","label_eff","shift_dx","shift_dy","noise","Δρ"]
println("Parameter    mean(u)   z(Lage)   Shrinkage  KS p       -> Lage erklaert Ablehnung?")
for (j,pp) in enumerate(r.sbc.per_param)
    u = (R[:,j] .+ 0.5) ./ (L+1)
    z = (mean(u)-0.5)/se
    @printf("%-11s  %.4f   %+6.1f    %.3f      %.2e   %s\n",
        labels[j], mean(u), z, pp.shrinkage, pp.ks_p,
        abs(z)>3 ? "JA (Lageverschiebung)" : (pp.ks_p<0.05 ? "NEIN -> Form, nicht Lage" : "-"))
end
