using JLD2, Printf, Statistics, HypothesisTests, Distributions, Random
d = JLD2.load(".planning/spikes/010-calibration-vs-imsize/mixture_replicate.jld2")
R,L,M,th,LAB = d["ranks"],d["L"],d["M"],d["theta"],d["labels"]
println("Atome pro Parameter (exakte Duplikate der Extremwerte):")
for p in 1:7
    x = th[:,p]; lo,hi = minimum(x),maximum(x)
    nlo,nhi = count(==(lo),x), count(==(hi),x)
    (nlo>1 || nhi>1) && @printf("  %-18s min-Atome=%3d  max-Atome=%3d  (%.1f%% der Ziehungen)\n",
                                LAB[p], nlo, nhi, 100*(nlo+nhi)/M)
end
println("\nKS-Test fuer rho_true:")
x = th[:,1]; lo,hi = minimum(x),maximum(x)
atom = (x .== lo) .| (x .== hi)
u = (R[:,1] .+ 0.5) ./ (L+1)
@printf("  ALLE Ziehungen        n=%4d  KS p = %.3e\n", M, pvalue(ExactOneSampleKSTest(u,Uniform(0,1))))
ui = u[.!atom]
@printf("  OHNE Atome (interior) n=%4d  KS p = %.3e   <-- gueltiger SBC-Test\n",
        length(ui), pvalue(ExactOneSampleKSTest(ui,Uniform(0,1))))
# randomisierter Rang (Standardkorrektur fuer gemischte Verteilungen)
Random.seed!(1234)
ur = copy(u); ur[atom] = rand(count(atom))
@printf("  randomisierte Raenge  n=%4d  KS p = %.3e\n", M, pvalue(ExactOneSampleKSTest(ur,Uniform(0,1))))
