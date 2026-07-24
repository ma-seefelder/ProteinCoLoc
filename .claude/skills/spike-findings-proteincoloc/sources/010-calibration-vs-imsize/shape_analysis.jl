# Spike 010d --- if the ranks are centred but still reject, the SHAPE must be wrong.
# U-shaped  => posterior TOO NARROW (overconfident); n-shaped => TOO WIDE (underconfident).
using JLD2, Printf, Statistics
d = JLD2.load(joinpath(@__DIR__, "mixture_replicate.jld2"))
R, L, M, LAB = d["ranks"], d["L"], d["M"], d["labels"]
println("Uniform expectation: tails(|u-0.5|>0.4)=0.20, middle(|u-0.5|<0.1)=0.20")
println("ratio = tails/middle;  >1.15 U-shaped (overconfident), <0.87 n-shaped (underconfident)\n")
@printf("%-18s %8s %8s %8s   %s\n","parameter","tails","middle","ratio","reading")
for p in 1:7
    u = (R[:,p] .+ 0.5) ./ (L+1)
    t = mean(abs.(u .- 0.5) .> 0.4); mid = mean(abs.(u .- 0.5) .< 0.1)
    r = t/mid
    rd = r > 1.15 ? "U -> POSTERIOR TOO NARROW (overconfident)" :
         r < 0.87 ? "n -> posterior too wide (underconfident)" : "flat"
    @printf("%-18s %8.3f %8.3f %8.2f   %s\n", LAB[p], t, mid, r, rd)
end
# binomial SE for the tail fraction at n=M
@printf("\n(SE on each fraction at M=%d is about %.3f)\n", M, sqrt(0.2*0.8/M))
