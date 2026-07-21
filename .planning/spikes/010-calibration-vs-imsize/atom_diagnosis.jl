using JLD2, Printf, Statistics
d = JLD2.load(".planning/spikes/010-calibration-vs-imsize/mixture_replicate.jld2")
R,L,M,th,psd = d["ranks"],d["L"],d["M"],d["theta"],d["psd"]
r = th[:,1]                          # theta*_rho_true
u = (R[:,1] .+ 0.5) ./ (L+1)
lo,hi = minimum(r), maximum(r)
@printf("rho_true* Spannweite: [%.4f, %.4f]\n", lo, hi)
for (nm,mask) in (("am unteren Rand (<lo+0.01)", r .< lo+0.01),
                  ("am oberen Rand (>hi-0.01)",  r .> hi-0.01),
                  ("INNEN (Rest)",               (r .>= lo+0.01) .& (r .<= hi-0.01)))
    n = count(mask); n==0 && continue
    uu = u[mask]
    @printf("%-26s n=%4d (%.1f%%)  Tails=%.3f  Mitte=%.3f  Verh.=%.2f\n",
            nm, n, 100n/M, mean(abs.(uu.-0.5).>0.4), mean(abs.(uu.-0.5).<0.1),
            mean(abs.(uu.-0.5).>0.4)/max(mean(abs.(uu.-0.5).<0.1),1e-9))
end
# Atome zaehlen: exakte Wiederholungen der Extremwerte
@printf("\nexakte Duplikate von min: %d, von max: %d  (Atome => >1)\n",
        count(==(lo), r), count(==(hi), r))
# Posterior-Breite in den Tail-Faellen vs Mitte (Ueberkonfidenz-Signatur)
tail = abs.(u .- 0.5) .> 0.4
@printf("median post_sd  Tail-Faelle %.4f   Mitte-Faelle %.4f\n",
        median(psd[tail,1]), median(psd[.!tail,1]))
