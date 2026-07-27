# ProteinCoLoc: Bayesian colocalization analysis of multi-channel fluorescence microscopy images
# Copyright (C) 2024  Manuel Seefelder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# spike/validation/run_p11_rhostrat.jl --- Phase-11: is the shift recoverable when COLOC IS STRONG?
#
# THE RATIONALE. Shift identifiability is structurally monotone in colocalization strength. When the
# two channels are independent, NO displacement makes them match better -- the correlation-vs-shift
# curve is flat and the shift is unidentifiable IN PRINCIPLE, not merely hard. When the channels
# resemble each other the curve has a genuine extremum and the shift is in principle findable.
# Every result so far pooled across the whole rho prior, so a uniform null could be masking a real
# dependence -- and the regime that matters practically, where colocalization actually exists, is
# exactly the strong end.
#
# ONE REFINEMENT TO THE STRATIFIER, AND IT IS NOT COSMETIC. The identifiability argument is about
# |rho|, NOT signed rho. At rho = -0.9 the channels are strongly ANTI-correlated: they are just as
# statistically dependent as at +0.9, and a displacement degrades that dependence just as
# systematically. It is rho = 0 -- not rho = -1 -- that is the unidentifiable case. The pool's rho
# prior spans [-0.99, +0.99] and is skewed negative (mean -0.106; only 10.8 % of samples exceed
# +0.5, but 30.1 % exceed |rho| = 0.5), so stratifying on SIGNED rho would put most of the
# genuinely-identifiable mass in the "low" bins and manufacture a null. Primary stratification is
# therefore on |rho|; signed bins are reported alongside so nothing is hidden.
#
# A CONSEQUENCE THAT AFFECTS CANDIDATE B'S VERDICT. `_subpixel_peak` takes an ARGMAX. For an
# anti-correlated pair the correlation-vs-offset curve has a MINIMUM at best alignment, so argmax
# locks onto noise for roughly half the pool. The paired feature was therefore tested under a
# handicap on the negative-rho half. This script re-extracts it BOTH ways -- as originally tested,
# and with a sign-aware extremum that tracks the correct direction -- so the verdict does not rest
# on that handicap. Both are reported.
#
# RE-ANALYSIS ONLY: no new pool, no new training, no threshold touched. The paired features are
# re-extracted (0.8 min) because the earlier run did not persist them, not because anything is
# being re-simulated differently.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p11_rhostrat.jl

using JLD2
using Dates
using Statistics
using LinearAlgebra

const P11_PAIRED_LOAD_ONLY = true
include(joinpath(@__DIR__, "run_p11_paired.jl"))
isdefined(@__MODULE__, :generate_p11_pool) || include(joinpath(@__DIR__, "..", "data", "p11_generate.jl"))

const RHOSTRAT_REPORT_PATH = joinpath(@__DIR__, "p11_rhostrat_report.jld2")

const RS_RIDGE_GRID = [1e-6, 1e-4, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0]
const RS_TEST_FRAC  = 0.25
const RS_N_PAIRED   = 6000
const RS_OFFSETS    = collect(-4.0:0.5:4.0)
const RS_IMSIZE     = P11_PROBE_COMPARABILITY_IMSIZE
const RS_ABS_EDGES  = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
const RS_MIN_N      = 100    # bins thinner than this are reported but NOT interpreted

_fit_std(X) = (vec(mean(X; dims = 2)),
               [(isfinite(s) && s > 1e-12) ? s : 1.0 for s in vec(std(X; dims = 2))])
_apply_std(X, mu, sd) = (X .- mu) ./ sd
function ridge_fit(Xs, y, p)
    b = mean(y); w = Symmetric(Xs * Xs' + p * I) \ (Xs * (y .- b)); return w, b
end
ridge_predict(Xs, w, b) = (Xs' * w) .+ b
_rmse(p, t) = sqrt(mean(abs2, p .- t))

"Fit ridge with the penalty chosen on validation, return held-out predictions."
function _fit_predict(Xf, yf, Xv, yv, Xt)
    bp, bv = RS_RIDGE_GRID[1], Inf
    for p in RS_RIDGE_GRID
        w, b = ridge_fit(Xf, yf, p)
        v = _rmse(ridge_predict(Xv, w, b), yv); v < bv && ((bp, bv) = (p, v))
    end
    w, b = ridge_fit(Xf, yf, bp)
    return ridge_predict(Xt, w, b), bp
end

"""
    _strat_table(label, pred, truth, strat, edges) -> NamedTuple

Ridge-vs-prior ratio per stratum. The prior baseline is the conditional prior mean (0), whose RMSE
inside a bin is computed EMPIRICALLY from the realised shifts there. Bins with fewer than
`RS_MIN_N` samples are reported with their n but flagged UNINTERPRETABLE rather than quoted.
"""
function _strat_table(label, pred, truth, strat, edges)
    println("\n  $label")
    println("  ", rpad("bin", 14), rpad("n", 8), rpad("ridge", 12), rpad("prior", 12),
            rpad("ratio", 10), "note")
    los, his, ns, rs, ps = Float64[], Float64[], Int[], Float64[], Float64[]
    for k in 1:(length(edges) - 1)
        lo, hi = edges[k], edges[k+1]
        sel = findall(x -> (x >= lo) & (x < hi + (k == length(edges)-1 ? 1e-9 : 0.0)), strat)
        n = length(sel)
        push!(los, lo); push!(his, hi); push!(ns, n)
        if n == 0
            push!(rs, NaN); push!(ps, NaN)
            println("  ", rpad("[$lo, $hi)", 14), rpad(0, 8), "-- empty --"); continue
        end
        r  = _rmse(pred[sel], truth[sel]); p0 = sqrt(mean(abs2, truth[sel]))
        push!(rs, r); push!(ps, p0)
        note = n < RS_MIN_N ? "UNINTERPRETABLE (n < $RS_MIN_N)" : ""
        println("  ", rpad("[$lo, $hi)", 14), rpad(n, 8), rpad(round(r; digits=5), 12),
                rpad(round(p0; digits=5), 12), rpad(round(r/p0; digits=5), 10), note)
    end
    return (lo = los, hi = his, n = ns, ridge = rs, prior = ps, ratio = rs ./ ps)
end

function main()
    println("="^78)
    println("PHASE-11 rho-STRATIFIED RECOVERABILITY  (re-analysis; no new pool, no training)")
    println("="^78)
    t0 = time()
    out = Dict{String,Any}()

    # =========================================================================================
    # PART 1 -- the ORIGINAL 128-row summary, stratified
    # =========================================================================================
    println("\nPART 1  128-row summary")
    pool = load_p11_pool(p11_pool_dir(P11_N_PAIRS))
    S  = Float64.(pool[:summary_min]); TH = Float64.(pool[:theta]); N = size(S, 2)
    RHO = TH[1, :]; DX = TH[5, :]; DY = TH[6, :]

    # CONFOUND CHECK, run BEFORE interpreting any gradient: does |rho| co-vary with the other
    # generative quantities, so that a rho gradient could actually be carried by something else?
    absr = abs.(RHO)
    names = ["rho_true","spillover","autofluor","label_eff","shift_dx","shift_dy","noise","chrom_eps"]
    println("  confound check -- cor(|rho|, .) across the pool:")
    confound = Float64[]
    for j in 2:8
        c = cor(absr, TH[j, :]); push!(confound, c)
        println("    ", rpad(names[j], 12), round(c; digits = 4))
    end
    c_lam = cor(absr, Float64.(pool[:lambda])); push!(confound, c_lam)
    println("    ", rpad("lambda", 12), round(c_lam; digits = 4))

    ntest = round(Int, RS_TEST_FRAC * N); tst = (N-ntest+1):N
    trall = 1:(N-ntest); nval = round(Int, 0.2*length(trall))
    val = (length(trall)-nval+1):length(trall); fit = 1:(length(trall)-nval)
    mu, sd = _fit_std(view(S,:,fit))
    Xf, Xv, Xt = _apply_std(view(S,:,fit),mu,sd), _apply_std(view(S,:,val),mu,sd), _apply_std(view(S,:,tst),mu,sd)

    pdx, _ = _fit_predict(Xf, DX[fit], Xv, DX[val], Xt)
    pdy, _ = _fit_predict(Xf, DY[fit], Xv, DY[val], Xt)
    out["summary_dx_abs"] = _strat_table("shift_dx by |rho|", pdx, DX[tst], absr[tst], RS_ABS_EDGES)
    out["summary_dy_abs"] = _strat_table("shift_dy by |rho|", pdy, DY[tst], absr[tst], RS_ABS_EDGES)
    out["summary_dx_signed"] = _strat_table("shift_dx by SIGNED rho", pdx, DX[tst], RHO[tst],
                                            [-1.0,-0.6,-0.2,0.2,0.6,1.0])

    prho, _ = _fit_predict(Xf, RHO[fit], Xv, RHO[val], Xt)
    ctrl = _rmse(prho, RHO[tst]) / sqrt(mean(abs2, RHO[tst] .- mean(RHO[fit])))
    println("\n  POSITIVE CONTROL rho_true (128-row): ratio = ", round(ctrl; digits = 5))

    # =========================================================================================
    # PART 2 -- the PAIRED probe-curve features, stratified, BOTH extremum conventions
    # =========================================================================================
    println("\nPART 2  paired probe-curve features (re-extracted, $RS_N_PAIRED samples)")
    Fa = zeros(Float64, 10, RS_N_PAIRED)   # as originally tested: argmax
    Fs = zeros(Float64, 10, RS_N_PAIRED)   # sign-aware extremum
    pDX = zeros(Float64, RS_N_PAIRED); pDY = zeros(Float64, RS_N_PAIRED)
    pRH = zeros(Float64, RS_N_PAIRED)
    done = Threads.Atomic{Int}(0)
    Threads.@threads for i in 1:RS_N_PAIRED
        rng = p11_datagen_rng(i)
        θ, _ = sample_p11_theta(rng)
        pr = simulate_pair(p11_datagen_rng(i), θ; imsize = RS_IMSIZE)
        cx, cy = _probe_curve(pr[1], pr[2], RS_OFFSETS)
        lx, vx, kx = _subpixel_peak(RS_OFFSETS, cx)
        ly, vy, ky = _subpixel_peak(RS_OFFSETS, cy)
        Fa[:, i] = [lx, ly, abs(kx), abs(ky), vx-minimum(cx), vy-minimum(cy), vx, vy,
                    maximum(cx)-minimum(cx), maximum(cy)-minimum(cy)]
        # SIGN-AWARE: for an anti-correlated pair the best alignment is the curve's MINIMUM, so
        # track the extremum in the direction the pair's own overall correlation dictates.
        sx = mean(cx) < 0 ? -1.0 : 1.0
        sy = mean(cy) < 0 ? -1.0 : 1.0
        lxs, vxs, kxs = _subpixel_peak(RS_OFFSETS, sx .* cx)
        lys, vys, kys = _subpixel_peak(RS_OFFSETS, sy .* cy)
        Fs[:, i] = [lxs, lys, abs(kxs), abs(kys), vxs-minimum(sx.*cx), vys-minimum(sy.*cy),
                    vxs, vys, maximum(sx.*cx)-minimum(sx.*cx), maximum(sy.*cy)-minimum(sy.*cy)]
        pDX[i], pDY[i], pRH[i] = θ.shift_dx, θ.shift_dy, θ.ρ_true
        n = Threads.atomic_add!(done,1)+1
        n % 400 == 0 && println("    $n/$RS_N_PAIRED")
    end

    pabs = abs.(pRH)
    nt = round(Int, RS_TEST_FRAC*RS_N_PAIRED); pt = (RS_N_PAIRED-nt+1):RS_N_PAIRED
    pa = 1:(RS_N_PAIRED-nt); nv = round(Int, 0.2*length(pa))
    pv = (length(pa)-nv+1):length(pa); pf = 1:(length(pa)-nv)

    for (tag, FF) in (("argmax (as tested)", Fa), ("sign-aware extremum", Fs))
        m2, s2 = _fit_std(view(FF,:,pf))
        Af, Av, At = _apply_std(view(FF,:,pf),m2,s2), _apply_std(view(FF,:,pv),m2,s2), _apply_std(view(FF,:,pt),m2,s2)
        qdx, _ = _fit_predict(Af, pDX[pf], Av, pDX[pv], At)
        qdy, _ = _fit_predict(Af, pDY[pf], Av, pDY[pv], At)
        println("\n  --- paired features, $tag ---")
        println("  POOLED  dx ratio = ",
                round(_rmse(qdx,pDX[pt])/sqrt(mean(abs2,pDX[pt])); digits=5),
                "   dy ratio = ", round(_rmse(qdy,pDY[pt])/sqrt(mean(abs2,pDY[pt])); digits=5))
        key = tag == "argmax (as tested)" ? "argmax" : "signaware"
        out["paired_dx_$key"] = _strat_table("shift_dx by |rho| ($tag)", qdx, pDX[pt], pabs[pt], RS_ABS_EDGES)
        out["paired_dy_$key"] = _strat_table("shift_dy by |rho| ($tag)", qdy, pDY[pt], pabs[pt], RS_ABS_EDGES)
        qr, _ = _fit_predict(Af, pRH[pf], Av, pRH[pv], At)
        println("  POSITIVE CONTROL rho_true ($tag): ratio = ",
                round(_rmse(qr,pRH[pt])/sqrt(mean(abs2,pRH[pt] .- mean(pRH[pf]))); digits=5))
    end

    elapsed = (time()-t0)/60
    tmp = RHOSTRAT_REPORT_PATH * ".tmp"
    jldsave(tmp; schema_version=1, abs_edges=RS_ABS_EDGES, min_n=RS_MIN_N,
        confound_cor=confound, confound_names=vcat(names[2:8], "lambda"),
        summary_dx_abs_ratio=out["summary_dx_abs"].ratio, summary_dx_abs_n=out["summary_dx_abs"].n,
        summary_dy_abs_ratio=out["summary_dy_abs"].ratio, summary_dy_abs_n=out["summary_dy_abs"].n,
        summary_dx_signed_ratio=out["summary_dx_signed"].ratio, summary_dx_signed_n=out["summary_dx_signed"].n,
        summary_control_ratio=ctrl,
        paired_dx_argmax_ratio=out["paired_dx_argmax"].ratio, paired_dx_argmax_n=out["paired_dx_argmax"].n,
        paired_dx_signaware_ratio=out["paired_dx_signaware"].ratio,
        paired_dy_signaware_ratio=out["paired_dy_signaware"].ratio,
        n_paired=RS_N_PAIRED, elapsed_min=elapsed,
        generated=string(Dates.now(Dates.UTC))*"Z",
        caption="Phase-11 rho-stratified recoverability. Primary stratifier |rho| (identifiability " *
                "is about dependence strength, not sign). Paired features re-extracted under both " *
                "argmax and a sign-aware extremum. Re-analysis only; gates nothing.")
    let d = JLD2.load(tmp); @assert haskey(d,"summary_dx_abs_ratio"); end
    mv(tmp, RHOSTRAT_REPORT_PATH; force=true)
    println("\nwrote ", RHOSTRAT_REPORT_PATH, "  (", round(elapsed;digits=2), " min)")
    return nothing
end

main()
