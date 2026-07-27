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

# spike/validation/run_p11_paired_ridge.jl --- Phase-11 redesign feasibility, RECOVERABILITY LEG.
#
# The second leg of the pre-committed verdict in `run_p11_paired.jl`. Leg one (SNR against the
# per-draw noise floor) came back 0.755 (x) / 1.227 (y) against 0.406 for the current summary --
# a real 2-3x improvement that straddles the SNR > 1 bar. This leg asks the recoverability
# question with the SAME ridge harness, split discipline and positive control that produced the
# 0.998 null for the current summary, so the two are directly comparable.
#
# WHAT IS DIFFERENT FROM `run_p11_recovery.jl`, AND WHAT IS DELIBERATELY IDENTICAL.
# Different: the predictors. Instead of the 128 raw summary rows, the features are read off the
# correlation-vs-probe-offset curve -- peak location per axis (itself a direct estimate of the
# misregistration), peak sharpness, and peak height relative to the curve's baseline.
# Identical, on purpose: the ridge solver, the train-only standardisation, the validation-selected
# penalty, the held-out evaluation, the conditional-prior baseline (predicting 0, whose RMSE inside
# a lambda bin is computed empirically), the per-rung breakdown, and the rho_true positive control.
# LAMBDA IS STILL EXCLUDED FROM THE PREDICTORS, for the same prior-echo reason.
#
# COST NOTE. Probe offsets step at 0.5 px here rather than the 0.25 px used for the SNR leg. That
# is not a corner cut: the feature's OWN per-draw precision was measured at 2.0-2.8 px, so 0.25 px
# resolution is an order of magnitude finer than the quantity being estimated and buys nothing but
# wall clock. The offset span is unchanged.
#
# Research lane only. Changes no shipped contract, touches no `src/` file, retrains nothing,
# alters no constant, spends no iteration allowance.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p11_paired_ridge.jl

using JLD2
using Dates
using Statistics
using LinearAlgebra

const P11_PAIRED_LOAD_ONLY = true
include(joinpath(@__DIR__, "run_p11_paired.jl"))
isdefined(@__MODULE__, :generate_p11_pool) || include(joinpath(@__DIR__, "..", "data", "p11_generate.jl"))

const PAIRED_RIDGE_REPORT_PATH = joinpath(@__DIR__, "p11_paired_ridge_report.jld2")

const RIDGE_N_SAMPLES  = 1200
const RIDGE_OFFSETS    = collect(-4.0:0.5:4.0)
const RIDGE_GRID       = [1e-6, 1e-4, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0]
const RIDGE_TEST_FRAC  = 0.25
const RIDGE_IMSIZE     = P11_PROBE_COMPARABILITY_IMSIZE

"""
    paired_features(ch1, ch2) -> Vector{Float64}

The candidate feature vector read off the two 1-D probe curves. Peak LOCATION is the feature
expected to carry the signal -- it is a direct estimate of the misregistration. Sharpness and
relative height are included because a flat or shallow peak signals an unreliable location
estimate, which a downstream estimator can in principle exploit.
"""
function paired_features(ch1::Matrix{Float64}, ch2::Matrix{Float64})
    cx, cy = _probe_curve(ch1, ch2, RIDGE_OFFSETS)
    lx, vx, kx = _subpixel_peak(RIDGE_OFFSETS, cx)
    ly, vy, ky = _subpixel_peak(RIDGE_OFFSETS, cy)
    return [lx, ly, abs(kx), abs(ky), vx - minimum(cx), vy - minimum(cy),
            vx, vy, maximum(cx) - minimum(cx), maximum(cy) - minimum(cy)]
end

_fit_std(X) = (vec(mean(X; dims = 2)),
               [(isfinite(s) && s > 1e-12) ? s : 1.0 for s in vec(std(X; dims = 2))])
_apply_std(X, mu, sd) = (X .- mu) ./ sd

function ridge_fit(Xs::AbstractMatrix, y::AbstractVector, penalty::Real)
    b  = mean(y)
    w  = Symmetric(Xs * Xs' + penalty * I) \ (Xs * (y .- b))
    return w, b
end
ridge_predict(Xs, w, b) = (Xs' * w) .+ b
_rmse(p, t) = sqrt(mean(abs2, p .- t))

function main()
    println("="^78)
    println("PHASE-11 PAIRED-FEATURE RECOVERABILITY  (lambda still EXCLUDED from predictors)")
    println("="^78)
    t0 = time()

    # Re-simulate a subset from the SAME keyed stream the pool uses, so these samples are drawn
    # from the identical joint. The pool stores summaries, not images, and the probe features need
    # the images -- but the Philox-per-index keying makes re-simulation exact, not approximate.
    println("simulating $RIDGE_N_SAMPLES samples for feature extraction...")
    F  = zeros(Float64, 10, RIDGE_N_SAMPLES)
    DX = zeros(Float64, RIDGE_N_SAMPLES)
    DY = zeros(Float64, RIDGE_N_SAMPLES)
    LA = zeros(Float64, RIDGE_N_SAMPLES)
    RH = zeros(Float64, RIDGE_N_SAMPLES)

    done = Threads.Atomic{Int}(0)
    Threads.@threads for i in 1:RIDGE_N_SAMPLES
        rng = p11_datagen_rng(i)
        θ, lam = sample_p11_theta(rng)
        pair = simulate_pair(p11_datagen_rng(i), θ; imsize = RIDGE_IMSIZE)
        F[:, i] = paired_features(pair[1], pair[2])
        DX[i], DY[i], LA[i], RH[i] = θ.shift_dx, θ.shift_dy, lam, θ.ρ_true
        n = Threads.atomic_add!(done, 1) + 1
        n % 200 == 0 && println("  $n/$RIDGE_N_SAMPLES ($(round((time()-t0)/60; digits=2)) min)")
    end

    ntest = round(Int, RIDGE_TEST_FRAC * RIDGE_N_SAMPLES)
    tst   = (RIDGE_N_SAMPLES - ntest + 1):RIDGE_N_SAMPLES
    trall = 1:(RIDGE_N_SAMPLES - ntest)
    nval  = round(Int, 0.2 * length(trall))
    val   = (length(trall) - nval + 1):length(trall)
    fit   = 1:(length(trall) - nval)
    println("split: fit = $(length(fit))  val = $(length(val))  test = $(length(tst))")

    mu, sd = _fit_std(view(F, :, fit))
    Xf, Xv, Xt = _apply_std(view(F,:,fit),mu,sd), _apply_std(view(F,:,val),mu,sd), _apply_std(view(F,:,tst),mu,sd)

    results = Dict{String,Any}()
    for (name, Y) in (("shift_dx", DX), ("shift_dy", DY))
        yf, yv, yt = Y[fit], Y[val], Y[tst]
        bp, bv = RIDGE_GRID[1], Inf
        for p in RIDGE_GRID
            w, b = ridge_fit(Xf, yf, p)
            v = _rmse(ridge_predict(Xv, w, b), yv); v < bv && ((bp, bv) = (p, v))
        end
        w, b = ridge_fit(Xf, yf, bp)
        pred = ridge_predict(Xt, w, b)
        r_all, p_all = _rmse(pred, yt), sqrt(mean(abs2, yt))

        println("\n" * "-"^78)
        println("$name   penalty = $bp")
        println("  POOLED  ridge = $(round(r_all;digits=5))  prior = $(round(p_all;digits=5))  " *
                "ratio = $(round(r_all/p_all;digits=5))")
        println(rpad("lambda bin",16), rpad("n",7), rpad("ridge",13), rpad("prior",13), rpad("ratio",9))

        edges = collect(Float64.(SC2_RUNGS)); lam_t = LA[tst]
        bl, bh, bn, br, bpp = Float64[], Float64[], Int[], Float64[], Float64[]
        for k in 1:length(edges)
            lo = k == 1 ? LAMBDA_MIN : edges[k-1]; hi = edges[k]
            sel = findall(x -> (x >= lo) & (x < hi + (k == length(edges) ? 1e-9 : 0.0)), lam_t)
            isempty(sel) && continue
            r = _rmse(pred[sel], yt[sel]); p0 = sqrt(mean(abs2, yt[sel]))
            push!(bl,lo); push!(bh,hi); push!(bn,length(sel)); push!(br,r); push!(bpp,p0)
            println(rpad("[$(round(lo;digits=2)), $(round(hi;digits=2)))",16), rpad(length(sel),7),
                    rpad(round(r;digits=5),13), rpad(round(p0;digits=5),13),
                    rpad(round(r/p0;digits=5),9))
        end
        results[name] = Dict("penalty"=>bp, "pooled_ridge"=>r_all, "pooled_prior"=>p_all,
            "pooled_ratio"=>r_all/p_all, "bin_lo"=>bl, "bin_hi"=>bh, "bin_n"=>bn,
            "bin_ridge"=>br, "bin_prior"=>bpp, "bin_ratio"=>br ./ bpp)
    end

    # Positive control on the SAME feature set. NOTE: these 10 curve features are NOT expected to
    # recover rho_true as well as the full 128-row summary did (ratio 0.157) -- they are a
    # 10-number registration-oriented digest, not a colocalization summary. The control is here to
    # show the harness works on this feature set, not to claim parity.
    yf, yv, yt = RH[fit], RH[val], RH[tst]
    bp, bv = RIDGE_GRID[1], Inf
    for p in RIDGE_GRID
        w, b = ridge_fit(Xf, yf, p); v = _rmse(ridge_predict(Xv,w,b), yv)
        v < bv && ((bp, bv) = (p, v))
    end
    w, b = ridge_fit(Xf, yf, bp)
    cr, cp = _rmse(ridge_predict(Xt,w,b), yt), sqrt(mean(abs2, yt .- mean(yf)))
    println("\n" * "-"^78)
    println("POSITIVE CONTROL rho_true on the 10 curve features")
    println("  ridge = $(round(cr;digits=5))  prior = $(round(cp;digits=5))  " *
            "ratio = $(round(cr/cp;digits=5))")

    elapsed = (time() - t0) / 60
    tmp = PAIRED_RIDGE_REPORT_PATH * ".tmp"
    jldsave(tmp; schema_version=1, n_samples=RIDGE_N_SAMPLES, offsets=RIDGE_OFFSETS,
        n_features=10, lambda_excluded=true, shift_dx=results["shift_dx"],
        shift_dy=results["shift_dy"], control_rho_ratio=cr/cp,
        control_rho_ridge=cr, control_rho_prior=cp,
        baseline_current_summary_ratio=0.998, imsize=RIDGE_IMSIZE,
        elapsed_min=elapsed, generated=string(Dates.now(Dates.UTC))*"Z",
        caption="Phase-11 paired-feature recoverability: ridge on 10 probe-curve features " *
                "(lambda excluded) vs the conditional-prior baseline, per lambda rung. " *
                "Research lane; gates nothing.")
    let d = JLD2.load(tmp); @assert haskey(d, "shift_dx"); end
    mv(tmp, PAIRED_RIDGE_REPORT_PATH; force = true)
    println("\nwrote ", PAIRED_RIDGE_REPORT_PATH, "  (", round(elapsed; digits=2), " min)")
    return nothing
end

main()
