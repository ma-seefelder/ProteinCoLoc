#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder
E-Mail: manuel.seefelder@uni-ulm.de
Postal address: Department of Gene Therapy, University of Ulm, Helmholzstr. 8/1, 89081 Ulm, Germany

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=#

# spike/validation/run_bf.jl --- Phase-5 REPORTED amortized Bayes factor (BF-01, BF-02).
#
# THE ACTUAL BF REPRODUCTION PROOF. Runs the amortized log Bayes factor over the FULL
# pre-registered Δρ sweep (BF_SWEEP_N=25 in [BF_SWEEP_LO, BF_SWEEP_HI]) on the FRESH
# RESERVED stream val_rng(VAL_MASTER_SEED), and gates it against the src/bayes.jl KDE
# baseline on BOTH pre-registered criteria (D-08): corr ≥ BF_CORR_MIN AND
# max|Δ log-BF| ≤ BF_LOGBF_TOL. It is a REAL GATE: it EXITS NONZERO if either fails.
#
# AMORTIZED SIDE = logratio ONLY (BF-02): the amortized log-BF is read in ONE forward
# pass via amortized_log_bf (logratio(m=1) − logratio(m=0) − measured log_prior_odds).
# NO KernelDensity / QuadGK / contrastive-shuffle appears on this side. The KDE baseline
# (the reproduction TARGET) is computed in the ISOLATED spike/baseline/ env by
# run_bf_baseline.jl (path b) on the SAME persisted Δρ draws (apples-to-apples, D-07);
# this script shells out to that env and never imports KDE/QuadGK into the lean spike env.
#
# ANTI-SNOOPING / HONESTY CONTRACT (D-08, CLAUDE.md HARD_GATE): BF_CORR_MIN/BF_LOGBF_TOL/
# BF_SWEEP_* were committed to consts.jl BEFORE this run. 05-02 reported an HONEST fixture
# finding — order-correct (corr=0.961) but magnitude-inflated (max|Δ|=3.18 > 0.5). Whether
# the reported-scale run clears D-08b is an honest empirical question; a shortfall is a
# Phase-6 finding, NEVER a reason to loosen a locked threshold or reseed.
#
# CPU-only (D-10). Run:  julia --project=spike spike/validation/run_bf.jl

using Test
using JLD2
using Statistics
using Dates

# Units under test: the amortized BF surface (pulls bf.jl → train_ratio.jl → harness.jl)
# and the figure surface (CairoMakie — figures are script artifacts, not gate assertions).
isdefined(@__MODULE__, :amortized_log_bf)  || include(joinpath(@__DIR__, "bf.jl"))
isdefined(@__MODULE__, :plot_bf_agreement) || include(joinpath(@__DIR__, "figures.jl"))

const BF_REPORT_PATH   = joinpath(@__DIR__, "bf_report.jld2")
const RATIO_NET_PATH    = joinpath(@__DIR__, "trained_ratio.jld2")
const BASELINE_PROJECT = joinpath(@__DIR__, "..", "baseline")
const BASELINE_SCRIPT  = joinpath(@__DIR__, "..", "baseline", "run_bf_baseline.jl")

# --- Atomic persistence of the reported BF artifact (save_npe idiom, T-05-13) ----------
function _save_bf_report(path; kwargs...)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "logbf_amortized") "_save_bf_report: integrity check failed ($tmp)"
    end
    mv(tmp, path; force = true)
    return path
end

function main()
    println("="^78)
    println("Phase-5 REPORTED amortized Bayes factor (BF-01/02) — LOCKED consts, RESERVED stream")
    println("  sweep: N=$BF_SWEEP_N  Δρ∈[$BF_SWEEP_LO, $BF_SWEEP_HI]  imsize=$SBC_IMSIZE")
    println("  gate:  corr ≥ BF_CORR_MIN=$BF_CORR_MIN  AND  max|ΔlogBF| ≤ BF_LOGBF_TOL=$BF_LOGBF_TOL")
    println("  stream=VAL_MASTER_SEED=$(repr(VAL_MASTER_SEED))")
    println("="^78)
    @assert VAL_MASTER_SEED != NPE_MASTER_SEED "reported stream must be disjoint from training stream (D-02)"
    @assert VAL_MASTER_SEED != VAL_FIX_SEED    "reported stream must be disjoint from fixture stream (D-02)"

    m   = load_frozen_model()
    rat = load_ratio(RATIO_NET_PATH)                 # frozen ratio net + MEASURED log_prior_odds
    println("loaded ratio net: log_prior_odds = $(round(rat.log_prior_odds; sigdigits=4)) (measured, Pitfall 5)")

    # --- (1) Amortized side: the FULL Δρ sweep, one forward pass per point (BF-01/BF-02).
    #     Persists the Δρ posterior + prior draws so the ISOLATED baseline scores IDENTICAL
    #     inputs (apples-to-apples, D-07).
    println("\n[1/3] generating reported amortized sweep (BF_SWEEP_N=$BF_SWEEP_N) …")
    t0 = time()
    sw = generate_bf_sweep(m, rat.estimator, rat.log_prior_odds;
                           n = BF_SWEEP_N, lo = BF_SWEEP_LO, hi = BF_SWEEP_HI,
                           rng = val_rng(VAL_MASTER_SEED), imsize = SBC_IMSIZE,
                           draws_path = BF_SWEEP_DRAWS_PATH)
    println("      amortized sweep done in $(round((time()-t0)/60; digits=2)) min; draws -> $BF_SWEEP_DRAWS_PATH")

    # --- (2) KDE baseline (path b): run_bf_baseline.jl in the ISOLATED spike/baseline/ env,
    #     scoring the SAME Δρ draws (the reproduction TARGET; the ONLY side with KDE/QuadGK).
    println("\n[2/3] running KDE baseline in the isolated spike/baseline/ env …")
    run(`julia --project=$BASELINE_PROJECT $BASELINE_SCRIPT`)

    # --- (3) The D-08 reproduction: corr + max|Δ log-BF| vs the KDE baseline ------------
    println("\n[3/3] scoring the reproduction …")
    rep = bf_reproduction(m, rat.estimator, rat.log_prior_odds; regenerate = false)

    # --- Persist the reported artifact BEFORE the gate (survives an honest fail) --------
    _save_bf_report(BF_REPORT_PATH;
        targets = collect(rep.targets),
        logbf_amortized = collect(rep.logbf_amortized),
        logbf_kde = collect(rep.logbf_kde),
        corr = rep.corr, max_abs_err = rep.max_abs_err,
        log_prior_odds = rat.log_prior_odds,
        BF_CORR_MIN = BF_CORR_MIN, BF_LOGBF_TOL = BF_LOGBF_TOL,
        BF_SWEEP_LO = BF_SWEEP_LO, BF_SWEEP_HI = BF_SWEEP_HI, BF_SWEEP_N = BF_SWEEP_N,
        VAL_MASTER_SEED = UInt64(VAL_MASTER_SEED),
        generated = string(Dates.now(Dates.UTC)) * "Z")
    println("persisted reported artifact -> $BF_REPORT_PATH")

    # --- Figure (script artifact — never a gate assertion) ------------------------------
    figpath = plot_bf_agreement(rep.logbf_amortized, rep.logbf_kde; filename = "bf_agreement.png")
    println("figure -> $figpath")

    # --- Headline (printed BEFORE the gate throws) --------------------------------------
    corr_pass = rep.corr >= BF_CORR_MIN
    err_pass  = rep.max_abs_err <= BF_LOGBF_TOL
    println("\n", "-"^78)
    println("REPORTED BF reproduction over $(length(rep.targets)) Δρ points:")
    println("  corr(amortized, KDE) = $(round(rep.corr; sigdigits=4))   " *
            "(≥ BF_CORR_MIN=$BF_CORR_MIN ? ", corr_pass ? "PASS" : "FAIL", ")")
    println("  max|Δ log-BF|        = $(round(rep.max_abs_err; sigdigits=4))   " *
            "(≤ BF_LOGBF_TOL=$BF_LOGBF_TOL ? ", err_pass ? "PASS" : "FAIL", ")")
    println("Overall reported BF gate: ", (corr_pass && err_pass) ? "PASS" : "FAIL")
    println("-"^78, "\n")

    # --- THE REAL GATE (D-08): exits nonzero if either pre-registered criterion fails.
    # Thresholds referenced from consts.jl (never tuned). An honest failure → Phase-6.
    @testset "Reported BF pre-registered gate (BF-01/02, D-08)" begin
        @test rep.corr        >= BF_CORR_MIN      # D-08a ordering
        @test rep.max_abs_err <= BF_LOGBF_TOL     # D-08b bounded magnitude
    end
    return nothing
end

main()
