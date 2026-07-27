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

# spike/validation/run_p11_recovery.jl --- Phase-11 shift-recoverability test (SC1g follow-up).
#
# THE QUESTION: can ANY estimator recover the sub-pixel shift from the 128-row summary?
#
# ROW 129 IS EXCLUDED FROM THE PREDICTORS, AND THIS IS THE WHOLE POINT OF THE TEST.
# Row 129 is `encode_lambda(lambda)` -- the uncertainty level handed to the network as an input.
# Because the training joint draws `shift ~ Uniform(-lambda, lambda)`, a predictor that sees
# lambda can reproduce the conditional prior spread perfectly WITHOUT extracting one bit of
# information from the image: the prior-SD ratio across the ladder is exactly
# LAMBDA_MAX / LAMBDA_MIN = 3.0 / 0.25 = 12.0 by construction, since Uniform(-lam, lam) has
# SD = lam / sqrt(3). Any diagnostic that admits row 129 therefore measures PRIOR ECHO, not
# learning, and is vacuous for this question. The pool stores the raw 128-row summary and lambda
# in separate fields, so the exclusion here is structural, not a filtering step that could slip.
#
# THE BASELINE: predicting the conditional prior mean (which is 0, by symmetry of
# Uniform(-lam, lam)). Its RMSE inside a lambda bin is sqrt(E[shift^2]) = sqrt(E[lambda^2] / 3).
# It is computed EMPIRICALLY from the realised shifts in each bin rather than from the analytic
# formula, so the comparison cannot be flattered by a mismatch between the nominal and realised
# lambda distribution. Ridge beating this baseline means the summary carries recoverable shift
# information; ridge merely matching it means it does not.
#
# Ridge is the right tool here precisely because it is weak and transparent: closed-form, no
# tuning theatre, no capacity to memorise. If a 128-predictor linear map cannot beat the prior,
# the honest reading is that the linear signal is absent -- and the companion script
# `run_p11_bottleneck.jl` independently measures WHY (the per-draw noise floor exceeds the whole
# effect). This is a DIAGNOSTIC, not a deliverable, and it gates nothing.
#
# DECOUPLING (S5): spike-local; reads the frozen pre-registration and the existing training pool;
# writes one new artifact. Trains no network, mutates no constant, touches no `src/` file.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p11_recovery.jl

using JLD2
using Dates
using Statistics
using LinearAlgebra

isdefined(@__MODULE__, :P11_DEV_SEED)        || include(joinpath(@__DIR__, "p11_consts.jl"))
isdefined(@__MODULE__, :generate_p11_pool)   || include(joinpath(@__DIR__, "..", "data", "p11_generate.jl"))

const RECOVERY_REPORT_PATH = joinpath(@__DIR__, "p11_recovery_report.jld2")

# Held-out fraction. A plain index split is sufficient and is the leak-free discipline that
# matters here: the transform and the ridge coefficients are fitted on TRAIN COLUMNS ONLY and
# applied to test, so no test-set statistic can leak into the fit.
const TEST_FRACTION = 0.25

# Ridge penalties swept on a validation slice carved out of TRAIN (never out of TEST).
const RIDGE_GRID = [1e-6, 1e-4, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0]

"""
    _fit_standardizer(X) -> (mu, sd)

Column-wise (per predictor row) mean and sd over TRAINING columns only. Rows with zero variance
(the near-constant present-mask block) get sd = 1 so they pass through as constants instead of
producing NaNs.
"""
function _fit_standardizer(X::AbstractMatrix)
    mu = vec(mean(X; dims = 2))
    sd = vec(std(X; dims = 2))
    @inbounds for i in eachindex(sd)
        (isfinite(sd[i]) && sd[i] > 1e-12) || (sd[i] = 1.0)
    end
    return mu, sd
end

_apply_standardizer(X, mu, sd) = (X .- mu) ./ sd

"""
    ridge_fit(Xs, y, penalty) -> (w, b)

Closed-form ridge on standardized predictors with an unpenalised intercept. `Xs` is
(features x samples). Solves (X X' + penalty I) w = X y_centered.
"""
function ridge_fit(Xs::AbstractMatrix, y::AbstractVector, penalty::Real)
    b  = mean(y)
    yc = y .- b
    G  = Symmetric(Xs * Xs' + penalty * I)
    w  = G \ (Xs * yc)
    return w, b
end

ridge_predict(Xs, w, b) = (Xs' * w) .+ b

_rmse(pred, truth) = sqrt(mean(abs2, pred .- truth))

function main()
    println("="^78)
    println("PHASE-11 SHIFT-RECOVERABILITY TEST  (row 129 EXCLUDED from predictors)")
    println("Ridge on the 128-row summary vs the conditional-prior baseline. Gates nothing.")
    println("="^78)
    t_start = time()

    dir  = p11_pool_dir(P11_N_PAIRS)
    pool = load_p11_pool(dir)
    S    = Float64.(pool[:summary_min])       # 128 x N -- the RAW summary; lambda is NOT in here
    TH   = Float64.(pool[:theta])             # 8 x N
    LAM  = Float64.(pool[:lambda])            # N
    N    = size(S, 2)

    @assert size(S, 1) == 128 "predictors must be exactly the 128 summary rows"
    @assert size(TH, 1) == P11_THETA_DIM
    @assert length(LAM) == N

    # Theta row indices for the shift, read from the frozen field order rather than hard-coded.
    fields = (:ρ_true, :spillover, :autofluorescence, :label_efficiency,
              :shift_dx, :shift_dy, :noise, :chromatic_eps)
    idx_dx = findfirst(==(:shift_dx), fields)
    idx_dy = findfirst(==(:shift_dy), fields)
    println("shift rows in theta: dx = $idx_dx, dy = $idx_dy   (N = $N)")

    # --- Leak-free split ---------------------------------------------------------------------
    ntest  = round(Int, TEST_FRACTION * N)
    test_i = (N - ntest + 1):N
    tr_all = 1:(N - ntest)
    nval   = round(Int, 0.2 * length(tr_all))
    val_i  = (length(tr_all) - nval + 1):length(tr_all)
    fit_i  = 1:(length(tr_all) - nval)
    println("split: fit = $(length(fit_i))  val = $(length(val_i))  test = $(length(test_i))")

    mu, sd = _fit_standardizer(view(S, :, fit_i))          # TRAIN-ONLY fit
    Xfit = _apply_standardizer(view(S, :, fit_i), mu, sd)
    Xval = _apply_standardizer(view(S, :, val_i), mu, sd)
    Xtst = _apply_standardizer(view(S, :, test_i), mu, sd)

    results = Dict{String,Any}()

    for (name, idx) in (("shift_dx", idx_dx), ("shift_dy", idx_dy))
        yfit = vec(TH[idx, fit_i])
        yval = vec(TH[idx, val_i])
        ytst = vec(TH[idx, test_i])

        # Select the penalty on VALIDATION, never on test.
        best_p, best_v = RIDGE_GRID[1], Inf
        for p in RIDGE_GRID
            w, b = ridge_fit(Xfit, yfit, p)
            v = _rmse(ridge_predict(Xval, w, b), yval)
            v < best_v && ((best_p, best_v) = (p, v))
        end
        w, b = ridge_fit(Xfit, yfit, best_p)
        pred = ridge_predict(Xtst, w, b)

        overall_ridge = _rmse(pred, ytst)
        overall_prior = sqrt(mean(abs2, ytst))     # predicting 0, the conditional prior mean

        println()
        println("-"^78)
        println("$name   best ridge penalty = $best_p  (chosen on validation)")
        println("-"^78)
        println("  POOLED   ridge RMSE = $(round(overall_ridge; digits = 5))   " *
                "prior RMSE = $(round(overall_prior; digits = 5))   " *
                "ratio = $(round(overall_ridge / overall_prior; digits = 5))")
        println()
        println(rpad("lambda bin", 16), rpad("n", 8), rpad("ridge RMSE", 14),
                rpad("prior RMSE", 14), rpad("ratio", 10))

        # Per-lambda-rung breakdown: a pooled number can hide a rung-dependent effect, so the
        # test set is binned by the sample's own lambda using the SC2 ladder edges.
        edges = collect(Float64.(SC2_RUNGS))
        lam_t = LAM[test_i]
        bin_lo, bin_hi, bin_n = Float64[], Float64[], Int[]
        bin_ridge, bin_prior  = Float64[], Float64[]
        for k in 1:(length(edges))
            lo = k == 1 ? LAMBDA_MIN : edges[k-1]
            hi = edges[k]
            sel = findall(x -> (x >= lo) & (x < hi + (k == length(edges) ? 1e-9 : 0.0)), lam_t)
            isempty(sel) && continue
            r = _rmse(pred[sel], ytst[sel])
            p0 = sqrt(mean(abs2, ytst[sel]))
            push!(bin_lo, lo); push!(bin_hi, hi); push!(bin_n, length(sel))
            push!(bin_ridge, r); push!(bin_prior, p0)
            println(rpad("[$(round(lo;digits=2)), $(round(hi;digits=2)))", 16),
                    rpad(length(sel), 8), rpad(round(r; digits = 5), 14),
                    rpad(round(p0; digits = 5), 14),
                    rpad(round(r / p0; digits = 5), 10))
        end

        results[name] = Dict(
            "best_penalty"   => best_p,
            "pooled_ridge"   => overall_ridge,
            "pooled_prior"   => overall_prior,
            "pooled_ratio"   => overall_ridge / overall_prior,
            "bin_lo"         => bin_lo,
            "bin_hi"         => bin_hi,
            "bin_n"          => bin_n,
            "bin_ridge_rmse" => bin_ridge,
            "bin_prior_rmse" => bin_prior,
            "bin_ratio"      => bin_ridge ./ bin_prior,
        )
    end

    # --- A positive control -------------------------------------------------------------------
    # rho_true IS strongly encoded in the patch-correlation summary (that is the whole premise of
    # the method), so the same ridge on the same predictors MUST recover it well. If this control
    # failed, the null on the shift would be uninformative -- it would just mean the harness is
    # broken. This is what separates "the summary lacks shift information" from "this script is
    # wrong".
    yfit = vec(TH[1, fit_i]); yval = vec(TH[1, val_i]); ytst = vec(TH[1, test_i])
    best_p, best_v = RIDGE_GRID[1], Inf
    for p in RIDGE_GRID
        w, b = ridge_fit(Xfit, yfit, p)
        v = _rmse(ridge_predict(Xval, w, b), yval)
        v < best_v && ((best_p, best_v) = (p, v))
    end
    w, b = ridge_fit(Xfit, yfit, best_p)
    ctrl_ridge = _rmse(ridge_predict(Xtst, w, b), ytst)
    ctrl_prior = sqrt(mean(abs2, ytst .- mean(yfit)))
    println()
    println("-"^78)
    println("POSITIVE CONTROL  rho_true (must be well recovered, else the harness is broken)")
    println("-"^78)
    println("  ridge RMSE = $(round(ctrl_ridge; digits = 5))   " *
            "prior RMSE = $(round(ctrl_prior; digits = 5))   " *
            "ratio = $(round(ctrl_ridge / ctrl_prior; digits = 5))")

    elapsed = (time() - t_start) / 60
    tmp = RECOVERY_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version   = 1,
        n_total          = N,
        n_fit            = length(fit_i),
        n_val            = length(val_i),
        n_test           = length(test_i),
        predictor_rows   = 128,
        row_129_excluded = true,
        ridge_grid       = RIDGE_GRID,
        shift_dx         = results["shift_dx"],
        shift_dy         = results["shift_dy"],
        control_rho_ridge_rmse = ctrl_ridge,
        control_rho_prior_rmse = ctrl_prior,
        control_rho_ratio      = ctrl_ridge / ctrl_prior,
        pool_dir         = dir,
        elapsed_min      = elapsed,
        generated        = string(Dates.now(Dates.UTC)) * "Z",
        caption          = "Phase-11 shift-recoverability: ridge on the 128-row summary with the " *
                           "lambda row EXCLUDED, against the conditional-prior baseline, per " *
                           "lambda rung, with a rho_true positive control. Diagnostic; gates nothing.",
    )
    let d = JLD2.load(tmp); @assert haskey(d, "shift_dx") && haskey(d, "control_rho_ratio"); end
    mv(tmp, RECOVERY_REPORT_PATH; force = true)
    println()
    println("wrote ", RECOVERY_REPORT_PATH, "  (", round(elapsed; digits = 2), " min)")
    return nothing
end

isdefined(@__MODULE__, :P11_RECOVERY_LOAD_ONLY) || main()
