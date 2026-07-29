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

# spike/validation/run_p12_eps_ridge.jl --- Phase-12 chromatic-eps identifiability probe (R-6).
#
# THE QUESTION: is `chromatic_eps` recoverable from the 128-row patch-correlation summary AT ALL?
# It matters because epsilon is the only nuisance in the model that CLEARS the per-draw noise
# floor: at the prior edge it moves the summary by ||ds||_2 = 1.2417 against a floor of 0.9765
# (ratio 1.272), roughly 3x more than the maximal registration shift. If it is identified, the net
# can explain the radial gradient away and S-4's chromatic confound is materially mitigated -- and
# it would be this project's FIRST identified nuisance, which changes how named limit #4 (vacuous
# nuisance columns) has to be written up. If it is vacuous, the confound is LIVE and the 12-20
# guards carry the weight. Either answer is publishable; not measuring it is not.
#
# THE PREDICTORS: ALL 128 SUMMARY ROWS. UNLIKE THE PHASE-11 SHIFT PROBE THERE IS NO CONDITIONING
# ROW TO EXCLUDE. `run_p11_recovery.jl` had to drop row 129 (`encode_lambda(lambda)`) because the
# Phase-11 training joint draws `shift ~ Uniform(-lambda, lambda)`, so a predictor that sees lambda
# reproduces the conditional prior spread as PRIOR ECHO rather than as learning. `chromatic_eps` is
# never conditioned on -- ruling Q3 of the Phase-11 datagen header is explicit that it is a FIXED
# marginalized nuisance drawn from CHROMATIC_PRIOR, deliberately NOT lambda-scaled and carrying no
# conditioning input of its own. There is therefore no echo channel to close here, and the 128 raw
# summary rows are the whole predictor set. A SECOND ARM uses only the 64 CONTINUOUS rows: if only
# the full arm works, the present-mask rows are carrying the signal, which is itself informative.
#
# THE BASELINE: predicting the prior mean of `CHROMATIC_PRIOR` (= Uniform(-0.02, 0.02), so the
# analytic mean is 0 by symmetry). It is computed EMPIRICALLY from the REALIZED epsilon in the FIT
# split rather than taken as the analytic 0, so the comparison cannot be flattered by a mismatch
# between the nominal and the realized draw distribution. Ridge beating this baseline means the
# summary carries recoverable chromatic information; ridge merely matching it means it does not.
#
# Ridge is the right tool here precisely because it is weak and transparent: closed-form, no tuning
# theatre, no capacity to memorise. If a 128-predictor linear map cannot beat the prior, the honest
# reading is that the LINEAR signal is absent.
#
# THE POSITIVE CONTROL: the same ridge, the same split, the same predictors, targeting `rho_true`
# (theta row 1). Phase 11 measured that control at ratio 0.157 -- an 84 % error reduction
# (11-DIAGNOSIS.md section 2). IF THAT CONTROL DOES NOT COME IN COMFORTABLY BELOW 1.0, THE EPSILON
# NULL IS UNINFORMATIVE and means only that this harness is broken, not that epsilon is vacuous.
#
# THE SHRINKAGE DIAGNOSTIC IS DELIBERATELY NOT NAMED `shrinkage` (F3). This file reports
# `ridge_residual_shrinkage` = residual_sd / prior_sd: the out-of-sample RMSE of a LINEAR POINT
# PREDICTOR over the realized prior spread. 12-18 reports a key called `shrinkage` = post_sd /
# prior_sd: the WIDTH of a neural flow's posterior. Same idea, different machinery; they coincide
# only under a well-specified linear-Gaussian model. One key name for two estimators is how a
# reader ends up comparing quantities that were never comparable, so they are reported apart and
# read as a pair, NEVER averaged. Both are compared against the one shared floor
# `P12_VACUOUS_SHRINKAGE_FLOOR`, which is a member of `P12_REPORTING_ONLY_CONSTANTS`.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; reads the frozen Phase-12 pre-registration
# and the EXISTING Phase-11 training pool; writes ONE new artifact. Trains no network, mutates no
# constant, appends no Tier-2 constant, touches no `src/` file, and -- the point of the assertion
# below -- DOES NOT REGENERATE OR OVERWRITE THE PHASE-11 POOL. The pool is opened through
# `p11_pool_dir` + `load_p11_pool` only; the pool GENERATOR and the shard WRITER are never called
# from this file, which is checked on the comment-stripped source rather than trusted from prose.
#
# INCLUDE ORDER, AND WHY IT IS NOT THE OBVIOUS ONE. `p11_consts.jl` guards its whole Tier-1 block
# on `isdefined(@__MODULE__, :P11_DEV_SEED)`, and `p12_consts.jl:108` legitimately BINDS
# `P11_DEV_SEED` (it forbids the Phase-11 stream by name). Loading the Phase-12 constants FIRST
# therefore makes the Phase-11 Tier-1 block a silent no-op -- `LAMBDA_MIN` never gets defined, and
# `p11_consts.jl`'s Tier-2 block then throws on `SC2_SPEARMAN_ATTENUATION`. Measured, not guessed.
# So the Phase-11 chain is loaded FIRST and the Phase-12 pre-registration SECOND. Both orders were
# run; only this one loads. The re-binding of the shared names is benign: they carry identical
# literals in both files, and `p12_consts.jl:424-428` already records the one deliberate exception
# (`Z_TWO_SIDED_90`, more digits of the same quantile) as legal under Julia 1.12.
#
# SEEDING. The split is a DETERMINISTIC tail block with no shuffle, exactly as the analog, so this
# run consumes NO random draw at all. `P12_EPS_RIDGE_COUNTER` is recorded in the artifact as the
# audit-trail tag for this activity rather than consumed as a key word -- the same distinction
# `spike/data/p11_generate.jl:92-96` draws for its own datagen counter. Had a shuffle been needed,
# `p12_rng(P12_EPS_RIDGE_COUNTER)` is the reserved stream for it.
#
# SCOPE: THIS IS A DIAGNOSTIC, NOT A DELIVERABLE, AND IT GATES NOTHING (R-6, reported-only). It
# changes no threshold, authorises no descope, and appends nothing to the pre-registration.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p12_eps_ridge.jl

using JLD2
using Dates
using Statistics
using LinearAlgebra

# ORDER IS LOAD-BEARING -- see the INCLUDE ORDER note in the header. The Phase-11 chain first
# (it pulls p11_consts.jl, the frozen prior and the sharded-cache layer), the Phase-12
# pre-registration second. The first guard is keyed on `load_p11_pool` -- the READER this file
# actually calls -- rather than on the pool generator the analog keys on, so that the name never
# has to appear outside a comment.
isdefined(@__MODULE__, :load_p11_pool) || include(joinpath(@__DIR__, "..", "data", "p11_generate.jl"))
isdefined(@__MODULE__, :P12_DEV_SEED)  || include(joinpath(@__DIR__, "p12_consts.jl"))

const EPS_RIDGE_REPORT_PATH = joinpath(@__DIR__, "p12_eps_ridge_report.jld2")

# Held-out fraction. A plain index split is sufficient and is the leak-free discipline that
# matters here: the transform and the ridge coefficients are fitted on TRAIN COLUMNS ONLY and
# applied to test, so no test-set statistic can leak into the fit.
const TEST_FRACTION = 0.25

# Ridge penalties swept on a validation slice carved out of TRAIN (never out of TEST).
const RIDGE_GRID = [1e-6, 1e-4, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0]

# The two reported predictor arms. `encode_d01` lays the summary out as rows 1:64 = patch values
# and rows 65:128 = the present-mask (src/amortized/summary.jl:72-76, column-major), so the
# continuous arm is exactly the first 64 rows.
const EPS_ARMS = (("all128", 1:128), ("continuous64", 1:64))

# Five |rho| bins. The epsilon perturbation grows with |rho| BY PHYSICS: a radial magnification
# mismatch degrades an EXISTING correlation, so at rho ~ 0 there is nothing for it to degrade and
# a pooled number could hide a strongly rho-dependent effect.
const RHO_BIN_EDGES = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

# theta row indices, read from the frozen prior field order rather than hard-coded.
const THETA_FIELDS = (:ρ_true, :spillover, :autofluorescence, :label_efficiency,
                      :shift_dx, :shift_dy, :noise, :chromatic_eps)

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

"""
    _select_and_fit(Xfit, yfit, Xval, yval, Xtst) -> (pred, best_p)

Sweep `RIDGE_GRID`, select the penalty on the VALIDATION block, refit on FIT, predict on TEST.
The test block is touched exactly once, at the end, by construction.
"""
function _select_and_fit(Xfit, yfit, Xval, yval, Xtst)
    best_p, best_v = RIDGE_GRID[1], Inf
    for p in RIDGE_GRID
        w, b = ridge_fit(Xfit, yfit, p)
        v = _rmse(ridge_predict(Xval, w, b), yval)
        v < best_v && ((best_p, best_v) = (p, v))
    end
    w, b = ridge_fit(Xfit, yfit, best_p)
    return ridge_predict(Xtst, w, b), best_p
end

function main()
    println("="^78)
    println("PHASE-12 CHROMATIC-EPS IDENTIFIABILITY PROBE  (theta row 8, R-6)")
    println("Ridge on the EXISTING Phase-11 pool vs the empirical prior-mean baseline.")
    println("READ-ONLY on the pool. Gates nothing.")
    println("="^78)
    t_start = time()

    # --- The pool: RESOLVE and READ. Never write. -------------------------------------------
    dir = p11_pool_dir(P11_N_PAIRS)
    @assert isdir(dir) "run_p12_eps_ridge: the resolved Phase-11 pool path is not a directory: $dir"
    # THE `isdir` LINE ABOVE CANNOT ACTUALLY FIRE, AND THAT IS WHY THIS ONE EXISTS.
    # `p11_pool_dir` -> `open_or_invalidate` (spike/data/cache.jl:162) ends in
    # `isdir(dir) || mkpath(dir)`, so the directory is CREATED when absent and the existence
    # assertion is satisfied by the resolver itself. The load-bearing check is therefore on the
    # SHARDS: an absent pool shows up as a freshly-created EMPTY hash directory, not as a missing
    # one. This check is what refuses to proceed in that case.
    shards = filter(f -> startswith(basename(f), "shard_") && endswith(f, ".jld2"),
                    readdir(dir; join = true))
    @assert !isempty(shards) """
        run_p12_eps_ridge: no shard files under the resolved Phase-11 pool directory
            $dir
        THIS SCRIPT WILL NOT CREATE THEM. Rebuilding that pool is a ~1 h job owned by the script
        `spike/data/p11_generate.jl`, and running it from here would silently substitute a fresh
        pool for the one 12-15 and 12-17 depend on. Report the missing pool as a BLOCKER and stop."""
    pool = load_p11_pool(dir)

    S  = Float64.(pool[:summary_min])       # 128 x N -- the RAW summary
    TH = Float64.(pool[:theta])             # 8 x N   -- prior field order
    N  = size(S, 2)

    @assert size(S, 1) == 128 "predictors must be exactly the 128 summary rows, got $(size(S, 1))"
    @assert size(TH, 1) == 8 """
        run_p12_eps_ridge: theta has $(size(TH, 1)) rows, expected 8. A PRE-`ca02b0e` POOL HAS NO
        chromatic_eps COLUMN AT ALL and cannot answer this question -- row 8 would be out of
        bounds, or (worse, on a 7-row read) the wrong column would be silently regressed. Commit
        `ca02b0e` is the D-09/D-11 edit that appended `chromatic_eps` as the 8th theta row. Do not
        regenerate the pool from here; report the schema mismatch as a BLOCKER."""

    idx_eps = findfirst(==(:chromatic_eps), THETA_FIELDS)
    idx_rho = findfirst(==(:ρ_true), THETA_FIELDS)
    @assert idx_eps == 8 && idx_rho == 1
    println("theta rows: chromatic_eps = $idx_eps, rho_true = $idx_rho   (N = $N)")
    println("pool dir (read-only) = $dir")

    # --- Leak-free split ---------------------------------------------------------------------
    ntest  = round(Int, TEST_FRACTION * N)
    test_i = (N - ntest + 1):N
    tr_all = 1:(N - ntest)
    nval   = round(Int, 0.2 * length(tr_all))
    val_i  = (length(tr_all) - nval + 1):length(tr_all)
    fit_i  = 1:(length(tr_all) - nval)
    println("split: fit = $(length(fit_i))  val = $(length(val_i))  test = $(length(test_i))")

    y_eps_fit = vec(TH[idx_eps, fit_i])
    y_eps_val = vec(TH[idx_eps, val_i])
    y_eps_tst = vec(TH[idx_eps, test_i])

    # THE BASELINE, computed EMPIRICALLY from the realized epsilon in the FIT split. The analytic
    # prior mean of CHROMATIC_PRIOR is 0 by symmetry; using the realized fit-split mean instead
    # means the comparison cannot be flattered by a nominal-vs-realized mismatch.
    eps_bar_fit = mean(y_eps_fit)
    # THE DENOMINATOR OF THE SHRINKAGE DIAGNOSTIC: the REALIZED prior spread, again from the fit
    # split only. The analytic value is recorded alongside it purely so the two are comparable.
    eps_prior_sd_fit      = std(y_eps_fit)
    eps_prior_sd_analytic = std(CHROMATIC_PRIOR)
    println("epsilon prior: empirical fit-split mean = $(round(eps_bar_fit; sigdigits = 4))  " *
            "sd = $(round(eps_prior_sd_fit; sigdigits = 6))   " *
            "(analytic sd = $(round(eps_prior_sd_analytic; sigdigits = 6)))")

    rho_tst = vec(TH[idx_rho, test_i])

    eps_ratio   = Dict{String,Float64}()
    ctrl_ratio  = Dict{String,Float64}()
    shrink_arm  = Dict{String,Float64}()
    vacuous_arm = Dict{String,Bool}()
    arm_details = Dict{String,Any}()

    for (arm_name, rows) in EPS_ARMS
        Srows = view(S, rows, :)
        mu, sd = _fit_standardizer(view(Srows, :, fit_i))          # TRAIN-ONLY fit
        Xfit = _apply_standardizer(view(Srows, :, fit_i), mu, sd)
        Xval = _apply_standardizer(view(Srows, :, val_i), mu, sd)
        Xtst = _apply_standardizer(view(Srows, :, test_i), mu, sd)

        pred, best_p = _select_and_fit(Xfit, y_eps_fit, Xval, y_eps_val, Xtst)

        ridge_rmse = _rmse(pred, y_eps_tst)
        prior_rmse = sqrt(mean(abs2, y_eps_tst .- eps_bar_fit))
        ratio      = ridge_rmse / prior_rmse
        shrink     = ridge_rmse / eps_prior_sd_fit

        eps_ratio[arm_name]   = ratio
        shrink_arm[arm_name]  = shrink
        vacuous_arm[arm_name] = shrink > P12_VACUOUS_SHRINKAGE_FLOOR

        println()
        println("-"^78)
        println("ARM $arm_name  ($(length(rows)) predictor rows)   " *
                "best ridge penalty = $best_p  (chosen on validation)")
        println("-"^78)
        println("  chromatic_eps   ridge RMSE = $(round(ridge_rmse; sigdigits = 6))   " *
                "prior RMSE = $(round(prior_rmse; sigdigits = 6))   " *
                "ratio = $(round(ratio; digits = 5))")
        println("  ridge_residual_shrinkage = $(round(shrink; digits = 5))   " *
                "(floor $(P12_VACUOUS_SHRINKAGE_FLOOR))  vacuous = $(vacuous_arm[arm_name])")
        println()
        println(rpad("|rho| bin", 16), rpad("n", 8), rpad("ridge RMSE", 14),
                rpad("prior RMSE", 14), rpad("ratio", 10))

        bin_lo, bin_hi, bin_n = Float64[], Float64[], Int[]
        bin_ridge, bin_prior  = Float64[], Float64[]
        for k in 1:(length(RHO_BIN_EDGES) - 1)
            lo, hi = RHO_BIN_EDGES[k], RHO_BIN_EDGES[k + 1]
            last_bin = k == length(RHO_BIN_EDGES) - 1
            sel = findall(x -> (abs(x) >= lo) & (abs(x) < hi + (last_bin ? 1e-9 : 0.0)), rho_tst)
            isempty(sel) && continue
            r  = _rmse(pred[sel], y_eps_tst[sel])
            p0 = sqrt(mean(abs2, y_eps_tst[sel] .- eps_bar_fit))
            push!(bin_lo, lo); push!(bin_hi, hi); push!(bin_n, length(sel))
            push!(bin_ridge, r); push!(bin_prior, p0)
            println(rpad("[$(round(lo; digits = 2)), $(round(hi; digits = 2)))", 16),
                    rpad(length(sel), 8), rpad(round(r; sigdigits = 6), 14),
                    rpad(round(p0; sigdigits = 6), 14),
                    rpad(round(r / p0; digits = 5), 10))
        end

        # The SAME ridge, the SAME split, the SAME predictors -- targeting rho_true. This is what
        # separates "the summary lacks chromatic information" from "this script is wrong".
        yc_fit = vec(TH[idx_rho, fit_i])
        yc_val = vec(TH[idx_rho, val_i])
        yc_tst = vec(TH[idx_rho, test_i])
        cpred, cbest_p = _select_and_fit(Xfit, yc_fit, Xval, yc_val, Xtst)
        cridge = _rmse(cpred, yc_tst)
        cprior = sqrt(mean(abs2, yc_tst .- mean(yc_fit)))
        ctrl_ratio[arm_name] = cridge / cprior
        println()
        println("  POSITIVE CONTROL rho_true   ridge RMSE = $(round(cridge; sigdigits = 6))   " *
                "prior RMSE = $(round(cprior; sigdigits = 6))   " *
                "ratio = $(round(cridge / cprior; digits = 5))   " *
                "(Phase 11 measured 0.157; penalty = $cbest_p)")

        arm_details[arm_name] = Dict(
            "predictor_rows"           => collect(rows),
            "n_predictor_rows"         => length(rows),
            "best_penalty"             => best_p,
            "eps_ridge_rmse"           => ridge_rmse,
            "eps_prior_rmse"           => prior_rmse,
            "eps_ratio"                => ratio,
            "ridge_residual_shrinkage" => shrink,
            "vacuous"                  => vacuous_arm[arm_name],
            "bin_lo"                   => bin_lo,
            "bin_hi"                   => bin_hi,
            "bin_n"                    => bin_n,
            "bin_ridge_rmse"           => bin_ridge,
            "bin_prior_rmse"           => bin_prior,
            "bin_ratio"                => bin_ridge ./ bin_prior,
            "control_rho_best_penalty" => cbest_p,
            "control_rho_ridge_rmse"   => cridge,
            "control_rho_prior_rmse"   => cprior,
            "control_rho_ratio"        => cridge / cprior,
        )
    end

    # THE PRIMARY ARM IS `all128` -- the full predictor set, the one the header's question is
    # asked about. The scalar top-level keys are that arm; the per-arm dictionaries carry both.
    primary = "all128"

    elapsed = (time() - t_start) / 60
    tmp = EPS_RIDGE_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version            = 1,
        n_total                   = N,
        n_fit                     = length(fit_i),
        n_val                     = length(val_i),
        n_test                    = length(test_i),
        ridge_grid                = RIDGE_GRID,
        test_fraction             = TEST_FRACTION,
        zscore_arm                = "per_row",
        theta_row_eps             = idx_eps,
        theta_row_rho             = idx_rho,
        eps_ratio                 = eps_ratio,
        eps_ratio_primary         = eps_ratio[primary],
        primary_arm               = primary,
        control_rho_ratio         = ctrl_ratio[primary],
        control_rho_ratio_by_arm  = ctrl_ratio,
        ridge_residual_shrinkage  = shrink_arm[primary],
        ridge_residual_shrinkage_by_arm = shrink_arm,
        vacuous                   = vacuous_arm[primary],
        vacuous_by_arm            = vacuous_arm,
        vacuous_shrinkage_floor   = P12_VACUOUS_SHRINKAGE_FLOOR,
        eps_prior_mean_fit        = eps_bar_fit,
        eps_prior_sd_fit          = eps_prior_sd_fit,
        eps_prior_sd_analytic     = eps_prior_sd_analytic,
        rho_bin_edges             = RHO_BIN_EDGES,
        arms                      = arm_details,
        pool_dir                  = dir,
        p11_pool_readonly         = true,
        counter                   = P12_EPS_RIDGE_COUNTER,
        elapsed_min               = elapsed,
        generated                 = string(Dates.now(Dates.UTC)) * "Z",
        caption                   = "Phase-12 chromatic-eps identifiability: ridge on the 128-row " *
                                    "summary of the EXISTING Phase-11 50k pool (read-only), " *
                                    "targeting theta row 8, against the empirically-computed " *
                                    "prior-mean baseline, in two predictor arms, broken down by " *
                                    "|rho|, with a rho_true positive control. Reported only; it " *
                                    "changes no threshold and appends no constant. " *
                                    "Diagnostic; gates nothing.",
    )
    let d = JLD2.load(tmp)
        @assert haskey(d, "eps_ratio")                "artifact integrity: eps_ratio missing"
        @assert haskey(d, "control_rho_ratio")        "artifact integrity: control_rho_ratio missing"
        @assert haskey(d, "ridge_residual_shrinkage") "artifact integrity: ridge_residual_shrinkage missing"
        @assert haskey(d, "vacuous")                  "artifact integrity: vacuous missing"
    end
    mv(tmp, EPS_RIDGE_REPORT_PATH; force = true)

    println()
    println("="^78)
    println("REPORTED, NOT GATED. eps_ratio = ", eps_ratio,
            "   control_rho_ratio = ", round(ctrl_ratio[primary]; digits = 5))
    println("wrote ", EPS_RIDGE_REPORT_PATH, "  (", round(elapsed; digits = 2), " min)")
    println("="^78)
    return nothing
end

isdefined(@__MODULE__, :P12_EPS_RIDGE_LOAD_ONLY) || main()
