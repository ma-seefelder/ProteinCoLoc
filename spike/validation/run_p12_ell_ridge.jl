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

# spike/validation/run_p12_ell_ridge.jl --- the D-08 CORRELATION-LENGTH IDENTIFIABILITY PROBE.
#
# THE QUESTION: is the spatial correlation length r1 recoverable from the 128-row patch-correlation
# summary AT ALL? D-08 hedged -- "expect this column may prove unidentifiable; report it as such if
# so" -- and S-1 is explicit that prose is not enough here. This turns the hedge into a measurement.
#
# NO ROW IS EXCLUDED, AND THAT IS THE DIFFERENCE FROM THE PHASE-11 PROBE. `run_p11_recovery.jl`
# excluded row 129 because lambda was a CONDITIONING INPUT: admitting it would have measured prior
# echo rather than recovery, and "row 129 is excluded and this is the whole point" is that file's
# most quotable line. D-08 INFERS r1 rather than conditioning on it, so there is no conditioning row
# here, nothing to exclude, and no prior-echo channel of the Phase-11 kind. Copying the analog's
# exclusion sentence into this header would be false, so it is not copied.
#
# THE TWO ARMS, AND WHY ONE WOULD NOT DO.
#   (a) `raw`        -- the 128 raw summary rows.
#   (b) `with_moran` -- the same 128 rows PLUS Moran's I of the 64 continuous rows, computed under
#                       the lattice's OWN 4-neighbour adjacency.
# r1 is a property of the SPATIAL AUTOCORRELATION of the 64 continuous rows, not of their levels. A
# linear map on raw levels can barely see a second-order quantity, so arm (a) alone is a WEAK PROBE:
# reporting it by itself would let an implementation limitation read as a physical null -- the same
# error class as Pitfall 3, run in the opposite direction. Arm (b) is the honest one and is the
# PRIMARY; arm (a) is reported beside it so a reader can see how much of the answer is second-order.
# Neither is derived from the other: both are fitted from scratch, on their own design.
#
# THE BASELINE: predict the prior mean of r1, computed EMPIRICALLY from the realized draws in the
# FIT split -- never from the analytic `P12_R1_PRIOR` mean. The comparison must not be flattered by
# any mismatch between the nominal and the realized draw distribution.
#
# THE POSITIVE CONTROL: the SAME ridge, the SAME split, the SAME predictors, targeting the global
# DCT coefficient c0 -- the global field level, whose Phase-11 analogue was recovered at ratio 0.157
# (12-CONTEXT S-2). WITHOUT IT, "r1 is unidentifiable" and "this script is wrong" are the same
# observation. If `control_c0_ratio` is not comfortably below 1.0, BOTH arms are uninformative and
# no identifiability verdict may be recorded from this run.
#
# THE POOL DRAWS r1; IT DOES NOT PIN IT. The target must VARY or there is nothing to regress. That
# also puts this pool in a DIFFERENT content-hash directory from 12-11's five pinned rung pools
# (`r1_pin` reaches the hash), which is asserted below rather than assumed -- a collision would
# silently serve one rung's constant-r1 data to a probe whose whole premise is that r1 varies.
#
# `imsize_tag` IS PASSED EXPLICITLY. The image-size sampler is a closure and cannot be hashed, so a
# pinned-size pool WITHOUT a tag hashes into the F5-mixture directory and is indistinguishable from
# one (12-09 section 6.1). The tag is DERIVED from `P12_STAGE1_IMSIZE`, never typed.
#
# THE SHRINKAGE DIAGNOSTIC, AND WHY IT IS NOT CALLED `shrinkage`.
# `ridge_residual_shrinkage` = residual_sd / prior_sd: the out-of-sample RMSE of a LINEAR POINT
# PREDICTOR over the realized prior sd. 12-18 stores a key called `shrinkage`: the WIDTH of a NEURAL
# FLOW's posterior over the same prior sd. Both estimate the same IDEA -- how much of the prior
# spread survives once the image is seen -- by DIFFERENT MACHINERY, and they coincide only under a
# well-specified linear-Gaussian model, which nothing here guarantees. Giving both the name
# `shrinkage` would ASSERT a commensurability that has not been ARGUED. They are named apart on
# purpose. Read them as a PAIR under the rule fixed in 12-13:131-152 BEFORE either number existed:
#   both near 1.0            -> CORROBORATION; the image contributed nothing, D-08 resolved negative
#   ridge ~1, posterior < 1  -> do NOT claim identifiability on the posterior alone (Phase 11
#                               measured coverage 0.723-0.767 against a nominal 0.90 -- intervals
#                               ~40 % too narrow -- and a linear probe would not have shown it)
#   ridge < 1, posterior ~1  -> the information IS there and the estimator is not using it: a
#                               CAPACITY limitation, so D-08's hedge is UNRESOLVED, not negative
#   they disagree            -> that is a FINDING, reported as one. Never averaged, never reconciled.
#
# NO PRIOR-ECHO CEILING IS RECORDED HERE, AND ITS ABSENCE IS DELIBERATE. An earlier draft of the
# plan recorded `prior_echo_ceiling_for_reference = P12_R1_MAX / P12_R1_MIN = 19.0`. It is WITHDRAWN
# rather than quietly dropped, because the reasoning is this milestone's most expensive recurring
# error: 19.0 is a ratio of PRIOR SUPPORT BOUNDS presented as a ratio of POSTERIOR WIDTHS. Phase
# 11's 12.0 was meaningful only because lambda was a CONDITIONING INPUT -- it bounded how much a
# conditional marginal could widen ACROSS lambda levels, and there were levels to compare. r1 is
# INFERRED, so there is no across-level width ratio to bound: a network ignoring the image returns
# the same prior for every dataset, which is a width ratio of 1, not 19. It is also not
# commensurable with 12-18's `post_sd / prior_sd`, so printing them together would put two different
# quantities under one caption -- in a reported artifact key.
#
# SCOPE: THIS IS A DIAGNOSTIC, NOT A DELIVERABLE, AND IT **gates nothing** (S-1, D-08). It appends
# no Tier-2 constant, moves no threshold and authorises no descope. It runs on BOTH the PROCEED and
# the DESCOPE route, because the identifiability answer shapes how named limit #4 is written either
# way.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; generates its own pool under
# `spike/data/cache/p12/` through the one datagen path. Touches no `src/` file, and does not read,
# write or regenerate Phase 11's pool.
#
# SEEDING. The split is a DETERMINISTIC tail block with no shuffle, so the ridge half of this run
# consumes NO random draw at all. `P12_ELL_RIDGE_COUNTER` is recorded in the artifact as the
# audit-trail tag for this activity rather than consumed as a key word -- the same distinction
# `spike/data/p12_generate.jl:129-132` draws for its own datagen counter. The pool rides the datagen
# stream `P12_DEV_SEED xor P12_DATAGEN_SALT`, keyed per global index.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p12_ell_ridge.jl

using JLD2
using Dates
using Statistics
using LinearAlgebra

# ORDER: the Tier-1 pre-registration first, then the datagen path (which transitively pulls the
# contract, the Phase-12 prior, the forward model, the encoder, the sharded cache AND
# `p12_theta_index` via p12_architecture.jl), then the lattice arithmetic for `lattice_adjacency`.
# Only Phase-12-local sentinel names are guarded on; no foreign pre-registration is loaded into this
# process, so the include-guard poisoning class (p12_consts.jl BINDS :P11_DEV_SEED and
# :P13_DEV_SEED) cannot arise here.
isdefined(@__MODULE__, :P12_DEV_SEED)       || include(joinpath(@__DIR__, "p12_consts.jl"))
isdefined(@__MODULE__, :generate_p12_pool)  || include(joinpath(@__DIR__, "..", "data", "p12_generate.jl"))
isdefined(@__MODULE__, :lattice_adjacency)  || include(joinpath(@__DIR__, "..", "simulator", "p12_lattice.jl"))

const ELL_RIDGE_REPORT_PATH = joinpath(@__DIR__, "p12_ell_ridge_report.jld2")

# Ridge penalties swept on a validation slice carved out of TRAIN (never out of TEST). These are
# REGULARIZATION STRENGTHS, not thresholds: nothing here decides a pass or a fail. Taken unchanged
# from `run_p11_recovery.jl:129-158` and `run_p12_stage1_ridge.jl:124` so the probes are comparable.
const ELL_RIDGE_GRID = [1e-6, 1e-4, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0]

# CAR is the phase's default arm. The CAR-vs-GP comparison is 12-15's question, not this probe's.
const ELL_ARM = :car

# The pinned-size tag that keeps this pool out of the F5-mixture directory. DERIVED from the frozen
# size, never typed, so the tag cannot name a size the pool does not have.
const ELL_IMSIZE_TAG = Symbol("pinned_", P12_STAGE1_IMSIZE[1], "x", P12_STAGE1_IMSIZE[2])

# Split fractions, not thresholds; taken unchanged from the analogs so the probes are comparable.
const ELL_TEST_FRACTION = P12_STAGE1_TEST_FRACTION
const ELL_VAL_FRACTION  = 0.2

# FIVE bins over the realized r1, so the expected non-uniformity is visible instead of averaged.
# Six edges make five bins. Pitfall 2's mechanism predicts the shape: short lengths give large
# within-image contrast, long lengths give a nearly constant field.
const ELL_BIN_PROBS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)

"""
    morans_i(rows, present; W = lattice_adjacency(P12_G)) -> Float64

The standard Moran's I of one sample's 64 RAW continuous summary rows arranged on the lattice:

    I = (n / sum(W)) * (z' * W * z) / (z' * z),   z = rows .- mean(rows)

`W` DEFAULTS TO `lattice_adjacency(P12_G)` AND IS NOT RE-DERIVED HERE. A probe that used a different
neighbour relation than the one defining the prior would be measuring a different question -- so the
adjacency comes from the same function `car_sigma`/`induced_lag1` use, and this file contains no
second neighbour loop. Callers pass the precomputed `W` to avoid rebuilding it per sample.

MASKED REGIONS CONTRIBUTE NO NEIGHBOUR TERM, STATED RATHER THAN LEFT TO A ZERO. `present` is the
sample's mask row as a Bool vector. A masked patch's continuous entry is not a measurement of
anything, and letting its stored zero flow into `z` would make an ABSENT region look like a region
sitting exactly at the mean -- a phantom observation that would bias I toward whatever the zero
happens to imply. Both the mean and every edge term are therefore restricted to present regions,
and an edge counts only when BOTH of its endpoints are present.

Returns `NaN` when fewer than two regions are present or the present values are constant: there is
no autocorrelation to measure, and a NaN is honest where a 0.0 would be a measurement.
"""
function morans_i(rows::AbstractVector, present::AbstractVector{Bool};
                  W = lattice_adjacency(P12_G))
    idx = findall(present)
    length(idx) < 2 && return NaN
    x = collect(float.(rows[idx]))
    z = x .- mean(x)
    denom = sum(abs2, z)
    denom > 0 || return NaN
    Wsub  = W[idx, idx]                  # an edge survives only if BOTH endpoints are present
    wsum  = sum(Wsub)
    wsum > 0 || return NaN
    return (length(idx) / wsum) * (z' * (Wsub * z)) / denom
end

"""
    _p12ell_fit_standardizer(X) -> (mu, sd)

Column-wise (per predictor row) mean and sd over FIT columns only. Rows with zero variance get
`sd = 1` so they pass through as constants instead of producing NaNs -- which the near-constant
present-mask block needs. Ported unchanged from `run_p11_recovery.jl:75-82`; the `1e-12` is that
file's numerical zero-variance guard, not a bar this file invents.
"""
function _p12ell_fit_standardizer(X::AbstractMatrix)
    mu = vec(mean(X; dims = 2))
    sd = vec(std(X; dims = 2))
    @inbounds for i in eachindex(sd)
        (isfinite(sd[i]) && sd[i] > 1e-12) || (sd[i] = 1.0)
    end
    return mu, sd
end

_p12ell_apply_standardizer(X, mu, sd) = (X .- mu) ./ sd
_p12ell_predict(Xs, w, b) = (Xs' * w) .+ b
_p12ell_rmse(pred, truth) = sqrt(mean(abs2, pred .- truth))
_p12ell_gram(Xfit::AbstractMatrix) = Symmetric(Xfit * Xfit')

"""
    _p12ell_select_and_fit(A, Xfit, yfit, Xval, yval, Xtst) -> (pred, best_penalty)

Sweep `ELL_RIDGE_GRID`, select the penalty on the VALIDATION block, predict on TEST with the
selected fit. The test block is touched exactly ONCE, at the end, by construction -- there is no
path through this function on which a test-set statistic can reach the selection. Every solve
carries a penalty, so a rank-deficient design cannot throw `SingularException`.
"""
function _p12ell_select_and_fit(A, Xfit, yfit, Xval, yval, Xtst)
    b   = mean(yfit)                       # unpenalised intercept, from FIT only
    rhs = Xfit * (yfit .- b)
    best_p = ELL_RIDGE_GRID[1]
    best_w = (A + best_p * I) \ rhs
    best_v = _p12ell_rmse(_p12ell_predict(Xval, best_w, b), yval)
    for p in ELL_RIDGE_GRID
        w = (A + p * I) \ rhs
        v = _p12ell_rmse(_p12ell_predict(Xval, w, b), yval)   # VALIDATION, never TEST
        if v < best_v
            best_p, best_v, best_w = p, v, w
        end
    end
    return _p12ell_predict(Xtst, best_w, b), best_p
end

function main()
    println("="^78)
    println("D-08 CORRELATION-LENGTH IDENTIFIABILITY PROBE -- two arms, live c0 control")
    println("NO row is excluded: r1 is INFERRED, not conditioned. This runner gates nothing.")
    println("="^78)
    t_start = time()

    G    = P12_G
    nreg = G^2
    n_ask = P12_STAGE1_N

    println("n asked      = $n_ask   arm = $ELL_ARM   imsize = $(P12_STAGE1_IMSIZE)  tag = $ELL_IMSIZE_TAG")
    println("r1           = DRAWN from $(P12_R1_PRIOR)   (pinning it would leave nothing to regress)")

    # --- The pool must NOT land in one of 12-11's pinned rung directories --------------------
    # `p12_pool_dir` resolves through `open_or_invalidate`, which CREATES the directory it
    # resolves, so this is a HASH PROBE and not a presence test (12-09 section 5;
    # `p12_pool_complete` is the presence test). What it proves is the only thing at stake: that
    # a DRAWN-r1 pool hashes somewhere no PINNED-r1 pool does.
    probe_dir = p12_pool_dir(n_ask; arm = ELL_ARM, r1 = nothing, imsize_tag = ELL_IMSIZE_TAG)
    rung_dirs = [p12_pool_dir(P12_STAGE1_N ÷ length(P12_R1_LADDER); arm = ELL_ARM, r1 = rung,
                              imsize_tag = ELL_IMSIZE_TAG) for rung in P12_R1_LADDER]
    @assert probe_dir ∉ rung_dirs """
        run_p12_ell_ridge: the DRAWN-r1 pool resolved to $probe_dir, which is one of 12-11's
        PINNED rung directories. r1 is not reaching the content hash, so this probe would be
        served a pool in which r1 is CONSTANT -- and a regression on a constant target is not a
        null result, it is a meaningless one. Refusing to generate."""

    # And against what 12-11 actually WROTE, not only against what the hash recomputes.
    let s1 = joinpath(@__DIR__, "p12_stage1_report.jld2")
        if isfile(s1)
            recorded = JLD2.load(s1, "pool_dirs")
            @assert probe_dir ∉ recorded "run_p12_ell_ridge: the resolved pool directory " *
                "$probe_dir is recorded in p12_stage1_report.jld2's pool_dirs. See above."
            println("checked against the $(length(recorded)) pool dirs recorded by 12-11: distinct")
        else
            println("NOTE: p12_stage1_report.jld2 absent; checked against recomputed rung hashes only")
        end
    end

    dir = generate_p12_pool(n_ask; arm = ELL_ARM,
                            imsize_sampler = _ -> P12_STAGE1_IMSIZE,
                            r1 = nothing, imsize_tag = ELL_IMSIZE_TAG)
    @assert dir == probe_dir "run_p12_ell_ridge: generated pool landed in $dir, the hash probe " *
        "resolved $probe_dir"
    @assert p12_pool_complete(dir, n_ask) "run_p12_ell_ridge: pool at $dir is not complete for n = $n_ask"

    pool = load_p12_pool(dir)
    S  = Float64.(pool.summary_min)          # 128 x N, RAW (never z-scored at rest)
    TH = Float64.(pool.theta)                # D x N, P12_THETA_ROWS order
    N  = size(S, 2)

    # --- STRUCTURAL aborts. None of these is a result. ---------------------------------------
    @assert N == n_ask "run_p12_ell_ridge: pool at $dir holds $N samples, expected $n_ask"
    @assert size(S, 1) == P12_SUMMARY_MIN_DIM "run_p12_ell_ridge: predictors must be the " *
        "$(P12_SUMMARY_MIN_DIM) encoded summary rows, got $(size(S, 1))"
    @assert size(TH, 1) == P12_D_MINISPIKE "run_p12_ell_ridge: theta has $(size(TH, 1)) rows, " *
        "expected P12_D_MINISPIKE = $(P12_D_MINISPIKE)"
    @assert all(==(P12_STAGE1_IMSIZE), pool.imsize) "run_p12_ell_ridge: pool at $dir carries " *
        "image sizes $(unique(pool.imsize)), expected the pinned $(P12_STAGE1_IMSIZE)"
    # THE TARGET MUST VARY. A drawn-r1 pool that arrived constant would make every ratio below
    # meaningless, and it is exactly what a hash collision with a pinned rung would produce.
    @assert length(unique(pool.r1)) > 1 """
        run_p12_ell_ridge: the loaded pool carries a SINGLE r1 value $(first(pool.r1)). The draw
        was requested but a pinned pool was served. Refusing to report a ratio."""

    idx_r1 = p12_theta_index(:r1; K = P12_K_DEV)
    idx_c0 = p12_theta_index(:c0; K = P12_K_DEV)
    println("theta rows: r1 = $idx_r1, c0 = $idx_c0   (N = $N)")
    println("pool dir = $dir")
    # The theta row and the pool's own realized record must be the SAME quantity. A silent
    # off-by-one in the row layout would regress a nuisance and report it as r1.
    @assert maximum(abs.(vec(TH[idx_r1, :]) .- collect(float.(pool.r1)))) < 1e-12 """
        run_p12_ell_ridge: theta row $idx_r1 does not equal the pool's realized r1 vector. The
        theta layout and the pool record disagree; refusing to regress a row that may not be r1."""

    # --- Leak-free split: deterministic tail block, no shuffle, exactly as the analogs --------
    ntest  = round(Int, ELL_TEST_FRACTION * N)
    test_i = (N - ntest + 1):N
    tr_all = 1:(N - ntest)
    nval   = round(Int, ELL_VAL_FRACTION * length(tr_all))
    val_i  = (length(tr_all) - nval + 1):length(tr_all)
    fit_i  = 1:(length(tr_all) - nval)
    println("split: fit = $(length(fit_i))  val = $(length(val_i))  test = $(length(test_i))")

    y_fit = vec(TH[idx_r1, fit_i]); y_val = vec(TH[idx_r1, val_i]); y_tst = vec(TH[idx_r1, test_i])
    c_fit = vec(TH[idx_c0, fit_i]); c_val = vec(TH[idx_c0, val_i]); c_tst = vec(TH[idx_c0, test_i])

    # THE BASELINE and THE SHRINKAGE DENOMINATOR, both EMPIRICAL and both from the FIT split only.
    r1_bar_fit          = mean(y_fit)
    r1_prior_sd_fit     = std(y_fit)
    r1_prior_sd_analytic = std(P12_R1_PRIOR)
    println("r1 prior: empirical fit-split mean = $(round(r1_bar_fit; sigdigits = 4))  " *
            "sd = $(round(r1_prior_sd_fit; sigdigits = 6))   " *
            "(analytic sd = $(round(r1_prior_sd_analytic; sigdigits = 6)))")

    # --- Moran's I of every sample's 64 continuous rows, under the PRIOR'S OWN adjacency -------
    W = lattice_adjacency(G)
    @assert size(W) == (nreg, nreg)
    moran = Vector{Float64}(undef, N)
    @inbounds for j in 1:N
        present = view(S, (nreg + 1):(2 * nreg), j) .>= 0.5     # the mask block, 1 = observed
        moran[j] = morans_i(view(S, 1:nreg, j), present; W = W)
    end
    n_moran_nan = count(isnan, moran)
    # A NaN predictor would poison the design. Samples with too few present regions to define I are
    # given the FIT-SPLIT mean of the finite values -- an explicit, stated imputation rather than a
    # silent zero, and recorded in the artifact so its size is visible.
    finite_fit = filter(isfinite, moran[fit_i])
    moran_fill = isempty(finite_fit) ? 0.0 : mean(finite_fit)
    @inbounds for j in 1:N
        isnan(moran[j]) && (moran[j] = moran_fill)
    end
    moran_i_summary = (mean = mean(moran), sd = std(moran),
                       min = minimum(moran), max = maximum(moran),
                       n_nan_imputed = n_moran_nan, fill_value = moran_fill,
                       cor_with_r1 = cor(moran, vec(TH[idx_r1, :])))
    println("Moran's I: mean = $(round(moran_i_summary.mean; digits = 5))  " *
            "sd = $(round(moran_i_summary.sd; digits = 5))  " *
            "range [$(round(moran_i_summary.min; digits = 4)), $(round(moran_i_summary.max; digits = 4))]  " *
            "NaN-imputed = $n_moran_nan")
    println("  raw corr(Moran's I, r1) = $(round(moran_i_summary.cor_with_r1; digits = 5))   " *
            "(REPORTED; the ridge below is the actual probe)")

    # --- The two designs. Neither is derived from the other. ----------------------------------
    S_moran = vcat(S, reshape(moran, 1, N))            # 129 rows: the 128 raw rows PLUS Moran's I
    ARMS = ("raw" => S, "with_moran" => S_moran)
    primary_arm_name = "with_moran"                     # arm (b) is the honest one; see the header

    ratio_arm    = Dict{String,Float64}()
    shrink_arm   = Dict{String,Float64}()
    vacuous_arm  = Dict{String,Bool}()
    control_arm  = Dict{String,Float64}()
    penalty_arm  = Dict{String,Float64}()
    bin_arm      = Dict{String,Vector{Float64}}()
    pred_arm     = Dict{String,Vector{Float64}}()

    # Bin edges from the REALIZED r1 in the FIT split, through `realized_r1_quantiles` rather than
    # from the intended prior -- and from FIT rather than TEST, so no test-set structure shapes the
    # reporting grid. Five bins, six edges.
    # NOTE `q.values`, not `values(q)`: `realized_r1_quantiles` returns the NamedTuple
    # `(probs = ..., values = ...)`, so `values(q)` would be `Base.values` over BOTH fields.
    q = realized_r1_quantiles(vec(TH[idx_r1, fit_i]); probs = ELL_BIN_PROBS)
    bin_edges = collect(float.(q.values))
    @assert issorted(bin_edges) "run_p12_ell_ridge: realized r1 quantiles are not monotone: $bin_edges"
    @assert length(bin_edges) == length(ELL_BIN_PROBS) "run_p12_ell_ridge: expected " *
        "$(length(ELL_BIN_PROBS)) realized bin edges, got $(length(bin_edges))"
    println("realized r1 bin edges (fit split, from realized_r1_quantiles) = " *
            "$(round.(bin_edges; digits = 4))")

    for (arm_name, X) in ARMS
        mu, sd = _p12ell_fit_standardizer(view(X, :, fit_i))       # FIT-ONLY fit
        Xfit = _p12ell_apply_standardizer(view(X, :, fit_i), mu, sd)
        Xval = _p12ell_apply_standardizer(view(X, :, val_i), mu, sd)
        Xtst = _p12ell_apply_standardizer(view(X, :, test_i), mu, sd)
        A    = _p12ell_gram(Xfit)

        pred, best_p = _p12ell_select_and_fit(A, Xfit, y_fit, Xval, y_val, Xtst)
        ridge_rmse = _p12ell_rmse(pred, y_tst)
        prior_rmse = sqrt(mean(abs2, y_tst .- r1_bar_fit))

        ratio_arm[arm_name]   = ridge_rmse / prior_rmse
        shrink_arm[arm_name]  = ridge_rmse / r1_prior_sd_fit
        vacuous_arm[arm_name] = shrink_arm[arm_name] > P12_VACUOUS_SHRINKAGE_FLOOR
        penalty_arm[arm_name] = best_p
        pred_arm[arm_name]    = pred

        # THE POSITIVE CONTROL: same design, same split, same penalty grid, target = c0.
        cpred, _ = _p12ell_select_and_fit(A, Xfit, c_fit, Xval, c_val, Xtst)
        control_arm[arm_name] =
            _p12ell_rmse(cpred, c_tst) / sqrt(mean(abs2, c_tst .- mean(c_fit)))

        println()
        println("-"^78)
        println("ARM $arm_name  ($(size(X, 1)) predictor rows)   " *
                "best ridge penalty = $best_p  (chosen on VALIDATION)")
        println("-"^78)
        println("  r1   ridge RMSE = $(round(ridge_rmse; sigdigits = 6))   " *
                "prior RMSE = $(round(prior_rmse; sigdigits = 6))   " *
                "ratio = $(round(ratio_arm[arm_name]; digits = 5))")
        println("  ridge_residual_shrinkage = $(round(shrink_arm[arm_name]; digits = 5))   " *
                "(floor $(P12_VACUOUS_SHRINKAGE_FLOOR))  vacuous = $(vacuous_arm[arm_name])")
        println("  c0 positive control ratio = $(round(control_arm[arm_name]; digits = 5))   " *
                "(Phase 11's global analogue recovered at 0.157)")

        # --- The per-r1-bin breakdown. Non-uniformity is the EXPECTED shape (Pitfall 2). -------
        println()
        println(rpad("r1 bin", 20), rpad("n", 8), rpad("ridge RMSE", 14),
                rpad("prior RMSE", 14), rpad("ratio", 10))
        r1_tst = y_tst
        ratios = Float64[]
        for k in 1:(length(bin_edges) - 1)
            lo, hi   = bin_edges[k], bin_edges[k + 1]
            last_bin = k == length(bin_edges) - 1
            sel = findall(x -> (x >= lo) & (x < hi + (last_bin ? 1e-9 : 0.0)), r1_tst)
            if isempty(sel)
                push!(ratios, NaN)
                continue
            end
            r  = _p12ell_rmse(pred[sel], y_tst[sel])
            p0 = sqrt(mean(abs2, y_tst[sel] .- r1_bar_fit))
            push!(ratios, p0 > 0 ? r / p0 : NaN)
            println(rpad("[$(round(lo; digits = 3)), $(round(hi; digits = 3)))", 20),
                    rpad(length(sel), 8), rpad(round(r; sigdigits = 6), 14),
                    rpad(round(p0; sigdigits = 6), 14),
                    rpad(round(last(ratios); digits = 5), 10))
        end
        bin_arm[arm_name] = ratios
    end

    # --- The reading, printed but NOT gated ---------------------------------------------------
    primary_ratio   = ratio_arm[primary_arm_name]
    primary_control = control_arm[primary_arm_name]
    control_live    = primary_control < 1.0
    println()
    println("="^78)
    println("D-08 IDENTIFIABILITY -- REPORTED, GATES NOTHING")
    println("="^78)
    println("  ratio_raw        = $(round(ratio_arm["raw"]; digits = 5))")
    println("  ratio_with_moran = $(round(primary_ratio; digits = 5))   (PRIMARY)")
    println("  control_c0_ratio = $(round(primary_control; digits = 5))   " *
            "live = $control_live")
    if !control_live
        println("  !! The c0 control is not below 1.0. BOTH arms are UNINFORMATIVE and no " *
                "identifiability verdict may be recorded from this run.")
    end
    println("  ridge_residual_shrinkage = $(round(shrink_arm[primary_arm_name]; digits = 5))   " *
            "vacuous = $(vacuous_arm[primary_arm_name])")
    println("  NOTE: read this beside 12-18's `shrinkage` (a flow's posterior width) as a PAIR, " *
            "under the rule fixed in 12-13:131-152. They are DIFFERENT ESTIMATORS. A disagreement " *
            "is a FINDING and is reported as one -- never averaged, never reconciled.")

    elapsed = (time() - t_start) / 60
    tmp = ELL_RIDGE_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version            = 1,
        generated                 = string(Dates.now(Dates.UTC)) * "Z",
        elapsed_min               = elapsed,
        ratio_raw                 = ratio_arm["raw"],
        ratio_with_moran          = ratio_arm["with_moran"],
        ratio_by_arm              = ratio_arm,
        primary_arm               = primary_arm_name,
        control_c0_ratio          = control_arm[primary_arm_name],
        control_c0_ratio_by_arm   = control_arm,
        control_live              = control_live,
        ridge_residual_shrinkage  = shrink_arm[primary_arm_name],
        ridge_residual_shrinkage_by_arm = shrink_arm,
        vacuous                   = vacuous_arm[primary_arm_name],
        vacuous_by_arm            = vacuous_arm,
        vacuous_shrinkage_floor   = P12_VACUOUS_SHRINKAGE_FLOOR,
        r1_prior_sd_fit           = r1_prior_sd_fit,
        r1_prior_sd_analytic      = r1_prior_sd_analytic,
        r1_prior_mean_fit         = r1_bar_fit,
        moran_i_summary           = moran_i_summary,
        r1_bin_ratios             = bin_arm,
        r1_bin_edges              = bin_edges,
        r1_bin_probs              = collect(ELL_BIN_PROBS),
        n_total                   = N,
        n_fit                     = length(fit_i),
        n_val                     = length(val_i),
        n_test                    = length(test_i),
        test_fraction             = ELL_TEST_FRACTION,
        val_fraction              = ELL_VAL_FRACTION,
        ridge_grid                = ELL_RIDGE_GRID,
        selected_penalty          = penalty_arm,
        pool_dir                  = dir,
        imsize                    = P12_STAGE1_IMSIZE,
        imsize_tag                = ELL_IMSIZE_TAG,
        arm                       = ELL_ARM,
        zscore_arm                = P12_ZSCORE_ARM,
        theta_row_r1              = idx_r1,
        theta_row_c0              = idx_c0,
        counter                   = P12_ELL_RIDGE_COUNTER,
        shrinkage_naming_note     = "This key is ridge_residual_shrinkage = residual_sd/prior_sd " *
                                    "for a LINEAR POINT PREDICTOR. 12-18 stores `shrinkage` = " *
                                    "post_sd/prior_sd for a NEURAL FLOW's posterior. Two " *
                                    "estimators of the same idea by different machinery, named " *
                                    "apart on purpose; they coincide only under a well-specified " *
                                    "linear-Gaussian model. Read as a PAIR under 12-13:131-152. " *
                                    "NEVER average them; a disagreement is a FINDING. No " *
                                    "prior-echo ceiling is recorded: r1 is INFERRED, not " *
                                    "conditioned, so P12_R1_MAX/P12_R1_MIN is a ratio of prior " *
                                    "SUPPORT BOUNDS and is not commensurable with either.",
        caption                   = "Phase-12 D-08 correlation-length identifiability: ridge on r1 " *
                                    "from the 128 raw summary rows (arm `raw`) and from those rows " *
                                    "PLUS Moran's I under the lattice's own 4-neighbour adjacency " *
                                    "(arm `with_moran`, PRIMARY), against the EMPIRICAL fit-split " *
                                    "prior-mean baseline, with a c0 global-level positive control. " *
                                    "NO row is excluded -- r1 is inferred, not conditioned. " *
                                    "Diagnostic; gates nothing.",
    )
    let d = JLD2.load(tmp)
        @assert haskey(d, "ratio_raw")                "artifact integrity: ratio_raw missing"
        @assert haskey(d, "ratio_with_moran")         "artifact integrity: ratio_with_moran missing"
        @assert haskey(d, "control_c0_ratio")         "artifact integrity: control_c0_ratio missing"
        @assert haskey(d, "ridge_residual_shrinkage") "artifact integrity: ridge_residual_shrinkage missing"
        @assert haskey(d, "vacuous")                  "artifact integrity: vacuous missing"
        @assert haskey(d, "moran_i_summary")          "artifact integrity: moran_i_summary missing"
        @assert haskey(d, "r1_bin_ratios")            "artifact integrity: r1_bin_ratios missing"
        @assert !haskey(d, "shrinkage") "artifact integrity: a key named `shrinkage` must NOT " *
            "exist here -- that name belongs to 12-18's posterior-width estimator"
        @assert endswith(d["caption"], "Diagnostic; gates nothing.")
    end
    mv(tmp, ELL_RIDGE_REPORT_PATH; force = true)

    println()
    println("wrote ", ELL_RIDGE_REPORT_PATH, "  (", round(elapsed; digits = 2), " min)")
    return nothing
end

isdefined(@__MODULE__, :P12_ELL_RIDGE_LOAD_ONLY) || main()
