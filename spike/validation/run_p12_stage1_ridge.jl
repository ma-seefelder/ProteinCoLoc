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

# spike/validation/run_p12_stage1_ridge.jl --- the D-12 STAGE-1 GATE: leave-region-out ridge
# across the r1 ladder, run BEFORE any network is trained.
#
# THE QUESTION: can ANY estimator recover a held-out region's drawn lattice value from the OTHER
# 63 regions' summary rows? The whole phase rests on the answer. If the other regions carry no
# information about the one held out, there is no spatial borrowing to learn, and a network
# trained on the same data would at best reproduce the per-region prior more expensively.
#
# THE EXCLUSION, AND IT IS THE WHOLE POINT OF THE TEST. Region r's own CONTINUOUS row AND its own
# MASK row are zeroed -- through `mask_regions`, the SAME encoding `encode_d01` already produces
# for a patch that falls below the >=15-surviving-pixel floor, so "region unobserved" is expressed
# once and in the vocabulary the shipped encoder already speaks. A probe that admits region r's
# own row measures the TRIVIAL IDENTITY (the summary row IS a noisy read of that region's rho),
# not borrowing, and would pass at every rung for a reason that has nothing to do with the lattice
# prior. This is `run_p11_recovery.jl`'s row-129 exclusion transposed onto the per-region question.
#
# THE BASELINE: predict the PER-REGION PRIOR MEAN, computed EMPIRICALLY from the realized draws in
# the FIT split of each rung -- never from an analytic marginal. Under the D-05 copula every
# region's marginal is the same by construction, so an analytic baseline would be defensible; it
# is deliberately not used, because the comparison must not be flattered by any mismatch between
# the nominal and the realized draw distribution. Ridge beating this baseline means the
# NEIGHBOURING regions carry recoverable information about the held-out one; ridge merely matching
# it means they do not.
#
# Ridge is the right tool here precisely because it is weak and transparent: closed-form, no
# tuning theatre, no capacity to memorise. If a 128-predictor linear map cannot beat the
# per-region prior, the honest reading is that the LINEAR borrowing signal is absent -- and the
# gate is deliberately asked in the weakest form that can still be believed.
#
# THE POSITIVE CONTROL: the SAME ridge, the SAME split, the SAME penalty grid, the SAME target,
# with row r INCLUDED. It must land below `P12_STAGE1_CONTROL_CEILING`. WITHOUT IT, "no borrowing"
# and "this script is wrong" are the same observation: a broken standardizer, a target/predictor
# misalignment or a pool carrying the wrong r1 would each produce a flat null that looks exactly
# like a scientific finding. If the control fails, the run is UNINFORMATIVE and NO verdict --
# neither PROCEED nor DESCOPE -- may be recorded from it.
#
# THIS RUNNER **GATES**. It is D-12 Stage 1, and it is the one Phase-12 runner whose header must
# not claim to be a diagnostic: `P12_STAGE1_RATIO_CEILING`, `P12_STAGE1_CONTROL_CEILING` and
# `P12_STAGE1_MIN_RUNGS` are members of `P12_GATING_CONSTANTS`, frozen in Tier 1 before any
# Phase-12 field existed. Its verdict is recorded in `12-STAGE1-VERDICT.md` BEFORE the ~1.8 h
# mini-spike and the ~2.5 h full iteration are spent -- that ordering IS the value of the gate,
# because a probe run after the spend is not a gate. It therefore does NOT throw on failure: a
# DESCOPE is a legitimate, pre-registered outcome, and an exception would make the artifact
# unwritable at exactly the moment it matters most.
#
# THE GATED TARGET IS THE **GAUSSIAN** FIELD VALUE `z_field[r]` (R-2) -- the sampled parameter,
# atom-free. The rho-space target `rho_field[r]` is run as a SECONDARY REPORTED ARM and gates
# nothing: in rho space the elementwise `ghat` clamp puts a measured 6.79 % atom mass on every
# region, which compresses the target's variance and therefore FLATTERS any predictor's ratio. A
# reader must not take the rho column as the gate.
#
# ONE POOL PER RUNG, WITH `r1` PINNED, AND THE PIN IS LOAD-BEARING TWICE OVER. The gate asks *at
# what correlation length does borrowing exist*, so a pooled r1 draw would AVERAGE the answer away:
# under the copula every region's marginal is fixed regardless of r1, so the within-image contrast
# is controlled entirely by the correlation length (12-RESEARCH Pitfall 2). And `r1` enters the
# CONTENT HASH as well as the draw (`p12_generating_config`'s `r1_pin`): without that, `n`, `arm`
# and the r1 prior bounds are identical across rungs, all five hash to ONE directory,
# resume-by-skip serves rung 1's data to rungs 2-5, and the realized-count assertion below still
# passes because the COUNT is right while the CONTENT is one rung repeated five times.
#
# `imsize_tag` IS PASSED EXPLICITLY, AND IT IS NOT OPTIONAL. The image-size sampler is a closure
# and cannot be hashed, so a pinned-size pool WITHOUT a tag hashes into the F5-mixture directory
# and is indistinguishable from one (12-09 §6.1). The tag is DERIVED from `P12_STAGE1_IMSIZE`
# rather than typed, so it cannot drift from the size it names.
#
# THE FIVE RUNGS SHARE THEIR RANDOM KEYS, BY CONSTRUCTION, AND THAT IS A FEATURE WORTH RECORDING.
# `p12_datagen_rng` is keyed by the GLOBAL INDEX alone, and pinning r1 SKIPS the r1 draw
# (`p12_prior.jl:187`), so sample i of every rung consumes the same stream from the field draw
# onwards. The rungs are therefore COMMON-RANDOM-NUMBER paired rather than independent: the ladder
# comparison is variance-reduced, and no rung-to-rung difference can be a seed artifact. It also
# means the five rungs are NOT five independent experiments, which is how their spread must be read.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; reads the frozen pre-registration and
# GENERATES ITS OWN pools under `spike/data/cache/p12/` through the one datagen path. It trains no
# network, mutates no constant, appends no Tier-2 constant, touches no `src/` file, and does not
# read, write or regenerate Phase 11's pool.
#
# SEEDING. The split is a DETERMINISTIC tail block with no shuffle, exactly as the analog, so the
# ridge half of this run consumes NO random draw at all. `P12_STAGE1_COUNTER` is recorded in the
# artifact as the audit-trail tag for this activity rather than consumed as a key word -- the same
# distinction `spike/data/p12_generate.jl:129-132` draws for its own datagen counter. Had a
# shuffle been needed, `p12_rng(P12_STAGE1_COUNTER)` is the reserved stream for it. The pools ride
# the datagen stream `P12_DEV_SEED xor P12_DATAGEN_SALT`, keyed per global index.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p12_stage1_ridge.jl

using JLD2
using Dates
using Statistics
using LinearAlgebra

# ORDER: the Tier-1 pre-registration first (it binds the ladder, the sizes and the three gate
# constants), then the datagen path (which transitively pulls the contract, the Phase-12 prior,
# the forward model, the lattice arithmetic, the encoder and the sharded cache), then the
# architecture file for `mask_regions`. Only Phase-12-local sentinel names are guarded on; no
# foreign pre-registration is loaded into this process, so the include-guard poisoning class
# (p12_consts.jl BINDS :P11_DEV_SEED and :P13_DEV_SEED) cannot arise here.
isdefined(@__MODULE__, :P12_DEV_SEED)     || include(joinpath(@__DIR__, "p12_consts.jl"))
isdefined(@__MODULE__, :generate_p12_pool) || include(joinpath(@__DIR__, "..", "data", "p12_generate.jl"))
isdefined(@__MODULE__, :mask_regions)      || include(joinpath(@__DIR__, "..", "npe", "p12_architecture.jl"))

const STAGE1_REPORT_PATH = joinpath(@__DIR__, "p12_stage1_report.jld2")

# Ridge penalties swept on a validation slice carved out of TRAIN (never out of TEST). These are
# REGULARIZATION STRENGTHS, not thresholds: nothing here decides a pass or a fail, and the value
# selected is reported per region so a reader can see which one the validation split chose.
const RIDGE_GRID = [1e-6, 1e-4, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0]

# The prior arm of the Stage-1 ladder. CAR is the phase's default arm; the CAR-vs-GP comparison is
# 12-15's question, not this gate's, and running the gate on both arms would double the datagen
# spend to answer a question the gate does not ask.
const STAGE1_ARM = :car

# The pinned-size tag that keeps these pools out of the F5-mixture directory. DERIVED from the
# frozen size, never typed, so the tag cannot name a size the pool does not have.
const STAGE1_IMSIZE_TAG =
    Symbol("pinned_", P12_STAGE1_IMSIZE[1], "x", P12_STAGE1_IMSIZE[2])

# The fraction of TRAIN held out for penalty selection. A split fraction, not a threshold; taken
# unchanged from `run_p11_recovery.jl:133` so the two probes' splits are directly comparable.
const STAGE1_VAL_FRACTION = 0.2

"""
    _p12s1_fit_standardizer(X) -> (mu, sd)

Column-wise (per predictor row) mean and sd over TRAINING columns only. Rows with zero variance
get `sd = 1` so they pass through as constants instead of producing NaNs -- which is exactly what
the near-constant present-mask block needs, AND what the two rows `mask_regions` has just zeroed
need: a zeroed row has zero variance by construction, and without this guard the masked design
would be all-NaN and the null would be an artifact of a division rather than a finding.

Ported unchanged from `spike/validation/run_p11_recovery.jl:75-82`; the `1e-12` is that file's
numerical zero-variance guard, not a bar this file invents.
"""
function _p12s1_fit_standardizer(X::AbstractMatrix)
    mu = vec(mean(X; dims = 2))
    sd = vec(std(X; dims = 2))
    @inbounds for i in eachindex(sd)
        (isfinite(sd[i]) && sd[i] > 1e-12) || (sd[i] = 1.0)
    end
    return mu, sd
end

_p12s1_apply_standardizer(X, mu, sd) = (X .- mu) ./ sd

_p12s1_predict(Xs, w, b) = (Xs' * w) .+ b

_p12s1_rmse(pred, truth) = sqrt(mean(abs2, pred .- truth))

"""
    _p12s1_gram(Xfit) -> Symmetric

`Xfit * Xfit'`, the penalty-independent half of the closed-form ridge normal equations, computed
ONCE per design instead of once per (penalty, target). `Xfit` is (features x samples).

THIS IS AN ARITHMETIC REARRANGEMENT OF `run_p11_recovery.jl:92-98`, NOT A DIFFERENT ESTIMATOR.
That file solves `Symmetric(Xs*Xs' + penalty*I) w = Xs*(y - mean(y))` inside the penalty loop; the
Gram term does not depend on the penalty or on the target, and this gate fits 64 regions x 2
targets x 2 designs x 8 penalties per rung, so recomputing it would be ~2000 redundant 128x128
accumulations per rung for identical numbers.
"""
_p12s1_gram(Xfit::AbstractMatrix) = Symmetric(Xfit * Xfit')

"""
    _p12s1_select_and_fit(A, Xfit, yfit, Xval, yval, Xtst) -> (pred, best_penalty)

Sweep `RIDGE_GRID`, select the penalty on the VALIDATION block, and predict on TEST with the
selected fit. The test block is touched exactly ONCE, at the end, by construction -- there is no
path through this function on which a test-set statistic can reach the selection.

The selected weights are the ones already computed for the winning penalty, which is identical to
refitting on FIT at that penalty (same closed form, same data), so the analog's explicit refit is
omitted rather than changed.

EVERY SOLVE CARRIES A PENALTY, AND THAT IS NOT AN OPTIMISATION DETAIL. In the masked design two
rows -- region r's continuous row and its mask row -- are EXACTLY zero by construction, so the
Gram matrix is EXACTLY rank-deficient (rank 126 of 128) no matter how many fit samples there are.
An unpenalised `A \\ rhs` anywhere in this function throws `SingularException` on the arm the gate
is actually about, so the seed of the search is the first grid penalty rather than the OLS fit.
"""
function _p12s1_select_and_fit(A, Xfit, yfit, Xval, yval, Xtst)
    b   = mean(yfit)                       # unpenalised intercept, from FIT only
    rhs = Xfit * (yfit .- b)
    best_p = RIDGE_GRID[1]
    best_w = (A + best_p * I) \ rhs
    best_v = _p12s1_rmse(_p12s1_predict(Xval, best_w, b), yval)
    for p in RIDGE_GRID
        w = (A + p * I) \ rhs
        v = _p12s1_rmse(_p12s1_predict(Xval, w, b), yval)   # VALIDATION, never TEST
        if v < best_v
            best_p, best_v, best_w = p, v, w
        end
    end
    return _p12s1_predict(Xtst, best_w, b), best_p
end

"""
    _p12s1_radial_r2(per_region) -> Float64

The R² of an OLS regression of the 64 per-region ratios on `radial_basis(P12_G)` = `[1 radius]`.

REPORTED, NEVER GATED. If per-region recoverability is itself strongly radial, S-4's chromatic
gradient is a live candidate for doing the borrowing rather than the lattice prior, and knowing
that BEFORE the mini-spike is worth the few lines it costs. `radial_basis` is reused rather than
re-derived: a second derivation of the lattice geometry is a second chance to transpose it.
"""
function _p12s1_radial_r2(per_region::AbstractVector)
    X   = radial_basis(P12_G)
    y   = collect(float.(per_region))
    β   = X \ y
    res = y .- X * β
    sst = sum(abs2, y .- mean(y))
    return sst <= 0 ? NaN : 1.0 - sum(abs2, res) / sst
end

"""
    _p12s1_mean_within_image_sd(M) -> Float64

Mean over datasets of the WITHIN-IMAGE across-region sd of a `64 x n` per-region quantity.

A MEASUREMENT, WITH NO BAR ATTACHED. Applied to the 64 continuous summary rows it is directly
comparable to 12-RESEARCH Pitfall 2's per-draw noise sd in raw correlation units; applied to the
true rho lattice it is the injected per-region contrast the field actually carries at that rung.
The comparison between the two and the research figure is made in `12-STAGE1-VERDICT.md`, in
prose, against a number quoted from the research with its provenance -- deliberately NOT as a
literal in this file, which carries no bar of its own beyond the three Tier-1 gate constants.
"""
_p12s1_mean_within_image_sd(M::AbstractMatrix) = mean(vec(std(M; dims = 1)))

# =============================================================================================
# THE LIVENESS DIAGNOSTIC -- ADDED AFTER THE FIRST RUN, REPORTED, AND IT DECIDES NOTHING
# =============================================================================================
#
# WHEN AND WHY IT WAS ADDED, RECORDED HERE BECAUSE THE TIMING IS THE POINT. The first run of this
# runner produced `control_live = false`: the positive control cleared the ratio ceiling at every
# rung but did not clear `P12_STAGE1_CONTROL_CEILING` at the three SHORTEST correlation lengths.
# The frozen interpretation rule says a failed control makes the run UNINFORMATIVE and names three
# likely causes -- a broken standardizer, a target/predictor misalignment, a pool with the wrong r1
# -- and directs the phase's single iteration allowance at FIXING THE HARNESS. That instruction
# rests on a premise, and the premise is checkable. These three quantities check it.
#
# THEY CHANGE NO THRESHOLD, DECIDE NOTHING, AND DO NOT ENTER `stage1_pass`, `borrowing_ok` OR
# `control_live`. The gate is evaluated from exactly the constants and the arms it was evaluated
# from before they existed, and the gate values are bit-identical across the two runs. Adding a
# REPORTED diagnostic after a result is legitimate; moving a bar after a result is not, and nothing
# below moves one.
#
#   `_p12s1_own_row_bound`  -- sqrt(1 - corr(row_r, z_r)^2), the RMSE ratio of the best possible
#       LINEAR predictor of z_field[r] from region r's OWN summary row. It is an information limit
#       of the DATA, computed with no ridge, no split and no standardizer, so it cannot be produced
#       by any defect in the estimator it is used to judge. If the measured control sits AT this
#       bound, the control is optimal and the harness is live; a broken harness lands ABOVE it.
#   `ownrow_ratio`          -- the same ridge machinery restricted to region r's own two rows. Its
#       agreement with the bound is the executable form of the same argument.
#   `global_control_ratio`  -- the full 128-row ridge predicting the GLOBAL field level, which is
#       the Phase-11 benchmark 12-CONTEXT S-2 names (global rho_true recovered at ratio 0.157, an
#       84 % error reduction). It is the one number in this phase directly comparable to a
#       measurement made on a harness already known to work.

"""
    _p12s1_own_row_bound(srow, y) -> Float64

`sqrt(1 - corr(srow, y)^2)`: the RMSE-over-prior-sd ratio of the BEST LINEAR predictor of `y` from
the single predictor `srow`. No ridge, no penalty, no split, no standardizer -- a property of the
DATA, so it is independent of every mechanism it is used to audit.
"""
function _p12s1_own_row_bound(srow::AbstractVector, y::AbstractVector)
    c = cor(collect(float.(srow)), collect(float.(y)))
    return isfinite(c) ? sqrt(max(0.0, 1 - c^2)) : NaN
end

function main()
    println("="^78)
    println("D-12 STAGE-1 GATE -- leave-region-out ridge across the r1 ladder")
    println("Region r's OWN continuous row AND mask row are EXCLUDED. This runner GATES.")
    println("="^78)
    t_start = time()

    rungs       = collect(float.(P12_R1_LADDER))
    n_rungs     = length(rungs)
    n_per_rung  = P12_STAGE1_N ÷ n_rungs
    G           = P12_G
    nreg        = G^2

    println("ladder      = $(rungs)")
    println("n per rung  = $n_per_rung   (P12_STAGE1_N = $(P12_STAGE1_N) over $n_rungs rungs)")
    println("arm         = $STAGE1_ARM   imsize = $(P12_STAGE1_IMSIZE)  tag = $STAGE1_IMSIZE_TAG")
    println("gate        = ratio <= $(P12_STAGE1_RATIO_CEILING) at >= $(P12_STAGE1_MIN_RUNGS) rung(s)" *
            "  AND  max control ratio <= $(P12_STAGE1_CONTROL_CEILING)")

    # --- The five pools must be five DIRECTORIES. Asserted BEFORE anything is generated. --------
    # `p12_pool_dir` resolves through `open_or_invalidate`, which CREATES the directory it
    # resolves, so this call is a hash probe and not a presence test (12-09 §5, `p12_pool_complete`
    # is the presence test). What it proves is the only thing at stake here: that `r1_pin` reaches
    # the content hash, so the five rungs cannot collide into one directory.
    pool_dirs = [p12_pool_dir(n_per_rung; arm = STAGE1_ARM, r1 = rung,
                              imsize_tag = STAGE1_IMSIZE_TAG) for rung in rungs]
    @assert length(unique(pool_dirs)) == n_rungs """
        run_p12_stage1_ridge: the $n_rungs r1 rungs resolved to $(length(unique(pool_dirs)))
        distinct pool directories. The r1 PIN is not reaching the content hash, so resume-by-skip
        would serve one rung's data for all five and this gate would report a ladder that is one
        rung repeated. Refusing to generate."""

    ratio_mean          = Vector{Float64}(undef, n_rungs)
    ratio_min           = Vector{Float64}(undef, n_rungs)
    ratio_max           = Vector{Float64}(undef, n_rungs)
    control_ratio_mean  = Vector{Float64}(undef, n_rungs)
    control_ratio_max   = Vector{Float64}(undef, n_rungs)
    rho_arm_ratio_mean  = Vector{Float64}(undef, n_rungs)
    rho_arm_ratio_max   = Vector{Float64}(undef, n_rungs)
    rho_control_mean    = Vector{Float64}(undef, n_rungs)
    radial_r2           = Vector{Float64}(undef, n_rungs)
    within_image_sd_rows     = Vector{Float64}(undef, n_rungs)
    within_image_sd_rho      = Vector{Float64}(undef, n_rungs)
    within_image_sd_zfield   = Vector{Float64}(undef, n_rungs)
    baseline_rmse_mean       = Vector{Float64}(undef, n_rungs)
    ownrow_ratio_mean        = Vector{Float64}(undef, n_rungs)   # REPORTED liveness diagnostic
    own_row_bound_mean       = Vector{Float64}(undef, n_rungs)   # REPORTED liveness diagnostic
    global_control_ratio     = Vector{Float64}(undef, n_rungs)   # REPORTED liveness diagnostic
    ownrow_ratio_per_region  = Matrix{Float64}(undef, nreg, n_rungs)
    own_row_bound_per_region = Matrix{Float64}(undef, nreg, n_rungs)
    ratio_per_region         = Matrix{Float64}(undef, nreg, n_rungs)
    control_ratio_per_region = Matrix{Float64}(undef, nreg, n_rungs)
    rho_ratio_per_region     = Matrix{Float64}(undef, nreg, n_rungs)
    selected_penalty         = Matrix{Float64}(undef, nreg, n_rungs)
    selected_penalty_control = Matrix{Float64}(undef, nreg, n_rungs)
    n_test_per_rung          = Vector{Int}(undef, n_rungs)
    realized_dirs            = Vector{String}(undef, n_rungs)

    for (k, rung) in enumerate(rungs)
        println()
        println("#"^78)
        println("RUNG $k/$n_rungs   r1 = $rung  (PINNED)")
        println("#"^78)

        # Resume-by-skip lives inside `generate_p12_pool` and costs nothing when the shards are
        # already complete, so this is both "make sure it is there" and "generate it".
        dir = generate_p12_pool(n_per_rung; arm = STAGE1_ARM,
                                imsize_sampler = _ -> P12_STAGE1_IMSIZE,
                                r1 = rung, imsize_tag = STAGE1_IMSIZE_TAG)
        @assert dir == pool_dirs[k] "run_p12_stage1_ridge: generated pool landed in $dir, the " *
            "hash probe resolved $(pool_dirs[k])"
        @assert p12_pool_complete(dir, n_per_rung) "run_p12_stage1_ridge: pool at $dir is not " *
            "complete for n = $n_per_rung"
        realized_dirs[k] = dir

        pool = load_p12_pool(dir)
        S    = Float64.(pool.summary_min)          # 128 x n, RAW (never z-scored at rest)
        n    = size(S, 2)

        # --- STRUCTURAL aborts. None of these is a gate outcome. --------------------------------
        @assert n == n_per_rung "run_p12_stage1_ridge: pool at $dir holds $n samples, expected $n_per_rung"
        @assert size(S, 1) == P12_SUMMARY_MIN_DIM "run_p12_stage1_ridge: predictors must be the " *
            "$(P12_SUMMARY_MIN_DIM) encoded summary rows, got $(size(S, 1))"
        @assert size(pool.theta, 1) == P12_D_MINISPIKE "run_p12_stage1_ridge: theta has " *
            "$(size(pool.theta, 1)) rows, expected P12_D_MINISPIKE = $(P12_D_MINISPIKE)"
        @assert all(==(rung), pool.r1) """
            run_p12_stage1_ridge: rung $k asked for r1 = $rung but the loaded pool carries
            $(length(unique(pool.r1))) distinct value(s) in $(extrema(pool.r1)). The pin was
            accepted and not applied, or a foreign rung's directory was served."""
        @assert all(==(P12_STAGE1_IMSIZE), pool.imsize) "run_p12_stage1_ridge: pool at $dir " *
            "carries image sizes $(unique(pool.imsize)), expected the pinned $(P12_STAGE1_IMSIZE)"

        # The lattices reshape column-major into G² x n, which is the SAME `p12_idx` ordering the
        # summary rows use (`encode_d01` vecs the G×G matrix column-major). Region r of the summary
        # and row r of these matrices are therefore the same region, with no transpose anywhere.
        ZF = reshape(pool.z_field,   nreg, n)      # the GAUSSIAN field -- the GATED target (R-2)
        RF = reshape(pool.rho_field, nreg, n)      # the rho field -- SECONDARY, reported only

        # --- Leak-free split: deterministic tail block, no shuffle, exactly as the analog --------
        ntest  = round(Int, P12_STAGE1_TEST_FRACTION * n)
        test_i = (n - ntest + 1):n
        tr_all = 1:(n - ntest)
        nval   = round(Int, STAGE1_VAL_FRACTION * length(tr_all))
        val_i  = (length(tr_all) - nval + 1):length(tr_all)
        fit_i  = 1:(length(tr_all) - nval)
        @assert length(test_i) == P12_STAGE1_N_TEST_PER_RUNG_DERIVED """
            run_p12_stage1_ridge: rung $k realized n_test = $(length(test_i)), but
            P12_STAGE1_N_TEST_PER_RUNG_DERIVED = $(P12_STAGE1_N_TEST_PER_RUNG_DERIVED).
            P12_STAGE1_RATIO_CEILING was DERIVED from that n (Tier 1 §6: sd ~ 1/sqrt(2*n_test)),
            so a different n silently changes this gate's power. Refusing to continue."""
        n_test_per_rung[k] = length(test_i)
        println("split: fit = $(length(fit_i))  val = $(length(val_i))  test = $(length(test_i))")

        # --- The CONTROL design: row r INCLUDED, so it is the same for every r ------------------
        muC, sdC = _p12s1_fit_standardizer(view(S, :, fit_i))       # FIT-ONLY fit
        XfitC = _p12s1_apply_standardizer(view(S, :, fit_i), muC, sdC)
        XvalC = _p12s1_apply_standardizer(view(S, :, val_i), muC, sdC)
        XtstC = _p12s1_apply_standardizer(view(S, :, test_i), muC, sdC)
        AC    = _p12s1_gram(XfitC)

        rr      = Vector{Float64}(undef, nreg)   # gated arm: z_field, row r EXCLUDED
        rc      = Vector{Float64}(undef, nreg)   # positive control: z_field, row r INCLUDED
        rrho    = Vector{Float64}(undef, nreg)   # secondary arm: rho_field, row r EXCLUDED
        rrhoc   = Vector{Float64}(undef, nreg)   # secondary control
        basev   = Vector{Float64}(undef, nreg)
        pens    = Vector{Float64}(undef, nreg)
        pensC   = Vector{Float64}(undef, nreg)
        rown    = Vector{Float64}(undef, nreg)   # REPORTED: own-two-rows ridge
        rbound  = Vector{Float64}(undef, nreg)   # REPORTED: own-row information limit

        for r in 1:nreg
            # THE EXCLUSION. Both of region r's rows -- continuous and mask -- are zeroed, in the
            # one place this phase expresses "region unobserved".
            Sm = mask_regions(S, [r])
            mu, sd = _p12s1_fit_standardizer(view(Sm, :, fit_i))    # FIT-ONLY fit
            Xfit = _p12s1_apply_standardizer(view(Sm, :, fit_i), mu, sd)
            Xval = _p12s1_apply_standardizer(view(Sm, :, val_i), mu, sd)
            Xtst = _p12s1_apply_standardizer(view(Sm, :, test_i), mu, sd)
            A    = _p12s1_gram(Xfit)

            # --- gated arm: the GAUSSIAN field value at the held-out region ---------------------
            yfit = vec(ZF[r, fit_i]); yval = vec(ZF[r, val_i]); ytst = vec(ZF[r, test_i])
            pred, bp = _p12s1_select_and_fit(A, Xfit, yfit, Xval, yval, Xtst)
            # THE BASELINE: the per-region prior mean, EMPIRICAL, from the FIT split only.
            base    = sqrt(mean(abs2, ytst .- mean(yfit)))
            rr[r]   = _p12s1_rmse(pred, ytst) / base
            basev[r] = base
            pens[r] = bp

            # --- positive control: the same everything, row r INCLUDED -------------------------
            cpred, cbp = _p12s1_select_and_fit(AC, XfitC, yfit, XvalC, yval, XtstC)
            rc[r]    = _p12s1_rmse(cpred, ytst) / base
            pensC[r] = cbp

            # --- secondary REPORTED arm: the rho-space value, clamp atoms and all ---------------
            yfitρ = vec(RF[r, fit_i]); yvalρ = vec(RF[r, val_i]); ytstρ = vec(RF[r, test_i])
            predρ, _ = _p12s1_select_and_fit(A, Xfit, yfitρ, Xval, yvalρ, Xtst)
            baseρ    = sqrt(mean(abs2, ytstρ .- mean(yfitρ)))
            rrho[r]  = _p12s1_rmse(predρ, ytstρ) / baseρ
            cpredρ, _ = _p12s1_select_and_fit(AC, XfitC, yfitρ, XvalC, yvalρ, XtstC)
            rrhoc[r]  = _p12s1_rmse(cpredρ, ytstρ) / baseρ

            # --- REPORTED liveness diagnostic. Gates nothing; see the block above. -------------
            own_rows = [r, nreg + r]
            muO, sdO = _p12s1_fit_standardizer(view(S, own_rows, fit_i))
            XfitO = _p12s1_apply_standardizer(view(S, own_rows, fit_i), muO, sdO)
            XvalO = _p12s1_apply_standardizer(view(S, own_rows, val_i), muO, sdO)
            XtstO = _p12s1_apply_standardizer(view(S, own_rows, test_i), muO, sdO)
            opred, _ = _p12s1_select_and_fit(_p12s1_gram(XfitO), XfitO, yfit, XvalO, yval, XtstO)
            rown[r]   = _p12s1_rmse(opred, ytst) / base
            rbound[r] = _p12s1_own_row_bound(view(S, r, test_i), ytst)
        end

        ratio_per_region[:, k]         = rr
        control_ratio_per_region[:, k] = rc
        rho_ratio_per_region[:, k]     = rrho
        selected_penalty[:, k]         = pens
        selected_penalty_control[:, k] = pensC

        ratio_mean[k]         = mean(rr)
        ratio_min[k]          = minimum(rr)
        ratio_max[k]          = maximum(rr)
        control_ratio_mean[k] = mean(rc)
        control_ratio_max[k]  = maximum(rc)
        rho_arm_ratio_mean[k] = mean(rrho)
        rho_arm_ratio_max[k]  = maximum(rrho)
        rho_control_mean[k]   = mean(rrhoc)
        radial_r2[k]          = _p12s1_radial_r2(rr)
        baseline_rmse_mean[k] = mean(basev)

        ownrow_ratio_per_region[:, k]  = rown
        own_row_bound_per_region[:, k] = rbound
        ownrow_ratio_mean[k]  = mean(rown)
        own_row_bound_mean[k] = mean(rbound)

        # The GLOBAL field level through the same full-128 design: the Phase-11-comparable control
        # (12-CONTEXT S-2 records global rho_true recovered at ratio 0.157). REPORTED, not gated.
        gtar  = vec(mean(ZF; dims = 1))
        gpred, _ = _p12s1_select_and_fit(AC, XfitC, gtar[fit_i], XvalC, gtar[val_i], XtstC)
        global_control_ratio[k] =
            _p12s1_rmse(gpred, gtar[test_i]) / sqrt(mean(abs2, gtar[test_i] .- mean(gtar[fit_i])))

        # Measured contrasts, for the noise-floor reading. No bar is applied to them here.
        within_image_sd_rows[k]   = _p12s1_mean_within_image_sd(view(S, 1:nreg, test_i))
        within_image_sd_rho[k]    = _p12s1_mean_within_image_sd(view(RF, :, test_i))
        within_image_sd_zfield[k] = _p12s1_mean_within_image_sd(view(ZF, :, test_i))

        println()
        println(rpad("r1", 8), rpad("ratio_mean", 13), rpad("ratio_min", 12), rpad("ratio_max", 12),
                rpad("control", 12), rpad("rho_arm", 12), rpad("radial_R2", 11))
        println(rpad(rung, 8), rpad(round(ratio_mean[k]; digits = 5), 13),
                rpad(round(ratio_min[k]; digits = 5), 12), rpad(round(ratio_max[k]; digits = 5), 12),
                rpad(round(control_ratio_mean[k]; digits = 5), 12),
                rpad(round(rho_arm_ratio_mean[k]; digits = 5), 12),
                rpad(round(radial_r2[k]; digits = 5), 11))
        println("  within-image across-region sd: summary rows = " *
                "$(round(within_image_sd_rows[k]; digits = 5))   " *
                "rho field = $(round(within_image_sd_rho[k]; digits = 5))   " *
                "z field = $(round(within_image_sd_zfield[k]; digits = 5))")
        println("  per-rung PASS against the ceiling: " *
                "$(ratio_mean[k] <= P12_STAGE1_RATIO_CEILING)")
        println("  liveness (REPORTED, decides nothing): own-row ridge = " *
                "$(round(ownrow_ratio_mean[k]; digits = 5))  vs its information limit " *
                "$(round(own_row_bound_mean[k]; digits = 5))   " *
                "global-level control = $(round(global_control_ratio[k]; digits = 5))")
    end

    # =============================================================================================
    # THE GATE -- evaluated in ONE place, from the frozen constants, and never thrown on
    # =============================================================================================
    borrowing_ok = count(k -> ratio_mean[k] <= P12_STAGE1_RATIO_CEILING, 1:n_rungs) >=
                   P12_STAGE1_MIN_RUNGS
    control_live = maximum(control_ratio_mean) <= P12_STAGE1_CONTROL_CEILING
    stage1_pass  = borrowing_ok && control_live
    rungs_passing = [rungs[k] for k in 1:n_rungs if ratio_mean[k] <= P12_STAGE1_RATIO_CEILING]

    # THE BANNER HAS THREE STATES, NOT TWO, AND THE THIRD IS NOT A HEDGE. `stage1_pass` is an AND
    # over two components, so `false` alone does not say WHICH half decided -- and the frozen
    # interpretation rule is explicit that a failed CONTROL makes the run UNINFORMATIVE, with
    # NEITHER a PROCEED nor a DESCOPE recordable from it. A two-state banner would therefore print
    # "DESCOPE" for a run from which a descope may not be read, and the banner is what a human sees.
    banner = !control_live ? "NO VERDICT — the positive control did not clear its ceiling" :
             stage1_pass   ? "PASS (PROCEED)" : "DESCOPE"
    println()
    println("="^78)
    println("D-12 STAGE-1 GATE — ", banner)
    println("="^78)
    println("  borrowing_ok = $borrowing_ok   ($(length(rungs_passing)) of $n_rungs rungs at or " *
            "below $(P12_STAGE1_RATIO_CEILING); $(P12_STAGE1_MIN_RUNGS) required; " *
            "passing rungs = $rungs_passing)")
    println("  control_live = $control_live   (max control ratio " *
            "$(round(maximum(control_ratio_mean); digits = 5)) against a ceiling of " *
            "$(P12_STAGE1_CONTROL_CEILING))")
    if !control_live
        println("  !! Neither a PROCEED nor a DESCOPE may be recorded from this run under the " *
                "frozen interpretation rule. Read the liveness diagnostic below BEFORE concluding " *
                "the harness is broken: if the control sits at its own-row information limit, the " *
                "ceiling is unreachable rather than the harness wrong, and which of those it is " *
                "is a pre-registration question, not a runtime one.")
        println("  own-row ridge  = $(round.(ownrow_ratio_mean; digits = 5))")
        println("  its info limit = $(round.(own_row_bound_mean; digits = 5))")
        println("  global control = $(round.(global_control_ratio; digits = 5))   " *
                "(Phase 11 measured 0.157 for global rho_true)")
    end
    println("  ratio_mean  = $ratio_mean")
    println("  control     = $control_ratio_mean")
    println("  rho arm     = $rho_arm_ratio_mean   (REPORTED, not gated)")
    println("  radial R2   = $radial_r2            (REPORTED, not gated)")

    elapsed = (time() - t_start) / 60
    tmp = STAGE1_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version           = 1,
        generated                = string(Dates.now(Dates.UTC)) * "Z",
        elapsed_min              = elapsed,
        stage1_pass              = stage1_pass,
        borrowing_ok             = borrowing_ok,
        control_live             = control_live,
        rungs_passing            = rungs_passing,
        r1_ladder                = rungs,
        n_per_rung               = n_per_rung,
        n_test_per_rung          = n_test_per_rung,
        n_test_per_rung_derived  = P12_STAGE1_N_TEST_PER_RUNG_DERIVED,
        test_fraction            = P12_STAGE1_TEST_FRACTION,
        val_fraction             = STAGE1_VAL_FRACTION,
        ridge_grid               = RIDGE_GRID,
        selected_penalty         = selected_penalty,
        selected_penalty_control = selected_penalty_control,
        imsize                   = P12_STAGE1_IMSIZE,
        imsize_tag               = STAGE1_IMSIZE_TAG,
        arm                      = STAGE1_ARM,
        zscore_arm               = P12_ZSCORE_ARM,
        target                   = "z_field (Gaussian space, R-2)",
        ratio_mean               = ratio_mean,
        ratio_min                = ratio_min,
        ratio_max                = ratio_max,
        ratio_per_region         = ratio_per_region,
        control_ratio_mean       = control_ratio_mean,
        control_ratio_max        = control_ratio_max,
        control_ratio_per_region = control_ratio_per_region,
        rho_arm_ratio_mean       = rho_arm_ratio_mean,
        rho_arm_ratio_max        = rho_arm_ratio_max,
        rho_arm_control_mean     = rho_control_mean,
        rho_ratio_per_region     = rho_ratio_per_region,
        baseline_rmse_mean       = baseline_rmse_mean,
        within_image_sd_rows     = within_image_sd_rows,
        within_image_sd_rho      = within_image_sd_rho,
        within_image_sd_zfield   = within_image_sd_zfield,
        radial_r2                = radial_r2,
        ownrow_ratio_mean        = ownrow_ratio_mean,
        ownrow_ratio_per_region  = ownrow_ratio_per_region,
        own_row_bound_mean       = own_row_bound_mean,
        own_row_bound_per_region = own_row_bound_per_region,
        global_control_ratio     = global_control_ratio,
        liveness_diagnostic_note = "REPORTED, decides nothing, added after the first run when the " *
                                   "positive control missed P12_STAGE1_CONTROL_CEILING. " *
                                   "own_row_bound = sqrt(1 - corr(row_r, z_r)^2) is an information " *
                                   "limit of the DATA (no ridge, no split, no standardizer), so a " *
                                   "control sitting AT it is optimal rather than broken. " *
                                   "global_control_ratio is the Phase-11-comparable benchmark " *
                                   "(12-CONTEXT S-2: global rho_true recovered at 0.157). Neither " *
                                   "enters stage1_pass, borrowing_ok or control_live, and the gate " *
                                   "values are bit-identical to the run that preceded them.",
        pool_dirs                = realized_dirs,
        ratio_ceiling            = P12_STAGE1_RATIO_CEILING,
        control_ceiling          = P12_STAGE1_CONTROL_CEILING,
        min_rungs                = P12_STAGE1_MIN_RUNGS,
        counter                  = P12_STAGE1_COUNTER,
        caption                  = "D-12 STAGE 1: THIS RUNNER **GATES**. Leave-region-out ridge " *
                                   "across the five pinned-r1 rungs of P12_R1_LADDER, region r's " *
                                   "own continuous AND mask rows excluded via mask_regions, " *
                                   "against the EMPIRICAL per-region prior-mean baseline, with a " *
                                   "row-r-included positive control. THE GATED TARGET IS THE " *
                                   "GAUSSIAN-SPACE LATTICE VALUE z_field[r] (R-2, atom-free); the " *
                                   "rho-space arm is REPORTED ONLY and its ratio is flattered by " *
                                   "the elementwise ghat clamp's 6.79 % per-region atom mass. The " *
                                   "verdict this artifact supports is recorded in " *
                                   "12-STAGE1-VERDICT.md.",
    )
    let d = JLD2.load(tmp)
        @assert haskey(d, "stage1_pass")        "artifact integrity: stage1_pass missing"
        @assert haskey(d, "ratio_mean")         "artifact integrity: ratio_mean missing"
        @assert haskey(d, "control_ratio_mean") "artifact integrity: control_ratio_mean missing"
        @assert haskey(d, "ratio_per_region")   "artifact integrity: ratio_per_region missing"
    end
    mv(tmp, STAGE1_REPORT_PATH; force = true)

    println()
    println("wrote ", STAGE1_REPORT_PATH, "  (", round(elapsed; digits = 2), " min)")
    return nothing
end

isdefined(@__MODULE__, :P12_STAGE1_LOAD_ONLY) || main()
