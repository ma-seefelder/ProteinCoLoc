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

# spike/validation/run_p12_minispike.jl --- the D-12 Stage-2 MINI-SPIKE: three arms, matched r1.
#
# THE QUESTION: at MATCHED induced lag-1 correlation, does a CAR or a GP lattice prior yield better
# per-region recovery -- and does either beat a NEUTRALIZED prior?
#
# THE HEAD: K = P12_K_DEV = 63, FULL RANK on the deviation space (D-04). Chosen so that a CAR loss
# cannot be a head artifact: if the head could not represent the deviation space, a arm difference
# might be reporting the truncation rather than the prior. Full rank removes that explanation before
# it can be offered.
#
# THE PARAMETRIZATION: both arms are drawn from `P12_R1_PRIOR` on INDUCED r1, never on alpha or ell
# (R-8). That is what "matched" means here -- the two kernels are compared at equal induced
# correlation rather than at equal nominal hyperparameter, which would compare different fields.
#
# THE SCORING SPACE: Gaussian-space DCT coefficients, ATOM-FREE (R-2). The rho-space arm carries the
# elementwise `ghat` clamp's atom mass and would flatter every ratio.
#
# SCOPE: THIS RUN **SELECTS** THE PRIOR under the rule fixed below, and **gates nothing** else --
# every other quantity it reports (coverage, c0 split, r1 shrinkage, the truncation curve) is
# REPORTED. It is NOT the D-12 Stage-1 gate; that already ran in 12-11 and its verdict is PROCEED.
#
# TWO WALL-CLOCK BUDGETS, AND THEY ARE **NOT** THE SAME BUDGET DESPITE CARRYING THE SAME NUMBER.
#   `P12_DATAGEN_WALLCLOCK_CEILING_MIN   = 150`  (p12_consts.jl:476) -- enforced INSIDE
#       `generate_p12_pool`, PER CALL, by projection. Each arm's pool is judged on its own; the
#       three pools' SUM is never a quantity anything compares to 150.
#   `P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150`  (p12_consts.jl:477) -- THIS RUN'S TOTAL, and until
#       this file existed NOTHING IN THE REPOSITORY READ IT ON EITHER AXIS: no runner applied it and
#       no test asserted its value. It was a number in a file. `_p12ms_check_budget` below is its
#       first enforcer, and a value assertion was added to `test_p12_consts.jl` in the same commit
#       so the constant is guarded on both axes like its siblings.
# Conflating the two is a live hazard precisely BECAUSE both read 150: the error is invisible at
# 150-vs-150 and would be obvious at 150-vs-90. Reported separately, each against its own constant.
#
# THE ESCAPES THAT ARE FORBIDDEN, NAMED SO THEY CANNOT BE TAKEN QUIETLY: shrinking
# `P12_MINISPIKE_N`, dropping an arm, or downgrading the image size to fit the budget. THREE ARMS IS
# THE PRE-REGISTERED DESIGN AND A TWO-ARM RUN IS A DIFFERENT EXPERIMENT. If the budget is exceeded
# the run throws with the literal `BLOCKER` and the shortfall is reported as one.
#
# THE HELD-OUT BLOCK IS CARVED FROM THE ONE POOL, NEVER GENERATED AS A SECOND POOL (CONVENTIONS
# C-04). `generate_p12_sample(idx)` keys its RNG on the GLOBAL INDEX ALONE and `generate_p12_pool(n)`
# always generates `1:n`, so two pools at the same arm/imsize share samples `1:min(n,m)`
# byte-identically -- a "separate scoring pool" is a SUPERSET of the training pool, never a held-out
# set. The scoring block is therefore the pool HEAD and the remainder is handed to training through
# `train_p12_npe`'s `pool_indices`, with disjointness ASSERTED against the bundle's own recorded
# `train_indices` AND `val_indices` -- both, because the val block drives early stopping and is
# contaminated for model selection while its indices genuinely are disjoint from training.
#
# PITFALL 5 -- THE ESTIMATOR EMITS **STANDARDIZED** theta, AND SCORING MUST UN-STANDARDIZE IT.
# `train_p12_npe` fits `theta_zt` and trains on `transform(theta_zt, theta)` (train_p12_npe.jl:458),
# then carries `theta_zt` in the bundle (:476) for exactly this purpose. EVERY read of a posterior
# draw therefore goes through `StatsBase.reconstruct(bundle.theta_zt, M)` FIRST. The contract is
# written down at `spike/npe/infer.jl:35-38` -- "reading a standardized draw as rho would be
# silently wrong" -- and honoured by `infer.jl:110,122` and `benchmark.jl:188`.
#
# THE FIRST RUN OF THIS FILE (2026-07-31T09:51Z) DID NOT DO IT, AND THE RESULT WAS HARVEST-READY
# GARBAGE: it scored standardized draws against the RAW `z_field`. The report was COMPLETE,
# internally consistent, and reconciled on every structural check -- index disjointness, per-region
# pooling, budget accounting -- because none of those checks look at the SCALE. What exposed it was
# arithmetic no assertion was making: the winning arm's RMSE (0.9898) sat 0.45 % below the RMSE of
# PREDICTING A CONSTANT ZERO (0.9943, the field's own marginal sd), the GP arm scored 60 % WORSE
# than that trivial predictor, and `r1_shrinkage` came to 4.298 -- which is just a standardized
# posterior sd of ~1.1 divided by the RAW prior sd of 0.2598. That run's `selected_prior`,
# `k_prod_recommended` and `n_low_recommended` were all invalid and were NEVER appended.
#
# PITFALL 6 -- THE theta ROWS ARE IN **SMOOTHNESS ORDER**, AND SCORING MUST INVERT THAT PERMUTATION.
# `p12_generate.jl:193` builds the coefficient rows as
#     c = p12_dct_vec(vec(draw.z_field); G = G)[p12_dct_order(G)]
# -- the flat DCT vector PERMUTED into `p12_dct_order` smoothness order, which is what makes "the
# first K modes" a defined object (p12_lattice.jl:400-416). Reconstructing a field therefore means
# SCATTERING each coefficient back to its flat mode index BEFORE the inverse DCT. `p12_region_field`
# (`spike/p12/result.jl:562`) is that inverse and is the ONE that exists; this file calls it and
# defines no second one.
#
# THE SECOND RUN OF THIS FILE APPLIED `p12_idct_vec` DIRECTLY TO THE PERMUTED VECTOR, and so scored
# every arm against a SCRAMBLED reconstruction of the truth. `p12_dct_order(8)` begins
# [1, 2, 9, 10, 3, 17, 11, 18, 4, 25, 19, 12, ...] and is emphatically NOT the identity. Measured on
# a real `sample_p12_prior` draw: the correct inverse recovers the field to 3.11e-15, while the
# naive one lands 3.72 away, with RMS error 1.263 against the field's OWN RMS of 0.970 -- a ratio of
# **1.302**. So A PERFECT POSTERIOR WOULD STILL HAVE SCORED ~1.30x THE TRIVIAL PREDICTOR, and the
# three arms were ranked on how well each net's coefficients happened to survive a scramble.
# (The ratio is 1.302 rather than sqrt(2) because the permutation maps low smoothness ranks onto
# lowish FLAT indices, so on a smooth CAR field it is a partial scramble rather than a full one.)
#
# WHY NO CHECK CAUGHT IT, AND WHY THAT IS THE SAME LESSON AS PITFALL 5 (12-STAGE1-VERDICT.md §7.1).
# A PERMUTATION IS ORTHOGONAL. Every structural property -- cardinality, per-region pooling to the
# reported scalar, index disjointness, budget reconciliation -- is preserved exactly, and Parseval
# means even the total energy is unchanged. `trivial_rmse` and `skill`, added after Pitfall 5, DID
# fire: skill came out negative for every arm. What they could not do is say WHICH wrong space.
#
# WHAT THE DEFECT DID AND DID NOT REACH, established by execution rather than by inspection:
#   AFFECTED  -- `rmse`, `skill`, `rmse_per_region`, `coverage`, `coverage_per_region`, and hence
#                ADMISSIBILITY, and hence the NONE-BEATS-ABLATION verdict.
#   UNAFFECTED -- `c0_rmse` / `c0_coverage`, because `p12_dct_order(G)[1] == 1` leaves the global
#                term at its own flat index (verified bit-identical on both sides);
#                `r1_shrinkage`, which reads row `end`; the truncation curve, which works in FLAT
#                index space throughout and never touches a theta row; and THE TRAINED BUNDLES,
#                because training consumes theta AS STORED and never reconstructs a field.
#
# SO THIS FILE NOW REPORTS `trivial_rmse` AND `skill` PER ARM. They are REPORTED and GATE NOTHING --
# no bar was added and no threshold moved, because inventing a floor after seeing a number is the
# move the two-tier pre-registration exists to forbid. But a skill number makes "the net learned
# nothing" VISIBLE in the artifact instead of invisible, which is the whole difference between the
# defect above being caught in three minutes and being appended to Tier 2 forever.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. Trains through 12-14's `train_p12_npe` and
# standardizes through its `standardize_p12`; it re-implements NEITHER, so no second training loop
# and no second standardizer fit can drift from the one surface. Touches no `src/` file.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p12_minispike.jl

using JLD2
using Dates
using Statistics
using LinearAlgebra
using StatsBase          # reconstruct -- the INVERSE theta-standardization (Pitfall 5, below)
using Random             # seed! -- the sampling seed, so the artifact is reproducible

isdefined(@__MODULE__, :P12_DEV_SEED)   || include(joinpath(@__DIR__, "p12_consts.jl"))
isdefined(@__MODULE__, :train_p12_npe)  || include(joinpath(@__DIR__, "..", "npe", "train_p12_npe.jl"))
# `p12_region_field` -- THE FIELD RECONSTRUCTION, REUSED RATHER THAN RE-IMPLEMENTED (Pitfall 6).
# This file previously built its own inverse inline and got it wrong; the repository already
# containing BOTH a right and a wrong reconstruction of the same object is how that happened.
isdefined(@__MODULE__, :p12_region_field) || include(joinpath(@__DIR__, "..", "p12", "result.jl"))

# THE ARTIFACT PATH IS OVERRIDABLE, AND THAT IS A SAFETY REQUIREMENT RATHER THAN A CONVENIENCE.
# Three PREVIOUS mini-spike artifacts are tracked in git and are cited BY NAME in
# `12-15-RERUN-VERDICT.md` (`p12_minispike_report.jld2`, `.rerun_control_18ep.jld2`,
# `.rerun_treatment_100ep.jld2`). A fixed const meant every re-run overwrote the unseeded baseline
# in place, so the only way to keep the earlier record was to remember to rename the file
# afterwards. Superseded numbers are recorded ADDITIVELY in this phase; an artifact that can only
# be superseded by being destroyed is the same rule broken on disk instead of in a document.
const MINISPIKE_REPORT_PATH =
    get(ENV, "P12_MINISPIKE_REPORT", joinpath(@__DIR__, "p12_minispike_report.jld2"))

"The three pre-registered arms, ablation LAST. Two arms would be a different experiment."
minispike_arms() = (:car, :gp, :none)

# The scoring block: the pool HEAD. Sized well above `P12_STAGE2_N_MIN = 271` so the coverage
# tolerance `select_prior` applies is used at the n it was derived for.
const MINISPIKE_N_HELDOUT = 1_000

# Posterior draws per scored dataset. `P12_SBC_L + 1` = 1000, the rank-bin-even count Tier 1 fixes.
const MINISPIKE_N_DRAWS = P12_SBC_L + 1

# Datasets scored per `sampleposterior` call. A CHUNK SIZE, not a threshold: 1000 datasets x 1000
# draws x 72 rows at once would be ~576 MB of posterior samples for no benefit.
const MINISPIKE_SCORE_CHUNK = 100

# Training epochs for each arm. Held IDENTICAL across arms -- the arm is the only thing that varies,
# which is what makes the comparison attributable (D-10).
const MINISPIKE_EPOCHS = parse(Int, get(ENV, "P12_MINISPIKE_EPOCHS", "18"))

# THE TRUNCATION FRACTIONS, NAMED HERE BEFORE THE CURVE IS COMPUTED (the plan requires the fraction
# be stated in the header before running).
#   `P12_K_PROD` = the smallest DCT rank whose TRUNCATION error falls to at most
#       MINISPIKE_KPROD_FRACTION of the winning arm's measured per-region POSTERIOR RMSE -- i.e. the
#       rank past which added modes are comfortably hidden by posterior noise.
#   `P12_N_LOW`  = the smallest rank whose truncation error falls BELOW the posterior noise at all
#       (fraction 1.0). This is the identified/vacuous boundary 12-18's SBC row classes need.
# `P12_N_LOW <= P12_K_PROD` HOLDS BY CONSTRUCTION, and that is the whole reason `n_low` was ruled a
# Tier-2 measurement rather than a Tier-1 guess: a STRICTER fraction requires MORE modes, so the
# stricter threshold's rank can never be the smaller of the two. The assertion at the append site
# CHECKS the construction rather than trusting it.
const MINISPIKE_KPROD_FRACTION = 0.5
const MINISPIKE_NLOW_FRACTION  = 1.0

# The pre-registered DCT ranks the curve is evaluated at.
const MINISPIKE_K_GRID = (3, 8, 15, 24, 35, 48, 63)

"""
    _p12ms_check_budget(elapsed_min, stage)

THE MINI-SPIKE WALL-CLOCK ENFORCER -- the first thing in this repository to read
`P12_MINISPIKE_WALLCLOCK_CEILING_MIN`. Throws with the literal `BLOCKER` on exceedance.

NOT the datagen ceiling. That one is enforced per call inside `generate_p12_pool`; this is THIS
RUN'S TOTAL across datagen, training and scoring for all three arms.
"""
function _p12ms_check_budget(elapsed_min::Real, stage::AbstractString)
    if elapsed_min > P12_MINISPIKE_WALLCLOCK_CEILING_MIN
        println("!"^78)
        println("PHASE-12 MINI-SPIKE WALL-CLOCK BLOCKER")
        println("!"^78)
        throw(ErrorException(
            "BLOCKER: the mini-spike total reached $(round(elapsed_min; digits = 2)) min at stage " *
            "'$stage', exceeding P12_MINISPIKE_WALLCLOCK_CEILING_MIN = " *
            "$(P12_MINISPIKE_WALLCLOCK_CEILING_MIN) min. THIS IS A SEPARATE BUDGET FROM " *
            "P12_DATAGEN_WALLCLOCK_CEILING_MIN, which happens to carry the same number and is " *
            "enforced per generate_p12_pool call. Record the shortfall as a BLOCKER naming the " *
            "measured rate. THE FOLLOWING ESCAPES ARE FORBIDDEN: shrinking P12_MINISPIKE_N, " *
            "dropping an arm, or downgrading the image size -- three arms at the pre-registered " *
            "size IS the design, and a two-arm run is a different experiment. Both ceilings are " *
            "APPEND-ONLY Tier 1 and may NOT be raised."))
    end
    return elapsed_min
end

"""
    select_prior(scores) -> (selected, admissible, beats_ablation, reason)

THE SELECTION RULE, FIXED IN CODE BEFORE ANY NUMBER EXISTS, so it cannot drift between the header's
prose and the run's behaviour. `scores` maps arm symbol -> NamedTuple with at least
`(coverage, rmse, n_test)`.

  1. ADMISSIBLE only if Gaussian-space per-region coverage is within `P12_STAGE2_COVERAGE_TOST_DELTA`
     of `P12_COVERAGE_NOMINAL`. A MISCALIBRATED ARM IS DISQUALIFIED, NOT "BETTER" -- the same rule
     the SC3 amendment pre-declares for the coverage gate. **`n_test >= P12_STAGE2_N_MIN` is
     ASSERTED before that bound is applied**, because the 0.03 tolerance was derived in 12-01 for
     N >= 271 DATASETS; applying it at an unasserted n is the unit/power mismatch this phase has
     already corrected three times.
  2. Among admissible arms, select the LOWER pooled per-region RMSE against the drawn lattice.
  3. The winner must BEAT the neutralized `:none` arm on that RMSE by any margin. If it does not,
     return `NONE-BEATS-ABLATION` -- a reportable outcome, never silently overwritten by picking a
     winner anyway.
  4. TIES (RMSE within 1 %) resolve to **CAR**, declared here IN ADVANCE because it is
     lattice-native and sparse and therefore the cheaper production choice. Recording the tie-break
     before the numbers exist is what stops it becoming a post-hoc preference.

A PURE FUNCTION OF ITS ARGUMENT: same input, same answer, no globals read but the frozen constants.
"""
function select_prior(scores)
    for a in (:car, :gp, :none)
        haskey(scores, a) || error("select_prior: missing arm :$a; three arms is the design")
        scores[a].n_test >= P12_STAGE2_N_MIN || error(
            "select_prior: arm :$a has n_test = $(scores[a].n_test) but " *
            "P12_STAGE2_COVERAGE_TOST_DELTA = $(P12_STAGE2_COVERAGE_TOST_DELTA) was derived for " *
            "N >= P12_STAGE2_N_MIN = $(P12_STAGE2_N_MIN) DATASETS. Applying a tolerance at an " *
            "unasserted n is the unit/power mismatch this phase has already corrected three times.")
    end

    adm = Symbol[]
    for a in (:car, :gp)
        abs(scores[a].coverage - P12_COVERAGE_NOMINAL) <= P12_STAGE2_COVERAGE_TOST_DELTA &&
            push!(adm, a)
    end

    isempty(adm) && return (selected = "NONE-BEATS-ABLATION", admissible = adm,
                            beats_ablation = false,
                            reason = "BOTH spatial arms are miscalibrated: coverage " *
                                     "$(round(scores[:car].coverage; digits=4)) (car) and " *
                                     "$(round(scores[:gp].coverage; digits=4)) (gp) against " *
                                     "$(P12_COVERAGE_NOMINAL) +/- $(P12_STAGE2_COVERAGE_TOST_DELTA). " *
                                     "A miscalibrated arm is DISQUALIFIED, not 'better' (rule 1).")

    # Rule 2, with rule 4's tie-break applied first so a tie cannot be decided by float noise.
    best = if length(adm) == 1
        only(adm)
    else
        rc, rg = scores[:car].rmse, scores[:gp].rmse
        abs(rc - rg) / min(rc, rg) <= 0.01 ? :car : (rc < rg ? :car : :gp)
    end

    if !(scores[best].rmse < scores[:none].rmse)
        return (selected = "NONE-BEATS-ABLATION", admissible = adm, beats_ablation = false,
                reason = "the best admissible arm :$best (RMSE $(round(scores[best].rmse; digits=5))) " *
                         "does NOT beat the neutralized ablation " *
                         "(RMSE $(round(scores[:none].rmse; digits=5))). The spatial premise is " *
                         "UNSUPPORTED at mini-spike scale (rule 3).")
    end

    tie = length(adm) == 2 &&
          abs(scores[:car].rmse - scores[:gp].rmse) / min(scores[:car].rmse, scores[:gp].rmse) <= 0.01
    return (selected = uppercase(string(best)), admissible = adm, beats_ablation = true,
            reason = tie ?
                "TIE (RMSE within 1 %) resolved to CAR by the tie-break declared in advance " *
                "(lattice-native and sparse, the cheaper production choice)." :
                "admissible = $adm; lower pooled per-region RMSE, and it beats the ablation.")
end

"""
    _p12ms_truncation_curve(Z) -> (ks, rmse)

For each rank `K` in `MINISPIKE_K_GRID`, the RMSE of reconstructing the TRUE field from only its
first `K` DCT modes in `p12_dct_order`.

A PROPERTY OF THE FIELD ENSEMBLE, NOT OF ANY NET. It answers "what rank is enough?" as ARITHMETIC,
so the production rank is not chosen by looking at which head trained best -- which would be
selection on the outcome.
"""
function _p12ms_truncation_curve(Ztrue::AbstractMatrix)
    order = p12_dct_order(P12_G)
    nreg  = size(Ztrue, 1)
    ks    = collect(MINISPIKE_K_GRID)
    out   = Vector{Float64}(undef, length(ks))
    coeffs = hcat([p12_dct_vec(vec(Ztrue[:, j])) for j in axes(Ztrue, 2)]...)   # nreg x n
    for (i, K) in enumerate(ks)
        keep = order[1:min(K + 1, length(order))]     # +1: rank K deviations PLUS the c0 term
        acc, cnt = 0.0, 0
        for j in axes(coeffs, 2)
            c = zeros(eltype(coeffs), nreg)
            c[keep] = coeffs[keep, j]
            r = p12_idct_vec(c) .- vec(Ztrue[:, j])
            acc += sum(abs2, r); cnt += length(r)
        end
        out[i] = sqrt(acc / cnt)
    end
    return ks, out
end

"`_p12ms_first_rank_below(ks, curve, bar)` -- the smallest rank whose truncation RMSE is <= `bar`."
function _p12ms_first_rank_below(ks, curve, bar)
    for (k, v) in zip(ks, curve)
        v <= bar && return k
    end
    return last(ks)      # the full-rank head is the honest answer when nothing clears the bar
end

"""
    _p12ms_score_arm(bundle, Zraw_held, Ztrue_held) -> NamedTuple

Score one trained arm on the HELD-OUT block. Reconstructs each posterior draw's field from its DCT
coefficients and compares to the drawn `z_field` in GAUSSIAN space (R-2, atom-free).

TWO INVERSES ARE REQUIRED AND THIS FUNCTION HAS GOT EACH OF THEM WRONG ONCE. They are independent
and neither substitutes for the other:

  * Pitfall 5 -- THE DRAWS ARE UN-STANDARDIZED FIRST. `sampleposterior` returns theta in the
    STANDARDIZED space the estimator was trained in; `bundle.theta_zt` is the frozen inverse.
    Scoring without it compares a standardized draw to raw truth and reports noise as recovery.
  * Pitfall 6 -- THE SMOOTHNESS PERMUTATION IS INVERTED BEFORE THE INVERSE DCT, through
    `p12_region_field`. theta row `k` is the `k`-th SMOOTHEST mode, not flat mode `k`; treating it
    as flat mode `k` reconstructs a scrambled field against which a PERFECT posterior still scores
    ~1.30x the trivial predictor.

Both were silent, both survived every structural check, and each invalidated a complete run. See
the header for the measured numbers and for exactly which reported quantities each one reached.
"""
function _p12ms_score_arm(bundle, Zraw_held::AbstractMatrix, Ztrue_held::AbstractMatrix)
    nreg  = P12_G^2
    ntest = size(Zraw_held, 2)

    # THE HEAD-WIDTH GUARD. Everything below reads rows `1:nreg` as "the field's G^2 DCT
    # coefficients", which is true ONLY at the full-rank head: `P12_THETA_ROWS(63)` puts `:c63` at
    # row 64 and `:spillover` at row 65. At a TRUNCATED production head (`K_dev = P12_K_PROD`, which
    # 12-17 trains) rows 1:64 would run PAST the deviation block and this function would reconstruct
    # a field out of `spillover` and `shift_dx` -- a wrong-space error one width away from the one
    # this file has already made twice, latent until exactly the run that matters most. It is
    # asserted rather than documented because no structural check can see it: the shapes all agree.
    bundle.K_dev == P12_K_DEV || error(
        "_p12ms_score_arm: this scorer reads theta rows 1:$(nreg) as the field's G^2 DCT " *
        "coefficients, which holds ONLY at the full-rank head K_dev = P12_K_DEV = $(P12_K_DEV). " *
        "The bundle carries K_dev = $(bundle.K_dev), so rows 1:$(nreg) run past the deviation " *
        "block into the nuisances and the reconstructed 'field' would be built from spillover and " *
        "shift_dx. Widen this function to read `1 + bundle.K_dev` rows and pad the truncated modes " *
        "with zeros before scoring a production head -- do NOT relax this check.")
    bundle.D == P12_D_MINISPIKE || error(
        "_p12ms_score_arm: expected the mini-spike layout D = P12_D_MINISPIKE = " *
        "$(P12_D_MINISPIKE), got D = $(bundle.D). `M[end, :]` is read as the r1 row and " *
        "`M[1:$(nreg), :]` as the field; neither holds at another width.")

    Zstd  = standardize_p12(Zraw_held, bundle.zt)          # 12-14's frozen transform, never refit

    # SEED THE SAMPLING. `sampleposterior` is STOCHASTIC, so seeding the weight init alone leaves
    # the artifact irreproducible: identical nets still give different scores. Counter + 3, off the
    # same reserved stream as the mask (+0/+1) and weight-init (+2) seeds, so still no new Tier-1
    # constant.
    #
    # DELIBERATELY THE SAME SEED FOR EVERY ARM -- COMMON RANDOM NUMBERS. The arms are scored on the
    # same held-out block, so giving them a common Monte-Carlo stream removes sampling noise as a
    # source of BETWEEN-ARM difference. That is variance reduction on the comparison, not a
    # thumb on it: it cannot favour any arm, and with 1000 datasets x 1000 draws the residual MC
    # noise is small either way. It matters because 12-15's decisive margin was 0.33 %.
    Random.seed!(rand(p12_rng(P12_MINISPIKE_COUNTER + 3), UInt64))

    se_sum   = zeros(Float64, nreg)      # per-region squared error of the posterior MEAN
    cov_hit  = zeros(Int, nreg)
    c0_se    = 0.0
    c0_hit   = 0
    r1_sds   = Float64[]
    a = (1 - P12_COVERAGE_NOMINAL) / 2

    for lo in 1:MINISPIKE_SCORE_CHUNK:ntest
        hi = min(lo + MINISPIKE_SCORE_CHUNK - 1, ntest)
        Zc = reshape_summary(Zstd[:, lo:hi], P12_G)
        smp = sampleposterior(bundle.estimator, Zc; N = MINISPIKE_N_DRAWS, use_gpu = false)
        raws = smp isa AbstractVector ? smp : [smp]
        # PITFALL 5: standardized theta -> RAW theta, through the bundle's own frozen transform.
        # This is the line whose absence invalidated the 2026-07-31T09:51Z run.
        mats = [StatsBase.reconstruct(bundle.theta_zt, Float64.(Mr)) for Mr in raws]
        for (t, M) in enumerate(mats)
            j = lo + t - 1
            # Rows 1:(1+K_dev) are [c0; the deviation coefficients] -- the field's DCT vector IN
            # SMOOTHNESS ORDER. PITFALL 6: the permutation must be inverted before the inverse DCT,
            # and `p12_region_field` is the ONE inverse that does it (`spike/p12/result.jl:562`).
            # Calling `p12_idct_vec` directly here -- which this file did until 2026-07-31 -- treats
            # theta row k as flat mode k and reconstructs a SCRAMBLED field.
            coeffs = @view M[1:nreg, :]
            fields = hcat([vec(p12_region_field(view(coeffs, :, d); G = P12_G))
                           for d in axes(coeffs, 2)]...)                              # nreg x N
            truth  = vec(Ztrue_held[:, j])
            mu     = vec(mean(fields; dims = 2))
            se_sum .+= abs2.(mu .- truth)
            for r in 1:nreg
                lo_r, hi_r = _p12_interval(view(fields, r, :), P12_COVERAGE_NOMINAL)
                (lo_r <= truth[r] <= hi_r) && (cov_hit[r] += 1)
            end
            # The GLOBAL term alone, so a reader can see how much of any arm difference is global
            # rather than spatial.
            c0_draws = vec(M[1, :])
            c0_true  = p12_dct_vec(truth)[1]
            c0_se   += abs2(mean(c0_draws) - c0_true)
            lo_c, hi_c = _p12_interval(c0_draws, P12_COVERAGE_NOMINAL)
            (lo_c <= c0_true <= hi_c) && (c0_hit += 1)
            push!(r1_sds, std(vec(M[end, :])))
        end
    end

    rmse_per_region = sqrt.(se_sum ./ ntest)
    cov_per_region  = cov_hit ./ ntest
    # REPORTING-ONLY (F3). NOT the same quantity as 12-13's ridge ratio: this is a posterior WIDTH,
    # that is a linear point predictor's residual sd. Named apart on purpose; read as a PAIR under
    # 12-13:131-152, and a disagreement is a FINDING, never averaged.
    r1_shrinkage = mean(r1_sds) / std(P12_R1_PRIOR)
    # REPORTED, GATES NOTHING. `trivial_rmse` is the RMSE of the null predictor "the field is zero"
    # -- i.e. the ensemble's own marginal sd -- and `skill` is the fraction of it the posterior mean
    # actually removes. No pre-registered bar reads either. They exist because an arm that has learned
    # NOTHING scores rmse ~= trivial_rmse, and without this line that fact is invisible in the report:
    # a pooled RMSE of 0.99 looks like a measurement rather than like the prior talking.
    trivial_rmse = sqrt(mean(abs2, Ztrue_held))
    return (rmse = sqrt(mean(se_sum) / ntest),
            trivial_rmse = trivial_rmse,
            skill = 1 - sqrt(mean(se_sum) / ntest) / trivial_rmse,
            rmse_per_region = rmse_per_region,
            coverage = mean(cov_per_region),
            coverage_per_region = cov_per_region,
            c0_rmse = sqrt(c0_se / ntest),
            c0_coverage = c0_hit / ntest,
            r1_shrinkage = r1_shrinkage,
            r1_vacuous = r1_shrinkage > P12_VACUOUS_SHRINKAGE_FLOOR,
            n_test = ntest)
end

function main()
    println("="^78)
    println("D-12 STAGE-2 MINI-SPIKE — three arms at matched r1, full-rank head (K = $(P12_K_DEV))")
    println("This run SELECTS the prior; every other quantity it reports gates nothing.")
    println("="^78)
    t0 = time()
    datagen_min = Dict{Symbol,Float64}()
    train_min   = Dict{Symbol,Float64}()
    score_min   = Dict{Symbol,Float64}()
    scores      = Dict{Symbol,Any}()
    bundles     = Dict{Symbol,Any}()
    pool_dirs   = Dict{Symbol,String}()

    held = 1:MINISPIKE_N_HELDOUT
    rest = (MINISPIKE_N_HELDOUT + 1):P12_MINISPIKE_N
    @assert length(held) >= P12_STAGE2_N_MIN
    println("held-out scoring block = $(first(held)):$(last(held))  ($(length(held)) datasets, " *
            "P12_STAGE2_N_MIN = $(P12_STAGE2_N_MIN))")
    println("training block         = $(first(rest)):$(last(rest))  ($(length(rest)) datasets)")
    println("budgets: datagen -> P12_DATAGEN_WALLCLOCK_CEILING_MIN = " *
            "$(P12_DATAGEN_WALLCLOCK_CEILING_MIN) (per pool call);  this run's TOTAL -> " *
            "P12_MINISPIKE_WALLCLOCK_CEILING_MIN = $(P12_MINISPIKE_WALLCLOCK_CEILING_MIN)")

    for arm in minispike_arms()
        println("\n" * "#"^78); println("ARM :$arm"); println("#"^78)

        td = time()
        # Pinned image size, NOT the F5 mixture: this is a CONTROLLED comparison between kernels,
        # and the image-size mixture belongs to the full run in 12-17. A mixture here would let an
        # arm difference be an image-size difference.
        dir = generate_p12_pool(P12_MINISPIKE_N; arm = arm,
                                imsize_sampler = _ -> P12_MINISPIKE_IMSIZE,
                                imsize_tag = Symbol("pinned_", P12_MINISPIKE_IMSIZE[1], "x",
                                                    P12_MINISPIKE_IMSIZE[2]))
        @assert p12_pool_complete(dir, P12_MINISPIKE_N)
        pool_dirs[arm]   = dir
        datagen_min[arm] = (time() - td) / 60
        _p12ms_check_budget((time() - t0) / 60, "datagen :$arm")

        tt = time()
        b = train_p12_npe(; arm = arm, n = P12_MINISPIKE_N, pool_dir = dir,
                          pool_indices = rest, epochs = MINISPIKE_EPOCHS, verbose = true)
        train_min[arm] = (time() - tt) / 60
        bundles[arm] = b
        _p12ms_check_budget((time() - t0) / 60, "training :$arm")

        # THE LEAK CHECK, asserted against the bundle's OWN recorded index sets rather than
        # re-derived arithmetic -- and against BOTH, because the val block drives early stopping.
        seen = vcat(b.train_indices, b.val_indices)
        @assert isempty(intersect(seen, collect(held))) """
            run_p12_minispike: arm :$arm trained on $(length(intersect(seen, collect(held)))) of the
            held-out scoring indices. Scoring here would report memorization as recovery."""
        @assert sort(seen) == collect(rest)

        ts = time()
        pool = load_p12_pool(dir)
        # THE COLUMN/GLOBAL-INDEX WITNESS. The leak argument above is stated in GLOBAL INDEX terms
        # ("generate_p12_sample keys its RNG on the global index alone") but ENFORCED in COLUMN
        # terms (`[:, held]`, `pool_indices = rest`). The two coincide because `generate_p12_pool`
        # always generates `1:n` and the shards concatenate in numeric order -- which is true, and
        # was unwitnessed until now. It becomes false the first time anyone gives the pool builder
        # an offset, and nothing else in this file would notice.
        @assert pool.global_index == collect(1:size(pool.summary_min, 2)) """
            run_p12_minispike: pool column j is not global index j in $dir. Every disjointness and
            held-out claim in this run is argued in global-index terms and applied by column."""
        Zraw_held  = Float64.(pool.summary_min)[:, held]
        Ztrue_held = reshape(pool.z_field, P12_G^2, size(pool.summary_min, 2))[:, held]
        scores[arm] = _p12ms_score_arm(b, Zraw_held, Ztrue_held)
        score_min[arm] = (time() - ts) / 60
        _p12ms_check_budget((time() - t0) / 60, "scoring :$arm")

        s = scores[arm]
        println("  pooled per-region RMSE = $(round(s.rmse; digits = 5))   " *
                "coverage = $(round(s.coverage; digits = 5)) " *
                "(nominal $(P12_COVERAGE_NOMINAL) +/- $(P12_STAGE2_COVERAGE_TOST_DELTA))")
        println("  trivial RMSE (predict 0) = $(round(s.trivial_rmse; digits = 5))   " *
                "skill = $(round(s.skill; digits = 5))   [REPORTED; gates nothing]")
        println("  c0 RMSE = $(round(s.c0_rmse; digits = 5))   " *
                "c0 coverage = $(round(s.c0_coverage; digits = 5))")
        println("  r1_shrinkage = $(round(s.r1_shrinkage; digits = 5))  " *
                "vacuous = $(s.r1_vacuous)   [REPORTED; a posterior WIDTH, NOT 12-13's ridge ratio]")
        println("  datagen $(round(datagen_min[arm]; digits=2)) min | " *
                "train $(round(train_min[arm]; digits=2)) min | " *
                "score $(round(score_min[arm]; digits=2)) min")
    end

    # --- THE SELECTION, through the one function ---------------------------------------------
    sel = select_prior(scores)

    # --- THE TRUNCATION CURVE: arithmetic on the field ensemble, not on any net ---------------
    #
    # TWO CURVES ARE COMPUTED AND BOTH ARE REPORTED, AND NEITHER IS PROMOTED. The curve has always
    # been taken from the `:car` pool while `winner_rmse` -- the BAR it is read against -- can come
    # from `:gp`, or, when no arm is selected, from whichever spatial arm scored lower. CAR, GP and
    # `:none` ensembles have genuinely different spectra (a `:none` field is white and needs every
    # mode), so a CAR-ensemble curve read against a non-CAR bar compares two different objects --
    # the Family-B shape this milestone has now met several times.
    #
    # THE FIX IS NOT MADE UNILATERALLY, BECAUSE IT WOULD CHANGE THE DEFINITION OF TWO TIER-2
    # CONSTANTS. `P12_K_PROD` and `P12_N_LOW` were authorised as "derive it from the measured
    # truncation curve" without a pool being named, so choosing one silently would settle a
    # definitional question by implementation. Both are therefore computed, both are written to the
    # artifact, and NEITHER IS APPENDED -- the choice is made with the two numbers side by side, by
    # the user, and only if a winner is ever actually selected.
    #
    # THE `:car` CURVE'S ROLE AS A CONTROL IS UNCHANGED AND MUST NOT BE LAUNDERED BY THIS CHANGE.
    # It was byte-identical across all three previous runs -- the one net-independent quantity that
    # did not move while every net-derived one reversed -- and that remains true and remains
    # evidence about the NET-side defect. It is NOT evidence that the curve was read against the
    # right bar. Both statements hold; neither licenses the other.
    winner_rmse = sel.beats_ablation ? scores[Symbol(lowercase(sel.selected))].rmse :
                                       minimum(scores[a].rmse for a in (:car, :gp))
    # The arm the BAR came from -- the selected one, or the lower-RMSE spatial arm on a no-selection
    # return. This is the pool the curve arguably ought to be taken from.
    bar_arm = sel.beats_ablation ? Symbol(lowercase(sel.selected)) :
              (scores[:car].rmse <= scores[:gp].rmse ? :car : :gp)

    curves    = Dict{Symbol,Vector{Float64}}()
    kprod_by  = Dict{Symbol,Int}()
    nlow_by   = Dict{Symbol,Int}()
    local ks
    for a in unique((:car, bar_arm))
        poolw = load_p12_pool(pool_dirs[a])
        Zt    = reshape(poolw.z_field, P12_G^2, size(poolw.summary_min, 2))[:, held]
        ks, cv = _p12ms_truncation_curve(Zt)
        curves[a]   = cv
        kprod_by[a] = _p12ms_first_rank_below(ks, cv, MINISPIKE_KPROD_FRACTION * winner_rmse)
        nlow_by[a]  = _p12ms_first_rank_below(ks, cv, MINISPIKE_NLOW_FRACTION  * winner_rmse)
        @assert nlow_by[a] <= kprod_by[a] "run_p12_minispike: n_low ($(nlow_by[a])) > k_prod " *
            "($(kprod_by[a])) on the :$a curve. Both come off ONE curve with n_low's fraction the " *
            "LOOSER, so this cannot happen unless the curve reader is wrong. Do NOT lower n_low " *
            "after the fact."
    end

    # The `:car` curve keeps the EXISTING key names, so the artifact stays directly comparable with
    # the three earlier runs and the control above stays readable. The bar-arm curve is ADDITIONAL.
    curve  = curves[:car]
    k_prod = kprod_by[:car]
    n_low  = nlow_by[:car]

    elapsed = _p12ms_check_budget((time() - t0) / 60, "final")

    println("\n" * "="^78)
    println("SELECTED PRIOR: ", sel.selected)
    println("PRODUCTION RANK: K = ", k_prod)
    println("="^78)
    println("  reason        : ", sel.reason)
    println("  admissible    : ", sel.admissible, "   beats_ablation = ", sel.beats_ablation)
    println("  truncation    : ", collect(zip(ks, round.(curve; digits = 5))), "   [:car pool]")
    println("  n_low         : ", n_low, "   (fraction $(MINISPIKE_NLOW_FRACTION) of winner RMSE)")
    println("  k_prod        : ", k_prod, "   (fraction $(MINISPIKE_KPROD_FRACTION))")
    println("  bar arm       : :", bar_arm, "   (the arm winner_rmse came from)")
    if bar_arm !== :car
        println("  truncation    : ", collect(zip(ks, round.(curves[bar_arm]; digits = 5))),
                "   [:$(bar_arm) pool -- the BAR's own ensemble]")
        println("  n_low/k_prod on the :$(bar_arm) curve : ", nlow_by[bar_arm], " / ",
                kprod_by[bar_arm])
    end
    println("  BOTH CURVES ARE REPORTED AND NEITHER IS APPENDED -- the pool the curve is taken " *
            "from is a definitional choice for the user, not for this run.")
    println("  TOTAL         : $(round(elapsed; digits = 2)) min against " *
            "P12_MINISPIKE_WALLCLOCK_CEILING_MIN = $(P12_MINISPIKE_WALLCLOCK_CEILING_MIN)")
    println("  (datagen was judged PER POOL against P12_DATAGEN_WALLCLOCK_CEILING_MIN = " *
            "$(P12_DATAGEN_WALLCLOCK_CEILING_MIN); the two budgets are NOT the same budget)")

    tmp = MINISPIKE_REPORT_PATH * ".tmp"
    jldsave(tmp;
        schema_version = 1,
        generated = string(Dates.now(Dates.UTC)) * "Z",
        selected_prior = sel.selected,
        admissible = collect(sel.admissible),
        beats_ablation = sel.beats_ablation,
        selection_reason = sel.reason,
        selection_rule = "1. ADMISSIBLE only if |coverage - P12_COVERAGE_NOMINAL| <= " *
                         "P12_STAGE2_COVERAGE_TOST_DELTA, asserted at n_test >= P12_STAGE2_N_MIN " *
                         "(a miscalibrated arm is DISQUALIFIED, not 'better'). 2. Among admissible, " *
                         "the LOWER pooled per-region RMSE. 3. The winner must BEAT the neutralized " *
                         ":none arm by any margin, else NONE-BEATS-ABLATION. 4. Ties within 1 % " *
                         "resolve to CAR, declared in advance. Fixed in code before any number.",
        arms = collect(minispike_arms()),
        rmse = Dict(string(a) => scores[a].rmse for a in keys(scores)),
        trivial_rmse = Dict(string(a) => scores[a].trivial_rmse for a in keys(scores)),
        skill = Dict(string(a) => scores[a].skill for a in keys(scores)),
        rmse_per_region = Dict(string(a) => scores[a].rmse_per_region for a in keys(scores)),
        coverage = Dict(string(a) => scores[a].coverage for a in keys(scores)),
        coverage_per_region = Dict(string(a) => scores[a].coverage_per_region for a in keys(scores)),
        c0_rmse = Dict(string(a) => scores[a].c0_rmse for a in keys(scores)),
        c0_coverage = Dict(string(a) => scores[a].c0_coverage for a in keys(scores)),
        r1_shrinkage = Dict(string(a) => scores[a].r1_shrinkage for a in keys(scores)),
        r1_vacuous = Dict(string(a) => scores[a].r1_vacuous for a in keys(scores)),
        truncation_curve = Dict("k" => ks, "rmse" => curve),
        k_prod_recommended = k_prod,
        n_low_recommended = n_low,
        # --- H2: the SAME quantities on the BAR'S OWN ensemble. ADDED keys, nothing renamed, so
        #     this artifact stays comparable with the three earlier runs. Identical to the `:car`
        #     entries whenever the bar came from the CAR arm.
        truncation_curve_bar_arm = Dict("k" => ks, "rmse" => curves[bar_arm]),
        k_prod_bar_arm = kprod_by[bar_arm],
        n_low_bar_arm = nlow_by[bar_arm],
        truncation_bar_arm = string(bar_arm),
        truncation_pool_note =
            "TWO CURVES, NEITHER PROMOTED. `truncation_curve` / `k_prod_recommended` / " *
            "`n_low_recommended` are taken from the :car pool, which is what all previous runs " *
            "reported and is why they remain directly comparable. `truncation_curve_bar_arm` / " *
            "`k_prod_bar_arm` / `n_low_bar_arm` are the same arithmetic on the ensemble the BAR " *
            "(`winner_rmse_used_for_fractions`) actually came from, named in " *
            "`truncation_bar_arm`. CAR, GP and :none ensembles have different spectra, so a " *
            "CAR-ensemble curve read against a non-CAR bar compares two different objects. " *
            "P12_K_PROD and P12_N_LOW were authorised as 'derive from the measured truncation " *
            "curve' WITHOUT a pool being named, so which of these two pairs is the constant is a " *
            "DEFINITIONAL choice for the user. NEITHER PAIR IS APPENDED BY THIS RUN. Separately: " *
            "the :car curve being byte-identical across runs is evidence about the NET-side " *
            "defect, not evidence that it was read against the right bar.",
        kprod_fraction = MINISPIKE_KPROD_FRACTION,
        nlow_fraction = MINISPIKE_NLOW_FRACTION,
        winner_rmse_used_for_fractions = winner_rmse,
        n_test = length(held),
        heldout_indices = collect(held),
        train_indices_by_arm = Dict(string(a) => bundles[a].train_indices for a in keys(bundles)),
        val_indices_by_arm = Dict(string(a) => bundles[a].val_indices for a in keys(bundles)),
        pool_dirs = Dict(string(a) => pool_dirs[a] for a in keys(pool_dirs)),
        epochs = MINISPIKE_EPOCHS,
        imsize = P12_MINISPIKE_IMSIZE,
        n_pool = P12_MINISPIKE_N,
        n_draws = MINISPIKE_N_DRAWS,
        datagen_min = Dict(string(a) => datagen_min[a] for a in keys(datagen_min)),
        train_min = Dict(string(a) => train_min[a] for a in keys(train_min)),
        score_min = Dict(string(a) => score_min[a] for a in keys(score_min)),
        elapsed_min = elapsed,
        minispike_ceiling_min = P12_MINISPIKE_WALLCLOCK_CEILING_MIN,
        datagen_ceiling_min = P12_DATAGEN_WALLCLOCK_CEILING_MIN,
        budget_note = "TWO SEPARATE BUDGETS THAT HAPPEN TO CARRY THE SAME NUMBER. " *
                      "elapsed_min is THIS RUN'S TOTAL against " *
                      "P12_MINISPIKE_WALLCLOCK_CEILING_MIN. datagen_min is per arm and each pool " *
                      "was judged INDEPENDENTLY inside generate_p12_pool against " *
                      "P12_DATAGEN_WALLCLOCK_CEILING_MIN; their SUM is not a quantity anything " *
                      "compares to that ceiling. Do not add them together or report one against " *
                      "the other's constant.",
        julia_version = string(VERSION),
        caption = "Phase-12 D-12 Stage-2 mini-spike: three arms (:car, :gp, :none) at MATCHED " *
                  "induced r1 through a full-rank K = $(P12_K_DEV) head, scored in Gaussian DCT " *
                  "space on a held-out block CARVED FROM THE POOL HEAD (never a second pool -- " *
                  "CONVENTIONS C-04). This run SELECTS the prior under the rule in " *
                  "`selection_rule`; every other quantity here gates nothing.",
    )
    let d = JLD2.load(tmp)
        for k in ("selected_prior", "truncation_curve", "k_prod_recommended")
            @assert haskey(d, k) "artifact integrity: $k missing"
        end
    end
    mv(tmp, MINISPIKE_REPORT_PATH; force = true)
    println("\nwrote ", MINISPIKE_REPORT_PATH)
    return nothing
end

isdefined(@__MODULE__, :P12_MINISPIKE_LOAD_ONLY) || main()
