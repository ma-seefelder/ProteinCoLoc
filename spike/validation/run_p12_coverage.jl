# spike/validation/run_p12_coverage.jl --- the D-09 leave-region-out coverage and scoring run.
#
# THIS RUN IS THE **SIMULATION** ARM ONLY. The `--real` switch is wired below and is left
# UNEXERCISED here. Scoring the two physical specimens before the construction has been validated
# on simulation would spend the only check that exists: on a real image there is no ground truth,
# so a wrong construction cannot be detected there at all. `--real` is 12-19's, and 12-19 does not
# run on this route.
#
# =============================================================================================
# WHAT THIS RUN CAN AND CANNOT CONCLUDE -- FIXED BEFORE IT LAUNCHED
# =============================================================================================
# The reading of this run's numbers was pre-declared, in its own commit, BEFORE any number
# existed: .planning/phases/12-spatial-colocalization-map/12-16-PREDECLARATION.md (commit df9ef91,
# 2026-07-31T22:56:50). The branch this run lands in is computed HERE, by
# `_p12cov_predeclared_branch`, so the pre-declaration is EXECUTABLE rather than prose that a
# reader has to apply by hand afterwards.
#
#   COVERAGE ALONE CANNOT SEPARATE *CALIBRATED AND INFORMATIVE* FROM *CALIBRATED BECAUSE
#   UNINFORMATIVE*. A predictive interval that is wide and centred near the marginal mean covers
#   ~90 % while carrying no information about which region is which. The DISCRIMINATOR is
#   `prior_only_floor`, and the margin is the EXISTING pre-registered P12_STAGE2_LOGSCORE_MIN
#   = 0.02 nats per held-out region -- reused, not invented, and not a new bar.
#
#   1. coverage in [0.87, 0.93] AND the floor beaten by >= 0.02  -> SPAT-05/06 supported.
#   2. coverage in band, floor NOT cleared -> "the uncertainty quantification is honest; the map
#      carries no information beyond the prior." NOT A PASS, and NOT to be written up as
#      "promising".
#   3. coverage outside the band -> SPAT-06 NOT MET. Reported, NOT TUNED.
#
# NO EPOCH, THRESHOLD OR N TUNING IF THE RESULT IS UNFAVOURABLE. P12_ITERATION_ALLOWANCE is
# already SPENT (12-D13-AUTHORISATION.md 6). ONE run at N = P12_STAGE2_N_MIN = 271.
#
# =============================================================================================
# THE INDEPENDENT UNIT IS THE **DATASET**. 271, NOT 17,344.
# =============================================================================================
# At 64 regions per image this run produces 271 * 64 = 17,344 region-draws but only 271
# INDEPENDENT units: the 64 regions of one image share every nuisance and the global term, so they
# are PSEUDO-REPLICATES (12-19-PLAN.md:148-150, 12-02-PLAN.md:235). An interval built on 17,344
# would be roughly 8x too tight and would manufacture significance out of pseudo-replication --
# which is the error this milestone has ALREADY MADE ONCE, applying a gate sized for N >= 271 to
# an arm whose effective_independent_n is 2. `effective_independent_n` is a mandatory artifact key
# beside `n_region_draws`, and the 17,344 figure never appears without the sentence saying the
# region-draws are pseudo-replicates.
#
# =============================================================================================
# THE TWO-SIDED GATE, AND WHY IT CANNOT FIRE HERE
# =============================================================================================
# `lro_pass` requires BOTH arms calibrated (|coverage - P12_COVERAGE_NOMINAL| <=
# P12_STAGE2_COVERAGE_TOST_DELTA) and then a spatial log-score advantage of at least
# P12_STAGE2_LOGSCORE_MIN. A coverage-only win is NOT a pass; a miscalibrated arm is
# DISQUALIFIED, not "worse".
#
# ON THIS ROUTE THERE IS NO SPATIAL ARM. The mini-spike returned NONE-BEATS-ABLATION and 12-17
# never ran, so `stage2_gate = :not_applicable_descope` BY CONSTRUCTION and the gate is
# UNEXERCISED. A gate reporting "fail" where there was no comparison would be read as evidence
# against the spatial prior, which would be wrong.
#
# =============================================================================================
# THE THREE ARMS, AND PITFALL 6's WORDING RULE
# =============================================================================================
#   ablation         -- the D-13 `arm = :none` net. Describe it ONLY as
#                       "no spatial borrowing, full nuisance and global borrowing". D-10's phrase
#                       "no mechanism to predict a held-out region beyond the prior" is FALSE as
#                       written: with r1 -> 0 the REGIONS become conditionally independent, but the
#                       other 63 observed regions still inform the shared nuisances (chromatic_eps
#                       above all, whose radial gradient is global) and the global term c0.
#   prior_only_floor -- the genuinely prior-only third arm, so a reader can SEE how much of the
#                       ablation's performance is nuisance/global borrowing rather than spatial.
#   spatial          -- ABSENT BY DESCOPE, not missing by error.
#
# =============================================================================================
# BUNDLE RESOLUTION, AND THE HASH THAT REFUSES
# =============================================================================================
#   p12_train_full_report.jld2 present -> score its two arms            (:full_two_arm)
#   else p12_ablation_report.jld2      -> score the single arm it names (:descope_ablation_only)
#   else                               -> THROW, naming both paths. A third, unhandled case does
#                                         not fall through silently.
# Each bundle's sha256 is re-verified against the report that actually exists ON THE ROUTE TAKEN
# (12-17 never ran on the descope route, so its report is absent and must not be consulted). A
# hash mismatch REFUSES TO SCORE -- that is the entire reason the hashes were recorded.
#
# =============================================================================================
# UN-STANDARDIZE THETA BEFORE THE PREDICTIVE INTERVAL IS BUILT
# =============================================================================================
# `sampleposterior` returns theta in the STANDARDIZED space the estimator was trained in, so every
# draw goes through `StatsBase.reconstruct(bundle.theta_zt, ...)` before it is compared to an
# observed correlation entry. The contract is `spike/npe/infer.jl:35-38`. 12-15's first run omitted
# it and produced a report that was COMPLETE, internally consistent, reconciled to 0.006 min
# against its own parts, and SCIENTIFICALLY VOID -- because no structural check looks at SCALE.
# (`infer.jl` calls this "Pitfall 5"; `12-RESEARCH.md:867` numbers a DIFFERENT Pitfall 5. Cite the
# file, never the number.) That un-standardization happens inside `lro_arm` and
# `_p12cov_single_stack` in `p12_coverage.jl`; this runner never reads a raw theta row itself.
#
# No `trivial_rmse` requirement attaches to this plan: it reports coverage and a log score, not an
# RMSE.
#
# =============================================================================================
# SPAT-07 AND SPAT-08 ARE DEFERRED TO v2.1 -- A DELIBERATE DECISION OF 2026-07-31
# =============================================================================================
# User ruling, recorded at 12-D13-AUTHORISATION.md 7 and 9. NOT DELIVERED, stated plainly:
#   * No Gaussian-space SBC, no randomized-rank rho-space SBC, no nuisance-appropriate equivalence
#     testing. `spat07_scope = :deferred_to_v2_1`, NOT `:reduced_descope`.
#   * THE PHASE MAKES NO PER-REGION CALIBRATION CLAIM BEYOND PREDICTIVE COVERAGE. What SPAT-06
#     establishes is that the leave-region-out predictive intervals cover what was actually
#     measured. It does NOT establish that the posterior is calibrated in the SBC sense, and the
#     two must not be conflated.
#   * SPAT-08's radial-energy, offset-grid and eps = 0 guards do not run, so the phase makes NO
#     radial-confound claim and the S-4 warning is left standing and unaddressed.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. Reaches `src/` only transitively and
# read-only. Adds no package: `SHA` is a stdlib on the default `@stdlib` LOAD_PATH entry, reached
# the same way `run_p12_ablation.jl:75` reaches it.
#
# Run:
#     julia --project=spike -t auto spike/validation/run_p12_coverage.jl

using JLD2
using Dates
using Statistics
using Random
using SHA

isdefined(@__MODULE__, :lro_scores) || include(joinpath(@__DIR__, "p12_coverage.jl"))

const P12COV_REPORT_PATH = joinpath(@__DIR__, "p12_coverage_sim_report.jld2")

# --- Run sizing, stated with its reason. NONE OF THESE IS A THRESHOLD --------------------------

# THE PRE-REGISTERED MINIMUM, NOT A CHOICE. `P12_STAGE2_N_MIN` = 271 is the N the coverage
# tolerance was DERIVED for (p12_consts.jl:409-422). It is not exceeded without a stated reason
# and it is never shortened to save time.
const P12COV_N_DATASETS = P12_STAGE2_N_MIN

# Posterior draws per held-out region. `P12_SBC_L + 1` = 1000 -- the rank-bin-even count Tier 1
# already fixes, reused rather than invented so no new number enters the run.
const P12COV_N_DRAWS = P12_SBC_L + 1

# The noise-model calibration: `P12COV_NEFF_THETA` independent parameters, each observed
# `P12COV_NEFF_OBS` times, at EVERY size in `P12_IMSIZE_SET`. 5 x 20 x 64 regions = 6,400
# per-region variance estimates per size, so the pooled Var(z) is precise to ~2 % while the
# between-theta spread stays measurable at 5 blocks.
const P12COV_NEFF_THETA = 5
const P12COV_NEFF_OBS   = 20

# The evaluation arm. THE F5 BINDING INVARIANT IS WHY THIS IS `:none` AND NOT `:car`: rank and
# coverage claims hold only under the joint the estimator was TRAINED on, and the D-13 bundle was
# trained on `arm = :none` pools. Scoring it against spatially correlated truth would measure
# MISSPECIFICATION, not calibration -- a different and unasked question. This is a NAMED LIMIT of
# the result, recorded in the artifact as `eval_arm` and in the caption, not left to be inferred.
const P12COV_EVAL_ARM = :none

# =============================================================================================
# Bundle resolution
# =============================================================================================

"""
    _p12cov_resolve_bundles(; dir) -> NamedTuple

Resolve which bundles this run scores, by FILE EXISTENCE, and verify each `sha256` against the
report that exists on the route actually taken. Refuses to score a bundle whose hash does not
match. Throws naming BOTH paths when neither report exists -- a third case does not fall through.
"""
function _p12cov_resolve_bundles(; dir = @__DIR__)
    full = joinpath(dir, "p12_train_full_report.jld2")
    abl  = joinpath(dir, "p12_ablation_report.jld2")

    if isfile(full)
        d = JLD2.load(full)
        arms = NamedTuple[]
        for k in ("car", "gp", "none", "spatial", "ablation")
            pk, hk = "bundle_path_$k", "bundle_sha256_$k"
            haskey(d, pk) && haskey(d, hk) &&
                push!(arms, (name = Symbol(k), path = d[pk], sha = d[hk]))
        end
        isempty(arms) && error(
            "run_p12_coverage: $full exists but names no per-arm bundle_path/bundle_sha256 " *
            "pair. Refusing to guess which net 12-17 trained.")
        return (source = :full_two_arm, report = full, arms = arms)
    elseif isfile(abl)
        d = JLD2.load(abl)
        for k in ("bundle_path", "bundle_sha256", "arm")
            haskey(d, k) || error("run_p12_coverage: $abl is missing the key `$k`.")
        end
        return (source = :descope_ablation_only, report = abl,
                arms = [(name = Symbol(d["arm"]), path = d["bundle_path"],
                         sha = d["bundle_sha256"])])
    end
    error("""
        run_p12_coverage: NEITHER bundle report exists. Both paths were checked:
            $full   (12-17, the PROCEED path)
            $abl    (12-14 / the D-13 descope deliverable)
        A third, unhandled case is not permitted to fall through silently. Produce one of them
        before scoring; this runner refuses to invent a bundle to score.""")
end

"""
    _p12cov_verify_sha(path, expected) -> String

Re-hash the bundle and REFUSE to score on a mismatch. Recording a hash and then not checking it
would make the record decorative; a mismatch means the net on disk is not the net the report
attributes the numbers to, and the comparison would be unattributable.
"""
function _p12cov_verify_sha(path, expected)
    isfile(path) || error("""
        run_p12_coverage: the report names a bundle that is not on disk:
            $path
        Trained bundles under `spike/artifacts/` are GITIGNORED BY DESIGN (durability is the
        P12_PRIMARY_CHECKOUT rule), so an absent bundle usually means this checkout is not the one
        that trained it. Refusing to score a bundle that does not exist.""")
    got = bytes2hex(open(sha256, path))
    got == expected || error("""
        run_p12_coverage: BUNDLE HASH MISMATCH -- refusing to score.
            path     : $path
            expected : $expected
            got      : $got
        The recorded hash exists precisely so that scoring a different net than the one the report
        attributes the numbers to is impossible. Do not relax this check.""")
    return got
end

# =============================================================================================
# The pre-declared reading, made executable
# =============================================================================================

"""
    _p12cov_predeclared_branch(coverage, logscore, floor_logscore) -> (branch, text)

Apply 12-16-PREDECLARATION.md section 3, which was committed BEFORE this run existed. The branch
is COMPUTED here rather than read off the numbers afterwards, because a reading fixed after the
numbers exist is not a pre-declaration.

`:calibrated_but_uninformative` carries its wording VERBATIM, because it is the most confusing
outcome available and the one most open to being softened. **It is NOT a pass and may not be
written up as "promising".**
"""
function _p12cov_predeclared_branch(coverage::Real, logscore::Real, floor_logscore::Real;
                                    nominal::Real = P12_COVERAGE_NOMINAL,
                                    tost_delta::Real = P12_STAGE2_COVERAGE_TOST_DELTA,
                                    logscore_min::Real = P12_STAGE2_LOGSCORE_MIN)
    in_band = abs(coverage - nominal) <= tost_delta
    margin  = logscore - floor_logscore
    if !in_band
        return (:miscalibrated,
            "SPAT-06 IS NOT MET: pooled leave-region-out coverage $(coverage) lies outside " *
            "[$(nominal - tost_delta), $(nominal + tost_delta)]. REPORTED, NOT TUNED -- moving " *
            "epochs, N, thresholds or the noise model to bring it into band would be tuning a " *
            "model until it passes its own calibration gate, and a gate that has been tuned to " *
            "has stopped measuring anything.")
    elseif margin >= logscore_min
        return (:supports_spat05_spat06,
            "SPAT-05/SPAT-06 SUPPORTED: coverage $(coverage) is in band AND the mean log score " *
            "beats prior_only_floor by $(margin) >= $(logscore_min) nats per held-out region. " *
            "The per-region map is calibrated AND informative.")
    else
        return (:calibrated_but_uninformative,
            "The uncertainty quantification is honest; the map carries no information beyond " *
            "the prior. THIS IS NOT A PASS, AND IT MAY NOT BE WRITTEN UP AS \"PROMISING\". " *
            "Coverage $(coverage) is in band, but the mean log score beats prior_only_floor by " *
            "only $(margin) < $(logscore_min) nats per held-out region. It is the direct " *
            "analogue of the 12-15 re-run's 2.3(a) -- the method works and the spatial prior " *
            "adds nothing -- one layer out.")
    end
end

# =============================================================================================
# Scoring accumulators
# =============================================================================================

"An arm's running tallies. `by_imsize` keeps the per-size breakdown 12-20's Guard 3 reads."
function _p12cov_new_arm(name::Symbol, G::Integer)
    return (name = name,
            inside = zeros(Int, G^2), n = zeros(Int, G^2),
            logscore = zeros(Float64, G^2), crps = zeros(Float64, G^2),
            width = zeros(Float64, G^2),
            per_dataset_cov = Float64[],
            ls_by_imsize = Dict{String,Vector{Float64}}(),
            cov_by_imsize = Dict{String,Vector{Float64}}())
end

function _p12cov_pool(a, G::Integer)
    used = findall(>(0), a.n)
    isempty(used) && error(
        "run_p12_coverage: arm :$(a.name) scored ZERO regions. Every region was degenerate or " *
        "masked, which is a defect in the run, not a coverage of 0.")
    tot   = sum(a.n[used])
    covp  = [a.inside[r] / a.n[r] for r in used]
    return (name = a.name,
            coverage = sum(a.inside[used]) / tot,
            coverage_per_region = [a.n[r] > 0 ? a.inside[r] / a.n[r] : NaN for r in 1:G^2],
            mean_logscore = sum(a.logscore[used]) / tot,
            logscore_per_region = [a.n[r] > 0 ? a.logscore[r] / a.n[r] : NaN for r in 1:G^2],
            mean_crps = sum(a.crps[used]) / tot,
            mean_width = sum(a.width[used]) / tot,
            per_dataset_coverage = a.per_dataset_cov,
            n_scored = tot,
            n_regions_used = length(used),
            coverage_spread_across_regions = length(covp) >= 2 ? std(covp) : NaN,
            logscore_by_imsize = Dict(k => mean(v) for (k, v) in a.ls_by_imsize),
            coverage_by_imsize = Dict(k => mean(v) for (k, v) in a.cov_by_imsize))
end

# =============================================================================================
# main
# =============================================================================================

"""
    main(args = ARGS; n_datasets, neff_theta, neff_obs, report_path) -> Dict

The reported run. **EVERY KEYWORD DEFAULTS TO THE PRE-REGISTERED VALUE**, so the command line in
this file's header reproduces the artifact exactly and no default can be reached by omission.

THE KEYWORDS EXIST FOR ONE PURPOSE: a SMOKE run that exercises the whole path in seconds before
~20 minutes of compute is committed to it. A run at any non-default sizing records
`smoke = true` in its payload and is written wherever `report_path` says -- it is NOT the
reported artifact and must never be mistaken for one. **They are not a tuning knob.** Shortening
`n_datasets` below `P12_STAGE2_N_MIN` for the reported run would apply a tolerance at an N it was
not derived for, which is the unit/power mismatch this phase has already corrected three times.
"""
function main(args = ARGS;
              n_datasets::Integer = P12COV_N_DATASETS,
              neff_theta::Integer = P12COV_NEFF_THETA,
              neff_obs::Integer   = P12COV_NEFF_OBS,
              report_path = P12COV_REPORT_PATH)
    smoke = (n_datasets != P12COV_N_DATASETS || neff_theta != P12COV_NEFF_THETA ||
             neff_obs != P12COV_NEFF_OBS || report_path != P12COV_REPORT_PATH)
    smoke && @warn "run_p12_coverage: SMOKE SIZING — this is NOT the reported artifact." *
                   " n_datasets=$n_datasets (pre-registered $(P12COV_N_DATASETS))," *
                   " neff=$(neff_theta)x$(neff_obs), path=$report_path"
    real_mode = "--real" in args
    real_mode && error("""
        run_p12_coverage: `--real` is 12-19's path and MUST NOT be run from plan 12-16.
        Scoring the two physical specimens before the D-09 construction has been validated on
        simulation would spend the only check that exists -- on a real image there is no ground
        truth, so a wrong construction cannot be detected there at all. 12-19 does not run on the
        D-13 descope route. The switch is wired and deliberately left UNEXERCISED.""")

    t0 = time()
    println("="^78)
    println("D-09 LEAVE-REGION-OUT COVERAGE AND SCORING — SIMULATION ARM")
    println("The reading of these numbers was PRE-DECLARED before this run launched:")
    println("  .planning/phases/12-spatial-colocalization-map/12-16-PREDECLARATION.md (df9ef91)")
    println("="^78)

    G    = P12_G
    nreg = G^2

    # --- 1. Resolve the bundles and verify each hash ------------------------------------------
    res = _p12cov_resolve_bundles()
    println("bundle_source = :$(res.source)   (from $(basename(res.report)))")
    bundles = Dict{Symbol,Any}()
    bundle_meta = Dict{String,Any}()
    for a in res.arms
        sha = _p12cov_verify_sha(a.path, a.sha)
        println("  arm :$(a.name)  sha256 VERIFIED  $(sha[1:16])…")
        bundles[a.name] = load_p12_npe(a.path)
        bundle_meta[string(a.name)] = Dict("path" => string(a.path), "sha256" => sha)
    end
    res.source === :descope_ablation_only &&
        println("  SPATIAL ARM: ABSENT BY DESCOPE, not missing by error. stage2_gate will be " *
                ":not_applicable_descope.")
    abl_name = res.source === :descope_ablation_only ? res.arms[1].name :
               (haskey(bundles, :none) ? :none : res.arms[end].name)
    bundle = bundles[abl_name]

    # --- 2. Calibrate the observation-noise model ---------------------------------------------
    println("\n[1/4] calibrating n_eff on the reserved coverage stream " *
            "($(neff_theta) theta x $(neff_obs) obs per size) ...")
    cal = calibrate_neff_all(; n_theta = neff_theta, n_obs = neff_obs,
                             rng = p12_rng(P12_COVERAGE_COUNTER), arm = P12COV_EVAL_ARM, G = G)
    for sz in P12_IMSIZE_SET
        k = string(sz)
        println("      $(rpad(k, 16)) n_eff = $(round(cal.by_imsize[k]; digits = 2))   " *
                "Var(z) = $(round(cal.varz_by_imsize[k]; digits = 6))   " *
                "between-theta sd = $(round(cal.varz_theta_spread_by_imsize[k]; digits = 6))")
    end
    n_eff = cal.n_eff
    println("      POOLED n_eff = $(round(n_eff; digits = 3))   " *
            "implied Var(z) = $(round(1 / (n_eff - 3); digits = 6))   " *
            "varz_fit_residual = $(round(cal.varz_fit_residual; digits = 6))  (REPORTED, gates nothing)")

    # THE PREDICTIVE USES THE **PER-SIZE** CONSTANT, NOT THE POOLED ONE. The image size is an
    # OBSERVABLE at read time, not a latent, so conditioning on it is available for free -- and
    # `P12_IMSIZE_SET` spans a 16x pixel-count range, over which n_eff moves by an order of
    # magnitude. Applying one pooled constant would make every 512-squared interval too narrow and
    # every 2048-squared interval too wide, and the pooled coverage would then be an average of two
    # opposite miscalibrations that could sit at nominal while neither size did. `12-16-PLAN.md:168`
    # specifies the per-size calibration for exactly this reason (Pitfall 2 measured the per-region
    # noise sd as image-size dependent), and `varz_fit_residual` is what makes the size dependence
    # visible rather than assumed. The POOLED value is still appended and reported, because it is
    # the single-number summary a reader will ask for -- it is simply not what the run applies.
    println("      the sweep applies the PER-SIZE n_eff (the image size is observed, not latent);" *
            " the pooled value is reported but NOT applied.")

    # --- 3. The LRO sweep ----------------------------------------------------------------------
    println("\n[2/4] leave-region-out sweep: $(n_datasets) datasets x $nreg regions " *
            "(arm = :$(P12COV_EVAL_ARM), N = $(P12COV_N_DRAWS) draws/region) ...")
    rng      = p12_rng(P12_COVERAGE_COUNTER + 1)   # a distinct sub-stream from the calibration
    rng_pred = p12_rng(P12_COVERAGE_COUNTER + 2)   # the observation-noise draws
    rng_flr  = p12_rng(P12_COVERAGE_COUNTER + 3)   # the prior-only floor

    arm_abl = _p12cov_new_arm(:ablation, G)
    arm_flr = _p12cov_new_arm(:prior_only_floor, G)
    param_inside = zeros(Int, nreg); param_n = zeros(Int, nreg)
    clamped_region_count = 0
    masked_region_count  = 0
    imsize_counts = Dict{String,Int}()

    for i in 1:n_datasets
        s   = draw_simulate_infer_p12(rng; arm = P12COV_EVAL_ARM, G = G)   # size from the F5 mixture
        key = string(s.imsize)
        imsize_counts[key] = get(imsize_counts, key, 0) + 1
        # The noise model conditioned on the OBSERVED image size -- see the note at [1/4].
        n_eff_i = cal.by_imsize[key]

        a   = lro_arm(bundle, s.Zraw, 1:nreg; N = P12COV_N_DRAWS, G = G)
        flr = prior_only_floor(; N = P12COV_N_DRAWS, rng = rng_flr, G = G)
        ztrue = vec(s.draw.z_field)

        # PER-DATASET hit counters. These are THIS dataset's own tallies, kept separate from the
        # arms' cumulative ones, because the per-dataset coverage FRACTION is the independent unit
        # this run's interval is built from. Deriving it from the cumulative tallies is exactly
        # how it would quietly become a region-draw statistic instead.
        hit_abl = 0; hit_flr = 0; nscored = 0
        for r in 1:nreg
            obs = s.Zraw[r]
            if s.Zraw[nreg + r] == 0.0                      # the patch was never scorable
                masked_region_count += 1
                continue
            elseif is_clamped_correlation(obs)              # a degenerate +-1, not a measurement
                clamped_region_count += 1
                continue
            end
            nscored += 1

            sc_abl = lro_scores(obs, predictive_draws(view(a.rho, r, :), n_eff_i; rng = rng_pred))
            sc_flr = lro_scores(obs, predictive_draws(view(flr,   r, :), n_eff_i; rng = rng_pred))
            for (arm, sc) in ((arm_abl, sc_abl), (arm_flr, sc_flr))
                arm.inside[r]   += sc.inside ? 1 : 0
                arm.n[r]        += 1
                arm.logscore[r] += sc.logscore
                arm.crps[r]     += sc.crps
                arm.width[r]    += sc.width
                push!(get!(arm.ls_by_imsize,  key, Float64[]), sc.logscore)
                push!(get!(arm.cov_by_imsize, key, Float64[]), sc.inside ? 1.0 : 0.0)
            end
            hit_abl += sc_abl.inside ? 1 : 0
            hit_flr += sc_flr.inside ? 1 : 0

            # THE TRUTH CROSS-CHECK ONLY SIMULATION PERMITS: parameter coverage of the DRAWN
            # latent value, with NO observation noise added, in atom-free Gaussian space (R-2).
            lo, hi = _p12_interval(view(a.z, r, :), P12_COVERAGE_NOMINAL)
            param_n[r] += 1
            (lo <= ztrue[r] <= hi) && (param_inside[r] += 1)
        end
        nscored > 0 || error(
            "run_p12_coverage: dataset $i had ZERO scorable regions. That is a defect in the " *
            "run, not a coverage of 0.")
        push!(arm_abl.per_dataset_cov, hit_abl / nscored)
        push!(arm_flr.per_dataset_cov, hit_flr / nscored)

        if i % 25 == 0 || i == n_datasets
            println("      $(lpad(i, 4))/$(n_datasets)   " *
                    "elapsed $(round((time() - t0) / 60; digits = 2)) min")
        end
    end

    pooled_abl = _p12cov_pool(arm_abl, G)
    pooled_flr = _p12cov_pool(arm_flr, G)
    param_used = findall(>(0), param_n)
    param_coverage = sum(param_inside[param_used]) / sum(param_n[param_used])
    param_coverage_per_region = [param_n[r] > 0 ? param_inside[r] / param_n[r] : NaN
                                 for r in 1:nreg]

    # --- 4. Intervals, the gate, and the pre-declared branch ----------------------------------
    ci_cluster = _p12cov_cluster_interval(pooled_abl.per_dataset_coverage)
    ci_wilson_271 = _p12cov_wilson(round(Int, pooled_abl.coverage * n_datasets),
                                   n_datasets)
    ci_invalid = _p12cov_wilson(sum(arm_abl.inside), pooled_abl.n_scored)

    stage2_gate = res.source === :descope_ablation_only ? :not_applicable_descope : :fail
    disqualified_arm = nothing
    headline_logscore_delta = nothing
    if res.source !== :descope_ablation_only && haskey(bundles, :car)
        error("run_p12_coverage: a two-arm route reached this runner, which the D-13 descope " *
              "route was built for. Refusing to guess how to pair the arms.")
    end

    discriminator_delta = pooled_abl.mean_logscore - pooled_flr.mean_logscore
    branch, branch_text = _p12cov_predeclared_branch(pooled_abl.coverage,
                                                     pooled_abl.mean_logscore,
                                                     pooled_flr.mean_logscore)

    # PITFALL 4: a shared miscalibration across arms that differ only in spatial correlation
    # cannot be caused by the thing that differs between them -- suspect the MASKING augmentation
    # before the prior. On this route there is only ONE model arm, so the "both arms" form of the
    # check cannot be evaluated; the single-arm form is recorded together with the rates a reader
    # needs, and the limitation is named rather than hidden.
    suspect_mask_ood = abs(pooled_abl.coverage - P12_COVERAGE_NOMINAL) >
                       2 * P12_STAGE2_COVERAGE_TOST_DELTA

    elapsed_min = (time() - t0) / 60

    println("\n[3/4] RESULTS")
    println("  " * "-"^74)
    for p in (pooled_abl, pooled_flr)
        println("  arm :$(rpad(string(p.name), 18)) coverage $(round(p.coverage; digits = 4))   " *
                "mean width $(round(p.mean_width; digits = 4))   " *
                "logscore $(round(p.mean_logscore; digits = 4))   " *
                "CRPS $(round(p.mean_crps; digits = 4))")
    end
    println("  " * "-"^74)
    println("  coverage 90% interval ON THE INDEPENDENT UNIT (n = $(n_datasets) DATASETS):")
    println("      [$(round(ci_cluster.lo; digits = 4)), $(round(ci_cluster.hi; digits = 4))]" *
            "   half-width $(round(ci_cluster.half_width; digits = 4))")
    println("  the SAME construction at the naive n = $(pooled_abl.n_scored) region-draws would be")
    println("      [$(round(ci_invalid.lo; digits = 4)), $(round(ci_invalid.hi; digits = 4))] " *
            "— INVALID. The 64 regions of one image share every nuisance and the global term, " *
            "so they are PSEUDO-REPLICATES; this interval is quoted ONLY to show the ~8x " *
            "tightening it would manufacture.")
    println("  parameter coverage (drawn z_field, NO observation noise, Gaussian space) = " *
            "$(round(param_coverage; digits = 4))")
    println("  predictive $(round(pooled_abl.coverage; digits = 4)) vs parameter " *
            "$(round(param_coverage; digits = 4)): predictive near nominal while parameter is " *
            "far from it would mean the observation-noise model is ABSORBING a miscalibrated " *
            "posterior — precisely the failure D-09 cannot detect on a real image.")
    println("  discriminator: ablation − prior_only_floor log score = " *
            "$(round(discriminator_delta; digits = 5)) nats/region " *
            "(bar: $(P12_STAGE2_LOGSCORE_MIN))")
    println("  stage2_gate = :$(stage2_gate)   (UNEXERCISED: there is no spatial arm)")
    println("  PRE-DECLARED BRANCH = :$branch")
    println("      $branch_text")

    # --- 5. Atomic write -----------------------------------------------------------------------
    println("\n[4/4] writing $(basename(report_path)) ...")
    payload = Dict{String,Any}(
        "schema_version" => 1,
        "bundle_source"  => res.source,
        "bundle_report"  => string(res.report),
        "bundles"        => bundle_meta,
        "bundle_arm"     => abl_name,
        "eval_arm"       => P12COV_EVAL_ARM,
        "eval_arm_note"  =>
            "NAMED LIMIT. The evaluation datasets are drawn at arm = :$(P12COV_EVAL_ARM), the " *
            "SAME joint the D-13 bundle was trained on (F5: train-joint must equal eval-joint). " *
            "These numbers are therefore a calibration statement UNDER THE MODEL'S OWN " *
            "GENERATIVE ASSUMPTIONS. They say nothing about behaviour under a spatially " *
            "correlated truth, which would be a misspecification question and is not asked here.",

        # --- per arm ---
        "coverage"             => pooled_abl.coverage,
        "coverage_per_region"  => pooled_abl.coverage_per_region,
        "logscore"             => pooled_abl.mean_logscore,
        "logscore_per_region"  => pooled_abl.logscore_per_region,
        "crps"                 => pooled_abl.mean_crps,
        "mean_width"           => pooled_abl.mean_width,
        "logscore_by_imsize"   => pooled_abl.logscore_by_imsize,
        "coverage_by_imsize"   => pooled_abl.coverage_by_imsize,
        "coverage_spread_across_regions" => pooled_abl.coverage_spread_across_regions,
        "per_dataset_coverage" => pooled_abl.per_dataset_coverage,

        "prior_only_floor_coverage"   => pooled_flr.coverage,
        "prior_only_floor_logscore"   => pooled_flr.mean_logscore,
        "prior_only_floor_crps"       => pooled_flr.mean_crps,
        "prior_only_floor_mean_width" => pooled_flr.mean_width,
        "prior_only_floor_logscore_by_imsize" => pooled_flr.logscore_by_imsize,

        "ablation_wording_rule" =>
            "The ablation arm is described ONLY as \"no spatial borrowing, full nuisance and " *
            "global borrowing\". D-10's phrase \"no mechanism to predict a held-out region " *
            "beyond the prior\" is FALSE as written: the other 63 observed regions still inform " *
            "the shared nuisances and the global term, so the ablation's predictive is strictly " *
            "narrower than the marginal prior. prior_only_floor is the genuinely prior-only arm.",

        # --- the truth cross-check ---
        "param_coverage"            => param_coverage,
        "param_coverage_per_region" => param_coverage_per_region,
        "param_coverage_note" =>
            "Coverage of the DRAWN z_field value at the held-out region against the posterior, " *
            "with NO observation noise added, in atom-free GAUSSIAN space (R-2). READ IT " *
            "BESIDE the predictive coverage: predictive near nominal while parameter is far " *
            "from nominal means the observation-noise model is absorbing a miscalibrated " *
            "posterior, which is exactly the failure D-09 cannot detect on a real image. This " *
            "comparison is the deliverable of this plan.",

        # --- the noise model ---
        "n_eff"                => n_eff,
        "n_eff_by_imsize"      => cal.by_imsize,
        "varz"                 => cal.varz,
        "varz_by_imsize"       => cal.varz_by_imsize,
        "varz_theta_spread_by_imsize" => cal.varz_theta_spread_by_imsize,
        "varz_fit_residual"    => cal.varz_fit_residual,
        "n_eff_applied" => :per_imsize,
        "n_eff_note" =>
            "P12_FISHERZ_NEFF is a DECLARED MODELLING ASSUMPTION, not a threshold. It " *
            "parametrizes the observation-noise model whose JOINT with the posterior D-09 " *
            "validates; it gates nothing by itself and no branch compares anything against it. " *
            "THE SWEEP APPLIED THE PER-SIZE CONSTANT (`n_eff_applied = :per_imsize`): the image " *
            "size is OBSERVED at read time, not latent, and n_eff moves by an order of magnitude " *
            "across P12_IMSIZE_SET's 16x pixel-count range. One pooled constant would make every " *
            "small-image interval too narrow and every large-image interval too wide, and the " *
            "pooled coverage could then sit at nominal as the average of two opposite " *
            "miscalibrations. The pooled value is appended and reported, but NOT applied.",

        # --- sizing and the independent unit ---
        "n_datasets"              => n_datasets,
        "smoke"                   => smoke,
        "n_draws_per_region"      => P12COV_N_DRAWS,
        "effective_independent_n" => n_datasets,
        "n_region_draws"          => pooled_abl.n_scored,
        "pseudoreplication_note" =>
            "THE INDEPENDENT UNIT IS THE DATASET. This run scored $(pooled_abl.n_scored) " *
            "region-draws from $(n_datasets) datasets, but the 64 regions of one image " *
            "share every nuisance and the global term, so the region-draws are " *
            "PSEUDO-REPLICATES. An interval built on them would be ~8x too tight and would " *
            "manufacture significance out of pseudo-replication. Never quote the region-draw " *
            "count without this sentence.",
        "coverage_ci_cluster" => [ci_cluster.lo, ci_cluster.hi],
        "coverage_ci_cluster_note" =>
            "THE INTERVAL THAT MAY BE QUOTED: mean +- z*sd/sqrt(n) over the PER-DATASET coverage " *
            "fractions, n = $(n_datasets) datasets. The between-dataset spread absorbs " *
            "the within-image correlation instead of assuming it away.",
        "coverage_ci_wilson_271" => [ci_wilson_271.lo, ci_wilson_271.hi],
        "coverage_ci_invalid_pseudoreplicated" => [ci_invalid.lo, ci_invalid.hi],
        "coverage_ci_invalid_note" =>
            "INVALID — RECORDED ONLY TO SHOW WHAT PSEUDO-REPLICATION WOULD MANUFACTURE. Wilson " *
            "at the region-draw count treats 64 regions of one image as 64 independent units. " *
            "DO NOT QUOTE IT.",

        # --- degeneracy accounting ---
        "clamped_region_count" => clamped_region_count,
        "masked_region_count"  => masked_region_count,
        "degeneracy_note" =>
            "A clamped +-1 correlation is a DEGENERATE region, not a measurement, and a " *
            "mask-row-zero region was never scorable at all. Both are EXCLUDED from every " *
            "coverage and score tally and counted here instead.",
        "realized_imsize_counts" => imsize_counts,

        # --- the gate ---
        "nominal"        => P12_COVERAGE_NOMINAL,
        "tost_delta"     => P12_STAGE2_COVERAGE_TOST_DELTA,
        "logscore_min"   => P12_STAGE2_LOGSCORE_MIN,
        "stage2_gate"    => stage2_gate,
        "stage2_gate_note" =>
            "UNEXERCISED BY CONSTRUCTION. lro_pass compares a spatial arm against the ablation " *
            "and there is no spatial arm on this route: the mini-spike returned " *
            "NONE-BEATS-ABLATION and 12-17 never ran. The spatial arm is ABSENT BY DESCOPE, not " *
            "missing by error, and a \"fail\" here would be read as evidence against the " *
            "spatial prior, which would be wrong. lro_pass is BUILT and unit-tested on " *
            "synthetic inputs; a green test is NOT the gate having run.",
        "disqualified_arm" => disqualified_arm,
        "headline_logscore_delta" => headline_logscore_delta,
        "headline_logscore_delta_note" =>
            "`nothing` BY CONSTRUCTION: this key is the SPATIAL-minus-ablation mean log-score " *
            "advantage and there is no spatial arm. 12-20's Guard 3 reads it, and 12-20 does " *
            "not run on this route. The comparison that WAS made is " *
            "`discriminator_logscore_delta` below, which is a different quantity and must not " *
            "be substituted for this one.",
        "discriminator_logscore_delta" => discriminator_delta,
        "discriminator_note" =>
            "ablation minus prior_only_floor, in nats per held-out region. This is the " *
            "PRE-DECLARED discriminator (12-16-PREDECLARATION.md section 2): coverage alone " *
            "cannot separate CALIBRATED AND INFORMATIVE from CALIBRATED BECAUSE UNINFORMATIVE, " *
            "and the margin over the floor is what does. The bar is the EXISTING " *
            "P12_STAGE2_LOGSCORE_MIN = $(P12_STAGE2_LOGSCORE_MIN), reused rather than invented.",

        # --- the pre-declared reading ---
        "predeclared_branch"      => branch,
        "predeclared_branch_text" => branch_text,
        "predeclaration_doc" =>
            ".planning/phases/12-spatial-colocalization-map/12-16-PREDECLARATION.md (df9ef91, " *
            "committed 2026-07-31T22:56:50, BEFORE this run existed)",

        # --- Pitfall 4 ---
        "suspect_mask_ood" => suspect_mask_ood,
        "suspect_mask_ood_note" =>
            "SINGLE-ARM FORM. Pitfall 4's check is that BOTH arms far off nominal implicates the " *
            "MCAR masking augmentation rather than the prior. There is only one model arm here, " *
            "so that form cannot be evaluated; this flag fires on the ablation alone at twice " *
            "the tolerance. The rates a reader needs are recorded beside it.",
        "training_realized_mask_rate" => bundle.realized_mask_rate,
        "d09_holdout_mask_rate"       => 1 / nreg,
        "mask_rate_note" =>
            "The net saw a realized mask rate of $(bundle.realized_mask_rate) in training; " *
            "D-09's held-out configuration masks exactly ONE of $nreg regions, a rate of " *
            "$(1 / nreg). If coverage is far off nominal, look at this pair BEFORE the prior.",

        # --- deferred scope, named ---
        "spat07_scope" => :deferred_to_v2_1,
        "spat08_scope" => :deferred_to_v2_1,
        "deferred_scope_note" =>
            "USER RULING, 2026-07-31 (12-D13-AUTHORISATION.md sections 7 and 9). NOT DELIVERED: " *
            "no Gaussian-space SBC, no randomized-rank rho-space SBC, no nuisance-appropriate " *
            "equivalence testing (SPAT-07); no radial-energy, offset-grid or eps = 0 guard " *
            "(SPAT-08). THE PHASE THEREFORE MAKES NO PER-REGION CALIBRATION CLAIM BEYOND " *
            "PREDICTIVE COVERAGE, and NO radial-confound claim -- the S-4 warning is left " *
            "standing and unaddressed. Predictive coverage and SBC calibration are different " *
            "statements and must not be conflated.",
        "scale_limit" =>
            "SPIKE SCALE. The bundle scored here is the D-13 fallback trained on " *
            "P12_MINISPIKE_N = 10,000 pairs, not the 50,000-pair version 12-17 would have " *
            "built. Any claim from this run inherits that limit.",

        "elapsed_min"   => elapsed_min,
        "generated"     => string(Dates.now(Dates.UTC)) * "Z",
        "julia_version" => string(VERSION),
        "caption" =>
            "D-09 leave-region-out predictive coverage and scoring, SIMULATION arm, at " *
            "N = P12_STAGE2_N_MIN = $(n_datasets) datasets x $nreg regions on the " *
            "reserved P12_COVERAGE_COUNTER stream. Two arms: the D-13 ablation (NO SPATIAL " *
            "BORROWING, FULL NUISANCE AND GLOBAL BORROWING) and the genuinely prior-only floor. " *
            "The spatial arm is ABSENT BY DESCOPE. Coverage is reported beside interval WIDTH " *
            "and the log score beside prior_only_floor, because neither can be read without " *
            "its baseline. The independent unit is the DATASET.")

    mkpath(dirname(report_path))
    tmp = report_path * ".tmp"
    jldsave(tmp; (Symbol(k) => v for (k, v) in payload)...)
    JLD2.jldopen(tmp, "r") do f
        for k in ("stage2_gate", "param_coverage", "spat07_scope", "logscore_by_imsize",
                  "headline_logscore_delta", "n_eff", "coverage", "mean_width",
                  "prior_only_floor_logscore", "effective_independent_n", "n_region_draws",
                  "clamped_region_count", "predeclared_branch")
            haskey(f, k) || error("run_p12_coverage: integrity check failed, $tmp missing `$k`")
        end
    end
    mv(tmp, report_path; force = true)
    println("      -> $(report_path)")
    println("\nDONE in $(round(elapsed_min; digits = 2)) min.")
    return payload
end

isdefined(@__MODULE__, :P12_COVERAGE_LOAD_ONLY) || main()
