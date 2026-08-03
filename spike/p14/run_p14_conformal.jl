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

# spike/p14/run_p14_conformal.jl --- REPORTED SC1-d: realized split-conformal marginal coverage,
# plus the SHARED evaluation pool three later runners consume.
#
# =============================================================================================
# ANTI-SNOOPING BANNER. READ BEFORE CHANGING ANY NUMBER IN THIS FILE.
# =============================================================================================
# EVERY THRESHOLD THIS RUNNER SCORES AGAINST WAS COMMITTED TO `spike/p14/consts.jl` BEFORE ANY
# PHASE-14 RESULT EXISTED. The miscoverage level, the two sample sizes, the two reserved stream
# counters, the vacuity floor and the iteration allowance are all READ from that frozen file by
# name. Not one of them is a literal here, and the acceptance criteria assert that, so a value
# cannot be "temporarily" inlined and left behind. The SC1-d band itself is not even a constant:
# it is DERIVED from the miscoverage level and the evaluation size by formula, which is precisely
# what makes it un-choosable.
#
# ROADMAP SC1 IS AMENDED BY D-02 and SC2 BY D-05, both frozen in `14-CONTEXT.md` while NO
# Phase-14 result of any kind existed. ANY REPORT THAT CITES THIS RESULT MUST CITE THE ORIGINAL
# ROADMAP CRITERIA ALONGSIDE IT. The amendment notice travels INSIDE the artifact so a reader
# who never opens a planning document still meets it.
#
# AN HONEST SHORTFALL IS A PHASE-14 FINDING, NEVER A LICENCE TO ACT. It is not a reason to
# loosen the miscoverage level, to redraw the calibration or evaluation set, to reseed, or to
# re-scope the band. `P14_ITERATION_ALLOWANCE = 1` has exactly ONE pre-declared trigger, written
# into `P14_ITERATION_TRIGGER` before any result existed, and it is a SCORE SWITCH (LAC -> APS,
# re-run once) -- never a threshold relaxation. Phase 7 was amended TWICE after seeing results
# and the credibility cost of that is why this paragraph is here. Phase 11 and Phase 12 each
# missed a pre-registered bar and REPORTED it; that is the precedent this runner follows.
#
# THE ORDER OF OPERATIONS IS THE POINT: read the thresholds, compute the statistics, PERSIST THE
# ARTIFACTS, print the headline, then assert. Both artifacts are written before anything can
# throw, so a FAILING gate still leaves a complete, self-describing report on disk.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY AND NOT A TEST. It composes surfaces that already exist
# and are already tested per file. Including it does nothing; it must be run deliberately:
#
#     julia --project=spike -t auto spike/p14/run_p14_conformal.jl
#
# WHAT THIS RUN'S GUARANTEE IS, IN FULL, AND WHERE THE INTENDED BOUND WENT.
# The calibration draws and the evaluation draws come from the SAME simulator the evidence net
# was trained on, at two reserved counters that are disjoint from the training stream. The
# training joint therefore equals the evaluation joint BY CONSTRUCTION, so every number this
# runner reports is a WELL-SPECIFIED-REGIME number that inherits the simulator's misspecification
# in full. D-03 intended a real-data bound beside it, drawn from the project's curated corpus
# collection. D-03a withdrew that half after checking it directly against corpus/manifest.csv:
# 30 of the 32 rows are `tier: simulated-secondary` (the CBS benchmark -- i.e. a SECOND
# simulator, not physical ground truth), the only 2 `tier: physical-primary` rows are the
# `sealed_holdout` reserved for Phase 16's blind evaluation, and every row records `bytes = 0`
# because `.gitignore:451-453` keeps the image bytes out of git and they are unfetched by design.
# So the bound D-03 intended is ABSENT, not merely loose -- and this file says so in the artifact
# rather than leaving a reader to infer that a simulator-derived number was checked against
# reality. The sanctioned real substrate is the six committed microscopy TIFFs under
# `test/test_images/` (D-03a, following `spike/p13/real_images.jl:71,145`), and that check belongs
# to the SC1-f runner, not to this one.
#
# THAT PARAGRAPH LIVES IN A `#` BLOCK ON PURPOSE. `spike/test/test_p14_decoupling.jl` strips
# whole-line comments and then bans the corpus tokens outright in every non-allowlisted Phase-14
# source, because a docstring or a string literal is CODE to that scan while a `#` block is not.
# This runner is not on the allowlist and must not be added to one -- an allowlist entry is a
# decision, not a fix. The persisted `circularity_note` therefore carries the same substance in
# the same detail without the banned tokens, and points here. This is the third recorded instance
# of "the guard must be allowed to name what it forbids" landing on this phase
# (14-02, 14-05 deviation 2, and now here).
#
# NO NEW DEPENDENCY, NO INSTALL, CPU-ONLY (CLAUDE.md, D-02). The hedge is `sort` + `ceil` +
# integer indexing; the draw, the net, the density null and the fusion all already exist and are
# consumed unchanged. `ConformalPrediction.jl` is NOT added (D-02, and SC1 is amended for it).
#
# DECOUPLING: spike-local. `src/` is reached only READ-ONLY and transitively, and step 0 PROVES
# at run time that it is byte-unchanged AND that no untracked file has appeared under it.
#
# RUN THE UNIT FILES PER FILE, never through the aggregate suite -- it exits 1 early at the
# Phase-4 speedup gate and masks every later include block.

using JLD2
using Statistics
using Dates

# --- Guarded includes, in dependency order --------------------------------------------------------
# `decide.jl` transitively pulls the pre-registration, the provenance loader, the posterior, the
# FDR rule, the hedge, the fusion, the result type, the Phase-13 read surface, the SHIPPED OOD
# channel and the classical comparator. `pools.jl` is the one surface it does not reach, and it is
# the one that owns the unstratified draw. Every include is guarded, so listing the modules this
# runner names directly costs nothing and documents the dependency set.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :P14_PROVENANCE_LOADED)   || include(joinpath(@__DIR__, "provenance.jl"))
isdefined(@__MODULE__, :P14_POSTERIOR_LOADED)    || include(joinpath(@__DIR__, "posterior.jl"))
isdefined(@__MODULE__, :p14_bayes_fdr)           || include(joinpath(@__DIR__, "fdr.jl"))
isdefined(@__MODULE__, :P14_CONFORMAL_LOADED)    || include(joinpath(@__DIR__, "conformal.jl"))
isdefined(@__MODULE__, :P14_FUSE_LOADED)         || include(joinpath(@__DIR__, "fuse.jl"))
isdefined(@__MODULE__, :P14Result)               || include(joinpath(@__DIR__, "result.jl"))
isdefined(@__MODULE__, :P14_DECIDE_LOADED)       || include(joinpath(@__DIR__, "decide.jl"))
isdefined(@__MODULE__, :P14_POOLS_LOADED)        || include(joinpath(@__DIR__, "pools.jl"))

if !isdefined(@__MODULE__, :P14_CONFORMAL_RUNNER_LOADED)
    "The reported SC1-d artifact. Written BEFORE the coverage assertion can throw."
    const P14_CONFORMAL_REPORT_PATH = joinpath(@__DIR__, "p14_conformal_report.jld2")

    """
    The SHARED evaluation pool. SC1-b (realized FDP), SC3-a/b/c (risk-coverage) and the SC3-d
    in-distribution arm all read THIS one set, so those numbers are computed on the SAME items at
    the SAME reserved counter rather than on three independent redraws. Three redraws would make a
    risk-coverage curve and the coverage it is scored at into two different experiments.
    """
    const P14_EVAL_POOL_PATH = joinpath(@__DIR__, "p14_eval_pool.jld2")

    # THE SC1-d BAND'S STANDARD-ERROR MULTIPLIER. This is the one number in this file that is not
    # read from `spike/p14/consts.jl`, and it is not a threshold: it is part of the BAND FORMULA
    # itself, transcribed from `14-VALIDATION.md` row SC1-d and restated in `14-09-PLAN.md`'s
    # objective, both of which were frozen before this runner existed. It is named here rather than
    # buried inside the expression so that it is greppable and so that a reader can see there is
    # exactly one such transcription.
    const P14_SC1D_BAND_SD_MULTIPLIER = 3

    """
    The placeholder hedge threshold used while the CALIBRATION posteriors are computed.

    The calibration pass has to reach the class posterior, and the only surface that produces one
    is the composed per-pair path -- which requires a threshold it does not yet have, because the
    threshold is what the calibration pass exists to fit. Rather than write a SECOND spelling of
    the encode/forward/posterior chain (two spellings of one computation is how a permutation bug
    survives, which cost this project a headline finding in Phase 12), the composed path is called
    once with this placeholder and its conformal, fusion and decision fields are DISCARDED.

    The class posterior does not depend on it, and that is PROVED at run time rather than argued:
    a sample of calibration items is recomputed at `P14_CAL_QHAT_PLACEHOLDER_ALT` and the
    posteriors are asserted identical.
    """
    const P14_CAL_QHAT_PLACEHOLDER     = 1.0
    "The second placeholder, used only to prove the calibration posterior is independent of it."
    const P14_CAL_QHAT_PLACEHOLDER_ALT = 0.5
    "How many calibration items the independence proof recomputes."
    const P14_CAL_INDEPENDENCE_SAMPLE  = 8

    """
    The value stored in the evaluation pool's `abstain_reason` column for an item that did NOT
    abstain. A per-pair result carries `nothing` there; a `Union{Nothing,Symbol}` column would be
    a deserialization surface for no benefit, and an EMPTY symbol would read as a missing reason
    rather than as an absent one. The sentinel is stored beside the column under
    `abstain_reason_sentinel` so a reader never has to guess it.
    """
    const P14_ABSTAIN_REASON_SENTINEL = :not_abstained

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else can
    # make this body skip.
    const P14_CONFORMAL_RUNNER_LOADED = true
end

# =============================================================================================
# Small local helpers
# =============================================================================================

"""
    p14_sc1d_band_lower(alpha, n) -> Float64

The pre-registered SC1-d lower band, `1 - alpha - k * sqrt(alpha * (1 - alpha) / n)`, computed
from the frozen miscoverage level and the frozen evaluation size.

# It is DERIVED, which is what makes it un-choosable

Split conformal's realized coverage on a fresh evaluation set is Binomial, so its standard error
is `sqrt(alpha * (1 - alpha) / n)` and the band is that error multiplied by
`P14_SC1D_BAND_SD_MULTIPLIER`. Nobody gets to pick the resulting number: fix the level and the
size and the band follows. That is exactly why `P14_ITERATION_TRIGGER` names the band as
non-adjustable -- the only way to move it is to move the level or the size, and both are frozen.
"""
p14_sc1d_band_lower(alpha::Real, n::Integer) =
    1.0 - Float64(alpha) -
    P14_SC1D_BAND_SD_MULTIPLIER * sqrt(Float64(alpha) * (1.0 - Float64(alpha)) / n)

"""
    _p14_class_key(c::ThreeWayClass) -> Symbol

The frozen class SYMBOL for a three-way label, by an explicit three-branch map.

# Why this is written out rather than derived from the enum's name

Five class orderings coexist in this codebase and every headline statistic -- coverage included --
is invariant to a CONSISTENT relabelling, which is what makes an inconsistent one produce numbers
that look right. A lowercased `string(c)` would be a silent coupling to the enum's spelling; an
explicit map is auditable in one glance. It is also not taken on trust: `_p14_assert_class_keys`
pins the map against the frozen `class_masses` naming on every draw, so a swapped pair fails on
the draw rather than travelling into a coverage number.
"""
function _p14_class_key(c::ThreeWayClass)
    c === COLOC     && return :coloc
    c === RANDOM    && return :random
    c === EXCLUSION && return :exclusion
    error("_p14_class_key: unhandled three-way class $(c). The class-symbol map must be total; " *
          "a fall-through would silently drop a class out of every rate this runner reports.")
end

"""
    _p14_assert_class_keys(masses, keys) -> nothing

Pin [`_p14_class_key`](@ref) to the frozen `class_masses` naming: the per-symbol frequency of the
mapped keys must equal the mass the label surface itself reports for that name.

This catches the one failure the map can produce and the one no downstream number would notice: a
SWAPPED pair. Coverage, the set-size rates and the class-mass tables would all still look
perfectly ordinary under a consistent swap.
"""
function _p14_assert_class_keys(masses::NamedTuple, keys::AbstractVector{Symbol})
    n = length(keys)
    for k in P14_MASS_KEYS
        got = count(==(k), keys) / n
        want = Float64(getproperty(masses, k))
        isapprox(got, want; atol = 1e-12) || error(
            "_p14_assert_class_keys: the class-symbol map puts $(round(got; digits = 6)) of the " *
            "draw on `$k` while the frozen label surface reports $(round(want; digits = 6)). The " *
            "map and `class_masses` disagree, which means a class name is attached to the wrong " *
            "label -- and a CONSISTENT swap is invisible to every rate this runner reports.")
    end
    return nothing
end

"""
    _p14_save_report(path, required; kwargs...) -> String

Atomically persist an artifact: write `path * ".tmp"`, REOPEN it read-only and integrity-check
every key in `required`, then `mv(...; force = true)`.

The `_p13_gate_save_report` idiom (`spike/p13/run_three_way_gate.jl:171-190`), with the required
keys passed in because this runner writes TWO artifacts with different contents and one shared
save path. A hard-coded key list would have had to name the union of both, which would either
pass vacuously on each or fail on both.

A crash mid-write leaves a discardable `.tmp` rather than a torn artifact a later reader would
happily believe.
"""
function _p14_save_report(path, required; kwargs...)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        for k in required
            @assert haskey(f, k) "_p14_save_report: integrity check failed -- $tmp is missing the required key `$k`"
        end
    end
    mv(tmp, path; force = true)
    return path
end

"""
    _p14_cross_tab(rows, row_levels, cols, col_levels) -> Matrix{Int}

A counts matrix over two symbol vectors at declared levels, so a cross-tabulation is never read
off an implicit ordering.

Built for the one cross-tabulation `14-RESEARCH` B.3 asks for by name: the EMPTY conformal set
against the OOD state. An empty set means the point is more nonconforming than `(1 - alpha)` of
the whole calibration set -- a distribution-free misspecification signal arriving on a completely
different channel from the density detector. If the two agree that is corroboration; if they
disagree that is the more interesting finding. Neither can be seen from marginal rates.
"""
function _p14_cross_tab(rows::AbstractVector{Symbol}, row_levels,
                        cols::AbstractVector{Symbol}, col_levels)
    length(rows) == length(cols) || throw(ArgumentError(
        "_p14_cross_tab: the two vectors must be the same length; got $(length(rows)) and $(length(cols))."))
    M = zeros(Int, length(row_levels), length(col_levels))
    for i in eachindex(rows)
        r = findfirst(==(rows[i]), row_levels)
        c = findfirst(==(cols[i]), col_levels)
        (r === nothing || c === nothing) && error(
            "_p14_cross_tab: the pair ($(rows[i]), $(cols[i])) is outside the declared levels " *
            "$(row_levels) x $(col_levels). A cell counted into nothing would make the table sum " *
            "to less than n with nothing saying why.")
        M[r, c] += 1
    end
    @assert sum(M) == length(rows) "_p14_cross_tab: the table does not sum to n"
    return M
end

# =============================================================================================
# The reported run
# =============================================================================================

"""
    main(; n_cal = P14_N_CAL, n_eval = P14_N_EVAL, ...) -> NamedTuple

The reported SC1-d measurement, end to end and once, plus the shared evaluation pool.

`n_cal` and `n_eval` default to the FROZEN sizes. Passing anything else puts the run in SMOKE
mode: it writes to `_smoke` artifact paths and prints a banner saying so, so a load-check can
never be mistaken for, or quoted as, a verdict, and can never overwrite the reported artifact.
"""
function main(; n_cal::Integer = P14_N_CAL,
                n_eval::Integer = P14_N_EVAL,
                report_path = P14_CONFORMAL_REPORT_PATH,
                pool_path   = P14_EVAL_POOL_PATH,
                verbose::Bool = true)

    t_start  = time()
    reported = (n_cal == P14_N_CAL) && (n_eval == P14_N_EVAL)
    if !reported
        report_path = replace(report_path, ".jld2" => "_smoke.jld2")
        pool_path   = replace(pool_path,   ".jld2" => "_smoke.jld2")
    end

    # --- 0. THE DECOUPLING PROOF, AT RUN TIME (D-01 / D-02) ---------------------------------
    # Asserted WHILE the run happens, not only in review: a long reported run that started on a
    # clean tree and finished on a dirty one would be a repudiation hole (T-14-11). The shared
    # lane guard is called first; the two commands are then re-run here so the runner's own
    # failure message names the breach, and a THIRD check is added that neither the guard nor a
    # diff can make -- an untracked NEW file under `src/` is invisible to `git diff HEAD` and
    # would be exactly the productionization-by-stealth D-01 forbids.
    lane = p14_assert_lane_clean(; verbose = verbose)
    if lane.checked
        src_clean = success(Cmd(`git diff --quiet HEAD -- src`; dir = P14_REPO_ROOT))
        env_clean = success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                                dir = P14_REPO_ROOT))
        src_untracked = readchomp(Cmd(`git status --porcelain -- src`; dir = P14_REPO_ROOT))
        @assert src_clean "D-01 decoupling breach: `git diff --quiet HEAD -- src` FAILED -- src/ is not byte-unchanged, so this run may not be reported"
        @assert env_clean "D-02 decoupling breach: spike/Project.toml or spike/Manifest.toml is modified -- the frozen spike environment moved under the run, and the hedge needs no dependency at all"
        @assert isempty(src_untracked) "D-01 decoupling breach: an UNTRACKED file appeared under src/ ($(src_untracked)) -- a diff against HEAD cannot see it, and shipping the decision layer by stealth is precisely what D-01 forbids"
    end

    # --- 1. THE INHERITANCE (D-07) AND THE DENSITY NULL --------------------------------------
    bundle = p14_load_bundle()
    ref    = p14_ood_reference()
    ood_in = p14_ood_input(ref)
    basis  = load_p13_basis()

    if !ref.available
        # NOT an abort. A null that was never fitted is the NOT-CHECKED state, every item then
        # abstains by default, and that is itself a reportable D-06 outcome rather than a crash.
        println("!"^78)
        println("OOD REFERENCE UNAVAILABLE -- every evaluation item will resolve to :not_checked")
        println("and ABSTAIN by default (D-06). This is a REPORTABLE OUTCOME, not a failure of")
        println("this runner, and the coverage number below is unaffected by it: coverage is a")
        println("property of the conformal set, not of the decision.")
        println("reason: ", ref.note)
        println("!"^78)
    end

    # --- 2. THE LOCKED-THRESHOLD BANNER, BEFORE ANYTHING IS COMPUTED -------------------------
    band_lower = p14_sc1d_band_lower(P14_ALPHA_CONFORMAL, n_eval)
    if verbose
        println("="^78)
        println("PHASE-14 REPORTED SC1-d: SPLIT-CONFORMAL MARGINAL COVERAGE")
        println("="^78)
        reported || println("*** SMOKE MODE -- NOT A REPORTED RUN, THE BAND IS NOT ASSERTED ***")
        println("AMENDED CRITERIA: ", P14_AMENDMENT_NOTICE)
        println("GATED (SC1-d):")
        println("  P14_ALPHA_CONFORMAL     = $P14_ALPHA_CONFORMAL          " *
                "(spike/p14/consts.jl section C; the hedge's miscoverage level,")
        println("                            a PRE-REGISTERED CONSTANT and NEVER alpha_FDR,")
        println("                            which is a per-call user parameter -- D-07)")
        println("  P14_N_CAL               = $P14_N_CAL           (consts.jl section B; requested n_cal = $n_cal)")
        println("  P14_N_EVAL              = $P14_N_EVAL           (consts.jl section B; requested n_eval = $n_eval)")
        println("  SC1-d band lower        = $band_lower")
        println("                            = 1 - alpha - $P14_SC1D_BAND_SD_MULTIPLIER * sqrt(alpha*(1-alpha)/n_eval),")
        println("                            DERIVED from the two frozen constants by formula")
        println("                            (14-VALIDATION.md row SC1-d). Not a chosen number.")
        println("REPORTED, NOT GATED:")
        println("  P14_AMBIGUOUS_RATE_FLOOR = $P14_AMBIGUOUS_RATE_FLOOR         (consts.jl section F) -- a VACUITY LABEL,")
        println("                            never a verdict. Abstaining rarely can also be the")
        println("                            CORRECT behaviour on well-specified data; what must")
        println("                            not happen is a coverage number quoted without it.")
        println("  q-hat, the realized set-size distribution, the class-mass tables and the")
        println("  conformal-status x OOD-state cross-tabulation are all reported beside the")
        println("  coverage number and none of them is gated.")
        println("STREAM (D-01, both draws UNSTRATIFIED -- 14-RESEARCH Pitfall 1):")
        println("  P14_DEV_SEED            = $(repr(UInt64(P14_DEV_SEED)))")
        println("  P14_CAL_COUNTER         = $P14_CAL_COUNTER   (calibration)")
        println("  P14_EVAL_COUNTER        = $P14_EVAL_COUNTER   (the SHARED reported evaluation set)")
        println("ITERATION DISCIPLINE:")
        println("  P14_ITERATION_ALLOWANCE = $P14_ITERATION_ALLOWANCE")
        println(P14_ITERATION_TRIGGER)
        println("INHERITED (D-07, LOADED and four-way agreed, never re-derived):")
        println("  tau                     = $(bundle.tau)")
        println("  class prior (MEASURED)  = $(bundle.prior)")
        println("  Phase-13 coloc-head calibration carried forward: ECE = $(bundle.calibration.ece)")
        println("OOD CHANNEL:")
        println("  available               = $(ref.available)   thr = $(ref.thr)")
        println("  channels wired          = $(ref.channels_wired)   " *
                "NOT wired = $(ref.channels_not_wired)")
        println("="^78)
    end

    # --- 3. THE CALIBRATION SET ---------------------------------------------------------------
    verbose && println("[1/6] drawing the UNSTRATIFIED calibration set (n_cal = $n_cal) …")
    cal        = p14_draw_pool(n_cal; counter = P14_CAL_COUNTER, basis = basis)
    cal_masses = p14_assert_unstratified([it.class for it in cal])
    cal_keys   = [_p14_class_key(it.class) for it in cal]
    _p14_assert_class_keys(cal_masses, cal_keys)
    verbose && println("      realized calibration class masses: $cal_masses")

    verbose && println("[2/6] scoring the calibration set and fitting q-hat …")
    cal_results    = [p14_decide_one(bundle, it.Zs, it.Zc, ood_in, P14_CAL_QHAT_PLACEHOLDER;
                                     lambda = it.lambda, idx = it.idx) for it in cal]
    cal_posteriors = [class_posterior(r) for r in cal_results]

    # THE PLACEHOLDER IS PROVED IRRELEVANT, NOT ASSERTED TO BE. See P14_CAL_QHAT_PLACEHOLDER.
    for j in 1:min(P14_CAL_INDEPENDENCE_SAMPLE, length(cal))
        it  = cal[j]
        alt = p14_decide_one(bundle, it.Zs, it.Zc, ood_in, P14_CAL_QHAT_PLACEHOLDER_ALT;
                             lambda = it.lambda, idx = it.idx)
        @assert class_posterior(alt) === cal_posteriors[j] "the calibration posterior moved when the placeholder hedge threshold moved (item $j) -- the composed path is not being used the way this runner claims"
    end

    fit  = p14_conformal_calibrate(cal_posteriors, cal_keys, P14_ALPHA_CONFORMAL)
    qhat = fit.qhat
    verbose && println("      q-hat = $qhat   scores (min/median/max) = $(fit.scores_summary)")

    # --- 4. THE SHARED EVALUATION SET ---------------------------------------------------------
    verbose && println("[3/6] drawing the UNSTRATIFIED evaluation set (n_eval = $n_eval) …")
    ev          = p14_draw_pool(n_eval; counter = P14_EVAL_COUNTER, basis = basis)
    eval_masses = p14_assert_unstratified([it.class for it in ev])
    eval_keys   = [_p14_class_key(it.class) for it in ev]
    _p14_assert_class_keys(eval_masses, eval_keys)
    verbose && println("      realized evaluation class masses: $eval_masses")

    verbose && println("[4/6] deciding every evaluation item at q-hat …")
    res = [p14_decide_one(bundle, it.Zs, it.Zc, ood_in, qhat;
                          lambda = it.lambda, idx = it.idx) for it in ev]

    # --- 5. THE MEASUREMENT --------------------------------------------------------------------
    covered  = [eval_keys[i] in conformal_set(res[i]) for i in eachindex(res)]
    coverage = count(covered) / length(covered)

    statuses = [conformal_status(r) for r in res]
    states   = [ood_state(r) for r in res]
    diag     = p14_hedge_diagnostics(statuses, qhat)
    set_size_counts = (singleton = count(==(:singleton), statuses),
                       ambiguous = count(==(:ambiguous), statuses),
                       empty     = count(==(:empty),     statuses))
    cross_tab = _p14_cross_tab(statuses, collect(P14_CONFORMAL_STATUSES),
                               states,   collect(P14_OOD_STATES))

    iteration_trigger_fired = coverage < band_lower

    # --- 6. PERSIST THE REPORT, BEFORE THE HEADLINE AND BEFORE THE ASSERTION -----------------
    verbose && println("[5/6] persisting the SC1-d report (BEFORE any verdict) …")
    _p14_save_report(report_path,
        ("qhat", "coverage", "band_lower", "cal_masses", "eval_masses", "vacuous_hedge",
         "guarantee_basis", "circularity_note", "amendment", "iteration_trigger_fired");
        # --- the measurement ---
        qhat = qhat, coverage = coverage, band_lower = band_lower,
        coverage_meets_band = !iteration_trigger_fired,
        coverage_gap = coverage - band_lower,
        n_covered = count(covered),
        n_cal = Int(n_cal), n_eval = Int(n_eval),
        cal_masses = cal_masses, eval_masses = eval_masses,
        scores_summary = fit.scores_summary,
        # --- the hedge, reported BESIDE the coverage number (14-RESEARCH Pitfall 9) ---
        set_size_counts = set_size_counts,
        singleton_rate = diag.singleton_rate,
        ambiguous_rate = diag.ambiguous_rate,
        empty_rate = diag.empty_rate,
        vacuous_hedge = diag.vacuous_hedge,
        cross_tab_status_by_ood = cross_tab,
        cross_tab_status_levels = collect(P14_CONFORMAL_STATUSES),
        cross_tab_ood_levels = collect(P14_OOD_STATES),
        # --- EVERY THRESHOLD THIS RUN WAS SCORED AGAINST, written INTO the artifact ---
        P14_ALPHA_CONFORMAL = P14_ALPHA_CONFORMAL,
        P14_N_CAL = P14_N_CAL, P14_N_EVAL = P14_N_EVAL,
        P14_CAL_COUNTER = Int(P14_CAL_COUNTER), P14_EVAL_COUNTER = Int(P14_EVAL_COUNTER),
        P14_AMBIGUOUS_RATE_FLOOR = P14_AMBIGUOUS_RATE_FLOOR,
        P14_SC1D_BAND_SD_MULTIPLIER = P14_SC1D_BAND_SD_MULTIPLIER,
        P14_ITERATION_ALLOWANCE = P14_ITERATION_ALLOWANCE,
        P14_ITERATION_TRIGGER = P14_ITERATION_TRIGGER,
        iteration_trigger_fired = iteration_trigger_fired,
        # --- what the guarantee actually rests on (D-03 / D-03a, T-14-32) ---
        guarantee_basis = :simulator_derived_held_out_draws,
        circularity_note =
            "SIMULATOR-DERIVED. The calibration draws and the evaluation draws come from the " *
            "SAME forward simulator the evidence net was trained on, at two reserved counters " *
            "disjoint from the training stream. The training joint therefore equals the " *
            "evaluation joint BY CONSTRUCTION, so the coverage number here is a " *
            "WELL-SPECIFIED-REGIME number that inherits the simulator's misspecification in " *
            "full. D-03 intended a real-data bound beside it; D-03a WITHDREW that half after " *
            "checking the intended substrate directly: 30 of its 32 rows are tier " *
            "`simulated-secondary` (a SECOND simulator, not physical ground truth), the only 2 " *
            "`physical-primary` rows are reserved unread for Phase 16's blind evaluation, and " *
            "every row records zero bytes because the image data is unfetched by design. So the " *
            "bound D-03 intended is ABSENT, NOT MERELY LOOSE. The sanctioned real substrate is " *
            "the six committed microscopy TIFFs under test/test_images/ (D-03a), and that is an " *
            "ILLUSTRATION run by the SC1-f runner, not a coverage claim and not run here. The " *
            "withdrawn feasibility figure D-03 quoted must not be carried forward. The full " *
            "statement, with the token names this artifact may not spell, is in the `#` banner " *
            "of spike/p14/run_p14_conformal.jl.",
        conformal_library = "none -- hand-rolled; SC1 AMENDED by D-02",
        amendment = P14_AMENDMENT_NOTICE,
        named_limits = p14_named_limits(),
        # --- provenance ---
        provenance = p14_provenance_record(bundle.prov),
        pool_provenance = p14_pool_provenance(ref),
        ood_available = ref.available,
        ood_note = ref.note,
        class_prior_used = bundle.prior,
        tau = bundle.tau,
        grid = bundle.grid,
        reported = reported,
        eval_pool_path = pool_path,
        julia_version = string(VERSION),
        nthreads = Threads.nthreads(),
        elapsed_s = time() - t_start,
        generated = string(Dates.now(Dates.UTC)) * "Z")
    verbose && println("      persisted (before any verdict) -> $report_path")

    # --- 7. PERSIST THE SHARED EVALUATION POOL ------------------------------------------------
    # Column-major, one array per field, so a later runner reads what it needs without
    # deserializing a vector of result objects. The raw summaries are deliberately NOT stored:
    # every item is a pure function of its index on this counter, so a runner that needs them
    # redraws them bitwise rather than carrying megabytes here.
    verbose && println("[6/6] persisting the SHARED evaluation pool …")
    _p14_save_report(pool_path,
        ("true_class", "class_posterior_coloc", "null_posterior", "conformal_status",
         "ood_state", "confidence", "decision", "abstain_reason", "idx");
        n = length(res),
        idx = [Int(it.idx) for it in ev],
        lambda = [Float64(it.lambda) for it in ev],
        rho_sample = [Float64(it.rho_s) for it in ev],
        rho_control = [Float64(it.rho_c) for it in ev],
        true_class = eval_keys,
        # the three-class posterior, BY NAME -- never a positional matrix
        class_posterior_coloc     = [class_posterior(r).coloc     for r in res],
        class_posterior_random    = [class_posterior(r).random    for r in res],
        class_posterior_exclusion = [class_posterior(r).exclusion for r in res],
        null_posterior = [null_posterior(r) for r in res],
        null_split_p_random    = [null_split(r).p_random    for r in res],
        null_split_p_exclusion = [null_split(r).p_exclusion for r in res],
        # the evidence, so a later runner can re-derive the posterior at another class prior
        # without re-running the net
        logbf_coloc     = [log_bf_vs_random(r).coloc     for r in res],
        logbf_exclusion = [log_bf_vs_random(r).exclusion for r in res],
        conformal_status = statuses,
        conformal_in_set_coloc     = [(:coloc     in conformal_set(r)) for r in res],
        conformal_in_set_random    = [(:random    in conformal_set(r)) for r in res],
        conformal_in_set_exclusion = [(:exclusion in conformal_set(r)) for r in res],
        covered = covered,
        ood_state = states,
        confidence = [p14_confidence(class_posterior(r)) for r in res],
        decision = [decision(r) for r in res],
        abstain_reason = [abstain_reason(r) === nothing ? P14_ABSTAIN_REASON_SENTINEL :
                          abstain_reason(r) for r in res],
        abstain_reason_sentinel = P14_ABSTAIN_REASON_SENTINEL,
        # the threshold these fields were produced at, so a set is never read without its cut
        qhat = qhat,
        P14_ALPHA_CONFORMAL = P14_ALPHA_CONFORMAL,
        P14_EVAL_COUNTER = Int(P14_EVAL_COUNTER),
        n_eval = Int(n_eval),
        eval_masses = eval_masses,
        class_prior_used = bundle.prior,
        tau = bundle.tau,
        ood_available = ref.available,
        stratified = false,
        redraw_recipe = "p14_draw_pool(n_eval; counter = P14_EVAL_COUNTER) reproduces these " *
                        "items bitwise on any thread count -- each item is a pure function of " *
                        "its index -- which is why the raw summaries are not stored here.",
        guarantee_basis = :simulator_derived_held_out_draws,
        amendment = P14_AMENDMENT_NOTICE,
        provenance = p14_provenance_record(bundle.prov),
        pool_provenance = p14_pool_provenance(ref),
        reported = reported,
        generated = string(Dates.now(Dates.UTC)) * "Z")
    verbose && println("      persisted -> $pool_path")

    # --- 8. THE HEADLINE, PRINTED BEFORE THE ASSERTION CAN THROW ------------------------------
    if verbose
        println("\n", "-"^78)
        println("SC1-d realized split-conformal marginal coverage: ",
                iteration_trigger_fired ? "BELOW THE BAND" : "WITHIN THE BAND")
        println("  q-hat                  = $qhat")
        println("  realized coverage      = $coverage   ($(count(covered)) of $(length(covered)))")
        println("  pre-registered band    = $band_lower  (lower bound)")
        println("  gap (coverage - band)  = $(coverage - band_lower)")
        println("  nominal 1 - alpha      = $(1.0 - P14_ALPHA_CONFORMAL)")
        println("  set sizes              singleton $(set_size_counts.singleton) / " *
                "ambiguous $(set_size_counts.ambiguous) / empty $(set_size_counts.empty)")
        println("  set-size rates         singleton $(diag.singleton_rate) / " *
                "ambiguous $(diag.ambiguous_rate) / empty $(diag.empty_rate)")
        println("  vacuous_hedge          = $(diag.vacuous_hedge)   " *
                "(LABEL, not a verdict; floor $P14_AMBIGUOUS_RATE_FLOOR)")
        println("  calibration masses     = $cal_masses")
        println("  evaluation masses      = $eval_masses")
        println("  expected class prior   = $(p14_pi_class_masses())")
        println("  conformal status x OOD state (rows $(P14_CONFORMAL_STATUSES), " *
                "cols $(P14_OOD_STATES)):")
        for (i, s) in enumerate(P14_CONFORMAL_STATUSES)
            println("    $s: ", cross_tab[i, :])
        end
        println("  guarantee basis        = simulator-derived held-out draws; the real-data " *
                "bound D-03 intended is ABSENT (D-03a), not merely loose")
        println("  iteration_trigger_fired = $iteration_trigger_fired")
        println("  elapsed                = $(round(time() - t_start; digits = 1)) s")
        println("-"^78)
        if !reported
            println("SMOKE MODE: the band is NOT asserted and these numbers are NOT a verdict.")
        end
    end

    # --- 9. THE ASSERTION, LAST ----------------------------------------------------------------
    if reported
        @assert coverage >= band_lower """
        SC1-d NOT MET. Realized split-conformal marginal coverage is $coverage on the reported
        evaluation set, BELOW the pre-registered lower band $band_lower (gap
        $(coverage - band_lower)). The band is DERIVED from P14_ALPHA_CONFORMAL and P14_N_EVAL by
        formula and is therefore not itself adjustable.

        AN HONEST SHORTFALL IS A PHASE-14 FINDING. It is REPORTED, not repaired. The ONLY
        authorised response is the single pre-declared iteration trigger, frozen in
        spike/p14/consts.jl before any result existed:

        $P14_ITERATION_TRIGGER

        i.e. switch the nonconformity score from LAC to APS and re-run ONCE, spending
        P14_ITERATION_ALLOWANCE. Relaxing P14_ALPHA_CONFORMAL, moving n_eval, widening the band,
        redrawing either set, reseeding or re-scoping the criterion are ALL forbidden, and
        spending the allowance is a pre-registration matter to put to the user rather than an
        executor's call.

        The full report was persisted BEFORE this assertion and is intact at:
          $report_path
          $pool_path
        """
    end

    return (qhat = qhat, coverage = coverage, band_lower = band_lower,
            iteration_trigger_fired = iteration_trigger_fired,
            cal_masses = cal_masses, eval_masses = eval_masses,
            set_size_counts = set_size_counts, vacuous_hedge = diag.vacuous_hedge,
            reported = reported, report_path = report_path, pool_path = pool_path)
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so including this file can NEVER trigger the reported run. Run it deliberately:
#
#     julia --project=spike -t auto spike/p14/run_p14_conformal.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
