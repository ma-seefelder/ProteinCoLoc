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

# spike/p14/run_p14_riskcoverage.jl --- GATED SC3-a / SC3-b / SC3-c: the risk-coverage curve of
# the scalar confidence score, its TWO reference curves computed on the SAME data, the selective
# skill normalized between them, the monotone-trend statistic on a PRE-REGISTERED coverage range,
# and the concentration of abstention on would-have-been-wrong items.
#
# =============================================================================================
# ANTI-SNOOPING BANNER. READ BEFORE CHANGING ANY NUMBER IN THIS FILE.
# =============================================================================================
# EVERY BAR THIS RUNNER SCORES AGAINST WAS COMMITTED TO `spike/p14/consts.jl` BEFORE ANY
# PHASE-14 RESULT EXISTED AND BEFORE THIS FILE EXISTED. The four are P14_SKILL_FLOOR,
# P14_SPEARMAN_FLOOR, P14_COVERAGE_FLOOR and P14_AUC_HARD_FLOOR; the evaluation size is
# P14_N_EVAL and the reserved stream counter is P14_EVAL_COUNTER. Every one of them is READ from
# that frozen file BY NAME. Not one of them is a literal here, and the acceptance criteria grep
# the comment-stripped source for exactly that, so a value cannot be "temporarily" inlined and
# left behind.
#
# ALL FOUR BARS AND THE COVERAGE FLOOR ARE JUDGEMENT CALLS WITH NO DERIVATION. `consts.jl`
# section E says so at length and `P14_JUDGEMENT_CALL_BARS` enumerates them as symbols so this
# runner can print and persist the label MECHANICALLY rather than by hand. They were put to the
# user as a BLOCKING checkpoint before `spike/p14/` existed and ruled accept-as-proposed, which
# makes them the user's bars rather than an agent's -- it does NOT make them derived.
#
# THIS PROJECT HAS RECORDED FOUR SEPARATE CASES OF AN UNDERIVED BAR MEASURING SOMETHING OTHER
# THAN WHAT IT NAMED (SC1g's component-vs-total, the n=2-against-271 real arm, the
# Wald-labelled-Wilson sizing, and the Phase-12 Stage-1 control ceiling). That is why the
# licensing standard for amending one of these is EVIDENCE, INDEPENDENT OF THE MACHINERY UNDER
# AUDIT, THAT THE BAR MEASURES SOMETHING OTHER THAN WHAT IT NAMES -- and why a shortfall here is
# REPORTED rather than repaired. `P14_ITERATION_ALLOWANCE = 1` cannot be spent on any of these:
# its ONE pre-declared trigger concerns split-conformal COVERAGE falling below the SC1-d band,
# and 14-09 measured that ABOVE the band, so the trigger did not fire.
#
# THE COVERAGE FLOOR IS THE MOST SNOOPABLE NUMBER IN THIS FILE AND IS TREATED AS SUCH. A floor
# chosen AFTER seeing where the empirical curve becomes monotone is invisible in the resulting
# statistic and total in the claim. `P14_COVERAGE_FLOOR` was frozen in commit
# a6c825867dbc786f7c3927df6760295ee0930c77 -- which contains no Phase-14 result of any kind --
# and the artifact records that sha next to `coverage_floor_frozen_before_run = true` and a
# `forbidden_action` string naming the two moves that are not available here.
#
# MONOTONICITY IS NOT ASSERTED POINTWISE, AND THAT IS DELIBERATE. The empirical risk-coverage
# curve is a step function whose denominator is the number of SELECTED items, so at low coverage
# a single hard item moves the selective risk by 1/k. SC3-b is therefore a Spearman statistic on
# a declared coverage range. The RAW curve over the FULL range is persisted and plotted anyway,
# including the noisy low-coverage part the statistic deliberately excludes -- shown, not hidden.
#
# A CURVE CANNOT BE DRAWN OVER THE D-05 FUSION, AND THIS RUNNER DOES NOT TRY. A risk-coverage
# sweep needs a single scalar to threshold. The fusion's inputs are a BINARY OOD flag, an INTEGER
# conformal set size and a cross-method record; sweeping "the fusion" would mean sweeping an
# object that has no total order, and the resulting picture would be meaningless. The sweep here
# is over kappa = max_y p-hat(y | Z) ONLY -- which is also the quantity the LAC hedge thresholds,
# so the curve and the conformal layer describe ONE ordering rather than two.
#
# THE ORDER OF OPERATIONS IS THE POINT: read the bars, compute the statistics, RENDER THE FIGURE,
# PERSIST THE ARTIFACT, print the headline, then assert. The figure and the artifact are written
# before anything can throw, so a FAILING gate still leaves a complete, self-describing report
# and a picture of the failure on disk. The evidence of an honest failure must survive the
# failure.
#
# THE EVALUATION SET IS CONSUMED, NEVER REDRAWN. `spike/p14/p14_eval_pool.jld2` was drawn ONCE by
# plan 14-09 at `P14_EVAL_COUNTER` and persisted precisely so that SC1-b, SC1-d and SC3-a/b/c are
# computed on the SAME 2000 items. If the file is missing this runner FAILS LOUDLY and names
# 14-09; it never redraws silently.
#
# NOTHING NUMERIC IS RE-IMPLEMENTED. The AUC is the in-repo tie-aware `roc_auc`
# (`spike/validation/ood.jl:319`) behind the same one-line wrapper Phase 13 used, and the
# Spearman is `StatsBase.corspearman` -- never hand-rolled (`spike/validation/run_p11_probe.jl:59`,
# PATTERNS S6). No dependency is added; StatsBase and CairoMakie are among the frozen sixteen.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY AND NOT A TEST. Including it does nothing; run it:
#
#     julia --project=spike -t auto spike/p14/run_p14_riskcoverage.jl
#
# DECOUPLING: spike-local. `src/` is reached only READ-ONLY and transitively, and step 0 PROVES
# at run time that it is byte-unchanged AND that no untracked file has appeared under it.

using JLD2
using Statistics
using Dates
using Random          # stdlib: `shuffle` for the RECORDED reference-curve permutation
import StatsBase      # corspearman, never hand-rolled

# --- Guarded includes, in dependency order --------------------------------------------------------
# `decide.jl` transitively pulls the pre-registration, the provenance loader, the posterior, the
# hedge, the fusion, the result type and the Phase-13 read surface. `ood.jl` is named explicitly
# because `roc_auc` is a DIRECT dependency of SC3-c and must not arrive by accident, and
# `figures.jl` is where the headless CairoMakie backend is activated (`spike/p13/run_three_way_gate.jl:85`).
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :P14_PROVENANCE_LOADED)   || include(joinpath(@__DIR__, "provenance.jl"))
isdefined(@__MODULE__, :P14_POSTERIOR_LOADED)    || include(joinpath(@__DIR__, "posterior.jl"))
isdefined(@__MODULE__, :P14_CONFORMAL_LOADED)    || include(joinpath(@__DIR__, "conformal.jl"))
isdefined(@__MODULE__, :P14_FUSE_LOADED)         || include(joinpath(@__DIR__, "fuse.jl"))
isdefined(@__MODULE__, :P14Result)               || include(joinpath(@__DIR__, "result.jl"))
isdefined(@__MODULE__, :P14_DECIDE_LOADED)       || include(joinpath(@__DIR__, "decide.jl"))
isdefined(@__MODULE__, :P14_POOLS_LOADED)        || include(joinpath(@__DIR__, "pools.jl"))
isdefined(@__MODULE__, :roc_auc) ||
    include(joinpath(@__DIR__, "..", "validation", "ood.jl"))
isdefined(@__MODULE__, :VAL_FIG_DIR) ||
    include(joinpath(@__DIR__, "..", "validation", "figures.jl"))

if !isdefined(@__MODULE__, :P14_RC_RUNNER_LOADED)
    "The SC3-a/b/c artifact. Written BEFORE any of the three assertions can throw."
    const P14_RC_REPORT_PATH = joinpath(@__DIR__, "p14_riskcoverage_report.jld2")

    """
    The SHARED evaluation pool, produced by plan 14-09.

    NAMED here rather than obtained by including `run_p14_conformal.jl`: that file defines its own
    `main`, and including one runner from another would put two zero-argument `main` methods in one
    namespace where the last include silently wins. The two spellings of the path are pinned to
    each other by the pool's own recorded `P14_EVAL_COUNTER` and `n_eval`, which this runner
    asserts, so a path pointing at some other file could not pass.
    """
    const P14_RC_EVAL_POOL_PATH = joinpath(@__DIR__, "p14_eval_pool.jld2")

    "The risk-coverage figure. Rendered BEFORE the gate; a gitignored regenerable artifact."
    const P14_RC_FIG_PATH = normpath(joinpath(@__DIR__, "..", "figures", "p14_riskcoverage.png"))

    "The columns this runner requires from the shared pool, so a schema drift fails by name."
    const P14_RC_REQUIRED_POOL_KEYS = ("n", "n_eval", "P14_EVAL_COUNTER", "true_class",
                                       "class_posterior_coloc", "class_posterior_random",
                                       "class_posterior_exclusion", "confidence",
                                       "conformal_status", "ood_state", "decision",
                                       "abstain_reason", "qhat", "class_prior_used",
                                       "eval_masses", "idx")

    """
    The tolerance at which the recomputed confidence must equal the column 14-09 persisted.

    Not a threshold on a result: a reproduction check. `kappa` is recomputed here BY NAME from the
    three stored class posteriors rather than read off the `confidence` column, because reading the
    column would make this runner trust a summary instead of the quantity it names; the two are
    then pinned to each other, so a drift in either is a named failure rather than a silent one.
    """
    const P14_RC_CONF_TOL = 1e-12

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else can
    # make this body skip.
    const P14_RC_RUNNER_LOADED = true
end

# =============================================================================================
# Small local helpers
# =============================================================================================

"""
    _p14_rc_save_report(path, required; kwargs...) -> String

Atomically persist an artifact: write `path * ".tmp"`, REOPEN it read-only and integrity-check
every key in `required`, then `mv(...; force = true)`.

The `_p13_gate_save_report` idiom (`spike/p13/run_three_way_gate.jl:171-190`), carried per runner
exactly as Phase 13 and `run_p14_fdr_check.jl` carry it. It is deliberately NOT imported from a
sibling runner: importing would mean including that runner and inheriting its `main`.

A crash mid-write leaves a discardable `.tmp` rather than a torn artifact a later reader would
happily believe.
"""
function _p14_rc_save_report(path, required; kwargs...)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        for k in required
            @assert haskey(f, k) "_p14_rc_save_report: integrity check failed -- $tmp is missing the required key `$k`"
        end
    end
    mv(tmp, path; force = true)
    return path
end

"""
    p14_rc_auc(neg_scores, pos_scores) -> Float64

The AUC only, from the in-repo hand-rolled tie-aware `roc_auc` (`spike/validation/ood.jl:319`).

Copied verbatim in shape from `p13_gate_auc` (`spike/p13/run_three_way_gate.jl:215`) rather than
re-indexing `[3]` at every call site, and rather than adding a ROC package. `roc_auc` returns the
Mann-Whitney U statistic exactly, which is what SC3-c wants: `P(score_pos > score_neg) + 0.5 *
P(=)`, robust to the threshold degeneracies a confident head produces.
"""
p14_rc_auc(neg_scores, pos_scores) = roc_auc(neg_scores, pos_scores)[3]

"""
    _p14_rc_class_from_key(s::Symbol) -> ThreeWayClass

The inverse of the class-symbol map the shared pool was written through, by an explicit
three-branch map.

Written out rather than derived from the enum's spelling because five class orderings coexist in
this codebase and every statistic here is invariant to a CONSISTENT relabelling -- which is
exactly what makes an inconsistent one produce numbers that look right. A fall-through would
silently drop a class out of the mass guard below, so the map is total and errors otherwise.
"""
function _p14_rc_class_from_key(s::Symbol)
    s === :coloc     && return COLOC
    s === :random    && return RANDOM
    s === :exclusion && return EXCLUSION
    error("_p14_rc_class_from_key: unhandled class symbol `$s`. The inverse class map must be " *
          "total; a fall-through would drop a class out of the mass re-derivation.")
end

"""
    _p14_rc_argmax(p::NamedTuple) -> Symbol

The argmax class of a three-way posterior, resolved BY NAME.

NEVER `argmax` over a positional vector: `P14_CLASS_KEYS` exists because five different class
orderings coexist here, and a positional argmax would silently inherit whichever one the caller
happened to build. Ties resolve in the frozen `P14_CLASS_KEYS` order, stated so the rule is
recorded rather than emergent; with a continuous head an exact three-way tie does not occur, and
if one did, any fixed rule is defensible so long as it is written down.
"""
function _p14_rc_argmax(p::NamedTuple)
    best_k = first(P14_CLASS_KEYS)
    best_v = getproperty(p, best_k)
    for k in P14_CLASS_KEYS
        v = getproperty(p, k)
        if v > best_v
            best_v = v
            best_k = k
        end
    end
    return best_k
end

"""
    _p14_rc_curve(kappa, loss) -> NamedTuple

The MEASURED risk-coverage curve, swept over the DESCENDING UNIQUE values of the confidence score.

The threshold convention is `roc_auc`'s (`spike/validation/ood.jl:325`, and its shipped twin at
`src/amortized/ood.jl:281`): descending unique attained values, `>=` at each. With a continuous
head ties are rare, but a sweep must define a rule and the rule must be the one already used
elsewhere in this repository, or two curves in one report would be swept differently.

At threshold `t`: `coverage = mean(kappa .>= t)` and `selective_risk = (number of selected items
whose argmax call is wrong) / (number selected)`. The FULL range is returned -- the first point
sits at coverage 1/n, where the denominator is one item and the risk is 0 or 1. That noise is the
point of `P14_COVERAGE_FLOOR`, and it is REPORTED rather than trimmed.
"""
function _p14_rc_curve(kappa::AbstractVector{<:Real}, loss::AbstractVector{<:Integer})
    n = length(kappa)
    n == length(loss) || throw(ArgumentError(
        "_p14_rc_curve: kappa holds $n entries and loss holds $(length(loss))"))
    thr = sort(unique(kappa); rev = true)
    cov  = Float64[]
    risk = Float64[]
    for t in thr
        sel   = 0
        nloss = 0
        for i in 1:n
            if kappa[i] >= t
                sel   += 1
                nloss += loss[i]
            end
        end
        sel == 0 && continue
        push!(cov, sel / n)
        push!(risk, nloss / sel)
    end
    return (coverage = cov, selective_risk = risk, thresholds = collect(Float64, thr))
end

"""
    _p14_rc_order_curve(loss, ord) -> NamedTuple

A reference risk-coverage curve from an EXPLICIT ORDERING of the items, evaluated at every
coverage `k/n`.

Both reference curves are orderings rather than scores, so the descending-unique-threshold sweep
does not apply to them: the oracle's natural score takes two values only, and a two-point curve
would under-resolve the very envelope it exists to draw. Evaluating at `k = 1 ... n` gives the
finest grid the data admits and spans the same `[1/n, 1]` coverage range the measured curve does,
which is what makes the three integrals comparable. The runner asserts that span agreement rather
than assuming it.
"""
function _p14_rc_order_curve(loss::AbstractVector{<:Integer}, ord::AbstractVector{<:Integer})
    n = length(ord)
    n == length(loss) || throw(ArgumentError(
        "_p14_rc_order_curve: the ordering holds $n entries and loss holds $(length(loss))"))
    cov  = Vector{Float64}(undef, n)
    risk = Vector{Float64}(undef, n)
    c = 0
    for k in 1:n
        c += loss[ord[k]]
        cov[k]  = k / n
        risk[k] = c / k
    end
    return (coverage = cov, selective_risk = risk)
end

"""
    _p14_rc_aurc(coverage, selective_risk) -> Float64

The area under the risk-coverage curve, by the trapezoidal rule against coverage.

Left as a RAW integral over the realized coverage span rather than divided by that span: the
selective-skill statistic is a ratio of DIFFERENCES of three such integrals over the same span, so
any common normalization cancels exactly. The span itself is recorded in the artifact so a reader
can convert to a mean selective risk without re-running anything.
"""
function _p14_rc_aurc(coverage::AbstractVector{<:Real}, selective_risk::AbstractVector{<:Real})
    length(coverage) == length(selective_risk) || throw(ArgumentError(
        "_p14_rc_aurc: $(length(coverage)) coverage points against $(length(selective_risk)) risks"))
    length(coverage) >= 2 || return NaN
    s = 0.0
    for i in 2:length(coverage)
        s += 0.5 * (selective_risk[i] + selective_risk[i - 1]) * (coverage[i] - coverage[i - 1])
    end
    return s
end

"""
    _p14_rc_load_pool(path; reported = true) -> NamedTuple

Load the SHARED evaluation pool written by plan 14-09 and prove it is the pool this runner is
entitled to score, rather than merely a file with the right name.

Four checks, each catching something the others cannot: the file EXISTS (and if not, this fails
loudly naming 14-09 and never redraws); every column read here is PRESENT by name; the recorded
`P14_EVAL_COUNTER` and `n_eval` are the FROZEN ones; and the class masses are RE-DERIVED from the
stored labels through the frozen `p14_assert_unstratified` guard rather than trusted from the
stored summary. The guard runs only on a reported pool -- 14-09 recorded that below roughly
n = 400 the equal-thirds detector has no power, so a smoke pool would trip it as a small-sample
artifact rather than as a finding.
"""
function _p14_rc_load_pool(path::AbstractString; reported::Bool = true)
    isfile(path) || error("""
        run_p14_riskcoverage: the SHARED evaluation pool is ABSENT at
            $path

        That pool is produced by PLAN 14-09 (`spike/p14/run_p14_conformal.jl`), which draws it
        ONCE at the reserved P14_EVAL_COUNTER and persists it so that SC1-b, SC1-d and SC3-a/b/c
        are scored on the SAME items. Run it first:

            julia --project=spike -t auto spike/p14/run_p14_conformal.jl

        THIS RUNNER WILL NOT REDRAW THE EVALUATION SET, silently or otherwise. A redraw would put
        the risk-coverage curve on a different sample from the coverage number it is quoted
        beside, and nothing in either artifact would say so.""")

    d = JLD2.load(path)
    for k in P14_RC_REQUIRED_POOL_KEYS
        haskey(d, k) || error(
            "run_p14_riskcoverage: the shared evaluation pool at $path carries no `$k` column. " *
            "The 14-09 pool schema moved under this runner; scoring SC3 against a partially " *
            "understood artifact is how a number ends up meaning something other than its name.")
    end

    n = Int(d["n"])
    @assert n == Int(d["n_eval"]) "run_p14_riskcoverage: the pool stores $(n) items but records n_eval = $(d["n_eval"]); the artifact disagrees with itself and neither number can be quoted"
    @assert Int(d["P14_EVAL_COUNTER"]) == Int(P14_EVAL_COUNTER) "run_p14_riskcoverage: the pool was drawn at counter $(d["P14_EVAL_COUNTER"]) but the frozen reserved counter is $(P14_EVAL_COUNTER); a pool from another stream is a different set of items and may not be scored as the reported evaluation set"
    if reported
        @assert n == Int(P14_N_EVAL) "run_p14_riskcoverage: the pool holds $n items but the frozen reported evaluation size is $(P14_N_EVAL)"
    end

    true_class = Vector{Symbol}(d["true_class"])
    @assert length(true_class) == n "run_p14_riskcoverage: the pool's true_class column holds $(length(true_class)) entries for $n items"

    realized = if reported
        p14_assert_unstratified([_p14_rc_class_from_key(s) for s in true_class])
    else
        class_masses([_p14_rc_class_from_key(s) for s in true_class])
    end
    stored = d["eval_masses"]
    for k in P14_MASS_KEYS
        @assert getproperty(realized, k) == getproperty(stored, k) "run_p14_riskcoverage: the `$k` mass re-derived from the pool's own labels is $(getproperty(realized, k)) while the pool records $(getproperty(stored, k)); the stored summary and the stored labels disagree"
    end

    return (n = n,
            path = path,
            true_class = true_class,
            p_coloc = Vector{Float64}(d["class_posterior_coloc"]),
            p_random = Vector{Float64}(d["class_posterior_random"]),
            p_exclusion = Vector{Float64}(d["class_posterior_exclusion"]),
            confidence_stored = Vector{Float64}(d["confidence"]),
            conformal_status = Vector{Symbol}(d["conformal_status"]),
            ood_state = Vector{Symbol}(d["ood_state"]),
            decision = Vector{Symbol}(d["decision"]),
            abstain_reason = Vector{Symbol}(d["abstain_reason"]),
            idx = Vector{Int}(d["idx"]),
            qhat = Float64(d["qhat"]),
            prior = d["class_prior_used"],
            eval_masses = realized,
            pool_provenance = get(d, "pool_provenance", nothing),
            pool_generated = get(d, "generated", ""),
            pool_reported = get(d, "reported", false))
end

"""
    p14_rc_figure(curve, oracle, rnd, base_error; path) -> String

The measured risk-coverage curve with BOTH reference curves, over the FULL coverage range, the
pre-registered statistic range shaded, and the base error rate annotated.

A SCRIPT ARTIFACT AND NEVER A GATE ASSERTION, and rendered BEFORE the gate for the same reason the
artifact is persisted before it: a failing assertion throws and everything after it is never
reached, so a figure emitted afterwards would exist only for a PASSING run. The evidence of an
honest failure must survive the failure.

The low-coverage part is DRAWN, not trimmed. A reader must be able to see the region the SC3-b
statistic excludes and judge the exclusion for themselves; a plot that started at the floor would
be an argument dressed as a picture.
"""
function p14_rc_figure(curve, oracle, rnd, base_error; path = P14_RC_FIG_PATH)
    mkpath(dirname(path))
    fig = Figure(size = (960, 560))
    ax = CairoMakie.Axis(fig[1, 1];
        title = "risk-coverage over kappa = max_y p-hat(y | Z) -- RAW, FULL RANGE",
        subtitle = "shaded = the PRE-REGISTERED SC3-b statistic range [P14_COVERAGE_FLOOR, 1]; " *
                   "floor frozen in a6c8258, before this runner existed",
        xlabel = "coverage (fraction of the batch retained)",
        ylabel = "selective risk (error rate among retained items)")
    vspan!(ax, P14_COVERAGE_FLOOR, 1.0; color = (:grey, 0.16),
           label = "SC3-b statistic range")
    hlines!(ax, [base_error]; color = :black, linestyle = :dot,
            label = "base error rate = $(round(base_error; digits = 4))")
    lines!(ax, rnd.coverage, rnd.selective_risk; color = :darkorange,
           label = "random reference (recorded shuffle)")
    lines!(ax, oracle.coverage, oracle.selective_risk; color = :seagreen,
           label = "oracle reference (lower envelope)")
    lines!(ax, curve.coverage, curve.selective_risk; color = :steelblue, linewidth = 2,
           label = "MEASURED")
    vlines!(ax, [P14_COVERAGE_FLOOR]; color = :grey, linestyle = :dash)
    axislegend(ax; position = :lt, framevisible = true)
    save(path, fig)
    return path
end

# =============================================================================================
# The reported run
# =============================================================================================

"""
    main(; pool_path = P14_RC_EVAL_POOL_PATH, report_path = P14_RC_REPORT_PATH,
           fig_path = P14_RC_FIG_PATH, verbose = true) -> NamedTuple

The GATED SC3-a / SC3-b / SC3-c measurement on the SHARED 14-09 evaluation pool.

`pool_path` defaults to that pool. Passing anything else puts the run in SMOKE mode: it writes to
`_smoke` artifact paths and prints a banner saying so, so a load check can never be mistaken for,
or quoted as, a verdict, and can never overwrite the reported artifact.
"""
function main(; pool_path = P14_RC_EVAL_POOL_PATH,
                report_path = P14_RC_REPORT_PATH,
                fig_path = P14_RC_FIG_PATH,
                verbose::Bool = true)

    t_start  = time()
    reported = (pool_path == P14_RC_EVAL_POOL_PATH)
    if !reported
        report_path = replace(report_path, ".jld2" => "_smoke.jld2")
        fig_path    = replace(fig_path, ".png" => "_smoke.png")
    end

    # --- 0. THE DECOUPLING PROOF, AT RUN TIME (D-01 / D-02) ---------------------------------
    # Asserted WHILE the run happens, not only in review: a run that started on a clean tree and
    # finished on a dirty one would be a repudiation hole (T-14-11). The shared lane guard is
    # called first; the two commands are then re-run here so this runner's own failure message
    # names the breach, and a THIRD check is added that neither the guard nor a diff can make --
    # an untracked NEW file under `src/` is invisible to `git diff HEAD` and would be exactly the
    # productionization-by-stealth D-01 forbids.
    lane = p14_assert_lane_clean(; verbose = verbose)
    if lane.checked
        src_clean = success(Cmd(`git diff --quiet HEAD -- src`; dir = P14_REPO_ROOT))
        env_clean = success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                                dir = P14_REPO_ROOT))
        src_untracked = readchomp(Cmd(`git status --porcelain -- src`; dir = P14_REPO_ROOT))
        @assert src_clean "D-01 decoupling breach: `git diff --quiet HEAD -- src` FAILED -- src/ is not byte-unchanged, so this run may not be reported"
        @assert env_clean "D-02 decoupling breach: spike/Project.toml or spike/Manifest.toml is modified -- the frozen spike environment moved under the run, and this runner needs no dependency at all"
        @assert isempty(src_untracked) "D-01 decoupling breach: an UNTRACKED file appeared under src/ ($(src_untracked)) -- a diff against HEAD cannot see it, and shipping the decision layer by stealth is precisely what D-01 forbids"
    end

    # --- 1. THE INHERITANCE (D-07) AND THE SHARED POOL --------------------------------------
    bundle = p14_load_bundle()
    pool   = _p14_rc_load_pool(pool_path; reported = reported)
    n_eval = pool.n
    for k in P14_CLASS_KEYS
        @assert getproperty(pool.prior, k) == getproperty(bundle.prior, k) "run_p14_riskcoverage: the pool records a `$k` prior mass of $(getproperty(pool.prior, k)) while the inherited bundle measures $(getproperty(bundle.prior, k)); the pool was scored under a different prior from the one this runner inherits"
    end

    # --- 2. THE LOCKED-THRESHOLD BANNER, BEFORE ANYTHING IS COMPUTED -------------------------
    if verbose
        println("="^78)
        println("PHASE-14 GATED SC3-a / SC3-b / SC3-c: SELECTIVE RISK-COVERAGE")
        println("="^78)
        reported || println("*** SMOKE MODE -- NOT A REPORTED RUN, THE GATES ARE NOT ASSERTED ***")
        println("AMENDED CRITERIA: ", P14_AMENDMENT_NOTICE)
        println("GATED BARS -- EVERY ONE A JUDGEMENT CALL WITH NO DERIVATION:")
        println("  P14_SKILL_FLOOR      = $P14_SKILL_FLOOR   (consts.jl section E; SC3-a)")
        println("                         GATED. JUDGEMENT CALL -- NO DERIVATION.")
        println("  P14_SPEARMAN_FLOOR   = $P14_SPEARMAN_FLOOR   (consts.jl section E; SC3-b)")
        println("                         GATED. JUDGEMENT CALL -- NO DERIVATION.")
        println("  P14_COVERAGE_FLOOR   = $P14_COVERAGE_FLOOR   (consts.jl section E; SC3-b range)")
        println("                         GATED. JUDGEMENT CALL -- NO DERIVATION.")
        println("  P14_AUC_HARD_FLOOR   = $P14_AUC_HARD_FLOOR   (consts.jl section E; SC3-c)")
        println("                         GATED. JUDGEMENT CALL -- NO DERIVATION.")
        println("  the judgement-call enumeration, printed MECHANICALLY from the frozen file so a")
        println("  bar cannot quietly stop being labelled one:")
        println("  P14_JUDGEMENT_CALL_BARS = $P14_JUDGEMENT_CALL_BARS")
        println("  (the fifth, P14_OOD_MARGIN_FLOOR, is SC3-d's and belongs to another runner.)")
        println("FREEZE PROVENANCE:")
        println("  all four were committed in a6c825867dbc786f7c3927df6760295ee0930c77, a commit")
        println("  containing no Phase-14 result of any kind. The coverage floor specifically was")
        println("  NOT chosen after seeing where the curve becomes monotone, and re-choosing it")
        println("  now -- or monotonizing the curve -- is the failure mode this project has paid")
        println("  for repeatedly and is not available here.")
        println("WHAT IS SWEPT (and what CANNOT be):")
        println("  kappa = max_y p-hat(y | Z), the SAME scalar the LAC hedge thresholds, so the")
        println("  curve and the conformal layer describe ONE ordering. A curve cannot be drawn")
        println("  over the D-05 fusion: an OOD flag is binary and a set size is an integer, so")
        println("  the fusion has no total order to sweep.")
        println("N_EVAL / STREAM:")
        println("  P14_N_EVAL           = $P14_N_EVAL           (pool holds $n_eval)")
        println("  P14_EVAL_COUNTER     = $P14_EVAL_COUNTER               (the SHARED reported evaluation set)")
        println("ITERATION DISCIPLINE:")
        println("  P14_ITERATION_ALLOWANCE = $P14_ITERATION_ALLOWANCE, and it CANNOT be spent here:")
        println(P14_ITERATION_TRIGGER)
        println("INHERITED (D-07, LOADED and four-way agreed, never re-derived):")
        println("  tau                     = $(bundle.tau)")
        println("  class prior (MEASURED)  = $(bundle.prior)")
        println("SHARED POOL (CONSUMED, never redrawn):")
        println("  path                    = $(pool.path)")
        println("  q-hat the sets were cut at = $(pool.qhat)")
        println("  realized class masses   = $(pool.eval_masses)")
        println("  generated               = $(pool.pool_generated)")
        println("="^78)
    end

    # --- 3. THE CONFIDENCE SCORE AND THE LOSS ------------------------------------------------
    # kappa is RECOMPUTED by name from the three stored class posteriors and then pinned to the
    # column 14-09 wrote, rather than read off that column. `p14_confidence` is the same function
    # the LAC hedge and the fusion consult, so the curve below is swept over the SAME ordering the
    # conformal layer cuts -- not a second, incidentally-similar score.
    verbose && println("[1/6] recomputing kappa by name and pinning it to the persisted column …")
    posteriors = [(coloc = pool.p_coloc[i], random = pool.p_random[i],
                   exclusion = pool.p_exclusion[i]) for i in 1:n_eval]
    kappa = [p14_confidence(p) for p in posteriors]
    dev_conf = maximum(abs.(kappa .- pool.confidence_stored))
    @assert dev_conf < P14_RC_CONF_TOL "run_p14_riskcoverage: the recomputed confidence deviates from the pool's stored column by up to $dev_conf (tolerance $(P14_RC_CONF_TOL)); the curve would then be swept over a different ordering from the one the hedge cut"

    # THE LOSS IS THE ARGMAX CALL'S CORRECTNESS AGAINST THE SIMULATOR'S TRUE LABEL. The simulator
    # knows the class of every item, so this is a MEASUREMENT and not an estimate -- and it is the
    # same vector SC3-c calls `hard`, which is why there is exactly one of it: `hard_i` is
    # precisely "the argmax class would have been wrong", so defining it twice could only create a
    # way for the two to disagree.
    argmax_class = [_p14_rc_argmax(p) for p in posteriors]
    loss = [argmax_class[i] === pool.true_class[i] ? 0 : 1 for i in 1:n_eval]
    base_error = mean(loss)
    verbose && println("      max |kappa_recomputed - kappa_stored| = $dev_conf   " *
                       "base error rate = $base_error")

    # --- 4. THE MEASURED CURVE, OVER THE FULL RANGE ------------------------------------------
    verbose && println("[2/6] sweeping the descending unique kappa thresholds …")
    curve = _p14_rc_curve(kappa, loss)
    @assert last(curve.coverage) == 1.0 "run_p14_riskcoverage: the measured curve does not reach full coverage; the descending sweep did not include the smallest attained kappa"
    @assert abs(last(curve.selective_risk) - base_error) < P14_RC_CONF_TOL "run_p14_riskcoverage: at full coverage the selective risk must BE the base error rate; got $(last(curve.selective_risk)) against $base_error"
    verbose && println("      $(length(curve.coverage)) curve points, coverage in " *
                       "[$(minimum(curve.coverage)), $(maximum(curve.coverage))]")

    # --- 5. THE TWO REFERENCE CURVES, ON THE SAME DATA ---------------------------------------
    # Normalizing between BOTH references is what makes SC3-a immune to the base error rate: a
    # near-perfect classifier and a mediocre one produce very different raw AURCs, and a bar on the
    # raw number would mostly be measuring how easy the batch was.
    verbose && println("[3/6] building the oracle and random reference curves …")
    oracle_ord = sortperm(loss)                       # stable: every correct item first
    oracle = _p14_rc_order_curve(loss, oracle_ord)
    # A RECORDED shuffle, not an unseeded one. The permutation is a reproducible function of the
    # frozen stream, so the random reference curve is regenerable bitwise; the seed and counter are
    # persisted beside it.
    rng = p14_rng(P14_EVAL_COUNTER)
    random_ord = shuffle(rng, collect(1:n_eval))
    rnd = _p14_rc_order_curve(loss, random_ord)
    @assert sort(random_ord) == collect(1:n_eval) "run_p14_riskcoverage: the recorded shuffle is not a permutation of the batch"
    @assert abs(last(oracle.selective_risk) - base_error) < P14_RC_CONF_TOL "run_p14_riskcoverage: the oracle curve must also end at the base error rate"
    @assert abs(last(rnd.selective_risk) - base_error) < P14_RC_CONF_TOL "run_p14_riskcoverage: the random curve must also end at the base error rate"

    # The three integrals are only comparable over a common coverage span. Asserted rather than
    # assumed: a tie at the top of the kappa ordering would start the measured curve above 1/n, and
    # a difference of that order would leak straight into the skill ratio.
    span_tol = 5.0 / n_eval
    cov_min_measured = minimum(curve.coverage)
    cov_min_ref = minimum(oracle.coverage)
    @assert abs(cov_min_measured - cov_min_ref) <= span_tol "run_p14_riskcoverage: the measured curve starts at coverage $cov_min_measured while the reference curves start at $cov_min_ref; the three AURCs would be integrals over different intervals and their difference would not be a skill"

    verbose && println("[4/6] integrating the three curves (trapezoidal, against coverage) …")
    aurc        = _p14_rc_aurc(curve.coverage, curve.selective_risk)
    aurc_oracle = _p14_rc_aurc(oracle.coverage, oracle.selective_risk)
    aurc_random = _p14_rc_aurc(rnd.coverage, rnd.selective_risk)
    e_aurc      = aurc - aurc_oracle
    oracle_is_lower_envelope = aurc_oracle <= aurc
    verbose && println("      AURC = $aurc   oracle = $aurc_oracle   random = $aurc_random   " *
                       "E-AURC = $e_aurc")

    # --- 6. THE FIGURE, RENDERED BEFORE ANY GATE IS EVALUATED --------------------------------
    verbose && println("[5/6] rendering the risk-coverage figure (BEFORE any verdict) …")
    figout = p14_rc_figure(curve, oracle, rnd, base_error; path = fig_path)
    verbose && println("      figure -> $figout")

    # --- 7. PERSIST THE REPORT, BEFORE THE HEADLINE AND BEFORE ANY ASSERTION -----------------
    verbose && println("[6/6] persisting the risk-coverage report (BEFORE any verdict) …")
    consts_path = joinpath(@__DIR__, "consts.jl")
    _p14_rc_save_report(report_path,
        ("curve_coverage", "curve_selective_risk", "oracle_coverage", "oracle_selective_risk",
         "random_coverage", "random_selective_risk", "aurc", "aurc_oracle", "aurc_random",
         "e_aurc", "base_error_rate", "n_eval", "shuffle_seed", "shuffle_counter", "provenance");
        # --- the RAW measured curve, over the FULL coverage range ---
        curve_coverage = curve.coverage,
        curve_selective_risk = curve.selective_risk,
        curve_thresholds = curve.thresholds,
        curve_is_raw_full_range = true,
        curve_n_points = length(curve.coverage),
        # --- the two reference curves, computed on the SAME data ---
        oracle_coverage = oracle.coverage,
        oracle_selective_risk = oracle.selective_risk,
        random_coverage = rnd.coverage,
        random_selective_risk = rnd.selective_risk,
        random_permutation = random_ord,
        # --- the integrals ---
        aurc = aurc, aurc_oracle = aurc_oracle, aurc_random = aurc_random, e_aurc = e_aurc,
        aurc_coverage_span_lo = cov_min_measured,
        aurc_coverage_span_hi = maximum(curve.coverage),
        aurc_reference_span_lo = cov_min_ref,
        aurc_span_tolerance = span_tol,
        oracle_is_lower_envelope = oracle_is_lower_envelope,
        # --- the loss the whole object rests on ---
        base_error_rate = base_error,
        n_errors = sum(loss),
        loss = loss,
        confidence = kappa,
        argmax_class = argmax_class,
        true_class = pool.true_class,
        max_abs_dev_confidence = dev_conf,
        confidence_tolerance = P14_RC_CONF_TOL,
        score_note =
            "The sweep is over kappa = max_y p-hat(y | Z) ONLY -- the SAME scalar the LAC hedge " *
            "thresholds -- so the curve and the conformal layer describe ONE ordering. A " *
            "risk-coverage curve CANNOT be drawn over the D-05 fusion: an OOD flag is binary and " *
            "a conformal set size is an integer, so the fusion admits no total order to sweep, " *
            "and the two objects must not be conflated.",
        # --- EVERY BAR THIS RUN WILL BE SCORED AGAINST, written INTO the artifact ---
        P14_SKILL_FLOOR = P14_SKILL_FLOOR,
        P14_SPEARMAN_FLOOR = P14_SPEARMAN_FLOOR,
        P14_COVERAGE_FLOOR = P14_COVERAGE_FLOOR,
        P14_AUC_HARD_FLOOR = P14_AUC_HARD_FLOOR,
        P14_JUDGEMENT_CALL_BARS = collect(P14_JUDGEMENT_CALL_BARS),
        P14_N_EVAL = P14_N_EVAL,
        P14_EVAL_COUNTER = Int(P14_EVAL_COUNTER),
        P14_ITERATION_ALLOWANCE = P14_ITERATION_ALLOWANCE,
        P14_ITERATION_TRIGGER = P14_ITERATION_TRIGGER,
        iteration_trigger_fired = false,
        iteration_allowance_applies_here = false,
        # --- the recorded shuffle, so the random reference is regenerable bitwise ---
        shuffle_seed = UInt64(P14_DEV_SEED),
        shuffle_salt = UInt64(P14_SALT),
        shuffle_counter = Int(P14_EVAL_COUNTER),
        shuffle_recipe = "shuffle(p14_rng(P14_EVAL_COUNTER), collect(1:n_eval))",
        # --- what the numbers actually rest on ---
        guarantee_basis = :simulator_derived_held_out_draws,
        named_limits = p14_named_limits(),
        named_limits_count = length(p14_named_limits()),
        amendment = P14_AMENDMENT_NOTICE,
        # --- provenance ---
        n_eval = n_eval,
        provenance = p14_provenance_record(bundle.prov),
        eval_pool_path = pool.path,
        eval_pool_sha = p14_blob_sha(pool.path),
        eval_pool_provenance = pool.pool_provenance,
        eval_pool_generated = pool.pool_generated,
        eval_pool_reported = pool.pool_reported,
        eval_masses = pool.eval_masses,
        qhat = pool.qhat,
        class_prior_used = bundle.prior,
        tau = bundle.tau,
        grid = bundle.grid,
        figure_path = figout,
        consts_sha = p14_consts_sha(),
        consts_git_blob_sha = p14_blob_sha(consts_path),
        reported = reported,
        julia_version = string(VERSION),
        nthreads = Threads.nthreads(),
        elapsed_s = time() - t_start,
        generated = string(Dates.now(Dates.UTC)) * "Z")
    verbose && println("      persisted (before any verdict) -> $report_path")

    # --- 8. THE HEADLINE, PRINTED BEFORE ANY ASSERTION CAN THROW -----------------------------
    if verbose
        println("\n", "-"^78)
        println("RISK-COVERAGE, RAW AND OVER THE FULL RANGE")
        println("  n_eval                 = $n_eval")
        println("  base error rate        = $base_error   ($(sum(loss)) of $n_eval argmax calls wrong)")
        println("  curve points           = $(length(curve.coverage))")
        println("  coverage range         = [$(minimum(curve.coverage)), $(maximum(curve.coverage))]")
        println("  AURC (measured)        = $aurc")
        println("  AURC (oracle)          = $aurc_oracle")
        println("  AURC (random)          = $aurc_random")
        println("  E-AURC                 = $e_aurc")
        println("  oracle below measured  = $oracle_is_lower_envelope")
        println("  figure                 = $figout")
        println("  elapsed                = $(round(time() - t_start; digits = 1)) s")
        println("-"^78)
        reported || println("SMOKE MODE: no gate was asserted and these numbers are NOT a verdict.")
    end

    return (curve = curve, oracle = oracle, random = rnd,
            aurc = aurc, aurc_oracle = aurc_oracle, aurc_random = aurc_random,
            e_aurc = e_aurc, base_error_rate = base_error,
            n_eval = n_eval, reported = reported,
            report_path = report_path, figure_path = figout)
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so including this file can NEVER trigger the reported run. Run it deliberately:
#
#     julia --project=spike -t auto spike/p14/run_p14_riskcoverage.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
