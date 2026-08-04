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

# spike/p14/run_p14_fdr_check.jl --- REPORTED SC1-b: the REALIZED false discovery proportion among
# decided-and-accepted items, against the rule's OWN predicted posterior expected FDP, at four
# pre-registered levels -- plus the REPORTED-NOT-GATED prior-sensitivity sweep (SC1-c).
#
# =============================================================================================
# ANTI-SNOOPING BANNER. READ BEFORE CHANGING ANY NUMBER IN THIS FILE.
# =============================================================================================
# EVERY THRESHOLD THIS RUNNER SCORES AGAINST WAS COMMITTED TO `spike/p14/consts.jl` BEFORE ANY
# PHASE-14 RESULT EXISTED. The four-level sweep is `P14_ALPHA_FDR_GRID`, the five-rung prior
# sweep is `P14_PI_COLOC_GRID`, the evaluation size is `P14_N_EVAL` and the reserved stream
# counter is `P14_EVAL_COUNTER`. Every one of them is READ from that frozen file BY NAME. Not one
# of them is a literal here, and the acceptance criteria grep for exactly that, so a value cannot
# be "temporarily" inlined and left behind. The comparison band is not a constant either: it is
# DERIVED from the rule's own predicted value and the realized number of discoveries by the
# binomial standard-error formula, which is what makes it un-choosable.
#
# ROADMAP SC1 IS AMENDED BY D-02 AND SC2 BY D-05, both frozen in `14-CONTEXT.md` while NO
# Phase-14 result of any kind existed. ANY REPORT THAT CITES THIS RESULT MUST CITE THE ORIGINAL
# ROADMAP CRITERIA ALONGSIDE IT. The amendment notice travels INSIDE the artifact so a reader who
# never opens a planning document still meets it.
#
# AN HONEST SHORTFALL IS A PHASE-14 FINDING, NEVER A LICENCE TO ACT. A realized FDP above the
# band is REPORTED. It is not a reason to move a level, to widen the band, to redraw or restratify
# the evaluation set, to reseed, or to re-scope the criterion to a friendlier subset.
# `P14_ITERATION_ALLOWANCE = 1` cannot be spent here at all: its ONE pre-declared trigger concerns
# split-conformal COVERAGE falling below the SC1-d band, and 14-09 measured that at 0.9115 against
# 0.8799 -- so the trigger did not fire and the allowance is not available to this runner under any
# reading. Phase 7 was amended TWICE after seeing results and the credibility cost of that is why
# this paragraph is here. Phase 11 and Phase 12 each missed a pre-registered bar and REPORTED it;
# that is the precedent this runner follows.
#
# THE ORDER OF OPERATIONS IS THE POINT: read the thresholds, compute the statistics, PERSIST THE
# ARTIFACT, print the headline, then assert. The artifact is written before anything can throw, so
# a FAILING gate still leaves a complete, self-describing report on disk.
#
# THE EVALUATION SET IS CONSUMED, NEVER REDRAWN. `spike/p14/p14_eval_pool.jld2` was drawn ONCE by
# plan 14-09 at `P14_EVAL_COUNTER` and persisted precisely so that SC1-b, SC1-d and SC3-a/b/c are
# computed on the SAME 2000 items rather than on three independent draws at the same counter --
# which would make a coverage number and the FDR number quoted beside it two different experiments.
# Reusing the persisted pool also costs roughly four minutes less than redrawing it. If the file is
# missing this runner FAILS LOUDLY and names 14-09; it never redraws silently.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY AND NOT A TEST. It composes surfaces that already exist and
# are already tested per file. Including it does nothing; it must be run deliberately:
#
#     julia --project=spike -t auto spike/p14/run_p14_fdr_check.jl
#
# WHAT SC1-b MEASURES, AND WHAT IT DOES NOT.
# The rule emits a PREDICTED number: the running mean of the accepted items' null posteriors at
# k*, i.e. the model's own posterior expected false discovery proportion. The simulator knows each
# item's TRUE class, so the REALIZED false discovery proportion among the accepted set is a direct
# measurement rather than an estimate. SC1-b asks whether the two agree. Because both are computed
# on simulator draws from the same forward model the evidence net was trained on, agreement here is
# a WELL-SPECIFIED-REGIME result: it says the arithmetic and the calibration are consistent with
# each other, not that either survives real data. The `guarantee_basis` key records that.
#
# NO NEW DEPENDENCY, NO INSTALL, CPU-ONLY (CLAUDE.md, D-02). Bayesian FDR is a sort over posterior
# probabilities, so `MultipleTesting` would be the wrong tool as well as a new dependency in an
# environment a running test asserts byte-frozen.
#
# DECOUPLING: spike-local. `src/` is reached only READ-ONLY and transitively, and step 0 PROVES at
# run time that it is byte-unchanged AND that no untracked file has appeared under it.
#
# THE UNIT IS THE IMAGE PAIR (D-04). No per-region and no per-tile false discovery number is
# produced anywhere in this file, and that is not an omission: Phase 12 returned NO on calibrated
# per-region uncertainty and the per-tile map carries no uncertainty field at all, so there is
# nothing to control against. The persisted `unit_note` says so in the artifact.
#
# RUN THE UNIT FILES PER FILE, never through the aggregate suite -- it exits 1 early at the
# Phase-4 speedup gate and masks every later include block.

using JLD2
using Statistics
using Dates

# --- Guarded includes, in dependency order --------------------------------------------------------
# `decide.jl` transitively pulls the pre-registration, the provenance loader, the posterior, the
# FDR rule, the hedge, the fusion, the result type and the Phase-13 read surface. `pools.jl` owns
# the equal-thirds detector this runner re-applies to the persisted pool. Every include is guarded,
# so naming the modules this runner reads directly costs nothing and documents the dependency set.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :P14_PROVENANCE_LOADED)   || include(joinpath(@__DIR__, "provenance.jl"))
isdefined(@__MODULE__, :P14_POSTERIOR_LOADED)    || include(joinpath(@__DIR__, "posterior.jl"))
isdefined(@__MODULE__, :p14_bayes_fdr)           || include(joinpath(@__DIR__, "fdr.jl"))
isdefined(@__MODULE__, :P14_CONFORMAL_LOADED)    || include(joinpath(@__DIR__, "conformal.jl"))
isdefined(@__MODULE__, :P14_FUSE_LOADED)         || include(joinpath(@__DIR__, "fuse.jl"))
isdefined(@__MODULE__, :P14Result)               || include(joinpath(@__DIR__, "result.jl"))
isdefined(@__MODULE__, :P14_DECIDE_LOADED)       || include(joinpath(@__DIR__, "decide.jl"))
isdefined(@__MODULE__, :P14_POOLS_LOADED)        || include(joinpath(@__DIR__, "pools.jl"))

if !isdefined(@__MODULE__, :P14_FDR_RUNNER_LOADED)
    "The reported SC1-b / SC1-c artifact. Written BEFORE the SC1-b assertion can throw."
    const P14_FDR_REPORT_PATH = joinpath(@__DIR__, "p14_fdr_report.jld2")

    """
    The SHARED evaluation pool, produced by plan 14-09.

    The path is NAMED here rather than obtained by including `run_p14_conformal.jl`, because that
    file defines its own `main` and including one runner from another would put two zero-argument
    `main` methods in one namespace, where the last include silently wins. The two spellings of the
    path are pinned to each other by the pool's own recorded `P14_EVAL_COUNTER` and `n_eval`, which
    this runner asserts, so a path that pointed at some other file could not pass.
    """
    const P14_FDR_EVAL_POOL_PATH = joinpath(@__DIR__, "p14_eval_pool.jld2")

    """
    THE SC1-b BAND'S STANDARD-ERROR MULTIPLIER: the two-sided 95 % normal quantile.

    This is the one number in this file that is not read from `spike/p14/consts.jl`, and it is not
    a threshold: it is part of the BAND FORMULA itself, transcribed from `14-VALIDATION.md` row
    SC1-b and restated verbatim in `14-10-PLAN.md`'s objective, both frozen before this runner
    existed. Named here rather than buried inside the expression so it is greppable and so a reader
    can see there is exactly one such transcription.

    Nobody gets to pick the resulting band: fix the rule's predicted value and the realized number
    of discoveries and the interval follows.
    """
    const P14_SC1B_BAND_Z = 1.96

    """
    The tolerance for the section-E.1 sanity assertion: the three class posteriors must sum to one,
    and the composite null must equal the sum of its two parts, to this many absolute units on
    EVERY item.

    Cheap, and it catches every arithmetic slip in the composite pooling -- including the one that
    matters most, pooling a composite null by summing Bayes factors or logits instead of
    probabilities, which is not exact and would leave a plausible-looking number behind.
    """
    const P14_FDR_SUM_TOL = 1e-12

    "The columns this runner requires from the shared pool, so a schema drift fails by name."
    const P14_FDR_REQUIRED_POOL_KEYS = ("n", "n_eval", "P14_EVAL_COUNTER", "true_class",
                                        "decision", "abstain_reason", "null_posterior",
                                        "null_split_p_random", "null_split_p_exclusion",
                                        "class_posterior_coloc", "class_posterior_random",
                                        "class_posterior_exclusion", "logbf_coloc",
                                        "logbf_exclusion", "conformal_status", "ood_state",
                                        "qhat", "class_prior_used", "eval_masses", "idx")

    """
    The ONE row type of the SC1-b table, declared so all four rows are the SAME concrete
    NamedTuple.

    This is the structural half of T-14-06: `decided_fraction` is a FIELD OF THE TYPE, so a row
    without it cannot be constructed at all, rather than being caught by a reviewer noticing its
    absence. The band fields are `Union{Missing, Float64}` because a level at which the rule
    accepts nothing has no binomial interval -- `missing` says "there is no band here", where a
    `NaN` would read as "the band was computed and came out undefined".
    """
    const P14_FDR_ROW_KEYS = (:alpha, :k_star, :n_decided, :decided_fraction, :n_accepted,
                              :n_false_discoveries, :predicted_fdp, :realized_fdp,
                              :band_lo, :band_hi, :within_band, :exceeds_upper_band,
                              :t_star, :cost_ratio, :mean_p_random, :mean_p_exclusion,
                              :share_random, :share_exclusion, :note)
    const P14FDRRow = NamedTuple{P14_FDR_ROW_KEYS,
                                 Tuple{Float64, Int, Int, Float64, Int, Int,
                                       Float64, Float64,
                                       Union{Missing, Float64}, Union{Missing, Float64},
                                       Bool, Bool,
                                       Float64, Float64, Float64, Float64,
                                       Float64, Float64, String}}

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else can
    # make this body skip.
    const P14_FDR_RUNNER_LOADED = true
end

# =============================================================================================
# Small local helpers
# =============================================================================================

"""
    _p14_fdr_save_report(path, required; kwargs...) -> String

Atomically persist an artifact: write `path * ".tmp"`, REOPEN it read-only and integrity-check
every key in `required`, then `mv(...; force = true)`.

The `_p13_gate_save_report` idiom (`spike/p13/run_three_way_gate.jl:171-190`), carried per runner
exactly as Phase 13 carried it per runner. It is deliberately NOT imported from
`run_p14_conformal.jl`: importing it would mean including a second runner and inheriting its
`main`, and a twelve-line save helper is a smaller cost than two `main` methods in one namespace.

A crash mid-write leaves a discardable `.tmp` rather than a torn artifact a later reader would
happily believe.
"""
function _p14_fdr_save_report(path, required; kwargs...)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        for k in required
            @assert haskey(f, k) "_p14_fdr_save_report: integrity check failed -- $tmp is missing the required key `$k`"
        end
    end
    mv(tmp, path; force = true)
    return path
end

"""
    p14_sc1b_band(predicted, k_star) -> (band_lo, band_hi)

The binomial 95 % interval around the rule's OWN predicted false discovery proportion, at a sample
size equal to the realized number of discoveries: `predicted +- z * sqrt(p(1-p)/k*)`.

Returns `(missing, missing)` when `k_star == 0`. THE DEGENERATE BRANCH COMES FIRST, BEFORE THE
ARITHMETIC: a level at which the rule accepts nothing has no realized proportion to compare and no
sample size to compute a standard error at, and reaching that case through a division by zero would
put a `NaN` band into the table that a later reader would have to reverse-engineer.

# Why the band is around the PREDICTED value rather than around the realized one

The question SC1-b asks is whether the realized proportion is consistent with what the rule
CLAIMED. That makes the rule's claim the null hypothesis and the realized count the observation, so
the sampling interval belongs around the claim. An interval around the observation would be
answering the reverse question, and at small `k*` the two are not interchangeable.
"""
function p14_sc1b_band(predicted::Real, k_star::Integer)
    k_star > 0 || return (missing, missing)
    p  = Float64(predicted)
    se = sqrt(max(p * (1.0 - p), 0.0) / k_star)
    return (p - P14_SC1B_BAND_Z * se, p + P14_SC1B_BAND_Z * se)
end

"""
    _p14_fdr_class_from_key(s::Symbol) -> ThreeWayClass

The inverse of the class-symbol map the shared pool was written through, by an explicit
three-branch map.

Written out rather than derived from the enum's spelling for the same reason `_p14_class_key` is:
five class orderings coexist in this codebase and every statistic here is invariant to a CONSISTENT
relabelling, which is exactly what makes an inconsistent one produce numbers that look right. A
fall-through would silently drop a class out of the equal-thirds guard below, so the map is total
and errors otherwise.
"""
function _p14_fdr_class_from_key(s::Symbol)
    s === :coloc     && return COLOC
    s === :random    && return RANDOM
    s === :exclusion && return EXCLUSION
    error("_p14_fdr_class_from_key: unhandled class symbol `$s`. The inverse class map must be " *
          "total; a fall-through would drop a class out of every mass this runner re-checks.")
end

"""
    _p14_fdr_load_pool(path; reported = true) -> NamedTuple

Load the SHARED evaluation pool written by plan 14-09 and prove it is the pool this runner is
entitled to score, rather than merely a file with the right name.

Four checks, each catching something the others cannot:

 1. the file EXISTS -- and if it does not, this fails loudly naming 14-09 and the command that
    produces it. It never redraws. A silent redraw at the same counter would still land on the same
    items, but a redraw is a second experiment as soon as anything about the draw changes, and the
    whole point of the persisted pool is that SC1-b, SC1-d and SC3-a/b/c share one;
 2. every column this runner reads is PRESENT, by name, so a schema drift fails with the missing
    column rather than with a `KeyError` deep inside a comprehension;
 3. the recorded `P14_EVAL_COUNTER` and `n_eval` are the FROZEN ones -- a pool drawn at another
    counter is a different set of items and may not be scored here;
 4. the recorded class masses are re-derived from the stored labels through the frozen
    `p14_assert_unstratified` guard, which rejects an equal-thirds draw and then checks every mass
    against the measured prior. Re-running the guard rather than trusting the stored
    `eval_masses` is what makes "the pool is unstratified" a checked property of THIS run.

The guard in (4) is applied only on a reported pool: 14-09 recorded that below roughly n = 400 the
equal-thirds detector has no power (three binomial standard errors around one third is 0.129 at
n = 120), so a smoke pool would trip it as a small-sample artifact rather than as a finding.
"""
function _p14_fdr_load_pool(path::AbstractString; reported::Bool = true)
    isfile(path) || error("""
        run_p14_fdr_check: the SHARED evaluation pool is ABSENT at
            $path

        That pool is produced by PLAN 14-09 (`spike/p14/run_p14_conformal.jl`), which draws it ONCE
        at the reserved P14_EVAL_COUNTER and persists it so that SC1-b, SC1-d and SC3-a/b/c are
        scored on the SAME items. Run it first:

            julia --project=spike -t auto spike/p14/run_p14_conformal.jl

        THIS RUNNER WILL NOT REDRAW THE EVALUATION SET, silently or otherwise. A redraw would put
        the realized-FDP number on a different sample from the coverage number it is quoted beside,
        and nothing in either artifact would say so.""")

    d = JLD2.load(path)
    for k in P14_FDR_REQUIRED_POOL_KEYS
        haskey(d, k) || error(
            "run_p14_fdr_check: the shared evaluation pool at $path carries no `$k` column. The " *
            "14-09 pool schema moved under this runner; scoring SC1-b against a partially " *
            "understood artifact is how a number ends up meaning something other than its name.")
    end

    n = Int(d["n"])
    @assert n == Int(d["n_eval"]) "run_p14_fdr_check: the pool stores $(n) items but records n_eval = $(d["n_eval"]); the artifact disagrees with itself and neither number can be quoted"
    @assert Int(d["P14_EVAL_COUNTER"]) == Int(P14_EVAL_COUNTER) "run_p14_fdr_check: the pool was drawn at counter $(d["P14_EVAL_COUNTER"]) but the frozen reserved counter is $(P14_EVAL_COUNTER); a pool from another stream is a different set of items and may not be scored as the reported evaluation set"
    if reported
        @assert n == Int(P14_N_EVAL) "run_p14_fdr_check: the pool holds $n items but the frozen reported evaluation size is $(P14_N_EVAL)"
    end

    true_class = Vector{Symbol}(d["true_class"])
    @assert length(true_class) == n "run_p14_fdr_check: the pool's true_class column holds $(length(true_class)) entries for $n items"

    # (4) THE EQUAL-THIRDS GUARD, RE-RUN rather than read off the artifact. `p14_assert_unstratified`
    # is the frozen surface; a local restatement of "is this balanced" would be a second spelling of
    # a check that already exists, which is how two versions of one guard silently diverge.
    realized = if reported
        p14_assert_unstratified([_p14_fdr_class_from_key(s) for s in true_class])
    else
        class_masses([_p14_fdr_class_from_key(s) for s in true_class])
    end
    stored = d["eval_masses"]
    for k in P14_MASS_KEYS
        @assert getproperty(realized, k) == getproperty(stored, k) "run_p14_fdr_check: the `$k` mass re-derived from the pool's own labels is $(getproperty(realized, k)) while the pool records $(getproperty(stored, k)); the stored summary and the stored labels disagree"
    end

    return (n = n,
            path = path,
            true_class = true_class,
            decision = Vector{Symbol}(d["decision"]),
            abstain_reason = Vector{Symbol}(d["abstain_reason"]),
            v = Vector{Float64}(d["null_posterior"]),
            p_random = Vector{Float64}(d["null_split_p_random"]),
            p_exclusion = Vector{Float64}(d["null_split_p_exclusion"]),
            p_coloc = Vector{Float64}(d["class_posterior_coloc"]),
            post_random = Vector{Float64}(d["class_posterior_random"]),
            post_exclusion = Vector{Float64}(d["class_posterior_exclusion"]),
            logbf_coloc = Vector{Float64}(d["logbf_coloc"]),
            logbf_exclusion = Vector{Float64}(d["logbf_exclusion"]),
            conformal_status = Vector{Symbol}(d["conformal_status"]),
            ood_state = Vector{Symbol}(d["ood_state"]),
            idx = Vector{Int}(d["idx"]),
            qhat = Float64(d["qhat"]),
            prior = d["class_prior_used"],
            eval_masses = realized,
            pool_provenance = get(d, "pool_provenance", nothing),
            pool_generated = get(d, "generated", ""),
            pool_reported = get(d, "reported", false))
end

"""
    _p14_fdr_partition_sanity(pool) -> NamedTuple

The section-E.1 sanity assertion, on EVERY item: the three class posteriors sum to one, and the
composite null equals the sum of its two parts.

Both are cheap and neither is decorative. The composite null `v = P(random | Z) + P(exclusion | Z)`
is exact only because posterior probability is additive over a genuine PARTITION; if the three
masses did not sum to one, the object being pooled would not be a partition and the FDR sort would
be running over a quantity with no interpretation. Returns the realized maximum deviations so the
artifact records HOW exact rather than merely that an assertion passed.
"""
function _p14_fdr_partition_sanity(pool::NamedTuple)
    dev_sum  = 0.0
    dev_pool = 0.0
    for i in 1:pool.n
        s = pool.p_coloc[i] + pool.post_random[i] + pool.post_exclusion[i]
        dev_sum = max(dev_sum, abs(s - 1.0))
        # The two split columns must reconstruct the scalar null the rule actually sorted on.
        dev_pool = max(dev_pool, abs((pool.p_random[i] + pool.p_exclusion[i]) - pool.v[i]))
    end
    @assert dev_sum < P14_FDR_SUM_TOL "run_p14_fdr_check (14-RESEARCH section E.1): the three class posteriors do not sum to 1 on every item (max absolute deviation $dev_sum, tolerance $(P14_FDR_SUM_TOL)). The composite null is additive ONLY over a genuine partition, so a sort over `v` would be a sort over a quantity with no interpretation."
    @assert dev_pool < P14_FDR_SUM_TOL "run_p14_fdr_check (14-RESEARCH section E.1): the stored composite null does not equal p_random + p_exclusion on every item (max absolute deviation $dev_pool, tolerance $(P14_FDR_SUM_TOL)). Either the pooling was done on the wrong object or the two columns came from different computations."
    return (max_abs_dev_sum_to_one = dev_sum, max_abs_dev_null_pooling = dev_pool)
end

"""
    _p14_fdr_row(; kwargs...) -> P14FDRRow

Build one SC1-b table row BY NAME and convert it to the single declared row type.

The conversion is what makes "every row carries its decided fraction" structural rather than
conventional: a row missing a key, or carrying an extra one, cannot be converted at all. The keys
are asserted equal to `P14_FDR_ROW_KEYS` first so the failure names the offending key set rather
than surfacing as a conversion `MethodError`.
"""
function _p14_fdr_row(nt::NamedTuple)
    keys(nt) === P14_FDR_ROW_KEYS || throw(ArgumentError(
        "_p14_fdr_row: a table row must carry exactly $(P14_FDR_ROW_KEYS) in that order; got " *
        "$(keys(nt)). The row type is the structural guarantee that no FDR number is recorded " *
        "without its decided fraction."))
    return convert(P14FDRRow, nt)
end

# =============================================================================================
# The prior-sensitivity sweep (SC1-c) -- REPORTED, NOT GATED
# =============================================================================================

"""
    _p14_fdr_recs(pool) -> Vector{NamedTuple}

Rebuild the rule inputs that `_p14_rule_at_prior` consumes, from the persisted pool columns.

Three fields, and each is exactly what the sweep needs:

- `lbf`: the evidence triple, rebuilt BY NAME from the two stored log Bayes factors with the D-08
  structural zero on the reference class. The evidence does not depend on the prior, which is why
  14-09 stored it -- so the sweep re-derives the posterior at another prior without re-running the
  net at all.
- `ood`: the three-valued OOD state. It does not depend on the prior either, so it is carried.
- `cm`: the cross-method record. The shared pool was drawn WITHOUT images -- `p14_draw_pool` yields
  summary vectors, not `MultiChannelImage`s -- so the classical channels were never consulted and
  the record is the explicit `:not_computed` one, carrying its note. **Never a bare `false`:** a
  `disagree_any = false` with no note reads as "the classics were consulted and concurred", which
  is not what happened.
"""
function _p14_fdr_recs(pool::NamedTuple)
    cm = _p14_no_cross_method(
        "the shared evaluation pool carries summary vectors, not image pairs, so no classical " *
        "comparison was made for any item in it; this is NOT agreement with the classics")
    return [(lbf = (coloc = pool.logbf_coloc[i], random = 0.0,
                    exclusion = pool.logbf_exclusion[i]),
             ood = pool.ood_state[i],
             cm  = cm)
            for i in 1:pool.n]
end

"""
    _p14_fdr_assert_reconstruction(recs, pool, prior, qhat) -> NamedTuple

Prove that re-running the WHOLE rule from the rebuilt inputs at the pool's OWN prior reproduces the
pool, before any swept rung is believed.

**This is what makes the sweep a re-run rather than a second experiment.** The sweep goes through
`_p14_rule_at_prior`, the same code path `decide_coloc` uses; the SC1-b table above goes through the
recorded per-item fields. Those are two routes to one set of numbers, and a sweep that quietly
disagreed with the table at the very rung where they should coincide would keep producing a
perfectly well-formed curve with nothing saying which route was wrong. Four things are pinned: the
three class posteriors, the composite null, the conformal status, and the abstain/decide partition.

Returns the realized maximum deviations, so the artifact records HOW exactly the two agree rather
than merely that an assertion passed.
"""
function _p14_fdr_assert_reconstruction(recs, pool::NamedTuple, prior::NamedTuple, qhat::Real)
    run = _p14_rule_at_prior(recs, qhat, prior, first(P14_ALPHA_FDR_GRID), false, pool.n)
    dev_post = 0.0
    dev_null = 0.0
    for i in 1:pool.n
        p = run.ps[i]
        dev_post = max(dev_post, abs(p.coloc - pool.p_coloc[i]),
                                 abs(p.random - pool.post_random[i]),
                                 abs(p.exclusion - pool.post_exclusion[i]))
        dev_null = max(dev_null, abs(p14_null_posterior(p) - pool.v[i]))
        @assert run.csets[i].status === pool.conformal_status[i] "run_p14_fdr_check: item $i re-runs to conformal status $(run.csets[i].status) but the pool recorded $(pool.conformal_status[i]); the sweep is not reproducing the run it sweeps around"
        @assert (run.acts[i].action === :abstain) == (pool.decision[i] === :abstain) "run_p14_fdr_check: item $i re-runs to action $(run.acts[i].action) but the pool recorded decision $(pool.decision[i]); the abstain/decide partition has drifted between the two routes"
    end
    @assert dev_post < P14_FDR_SUM_TOL "run_p14_fdr_check: the re-run class posterior deviates from the pool by up to $dev_post (tolerance $(P14_FDR_SUM_TOL)); the prior sweep would be sweeping around a different point from the one SC1-b was measured at"
    @assert dev_null < P14_FDR_SUM_TOL "run_p14_fdr_check: the re-run composite null deviates from the pool by up to $dev_null (tolerance $(P14_FDR_SUM_TOL))"
    return (max_abs_dev_posterior = dev_post, max_abs_dev_null = dev_null,
            n_decided = run.fdr.n_decided)
end

"""
    _p14_fdr_claims() -> NamedTuple

The honesty block, written INTO the artifact as string keys so it travels WITH the numbers and
cannot be dropped by a report writer.

Every one of these is a sentence a reader of a bare FDR number would otherwise have to be told by
someone who remembered to tell them. A `.jld2` that carries the number and not the scope is a
repudiation hole (T-14-33), which is why these are REQUIRED keys of the save rather than optional
extras.
"""
function _p14_fdr_claims()
    return (
        fdr_claim =
            "The posterior expected false discovery proportion is controlled at alpha_FDR over " *
            "the DECIDED SUBSET ONLY -- the batch minus abstentions. It is NOT controlled over " *
            "the whole batch: an abstained item is neither a discovery nor a non-discovery and " *
            "appears in neither the numerator nor the denominator. The guarantee is conditional " *
            "on the model -- calibrated posteriors and a correct class prior -- and is a " *
            "POSTERIOR EXPECTATION, not a frequentist long-run rate. A rule that abstains on " *
            "most of a batch meets any level trivially, so the level and the decided fraction " *
            "are ONE quantity and are never quoted apart.",

        assumption_a2_statement =
            "ASSUMPTION A2, as a derivation a referee can check in two lines rather than an " *
            "assertion to be believed. (1) The controlled quantity is E[FDP | data] = (1/|R|) * " *
            "sum over i in R of P(H0_i | data), where R is the rejection set. (2) R is a " *
            "DETERMINISTIC FUNCTION OF THE OBSERVED DATA -- both the abstention filter and the " *
            "prefix sort are computed from the data alone -- hence R is sigma(data)-measurable " *
            "and pulls straight out of the conditional expectation, so the identity holds for " *
            "ANY data-dependent selection of R, including one that filtered on OOD status or " *
            "conformal ambiguity. Unlike the frequentist case, where a data-dependent " *
            "pre-selection is a genuine selective-inference problem requiring correction, the " *
            "Bayesian posterior quantity is immune because the conditioning has already happened. " *
            "This is LOAD-BEARING: it is the entire justification for abstain-then-sort and for " *
            "the decided-subset scope.",

        ordering_note =
            "ABSTAIN FIRST, THEN SORT. The reverse order -- sort the whole batch, then abstain " *
            "from some of the accepted items -- is WRONG and is recorded here so it is not " *
            "re-derived: removing items from a set whose running mean was computed over all of " *
            "them changes both the numerator and the denominator in an uncontrolled way, so the " *
            "residual accepted set's running mean can EXCEED the level. Abstaining first also " *
            "improves the premise, because it routes away precisely the items whose posteriors " *
            "are least likely to be right.",

        independence_note =
            "NO INDEPENDENCE ACROSS THE BATCH IS ASSUMED. E[sum of 1{H0_i}] = sum of v_i follows " *
            "from LINEARITY OF EXPECTATION alone, which holds whether or not the items are " *
            "independent. This is the standard reviewer question about a batch-level error rate " *
            "and the answer is clean; it is recorded here so the answer travels with the number " *
            "instead of having to be reconstructed.",

        prior_note =
            "THE PRIOR IS THE SIMULATOR'S CLASS MIX, NOT ANY REAL BATCH'S PREVALENCE. " *
            "pi_class_masses is the realized class mix of the prior draws the evidence net was " *
            "trained under -- measured, never chosen -- and every null posterior v_i depends on " *
            "it. A batch that is 90 % coloc makes every v_i too large, so the rule is merely too " *
            "conservative; a batch that is 1 % coloc makes them too small and THE FDR CLAIM IS " *
            "THEN FALSE. The pi-sensitivity table is the honest answer to that and is the reason " *
            "it is required output rather than an optional extra. Empirical-Bayes estimation of " *
            "the prior from the batch is DELIBERATELY NOT the default: estimating the prior from " *
            "the same batch the guarantee is quoted over would make that guarantee CIRCULAR, and " *
            "it would add a second estimated quantity with its own failure modes.",

        prior_atom_note =
            "THE PRIOR CARRIES CLAMP ATOMS, AND THEY ARE ASYMMETRIC. The theta prior places " *
            "roughly 4.85 % of its mass at the clamped endpoint rho = -0.99 against roughly " *
            "1.91 % at rho = +0.99, a ratio of about 2.5 to 1 AGAINST the exclusion end. The " *
            "measured pi_class_masses inherit those atoms, and in this phase the prior enters " *
            "the PRIOR TERM of the class posterior directly -- a term Phase 13 never used, " *
            "because Phase 13 published log Bayes factors and deliberately declined to introduce " *
            "a prior at all. So an asymmetry that was inert upstream is load-bearing here. " *
            "Recorded as a stated property of the prior, and as one more reason the " *
            "pi-sensitivity table is required rather than optional.",

        unit_note =
            "FDR IS CONTROLLED PER IMAGE PAIR (D-04). One call, one pair, one decision, one " *
            "entry in the batch the prefix rule sorts. A per-region or per-tile false discovery " *
            "rate is OUT OF SCOPE and is not merely unimplemented: the per-tile map type carries " *
            "no uncertainty field at all, and Phase 12 returned NO on calibrated per-tile " *
            "uncertainty, so there is nothing to control against. A tile map may be DISPLAYED " *
            "beside a decision; it is never FDR-controlled, and claiming otherwise would repeat " *
            "the counting-a-closed-negative-as-a-success failure this project has already " *
            "corrected once.",

        amendment = P14_AMENDMENT_NOTICE,
    )
end

# =============================================================================================
# The reported run
# =============================================================================================

"""
    main(; pool_path = P14_FDR_EVAL_POOL_PATH, report_path = P14_FDR_REPORT_PATH,
           verbose = true) -> NamedTuple

The reported SC1-b measurement: realized versus predicted false discovery proportion at every
level of the frozen `P14_ALPHA_FDR_GRID`, each row carrying its decided fraction.

`pool_path` defaults to the SHARED evaluation pool 14-09 persisted. Passing anything else puts the
run in SMOKE mode: it writes to a `_smoke` artifact path and prints a banner saying so, so a load
check can never be mistaken for, or quoted as, a verdict, and can never overwrite the reported
artifact.
"""
function main(; pool_path = P14_FDR_EVAL_POOL_PATH,
                report_path = P14_FDR_REPORT_PATH,
                verbose::Bool = true)

    t_start  = time()
    reported = (pool_path == P14_FDR_EVAL_POOL_PATH)
    reported || (report_path = replace(report_path, ".jld2" => "_smoke.jld2"))

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
        @assert env_clean "D-02 decoupling breach: spike/Project.toml or spike/Manifest.toml is modified -- the frozen spike environment moved under the run, and the FDR rule needs no dependency at all"
        @assert isempty(src_untracked) "D-01 decoupling breach: an UNTRACKED file appeared under src/ ($(src_untracked)) -- a diff against HEAD cannot see it, and shipping the decision layer by stealth is precisely what D-01 forbids"
    end

    # --- 1. THE INHERITANCE (D-07) AND THE SHARED POOL --------------------------------------
    bundle = p14_load_bundle()
    pool   = _p14_fdr_load_pool(pool_path; reported = reported)
    n_eval = pool.n

    # The pool was scored at the MEASURED prior; the bundle is where that prior comes from. Pinned
    # rather than assumed: a pool scored under some other prior would make every `v` in it a
    # different quantity from the one this runner's sweep recomputes.
    for k in P14_CLASS_KEYS
        @assert getproperty(pool.prior, k) == getproperty(bundle.prior, k) "run_p14_fdr_check: the pool records a `$k` prior mass of $(getproperty(pool.prior, k)) while the inherited bundle measures $(getproperty(bundle.prior, k)); the pool was scored under a different prior from the one this runner would sweep around"
    end

    # --- 2. THE LOCKED-THRESHOLD BANNER, BEFORE ANYTHING IS COMPUTED -------------------------
    if verbose
        println("="^78)
        println("PHASE-14 REPORTED SC1-b: REALIZED vs PREDICTED FALSE DISCOVERY PROPORTION")
        println("="^78)
        reported || println("*** SMOKE MODE -- NOT A REPORTED RUN, THE BAND IS NOT ASSERTED ***")
        println("AMENDED CRITERIA: ", P14_AMENDMENT_NOTICE)
        println("GATED (SC1-b):")
        println("  P14_ALPHA_FDR_GRID      = $P14_ALPHA_FDR_GRID")
        println("                            (consts.jl section C. alpha_FDR is a per-call USER")
        println("                            parameter -- SC1 -- and is NEVER the conformal")
        println("                            miscoverage level, which is a different quantity that")
        println("                            happens to share one grid point (D-07). What is frozen")
        println("                            is the SWEEP, so `we controlled FDR' is not a claim")
        println("                            about a grid chosen after the fact.)")
        println("  P14_N_EVAL              = $P14_N_EVAL           (consts.jl section B; pool holds $n_eval)")
        println("  P14_EVAL_COUNTER        = $P14_EVAL_COUNTER               (the SHARED reported evaluation set)")
        println("  band                    = predicted +- $P14_SC1B_BAND_Z * sqrt(predicted*(1-predicted)/k*)")
        println("                            DERIVED from the rule's OWN predicted value and the")
        println("                            realized k* (14-VALIDATION.md row SC1-b). Not a chosen")
        println("                            number. FALSIFIED if the realized FDP EXCEEDS the")
        println("                            upper band at any grid level.")
        println("REPORTED, NOT GATED:")
        println("  the exclusion/random split of the accepted sets' null mass, the four cost")
        println("  ratios and the four acceptance thresholds are reported beside the verdict and")
        println("  none of them is gated.")
        println("SCOPE (D-04, and it is the whole claim):")
        println("  FDR is controlled over the DECIDED SUBSET ONLY, per IMAGE PAIR. No per-region")
        println("  and no per-tile number is produced anywhere in this runner.")
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

    # --- 3. THE SECTION-E.1 SANITY ASSERTION, ON EVERY ITEM ----------------------------------
    verbose && println("[1/5] checking the composite-null partition on every item …")
    sanity = _p14_fdr_partition_sanity(pool)
    verbose && println("      max |sum(p) - 1| = $(sanity.max_abs_dev_sum_to_one)   " *
                       "max |p_random + p_exclusion - v| = $(sanity.max_abs_dev_null_pooling)")

    # --- 4. ABSTAIN FIRST. The partition happens before anything is sorted. -------------------
    # The recorded per-item fused action IS the partition: `p14_decide_one` wrote `:abstain` when
    # the D-05 fusion abstained and an evidence call otherwise. Reading it back rather than
    # recomputing it is what makes SC1-b a measurement ON the shared pool rather than beside it.
    verbose && println("[2/5] partitioning the pool -- ABSTAIN FIRST, then sort …")
    decided_idx = [i for i in 1:n_eval if pool.decision[i] !== :abstain]
    n_decided   = length(decided_idx)
    abstained   = n_eval - n_decided
    reasons     = Dict{Symbol, Int}()
    for i in 1:n_eval
        pool.decision[i] === :abstain || continue
        reasons[pool.abstain_reason[i]] = get(reasons, pool.abstain_reason[i], 0) + 1
    end
    @assert n_decided + abstained == n_eval "run_p14_fdr_check: the decided/abstained partition does not cover the batch"
    verbose && println("      decided $n_decided / $n_eval   abstained $abstained   reasons: $reasons")

    v_decided = pool.v[decided_idx]

    # --- 5. THE SWEEP OVER THE FROZEN LEVELS -------------------------------------------------
    verbose && println("[3/5] sweeping the frozen alpha grid …")
    rows = P14FDRRow[]
    for alpha in P14_ALPHA_FDR_GRID
        # `p14_fdr_over_decided` REQUIRES `n_total`, so a row without a decided fraction is
        # structurally impossible: the fraction comes out of the rule, not out of a later merge.
        res = p14_fdr_over_decided(v_decided, alpha; n_total = n_eval)
        @assert res.fdr_scope === :decided_subset_only "run_p14_fdr_check: the rule returned scope $(res.fdr_scope)"
        @assert res.n_decided == n_decided "run_p14_fdr_check: the rule saw $(res.n_decided) decided items but the partition holds $n_decided"

        # `res.accepted` holds positions into the DECIDED SUBSET, never into the batch. The mapping
        # back is explicit and then asserted: an off-by-one here is silent and would report the
        # wrong items as discoveries.
        accepted_original = [decided_idx[j] for j in res.accepted]
        @assert length(accepted_original) == res.k_star "run_p14_fdr_check: the mapped accepted set holds $(length(accepted_original)) items but k_star is $(res.k_star)"
        @assert all(i -> pool.decision[i] !== :abstain, accepted_original) "run_p14_fdr_check: an ABSTAINED item reached the accepted set at level $alpha; the decided-subset index mapping is off"

        # THE REALIZED FDP IS A DIRECT MEASUREMENT. The simulator's true label exists for every
        # item, so this is a count of accepted items whose true class is not the coloc class --
        # not an estimate of one.
        n_false   = count(i -> pool.true_class[i] !== :coloc, accepted_original)
        realized  = res.k_star == 0 ? NaN : n_false / res.k_star
        predicted = res.fdp
        band_lo, band_hi = p14_sc1b_band(predicted, res.k_star)

        # The composite-null decomposition of the ACCEPTED set: a batch whose false discoveries sit
        # mostly on exclusion mass is confusing segregation with colocalization, one whose sit
        # mostly on random mass is confusing noise with signal, and the total hides which.
        mean_pr = res.k_star == 0 ? NaN : mean(pool.p_random[accepted_original])
        mean_pe = res.k_star == 0 ? NaN : mean(pool.p_exclusion[accepted_original])
        null_mass = res.k_star == 0 ? NaN : sum(pool.v[accepted_original])
        share_r = (res.k_star == 0 || null_mass == 0) ? NaN :
                  sum(pool.p_random[accepted_original]) / null_mass
        share_e = (res.k_star == 0 || null_mass == 0) ? NaN :
                  sum(pool.p_exclusion[accepted_original]) / null_mass

        exceeds = res.k_star > 0 && realized > band_hi
        within  = res.k_star > 0 && realized >= band_lo && realized <= band_hi

        push!(rows, _p14_fdr_row((
            alpha               = Float64(alpha),
            k_star              = Int(res.k_star),
            n_decided           = Int(res.n_decided),
            decided_fraction    = Float64(res.decided_fraction),
            n_accepted          = length(accepted_original),
            n_false_discoveries = Int(n_false),
            predicted_fdp       = Float64(predicted),
            realized_fdp        = Float64(realized),
            band_lo             = band_lo,
            band_hi             = band_hi,
            within_band         = within,
            exceeds_upper_band  = exceeds,
            t_star              = Float64(res.t_star),
            cost_ratio          = Float64(res.cost_ratio),
            mean_p_random       = Float64(mean_pr),
            mean_p_exclusion    = Float64(mean_pe),
            share_random        = Float64(share_r),
            share_exclusion     = Float64(share_e),
            note                = res.k_star == 0 ? "no discoveries at this alpha" : "")))
    end
    @assert length(rows) == length(P14_ALPHA_FDR_GRID) "run_p14_fdr_check: the table has $(length(rows)) rows for $(length(P14_ALPHA_FDR_GRID)) frozen levels"

    offenders = [r for r in rows if r.exceeds_upper_band]
    sc1b_met  = isempty(offenders)

    # --- 6. THE PRIOR-SENSITIVITY SWEEP (SC1-c) -- REPORTED, NOT GATED -----------------------
    # There is NO BAR anywhere below and there is not meant to be one. SC1-c's failure condition
    # is FAILING TO REPORT IT; a threshold attached here would quietly convert a diagnostic into
    # a gate, which is the thing the pre-registration forbids.
    verbose && println("[4/5] re-running the whole rule at every rung of the frozen prior grid …")
    recs  = _p14_fdr_recs(pool)
    recon = _p14_fdr_assert_reconstruction(recs, pool, bundle.prior, pool.qhat)
    @assert recon.n_decided == n_decided "run_p14_fdr_check: the re-run decides $(recon.n_decided) items while the recorded partition decides $n_decided"
    verbose && println("      reconstruction verified: max |dp| = $(recon.max_abs_dev_posterior), " *
                       "max |dv| = $(recon.max_abs_dev_null), decided $(recon.n_decided)")

    pi_rows = NamedTuple[]
    for pi_c in P14_PI_COLOC_GRID
        # The coloc mass is replaced and the other two are renormalised IN THEIR MEASURED RATIO,
        # so the sweep asks "what if this batch's coloc prevalence were different" rather than
        # "what if the whole class structure were different" -- two questions, and only one of
        # them was asked.
        prior_c = _p14_prior_at(bundle.prior, pi_c)
        for alpha in P14_ALPHA_FDR_GRID
            # THROUGH THE SAME CODE PATH the headline run uses. The prior moves the class
            # posterior, hence the conformal set, hence the abstention, hence the amortized call
            # and the fusion -- so every rung is a full re-run and `n_decided` is free to move.
            sr = _p14_rule_at_prior(recs, pool.qhat, prior_c, alpha, false, n_eval)
            acc_orig = [sr.decided_idx[j] for j in sr.fdr.accepted]
            @assert length(acc_orig) == sr.fdr.k_star "run_p14_fdr_check: the swept accepted set at pi_c = $pi_c, alpha = $alpha holds $(length(acc_orig)) items but k_star is $(sr.fdr.k_star)"
            n_false_c = count(i -> pool.true_class[i] !== :coloc, acc_orig)
            push!(pi_rows, (pi_coloc         = Float64(pi_c),
                            alpha            = Float64(alpha),
                            k_star           = Int(sr.fdr.k_star),
                            predicted_fdp    = Float64(sr.fdr.fdp),
                            realized_fdp     = sr.fdr.k_star == 0 ? NaN :
                                               n_false_c / sr.fdr.k_star,
                            n_decided        = Int(sr.fdr.n_decided),
                            decided_fraction = Float64(sr.fdr.decided_fraction)))
        end
    end
    @assert length(pi_rows) == length(P14_PI_COLOC_GRID) * length(P14_ALPHA_FDR_GRID) "run_p14_fdr_check: the sensitivity table has $(length(pi_rows)) rows for a $(length(P14_PI_COLOC_GRID)) x $(length(P14_ALPHA_FDR_GRID)) grid"
    pi_sensitivity = [NamedTuple{(:pi_coloc, :alpha, :k_star, :predicted_fdp, :realized_fdp,
                                  :n_decided, :decided_fraction),
                                 Tuple{Float64, Float64, Int, Float64, Float64, Int, Float64}}(r)
                      for r in pi_rows]

    # --- 7. PERSIST THE REPORT, BEFORE THE HEADLINE AND BEFORE THE ASSERTION -----------------
    verbose && println("[5/5] persisting the SC1-b / SC1-c report (BEFORE any verdict) …")
    consts_path = joinpath(@__DIR__, "consts.jl")
    _p14_fdr_save_report(report_path,
        ("alpha_table", "fdr_scope", "sc1b_met", "n_eval", "n_decided", "decided_fraction",
         "amendment", "guarantee_basis", "pi_sensitivity", "pi_sensitivity_gated",
         "fdr_claim", "assumption_a2_statement", "ordering_note", "independence_note",
         "prior_note", "prior_atom_note", "unit_note");
        # --- the honesty block, splatted in so it CANNOT be dropped: every one of its keys is
        #     also a REQUIRED key of the integrity check above, so an artifact missing one is
        #     never written at all (T-14-33).
        _p14_fdr_claims()...,
        # --- the measurement ---
        alpha_table = rows,
        alpha_table_keys = collect(P14_FDR_ROW_KEYS),
        sc1b_met = sc1b_met,
        n_eval = n_eval, n_decided = n_decided, n_abstained = abstained,
        decided_fraction = n_decided / n_eval,
        abstain_reason_counts = NamedTuple{Tuple(keys(reasons))}(Tuple(values(reasons))),
        # --- the scope label, machine-readable, so no reader can take this for batch-wide control
        fdr_scope = :decided_subset_only,
        # --- SC1-c: REPORTED, NOT GATED. The flag is persisted rather than left to prose so a
        #     later reader cannot mistake the table for a criterion that was met.
        pi_sensitivity = pi_sensitivity,
        pi_sensitivity_gated = false,
        pi_sensitivity_reconstruction = recon,
        # --- the section-E.1 sanity numbers, recorded rather than merely asserted ---
        max_abs_dev_sum_to_one = sanity.max_abs_dev_sum_to_one,
        max_abs_dev_null_pooling = sanity.max_abs_dev_null_pooling,
        partition_tolerance = P14_FDR_SUM_TOL,
        # --- EVERY THRESHOLD THIS RUN WAS SCORED AGAINST, written INTO the artifact ---
        P14_ALPHA_FDR_GRID = P14_ALPHA_FDR_GRID,
        P14_PI_COLOC_GRID = P14_PI_COLOC_GRID,
        P14_N_EVAL = P14_N_EVAL,
        P14_EVAL_COUNTER = Int(P14_EVAL_COUNTER),
        P14_SC1B_BAND_Z = P14_SC1B_BAND_Z,
        P14_ITERATION_ALLOWANCE = P14_ITERATION_ALLOWANCE,
        P14_ITERATION_TRIGGER = P14_ITERATION_TRIGGER,
        iteration_trigger_fired = false,
        iteration_allowance_applies_here = false,
        # --- what the numbers actually rest on ---
        guarantee_basis = :simulator_derived_held_out_draws,
        named_limits = p14_named_limits(),
        # `amendment` is NOT repeated here: it is one of the honesty keys splatted in above, so
        # there is exactly one spelling of it and it cannot drift from the others.
        # --- provenance ---
        named_limits_count = length(p14_named_limits()),
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
        consts_sha = p14_consts_sha(),
        consts_git_blob_sha = p14_blob_sha(consts_path),
        reported = reported,
        julia_version = string(VERSION),
        nthreads = Threads.nthreads(),
        elapsed_s = time() - t_start,
        generated = string(Dates.now(Dates.UTC)) * "Z")
    verbose && println("      persisted (before any verdict) -> $report_path")

    # --- 7. THE HEADLINE, PRINTED BEFORE THE ASSERTION CAN THROW ------------------------------
    if verbose
        println("\n", "-"^78)
        println("SC1-b realized vs predicted false discovery proportion: ",
                sc1b_met ? "WITHIN BAND AT EVERY LEVEL" : "ABOVE THE UPPER BAND AT $(length(offenders)) LEVEL(S)")
        for r in rows
            # ONE line per level carrying the level, the rule's predicted FDP, the decided fraction
            # as both a ratio and a percentage and the scope token -- emitted through the SAME
            # `_p14_headline_line` a batch object would print, never a second format string -- and
            # the MEASURED realized FDP appended to it, so no FDR number in this printout can be
            # read without its decided fraction.
            println("  ", _p14_headline_line(r.alpha, r.predicted_fdp, r.n_decided, n_eval),
                    " | MEASURED realized FDP = ", r.realized_fdp,
                    " (", r.n_false_discoveries, "/", r.k_star, ")",
                    " | band [", r.band_lo, ", ", r.band_hi, "]",
                    isempty(r.note) ? "" : " | " * r.note)
        end
        println()
        println("  full table:")
        println("    alpha     k*    decided(frac)   predicted     realized      band_hi      t*        cost")
        for r in rows
            println("    ", rpad(r.alpha, 8), "  ", rpad(r.k_star, 5), " ",
                    rpad(string(r.n_decided, "(", round(r.decided_fraction; digits = 4), ")"), 15),
                    " ", rpad(round(r.predicted_fdp; digits = 6), 12),
                    " ", rpad(round(r.realized_fdp; digits = 6), 12),
                    " ", rpad(r.band_hi === missing ? "n/a" : string(round(r.band_hi; digits = 6)), 12),
                    " ", rpad(round(r.t_star; digits = 5), 9),
                    " ", round(r.cost_ratio; digits = 5))
        end
        println()
        println("  composite-null decomposition of the ACCEPTED sets (REPORTED, not gated):")
        for r in rows
            println("    alpha ", rpad(r.alpha, 8), "  mean p_random = ",
                    rpad(round(r.mean_p_random; digits = 6), 12), "  mean p_exclusion = ",
                    rpad(round(r.mean_p_exclusion; digits = 6), 12), "  share random/exclusion = ",
                    round(r.share_random; digits = 4), " / ", round(r.share_exclusion; digits = 4))
        end
        println()
        println("  SC1-c PRIOR SENSITIVITY -- REPORTED, NOT GATED -- NO BAR IS ATTACHED TO THIS")
        println("  TABLE. Failing to report it is the failure; nothing here can be missed.")
        println("  frozen grid P14_PI_COLOC_GRID = $P14_PI_COLOC_GRID")
        println("    pi_coloc  alpha     k*     n_decided(frac)  predicted     realized")
        for r in pi_sensitivity
            println("    ", rpad(r.pi_coloc, 9), " ", rpad(r.alpha, 8), "  ", rpad(r.k_star, 6),
                    " ", rpad(string(r.n_decided, "(", round(r.decided_fraction; digits = 4), ")"), 16),
                    " ", rpad(round(r.predicted_fdp; digits = 6), 12),
                    " ", round(r.realized_fdp; digits = 6))
        end
        println()
        println("  THE HONESTY BLOCK IS PERSISTED INSIDE THE ARTIFACT, as the string keys")
        println("  fdr_claim, assumption_a2_statement, ordering_note, independence_note,")
        println("  prior_note, prior_atom_note, unit_note and amendment, so it travels with the")
        println("  numbers rather than depending on a report writer's memory.")
        println("  unit: ", _p14_fdr_claims().unit_note)
        println()
        println("  decided $n_decided / $n_eval   abstained $abstained   reasons: $reasons")
        println("  guarantee basis        = simulator-derived held-out draws; both the predicted")
        println("                           and the realized number come from the same forward")
        println("                           model the evidence net was trained on")
        println("  elapsed                = $(round(time() - t_start; digits = 1)) s")
        println("-"^78)
        reported || println("SMOKE MODE: the band is NOT asserted and these numbers are NOT a verdict.")
    end

    # --- 8. THE ASSERTION, LAST ----------------------------------------------------------------
    if reported
        @assert sc1b_met """
        SC1-b NOT MET. The realized false discovery proportion EXCEEDS the upper binomial 95 %
        band of the rule's own predicted value at $(length(offenders)) of
        $(length(P14_ALPHA_FDR_GRID)) pre-registered levels.

        Offending rows (each with its decided fraction, which may never be quoted separately):
        $(join(["  alpha = $(r.alpha): realized $(r.realized_fdp) > band_hi $(r.band_hi) " *
                "(predicted $(r.predicted_fdp), k* = $(r.k_star), decided $(r.n_decided)/$n_eval " *
                "= $(r.decided_fraction))" for r in offenders], "\n"))

        AN HONEST SHORTFALL IS A PHASE-14 FINDING. It is REPORTED, not repaired. Relaxing
        P14_ALPHA_FDR_GRID, widening the band, moving n_eval, redrawing or restratifying the
        evaluation set, reseeding, or re-scoping the criterion to a friendlier subset are ALL
        forbidden.

        P14_ITERATION_ALLOWANCE CANNOT BE SPENT HERE. Its one pre-declared trigger concerns
        split-conformal coverage falling below the SC1-d band, and 14-09 measured that coverage
        ABOVE the band, so the trigger did not fire and this is not the condition it names:

        $P14_ITERATION_TRIGGER

        The full report was persisted BEFORE this assertion and is intact at:
          $report_path
        """
    end

    return (alpha_table = rows, sc1b_met = sc1b_met,
            n_eval = n_eval, n_decided = n_decided,
            decided_fraction = n_decided / n_eval,
            reported = reported, report_path = report_path)
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so including this file can NEVER trigger the reported run. Run it deliberately:
#
#     julia --project=spike -t auto spike/p14/run_p14_fdr_check.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
