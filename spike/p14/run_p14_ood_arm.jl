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

# spike/p14/run_p14_ood_arm.jl --- GATED SC3-d: does abstention concentrate on OUT-OF-DISTRIBUTION
# data? Five matched arms -- one in-distribution and one per in-repo misspecification family at the
# strongest existing grid rung -- each scored through the SAME composed decision path, at the SAME
# conformal q-hat, against the SAME density null.
#
# =============================================================================================
# ANTI-SNOOPING BANNER. READ BEFORE CHANGING ANY NUMBER IN THIS FILE.
# =============================================================================================
# THE BAR THIS RUNNER SCORES AGAINST WAS COMMITTED TO `spike/p14/consts.jl` BEFORE ANY PHASE-14
# RESULT EXISTED AND BEFORE THIS FILE EXISTED. It is P14_OOD_MARGIN_FLOOR; the matched arm size is
# P14_N_OOD and the reserved stream counter is P14_OOD_COUNTER. Every one of them is READ from that
# frozen file BY NAME. Not one of them is a literal here, and the acceptance criteria grep the
# comment-stripped source for exactly that, so a value cannot be "temporarily" inlined and left
# behind.
#
# P14_OOD_MARGIN_FLOOR IS A JUDGEMENT CALL WITH NO DERIVATION. It is the fifth member of
# `P14_JUDGEMENT_CALL_BARS`, printed and persisted MECHANICALLY from that frozen tuple rather than
# by hand, so it cannot quietly stop being labelled one. It was frozen in a commit that precedes
# every result-producing Phase-14 commit and it may NEVER be relaxed after a result.
# `P14_ITERATION_ALLOWANCE = 1` cannot be spent on it: its ONE pre-declared trigger concerns
# split-conformal COVERAGE falling below the SC1-d band, and 14-09 measured that ABOVE the band, so
# the trigger did not fire.
#
# SC3-d IS A DIFFERENT CLAIM FROM SC3-c AND BOTH ARE NEEDED. SC3's "hard/OOD" is a slash, not a
# synonym. SC3-c (14-11) asks whether abstention tracks WOULD-HAVE-BEEN-WRONG; SC3-d asks whether
# it tracks THE DATA IS NOT WHAT THE MODEL WAS TRAINED ON. A layer can do either without the other.
#
# NO MISSPECIFICATION IS AUTHORED HERE. The four families are the in-repo `OOD_FAMILIES`
# (`spike/validation/ood.jl:499`) at the strongest rung of the EXISTING grid, consumed through
# `p14_misspec_pool`. Authoring the stressor and the bar in the same phase is tuning, and an
# acceptance criterion greps this file to prove no family is defined in it.
#
# A NOT-CHECKED ABSTENTION IS NOT DETECTOR PERFORMANCE, AND THE BREAKDOWN EXISTS TO SAY SO.
# Under D-06 an item whose OOD state is `:not_checked` ABSTAINS by default. An arm that abstained
# because nothing ever checked it has demonstrated nothing about OOD sensitivity, so `ood_fired`
# and `ood_not_checked` are counted SEPARATELY per arm and the printout states which drove the
# margin. That is also why this runner HARD-STOPS when the null is unavailable: unlike every other
# Phase-14 runner, SC3-d is meaningless without a live OOD channel.
#
# THE HEADLINE IS THE MINIMUM OVER FAMILIES, NOT THE MEAN. SC3-d is falsified if abstention is no
# more frequent under A strong misspecification, so a single weak family must stay visible rather
# than be averaged away by three strong ones. The per-family table is persisted beside it.
#
# THE ORDER OF OPERATIONS IS THE POINT: read the bar, build the arms, compute the statistics,
# PERSIST THE ARTIFACT, print the headline, then assert. The artifact is written before anything can
# throw, so a FAILING gate still leaves a complete, self-describing report on disk. The evidence of
# an honest failure must survive the failure.
#
# THE Q-HAT IS LOADED, NEVER RE-CALIBRATED. It is the quantile 14-09 cut the conformal sets at. A
# fresh quantile here would make this runner's abstention rates incomparable with 14-09's, and
# nothing in either artifact would say so. If that artifact is absent this runner FAILS LOUDLY and
# names 14-09.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY AND NOT A TEST. Including it does nothing; run it:
#
#     julia --project=spike -t auto spike/p14/run_p14_ood_arm.jl
#
# DECOUPLING: spike-local. `src/` is reached only READ-ONLY and transitively, and step 0 PROVES at
# run time that it is byte-unchanged AND that no untracked file has appeared under it. Nothing here
# writes to the gitignored Phase-13 pool; it is read-only for the whole phase.

using JLD2
using Statistics
using Dates

# --- Guarded includes, in dependency order --------------------------------------------------------
# `decide.jl` transitively pulls the pre-registration, the provenance loader, the posterior, the
# hedge, the fusion, the result type and the Phase-13 read surface. `pools.jl` carries the draw,
# the density null and the MATCHED misspecification arms.
#
# THE MISSPECIFICATION FAMILIES COME FIRST, AND THE ORDER IS LOAD-BEARING RATHER THAN TIDY.
# `OOD_FAMILIES` lives ONLY in `spike/validation/ood.jl`, but that file also carries spike-lane
# twins of `fit_ood_nulls`, `maha_score` and `roc_auc` whose SHIPPED counterparts in
# `src/amortized/ood.jl` arrive transitively through `decide.jl`. In Julia the LAST definition
# wins, so including the validation file AFTER `decide.jl` would silently swap the implementation
# the density null is fitted with -- and plans 14-09, 14-10 and 14-11 all fitted theirs with the
# SHIPPED one (their guards named `fit_ood_nulls` / `roc_auc`, which `src/` had already defined, so
# the validation file was never reached). Loading it FIRST leaves the shipped definitions on top,
# exactly as in 14-09, while making the four families available. The operating point is then PINNED
# to the one 14-09 recorded, so this ordering argument is proved at run time and not merely argued.
#
# The guard names `OOD_FAMILIES` itself -- the symbol this runner actually needs -- rather than a
# name `src/amortized/ood.jl` also defines, which is precisely the guard that would silently skip.
isdefined(@__MODULE__, :OOD_FAMILIES) ||
    include(joinpath(@__DIR__, "..", "validation", "ood.jl"))
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :P14_PROVENANCE_LOADED)   || include(joinpath(@__DIR__, "provenance.jl"))
isdefined(@__MODULE__, :P14_POSTERIOR_LOADED)    || include(joinpath(@__DIR__, "posterior.jl"))
isdefined(@__MODULE__, :P14_CONFORMAL_LOADED)    || include(joinpath(@__DIR__, "conformal.jl"))
isdefined(@__MODULE__, :P14_FUSE_LOADED)         || include(joinpath(@__DIR__, "fuse.jl"))
isdefined(@__MODULE__, :P14Result)               || include(joinpath(@__DIR__, "result.jl"))
isdefined(@__MODULE__, :P14_DECIDE_LOADED)       || include(joinpath(@__DIR__, "decide.jl"))
isdefined(@__MODULE__, :P14_POOLS_LOADED)        || include(joinpath(@__DIR__, "pools.jl"))

if !isdefined(@__MODULE__, :P14_OOD_ARM_RUNNER_LOADED)
    "The SC3-d artifact. Written BEFORE the assertion can throw."
    const P14_OOD_ARM_REPORT_PATH = joinpath(@__DIR__, "p14_ood_arm_report.jld2")

    """
    The 14-09 artifact carrying the q-hat this runner cuts every conformal set at.

    NAMED here rather than obtained by including `run_p14_conformal.jl`: that file defines its own
    `main`, and including one runner from another would put two zero-argument `main` methods in one
    namespace where the last include silently wins.
    """
    const P14_OOD_ARM_QHAT_SOURCE = joinpath(@__DIR__, "p14_conformal_report.jld2")

    "The keys this runner requires from the 14-09 artifact, so a schema drift fails BY NAME."
    const P14_OOD_ARM_REQUIRED_QHAT_KEYS = ("qhat", "coverage", "band_lower", "n_eval",
                                            "P14_ALPHA_CONFORMAL", "reported", "pool_provenance")

    """
    The IN-DISTRIBUTION arm's name, used as the first key of every per-arm table.

    A Symbol rather than a position, for the same reason every class read in this lane is by name:
    the margin is `rate(family) - rate(ID)` and getting that subtraction backwards would produce a
    sign-flipped number that still looks like a rate difference.
    """
    const P14_OOD_ARM_ID_KEY = :in_distribution

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else can
    # make this body skip.
    const P14_OOD_ARM_RUNNER_LOADED = true
end

# =============================================================================================
# Small local helpers
# =============================================================================================

"""
    _p14_ood_save_report(path, required; kwargs...) -> String

Atomically persist an artifact: write `path * ".tmp"`, REOPEN it read-only and integrity-check
every key in `required`, then `mv(...; force = true)`.

The `_p13_gate_save_report` idiom (`spike/p13/run_three_way_gate.jl:171-190`), carried per runner
exactly as `run_p14_conformal.jl`, `run_p14_fdr_check.jl` and `run_p14_riskcoverage.jl` carry it.
It is deliberately NOT imported from a sibling runner: importing would mean including that runner
and inheriting its `main`.

A crash mid-write leaves a discardable `.tmp` rather than a torn artifact a later reader would
happily believe.
"""
function _p14_ood_save_report(path, required; kwargs...)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        for k in required
            @assert haskey(f, k) "_p14_ood_save_report: integrity check failed -- $tmp is missing the required key `$k`"
        end
    end
    mv(tmp, path; force = true)
    return path
end

"""
    _p14_ood_load_qhat(path) -> NamedTuple

Load the split-conformal threshold PLAN 14-09 CALIBRATED, and prove the artifact is the one this
runner is entitled to read rather than merely a file with the right name.

Three checks, each catching something the others cannot: the file EXISTS (and if not, this fails
loudly naming 14-09 and never re-calibrates); every key read here is PRESENT by name; and the
stored q-hat is FINITE, because a non-finite conformal threshold silently admits every class or
none and both look like a legitimate regime.

RE-CALIBRATING HERE WOULD BE THE SILENT FAILURE. The abstention rates below are quoted beside
14-09's; a fresh quantile would put them on a different hedge, and no reader of either artifact
could tell.
"""
function _p14_ood_load_qhat(path::AbstractString)
    isfile(path) || error("""
        run_p14_ood_arm: the 14-09 conformal report is ABSENT at
            $path

        That artifact carries the split-conformal q-hat this runner cuts every conformal set at.
        It is produced by PLAN 14-09 (`spike/p14/run_p14_conformal.jl`). Run it first:

            julia --project=spike -t auto spike/p14/run_p14_conformal.jl

        THIS RUNNER WILL NOT RE-CALIBRATE THE HEDGE, silently or otherwise. A fresh quantile would
        make its abstention rates incomparable with the ones 14-09 reported, and nothing in either
        artifact would say so.""")

    d = JLD2.load(path)
    for k in P14_OOD_ARM_REQUIRED_QHAT_KEYS
        haskey(d, k) || error(
            "run_p14_ood_arm: the 14-09 conformal report at $path carries no `$k` key. The 14-09 " *
            "artifact schema moved under this runner; reading a q-hat out of a partially " *
            "understood artifact is how a number ends up meaning something other than its name.")
    end
    q = Float64(d["qhat"])
    isfinite(q) || error(
        "run_p14_ood_arm: the q-hat stored by 14-09 is $q, which is not finite. A non-finite " *
        "conformal threshold admits every class or none, and both look like a legitimate regime.")

    return (qhat = q,
            path = path,
            coverage = Float64(d["coverage"]),
            band_lower = Float64(d["band_lower"]),
            n_eval = Int(d["n_eval"]),
            alpha = Float64(d["P14_ALPHA_CONFORMAL"]),
            source_reported = d["reported"],
            pool_provenance = d["pool_provenance"],
            sha = p14_blob_sha(path))
end

"""
    _p14_ood_breakdown(results) -> NamedTuple

The count of every member of the CLOSED reason set `p14_abstain_reasons()`, plus the decided count,
over one arm.

ENUMERATED, NEVER COLLECTED. Counting whatever reasons turn up would silently drop a reason that
stopped being reachable and silently admit one that was never declared; both are exactly the kind
of drift the closed set exists to prevent. `:ood_fired` and `:ood_not_checked` are therefore
SEPARATE keys and are never summed into an "OOD" total anywhere in this file: an arm that abstained
because nothing ever checked it has demonstrated no OOD sensitivity at all (T-14-39, D-06).
"""
function _p14_ood_breakdown(results)
    reasons = p14_abstain_reasons()
    counts  = map(r -> count(x -> abstain_reason(x) === r, results), reasons)
    bd = NamedTuple{reasons}(counts)
    n_abstain = count(x -> decision(x) === :abstain, results)
    @assert sum(counts) == n_abstain "_p14_ood_breakdown: $(sum(counts)) reasons were counted over $(n_abstain) abstentions -- an abstention carrying a reason outside the closed set p14_abstain_reasons() would be invisible in this table"
    return bd
end

"""
    _p14_ood_cross_tab(rows, row_levels, cols, col_levels) -> Matrix{Int}

A dense contingency table over DECLARED levels, erroring on any pair outside them.

The same shape `run_p14_conformal.jl:275` carries, re-stated per runner for the reason every
helper in this lane is: importing it would mean including that runner and inheriting its `main`. A
cell counted into nothing would make the table sum to less than n with nothing saying why, which is
why the fall-through errors rather than dropping the pair.
"""
function _p14_ood_cross_tab(rows::AbstractVector{Symbol}, row_levels,
                            cols::AbstractVector{Symbol}, col_levels)
    length(rows) == length(cols) || throw(ArgumentError(
        "_p14_ood_cross_tab: the two vectors must be the same length; got $(length(rows)) and $(length(cols))."))
    M = zeros(Int, length(row_levels), length(col_levels))
    for i in eachindex(rows)
        r = findfirst(==(rows[i]), row_levels)
        c = findfirst(==(cols[i]), col_levels)
        (r === nothing || c === nothing) && error(
            "_p14_ood_cross_tab: the pair ($(rows[i]), $(cols[i])) is outside the declared levels " *
            "$(row_levels) x $(col_levels). A cell counted into nothing would make the table sum " *
            "to less than n with nothing saying why.")
        M[r, c] += 1
    end
    @assert sum(M) == length(rows) "_p14_ood_cross_tab: the table does not sum to n"
    return M
end

"""
    _p14_ood_arm_record(results) -> NamedTuple

One arm's reported record: `n`, `abstain_rate`, the trigger `breakdown`, the state and status
vectors it was computed from, and the two rate columns SC3-d's honesty rests on.

`ood_fired_rate` and `ood_not_checked_rate` are carried as first-class fields rather than left to
be divided out of the breakdown by a reader, because the whole question "was the margin driven by
the detector or by the detector never running?" is answered by comparing exactly those two columns
across arms.
"""
function _p14_ood_arm_record(results)
    n  = length(results)
    bd = _p14_ood_breakdown(results)
    statuses = Symbol[conformal_status(r) for r in results]
    states   = Symbol[ood_state(r) for r in results]
    n_abstain = count(x -> decision(x) === :abstain, results)
    return (n = n,
            n_abstain = n_abstain,
            abstain_rate = n_abstain / n,
            trigger_breakdown = bd,
            ood_fired_rate = bd.ood_fired / n,
            ood_not_checked_rate = bd.ood_not_checked / n,
            conformal_empty_rate = bd.conformal_empty / n,
            conformal_ambiguous_rate = bd.conformal_ambiguous / n,
            disagreement_and_ood_rate = bd.disagreement_and_ood / n,
            conformal_status = statuses,
            ood_state = states,
            decision = Symbol[decision(r) for r in results],
            set_size_counts = (singleton = count(==(:singleton), statuses),
                               ambiguous = count(==(:ambiguous), statuses),
                               empty     = count(==(:empty),     statuses)),
            ood_state_counts = (fired       = count(==(:fired),       states),
                                clear       = count(==(:clear),       states),
                                not_checked = count(==(:not_checked), states)))
end

# =============================================================================================
# The reported run
# =============================================================================================

"""
    main(; n_ood = P14_N_OOD, report_path = P14_OOD_ARM_REPORT_PATH, verbose = true)
        -> NamedTuple

The GATED SC3-d measurement: five matched arms, an abstention rate and a trigger breakdown for
each, a per-family margin table and the conservative headline.

`n_ood` defaults to the FROZEN per-arm size. Passing anything else puts the run in SMOKE mode: it
writes to a `_smoke` artifact path and prints a banner saying so, so a load check can never be
mistaken for, or quoted as, a verdict, and can never overwrite the reported artifact.
"""
function main(; n_ood::Integer = P14_N_OOD,
                report_path = P14_OOD_ARM_REPORT_PATH,
                verbose::Bool = true)

    t_start  = time()
    reported = (n_ood == P14_N_OOD)
    reported || (report_path = replace(report_path, ".jld2" => "_smoke.jld2"))

    # --- 0. THE DECOUPLING PROOF, AT RUN TIME (D-01 / D-02) ---------------------------------
    # Asserted WHILE the run happens, not only in review: a long reported run that started on a
    # clean tree and finished on a dirty one would be a repudiation hole (T-14-11). The shared lane
    # guard is called first; the two commands are then re-run here so this runner's own failure
    # message names the breach, and a THIRD check is added that neither the guard nor a diff can
    # make -- an untracked NEW file under `src/` is invisible to `git diff HEAD` and would be
    # exactly the productionization-by-stealth D-01 forbids.
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

    # --- 1. THE INHERITANCE (D-07) AND THE DENSITY NULL --------------------------------------
    # UNLIKE EVERY OTHER PHASE-14 RUNNER, A MISSING NULL IS A HARD STOP HERE. Elsewhere an
    # unavailable OOD channel is a reportable D-06 outcome: every item resolves to :not_checked and
    # abstains, and the coverage or FDR number is unaffected because it is a property of the
    # conformal set rather than of the decision. SC3-d is the one claim that IS about the detector.
    # With no null, both arms would abstain at rate 1 by default, the margin would be exactly zero,
    # and that zero would be a statement about this runner rather than about the layer.
    bundle = p14_load_bundle()
    ref    = p14_ood_reference()
    ood_in = p14_ood_input(ref)
    basis  = load_p13_basis()

    ref.available || error("""
        run_p14_ood_arm: THE OOD REFERENCE IS UNAVAILABLE, so SC3-d cannot be measured at all.

        reason:   $(ref.note)
        pool_dir: $(repr(ref.pool_dir))

        The density null is fitted on the PHASE-13 NET'S OWN training pool, resolved from the
        trained net's `h.meta.pool_dir`. That pool is ~51 MB of GITIGNORED bulk data, so it exists
        only in a working tree that has actually generated it. THIS PLAN MUST THEREFORE RUN ON THE
        MAIN WORKING TREE AND NOT IN A GIT WORKTREE (`execution_environment:
        main-working-tree-no-worktrees`): a fresh worktree has no pool, and a `git worktree remove
        --force` has already destroyed a 54 MB gitignored pool once in this milestone.

        This is a HARD STOP rather than the reportable :not_checked outcome the other Phase-14
        runners carry, because SC3-d is the one claim that is ABOUT the detector: with no null both
        arms abstain at rate 1 under D-06, the margin is exactly zero, and that zero would be a
        statement about this runner rather than about the decision layer.""")

    isfinite(ref.thr) || error(
        "run_p14_ood_arm: the fitted operating point is $(ref.thr), which is not finite. Under the " *
        "shipped semantics (`src/amortized/local_map.jl:107`) a non-finite threshold is NOT " *
        "CHECKED, so every item would abstain by default and the SC3-d margin would be zero for a " *
        "reason that has nothing to do with the detector.")

    qh = _p14_ood_load_qhat(P14_OOD_ARM_QHAT_SOURCE)

    # THE OPERATING POINT IS PINNED TO THE ONE 14-09 SCORED AGAINST, and this is the assertion that
    # PROVES the include-ordering argument at the head of this file rather than leaving it argued.
    # Both the density fit and the threshold are re-derived here from the same pool; if a different
    # implementation of `fit_ood_nulls` / `maha_score` had won the load order, the quantile would
    # move and the in-distribution abstain rate below would no longer be comparable with 14-09's --
    # silently, because a slightly different threshold is still a perfectly plausible threshold.
    @assert qh.pool_provenance.ood_thr == ref.thr """
    run_p14_ood_arm: the density-null operating point fitted here is $(ref.thr), while plan 14-09
    recorded $(qh.pool_provenance.ood_thr) for the same pool. The two runners' abstention rates are
    quoted beside each other, so they must be decided against the SAME operating point. A moved
    threshold means the null was fitted on a different pool, at a different quantile, on a
    different standardizer basis, or with a different implementation of `fit_ood_nulls` /
    `maha_score` winning the include order.
      this run: pool_dir = $(ref.pool_dir), q = $(ref.q), basis = $(ref.basis_provenance)
      14-09   : pool_dir = $(qh.pool_provenance.pool_dir), q = $(qh.pool_provenance.ood_id_quantile), basis = $(qh.pool_provenance.ood_basis_provenance)
    """

    # --- 2. THE LOCKED-THRESHOLD BANNER, BEFORE ANYTHING IS COMPUTED -------------------------
    families = collect(keys(OOD_FAMILIES))
    if verbose
        println("="^78)
        println("PHASE-14 GATED SC3-d: DOES ABSTENTION CONCENTRATE ON OUT-OF-DISTRIBUTION DATA?")
        println("="^78)
        reported || println("*** SMOKE MODE -- NOT A REPORTED RUN, THE GATE IS NOT ASSERTED ***")
        println("AMENDED CRITERIA: ", P14_AMENDMENT_NOTICE)
        println("GATED (SC3-d):")
        println("  P14_OOD_MARGIN_FLOOR = $P14_OOD_MARGIN_FLOOR   (consts.jl section E)")
        println("                         GATED. JUDGEMENT CALL -- NO DERIVATION.")
        println("  the judgement-call enumeration, printed MECHANICALLY from the frozen file so a")
        println("  bar cannot quietly stop being labelled one:")
        println("  P14_JUDGEMENT_CALL_BARS = $P14_JUDGEMENT_CALL_BARS")
        println("  (the other four are SC3-a/b/c's and belong to run_p14_riskcoverage.jl.)")
        println("MATCHED ARMS / STREAM:")
        println("  P14_N_OOD            = $P14_N_OOD   PER ARM (requested n_ood = $n_ood)")
        println("  P14_OOD_COUNTER      = $P14_OOD_COUNTER      (the reserved SC3-d stream)")
        println("  arms                 = 1 in-distribution + $(length(families)) misspecified")
        println("WHAT IS INJECTED (consumed, NEVER authored here):")
        println("  families (spike/validation/ood.jl, the in-repo positive controls): ",
                join(families, ", "))
        println("  level                = :strongest -> rung $OOD_GRID_LEVELS of the EXISTING grid")
        println("  No new misspecification was written for this phase, so the arm cannot have been")
        println("  tuned to produce a margin.")
        println("THE HEDGE IS LOADED, NEVER RE-CALIBRATED:")
        println("  q-hat                = $(qh.qhat)   (from $(basename(qh.path)), plan 14-09)")
        println("  14-09 coverage       = $(qh.coverage) against band $(qh.band_lower)")
        println("ITERATION DISCIPLINE:")
        println("  P14_ITERATION_ALLOWANCE = $P14_ITERATION_ALLOWANCE, and it CANNOT be spent here:")
        println(P14_ITERATION_TRIGGER)
        println("INHERITED (D-07, LOADED and four-way agreed, never re-derived):")
        println("  tau                  = $(bundle.tau)")
        println("  class prior (MEASURED) = $(bundle.prior)")
        println("OOD CHANNEL -- A NAMED LIMIT, RECORDED RATHER THAN PAPERED OVER:")
        println("  available            = $(ref.available)   thr = $(ref.thr)")
        println("  channels wired       = $(ref.channels_wired)")
        println("  channels NOT wired   = $(ref.channels_not_wired)")
        println("  basis                = $(ref.basis_provenance)   pool = $(ref.pool_dir)")
        println("  a `:clear` here means ONE channel says clear, not three.")
        println("D-06:")
        println("  an item whose OOD state is :not_checked ABSTAINS by default. `ood_fired` and")
        println("  `ood_not_checked` are counted SEPARATELY per arm, because an arm that abstained")
        println("  because nothing checked it has demonstrated NO OOD sensitivity.")
        println("="^78)
    end

    # --- 3. THE IN-DISTRIBUTION ARM -----------------------------------------------------------
    verbose && println("[1/6] drawing the UNSTRATIFIED in-distribution arm (n = $n_ood) …")
    id_items = p14_draw_pool(n_ood; counter = P14_OOD_COUNTER, basis = basis)
    id_masses = reported ? p14_assert_unstratified([it.class for it in id_items]) :
                           class_masses([it.class for it in id_items])
    verbose && println("      realized in-distribution class masses: $id_masses")

    # --- 4. THE FOUR MISSPECIFIED ARMS, MATCHED ITEM BY ITEM ----------------------------------
    # `p14_misspec_pool` returns BOTH halves of a matched pair built from ONE parameter draw per
    # index: identical theta, lambda and image size, differing only in which generator produced the
    # pixels. Its in-distribution half rides the SAME stream and the SAME indices as the arm drawn
    # above, so the two are asserted EQUAL summary by summary rather than assumed to be. If they
    # ever diverged, the reported margin would be partly a difference between two draws rather than
    # between two generators, and nothing would say so.
    arms_raw = Dict{Symbol, Any}()
    levels   = Dict{Symbol, Int}()
    for fam in families
        verbose && println("[2/6] building the matched `$fam` arm at the strongest rung (n = $n_ood) …")
        mp = p14_misspec_pool(n_ood; counter = P14_OOD_COUNTER, family = fam,
                              level = :strongest, basis = basis)
        @assert length(mp.ood) == length(id_items) "run_p14_ood_arm: the `$fam` arm holds $(length(mp.ood)) items against the in-distribution arm's $(length(id_items)). A rate difference between arms of unequal n is partly an artefact of sample size."
        @assert length(mp.id) == length(id_items) "run_p14_ood_arm: the `$fam` matched in-distribution half holds $(length(mp.id)) items against the drawn in-distribution arm's $(length(id_items))"
        for j in eachindex(id_items)
            @assert mp.id[j].Zs == id_items[j].Zs "run_p14_ood_arm: the `$fam` arm's own in-distribution half disagrees with the drawn in-distribution arm at item $j. The two arms would then differ by the DRAW as well as by the misspecification, and the reported margin would be partly a stream artefact."
            @assert mp.id[j].class === id_items[j].class "run_p14_ood_arm: the `$fam` arm's matched class label disagrees with the drawn in-distribution arm at item $j"
        end
        arms_raw[fam] = mp
        levels[fam]   = mp.level
    end
    @assert length(Set(values(levels))) == 1 "run_p14_ood_arm: the four families resolved to different grid rungs $(levels); `:strongest` must resolve to one rung of the existing grid for every family"

    # --- 5. THE DECISIONS, ONE COMPOSED PATH, ONE Q-HAT, ONE NULL -----------------------------
    # `allow_unchecked_ood = false` is the D-06 default and is passed EXPLICITLY, because the whole
    # SC3-d reading depends on it: with it true a :not_checked item would DECIDE rather than
    # abstain, and the not-checked column below -- the column that keeps a vacuous abstention from
    # being read as detector performance -- would be empty by construction.
    verbose && println("[3/6] deciding every item of every arm at the LOADED q-hat …")
    id_results = [p14_decide_one(bundle, it.Zs, it.Zc, ood_in, qh.qhat;
                                 lambda = it.lambda, idx = it.idx,
                                 allow_unchecked_ood = false) for it in id_items]
    id_record  = _p14_ood_arm_record(id_results)

    fam_records = Dict{Symbol, Any}()
    for fam in families
        mp = arms_raw[fam]
        rs = [p14_decide_one(bundle, it.Zs, it.Zc, ood_in, qh.qhat;
                             lambda = it.lambda, idx = it.idx,
                             allow_unchecked_ood = false) for it in mp.ood]
        fam_records[fam] = _p14_ood_arm_record(rs)
    end

    # --- 6. THE PER-FAMILY MARGINS AND THE CONSERVATIVE HEADLINE ------------------------------
    verbose && println("[4/6] computing the per-family margins …")
    margin_per_family = NamedTuple{Tuple(families)}(
        Tuple(fam_records[f].abstain_rate - id_record.abstain_rate for f in families))
    margin = minimum(values(margin_per_family))
    weakest_family = families[argmin([getproperty(margin_per_family, f) for f in families])]

    arms = merge(NamedTuple{(P14_OOD_ARM_ID_KEY,)}((id_record,)),
                 NamedTuple{Tuple(families)}(Tuple(fam_records[f] for f in families)))
    @assert length(arms) == length(families) + 1 "run_p14_ood_arm: the arm table holds $(length(arms)) entries for $(length(families)) families plus the in-distribution arm"
    for k in keys(arms)
        @assert getproperty(arms, k).n == length(id_items) "run_p14_ood_arm: arm `$k` holds $(getproperty(arms, k).n) items against $(length(id_items)); the arms are not matched in n"
    end

    return (arms = arms, margin = margin, margin_per_family = margin_per_family,
            weakest_family = weakest_family, id_record = id_record,
            qhat = qh.qhat, n_ood = Int(n_ood), reported = reported,
            report_path = report_path, elapsed_s = time() - t_start)
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so including this file can NEVER trigger the reported run. Run it deliberately:
#
#     julia --project=spike -t auto spike/p14/run_p14_ood_arm.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
