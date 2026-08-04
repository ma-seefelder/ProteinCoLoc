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

"""
    _p14_ood_claims(ref, families, level) -> NamedTuple

The honesty block, written INTO the artifact as string keys so it travels WITH the numbers and
cannot be dropped by a report writer.

Every one of these is a sentence a reader of a bare SC3-d margin would otherwise have to be told by
someone who remembered to tell them. A `.jld2` that carries the number and not the scope is a
repudiation hole (T-14-39, T-14-27, T-14-40, T-14-02), which is why each of these is a REQUIRED key
of the save rather than an optional extra: an artifact missing one is never written at all.
"""
function _p14_ood_claims(ref::NamedTuple, families, level::Integer)
    fam_list = join(families, ", ")
    return (
        bar_note =
            "P14_OOD_MARGIN_FLOOR = $(P14_OOD_MARGIN_FLOOR) IS A JUDGEMENT CALL WITH NO " *
            "DERIVATION. Nothing derives an absolute abstain-rate margin of that size from the " *
            "detector's operating characteristic, from the misspecification grid, or from " *
            "anything else; it is a number a human chose as the level at which `abstention " *
            "concentrates on OOD data' would count as demonstrated. It is the fifth member of " *
            "P14_JUDGEMENT_CALL_BARS and is printed and persisted from that frozen tuple " *
            "MECHANICALLY, so it cannot quietly stop being labelled one. It was frozen in " *
            "`spike/p14/consts.jl` in a commit that precedes every result-producing Phase-14 " *
            "commit, and it MAY NEVER BE RELAXED AFTER A RESULT. This project has recorded four " *
            "separate cases of an underived bar measuring something other than what it named, so " *
            "the licensing standard for amending one is EVIDENCE, INDEPENDENT OF THE MACHINERY " *
            "UNDER AUDIT, that the bar measures something other than what it names. A shortfall " *
            "is a PHASE-14 FINDING and is REPORTED, not re-tuned. P14_ITERATION_ALLOWANCE has " *
            "exactly ONE pre-declared trigger and it is the split-conformal coverage band, which " *
            "14-09 measured ABOVE the band, so the allowance is unspent and is not available here.",

        channels_note =
            "ONLY THE SUMMARY-DENSITY CHANNEL IS WIRED IN THIS LANE. channels_wired = " *
            "$(ref.channels_wired), channels_not_wired = $(ref.channels_not_wired). The shipped " *
            "flag is an OR-fusion over three channels; here a `:clear` means ONE channel says " *
            "clear, not three, and a family the density channel happens to be blind to will show " *
            "a small margin THAT THE UNWIRED CHANNELS MIGHT HAVE CAUGHT. That is recorded rather " *
            "than papered over, and it cuts BOTH ways: it is not an excuse for a shortfall, " *
            "because the layer as measured is the layer as it exists, but a reader comparing this " *
            "number to a three-channel detector's is not comparing like with like.",

        basis_note =
            "THE NULL IS FITTED ON THE PHASE-13 NET'S OWN POOL ($(ref.pool_dir)), at the " *
            "pre-registered in-distribution quantile $(ref.q), basis_provenance = " *
            "$(ref.basis_provenance). PITFALL 4 is the reason: the shipped bundle at " *
            "artifacts/amended_v2/grid_8/ood_nulls_8.jld2 rides the PHASE-7 standardizer while " *
            "this net rides PHASE 11's, inherited and never re-fit, so scoring a " *
            "Phase-11-standardized summary against a Phase-7-fitted Mahalanobis null would " *
            "produce a number with NO INTERPRETATION -- not a wrong number, an uninterpretable " *
            "one. THE SHIPPED BUNDLE WAS NOT READ BY THIS RUNNER AT ALL. If any runner ever reads " *
            "it, it is a labelled COMPARISON REFERENCE and the two thresholds are NEVER fused. " *
            "The operating point used here is additionally PINNED, by an executable assertion, to " *
            "the one plan 14-09 recorded, so the two runners' abstain rates are decided against " *
            "the same threshold rather than merely believed to be.",

        misspecification_note =
            "NO MISSPECIFICATION WAS AUTHORED FOR THIS PHASE. The families are the in-repo " *
            "OOD_FAMILIES ($fam_list) defined in `spike/validation/ood.jl`, and " *
            "the level is `:strongest`, which resolves to rung $level of the EXISTING grid. Each " *
            "is a pre-registered positive control producing structure the shared-latent " *
            "smooth-field simulator cannot produce. Authoring the stressor and the bar in the " *
            "same phase is tuning, and an acceptance criterion greps this runner to prove no " *
            "family is defined in it (T-14-40). The arms are MATCHED item by item -- identical " *
            "theta, lambda, image size and stream position, differing only in which generator " *
            "produced the pixels -- and the matched in-distribution half of every family arm is " *
            "asserted EQUAL, summary by summary, to the separately drawn in-distribution arm.",

        d06_note =
            "AN ABSTENTION DRIVEN BY `ood_not_checked` IS NOT DETECTOR PERFORMANCE. Under D-06 an " *
            "item whose OOD state is `:not_checked` ABSTAINS by default, and that default exists " *
            "because a `false` flag from an unfitted null means NOT CHECKED rather than IN " *
            "DISTRIBUTION (`src/amortized/local_map.jl:70-72`). An arm that abstained because " *
            "nothing ever checked it has demonstrated NOTHING about OOD sensitivity. So " *
            "`ood_fired` and `ood_not_checked` are counted SEPARATELY in every arm's " *
            "trigger_breakdown, are never summed into an `OOD' total anywhere, and the headline " *
            "records `margin_driver` per family. This runner additionally HARD-STOPS when the " *
            "null is unavailable, which is the one case in which every arm would abstain at rate " *
            "one for a reason that has nothing to do with the detector (T-14-39).",

        crosstab_note =
            "THE EMPTY CONFORMAL SET AND THE DENSITY FLAG ARE TWO INDEPENDENT MISSPECIFICATION " *
            "SIGNALS, AND THE CROSS-TAB IS WHERE THAT IS VISIBLE. They arrive through completely " *
            "different channels: the density null depends on the simulator being right about " *
            "densities, whereas an empty conformal set depends only on EXCHANGEABILITY -- which " *
            "arguably makes it the most valuable single output of the conformal hedge, because it " *
            "is the one signal in this phase that does not depend on the simulator being right " *
            "about densities. BOTH READINGS ARE REPORTED AND NEITHER IS EDITORIALISED INTO THE " *
            "OTHER: where the two AGREE (an empty set on an item the density null also flagged) " *
            "that is corroboration worth reporting, and where they DISAGREE (an empty set on an " *
            "item the density null called clear, or vice versa) that is the more interesting " *
            "finding, because one channel is seeing something the other cannot. The 3x3 table " *
            "keeps `:empty` distinct from `:ambiguous` throughout (T-14-22).",

        gated_reading_note =
            "THE GATED READING IS THE MINIMUM OVER FAMILIES, NOT THE MEAN, AND THE PER-FAMILY " *
            "TABLE IS PERSISTED BESIDE IT. SC3-d is falsified if abstention is no more frequent " *
            "under A strong, deliberately injected misspecification than in-distribution, so a " *
            "single weak family must stay VISIBLE rather than be averaged away by three strong " *
            "ones (T-14-41). `margin` is therefore min over families of " *
            "(abstain_rate(family) - abstain_rate(in-distribution)), and `margin_per_family` " *
            "carries every one of them. A report quoting a mean would be quoting a different " *
            "statistic under this statistic's name.",

        circularity_note =
            "THE IN-DISTRIBUTION ARM IS A WELL-SPECIFIED-REGIME ARM (D-03 / D-03a). It is drawn " *
            "from the SAME simulator the evidence net was trained on and the density null was " *
            "fitted on, so its abstain rate is a statement about internal consistency and NOT " *
            "about real microscopy. The misspecified arms are that same simulator perturbed by " *
            "the in-repo positive controls, which is a bound on nothing physical: they say the " *
            "detector separates the simulator from a deliberately corrupted simulator. The " *
            "real-data bounding evidence D-03 intended is ABSENT, NOT MERELY LOOSE.",

        amendment = P14_AMENDMENT_NOTICE,
    )
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

    # --- 7. THE EMPTY-SET x OOD-STATE CROSS-TAB, PER ARM --------------------------------------
    # Two INDEPENDENT misspecification signals, tabulated against each other rather than merged.
    # `:empty` is kept distinct from `:ambiguous` throughout (T-14-22): collapsing them would
    # destroy exactly the channel this table exists to expose.
    verbose && println("[5/6] cross-tabulating the conformal status against the OOD state …")
    status_levels = collect(P14_CONFORMAL_STATUSES)
    ood_levels    = collect(P14_OOD_STATES)
    crosstab = NamedTuple{keys(arms)}(Tuple(
        _p14_ood_cross_tab(getproperty(arms, k).conformal_status, status_levels,
                           getproperty(arms, k).ood_state, ood_levels) for k in keys(arms)))
    for k in keys(crosstab)
        @assert sum(getproperty(crosstab, k)) == n_ood "run_p14_ood_arm: the `$k` cross-tab sums to $(sum(getproperty(crosstab, k))) rather than $n_ood; a cell counted into nothing would make the table understate one of the two channels"
    end

    # The two readings, computed rather than asserted, and BOTH reported. `empty_row` is the index
    # of `:empty`; `fired_col` / `clear_col` index the density verdict. Agreement is an empty set on
    # an item the null also flagged; disagreement is an empty set the null called clear, or a
    # flagged item whose conformal set was a perfectly ordinary singleton.
    empty_row = findfirst(==(:empty),  status_levels)
    fired_col = findfirst(==(:fired),  ood_levels)
    clear_col = findfirst(==(:clear),  ood_levels)
    agree_disagree = NamedTuple{keys(crosstab)}(Tuple(
        let M = getproperty(crosstab, k)
            (empty_and_fired  = M[empty_row, fired_col],
             empty_but_clear  = M[empty_row, clear_col],
             fired_but_not_empty = sum(M[:, fired_col]) - M[empty_row, fired_col],
             n_empty = sum(M[empty_row, :]),
             n_fired = sum(M[:, fired_col]))
        end for k in keys(crosstab)))

    # WHICH CHANNEL DROVE THE MARGIN, per family, computed from the two columns D-06 keeps apart.
    # A margin driven by `ood_not_checked` would be an arm abstaining because nothing checked it,
    # which demonstrates no OOD sensitivity at all (T-14-39).
    margin_driver = NamedTuple{Tuple(families)}(Tuple(
        let d_fired = fam_records[f].ood_fired_rate - id_record.ood_fired_rate,
            d_unchecked = fam_records[f].ood_not_checked_rate - id_record.ood_not_checked_rate,
            d_empty = fam_records[f].conformal_empty_rate - id_record.conformal_empty_rate,
            d_ambig = fam_records[f].conformal_ambiguous_rate - id_record.conformal_ambiguous_rate,
            d_disag = fam_records[f].disagreement_and_ood_rate - id_record.disagreement_and_ood_rate
            deltas = (ood_fired = d_fired, ood_not_checked = d_unchecked,
                      conformal_empty = d_empty, conformal_ambiguous = d_ambig,
                      disagreement_and_ood = d_disag)
            ks = keys(deltas)
            (deltas = deltas,
             dominant = ks[argmax([getproperty(deltas, k) for k in ks])],
             driven_by_not_checked = d_unchecked > d_fired)
        end for f in families))
    any_not_checked = any(k -> getproperty(arms, k).trigger_breakdown.ood_not_checked > 0,
                          keys(arms))

    margin_meets_floor = margin >= P14_OOD_MARGIN_FLOOR

    # --- 8. PERSIST THE REPORT, BEFORE THE HEADLINE AND BEFORE THE ASSERTION ------------------
    verbose && println("[6/6] persisting the SC3-d report (BEFORE any verdict) …")
    consts_path = joinpath(@__DIR__, "consts.jl")
    _p14_ood_save_report(report_path,
        ("arms", "margin", "margin_per_family", "margin_floor", "bar_is_judgement_call",
         "bar_note", "empty_set_crosstab", "channels_wired", "channels_not_wired", "basis_note",
         "misspecification_note", "d06_note", "amendment", "crosstab_note", "channels_note",
         "gated_reading_note", "circularity_note", "margin_driver", "sc3d_met",
         "abstain_rate_in_distribution", "abstain_rate_per_family", "provenance");
        # --- the honesty block, splatted in so it CANNOT be dropped: every one of its keys is
        #     also a REQUIRED key of the integrity check above, so an artifact missing one is
        #     never written at all.
        _p14_ood_claims(ref, families, first(values(levels)))...,
        # --- THE MEASUREMENT (SC3-d, GATED on the MINIMUM over families) ---
        arms = arms,
        margin = margin,
        margin_per_family = margin_per_family,
        weakest_family = weakest_family,
        gated_reading = :minimum_over_families,
        sc3d_met = margin_meets_floor,
        abstain_rate_in_distribution = id_record.abstain_rate,
        abstain_rate_per_family = NamedTuple{Tuple(families)}(
            Tuple(fam_records[f].abstain_rate for f in families)),
        # --- WHICH CHANNEL FIRED, kept apart rather than summed (D-06, T-14-39) ---
        margin_driver = margin_driver,
        any_not_checked = any_not_checked,
        ood_fired_rate_per_arm = NamedTuple{keys(arms)}(
            Tuple(getproperty(arms, k).ood_fired_rate for k in keys(arms))),
        ood_not_checked_rate_per_arm = NamedTuple{keys(arms)}(
            Tuple(getproperty(arms, k).ood_not_checked_rate for k in keys(arms))),
        trigger_breakdown_per_arm = NamedTuple{keys(arms)}(
            Tuple(getproperty(arms, k).trigger_breakdown for k in keys(arms))),
        abstain_reasons = collect(p14_abstain_reasons()),
        # --- THE TWO INDEPENDENT MISSPECIFICATION SIGNALS, CROSS-TABULATED (T-14-22) ---
        empty_set_crosstab = crosstab,
        empty_set_crosstab_status_levels = status_levels,
        empty_set_crosstab_ood_levels = ood_levels,
        empty_set_agreement = agree_disagree,
        set_size_counts_per_arm = NamedTuple{keys(arms)}(
            Tuple(getproperty(arms, k).set_size_counts for k in keys(arms))),
        ood_state_counts_per_arm = NamedTuple{keys(arms)}(
            Tuple(getproperty(arms, k).ood_state_counts for k in keys(arms))),
        # --- THE BAR, written INTO the artifact and labelled MECHANICALLY ---
        margin_floor = Float64(P14_OOD_MARGIN_FLOOR),
        bar_is_judgement_call = true,
        P14_OOD_MARGIN_FLOOR = P14_OOD_MARGIN_FLOOR,
        P14_JUDGEMENT_CALL_BARS = collect(P14_JUDGEMENT_CALL_BARS),
        bar_is_in_judgement_call_enumeration =
            :P14_OOD_MARGIN_FLOOR in P14_JUDGEMENT_CALL_BARS,
        P14_N_OOD = P14_N_OOD,
        P14_OOD_COUNTER = Int(P14_OOD_COUNTER),
        P14_ITERATION_ALLOWANCE = P14_ITERATION_ALLOWANCE,
        P14_ITERATION_TRIGGER = P14_ITERATION_TRIGGER,
        iteration_trigger_fired = false,
        iteration_allowance_applies_here = false,
        # --- THE ARMS AND WHAT MADE THEM ---
        n_ood = Int(n_ood),
        n_arms = length(arms),
        arms_matched_n = true,
        families = collect(families),
        family_level = first(values(levels)),
        family_level_is_strongest_existing_rung = true,
        ood_grid_levels = OOD_GRID_LEVELS,
        id_class_masses = id_masses,
        allow_unchecked_ood = false,
        # --- THE HEDGE, LOADED not re-calibrated ---
        qhat = qh.qhat,
        qhat_source = qh.path,
        qhat_source_sha = qh.sha,
        qhat_source_coverage = qh.coverage,
        qhat_source_band_lower = qh.band_lower,
        qhat_source_alpha = qh.alpha,
        qhat_source_reported = qh.source_reported,
        qhat_recalibrated_here = false,
        # --- THE DETECTOR ---
        channels_wired = ref.channels_wired,
        channels_not_wired = ref.channels_not_wired,
        ood_threshold = ref.thr,
        ood_id_quantile = ref.q,
        ood_basis_provenance = ref.basis_provenance,
        ood_score_quantiles = ref.score_quantiles,
        ood_threshold_pinned_to_14_09 = true,
        shipped_bundle_read = false,
        pool_provenance = p14_pool_provenance(ref),
        # --- what the numbers rest on ---
        guarantee_basis = :simulator_derived_matched_arms,
        named_limits = p14_named_limits(),
        named_limits_count = length(p14_named_limits()),
        # --- provenance ---
        provenance = p14_provenance_record(bundle.prov),
        class_prior_used = bundle.prior,
        tau = bundle.tau,
        grid = bundle.grid,
        master_seed = UInt64(P14_DEV_SEED),
        salt = UInt64(P14_SALT),
        consts_sha = p14_consts_sha(),
        consts_git_blob_sha = p14_blob_sha(consts_path),
        reported = reported,
        julia_version = string(VERSION),
        nthreads = Threads.nthreads(),
        elapsed_s = time() - t_start,
        generated = string(Dates.now(Dates.UTC)) * "Z")
    verbose && println("      persisted (before any verdict) -> $report_path")

    # --- 9. THE HEADLINE, PRINTED BEFORE ANY ASSERTION CAN THROW ------------------------------
    if verbose
        println("\n", "-"^78)
        println("SC3-d: ABSTENTION UNDER STRONG MISSPECIFICATION vs IN-DISTRIBUTION")
        println("  matched n per arm      = $n_ood        arms = $(length(arms))")
        println("  in-distribution abstain rate = $(id_record.abstain_rate)")
        println()
        println("  PER-ARM ABSTAIN RATE AND TRIGGER BREAKDOWN")
        println("  (ood_fired and ood_not_checked are SEPARATE columns and are never summed)")
        for k in keys(arms)
            a = getproperty(arms, k)
            println("    $(rpad(string(k), 16)) rate = $(rpad(round(a.abstain_rate; digits = 4), 8)) " *
                    "n_abstain = $(a.n_abstain)")
            println("      $(a.trigger_breakdown)")
            println("      set sizes $(a.set_size_counts)   ood states $(a.ood_state_counts)")
        end
        println()
        println("  PER-FAMILY MARGIN = abstain_rate(family) - abstain_rate(in-distribution)")
        for f in families
            m  = getproperty(margin_per_family, f)
            dr = getproperty(margin_driver, f)
            println("    $(rpad(string(f), 16)) margin = $(rpad(round(m; digits = 4), 8)) " *
                    "$(m >= P14_OOD_MARGIN_FLOOR ? ">=" : "<") floor   " *
                    "dominant channel = $(dr.dominant)   " *
                    "driven by NOT-CHECKED = $(dr.driven_by_not_checked)")
            println("      deltas $(dr.deltas)")
        end
        println()
        println("  WHICH CHANNEL DROVE THE MARGIN")
        println("    any item in any arm resolved :not_checked = $any_not_checked")
        if any_not_checked
            println("    SOME ABSTENTIONS WERE NOT-CHECKED ABSTENTIONS. An arm that abstains")
            println("    because the detector never ran has demonstrated NO OOD sensitivity, so")
            println("    the ood_not_checked column above must be read before the margin (D-06).")
        else
            println("    NO item in any arm resolved :not_checked, so no part of any margin is a")
            println("    default abstention: the density detector ran on every item of every arm,")
            println("    and every OOD-triggered abstention above is a FIRED detector (D-06).")
        end
        println()
        println("  EMPTY CONFORMAL SET x OOD STATE -- TWO INDEPENDENT MISSPECIFICATION SIGNALS")
        println("    rows = $(status_levels)   cols = $(ood_levels)")
        for k in keys(crosstab)
            println("    $(k):")
            M = getproperty(crosstab, k)
            for (i, s) in enumerate(status_levels)
                println("      $(rpad(string(s), 10)) $(M[i, :])")
            end
            g = getproperty(agree_disagree, k)
            println("      AGREEMENT   : empty AND fired = $(g.empty_and_fired)")
            println("      DISAGREEMENT: empty but CLEAR = $(g.empty_but_clear)   " *
                    "fired but NOT empty = $(g.fired_but_not_empty)")
        end
        println("    BOTH READINGS ARE REPORTED AND NEITHER IS EDITORIALISED INTO THE OTHER.")
        println("    Agreement is corroboration through two channels that share no assumption --")
        println("    the density null needs the simulator to be right about densities, the empty")
        println("    conformal set needs only exchangeability. Disagreement is the MORE")
        println("    INTERESTING finding: one channel is seeing something the other cannot.")
        println()
        println("  GATED READING: MINIMUM over families (never the mean -- a weak family must stay")
        println("  visible rather than be averaged away by three strong ones).")
        println("    weakest family        = $weakest_family")
        println("    SC3-d margin          = $margin  " *
                "$(margin_meets_floor ? ">=" : "<") P14_OOD_MARGIN_FLOOR = $P14_OOD_MARGIN_FLOOR" *
                "   -> $(margin_meets_floor ? "MET" : "NOT MET")")
        println("    THE FLOOR IS A JUDGEMENT CALL WITH NO DERIVATION, frozen in consts.jl before")
        println("    this runner existed, and listed in P14_JUDGEMENT_CALL_BARS =")
        println("    $P14_JUDGEMENT_CALL_BARS")
        println()
        println("  NAMED LIMITS CARRIED BY THIS NUMBER")
        println("    channels wired = $(ref.channels_wired), NOT wired = $(ref.channels_not_wired)")
        println("    a family the density channel is blind to shows a small margin here that the")
        println("    unwired noise / posterior-predictive channels might have caught; that is")
        println("    RECORDED, and it is not an excuse -- the layer as measured is the layer as it")
        println("    exists.")
        println("    the misspecification families are the in-repo OOD_FAMILIES at rung " *
                "$(first(values(levels))); no new one was authored, so the arm cannot have been")
        println("    tuned to produce a margin.")
        println("    elapsed = $(round(time() - t_start; digits = 1)) s")
        println("-"^78)
        reported || println("SMOKE MODE: no gate was asserted and these numbers are NOT a verdict.")
    end

    # --- 10. THE ASSERTION, LAST ---------------------------------------------------------------
    # PERSISTED FIRST, PRINTED SECOND, ASSERTED THIRD. The artifact is already on disk, so a
    # failing bar below leaves complete evidence of the failure behind it (T-14-31).
    if reported
        @assert margin_meets_floor """
        SC3-d NOT MET. The MINIMUM per-family abstain-rate margin is $margin against the frozen
        floor P14_OOD_MARGIN_FLOOR = $P14_OOD_MARGIN_FLOOR.

          in-distribution abstain rate : $(id_record.abstain_rate)
          per-family margins           : $margin_per_family
          weakest family               : $weakest_family
          margin drivers               : $(NamedTuple{Tuple(families)}(Tuple(getproperty(margin_driver, f).dominant for f in families)))
          any :not_checked item        : $any_not_checked
          channels wired               : $(ref.channels_wired)   NOT wired: $(ref.channels_not_wired)

        AN HONEST SHORTFALL IS A PHASE-14 FINDING, NEVER A LICENCE TO ACT. The floor is a
        JUDGEMENT CALL WITH NO DERIVATION, and that is a reason to REPORT it as such -- not a
        licence to relax it after seeing this number, nor to add a family, nor to change the grid
        rung, nor to wire a further OOD channel in order to clear a bar. Every one of those moves
        would be authoring the experiment after seeing its result.

        P14_ITERATION_ALLOWANCE has exactly ONE pre-declared trigger and it is the split-conformal
        coverage band, not this bar:

        $P14_ITERATION_TRIGGER

        The gated reading is the MINIMUM over families, deliberately, so that a single family the
        wired channel is blind to stays visible rather than being averaged away. Read
        `margin_per_family`, `margin_driver` and `empty_set_crosstab` in the artifact before
        concluding anything about the layer as a whole.

        The full report was persisted BEFORE this assertion and is intact at:
          $report_path
        """
    end

    return (arms = arms, margin = margin, margin_per_family = margin_per_family,
            margin_driver = margin_driver, weakest_family = weakest_family,
            empty_set_crosstab = crosstab, empty_set_agreement = agree_disagree,
            any_not_checked = any_not_checked, sc3d_met = margin_meets_floor,
            id_record = id_record,
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
