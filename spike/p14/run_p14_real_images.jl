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

# spike/p14/run_p14_real_images.jl --- SC1-f (D-03a): the decision layer, run on the ONLY real
# microscopy this project may touch.
#
# =============================================================================================
# THIS RUNNER PRODUCES AN ILLUSTRATION, NOT A COVERAGE CLAIM. READ THIS FIRST.
# =============================================================================================
# SIX TIFFS IN TWO CONDITIONS, WITH NO COLOCALIZATION GROUND-TRUTH LABELS, BOUNDS NOTHING
# TIGHTLY. No rate, fraction or percentage is computed over these six images anywhere in this
# file, `is_coverage_claim = false` and `illustration_note` are REQUIRED keys of the artifact, and
# the integrity check refuses to write a report that lacks either. A reader who quotes a number
# from here as a coverage number is quoting a statistic this runner deliberately never computed.
#
# READ-ONLY DISCIPLINE. test/test_images/ IS READ-ONLY INPUT. THIS FILE OPENS THE TIFFS AND
# NOTHING ELSE: IT NEVER WRITES, MOVES OR MODIFIES ANYTHING UNDER test/, AND NO test/ PATH APPEARS
# IN THIS PLAN'S files_modified LIST. The `git ls-files -s` digest of the fixture tree is recorded
# BEFORE and AFTER the run and asserted equal, and the sha256 of each of the six files is recorded
# beside its row, so a later change to a fixture is visible rather than silent.
#
# WHY THESE SIX AND NOT THE PROVENANCE MANIFEST'S PHYSICAL ANCHORS (D-03a). The manifest holds no
# unsealed physical ground truth: 30 of its 32 rows are `tier: simulated-secondary` -- the CBS
# benchmark, i.e. ANOTHER SIMULATOR -- and the only 2 `physical-primary` rows are the Phase-16
# sealed holdout. Its image bytes are also unfetched by design. Reading one of those two rows here
# would irreversibly burn Phase 16 on the very hypothesis Phase 16 exists to evaluate, and there is
# no second blind set. Phase 13 hit this exact problem and resolved it the same way
# (`spike/p13/real_images.jl:65-73`): THE SANCTIONED SUBSTITUTE IS THE SIX COMMITTED MICROSCOPY
# TIFFS UNDER test/test_images/, adopted precisely so that Phase 16's blind evaluation stays blind.
# Task 2 of this plan appends the METADATA-ONLY availability record that says all of that with
# numbers computed from the manifest rather than copied from prose.
#
# THE WITHDRAWN FIGURE APPEARS NOWHERE. D-03's corpus-derived miscoverage floor -- the ratio one
# over thirty-one, roughly three per cent -- is WITHDRAWN, not merely superseded: it was computed
# from the wrong n and on rows that are not physical truth. The numeral is banned from Phase-14
# executable code by `spike/test/test_p14_decoupling.jl`, and it is expressed in prose here so the
# ban has nothing to catch.
#
# A PREDICTION, RECORDED AS A PREDICTION AND NOT AS AN INPUT (Assumptions Log A6). 13-REPORT limit
# D records both fixtures scoring a summary density of 703.30 against their own threshold of
# 167.54 -- a factor of 4.198 -- so under D-05/D-06 the correct output for both specimens is
# ABSTAIN. That is a result worth reporting: on the only real microscopy this project holds, the
# decision layer abstains, and that is the DESIGNED behaviour rather than a failure. The prediction
# is persisted BESIDE the observed decisions with `prediction_held`, and if the images do NOT
# abstain that is recorded loudly as a finding about the OOD wiring -- nothing is adjusted here.
#
# THE Q-HAT IS LOADED, NEVER RE-CALIBRATED. It is the quantile 14-09 cut its conformal sets at. A
# fresh quantile here would put these decisions on a different hedge from every other Phase-14
# number, and nothing in either artifact would say so. If that artifact is absent this runner FAILS
# LOUDLY and names 14-09.
#
# THE OPERATIVE CHANNEL PAIR IS PASSED EXPLICITLY AT EVERY CALL SITE, never inherited from a
# keyword default. That is the Phase-13 correction recorded at `spike/p13/real_images.jl:59-63`:
# `P13_REAL_CHANNEL_PAIR` was amended from c1/c2 to the c2/c3 green/red pair because c1 is the
# DAPI/Hoechst nuclear counterstain, so any pair containing it measures counterstain-versus-protein
# overlap and NOT colocalization. All three channels of both conditions are opened and digested --
# that is what makes this a six-TIFF record -- but only the c2/c3 pair enters a decision, and every
# row says which of the two it is.
#
# THERE IS NO GATE HERE. SC1-f's failure conditions are all structural rather than numeric: a
# coverage number was emitted, the sealed holdout was read, the manifest was reported as gone
# rather than unfetched, or the six-TIFF result was presented as a coverage claim. Each is
# prevented by construction and asserted by `spike/test/test_p14_decoupling.jl`, so this runner
# persists and prints, and asserts no threshold.
#
# THE ORDER OF OPERATIONS IS THE POINT: read the inheritance, decide, PERSIST THE ARTIFACT, then
# print. The artifact is written before anything can throw, so even a run that dies in the printout
# leaves a complete, self-describing report on disk.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY AND NOT A TEST. Including it does nothing; run it:
#
#     julia --project=spike -t auto spike/p14/run_p14_real_images.jl
#
# DECOUPLING: spike-local. `src/` is reached only READ-ONLY and transitively, and step 0 PROVES at
# run time that it is byte-unchanged and that no untracked file has appeared under it.

using JLD2
using Statistics
using Dates
using SHA

# --- Guarded includes, in dependency order --------------------------------------------------------
# THE ORDER IS 14-09's, COPIED RATHER THAN REDERIVED. `decide.jl` pulls the SHIPPED
# `src/amortized/ood.jl` (hence `fit_ood_nulls`, `maha_score`, `roc_auc`) before `pools.jl` is
# reached, so `pools.jl`'s own guard on `fit_ood_nulls` is already satisfied and the spike-lane
# twins in `spike/validation/ood.jl` are never loaded. In Julia the last definition wins, so this
# ordering is what puts these decisions on the SAME density-null implementation plans 14-09, 14-10
# and 14-11 used. The real-image ingestion comes last: it is Phase 13's frozen loader, reused
# verbatim rather than re-authored, so there is exactly one real-image read path in this project.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :P14_PROVENANCE_LOADED)   || include(joinpath(@__DIR__, "provenance.jl"))
isdefined(@__MODULE__, :P14_POSTERIOR_LOADED)    || include(joinpath(@__DIR__, "posterior.jl"))
isdefined(@__MODULE__, :p14_bayes_fdr)           || include(joinpath(@__DIR__, "fdr.jl"))
isdefined(@__MODULE__, :P14_CONFORMAL_LOADED)    || include(joinpath(@__DIR__, "conformal.jl"))
isdefined(@__MODULE__, :P14_FUSE_LOADED)         || include(joinpath(@__DIR__, "fuse.jl"))
isdefined(@__MODULE__, :P14Result)               || include(joinpath(@__DIR__, "result.jl"))
isdefined(@__MODULE__, :P14_DECIDE_LOADED)       || include(joinpath(@__DIR__, "decide.jl"))
isdefined(@__MODULE__, :P14_POOLS_LOADED)        || include(joinpath(@__DIR__, "pools.jl"))
isdefined(@__MODULE__, :load_real) ||
    include(joinpath(@__DIR__, "..", "p13", "real_images.jl"))

if !isdefined(@__MODULE__, :P14_REAL_RUNNER_LOADED)
    "The SC1-f artifact. Written BEFORE anything is printed and BEFORE anything can throw."
    const P14_REAL_REPORT_PATH = joinpath(@__DIR__, "p14_real_images_report.jld2")

    """
    The 14-09 artifact carrying the q-hat every conformal set below is cut at.

    NAMED here rather than obtained by including `run_p14_conformal.jl`: that file defines its own
    `main`, and including one runner from another would put two zero-argument `main` methods in one
    namespace where the last include silently wins.
    """
    const P14_REAL_QHAT_SOURCE = joinpath(@__DIR__, "p14_conformal_report.jld2")

    "The keys this runner requires from the 14-09 artifact, so a schema drift fails BY NAME."
    const P14_REAL_REQUIRED_QHAT_KEYS = ("qhat", "coverage", "band_lower", "n_eval",
                                         "P14_ALPHA_CONFORMAL", "reported", "pool_provenance")

    """
    The repo-relative root of the D-03a substrate, spelled as ONE literal.

    It is a literal on purpose. `spike/test/test_p14_decoupling.jl` asserts the sealed-holdout
    negatives POSITIVELY -- it requires some Phase-14 source to name this path -- precisely so that
    "no Phase-14 file reads a sealed image" cannot be satisfied vacuously by a phase that reads no
    real images at all. Removing this string would silently turn three live assertions into three
    tautologies.
    """
    const P14_REAL_SUBSTRATE_ROOT = "test/test_images"

    """
    All three committed channels, opened per condition.

    The DECISION rides `P13_REAL_CHANNEL_PAIR` alone (c2/c3, the two target proteins). c1 is the
    DAPI/Hoechst nuclear counterstain and enters no comparison. Every channel is nevertheless
    opened and digested, because SC1-f's unit of record is THE SIX COMMITTED TIFFS and a file that
    was never opened cannot carry a digest that proves it was not modified.
    """
    const P14_REAL_ALL_CHANNELS = (1, 2, 3)

    """
    The two decision units: each condition once as the sample, with the other as its control.

    Both orientations are run because the pairing is not symmetric -- the three-way rule reads
    `rho_sample` and `rho_sample - rho_control` -- and because Phase 13's real arm ran exactly
    these two (`spike/p13/run_p13_realimage.jl:729-730`). Naming them here rather than deriving
    them from a loop over conditions keeps the sample/control roles visible at the point of use.
    """
    const P14_REAL_UNITS = ((name = :positive_as_sample,
                             sample = P13_REAL_SAMPLE, control = P13_REAL_CONTROL),
                            (name = :negative_as_sample,
                             sample = P13_REAL_CONTROL, control = P13_REAL_SAMPLE))

    """
    The prediction, recorded as a PREDICTION and never as an input (Assumptions Log A6).

    Stated in `spike/p14/run_p14_real_images.jl` before the run and persisted verbatim beside the
    observed decisions, so a reader can see that the expectation was written down first rather
    than fitted afterwards.
    """
    const P14_REAL_PREDICTED_BEHAVIOUR =
        "ABSTAIN on both specimens (13-REPORT limit D: density 703.30 vs threshold 167.54, 4.198x)"

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else can
    # make this body skip.
    const P14_REAL_RUNNER_LOADED = true
end

# =============================================================================================
# Small local helpers
# =============================================================================================

"""
    _p14_real_save_report(path, required; kwargs...) -> String

Atomically persist an artifact: write `path * ".tmp"`, REOPEN it read-only and integrity-check
every key in `required`, then `mv(...; force = true)`.

The `_p13_gate_save_report` idiom (`spike/p13/run_three_way_gate.jl:171-190`), carried per runner
exactly as `run_p14_conformal.jl`, `run_p14_fdr_check.jl`, `run_p14_riskcoverage.jl` and
`run_p14_ood_arm.jl` carry it. It is deliberately NOT imported from a sibling runner: importing
would mean including that runner and inheriting its `main`.

`is_coverage_claim` and `illustration_note` are among the required keys at every call site, so a
report that failed to say what it is not is never written at all.
"""
function _p14_real_save_report(path, required; kwargs...)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        for k in required
            @assert haskey(f, k) "_p14_real_save_report: integrity check failed -- $tmp is missing the required key `$k`"
        end
    end
    mv(tmp, path; force = true)
    return path
end

"""
    _p14_real_load_qhat(path) -> NamedTuple

Load the split-conformal threshold PLAN 14-09 CALIBRATED, and prove the artifact is the one this
runner is entitled to read rather than merely a file with the right name.

Three checks, each catching something the others cannot: the file EXISTS (and if not, this fails
loudly naming 14-09 and never re-calibrates); every key read here is PRESENT by name; and the
stored q-hat is FINITE, because a non-finite conformal threshold silently admits every class or
none and both look like a legitimate regime.
"""
function _p14_real_load_qhat(path::AbstractString)
    isfile(path) || error("""
        run_p14_real_images: the 14-09 conformal report is ABSENT at
            $path

        That artifact carries the split-conformal q-hat every set below is cut at. It is produced
        by PLAN 14-09 (`spike/p14/run_p14_conformal.jl`). Run it first:

            julia --project=spike -t auto spike/p14/run_p14_conformal.jl

        THIS RUNNER WILL NOT RE-CALIBRATE THE HEDGE, silently or otherwise. A fresh quantile would
        put these six decisions on a different hedge from every other Phase-14 number, and nothing
        in either artifact would say so.""")

    d = JLD2.load(path)
    for k in P14_REAL_REQUIRED_QHAT_KEYS
        haskey(d, k) || error(
            "run_p14_real_images: the 14-09 conformal report at $path carries no `$k` key. The " *
            "14-09 artifact schema moved under this runner; reading a q-hat out of a partially " *
            "understood artifact is how a number ends up meaning something other than its name.")
    end
    q = Float64(d["qhat"])
    isfinite(q) || error(
        "run_p14_real_images: the q-hat stored by 14-09 is $q, which is not finite. A non-finite " *
        "conformal threshold admits every class or none, and both look like a legitimate regime.")

    return (qhat = q, path = path,
            coverage = Float64(d["coverage"]), band_lower = Float64(d["band_lower"]),
            n_eval = Int(d["n_eval"]), alpha = Float64(d["P14_ALPHA_CONFORMAL"]),
            source_reported = d["reported"], pool_provenance = d["pool_provenance"],
            sha = p14_blob_sha(path))
end

"""
    _p14_real_costes_idx(i) -> Int

The Costes block-permutation stream index for decision unit `i`, derived from the RESERVED SC1-f
counter rather than from a bare position.

`P14_REAL_COUNTER` is the counter the pre-registration set aside for this arm. This runner draws no
pool and therefore consumes no Philox stream at that counter directly; the ONE random stream it
does consume is the Costes null, keyed by `costes_rng(P14_DEV_SEED, idx)`. Folding the reserved
counter into that index is what keeps the reservation meaningful: the resulting indices lie far
above any pool item index a sibling runner could pass, so no two Phase-14 Costes nulls can share a
key by accident.
"""
_p14_real_costes_idx(i::Integer) = 1000 * Int(P14_REAL_COUNTER) + Int(i)

"""
    _p14_real_digest(path) -> String

The sha256 hex digest of one committed fixture, recorded beside its row.

Read-only, and the point of it: the tree-level `git ls-files -s` digest proves nothing under
`test/` was written during the run, while these six per-file digests make a change to any single
fixture visible in the artifact itself, in a later checkout, without git.
"""
_p14_real_digest(path::AbstractString) = bytes2hex(SHA.sha256(read(path)))

"""
    _p14_real_row(r) -> NamedTuple

The decision columns of one unit, flattened into plain data.

The `ThreeWayClass` enum the classical channel returns is converted to a `Symbol` and the conformal
set to a `Vector{Symbol}` so every row of the persisted table has the SAME concrete type: a table
whose rows differ in type deserialises as `Vector{Any}` and a later reader cannot tell a schema
drift from a legal absence. The cross-method NamedTuple is otherwise carried whole -- `basis` and
`costes_seed` included -- because a comparison whose scale was not recorded is one nobody can
reconstruct later.
"""
function _p14_real_row(r::P14Result)
    cm = cross_method(r)
    cm_flat = merge(cm, (classical_call = cm.classical_call === nothing ? :not_computed :
                                          Symbol(cm.classical_call),))
    return (decision          = decision(r),
            abstain_reason    = abstain_reason(r) === nothing ? :none : abstain_reason(r),
            ood_state         = ood_state(r),
            ood_score         = Float64(r.three_way.ood.score),
            conformal_status  = conformal_status(r),
            conformal_set     = collect(Symbol, conformal_set(r)),
            class_posterior   = class_posterior(r),
            null_posterior    = Float64(null_posterior(r)),
            null_split        = null_split(r),
            cross_method      = cm_flat)
end

# =============================================================================================
# The reported run
# =============================================================================================

"""
    main(; channels = P13_REAL_CHANNEL_PAIR, report_path = P14_REAL_REPORT_PATH, verbose = true)
        -> NamedTuple

The SC1-f illustration: the six committed microscopy TIFFs, opened read-only, decided through the
composed Phase-14 path at the LOADED 14-09 q-hat, with the pre-stated ABSTAIN prediction persisted
beside the observed outcome.

`channels` defaults to the AMENDED operative pair and is passed EXPLICITLY to every call site.
Passing anything else puts the run in SMOKE mode: it writes to a `_smoke` artifact path and prints
a banner saying so, so a load check can never be mistaken for the reported record and can never
overwrite it.

**No rate, fraction or percentage is computed over these six images.** See the banner.
"""
function main(; channels = P13_REAL_CHANNEL_PAIR,
                report_path = P14_REAL_REPORT_PATH,
                verbose::Bool = true)

    t_start  = time()
    reported = Tuple(channels) == Tuple(P13_REAL_CHANNEL_PAIR)
    reported || (report_path = replace(report_path, ".jld2" => "_smoke.jld2"))
    chans = collect(Int, channels)

    # --- 0. THE DECOUPLING PROOF, AT RUN TIME (D-01 / D-02) ---------------------------------
    # Asserted WHILE the run happens, not only in review: a long reported run that started on a
    # clean tree and finished on a dirty one would be a repudiation hole (T-14-11). A THIRD check
    # is added that neither the guard nor a diff can make -- an untracked NEW file under `src/` is
    # invisible to `git diff HEAD` and would be exactly the productionization-by-stealth D-01
    # forbids -- and a FOURTH that is specific to this runner: `test/` must be clean before the
    # fixtures are opened, so a dirty fixture tree is caught BEFORE it is read rather than after.
    lane = p14_assert_lane_clean(; verbose = verbose)
    if lane.checked
        src_untracked  = readchomp(Cmd(`git status --porcelain -- src`; dir = P14_REPO_ROOT))
        test_untracked = readchomp(Cmd(`git status --porcelain -- test`; dir = P14_REPO_ROOT))
        @assert isempty(src_untracked) "D-01 decoupling breach: an UNTRACKED file appeared under src/ ($(src_untracked)) -- a diff against HEAD cannot see it, and shipping the decision layer by stealth is precisely what D-01 forbids"
        @assert isempty(test_untracked) "T-14-45: the fixture tree under test/ is not clean before the run ($(test_untracked)). The six TIFFs are READ-ONLY INPUT; a dirty tree means this is not the shipped substrate."
    end

    # --- 1. THE STRUCTURAL READ-ONLY DIGEST, BEFORE ANY FILE IS OPENED -----------------------
    digest_before = verify_real_readonly_digest()

    # --- 2. THE INHERITANCE (D-07), THE DENSITY NULL AND THE LOADED HEDGE --------------------
    # A missing null is NOT a hard stop here. Unlike SC3-d, SC1-f is not a claim about the
    # detector: with no null every specimen resolves to :not_checked and ABSTAINS by default, which
    # is a reportable D-06 outcome and is recorded as such rather than crashing the arm.
    bundle = p14_load_bundle()
    ref    = p14_ood_reference()
    ood_in = p14_ood_input(ref)
    basis  = load_p13_basis()
    qh     = _p14_real_load_qhat(P14_REAL_QHAT_SOURCE)

    # --- 3. THE BANNER, BEFORE ANYTHING IS COMPUTED ------------------------------------------
    if verbose
        println("="^78)
        println("PHASE-14 SC1-f (D-03a): THE DECISION LAYER ON THE SIX COMMITTED MICROSCOPY TIFFS")
        println("="^78)
        reported || println("*** SMOKE MODE -- NOT THE REPORTED RECORD ***")
        println("THIS IS AN ILLUSTRATION, NOT A COVERAGE CLAIM.")
        println("  Six TIFFs in two conditions, with NO colocalization ground-truth labels, is a")
        println("  very small real-data check. It bounds nothing tightly. No rate, fraction or")
        println("  percentage is computed over these six images anywhere in this runner.")
        println("READ-ONLY DISCIPLINE:")
        println("  substrate            = $P14_REAL_SUBSTRATE_ROOT   (committed fixtures)")
        println("  this runner opens the TIFFs and nothing else; it never writes, moves or")
        println("  modifies anything under test/, and no test/ path is in files_modified.")
        println("  git ls-files digest BEFORE = $digest_before")
        println("AMENDED CRITERIA: ", P14_AMENDMENT_NOTICE)
        println("THE OPERATIVE PAIR, PASSED EXPLICITLY AT EVERY CALL SITE:")
        println("  P13_REAL_CHANNEL_PAIR = $P13_REAL_CHANNEL_PAIR  (c2/c3 = green/red, the two")
        println("    target proteins). c1 is the DAPI/Hoechst nuclear counterstain and enters no")
        println("    comparison; all three channels are still opened and digested.")
        println("  requested channels    = $chans")
        println("THE HEDGE IS LOADED, NEVER RE-CALIBRATED:")
        println("  q-hat                = $(qh.qhat)   (from $(basename(qh.path)), plan 14-09)")
        println("  14-09 coverage       = $(qh.coverage) against band $(qh.band_lower)")
        println("INHERITED (D-07, LOADED and four-way agreed, never re-derived):")
        println("  tau                  = $(bundle.tau)")
        println("  class prior (MEASURED) = $(bundle.prior)")
        println("OOD CHANNEL -- A NAMED LIMIT, RECORDED RATHER THAN PAPERED OVER:")
        println("  available            = $(ref.available)   thr = $(ref.thr)")
        println("  channels wired       = $(ref.channels_wired)")
        println("  channels NOT wired   = $(ref.channels_not_wired)")
        println("  a `:clear` here means ONE channel says clear, not three.")
        println("THE PREDICTION, STATED BEFORE THE RUN (Assumptions Log A6):")
        println("  $P14_REAL_PREDICTED_BEHAVIOUR")
        println("  If the specimens do NOT abstain, that is recorded LOUDLY as a finding about the")
        println("  OOD wiring. Nothing is adjusted here either way.")
        println("="^78)
    end

    # --- 4. THE SIX FILES, ASSERTED PRESENT AND DIGESTED BEFORE ANYTHING IS DECIDED ----------
    # `real_tif` validates the condition and the channel index against the pre-registration and
    # raises naming the expected absolute path if a fixture is missing, so the failure says WHICH
    # file rather than surfacing a bare SystemError from inside the image loader.
    verbose && println("[1/5] asserting and digesting the six committed TIFFs …")
    files = NamedTuple[]
    for cond in P13_REAL_CONDITIONS, ch in P14_REAL_ALL_CHANNELS
        p = real_tif(cond, ch)
        push!(files, (condition = String(cond), channel = Int(ch),
                      rel_path = P14_REAL_SUBSTRATE_ROOT * "/" * cond * "/$(cond)_c$(ch).tif",
                      sha256 = _p14_real_digest(p),
                      in_operative_pair = ch in Tuple(channels)))
    end
    @assert length(files) == 6 "run_p14_real_images: the substrate is the SIX committed TIFFs; $(length(files)) were enumerated"
    @assert length(unique(f.sha256 for f in files)) == 6 "run_p14_real_images: two of the six fixtures have the same sha256; the six-row record would then be describing fewer than six distinct files"

    # --- 5. THE TWO DECISION UNITS ------------------------------------------------------------
    # `decide_coloc` is the src-SHAPED per-pair entry point -- a bundle plus two images and the two
    # channels to compare, in; an AbstractColocResult subtype, out -- and "src-shaped" is NOT "in
    # src": no byte of src/ is edited by this phase and the decision layer is NOT shipped in v2.0.
    # `allow_unchecked_ood = false` is the D-06 default and is passed EXPLICITLY, because the whole
    # reading depends on it: with it true a :not_checked specimen would DECIDE rather than abstain.
    verbose && println("[2/5] loading both conditions and deciding at the LOADED q-hat …")
    results = P14Result[]
    unit_rows = NamedTuple[]
    for (i, u) in enumerate(P14_REAL_UNITS)
        img = load_real(u.sample;  channels = P14_REAL_ALL_CHANNELS)
        ctl = load_real(u.control; channels = P14_REAL_ALL_CHANNELS)
        r = decide_coloc(bundle, img, ctl, chans;
                         ood_nulls = ood_in, qhat = qh.qhat,
                         lambda = P13_PHASE11_REFERENCE_LAMBDA,
                         idx = _p14_real_costes_idx(i),
                         basis = basis,
                         allow_unchecked_ood = false,
                         meta = (unit = u.name, substrate = :real_committed_tiffs))
        push!(results, r)
        push!(unit_rows, merge((unit = u.name, sample = u.sample, control = u.control,
                                channels = chans,
                                ood_threshold = ref.thr === nothing ? NaN : Float64(ref.thr),
                                costes_idx = _p14_real_costes_idx(i)),
                               _p14_real_row(r)))
        verbose && println("      $(u.name): $(decision(r))  " *
                           "(ood_state = $(ood_state(r)), reason = $(abstain_reason(r)))")
    end

    # --- 6. THE SIX-ROW RECORD ----------------------------------------------------------------
    # THE UNIT OF RECORD IS THE FILE; THE UNIT OF DECISION IS THE ACQUISITION. A decision is a
    # property of a two-channel acquisition compared against a control, not of a single TIFF, so
    # each file row carries the decision of the unit in which ITS CONDITION was the sample, and
    # `in_operative_pair` says whether that channel actually entered the comparison. Stating that
    # rather than quietly duplicating a number is the difference between a table a reader can use
    # and one that invites a wrong reading.
    unit_of = Dict(String(u.sample) => unit_rows[i] for (i, u) in enumerate(P14_REAL_UNITS))
    images = [merge(f, (decision_unit = unit_of[f.condition].unit,
                        decision_unit_control = unit_of[f.condition].control,
                        decision = unit_of[f.condition].decision,
                        abstain_reason = unit_of[f.condition].abstain_reason,
                        ood_state = unit_of[f.condition].ood_state,
                        ood_score = unit_of[f.condition].ood_score,
                        ood_threshold = unit_of[f.condition].ood_threshold,
                        conformal_status = unit_of[f.condition].conformal_status,
                        conformal_set = unit_of[f.condition].conformal_set,
                        class_posterior = unit_of[f.condition].class_posterior,
                        null_posterior = unit_of[f.condition].null_posterior,
                        cross_method = unit_of[f.condition].cross_method)) for f in files]
    @assert length(images) == 6 "run_p14_real_images: the images table must hold exactly the six committed TIFFs"

    # --- 7. THE PREDICTION AND ITS OUTCOME, SIDE BY SIDE --------------------------------------
    observed = Symbol[decision(r) for r in results]
    prediction_held = all(d -> d === :abstain, observed)

    # --- 8. THE READ-ONLY PROOF, RE-TAKEN -----------------------------------------------------
    digest_after = verify_real_readonly_digest()

    # --- 9. PERSIST THE REPORT, BEFORE ANY TABLE IS PRINTED -----------------------------------
    verbose && println("[3/5] persisting the SC1-f report (BEFORE anything is printed) …")
    consts_path = joinpath(@__DIR__, "consts.jl")
    _p14_real_save_report(report_path,
        ("images", "decision_units", "predicted_behaviour", "prediction_held", "observed_decisions",
         "is_coverage_claim", "illustration_note", "readonly_digest_before",
         "readonly_digest_after", "substrate_root", "provenance");
        # --- WHAT THIS IS NOT. Required keys, so a report that failed to say so is never written.
        is_coverage_claim = false,
        illustration_note =
            "Six TIFFs in two conditions, with no colocalization ground-truth labels, is a very " *
            "small real-data check. It bounds nothing tightly and is reported as an " *
            "illustration, not a coverage claim. No rate, fraction or percentage is computed " *
            "over these six images anywhere in the runner that produced this artifact.",
        no_rate_computed = true,
        n_images = length(images),
        n_decision_units = length(unit_rows),
        decision_unit_note =
            "THE UNIT OF RECORD IS THE FILE; THE UNIT OF DECISION IS THE ACQUISITION. A decision " *
            "is a property of a two-channel acquisition compared against a control, so each of " *
            "the six file rows carries the decision of the unit in which its condition was the " *
            "SAMPLE, and `in_operative_pair` records whether that channel entered the comparison " *
            "at all. c1 is the DAPI/Hoechst nuclear counterstain and enters none of them.",
        # --- THE SIX FILES AND THE TWO DECISIONS ---
        images = images,
        decision_units = unit_rows,
        observed_decisions = observed,
        substrate_root = P14_REAL_SUBSTRATE_ROOT,
        conditions = collect(String, P13_REAL_CONDITIONS),
        channels_used = chans,
        P13_REAL_CHANNEL_PAIR = collect(Int, P13_REAL_CHANNEL_PAIR),
        channel_pair_note =
            "P13_REAL_CHANNEL_PAIR was AMENDED from c1/c2 to the c2/c3 green/red pair by " *
            "13-D15-AMENDMENT.md: test/runtests records the fixture channels as blue/green/red, " *
            "so c1 is the nuclear counterstain and any pair containing it measures " *
            "counterstain-versus-protein overlap and NOT colocalization. The pair is passed " *
            "EXPLICITLY at every call site so the operative pair is visible at the point of use " *
            "rather than inherited from a keyword default.",
        naming_correction_note =
            "THE FOLDER NAMES positive/ AND negative/ ARE THE ORIGINAL PACKAGE'S BIOLOGICAL TEST " *
            "CONDITIONS, NOT COLOCALIZATION LABELS. The `negative` pair is a positively " *
            "correlated pair, not an anti-correlated one; nothing here may treat it as an " *
            "exclusion example (spike/p13/real_images.jl, correction (b)).",
        # --- THE PREDICTION, WRITTEN DOWN FIRST ---
        predicted_behaviour = P14_REAL_PREDICTED_BEHAVIOUR,
        prediction_held = prediction_held,
        prediction_note =
            "The prediction was stated in the runner BEFORE the run and is persisted verbatim " *
            "beside the observed decisions. If it did NOT hold, that is a finding about the OOD " *
            "wiring and is reported as one; no threshold, constant or channel was adjusted here " *
            "in either case.",
        # --- THE READ-ONLY PROOF ---
        readonly_digest_before = digest_before,
        readonly_digest_after = digest_after,
        readonly_digest_equal = digest_before == digest_after,
        read_only = P13_REAL_READ_ONLY,
        file_sha256 = [f.sha256 for f in files],
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
        ood_available = ref.available,
        ood_threshold = ref.thr === nothing ? NaN : Float64(ref.thr),
        ood_id_quantile = ref.q,
        ood_basis_provenance = ref.basis_provenance,
        channels_wired = ref.channels_wired,
        channels_not_wired = ref.channels_not_wired,
        shipped_bundle_read = false,
        pool_provenance = p14_pool_provenance(ref),
        allow_unchecked_ood = false,
        # --- what the numbers rest on ---
        guarantee_basis = :simulator_derived_illustrated_on_six_committed_tiffs,
        named_limits = p14_named_limits(),
        named_limits_count = length(p14_named_limits()),
        lambda_used = P13_PHASE11_REFERENCE_LAMBDA,
        lambda_note =
            "Phase 11 closed NEGATIVE on inferring registration from the 8x8 summary, so lambda " *
            "is a STATED REFERENCE VALUE and not an estimate recovered from these images.",
        # --- provenance ---
        provenance = p14_provenance_record(bundle.prov),
        class_prior_used = bundle.prior,
        tau = bundle.tau,
        grid = bundle.grid,
        master_seed = UInt64(P14_DEV_SEED),
        salt = UInt64(P14_SALT),
        P14_REAL_COUNTER = Int(P14_REAL_COUNTER),
        costes_seed = UInt64(P14_DEV_SEED),
        consts_sha = p14_consts_sha(),
        consts_git_blob_sha = p14_blob_sha(consts_path),
        reported = reported,
        julia_version = string(VERSION),
        nthreads = Threads.nthreads(),
        elapsed_s = time() - t_start,
        generated = string(Dates.now(Dates.UTC)) * "Z")
    verbose && println("      persisted (before anything is printed) -> $report_path")

    # --- 10. THE SIX-ROW TABLE, PRINTED AFTER THE ARTIFACT IS ON DISK -------------------------
    if verbose
        println("[4/5] the six committed TIFFs")
        println("\n", "-"^78)
        println("SC1-f: THE DECISION LAYER ON THE SIX COMMITTED MICROSCOPY TIFFS (D-03a)")
        println("  substrate = $P14_REAL_SUBSTRATE_ROOT   operative pair = $chans")
        println()
        for row in images
            println("  $(rpad(row.rel_path, 42)) c$(row.channel)  " *
                    "in_pair = $(row.in_operative_pair)")
            println("    sha256      $(row.sha256)")
            println("    decision    $(row.decision)   reason = $(row.abstain_reason)   " *
                    "unit = $(row.decision_unit)")
            println("    ood         state = $(row.ood_state)   score = " *
                    "$(round(row.ood_score; digits = 4))   threshold = " *
                    "$(round(row.ood_threshold; digits = 4))")
            println("    conformal   status = $(row.conformal_status)   set = $(row.conformal_set)")
            println("    posterior   $(row.class_posterior)   null = " *
                    "$(round(row.null_posterior; digits = 6))")
            println("    classical   call = $(row.cross_method.classical_call)   costes_p = " *
                    "$(row.cross_method.costes_p_sample)   basis = $(row.cross_method.basis)")
        end
        println()
        println("[5/5] the prediction and its outcome, side by side")
        println("  PREDICTED : $P14_REAL_PREDICTED_BEHAVIOUR")
        println("  OBSERVED  : $observed")
        println("  HELD      : $prediction_held")
        if !prediction_held
            println("  THE PREDICTION DID NOT HOLD. That is a FINDING about the OOD wiring, not a")
            println("  licence to adjust anything: the specimens score far above the fitted")
            println("  operating point, so a non-abstaining outcome means the state the decision")
            println("  read was not the state the detector computed. Chase the wiring; do not")
            println("  move a threshold.")
        else
            println("  On the only real microscopy this project holds, the decision layer")
            println("  ABSTAINS -- and that is the DESIGNED behaviour (D-05/D-06), not a failure.")
        end
        println()
        println("  READ-ONLY PROOF")
        println("    git ls-files -s digest before = $digest_before")
        println("    git ls-files -s digest after  = $digest_after")
        println("    equal = $(digest_before == digest_after)")
        println()
        println("  THIS IS AN ILLUSTRATION, NOT A COVERAGE CLAIM. No rate, fraction or percentage")
        println("  was computed over these six images. Six TIFFs in two conditions, with no")
        println("  colocalization ground-truth labels, bounds nothing tightly.")
        println("    elapsed = $(round(time() - t_start; digits = 1)) s")
        println("-"^78)
        reported || println("SMOKE MODE: these numbers are NOT the reported record.")
    end

    # --- 11. THE STRUCTURAL ASSERTIONS, LAST --------------------------------------------------
    # PERSISTED FIRST, PRINTED SECOND, ASSERTED THIRD. There is no THRESHOLD to assert -- SC1-f has
    # none -- but the read-only discipline is a property this run must not have broken, and it is
    # checked after the artifact is safely on disk so the evidence survives the failure.
    @assert digest_before == digest_after """
    T-14-45: the `git ls-files -s test/test_images` digest CHANGED across this run.

      before : $digest_before
      after  : $digest_after

    test/test_images/ is READ-ONLY INPUT. This runner opens the TIFFs and nothing else; a moved
    digest means a fixture's content, mode or name changed while the run was in progress. The
    report was persisted BEFORE this assertion and is intact at:
      $report_path
    """

    return (images = images, decision_units = unit_rows, results = results,
            observed_decisions = observed, prediction_held = prediction_held,
            is_coverage_claim = false, qhat = qh.qhat, reported = reported,
            readonly_digest = digest_after, report_path = report_path,
            elapsed_s = time() - t_start)
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so including this file can NEVER trigger the reported run. Run it deliberately:
#
#     julia --project=spike -t auto spike/p14/run_p14_real_images.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
