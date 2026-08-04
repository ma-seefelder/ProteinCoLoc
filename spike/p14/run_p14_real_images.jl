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
#
# THE AVAILABILITY RECORD IS METADATA ONLY, AND IT IS COMPUTED. This runner reads the tier, split
# and byte columns of the committed provenance manifest -- a text file whose whole purpose is
# provenance -- and NOTHING ELSE. It opens no image row, lists no directory of image bytes, and
# neither calls nor names the sole reachable accessor for the sealed rows. That is the narrower
# property `spike/test/test_p14_decoupling.jl` grants the files on its metadata allowlist, and
# being on that allowlist is a DECISION rather than a fix. Every count in the record -- the tier
# breakdown, the split partition, the byte accounting -- is derived from the file, so a manifest
# that changed would change the record instead of silently outliving it. The record says
# UNFETCHED BY DESIGN rather than missing, because the bytes were never expected to be in a
# checkout, and it says the D-03 bound is ABSENT NOT MERELY LOOSE, because `loose` would invite a
# reader to discount a weak bound where in fact there is none.
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
# THE MANIFEST IS A COMMITTED TEXT FILE AND IS READ AS METADATA ONLY. `CSV` and `DataFrames` are
# both among the frozen sixteen, so this adds no dependency and forces no `Pkg.resolve` (D-02).
# A hand-rolled comma split would be WRONG here as well as unnecessary: two rows of the manifest
# carry quoted fields containing commas, and splitting on the delimiter shears them.
using CSV
using DataFrames

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

    """
    The provenance manifest, read as METADATA ONLY.

    This runner is on the `_P14_CORPUS_METADATA_ALLOWED` allowlist of
    `spike/test/test_p14_decoupling.jl`, and that allowlist grants a NARROWER property, not a pass:
    metadata yes, the sealed image directory NEVER -- neither spelled as a path nor reached through
    a constant nor opened by the sole reachable accessor. Adding a file to that allowlist is a
    DECISION, not a fix. The manifest itself is a committed text file whose entire purpose is
    provenance; reading its tier, split and byte columns opens no image.
    """
    const P14_MANIFEST_PATH = joinpath(P14_REPO_ROOT, "corpus", "manifest.csv")

    """
    The ignore file whose rules are the EVIDENCE that the image bytes are unfetched BY DESIGN
    rather than lost.

    The matching lines and their line numbers are read at RUN TIME rather than transcribed, because
    a transcribed line number goes stale the first time an unrelated rule is inserted above it --
    and this one already had: 14-CONTEXT D-03a cites a line range that no longer holds.
    """
    const P14_GITIGNORE_PATH = joinpath(P14_REPO_ROOT, ".gitignore")

    "The needle the ignore-rule evidence is selected by. Deliberately the ROOT name and nothing more."
    const P14_MANIFEST_ROOT_TOKEN = "corpus"

    "The two tier labels the record reports by name, so a renamed tier shows up as a zero count."
    const P14_TIER_LABELS = ("simulated-secondary", "physical-primary")

    "The three split labels, likewise. `sealed_holdout` is Phase 16's, and is never read here."
    const P14_SPLIT_LABELS = ("dev", "eval", "sealed_holdout")

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

"""
    _p14_gitignore_evidence() -> NamedTuple

The ignore rules that keep the downloaded image bytes out of git, read at RUN TIME with their line
numbers, as `(available, line_numbers, lines, note)`.

**Computed, never transcribed.** 14-CONTEXT D-03a cites a line range for these rules that no longer
holds -- unrelated rules were inserted above them -- which is exactly the failure mode a quoted
line number has. Selecting the lines by the manifest's own root name means the record moves with
the file instead of drifting away from it, and a reader gets the rule text itself rather than a
pointer that may or may not still land on it.

This is the evidence for the word BY DESIGN in `corpus_status`. Without it, "the bytes are not on
disk" is a claim; with it, the reason is on the record and checkable.
"""
function _p14_gitignore_evidence()
    isfile(P14_GITIGNORE_PATH) || return (
        available = false, line_numbers = Int[], lines = String[],
        note = "no ignore file was found at $(P14_GITIGNORE_PATH); the by-design reading rests " *
               "on the manifest's own byte accounting alone")
    nums = Int[]
    lns  = String[]
    for (i, ln) in enumerate(eachline(P14_GITIGNORE_PATH))
        occursin(P14_MANIFEST_ROOT_TOKEN, ln) || continue
        push!(nums, i)
        push!(lns, rstrip(ln))
    end
    return (available = !isempty(nums), line_numbers = nums, lines = lns,
            note = "the ignore rules were located by name at run time, not transcribed from a " *
                   "cited line range; the range 14-CONTEXT D-03a quotes has since drifted")
end

"""
    _p14_manifest_metadata() -> NamedTuple

The provenance manifest's tier, split and byte accounting, COMPUTED FROM THE FILE.

# Metadata only, and what that excludes

Three columns are read -- `tier`, `split`, `bytes` -- from a committed text file. No image row is
opened, no directory under the manifest's root is listed, and the sole reachable accessor for the
sealed rows is neither called nor named. The narrower allowlisted property
`spike/test/test_p14_decoupling.jl` grants this runner is exactly that one.

# Why every count is derived rather than written down

The tier and split breakdowns are the whole substance of D-03a's first claim, and a hard-coded
"30 / 2" would keep reading correctly on the day the manifest changed and the claim stopped being
true. Deriving them means the record changes with the file. The two tier labels and the three split
labels ARE named -- so a renamed label shows up as a zero count rather than vanishing from a
collected set -- and the per-label counts are then asserted to exhaust the table.

`CSV.read` is used rather than a comma split because two rows carry quoted fields containing
commas; splitting on the delimiter shears them and silently mis-assigns every later column.
"""
function _p14_manifest_metadata()
    isfile(P14_MANIFEST_PATH) || error(
        "run_p14_real_images: the provenance manifest is absent at $(P14_MANIFEST_PATH). It is a " *
        "COMMITTED text file, so its absence means this is not the shipped tree -- and the SC1-f " *
        "availability record must be computed from it rather than asserted from prose.")
    df = CSV.read(P14_MANIFEST_PATH, DataFrame; comment = "#")
    for col in (:tier, :split, :bytes)
        hasproperty(df, col) || error(
            "run_p14_real_images: the manifest carries no `$col` column. The availability record " *
            "is computed from the file; a schema drift must fail BY NAME rather than produce a " *
            "record about columns that are not there.")
    end

    tiers  = String.(df.tier)
    splits = String.(df.split)
    bytes  = Int.(df.bytes)
    n_rows = length(tiers)

    tier_breakdown  = [(tier = t,  n = count(==(t), tiers))  for t in P14_TIER_LABELS]
    split_breakdown = [(split = s, n = count(==(s), splits)) for s in P14_SPLIT_LABELS]
    @assert sum(r.n for r in tier_breakdown) == n_rows "run_p14_real_images: the named tier labels $(P14_TIER_LABELS) account for $(sum(r.n for r in tier_breakdown)) of $n_rows manifest rows; a row in an unnamed tier would be invisible in the breakdown"
    @assert sum(r.n for r in split_breakdown) == n_rows "run_p14_real_images: the named split labels $(P14_SPLIT_LABELS) account for $(sum(r.n for r in split_breakdown)) of $n_rows manifest rows"

    is_physical = tiers .== "physical-primary"
    is_sealed   = splits .== "sealed_holdout"
    has_bytes   = bytes .> 0

    return (n_rows = n_rows,
            tier_breakdown = tier_breakdown,
            split_breakdown = split_breakdown,
            n_simulated_secondary = count(==("simulated-secondary"), tiers),
            n_physical_primary = count(is_physical),
            n_dev = count(==("dev"), splits),
            n_eval = count(==("eval"), splits),
            n_sealed_holdout = count(is_sealed),
            images_available = count(has_bytes),
            bytes_total = sum(bytes),
            n_physical_unsealed_with_bytes = count(is_physical .& (.!is_sealed) .& has_bytes),
            n_physical_sealed = count(is_physical .& is_sealed),
            path = P14_MANIFEST_PATH,
            sha256 = _p14_real_digest(P14_MANIFEST_PATH))
end

"""
    _p14_corpus_record() -> NamedTuple

The SC1-f availability record: what the provenance manifest holds, why none of it is read here, and
what that costs the claim.

# The four sentences this record exists to make un-droppable

**UNFETCHED, NOT GONE.** Every row records zero bytes and the ignore rules keep downloaded image
bytes out of git, so the data was never expected to be in a checkout. Reporting it as missing would
be false in the direction that flatters this phase -- it would make an unmade decision look like an
unavoidable circumstance. A fetch is a licence and bandwidth decision, explicitly a human one, and
out of Phase-14 scope.

**NO UNSEALED PHYSICAL TRUTH.** The overwhelming majority of rows are `simulated-secondary` -- the
CBS benchmark, i.e. another simulator -- and every `physical-primary` row is Phase 16's sealed
holdout. So even a completed fetch would not have produced the physical check D-03 wanted.

**THE SEALED HOLDOUT IS NOT READ.** Not in a script, not in a test, not behind a flag. There is no
second blind set, so reading one row would irreversibly burn the evaluation it exists for.

**THE MISSING BOUND IS ABSENT, NOT LOOSE.** `d03_bound_status = :absent_not_loose`, in those words,
because "loose" implies a weak bound exists and a reader may discount it accordingly. Nothing here
bounds the real-data behaviour at all; the six committed TIFFs illustrate it.
"""
function _p14_corpus_record()
    m  = _p14_manifest_metadata()
    gi = _p14_gitignore_evidence()

    return (manifest_metadata = m,
            gitignore_evidence = gi,
            corpus_images_available = m.images_available,
            corpus_status = :unfetched_by_design,
            corpus_status_note =
                "UNFETCHED BY DESIGN. Every one of the $(m.n_rows) manifest rows records zero " *
                "bytes (total $(m.bytes_total)), and the ignore rules recorded verbatim in " *
                "`gitignore_evidence` keep downloaded image bytes out of git while " *
                "`fetch.jl` beside the manifest is the fetcher. The data is therefore " *
                "unfetched: it was never expected to be present in a checkout, and it has not " *
                "been deleted or lost. Fetching it is a licence and bandwidth decision, " *
                "explicitly flagged as a human one, and it is OUT OF PHASE-14 SCOPE. The byte " *
                "accounting above is COMPUTED from the manifest; the sealed image directory is " *
                "never listed, opened or named by this runner.",
            tier_breakdown = m.tier_breakdown,
            split_breakdown = m.split_breakdown,
            physical_truth_available = m.n_physical_unsealed_with_bytes > 0,
            physical_truth_note =
                "The only $(m.n_physical_primary) `physical-primary` rows are Phase 16's sealed " *
                "holdout ($(m.n_physical_sealed) of them), and they are NOT read here -- not in " *
                "a script, not in a test, not behind a flag. There is no second blind set, so " *
                "consuming one would irreversibly burn the blind evaluation on the very " *
                "hypothesis it exists to evaluate. The remaining " *
                "$(m.n_simulated_secondary) rows are `simulated-secondary` (the CBS benchmark, " *
                "i.e. another simulator) and are not physical ground truth either, so a " *
                "completed fetch would still not have produced the physical check D-03 wanted.",
            substitution_record =
                (from = "the provenance manifest's physical-primary anchors, both of which are " *
                        "`split = sealed_holdout` and reserved for the Phase-16 blind evaluation",
                 to = "the six committed microscopy TIFFs under " * P14_REAL_SUBSTRATE_ROOT,
                 precedent = "spike/p13/real_images.jl:71,145",
                 amendment = ".planning/phases/13-three-hypothesis-amortized-bayes-factor/" *
                             "13-D15-AMENDMENT.md",
                 rationale = "precisely so that Phase 16's blind corpus stays blind",
                 decision = "D-03a, which supersedes D-03's corpus half and not its simulator half",
                 sealed_holdout_read = false,
                 metadata_only = true),
            withdrawn_figure_note =
                "WITHDRAWN, NOT MERELY SUPERSEDED. D-03 recorded a miscoverage floor derived " *
                "from the manifest -- the ratio of one to thirty-one, roughly three per cent -- " *
                "and it must not appear in any Phase-14 artifact or report. It was computed from " *
                "the wrong n (the manifest is partitioned dev = $(m.n_dev) / eval = $(m.n_eval) " *
                "/ sealed_holdout = $(m.n_sealed_holdout), not `30 open`) AND on rows that are " *
                "not physical truth ($(m.n_simulated_secondary) of $(m.n_rows) are " *
                "`simulated-secondary`). The numeral is banned from Phase-14 executable code by " *
                "`spike/test/test_p14_decoupling.jl`, so it is written here in prose rather " *
                "than as a literal -- which is also why a reader will not find it as a number " *
                "anywhere in this artifact.",
            d03_bound_status = :absent_not_loose,
            d03_bound_note =
                "The bounding evidence D-03 intended is UNAVAILABLE, so the bound is ABSENT, NOT " *
                "MERELY LOOSE. Those are the words, and the difference matters: `loose` invites " *
                "a reader to discount a weak bound, whereas nothing in this phase bounds the " *
                "layer's behaviour on real microscopy at all. The conformal guarantee remains " *
                "SIMULATOR-DERIVED and inherits the simulator's misspecification in full; the " *
                "six committed TIFFs ILLUSTRATE the behaviour and bound nothing.",
            corpus_metadata_only = true,
            corpus_images_opened = 0,
            sealed_holdout_read = false)
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
        # T-14-04, asserted on BOTH sides of the metadata read: the provenance tree is byte-
        # unchanged and untracked-clean. A metadata read cannot dirty it, which is precisely why
        # an assertion here is cheap and why its failure would mean something other than a read
        # happened.
        prov_clean     = success(Cmd(`git diff --quiet HEAD -- corpus`; dir = P14_REPO_ROOT))
        prov_untracked = readchomp(Cmd(`git status --porcelain -- corpus`; dir = P14_REPO_ROOT))
        @assert prov_clean "T-14-04: the provenance tree is not byte-unchanged against HEAD before the run; this runner reads its manifest as METADATA ONLY and may not run against a modified one"
        @assert isempty(prov_untracked) "T-14-04: an untracked file appeared in the provenance tree ($(prov_untracked)); a fetch or a stray write would make the availability record describe something other than the committed manifest"
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
    verbose && println("[1/6] asserting and digesting the six committed TIFFs …")
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
    verbose && println("[2/6] loading both conditions and deciding at the LOADED q-hat …")
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

    # --- 8. THE CORPUS AVAILABILITY RECORD, METADATA ONLY -------------------------------------
    # Computed from the committed manifest text, never copied from prose and never hard-coded, so
    # the record changes if the manifest does. No image row is opened; the sealed holdout is not
    # read -- not in a script, not in a test, not behind a flag.
    verbose && println("[3/6] computing the corpus availability record (METADATA ONLY) …")
    corpus = _p14_corpus_record()

    # --- 9. THE READ-ONLY PROOF, RE-TAKEN -----------------------------------------------------
    digest_after = verify_real_readonly_digest()

    # --- 10. PERSIST THE REPORT, BEFORE ANY TABLE IS PRINTED ----------------------------------
    verbose && println("[4/6] persisting the SC1-f report (BEFORE anything is printed) …")
    consts_path = joinpath(@__DIR__, "consts.jl")
    _p14_real_save_report(report_path,
        ("images", "decision_units", "predicted_behaviour", "prediction_held", "observed_decisions",
         "is_coverage_claim", "illustration_note", "readonly_digest_before",
         "readonly_digest_after", "substrate_root", "provenance",
         # --- the availability record: every key of it REQUIRED, so a report that dropped the
         #     honest account of what was NOT read is never written at all.
         "corpus_images_available", "corpus_status", "corpus_status_note", "tier_breakdown",
         "split_breakdown", "physical_truth_available", "substitution_record",
         "withdrawn_figure_note", "d03_bound_status");
        # --- THE AVAILABILITY RECORD, splatted in so it CANNOT be dropped ---
        corpus...,
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

    # --- 11. THE SIX-ROW TABLE, PRINTED AFTER THE ARTIFACT IS ON DISK -------------------------
    if verbose
        println("[5/6] the six committed TIFFs")
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
        println("[6/6] the availability record -- METADATA ONLY, no image row was opened")
        println("  manifest             = $(corpus.manifest_metadata.path)")
        println("  rows                 = $(corpus.manifest_metadata.n_rows)   " *
                "sha256 = $(corpus.manifest_metadata.sha256)")
        println("  TIER BREAKDOWN (computed from the file, never hard-coded)")
        for r in corpus.tier_breakdown
            println("    $(rpad(r.tier, 22)) $(r.n)")
        end
        println("  SPLIT BREAKDOWN")
        for r in corpus.split_breakdown
            println("    $(rpad(r.split, 22)) $(r.n)")
        end
        println("  images available     = $(corpus.corpus_images_available)   " *
                "total bytes = $(corpus.manifest_metadata.bytes_total)")
        println("  status               = $(corpus.corpus_status)")
        println("  physical truth available = $(corpus.physical_truth_available)")
        println("  sealed holdout read  = $(corpus.sealed_holdout_read)   " *
                "images opened = $(corpus.corpus_images_opened)")
        println("  D-03 bound status    = $(corpus.d03_bound_status)")
        println("  IGNORE-RULE EVIDENCE (located at run time, not transcribed)")
        if corpus.gitignore_evidence.available
            for (n, l) in zip(corpus.gitignore_evidence.line_numbers,
                              corpus.gitignore_evidence.lines)
                println("    .gitignore:$(rpad(n, 5)) $l")
            end
        else
            println("    $(corpus.gitignore_evidence.note)")
        end
        println("  UNFETCHED, NOT MISSING. The bytes were never expected to be in a checkout; a")
        println("  fetch is a human licence-and-bandwidth decision and is out of Phase-14 scope.")
        println("  THE MISSING D-03 BOUND IS ABSENT, NOT MERELY LOOSE.")
        println("  THE WITHDRAWN D-03 FIGURE IS RECORDED AS WITHDRAWN AND APPEARS AS NO NUMERAL.")
        println()
        println("THE PREDICTION AND ITS OUTCOME, SIDE BY SIDE")
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

    # --- 12. THE STRUCTURAL ASSERTIONS, LAST --------------------------------------------------
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
            corpus = corpus, observed_decisions = observed, prediction_held = prediction_held,
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
