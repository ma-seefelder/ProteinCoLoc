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

# spike/p13/run_p13_realimage.jl --- D-15 (AMENDED) qualitative real-image check. REPORTED, NOT GATED.
#
# =============================================================================================
# (a) CHARACTERIZATION, NOT A GATE.
# =============================================================================================
# THIS RUN CANNOT FAIL THE PHASE. THE PHASE-13 GATE IS SIMULATOR GROUND TRUTH -- the D-12
# discrimination floors and the D-13 calibration band, evaluated on drawn theta by
# spike/p13/run_three_way_gate.jl and reported there. NOTHING THIS FILE PRINTS FEEDS THAT GATE.
#
# THE SIX COMMITTED test/test_images/ TIFFs CARRY NO COLOCALIZATION GROUND-TRUTH LABEL. There is
# no recorded rho, no recorded Delta-rho, no segregation degree and no provenance metadata
# asserting one. SO THIS ARM CAN SHOW THAT THE THREE-WAY VERDICTS BEHAVE SENSIBLY ON REAL
# MICROSCOPY; IT CANNOT SHOW THAT THEY ARE CORRECT, because correctness requires a label the data
# does not have. NO PASS/FAIL THRESHOLD IS DEFINED FOR ANY QUANTITY THIS RUNNER PRINTS
# (P13_REAL_IS_GATED = false, P13_REAL_QUALITATIVE_ONLY = true). A DISAPPOINTING, FLAT OR
# CONTRADICTORY READING IS A REPORTABLE FINDING, NOT A FAILURE, AND IS NEVER A REASON TO RETRAIN,
# RESEED, RETUNE OR WIDEN THE SUBSTRATE. Labelled real segregation validation is DEFERRED TO THE
# PHASE-16 BLIND EVALUATION.
#
# =============================================================================================
# (b) THE NAMING CORRECTION (P13_REAL_NAMING_CORRECTION -- carried here, not only in the report).
# =============================================================================================
# THE FOLDER NAMES positive/ AND negative/ ARE THE ORIGINAL PACKAGE'S *BIOLOGICAL* TEST
# CONDITIONS, NOT COLOCALIZATION LABELS. The "negative" pair measures mean patch correlation
# +0.2481 (rho_true ~ +0.215) -- A POSITIVELY CORRELATED PAIR, NOT AN ANTI-CORRELATED ONE.
# NOTHING HERE TREATS negative/ AS AN EXCLUSION EXAMPLE, and a reader who assumes otherwise will
# misread every real-image figure in the phase.
#
# THE CHANNEL-PAIR CAVEAT, RESTATED BECAUSE IT SHARPENS (b) RATHER THAN SOFTENING IT. The Tier-1
# pre-registration freezes P13_REAL_CHANNEL_PAIR = (1, 2) and P13_REAL_REDUNDANCY_PAIR = (1, 3),
# and test/runtests.jl:105 records the fixture channels as ["blue", "green", "red"] -- so channel
# 1 is the DAPI/Hoechst NUCLEAR COUNTERSTAIN, not a target protein. BOTH pre-registered pairs
# therefore involve the counterstain. spike/simulator/ghat.jl records (second hand-patch entry,
# 2026-07-25) that the c1/c2 figures are SUPERSEDED FOR COLOCALIZATION PURPOSES by the c2/c3
# green/red pair. THIS RUNNER DOES NOT ADD A THIRD, UN-PRE-REGISTERED PAIR: choosing a new
# channel pair after the fixtures had been measured is exactly the move the D-04 anti-snooping
# contract forbids. The frozen pairs are read, and the limit is NAMED rather than repaired.
#
# =============================================================================================
# (c) THE TARGET-SUBSTITUTION RECORD (P13_REAL_SUBSTITUTION_RECORD, 11-10-PLAN.md house shape).
# =============================================================================================
# D-15 ORIGINALLY NAMED THE PHYSICAL MITOCHONDRIA ANCHOR OF THE FROZEN PROVENANCE MANIFEST as the
# real-data check. Verification showed both physical-primary rows are sha256 = "PENDING-FETCH",
# bytes = 0 and split = sealed_holdout -- UNFETCHED, and reserved behind the anti-snooping control
# for the PHASE-16 BLIND EVALUATION. Consuming one here would irreversibly burn Phase 16 on the
# very hypothesis Phase 16 exists to evaluate. THE SANCTIONED SUBSTITUTE is test/test_images/ --
# six real committed microscopy TIFFs already exercised by test/runtests.jl and already the frozen
# anchor of this project's own simulator calibration.
#
# THE ANTI-SNOOPING SEAL IS LEFT INTACT. The sealed-holdout accessor -- `open_sealed_holdout(df;
# reason)` -- IS NOT CALLED: not in a script, not in a test, not behind a flag, and no seal-break
# escape hatch exists or may be re-introduced. Naming it here, inside a whole-line comment, is
# provenance prose, which is the only form the ruling permits; the machine-checked assertion below
# proves the name appears nowhere in this file's executable text.
#
# =============================================================================================
# (d) THE OOD FINDING, STATED UP FRONT.
# =============================================================================================
# THE SHIPPED DETECTOR FLAGS BOTH FIXTURES: Mahalanobis summary density
# P13_REAL_OOD_SHIPPED_DENSITY = 433.69 against the frozen in-distribution threshold
# P13_REAL_OOD_SHIPPED_THRESHOLD = 179.14 -- 2.4x over, in BOTH read directions. THAT IS THE
# NET'S OWN MISSPECIFICATION CHANNEL DOING ITS JOB, and it is arguably the most honest single
# number this arm produces. This run measures the PHASE-13 net's verdict on the same inputs and
# prints it BESIDE the shipped reference. NEITHER VERDICT IS SUPPRESSED, SOFTENED OR "FIXED"
# (Pitfall 13). AGREEMENT CORROBORATES THE NAMED LIMIT; DIVERGENCE IS ITSELF A RESULT about what
# the registration-aware basis buys on real data. EVERY REAL-IMAGE QUANTITY IS PRINTED AND
# PERSISTED ON THE SAME ROW AS ITS OOD DENSITY, ITS OOD THRESHOLD, ITS BOOLEAN VERDICT AND ITS
# LAMBDA.
#
# THE SHIPPED BUNDLE IS READ AS A FROZEN COMPARISON REFERENCE ONLY AND IS NEVER AN EVIDENCE
# BASIS. D-02 forbids falling back to the shipped grid-8 basis for the EVIDENCE read; it does not
# forbid quoting that bundle's own OOD numbers for comparison, which is what
# P13_REAL_OOD_COMPARISON = true authorises.
#
# =============================================================================================
# (e) READ-ONLY DISCIPLINE.
# =============================================================================================
# test/test_images/ IS READ-ONLY INPUT. The `git ls-files -s` digest is recorded BEFORE and AFTER
# the run and asserted identical, so "we did not write there" is a checked property rather than a
# stated intention. Nothing under test/ is written, moved or modified by any line below.
#
# =============================================================================================
# (f) DECOUPLING AND HOW TO RUN IT.
# =============================================================================================
# DECOUPLING (hard constraint, CLAUDE.md, D-01): spike-local. src/ is reached ONLY read-only and
# transitively through the summary-contract include chain, and step 0 PROVES at run time that
# src/ and both spike manifests are byte-unchanged. CPU-only; CUDA is never imported; no package
# is installed and the figure reuses the CairoMakie surface already resolved in
# spike/validation/figures.jl. THIS ARM CONSUMES NO RNG STREAM AT ALL: the substrate is six
# committed files, the alpha transform is a deterministic function of its input, its mask and its
# floor, and a forward pass through a trained net is deterministic -- so there is no seed to burn
# and no reserved counter this file could pre-observe.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY. Including it does nothing; run it deliberately:
#
#     julia --project=spike spike/p13/run_p13_realimage.jl
#
# ===== END OF EXPLANATORY HEADER -- SEAL SELF-CHECK SENTINEL, DO NOT REWORD =====

using JLD2
using Statistics
using LinearAlgebra
using StatsBase
using Dates
using SHA

# ORDER MATTERS. datagen.jl first -- it carries the whole Phase-11 precondition preamble
# (preconditions.jl -> p11_architecture.jl, consts.jl, labels.jl, net.jl, the simulator halves and
# the src/ summary contract) AND the labelled-pool reader this arm's in-distribution reference is
# read from. Then the alpha transform, then the real-image ingestion, then the figure surface.
# All guarded, the house guarded-include idiom.
isdefined(@__MODULE__, :load_pool)               || include(joinpath(@__DIR__, "datagen.jl"))
isdefined(@__MODULE__, :verify_alpha_invariants) || include(joinpath(@__DIR__, "alpha_series.jl"))
isdefined(@__MODULE__, :real_alpha_ladder)       || include(joinpath(@__DIR__, "real_images.jl"))
isdefined(@__MODULE__, :plot_roc)                ||
    include(joinpath(@__DIR__, "..", "validation", "figures.jl"))

# --- (0) The MACHINE-CHECKED seal assertion ---------------------------------------------------
# The 11-10-PLAN.md:191-195 mechanism, applied one phase later and asserted on this file's OWN
# text: strip the explanatory header above (everything up to and including the sentinel line) and
# assert the remainder names neither the sealed provenance directory nor the sealed-holdout
# accessor. That turns "we intended not to reach the seal" into "this file FAILS TO LOAD if we
# did".
#
# BOTH FORBIDDEN NAMES ARE ASSEMBLED FROM CHARACTER FRAGMENTS AT THE ASSERTION SITE, ON PURPOSE,
# so the assertion does not itself embed the literals it forbids -- otherwise the check would trip
# on its own text and would have to be weakened into uselessness. The same construction guards
# spike/p13/real_images.jl and spike/test/test_p13_real.jl.
let src = read(@__FILE__, String)
    marker = "END OF EXPLANATORY HEADER"
    hit    = findfirst(marker, src)
    body   = hit === nothing ? src : src[(last(hit) + 1):end]
    forbidden_dir = join(('c', 'o', 'r', 'p', 'u', 's'))
    forbidden_fn  = "open" * "_sealed" * "_holdout"
    occursin(forbidden_dir, body) && error(
        "run_p13_realimage.jl: the sealed provenance directory is named in executable code. The " *
        "amended D-15 substitutes test/test_images/ precisely so the Phase 16 blind evaluation " *
        "is not burned here; remove the reference.")
    occursin(forbidden_fn, body) && error(
        "run_p13_realimage.jl: the sealed-holdout accessor is named in executable code. It is " *
        "never called -- not in a script, not in a test, not behind a flag -- and no seal-break " *
        "escape hatch may be re-introduced.")
end

if !isdefined(@__MODULE__, :P13_REALIMAGE_REPORT_PATH)
    "The reported real-image artifact. Written BEFORE the read-only digest comparison can throw."
    const P13_REALIMAGE_REPORT_PATH = joinpath(@__DIR__, "realimage_report.jld2")

    "The trained three-way evidence net (plan 13-11). There is no substitute and no fallback."
    const P13_REALIMAGE_NET_PATH = joinpath(@__DIR__, "three_way_net.jld2")

    "The two-panel real-image figure. A SCRIPT ARTIFACT, never a gate assertion."
    const P13_REALIMAGE_FIG_PATH = joinpath(@__DIR__, "..", "figures", "p13_realimage.png")

    """
        P13_REALIMAGE_SHIPPED_CANDIDATES

    The frozen shipped-bundle directories, RELATIVE to the repository artifacts root, searched in
    this declaration order for the OOD COMPARISON REFERENCE (P13_REAL_OOD_COMPARISON = true).

    WHICH ONE IS THE REFERENCE IS *DERIVED, NOT CHOSEN*. The selected bundle is the first whose
    OWN recorded `report.ood.id_threshold` reproduces the frozen constant
    `P13_REAL_OOD_SHIPPED_THRESHOLD`. That makes the pre-registration itself the selector: the
    reference is the bundle the frozen number came from, and if no bundle on disk reproduces it,
    the run says so and quotes the frozen constants alone rather than silently comparing against
    some other net.

    Every candidate is opened READ-ONLY, for its standardizer, its density null and its recorded
    operating point. Nothing under the artifacts tree is written, and the Phase-7 ship gate is not
    re-run, not re-scored and not reopened.
    """
    const P13_REALIMAGE_SHIPPED_CANDIDATES = (joinpath("amended_v2", "grid_8"), "grid_8")

    """
        P13_REALIMAGE_OOD_ID_QUANTILE

    ATTRIBUTED COPY of the shipped `OOD_ID_QUANTILE` (`src/amortized/ood.jl:61`), the
    pre-registered in-distribution operating quantile (~5% ID false-positive rate) that the
    shipped detector's own threshold was set at.

    IT IS COPIED, NOT CHOSEN, AND IT SCORES NOTHING IN THIS ARM. It exists so the Phase-13 OOD
    channel is constructed by the IDENTICAL recipe as the shipped one -- same continuous-rows-only
    Mahalanobis fit, same ridge, same ID quantile -- because a side-by-side comparison of two
    detectors built by different recipes would be a statement about the recipes rather than about
    the nets. No Phase-13 verdict turns on the resulting number: P13_REAL_IS_GATED = false, and an
    OOD flag here is a REPORTED DIAGNOSTIC that travels beside every quantity, never a bar a
    result is scored against.
    """
    const P13_REALIMAGE_OOD_ID_QUANTILE = 0.95

    "Mahalanobis covariance ridge, copied verbatim from the shipped `fit_ood_nulls` (Pitfall 2)."
    const P13_REALIMAGE_OOD_RIDGE = 1e-6
end

# =============================================================================================
# Small local helpers
# =============================================================================================

"""
    _p13_real_save_report(path; kwargs...) -> String

Atomically persist the reported artifact: write `path * ".tmp"`, REOPEN to integrity-check the
required keys, then `mv(...; force = true)` -- the `_save_bf_report` idiom
(`spike/validation/run_bf.jl:59-69`). A crash mid-write leaves a discardable `.tmp`, never a torn
artifact a later reader would happily believe.
"""
function _p13_real_save_report(path; kwargs...)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        for k in ("is_gated", "lambda_sweep", "ood_phase13", "ood_shipped",
                  "alpha_ladder_real", "digest_before", "digest_after")
            @assert haskey(f, k) "_p13_real_save_report: integrity check failed, $tmp missing $k"
        end
    end
    mv(tmp, path; force = true)
    return path
end

"""
    p13_real_blob_sha(path) -> String

The GIT BLOB sha1 of `path`, i.e. exactly what `git hash-object <path>` prints, so the artifact
points INTO the audit trail: a reader can `git cat-file -p <sha>` the frozen pre-registration text
out of history, or `git log --find-object=<sha>` the commit that introduced it.
"""
function p13_real_blob_sha(path)
    bytes = read(path)
    ctx   = SHA.SHA1_CTX()
    SHA.update!(ctx, Vector{UInt8}("blob $(length(bytes))\0"))
    SHA.update!(ctx, bytes)
    return bytes2hex(SHA.digest!(ctx))
end

"""
    p13_real_cont_rows(nrows) -> UnitRange{Int}

The CONTINUOUS rows of a `:min` summary, `1:(nrows ÷ 2)`.

ATTRIBUTED COPY of `_summary_row_partition(:min, nrows)` (`src/amortized/summary.jl`), reproduced
here rather than reached through an `include` because that file's method signatures mention
`MultiChannelImage` and it is not loadable standalone. `src/` stays byte-unchanged.

THE MASK ROWS ARE NEVER READ. They are near-constant 0/1 indicators, so including them makes the
covariance singular -- the Pitfall-2 the shipped fit documents. The split is derived from the
realized row count, never written as a literal.
"""
function p13_real_cont_rows(nrows::Integer)
    iseven(nrows) || error(
        "p13_real_cont_rows: a :min summary has an even row count (2*G^2), got $nrows")
    return 1:(nrows ÷ 2)
end

"""
    p13_real_fit_density(Z) -> NamedTuple

Fit the summary-density Mahalanobis null on a `d x N` matrix of STANDARDIZED summaries, returning
`(cont, μS, C)` with `C = cholesky(Symmetric(ΣS + ridge*I))`.

ATTRIBUTED COPY of the shipped `fit_ood_nulls` (`src/amortized/ood.jl:77`, twin at
`spike/validation/ood.jl:112`), reproduced under a `p13_` name for the same reason `net.jl`
reproduces `pair_encode`: the shipped file is not loadable standalone and the spike twin drags the
whole Phase-5 validation harness (a frozen Phase-4 model load) into a file that needs one
covariance. Copying rather than re-deriving is what makes the two detectors comparable.
"""
function p13_real_fit_density(Z::AbstractMatrix; ridge::Real = P13_REALIMAGE_OOD_RIDGE)
    cont = p13_real_cont_rows(size(Z, 1))
    Zc   = Float64.(Z[cont, :])
    μS   = vec(mean(Zc; dims = 2))
    ΣS   = cov(Zc; dims = 2)
    C    = cholesky(Symmetric(ΣS + ridge * I))
    return (cont = collect(cont), μS = μS, C = C)
end

"""
    p13_real_maha(nulls, z) -> Float64

Squared Mahalanobis distance of a standardized summary `z` to a density null, on the CONTINUOUS
rows only. ATTRIBUTED COPY of the shipped `maha_score` (`src/amortized/ood.jl`); an
in-distribution summary scores LOW, a summary shifted off the training manifold scores HIGH.
"""
function p13_real_maha(nulls, z::AbstractVector)
    r = Float64.(z[nulls.cont]) .- nulls.μS
    return sum(abs2, nulls.C.L \ r)
end

"""
    p13_real_pair_density(nulls, Zs, Zc) -> Float64

The PAIR density score `max(maha(Zs), maha(Zc))`, exactly as the shipped `ood_verdict` fuses a
(sample, control) analysis on the density channel. Reported alongside the two per-acquisition
scores so the asymmetry between the two fixtures stays visible rather than being hidden by the max.
"""
p13_real_pair_density(nulls, Zs, Zc) = max(p13_real_maha(nulls, Zs), p13_real_maha(nulls, Zc))

"""
    p13_real_shipped_reference() -> NamedTuple

Resolve and read the frozen SHIPPED grid-8 comparison reference, READ-ONLY.

Returns `(available, bundle, zt, density, threshold, reproduced_frozen, note)`. The bundle is the
first `P13_REALIMAGE_SHIPPED_CANDIDATES` entry whose recorded `report.ood.id_threshold` reproduces
the frozen `P13_REAL_OOD_SHIPPED_THRESHOLD` -- DERIVED FROM THE PRE-REGISTRATION, never picked by
hand. When no bundle on disk reproduces it the first readable bundle is reported with
`reproduced_frozen = false`, and when none is readable the arm reports the frozen constants alone.

The operating point lives in `gate_report_8.jld2` rather than beside the null fit, exactly as
`src/amortized/local_map.jl:93-152` reads it. Only the `zt` key of `npe_8.jld2` is opened, so no
estimator object and no theta transform is reconstructed.
"""
function p13_real_shipped_reference(; root = joinpath(P13_REPO_ROOT, "artifacts"))
    first_readable = nothing
    for rel in P13_REALIMAGE_SHIPPED_CANDIDATES
        dir  = joinpath(root, rel)
        npe  = joinpath(dir, "npe_8.jld2")
        ood  = joinpath(dir, "ood_nulls_8.jld2")
        rep  = joinpath(dir, "gate_report_8.jld2")
        (isfile(npe) && isfile(ood) && isfile(rep)) || continue
        thr = nothing
        try
            r   = JLD2.jldopen(f -> f["report"], rep, "r")
            thr = Float64(r.ood.id_threshold)
        catch err
            @warn "p13_real_shipped_reference: could not read the recorded operating point" dir err
            continue
        end
        zt    = JLD2.jldopen(f -> f["zt"], npe, "r")
        nulls = JLD2.jldopen(f -> f["ood_nulls"], ood, "r").density
        hit   = (available = true, bundle = rel, zt = zt, density = nulls, threshold = thr,
                 reproduced_frozen = isapprox(thr, P13_REAL_OOD_SHIPPED_THRESHOLD; atol = 0.005),
                 note = "READ-ONLY frozen COMPARISON REFERENCE, never an evidence basis (D-02).")
        hit.reproduced_frozen && return hit
        first_readable === nothing && (first_readable = hit)
    end
    first_readable === nothing && return (
        available = false, bundle = nothing, zt = nothing, density = nothing,
        threshold = P13_REAL_OOD_SHIPPED_THRESHOLD, reproduced_frozen = false,
        note = "No shipped bundle on disk. The frozen constants " *
               "P13_REAL_OOD_SHIPPED_DENSITY / P13_REAL_OOD_SHIPPED_THRESHOLD are quoted alone.")
    return first_readable
end

"""
    p13_real_id_reference(pool_dir) -> NamedTuple

Build the PHASE-13 in-distribution reference: fit the density null on the standardized summaries
of the net's OWN labelled training pool, and take the operating point at the copied shipped ID
quantile.

`pool_dir` is READ OFF THE TRAINED NET'S PERSISTED METADATA (`meta.pool_dir`), never guessed and
never re-derived through the cache layer -- so this function cannot create, invalidate or
regenerate a pool. Every shard is opened READ-ONLY through `load_pool`; nothing is written.

BOTH HALVES OF EVERY PAIR ARE USED. `Zs` and `Zc` are both simulator acquisitions standardized
through the SAME frozen Phase-11 `zt`, so the union is precisely "the acquisitions this net was
trained on", which is what an in-distribution reference should be referenced to.

DISCLOSED, NOT BURIED: that pool is CLASS-FREQUENCY STRATIFIED (D-07, `P13_STRATIFICATION =
:class_frequency`), so its rho marginal is not the prior's. The shipped null was fit on
prior-distributed train summaries. The two references are therefore not identical constructions of
the same distribution -- they are each net's own training distribution, which is the comparison the
arm is actually making. Recorded in the artifact so nobody has to infer it.

Returns `(available, nulls, threshold, n_id, scores_quantiles, pool_dir, n_shards, note)`. An
absent, unreadable or schema-shifted pool yields `available = false` with the reason recorded: the
arm is not gated, so a missing in-distribution reference is reported rather than fatal.
"""
function p13_real_id_reference(pool_dir)
    bad(note) = (available = false, nulls = nothing, threshold = nothing, n_id = 0,
                 scores_quantiles = nothing, pool_dir = pool_dir, n_shards = 0, note = note)
    (pool_dir isa AbstractString && isdir(pool_dir)) ||
        return bad("the trained net records no readable pool directory, so the Phase-13 " *
                   "in-distribution reference could not be built")
    shards = sort(filter(p -> startswith(basename(p), "shard_") && endswith(p, ".jld2"),
                         readdir(pool_dir; join = true)))
    isempty(shards) && return bad("the recorded pool directory holds no shard file")
    cols = Any[]
    n_sh = 0
    try
        for s in shards
            p = load_pool(s)
            push!(cols, p.Zs)
            push!(cols, p.Zc)
            n_sh += 1
        end
    catch err
        return bad("a pool shard could not be read: $(sprint(showerror, err))")
    end
    Z      = reduce(hcat, cols)
    nulls  = p13_real_fit_density(Z)
    scores = [p13_real_maha(nulls, Float64.(@view Z[:, j])) for j in axes(Z, 2)]
    thr    = quantile(scores, P13_REALIMAGE_OOD_ID_QUANTILE)
    qs     = (q50 = quantile(scores, 0.50), q90 = quantile(scores, 0.90),
              q95 = quantile(scores, 0.95), q99 = quantile(scores, 0.99),
              qmax = maximum(scores))
    return (available = true, nulls = nulls, threshold = thr, n_id = size(Z, 2),
            scores_quantiles = qs, pool_dir = pool_dir, n_shards = n_sh,
            note = "Density null fit on the Phase-13 labelled training pool (both acquisitions " *
                   "of every pair), standardized through the FROZEN Phase-11 zt. Operating " *
                   "point at the COPIED shipped ID quantile. Scores nothing.")
end

"""
    p13_real_figure(sweep, ladders, headline; path) -> String

The two-panel real-image figure: the lambda sweep of both read directions on one panel, the real
alpha ladders of both conditions on the other, with OOD-flagged points VISUALLY MARKED.

A SCRIPT ARTIFACT AND NEVER A GATE ASSERTION -- the Phase-4 rule `spike/validation/figures.jl`
states at its head, honoured here by rendering the figure before anything else can throw. Composed
in this runner rather than through a `figures.jl` helper for the two reasons plans 13-12 and 13-13
already recorded for their own figures: that file has no sweep or ladder helper, and its
`_val_fig_path` hard-routes every output into `spike/validation/figures/`, which is not where this
deliverable lives. The CairoMakie surface, the headless backend activation and the caption
discipline are all inherited from it. No plotting dependency is added.
"""
function p13_real_figure(sweep, ladders, headline; path = P13_REALIMAGE_FIG_PATH)
    mkpath(dirname(path))
    fig = Figure(size = (1180, 500))

    ax1 = CairoMakie.Axis(fig[1, 1];
        title    = "real-image three-way evidence across the Phase-11 lambda ladder",
        subtitle = "QUALITATIVE, REPORTED, NOT A GATE -- open circles mark OOD-FLAGGED reads",
        xlabel   = "lambda (registration uncertainty, px -- a real image has NO KNOWN lambda)",
        ylabel   = "log Bayes factor (nats)")
    hlines!(ax1, [0.0]; color = :gray, linestyle = :dash)
    dirs = unique([r.direction for r in sweep])
    cols = (:firebrick, :steelblue)
    for (k, dname) in enumerate(dirs)
        rows = filter(r -> r.direction == dname, sweep)
        xs   = [r.lambda for r in rows]
        lines!(ax1, xs, [r.logbf_exclusion for r in rows];
               color = cols[k], label = "log BF(E:R)  $dname")
        lines!(ax1, xs, [r.logbf_coloc for r in rows];
               color = cols[k], linestyle = :dot, label = "log BF(C:R)  $dname")
        # AN OOD-FLAGGED READ IS DRAWN HOLLOW, AN IN-DISTRIBUTION READ FILLED, so no point on
        # this figure can be read without its misspecification verdict (T-13-56). The two states
        # are drawn as two INDEX-SPLIT scatter calls rather than as one call with a per-point
        # alpha vector: the pinned Makie resolves a scalar alpha only, and a per-point vector
        # fails the whole figure at draw time.
        flagged = [r.ood_is_ood for r in rows]
        for (yvals, mk, ms) in ((Float64[r.logbf_exclusion for r in rows], :circle, 11),
                                (Float64[r.logbf_coloc for r in rows], :utriangle, 10))
            hol = findall(flagged)
            fil = findall(.!flagged)
            isempty(hol) || scatter!(ax1, xs[hol], yvals[hol]; color = (:white, 0.0),
                                     marker = mk, markersize = ms,
                                     strokewidth = 2, strokecolor = cols[k])
            isempty(fil) || scatter!(ax1, xs[fil], yvals[fil]; color = cols[k],
                                     marker = mk, markersize = ms)
        end
    end
    vlines!(ax1, [headline]; color = :seagreen, linestyle = :dashdot, linewidth = 2)
    axislegend(ax1; position = :lb, labelsize = 9)

    ax2 = CairoMakie.Axis(fig[1, 2];
        title    = "real-substrate alpha ladder through the trained net",
        subtitle = "alpha = 0 is MODERATELY COLOCALIZED here; NEVER averaged with 13-13",
        xlabel   = "alpha (construction parameter, not a physical quantity)",
        ylabel   = "log Bayes factor (nats)")
    hlines!(ax2, [0.0]; color = :gray, linestyle = :dash)
    for (k, L) in enumerate(ladders)
        lines!(ax2, L.alpha, L.logbf_exclusion; color = cols[k],
               label = "log BF(E:R)  $(L.condition)")
        scatter!(ax2, L.alpha, L.logbf_exclusion; color = cols[k])
        lines!(ax2, L.alpha, L.logbf_coloc; color = cols[k], linestyle = :dot,
               label = "log BF(C:R)  $(L.condition)")
        L.alpha_star === nothing ||
            vlines!(ax2, [L.alpha_star]; color = cols[k], linestyle = :dashdot)
    end
    axislegend(ax2; position = :lt, labelsize = 9)

    save(path, fig)
    return path
end

# =============================================================================================
# The reported run
# =============================================================================================

"""
    main(; ...) -> NamedTuple

The amended D-15 real-image arm, end to end and once. There is no gate and no smoke/reported mode
split: every number this function produces is reported, none of it is scored against anything, and
nothing it prints can fail the phase.
"""
function main(; net_path = P13_REALIMAGE_NET_PATH,
              report_path = P13_REALIMAGE_REPORT_PATH,
              fig_path = P13_REALIMAGE_FIG_PATH,
              verbose::Bool = true)

    t_start = time()

    # --- 0. THE DECOUPLING PROOF, AT RUN TIME (D-01) ---------------------------------------
    # Asserted while the run happens, not only in review: a reported run that started on a clean
    # tree and finished on a dirty one would be a repudiation hole (T-13-39).
    src_clean = success(Cmd(`git diff --quiet HEAD -- src`; dir = P13_REPO_ROOT))
    env_clean = success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                            dir = P13_REPO_ROOT))
    @assert src_clean "D-01 decoupling breach: `git diff --quiet HEAD -- src` FAILED -- src/ is not byte-unchanged, so this run may not be reported"
    @assert env_clean "D-01 decoupling breach: spike/Project.toml or spike/Manifest.toml is modified -- the spike environment moved under the run"
    verbose && println("STEP 0 (D-01): src/ byte-unchanged = $src_clean   " *
                       "spike Project/Manifest byte-unchanged = $env_clean")

    # --- 1. THE READ-ONLY DIGEST, BEFORE ANYTHING IS OPENED --------------------------------
    digest_before = verify_real_readonly_digest()
    verbose && println("STEP 1 (T-13-50): read-only digest BEFORE = $digest_before")

    # --- 2. THE HARD PRECONDITION (D-02) AND THE TRAINED NET -------------------------------
    p13_require_phase11()
    basis = load_p13_basis()
    G     = isqrt(length(basis.zt.mean))

    h   = load_three_way(net_path)
    net = h.net
    hlo = h.head_log_odds     # the PERSISTED, MEASURED per-head correction. NEVER re-measured
                              # here: measuring it on this substrate would tune the reported
                              # evidence scale to the very data it is read on (D-07).

    # The conditioning width is DERIVED and asserted, never compared against a literal: Phase 11
    # may deliver a conditioning vector whose length is not one, and a hard-coded width would be
    # silently wrong with nothing thrown (T-13-51).
    n_cond   = p13_conditioning_length()
    expected = three_way_input_dim(G, n_cond)
    h.input_dim == expected || error(
        "the trained net's input width $(h.input_dim) does not equal the DERIVED " *
        "three_way_input_dim($G, $n_cond) = $expected -- the conditioning encoder is not the one " *
        "this net was trained with, or the summary row partition changed")

    # --- 3. THE PRE-REGISTERED LAMBDA SWEEP, RESOLVED FROM ITS RULES -----------------------
    P13_REAL_LAMBDA_READS_RULE === :full_phase11_ladder || error(
        "P13_REAL_LAMBDA_READS_RULE is $(P13_REAL_LAMBDA_READS_RULE), not :full_phase11_ladder")
    P13_REAL_LAMBDA_HEADLINE_RULE === :widest_rung || error(
        "P13_REAL_LAMBDA_HEADLINE_RULE is $(P13_REAL_LAMBDA_HEADLINE_RULE), not :widest_rung")
    lambda_rungs = collect(Float64, SC2_RUNGS)
    headline     = maximum(lambda_rungs)
    headline == Float64(LAMBDA_MAX) || error(
        "the realized widest rung $(headline) does not equal Phase-11's LAMBDA_MAX " *
        "$(LAMBDA_MAX). :widest_rung would resolve to a rung Phase 11 never trained at. See " *
        "spike/validation/p11_consts.jl (SC2_RUNGS / LAMBDA_MAX) before touching anything here.")
    headline == Float64(P13_REAL_LAMBDA_HEADLINE_EXPECTED) || error(
        "the realized headline rung $(headline) does not equal the frozen " *
        "P13_REAL_LAMBDA_HEADLINE_EXPECTED $(P13_REAL_LAMBDA_HEADLINE_EXPECTED) -- the " *
        "pre-registration and Phase 11's ladder have drifted apart")

    # --- 4. THE LOCKED-VALUE BANNER ---------------------------------------------------------
    if verbose
        println("="^78)
        println("PHASE-13 REAL-IMAGE ARM (D-15 AMENDED) -- QUALITATIVE, REPORTED, NOT A GATE")
        println("="^78)
        println("POSTURE (pre-registered sentinels, machine-checked):")
        println("  P13_REAL_IS_GATED         = $P13_REAL_IS_GATED")
        println("  P13_REAL_QUALITATIVE_ONLY = $P13_REAL_QUALITATIVE_ONLY")
        println("  P13_REAL_READ_ONLY        = $P13_REAL_READ_ONLY")
        println("  NOTHING BELOW CAN FAIL THE PHASE. The six committed TIFFs carry NO")
        println("  colocalization ground-truth label, so this arm shows BEHAVIOUR and never")
        println("  CORRECTNESS, and NO PASS/FAIL THRESHOLD IS DEFINED FOR ANY REAL-IMAGE")
        println("  QUANTITY. Labelled real segregation validation is deferred to the Phase 16")
        println("  blind evaluation.")
        println("SUBSTRATE (all READ from the frozen pre-registration):")
        println("  P13_REAL_CONDITIONS     = $P13_REAL_CONDITIONS")
        println("  P13_REAL_SAMPLE         = $P13_REAL_SAMPLE")
        println("  P13_REAL_CONTROL        = $P13_REAL_CONTROL")
        println("  P13_REAL_CHANNEL_PAIR   = $P13_REAL_CHANNEL_PAIR   (primary arm)")
        println("  P13_REAL_REDUNDANCY_PAIR= $P13_REAL_REDUNDANCY_PAIR   (REDUNDANCY, not a claim)")
        println("  P13_REAL_IMSIZE         = $P13_REAL_IMSIZE")
        println("  P13_REAL_ALPHA_GRID     = $P13_REAL_ALPHA_GRID")
        println("  grid                    = $G x $G   n_cond = $n_cond   input width = $expected")
        println("HYPOTHESIS BOUNDARY (Tier 2, MEASURED not chosen):")
        println("  P13_TAU                 = $P13_TAU  at reference lambda " *
                "$P13_TAU_REFERENCE_LAMBDA")
        println("LAMBDA (a real image has NO KNOWN lambda -- the rule, not a guess):")
        println("  P13_REAL_LAMBDA_READS_RULE    = $P13_REAL_LAMBDA_READS_RULE")
        println("  P13_REAL_LAMBDA_HEADLINE_RULE = $P13_REAL_LAMBDA_HEADLINE_RULE")
        println("  resolved sweep                = $(lambda_rungs)")
        println("  resolved headline rung        = $headline  (= Phase-11 LAMBDA_MAX = " *
                "$(LAMBDA_MAX); the most conservative reading, which WEAKENS the evidence, so a")
        println("                                  verdict that survives it is the credible one)")
        println("  [ASSUMED] (A11): widefield chromatic aberration plus filter-cube/stage")
        println("  repeatability on these particular microscopes is at least 1 px, which is why")
        println("  LAMBDA_MIN = $(LAMBDA_MIN) is not a credible operating point for them. THIS IS A")
        println("  DOMAIN JUDGEMENT ABOUT TYPICAL WIDEFIELD SYSTEMS, NOT A MEASUREMENT OF THESE")
        println("  MICROSCOPES, AND IT GATES NOTHING. Estimating lambda from the images is")
        println("  explicitly forbidden: that is a new estimator with its own validation burden,")
        println("  inside a phase that is not about registration.")
        println("OOD COMPARISON (P13_REAL_OOD_COMPARISON = $P13_REAL_OOD_COMPARISON):")
        println("  frozen shipped density   = $P13_REAL_OOD_SHIPPED_DENSITY")
        println("  frozen shipped threshold = $P13_REAL_OOD_SHIPPED_THRESHOLD")
        println("  The shipped bundle is a FROZEN COMPARISON REFERENCE ONLY and is NEVER an")
        println("  evidence basis; D-02's no-fallback rule governs the EVIDENCE read.")
        println("NET UNDER TEST:")
        println("  artifact                = $net_path")
        println("  head_log_odds (MEASURED, persisted) = $hlo")
        println("  consts sha256 at train time         = $(h.consts_sha)   " *
                "(P13_CONSTS_SHA.pre_amendment)")
        println("  consts sha256 of the file NOW        = $(p13_consts_sha())   " *
                "(P13_CONSTS_SHA.post_amendment)")
        println("  THE TWO consts_sha VALUES DIFFER BY DESIGN: 13-D15-AMENDMENT.md amended")
        println("  consts.jl AFTER this net was trained, and training reads no P13_REAL_*")
        println("  constant, so the amended values are ones this net never consumed. Both")
        println("  sides are asserted below against the NAMED, DATED P13_CONSTS_SHA pair --")
        println("  a third consts.jl sha256 is drift, not an amendment.")
        println("  tau at train time                   = $(h.tau)   cut = $(h.cut_variant)")
        println("STREAM (D-01): THIS ARM CONSUMES NO RNG STREAM. Six committed files, a")
        println("  deterministic transform and a deterministic forward pass -- no seed is burned")
        println("  and no reserved counter is pre-observed.")
        println("THE NAMING CORRECTION, VERBATIM FROM THE PRE-REGISTRATION:")
        println(P13_REAL_NAMING_CORRECTION)
        println("THE TARGET-SUBSTITUTION RECORD, VERBATIM FROM THE PRE-REGISTRATION:")
        println(P13_REAL_SUBSTITUTION_RECORD)
        println("  The anti-snooping seal is NOT opened by this runner: the sealed-holdout")
        println("  accessor is never called, no path under the sealed provenance tree is")
        println("  constructed, and this file asserts that about its own source text at load.")
        println("="^78)
    end

    @assert P13_REAL_IS_GATED == false "P13_REAL_IS_GATED is not false -- this runner may not become a gate"
    @assert P13_REAL_QUALITATIVE_ONLY == true "P13_REAL_QUALITATIVE_ONLY is not true"
    # THE TRAIN-TIME PROVENANCE GUARD, RE-DERIVED BY 13-D15-AMENDMENT.md section 7.3. The two
    # sides are now pinned to NAMED, DATED literals instead of merely to each other, which is
    # STRICTLY STRONGER: the old single equality would have passed silently after a retrain plus
    # an undisclosed edit. They are EXPECTED to differ here -- the net was trained under the
    # pre-amendment pre-registration and consts.jl now carries the amendment -- and that
    # divergence is asserted rather than tolerated. Training reads no P13_REAL_* constant, so the
    # amended values are ones this net never consumed. The widened `||` form is REJECTED: it
    # would accept any future drift silently. See P13_CONSTS_SHA in spike/p13/net.jl.
    @assert h.consts_sha     == P13_CONSTS_SHA.pre_amendment  "the net was NOT trained under the frozen pre-registration: its recorded consts_sha is not the pre-amendment sha256 (spike/p13/net.jl P13_CONSTS_SHA.pre_amendment)"
    @assert p13_consts_sha() == P13_CONSTS_SHA.post_amendment "spike/p13/consts.jl is at NEITHER sha256 this amendment authorises (13-D15-AMENDMENT.md section 7.3): a third value is drift, not an amendment"

    # --- 5. THE TWO OOD REFERENCES, BUILT BEFORE ANY VERDICT IS PRINTED --------------------
    verbose && println("[1/6] building the OOD references (Phase-13 in-distribution pool, and " *
                       "the frozen shipped bundle) ...")
    pool_dir = (h.meta isa NamedTuple && haskey(h.meta, :pool_dir)) ? h.meta.pool_dir : nothing
    idref    = p13_real_id_reference(pool_dir)
    shipped  = p13_real_shipped_reference()
    if verbose
        if idref.available
            println("      Phase-13 ID reference: $(idref.n_id) standardized acquisitions from " *
                    "$(idref.n_shards) shard(s)")
            println("      ID score quantiles: $(idref.scores_quantiles)")
            println("      Phase-13 operating point (COPIED shipped ID quantile " *
                    "$(P13_REALIMAGE_OOD_ID_QUANTILE)) = $(idref.threshold)")
        else
            println("      Phase-13 ID reference UNAVAILABLE: $(idref.note)")
            println("      The arm continues -- it is not a gate -- and the absence is recorded.")
        end
        if shipped.available
            println("      shipped comparison bundle = $(shipped.bundle)   recorded threshold " *
                    "= $(shipped.threshold)   reproduces the frozen constant = " *
                    "$(shipped.reproduced_frozen)")
        else
            println("      shipped comparison bundle UNAVAILABLE: $(shipped.note)")
        end
    end

    # --- helper: one (direction, lambda) row ------------------------------------------------
    # The evidence read and the OOD read are computed TOGETHER and returned on ONE row, so no
    # printed or persisted quantity can ever be separated from its OOD verdict (T-13-56).
    function sweep_rows(arm, channels)
        pair_pos = real_pair(P13_REAL_SAMPLE;  channels = channels)
        pair_neg = real_pair(P13_REAL_CONTROL; channels = channels)
        enc_pos  = encode_d01(patch_summary(build_mci(pair_pos)))
        enc_neg  = encode_d01(patch_summary(build_mci(pair_neg)))
        Zpos     = standardize_summary(enc_pos, basis.zt, basis.variant)
        Zneg     = standardize_summary(enc_neg, basis.zt, basis.variant)

        d13_pos = idref.available ? p13_real_maha(idref.nulls, Zpos) : NaN
        d13_neg = idref.available ? p13_real_maha(idref.nulls, Zneg) : NaN
        d13     = max(d13_pos, d13_neg)
        thr13   = idref.available ? Float64(idref.threshold) : NaN

        rows = NamedTuple[]
        for (dname, Zs, Zc, sname, cname) in
                (("positive as sample", Zpos, Zneg, P13_REAL_SAMPLE, P13_REAL_CONTROL),
                 ("negative as sample", Zneg, Zpos, P13_REAL_CONTROL, P13_REAL_SAMPLE))
            for lam in lambda_rungs
                b = three_way_log_bf(net, p13_encode_pair(Zs, Zc, lam), hlo)
                push!(rows, (arm = String(arm), direction = dname,
                             sample = String(sname), control = String(cname),
                             channels = Tuple(channels), lambda = Float64(lam),
                             is_headline = lam == headline,
                             logbf_coloc = b.coloc, logbf_exclusion = b.exclusion,
                             logbf_random = b.random,
                             ood_density = d13, ood_threshold = thr13,
                             ood_is_ood = idref.available ? d13 > thr13 : false,
                             ood_ratio = idref.available ? d13 / thr13 : NaN,
                             ood_density_sample = d13_pos, ood_density_control = d13_neg,
                             ood_available = idref.available))
            end
        end
        return (rows = [r for r in rows], Zpos = Zpos, Zneg = Zneg)
    end

    function print_sweep(label, rows)
        println("-"^78)
        println(label)
        println("  Every log Bayes factor is printed ON THE SAME ROW as its OOD density, its OOD")
        println("  threshold and its boolean verdict. That pairing is the whole point: a verdict")
        println("  quoted without its OOD flag would read as if the net were operating inside its")
        println("  training distribution.")
        println("  direction            lambda   logBF(C:R)   logBF(E:R)    OOD density    " *
                "OOD thr     is_ood")
        for r in rows
            println("  $(rpad(r.direction, 20)) $(rpad(r.lambda, 8)) " *
                    "$(lpad(round(r.logbf_coloc; sigdigits = 6), 12)) " *
                    "$(lpad(round(r.logbf_exclusion; sigdigits = 6), 12)) " *
                    "$(lpad(round(r.ood_density; sigdigits = 6), 14)) " *
                    "$(lpad(round(r.ood_threshold; sigdigits = 6), 11)) " *
                    "$(lpad(string(r.ood_is_ood), 10))" *
                    (r.is_headline ? "   <-- HEADLINE RUNG" : ""))
        end
        println("-"^78)
    end

    # --- 6. SECTION 1 -- THE UNMODIFIED READ, BOTH DIRECTIONS, EVERY RUNG ------------------
    verbose && println("[2/6] SECTION 1 -- the unmodified read at every pre-registered lambda " *
                       "rung, both directions ...")
    prim  = sweep_rows("primary", P13_REAL_CHANNEL_PAIR)
    sweep = prim.rows
    verbose && print_sweep("SECTION 1 -- THE UNMODIFIED REAL PAIRS (primary arm, channels " *
                           "$(P13_REAL_CHANNEL_PAIR)).", sweep)

    # --- 7. SECTION 2 -- THE SHIPPED-NET OOD COMPARISON ------------------------------------
    verbose && println("[3/6] SECTION 2 -- the shipped-net OOD comparison ...")
    ship_pos = NaN
    ship_neg = NaN
    if shipped.available
        sp = real_pair(P13_REAL_SAMPLE;  channels = P13_REAL_CHANNEL_PAIR)
        sn = real_pair(P13_REAL_CONTROL; channels = P13_REAL_CHANNEL_PAIR)
        Sp = standardize_summary(encode_d01(patch_summary(build_mci(sp))), shipped.zt, :min)
        Sn = standardize_summary(encode_d01(patch_summary(build_mci(sn))), shipped.zt, :min)
        ship_pos = p13_real_maha(shipped.density, Sp)
        ship_neg = p13_real_maha(shipped.density, Sn)
    end
    ship_pair = shipped.available ? max(ship_pos, ship_neg) : P13_REAL_OOD_SHIPPED_DENSITY
    ship_thr  = Float64(shipped.threshold)
    ship_flag = ship_pair > ship_thr
    d13_pair  = idref.available ? maximum(r.ood_density for r in sweep) : NaN
    d13_thr   = idref.available ? Float64(idref.threshold) : NaN
    d13_flag  = idref.available && d13_pair > d13_thr
    a12       = if !(idref.available && shipped.available)
                    :indeterminate
                elseif d13_flag == ship_flag
                    :agreement
                else
                    :divergence
                end
    if verbose
        println("-"^78)
        println("SECTION 2 -- THE TWO DETECTORS, SIDE BY SIDE (P13_REAL_OOD_COMPARISON = true).")
        println("  detector      density        threshold      ratio to thr   is_ood")
        println("  Phase-13      $(lpad(round(d13_pair; sigdigits = 6), 12)) " *
                "$(lpad(round(d13_thr; sigdigits = 6), 14)) " *
                "$(lpad(round(d13_pair / d13_thr; sigdigits = 4), 14)) " *
                "$(lpad(string(d13_flag), 8))")
        println("  shipped       $(lpad(round(ship_pair; sigdigits = 6), 12)) " *
                "$(lpad(round(ship_thr; sigdigits = 6), 14)) " *
                "$(lpad(round(ship_pair / ship_thr; sigdigits = 4), 14)) " *
                "$(lpad(string(ship_flag), 8))")
        println("  frozen shipped reference (pre-registration): density " *
                "$P13_REAL_OOD_SHIPPED_DENSITY vs threshold $P13_REAL_OOD_SHIPPED_THRESHOLD")
        println("  shipped bundle read = $(shipped.bundle)   reproduces the frozen constant = " *
                "$(shipped.reproduced_frozen)")
        println("  THE SHIPPED BUNDLE IS USED HERE AS A FROZEN COMPARISON REFERENCE ONLY AND IS")
        println("  NEVER AN EVIDENCE BASIS. D-02 forbids falling back to the shipped grid-8")
        println("  basis for the EVIDENCE read; every log Bayes factor above came from the")
        println("  Phase-13 net on Phase 11's registration-aware basis, with no fallback path.")
        println("  A12 OUTCOME OBSERVED: $a12")
        println("    :agreement    -- both detectors flag the same way, which CORROBORATES the")
        println("                     named limit that real microscopy at this size and these")
        println("                     statistics sits outside the simulator's training joint.")
        println("    :divergence   -- the two nets disagree, which is ITSELF A RESULT about what")
        println("                     the registration-aware basis buys on real data.")
        println("  NEITHER VERDICT IS SUPPRESSED, SOFTENED OR FIXED. An OOD flag here is the")
        println("  misspecification channel doing its job -- not a bug, and not a reason to skip")
        println("  the check.")
        println("-"^78)
    end

    # --- 8. SECTION 3 -- THE REAL ALPHA LADDER THROUGH THE NET -----------------------------
    verbose && println("[4/6] SECTION 3 -- the real alpha ladder, invariants FIRST ...")
    ladder_rows = NamedTuple[]
    ladder_figs = NamedTuple[]
    inv_records = NamedTuple[]
    alpha_star  = Dict{String,Any}()
    inv_fields  = (:bitwise_alpha0, :no_new_zeros, :intensity_conserved, :mask_invariant,
                   :mask_fraction_ok, :max_value_ok, :mbar_monotone, :missing_nonincreasing)
    invariants_ok = true

    for cond in P13_REAL_CONDITIONS
        pair = real_pair(cond; channels = P13_REAL_CHANNEL_PAIR)
        v    = verify_alpha_invariants(pair; ladder = P13_REAL_ALPHA_GRID,
                                       max_value = P13_ALPHA_MAX_VALUE_BOUND)
        ok   = all(getproperty(v, f) for f in inv_fields)
        invariants_ok &= ok
        push!(inv_records, merge((condition = String(cond), all_ok = ok), v))
        verbose && println("      $(rpad(cond, 9)) invariants $(ok ? "PASS" : "FAIL")  " *
                           "mask frac $(round(v.realized_mask_fraction; digits = 5))  " *
                           "src zeros $(v.source_zero_count)  " *
                           "max rel intensity err $(round(v.max_rel_intensity_err; sigdigits = 3))" *
                           "  realized max $(round(v.realized_max_value; sigdigits = 5))")
    end

    if !invariants_ok
        println("\n", "!"^78)
        println("INVARIANT VIOLATION -- THE LADDER CURVES ARE NOT READ AND MUST NOT BE READ.")
        println("!"^78)
        for rec in inv_records, f in inv_fields
            getproperty(rec, f) || println("  $(rec.condition): $(String(f)) FAILED")
        end
        println("")
        println("  WHY THIS STOPS THE LADDER INSTEAD OF DEGRADING IT. A FLAT LADDER CAUSED BY")
        println("  DELETED PIXELS IS INDISTINGUISHABLE FROM A NULL RESULT. `_exclude_zero` drops")
        println("  any pixel where either channel is zero, so a transform that pushed masked ch2")
        println("  pixels to hard zero would REMOVE them from the patch correlation rather than")
        println("  anti-correlate them, and the ladder would read flat for a reason that has")
        println("  nothing whatsoever to do with segregation. What HAS been computed is persisted")
        println("  below so the violation itself is on the record.")
        println("!"^78, "\n")
    else
        verbose && println("      ALL invariants hold on BOTH real conditions; the curves may be " *
                           "read.")
        for cond in P13_REAL_CONDITIONS
            RL   = real_alpha_ladder(cond; channels = P13_REAL_CHANNEL_PAIR,
                                     grid = P13_REAL_ALPHA_GRID)
            L    = alpha_ladder(real_pair(cond; channels = P13_REAL_CHANNEL_PAIR);
                                ladder = P13_REAL_ALPHA_GRID,
                                max_value = P13_ALPHA_MAX_VALUE_BOUND)
            ctrl = real_pair(cond == P13_REAL_SAMPLE ? P13_REAL_CONTROL : P13_REAL_SAMPLE;
                             channels = P13_REAL_CHANNEL_PAIR)
            Zc   = standardize_summary(encode_d01(patch_summary(build_mci(ctrl))),
                                       basis.zt, basis.variant)
            d_c  = idref.available ? p13_real_maha(idref.nulls, Zc) : NaN

            xs, lbe, lbc = Float64[], Float64[], Float64[]
            for (r, rung) in enumerate(L.pairs)
                enc  = encode_d01(patch_summary(build_mci(rung)))
                Zs   = standardize_summary(enc, basis.zt, basis.variant)
                b    = three_way_log_bf(net, p13_encode_pair(Zs, Zc, headline), hlo)
                d_s  = idref.available ? p13_real_maha(idref.nulls, Zs) : NaN
                dens = idref.available ? max(d_s, d_c) : NaN
                thr  = idref.available ? Float64(idref.threshold) : NaN
                s    = RL.rungs[r]
                push!(xs, Float64(s.alpha)); push!(lbe, b.exclusion); push!(lbc, b.coloc)
                push!(ladder_rows, (condition = String(cond),
                                    channels = Tuple(P13_REAL_CHANNEL_PAIR),
                                    lambda = headline, alpha = Float64(s.alpha),
                                    logbf_coloc = b.coloc, logbf_exclusion = b.exclusion,
                                    mbar = s.mbar, n_missing = s.n_missing,
                                    max_value = s.max_value, sum_ratio = s.sum_ratio,
                                    n_zero = s.n_zero,
                                    ood_density = dens, ood_threshold = thr,
                                    ood_is_ood = idref.available ? dens > thr : false,
                                    ood_available = idref.available))
            end
            # alpha_star: the SMALLEST rung at which the mean log BF(E:R) first exceeds 0.
            # `nothing` is a legitimate, reportable outcome and is NEVER interpolated away -- the
            # ladder's own resolution IS the resolution of this statistic.
            star = nothing
            for i in eachindex(lbe)
                if isfinite(lbe[i]) && lbe[i] > 0.0
                    star = xs[i]
                    break
                end
            end
            alpha_star[String(cond)] = star
            push!(ladder_figs, (condition = String(cond), alpha = xs,
                                logbf_exclusion = lbe, logbf_coloc = lbc, alpha_star = star))
        end

        if verbose
            for cond in P13_REAL_CONDITIONS
                rows = filter(r -> r.condition == String(cond), ladder_rows)
                println("-"^78)
                println("SECTION 3 -- REAL ALPHA LADDER, condition \"$cond\", at the headline " *
                        "lambda $headline.")
                println("  alpha    logBF(E:R)    logBF(C:R)       m-bar   n_missing   " *
                        "max_value     is_ood")
                for r in rows
                    println("  $(rpad(r.alpha, 8)) " *
                            "$(lpad(round(r.logbf_exclusion; sigdigits = 6), 12)) " *
                            "$(lpad(round(r.logbf_coloc; sigdigits = 6), 13)) " *
                            "$(lpad(round(r.mbar; sigdigits = 5), 11)) " *
                            "$(lpad(r.n_missing, 11)) " *
                            "$(lpad(round(r.max_value; sigdigits = 5), 11)) " *
                            "$(lpad(string(r.ood_is_ood), 10))")
                end
                st = alpha_star[String(cond)]
                println("  alpha_star_real(\"$cond\") = " *
                        (st === nothing ? "NOTHING -- the ladder never crosses zero, which is a " *
                                          "reportable outcome and is not interpolated away" :
                                          "$st  (a RUNG of the frozen ladder, never interpolated)"))
                println("-"^78)
            end
            println("THIS IS NOT THE SIMULATED ARM, AND THE TWO MUST NEVER BE AVERAGED.")
            println("  On SIMULATED substrate the ladder runs RANDOM -> EXCLUSION by")
            println("  construction: both acquisitions are drawn at rho_true = 0. On THIS")
            println("  substrate alpha = 0 is a MODERATELY COLOCALIZED pair -- measured m-bar")
            println("  $(P13_REAL_ANCHOR_MBAR) -- so the ladder runs MODERATELY COLOCALIZED ->")
            println("  NEAR-EXCLUSION. Plan 13-13 reports the simulated arm separately. They are")
            println("  different experiments, both informative, and a merged curve would be a")
            println("  statement about neither.")
            println("-"^78)
            # --- 9. SECTION 4 -- THE SCOPE LIMIT ON THE LADDER ---------------------------------
            println("SECTION 4 -- THE SCOPE LIMIT OF THIS LADDER, STATED SO IT IS NOT OVER-READ.")
            println("  The real ladder spans rho_true in roughly [-0.44, +0.33] -- a DENSELY")
            println("  COVERED INTERIOR region of GHAT_RHO_KNOTS, comfortably away from the")
            println("  +/-0.99 clamp atoms. It therefore probes THE SIGN TRANSITION AND THE")
            println("  NEAR-EXCLUSION REGIME, and it DOES NOT PROBE THE DEEP-EXCLUSION TAIL where")
            println("  the prior-atom asymmetry lands hardest. A reader who assumes otherwise")
            println("  will over-read the result.")
            println("  HONEST CAVEAT ON THE ENDPOINT: at alpha = 1 the mask region is not EMPTY,")
            println("  it is AT the background floor. The endpoint is 'ch2 reduced to background")
            println("  wherever ch1 has objects', not 'ch2 absent' -- do not write 'fully")
            println("  disjoint' without that qualifier.")
            println("-"^78)
        end
    end

    # --- 10. SECTION 5 -- THE REDUNDANCY ARM -----------------------------------------------
    verbose && println("[5/6] SECTION 5 -- the REDUNDANCY arm at " *
                       "P13_REAL_REDUNDANCY_PAIR = $(P13_REAL_REDUNDANCY_PAIR) ...")
    red = sweep_rows("redundancy", P13_REAL_REDUNDANCY_PAIR)
    if verbose
        print_sweep("SECTION 5 -- REDUNDANCY (channels $(P13_REAL_REDUNDANCY_PAIR)). THIS IS " *
                    "NOT A SEPARATE CLAIM.", red.rows)
        println("  Labelled REDUNDANCY on purpose: it re-reads the SAME two specimens through a")
        println("  second pre-registered channel pair, so it can corroborate or contradict the")
        println("  primary arm but adds no independent specimen and no new evidence about")
        println("  colocalization. Note that BOTH pre-registered pairs include channel 1, the")
        println("  nuclear counterstain -- a limit of the frozen pre-registration, named here")
        println("  rather than repaired by choosing a new pair after the fact.")
        println("-"^78)
    end

    # --- 11. THE READ-ONLY DIGEST, AFTER ---------------------------------------------------
    digest_after = verify_real_readonly_digest()
    verbose && println("[6/6] read-only digest AFTER  = $digest_after")

    # --- 12. PERSIST ATOMICALLY, BEFORE THE DIGEST COMPARISON CAN THROW --------------------
    ood_phase13 = (available = idref.available,
                   source = "Phase-13 labelled training pool, standardized through the FROZEN " *
                            "Phase-11 zt (D-02); density-Mahalanobis on the continuous rows only",
                   n_id = idref.n_id, n_shards = idref.n_shards,
                   pool_dir = idref.pool_dir === nothing ? "" : String(idref.pool_dir),
                   id_quantile = P13_REALIMAGE_OOD_ID_QUANTILE,
                   ridge = P13_REALIMAGE_OOD_RIDGE,
                   density_sample = idref.available ? sweep[1].ood_density_sample : NaN,
                   density_control = idref.available ? sweep[1].ood_density_control : NaN,
                   density_pairmax = d13_pair, threshold = d13_thr,
                   ratio_to_threshold = d13_pair / d13_thr, is_ood = d13_flag,
                   score_quantiles = idref.scores_quantiles === nothing ?
                                     (q50 = NaN, q90 = NaN, q95 = NaN, q99 = NaN, qmax = NaN) :
                                     idref.scores_quantiles,
                   note = idref.note)
    ood_shipped = (available = shipped.available,
                   bundle = shipped.bundle === nothing ? "" : String(shipped.bundle),
                   density_sample = ship_pos, density_control = ship_neg,
                   density_pairmax = ship_pair, threshold = ship_thr,
                   ratio_to_threshold = ship_pair / ship_thr, is_ood = ship_flag,
                   frozen_density_const = P13_REAL_OOD_SHIPPED_DENSITY,
                   frozen_threshold_const = P13_REAL_OOD_SHIPPED_THRESHOLD,
                   reproduced_frozen = shipped.reproduced_frozen,
                   role = "FROZEN COMPARISON REFERENCE ONLY, never an evidence basis (D-02)",
                   note = shipped.note)

    verbose && println("      persisting the reported artifact ...")
    _p13_real_save_report(report_path;
        # --- the posture, machine-checked and carried with every quotable number ---
        is_gated = P13_REAL_IS_GATED,
        qualitative_only = P13_REAL_QUALITATIVE_ONLY,
        read_only = P13_REAL_READ_ONLY,
        naming_correction = P13_REAL_NAMING_CORRECTION,
        substitution_record = P13_REAL_SUBSTITUTION_RECORD,
        physical_anchor_deferred_to = P13_PHYSICAL_ANCHOR_DEFERRED_TO,
        # --- Section 1: one row per (direction, lambda), each carrying its own OOD verdict ---
        lambda_sweep = [r for r in sweep],
        lambda_rungs = lambda_rungs,
        lambda_headline = headline,
        lambda_reads_rule = P13_REAL_LAMBDA_READS_RULE,
        lambda_headline_rule = P13_REAL_LAMBDA_HEADLINE_RULE,
        lambda_headline_expected = P13_REAL_LAMBDA_HEADLINE_EXPECTED,
        phase11_lambda_range = collect(Float64, P13_PHASE11_LAMBDA_RANGE),
        # --- Section 2: the two detectors side by side ---
        ood_phase13 = ood_phase13,
        ood_shipped = ood_shipped,
        ood_a12_outcome = String(a12),
        # --- Section 3: the real alpha ladder, its invariants and its crossing points ---
        alpha_ladder_real = [r for r in ladder_rows],
        alpha_star_real = (positive = get(alpha_star, "positive", nothing),
                           negative = get(alpha_star, "negative", nothing)),
        alpha_invariants_real = [r for r in inv_records],
        alpha_invariants_ok = invariants_ok,
        alpha_grid = collect(Float64, P13_REAL_ALPHA_GRID),
        alpha_never_averaged_with = "the :simulated arm of plan 13-13 -- alpha = 0 means a " *
                                    "RANDOM pair there and a MODERATELY COLOCALIZED pair here",
        # --- Section 5: the redundancy arm, labelled, never a separate claim ---
        redundancy_sweep = [r for r in red.rows],
        redundancy_channels = collect(Int, P13_REAL_REDUNDANCY_PAIR),
        # --- the substrate provenance of both fixtures ---
        provenance_sample = real_provenance(P13_REAL_SAMPLE; channels = P13_REAL_CHANNEL_PAIR),
        provenance_control = real_provenance(P13_REAL_CONTROL; channels = P13_REAL_CHANNEL_PAIR),
        grid_provenance = real_grid_provenance(; G = G),
        anchor_mbar = P13_REAL_ANCHOR_MBAR,
        # --- the read-only proof ---
        digest_before = digest_before,
        digest_after = digest_after,
        digest_source = "sha256 of `git ls-files -s test/test_images`",
        # --- EVERY pre-registered value this run READ, written INTO the artifact ---
        P13_REAL_CONDITIONS = collect(String.(P13_REAL_CONDITIONS)),
        P13_REAL_SAMPLE = P13_REAL_SAMPLE, P13_REAL_CONTROL = P13_REAL_CONTROL,
        P13_REAL_CHANNEL_PAIR = collect(Int, P13_REAL_CHANNEL_PAIR),
        P13_REAL_IMSIZE = collect(Int, P13_REAL_IMSIZE),
        P13_REAL_OOD_COMPARISON = P13_REAL_OOD_COMPARISON,
        P13_TAU = P13_TAU, P13_TAU_REFERENCE_LAMBDA = P13_TAU_REFERENCE_LAMBDA,
        P13_ALPHA_GATED = P13_ALPHA_GATED,
        P13_ITERATION_ALLOWANCE = P13_ITERATION_ALLOWANCE,
        # --- provenance: the net, the basis, the frozen pre-registration ---
        rng_stream = "NONE -- this arm consumes no RNG stream (deterministic substrate, " *
                     "deterministic transform, deterministic forward pass)",
        net_path = net_path, net_head_log_odds = hlo,
        net_consts_sha = h.consts_sha, net_meta = h.meta,
        net_input_dim = h.input_dim, derived_input_dim = expected,
        grid = G, n_cond = n_cond,
        consts_sha = p13_consts_sha(),
        consts_git_blob_sha = p13_real_blob_sha(joinpath(@__DIR__, "consts.jl")),
        phase11_net = P13_PHASE11_NET, phase11_sha = p13_phase11_sha(),
        assumed_a11 = "[ASSUMED] widefield chromatic aberration plus filter-cube/stage " *
                      "repeatability on these particular microscopes is at least 1 px. A DOMAIN " *
                      "JUDGEMENT ABOUT TYPICAL WIDEFIELD SYSTEMS, NOT A MEASUREMENT OF THESE " *
                      "MICROSCOPES, AND IT GATES NOTHING.",
        note = "QUALITATIVE, REPORTED, NOT GATED (P13_REAL_IS_GATED = false, " *
               "P13_REAL_QUALITATIVE_ONLY = true). The six committed TIFFs carry NO " *
               "colocalization ground-truth label, so this arm shows BEHAVIOUR and never " *
               "CORRECTNESS. No pass/fail threshold is defined for any real-image quantity. " *
               "Labelled real segregation validation is deferred to the Phase 16 blind " *
               "evaluation, which is what the sealed holdout exists for and why Phase 13 did " *
               "not open it.",
        generated = string(Dates.now(Dates.UTC)) * "Z")
    verbose && println("      persisted reported artifact -> $report_path")

    # --- 13. THE FIGURE. A SCRIPT ARTIFACT, NEVER A GATE ASSERTION -------------------------
    figout = nothing
    try
        figout = p13_real_figure(sweep, ladder_figs, headline; path = fig_path)
        verbose && println("      figure -> $figout")
    catch err
        @warn "the figure could not be rendered; the reported artifact is already persisted" err
    end

    # --- 14. THE READ-ONLY COMPARISON, AFTER THE ARTIFACT IS SAFELY ON DISK ----------------
    # Deliberately AFTER the persist: the digest pair is a property of the run and belongs in the
    # record whatever it says, and an artifact lost to an assertion would be the worse failure.
    @assert digest_before == digest_after "READ-ONLY BREACH (T-13-50): the git-index digest of test/test_images changed during the run.\n  before = $digest_before\n  after  = $digest_after"
    verbose && println("READ-ONLY PROOF: the digests match; the committed fixtures are " *
                       "byte-identical before and after.")
    verbose && println("total wall clock: $(round((time() - t_start) / 60; digits = 2)) min")

    return (is_gated = false, qualitative_only = true,
            lambda_sweep = sweep, lambda_headline = headline,
            ood_phase13 = ood_phase13, ood_shipped = ood_shipped, a12 = a12,
            alpha_ladder_real = ladder_rows, alpha_star_real = alpha_star,
            invariants_ok = invariants_ok, redundancy_sweep = red.rows,
            digest_before = digest_before, digest_after = digest_after,
            report_path = report_path, figure = figout)
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so `include`-ing this file can NEVER trigger the reported run. Run it deliberately:
#
#     julia --project=spike spike/p13/run_p13_realimage.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
