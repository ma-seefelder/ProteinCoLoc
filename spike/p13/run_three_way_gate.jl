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

# spike/p13/run_three_way_gate.jl --- REPORTED D-12/D-13 gate: LOCKED consts, RESERVED stream.
#
# =============================================================================================
# ANTI-SNOOPING BANNER. READ BEFORE CHANGING ANY NUMBER IN THIS FILE.
# =============================================================================================
# EVERY THRESHOLD THIS RUNNER SCORES AGAINST WAS COMMITTED TO `spike/p13/consts.jl` BEFORE THE
# NET EXISTED. The AUC floors, the ECE band and its bin count, the evaluation size, the minimum
# per-head evaluation count, the vacuous-pass AUC floor, the gate statistic, the confusion rule
# and the not-gated status of the continuity check are all READ from that frozen file by name.
# Not one of them is a literal here, and the test-side acceptance criteria assert that, so a
# value cannot be "temporarily" inlined and left behind.
#
# ROADMAP SC2 IS AMENDED BY `.planning/phases/13-three-hypothesis-amortized-bayes-factor/
# 13-SC2-AMENDMENT.md`, which was written and frozen while NO Phase-13 result of any kind
# existed. ANY REPORT THAT CITES THIS RESULT MUST CITE THE ORIGINAL ROADMAP CRITERION ALONGSIDE
# IT -- reporting only the amended criterion is exactly the failure that document exists to
# prevent, and the phase report (plan 13-14) is required to print both side by side.
#
# AN HONEST SHORTFALL IS A PHASE-13 FINDING, NEVER A LICENCE TO ACT. It is not a reason to
# loosen a locked threshold, to reseed, to redraw the evaluation set, to retrain, or to
# re-render the figure with different scaling. `P13_ITERATION_ALLOWANCE = 1` has exactly ONE
# pre-declared trigger, written into `P13_ITERATION_TRIGGER` before any result existed, and it
# is a stratification switch -- never a threshold relaxation. Phase 7 was amended TWICE after
# seeing results, and the credibility cost of that is why this paragraph is here.
#
# THE ORDER OF OPERATIONS IS THE POINT: read the thresholds, compute the statistics, PERSIST THE
# ARTIFACT, print the headline, then assert. The artifact is written before anything can throw,
# so a FAILING gate still leaves a complete, self-describing report on disk.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY AND NOT A TEST. It composes surfaces that already exist
# and are already tested. Including it does nothing; it must be run deliberately:
#
#     julia --project=spike -t auto spike/p13/run_three_way_gate.jl
#
# NO KDE, NO QUADGK. The retired KDE Bayes-factor baseline is invalid past |logBF| ~ log(L) ~ 6.9
# (`docs/amortized.md` named limit 3), its attrition is 100% baseline-side, and it is not even
# reachable from the lean spike environment. Every AUC here is the in-repo hand-rolled tie-aware
# `roc_auc` (`spike/validation/ood.jl`); no package is installed.
#
# CPU-ONLY (CLAUDE.md, D-01). CUDA is never imported. DECOUPLING: spike-local; `src/` is reached
# only READ-ONLY and transitively, and step 0 PROVES at run time that it is byte-unchanged.

using Test
using JLD2
using Statistics
using Dates
using SHA
import StatsBase
import Random123: Philox4x

# ORDER MATTERS. The datagen surface first -- it carries the whole Phase-11 precondition preamble
# (preconditions.jl), the pre-registration (consts.jl), the label surface (labels.jl) and the net
# (net.jl). Then the Phase-13 result/calibration surface, then the hand-rolled AUC, then the
# SHIPPED-analog binary NRE read surface (continuity, reported only), then the frozen Phase-11
# rung ladder, then the figure surface. All guarded, the house guarded-include idiom.
isdefined(@__MODULE__, :stratify_by_class)    || include(joinpath(@__DIR__, "datagen.jl"))
isdefined(@__MODULE__, :p13_calibration_meta) || include(joinpath(@__DIR__, "result.jl"))
isdefined(@__MODULE__, :roc_auc)              ||
    include(joinpath(@__DIR__, "..", "validation", "ood.jl"))
isdefined(@__MODULE__, :amortized_log_bf)     ||
    include(joinpath(@__DIR__, "..", "validation", "bf.jl"))
isdefined(@__MODULE__, :SC2_RUNGS)            ||
    include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
isdefined(@__MODULE__, :plot_roc)             ||
    include(joinpath(@__DIR__, "..", "validation", "figures.jl"))

if !isdefined(@__MODULE__, :P13_GATE_REPORT_PATH)
    # NOTE THE PATH DEPTH: from spike/p13/ the repo root is TWO levels up, not three.
    "Repository root, for the run-time decoupling proof (step 0)."
    const P13_GATE_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

    "The reported D-12/D-13 gate artifact. Written BEFORE the gate can throw."
    const P13_GATE_REPORT_PATH = joinpath(@__DIR__, "three_way_gate_report.jld2")

    "The trained three-way evidence net (plan 13-11). There is no substitute and no fallback."
    const P13_GATE_NET_PATH = joinpath(@__DIR__, "three_way_net.jld2")

    """
        P13_GATE_BINARY_NET

    The frozen SPIKE-LANE binary NRE read for the CONTINUITY OBSERVATION ONLY (reported, never
    gated). It carries its OWN frozen `zt`, it is the artifact `spike/validation/run_bf.jl`
    reports against, and it is the same shuffled-theta `RatioEstimator` construction as the
    shipped bundle.

    THE SHIPPED `artifacts/grid_8/ratio_8.jld2` IS DELIBERATELY NOT READ HERE. `consts.jl`
    section I2 sanctions exactly ONE read of the shipped grid-8 bundle -- as a frozen OOD
    COMPARISON REFERENCE -- and this is not that read. `src/amortized/bf.jl` is likewise not
    included: it imports KernelDensity and QuadGK, neither of which is in the lean spike
    environment, which is the same reason `spike/validation/bf.jl` exists.
    """
    const P13_GATE_BINARY_NET = joinpath(@__DIR__, "..", "validation", "trained_ratio.jld2")

    "The confusion/ROC figure. A SCRIPT ARTIFACT, never a gate assertion."
    const P13_GATE_FIG_PATH = joinpath(@__DIR__, "..", "figures", "p13_confusion.png")

    """
        P13_GATE_RHO_BANDS

    Upper edges of the |rho_sample| bands the EXCLUSION-class AUC is ADDITIONALLY reported over.

    THESE ARE REPORTING BANDS, NOT THRESHOLDS. No bar of any kind is attached to a band; the only
    bar in play is the frozen `P13_AUC_FLOOR_EXCLUSION`, applied unchanged inside each band. They
    exist because `P13_ITERATION_TRIGGER` is worded specifically about the deep-exclusion tail
    near the -0.99 atom, so the trigger can only be evaluated MECHANICALLY if the exclusion AUC
    is resolved as a function of |rho| rather than argued from the pooled number.
    """
    const P13_GATE_RHO_BANDS = (0.5, 0.75, 0.9, 1.0)

    "The class order every table, matrix row and matrix column in this runner uses."
    const P13_GATE_CLASS_ORDER = (:exclusion, :random, :coloc)
end

# =============================================================================================
# The RESERVED evaluation stream (D-01)
# =============================================================================================

"""
    p13_gate_rng(idx::Integer) -> Philox4x

The per-item keyed RNG of the REPORTED evaluation set:
`Philox4x((P13_DEV_SEED xor P13_SALT xor P13_GATE_COUNTER, idx))`.

IT IS THE `p13_datagen_rng` CONSTRUCTION AT THE GATE COUNTER, AND THAT IS THE WHOLE POINT. The
counter goes in the KEY word and the per-item GLOBAL INDEX in the counter word, so a column is a
pure function of its index and the evaluation set is byte-identical for any thread count. Because
`P13_GATE_COUNTER != P13_DATAGEN_COUNTER`, the resulting Philox key differs from the training
pool's, so the evaluation set is DISJOINT from the 48,000-item pool the net was trained on. An
evaluation set drawn from the training counter would report memorization as discrimination; the
assertion in `main` makes that structural rather than a matter of trust.

The evaluation set therefore resolves to its OWN content-hash directory if it is ever cached.
That is the CORRECT outcome, not a cache miss.
"""
p13_gate_rng(idx::Integer) =
    Philox4x(UInt64, (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_GATE_COUNTER),
                      UInt64(idx)))

# The key words that must differ, asserted rather than assumed. "Obviously distinct" is how two
# lanes end up sharing a Philox key, and a shared key between the training pool and the reported
# evaluation set is the single most damaging thing that could silently happen to this gate.
@assert (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_GATE_COUNTER)) !=
        (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_DATAGEN_COUNTER)) "the gate key word collides with the training-pool key word"
@assert (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_GATE_COUNTER)) !=
        (UInt64(P13_FIX_SEED) ⊻ P13_SALT ⊻ UInt64(P13_FIXTURE_COUNTER)) "the gate key word collides with the fixture family"

# =============================================================================================
# Small local helpers
# =============================================================================================

"""
    _p13_gate_save_report(path; kwargs...) -> String

Atomically persist the reported gate artifact: write `path * ".tmp"`, REOPEN to integrity-check a
required key, then `mv(...; force = true)`. The `_save_bf_report` idiom
(`spike/validation/run_bf.jl:59-69`). A crash mid-write leaves a discardable `.tmp`, never a torn
artifact a later reader would happily believe.
"""
function _p13_gate_save_report(path; kwargs...)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "auc_coloc") "_p13_gate_save_report: integrity check failed ($tmp)"
        @assert haskey(f, "confusion") "_p13_gate_save_report: integrity check failed ($tmp)"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    p13_blob_sha(path) -> String

The GIT BLOB sha1 of `path`, i.e. exactly what `git hash-object <path>` prints. Recorded so the
artifact points INTO the audit trail: a reader can run `git cat-file -p <sha>` to get the frozen
pre-registration text out of history, or `git log --find-object=<sha>` to get the commit that
introduced it.
"""
function p13_blob_sha(path)
    bytes = read(path)
    ctx   = SHA.SHA1_CTX()
    SHA.update!(ctx, Vector{UInt8}("blob $(length(bytes))\0"))
    SHA.update!(ctx, bytes)
    return bytes2hex(SHA.digest!(ctx))
end

"""
    p13_gate_auc(neg_scores, pos_scores) -> Float64

The AUC only, from the in-repo hand-rolled tie-aware `roc_auc` (`spike/validation/ood.jl:309`).
No package is installed and no ROC library is added; `roc_auc` returns the Mann-Whitney U
statistic exactly, which is robust to the threshold degeneracies a well-separated head produces.
"""
p13_gate_auc(neg_scores, pos_scores) = roc_auc(neg_scores, pos_scores)[3]

"""
    p13_to_binary_basis(Z_p11, basis, zt_bin) -> Vector{Float32}

Move ONE Phase-11-standardized summary into the SHIPPED-ANALOG BINARY NET'S OWN frozen basis, for
the continuity observation only.

`standardize_summary` (`spike/npe/infer.jl:76-84`) is an AFFINE map on the CONTINUOUS rows and
leaves the mask rows untouched, so it is exactly invertible: `StatsBase.reconstruct` recovers the
raw summary from the Phase-11 rows, and the binary net's own frozen `zt` then puts that raw
summary into ITS input space. Nothing is re-fitted, nothing is estimated, and the mask block is
carried through unchanged by both transforms.

THIS IS WHY THE CONTINUITY NUMBER IS AN OBSERVATION AND NOT A CRITERION: the round trip proves the
two nets can be shown the same acquisition, and it does NOT make their outputs the same object.
See the continuity section's printed reasons.
"""
function p13_to_binary_basis(Z_p11::AbstractVector, basis, zt_bin)
    n  = length(Z_p11)
    nc = n ÷ 2
    raw = Vector{Float64}(undef, n)
    raw[1:nc]      = vec(StatsBase.reconstruct(basis.zt,
                                               reshape(Float64.(Z_p11[1:nc]), :, 1)))
    raw[nc+1:end]  = Float64.(Z_p11[nc+1:end])      # mask rows: bypassed by every transform
    return standardize_summary(raw, zt_bin, basis.variant)
end

# =============================================================================================
# The reported run
# =============================================================================================

"""
    main(; m = P13_GATE_M, min_per_head = P13_MIN_EVAL_PER_HEAD, ...) -> NamedTuple

The reported Phase-13 D-12/D-13 gate, end to end and once.

`m` and `min_per_head` default to the FROZEN constants. Passing anything else puts the run in
SMOKE mode: it writes to a `_smoke` artifact path, prints a banner saying so, and DOES NOT RUN
THE GATE ASSERTIONS -- so a load-check can never be mistaken for, or quoted as, a verdict.
"""
function main(; m::Integer = P13_GATE_M,
              min_per_head::Integer = P13_MIN_EVAL_PER_HEAD,
              net_path = P13_GATE_NET_PATH,
              report_path = P13_GATE_REPORT_PATH,
              fig_path = P13_GATE_FIG_PATH,
              binary_net_path = P13_GATE_BINARY_NET,
              verbose::Bool = true)

    reported = (m == P13_GATE_M) && (min_per_head == P13_MIN_EVAL_PER_HEAD)
    if !reported
        report_path = replace(report_path, ".jld2" => "_smoke.jld2")
        fig_path    = replace(fig_path, ".png" => "_smoke.png")
    end

    # --- 0. THE DECOUPLING PROOF, AT RUN TIME (D-01) ---------------------------------------
    # Asserted while the run happens, not only in review: a long reported run that started on a
    # clean tree and finished on a dirty one would be a repudiation hole (T-13-39). The
    # `spike/test/test_p13_result.jl:203-204` precedent, promoted to a runner precondition.
    src_clean = success(Cmd(`git diff --quiet HEAD -- src`; dir = P13_GATE_REPO_ROOT))
    env_clean = success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                            dir = P13_GATE_REPO_ROOT))
    @assert src_clean "D-01 decoupling breach: `git diff --quiet HEAD -- src` FAILED -- src/ is not byte-unchanged, so this run may not be reported"
    @assert env_clean "D-01 decoupling breach: spike/Project.toml or spike/Manifest.toml is modified -- the spike environment moved under the run"
    verbose && println("STEP 0 (D-01): src/ byte-unchanged = $src_clean   " *
                       "spike Project/Manifest byte-unchanged = $env_clean")

    # --- 1. THE HARD PRECONDITION (D-02) ----------------------------------------------------
    # A missing Phase-11 research NPE is a BLOCK on the phase, not a skip and not a fallback to
    # the shipped grid-8 basis.
    p13_require_phase11()
    basis = load_p13_basis()
    G     = isqrt(length(basis.zt.mean))

    h   = load_three_way(net_path)
    net = h.net
    hlo = h.head_log_odds     # the PERSISTED, MEASURED per-head correction. Never re-measured
                              # here: measuring it on the evaluation set would tune the reported
                              # evidence scale to the data it is scored against (D-07).

    # --- 2. THE LOCKED-THRESHOLD BANNER -----------------------------------------------------
    if verbose
        println("="^78)
        println("PHASE-13 REPORTED THREE-WAY GATE (D-12/D-13) -- LOCKED consts, RESERVED stream")
        println("="^78)
        reported || println("*** SMOKE MODE -- NOT A REPORTED RUN, THE GATE IS NOT ASSERTED ***")
        println("AMENDED CRITERION: ROADMAP Phase-13 SC2 is amended by 13-SC2-AMENDMENT.md,")
        println("  frozen while NO Phase-13 result existed. Any report citing this result MUST")
        println("  cite the original ROADMAP criterion alongside it (amendment section 6).")
        println("GATED (D-12, discrimination):")
        println("  P13_GATE_M              = $P13_GATE_M            (requested m = $m)")
        println("  P13_AUC_FLOOR_COLOC     = $P13_AUC_FLOOR_COLOC")
        println("  P13_AUC_FLOOR_EXCLUSION = $P13_AUC_FLOOR_EXCLUSION")
        println("  P13_MIN_EVAL_PER_HEAD   = $P13_MIN_EVAL_PER_HEAD          " *
                "(requested min_per_head = $min_per_head)")
        println("GATED (D-13, calibration):")
        println("  P13_GATE_STATISTIC      = $P13_GATE_STATISTIC          " *
                "(ECE GATES; MCE is REPORTED and NEVER gated)")
        println("  P13_ECE_GREEN           = $P13_ECE_GREEN")
        println("  P13_ECE_YELLOW          = $P13_ECE_YELLOW")
        println("  P13_ECE_NBINS           = $P13_ECE_NBINS")
        println("  P13_VACUOUS_AUC_FLOOR   = $P13_VACUOUS_AUC_FLOOR")
        println("  P13_TOST_REQUIRED       = $P13_TOST_REQUIRED        " *
                "(the traffic light is a BAND, not a p-value: D-14)")
        println("REPORTED, NOT GATED:")
        println("  P13_CONFUSION_RULE      = $P13_CONFUSION_RULE")
        println("  P13_CONTINUITY_GATED    = $P13_CONTINUITY_GATED")
        println("  lambda ladder           = $SC2_RUNGS  (Phase-11 SC2 rungs, READ not chosen)")
        println("HYPOTHESIS BOUNDARY (Tier 2, MEASURED not chosen):")
        println("  P13_TAU                 = $P13_TAU  at reference lambda " *
                "$P13_TAU_REFERENCE_LAMBDA")
        println("  P13_TAU_MEASURED_AUC    = $P13_TAU_MEASURED_AUC  (bar $P13_TAU_AUC)")
        println("  P13_TAU_PROBE_SHA       = $P13_TAU_PROBE_SHA")
        println("  P13_TAU_PROBE_ARTIFACT  = $P13_TAU_PROBE_ARTIFACT")
        println("  P13_CUT_VARIANT         = $P13_CUT_VARIANT   (labels from the D-05 two-factor")
        println("                            cut on the rho_sample LEVEL x the control contrast,")
        println("                            NEVER from delta_rho)")
        println("ITERATION DISCIPLINE:")
        println("  P13_ITERATION_ALLOWANCE = $P13_ITERATION_ALLOWANCE  " *
                "(one pre-declared trigger; this run does NOT spend it)")
        println("NET UNDER TEST:")
        println("  artifact                = $net_path")
        println("  head_log_odds (MEASURED, persisted) = $hlo")
        println("  consts sha256 at train time         = $(h.consts_sha)")
        println("  tau at train time                   = $(h.tau)   " *
                "cut = $(h.cut_variant)")
        println("STREAM (D-01):")
        println("  P13_DEV_SEED            = $(repr(UInt64(P13_DEV_SEED))) " *
                "at counter P13_GATE_COUNTER = $P13_GATE_COUNTER")
        println("  training pool counter   = P13_DATAGEN_COUNTER = $P13_DATAGEN_COUNTER " *
                "(DISJOINT)")
        println("="^78)
    end

    # --- 3. STREAM HYGIENE, ASSERTED BEFORE ANYTHING IS DRAWN (T-13-37) --------------------
    @assert P13_GATE_COUNTER != P13_DATAGEN_COUNTER "the reported evaluation set must not ride the training pool's counter"
    @assert !(UInt64(P13_DEV_SEED) in _p13_forbidden()) "P13_DEV_SEED is a forbidden (pre-observed) seed"
    @assert UInt64(P13_DEV_SEED) != UInt64(P13_FIX_SEED) "the reported stream must differ from the fixture stream"
    @assert h.consts_sha == p13_consts_sha() "the frozen pre-registration changed since the net was trained: consts.jl sha256 does not match the artifact's"

    # --- 4. THE FRESH EVALUATION SET --------------------------------------------------------
    verbose && println("[1/6] generating the FRESH evaluation set: $m labelled pairs at " *
                       "P13_GATE_COUNTER …")
    t0  = time()
    strat = stratify_by_class(m, basis; rng_for = p13_gate_rng, verbose = verbose)
    items = strat.items
    verbose && println("      done in $(round((time() - t0) / 60; digits = 2)) min " *
                       "(overhead $(round(strat.n_draws / m; digits = 3))x)")

    classes = [it.class for it in items]
    n_excl  = count(==(EXCLUSION), classes)
    n_rand  = count(==(RANDOM), classes)
    n_col   = count(==(COLOC), classes)
    realized = class_masses(classes)
    if verbose
        println("      realized class counts  E=$n_excl  R=$n_rand  C=$n_col")
        println("      realized frequencies   $realized  (target $P13_TARGET_CLASS_FREQ)")
    end

    # Each head is scored ONLY on its own class pair (D-11), so its restricted evaluation set is
    # its own positives plus the shared RANDOM negatives.
    n_head_coloc = n_col + n_rand
    n_head_excl  = n_excl + n_rand
    for (name, nh) in ((:coloc, n_head_coloc), (:exclusion, n_head_excl))
        nh >= min_per_head || error("""
            run_three_way_gate: the $name head's RESTRICTED evaluation set holds $nh items,
            below the pre-registered P13_MIN_EVAL_PER_HEAD = $min_per_head.

            ECE IS BIASED DOWNWARD AT SMALL n AND WHEN BINS ARE SPARSE, so a green ECE on a small
            evaluation set is not evidence of calibration -- it is an artifact of having too few
            samples per bin. The minimum was pre-registered precisely so this cannot be waved
            through. This is a diagnosis, not a bar to lower.""")
    end

    # --- 5. SCORE EVERY PAIR THROUGH THE DOCUMENTED READ SURFACE ---------------------------
    verbose && println("[2/6] scoring $m pairs through three_way_log_bf / three_way_probs …")
    lbf_c   = Vector{Float64}(undef, m)   # CORRECTED log BF(coloc : random)
    lbf_e   = Vector{Float64}(undef, m)   # CORRECTED log BF(exclusion : random)
    prob_c  = Vector{Float64}(undef, m)   # UNCORRECTED head probability, coloc head
    prob_e  = Vector{Float64}(undef, m)   # UNCORRECTED head probability, exclusion head
    lambdas = Vector{Float64}(undef, m)
    rho_s   = Vector{Float64}(undef, m)
    rho_c   = Vector{Float64}(undef, m)
    for j in 1:m
        it       = items[j]
        Zj       = p13_encode_pair(it.Zs, it.Zc, it.lambda)
        b        = three_way_log_bf(net, Zj, hlo)
        p        = three_way_probs(net, Zj)
        lbf_c[j] = b.coloc
        lbf_e[j] = b.exclusion
        prob_c[j] = p.coloc
        prob_e[j] = p.exclusion
        lambdas[j] = Float64(it.lambda)
        rho_s[j]   = Float64(it.rho_s)
        rho_c[j]   = Float64(it.rho_c)
    end

    is_c = classes .== COLOC
    is_e = classes .== EXCLUSION
    is_r = classes .== RANDOM

    # --- 6. DISCRIMINATION -- THIS IS THE GATE (D-12) ---------------------------------------
    # Per-class one-vs-random AUC at SIMULATOR GROUND-TRUTH labels, each head scored on its OWN
    # restricted set. Threshold-free by construction: an AUC is a ranking statistic, so no
    # decision rule and no operating point enters the gated number.
    fpr_c, tpr_c, auc_coloc     = roc_auc(lbf_c[is_r], lbf_c[is_c])
    fpr_e, tpr_e, auc_exclusion = roc_auc(lbf_e[is_r], lbf_e[is_e])
    auc_c_pass = auc_coloc     >= P13_AUC_FLOOR_COLOC
    auc_e_pass = auc_exclusion >= P13_AUC_FLOOR_EXCLUSION

    if verbose
        println("-"^78)
        println("DISCRIMINATION (D-12) -- THE GATE. One-vs-random AUC at simulator ground truth.")
        println("  coloc     vs random   AUC = $(round(auc_coloc; sigdigits = 6))   " *
                "(>= P13_AUC_FLOOR_COLOC=$P13_AUC_FLOOR_COLOC ? " *
                "$(auc_c_pass ? "PASS" : "FAIL"))   n = $n_col vs $n_rand")
        println("  exclusion vs random   AUC = $(round(auc_exclusion; sigdigits = 6))   " *
                "(>= P13_AUC_FLOOR_EXCLUSION=$P13_AUC_FLOOR_EXCLUSION ? " *
                "$(auc_e_pass ? "PASS" : "FAIL"))   n = $n_excl vs $n_rand")
    end

    # The exclusion AUC resolved over |rho_sample| bands. REPORTED, and it is the number the
    # D-04 iteration trigger is evaluated against -- the trigger is worded specifically about
    # the deep-exclusion tail near the -0.99 atom, so it cannot be evaluated from the pooled
    # AUC alone. The bands carry NO bar of their own; the frozen floor is applied unchanged.
    band_lo   = Float64[]
    band_hi   = Float64[]
    band_n    = Int[]
    band_auc  = Float64[]
    lo = Float64(P13_TAU)
    for hi in P13_GATE_RHO_BANDS
        sel = is_e .& (abs.(rho_s) .>= lo) .& (abs.(rho_s) .<= hi)
        push!(band_lo, lo); push!(band_hi, Float64(hi)); push!(band_n, count(sel))
        push!(band_auc, count(sel) > 0 ? p13_gate_auc(lbf_e[is_r], lbf_e[sel]) : NaN)
        lo = Float64(hi)
    end
    # The deep tail is the LAST band; the mid-range is everything below it.
    deep_sel     = is_e .& (abs.(rho_s) .> Float64(P13_GATE_RHO_BANDS[end - 1]))
    mid_sel      = is_e .& (abs.(rho_s) .<= Float64(P13_GATE_RHO_BANDS[end - 1]))
    auc_deep     = count(deep_sel) > 0 ? p13_gate_auc(lbf_e[is_r], lbf_e[deep_sel]) : NaN
    auc_mid      = count(mid_sel)  > 0 ? p13_gate_auc(lbf_e[is_r], lbf_e[mid_sel])  : NaN
    # MECHANICAL evaluation of the ONE pre-declared trigger (P13_ITERATION_TRIGGER): the
    # exclusion AUC is below its floor SPECIFICALLY AT HIGH |rho| WHILE THE MID-RANGE IS
    # RESOLVED. Both halves must hold; either alone does not fire it.
    trigger_fired = isfinite(auc_deep) && isfinite(auc_mid) &&
                    (auc_deep < P13_AUC_FLOOR_EXCLUSION) &&
                    (auc_mid  >= P13_AUC_FLOOR_EXCLUSION)

    if verbose
        println("  exclusion AUC resolved over |rho_sample| bands (REPORTING BANDS, no bar of")
        println("  their own -- the frozen floor is applied unchanged inside each):")
        for i in eachindex(band_lo)
            println("    |rho| in ($(band_lo[i]), $(band_hi[i])]   n = $(band_n[i])   " *
                    "AUC = $(round(band_auc[i]; sigdigits = 6))")
        end
        println("  deep tail |rho| > $(P13_GATE_RHO_BANDS[end - 1]):  AUC = " *
                "$(round(auc_deep; sigdigits = 6))   n = $(count(deep_sel))")
        println("  mid range |rho| <= $(P13_GATE_RHO_BANDS[end - 1]): AUC = " *
                "$(round(auc_mid; sigdigits = 6))   n = $(count(mid_sel))")
        println("  D-04 ITERATION TRIGGER (deep tail below floor AND mid-range resolved): " *
                "$(trigger_fired ? "FIRED" : "NOT FIRED")")
    end

    # --- 7. THE CONFUSION MATRIX -- DESCRIPTIVE, NOT A DECISION RULE (D-12) -----------------
    conf = zeros(Int, 3, 3)      # rows = TRUE class, cols = PREDICTED, both in class order
    class_index(cl) = cl === EXCLUSION ? 1 : (cl === RANDOM ? 2 : 3)
    for j in 1:m
        # argmax over (log BF(C:R), 0, log BF(E:R)) -- the random entry is the STRUCTURAL zero.
        trip = (lbf_c[j], 0.0, lbf_e[j])
        k    = argmax(trip)
        predicted = k == 1 ? 3 : (k == 2 ? 2 : 1)     # (coloc, random, exclusion) -> class order
        conf[class_index(classes[j]), predicted] += 1
    end

    if verbose
        println("-"^78)
        println("CONFUSION MATRIX AT argmax -- DESCRIPTIVE, NOT A DECISION RULE.")
        println("  The argmax can only declare `random` when BOTH log Bayes factors are")
        println("  negative, which is a DECISION rule; decisions and abstention are Phase 14's")
        println("  scope, and a descriptive argmax quoted as a classifier would pre-empt that")
        println("  phase's design. Nothing below carries an assertion.")
        println("  rows = TRUE class, cols = PREDICTED     " *
                "$(P13_GATE_CLASS_ORDER)")
        for i in 1:3
            println("    $(rpad(P13_GATE_CLASS_ORDER[i], 10)) " *
                    join([lpad(conf[i, j], 8) for j in 1:3], " "))
        end
    end

    # --- 8. CALIBRATION -- THE OTHER HALF OF THE GATE (D-13) --------------------------------
    # Per head, on that head's OWN restricted class pair, on the UNCORRECTED probabilities
    # against the binary label. Calibration is a property of the classifier UNDER ITS OWN
    # TRAINING FREQUENCIES, so subtracting the evidence-scale correction first would score a
    # quantity that is not a probability of anything (three_way_probs' docstring).
    sel_head_c = is_c .| is_r
    sel_head_e = is_e .| is_r
    cal_c = _bin_calibration(prob_c[sel_head_c], Vector{Bool}(is_c[sel_head_c]);
                             n_bins = P13_ECE_NBINS)
    cal_e = _bin_calibration(prob_e[sel_head_e], Vector{Bool}(is_e[sel_head_e]);
                             n_bins = P13_ECE_NBINS)

    # The head's OWN discrimination number, computed on the very quantity the ECE was computed
    # on. It coincides with the gated AUC above because subtracting a constant is rank-preserving
    # and an AUC is a ranking statistic -- both are reported so that identity is visible rather
    # than assumed.
    head_auc_c = p13_gate_auc(prob_c[is_r], prob_c[is_c])
    head_auc_e = p13_gate_auc(prob_e[is_r], prob_e[is_e])

    meta_c = p13_calibration_meta(cal_c; grid = G, auc = head_auc_c)
    meta_e = p13_calibration_meta(cal_e; grid = G, auc = head_auc_e)
    empty_c = meta_c.gate.empty_bins
    empty_e = meta_e.gate.empty_bins
    vac_c   = meta_c.gate.vacuous
    vac_e   = meta_e.gate.vacuous
    ece_c_pass = cal_c.ece <= P13_ECE_GREEN
    ece_e_pass = cal_e.ece <= P13_ECE_GREEN

    if verbose
        println("-"^78)
        println("CALIBRATION (D-13) -- ECE GATES, MCE DOES NOT. Every ECE carries its head AUC")
        println("and its empty-bin count, because a green ECE on a head that learned nothing is")
        println("a VACUOUS pass, and because _bin_calibration scores an EMPTY bin at ZERO ECE")
        println("weight but FULL MCE weight -- so MCE is dominated by empty bins on a")
        println("well-separated three-way problem and is not a usable gate statistic here.")
        for (nm, cal, mt, ab, ha, ep, vc) in
                ((:coloc, cal_c, meta_c, ece_c_pass, head_auc_c, empty_c, vac_c),
                 (:exclusion, cal_e, meta_e, ece_e_pass, head_auc_e, empty_e, vac_e))
            println("  $(rpad(nm, 10)) ECE = $(round(cal.ece; sigdigits = 6)) " *
                    "[$(mt.gate.verdict)] (<= P13_ECE_GREEN=$P13_ECE_GREEN ? " *
                    "$(ab ? "PASS" : "FAIL"))   " *
                    "MCE = $(round(cal.mce; sigdigits = 6)) (REPORTED, NEVER GATED; " *
                    "empty bins = $ep of $P13_ECE_NBINS)   " *
                    "head AUC = $(round(ha; sigdigits = 6))" *
                    (vc ? "   *** VACUOUS PASS ***" : ""))
        end
        if vac_c || vac_e
            println("  *** VACUOUS: a green ECE on a head whose AUC is at or below")
            println("      P13_VACUOUS_AUC_FLOOR = $P13_VACUOUS_AUC_FLOOR is NOT calibration")
            println("      evidence and must never be quoted as one. A head that says")
            println("      \"one half, always\" is perfectly calibrated and perfectly useless.")
        end
        println("  CONTEXT, NOT A THRESHOLD: the shipped BINARY net's decision-calibration ECE")
        println("  in spike 014 was 0.0188. It is quoted beside the coloc head's ECE " *
                "($(round(cal_c.ece; sigdigits = 6))) as a")
        println("  scale reference only -- it is a different net on a different input surface")
        println("  and no Phase-13 verdict turns on it.")
    end

    # --- 9. THE LAMBDA RESPONSE (D-03) -- REPORTED --------------------------------------------
    # The SAME summaries at OTHERWISE MATCHED inputs, re-encoded at each rung of Phase 11's own
    # frozen SC2 ladder. Only the appended conditioning value moves.
    verbose && println("[3/6] measuring the lambda response over the Phase-11 SC2 ladder …")
    rungs      = collect(Float64, SC2_RUNGS)
    mean_abs_c = Float64[]; mean_abs_e = Float64[]
    mean_c_on_c = Float64[]; mean_e_on_e = Float64[]
    mean_c_on_r = Float64[]; mean_e_on_r = Float64[]
    for lam in rungs
        bc = Vector{Float64}(undef, m); be = Vector{Float64}(undef, m)
        for j in 1:m
            Zj    = p13_encode_pair(items[j].Zs, items[j].Zc, lam)
            b     = three_way_log_bf(net, Zj, hlo)
            bc[j] = b.coloc; be[j] = b.exclusion
        end
        push!(mean_abs_c, mean(abs.(bc)));      push!(mean_abs_e, mean(abs.(be)))
        push!(mean_c_on_c, mean(bc[is_c]));     push!(mean_e_on_e, mean(be[is_e]))
        push!(mean_c_on_r, mean(bc[is_r]));     push!(mean_e_on_r, mean(be[is_r]))
    end
    span_c = maximum(mean_abs_c) - minimum(mean_abs_c)
    span_e = maximum(mean_abs_e) - minimum(mean_abs_e)
    spear_c = StatsBase.corspearman(rungs, mean_abs_c)
    spear_e = StatsBase.corspearman(rungs, mean_abs_e)

    if verbose
        println("-"^78)
        println("LAMBDA RESPONSE (D-03) -- REPORTED, NOT GATED. Identical summaries, only the")
        println("appended registration-uncertainty conditioning value moves.")
        println("  lambda   mean|logBF_C|   mean|logBF_E|   mean logBF_C|coloc   " *
                "mean logBF_E|exclusion")
        for i in eachindex(rungs)
            println("  $(rpad(rungs[i], 7)) $(lpad(round(mean_abs_c[i]; sigdigits = 5), 14)) " *
                    "$(lpad(round(mean_abs_e[i]; sigdigits = 5), 15)) " *
                    "$(lpad(round(mean_c_on_c[i]; sigdigits = 5), 21)) " *
                    "$(lpad(round(mean_e_on_e[i]; sigdigits = 5), 22))")
        end
        println("  span of mean|logBF| across the ladder: coloc head " *
                "$(round(span_c; sigdigits = 5)), exclusion head $(round(span_e; sigdigits = 5))")
        println("  Spearman(rung, mean|logBF|): coloc $(round(spear_c; sigdigits = 4)), " *
                "exclusion $(round(spear_e; sigdigits = 4))")
        println("  A FLAT RESPONSE IS AN HONEST FINDING ABOUT THE CONDITIONING INPUT, NOT A BUG")
        println("  AND NOT SOMETHING TO EXPLAIN AWAY. Phase 11 established that registration")
        println("  uncertainty at this scale is not recoverable from the fixed 8x8 summary and")
        println("  that the posterior width is essentially flat in lambda; a flat evidence")
        println("  response here is the same finding one level up. It is recorded either way.")
    end

    # --- 10. BINARY-NRE CONTINUITY -- REPORTED, NEVER GATED ---------------------------------
    verbose && println("[4/6] continuity against the frozen binary NRE (reported, not gated) …")
    cont_n = 0; cont_corr = NaN; cont_maxabs = NaN; cont_mean = NaN
    cont_three = Float64[]; cont_bin = Float64[]; cont_note = ""
    cont_lambda = Float64(P13_PHASE11_REFERENCE_LAMBDA)
    if isfile(binary_net_path)
        rat = load_ratio(binary_net_path)
        idx = findall(sel_head_c)                    # the overlapping two-way regime: C and R
        cont_three = Vector{Float64}(undef, length(idx))
        cont_bin   = Vector{Float64}(undef, length(idx))
        for (k, j) in enumerate(idx)
            Zj  = p13_encode_pair(items[j].Zs, items[j].Zc, cont_lambda)
            cont_three[k] = three_way_log_bf(net, Zj, hlo).coloc
            Zs_b = p13_to_binary_basis(collect(items[j].Zs), basis, rat.zt)
            Zc_b = p13_to_binary_basis(collect(items[j].Zc), basis, rat.zt)
            cont_bin[k] = amortized_log_bf(rat.estimator, pair_encode(Zs_b, Zc_b),
                                           rat.log_prior_odds)
        end
        cont_n       = length(idx)
        cont_corr    = cor(cont_three, cont_bin)
        cont_maxabs  = maximum(abs.(cont_three .- cont_bin))
        cont_mean    = mean(cont_three .- cont_bin)
        cont_note    = "spike-lane frozen binary NRE (Phase-5 RatioEstimator on the frozen " *
                       "Phase-4 basis); the shipped grid-8 bundle is reserved by consts.jl " *
                       "section I2 as an OOD comparison reference only"
    else
        cont_note = "NOT COMPUTED: $binary_net_path is absent (gitignored, regenerable). " *
                    "Continuity is reported-not-gated, so no verdict depends on it."
    end

    if verbose
        println("-"^78)
        println("BINARY-NRE CONTINUITY -- REPORTED, NEVER GATED (P13_CONTINUITY_GATED = " *
                "$P13_CONTINUITY_GATED).")
        println("  TWO ARCHITECTURAL REASONS, BOTH DECIDABLE BEFORE A SINGLE WEIGHT EXISTED:")
        println("   (1) DIFFERENT INPUT SURFACES (D-02). This net reads Phase-11's")
        println("       registration-aware frozen basis plus an appended lambda conditioning")
        println("       row; the binary net reads a different frozen basis and has no")
        println("       conditioning input at all. A disagreement is then a statement about the")
        println("       basis change and the conditioning input, not about either net's")
        println("       evidence.")
        println("   (2) THE TWO LOGITS ARE DIFFERENT OBJECTS. The binary net's logit is a")
        println("       shuffled-theta likelihood-to-evidence ratio, so its model-index")
        println("       difference IS the log Bayes factor with no correction term. The D-09")
        println("       heads are plain BCE classifiers on class labels, whose Bayes-optimal")
        println("       logit carries log(q_A/q_B) and therefore NEEDS a correction the binary")
        println("       net does not. They are not numerically interchangeable even at")
        println("       identical inputs.")
        println("  read at lambda = $cont_lambda over $cont_n coloc/random pairs")
        println("  corr(three-way logBF(C:R), binary logBF) = $(round(cont_corr; sigdigits = 5))")
        println("  max|difference|                          = " *
                "$(round(cont_maxabs; sigdigits = 5))")
        println("  mean difference                          = $(round(cont_mean; sigdigits = 5))")
        println("  source: $cont_note")
        println("  NO ASSERTION IS ATTACHED TO ANY NUMBER IN THIS SECTION.")
    end

    # --- 11. PERSIST THE ARTIFACT BEFORE ANY VERDICT (T-13-03, D-04) ------------------------
    verbose && println("[5/6] persisting the reported artifact BEFORE the gate …")
    consts_path = joinpath(@__DIR__, "consts.jl")
    lambda_response = (lambda = rungs,
                       mean_abs_logbf_coloc = mean_abs_c,
                       mean_abs_logbf_exclusion = mean_abs_e,
                       mean_logbf_coloc_on_coloc = mean_c_on_c,
                       mean_logbf_exclusion_on_exclusion = mean_e_on_e,
                       mean_logbf_coloc_on_random = mean_c_on_r,
                       mean_logbf_exclusion_on_random = mean_e_on_r,
                       span_coloc = span_c, span_exclusion = span_e,
                       spearman_coloc = spear_c, spearman_exclusion = spear_e,
                       ladder_source = "SC2_RUNGS (spike/validation/p11_consts.jl), READ",
                       gated = false)
    continuity = (n = cont_n, lambda = cont_lambda, corr = cont_corr,
                  max_abs_diff = cont_maxabs, mean_diff = cont_mean,
                  logbf_three_way = cont_three, logbf_binary = cont_bin,
                  binary_net = binary_net_path, note = cont_note,
                  gated = P13_CONTINUITY_GATED)

    _p13_gate_save_report(report_path;
        # --- the GATED statistics ---
        auc_coloc = auc_coloc, auc_exclusion = auc_exclusion,
        ece_coloc = cal_c.ece, ece_exclusion = cal_e.ece,
        vacuous_coloc = vac_c, vacuous_exclusion = vac_e,
        # --- reported beside every gated number ---
        mce_coloc = cal_c.mce, mce_exclusion = cal_e.mce,
        empty_bins_coloc = empty_c, empty_bins_exclusion = empty_e,
        head_auc_coloc = head_auc_c, head_auc_exclusion = head_auc_e,
        verdict_coloc = meta_c.gate.verdict, verdict_exclusion = meta_e.gate.verdict,
        confusion = conf, confusion_class_order = collect(String.(P13_GATE_CLASS_ORDER)),
        lambda_response = lambda_response, continuity = continuity,
        # --- the reliability curves themselves ---
        bin_midpoints_coloc = cal_c.bin_midpoints, predicted_rate_coloc = cal_c.predicted_rate,
        observed_rate_coloc = cal_c.observed_rate, bin_counts_coloc = cal_c.bin_counts,
        bin_midpoints_exclusion = cal_e.bin_midpoints,
        predicted_rate_exclusion = cal_e.predicted_rate,
        observed_rate_exclusion = cal_e.observed_rate, bin_counts_exclusion = cal_e.bin_counts,
        roc_fpr_coloc = fpr_c, roc_tpr_coloc = tpr_c,
        roc_fpr_exclusion = fpr_e, roc_tpr_exclusion = tpr_e,
        # --- the exclusion AUC resolved over |rho| and the MECHANICAL D-04 trigger ---
        band_lo = band_lo, band_hi = band_hi, band_n = band_n, band_auc = band_auc,
        auc_exclusion_deep_tail = auc_deep, auc_exclusion_mid_range = auc_mid,
        n_exclusion_deep_tail = count(deep_sel), n_exclusion_mid_range = count(mid_sel),
        iteration_trigger_fired = trigger_fired,
        iteration_trigger_text = P13_ITERATION_TRIGGER,
        # --- the per-item reported scores, so a later plan need not re-run the gate ---
        logbf_coloc = lbf_c, logbf_exclusion = lbf_e,
        prob_coloc = prob_c, prob_exclusion = prob_e,
        class = Int.(classes), rho_sample = rho_s, rho_control = rho_c, lambda = lambdas,
        # --- the realized evaluation set ---
        m = m, n_exclusion = n_excl, n_random = n_rand, n_coloc = n_col,
        n_head_coloc = n_head_coloc, n_head_exclusion = n_head_excl,
        realized_class_freq = realized, n_draws = strat.n_draws,
        # --- EVERY THRESHOLD THIS RUN WAS SCORED AGAINST, written INTO the artifact ---
        P13_GATE_M = P13_GATE_M,
        P13_AUC_FLOOR_COLOC = P13_AUC_FLOOR_COLOC,
        P13_AUC_FLOOR_EXCLUSION = P13_AUC_FLOOR_EXCLUSION,
        P13_MIN_EVAL_PER_HEAD = P13_MIN_EVAL_PER_HEAD,
        P13_ECE_GREEN = P13_ECE_GREEN, P13_ECE_YELLOW = P13_ECE_YELLOW,
        P13_ECE_NBINS = P13_ECE_NBINS,
        P13_VACUOUS_AUC_FLOOR = P13_VACUOUS_AUC_FLOOR,
        P13_GATE_STATISTIC = P13_GATE_STATISTIC,
        P13_CONFUSION_RULE = P13_CONFUSION_RULE,
        P13_CONTINUITY_GATED = P13_CONTINUITY_GATED,
        P13_TOST_REQUIRED = P13_TOST_REQUIRED,
        P13_TARGET_CLASS_FREQ = P13_TARGET_CLASS_FREQ,
        P13_ITERATION_ALLOWANCE = P13_ITERATION_ALLOWANCE,
        P13_TAU = P13_TAU, P13_TAU_REFERENCE_LAMBDA = P13_TAU_REFERENCE_LAMBDA,
        P13_CUT_VARIANT = P13_CUT_VARIANT,
        rho_bands = collect(Float64, P13_GATE_RHO_BANDS),
        # --- provenance: the stream, the net, the frozen pre-registration ---
        master_seed = UInt64(P13_DEV_SEED), salt = UInt64(P13_SALT),
        gate_counter = Int(P13_GATE_COUNTER), datagen_counter = Int(P13_DATAGEN_COUNTER),
        net_path = net_path, net_head_log_odds = hlo,
        net_consts_sha = h.consts_sha, net_meta = h.meta,
        consts_sha = p13_consts_sha(), consts_git_blob_sha = p13_blob_sha(consts_path),
        phase11_net = P13_PHASE11_NET, phase11_sha = p13_phase11_sha(),
        grid = G, reported = reported,
        amendment = "13-SC2-AMENDMENT.md -- any report citing this result must cite the " *
                    "original ROADMAP SC1/SC2/SC3 alongside it",
        generated = string(Dates.now(Dates.UTC)) * "Z")
    verbose && println("      persisted (before any verdict) -> $report_path")

    # --- 12. THE FIGURE. A SCRIPT ARTIFACT, NEVER A GATE ASSERTION --------------------------
    # RENDERED BEFORE THE GATE, following run_bf.jl:118-120, for the same reason the artifact is:
    # a failing @testset throws and everything after it is never reached, so a figure emitted
    # after the assertions would exist only for a PASSING run. The evidence of an honest failure
    # must survive the failure.
    verbose && println("[6/6] rendering the confusion / ROC figure …")
    figout = p13_gate_figure(conf, fpr_c, tpr_c, auc_coloc, fpr_e, tpr_e, auc_exclusion;
                             path = fig_path)
    verbose && println("      figure -> $figout")

    # --- 13. THE HEADLINE, PRINTED BEFORE THE GATE CAN THROW --------------------------------
    overall = auc_c_pass && auc_e_pass && ece_c_pass && ece_e_pass && !vac_c && !vac_e
    if verbose
        println("\n", "-"^78)
        println("Overall reported three-way gate: ", overall ? "PASS" : "FAIL")
        println("  AUC coloc     $(round(auc_coloc; sigdigits = 6)) " *
                "$(auc_c_pass ? ">=" : "<") $P13_AUC_FLOOR_COLOC")
        println("  AUC exclusion $(round(auc_exclusion; sigdigits = 6)) " *
                "$(auc_e_pass ? ">=" : "<") $P13_AUC_FLOOR_EXCLUSION")
        println("  ECE coloc     $(round(cal_c.ece; sigdigits = 6)) " *
                "$(ece_c_pass ? "<=" : ">") $P13_ECE_GREEN")
        println("  ECE exclusion $(round(cal_e.ece; sigdigits = 6)) " *
                "$(ece_e_pass ? "<=" : ">") $P13_ECE_GREEN")
        println("  vacuous coloc $vac_c   vacuous exclusion $vac_e")
        println("  D-04 iteration trigger: $(trigger_fired ? "FIRED" : "NOT FIRED")  " *
                "(allowance $P13_ITERATION_ALLOWANCE, evaluated mechanically)")
        println("-"^78, "\n")
    end

    result = (auc_coloc = auc_coloc, auc_exclusion = auc_exclusion,
              ece_coloc = cal_c.ece, ece_exclusion = cal_e.ece,
              mce_coloc = cal_c.mce, mce_exclusion = cal_e.mce,
              empty_bins_coloc = empty_c, empty_bins_exclusion = empty_e,
              head_auc_coloc = head_auc_c, head_auc_exclusion = head_auc_e,
              vacuous_coloc = vac_c, vacuous_exclusion = vac_e,
              confusion = conf, lambda_response = lambda_response, continuity = continuity,
              auc_exclusion_deep_tail = auc_deep, auc_exclusion_mid_range = auc_mid,
              iteration_trigger_fired = trigger_fired,
              overall = overall, report_path = report_path, figure = figout,
              reported = reported)

    # --- 14. THE REAL GATE: exits nonzero if any pre-registered criterion fails --------------
    # SIX assertions, and only six. NO MCE assertion (the empty-bin trap: an empty bin carries
    # ZERO ece weight and FULL mce weight, so it would fail for a reason unrelated to
    # calibration). NO continuity assertion (P13_CONTINUITY_GATED is false and the two nets do
    # not share an input surface). NO confusion-matrix assertion (it is descriptive). Every
    # threshold is REFERENCED from the frozen pre-registration, never inlined.
    if reported
        @testset "Reported three-way pre-registered gate (D-12/D-13)" begin
            @test auc_coloc     >= P13_AUC_FLOOR_COLOC
            @test auc_exclusion >= P13_AUC_FLOOR_EXCLUSION
            @test cal_c.ece     <= P13_ECE_GREEN
            @test cal_e.ece     <= P13_ECE_GREEN
            @test !vac_c
            @test !vac_e
        end
    else
        println("SMOKE MODE: the gate assertions were NOT run and this output is NOT a verdict.")
    end
    return result
end

"""
    p13_gate_figure(conf, fpr_c, tpr_c, auc_c, fpr_e, tpr_e, auc_e; path) -> String

The confusion-matrix heatmap and the two one-vs-random ROC curves, in one PNG.

A SCRIPT ARTIFACT AND NEVER A GATE ASSERTION (the Phase-4 rule `spike/validation/figures.jl`
states at its head, and which this file honours by rendering before the gate rather than after).
Composed here rather than through a `figures.jl` helper for two reasons: that file has no
confusion-matrix helper, and its `_val_fig_path` hard-routes every output into
`spike/validation/figures/`, which is not where this deliverable lives. The CairoMakie surface,
the headless backend activation and the caption discipline are all inherited from it.
"""
function p13_gate_figure(conf, fpr_c, tpr_c, auc_c, fpr_e, tpr_e, auc_e;
                         path = P13_GATE_FIG_PATH)
    mkpath(dirname(path))
    labels = String.(collect(P13_GATE_CLASS_ORDER))
    fig = Figure(size = (1000, 460))
    ax1 = CairoMakie.Axis(fig[1, 1];
        title = "confusion at argmax -- DESCRIPTIVE, NOT A DECISION RULE",
        subtitle = "rows = simulator ground truth, cols = argmax(logBF_C, 0, logBF_E)",
        xlabel = "predicted", ylabel = "true",
        xticks = (1:3, labels), yticks = (1:3, labels))
    heatmap!(ax1, 1:3, 1:3, permutedims(Float64.(conf)); colormap = :Blues)
    for i in 1:3, j in 1:3
        text!(ax1, j, i; text = string(conf[i, j]), align = (:center, :center),
              color = :black)
    end
    ax2 = CairoMakie.Axis(fig[1, 2];
        title = "one-vs-random ROC (D-12 gate statistic)",
        subtitle = "coloc AUC=$(round(auc_c; digits = 4))  " *
                   "exclusion AUC=$(round(auc_e; digits = 4))",
        xlabel = "false positive rate", ylabel = "true positive rate")
    lines!(ax2, [0, 1], [0, 1]; color = :gray, linestyle = :dash)
    lines!(ax2, fpr_c, tpr_c; color = :steelblue, label = "coloc vs random")
    lines!(ax2, fpr_e, tpr_e; color = :firebrick, label = "exclusion vs random")
    axislegend(ax2; position = :rb)
    save(path, fig)
    return path
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so `include`-ing this file can NEVER trigger the reported gate. Run it deliberately:
#
#     julia --project=spike -t auto spike/p13/run_three_way_gate.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
