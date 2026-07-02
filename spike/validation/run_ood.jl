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

# spike/validation/run_ood.jl --- Phase-5 REPORTED OOD ROC (OOD-01, OOD-02).
#
# THE ACTUAL OOD DETECTION PROOF. Runs the controlled ROC experiment over the FULL
# 4-family × OOD_GRID_LEVELS misspecification grid on the FRESH RESERVED stream
# val_rng(VAL_MASTER_SEED), at the PRE-REGISTERED ID-quantile operating point
# (OOD_ID_QUANTILE), and gates it: the summary-density (Mahalanobis) AUC — the OR-flag's
# primary channel — (combined AND each family's strongest-level) ≥ OOD_AUC_MIN AND each
# summary-orthogonal negative control KS-invariant (< OOD_KS_EPS) with the flag QUIET. It
# is a REAL GATE: it EXITS NONZERO on any pre-registered failure.
#
# TRAIN-ONLY FREEZE (SC5 / D-06): the OOD nulls are fit on FRESH IN-DISTRIBUTION draws
# through the FROZEN read chain (id_summary → the training distribution in the frozen
# m.zt space) — NEVER on eval/misspec data, which is fully external and only ever SCORED.
# The Mahalanobis fit touches the continuous rows only (Pitfall 2); the ID-quantile
# operating point is committed on a HELD-OUT ID pool (out-of-sample, honest ~5% FPR).
#
# POSTERIOR-PREDICTIVE (PP) CHANNEL — REPORTED NON-VIABLE ON THIS FROZEN NET (honest
# Rule-1 finding, 05-04): on strongly-OOD inputs the frozen Phase-4 NPE returns a
# NON-FINITE posterior mean θ̂ (the flow extrapolates to NaN out of distribution — the same
# overconfidence/miscalibration the reported SBC gate independently exposed), so the PP
# re-simulation throws in simulate_pair. The PP channel is therefore excluded from the
# reported OR-fusion (with_pp=false); the density channel is the reported OOD detector and
# the PP non-viability is EVIDENCED (a finiteness probe) and reported, NOT hidden. The fix
# belongs in ood.jl's pp_mismatch_score (guard non-finite θ̂), which this plan is not
# permitted to edit (shared-file isolation) — surfaced to Phase-6 / gap-closure.
#
# NAMED BLIND SPOT (D-04, honesty capstone): the fixed 8×8 patch-Pearson summary is a
# CORRELATION statistic, so correlation-preserving transforms (affine / rotate-flip /
# block-permute) are provably KS-invariant on the summary and the flag stays quiet on them
# BY CONSTRUCTION — the structural blind spot MEASURED and NAMED, paired with the SBC
# "under the simulator" caveat (SBC-04). A per-family AUC shortfall (e.g. the documented
# detector-noise blind spot, 05-03) is an HONEST finding for Phase-6, never tuned away by
# editing consts.jl or reseeding.
#
# CPU-only (D-10), read-only src/ (via the harness). Run:
#     julia --project=spike spike/validation/run_ood.jl

using Test
using JLD2
using Statistics
using Dates

# Units under test: the OOD surface (pulls ood.jl → harness.jl → consts.jl transitively)
# and the figure surface (CairoMakie — figures are script artifacts, not gate assertions).
isdefined(@__MODULE__, :ood_roc_over_grid) || include(joinpath(@__DIR__, "ood.jl"))
isdefined(@__MODULE__, :plot_roc)          || include(joinpath(@__DIR__, "figures.jl"))

const OOD_REPORT_PATH = joinpath(@__DIR__, "ood_report.jld2")

# Monte-Carlo resolution knobs for the reported grid (NOT anti-snooping thresholds —
# the pass/fail constants OOD_AUC_MIN/OOD_KS_EPS/OOD_ID_QUANTILE/OOD_GRID_LEVELS are
# locked in consts.jl). Sized for a stable reported AUC.
const OOD_N_ID  = 100     # fresh ID negatives (null-fit + held-out operating point + ROC)
const OOD_N_POS = 40      # external positives per family × level
const OOD_NEG_REPS = 30   # negative-control replicates (KS-invariance + flag-quiet rate)

# --- Atomic persistence of the reported OOD artifact (save_npe idiom, T-05-13) ---------
function _save_ood_report(path; kwargs...)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "combined_auc") "_save_ood_report: integrity check failed ($tmp)"
    end
    mv(tmp, path; force = true)
    return path
end

# --- PP-viability probe: EVIDENCE the reported PP-channel non-viability (honest finding).
# For each family at the strongest level, infer the NPE posterior draws on ONE external
# misspec summary and record whether θ̂ is finite. A non-finite θ̂ is exactly what makes
# pp_mismatch_score throw in simulate_pair (the reason PP is excluded from the reported
# OR-fusion). Cheap (one posterior pass per family; no PP re-simulation).
function _pp_viability_probe(m, rng)
    viable = Dict{Symbol,Bool}()
    for (fam, gen) in pairs(OOD_FAMILIES)
        θ  = sample_prior(rng)
        Z  = frozen_summary(m, gen(rng, θ; imsize = SBC_IMSIZE, level = OOD_GRID_LEVELS))
        draws = posterior_for(m.estimator, Z; N = 64, use_gpu = false)
        θhat  = vec(mean(StatsBase.reconstruct(m.θzt, draws); dims = 2))
        viable[fam] = all(isfinite, θhat)
    end
    return viable
end

function main()
    println("="^78)
    println("Phase-5 REPORTED OOD ROC (OOD-01/02) — LOCKED consts, RESERVED stream")
    println("  grid: $(length(OOD_FAMILIES)) families × $OOD_GRID_LEVELS levels  imsize=$SBC_IMSIZE")
    println("  n_id=$OOD_N_ID  n_pos=$OOD_N_POS  operating point q=$OOD_ID_QUANTILE")
    println("  gate: AUC ≥ OOD_AUC_MIN=$OOD_AUC_MIN  AND  neg-ctrl KS < OOD_KS_EPS=$OOD_KS_EPS (flag quiet)")
    println("  stream=VAL_MASTER_SEED=$(repr(VAL_MASTER_SEED))")
    println("  NOTE: PP channel excluded from OR-fusion (non-viable on this net — see banner).")
    println("="^78)
    @assert VAL_MASTER_SEED != NPE_MASTER_SEED "reported stream must be disjoint from training stream (D-02)"
    @assert VAL_MASTER_SEED != VAL_FIX_SEED    "reported stream must be disjoint from fixture stream (D-02)"

    m = load_frozen_model()

    # --- Nulls fit on TRAIN-DISTRIBUTION only (SC5/D-06): fresh in-distribution draws
    #     through the frozen chain. NEVER fit on eval/misspec data (which is external). ---
    println("\nfitting OOD nulls on fresh in-distribution (train-distribution) draws …")
    rng_fit = val_rng(VAL_MASTER_SEED)
    Ztrain = reduce(hcat, [id_summary(m, rng_fit; imsize = SBC_IMSIZE) for _ in 1:OOD_N_ID])
    nulls  = fit_ood_nulls(Ztrain)
    @assert nulls.cont == collect(1:64) "nulls must fit continuous rows only (Pitfall 2)"

    # --- PP-channel viability, evidenced (honest finding) -------------------------------
    println("probing PP-channel viability (posterior θ̂ finiteness on strong misspec) …")
    pp_viable = _pp_viability_probe(m, val_rng(VAL_MASTER_SEED))
    println("      PP θ̂ finite per family: ", pp_viable)

    # --- (1) The controlled ROC over the FULL 4-family grid, DENSITY channel (OOD-02).
    #     with_pp=false: the PP channel is non-viable on this net (see banner), so it is
    #     excluded from the reported OR-fusion; the summary-density channel is the reported
    #     OOD detector. n_id/n_pos are resolution knobs; the gate constant is OOD_AUC_MIN.
    println("[1/2] running the density-channel ROC over the full grid …")
    t0 = time()
    report = ood_roc_over_grid(m, nulls; families = OOD_FAMILIES, levels = OOD_GRID_LEVELS,
                               n_id = OOD_N_ID, n_pos = OOD_N_POS,
                               rng = val_rng(VAL_MASTER_SEED), imsize = SBC_IMSIZE,
                               with_pp = false)
    println("      grid ROC done in $(round((time()-t0)/60; digits=2)) min")

    fam_names    = collect(report.families)
    maha_thr     = report.maha_thr
    fam_best_auc = Dict(fam => maximum(report.maha_auc[fam]) for fam in fam_names)

    # --- (2) Summary-orthogonal NEGATIVE control (D-04): KS-invariance + density-flag
    #     quiet-rate at the pre-registered ID-quantile operating point (density channel).
    println("[2/2] evaluating the summary-orthogonal negative control (named blind spot) …")
    rng_neg = val_rng(VAL_MASTER_SEED)
    negctrls = (
        affine = (p -> negctrl_affine(p; a = 2.0, b = 0.5)),
        rotate = (p -> negctrl_rotate_flip(p; k = 1)),
        block  = (p -> negctrl_block_permute(val_rng(VAL_MASTER_SEED), p; blocks = 8)),
    )
    neg_ks        = Dict{Symbol,Vector{Float64}}()
    neg_fire_rate = Dict{Symbol,Float64}()
    for (name, tf) in pairs(negctrls)
        ksvals = Float64[]
        fires  = 0
        for _ in 1:OOD_NEG_REPS
            θ    = sample_prior(rng_neg)
            base = simulate_pair(rng_neg, θ; imsize = SBC_IMSIZE)
            push!(ksvals, verify_summary_invariance(base, tf))
            # Density-channel flag (PP excluded): fires iff maha exceeds the ID-quantile thr.
            fired = maha_score(nulls, frozen_summary(m, tf(base))) > maha_thr
            fires += fired ? 1 : 0
        end
        neg_ks[name]        = ksvals
        neg_fire_rate[name] = fires / OOD_NEG_REPS
    end

    # --- Persist the reported artifact BEFORE the gate (survives an honest fail) --------
    _save_ood_report(OOD_REPORT_PATH;
        families = String.(fam_names),
        levels = report.levels,
        maha_auc = Dict(String(k) => v for (k, v) in report.maha_auc),
        fam_best_auc = Dict(String(k) => v for (k, v) in fam_best_auc),
        fire_rate = Dict(String(k) => v for (k, v) in report.fire_rate),
        combined_auc = report.combined_auc,
        roc_fpr = report.roc_fpr, roc_tpr = report.roc_tpr,
        maha_thr = maha_thr,
        id_fire_rate = report.id_fire_rate,
        youden_j = report.youden.j, youden_fpr = report.youden.fpr, youden_tpr = report.youden.tpr,
        pp_channel_viable = Dict(String(k) => v for (k, v) in pp_viable),
        pp_channel_note = "PP excluded from OR-fusion: non-finite posterior θ̂ on strong misspec (NPE overconfidence)",
        neg_ks = Dict(String(k) => v for (k, v) in neg_ks),
        neg_ks_max = Dict(String(k) => maximum(v) for (k, v) in neg_ks),
        neg_fire_rate = Dict(String(k) => v for (k, v) in neg_fire_rate),
        OOD_AUC_MIN = OOD_AUC_MIN, OOD_KS_EPS = OOD_KS_EPS,
        OOD_ID_QUANTILE = OOD_ID_QUANTILE, OOD_GRID_LEVELS = OOD_GRID_LEVELS,
        OOD_PP_REPS = OOD_PP_REPS, OOD_IMSIZE = collect(SBC_IMSIZE),
        n_id = OOD_N_ID, n_pos = OOD_N_POS,
        VAL_MASTER_SEED = UInt64(VAL_MASTER_SEED),
        generated = string(Dates.now(Dates.UTC)) * "Z")
    println("persisted reported artifact -> $OOD_REPORT_PATH")

    # --- Figure: the density-channel pooled-strongest-level ROC (reported headline curve) ---
    plot_roc(report.roc_fpr, report.roc_tpr, report.combined_auc; filename = "ood_roc_combined.png")
    println("figure -> $(joinpath(@__DIR__, "figures"))/ood_roc_combined.png")

    # --- Reported tables (printed BEFORE the gate throws) -------------------------------
    println("\n", "-"^78)
    println("Per-family summary-density (Mahalanobis) AUC over the misspec grid:")
    println(rpad("family", 14), [rpad("L$l", 9) for l in report.levels]..., rpad("best", 9), "gate")
    fam_pass = Dict{Symbol,Bool}()
    for fam in fam_names
        pass = fam_best_auc[fam] >= OOD_AUC_MIN
        fam_pass[fam] = pass
        println(rpad(String(fam), 14),
                [rpad(round(report.maha_auc[fam][l]; sigdigits=3), 9) for l in report.levels]...,
                rpad(round(fam_best_auc[fam]; sigdigits=3), 9),
                pass ? "PASS" : "FAIL")
    end
    println("\ndensity pooled strongest-level AUC = $(round(report.combined_auc; sigdigits=4)) " *
            "(≥ OOD_AUC_MIN=$OOD_AUC_MIN ? ", report.combined_auc >= OOD_AUC_MIN ? "PASS" : "FAIL", ")")
    println("held-out ID fire-rate (should ≈ 1−q = $(round(1-OOD_ID_QUANTILE; digits=2))): " *
            "$(round(report.id_fire_rate; sigdigits=3))")
    println("post-hoc Youden-J (reference only, NOT the gate): j=$(round(report.youden.j; sigdigits=3))")
    println("PP-channel viability (θ̂ finite on strong misspec): $pp_viable  -> excluded from OR-fusion")

    println("\nSummary-orthogonal negative control (D-04 named blind spot):")
    println(rpad("transform", 14), rpad("max KS", 12), rpad("< OOD_KS_EPS", 14),
            rpad("flag-fire", 12), "quiet")
    neg_pass = Dict{Symbol,Bool}()
    for name in keys(negctrls)
        ksmax = maximum(neg_ks[name])
        fr    = neg_fire_rate[name]
        ks_ok = ksmax < OOD_KS_EPS
        quiet = fr <= (1 - OOD_ID_QUANTILE) + 1e-9   # ≤ ID false-positive rate
        neg_pass[name] = ks_ok && quiet
        println(rpad(String(name), 14), rpad(round(ksmax; sigdigits=3), 12),
                rpad(ks_ok ? "yes" : "NO", 14), rpad(round(fr; sigdigits=3), 12),
                quiet ? "yes" : "NO")
    end
    overall = all(values(fam_pass)) && (report.combined_auc >= OOD_AUC_MIN) && all(values(neg_pass))
    println("-"^78)
    println("NAMED BLIND SPOT: the fixed 8×8 patch-Pearson summary is a CORRELATION statistic;")
    println("  summary-orthogonal discrepancies it provably cannot see are measured above and")
    println("  paired with the SBC \"under the simulator\" caveat (SBC-04).")
    println("Overall reported OOD gate: ", overall ? "PASS" : "FAIL")
    println("-"^78, "\n")

    # --- THE REAL GATE (OOD-01/02): exits nonzero on any pre-registered failure.
    # Thresholds referenced from consts.jl (never tuned). An honest failure → Phase-6.
    @testset "Reported OOD pre-registered gate (OOD-01/02, D-03..06)" begin
        @testset "positive controls: AUC ≥ OOD_AUC_MIN" begin
            @test report.combined_auc >= OOD_AUC_MIN
            for fam in fam_names
                @testset "$(fam)" begin
                    @test fam_best_auc[fam] >= OOD_AUC_MIN
                end
            end
        end
        @testset "negative control: KS-invariant (< OOD_KS_EPS) + flag quiet" begin
            for name in keys(negctrls)
                @testset "$(name)" begin
                    @test maximum(neg_ks[name]) < OOD_KS_EPS
                    @test neg_fire_rate[name] <= (1 - OOD_ID_QUANTILE) + 1e-9
                end
            end
        end
    end
    return nothing
end

main()
