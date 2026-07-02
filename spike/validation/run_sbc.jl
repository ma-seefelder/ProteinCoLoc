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

# spike/validation/run_sbc.jl --- Phase-5 REPORTED SBC run (SBC-01..04).
#
# THE ACTUAL CALIBRATION PROOF. This script runs the Simulation-Based-Calibration
# experiment at the LOCKED, PRE-REGISTERED constants (SBC_M=2000, SBC_L=999,
# SBC_BINS=50; SBC_KS_ALPHA/SBC_CHI2_ALPHA/SBC_ECE_GREEN thresholds) on the FRESH,
# RESERVED, never-tuned-against stream val_rng(VAL_MASTER_SEED) — disjoint from the
# fixture stream (VAL_FIX_SEED) and the training stream (NPE_MASTER_SEED). It is a
# REAL GATE: it EXITS NONZERO if any pre-registered rank-uniformity / calibration
# threshold fails.
#
# ANTI-SNOOPING / HONESTY CONTRACT (SBC-03 / D-02, CLAUDE.md HARD_GATE): every M/L/bin
# and every pass/fail threshold was committed to consts.jl BEFORE this run. A genuine
# calibration shortfall (e.g. the documented negative-Δρ tail, RESEARCH Pitfall 6) is an
# HONEST SCIENTIFIC FINDING surfaced to the Phase-6 Go/No-Go memo, NOT tuned away by
# editing consts.jl or reseeding. This script NEVER changes a locked value.
#
# SBC-04 CAPTION: the calibration proof is "calibrated under the simulator; pair with
# the OOD result" — calibration is conditional on the forward model, not unconditional.
#
# CPU-only (D-10), read-only src/ (via the harness). Run:
#     julia --project=spike spike/validation/run_sbc.jl

using Test
using JLD2
using Statistics
using Dates

# Units under test: the SBC diagnostics (pulls harness.jl + consts.jl transitively) and
# the figure surface (CairoMakie, loaded here — figures are script artifacts, not gates).
isdefined(@__MODULE__, :sbc_ranks)           || include(joinpath(@__DIR__, "sbc.jl"))
isdefined(@__MODULE__, :plot_rank_histogram) || include(joinpath(@__DIR__, "figures.jl"))

const SBC_REPORT_PATH = joinpath(@__DIR__, "sbc_report.jld2")

# --- Atomic persistence of the reported SBC artifact (save_npe idiom, T-05-13) ---------
function _save_sbc_report(path; kwargs...)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "rank_table") "_save_sbc_report: integrity check failed ($tmp)"
    end
    mv(tmp, path; force = true)
    return path
end

function main()
    println("="^78)
    println("Phase-5 REPORTED SBC run (SBC-01..04) — LOCKED consts, RESERVED VAL stream")
    println("  M=$SBC_M  L=$SBC_L  bins=$SBC_BINS  imsize=$SBC_IMSIZE")
    println("  KS_α=$SBC_KS_ALPHA  χ²_α=$SBC_CHI2_ALPHA  ECE_green=$SBC_ECE_GREEN")
    println("  stream=VAL_MASTER_SEED=$(repr(VAL_MASTER_SEED))  (≠ NPE_MASTER_SEED=$(repr(NPE_MASTER_SEED)))")
    println("="^78)
    @assert VAL_MASTER_SEED != NPE_MASTER_SEED "reported stream must be disjoint from training stream (D-02)"
    @assert VAL_MASTER_SEED != VAL_FIX_SEED    "reported stream must be disjoint from fixture stream (D-02)"

    m = load_frozen_model()

    # --- The reported M×8 rank table on the RESERVED fresh stream (SBC-01) --------------
    t0 = time()
    rank_table = sbc_ranks(m; M = SBC_M, L = SBC_L, imsize = SBC_IMSIZE,
                           rng = val_rng(VAL_MASTER_SEED))
    elapsed = time() - t0
    println("sbc_ranks: $(size(rank_table,1))×$(size(rank_table,2)) rank table in $(round(elapsed/60; digits=2)) min")

    labels = SBC_PARAM_LABELS
    P = length(labels)                        # 7 θ + Δρ = 8

    # --- Per-parameter uniformity, coverage, ECE/MCE, traffic-light (SBC-02) ------------
    ks_p    = Vector{Float64}(undef, P)
    chi2_p  = Vector{Float64}(undef, P)
    ece     = Vector{Float64}(undef, P)
    mce     = Vector{Float64}(undef, P)
    verdict = Vector{Symbol}(undef, P)
    cov_nominal   = Vector{Vector{Float64}}(undef, P)
    cov_empirical = Vector{Vector{Float64}}(undef, P)

    for p in 1:P
        ranks = @view rank_table[:, p]
        u = sbc_uniformity(ranks; L = SBC_L, bins = SBC_BINS)
        ks_p[p]   = u.ks_p
        chi2_p[p] = u.chi2_p
        cov = sbc_coverage(ranks; L = SBC_L)
        cov_nominal[p]   = cov.nominal
        cov_empirical[p] = cov.empirical
        cal = sbc_calibration(ranks; L = SBC_L, n_bins = SBC_BINS)
        ece[p]     = cal.ece
        mce[p]     = cal.mce
        verdict[p] = sbc_traffic_light(cal.ece)
    end

    # --- Persist the reported artifact BEFORE the gate (so it survives an honest fail) ---
    _save_sbc_report(SBC_REPORT_PATH;
        rank_table = rank_table, labels = collect(labels),
        ks_p = ks_p, chi2_p = chi2_p, ece = ece, mce = mce,
        verdict = String.(verdict),
        cov_nominal = cov_nominal, cov_empirical = cov_empirical,
        SBC_M = SBC_M, SBC_L = SBC_L, SBC_BINS = SBC_BINS,
        SBC_KS_ALPHA = SBC_KS_ALPHA, SBC_CHI2_ALPHA = SBC_CHI2_ALPHA,
        SBC_ECE_GREEN = SBC_ECE_GREEN, SBC_ECE_YELLOW = SBC_ECE_YELLOW,
        SBC_IMSIZE = collect(SBC_IMSIZE),
        VAL_MASTER_SEED = UInt64(VAL_MASTER_SEED),
        caption = SBC_CAPTION,
        elapsed_min = elapsed / 60,
        generated = string(Dates.now(Dates.UTC)) * "Z")
    println("persisted reported artifact -> $SBC_REPORT_PATH")

    # --- Figures (script artifacts — never gate assertions) -----------------------------
    for p in 1:P
        lbl = labels[p]
        plot_rank_histogram(rank_table[:, p]; bins = SBC_BINS, L = SBC_L, label = lbl,
                            filename = "sbc_rank_hist_$(lbl).png")
        plot_coverage_curve(cov_nominal[p], cov_empirical[p];
                            filename = "sbc_coverage_$(lbl).png")
        plot_reliability_ece(sbc_calibration(rank_table[:, p]; L = SBC_L, n_bins = SBC_BINS);
                             filename = "sbc_reliability_$(lbl).png")
    end
    println("figures -> $(joinpath(@__DIR__, "figures"))/sbc_*")

    # --- Per-parameter PASS/FAIL table (printed BEFORE the gate throws) ------------------
    println("\n", "-"^78)
    println(rpad("parameter", 18), rpad("KS p", 12), rpad("χ² p", 12),
            rpad("ECE", 10), rpad("MCE", 10), rpad("verdict", 9), "gate")
    println("-"^78)
    param_pass = Vector{Bool}(undef, P)
    for p in 1:P
        pass = (ks_p[p] > SBC_KS_ALPHA) && (chi2_p[p] > SBC_CHI2_ALPHA) && (ece[p] < SBC_ECE_GREEN)
        param_pass[p] = pass
        println(rpad(labels[p], 18),
                rpad(round(ks_p[p];   sigdigits=3), 12),
                rpad(round(chi2_p[p]; sigdigits=3), 12),
                rpad(round(ece[p];    sigdigits=3), 10),
                rpad(round(mce[p];    sigdigits=3), 10),
                rpad(String(verdict[p]), 9),
                pass ? "PASS" : "FAIL")
    end
    println("-"^78)
    println("SBC caption: \"$SBC_CAPTION\"")
    println("Overall reported SBC gate: ", all(param_pass) ? "PASS" : "FAIL")
    println("-"^78, "\n")

    # --- THE REAL GATE (SBC-01..04): exits nonzero if any pre-registered threshold fails.
    # Thresholds referenced from consts.jl (never tuned). An honest failure is a Phase-6
    # finding, NOT a reason to loosen SBC_KS_ALPHA / SBC_CHI2_ALPHA / SBC_ECE_GREEN.
    @testset "Reported SBC pre-registered gate (SBC-01..04)" begin
        for p in 1:P
            @testset "$(labels[p])" begin
                @test ks_p[p]   > SBC_KS_ALPHA
                @test chi2_p[p] > SBC_CHI2_ALPHA
                @test ece[p]    < SBC_ECE_GREEN
            end
        end
    end
    return nothing
end

main()
