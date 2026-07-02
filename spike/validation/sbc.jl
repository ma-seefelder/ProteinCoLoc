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

# spike/validation/sbc.jl --- Phase-5 Simulation-Based Calibration (SBC-01..04).
#
# The highest-scrutiny scientific-honesty deliverable: the calibration PROOF for the
# frozen Phase-4 net "under the simulator". It rides the shared harness (no re-fitting)
# and computes, from a pre-registered M×L SBC run:
#   • per-parameter rank histograms for all 7 θ PLUS a dedicated paired-draw Δρ SBC
#     (D-01) -- an M×8 integer rank table, every rank ∈ 0:L (SBC-01);
#   • KS + χ² rank-uniformity p-values via HypothesisTests (never hand-rolled), a
#     nominal-vs-empirical coverage curve, and ECE/MCE with a traffic-light verdict
#     against the pre-registered SBC_ECE_* cutoffs (SBC-02);
#   • the SBC_CAPTION contract: results are "calibrated UNDER THE SIMULATOR" and must
#     be paired with the OOD result (SBC-04) -- calibration is conditional on the
#     forward model, not an unconditional guarantee.
#
# ANTI-SNOOPING (SBC-03 / D-02): M, L, bin count and all thresholds live in the
# committed consts.jl; the reported run consumes the reserved disjoint VAL_MASTER_SEED
# stream. This file NEVER tunes those values.
#
# ECE/MCE: `_bin_calibration` + `CalibrationResult` are ported from BayesInteractomics
# (calibration.jl / types.jl) -- the PATTERN is copied (~50 lines), NOT the dependency.
#
# CPU-only, read-only src/ (via the harness). Flat top-level functions; guarded includes.

using HypothesisTests    # ExactOneSampleKSTest, ChisqTest, pvalue (rank uniformity)
using Distributions      # Uniform (the KS reference)
using Statistics         # mean
using StatsBase          # (reconstruct etc. come through the harness)

isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "harness.jl"))

# --- SBC-04 caption contract ----------------------------------------------------
# The EXACT string every SBC figure/report is captioned with: calibration is proven
# CONDITIONAL on the simulator and must be read alongside the OOD misspecification
# flag (a well-specified calibration proof says nothing about out-of-distribution data).
const SBC_CAPTION = "calibrated under the simulator; pair with the OOD result"

# --- Ported reliability-diagram struct + binning (BayesInteractomics pattern) ----
"""
    CalibrationResult

Binned reliability-diagram result: per-bin midpoints, mean predicted probability,
observed positive rate, bin counts, plus the Expected (ECE) and Maximum (MCE)
Calibration Errors. Ported from BayesInteractomics `types.jl` (pattern, not a dep).
"""
struct CalibrationResult
    bin_midpoints::Vector{Float64}
    predicted_rate::Vector{Float64}
    observed_rate::Vector{Float64}
    bin_counts::Vector{Int}
    ece::Float64
    mce::Float64
end

"""
    _bin_calibration(posterior_probs, empirical_positive; n_bins = 10) -> CalibrationResult

Shared binning logic (ported verbatim from BayesInteractomics `calibration.jl`, the
PATTERN not the dependency): bin samples by predicted probability into `n_bins`
equal-width bins over [0,1] and compare each bin's mean predicted probability to its
observed positive rate. ECE = Σ (bin_count/total)·|pred−obs|; MCE = max gap.
"""
function _bin_calibration(posterior_probs::AbstractVector{<:Real},
                          empirical_positive::AbstractVector{Bool};
                          n_bins::Int = 10)
    bin_edges = range(0.0, 1.0, length = n_bins + 1)

    bin_midpoints  = Float64[]
    predicted_rate = Float64[]
    observed_rate  = Float64[]
    bin_counts     = Int[]

    for i in 1:n_bins
        lo = bin_edges[i]
        hi = bin_edges[i + 1]

        # Include the right endpoint for the last bin.
        if i == n_bins
            mask = (posterior_probs .>= lo) .& (posterior_probs .<= hi)
        else
            mask = (posterior_probs .>= lo) .& (posterior_probs .< hi)
        end

        count_in_bin = sum(mask)
        push!(bin_counts, count_in_bin)

        midpoint = (lo + hi) / 2.0
        push!(bin_midpoints, midpoint)

        if count_in_bin > 0
            push!(predicted_rate, mean(posterior_probs[mask]))
            push!(observed_rate, mean(empirical_positive[mask]))
        else
            push!(predicted_rate, midpoint)
            push!(observed_rate, 0.0)
        end
    end

    total = sum(bin_counts)
    ece = 0.0
    mce = 0.0
    for i in eachindex(bin_midpoints)
        gap    = abs(predicted_rate[i] - observed_rate[i])
        weight = bin_counts[i] / max(total, 1)
        ece   += weight * gap
        mce    = max(mce, gap)
    end

    return CalibrationResult(bin_midpoints, predicted_rate, observed_rate,
                             bin_counts, ece, mce)
end

"""
    sbc_traffic_light(ece::Real) -> Symbol

Map an ECE to the pre-registered traffic-light verdict: `:green` if
`ece ≤ SBC_ECE_GREEN`, `:yellow` if `ece ≤ SBC_ECE_YELLOW`, else `:red`. The cutoffs
are locked in consts.jl (SBC-03).
"""
function sbc_traffic_light(ece::Real)
    ece <= SBC_ECE_GREEN  && return :green
    ece <= SBC_ECE_YELLOW && return :yellow
    return :red
end

"""
    sbc_ranks(m; M = SBC_M, L = SBC_L, imsize = SBC_IMSIZE, rng = val_rng()) -> Matrix{Int}

The SBC rank table (SBC-01): an `M×8` integer matrix. For each of `M` prior draws it
runs the shared draw→simulate→infer chain and records, for every one of the 7 θ
parameters, the rank of θ*ₚ among its `L` posterior draws (`count(<(θ*ₚ), draws[p,:])`,
∈ `0:L`). Column 8 is the DEDICATED paired-draw Δρ SBC (D-01): two independent draws
give ρ_true draw vectors ρs, ρc; the Δρ rank is `count(<(Δρ*), ρs .- ρc)` with
`Δρ* = θs.ρ_true − θc.ρ_true` -- NOT inherited from the ρ_true marginal ranks.
CPU-only; consumes the passed `rng` sequentially (a single disjoint stream, D-02).
"""
function sbc_ranks(m; M::Integer = SBC_M, L::Integer = SBC_L,
                   imsize = SBC_IMSIZE, rng = val_rng())
    rank_table = Matrix{Int}(undef, M, 8)
    for i in 1:M
        # 7 marginal θ ranks from one draw→simulate→infer.
        t = draw_simulate_infer(m, rng; imsize = imsize, N = L)
        for p in 1:7
            rank_table[i, p] = count(<(t.θ[p]), @view t.draws[p, :])
        end
        # Δρ rank from a paired draw (D-01), independent of the marginal ρ_true rank.
        pr        = draw_simulate_infer_paired(m, rng; imsize = imsize, N = L)
        Δρ_draws  = pr.ρs .- pr.ρc
        Δρ_star   = pr.θs.ρ_true - pr.θc.ρ_true
        rank_table[i, 8] = count(<(Δρ_star), Δρ_draws)
    end
    return rank_table
end

"""
    sbc_uniformity(ranks; L = SBC_L, bins = SBC_BINS) -> (ks_p, chi2_p)

Rank-uniformity p-values for one parameter's SBC ranks (SBC-02). Maps ranks to PIT
values `u = (ranks .+ 0.5) ./ (L+1)` and runs `ExactOneSampleKSTest(u, Uniform(0,1))`;
independently bins the ranks into `bins` equal-width rank bins (L+1 divisible by bins,
Pitfall 4) and runs the one-way `ChisqTest(counts)` (goodness-of-fit to uniform bin
counts). A well-calibrated net gives both p > the pre-registered alphas; a biased rank
vector (e.g. all zeros) drives both p → 0. Never hand-rolled (RESEARCH Don't-Hand-Roll).
"""
function sbc_uniformity(ranks::AbstractVector{<:Integer};
                        L::Integer = SBC_L, bins::Integer = SBC_BINS)
    u = (ranks .+ 0.5) ./ (L + 1)
    ks_p = pvalue(ExactOneSampleKSTest(u, Uniform(0.0, 1.0)))

    counts = zeros(Int, bins)
    for r in ranks
        b = min(bins, floor(Int, (r / (L + 1)) * bins) + 1)
        counts[b] += 1
    end
    chi2_p = pvalue(ChisqTest(counts))
    return (ks_p = ks_p, chi2_p = chi2_p)
end

# Default nominal central-credible levels for the coverage curve.
const SBC_COVERAGE_LEVELS = collect(0.05:0.05:0.95)

"""
    sbc_coverage(ranks; levels = SBC_COVERAGE_LEVELS, L = SBC_L) -> (nominal, empirical)

The nominal-vs-empirical coverage curve (SBC-02). For each nominal central-credible
level α, the empirical coverage is the fraction of draws whose PIT
`u = (rank+0.5)/(L+1)` falls in the central α interval `|u − 0.5| ≤ α/2`. Under a
well-calibrated net the empirical curve tracks the nominal diagonal (monotone,
near-diagonal). Returns parallel `nominal`/`empirical` vectors.
"""
function sbc_coverage(ranks::AbstractVector{<:Integer};
                      levels = SBC_COVERAGE_LEVELS, L::Integer = SBC_L)
    u = (ranks .+ 0.5) ./ (L + 1)
    nominal   = collect(levels)
    empirical = [mean(abs.(u .- 0.5) .<= (α / 2)) for α in nominal]
    return (nominal = nominal, empirical = empirical)
end

"""
    sbc_calibration(ranks; levels = SBC_COVERAGE_LEVELS, L = SBC_L,
                    n_bins = SBC_BINS) -> CalibrationResult

Route the SBC central-credible-interval coverage into the reliability-diagram binning
(the ECE/MCE calibration metric). For every (draw, nominal-level) pair the predicted
probability is the nominal level α and the empirical positive is the in-interval
indicator `|u − 0.5| ≤ α/2`; `_bin_calibration` then bins by predicted level and
compares to observed coverage. The resulting ECE is the coverage miscalibration,
scored by `sbc_traffic_light`.
"""
function sbc_calibration(ranks::AbstractVector{<:Integer};
                         levels = SBC_COVERAGE_LEVELS, L::Integer = SBC_L,
                         n_bins::Integer = SBC_BINS)
    u = (ranks .+ 0.5) ./ (L + 1)
    probs = Float64[]
    pos   = Bool[]
    for α in levels
        for ui in u
            push!(probs, α)
            push!(pos, abs(ui - 0.5) <= (α / 2))
        end
    end
    return _bin_calibration(probs, pos; n_bins = n_bins)
end

# --- Parameter labels for the 8 SBC rank columns (7 θ + Δρ) ---------------------
const SBC_PARAM_LABELS = ("ρ_true", "spillover", "autofluorescence",
                          "label_efficiency", "shift_dx", "shift_dy", "noise", "Δρ")
