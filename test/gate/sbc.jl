#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/sbc.jl --- per-grid Simulation-Based Calibration for the D-05 ship-gate.
#
# The calibration PROOF for a shipped grid's frozen net "under the simulator". It rides the
# grid-parametrized CPU harness (harness.jl, use_gpu = false) with NO re-fitting and computes,
# from a fresh-seeded (prod_rng[G]) M×L SBC run:
#   • per-parameter rank histograms for all 7 θ PLUS a dedicated paired-draw Δρ SBC (D-01) — an
#     M×8 integer rank table, every rank ∈ 0:L;
#   • KS + χ² rank-uniformity p-values via HypothesisTests (never hand-rolled), a nominal-vs-
#     empirical coverage curve, and ECE/MCE with a traffic-light verdict against the pre-registered
#     SBC_ECE_* cutoffs.
#
# Promoted (grid-generalized, fresh-seeded) from spike/validation/sbc.jl:65-196. The ported
# CalibrationResult / _bin_calibration ECE-MCE pattern is the BayesInteractomics reliability
# diagram (pattern copied, NOT the dependency).
#
# ANTI-SNOOPING (Pitfall 3 / T-7-08): M, L, bin count and all thresholds live in the committed
# gate_consts_<G>.jl; the reported run consumes the FRESH disjoint prod_rng(G) stream (NEVER the
# spike VAL_MASTER_SEED / the training NPE_MASTER_SEED). This file NEVER tunes those values.

using HypothesisTests    # ExactOneSampleKSTest, ChisqTest, pvalue (rank uniformity)
using Distributions      # Uniform (the KS reference)
using Statistics         # mean

isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "harness.jl"))

# --- SBC caption contract --------------------------------------------------------------------
# Calibration is proven CONDITIONAL on the simulator and must be read alongside the OOD flag.
const SBC_CAPTION = "calibrated under the simulator; pair with the OOD result"

# --- Ported reliability-diagram struct + binning (BayesInteractomics pattern) ----------------
"""
    CalibrationResult

Binned reliability-diagram result: per-bin midpoints, mean predicted probability, observed
positive rate, bin counts, plus the Expected (ECE) and Maximum (MCE) Calibration Errors. Ported
from BayesInteractomics `types.jl` (pattern, not a dependency).
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

Bin samples by predicted probability into `n_bins` equal-width bins over [0,1] and compare each
bin's mean predicted probability to its observed positive rate. ECE = Σ (bin_count/total)·|pred−obs|;
MCE = max gap. Ported verbatim from BayesInteractomics `calibration.jl` (pattern, not a dependency).
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

Map an ECE to the pre-registered verdict: `:green` if `ece ≤ SBC_ECE_GREEN`, `:yellow` if
`ece ≤ SBC_ECE_YELLOW`, else `:red`. The cutoffs are locked in the gate consts.
"""
function sbc_traffic_light(ece::Real)
    ece <= SBC_ECE_GREEN  && return :green
    ece <= SBC_ECE_YELLOW && return :yellow
    return :red
end

# --- Parameter labels for the 8 SBC rank columns (7 θ + Δρ) ----------------------------------
const SBC_PARAM_LABELS = ("ρ_true", "spillover", "autofluorescence",
                          "label_efficiency", "shift_dx", "shift_dy", "noise", "Δρ")

"""
    sbc_ranks(m; G, M = SBC_M, L = SBC_L, imsize = SBC_IMSIZE,
              sim = default_simulator(), rng = prod_rng(G)) -> Matrix{Int}

The SBC rank table: an `M×8` integer matrix. For each of `M` prior draws it runs the grid-`G`
draw→simulate→CPU-infer chain and records, for every one of the 7 θ parameters, the rank of θ*ₚ
among its `L` posterior draws (`count(<(θ*ₚ), draws[p,:])`, ∈ `0:L`). Column 8 is the DEDICATED
paired-draw Δρ SBC (D-01): two independent draws give ρ_true draw vectors ρs, ρc; the Δρ rank is
`count(<(Δρ*), ρs .- ρc)` with `Δρ* = θs.ρ_true − θc.ρ_true`. CPU-only; consumes the passed FRESH
disjoint `rng` (prod_rng[G]) sequentially.
"""
function sbc_ranks(m; G::Integer, M::Integer = SBC_M, L::Integer = SBC_L,
                   imsize = SBC_IMSIZE, imsize_set = GATE_IMSIZE_SET,
                   imsize_weights = GATE_IMSIZE_WEIGHTS,
                   sim = default_simulator(), rng = prod_rng(G))
    return sbc_ranks_and_spread(m; G = G, M = M, L = L, imsize = imsize,
                                imsize_set = imsize_set, imsize_weights = imsize_weights,
                                sim = sim, rng = rng).ranks
end

"""
    sbc_ranks_and_spread(m; G, M = SBC_M, L = SBC_L, imsize = SBC_IMSIZE,
                         sim = default_simulator(), rng = prod_rng(G))
        -> (ranks, post_sd, prior_draws)

The SBC rank table PLUS the spread quantities the shrinkage diagnostic needs, from the SAME
draw→simulate→infer chain (`sbc_ranks` is a thin wrapper on this):

  • `ranks`       — the `M×8` rank table (identical to `sbc_ranks`; see its docstring);
  • `post_sd`     — `M×8`, the posterior standard deviation of the `L` draws behind each rank;
  • `prior_draws` — `M×8`, the θ*ₚ actually drawn from π (column 8 = the paired Δρ*);
  • `imsizes`     — length-`M`, the REALISED image size of each marginal (7-θ) draw;
  • `imsizes_delta` — length-`M`, the REALISED (shared) image size of each paired Δρ draw.

Consumes the passed `rng` in EXACTLY the same order and amount as `sbc_ranks` did — the extra
outputs are `std`/bookkeeping over samples that were already drawn — so the rank table, and every
verdict computed from it, is unchanged.

F5 MIXTURE (07-GATE-AMENDMENT §4, required change 3): under the amended pre-registration the image
size is drawn PER DRAW from `imsize_set`/`imsize_weights` off THIS `rng` (so it is part of the
reproducible stream), and the realised size is recorded alongside the rank table so the realised
mixture is auditable post-hoc rather than reconstructed from the weights.
"""
function sbc_ranks_and_spread(m; G::Integer, M::Integer = SBC_M, L::Integer = SBC_L,
                              imsize = SBC_IMSIZE, imsize_set = GATE_IMSIZE_SET,
                              imsize_weights = GATE_IMSIZE_WEIGHTS,
                              sim = default_simulator(),
                              rng = prod_rng(G))
    # REPRODUCIBILITY (07-10): pin BOTH streams. `rng` (prod_rng(G)) drives the prior draw + the
    # forward simulation; the flow's posterior draws come from the GLOBAL RNG because
    # `sampleposterior` threads no rng (see the GLOBAL-RNG SEEDING block in harness.jl). Seeding
    # here — deterministically from the SAME frozen PROD_SEED[G] — is what makes this rank table
    # reproducible bit-for-bit on a re-run. NOTE: the grid-4/8/16 tables already recorded in
    # artifacts/ predate this fix and were drawn from an UNSEEDED global stream; they are NOT
    # retroactively reproducible, and a re-run will legitimately differ from them.
    seed_gate_global!(G)
    rank_table  = Matrix{Int}(undef, M, 8)
    post_sd     = Matrix{Float64}(undef, M, 8)
    prior_draws = Matrix{Float64}(undef, M, 8)
    imsizes       = Vector{Tuple{Int,Int}}(undef, M)   # realised size, marginal draws
    imsizes_delta = Vector{Tuple{Int,Int}}(undef, M)   # realised (shared) size, paired Δρ draws
    for i in 1:M
        # 7 marginal θ ranks from one draw→simulate→infer.
        t = draw_simulate_infer(m, rng; G = G, imsize = imsize, imsize_set = imsize_set,
                                imsize_weights = imsize_weights, N = L, sim = sim)
        imsizes[i] = t.imsize
        for p in 1:7
            rank_table[i, p]  = count(<(t.θ[p]), @view t.draws[p, :])
            post_sd[i, p]     = Statistics.std(@view t.draws[p, :])
            prior_draws[i, p] = t.θ[p]
        end
        # Δρ rank from a paired draw (D-01), independent of the marginal ρ_true rank.
        pr       = draw_simulate_infer_paired(m, rng; G = G, imsize = imsize,
                                              imsize_set = imsize_set,
                                              imsize_weights = imsize_weights, N = L, sim = sim)
        imsizes_delta[i] = pr.imsize
        Δρ_draws = pr.ρs .- pr.ρc
        Δρ_star  = pr.θs.ρ_true - pr.θc.ρ_true
        rank_table[i, 8]  = count(<(Δρ_star), Δρ_draws)
        post_sd[i, 8]     = Statistics.std(Δρ_draws)
        prior_draws[i, 8] = Δρ_star
    end
    return (ranks = rank_table, post_sd = post_sd, prior_draws = prior_draws,
            imsizes = imsizes, imsizes_delta = imsizes_delta)
end

# --- Shrinkage (vacuous-pass) diagnostic — REPORTING ONLY ------------------------------------
#
# WHY (07-CALIBRATION-FINDINGS F3): SBC rank uniformity is a NECESSARY but not SUFFICIENT
# calibration signal. A posterior that simply REPRODUCES THE PRIOR — i.e. the net learned nothing
# about that parameter — yields perfectly uniform ranks BY CONSTRUCTION. The grid-4
# `label_efficiency` column is exactly this case (post_sd/prior_sd = 0.994, shrinkage centre =
# the prior mean), and its KS pass must never be read as evidence of calibration.
#
# The fix is to make the vacuous case VISIBLE, not to change the verdict: `sbc_gate` now reports a
# per-parameter `shrinkage = mean(post_sd) / prior_sd` alongside every p-value, plus a boolean
# `vacuous` annotation. Both are DIAGNOSTIC ONLY — neither enters `ks_pass`, `ece_pass` or
# `passed`, and no pre-registered `gate_consts_<G>.jl` threshold is touched (that file is frozen
# pre-registration; adding a pass condition here would be post-hoc re-registration).
#
# `SBC_VACUOUS_SHRINKAGE` is therefore a REPORTING LABEL CUTOFF, deliberately defined here and NOT
# in the frozen consts: shrinkage ≥ 0.95 means the posterior retains ≥95% of the prior spread, i.e.
# essentially no information was gained about that parameter.
const SBC_VACUOUS_SHRINKAGE = 0.95

"""
    sbc_shrinkage(post_sd, prior_draws) -> (shrinkage, post_sd_mean, prior_sd, vacuous)

Per-parameter posterior-to-prior spread ratio for the 8 SBC columns. `prior_sd` is the standard
deviation of the M θ* values ACTUALLY drawn from π in this run (an empirical prior sample — no
analytic form is assumed, which matters because ρ_true = `ghat(μ*)` has no closed-form sd);
`post_sd_mean` is the mean over the M draws of the posterior sd; `shrinkage` is their ratio.

`shrinkage ≈ 1` ⇒ the posterior is the prior ⇒ that column's rank uniformity is VACUOUS and is
flagged `vacuous = true` (at the `SBC_VACUOUS_SHRINKAGE` reporting cutoff). Diagnostic only —
never a pass/fail input.
"""
function sbc_shrinkage(post_sd::AbstractMatrix, prior_draws::AbstractMatrix)
    n = size(post_sd, 2)
    psd_mean = [Statistics.mean(view(post_sd, :, p)) for p in 1:n]
    prior_sd = [Statistics.std(view(prior_draws, :, p)) for p in 1:n]
    shrink   = [prior_sd[p] > 0 ? psd_mean[p] / prior_sd[p] : NaN for p in 1:n]
    vacuous  = [isfinite(s) && s >= SBC_VACUOUS_SHRINKAGE for s in shrink]
    return (shrinkage = shrink, post_sd_mean = psd_mean, prior_sd = prior_sd, vacuous = vacuous)
end

"""
    sbc_uniformity(ranks; L = SBC_L, bins = SBC_BINS) -> (ks_p, chi2_p)

Rank-uniformity p-values for one parameter's SBC ranks. Maps ranks to PIT values
`u = (ranks .+ 0.5) ./ (L+1)` and runs `ExactOneSampleKSTest(u, Uniform(0,1))`; independently bins
the ranks into `bins` equal-width rank bins (L+1 divisible by bins) and runs the one-way
`ChisqTest(counts)`. A well-calibrated net gives both p > the pre-registered alphas. Never
hand-rolled.
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

The nominal-vs-empirical coverage curve. For each nominal central-credible level α, the empirical
coverage is the fraction of draws whose PIT `u = (rank+0.5)/(L+1)` falls in the central α interval
`|u − 0.5| ≤ α/2`. Returns parallel `nominal`/`empirical` vectors.
"""
function sbc_coverage(ranks::AbstractVector{<:Integer};
                      levels = SBC_COVERAGE_LEVELS, L::Integer = SBC_L)
    u = (ranks .+ 0.5) ./ (L + 1)
    nominal   = collect(levels)
    empirical = [mean(abs.(u .- 0.5) .<= (α / 2)) for α in nominal]
    return (nominal = nominal, empirical = empirical)
end

"""
    sbc_calibration(ranks; levels = SBC_COVERAGE_LEVELS, L = SBC_L, n_bins = SBC_BINS) -> CalibrationResult

Route the SBC central-credible-interval coverage into the reliability-diagram binning (the ECE/MCE
calibration metric). For every (draw, nominal-level) pair the predicted probability is the nominal
level α and the empirical positive is the in-interval indicator `|u − 0.5| ≤ α/2`; `_bin_calibration`
bins by predicted level and compares to observed coverage. The resulting ECE is scored by
`sbc_traffic_light`.
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

# =============================================================================================
# A1 — MULTIPLICITY CONTROL (07-GATE-AMENDMENT §1) and the v1↔v2 rules-change attribution
# =============================================================================================
#
# THE DEFECT. v1's `ks_pass = all(x -> x.ks_p > SBC_KS_ALPHA, per)` runs 8 level-0.05 tests and
# requires ALL to survive. A PERFECTLY CALIBRATED net fails that conjunction with probability
# 1 − 0.95^8 = 33.7 % — the gate's true α was 0.34, not the 0.05 it claimed. This is arithmetic on
# the gate's own structure (8 tests, α = 0.05), computable with the artifacts directory deleted.
#
# THE RULE. Holm–Bonferroni step-down at a FAMILY-WISE error rate of 0.05, applied SEPARATELY to
# the 8 KS p-values and to the 8 χ² p-values. Holm is valid under ARBITRARY dependence (the 8
# columns share the same M datasets and the same posterior draws, so they are dependent by
# construction with an unknown structure) and is uniformly at least as powerful as plain Bonferroni.
#
# ATTRIBUTION (the v1 verdict is ALWAYS co-computed). The amended run changes TWO things at once —
# the gate rules AND the model (the bounded-θ retrain; amendment §5, "the bounded-θ retrain is a
# confound"). Reporting only the amended verdict would make the rules-change contribution
# unrecoverable. Every `sbc_gate` result therefore carries BOTH verdicts computed on the SAME rank
# table: `v1_verdict` (the original uncorrected all-8 α = 0.05 conjunction) and `v2_verdict`
# (Holm). The BINDING verdict is whichever the LOADED pre-registration prescribes, and is reported
# as `rules_version`. The v1 numbers are free — they are a second reduction of one rank table, not
# a second run — and they cost nothing in seed budget.

"""
    sbc_holm_adjusted(p) -> Vector{Float64}

Holm–Bonferroni adjusted p-values for a family of `m = length(p)` tests, in INPUT order: sort
ascending, multiply p₍ⱼ₎ by `(m − j + 1)`, take the running maximum (monotonicity) and cap at 1.
Valid under arbitrary dependence.

The FROZEN `holm_adjusted` in `gate_consts_8_v2.jl` is used instead whenever the amended
pre-registration is loaded (`sbc_holm` below dispatches to it); this implementation is the fallback
that lets the v1 consts still compute the amended verdict for the side-by-side attribution. A unit
test asserts the two agree.
"""
function sbc_holm_adjusted(p::AbstractVector{<:Real})
    m   = length(p)
    ord = sortperm(collect(float.(p)))
    adj = Vector{Float64}(undef, m)
    running = 0.0
    for (j, i) in enumerate(ord)
        running = max(running, (m - j + 1) * float(p[i]))
        adj[i]  = min(1.0, running)
    end
    return adj
end

"Holm adjusted p-values, PREFERRING the frozen `holm_adjusted` of the loaded pre-registration."
sbc_holm(p) = isdefined(@__MODULE__, :holm_adjusted) ?
    getfield(@__MODULE__, :holm_adjusted)(p) : sbc_holm_adjusted(p)

"""
    sbc_holm_pass(p; fwer = 0.05) -> Bool

`true` iff Holm–Bonferroni at family-wise error rate `fwer` rejects NO null in the family — i.e.
every adjusted p-value exceeds `fwer`. This is the amended (A1) uniformity verdict.
"""
sbc_holm_pass(p::AbstractVector{<:Real}; fwer::Real = 0.05) = all(>(fwer), sbc_holm(p))

# The v1 per-test α, used ONLY to co-report what the ORIGINAL rules would have concluded. Under a
# v1 consts file this is the frozen `SBC_KS_ALPHA` itself; under v2 that constant is deliberately
# RETIRED (undefined so that any code still gating on it errors loudly), so the historical value
# 0.05 is named here — as a REPORTING constant for the attribution, never as a gate input.
const SBC_V1_ALPHA_ATTRIBUTION = 0.05
_sbc_v1_alpha() = isdefined(@__MODULE__, :SBC_KS_ALPHA) ?
    getfield(@__MODULE__, :SBC_KS_ALPHA) : SBC_V1_ALPHA_ATTRIBUTION

"The family-wise α of the loaded pre-registration (v2), or the v1 default when it is not loaded."
_sbc_ks_fwer()   = isdefined(@__MODULE__, :SBC_KS_FWER)   ? getfield(@__MODULE__, :SBC_KS_FWER)   : 0.05
_sbc_chi2_fwer() = isdefined(@__MODULE__, :SBC_CHI2_FWER) ? getfield(@__MODULE__, :SBC_CHI2_FWER) : 0.05

"`2` when the amended pre-registration (`SBC_MULTIPLICITY = :holm`) is loaded, else `1`."
const SBC_RULES_VERSION = isdefined(@__MODULE__, :SBC_MULTIPLICITY) &&
                          SBC_MULTIPLICITY === :holm ? 2 : 1

# A1b — the χ² rank-bin count, DECOUPLED from the ECE reliability binning. v1 overloaded ONE
# constant (SBC_BINS = 50) for both, so neither could be set on its own merits. Under v2,
# `SBC_CHI2_BINS = 20` drives `sbc_uniformity` while `SBC_BINS = 50` continues to drive
# `sbc_calibration` — the ECE arm is therefore NUMERICALLY IDENTICAL to v1. Under a v1 consts file
# there is no separate constant and the χ² binning falls back to `bins`, i.e. v1 behaviour.
const GATE_CHI2_BINS = isdefined(@__MODULE__, :SBC_CHI2_BINS) ? SBC_CHI2_BINS : nothing
_gate_chi2_bins(bins) = GATE_CHI2_BINS === nothing ? bins : GATE_CHI2_BINS

"""
    sbc_gate(m; G, M = SBC_M, L = SBC_L, bins = SBC_BINS, chi2_bins = _gate_chi2_bins(bins),
             imsize = SBC_IMSIZE, imsize_set = GATE_IMSIZE_SET,
             imsize_weights = GATE_IMSIZE_WEIGHTS, sim = default_simulator(),
             rng = prod_rng(G)) -> NamedTuple

The aggregated per-grid SBC ship-gate verdict: runs the fresh-seeded M×L rank table, then for
each of the 8 columns (7 θ + Δρ) computes KS/χ² uniformity, the ECE/MCE `CalibrationResult`, and
its traffic-light. CPU-only, pre-registered thresholds, fresh disjoint stream.

**BINNING (A1b).** `chi2_bins` drives `sbc_uniformity` (the χ² rank bins) and `bins` drives
`sbc_calibration` (the ECE/MCE reliability bins). Under the amended pre-registration these are 20
and 50; under a v1 consts file `chi2_bins` defaults to `bins`, exactly reproducing v1.

**MULTIPLICITY (A1).** BOTH verdicts are computed on the SAME rank table and reported side by side:

  • `v1_verdict` — the ORIGINAL uncorrected conjunction `all(ks_p > 0.05)`, whose family-wise
    false-failure rate for a perfectly calibrated net is 1 − 0.95^8 = 33.7 %;
  • `v2_verdict` — Holm–Bonferroni step-down at FWER 0.05, applied separately to the KS and χ²
    families.

`rules_version` names which one is BINDING (the one the loaded pre-registration prescribes), and
`passed`/`ks_pass`/`ece_pass` are that one's fields. This co-reporting is a free attribution of the
rules change: the amended run alters both the rules and the model (bounded-θ retrain), and without
the v1 numbers on the same table the two contributions could not be separated (amendment §5).

`passed` retains v1's STRUCTURE — `ks_pass && ece_pass`. The χ² family verdict is computed and
reported (`chi2_pass`) under the same Holm correction, but is NOT folded into `passed`: v1 never
gated on χ², and the amendment specifies a multiplicity correction for the existing arms, not a new
gating arm. Adding one would be an unspecified tightening.

Each `per_param` entry ALSO carries the diagnostic `shrinkage = post_sd/prior_sd` (with its
`post_sd`/`prior_sd` components) and a `vacuous` flag, and the report carries the list of
`vacuous_params`. These make a "posterior == prior" column — which passes rank uniformity BY
CONSTRUCTION while having learned nothing (F3) — impossible to mis-read as calibration evidence.
They are REPORTING ONLY and feed no pass/fail decision.
"""
function sbc_gate(m; G::Integer, M::Integer = SBC_M, L::Integer = SBC_L,
                  bins::Integer = SBC_BINS, chi2_bins::Integer = _gate_chi2_bins(bins),
                  imsize = SBC_IMSIZE, imsize_set = GATE_IMSIZE_SET,
                  imsize_weights = GATE_IMSIZE_WEIGHTS,
                  sim = default_simulator(), rng = prod_rng(G))
    sr    = sbc_ranks_and_spread(m; G = G, M = M, L = L, imsize = imsize,
                                 imsize_set = imsize_set, imsize_weights = imsize_weights,
                                 sim = sim, rng = rng)
    ranks = sr.ranks
    shr   = sbc_shrinkage(sr.post_sd, sr.prior_draws)
    per = NamedTuple[]
    for p in 1:8
        # A1b: χ² gets its OWN bin count; the ECE reliability binning is untouched.
        u   = sbc_uniformity(view(ranks, :, p); L = L, bins = chi2_bins)
        cal = sbc_calibration(collect(view(ranks, :, p)); L = L, n_bins = bins)
        push!(per, (label = SBC_PARAM_LABELS[p], ks_p = u.ks_p, chi2_p = u.chi2_p,
                    ece = cal.ece, mce = cal.mce, light = sbc_traffic_light(cal.ece),
                    # --- diagnostic only; NOT inputs to any pass/fail decision (F3) ---
                    shrinkage = shr.shrinkage[p], post_sd = shr.post_sd_mean[p],
                    prior_sd = shr.prior_sd[p], vacuous = shr.vacuous[p]))
    end

    ks_p     = [x.ks_p   for x in per]
    chi2_p   = [x.chi2_p for x in per]
    ece_pass = all(x -> x.ece <= SBC_ECE_GREEN, per)     # UNCHANGED in both rule sets

    # --- v1 (original) rules, recomputed on THIS rank table purely for attribution -------------
    α  = _sbc_v1_alpha()
    v1 = (rules = 1, alpha = α,
          ks_pass   = all(>(α), ks_p),
          chi2_pass = all(>(α), chi2_p),
          ece_pass  = ece_pass)
    v1 = merge(v1, (passed = v1.ks_pass && v1.ece_pass,))

    # --- v2 (amended) rules: Holm–Bonferroni at FWER 0.05 on each family -----------------------
    ks_adj   = sbc_holm(ks_p)
    chi2_adj = sbc_holm(chi2_p)
    v2 = (rules = 2, multiplicity = :holm, ks_fwer = _sbc_ks_fwer(), chi2_fwer = _sbc_chi2_fwer(),
          ks_adj_p = ks_adj, chi2_adj_p = chi2_adj,
          ks_pass   = all(>(_sbc_ks_fwer()),   ks_adj),
          chi2_pass = all(>(_sbc_chi2_fwer()), chi2_adj),
          ece_pass  = ece_pass)
    v2 = merge(v2, (passed = v2.ks_pass && v2.ece_pass,))

    binding = SBC_RULES_VERSION == 2 ? v2 : v1

    return (grid = Int(G), M = Int(M), L = Int(L), ranks = ranks, per_param = per,
            post_sd = sr.post_sd, prior_draws = sr.prior_draws,
            bins = Int(bins), chi2_bins = Int(chi2_bins),
            imsize_set = imsize_set, imsize_weights = imsize_weights,
            realised_imsize_counts =
                realised_imsize_counts(vcat(sr.imsizes, sr.imsizes_delta)),
            imsizes = sr.imsizes, imsizes_delta = sr.imsizes_delta,
            vacuous_cutoff = SBC_VACUOUS_SHRINKAGE,
            vacuous_params = [SBC_PARAM_LABELS[p] for p in 1:8 if shr.vacuous[p]],
            rules_version = SBC_RULES_VERSION, v1_verdict = v1, v2_verdict = v2,
            ks_pass = binding.ks_pass, chi2_pass = binding.chi2_pass,
            ece_pass = ece_pass, passed = binding.passed,
            caption = SBC_CAPTION)
end
