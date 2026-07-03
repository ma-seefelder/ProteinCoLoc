#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder

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

#############################################################################################
# src/amortized/bf.jl --- the amortized Bayes-factor READ surface + non-clamped KDE baseline.
#
# Promoted (grid-generalized, CPU-default) from the proven spike:
#   spike/validation/bf.jl:94-100          amortized_log_bf (one-pass NRE log Bayes factor)
#   spike/validation/train_ratio.jl:185-188 pair_encode (research A7 difference encoding)
#
# THE AMORTIZED LOG BAYES FACTOR (BF-01/BF-02): read in ONE forward pass from the trained
# model-comparison `RatioEstimator`:
#
#     log-BF(Z) = logratio(Z; m=1) − logratio(Z; m=0) − log_prior_odds
#
# The NRE ratio identity `logratio(Z,m) = log[p(m|Z)/p(m)]` makes the difference already equal
# `log(posterior-odds) − log(prior-odds) = log-BF`; the MEASURED `log_prior_odds` is subtracted
# (Pitfall 5 — measured, never assumed 0). NO quadgk / KDE / shuffle on this amortized side.
#
# HONESTY NOTE (05-RESEARCH Pitfall 3): NeuralEstimators v0.2.1's `RatioEstimator` loss is
# hard-coded `logitbinarycrossentropy`; a custom l-POP loss passed to `train` is SILENTLY
# IGNORED. This is the standard Hermans-2020 NRE used as a model-comparison reduction — do NOT
# attempt an l-POP loss here.
#
# GRID COUPLING: `pair_encode` derives the continuous-row count from the summary length
# (`nc = length(Zs) ÷ 2 = G²`), so the paired encoding is `ratio_input_dim(G) = 5·G²`-long for
# any shipped grid, matching `summary.jl`'s single-source `ratio_input_dim`.
#
# NON-CLAMPED BF BASELINE (Memo §5 hardening / T-7-06): `kde_log_bf_unclamped` mirrors the
# demoted `compute_BayesFactor` KDE math but WITHOUT the spike baseline's `_clampp = 1e-8`
# probability floor (spike/baseline/run_bf_baseline.jl:57). The clamp saturates `logBF` to
# ±log((1−1e-8)/1e-8) ≈ ±18.4 at any sweep tail where the sharpened NPE Δρ posterior is fully
# one-sided, so `max|Δ logBF|` measured against the CLAMPED baseline reflects the FLOOR, not the
# estimator (STATE Phase-5-iter2). This non-clamped reference lets the ship-gate BF comparison
# measure `max|Δ logBF|` free of the clamp artifact. Uses KernelDensity/QuadGK/Distributions —
# all already root `[deps]` (no re-resolve of the Wave-0 Manifest). CPU-only, gate-side only.
#############################################################################################

import NeuralEstimators: logratio
import KernelDensity: kde
import QuadGK: quadgk
import Distributions: pdf

# ============================ amortized (NRE) log Bayes factor ==============================

"""
    amortized_log_bf(est_bf, Z_pair, log_prior_odds; use_gpu::Bool = false) -> Float64

The amortized log Bayes factor for one paired summary `Z_pair` (a `ratio_input_dim(G)`-vector or
`ratio_input_dim(G)×1` matrix), in ONE forward pass (BF-01/BF-02):

    logratio(Z; m=1) − logratio(Z; m=0) − log_prior_odds

`grid = [0 1]` evaluates the ratio at the null (m=0) and coloc (m=1) model indices; the MEASURED
`log_prior_odds` is subtracted (Pitfall 5). Uses NeuralEstimators `logratio` ONLY — no quadgk /
KDE / shuffle. `use_gpu = false` by default (CPU-reproducible shipped path, D-06). Returns a
finite scalar.
"""
function amortized_log_bf(est_bf, Z_pair, log_prior_odds; use_gpu::Bool = false)
    Zc   = Z_pair isa AbstractVector ? reshape(Z_pair, :, 1) : Z_pair
    grid = reshape(Float32[0.0 1.0], 1, 2)                       # null (m=0), coloc (m=1)
    lr   = logratio(est_bf, Zc; grid = grid, use_gpu = use_gpu)  # 1×2: [logr(m=0) logr(m=1)]
    return (lr[1, 2] - lr[1, 1]) - log_prior_odds
end

"""
    pair_encode(Zs, Zc) -> Vector{Float32}

The paired-summary encoding for the model-comparison ratio net (research A7 difference
encoding): `vcat(Zs, Zc, Zs[1:G²] − Zc[1:G²])` — the two `2·G²`-dim frozen-`zt` summaries
concatenated PLUS the `G²`-dim continuous-row (correlation) CONTRAST. The continuous-row count
`nc = G²` is derived from the summary length (`length(Zs) ÷ 2`), so the encoding is grid-general
and its length is `ratio_input_dim(G) = 5·G²`. Used IDENTICALLY by the trainer and the amortized
BF read, so training and inference share one input space.
"""
function pair_encode(Zs::AbstractVector, Zc::AbstractVector)
    length(Zs) == length(Zc) ||
        throw(DimensionMismatch("pair_encode: Zs ($(length(Zs))) and Zc ($(length(Zc))) " *
                                "must be the same length"))
    iseven(length(Zs)) ||
        throw(ArgumentError("pair_encode: summary length $(length(Zs)) must be 2·G² (even)"))
    nc = length(Zs) ÷ 2                                   # continuous rows = G² (cont_rows(G))
    d  = @view(Zs[1:nc]) .- @view(Zc[1:nc])
    return Float32.(vcat(Zs, Zc, d))
end

# ============================ non-clamped KDE BF baseline (gate only) =======================

"""
    _p_gt_threshold_unclamped(draws; threshold = 0.0) -> Float64

`P(Δρ > threshold)` via the `compute_BayesFactor` recipe (`1 − ∫_{-∞}^{threshold} kde(draws)`)
WITHOUT the spike baseline's `_clampp` 1e-8 floor: the raw KDE tail probability is returned even
when it underflows toward 0 or 1. This is the whole point of the non-clamped baseline — a fully
one-sided posterior yields a probability at the KDE's own numerical floor, not an artificial
1e-8 clamp, so a downstream `logBF` reflects the estimator rather than the clamp (T-7-06).
"""
function _p_gt_threshold_unclamped(draws; threshold::Real = 0.0)
    dist = kde(collect(float.(draws)))
    p_le, _ = quadgk(x -> pdf(dist, x), -Inf, threshold)
    return 1 - p_le                                       # NO _clampp floor (Memo §5)
end

"""
    kde_log_bf_unclamped(post_draws, prior_draws; threshold = 0.0) -> Vector{Float64}

The KDE log Bayes factor per sweep point mirroring the demoted `compute_BayesFactor` math
(`BF = (p_post/(1−p_post)) / (p_prior/(1−p_prior))`, `logBF = log(BF)`) but with NO `1e-8`
probability floor (Memo §5 / T-7-06). `post_draws` is a vector of per-sweep-point Δρ posterior
draw vectors; `prior_draws` is the fixed Δρ prior sample. Because the tail probabilities are not
clamped, `max|Δ logBF|` measured against THIS baseline is free of the KDE clamp/tail artifact
(STATE Phase-5-iter2) — it reflects the amortized estimator, not the floor. Gate-side only
(the amortized read path never touches KDE). May return `±Inf` for a perfectly one-sided
posterior; that is the honest, un-floored value the gate is meant to see.
"""
function kde_log_bf_unclamped(post_draws, prior_draws; threshold::Real = 0.0)
    p_prior    = _p_gt_threshold_unclamped(prior_draws; threshold = threshold)
    prior_odds = p_prior / (1 - p_prior)
    logbf = Vector{Float64}(undef, length(post_draws))
    for i in eachindex(post_draws)
        p_post         = _p_gt_threshold_unclamped(post_draws[i]; threshold = threshold)
        posterior_odds = p_post / (1 - p_post)
        logbf[i]       = log(posterior_odds / prior_odds)
    end
    return logbf
end
