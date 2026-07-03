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
# src/amortized/ood.jl --- the amortized OOD / misspecification FLAG (PROD-01).
#
# Promoted (grid-generalized, zero new deps — hand-rolled Mahalanobis + ROC via
# LinearAlgebra/Statistics) from the proven spike:
#   spike/validation/ood.jl:112-132  fit_ood_nulls / maha_score       (density channel)
#   spike/validation/ood.jl:145-159  _finite_or / _theta_tuple        (PP finite-guard, iter1)
#   spike/validation/ood.jl:182-203  pp_mismatch_score                (posterior-predictive)
#   spike/validation/ood.jl:233-303  noise_features / fit_noise_null / noise_score (noise chan.)
#   spike/validation/ood.jl:319-357  roc_auc / id_threshold / youden_j
#   spike/validation/ood.jl:375-380  ood_flag (OR-fusion)
#
# THREE-CHANNEL OR-FUSED detector (D-05):
#   • density  — Mahalanobis on the CONTINUOUS summary rows only (the binary mask rows are
#                near-constant and would make Σ singular; a 1e-6·I ridge stabilizes it).
#   • noise    — Mahalanobis on scale/rotation/permutation-INVARIANT image-noise features
#                (iter2), closing the detector-noise blind spot of the correlation-only summary.
#   • pp       — posterior-predictive per-feature discrepancy: infer θ̂, re-simulate a summary
#                cloud, score the observed summary against it.
#
# MEMO §5 HARDENING (T-7-04): the reported OR-fusion now includes the posterior-predictive
# channel (`with_pp = true` by default — the spike ran `with_pp = false`). The iter1 finite-guard
# (`_finite_or`/`_theta_tuple`) maps a non-finite posterior θ̂ (which a STRONG misspecification
# can drive to NaN/Inf) to an in-range, prior-valid fallback, so PP re-simulation stays
# COMPUTABLE on misspecified inputs instead of crashing (the misspecified summary then sits far
# from the typical-θ cloud → a HIGH, correct OOD score, not a crash).
#
# GRID-AGNOSTIC: the continuous/mask row split is delegated to `_summary_row_partition`
# (summary.jl), so a new shipped grid flows through here with no change.
#
# SCOPE: this file is the OOD READ SURFACE (the channels + ROC helpers + fused `ood_verdict`).
# The controlled misspecification-grid ROC EXPERIMENT (misspec_* families, negative controls,
# `ood_roc_over_grid`) is the per-grid D-05 SHIP-GATE's job and is added by a later Phase-7 plan
# (it needs the promoted simulator + ImageFiltering); it is deliberately NOT promoted here.
#############################################################################################

import LinearAlgebra: cholesky, Symmetric, I
import Statistics: mean, median, var, std, cor, quantile
import StatsBase

# --- OOD read-surface defaults (the per-grid D-05 gate re-expresses these fresh) ------------
const OOD_PP_REPS     = 50      # posterior-predictive re-simulations per PP score (spike value)
const OOD_ID_QUANTILE = 0.95    # pre-registered ID operating quantile (~5% ID false-positive)

# ============================================================================
# Channel 1 — summary-density Mahalanobis (fit on TRAIN continuous rows only)
# ============================================================================

"""
    fit_ood_nulls(Ztrain; variant = :min) -> NamedTuple

Fit the Mahalanobis summary-density null on the TRAIN split ONLY (frozen-stats discipline).
`Ztrain` is a `d×Ntrain` matrix of STANDARDIZED train summaries. Using the grid-general
`_summary_row_partition` (summary.jl), the fit touches the CONTINUOUS rows ONLY — the binary
mask rows are near-constant and would make Σ singular. Returns `(cont, μS, C)` with
`C = cholesky(Symmetric(ΣS + 1e-6·I))` (a numerical ridge). The returned scorer never touches
the mask rows. Fit inputs are TRAIN summaries only — misspec/eval data is never passed here.
"""
function fit_ood_nulls(Ztrain::AbstractMatrix; variant::Symbol = :min)
    cont, _ = _summary_row_partition(variant, size(Ztrain, 1))
    Zc = Float64.(Ztrain[cont, :])
    μS = vec(mean(Zc; dims = 2))
    ΣS = cov_cols(Zc)
    C  = cholesky(Symmetric(ΣS + 1e-6 * I))     # ridge 1e-6 (Σ-singularity guard)
    return (cont = cont, μS = μS, C = C)
end

# Column-wise covariance (Statistics.cov is not imported into core; hand-roll over columns so
# the OOD channel keeps its LinearAlgebra/Statistics-only footprint with no new import).
function cov_cols(Zc::AbstractMatrix)
    n = size(Zc, 2)
    μ = vec(mean(Zc; dims = 2))
    X = Zc .- μ
    return (X * X') ./ (n - 1)
end

"""
    maha_score(nulls, z) -> Float64

Squared Mahalanobis distance of a standardized summary `z` to the train density null, on the
CONTINUOUS rows only: `r = z[cont] − μS; ‖C.L \\ r‖²`. Non-negative and finite; an
in-distribution `z` scores LOW, a summary off the train manifold scores HIGH. The mask rows are
never read.
"""
function maha_score(nulls, z::AbstractVector)
    r = Float64.(z[nulls.cont]) .- nulls.μS
    return sum(abs2, nulls.C.L \ r)
end

# ============================================================================
# Channel 2 — posterior-predictive mismatch (with the iter1 finite-guard)
# ============================================================================

# Finite-guard (iter1): a strongly MISSPECIFIED input can drive the flow's posterior mean to a
# NON-FINITE θ̂ (NaN/Inf). Since `clamp(NaN, …) == NaN` and `max(NaN, 0) == NaN`, an unguarded
# NaN would propagate into `simulate_pair` and CRASH the PP channel. Mapping a non-finite
# component to an in-range fallback keeps PP COMPUTABLE: the re-simulated cloud is a valid
# (typical-θ) reference, and the misspecified `Z_obs` sits far from it, so the PP score comes out
# HIGH — the correct OOD signal, not a crash (Memo §5 / T-7-04).
_finite_or(x, default) = isfinite(x) ? float(x) : default

# Reconstruct a simulate_pair-valid θ NamedTuple from a physical-θ posterior-mean vector, first
# mapping any non-finite component to an in-range fallback, then clamping into the prior-valid
# ranges (mirrors simulate_pair's entry guard: ρ∈[-1,1], spillover∈[0,1], autofluorescence≥0,
# label_efficiency∈[0,1], noise≥0).
_theta_tuple(v) = (
    ρ_true           = clamp(_finite_or(v[1], 0.0), -1.0, 1.0),
    spillover        = clamp(_finite_or(v[2], 0.0),  0.0, 1.0),
    autofluorescence = max(_finite_or(v[3], 0.0), 0.0),
    label_efficiency = clamp(_finite_or(v[4], 1.0),  0.0, 1.0),
    shift_dx         = _finite_or(v[5], 0.0),
    shift_dy         = _finite_or(v[6], 0.0),
    noise            = max(_finite_or(v[7], 0.0), 0.0),
)

"""
    frozen_summary(m, pair) -> Vector{Float32}

Put a 2-channel external image `pair::Vector{Matrix{Float64}}` into the estimator's input space
through the FROZEN read chain: `build_mci → patch_summary(·, m.grid) → encode_d01 →
standardize_summary(·, m.zt, :min)`. The SAME path the read surface rides, so a misspecified
image is scored byte-identically to training data with NO re-fitting. `m` is a NamedTuple
`(estimator, θzt, zt, grid)` (the frozen model handle).
"""
frozen_summary(m, pair) =
    standardize_summary(encode_d01(patch_summary(build_mci(pair), m.grid)), m.zt, :min)

"""
    pp_mismatch_score(m, Z_obs; reps = OOD_PP_REPS, rng, imsize, N = 200) -> Float64

Posterior-predictive mismatch (Memo §5, re-enabled): infer θ̂ = posterior mean from `Z_obs` (one
`posterior_for` pass, un-standardized with the frozen `m.θzt`, then finite-guarded via
`_theta_tuple`), re-simulate `reps` summaries from θ̂ through `simulate_pair → frozen_summary`,
and return the mean squared PER-FEATURE z-score of `Z_obs`'s continuous rows against that
posterior-predictive summary cloud. A DIAGONAL (per-feature) discrepancy is used deliberately
(the PP cloud has only `reps` samples; a full continuous-row covariance would be rank-deficient
for `reps < G²`). CPU-only. `m` is the frozen model handle `(estimator, θzt, zt, grid)`.
"""
function pp_mismatch_score(m, Z_obs::AbstractVector; reps::Integer = OOD_PP_REPS,
                           rng, imsize, N::Integer = 200)
    cont, _ = _summary_row_partition(:min, length(Z_obs))

    # Infer θ̂ (posterior mean, physical units) in ONE amortized pass, then FINITE-GUARD it.
    draws = posterior_for(m.estimator, Z_obs; N = N, use_gpu = false)
    orig  = StatsBase.reconstruct(m.θzt, draws)
    θhat  = _theta_tuple(vec(mean(orig; dims = 2)))

    # Posterior-predictive summary cloud (continuous rows) from θ̂.
    cloud = Matrix{Float64}(undef, length(cont), reps)
    for r in 1:reps
        Zr = frozen_summary(m, simulate_pair(rng, θhat; imsize = imsize))
        cloud[:, r] = Float64.(Zr[cont])
    end

    # Diagonal (per-feature) discrepancy — well-conditioned for reps < dim (ridge on σ²).
    μc = vec(mean(cloud; dims = 2))
    σ2 = vec(var(cloud; dims = 2)) .+ 1e-6
    z  = (Float64.(Z_obs[cont]) .- μc) .^ 2 ./ σ2
    return mean(z)
end

# ============================================================================
# Channel 3 — image-noise Mahalanobis (narrows the detector-noise blind spot)
# ============================================================================
#
# Auxiliary noise-sensitive features computed DIRECTLY from the images (NOT from the frozen
# summary, which is the UNCHANGED NPE/NRE input and must not change), OR-fused into the detector.
# Every feature is a RATIO / normalized moment / correlation, hence INVARIANT to the affine
# (a·x+b, a>0) and rotation/permutation negative-control transforms, so the noise channel stays
# quiet on the summary-orthogonal blind-spot transforms BY CONSTRUCTION.

# Median absolute deviation (robust scale). Hand-rolled (no new dependency).
_mad(v) = (mv = median(v); median(abs.(v .- mv)))

"""
    noise_features(pair) -> Vector{Float64}

Per-channel scale/rotation/permutation-INVARIANT noise-sensitivity features for a 2-channel
image `pair::Vector{Matrix{Float64}}` (10 features = 5 per channel): high-frequency energy ratio
`var(∇²x)/var(x)`, robust HF scale ratio `MAD(∇²x)/MAD(x)`, pixel-outlier fraction, HF excess
kurtosis, and lag-1 spatial autocorrelation. `∇²x` is the interior 5-point finite-difference
Laplacian (no FFT/dependency). Non-finite guards map degenerate values to 0.
"""
function noise_features(pair)
    feats = Float64[]
    for ch in pair
        X = Float64.(ch)
        H, W = size(X)
        lap = @views 4.0 .* X[2:H-1, 2:W-1] .- X[1:H-2, 2:W-1] .- X[3:H, 2:W-1] .-
                     X[2:H-1, 1:W-2] .- X[2:H-1, 3:W]
        lv  = vec(lap)
        xv  = vec(X)
        madl = _mad(lv); madx = _mad(xv)
        f1 = var(lv) / (var(xv) + 1e-12)
        f2 = madl / (madx + 1e-12)
        f3 = mean(abs.(lv) .> 6.0 * madl + 1e-12)
        f4 = StatsBase.kurtosis(lv)                                   # excess kurtosis
        ah = cor(vec(@view X[1:H-1, :]), vec(@view X[2:H, :]))
        av = cor(vec(@view X[:, 1:W-1]), vec(@view X[:, 2:W]))
        f5 = 0.5 * (ah + av)
        for f in (f1, f2, f3, f4, f5)
            push!(feats, isfinite(f) ? f : 0.0)
        end
    end
    return feats
end

"""
    fit_noise_null(F) -> NamedTuple

Fit the image-noise Mahalanobis null on the TRAIN-ONLY feature matrix `F` (`nfeat×Ntrain`, one
column per fit image). Features are z-standardized (a constant feature guarded to scale 1) so
the covariance ridge is meaningful across heterogeneous feature scales; returns `(μ, σ, C)` with
`C = cholesky(Symmetric(Σ_std + 1e-6·I))`.
"""
function fit_noise_null(F::AbstractMatrix)
    μ = vec(mean(F; dims = 2))
    σ = vec(std(F; dims = 2)); @inbounds for i in eachindex(σ); (σ[i] < 1e-8) && (σ[i] = 1.0); end
    Fs = (F .- μ) ./ σ
    Σ  = cov_cols(Fs)
    C  = cholesky(Symmetric(Σ + 1e-6 * I))
    return (μ = μ, σ = σ, C = C)
end

"""
    noise_score(nn, pair) -> Float64

Squared Mahalanobis distance of an image pair's `noise_features` to the train noise null `nn`
(`fit_noise_null`), in the standardized feature space. In-distribution (smooth) images score
LOW; salt-and-pepper / heavy-tailed detector-noise mismatch scores HIGH.
"""
function noise_score(nn, pair)
    f = (noise_features(pair) .- nn.μ) ./ nn.σ
    return sum(abs2, nn.C.L \ f)
end

# ============================================================================
# ROC / AUC (hand-rolled — NO new package) + operating points
# ============================================================================

"""
    roc_auc(neg_scores, pos_scores) -> (fpr, tpr, auc)

Hand-rolled ROC curve + AUC for an OOD detector where HIGHER score ⇒ more OOD. `neg_scores` are
in-distribution (label 0), `pos_scores` misspecified (label 1). The curve `(fpr, tpr)` is swept
over descending unique thresholds; the AUC is the tie-aware Mann–Whitney U statistic (exact,
robust to threshold degeneracies) `= P(score_pos > score_neg) + ½·P(=)`. `auc ≈ 1.0` for a
perfectly separable input, `≈ 0.5` for identical distributions. No dependency.
"""
function roc_auc(neg_scores::AbstractVector, pos_scores::AbstractVector)
    P = length(pos_scores); N = length(neg_scores)
    (P == 0 || N == 0) && return (Float64[0.0, 1.0], Float64[0.0, 1.0], NaN)

    thr = sort(unique(vcat(collect(neg_scores), collect(pos_scores))); rev = true)
    fpr = Float64[0.0]; tpr = Float64[0.0]
    for t in thr
        push!(tpr, count(>=(t), pos_scores) / P)
        push!(fpr, count(>=(t), neg_scores) / N)
    end

    u = 0.0
    for p in pos_scores, n in neg_scores
        u += p > n ? 1.0 : (p == n ? 0.5 : 0.0)
    end
    return (fpr, tpr, u / (P * N))
end

"""
    id_threshold(id_scores; q = OOD_ID_QUANTILE) -> Float64

The PRE-REGISTERED operating point: the `q`-quantile of the IN-DISTRIBUTION scores (default
`OOD_ID_QUANTILE`, ~5% ID false-positive rate). Committed BEFORE seeing misspec data — this, NOT
Youden-J, is the gate.
"""
id_threshold(id_scores::AbstractVector; q::Real = OOD_ID_QUANTILE) = quantile(id_scores, q)

"""
    youden_j(fpr, tpr) -> NamedTuple

The POST-HOC Youden-J operating point (`max TPR − FPR`) along a ROC curve, returned as
`(j, fpr, tpr, idx)`. A LABELED REFERENCE only — NEVER the gate (the gate is `id_threshold`).
"""
function youden_j(fpr::AbstractVector, tpr::AbstractVector)
    j    = tpr .- fpr
    idx  = argmax(j)
    return (j = j[idx], fpr = fpr[idx], tpr = tpr[idx], idx = idx)
end

# ============================================================================
# OR-fusion flag + productionized fused verdict (D-05 / Memo §5)
# ============================================================================

"""
    ood_flag(nulls, m, Z; maha_thr, pp_thr, rng, reps = OOD_PP_REPS, imsize, N = 200) -> Bool

The OR-combined OOD flag (D-05): fires iff EITHER the summary-density channel
(`maha_score > maha_thr`) OR the posterior-predictive channel (`pp_mismatch_score > pp_thr`)
exceeds its PRE-REGISTERED ID-quantile threshold. Both thresholds come from `id_threshold` fit on
train/ID data BEFORE misspec data. `reps`/`N`/`imsize` MUST match the values `pp_thr` was
computed with. CPU-only.
"""
function ood_flag(nulls, m, Z::AbstractVector; maha_thr::Real, pp_thr::Real,
                  rng, reps::Integer = OOD_PP_REPS, imsize, N::Integer = 200)
    (maha_score(nulls, Z) > maha_thr) && return true
    return pp_mismatch_score(m, Z; reps = reps, rng = rng, imsize = imsize, N = N) > pp_thr
end

# robust-z of a raw channel score against its train reference (median / MAD-scale), the
# continuous analogue of "fire iff the channel exceeds its ID-quantile point" (D-05 OR-fusion).
function _robust_z(zref, name::Symbol, s::Real)
    haskey(zref, name) || return float(s)
    r = zref[name]
    return (float(s) - r.med) / r.scl
end

"""
    ood_verdict(ood_nulls, Zs, Zc = nothing; pair = nothing, with_pp = true,
                rng = nothing, reps = OOD_PP_REPS, imsize = (256, 256)) -> OODVerdict

The productionized fused OOD verdict for a (sample, control) analysis. Fuses (OR, via robust-z
against the train reference `ood_nulls.zref`) every channel its inputs support:

  • density — `maha_score(ood_nulls.density, Zs)` (and `Zc` when given; the pair max),
  • noise   — `noise_score(ood_nulls.noise, pair)` when the image `pair` and a fitted
              `ood_nulls.noise` null are available,
  • pp      — `pp_mismatch_score(ood_nulls.model, Zs; …)` when `with_pp` and a frozen
              `ood_nulls.model` handle are available (Memo §5: `with_pp = true` by DEFAULT — the
              re-enabled posterior-predictive channel, finite-guarded).

Returns an `OODVerdict(score, flag, per_channel)`: `score` is the fused (max robust-z) statistic,
`flag = score > ood_nulls.thr` (the pre-registered fused operating point), and `per_channel`
carries each active channel's RAW score. `ood_nulls` is a NamedTuple carrying at least
`:density`; `:noise`, `:model`, `:zref`, `:thr` are populated by the per-grid D-05 ship-gate.
CPU-only.
"""
function ood_verdict(ood_nulls, Zs, Zc = nothing; pair = nothing, with_pp::Bool = true,
                     rng = nothing, reps::Integer = OOD_PP_REPS, imsize = (256, 256))
    channels = Tuple{Symbol,Float64}[]

    # Channel 1 — summary-density Mahalanobis (always available from the standardized summary).
    dens = maha_score(ood_nulls.density, Zs)
    if Zc !== nothing
        dens = max(dens, maha_score(ood_nulls.density, Zc))
    end
    push!(channels, (:density, dens))

    # Channel 3 — image-noise Mahalanobis (needs the image pair + a fitted noise null).
    if pair !== nothing && haskey(ood_nulls, :noise) && ood_nulls.noise !== nothing
        push!(channels, (:noise, noise_score(ood_nulls.noise, pair)))
    end

    # Channel 2 — posterior-predictive mismatch (Memo §5: with_pp = true by default).
    if with_pp && haskey(ood_nulls, :model) && ood_nulls.model !== nothing
        push!(channels, (:pp, pp_mismatch_score(ood_nulls.model, Zs;
                                                reps = reps, rng = rng, imsize = imsize)))
    end

    zref  = haskey(ood_nulls, :zref) ? ood_nulls.zref : nothing
    fused = -Inf
    per   = NamedTuple()
    for (name, s) in channels
        z     = zref === nothing ? float(s) : _robust_z(zref, name, s)
        fused = max(fused, z)
        per   = merge(per, NamedTuple{(name,)}((s,)))
    end

    thr  = haskey(ood_nulls, :thr) ? ood_nulls.thr : nothing
    flag = thr === nothing ? false : fused > thr
    return OODVerdict(fused, flag, per)
end
