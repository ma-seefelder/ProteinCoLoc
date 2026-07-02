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

# spike/validation/ood.jl --- Phase-5 OOD / misspecification flag (OOD-01, OOD-02).
#
# The scientific-honesty capstone: an amortized out-of-distribution flag for the FROZEN
# Phase-4 net, built as TWO channels fused by OR (D-05) and validated as a controlled
# ROC experiment over a misspecification grid (D-03) with a summary-orthogonal negative
# control whose blind spot is MEASURED and NAMED (D-04), not hidden.
#
# CENTRAL INVARIANTS (honesty-critical -- see 05-CONTEXT D-03..06):
#   • TRAIN-ONLY FREEZE (SC5 / D-06): every OOD null (the Mahalanobis μ_S/Σ_S density
#     model, and the posterior-predictive reference cloud) is fit on the TRAIN split
#     ONLY. `fit_ood_nulls(Ztrain)` never sees eval/misspec data; misspecified test
#     images are FULLY EXTERNAL and only ever SCORED, never fit against.
#   • CONTINUOUS-ROWS-ONLY MAHALANOBIS (Pitfall 2 / RESEARCH §OOD): the 128-dim :min
#     summary has 64 continuous correlation rows (1:64) + 64 binary MASK rows (65:128).
#     The near-constant mask rows make Σ singular, so Mahalanobis fits and scores the
#     CONTINUOUS rows ONLY (via the loader's authoritative Loader._row_partition), with
#     a ridge (1e-6) on the covariance for numerical stability.
#   • OR-FUSION (D-05): report each channel's own ROC/AUC AND an OR-combined flag that
#     fires iff EITHER channel exceeds its pre-registered ID-quantile threshold.
#   • PRE-REGISTERED OPERATING POINT (D-06): the gate is the ID-quantile threshold
#     (OOD_ID_QUANTILE, committed in consts.jl BEFORE seeing misspec data). Youden-J is
#     a labeled POST-HOC reference only -- NEVER the gate.
#   • NAMED BLIND SPOT (D-04): the negative control uses correlation-preserving
#     transforms that are empirically KS-invariant on the 8×8 summary (statistic <
#     OOD_KS_EPS), so the flag stays quiet there BY CONSTRUCTION -- the structural blind
#     spot of the fixed summary is measured, paired with the SBC "under the simulator"
#     caveat (SBC-04).
#   • NO NEW PACKAGE (Pitfall 7): the ROC/AUC is hand-rolled (LinearAlgebra/Statistics);
#     no ROC/normalizing-flow dependency enters the pinned env. The optional
#     NormalisingFlow density escalation (D-05) is out of scope unless the reported grid
#     shows Mahalanobis AUC < OOD_AUC_MIN -- documented as a follow-up, not added here.
#   • CPU-ONLY (D-10): use_gpu = false on every NeuralEstimators call.
#   • READ-ONLY src/ (CLAUDE.md): reached only transitively through the harness chain.
#
# Flat top-level functions (sibling style of sbc.jl/bf.jl). Guarded includes keep the
# file loadable standalone AND idempotent under runtests.jl.

using LinearAlgebra      # cholesky, Symmetric, UniformScaling I (Mahalanobis)
using Statistics         # mean, cov, quantile
using StatsBase          # reconstruct (inverse ZScoreTransform; via the harness too)
using HypothesisTests    # ApproximateTwoSampleKSTest (D-04 summary-invariance verification)

# The shared harness pulls consts.jl (OOD_* thresholds), the frozen model surface
# (load_frozen_model / posterior_for / standardize_summary), the simulator
# (sample_prior / simulate_pair), the summary contract (build_mci / patch_summary),
# encode_d01, val_rng, AND the Loader module (_row_partition) -- all transitively.
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "harness.jl"))

# ============================================================================
# Frozen-chain image → standardized summary (the external-image scorer)
# ============================================================================

"""
    frozen_summary(m, pair) -> Vector{Float32}

Put a 2-channel external image `pair::Vector{Matrix{Float64}}` into the estimator's
input space through the FROZEN read chain: `build_mci → patch_summary (UNCHANGED src
8×8 Pearson) → encode_d01 (128-dim) → standardize_summary(·, m.zt, :min)`. The SAME
path SBC/BF ride, so a misspecified image is scored byte-identically to training data
with NO re-fitting (T-05-01 / SC5). Returns the 128-dim standardized summary vector.
"""
frozen_summary(m, pair) =
    standardize_summary(encode_d01(patch_summary(build_mci(pair))), m.zt, :min)

"""
    id_summary(m, rng; imsize = SBC_IMSIZE) -> Vector{Float32}

One IN-DISTRIBUTION standardized summary: draw θ ~ π, `simulate_pair`, and push through
`frozen_summary`. Skips the posterior pass (the density channel needs only the summary),
so building a train/ID pool is cheap. Consumes `rng` sequentially (a disjoint stream).
"""
id_summary(m, rng; imsize = SBC_IMSIZE) =
    frozen_summary(m, simulate_pair(rng, sample_prior(rng); imsize = imsize))

# ============================================================================
# Channel 1 — summary-density Mahalanobis (fit on TRAIN continuous rows only)
# ============================================================================

"""
    fit_ood_nulls(Ztrain; variant = :min) -> NamedTuple

Fit the Mahalanobis summary-density null on the TRAIN split ONLY (SC5/D-06). `Ztrain`
is a `d×Ntrain` matrix of STANDARDIZED train summaries (e.g. `load_fold(...).Ztr` in the
reported run, or a train-distribution fixture pool in the fast gate). Using the loader's
authoritative `Loader._row_partition`, the fit touches the CONTINUOUS rows (1:64 for
:min) ONLY -- the binary mask rows (65:128) are near-constant and would make Σ singular
(Pitfall 2). Returns `(cont, μS, C)` where `C = cholesky(Symmetric(ΣS + 1e-6·I))` (a
ridge for numerical stability). The returned scorer NEVER touches the mask rows.

Fit inputs are TRAIN summaries only -- misspec/eval data is never passed here (T-05-07).
"""
function fit_ood_nulls(Ztrain::AbstractMatrix; variant::Symbol = :min)
    cont, _ = Loader._row_partition(variant, size(Ztrain, 1))
    Zc = Float64.(Ztrain[cont, :])
    μS = vec(mean(Zc; dims = 2))
    ΣS = cov(Zc; dims = 2)
    C  = cholesky(Symmetric(ΣS + 1e-6 * I))     # ridge 1e-6 (Pitfall 2 / RESEARCH)
    return (cont = cont, μS = μS, C = C)
end

"""
    maha_score(nulls, z) -> Float64

Squared Mahalanobis distance of a standardized summary `z` to the train density null,
on the CONTINUOUS rows only: `r = z[cont] − μS; ‖C.L \\ r‖²`. Non-negative and finite;
an in-distribution `z` scores LOW, a summary shifted off the train manifold scores HIGH.
The mask rows (65:128) are never read (Pitfall 2).
"""
function maha_score(nulls, z::AbstractVector)
    r = Float64.(z[nulls.cont]) .- nulls.μS
    return sum(abs2, nulls.C.L \ r)
end

# ============================================================================
# Channel 2 — posterior-predictive mismatch
# ============================================================================

# Reconstruct a simulate_pair-valid θ NamedTuple from a physical-θ posterior-mean
# vector, clamping into the prior-valid ranges (mirrors simulate_pair's entry guard:
# ρ∈[-1,1], spillover∈[0,1], autofluorescence≥0, label_efficiency∈[0,1], noise≥0).
_theta_tuple(v) = (
    ρ_true           = clamp(v[1], -1.0, 1.0),
    spillover        = clamp(v[2],  0.0, 1.0),
    autofluorescence = max(v[3], 0.0),
    label_efficiency = clamp(v[4],  0.0, 1.0),
    shift_dx         = v[5],
    shift_dy         = v[6],
    noise            = max(v[7], 0.0),
)

"""
    pp_mismatch_score(m, Z_obs; reps = OOD_PP_REPS, rng = val_rng(),
                      imsize = SBC_IMSIZE, N = 200) -> Float64

Posterior-predictive mismatch (D-05): infer θ̂ = posterior mean from `Z_obs` (one
`posterior_for` pass, un-standardized with the frozen `m.θzt`), re-simulate `reps`
summaries from θ̂ through `simulate_pair → frozen_summary`, and return the Mahalanobis
distance of `Z_obs`'s continuous rows to that posterior-predictive summary cloud (a
LOCAL mean/cov fit on the cloud, ridge 1e-6). Grows when the inferred θ̂ regenerates
summaries INCONSISTENT with the observed summary -- reaching some summary-orthogonal
cases the density channel misses. CPU-only (`use_gpu = false`). The re-simulation cloud
is a POSTERIOR-PREDICTIVE reference for THIS observation, not a train null, so it is not
subject to the train-only-freeze (it depends only on θ̂ and the frozen simulator).
"""
function pp_mismatch_score(m, Z_obs::AbstractVector; reps::Integer = OOD_PP_REPS,
                           rng = val_rng(), imsize = SBC_IMSIZE, N::Integer = 200)
    cont, _ = Loader._row_partition(:min, length(Z_obs))

    # Infer θ̂ (posterior mean, physical units) in ONE amortized pass.
    draws = posterior_for(m.estimator, Z_obs; N = N, use_gpu = false)
    orig  = StatsBase.reconstruct(m.θzt, draws)
    θhat  = _theta_tuple(vec(mean(orig; dims = 2)))

    # Posterior-predictive summary cloud (continuous rows) from θ̂.
    cloud = Matrix{Float64}(undef, length(cont), reps)
    for r in 1:reps
        Zr = frozen_summary(m, simulate_pair(rng, θhat; imsize = imsize))
        cloud[:, r] = Float64.(Zr[cont])
    end

    μc = vec(mean(cloud; dims = 2))
    Σc = cov(cloud; dims = 2)
    Cc = cholesky(Symmetric(Σc + 1e-6 * I))
    rr = Float64.(Z_obs[cont]) .- μc
    return sum(abs2, Cc.L \ rr)
end

# ============================================================================
# ROC / AUC (hand-rolled — NO new package, Pitfall 7) + operating points
# ============================================================================

"""
    roc_auc(neg_scores, pos_scores) -> (fpr, tpr, auc)

Hand-rolled ROC curve + AUC for an OOD detector where HIGHER score ⇒ more OOD.
`neg_scores` are in-distribution (label 0), `pos_scores` misspecified (label 1). The
curve `(fpr, tpr)` is swept over descending unique thresholds (monotone, endpoints
(0,0)→(1,1)); the AUC is the tie-aware Mann–Whitney U statistic (exact, robust to
threshold degeneracies) `= P(score_pos > score_neg) + ½·P(=)`. `auc ≈ 1.0` for a
perfectly separable input, `≈ 0.5` for identical distributions. No dependency (~15 lines).
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

The PRE-REGISTERED operating point (D-06): the `q`-quantile of the IN-DISTRIBUTION
scores (default `OOD_ID_QUANTILE`, ~5% ID false-positive rate). Committed in consts.jl
BEFORE seeing misspec data -- this, NOT Youden-J, is the gate.
"""
id_threshold(id_scores::AbstractVector; q::Real = OOD_ID_QUANTILE) =
    quantile(id_scores, q)

"""
    youden_j(fpr, tpr) -> NamedTuple

The POST-HOC Youden-J operating point (`max TPR − FPR`) along a ROC curve. Returned as
`(j, fpr, tpr, idx)` and shown as a LABELED REFERENCE on the ROC only -- it is NEVER the
gate (D-06). The gate is always the pre-registered `id_threshold`.
"""
function youden_j(fpr::AbstractVector, tpr::AbstractVector)
    j    = tpr .- fpr
    idx  = argmax(j)
    return (j = j[idx], fpr = fpr[idx], tpr = tpr[idx], idx = idx)
end

# ============================================================================
# OR-fusion flag (D-05)
# ============================================================================

"""
    ood_flag(nulls, m, Z; maha_thr, pp_thr, rng = val_rng(),
             reps = OOD_PP_REPS, imsize = SBC_IMSIZE) -> Bool

The OR-combined OOD flag (D-05): fires iff EITHER the summary-density channel
(`maha_score > maha_thr`) OR the posterior-predictive channel
(`pp_mismatch_score > pp_thr`) exceeds its PRE-REGISTERED ID-quantile threshold. Both
thresholds come from `id_threshold(id_scores; q = OOD_ID_QUANTILE)` fit on train/ID data
BEFORE misspec data (D-06). CPU-only.
"""
function ood_flag(nulls, m, Z::AbstractVector; maha_thr::Real, pp_thr::Real,
                  rng = val_rng(), reps::Integer = OOD_PP_REPS, imsize = SBC_IMSIZE)
    (maha_score(nulls, Z) > maha_thr) && return true
    return pp_mismatch_score(m, Z; reps = reps, rng = rng, imsize = imsize) > pp_thr
end
