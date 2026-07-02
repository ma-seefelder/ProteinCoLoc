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

# Finite-guard: a strongly MISSPECIFIED input can drive the flow's posterior mean to a
# NON-FINITE θ̂ (NaN/Inf). Since `clamp(NaN, …) == NaN` and `max(NaN, 0) == NaN`, an
# unguarded NaN would propagate into `simulate_pair` and CRASH the PP channel (STATE
# Phase-5 blocker: "PP channel non-viable on frozen net"). Replacing a non-finite
# component with an in-range fallback keeps the PP channel COMPUTABLE: the re-simulated
# cloud is then a valid (typical-θ) reference, and the misspecified `Z_obs` sits far from
# it, so the PP score comes out HIGH — the correct OOD signal, not a crash.
_finite_or(x, default) = isfinite(x) ? float(x) : default

# Reconstruct a simulate_pair-valid θ NamedTuple from a physical-θ posterior-mean
# vector, first mapping any non-finite component to an in-range fallback, then clamping
# into the prior-valid ranges (mirrors simulate_pair's entry guard: ρ∈[-1,1],
# spillover∈[0,1], autofluorescence≥0, label_efficiency∈[0,1], noise≥0).
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
    pp_mismatch_score(m, Z_obs; reps = OOD_PP_REPS, rng = val_rng(),
                      imsize = SBC_IMSIZE, N = 200) -> Float64

Posterior-predictive mismatch (D-05): infer θ̂ = posterior mean from `Z_obs` (one
`posterior_for` pass, un-standardized with the frozen `m.θzt`), re-simulate `reps`
summaries from θ̂ through `simulate_pair → frozen_summary`, and return the mean squared
PER-FEATURE z-score of `Z_obs`'s continuous rows against that posterior-predictive
summary cloud: `mean_f ((z_f − μc_f) / σc_f)²`. Grows when the inferred θ̂ regenerates
summaries INCONSISTENT with the observed summary -- reaching some summary-orthogonal
cases the density channel misses.

A DIAGONAL (per-feature) discrepancy is used deliberately, NOT a full-covariance
Mahalanobis (RESEARCH §PP "per-feature z-score" option): the PP cloud has only `reps`
samples, and with `reps < 64` (the continuous-row dimension -- true even at the reported
OOD_PP_REPS=50) a full 64×64 cloud covariance is RANK-DEFICIENT and ridge-dominated,
producing meaningless ~1e7 scores. The per-feature variance is well-conditioned at any
`reps ≥ 2`. CPU-only (`use_gpu = false`). The PP cloud is a posterior-predictive
reference for THIS observation, not a train null, so it is not subject to the
train-only-freeze (it depends only on θ̂ and the frozen simulator).
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

    # Diagonal (per-feature) discrepancy -- well-conditioned for reps < dim (ridge on σ²).
    μc = vec(mean(cloud; dims = 2))
    σ2 = vec(var(cloud; dims = 2)) .+ 1e-6
    z  = (Float64.(Z_obs[cont]) .- μc) .^ 2 ./ σ2
    return mean(z)
end

# ============================================================================
# Channel 3 — image-noise Mahalanobis (narrows the detector-noise blind spot)
# ============================================================================
#
# WHY (05-03 named the detector-NOISE mismatch as a MEASURED blind spot of the fixed 8×8
# patch-Pearson summary: salt-and-pepper + heavy-tailed corruption DECORRELATES patches
# toward a central summary, so the correlation-only summary — and every channel built on
# it (density, PP) — is structurally near-blind to it, AUC≈0). This channel adds AUXILIARY
# noise-sensitive features computed DIRECTLY FROM THE IMAGES (NOT from the frozen 8×8
# summary, which is the UNCHANGED NPE/NRE inference input and MUST NOT change) and OR-fuses
# them into the OOD detector.
#
# FIRST-PRINCIPLES DESIGN (SC5 — designed from the TRAINING distribution, NOT from the
# held-out noise test set): the simulator (forward.jl) produces SMOOTH images — STRUCT_σ=6
# Gaussian fields ≫ the σ_psf=1.3 PSF — with only bounded Poisson(GAIN=50) + tiny
# Gaussian(READ_NOISE·θ.noise, θ.noise∈[0,1]) noise. So an in-distribution image has LOW
# high-frequency energy, near-Gaussian local statistics, few pixel-scale outliers, and HIGH
# nearest-neighbour autocorrelation. The misspec_noise family injects pixel-scale
# salt-and-pepper spikes + heavy-tailed (cubic) additive noise — exactly the high-frequency,
# heavy-tailed, low-autocorrelation signature the training distribution lacks.
#
# INVARIANCE (honesty-critical): every feature is a RATIO or a normalized moment or a
# correlation, so it is INVARIANT to the affine (a·x+b, a>0) and rotation/permutation
# negative-control transforms (D-04) — the noise channel therefore stays quiet on the
# summary-orthogonal blind-spot transforms BY CONSTRUCTION, exactly as the density channel
# does, so the OR-fused detector does not false-fire on them (verified on DEV).

# Median absolute deviation (robust scale). Hand-rolled (no new dependency).
_mad(v) = (mv = median(v); median(abs.(v .- mv)))

"""
    noise_features(pair) -> Vector{Float64}

Per-channel scale/rotation/permutation-INVARIANT noise-sensitivity features for a
2-channel image `pair::Vector{Matrix{Float64}}` (10 features = 5 per channel):

  f1 = high-frequency ENERGY ratio  var(∇²x) / var(x)                (smooth⇒low, noisy⇒high)
  f2 = robust HF scale ratio        MAD(∇²x) / MAD(x)                (robust twin of f1)
  f3 = pixel-outlier fraction       mean(|∇²x| > 6·MAD(∇²x))         (salt-and-pepper⇒high)
  f4 = HF excess kurtosis           kurtosis(∇²x)                    (heavy tails⇒high)
  f5 = lag-1 spatial autocorrelation ½(corr_h + corr_v) of x        (smooth⇒~1, noisy⇒low)

`∇²x` is the interior 5-point finite-difference Laplacian (no FFT/dependency). Every feature
is a ratio / normalized moment / correlation, hence invariant to `x → a·x + b` (a>0) and to
grid rotation/permutation. Non-finite guards map degenerate (constant-slice) values to 0.
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

Fit the image-noise Mahalanobis null on the TRAIN-ONLY feature matrix `F` (`nfeat×Ntrain`,
one column per fit image; SC5 — never misspec data). Features are first z-standardized
(per-feature mean/std, a constant feature guarded to scale 1) so the covariance ridge is
meaningful across the heterogeneous feature scales; returns `(μ, σ, C)` with
`C = cholesky(Symmetric(Σ_std + 1e-6·I))`.
"""
function fit_noise_null(F::AbstractMatrix)
    μ = vec(mean(F; dims = 2))
    σ = vec(std(F; dims = 2)); @inbounds for i in eachindex(σ); (σ[i] < 1e-8) && (σ[i] = 1.0); end
    Fs = (F .- μ) ./ σ
    Σ  = cov(Fs; dims = 2)
    C  = cholesky(Symmetric(Σ + 1e-6 * I))
    return (μ = μ, σ = σ, C = C)
end

"""
    noise_score(nn, pair) -> Float64

Squared Mahalanobis distance of an image pair's `noise_features` to the train noise null
`nn` (fit_noise_null), in the standardized feature space. In-distribution (smooth) images
score LOW; salt-and-pepper / heavy-tailed detector-noise mismatch scores HIGH.
"""
function noise_score(nn, pair)
    f = (noise_features(pair) .- nn.μ) ./ nn.σ
    return sum(abs2, nn.C.L \ f)
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
             reps = OOD_PP_REPS, imsize = SBC_IMSIZE, N = 200) -> Bool

The OR-combined OOD flag (D-05): fires iff EITHER the summary-density channel
(`maha_score > maha_thr`) OR the posterior-predictive channel
(`pp_mismatch_score > pp_thr`) exceeds its PRE-REGISTERED ID-quantile threshold. Both
thresholds come from `id_threshold(id_scores; q = OOD_ID_QUANTILE)` fit on train/ID data
BEFORE misspec data (D-06). `reps`/`N`/`imsize` MUST match the values the `pp_thr` was
computed with, so the PP channel is scored on the same footing as its threshold. CPU-only.
"""
function ood_flag(nulls, m, Z::AbstractVector; maha_thr::Real, pp_thr::Real,
                  rng = val_rng(), reps::Integer = OOD_PP_REPS, imsize = SBC_IMSIZE,
                  N::Integer = 200)
    (maha_score(nulls, Z) > maha_thr) && return true
    return pp_mismatch_score(m, Z; reps = reps, rng = rng, imsize = imsize, N = N) > pp_thr
end

# ============================================================================
# Misspecification grid — POSITIVE controls (D-03, four families)
# ============================================================================
#
# All four families produce FULLY EXTERNAL image pairs (SC5) that the shared-latent
# smooth-Gaussian-field simulator CANNOT produce, each swept over OOD_GRID_LEVELS
# magnitudes for a dose-response ROC. Every generator mirrors simulate_pair's entry
# guard (imsize dims ≥ 64, T-05-09) and returns finite, non-negative channels.

# Shared imsize≥64 guard (mirrors simulate_pair's entry validation; T-05-09 / CR-02).
function _guard_misspec_imsize(imsize)
    (imsize[1] ≥ 64 && imsize[2] ≥ 64) ||
        throw(ArgumentError("misspec imsize dims must be ≥ 64 for the 8×8 patch grid, " *
                            "got $imsize"))
    return nothing
end

# Sparse Gaussian-puncta field: `nspots` bright impulses smoothed at scale σ. The
# building block of the texture family — sharp, sparse spots the SMOOTH-field simulator
# structurally cannot generate (Phase-2 rejected puncta, D-15 — the canonical real-world
# OOD since real immunofluorescence is punctate).
function _puncta_field(rng, imsize, nspots::Integer, σ::Real)
    f = zeros(Float64, imsize...)
    H, W = imsize
    for _ in 1:nspots
        f[rand(rng, 1:H), rand(rng, 1:W)] += 0.5 + rand(rng)
    end
    return imfilter(f, Kernel.gaussian(σ))
end

"""
    misspec_texture(rng, θ; imsize = SBC_IMSIZE, level = OOD_GRID_LEVELS)
        -> Vector{Matrix{Float64}}

Family 1 (D-03) — TEXTURE-model mismatch. A NEW puncta/granular generator (a shared +
private Gaussian-spot mixture, NOT `simulate_pair`): higher `level` ⇒ more, sharper
spots ⇒ further from the smooth-field training manifold. θ.ρ_true drives the shared spot
weight (with `sign(ρ)` on channel 2, mirroring the simulator's anti-correlation path) so
the OOD image still carries a colocalization structure. The canonical real-world OOD
positive control.
"""
function misspec_texture(rng, θ; imsize = SBC_IMSIZE, level::Integer = OOD_GRID_LEVELS)
    _guard_misspec_imsize(imsize)
    ρ = clamp(θ.ρ_true, -1.0, 1.0)
    σ = max(0.8, 3.0 - 0.5 * level)              # sharper (more un-smooth) at higher level
    a = sqrt(abs(ρ)); b = sqrt(1.0 - abs(ρ))
    shared = _puncta_field(rng, imsize, 30 * level, σ)
    p1     = _puncta_field(rng, imsize, 20 * level, σ)
    p2     = _puncta_field(rng, imsize, 20 * level, σ)
    ch1 = BG_FLOOR .+ a .* shared .+ b .* p1
    ch2 = BG_FLOOR .+ sign(ρ) .* a .* shared .+ b .* p2
    return [Matrix{Float64}(max.(ch1, 0.0)), Matrix{Float64}(max.(ch2, 0.0))]
end

"""
    misspec_noise(rng, θ; imsize = SBC_IMSIZE, level = OOD_GRID_LEVELS)
        -> Vector{Matrix{Float64}}

Family 2 (D-03) — NOISE-model mismatch. Corrupts a `simulate_pair` base with
salt-and-pepper spikes plus heavy-tailed (cubic) additive noise BEYOND the simulator's
Poisson+Gaussian model; `level` scales both the spike fraction and the amplitude.
"""
function misspec_noise(rng, θ; imsize = SBC_IMSIZE, level::Integer = OOD_GRID_LEVELS)
    _guard_misspec_imsize(imsize)
    base  = simulate_pair(rng, θ; imsize = imsize)
    frac  = 0.02 * level
    scale = 2.0 * level
    out = map(base) do ch
        c  = copy(ch)
        mx = maximum(c)
        @inbounds for idx in eachindex(c)
            rand(rng) < frac && (c[idx] = rand(rng) < 0.5 ? 0.0 : mx * (1 + scale))
        end
        c .+= scale .* mx .* (rand(rng, size(c)...) .- 0.5) .^ 3   # heavy-tailed, non-Gaussian
        Matrix{Float64}(max.(c, 0.0))
    end
    return collect(out)
end

"""
    misspec_optics(rng, θ; imsize = SBC_IMSIZE, level = OOD_GRID_LEVELS)
        -> Vector{Matrix{Float64}}

Family 3 (D-03) — OPTICS/PSF mismatch. Convolves a `simulate_pair` base with a strongly
ANISOTROPIC/aberrated Gaussian PSF (σy ≫ σx) far beyond the fixed isotropic `σ_psf`;
`level` scales the anisotropy the simulator cannot produce.
"""
function misspec_optics(rng, θ; imsize = SBC_IMSIZE, level::Integer = OOD_GRID_LEVELS)
    _guard_misspec_imsize(imsize)
    base = simulate_pair(rng, θ; imsize = imsize)
    σy   = σ_psf + 2.0 * level                    # ≫ σ_psf — aberrated/out-of-focus
    σx   = σ_psf + 0.2 * level
    return [Matrix{Float64}(max.(imfilter(ch, Kernel.gaussian((σy, σx))), 0.0)) for ch in base]
end

"""
    misspec_background(rng, θ; imsize = SBC_IMSIZE, level = OOD_GRID_LEVELS)
        -> Vector{Matrix{Float64}}

Family 4 (D-03) — BACKGROUND/illumination mismatch. Applies a multiplicative
illumination gradient + radial vignette and an autofluorescence bleed BEYOND the
simulator's `Uniform(0,0.1)` offset; `level` scales the gradient, vignette and bleed.
"""
function misspec_background(rng, θ; imsize = SBC_IMSIZE, level::Integer = OOD_GRID_LEVELS)
    _guard_misspec_imsize(imsize)
    base = simulate_pair(rng, θ; imsize = imsize)
    H, W = imsize
    gy   = reshape(range(-1.0, 1.0; length = H), H, 1)   # H×1
    gx   = reshape(range(-1.0, 1.0; length = W), 1, W)   # 1×W
    grad  = 1.0 .+ (0.5 * level) .* gy .+ (0.3 * level) .* gx           # H×W
    vign  = 1.0 .- (0.2 * level) .* (gy .^ 2 .+ gx .^ 2)               # H×W
    bleed = 0.1 + 0.3 * level                                          # beyond U(0,0.1)
    return [Matrix{Float64}(max.(ch .* max.(grad, 0.0) .* max.(vign, 0.0) .+ bleed, 0.0))
            for ch in base]
end

# The four positive-control families, in fixed order (the reported ROC grid, D-03).
const OOD_FAMILIES = (texture    = misspec_texture,
                      noise      = misspec_noise,
                      optics     = misspec_optics,
                      background = misspec_background)

# ============================================================================
# Summary-orthogonal NEGATIVE control (D-04) — the named blind spot
# ============================================================================
#
# Correlation-preserving transforms of IN-DISTRIBUTION images. The 8×8 patch-Pearson
# summary is (exactly, for affine; as a value-multiset, for rotate/permute) UNCHANGED, so
# the flag stays quiet BY CONSTRUCTION — measuring and NAMING the structural blind spot of
# the fixed summary (paired with SBC-04's "under the simulator" caveat).

"""
    negctrl_affine(pair; a = 2.0, b = 0.5) -> Vector{Matrix{Float64}}

Per-channel positive affine `a·x + b` (a > 0) applied to the SIGNAL pixels (x > 0),
leaving the exact-zero (dead/undershoot-clamped) pixels at zero. Pearson is exactly
scale/shift invariant, so the 8×8 summary is UNCHANGED — the strongest theoretical
guarantee. Preserving the zero set is load-bearing: a blanket `+b` would turn the
`max(·,0)` read-noise zeros into signal and change which pixels the src `_exclude_zero`
drops per patch, perturbing a few patch correlations (empirically ~3/64). Gating the
shift on `x > 0` keeps the excluded set identical, so per-patch Pearson is invariant and
KS ≈ 0 robustly. `b ≥ 0` keeps signal pixels positive.
"""
function negctrl_affine(pair; a::Real = 2.0, b::Real = 0.5)
    a > 0 || throw(ArgumentError("negctrl_affine: slope a must be > 0, got $a"))
    return [Matrix{Float64}(map(v -> v > 0 ? a * v + b : v, ch)) for ch in pair]
end

"""
    negctrl_rotate_flip(pair; k = 1) -> Vector{Matrix{Float64}}

Global 90°·k rotation (no interpolation). On a square image this permutes WHICH patch is
where but preserves the SET of 8×8 per-patch correlation values, so the summary
DISTRIBUTION is unchanged (KS on the 64 values passes).
"""
function negctrl_rotate_flip(pair; k::Integer = 1)
    rot(x) = (kk = mod(k, 4); kk == 0 ? x : kk == 1 ? rotl90(x) : kk == 2 ? rot180(x) : rotr90(x))
    return [Matrix{Float64}(rot(ch)) for ch in pair]
end

"""
    negctrl_block_permute(rng, pair; blocks = 8) -> Vector{Matrix{Float64}}

Distribution-preserving spatial rearrangement: partition into `blocks × blocks` tiles and
apply the SAME random tile permutation to BOTH channels. Each patch's paired pixels move
together, so the per-patch correlation values are preserved as a SET (KS-invariant). With
`blocks = 8` the tiles align with the 8×8 patch grid.
"""
function negctrl_block_permute(rng, pair; blocks::Integer = 8)
    H, W = size(pair[1])
    (H % blocks == 0 && W % blocks == 0) ||
        throw(ArgumentError("negctrl_block_permute: $((H, W)) not divisible by blocks=$blocks"))
    bh = H ÷ blocks; bw = W ÷ blocks
    perm = randperm(rng, blocks * blocks)             # SAME permutation for both channels
    function permute(x)
        out = similar(x)
        for (dst, src) in enumerate(perm)
            di, dj = fldmod1(dst, blocks); si, sj = fldmod1(src, blocks)
            out[(di-1)*bh+1:di*bh, (dj-1)*bw+1:dj*bw] = x[(si-1)*bh+1:si*bh, (sj-1)*bw+1:sj*bw]
        end
        return out
    end
    return [Matrix{Float64}(permute(ch)) for ch in pair]
end

"""
    verify_summary_invariance(pair, transform) -> Float64

Empirically verify (D-04) that `transform` leaves the 8×8 patch-Pearson summary
statistically UNCHANGED: compute the summary before and after, then return the
two-sample KS statistic (`ApproximateTwoSampleKSTest.δ`) between the two multisets of
non-missing correlation values. A correlation-preserving transform gives δ ≈ 0; call
sites assert `δ < OOD_KS_EPS`. Returns `NaN` if either summary is fully missing.
"""
function verify_summary_invariance(pair, transform)
    v0 = Float64.(collect(skipmissing(patch_summary(build_mci(pair)))))
    v1 = Float64.(collect(skipmissing(patch_summary(build_mci(transform(pair))))))
    (isempty(v0) || isempty(v1)) && return NaN
    return ApproximateTwoSampleKSTest(v0, v1).δ
end

# ============================================================================
# ROC over the misspecification grid (D-03/D-04/D-06) — reported by 05-04
# ============================================================================

"""
    ood_roc_over_grid(m, nulls; families = OOD_FAMILIES, levels = OOD_GRID_LEVELS,
                      n_id = 100, n_pos = 40, rng = val_rng(), imsize = SBC_IMSIZE,
                      with_pp = true, pp_reps = OOD_PP_REPS) -> NamedTuple

The controlled ROC/AUC experiment (OOD-02). Draws `n_fit` FRESH train-only ID images to fit
the image-noise null (Channel 3) and the per-channel robust-z fusion reference, `n_id` FRESH
in-distribution negatives, and, for each positive-control family × magnitude level, `n_pos`
fully-external positives (SC5). For every sample it scores BOTH the summary-density
(Mahalanobis) channel and the image-noise (Mahalanobis) channel and OR-FUSES them into one
detector (per-sample max of robust-z scores, D-05); when `with_pp` the posterior-predictive
channel additionally contributes to the flag.

`maha_auc[fam]` is the REPORTED (gated) per-level AUC of the OR-FUSED density∨noise detector
— renamed in spirit from the pure-density channel because 05-03 named the detector-noise
family a MEASURED blind spot of the correlation-only summary, which this fused detector now
NARROWS. `density_auc[fam]` / `noise_auc[fam]` are the two standalone channels, returned for
transparency. Also returns the OR-flag fire-rate at the PRE-REGISTERED ID-quantile operating
point (`id_threshold` on the FUSED ID scores, committed before misspec data, D-06); the
DENSITY threshold `maha_thr` (still used by the run_ood negative-control flag); the pooled
strongest-level combined ROC curve (for `plot_ood_roc`); the in-distribution fire-rate; the
`noise_null` + `zref` (so a caller can reconstruct the fused flag on arbitrary images); and
the POST-HOC `youden_j` reference (never the gate). The full OOD_GRID_LEVELS reported run is
owned by Wave-3 (05-04); the fast gate calls this with small `n_id`/`n_pos`/`with_pp=false`.
"""
function ood_roc_over_grid(m, nulls; families = OOD_FAMILIES, levels::Integer = OOD_GRID_LEVELS,
                           n_id::Integer = 100, n_pos::Integer = 40, n_fit::Integer = n_id,
                           rng = val_rng(), imsize = SBC_IMSIZE,
                           with_pp::Bool = true, pp_reps::Integer = OOD_PP_REPS)
    # --- (A) TRAIN-ONLY noise null + per-channel robust-z reference (SC5/D-06) -----
    #     A fresh ID image pool, DISJOINT from the ROC negatives drawn next off the same
    #     sequential stream. It fits the image-noise Mahalanobis null AND the density/noise
    #     robust-z medians/scales used to put the two channels on a comparable footing for
    #     OR-fusion. Fit on ID (training-distribution) draws ONLY — never on misspec data.
    fit_pairs = [simulate_pair(rng, sample_prior(rng); imsize = imsize) for _ in 1:n_fit]
    Ffit  = reduce(hcat, [noise_features(p) for p in fit_pairs])
    nnull = fit_noise_null(Ffit)
    fit_maha  = [maha_score(nulls, frozen_summary(m, p)) for p in fit_pairs]
    fit_noise = [noise_score(nnull, p) for p in fit_pairs]
    med_d = median(fit_maha);  scl_d = _mad(fit_maha)  + 1e-12
    med_n = median(fit_noise); scl_n = _mad(fit_noise) + 1e-12
    # OR-fusion score (D-05): the per-sample max of the two channels' robust-z scores —
    # the continuous analogue of "fire iff EITHER channel exceeds its ID-quantile point".
    zfuse(sd, sn) = max((sd - med_d) / scl_d, (sn - med_n) / scl_n)

    # --- (B) ID negatives (fresh; the ROC negatives + the pre-registered operating point) -
    id_pairs = [simulate_pair(rng, sample_prior(rng); imsize = imsize) for _ in 1:n_id]
    Zid      = [frozen_summary(m, p) for p in id_pairs]
    id_maha  = [maha_score(nulls, z) for z in Zid]
    maha_thr = id_threshold(id_maha)                   # density threshold (run_ood neg-ctrl)
    id_noise = [noise_score(nnull, p) for p in id_pairs]
    id_fused = [zfuse(id_maha[k], id_noise[k]) for k in 1:n_id]
    fused_thr = id_threshold(id_fused)                 # OR-fused operating point (D-06)
    id_pp = with_pp ?
        [pp_mismatch_score(m, z; reps = pp_reps, rng = val_rng(), imsize = imsize) for z in Zid] :
        Float64[]
    pp_thr = with_pp ? id_threshold(id_pp) : Inf

    fam_names   = keys(families)
    maha_auc    = Dict{Symbol,Vector{Float64}}()       # REPORTED detector = density∨noise fusion
    density_auc = Dict{Symbol,Vector{Float64}}()       # standalone density channel (transparency)
    noise_auc   = Dict{Symbol,Vector{Float64}}()       # standalone image-noise channel (transparency)
    pp_auc      = Dict{Symbol,Vector{Float64}}()
    fire_rate   = Dict{Symbol,Vector{Float64}}()

    pooled_pos_fused = Float64[]                        # strongest level, all families (ROC curve)
    for fam in fam_names
        gen = families[fam]
        maha_auc[fam]    = Float64[]
        density_auc[fam] = Float64[]
        noise_auc[fam]   = Float64[]
        pp_auc[fam]      = Float64[]
        fire_rate[fam]   = Float64[]
        for lvl in 1:levels
            pos_maha = Float64[]; pos_noise = Float64[]; pos_fused = Float64[]
            pos_pp = Float64[]; fires = 0
            for _ in 1:n_pos
                θ    = sample_prior(rng)
                img  = gen(rng, θ; imsize = imsize, level = lvl)
                Zp   = frozen_summary(m, img)
                ms   = maha_score(nulls, Zp);   push!(pos_maha,  ms)
                ns   = noise_score(nnull, img); push!(pos_noise, ns)
                fs   = zfuse(ms, ns);           push!(pos_fused, fs)
                fired = fs > fused_thr
                if with_pp
                    ps = pp_mismatch_score(m, Zp; reps = pp_reps, rng = val_rng(), imsize = imsize)
                    push!(pos_pp, ps)
                    fired = fired || (ps > pp_thr)
                end
                fires += fired ? 1 : 0
                lvl == levels && push!(pooled_pos_fused, fs)
            end
            push!(maha_auc[fam],    roc_auc(id_fused, pos_fused)[3])   # FUSED — the gated AUC
            push!(density_auc[fam], roc_auc(id_maha,  pos_maha)[3])
            push!(noise_auc[fam],   roc_auc(id_noise, pos_noise)[3])
            with_pp && push!(pp_auc[fam], roc_auc(id_pp, pos_pp)[3])
            push!(fire_rate[fam], fires / n_pos)
        end
    end

    id_fire_rate = mean(id_fused .> fused_thr)
    fpr, tpr, or_auc = roc_auc(id_fused, pooled_pos_fused)
    yj = youden_j(fpr, tpr)                            # POST-HOC reference only (D-06)

    return (levels = collect(1:levels), families = fam_names,
            maha_thr = maha_thr, pp_thr = pp_thr, fused_thr = fused_thr,
            maha_auc = maha_auc, density_auc = density_auc, noise_auc = noise_auc,
            pp_auc = pp_auc, fire_rate = fire_rate,
            id_fire_rate = id_fire_rate, noise_null = nnull,
            zref = (med_d = med_d, scl_d = scl_d, med_n = med_n, scl_n = scl_n),
            roc_fpr = fpr, roc_tpr = tpr, combined_auc = or_auc,
            youden = yj, auc_min = OOD_AUC_MIN, id_quantile = OOD_ID_QUANTILE)
end

"""
    plot_ood_roc(fpr, tpr, auc; filename = "ood_roc.png") -> String

Write the OOD ROC PNG by CALLING figures.jl's owned `plot_roc` (D-05 figure). figures.jl
is loaded LAZILY here so the fast test gate never pulls the CairoMakie stack (figures are
script artifacts, not gate assertions — HARD RULE in figures.jl). Called by the reported
Wave-3 run (05-04), never by test_ood.jl.
"""
function plot_ood_roc(fpr, tpr, auc; filename::AbstractString = "ood_roc.png")
    isdefined(@__MODULE__, :plot_roc) || include(joinpath(@__DIR__, "figures.jl"))
    return plot_roc(fpr, tpr, auc; filename = filename)
end
