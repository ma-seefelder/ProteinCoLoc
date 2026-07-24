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
# src/amortized/api.jl --- the PUBLIC amortized colocalization entry point (PROD-01, D-01).
#
# `colocalization_amortized(img, control, channels; num_patches, N)` is the v2.0 primary public
# API: ONE registry lookup + a few frozen CPU forward passes yield an `AmortizedColocResult`
# (posterior draws, Δρ cloud, amortized log Bayes factor, OOD verdict, calibration provenance).
# It composes the already-promoted read surfaces — it defines NO new inference math:
#
#   estimator_for(G)                         # registry: lazy-load the shipped 8×8 bundle (PROD-02)
#   patch_summary → encode_d01               # grid-parametric summary (summary.jl)
#   standardize_summary(·, b.zt, :min)       # FROZEN zt, mask rows bypassed (infer.jl, Pitfall 5)
#   posterior_for / rho_draws                # single-pass NPE posterior, un-standardized (infer.jl)
#   amortized_log_bf(pair_encode(Zs, Zc))    # one-pass NRE log Bayes factor (bf.jl)
#   ood_verdict(b.ood_nulls, Zs, Zc)         # density-channel OOD flag at the recorded op-point (ood.jl)
#
# CPU-DEFAULT (D-06): `use_gpu = false` on EVERY inference call — the shipped/reproducible path is
# deterministic regardless of the hardware the net was TRAINED on. `use_gpu` stays a plumbed keyword.
#
# GATE-VERDICT HONESTY: the returned `calibration` is the grid's recorded ship-gate report; the
# 8×8 GO carries NAMED LIMITS (twice-amended gate; ρ_true randomized-rank atom handling;
# simulation-based BF validation, not a per-pair KDE reference; ~0.08-SD nuisance-marginal drift;
# ±1.5-z one-shot seed variance). See `docs/amortized.md` and `07-GO-NO-GO-UPDATE.md §5`.
#############################################################################################

import StatsBase

# Build a 2-channel MultiChannelImage from the two SELECTED channels of a parent image, carrying the
# parent's Otsu thresholds (mirrors the training/gate `build_mci` and the 07-08 per-tile `_tile_mci`
# construction, so the public summary rides the SAME frozen read chain training did). `patch_summary`
# reads `.data[1]`/`.data[2]`, i.e. exactly the selected pair.
function _pair_mci(mci::MultiChannelImage, channels::AbstractVector{<:Integer})
    c1, c2 = channels[1], channels[2]
    x = mci.data[c1]; y = mci.data[c2]
    return MultiChannelImage(
        [x, y],
        [mci.channels[c1], mci.channels[c2]],
        mci.name,
        ["", ""],
        size(x),
        [mci.otsu_threshold[c1], mci.otsu_threshold[c2]],
    )
end

# Frozen, standardized summary vector for one image + its two selected channels (raw encode → frozen
# zt; mask rows bypassed). The one place the public path touches the summary chain.
function _frozen_summary(mci::MultiChannelImage, channels::AbstractVector{<:Integer}, b, G::Integer)
    pair = _pair_mci(mci, channels)
    return standardize_summary(encode_d01(patch_summary(pair, G)), b.zt, :min)
end

"""
    colocalization_amortized(img, control, channels; num_patches = 8, N = 2000, use_gpu = false)
        -> AmortizedColocResult

The v2.0 PUBLIC amortized colocalization entry point (D-01): infer a calibrated colocalization
posterior, a Δρ (sample − control) contrast cloud, an amortized log Bayes factor, and an OOD /
misspecification flag for a `(img, control)` pair, in a handful of frozen CPU forward passes through
the `num_patches × num_patches` estimator bundle.

# Arguments
- `img::MultiChannelImage`: the sample image.
- `control::MultiChannelImage`: the control image (Δρ is `ρ(img) − ρ(control)`, D-03).
- `channels::AbstractVector{<:Integer}`: the TWO channels to compare (e.g. `[1, 2]`).

# Keywords
- `num_patches::Integer = 8`: the registry key (PROD-02). Only shipped grids resolve; an
  unregistered grid raises an `ArgumentError` that names the shipped grids and `train_and_register`
  (T-7-05 — never a silent default). v2.0 ships **only 8×8** (`07-GO-NO-GO-UPDATE.md`, Option A).
- `N::Integer = 2000`: posterior draws per stack.
- `use_gpu::Bool = false`: CPU-reproducible shipped default (D-06); a plumbed opt-in only.

# Returns
An [`AmortizedColocResult`](@ref) with the 7×N physical-θ posterior (row 1 = ρ_true), the length-N
Δρ draw cloud, the scalar log Bayes factor, the [`OODVerdict`](@ref), and the grid's recorded
[`CalibrationMeta`](@ref) ship-gate provenance. Read it through the shared accessors
(`posterior_draws`, `delta_rho`, `bayes_factor`, `is_ood`).

# Named limits (non-negotiable honesty — carried into the manuscript)
The 8×8 calibration rests on a POST-HOC re-reading: the ship-gate was amended twice; ρ_true SBC
uses randomized-rank atom handling (the prior places mass on non-realizable μ); the Bayes factor is
validated by simulation-based DISCRIMINATION (spike 014, AUC 0.994), NOT magnitude agreement with a
per-pair KDE reference (shown invalid past |logBF| ≈ log(L)); nuisance-parameter marginals carry a
~0.08-SD residual drift not certifiable as negligible at M=2000; and the one-shot gate carries
±1.5-z per-parameter seed variance. See `docs/amortized.md` and `07-GO-NO-GO-UPDATE.md §5`. The
shipped OOD flag is the summary-density (Mahalanobis) channel at the recorded ID operating point;
the fused noise/PP channels that reached gate AUC 1.0 are gate-time constructs (they need the raw
image pair / re-simulation and are not carried in the frozen bundle).
"""
function colocalization_amortized(img::MultiChannelImage, control::MultiChannelImage,
                                  channels::AbstractVector{<:Integer};
                                  num_patches::Integer = 8, N::Integer = 2000,
                                  use_gpu::Bool = false)
    # --- input validation (V5) --------------------------------------------------------------
    length(channels) == 2 || throw(ArgumentError(
        "colocalization_amortized: `channels` must select exactly 2 channels, got $(length(channels))."))
    for (nm, m) in (("img", img), ("control", control))
        for c in channels
            (1 <= c <= length(m.data)) || throw(ArgumentError(
                "colocalization_amortized: channel $c out of range for $nm (has $(length(m.data)) channels)."))
        end
    end
    N >= 1 || throw(ArgumentError("colocalization_amortized: N must be ≥ 1, got $N."))

    G = Int(num_patches)
    b = estimator_for(G)                       # registry lookup / lazy artifact load (T-7-05)

    # --- frozen summaries (raw encode → frozen zt; mask rows bypassed) -----------------------
    Zs = _frozen_summary(img,     channels, b, G)
    Zc = _frozen_summary(control, channels, b, G)
    (all(isfinite, Zs) && all(isfinite, Zc)) || throw(ArgumentError(
        "colocalization_amortized: non-finite standardized summary — the image/channel data is " *
        "degenerate for a $(G)×$(G) grid (check for NaN/Inf pixels or an all-background stack)."))

    # --- amortized reads (all CPU-default, D-06) --------------------------------------------
    post_std = posterior_for(b.npe, Zs; N = N, use_gpu = use_gpu)          # 7×N standardized θ
    post     = Matrix{Float64}(StatsBase.reconstruct(b.θzt, post_std))     # 7×N physical θ (Pitfall 5)

    ρs = rho_draws(b.npe, Zs, b.θzt; N = N, use_gpu = use_gpu)
    ρc = rho_draws(b.npe, Zc, b.θzt; N = N, use_gpu = use_gpu)
    Δρ = collect(Float64, ρs .- ρc)                                        # D-03 Δρ draw cloud

    logbf = Float64(amortized_log_bf(b.ratio.estimator, pair_encode(Zs, Zc),
                                     b.ratio.log_prior_odds; use_gpu = use_gpu))

    ood = ood_verdict(b.ood_nulls, Zs, Zc)                                 # density-channel OOD flag

    return AmortizedColocResult(G, post, Δρ, logbf, ood, b.calibration,
                                (; N = Int(N), channels = collect(Int, channels),
                                   num_patches = G, use_gpu = use_gpu))
end
