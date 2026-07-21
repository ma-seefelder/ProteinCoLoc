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
# src/amortized/train_npe.jl --- amortized-NPE training (PROD-01, D-06 GPU-plumbed).
#
# Promoted (grid-general, GPU-plumbed) from the proven spike:
#   spike/npe/train_npe.jl:80        fit_theta_transform (leak-free θ ZScoreTransform)
#   spike/npe/train_npe.jl:102-149   train_fold          (fixed-data train, use_gpu delta)
#
# THREE DELTAS on promotion (07-PATTERNS §use_gpu split):
#   (1) D-06 GPU: `use_gpu` DEFAULTS to `has_cuda_device()` (GPU when a device is present,
#       graceful CPU fallback). The spike's `use_gpu && throw(...)` CPU-only guard is REMOVED
#       from the training path. Persisted models + the ship-gate + shipped inference stay
#       CPU-reproducible (persist.jl / infer.jl default `use_gpu=false`), so training on GPU
#       never leaks device non-determinism into the frozen, pre-registered object.
#   (2) The leak-free summary `zt` and the fixed-data train form are decoupled from the spike's
#       cache `loader.jl`: `train_npe` takes ALREADY-STANDARDIZED train/val summaries (the
#       pipeline fits `zt` on TRAIN continuous rows only, via `fit_summary_transform`) plus the
#       RAW train/val θ, and fits the θ transform itself (leak-free).
#   (3) src module convention (no standalone `isdefined || include` guards; the module's ordered
#       includes provide `build_estimator`, `_summary_row_partition`, `has_cuda_device`).
#
# AdamW-Float64 GOTCHA (train_npe.jl:129-131): the learning-rate / weight-decay stay Float64 to
# match NeuralEstimators' Float64 CosAnneal `lr_schedule` — a Float32-parameterised optimiser
# state trips `Optimisers.adjust!` when the schedule feeds a Float64 eta.
#############################################################################################

import NeuralEstimators: train
import Flux
import StatsBase

"""
    fit_theta_transform(θtr::AbstractMatrix; bounds = theta_prior_bounds()) -> BoundedThetaTransform

Fit the leak-free FROZEN θ transform on the TRAIN θ columns ONLY (`dims = 2`).

**BOUNDED θ-SPACE (07-CALIBRATION-FINDINGS F2).** θ is FIRST mapped out of the truncated prior box
onto ℝ by a per-parameter logit bijection (`theta_to_unbounded`, bounds from the single-source
`theta_prior_bounds()`), and the `ZScoreTransform` is fitted in THAT unconstrained space. The
`NormalisingFlow` therefore models an unbounded quantity — which is what it can actually
represent — and `StatsBase.reconstruct` maps every posterior draw back STRICTLY INSIDE the prior
support. Previously this was a plain `ZScoreTransform` on raw θ, which left the flow unable to
represent the box: measured out-of-support leakage plus a location-biased `label_efficiency`
conditional (SBC z = +7.4 at G = 8).

`bounds` is a keyword ONLY so tests can exercise a different box; production callers must let it
default so the prior definition stays the single source of truth. Rows whose logit-space spread is
degenerate (a synthetic all-identical θ row) get a unit scale instead of a 0 divisor.
"""
function fit_theta_transform(θtr::AbstractMatrix; bounds = theta_prior_bounds())
    lo = Float64[b[1] for b in bounds]
    hi = Float64[b[2] for b in bounds]
    size(θtr, 1) == length(lo) || throw(DimensionMismatch(
        "fit_theta_transform: θ has $(size(θtr,1)) rows but the prior box has $(length(lo))"))
    U  = theta_to_unbounded(θtr, lo, hi)
    zt = StatsBase.fit(StatsBase.ZScoreTransform, U; dims = 2)
    @inbounds for i in eachindex(zt.scale)      # degenerate-row guard (never a 0 divisor)
        (isfinite(zt.scale[i]) && zt.scale[i] > 0) || (zt.scale[i] = 1.0)
    end
    return BoundedThetaTransform(lo, hi, zt)
end

"""
    fit_summary_transform(Ztr_raw::AbstractMatrix; variant::Symbol = :min) -> ZScoreTransform

Fit the leak-free summary `ZScoreTransform` on the TRAIN-only RAW summaries `Ztr_raw` (a
`2·G²×Ntrain` matrix), restricted to the CONTINUOUS rows via `_summary_row_partition` (the binary
mask rows are NEVER z-scored — z-scoring a 0/1 mask would re-couple folds through the mask mean).
The returned `zt` is the frozen standardizer `standardize_summary` (infer.jl) applies at
inference and the ratio/OOD training see, so training and deployment share ONE input space.
"""
function fit_summary_transform(Ztr_raw::AbstractMatrix; variant::Symbol = :min)
    cont, _ = _summary_row_partition(variant, size(Ztr_raw, 1))
    return StatsBase.fit(StatsBase.ZScoreTransform, Ztr_raw[cont, :]; dims = 2)
end

"""
    train_npe(Ztr_std, Zva_std, θtr, θva; use_gpu = has_cuda_device(), epochs = 300,
              batchsize = 128, learning_rate = 2.5e-4, weight_decay = 1e-4,
              dstar = NPE_DSTAR, depth = NPE_DEPTH, width = NPE_WIDTH,
              num_coupling_layers = NPE_COUPLING, flow_depth = NPE_FLOW_DEPTH,
              flow_width = NPE_FLOW_WIDTH, stopping_epochs = 40, verbose = false) -> NamedTuple

Train the amortized NPE on the ALREADY-STANDARDIZED train/val summaries `Ztr_std`/`Zva_std`
(`d_in×n`, frozen-`zt` space; the pipeline standardizes) and the RAW train/val θ `θtr`/`θva`
(`7×n`). Fits the leak-free BOUNDED θ transform on `θtr` (logit-of-prior-box ∘ z-score, F2), builds
`build_estimator(size(Ztr_std,1), ...)`, and calls the FIXED-DATA
`train(est, θtr_std, θva_std, Ztr_std, Zva_std; use_gpu, ...)` with AdamW weight decay and early
stopping.

**D-06:** `use_gpu` defaults to `has_cuda_device()` (GPU when available, graceful CPU fallback);
NO `use_gpu && throw` guard. LR / weight-decay stay Float64 (AdamW-CosAnneal gotcha).

Returns `(estimator, θzt, d_in, variant, arch)` where `θzt` is the frozen θ-transform and `arch`
is the architecture metadata `persist.save_estimator` needs to rebuild the net for
`Flux.loadmodel!` (Pitfall 4).
"""
function train_npe(Ztr_std::AbstractMatrix, Zva_std::AbstractMatrix,
                   θtr::AbstractMatrix, θva::AbstractMatrix;
                   use_gpu::Bool = has_cuda_device(),
                   epochs::Integer = 300, batchsize::Integer = 128,
                   learning_rate::Real = 2.5e-4, weight_decay::Real = 1e-4,
                   dstar::Integer = NPE_DSTAR, depth::Integer = NPE_DEPTH,
                   width::Integer = NPE_WIDTH,
                   num_coupling_layers::Integer = NPE_COUPLING,
                   flow_depth::Integer = NPE_FLOW_DEPTH,
                   flow_width::Integer = NPE_FLOW_WIDTH,
                   stopping_epochs::Integer = 40, verbose::Bool = false)
    size(Ztr_std, 2) == size(θtr, 2) ||
        throw(DimensionMismatch("train_npe: Ztr_std cols $(size(Ztr_std,2)) != θtr cols $(size(θtr,2))"))
    size(Zva_std, 2) == size(θva, 2) ||
        throw(DimensionMismatch("train_npe: Zva_std cols $(size(Zva_std,2)) != θva cols $(size(θva,2))"))

    # Leak-free θ-standardization (Pitfall 5): fit on train θ only, apply to both.
    θzt     = fit_theta_transform(θtr)
    θtr_std = Float32.(StatsBase.transform(θzt, θtr))
    θva_std = Float32.(StatsBase.transform(θzt, θva))

    d_in = size(Ztr_std, 1)
    est  = build_estimator(d_in; dstar = dstar, depth = depth, width = width,
                           num_coupling_layers = num_coupling_layers,
                           flow_depth = flow_depth, flow_width = flow_width)

    # FIXED-DATA train form. D-06: `use_gpu` plumbed (GPU-train allowed). AdamW args are Float64
    # to match NeuralEstimators' Float64 CosAnneal lr_schedule (Optimisers.adjust! gotcha).
    est = train(est, θtr_std, θva_std, Float32.(Ztr_std), Float32.(Zva_std);
                epochs = epochs, batchsize = batchsize, use_gpu = use_gpu,
                optimiser = Flux.Optimisers.AdamW(learning_rate, (0.9, 0.999), weight_decay),
                stopping_epochs = stopping_epochs, verbose = verbose)

    arch = (d_in = d_in, D = NPE_D, dstar = dstar, depth = depth, width = width,
            num_coupling_layers = num_coupling_layers,
            flow_depth = flow_depth, flow_width = flow_width)
    return (estimator = est, θzt = θzt, d_in = d_in, variant = :min, arch = arch)
end
