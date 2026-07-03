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
# src/amortized/train_ratio.jl --- amortized-BF NRE training (PROD-01, D-06 GPU-plumbed).
#
# Promoted (grid-general, GPU-plumbed) from the proven spike:
#   spike/validation/train_ratio.jl:208-229  assemble_ratio_data_cache (cache-paired data)
#   spike/validation/train_ratio.jl:240-244  measure_log_prior_odds
#   spike/validation/train_ratio.jl:278-321  train_ratio (fixed-data train, use_gpu delta)
#
# Trains a NeuralEstimators `RatioEstimator` as a model-comparison / evidence net over the binary
# model index m ∈ {0,1}: m = 1 ⟺ coloc (Δρ > 0), m = 0 ⟺ null. Reading `logratio(Z; m=1) −
# logratio(Z; m=0)` in ONE forward pass is the amortized log Bayes factor (bf.jl).
#
# THREE DELTAS on promotion:
#   (1) D-06 GPU: `use_gpu` DEFAULTS to `has_cuda_device()`; the spike's `use_gpu && throw(...)`
#       CPU-only guard is REMOVED. Shipped BF inference stays CPU (bf.jl default `use_gpu=false`).
#   (2) GRID-GENERAL: the conditioner input width derives from `ratio_input_dim(G) = 5·G²`
#       (summary.jl), NOT the literal 320. `assemble_ratio_pairs` pairs columns of an
#       already-`zt`-standardized summary matrix (the pipeline standardizes) via the grid-general
#       `pair_encode` (bf.jl) — no cache/loader dependency.
#   (3) HARDCODED LOSS: NO custom loss is passed to `train` — `RatioEstimator`'s loss is
#       hard-coded `logitbinarycrossentropy` and a passed loss is SILENTLY IGNORED (v0.2.1).
#
# AdamW-Float64 gotcha (mirrors train_npe): LR / weight-decay stay Float64.
#############################################################################################

import NeuralEstimators: RatioEstimator, train
import Flux
import Flux: Chain, Dense, gelu
import Statistics: mean
import Random123: Philox4x
import Random: rand

# --- Ratio-net constants (Phase-5 BF-iteration values) --------------------------------------
const RATIO_MODEL_SCHEMA    = 1
const RATIO_NUM_SUMMARIES   = 64      # learned-summary bottleneck feeding the model-index head
const RATIO_SUMMARY_WIDTH   = 256     # summary-conditioner hidden width (mirrors the NPE lever)
const RATIO_SPLIT_THRESHOLD = 0.0     # D-07 coloc/null split on the ρ_true CONTRAST (≡ {Δμ>0})

# Deterministic, disjoint stream for CACHE/POOL PAIRING (which two summary columns form a
# training pair + the train/val split). Distinct from the datagen master seed and the eval
# streams, so the assembled training pairs never coincide with the reported BF eval stream.
const RATIO_PAIR_SEED = 0x00000000004A7107

"""
    build_ratio_estimator(input_dim::Integer; num_summaries = RATIO_NUM_SUMMARIES,
                          width = RATIO_SUMMARY_WIDTH) -> RatioEstimator

Build the model-comparison `RatioEstimator`: a summary conditioner
`Chain(Dense(input_dim, W, gelu), Dense(W, W, gelu), Dense(W, W, gelu), Dense(W, num_summaries))`
wrapped as `RatioEstimator(net, 1; num_summaries)` (num_parameters = 1 = the binary model index).
`input_dim` is `ratio_input_dim(G) = 5·G²` (grid-general) — NOT the literal 320. Shared by the
trainer and by `persist.load_ratio` so training and reload build the identical topology.
"""
function build_ratio_estimator(input_dim::Integer;
                               num_summaries::Integer = RATIO_NUM_SUMMARIES,
                               width::Integer = RATIO_SUMMARY_WIDTH)
    input_dim >= 1 || throw(ArgumentError("build_ratio_estimator: input_dim must be ≥ 1"))
    W = width
    net = Chain(Dense(input_dim, W, gelu), Dense(W, W, gelu),
                Dense(W, W, gelu), Dense(W, num_summaries))
    return RatioEstimator(net, 1; num_summaries = num_summaries)
end

"""
    measure_log_prior_odds(model_index) -> Float64

The MEASURED log prior-odds of the coloc event (Pitfall 5 — measured, NEVER assumed 0):
`log(p / (1 − p))` with `p = mean(model_index)` the prior fraction of coloc labels. Because the
paired θ are i.i.d. from the SAME prior, `p ≈ 0.5` by symmetry — but it is MEASURED (not assumed)
to absorb finite-sample/tie imbalance (D-07). Returns a finite scalar; errors on degenerate
balance.
"""
function measure_log_prior_odds(model_index::AbstractVector{<:Real})
    p = mean(model_index)
    (0.0 < p < 1.0) ||
        error("measure_log_prior_odds: degenerate label balance p=$p (need 0<p<1)")
    return log(p / (1 - p))
end

"""
    assemble_ratio_pairs(Zstd, ρ, n; rng = Philox4x(UInt64, (RATIO_PAIR_SEED, UInt64(0))),
                         threshold = RATIO_SPLIT_THRESHOLD) -> (Z_pair, model_index)

Build `n` paired training columns for the model-comparison ratio net (BF-01, D-07) by PAIRING
distinct columns of an already-FROZEN-`zt`-standardized summary matrix `Zstd` (`2·G²×N`, the
deployed NPE input space) with its ρ_true row `ρ` (`θ[1, :]`). For each column: draw two DISTINCT
pool indices `(i, j)` from `rng`, set the column to `pair_encode(Zstd[:,i], Zstd[:,j])` (length
`ratio_input_dim(G) = 5·G²`, grid-general) and the label to
`Float32((ρ[i] − ρ[j]) > threshold)` — the D-07 coloc/null split on the ρ_true contrast. Because
`(i, j)` are i.i.d. over the pool, `P(ρ_i > ρ_j) ≈ 0.5`, so labels are balanced by construction.

Returns `(Z_pair::Matrix{Float32} 5G²×n, model_index::Vector{Float32} length n)`.
"""
function assemble_ratio_pairs(Zstd::AbstractMatrix, ρ::AbstractVector, n::Integer;
                              rng = Philox4x(UInt64, (RATIO_PAIR_SEED, UInt64(0))),
                              threshold::Real = RATIO_SPLIT_THRESHOLD)
    N = size(Zstd, 2)
    N >= 2 || error("assemble_ratio_pairs: need ≥2 summary columns, got $N")
    length(ρ) == N ||
        throw(DimensionMismatch("assemble_ratio_pairs: ρ length $(length(ρ)) != Zstd cols $N"))

    G   = isqrt(size(Zstd, 1) ÷ 2)                 # 2·G² rows ⇒ grid G (single-source coupling)
    idm = ratio_input_dim(G)                        # 5·G² — grid-general conditioner width
    Z_pair      = Matrix{Float32}(undef, idm, n)
    model_index = Vector{Float32}(undef, n)
    for k in 1:n
        i1 = rand(rng, 1:N)
        i2 = rand(rng, 1:N)
        while i2 == i1
            i2 = rand(rng, 1:N)
        end
        Z_pair[:, k]      = pair_encode(view(Zstd, :, i1), view(Zstd, :, i2))
        model_index[k]    = Float32((ρ[i1] - ρ[i2]) > threshold)
    end
    return Z_pair, model_index
end

"""
    train_ratio(Zstd, ρ; n = 48_000, use_gpu = has_cuda_device(), epochs = 300,
                batchsize = 128, learning_rate = 2.5e-4, weight_decay = 1e-4, val_frac = 0.15,
                stopping_epochs = 40, num_summaries = RATIO_NUM_SUMMARIES,
                summary_width = RATIO_SUMMARY_WIDTH,
                pair_rng = Philox4x(UInt64, (RATIO_PAIR_SEED, UInt64(0))),
                verbose = false) -> NamedTuple

Train the model-comparison `RatioEstimator` on `n` pairs assembled from the FROZEN-`zt`
standardized summary pool `Zstd` (`2·G²×N`) and its ρ_true row `ρ`. Builds
`build_ratio_estimator(ratio_input_dim(G))`, measures `log_prior_odds` from the assembled labels,
splits train/val, and calls the FIXED-DATA `train(est, midx_tr, midx_va, Z_tr, Z_va; ...)` with
AdamW weight decay and early stopping.

**D-06:** `use_gpu` defaults to `has_cuda_device()` (GPU-train allowed, graceful CPU fallback);
NO `use_gpu && throw` guard. **No custom loss** is passed (the ratio loss is the hard-coded
`logitbinarycrossentropy`; a passed loss is silently ignored). Returns the ratio HANDLE
`(estimator, log_prior_odds, num_summaries, summary_width, input_dim)`.
"""
function train_ratio(Zstd::AbstractMatrix, ρ::AbstractVector;
                     n::Integer = 48_000, use_gpu::Bool = has_cuda_device(),
                     epochs::Integer = 300, batchsize::Integer = 128,
                     learning_rate::Real = 2.5e-4, weight_decay::Real = 1e-4,
                     val_frac::Real = 0.15, stopping_epochs::Integer = 40,
                     num_summaries::Integer = RATIO_NUM_SUMMARIES,
                     summary_width::Integer = RATIO_SUMMARY_WIDTH,
                     pair_rng = Philox4x(UInt64, (RATIO_PAIR_SEED, UInt64(0))),
                     verbose::Bool = false)
    Z_pair, model_index = assemble_ratio_pairs(Zstd, ρ, n; rng = pair_rng)
    log_prior_odds = measure_log_prior_odds(model_index)

    # Train/val split (contiguous — the stream is already i.i.d. across columns).
    nval = clamp(round(Int, val_frac * n), 1, n - 1)
    ntr  = n - nval
    Z_tr = Z_pair[:, 1:ntr]
    Z_va = Z_pair[:, ntr+1:end]
    midx_tr = reshape(model_index[1:ntr], 1, :)         # 1×ntr (num_parameters = 1)
    midx_va = reshape(model_index[ntr+1:end], 1, :)

    input_dim = size(Z_pair, 1)                          # == ratio_input_dim(G)
    est = build_ratio_estimator(input_dim; num_summaries = num_summaries, width = summary_width)

    # FIXED-DATA train form. D-06: `use_gpu` plumbed. NO custom loss (hard-coded logit-BCE).
    # AdamW args Float64 to match NeuralEstimators' Float64 CosAnneal lr_schedule.
    est = train(est, midx_tr, midx_va, Z_tr, Z_va;
                epochs = epochs, batchsize = batchsize, use_gpu = use_gpu,
                optimiser = Flux.Optimisers.AdamW(learning_rate, (0.9, 0.999), weight_decay),
                stopping_epochs = stopping_epochs, verbose = verbose)

    return (estimator = est, log_prior_odds = log_prior_odds,
            num_summaries = num_summaries, summary_width = summary_width,
            input_dim = input_dim)
end
