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

# spike/validation/train_ratio.jl --- Phase-5 amortized-BF ratio net (BF-01).
#
# Trains a NeuralEstimators `RatioEstimator` as a neural MODEL-COMPARISON / Evidence
# Network over a BINARY model index m ∈ {0,1}: m = 1 ⟺ coloc (Δρ > 0), m = 0 ⟺ null
# (Δρ ≤ 0). Reading `logratio(Z; m=1) − logratio(Z; m=0)` in ONE forward pass is the
# amortized log Bayes factor (bf.jl) -- provably the same posterior-odds/prior-odds
# quantity `compute_BayesFactor` computes (proof: 05-RESEARCH.md §Amortized Bayes Factor).
#
# HONESTY NOTE (05-RESEARCH.md, Pitfall 3): the pinned NeuralEstimators v0.2.1 has NO
# l-POP / Evidence-Networks loss. `RatioEstimator`'s loss is HARD-CODED
# `logitbinarycrossentropy` (RatioEstimator.jl:_loss) and a custom loss passed to
# `train` is SILENTLY IGNORED. The CLAUDE.md/CONTEXT.md "l-POP-style loss" is therefore
# re-labeled honestly as the standard Hermans-2020 NRE used as a model-comparison
# reduction. Literal l-POP (a hand-rolled Evidence Network) is a documented fallback only.
#
# D-07 (the coloc/null split mirrors compute_BayesFactor EXACTLY): the baseline event is
# {Δμ > 0} on the induced-mean contrast Δμ = μ_sample − μ_control at threshold 0. Because
# the simulator sets ρ_true = ghat(μ*) with ghat STRICTLY MONOTONE INCREASING (SIM-02,
# ghat.jl), sign(ρ_true_s − ρ_true_c) == sign(μ*_s − μ*_c) for every non-clamped pair, so
# the model index defined by the ρ_true CONTRAST at threshold 0 reproduces {Δμ>0} EXACTLY,
# INDEPENDENT of ghat(0). (Numerically ghat(0) ≈ −0.0375 ≠ 0, but that only matters for a
# POINT null on a single μ — this is a two-stack CONTRAST, so the threshold on the ρ-contrast
# is 0, matching the baseline's Δμ threshold of 0. See RATIO_SPLIT_THRESHOLD below.)
#
# FROZEN STATS (SC5 / T-05-01): BOTH summaries in every pair are standardized with the
# FROZEN train `m.zt` (via standardize_summary, the harness surface). NO ZScoreTransform is
# re-fit in this file -- the ratio net must see summaries in the deployed NPE's input space.
#
# CPU-ONLY (D-10 / Pitfall 1): `use_gpu` DEFAULTS TO TRUE in NeuralEstimators, so
# `train_ratio` hard-throws on use_gpu=true and every NeuralEstimators call passes
# use_gpu=false. CUDA is never imported.
#
# PERSISTENCE (save_npe idiom, train_npe.jl:147-186): `save_ratio` writes {estimator, zt,
# log_prior_odds, ...} to a `.tmp`, reopens to integrity-check the "estimator" key, then
# `mv(...; force=true)` (filesystem-atomic). `load_ratio` is the paired reader.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; touches no src/. Guarded includes
# keep the file loadable standalone AND idempotent under runtests.jl.

using Flux               # Chain, Dense, gelu, Flux.Optimisers.AdamW
using NeuralEstimators   # RatioEstimator, train
using JLD2               # atomic model persistence
using StatsBase          # (mean etc. come through Statistics/harness)
using Statistics         # mean
using Dates              # UTC-labeled artifact timestamp
using Random123          # Philox4x (deterministic, disjoint cache-pairing stream)
using Random             # rand over ranges / randperm on an AbstractRNG

# --- ORDER MATTERS: consts first, then the shared harness (load_frozen_model +
#     sample_prior/simulate_pair/build_mci/patch_summary/encode_d01 + standardize_summary
#     + val_rng). Guarded for idempotency under runtests.jl.
isdefined(@__MODULE__, :VAL_MASTER_SEED)   || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :load_frozen_model) || include(joinpath(@__DIR__, "harness.jl"))

# Persistence schema for trained_ratio.jld2 (bump on any breaking layout change).
if !isdefined(@__MODULE__, :RATIO_MODEL_SCHEMA)
    const RATIO_MODEL_SCHEMA = 1

    # Summary-network output width (num_summaries): the 256-dim paired summary is
    # compressed to this many learned summaries before the inference MLP. PHASE-5
    # BF ITERATION: 32→64 to mirror the NPE lever (dstar 32→64) — more capacity in the
    # learned-summary bottleneck to sharpen the model-index classifier.
    const RATIO_NUM_SUMMARIES = 64

    # Summary-network hidden width. PHASE-5 BF ITERATION: 64→256 (the NPE's 2×128→3×256
    # conditioner scale-up, mirrored on the NRE side).
    const RATIO_SUMMARY_WIDTH = 256

    # The paired-summary input dimension (research A7 difference encoding): the two
    # 128-dim :min encodings concatenated PLUS the 64-dim CONTINUOUS-ROW contrast
    # (Zs−Zc on rows 1:64) — 128+128+64 = 320. The contrast is the exact per-patch
    # correlation-difference signal the D-07 model index depends on, so exposing it
    # directly sharpens the near-threshold classifier (PHASE-5 BF ITERATION).
    const RATIO_CONT_ROWS = 64
    const RATIO_INPUT_DIM  = 320

    # D-07 coloc/null split threshold on the ρ_true CONTRAST Δρ_ρ = ρ_true_s − ρ_true_c.
    # It is 0.0 (NOT ghat(0)) because the baseline threshold on the induced-mean contrast
    # Δμ is 0 and ghat is strictly monotone through both stacks (see file banner).
    const RATIO_SPLIT_THRESHOLD = 0.0

    # Reported training defaults (NOT anti-snooping thresholds — training hyperparameters,
    # tunable; the pass/fail BF thresholds live in consts.jl). PHASE-5 BF ITERATION: the
    # net is now trained on PAIRS assembled from the SAME 50k (θ, summary) cache the NPE
    # trained on (assemble_ratio_data_cache) — leak-free (cache gen-seed 0x134d8f3 is
    # DISJOINT from VAL_MASTER_SEED and NPE_MASTER_SEED) and 16× more data than the prior
    # 3k fresh-sim pairs — with the NPE-mirrored stabler recipe (LR 2.5e-4, batch 128,
    # 300 epochs / patience 40).
    const RATIO_TRAIN_N = 48000    # cache-paired (Z_s,Z_c) training columns
    const RATIO_EPOCHS  = 300      # max epochs (early-stopped at stopping_epochs=40)
    const RATIO_LPO_N   = 20000    # prior label draws for the measured log_prior_odds

    # Deterministic, disjoint stream for CACHE PAIRING (which two cache columns form a
    # training pair + the train/val split). Distinct from VAL_MASTER_SEED / VAL_FIX_SEED /
    # NPE_MASTER_SEED and from the cache generation seed, so the assembled training pairs
    # never coincide with the reported BF eval stream (D-02 spirit).
    const RATIO_PAIR_SEED = 0x00000000004A7107

    # Default on-disk path for the persisted trained ratio net (large binary, gitignored).
    const RATIO_MODEL_PATH = joinpath(@__DIR__, "trained_ratio.jld2")
end

"""
    resolve_cache_dir() -> String

Locate the single 50k main-pool cache directory (the one the NPE trained on): the
`spike/data/cache/<hash>/` subdir that contains both `meta.jld2` and `shard_0001.jld2`
(the tiny `fixture/` cache and any stray dir are excluded by requiring a real shard +
meta). Errors if zero or more than one such directory exists, so the ratio net can never
silently train on the wrong pool.
"""
function resolve_cache_dir()
    root = joinpath(@__DIR__, "..", "data", "cache")
    isdir(root) || error("resolve_cache_dir: no cache root at $root")
    hits = String[]
    for d in readdir(root; join = true)
        isdir(d) || continue
        (isfile(joinpath(d, "meta.jld2")) && isfile(joinpath(d, "shard_0001.jld2"))) &&
            push!(hits, d)
    end
    length(hits) == 1 ||
        error("resolve_cache_dir: expected exactly one main-pool cache dir, found $(length(hits)): $hits")
    return hits[1]
end

"""
    assemble_ratio_data(m, n; rng = val_rng(), imsize = SBC_IMSIZE) -> (Z_pair, model_index)

Build `n` paired training columns for the model-comparison ratio net (BF-01, D-07).

For each pair: draw two INDEPENDENT priors θ_s, θ_c (`sample_prior`), simulate both
(`simulate_pair`), extract the FROZEN 8×8 patch-correlation summary and `encode_d01`
each to 128-dim, then standardize BOTH with the FROZEN `m.zt` (`standardize_summary`,
never re-fit). The column is `vcat(Z_s, Z_c)` (256-dim). The label is
`Float32((θ_s.ρ_true − θ_c.ρ_true) > RATIO_SPLIT_THRESHOLD)` ∈ {0f0,1f0} — the D-07
coloc/null split on the ρ_true contrast (≡ the baseline {Δμ>0} event, monotone ghat).

Returns `(Z_pair::Matrix{Float32} 256×n, model_index::Vector{Float32} length n)`.
CPU-only; consumes `rng` sequentially (a single disjoint stream, D-02).
"""
function assemble_ratio_data(m, n::Integer; rng = val_rng(), imsize = SBC_IMSIZE)
    Z_pair = Matrix{Float32}(undef, RATIO_INPUT_DIM, n)
    model_index = Vector{Float32}(undef, n)
    for j in 1:n
        θs = sample_prior(rng)
        θc = sample_prior(rng)
        Zs = standardize_summary(encode_d01(patch_summary(build_mci(simulate_pair(rng, θs; imsize = imsize)))), m.zt, :min)
        Zc = standardize_summary(encode_d01(patch_summary(build_mci(simulate_pair(rng, θc; imsize = imsize)))), m.zt, :min)
        Z_pair[:, j] = pair_encode(Zs, Zc)
        model_index[j] = Float32((θs.ρ_true - θc.ρ_true) > RATIO_SPLIT_THRESHOLD)
    end
    return Z_pair, model_index
end

"""
    pair_encode(Zs, Zc) -> Vector{Float32}

The paired-summary encoding for the model-comparison ratio net (research A7 difference
encoding): `vcat(Zs, Zc, Zs[1:RATIO_CONT_ROWS] − Zc[1:RATIO_CONT_ROWS])` — the two 128-dim
frozen-`m.zt` summaries concatenated PLUS the 64-dim continuous-row (correlation) CONTRAST.
The contrast row block is the direct per-patch Δcorrelation signal the D-07 coloc/null model
index is a function of. Used IDENTICALLY by the trainer (`assemble_ratio_data*`) and the
amortized BF read (`bf.jl::build_bf_pair`), so training and inference share one input space.
"""
function pair_encode(Zs::AbstractVector, Zc::AbstractVector)
    d = @view(Zs[1:RATIO_CONT_ROWS]) .- @view(Zc[1:RATIO_CONT_ROWS])
    return Float32.(vcat(Zs, Zc, d))
end

"""
    assemble_ratio_data_cache(m, n; cache_dir = resolve_cache_dir(),
                              rng = Philox4x(UInt64, (RATIO_PAIR_SEED, UInt64(0))),
                              variant = :min) -> (Z_pair, model_index)

Build `n` paired training columns for the model-comparison ratio net (BF-01, D-07) by
PAIRING the SAME 50k (θ, summary) cache the NPE trained on — NOT re-simulating (16× more
data, no simulation cost, and provably leak-free: the cache generation seed 0x134d8f3 is
disjoint from VAL_MASTER_SEED/NPE_MASTER_SEED, so these training pairs never coincide with
the reported BF eval stream). Loads the RAW main pool once, standardizes ALL summaries with
the FROZEN `m.zt` (`standardize_summary`, the deployed NPE input space, never re-fit), then
for each of `n` columns draws two DISTINCT cache indices `(i,j)` from `rng`, sets the column
to `vcat(Z_i, Z_j)` (256-dim) and the label to `Float32((ρ_i − ρ_j) > RATIO_SPLIT_THRESHOLD)`
— the D-07 coloc/null split on the ρ_true contrast (row 1 of θ). Because `(i,j)` are i.i.d.
over the pool, P(ρ_i > ρ_j) ≈ 0.5, so labels are balanced by construction.

Returns `(Z_pair::Matrix{Float32} 256×n, model_index::Vector{Float32} length n)`. CPU-only.
"""
function assemble_ratio_data_cache(m, n::Integer; cache_dir = resolve_cache_dir(),
                                   rng = Philox4x(UInt64, (RATIO_PAIR_SEED, UInt64(0))),
                                   variant::Symbol = :min)
    pool = Loader.load_main_pool(cache_dir)
    Zstd = standardize_summary(pool.summary_min, m.zt, variant)   # 128×N in the frozen space
    ρ    = Float64.(pool.theta[1, :])                             # ρ_true row (D-07 contrast)
    N    = size(Zstd, 2)
    N ≥ 2 || error("assemble_ratio_data_cache: need ≥2 cache columns, got $N")

    Z_pair = Matrix{Float32}(undef, RATIO_INPUT_DIM, n)
    model_index = Vector{Float32}(undef, n)
    for j in 1:n
        i1 = rand(rng, 1:N)
        i2 = rand(rng, 1:N)
        while i2 == i1
            i2 = rand(rng, 1:N)
        end
        Z_pair[:, j] = pair_encode(view(Zstd, :, i1), view(Zstd, :, i2))
        model_index[j] = Float32((ρ[i1] - ρ[i2]) > RATIO_SPLIT_THRESHOLD)
    end
    return Z_pair, model_index
end

"""
    measure_log_prior_odds(model_index) -> Float64

The MEASURED log prior-odds of the coloc event (Pitfall 5 — measured, NEVER assumed 0):
`log(p / (1 − p))` with `p = mean(model_index)` the prior fraction of coloc labels.
Because θ_s, θ_c are drawn i.i.d. from the SAME prior, P(ρ_s > ρ_c) = 0.5 by symmetry,
so this is ≈ 0 by construction — but it is MEASURED (not assumed) to absorb any
finite-sample/tie imbalance, per D-07. Returns a finite scalar.
"""
function measure_log_prior_odds(model_index::AbstractVector{<:Real})
    p = mean(model_index)
    (0.0 < p < 1.0) || error("measure_log_prior_odds: degenerate label balance p=$p (need 0<p<1)")
    return log(p / (1 - p))
end

"""
    prior_label_draws(n; rng = val_rng(), threshold = RATIO_SPLIT_THRESHOLD) -> Vector{Float32}

CHEAP model-index labels drawn from the prior WITHOUT simulating images: the D-07 label
depends only on the ρ_true contrast of two `sample_prior` draws, so the prior-odds can be
measured from `n` prior pairs at negligible cost. Used to measure `log_prior_odds` at a
large `n` when persisting the trained net.
"""
function prior_label_draws(n::Integer; rng = val_rng(), threshold = RATIO_SPLIT_THRESHOLD)
    labels = Vector{Float32}(undef, n)
    for j in 1:n
        ρs = sample_prior(rng).ρ_true
        ρc = sample_prior(rng).ρ_true
        labels[j] = Float32((ρs - ρc) > threshold)
    end
    return labels
end

"""
    train_ratio(m; n = RATIO_TRAIN_N, epochs = RATIO_EPOCHS, batchsize = 64,
                use_gpu = false, imsize = SBC_IMSIZE, val_frac = 0.2,
                rng = val_rng(), verbose = false) -> RatioEstimator

Train the model-comparison `RatioEstimator` CPU-only on `n` paired frozen-`m.zt`
summaries (BF-01). Builds the summary Chain `Dense(256,64,gelu)→Dense(64,64,gelu)→
Dense(64,RATIO_NUM_SUMMARIES)`, wraps it as `RatioEstimator(summary_network, 1;
num_summaries=RATIO_NUM_SUMMARIES)` (num_parameters = 1 = the binary model index), and
calls the FIXED-DATA `train(est, midx_tr, midx_va, Z_tr, Z_va; ...)` with AdamW weight
decay and early stopping. A custom loss is NEVER passed (it is silently ignored — the
loss is hard-coded `logitbinarycrossentropy`). `use_gpu=true` HARD-THROWS (CPU-only gate);
every NeuralEstimators call passes use_gpu=false.
"""
function train_ratio(m; n::Integer = RATIO_TRAIN_N, epochs::Integer = RATIO_EPOCHS,
                     batchsize::Integer = 128, use_gpu::Bool = false,
                     imsize = SBC_IMSIZE, val_frac::Real = 0.15,
                     learning_rate::Real = 2.5e-4, weight_decay::Real = 1e-4,
                     stopping_epochs::Integer = 40,
                     from_cache::Bool = true, cache_dir = nothing,
                     pair_rng = Philox4x(UInt64, (RATIO_PAIR_SEED, UInt64(0))),
                     rng = val_rng(), verbose::Bool = false)
    use_gpu && throw(ArgumentError("train_ratio: use_gpu=true out of scope (CPU-only gate, D-10)"))

    # PHASE-5 BF ITERATION: default to cache-paired data (leak-free, 16× the prior 3k
    # fresh-sim pairs). `from_cache=false` restores the original online-simulation path.
    Z_pair, model_index = from_cache ?
        assemble_ratio_data_cache(m, n;
            cache_dir = cache_dir === nothing ? resolve_cache_dir() : cache_dir,
            rng = pair_rng) :
        assemble_ratio_data(m, n; rng = rng, imsize = imsize)

    # Train/val split (contiguous — the stream is already i.i.d. across columns).
    nval = clamp(round(Int, val_frac * n), 1, n - 1)
    ntr  = n - nval
    Z_tr = Z_pair[:, 1:ntr]
    Z_va = Z_pair[:, ntr+1:end]
    midx_tr = reshape(model_index[1:ntr], 1, :)         # 1×ntr (num_parameters = 1)
    midx_va = reshape(model_index[ntr+1:end], 1, :)

    # PHASE-5 BF ITERATION: higher-capacity summary conditioner mirroring the NPE lever
    # (64-wide → 3×256-wide + a 64-dim learned-summary bottleneck).
    W = RATIO_SUMMARY_WIDTH
    summary_network = Chain(Dense(RATIO_INPUT_DIM, W, gelu),
                            Dense(W, W, gelu),
                            Dense(W, W, gelu),
                            Dense(W, RATIO_NUM_SUMMARIES))
    est = RatioEstimator(summary_network, 1; num_summaries = RATIO_NUM_SUMMARIES)

    # FIXED-DATA train form (mirrors train_npe.jl). AdamW args are Float64 to match
    # NeuralEstimators' Float64 CosAnneal lr_schedule (Optimisers.adjust! gotcha). Stabler
    # NPE-mirrored recipe: lower LR (2.5e-4), larger batch (128), 300 epochs / patience 40.
    est = train(est, midx_tr, midx_va, Z_tr, Z_va;
                epochs = epochs, batchsize = batchsize, use_gpu = false,
                optimiser = Flux.Optimisers.AdamW(learning_rate, (0.9, 0.999), weight_decay),
                stopping_epochs = stopping_epochs, verbose = verbose)
    return est
end

"""
    save_ratio(path, est, zt, log_prior_odds; meta = (;)) -> String

Atomically persist the trained ratio net + the FROZEN zt reference + the MEASURED
`log_prior_odds` (Pitfall 5): `jldsave` {schema_version, estimator, zt, log_prior_odds,
num_summaries, meta} to a `.tmp`, reopen to integrity-check the estimator key, then
`mv(...; force=true)` (the filesystem-atomic save_npe idiom). Returns `path`.
"""
function save_ratio(path, est, zt, log_prior_odds; meta = (;))
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp;
        schema_version = RATIO_MODEL_SCHEMA,
        estimator      = est,
        zt             = zt,
        log_prior_odds = log_prior_odds,
        num_summaries  = RATIO_NUM_SUMMARIES,
        split_threshold = RATIO_SPLIT_THRESHOLD,
        meta           = merge((generated = string(Dates.now(Dates.UTC)) * "Z",), meta))
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "estimator") "save_ratio: integrity check failed, $tmp missing estimator"
    end
    mv(tmp, path; force = true)   # atomic commit
    return path
end

"""
    load_ratio(path) -> NamedTuple

Reload a trained ratio net persisted by `save_ratio`. Returns
`(estimator, zt, log_prior_odds, num_summaries, meta)`. Integrity-checks the estimator
key on open (the frozen-stats + measured-prior-odds reproducibility contract bf.jl consumes).
"""
function load_ratio(path)
    isfile(path) || error("load_ratio: no model file at $path")
    d = JLD2.load(path)
    haskey(d, "estimator") || error("load_ratio: $path missing estimator key")
    return (estimator = d["estimator"], zt = d["zt"],
            log_prior_odds = d["log_prior_odds"],
            num_summaries = get(d, "num_summaries", RATIO_NUM_SUMMARIES),
            meta = get(d, "meta", nothing))
end

# --- Script entry point: persist a reported small training to trained_ratio.jld2 -----
# GUARDED so `include`-ing this file (bf.jl / test_bf.jl) never triggers training. Run
# `julia --project=spike spike/validation/train_ratio.jl` to (re)produce the net.
if abspath(PROGRAM_FILE) == @__FILE__
    m   = load_frozen_model()
    est = train_ratio(m; n = RATIO_TRAIN_N, epochs = RATIO_EPOCHS, verbose = true)
    lpo = measure_log_prior_odds(prior_label_draws(RATIO_LPO_N; rng = val_rng()))
    save_ratio(RATIO_MODEL_PATH, est, m.zt, lpo;
               meta = (n = RATIO_TRAIN_N, epochs = RATIO_EPOCHS,
                       num_summaries = RATIO_NUM_SUMMARIES,
                       split_threshold = RATIO_SPLIT_THRESHOLD))
    @info "train_ratio: persisted ratio net" path = RATIO_MODEL_PATH log_prior_odds = lpo
end
