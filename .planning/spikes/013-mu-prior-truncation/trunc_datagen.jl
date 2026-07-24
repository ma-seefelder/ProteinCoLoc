#########################################################################################
# Spike 013 --- the SIMULATOR-side μ-prior TRUNCATION lever (shared datagen).
#
# QUESTION (Spike 010 addendum). `ProteinCoLoc.sample_prior` draws μ* ~ Truncated(Cauchy(0,0.3),
# -1, 1) (the Turing μ-prior) and sets ρ_true = ghat(μ*). `ghat` is a CLAMPED piecewise-linear
# inverse over the physically-realized sweep range μ ∈ [GHAT_MU_MIN, GHAT_MU_MAX] =
# [-0.67976, 0.847149]; Cauchy tails past that map to EXACTLY ±0.99 (~6.5% of draws = prior atoms).
# SBC assumes a continuous prior, so those atoms make the ρ_true rank test invalid without a
# randomized-rank correction.
#
# THE LEVER. TRUNCATE the μ-prior to the ghat-achievable support [GHAT_MU_MIN, GHAT_MU_MAX]. Then
# ρ_true = ghat(μ*) is CONTINUOUS over [-0.99, 0.99] (the endpoints are hit only in the measure-zero
# limit), so there are NO clamp atoms and ρ_true SBC is valid WITHOUT randomized ranks.
#
# COST (the reason this is a spike, not a fix). This DEVIATES from the CLAUDE.md constraint that the
# simulator prior match the Turing μ-prior for ADVI comparability, and it removes ALL training mass
# beyond μ ∈ [-0.68, 0.847] — i.e. the ±0.99 correlation regime the clamp used to populate. The
# spike exists to MEASURE whether removing that mass degrades high-|ρ| inference.
#
# DECOUPLING. src/ is NOT edited. This file DEFINES a truncated-prior variant of `sample_prior` /
# `generate_sample` / `generate_samples` by REUSING every src/ constant and physics routine
# (SPILLOVER_PRIOR, ..., ghat, sample_rng, sample_imsize, simulate_pair, build_mci, patch_summary,
# encode_d01, summary_dim). The ONLY change vs the src/ path is which μ-prior μ* is drawn from.
# It is injected into `_train_grid_pipeline` via the `datagen` seam — src/ stays byte-untouched.
#########################################################################################
using ProteinCoLoc
import Distributions: Truncated, Cauchy
import Random: AbstractRNG

# The truncated μ-prior: the SAME Cauchy(0,0.3) location/scale as ProteinCoLoc.MU_PRIOR, truncated
# to the ghat-achievable support instead of [-1, 1]. ρ_true = ghat(μ*) is then continuous, no atoms.
const MU_PRIOR_TRUNC = Truncated(Cauchy(0.0, 0.3),
                                 ProteinCoLoc.GHAT_MU_MIN, ProteinCoLoc.GHAT_MU_MAX)

"""
    sample_prior_trunc(rng) -> NamedTuple

Identical to `ProteinCoLoc.sample_prior` EXCEPT μ* is drawn from `MU_PRIOR_TRUNC` (the
ghat-achievable support) so `ρ_true = ghat(μ*)` is CONTINUOUS over [-0.99, 0.99] — no ±0.99 clamp
atoms. The 6 nuisances are drawn from the SAME src/ prior constants (`ProteinCoLoc.SPILLOVER_PRIOR`
etc — NOT redefined here) in the SAME field order, so the per-index rng draw STRUCTURE matches
`sample_prior` exactly (7 draws before the imsize draw).
"""
function sample_prior_trunc(rng::AbstractRNG)
    μ_star = rand(rng, MU_PRIOR_TRUNC)
    return (
        ρ_true           = ProteinCoLoc.ghat(μ_star),
        spillover        = rand(rng, ProteinCoLoc.SPILLOVER_PRIOR),
        autofluorescence = rand(rng, ProteinCoLoc.AUTOFLUORESCENCE_PRIOR),
        label_efficiency = rand(rng, ProteinCoLoc.LABEL_EFFICIENCY_PRIOR),
        shift_dx         = rand(rng, ProteinCoLoc.SHIFT_PRIOR),
        shift_dy         = rand(rng, ProteinCoLoc.SHIFT_PRIOR),
        noise            = rand(rng, ProteinCoLoc.NOISE_PRIOR),
    )
end

"""
    generate_sample_trunc(master_seed, idx, grid; imsize_set, imsize_weights) -> NamedTuple

Truncated-prior mirror of `ProteinCoLoc.generate_sample`: identical Philox-per-index keying
(`ProteinCoLoc.sample_rng(master_seed, idx)`), identical imsize draw, forward model, patch summary
and D-01 encoding — the ONLY change is `sample_prior_trunc` in place of `sample_prior`.
"""
function generate_sample_trunc(master_seed::Integer, idx::Integer, grid::Integer;
                               imsize_set, imsize_weights)
    rng = ProteinCoLoc.sample_rng(master_seed, idx)
    θ   = sample_prior_trunc(rng)
    isz = ProteinCoLoc.sample_imsize(rng; imsize_set = imsize_set, imsize_weights = imsize_weights)
    mci = ProteinCoLoc.build_mci(ProteinCoLoc.simulate_pair(rng, θ; imsize = isz))
    M   = ProteinCoLoc.patch_summary(mci, grid)
    return (theta = collect(values(θ)), s_min = ProteinCoLoc.encode_d01(M),
            idx = Int(idx), imsize = isz)
end

"""
    generate_samples_trunc(N; grid, master_seed, indices, imsize_set, imsize_weights, parallel)
        -> (; theta, summary_min, global_index, imsize)

Truncated-prior mirror of `ProteinCoLoc.generate_samples`: column-major pre-allocated buffers, one
sample per column, each a pure function of `(master_seed, grid, indices[j])` so the parallel and
serial paths are byte-identical. Returned shape `(; theta::7×N, summary_min::summary_dim(grid)×N,
...)` — exactly what the `_train_grid_pipeline` `datagen` seam consumes.
"""
function generate_samples_trunc(N::Integer; grid::Integer, master_seed::Integer,
                                indices = 1:N, imsize_set, imsize_weights,
                                parallel::Bool = Threads.nthreads() > 1)
    idxv = collect(indices)
    @assert length(idxv) == N "indices length $(length(idxv)) != N=$N"
    d            = ProteinCoLoc.summary_dim(grid)
    theta        = Matrix{Float64}(undef, 7, N)
    summary_min  = Matrix{Float64}(undef, d, N)
    global_index = Vector{Int}(undef, N)
    imsize       = Vector{Tuple{Int,Int}}(undef, N)
    fill_col! = function (j::Int)
        s = generate_sample_trunc(master_seed, idxv[j], grid;
                                  imsize_set = imsize_set, imsize_weights = imsize_weights)
        @inbounds theta[:, j]       = s.theta
        @inbounds summary_min[:, j] = s.s_min
        @inbounds global_index[j]   = s.idx
        @inbounds imsize[j]         = s.imsize
        return nothing
    end
    if parallel
        Threads.@threads for j in 1:N
            fill_col!(j)
        end
    else
        for j in 1:N
            fill_col!(j)
        end
    end
    return (theta = theta, summary_min = summary_min,
            global_index = global_index, imsize = imsize)
end
