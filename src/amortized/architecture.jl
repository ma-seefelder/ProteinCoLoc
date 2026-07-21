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
# src/amortized/architecture.jl --- the NPE `PosteriorEstimator` topology (PROD-01).
#
# Promoted (input-width-agnostic, unchanged topology) from the proven spike:
#   spike/npe/architecture.jl:60-66   NPE_* architecture constants (Phase-5 calibrated)
#   spike/npe/architecture.jl:87-115  build_estimator (MLP summary net → NormalisingFlow)
#
# `build_estimator(d_in; ...)` reads its input width from `d_in = size(Ztr, 1)`, so it is
# ALREADY grid-general: a new shipped grid changes `d_in = 2·G²` but NOT this file. The
# Phase-5-calibrated higher-capacity defaults (dstar=64, depth=3/width=256, 10 coupling layers)
# are carried verbatim so the productionized estimator matches the spike-validated topology.
#
# API GOTCHA (preserved from the ONLY verified-working v0.2.1 construction): `q` is a
# CONSTRUCTED `NormalisingFlow` INSTANCE passed POSITIONALLY — the `q=` keyword convenience form
# of `PosteriorEstimator` expects a TYPE, not an instance. Do NOT "simplify" to `q=Type`.
#
# DEVICE-AGNOSTIC: the architecture itself carries no device; `use_gpu` is decided at the
# training call site (train_npe.jl, D-06). Persisted models are always CPU-resident (persist.jl).
#############################################################################################

import Flux: Chain, Dense, gelu
import NeuralEstimators: PosteriorEstimator, NormalisingFlow
import StatsBase

# =============================================================================================
# BOUNDED θ-SPACE (07-CALIBRATION-FINDINGS F2)
# =============================================================================================
#
# THE DEFECT: π(θ) is a hard-truncated BOX (e.g. `label_efficiency ~ Uniform(0.6, 1.0)`), but
# `NeuralEstimators.NormalisingFlow` models θ in an UNBOUNDED space after a plain z-score. It
# cannot represent the truncated support, and NPE's forward-KL objective constrains nothing about
# the implied marginal ∫q(θ|Z)p(Z)dZ = π(θ). Measured consequence (grid 8, `label_efficiency`):
# SBC mean rank 0.548 (z = +7.4) with ~2% of posterior mass leaking outside [0.6, 1.0], the
# leakage sign-matched to the bias.
#
# THE FIX: interpose a per-parameter LOGIT bijection from the prior box onto ℝ BEFORE the
# z-score, and invert it on read-out. The flow then works in a space where the whole real line is
# legal, and every posterior draw maps back STRICTLY INSIDE the prior box — out-of-support mass is
# 0 BY CONSTRUCTION, not by learning.
#
#   forward :  u = (θ − lo)/(hi − lo) ∈ [0,1]  →  y = logit(clamp(u, ε, 1−ε))  →  z-score
#   inverse :  un-z-score  →  θ = lo + (hi − lo)·logistic(y)                    ∈ [lo, hi]
#
# The map is STRICTLY MONOTONE per component, so SBC ranks are identical whether computed in
# θ-space or in flow-space — the fix changes what the flow must represent, not what is measured.
#
# BOUNDS PROVENANCE: `lo`/`hi` come from `theta_prior_bounds()` (simulator.jl), which DERIVES them
# from the prior objects `sample_prior` draws from. No bound is ever hardcoded here.
#
# ε GUARD: ρ_true = `ghat(μ*)` is a CLAMPED map, so it puts genuine POINT MASSES on both endpoints
# (μ* outside the swept knot range). `logit` of an exact endpoint is ±Inf, so u is clamped into
# [ε, 1−ε] — the atoms become the large-but-finite values ±logit(ε). ε is a numerical guard only;
# it never widens the support (the INVERSE is exact and always lands in [lo, hi]).
#
# DISPATCH-COMPATIBLE: `StatsBase.transform`/`StatsBase.reconstruct` methods are defined for this
# type, so EVERY existing call site (`infer.jl` ρ reads, `ood.jl`, `local_map.jl`, the gate
# harness) is unchanged — they already call `reconstruct(θzt, draws)` generically. A previously
# frozen artifact carrying a plain `ZScoreTransform` keeps working through StatsBase's own methods.

"""
    THETA_LOGIT_EPS

Numerical clamp for the logit forward map (`u ∈ [ε, 1−ε]`). Guards the prior's endpoint ATOMS
(ρ_true = `ghat(μ*)` clamps) against `logit(0) = -Inf`. It bounds the magnitude of the mapped
atoms (`|y| ≤ logit(ε) ≈ 13.8`); it does NOT widen the support, because the inverse map is exact.
"""
const THETA_LOGIT_EPS = 1e-6

"""
    BoundedThetaTransform

The frozen θ transform for a bounded-support prior: a per-parameter logit bijection from the prior
box `[lo, hi]` onto ℝ, composed with a `ZScoreTransform` fitted in that unconstrained space.

Fields: `lo`/`hi` (the prior box, from `theta_prior_bounds()`) and `zt` (the leak-free
train-fitted `ZScoreTransform` in logit space).

Implements `StatsBase.transform` (θ → flow space) and `StatsBase.reconstruct` (flow space → θ),
so it is a DROP-IN for the plain `ZScoreTransform` at every existing call site. `reconstruct`
output is guaranteed to lie in `[lo, hi]`.
"""
struct BoundedThetaTransform{T<:Real,Z}
    lo::Vector{T}
    hi::Vector{T}
    zt::Z
end

_logit_eps(u::Real, ε::Real) = (c = clamp(u, ε, one(u) - ε); log(c / (one(c) - c)))
_logistic(y::Real) = inv(one(y) + exp(-y))

"""
    theta_to_unbounded(θ, lo, hi; eps = THETA_LOGIT_EPS) -> Matrix

Per-row logit map of a `d×n` θ matrix from the prior box onto ℝ. Rows are parameters. Inputs are
`clamp`ed into `[lo, hi]` first, so a caller that hands over synthetic (out-of-box) θ degrades to
the boundary rather than producing `NaN`.
"""
function theta_to_unbounded(θ::AbstractMatrix, lo::AbstractVector, hi::AbstractVector;
                            eps::Real = THETA_LOGIT_EPS)
    size(θ, 1) == length(lo) == length(hi) ||
        throw(DimensionMismatch("theta_to_unbounded: θ rows $(size(θ,1)) vs bounds $(length(lo))"))
    out = Matrix{Float64}(undef, size(θ))
    @inbounds for i in axes(θ, 1)
        w = hi[i] - lo[i]
        w > 0 || throw(ArgumentError("theta_to_unbounded: degenerate bound on row $i"))
        for j in axes(θ, 2)
            out[i, j] = _logit_eps((clamp(θ[i, j], lo[i], hi[i]) - lo[i]) / w, eps)
        end
    end
    return out
end

"""
    theta_from_unbounded(Y, lo, hi) -> Matrix

Inverse of `theta_to_unbounded`: `θ = lo + (hi − lo)·logistic(y)`, row-wise. The result is
GUARANTEED inside `[lo, hi]` (a final `clamp` absorbs float round-off at the extremes), which is
the whole point — no posterior draw can leave the prior support.
"""
function theta_from_unbounded(Y::AbstractMatrix, lo::AbstractVector, hi::AbstractVector)
    size(Y, 1) == length(lo) == length(hi) ||
        throw(DimensionMismatch("theta_from_unbounded: rows $(size(Y,1)) vs bounds $(length(lo))"))
    out = similar(Y, float(eltype(Y)))
    @inbounds for i in axes(Y, 1)
        w = hi[i] - lo[i]
        for j in axes(Y, 2)
            out[i, j] = clamp(lo[i] + w * _logistic(Y[i, j]), lo[i], hi[i])
        end
    end
    return out
end

function StatsBase.transform(t::BoundedThetaTransform, θ::AbstractMatrix)
    return StatsBase.transform(t.zt, theta_to_unbounded(θ, t.lo, t.hi))
end
StatsBase.transform(t::BoundedThetaTransform, v::AbstractVector) =
    vec(StatsBase.transform(t, reshape(v, :, 1)))

function StatsBase.reconstruct(t::BoundedThetaTransform, Y::AbstractMatrix)
    return theta_from_unbounded(StatsBase.reconstruct(t.zt, Y), t.lo, t.hi)
end
StatsBase.reconstruct(t::BoundedThetaTransform, v::AbstractVector) =
    vec(StatsBase.reconstruct(t, reshape(v, :, 1)))

# --- Architecture constants (Phase-5 calibrated; every consumer reuses these) ---------------
# The frozen Phase-4 topology (dstar=32, depth=2/width=128, 6 coupling layers) produced
# OVER-DISPERSED ρ_true/Δρ posteriors that failed the Phase-5 SBC gates; the flow lacked
# capacity to sharpen the conditionals. These raised defaults give the SAME `build_estimator`
# a higher-capacity flow. ρ_true stays θ ROW 1 and D stays 7, so downstream index/axis
# assumptions are unchanged.
const NPE_D          = 7    # θ dimension (ρ_true is ROW 1; prior field order)
const NPE_DSTAR      = 64   # learned-summary dimension feeding the flow conditioner
const NPE_DEPTH      = 3    # hidden Dense layers in the MLP summary net
const NPE_WIDTH      = 256  # hidden width of the MLP summary net
const NPE_COUPLING   = 10   # NormalisingFlow coupling layers (was library default 6)
const NPE_FLOW_DEPTH = 2    # conditioner-MLP hidden layers per affine coupling block
const NPE_FLOW_WIDTH = 128  # conditioner-MLP hidden width per affine coupling block

"""
    build_estimator(d_in::Integer, D::Integer = NPE_D;
                    dstar::Integer = NPE_DSTAR, depth::Integer = NPE_DEPTH,
                    width::Integer = NPE_WIDTH,
                    num_coupling_layers::Integer = NPE_COUPLING,
                    flow_depth::Integer = NPE_FLOW_DEPTH,
                    flow_width::Integer = NPE_FLOW_WIDTH) -> PosteriorEstimator

Build the amortized-NPE `PosteriorEstimator`: an MLP summary net mapping a `d_in`-vector
standardized summary (`d_in = 2·G²` for a `G×G` grid) to `dstar` learned summaries, feeding a
NeuralEstimators `NormalisingFlow` over the `D`-dim (default 7) simulator θ. The flow therefore
has `D` marginals; marginal 1 is ρ_true.

`q` is constructed as a `NormalisingFlow` INSTANCE and passed POSITIONALLY (the v0.2.1 API
gotcha). `num_coupling_layers` sets the flow's coupling-stack depth; `flow_depth`/`flow_width`
are forwarded to the per-block conditioner MLPs. All knobs default to the module `const`s.
`build_estimator` is INPUT-WIDTH-AGNOSTIC (reads `d_in`), so it needs no per-grid change.
"""
function build_estimator(d_in::Integer, D::Integer = NPE_D;
                         dstar::Integer = NPE_DSTAR, depth::Integer = NPE_DEPTH,
                         width::Integer = NPE_WIDTH,
                         num_coupling_layers::Integer = NPE_COUPLING,
                         flow_depth::Integer = NPE_FLOW_DEPTH,
                         flow_width::Integer = NPE_FLOW_WIDTH)
    d_in  >= 1 || throw(ArgumentError("build_estimator: d_in must be ≥ 1, got $d_in"))
    D     >= 1 || throw(ArgumentError("build_estimator: D must be ≥ 1, got $D"))
    dstar >= D || throw(ArgumentError("build_estimator: dstar ($dstar) must be ≥ D ($D)"))
    depth >= 1 || throw(ArgumentError("build_estimator: depth must be ≥ 1, got $depth"))
    num_coupling_layers >= 1 ||
        throw(ArgumentError("build_estimator: num_coupling_layers must be ≥ 1"))

    # MLP summary net: Dense(d_in→width, gelu), (depth-1) width→width gelu blocks, then a
    # linear width→dstar projection to the learned summaries. Concretely typed (not Any[...]).
    layers = Dense[Dense(d_in, width, gelu)]
    for _ in 2:depth
        push!(layers, Dense(width, width, gelu))
    end
    push!(layers, Dense(width, dstar))
    network = Chain(layers...)

    # NormalisingFlow INSTANCE, passed POSITIONALLY. `num_coupling_layers` deepens the coupling
    # stack; `depth`/`width` widen the per-block conditioner nets (the capacity that sharpens
    # the flow).
    q = NormalisingFlow(D; num_summaries = dstar,
                        num_coupling_layers = num_coupling_layers,
                        depth = flow_depth, width = flow_width)
    return PosteriorEstimator(network, q)
end
