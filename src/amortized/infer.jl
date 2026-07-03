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
# src/amortized/infer.jl --- the amortized NPE READ surface (PROD-01).
#
# Promoted (grid-generalized, CPU-default) from the proven spike read side:
#   spike/npe/infer.jl:76-84   standardize_summary  (frozen `zt`; mask rows bypassed)
#   spike/npe/infer.jl:96-138  posterior_for / rho_hat / rho_draws / delta_rho
#   spike/npe/infer.jl:140-151 interval_width
#
# This is the request-response service `colocalization_amortized` composes: a single CPU
# forward pass through the FROZEN `PosteriorEstimator` yields posterior draws, ρ_true draws,
# and the Δρ Monte-Carlo difference (D-03). The estimator is trained in STANDARDIZED θ-space,
# so EVERY ρ read first un-standardizes the draws with `StatsBase.reconstruct(θzt, draws)`
# before touching row 1 (Pitfall 5) — reading a standardized draw as ρ would be silently wrong.
#
# CPU-DEFAULT (D-06 reproducibility split): the shipped/reproducible read path defaults
# `use_gpu = false` on EVERY call (the pre-registered ship-gate numbers must be deterministic
# regardless of the hardware the net was TRAINED on). `use_gpu` stays a plumbed keyword so a
# caller can opt into GPU inference, but the shipped default is CPU.
#
# GRID-AGNOSTIC: the continuous/mask row split is delegated to `_summary_row_partition`
# (src/amortized/summary.jl) — the single grid-coupling site — so a new shipped grid flows
# through here with no change.
#############################################################################################

import NeuralEstimators: sampleposterior, posteriormean, interval
import StatsBase
import Statistics: mean

"""
    standardize_summary(S, zt, variant::Symbol = :min) -> AbstractVecOrMat{Float32}

Apply the FROZEN train-fit `ZScoreTransform` `zt` to a RAW summary (a `nrows`-vector for one
stack, or `nrows×K` for many) exactly as training did: the CONTINUOUS rows are z-scored with
`zt`, the binary MASK rows pass through UNCHANGED (never re-z-scored — z-scoring a 0/1 mask
would re-couple folds through the mask mean). This puts a holdout summary into the estimator's
input space with NO re-fitting (frozen-stats discipline). The row split is grid-general via
`_summary_row_partition` (summary.jl).
"""
function standardize_summary(S::AbstractMatrix, zt, variant::Symbol = :min)
    cont, mask = _summary_row_partition(variant, size(S, 1))
    out = Matrix{Float32}(undef, size(S))
    out[cont, :] = Float32.(StatsBase.transform(zt, S[cont, :]))
    out[mask, :] = Float32.(S[mask, :])
    return out
end
standardize_summary(v::AbstractVector, zt, variant::Symbol = :min) =
    vec(standardize_summary(reshape(v, :, 1), zt, variant))

"""
    posterior_for(est, Z; N::Integer = 2000, use_gpu::Bool = false) -> Matrix

Single-pass posterior for ONE stack: `Z` is a standardized summary VECTOR (or a single-column
`d_in×1` matrix). A vector is coerced to a `d_in×1` matrix so NeuralEstimators treats it as ONE
dataset (a bare vector is otherwise read as many 1-element datasets). Returns the `7×N` matrix
of posterior draws from one `sampleposterior` pass, in STANDARDIZED θ-space — un-standardize
with `θzt` before reading ρ. `use_gpu = false` by default (CPU-reproducible shipped path, D-06).
"""
function posterior_for(est, Z; N::Integer = 2000, use_gpu::Bool = false)
    Zc = Z isa AbstractVector ? reshape(Z, :, 1) : Z   # coerce to d_in×1 = ONE dataset
    return sampleposterior(est, Zc; N = N, use_gpu = use_gpu)
end

"""
    rho_hat(est, Z, θzt; N::Integer = 2000, use_gpu::Bool = false) -> Float64

ρ̂ for one stack: draw a single-pass posterior, UN-STANDARDIZE it with the frozen θ transform
`θzt` (`StatsBase.reconstruct`), then return the row-1 (ρ_true) posterior mean. On an
in-distribution stack the result lies in [-1, 1] (ρ is a correlation). CPU-default (D-06).
"""
function rho_hat(est, Z, θzt; N::Integer = 2000, use_gpu::Bool = false)
    draws = posterior_for(est, Z; N = N, use_gpu = use_gpu)   # 7×N standardized
    orig  = StatsBase.reconstruct(θzt, draws)                 # 7×N un-standardized (Pitfall 5)
    return posteriormean(orig)[1]                             # row 1 = ρ_true
end

"""
    rho_draws(est, Z, θzt; N::Integer = 2000, use_gpu::Bool = false) -> Vector

The `N` un-standardized ρ_true posterior draws for one stack (row 1 of the un-standardized
posterior). The building block of the Δρ MC difference. `StatsBase.reconstruct(θzt, ...)` is
applied BEFORE reading row 1 (Pitfall 5). CPU-default (D-06).
"""
function rho_draws(est, Z, θzt; N::Integer = 2000, use_gpu::Bool = false)
    draws = posterior_for(est, Z; N = N, use_gpu = use_gpu)
    orig  = StatsBase.reconstruct(θzt, draws)
    return vec(orig[1, :])                                    # row 1 = ρ_true
end

"""
    delta_rho(est, Zsample, Zcontrol, θzt; N::Integer = 2000, use_gpu::Bool = false) -> Float64

Δρ for a (sample, control) pair as the Monte-Carlo difference of TWO INDEPENDENT single-stack
posterior passes (D-03): `mean(ρ_sample_draws .- ρ_control_draws)`. Two independent forward
passes, not a joint estimator — the fully-amortized Δρ. CPU-default (D-06).

This is a 4-argument METHOD of the `delta_rho` accessor generic (results.jl also defines the
1-argument `delta_rho(::AbstractColocResult)` interface accessor); the two coexist by dispatch.
"""
function delta_rho(est, Zsample, Zcontrol, θzt; N::Integer = 2000, use_gpu::Bool = false)
    ρs = rho_draws(est, Zsample,  θzt; N = N, use_gpu = use_gpu)
    ρc = rho_draws(est, Zcontrol, θzt; N = N, use_gpu = use_gpu)
    return mean(ρs .- ρc)
end

"""
    interval_width(draws; probs = [0.05, 0.95]) -> Vector

Per-parameter credible-interval width from a `d×N` posterior-draw matrix, via the
NeuralEstimators `interval` (column 2 upper − column 1 lower). Pass UN-STANDARDIZED draws for
physical-unit widths; row 1 is the ρ_true interval width. Never hand-roll quantiles.
"""
function interval_width(draws; probs = [0.05, 0.95])
    I = interval(draws; probs = probs)   # d×2 (lower, upper)
    return collect(I[:, 2] .- I[:, 1])
end
