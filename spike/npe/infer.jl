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

# spike/npe/infer.jl --- Phase-4 NPE inference surface (NPE-01/02, D-03, RESEARCH Pattern 2).
#
# The amortized read side of the trained PosteriorEstimator: single-pass posteriors,
# ρ̂ (row-1 posterior mean, un-standardized), the Δρ Monte-Carlo difference (D-03),
# NeuralEstimators interval widths, and the μ→ρ `ghat` mapping the benchmark (04-05)
# uses to put the ADVI baseline on the NPE's known ρ_true axis.
#
# CENTRAL BENCHMARK-AXIS ASSUMPTION (RESEARCH A1 / Pattern 2): the NPE infers the
# SIMULATOR θ (ρ_true, the known coloc ground truth); the Turing/ADVI baseline infers
# the INDUCED per-condition mean μ. They live in different spaces. The whole
# NPE-vs-ADVI comparison is done in ρ-SPACE by pushing ADVI's μ posterior through the
# FROZEN `ghat` (μ→ρ) via `rho_from_mu`, so both methods score against each stack's
# known ρ_true. `ghat` is auto-generated/frozen (SIM-02); it is CALLED, never recomputed.
#
# θ-UN-STANDARDIZATION (Pitfall 5): the estimator is trained in standardized θ-space
# (train_npe.jl fits `θzt` on θtr only). EVERY ρ read here first un-standardizes the
# draws with `StatsBase.reconstruct(θzt, draws)` before touching row 1 -- reading a
# standardized draw as ρ would be silently wrong.
#
# CPU-ONLY (D-10 / Pitfall 1): `use_gpu = false` is passed on every `sampleposterior`
# call. Per-parameter RMSE, when needed downstream (04-05), MUST go through
# NeuralEstimators `assess`/`rmse` (RESEARCH "Don't Hand-Roll") -- never a manual
# error loop -- and be scaled back to physical units by `θzt.scale`.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; reaches src/ only transitively
# through train_npe.jl's guarded chain. Guarded includes keep the file loadable
# standalone AND inside runtests.jl.

using NeuralEstimators   # sampleposterior, posteriormean, interval
using StatsBase          # reconstruct (inverse ZScoreTransform)
using Statistics         # mean

# --- ORDER MATTERS: trainer first (build_estimator/train_fold/load_npe + loader +
#     architecture), then the FROZEN ghat μ→ρ map. Guarded for idempotency.
isdefined(@__MODULE__, :train_fold) || include(joinpath(@__DIR__, "train_npe.jl"))
isdefined(@__MODULE__, :ghat)       || include(joinpath(@__DIR__, "..", "simulator", "ghat.jl"))

"""
    _summary_row_partition(variant::Symbol, nrows::Int) -> (cont_rows, mask_rows)

The SAME continuous/mask row split the loader uses (delegates to the authoritative
`Loader._row_partition`) so a holdout summary is standardized byte-identically to the
training Z: continuous rows z-scored, binary mask rows (65:128) passed through.
"""
_summary_row_partition(variant::Symbol, nrows::Int) = Loader._row_partition(variant, nrows)

"""
    standardize_summary(S, zt, variant::Symbol = :min) -> AbstractVecOrMat{Float32}

Apply the FROZEN loader train-fit `ZScoreTransform` `zt` to a RAW summary (a
`nrows`-vector for one stack, or `nrows×K` for many) exactly as `load_fold` did:
continuous rows are z-scored with `zt`, the binary mask rows (65:128) pass through
UNCHANGED (never re-z-scored; Phase-3 D-08 / Anti-Pattern). This is how a holdout
`summary_min` is put into the estimator's input space with no re-fitting.
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

Single-pass posterior for ONE stack: `Z` is a standardized summary VECTOR (or a
single-column `d_in×1` matrix). A vector is coerced to a `d_in×1` matrix so
NeuralEstimators treats it as ONE dataset (a bare vector is otherwise read as many
1-element datasets). Returns the `7×N` matrix of posterior draws from one
`sampleposterior` pass (STANDARDIZED θ-space -- un-standardize before reading ρ).
`use_gpu=false` on every call (D-10).
"""
function posterior_for(est, Z; N::Integer = 2000, use_gpu::Bool = false)
    Zc = Z isa AbstractVector ? reshape(Z, :, 1) : Z
    return sampleposterior(est, Zc; N = N, use_gpu = use_gpu)
end

"""
    rho_hat(est, Z, θzt; N::Integer = 2000, use_gpu::Bool = false) -> Float64

ρ̂ for one stack: draw a single-pass posterior, UN-STANDARDIZE it with the frozen θ
transform `θzt` (`StatsBase.reconstruct`), then return the row-1 (ρ_true) posterior
mean. On an in-distribution stack the result lies in [-1,1] (ρ is a correlation).
"""
function rho_hat(est, Z, θzt; N::Integer = 2000, use_gpu::Bool = false)
    draws = posterior_for(est, Z; N = N, use_gpu = use_gpu)   # 7×N standardized
    orig  = StatsBase.reconstruct(θzt, draws)                 # 7×N un-standardized
    return posteriormean(orig)[1]                             # row 1 = ρ_true
end

"""
    rho_draws(est, Z, θzt; N, use_gpu) -> Vector

The N un-standardized ρ_true posterior draws for one stack (row 1 of the
un-standardized posterior). The building block of the Δρ MC difference.
"""
function rho_draws(est, Z, θzt; N::Integer = 2000, use_gpu::Bool = false)
    draws = posterior_for(est, Z; N = N, use_gpu = use_gpu)
    orig  = StatsBase.reconstruct(θzt, draws)
    return vec(orig[1, :])
end

"""
    delta_rho(est, Zsample, Zcontrol, θzt; N::Integer = 2000, use_gpu = false) -> Float64

Δρ for a (sample, control) pair as the Monte-Carlo difference of TWO INDEPENDENT
single-stack posterior passes (D-03): `mean(ρ_sample_draws .- ρ_control_draws)`.
Two independent forward passes, not a joint estimator -- the fully-amortized Δρ the
benchmark (04-05) times against the ADVI joint `vi()` run (D-04).
"""
function delta_rho(est, Zsample, Zcontrol, θzt; N::Integer = 2000, use_gpu::Bool = false)
    ρs = rho_draws(est, Zsample,  θzt; N = N, use_gpu = use_gpu)
    ρc = rho_draws(est, Zcontrol, θzt; N = N, use_gpu = use_gpu)
    return mean(ρs .- ρc)
end

"""
    interval_width(draws; probs = [0.05, 0.95]) -> Vector

Per-parameter credible-interval width from a `d×N` posterior-draw matrix, via the
NeuralEstimators `interval` (column 2 upper − column 1 lower). Pass UN-STANDARDIZED
draws for physical-unit widths. Row 1 is the ρ_true interval width. Never hand-roll
quantiles (RESEARCH "Don't Hand-Roll").
"""
function interval_width(draws; probs = [0.05, 0.95])
    I = interval(draws; probs = probs)   # NamedArray d×2 (lower, upper)
    return collect(I[:, 2] .- I[:, 1])
end

"""
    rho_from_mu(μ_draws) -> Float64

Map an ADVI μ posterior onto the NPE's known ρ_true axis: `mean(ghat.(μ_draws))`
(RESEARCH Pattern 2 / A1). `ghat` is the FROZEN monotone induced-μ inverse
(SIM-02) -- called per draw, NEVER recomputed. This is the cross-method common axis
the 04-05 benchmark scores both methods on.
"""
rho_from_mu(μ_draws) = mean(ghat.(μ_draws))
