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

# spike/npe/p11_architecture.jl --- Phase-11 research-net model surface (D-02, D-03, ruling Q4).
#
# DECLARED DEVIATIONS (D-04). The Phase-11 research net differs from the shipped 8x8 net in
# EXACTLY FOUR ways, and those four are the ones pre-registered in
# `P11_RESEARCH_NET_DEVIATIONS` (spike/validation/p11_consts.jl, Tier 1, committed before
# anything ran):
#   (1) chromatic_eps appended as the 8th theta column (D-09), so D = 8 rather than 7;
#   (2) SHIFT_PRIOR widened to Uniform(-3, 3) px -- SPIKE ONLY, deliberately NOT mirrored
#       into src/ (D-02);
#   (3) lambda appended as a 129th input row, so d_in = 129 rather than 128 (D-03);
#   (4) BoundedThetaTransform ported from src/amortized/architecture.jl into the spike
#       trainer (the 07-CALIBRATION-FINDINGS F2 remedy).
# THIS FILE IMPLEMENTS TWO OF THEM -- (4) the ported theta transform and (3) the 129th lambda
# row. The other two live in `spike/simulator/prior.jl` ((1) chromatic_eps as an 8th prior
# field, (2) the widened SHIFT_PRIOR) and are consumed here only through
# `p11_theta_prior_bounds()`; the lambda-conditional shift draw that feeds them is
# `spike/data/p11_generate.jl`'s job, not this file's.
#
# CONSEQUENCE, STATED PLAINLY: BECAUSE OF (1) AND (3) THIS RESEARCH NET IS **NOT DROP-IN
# COMPARABLE TO THE SHIPPED READ SURFACE**. The shipped `EstimatorBundle` reads a 128-row input
# and a 7-marginal flow; nothing in this file may propagate into it. `NPE_D` stays 7 and
# `build_estimator`'s default D stays `NPE_D` (D-16: no public API change, no artifact change).
# Any number produced through this file must be attributed to the RESEARCH net by name.
#
# NAME COLLISION WARNING. `P11_THETA_LOGIT_EPS` (the logit clamp guard, ported from src) and
# D-09's `chromatic_eps` (a theta field, a radial magnification difference between channels)
# are DIFFERENT THINGS that share a Greek letter in the literature. No bare single-letter
# epsilon identifier is written anywhere in this file; every one is spelled out.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local; reaches src/ only READ-ONLY -- the
# bounded theta transform below is a COPY of src bytes with attribution, not a dependency, and
# nothing here writes to or imports from the ProteinCoLoc module. Guarded includes keep the
# file loadable standalone AND inside a test harness.

using StatsBase          # ZScoreTransform, fit, transform, reconstruct
using Distributions      # minimum/maximum of the spike prior objects (bounds derivation)

# --- ORDER MATTERS: the Tier-1 pre-registration first (LAMBDA_MIN/LAMBDA_MAX and the declared
#     deviation list), then the estimator builder (build_estimator + the NPE_* consts), then the
#     simulator prior (the objects `p11_theta_prior_bounds()` DERIVES its box from). Guarded for
#     idempotency: a re-include under a test file that already loaded them is a silent no-op.
isdefined(@__MODULE__, :LAMBDA_MIN)      || include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
isdefined(@__MODULE__, :build_estimator) || include(joinpath(@__DIR__, "architecture.jl"))
isdefined(@__MODULE__, :sample_prior)    || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))

# =============================================================================================
# 1. BOUNDED theta-SPACE (07-CALIBRATION-FINDINGS F2) -- PORTED FROM src/
# =============================================================================================
#
# ATTRIBUTION. The struct, its two `StatsBase` methods and the rationale block below are COPIED
# from `src/amortized/architecture.jl:43-104,157-159` (and the fit-side construction from
# `src/amortized/train_npe.jl:66-77`). They are COPIED RATHER THAN DEPENDED ON because the spike
# environment does not load the `ProteinCoLoc` module -- `spike/Project.toml` has no path
# dependency on the root package, and adding one would violate the CLAUDE.md decoupling
# constraint and the D-01 no-new-dependency rule. The port is AUTHORISED BY RULING Q4 and is
# entry (4) of the pre-registered `P11_RESEARCH_NET_DEVIATIONS`.
#
# WHY THE PORT IS NECESSARY AND NOT A CONVENIENCE. `spike/npe/train_npe.jl:80` fits a plain
# `ZScoreTransform` on raw theta. D-02 requires the rho_true attenuation under the widened
# nuisance prior to be MEASURED against the frozen shipped net -- and the shipped net HAS the
# bounded transform. Training the research net without it would confound "the widened nuisance
# prior attenuates rho_true" (the quantity D-02 asks for) with "a plain z-score leaks ~2% of
# posterior mass outside the prior box" (a known, already-diagnosed defect), on exactly the one
# number the comparison reports.
#
# THE DEFECT (verbatim rationale from the src source). pi(theta) is a hard-truncated BOX (e.g.
# `label_efficiency ~ Uniform(0.6, 1.0)`), but `NeuralEstimators.NormalisingFlow` models theta in
# an UNBOUNDED space after a plain z-score. It cannot represent the truncated support, and NPE's
# forward-KL objective constrains nothing about the implied marginal
# int q(theta|Z) p(Z) dZ = pi(theta). Measured consequence (grid 8, `label_efficiency`): SBC mean
# rank 0.548 (z = +7.4) with ~2% of posterior mass leaking outside [0.6, 1.0], the leakage
# sign-matched to the bias.
#
# THE FIX. Interpose a per-parameter LOGIT bijection from the prior box onto R BEFORE the
# z-score, and invert it on read-out. The flow then works in a space where the whole real line is
# legal, and every posterior draw maps back STRICTLY INSIDE the prior box -- out-of-support mass
# is 0 BY CONSTRUCTION, not by learning.
#
#   forward :  u = (theta - lo)/(hi - lo) in [0,1]
#              -> y = logit(clamp(u, P11_THETA_LOGIT_EPS, 1 - P11_THETA_LOGIT_EPS))
#              -> z-score
#   inverse :  un-z-score  ->  theta = lo + (hi - lo) * logistic(y)   in [lo, hi]
#
# The map is STRICTLY MONOTONE per component, so SBC ranks are identical whether computed in
# theta-space or in flow-space -- the fix changes what the flow must represent, not what is
# measured.
#
# BOUNDS PROVENANCE: `lo`/`hi` come from `p11_theta_prior_bounds()` below, which DERIVES them
# from the spike prior objects `sample_prior` draws from. No bound is ever hardcoded here.
#
# LOGIT-CLAMP GUARD, AND THE ONE PIECE OF src PROSE THAT DOES **NOT** CARRY OVER. In src the
# clamp is justified by rho_true = `ghat(mu*)` being a CLAMPED map, which puts genuine POINT
# MASSES (atoms) on both endpoints of ROW 1; `logit` of an exact endpoint is +-Inf, so u is
# clamped and the atoms become the large-but-finite values +-logit(clamp). That argument applies
# to ROW 1 ONLY and is retained for it. It does NOT apply to the new 8th row: `chromatic_eps` is
# drawn straight from `CHROMATIC_PRIOR = Uniform(-0.02, 0.02)` with no clamped map in front of
# it, so THE chromatic_eps COLUMN HAS NO ATOMS and the clamp there is purely a numerical guard
# against float round-off at the box edge. Do not read the atom prose onto it.
# In neither case does the clamp widen the support: the INVERSE is exact and always lands in
# [lo, hi].
#
# DISPATCH-COMPATIBLE: `StatsBase.transform`/`StatsBase.reconstruct` methods are defined for the
# ported type, so every existing spike call site (`infer.jl`'s rho reads, which already call
# `reconstruct(theta_zt, draws)` generically) is unchanged.

"""
    P11_THETA_LOGIT_EPS

Numerical clamp for the logit forward map (`u` is confined to `[c, 1-c]`). Guards ROW 1's
endpoint ATOMS (`rho_true = ghat(mu*)` clamps) against `logit(0) = -Inf`, and guards every other
row against float round-off exactly at the box edge. It bounds the magnitude of the mapped atoms
(`|y| <= logit(c) ~ 13.8`); it does NOT widen the support, because the inverse map is exact.

Value copied from `src/amortized/architecture.jl:85` (`THETA_LOGIT_EPS`) so the ported transform
is numerically identical to the shipped one; the spike name is prefixed so the two can never be
confused in a stack trace.
"""
const P11_THETA_LOGIT_EPS = 1e-6

"""
    P11BoundedThetaTransform

The frozen theta transform for a bounded-support prior: a per-parameter logit bijection from the
prior box `[lo, hi]` onto R, composed with a `ZScoreTransform` fitted in that unconstrained
space.

Fields: `lo`/`hi` (the prior box, from `p11_theta_prior_bounds()`) and `zt` (the leak-free
train-fitted `ZScoreTransform` in logit space).

Implements `StatsBase.transform` (theta -> flow space) and `StatsBase.reconstruct` (flow space ->
theta), so it is a DROP-IN for the plain `ZScoreTransform` at every existing spike call site.
`reconstruct` output is guaranteed to lie in `[lo, hi]`.

PORTED from `src/amortized/architecture.jl:88-104` (ruling Q4). The name carries the `P11`
prefix ON PURPOSE: it must never be mistaken for the `src` type in a stack trace, because the
two are structurally identical but belong to different nets -- the shipped one and the Phase-11
research one.
"""
struct P11BoundedThetaTransform{T<:Real,Z}
    lo::Vector{T}
    hi::Vector{T}
    zt::Z
end

_p11_logit_clamped(u::Real, logit_eps::Real) =
    (c = clamp(u, logit_eps, one(u) - logit_eps); log(c / (one(c) - c)))
_p11_logistic(y::Real) = inv(one(y) + exp(-y))

"""
    p11_theta_to_unbounded(theta, lo, hi; logit_eps = P11_THETA_LOGIT_EPS) -> Matrix

Per-row logit map of a `d x n` theta matrix from the prior box onto R. Rows are parameters.
Inputs are `clamp`ed into `[lo, hi]` first, so a caller that hands over synthetic (out-of-box)
theta degrades to the boundary rather than producing `NaN`.

Ported from `src/amortized/architecture.jl:116-129` (`theta_to_unbounded`).
"""
function p11_theta_to_unbounded(θ::AbstractMatrix, lo::AbstractVector, hi::AbstractVector;
                                logit_eps::Real = P11_THETA_LOGIT_EPS)
    size(θ, 1) == length(lo) == length(hi) || throw(DimensionMismatch(
        "p11_theta_to_unbounded: θ rows $(size(θ,1)) vs bounds $(length(lo))"))
    out = Matrix{Float64}(undef, size(θ))
    @inbounds for i in axes(θ, 1)
        w = hi[i] - lo[i]
        w > 0 || throw(ArgumentError("p11_theta_to_unbounded: degenerate bound on row $i"))
        for j in axes(θ, 2)
            out[i, j] = _p11_logit_clamped((clamp(θ[i, j], lo[i], hi[i]) - lo[i]) / w, logit_eps)
        end
    end
    return out
end

"""
    p11_theta_from_unbounded(Y, lo, hi) -> Matrix

Inverse of `p11_theta_to_unbounded`: `theta = lo + (hi - lo) * logistic(y)`, row-wise. The result
is GUARANTEED inside `[lo, hi]` (a final `clamp` absorbs float round-off at the extremes), which
is the whole point -- no posterior draw can leave the prior support.

Ported from `src/amortized/architecture.jl:138-149` (`theta_from_unbounded`).
"""
function p11_theta_from_unbounded(Y::AbstractMatrix, lo::AbstractVector, hi::AbstractVector)
    size(Y, 1) == length(lo) == length(hi) || throw(DimensionMismatch(
        "p11_theta_from_unbounded: rows $(size(Y,1)) vs bounds $(length(lo))"))
    out = similar(Y, float(eltype(Y)))
    @inbounds for i in axes(Y, 1)
        w = hi[i] - lo[i]
        for j in axes(Y, 2)
            out[i, j] = clamp(lo[i] + w * _p11_logistic(Y[i, j]), lo[i], hi[i])
        end
    end
    return out
end

# The two `StatsBase` methods, ported from `src/amortized/architecture.jl:151-161`.
function StatsBase.transform(t::P11BoundedThetaTransform, θ::AbstractMatrix)
    return StatsBase.transform(t.zt, p11_theta_to_unbounded(θ, t.lo, t.hi))
end
StatsBase.transform(t::P11BoundedThetaTransform, v::AbstractVector) =
    vec(StatsBase.transform(t, reshape(v, :, 1)))

function StatsBase.reconstruct(t::P11BoundedThetaTransform, Y::AbstractMatrix)
    return p11_theta_from_unbounded(StatsBase.reconstruct(t.zt, Y), t.lo, t.hi)
end
StatsBase.reconstruct(t::P11BoundedThetaTransform, v::AbstractVector) =
    vec(StatsBase.reconstruct(t, reshape(v, :, 1)))

"""
    p11_theta_prior_bounds() -> NTuple{8,Tuple{Float64,Float64}}

The Phase-11 research prior box, in `sample_prior` FIELD ORDER (`spike/simulator/prior.jl:93`):
`rho_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise,
chromatic_eps`.

**SINGLE SOURCE OF TRUTH.** Every bound is DERIVED from the prior objects `sample_prior` actually
draws from -- nothing here is a duplicated literal, so changing a prior in `prior.jl` moves this
box with it. This is the `src/amortized/simulator.jl:166-175` `theta_prior_bounds()` rule applied
to the spike's own priors; the two boxes DELIBERATELY DISAGREE on rows 5-6, because D-02 widened
`SHIFT_PRIOR` in the spike only (deviation (2)).
"""
p11_theta_prior_bounds() = (
    (first(GHAT_RHO_KNOTS), last(GHAT_RHO_KNOTS)),                       # rho_true = ghat(mu*)
    (minimum(SPILLOVER_PRIOR),        maximum(SPILLOVER_PRIOR)),
    (minimum(AUTOFLUORESCENCE_PRIOR), maximum(AUTOFLUORESCENCE_PRIOR)),
    (minimum(LABEL_EFFICIENCY_PRIOR), maximum(LABEL_EFFICIENCY_PRIOR)),
    (minimum(SHIFT_PRIOR),            maximum(SHIFT_PRIOR)),             # shift_dx (D-02 widened)
    (minimum(SHIFT_PRIOR),            maximum(SHIFT_PRIOR)),             # shift_dy (D-02 widened)
    (minimum(NOISE_PRIOR),            maximum(NOISE_PRIOR)),
    (minimum(CHROMATIC_PRIOR),        maximum(CHROMATIC_PRIOR)),         # chromatic_eps (D-09)
)

"""
    fit_p11_theta_transform(theta_train; bounds = p11_theta_prior_bounds())
        -> P11BoundedThetaTransform

Fit the leak-free FROZEN theta transform on the TRAIN theta columns ONLY (`dims = 2`): map out of
the truncated prior box onto R with the per-parameter logit bijection, then fit the
`ZScoreTransform` in THAT unconstrained space.

`bounds` is a keyword ONLY so tests can exercise a different box; real callers must let it
default so `prior.jl` stays the single source of truth. Rows whose logit-space spread is
degenerate (a synthetic all-identical theta row) get a unit scale instead of a 0 divisor.

Ported from `src/amortized/train_npe.jl:66-77` (ruling Q4).
"""
function fit_p11_theta_transform(θtr::AbstractMatrix; bounds = p11_theta_prior_bounds())
    lo = Float64[b[1] for b in bounds]
    hi = Float64[b[2] for b in bounds]
    size(θtr, 1) == length(lo) || throw(DimensionMismatch(
        "fit_p11_theta_transform: θ has $(size(θtr,1)) rows but the prior box has $(length(lo))"))
    U  = p11_theta_to_unbounded(θtr, lo, hi)
    zt = StatsBase.fit(StatsBase.ZScoreTransform, U; dims = 2)
    @inbounds for i in eachindex(zt.scale)      # degenerate-row guard (never a 0 divisor)
        (isfinite(zt.scale[i]) && zt.scale[i] > 0) || (zt.scale[i] = 1.0)
    end
    return P11BoundedThetaTransform(lo, hi, zt)
end

# =============================================================================================
# 2. The lambda conditioning encoder (D-03) -- spike-local, NEVER crosses into src/
# =============================================================================================

"""
    encode_lambda(lam) -> Float64

Map the registration-uncertainty level `lam` (the shift prior's half-width) onto `[0, 1]`:
`(lam - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN)`.

This is a FIXED, INVERTIBLE, **DATA-INDEPENDENT** affine map -- NOT a fitted transform. Nothing
about it is estimated from the training pool, so there is no leak to be leak-free about, and
there is no frozen statistic to persist alongside the net beyond the two constants
`LAMBDA_MIN`/`LAMBDA_MAX`, which are themselves Tier-1 pre-registered
(`spike/validation/p11_consts.jl`). Values outside `[LAMBDA_MIN, LAMBDA_MAX]` are NOT clamped:
reading the net outside its trained lambda range is an extrapolation the caller must own, and
silently clamping would hide it.
"""
encode_lambda(lam) = (lam - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN)

"""
    decode_lambda(u) -> Float64

Inverse of `encode_lambda`: `LAMBDA_MIN + u * (LAMBDA_MAX - LAMBDA_MIN)`. Kept for symmetry and
for the figures, which label ladder rungs in PHYSICAL px on the axis while the net sees the
encoded value.
"""
decode_lambda(u) = LAMBDA_MIN + u * (LAMBDA_MAX - LAMBDA_MIN)

"""
    augment_input(Z128, lam) -> Matrix{Float32}

Append the encoded lambda as a 129th input row to an already-standardized 128-row summary matrix
(D-03, deviation (3)). Returns a `129 x size(Z128, 2)` matrix whose last row is constant across
columns.

**ORDER INVARIANT -- THIS IS A CORRECTNESS CONTRACT, NOT A STYLE PREFERENCE (Pitfall 2).**
`Z128` MUST ALREADY BE THE OUTPUT OF THE FROZEN PATH
`standardize_summary(encode_d01(M), zt, :min)`. Standardize the 128 rows FIRST, then append.

Appending BEFORE standardizing fails one of two ways, and the second is the dangerous one:
  * loudly -- `Loader._row_partition` (`spike/data/loader.jl:126`) requires EXACTLY 128 rows for
    `:min` and errors on 129; the `src` twin errors on an odd count. Both errors are GOOD and
    must never be "fixed" to accept 129 rows; and
  * silently -- if someone does "fix" the partition, the conditioning variable gets z-scored
    against the TRAINING POOL's lambda distribution, coupling a user-declared input to the
    training data and destroying the read-at-a-chosen-lambda sweep that is the whole point of
    D-03.
The `@assert` below enforces the invariant AT RUNTIME so it is checked, not merely documented;
it also catches a double-append (129 in) and a raw-64-row summary.
"""
function augment_input(Z128::AbstractMatrix, lam)
    @assert size(Z128, 1) == 128 "augment_input: expected the frozen 128-row standardized " *
        "summary, got $(size(Z128, 1)) rows — standardize FIRST, then append λ (Pitfall 2)"
    return vcat(Float32.(Z128), fill(Float32(encode_lambda(lam)), 1, size(Z128, 2)))
end
augment_input(v::AbstractVector, lam) = augment_input(reshape(v, :, 1), lam)

# =============================================================================================
# 3. The research estimator (D-03) -- d_in = 129, D = 8
# =============================================================================================

"""
    build_p11_estimator(; d_in = 129, D = 8, dstar = 64, depth = 3, width = 256,
                        num_coupling_layers = 10, flow_depth = 2, flow_width = 128)
        -> PosteriorEstimator

The Phase-11 research estimator: a thin wrapper over `build_estimator` at `d_in = 129`
(128 frozen summary rows + the lambda row, deviation (3)) and `D = 8` (the 7 D-04 theta fields
plus `chromatic_eps`, deviation (1)).

**NO ARCHITECTURE CODE CHANGES.** `build_estimator` is input-width-agnostic and its `D` is
already POSITIONAL (`spike/npe/architecture.jl:87`), so `build_estimator(129, 8; ...)` needs no
signature change, and the `dstar >= D` guard (`architecture.jl:95`) holds comfortably at
64 >= 8.

**CAPACITY IS HELD FIXED ON PURPOSE.** Every knob defaults to the unchanged Phase-5-calibrated
value (`NPE_DSTAR`/`NPE_DEPTH`/`NPE_WIDTH`/`NPE_COUPLING`/`NPE_FLOW_DEPTH`/`NPE_FLOW_WIDTH`), so
the comparison against the shipped net is CAPACITY-CONTROLLED and the only things that moved are
the four declared deviations. Spike 012 additionally found this architecture near its useful
capacity ceiling, so raising capacity is contraindicated as well as off-protocol.

`NPE_D` (7) and `build_estimator`'s default `D` are NOT touched by this wrapper -- the shipped
read surface keeps its 7 marginals (D-16).
"""
function build_p11_estimator(; d_in::Integer = 129, D::Integer = 8,
                             dstar::Integer = NPE_DSTAR, depth::Integer = NPE_DEPTH,
                             width::Integer = NPE_WIDTH,
                             num_coupling_layers::Integer = NPE_COUPLING,
                             flow_depth::Integer = NPE_FLOW_DEPTH,
                             flow_width::Integer = NPE_FLOW_WIDTH)
    return build_estimator(d_in, D; dstar = dstar, depth = depth, width = width,
                           num_coupling_layers = num_coupling_layers,
                           flow_depth = flow_depth, flow_width = flow_width)
end
