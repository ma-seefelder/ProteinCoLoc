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

# spike/simulator/forward.jl --- SIM-01: the 2-channel 2D forward physics
# simulator `simulate_pair(rng, θ; imsize)`.
#
# Runs the seven D-06 stages IN ORDER on an explicitly-threaded rng (D-14):
#   (1) correlated bivariate intensity fields  -- the SPIKE-VALIDATED
#       shared-latent smooth Gaussian-field generator (D-15), NOT puncta:
#       a smooth latent field L plus two independent smoothed fields ε₁,ε₂
#       are mixed and pushed through softplus so ρ_true monotonically (and
#       signed) drives the induced per-patch correlation through the REAL
#       summary -- the de-risk experiment reached induced-μ ∈ [−0.76,+0.76]
#       (Spearman 1.0) only when the field-smoothing scale STRUCT_σ ≫ the PSF
#       width σ_psf (a puncta/spot mixture saturated ~0.21 and could not go
#       negative -- rejected, D-15).
#   (2) Bernoulli thinning      (label_efficiency = keep-probability)
#   (3) fixed Gaussian PSF      imfilter(·, Kernel.gaussian(σ_psf))  [D-05 nuisance, NOT in θ]
#   (4) 2×2 directional spillover  ch1 += spillover·ch2  (directional, NOT symmetric)
#   (5) autofluorescence offset + BG_FLOOR  (background stays small-positive, never hard 0.0)
#   (6) misregistration of ch2 only -- ONE composed affine (D-10): a radial
#       chromatic magnification difference about the image centre (chromatic_eps)
#       COMPOSED with the sub-pixel registration shift (dy, dx), applied in a
#       SINGLE resampling pass with fillvalue = BG_FLOOR. One pass is a
#       correctness requirement, not a style choice: every interpolation pass
#       smooths and smoothing decorrelates, so a second pass would widen the
#       posterior for a reason unrelated to misalignment.
#   (7) Poisson(shot) + Gaussian(read) noise   (scaled by θ.noise)  [D-07]
#
# Output is a 2-element Vector{Matrix{Float64}} (all finite, ≥ 0, mostly-non-zero)
# that flows UNCHANGED through spike/contract.jl's build_mci/summary (SIM-03).
#
# θ is the 8-field NamedTuple (the 7 D-04 fields plus chromatic_eps, D-09):
#   (ρ_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy,
#    noise, chromatic_eps)
# chromatic_eps is APPENDED LAST (downstream code indexes θ positionally) and is
# read DEFENSIVELY, so a legacy 7-field θ remains a valid input meaning
# "no chromatic aberration" and reproduces the pre-D-09 output bit for bit.

using Distributions
using ImageFiltering                 # imfilter, Kernel.gaussian (PSF + field smoothing)
using ImageTransformations           # warp (sub-pixel shift)
using CoordinateTransformations      # Translation, LinearMap, recenter (the composed affine)
using Interpolations                 # BSpline, Linear (warp interpolation scheme)
using Statistics                     # mean, std (field standardization)
using Random                         # AbstractRNG, Xoshiro

# --- Fixed nuisance constants (Claude's discretion within RESEARCH ranges) ------
# σ_psf: diffraction-limited fixed PSF (D-05), ≈1.0-1.5 px (RESEARCH Pattern 2).
const σ_psf      = 1.3
# STRUCT_σ: field-smoothing scale. MUST be ≫ σ_psf so the induced spatial
# correlation survives PSF+noise (D-15: the de-risk used ≈6 px vs 1.5 px PSF).
const STRUCT_σ   = 6.0
# Poisson shot gain + Gaussian read-noise sd (plausible microscopy ranges, D-07).
const GAIN       = 50.0
const READ_NOISE = 0.01
# Small-positive background floor: keeps every pixel > 0 so _exclude_zero does
# NOT thin patches below the 15-px floor (SIM-03 background-not-zero trap).
const BG_FLOOR   = 0.02

# Numerically-stable softplus = log(1 + exp(x)); the nonneg intensity transform.
_softplus(x::Real) = log1p(exp(-abs(x))) + max(x, zero(x))

"""
    _smooth_field(rng, imsize) -> Matrix{Float64}

A single smooth latent Gaussian field: white noise low-pass-filtered at the
STRUCT_σ structural scale, then **standardized to zero-mean / unit-variance**.
The building block of the shared-latent generator.

Standardization is load-bearing (D-15): Gaussian smoothing at STRUCT_σ collapses
the raw field variance to ≈1/(4π·STRUCT_σ²), which would bury the correlation
signal under Poisson/read noise and flatten induced-μ to ≈0. Renormalizing to
unit variance restores full amplitude AND makes the linear mix `a·L + b·ε`
carry correlation exactly `sign(ρ)·a² = ρ` into the softplus transform.
"""
function _smooth_field(rng::AbstractRNG, imsize::Tuple{Int,Int})
    f = imfilter(randn(rng, imsize...), Kernel.gaussian(STRUCT_σ))
    # IN-03: floor the denominator at 1e-10 instead of adding eps()≈2.2e-16; a near-zero
    # std divided by eps() would explode to ~1e15, whereas max(·,1e-10) caps the blow-up.
    return (f .- mean(f)) ./ max(std(f), 1e-10)
end

"""
    simulate_pair(rng::AbstractRNG, θ; imsize=(256,256)) -> Vector{Matrix{Float64}}

Generate a 2-channel microscopy image pair from the 8-field θ NamedTuple (the 7
D-04 fields plus `chromatic_eps`, D-09) by running the seven D-06 stages in order
on the explicitly-threaded `rng` (D-14). A legacy 7-field θ is still accepted and
is treated as `chromatic_eps = 0`, which reproduces the pre-D-09 output exactly.

Stage 1 is the spike-validated shared-latent correlated smooth Gaussian-field
generator (D-15): `c1 = softplus(√|ρ|·L + √(1−|ρ|)·ε₁)`,
`c2 = softplus(sign(ρ)·√|ρ|·L + √(1−|ρ|)·ε₂)` -- the `sign(ρ)` on channel-2's
shared component makes anti-correlation reachable at negative ρ_true. Higher
ρ_true ⇒ monotonically higher induced cross-channel patch correlation.

Returns `[ch1, ch2] :: Vector{Matrix{Float64}}`, each `imsize`, all entries finite
and ≥ 0 with a small-positive background (no hard zeros), ready for `build_mci`.

Throws `ArgumentError` for `ρ_true ∉ [-1,1]`, any out-of-range / non-finite θ
field, or an `imsize` with a dimension < 64 (CR-02: the 8×8 patch grid needs
≥ 8 px/patch-side so patches clear the ≥15-px floor).
"""
function simulate_pair(rng::AbstractRNG, θ; imsize::Tuple{Int,Int} = (256, 256))
    # --- θ / imsize validation at entry (T-02-IV: untrusted parameter vector) ---
    # CR-02: the 8×8 patch grid divides each axis into 8 patches; src/colocalization.jl
    # drops a patch to `missing` when it has ≤ 15 surviving pixels. At imsize < 64 a
    # patch side is < 8 px (< 64 px/patch) and (8,8) gives 1-px patches → an all-missing
    # summary that crashed induced_mu. Require ≥ 64 so each patch is ≥ 8×8 = 64 px.
    (imsize[1] ≥ 64 && imsize[2] ≥ 64) ||
        throw(ArgumentError("imsize dims must be ≥ 64 for the 8×8 patch grid to " *
                            "produce non-missing patches (>15 px each), got $imsize"))
    # D-09 backward compatibility (SC1e): chromatic_eps is read through ONE defensive local
    # binding, computed here so both the guard below and stage 6 share it. A legacy 7-field θ
    # (e.g. the θ reconstructed by src/amortized/ood.jl `_theta_tuple` from a 7-row posterior
    # mean) is therefore a valid input meaning "no chromatic aberration".
    chromatic_eps_val = hasproperty(θ, :chromatic_eps) ? θ.chromatic_eps : 0.0
    all(isfinite, (θ.ρ_true, θ.spillover, θ.autofluorescence, θ.label_efficiency,
                   θ.shift_dx, θ.shift_dy, θ.noise, chromatic_eps_val)) ||
        throw(ArgumentError("all θ fields must be finite, got $θ"))
    (-1.0 ≤ θ.ρ_true ≤ 1.0) ||
        throw(ArgumentError("ρ_true must be in [-1, 1], got $(θ.ρ_true)"))
    # WR-07: reject physically out-of-range nuisances rather than silently clamping
    # (label_efficiency) or corrupting downstream stages (negative spillover would
    # SUBTRACT signal in stage 4; negative autofluorescence/noise are non-physical).
    (0.0 ≤ θ.spillover ≤ 1.0) ||
        throw(ArgumentError("spillover must be in [0, 1], got $(θ.spillover)"))
    (0.0 ≤ θ.autofluorescence) ||
        throw(ArgumentError("autofluorescence must be ≥ 0, got $(θ.autofluorescence)"))
    (0.0 ≤ θ.label_efficiency ≤ 1.0) ||
        throw(ArgumentError("label_efficiency must be in [0, 1], got $(θ.label_efficiency)"))
    (0.0 ≤ θ.noise) ||
        throw(ArgumentError("noise must be ≥ 0, got $(θ.noise)"))
    # WR-07 (D-09): the stage-6 backward scale is s = 1/(1 + chromatic_eps), so 1 + chromatic_eps
    # must stay STRICTLY POSITIVE. At chromatic_eps = -1 the scale is a division by zero; below
    # -1 it goes negative and the affine map MIRRORS channel 2 instead of magnifying it. This is
    # a correctness bound, not a style choice (T-11-07).
    (-1.0 < chromatic_eps_val) ||
        throw(ArgumentError("chromatic_eps must be > -1 (so 1 + chromatic_eps stays positive), " *
                            "got $(chromatic_eps_val)"))

    ρ = θ.ρ_true

    # --- (1) shared-latent correlated smooth Gaussian fields (D-15) -------------
    # Smooth latent L shared across channels; ε₁,ε₂ independent per-channel fields.
    L  = _smooth_field(rng, imsize)
    ε1 = _smooth_field(rng, imsize)
    ε2 = _smooth_field(rng, imsize)
    a  = sqrt(abs(ρ))           # shared-component weight
    b  = sqrt(1.0 - abs(ρ))     # private-component weight  (a² + b² = 1)
    # sign(ρ) flips channel-2's shared component so negative ρ ⇒ anti-correlation.
    ch1 = _softplus.(a .* L .+ b .* ε1)
    ch2 = _softplus.(sign(ρ) .* a .* L .+ b .* ε2)

    # --- (2) Bernoulli thinning (label_efficiency = keep-probability) -----------
    p     = clamp(θ.label_efficiency, 0.0, 1.0)
    keep1 = rand(rng, Bernoulli(p), imsize...)
    keep2 = rand(rng, Bernoulli(p), imsize...)
    ch1 .= ch1 .* keep1
    ch2 .= ch2 .* keep2

    # --- (3) fixed Gaussian PSF (D-05 nuisance; width is NOT read from θ) --------
    psf = Kernel.gaussian(σ_psf)
    ch1 = imfilter(ch1, psf)            # PSF also re-smooths the thinning zeros away
    ch2 = imfilter(ch2, psf)

    # --- (4) 2×2 directional spillover (bleed-through, directional ≠ symmetric) --
    ch1 = ch1 .+ θ.spillover .* ch2     # only ch1 receives ch2's bleed-through

    # --- (5) autofluorescence offset + small-positive background floor ----------
    ch1 = ch1 .+ θ.autofluorescence .+ BG_FLOOR
    ch2 = ch2 .+ θ.autofluorescence .+ BG_FLOOR

    # --- (6) ONE composed affine on channel 2 only: radial chromatic scale ∘ registration
    #         shift (D-10; fillvalue > 0) ------------------------------------------------
    # WR-06 (RE-DERIVED for the composed map, D-10). `warp` is BACKWARD-mode:
    #     out[I] = ch2[A(I)]        (ImageTransformations warp.jl:4 and :167-176)
    # so `A` maps a DESTINATION index to the SOURCE index it samples, and the image CONTENT
    # therefore undergoes A⁻¹. Two consequences:
    #   • The Translation slot order is UNCHANGED by the composition: slot 1 is the FIRST array
    #     axis (rows = vertical = dy), slot 2 the SECOND (columns = horizontal = dx). Pass
    #     (dy, dx) exactly as before, so the named fields keep their conventional physical axes.
    #   • To MAGNIFY channel 2's content by (1 + chromatic_eps) about the image centre, the
    #     BACKWARD map must SHRINK coordinates by that factor — hence s = 1/(1 + chromatic_eps).
    #     The guard above keeps 1 + chromatic_eps strictly positive, so s is finite and positive.
    # `recenter(t, c) == Translation(c) ∘ t ∘ Translation(-c)` (CoordinateTransformations
    # core.jl:103-105), which is what applies the scale about `c`, the geometric centre of
    # axes(ch2) (128.5, 128.5 at 256²) rather than about the origin. `c` is a plain Tuple —
    # recenter's Tuple method converts it internally, so no StaticArrays import is needed.
    # The whole composition COLLAPSES to a single AffineMap (compose methods, affine.jl:139-165),
    # so exactly ONE interpolation pass runs — the D-10 correctness requirement, because every
    # resampling pass smooths and smoothing decorrelates, which would contaminate the very
    # posterior width SC2 measures. At chromatic_eps = 0 the linear part is the exact identity
    # and the result is bit-identical to the previous translation-only stage (regression-tested
    # against spike/test/fixtures/p11_stage6_golden.jld2).
    s = 1.0 / (1.0 + chromatic_eps_val)
    c = map(ax -> (first(ax) + last(ax)) / 2, axes(ch2))
    A = Translation(θ.shift_dy, θ.shift_dx) ∘ recenter(LinearMap([s 0.0; 0.0 s]), c)
    shifted = warp(ch2, A, axes(ch2); method = BSpline(Linear()), fillvalue = BG_FLOOR)
    ch2 = Matrix{Float64}(collect(shifted))

    # --- (7) Poisson(shot) + Gaussian(read) noise (scaled by θ.noise) -----------
    read_sd = max(READ_NOISE * θ.noise, 0.0)
    ch1 = Float64.(rand.(Ref(rng), Poisson.(GAIN .* max.(ch1, 0.0)))) ./ GAIN .+
          rand(rng, Normal(0.0, read_sd), imsize)
    ch2 = Float64.(rand.(Ref(rng), Poisson.(GAIN .* max.(ch2, 0.0)))) ./ GAIN .+
          rand(rng, Normal(0.0, read_sd), imsize)

    # Clamp to non-negative; the softplus baseline keeps the background well above
    # zero, so this only trims rare read-noise undershoots (no hard-zero regions).
    return [Matrix{Float64}(max.(ch1, 0.0)), Matrix{Float64}(max.(ch2, 0.0))]
end
