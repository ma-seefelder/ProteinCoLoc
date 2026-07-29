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

# spike/simulator/p12_prior.jl --- D-05: the MARGINAL-PRESERVING copula that turns the lattice
# arithmetic of p12_lattice.jl into a prior over rho FIELDS.
#
# THE CLAIM, AND THE REASON IT IS SAFE. The lattice prior contributes ONLY SPATIAL DEPENDENCE.
# Each region's marginal is pushed through Phi, then through the MU_PRIOR quantile, then through
# the frozen ghat -- so every one of the G^2 regions carries EXACTLY the prior `prior.jl:40` puts
# on the scalar knob, and the CLAUDE.md constraint that pi(theta) match the Turing mu/nu/sigma/tau
# ranges is untouched by anything in this file. MU_PRIOR and all seven nuisance priors are REACHED
# through the guarded include of prior.jl and are never restated here: a re-declaration would fork
# pi(theta) from the ADVI baseline silently, while looking like a convenience.
#
# THE SAMPLED OBJECT IS THE GAUSSIAN FIELD; rho IS ITS MONOTONE IMAGE (R-2 / Pitfall 7). That is
# why `z` is part of the return value and not an implementation detail. In rho space the ELEMENTWISE
# ghat clamp puts a measured 6.79 % atom mass on EVERY region (12-RESEARCH.md Pattern 2, measured
# 2026-07-27, against a global 1-D reference of 6.78 %); in Gaussian space there are no atoms at
# all. The theta rows downstream (12-09) are built from `z`, and 12-18's clean SBC arm has no other
# input, so `z` must not be dropped because "the deliverable is a Delta-rho map" -- the deliverable
# is DERIVED from `z`.
#
# THE SAME TRAP AS prior.jl:30-32, NOW PER REGION. Do NOT assert `rho_field[r] == mu_field[r]`:
# the nonneg softplus mixing + PSF + noise attenuate rho into the induced mu nonlinearly, and ghat
# inverts exactly that measured map. What holds per region is the MARGINAL claim -- `mu_field[r]`
# is distributed as MU_PRIOR -- not a pointwise identity.
#
# WHAT THIS FILE DOES NOT KNOW ABOUT. theta vectors, DCT rows, networks. It draws fields and
# nuisances. The 72-row theta assembly lives in `spike/data/p12_generate.jl` (12-09), the only file
# that includes both this and `p12_architecture.jl`, so the row layout is defined exactly once at
# `P12_THETA_ROWS`.
#
# ZERO NEW DEPENDENCIES. `Distributions` and `Random` are already in spike/Project.toml;
# `Statistics` is a stdlib reached through Julia's default `@stdlib` LOAD_PATH entry exactly as
# `run_p11_recovery.jl:50-53` reaches it, and is deliberately NOT added to spike/Project.toml
# (`spike/test/runtests.jl:123-126` asserts the exact 16-name dependency set).
#
# Run (self-check):
#     julia --project=spike -e 'include("spike/simulator/p12_prior.jl")'

using Distributions          # Normal, cdf, quantile (the copula); MU_PRIOR's own type
using Random                 # AbstractRNG
using Statistics             # mean (the derived read-time scalar); stdlib, see header

# ORDER MATTERS, AND prior.jl IS FIRST. It binds MU_PRIOR, SIM02_W1_TOL and every nuisance prior
# this file reuses VERBATIM -- Phase 12 introduces NO new nuisance prior and must not restate one.
# It also transitively includes the frozen `ghat.jl` (prior.jl:37, an UNGUARDED include), so the
# guard below is what keeps ghat's knot vectors from being re-bound on a second load.
# p12_lattice.jl second: it supplies `lattice_sigma` / `field_sampler` and transitively guard-
# includes the Tier-1 pre-registration, which is where P12_G and P12_R1_PRIOR come from.
isdefined(@__MODULE__, :MU_PRIOR)      || include(joinpath(@__DIR__, "prior.jl"))
isdefined(@__MODULE__, :lattice_sigma) || include(joinpath(@__DIR__, "p12_lattice.jl"))

# =====================================================================================
# 0. THE r1 PRIOR IS READ, NEVER RE-DECLARED
# =====================================================================================
# `P12_R1_PRIOR` arrives through the guarded includes above (p12_lattice.jl -> p12_consts.jl §4).
# THERE IS DELIBERATELY NO `const P12_R1_PRIOR = ...` IN THIS FILE. 12-01 binds it in Tier 1
# expressly so that no plan constructs its own and two arms quietly compare different priors --
# the R-8 trap again. A local `const` carrying the SAME value would be benign in Julia and is
# therefore WORSE than a wrong one: it silently defeats the guard, and the day Tier 1 is read
# differently the two copies diverge with nothing to catch it.
#
# WHY THE PRIOR IS ON THE INDUCED LAG-1 CORRELATION AND NOT ON THE CAR alpha (R-8). Measured lag-1
# runs 0.136 at alpha = 0.5 to 0.947 at alpha = 0.999, so a uniform prior on alpha would put
# roughly 90 % of its mass on "no spatial structure". That walks D-08 straight into the S-1
# prior-echo trap -- the correlation length would be unidentifiable BY PARAMETRIZATION rather than
# BY PHYSICS -- while appearing to test a prior SHAPE. Both arms are addressed by induced r1 so
# the CAR-vs-GP comparison is a question about the kernels, not about their coordinates.

# =====================================================================================
# 1. ASVS V5: a field must validate itself
# =====================================================================================

"""
    validate_rho_field(F) -> F

Throw `ArgumentError` naming the offending value when `F` is not square, carries a non-finite
entry, or carries an entry outside `[-1, 1]`.

THE SCALAR GUARD IN `simulate_pair` CANNOT SEE INSIDE A FIELD. `forward.jl:136-158` validates
`theta.rho_true` as a single number; once rho is a `G x G` field that check is structurally blind
to 63 of the 64 values it is supposed to cover. So the field validates itself HERE, before it can
ever reach the forward model, mirroring that block's discipline: reject rather than clamp, and
name the failing value in the message.
"""
function validate_rho_field(F::AbstractMatrix)
    size(F, 1) == size(F, 2) || throw(ArgumentError(
        "validate_rho_field: expected a square G×G rho field, got $(size(F,1))×$(size(F,2))"))
    @inbounds for k in eachindex(F)
        v = F[k]
        isfinite(v) || throw(ArgumentError(
            "validate_rho_field: non-finite rho at flat index $k, got $v"))
        (-1.0 <= v <= 1.0) || throw(ArgumentError(
            "validate_rho_field: rho must be in [-1, 1], got $v at flat index $k"))
    end
    return F
end

# =====================================================================================
# 2. The copula (D-05)
# =====================================================================================

"""
    rho_field(rng, sampler; G = P12_G) -> (rho = G×G, mu = G×G, z = G×G)

Draw one spatial rho field through the D-05 copula.

    z = sampler(rng)                  # marginally N(0,1) PER CELL -- field_sampler already
                                      # divides by sqrt.(diag(Sigma)), so do NOT re-standardize
    u = cdf.(Normal(), z)             # Uniform(0,1) per region
    mu = quantile.(Ref(MU_PRIOR), u)  # MU_PRIOR per region, EXACTLY
    rho = ghat.(mu)                   # the frozen SIM-02 inverse, ELEMENTWISE

`ghat` IS APPLIED ELEMENTWISE AND IS OTHERWISE UNCHANGED -- no knot is touched, no endpoint is
moved. The elementwise application is also what makes the clamp's atoms PER REGION: a measured
6.79 % of entries land on +/-0.99, per region, which is the SAME RATE as the current global 1-D
prior (6.78 %) but 64x more frequent per draw. That is why `z` -- which is atom-free -- is the row
the head is parametrized on (R-2), and why the randomized-rank budget is planned rather than
discovered.

The returned `rho` is passed through [`validate_rho_field`](@ref) before it is returned.
"""
function rho_field(rng::AbstractRNG, sampler; G::Integer = P12_G)
    z = sampler(rng)
    length(z) == G * G || throw(ArgumentError(
        "rho_field: sampler returned $(length(z)) values; expected G² = $(G*G) for G = $G"))
    u = cdf.(Normal(), z)
    μ = quantile.(Ref(MU_PRIOR), u)
    ρ = ghat.(μ)
    F = reshape(ρ, G, G)
    validate_rho_field(F)
    return (rho = F, mu = reshape(μ, G, G), z = reshape(collect(float.(z)), G, G))
end

# =====================================================================================
# 3. The field-aware prior draw -- AND ITS ORDER
# =====================================================================================

"""
    sample_p12_prior(rng; G = P12_G, arm = :car, r1 = nothing) -> NamedTuple

Draw one Phase-12 prior sample: a spatial rho field plus the seven non-rho nuisances.

THE GENERATIVE STORY, IN EXACTLY THIS ORDER, AND NO OTHER:

    r1_i    ~ P12_R1_PRIOR                       <- THE CORRELATION LENGTH, DRAWN **FIRST**
    Sigma_i  = lattice_sigma(arm, r1_i)          <- the covariance AT that induced lag-1
    field_i ~ field_sampler(Sigma_i) -> copula   <- THE FIELD, DRAWN **CONDITIONAL ON IT**
    spillover_i, autofluorescence_i, label_efficiency_i, shift_dx_i, shift_dy_i,
    noise_i, chromatic_eps_i ~ their EXISTING priors, in their EXISTING relative order

THIS ORDER IS THE ENTIRE POINT OF THE FUNCTION, and it is the `p11_generate.jl:20-33` discipline
transposed onto a field. Getting it backwards -- drawing a field and then reporting a correlation
length for it -- manufactures a null that looks EXACTLY like "the correlation length is below the
8x8 summary's resolution" while actually being an implementation bug. Phase 11's header records
what that costs; D-08 is a question about physics only if the physics was drawn first.

`arm` IS CARRIED IN THE RETURN VALUE so a persisted sample records which prior arm produced it.
An unrecorded arm makes the 12-15 CAR-vs-GP mini-spike unreadable after the fact: two pools that
differ only in a keyword nobody wrote down are two pools that cannot be told apart.

`r1` may be PINNED (a number) instead of drawn, which is what the r1 ladder and the D-10 ablation
use. Pinning consumes less of `rng` than drawing does, by construction -- see the note in
`test_p12_prior.jl` testset 4, which asserts exactly that as the behavioural proof of the order.
"""
function sample_p12_prior(rng::AbstractRNG; G::Integer = P12_G,
                          arm::Symbol = :car, r1 = nothing)
    # (1) the correlation length FIRST
    r1v = r1 === nothing ? rand(rng, P12_R1_PRIOR) : Float64(r1)
    # (2) the field CONDITIONAL on it
    Σ   = lattice_sigma(arm, r1v; G = G)
    f   = field_sampler(Σ)
    fld = rho_field(rng, f; G = G)
    # (3) the seven non-rho nuisances, from prior.jl's OWN priors, in prior.jl's OWN order
    return (
        rho_field        = fld.rho,
        mu_field         = fld.mu,
        z_field          = fld.z,
        r1               = r1v,
        arm              = arm,
        spillover        = rand(rng, SPILLOVER_PRIOR),
        autofluorescence = rand(rng, AUTOFLUORESCENCE_PRIOR),
        label_efficiency = rand(rng, LABEL_EFFICIENCY_PRIOR),
        shift_dx         = rand(rng, SHIFT_PRIOR),
        shift_dy         = rand(rng, SHIFT_PRIOR),
        noise            = rand(rng, NOISE_PRIOR),
        chromatic_eps    = rand(rng, CHROMATIC_PRIOR),   # D-09: appended LAST (positional θ)
    )
end

"""
    theta_scalar_view(draw) -> NamedTuple

The LEGACY 8-field theta NamedTuple `simulate_pair` already accepts, built from a Phase-12 draw.

`rho_true` HERE IS A **DERIVED** READ-TIME SCALAR (R-1), NOT THE OBJECT THE HEAD ESTIMATES. It is
`ghat(quantile(MU_PRIOR, cdf(Normal(), mean(z_field))))` -- the copula applied to the MEAN of the
Gaussian field, i.e. the one number a scalar-rho simulator can consume. The head is parametrized
on the field's DCT coefficients; collapsing 64 regions to one number is a lossy projection, and
D-07 forbids SCORING such a projection. This view therefore exists for exactly two purposes:
comparability against the scalar lane, and validation plumbing (feeding the unmodified
`simulate_pair` so the legacy path can be shown to still accept a Phase-12 draw). Do not rank it,
do not report coverage on it, and do not treat it as "the" rho of the draw.
"""
function theta_scalar_view(draw)
    ρ_scalar = ghat(quantile(MU_PRIOR, cdf(Normal(), mean(draw.z_field))))
    return (
        ρ_true           = ρ_scalar,
        spillover        = draw.spillover,
        autofluorescence = draw.autofluorescence,
        label_efficiency = draw.label_efficiency,
        shift_dx         = draw.shift_dx,
        shift_dy         = draw.shift_dy,
        noise            = draw.noise,
        chromatic_eps    = draw.chromatic_eps,
    )
end

# =====================================================================================
# 4. G×G field -> pixel grid, in two matmuls (D-06)
# =====================================================================================
# THE COST IS THE REASON THIS IS SEPARABLE AND NOT A PER-PIXEL INTERPOLANT. Measured
# (12-RESEARCH.md Pattern 3, 2026-07-27) at 1376x1028: a scalar interpolant evaluated per output
# pixel costs 0.3622 s against 0.0252 s for `A * F * B'` -- 14x. Over a 50 000-pair pool that is
# ~5 HOURS added to datagen against ~8 MINUTES, i.e. the difference between clearing and blowing
# `P12_DATAGEN_WALLCLOCK_CEILING_MIN`. A plan that says "bilinear upsample" without naming the
# implementation silently produces the 5-hour version.

"""
    p12_bilinear_op(G, m; offset = false) -> Matrix{Float64}

The `m x G` bilinear interpolation operator taking a length-`G` lattice axis to `m` pixels: row
`r` carries the two weights `1-t` and `t` of the cell pair the r-th sample position falls between.
Every row sums to 1, so a constant field renders to that constant.

`offset = true` **IS** THE D-06 HALF-CELL GUARD ARM. It samples `range(1 - 0.5, G + 0.5;
length = m)` instead of `range(1, G; length = m)`, offsetting the field lattice by half a cell
against the summary's patch boundaries. ONE KEYWORD, NOT A SECOND CODE PATH -- so the guard
cannot drift from the thing it guards, which is the whole failure mode a separately-written
"offset version" would reintroduce.

At `offset = true` the two half-cell margins fall OUTSIDE `[1, G]`, so the end rows carry weights
`(1.5, -0.5)`: linear EXTRAPOLATION, not interpolation. That is deliberate and is what "the
lattice is offset" means; the row still sums to 1, so the level of the field is untouched and
only the sampling POSITIONS move.
"""
function p12_bilinear_op(G::Integer, m::Integer; offset::Bool = false)
    G >= 2 || throw(ArgumentError("p12_bilinear_op: G must be ≥ 2, got $G"))
    m >= 1 || throw(ArgumentError("p12_bilinear_op: m must be ≥ 1, got $m"))
    A  = zeros(Float64, m, G)
    xs = offset ? range(1 - 0.5, G + 0.5; length = m) : range(1, G; length = m)
    @inbounds for (r, x) in enumerate(xs)
        i = clamp(floor(Int, x), 1, G - 1)
        t = x - i
        A[r, i]     = 1 - t
        A[r, i + 1] = t
    end
    return A
end

"""
    p12_upsample_ops(G, imsize; offset = false) -> (A, B)

The two separable operators for a `G x G` field rendered to `imsize` pixels, so a datagen loop can
HOIST their construction out of the inner loop: they depend only on `(G, imsize, offset)`, never
on the draw.
"""
p12_upsample_ops(G::Integer, imsize::Tuple{Int,Int}; offset::Bool = false) =
    (p12_bilinear_op(G, imsize[1]; offset = offset),
     p12_bilinear_op(G, imsize[2]; offset = offset))

"""
    p12_upsample(F, imsize; offset = false) -> Matrix{Float64}

Render the lattice field `F` to `imsize` pixels as `A * F * B'` -- two matmuls, no per-pixel loop
and no per-pixel interpolant object. See the cost note above the operator: 0.0252 s vs 0.3622 s at
1376x1028.

`B'` (the TRANSPOSE) is not cosmetic: `B` is `n x size(F,2)`, so `F * B'` is what contracts the
field's second axis onto the output columns. `A * F * B` is a different (and usually
dimension-invalid) expression, which `test_p12_prior.jl` testset 6 exists to catch.
"""
function p12_upsample(F::AbstractMatrix, imsize::Tuple{Int,Int}; offset::Bool = false)
    A = p12_bilinear_op(size(F, 1), imsize[1]; offset = offset)
    B = p12_bilinear_op(size(F, 2), imsize[2]; offset = offset)
    return A * F * B'
end
