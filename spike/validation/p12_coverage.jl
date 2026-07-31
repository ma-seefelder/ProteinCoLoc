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

# spike/validation/p12_coverage.jl --- D-09's leave-region-out predictive construction, the
# Fisher-z observation-noise model, the two-sided scoring rule, and `p12_coloc_map` (SPAT-05).
#
# THIS FILE IS A LIBRARY. It defines functions and runs nothing.
#
# =============================================================================================
# THE D-09 CONSTRUCTION, THE SEVEN STEPS, VERBATIM
# =============================================================================================
#
#   1. Encode the observed image ONCE:  Z_full = encode_d01(patch_summary(mci)).  RAW 128 rows.
#   2. For each held-out region r, mask r on the RAW rows -- value row and mask row both to 0 --
#      via `mask_regions`. That is byte-identically the encoding `encode_d01` already produces
#      for a patch below the >=15-surviving-pixel floor, so a held-out region is IN-VOCABULARY
#      rather than a new sentinel the net was never taught (Pitfall 4).
#   3. Standardize the masked rows with the FROZEN `bundle.zt` (`standardize_p12`), then
#      `reshape_summary` to the (G, G, 2, K) spatial view. MASK ON RAW, STANDARDIZE SECOND --
#      the ordering `train_p12_npe.jl:30-42` establishes and 12-14's testset 3 pins.
#   4. `sampleposterior` -> UN-STANDARDIZE with `bundle.theta_zt` -> reconstruct region r's
#      latent Gaussian value from the DCT rows -> push through the frozen read-time copula into
#      rho units. That is the posterior for region r GIVEN THE OTHER 63 REGIONS.
#   5. Convolve with the observation-noise model: map each rho draw through `fisherz`, add
#      N(0, 1/(n_eff - 3)), map back with `inv_fisherz`. The result is a PREDICTIVE for the
#      quantity actually observed -- region r's correlation entry -- not for its latent rho.
#   6a. Coverage: is the OBSERVED entry inside the equal-tailed `nominal` interval of those
#       predictive draws (`_p12_interval`)?
#   6b. Sharpness given calibration: a proper scoring rule -- the log predictive density, and
#       CRPS beside it.
#   7. Pool over regions and datasets, and compare arms under `lro_pass`.
#
# =============================================================================================
# THE HONEST SCOPE SENTENCE. D-09 VALIDATES A JOINT, NOT A POSTERIOR.
# =============================================================================================
#
# What step 6a measures is the calibration of the JOINT of (the posterior from step 4) AND (the
# observation-noise model of step 5). It is NOT a calibration statement about the posterior
# alone. A generous noise model widens every interval and can carry a MISCALIBRATED posterior to
# nominal coverage; nothing in the construction can see that, which is precisely why the
# simulated arm additionally reports PARAMETER coverage of the drawn latent value with no
# observation noise added. This is the same class of caveat Phase 11 attached to its hybrid
# quadrature column, and it is stated here rather than in a report because the report is where a
# caveat goes to be skipped.
#
# =============================================================================================
# PITFALL 6 -- A WORDING RULE WITH TEETH
# =============================================================================================
#
# D-10's phrase "independent pooling has no mechanism to predict a held-out region beyond the
# prior" is FALSE AS WRITTEN, and putting it in a report would be a falsifiable overclaim. With
# r1 -> 0 the REGIONS become conditionally independent, but the other 63 OBSERVED regions still
# inform the shared nuisances -- `chromatic_eps` above all, whose radial gradient is global --
# and the global DCT term c0. The ablation's predictive is therefore the prior conditional on the
# shared nuisances and the global term, which is STRICTLY NARROWER than the marginal prior.
#
# THE REQUIRED WORDING, EVERYWHERE, IS:
#     no spatial borrowing, full nuisance and global borrowing
#
# and a third, genuinely prior-only arm (`prior_only_floor`) is reported beside it so a reader
# can SEE how much of the ablation's performance is nuisance/global borrowing rather than spatial
# borrowing. The frozen `12-SC3-AMENDMENT.md` §4(c) already commits the phase to this wording, so
# a deviation here contradicts a committed document.
#
# =============================================================================================
# PITFALL 5 (12-RESEARCH.md:867 -- NOT the theta-un-standardization one) -- THE GATE IS TWO-SIDED
# =============================================================================================
#
# The predictive interval convolves posterior uncertainty about a region's rho with per-region
# OBSERVATION noise, measured at 0.0662 in correlation units against an across-region spread of
# only 0.0847 at constant rho -- a ratio of 1.28. When observation noise dominates, BOTH arms
# show approximately nominal coverage and "beats independent pooling in coverage" becomes a coin
# flip. Hence:
#
#   * CALIBRATION is a QUALIFYING condition, required of both arms. A miscalibrated arm is
#     DISQUALIFIED, not "worse".
#   * SHARPNESS GIVEN CALIBRATION -- a proper scoring rule -- is the DISCRIMINATING one, at the
#     pre-registered `P12_STAGE2_LOGSCORE_MIN` = 0.02 nats per held-out region.
#   * A COVERAGE-ONLY WIN WITH NO SCORING-RULE IMPROVEMENT IS EXPLICITLY NOT A PASS.
#
# `lro_pass` below is the ONE implementation of that rule, so it cannot drift between the runner
# and the report.
#
# THE OTHER PITFALL 5 IS ALSO LIVE HERE AND IS CITED BY FILE, NEVER BY NUMBER. The
# theta-un-standardization contract is `spike/npe/infer.jl:35-38`: `sampleposterior` returns theta
# in the STANDARDIZED space the estimator was trained in, and every read goes through
# `StatsBase.reconstruct(bundle.theta_zt, ...)` FIRST. 12-15's first run omitted it and produced a
# report that was complete, internally consistent, reconciled against its own parts, and
# scientifically void -- because no structural check looks at SCALE.
#
# =============================================================================================
# THE FIELD RECONSTRUCTION GOES THROUGH `p12_region_field`, AND THAT IS LOAD-BEARING
# =============================================================================================
#
# The theta rows are NOT the flat DCT coefficient vector. `p12_generate.jl:193` stores
# `p12_dct_vec(vec(z_field))[p12_dct_order(G)]` -- the coefficients PERMUTED into smoothness
# order -- and `p12_dct_order(8)` is not the identity (`[1, 2, 9, 10, 3, 17, 11, 18, ...]`).
# Applying `p12_idct_vec` directly to a theta column therefore reconstructs a DIFFERENT FIELD:
# verified, a round-trip that is exact to 3.1e-15 through `p12_region_field` is wrong by 3.17 in
# absolute value without the inverse permutation. `p12_region_field` (`spike/p12/result.jl:562`)
# scatters each coefficient back through `p12_dct_order` BEFORE the inverse transform, and it is
# the only reconstruction this file uses. Because a permutation is orthogonal, this error is not
# a small perturbation and no structural check can see it -- the same shape as the wrong-space
# episode above, one basis out.
#
# =============================================================================================
# WHY THIS FILE ADDS `draw_simulate_infer_p12` RATHER THAN CALLING `draw_simulate_infer`
# =============================================================================================
#
# `harness.jl:111`'s `draw_simulate_infer` is hardwired to the Phase-5 SCALAR lane: it calls
# `sample_prior` (7-row theta, a scalar rho_true), `standardize_summary(..., m.zt, :min)` and
# `posterior_for`. A Phase-12 draw is a FIELD plus 72 theta rows and standardizes through
# `standardize_p12`, so that function cannot express it and widening it would change a frozen
# Phase-5 surface. THE CHAIN IS NOT FORKED: `draw_simulate_infer_p12` calls exactly the same
# primitives, in the same order, reached through the same guarded includes --
# `simulate_pair` -> `build_mci` -> `patch_summary` -> `encode_d01` -- and adds only the field
# draw, the mask step and the reshape. There is still one simulate-and-encode path in the spike.
#
# THE STREAM IS THE **VALIDATION** ONE, AND THAT IS A POSITIVE PROPERTY WORTH CLAIMING. This file
# never reads a cached pool. It draws fresh through `p12_rng(P12_COVERAGE_COUNTER)`, which rides
# `P12_SALT`; pools ride `P12_DATAGEN_SALT`. The two are asserted disjoint at load time in both
# `p12_consts.jl` and `p12_generate.jl`, so a coverage number computed on training data is not
# merely avoided here -- THAT PATH DOES NOT EXIST.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. `src/` is reached only transitively and
# READ-ONLY, through `spike/p12/result.jl`'s existing include of `src/results.jl`. No package is
# added; every `using` below is already in `spike/Project.toml` or is a stdlib on the default
# LOAD_PATH.
#
# Run (self-check):
#     julia --project=spike -e 'include("spike/validation/p12_coverage.jl")'

using Statistics         # mean, std, var, quantile
using StatsBase          # reconstruct -- the INVERSE theta-standardization (infer.jl:35-38)
using Random             # AbstractRNG, randn
using Distributions      # Normal, cdf, quantile -- the frozen read-time copula

# --- ORDER MATTERS, AND THE FOREIGN PRE-REGISTRATION GOES FIRST -------------------------------
# `spike/p12/result.jl` is included FIRST because it already solves the include-guard poisoning
# problem for this dependency set (`result.jl:74-81`): it reaches `spike/validation/sbc.jl` ->
# `harness.jl` -> `consts.jl` (guarded on `:SBC_M`, which `consts.jl` alone owns) BEFORE it
# reaches the Phase-12 pre-registration. `p12_consts.jl` BINDS `:P11_DEV_SEED` and
# `:P13_DEV_SEED`, the sentinels the Phase-11 and Phase-13 blocks guard themselves on, so loading
# it first can silently skip a foreign one. Inheriting result.jl's order rather than restating it
# is what keeps this file out of that class -- which the executor briefing forbids fixing here,
# because it is coupled to an open user ruling.
#
# That one include transitively supplies: `SpatialColocResultSpike` and the whole 12-12 surface
# (`p12_region_maps`, `p12_mark_unscorable!`, `p12_result_meta`, `p12_region_field`,
# `P12_SPATIAL_SENTINEL`), the simulate/encode primitives (`simulate_pair`, `build_mci`,
# `patch_summary`, `encode_d01`), the frozen scalar priors (`MU_PRIOR`, `ghat`), the Tier-1
# pre-registration and the lattice arithmetic.
isdefined(@__MODULE__, :SpatialColocResultSpike) ||
    include(joinpath(@__DIR__, "..", "p12", "result.jl"))
# The Phase-12 training surface: `load_p12_npe` and the FROZEN `standardize_p12` (never refit).
# It transitively supplies `spike/data/p12_generate.jl` (`sample_p12_imsize`, `p12_theta_column`).
isdefined(@__MODULE__, :load_p12_npe) ||
    include(joinpath(@__DIR__, "..", "npe", "train_p12_npe.jl"))
# The field-aware prior draw: `sample_p12_prior`, `theta_scalar_view`, `rho_field`.
isdefined(@__MODULE__, :sample_p12_prior) ||
    include(joinpath(@__DIR__, "..", "simulator", "p12_prior.jl"))

# =============================================================================================
# 1. The Fisher-z transform, and the clamp that is a GUARD rather than a modelling choice
# =============================================================================================

"""
    P12_FISHERZ_CLAMP :: Float64

The magnitude `correlation()` output is clamped to before `atanh`.

**THIS IS A GUARD, NOT A MODELLING CHOICE, AND IT IS NOT A THRESHOLD.** `correlation()` can
return exactly `±1.0` on a degenerate patch -- two constant channels, or a patch reduced to a
handful of surviving pixels -- and `atanh(±1)` is infinite, which would poison every pooled
variance, log score and CRPS downstream with a `NaN`. Clamping keeps the arithmetic finite.

It is deliberately NOT in `p12_consts.jl`: only `P12_FISHERZ_NEFF` is authorised to be appended
there, and a numerical guard against a floating-point infinity is not a pre-registered bar. It
gates nothing.

**A CLAMPED ENTRY IS A DEGENERATE REGION, NOT A MEASUREMENT.** Every clamp is counted
([`count_clamped`](@ref)) and the count is a mandatory report key, because a silently clamped
`±1` correlation entering a coverage tally would be an unscorable region reported as a scored
one. `1 - 1e-6` puts `|z| = 7.25`, far outside any predictive interval this construction
produces, so a clamped entry cannot be quietly covered either.
"""
const P12_FISHERZ_CLAMP = 0.999999

"""
    fisherz(r) -> Float64

`atanh(clamp(r, -P12_FISHERZ_CLAMP, P12_FISHERZ_CLAMP))` -- the variance-stabilizing transform
of a correlation. See [`P12_FISHERZ_CLAMP`](@ref) for why the clamp exists and why a clamped
entry must be counted rather than swallowed.

**WHY FISHER-z RATHER THAN AN EMPIRICAL NOISE SURFACE.** Three reasons, recorded here because the
alternative is the obvious one and would be defensible if the reasons were not written down:

  1. **ONE calibrated constant instead of a two-dimensional surface.** An empirical model needs
     the observation-noise sd as a function of (image size, |rho|); Fisher-z needs a single
     `n_eff` per image size, and nothing else.
  2. **It is the standard treatment** of a correlation's sampling distribution, so it is a
     modelling assumption a referee can recognize rather than one this project invented.
  3. **The |rho| dependence becomes AUTOMATIC.** Research measured the per-region noise sd at
     0.0799 at rho ~= 0.16 against 0.0393 at rho ~= -0.82 -- a factor of 2.03, against the
     `(1 - r^2)` shape Fisher-z absorbs exactly (`(1-0.16^2)/(1-0.82^2) = 2.99` in variance,
     `1.73` in sd). An empirical surface would have to LEARN that shape from finite draws; the
     transform gets it for free and therefore extrapolates honestly to |rho| values the
     calibration draws did not visit.
"""
fisherz(r::Real) = atanh(clamp(float(r), -P12_FISHERZ_CLAMP, P12_FISHERZ_CLAMP))

"""
    inv_fisherz(z) -> Float64

`tanh(z)`: the exact inverse of [`fisherz`](@ref) on the unclamped range, mapping back into
correlation units -- the units the observed summary entry is in.
"""
inv_fisherz(z::Real) = tanh(float(z))

"""
    is_clamped_correlation(r) -> Bool

Whether `r` would be (or already sits at) the [`P12_FISHERZ_CLAMP`](@ref) boundary, i.e. whether
this entry is a DEGENERATE region rather than a measurement.
"""
is_clamped_correlation(r::Real) = abs(float(r)) >= P12_FISHERZ_CLAMP

"""
    count_clamped(v) -> Int

How many entries of `v` are degenerate by [`is_clamped_correlation`](@ref). Reported, never
swallowed: `clamped_region_count` is a mandatory artifact key.
"""
count_clamped(v) = count(is_clamped_correlation, v)

# =============================================================================================
# 2. The frozen read-time copula: latent Gaussian -> rho
# =============================================================================================

"""
    p12_z_to_rho(z) -> Float64

The FROZEN read-time map from the latent Gaussian field value to rho units:

    rho = ghat(quantile(MU_PRIOR, cdf(Normal(), z)))

**THIS IS NOT A NEW MAP.** It is the same composition, applied elementwise, that
`rho_field` uses to build the simulated field (`spike/simulator/p12_prior.jl:145-147`) and that
`theta_scalar_view` uses for the derived scalar (`:224`). It is written once here because this
file needs it per REGION per DRAW, and re-deriving the composition at each call site is how the
copula would silently fork from the one the simulator drew through.

`ghat` clamps rho to `±0.99`, which places a measured 6.79 % atom mass on every region
(12-RESEARCH.md Pattern 2). Those atoms are a property of the MODEL, not a degeneracy of the
data, and they are far inside [`P12_FISHERZ_CLAMP`](@ref) -- so a `ghat` atom never triggers the
degenerate-region path, and the two must not be conflated.
"""
p12_z_to_rho(z::Real) = ghat(quantile(MU_PRIOR, cdf(Normal(), float(z))))

"""
    p12_field_rho_draws(M; G = P12_G) -> Matrix{Float64}

Turn an UN-STANDARDIZED theta draw matrix `M` (`D x N`, physical units) into a `G^2 x N` matrix
of per-region rho draws.

Each column's first `G^2` rows are `[c0; the deviation coefficients in p12_dct_order]`. They go
through [`p12_region_field`](@ref) -- which INVERTS the smoothness permutation before the inverse
DCT, and is the only correct reconstruction in this repository (see the header) -- and then
elementwise through [`p12_z_to_rho`](@ref).

**`M` MUST ALREADY BE UN-STANDARDIZED.** This function cannot check that: a standardized theta
matrix has exactly the same shape and finiteness, and reconstructs into a perfectly plausible
field. The caller passes `StatsBase.reconstruct(bundle.theta_zt, ...)`; the contract is
`spike/npe/infer.jl:35-38`.
"""
function p12_field_rho_draws(M::AbstractMatrix; G::Integer = P12_G)
    nreg = G^2
    size(M, 1) >= nreg || throw(DimensionMismatch(
        "p12_field_rho_draws: expected at least $nreg theta rows (the c0 + deviation block at " *
        "G = $G), got $(size(M, 1))"))
    N   = size(M, 2)
    out = Matrix{Float64}(undef, nreg, N)
    @inbounds for d in 1:N
        F = p12_region_field(Float64.(view(M, 1:nreg, d)); G = G)   # G x G, GAUSSIAN space
        for r in 1:nreg
            out[r, d] = p12_z_to_rho(F[r])
        end
    end
    return out
end

"""
    p12_field_z_draws(M; G = P12_G) -> Matrix{Float64}

As [`p12_field_rho_draws`](@ref) but stopping in GAUSSIAN space -- the `G^2 x N` latent field
draws, before the copula.

This is what the simulated arm's PARAMETER-coverage cross-check scores against the drawn
`z_field`: Gaussian space is ATOM-FREE (R-2), while the elementwise `ghat` clamp puts 6.79 % atom
mass on every region in rho space, so a coverage statement made in rho space would be partly a
statement about the clamp.
"""
function p12_field_z_draws(M::AbstractMatrix; G::Integer = P12_G)
    nreg = G^2
    size(M, 1) >= nreg || throw(DimensionMismatch(
        "p12_field_z_draws: expected at least $nreg theta rows at G = $G, got $(size(M, 1))"))
    N   = size(M, 2)
    out = Matrix{Float64}(undef, nreg, N)
    @inbounds for d in 1:N
        F = p12_region_field(Float64.(view(M, 1:nreg, d)); G = G)
        for r in 1:nreg
            out[r, d] = F[r]
        end
    end
    return out
end

# =============================================================================================
# 3. `calibrate_neff` -- THE ONE CALIBRATION THIS FILE PERFORMS
# =============================================================================================

"""
    draw_simulate_infer_p12(rng; imsize, arm = :none, r1 = nothing, G = P12_G) -> NamedTuple

ONE Phase-12 draw-and-simulate step on the VALIDATION stream: draw a field-aware prior sample,
render it through the UNCHANGED forward model, and encode it through the UNCHANGED summary
contract. Returns `(draw, Zraw, imsize)` with `Zraw` the RAW `2G^2`-row encoding.

See the header for why this exists beside `harness.jl:111`'s `draw_simulate_infer` rather than
calling it: the Phase-5 function is hardwired to the scalar lane's 7-row theta and to
`standardize_summary`, so it cannot express a field draw. **The chain is not forked** -- the same
`simulate_pair` -> `build_mci` -> `patch_summary` -> `encode_d01` primitives are called here, in
the same order, reached through the same guarded includes. Only the field draw is added.

The body mirrors `generate_p12_sample` (`spike/data/p12_generate.jl:286`) with ONE deliberate
difference: that function keys its RNG on the global pool index and rides the DATAGEN salt, while
this one takes an `rng` off the VALIDATION stream. That difference is the whole reason this
plan's numbers cannot be computed on training data.
"""
function draw_simulate_infer_p12(rng::AbstractRNG; imsize, arm::Symbol = :none,
                                 r1 = nothing, G::Integer = P12_G)
    draw = sample_p12_prior(rng; G = G, arm = arm, r1 = r1)
    θs   = merge(theta_scalar_view(draw), (rho_field = draw.rho_field,))
    mci  = build_mci(simulate_pair(rng, θs; imsize = imsize))
    Zraw = encode_d01(patch_summary(mci))
    return (draw = draw, Zraw = Zraw, imsize = imsize)
end

"""
    calibrate_neff(; imsize, n_theta, n_obs, rng, arm = :none, G = P12_G) -> NamedTuple

Calibrate the ONE constant of the observation-noise model at ONE image size.

# THE DESIGN: A FIXED PARAMETER, REPEATEDLY OBSERVED -- REPLICATED OVER SEVERAL PARAMETERS

`n_theta` prior draws are taken. Each is ONE rho field plus one set of nuisances, and each is then
observed `n_obs` times through INDEPENDENT channel-pair simulations. **Within a theta block the
parameter is held constant and only the observation varies**, so the spread of the resulting
per-region correlation entries is OBSERVATION NOISE and nothing else.

**A DESIGN THAT REDREW THE FIELD EACH TIME WOULD BE WRONG, NOT MERELY NOISIER.** It would measure
the prior's spread CONVOLVED with the noise -- a different and much larger quantity -- and the
fitted `n_eff` would come out far too small, which would widen every predictive interval and buy
nominal coverage that means nothing. The conditioning is the whole design.

**WHY `n_theta > 1`, AND WHAT IT BUYS BEYOND PRECISION.** A single theta would make the constant
conditional on one particular rho field. Replicating over `n_theta` independent fields turns the
BETWEEN-THETA spread of the fitted variance into a direct test of the modelling claim that
motivates Fisher-z in the first place: if the transform really absorbs the `(1 - r^2)` dependence,
`Var(z)` should be approximately the same whatever field was drawn. `varz_by_theta` is returned
so that claim is CHECKABLE rather than asserted -- a large spread there would say the transform is
not doing its job at this image size, and that is a finding, not a nuisance.

Under the Fisher-z model `Var(z) = 1/(n_eff - 3)` per region, so `n_eff = 3 + 1/Var(z)`, pooled
over the `G^2` regions and the `n_theta` blocks. Regions whose MASK row reads 0 in a draw are
EXCLUDED from that draw -- an absent patch is not a noisy measurement of anything -- and clamped
entries are counted and reported rather than swallowed.

Returns `(n_eff, varz, varz_by_theta, varz_per_region, n_used, n_clamped, n_masked, imsize,
n_theta, n_obs)`.

**`n_eff` IS A DECLARED MODELLING ASSUMPTION, NOT A THRESHOLD.** It parametrizes the noise model
whose JOINT with the posterior D-09 validates. It gates nothing by itself, and no branch anywhere
compares a measured quantity against it.
"""
function calibrate_neff(; imsize, n_theta::Integer = 5, n_obs::Integer, rng::AbstractRNG,
                        arm::Symbol = :none, G::Integer = P12_G)
    n_obs >= 2 || throw(ArgumentError(
        "calibrate_neff: n_obs must be >= 2 to have a within-theta variance at all, got $n_obs"))
    n_theta >= 1 || throw(ArgumentError("calibrate_neff: n_theta must be >= 1, got $n_theta"))
    nreg = G^2

    varz_blocks = Float64[]              # one pooled Var(z) per theta block
    varz_all    = Float64[]              # every per-region, per-theta variance
    varz_per_region = zeros(Float64, nreg)
    varz_per_region_n = zeros(Int, nreg)
    n_clamped = 0
    n_masked  = 0

    for _ in 1:n_theta
        # ONE draw: one field, one set of nuisances. The block below re-observes THIS parameter.
        draw = sample_p12_prior(rng; G = G, arm = arm)
        θs   = merge(theta_scalar_view(draw), (rho_field = draw.rho_field,))

        zvals   = Matrix{Float64}(undef, nreg, n_obs)
        present = falses(nreg, n_obs)
        for d in 1:n_obs
            Zraw = encode_d01(patch_summary(build_mci(simulate_pair(rng, θs; imsize = imsize))))
            @inbounds for r in 1:nreg
                if Zraw[nreg + r] != 0.0                   # the binary present-mask row
                    present[r, d] = true
                    is_clamped_correlation(Zraw[r]) && (n_clamped += 1)
                    zvals[r, d] = fisherz(Zraw[r])
                end
            end
        end
        n_masked += count(!, present)

        block = Float64[]
        @inbounds for r in 1:nreg
            v = @view zvals[r, present[r, :]]
            length(v) >= 2 || continue
            vr = var(v)
            isfinite(vr) || continue
            push!(block, vr); push!(varz_all, vr)
            varz_per_region[r] += vr; varz_per_region_n[r] += 1
        end
        isempty(block) || push!(varz_blocks, mean(block))
    end

    isempty(varz_all) && error(
        "calibrate_neff: no region had >= 2 present observations at imsize = $imsize across " *
        "$n_theta theta block(s). The observation-noise model cannot be fitted from this " *
        "configuration; report the shortfall rather than substituting a default n_eff.")
    varz = mean(varz_all)
    varz > 0 || error(
        "calibrate_neff: pooled Var(z) = $varz is not positive at imsize = $imsize -- every " *
        "observation produced an identical correlation entry, which means the simulation was " *
        "not re-randomized between draws.")

    vpr = [varz_per_region_n[r] == 0 ? NaN : varz_per_region[r] / varz_per_region_n[r]
           for r in 1:nreg]
    return (n_eff           = 3 + 1 / varz,
            varz            = varz,
            varz_by_theta   = varz_blocks,
            varz_per_region = vpr,
            n_used          = length(varz_all),
            n_clamped       = n_clamped,
            n_masked        = n_masked,
            imsize          = imsize,
            n_theta         = Int(n_theta),
            n_obs           = Int(n_obs))
end

"""
    calibrate_neff_all(; n_theta, n_obs, rng, sizes = P12_IMSIZE_SET, arm = :none, G = P12_G)
        -> NamedTuple

Run [`calibrate_neff`](@ref) at every image size in `sizes` and pool.

**THE POOLED `n_eff` IS FITTED ON THE POOLED VARIANCE, NOT AVERAGED OVER THE PER-SIZE `n_eff`s.**
`n_eff` is a nonlinear function of the variance, so a mean of `n_eff`s is not the `n_eff` of the
mean variance, and the pooled value has to be the one whose implied variance a reader can check
against the realized per-size variances. `varz_fit_residual` is exactly that check: the largest
absolute discrepancy between a size's REALIZED `Var(z)` and the variance the POOLED constant
implies. It is REPORTED so the fit quality is visible rather than trusted -- it gates nothing, and
no branch reads it.

The pooling is unweighted across sizes because `P12_IMSIZE_SET` is the joint the estimator was
trained on and every size is scored; weighting by `P12_IMSIZE_WEIGHTS` would make the constant
track the mixture rather than the model, and would have to be re-derived if the mixture moved.
"""
function calibrate_neff_all(; n_theta::Integer = 5, n_obs::Integer, rng::AbstractRNG,
                            sizes = P12_IMSIZE_SET, arm::Symbol = :none, G::Integer = P12_G)
    per = [calibrate_neff(; imsize = sz, n_theta = n_theta, n_obs = n_obs, rng = rng,
                          arm = arm, G = G)
           for sz in sizes]
    varz_pooled  = mean(p.varz for p in per)
    n_eff_pooled = 3 + 1 / varz_pooled
    implied      = 1 / (n_eff_pooled - 3)
    return (n_eff        = n_eff_pooled,
            varz         = varz_pooled,
            by_imsize    = Dict(string(p.imsize) => p.n_eff for p in per),
            varz_by_imsize = Dict(string(p.imsize) => p.varz for p in per),
            # REPORTED, GATES NOTHING: the largest gap between a size's REALIZED Var(z) and the
            # variance the POOLED constant implies. The fit quality is visible, not trusted.
            varz_fit_residual = maximum(abs(p.varz - implied) for p in per),
            # The between-theta spread at each size: the direct check on Fisher-z's |rho|
            # -absorption claim (see `calibrate_neff`). A large value is a FINDING.
            varz_theta_spread_by_imsize =
                Dict(string(p.imsize) => (length(p.varz_by_theta) >= 2 ?
                                          std(p.varz_by_theta) : NaN) for p in per),
            n_clamped    = sum(p.n_clamped for p in per),
            n_masked     = sum(p.n_masked  for p in per),
            per_size     = per)
end

# =============================================================================================
# 4. Step 5: the predictive draws
# =============================================================================================

"""
    predictive_draws(rho_draws, n_eff; rng = Random.default_rng()) -> Vector{Float64}

Step 5 of the D-09 construction. Map each posterior rho draw for the held-out region through
[`fisherz`](@ref), add `randn * sqrt(1 / (n_eff - 3))`, and map back with [`inv_fisherz`](@ref).

Returns draws **in correlation units** -- the units the observed summary entry is in, which is
what makes the coverage comparison in step 6a a like-for-like one.

This is a CONVOLUTION, not a widening factor: the posterior spread and the observation noise
enter as independent contributions in `z` space, so a sharp posterior under a noisy observation
model and a vague posterior under a precise one are distinguishable here even though they can
produce the same interval width. That distinction is what the scoring rule in
[`lro_scores`](@ref) reads and raw coverage cannot.
"""
function predictive_draws(rho_draws, n_eff::Real; rng::AbstractRNG = Random.default_rng())
    n_eff > 3 || throw(ArgumentError(
        "predictive_draws: n_eff must exceed 3 for the Fisher-z variance 1/(n_eff-3) to be " *
        "positive, got $n_eff"))
    s = sqrt(1 / (n_eff - 3))
    return [inv_fisherz(fisherz(r) + s * randn(rng)) for r in rho_draws]
end

# =============================================================================================
# 5. Steps 6a and 6b: coverage, log score, CRPS -- and the WIDTH that makes coverage readable
# =============================================================================================

"""
    _p12cov_crps(pred_draws, observed) -> Float64

The empirical CRPS of a predictive SAMPLE against a scalar observation, in correlation units:

    CRPS = E|X - y| - 0.5 * E|X - X'|

with the second term computed exactly from the sorted sample as
`(2/n^2) * sum_i (2i - n - 1) * x_(i)`, which is `O(n log n)` rather than the `O(n^2)` double
loop. LOWER IS BETTER, so it is reported with the opposite sign convention from the log score --
stated here because reading the two in the same direction is an easy and silent error.
"""
function _p12cov_crps(pred_draws::AbstractVector, observed::Real)
    x = sort(collect(float.(pred_draws)))
    n = length(x)
    y = float(observed)
    t1 = mean(abs.(x .- y))
    s  = 0.0
    @inbounds for i in 1:n
        s += (2i - n - 1) * x[i]
    end
    t2 = (2 / n^2) * s
    return t1 - 0.5 * t2
end

"""
    lro_scores(observed, pred_draws; nominal = P12_COVERAGE_NOMINAL) -> NamedTuple

Steps 6a and 6b for ONE held-out region: score the observed correlation entry against that
region's predictive draws.

Returns `(inside, lo, hi, width, logscore, crps, z_mean, z_sd, clamped, n_draws)`.

# THE WIDTH IS RETURNED BESIDE THE COVERAGE, AND THAT IS MANDATORY

**Coverage is a property of interval WIDTH as much as of centre.** A predictive interval that is
wide and centred near the marginal mean covers ~= 90 % *while carrying no information about which
region is which*. So coverage alone cannot separate CALIBRATED AND INFORMATIVE from CALIBRATED
BECAUSE UNINFORMATIVE, and a coverage figure reported without its width cannot be read at all.
The same rule is why the runner reports `prior_only_floor` beside every log score: a number whose
null baseline is absent is uninterpretable, and this phase has been bitten by that twice.

# `inside` USES `_p12_interval`, WHICH IS NOT DEFINED HERE

The equal-tailed interval helper is `_p12_interval` from `spike/validation/p12_consts.jl` (12-01,
wave 1), which carries the three-line port and names `run_p11_coverage.jl:93-95` as its source.
**`run_p11_coverage.jl` IS NOT INCLUDED** to reach the original: it has no `*_LOAD_ONLY` guard and
ends in a bare `main()`, so an include to borrow a helper would EXECUTE the entire Phase-11
coverage run. There is exactly one definition, and it lives in the wave-1 file because 12-15
needed it two waves before this file existed.

# THE LOG SCORE IS CLOSED-FORM IN FISHER-z SPACE, AND THAT IS A CHOICE ABOUT INTEGRITY

The predictive draws are approximately Gaussian in `z` by construction -- step 5 adds Gaussian
noise in exactly that space -- so the log predictive density is evaluated as
`logpdf(Normal(mean(z), sd(z)), fisherz(observed))`, in closed form.

**THE ALTERNATIVE WAS A KDE, AND IT WAS REJECTED BECAUSE A BANDWIDTH IS A FREE PARAMETER.** A
bandwidth rule can be chosen -- consciously or not -- to favour an arm, and no artifact would
show it. The closed form has nothing to tune. (Had a KDE been used, `12-16-PLAN.md` requires its
bandwidth rule to be a pre-registered constant; it is not used, so no such constant exists.)

**THE CHANGE-OF-VARIABLES JACOBIAN IS OMITTED, DELIBERATELY AND VISIBLY.** The density is in
Fisher-z units, not correlation units; converting would add `-log(1 - observed^2)`. That term
depends ONLY on the observation, so it is IDENTICAL across every arm scored at the same held-out
region and cancels exactly in every difference. `P12_STAGE2_LOGSCORE_MIN` reads a DIFFERENCE, so
the gate is unaffected; but the absolute `logscore` figure is a Fisher-z-space density and must
not be quoted as a correlation-space one. `crps` is reported alongside precisely because it IS in
correlation units, is bounded, and is less sensitive to a single unlucky tail draw.
"""
function lro_scores(observed::Real, pred_draws::AbstractVector;
                    nominal::Real = P12_COVERAGE_NOMINAL)
    @assert length(pred_draws) >= 2 "lro_scores: need at least 2 predictive draws to form an " *
        "interval and a spread; got $(length(pred_draws))"
    y      = float(observed)
    lo, hi = _p12_interval(pred_draws, nominal)
    zs     = fisherz.(pred_draws)
    m      = mean(zs)
    s      = std(zs)
    s > 0 || (s = eps(Float64))     # a degenerate predictive is a defect, not a certainty; the
                                    # floor keeps the score finite so it is REPORTED, not thrown
    zy = fisherz(y)
    ls = -0.5 * log(2π) - log(s) - 0.5 * ((zy - m) / s)^2
    return (inside   = (lo <= y <= hi),
            lo       = lo,
            hi       = hi,
            width    = hi - lo,
            logscore = ls,
            crps     = _p12cov_crps(pred_draws, y),
            z_mean   = m,
            z_sd     = s,
            clamped  = is_clamped_correlation(y),
            n_draws  = length(pred_draws))
end

# =============================================================================================
# 6. `lro_pass` -- THE TWO-SIDED RULE, IMPLEMENTED ONCE
# =============================================================================================

"""
    lro_pass(spatial, ablation; nominal = P12_COVERAGE_NOMINAL,
             tost_delta = P12_STAGE2_COVERAGE_TOST_DELTA,
             logscore_min = P12_STAGE2_LOGSCORE_MIN) -> (verdict, details)

The pre-registered two-sided Stage-2 rule, as ONE function so it cannot drift between the runner
and the report. `spatial` and `ablation` are NamedTuples carrying at least `coverage` and
`mean_logscore`.

    1. BOTH arms must be CALIBRATED: `|coverage - nominal| <= tost_delta`.
       A miscalibrated arm is **DISQUALIFIED, not "worse"** -- returns `(:disqualified, ...)`
       naming the offending arm. The spatial arm is named first when both fail, and `details`
       records both flags so "both were miscalibrated" is never lost.
    2. Given both calibrated, require
       `mean_logscore(spatial) - mean_logscore(ablation) >= logscore_min`.
    3. Verdict is one of `:pass`, `:fail`, `:disqualified` -- a SYMBOL, never a bare boolean.

**A coverage-only win returns `:fail`.** That is the entire point of the rule: when observation
noise dominates, both arms sit near nominal and a coverage comparison becomes a coin flip, so
calibration is the QUALIFYING condition and the proper scoring rule is the DISCRIMINATING one.
The four branches of this function ARE the frozen SC3 criterion; a change that breaks one of them
is a change to a frozen success criterion, not a refactor.

# ON THE D-13 DESCOPE ROUTE THIS FUNCTION CANNOT FIRE, BY CONSTRUCTION

It compares a spatial arm against the ablation, and on the route this phase took **there is no
spatial arm**: the mini-spike returned NONE-BEATS-ABLATION and 12-17 never ran. The runner
therefore records `stage2_gate = :not_applicable_descope` and never calls this function on real
data. It is built anyway -- it is small next to `p12_coloc_map`, it keeps the artifact contract
intact for 12-20's Guard 3, and it is the first thing a v2.1 spatial arm needs. **Its tests
exercise the BRANCH LOGIC on synthetic inputs, not the gate. A green test here must never be read
as the gate having run.**
"""
function lro_pass(spatial, ablation;
                  nominal::Real = P12_COVERAGE_NOMINAL,
                  tost_delta::Real = P12_STAGE2_COVERAGE_TOST_DELTA,
                  logscore_min::Real = P12_STAGE2_LOGSCORE_MIN)
    cal_s = abs(spatial.coverage  - nominal) <= tost_delta
    cal_a = abs(ablation.coverage - nominal) <= tost_delta
    delta = spatial.mean_logscore - ablation.mean_logscore
    base  = (spatial_calibrated = cal_s, ablation_calibrated = cal_a,
             spatial_coverage = spatial.coverage, ablation_coverage = ablation.coverage,
             logscore_delta = delta, nominal = nominal, tost_delta = tost_delta,
             logscore_min = logscore_min)

    if !cal_s || !cal_a
        # The spatial arm is named first when BOTH fail; `details` keeps both flags, so a
        # both-miscalibrated run is never reported as a single-arm problem.
        which = !cal_s ? :spatial : :ablation
        return (:disqualified, merge(base, (which_arm = which,
            reason = "arm :$which is MISCALIBRATED (|coverage - $nominal| > $tost_delta). A " *
                     "miscalibrated arm is DISQUALIFIED, not 'worse' -- its scoring-rule value " *
                     "is not comparable, because sharpness is only meaningful GIVEN calibration.")))
    end
    if delta >= logscore_min
        return (:pass, merge(base, (which_arm = nothing,
            reason = "both arms calibrated AND the spatial arm beats the ablation by " *
                     "$(delta) >= $logscore_min nats per held-out region.")))
    end
    return (:fail, merge(base, (which_arm = nothing,
        reason = "both arms calibrated but the spatial log-score advantage is $(delta) < " *
                 "$logscore_min nats per held-out region. A COVERAGE-ONLY WIN IS NOT A PASS.")))
end

# =============================================================================================
# 7. `prior_only_floor` -- the third arm Pitfall 6 requires
# =============================================================================================

"""
    prior_only_floor(; N, rng, G = P12_G) -> Matrix{Float64}

The THIRD arm: predict the held-out region from its MARGINAL prior alone -- draws from
`MU_PRIOR` pushed through the frozen `ghat`, **with no conditioning on the image at all**.
Returns a `G^2 x N` matrix so it plugs into the same scoring path as the model arms; every region
carries an independent draw set, because the marginal prior is identical across regions but the
Monte-Carlo noise should not be shared.

**THIS IS THE FLOOR THAT MAKES THE ABLATION READABLE.** The D-10 ablation is *no spatial
borrowing, full nuisance and global borrowing* -- its predictive is the prior CONDITIONAL on the
shared nuisances and the global term, which is strictly narrower than the marginal prior. Without
this arm a reader cannot tell how much of the ablation's performance is nuisance/global borrowing
rather than spatial borrowing, and the Pitfall-6 overclaim becomes unfalsifiable rather than
merely wrong.

It is also the DISCRIMINATOR the 12-16 pre-declaration fixed before any number existed: coverage
alone cannot separate *calibrated and informative* from *calibrated because uninformative*, and
the margin over this floor is what does.

**WHY `MU_PRIOR . ghat` IS EXACTLY THE PER-REGION MARGINAL.** The D-05 copula pushes each
region's latent Gaussian through `Phi`, then `quantile(MU_PRIOR, .)`, then `ghat` -- so every one
of the `G^2` regions carries EXACTLY the prior `spike/simulator/prior.jl:40` puts on the scalar
knob. The lattice prior contributes only spatial DEPENDENCE, never a different marginal. Drawing
`ghat.(rand(rng, MU_PRIOR, N))` is therefore the region's marginal prior, not an approximation of
it.
"""
function prior_only_floor(; N::Integer, rng::AbstractRNG, G::Integer = P12_G)
    N >= 2 || throw(ArgumentError("prior_only_floor: N must be >= 2, got $N"))
    nreg = G^2
    out  = Matrix{Float64}(undef, nreg, N)
    @inbounds for r in 1:nreg, d in 1:N
        out[r, d] = ghat(rand(rng, MU_PRIOR))
    end
    return out
end

# =============================================================================================
# 8. `lro_arm` -- steps 1 to 4 for one arm
# =============================================================================================

"""
    lro_arm(bundle, Zraw, region_range; N, G = P12_G) -> NamedTuple

Steps 1-4 of the D-09 construction for ONE arm on ONE dataset: for every region `r` in
`region_range`, mask `r`, read the posterior, and return region `r`'s own posterior draws.

Returns `(rho = length(region_range) x N, z = length(region_range) x N, regions = ...)`: the
held-out region's draws in rho units AND in latent Gaussian space, because the simulated arm's
parameter-coverage cross-check is scored in the atom-free Gaussian space (R-2) while the
predictive comparison is scored in correlation units.

# MASKING IS APPLIED TO THE **RAW** ROWS

`mask_regions` runs on `Zraw` BEFORE `standardize_p12`, which is the ordering
`train_p12_npe.jl:30-42` establishes and 12-14's testset 3 pins. The reason is the whole of
Pitfall 4: a naturally unusable patch reaches the net as raw value `0.0` with mask `0`, which
after per-row z-scoring becomes `(0 - mu_r)/sd_r` -- NOT zero. Masking after standardization would
instead place the ROW MEAN at the held-out region, so the net would read *"this region is
perfectly average"* exactly where read time says *"this region is absent"*. That is an
out-of-distribution encoding and a silent one. It is asserted again here, at the CONSUMER, in
this file's testset 6 -- deliberately twice, because the producer and the consumer can drift
apart independently.

# ONE BATCHED `sampleposterior` CALL, NOT ONE PER REGION

The `length(region_range)` masked variants are assembled into a single `(G, G, 2, K)` batch, so
the whole leave-one-region-out sweep of a dataset is ONE forward pass set rather than 64. That is
a performance property with a correctness consequence worth naming: every region of a dataset is
scored by the same net in the same call, so no between-region difference can come from batching.

# EVERY DRAW IS UN-STANDARDIZED FIRST

`StatsBase.reconstruct(bundle.theta_zt, ...)` is applied before any row is touched
(`spike/npe/infer.jl:35-38`). CPU-only: `use_gpu = false` is passed EXPLICITLY, because
NeuralEstimators v0.2.1 defaults it to `true`.
"""
function lro_arm(bundle, Zraw::AbstractVector, region_range; N::Integer, G::Integer = P12_G)
    nreg = G^2
    regs = collect(region_range)
    all(r -> 1 <= r <= nreg, regs) || throw(ArgumentError(
        "lro_arm: every held-out region must lie in 1:$nreg (G = $G), got $(extrema(regs))"))
    length(Zraw) == 2nreg || throw(DimensionMismatch(
        "lro_arm: expected $(2nreg) RAW summary rows, got $(length(Zraw))"))

    # STEP 2, ON THE RAW ROWS -- one masked column per held-out region.
    Zmasked = Matrix{Float64}(undef, 2nreg, length(regs))
    @inbounds for (t, r) in enumerate(regs)
        Zmasked[:, t] = mask_regions(Float64.(Zraw), r; G = G)
    end
    # STEP 3: standardize with the FROZEN transform, THEN reshape. Never refit.
    Zb = reshape_summary(standardize_p12(Zmasked, bundle.zt; G = G), G)
    # STEP 4: one batched posterior read, un-standardized before anything reads a row.
    smp  = sampleposterior(bundle.estimator, Zb; N = N, use_gpu = false)
    raws = smp isa AbstractVector ? smp : [smp]
    length(raws) == length(regs) || error(
        "lro_arm: sampleposterior returned $(length(raws)) draw sets for $(length(regs)) " *
        "masked datasets. The batch axis and the region axis must correspond one-to-one.")

    rho = Matrix{Float64}(undef, length(regs), N)
    z   = Matrix{Float64}(undef, length(regs), N)
    for (t, r) in enumerate(regs)
        M = StatsBase.reconstruct(bundle.theta_zt, Float64.(raws[t]))   # infer.jl:35-38
        zd = p12_field_z_draws(M; G = G)
        @inbounds for d in 1:N
            z[t, d]   = zd[r, d]
            rho[t, d] = p12_z_to_rho(zd[r, d])
        end
    end
    return (rho = rho, z = z, regions = regs)
end

# =============================================================================================
# 9. `p12_coloc_map` -- SPAT-05, the read surface 12-12 deferred to this plan
# =============================================================================================

"""
    p12_coloc_map(bundle, mci, control_mci::Union{Nothing,MultiChannelImage};
                  N, nominal, n_eff) -> SpatialColocResultSpike

Produce the per-region colocalization map: **all three** maps (`region_rho_sample`,
`region_rho_control` and the PRIMARY `region_delta_rho`) each with its per-region uncertainty,
assembled into the 12-12 result type.

`control_mci` IS POSITIONAL AND MANDATORY. `N`, `nominal` and `n_eff` are keywords with NO
DEFAULTS -- every call site states them.

# TWO INDEPENDENT SINGLE-STACK PASSES, NOT A JOINT ESTIMATOR

For each of `mci` and `control_mci`, separately and identically:

    encode_d01(patch_summary(.)) -> standardize_p12 with the FROZEN `bundle.zt`
    -> reshape_summary -> sampleposterior -> reconstruct(bundle.theta_zt, .)
    -> p12_region_field (which INVERTS the smoothness permutation) -> p12_z_to_rho

and then the per-region MONTE-CARLO DIFFERENCE of the two draw sets. This mirrors
`src/amortized/infer.jl:108` rather than inventing a joint estimator (D-03): there are TWO
`sampleposterior` calls in this function and no third path.

**EXPOSING THE TWO SINGLE-STACK MAPS COSTS ZERO EXTRA FORWARD PASSES.** Both draw sets already
exist before the subtraction -- a Delta-rho read performs both passes internally regardless -- so
the marginal cost of shipping them is artifact keys, not compute. A reader should not conclude
the scope doubled.

# `control_mci === nothing` IS A SANCTIONED, EXPLICIT SINGLE-STACK MODE

It means exactly one thing: *there is no control stack, therefore there is no Delta-rho here.*
The real-image arm is why it exists -- 12-19 must report `region_rho_sample` on two physical
specimens and is simultaneously forbidden to pair them, because the simulated pair is two
exchangeable draws from one prior (a WITHIN-population difference) while two real anchors are
deliberately non-exchangeable (a BETWEEN-population one). With `control_mci` mandatory those two
requirements are unsatisfiable together.

In that mode the sample pass runs alone and the result is built with `region_rho_control`
entirely `P12_SPATIAL_SENTINEL`, `region_sd_control` entirely `NaN` and `ood_control` entirely
`true`; 12-12's union rule then forces `ood` all-true and `region_delta_rho` all-sentinel. **That
is the honest encoding: the named Delta-rho limit becomes a machine-checkable property of the
returned object rather than a sentence in a report.** No constructor change is needed -- the
structural identity is scoped to regions where no map is a sentinel, so it holds vacuously, and
`P12_SPATIAL_SENTINEL = 0.0` sits inside both sanity ranges. `meta.delta_rho_available` records
the mode by name so a consumer never infers it from a sentinel pattern. **A single-stack result
is NOT a Delta-rho of zero and must never be summarized as one.**

`control_mci` stays MANDATORY AS AN ARGUMENT: passing `nothing` must be a written decision at the
call site, never a default a caller reaches by omission.

# DEGENERATE REGIONS GO THROUGH 12-12's SINGLE SENTINEL WRITER

A region is unscorable in a stack when its MASK row reads 0 (the patch fell below the
>=15-surviving-pixel floor) **or** its correlation entry is clamped at `±1`
([`is_clamped_correlation`](@ref)). Those regions are marked with `p12_mark_unscorable!`, which
writes the sentinel TRIPLE -- value, `NaN` sd, `true` flag -- and propagates through the union
rule in one call. **A sentinel is therefore never written without its OOD flag**, and a genuine
measured Delta-rho of exactly `0.0` on a non-degenerate region carries `ood = false` and remains
distinguishable from a refusal (T-7-04).

# THE SCALAR `rho` THIS RESULT REPORTS IS **DERIVED** (R-1)

`meta.derived_rho` is `ghat(quantile(MU_PRIOR, Phi(c0 / G)))` computed from the SAMPLE stack's
posterior-mean global DCT term. It is NOT a sampled parameter, and the label is in the field name,
in this docstring and in `delta_rho`'s own docstring. **The per-region map is the deliverable, not
this scalar.**

# THE OOD OPERATING POINT

`meta.ood_threshold` is left `nothing` unless one is supplied, so a `false` flag reads as **NOT
CHECKED** rather than as "in distribution" (`src/amortized/local_map.jl:88-91`). Branch on
`p12_ood_checked` before reading an all-`false` OOD map as a clean bill of health.

# WHAT `nominal` AND `n_eff` DO HERE, STATED PLAINLY RATHER THAN IMPLIED

`12-16-PLAN.md:201` mandates this signature, so both are required and neither has a default. **But
this function reports the posterior sd per region, not a predictive interval, so it consumes them
only to VALIDATE and to refuse an inconsistent call.** They are the observation-noise model's
settings, and the predictive construction that uses them is [`predictive_draws`](@ref) /
[`lro_scores`](@ref) in the coverage path. Requiring them here means a map and the coverage
numbers computed beside it cannot be built under silently different settings -- but a reader
should not infer from the signature that a predictive interval is stored in the returned object.
It is not. Recorded rather than quietly dropped: the arguments are the plan's, the honesty about
what they do is this file's.
"""
function p12_coloc_map(bundle, mci, control_mci;
                       N::Integer, nominal::Real, n_eff::Real,
                       G::Integer = P12_G, ood_threshold = nothing, net_path = nothing)
    N >= 2 || throw(ArgumentError("p12_coloc_map: N must be >= 2, got $N"))
    0 < nominal < 1 || throw(ArgumentError(
        "p12_coloc_map: nominal must lie in (0,1), got $nominal"))
    n_eff > 3 || throw(ArgumentError("p12_coloc_map: n_eff must exceed 3, got $n_eff"))
    nreg = G^2

    # --- The SAMPLE stack: pass 1 of 2 ---------------------------------------------------------
    smp = _p12cov_single_stack(bundle, mci; N = N, G = G)
    # --- The CONTROL stack: pass 2 of 2, or the sanctioned refusal -----------------------------
    ctl = control_mci === nothing ? nothing : _p12cov_single_stack(bundle, control_mci; N = N, G = G)
    delta_rho_available = ctl !== nothing

    # FAIL-CLOSED: every region is a refusal until something measures it (12-12's own words).
    maps = p12_region_maps(G)

    for r in 1:nreg
        s_ok = !smp.degenerate[r]
        c_ok = ctl !== nothing && !ctl.degenerate[r]
        if s_ok
            maps.region_rho_sample[r] = mean(view(smp.rho, r, :))
            maps.region_sd_sample[r]  = std(view(smp.rho, r, :))
            maps.ood_sample[r]        = false
        end
        if c_ok
            maps.region_rho_control[r] = mean(view(ctl.rho, r, :))
            maps.region_sd_control[r]  = std(view(ctl.rho, r, :))
            maps.ood_control[r]        = false
        end
        if s_ok && c_ok
            d = view(smp.rho, r, :) .- view(ctl.rho, r, :)   # the per-region MC difference
            maps.region_delta_rho[r] = mean(d)
            maps.region_sd[r]        = std(d)
            maps.ood[r]              = false
        end
    end
    # The regions that were ATTEMPTED and FAILED go through the ONE sentinel writer, so the
    # sentinel/NaN/flag triple and the union propagation cannot come apart. Regions never
    # measured at all are already in that state from `p12_region_maps`.
    p12_mark_unscorable!(maps, findall(smp.degenerate); stack = :sample)
    ctl === nothing ? p12_mark_unscorable!(maps, 1:nreg; stack = :control) :
                      p12_mark_unscorable!(maps, findall(ctl.degenerate); stack = :control)

    # A zero-sd region is a defect the constructor would reject with a confusing message; catch
    # it here and say what it means. It can only arise from a degenerate posterior.
    for (nm, sd, fl) in (("region_sd_sample",  maps.region_sd_sample,  maps.ood_sample),
                         ("region_sd_control", maps.region_sd_control, maps.ood_control),
                         ("region_sd",         maps.region_sd,         maps.ood))
        for i in eachindex(sd)
            fl[i] && continue
            (isfinite(sd[i]) && sd[i] > 0) || error(
                "p12_coloc_map: `$nm` is $(sd[i]) at a region flagged SCORABLE. Every posterior " *
                "draw collapsed to one value there, which is a degenerate read, not a confident " *
                "one. Refusing to construct a result that would report zero uncertainty.")
        end
    end

    # The G x G x N stack of DIFFERENCE draws. `nothing` in single-stack mode, because there is
    # no difference to carry -- an all-NaN stack would look like a computed-then-lost quantity.
    # Unscorable regions carry NaN rather than a plausible number, matching their NaN sd.
    draws = nothing
    if delta_rho_available
        draws = Array{Float64,3}(undef, G, G, N)
        @inbounds for d in 1:N, r in 1:nreg
            draws[r + (d - 1) * nreg] =           # column-major: flat r IS the (i,j) of region r
                maps.ood[r] ? NaN : (smp.rho[r, d] - ctl.rho[r, d])
        end
    end

    # The DERIVED (R-1) scalar, from the SAMPLE stack's posterior-mean global term.
    c0_mean = mean(view(smp.theta_c0, :))
    derived_rho = ghat(quantile(MU_PRIOR, cdf(Normal(), c0_mean / G)))
    r1_median = quantile(vec(smp.theta_r1), 0.5)

    meta = p12_result_meta(; arm = bundle.arm,
                           r1_posterior_median = r1_median,
                           derived_rho = derived_rho,
                           ood_threshold = ood_threshold,
                           net_path = net_path,
                           delta_rho_available = delta_rho_available)
    return SpatialColocResultSpike(G, maps, draws, nothing, meta)
end

"""
    _p12cov_single_stack(bundle, mci; N, G) -> NamedTuple

ONE single-stack posterior read: `encode_d01(patch_summary(mci))` -> `standardize_p12` with the
FROZEN `bundle.zt` -> `reshape_summary` -> `sampleposterior` -> `reconstruct(bundle.theta_zt, .)`
-> `p12_region_field` -> `p12_z_to_rho`.

Returns `(rho = G^2 x N, degenerate = G^2 Bool, theta_c0 = N, theta_r1 = N, Zraw)`.

`degenerate[r]` is true when region `r`'s MASK row reads 0 (an unusable patch) or its correlation
entry is clamped at `±1`. Those are the two ways a region can be present in the array but absent
as a measurement, and both must reach the sentinel writer.

CPU-only: `use_gpu = false` is EXPLICIT (v0.2.1 defaults it to `true`).
"""
function _p12cov_single_stack(bundle, mci; N::Integer, G::Integer = P12_G)
    nreg = G^2
    Zraw = Float64.(encode_d01(patch_summary(mci)))
    length(Zraw) == 2nreg || throw(DimensionMismatch(
        "_p12cov_single_stack: expected $(2nreg) encoded rows at G = $G, got $(length(Zraw))"))

    Z    = reshape_summary(standardize_p12(reshape(Zraw, :, 1), bundle.zt; G = G), G)
    smp  = sampleposterior(bundle.estimator, Z; N = N, use_gpu = false)
    raws = smp isa AbstractVector ? smp : [smp]
    M    = StatsBase.reconstruct(bundle.theta_zt, Float64.(raws[1]))   # infer.jl:35-38

    degenerate = falses(nreg)
    @inbounds for r in 1:nreg
        degenerate[r] = (Zraw[nreg + r] == 0.0) || is_clamped_correlation(Zraw[r])
    end
    return (rho        = p12_field_rho_draws(M; G = G),
            degenerate = degenerate,
            theta_c0   = vec(M[1, :]),
            theta_r1   = vec(M[end, :]),
            Zraw       = Zraw)
end

# =============================================================================================
# 10. Pooling helpers: the coverage interval, ON THE INDEPENDENT UNIT
# =============================================================================================

"""
    _p12cov_wilson(k, n; z = Z_TWO_SIDED_90) -> (lo, hi)

Wilson score interval for a binomial proportion `k/n`.

**COPIED WITH ATTRIBUTION from `spike/validation/p11_stats.jl:49-64` rather than included**, which
is the same discipline `p11_stats.jl` itself records for the two functions it copied out of
`test/gate/`. `p11_stats.jl` guards its own `p11_consts.jl` include on `:P11_DEV_SEED` -- a
sentinel `p12_consts.jl` BINDS -- so including it from a Phase-12 runner is a textbook instance of
the include-guard poisoning the executor briefing forbids fixing here. Four lines of arithmetic
are cheaper than a poisoned load, and `Z_TWO_SIDED_90` is read from `p12_consts.jl`.

Wald is deliberately not used for the INTERVAL: near `p = 0.90` with `n` in the low hundreds its
actual coverage is poor and it can extend past 1. (`P12_STAGE2_N_MIN` was nonetheless SIZED with
the Wald half-width, which `p12_consts.jl:416-422` records verbatim and labels precisely; sizing
and reporting are two different constructions and this file does not collapse them.)
"""
function _p12cov_wilson(k::Integer, n::Integer; z::Real = Z_TWO_SIDED_90)
    p̂ = k / n
    d = 1 + z^2 / n
    c = (p̂ + z^2 / (2n)) / d
    h = (z / d) * sqrt(p̂ * (1 - p̂) / n + z^2 / (4n^2))
    return (lo = c - h, hi = c + h)
end

"""
    _p12cov_cluster_interval(per_dataset_fractions; z = Z_TWO_SIDED_90) -> (lo, hi, half_width)

**THE INTERVAL THAT MAY BE QUOTED.** `per_dataset_fractions[i]` is the fraction of dataset `i`'s
regions that were covered. The interval is `mean +/- z * sd / sqrt(n_datasets)` -- a
BETWEEN-CLUSTER interval whose `n` is the number of DATASETS.

**THE INDEPENDENT UNIT IS THE DATASET, AND THIS MILESTONE HAS ALREADY MADE THE OTHER MISTAKE
ONCE.** At 64 regions per image a run at `N = P12_STAGE2_N_MIN = 271` produces `271 * 64 = 17,344`
region-draws, **but only 271 independent units**: the 64 regions of one image share every nuisance
and the global term, so they are PSEUDO-REPLICATES (`12-19-PLAN.md:148-150`, `12-02-PLAN.md:235`).
An interval built on 17,344 would be roughly 8x too tight and would manufacture significance out
of pseudo-replication -- which is exactly what the real-image arm did with a gate sized for
`N >= 271` applied to an arm whose `effective_independent_n` is 2.

This construction takes the between-dataset spread as the standard error, so it absorbs the
within-image correlation automatically instead of assuming it away. `effective_independent_n` is a
mandatory artifact key beside `n_region_draws`, and the 17,344 figure is NEVER reported alone.
"""
function _p12cov_cluster_interval(per_dataset_fractions::AbstractVector; z::Real = Z_TWO_SIDED_90)
    n = length(per_dataset_fractions)
    n >= 2 || throw(ArgumentError(
        "_p12cov_cluster_interval: need at least 2 datasets to have a between-cluster spread, " *
        "got $n. The independent unit is the DATASET, never the region-draw."))
    m = mean(per_dataset_fractions)
    h = z * std(per_dataset_fractions) / sqrt(n)
    return (lo = m - h, hi = m + h, half_width = h, n_independent = n)
end
