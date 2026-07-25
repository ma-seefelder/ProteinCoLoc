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

# spike/simulator/prior.jl --- SIM-02: the simulator prior π(θ), made provably
# consistent with the Turing μ-prior via the frozen induced-μ inverse ĝ.
#
# sample_prior draws μ*~Truncated(Cauchy(0,0.3),-1,1) (the VERBATIM Turing μ-prior,
# src/bayes.jl ~274-291) and sets ρ_true = ghat(μ*) so the per-pair INDUCED μ
# (measured through the UNCHANGED summary) matches that target by construction
# (transform sampling). The 6 nuisances are drawn from their own prior ranges, so
# the induced μ carries the real generative spread the calibration averaged over.
#
# Do NOT assert ρ_true == μ (02-PATTERNS Anti-Pattern / Pitfall 2): the nonneg
# softplus mixing + PSF + noise attenuate ρ_true into μ nonlinearly -- ĝ inverts
# exactly that measured map.

using Distributions          # Truncated, Cauchy, Uniform
using Random                 # AbstractRNG

include(joinpath(@__DIR__, "ghat.jl"))   # ghat, GHAT_MU_KNOTS/RHO_KNOTS, GHAT_MU_MIN/MAX

# --- The Turing μ-prior (src/bayes.jl ~274-291), the SIM-02 consistency target ---
const MU_PRIOR = Truncated(Cauchy(0.0, 0.3), -1.0, 1.0)

# --- PRE-DECLARED SIM-02 pass bar (T-02-CAL) -------------------------------------
# Identical to calibration.jl's tolerance; fixed BEFORE the calibration run, NOT
# tuned to pass. The induced-μ Wasserstein-1 distance to MU_PRIOR (scoped to the
# physically-realized range [GHAT_MU_MIN, GHAT_MU_MAX], D-16) must clear this.
const SIM02_W1_TOL = 0.10

# --- Nuisance prior ranges (Claude's discretion, D-03) ---------------------------
# MIRROR of calibration.jl's `_sample_nuisances` (the calibration measured E[μ|ρ]
# averaging over EXACTLY these ranges; documented in spike/NOTES.md §3). Chosen so
# the induced per-patch-correlation spread (→ σ/τ) stays inside the Turing scale
# prior Truncated(Cauchy(0.1,0.3),1e-4,1) -- verified: pooled SD ≈ 0.32.
const SPILLOVER_PRIOR        = Uniform(0.0, 0.2)   # directional bleed-through (modest)
const AUTOFLUORESCENCE_PRIOR = Uniform(0.0, 0.1)   # additive background offset
const LABEL_EFFICIENCY_PRIOR = Uniform(0.6, 1.0)   # Bernoulli keep-probability
const SHIFT_PRIOR            = Uniform(-3.0, 3.0)  # registration error (dx, dy); D-02 WIDENED, see below
const CHROMATIC_PRIOR        = Uniform(-0.02, 0.02) # radial chromatic scale ε (D-09), see below
const NOISE_PRIOR            = Uniform(0.0, 1.0)   # Poisson+Gaussian noise scale

# SHIFT_PRIOR (D-02) DELIBERATELY DIVERGES from src/amortized/simulator.jl, which STAYS at
# Uniform(-1.0, 1.0). D-02 scopes the ±3 px widening to the RESEARCH net only, and D-11 mirrors
# ONLY the chromatic term into src/. Mirroring the widening would desynchronise
# theta_prior_bounds()[5:6] from the shipped bundle's FROZEN theta_zt.lo/.hi — a silent semantic
# drift on the reconstruct path. ±3 px covers realistic multi-channel error (chromatic aberration
# plus stage/filter-cube drift) and gives the SC2 uncertainty ladder real headroom (0.25 → 3 px).
#
# CHROMATIC_PRIOR (D-09) is a 1-parameter RADIAL MAGNIFICATION difference between the channels:
# zero at the image centre, maximal at the corners, so a global dx/dy translation cannot represent
# it — a distinct failure mode, not a reparameterisation. On a 256² field the half-diagonal is
# ≈181 px, so corner displacement ≈ ε·181: 0.005 → ≈0.9 px (well-corrected apochromat), 0.01 →
# ≈1.8 px (typical achromat), 0.02 → ≈3.6 px (poorly corrected) — magnitude-comparable to the ±3 px
# translation headroom. SYMMETRIC because channel ordering is the USER's choice, not the model's:
# a one-sided prior would silently assume ch2 is always the longer wavelength.

"""
    sample_prior(rng::AbstractRNG) -> NamedTuple

Draw one θ from the simulator prior π(θ) (the 8-field NamedTuple: the 7 D-04
fields plus `chromatic_eps`, D-09). The colocalization knob is set through the
frozen induced-μ inverse:

  μ* ~ Truncated(Cauchy(0,0.3),-1,1)      # the Turing μ-prior (MU_PRIOR)
  ρ_true = ghat(μ*)                        # so E[induced μ | ρ_true] = μ*

and the 7 nuisances are drawn from their own priors. Deterministic under a fixed
`rng` (D-14). The result flows straight into `simulate_pair(rng, θ; imsize)`.

`chromatic_eps` is APPENDED AT THE END of the tuple, never slotted next to
`shift_dx`/`shift_dy`, because downstream code indexes θ POSITIONALLY
(`collect(values(θ))` row order; `src/amortized/ood.jl` `_theta_tuple` reads
`v[1]..v[7]`). Appending keeps every existing index valid.
"""
function sample_prior(rng::AbstractRNG)
    μ_star = rand(rng, MU_PRIOR)
    return (
        ρ_true           = ghat(μ_star),
        spillover        = rand(rng, SPILLOVER_PRIOR),
        autofluorescence = rand(rng, AUTOFLUORESCENCE_PRIOR),
        label_efficiency = rand(rng, LABEL_EFFICIENCY_PRIOR),
        shift_dx         = rand(rng, SHIFT_PRIOR),
        shift_dy         = rand(rng, SHIFT_PRIOR),
        noise            = rand(rng, NOISE_PRIOR),
        chromatic_eps    = rand(rng, CHROMATIC_PRIOR),   # D-09: appended LAST (positional θ)
    )
end
