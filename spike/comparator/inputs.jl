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

# spike/comparator/inputs.jl --- CMP-03: seeded shared-input builder (D-02/D-03).
#
# The comparator's shared-input SOURCE: a seeded θ grid spanning the three
# colocalization regimes (:coloc / :random / :exclusion) that reuses the FROZEN
# Phase-2 chain (sample_prior -> simulate_pair -> build_mci) to emit a reproducible
# `Vector{MultiChannelImage}`, each carrying its simulator GROUND-TRUTH regime
# label (D-02). The ground-truth label is the anchor that lets the comparison
# table show *when* a classical estimator is wrong.
#
# The SAME entry point also accepts a caller-supplied `Vector{MultiChannelImage}`
# (D-03) so a future external corpus (Phase 8) can be fed WITHOUT redesign and
# WITHOUT taking a Phase-8 dependency now. External images carry NO trusted ground
# truth, so their regime is labelled `:unknown` (T-09-07: never fabricate labels).
#
# REPRODUCIBILITY (D-12 / T-09-08): every generated input is a pure function of
# `(master_seed, global_index)` via the Philox `sample_rng` stream, exactly like
# the Phase-3 cache. Two builds under the same master seed are byte-identical in
# both images and labels, independent of thread count or scan order.
#
# ANTI-DATA-SNOOPING (D-14): the regime band boundaries are FIXED `const`s here,
# declared once and never tuned to make a table look good.
#
# DECOUPLING (hard constraint, CLAUDE.md): reaches src/ ONLY transitively through
# contract.jl's read-only include(). All code lives under spike/. The includes
# below are GUARDED (isdefined) so this file loads both standalone and after a
# sibling already pulled the frozen chain into scope.

# --- ORDER MATTERS: contract.jl FIRST (brings build_mci + the read-only src/
#     coupling into scope), then the simulator, then the seeding util. Mirrors
#     spike/data/generate.jl:44-47 exactly so the frozen chain is reached read-only.
isdefined(@__MODULE__, :build_mci)     || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :sample_prior)  || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair) || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :HOLDOUT_SALT)  || include(joinpath(@__DIR__, "..", "data", "seeding.jl"))
# Sibling config.jl for MASTER_SEED (D-12); guarded so inputs.jl stays loadable
# standalone as well as after run_comparator.jl already loaded the config.
isdefined(@__MODULE__, :MASTER_SEED)   || include(joinpath(@__DIR__, "config.jl"))

# --- PRE-DECLARED regime band boundaries (D-14; anti-data-snooping) --------------
# The ground-truth regime is derived from the simulator's ρ_true knob (which spans
# [-0.99, 0.99] via GHAT_RHO_KNOTS). Bands are FIXED here, BEFORE any comparison
# table exists, and are NOT tuned post-hoc:
#   ρ_true ≥ REGIME_COLOC_BAND      ⇒ :coloc      (strong positive correlation)
#   ρ_true ≤ REGIME_EXCLUSION_BAND  ⇒ :exclusion  (anti-correlation / mutual exclusion)
#   otherwise                        ⇒ :random     (near-zero correlation)
const REGIME_COLOC_BAND     = 0.30
const REGIME_EXCLUSION_BAND = -0.30

# Safety cap on the reject-and-redraw scan: at most SCAN_CAP_FACTOR * (3*n_per_regime)
# global indices are inspected before we declare a regime unfillable. Generous
# because the Cauchy μ-prior places ≥20% mass in each band, so buckets fill fast.
const REGIME_SCAN_CAP_FACTOR = 1000

"""
    regime_of(ρ_true::Real) -> Symbol

Map a simulator ground-truth `ρ_true` to its regime label using the pre-declared
(D-14) fixed band boundaries: `:coloc` above `REGIME_COLOC_BAND`, `:exclusion`
below `REGIME_EXCLUSION_BAND`, `:random` in between.
"""
function regime_of(ρ_true::Real)
    ρ_true ≥ REGIME_COLOC_BAND     && return :coloc
    ρ_true ≤ REGIME_EXCLUSION_BAND && return :exclusion
    return :random
end

"""
    build_shared_inputs(; master_seed=MASTER_SEED, n_per_regime=8,
                        imsize=(128,128)) -> NamedTuple

Build the comparator's seeded, regime-labelled shared-input set (D-02). Returns

    (images::Vector{MultiChannelImage},   # n_per_regime per regime, regime-ordered
     regime::Vector{Symbol},              # :coloc / :random / :exclusion
     rho_true::Vector{Float64},           # simulator ground-truth ρ_true per input
     master_seed)                         # the seed used (for the artifact header)

grouped in the fixed regime order `(:coloc, :random, :exclusion)`, with
`n_per_regime` inputs per regime.

Each input is a pure function of `(master_seed, global_index)` via the frozen
chain `sample_rng -> sample_prior -> simulate_pair -> build_mci` (exactly the
`generate_sample` composition, contract.jl reached READ-ONLY). To span the three
regimes despite the Cauchy μ-prior concentrating near ρ_true≈0, global indices are
scanned in increasing order and each draw is REJECT-AND-REDRAW routed into its
target regime bucket until every bucket holds `n_per_regime` inputs — still fully
seeded and order-independent (D-12), because acceptance depends only on the draw's
own `(master_seed, idx)`. Two calls with the same `master_seed` therefore produce
BYTE-IDENTICAL images and labels (T-09-08).

`imsize` dims must be ≥ 64 (the 8×8 patch grid floor enforced by `simulate_pair`).
"""
function build_shared_inputs(; master_seed::Integer = MASTER_SEED,
                             n_per_regime::Int = 8,
                             imsize::Tuple{Int,Int} = (128, 128))
    n_per_regime ≥ 1 ||
        throw(ArgumentError("n_per_regime must be ≥ 1, got $n_per_regime"))
    (imsize[1] ≥ 64 && imsize[2] ≥ 64) ||
        throw(ArgumentError("imsize dims must be ≥ 64 for the 8×8 patch grid, got $imsize"))

    targets = (:coloc, :random, :exclusion)
    # Per-regime buckets of accepted global indices (θ is a pure fn of the index,
    # regenerated at assembly, so we only need to remember which indices to keep).
    buckets   = Dict{Symbol,Vector{Int}}(r => Int[] for r in targets)
    remaining = length(targets) * n_per_regime
    max_scan  = REGIME_SCAN_CAP_FACTOR * remaining

    idx = 0
    while remaining > 0
        idx += 1
        idx > max_scan && throw(ErrorException(
            "regime scan exceeded $max_scan indices before filling all regimes " *
            "(filled: $(Dict(r => length(v) for (r, v) in buckets)))"))
        # Fresh keyed stream for regime routing only; discarded after (does not
        # perturb the assembly stream, which is reconstructed from the same key).
        θ = sample_prior(sample_rng(master_seed, idx))
        r = regime_of(θ.ρ_true)
        if length(buckets[r]) < n_per_regime
            push!(buckets[r], idx)
            remaining -= 1
        end
    end

    # Assemble in the fixed regime order. Reconstruct the keyed rng FRESH per index
    # so the sample_prior -> simulate_pair continuation matches generate_sample
    # byte-for-byte (the routing draw above used a separate rng object).
    images   = MultiChannelImage[]
    regime   = Symbol[]
    rho_true = Float64[]
    for r in targets
        for i in buckets[r]
            rng = sample_rng(master_seed, i)
            θ   = sample_prior(rng)
            mci = build_mci(simulate_pair(rng, θ; imsize = imsize); name = "sim-$(r)-$(i)")
            push!(images, mci)
            push!(regime, r)
            push!(rho_true, θ.ρ_true)
        end
    end

    return (images = images, regime = regime, rho_true = rho_true,
            master_seed = master_seed)
end

"""
    build_shared_inputs(images::Vector{<:MultiChannelImage}) -> NamedTuple

D-03 passthrough: accept a caller-supplied vector of `MultiChannelImage`s (a
future Phase-8 external corpus) and return it UNCHANGED, taking NO Phase-8
dependency. External inputs carry no trusted ground truth, so every regime label
is `:unknown` and `rho_true` is `missing` (T-09-07: never fabricate ground-truth
labels for untrusted inputs). Same NamedTuple shape as the keyword method.
"""
function build_shared_inputs(images::Vector{<:MultiChannelImage})
    n = length(images)
    return (images      = images,
            regime      = fill(:unknown, n),
            rho_true    = fill(missing, n),
            master_seed = missing)
end
