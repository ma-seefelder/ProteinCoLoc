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

# spike/validation/harness.jl --- Phase-5 shared harness (SBC-01, BF, OOD).
#
# THE ONE θ*~π → simulate → infer chain all three Phase-5 deliverables ride. It is
# pure COMPOSITION of the frozen amortized read surface (infer.jl) over the FROZEN
# Phase-4 net -- it adds no inference machinery, only the outer draw loop and the
# paired-Δρ path (D-01).
#
# CENTRAL INVARIANTS (every consumer inherits them):
#   • FROZEN STATS (T-05-01): the net + BOTH transforms (zt, θzt) are loaded ONCE
#     via load_npe; standardization is NEVER re-fit here (no fit(ZScoreTransform)).
#     Every summary is put into the estimator's input space with the frozen m.zt,
#     and every ρ read un-standardizes with m.θzt (Pitfall 5).
#   • CPU-ONLY (D-10 / Pitfall 1): use_gpu = false on EVERY NeuralEstimators call.
#   • DISJOINT SEED (D-02): val_rng threads a Random123 stream keyed off
#     VAL_MASTER_SEED ⊻ VAL_SALT -- provably disjoint from the training/holdout/fold
#     streams (seeding.jl) AND from NPE_MASTER_SEED (VAL_MASTER_SEED ≠ 0xC0FFEE), so
#     the reported numbers never reuse the stream the net was trained on.
#   • READ-ONLY src/ (hard constraint, CLAUDE.md): src/ is reached only transitively
#     through contract.jl's read-only include(); this file lives entirely under spike/.
#
# Flat top-level functions (sibling style of infer.jl -- no module wrapper). Guarded,
# order-dependent includes keep it loadable standalone AND idempotent under runtests.jl.

using StatsBase          # reconstruct (inverse ZScoreTransform)
using Statistics         # mean
using Random123          # Philox4x (counter/key-based AbstractRNG for the VAL stream)

# --- ORDER MATTERS: consts first, then the trainer/inference surface (load_npe +
#     posterior_for + rho_draws), then the simulator halves, contract, encode, and
#     the seeding primitives whose salt idiom val_rng mirrors. Guarded for idempotency.
# SENTINEL IS :SBC_M, NOT :VAL_MASTER_SEED. `consts.jl:38-42` guards its own block on :SBC_M and
# says so; but :VAL_MASTER_SEED is ALSO declared by p11_consts.jl:87, p12_consts.jl:103 and
# p13/consts.jl:109, so guarding on that name made this include a silent no-op whenever any Phase
# 11/12/13 pre-registration loaded first -- and then SBC_M / SBC_FIX_M never existed at all.
# :SBC_M is owned by consts.jl alone, so the caller guard and the file's internal guard agree.
isdefined(@__MODULE__, :SBC_M) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :load_npe)        || include(joinpath(@__DIR__, "..", "npe", "train_npe.jl"))
isdefined(@__MODULE__, :posterior_for)   || include(joinpath(@__DIR__, "..", "npe", "infer.jl"))
isdefined(@__MODULE__, :sample_prior)    || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair)   || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :build_mci)       || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :encode_d01)      || include(joinpath(@__DIR__, "..", "data", "encode.jl"))
isdefined(@__MODULE__, :sample_rng)      || include(joinpath(@__DIR__, "..", "data", "seeding.jl"))

# --- Disjoint VAL key-namespace salt (D-02) -------------------------------------
# A fixed nonzero UInt64 XOR-ed into the first key word to carve the reported
# validation draw stream PROVABLY DISJOINT from the generation salts in seeding.jl
# (HOLDOUT_SALT, FOLD_SALT) and from the un-salted main-pool stream. Value is a
# SplitMix64 mixing constant (good avalanche, distinct from the other two salts);
# only nonzero-and-distinct matters.
if !isdefined(@__MODULE__, :VAL_SALT)
    const VAL_SALT = 0xBF58476D1CE4E5B9   # distinct from HOLDOUT_SALT / FOLD_SALT
end

"""
    load_frozen_model(path = joinpath(@__DIR__, "..", "npe", "trained_npe.jld2")) -> NamedTuple

Load the FROZEN Phase-4 NPE ONCE: `(estimator, θzt, zt, variant, d_in, meta)` via
`load_npe` (which integrity-checks the "estimator" key on open). This is the single
frozen-stats source the whole Phase-5 bundle infers against -- the net is NEVER
re-fit and the transforms `zt`/`θzt` are the deployed ones (T-05-01).
"""
function load_frozen_model(path = joinpath(@__DIR__, "..", "npe", "trained_npe.jld2"))
    return load_npe(path)
end

"""
    val_rng(master_seed = VAL_MASTER_SEED) -> Philox4x

The Phase-5 validation RNG: a Random123 `Philox4x` counter-based `AbstractRNG` keyed
by `(master_seed ⊻ VAL_SALT, 0)`. Its draw stream is PROVABLY DISJOINT from the
training/holdout/fold streams (seeding.jl uses different salts) and, at the default
`VAL_MASTER_SEED = 0x5BC0FFEE ≠ NPE_MASTER_SEED = 0xC0FFEE`, from the stream the net
was trained on (D-02 fresh disjoint stream). Threaded sequentially through the SBC/
BF/OOD draw loops; two constructions with the same `master_seed` are bit-identical.
Pass `VAL_FIX_SEED` for the fast fixture gates so they never consume the reported stream.
"""
val_rng(master_seed = VAL_MASTER_SEED) =
    Philox4x(UInt64, (UInt64(master_seed) ⊻ VAL_SALT, UInt64(0)))

"""
    draw_simulate_infer(m, rng; imsize = SBC_IMSIZE, N = SBC_L) -> NamedTuple

ONE θ*~π → simulate → infer step (RESEARCH Pattern 1). Draws a prior θ*, simulates a
channel pair, runs it through the UNCHANGED src/ summary contract, standardizes with
the FROZEN `m.zt`, takes ONE `sampleposterior` pass, and un-standardizes with the
frozen `m.θzt`. Returns `(θ = θ*, Z = standardized-summary, draws = 7×N physical-θ)`.
`draws` is in PHYSICAL θ units (row 1 = ρ_true); no standardization is re-fit here.
CPU-only (`use_gpu = false`).
"""
function draw_simulate_infer(m, rng; imsize = SBC_IMSIZE, N = SBC_L)
    θ   = sample_prior(rng)                                          # θ* ~ π(θ)
    mci = build_mci(simulate_pair(rng, θ; imsize = imsize))          # forward sim → MCI
    Z   = standardize_summary(encode_d01(patch_summary(mci)), m.zt, :min)  # FROZEN m.zt
    draws_std = posterior_for(m.estimator, Z; N = N, use_gpu = false)      # 7×N standardized
    draws     = StatsBase.reconstruct(m.θzt, draws_std)             # 7×N physical (Pitfall 5)
    return (θ = θ, Z = Z, draws = draws)
end

"""
    draw_simulate_infer_paired(m, rng; imsize = SBC_IMSIZE, N = SBC_L) -> NamedTuple

The paired-draw path for the Δρ SBC (D-01): draw TWO INDEPENDENT priors θ*_s, θ*_c,
simulate/encode/standardize both against the FROZEN `m.zt`, and return the N
un-standardized ρ_true draw VECTORS for each via `rho_draws` (reused directly, NOT
re-implemented). Returns `(θs, θc, ρs, ρc)`; consumers form the MC difference
`Δρ_draws = ρs .- ρc` and the truth `Δρ* = θs.ρ_true - θc.ρ_true`. CPU-only.
"""
function draw_simulate_infer_paired(m, rng; imsize = SBC_IMSIZE, N = SBC_L)
    θs   = sample_prior(rng)
    mcis = build_mci(simulate_pair(rng, θs; imsize = imsize))
    Zs   = standardize_summary(encode_d01(patch_summary(mcis)), m.zt, :min)

    θc   = sample_prior(rng)
    mcic = build_mci(simulate_pair(rng, θc; imsize = imsize))
    Zc   = standardize_summary(encode_d01(patch_summary(mcic)), m.zt, :min)

    ρs = rho_draws(m.estimator, Zs, m.θzt; N = N, use_gpu = false)   # N ρ_true draws (sample)
    ρc = rho_draws(m.estimator, Zc, m.θzt; N = N, use_gpu = false)   # N ρ_true draws (control)
    return (θs = θs, θc = θc, ρs = ρs, ρc = ρc)
end
