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

# spike/data/seeding.jl --- DATA-01 reproducibility primitives (D-11 / D-10 / D-03).
#
# Each sample owns a Random123 `Philox4x` RNG KEYED BY ITS GLOBAL INDEX, so the
# draw stream is a pure function of `(master_seed, global_index)` -- independent
# of thread count and execution order (D-11). `Philox4x` is an `AbstractRNG`, so
# it drops straight into `sample_prior(rng)` (prior.jl:71) and
# `simulate_pair(rng, θ; imsize)` (forward.jl:110) with ZERO simulator change.
#
# The reserved ADVI holdout (D-10) is drawn from a PROVABLY DISJOINT key namespace
# (first key word XOR-salted) so it can never collide with a main-pool index, and
# the fold permutation (Wave-4 loader) uses a third disjoint salt.
#
# `sample_imsize` (D-03) draws the per-sample image size from a fixed cost-aware
# categorical so size-robustness is trained in without blowing the CPU budget.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local utility; touches no src/.

using Random123              # Philox4x (counter/key-based AbstractRNG)
using Random                 # AbstractRNG

# --- Disjoint key-namespace salts (D-10) ----------------------------------------
# Fixed nonzero UInt64s XOR-ed into the first key word to carve provably disjoint
# draw streams. Values are the SplitMix64 / golden-ratio mixing constants (good
# avalanche, well away from 0); only nonzero-and-distinct matters here.
const HOLDOUT_SALT = 0x9E3779B97F4A7C15   # reserved ADVI holdout stream (D-10)
const FOLD_SALT    = 0xD1B54A32D192ED03   # loader k-fold permutation stream (Wave-4)

# --- imsize categorical (D-03 / D-09) -------------------------------------------
# Discrete size set: all dims ≥ 64 (forward.jl validation) and 8-divisible-friendly.
# 1376×1028 is the real-data anchor (Phase-2 D-08). cost ∝ W·H, so the relative
# per-sample cost vs 256² is {1×, 4×, 16×, ~21.6×, 64×}.
const IMSIZE_SET = (
    (256, 256),     # 1×    -- the cheap workhorse
    (512, 512),     # 4×
    (1024, 1024),   # 16×
    (1376, 1028),   # ~21.6× -- real-data anchor (D-08)
    (2048, 2048),   # 64×   -- the expensive tail
)
# Cost-aware weights: heavily favor the small sizes; the ≥1024² fraction is capped
# at 0.10 so the 50k-default budget (D-06) stays CPU-tractable AND the per-thread
# peak memory (a 2048² in-flight pair ≈ 170 MB) is hit only ~2% of the time (D-12).
const IMSIZE_WEIGHTS = (0.55, 0.35, 0.05, 0.03, 0.02)   # Σ = 1.0; ≥1024² fraction = 0.10
# Expected per-sample cost (Σ wᵢ·costᵢ, units of the 256² cost), for budget sizing:
#   0.55·1 + 0.35·4 + 0.05·16 + 0.03·21.6 + 0.02·64 ≈ 4.68× the 256² baseline.
const IMSIZE_EXPECTED_COST = 4.678

"""
    sample_rng(master_seed::Integer, idx::Integer) -> Philox4x

The main-pool per-sample keyed RNG (D-11). The draw stream depends ONLY on
`(master_seed, idx)` -- two constructions with the same args produce identical
`rand` sequences regardless of which thread or in what order they run. Returns an
`AbstractRNG` usable directly in `sample_prior` / `simulate_pair`.
"""
sample_rng(master_seed::Integer, idx::Integer) =
    Philox4x(UInt64, (UInt64(master_seed), UInt64(idx)))

"""
    holdout_rng(master_seed::Integer, idx::Integer) -> Philox4x

The reserved-holdout per-sample keyed RNG (D-10). Identical construction to
`sample_rng` but with `HOLDOUT_SALT` XOR-ed into the first key word, so its draw
stream is PROVABLY DISJOINT from `sample_rng(master_seed, idx)` for the same idx
-- the reserved set can never collide with a main-pool index.
"""
holdout_rng(master_seed::Integer, idx::Integer) =
    Philox4x(UInt64, (UInt64(master_seed) ⊻ HOLDOUT_SALT, UInt64(idx)))

"""
    fold_rng(master_seed::Integer) -> Philox4x

The Wave-4 loader's k-fold permutation RNG, keyed by `(master_seed ⊻ FOLD_SALT, 0)`
so the fold split is reproducible AND independent of the generation key stream.
"""
fold_rng(master_seed::Integer) =
    Philox4x(UInt64, (UInt64(master_seed) ⊻ FOLD_SALT, UInt64(0)))

"""
    sample_imsize(rng::AbstractRNG) -> Tuple{Int,Int}

Draw one image size from `IMSIZE_SET` under the cost-aware `IMSIZE_WEIGHTS` (D-03).
Every returned dim is ≥ 64 (so the 8×8 patch grid clears the ≥15-px floor) and the
draw is reproducible under a fixed `rng`. Hand-rolled inverse-CDF over the fixed
5-element categorical (no extra dependency).
"""
function sample_imsize(rng::AbstractRNG)
    r = rand(rng)
    c = 0.0
    @inbounds for i in 1:length(IMSIZE_SET)
        c += IMSIZE_WEIGHTS[i]
        r <= c && return IMSIZE_SET[i]
    end
    return IMSIZE_SET[end]   # float-rounding guard: r just over Σw -> the last bin
end
