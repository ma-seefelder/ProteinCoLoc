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

# spike/comparator/classical.jl --- CMP-01/CMP-02: classical estimator battery (D-04/D-06)
#
# The columns of the Phase-9 cross-method comparison table: the classical
# fluorescence-colocalization estimators (Manders M1/M2, whole-image Pearson,
# whole-image Spearman, patch-grid correlation) plus the one genuinely NEW
# algorithm — the Costes randomization block-scramble p-value (D-05, appended
# by Task 2). Every estimator is a NaN-safe callable over a `MultiChannelImage`
# that NEVER throws on degenerate (all-zero / all-missing) input (D-06).
#
# DECOUPLING (hard constraint, CLAUDE.md):
#   * `spike/data/encode.jl` is NOT edited here. It is entry 5 of
#     `hashguard.jl:HASH_SRC_FILES`; any byte change flips `cache_hash` and
#     auto-invalidates the completed Phase-3 ≥50k-sample cache + the trained
#     Phase-4 NPE. `manders`/`pearson_whole` below therefore REPRODUCE the
#     `encode_aug` formula byte-for-byte (encode.jl:102-109) as a new callable;
#     D-04's "stay consistent" intent is enforced by an equality test (09-06),
#     NOT by refactoring the hash-guarded file.
#   * The frozen `src/` math (`patch`, `correlation`, `MultiChannelImage`) is
#     reached ONLY through the read-only `include()` boundary in `spike/contract.jl`
#     (already pulled into scope by the harness). This file adds NO second
#     `include(...src/...)` — that would double-define the frozen functions.

using Statistics     # cor, mean (whole-image Pearson, patch-grid reduction)
using StatsBase      # corspearman (whole-image Spearman)
using Random123      # Philox4x — disjoint COSTES_SALT scramble stream (mirror seeding.jl)
using Random         # randperm(rng, n) — reproducible block-position permutation

# Empty-collection-safe reduction: mirror the `_safe` NaN-not-throw idiom from
# encode.jl:80 (itself mirroring induced_mu, contract.jl:104). A reduction over
# an empty vector returns NaN rather than throwing, so every estimator stays
# total on degenerate input (the caller/table filters on isfinite).
_safe(f, v::AbstractVector) = isempty(v) ? NaN : f(v)

"""
    manders(mci) -> (M1, M2)

The Manders colocalization coefficients, computed with the SAME Otsu thresholds
(`mci.otsu_threshold[1]/[2]`) `encode_aug` uses — reproducing the encode.jl:102-106
formula byte-for-byte so the comparator column and the D-02 augmented encoder can
NEVER silently diverge (D-04). `M1` is the fraction of channel-1 signal where
channel-2 is above its threshold; `M2` the symmetric fraction for channel-2.

Degenerate guard (D-06): a channel whose total signal is 0 yields `NaN` for that
coefficient (never a divide-by-zero throw).
"""
function manders(mci)
    ch1 = mci.data[1]
    ch2 = mci.data[2]
    t1  = mci.otsu_threshold[1]
    t2  = mci.otsu_threshold[2]
    s1  = sum(ch1)
    s2  = sum(ch2)
    M1 = s1 > 0 ? sum(ch1 .* (ch2 .> t2)) / s1 : NaN   # frac. ch1 signal where ch2 present
    M2 = s2 > 0 ? sum(ch2 .* (ch1 .> t1)) / s2 : NaN   # frac. ch2 signal where ch1 present
    return (M1, M2)
end

"""
    pearson_whole(mci) -> Float64

Whole-image Pearson correlation over all pixels — identical to `encode_aug`'s
`pearson_whole` moment (encode.jl:109), so this comparator column and the
augmented-summary moment are the SAME statistic (D-04).
"""
pearson_whole(mci) = cor(vec(mci.data[1]), vec(mci.data[2]))

"""
    spearman_whole(mci) -> Float64

Whole-image Spearman rank correlation over all pixels (the rank-based classical
companion to `pearson_whole`). Uses `StatsBase.corspearman`, the same rank
correlation the frozen `src/colocalization.jl:correlation` dispatches on for
`method = :spearman`.
"""
spearman_whole(mci) = corspearman(vec(mci.data[1]), vec(mci.data[2]))

"""
    patch_correlation(mci; method::Symbol = :spearman) -> Float64

The mean per-patch correlation over the FROZEN 8×8 patch grid: calls the
UNCHANGED `patch(...,8)` and `correlation(...; method)` from
`src/colocalization.jl` reached through `contract.jl` (no second `include` of
src/), then reduces the non-missing entries via `_safe(mean, ...)`. Per patch,
`_exclude_zero` drops 0.0/NaN/missing and any patch with ≤15 survivors becomes
`missing` (the frozen contract traps); if ALL patches are missing the reduction
returns `NaN` (D-06).
"""
function patch_correlation(mci; method::Symbol = :spearman)
    x = mci.data[1]
    y = mci.data[2]
    xp, yp = patch.([x, y], 8)                    # src/colocalization.jl:37 (frozen, via contract.jl)
    ρ = correlation(xp, yp; method = method)      # src/colocalization.jl:221 (frozen, via contract.jl)
    return _safe(mean, collect(skipmissing(ρ)))
end

# --- Costes randomization significance p-value (NEW, D-05 / CMP-02) -------------
# Costes et al. 2004, Biophys J 86:3993-4003. The canonical colocalization
# significance null: scramble one channel in BLOCKS the size of the resolution
# element (NOT individual pixels), recompute the correlation, repeat n times, and
# report the fraction of scrambles whose correlation ≥ the observed one. BLOCK
# scrambling is essential — pixel scrambling destroys spatial autocorrelation and
# yields falsely tiny (over-significant) p-values (RESEARCH Pitfall 2). Statistic
# is whole-image Pearson, matching `pearson_whole` so the observed and null
# statistics are identical in kind.

"""
    costes_rng(master_seed, idx) -> Philox4x

The Costes scramble RNG, keyed by `(master_seed ⊻ COSTES_SALT, idx)` — a Philox
stream PROVABLY DISJOINT from `sample_rng`/`holdout_rng`/`fold_rng` (seeding.jl),
so the Costes null draws never collide with the data/CV streams (D-10). `COSTES_SALT`
is declared in comparator/config.jl (not redeclared here). Two constructions with
the same args produce identical draw sequences ⇒ `costes_p` is bit-reproducible.
"""
costes_rng(master_seed, idx) =
    Philox4x(UInt64, (UInt64(master_seed) ⊻ COSTES_SALT, UInt64(idx)))

# Partition `ch2` into a grid of square blocks of side `block` and randomly permute
# the block POSITIONS on that grid (a grid-cell position permutation), preserving
# each block's within-block spatial autocorrelation — the canonical Costes null.
# This is a BLOCK scramble, NOT a pixel scramble. Edge blocks are partial; the copy
# is clamped to the min overlap of source/destination cells (image bounds), and the
# output is seeded from a copy of ch2 so no pixel is ever left unassigned.
function _block_permute(rng, ch2, block)
    W, H = size(ch2)
    nbx  = cld(W, block)
    nby  = cld(H, block)
    ncell = nbx * nby
    rowrange(bi) = ((bi - 1) * block + 1):min(bi * block, W)
    colrange(bj) = ((bj - 1) * block + 1):min(bj * block, H)
    perm = randperm(rng, ncell)                    # reproducible block-position permutation
    scr  = copy(ch2)
    @inbounds for dst in 1:ncell
        src = perm[dst]
        dbi, dbj = fldmod1(dst, nby)               # dst linear -> (block-row, block-col)
        sbi, sbj = fldmod1(src, nby)               # src linear -> (block-row, block-col)
        dr = rowrange(dbi); dc = colrange(dbj)
        sr = rowrange(sbi); sc = colrange(sbj)
        nr = min(length(dr), length(sr))           # clamp partial edge blocks to bounds
        nc = min(length(dc), length(sc))
        scr[dr[1:nr], dc[1:nc]] = ch2[sr[1:nr], sc[1:nc]]
    end
    return scr
end

"""
    costes_p(master_seed, idx, mci; n=COSTES_N_SCRAMBLE, block=COSTES_BLOCK_PX) -> Float64

The Costes randomization significance p-value for `mci` (Costes et al. 2004). Draws
`n` BLOCK-scramble nulls of channel 2 (block side `block` ≈ one resolution element —
a grid-cell position permutation, NOT a pixel scramble), recomputes whole-image
Pearson on each, and returns

    p = (#{ r_scramble ≥ r_obs } + 1) / (n + 1)

with the Davison–Hinkley `+1/+1` Monte-Carlo correction so `p` is never 0 and always
lies in `[1/(n+1), 1]`. Seeded via `costes_rng(master_seed, idx)` (disjoint Philox
stream) ⇒ bit-reproducible. Returns `NaN` when the observed correlation is non-finite
(degenerate input, D-06). Non-finite scramble correlations are treated as not-exceeding.
"""
function costes_p(master_seed, idx, mci; n = COSTES_N_SCRAMBLE, block = COSTES_BLOCK_PX)
    ch1 = mci.data[1]
    ch2 = mci.data[2]
    r_obs = cor(vec(ch1), vec(ch2))
    isfinite(r_obs) || return NaN
    rng = costes_rng(master_seed, idx)
    c1 = vec(ch1)
    ge = 0
    for _ in 1:n
        scr = _block_permute(rng, ch2, block)      # block-position permutation (not pixel)
        r   = cor(c1, vec(scr))
        (isfinite(r) && r >= r_obs) && (ge += 1)
    end
    return (ge + 1) / (n + 1)                       # Davison–Hinkley +1/+1
end
