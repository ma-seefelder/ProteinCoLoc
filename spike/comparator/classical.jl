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
