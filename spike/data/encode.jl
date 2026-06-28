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

# spike/data/encode.jl --- DATA-01 fixed-dimension feature encoding (D-01 / D-02).
#
# Pure transforms of the FROZEN 8x8 patch-correlation matrix `patch_summary(mci)`
# (contract.jl:85) into the training-vector layouts the amortized nets consume:
#
#   encode_d01(M) -> 128-dim  : the LOCKED D-01 encoding. Rows 1:64 are the
#       imputed (missing->0) per-patch correlations in column-major vec() order;
#       rows 65:128 are the parallel binary present/absent mask (1.0 present,
#       0.0 missing). The mask makes degeneracy self-describing, so Phase-4's
#       summary-net input dim is unambiguously 128.
#   encode_aug(mci, M) -> AUG_DIM : the D-02 augmented SUPERSET = the 128 D-01
#       dims PLUS N_AUG_MOMENTS scalar moments computed over the two channel
#       matrices + the 8x8 matrix, so Phase-4 (ABL-01) can slice extra features
#       with zero re-simulation.
#
# DEGENERACY (D-13, mirrors induced_mu's NaN-not-throw choice, contract.jl:104):
# a fully-missing 8x8 is KEPT (vals all 0, mask all 0) -- NEVER dropped. Dropping
# degenerate samples would distort the pi(theta)-faithful training distribution.
#
# DECOUPLING (hard constraint, CLAUDE.md): pure spike-local helpers over the
# matrix `patch_summary` already returns; never re-computes or edits the src/
# summary. Reached only transitively through contract.jl's read-only include().

using StatsBase      # skewness, kurtosis (the patch-grid tail moments)
using Statistics     # mean, median, std, quantile, cor

# --- Augmented-moment layout contract (D-02) ------------------------------------
# A NAMED constant so the cache layout is a contract, not a magic number. The 14
# moments (in append order) are: Manders M1, M2; whole-image Pearson; patch-grid
# median, IQR, mean (=induced_mu), std, skewness, excess kurtosis; fraction
# missing; per-channel intensity median + IQR (ch1, ch2). Exact composition is
# Claude's discretion (D-02); only "cache D-01 + an augmented superset" is locked.
const N_AUG_MOMENTS = 14
const AUG_DIM       = 128 + N_AUG_MOMENTS   # = 142

"""
    encode_d01(M::AbstractMatrix) -> Vector{Float64}

The LOCKED 128-dim D-01 encoding of the 8x8 `Matrix{Union{Float64,Missing}}`
returned by `patch_summary` (contract.jl:85).

  vals = vec(coalesce.(M, 0.0))          # 64-dim, missing -> 0, column-major
  mask = vec(Float64.(.!ismissing.(M)))  # 64-dim, 1.0 present / 0.0 missing
  return vcat(vals, mask)                 # 128-dim; rows 1:64 vals, 65:128 mask

Rows 1:64 are the imputed correlations in column-major `vec()` order of the 8x8
grid; rows 65:128 are the parallel binary mask in the same ordering. A fully-
missing 8x8 returns a 128-vector of zeros (sample KEPT, not dropped, D-13) --
the mask channel encodes the degeneracy for the net.
"""
function encode_d01(M::AbstractMatrix)
    vals = vec(coalesce.(M, 0.0))                 # 64-dim, missing->0 (column-major)
    mask = vec(Float64.(.!ismissing.(M)))         # 64-dim, 1=present / 0=missing
    return vcat(vals, mask)                        # 128-dim; rows 1:64 vals, 65:128 mask
end

# Empty-collection-safe reduction: mirror induced_mu's NaN-not-throw degeneracy
# idiom (contract.jl:104). A reduction over a fully-missing patch grid returns NaN
# rather than throwing, so encode_aug stays total (the mask already flags it).
_safe(f, v::AbstractVector) = isempty(v) ? NaN : f(v)

"""
    encode_aug(mci::MultiChannelImage, M::AbstractMatrix) -> Vector{Float64}

The D-02 augmented SUPERSET: `vcat(encode_d01(M), moments)` of length `AUG_DIM`,
whose first 128 rows are exactly `encode_d01(M)`. `moments` (length
`N_AUG_MOMENTS`) are per-pair scalar descriptors over the two channel matrices
and the 8x8 matrix -- Manders M1/M2 (thresholds = `mci.otsu_threshold`, already
on the struct, contract.jl:68), whole-image Pearson, the patch-grid distribution
moments (median/IQR/mean/std/skewness/excess-kurtosis), fraction-missing, and
per-channel intensity median+IQR. Each `skipmissing`/empty reduction is guarded
the way `induced_mu` is, so a degenerate (fully-missing) input yields documented
NaN-safe moments rather than an error. On a non-degenerate fixture all entries
are finite.
"""
function encode_aug(mci::MultiChannelImage, M::AbstractMatrix)
    ch1 = mci.data[1]
    ch2 = mci.data[2]
    t1  = mci.otsu_threshold[1]
    t2  = mci.otsu_threshold[2]

    # --- Manders colocalization coefficients (thresholds from the MCI) ----------
    s1 = sum(ch1)
    s2 = sum(ch2)
    M1 = s1 > 0 ? sum(ch1 .* (ch2 .> t2)) / s1 : NaN   # frac. ch1 signal where ch2 present
    M2 = s2 > 0 ? sum(ch2 .* (ch1 .> t1)) / s2 : NaN   # frac. ch2 signal where ch1 present

    # --- whole-image Pearson over all pixels ------------------------------------
    pearson_whole = cor(vec(ch1), vec(ch2))

    # --- patch-grid distribution moments over the present (non-missing) entries -
    v = collect(skipmissing(M))
    grid_median   = _safe(median, v)
    grid_iqr      = _safe(x -> quantile(x, 0.75) - quantile(x, 0.25), v)
    grid_mean     = _safe(mean, v)                       # = induced_mu
    grid_std      = _safe(std, v)
    grid_skew     = _safe(StatsBase.skewness, v)
    grid_kurt     = _safe(StatsBase.kurtosis, v)         # excess kurtosis
    frac_missing  = count(ismissing, M) / length(M)      # /64

    # --- per-channel intensity descriptors --------------------------------------
    c1 = vec(ch1); c2 = vec(ch2)
    ch1_median = median(c1)
    ch1_iqr    = quantile(c1, 0.75) - quantile(c1, 0.25)
    ch2_median = median(c2)
    ch2_iqr    = quantile(c2, 0.75) - quantile(c2, 0.25)

    moments = Float64[
        M1, M2,
        pearson_whole,
        grid_median, grid_iqr, grid_mean, grid_std, grid_skew, grid_kurt,
        frac_missing,
        ch1_median, ch1_iqr, ch2_median, ch2_iqr,
    ]
    @assert length(moments) == N_AUG_MOMENTS "moment vector length $(length(moments)) != N_AUG_MOMENTS"
    return vcat(encode_d01(M), moments)
end
