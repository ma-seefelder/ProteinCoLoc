#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder

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

#############################################################################################
# src/amortized/summary.jl --- grid-parametric patch-correlation summary + encoder (PROD-02).
#
# The ONE shared coupling site between the patch grid `G` and the amortized summary vector.
# D-04 couples the summary-vector dimension to `G`: an `G×G` per-patch Pearson correlation
# matrix encodes to a length-`2·G²` vector (imputed values + a binary present-mask). This file
# is the SINGLE SOURCE OF TRUTH for that coupling — every downstream amortized module (datagen,
# training, inference, OOD, registry) derives its dimensions from the helpers here, so a shipped
# grid is added by instantiating the SAME parametrized encoder, never by re-hardcoding sizes.
#
# Promoted (grid-generalized) from the proven spike:
#   spike/contract.jl:85-90    patch_summary(mci)          -> patch_summary(mci, G)
#   spike/data/encode.jl:71-75 encode_d01(M)               -> body unchanged (vec handles any G)
#   spike/data/loader.jl:126-136 _row_partition(:min, 128) -> _summary_row_partition(:min, 2G²)
#
# The underlying `patch()`/`correlation()` in src/colocalization.jl are ALREADY grid-general
# (they read `num_patches`), so `G` threads straight through; the ≥15-survivor floor that maps a
# sparse patch to `missing` lives unchanged in `correlation`/`_exclude_zero`.
#############################################################################################

"""
    patch_summary(mci::MultiChannelImage, G::Integer) -> Matrix{Union{Float64,Missing}}

The grid-parametric SIM-03 summary statistic: the `G×G` per-patch Pearson correlation matrix of
the first two channels of `mci`, computed by the UNCHANGED `src/colocalization.jl` `patch()` and
`correlation()`. Each patch with `≤15` non-zero survivors becomes `missing` (the sparsity floor
in `correlation`/`_exclude_zero`). `G` is the registry key; `G = 8` reproduces the spike summary.

This is a plain module-local function (NOT an extension of `Base.summary`, whose contract is
`(io, x)::String` — returning a Matrix there would poison Julia's `show`/REPL machinery).
"""
function patch_summary(mci::MultiChannelImage, G::Integer)
    x = mci.data[1]
    y = mci.data[2]
    xp, yp = patch.([x, y], G)                        # src/colocalization.jl:37 — grid-general
    return correlation(xp, yp; method = :pearson)     # src/colocalization.jl:221 — grid-general
end

"""
    encode_d01(M::AbstractMatrix) -> Vector{Float64}

The `2·G²`-dim encoding of the `G×G` `Matrix{Union{Float64,Missing}}` returned by
`patch_summary`:

  vals = vec(coalesce.(M, 0.0))          # G²-dim, missing -> 0, column-major
  mask = vec(Float64.(.!ismissing.(M)))  # G²-dim, 1.0 present / 0.0 missing
  return vcat(vals, mask)                # 2·G²-dim; rows 1:G² vals, G²+1:2G² mask

Rows `1:G²` are the imputed correlations in column-major `vec()` order; rows `G²+1:2G²` are the
parallel binary present-mask in the same ordering (the mask makes degeneracy self-describing to
the summary net). A fully-missing grid returns a `2·G²`-vector of zeros (sample KEPT, not
dropped). The body is grid-agnostic: `vec` handles any `G`.
"""
function encode_d01(M::AbstractMatrix)
    vals = vec(coalesce.(M, 0.0))                 # G²-dim, missing->0 (column-major)
    mask = vec(Float64.(.!ismissing.(M)))         # G²-dim, 1=present / 0=missing
    return vcat(vals, mask)                        # 2·G²-dim; rows 1:G² vals, G²+1:2G² mask
end

"""
    _summary_row_partition(variant::Symbol, nrows::Integer) -> (cont_rows, mask_rows)

Split the `nrows` feature rows of a summary vector into the CONTINUOUS rows (z-scored at
standardization) and the binary MASK rows (passed through unchanged). For `:min` — the shipped
D-04 encoding — the split is grid-general: `cont = 1:G²`, `mask = G²+1:2G²` where `G² = nrows÷2`.

z-scoring a 0/1 mask would re-couple folds through the mask mean, so the mask rows are NEVER
standardized. The hardcoded 128-row guard of the spike loader is replaced by an `iseven(nrows)`
check, so the split is correct for every `G` (`nrows = 2·G²` is always even).
"""
function _summary_row_partition(variant::Symbol, nrows::Integer)
    if variant === :min
        iseven(nrows) || error(
            "_summary_row_partition: :min expects an even row count (2·G²), got $nrows")
        half = nrows ÷ 2
        return (collect(1:half), collect(half+1:nrows))
    else
        error("_summary_row_partition: unknown variant $variant (shipped path uses :min)")
    end
end

# --- Single-source dimension helpers (grid G -> summary/encoder/NRE dimensions) -------------
# Every amortized module derives its buffer/architecture widths from these, so the grid coupling
# is stated exactly once. `G = 8` reproduces the proven spike constants (128 / 64 / 320).

"""
    summary_dim(G::Integer) -> Int

Length of the D-01 encoded summary vector for a `G×G` patch grid: `2·G²` (imputed values +
present-mask). `summary_dim(8) == 128`.
"""
summary_dim(G::Integer) = 2 * G^2

"""
    cont_rows(G::Integer) -> Int

Number of CONTINUOUS (z-scored) summary rows for a `G×G` grid: `G²` (the imputed correlations;
the other `G²` rows are the binary mask). `cont_rows(8) == 64`.
"""
cont_rows(G::Integer) = G^2

"""
    ratio_input_dim(G::Integer) -> Int

Input dimension of the NRE pair-encoding for a `G×G` grid: `5·G²`. The A7 difference encoding
concatenates the two `2·G²` summaries plus their `G²`-dim continuous-row contrast
(`2·(2G²) + G² = 5G²`). `ratio_input_dim(8) == 320`.
"""
ratio_input_dim(G::Integer) = 5 * G^2
