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

# spike/npe/p12_architecture.jl --- the SPATIAL VIEW of an UNCHANGED summary, plus the lean
# convolutional head, the wide flow, and the 72-row theta layout every later Phase-12 plan
# indexes by name.
#
# WHAT THIS FILE DOES NOT DO: it does not change the summary statistic. The 8x8
# patch-correlation grid and its `encode_d01` encoding are byte-for-byte what Phase 4 shipped.
# D-02's claim is about the INDUCTIVE BIAS, not the information: the per-region content is
# already present in the 128 rows, and flattening them into an MLP discards the ADJACENCY, not
# the data. `reshape_summary` below is a pure re-view of the same numbers.
#
# THE COLUMN-MAJOR CONTRACT IS INHERITED, NOT CHOSEN, AND MUST NEVER BE TRANSPOSED.
# `src/amortized/summary.jl:72-76` (`encode_d01`, READ-ONLY reference) builds
# `vals = vec(coalesce.(M, 0.0))` COLUMN-MAJOR from the GxG matrix `patch_summary` returns and
# `mask = vec(Float64.(.!ismissing.(M)))` in the same ordering, then `vcat`s them -- rows 1:G^2
# are the imputed correlations, rows G^2+1:2G^2 the binary present-mask. So
# `reshape(rows[1:G^2], G, G)` recovers that matrix EXACTLY, with NO transpose. A transposed
# implementation rotates the lattice relative to the image, survives visual inspection, and
# passes every symmetric fixture -- which is why testset 1 of `test_p12_architecture.jl`
# asserts exact equality against a locally reconstructed encoding rather than eyeballing it
# (T-12-14).
#
# THE MASK ROWS ARE NEVER STANDARDIZED, AND THEY ALWAYS LAND IN CHANNEL 2.
# `src/amortized/summary.jl:85-87` gives the reason in one line: "z-scoring a 0/1 mask would
# re-couple folds through the mask mean". The frozen `zt` therefore z-scores rows 1:G^2 only and
# passes rows G^2+1:2G^2 through unchanged; `reshape_summary` preserves that split as channel 1
# (continuous) and channel 2 (binary), so the convolution sees a genuine 0/1 occupancy plane.
#
# THE ANTI-PATTERN, STATED SO IT CANNOT BE RE-INTRODUCED BY ACCIDENT: there is NO
# `GlobalMeanPool`, NO `DeepSet` and NO permutation-invariant pooling of any kind anywhere in
# this net. Permutation-invariance over patches IS the exchangeable pooling this phase exists to
# replace (12-RESEARCH.md Anti-Patterns; 12-CONTEXT records it as an SC1 correction). A
# flatten-then-Dense head after the convolutions is the entire point -- it is what makes region
# (2,3) and region (7,7) distinguishable to the flow. `test_p12_architecture.jl` testset 3
# asserts the absence on the COMMENT-STRIPPED source, so this paragraph documents the rule
# without weakening the check on the code (T-12-18).
#
# WHY THERE IS NO IMAGE-TO-IMAGE HEAD. The NeuralEstimators docs' spatially-varying-parameters
# example has exactly this phase's shape (image in, field out) but uses `PointEstimator` with a
# `Unet` -- a POINT estimator with no posterior. `Unet` is not exported by the pinned v0.2.1
# (12-RESEARCH.md Pattern 4, verified against the local depot). D-01's finite-dimensional
# vector head is therefore not merely a calibration convenience: it is the only route the
# pinned library supports.
#
# ZERO NEW DEPENDENCIES (T-12-09). `Flux` and `NeuralEstimators` are already the two pinned
# neural deps; `Random` is a stdlib reached through the default LOAD_PATH `@stdlib` entry, the
# same mechanism `spike/validation/run_p11_recovery.jl:50-53` uses for `LinearAlgebra`. Nothing
# below installs a package or mutates the load path -- `spike/test/runtests.jl:123-126` asserts
# the exact 16-name dependency set and a resolve could dislodge the NeuralEstimators 0.2.1 pin.
#
# Run (self-check):
#     julia --project=spike -e 'include("spike/npe/p12_architecture.jl"); build_p12_estimator()'

using Flux               # Chain, Conv, Dense, SamePad, gelu, flatten
using NeuralEstimators   # PosteriorEstimator, NormalisingFlow
using Random             # AbstractRNG, randperm -- MCAR masking augmentation

# ORDER MATTERS: the Tier-1 pre-registration (P12_G, P12_K_DEV, P12_D_MINISPIKE,
# P12_MASK_K_SET, P12_RESEARCH_NET_DEVIATIONS) before anything that consumes it, then the
# lattice arithmetic that owns the shared indexing contract `p12_idx`. Both guarded for
# idempotency, the `run_p11_recovery.jl:55-56` idiom.
isdefined(@__MODULE__, :P12_G)   || include(joinpath(@__DIR__, "..", "validation", "p12_consts.jl"))
isdefined(@__MODULE__, :p12_idx) || include(joinpath(@__DIR__, "..", "simulator", "p12_lattice.jl"))

# =============================================================================================
# 1. The spatial view of the unchanged summary
# =============================================================================================

"""
    reshape_summary(Zraw::AbstractMatrix, G::Integer = P12_G) -> Array{Float32,4}

Re-view the `(2G²) × K` encoded summary -- exactly what `encode_d01` produces, `K` data sets
along the LAST dimension -- as a `(G, G, 2, K)` array. Channel 1 is the correlation lattice
(rows `1:G²`), channel 2 is the binary present-mask lattice (rows `G²+1:2G²`).

`reshape(view, G, G)` recovers the original `patch_summary` matrix EXACTLY, with **no
transpose**: `encode_d01` `vec`s column-major, so summary row `(j-1)*G + i` is patch `(i, j)`
(the `p12_idx` contract). Throws `DimensionMismatch` unless `size(Zraw, 1) == 2G²`.

**ORDER INVARIANT -- A CORRECTNESS CONTRACT, NOT A STYLE PREFERENCE.** `Zraw` must ALREADY have
been through the frozen standardization path (`standardize_summary(..., zt, :min)`, which
z-scores the `G²` continuous rows and passes the `G²` mask rows through untouched). Standardize
FIRST, then reshape. A reshape-then-standardize ordering re-couples the mask: once the mask rows
live in a channel of a 4-D array there is no row partition left to exempt them, so a
whole-array z-score would standardize a 0/1 plane against its own mean -- the exact fold
re-coupling `src/amortized/summary.jl:85-87` forbids, and it fails silently rather than loudly.
This mirrors `spike/npe/p11_architecture.jl:330-333`, where the same ordering trap is asserted
at runtime for the λ row.

The vector method is the single-data-set case, `reshape_summary(v, G) == reshape_summary(reshape(v, :, 1), G)`.
"""
function reshape_summary(Zraw::AbstractMatrix, G::Integer = P12_G)
    G >= 2 || throw(ArgumentError("reshape_summary: G must be ≥ 2, got $G"))
    n = G^2
    size(Zraw, 1) == 2n || throw(DimensionMismatch(
        "reshape_summary: expected $(2n) rows (2·G² at G = $G), got $(size(Zraw, 1))"))
    K = size(Zraw, 2)
    out = Array{Float32,4}(undef, G, G, 2, K)
    @views out[:, :, 1, :] .= reshape(Zraw[1:n,        :], G, G, K)
    @views out[:, :, 2, :] .= reshape(Zraw[n+1:2n,     :], G, G, K)
    return out
end

reshape_summary(v::AbstractVector, G::Integer = P12_G) = reshape_summary(reshape(v, :, 1), G)

# =============================================================================================
# 2. Region masking: the read-time construction and the training-time augmentation
# =============================================================================================

"""
    mask_regions(Zraw, idxs; G = P12_G) -> copy of Zraw

Return a COPY of the encoded summary with the continuous row AND the mask row of every region
in `idxs` set to zero. `idxs` may be a single region index or any iterable of them, each in
`1:G²` in the `p12_idx` column-major ordering.

**This is the SAME encoding `encode_d01` already produces for a patch that falls below the
≥15-surviving-pixel floor** -- value 0, mask 0. "Region unobserved" is therefore expressed in
the vocabulary the shipped encoder already speaks, rather than in a new sentinel that the net
would have to be taught. That is what makes D-09's leave-region-out construction a read of the
trained model rather than a read of an out-of-vocabulary input.
"""
function mask_regions(Zraw::AbstractMatrix, idxs; G::Integer = P12_G)
    G >= 2 || throw(ArgumentError("mask_regions: G must be ≥ 2, got $G"))
    n = G^2
    size(Zraw, 1) == 2n || throw(DimensionMismatch(
        "mask_regions: expected $(2n) rows (2·G² at G = $G), got $(size(Zraw, 1))"))
    out = copy(Zraw)
    z = zero(eltype(out))
    for r in idxs
        (1 <= r <= n) || throw(ArgumentError(
            "mask_regions: region index $r is outside 1:$n (G = $G)"))
        @views out[r,     :] .= z
        @views out[n + r, :] .= z
    end
    return out
end

mask_regions(v::AbstractVector, idxs; G::Integer = P12_G) =
    vec(mask_regions(reshape(v, :, 1), idxs; G = G))

"""
    sample_mask_k(rng; k_set = P12_MASK_K_SET) -> Int

Draw the number of regions to mask for one training sample, uniformly from the pre-registered
`P12_MASK_K_SET` (Tier 1: `0:8`, i.e. 0 – 12.5 % of the 64 regions).
"""
sample_mask_k(rng::AbstractRNG; k_set = P12_MASK_K_SET) = rand(rng, k_set)

"""
    augment_mask!(Zraw, rng; G = P12_G, k_set = P12_MASK_K_SET) -> Zraw

MCAR training-time region-masking augmentation, IN PLACE, one independent draw of `k` per
column (per data set). For each column: draw `k` from `k_set`, choose `k` distinct regions
uniformly without replacement, and zero both the continuous and the mask row of each.

**WHY THIS EXISTS (Pitfall 4).** On the frozen F5 image-size mixture with `BG_FLOOR = 0.02`, a
patch almost never falls below the ≥15-surviving-pixel floor, so masked regions are RARE in
training. A net that has essentially never seen a masked region would turn D-09's
leave-region-out coverage measurement into a measurement of OUT-OF-DISTRIBUTION behaviour --
covariate shift -- rather than of the spatial borrowing D-09 is asking about. This augmentation
makes the held-out-region configuration in-distribution at read time. It is a TRAINING-TIME
transformation only; nothing at read time calls it.

**THE WARNING SIGN, RECORDED IN ADVANCE.** If leave-region-out coverage comes out wildly off
nominal for BOTH the spatial model AND the D-10 matched ablation, suspect this augmentation
(too little masking seen, or masked at the wrong rate) BEFORE suspecting the prior. A shared
failure across two arms that differ only in spatial correlation cannot be caused by the thing
that differs between them.
"""
function augment_mask!(Zraw::AbstractMatrix, rng::AbstractRNG;
                       G::Integer = P12_G, k_set = P12_MASK_K_SET)
    G >= 2 || throw(ArgumentError("augment_mask!: G must be ≥ 2, got $G"))
    n = G^2
    size(Zraw, 1) == 2n || throw(DimensionMismatch(
        "augment_mask!: expected $(2n) rows (2·G² at G = $G), got $(size(Zraw, 1))"))
    z = zero(eltype(Zraw))
    for c in axes(Zraw, 2)
        k = sample_mask_k(rng; k_set = k_set)
        k == 0 && continue
        perm = randperm(rng, n)
        @inbounds for t in 1:k
            r = perm[t]
            Zraw[r,     c] = z
            Zraw[n + r, c] = z
        end
    end
    return Zraw
end

# =============================================================================================
# 3. The LEAN summary net (R-7) -- written literally, with the measured cost that justifies it
# =============================================================================================
#
# THE TOPOLOGY BELOW IS PRE-REGISTERED AND THE MEASUREMENT THAT PICKED IT IS RECORDED HERE, SO
# THAT NOBODY REDISCOVERS THE EXPENSIVE ONE (R-7; 12-RESEARCH.md Compute Budget, measured in the
# live spike environment 2026-07-27, D = 72 over a 50 000-pair pool):
#
#   LEAN (this file)  70 872 summary parameters ->  1.8 h training
#   FAT (REJECTED)    a 2 -> 32 -> 64 -> 64 convolution stack with no channel bottleneck,
#                     flattening 8*8*64 = 4096 into a 4096 -> 256 Dense:
#                     1 137 760 summary parameters -> 12.5 h training
#
# That is a 7x training-cost multiplier for no representational argument at G = 8. The single
# wide Dense at the head of the fat variant carries 92 % of its summary parameters ALL BY
# ITSELF; the 1x1 convolution below is a channel bottleneck placed specifically to make that
# Dense small. A task specification that says only "a CNN" gets the 12.5-hour version, which is
# why the layer stack is transcribed literally rather than described.
#
# The layer-count and channel-width assertions live in `test_p12_architecture.jl` testset 3,
# with the rejected 1 137 760 named there too, so a later "small" widening announces itself.

"""
    build_p12_summary_net(; G = P12_G, dstar = 128) -> Chain

The LEAN pre-registered convolutional summary network (R-7). Maps a `(G, G, 2, K)` spatial view
of the encoded summary to a `dstar × K` matrix of learned summaries:

```
Conv((3,3), 2  => 16, gelu; pad = SamePad())   #  G×G×2  → G×G×16
Conv((3,3), 16 => 32, gelu; pad = SamePad())   #  G×G×16 → G×G×32
Conv((1,1), 32 =>  8, gelu)                    #  G×G×32 → G×G×8   channel bottleneck
Flux.flatten                                   #  → 8G² (= 512 at G = 8)
Dense(8G², dstar)                              #  → dstar (= 128)
```

At `G = 8, dstar = 128` this is 70 872 trainable parameters. The head width is DERIVED from
`G`, never hard-coded, so the same builder is correct at any grid size.

The `SamePad()` on the two 3×3 convolutions keeps the lattice at `G × G` throughout, so every
region retains a one-to-one correspondence with a spatial position right up to the flatten --
which is what makes the head able to tell region `(2,3)` from region `(7,7)`.
"""
function build_p12_summary_net(; G::Integer = P12_G, dstar::Integer = 128)
    G     >= 2 || throw(ArgumentError("build_p12_summary_net: G must be ≥ 2, got $G"))
    dstar >= 1 || throw(ArgumentError("build_p12_summary_net: dstar must be ≥ 1, got $dstar"))
    return Chain(
        Conv((3,3), 2  => 16, gelu; pad = SamePad()),   #  G×G×2  → G×G×16
        Conv((3,3), 16 => 32, gelu; pad = SamePad()),   #  G×G×16 → G×G×32
        Conv((1,1), 32 =>  8, gelu),                    #  G×G×32 → G×G×8   channel bottleneck
        Flux.flatten,                                   #  → 8G²  (512 at G = 8)
        Dense(8 * G^2, dstar),                          #  → dstar
    )
end

# =============================================================================================
# 4. The estimator: lean CNN + a wide NormalisingFlow over 72 marginals
# =============================================================================================

"""
    build_p12_estimator(; G = P12_G, D = P12_D_MINISPIKE, dstar = 128,
                        num_coupling_layers = 10, flow_depth = 2, flow_width = 128)
        -> PosteriorEstimator

The Phase-12 mini-spike estimator: the lean convolutional summary net over the `(G, G, 2, K)`
spatial view, feeding a NeuralEstimators `NormalisingFlow` over the `D`-dimensional θ
(`P12_D_MINISPIKE = 72`; see `P12_THETA_ROWS`).

**VERIFIED v0.2.1 FACTS THIS CONSTRUCTION RESTS ON** (12-RESEARCH.md Pattern 4, read from the
pinned local depot and run live):

  * The summary network is CALLED DIRECTLY -- `_summarystatistics(estimator, Z) =
    estimator.summary_network(Z)` (`Estimators.jl:76`). There is no `DeepSet` requirement and no
    exchangeable-pooling wrapper: ANY Flux model that accepts the 4-D `Z` works, which is what
    lets a CNN drop in where the Phase-4/5 MLP was.
  * `NormalisingFlow` is DIMENSION-GENERIC (`CouplingLayer` splits `d` into ⌊d/2⌋ / ⌈d/2⌉ and
    works for any `d ≥ 2`), so widening θ from 7 to 72 rows needs **no custom
    `ApproximateDistribution`**. The flow just gets a wider θ vector.
  * `Unet` is NOT exported by v0.2.1, so the docs' image-to-image example is unreachable here
    and D-01's finite-dimensional vector head is the only route the pinned library supports.
  * `q` is a CONSTRUCTED `NormalisingFlow` INSTANCE passed POSITIONALLY -- the `q =` keyword
    convenience form of `PosteriorEstimator` expects a TYPE, not an instance. This is the
    Phase-4 API gotcha preserved verbatim from `spike/npe/architecture.jl:108-114`.

CPU-only is enforced at the CALL SITES (`use_gpu = false` on every NeuralEstimators call); the
architecture itself is device-agnostic.
"""
function build_p12_estimator(; G::Integer = P12_G, D::Integer = P12_D_MINISPIKE,
                             dstar::Integer = 128,
                             num_coupling_layers::Integer = 10,
                             flow_depth::Integer = 2, flow_width::Integer = 128)
    # `D >= 2` is the flow's own requirement, not a taste call: a CouplingLayer must have a
    # non-empty half on each side of the split.
    D     >= 2 || throw(ArgumentError("build_p12_estimator: D must be ≥ 2, got $D"))
    dstar >= 1 || throw(ArgumentError("build_p12_estimator: dstar must be ≥ 1, got $dstar"))
    G     >= 2 || throw(ArgumentError("build_p12_estimator: G must be ≥ 2, got $G"))
    num_coupling_layers >= 1 || throw(ArgumentError(
        "build_p12_estimator: num_coupling_layers must be ≥ 1, got $num_coupling_layers"))

    network = build_p12_summary_net(; G = G, dstar = dstar)
    q = NormalisingFlow(D; num_summaries = dstar,
                        num_coupling_layers = num_coupling_layers,
                        depth = flow_depth, width = flow_width)
    return PosteriorEstimator(network, q)
end

# =============================================================================================
# 5. The theta row layout (D-01 / D-08 / R-1) -- fixed here, consumed by every later plan
# =============================================================================================
#
# ROW ORDER, AND THE APPENDED-LAST DISCIPLINE `spike/simulator/prior.jl:87-92` established:
#
#   1                c0            global DCT DC coefficient (R-1: a genuinely SAMPLED parameter)
#   2 .. 1+K         c1..cK        deviation DCT coefficients in `p12_dct_order`, K = P12_K_DEV = 63
#   2+K .. 8+K       spillover, autofluorescence, label_efficiency, shift_dx, shift_dy,
#                    noise, chromatic_eps   (the 7 non-rho nuisances, EXISTING relative order)
#   9+K              r1            induced lag-1 correlation, D-08, APPENDED LAST
#
# D = 1 + K + 7 + 1 = 72 = P12_D_MINISPIKE.
#
# rho_true IS NOT A THETA ROW IN THE SPIKE LANE. It is replaced by the field, and the
# shipped-comparable scalar rho is a DERIVED read-time quantity,
# `ghat(quantile(MU_PRIOR, Phi(c0 / G)))`, which must be LABELLED AS DERIVED wherever it is
# reported (R-1). Naming the rows here, once, is what stops each downstream reader re-deriving
# a positional guess.

"""
    P12_THETA_ROWS(K = P12_K_DEV) -> Vector{Symbol}

The θ row order, in full: `:c0`, then `:c1 … :cK` (the deviation DCT coefficients in
`p12_dct_order`), then the seven non-ρ nuisances in their existing relative order, then `:r1`
appended LAST. `length(P12_THETA_ROWS()) == P12_D_MINISPIKE == 72`.
"""
function P12_THETA_ROWS(K::Integer = P12_K_DEV)
    K >= 0 || throw(ArgumentError("P12_THETA_ROWS: K must be ≥ 0, got $K"))
    rows = Symbol[:c0]
    for k in 1:K
        push!(rows, Symbol("c", k))
    end
    append!(rows, (:spillover, :autofluorescence, :label_efficiency,
                   :shift_dx, :shift_dy, :noise, :chromatic_eps))
    push!(rows, :r1)                       # APPENDED LAST (D-08)
    return rows
end

"""
    p12_theta_index(sym::Symbol; K = P12_K_DEV) -> Int

The θ row index of `sym`, derived from `P12_THETA_ROWS` rather than from a re-derived
arithmetic guess. Throws `ArgumentError` on an unknown name.
"""
function p12_theta_index(sym::Symbol; K::Integer = P12_K_DEV)
    rows = P12_THETA_ROWS(K)
    i = findfirst(==(sym), rows)
    i === nothing && throw(ArgumentError(
        "p12_theta_index: $sym is not a θ row at K = $K (rows are :c0, :c1..:c$K, the 7 " *
        "nuisances, :r1)"))
    return i
end

# The convenience positions, at the pre-registered K = P12_K_DEV. Every downstream reader takes
# these rather than writing `2:64` or `65:71` inline, so a K change is one edit, not a search.
const p12_row_c0        = 1                                  # the global DCT term
const p12_rows_dev      = 2:(1 + P12_K_DEV)                   # 2:64   -- 63 deviation coefficients
const p12_rows_nuisance = (2 + P12_K_DEV):(8 + P12_K_DEV)     # 65:71  -- the 7 non-rho nuisances
const p12_row_r1        = 9 + P12_K_DEV                       # 72     -- appended last

@assert length(P12_THETA_ROWS()) == P12_D_MINISPIKE
@assert first(P12_THETA_ROWS()) === :c0 && last(P12_THETA_ROWS()) === :r1
@assert length(p12_rows_dev) == P12_K_DEV
@assert length(p12_rows_nuisance) == 7
@assert p12_row_r1 == P12_D_MINISPIKE

# =============================================================================================
# 6. The recorded z-scoring arm (R-3, T-12-19)
# =============================================================================================
#
# PER-ROW Z-SCORING IS THE DEFAULT AND IT IS THE ARM IN FORCE. It matches the existing frozen
# `zt` and the existing loader exactly, so nothing about the standardization moved when the
# summary gained a spatial view.
#
# The declared ALTERNATIVE is `:shared_scalar` -- one mean/sd shared across all G^2 continuous
# rows, which would preserve translation equivariance instead of applying a per-pixel affine
# map that fights a convolution's weight sharing. It is taken ONLY if the CNN underperforms the
# MLP control, and taking it is a recorded deviation.
#
# WHICHEVER IS USED MUST BE RECORDED. An unrecorded choice here makes every later comparison
# unreadable, because the two arms differ in a way no downstream artifact would show. EVERY
# PERSISTED PHASE-12 ARTIFACT MUST CARRY THIS SYMBOL.
const P12_ZSCORE_ARM = :per_row
