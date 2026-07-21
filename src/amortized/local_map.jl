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
# src/amortized/local_map.jl --- windowed sub-tile local colocalization map (PROD-02, D-04).
#
# The HONEST Phase-7 answer to "local localisation": run the ALREADY-GATED estimator of one
# fixed grid over image SUB-TILES and read a per-tile Δρ. Nothing is retrained and no new gate
# is run — every tile is scored by the SAME frozen bundle, so the map inherits that grid's
# ship-gate exactly.
#
# WHAT THIS IS NOT: a calibrated per-region posterior. A finer patch grid lengthens the summary
# vector but the NPE output stays a single GLOBAL ρ (07-RESEARCH Finding 2), and windowing does
# not add a spatial prior, borrow strength between neighbours, or propagate per-region
# uncertainty. The true per-region GP/CAR map — `SpatialColocResult` / `delta_rho_map` /
# `uncertainty_map` — is Phase 12 and is DELIBERATELY not implemented here (the extension-point
# sketch in src/results.jl stays untouched).
#
# WAVE-6 REGISTRATION BRIDGE: `_lazy_load_from_artifact!` (Artifacts.toml, content-hashed) is
# wired by a LATER Phase-7 plan, so `_REGISTRY` is empty at this point in the build order.
# `_ensure_grid_registered` therefore loads the bundle IN-PROCESS from the local, in-repo
# `artifacts/grid_G/` files produced by that grid's own gate pipeline (a trusted, dev-time
# source — T-7-01) and `register!`s it, so `estimator_for(G)` works without touching the
# Artifacts path.
#############################################################################################

"""
    LOCAL_MAP_SENTINEL :: Float64

The documented finite sentinel written into `LocalColocMap.delta_rho` for a sub-tile that could
NOT be scored — a tile smaller than the patch grid, or one whose patches all fall below the
`≥15`-survivor floor of `correlation`, or one whose posterior read came back non-finite. The
value is `0.0` ("no evidence of a sample-vs-control difference"), and EVERY sentinel tile also
carries `ood_flag = true`, so a sentinel is never silently indistinguishable from a measured
Δρ of zero (T-7-04: a degenerate tile must never crash the map, and must never masquerade as a
finding).
"""
const LOCAL_MAP_SENTINEL = 0.0

"""
    LocalColocMap

A COARSE windowed colocalization readout: one Δρ + one OOD flag per image sub-tile, all produced
by a single frozen, already-gated estimator bundle. NOT a `SpatialColocResult` (Phase 12) and
deliberately NOT an `AbstractColocResult` — it carries no posterior draws, no Bayes factor and no
per-region uncertainty.

# Fields
- `grid::Int`: the patch grid `G` of the frozen bundle each tile was scored with (the registry key).
- `tiles::Tuple{Int,Int}`: the `(rows, cols)` sub-tile layout of the map.
- `delta_rho::Matrix{Float64}`: `rows×cols` per-tile Δρ; unscorable tiles hold `LOCAL_MAP_SENTINEL`.
- `ood_flag::Matrix{Bool}`: `rows×cols` per-tile flag — `true` for a fired OOD channel OR a sentinel tile.
- `meta::NamedTuple`: provenance (`N`, `ood_score`, `n_sentinel`, `tile_size`, `channels`,
  `sentinel`, `ood_thr`). `ood_thr` is the recorded pre-registered operating point the flags were
  decided against; `nothing` means none was available, so a `false` flag on a non-sentinel tile
  means "NOT CHECKED", not "in distribution".
"""
struct LocalColocMap
    grid      :: Int
    tiles     :: Tuple{Int,Int}
    delta_rho :: Matrix{Float64}
    ood_flag  :: Matrix{Bool}
    meta      :: NamedTuple
end

"""
    _recorded_ood_threshold(artifact_dir, grid) -> Union{Float64, Nothing}

Read the grid's ALREADY-PRE-REGISTERED OOD operating point out of the gate report
(`artifact_dir/gate_report_G.jld2` → `report.ood.id_threshold`, as written by the D-05 ship-gate
through `write_gate_report`). Nothing is re-run, re-scored or re-chosen here: this only recovers
the number the gate committed BEFORE seeing misspecified data.

Returns `nothing` — never a substitute value — when the report is absent, unreadable, carries no
OOD section, or records a non-finite threshold. A missing operating point must leave the flag
HONESTLY INERT rather than silently default to something arbitrary.
"""
function _recorded_ood_threshold(artifact_dir, grid)
    path = joinpath(artifact_dir, "gate_report_$(Int(grid)).jld2")
    isfile(path) || return nothing
    rep = try
        d = JLD2.load(path)
        get(d, "report", nothing)
    catch e
        e isa InterruptException && rethrow()
        nothing
    end
    (rep isa NamedTuple && haskey(rep, :ood)) || return nothing
    o = rep.ood
    (o isa NamedTuple && haskey(o, :id_threshold)) || return nothing
    t = o.id_threshold
    return (t isa Real && isfinite(t)) ? float(t) : nothing
end

"""
    _activate_recorded_ood_threshold(b::EstimatorBundle, artifact_dir) -> EstimatorBundle

Wire the grid's RECORDED pre-registered OOD operating point into the bundle's `ood_nulls` as
`:thr`, so `ood_verdict` actually decides instead of returning `flag = false` for every tile.

The persisted `ood_nulls_G.jld2` carries only the frozen `:density` fit (`_train_grid_pipeline`
step 5) — no `:thr`. Without one, `ood_verdict` short-circuits to `flag = false` and the per-tile
flag is VACUOUS: it then fires only structurally (sentinel tiles), never as a detector, which
silently understates risk.

**SCALE GUARD (mandatory).** The gate's `id_threshold` is the `OOD_ID_QUANTILE` quantile of RAW
`maha_score(density, ·)` values (`ood_gate` in `test/gate/run_gate.jl`). `ood_verdict` compares
`thr` against its FUSED statistic, which equals that same raw density Mahalanobis score only when
(a) there is no `:zref` (else channels are robust-z rescaled) and (b) no other channel is active
(`:noise`, `:model` would enter the `max` on their own raw scales). If either holds, the recorded
number is NOT comparable, so the threshold is NOT installed and the flag stays inert with a loud
warning — a mis-scaled threshold is worse than an inert flag.

A bundle that already carries a `:thr` is returned untouched (never overwritten).
"""
function _activate_recorded_ood_threshold(b::EstimatorBundle, artifact_dir)
    nulls = b.ood_nulls
    _has(k) = haskey(nulls, k) && getfield(nulls, k) !== nothing
    _has(:thr) && return b                       # already carries its own operating point

    if _has(:zref) || _has(:noise) || _has(:model)
        @warn "_activate_recorded_ood_threshold: ood_nulls carries extra channels/zref, so the " *
              "gate's raw density-Mahalanobis `id_threshold` is NOT on the same scale as the " *
              "fused ood_verdict statistic — the recorded operating point was NOT installed and " *
              "the per-tile OOD flag stays INERT (structural/sentinel firing only)." grid = b.grid
        return b
    end

    thr = _recorded_ood_threshold(artifact_dir, b.grid)
    if thr === nothing
        @warn "_activate_recorded_ood_threshold: no recorded OOD operating point in " *
              "$(joinpath(artifact_dir, "gate_report_$(b.grid).jld2")) — the per-tile OOD flag " *
              "is INERT (it can only fire structurally, on sentinel tiles)." grid = b.grid maxlog = 1
        return b
    end

    return EstimatorBundle(b.grid, b.npe, b.ratio, merge(nulls, (; thr = thr)),
                           b.zt, b.θzt, b.calibration)
end

"""
    _ensure_grid_registered(grid::Integer; artifact_dir = ...) -> EstimatorBundle

Make `estimator_for(grid)` work IN-PROCESS by loading the grid's frozen bundle straight from the
local `artifact_dir` (`npe_G.jld2` / `ratio_G.jld2` / `ood_nulls_G.jld2`, as written by
`_train_grid_pipeline`) and `register!`ing it. A no-op returning the registered bundle when
`grid` is already in `_REGISTRY`.

This is the WAVE-6 bridge over the not-yet-wired Artifacts path: it does NOT call
`_lazy_load_from_artifact!` and reads no `Artifacts.toml`. The files are in-repo outputs of this
package's own gate pipeline (trusted dev-time inputs, T-7-01); the public, content-hashed
distribution path is wired by a later Phase-7 plan. A missing artifact errors loudly and names
`train_and_register(grid)`.

The registered bundle also gets the grid's RECORDED, pre-registered OOD operating point wired in
via `_activate_recorded_ood_threshold` (read from `gate_report_G.jld2`, scale-guarded), so the
per-tile `ood_verdict` flag is operative rather than vacuously `false`.
"""
function _ensure_grid_registered(grid::Integer;
                                 artifact_dir = _grid_dir(default_artifacts_root(), grid))
    g = Int(grid)
    haskey(_REGISTRY, g) && return _REGISTRY[g]
    npe_path   = joinpath(artifact_dir, "npe_$(g).jld2")
    ratio_path = joinpath(artifact_dir, "ratio_$(g).jld2")
    ood_path   = joinpath(artifact_dir, "ood_nulls_$(g).jld2")
    for p in (npe_path, ratio_path, ood_path)
        isfile(p) || error("_ensure_grid_registered: no artifact at $p — the $(g)×$(g) bundle " *
                           "has not been produced yet. Build one with `train_and_register($g)`.")
    end
    b = _bundle_from_artifacts(g, npe_path, ratio_path, ood_path)
    return register!(_activate_recorded_ood_threshold(b, artifact_dir))
end

# Split one channel matrix into an r×c grid of sub-tile matrices, reusing the package's own
# `patch(img, num_patches_x, num_patches_y)` tiler (src/colocalization.jl) — same trimming rule
# as the summary path, so a sub-tile is exactly the block the correlation summary would see.
function _sub_tiles(m::AbstractMatrix, r::Integer, c::Integer)
    P = patch(m, Int(r), Int(c))                       # r×c×px×py
    return [P[i, j, :, :] for i in 1:Int(r), j in 1:Int(c)]
end

# Build a 2-channel MultiChannelImage for ONE sub-tile, carrying the PARENT image's Otsu
# thresholds (a tile-local Otsu would silently re-threshold each window and make tiles
# incomparable). Mirrors the `build_mci` convenience shape.
function _tile_mci(x::AbstractMatrix, y::AbstractMatrix, parent::MultiChannelImage,
                   ch::AbstractVector{<:Integer}, name::AbstractString)
    return MultiChannelImage(
        [x, y],
        [parent.channels[ch[1]], parent.channels[ch[2]]],
        String(name),
        ["", ""],
        size(x),
        [parent.otsu_threshold[ch[1]], parent.otsu_threshold[ch[2]]],
    )
end

# Encoded (raw) summary for one sub-tile + whether ANY patch survived the ≥15-survivor floor.
function _tile_summary(mci::MultiChannelImage, grid::Int)
    Zraw = encode_d01(patch_summary(mci, grid))
    _, mask = _summary_row_partition(:min, length(Zraw))
    return Zraw, any(>(0.5), @view Zraw[mask])
end

"""
    local_coloc_map(img, control, channels; grid = 8, tiles = (4, 4), N = 2000,
                    artifact_dir = ..., use_gpu = false) -> LocalColocMap

COARSE WINDOWED READOUT — **not** a calibrated per-region colocalization map.

Partition `img` and `control` (two `MultiChannelImage`s) into a `tiles = (rows, cols)` grid of
sub-tiles, and score each sample tile against the CORRESPONDING control tile with the frozen,
already-gated `grid×grid` estimator bundle: the tile's `patch_summary`/`encode_d01` summary is
standardized with the bundle's FROZEN `zt`, read through `delta_rho` (the D-03 Monte-Carlo
difference of two independent posterior passes), and passed to `ood_verdict`. `channels` selects
the two channels (e.g. `[1, 2]`) compared in every tile.

No training and no new ship-gate happen here: `_ensure_grid_registered(grid)` registers the
bundle in-process from the local artifacts, `estimator_for(grid)` retrieves it, and the map
therefore inherits that grid's recorded ship-gate verdict verbatim — including its caveats.

## What this output is, and is not

Each tile is an INDEPENDENT global-ρ read on a smaller window. There is no spatial prior, no
borrowing of strength between neighbouring tiles, and no per-region uncertainty: the returned
`delta_rho` is a matrix of point estimates, not a posterior field. Smaller tiles also carry fewer
pixels per patch, so tiles drift toward — and can fall below — the `≥15`-survivor floor, where
the tile is marked with `LOCAL_MAP_SENTINEL` and flagged rather than reported. The calibrated
per-region map (GP/CAR lattice, `SpatialColocResult` with `delta_rho_map` / `uncertainty_map`) is
Phase 12 and is not implemented by this function.

## Returns

A [`LocalColocMap`](@ref) with exactly `rows×cols` entries in both `delta_rho` and `ood_flag`;
every `delta_rho` entry is finite (a measured Δρ or `LOCAL_MAP_SENTINEL`) — a degenerate tile
never throws (T-7-04).
"""
function local_coloc_map(img::MultiChannelImage, control::MultiChannelImage,
                         channels::AbstractVector{<:Integer};
                         grid::Integer = 8,
                         tiles::Tuple{Int,Int} = (4, 4),
                         N::Integer = 2000,
                         artifact_dir = _grid_dir(default_artifacts_root(), grid),
                         use_gpu::Bool = false)
    length(channels) == 2 ||
        throw(ArgumentError("local_coloc_map: `channels` must select exactly 2 channels, got $(length(channels))"))
    r, c = tiles
    (r >= 1 && c >= 1) ||
        throw(ArgumentError("local_coloc_map: `tiles` must be ≥1 in both directions, got $tiles"))

    G = Int(grid)
    _ensure_grid_registered(G; artifact_dir = artifact_dir)
    b = estimator_for(G)

    xs = _sub_tiles(img.data[channels[1]],     r, c)
    ys = _sub_tiles(img.data[channels[2]],     r, c)
    xc = _sub_tiles(control.data[channels[1]], r, c)
    yc = _sub_tiles(control.data[channels[2]], r, c)

    Δ     = fill(LOCAL_MAP_SENTINEL, r, c)
    flag  = falses(r, c)
    score = fill(NaN, r, c)
    n_sentinel = 0

    for j in 1:c, i in 1:r
        tile_s = _tile_mci(xs[i, j], ys[i, j], img,     channels, "$(img.name)_tile_$(i)_$(j)")
        tile_c = _tile_mci(xc[i, j], yc[i, j], control, channels, "$(control.name)_tile_$(i)_$(j)")

        # Sentinel path 1: the sub-tile is smaller than the patch grid, so `patch` cannot even
        # form G×G patches. Never attempt the read.
        if minimum(size(xs[i, j])) < G || minimum(size(xc[i, j])) < G
            flag[i, j] = true; n_sentinel += 1
            continue
        end

        Zs_raw, ok_s = _tile_summary(tile_s, G)
        Zc_raw, ok_c = _tile_summary(tile_c, G)

        # Sentinel path 2: below the ≥15-survivor floor — no patch in the tile survived, so the
        # summary carries no correlation information at all.
        if !(ok_s && ok_c)
            flag[i, j] = true; n_sentinel += 1
            continue
        end

        Zs = standardize_summary(Zs_raw, b.zt, :min)
        Zc = standardize_summary(Zc_raw, b.zt, :min)
        d  = delta_rho(b.npe, Zs, Zc, b.θzt; N = N, use_gpu = use_gpu)

        # Sentinel path 3: a non-finite posterior read (finite-guard, T-7-04).
        if !isfinite(d)
            flag[i, j] = true; n_sentinel += 1
            continue
        end

        Δ[i, j] = d
        v = ood_verdict(b.ood_nulls, Zs, Zc)
        score[i, j] = v.score
        flag[i, j]  = v.flag
    end

    # `ood_thr` is the operating point the flags were actually decided against — `nothing` means
    # NO threshold was available, i.e. the non-sentinel flags are inert, not "in distribution".
    thr = haskey(b.ood_nulls, :thr) ? b.ood_nulls.thr : nothing

    return LocalColocMap(G, (r, c), Δ, flag,
                         (; N = Int(N), ood_score = score, n_sentinel = n_sentinel,
                            tile_size = size(xs[1, 1]), channels = collect(Int, channels),
                            sentinel = LOCAL_MAP_SENTINEL, ood_thr = thr))
end
