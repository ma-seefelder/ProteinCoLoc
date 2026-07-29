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

# spike/p12/result.jl --- SpatialColocResultSpike, the spike-lane per-region result type (SPAT-05).
#
# (a) WHAT THIS FILE EXISTS TO FILL. src/results.jl:180-191 has carried a COMMENT-ONLY
#     `SpatialColocResult` sketch since Phase 7, and src/amortized/local_map.jl:56-79 states in its
#     own docstring that the shipped `LocalColocMap` "carries no posterior draws, no Bayes factor
#     and no per-region uncertainty". A per-region Delta-rho map WITH a per-region uncertainty map
#     is this phase's product, and it is also the D-13 descope deliverable -- so the type must
#     exist and be tested regardless of how the Stage-1 and Stage-2 gates land.
#
# (b) WHERE THE EXECUTABLE TYPE LIVES, AND WHY THE SKETCH IS NOT EDITED. D-01 and CLAUDE.md scope
#     all Phase-12 work to spike/, and CLAUDE.md requires src/ to stay provably untouched during
#     the spike. So the EXECUTABLE subtype is defined HERE, in spike/p12/, reached to its supertype
#     by a READ-ONLY include of src/results.jl -- the same read-only reach spike/contract.jl:46-47
#     already uses for src/LoadImages.jl and src/colocalization.jl. The comment-only sketch block
#     at src/results.jl:180-191 stays BYTE-UNCHANGED, and its identifier `SpatialColocResult` stays
#     UNCLAIMED, so the future productionization phase inherits a free name rather than a conflict.
#     spike/p13/result.jl solved exactly this problem for Phase 13 and is the precedent mirrored
#     here structurally.
#
# (c) PHASE 7 D-02 IS STILL SATISFIED. D-02 requires every new result variant to slot in as a NEW
#     subtype of the real `AbstractColocResult`, never as bolted-on fields on an existing struct.
#     That is structurally what happens here: `SpatialColocResultSpike <: AbstractColocResult`
#     against the ACTUAL supertype from src/results.jl, and the four shared accessors are
#     implemented for it. Defining it in spike/ changes nothing about that relationship.
#
# (d) DECOUPLING (hard constraint, CLAUDE.md): spike-local code. The reach into src/ is a read-only
#     `include`; no src/ byte changes, no package is added, and nothing here is imported by the
#     shipped package.
#
# WHAT THIS FILE DELIBERATELY DOES NOT CONTAIN. The `p12_coloc_map` FUNCTION that produces one of
# these from a trained estimator and an image needs the posterior read path and lands with the
# coverage machinery; this file defines the CONTAINER only. Defining it early is deliberate: two
# later plans write into it, and a container invented mid-run acquires fields by accretion -- which
# is exactly what D-02's "new subtypes, never bolted-on fields" rule exists to prevent.
#
# ORDER MATTERS (spike/contract.jl:34-39, the documented load order): StatsBase + Statistics must
# be in scope before src/colocalization.jl is reached transitively, and Images before
# src/LoadImages.jl (its convenience constructor calls Images.otsu_threshold). Both are reached
# through spike/validation/sbc.jl -> harness.jl -> contract.jl below, so the `using` lines come
# FIRST and the includes after.
using StatsBase      # corspearman/corkendall for the transitively-reached correlation() Dict
using Statistics     # mean, quantile
using Images         # otsu_threshold (transitive src/LoadImages.jl requirement)

# --- Guarded includes, in dependency order --------------------------------------------------
# READ-ONLY reach into src/ for the supertype and `_iface_error`.
# NOTE THE PATH DEPTH: from spike/p12/ the repo root is TWO levels up, not three.
isdefined(@__MODULE__, :AbstractColocResult) ||
    include(joinpath(@__DIR__, "..", "..", "src", "results.jl"))
# The shared, gate-lineage calibration surface: `CalibrationResult` + `_bin_calibration`. Guarded
# on the FUNCTION rather than the struct, matching spike/p13/result.jl:83-84. A SECOND calibration
# container is NEVER defined here -- the `calibration` field carries the shared one.
#
# THIS INCLUDE COMES BEFORE THE PHASE-12 PRE-REGISTRATION ON PURPOSE (DEF-12-02 / DEF-12-03).
# `spike/validation/p12_consts.jl` binds `:P11_DEV_SEED` and `:P13_DEV_SEED` -- the very sentinels
# `p11_consts.jl` and `spike/p13/consts.jl` guard their own Tier-1 blocks on -- so a process that
# loads the Phase-12 pre-registration first can silently skip a FOREIGN one. Loading the foreign
# surfaces first is the standing workaround. `spike/validation/harness.jl:51-56` records the same
# hazard for `:VAL_MASTER_SEED` and is why its own guard is keyed on `:SBC_M`.
isdefined(@__MODULE__, :_bin_calibration) ||
    include(joinpath(@__DIR__, "..", "validation", "sbc.jl"))
# The Phase-12 model surface: `P12_ZSCORE_ARM` (R-3, the recorded z-scoring arm) and the theta row
# map. It pulls the Tier-1 pre-registration (`P12_G`, `P12_K_DEV`, ...) transitively.
isdefined(@__MODULE__, :P12_ZSCORE_ARM) ||
    include(joinpath(@__DIR__, "..", "npe", "p12_architecture.jl"))
# The lattice arithmetic: `p12_idct_vec` + `p12_dct_order`, for the field reconstruction below.
isdefined(@__MODULE__, :p12_idx) ||
    include(joinpath(@__DIR__, "..", "simulator", "p12_lattice.jl"))

# THE GUARD BLOCK COVERS ONLY THE `const`s AND THE `struct`, NOT THE FUNCTIONS, AND THAT SPLIT IS
# DELIBERATE. Julia 1.12 DROPS every docstring written inside an `if ... end` block -- the parser
# emits the `Core.@doc` call but the docsystem never registers it, so `@doc f` returns `nothing`
# for a function documented inside a guard (verified on 1.12.6, spike/p13/result.jl:86-93). R-1
# requires the derived scalar rho to be LABELLED DERIVED where a reader will see it, and the
# accompanying test asserts that label is retrievable from the RUNNING system rather than from the
# source text, so the functions live at TOP LEVEL where the docsystem can see them. Method
# redefinition under a re-include is silent and harmless in Julia; only `const` and `struct`
# redefinition would warn or throw, and those are what the guard protects.
if !isdefined(@__MODULE__, :SpatialColocResultSpike)

# --- The sentinel, and the promise that it does not drift ------------------------------------

"""
    P12_SPATIAL_SENTINEL :: Float64

The finite sentinel written into a region of any of this type's three maps when that region could
NOT be scored. It is the SAME VALUE and the SAME DISCIPLINE as `LOCAL_MAP_SENTINEL`
(`src/amortized/local_map.jl:42-53`): the value is `0.0` and EVERY sentinel region also carries a
`true` OOD flag, so a sentinel is never silently indistinguishable from a measured Delta-rho of
zero (T-7-04).

It is REPRODUCED here rather than imported because `src/` is read-only in this phase and
`local_map.jl` is not loaded by the spike lane. The two must not drift, so a testset reads the
literal out of `local_map.jl` AS TEXT and asserts it parses to this value -- trusting a comment is
exactly what this project's seed-disjointness rule already rejects.
"""
const P12_SPATIAL_SENTINEL = 0.0

"""
    P12_IDENTITY_ATOL :: Float64

The stated absolute tolerance of the STRUCTURAL IDENTITY
`region_delta_rho ~= region_rho_sample - region_rho_control`, checked elementwise by the
constructor. It is a float-rounding allowance on a difference of two O(1) correlations, not a
modelling slack: anything a Monte-Carlo difference of the two carried maps produces agrees to far
better than this.
"""
const P12_IDENTITY_ATOL = 1e-9

"""
    P12_RESULT_META_KEYS

The EXACT key set a [`SpatialColocResultSpike`](@ref)'s `meta` must carry. Guarded as an exact
match (not a subset) by [`_p12_check_result_meta`](@ref), so a missing key AND an extra key both
throw and name D-02 rather than being stored silently.

`:delta_rho_available` is in this list, and the list is EXACT-MATCH, so adding the flag to
[`p12_result_meta`](@ref) without adding it here would make the guard reject every result the
constructor is handed -- a missing key on one side and an extra key on the other, from the same
change. That coupling is the point: the two declaration sites cannot drift apart quietly.
"""
const P12_RESULT_META_KEYS = (:zscore_arm, :arm, :r1_posterior_median, :derived_rho,
                              :ood_threshold, :net_path, :delta_rho_available, :schema_version)

# --- The result type --------------------------------------------------------------------------

"""
    SpatialColocResultSpike(grid, region_delta_rho, region_sd, region_rho_sample, region_sd_sample,
                            region_rho_control, region_sd_control, ood_sample, ood_control, ood,
                            region_draws, calibration, meta)
    SpatialColocResultSpike(grid, maps, region_draws, calibration, meta)

The Phase-12 per-region colocalization map: a per-region Delta-rho lattice WITH a per-region
uncertainty lattice, as a new subtype of the real `AbstractColocResult` (Phase 7 D-02), defined in
the spike lane per D-01. `src/results.jl` is byte-unchanged and its comment-only sketch is
superseded by this file rather than edited.

# THE RESULT CARRIES THREE MAPS

`region_delta_rho` is the PRIMARY named deliverable (the `src/results.jl:180-191` sketch names it),
computed as the per-region Monte-Carlo difference mirroring `src/amortized/infer.jl:117-121`.
`region_rho_sample` and `region_rho_control` are ADDITIONAL exposed artifacts, never a replacement:
a reader benefits from seeing WHERE a difference comes from rather than being handed a contrast
they cannot interrogate. Exposing them costs ZERO extra forward passes -- a Delta-rho read already
performs both single-stack passes internally, so these are values the computation already held and
threw away.

# THE GUARD IS A STRUCTURAL IDENTITY, NOT A RANGE

The constructor asserts, elementwise and to `P12_IDENTITY_ATOL`,

    region_delta_rho ~= region_rho_sample - region_rho_control

over every region where none of the three is a sentinel. This REPLACED -- did not tighten -- an
earlier `[-2, 2]` range check on the difference. The old bound admitted `rho` in `[-1, 1]`
silently, so the one check that could have caught a rho map mislabelled as a Delta-rho map is
precisely what hid it; and tightening is no repair either, because with a TRUE Delta-rho the range
genuinely IS `[-2, 2]`. A range cannot distinguish the two quantities at all. An identity can, and
it becomes available only because all three maps are now shipped.

The `[-1, 1]` and `[-2, 2]` bounds are RETAINED as SANITY checks, explicitly labelled as such in
the code and in the error messages, and are never the guard.

# SENTINEL SEMANTICS, DEFINED PER MAP

`local_map.jl:48`'s wording is about the DIFFERENCE and is meaningless for a single-stack map:

| Map | A sentinel there means |
|---|---|
| `region_delta_rho` | "no evidence of a sample-vs-control difference" -- the inherited meaning |
| `region_rho_sample` | "this region was not scorable in THIS stack" -- measurability, not zero correlation |
| `region_rho_control` | as above, for the control stack |

The invariant holds in all three: a sentinel region carries `ood = true` and a `NaN` `sd`, so no
sentinel is ever silently indistinguishable from a measured zero. A region unscorable in EITHER
stack is unscorable in the difference, so the difference's sentinel set is the UNION of the two
single-stack sets -- `ood == ood_sample .| ood_control`, enforced here rather than hoped for at the
call site.

# THE SINGLE-STACK CONSTRUCTION IS A LEGAL, MEANINGFUL STATE

When there is no control stack at all -- the real-image arm, where two physical specimens are
NOT exchangeable and so must not be paired -- the result is built with `region_rho_control`
entirely `P12_SPATIAL_SENTINEL`, `region_sd_control` entirely `NaN`, `ood_control` entirely `true`,
and therefore by the union rule `ood` entirely true and `region_delta_rho` entirely sentinel. The
constructor needs no special case to admit this: the structural identity is scoped to regions where
no map is a sentinel and so holds vacuously, `P12_SPATIAL_SENTINEL = 0.0` lies inside both sanity
ranges, and the sentinel/OOD pairing is exact in all three maps.

That is the point of writing it down: this is how *"Delta-rho is not computable on this data"*
stops being a sentence in a report and becomes a PROPERTY OF THE RETURNED OBJECT. `meta` carries
`delta_rho_available::Bool` so a consumer reads the mode directly instead of inferring it from a
sentinel pattern. **A single-stack result is NOT a Delta-rho of zero and must never be summarized
as one.**

# Fields
- `grid::Int`: patch grid `G`; every map is `G x G`.
- `region_delta_rho::Matrix{Float64}`: the PRIMARY map, in rho units; sentinel where unscorable.
- `region_sd::Matrix{Float64}`: per-region sd OF THE DIFFERENCE; `NaN` exactly where `ood` is true.
- `region_rho_sample::Matrix{Float64}` / `region_sd_sample::Matrix{Float64}`: the sample stack.
- `region_rho_control::Matrix{Float64}` / `region_sd_control::Matrix{Float64}`: the control stack.
- `ood_sample` / `ood_control` / `ood::Matrix{Bool}`: per-map flags; `ood` is the UNION.
- `region_draws::Union{Nothing,Array{Float64,3}}`: `G x G x N` draws OF THE DIFFERENCE, or nothing.
- `calibration::Union{Nothing,CalibrationResult}`: the shared calibration container, or nothing.
- `meta::NamedTuple`: exactly `P12_RESULT_META_KEYS`; build it with [`p12_result_meta`](@ref).

The convenience constructor takes the nine-matrix working bundle from
[`p12_region_maps`](@ref) in place of the nine positional matrices, so the argument order is
written down once instead of at every call site.
"""
struct SpatialColocResultSpike <: AbstractColocResult
    grid               :: Int
    # --- the PRIMARY named deliverable (the src/results.jl:180-191 sketch) ---
    region_delta_rho   :: Matrix{Float64}   # GxG, rho units, sentinel where unscorable
    region_sd          :: Matrix{Float64}   # GxG, sd of the DIFFERENCE, NaN where unscorable
    # --- the two single-stack maps, exposed so a reader can see WHERE the difference comes from --
    region_rho_sample  :: Matrix{Float64}   # GxG, sentinel where the SAMPLE stack was unscorable
    region_sd_sample   :: Matrix{Float64}
    region_rho_control :: Matrix{Float64}   # GxG, sentinel where the CONTROL stack was unscorable
    region_sd_control  :: Matrix{Float64}
    # --- per-map OOD, because "unscorable" is a per-stack fact ---
    ood_sample         :: Matrix{Bool}
    ood_control        :: Matrix{Bool}
    ood                :: Matrix{Bool}      # the DIFFERENCE's flag == ood_sample .| ood_control
    region_draws       :: Union{Nothing, Array{Float64,3}}
    calibration        :: Union{Nothing, CalibrationResult}
    meta               :: NamedTuple

    function SpatialColocResultSpike(grid::Integer,
                                     region_delta_rho::AbstractMatrix{<:Real},
                                     region_sd::AbstractMatrix{<:Real},
                                     region_rho_sample::AbstractMatrix{<:Real},
                                     region_sd_sample::AbstractMatrix{<:Real},
                                     region_rho_control::AbstractMatrix{<:Real},
                                     region_sd_control::AbstractMatrix{<:Real},
                                     ood_sample::AbstractMatrix{Bool},
                                     ood_control::AbstractMatrix{Bool},
                                     ood::AbstractMatrix{Bool},
                                     region_draws::Union{Nothing,AbstractArray{<:Real,3}},
                                     calibration::Union{Nothing,CalibrationResult},
                                     meta::NamedTuple)
        G = Int(grid)
        G >= 2 || throw(ArgumentError(
            "SpatialColocResultSpike: `grid` must be >= 2 -- a 1x1 lattice carries no spatial " *
            "structure for a per-region map to describe. Got grid = $G."))

        dlt = Matrix{Float64}(region_delta_rho)
        sd  = Matrix{Float64}(region_sd)
        smp = Matrix{Float64}(region_rho_sample)
        sds = Matrix{Float64}(region_sd_sample)
        ctl = Matrix{Float64}(region_rho_control)
        sdc = Matrix{Float64}(region_sd_control)
        osm = Matrix{Bool}(ood_sample)
        oct = Matrix{Bool}(ood_control)
        odf = Matrix{Bool}(ood)

        # (1) Shape. Checked first so every later index is known to be in range.
        for (nm, m) in (("region_delta_rho",   size(dlt)), ("region_sd",          size(sd)),
                        ("region_rho_sample",  size(smp)), ("region_sd_sample",   size(sds)),
                        ("region_rho_control", size(ctl)), ("region_sd_control",  size(sdc)),
                        ("ood_sample",         size(osm)), ("ood_control",        size(oct)),
                        ("ood",                size(odf)))
            m == (G, G) || throw(ArgumentError(
                "SpatialColocResultSpike: `$nm` must be exactly grid x grid = $G x $G, got $m."))
        end

        # (2) ==== THE STRUCTURAL IDENTITY -- THIS CONSTRUCTOR'S MOST IMPORTANT CHECK ==========
        # Scoped to regions where NONE of the three is a sentinel: a sentinel is a refusal to
        # measure, not a measurement, so differencing two refusals asserts nothing.
        for i in eachindex(dlt)
            (dlt[i] == P12_SPATIAL_SENTINEL || smp[i] == P12_SPATIAL_SENTINEL ||
             ctl[i] == P12_SPATIAL_SENTINEL) && continue
            isapprox(dlt[i], smp[i] - ctl[i]; atol = P12_IDENTITY_ATOL, rtol = 0.0) || throw(
                ArgumentError(
                "SpatialColocResultSpike: THE STRUCTURAL IDENTITY FAILED at region " *
                "$(Tuple(CartesianIndices(dlt)[i])) -- region_delta_rho = $(dlt[i]) but " *
                "region_rho_sample - region_rho_control = $(smp[i] - ctl[i]) " *
                "(atol = $(P12_IDENTITY_ATOL)). " *
                "THIS CHECK REPLACED a [-2, 2] range check on the difference, and reinstating a " *
                "range in its place is a regression: the old bound admitted rho in [-1, 1] " *
                "SILENTLY, so the one check that could have caught a rho map mislabelled as a " *
                "Delta-rho map is exactly what hid it -- and tightening is no repair either, " *
                "because with a TRUE Delta-rho the range genuinely IS [-2, 2]. A range cannot " *
                "distinguish a Delta-rho map from a rho map at all; this identity can, and only " *
                "because all three maps are shipped."))
        end

        # (3) SANITY ranges. EXPLICITLY NOT THE GUARD -- see (2). They exist to catch a garbled
        # array, not to identify a quantity.
        for (nm, m, lo, hi, why) in
            (("region_rho_sample",  smp, -1.0, 1.0, "a correlation"),
             ("region_rho_control", ctl, -1.0, 1.0, "a correlation"),
             ("region_delta_rho",   dlt, -2.0, 2.0,
              "a DIFFERENCE OF TWO CORRELATIONS, which genuinely spans [-2, 2] -- do NOT " *
              "\"fix\" this bound to [-1, 1]"))
            for i in eachindex(m)
                (isfinite(m[i]) && lo <= m[i] <= hi) || throw(ArgumentError(
                    "SpatialColocResultSpike: SANITY range violated -- `$nm` at region " *
                    "$(Tuple(CartesianIndices(m)[i])) is $(m[i]), outside [$lo, $hi]. It is " *
                    "$why. This is a sanity bound, NOT the Delta-rho guard; the guard is the " *
                    "structural identity."))
            end
        end

        # (4) The sd/OOD pairing, PER MAP. `NaN` EXACTLY where the map's own flag is true. This is
        # half of what makes an unscorable region unmistakable, and it is checked here rather than
        # hoped for at the call site: a small sd on a failed region reads as a confident
        # measurement.
        for (snm, s, fnm, f) in (("region_sd",         sd,  "ood",         odf),
                                 ("region_sd_sample",  sds, "ood_sample",  osm),
                                 ("region_sd_control", sdc, "ood_control", oct))
            for i in eachindex(s)
                pos = Tuple(CartesianIndices(s)[i])
                if f[i]
                    isnan(s[i]) || throw(ArgumentError(
                        "SpatialColocResultSpike: `$snm` must be NaN EXACTLY where `$fnm` is " *
                        "true; at region $pos the flag is true but the sd is $(s[i]). An " *
                        "unscorable region must not carry a finite uncertainty."))
                else
                    (isfinite(s[i]) && s[i] > 0.0) || throw(ArgumentError(
                        "SpatialColocResultSpike: `$snm` must be finite and > 0 wherever `$fnm` " *
                        "is false; at region $pos the flag is false but the sd is $(s[i])."))
                end
            end
        end

        # (5) The sentinel/OOD pairing, PER MAP -- the other half. A near-zero value on a failed
        # region would otherwise read as "measured no difference" (T-7-04).
        for (mnm, m, fnm, f) in (("region_delta_rho",   dlt, "ood",         odf),
                                 ("region_rho_sample",  smp, "ood_sample",  osm),
                                 ("region_rho_control", ctl, "ood_control", oct))
            for i in eachindex(m)
                (!f[i] || m[i] == P12_SPATIAL_SENTINEL) || throw(ArgumentError(
                    "SpatialColocResultSpike: `$mnm` must be P12_SPATIAL_SENTINEL " *
                    "($(P12_SPATIAL_SENTINEL)) EXACTLY where `$fnm` is true; at region " *
                    "$(Tuple(CartesianIndices(m)[i])) the flag is true but the value is $(m[i]). " *
                    "A sentinel must never be silently indistinguishable from a measured zero."))
            end
        end

        # (6) THE UNION RULE. A difference is unscorable whenever EITHER side is, so the
        # difference's sentinel set is the union of the two single-stack sets. Checked here rather
        # than hoped for, because it is the one invariant that ties the three maps' degeneracy
        # together -- and it is what makes the no-control-stack construction legal by consequence
        # rather than by exception.
        odf == (osm .| oct) || throw(ArgumentError(
            "SpatialColocResultSpike: THE UNION RULE FAILED -- `ood` must equal " *
            "`ood_sample .| ood_control` exactly. Regions that disagree: " *
            "$(Tuple.(findall(odf .!= (osm .| oct))))."))

        # (7) Draws, if any, are a GxGxN stack of the DIFFERENCE.
        draws = if region_draws === nothing
            nothing
        else
            size(region_draws)[1:2] == (G, G) || throw(ArgumentError(
                "SpatialColocResultSpike: `region_draws` must be grid x grid x N = $G x $G x N, " *
                "got $(size(region_draws))."))
            Array{Float64,3}(region_draws)
        end

        # (8) The meta key set, exact-match.
        _p12_check_result_meta(meta)

        return new(G, dlt, sd, smp, sds, ctl, sdc, osm, oct, odf, draws, calibration, meta)
    end
end

end # if !isdefined(@__MODULE__, :SpatialColocResultSpike) -- consts + struct only

"""
    _p12_check_result_meta(meta::NamedTuple) -> NamedTuple

Throw an `ArgumentError` naming D-02 unless `meta` carries EXACTLY `P12_RESULT_META_KEYS` -- no
missing key, no extra key, in any order.

Modelled on the keyed-NamedTuple guard `_p13_as_three_way_logbf` (`spike/p13/result.jl:179`). An
exact-match key set is what makes the `meta` a contract rather than a bag: a caller who forgets
`:delta_rho_available` is told so at construction instead of producing a result whose mode cannot
be read, and a caller who invents a key is told so instead of having it silently stored where no
consumer will look for it.
"""
function _p12_check_result_meta(meta::NamedTuple)
    ks = Set(keys(meta))
    want = Set(P12_RESULT_META_KEYS)
    ks == want || throw(ArgumentError(
        "SpatialColocResultSpike (Phase-7 D-02): `meta` must carry EXACTLY the keys " *
        "$(P12_RESULT_META_KEYS); got $(keys(meta)). Missing: " *
        "$(Tuple(sort(collect(setdiff(want, ks))))). Extra: " *
        "$(Tuple(sort(collect(setdiff(ks, want))))). Build it with `p12_result_meta`."))
    return meta
end

"""
    SpatialColocResultSpike(grid, maps::NamedTuple, region_draws, calibration, meta)

Convenience constructor that splats the nine-matrix working bundle from [`p12_region_maps`](@ref)
into the positional constructor, in the bundle's declared order.

It exists so the nine-matrix argument order is written down ONCE instead of at every call site. A
seven-argument call is the OLD single-map shape and no longer type-checks -- which is the intended
outcome, because the single-map shape cannot express the structural identity.
"""
SpatialColocResultSpike(grid::Integer, maps::NamedTuple,
                        region_draws::Union{Nothing,AbstractArray{<:Real,3}},
                        calibration::Union{Nothing,CalibrationResult},
                        meta::NamedTuple) =
    SpatialColocResultSpike(grid,
                            maps.region_delta_rho, maps.region_sd,
                            maps.region_rho_sample, maps.region_sd_sample,
                            maps.region_rho_control, maps.region_sd_control,
                            maps.ood_sample, maps.ood_control, maps.ood,
                            region_draws, calibration, meta)

# --- The four shared accessors (the AbstractColocResult interface) ----------------------------

"""
    posterior_draws(r::SpatialColocResultSpike)

The per-region posterior draws backing the map: a `G x G x N` stack of DIFFERENCE draws, or
`nothing` when the result was built without them.

Note the layout difference from the shipped `AmortizedColocResult`, whose `posterior` is a
parameter-by-draw MATRIX: here the leading two axes are the LATTICE, because the quantity is a
field rather than a parameter vector.
"""
posterior_draws(r::SpatialColocResultSpike) = r.region_draws

"""
    is_ood(r::SpatialColocResultSpike) -> Bool

`true` if ANY region carries an OOD flag, i.e. `any(r.ood)`.

THE CONSERVATIVE DIRECTION IS THE ONLY HONEST ONE for a single-boolean summary of a lattice: a map
with one unscorable region is not an in-distribution map, and reporting it as one would be the
sentinel-reads-as-a-measurement failure at map scale. A consumer that needs the per-region picture
must read `r.ood` (and `r.ood_sample` / `r.ood_control`) directly.

**A `false` here is only meaningful if the operating point was actually checked.** Per
`src/amortized/local_map.jl:88-91`, a `false` flag under a `nothing` threshold means NOT CHECKED,
not "in distribution" -- branch on [`p12_ood_checked`](@ref) first.
"""
is_ood(r::SpatialColocResultSpike) = any(r.ood)

"""
    bayes_factor(r::SpatialColocResultSpike)

Falls through to the documented interface error: this result carries NO Bayes factor.

DO NOT FABRICATE ONE. The Phase-12 read is a posterior over a spatial field; no ratio estimator is
trained against it and no evidence quantity is computed anywhere in this phase. Synthesising a
plausible-looking number -- a threshold crossing, a posterior mass, a per-region tail probability
relabelled as evidence -- would put a figure in front of a reader that no measurement backs.
Letting `_iface_error` (`src/results.jl:47-49`) fire is the honest answer, and it names the
accessor, the type and the four accessors the interface requires. This mirrors
`spike/p13/result.jl`'s treatment of the Delta-rho it does not carry.
"""
bayes_factor(r::SpatialColocResultSpike) = _iface_error(r, :bayes_factor)

"""
    delta_rho(r::SpatialColocResultSpike) -> Real

Return the single scalar rho this result carries, which is **DERIVED** (R-1) from the global DCT
coefficient `c0` and is **not a sampled parameter**: it is reported only for continuity with the
shipped estimator, whose `rho_true` row this phase's theta does not have.

A reader who takes this number as the estimand has been misled, so the label goes where they will
read it -- in the first sentence, not in a footnote. In the spike lane there IS no `rho_true`
theta row; it is replaced by the field, and the shipped-comparable scalar is the read-time
quantity `ghat(quantile(MU_PRIOR, Phi(c0 / G)))` recorded in `meta.derived_rho`.

**THE PER-REGION MAP IS THE DELIVERABLE, NOT THIS SCALAR.** Use [`delta_rho_map`](@ref) and
[`uncertainty_map`](@ref); this accessor exists because the shared `AbstractColocResult` interface
requires a `delta_rho` method, and refusing it outright would break every type-agnostic consumer
for the sake of a number the result genuinely has.
"""
delta_rho(r::SpatialColocResultSpike) = r.meta.derived_rho

# --- The two names the src/ sketch specifies ---------------------------------------------------

"""
    delta_rho_map(r::SpatialColocResultSpike) -> Matrix{Float64}

The PRIMARY deliverable: the `G x G` per-region Delta-rho lattice, in rho units, carrying
`P12_SPATIAL_SENTINEL` at every region that could not be scored (each of which also carries
`r.ood[i] == true` and `r.region_sd[i] === NaN`).

One of the two accessor names the `src/results.jl:180-191` sketch specifies. Returns the field
itself, not a copy and not a transpose.

**Read `r.meta.delta_rho_available` before summarizing this map.** `false` means no control stack
existed, so the map is a REFUSAL rather than a measurement, and an all-sentinel map must never be
reported as a Delta-rho of zero.
"""
delta_rho_map(r::SpatialColocResultSpike) = r.region_delta_rho

"""
    uncertainty_map(r::SpatialColocResultSpike) -> Matrix{Float64}

The per-region posterior sd OF THE DIFFERENCE: the `G x G` uncertainty lattice, `NaN` exactly where
`r.ood` is true.

This is the field `src/amortized/local_map.jl:56-79` says the shipped `LocalColocMap` does not
have ("carries no posterior draws, no Bayes factor and no per-region uncertainty"). Supplying it is
the point of this type. The single-stack uncertainties are `r.region_sd_sample` and
`r.region_sd_control`.
"""
uncertainty_map(r::SpatialColocResultSpike) = r.region_sd

# --- The OOD "was it even checked?" predicate --------------------------------------------------

"""
    p12_ood_checked(r::SpatialColocResultSpike) -> Bool

Whether an OOD operating point was recorded at all, i.e. `r.meta.ood_threshold !== nothing`.

Carries `src/amortized/local_map.jl:88-91`'s rule into this type: **a `false` OOD flag under a
`nothing` threshold means NOT CHECKED, never "in distribution"**. A missing operating point must
leave the flag honestly inert rather than silently default to something arbitrary, and downstream
code must branch on THIS predicate before interpreting an all-`false` OOD map as a clean bill of
health.
"""
p12_ood_checked(r::SpatialColocResultSpike) = r.meta.ood_threshold !== nothing

# --- Field reconstruction from the theta DCT rows ----------------------------------------------

"""
    p12_region_field(coefficients; G = isqrt(length(coefficients))) -> Matrix{Float64}

Reconstruct the `G x G` lattice field from a theta column's DCT coefficient rows.

`coefficients` are in SMOOTHNESS ORDER -- `c0` first, then the deviation coefficients in
`p12_dct_order(G)[2:end]` -- exactly the order `P12_THETA_ROWS` declares. The permutation is
INVERTED here (scattering each coefficient back to its flat mode index) before
`p12_idct_vec` is applied.

THIS IS AN EXACT ORTHONORMAL INVERSE, NOT AN APPROXIMATION. `p12_dct_matrix` is orthonormal
(`C'C == I` to float tolerance), so `p12_idct_vec . p12_dct_vec` is the identity up to rounding and
Parseval holds exactly. That is why D-07's no-lossy-transform rule is NOT engaged by this step: no
information is discarded here. (Truncating the coefficient set to fewer than `G^2` rows WOULD be
lossy -- but that is a choice made upstream, at the theta row budget, and it is visible there.)
"""
function p12_region_field(coefficients::AbstractVector{<:Real};
                          G::Integer = isqrt(length(coefficients)))
    n = length(coefficients)
    G * G == n || throw(ArgumentError(
        "p12_region_field: length(coefficients) = $n is not G^2 for G = $G. Pass the full " *
        "G^2-row coefficient column (c0 plus the G^2-1 deviation rows)."))
    order = p12_dct_order(G)
    c = zeros(Float64, n)
    @inbounds for k in 1:n
        c[order[k]] = float(coefficients[k])   # INVERT the smoothness permutation
    end
    return reshape(p12_idct_vec(c; G = G), G, G)
end

# --- The pre-construction working bundle, and its single sentinel writer ------------------------

"""
    p12_region_maps(grid::Int) -> NamedTuple

Allocate the PRE-CONSTRUCTION working bundle `maps`: the nine mutable `grid x grid` matrices

    (region_delta_rho, region_sd, region_rho_sample, region_sd_sample,
     region_rho_control, region_sd_control, ood_sample, ood_control, ood)

in the constructor's positional order after `grid`. It is named HERE because the bundle is a real
object the write path passes around, and because the marking must happen BEFORE construction: the
constructor validates the sentinel/OOD/`NaN` triple and throws, so a half-marked set can never
legally exist as a `SpatialColocResultSpike`.

**IT IS INITIALISED FULLY UNSCORABLE** -- every rho map at `P12_SPATIAL_SENTINEL`, every sd at
`NaN`, every flag `true`. That is fail-closed: a region is a refusal until something measures it,
so a write path that skips a region yields an honest sentinel rather than a confident zero. The
initial bundle is itself a legal construction (the identity holds vacuously, the union rule holds,
every pairing is exact) -- it is precisely the no-control-stack shape, with the sample stack
unmeasured too.

Use [`p12_mark_unscorable!`](@ref) to mark a region that was ATTEMPTED and failed; write measured
values into the matrices directly, clearing the corresponding flag in the same step.
"""
function p12_region_maps(grid::Integer)
    G = Int(grid)
    G >= 2 || throw(ArgumentError("p12_region_maps: `grid` must be >= 2, got $G."))
    return (region_delta_rho   = fill(P12_SPATIAL_SENTINEL, G, G),
            region_sd          = fill(NaN, G, G),
            region_rho_sample  = fill(P12_SPATIAL_SENTINEL, G, G),
            region_sd_sample   = fill(NaN, G, G),
            region_rho_control = fill(P12_SPATIAL_SENTINEL, G, G),
            region_sd_control  = fill(NaN, G, G),
            ood_sample         = fill(true, G, G),
            ood_control        = fill(true, G, G),
            ood                = fill(true, G, G))
end

"""
    p12_mark_unscorable!(maps, idxs; stack::Symbol) -> maps

Write the sentinel TRIPLE -- `P12_SPATIAL_SENTINEL`, `NaN`, `true` -- for every region index in
`idxs`, in the stack named by `stack` (`:sample`, `:control` or `:both`).

**`stack` HAS NO DEFAULT.** Which stack failed is a fact about the data, not a convention, and a
defaulted answer would be a guess recorded as a measurement. Every call site names it.

**IT WRITES THROUGH THE UNION RULE.** Marking `:sample` sets the sample triple AND the DIFFERENCE
triple -- because a difference is unscorable whenever either side is -- and leaves the control
untouched. `maps` is the bundle of all six value matrices plus all three flag matrices, so ONE call
cannot leave them inconsistent.

Having exactly one writer of the sentinel triple is what keeps the constructor's invariants
satisfiable. A call site that set only one map, or that forgot the union propagation, is exactly
the bug this signature exists to make impossible: the constructor would throw on the half-marked
bundle, but only after the write path had already decided the region was fine.
"""
function p12_mark_unscorable!(maps::NamedTuple, idxs; stack::Symbol)
    stack in (:sample, :control, :both) || throw(ArgumentError(
        "p12_mark_unscorable!: `stack` must be :sample, :control or :both, got :$stack."))
    for i in idxs
        if stack === :sample || stack === :both
            maps.region_rho_sample[i] = P12_SPATIAL_SENTINEL
            maps.region_sd_sample[i]  = NaN
            maps.ood_sample[i]        = true
        end
        if stack === :control || stack === :both
            maps.region_rho_control[i] = P12_SPATIAL_SENTINEL
            maps.region_sd_control[i]  = NaN
            maps.ood_control[i]        = true
        end
        # THE UNION PROPAGATION, ALWAYS -- never conditional on `stack`.
        maps.region_delta_rho[i] = P12_SPATIAL_SENTINEL
        maps.region_sd[i]        = NaN
        maps.ood[i]              = true
    end
    return maps
end

# --- The meta constructor ----------------------------------------------------------------------

"""
    p12_result_meta(; zscore_arm = P12_ZSCORE_ARM, arm, r1_posterior_median, derived_rho,
                      ood_threshold, net_path, delta_rho_available::Bool,
                      schema_version = 1) -> NamedTuple

Build the `meta` a [`SpatialColocResultSpike`](@ref) carries, with exactly
`P12_RESULT_META_KEYS` in that order. Mirrors `p13_calibration_meta` (`spike/p13/result.jl:337`).

`ood_threshold = nothing` MEANS THE OPERATING POINT WAS **NOT CHECKED**, and a `false` OOD flag
under a `nothing` threshold must NEVER be read as "in distribution" -- the rule is
`src/amortized/local_map.jl:88-91`'s, carried into this type verbatim in spirit. A missing
operating point leaves the flag honestly inert rather than silently defaulting to something
arbitrary. Read it through [`p12_ood_checked`](@ref).

`delta_rho_available` HAS NO DEFAULT, DELIBERATELY. It says whether a control stack existed at all,
and every construction must STATE it rather than inherit it. `false` means the single-stack mode:
the Delta-rho map is a REFUSAL, not a measurement. A summary that reports a Delta-rho figure from a
`delta_rho_available = false` result is a defect the constructor cannot catch -- which is exactly
why the flag is mandatory and NAMED rather than inferred from the sentinel pattern.

`zscore_arm` defaults to `P12_ZSCORE_ARM` (R-3): the arm actually in force must be recorded on
every persisted artifact, because the two arms differ in a way no downstream artifact would
otherwise show.
"""
function p12_result_meta(; zscore_arm::Symbol = P12_ZSCORE_ARM,
                         arm,
                         r1_posterior_median,
                         derived_rho,
                         ood_threshold,
                         net_path,
                         delta_rho_available::Bool,
                         schema_version::Integer = 1)
    return (zscore_arm          = zscore_arm,
            arm                 = arm,
            r1_posterior_median = r1_posterior_median,
            derived_rho         = derived_rho,
            ood_threshold       = ood_threshold,
            net_path            = net_path,
            delta_rho_available = delta_rho_available,
            schema_version      = Int(schema_version))
end
