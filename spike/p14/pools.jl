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

# spike/p14/pools.jl --- the shared draw and null-fitting machinery, with TWO SILENT FAILURE
# MODES closed by executable checks.
#
# THIS FILE EXISTS SO THE CONSTRUCTION IS WRITTEN ONCE. Five Phase-14 runners need a draw and a
# density null. Copied into each of them, the two properties below would drift, and both of them
# fail SILENTLY: the code runs, the quantile computes, the number looks fine.
#
# =================================================================================================
# PITFALL 1 -- THE CALIBRATION SET MUST BE pi-DISTRIBUTED, NOT CLASS-BALANCED.
# =================================================================================================
# The Phase-13 gate set is balanced to equal thirds (realized E 1334 / R 1333 / C 1333) while pi at
# the inherited cut is E 0.3119 / R 0.4713 / C 0.2168. SPLIT CONFORMAL'S GUARANTEE REQUIRES THE
# CALIBRATION POINTS TO BE EXCHANGEABLE WITH THE TEST POINT, and a class-balanced calibration set
# against a pi-distributed deployment stream is not exchangeable: the realized coverage then misses
# 1 - alpha by an amount driven by the class-conditional score distributions. The Bayesian-FDR rule
# is hit twice over, because a false-discovery proportion depends on PREVALENCE and a balanced set
# has had its prevalence replaced.
#
# The forbidden route is the Phase-13 accept/reject balancer at spike/p13/datagen.jl:583 -- named
# HERE, IN A COMMENT, and nowhere in the executable text of this file, because the guard that
# forbids it greps the comment-stripped source and would otherwise flag the file that forbids it.
# Every draw below goes through `p13_simulate_indices`, the unstratified path.
#
# =================================================================================================
# PITFALL 4 -- THE SHIPPED OOD NULL IS ON A DIFFERENT BASIS THAN THE PHASE-13 NET.
# =================================================================================================
# The shipped bundle's standardizer is the Phase-7 one. The Phase-13 net rides PHASE 11's, inherited
# and never re-fit (`meta.zt_provenance` = "INHERITED from the Phase-11 research NPE, never re-fit
# (D-02)"). Scoring a Phase-11-standardized summary against a Phase-7-fitted Mahalanobis null
# produces a number with NO INTERPRETATION -- not a wrong number, an uninterpretable one. So the
# null here is fitted on the Phase-13 net's OWN pool, `meta.pool_dir`, and the fact is recorded on
# every return as `basis_provenance = :p13_pool_phase11_zt` rather than left to be inferred.
#
# The shipped comparison bundle lives at artifacts/amended_v2/grid_8/ood_nulls_8.jld2. That path is
# written HERE, IN A COMMENT, for the same reason as above, and this file does not read it: if a
# runner ever does, it is a COMPARISON REFERENCE and the two thresholds are NEVER fused.
#
# =================================================================================================
# THE POOL-OVERLAP HAZARD (.planning/CONVENTIONS.md C-04) AND WHY IT DOES NOT APPLY.
# =================================================================================================
# Index-keyed generation means a "separate scoring pool" at the SAME configuration is a SUPERSET of
# the training pool, never a held-out set. Phase 14's draws therefore ride a DIFFERENT Philox key
# family (P14_DEV_SEED / P14_FIX_SEED xor P14_SALT), and that is ASSERTED rather than assumed:
# `p14_datagen_key_word` proves its key word is a member of the Phase-14 cross-product
# `_p14_key_words()` and a non-member of the whole Phase-13 family `_p14_p13_key_words()`, both of
# which spike/test/test_p14_consts.jl already exhausts. A NEW SEED PATH THAT BYPASSES THAT PROOF IS
# THE FAILURE THIS FUNCTION EXISTS TO PREVENT.
#
# READ-ONLY, AND THE 51 MB POOL IS WHY. `spike/data/cache/p13` is gitignored bulk data that this
# file READS and never writes; this milestone has already destroyed a 54 MB gitignored pool once.
# Phase-14 caches live under a NEW root, `p14_cache_root()`. Nothing here writes anything at all:
# no mkpath, no save, no artifact.
#
# DECOUPLING (hard constraint, CLAUDE.md and D-01): spike-local. Reaches spike/p13/ and
# spike/validation/ READ-ONLY through guarded includes, touches no src/ byte, adds no dependency.
#
# NOTE THE PATH DEPTH: from spike/p14/ the repo root is TWO levels up, not three.
#
# EVERY CALLER USES THE GUARDED INCLUDE:
#     isdefined(@__MODULE__, :P14_POOLS_LOADED) || include(joinpath(@__DIR__, "pools.jl"))

# --- Guarded includes, in dependency order ------------------------------------------------------
# The Phase-14 Tier-1 pre-registration: the seeds, the salt, the reserved counters, the sample sizes.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
# The D-07 threshold loader -- the ONLY sanctioned route to the inherited three-way cut.
isdefined(@__MODULE__, :P14_PROVENANCE_LOADED) || include(joinpath(@__DIR__, "provenance.jl"))
# The three-class posterior surface (and, transitively, the Phase-13 net and result types).
isdefined(@__MODULE__, :p14_class_posterior) || include(joinpath(@__DIR__, "posterior.jl"))
# The Phase-13 generator: `p13_simulate_indices`, `p13_draw_labelled_item`, `load_pool`.
isdefined(@__MODULE__, :load_pool) || include(joinpath(@__DIR__, "..", "p13", "datagen.jl"))
# The Phase-13 label surface: `class_masses`, and the EXCLUSION / RANDOM / COLOC enum.
isdefined(@__MODULE__, :three_way_label) || include(joinpath(@__DIR__, "..", "p13", "labels.jl"))
# The frozen Phase-11 basis handle, memoized so every Phase-13/14 file shares ONE `zt` object.
isdefined(@__MODULE__, :load_p13_basis) ||
    include(joinpath(@__DIR__, "..", "p13", "preconditions.jl"))
# The in-repo OOD machinery, REUSED rather than re-implemented: `fit_ood_nulls`, `maha_score`,
# `id_threshold`, `OOD_ID_QUANTILE`, `OOD_FAMILIES`, `OOD_GRID_LEVELS`.
isdefined(@__MODULE__, :fit_ood_nulls) ||
    include(joinpath(@__DIR__, "..", "validation", "ood.jl"))

import Random123: Philox4x     # counter-based RNG; re-imported so this file is self-sufficient

# THE GUARD COVERS ONLY THE `const`s. Julia 1.12 DROPS every docstring written inside an
# `if ... end` block (the parser emits the `Core.@doc` call but the docsystem never registers it),
# so the documented functions sit at top level; method redefinition under a re-include is silent
# and harmless, while `const` redefinition is what would warn. Same split, same reason, as
# spike/p14/provenance.jl and spike/p14/posterior.jl.
if !isdefined(@__MODULE__, :P14_POOLS_LOADED)
    # The three class names, in the order `class_masses` returns them. Every read in this file is
    # BY NAME through this tuple; there is no positional class index anywhere below, for the same
    # reason posterior.jl has none -- five class orderings coexist in this codebase and an
    # inconsistent relabelling produces numbers that look right.
    const P14_MASS_KEYS = (:exclusion, :random, :coloc)

    # Memoized so a runner that calls the loader in a loop pays for the artifact reads ONCE. The
    # loader itself is unchanged and still asserts the full four-way agreement on first call.
    const _P14_TAU_CACHE  = Ref{Any}(nothing)
    const _P14_META_CACHE = Ref{Any}(nothing)

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else can
    # make this body skip.
    const P14_POOLS_LOADED = true
end

# --- The inherited numbers, loaded once ---------------------------------------------------------

"""
    p14_tau() -> Float64

The inherited three-way cut, memoized. Obtained ONLY through `p14_load_tau()` (D-07), which asserts
the four-way agreement between the net, the probe report, the gate report and the Phase-13
pre-registration before returning. No threshold literal appears in this file under any name.
"""
function p14_tau()
    _P14_TAU_CACHE[] === nothing && (_P14_TAU_CACHE[] = p14_load_tau().tau)
    return _P14_TAU_CACHE[]::Float64
end

"""
    p14_net_meta() -> NamedTuple

The Phase-13 net's persisted metadata, memoized: the source of `pi_class_masses` (the MEASURED
class prior), `pool_dir` (the pool the OOD null is fitted on) and `zt_provenance`.

READ, NEVER GUESSED. `pool_dir` in particular is read off the trained net rather than re-derived
through the cache layer, so nothing here can create, invalidate or regenerate a pool.
"""
function p14_net_meta(; net_path = joinpath(P14_P13_DIR, "three_way_net.jld2"))
    if _P14_META_CACHE[] === nothing
        isfile(net_path) || error("p14_net_meta: no three-way net artifact at $net_path")
        _P14_META_CACHE[] = load_three_way(net_path).meta
    end
    return _P14_META_CACHE[]
end

"""
    p14_pi_class_masses() -> NamedTuple

`meta.pi_class_masses`, the MEASURED class prior: the realized class fractions of 200,000 prior
draws under the inherited cut, `(exclusion = 0.31272, random = 0.47073, coloc = 0.21655)`. Read off
the net's own metadata and re-keyed by nothing -- it is returned exactly as stored, and every
consumer below reads it BY NAME.
"""
p14_pi_class_masses() = p14_net_meta().pi_class_masses

# --- The Phase-14 per-item stream ---------------------------------------------------------------

"""
    p14_datagen_key_word(counter; rng_for = p14_rng) -> UInt64

The FIRST Philox key word of a Phase-14 per-item draw stream: `seed ⊻ P14_SALT ⊻ counter`, with the
seed taken from the family `rng_for` names (`p14_rng` = reported, `p14_fix_rng` = fixture).

# Why the counter goes in the KEY word and the index in the second word

Mirrors `p13_datagen_rng` (`spike/p13/datagen.jl:180`) exactly. `p14_rng(counter)` puts the activity
counter in the second word, which is right for an activity consuming ONE sequential stream. A POOL
cannot do that: the second word must carry the per-item GLOBAL INDEX, because that is what makes a
column a pure function of its index and therefore the draw byte-identical for any thread count. So
the reserved counter is folded into the key word instead.

# The disjointness is PROVED here, not assumed

Two assertions, and they are the whole point of this function existing rather than a lambda at each
call site. `spike/p14/consts.jl` proves disjointness on the FULL Philox key-word cross-product --
`seed ⊻ salt ⊻ counter` against every Phase-13 key word, not merely seed inequality. This function
asserts that the key word it is about to hand out is (a) a MEMBER of that proved Phase-14 set and
(b) a NON-MEMBER of the whole Phase-13 family. A new seed path that bypassed the cross-product would
therefore fail HERE, on the first draw, instead of silently sharing a stream with a Phase-13 pool.

That matters because of the pool-overlap hazard (`.planning/CONVENTIONS.md` C-04): index-keyed
generation makes a same-configuration pool a SUPERSET of the training pool, never a held-out set. A
different key family is the only thing that makes a Phase-14 draw genuinely fresh.
"""
function p14_datagen_key_word(counter::Integer; rng_for = p14_rng)
    base = rng_for === p14_rng     ? UInt64(P14_DEV_SEED) ⊻ P14_SALT :
           rng_for === p14_fix_rng ? UInt64(P14_FIX_SEED) ⊻ P14_SALT :
           throw(ArgumentError(
               "p14_datagen_key_word: `rng_for` must be `p14_rng` (reported) or `p14_fix_rng` " *
               "(fixture). A foreign RNG family would key a stream whose disjointness from the " *
               "Phase-13 lane was never proved by the key-word cross-product in " *
               "spike/p14/consts.jl, which is exactly the failure this argument is restricted to " *
               "prevent."))
    counter in P14_ALL_COUNTERS || throw(ArgumentError(
        "p14_datagen_key_word: counter $counter is not one of the reserved Phase-14 counters " *
        "$(P14_ALL_COUNTERS). Only reserved counters were carried through the pre-registration's " *
        "key-word cross-product, so an ad-hoc counter is an unproved stream."))
    kw = base ⊻ UInt64(counter)
    kw in Set(_p14_key_words()) || error(
        "p14_datagen_key_word: the key word 0x$(string(kw, base = 16, pad = 16)) is NOT a member " *
        "of the Phase-14 key-word cross-product asserted in spike/p14/consts.jl. The stream would " *
        "then be one whose disjointness was never proved.")
    kw in Set(_p14_p13_key_words()) && error(
        "p14_datagen_key_word: the key word 0x$(string(kw, base = 16, pad = 16)) collides with a " *
        "Phase-13 key word. A Phase-14 evaluation would then ride draws a Phase-13 number already " *
        "rode -- silently, and undetectably from the artifacts alone.")
    return kw
end

"""
    p14_datagen_rng(counter; rng_for = p14_rng) -> Function

The per-item keyed RNG CONSTRUCTOR of a Phase-14 draw: `idx -> Philox4x((key_word, idx))`, with the
key word from [`p14_datagen_key_word`](@ref).

Returned as a function of the index because that is the calling convention
`p13_simulate_indices(indices, basis; rng_for = ...)` expects (`spike/p13/datagen.jl:502`, through
`p13_draw_labelled_item`'s `rng = rng_for(idx)`). Passing a constructor that IGNORED the index would
hand every item the same stream and produce `n` identical items -- which is silent, because
identical items are perfectly valid items.
"""
function p14_datagen_rng(counter::Integer; rng_for = p14_rng)
    kw = p14_datagen_key_word(counter; rng_for = rng_for)
    return idx -> Philox4x(UInt64, (kw, UInt64(idx)))
end

# --- The draw -----------------------------------------------------------------------------------

"""
    p14_draw_classes(n; counter, rng_for = p14_rng, tau = p14_tau()) -> Vector{ThreeWayClass}

The CLASS LABELS of the first `n` items of the very draw [`p14_draw_pool`](@ref) would produce, at a
cost of two prior draws each and NO simulation.

Same stream, same index map, same label rule: the label is a function of theta alone
(`p13_draw_labelled_item`, `spike/p13/datagen.jl:308`), so it can be read before any image exists.
That is not an approximation of the pool's classes, it is the same computation the pool performs --
and `p14_draw_pool` asserts the agreement item by item through `p13_assert_classes_agree`, so the
two cannot drift.

It exists because the pi-distribution check wants thousands of labels while an image draw costs
about a second each. Never returns summaries: a caller that needs `Z` needs the pool.
"""
function p14_draw_classes(n::Integer; counter::Integer, rng_for = p14_rng, tau::Real = p14_tau())
    n >= 1 || throw(ArgumentError("p14_draw_classes: n must be >= 1 (got $n)"))
    item_rng = p14_datagen_rng(counter; rng_for = rng_for)
    return [p13_draw_labelled_item(i; tau = tau, rng_for = item_rng).class for i in 1:n]
end

"""
    p14_draw_pool(n; counter, basis = load_p13_basis(), rng_for = p14_rng, tau = p14_tau(),
                  parallel = Threads.nthreads() > 1) -> Vector{<:NamedTuple}

Draw `n` labelled items UNSTRATIFIED on the Phase-14 stream at `counter`. Each item carries
`(Zs, Zc, lambda, theta, class)` plus the Phase-13 fields `(rho_s, rho_c, imsize, idx)`, where
`theta = (sample = ..., control = ...)` is the pair of drawn parameter vectors and `class` is the
three-way label of `(rho_s, rho_c)` under the INHERITED cut (D-07, obtained through `p14_tau()`).

# THE DRAW IS UNSTRATIFIED, AND THAT IS A REQUIREMENT, NOT A DEFAULT

Every reported Phase-14 draw is pi-distributed. The reason is not aesthetic:

 1. **Split conformal's guarantee requires the calibration points to be EXCHANGEABLE with the test
    point.** A class-balanced calibration set scored against a pi-distributed deployment stream is
    not exchangeable, and the realized coverage then misses `1 - alpha` by an amount driven by the
    class-conditional score distributions. Nothing throws; the quantile computes; the number looks
    fine.
 2. **A false-discovery proportion depends on PREVALENCE.** Balancing the classes replaces the very
    prevalence the Bayesian-FDR rule is controlling against.

So the draw goes through `p13_simulate_indices`, the unstratified path, and the balancing helper is
never called -- asserted by a comment-stripped source grep in `spike/test/test_p14_pools.jl`, and
the realized masses are checked by [`p14_assert_unstratified`](@ref) and persisted into every
artifact rather than assumed.

# The stream, and why CONVENTIONS C-04 does not bite

The draw rides `P14_DEV_SEED ⊻ P14_SALT ⊻ counter` (or the fixture seed), a DIFFERENT Philox key
family from every Phase-13 stream. That matters because index-keyed generation means a "separate
scoring pool" at the same configuration is a SUPERSET of the training pool, never a held-out set
(`.planning/CONVENTIONS.md` C-04). The freshness is ASSERTED, not assumed: `p14_datagen_key_word`
checks membership in the Phase-14 cross-product and non-membership in the Phase-13 family, both of
which `spike/test/test_p14_consts.jl` exhausts over every seed/salt/counter combination.

# Determinism

Each item is a pure function of its index, so `(n, counter, rng_for)` reproduces the same items
bitwise on any thread count and on any resume; `parallel` changes the wall clock and nothing else.
"""
function p14_draw_pool(n::Integer; counter::Integer, basis = load_p13_basis(),
                       rng_for = p14_rng, tau::Real = p14_tau(),
                       parallel::Bool = Threads.nthreads() > 1)
    n >= 1 || throw(ArgumentError("p14_draw_pool: n must be >= 1 (got $n)"))
    item_rng = p14_datagen_rng(counter; rng_for = rng_for)
    idxs  = collect(1:Int(n))
    items = p13_simulate_indices(idxs, basis; tau = tau, rng_for = item_rng, parallel = parallel)

    # The cheap theta-only half, re-derived at the SAME indices on the SAME stream, so the returned
    # item can carry its theta. `p13_draw_labelled_item` is deterministic in the index, so this is
    # the same draw the simulation performed -- and the agreement is asserted rather than trusted,
    # which is the `p13_assert_classes_agree` contract (spike/p13/datagen.jl:530).
    draws = [p13_draw_labelled_item(i; tau = tau, rng_for = item_rng) for i in idxs]
    p13_assert_classes_agree(items, [d.class for d in draws], idxs)

    return [merge(items[j], (theta = (sample = draws[j].theta_s, control = draws[j].theta_c),))
            for j in eachindex(idxs)]
end

# --- The equal-thirds detector ------------------------------------------------------------------

"""
    p14_assert_unstratified(classes; n = length(classes), tol_sd = 3.0) -> NamedTuple

Assert that a realized class vector is pi-DISTRIBUTED and is NOT class-balanced, and return the
realized masses so the caller can persist them into the artifact.

Two separate checks, in this order and for two different reasons:

 1. **The balance check, first, because it is the more specific diagnosis.** If all three realized
    masses sit within `tol_sd` binomial standard errors of `1/3` SIMULTANEOUSLY, this throws and
    names the pitfall. A balanced set would also fail check 2, but it would fail it as "a mass is
    off pi", which is the symptom rather than the cause.
 2. **The pi-agreement check.** Each class must be within `tol_sd` binomial standard errors
    `sqrt(p(1-p)/n)` of its value in `meta.pi_class_masses`. This catches the balanced set, a
    truncated set, a filtered set, and a draw that silently rode the wrong stream.

The two are not redundant. At the sample sizes this phase draws, the EXCLUSION mass alone is within
3 standard errors of `1/3` -- so a check for "some mass is near a third" would fire on a perfectly
good pi-distributed draw. Only the SIMULTANEOUS condition distinguishes a balanced set, and only
check 2 catches a set that is neither balanced nor pi.

`n` is separable from `length(classes)` so a caller can state the sample size the standard error
should be computed at when it summarizes a larger draw. Every mass is read BY NAME.
"""
function p14_assert_unstratified(classes; n::Integer = length(classes), tol_sd::Real = 3.0)
    n >= 1 || throw(ArgumentError("p14_assert_unstratified: n must be >= 1 (got $n)"))
    realized = class_masses(classes)
    expected = p14_pi_class_masses()

    third    = 1.0 / 3.0
    se_third = sqrt(third * (1.0 - third) / n)
    balanced = all(k -> abs(getproperty(realized, k) - third) <= tol_sd * se_third, P14_MASS_KEYS)
    balanced && error("""
        p14_assert_unstratified: THE DRAW IS CLASS-BALANCED, NOT pi-DISTRIBUTED.

        realized: $(realized)
        expected: $(expected)   (meta.pi_class_masses, the MEASURED class prior)
        all three masses are within $(tol_sd) binomial SE ($(round(tol_sd * se_third; digits = 4)))
        of 1/3 at n = $n.

        This is the phase's first silent-failure mode. Split conformal's guarantee requires the
        calibration points to be EXCHANGEABLE with the test point, and a class-balanced set against
        a pi-distributed deployment stream is not exchangeable -- the realized coverage then misses
        1 - alpha by an amount driven by the class-conditional score distributions. A Bayesian false
        discovery proportion depends on PREVALENCE, which balancing replaces outright. Nothing
        throws downstream: the code runs, the quantile computes, the number looks plausible.

        The Phase-13 gate set IS balanced by construction and is DEVELOPMENT SUBSTRATE ONLY -- see
        `p14_gate_dev_substrate`. Reported Phase-14 draws come from `p14_draw_pool`.
        """)

    for k in P14_MASS_KEYS
        p  = Float64(getproperty(expected, k))
        r  = Float64(getproperty(realized, k))
        se = sqrt(p * (1.0 - p) / n)
        abs(r - p) <= tol_sd * se || error("""
            p14_assert_unstratified: the realized `$k` mass is $(round(r; digits = 5)), which is
            $(round(abs(r - p) / se; digits = 2)) binomial SE from its expected $(round(p; digits = 5))
            (bar: $(tol_sd) SE, se = $(round(se; digits = 5)) at n = $n).

            realized: $(realized)
            expected: $(expected)

            The draw is not the prior's class mix. Either it was filtered, truncated or balanced, or
            it rode a stream or a cut it was not supposed to ride.
            """)
    end
    return realized
end

# --- The Phase-13 gate set: development substrate, never a reported number ------------------------

"""
    p14_gate_dev_substrate(; path = ...) -> NamedTuple

Load the Phase-13 three-way gate report and return it LABELLED
`(role = :development_and_crosscheck_only, reported = false, ...)`.

# It may never be the source of a reported Phase-14 number, for TWO independent reasons

 1. **It is class-balanced** (realized E 1334 / R 1333 / C 1333), so every argument in
    [`p14_assert_unstratified`](@ref) applies to it directly.
 2. **It is the very set the Phase-13 ECE gate was reported on.** Calibrating a Phase-14 hedge on it
    would calibrate the hedge on the set whose calibration was the headline claim -- the coverage
    number would then be partly a restatement of a result already banked, and no reader could tell
    which part.

# What it IS good for, and this is not a grudging concession

Getting the code right, sanity plots, fixture construction, and a REPORTED CROSS-CHECK beside a
number derived from a fresh draw. A cross-check that agrees is informative and costs nothing; what
is forbidden is sourcing the number itself.

`role` and `reported` are FIELDS rather than documentation, so a runner cannot consume this set by
accident and silently: the labels travel with the data into whatever record embeds it.
"""
function p14_gate_dev_substrate(; path = joinpath(P14_P13_DIR, "three_way_gate_report.jld2"))
    isfile(path) || error("p14_gate_dev_substrate: no Phase-13 gate report at $path")
    report = JLD2.load(path)
    return (report = report, path = path,
            role = :development_and_crosscheck_only,
            reported = false,
            note = "The Phase-13 three-way gate set is CLASS-BALANCED and is the set the Phase-13 " *
                   "ECE gate was itself reported on. Development, fixtures and reported " *
                   "CROSS-CHECKS only; never the source of a reported Phase-14 number.")
end

# --- Cache roots ----------------------------------------------------------------------------------

# THE PHASE-13 CACHE ROOT IS READ-ONLY. Phase 14 reads the ~51 MB gitignored pool under it (the OOD
# null fit below) and writes NOTHING there -- no shard, no meta, no temporary file. A datagen run
# that reached for the same root would overwrite it silently, because a shard writer's job is to
# write shards, and this milestone has already destroyed a 54 MB gitignored pool once.
#
# The root separation mirrors the p11/p12 separation asserted at
# spike/test/test_p12_decoupling.jl:321-349, and is asserted again for this phase at
# spike/test/test_p14_decoupling.jl's "Phase-14 caches use a NEW root" testset.

"""
    p14_cache_root() -> String

The Phase-14 cache root, `spike/data/cache/p14`. A PATH, not a directory: this function creates
nothing, and no function in this file creates anything. The Phase-13 root beside it is read-only for
the whole phase (see the comment block above this definition).
"""
p14_cache_root() = normpath(joinpath(@__DIR__, "..", "data", "cache", "p14"))

# --- The OOD reference: the density null on the PHASE-13 net's own basis ---------------------------

"""
    p14_ood_reference(; pool_dir = nothing, q = OOD_ID_QUANTILE, basis = load_p13_basis())
        -> NamedTuple

Fit the summary-density Mahalanobis null on the Phase-13 net's OWN training pool and take the
operating point at the pre-registered in-distribution quantile.

Returns `(available, thr, nulls, ood_nulls, n_pool, n_shards, pool_dir, q, basis_provenance,
channels_wired, channels_not_wired, score_quantiles, note)`.

# The construction, and what it REUSES

The same construction as `p13_real_id_reference` (`spike/p13/run_p13_realimage.jl:412`): resolve
`pool_dir` from the trained net's `meta.pool_dir` (**never guessed, never re-derived through the
cache layer**, so this function cannot create, invalidate or regenerate a pool), read every
`shard_*.jld2`, stack BOTH the `Zs` and `Zc` halves -- both are simulator acquisitions standardized
through the SAME frozen basis, so their union is precisely "the acquisitions this net was trained
on" -- fit the null, and take `id_threshold(scores; q = OOD_ID_QUANTILE)`.

It is expressed through the in-repo `fit_ood_nulls` / `maha_score` / `id_threshold`
(`spike/validation/ood.jl`) rather than through a local copy, so there is one implementation of the
covariance and one of the score. The operating point is the PRE-REGISTERED in-distribution quantile
and never the post-hoc maximum of `TPR - FPR` along a ROC, which is a labelled-data quantity and is
therefore not an operating point a detector can be shipped with.

# PITFALL 4, which is the reason this function exists at all

The shipped bundle's null was fitted on the PHASE-7 standardizer, while the Phase-13 net rides
PHASE 11's, inherited and never re-fit (`meta.zt_provenance` = "INHERITED from the Phase-11 research
NPE, never re-fit (D-02)"). Scoring a Phase-11-standardized summary against the Phase-7-fitted null
produces a number with NO INTERPRETATION. So the null is fitted HERE, on this net's own pool, and
`basis_provenance = :p13_pool_phase11_zt` is recorded on every successful return. If the shipped
bundle is ever read by a runner it is a COMPARISON REFERENCE, labelled as such, and the two
thresholds are NEVER fused into one.

# A NAMED LIMIT FOR THE REPORT: only ONE channel is wired

`channels_wired = (:density,)` and `channels_not_wired = (:pp, :noise)`. The shipped flag is an
OR-fusion over three channels; this lane wires the density channel only. A `:clear` resolved
downstream therefore means "one channel says clear", not "three channels say clear". Recorded as a
field so it travels into the artifact instead of being papered over.

# It NEVER throws and NEVER substitutes a threshold

Every failure path returns `available = false`, `thr = nothing` and a `note` explaining why -- the
`bad(note)` idiom, copied from `p13_real_id_reference`. A missing operating point must leave the
flag HONESTLY INERT rather than default to an arbitrary value: `p14_ood_state` then reads
`:not_checked` and `p14_fuse` ABSTAINS (D-06). Substituting a plausible number would convert a
loud, correct abstention into a silent, confident answer.
"""
function p14_ood_reference(; pool_dir = nothing, q::Real = OOD_ID_QUANTILE,
                             basis = load_p13_basis())
    bad(note, pd) = (available = false, thr = nothing, nulls = nothing,
                     ood_nulls = NamedTuple(), n_pool = 0, n_shards = 0,
                     pool_dir = pd, q = q,
                     basis_provenance = :unavailable,
                     channels_wired = (:density,), channels_not_wired = (:pp, :noise),
                     score_quantiles = nothing, note = note)

    pd = pool_dir
    if pd === nothing
        # RESOLVED FROM THE TRAINED NET, never guessed. A net that records no pool is a net whose
        # in-distribution reference cannot be reconstructed, which is a reportable fact, not an
        # excuse to pick a directory.
        try
            pd = p14_net_meta().pool_dir
        catch err
            return bad("the Phase-13 net metadata could not be read, so no pool directory could " *
                       "be resolved: $(sprint(showerror, err))", nothing)
        end
    end
    (pd isa AbstractString && isdir(pd)) ||
        return bad("no readable pool directory at $(repr(pd)); the in-distribution reference " *
                   "could not be built, so no operating point exists", pd)

    shards = sort(filter(p -> startswith(basename(p), "shard_") && endswith(p, ".jld2"),
                         readdir(pd; join = true)))
    isempty(shards) && return bad("the pool directory holds no shard file", pd)

    cols = Any[]
    n_sh = 0
    try
        for s in shards
            p = load_pool(s)              # READ-ONLY, schema-checked by the Phase-13 reader
            push!(cols, p.Zs)
            push!(cols, p.Zc)
            n_sh += 1
        end
    catch err
        return bad("a pool shard could not be read: $(sprint(showerror, err))", pd)
    end

    Z      = reduce(hcat, cols)
    nulls  = fit_ood_nulls(Z; variant = basis.variant)
    scores = [maha_score(nulls, Float64.(@view Z[:, j])) for j in axes(Z, 2)]
    thr    = id_threshold(scores; q = q)
    isfinite(thr) ||
        return bad("the fitted operating point is not finite ($thr); a non-finite threshold is " *
                   "NOT CHECKED under the shipped semantics, so it is reported as unavailable " *
                   "rather than carried", pd)

    return (available = true,
            thr = thr,
            nulls = nulls,
            # The shape `p14_ood_state` and the shipped verdict both read: the fit PLUS the
            # operating point. `:thr` is present ONLY on this path, which is what makes the
            # three-valued resolver able to tell a fitted null from an unfitted one.
            ood_nulls = merge(nulls, (thr = thr,)),
            n_pool = size(Z, 2),
            n_shards = n_sh,
            pool_dir = pd,
            q = q,
            basis_provenance = :p13_pool_phase11_zt,
            channels_wired = (:density,),
            channels_not_wired = (:pp, :noise),
            score_quantiles = (q50 = quantile(scores, 0.50), q90 = quantile(scores, 0.90),
                               q95 = quantile(scores, 0.95), q99 = quantile(scores, 0.99),
                               qmax = maximum(scores)),
            note = "Density null fit on the Phase-13 labelled training pool (both acquisitions of " *
                   "every pair), standardized through the net's OWN inherited Phase-11 basis. " *
                   "Operating point at the pre-registered in-distribution quantile. The noise and " *
                   "posterior-predictive channels are NOT wired in this lane.")
end

"""
    p14_pool_provenance(ref) -> NamedTuple

The pool block every Phase-14 runner embeds beside its result, built from what
[`p14_ood_reference`](@ref) returned: which directory, how many shards, how many bytes, how many
columns, at which quantile, on which basis.

`pool_listing_sha256` digests the shard NAMES AND SIZES, not their contents, and is named so that it
cannot be mistaken for a content hash. Hashing 51 MB on every run to detect a change that a schema
check and a column count already catch would be a cost paid for a stronger-sounding word; a listing
digest catches the failures that actually happen here -- a shard added, removed, truncated or
regenerated at a different size.
"""
function p14_pool_provenance(ref::NamedTuple)
    pd = ref.pool_dir
    listing = String[]
    bytes = 0
    if pd isa AbstractString && isdir(pd)
        for f in sort(filter(p -> startswith(p, "shard_") && endswith(p, ".jld2"), readdir(pd)))
            sz = filesize(joinpath(pd, f))
            bytes += sz
            push!(listing, "$f:$sz")
        end
    end
    return (pool_dir = pd,
            pool_n_shards = ref.n_shards,
            pool_bytes = bytes,
            pool_listing_sha256 = bytes2hex(SHA.sha256(join(listing, "\n"))),
            pool_n_columns = ref.n_pool,
            ood_available = ref.available,
            ood_thr = ref.thr,
            ood_id_quantile = ref.q,
            ood_basis_provenance = ref.basis_provenance,
            ood_channels_wired = ref.channels_wired,
            ood_channels_not_wired = ref.channels_not_wired)
end

# --- The misspecification arm ----------------------------------------------------------------------

"""
    p14_misspec_pool(n; counter = P14_OOD_COUNTER, family, level = :strongest,
                     basis = load_p13_basis(), rng_for = p14_rng, tau = p14_tau()) -> NamedTuple

Draw a MATCHED pair of arms: `n` in-distribution items and `n` misspecified items generated from the
SAME drawn parameters, at the same image sizes, on the same stream.

`family` is a key of the in-repo `OOD_FAMILIES` (`:texture`, `:noise`, `:optics`, `:background`);
`level = :strongest` resolves to the top rung of the existing grid. **No new misspecification is
authored here.** These four families are the pre-registered positive controls, each one producing
structure the shared-latent smooth-field simulator cannot produce, and inventing a fifth in this
phase would be a new experiment wearing an old experiment's name.

# Why the arms are MATCHED item by item

Both arms are built from one `p13_draw_labelled_item(idx)`: same lambda, same theta pair, same image
size, same RNG position. The only difference between the two arms is which generator produced the
pixels. An unmatched comparison would confound the misspecification with the parameter draw, and an
abstention-rate margin between two differently-parameterized arms measures partly the parameters.

# The summariser is asserted to be the canonical one

`_p14_summarize` composes the same four calls the Phase-13 generator composes
(`spike/p13/datagen.jl:363-368`). Because it is a second spelling of one composition, the first item
of every call re-derives the in-distribution summary THROUGH IT and asserts equality with the
canonical `p13_simulate_labelled` output. If the two ever diverge, this throws on item one rather
than reporting an OOD margin that is partly a summarization difference.

Returns `(id, ood, family, level, n, counter, note)` with both arms of equal length.
"""
function p14_misspec_pool(n::Integer; counter::Integer = P14_OOD_COUNTER,
                          family::Symbol, level = :strongest,
                          basis = load_p13_basis(), rng_for = p14_rng, tau::Real = p14_tau())
    n >= 1 || throw(ArgumentError("p14_misspec_pool: n must be >= 1 (got $n)"))
    haskey(OOD_FAMILIES, family) || throw(ArgumentError(
        "p14_misspec_pool: family must be one of $(keys(OOD_FAMILIES)) -- the pre-registered " *
        "positive controls. Authoring a new misspecification here would be a new experiment " *
        "under an existing name. Got $family."))
    lvl = level === :strongest ? OOD_GRID_LEVELS : Int(level)
    (1 <= lvl <= OOD_GRID_LEVELS) || throw(ArgumentError(
        "p14_misspec_pool: level must be :strongest or an integer in 1:$(OOD_GRID_LEVELS); " *
        "got $level"))

    gen      = OOD_FAMILIES[family]
    item_rng = p14_datagen_rng(counter; rng_for = rng_for)

    # THE SUMMARISER AGREEMENT CHECK, ON ITEM ONE, BEFORE ANY ARM IS BUILT. See the docstring.
    let d = p13_draw_labelled_item(1; tau = tau, rng_for = item_rng),
        canonical = p13_simulate_labelled(1, basis; tau = tau, rng_for = item_rng)
        mine = _p14_summarize(simulate_pair(d.rng, d.theta_s; imsize = d.imsize), basis)
        @assert mine == canonical.Zs "p14_misspec_pool: the local summariser disagrees with the Phase-13 generator on the in-distribution half of item 1. The two arms would then differ by a summarization step as well as by the misspecification, and the reported margin would be partly an artifact of this file."
    end

    id_arm  = Vector{Any}(undef, Int(n))
    ood_arm = Vector{Any}(undef, Int(n))
    for j in 1:Int(n)
        canonical = p13_simulate_labelled(j, basis; tau = tau, rng_for = item_rng)
        id_arm[j] = canonical
        d = p13_draw_labelled_item(j; tau = tau, rng_for = item_rng)
        Zs = _p14_summarize(gen(d.rng, d.theta_s; imsize = d.imsize, level = lvl), basis)
        Zc = _p14_summarize(gen(d.rng, d.theta_c; imsize = d.imsize, level = lvl), basis)
        ood_arm[j] = (Zs = Zs, Zc = Zc, lambda = d.lambda, class = d.class,
                      imsize = d.imsize, idx = d.idx)
    end

    return (id = [x for x in id_arm], ood = [x for x in ood_arm],
            family = family, level = lvl, n = Int(n), counter = Int(counter),
            note = "Matched arms: identical theta, lambda and image size in both, differing only " *
                   "in which generator produced the pixels. Misspecification families are the " *
                   "pre-registered in-repo positive controls at the strongest grid rung.")
end

"""
    _p14_summarize(pair, basis) -> Vector

One image pair to one standardized summary, through the SAME four calls the Phase-13 generator uses
(`spike/p13/datagen.jl:363-368`): `standardize_summary(encode_d01(patch_summary(build_mci(pair))),
basis.zt, basis.variant)`.

It exists only so that both arms of [`p14_misspec_pool`](@ref) pass through ONE code path -- the
misspecified arm cannot call `p13_simulate_labelled`, because that function owns its own generator.
It is a second spelling of an existing composition and is therefore asserted equal to it on the
first item of every call, which is the check that keeps "second spelling" from becoming "second
implementation".
"""
_p14_summarize(pair, basis) =
    standardize_summary(encode_d01(patch_summary(build_mci(pair))), basis.zt, basis.variant)
