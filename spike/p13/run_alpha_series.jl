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

# spike/p13/run_alpha_series.jl --- REPORTED D-15/D-16 alpha-graded segregation series (NOT a gate).
#
# =============================================================================================
# (a) THIS SERIES CARRIES NO GATE, AND THAT IS A PRE-REGISTERED PROPERTY, NOT A CONCESSION.
# =============================================================================================
# `P13_ALPHA_GATED = false` was frozen in `spike/p13/consts.jl` before any Phase-13 number
# existed. The D-15 discrimination and D-13 calibration gates rest on SIMULATOR GROUND TRUTH --
# the two-factor cut on drawn theta -- and they were run once, and reported, by
# `spike/p13/run_three_way_gate.jl`. This runner scores a CONSTRUCTED ladder whose alpha is a
# knob, not a measured physical quantity, so there is nothing here a pass/fail bar could
# honestly be attached to. This file therefore contains no gate assertion of any kind, and the
# acceptance criteria assert that by source grep so it cannot quietly acquire one later.
#
# A FLAT, NON-MONOTONE OR OTHERWISE DISAPPOINTING LADDER IS A FINDING AND IS REPORTED AS ONE.
# It is not a reason to move a rung, change the transform, retrain, reseed or redraw. The
# pre-registration authorises exactly ONE iteration on exactly ONE pre-declared trigger
# (`P13_ITERATION_TRIGGER`), that trigger is about the D-12 gate and not about this series, and
# plan 13-12 already evaluated it mechanically and recorded that it did NOT fire.
#
# =============================================================================================
# (b) WHAT THIS SERIES IS FOR: the correlation-versus-localization gap.
# =============================================================================================
# The simulator's "exclusion" mechanism IS NEGATIVE INTENSITY CORRELATION across patches (the
# `sign(rho)` flip on channel 2's shared latent component, spike/simulator/forward.jl stage 1).
# A biologist's "mutually exclusive" is DISJOINT SPATIAL LOCALIZATION -- protein B is not where
# protein A is. Those two notions coincide often but not necessarily, and NOTHING ELSE IN THIS
# PHASE PROBES THE DIFFERENCE: the provenance corpus holds exactly one real segregated anchor
# and no graded series at all.
#
# So the ladder is built SPATIALLY (`alpha_segregate` moves ch2 mass out of the ch1 object mask
# and returns it to the complement) and scored by a net trained ONLY on the anti-correlation
# notion. The interesting output is therefore not "is exclusion detected" but AT WHAT ALPHA THE
# INTENSITY-CORRELATION NOTION STARTS TO TRACK DISJOINT LOCALIZATION -- reported here as
# `alpha_star`, the smallest rung at which the mean log BF(exclusion : random) first exceeds 0.
#
# Building the ladder by anti-correlation instead would test the net against its own training
# assumption, which is circular, and would leave the actual question untouched. Displacing the
# ch2 objects instead is likewise rejected: a global displacement IS the mis-registration
# nuisance Phase 11 marginalizes over, so a null would be uninterpretable and a positive result
# unattributable.
#
# =============================================================================================
# (c) THE PHASE-16 DEFERRAL OF THE SINGLE PHYSICAL SEGREGATED ANCHOR.
# =============================================================================================
# The provenance manifest contains exactly ONE real segregated anchor -- `neg-lightmycells-01`,
# physical-primary, negative, segregated, mitochondria. Phase 13 does not read it, and it is
# DEFERRED TO PHASE 16 for three independent reasons, every one of them verified:
#
#   1. its recorded `sha256` is the unfilled sentinel "PENDING-FETCH" -- the bootstrap fetch was
#      deliberately never run, so there is no hash to verify anything against;
#   2. the data directory is EMPTY and carries no tracked bytes -- the file is simply not on
#      disk, and `bytes = 0` on the row records exactly that;
#   3. the row's split is `sealed_holdout`, reachable only through the anti-snooping accessor
#      `open_sealed_holdout(df; reason)`, which is a deliberate access control.
#
# THAT ACCESSOR IS NOT CALLED BY THIS RUNNER -- not directly, not in a test, not behind a flag,
# and no seal-break escape hatch exists anywhere in Phase 13. THE SEAL IS LEFT INTACT ON
# PURPOSE: it exists for Phase 16's BLIND evaluation, and consuming the one segregated anchor
# here would irreversibly burn Phase 16 on the very hypothesis Phase 16 exists to evaluate. The
# n = 1 physical check is therefore deferred, and this file does not reference the corpus
# directory at all -- the paragraph above is provenance prose in a comment, which is the only
# form the ruling permits.
#
# The download-size question that earlier research attached to this deferral is RETIRED
# (assumption A6, withdrawn): Phase 13 fetches nothing, so the size of an unfetched file is moot.
# `bytes = 0` records that nothing was fetched and nothing more.
#
# =============================================================================================
# (d) SCOPE: this run is the `:simulated` arm ONLY, and the two arms are never averaged.
# =============================================================================================
# `P13_ALPHA_SUBSTRATE = (:simulated, :real)` is two-valued and pre-registered. The `:real` arm
# runs the SAME `alpha_segregate` code path over the six committed `test/test_images/` TIFFs and
# is reported separately by plan 13-16's `run_p13_realimage.jl`. THEIR CURVES MUST NEVER BE
# AVERAGED OR MERGED, because `alpha = 0` means different things on the two substrates: here it
# is a GENUINELY RANDOM pair by construction (both acquisitions are drawn at rho_true = 0),
# while on real substrate the unmodified pairs are MODERATELY COLOCALIZED (measured m-bar
# +0.3292 positive / +0.2481 negative). Both are informative; they are not the same experiment.
#
# The CBS `cbs-RG-000` arm is DROPPED FROM SCOPE, not optional: it would need the fetch path,
# and the amended D-15 supplies a real substrate that needs no fetch, no hash gate and no seal.
#
# DECOUPLING (hard constraint, CLAUDE.md, D-01): spike-local. `src/` is reached only READ-ONLY
# and transitively through the simulator/summary contract chain, and step 0 PROVES at run time
# that `src/` and both spike manifests are byte-unchanged. CPU-only; CUDA is never imported; no
# package is installed and the figure reuses the CairoMakie surface already resolved in
# `spike/validation/figures.jl`.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY. Including it does nothing; run it deliberately:
#
#     julia --project=spike spike/p13/run_alpha_series.jl

using JLD2
using Statistics
using Dates
using SHA
import Random123: Philox4x

# ORDER MATTERS. The datagen surface first -- it carries the whole Phase-11 precondition preamble
# (preconditions.jl), the pre-registration (consts.jl), the label surface (labels.jl), the net
# (net.jl), the simulator halves and the src/ summary contract. Then the alpha transform, then
# the figure surface. All guarded, the house guarded-include idiom.
isdefined(@__MODULE__, :p13_sample_theta_given_lambda) || include(joinpath(@__DIR__, "datagen.jl"))
isdefined(@__MODULE__, :verify_alpha_invariants)       || include(joinpath(@__DIR__, "alpha_series.jl"))
isdefined(@__MODULE__, :plot_roc)                      ||
    include(joinpath(@__DIR__, "..", "validation", "figures.jl"))

if !isdefined(@__MODULE__, :P13_ALPHA_REPORT_PATH)
    # NOTE THE PATH DEPTH: from spike/p13/ the repo root is TWO levels up, not three.
    "Repository root, for the run-time decoupling proof (step 0)."
    const P13_ALPHA_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

    "The reported alpha-series artifact. Written BEFORE anything can throw."
    const P13_ALPHA_REPORT_PATH = joinpath(@__DIR__, "alpha_report.jld2")

    "The trained three-way evidence net (plan 13-11). There is no substitute and no fallback."
    const P13_ALPHA_NET_PATH = joinpath(@__DIR__, "three_way_net.jld2")

    "The three-curve ladder figure. A SCRIPT ARTIFACT, never an assertion."
    const P13_ALPHA_FIG_PATH = joinpath(@__DIR__, "..", "figures", "p13_alpha_ladder.png")

    """
        P13_ALPHA_SIM_MAX_VALUE

    The DECLARED dynamic range of the SIMULATED substrate, handed to `alpha_segregate` and to
    `verify_alpha_invariants` as their `max_value` keyword.

    IT IS DELIBERATELY UNBOUNDED, AND THAT IS THE CHOICE-FREE OPTION RATHER THAN A LOOSENED BAR.
    `P13_ALPHA_MAX_VALUE_BOUND = 1.0` is the REAL substrate's range: the committed TIFFs are
    Gray values in `[0, 1]`. The simulator's softplus intensities are unnormalized and routinely
    exceed it (measured realized maxima of 5-10 on this substrate, before any alpha), which is
    exactly why `alpha_segregate` takes the range as a CALLER-SUPPLIED knob instead of deciding
    it internally. Declaring some other finite number here would be inventing a bar this series
    could then be scored against, in a runner whose entire posture is that it scores nothing.
    So the range is declared unbounded and the REALIZED maximum is RECORDED per rung instead --
    which is what the pre-registration actually asks for ("ASSERTED AND RECORDED, NEVER SILENTLY
    CLAMPED"). Nothing is ever clamped, and `max_value_ok` is consequently vacuously true on
    this arm; the report says so rather than quoting it as evidence.
    """
    const P13_ALPHA_SIM_MAX_VALUE = Inf

    """
        P13_ALPHA_RHO

    The true rho BOTH acquisitions of every reported item are drawn at: exactly zero.

    This is what makes the `alpha = 0` rung a GENUINELY RANDOM PAIR BY CONSTRUCTION rather than
    by measurement, which is the property the simulated arm exists to have and the property the
    real arm cannot have. It also lands the item squarely inside the D-05 `RANDOM` class at
    alpha = 0: `|rho_s| < tau` and `|rho_s - rho_c| < tau` both hold with the full tau margin.
    """
    const P13_ALPHA_RHO = 0.0
end

# =============================================================================================
# The RESERVED alpha stream (D-01)
# =============================================================================================

"""
    p13_alpha_rng(idx::Integer) -> Philox4x

The per-item keyed RNG of the reported alpha substrate:
`Philox4x((P13_DEV_SEED xor P13_SALT xor P13_ALPHA_COUNTER, idx))`.

IT IS THE `p13_datagen_rng` CONSTRUCTION AT THE ALPHA COUNTER. The reserved activity counter
goes in the KEY word and the per-item index in the COUNTER word, so an item is a pure function
of its index and the substrate is byte-identical for any thread count and any resumed run.
Because `P13_ALPHA_COUNTER` differs from `P13_DATAGEN_COUNTER` and from `P13_GATE_COUNTER`, the
resulting Philox key differs from both, so this substrate is DISJOINT from the 48,000-item
training pool AND from the 4,000-pair reported evaluation set. The assertions below make that
structural rather than a matter of trust.
"""
p13_alpha_rng(idx::Integer) =
    Philox4x(UInt64, (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_ALPHA_COUNTER),
                      UInt64(idx)))

# The key words that must differ, asserted rather than assumed. "Obviously distinct" is how two
# lanes end up sharing a Philox key.
@assert (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_ALPHA_COUNTER)) !=
        (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_DATAGEN_COUNTER)) "the alpha key word collides with the training-pool key word"
@assert (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_ALPHA_COUNTER)) !=
        (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_GATE_COUNTER)) "the alpha key word collides with the reported gate key word"
@assert (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_ALPHA_COUNTER)) !=
        (UInt64(P13_FIX_SEED) ⊻ P13_SALT ⊻ UInt64(P13_FIXTURE_COUNTER)) "the alpha key word collides with the fixture family"

# =============================================================================================
# Small local helpers
# =============================================================================================

"""
    _p13_alpha_save_report(path; kwargs...) -> String

Atomically persist the reported artifact: write `path * ".tmp"`, REOPEN to integrity-check the
required keys, then `mv(...; force = true)`. The `_save_bf_report` idiom
(`spike/validation/run_bf.jl:59-69`). A crash mid-write leaves a discardable `.tmp`, never a
torn artifact a later reader would happily believe.
"""
function _p13_alpha_save_report(path; kwargs...)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        for k in ("alpha_star", "logbf_exclusion", "logbf_coloc", "mbar", "n_missing", "ladder")
            @assert haskey(f, k) "_p13_alpha_save_report: integrity check failed, $tmp missing $k"
        end
    end
    mv(tmp, path; force = true)
    return path
end

"""
    p13_alpha_blob_sha(path) -> String

The GIT BLOB sha1 of `path`, i.e. exactly what `git hash-object <path>` prints. Recorded so the
artifact points INTO the audit trail: a reader can run `git cat-file -p <sha>` to get the frozen
pre-registration text out of history, or `git log --find-object=<sha>` to get the commit that
introduced it.
"""
function p13_alpha_blob_sha(path)
    bytes = read(path)
    ctx   = SHA.SHA1_CTX()
    SHA.update!(ctx, Vector{UInt8}("blob $(length(bytes))\0"))
    SHA.update!(ctx, bytes)
    return bytes2hex(SHA.digest!(ctx))
end

"""
    p13_alpha_draw_substrate(idx, lambda) -> NamedTuple

Draw ONE reported substrate item: `(pair_s, pair_c, imsize, theta_s, theta_c, lambda, idx)`.

TWO ACQUISITIONS, ONE OF WHICH IS THE CONTROL AND IS NEVER TRANSFORMED. The evidence net reads
a PAIR -- a sample stack and its control stack -- so a ladder over the sample alone is what the
D-05 cut is written in terms of: as alpha grows the SAMPLE moves toward exclusion while the
CONTROL stays random, which is exactly `rho_s < -tau AND rho_s - rho_c < -tau`. Grading BOTH
halves would move the contrast and the level together and the ladder would measure nothing.

BOTH acquisitions are drawn at `P13_ALPHA_RHO = 0.0` and share ONE lambda and ONE image size,
because they are two channels of the SAME acquisition (D-03). The nuisances are drawn from
`p13_sample_theta_given_lambda` -- the SAME conditional draw the training pool used, at the SAME
reference lambda the pair is then encoded at, so the conditioning row never states a level the
column's own shift was not drawn under. The image size comes from the frozen F5 mixture, so the
substrate's joint is the joint the net was trained on except for the one thing deliberately
pinned: rho.
"""
function p13_alpha_draw_substrate(idx::Integer, lambda::Real)
    rng     = p13_alpha_rng(idx)
    theta_s = merge(p13_sample_theta_given_lambda(rng, lambda), (ρ_true = P13_ALPHA_RHO,))
    theta_c = merge(p13_sample_theta_given_lambda(rng, lambda), (ρ_true = P13_ALPHA_RHO,))
    isz     = p13_sample_imsize(rng)
    pair_s  = simulate_pair(rng, theta_s; imsize = isz)
    pair_c  = simulate_pair(rng, theta_c; imsize = isz)
    return (pair_s = pair_s, pair_c = pair_c, imsize = isz,
            theta_s = theta_s, theta_c = theta_c, lambda = Float64(lambda), idx = Int(idx))
end

"""
    p13_alpha_crossing(ladder, curve) -> Union{Float64,Nothing}

`alpha_star`: the SMALLEST ladder value at which `curve` first exceeds zero, or `nothing` when
the ladder never crosses.

`nothing` IS A LEGITIMATE, REPORTABLE OUTCOME AND MUST NOT BE INTERPOLATED AWAY. The crossing
point is reported as a RUNG of the frozen `P13_ALPHA_LADDER`, never as an interpolated value
between rungs: an interpolated crossing would be a number the pre-registration does not contain,
computed from a curve nobody committed to a functional form for. The ladder's resolution IS the
resolution of this statistic, and saying so is part of the result.
"""
function p13_alpha_crossing(ladder, curve)
    for i in eachindex(curve)
        isfinite(curve[i]) && curve[i] > 0.0 && return Float64(ladder[i])
    end
    return nothing
end

"""
    p13_alpha_figure(ladder, lbf_e, lbf_c, mbar, alpha_star; path) -> String

The three-curve ladder figure: `log BF(E:R)(alpha)` and `log BF(C:R)(alpha)` on a shared nats
axis, `m-bar(alpha)` on a second axis, `alpha_star` marked, over one shared alpha axis.

A SCRIPT ARTIFACT AND NEVER AN ASSERTION -- the Phase-4 rule `spike/validation/figures.jl`
states at its head, honoured here by rendering the figure BEFORE the headline rather than after.
Composed in this runner rather than through a `figures.jl` helper for the two reasons plan 13-12
already recorded for its own figure: that file has no ladder helper, and its `_val_fig_path`
hard-routes every output into `spike/validation/figures/`, which is not where this deliverable
lives. The CairoMakie surface, the headless backend activation and the caption discipline are
all inherited from it.
"""
function p13_alpha_figure(ladder, lbf_e, lbf_c, mbar, alpha_star;
                          path = P13_ALPHA_FIG_PATH)
    mkpath(dirname(path))
    a = collect(Float64, ladder)
    star_txt = alpha_star === nothing ? "never crosses 0" : "alpha* = $(alpha_star)"
    fig = Figure(size = (1000, 460))

    ax1 = CairoMakie.Axis(fig[1, 1];
        title    = "alpha-graded segregation, SIMULATED arm -- REPORTED, NOT A GATE",
        subtitle = "mean log Bayes factor per rung over the reported substrate; $star_txt",
        xlabel = "alpha (construction parameter, not a physical quantity)",
        ylabel = "mean log Bayes factor (nats)")
    hlines!(ax1, [0.0]; color = :gray, linestyle = :dash)
    lines!(ax1, a, lbf_e; color = :firebrick, label = "log BF(exclusion : random)")
    scatter!(ax1, a, lbf_e; color = :firebrick)
    lines!(ax1, a, lbf_c; color = :steelblue, label = "log BF(coloc : random)")
    scatter!(ax1, a, lbf_c; color = :steelblue)
    alpha_star === nothing ||
        vlines!(ax1, [alpha_star]; color = :seagreen, linestyle = :dot, linewidth = 2)
    axislegend(ax1; position = :rc)

    ax2 = CairoMakie.Axis(fig[1, 2];
        title    = "summary-level response m-bar(alpha)",
        subtitle = "mask-weighted mean of the PRESENT continuous summary rows",
        xlabel = "alpha", ylabel = "m-bar (mean patch correlation)")
    hlines!(ax2, [0.0]; color = :gray, linestyle = :dash)
    lines!(ax2, a, mbar; color = :darkorange)
    scatter!(ax2, a, mbar; color = :darkorange)
    alpha_star === nothing ||
        vlines!(ax2, [alpha_star]; color = :seagreen, linestyle = :dot, linewidth = 2)

    save(path, fig)
    return path
end

# =============================================================================================
# The reported run
# =============================================================================================

"""
    main(; n_images = P13_ALPHA_N_IMAGES, ...) -> NamedTuple

The reported D-15/D-16 alpha-graded series, end to end and once. There is no gate, so there is
no smoke/reported mode split of the kind the D-12 gate runner needs: every number this function
produces is reported, and none of it is scored against anything.
"""
function main(; n_images::Integer = P13_ALPHA_N_IMAGES,
              net_path = P13_ALPHA_NET_PATH,
              report_path = P13_ALPHA_REPORT_PATH,
              fig_path = P13_ALPHA_FIG_PATH,
              max_value::Real = P13_ALPHA_SIM_MAX_VALUE,
              verbose::Bool = true)

    # --- 0. THE DECOUPLING PROOF, AT RUN TIME (D-01) ---------------------------------------
    # Asserted while the run happens, not only in review: a long reported run that started on a
    # clean tree and finished on a dirty one would be a repudiation hole (T-13-39).
    src_clean = success(Cmd(`git diff --quiet HEAD -- src`; dir = P13_ALPHA_REPO_ROOT))
    env_clean = success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                            dir = P13_ALPHA_REPO_ROOT))
    @assert src_clean "D-01 decoupling breach: `git diff --quiet HEAD -- src` FAILED -- src/ is not byte-unchanged, so this run may not be reported"
    @assert env_clean "D-01 decoupling breach: spike/Project.toml or spike/Manifest.toml is modified -- the spike environment moved under the run"
    verbose && println("STEP 0 (D-01): src/ byte-unchanged = $src_clean   " *
                       "spike Project/Manifest byte-unchanged = $env_clean")

    # --- 1. THE HARD PRECONDITION (D-02) ----------------------------------------------------
    p13_require_phase11()
    basis = load_p13_basis()
    G     = isqrt(length(basis.zt.mean))

    h   = load_three_way(net_path)
    net = h.net
    hlo = h.head_log_odds     # the PERSISTED, MEASURED per-head correction. NEVER re-measured
                              # here: measuring it on this substrate would tune the reported
                              # evidence scale to the very data it is read on (D-07).

    lambda = Float64(P13_PHASE11_REFERENCE_LAMBDA)

    # --- 2. THE LOCKED-VALUE BANNER ---------------------------------------------------------
    if verbose
        println("="^78)
        println("PHASE-13 REPORTED ALPHA SERIES (D-15/D-16) -- SUPPORTING EVIDENCE, NO GATE")
        println("="^78)
        println("SUBSTRATE SCOPE:")
        println("  P13_ALPHA_SUBSTRATE     = $P13_ALPHA_SUBSTRATE")
        println("  THIS RUN COVERS         = :simulated   (alpha = 0 is a GENUINELY RANDOM pair")
        println("                            by construction: both acquisitions drawn at")
        println("                            rho_true = $P13_ALPHA_RHO)")
        println("  THE :real ARM           = plan 13-16's spike/p13/run_p13_realimage.jl, over")
        println("                            the six committed test/test_images/ TIFFs, where")
        println("                            alpha = 0 is MODERATELY COLOCALIZED (m-bar +0.3292")
        println("                            / +0.2481). THE TWO ARMS ARE NEVER AVERAGED.")
        println("LADDER AND TRANSFORM (all READ from the frozen pre-registration):")
        println("  P13_ALPHA_LADDER        = $P13_ALPHA_LADDER")
        println("  P13_ALPHA_BG_QUANTILE   = $P13_ALPHA_BG_QUANTILE   " *
                "(the strictly-positive background floor rule)")
        println("  P13_ALPHA_MASK_RULE     = $P13_ALPHA_MASK_RULE")
        println("  P13_ALPHA_MASK_FRACTION_BOUNDS = $P13_ALPHA_MASK_FRACTION_BOUNDS")
        println("  P13_ALPHA_ZERO_INVARIANT       = $P13_ALPHA_ZERO_INVARIANT")
        println("  P13_ALPHA_INVARIANT_RTOL       = $P13_ALPHA_INVARIANT_RTOL")
        println("  P13_ALPHA_N_IMAGES      = $P13_ALPHA_N_IMAGES        (requested $n_images)")
        println("  declared dynamic range  = $max_value   (the SIMULATED substrate declares an")
        println("                            UNBOUNDED range; the REAL substrate's frozen")
        println("                            P13_ALPHA_MAX_VALUE_BOUND = $P13_ALPHA_MAX_VALUE_BOUND is a Gray-TIFF")
        println("                            property and does not apply here. Realized maxima")
        println("                            are RECORDED per rung, never clamped.)")
        println("NOT A GATE:")
        println("  P13_ALPHA_GATED         = $P13_ALPHA_GATED   " *
                "(no criterion, no bar, no assertion on any curve)")
        println("HYPOTHESIS BOUNDARY (Tier 2, MEASURED not chosen):")
        println("  P13_TAU                 = $P13_TAU  at reference lambda " *
                "$P13_TAU_REFERENCE_LAMBDA")
        println("  encoding lambda         = $lambda   " *
                "(:widest_rung -- the most conservative reading)")
        println("NET UNDER TEST:")
        println("  artifact                = $net_path")
        println("  head_log_odds (MEASURED, persisted) = $hlo")
        println("  consts sha256 at train time         = $(h.consts_sha)")
        println("  tau at train time                   = $(h.tau)   cut = $(h.cut_variant)")
        println("STREAM (D-01):")
        println("  P13_DEV_SEED            = $(repr(UInt64(P13_DEV_SEED))) " *
                "at counter P13_ALPHA_COUNTER = $P13_ALPHA_COUNTER")
        println("  training pool counter   = P13_DATAGEN_COUNTER = $P13_DATAGEN_COUNTER " *
                "(DISJOINT)")
        println("  reported gate counter   = P13_GATE_COUNTER    = $P13_GATE_COUNTER " *
                "(DISJOINT)")
        println("PHYSICAL ANCHOR:")
        println("  P13_PHYSICAL_ANCHOR_DEFERRED_TO = $P13_PHYSICAL_ANCHOR_DEFERRED_TO")
        println("  The one real segregated anchor is triply unavailable (unfilled hash, no")
        println("  bytes on disk, sealed split). The sealed-holdout accessor is NOT called by")
        println("  this runner and the anti-snooping seal is left INTACT for Phase 16's blind")
        println("  evaluation. Nothing under the provenance data tree is checked, fetched or")
        println("  referenced by any line of executable code in this file.")
        println("="^78)
    end

    # --- 3. STREAM AND SCOPE HYGIENE, ASSERTED BEFORE ANYTHING IS DRAWN --------------------
    @assert :simulated in P13_ALPHA_SUBSTRATE "the :simulated substrate is not a member of the pre-registered P13_ALPHA_SUBSTRATE"
    @assert P13_ALPHA_COUNTER != P13_DATAGEN_COUNTER "the reported alpha substrate must not ride the training pool's counter"
    @assert P13_ALPHA_COUNTER != P13_GATE_COUNTER "the reported alpha substrate must not ride the reported gate's counter"
    @assert P13_ALPHA_GATED == false "P13_ALPHA_GATED is not false -- this runner may not become a gate"
    @assert !(UInt64(P13_DEV_SEED) in _p13_forbidden()) "P13_DEV_SEED is a forbidden (pre-observed) seed"
    @assert h.consts_sha == p13_consts_sha() "the frozen pre-registration changed since the net was trained: consts.jl sha256 does not match the artifact's"

    ladder = collect(Float64, P13_ALPHA_LADDER)
    R      = length(ladder)

    # --- 4. THE REPORTED SUBSTRATE ----------------------------------------------------------
    verbose && println("[1/5] drawing $n_images simulated substrate items at " *
                       "P13_ALPHA_COUNTER (rho_true = $P13_ALPHA_RHO, lambda = $lambda) ...")
    t0 = time()
    items = [p13_alpha_draw_substrate(i, lambda) for i in 1:n_images]
    verbose && println("      done in $(round((time() - t0) / 60; digits = 2)) min; " *
                       "image sizes $(sort(unique([it.imsize for it in items]), by = prod))")

    # --- 5. INVARIANTS FIRST, CURVES SECOND -------------------------------------------------
    # A ladder that is flat because pixels were DELETED looks identical to a ladder that is flat
    # because the net cannot see spatial segregation. Only the invariants separate the two, so
    # they are verified on the REPORTED substrate -- not merely on a fixture -- before a single
    # curve is read (T-13-14).
    verbose && println("[2/5] re-verifying the D-16 transform invariants on the REPORTED " *
                       "substrate ...")
    invs = Vector{Any}(undef, n_images)
    all_ok = true
    for i in 1:n_images
        v = verify_alpha_invariants(items[i].pair_s; max_value = max_value)
        invs[i] = v
        ok = v.bitwise_alpha0 && v.no_new_zeros && v.intensity_conserved && v.mask_invariant &&
             v.mask_fraction_ok && v.max_value_ok && v.mbar_monotone && v.missing_nonincreasing
        all_ok &= ok
        verbose && println("      image $(lpad(i, 3))  $(ok ? "PASS" : "FAIL")  " *
                           "mask frac $(round(v.realized_mask_fraction; digits = 4))  " *
                           "src zeros $(v.source_zero_count)  " *
                           "max rel intensity err $(round(v.max_rel_intensity_err; sigdigits = 3))  " *
                           "realized max $(round(v.realized_max_value; sigdigits = 5))")
    end

    # The per-invariant failure counts, so a violation names ITSELF rather than being hunted.
    inv_fields = (:bitwise_alpha0, :no_new_zeros, :intensity_conserved, :mask_invariant,
                  :mask_fraction_ok, :max_value_ok, :mbar_monotone, :missing_nonincreasing)
    inv_fail = Dict(String(f) => count(v -> !getproperty(v, f), invs) for f in inv_fields)

    if !all_ok
        println("\n", "!"^78)
        println("INVARIANT VIOLATION -- THE CURVES ARE NOT READ AND MUST NOT BE READ.")
        println("!"^78)
        for f in inv_fields
            n_bad = inv_fail[String(f)]
            n_bad == 0 || println("  $(rpad(String(f), 24)) FAILED on $n_bad of $n_images images")
        end
        println("")
        println("  WHY THIS STOPS THE RUN INSTEAD OF DEGRADING IT. A FLAT LADDER CAUSED BY")
        println("  DELETED PIXELS IS INDISTINGUISHABLE FROM A NULL RESULT. `_exclude_zero`")
        println("  drops any pixel where either channel is zero, so a transform that pushed")
        println("  masked ch2 pixels to hard zero would REMOVE them from the patch correlation")
        println("  rather than anti-correlate them: the measured correlation would then be over")
        println("  the untouched complement only -- approximately unchanged, i.e. random -- and")
        println("  object-dense patches could fall below the 15-survivor floor and go missing.")
        println("  The ladder would read flat for a reason that has nothing whatsoever to do")
        println("  with segregation, and the phase would conclude the net cannot see spatial")
        println("  exclusion when it was never shown any. Ruling that out is the ONLY thing")
        println("  these invariants are for, so a curve read past a violation would be")
        println("  uninterpretable rather than merely uncertain.")
        println("")
        println("  What HAS been computed is persisted below so the violation itself is on the")
        println("  record; no curve, no crossing point and no figure is produced.")
        println("!"^78, "\n")

        _p13_alpha_save_report(report_path;
            invariants_ok = false, invariant_failures = inv_fail,
            alpha_star = nothing,
            logbf_exclusion = fill(NaN, R), logbf_coloc = fill(NaN, R),
            mbar = fill(NaN, R), n_missing = fill(-1, R), ladder = ladder,
            n_images = n_images, substrate = "simulated",
            invariant_bitwise_alpha0 = [v.bitwise_alpha0 for v in invs],
            invariant_no_new_zeros = [v.no_new_zeros for v in invs],
            invariant_intensity_conserved = [v.intensity_conserved for v in invs],
            invariant_mask_invariant = [v.mask_invariant for v in invs],
            invariant_mask_fraction_ok = [v.mask_fraction_ok for v in invs],
            invariant_max_value_ok = [v.max_value_ok for v in invs],
            invariant_mbar_monotone = [v.mbar_monotone for v in invs],
            invariant_missing_nonincreasing = [v.missing_nonincreasing for v in invs],
            max_rel_intensity_err = [v.max_rel_intensity_err for v in invs],
            source_zero_count = [v.source_zero_count for v in invs],
            realized_mask_fraction = [v.realized_mask_fraction for v in invs],
            realized_max_value = [v.realized_max_value for v in invs],
            master_seed = UInt64(P13_DEV_SEED), salt = UInt64(P13_SALT),
            alpha_counter = Int(P13_ALPHA_COUNTER),
            datagen_counter = Int(P13_DATAGEN_COUNTER),
            gate_counter = Int(P13_GATE_COUNTER),
            consts_sha = p13_consts_sha(),
            consts_git_blob_sha = p13_alpha_blob_sha(joinpath(@__DIR__, "consts.jl")),
            generated = string(Dates.now(Dates.UTC)) * "Z")
        println("persisted the violation record -> $report_path")
        return (invariants_ok = false, invariant_failures = inv_fail,
                alpha_star = nothing, report_path = report_path, figure = nothing)
    end

    verbose && println("      ALL invariants hold on all $n_images reported images; " *
                       "the curves may be read.")

    # --- 6. SCORE EVERY IMAGE AT EVERY RUNG THROUGH THE DOCUMENTED READ SURFACE ------------
    # The CONTROL half is summarized ONCE per item and held fixed across the whole ladder: only
    # the sample is graded, which is what makes the ladder move the D-05 two-factor cut rather
    # than move both of its arguments together.
    verbose && println("[3/5] scoring $(n_images) x $R rungs through three_way_log_bf ...")
    t1 = time()
    LBE  = Matrix{Float64}(undef, R, n_images)   # log BF(exclusion : random)
    LBC  = Matrix{Float64}(undef, R, n_images)   # log BF(coloc : random)
    MB   = Matrix{Float64}(undef, R, n_images)   # m-bar
    NMIS = Matrix{Int}(undef, R, n_images)       # absent-patch count
    MASKROW_STABLE = Vector{Bool}(undef, n_images)   # the summary's mask rows, alpha-invariant?

    for i in 1:n_images
        it = items[i]
        Zc = standardize_summary(encode_d01(patch_summary(build_mci(it.pair_c))),
                                 basis.zt, basis.variant)
        L  = alpha_ladder(it.pair_s; ladder = P13_ALPHA_LADDER, max_value = max_value)
        mask_rows_ref = nothing
        stable = true
        for (r, rung) in enumerate(L.pairs)
            enc = encode_d01(patch_summary(build_mci(rung)))
            nc  = length(enc) ÷ 2
            vals = @view enc[1:nc]
            msk  = @view enc[(nc + 1):(2nc)]
            wsum = sum(msk)
            MB[r, i]   = wsum == 0.0 ? NaN : sum(vals .* msk) / wsum
            NMIS[r, i] = Int(nc - wsum)
            mask_rows_ref === nothing ? (mask_rows_ref = collect(msk)) :
                                        (stable &= collect(msk) == mask_rows_ref)
            Zs = standardize_summary(enc, basis.zt, basis.variant)
            b  = three_way_log_bf(net, p13_encode_pair(Zs, Zc, lambda), hlo)
            LBE[r, i] = b.exclusion
            LBC[r, i] = b.coloc
        end
        MASKROW_STABLE[i] = stable
        verbose && (i % 8 == 0 || i == n_images) &&
            println("      scored $i / $n_images items " *
                    "($(round((time() - t1) / 60; digits = 2)) min elapsed)")
    end

    mean_lbe   = [mean(@view LBE[r, :])   for r in 1:R]
    median_lbe = [median(@view LBE[r, :]) for r in 1:R]
    mean_lbc   = [mean(@view LBC[r, :])   for r in 1:R]
    median_lbc = [median(@view LBC[r, :]) for r in 1:R]
    mean_mbar  = [mean(@view MB[r, :])    for r in 1:R]
    tot_missing = [sum(@view NMIS[r, :])  for r in 1:R]

    # --- 7. THE CROSSING POINT, THE HEADLINE NUMBER -----------------------------------------
    alpha_star = p13_alpha_crossing(ladder, mean_lbe)

    if verbose
        println("-"^78)
        println("THE LADDER (REPORTED; no bar is attached to any number below).")
        println("  alpha    mean logBF(E:R)  median         mean logBF(C:R)  median         " *
                "m-bar      missing")
        for r in 1:R
            println("  $(rpad(ladder[r], 8)) " *
                    "$(lpad(round(mean_lbe[r]; sigdigits = 6), 15)) " *
                    "$(lpad(round(median_lbe[r]; sigdigits = 6), 14)) " *
                    "$(lpad(round(mean_lbc[r]; sigdigits = 6), 16)) " *
                    "$(lpad(round(median_lbc[r]; sigdigits = 6), 14)) " *
                    "$(lpad(round(mean_mbar[r]; sigdigits = 5), 10)) " *
                    "$(lpad(tot_missing[r], 8))")
        end
        println("-"^78)
        if alpha_star === nothing
            println("alpha* : THE LADDER NEVER CROSSES ZERO. The mean log BF(exclusion : random)")
            println("         does not become positive at any rung of the frozen ladder, so the")
            println("         intensity-correlation notion of exclusion does not begin to track")
            println("         disjoint spatial localization anywhere on this ladder. That is the")
            println("         result; it is not interpolated, extended or explained away.")
        else
            println("alpha* = $alpha_star   <-- THE HEADLINE NUMBER.")
            println("         The smallest rung of the frozen ladder at which the mean")
            println("         log BF(exclusion : random) first exceeds 0, i.e. WHERE THE")
            println("         INTENSITY-CORRELATION NOTION OF EXCLUSION BEGINS TO TRACK")
            println("         DISJOINT SPATIAL LOCALIZATION. It is a RUNG, never an")
            println("         interpolated value, so the ladder's own resolution is the")
            println("         resolution of this statistic.")
            println("         alpha IS A CONSTRUCTION PARAMETER, NOT A PHYSICAL QUANTITY, so")
            println("         alpha* CALIBRATES the correlation-versus-localization gap rather")
            println("         than MEASURING it.")
        end
        println("-"^78)
    end

    # --- 8. THE PITFALL-2 WARNING SIGNS, PRINTED EXPLICITLY ---------------------------------
    miss_rose   = tot_missing[end] > tot_missing[1]
    miss_any    = any(tot_missing .> 0)
    maskrows_ok = all(MASKROW_STABLE)
    if verbose
        println("PITFALL-2 WARNING SIGNS (the mechanical explanations a flat ladder could have).")
        println("  absent patches at the FIRST rung (alpha = $(ladder[1])): $(tot_missing[1]) " *
                "of $(n_images * G * G)")
        println("  absent patches at the LAST  rung (alpha = $(ladder[end])): $(tot_missing[end]) " *
                "of $(n_images * G * G)")
        println("  did the absent-patch count RISE across the ladder? " *
                "$(miss_rose ? "YES -- pixels are being DELETED, not moved" : "NO")")
        println("  any absent patch at any rung on any image? $(miss_any ? "YES" : "NO")")
        println("  did the summary's MASK ROWS change with alpha on any image? " *
                "$(maskrows_ok ? "NO -- byte-identical across every rung" : "YES on $(count(!, MASKROW_STABLE)) image(s)")")
        println("  Both signs are the ones the strictly-positive background floor exists to")
        println("  prevent: masked ch2 pixels are pushed DOWN TO a positive floor, never to")
        println("  zero, so the zero SET is untouched and every pixel keeps participating in")
        println("  its patch correlation.")
        println("  NOTE, STATED RATHER THAN QUOTED AS EVIDENCE: because the simulated arm")
        println("  declares an UNBOUNDED dynamic range, the invariant `max_value_ok` is")
        println("  vacuously true on this arm. The realized maxima are reported instead.")
        println("-"^78)
    end

    # --- 9. THE CROSS-SUBSTRATE POINTER AND THE DECLARED SCOPE DROP -------------------------
    if verbose
        println("CROSS-SUBSTRATE POINTER -- READ THIS BEFORE QUOTING ANY NUMBER ABOVE.")
        println("  This run is the :simulated member of the pre-registered two-valued")
        println("  P13_ALPHA_SUBSTRATE. The :real member is reported SEPARATELY by")
        println("  spike/p13/run_p13_realimage.jl (plan 13-16) over the six committed")
        println("  test/test_images/ TIFFs.")
        println("  ONE CODE PATH, TWO EXPERIMENTS. Both arms call the SAME alpha_segregate,")
        println("  typed on AbstractMatrix{Float64} with no substrate branch anywhere, so the")
        println("  rungs are comparable rung for rung. What DIFFERS is what alpha = 0 MEANS:")
        println("  here it is a genuinely random pair by construction (rho_true = " *
                "$P13_ALPHA_RHO on")
        println("  both acquisitions); there it is moderately COLOCALIZED (measured m-bar")
        println("  +0.3292 positive / +0.2481 negative). So the two curves are reported side by")
        println("  side and are NEVER AVERAGED OR MERGED -- they are different experiments,")
        println("  and both are informative.")
        println("  DECLARED SCOPE DECISION: the CBS cbs-RG-000 arm is DROPPED, not optional --")
        println("  it would require the fetch path, and the amended D-15 supplies a real")
        println("  substrate that needs no fetch, no hash gate and no seal. Nothing under the")
        println("  provenance data tree was checked, fetched or referenced by this run.")
        println("-"^78)
    end

    # --- 10. PERSIST ATOMICALLY -------------------------------------------------------------
    verbose && println("[4/5] persisting the reported artifact ...")
    _p13_alpha_save_report(report_path;
        # --- the three reported curves, one entry per ladder rung ---
        ladder = ladder,
        logbf_exclusion = mean_lbe, logbf_exclusion_median = median_lbe,
        logbf_coloc = mean_lbc, logbf_coloc_median = median_lbc,
        mbar = mean_mbar, n_missing = tot_missing,
        alpha_star = alpha_star,
        # --- the per-image, per-rung scores, so no later plan need re-run this ---
        logbf_exclusion_all = LBE, logbf_coloc_all = LBC,
        mbar_all = MB, n_missing_all = NMIS,
        # --- the invariant verdicts and residuals on the REPORTED substrate ---
        invariants_ok = true, invariant_failures = inv_fail,
        invariant_bitwise_alpha0 = [v.bitwise_alpha0 for v in invs],
        invariant_no_new_zeros = [v.no_new_zeros for v in invs],
        invariant_intensity_conserved = [v.intensity_conserved for v in invs],
        invariant_mask_invariant = [v.mask_invariant for v in invs],
        invariant_mask_fraction_ok = [v.mask_fraction_ok for v in invs],
        invariant_max_value_ok = [v.max_value_ok for v in invs],
        invariant_mbar_monotone = [v.mbar_monotone for v in invs],
        invariant_missing_nonincreasing = [v.missing_nonincreasing for v in invs],
        max_rel_intensity_err = [v.max_rel_intensity_err for v in invs],
        source_zero_count = [v.source_zero_count for v in invs],
        realized_mask_fraction = [v.realized_mask_fraction for v in invs],
        realized_max_value = [v.realized_max_value for v in invs],
        # --- the Pitfall-2 warning signs ---
        missing_rose = miss_rose, missing_any = miss_any,
        mask_rows_stable = maskrows_ok, mask_rows_stable_per_image = MASKROW_STABLE,
        # --- the realized substrate ---
        substrate = "simulated", n_images = n_images,
        imsize = [it.imsize for it in items],
        rho_true = P13_ALPHA_RHO, lambda = lambda, grid = G,
        # --- EVERY pre-registered value this run READ, written INTO the artifact ---
        P13_ALPHA_LADDER = collect(Float64, P13_ALPHA_LADDER),
        P13_ALPHA_BG_QUANTILE = P13_ALPHA_BG_QUANTILE,
        P13_ALPHA_MASK_RULE = P13_ALPHA_MASK_RULE,
        P13_ALPHA_MASK_FRACTION_BOUNDS = collect(Float64, P13_ALPHA_MASK_FRACTION_BOUNDS),
        P13_ALPHA_ZERO_INVARIANT = P13_ALPHA_ZERO_INVARIANT,
        P13_ALPHA_INVARIANT_RTOL = P13_ALPHA_INVARIANT_RTOL,
        P13_ALPHA_N_IMAGES = P13_ALPHA_N_IMAGES,
        P13_ALPHA_GATED = P13_ALPHA_GATED,
        P13_ALPHA_SUBSTRATE = collect(String.(P13_ALPHA_SUBSTRATE)),
        P13_ALPHA_MAX_VALUE_BOUND = P13_ALPHA_MAX_VALUE_BOUND,
        declared_max_value = Float64(max_value),
        P13_TAU = P13_TAU, P13_TAU_REFERENCE_LAMBDA = P13_TAU_REFERENCE_LAMBDA,
        P13_PHYSICAL_ANCHOR_DEFERRED_TO = P13_PHYSICAL_ANCHOR_DEFERRED_TO,
        # --- provenance: the stream, the net, the frozen pre-registration ---
        master_seed = UInt64(P13_DEV_SEED), salt = UInt64(P13_SALT),
        alpha_counter = Int(P13_ALPHA_COUNTER),
        datagen_counter = Int(P13_DATAGEN_COUNTER),
        gate_counter = Int(P13_GATE_COUNTER),
        net_path = net_path, net_head_log_odds = hlo,
        net_consts_sha = h.consts_sha, net_meta = h.meta,
        consts_sha = p13_consts_sha(),
        consts_git_blob_sha = p13_alpha_blob_sha(joinpath(@__DIR__, "consts.jl")),
        phase11_net = P13_PHASE11_NET, phase11_sha = p13_phase11_sha(),
        note = "REPORTED, NOT GATED (P13_ALPHA_GATED = false). The :simulated arm of the " *
               "two-valued P13_ALPHA_SUBSTRATE; the :real arm is plan 13-16's and the two " *
               "are never averaged. alpha is a construction parameter, not a physical " *
               "quantity, so alpha_star calibrates the correlation-versus-localization gap " *
               "rather than measuring it.",
        generated = string(Dates.now(Dates.UTC)) * "Z")
    verbose && println("      persisted reported artifact -> $report_path")

    # --- 11. THE FIGURE. A SCRIPT ARTIFACT, NEVER AN ASSERTION ------------------------------
    verbose && println("[5/5] rendering the three-curve ladder figure ...")
    figout = p13_alpha_figure(ladder, mean_lbe, mean_lbc, mean_mbar, alpha_star;
                              path = fig_path)
    verbose && println("      figure -> $figout")

    return (ladder = ladder, alpha_star = alpha_star,
            logbf_exclusion = mean_lbe, logbf_exclusion_median = median_lbe,
            logbf_coloc = mean_lbc, logbf_coloc_median = median_lbc,
            mbar = mean_mbar, n_missing = tot_missing,
            invariants_ok = true, invariant_failures = inv_fail,
            missing_rose = miss_rose, mask_rows_stable = maskrows_ok,
            n_images = n_images, substrate = "simulated",
            report_path = report_path, figure = figout)
end

# --- Script entry point ----------------------------------------------------------------------
# GUARDED so `include`-ing this file can NEVER trigger the reported run. Run it deliberately:
#
#     julia --project=spike spike/p13/run_alpha_series.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
