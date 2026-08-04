#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/p15_misspec.jl --- the THREE Phase-15 sweep axes, and the seven-axis family table.
#
# `.planning/phases/15-calibration-operating-envelope-and-ci-gate/15-CONTEXT.md` names spillover,
# PSF, autofluorescence and registration as the misspecification axes the calibration operating
# envelope is measured over. `test/gate/misspec.jl:181-184` ships four families —
# `(texture, noise, optics, background)` — so `spillover` and `registration` have no generator at
# all, and `background` is not the autofluorescence axis it looks like (see D-02a below). This file
# supplies the three missing generators and the seven-axis table, and CHANGES NOTHING ELSE.
#
# ================== THIS FILE ADDS AXES AND MUTATES NOTHING — WHY THAT MATTERS ==================
# `OOD_FAMILIES` is an INPUT TO A RECORDED VERDICT. `gate_ood_roc` builds `pooled_pos` from the
# strongest-level scores over ALL families, and `ood_gate` returns
# `passed = auc_ok && all(values(fam_pass))` (`misspec.jl:316-396`). Appending three families to
# that table would therefore change the pooled AUC and add three terms to the pass conjunction —
# so the shipped, recorded `OOD — PASS (pooled AUC 1.0)` result would stop being reproducible from
# the code that claims it. Not "would need re-running": would silently become a DIFFERENT number
# reported under the old one's name.
#
# So the seven axes ride a LOCAL `P15_FAMILIES = merge(OOD_FAMILIES, ...)`, passed through
# `gate_ood_roc`'s EXISTING `families =` keyword. `misspec.jl` stays byte-unchanged and its
# LF-normalized SHA-256 is pinned in `p15_consts.jl`'s `P15_FROZEN_FILE_SHA256.misspec`, asserted
# on every test run by `test/gate/test_p15_families.jl`. There is no assignment to `OOD_FAMILIES`
# anywhere below, and there must never be one.
#
# ============================ THE MECHANISM TYPING IS A DELIVERABLE =============================
# Each axis carries a mechanism label from `P15_AXIS_MECHANISM`, and the label decides how a break
# READS. An IN-PRIOR axis that breaks coverage is a TRAINING failure: the net saw those θ during
# training and still miscalibrates. An OUT-OF-MODEL axis that breaks is a SCOPE limit — expected,
# honest, and exactly what the OOD flag exists for. A domain map that reported one number for both
# would be unactionable. The three axes here are the only IN-PRIOR ones, which is precisely why
# they had to be built: without them the envelope would have no in-prior half at all.
#
# RUNG 0 IS AN IDENTITY, NOT AN APPROXIMATION. Every generator here delegates to
# `ProteinCoLoc.simulate_pair` and at `level = 0` passes θ through untouched, so rung 0 is
# byte-identical to the unmisspecified forward model on an identically-seeded stream. That identity
# is what makes rung 0 usable as the break anchor (D-03/D-04a) and it is asserted, not asserted-to.
#
# STANDALONE-LOADABLE, in the established gate-file style: the guarded includes below let this file
# be loaded on its own, and are idempotent under `run_gate.jl` and `runtests.jl`.

using ProteinCoLoc

isdefined(@__MODULE__, :SBC_M)               || include(joinpath(@__DIR__, "p15_consts.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "harness.jl"))
isdefined(@__MODULE__, :OOD_FAMILIES)        || include(joinpath(@__DIR__, "misspec.jl"))

# ============================================================================
# The three IN-PRIOR sweep axes (D-02, D-02a, D-03)
# ============================================================================
#
# All three are PURE θ-OVERRIDES on the unmodified forward model: they change one field of θ to a
# frozen ladder magnitude and delegate. Nothing about `simulate_pair` is reimplemented here — a
# forked copy of the forward model would make the rung-0 identity unverifiable by construction,
# which is the one property these axes exist to have.
#
# WHAT THE SHARED STREAM DOES AND DOES NOT BUY. `simulate_pair` consumes rng in stage 1 (three
# latent Gaussian fields), stage 2 (two Bernoulli thinning masks) and stage 7 (Poisson + Gaussian
# noise). Stages 4 (spillover mix), 5 (autofluorescence offset) and 6 (sub-pixel shift) consume
# NONE. So at a fixed seed, rung k and rung 0 share θ*, both latent triples and both masks — the
# matched-pairs design, and it is asserted as an exact image equality in `test_p15_families.jl`.
# It does NOT extend past the end of the call: stage 7's Poisson consumption is RATE-dependent and
# the override is a change of rate, so the stream state after the call differs by rung. MEASURED,
# not assumed (same testset): `registration` and `autofluorescence` diverge from rung 0 at all five
# rungs for every fixture θ tried, while `spillover` — which perturbs channel 1 only — is
# θ-dependent, diverging at 0, 1, 2 or 5 of 5 rungs depending on the draw. Since `sbc_gate` threads
# ONE rng through all M draws (`harness.jl:188-204`), two rungs agree on draw 1 and then see
# different θ* sequences from draw 2 on. Cross-rung ECEs are therefore neither independent
# replicates nor exactly paired ones; a report must claim neither.

"""
    misspec_spillover(rng, θ; imsize = SBC_IMSIZE, level = 0, G = 8)

Axis SPILLOVER — mechanism **IN-PRIOR EXTRAPOLATION**.

`level = 0` is the pass-through rung: byte-identical to `ProteinCoLoc.simulate_pair(rng, θ)` on an
identically-seeded stream. `level = k` in `1:length(P15_SPILLOVER_LADDER)` overrides `θ.spillover`
with `P15_SPILLOVER_LADDER[k]` and changes nothing else.

`P15_SPILLOVER_LADDER[2]` **is** `maximum(ProteinCoLoc.SPILLOVER_PRIOR)` (= 0.2), so rung 2 is the
MARKED PRIOR BOUNDARY rung (D-03): a break at or below it is inside the training distribution and
reads as a training failure, a break above it reads as extrapolation. Rung 1 sits strictly inside
the prior so the two cases are distinguishable at all. The ladder's top rung is the SIMULATOR'S OWN
limit — `simulate_pair` throws above `spillover = 1.0` (`simulator.jl:232`) — not a chosen stopping
point.

Stage 4 (`ch1 .+= θ.spillover .* ch2`) consumes no rng, so every rung shares the latent fields and
thinning masks of rung 0 at the same seed (see the block comment above for where that stops).
"""
function misspec_spillover(rng, θ; imsize = SBC_IMSIZE, level::Integer = 0, G::Integer = 8)
    # Guarded FIRST and UNCONDITIONALLY, including at level 0: an imsize the patch grid cannot
    # summarize is an error at every rung, and finding that out only at rung 1 would waste the
    # anchor draw.
    _guard_misspec_imsize(imsize, G)
    level == 0 && return ProteinCoLoc.simulate_pair(rng, θ; imsize = imsize)
    return ProteinCoLoc.simulate_pair(rng, merge(θ, (; spillover = P15_SPILLOVER_LADDER[level]));
                                      imsize = imsize)
end

"""
    misspec_registration(rng, θ; imsize = SBC_IMSIZE, level = 0, G = 8)

Axis REGISTRATION — mechanism **IN-PRIOR EXTRAPOLATION**.

`level = 0` is the pass-through rung. `level = k` overrides BOTH `θ.shift_dx` and `θ.shift_dy` with
`P15_SHIFT_LADDER[k] / √2` — a fixed 45° direction chosen so the ladder is defined on the shift
NORM rather than on a per-axis component, which is the quantity `maximum(ProteinCoLoc.SHIFT_PRIOR)`
is expressed in. `P15_SHIFT_LADDER[2] = 1.0 px` is that prior boundary, so rung 2 is the MARKED
BOUNDARY rung (D-03).

The norm reproduces the ladder value to within ONE ULP, not exactly: `r/√2` followed by `hypot` is
not an exact Float64 round trip (0.5 → 0.49999999999999994). Sub-ulp on a sub-pixel shift is
physically nothing, but it is stated rather than claimed exact.

HONEST NOTE ON THE UPPER RUNGS: `simulate_pair`'s stage-6 `warp` fills the vacated border with
`ProteinCoLoc.BG_FLOOR` (`simulator.jl:290`). At rung 5 (8 px on each axis) that fill is a visible
frame around the image. **The border fill is part of the misspecification, not an artifact to be
corrected** — a real misregistered acquisition has no data there either, and cropping it away would
measure a different, easier problem than the one the envelope is about. The ladder stops at 8 px
because beyond it the fill, rather than the registration error, would dominate what is measured.

Stage 6 consumes no rng, so the matched-pairs property of the block comment applies here too.
`simulate_pair` imposes no range guard on the shift beyond finiteness, so this ladder's endpoint is
a stated choice, unlike spillover's.
"""
function misspec_registration(rng, θ; imsize = SBC_IMSIZE, level::Integer = 0, G::Integer = 8)
    _guard_misspec_imsize(imsize, G)
    level == 0 && return ProteinCoLoc.simulate_pair(rng, θ; imsize = imsize)
    d = P15_SHIFT_LADDER[level] / sqrt(2.0)          # 45°: |(d, d)| == the ladder value (to 1 ulp)
    return ProteinCoLoc.simulate_pair(rng, merge(θ, (; shift_dx = d, shift_dy = d));
                                      imsize = imsize)
end

"""
    misspec_autofluorescence(rng, θ; imsize = SBC_IMSIZE, level = 0, G = 8)

Axis AUTOFLUORESCENCE — mechanism **IN-PRIOR EXTRAPOLATION**.

`level = 0` is the pass-through rung. `level = k` overrides `θ.autofluorescence` with
`P15_AF_LADDER[k]`. `P15_AF_LADDER[2]` **is** `maximum(ProteinCoLoc.AUTOFLUORESCENCE_PRIOR)`
(= 0.1), the MARKED PRIOR BOUNDARY rung (D-03).

WHY THIS AXIS EXISTS AT ALL (D-02a, in one sentence): this is a PURE ADDITIVE OFFSET axis, and it
is distinct from the shipped `background` family, which composes an illumination gradient × a
radial vignette × a bleed the simulator has NO parameter for at any θ and whose level-1 bleed is
already 0.4 — four times the prior maximum before its ladder even starts. `background` is therefore
reclassified OUT-OF-MODEL and carries no boundary rung; adding this seventh, pure-offset axis is
what makes D-03's marked boundary rung attainable on every in-prior axis without an exception.

Stage 5 adds the offset rather than sampling it, so no rng is consumed there and the matched-pairs
property of the block comment applies. `simulate_pair` guards only `autofluorescence ≥ 0`
(`simulator.jl:234`), so this ladder's endpoint, like registration's, is a stated choice.
"""
function misspec_autofluorescence(rng, θ; imsize = SBC_IMSIZE, level::Integer = 0, G::Integer = 8)
    _guard_misspec_imsize(imsize, G)
    level == 0 && return ProteinCoLoc.simulate_pair(rng, θ; imsize = imsize)
    return ProteinCoLoc.simulate_pair(rng,
                                      merge(θ, (; autofluorescence = P15_AF_LADDER[level]));
                                      imsize = imsize)
end

# ============================================================================
# The seven-axis table (D-02, D-02a)
# ============================================================================

"""
    P15_FAMILIES :: NamedTuple

The SEVEN Phase-15 sweep axes: the four shipped `OOD_FAMILIES` entries FIRST, in their fixed
Phase-5 order so the existing report ordering is preserved, then the three in-prior axes above.
Key order matches `P15_AXES`.

This is a `merge` — a NEW NamedTuple — precisely so that `OOD_FAMILIES` is not mutated; see the
header. Pass it through `gate_ood_roc`'s existing `families =` keyword; never append to the frozen
table.
"""
const P15_FAMILIES = merge(OOD_FAMILIES, (spillover        = misspec_spillover,
                                          registration     = misspec_registration,
                                          autofluorescence = misspec_autofluorescence))

"""
    p15_axis_mechanism(axis) -> Symbol

`:in_prior` or `:out_of_model` for a Phase-15 axis. The label is a DELIVERABLE of the domain map,
not bookkeeping: it is what turns "this axis broke at rung 3" into either "the net miscalibrates on
data it was trained for" or "this is outside the forward model's scope, which is what the OOD flag
is for".
"""
p15_axis_mechanism(axis) = P15_AXIS_MECHANISM[axis]

"""
    p15_axis_boundary_rung(axis) -> Union{Int,Nothing}

The rung at which an axis's ladder sits exactly at the simulator prior's boundary, or `nothing`.

`nothing` means "THIS AXIS HAS NO PRIOR SUPPORT TO MARK" — the out-of-model axes are not
expressible at any θ, so there is no boundary to place. It never means "unmeasured".
"""
p15_axis_boundary_rung(axis) = P15_PRIOR_BOUNDARY_RUNG[axis]

# NOTE FOR THE SWEEP RUNNER (15-05): there is deliberately NO level-0 entry point for the four
# SHIPPED families here. Rung 0 is the unmisspecified model, so it is the SAME measurement for
# every axis; it is measured ONCE per net through `default_simulator()` and shared. Giving each
# shipped family its own level-0 path would spend seven anchor arms to compute one number, and
# would invite seven slightly different anchors to be compared against each other.
