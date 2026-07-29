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

# spike/p13/preconditions.jl --- the Phase-11 hard precondition and the p13 include preamble (D-02, D-03).
#
# =============================================================================================
# PHASE 13 TRAINS ON PHASE 11's REGISTRATION-AWARE **FROZEN** `zt`, NOT THE SHIPPED GRID-8 `zt`.
# THERE IS NO FALLBACK, ANYWHERE, IN ANY PHASE-13 FILE.
# A MISSING PHASE-11 NET IS A **BLOCK**, NOT A SKIP: THIS FILE THROWS.
# THIS IS THE SINGLE FILE THAT BINDS PHASE-11 NAMES. EVERY LATER `spike/p13/*.jl` INCLUDES THIS
# FILE AND IMPORTS THE PHASE-11 NAMES FROM HERE, RATHER THAN REACHING INTO PHASE 11 DIRECTLY.
# =============================================================================================
#
# WHY A WHOLE FILE FOR A PRECONDITION. The failure this guards against is SILENT. A Phase-13
# evidence net trained on the SHIPPED grid-8 standardizer would build, train, converge and score
# exactly as if nothing were wrong -- and it would mean something different from what the phase
# claims. D-02 makes the Phase-13 net train on Phase 11's registration-aware input surface so the
# two-way-versus-three-way comparison is like-for-like; D-03 makes the registration-uncertainty
# level a conditioning input read in EXACTLY the encoding the Phase-11 net was trained with. Both
# of those are properties of WHICH ARTIFACT WAS LOADED, and nothing downstream can recover them
# from a trained net. So the check has to happen at the load boundary, once, loudly.
#
# ONE LAMBDA, NOT TWO (pre-registered modelling choice; do not reopen). The sample stack and the
# control stack are two channels of the SAME acquisition under the SAME stated registration
# uncertainty. A second lambda would encode the belief that the user holds DIFFERENT beliefs about
# the alignment of the two stacks, which is not a claim the read surface expresses and not a claim
# any acquisition supports. `p13_encode_pair` therefore takes ONE lambda for the pair.
#
# THE CONDITIONING LENGTH IS DERIVED, NEVER WRITTEN. `n_cond = length(encode_lambda(lambda))`,
# computed once from the Phase-11 encoder itself. The Phase-13 input width follows as
# `three_way_input_dim(G, p13_conditioning_length())`. Writing the realized width as a literal
# would be silently wrong the moment Phase 11 delivered a conditioning vector of another length,
# and nothing would throw.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. This file reaches `src/` only
# transitively and READ-ONLY through the simulator/contract chain below; it writes nothing, it
# imports no `ProteinCoLoc` symbol, and it adds no package. `src/` stays byte-unchanged.
#
# =============================================================================================
# PHASE-11 STATUS RULING (recorded here on purpose -- read this before you conclude the
# precondition binds an artifact it should not)
# =============================================================================================
# A reader running `git log` will find that **Phase 11 is CLOSED (2026-07-27), NOT COMPLETE**
# (`.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-CLOSURE.md`), that
# `11-07-SUMMARY.md` carries `status: blocked`, and that plans 11-08..11-11 are SUPERSEDED. Taken
# literally, this plan's Task-1 item "confirm Phase 11 is COMPLETE" FAILS. Phase 13 binds the
# Phase-11 net anyway, and the reasoning is written here rather than in a planning document so it
# can stand or fall in the open:
#
#   1. THE GATE THAT FAILED WAS MIS-SPECIFIED, NOT FAILED. `P11_LAMBDA_ABLATION_FACTOR = 2.502`
#      was derived as half a ratio of the REGISTRATION-INDUCED COMPONENT of the Delta-rho
#      posterior spread, but the tripwire applied it to the TOTAL posterior SD. Phase 11's own
#      independent coverage evidence puts the true width ratio at ~1.00 (RMSE 0.13334 at
#      lambda_min versus 0.13338 at lambda_max, ratio 1.0003), against measured net ratios of
#      1.109 / 0.987 / 0.936. 11-CLOSURE states it plainly: the correct answer is ~1.00, the net
#      was right, and the gate failed it anyway. NO THRESHOLD WAS EDITED to reach that reading.
#   2. THE CONDITIONING INPUT IS ALIVE AND MEASURED. 11-07-SUMMARY records shift marginals whose
#      SDs track lambda with ratios 4.84-9.05 against an ideal prior ratio of 12.0, plus a fourth
#      decisive check that the 129th input row is read and used. What does NOT track lambda is the
#      rho_true posterior WIDTH -- and that is the finding, not a defect in the artifact.
#   3. WHAT PHASE 13 ACTUALLY NEEDS IS THE SIX-ITEM ARTIFACT CONTRACT ASSERTED BELOW: a trained
#      net, a frozen summary transform, a frozen theta transform, the lambda encoder plus its
#      range, the imsize and prior provenance, and the theta arity and field order. All six are
#      present and are asserted here at runtime. The blocking checkpoint that preceded this file
#      existed to stop the artifact NAMES being GUESSED; nothing below is guessed.
#   4. PHASE 13 IS ALREADY DESIGNED TO ABSORB THIS HONESTLY. Plan 13-12 carries the must_have
#      "(D-03) the log-BF pair's response to the registration-uncertainty level is measured and
#      reported honestly, INCLUDING A FLAT RESPONSE AS A FINDING". Phase 11 having established
#      that the width is genuinely flat in lambda makes a flat Phase-13 lambda response a
#      PRE-AUTHORIZED, EXPECTED outcome rather than something to be explained away. Proceeding
#      therefore does not launder a negative result.
#   5. NOTE PRECISELY WHAT D-02 CLAIMS. It claims Phase 13 trains on Phase 11's registration-aware
#      INPUT SURFACE rather than the shipped grid-8 basis, so the comparison is like-for-like. It
#      does NOT claim registration is inferable from the 8x8 summary. Phase 11's negative result
#      leaves D-02 intact.
#
# NOTHING IN THIS FILE ASSERTS "PHASE 11 IS COMPLETE", and nothing in this file weakens the
# no-fallback rule. A missing net is still a hard `error`.
#
# =============================================================================================
# THREE DEVIATIONS FROM THE 13-09 PLAN SKELETON, ALL FORCED BY THE ARTIFACT AS BUILT
# =============================================================================================
# (a) `p11_architecture.jl` MUST BE IN SCOPE BEFORE `load_npe`. `P11BoundedThetaTransform` is
#     defined there. Load the artifact without it and JLD2 emits a mere WARNING -- "type
#     Main.P11BoundedThetaTransform{...} does not exist in workspace; reconstructing" -- and hands
#     back a `JLD2.ReconstructedMutable` PLACEHOLDER. The load appears to succeed, `h.theta-zt`
#     has the right field names, and the frozen theta transform is silently not the real one.
#     That is exactly the class of silent-wrong-basis failure this file exists to prevent, so the
#     include preamble orders `p11_architecture.jl` first AND `p13_require_phase11` asserts the
#     CONCRETE TYPE, so a future regression is loud rather than silent.
# (b) THE HANDLE CARRIES NEITHER `lambda_range` NOR `training_imsize_provenance`. The plan
#     skeleton asserted `h.lambda_range isa Tuple` and `h.training_imsize_provenance.recorded ==
#     true`; NEITHER PROPERTY EXISTS on this handle and both asserts would throw. Both are bound
#     instead from Phase 11's Tier-1 PRE-REGISTRATION (`spike/validation/p11_consts.jl`):
#     `LAMBDA_MIN`/`LAMBDA_MAX` for the range, `P11_IMSIZE_SET`/`P11_IMSIZE_WEIGHTS` for the
#     provenance. THAT IS A STRONGER SOURCE THAN NET METADATA, not a weaker one: a frozen,
#     pre-registered consts file cannot drift with a retrain, and it sits inside the
#     pre-registration audit trail, whereas a metadata field is whatever the trainer happened to
#     write. What CAN be read off the handle is read off the handle -- the theta arity, the
#     positional prior box including the +/-3.0 shift half-width, the row arithmetic and the
#     research-lane flags -- so the binding is not merely documentary.
# (c) THE TWO PRE-REGISTRATIONS COLLIDE ON A GUARD SENTINEL, so Phase 11's Tier-1 constants are
#     read through an ISOLATED MODULE instead of through `p11_consts.jl`'s own guarded include.
#     `p13/consts.jl` RESERVES the name `P11_DEV_SEED` in its forbidden-seed block (so two
#     research lanes provably cannot share a Philox key) -- and that name is PRECISELY the guard
#     sentinel `p11_consts.jl` keys its whole Tier-1 block on. So once the Phase-13
#     pre-registration is loaded into a module, `p11_consts.jl` believes it has already installed
#     itself there, skips its Tier-1 block, and then dies on its own out-of-guard self-check with
#     `UndefVarError: SC2_SPEARMAN_ATTENUATION`. That is not hypothetical: `runtests.jl` loads
#     `test_p13_consts.jl` before this file, so the naive ordering fails in the suite while
#     passing standalone -- the worst shape of failure to leave in place.
#     NEITHER FILE MAY BE EDITED: both are frozen pre-registrations, and `P11_DEV_SEED` is
#     reserved in `p13/consts.jl` for a good reason. So the repair is entirely local to this
#     file: read `p11_consts.jl` in a private module and bind what is needed from there. The two
#     values that matter (`LAMBDA_MIN`, `LAMBDA_MAX`) are then already in scope when
#     `p11_architecture.jl` is included, so its own guard makes it skip `p11_consts.jl` and the
#     shadowed-sentinel path is never taken. THE POINT OF THE REPAIR IS ORDER-INDEPENDENCE: this
#     file now loads correctly whether the caller reached `p13/consts.jl` first or not, which is
#     the property the guarded-include idiom is supposed to buy everywhere else in the spike.
#
# Flat top-level functions (sibling style of `net.jl` / `harness.jl` -- no module wrapper).
# Guarded, order-dependent includes keep the file loadable standalone AND idempotent under
# `runtests.jl`.

using StatsBase          # ZScoreTransform (the frozen summary standardizer's concrete type)
using Statistics         # mean, std (the inherited-versus-re-fit moment check)

# --- PHASE 11's TIER-1 PRE-REGISTRATION, READ THROUGH AN ISOLATED MODULE ----------------------
# This is the deviation-(c) repair: `p13/consts.jl` reserves `P11_DEV_SEED`, which is exactly the
# guard sentinel `p11_consts.jl` keys its Tier-1 block on, so `p11_consts.jl` can no longer
# install itself in a module that already holds the Phase-13 pre-registration. Reading it in a
# private module sidesteps the collision without editing either frozen file: the module gets a
# clean namespace, so the guard opens and the whole block runs. NO GATE IS RUN AND THE FILE IS
# NOT MODIFIED -- this is a read. (Same precedent as `module _GC` in `p13/consts.jl` and
# `module GateV2` in `p11_consts.jl` itself.)
#
# THIS SITS OUTSIDE ANY GUARD BLOCK ON PURPOSE: Julia rejects a `module` expression that is not
# at top level, so it cannot live inside an `if`. Idempotency comes from the guarded-include
# idiom every caller uses (`isdefined(@__MODULE__, :p13_require_phase11) || include(...)`).
#
# NOTE THE PATH DEPTH: from `spike/p13/` the validation directory is one level up, not two.
module _P11C
    include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
end

# Bind, from that isolated read, the Phase-11 Tier-1 constants this lane needs -- NEVER retyped,
# so they cannot drift from Phase 11's own frozen file. Guarded, so a module that genuinely did
# manage to load `p11_consts.jl` keeps its own bindings and nothing is redefined.
#   LAMBDA_MIN / LAMBDA_MAX  the trained lambda range, and the two values `encode_lambda` closes
#                            over -- putting them in scope HERE is what makes
#                            `p11_architecture.jl` skip its own `p11_consts.jl` include below
#   SC2_RUNGS                the frozen lambda ladder, for the reads plan 13-12 sweeps over
#   P11_IMSIZE_* the F5 training-joint provenance asserted against the Phase-13 mixture
if !isdefined(@__MODULE__, :LAMBDA_MIN)
    const LAMBDA_MIN = _P11C.LAMBDA_MIN
    const LAMBDA_MAX = _P11C.LAMBDA_MAX
    const SC2_RUNGS  = _P11C.SC2_RUNGS
end
if !isdefined(@__MODULE__, :P11_IMSIZE_SET)
    const P11_IMSIZE_SET     = _P11C.P11_IMSIZE_SET
    const P11_IMSIZE_WEIGHTS = _P11C.P11_IMSIZE_WEIGHTS
end
# Any OTHER Phase-11 Tier-1 constant a later Phase-13 file needs is reachable as `_P11C.NAME`
# rather than by adding a second include of `p11_consts.jl` somewhere else in the lane.

# --- ORDER MATTERS. See deviation (a) above: PHASE 11's ARCHITECTURE FIRST, because it owns the
#     concrete theta-transform type the persisted artifact stores, and a `load_npe` without it in
#     scope only WARNS and substitutes a placeholder. Then the Phase-13 pre-registration, then
#     the loader whose `load_npe` reconstructs the handle, then the standardization/read surface,
#     then the encoding surface that derives the input width, then the simulator halves, the src/
#     summary contract, the encoder and the seeding primitives. Guarded for idempotency.
isdefined(@__MODULE__, :P11BoundedThetaTransform) ||
    include(joinpath(@__DIR__, "..", "npe", "p11_architecture.jl"))
isdefined(@__MODULE__, :P13_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :load_npe)          || include(joinpath(@__DIR__, "..", "npe", "train_npe.jl"))
isdefined(@__MODULE__, :standardize_summary) || include(joinpath(@__DIR__, "..", "npe", "infer.jl"))
isdefined(@__MODULE__, :ThreeWayEvidenceNet) || include(joinpath(@__DIR__, "net.jl"))
isdefined(@__MODULE__, :sample_prior)      || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair)     || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :build_mci)         || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :encode_d01)        || include(joinpath(@__DIR__, "..", "data", "encode.jl"))
isdefined(@__MODULE__, :sample_rng)        || include(joinpath(@__DIR__, "..", "data", "seeding.jl"))

if !isdefined(@__MODULE__, :P13_PHASE11_NET)

    # =========================================================================================
    # THE BOUND PHASE-11 ARTIFACT
    # =========================================================================================
    # NAME BOUND BY THE 13-09 TASK-1 BLOCKING CHECKPOINT, resolved 2026-07-28 by inspecting the
    # artifact on disk rather than by reading a planning document. It is the RESEARCH net trained
    # by plan 11-07 (CPU-only, d_in = 129, D = 8, 50000 pairs, research_lane = true,
    # shipped = false). It is NOT the pre-Phase-11 net that sits beside it in the same directory,
    # and it is NOT any shipped bundle: those are a DIFFERENT input surface (128 rows, no
    # conditioning row, the narrow shift prior, 7 theta fields) and binding one would invalidate
    # D-02, D-03 and the D-12 continuity framing.
    const P13_PHASE11_NET = joinpath(@__DIR__, "..", "npe", "p11_research_npe.jld2")

    # The Phase-11 lambda range, bound from Phase 11's Tier-1 pre-registration rather than from
    # net metadata -- see deviation (b). `LAMBDA_MIN` is deliberately never 0 (the 0 -> 0.25 step
    # is interpolation onset, a qualitatively different regime, not misalignment sensitivity).
    const P13_PHASE11_LAMBDA_RANGE = (LAMBDA_MIN, LAMBDA_MAX)

    # The reference lambda that `P13_TAU_REFERENCE_LAMBDA_RULE = :widest_rung` resolves to. Read
    # through the frozen Tier-2 constant so this file never re-derives it, and asserted below to
    # equal Phase 11's own `LAMBDA_MAX` so the rule cannot quietly point at a different rung.
    const P13_PHASE11_REFERENCE_LAMBDA = P13_TAU_REFERENCE_LAMBDA

    # --- Diagnostic tolerances for the inherited-versus-re-fit moment check -------------------
    # LOCAL DIAGNOSTIC KNOBS, NOT PRE-REGISTERED BARS. They do not gate any reported number; they
    # set the sensitivity of a warning-sign heuristic (see `assert_frozen_zt`). They live here
    # rather than in `consts.jl` because `consts.jl` is the pre-registration and must not grow a
    # knob that no claim depends on.
    const P13_FROZEN_ZT_MEAN_TOL = 0.05   # |per-row mean| below this counts as "centred"
    const P13_FROZEN_ZT_SD_TOL   = 0.05   # |per-row sd - 1| below this counts as "unit-scaled"
    const P13_FROZEN_ZT_MIN_COLS = 200    # below this the per-row moments are too noisy to judge

    # Memoized single load (the artifact is ~5 MB; every later p13 file shares ONE read) and the
    # memoized derived conditioning length.
    const _P13_BASIS  = Ref{Any}(nothing)
    const _P13_NCOND  = Ref{Int}(0)

end

"""
    p13_require_phase11(path = P13_PHASE11_NET) -> NamedTuple

Bind the Phase-11 research NPE, or BLOCK.

When `path` is absent this throws the D-02 block: an instructional, multi-line error naming
Phase 11, naming the expected artifact, and stating that Phase 13 must NOT fall back to the
shipped grid-8 basis. IT IS A HARD `error`, NEVER A STATUS RETURN AND NEVER A SKIP. The
`:not_trained` status shape used by the Phase-5 runner scripts is right for a RUNNER, which may
legitimately decline to score; it is wrong here, because D-02 makes a missing Phase-11 net a
block on the phase rather than a gate that reports "not run".

When `path` is present the handle is loaded through `load_npe` and the SIX-ITEM PHASE-11
CONTRACT is asserted item by item, each with its own message naming what it protects:

  1. the trained net is present, and its metadata marks it the RESEARCH lane, not a shipped one;
  2. the FROZEN summary transform is present and is a real `ZScoreTransform`;
  3. the FROZEN theta transform is present AND IS THE CONCRETE `P11BoundedThetaTransform` --
     see deviation (a) at the head of this file: a JLD2-reconstructed placeholder passes a
     field-name check and is not the real transform;
  4. the lambda range is a two-element real range whose upper end equals the reference lambda
     `P13_TAU_REFERENCE_LAMBDA_RULE = :widest_rung` resolves to, so the frozen tau is
     interpretable against this net;
  5. the imsize provenance equals `P13_IMSIZE_SET` / `P13_IMSIZE_WEIGHTS` -- the F5
     covariate-shift check, in the assert-then-print shape of the spike-014 precedent;
  6. the theta arity and FIELD ORDER, asserted POSITIONALLY against the live prior box
     `p11_theta_prior_bounds()`, so an arity change or a reordered row cannot silently shift
     every theta row underneath the phase.

Plus the ROW ARITHMETIC, which closes the input width without a literal:
`d_in == 2 * length(zt.mean) + p13_conditioning_length()`.

Returns the loaded handle.
"""
function p13_require_phase11(path = P13_PHASE11_NET)
    isfile(path) || error("""
      Phase 13 is BLOCKED on Phase 11 (13-CONTEXT D-02).
      Expected the Phase-11 registration-aware research NPE at:
          $path
      Phase 13 must NOT fall back to the shipped grid-8 basis -- that is a different input
      surface and would invalidate D-02, D-03 and the D-12 continuity framing.

      There is no substitute and no fallback path: not the pre-Phase-11 net beside it, not any
      shipped bundle, not a re-standardized Phase-13 pool. Phase 11 must deliver the trained
      research net together with its FROZEN summary transform, its FROZEN theta transform, its
      lambda encoder and range, its recorded image-size and prior provenance, and its 8-field
      theta arity. Until it does, Phase 13 does not run.""")

    h = load_npe(path)

    # --- item 1: the trained net, and the research lane ------------------------------------
    hasproperty(h, :estimator) && h.estimator !== nothing ||
        error("Phase-11 contract item 1: the loaded handle carries no trained estimator ($path)")
    hasproperty(h, :meta) && h.meta !== nothing ||
        error("Phase-11 contract item 1: the loaded handle carries no metadata, so the research " *
              "lane cannot be verified ($path)")
    haskey(h.meta, :research_lane) && h.meta.research_lane === true ||
        error("Phase-11 contract item 1: meta.research_lane is not true -- this is not the " *
              "Phase-11 RESEARCH net, and Phase 13 must not train on any other basis (D-02)")
    haskey(h.meta, :shipped) && h.meta.shipped === false ||
        error("Phase-11 contract item 1: meta.shipped is not false -- a shipped bundle is a " *
              "DIFFERENT input surface and binding it would invalidate D-02/D-03/D-12")

    # --- items 2 and 3: the two FROZEN transforms -------------------------------------------
    hasproperty(h, :zt) && hasproperty(h, :θzt) ||
        error("Phase-11 contract items 2-3: the handle is missing a frozen transform " *
              "(properties present: $(propertynames(h)))")
    h.zt isa StatsBase.ZScoreTransform ||
        error("Phase-11 contract item 2: the frozen summary transform is a " *
              "$(typeof(h.zt)), not a ZScoreTransform -- the standardizer is not the " *
              "deployed one")
    # DEVIATION (a) REGRESSION GUARD. Without `p11_architecture.jl` in scope, JLD2 WARNS and
    # substitutes a reconstructed placeholder whose field names match. This assertion is the
    # difference between that failure being loud and being invisible.
    h.θzt isa P11BoundedThetaTransform ||
        error("Phase-11 contract item 3: the frozen theta transform loaded as a " *
              "$(typeof(h.θzt)), not a P11BoundedThetaTransform. p11_architecture.jl was not " *
              "in scope at load time, so JLD2 reconstructed a PLACEHOLDER and the frozen " *
              "theta basis is silently wrong. Include spike/npe/p11_architecture.jl first.")

    # --- item 6: theta arity and POSITIONAL field order -------------------------------------
    # Asserted against the LIVE prior box rather than against copied literals, so the check moves
    # with `prior.jl` and a reordered field cannot pass by coincidence. Field order is
    # rho_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise,
    # chromatic_eps (prior.jl:93; chromatic_eps appended LAST so every existing index stays
    # valid).
    bounds = p11_theta_prior_bounds()
    length(h.θzt.lo) == length(bounds) ||
        error("Phase-11 contract item 6: the frozen theta transform has " *
              "$(length(h.θzt.lo)) rows but the live prior box has $(length(bounds)) -- the " *
              "theta arity moved under the artifact")
    h.θzt.lo == Float64[b[1] for b in bounds] && h.θzt.hi == Float64[b[2] for b in bounds] ||
        error("Phase-11 contract item 6: the frozen theta box disagrees with the live prior " *
              "box POSITIONALLY. Persisted lo=$(h.θzt.lo) hi=$(h.θzt.hi); live " *
              "lo=$(Float64[b[1] for b in bounds]) hi=$(Float64[b[2] for b in bounds]). A " *
              "reordered or re-widened theta row shifts every row of the frozen basis.")
    # The widened shift half-width IS readable off the handle, so it is asserted off the handle
    # (see deviation (b)): rows 5-6 are shift_dx/shift_dy at Phase 11's widened +/-3.0.
    h.θzt.hi[5] == h.θzt.hi[6] == LAMBDA_MAX && h.θzt.lo[5] == h.θzt.lo[6] == -LAMBDA_MAX ||
        error("Phase-11 contract item 6: the persisted shift half-width " *
              "($(h.θzt.lo[5]), $(h.θzt.hi[5])) is not the widened +/-$(LAMBDA_MAX) that " *
              "D-02 requires, so this net did not see the widened registration prior")

    # --- item 4: the lambda range and the :widest_rung reference ----------------------------
    P13_PHASE11_LAMBDA_RANGE isa Tuple{<:Real,<:Real} ||
        error("Phase-11 contract item 4: the lambda range is not a two-element real range")
    P13_PHASE11_LAMBDA_RANGE[1] < P13_PHASE11_LAMBDA_RANGE[2] ||
        error("Phase-11 contract item 4: the lambda range is not ordered low-to-high")
    P13_PHASE11_LAMBDA_RANGE[1] > 0 ||
        error("Phase-11 contract item 4: LAMBDA_MIN is not strictly positive -- a lambda = 0 " *
              "rung is the no-interpolation regime, which Phase 11 excluded by design")
    P13_PHASE11_LAMBDA_RANGE[2] == P13_PHASE11_REFERENCE_LAMBDA ||
        error("Phase-11 contract item 4: LAMBDA_MAX = $(P13_PHASE11_LAMBDA_RANGE[2]) but the " *
              "frozen P13_TAU_REFERENCE_LAMBDA = $(P13_PHASE11_REFERENCE_LAMBDA). " *
              ":widest_rung would resolve to a lambda the measured tau was not measured at, so " *
              "tau would no longer be interpretable against this net.")

    # --- item 5: the F5 imsize provenance (assert, THEN print) -----------------------------
    P11_IMSIZE_SET == P13_IMSIZE_SET ||
        error("Phase-11 contract item 5: training imsize_set $(P11_IMSIZE_SET) != the " *
              "Phase-13 F5 mixture $(P13_IMSIZE_SET) -- covariate shift between the training " *
              "joint and the evaluation joint")
    P11_IMSIZE_WEIGHTS == P13_IMSIZE_WEIGHTS ||
        error("Phase-11 contract item 5: training imsize_weights $(P11_IMSIZE_WEIGHTS) != the " *
              "Phase-13 F5 weights $(P13_IMSIZE_WEIGHTS)")

    # --- the row arithmetic: the input width closes WITHOUT a literal ----------------------
    # `zt` standardizes the CONTINUOUS rows only (the mask rows are bypassed by
    # `standardize_summary`), so the summary is 2 * length(zt.mean) rows and the net's input is
    # that plus the conditioning rows. Everything here is derived; nothing is typed.
    n_cont  = length(h.zt.mean)
    n_cond  = p13_conditioning_length()
    n_sum   = 2 * n_cont
    h.variant === :min ||
        error("Phase-11 contract: the persisted summary variant is $(h.variant), not :min, so " *
              "the continuous/mask row partition this arithmetic assumes does not hold")
    h.d_in == n_sum + n_cond ||
        error("Phase-11 contract: d_in = $(h.d_in) does not equal 2 * $(n_cont) continuous " *
              "summary rows + $(n_cond) conditioning row(s) = $(n_sum + n_cond). Either the " *
              "conditioning encoder is not the one this net was trained with, or the summary " *
              "row partition changed.")
    G = isqrt(n_cont)
    G * G == n_cont ||
        error("Phase-11 contract: $(n_cont) continuous rows is not a square patch grid")

    println("Phase-11 basis bound OK: ", basename(path),
            " | research_lane=", h.meta.research_lane, " shipped=", h.meta.shipped,
            " | d_in=", h.d_in, " = 2*", n_cont, " + ", n_cond,
            " | G=", G, " | theta arity=", length(h.θzt.lo),
            " | lambda in ", P13_PHASE11_LAMBDA_RANGE,
            " | imsize provenance ", P11_IMSIZE_SET, " w=", P11_IMSIZE_WEIGHTS)

    return h
end

"""
    load_p13_basis(; path = P13_PHASE11_NET) -> NamedTuple

The SINGLE-LOAD accessor every later Phase-13 file uses. Loads and contract-checks the Phase-11
basis on first call through `p13_require_phase11`, memoizes the handle in a module-level `Ref`,
and returns the same object on every later call -- so the ~5 MB artifact is read ONCE per session
and, more importantly, every Phase-13 file demonstrably shares ONE `zt` object.

The analog is `spike/validation/harness.jl:70-80`'s `load_frozen_model`, pointed at the Phase-11
research artifact instead of the shipped Phase-4 net. The difference matters and is the whole
point of D-02: same shape of accessor, different basis.

Passing an explicit `path` that differs from the memoized one RELOADS and re-checks, so a test
can bind a temporary artifact without poisoning the session's memoized basis.
"""
function load_p13_basis(; path = P13_PHASE11_NET)
    if _P13_BASIS[] === nothing || path != P13_PHASE11_NET
        h = p13_require_phase11(path)
        path == P13_PHASE11_NET && (_P13_BASIS[] = h)
        return h
    end
    return _P13_BASIS[]
end

"""
    p13_conditioning_length() -> Int

`n_cond`, DERIVED as `length(encode_lambda(lambda))` from Phase 11's OWN encoder, evaluated once
at the reference lambda and memoized.

**THE PHASE-13 INPUT WIDTH IS `three_way_input_dim(G, p13_conditioning_length())`. NEVER WRITE
THE REALIZED WIDTH AS A LITERAL.** Phase 11 delivers a scalar encoding today, so `n_cond` is 1
and the width at `G = 8` comes out at `5*G^2 + 1`; a Phase 11 that delivered a two-component
encoding would change the width and a hard-coded number would be silently wrong with nothing
thrown. `length` of a scalar is 1 in Julia, so the same expression covers both cases without a
branch.

Computed from `encode_lambda` DIRECTLY rather than from `p13_encode_lambda`, because
`p13_encode_lambda` asserts its own return length against this function -- routing it through the
wrapper would make the two mutually recursive.
"""
function p13_conditioning_length()
    if _P13_NCOND[] == 0
        _P13_NCOND[] = length(encode_lambda(P13_PHASE11_REFERENCE_LAMBDA))
    end
    return _P13_NCOND[]
end

"""
    p13_encode_lambda(lambda) -> Vector{Float64}

Thin, documented PASS-THROUGH to Phase 11's `encode_lambda` (`spike/npe/p11_architecture.jl`),
which maps the registration-uncertainty level onto `[0, 1]` as
`(lambda - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN)`.

**THE ENCODING IS INHERITED, NEVER RE-DERIVED.** The Phase-11 net was trained with that exact
affine map; re-deriving it here -- rescaling, clamping, standardizing, or "improving" it -- would
feed the net a conditioning value it never saw at a row it reads, and nothing would throw. This
function exists ONLY to normalize the container (Phase 11 returns a scalar; the pair encoder
needs an iterable) and to assert the length against `p13_conditioning_length`.

Values outside the trained range are NOT clamped, exactly as Phase 11 declines to clamp them:
reading the net outside its trained lambda range is an extrapolation the caller must own, and
silently clamping would hide it.
"""
function p13_encode_lambda(lambda)
    u = encode_lambda(lambda)
    v = u isa AbstractVector ? collect(Float64, u) : Float64[u]
    length(v) == p13_conditioning_length() ||
        error("p13_encode_lambda: the Phase-11 encoder returned $(length(v)) component(s) at " *
              "lambda = $lambda but p13_conditioning_length() is " *
              "$(p13_conditioning_length()) -- the conditioning width is not stable, so the " *
              "net's input width is not well defined")
    return v
end

"""
    p13_encode_pair(Zs, Zc, lambda) -> Vector{Float32}

The Phase-13 net's input for ONE acquisition: `encode_conditioned_pair(Zs, Zc,
p13_encode_lambda(lambda))`. `Zs` and `Zc` are the sample and control summaries ALREADY
standardized through the FROZEN Phase-11 `zt` (`2*G^2` rows each, mask rows bypassed); ONE lambda
governs the pair, for the reason stated in this file's banner.

THE CORRECT FORM AND BOTH WRONG ONES (D-03, `P13_LAMBDA_PLACEMENT = :append_after_pair_encode`):

    CORRECT     vcat(p13_pair_encode(Zs_2G2, Zc_2G2), cond)   # the conditioning goes LAST
    WRONG       p13_pair_encode(Zs_2G2p1, Zc_2G2p1)           # the lambda folded INTO each
                                                              # summary, so each is odd-length
    ALSO WRONG  vcat(p13_pair_encode(Zs, Zc), cond, cond)     # the conditioning counted twice

The FIRST wrong form fails LOUDLY on `p13_pair_encode`'s `iseven` guard, and THAT IS A GIFT: the
guard is what turns a mis-shaped input into an exception instead of a net that trains happily on
a conditioning value smeared through the contrast block. Never weaken that guard to accommodate a
caller. The SECOND wrong form is SILENT, which is why `encode_conditioned_pair` re-asserts the
realized length against `three_way_input_dim(G, length(cond))` rather than trusting it.
"""
p13_encode_pair(Zs::AbstractVector, Zc::AbstractVector, lambda) =
    encode_conditioned_pair(Zs, Zc, p13_encode_lambda(lambda))

"""
    assert_frozen_zt(basis, Zstd_sample) -> NamedTuple

Evidence that the standardizer was INHERITED from Phase 11 rather than RE-FIT on the Phase-13
pool. Returns `(verdict, max_abs_mean, max_abs_sd_dev, n_cols, n_rows)` where `verdict` is one of
`:inherited`, `:refit_suspected` or `:indeterminate`.

PITFALL 5, IN FULL. The spike's frozen-stats discipline is that the net and BOTH transforms are
loaded ONCE and standardization is NEVER re-fit downstream: every summary is put into the
estimator's input space with the frozen `zt`, and every theta read un-standardizes with the
frozen theta transform (`spike/validation/train_ratio.jl:45-47,212` -- the ratio net always
standardizes through the FROZEN `m.zt` and no `ZScoreTransform` is ever fitted in the validation
lane). Re-fitting is catastrophic AND INVISIBLE: the net still receives numbers of the right
shape and roughly the right scale, it trains, it converges, and every claim about it silently
refers to a basis the net was not trained on. THE WARNING SIGN is standardized Phase-13 summaries
whose per-row moments come out at mean ~0 and sd ~1 ON THE PHASE-13 POOL -- because a transform
frozen on the PHASE-11 pool, applied to a DIFFERENT pool, should NOT reproduce exactly that.

TWO CHECKS, ONE HARD AND ONE SOFT, AND THE DISTINCTION IS DELIBERATE:

  * HARD: `basis.zt` must be the very object the memoized single load returned (`===`). A caller
    who fitted their own standardizer and packed it into a look-alike handle fails here, and
    that is an outright error because it is unambiguous. NOTE WHAT THIS IMPLIES: `basis` must
    have come from `load_p13_basis()`. A handle obtained from a SEPARATE `p13_require_phase11`
    call is a second read of the same bytes and therefore a DIFFERENT object, and it fails this
    check on purpose -- the discipline being enforced is one load per session, shared by every
    Phase-13 file, not merely one transform that happens to hold equal numbers.
  * SOFT: the moment check is a HEURISTIC. Perfect centring on a different pool is strong
    evidence of a re-fit but not proof -- two pools drawn from similar joints can legitimately
    give similar moments, and a small sample gives noisy ones. So a soft violation is REPORTED
    (a `@warn` plus `verdict = :refit_suspected` in the return value) and NEVER silently
    accepted; the caller decides. With fewer than `P13_FROZEN_ZT_MIN_COLS` columns the moments
    are too noisy to judge and the verdict is `:indeterminate` rather than a false clean bill.

Only the CONTINUOUS rows are examined, because `standardize_summary` bypasses the mask rows: the
mask block is untouched by any transform and its moments say nothing about either hypothesis.
"""
function assert_frozen_zt(basis, Zstd_sample::AbstractMatrix;
                          tol_mean::Real = P13_FROZEN_ZT_MEAN_TOL,
                          tol_sd::Real   = P13_FROZEN_ZT_SD_TOL,
                          min_cols::Integer = P13_FROZEN_ZT_MIN_COLS)
    basis.zt === load_p13_basis().zt ||
        error("assert_frozen_zt: basis.zt is NOT the object the Phase-11 artifact load " *
              "returned. The standardizer must be INHERITED through load_p13_basis(), never " *
              "re-fit and never rebuilt (D-02).")

    n_cont = length(basis.zt.mean)
    size(Zstd_sample, 1) >= n_cont || throw(DimensionMismatch(
        "assert_frozen_zt: the standardized sample has $(size(Zstd_sample, 1)) rows but the " *
        "frozen transform covers $(n_cont) continuous rows"))
    C = @view Zstd_sample[1:n_cont, :]
    n = size(C, 2)

    if n < min_cols
        return (verdict = :indeterminate, max_abs_mean = NaN, max_abs_sd_dev = NaN,
                n_cols = n, n_rows = n_cont)
    end

    m = vec(mean(C, dims = 2))
    s = vec(std(C, dims = 2))
    max_abs_mean   = maximum(abs, m)
    max_abs_sd_dev = maximum(abs, s .- 1)
    looks_refit    = max_abs_mean <= tol_mean && max_abs_sd_dev <= tol_sd

    if looks_refit
        @warn """
        assert_frozen_zt: RE-FIT WARNING SIGN (Pitfall 5). The standardized Phase-13 summaries
        are centred and unit-scaled to within the diagnostic tolerance on THIS pool
        (max |mean| = $max_abs_mean <= $tol_mean, max |sd - 1| = $max_abs_sd_dev <= $tol_sd over
        $n columns). A transform FROZEN on the Phase-11 pool should NOT reproduce that on a
        different pool. Confirm the summaries came through the frozen basis.zt and that no
        ZScoreTransform was fitted anywhere in the Phase-13 data path. Reported, NOT accepted.
        """
    end

    return (verdict = looks_refit ? :refit_suspected : :inherited,
            max_abs_mean = max_abs_mean, max_abs_sd_dev = max_abs_sd_dev,
            n_cols = n, n_rows = n_cont)
end

assert_frozen_zt(basis, Zstd_sample::AbstractVector; kwargs...) =
    assert_frozen_zt(basis, reshape(Zstd_sample, :, 1); kwargs...)
