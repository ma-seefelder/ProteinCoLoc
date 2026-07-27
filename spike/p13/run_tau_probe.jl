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

# spike/p13/run_tau_probe.jl --- REPORTED D-06 resolution probe: LOCKED spec, RESERVED stream.
#
# ANTI-SNOOPING / HONESTY CONTRACT (D-04, D-06, CLAUDE.md HARD_GATE).
# EVERY knob this script obeys -- the delta grid, the AUC bar, R, the statistic, the
# unpaired design, the both-directions rule, the bootstrap count, the reference-lambda
# rule, the abort criterion and the two Phase-13 streams -- was committed to
# `spike/p13/consts.jl` BEFORE this file existed and BEFORE a single reported draw was
# taken. NOTHING below selects, tunes or derives a threshold: every bar is REFERENCED
# from the frozen file by name and never inlined as a literal, so a reader diffing this
# script against the pre-registration can see at a glance that no number moved.
#
# THE GRID IS NEVER EXTENDED. If no delta on `P13_TAU_DELTA_GRID` reaches
# `P13_TAU_AUC`, the run reports `:below_resolution` and STOPS. That is an HONEST
# FINDING about the fixed 8x8 patch-correlation summary -- it says the summary cannot
# reliably separate ANY pre-registered rho difference, and therefore that a defensible
# three-way cut cannot be drawn at this summary's resolution. It is NOT a reason to add
# a grid point, lower the bar, raise R, or reroll a seed. Adding a larger delta after
# SEEING a miss would make "below resolution" unfalsifiable, which is exactly why
# `P13_TAU_ABORT_EXTEND_GRID = false` is frozen and is READ rather than described.
#
# ONE RUN. This script consumes `p13_rng(P13_TAU_COUNTER)` and the two arm keys carved
# from `P13_DEV_SEED` ONCE. A SECOND reported run on a fresh counter would be a second
# measurement of a pre-registered quantity and would have to be declared as spending
# `P13_ITERATION_ALLOWANCE` -- the single documented iteration, which is RESERVED for
# the D-07 stratification switch and for nothing else.
#
# BLOCKED ON THE PHASE-11 SIMULATOR MERGE, NOT ON PHASE 11'S TRAINED NET. The probe
# stops at `encode_d01(patch_summary(...))`, so it needs no `zt`, no conditioning
# encoder and no network. But a tau measured against the PRE-Phase-11 simulator would
# describe an EASIER world than the net will train in (narrower shift prior, no
# chromatic term), so it would be OPTIMISTIC and would starve the random class
# (13-RESEARCH Pitfall 10). `simulator_provenance_guard()` is therefore turned into a
# HARD STOP below, before the reserved stream is spent.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local. `src/` is reached ONLY
# read-only and transitively through `spike/contract.jl` for the UNCHANGED
# patch/correlation summary. No src/ byte, no manifest byte, no new dependency.
#
# CPU-only (CLAUDE.md). Run:
#     julia --project=spike -t auto spike/p13/run_tau_probe.jl

using JLD2
using Dates
using SHA          # Julia standard library: the pre-registration fingerprint. No dependency added.

# The probe machinery (plan 13-05). Guarded, and including it simulates nothing.
isdefined(@__MODULE__, :tau_curve) || include(joinpath(@__DIR__, "tau_probe.jl"))

# Phase-11's Tier-1 pre-registration, loaded into an ISOLATED module purely to READ
# `LAMBDA_MAX` and the SC2 rung ladder so `P13_TAU_REFERENCE_LAMBDA_RULE = :widest_rung`
# can be RESOLVED against Phase 11's own frozen numbers instead of retyped. Isolated
# because that file re-declares names this scope already holds (`P11_DEV_SEED` among
# them) and because nothing of Phase 11 is RUN here. Same idiom as `_GC` in
# `spike/p13/consts.jl`, and it must sit at top level for the same reason.
module _P11C
    include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
end

const TAU_REPORT_PATH = joinpath(@__DIR__, "tau_probe_report.jld2")
const REPO_ROOT       = normpath(joinpath(@__DIR__, "..", ".."))

# --- Provenance helpers ---------------------------------------------------------------
# `git -C REPO_ROOT` rather than a bare `git`, so the recorded shas describe THIS
# checkout whatever directory the script was launched from.
_repo_head_sha()  = readchomp(`git -C $REPO_ROOT rev-parse HEAD`)
_simulator_sha()  = readchomp(`git -C $REPO_ROOT log -1 --format=%H -- spike/simulator`)
# Recorded because a HEAD sha alone is NOT complete provenance when the working tree
# carries uncommitted changes -- and at probe time this runner itself is necessarily
# still uncommitted (the two-commit discipline puts the artifact commit AFTER the run).
# `runner_sha256` below is what actually pins the code that produced the artifact.
_repo_dirty()     = !isempty(readchomp(`git -C $REPO_ROOT status --porcelain`))
_sha256_of(path)  = bytes2hex(SHA.sha256(read(path)))
# The same convention as `p13_consts_sha()` (spike/p13/net.jl:524), recomputed locally
# rather than included: `net.jl` pulls the Flux / NeuralEstimators stack, and a
# simulator-only probe must not load a training stack to fingerprint a text file.
_consts_sha256()  = _sha256_of(joinpath(@__DIR__, "consts.jl"))
_runner_sha256()  = _sha256_of(@__FILE__)

# --- Atomic persistence (the `_save_bf_report` idiom, spike/validation/run_bf.jl:60) ---
# Write to `.tmp`, REOPEN and assert the load-bearing keys are really there, only then
# `mv(force = true)`. A crash mid-write leaves the `.tmp`, never a torn report at the
# reported path (T-13-03).
function _save_tau_report(path; kwargs...)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "curve")              "_save_tau_report: integrity check failed ($tmp)"
        @assert haskey(f, "P13_TAU_DELTA_GRID") "_save_tau_report: integrity check failed ($tmp)"
        @assert haskey(f, "P13_TAU_AUC")        "_save_tau_report: integrity check failed ($tmp)"
        @assert haskey(f, "simulator_sha")      "_save_tau_report: integrity check failed ($tmp)"
    end
    mv(tmp, path; force = true)
    return path
end

function main()
    t_start = time()

    # --- (1) The banner: every locked knob echoed BY NAME and VALUE ---------------------
    println("="^78)
    println("Phase-13 REPORTED D-06 resolution probe (tau) — LOCKED consts, RESERVED stream")
    println("  grid:      P13_TAU_DELTA_GRID = $P13_TAU_DELTA_GRID")
    println("  bar:       P13_TAU_AUC = $P13_TAU_AUC   (tau = smallest delta with A(delta) >= bar)")
    println("  draws:     P13_TAU_R = $P13_TAU_R per arm")
    println("  statistic: P13_TAU_STATISTIC = $(repr(P13_TAU_STATISTIC))")
    println("  design:    P13_TAU_DESIGN = $(repr(P13_TAU_DESIGN))   " *
            "both directions: P13_TAU_BOTH_DIRECTIONS = $P13_TAU_BOTH_DIRECTIONS")
    println("  bootstrap: P13_TAU_BOOTSTRAP_B = $P13_TAU_BOOTSTRAP_B (median reported)")
    println("  lambda:    P13_TAU_REFERENCE_LAMBDA_RULE = $(repr(P13_TAU_REFERENCE_LAMBDA_RULE)) " *
            "(expected $P13_TAU_REFERENCE_LAMBDA_EXPECTED)")
    println("  abort:     P13_TAU_ABORT_EXTEND_GRID = $P13_TAU_ABORT_EXTEND_GRID  " *
            "(a miss is REPORTED, the grid is NEVER extended)")
    println("  simulator: P13_TAU_REQUIRE_POST_P11_SIMULATOR = $P13_TAU_REQUIRE_POST_P11_SIMULATOR")
    println("  frames:    P13_IMSIZE_SET = $P13_IMSIZE_SET")
    println("             P13_IMSIZE_WEIGHTS = $P13_IMSIZE_WEIGHTS")
    println("  stream:    P13_DEV_SEED = $(repr(P13_DEV_SEED)) at counter P13_TAU_COUNTER = $P13_TAU_COUNTER")
    println("="^78)

    # --- (2) Stream disjointness, asserted BEFORE anything is drawn (T-13-01) -----------
    @assert UInt64(P13_DEV_SEED) != UInt64(P13_FIX_SEED) "the reported stream must be disjoint from the fixture stream (D-01)"
    @assert !(UInt64(P13_DEV_SEED) in _p13_forbidden())  "the reported stream collides with a reserved or burned seed (D-01)"
    @assert P13_TAU_COUNTER != P13_FIXTURE_COUNTER       "the reported counter must not be the fixture counter (D-01)"

    # --- (3) Simulator provenance: the RUNTIME enforcement of the Task-1 checkpoint -----
    # `simulator_provenance_guard` REPORTS; this is the caller that turns a false into a
    # hard stop, and it runs BEFORE the curve is computed so the reserved stream is never
    # spent against the wrong simulator.
    guard = simulator_provenance_guard()
    guard.post_p11 || error("""
        REFUSING TO SPEND THE RESERVED STREAM: the Phase-11 simulator surgery has NOT merged.

        $(guard.reason)

        13-RESEARCH Pitfall 10: the POST-Phase-11 deployment distribution is strictly WIDER
        (widened SHIFT_PRIOR, added chromatic term), and a wider nuisance joint can only make
        two rho values HARDER to tell apart. A tau measured against the pre-Phase-11 simulator
        is therefore OPTIMISTIC: it would set the random band too narrow, starve the random
        class, and hand the D-12 confusion matrix a flattering result that came from measuring
        an easier world than the net trains in. Freezing such a tau is a pre-registration error.
        """)
    println("\nsimulator provenance: ", guard.reason)

    # `:widest_rung` RESOLVED, not retyped: the widest SC2 rung must be Phase-11's
    # LAMBDA_MAX, which must equal the frozen expectation AND the half-width the guard
    # actually read off the live SHIFT_PRIOR. Resolution degrades with registration
    # uncertainty, so tau is really tau-of-lambda and the reference must travel with it.
    lambda_ref = float(last(_P11C.SC2_RUNGS))
    @assert isapprox(lambda_ref, float(_P11C.LAMBDA_MAX); rtol = 1e-9) "the widest SC2 rung is not Phase-11's LAMBDA_MAX"
    @assert isapprox(lambda_ref, float(P13_TAU_REFERENCE_LAMBDA_EXPECTED); rtol = 1e-9) "the realized widest rung disagrees with P13_TAU_REFERENCE_LAMBDA_EXPECTED"
    @assert isapprox(guard.shift_halfwidth, lambda_ref; rtol = 1e-9) "the live shift-prior half-width disagrees with the widest rung"
    println("reference lambda ($(repr(P13_TAU_REFERENCE_LAMBDA_RULE))) resolved to $lambda_ref " *
            "= Phase-11 LAMBDA_MAX = the live shift-prior half-width")
    # The probe MARGINALIZES over the full SHIFT_PRIOR rather than conditioning on a rung,
    # and that prior's half-width IS the widest rung -- so the measured tau is already the
    # tau at the widest registration uncertainty, which is the conservative reading the
    # rule asks for.

    # --- (4) Provenance shas, recorded before the measurement --------------------------
    repo_sha   = _repo_head_sha()
    sim_sha    = _simulator_sha()
    dirty      = _repo_dirty()
    consts_sha = _consts_sha256()
    runner_sha = _runner_sha256()
    println("repo HEAD        = $repo_sha", dirty ? "  (working tree DIRTY at probe time)" : "")
    println("simulator commit = $sim_sha")
    println("consts.jl sha256 = $consts_sha")
    println("runner    sha256 = $runner_sha")

    # --- (5) The WHOLE curve, never short-circuited at the first clearing delta ---------
    # The curve is the reported artifact: it lets a reader re-read tau at ANY bar without
    # re-running the probe, and a FLAT or NON-MONOTONE curve is a different and more
    # serious finding than "tau is large". Every other knob is left at its frozen default
    # -- `arm_keys` is passed only so the exact Philox first-key-words can be persisted,
    # and it is precisely the value the default expression produces.
    arm = tau_probe_arm_keys()
    n_draws = P13_TAU_R * (1 + length(P13_TAU_DELTA_GRID) * (P13_TAU_BOTH_DIRECTIONS ? 2 : 1))
    println("\ncomputing the whole A(delta) curve: $(length(P13_TAU_DELTA_GRID)) deltas, " *
            "$n_draws full simulate-and-summarize draws over the F5 frame mixture …")
    t0 = time()
    curve = tau_curve(; arm_keys = arm)
    t_curve = time() - t0
    println("curve done in $(round(t_curve / 60; digits = 2)) min")

    # --- (6) PERSIST BEFORE ANY VERDICT ------------------------------------------------
    # The measurement is written and integrity-checked before `tau_from_curve` is ever
    # called, so the curve survives whatever the verdict turns out to be.
    cols = (delta            = Float64[r.delta            for r in curve],
            auc_neg          = Float64[r.auc_neg          for r in curve],
            auc_pos          = Float64[r.auc_pos          for r in curve],
            auc              = Float64[r.auc              for r in curve],
            auc_median_boot  = Float64[r.auc_median_boot  for r in curve],
            n_reference      = Int[r.n_reference          for r in curve],
            n_neg            = Int[r.n_neg                for r in curve],
            n_pos            = Int[r.n_pos                for r in curve],
            n_degenerate     = Int[r.n_degenerate         for r in curve])
    spec = (; curve = identity.(curve),
            delta = cols.delta, auc_neg = cols.auc_neg, auc_pos = cols.auc_pos,
            auc = cols.auc, auc_median_boot = cols.auc_median_boot,
            n_reference = cols.n_reference, n_neg = cols.n_neg, n_pos = cols.n_pos,
            n_degenerate = cols.n_degenerate,
            # The locked spec, each knob its own key, so the artifact is self-describing.
            P13_TAU_DELTA_GRID = P13_TAU_DELTA_GRID,
            P13_TAU_AUC = P13_TAU_AUC,
            P13_TAU_R = P13_TAU_R,
            P13_TAU_STATISTIC = P13_TAU_STATISTIC,
            P13_TAU_DESIGN = P13_TAU_DESIGN,
            P13_TAU_BOTH_DIRECTIONS = P13_TAU_BOTH_DIRECTIONS,
            P13_TAU_BOOTSTRAP_B = P13_TAU_BOOTSTRAP_B,
            P13_TAU_REFERENCE_LAMBDA_RULE = P13_TAU_REFERENCE_LAMBDA_RULE,
            P13_TAU_REFERENCE_LAMBDA_EXPECTED = P13_TAU_REFERENCE_LAMBDA_EXPECTED,
            P13_TAU_ABORT_EXTEND_GRID = P13_TAU_ABORT_EXTEND_GRID,
            P13_TAU_REQUIRE_POST_P11_SIMULATOR = P13_TAU_REQUIRE_POST_P11_SIMULATOR,
            P13_IMSIZE_SET = P13_IMSIZE_SET,
            P13_IMSIZE_WEIGHTS = P13_IMSIZE_WEIGHTS,
            # The stream, in full, so the run is re-derivable to the bit.
            P13_DEV_SEED = UInt64(P13_DEV_SEED),
            P13_SALT = UInt64(P13_SALT),
            P13_TAU_COUNTER = P13_TAU_COUNTER,
            arm_key_reference = UInt64(arm.reference),
            arm_key_contrast = UInt64(arm.contrast),
            # Provenance.
            repo_sha = repo_sha,
            simulator_sha = sim_sha,
            repo_dirty_at_run = dirty,
            consts_sha256 = consts_sha,
            runner_sha256 = runner_sha,
            simulator_guard_reason = guard.reason,
            simulator_theta_arity = guard.theta_arity,
            simulator_shift_halfwidth = guard.shift_halfwidth,
            reference_lambda = lambda_ref,
            curve_seconds = t_curve,
            generated = string(Dates.now(Dates.UTC)) * "Z")
    _save_tau_report(TAU_REPORT_PATH; spec...)
    measurement_sha = _sha256_of(TAU_REPORT_PATH)
    println("persisted the MEASUREMENT (pre-verdict) -> $TAU_REPORT_PATH")
    println("pre-verdict artifact sha256 = $measurement_sha")

    # --- The reported table (printed before the verdict, like the BF runner) ------------
    println("\n", "-"^78)
    println("A(delta) over the frozen grid — magnitudes max(A, 1-A), both directions:")
    println(rpad("delta", 10), rpad("A_neg", 10), rpad("A_pos", 10), rpad("A = max", 10),
            rpad("boot med", 10), rpad("n_ref", 8), rpad("n_neg", 8), rpad("n_pos", 8), "n_degen")
    for r in curve
        println(rpad(r.delta, 10),
                rpad(round(r.auc_neg; digits = 4), 10),
                rpad(round(r.auc_pos; digits = 4), 10),
                rpad(round(r.auc; digits = 4), 10),
                rpad(round(r.auc_median_boot; digits = 4), 10),
                rpad(r.n_reference, 8), rpad(r.n_neg, 8), rpad(r.n_pos, 8), r.n_degenerate)
    end
    println("-"^78)

    # --- (7) The verdict ---------------------------------------------------------------
    v = tau_from_curve(curve)
    spacing = ghat_knot_spacing_near_zero()
    println()
    if v.status === :measured
        println("MEASURED tau = $(v.tau)   (grid position $(v.delta_index) of " *
                "$(length(P13_TAU_DELTA_GRID)); A(tau) = $(round(curve[v.delta_index].auc; digits = 4)) " *
                ">= P13_TAU_AUC = $(v.bar))")
        println("  bootstrap median at tau = $(round(curve[v.delta_index].auc_median_boot; digits = 4)) " *
                "over P13_TAU_BOOTSTRAP_B = $P13_TAU_BOOTSTRAP_B resamples")
        println("  D-05 random band = (-tau, +tau) at reference lambda $lambda_ref")
        # --- (8) The sub-knot caution: PRINTED, never an error --------------------------
        if v.sub_knot
            println("  CAUTION (13-RESEARCH Pitfall 8): tau = $(v.tau) is BELOW the ghat knot " *
                    "spacing near zero ($spacing). The frozen calibration map is a " *
                    "piecewise-linear interpolation of a measured sweep, so a tau finer than " *
                    "its knot spacing claims a resolution the map itself does not carry. " *
                    "REPORTED, NOT CORRECTED — state it in the report; do NOT round tau up.")
        else
            println("  tau is at or above the ghat knot spacing near zero ($spacing): no " *
                    "Pitfall-8 sub-knot caution applies.")
        end
    else
        println("!"^78)
        println("BELOW RESOLUTION — NO tau WAS MEASURED.")
        println("No delta on the pre-registered P13_TAU_DELTA_GRID = $P13_TAU_DELTA_GRID reached")
        println("the frozen bar P13_TAU_AUC = $(v.bar). The abort criterion has fired.")
        println()
        println("THE GRID WAS NOT EXTENDED (P13_TAU_ABORT_EXTEND_GRID = $P13_TAU_ABORT_EXTEND_GRID),")
        println("the bar was not lowered, R was not raised and no seed was rerolled. This is a")
        println("REPORTABLE PHASE-LEVEL FINDING: the fixed 8x8 patch-correlation summary cannot")
        println("reliably separate ANY pre-registered rho difference, so Phase 13 cannot draw a")
        println("defensible three-way cut at this summary's resolution and NO labelled data may")
        println("be generated. No Tier-2 P13_TAU is written.")
        println("!"^78)
    end

    # --- Re-persist WITH the verdict ----------------------------------------------------
    # The verdict is a DETERMINISTIC PURE FUNCTION of the curve just persisted, so
    # appending it cannot influence what was measured. `measurement_sha256` pins the
    # pre-verdict file, so the shipped artifact provably carries the same curve that was
    # written before `tau_from_curve` was ever called.
    _save_tau_report(TAU_REPORT_PATH; spec...,
        tau = v.tau === nothing ? NaN : float(v.tau),
        tau_status = string(v.status),
        tau_measured = v.status === :measured,
        tau_delta_index = v.delta_index === nothing ? 0 : v.delta_index,
        tau_bar = float(v.bar),
        tau_sub_knot = v.sub_knot,
        ghat_knot_spacing_near_zero = spacing,
        measurement_sha256 = measurement_sha)
    println("\npersisted the verdict alongside the measurement -> $TAU_REPORT_PATH")
    println("total wall clock: $(round((time() - t_start) / 60; digits = 2)) min")
    return nothing
end

main()
