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

# spike/test/capture_p12_golden.jl --- PRE-EDIT golden capture for the D-06 stage-1 regression.
#
# TIMING DISCIPLINE. THIS SCRIPT MUST RUN, AND ITS ARTIFACT MUST BE COMMITTED, BEFORE A SINGLE
# BYTE OF `spike/simulator/forward.jl` STAGE 1 CHANGES. The whole value of the fixture is that it
# records what the simulator produced BEFORE the scalar rho mixer learned to accept a spatially-
# varying rho FIELD; a fixture regenerated AFTERWARDS would compare the new code against itself
# and pass vacuously. That failure mode is silent, so it is made DETECTABLE two ways: the
# capture-time commit sha is embedded INSIDE the JLD2 as `p12_base_sha`, and the script REFUSES
# TO RUN if the working tree is dirty for the simulator files whose pre-edit behaviour it is
# claiming to preserve.
#
# ONE-SHOT. This is a capture script, not a gate: it writes an artifact and asserts nothing about
# the simulator's values. The `==` regression assertion belongs in `test_p12_prior.jl`, which
# reads this fixture back after the edit.
#
# CPU-only, spike-local, touches no `src/`. Run:
#     julia --project=spike spike/test/capture_p12_golden.jl

using JLD2
using Dates
using Random123

# --- ORDER MATTERS: the pre-registered stream (p12_fix_rng, P12_FIXTURE_COUNTER) must exist
#     before the simulator is loaded and driven by it. Both guarded for idempotency (S2).
isdefined(@__MODULE__, :P12_DEV_SEED) ||
    include(joinpath(@__DIR__, "..", "validation", "p12_consts.jl"))
isdefined(@__MODULE__, :simulate_pair) ||
    include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))

const P12_GOLDEN_PATH = joinpath(@__DIR__, "fixtures", "p12_stage1_golden.jld2")

# FIXTURE GEOMETRY, BOUND HERE AND DELIBERATELY NOT IN `p12_consts.jl`. These two are not
# pre-registered thresholds -- nothing gates on them and no reported number reads them -- and
# Tier 1 is CLOSED (append-only, three reserved Tier-2 sentinels, none of them this). 128^2
# clears `simulate_pair`'s >= 64 guard with margin, keeps the committed fixture small, and still
# exercises every one of the seven stages end to end.
const P12_GOLDEN_IMSIZE = (128, 128)
const P12_GOLDEN_KEYS   = 1:6

# An EXPLICIT 8-field theta literal. It is deliberately NOT drawn from any prior sampler, and the
# Phase-11 reason transfers with an extra edge: `sample_p12_prior`'s RNG consumption differs from
# `sample_prior`'s BY CONSTRUCTION (it draws r1 and a 64-cell field first), so a golden captured
# through a sampler would differ after the edit for a reason that has nothing to do with stage 1,
# and the regression check would be meaningless. This is the Phase-11 literal plus an EXPLICIT
# `chromatic_eps = 0.0`, so the fixture is about stage 1 ALONE and the stage-6 term is pinned
# inert. Every field sits well inside its prior range, so `simulate_pair`'s entry guards accept it
# unchanged before and after the edit.
const THETA_GOLDEN_P12 = (ρ_true           = 0.7,
                          spillover        = 0.1,
                          autofluorescence = 0.05,
                          label_efficiency = 0.85,
                          shift_dx         = 0.37,
                          shift_dy         = -0.82,
                          noise            = 0.5,
                          chromatic_eps    = 0.0)

# THE GOLDEN theta CARRIES NO FIELD, AND THAT IS THE CLAIM THE CONSUMER REFUSES A REGENERATED
# FIXTURE ON. A fixture captured after the edit through a Phase-12 draw would carry `rho_field`.
@assert !hasproperty(THETA_GOLDEN_P12, :rho_field) """
    capture_p12_golden: THETA_GOLDEN_P12 must NOT carry a rho_field -- the exact-equality arm is
    the LEGACY SCALAR path, which never enters the upsample."""
@assert length(THETA_GOLDEN_P12) == 8 "capture_p12_golden: expected the 8-field θ, got $(length(THETA_GOLDEN_P12))"

"""
    capture_p12_golden(path = P12_GOLDEN_PATH) -> String

Persist the PRE-EDIT `simulate_pair` output at six Philox keys off the reserved FIXTURE counter,
together with the commit sha the capture was taken at. Atomic: `mkpath` -> `.tmp` -> `jldsave` ->
reopen and integrity-check -> `mv(...; force = true)`, so a crash leaves only a discardable
`.tmp`, never a half-written golden.
"""
function capture_p12_golden(path::AbstractString = P12_GOLDEN_PATH)
    # --- The pre-edit claim must be TRUE at capture time, not merely intended ---------------
    dirty = readchomp(`git status --porcelain -- spike/simulator`)
    @assert isempty(dirty) """
        capture_p12_golden: the working tree is DIRTY for the simulator files this fixture
        claims to precede:
        $dirty
        A golden captured against uncommitted simulator edits makes the pre-edit claim already
        false. Commit or revert those files, then re-run."""

    base_sha       = readchomp(`git rev-parse HEAD`)
    base_sha_short = readchomp(`git rev-parse --short HEAD`)
    @assert length(base_sha) == 40 "capture_p12_golden: expected a 40-char sha, got $base_sha"

    # --- The golden bytes ------------------------------------------------------------------
    # Key k rides p12_fix_rng(P12_FIXTURE_COUNTER + k): the FIXTURE stream on the FIXTURE
    # counter (p12_consts.jl §2), never one a reported number consumes, and offset per key so the
    # six draws are independent sub-streams.
    channels = [simulate_pair(p12_fix_rng(P12_FIXTURE_COUNTER + k), THETA_GOLDEN_P12;
                              imsize = P12_GOLDEN_IMSIZE) for k in P12_GOLDEN_KEYS]

    # --- Atomic write (S4) -----------------------------------------------------------------
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp;
        schema_version  = 1,
        theta           = THETA_GOLDEN_P12,
        channels        = channels,
        imsize          = P12_GOLDEN_IMSIZE,
        golden_keys     = collect(P12_GOLDEN_KEYS),
        p12_base_sha    = base_sha,
        p12_base_sha_short = base_sha_short,
        seed            = UInt64(P12_FIX_SEED),
        salt            = UInt64(P12_SALT),
        fixture_counter = P12_FIXTURE_COUNTER,
        generated       = string(Dates.now(Dates.UTC)) * "Z",
        caption         = "PRE-EDIT simulate_pair golden for the D-06 stage-1 regression: a θ " *
                          "carrying NO rho_field must reproduce these bytes EXACTLY after stage " *
                          "1 learns to accept a spatially-varying ρ field")
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "channels") "capture_p12_golden: integrity check failed, $tmp missing channels"
        @assert haskey(f, "p12_base_sha") "capture_p12_golden: integrity check failed, $tmp missing p12_base_sha"
    end
    mv(tmp, path; force = true)   # atomic commit

    println("="^78)
    println("PRE-EDIT stage-1 golden captured (D-06)")
    println("  path            : $path")
    println("  keys            : $(collect(P12_GOLDEN_KEYS)) off counter $P12_FIXTURE_COUNTER")
    println("  imsize          : $P12_GOLDEN_IMSIZE")
    println("  channels        : $(length(channels)) × $(length(first(channels))) matrices")
    println("  stream          : P12_FIX_SEED=$(repr(P12_FIX_SEED)) ⊻ P12_SALT=$(repr(P12_SALT))")
    println()
    println("  P12_BASE_SHA (full ) : $base_sha")
    println("  P12_BASE_SHA (short) : $base_sha_short")
    println()
    println("  Record the full sha in the plan summary -- testset 12 asserts it is an ANCESTOR")
    println("  of HEAD, which is what makes 'the fixture precedes the edit' checkable.")
    println("="^78)
    return path
end

capture_p12_golden()
