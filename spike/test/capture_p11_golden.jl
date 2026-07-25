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

# spike/test/capture_p11_golden.jl --- PRE-EDIT golden capture for the D-10 stage-6 regression.
#
# TIMING DISCIPLINE. THIS SCRIPT MUST RUN, AND ITS ARTIFACT MUST BE COMMITTED, BEFORE A SINGLE
# BYTE OF `spike/simulator/forward.jl` STAGE 6 OR `src/amortized/simulator.jl` CHANGES. The
# whole value of the fixture is that it records what the simulator produced BEFORE the composed
# affine warp replaced the translation-only warp; a fixture regenerated AFTERWARDS would compare
# the new code against itself and pass vacuously. That failure mode is silent, so it is made
# DETECTABLE two ways: the capture-time commit sha is embedded INSIDE the JLD2 as
# `phase11_base_sha`, and the script REFUSES TO RUN if the working tree is dirty for either of
# the two files whose pre-edit behaviour it is claiming to preserve.
#
# ONE-SHOT. This is a capture script, not a gate: it writes an artifact and asserts nothing
# about the simulator's values. The `==` regression assertion belongs in the stage-6 regression
# test, which reads this fixture back after the edit.
#
# CPU-only (D-10), spike-local, touches no `src/`. Run:
#     julia --project=spike spike/test/capture_p11_golden.jl

using JLD2
using Dates
using Random123

# --- ORDER MATTERS: the pre-registered stream (p11_rng, P11_FIXTURE_COUNTER) must exist
#     before the simulator is loaded and driven by it. Both guarded for idempotency (S2).
isdefined(@__MODULE__, :P11_DEV_SEED) ||
    include(joinpath(@__DIR__, "..", "validation", "p11_consts.jl"))
isdefined(@__MODULE__, :simulate_pair) ||
    include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))

const P11_GOLDEN_PATH   = joinpath(@__DIR__, "fixtures", "p11_stage6_golden.jld2")
const P11_GOLDEN_IMSIZE = (256, 256)
const P11_GOLDEN_KEYS   = 1:4

# An EXPLICIT 7-field θ literal. It is deliberately NOT drawn from the prior sampler: that
# sampler's RNG consumption changes the instant the 8th (chromatic) field is added, so a golden
# captured through it would differ after the edit for a reason that has nothing to do with
# stage 6, and the regression check would be meaningless. Every field sits well inside its
# prior range, so `simulate_pair`'s entry guards accept it unchanged before and after the edit.
const THETA_GOLDEN = (ρ_true           = 0.7,
                      spillover        = 0.1,
                      autofluorescence = 0.05,
                      label_efficiency = 0.85,
                      shift_dx         = 0.37,
                      shift_dy         = -0.82,
                      noise            = 0.5)

"""
    capture_p11_golden(path = P11_GOLDEN_PATH) -> String

Persist the PRE-EDIT `simulate_pair` output at four Philox keys off the reserved FIXTURE
counter, together with the commit sha the capture was taken at. Atomic: `mkpath` → `.tmp` →
`jldsave` → reopen and integrity-check → `mv(...; force = true)`, so a crash leaves only a
discardable `.tmp`, never a half-written golden.
"""
function capture_p11_golden(path::AbstractString = P11_GOLDEN_PATH)
    # --- The pre-edit claim must be TRUE at capture time, not merely intended ---------------
    dirty = readchomp(`git status --porcelain -- spike/simulator src/amortized/simulator.jl`)
    @assert isempty(dirty) """
        capture_p11_golden: the working tree is DIRTY for the simulator files this fixture
        claims to precede:
        $dirty
        A golden captured against uncommitted simulator edits makes the pre-edit claim already
        false. Commit or revert those files, then re-run."""

    base_sha       = readchomp(`git rev-parse HEAD`)
    base_sha_short = readchomp(`git rev-parse --short HEAD`)
    @assert length(base_sha) == 40 "capture_p11_golden: expected a 40-char sha, got $base_sha"

    # --- The golden bytes ------------------------------------------------------------------
    # Key k rides p11_rng(P11_FIXTURE_COUNTER + k): the FIXTURE counter, never one a reported
    # number consumes (D-01), and offset per key so the four draws are independent sub-streams.
    channels = [simulate_pair(p11_rng(P11_FIXTURE_COUNTER + k), THETA_GOLDEN;
                              imsize = P11_GOLDEN_IMSIZE) for k in P11_GOLDEN_KEYS]

    # --- Atomic write (S4) -----------------------------------------------------------------
    mkpath(dirname(path))
    tmp = path * ".tmp"
    jldsave(tmp;
        schema_version         = 1,
        theta                  = THETA_GOLDEN,
        keys_used              = collect(P11_GOLDEN_KEYS),
        imsize                 = P11_GOLDEN_IMSIZE,
        channels               = channels,
        phase11_base_sha       = base_sha,
        phase11_base_sha_short = base_sha_short,
        seed                   = UInt64(P11_DEV_SEED),
        salt                   = UInt64(P11_SALT),
        fixture_counter        = P11_FIXTURE_COUNTER,
        generated              = string(Dates.now(Dates.UTC)) * "Z",
        caption                = "PRE-EDIT simulate_pair golden for the D-10 chromatic-zero " *
                                 "regression check: at chromatic_eps = 0 the composed affine " *
                                 "warp must reproduce these bytes exactly")
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "channels") "capture_p11_golden: integrity check failed, $tmp missing channels"
        @assert haskey(f, "phase11_base_sha") "capture_p11_golden: integrity check failed, $tmp missing phase11_base_sha"
    end
    mv(tmp, path; force = true)   # atomic commit

    println("="^78)
    println("PRE-EDIT stage-6 golden captured (D-10)")
    println("  path            : $path")
    println("  keys            : $(collect(P11_GOLDEN_KEYS)) off counter $P11_FIXTURE_COUNTER")
    println("  imsize          : $P11_GOLDEN_IMSIZE")
    println("  channels        : $(length(channels)) × $(length(first(channels))) matrices")
    println("  stream          : P11_DEV_SEED=$(repr(P11_DEV_SEED)) ⊻ P11_SALT=$(repr(P11_SALT))")
    println()
    println("  PHASE11_BASE_SHA (full ) : $base_sha")
    println("  PHASE11_BASE_SHA (short) : $base_sha_short")
    println()
    println("  Record the full sha in the plan summary — the decoupling/provenance diff")
    println("  commands consume it as `git diff <PHASE11_BASE_SHA>..HEAD -- <paths>`.")
    println("="^78)
    return path
end

capture_p11_golden()
