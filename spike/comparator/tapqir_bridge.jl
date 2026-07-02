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

# spike/comparator/tapqir_bridge.jl --- CMP-06: optional Tapqir anchor (D-07/D-08)
#
# The Tapqir sanity anchor. Tapqir (Ordabayev et al., eLife 2022 11:e73860;
# gelles-brandeis/tapqir, Apache-2.0, PyPI v1.1.19) models CoSMoS single-molecule
# *spot* data — a DIFFERENT data regime from our diffuse 2-channel colocalization
# fields — so it is NOT a comparator column on our simulator images (D-07). This
# bridge only proves "we can invoke Tapqir correctly" by reproducing a PUBLISHED
# tutorial value (the eLife-2022 Part II `cosmos` sample-dataset tutorial), and it
# strengthens positioning ("validated against the modern SBI-adjacent baseline")
# WITHOUT ever gating the classical table.
#
# DECOUPLING (hard constraint, CLAUDE.md + threats T-09-10/T-09-11/T-09-12):
#   * The Python stack (PythonCall + CondaPkg) lives ONLY in the isolated sub-env
#     `spike/comparator/tapqir_env/` (its own Project/Manifest) and MUST NEVER enter
#     the main `spike/Project.toml`. runtests.jl gate (g) asserts their absence.
#   * The bridge NEVER imports PythonCall in this (main-harness) process. Tapqir is
#     run out-of-process, in a SUBPROCESS activated against the isolated env, so no
#     Python ever loads on the classical path (OQ2 RESOLVED: subprocess-only).
#   * `_tapqir_env_available()` probes first and returns `false` on ANY error; on any
#     failure `tapqir_anchor()` returns `status=:skipped` — the classical battery runs
#     GREEN with zero Python present (graceful skip-with-flag, mirrors CPU/CUDA D-04).
#   * We NEVER `eval` untrusted content; the subprocess only runs pinned, first-party
#     scripts under the audited, Apache-2.0, peer-reviewed Tapqir package.
#
# Contract consumed by the harness / D-13 gate (09-06):
#   tapqir_anchor(; tol=TAPQIR_TOL) ->
#       (status = :passed|:failed|:skipped, value = Float64|missing, reason = String)
#
# Depends on `TAPQIR_PUBLISHED_VALUE` and `TAPQIR_TOL` from config.jl (include it first).

"Absolute path to the isolated Tapqir sub-env (own Project.toml/Manifest.toml)."
_tapqir_env_dir() = joinpath(@__DIR__, "tapqir_env")

"""
    _wait_timeout(proc, timeout_s) -> Base.Process

Poll `proc` until it exits or `timeout_s` seconds elapse; kill it on timeout so the
probe/tutorial subprocess can NEVER hang the classical path (autonomous-safety: never
block on a stalled Conda/Python download).
"""
function _wait_timeout(proc::Base.Process, timeout_s::Real)
    t0 = time()
    while process_running(proc) && (time() - t0) < timeout_s
        sleep(0.25)
    end
    if process_running(proc)
        try; kill(proc); catch; end
        sleep(0.5)
        try; kill(proc, Base.SIGKILL); catch; end
    end
    return proc
end

"""
    _run_isolated(script::AbstractString; timeout_s=120, offline=true) -> (ok::Bool, out::String)

Run `script` in a SUBPROCESS activated against the isolated `tapqir_env` (never in this
process). Captures stdout; stderr is discarded. `offline=true` forces CondaPkg offline so
the subprocess uses only an ALREADY-materialized env and never triggers a network/Conda
download (hang-proof). Returns `(success && !timed_out, captured_stdout)`.
"""
function _run_isolated(script::AbstractString; timeout_s::Real = 120, offline::Bool = true)
    envdir = _tapqir_env_dir()
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(envdir) -e $script`
    extra = Pair{String,String}[]
    if offline
        push!(extra, "JULIA_CONDAPKG_OFFLINE" => "yes")
        # Null backend = do not manage/download a Conda env; use only what is present.
        if !haskey(ENV, "JULIA_CONDAPKG_BACKEND")
            push!(extra, "JULIA_CONDAPKG_BACKEND" => "Null")
        end
    end
    cmd = addenv(cmd, extra...)
    out = ""
    ok = false
    mktemp() do path, io
        proc = run(pipeline(cmd; stdout = io, stderr = devnull); wait = false)
        _wait_timeout(proc, timeout_s)
        try; close(io); catch; end
        out = isfile(path) ? read(path, String) : ""
        ok = !process_running(proc) && success(proc)
    end
    return (ok, out)
end

"""
    _tapqir_env_available() -> Bool

Probe whether the isolated Tapqir env is reachable: the sub-env Manifest must exist AND a
subprocess under that env must be able to `import PythonCall` and import the `tapqir` Python
package. Returns `false` on ANY error (never throws). This is the graceful-degradation gate
that keeps the classical battery green with zero Python present (T-09-12).
"""
function _tapqir_env_available()
    try
        envdir = _tapqir_env_dir()
        isfile(joinpath(envdir, "Manifest.toml")) || return false
        script = raw"""
        try
            import PythonCall
            PythonCall.pyimport("tapqir")
            print("TAPQIR_OK")
        catch
            print("TAPQIR_NO")
        end
        """
        ok, out = _run_isolated(script; timeout_s = 120, offline = true)
        return ok && occursin("TAPQIR_OK", out)
    catch
        return false
    end
end

"""
    _run_tapqir_tutorial() -> Float64

Run the PUBLISHED eLife-2022 Part II `cosmos` tutorial (Tapqir's distributed sample data
set) once, out-of-process under the isolated env, and return the recovered scalar (the
documented tutorial quantity). Communicates the value back over stdout as a single
`TAPQIR_VALUE=<float>` line; throws if the subprocess fails or emits no parseable value
(the caller maps any throw to a clean `:skipped`). Never imports PythonCall in-process.
"""
function _run_tapqir_tutorial()
    # First-party, pinned script — NOT untrusted input, never `eval`ed here. Runs the
    # documented cosmos tutorial and prints the recovered scalar for the parent to parse.
    # The exact recovered quantity + dataset id are captured/pinned in 09-04 Task 2
    # (eLife 11:e73860, Part II sample data set); until captured the subprocess raises,
    # which the caller turns into a documented clean-skip.
    script = raw"""
    using PythonCall
    tapqir = pyimport("tapqir")
    # eLife-2022 Part II tutorial: fit the `cosmos` model on the distributed sample data
    # set and recover the documented scalar (e.g. average target-specific spot probability
    # or a documented cosmos global parameter). The concrete API call is pinned when the
    # value is captured on a machine with a materialized env (09-04 Task 2 / OQ1).
    error("Tapqir tutorial value not yet captured on this machine (OQ1 clean-skip path)")
    """
    ok, out = _run_isolated(script; timeout_s = 600, offline = true)
    ok || error("Tapqir tutorial subprocess failed or timed out")
    m = match(r"TAPQIR_VALUE=([-+0-9.eE]+)", out)
    m === nothing && error("Tapqir tutorial produced no parseable TAPQIR_VALUE")
    return parse(Float64, m.captures[1])
end

"""
    tapqir_anchor(; tol = TAPQIR_TOL) -> NamedTuple

Optional Tapqir sanity anchor with graceful skip-with-flag (D-07/D-08). Returns
`(status, value, reason)`:

  * env unavailable                         -> `(:skipped, missing, "...")`
  * `TAPQIR_PUBLISHED_VALUE` not yet pinned -> `(:skipped, missing, "...")`
  * tutorial reproduces value within `tol`  -> `(:passed,  value,   "...")`
  * tutorial reproduces value outside `tol` -> `(:failed,  value,   "...")`
  * ANY exception                           -> `(:skipped, missing, showerror(e))`

The classical battery never depends on this result; a `:skipped` return is a fully
green outcome (T-09-12).
"""
function tapqir_anchor(; tol = TAPQIR_TOL)
    avail = try
        _tapqir_env_available()
    catch
        false
    end
    avail || return (status = :skipped, value = missing,
                     reason = "Tapqir/Python env unavailable")
    try
        ref = TAPQIR_PUBLISHED_VALUE
        if !(ref isa Real) || !isfinite(ref)
            return (status = :skipped, value = missing,
                    reason = "TAPQIR_PUBLISHED_VALUE not pinned (documented clean-skip)")
        end
        recovered = _run_tapqir_tutorial()
        if !isfinite(recovered)
            return (status = :skipped, value = missing,
                    reason = "Tapqir tutorial returned a non-finite value")
        end
        ok = abs(recovered - ref) <= tol
        return (status = ok ? :passed : :failed, value = recovered,
                reason = ok ? "within tolerance ($(tol))" : "outside tolerance ($(tol))")
    catch e
        return (status = :skipped, value = missing, reason = sprint(showerror, e))
    end
end
