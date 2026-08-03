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

# spike/p14/provenance.jl --- the D-07 threshold loader and the Phase-14 provenance block.
#
# TAU IS INHERITED, NOT OWNED. D-07: "tau is INHERITED from Phase 13's artifact by loading it,
# with provenance asserted. Copying the value is how two numbers silently diverge; re-deriving
# guarantees a second, different cut on the same hypothesis space." Phase 13 MEASURED that
# number against a pre-registered bar on a frozen grid; Phase 14 has no licence to have an
# opinion about it. This file is the ONLY route by which any Phase-14 code may obtain it, and
# the route ERRORS rather than defaults whenever the provenance does not check out.
#
# WHY THIS FILE EXISTS AT ALL, RATHER THAN A TIER-2 APPEND TO spike/p14/consts.jl. Phase 13's
# Tier 2 was a value it would MEASURE and APPEND to its own pre-registration. Phase 14's
# equivalent is a value it will LOAD from someone else's file, so there is no second guard
# block in spike/p14/consts.jl and there must never be one -- that file asserts tau's own
# absence three separate ways. The loader needed a home; this is it.
#
# THE FOUR-WAY AGREEMENT, AND WHY IT IS ASSERTED RATHER THAN ASSUMED. The measured threshold
# appears in FOUR places on disk (14-RESEARCH section D.3, every one loaded and verified
# 2026-08-03):
#
#   1. spike/p13/three_way_net.jld2         key "tau"        <-- THE SOURCE OF TRUTH
#   2. spike/p13/tau_probe_report.jld2      key "tau"
#   3. spike/p13/three_way_gate_report.jld2 key "P13_TAU"
#   4. spike/p13/consts.jl Tier-2           const P13_TAU  (in scope here via the include below)
#
# (1) is primary because the threshold travels with the NET THAT WAS TRAINED UNDER IT: a net and
# a cut that disagree is the failure that matters, and reading the constant as primary would not
# catch it. All four agree TODAY. This file asserts that agreement at load time so that a future
# divergence is LOUD rather than silent -- which is the whole of D-07's reasoning.
#
# THE PROBE'S THIRD consts_sha256 IS ASSERTED **DIVERGENT**, AND THAT IS NOT A BUG. See the long
# note at the assertion itself. Widening it to a disjunction is itself the failure.
#
# DECOUPLING (hard constraint, CLAUDE.md and D-01): spike-local. Reaches spike/p13/ READ-ONLY
# through guarded includes and JLD2 loads, writes nothing anywhere, touches no src/ byte and
# adds no dependency -- JLD2 is already among the frozen sixteen; SHA and Dates are stdlib.
#
# NOTE THE PATH DEPTH: from spike/p14/ the repo root is TWO levels up, not three. The same
# reminder is carried in three separate spike/p13/ files because the mistake is silent.
#
# EVERY CALLER USES THE GUARDED INCLUDE:
#     isdefined(@__MODULE__, :P14_PROVENANCE_LOADED) || include(joinpath(@__DIR__, "provenance.jl"))

using JLD2      # already among the frozen sixteen; the artifact format Phase 13 persisted in
using SHA       # stdlib: the sha256 of the frozen pre-registration and the git blob sha1
using Dates     # stdlib: the UTC stamp every persisted Phase-14 artifact carries

# The Phase-14 Tier-1 pre-registration: P14_DEV_SEED, P14_SALT and the frozen bars.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
# The Phase-13 net surface: `load_three_way`, `p13_consts_sha()` and `P13_CONSTS_SHA`.
# INCLUDING THIS FILE TRANSITIVELY BRINGS THE PHASE-13 PRE-REGISTRATION INTO SCOPE -- P13_TAU,
# P13_TAU_AUC, P13_TAU_MEASURED_AUC, P13_TAU_PROBE_SHA, P13_TAU_SIMULATOR_SHA and the accessor
# `p13_tau()`. That is what makes the FOURTH cross-check free: the constant is read at the one
# place it was measured, never mirrored into the Phase-14 namespace.
isdefined(@__MODULE__, :ThreeWayEvidenceNet) || include(joinpath(@__DIR__, "..", "p13", "net.jl"))

# THE GUARD COVERS ONLY THE `const`s. Julia 1.12 DROPS every docstring written inside an
# `if ... end` block (the parser emits the `Core.@doc` call but the docsystem never registers
# it), so the documented functions sit at top level; method redefinition under a re-include is
# silent and harmless, while `const` redefinition is what would warn. Same split, same reason,
# as spike/p13/result.jl.
if !isdefined(@__MODULE__, :P14_PROVENANCE_LOADED)
    const P14_P13_DIR    = normpath(joinpath(@__DIR__, "..", "p13"))
    const P14_REPO_ROOT  = normpath(joinpath(@__DIR__, "..", ".."))   # TWO levels up, not three
    const P14_CONSTS_PATH = joinpath(@__DIR__, "consts.jl")

    # The SC1/SC2 amendment citation, carried in CODE so it travels with every artifact rather
    # than living only in a document a reader may not open (14-PATTERNS section 4.6). Copied in
    # shape from the `amendment =` key of spike/p13/run_three_way_gate.jl:763-764.
    const P14_AMENDMENT_NOTICE =
        "SC1 is AMENDED by D-02 (ConformalPrediction.jl dropped; conformal met in substance, " *
        "hand-rolled). SC2 is AMENDED by D-05 (asymmetric fusion, not the literal OR). Any report citing " *
        "this result must cite the original ROADMAP SC1/SC2 alongside it."

    # Declared UNCONDITIONALLY at the foot of the block, so it is the sentinel and nothing else
    # can make this body skip.
    const P14_PROVENANCE_LOADED = true
end

# --- Provenance primitives ---------------------------------------------------------------------

"""
    p14_blob_sha(path) -> String

The GIT BLOB sha1 of `path`, i.e. exactly what `git hash-object <path>` prints. ATTRIBUTED COPY
of `p13_blob_sha` (`spike/p13/run_three_way_gate.jl:200-206`), reproduced here rather than
reached through an `include` because that file is a RUNNER, not a library, and including it to
borrow one helper would be a far larger surface than copying six lines.

Recorded so a Phase-14 artifact points INTO the audit trail: a reader can run
`git cat-file -p <sha>` to get the frozen pre-registration text out of history, or
`git log --find-object=<sha>` to get the commit that introduced it.
"""
function p14_blob_sha(path)
    bytes = read(path)
    ctx   = SHA.SHA1_CTX()
    SHA.update!(ctx, Vector{UInt8}("blob $(length(bytes))\0"))
    SHA.update!(ctx, bytes)
    return bytes2hex(SHA.digest!(ctx))
end

"""
    p14_consts_sha() -> String

The sha256 of `spike/p14/consts.jl`, hex-encoded, recorded on every persisted Phase-14 artifact.
The Phase-14 analogue of `p13_consts_sha()` (`spike/p13/net.jl:524`), and it carries the same
audit meaning: the pre-registration forbids EDITING a constant, so a changed sha on a file that
should only ever have been APPENDED to is itself the signal.
"""
p14_consts_sha() = bytes2hex(SHA.sha256(read(joinpath(@__DIR__, "consts.jl"))))

"""
    _p14_load_report(path, what) -> Dict

Read one Phase-13 `.jld2` provenance report. A HARD PRECONDITION, unlike the deliberately inert
`_recorded_ood_threshold` path (`src/amortized/local_map.jl:93`) which returns `nothing` when its
gate report is absent: there, a missing operating point must leave a flag honestly NOT CHECKED;
here, a missing report means the threshold cannot be inherited at all, and D-07 forbids
substituting a literal for it. So this errors and never returns `nothing`.
"""
function _p14_load_report(path, what)
    isfile(path) || error("""
        p14_load_tau: no $what artifact at $path.

        D-07 requires tau to be INHERITED by LOADING the Phase-13 artifact with its provenance
        asserted. Substituting a literal here -- or defaulting, or skipping the check -- is
        exactly the divergence D-07 exists to prevent. Restore the artifact or do not run.
        """)
    try
        return JLD2.load(path)
    catch e
        e isa InterruptException && rethrow()
        error("""
            p14_load_tau: the $what artifact at $path exists but could not be read ($e).

            This is a hard precondition, not a skippable check: an unreadable provenance report
            is indistinguishable from an absent one for D-07's purposes.
            """)
    end
end

"Fetch `key` from a loaded Phase-13 report, naming the artifact and the key when it is absent."
function _p14_report_key(d, key, what, path)
    haskey(d, key) || error("p14_load_tau: the $what artifact at $path has no \"$key\" key; " *
                            "its layout is not the one D-07's provenance table was written against")
    return d[key]
end

# --- The loader --------------------------------------------------------------------------------

"""
    p14_load_tau(; net_path, probe_path, gate_path) -> NamedTuple

LOAD the Phase-13 three-way cut threshold and ASSERT its provenance. This is the only sanctioned
route by which Phase-14 code may obtain that number (D-07); no Phase-14 file writes it as a
literal, and `spike/test/test_p14_provenance.jl` greps the lane to keep it that way.

Order of operations:

  1. all three Phase-13 artifacts must exist -- otherwise `error`, never a default;
  2. the NET is loaded and is the SOURCE OF TRUTH, because the threshold travels with the net
     that was trained under it (14-RESEARCH section D.3);
  3. the probe and gate reports are loaded (hard preconditions, see `_p14_load_report`);
  4. the FOUR on-disk values are asserted equal, one assertion each so a failure names which
     site diverged;
  5. the sha table is pinned on BOTH sides against named, dated literals -- strictly stronger
     than pinning the two sides merely to each other, which a retrain plus an undisclosed edit
     would satisfy silently (`spike/p13/run_three_way_gate.jl:361-370`);
  6. the probe's THIRD `consts_sha256` is asserted DIVERGENT rather than tolerated.

RECORDED HONESTLY, NOT ASSERTED AWAY: `probe_repo_dirty_at_run` is `true` on disk. The Phase-13
tau probe ran on a dirty working tree. Phase 14 SURFACES that fact on every artifact it writes
and does NOT assert it false; a provenance field that only ever reports the convenient value is
not provenance.

Returns a NamedTuple carrying the threshold, the artifact paths it came from, every sha in the
table and the two measured AUC numbers, so a caller can embed the whole block without re-reading
any file. `four_way_agreed` is `true` by construction: this function cannot return with it false,
because a disagreement throws.
"""
function p14_load_tau(; net_path   = joinpath(P14_P13_DIR, "three_way_net.jld2"),
                        probe_path = joinpath(P14_P13_DIR, "tau_probe_report.jld2"),
                        gate_path  = joinpath(P14_P13_DIR, "three_way_gate_report.jld2"))

    # 1. Existence first, so a missing artifact is reported as a missing artifact rather than as
    #    a downstream key error.
    isfile(net_path) || error("""
        p14_load_tau: no three-way net artifact at $net_path.

        The net is the SOURCE OF TRUTH for the inherited threshold (D-07, 14-RESEARCH D.3): the
        cut travels with the net that was trained under it. D-07 forbids substituting a literal
        when it is absent, so this is a block on the phase, not a skip.
        """)

    # 2. The net.
    h = load_three_way(net_path)
    tau = h.tau
    tau === nothing && error("p14_load_tau: $net_path records no tau; it was persisted before " *
                             "the Phase-13 Tier-2 measurement existed and cannot be inherited from")

    # 3. The two provenance reports.
    probe = _p14_load_report(probe_path, "tau probe")
    gate  = _p14_load_report(gate_path, "three-way gate")

    probe_tau         = _p14_report_key(probe, "tau", "tau probe", probe_path)
    probe_bar         = _p14_report_key(probe, "tau_bar", "tau probe", probe_path)
    probe_repo_sha    = _p14_report_key(probe, "repo_sha", "tau probe", probe_path)
    probe_sim_sha     = _p14_report_key(probe, "simulator_sha", "tau probe", probe_path)
    probe_consts_sha  = _p14_report_key(probe, "consts_sha256", "tau probe", probe_path)
    probe_dirty       = _p14_report_key(probe, "repo_dirty_at_run", "tau probe", probe_path)
    gate_tau          = _p14_report_key(gate, "P13_TAU", "three-way gate", gate_path)
    gate_net_sha      = _p14_report_key(gate, "net_consts_sha", "three-way gate", gate_path)

    # 4. THE FOUR-WAY AGREEMENT. One assertion per equality, so the message names WHICH of the
    #    four sites moved. A single fused conjunction would report only that "something" drifted.
    @assert tau == probe_tau "D-07 breach: the net artifact and the tau probe report disagree about the measured three-way cut ($tau vs $probe_tau)"
    @assert tau == gate_tau "D-07 breach: the net artifact and the Phase-13 gate report disagree about the measured three-way cut ($tau vs $gate_tau)"
    @assert tau == P13_TAU "D-07 breach: the net artifact and the Phase-13 Tier-2 pre-registration disagree about the measured three-way cut ($tau vs $P13_TAU)"
    @assert tau == p13_tau() "D-07 breach: the net artifact and the Phase-13 Tier-2 accessor disagree about the measured three-way cut ($tau vs $(p13_tau()))"

    # 5. THE SHA TABLE, PINNED ON BOTH SIDES to named, dated literals. The two consts sha256
    #    values are EXPECTED to differ: the net was trained under the pre-amendment
    #    pre-registration and spike/p13/consts.jl now carries 13-D15-AMENDMENT.md. Training reads
    #    no P13_REAL_* constant, so the amended values are ones this net never consumed. Pinning
    #    each side to its own literal catches a retrain plus an undisclosed edit, which the older
    #    "the two sides agree with each other" form would have passed silently.
    @assert h.consts_sha == P13_CONSTS_SHA.pre_amendment "D-07 breach: the inherited net was NOT trained under the frozen Phase-13 pre-registration -- its recorded consts_sha is not P13_CONSTS_SHA.pre_amendment (spike/p13/net.jl:556-559)"
    @assert p13_consts_sha() == P13_CONSTS_SHA.post_amendment "D-07 breach: spike/p13/consts.jl is at NEITHER sha256 the Phase-13 amendment authorises -- a third value is drift, not an amendment"
    @assert gate_net_sha == h.consts_sha "D-07 breach: the Phase-13 gate report scored a net whose consts_sha differs from the net this loader read ($gate_net_sha vs $(h.consts_sha))"
    @assert probe_repo_sha == P13_TAU_PROBE_SHA "D-07 breach: the tau probe report records repo sha $probe_repo_sha, not the pinned P13_TAU_PROBE_SHA"
    @assert probe_sim_sha == P13_TAU_SIMULATOR_SHA "D-07 breach: the tau probe ran against simulator $probe_sim_sha, not the pinned P13_TAU_SIMULATOR_SHA"
    @assert probe_bar == P13_TAU_AUC "D-07 breach: the tau probe report's bar ($probe_bar) is not the frozen P13_TAU_AUC ($P13_TAU_AUC); the threshold was measured against a different bar than the one now on record"
    @assert P13_TAU_MEASURED_AUC >= P13_TAU_AUC "D-07 breach: the realized AUC at the inherited threshold does not clear the frozen P13_TAU_AUC bar"

    # 6. ASSERT THE DIVERGENCE. DO NOT TOLERATE IT, AND DO NOT WIDEN IT.
    #
    #    tau_probe_report.jld2 records a THIRD consts sha256, matching NEITHER literal of
    #    P13_CONSTS_SHA. That is CORRECT AND EXPECTED: the probe ran on 2026-07-27, BEFORE the
    #    Tier-2 tau block was appended to spike/p13/consts.jl. That ordering is precisely what
    #    spike/p13/consts.jl:842-846 cites as the auditable evidence that the threshold was
    #    MEASURED BEFORE IT WAS LOCKED -- the measurement exists in git history at a commit
    #    carrying no Tier-2 value, so no reader has to take the ordering on trust. A probe sha
    #    that MATCHED one of the two literals would mean the probe had seen the locked value.
    #
    #    REJECTED, and recorded here so a later reader knows it was considered rather than
    #    missed: widening this to a disjunction that ACCEPTS the third value -- the shape
    #    spike/p13/net.jl:549-554 already rejects for the net-side guard. Once such a
    #    disjunction exists it accepts ANY future drift silently, converting a real integrity
    #    check into decoration. That is forbidden. If this assertion ever fires, the probe
    #    artifact has been regenerated or replaced, and THAT is the finding -- not this line.
    @assert probe_consts_sha ∉ values(P13_CONSTS_SHA) "D-07 breach: the tau probe report's consts sha256 now MATCHES a pinned Phase-13 literal. It must diverge from both -- the probe ran before the Tier-2 append, which is the ordering evidence that the threshold was measured before it was locked. Do NOT repair this by widening the check to a disjunction; that destroys the integrity argument and is forbidden."

    hasproperty(h.meta, :phase11_sha) ||
        error("p14_load_tau: $net_path carries no meta.phase11_sha; the basis the net rode " *
              "cannot be pinned, so its provenance is incomplete")

    return (tau = tau,
            source = :three_way_net_jld2,
            net_path = net_path, probe_path = probe_path, gate_path = gate_path,
            net_consts_sha = h.consts_sha,
            file_consts_sha = p13_consts_sha(),
            probe_consts_sha = probe_consts_sha,
            probe_repo_sha = probe_repo_sha,
            probe_simulator_sha = probe_sim_sha,
            probe_repo_dirty_at_run = probe_dirty,
            tau_measured_auc = P13_TAU_MEASURED_AUC,
            tau_bar = P13_TAU_AUC,
            phase11_sha = h.meta.phase11_sha,
            four_way_agreed = true)
end

# --- The block every Phase-14 artifact embeds ---------------------------------------------------

"""
    p14_provenance_record(prov; extra = NamedTuple()) -> NamedTuple

The provenance block every Phase-14 `.jld2` artifact embeds, built from the NamedTuple
`p14_load_tau` returned.

Any caller-supplied `extra` (run-specific provenance such as `(; counter, n, passed)`) is MERGED
**UNDER** the required keys -- the `p13_calibration_meta` idiom (`spike/p13/result.jl:337-352`) --
so caller provenance is preserved but a required key can never be overwritten or omitted. A run
that "helpfully" supplied its own `tau` or its own `amendment` would otherwise be able to
overwrite the inherited ones, which is the single thing this record exists to prevent.

The amendment citation travels IN THE ARTIFACT, not only in a document: SC1 and SC2 are both
amended for Phase 14, and a number quoted without its amendment is a number quoted against a
criterion nobody applied.
"""
function p14_provenance_record(prov::NamedTuple; extra::NamedTuple = NamedTuple())
    required = merge(prov, (
        p14_consts_sha          = p14_consts_sha(),
        p14_consts_git_blob_sha = p14_blob_sha(P14_CONSTS_PATH),
        master_seed             = UInt64(P14_DEV_SEED),
        salt                    = UInt64(P14_SALT),
        generated               = string(Dates.now(Dates.UTC)) * "Z",
        amendment               = P14_AMENDMENT_NOTICE,
    ))
    return merge(extra, required)
end

# --- The D-01 lane guard, asserted WHILE the run happens ----------------------------------------

"""
    p14_assert_lane_clean(; repo_root = P14_REPO_ROOT, verbose = false) -> NamedTuple

STEP 0 of every Phase-14 runner (D-01, threat T-14-11). Asserts that `src/` and the frozen spike
environment are byte-unchanged against `HEAD`, WHILE the run happens rather than only in review:
a long reported run that started on a clean tree and finished on a dirty one is a repudiation
hole. Promoted from the `spike/p13/run_three_way_gate.jl:270-280` precondition, which is itself
the `spike/test/test_p13_result.jl:203-204` test promoted into a runner.

Degrades to an `@info` skip when `git` is unavailable, rather than failing a machine that simply
has no git -- the same degradation the Phase-13 test half uses. It never degrades on a git that
is present and reports a diff.
"""
function p14_assert_lane_clean(; repo_root = P14_REPO_ROOT, verbose::Bool = false)
    if Sys.which("git") === nothing
        @info "p14_assert_lane_clean: git unavailable -- skipping the executable D-01 lane check"
        return (checked = false, src_clean = missing, env_clean = missing)
    end
    src_clean = success(Cmd(`git diff --quiet HEAD -- src`; dir = repo_root))
    env_clean = success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                            dir = repo_root))
    @assert src_clean "D-01 decoupling breach: `git diff --quiet HEAD -- src` FAILED -- src/ is not byte-unchanged, so this run may not be reported"
    @assert env_clean "D-01 decoupling breach: spike/Project.toml or spike/Manifest.toml is modified -- the frozen spike environment moved under the run"
    verbose && println("STEP 0 (D-01): src/ byte-unchanged = $src_clean   " *
                       "spike Project/Manifest byte-unchanged = $env_clean")
    return (checked = true, src_clean = src_clean, env_clean = env_clean)
end
