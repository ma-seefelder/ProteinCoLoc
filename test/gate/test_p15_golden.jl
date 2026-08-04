#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/test_p15_golden.jl --- the assertions that make the D-09 golden a CONTROL rather
# than a rubber stamp.
#
# `ci_golden.jl` proves the golden MATCHES. That is the easy half. The hard half is that the
# golden cannot be quietly regenerated to make a red run green, and THAT is what this file
# asserts, from outside the script:
#
#   * `--rebless` REFUSES without a reason, and refuses WITHOUT WRITING (testset 3, checked by
#     byte-comparing the whole pre-registration file before and after);
#   * re-blessing APPENDS a new sentinel-keyed block and leaves the previous one byte-identical
#     (testset 4, on a temp COPY — never on the committed file);
#   * a MODIFIED `const` declaration line anywhere in `p15_consts.jl` fails the suite (testset 5,
#     the CONVENTIONS C-03 `-U0` declaration-line diff), and the guard is SHOWN TO FIRE on a
#     synthetic diff (testset 6, CONVENTIONS C-04 rule 3 — an assurance never shown non-empty is
#     not a test).
#
# ISOLATED MODULE IS MANDATORY, NOT STYLE. `ci_golden.jl` loads `p15_consts.jl`, which defines the
# SAME const names (`SBC_M`, `prod_seed`, …) that `runtests.jl` has already loaded from another
# gate pre-registration into the test module. Including it directly would hit `p15_consts.jl`'s
# `:SBC_M` guard, silently SKIP the load, and leave every bar below reading the wrong file — so
# `ci_golden.jl`'s own `GATE_CONSTS_VERSION == 15` check would then fire. `P15Golden` is a fresh
# module, so it gets the Phase-15 pre-registration and nothing else.
#
# WHAT THIS FILE DOES NOT NEED. Only testsets 1 and 2 need the shipped `grid_8` net. The re-bless
# discipline (3, 4) and the append-only guard (5, 6) are asserted against a SYNTHETIC measurement,
# so they run — and can fail — on a machine or a runner that cannot resolve the artifact at all.
# When the net is unresolvable, testsets 1 and 2 report an honest `@info` skip rather than a pass
# (the `:not_trained` pattern at `run_gate.jl:16-19`); a skipped calibration check must never look
# like a passing one.
#
# RUNS BOTH WAYS: standalone (`julia --project=. --threads=auto test/gate/test_p15_golden.jl`) and
# `include`d from `test/runtests.jl` (15-01 already wired the include list). Every path is built
# from `@__DIR__`, never from the process working directory.

using Test

module P15Golden
    include(joinpath(@__DIR__, "ci_golden.jl"))
end

"The pre-registration file the golden blocks live in. Read-only in every testset below."
const P15G_CONSTS_PATH = joinpath(@__DIR__, "p15_consts.jl")

"The repository root, resolved from this file rather than from the process working directory."
const P15G_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# =============================================================================================
# Net availability — probed ONCE, honestly
# =============================================================================================
# `estimator_for(8)` registers the bundle process-locally on first success, so this probe costs
# the resolution once and every later call is a dict lookup.

const P15G_NET_ERROR = Ref{Any}(nothing)

const P15G_NET_OK = try
    P15Golden.p15_golden_net()
    true
catch err
    P15G_NET_ERROR[] = err
    false
end

if !P15G_NET_OK
    @info """
    test_p15_golden: the shipped grid_8 bundle could NOT be resolved, so the two MEASUREMENT
    testsets are SKIPPED (not passed). `artifacts/` is git-ignored and `Artifacts.toml` points
    grid_8 at a GitHub release asset; on a clone or a runner without that release the net is
    genuinely unavailable and the golden cannot be evaluated. The re-bless discipline and the
    append-only guard below still run in full.""" exception = P15G_NET_ERROR[]
end

"Measured once and reused: two independent measurements for the determinism testset."
const P15G_MEASURED = P15G_NET_OK ?
    (P15Golden.p15_golden_measure(), P15Golden.p15_golden_measure()) : nothing

# =============================================================================================
# A synthetic measurement — so the DISCIPLINE testsets need no net
# =============================================================================================

"""
A stand-in `p15_golden_measure()` result with obviously-fake values. Used only to exercise block
GENERATION and APPENDING; it is never compared against a committed golden and never written to
the committed file.
"""
const P15G_FAKE = (artifact_hash      = "0123456789abcdef0123456789abcdef01234567",
                   net_content_digest = "11111111111111111111",
                   resolved_via       = :test_fixture,
                   rank_digest        = "22222222222222222222",
                   ece                = collect(range(0.1; step = 0.01, length = 8)),
                   labels             = ["c$i" for i in 1:8],
                   julia_version      = string(VERSION))

# =============================================================================================
# The CONVENTIONS C-03 append-only predicate — a PURE function, so it can be shown to fire
# =============================================================================================

"""
    p15g_appendonly_violations(diff) -> Vector{String}

The changed-DECLARATION-line rule of CONVENTIONS C-03, as a pure function of a `git diff -U0`
string so testset 6 can feed it a synthetic diff and SEE it report a violation.

Keeps only added/removed lines (`^[+-]`), drops the `+++`/`---` file headers, narrows to lines
whose remainder after optional whitespace begins with `const `, and returns every REMOVED one.
`-U0` is what makes this correct: with default context, a `const` line sitting three lines away
from an unrelated edit would be counted even though nothing about it changed, and a criterion that
fails correct work trains an executor to disable it.
"""
function p15g_appendonly_violations(diff::AbstractString)
    out = String[]
    for line in split(diff, '\n')
        isempty(line) && continue
        c = line[1]
        (c == '+' || c == '-') || continue
        length(line) >= 2 && (line[2] == '+' || line[2] == '-') && continue   # +++ / --- headers
        rest = lstrip(line[2:end])
        startswith(rest, "const ") || continue
        c == '-' && push!(out, line)
    end
    return out
end

"""
    p15g_base_ref() -> Union{String, Nothing}

The first base ref that actually resolves in THIS checkout, or `nothing`. A shallow CI checkout
has no `origin/main`, and a guard that fails because the base ref is missing is a guard that gets
disabled — so an unavailable base is an honest skip, never a failure.
"""
function p15g_base_ref()
    for r in ("origin/main", "origin/master", "main", "master")
        try
            s = read(Cmd(`git rev-parse --verify -q $r`; dir = P15G_REPO_ROOT), String)
            isempty(strip(s)) || return r
        catch
        end
    end
    return nothing
end

# =============================================================================================
@testset "P15 D-09 golden (15-06)" begin

    @testset "1. golden is deterministic in-process" begin
        if P15G_MEASURED === nothing
            @info "SKIPPED (1): no resolvable grid_8 net — see the @info above."
            @test_skip false
        else
            a, b = P15G_MEASURED
            # The rank digest is the strong claim: identical draw order, summaries and flow pass.
            @test a.rank_digest == b.rank_digest
            @test a.artifact_hash == b.artifact_hash
            @test a.net_content_digest == b.net_content_digest
            @test length(a.ece) == 8
            @test a.ece == b.ece          # ELEMENTWISE identical, not merely within tolerance
            @test a.ranks == b.ranks
            @test a.global_seed == b.global_seed
        end
    end

    @testset "2. golden matches the committed values" begin
        if P15G_MEASURED === nothing
            @info "SKIPPED (2): no resolvable grid_8 net — see the @info above."
            @test_skip false
        else
            golden = P15Golden.p15_golden_blessed()
            @test golden !== nothing
            r = P15Golden.p15_golden_check(P15G_MEASURED[1], golden)
            # Print every failure before asserting, so a red run is diagnosable from the log alone.
            for f in r.failures
                @warn "golden mismatch" detail = f
            end
            @test r.ok
            @test isempty(r.failures)
        end
    end

    @testset "3. rebless refuses without a reason (D-09)" begin
        before = read(P15G_CONSTS_PATH)                       # RAW BYTES, not a digest
        script = joinpath(@__DIR__, "ci_golden.jl")
        proj   = Base.active_project()
        cmd    = `$(Base.julia_cmd()) --project=$(proj) --startup-file=no $script --rebless`
        # The env var IS a legitimate way to supply a reason, so it must be absent for the
        # refusal to be the thing under test.
        cmd = addenv(cmd, "P15_REBLESS_REASON" => nothing)
        p   = run(pipeline(ignorestatus(cmd); stdout = devnull, stderr = devnull))
        @test !success(p)
        @test p.exitcode == 3
        after = read(P15G_CONSTS_PATH)
        # Byte equality is STRICTLY STRONGER than the digest equality D-09 asks for, and needs no
        # crypto: the `SHA` stdlib is not a [dep] and importing it aborts the suite under
        # `Pkg.test()`'s sandbox (15-01, deviation 2).
        @test after == before
        @test hash(after) == hash(before)
    end

    @testset "4. rebless appends, never edits (D-09)" begin
        mktempdir() do dir
            copy = joinpath(dir, "p15_consts.jl")
            cp(P15G_CONSTS_PATH, copy)
            before_lines = readlines(copy)
            before_bytes = read(copy)

            nextblock = P15Golden.p15_golden_next_block()
            @test nextblock >= 2                       # block 1 is already committed
            text = P15Golden.p15_golden_block_text(P15G_FAKE; block = nextblock,
                                                   reason = "testset 4 fixture — not a real bless")
            P15Golden.p15_golden_append!(copy, text)

            after_lines = readlines(copy)
            after_bytes = read(copy)

            # It GREW …
            @test length(after_bytes) > length(before_bytes)
            @test length(after_lines) > length(before_lines)
            # … and every pre-existing line — including every line of the committed golden block —
            # survives byte-identically, in order, as a strict PREFIX of the new file.
            @test after_lines[1:length(before_lines)] == before_lines
            @test after_bytes[1:length(before_bytes)] == before_bytes
            # The new block carries a DISTINCT sentinel; the previous one is still declared.
            newtext = join(after_lines[(length(before_lines) + 1):end], "\n")
            @test occursin("P15_GOLDEN_RANK_DIGEST_V$(nextblock)", newtext)
            @test any(l -> occursin("const P15_GOLDEN_RANK_DIGEST ", l), before_lines)
            @test any(l -> occursin("const P15_GOLDEN_RANK_DIGEST ", l), after_lines)
            # The reason is echoed into the appended block, so the audit trail is in the file too.
            @test occursin("testset 4 fixture", newtext)
            # And the COMMITTED file was not touched by any of this.
            @test read(P15G_CONSTS_PATH) == before_bytes
        end
    end

    @testset "5. append-only guard (CONVENTIONS C-03)" begin
        base = p15g_base_ref()
        if base === nothing
            @info "SKIPPED (5): no base ref resolves in this checkout (a shallow CI clone has " *
                  "no origin/main). The guard is exercised on a synthetic diff in testset 6."
            @test_skip false
        else
            diff = try
                read(Cmd(`git diff -U0 $base -- test/gate/p15_consts.jl`;
                         dir = P15G_REPO_ROOT), String)
            catch err
                @info "SKIPPED (5): git diff unavailable" exception = err
                nothing
            end
            if diff === nothing
                @test_skip false
            else
                v = p15g_appendonly_violations(diff)
                for l in v
                    @warn "PRE-REGISTRATION BREACH: a `const` declaration in p15_consts.jl was " *
                          "MODIFIED rather than appended. Constants are only ever APPENDED, so " *
                          "the git history of that file IS the audit trail (D-09)." line = l
                end
                @test isempty(v)
            end
        end
    end

    @testset "6. the guard fires (CONVENTIONS C-04 rule 3)" begin
        # A synthetic `-U0` diff carrying one MODIFIED declaration line. An assurance never shown
        # non-empty is an assurance, not a test.
        fired = """
        --- a/test/gate/p15_consts.jl
        +++ b/test/gate/p15_consts.jl
        @@ -1004 +1004 @@
        -    const P15_GOLDEN_ECE = [0.1, 0.2]
        +    const P15_GOLDEN_ECE = [0.3, 0.4]
        """
        v = p15g_appendonly_violations(fired)
        @test length(v) == 1
        @test occursin("P15_GOLDEN_ECE", v[1])

        # A pure APPEND must NOT fire — otherwise the guard would fail correct work.
        clean = """
        --- a/test/gate/p15_consts.jl
        +++ b/test/gate/p15_consts.jl
        @@ -1142,0 +1143,2 @@
        +    const P15_GOLDEN_RANK_DIGEST_V2 = "abc"
        +    const P15_GOLDEN_ECE_V2 = [0.3, 0.4]
        """
        @test isempty(p15g_appendonly_violations(clean))

        # And the `---`/`+++` headers must never be mistaken for removed declarations, nor a
        # removed NON-declaration line (a comment, a blank) counted as a breach.
        headers_only = """
        --- a/test/gate/p15_consts.jl
        +++ b/test/gate/p15_consts.jl
        @@ -5 +5 @@
        -# a comment mentioning const P15_GOLDEN_ECE
        +# a different comment
        """
        @test isempty(p15g_appendonly_violations(headers_only))
    end
end
