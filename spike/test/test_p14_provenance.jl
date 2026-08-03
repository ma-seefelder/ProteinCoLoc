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

# spike/test/test_p14_provenance.jl --- the D-07 threshold-provenance gate.
#
# WHAT THIS FILE IS FOR. D-07 says the Phase-13 three-way cut is INHERITED: Phase 14 LOADS it
# from the Phase-13 artifact with its provenance asserted, and never copies the literal and never
# re-derives it. "Copying the value is how two numbers silently diverge; re-deriving guarantees a
# second, different cut on the same hypothesis space." That decision is only worth anything if a
# future divergence is LOUD, so this file turns the whole D-07 table into a running test that
# finishes in seconds:
#
#   1. the threshold is LOADED, and no Phase-14 source writes it as a literal;
#   2. the FOUR on-disk sites agree -- checked by loading all four INDEPENDENTLY here, so this
#      file is a check ON the loader rather than a restatement OF it;
#   3. the sha table is pinned on BOTH sides to named, dated literals;
#   4. the tau probe's THIRD consts_sha256 is asserted DIVERGENT -- see the long note in
#      testset 4; widening that to a disjunction is itself the failure;
#   5. the probe's dirty-tree fact is SURFACED, not asserted away;
#   6. the Phase-13 lane and src/ are byte-unchanged (D-01).
#
# Mirrors test_p13_consts.jl / test_p13_result.jl: license header -> using Test -> guarded
# include of the unit under test -> sources read ONCE at module level -> one outer @testset.
# No training, no simulation, no figure call.
#
# THIS FILE IS RUN PER-FILE, NEVER THROUGH runtests.jl, which exits 1 early at the Phase-4
# SPEEDUP_GATE and masks every later include block:
#
#     julia --project=spike spike/test/test_p14_provenance.jl

using Test
using JLD2

# Unit under test. Guarded so a re-include under a runner is a silent no-op; the include
# transitively brings the Phase-13 pre-registration (P13_TAU, P13_CONSTS_SHA, p13_tau(),
# p13_consts_sha()) and the Phase-14 Tier-1 block into this flat top-level namespace.
isdefined(@__MODULE__, :P14_PROVENANCE_LOADED) ||
    include(joinpath(@__DIR__, "..", "p14", "provenance.jl"))

# The loader is called ONCE at module level (the test_bf.jl:67 read-once idiom): every testset
# below reads the same returned block, so a testset cannot accidentally pass against a second,
# differently-loaded copy.
const P14_PROV = p14_load_tau()

# The EXECUTABLE source of the loader, with whole-line comments stripped BEFORE any grep. The
# file necessarily NAMES the forbidden widened form in order to forbid it, so an un-stripped scan
# would be red by construction -- the same reason test_p13_consts.jl:46-48 strips first.
const P14_PROV_SRC  = read(joinpath(@__DIR__, "..", "p14", "provenance.jl"), String)
const P14_PROV_CODE = join(filter(l -> !startswith(strip(l), "#"),
                                  split(P14_PROV_SRC, '\n')), '\n')

# Every Phase-14 source that exists RIGHT NOW. Globbed rather than enumerated so a file added in
# a later wave is scanned automatically instead of needing this list edited -- an unedited list
# is how a new file quietly escapes a lane-wide guard.
const P14_LANE_DIR   = normpath(joinpath(@__DIR__, "..", "p14"))
const P14_LANE_FILES = sort(filter(f -> endswith(f, ".jl"), readdir(P14_LANE_DIR)))

# The forbidden shape: any name containing `tau`, in any casing, assigned the Phase-13 literal.
# Built as a Regex over the COMMENT-STRIPPED source of each lane file. A legitimate
# `tau = prov.tau` is not matched; only a hardcoded copy is.
const P14_TAU_LITERAL_RE = r"\w*[Tt][Aa][Uu]\w*\s*=\s*0\.15"

"Comment-stripped source of one `spike/p14/` file."
_p14_lane_code(f) = join(filter(l -> !startswith(strip(l), "#"),
                                split(read(joinpath(P14_LANE_DIR, f), String), '\n')), '\n')

@testset "P14 tau provenance (D-07)" verbose = true begin

    @testset "tau is loaded, not written here" begin
        # The value, its declared source, and the loader's own statement that all four sites
        # agreed. `four_way_agreed` cannot be false: a disagreement throws inside the loader.
        @test P14_PROV.tau === 0.15
        @test P14_PROV.source === :three_way_net_jld2
        @test P14_PROV.four_way_agreed

        # NO PHASE-14 SOURCE MAY WRITE THE THRESHOLD AS A LITERAL. This is the executable half of
        # D-07: the decision is only real if a copied literal FAILS rather than merely being
        # discouraged in a comment.
        @test !isempty(P14_LANE_FILES)
        for f in P14_LANE_FILES
            @test !occursin(P14_TAU_LITERAL_RE, _p14_lane_code(f))
        end
        # spike/p14/consts.jl additionally has no Phase-14 tau binding at all -- the Tier-2 slot
        # is INHERITED, never appended (needle assembled at run time so this line does not itself
        # plant the forbidden name in a scanned source; the test_p12_decoupling.jl:34-45 idiom).
        @test !isdefined(@__MODULE__, Symbol("P14_", "TAU"))

        # ARMS ITSELF LATER. These Phase-14 modules land in later waves; when they do, the glob
        # above picks them up automatically. Named here so a reader can see the guard is
        # lane-wide by design rather than by accident.
        for later in ("result.jl", "fuse.jl", "fdr.jl", "conformal.jl", "posterior.jl",
                      "decide.jl")
            if isfile(joinpath(P14_LANE_DIR, later))
                @test !occursin(P14_TAU_LITERAL_RE, _p14_lane_code(later))
            else
                @info "spike/p14/$later does not exist yet — the literal guard arms itself when it lands"
            end
        end
    end

    @testset "the four on-disk taus agree" begin
        # LOADED INDEPENDENTLY HERE, not read back off P14_PROV. Asserting the loader's own
        # return value against itself would be a tautology; the point of D-07 is that the four
        # SITES agree, so this testset goes to the four sites.
        h     = load_three_way(P14_PROV.net_path)
        probe = JLD2.load(P14_PROV.probe_path)
        gate  = JLD2.load(P14_PROV.gate_path)

        @test h.tau == P14_PROV.tau
        @test probe["tau"] == P14_PROV.tau
        @test gate["P13_TAU"] == P14_PROV.tau
        @test P13_TAU == P14_PROV.tau
        @test p13_tau() == P14_PROV.tau
        # The measured value is a point ON the frozen Phase-13 grid, never interpolated: the
        # inherited number is the one the reading rule produced.
        @test P14_PROV.tau in P13_TAU_DELTA_GRID
        # And it cleared the bar it was measured against, which is the reason it exists.
        @test P14_PROV.tau_bar == P13_TAU_AUC
        @test probe["tau_bar"] == P14_PROV.tau_bar
        @test P14_PROV.tau_measured_auc >= P14_PROV.tau_bar
    end

    @testset "the sha table is pinned on both sides" begin
        # THE TWO consts sha256 VALUES ARE EXPECTED TO DIFFER. The net was trained under the
        # pre-amendment pre-registration; spike/p13/consts.jl now carries 13-D15-AMENDMENT.md.
        # Pinning EACH SIDE to its own named, dated literal is strictly stronger than the older
        # "the two sides agree with each other" form, which a retrain plus an undisclosed edit
        # would have satisfied silently (spike/p13/run_three_way_gate.jl:361-370).
        @test P14_PROV.net_consts_sha == P13_CONSTS_SHA.pre_amendment
        @test P14_PROV.file_consts_sha == P13_CONSTS_SHA.post_amendment
        @test P14_PROV.net_consts_sha != P14_PROV.file_consts_sha
        @test occursin(r"^[0-9a-f]{64}$", P14_PROV.net_consts_sha)
        @test occursin(r"^[0-9a-f]{64}$", P14_PROV.file_consts_sha)

        # The tree the probe ran on and the simulator it measured against are both PINNED as
        # full 40-character git object names rather than described in prose.
        @test P14_PROV.probe_repo_sha == P13_TAU_PROBE_SHA
        @test P14_PROV.probe_simulator_sha == P13_TAU_SIMULATOR_SHA
        @test occursin(r"^[0-9a-f]{40}$", P14_PROV.probe_repo_sha)
        @test occursin(r"^[0-9a-f]{40}$", P14_PROV.probe_simulator_sha)

        # The Phase-11 basis the net rode, pinned so a re-based net cannot pass unnoticed.
        @test P14_PROV.phase11_sha ==
              "659cda09112c4f20d82e7fd1bc6dc22d57f3853ee245f3c6f781407eb68eb571"

        # The gate report scored the SAME net this loader read.
        @test JLD2.load(P14_PROV.gate_path)["net_consts_sha"] == P14_PROV.net_consts_sha
    end

    @testset "the probe's third sha is asserted DIVERGENT (the trap)" begin
        # THIS DIVERGENCE IS CORRECT, AND IT IS THE POINT.
        #
        # tau_probe_report.jld2 records a THIRD consts sha256 matching NEITHER literal of
        # P13_CONSTS_SHA, because the probe ran on 2026-07-27 BEFORE the Tier-2 tau block was
        # appended to spike/p13/consts.jl. That ordering is exactly what spike/p13/consts.jl:842-846
        # cites as the auditable evidence that the threshold was MEASURED BEFORE IT WAS LOCKED.
        # A probe sha that MATCHED a pinned literal would mean the probe had seen the locked
        # value -- which is why the assertion runs in this direction.
        #
        # WIDENING THIS TO A DISJUNCTION IS ITSELF THE FAILURE, not a repair. The rejected form
        # is recorded at spike/p13/net.jl:549-554 for the net-side guard, for the same reason:
        # once such a disjunction exists it accepts ANY future drift silently, converting a real
        # integrity check into decoration.
        @test !(P14_PROV.probe_consts_sha in values(P13_CONSTS_SHA))
        @test P14_PROV.probe_consts_sha ==
              "640ec90edb80ef62eacbbf64c7b40cb9ca87c59d28a3a7438ef9d382b0aca4cf"
        @test occursin(r"^[0-9a-f]{64}$", P14_PROV.probe_consts_sha)

        # THE LOADER MUST ASSERT IT IN THE SAME DIRECTION. Grepped over the comment-stripped
        # source, so the loader's own explanatory prose -- which names the forbidden form in
        # order to forbid it -- cannot satisfy or defeat this check.
        @test occursin("∉ values(P13_CONSTS_SHA)", P14_PROV_CODE)
        @test !occursin(" in values(P13_CONSTS_SHA)", P14_PROV_CODE)
        @test !occursin("∈ values(P13_CONSTS_SHA)", P14_PROV_CODE)
        # ...and its CONDITION carries no disjunction. Taken over the part of the line before the
        # assertion message, so the prose in the message cannot trip or excuse this.
        _sha_lines = filter(l -> occursin("values(P13_CONSTS_SHA)", l),
                            split(P14_PROV_CODE, '\n'))
        @test length(_sha_lines) == 1
        @test !occursin("||", first(split(_sha_lines[1], '"')))
    end

    @testset "the dirty-tree fact is surfaced, not asserted away" begin
        # The Phase-13 tau probe ran on a DIRTY working tree, and Phase 14 reports that on every
        # artifact it writes rather than pretending otherwise. A provenance field that only ever
        # carries the convenient value is not provenance. 14-RESEARCH D.3: "Surface it in the
        # Phase-14 report; do not assert it false."
        @test P14_PROV.probe_repo_dirty_at_run === true
        @test JLD2.load(P14_PROV.probe_path)["repo_dirty_at_run"] === true
        # And it travels: the embedded block a Phase-14 artifact carries keeps the field.
        @test p14_provenance_record(P14_PROV).probe_repo_dirty_at_run === true
    end

    @testset "the provenance record cannot be overwritten by a caller" begin
        # The p13_calibration_meta idiom (spike/p13/result.jl:337-352): caller extras are merged
        # UNDER the required keys. A run that "helpfully" supplied its own tau or its own
        # amendment would otherwise overwrite the inherited ones, which is the single thing this
        # record exists to prevent.
        rec = p14_provenance_record(P14_PROV;
                                    extra = (tau = -1.0, amendment = "none", run_label = :smoke))
        @test rec.tau === P14_PROV.tau
        @test occursin("AMENDED by D-02", rec.amendment)
        @test occursin("AMENDED by D-05", rec.amendment)
        @test rec.run_label === :smoke          # caller provenance survives
        @test rec.master_seed == UInt64(P14_DEV_SEED)
        @test rec.salt == UInt64(P14_SALT)
        @test occursin(r"^[0-9a-f]{64}$", rec.p14_consts_sha)
        @test occursin(r"^[0-9a-f]{40}$", rec.p14_consts_git_blob_sha)
        @test endswith(rec.generated, "Z")
    end

    @testset "a missing artifact is a block, never a default" begin
        # D-07 forbids substituting a literal when the artifact is absent, so the loader must
        # THROW rather than fall back. Non-existent paths only -- nothing is written or removed.
        @test_throws Exception p14_load_tau(net_path = joinpath(P14_LANE_DIR, "no_such_net.jld2"))
        @test_throws Exception p14_load_tau(probe_path = joinpath(P14_LANE_DIR, "no_such_probe.jld2"))
        @test_throws Exception p14_load_tau(gate_path = joinpath(P14_LANE_DIR, "no_such_gate.jld2"))
    end

    @testset "the Phase-13 lane is read-only (D-01)" begin
        # Phase 14 CONSUMES Phase 13's artifacts and must not write to them; and src/ stays
        # byte-frozen for the whole spike. This is what makes "read-only" checkable rather than
        # asserted in prose. Guarded: a sandbox without git skips the executable half rather
        # than failing for an unrelated reason. `isdir(.git) || isfile(.git)` because inside a
        # linked worktree `.git` is a FILE (test_p12_decoupling.jl:64).
        @test P14_SRC_UNTOUCHED == true
        _has_git = Sys.which("git") !== nothing &&
                   (isdir(joinpath(P14_REPO_ROOT, ".git")) ||
                    isfile(joinpath(P14_REPO_ROOT, ".git")))
        if _has_git
            @test success(Cmd(`git diff --quiet HEAD -- spike/p13`; dir = P14_REPO_ROOT))
            @test isempty(readchomp(Cmd(`git status --porcelain -- spike/p13`;
                                        dir = P14_REPO_ROOT)))
            @test success(Cmd(`git diff --quiet HEAD -- src`; dir = P14_REPO_ROOT))
            @test isempty(readchomp(Cmd(`git status --porcelain -- src`; dir = P14_REPO_ROOT)))
            # The frozen spike environment, the constraint test_p12_decoupling.jl:176-178 owns.
            @test success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                              dir = P14_REPO_ROOT))
            # The runner-side guard promoted from spike/p13/run_three_way_gate.jl:270-280.
            _lane = p14_assert_lane_clean()
            @test _lane.checked
            @test _lane.src_clean && _lane.env_clean
        else
            @info "git unavailable — skipping the executable read-only assertions"
        end
    end

    @testset "P14 provenance ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
