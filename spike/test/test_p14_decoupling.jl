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

# spike/test/test_p14_decoupling.jl --- Phase 14's hard constraints, made EXECUTABLE.
#
# WHY THIS FILE EXISTS, AND WHY IT EXISTS *NOW*. 14-CONTEXT.md's "Constraints that bind this
# phase" section, CLAUDE.md's decoupling rule and D-01/D-02/D-03a/D-06/D-07 are all stated in
# PROSE, across several documents. Prose does not fail a build. Every rule below is an assertion
# that turns red the moment it stops being true. It runs in WAVE 1, before the modules it guards
# exist, so a violation fails on the commit that INTRODUCES it rather than at review. An assertion
# added at phase close proves only that nothing broke SINCE the assertion was added.
#
# IT RESTATES test_p12_decoupling.jl RATHER THAN EXTENDING IT. That file's own banner records the
# house rule -- "each phase states its own claim in its own file" (:187-188) -- and its scanned
# surface is `spike/p12/**` plus `test_p12_*.jl` (:78-88), which CANNOT see a Phase-14 file. A
# Phase-14 source is therefore invisible to the Phase-12 guard; that gap is exactly what this file
# closes. `test_p12_decoupling.jl` is NOT edited by this plan.
#
# THE SOURCE SCAN SCANS *ITSELF*, AND THAT IS WHY THE NEEDLES ARE BUILT BY CONCATENATION.
# This file matches `spike/test/test_p14_*.jl`, so it is one of the Phase-14 sources the scans
# below enumerate. A naive token ban is therefore RED BY CONSTRUCTION: the file has to name the
# thing it forbids in order to forbid it. Two different answers are used for two different halves
# of the problem:
#   * the ALLOWLIST (`_P14_CORPUS_METADATA_ALLOWED`) exempts the files that legitimately mention
#     corpus METADATA -- see testset 4's comment for why each of them must;
#   * every needle that would still match this file even under the allowlist is assembled at run
#     time from fragments, so the contiguous literal never appears in this source at all. That is
#     what lets the scan INCLUDE this file rather than exempt it.
# Neither trick weakens anything: the needle a run searches for is the full token either way.
#
# THIS FILE IS RUN PER-FILE, NEVER THROUGH runtests.jl. The aggregate suite exits 1 early at the
# Phase-4 SPEEDUP_GATE and masks every later include block, so a green aggregate run would be no
# evidence at all:
#
#     julia --project=spike spike/test/test_p14_decoupling.jl
#
# DECOUPLING (CLAUDE.md): spike-local. Shells out to `git` READ-ONLY, reads Phase-14 source TEXT
# and queries the active Pkg environment. Writes nothing, adds no package, touches no `src/` file,
# opens no image and loads no artifact.

using Test
using Pkg

# Resolved statically from this file's own location: from `spike/test/` the repo root is TWO
# levels up. Deliberately NOT obtained by including `spike/p14/provenance.jl` -- this file must be
# runnable even if the loader is broken, since "the loader is broken" is one of the things a
# decoupling failure looks like. The value is identical to `P14_REPO_ROOT` there, so a re-include
# in one module is a same-value redefinition rather than a conflict.
const P14_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# The ONE environmental question this file cannot pre-register: is `git` usable here at all? Same
# predicate as `test_p12_decoupling.jl:65-67` -- a source tree that is not a checkout skips the
# executable half rather than failing for an unrelated reason. `isdir(...) || isfile(...)` because
# inside a linked worktree `.git` is a FILE, not a directory.
const _HAS_GIT = Sys.which("git") !== nothing &&
                 (isdir(joinpath(P14_REPO_ROOT, ".git")) ||
                  isfile(joinpath(P14_REPO_ROOT, ".git")))

# Source text with whole-line comments stripped, so a path MENTIONED IN PROSE cannot be mistaken
# for a call. Ported from `spike/test/test_p12_decoupling.jl:71-73`.
_strip_comment_lines(src::AbstractString) =
    join(filter(l -> !startswith(strip(l), "#"), split(src, '\n')), '\n')

# --- The Phase-14 source surface ----------------------------------------------------------------
# Every file this phase may add, as (directory, filename pattern). A GLOB, never a hand-maintained
# list: an enumerated list is precisely how a new file escapes a lane-wide guard (14-02-SUMMARY
# "Assumption Drift" 2). Files that do not exist yet are simply absent from the enumeration -- this
# guard runs at wave 1 and six of the seven modules land later -- so the scan GROWS with the phase
# instead of asserting against a frozen list that goes stale.
const _P14_SOURCE_PATTERNS = (("spike/p14", r"^.*\.jl$"),
                              ("spike/test", r"^test_p14_.*\.jl$"))

"Repo-relative, forward-slashed paths of every Phase-14 source that EXISTS right now."
function _p14_sources()
    out = String[]
    for (rel, pat) in _P14_SOURCE_PATTERNS
        d = joinpath(P14_REPO_ROOT, rel)
        isdir(d) || continue
        for f in readdir(d)
            occursin(pat, f) || continue
            p = rel * "/" * f
            p in out || push!(out, p)
        end
    end
    return sort(out)
end

_p14_body(relpath::AbstractString) =
    _strip_comment_lines(read(joinpath(P14_REPO_ROOT, relpath), String))

# --- The corpus needles ---------------------------------------------------------------------------
# THE EXEMPTION SURFACE, DECLARED AS DATA SO A READER SEES ALL OF IT AT ONCE. It is deliberately
# short. ADDING AN ENTRY HERE IS A DECISION, NOT A FIX -- if a new Phase-14 file needs to say
# `corpus`, the question to answer first is why.
const _P14_CORPUS_METADATA_ALLOWED = ("spike/p14/consts.jl",
                                      "spike/p14/run_p14_real_images.jl",
                                      "spike/test/test_p14_consts.jl",
                                      "spike/test/test_p14_decoupling.jl")

# Banned outright in every NON-allowlisted Phase-14 source.
const _P14_CORPUS_TOKENS = ("corpus", "sealed_holdout", "manifest.csv", "CORPUS_MASTER_SEED")
const _CORPUS_PATH_TOKEN = "corpus"
# Image-loading calls. `build_mci` is the spike's own MultiChannelImage constructor
# (`spike/contract.jl`), so it is the third way a byte of image data enters this phase.
const _IMAGE_LOADERS = ("load(", "FileIO.load", "build_mci")

# ASSEMBLED AT RUN TIME (see the banner): the contiguous literals must not appear in this source,
# because this source is itself scanned. The value searched for is the full token.
const _SEALED_OPENER            = "open_" * "sealed_holdout"
const _SEALED_IMAGE_DIR_LITERAL = "cor" * "pus/data"
const _SEALED_IMAGE_DIR_SYMBOL  = "CORPUS_" * "DATA_DIR"

# --- The write-path needles -----------------------------------------------------------------------
# Phase 14 READS Phase 13's artifacts and pool and WRITES only under its own roots. Both checks
# below are conjunctions (a path token AND a write call) precisely so that a read stays legal:
# banning the path outright would forbid `h.meta.pool_dir`, which the OOD null fit requires.
const _ARTIFACT_WRITERS = ("jldsave", "mv(")
const _CACHE_WRITERS    = ("jldsave", "mkpath")
# Assembled from fragments for the same reason as the sealed needles: this file is scanned too.
const _P13_CACHE_LITERAL = "\"cache\", " * "\"p13\""

"""
Lines of a comment-stripped body that would read a sealed image: the sole reachable opener
(`corpus/manifest.jl` D-09), or an image-loading call on a path expression naming the corpus.
Returns the offending LINES so a failure shows the code, not just a boolean.
"""
function _sealed_image_load_hits(body::AbstractString)
    hits = String[]
    for ln in split(body, '\n')
        s = strip(ln)
        isempty(s) && continue
        if occursin(_SEALED_OPENER, s)
            push!(hits, s)
            continue
        end
        occursin(_CORPUS_PATH_TOKEN, s) || continue
        any(t -> occursin(t, s), _IMAGE_LOADERS) && push!(hits, s)
    end
    return hits
end

@testset "P14 decoupling: src/, deps, the sealed holdout and the phase's own greps (14-03)" verbose = true begin

    @testset "src/ is byte-unchanged (CLAUDE.md hard constraint, D-01)" begin
        # D-01 SAYS THE DECISION LAYER IS BUILT WITH A `src/`-SHAPED SIGNATURE. "src-shaped" is
        # NOT "in src", and this testset is what makes that distinction EXECUTABLE rather than a
        # sentence in a context document: `decide_coloc` may take the arguments it would take in
        # `src/` and return an `AbstractColocResult` subtype, and it may still not be written
        # there. 14-CONTEXT D-01: "the decision layer is NOT shipped in v2.0. Say so plainly in
        # any report; do not let 'src-shaped' be read as 'in src'."
        if _HAS_GIT
            # Whole tree first, then the four files 14-PATTERNS names as READ-ONLY references, so
            # a failure NAMES the file instead of saying "something under src/".
            @test success(Cmd(`git diff --quiet HEAD -- src`; dir = P14_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/results.jl`; dir = P14_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/amortized/ood.jl`; dir = P14_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/amortized/local_map.jl`; dir = P14_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/amortized/api.jl`; dir = P14_REPO_ROOT))
            # THE UNTRACKED CASE, WHICH `git diff HEAD` CANNOT SEE. A NEW file dropped under
            # `src/` is invisible to a diff against HEAD and would pass every assertion above
            # while being exactly the productionization-by-stealth this constraint forbids.
            @test isempty(readchomp(Cmd(`git status --porcelain -- src`; dir = P14_REPO_ROOT)))
        else
            @info "git unavailable — skipping the executable src/ byte-equality assertions"
        end
    end

    @testset "the spike environment is byte-frozen (D-02)" begin
        # 14-PATTERNS §0: every Phase-14 pattern was checked against the frozen sixteen, and the
        # stdlib reach (Test, Pkg, Statistics, Dates, SHA, LinearAlgebra) rides the default
        # LOAD_PATH `@stdlib` entry and needs no Project.toml entry. So ANY movement in these two
        # files means a `Pkg.add`/`Pkg.resolve` happened -- which is the thing the frozen manifest
        # exists to prevent, and the concrete reason D-02 declines `ConformalPrediction.jl`.
        if _HAS_GIT
            @test success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                              dir = P14_REPO_ROOT))
            @test isempty(readchomp(Cmd(`git status --porcelain -- spike/Project.toml spike/Manifest.toml`;
                                        dir = P14_REPO_ROOT)))
        else
            @info "git unavailable — skipping the executable spike-environment freeze assertions"
        end
    end

    @testset "the dependency name set is exactly the frozen sixteen (D-02)" begin
        # Re-asserted FROM THIS FILE rather than by editing `runtests.jl` or the Phase-12 guard:
        # each phase states its own claim in its own file (`test_p12_decoupling.jl:187-188`).
        @test Set(keys(Pkg.project().dependencies)) == Set([
            "BenchmarkTools", "CSV", "CairoMakie", "CoordinateTransformations", "DataFrames",
            "Distributions", "Flux", "HypothesisTests", "ImageFiltering", "ImageTransformations",
            "Images", "Interpolations", "JLD2", "NeuralEstimators", "Random123", "StatsBase"])
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"

        # THE PHASE-14-SPECIFIC HALF, ASSERTED BY NAME -- one @test each, so the failure says
        # WHICH package appeared.
        #
        # `ConformalPrediction` is named in the ORIGINAL SC1, and SC1 IS AMENDED BY D-02
        # (14-CONTEXT.md, "Decision machinery"): the library name is DROPPED and conformal is met
        # IN SUBSTANCE, hand-rolled. Split conformal is a sorted-score order statistic; the
        # library wraps MLJ models rather than a NeuralEstimators posterior, so it is the one
        # option with a concrete verifiable cost and the least payoff. Adding it drags in the MLJ
        # tree and forces a `Pkg.resolve()` against the frozen manifest -- and this milestone has
        # already paid for that once, when the Phase-7 Wave-0 co-resolution gate FAILED on a
        # transitive Makie conflict that pushed NeuralEstimators back below the v0.2.1 pin the
        # whole spike is built on.
        @test !haskey(Pkg.project().dependencies, "ConformalPrediction")
        @test !haskey(Pkg.project().dependencies, "MLJ")
        # `roc_auc` and the risk-coverage curve are computed with the in-repo hand-rolled helper
        # (`spike/validation/ood.jl:319`), never by adding an ROC package.
        @test !haskey(Pkg.project().dependencies, "ROCAnalysis")
        # Bayesian FDR is a sort over posterior probabilities, not a frequentist p-value
        # adjustment; a multiple-testing package would be the wrong tool as well as a new dep.
        @test !haskey(Pkg.project().dependencies, "MultipleTesting")
        @test !haskey(Pkg.project().dependencies, "Distances")
        # CPU-only is a CLAUDE.md constraint, not merely a preference.
        @test !haskey(Pkg.project().dependencies, "CUDA")
    end

    @testset "the sealed holdout is never consumed (D-03a)" begin
        # THE PROPERTY, STATED BEFORE THE CODE IT BINDS EXISTS: no byte of a sealed corpus IMAGE
        # is ever read by Phase 14. It is NOT "the string `corpus` appears nowhere", and the
        # difference is load-bearing in four places this phase actually requires:
        #   * `spike/p14/consts.jl` records the D-03a substitution as a DECLARED DEVIATION, and
        #     that record has to name what was substituted away FROM;
        #   * `spike/test/test_p14_consts.jl` MUST name CORPUS_MASTER_SEED in the seed-
        #     disjointness assertion -- a token ban would forbid the pre-registration's own test
        #     from listing the very seed it exists to forbid;
        #   * THIS file must contain `corpus` and `sealed_holdout` to express the rule at all;
        #   * `run_p14_real_images.jl` (SC1-f) must READ `corpus/manifest.csv` to record
        #     `corpus_images_available == 0` as UNFETCHED BY DESIGN with its tier breakdown. That
        #     is a METADATA read of a committed text file whose whole purpose is provenance; it
        #     opens no image.
        #
        # WHY THE EXCLUSION IS ABSOLUTE. D-03a, verified directly against `corpus/manifest.csv`:
        # 30 of 32 rows are `tier: simulated-secondary` (the CBS benchmark -- another simulator),
        # and the only 2 `physical-primary` rows are `split = sealed_holdout`, `bytes = 0`,
        # `sha256 = PENDING-FETCH`, RESERVED for the Phase-16 blind evaluation. Scoring one here
        # would irreversibly burn Phase 16 -- there is no second blind set. Phase 13 hit this
        # exact problem and resolved it the same way (`spike/p13/real_images.jl:65-73`): the
        # SANCTIONED SUBSTITUTE IS THE SIX COMMITTED MICROSCOPY TIFFS UNDER test/test_images/.
        sources = _p14_sources()
        @test !isempty(sources)                       # the scan has something to scan
        # The allowlist must describe files that really are in the scanned surface, or an
        # exemption could silently protect nothing (a typo, or a renamed file).
        @test all(a -> !isfile(joinpath(P14_REPO_ROOT, a)) || a in sources,
                  _P14_CORPUS_METADATA_ALLOWED)

        # (a) NO IMAGE LOAD RESOLVES ONTO THE CORPUS, in ANY Phase-14 source -- allowlisted or
        #     not. Metadata is exempt; opening a row never is.
        load_hits = String[]
        for p in sources
            for ln in _sealed_image_load_hits(_p14_body(p))
                push!(load_hits, "$p: $ln")
            end
        end
        @test load_hits == String[]

        # (b) NON-allowlisted sources may not mention the corpus at all. Comment-stripped first,
        #     so prose about the holdout is free and a CODE reference is not.
        token_hits = String[]
        for p in sources
            p in _P14_CORPUS_METADATA_ALLOWED && continue
            body = _p14_body(p)
            for tok in _P14_CORPUS_TOKENS
                occursin(tok, body) && push!(token_hits, "$p: $tok")
            end
        end
        @test token_hits == String[]

        # (c) The allowlisted files get the NARROWER property: metadata yes, the sealed IMAGE
        #     DIRECTORY never -- neither spelled as a path nor reached through the constant
        #     `corpus/config.jl` binds to it, and never through the sole reachable opener.
        image_dir_hits = String[]
        for p in _P14_CORPUS_METADATA_ALLOWED
            isfile(joinpath(P14_REPO_ROOT, p)) || continue
            body = _p14_body(p)
            occursin(_SEALED_IMAGE_DIR_LITERAL, body) &&
                push!(image_dir_hits, "$p: $(_SEALED_IMAGE_DIR_LITERAL)")
            occursin(_SEALED_IMAGE_DIR_SYMBOL, body) &&
                push!(image_dir_hits, "$p: $(_SEALED_IMAGE_DIR_SYMBOL)")
            occursin(_SEALED_OPENER, body) &&
                push!(image_dir_hits, "$p: $(_SEALED_OPENER)")
        end
        @test image_dir_hits == String[]

        # (d) The corpus tree is byte-unchanged and untracked-clean, AND the PERMITTED real-image
        #     path is asserted POSITIVELY so the three negatives above cannot be satisfied
        #     vacuously by a phase that reads no images at all. The six committed TIFFs are that
        #     path (D-03a). The `if isfile(...)` half is the arms-itself-later idiom: the check is
        #     derived from the FILE rather than from a fixed list precisely so it arms itself when
        #     the SC1-f runner lands, instead of being faked now or forgotten later.
        for cond in ("positive", "negative"), ch in ("c1", "c2", "c3")
            @test isfile(joinpath(P14_REPO_ROOT, "test", "test_images", cond, "$(cond)_$(ch).tif"))
        end
        _real_runner = joinpath(P14_REPO_ROOT, "spike", "p14", "run_p14_real_images.jl")
        if isfile(_real_runner)
            @test any(p -> occursin("test/test_images", _p14_body(p)), sources)
        else
            @info "spike/p14/run_p14_real_images.jl does not exist yet (SC1-f) — skipping the " *
                  "positive test/test_images assertion; it arms itself when the runner lands"
        end
        if _HAS_GIT
            @test success(Cmd(`git diff --quiet HEAD -- corpus`; dir = P14_REPO_ROOT))
            @test isempty(readchomp(Cmd(`git status --porcelain -- corpus`; dir = P14_REPO_ROOT)))
        else
            @info "git unavailable — skipping the executable corpus/ byte-equality assertion"
        end
    end

    @testset "the shipped bundle and the Phase-13 artifacts are read-only" begin
        # Phase 14 is RESEARCH-LANE (D-01): no reship, no new shipped artifact, no ship gate
        # reopened. The Phase-7 GO rests on the `amended_v2/grid_8` bundle, which is
        # content-hashed and Release-hosted rather than committed, and `artifacts/` is gitignored
        # here -- so the git half is SKIPPED rather than passing vacuously against a path git was
        # never watching.
        if _HAS_GIT
            # Phase 13 is this phase's INPUT: the net, the tau probe report and the gate report
            # are read by `p14_load_tau()` and must not move underneath it.
            @test isempty(readchomp(Cmd(`git status --porcelain -- spike/p13`; dir = P14_REPO_ROOT)))
            _tracked = readchomp(Cmd(`git ls-files -- artifacts/amended_v2/grid_8`;
                                     dir = P14_REPO_ROOT))
            if isempty(_tracked)
                @info "artifacts/amended_v2/grid_8 is not tracked in this checkout — skipping " *
                      "the executable shipped-bundle assertion (the bundle is Release-hosted)"
            else
                @test success(Cmd(`git diff --quiet HEAD -- artifacts`; dir = P14_REPO_ROOT))
                @test isempty(readchomp(Cmd(`git status --porcelain -- artifacts`;
                                            dir = P14_REPO_ROOT)))
            end
        else
            @info "git unavailable — skipping the shipped-bundle and Phase-13 read-only assertions"
        end

        # No Phase-14 source WRITES under `artifacts/`. Reading the shipped bundle as a comparison
        # reference is permitted; adding to it is a ship-gate act.
        artifact_write_hits = String[]
        for p in _p14_sources()
            for ln in split(_p14_body(p), '\n')
                s = strip(ln)
                occursin("artifacts", s) || continue
                any(t -> occursin(t, s), _ARTIFACT_WRITERS) && push!(artifact_write_hits, "$p: $s")
            end
        end
        @test artifact_write_hits == String[]
    end

    @testset "Phase-14 caches use a NEW root, never Phase 13's" begin
        # Phase 13's pool is gitignored bulk data that the Phase-14 OOD null fit READS
        # (`h.meta.pool_dir`), and this milestone has already destroyed a 54 MB gitignored pool
        # once. A Phase-14 datagen run that reached for the same cache root would overwrite it
        # SILENTLY -- a shard writer's job is to write shards.
        p13_root = normpath(joinpath(P14_REPO_ROOT, "spike", "data", "cache", "p13"))
        p14_root = normpath(joinpath(P14_REPO_ROOT, "spike", "data", "cache", "p14"))
        @test p13_root != p14_root
        @test !startswith(p13_root, p14_root * Base.Filesystem.path_separator)

        # READING the Phase-13 pool is permitted and required; WRITING into it is not. Only the
        # conjunction is banned, so `h.meta.pool_dir` stays reachable.
        p13_write_hits = String[]
        for p in _p14_sources()
            for ln in split(_p14_body(p), '\n')
                s = strip(ln)
                occursin(_P13_CACHE_LITERAL, s) || continue
                (any(t -> occursin(t, s), _CACHE_WRITERS) ||
                 (occursin("open(", s) && occursin("\"w\"", s))) &&
                    push!(p13_write_hits, "$p: $s")
            end
        end
        @test p13_write_hits == String[]
    end

    @testset "decoupling checks ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
