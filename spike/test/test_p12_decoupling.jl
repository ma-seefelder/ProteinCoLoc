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

# spike/test/test_p12_decoupling.jl --- the phase's three hard constraints, made EXECUTABLE.
#
# WHY THIS FILE EXISTS. CLAUDE.md's decoupling rule ("the main package, `src/`, and both
# manuscript pipelines must remain provably untouched during the spike"; "spike/Project.toml /
# spike/Manifest.toml byte-unchanged") and D-11's sealed-holdout exclusion are both stated in
# PROSE, in several documents. Prose does not fail a build. Each rule below is an assertion that
# turns the suite red the moment it stops being true -- including the one nobody would notice, a
# Phase-12 code path quietly reading a `corpus/` row.
#
# IT RUNS EARLY IN THE PHASE ON PURPOSE. An assertion added at phase close proves only that
# nothing broke SINCE the assertion was added. These assert the state that holds now, in the wave
# where it still holds, so every later wave is measured against it.
#
# THE SOURCE SCAN SCANS *ITSELF*, AND THAT IS WHY THREE NEEDLES ARE BUILT BY CONCATENATION.
# This file matches `spike/test/test_p12_*.jl`, so it is one of the Phase-12 sources testset 4
# enumerates. A naive token ban is therefore RED BY CONSTRUCTION: the file has to name the thing
# it forbids in order to forbid it. Two answers are used, and they are different answers to
# different halves of the problem:
#   * the ALLOWLIST (`_P12_CORPUS_METADATA_ALLOWED`) exempts the four files that legitimately
#     mention corpus METADATA -- see testset 4's comment for why each of the four must;
#   * the three needles that would still match this file even under the allowlist
#     (`open_` + `sealed_holdout`, and the two spellings of the sealed IMAGE directory) are
#     assembled at run time from fragments, so the contiguous literal never appears in this
#     source at all. That is what lets the scan include this file rather than exempt it.
# Neither trick weakens anything: the needle a run searches for is the full token either way.
#
# DECOUPLING (CLAUDE.md): spike-local. Reads the Tier-1 pre-registration, shells out to `git`
# read-only, and reads Phase-12 source text. Writes nothing, adds no package, touches no `src/`
# file and opens no image.

using Test
using Pkg

# The repo root is READ from the Tier-1 pre-registration, never re-derived here: p12_consts.jl
# §12 is explicit that "12-05 READS these bindings rather than re-deriving its own", because two
# derivations can drift and only one of them would be the pre-registered one. Guarded for
# idempotency (S2), the same idiom `test_p12_consts.jl:41-42` uses.
isdefined(@__MODULE__, :P12_DEV_SEED) ||
    include(joinpath(@__DIR__, "..", "validation", "p12_consts.jl"))

# The ONE thing Tier 1 does not answer, so it stays local: is `git` usable in this sandbox at
# all? Exactly the predicate `spike/test/test_p13_result.jl:199-201` uses -- a source tree that
# is not a checkout skips the executable half rather than failing for an unrelated reason.
# `isdir(...) || isfile(...)` because inside a linked worktree `.git` is a FILE.
const _HAS_GIT = Sys.which("git") !== nothing &&
                 (isdir(joinpath(P12_REPO_ROOT, ".git")) ||
                  isfile(joinpath(P12_REPO_ROOT, ".git")))

# Source text with whole-line comments stripped, so a path MENTIONED IN PROSE cannot be mistaken
# for a call. Ported from `spike/test/test_stage6_regression.jl:79-80`.
_strip_comment_lines(src::AbstractString) =
    join(filter(l -> !startswith(strip(l), "#"), split(src, '\n')), '\n')

# --- The Phase-12 source surface ---------------------------------------------------------------
# Every file this phase may add, as (directory, filename pattern). Files that do not exist yet
# are simply absent from the enumeration: this plan runs at wave 2 and most of these land later,
# so the scan grows with the phase instead of asserting against a frozen list that goes stale.
const _P12_SOURCE_PATTERNS = (
    ("spike/validation", r"^p12_.*\.jl$"),
    ("spike/validation", r"^run_p12_.*\.jl$"),
    ("spike/simulator",  r"^p12_.*\.jl$"),
    ("spike/npe",        r"^p12_architecture\.jl$"),
    ("spike/npe",        r"^train_p12_npe\.jl$"),
    ("spike/data",       r"^p12_generate\.jl$"),
    ("spike/p12",        r"^.*\.jl$"),
    ("spike/test",       r"^test_p12_.*\.jl$"),
    ("spike/test",       r"^capture_p12_golden\.jl$"),
)

"Repo-relative, forward-slashed paths of every Phase-12 source that EXISTS right now."
function _p12_sources()
    out = String[]
    for (rel, pat) in _P12_SOURCE_PATTERNS
        d = joinpath(P12_REPO_ROOT, rel)
        isdir(d) || continue
        for f in readdir(d)
            occursin(pat, f) || continue
            p = rel * "/" * f
            p in out || push!(out, p)
        end
    end
    return sort(out)
end

_p12_body(relpath::AbstractString) =
    _strip_comment_lines(read(joinpath(P12_REPO_ROOT, relpath), String))

# --- The corpus needles -------------------------------------------------------------------------
# THE EXEMPTION SURFACE, DECLARED AS DATA SO A READER SEES ALL OF IT AT ONCE. It is deliberately
# short. Adding a fifth entry is a DECISION, not a fix -- if a new Phase-12 file needs to say
# `corpus`, the question to answer first is why.
const _P12_CORPUS_METADATA_ALLOWED = ("spike/validation/p12_consts.jl",
                                      "spike/test/test_p12_consts.jl",
                                      "spike/test/test_p12_decoupling.jl",
                                      "spike/validation/run_p12_coverage.jl")

# Banned outright in every NON-allowlisted Phase-12 source.
const _P12_CORPUS_TOKENS = ("corpus", "sealed_holdout", "manifest.csv", "CORPUS_MASTER_SEED")
const _CORPUS_PATH_TOKEN = "corpus"
# Image-loading calls. `build_mci` is the spike's own MultiChannelImage constructor
# (`spike/simulator/*`), so it is the third way a byte of image data enters this phase.
const _IMAGE_LOADERS = ("load(", "FileIO.load", "build_mci")

# ASSEMBLED AT RUN TIME (see the banner): the contiguous literals must not appear in this source,
# because this source is itself scanned. The value searched for is the full token.
const _SEALED_OPENER            = "open_" * "sealed_holdout"
const _SEALED_IMAGE_DIR_LITERAL = "cor" * "pus/data"
const _SEALED_IMAGE_DIR_SYMBOL  = "CORPUS_" * "DATA_DIR"

"""
Lines of a comment-stripped body that would read a sealed image: the sole reachable opener
(`corpus/manifest.jl` D-09), or an image-loading call on a path expression naming the corpus.
Returns the offending lines so a failure shows the code, not just a boolean.
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

@testset "P12 decoupling: src/, deps and the sealed holdout (12-05)" verbose = true begin

    @testset "src/ is byte-unchanged (CLAUDE.md hard constraint)" begin
        if _HAS_GIT
            # Whole-tree first, then the four files 12-PATTERNS §0.2 names as READ-ONLY
            # references, so a failure NAMES the file instead of saying "something under src/".
            @test success(Cmd(`git diff --quiet HEAD -- src`; dir = P12_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/results.jl`; dir = P12_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/amortized/summary.jl`; dir = P12_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/amortized/architecture.jl`; dir = P12_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/amortized/local_map.jl`; dir = P12_REPO_ROOT))
            # THE UNTRACKED CASE, WHICH `git diff HEAD` CANNOT SEE. A NEW file dropped under
            # `src/` is invisible to a diff against HEAD and would pass every assertion above
            # while being exactly the productionization-by-stealth this constraint forbids.
            @test isempty(readchomp(Cmd(`git status --porcelain -- src`; dir = P12_REPO_ROOT)))
        else
            @info "git unavailable — skipping the executable src/ byte-equality assertions"
        end
    end

    @testset "the spike environment is byte-frozen" begin
        # 12-PATTERNS §0.1: the spatial map's stdlib reach (LinearAlgebra, SparseArrays) rides
        # the default LOAD_PATH `@stdlib` entry and needs no Project.toml entry, so ANY movement
        # in these two files means a `Pkg.add`/`Pkg.resolve` happened -- which is the thing the
        # frozen manifest exists to prevent.
        if _HAS_GIT
            @test success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                              dir = P12_REPO_ROOT))
            @test isempty(readchomp(Cmd(`git status --porcelain -- spike/Project.toml spike/Manifest.toml`;
                                        dir = P12_REPO_ROOT)))
        else
            @info "git unavailable — skipping the executable spike-environment freeze assertions"
        end
    end

    @testset "the dependency name set is exactly the frozen sixteen (resolve-risk clause k)" begin
        # (k) RESOLVE-RISK GATE (Phase 12). Re-asserted FROM THIS FILE rather than by editing
        # `runtests.jl:114-127`: the Phase-11 block there is a committed record of what Phase 11
        # froze, and each phase states its own claim in its own file.
        @test Set(keys(Pkg.project().dependencies)) == Set([
            "BenchmarkTools", "CSV", "CairoMakie", "CoordinateTransformations", "DataFrames",
            "Distributions", "Flux", "HypothesisTests", "ImageFiltering", "ImageTransformations",
            "Images", "Interpolations", "JLD2", "NeuralEstimators", "Random123", "StatsBase"])
        @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"

        # THE PHASE-12-SPECIFIC HALF. The spatial map adds NO new package: at G = 8 the lattice
        # is 64 cells, so the CAR precision and the GP covariance are DENSE 64x64 matrices under
        # stdlib `LinearAlgebra` -- there is nothing for a Gaussian-process package to do that a
        # 64x64 Cholesky does not already do in three lines.
        #
        # WHY THE BAN IS NOT MERELY TIDINESS. Adding any of the three below forces a
        # `Pkg.resolve()` against the frozen manifest, and this milestone has already paid for
        # that once: the Phase-7 Wave-0 co-resolution gate FAILED on a transitive Makie conflict
        # that pushed NeuralEstimators back toward 0.1.4, below the v0.2.1 pin the whole spike is
        # built on. Asserted BY NAME, one test each, so the failure says WHICH one appeared.
        @test !haskey(Pkg.project().dependencies, "AbstractGPs")
        @test !haskey(Pkg.project().dependencies, "KernelFunctions")
        @test !haskey(Pkg.project().dependencies, "GaussianRandomFields")
    end

    @testset "the Phase-16 sealed holdout is not consumed (D-11 AMENDED)" begin
        # THE PROPERTY, STATED BEFORE IT IS IMPLEMENTED: no byte of a sealed corpus IMAGE is ever
        # read by Phase 12. It is NOT "the string `corpus` appears nowhere", and the difference is
        # load-bearing in three places this phase actually requires:
        #   * `p12_consts.jl` MUST define CORPUS_MASTER_SEED -- it belongs in the forbidden-seed
        #     inventory (mirroring `p11_consts.jl:90`), and `test_p12_consts.jl` MUST name it in
        #     the disjointness assertion. A token ban would forbid the pre-registration file from
        #     listing the very seed it exists to forbid.
        #   * THIS file must contain `corpus` and `sealed_holdout` to express the rule at all.
        #   * 12-19 must READ `corpus/manifest.csv` to re-derive `corpus_tier_counts` and
        #     `corpus_images_available = 0`. That is a METADATA read of a committed text file whose
        #     whole purpose is provenance; it opens no image, and it is what makes the
        #     reported-not-gated demotion of the real-image arm checkable.
        #
        # WHY THE EXCLUSION IS ABSOLUTE. Both `physical-primary` rows of `corpus/manifest.csv`
        # (`pos-tetraspeck-01`, `neg-lightmycells-01`) carry `split = sealed_holdout`,
        # `sha256 = PENDING-FETCH` and `bytes = 0`: nothing is fetched, and they are RESERVED for
        # the Phase-16 blind evaluation. Scoring them here would irreversibly burn Phase 16 --
        # there is no second blind set.
        sources = _p12_sources()
        @test !isempty(sources)                       # the scan has something to scan
        # The allowlist must describe files that really are in the scanned surface, or an
        # exemption could silently protect nothing (a typo, or a renamed file).
        @test all(a -> !isfile(joinpath(P12_REPO_ROOT, a)) || a in sources,
                  _P12_CORPUS_METADATA_ALLOWED)

        # (a) NO IMAGE LOAD RESOLVES ONTO THE CORPUS, in ANY Phase-12 source -- allowlisted or
        #     not. Metadata is exempt; opening a row never is.
        load_hits = String[]
        for p in sources
            for ln in _sealed_image_load_hits(_p12_body(p))
                push!(load_hits, "$p: $ln")
            end
        end
        @test load_hits == String[]

        # (b) NON-allowlisted sources may not mention the corpus at all. Comment-stripped first,
        #     so prose about the holdout is free and a CODE reference is not.
        token_hits = String[]
        for p in sources
            p in _P12_CORPUS_METADATA_ALLOWED && continue
            body = _p12_body(p)
            for tok in _P12_CORPUS_TOKENS
                occursin(tok, body) && push!(token_hits, "$p: $tok")
            end
        end
        @test token_hits == String[]

        #     The allowlisted four get the NARROWER property: metadata yes, the sealed IMAGE
        #     DIRECTORY never -- neither spelled as a path nor reached through the constant
        #     `corpus/config.jl` binds to it.
        image_dir_hits = String[]
        for p in _P12_CORPUS_METADATA_ALLOWED
            isfile(joinpath(P12_REPO_ROOT, p)) || continue
            body = _p12_body(p)
            occursin(_SEALED_IMAGE_DIR_LITERAL, body) &&
                push!(image_dir_hits, "$p: $(_SEALED_IMAGE_DIR_LITERAL)")
            occursin(_SEALED_IMAGE_DIR_SYMBOL, body) &&
                push!(image_dir_hits, "$p: $(_SEALED_IMAGE_DIR_SYMBOL)")
        end
        @test image_dir_hits == String[]

        # (c) The corpus tree is byte-unchanged, and the PERMITTED real-image path is asserted
        #     POSITIVELY so the three negatives above cannot be satisfied by a phase that reads no
        #     images at all. `test/test_images/` (tracked, `positive/` + `negative/`) is that path.
        #     No coverage runner exists yet, so this half is skipped rather than faked -- 12-19 is
        #     the plan that lands it.
        _coverage_runner = joinpath(P12_REPO_ROOT, "spike", "validation", "run_p12_coverage.jl")
        if isfile(_coverage_runner)
            @test any(p -> occursin("test/test_images", _p12_body(p)), sources)
        else
            @info "run_p12_coverage.jl does not exist yet (12-19) — skipping the positive " *
                  "test/test_images assertion"
        end
        if _HAS_GIT
            @test success(Cmd(`git diff --quiet HEAD -- corpus`; dir = P12_REPO_ROOT))
            @test isempty(readchomp(Cmd(`git status --porcelain -- corpus`; dir = P12_REPO_ROOT)))
        else
            @info "git unavailable — skipping the executable corpus/ byte-equality assertion"
        end
    end

    @testset "the shipped bundle is untouched" begin
        # The Phase-7 GO rests on the `amended_v2/grid_8` bundle: content-hashed
        # (`git-tree-sha1 = 90e6b63a...`) and Release-hosted, NOT committed into this tree.
        # Phase 12 is RESEARCH-LANE: no reship, no new shipped artifact, no ship gate reopened.
        # `artifacts/` is gitignored here (`.gitignore:63`), so the executable half is skipped
        # rather than passing vacuously against a path git was never watching.
        if _HAS_GIT
            _tracked = readchomp(Cmd(`git ls-files -- artifacts/amended_v2/grid_8`;
                                     dir = P12_REPO_ROOT))
            if isempty(_tracked)
                @info "artifacts/amended_v2/grid_8 is not tracked in this checkout — skipping " *
                      "the executable shipped-bundle assertion (the bundle is Release-hosted)"
            else
                @test success(Cmd(`git diff --quiet HEAD -- artifacts`; dir = P12_REPO_ROOT))
                @test isempty(readchomp(Cmd(`git status --porcelain -- artifacts`;
                                            dir = P12_REPO_ROOT)))
            end
        else
            @info "git unavailable — skipping the executable shipped-bundle assertion"
        end
    end

    @testset "the Phase-11 pool is read-only to Phase 12" begin
        # Phase 11's 50 000-pair pool is 54 MB of gitignored bulk data that 12-15 and 12-17 both
        # depend on, and it has already been destroyed once in this milestone (by a worktree
        # cleanup). A Phase-12 datagen run that reached for the same cache root would overwrite it
        # SILENTLY -- the shard writer's job is to write shards.
        p11_root = normpath(joinpath(P12_REPO_ROOT, "spike", "data", "cache", "p11"))
        p12_root = normpath(joinpath(P12_REPO_ROOT, "spike", "data", "cache", "p12"))
        @test p11_root != p12_root
        @test !startswith(p11_root, p12_root * Base.Filesystem.path_separator)

        # The two literals above are only as good as their agreement with the DECLARING sites, so
        # each is checked against its own file where that file exists. `p11_generate.jl` exists
        # today (this half is live); `p12_generate.jl` lands at 12-09 (this half arms itself then,
        # which is the point of deriving the check from the file rather than from a fixed list).
        _p11_gen = joinpath(P12_REPO_ROOT, "spike", "data", "p11_generate.jl")
        if isfile(_p11_gen)
            @test occursin("joinpath(@__DIR__, \"cache\", \"p11\")",
                           _strip_comment_lines(read(_p11_gen, String)))
        end
        _p12_gen = joinpath(P12_REPO_ROOT, "spike", "data", "p12_generate.jl")
        if isfile(_p12_gen)
            _p12_gen_body = _strip_comment_lines(read(_p12_gen, String))
            @test occursin("joinpath(@__DIR__, \"cache\", \"p12\")", _p12_gen_body)
            @test !occursin("\"p11\"", _p12_gen_body)
        else
            @info "spike/data/p12_generate.jl does not exist yet (12-09) — its cache-root " *
                  "declaration is checked once it lands"
        end
    end

    @testset "decoupling checks ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
