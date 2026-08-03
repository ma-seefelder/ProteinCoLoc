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

# --- The per-file source table (the `test_p13_tau.jl:102-104` triple, generalized) ----------------
# That file reads ONE unit under test as `TAU_SRC` / `TAU_LINES` / `TAU_CODE`, because it guards one
# probe. Phase 14's integrity greps are LANE-WIDE, so the same triple is built once PER FILE and the
# line index is carried with it -- a failure then reports `path:line: the offending code` rather than
# a boolean, which is the difference between a report you can act on and one you have to reproduce.
#
# The `lines` here are the SAME comment-stripped lines `_p14_body` joins, so an index into `lines` is
# meaningful for the positional SC2-c exemption below. It is NOT a line number in the file on disk.
const _P14_TABLE = [begin
                        raw = read(joinpath(P14_REPO_ROOT, p), String)
                        lns = String[l for l in split(raw, '\n') if !startswith(strip(l), "#")]
                        (path = p, lines = lns, code = join(lns, '\n'))
                    end for p in _p14_sources()]

_p14_entry(path::AbstractString) =
    (i = findfirst(e -> e.path == path, _P14_TABLE); i === nothing ? nothing : _P14_TABLE[i])

"Offending `path:index: line` strings for every comment-stripped line satisfying `pred`."
function _p14_line_hits(pred)
    hits = String[]
    for e in _P14_TABLE, (i, ln) in enumerate(e.lines)
        s = strip(ln)
        isempty(s) && continue
        pred(s) && push!(hits, "$(e.path):$i: $s")
    end
    return hits
end

# --- SC2-c: the D-06 anti-regression needles ------------------------------------------------------
# D-06: the OOD input is THREE-VALUED -- fired / clear / not-checked. `src/amortized/local_map.jl`
# states outright that "a `false` flag on a non-sentinel tile means NOT CHECKED, not in
# distribution", and `_recorded_ood_threshold` (:93) returns `nothing` when the gate report is
# absent, unreadable, or records a non-finite threshold. Reading a `false` flag as "in distribution"
# is therefore a SILENT scientific failure: the tool speaks confidently in exactly the regime where
# it has no evidence it should. Every spelling of that read is banned.
#
# Each needle is an ESCAPED regex, and that is also what keeps this file out of its own scan: in
# `r"\.flag\s*==\s*false"` the characters following `.flag` are `\`, `s`, `*` -- not whitespace and
# not `=` -- so the pattern cannot match its own declaration.
const _SC2C_PATTERNS = (r"\.flag\s*==\s*false",
                        r"\.flag\s*===\s*false",
                        r"!\s*\w*\.flag",
                        r"verdict\.flag\s*\?")

"""
The comment-stripped line range of `function p14_ood_state ... end` in a Phase-14 source, or
`nothing` if this file does not define it. The closing `end` is located POSITIONALLY: the first
later line whose stripped text is exactly `end` at the SAME indentation as the `function` line.

The LONG form is required, and deliberately so: `p14_ood_state` is the ONE place allowed to inspect
a raw flag, because deriving the three-valued state is its whole job. A short-form definition would
be un-delimitable, so the exemption would have to be either a whole-file pass (which is not an
exemption, it is a hole) or nothing at all. A missing `end` raises rather than returning `nothing`,
so a malformed exemption fails loudly instead of silently protecting nothing.
"""
function _p14_ood_state_range(entry)
    i0 = findfirst(l -> occursin("function p14_ood_state", l), entry.lines)
    i0 === nothing && return nothing
    indent = length(entry.lines[i0]) - length(lstrip(entry.lines[i0]))
    for j in (i0 + 1):length(entry.lines)
        l = entry.lines[j]
        strip(l) == "end" || continue
        (length(l) - length(lstrip(l))) == indent || continue
        return i0:j
    end
    error("test_p14_decoupling.jl: `function p14_ood_state` was found in $(entry.path) but no " *
          "matching `end` at the same indentation. The SC2-c exemption cannot be located, and an " *
          "exemption that cannot be located must fail loudly rather than protect nothing.")
end

# --- SC2-d, the withdrawn figure, the forbidden seed, the helper bans ----------------------------
# SC2-d. tau = 0.15 is defined on rho_true -- measured by a probe that stops at
# `encode_d01(patch_summary(...))` and never standardizes -- while `patch_correlation` returns the
# INDUCED MEAN PER-PATCH CORRELATION mu. Comparing mu to tau directly is a criterion applied in a
# unit it was not derived for, which is one of the four recorded "the bar was wrong" families
# (14-RESEARCH Pitfall 5; the same failure `P12_STAGE1_RATIO_CEILING` carries a recorded warning
# about). `ghat` is the frozen, monotone, clamped map from mu to rho_true and is the ONLY in-repo
# bridge between the two scales.
const _PC_CALL       = "patch_" * "correlation("
const _TAU_MENTIONS  = ("tau", "τ")
const _GHAT_CALL     = "ghat("

# THE WITHDRAWN FIGURE. D-03a withdraws "alpha >= 1/31 ~ 0.032": it was computed from the wrong n
# (the corpus is partitioned dev = 14 / eval = 16 / sealed = 2, not "30 open") AND on rows that are
# not physical truth (30 of 32 are `simulated-secondary`). This check exists so a corpus-derived
# floor can never become an EXECUTABLE DEFAULT -- which would be the corpus silently calibrating the
# hedge that D-03a just removed it from.
#
# THE EXEMPTION IS DELIBERATE AND NARROW: `p14_conformal_quantile`'s GENERAL feasibility expression
# `1/(n + 1)`, computed from its own `n` argument, is permitted and expected. Only a HARD-CODED 31
# or the decimal 0.032 is flagged. Needles assembled from fragments so this file passes its own scan.
const _WITHDRAWN_FIGURE_NEEDLES = ("0." * "032", "1/" * "31", "1 / " * "31", "1/(30" * "+1)")

# PITFALL 7. `spike/comparator/config.jl:56` sets `MASTER_SEED = 0x00000000_00C0FFEE`, which IS
# `NPE_MASTER_SEED` in the forbidden-seed inventory. `costes_p` takes `master_seed` as its first
# positional argument (`classical.jl:181`), so an executor copying the Phase-9 call site verbatim
# would consume a forbidden stream. `test_p14_consts.jl` CANNOT catch this: it checks CONSTANTS, not
# CALL SITES, and the constant it asserts disjoint is exactly the one that would be passed here.
#
# The banned call site, spelled out once so a reader can see exactly what is looked for -- and
# spelled out HERE, in a comment, because the scan strips comments and would otherwise flag this
# file for naming the thing it forbids:
#
#     costes_p(MASTER_SEED, idx, mci)     # FORBIDDEN -- pass P14_DEV_SEED as the first argument
#
# and likewise any `costes_p(COSTES...` spelling that reaches for the comparator's own salt.
const _FORBIDDEN_COSTES_CALLS = ("costes_p(" * "MASTER_SEED", "costes_p(" * "COSTES")

# No dependency may be added (D-02), so every shared helper is REUSED. A local re-implementation is
# how two versions of one statistic silently diverge -- and `ghat` in particular is FROZEN
# calibration, not a formula to be retyped.
const _NO_REIMPLEMENT = ("roc_auc", "corspearman", "_bin_calibration", "three_way_label", "ghat")

# D-07: alpha_FDR is a per-call USER parameter, alpha_conformal is a pre-registered miscoverage
# level. They happen to share the value 0.10 at one grid point, which is exactly why neither may
# ever be ASSIGNED from the other: the day one moves, a derived assignment moves the other silently
# and the two stop being distinct quantities at all.
const _ALPHA_CROSS_ASSIGN = (r"alpha_fdr\s*=\s*.*ALPHA_CONFORMAL"i,
                             r"alpha_conformal\s*=\s*.*alpha_fdr"i)

# The aggregate suite exits 1 early at the Phase-4 SPEEDUP_GATE, masking every later include block,
# so a Phase-14 source that reached for runtests.jl would be wiring its evidence to a runner that
# never reaches it. Per-file runs are the only reliable signal (14-PATTERNS §4.8).
const _RUNTESTS_NEEDLE = "runtests" * ".jl"

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

    @testset "SC2-c: no `false` OOD flag is read as in-distribution (D-06)" begin
        # See `_SC2C_PATTERNS` above for why this is a scientific failure rather than a style one.
        # The scan is lane-wide and the exemption is POSITIONAL: only the body of `p14_ood_state`
        # may inspect a raw flag, because turning a raw flag into the three-valued state is that
        # function's entire job.
        sc2c_hits = String[]
        for e in _P14_TABLE
            exempt = 0:-1
            if e.path == "spike/p14/fuse.jl"
                r = _p14_ood_state_range(e)
                r === nothing || (exempt = r)
            end
            for (i, ln) in enumerate(e.lines)
                i in exempt && continue
                s = strip(ln)
                isempty(s) && continue
                any(p -> occursin(p, s), _SC2C_PATTERNS) && push!(sc2c_hits, "$(e.path):$i: $s")
            end
        end
        @test sc2c_hits == String[]

        # THE POSITIVE HALF, so the guard cannot be satisfied VACUOUSLY by a file that simply never
        # mentions the three-valued state. `fuse.jl` lands in wave 2; this arms itself then.
        _fuse = _p14_entry("spike/p14/fuse.jl")
        if _fuse === nothing
            @info "spike/p14/fuse.jl does not exist yet (wave 2) — the positive `:not_checked` " *
                  "assertion arms itself when it lands"
        else
            @test occursin(":not_checked", _fuse.code)
            # The long form is REQUIRED so the exemption above is delimitable at all.
            @test occursin("function p14_ood_state", _fuse.code)
            @test _p14_ood_state_range(_fuse) !== nothing
        end
    end

    @testset "SC2-d: every classical statistic compared to tau passes through ghat" begin
        # PER-LINE, not whole-file: the question is whether THIS comparison was rescaled, and a
        # `ghat` call somewhere else in the file answers a different question.
        unit_hits = _p14_line_hits(s -> occursin(_PC_CALL, s) &&
                                        any(t -> occursin(t, s), _TAU_MENTIONS) &&
                                        !occursin(_GHAT_CALL, s))
        @test unit_hits == String[]

        # THE POSITIVE HALF, scoped to `decide.jl` rather than to the lane: a lane-wide `occursin`
        # would be satisfied by THIS file's own needles, which is not evidence about anything.
        _decide = _p14_entry("spike/p14/decide.jl")
        if _decide === nothing
            @info "spike/p14/decide.jl does not exist yet (wave 2) — the positive " *
                  "ghat-wrapping assertion arms itself when it lands"
        else
            @test occursin("ghat(patch_correlation(", _decide.code)
            # And the basis is RECORDED on the result, so the rescaling is auditable rather than
            # implicit (14-RESEARCH §D.5).
            @test occursin(":ghat_rho_true_scale", _decide.code)
        end
    end

    @testset "the WITHDRAWN corpus-derived alpha floor is not an executable default (D-03a)" begin
        # 0.032 and 1/31 are withdrawn figures, not merely superseded ones -- see the note at
        # `_WITHDRAWN_FIGURE_NEEDLES`. The general `1/(n + 1)` feasibility expression is exempt.
        figure_hits = _p14_line_hits(s -> any(n -> occursin(n, s), _WITHDRAWN_FIGURE_NEEDLES))
        @test figure_hits == String[]
    end

    @testset "costes_p never consumes the comparator MASTER_SEED (Pitfall 7)" begin
        seed_hits = _p14_line_hits(s -> any(n -> occursin(n, s), _FORBIDDEN_COSTES_CALLS))
        @test seed_hits == String[]

        # THE POSITIVE HALF: once the composition point exists, EVERY `costes_p(` call site must
        # name the Phase-14 dev stream on the same line. Scoped to `decide.jl`, and collecting the
        # offending call sites rather than a boolean.
        _decide = _p14_entry("spike/p14/decide.jl")
        if _decide === nothing
            @info "spike/p14/decide.jl does not exist yet (wave 2) — the positive " *
                  "P14_DEV_SEED-at-every-costes-call-site assertion arms itself when it lands"
        else
            bad = String[]
            for (i, ln) in enumerate(_decide.lines)
                s = strip(ln)
                occursin("costes_p(", s) || continue
                occursin("P14_DEV_SEED", s) || push!(bad, "spike/p14/decide.jl:$i: $s")
            end
            @test bad == String[]
        end
    end

    @testset "shared helpers are REUSED, never re-implemented (D-02)" begin
        reimpl_hits = String[]
        for name in _NO_REIMPLEMENT
            needle = "function " * name
            append!(reimpl_hits, _p14_line_hits(s -> occursin(needle, s)))
        end
        @test reimpl_hits == String[]

        # The positive twin (`test_p13_tau.jl:275-315`): a ban alone is satisfied by never using
        # the helper at all, which is not what "reuse" means.
        _rc = _p14_entry("spike/p14/run_p14_riskcoverage.jl")
        if _rc === nothing
            @info "spike/p14/run_p14_riskcoverage.jl does not exist yet (SC3-a/b/c) — the " *
                  "positive roc_auc / corspearman reuse assertions arm themselves when it lands"
        else
            @test occursin("roc_auc(", _rc.code)
            @test occursin("corspearman", _rc.code)
        end
    end

    @testset "the two alphas are never derived from one another (D-07)" begin
        alpha_hits = _p14_line_hits(s -> any(p -> occursin(p, s), _ALPHA_CROSS_ASSIGN))
        @test alpha_hits == String[]
    end

    @testset "the aggregate suite is NOT a Phase-14 gate" begin
        # THIS TESTSET'S NAME IS DELIBERATELY PARAPHRASED, and the first run is why: a `@testset`
        # title is a STRING IN CODE, not a comment, so naming the runner in the title made this
        # very check report itself as the lane's only offender. Recorded rather than quietly
        # renamed -- it is the cheapest available demonstration that the scan is live and that
        # comment-stripping is not the same thing as string-stripping.
        runner_hits = _p14_line_hits(s -> occursin(_RUNTESTS_NEEDLE, s))
        @test runner_hits == String[]
    end

    @testset "decoupling checks ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
