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

# spike/test/test_p13_real.jl --- the amended D-15 real-image arm (D-15 AMENDED, D-16).
#
# WHAT THIS FILE ASSERTS, AND -- MORE IMPORTANTLY -- WHAT IT DOES NOT.
# The amended D-15 real arm is QUALITATIVE AND IS NOT A GATE. The six committed microscopy TIFFs
# under test/test_images/ carry NO colocalization ground-truth label, so no real-image quantity
# has a pass/fail threshold anywhere in this phase. This file therefore asserts the arm's
# MECHANICS -- that the ingestion is the frozen read chain, that the shared alpha transform's
# invariants hold on real input, that the seal stays sealed, that the fixtures stay unmodified,
# and that the not-a-gate posture is machine-checked -- and it NEVER asserts the CORRECTNESS of
# any verdict, because there is no label against which correctness could be defined. Labelled
# real segregation validation is deferred to the Phase-16 blind evaluation.
#
# THE NAMING CORRECTION APPLIES HERE TOO. positive/ and negative/ are the original package's
# BIOLOGICAL test conditions, not colocalization labels; the "negative" pair measures mean patch
# correlation +0.3815 on the OPERATIVE c2/c3 pair (+0.2481 on the SUPERSEDED c1/c2 pair;
# 13-D15-AMENDMENT.md CHANGE A) and is a POSITIVELY correlated pair either way -- MORE so on the
# corrected pair, not less. Nothing below treats negative/ as an exclusion example.
#
# PITFALL 2b, THE ONE THAT WOULD SILENTLY BREAK THIS SUITE. The naive repair of the _exclude_zero
# trap is to assert that every output pixel is strictly positive. That property is FALSE on this
# substrate BEFORE any alpha is applied -- positive_c2 carries 2 exact-zero pixels and
# negative_c2 carries 3, of 1 414 528, in the source files -- so at alpha = 0, where the output
# IS the input, the naive form would report a CORRECT algorithm as broken. The pre-registered
# rule is P13_ALPHA_ZERO_INVARIANT = :count_iszero_preserved and the assertion below compares
# count(iszero, ...) against the count on the UNMODIFIED channel. The forbidden form is not even
# quoted anywhere in this file, so a source-grep gate cannot be defeated by prose.
#
# THIS SUITE CONSUMES NO RNG STREAM AT ALL, AND THAT IS ITSELF ASSERTED BELOW. A reader will
# assume randomness is involved because every other Phase-13 suite rides a fixture stream. It is
# not: the substrate is six committed files, and the alpha transform is a deterministic function
# of its input, its mask and its floor. There is no seed to burn, no fixture stream to protect,
# and therefore no way for this suite to pre-observe a REPORTED Phase-13 stream (D-01).
#
# Mirrors test_bf.jl / test_p13_alpha.jl: license header -> using Test -> guarded include of the
# unit under test -> fixtures computed ONCE at module level -> one outer @testset. No training,
# no figure call (figures are script artifacts and belong to plan 13-16). CPU-only (D-10),
# spike-local, touches no src/ byte, writes nothing under test/, and never opens the seal.

using Test
using Statistics
using SHA

# Unit under test: the real-image ingestion (pulls p13/consts.jl -> contract.jl ->
# data/encode.jl -> p13/alpha_series.jl).
isdefined(@__MODULE__, :real_tif) ||
    include(joinpath(@__DIR__, "..", "p13", "real_images.jl"))

# --- Fixtures computed ONCE (the test_bf.jl:49-65 idiom) --------------------------------------
# Each real summary is one patch_summary over 1.4 M pixels, so every ladder is built once here
# and the testsets assert against the stored results rather than recomputing.
const REAL_PAIRS  = Dict(c => real_pair(c) for c in P13_REAL_CONDITIONS)
const REAL_MASKS  = Dict(c => ch1_object_mask(REAL_PAIRS[c]) for c in P13_REAL_CONDITIONS)
const REAL_FLOORS = Dict(c => alpha_background_floor(REAL_PAIRS[c][2]) for c in P13_REAL_CONDITIONS)
const REAL_ZEROS  = Dict(c => count(iszero, REAL_PAIRS[c][2]) for c in P13_REAL_CONDITIONS)
const REAL_PS     = Dict(c => patch_summary(build_mci(REAL_PAIRS[c])) for c in P13_REAL_CONDITIONS)
const REAL_MBARS  = Dict(c => real_mbar_of(REAL_PAIRS[c]) for c in P13_REAL_CONDITIONS)
const REAL_LADDER = Dict(c => alpha_ladder(REAL_PAIRS[c]; mask = REAL_MASKS[c],
                                           floor = REAL_FLOORS[c])
                         for c in P13_REAL_CONDITIONS)
const REAL_SUMS   = Dict(c => alpha_ladder_summaries(REAL_PAIRS[c]) for c in P13_REAL_CONDITIONS)
const REAL_PROV   = real_provenance(P13_REAL_SAMPLE)

# Sources with WHOLE-LINE COMMENTS STRIPPED, so a file's own prose -- which must NAME the sealed
# holdout in order to state that it is never opened -- cannot invalidate its own gate (the
# test_p13_alpha.jl:123-124 / test_p13_consts.jl:49-50 idiom).
const P13_DIR    = joinpath(@__DIR__, "..", "p13")
const P13_JL     = sort(filter(f -> endswith(f, ".jl"), readdir(P13_DIR)))
_strip_comments(s) = join(filter(l -> !startswith(strip(l), "#"), split(s, '\n')), '\n')
const P13_CODE   = Dict(f => _strip_comments(read(joinpath(P13_DIR, f), String)) for f in P13_JL)
const REAL_CODE  = P13_CODE["real_images.jl"]
const CONSTS_CODE = P13_CODE["consts.jl"]
# FUTURE-PROOFING THE NOT-A-GATE SOURCE GATE, and a NO-OP today. Under the M1 mechanism
# (13-D15-AMENDMENT.md section 7.1) the operative constants stay in consts.jl, so this `get`
# returns "" and changes nothing. If a LATER amendment ever adds a sibling constants file, the
# bar-name set equality below would otherwise stop seeing it. The `readdir` discovery at :86-88
# is NOT touched: a discovered gate must never be swapped for an enumerated one.
const AMENDED_CODE = get(P13_CODE, "consts_d15_amendment.jl", "")

# THE TWO FORBIDDEN NAMES ARE ASSEMBLED FROM FRAGMENTS AT THE ASSERTION SITE, ON PURPOSE, so this
# test file does not itself embed the literals it forbids -- otherwise the gate would trip on its
# own text and would have to be weakened into uselessness.
const SEALED_DIR = join(('c', 'o', 'r', 'p', 'u', 's'))
const SEALED_FN  = "open" * "_sealed" * "_holdout"

@testset "P13 real-image arm (D-15 AMENDED, D-16)" verbose = true begin

    @testset "P13 real image ingestion" begin
        # All six committed fixtures resolve and exist. real_tif throws if one is missing, so
        # reaching a path at all is already the isfile assertion; it is restated explicitly.
        for c in P13_REAL_CONDITIONS, i in 1:3
            @test isfile(real_tif(c, i))
        end
        for c in P13_REAL_CONDITIONS
            mci = load_real(c)
            @test length(mci.data) == 2
            for d in mci.data
                @test size(d) == P13_REAL_IMSIZE
                @test d isa Matrix{Float64}
                # The frozen loader is Float64.(Images.Gray.(Images.load(path))): Rec-601 luma on
                # the on-disk RGB{N0f16}, so values are in [0, 1], neither z-scored nor
                # background-subtracted. The frames are dim; nothing sits near the top.
                @test minimum(d) >= 0.0
                @test maximum(d) <= 1.0
            end
        end

        # The 1028 x 1376 frame feeds the fixed grid UNMODIFIED. patch() truncates to the largest
        # exact multiple, so the shape is DERIVED here rather than transcribed.
        g = real_grid_provenance()
        @test g.rows_dropped == P13_REAL_GRID_TRUNCATION_ROWS
        @test g.cols_dropped == 0
        @test g.px_per_patch == prod(g.patch_size)
        @test occursin("truncates", g.truncation_note)
        for c in P13_REAL_CONDITIONS
            @test size(patch(REAL_PAIRS[c][1], 8)) == (8, 8, g.patch_size[1], g.patch_size[2])
            # Every patch holds three orders of magnitude more pixels than the 15-survivor floor
            # of _exclude_zero, so zero patches go missing on either fixture.
            @test count(ismissing, REAL_PS[c]) == 0
        end

        # THE LINEAGE CLAIM, and the reason this assertion is here rather than in a report: it
        # proves this ingestion is the SAME read chain the frozen ghat calibration used, not a
        # lookalike that happens to land nearby. spike/simulator/ghat.jl:58-64 records the
        # unmasked c1/c2 anchors as positive mu = 0.3292, negative mu = 0.2481.
        #
        # HONEST CAVEAT (recorded in real_images.jl header (b2)): ghat.jl now also records that
        # the c1/c2 figures are SUPERSEDED *for colocalization purposes* by the c2/c3 green/red
        # pair, because c1 is the DAPI/Hoechst nuclear counterstain. That does not weaken the
        # read-chain identity claim asserted here, and it does not license reinterpreting the
        # c1/c2 number as a colocalization measurement.
        #
        # BOTH SIDES ARE PINNED TO LITERALS ON PURPOSE, AND THAT IS THE POINT OF THIS BLOCK.
        # 13-D15-AMENDMENT.md CHANGE A moved P13_REAL_CHANNEL_PAIR to (2, 3) AND
        # P13_REAL_ANCHOR_MBAR to the c2/c3 values TOGETHER. Had this assertion gone on riding
        # the default pair and comparing against the constant, both sides would have moved
        # together, it would STILL HAVE PASSED, and the c1/c2 read-chain identity proof would
        # have been silently destroyed while this very comment went on claiming to test it. So
        # the pair is named explicitly and the expected values are literals.
        #
        # THIS IS A READ-CHAIN IDENTITY CHECK ON A SUPERSEDED PAIR. It is NEVER a colocalization
        # measurement and NEVER a pass condition for a colocalization claim.
        @test isapprox(real_mbar("positive"; channels = (1, 2)), 0.3292;
                       atol = P13_REAL_ANCHOR_TOL)
        @test isapprox(real_mbar("negative"; channels = (1, 2)), 0.2481;
                       atol = P13_REAL_ANCHOR_TOL)

        # THE OPERATIVE REGRESSION, on the AMENDED pair, against the amended anchors. This one
        # IS the colocalization read chain: c2/c3 = green/red, the two target proteins.
        for cond in P13_REAL_CONDITIONS
            @test isapprox(real_mbar(cond; channels = P13_REAL_CHANNEL_PAIR),
                           getproperty(P13_REAL_ANCHOR_MBAR, Symbol(cond));
                           atol = P13_REAL_ANCHOR_TOL)
        end

        # Boundary validation happens BEFORE the filesystem is touched, and the sealed substrate
        # is not reachable by asking for it under another name.
        @test_throws ArgumentError real_tif("sealed", 1)
        @test_throws ArgumentError real_tif(P13_REAL_SAMPLE, 0)
        @test_throws ArgumentError real_mbar(P13_REAL_SAMPLE; G = 4)

        # The provenance record carries the substrate-dependent knobs, not the algorithm.
        @test REAL_PROV.substrate === :real
        @test REAL_PROV.imsize == P13_REAL_IMSIZE
        @test REAL_PROV.source_zero_count == REAL_ZEROS[P13_REAL_SAMPLE]
        @test isapprox(REAL_PROV.alpha0_mbar, REAL_MBARS[P13_REAL_SAMPLE]; atol = 0.0)
    end

    @testset "P13 real alpha ladder" begin
        for c in P13_REAL_CONDITIONS
            x, y = REAL_PAIRS[c][1], REAL_PAIRS[c][2]
            M, b = REAL_MASKS[c], REAL_FLOORS[c]
            rungs = REAL_SUMS[c]
            built = REAL_LADDER[c].pairs

            # The ladder ran on the pre-registered grid, shared rung-for-rung with the simulated
            # arm, through the SAME alpha_segregate code path -- there is no substrate branch.
            @test REAL_LADDER[c].ladder == P13_REAL_ALPHA_GRID
            @test length(rungs) == length(P13_REAL_ALPHA_GRID)

            # (1) alpha = 0 reproduces the input BITWISE. Exact ==, never isapprox: an
            # approximate check would pass a transform that perturbs every pixel by a rounding
            # step, and a perturbed zero set is a CHANGED summary, because _exclude_zero keys on
            # which pixels are exactly zero.
            @test alpha_segregate(x, y, M, 0.0; b = b) == y
            @test built[1][2] == y

            # (2) alpha introduces NO NEW zero. THE ASSERTION IS DELIBERATELY NOT AN
            # ALL-PIXELS-STRICTLY-POSITIVE CHECK, and that naive form is not written anywhere in
            # this file. It is FALSE on this substrate: positive_c2 carries 2 exact-zero pixels
            # and negative_c2 carries 3, of 1 414 528, in the SOURCE files before any alpha is
            # applied -- so at alpha = 0, where the transformed channel IS the input channel, the
            # naive form would fail a CORRECT algorithm and turn this suite red for a reason that
            # has nothing to do with the transform. The correct reference is the count on the
            # UNMODIFIED channel, which is why it is recorded in provenance rather than assumed.
            #
            # MEASURED MISS AFTER 13-D15-AMENDMENT.md CHANGE A, LEFT FAILING DELIBERATELY.
            # The 2-and-3 source zeros quoted above are a property of c2, which was the
            # TRANSFORMED channel while the pair was c1/c2. On the amended pair c2/c3 the
            # transformed channel is c3, and c3 carries ZERO exact-zero pixels on one condition,
            # so this line now reads `0 > 0` and FAILS. The invariant it guards is UNAFFECTED:
            # the zero-count-preserved checks immediately below still pass on both conditions.
            # What has failed is a non-degeneracy guard asserting that the check is being
            # exercised against a non-trivial reference. It is NOT rewritten to the c3 counts,
            # because rewriting a substrate expectation to match what was measured after the
            # measurement is the act the amendment exists not to be. Reported, not adjusted.
            @test REAL_ZEROS[c] > 0
            for (k, p) in enumerate(built)
                @test count(iszero, p[2]) == REAL_ZEROS[c]
                @test rungs[k].n_zero == REAL_ZEROS[c]
            end

            # (3) total ch2 intensity is preserved, so an alpha response can never be a
            # brightness change wearing a segregation costume.
            for p in built
                @test isapprox(sum(p[2]), sum(y); rtol = P13_ALPHA_INVARIANT_RTOL)
            end

            # (4) the mask is alpha-INVARIANT: computed once from the UNMODIFIED ch1 and held
            # fixed. Recomputing it on each transformed pair is a real leak check, not a
            # tautology -- it fails the moment a future edit writes to channel 1.
            for p in built
                @test p[1] == x
                @test ch1_object_mask(p) == M
            end

            # (5) the mask fraction is inside the pre-registered band. Outside it the complement
            # cannot absorb the redistributed mass and the boost blows up.
            @test first(P13_ALPHA_MASK_FRACTION_BOUNDS) <= mean(M) <=
                  last(P13_ALPHA_MASK_FRACTION_BOUNDS)

            # (6) the realized maximum stays inside the declared dynamic range -- RECORDED here,
            # never silently clamped.
            for r in rungs
                @test r.max_value <= P13_ALPHA_MAX_VALUE_BOUND
            end

            # (7) the missing count is alpha-INVARIANT AT ZERO on this substrate, not merely
            # non-increasing. A rising count is the Pitfall-2 warning sign that pixels are being
            # DELETED rather than moved.
            for r in rungs
                @test r.n_missing == 0
            end

            # (8) m-bar is STRICTLY monotone decreasing, and the ladder crosses zero strictly
            # INSIDE the grid rather than pinned at an endpoint -- which is what makes the
            # crossing point informative. On this substrate alpha = 0 is a MODERATELY COLOCALIZED
            # pair, not a random one; the real and simulated arms are reported separately and are
            # never averaged, because their alpha = 0 semantics differ.
            #
            # MEASURED MISS AFTER 13-D15-AMENDMENT.md CHANGE A, LEFT FAILING DELIBERATELY, AND
            # THIS ONE IS A SUBSTANTIVE SCIENTIFIC FINDING RATHER THAN A BOOKKEEPING BREAK.
            # On the SUPERSEDED c1/c2 pair the ladder crossed zero inside the grid and reached
            # m-bar -0.16787 / -0.25198 at alpha = 1 -- the basis of the claim that negative
            # induced mu is CONSTRUCTIBLE from real microscopy pixels. On the AMENDED c2/c3
            # pair it does NOT cross zero at all: measured +0.460266 -> +0.138056 (positive)
            # and +0.381475 -> +0.204141 (negative). Strict monotone decrease (the line above)
            # still HOLDS; only the reach into negative m-bar is gone.
            #
            # Physically this is unsurprising: the mask is the Otsu mask of the pair's FIRST
            # channel, which was the nuclear counterstain c1 and is now the green target c2.
            # Redistributing red away from green objects, when green and red co-occur and share
            # a bright background, cannot drive the correlation below zero the way
            # redistributing green away from nuclei could.
            #
            # THIS LINE IS NOT REWRITTEN OR RELAXED. It is left failing so the retraction is
            # visible in the suite rather than only in a document. The disposition of the
            # superseded "negative induced mu is constructible from real pixels" claim is a
            # pre-registration decision, escalated by 13-17 rather than resolved in execution.
            mbars = [r.mbar for r in rungs]
            @test all(isfinite, mbars)
            @test all(mbars[i] < mbars[i - 1] for i in 2:length(mbars))
            @test first(mbars) > 0.0 > last(mbars)
            @test isapprox(first(mbars), REAL_MBARS[c]; atol = 0.0)
        end
    end

    @testset "P13 corpus seal untouched" begin
        # A source-grep gate over EVERY .jl file under spike/p13/, DISCOVERED with readdir rather
        # than enumerated, so a file added later cannot slip past the gate by not being listed.
        @test !isempty(P13_JL)
        @test "real_images.jl" in P13_JL
        for f in P13_JL
            @test !occursin(SEALED_DIR, P13_CODE[f])
            @test !occursin(SEALED_FN, P13_CODE[f])
        end

        # real_images.jl carries its own include-time self-check on its own text, so the property
        # is enforced at load, not only at test time.
        @test occursin("@__FILE__", REAL_CODE)

        # The reported runner is plan 13-16's artifact and does not exist yet. When it lands it
        # must carry the same read(@__FILE__, String) self-check, so the two phases' honesty
        # claims are verifiably consistent (the 11-10-PLAN.md:191-195 mechanism, applied one
        # phase later).
        runner = joinpath(P13_DIR, "run_p13_realimage.jl")
        if isfile(runner)
            @test occursin("@__FILE__", _strip_comments(read(runner, String)))
        else
            @test_skip "run_p13_realimage.jl is plan 13-16's artifact and has not landed yet"
        end
    end

    @testset "P13 real images are read-only" begin
        # test/test_images/ is READ-ONLY INPUT and appears in no Phase-13 files_modified list.
        # The subprocess gets its working directory via Cmd(...; dir = ...), so the running
        # process's own directory is never mutated.
        status = read(Cmd(`git status --porcelain test/test_images`; dir = P13_REPO_ROOT), String)
        @test isempty(strip(status))

        # The structural proof: git ls-files -s prints mode, blob hash, stage and path, so the
        # digest moves if any fixture's CONTENT, mode or name moves. It is not a timestamp check.
        d1 = verify_real_readonly_digest()
        d2 = verify_real_readonly_digest()
        @test d1 == d2
        @test length(d1) == 64
        @test occursin(r"^[0-9a-f]{64}$", d1)

        # No Phase-13 file mutates anything named test_images. Comment lines are stripped first,
        # so the prose that must NAME the fixtures in order to say they are read-only cannot
        # invalidate its own gate.
        mutators = r"\b(write|mv|rm|cp|touch|mkpath|chmod)\("
        for f in P13_JL, line in split(P13_CODE[f], '\n')
            if occursin("test_images", line)
                @test !occursin(mutators, line)
                @test !occursin(r"open\([^)]*\"[aw]", line)
            end
        end
    end

    @testset "P13 real arm is not a gate" begin
        # The three sentinels that make the posture machine-checked rather than merely written.
        @test P13_REAL_IS_GATED == false
        @test P13_REAL_QUALITATIVE_ONLY == true
        @test P13_REAL_READ_ONLY == true
        @test P13_ALPHA_GATED == false

        # Two substrates, one shared ladder, one shared code path.
        @test P13_ALPHA_SUBSTRATE == (:simulated, :real)
        @test P13_REAL_ALPHA_GRID == P13_ALPHA_LADDER

        # A real image has NO KNOWN registration uncertainty, so the arm reports the whole sweep
        # and headlines the widest, most conservative rung: a conservative lambda WEAKENS the
        # evidence, so a verdict that survives it is the credible one.
        @test P13_REAL_LAMBDA_HEADLINE_RULE === :widest_rung
        @test P13_REAL_LAMBDA_READS_RULE === :full_phase11_ladder
        @test P13_REAL_OOD_COMPARISON == true

        # The two honesty strings are carried verbatim, not paraphrased.
        @test !isempty(strip(P13_REAL_NAMING_CORRECTION))
        @test occursin("biological", P13_REAL_NAMING_CORRECTION)
        @test occursin("+0.3815", P13_REAL_NAMING_CORRECTION)
        @test !isempty(strip(P13_REAL_SUBSTITUTION_RECORD))
        @test occursin("Phase 16", P13_REAL_SUBSTITUTION_RECORD)

        # THE STRONGEST FORM: no real-image quantity has a pass/fail threshold at all. If no
        # threshold exists, the arm cannot accidentally become a gate. Asserted as a SOURCE
        # property of the pre-registration (comment-stripped) rather than as a list of values,
        # so a bar added later is caught even if nothing else references it.
        #
        # EXACTLY TWO EXEMPTIONS ARE ALLOWED, BOTH NAMED AND REASONED SO THE LIST CANNOT
        # SILENTLY GROW:
        #   P13_REAL_ANCHOR_TOL           -- an INGESTION REGRESSION tolerance on the frozen ghat
        #                                    read-chain identity check, not a Phase-13 bar.
        #   P13_REAL_OOD_SHIPPED_THRESHOLD -- a RECORDED reference value belonging to the SHIPPED
        #                                    net, quoted for comparison, not a Phase-13 bar.
        # MATCHED PER UNDERSCORE-SEGMENT, NOT AS A SUBSTRING, and that is a correctness fix and
        # not a loosening: a substring match reports P13_REAL_NAMING_CORRECTION as a bar because
        # "naMINg" contains MIN, which would force the honesty string onto the exemption list and
        # so make the exemption list -- the very thing that must not grow -- grow. A threshold
        # constant names its role as a whole segment.
        barwords = ("FLOOR", "MIN", "MAX", "THRESHOLD", "TOL", "BOUND", "CUTOFF")
        allowed = Set(["P13_REAL_ANCHOR_TOL", "P13_REAL_OOD_SHIPPED_THRESHOLD"])
        names = Set(String[m.match for m in eachmatch(r"P13_REAL_[A-Z0-9_]+",
                                                      CONSTS_CODE * "\n" * AMENDED_CODE)])
        barlike = filter(n -> any(seg -> seg in barwords, split(n, '_')), names)
        @test barlike == allowed
    end

    @testset "A14 transposition falsifier" begin
        # IMSIZE_SET (spike/data/seeding.jl:52-58) carries this frame size WITH THE AXES SWAPPED
        # at weight 0.03, while load_tiff yields the untransposed orientation. This one line is
        # the cheap falsifier that the axis order is statistically immaterial: the per-patch pixel
        # count and the truncation are identical either way, and only the column-major vec()
        # ordering of the grid corresponds to a transposed spatial layout. A failure here would
        # mean the real summaries sit on a permuted basis.
        for c in P13_REAL_CONDITIONS
            flipped = Vector{Matrix{Float64}}([permutedims(REAL_PAIRS[c][1]),
                                               permutedims(REAL_PAIRS[c][2])])
            @test isapprox(real_mbar_of(flipped), REAL_MBARS[c]; atol = P13_REAL_ANCHOR_TOL)
        end
    end

    @testset "P13 real arm consumes no RNG stream" begin
        # Worth asserting precisely because a reader will assume randomness is involved: every
        # other Phase-13 suite rides a fixture stream, and this one has nothing to ride. The
        # substrate is six committed files and the transform is a deterministic function of its
        # input, its mask and its floor, so there is no seed to burn and no reported stream this
        # suite could pre-observe (D-01).
        @test !occursin("rand(", REAL_CODE)
        @test !occursin("rng", REAL_CODE)
        c = P13_REAL_SAMPLE
        a = alpha_segregate(REAL_PAIRS[c][1], REAL_PAIRS[c][2], REAL_MASKS[c], 0.5;
                            b = REAL_FLOORS[c])
        b2 = alpha_segregate(REAL_PAIRS[c][1], REAL_PAIRS[c][2], REAL_MASKS[c], 0.5;
                             b = REAL_FLOORS[c])
        @test a == b2
    end

    @testset "P13 real arm ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
