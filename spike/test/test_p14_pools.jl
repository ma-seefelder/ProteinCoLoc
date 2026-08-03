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

# spike/test/test_p14_pools.jl --- the two silent failure modes, closed by RUNNING assertions.
#
# BOTH FAILURES THIS FILE GUARDS AGAINST PRODUCE PLAUSIBLE-LOOKING NUMBERS, WHICH IS WHY PROSE IS
# NOT ENOUGH AND A COMMENT IN THE MODULE IS NOT ENOUGH.
#
#   PITFALL 1 -- a CLASS-BALANCED calibration set. Split conformal's guarantee requires the
#   calibration points to be EXCHANGEABLE with the test point, and a balanced set against a
#   pi-distributed deployment stream is not exchangeable; a Bayesian false-discovery proportion
#   additionally depends on the PREVALENCE that balancing replaces. Nothing throws: the code runs,
#   the quantile computes, the coverage number looks fine and is wrong by an amount nobody can see.
#   Guarded here by (a) a comment-stripped source grep proving the balancing path is unreachable
#   from the module, (b) a realized-mass check against the MEASURED prior, and (c) feeding the
#   detector a synthetic equal-thirds vector and requiring it to throw.
#
#   PITFALL 4 -- an OOD null on the WRONG STANDARDIZER BASIS. The shipped null was fitted on the
#   Phase-7 basis; the Phase-13 net rides Phase 11's, inherited and never re-fit. A cross-basis
#   Mahalanobis score is not a wrong number, it is an uninterpretable one. Guarded here by pinning
#   `basis_provenance` and by asserting the pool the null was fitted on IS the one the trained net
#   records.
#
# AND D-06 END TO END. `p14_ood_reference -> p14_ood_state -> p14_fuse` is asserted from a MISSING
# pool all the way to `(:abstain, :ood_not_checked)`. 14-CONTEXT.md asks for the not-checked rule to
# be executable rather than prose; this is where it becomes executable, and the positive half (a
# present pool resolving to `:clear` and DECIDING) is asserted beside it so the chain cannot pass by
# abstaining unconditionally.
#
# FIXTURE SCALE ONLY, AND ON THE FIXTURE STREAM. Every draw here rides `p14_fix_rng`, never a
# reported stream: a test that pre-observed a reported draw would have consumed the very evidence a
# later runner reports. The image draws are deliberately tiny (an item costs about a second) and the
# thousands-of-items check runs on the theta-only twin, which simulates nothing.
#
# RUN IT PER FILE. The aggregate spike suite exits 1 early at the Phase-4 speedup gate and masks
# every later include block, so a green aggregate run would be no evidence at all:
#
#     julia --project=spike spike/test/test_p14_pools.jl
#
# DECOUPLING (CLAUDE.md): spike-local. Reads the Phase-13 pool and artifacts READ-ONLY, writes
# nothing, adds no package, touches no src/ byte.

using Test

# Unit under test (pulls the Phase-13 generator, the frozen basis and the in-repo OOD helpers).
isdefined(@__MODULE__, :P14_POOLS_LOADED) || include(joinpath(@__DIR__, "..", "p14", "pools.jl"))
# The D-06 chain's second and third links.
isdefined(@__MODULE__, :p14_fuse) || include(joinpath(@__DIR__, "..", "p14", "fuse.jl"))

# --- Fixtures built ONCE at module level ---------------------------------------------------------
# The expensive ones are the image draws. Everything below reuses these.

const POOLS_PATH = normpath(joinpath(@__DIR__, "..", "p14", "pools.jl"))
const P14_LANE   = normpath(joinpath(@__DIR__, "..", "p14"))

"Comment-stripped source of one file: a path MENTIONED IN PROSE must not read as a call."
_code(path) = join(filter(l -> !startswith(strip(l), "#"),
                          split(read(path, String), '\n')), '\n')

const POOLS_CODE = _code(POOLS_PATH)

# THE NEEDLES ARE ASSEMBLED FROM FRAGMENTS, the `test_p12_decoupling.jl:34-45` idiom. This file has
# to NAME the balancing helper in order to forbid it, and Phase-14 sources are themselves scanned by
# lane-wide guards; assembling the token at run time keeps the contiguous literal out of the scanned
# text while the value searched for is the full token either way.
const BALANCER_NEEDLE = "stratify" * "_by_class"
const UNSTRAT_NEEDLE  = "p13_simulate" * "_indices"
const SHIPPED_NEEDLES = ("ood_nulls_8" * ".jld2", "amended" * "_v2")
const OPERATING_POINT_NEEDLE = "youden" * "_j"

# The theta-only twin: 600 labels, no image, milliseconds. This is what the pi-distribution check
# runs on, and its agreement with the image draw is asserted in testset 1.
const FIX_CLASSES = p14_draw_classes(600; counter = P14_FIXTURE_COUNTER, rng_for = p14_fix_rng)

# Two draws at the SAME (n, counter) -- must be identical -- and one at a DIFFERENT counter on the
# same fixture seed, which must not be.
const FIX_A  = p14_draw_pool(4; counter = P14_FIXTURE_COUNTER, rng_for = p14_fix_rng)
const FIX_B  = p14_draw_pool(4; counter = P14_FIXTURE_COUNTER, rng_for = p14_fix_rng)
const FIX_C2 = p14_draw_pool(1; counter = P14_CAL_COUNTER,     rng_for = p14_fix_rng)

# The net's own record, read independently of the unit under test.
const NET_META = p14_net_meta()

# The OOD reference on the real pool, and the same call against a directory that does not exist.
const OOD_REF   = p14_ood_reference()
const OOD_BOGUS = p14_ood_reference(pool_dir = joinpath(@__DIR__, "definitely-not-a-pool-dir"))

# The read-only proof: the pool listing digest, taken BEFORE the misspecification draw and compared
# again at the end. The `digest_before`/`digest_after` idiom of the Phase-13 real-image runner,
# applied to a gitignored tree that `git status` cannot speak about at all.
const POOL_DIGEST_BEFORE = p14_pool_provenance(OOD_REF).pool_listing_sha256

# One matched misspecification pair, so the arm is actually EXECUTED rather than merely defined.
const MIS = p14_misspec_pool(1; counter = P14_FIXTURE_COUNTER, rng_for = p14_fix_rng,
                             family = :noise)

# A verdict object for the D-06 chain. Built by hand: this file never runs a detector, it only
# resolves a state from a null that either does or does not carry an operating point.
const VERDICT_QUIET = OODVerdict(0.0, false, NamedTuple())

# Binomial standard error, written out HERE rather than borrowed from the unit under test: a test
# that computes its expectation with the code under test proves only that one function agrees with
# itself.
_se(p, n) = sqrt(p * (1.0 - p) / n)

@testset "P14 pools: unstratified draws, the Phase-13 basis null, and D-06 end to end" verbose = true begin

    @testset "draws are reproducible and counter-separated" begin
        # Counter-based Philox: the same (n, counter, family) is the same draw, bitwise, on any
        # thread count and any resume. Compared field by field rather than by `===`, because two
        # separately-allocated identical items are not the same object.
        @test length(FIX_A) == 4 == length(FIX_B)
        @test all(j -> FIX_A[j] == FIX_B[j], eachindex(FIX_A))
        @test all(j -> FIX_A[j].Zs == FIX_B[j].Zs && FIX_A[j].Zc == FIX_B[j].Zc, eachindex(FIX_A))

        # A different counter is a different sub-stream. If this ever passes trivially, two
        # activities are sharing draws.
        @test FIX_C2[1].Zs != FIX_A[1].Zs

        # AND THE ITEMS WITHIN ONE DRAW ARE DISTINCT. A stream constructor that ignored the item
        # index would hand every item the same RNG and return `n` copies of one item -- which is
        # silent, because identical items are perfectly valid items and every reproducibility and
        # counter-separation assertion above still passes on them.
        @test length(unique(x -> x.Zs, FIX_A)) == length(FIX_A)
        @test length(unique(x -> x.rho_s, FIX_A)) == length(FIX_A)

        # THE FIXTURE STREAM DIFFERS FROM THE REPORTED STREAM AT THE SAME COUNTER. Asserted on the
        # Philox KEY WORD, which is the quantity 14-01 proved disjoint -- two seeds that differ are
        # not enough if `seed xor salt xor counter` can still collide -- and then confirmed on an
        # actual draw.
        kw_fix = p14_datagen_key_word(P14_FIXTURE_COUNTER; rng_for = p14_fix_rng)
        kw_dev = p14_datagen_key_word(P14_FIXTURE_COUNTER; rng_for = p14_rng)
        @test kw_fix != kw_dev
        @test rand(p14_datagen_rng(P14_FIXTURE_COUNTER; rng_for = p14_fix_rng)(1)) !=
              rand(p14_datagen_rng(P14_FIXTURE_COUNTER; rng_for = p14_rng)(1))

        # THE C-04 FRESHNESS PROPERTY, ASSERTED RATHER THAN ASSUMED. Index-keyed generation makes a
        # same-configuration pool a SUPERSET of the training pool, so the ONLY thing that makes a
        # Phase-14 draw genuinely fresh is a different key family.
        for c in P14_ALL_COUNTERS, fam in (p14_rng, p14_fix_rng)
            kw = p14_datagen_key_word(c; rng_for = fam)
            @test kw in Set(_p14_key_words())
            @test !(kw in Set(_p14_p13_key_words()))
        end

        # And the guard FALSIFIES: a foreign RNG family and an unreserved counter are both refused,
        # because either one would key a stream whose disjointness was never proved.
        @test_throws ArgumentError p14_datagen_key_word(P14_FIXTURE_COUNTER; rng_for = identity)
        @test_throws ArgumentError p14_datagen_key_word(777; rng_for = p14_fix_rng)

        # The theta-only twin is the SAME draw at the same indices, not an approximation of it.
        @test [x.class for x in FIX_A] == FIX_CLASSES[1:length(FIX_A)]
        # And the item really carries what the plan says it carries.
        @test all(x -> haskey(x, :Zs) && haskey(x, :Zc) && haskey(x, :lambda) &&
                       haskey(x, :theta) && haskey(x, :class), FIX_A)
        @test all(x -> haskey(x.theta, :sample) && haskey(x.theta, :control), FIX_A)
        @test all(x -> x.theta.sample.ρ_true == x.rho_s && x.theta.control.ρ_true == x.rho_c, FIX_A)
    end

    @testset "draws are pi-distributed, not equal thirds (Pitfall 1)" begin
        expected = NET_META.pi_class_masses
        realized = class_masses(FIX_CLASSES)
        n = length(FIX_CLASSES)

        # The expectation is computed HERE, from the net's own recorded prior, and every mass is
        # read BY NAME -- five class orderings coexist in this codebase and an inconsistent
        # relabelling produces numbers that look right.
        for k in (:exclusion, :random, :coloc)
            @test abs(getproperty(realized, k) - getproperty(expected, k)) <=
                  3.0 * _se(getproperty(expected, k), n)
        end

        # THE DRAW IS NOT BALANCED, and this is the assertion that would fail if a future edit
        # routed the draw through the accept/reject balancer.
        @test !all(k -> abs(getproperty(realized, k) - 1/3) <= 3.0 * _se(1/3, n),
                   (:exclusion, :random, :coloc))
        @test abs(realized.coloc - 1/3) > 3.0 * _se(1/3, n)
        @test abs(realized.random - 1/3) > 3.0 * _se(1/3, n)

        # The unit under test agrees, and returns the masses so a runner can persist them.
        m = p14_assert_unstratified(FIX_CLASSES)
        @test m == realized

        # AND IT FAILS LOUDLY ON A BALANCED SET. This is the falsification: without it, the
        # detector could be a function that returns its argument.
        balanced = vcat(fill(COLOC, 334), fill(RANDOM, 333), fill(EXCLUSION, 333))
        @test_throws ErrorException p14_assert_unstratified(balanced)
        # AND IT REPORTS THE CAUSE, NOT THE SYMPTOM. A balanced set also deviates from pi, so a
        # detector that only checked pi would still throw -- with a message that sends the reader
        # looking for a broken stream instead of a balanced one.
        msg = try
            p14_assert_unstratified(balanced)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("CLASS-BALANCED", msg)

        # The pi-agreement check is INDEPENDENT of the balance check: a set that is neither
        # balanced nor pi-distributed must also fail, or the second check is decoration.
        @test_throws ErrorException p14_assert_unstratified(fill(COLOC, 300))
        @test_throws ErrorException p14_assert_unstratified(
            vcat(fill(COLOC, 200), fill(RANDOM, 100), fill(EXCLUSION, 300)))
    end

    @testset "no stratification path is reachable" begin
        # Comment-stripped, so the module may EXPLAIN the forbidden path in prose while the
        # executable text may not contain it. A docstring is CODE to this grep, which is exactly
        # the trap this project has now hit repeatedly.
        @test !occursin(BALANCER_NEEDLE, POOLS_CODE)
        # The positive half: the ban is not satisfied vacuously by a module that draws nothing.
        @test occursin(UNSTRAT_NEEDLE, POOLS_CODE)
        # The operating point is the pre-registered in-distribution quantile, never the post-hoc
        # labelled-data maximum along a ROC.
        @test !occursin(OPERATING_POINT_NEEDLE, POOLS_CODE)
        @test occursin("id_threshold(", POOLS_CODE)
        # The shipped bundle is not read here at all; if a runner ever reads it, it is a labelled
        # comparison reference and the two thresholds are never fused.
        @test all(nd -> !occursin(nd, POOLS_CODE), SHIPPED_NEEDLES)
    end

    @testset "the OOD reference is on the Phase-13 basis (Pitfall 4)" begin
        @test OOD_REF.available
        @test OOD_REF.thr isa Real && isfinite(OOD_REF.thr)
        @test OOD_REF.n_pool > 0
        @test OOD_REF.n_shards > 0
        # THE POOL IS THE NET'S OWN. This is the whole of Pitfall 4: the null must be fitted on the
        # basis the net rides, and the only way to know which that is, is to read it off the net.
        @test normpath(OOD_REF.pool_dir) == normpath(NET_META.pool_dir)
        @test OOD_REF.basis_provenance === :p13_pool_phase11_zt
        @test occursin("INHERITED", NET_META.zt_provenance)
        # Both halves of every pair were used: the pool holds n_pool/2 pairs.
        @test OOD_REF.n_pool == 2 * (OOD_REF.n_pool ÷ 2)
        # The pre-registered operating quantile, not a chosen one.
        @test OOD_REF.q === OOD_ID_QUANTILE
        # The unwired channels are RECORDED rather than implied, so `:clear` downstream cannot be
        # read as "three channels say clear".
        @test OOD_REF.channels_wired == (:density,)
        @test OOD_REF.channels_not_wired == (:pp, :noise)
        # `:thr` is present exactly when the fit succeeded.
        @test haskey(OOD_REF.ood_nulls, :thr)
        @test OOD_REF.ood_nulls.thr === OOD_REF.thr
        # The provenance block a runner embeds points at the pool it actually read.
        prov = p14_pool_provenance(OOD_REF)
        @test prov.pool_n_columns == OOD_REF.n_pool
        @test prov.pool_bytes > 0
        @test prov.ood_basis_provenance === :p13_pool_phase11_zt
    end

    @testset "a missing pool degrades honestly, and D-06 fires end to end" begin
        # NO THROW, NO SUBSTITUTE. A missing operating point must leave the flag honestly inert.
        @test OOD_BOGUS.available == false
        @test OOD_BOGUS.thr === nothing
        @test OOD_BOGUS.nulls === nothing
        @test !haskey(OOD_BOGUS.ood_nulls, :thr)
        @test OOD_BOGUS.note isa AbstractString && !isempty(strip(OOD_BOGUS.note))
        # A null that was never fitted has NO basis, and claiming one would be Pitfall 4 in
        # miniature.
        @test OOD_BOGUS.basis_provenance === :unavailable

        # THE FULL D-06 CHAIN, from an absent pool to a silence with a reason attached.
        state = p14_ood_state(OOD_BOGUS.ood_nulls, VERDICT_QUIET)
        @test state === :not_checked
        fused = p14_fuse(ood_state = state, conformal_status = :singleton, disagree = false)
        @test fused.action === :abstain
        @test fused.reason === :ood_not_checked
        @test (fused.action, fused.reason) === (:abstain, :ood_not_checked)

        # THE POSITIVE HALF, so the chain cannot pass by abstaining unconditionally: with the real
        # reference the same quiet verdict resolves to `:clear` and the rule DECIDES.
        @test p14_ood_state(OOD_REF.ood_nulls, VERDICT_QUIET) === :clear
        decided = p14_fuse(ood_state = p14_ood_state(OOD_REF.ood_nulls, VERDICT_QUIET),
                           conformal_status = :singleton, disagree = false)
        @test decided.action === :decide
    end

    @testset "the gate set is labelled development-only (Pitfall 3)" begin
        g = p14_gate_dev_substrate()
        @test g.reported === false
        @test g.role === :development_and_crosscheck_only
        @test g.note isa AbstractString && !isempty(strip(g.note))
        @test isfile(g.path)
        # The labels are FIELDS, which is what stops a runner from consuming the set silently.
        @test haskey(g, :reported) && haskey(g, :role)
        # It really is the Phase-13 gate report and not an empty stand-in.
        @test haskey(g.report, "P13_TAU")
    end

    @testset "cache roots are separated and nothing was written to the Phase-13 pool" begin
        p13_root = normpath(joinpath(@__DIR__, "..", "data", "cache", "p13"))
        p14_root = normpath(p14_cache_root())
        @test p13_root != p14_root
        # Boundary-safe both ways: neither root may CONTAIN the other, so a recursive write under
        # one can never land in the other.
        sep = Base.Filesystem.path_separator
        @test !startswith(p13_root, p14_root * sep)
        @test !startswith(p14_root, p13_root * sep)
        @test basename(p14_root) == "p14"

        # NO PHASE-14 SOURCE WRITES UNDER THE PHASE-13 CACHE. Lane-wide (a glob, never a
        # hand-maintained list) and a CONJUNCTION of a path token and a write call, so that READING
        # `meta.pool_dir` -- which the null fit requires -- stays legal.
        cache_token = "\"cache\", " * "\"p13\""
        writers     = ("jldsave", "mkpath", "save_pool")
        hits = String[]
        for f in sort(filter(x -> endswith(x, ".jl"), readdir(P14_LANE)))
            for ln in split(_code(joinpath(P14_LANE, f)), '\n')
                s = strip(ln)
                occursin(cache_token, s) || continue
                (any(w -> occursin(w, s), writers) ||
                 (occursin("open(", s) && occursin("\"w\"", s))) && push!(hits, "$f: $s")
            end
        end
        @test hits == String[]
        # This module writes nothing AT ALL -- not into p13's root, not into its own.
        @test !any(w -> occursin(w, POOLS_CODE), ("jldsave", "mkpath", "save_pool"))

        # THE READ-ONLY PROOF, on a gitignored tree that `git status` cannot speak about: the
        # digest is recomputed FROM DISK here, after the null fit and the misspecification draw,
        # and every shard is still present at its original size.
        @test p14_pool_provenance(OOD_REF).pool_listing_sha256 == POOL_DIGEST_BEFORE
    end

    @testset "the misspecification arms are matched and pre-registered" begin
        @test length(MIS.id) == length(MIS.ood) == MIS.n == 1
        @test MIS.family === :noise
        # `:strongest` is the top rung of the EXISTING grid, not a new one.
        @test MIS.level == OOD_GRID_LEVELS
        # MATCHED: same parameters, same size, same label. Only the generator differs.
        @test MIS.id[1].lambda == MIS.ood[1].lambda
        @test MIS.id[1].imsize == MIS.ood[1].imsize
        @test MIS.id[1].class === MIS.ood[1].class
        @test MIS.id[1].Zs != MIS.ood[1].Zs
        @test all(isfinite, MIS.ood[1].Zs) && all(isfinite, MIS.ood[1].Zc)
        @test length(MIS.ood[1].Zs) == length(MIS.id[1].Zs)
        # No new misspecification may be authored under an existing name.
        @test_throws ArgumentError p14_misspec_pool(1; counter = P14_FIXTURE_COUNTER,
                                                    rng_for = p14_fix_rng, family = :not_a_family)
    end

    @testset "P14 pool machinery ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
