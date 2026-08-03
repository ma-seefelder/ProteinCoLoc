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

# spike/test/test_p14_consts.jl --- the Phase-14 Tier-1 pre-registration gate (D-07).
#
# WHAT THIS FILE IS FOR. `spike/p14/consts.jl` is the scientific-integrity artifact of
# Phase 14: every seed, counter, sample size, alpha and bar that a reported Phase-14
# number will be scored against. Its value is entirely in the fact that it was committed
# FIRST. This suite asserts every locked value AS A LITERAL, so a post-hoc edit -- the one
# failure mode that would silently invalidate every downstream number -- breaks the suite
# loudly instead of passing quietly.
#
# It also re-runs the seed-disjointness proof against the RECOMPUTED Phase-13 inventory (a
# comment naming a derived seed is not evidence), asserts the FULL Philox key-word
# cross-product against the whole Phase-13 family, and proves that tau is NOT here -- D-07
# requires it to be inherited from the Phase-13 artifact by loading, and a Phase-14 literal
# for it is the exact divergence D-07 exists to prevent.
#
# Mirrors test_p13_consts.jl / test_sbc.jl / test_bf.jl: license header -> using Test ->
# guarded include of the unit under test -> source read ONCE at module level -> one outer
# @testset. No training, no simulation, no figure call -- this gate is milliseconds.
#
# RUN IT DIRECTLY. `spike/test/runtests.jl` is NOT a gate for this phase: it exits 1 early
# at the Phase-4 SPEEDUP_GATE, masking every later include block.
#     julia --project=spike spike/test/test_p14_consts.jl

using Test

# Unit under test: the Tier-1 pre-registration. Guarded so a re-include under the harness
# is a silent no-op. Including it transitively includes spike/p13/consts.jl (also guarded),
# which is what puts the Phase-13 seeds, counters and the isolated `_GC` gate module in
# scope for the disjointness assertions below.
isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) ||
    include(joinpath(@__DIR__, "..", "p14", "consts.jl"))

# The pre-registration source, read ONCE (the test_bf.jl:67 / test_p13_consts.jl:44 idiom),
# for the source-grep assertions below. COMMENT LINES ARE STRIPPED before any count or
# `occursin` assertion, so the file's own explanatory prose -- which necessarily NAMES the
# forbidden forms in order to forbid them -- cannot invalidate its own gate. The RAW source
# is kept too, because two assertions below are deliberately about the prose itself.
const P14_CONSTS_SRC = read(joinpath(@__DIR__, "..", "p14", "consts.jl"), String)
const P14_CONSTS_CODE = join(
    filter(l -> !startswith(strip(l), "#"), split(P14_CONSTS_SRC, '\n')), '\n')

@testset "P14 Tier-1 pre-registration (D-07)" verbose = true begin

    @testset "the locked literals" begin
        # A CHANGED CONSTANT IS A CHANGED EXPERIMENT. Every frozen value is asserted BY
        # VALUE, not merely by existence, so an edit after a result breaks this file rather
        # than quietly re-scoring the phase.
        @test P14_DEV_SEED == 0x0000_0000_0B14_DE71
        @test P14_FIX_SEED == 0x0000_0000_0B14_F1F7
        @test P14_SALT     == 0x5851_F42D_4C95_7F2D

        # Counters, one per activity.
        @test P14_CAL_COUNTER     == 1
        @test P14_EVAL_COUNTER    == 2
        @test P14_OOD_COUNTER     == 3
        @test P14_REAL_COUNTER    == 4
        @test P14_FIXTURE_COUNTER == 99

        # Sample sizes. At n = 2000 and alpha = 0.10 the realized conformal coverage has
        # sd = sqrt(0.10 * 0.90 / 2000) = 0.0067; that arithmetic is exact and is what pins
        # the SC1-d band, so n is not a free knob after the fact.
        @test P14_N_CAL  == 2000
        @test P14_N_EVAL == 2000
        @test P14_N_OOD  == 1000
        @test isapprox(sqrt(P14_ALPHA_CONFORMAL * (1 - P14_ALPHA_CONFORMAL) / P14_N_EVAL),
                       0.0067; atol = 1e-4)

        # The two alphas and the two swept grids.
        @test P14_ALPHA_CONFORMAL == 0.10
        @test P14_ALPHA_FDR_GRID == (0.01, 0.05, 0.10, 0.20)
        @test P14_PI_COLOC_GRID == (0.05, 0.10, 0.217, 0.40, 0.70)

        # THE FIVE JUDGEMENT-CALL BARS, at the values the user ruled on 2026-08-03
        # (accept-as-proposed). These are the numbers the ruling froze; if any one of them
        # is missed by a result, it is REPORTED and NOT re-tuned, and this test is what
        # makes a quiet re-tune impossible.
        @test P14_SKILL_FLOOR      == 0.60
        @test P14_SPEARMAN_FLOOR   == 0.95
        @test P14_COVERAGE_FLOOR   == 0.20
        @test P14_AUC_HARD_FLOOR   == 0.80
        @test P14_OOD_MARGIN_FLOOR == 0.50

        # A LABEL, never a verdict (14-RESEARCH Pitfall 9).
        @test P14_AMBIGUOUS_RATE_FLOOR == 0.01

        # The class-order contract. `===` because the identity of the tuple of Symbols is
        # the contract: five different class orderings coexist in this codebase and every
        # Phase-14 read site must be BY NAME.
        @test P14_CLASS_KEYS === (:coloc, :random, :exclusion)
    end

    @testset "seeds are fresh and disjoint" begin
        @test !(UInt64(P14_DEV_SEED) in _p14_forbidden())
        @test !(UInt64(P14_FIX_SEED) in _p14_forbidden())
        @test UInt64(P14_DEV_SEED) != UInt64(P14_FIX_SEED)

        # A repeated (seed, salt) pair is a repeated Philox key, hence a repeated stream.
        @test !(P14_SALT in P14_REPO_SALTS)
        # P13_SALT is in the Phase-14 salt inventory even though it is not in P13_REPO_SALTS:
        # Phase 13 could not forbid its own salt to itself.
        @test P13_SALT in P14_REPO_SALTS
        @test length(P14_REPO_SALTS) == length(P13_REPO_SALTS) + 1

        # ONE @test PER RESERVED SEED, written out rather than looped, so a failure names
        # WHICH stream was collided with instead of reporting an anonymous iteration.
        @test UInt64(P14_DEV_SEED) != UInt64(VAL_MASTER_SEED)
        @test UInt64(P14_DEV_SEED) != UInt64(NPE_MASTER_SEED)
        @test UInt64(P14_DEV_SEED) != UInt64(VAL_FIX_SEED)
        @test UInt64(P14_DEV_SEED) != UInt64(RATIO_PAIR_SEED)
        @test UInt64(P14_DEV_SEED) != UInt64(CORPUS_MASTER_SEED)
        @test UInt64(P14_DEV_SEED) != UInt64(DEFAULT_MASTER_SEED)
        @test UInt64(P14_DEV_SEED) != UInt64(P11_DEV_SEED)
        @test UInt64(P14_DEV_SEED) != UInt64(P13_DEV_SEED)
        @test UInt64(P14_DEV_SEED) != UInt64(P13_FIX_SEED)

        @test UInt64(P14_FIX_SEED) != UInt64(VAL_MASTER_SEED)
        @test UInt64(P14_FIX_SEED) != UInt64(NPE_MASTER_SEED)
        @test UInt64(P14_FIX_SEED) != UInt64(VAL_FIX_SEED)
        @test UInt64(P14_FIX_SEED) != UInt64(RATIO_PAIR_SEED)
        @test UInt64(P14_FIX_SEED) != UInt64(CORPUS_MASTER_SEED)
        @test UInt64(P14_FIX_SEED) != UInt64(DEFAULT_MASTER_SEED)
        @test UInt64(P14_FIX_SEED) != UInt64(P11_DEV_SEED)
        @test UInt64(P14_FIX_SEED) != UInt64(P13_DEV_SEED)
        @test UInt64(P14_FIX_SEED) != UInt64(P13_FIX_SEED)

        # PROOF THAT THE RECOMPUTE PATH RAN, not that a comment claims it did. The Phase-14
        # forbidden set is DERIVED from `_p13_forbidden()`, which itself RECOMPUTES both
        # frozen ship-gate seed families rather than trusting a literal.
        @test any(v -> UInt64(v) in _p14_forbidden(), values(_GC.PROD_SEED_V2))
        @test any(v -> UInt64(v) in _p14_forbidden(), values(_GC.PROD_SEED))
        # Phase 13's own two seeds are in the set BECAUSE Phase 13 could not forbid itself.
        @test UInt64(P13_DEV_SEED) in _p14_forbidden()
        @test UInt64(P13_FIX_SEED) in _p14_forbidden()
        # The Phase-14 set is a strict superset of the Phase-13 one.
        @test issubset(_p13_forbidden(), _p14_forbidden())
        # A duplicate would silently shrink the forbidden set.
        @test length(_p14_forbidden()) == length(_p14_forbidden_list())
        # The forbidden inventory is DERIVED, not retyped: the include is what supplies it.
        @test occursin("p13", P14_CONSTS_CODE)
    end

    @testset "streams are counter-separated" begin
        # The draw-a-number form is the executable proof of stream separation: distinct
        # reserved counters must not produce the same first draw.
        @test rand(p14_rng(P14_CAL_COUNTER), UInt64) != rand(p14_rng(P14_EVAL_COUNTER), UInt64)
        @test rand(p14_rng(P14_EVAL_COUNTER), UInt64) != rand(p14_rng(P14_OOD_COUNTER), UInt64)
        @test rand(p14_rng(P14_OOD_COUNTER), UInt64) != rand(p14_rng(P14_REAL_COUNTER), UInt64)
        # The fixture stream never overlaps the reported one, EVEN AT THE SAME COUNTER --
        # which is the only case that matters, because it is the one a counter convention
        # cannot protect against.
        @test rand(p14_fix_rng(P14_FIXTURE_COUNTER), UInt64) !=
              rand(p14_rng(P14_FIXTURE_COUNTER), UInt64)
        # Reproducibility: same construction, same draw. Without this the disjointness
        # assertions above would be satisfiable by an RNG that is merely nondeterministic.
        @test rand(p14_rng(P14_CAL_COUNTER), UInt64) == rand(p14_rng(P14_CAL_COUNTER), UInt64)
        @test rand(p14_fix_rng(P14_FIXTURE_COUNTER), UInt64) ==
              rand(p14_fix_rng(P14_FIXTURE_COUNTER), UInt64)
        # Fixtures must never ride a reported counter.
        @test P14_FIXTURE_COUNTER ∉ (P14_CAL_COUNTER, P14_EVAL_COUNTER,
                                     P14_OOD_COUNTER, P14_REAL_COUNTER)
        @test P14_FIXTURE_COUNTER ∉ P14_REPORTED_COUNTERS
        @test length(Set(P14_ALL_COUNTERS)) == length(P14_ALL_COUNTERS)
    end

    @testset "the Philox key-word cross-product is pairwise distinct vs the whole P13 family" begin
        # TWO SEEDS THAT DIFFER ARE NOT ENOUGH. What actually keys Philox in this project's
        # runners is `seed xor salt xor counter` (spike/p13/run_three_way_gate.jl:155-165),
        # so THAT is the quantity that must be proven distinct -- across the full
        # cross-product, including both bare key words on each side. A shared key between a
        # Phase-13 reported stream and the Phase-14 evaluation set would mean Phase 14 scored
        # itself on draws a Phase-13 number already rode, undetectably from the artifacts.
        kw14 = _p14_key_words()
        kw13 = _p14_p13_key_words()
        @test length(kw14) == 2 * length(P14_ALL_COUNTERS) + 2
        @test length(kw13) == 2 * length(P14_P13_COUNTERS) + 2
        @test isempty(intersect(Set(kw14), Set(kw13)))
        for a in kw14, b in kw13
            @test a != b
        end
        # ...and pairwise distinct WITHIN Phase 14, so two Phase-14 activities cannot share
        # a key either.
        @test length(Set(kw14)) == length(kw14)

        # The draw-a-number cross-check against the Phase-13 REPORTED gate stream.
        @test rand(p14_rng(P14_EVAL_COUNTER), UInt64) != rand(p13_rng(P13_GATE_COUNTER), UInt64)
        @test rand(p14_rng(P14_CAL_COUNTER), UInt64) != rand(p13_rng(P13_DATAGEN_COUNTER), UInt64)
        @test rand(p14_fix_rng(P14_FIXTURE_COUNTER), UInt64) !=
              rand(p13_fix_rng(P13_FIXTURE_COUNTER), UInt64)
    end

    @testset "tau is NOT Tier-1 here (D-07)" begin
        # TAU IS INHERITED BY LOADING, NOT DECLARED. D-07: Phase 14 loads tau from the
        # Phase-13 artifact with its provenance asserted, and never re-derives it and never
        # copies its literal into the Phase-14 namespace. Copying the value is how two
        # numbers silently diverge; re-deriving guarantees a second, different cut on the
        # same hypothesis space. A Phase-14 literal for tau is THE divergence D-07 exists to
        # prevent, so it is forbidden twice: by binding and by source text.
        @test !isdefined(@__MODULE__, :P14_TAU)
        @test !occursin("P14_TAU", P14_CONSTS_CODE)
        # The Phase-13 value IS in scope, and that is the point -- it is read from the one
        # place it was measured rather than mirrored.
        @test isdefined(@__MODULE__, :P13_TAU)
        @test !occursin("const P13_TAU", P14_CONSTS_CODE)
        # Phase 14 must not have re-declared the Phase-13 probe provenance either; the
        # provenance travels with the artifact, not with a copy of it.
        @test !occursin("P14_TAU_PROBE", P14_CONSTS_SRC)
    end

    @testset "the two alphas are distinct bindings (D-07)" begin
        # They SHARE the value 0.10, which is precisely why neither may ever be assigned
        # from the other: such an assignment is invisible to every test that checks values.
        @test isdefined(@__MODULE__, :P14_ALPHA_CONFORMAL)
        @test P14_ALPHA_CONFORMAL in P14_ALPHA_FDR_GRID
        # Neither spelling of alpha_FDR may be assigned FROM the conformal level...
        @test !occursin(r"\balpha_fdr\b\s*=\s*P14_ALPHA_CONFORMAL"i, P14_CONSTS_CODE)
        @test !occursin(r"\bP14_ALPHA_FDR\w*\s*=\s*P14_ALPHA_CONFORMAL", P14_CONSTS_CODE)
        # ...nor the conformal level from an FDR name, in either direction.
        @test !occursin(r"\bP14_ALPHA_CONFORMAL\s*=\s*\w*ALPHA_FDR"i, P14_CONSTS_CODE)
        # alpha_FDR is a USER PARAMETER, so what is frozen is only the SWEEP grid, and that
        # grid must still span more than one level or it cannot demonstrate control ACROSS
        # levels.
        @test length(P14_ALPHA_FDR_GRID) > 1
        @test !isdefined(@__MODULE__, :P14_ALPHA_FDR)
    end

    @testset "the judgement-call bars are declared as such" begin
        @test length(P14_JUDGEMENT_CALL_BARS) == 5
        @test all(s -> isdefined(@__MODULE__, s), P14_JUDGEMENT_CALL_BARS)
        @test all(s -> getfield(@__MODULE__, s) isa Real, P14_JUDGEMENT_CALL_BARS)
        @test Set(P14_JUDGEMENT_CALL_BARS) == Set((:P14_SKILL_FLOOR, :P14_SPEARMAN_FLOOR,
                                                   :P14_COVERAGE_FLOOR, :P14_AUC_HARD_FLOOR,
                                                   :P14_OOD_MARGIN_FLOOR))
        # THE LABEL MUST TRAVEL WITH THE FILE. Asserted on the RAW source, because the label
        # lives in the prose the comment-stripping removes -- which is the whole point: a
        # reader who opens this file must meet the words before the numbers.
        @test occursin("JUDGEMENT CALL", P14_CONSTS_SRC)
        # The user's ruling and its date are on the record at the point of definition.
        @test occursin("USER RULING", P14_CONSTS_SRC)
        @test occursin("2026-08-03", P14_CONSTS_SRC)
        # ...and so is the binding consequence, un-softened.
        @test occursin("NOT RE-TUNED", uppercase(P14_CONSTS_SRC))
    end

    @testset "the iteration allowance is pre-declared and non-trivial" begin
        # ONE documented iteration, authorised in advance, with ONE named trigger. A second
        # is not authorised by that file and cannot be authorised by amending it.
        @test P14_ITERATION_ALLOWANCE == 1
        @test P14_ITERATION_TRIGGER isa AbstractString
        @test !isempty(strip(P14_ITERATION_TRIGGER))
        # The remedy is NAMED, so "we iterated" cannot later mean something else.
        @test occursin("APS", P14_ITERATION_TRIGGER)
        @test occursin("LAC", P14_ITERATION_TRIGGER)
        # The allowance concerns conformal coverage and explicitly NOT the five bars.
        @test occursin("P14_JUDGEMENT_CALL_BARS", P14_ITERATION_TRIGGER)
    end

    @testset "declared deviations and decoupling are on the record" begin
        @test P14_SRC_UNTOUCHED == true
        @test length(P14_DECLARED_DEVIATIONS) == 4
        @test all(d -> d isa AbstractString && !isempty(strip(d)), P14_DECLARED_DEVIATIONS)
        # The four amendments this phase runs under, named before any result exists.
        @test any(d -> occursin("D-02", d), P14_DECLARED_DEVIATIONS)
        @test any(d -> occursin("D-05", d), P14_DECLARED_DEVIATIONS)
        @test any(d -> occursin("D-07", d), P14_DECLARED_DEVIATIONS)
        @test any(d -> occursin("D-03a", d), P14_DECLARED_DEVIATIONS)
        # The guard sentinel is the declared-deviation tuple, NOT a seed. A seed is the
        # worst possible sentinel here: every Phase-N pre-registration must re-declare prior
        # phases' seeds to assert disjointness, and this file mirrors P13_DEV_SEED below.
        @test occursin("isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS)", P14_CONSTS_CODE)
        @test occursin("P13_DEV_SEED", P14_CONSTS_CODE)
    end

    @testset "P14 ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
