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

# spike/test/test_p14_fuse.jl --- the D-05 asymmetric fusion and D-06's three-valued OOD
# (SC2-a, SC2-b).
#
# THE EXPECTATION TABLE BELOW IS WRITTEN OUT LITERALLY, ALL 36 ROWS, AND THAT IS THE WHOLE POINT.
# A test that computes its expectation by re-implementing the rule proves only that two copies of
# the same mistake agree. Every cell here was worked out by hand from D-05 and D-06 and typed out;
# the unit under test is never consulted to produce an expected value. The table is then asserted
# to BE the complete 3 x 3 x 2 x 2 product -- no cell missing, no cell twice -- so it cannot pass by
# quietly omitting the row that hurts.
#
# THE ONE CELL THAT IS THE WHOLE AMENDMENT:
#
#     (:clear, :singleton, disagree = true, allow_unchecked_ood = false)
#         =>  (:decide, :decided_contra_classical)
#
# ROADMAP SC2 as originally written is `OOD OR cross-method disagreement OR ambiguous set` => be
# silent. Under that literal OR this cell would ABSTAIN. D-05 AMENDS it, because Phase 9 exists so
# that v2.0 can be positioned as "knows when the classics are wrong, not merely agrees with them" --
# so the literal OR would silence the tool in exactly the cases that demonstrate its headline
# differentiator, and three OR'd triggers compound into a tool that is usually silent. Disagreeing
# with Costes/Manders is the PRODUCT THESIS, not a failure mode. It therefore DECIDES, and the
# disagreement is carried out on the result as a named reason rather than being discarded.
#
# THAT AMENDMENT IS ASSERTED, NOT DESCRIBED. Reverting branch 5 of `p14_fuse` to SC2's literal
# `disagree` alone must make this file exit 1, and the plan requires that scratch edit to be RUN.
#
# THE SECOND HALF IS D-06, AND THE FAILURE IT PREVENTS IS SILENT AND SCIENTIFIC. `ood_verdict`
# (src/amortized/ood.jl:394-396) computes `flag = thr === nothing ? false : fused > thr`, so when
# the threshold is absent it returns a FALSE flag FOR EVERY INPUT, UNCONDITIONALLY. A false flag is
# genuinely two-valued at the source, and `src/amortized/local_map.jl:70-72` states outright that
# such a flag means NOT CHECKED, not in distribution. Reading it as "in distribution" makes the tool
# speak confidently in exactly the regime where it has no evidence it should. So `p14_ood_state` is
# THREE-valued, never a Bool, `:not_checked` ABSTAINS by default, and the opt-out is a keyword a
# caller has to set deliberately.
#
# Mirrors spike/test/test_p14_fdr.jl: license header -> using Test -> guarded include of the unit
# under test -> fixtures computed ONCE at module level -> one outer @testset. No training, no
# simulation, no random draw at all: every input here is a symbol or a hand-built struct.
#
# RUN IT PER FILE. The aggregate spike suite exits 1 early at the Phase-4 speedup gate and masks
# every later include block, so a green aggregate run would be no evidence at all:
#
#     julia --project=spike spike/test/test_p14_fuse.jl

using Test

# Unit under test (pulls spike/p14/consts.jl and, read-only, src/results.jl transitively).
isdefined(@__MODULE__, :p14_fuse) || include(joinpath(@__DIR__, "..", "p14", "fuse.jl"))

# --- Fixtures built ONCE at module level ---------------------------------------------------

const FUSE_FIX_OOD_STATES = (:fired, :clear, :not_checked)
const FUSE_FIX_CONFORMAL  = (:singleton, :ambiguous, :empty)
const FUSE_FIX_DISAGREE   = (false, true)
const FUSE_FIX_ALLOW      = (false, true)

# THE INDEPENDENT EXPECTATION TABLE. Columns:
#   ood_state | conformal_status | disagree | allow_unchecked_ood | expected action | expected reason
# Derived BY HAND from D-05 and D-06; `p14_fuse` is never called to produce a value in this table.
const FUSE_TRUTH_TABLE = (
    # --- OOD FIRED: abstains regardless of everything else (D-05). 12 cells. ---------------
    (:fired,       :singleton, false, false, :abstain, :ood_fired),
    (:fired,       :singleton, false, true,  :abstain, :ood_fired),
    (:fired,       :singleton, true,  false, :abstain, :ood_fired),
    (:fired,       :singleton, true,  true,  :abstain, :ood_fired),
    (:fired,       :ambiguous, false, false, :abstain, :ood_fired),
    (:fired,       :ambiguous, false, true,  :abstain, :ood_fired),
    (:fired,       :ambiguous, true,  false, :abstain, :ood_fired),
    (:fired,       :ambiguous, true,  true,  :abstain, :ood_fired),
    (:fired,       :empty,     false, false, :abstain, :ood_fired),
    (:fired,       :empty,     false, true,  :abstain, :ood_fired),
    (:fired,       :empty,     true,  false, :abstain, :ood_fired),
    (:fired,       :empty,     true,  true,  :abstain, :ood_fired),

    # --- OOD CLEAR. The conformal set decides, and DISAGREEMENT ALONE DOES NOT SILENCE. -----
    # The four :singleton rows are D-05's amendment; the third of them is the named cell.
    (:clear,       :singleton, false, false, :decide,  :decided),
    (:clear,       :singleton, false, true,  :decide,  :decided),
    (:clear,       :singleton, true,  false, :decide,  :decided_contra_classical),
    (:clear,       :singleton, true,  true,  :decide,  :decided_contra_classical),
    (:clear,       :ambiguous, false, false, :abstain, :conformal_ambiguous),
    (:clear,       :ambiguous, false, true,  :abstain, :conformal_ambiguous),
    (:clear,       :ambiguous, true,  false, :abstain, :conformal_ambiguous),
    (:clear,       :ambiguous, true,  true,  :abstain, :conformal_ambiguous),
    (:clear,       :empty,     false, false, :abstain, :conformal_empty),
    (:clear,       :empty,     false, true,  :abstain, :conformal_empty),
    (:clear,       :empty,     true,  false, :abstain, :conformal_empty),
    (:clear,       :empty,     true,  true,  :abstain, :conformal_empty),

    # --- OOD NOT CHECKED. Abstains by DEFAULT (D-06); the opt-out is deliberate and explicit.
    # Under the opt-out the remaining triggers apply in order, and disagreement is no longer
    # "alone" -- there is an independent reason to distrust ourselves, so it abstains.
    (:not_checked, :singleton, false, false, :abstain, :ood_not_checked),
    (:not_checked, :singleton, false, true,  :decide,  :decided),
    (:not_checked, :singleton, true,  false, :abstain, :ood_not_checked),
    (:not_checked, :singleton, true,  true,  :abstain, :disagreement_and_ood),
    (:not_checked, :ambiguous, false, false, :abstain, :ood_not_checked),
    (:not_checked, :ambiguous, false, true,  :abstain, :conformal_ambiguous),
    (:not_checked, :ambiguous, true,  false, :abstain, :ood_not_checked),
    (:not_checked, :ambiguous, true,  true,  :abstain, :conformal_ambiguous),
    (:not_checked, :empty,     false, false, :abstain, :ood_not_checked),
    (:not_checked, :empty,     false, true,  :abstain, :conformal_empty),
    (:not_checked, :empty,     true,  false, :abstain, :ood_not_checked),
    (:not_checked, :empty,     true,  true,  :abstain, :conformal_empty),
)

# The two DECIDE reasons, named here so the closed-set assertion below reads as a union of two
# enumerations rather than as a magic pair.
const FUSE_FIX_DECIDE_REASONS = (:decided, :decided_contra_classical)

# --- The `ood_nulls` / `OODVerdict` fixtures for the D-06 resolver -------------------------------
# A verdict is built for BOTH raw outcomes so that the resolver's answer can be shown to depend on
# the THRESHOLD's presence and finiteness rather than on the raw outcome alone.
const FUSE_FIX_VERDICT_OVER  = OODVerdict(9.9, true,  (density = 9.9,))
const FUSE_FIX_VERDICT_UNDER = OODVerdict(0.4, false, (density = 0.4,))

# `:thr` ABSENT is the shipped default: `ood_nulls_G.jld2` persists only the frozen `:density` fit.
const FUSE_FIX_NULLS_ABSENT  = (density = :fitted,)
# Non-finite and non-Real thresholds, each of which `_recorded_ood_threshold` also treats as
# "no operating point available" (src/amortized/local_map.jl:107).
const FUSE_FIX_NULLS_INF     = (density = :fitted, thr = Inf)
const FUSE_FIX_NULLS_NEGINF  = (density = :fitted, thr = -Inf)
const FUSE_FIX_NULLS_NAN     = (density = :fitted, thr = NaN)
const FUSE_FIX_NULLS_NOTHING = (density = :fitted, thr = nothing)
const FUSE_FIX_NULLS_STRING  = (density = :fitted, thr = "3.5")
# The only shape that yields a two-valued answer.
const FUSE_FIX_NULLS_FINITE  = (density = :fitted, thr = 3.5)

const FUSE_FIX_NULLS_UNCHECKED = (FUSE_FIX_NULLS_ABSENT, FUSE_FIX_NULLS_INF, FUSE_FIX_NULLS_NEGINF,
                                  FUSE_FIX_NULLS_NAN, FUSE_FIX_NULLS_NOTHING, FUSE_FIX_NULLS_STRING)

@testset "P14 D-05 asymmetric fusion and D-06 three-valued OOD (SC2-a, SC2-b)" verbose = true begin

    @testset "the expectation table IS the complete product, once each" begin
        @test length(FUSE_TRUTH_TABLE) == 36
        keys_seen = [(r[1], r[2], r[3], r[4]) for r in FUSE_TRUTH_TABLE]
        @test length(unique(keys_seen)) == 36
        full = [(o, c, d, a) for o in FUSE_FIX_OOD_STATES, c in FUSE_FIX_CONFORMAL,
                                 d in FUSE_FIX_DISAGREE, a in FUSE_FIX_ALLOW]
        @test length(full) == 36
        @test Set(keys_seen) == Set(full)
    end

    @testset "SC2-a: the full truth table, cell by cell" begin
        for row in FUSE_TRUTH_TABLE
            ood, conf, dis, allow, want_action, want_reason = row
            @testset "($(ood), $(conf), disagree=$(dis), allow_unchecked_ood=$(allow))" begin
                r = p14_fuse(; ood_state = ood, conformal_status = conf, disagree = dis,
                               allow_unchecked_ood = allow)
                @test (r.action, r.reason) === (want_action, want_reason)
            end
        end
    end

    @testset "the ONE cell that is D-05's amendment" begin
        # Under ROADMAP SC2's literal OR this would ABSTAIN. It DECIDES, and the disagreement
        # survives onto the result as a named reason instead of being thrown away.
        r = p14_fuse(; ood_state = :clear, conformal_status = :singleton, disagree = true,
                       allow_unchecked_ood = false)
        @test r.action === :decide
        @test r.reason === :decided_contra_classical
        # Its twin without the disagreement, so the field is shown to be CARRYING the disagreement
        # rather than being the constant answer for a clear singleton.
        q = p14_fuse(; ood_state = :clear, conformal_status = :singleton, disagree = false,
                       allow_unchecked_ood = false)
        @test q.action === :decide
        @test q.reason === :decided
        @test r.reason !== q.reason
        # The opt-out is irrelevant when the OOD channel actually answered.
        @test p14_fuse(; ood_state = :clear, conformal_status = :singleton, disagree = true,
                         allow_unchecked_ood = true).reason === :decided_contra_classical
    end

    @testset "SC2-b: not-checked abstains by default (D-06)" begin
        for nulls in FUSE_FIX_NULLS_UNCHECKED, v in (FUSE_FIX_VERDICT_OVER, FUSE_FIX_VERDICT_UNDER)
            st = p14_ood_state(nulls, v)
            @test st === :not_checked
            # ...and it routes to a NAMED abstention, not a silent one.
            @test (p14_fuse(; ood_state = st, conformal_status = :singleton,
                              disagree = false).action,
                   p14_fuse(; ood_state = st, conformal_status = :singleton,
                              disagree = false).reason) === (:abstain, :ood_not_checked)
            # ...and the opt-out, set deliberately, flips it to a decision.
            @test p14_fuse(; ood_state = st, conformal_status = :singleton, disagree = false,
                             allow_unchecked_ood = true).action === :decide
        end

        # A FINITE threshold is the only shape that yields a two-valued answer, and it tracks the
        # verdict rather than the threshold's mere presence.
        @test p14_ood_state(FUSE_FIX_NULLS_FINITE, FUSE_FIX_VERDICT_OVER) === :fired
        @test p14_ood_state(FUSE_FIX_NULLS_FINITE, FUSE_FIX_VERDICT_UNDER) === :clear
        # `ood_nulls` IS NAMEDTUPLE-SHAPED, and that is not an accident of this file: the shipped
        # call site itself does `haskey(ood_nulls, :thr) ? ood_nulls.thr : nothing`
        # (src/amortized/ood.jl:394-395). Broadening the resolver to `get(ood_nulls, :thr, nothing)`
        # so that a Dict would also work would be a SILENT divergence from the shipped semantics --
        # exactly the class of divergence D-06 exists to prevent -- so a Dict carrying a threshold
        # fails LOUDLY here instead of being quietly accommodated.
        @test_throws Exception p14_ood_state(Dict(:thr => 3.5), FUSE_FIX_VERDICT_OVER)
        # A container with no threshold at all is safe on anything `haskey` answers for, and
        # resolves to the conservative state.
        @test p14_ood_state(Dict{Symbol, Any}(), FUSE_FIX_VERDICT_UNDER) === :not_checked
        # An INTEGER threshold is Real and finite, so it is a genuine operating point.
        @test p14_ood_state((thr = 3,), FUSE_FIX_VERDICT_UNDER) === :clear
    end

    @testset "a Bool never leaks out of the resolver" begin
        # The whole of D-06 is that a two-valued answer cannot represent a three-valued question.
        for nulls in (FUSE_FIX_NULLS_UNCHECKED..., FUSE_FIX_NULLS_FINITE),
            v in (FUSE_FIX_VERDICT_OVER, FUSE_FIX_VERDICT_UNDER)
            st = p14_ood_state(nulls, v)
            @test st isa Symbol
            @test !(st isa Bool)
            @test st in FUSE_FIX_OOD_STATES
        end
    end

    @testset "validation names the function, the decision and the observed value" begin
        @test_throws ArgumentError p14_fuse(; ood_state = :maybe, conformal_status = :singleton,
                                              disagree = false)
        @test_throws ArgumentError p14_fuse(; ood_state = :clear, conformal_status = :huge,
                                              disagree = false)
        # A Bool smuggled in where the three-valued state belongs is the exact D-06 regression, and
        # it is refused by the SIGNATURE rather than by the body: `ood_state::Symbol` on the keyword
        # makes it a TypeError before a single branch is evaluated, which is a stronger guarantee
        # than any runtime check in the body could give.
        @test_throws TypeError p14_fuse(; ood_state = false, conformal_status = :singleton,
                                          disagree = false)
        @test_throws TypeError p14_fuse(; ood_state = :clear, conformal_status = :singleton,
                                          disagree = :yes)
        err = try
            p14_fuse(; ood_state = :maybe, conformal_status = :singleton, disagree = false)
        catch e
            e
        end
        @test occursin("p14_fuse", err.msg)
        @test occursin("D-06", err.msg)
        @test occursin("maybe", err.msg)
        err2 = try
            p14_fuse(; ood_state = :clear, conformal_status = :huge, disagree = false)
        catch e
            e
        end
        @test occursin("p14_fuse", err2.msg)
        @test occursin("huge", err2.msg)
        # Every argument is REQUIRED. There is no default abstention and no default disagreement:
        # a caller who has not measured one of them has not asked this question.
        @test_throws UndefKeywordError p14_fuse(; conformal_status = :singleton, disagree = false)
        @test_throws UndefKeywordError p14_fuse(; ood_state = :clear, disagree = false)
        @test_throws UndefKeywordError p14_fuse(; ood_state = :clear, conformal_status = :singleton)
    end

    @testset "the reason set is CLOSED" begin
        @test p14_abstain_reasons() === (:ood_fired, :ood_not_checked, :conformal_ambiguous,
                                         :conformal_empty, :disagreement_and_ood)
        @test length(p14_abstain_reasons()) == 5
        @test length(unique(p14_abstain_reasons())) == 5
        # No reason may appear on both sides: an abstention reason that could also be a decision
        # reason would make the artifact unreadable.
        @test isempty(intersect(Set(p14_abstain_reasons()), Set(FUSE_FIX_DECIDE_REASONS)))

        # Every reason actually produced across all 36 cells is in one of the two enumerations,
        # and the action agrees with which enumeration it came from.
        for row in FUSE_TRUTH_TABLE
            r = p14_fuse(; ood_state = row[1], conformal_status = row[2], disagree = row[3],
                           allow_unchecked_ood = row[4])
            @test r.action in (:decide, :abstain)
            if r.action === :abstain
                @test r.reason in p14_abstain_reasons()
            else
                @test r.reason in FUSE_FIX_DECIDE_REASONS
            end
        end
        # And every enumerated abstention reason is REACHABLE -- a closed set with a dead member is
        # a set that stopped describing the rule.
        produced = Set(p14_fuse(; ood_state = row[1], conformal_status = row[2],
                                  disagree = row[3], allow_unchecked_ood = row[4]).reason
                       for row in FUSE_TRUTH_TABLE)
        @test Set(p14_abstain_reasons()) ⊆ produced
        @test Set(FUSE_FIX_DECIDE_REASONS) ⊆ produced
    end

    @testset "the status vocabulary is single-sourced with the conformal hedge" begin
        # spike/p14/conformal.jl declares the same three symbols under the same name behind the
        # same guard, and whichever file loads second ASSERTS agreement rather than silently
        # losing. This pins the copy that lives here.
        @test P14_CONFORMAL_STATUSES === (:singleton, :ambiguous, :empty)
        @test Set(P14_CONFORMAL_STATUSES) == Set(FUSE_FIX_CONFORMAL)
    end

    @testset "P14 fusion ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
