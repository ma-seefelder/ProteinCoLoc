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

# spike/test/test_p14_result.jl --- the Phase-14 result-type contract.
#
# WHAT THIS FILE IS FOR. spike/p14/result.jl is where this phase's honesty commitments stop being
# prose and become fields. A field can be renamed, a docstring can be deleted and an invariant can
# be relaxed by whoever is in a hurry; each one below is an assertion that turns red when that
# happens.
#
#   1. THE TYPE IS A REAL SUBTYPE, AND NOTHING WAS BOLTED ONTO A SHIPPED STRUCT. Phase-7 D-02
#      requires every new result variant to slot in as a new `AbstractColocResult` subtype.
#      Defining the executable type in spike/ (D-01) must not weaken that, so the subtype relation
#      is asserted against the ACTUAL supertype from src/results.jl and the shipped
#      `AmortizedColocResult` field list is asserted UNCHANGED.
#   2. THE TWO LOSSY ACCESSORS ARE LOSSY IN THE DOCUMENTED DIRECTION. `is_ood` collapses `:clear`
#      and `:not_checked` -- which are DIFFERENT STATES (D-06) -- and the test asserts both that it
#      collapses them and that `ood_state` still distinguishes them. `delta_rho` must fall through
#      to the interface error rather than fabricate a contrast no control posterior backs.
#   3. THE PER-TILE MAP CANNOT BE READ AS FDR-CONTROLLED (D-04). Asserted, not merely documented:
#      the FIELD NAME carries the label, `p14_local_map_is_controlled` is executable, and the
#      shipped `LocalColocMap` is asserted to carry NO uncertainty field -- which is the evidence
#      D-04 rests on rather than an opinion about it.
#   4. THE BATCH OBJECT CANNOT OMIT ITS DENOMINATOR. `(alpha_FDR, n_decided/n_total)` is one
#      quantity; a rule that abstains on 95 % of a batch hits any alpha trivially.
#
# WHY THE `fdr` INPUT IS A HAND-BUILT FIXTURE RATHER THAN A CALL INTO spike/p14/fdr.jl. Loading
# fdr.jl pulls posterior.jl, which pulls the Phase-13 net and with it Flux and NeuralEstimators --
# roughly sixteen seconds of package load -- to produce a NamedTuple whose SHAPE is the only thing
# this file is testing. The arithmetic inside it belongs to spike/test/test_p14_fdr.jl and is
# tested there. The fixture is kept honest by a SOURCE GREP against fdr.jl's own returned key
# names (testset 7), so a rename there fails here instead of drifting silently.
#
# Mirrors spike/test/test_p13_result.jl: license header -> using Test -> guarded include of the
# unit under test -> fixtures computed ONCE at module level -> one outer @testset. No training, no
# simulation, no figure call, and no random number is drawn at all.
#
# RUN IT PER FILE. The aggregate spike suite exits 1 early at the Phase-4 speedup gate and masks
# every later include block, so a green aggregate run would be no evidence at all:
#
#     julia --project=spike spike/test/test_p14_result.jl

using Test

# Unit under test. Pulls src/results.jl, src/registry.jl and src/amortized/local_map.jl READ-ONLY,
# plus spike/p13/result.jl and spike/p14/{consts,fuse}.jl.
isdefined(@__MODULE__, :P14Result) || include(joinpath(@__DIR__, "..", "p14", "result.jl"))

# GUARDED, NOT REDEFINED. `P14_REPO_ROOT` is already a `const` in spike/p14/provenance.jl:83 and in
# spike/test/test_p14_decoupling.jl:66, and the spike lane is a FLAT top-level namespace: when
# several files run in one process the definitions land in the same module. The values are
# identical (all normalise to the repo root), so this file reuses an existing binding when one is
# present and defines it only when running standalone.
if !isdefined(@__MODULE__, :P14_REPO_ROOT)
    const P14_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
end

# --- Fixtures built ONCE at module level ---------------------------------------------------------
# A perfectly calibrated toy reliability input, only so the composed Phase-13 result has a real
# CalibrationMeta to carry. Nothing here gates on its numbers.
const RES14_CAL  = _bin_calibration(fill(0.5, 100), vcat(trues(50), falses(50));
                                    n_bins = P13_ECE_NBINS)
const RES14_META = p13_calibration_meta(RES14_CAL; grid = 8, auc = 0.97)
const RES14_OOD  = OODVerdict(1.23, false, (; density = 1.23))
const RES14_TW   = ThreeHypothesisColocResult(8, zeros(7, 5),
                                              (coloc = 2.5, random = 0.0, exclusion = -3.25),
                                              RES14_OOD, RES14_META, NamedTuple())

# The DELIBERATELY UNEQUAL class posterior: 0.7 / 0.2 / 0.1. Equal masses would let any of the five
# class orderings that coexist in this codebase pass testset 8 by coincidence.
const RES14_P = (coloc = 0.7, random = 0.2, exclusion = 0.1)

# The per-tile map fixture. Its field list is asserted in testset 5 -- it is the evidence for D-04,
# not a convenience.
const RES14_MAP = LocalColocMap(8, (2, 2), zeros(2, 2), falses(2, 2), (; n_sentinel = 0))

"Build a valid P14Result; every field is overridable so a single invariant can be broken at a time."
_p14res(; three_way = RES14_TW, class_posterior = RES14_P,
          null_posterior = 1.0 - RES14_P.coloc,
          null_split = (p_random = RES14_P.random, p_exclusion = RES14_P.exclusion,
                        v = 1.0 - RES14_P.coloc),
          decision = :coloc, abstain_reason = nothing, ood_state = :clear,
          conformal = (set = (:coloc,), status = :singleton, qhat = 0.4),
          cross_method = (basis = :ghat_rho_true_scale, disagree_any = false),
          local_map_uncontrolled = nothing, meta = NamedTuple()) =
    P14Result(three_way, class_posterior, null_posterior, null_split, decision, abstain_reason,
              ood_state, conformal, cross_method, local_map_uncontrolled, meta)

# --- Batch fixtures -------------------------------------------------------------------------------
# The rule output, HAND-BUILT (see the banner) and kept honest by the source grep in testset 7.
const RES14_FDR = (accepted = [1], k_star = 1, fdp = 0.02, t_star = 0.02, cost_ratio = 0.0204,
                   n = 2, n_total = 4, n_decided = 2, decided_fraction = 0.5,
                   fdr_scope = :decided_subset_only)
const RES14_SENS = [(pi_coloc = pc, k_star = 1, fdp = 0.02, n_accepted = 1)
                    for pc in P14_PI_COLOC_GRID]
const RES14_PRIOR = (coloc = 0.217, random = 0.471, exclusion = 0.312)
const RES14_BATCH_RESULTS = [_p14res(decision = :coloc),
                             _p14res(decision = :not_coloc),
                             _p14res(decision = :abstain, abstain_reason = :ood_not_checked,
                                     ood_state = :not_checked),
                             _p14res(decision = :abstain, abstain_reason = :conformal_empty,
                                     conformal = (set = (), status = :empty, qhat = 0.4))]
const RES14_COUNTS = (ood_not_checked = 1, conformal_empty = 1)

"Build a valid P14BatchDecision; `alpha_fdr` and the ledger are the parts worth varying."
_p14batch(; results = RES14_BATCH_RESULTS, alpha = 0.10, fdr = RES14_FDR, n_total = 4,
            counts = RES14_COUNTS, sens = RES14_SENS) =
    p14_batch_decision(results;
                       alpha_fdr = alpha,
                       fdr = fdr,
                       n_total = n_total,
                       abstain_reason_counts = counts,
                       class_prior_used = RES14_PRIOR,
                       class_prior_source = :measured_pi_class_masses,
                       pi_sensitivity = sens,
                       conformal = (qhat = 0.4, vacuous_hedge = false),
                       meta = (; note = "fixture"))

# The executable source of the two files this test greps, comment-stripped BEFORE any assertion so
# a grep cannot be satisfied (or defeated) by prose.
_p14_code(path) = join(filter(l -> !startswith(strip(l), "#"),
                              split(read(path, String), '\n')), '\n')
const RES14_FDR_CODE = _p14_code(joinpath(P14_REPO_ROOT, "spike", "p14", "fdr.jl"))

@testset "P14 result and batch types (D-01, D-04, D-05, D-06)" verbose = true begin

    @testset "is a new AbstractColocResult subtype (Phase-7 D-02, D-01)" begin
        # The supertype is the REAL one from src/results.jl, reached read-only -- not a spike-local
        # look-alike. If that read-only include ever broke, this assertion is what would notice.
        @test P14Result <: AbstractColocResult
        @test isabstracttype(AbstractColocResult)
        @test isconcretetype(P14Result)
        @test _p14res() isa AbstractColocResult
        # D-02's other half: the new variant does NOT bolt fields onto the shipped struct.
        @test fieldnames(AmortizedColocResult) == (:grid, :posterior, :delta_rho_draws,
                                                   :log_bayes_factor, :ood, :calibration, :meta)
        # COMPOSITION, NOT RESTATEMENT (14-RESEARCH D.2 (a)): the Phase-13 result is carried as a
        # FIELD, so D-08's structural zero is enforced once, at Phase 13's own constructor.
        @test :three_way in fieldnames(P14Result)
        @test fieldtype(P14Result, :three_way) === ThreeHypothesisColocResult
        @test three_way(_p14res()) === RES14_TW
        @test log_bf_vs_random(_p14res()).random === 0.0
    end

    @testset "the four required accessors" begin
        r = _p14res()
        @test posterior_draws(r) == zeros(7, 5)
        @test bayes_factor(r) == 2.5
        @test bayes_factor(r) == log_bf_vs_random(r).coloc
        @test is_ood(r) === false
        # DO NOT FABRICATE A DELTA RHO. No control posterior was carried, so there is nothing to
        # return; the documented interface error is the honest answer.
        @test_throws ErrorException delta_rho(r)
        err = try; delta_rho(r); catch e; e; end
        @test occursin("delta_rho", err.msg)
        @test occursin("P14Result", err.msg)
        @test occursin("accessor interface", err.msg)
        # When the run DID carry control draws they are returned unchanged -- the guard is about
        # absence, not a blanket refusal.
        draws = [0.11, -0.04, 0.37]
        @test delta_rho(_p14res(meta = (; delta_rho_draws = draws))) == draws
    end

    @testset "is_ood is lossy, and lossy in the documented direction (D-06)" begin
        # THE COLLAPSE, ASSERTED IN BOTH DIRECTIONS. `:clear` means a detector answered;
        # `:not_checked` means none ran, because `ood_verdict` returns a false verdict for EVERY
        # input when no operating point was recorded. Reading the second as the first makes the
        # tool speak confidently in exactly the regime where it has no evidence it should.
        @test is_ood(_p14res(ood_state = :clear)) === false
        @test is_ood(_p14res(ood_state = :not_checked,
                             decision = :abstain,
                             abstain_reason = :ood_not_checked)) === false
        @test is_ood(_p14res(ood_state = :fired,
                             decision = :abstain, abstain_reason = :ood_fired)) === true
        # ... and the accessor that does NOT lose the distinction.
        @test ood_state(_p14res(ood_state = :clear)) === :clear
        @test ood_state(_p14res(ood_state = :not_checked, decision = :abstain,
                                abstain_reason = :ood_not_checked)) === :not_checked
        @test ood_state(_p14res(ood_state = :fired, decision = :abstain,
                                abstain_reason = :ood_fired)) === :fired
        # A fourth state is not storable at all...
        @test_throws ArgumentError _p14res(ood_state = :unchecked)
        @test_throws ArgumentError _p14res(ood_state = :in_distribution)
        # ... and neither is a Bool. THAT IS THE WHOLE OF D-06: the OOD input is three-valued, so a
        # `false` cannot be stored here and then be read back as "in distribution".
        @test_throws MethodError _p14res(ood_state = false)
        # THE LOSSINESS MUST BE RETRIEVABLE FROM THE RUNNING SYSTEM, not just readable in the
        # source: Julia 1.12 drops every docstring written inside an `if ... end` block, so this is
        # what proves the guard-block split was copied correctly.
        doc = string(@doc(is_ood))
        @test occursin("not_checked", doc)
        @test occursin("LOSSY", doc)
        # THE SOURCE-GREP TWIN. The decision layer must never route through the Bool accessor. It
        # arms itself when decide.jl lands.
        _decide = joinpath(P14_REPO_ROOT, "spike", "p14", "decide.jl")
        if isfile(_decide)
            @test !occursin("is_ood(", _p14_code(_decide))
        else
            @info "spike/p14/decide.jl does not exist yet — the never-call-is_ood grep arms " *
                  "itself when it lands"
        end
    end

    @testset "abstention is explainable" begin
        # Every ABSTAIN carries which trigger fired: a silent abstention is nearly as bad as a
        # wrong call, because a reader cannot tell a tool that DECLINED from one that BROKE.
        for reason in p14_abstain_reasons()
            r = _p14res(decision = :abstain, abstain_reason = reason)
            @test decision(r) === :abstain
            @test abstain_reason(r) === reason
            @test abstain_reason(r) in p14_abstain_reasons()
        end
        # A silence with no recorded trigger is not constructible.
        @test_throws ArgumentError _p14res(decision = :abstain, abstain_reason = nothing)
        # ... and neither is a reason recorded on an item that was never silenced.
        @test_throws ArgumentError _p14res(decision = :coloc, abstain_reason = :ood_fired)
        @test_throws ArgumentError _p14res(decision = :not_coloc, abstain_reason = :ood_fired)
        # A reason from outside the closed set cannot appear untracked.
        @test_throws ArgumentError _p14res(decision = :abstain, abstain_reason = :because)
        # And the decision vocabulary is closed too.
        @test_throws ArgumentError _p14res(decision = :maybe)
        @test abstain_reason(_p14res()) === nothing
        @test Set(P14_DECISIONS) == Set((:coloc, :not_coloc, :abstain))
    end

    @testset "the per-tile map is labelled uncontrolled (D-04)" begin
        r = _p14res(local_map_uncontrolled = RES14_MAP)
        # THE LABEL IS IN THE TYPE, not only in a docstring: a caller reading a field list cannot
        # mistake this for a controlled call.
        @test :local_map_uncontrolled in fieldnames(P14Result)
        @test local_map(r) === RES14_MAP
        @test local_map(_p14res()) === nothing
        # EXECUTABLE, GREPPABLE, AND FALSE FOR EVERY RESULT -- including one carrying a map.
        @test p14_local_map_is_controlled(r) === false
        @test p14_local_map_is_controlled(_p14res()) === false
        # THE EVIDENCE D-04 RESTS ON, asserted rather than asserted-about: the shipped map carries
        # a per-tile Delta rho and a per-tile flag and NO uncertainty field, so there is nothing
        # per-tile to control against. Phase 12 returned NO on calibrated per-region uncertainty.
        @test fieldnames(LocalColocMap) == (:grid, :tiles, :delta_rho, :ood_flag, :meta)
        @test !(:uncertainty in fieldnames(LocalColocMap))
        @test !(:region_sd in fieldnames(LocalColocMap))
        # The docstring must keep saying so in the words D-04 uses.
        doc = string(@doc(local_map))
        @test occursin("UNCONTROLLED DISPLAY", doc)
        @test occursin("D-04", doc)
    end

    @testset "cross-method disagreement is recorded, not suppressed (D-05)" begin
        # THE CELL THAT IS THE WHOLE AMENDMENT: a disagreeing item that still DECIDES. Under SC2 as
        # originally written this would have been silenced -- and silenced exactly where Phase 9
        # says the tool should speak, since disagreeing with Costes/Manders is the product thesis.
        r = _p14res(decision = :coloc,
                    cross_method = (basis = :ghat_rho_true_scale, classical_call = :random,
                                    disagree_label = true, disagree_costes = false,
                                    disagree_any = true))
        @test decision(r) === :coloc
        @test abstain_reason(r) === nothing
        @test cross_method(r).disagree_any === true
        # `basis` records WHICH SCALE the comparison was made on -- the thing Pitfall 5 says is
        # otherwise forgotten -- so a record without it is not constructible.
        @test cross_method(r).basis === :ghat_rho_true_scale
        @test_throws ArgumentError _p14res(cross_method = (disagree_any = true,))
        @test_throws ArgumentError _p14res(cross_method = NamedTuple())
    end

    @testset "the batch object cannot omit its denominator" begin
        b = _p14batch()
        # THE SCOPE, MACHINE-READABLE. An artifact carrying this cannot be read as claiming
        # batch-wide control.
        @test b.fdr_scope === :decided_subset_only
        @test b.n_decided + b.n_abstained == b.n_total
        @test b.decided_fraction == b.n_decided / b.n_total
        @test length(b.pi_sensitivity) == length(P14_PI_COLOC_GRID)
        @test Tuple(e.pi_coloc for e in b.pi_sensitivity) == P14_PI_COLOC_GRID
        # THE DENOMINATOR IS A REQUIRED INPUT. `(alpha_FDR, n_decided/n_total)` is ONE quantity: a
        # rule that abstains on 95 % of a batch hits any alpha trivially.
        @test_throws UndefKeywordError p14_batch_decision(
            RES14_BATCH_RESULTS;
            alpha_fdr = 0.10,
            fdr = RES14_FDR,
            abstain_reason_counts = RES14_COUNTS,
            class_prior_used = RES14_PRIOR,
            class_prior_source = :measured_pi_class_masses,
            pi_sensitivity = RES14_SENS,
            conformal = (qhat = 0.4,))
        # ... and it is checked, not merely demanded: a misstated denominator is rejected.
        @test_throws ArgumentError _p14batch(n_total = 5)
        # The abstention ledger must balance -- a batch whose reasons do not add up has silenced
        # items nobody can explain.
        @test_throws ArgumentError _p14batch(counts = (ood_not_checked = 1,))
        @test_throws ArgumentError _p14batch(counts = (ood_not_checked = 1, because = 1))
        # The sensitivity curve is REQUIRED output, never optional (SC1-c: reported, not gated).
        @test_throws ArgumentError _p14batch(sens = RES14_SENS[1:2])
        @test_throws ArgumentError _p14batch(
            sens = [(pi_coloc = pc, k_star = 1) for pc in P14_PI_COLOC_GRID])
        # No other scope is storable, by any entry point.
        @test_throws ArgumentError P14BatchDecision(
            RES14_BATCH_RESULTS, 0.10, 1, 0.02, 0.0204, 0.02, [1], 4, 2, 2, 0.5, RES14_COUNTS,
            RES14_PRIOR, :measured_pi_class_masses, RES14_SENS, NamedTuple(), :whole_batch,
            NamedTuple())
        # THE HEADLINE CARRIES BOTH NUMBERS IN ONE STRING, so a copy-paste cannot separate them.
        h = p14_headline(b)
        @test occursin("0.1", h)
        @test occursin("2/4", h)
        @test occursin("FDR SCOPE: DECIDED SUBSET ONLY", h)
        @test occursin("FDR SCOPE: DECIDED SUBSET ONLY", p14_headline_template_probe())
        # THE HAND-BUILT `fdr` FIXTURE IS KEPT HONEST BY A SOURCE GREP against fdr.jl's own
        # returned key names, so a rename there fails HERE instead of drifting silently.
        for k in keys(RES14_FDR)
            @test occursin(Regex("\\b$(k)\\s*="), RES14_FDR_CODE)
        end
        @test occursin(r"fdr_scope\s*=\s*:decided_subset_only", RES14_FDR_CODE)
        # The named limits are enumerable FROM CODE, so a report that drops one fails to iterate a
        # tuple rather than forgetting a sentence.
        @test length(p14_named_limits()) >= 8
        @test all(s -> s isa String && !isempty(s), p14_named_limits())
        @test any(s -> occursin("NOT shipped in v2.0", s), p14_named_limits())
        @test any(s -> occursin("DECIDED SUBSET ONLY", s), p14_named_limits())
        @test any(s -> occursin("UNCONTROLLED DISPLAY", s), p14_named_limits())
    end

    @testset "class order survives the round trip" begin
        # FIVE class orderings coexist in this codebase and every one of ECE, AUC and FDR is
        # invariant to a CONSISTENT relabelling -- which is exactly what makes an inconsistent one
        # so dangerous: it produces numbers that look right. The masses are deliberately unequal,
        # and each is asserted separately so a failure NAMES the class that moved.
        r = _p14res()
        @test class_posterior(r).coloc     == 0.7
        @test class_posterior(r).random    == 0.2
        @test class_posterior(r).exclusion == 0.1
        # A posterior declared in a different order lands on the same names.
        r2 = _p14res(class_posterior = (exclusion = 0.1, random = 0.2, coloc = 0.7))
        @test class_posterior(r2).coloc     == 0.7
        @test class_posterior(r2).exclusion == 0.1
        # The composite null is the exact complement of coloc, and a stored value that disagrees
        # with its own class posterior is rejected.
        @test null_posterior(r) == 1.0 - 0.7
        @test null_split(r).p_random    == 0.2
        @test null_split(r).p_exclusion == 0.1
        @test_throws ArgumentError _p14res(null_posterior = 0.5)
        # A posterior that does not sum to 1 makes every null posterior downstream of it wrong.
        @test_throws ArgumentError _p14res(
            class_posterior = (coloc = 0.7, random = 0.2, exclusion = 0.2), null_posterior = 0.3)
        @test_throws ArgumentError _p14res(class_posterior = (coloc = 0.7, random = 0.3))
        # The hedge vocabulary is the one spike/p14/conformal.jl and fuse.jl share.
        @test_throws ArgumentError _p14res(
            conformal = (set = (:coloc,), status = :one, qhat = 0.4))
        @test_throws ArgumentError _p14res(conformal = (status = :singleton,))
        @test conformal_set(r) == (:coloc,)
        @test conformal_status(r) === :singleton
    end

    @testset "src/ is byte-unchanged (D-01)" begin
        # D-01 says the decision layer is built with a `src/`-shaped signature. "src-shaped" is NOT
        # "in src": this type may subtype the shipped supertype and may still not be written there.
        # The decision layer is NOT shipped in v2.0, and this assertion is what makes "provably
        # untouched" mean something in this file as well as in the lane guard.
        @test P14_SRC_UNTOUCHED == true
        _has_git = Sys.which("git") !== nothing &&
                   (isdir(joinpath(P14_REPO_ROOT, ".git")) ||
                    isfile(joinpath(P14_REPO_ROOT, ".git")))
        if _has_git
            @test success(Cmd(`git diff --quiet HEAD -- src`; dir = P14_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/results.jl`; dir = P14_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/amortized/local_map.jl`;
                              dir = P14_REPO_ROOT))
            @test isempty(readchomp(Cmd(`git status --porcelain -- src`; dir = P14_REPO_ROOT)))
        else
            @info "git unavailable — skipping the executable src/ byte-equality assertions"
        end
    end

    @testset "the result types ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
