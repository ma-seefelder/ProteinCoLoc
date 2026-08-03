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

# spike/test/test_p14_decide.jl --- the composed decision pipeline, asserted end to end at
# FIXTURE SCALE.
#
# THE FOUR CLAIMS THIS FILE BINDS, in the order they would fail silently:
#
#   1. THE ORDERING. ABSTAIN FIRST, THEN SORT. Testset 4 does not merely check that no abstained
#      item appears in the accepted set -- that is true under BOTH orders, which is exactly why it
#      is not enough on its own. It recomputes the rule BOTH WAYS in the test, SEARCHES for a
#      level at which the two orders give different answers, asserts such a level exists, and then
#      asserts `decide_coloc` follows the abstain-first one. A scratch reorder of the
#      implementation therefore turns this file red rather than leaving it green.
#   2. THE THRESHOLD IS LOADED, NEVER TYPED. Testset 1 asserts the inherited cut arrives with its
#      four-way provenance agreement AND that there is no keyword by which a caller could
#      substitute one (D-07).
#   3. NOT-CHECKED IS NOT "IN DISTRIBUTION". Testset 2 asserts the default abstains and names its
#      trigger, and that the opt-out is something a caller must set deliberately (D-06).
#   4. THE SCALE AND THE SEED. Testset 9 runs the src-shaped entry on a synthetic pair and asserts
#      the classical comparison was made on the rescaled scale and consumed the Phase-14 dev
#      stream -- on a LIVE result, not only in a source grep (SC2-d, Pitfall 7).
#
# NO SIMULATION AND NO IMAGE LOADING. Every summary is a synthetic draw from the FIXTURE stream at
# the fixture counter -- never a reported counter -- and the two synthetic images are random
# matrices built through the lane's own container. Loading the net and the frozen standardizer IS
# permitted: they are the objects under test, a forward pass is milliseconds, and a decision layer
# tested without its net would be testing nothing.
#
# THE FIXTURE IS DETERMINISTIC WITHOUT PINNED MAGIC NUMBERS. The Philox stream is counter-based, so
# the same draws come back on every machine; every threshold the fixture needs (the operating point
# that splits the batch, the level at which the two orders disagree) is DERIVED FROM THE DRAWS
# INSIDE THE TEST rather than transcribed. A hard-coded number here would go stale silently the
# first time anything upstream moved.
#
# THE OPERATING POINT BELOW IS A FIXTURE, NOT AN OPERATING POINT. It is chosen so that a known
# half of the batch fires; nothing about it is pre-registered and no number computed here is
# reportable.
#
# RUN IT PER FILE. The aggregate spike suite exits 1 early at the Phase-4 speedup gate and masks
# every later include block, so a green aggregate run would be no evidence at all:
#
#     julia --project=spike spike/test/test_p14_decide.jl

using Test
using Random

# Unit under test. Pulls the whole Phase-14 lane plus the Phase-13 read surface, the classical
# comparator and -- READ-ONLY -- the shipped result types and OOD channel.
isdefined(@__MODULE__, :decide_coloc) || include(joinpath(@__DIR__, "..", "p14", "decide.jl"))

# GUARDED, NOT REDEFINED. The spike lane is a FLAT top-level namespace: when several files run in
# one process every definition lands in the same module, so re-`const`-ing a name is a real
# collision. The value is identical wherever it is defined.
if !isdefined(@__MODULE__, :P14_REPO_ROOT)
    const P14_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
end

# --- Fixtures built ONCE at module level ----------------------------------------------------------

const DEC_B     = p14_load_bundle()
const DEC_N     = p14_summary_length(DEC_B)          # DERIVED from the net, never a literal
const DEC_BASIS = load_p13_basis()
const DEC_RNG   = p14_fix_rng(P14_FIXTURE_COUNTER)   # the FIXTURE stream, never a reported one

# The density null the OOD channel scores against, fitted on synthetic draws. Three hundred columns
# so the covariance over the continuous rows is not rank-deficient; the ridge in `fit_ood_nulls`
# would hide that rather than fix it.
const DEC_POOL  = reduce(hcat, [randn(DEC_RNG, DEC_N) for _ in 1:300])
const DEC_NULLS = fit_ood_nulls(DEC_POOL; variant = DEC_BASIS.variant)

# Twelve synthetic pairs. Small enough to run in seconds, large enough that the prefix rule has
# something to cut.
const DEC_ITEMS = [(Zs = randn(DEC_RNG, DEC_N), Zc = randn(DEC_RNG, DEC_N)) for _ in 1:12]

# The scores the shipped channel would compute, and a fixture split that fires the upper half.
# DERIVED from the draws: the sixth-largest score is the cut, so exactly six items fire.
const DEC_SCORES = [max(maha_score(DEC_NULLS, it.Zs), maha_score(DEC_NULLS, it.Zc))
                    for it in DEC_ITEMS]
const DEC_THR    = sort(DEC_SCORES)[6]

# THE TWO OOD INPUT SHAPES THE PHASE DISTINGUISHES (D-06). With the operating point, the channel
# answers; without it, `ood_verdict` returns a false flag for EVERY input unconditionally, which
# means NOT CHECKED and must never be read as in-distribution.
const DEC_OOD           = (density = DEC_NULLS, thr = DEC_THR)
const DEC_OOD_UNCHECKED = (density = DEC_NULLS,)

# The conformal threshold. A FIXTURE value chosen so most sets come out singletons and the OOD
# channel is what drives the partition; it is not the pre-registered miscoverage level and no
# number here is a coverage claim.
const DEC_QHAT = 0.55

# Two synthetic images, for the src-shaped entry point and its validation. Random matrices through
# the lane's own container -- no file is opened.
_dec_img(seed_matrices) = MultiChannelImage(seed_matrices, ["c1", "c2", "c3"], "fixture",
                                            ["", "", ""], size(seed_matrices[1]),
                                            [0.5, 0.5, 0.5])
const DEC_IMG_S = _dec_img([abs.(randn(DEC_RNG, 64, 64)) for _ in 1:3])
const DEC_IMG_C = _dec_img([abs.(randn(DEC_RNG, 64, 64)) for _ in 1:3])

"Run the batch rule at one level, with the fixture operating point and the deliberate opt-out."
_dec_batch(alpha; ood = DEC_OOD, allow = true, items = DEC_ITEMS) =
    decide_coloc(DEC_B, items, ood, DEC_QHAT; alpha_fdr = alpha, allow_unchecked_ood = allow)

# ONE batch at a level chosen below, reused by several testsets so the net is not re-run per
# assertion.
const DEC_BATCH = _dec_batch(0.70)

@testset "P14 decide: the composed pipeline, abstain-first, at fixture scale (14-07)" verbose = true begin

    @testset "the bundle LOADS the inherited cut and cannot be handed another (D-07)" begin
        @test DEC_B.tau === 0.15
        @test DEC_B.prov.four_way_agreed
        @test DEC_B.prov.source === :three_way_net_jld2

        # THERE IS NO KEYWORD BY WHICH A CALLER COULD SUBSTITUTE THE CUT. Copying the value is how
        # two numbers silently diverge and re-deriving it guarantees a second different cut on the
        # same hypothesis space, so the loader offers neither.
        kws = Base.kwarg_decl(first(methods(p14_load_bundle)))
        @test !(:tau in kws)
        @test Set(kws) == Set([:net_path, :probe_path, :gate_path])

        # The class prior is the MEASURED class mix of the prior the net was trained under, on the
        # frozen class keys and summing to 1. Measured rather than chosen, which is why it is the
        # default.
        @test Set(keys(DEC_B.prior)) == Set(P14_CLASS_KEYS)
        @test isapprox(DEC_B.prior.coloc + DEC_B.prior.random + DEC_B.prior.exclusion, 1.0;
                       atol = 1e-9)

        # PHASE 13's CALIBRATION IS CARRIED, NOT RECOMPUTED. Re-running the binning on Phase-14
        # draws would produce a second, differently-sourced ECE that a reader would conflate with
        # the gated one. Asserted against the gate report the number was banked in.
        rep = _p14_load_report(DEC_B.prov.gate_path, "three-way gate")
        @test DEC_B.calibration.ece == rep["ece_coloc"]
        @test DEC_B.calibration.mce == rep["mce_coloc"]
        @test DEC_B.calibration.gate.auc == rep["head_auc_coloc"]
        @test DEC_B.calibration.gate.carried_forward === true
        @test DEC_B.calibration.gate.recomputed_in_phase14 === false
        # The summary width is DERIVED from the net, so a grid change moves it rather than
        # silently mismatching a literal.
        @test DEC_N == 2 * Int(DEC_B.h.num_summaries)
    end

    @testset "a NOT-CHECKED OOD channel abstains by default (D-06)" begin
        it = DEC_ITEMS[1]
        r = p14_decide_one(DEC_B, it.Zs, it.Zc, DEC_OOD_UNCHECKED, DEC_QHAT;
                           lambda = P13_PHASE11_REFERENCE_LAMBDA, idx = 1)
        @test ood_state(r) === :not_checked
        @test decision(r) === :abstain
        @test abstain_reason(r) === :ood_not_checked
        # EVERY ABSTAIN CARRIES ITS TRIGGER, and the trigger is from the closed set.
        @test abstain_reason(r) in p14_abstain_reasons()

        # The opt-out is something a caller must set DELIBERATELY, and it opts out of the
        # not-checked branch only -- every other trigger still applies.
        r2 = p14_decide_one(DEC_B, it.Zs, it.Zc, DEC_OOD_UNCHECKED, DEC_QHAT;
                            lambda = P13_PHASE11_REFERENCE_LAMBDA, idx = 1,
                            allow_unchecked_ood = true)
        @test ood_state(r2) === :not_checked
        @test conformal_status(r2) === :singleton      # nothing else is firing on this item
        @test cross_method(r2).disagree_any === false
        @test decision(r2) !== :abstain
        @test abstain_reason(r2) === nothing

        # With the operating point present the channel actually answers, and the state is one of
        # the two ANSWERS rather than the third non-answer.
        r3 = p14_decide_one(DEC_B, it.Zs, it.Zc, DEC_OOD, DEC_QHAT;
                            lambda = P13_PHASE11_REFERENCE_LAMBDA, idx = 1)
        @test ood_state(r3) in (:fired, :clear)

        # NO IMAGES SUPPLIED IS NOT AGREEMENT WITH THE CLASSICS. A bare `false` here would read as
        # "the classics were consulted and concurred".
        @test cross_method(r2).basis === :not_computed
        @test haskey(cross_method(r2), :note)
    end

    @testset "abstained items never enter the accepted set" begin
        b = DEC_BATCH
        @test b.n_decided + b.n_abstained == b.n_total
        @test b.n_total == length(DEC_ITEMS)
        @test 0 < b.n_decided < b.n_total                    # the fixture really is mixed
        @test sum(values(b.abstain_reason_counts)) == b.n_abstained
        @test all(k -> k in p14_abstain_reasons(), keys(b.abstain_reason_counts))

        # `accepted` holds positions INTO THE DECIDED SUBSET; the batch-index mapping is carried in
        # meta, and it is the mapping an off-by-one would corrupt SILENTLY.
        acc = b.meta.accepted_batch_indices
        @test length(acc) == b.k_star
        @test all(i -> b.results[i].decision !== :abstain, acc)
        @test all(i -> b.results[i].decision === :coloc, acc)
        @test all(i -> 1 <= i <= b.n_decided, b.accepted)
        @test b.fdr_scope === :decided_subset_only
        @test b.decided_fraction == b.n_decided / b.n_total

        # A decided item that the rule did NOT accept is :not_coloc, and its provisional call is
        # preserved so a divergence between the argmax and the rule stays visible.
        for i in b.meta.decided_batch_indices
            @test b.results[i].decision === (i in acc ? :coloc : :not_coloc)
            @test haskey(b.results[i].meta, :provisional_decision)
        end
    end

    @testset "ABSTAIN FIRST, THEN SORT -- and the two orders really do differ here" begin
        b = DEC_BATCH
        v_all     = [null_posterior(r) for r in b.results]
        dec_idx   = b.meta.decided_batch_indices
        v_decided = v_all[dec_idx]

        # (a) THE RULE WAS RUN OVER THE DECIDED SUBSET ONLY, recomputed here from the results
        #     rather than read back off the object under test.
        indep = p14_bayes_fdr(v_decided, b.alpha_fdr)
        @test indep.k_star == b.k_star
        @test indep.fdp == b.posterior_expected_fdp
        @test sort([dec_idx[j] for j in indep.accepted]) == sort(b.meta.accepted_batch_indices)

        # (b) THE REALIZED FDP IS THE MEAN OF THE ACCEPTED NULL POSTERIORS, recomputed from
        #     scratch. Zero accepted is a legal answer and is not a mean.
        if b.k_star > 0
            @test isapprox(b.posterior_expected_fdp,
                           sum(v_all[b.meta.accepted_batch_indices]) / b.k_star; atol = 1e-12)
            @test b.posterior_expected_fdp <= b.alpha_fdr + 1e-12
        end

        # (c) THE DISCRIMINATING HALF. Sorting the whole batch first and abstaining afterwards
        #     also yields an accepted set containing no abstained item -- so "no abstained item was
        #     accepted" is true under BOTH orders and proves nothing on its own. What separates
        #     them is WHERE THE CUT FALLS: the whole-batch running mean is computed over items the
        #     decided-subset mean never sees, so the two stop at different k. This searches for a
        #     level at which they disagree, asserts one exists, and then asserts the implementation
        #     follows the abstain-first answer at every such level.
        abst = [r.decision === :abstain for r in b.results]
        disagreeing = Tuple{Float64, Int, Int}[]
        for a in 0.01:0.01:0.99
            k_abstain_first  = p14_bayes_fdr(v_decided, a).k_star
            k_sort_then_drop = count(i -> !abst[i], p14_bayes_fdr(v_all, a).accepted)
            k_abstain_first == k_sort_then_drop || push!(disagreeing, (a, k_abstain_first,
                                                                       k_sort_then_drop))
        end
        @test !isempty(disagreeing)
        for (a, k_abstain_first, _) in disagreeing
            @test _dec_batch(a).k_star == k_abstain_first
        end

        # (d) An abstained item is neither a discovery nor a non-discovery: it appears in NEITHER
        #     the numerator nor the denominator of the quoted quantity.
        @test length(v_decided) == b.n_decided
        @test b.fdr_scope === :decided_subset_only
    end

    @testset "the all-abstain batch is a real answer, not an error" begin
        b0 = _dec_batch(0.10; ood = DEC_OOD_UNCHECKED, allow = false)
        @test b0.n_abstained == b0.n_total
        @test b0.n_decided == 0
        @test b0.k_star == 0
        @test b0.decided_fraction == 0.0
        @test isempty(b0.accepted)
        @test b0.posterior_expected_fdp == 0.0
        @test b0.abstain_reason_counts == (ood_not_checked = length(DEC_ITEMS),)
        # AND THE HEADLINE STILL CARRIES BOTH NUMBERS, which is the whole point: a rule that
        # abstains on everything hits any level trivially.
        @test occursin("0/$(b0.n_total)", p14_headline(b0))
    end

    @testset "the prior-sensitivity curve is required output and it varies" begin
        b = DEC_BATCH
        @test length(b.pi_sensitivity) == length(P14_PI_COLOC_GRID)
        @test Tuple(e.pi_coloc for e in b.pi_sensitivity) == P14_PI_COLOC_GRID
        @test all(e -> all(k -> haskey(e, k), P14_PI_SENSITIVITY_KEYS), b.pi_sensitivity)
        @test all(e -> haskey(e, :decided_fraction), b.pi_sensitivity)

        # PREVALENCE MISMATCH IS THE BIGGEST WAY THE GUARANTEE FAILS SILENTLY: the default prior is
        # the SIMULATOR's class mix, not any batch's prevalence, and the null posteriors move with
        # it. A curve whose rungs were all identical would mean the sweep was not reaching the
        # arithmetic -- which is a FINDING, so it is recorded rather than passed over.
        ks = [e.k_star for e in b.pi_sensitivity]
        if length(unique(ks)) == 1
            @info "the prior-sensitivity curve is flat on this fixture (k_star = $(ks[1]) at " *
                  "every rung); that is a finding about the fixture, not a pass"
        end
        @test length(unique(ks)) > 1
        @test b.class_prior_source === :measured_pi_class_masses
        @test b.class_prior_used == DEC_B.prior

        # A caller-supplied prior is RECORDED as such, so a silently switched prior is not
        # possible.
        b2 = decide_coloc(DEC_B, DEC_ITEMS, DEC_OOD, DEC_QHAT; alpha_fdr = 0.70,
                          allow_unchecked_ood = true,
                          prior = (coloc = 0.5, random = 0.3, exclusion = 0.2))
        @test b2.class_prior_source === :caller_supplied
    end

    @testset "validation names the function it failed in" begin
        # The level is a per-call USER parameter (SC1) and the open interval is the whole of it.
        for bad in (0.0, 1.0, -0.1, 1.5, NaN)
            e = try; _dec_batch(bad); nothing; catch err; err; end
            @test e isa ArgumentError
            @test occursin("decide_coloc", sprint(showerror, e))
        end
        # It is a REQUIRED keyword: there is no constant it could quietly default to.
        @test_throws UndefKeywordError decide_coloc(DEC_B, DEC_ITEMS, DEC_OOD, DEC_QHAT)
        @test :alpha_fdr in Base.kwarg_decl(
            first(m for m in methods(decide_coloc) if m.sig.parameters[3] === AbstractVector))

        # An empty batch is not a decided fraction.
        @test_throws ArgumentError decide_coloc(DEC_B, NamedTuple[], DEC_OOD, DEC_QHAT;
                                                alpha_fdr = 0.1)
        # A misspelled optional key would silently disable a channel, and a disabled cross-method
        # channel reports no disagreement, which reads as AGREEMENT.
        @test_throws ArgumentError decide_coloc(
            DEC_B, [merge(DEC_ITEMS[1], (mcis = nothing,))], DEC_OOD, DEC_QHAT; alpha_fdr = 0.1)

        # The src-shaped entry validates its channel arity FIRST, before any artifact is opened.
        for chans in ([1], [1, 2, 3], Int[])
            e = try
                decide_coloc(DEC_B, DEC_IMG_S, DEC_IMG_C, chans;
                             ood_nulls = DEC_OOD, qhat = DEC_QHAT)
                nothing
            catch err; err; end
            @test e isa ArgumentError
            @test occursin("decide_coloc", sprint(showerror, e))
            @test occursin("exactly 2 channels", sprint(showerror, e))
        end
        # Out of range, and the degenerate self-comparison.
        @test_throws ArgumentError decide_coloc(DEC_B, DEC_IMG_S, DEC_IMG_C, [1, 9];
                                                ood_nulls = DEC_OOD, qhat = DEC_QHAT)
        @test_throws ArgumentError decide_coloc(DEC_B, DEC_IMG_S, DEC_IMG_C, [2, 2];
                                                ood_nulls = DEC_OOD, qhat = DEC_QHAT)
    end

    @testset "the headline carries the level and its denominator in ONE string" begin
        b = DEC_BATCH
        h = p14_headline(b)
        @test occursin("DECIDED SUBSET ONLY", h)
        @test occursin(string(b.alpha_fdr), h)
        @test occursin("$(b.n_decided)/$(b.n_total)", h)
        @test occursin(string(round(100 * b.decided_fraction; digits = 1)), h)
        # The eight named limits travel with the batch rather than being remembered.
        @test length(b.meta.named_limits) == length(p14_named_limits())
        @test any(s -> occursin("NOT shipped in v2.0", s), b.meta.named_limits)
        @test any(s -> occursin("DECIDED SUBSET ONLY", s), b.meta.named_limits)
    end

    @testset "the signature is SRC-SHAPED, asserted rather than described" begin
        # D-01: a bundle plus images in, an `AbstractColocResult` subtype out -- the shape
        # `colocalization_amortized(img, control, channels; ...)` has, so a later promotion would
        # be mechanical. THE ASSERTION IS THE CLAIM; the prose in the docstring is not.
        @test hasmethod(decide_coloc,
                        Tuple{NamedTuple, MultiChannelImage, MultiChannelImage, Vector{Int}})
        @test P14Result <: AbstractColocResult
        @test AbstractColocResult === supertype(P14Result)

        # ... and "src-shaped" is NOT "in src". The decision layer is NOT shipped in v2.0.
        @test P14_SRC_UNTOUCHED == true
        @test isfile(joinpath(P14_REPO_ROOT, "spike", "p14", "decide.jl"))
        @test !isfile(joinpath(P14_REPO_ROOT, "src", "amortized", "decide.jl"))

        # ONE LIVE RUN THROUGH THE WHOLE SRC-SHAPED PATH, so SC2-d and Pitfall 7 are asserted on a
        # RESULT and not only in a source grep: the frozen ingestion, the net, the classical
        # comparison rescaled onto the cut's own scale, and the Costes channel on the Phase-14 dev
        # stream.
        r = decide_coloc(DEC_B, DEC_IMG_S, DEC_IMG_C, [1, 2];
                         ood_nulls = DEC_OOD_UNCHECKED, qhat = DEC_QHAT, idx = 7,
                         allow_unchecked_ood = true)
        @test r isa P14Result
        @test r isa AbstractColocResult
        @test decision(r) in P14_DECISIONS
        @test cross_method(r).basis === :ghat_rho_true_scale
        @test cross_method(r).costes_seed == UInt64(P14_DEV_SEED)
        @test cross_method(r).tau == DEC_B.tau
        @test cross_method(r).classical_call isa ThreeWayClass

        # SC2-d ASSERTED ON THE VALUE, NOT ON THE LABEL. Recorded because it was FALSIFIED during
        # execution: dropping the rescaling from ONE of the two channels is caught by NEITHER the
        # lane grep (its per-line rule fires only on a line that mentions both the statistic and
        # the cut, and the two live on different lines here) NOR by a range check (a raw induced
        # mean also lies in [-1, 1]) NOR by the recorded basis (a literal, which a wrong
        # computation carries just as happily). So the rescaled value is recomputed here from the
        # same frozen ingestion and compared, on BOTH channels, and the rescaling is asserted to
        # have actually moved the number -- an identity-looking pass would prove nothing.
        pair_s = _p14_channel_pair(DEC_IMG_S, [1, 2])
        pair_c = _p14_channel_pair(DEC_IMG_C, [1, 2])
        mu_s = patch_correlation(pair_s)
        mu_c = patch_correlation(pair_c)
        @test cross_method(r).rho_hat_sample  == ghat(mu_s)
        @test cross_method(r).rho_hat_control == ghat(mu_c)
        @test abs(cross_method(r).rho_hat_sample  - mu_s) > 1e-6
        @test abs(cross_method(r).rho_hat_control - mu_c) > 1e-6
        @test -1.0 <= cross_method(r).rho_hat_sample <= 1.0
        @test -1.0 <= cross_method(r).rho_hat_control <= 1.0
        # And the classical call really is the frozen rule applied at the INHERITED cut, on that
        # rescaled scale -- not a second cut invented here.
        @test cross_method(r).classical_call ===
              three_way_label(ghat(mu_s), ghat(mu_c); tau = DEC_B.tau)
        # Manders travels as a DESCRIPTIVE column; no call is derived from it.
        @test cross_method(r).manders_sample isa NamedTuple
        @test !haskey(cross_method(r), :manders_call)
        @test r.meta.channels == [1, 2]
        @test r.meta.basis_provenance === :p13_inherited_phase11_zt

        # D-04: the per-tile map is UNCONTROLLED DISPLAY, never an FDR-controlled call.
        @test p14_local_map_is_controlled(r) === false
        @test local_map(r) === nothing

        # A forbidden stream is refused at RUN TIME, not only at review time -- the grep checks
        # spellings and this checks values.
        @test_throws ArgumentError p14_cross_method(DEC_IMG_S, DEC_IMG_C, COLOC, DEC_B.tau;
                                                    idx = 1, seed = first(_p13_forbidden()))
    end

    @testset "src/ and the frozen environment are unchanged (D-01)" begin
        _has_git = Sys.which("git") !== nothing &&
                   (isdir(joinpath(P14_REPO_ROOT, ".git")) ||
                    isfile(joinpath(P14_REPO_ROOT, ".git")))
        if _has_git
            @test success(Cmd(`git diff --quiet HEAD -- src`; dir = P14_REPO_ROOT))
            @test isempty(readchomp(Cmd(`git status --porcelain -- src`; dir = P14_REPO_ROOT)))
            @test success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                              dir = P14_REPO_ROOT))
        else
            @info "git unavailable — skipping the executable src/ byte-equality assertions"
        end
    end

    @testset "the decision layer ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
