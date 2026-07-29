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

# spike/test/test_p12_result.jl --- the SPAT-05 per-region result-type gate.
#
# WHAT THIS FILE IS FOR. The type under test makes four claims that are cheap to state and easy to
# lose, so each is asserted here rather than left to a docstring:
#
#   1. THE TYPE IS A REAL SUBTYPE, AND THE src/ NAME IS STILL FREE. Phase 7 D-02 requires every new
#      result variant to slot in as a new `AbstractColocResult` subtype, never as bolted-on fields.
#      Defining the executable type in spike/ (D-01) must not weaken that -- so the subtype relation
#      is asserted against the ACTUAL supertype from src/results.jl, reached read-only, and the
#      sketch's own identifier `SpatialColocResult` is asserted to be STILL UNCLAIMED.
#   2. AN UNSCORABLE REGION IS NEVER A MEASURED ZERO. The sentinel / `NaN` / `ood` triple is
#      enforced by the constructor, in ALL THREE maps, and each half is falsified independently.
#      The case that matters in practice is a near-zero Delta-rho on a failed region: without the
#      pairing it reads as "measured no difference", which is a false negative presented as a
#      finding (T-7-04).
#   3. THE Delta-rho GUARD IS A STRUCTURAL IDENTITY, NOT A RANGE. Testset 6 REPLACES a former
#      range-guard testset, and the replacement is the point: it shows a rho map mislabelled as a
#      Delta-rho map being REJECTED, and shows in the same breath that the retired `[-2, 2]` bound
#      would have admitted it silently.
#   4. THE DERIVED SCALAR SAYS SO, AND AN UNCHECKED OOD FLAG SAYS SO. Both are asserted through the
#      RUNNING system (`@doc`, a predicate), not by reading the source.
#
# Mirrors test_p13_result.jl: license header -> using Test -> guarded include of the unit under test
# -> fixtures built ONCE -> one outer @testset. No training, no simulation, no figure call, and no
# reported stream is consumed (nothing here draws random numbers at all).

using Test

# Unit under test: the spike-lane per-region result type. It pulls src/results.jl read-only, the
# shared calibration surface, the Phase-12 model surface and the lattice arithmetic transitively.
isdefined(@__MODULE__, :SpatialColocResultSpike) ||
    include(joinpath(@__DIR__, "..", "p12", "result.jl"))

# Source text with whole-line comments stripped, so a value MENTIONED IN PROSE cannot be mistaken
# for a declaration. Ported (one line) from `spike/test/test_p12_decoupling.jl:71-72`; re-declaring
# a function is silent in Julia, and the two definitions are byte-identical.
_strip_comment_lines(src::AbstractString) =
    join(filter(l -> !startswith(strip(l), "#"), split(src, '\n')), '\n')

# Guarded: a sandbox without git, or a source tree that is not a checkout, skips the executable
# half rather than failing for an unrelated reason. `isdir(...) || isfile(...)` because inside a
# linked worktree `.git` is a FILE. `P12_REPO_ROOT` is READ from the Tier-1 pre-registration rather
# than re-derived (p12_consts.jl section 12: two derivations can drift, and only one of them is the
# pre-registered one). The name is prefixed so it cannot collide with test_p12_decoupling.jl's own
# `_HAS_GIT` when the whole suite runs in one flat top-level namespace.
const _P12RES_HAS_GIT = Sys.which("git") !== nothing &&
                        (isdir(joinpath(P12_REPO_ROOT, ".git")) ||
                         isfile(joinpath(P12_REPO_ROOT, ".git")))

const _P12RES_G = 8

# --- Fixtures built ONCE ------------------------------------------------------------------------
#
# EVERY HELPER BELOW IS PHASE-PREFIXED, AND THAT IS NOT A STYLE CHOICE. The spike lane is a FLAT
# top-level namespace: when the whole suite runs in one process every test file's helpers land in
# the same module. `test_p13_result.jl:77` already binds a helper called `_res`, declared as
# `_res(; coloc = 2.5, ...)`. A Phase-12 `_res(maps = default; ...)` would generate a
# zero-positional method with the IDENTICAL signature `Tuple{typeof(_res)}`, so the later include
# would SILENTLY OVERWRITE the earlier one -- verified, not assumed. Today's include order
# (`test_p12_suite.jl` runs first) happens to hide it, which is exactly what makes it worth
# removing: the defect would surface only when someone reorders the suite, and it would then look
# like a Phase-13 failure. This is the same collision class `test_p13_result.jl:58-66` guards
# `P13_REPO_ROOT` against.

"""
A fresh nine-matrix working bundle in which EVERY region is measured: sample `s`, control `c`,
difference `s - c`, a finite positive sd everywhere and every flag `false`.

`p12_region_maps` deliberately hands back a FULLY UNSCORABLE bundle (fail-closed: a region is a
refusal until something measures it), so a test that wants measured regions has to write them --
which is exactly what a real write path does.
"""
function _p12_measured_maps(; G::Int = _P12RES_G, s::Float64 = 0.6, c::Float64 = 0.2,
                            sd::Float64 = 0.1)
    maps = p12_region_maps(G)
    fill!(maps.region_rho_sample,  s)
    fill!(maps.region_rho_control, c)
    fill!(maps.region_delta_rho,   s - c)
    fill!(maps.region_sd,          sd)
    fill!(maps.region_sd_sample,   sd)
    fill!(maps.region_sd_control,  sd)
    fill!(maps.ood_sample,  false)
    fill!(maps.ood_control, false)
    fill!(maps.ood,         false)
    return maps
end

"A valid `meta`; every keyword is overridable so a testset can vary exactly one thing."
_p12meta(; arm = :car, r1_posterior_median = 0.5, derived_rho = 0.3, ood_threshold = nothing,
         net_path = "spike/npe/p12_minispike.jld2", delta_rho_available = true) =
    p12_result_meta(; arm = arm, r1_posterior_median = r1_posterior_median,
                    derived_rho = derived_rho, ood_threshold = ood_threshold,
                    net_path = net_path, delta_rho_available = delta_rho_available)

"Build a result from a bundle."
_p12res(maps = _p12_measured_maps(); draws = nothing, calibration = nothing,
        meta = _p12meta(), G::Int = _P12RES_G) =
    SpatialColocResultSpike(G, maps, draws, calibration, meta)

const P12_RES_FIX = _p12res()

@testset "P12 spatial result type and its three maps (SPAT-05)" verbose = true begin

    @testset "SpatialColocResultSpike is a real AbstractColocResult subtype (Phase-7 D-02)" begin
        # The supertype is the REAL one from src/results.jl, reached read-only -- not a spike-local
        # look-alike. If the read-only include ever broke, this assertion is what would notice.
        @test SpatialColocResultSpike <: AbstractColocResult
        @test isabstracttype(AbstractColocResult)
        @test isconcretetype(SpatialColocResultSpike)
        @test P12_RES_FIX isa AbstractColocResult

        # The four shared accessors are all DEFINED for it ...
        for f in (delta_rho, bayes_factor, is_ood, posterior_draws)
            @test hasmethod(f, Tuple{SpatialColocResultSpike})
        end
        # ... and, the half that `hasmethod` alone cannot show: each resolves to a method written
        # FOR THIS TYPE rather than falling through to the supertype's `_iface_error` stub. Without
        # this, `hasmethod` would be `true` for a type that implements nothing at all, and the
        # "the interface is genuinely implemented" claim would be vacuous.
        for f in (delta_rho, bayes_factor, is_ood, posterior_draws)
            @test which(f, Tuple{SpatialColocResultSpike}).sig.parameters[2] ===
                  SpatialColocResultSpike
        end
        # `bayes_factor` is the one that deliberately routes BACK to the documented interface
        # error: this phase computes no evidence quantity, and fabricating one -- a threshold
        # crossing, a tail probability relabelled as evidence -- would put a figure in front of a
        # reader that no measurement backs. The refusal is a method, not an omission.
        @test_throws ErrorException bayes_factor(P12_RES_FIX)
        let err = try; bayes_factor(P12_RES_FIX); catch e; e; end
            @test occursin("accessor interface", err.msg)
            @test occursin("SpatialColocResultSpike", err.msg)
        end

        # THE src/ SKETCH NAME IS STILL UNCLAIMED. Asserted as "the identifier is not bound in this
        # module" rather than as `SpatialColocResultSpike !== SpatialColocResult`, which could not
        # even be written without binding the name it is trying to prove free.
        @test !isdefined(@__MODULE__, :SpatialColocResult)

        # D-02's other half: the new variant does NOT bolt fields onto the shipped struct.
        @test fieldnames(AmortizedColocResult) == (:grid, :posterior, :delta_rho_draws,
                                                   :log_bayes_factor, :ood, :calibration, :meta)
    end

    @testset "the accessors return the maps the src sketch names" begin
        r = P12_RES_FIX
        # `===`, not `==`: the sketch names these as the maps themselves. A copy would silently
        # decouple a consumer's view from the result, and a transpose would silently rotate the
        # lattice -- both are invisible to an `==` on a symmetric fixture.
        @test delta_rho_map(r)   === r.region_delta_rho
        @test uncertainty_map(r) === r.region_sd
        @test size(delta_rho_map(r))   == (_P12RES_G, _P12RES_G)
        @test size(uncertainty_map(r)) == (_P12RES_G, _P12RES_G)
        # The two single-stack maps are exposed too -- that is the whole reason the identity guard
        # of testset 6 is available at all.
        @test size(r.region_rho_sample)  == (_P12RES_G, _P12RES_G)
        @test size(r.region_rho_control) == (_P12RES_G, _P12RES_G)
        @test size(r.region_sd_sample)   == (_P12RES_G, _P12RES_G)
        @test size(r.region_sd_control)  == (_P12RES_G, _P12RES_G)
    end

    @testset "an unscorable region is unmistakable (T-7-04), in ALL THREE maps" begin
        idxs = [3, 11, 27]
        maps = _p12_measured_maps()
        # `stack` HAS NO DEFAULT and every call site must name it: which stack failed is a fact
        # about the data, not a convention, and a defaulted answer would be a guess recorded as a
        # measurement.
        @test_throws UndefKeywordError p12_mark_unscorable!(_p12_measured_maps(), idxs)
        @test_throws ArgumentError p12_mark_unscorable!(_p12_measured_maps(), idxs; stack = :bogus)

        p12_mark_unscorable!(maps, idxs; stack = :both)
        r = _p12res(maps)
        for i in idxs
            @test r.region_delta_rho[i]   == P12_SPATIAL_SENTINEL
            @test r.region_rho_sample[i]  == P12_SPATIAL_SENTINEL
            @test r.region_rho_control[i] == P12_SPATIAL_SENTINEL
            @test isnan(r.region_sd[i])
            @test isnan(r.region_sd_sample[i])
            @test isnan(r.region_sd_control[i])
            @test r.ood[i] && r.ood_sample[i] && r.ood_control[i]
        end
        # A region NOT marked is untouched -- the marking is surgical, not global. (`≈` on the
        # difference, because `0.6 - 0.2` is not the literal `0.4` in binary floating point; the
        # two single-stack values below ARE exact, since they were filled rather than computed.)
        @test r.region_delta_rho[1] ≈ 0.6 - 0.2
        @test !r.ood[1] && r.region_sd[1] == 0.1
        @test r.region_rho_sample[1] == 0.6 && r.region_rho_control[1] == 0.2

        # THE UNION RULE, exercised one-sided: a difference is unscorable whenever EITHER side is.
        maps2 = _p12_measured_maps()
        j = 42
        p12_mark_unscorable!(maps2, [j]; stack = :sample)
        r2 = _p12res(maps2)
        @test r2.region_rho_sample[j] == P12_SPATIAL_SENTINEL
        @test isnan(r2.region_sd_sample[j])
        @test r2.ood_sample[j]
        @test r2.region_delta_rho[j] == P12_SPATIAL_SENTINEL   # propagated
        @test isnan(r2.region_sd[j])
        @test r2.ood[j]
        # ... while the CONTROL stack at that region remains a normal measured value. That is the
        # asymmetry the union rule encodes and the reason the two single-stack maps are exposed.
        @test r2.region_rho_control[j] == 0.2
        @test r2.region_sd_control[j]  == 0.1
        @test !r2.ood_control[j]
        @test r2.ood == (r2.ood_sample .| r2.ood_control)

        # --- THE FALSIFIERS, which are the point of this testset -----------------------------
        # (a) `ood` true where `region_sd` is finite: an unscorable region reported with a
        #     confident uncertainty.
        let m = _p12_measured_maps()
            p12_mark_unscorable!(m, [5]; stack = :both)
            m.region_sd[5] = 0.1
            @test_throws ArgumentError _p12res(m)
        end
        # (b) `ood` true where `region_delta_rho` is a small NON-sentinel number. THIS IS THE CASE
        #     THAT MATTERS IN PRACTICE: a near-zero Delta-rho from a failed region reads as
        #     "measured no difference" unless the pairing is enforced -- a false negative presented
        #     as a finding.
        let m = _p12_measured_maps()
            p12_mark_unscorable!(m, [5]; stack = :both)
            m.region_delta_rho[5] = 1e-9
            @test_throws ArgumentError _p12res(m)
        end
        # (c) `region_sd` carrying a `NaN` where `ood` is false: an uncertainty that silently went
        #     missing on a region still being reported as measured.
        let m = _p12_measured_maps()
            m.region_sd[7] = NaN
            @test_throws ArgumentError _p12res(m)
        end
    end

    @testset "P12_SPATIAL_SENTINEL has not drifted from LOCAL_MAP_SENTINEL" begin
        # READ AS TEXT, NOT INCLUDED. src/ must stay UNLOADED by the spike lane's result surface,
        # so the shipped literal is recovered from the source itself. The alternative -- trusting
        # the comment in spike/p12/result.jl that says the two agree -- is exactly the
        # "a comment is not evidence" pattern this project's seed-disjointness rule already
        # rejects, and it is why the ship-gate seeds are recomputed rather than quoted.
        path = joinpath(P12_REPO_ROOT, "src", "amortized", "local_map.jl")
        @test isfile(path)
        body = _strip_comment_lines(read(path, String))
        ms = collect(eachmatch(r"const\s+LOCAL_MAP_SENTINEL\s*=\s*([-+0-9.eE]+)", body))
        @test length(ms) == 1                       # exactly one declaration, so the read is unambiguous
        @test parse(Float64, ms[1].captures[1]) == P12_SPATIAL_SENTINEL
        # ... and the shipped docstring's rationale is still there, so a future edit cannot delete
        # the discipline and leave a bare literal for this test to agree with vacuously.
        raw = read(path, String)
        @test occursin("T-7-04", raw)
        @test occursin("never silently indistinguishable", raw)
    end

    @testset "constructor guards" begin
        G = _P12RES_G
        ok = _p12_measured_maps()
        m  = _p12meta()

        # grid must be >= 2: a 1x1 lattice carries no spatial structure to map.
        @test_throws ArgumentError SpatialColocResultSpike(
            1, fill(0.4, 1, 1), fill(0.1, 1, 1), fill(0.6, 1, 1), fill(0.1, 1, 1),
            fill(0.2, 1, 1), fill(0.1, 1, 1), falses(1, 1), falses(1, 1), falses(1, 1),
            nothing, nothing, m)

        # A non-square `region_delta_rho`.
        @test_throws ArgumentError SpatialColocResultSpike(
            G, fill(0.4, G, G - 1), ok.region_sd, ok.region_rho_sample, ok.region_sd_sample,
            ok.region_rho_control, ok.region_sd_control, ok.ood_sample, ok.ood_control, ok.ood,
            nothing, nothing, m)

        # A `region_rho_sample` entry OUTSIDE the [-1, 1] SANITY range -- with `region_delta_rho`
        # ADJUSTED so the structural identity still holds. Without that adjustment the identity
        # check fires first and this test would pass for the wrong reason, telling us nothing about
        # the sanity bound it claims to exercise.
        let mm = _p12_measured_maps()
            mm.region_rho_sample[3] = 1.5
            mm.region_delta_rho[3]  = 1.5 - mm.region_rho_control[3]   # identity preserved
            @test isapprox(mm.region_delta_rho[3],
                           mm.region_rho_sample[3] - mm.region_rho_control[3]; atol = 1e-12)
            @test_throws ArgumentError _p12res(mm)
        end

        # A `region_sd` entry of 0.0 on a region flagged as measured: an uncertainty of exactly
        # zero is not an uncertainty.
        let mm = _p12_measured_maps(); mm.region_sd[2] = 0.0; @test_throws ArgumentError _p12res(mm); end

        # `region_draws` whose leading two dimensions disagree with `grid`.
        @test_throws ArgumentError _p12res(ok; draws = zeros(G - 1, G, 5))
        # ... and the correct shape constructs and round-trips through the shared accessor.
        let r = _p12res(ok; draws = zeros(G, G, 5))
            @test size(posterior_draws(r)) == (G, G, 5)
            @test posterior_draws(_p12res(ok)) === nothing
        end

        # A `meta` missing a required key. The key set is EXACT-MATCH, so an extra key is rejected
        # too -- both halves, because a bag that silently accepts either is not a contract.
        @test_throws ArgumentError _p12res(ok; meta = Base.structdiff(m, NamedTuple{(:zscore_arm,)}))
        @test_throws ArgumentError _p12res(ok; meta = merge(m, (; unexpected = 1)))

        # DO NOT TEST THE [-2, 2] BOUND ON `region_delta_rho` AS A FALSIFIER, AND HERE IS WHY.
        # Given the identity plus [-1, 1] on BOTH single-stack maps, |Delta-rho| > 2 is
        # UNCONSTRUCTIBLE: any input that violates the difference's range violates the identity or
        # a single-stack range FIRST, so a `@test_throws` there would pass just as happily on a
        # build with the [-2, 2] check deleted -- a test that cannot fail. The bound is a
        # DELIBERATE, DOCUMENTED REDUNDANCY, kept so a future reader who meets [-2, 2] learns that
        # a difference of two correlations genuinely spans it and does not "fix" it to [-1, 1]. It
        # is not a live guard and this testset must not imply otherwise.
    end

    @testset "the three maps satisfy the structural identity" begin
        # THIS TESTSET REPLACES THE FORMER RANGE-GUARD TESTSET, AND THE REPLACEMENT IS THE POINT.
        # A range could never have done this job: the old [-2, 2] bound admitted rho in [-1, 1]
        # SILENTLY, so the one check that could have caught a rho map mislabelled as a Delta-rho
        # map is exactly what hid it -- and tightening was no fix either, because with a TRUE
        # Delta-rho the range genuinely IS [-2, 2]. A range cannot distinguish the two quantities;
        # an identity can, and only because all three maps are now shipped.
        r = P12_RES_FIX
        @test maximum(abs, r.region_delta_rho .- (r.region_rho_sample .- r.region_rho_control)) <=
              P12_IDENTITY_ATOL

        # Perturb the CONTROL map alone: the difference no longer matches, and construction fails.
        let mm = _p12_measured_maps()
            mm.region_rho_control[9] = 0.25
            @test_throws ArgumentError _p12res(mm)
        end

        # THE MISLABEL, WHICH IS THE FAILURE THE GUARD EXISTS FOR: the SAMPLE map handed over where
        # the DIFFERENCE was meant.
        let mm = _p12_measured_maps()
            mm.region_delta_rho .= mm.region_rho_sample
            @test_throws ArgumentError _p12res(mm)
            # ... and, in the same breath, the proof that the retired range could not have caught
            # it: every entry of that mislabelled map sits comfortably inside [-2, 2], and inside
            # [-1, 1] as well.
            @test all(v -> -2.0 <= v <= 2.0, mm.region_delta_rho)
            @test all(v -> -1.0 <= v <= 1.0, mm.region_delta_rho)
        end
        # The sign error, which a range is equally blind to.
        let mm = _p12_measured_maps()
            mm.region_delta_rho .= mm.region_rho_control .- mm.region_rho_sample
            @test_throws ArgumentError _p12res(mm)
            @test all(v -> -2.0 <= v <= 2.0, mm.region_delta_rho)
        end

        # The SANITY checks, beside the guard and explicitly NOT it. A result whose difference sits
        # at -1.9 / +1.9 with CONSISTENT single-stack maps still constructs -- so a future reader
        # who "corrects" the difference's bound to [-1, 1] breaks a test that explains why.
        let mm = _p12_measured_maps(s = 0.95, c = -0.95)      # difference = +1.9
            @test _p12res(mm).region_delta_rho[1] ≈ 1.9
        end
        let mm = _p12_measured_maps(s = -0.95, c = 0.95)      # difference = -1.9
            @test _p12res(mm).region_delta_rho[1] ≈ -1.9
        end
        # And the [-1, 1] sanity bound on a single-stack map is live in the other direction too.
        let mm = _p12_measured_maps()
            mm.region_rho_control[4] = -1.5
            mm.region_delta_rho[4]   = mm.region_rho_sample[4] - (-1.5)   # identity preserved
            @test_throws ArgumentError _p12res(mm)
        end

        # THE NO-CONTROL-STACK MODE IS A LEGAL STATE, and it announces itself. Every invariant
        # holds: the identity is vacuous (every region is a sentinel somewhere), the sentinel value
        # lies inside both sanity ranges, and the pairing is exact in all three maps.
        let mm = _p12_measured_maps()
            p12_mark_unscorable!(mm, eachindex(mm.region_rho_control); stack = :control)
            rs = _p12res(mm; meta = _p12meta(delta_rho_available = false))
            @test all(==(P12_SPATIAL_SENTINEL), rs.region_rho_control)
            @test all(isnan, rs.region_sd_control)
            @test all(rs.ood_control)
            @test all(rs.ood)                                  # by the union rule
            @test all(==(P12_SPATIAL_SENTINEL), rs.region_delta_rho)
            @test all(isnan, rs.region_sd)
            @test rs.meta.delta_rho_available === false
            # ... and the sample stack is STILL fully measured, which is what makes this a refusal
            # to difference rather than an absence of data.
            @test all(==(0.6), rs.region_rho_sample)
            @test !any(rs.ood_sample)
            # A single-stack result is NOT a Delta-rho of zero. The flag is what a consumer must
            # branch on; the all-sentinel map alone would be indistinguishable from a real map of
            # regions that all failed, which is a different statement.
            @test P12_RES_FIX.meta.delta_rho_available === true
        end
    end

    @testset "the scalar rho is labelled DERIVED (R-1)" begin
        # THE LABEL MUST BE RETRIEVABLE FROM THE RUNNING SYSTEM, not just present in the file: a
        # reader who takes this scalar as the estimand has been misled, and the docstring is where
        # they will look. (This is also why the accessors live at top level rather than inside the
        # file's guard block -- Julia 1.12 drops docstrings written inside an `if` block.)
        doc = string(@doc(delta_rho))
        @test !isempty(doc)
        @test occursin("DERIVED", doc)
        @test occursin("not a sampled parameter", doc)
        @test delta_rho(P12_RES_FIX) == P12_RES_FIX.meta.derived_rho
        @test delta_rho(_p12res(; meta = _p12meta(derived_rho = -0.11))) == -0.11
    end

    @testset "a missing operating point reads as NOT CHECKED" begin
        # `local_map.jl:88-91`: a `false` OOD flag with NO recorded threshold means NOT CHECKED,
        # never "in distribution". Downstream code must branch on `p12_ood_checked` BEFORE
        # interpreting an all-`false` OOD map as a clean bill of health -- an unchecked flag read
        # as a pass is an unearned safety claim, not a conservative default.
        r = _p12res(; meta = _p12meta(ood_threshold = nothing))
        @test !any(r.ood)                       # the map is all-clear ...
        @test is_ood(r) === false               # ... and so is the single-boolean summary ...
        @test p12_ood_checked(r) === false      # ... but nothing was ever checked.
        # With an operating point recorded, the same all-`false` map means something.
        @test p12_ood_checked(_p12res(; meta = _p12meta(ood_threshold = 3.5))) === true
        # And the conservative direction of the summary: one unscorable region makes the map OOD.
        let mm = _p12_measured_maps()
            p12_mark_unscorable!(mm, [17]; stack = :both)
            @test is_ood(_p12res(mm)) === true
        end
    end

    @testset "the meta records the z-scoring arm (R-3)" begin
        # An unrecorded standardization arm makes every later comparison unreadable, because the
        # two arms differ in a way no downstream artifact would otherwise show.
        @test P12_RES_FIX.meta.zscore_arm === P12_ZSCORE_ARM === :per_row
        @test P12_RES_FIX.meta.schema_version === 1
        @test keys(P12_RES_FIX.meta) === P12_RESULT_META_KEYS
    end

    @testset "src/ is byte-unchanged" begin
        # CLAUDE.md requires src/ to stay provably untouched during the spike, and the executable
        # subtype under test exists precisely so the sketch does not have to be edited. This
        # assertion is what makes "provably" mean something.
        if _P12RES_HAS_GIT
            @test success(Cmd(`git diff --quiet HEAD -- src`; dir = P12_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/results.jl`; dir = P12_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src/amortized/local_map.jl`;
                              dir = P12_REPO_ROOT))
            # THE UNTRACKED CASE, WHICH `git diff HEAD` CANNOT SEE. A NEW file dropped under src/
            # is invisible to a diff against HEAD and would pass every assertion above while being
            # exactly the productionization-by-stealth the constraint forbids.
            @test isempty(readchomp(Cmd(`git status --porcelain -- src`; dir = P12_REPO_ROOT)))
        else
            @info "git unavailable — skipping the executable src/ byte-equality assertions"
        end
        # The sketch block is still there, still comments, still carrying the unclaimed name -- so
        # the "superseded, not edited" claim is checkable from the running suite.
        let sketch = read(joinpath(P12_REPO_ROOT, "src", "results.jl"), String)
            @test occursin("SKETCH ONLY", sketch)
            @test occursin("delta_rho_map", sketch)
            @test occursin("uncertainty_map", sketch)
        end
    end

    @testset "result surface ran CPU-only" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
