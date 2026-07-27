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

# spike/test/test_p13_tau.jl --- the Phase-13 tau resolution-probe gate (D-06).
#
# WHAT THIS FILE IS FOR. tau is the hinge of the whole phase: it sets the three-way class
# boundary, hence every label, hence every head's training set, hence every gate number. So
# the probe that measures it has to be right BEFORE it is run, because the reported run
# spends a reserved stream and its result cannot be honestly re-taken. Every assertion below
# checks a property that would be SILENT if it were wrong -- the probe would still return a
# plausible curve and a plausible tau.
#
#   13-RESEARCH Pitfall 8 -- A PAIRED DESIGN. If the two arms shared a Philox key, every
#     nuisance and the latent field would be held fixed across a comparison, tau would come
#     out several-fold too small, the random band would be razor-thin, the random class would
#     be starved, and the confusion matrix would then look excellent for a reason that has
#     nothing to do with the network. Nothing about the output SHAPE would change. Testset 3
#     asserts the arm keys differ AND that the arms genuinely draw different images.
#   13-RESEARCH Pitfall 10 -- THE WRONG SIMULATOR. A tau measured against the pre-Phase-11
#     simulator describes an easier world than the net trains in, so freezing it would be a
#     pre-registration error. Testset 8 asserts the provenance guard REPORTS rather than
#     throws, which is what lets plan 13-10's runner quote it in a hard stop.
#   T-13-17 -- GRID EXTENSION AFTER A MISS. Testset 6 asserts the stopping rule on SYNTHETIC
#     curves and greps the executable source for any path that could widen the frozen grid.
#     Extending it after seeing a miss would make :below_resolution unfalsifiable.
#
# FIXTURE SCALE IS MANDATORY, AND EVERY NUMBER BELOW IS AT FIXTURE POWER. The suite runs at
# R = TAU_FIX_R (16 draws per arm) on a single small 128-squared frame off the FIXTURE stream
# P13_FIX_SEED, so it stays far under the latency budget and NEVER consumes -- and so never
# pre-observes -- the reported P13_DEV_SEED stream (D-01, T-13-01). The REPORTED curve is
# produced by plan 13-10 at P13_TAU_R = 400 over the frozen F5 image-size mixture; nothing
# measured here is a reported number, and the AUC values seen here must never be quoted as
# the summary's resolution.
#
# Mirrors test_bf.jl / test_p13_labels.jl / test_p13_alpha.jl: license header -> using Test
# -> guarded include of the unit under test -> fixtures computed ONCE at module level -> one
# outer @testset -> the CPU-only postscript. No training, no figure call.

using Test
using Statistics

# Unit under test: the D-06 probe (pulls p13/consts.jl -> simulator -> contract -> the
# read-only reuse of the hand-rolled tie-aware AUC).
isdefined(@__MODULE__, :tau_curve) || include(joinpath(@__DIR__, "..", "p13", "tau_probe.jl"))

# --- Fixture knobs. NOT pre-registered thresholds; see the banner ---------------------------
# R = 16 gives an AUC standard error of roughly 0.09, i.e. an order of magnitude coarser than
# the reported R = 400 (~0.03). That is why testset 5 is a directional SANITY GUARD with a
# wide allowance and not a monotonicity assertion.
const TAU_FIX_R = 16
# ONE small frame, deliberately NOT the frozen F5 mixture: the frozen set starts at 512
# squared and its 2048-squared tail would blow the suite's latency budget many times over.
# The probe takes the mixture positionally exactly so a fixture can substitute one.
const TAU_FIX_IMSIZE_SET     = ((128, 128),)
const TAU_FIX_IMSIZE_WEIGHTS = (1.0,)
const TAU_FIX_BOOTSTRAP      = 64

# BOTH arms are carved off P13_FIX_SEED, so the reported P13_DEV_SEED stream stays untouched.
const TAU_FIX_KEYS = tau_probe_arm_keys(; seed = P13_FIX_SEED)

# The full fixture keyword bundle, so no testset can accidentally call the probe at reported
# scale (R = 400 over the F5 mixture) and spend minutes plus the wrong stream.
const TAU_FIX = (; R = TAU_FIX_R, arm_keys = TAU_FIX_KEYS,
                 imsize_set = TAU_FIX_IMSIZE_SET, imsize_weights = TAU_FIX_IMSIZE_WEIGHTS,
                 bootstrap = TAU_FIX_BOOTSTRAP,
                 boot_rng = p13_fix_rng(P13_FIXTURE_COUNTER))

# --- Fixtures computed ONCE (the test_bf.jl:49-65 idiom) -----------------------------------
# The rho = 0 reference arm is drawn ONCE and threaded into all three tau_auc_at calls below,
# which both halves the simulation cost and exercises the `reference` keyword the reported
# runner relies on. The arms stay UNPAIRED within each comparison -- the reference arm rides
# keys.reference and the contrast arms ride keys.contrast, two disjoint keys.
const TAU_FIX_REF = summary_draws(0.0, TAU_FIX_R, TAU_FIX_KEYS.reference;
                                  imsize_set = TAU_FIX_IMSIZE_SET,
                                  imsize_weights = TAU_FIX_IMSIZE_WEIGHTS)

_tau_fix_at(delta) = tau_auc_at(delta; reference = TAU_FIX_REF, TAU_FIX...)

const TAU_FIX_BIG   = _tau_fix_at(0.5)
const TAU_FIX_TINY  = _tau_fix_at(0.02)
const TAU_FIX_BIGGER = _tau_fix_at(0.6)

# Executable source of the unit under test (the test_p13_alpha.jl:123-124 idiom): comment
# lines stripped, so a grep gate cannot be satisfied -- or defeated -- by prose.
const TAU_SRC   = read(joinpath(@__DIR__, "..", "p13", "tau_probe.jl"), String)
const TAU_LINES = filter(l -> !startswith(strip(l), "#"), split(TAU_SRC, '\n'))
const TAU_CODE  = join(TAU_LINES, '\n')

@testset "P13 tau resolution probe (D-06)" verbose = true begin

    @testset "fixtures never consume a reported stream (D-01)" begin
        # The explicit statement of the fixture-stream rule: a gate that drew from the
        # reported stream would pre-observe the very draws the reported tau is measured on.
        @test UInt64(P13_FIX_SEED) != UInt64(P13_DEV_SEED)
        @test TAU_FIX_KEYS.reference != tau_probe_arm_keys().reference
        @test TAU_FIX_KEYS.contrast  != tau_probe_arm_keys().contrast
        # P13_TAU is Tier-2 and must still be ABSENT: this whole suite has to be runnable
        # BEFORE the probe has measured anything (D-06).
        @test !isdefined(@__MODULE__, :P13_TAU)
    end

    @testset "m-bar ignores absent patches" begin
        # Rows 1:4 are the per-patch correlations, rows 5:8 the present/absent indicator.
        # The fourth patch is ABSENT, and encode_d01 imputes it as 0.0 -- so the value row
        # carries a 0.0 that is not a measurement.
        z = vcat([0.5, 0.1, 0.9, 0.0], [1.0, 1.0, 1.0, 0.0])
        @test isapprox(mbar_from_summary(z), 0.5)          # (0.5 + 0.1 + 0.9) / 3

        # An UNWEIGHTED mean would answer 0.375: the missing-to-zero imputation would pull
        # the statistic toward 0 by 0.125 on a four-patch grid. That bias has a DIRECTION --
        # it drags both probe arms toward the same point, deflating discriminability and so
        # inflating tau -- which is why the weighting is a correctness requirement, not a
        # refinement.
        @test isapprox(mean(z[1:4]), 0.375)
        @test !isapprox(mean(z[1:4]), mbar_from_summary(z))

        # An all-absent grid reports the degeneracy instead of a silent zero.
        @test isnan(mbar_from_summary(vcat([0.0, 0.0], [0.0, 0.0])))
        # A single present patch is that patch, with no imputed rows averaged in.
        @test isapprox(mbar_from_summary(vcat([0.7, 0.0], [1.0, 0.0])), 0.7)
        # An odd row count is not a valid encoding and must fail loudly.
        @test_throws ArgumentError mbar_from_summary([0.1, 0.2, 0.3])
    end

    @testset "determinism under a fixed key" begin
        a = summary_draws(0.0, 6, TAU_FIX_KEYS.reference;
                          imsize_set = TAU_FIX_IMSIZE_SET,
                          imsize_weights = TAU_FIX_IMSIZE_WEIGHTS)
        b = summary_draws(0.0, 6, TAU_FIX_KEYS.reference;
                          imsize_set = TAU_FIX_IMSIZE_SET,
                          imsize_weights = TAU_FIX_IMSIZE_WEIGHTS)
        # BITWISE identical, never isapprox: the counter-based key is the whole reason the
        # reported probe is reproducible across processes and thread counts, and an
        # approximate check would pass a version that had quietly started drawing from the
        # global RNG somewhere in the chain.
        @test a == b
        @test all(isfinite, a)

        c = summary_draws(0.0, 6, TAU_FIX_KEYS.contrast;
                          imsize_set = TAU_FIX_IMSIZE_SET,
                          imsize_weights = TAU_FIX_IMSIZE_WEIGHTS)
        @test a != c

        # A different counter block on the SAME key is also a different sub-stream -- this is
        # what keeps the two directional contrast arms from sharing a draw.
        d = summary_draws(0.0, 6, TAU_FIX_KEYS.reference; counter_base = 6,
                          imsize_set = TAU_FIX_IMSIZE_SET,
                          imsize_weights = TAU_FIX_IMSIZE_WEIGHTS)
        @test a != d

        # rho is PINNED, not drawn: a merge under the wrong field spelling would APPEND a
        # ninth field, leave the knob at its prior draw, and make both arms identically
        # distributed -- the curve would flatten to chance for a reason having nothing to do
        # with the summary's resolution.
        @test length(merge(sample_prior(p13_fix_rng(P13_FIXTURE_COUNTER)),
                           (; ρ_true = -0.5))) == 8
        @test merge(sample_prior(p13_fix_rng(P13_FIXTURE_COUNTER)),
                    (; ρ_true = -0.5)).ρ_true == -0.5
        @test_throws ArgumentError summary_draws(0.0, 0, TAU_FIX_KEYS.reference)
        @test_throws ArgumentError summary_draws(1.5, 2, TAU_FIX_KEYS.reference)
    end

    @testset "the two arms are unpaired (disjoint streams)" begin
        # THE SINGLE MOST IMPORTANT PROPERTY OF THIS PROBE. In a PAIRED design the two arms
        # would share a key, so at the SAME rho they would be elementwise EQUAL by
        # construction -- every nuisance and the whole latent field held fixed. That would
        # report the resolution of a counterfactual the deployed tool never has and would
        # UNDERSTATE tau by a large factor (13-RESEARCH Pitfall 8). Phase 11's D-06 probe IS
        # paired, deliberately, because it measures a displacement; this one must not be.
        @test TAU_FIX_KEYS.reference != TAU_FIX_KEYS.contrast
        arm_a = summary_draws(0.0, 6, TAU_FIX_KEYS.reference;
                              imsize_set = TAU_FIX_IMSIZE_SET,
                              imsize_weights = TAU_FIX_IMSIZE_WEIGHTS)
        arm_b = summary_draws(0.0, 6, TAU_FIX_KEYS.contrast;
                              imsize_set = TAU_FIX_IMSIZE_SET,
                              imsize_weights = TAU_FIX_IMSIZE_WEIGHTS)
        @test !all(arm_a .== arm_b)          # NOT elementwise equal: unpaired
        @test all(arm_a[i] != arm_b[i] for i in eachindex(arm_a))
        # The pre-registration records the design, and it must agree with the code.
        @test P13_TAU_DESIGN === :unpaired
        @test P13_TAU_STATISTIC === :mbar_mask_weighted
    end

    @testset "tau_auc_at returns a discriminability magnitude" begin
        @test TAU_FIX_BIG.auc >= 0.5 && TAU_FIX_BIG.auc <= 1.0
        # The MAX over the two directions: the atoms are asymmetric by about 2.5x, so the
        # resolution may be too, and D-05 assumes ONE symmetric tau.
        @test TAU_FIX_BIG.auc >= TAU_FIX_BIG.neg
        @test TAU_FIX_BIG.auc >= TAU_FIX_BIG.pos
        # Each direction is reported as max(A, 1 - A). m-bar DECREASES with rho on the
        # negative side, so the RAW AUC of the rho = -delta arm runs toward 0, not 1: a bar
        # applied to the raw value would score perfect separation as "worse than chance".
        @test TAU_FIX_BIG.neg >= 0.5 && TAU_FIX_BIG.pos >= 0.5
        @test isfinite(TAU_FIX_BIG.auc_median_boot)
        @test TAU_FIX_BIG.n_reference == TAU_FIX_R
        @test TAU_FIX_BIG.n_neg == TAU_FIX_R && TAU_FIX_BIG.n_pos == TAU_FIX_R
        @test TAU_FIX_BIG.n_degenerate == 0
        # A(0) is chance by construction; measuring it would report "the summary cannot
        # resolve zero from zero" as though it were a finding.
        @test_throws ArgumentError tau_auc_at(0.0)
        @test_throws ArgumentError tau_auc_at(-0.1)
    end

    @testset "a large delta is more separable than a tiny one" begin
        # A SANITY GUARD AT FIXTURE POWER, NOT THE REPORTED MEASUREMENT. R = 16 gives an AUC
        # standard error of roughly 0.09, so the 0.15 allowance is measurement slop at this
        # n and nothing else; strict monotonicity is NOT asserted, because a non-monotone
        # rung at this power would be noise rather than a finding. The REPORTED curve comes
        # from plan 13-10 at P13_TAU_R over the frozen F5 mixture, and only that curve may be
        # quoted as the summary's resolution. To say it once in plain lower case so a grep
        # gate can find it: what this testset checks is a sanity guard at fixture power, and
        # it is not the reported measurement.
        @test TAU_FIX_BIGGER.auc >= TAU_FIX_TINY.auc - 0.15
        @test TAU_FIX_BIG.auc >= TAU_FIX_TINY.auc - 0.15
        @test isfinite(TAU_FIX_TINY.auc)
        @test TAU_FIX_TINY.auc >= 0.5 && TAU_FIX_TINY.auc <= 1.0
    end

    @testset "the abort criterion is immune to its own result" begin
        # SYNTHETIC curves, so the stopping rule is verified without spending one simulated
        # draw -- let alone the reserved reported stream. Only `delta` and `auc` are needed.
        curve = [(delta = 0.02, auc = 0.5), (delta = 0.05, auc = 0.95)]
        r = tau_from_curve(curve; bar = 0.90)
        @test r.tau == 0.05
        @test r.status === :measured
        @test r.delta_index == 2
        # 0.05 is below the ghat knot spacing near zero (0.0825), so the sub-knot caution
        # fires -- as a PRINTED note, never an error: rounding tau up to a knot would hide a
        # possibly real finding (13-RESEARCH Pitfall 8's warning sign).
        @test r.sub_knot === true
        @test isapprox(ghat_knot_spacing_near_zero(), 0.0825; atol = 1e-12)

        # THE FIRST clearing delta, not the largest, and not the best-looking one.
        multi = [(delta = 0.02, auc = 0.91), (delta = 0.05, auc = 0.99)]
        @test tau_from_curve(multi; bar = 0.90).tau == 0.02

        # NO delta clears the bar -> a reportable finding, NOT an extended grid.
        r2 = tau_from_curve(curve; bar = 0.99)
        @test r2.tau === nothing
        @test r2.status === :below_resolution
        @test r2.bar == 0.99
        @test_throws ArgumentError tau_from_curve(NamedTuple[])

        # The frozen constant that forbids the workaround, asserted as a literal.
        @test P13_TAU_ABORT_EXTEND_GRID === false

        # SOURCE GREP: no executable line may mention the frozen grid together with a
        # growth operation. Extending the grid after seeing a miss would make
        # :below_resolution unfalsifiable, so the absence of such a path is checked
        # mechanically rather than promised in a comment (T-13-17).
        for op in ("push!", "append!", "vcat", "prepend!", "insert!", "resize!")
            @test !any(l -> occursin("P13_TAU_DELTA_GRID", l) && occursin(op, l), TAU_LINES)
        end
        # Nor may the grid be rebuilt from a range or a comprehension.
        @test !any(l -> occursin("P13_TAU_DELTA_GRID", l) && occursin("range(", l), TAU_LINES)
        @test !occursin("P13_TAU_DELTA_GRID...", TAU_CODE)
        # The abort constant is READ by the code, not merely described in prose.
        @test occursin("P13_TAU_ABORT_EXTEND_GRID", TAU_CODE)
    end

    @testset "the frozen knobs are the defaults" begin
        # Keyword defaults are not directly introspectable without reflection gymnastics, so
        # the check is that each frozen knob is REFERENCED in the executable source: a probe
        # that had drifted onto private defaults would stop mentioning the pre-registration
        # it is supposed to be scored against.
        for knob in ("P13_TAU_DELTA_GRID", "P13_TAU_R", "P13_TAU_AUC", "P13_TAU_BOOTSTRAP_B",
                     "P13_TAU_BOTH_DIRECTIONS", "P13_IMSIZE_SET", "P13_IMSIZE_WEIGHTS",
                     "P13_TAU_COUNTER", "P13_DEV_SEED", "P13_SALT")
            @test occursin(knob, TAU_CODE)
        end
        # The frozen values themselves, asserted here so a silent edit to consts.jl shows up
        # in the probe's own gate as well as in test_p13_consts.jl.
        @test P13_TAU_DELTA_GRID == (0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20)
        @test P13_TAU_AUC == 0.90
        @test P13_TAU_R == 400
        @test P13_TAU_BOOTSTRAP_B == 1000
        @test P13_TAU_BOTH_DIRECTIONS === true
        @test issorted(P13_TAU_DELTA_GRID) && first(P13_TAU_DELTA_GRID) > 0.0

        # NO trained network and NO frozen standardization anywhere on the executable path:
        # that is the sequencing win (13-RESEARCH E4) that lets this file exist before the
        # Phase-11 net does. The token is checked against the comment-stripped source.
        @test !occursin("standardize_summary", TAU_CODE)
        @test !occursin("load_npe", TAU_CODE)
        @test !occursin(".zt", TAU_CODE)
        # NO new dependency: the tie-aware AUC is REUSED, never re-implemented.
        @test occursin("roc_auc(", TAU_CODE)
        @test !occursin("function roc_auc", TAU_CODE)
        # Pitfall 11: no Phase-13 executable line may reach the sealed corpus.
        @test !occursin("open_sealed_holdout", TAU_CODE)
        @test !occursin("corpus/", TAU_CODE)
        @test !occursin("corpus\\", TAU_CODE)

        # The whole curve is computed and reported, never short-circuited: a reader must be
        # able to re-read tau at any bar from the reported artifact.
        cv = tau_curve(; grid = (0.02, 0.5), R = 6, arm_keys = TAU_FIX_KEYS,
                       bootstrap = 8, boot_rng = p13_fix_rng(P13_FIXTURE_COUNTER),
                       imsize_set = TAU_FIX_IMSIZE_SET,
                       imsize_weights = TAU_FIX_IMSIZE_WEIGHTS)
        @test length(cv) == 2                              # BOTH rungs, including the miss
        @test [r.delta for r in cv] == [0.02, 0.5]
        @test all(haskey(r, :auc_neg) && haskey(r, :auc_pos) && haskey(r, :auc_median_boot)
                  for r in cv)
        @test all(isfinite(r.auc) for r in cv)
        @test_throws ArgumentError tau_curve(; grid = ())
    end

    @testset "simulator provenance guard reports, never errors" begin
        # REPORTS, never throws. Plan 13-10's reported runner is the caller that turns a
        # `false` into a HARD STOP immediately before it spends the reserved stream; keeping
        # the check reportable is what lets this suite exercise it at all and lets the runner
        # quote the reason in its own error message (13-RESEARCH Pitfall 10).
        g = simulator_provenance_guard()
        @test g isa NamedTuple
        @test g.post_p11 isa Bool
        @test g.reason isa AbstractString
        @test !isempty(g.reason)
        @test g.theta_arity isa Integer
        @test g.has_chromatic_eps isa Bool
        # The two markers are READ from the simulator, not asserted to a literal, so this
        # gate stays green on either side of the Phase-11 simulator merge. On the merged
        # simulator the guard reports true, the theta arity is 8 and the shift-prior
        # half-width equals the frozen reference-lambda value.
        @test g.expected_halfwidth == float(P13_TAU_REFERENCE_LAMBDA_EXPECTED)
        if g.post_p11
            @test g.has_chromatic_eps
            @test g.theta_arity == 8
            @test isapprox(g.shift_halfwidth, g.expected_halfwidth; rtol = 1e-9)
        end
        # The pre-registration insists the probe be run against the post-merge simulator.
        @test P13_TAU_REQUIRE_POST_P11_SIMULATOR === true
    end

    @testset "P13 tau probe ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
