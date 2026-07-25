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

# spike/validation/p11_consts.jl --- Phase-11 Tier-1 pre-registration (D-04).
#
# THE ANTI-SNOOPING CONTRACT (D-04). Every ladder, every N, every pass/fail
# threshold and the DEV seed that governs the registration/chromatic-uncertainty
# work is LOCKED in this ONE file and committed BEFORE anything runs -- before the
# simulator edit, before the D-06 pre-flight probe, and before any training.
# Changing a value below after a run would be "tune until it passes" data-snooping;
# these constants ARE the committed pre-registration.
#
# TWO-TIER STRUCTURE -- READ THIS BEFORE EDITING ANYTHING.
#   TIER 1 is the single guard block in THIS file, committed before the probe.
#   TIER 2 (SC2_SPEARMAN_FLOOR and the Δρ_eq calibration coefficients) is APPENDED
#     LATER as a SECOND guard block keyed on a DIFFERENT sentinel const, carrying a
#     one-line provenance comment per constant naming the probe artifact it came from.
#   NO CONSTANT IN EITHER BLOCK IS EVER *EDITED*. Constants are only ever APPENDED,
#   so the git history of this file is itself the pre-registration audit trail. A
#   diff that MODIFIES a line below is, by construction, a pre-registration breach.
#
# PROVENANCE DISCLOSURE (the gate_consts_8_v2.jl:15-26 voice). An indicative,
# throw-away probe was run while researching this phase, so some magnitudes below
# (the ladders' ranges, the ±3 pp band arithmetic) were chosen with an approximate
# sense of scale already in hand. That is precisely why every value here is argued
# ONLY from properties of the DESIGN -- interpolation-onset physics, Wilson-interval
# arithmetic, the coverage-to-width algebra, conventional α levels -- and NEVER from
# what the indicative run scored or from what any net would pass. The keys that
# indicative run consumed are listed as BURNED below and are forbidden here.
# NOTHING REPORTED HAS BEEN RUN AGAINST THIS FILE. It is a specification, frozen
# before use.
#
# SEED DISCIPLINE (D-01). P11_DEV_SEED is a FRESH Random123 stream, EXECUTABLY proven
# disjoint from every reserved stream in the repository: the spike training stream
# (NPE_MASTER_SEED), the spike validation stream (VAL_MASTER_SEED), the fixture stream
# (VAL_FIX_SEED), the productionization datagen stream (DEFAULT_MASTER_SEED), the
# corpus stream (CORPUS_MASTER_SEED), every prior DEV seed, the keys this phase's
# research probe burned, and all four seeds of BOTH frozen gate families (PROD_SEED,
# PROD_SEED_V2) -- which are DERIVED, not literal, and so are RECOMPUTED here rather
# than trusted from a comment.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local constants; reaches src/ not at
# all and test/gate/ only READ-ONLY, through an isolated module.
#
# Guarded as ONE Tier-1 block keyed on :P11_DEV_SEED so a re-include under a test file
# that also loads p11_stats.jl is a silent no-op -- a redefinition to the same value
# would otherwise warn on a `const`.

# The frozen amended pre-registration, loaded into an ISOLATED module: it defines the
# same const names as `gate_consts_8.jl`, and we want ONLY its imsize mixture (and, for
# the unit tests, its frozen Holm adjuster) from it. Nothing is run against it.
#
# THIS SITS OUTSIDE THE GUARD BLOCK ON PURPOSE: Julia rejects a `module` expression
# that is not at top level ("syntax: \"module\" expression not at top level"), so it
# cannot live inside the guard block below. Idempotency is instead provided by the
# guarded-include idiom every caller uses:
#     isdefined(@__MODULE__, :P11_DEV_SEED) || include(".../p11_consts.jl")
module GateV2
    include(joinpath(@__DIR__, "..", "..", "test", "gate", "gate_consts_8_v2.jl"))
end

if !isdefined(@__MODULE__, :P11_DEV_SEED)
    import Random123: Philox4x     # counter-based RNG; the fresh disjoint DEV stream

    # =====================================================================================
    # 1. FORBIDDEN seeds (D-01) -- every stream a Phase-11 run must NOT consume
    # =====================================================================================
    # Literals below were read from their defining sites, not from any summary document.
    const NPE_MASTER_SEED     = 0x0000_0000_00C0_FFEE  # spike TRAINING stream (npe/train_npe.jl:65)
    const VAL_MASTER_SEED     = 0x0000_0000_5BC0_FFEE  # spike VALIDATION stream (consts.jl:76)
    const VAL_FIX_SEED        = 0x0000_0000_00F1_F7ED  # spike FIXTURE stream (consts.jl:79 = 0xF1F7ED)
    const DEFAULT_MASTER_SEED = 0x0000_0000_0000_0001  # productionization DATAGEN stream
    const CORPUS_MASTER_SEED  = 0x0000_0000_00C0_5EED  # corpus stream (corpus/manifest.csv, D-17)

    # DEV seeds burned by the F2 bounded-θ remedy diagnostic (07-CALIBRATION-FINDINGS).
    const F2_DEV_SEEDS    = (0x0000_0000_0DE7_C0DE, 0x0000_0000_DE7C_0DE2)
    # DEV seeds burned by earlier spike work (de-risk / smoke / ablation streams).
    const SPIKE_DEV_SEEDS = (0x0000_0000_00CA_11B0, 0x0000_0000_00CA_11B1,
                             0x0000_0000_00CA_11B2, 0x0000_0000_00CA_11B3,
                             0x0000_0000_00DE_C0DE, 0x0000_0000_0DE7_C0D3)
    # Keys consumed by this phase's own INDICATIVE research probe (2026-07-25). They have
    # already been observed against the simulator, so they are BURNED and not reusable --
    # reusing one would let a pre-observed stream govern a reported number (D-01).
    const P11_RESEARCH_BURNED = (0x0000_0000_0000_BEEF, 0x0000_0000_0000_CAFE,
                                 0x0000_0000_0000_FEED, 0x0000_0000_0000_DEAD)

    # The v1 and v2 ship-gate seeds are DERIVED, not literal. RECOMPUTE them from the
    # frozen mixing constants so they can be forbidden -- a comment naming a derived seed
    # is not evidence. Nothing here seeds anything; these values exist only to be excluded.
    const PROD_SALT    = 0x94D0_49BB_1331_11EB  # v1 gate salt (gate_consts_8_v2.jl:266)
    const PROD_MASTER  = 0x0000_0000_09E3_779B  # v1 gate master (gate_consts_8_v2.jl:267)
    const AMEND_SALT   = 0xC4CE_B9FE_1A85_EC53  # v2 gate salt (gate_consts_8_v2.jl:295)
    const AMEND_MASTER = 0x0000_0000_C0DE_2026  # v2 gate master (gate_consts_8_v2.jl:296)

    "One salted-Philox gate draw, the frozen per-grid seed rule (gate_consts_8_v2.jl:275-282)."
    _derive(master, salt, G) =
        rand(Philox4x(UInt64, (UInt64(master) ⊻ UInt64(salt), UInt64(G))), UInt64)

    "The v1 per-grid ship-gate seeds, recomputed purely so they can be forbidden."
    const PROD_SEED    = Dict{Int,UInt64}(G => _derive(PROD_MASTER, PROD_SALT, G)
                                          for G in (4, 8, 16, 32))
    "The v2 (amended) per-grid ship-gate seeds, recomputed purely so they can be forbidden."
    const PROD_SEED_V2 = Dict{Int,UInt64}(G => _derive(AMEND_MASTER, AMEND_SALT, G)
                                          for G in (4, 8, 16, 32))

    """
        _p11_forbidden() -> Tuple

    Every seed a Phase-11 stream must not be: zero, the four named reserved streams, the
    corpus stream, every prior DEV seed, the four keys the research probe burned, and all
    eight derived ship-gate seeds. `P11_DEV_SEED` is asserted absent from this set below
    AND, independently, in `spike/test/test_p11_consts.jl`.
    """
    _p11_forbidden() = (UInt64(0),
                        NPE_MASTER_SEED, VAL_MASTER_SEED, VAL_FIX_SEED,
                        DEFAULT_MASTER_SEED, CORPUS_MASTER_SEED,
                        F2_DEV_SEEDS..., SPIKE_DEV_SEEDS..., P11_RESEARCH_BURNED...,
                        values(PROD_SEED)..., values(PROD_SEED_V2)...)

    # =====================================================================================
    # 2. The FRESH Phase-11 DEV stream (D-01)
    # =====================================================================================
    # A distinct (seed, salt) pair gives a different Philox key, hence a disjoint stream;
    # the assertions at the foot of this block make the disjointness explicit rather than
    # merely overwhelmingly probable.
    const P11_SALT     = 0xA24B_AED4_663E_E121  # mx3 finalizer constant; ≠ HOLDOUT/FOLD/VAL/PROD/AMEND salts
    const P11_DEV_SEED = 0x0000_0000_0B11_DE71  # "0B11" = phase 11, "DE71" = DEV-1 (D-01)

    "The Phase-11 DEV RNG: `Philox4x` keyed by `(P11_DEV_SEED ⊻ P11_SALT, counter)`."
    p11_rng(counter::Integer = 0) =
        Philox4x(UInt64, (UInt64(P11_DEV_SEED) ⊻ P11_SALT, UInt64(counter)))

    # RESERVED counters, one per activity, so two activities can never silently share a
    # sub-stream. Mirrors the VAL_FIX_SEED discipline of `spike/validation/consts.jl:77-79`:
    # fixtures ride their own counter and NEVER consume a counter a reported number rides on.
    const P11_PROBE_COUNTER       = 1   # D-06 simulator-only pre-flight probe
    const P11_DATAGEN_COUNTER     = 2   # research-net training-pair generation
    const P11_LADDER_COUNTER      = 3   # SC2 monotone-widening / per-rung coverage ladder
    const P11_BREAKDOWN_COUNTER   = 4   # SC3 beyond-prior breakdown curve (reported, D-14)
    const P11_ATTENUATION_COUNTER = 5   # D-02 ρ_true attenuation A/B against the frozen net
    const P11_REALIMAGE_COUNTER   = 6   # D-17 qualitative real-image width check
    const P11_FIXTURE_COUNTER     = 99  # FIXTURES ONLY -- never a reported counter

    # =====================================================================================
    # 3. λ, the amortized registration-uncertainty level (D-03)
    # =====================================================================================
    # λ is the half-width of the shift prior the net is conditioned on, appended as an extra
    # input row so one training run yields a continuous, controlled sweep.
    const LAMBDA_MIN = 0.25  # D-03: never 0 -- the 0 → 0.25 step is INTERPOLATION ONSET, not
                             # misalignment sensitivity. At exactly integer offsets the backward
                             # map samples on-grid and BSpline(Linear()) is an identity lookup;
                             # any sub-pixel offset engages the kernel. A λ = 0 rung would be a
                             # qualitatively different (no-interpolation) regime.
    const LAMBDA_MAX = 3.0   # D-02: equals the widened SHIFT_PRIOR half-width, so the widest
                             # rung IS the stated prior -- the ladder never leaves it.

    # =====================================================================================
    # 4. Ladders (D-07, D-14)
    # =====================================================================================
    # SC2 / SC3-in-prior rungs. Seven rungs spanning exactly [LAMBDA_MIN, LAMBDA_MAX]; each
    # rung draws shift ~ Uniform(-λ_r, λ_r). No λ = 0 rung (see LAMBDA_MIN).
    const SC2_RUNGS = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)

    # SC3 BEYOND-PRIOR ladder (D-14): fixed-magnitude injection (D-15), λ pinned at read time.
    # The shift injection is along the DIAGONAL, i.e. dx = dy = x/√2, so |(dx,dy)| = x px.
    const SC3_SHIFT_RUNGS = (3.5, 4.0, 5.0, 6.0, 8.0)      # px, beyond the ±3 px training prior
    const SC3_EPS_RUNGS   = (0.025, 0.03, 0.04, 0.05)      # |chromatic_eps|, beyond the ±0.02 prior

    # The λ AT READ TIME for the beyond-prior curve. λ = 3.0 is the maximal-declared-uncertainty
    # reading: the user already declares the widest registration uncertainty the net knows and
    # STILL gets mis-registered data, so a breakdown there is the most damning honest number.
    # λ = 1.0 is the realistic under-declaration user story and will break earlier.
    # NEITHER IS GATED (D-14): the breakdown point is REPORTED, never a pass/fail input.
    const SC3_READ_LAMBDA_PRIMARY   = 3.0
    const SC3_READ_LAMBDA_SECONDARY = 1.0

    # =====================================================================================
    # 5. Coverage and equivalence (D-07, D-08)
    # =====================================================================================
    # D-07 is a TWO-PART criterion and is expressed here as constants, not as prose:
    #   (a) monotone non-decreasing posterior width across SC2_RUNGS -- gated by a one-sided
    #       permutation test on Spearman(λ, width) at P11_PERM_ALPHA plus the Tier-2 floor;
    #   (b) per-rung 90 % credible-interval coverage equivalent to nominal -- gated by TOST
    #       at SC2_TOST_ALPHA per side, Holm-corrected across rungs at P11_HOLM_FWER.
    # (b) is what makes the widening HONEST rather than merely larger.
    const SC2_COVERAGE_NOMINAL = 0.90   # the credible level whose coverage is under test
    const SC2_TOST_DELTA       = 0.03   # D-08 equivalence band, ±3 pp (derivation below)
    const SC2_TOST_ALPHA       = 0.05   # per-side α of the two one-sided tests
    const P11_HOLM_FWER        = 0.05   # family-wise error rate of the Holm step-down over rungs
    const P11_PERM_ALPHA       = 0.05   # one-sided permutation-test α for SC2(a)
    const PERM_B               = 10_000 # permutation replicates (rung-label shuffles)
    const Z_TWO_SIDED_90       = 1.644853627  # Φ⁻¹(0.95); the 90 % two-sided limb = TOST at α = 0.05/side
    const N_PER_RUNG           = 500    # paired draws per rung (pre-registered; see N_MIN_DERIVED)
    const N_MIN_DERIVED        = 271    # the knife-edge minimum the band arithmetic implies

    # DERIVATION OF N_MIN (recorded verbatim, the way gate_consts_8_v2.jl:189-198 records the
    # Fisher-z derivation). Require the 90 % Wilson interval at p̂ = SC2_COVERAGE_NOMINAL to fit
    # inside a ±SC2_TOST_DELTA band:
    #     half-width ≈ z·√(p(1−p)/N) = 1.644854·√(0.09/N) = 0.49346/√N ≤ 0.03
    #     ⇒ √N ≥ 16.449 ⇒ N ≥ 270.6 ⇒ N_min = 271
    # That is the KNIFE-EDGE case (p̂ exactly nominal). N_PER_RUNG = 500 buys slack: half-width
    # 0.0221, so equivalence stays DECIDABLE for p̂ ∈ [0.878, 0.922].
    #
    # DERIVATION OF THE ±3 pp BAND (D-08 requires an outcome-independent justification).
    # Step 1 -- convert a coverage band into a reported-width-error band. For a symmetric
    # posterior, if the reported 90 % interval is ±1.644854·c·σ while the true scale is σ, the
    # empirical coverage is 2Φ(1.644854·c) − 1:
    #     reported-width factor c | coverage | deviation from 0.90
    #     0.92                    | 0.8698   | −3.0 pp
    #     1.00                    | 0.9000   |  0
    #     1.10                    | 0.9260   | +2.6 pp
    # So a ±3 pp COVERAGE band ≈ a −8 % / +10 % error in the reported INTERVAL WIDTH.
    # Step 2 -- check that against the effect being demonstrated. The SC2 ladder is designed to
    # move the posterior width by a factor of order 3.7× (≈ +270 %) across [LAMBDA_MIN, LAMBDA_MAX].
    # The equivalence band is therefore ≈ 1/30 of the claimed effect -- small enough that
    # "coverage is provably close to nominal" cannot be confused with "coverage moved with the
    # ladder". Both steps are arithmetic on the design; neither reads a result.
    #
    # WHY TOST AND NOT A POINT NULL (the M = 2000 over-power lesson, 07-GATE-AMENDMENT). A
    # point-null test at large N rejects on drift too small to matter. TOST INVERTS the
    # incentive: larger N makes equivalence EASIER to establish, because the interval shrinks
    # toward the point estimate. There is no over-power pathology for an equivalence test, only
    # a cost ceiling -- so the Phase-7 lesson does not cap N here; it is the REASON D-08 chose
    # TOST in the first place.

    # =====================================================================================
    # 6. The D-06 pre-flight probe, and the abort criterion that must outrank it
    # =====================================================================================
    # The probe is simulator-only (no training): sweep shift and chromatic_eps through the
    # forward model at fixed θ and measure how far the 128-dim summary actually moves,
    # relative to the movement a known Δρ induces. Minutes, not hours.
    const P11_PROBE_THETA_BASE_N  = 5     # distinct base θ draws the sweep is repeated over
    const P11_PROBE_R             = 32    # paired replicates per (θ base, rung) cell
    const P11_PROBE_SHIFT_RUNGS   = (0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0)  # px, diagonal
    const P11_PROBE_SHIFT_EXT     = (4.0, 5.0, 6.0, 8.0)                  # px, D-14 range
    const P11_PROBE_EPS_RUNGS     = (0.0, 0.0025, 0.005, 0.01, 0.02)      # inside the ε prior
    const P11_PROBE_EPS_EXT       = (0.03, 0.05)                          # D-14 range
    const P11_PROBE_DRHO_RUNGS    = (0.02, 0.05, 0.10, 0.20)  # the Δρ reference scale
    # The metric is FIXED HERE so it cannot be selected to make a threshold convenient: the
    # mean AND median PAIRED 2-norm of the summary difference, over the 64 CONTINUOUS rows only.
    # NEVER the mask rows 65:128 -- they are near-constant and are exactly the rows
    # `_summary_row_partition` excludes; including them would dilute the norm with structural zeros.
    const P11_PROBE_METRIC        = :paired_l2_rows_1_64

    # THE ABORT CRITERION IS TIER 1 AND MUST BE IMMUNE TO THE PROBE'S OWN RESULT. If it were
    # set after seeing the probe, the "SC2 is below the fixed-summary resolution" branch would be
    # unfalsifiable. Declare "SC2 below resolution" -- do NOT train -- if either holds:
    #   (i)  S_probe = corspearman(λ rungs, mean Δρ_eq per rung) < P11_PROBE_S_FLOOR, or
    #   (ii) Δρ_eq(λ_max) − Δρ_eq(λ_min) < P11_PROBE_SPAN_FLOOR.
    # The span criterion is on the LADDER SPAN, not on every adjacent pair, because a local
    # plateau at the interpolation onset is a known simulator artifact, not a design failure.
    const P11_PROBE_S_FLOOR    = 0.9    # rank-correlation floor across the λ ladder
    const P11_PROBE_SPAN_FLOOR = 0.02   # Δρ_eq units; = the smallest calibrated Δρ step

    # The SC2(a) Spearman floor is SC2_SPEARMAN_ATTENUATION · S_probe. The ATTENUATION FACTOR is
    # Tier 1 -- fixed here, before the run; only the measured S_probe is Tier 2. Halving is a
    # round, conventional allowance for the fact that the width-vs-λ response is a COMPOSITION of
    # (summary displacement ← λ) with (posterior width ← summary), the second of which is a
    # learned, lossy stage. It is the same move 07-NUISANCE-SBC-SPEC-DRAFT.md:91-94 makes for
    # δ = 0.10 SD: a judgment allowance named in advance rather than discovered afterwards.
    const SC2_SPEARMAN_ATTENUATION = 0.5

    # =====================================================================================
    # 7. Image-size arm (F5): the research net TRAINS and the ladder EVALUATES on the mixture
    # =====================================================================================
    # BINDING INVARIANT (F5): SBC-style rank/coverage claims hold only under the joint the
    # estimator was TRAINED on, so train-joint MUST equal eval-joint. Both are therefore READ
    # from the frozen amended pre-registration through the isolated `GateV2` module above and
    # are NEVER retyped here -- a retyped mixture is a mixture that can silently drift.
    const P11_IMSIZE_SET     = GateV2.SBC_IMSIZE_SET
    const P11_IMSIZE_WEIGHTS = GateV2.SBC_IMSIZE_WEIGHTS

    # The one extra probe arm, kept ONLY so the pre-registered probe is comparable with the
    # indicative research run. It is a probe arm, never a training or ladder arm.
    const P11_PROBE_COMPARABILITY_IMSIZE = (256, 256)

    # =====================================================================================
    # 8. Budget
    # =====================================================================================
    const P11_N_PAIRS = 50_000  # matches every prior run (train_grid8_amended_v2.jl:11), so the
                                # comparison against the shipped net is CAPACITY-CONTROLLED.
    # Datagen wall-clock ceiling. Expected ≈ 1.2 h at 32 threads from the frozen mixture's own
    # cost model (gate_consts_8_v2.jl:148, E[cost] = 0.08346 s/pair ⇒ 50 000 pairs ≈ 4173 s).
    # EXCEEDING THIS IS A BLOCKER, RECORDED AS SUCH -- it is NOT licence to downgrade the arm to
    # 256², which would recreate the exact compute-budget artifact the F5 amendment
    # (gate_consts_8_v2.jl:126-151) exists to correct.
    const P11_DATAGEN_WALLCLOCK_CEILING_MIN = 150

    # =====================================================================================
    # 9. Stage-6 regression tolerance (D-10)
    # =====================================================================================
    # The composed affine warp must reproduce today's translation-only stage 6 EXACTLY at
    # chromatic_eps = 0. Bit-for-bit equality is MEASURED to hold (twelve Philox keys plus four
    # shift pairs), not hoped for: recenter ∘ LinearMap(I) composes to an exact identity linear
    # part with a zero translation, so the source coordinates -- and hence the B-spline
    # interpolation -- are identical.
    const P11_STAGE6_EXACT             = true   # `==` is the contract; do NOT soften to isapprox
    const P11_STAGE6_TOLERANCE_FALLBACK = 1e-13 # a few ULP; ANY use of this escape hatch must be
                                                # recorded as a deviation in the report, because it
                                                # means an upstream numerics change, not a pass.

    # =====================================================================================
    # 10. The D-04 novelty: the iteration allowance, declared in advance
    # =====================================================================================
    # Phase 7 was amended TWICE after seeing results, and the credibility cost of that is the
    # reason this constant exists. ONE documented iteration is authorised, and it is authorised
    # HERE, before any result exists. A SECOND iteration is NOT authorised by this file and
    # cannot be authorised by amending this file. The single allowance is cheapest spent on a
    # RE-PARAMETERISED PROBE (a wider λ ladder and/or the realistic image sizes) BEFORE training
    # -- never on relaxing a threshold after training.
    const P11_ITERATION_ALLOWANCE = 1

    # =====================================================================================
    # 11. Declared research-net deviations
    # =====================================================================================
    # Naming the model changes IN ADVANCE is the D-04 discipline: the Phase-7 precedent is that
    # UNNAMED model changes are how credibility gets spent. The report must enumerate exactly
    # these four and no others; anything else that changes is an undeclared deviation.
    const P11_RESEARCH_NET_DEVIATIONS = (
        "chromatic_eps appended as the 8th theta column (D-09), so D = 8 rather than 7",
        "SHIFT_PRIOR widened to Uniform(-3, 3) px -- SPIKE ONLY, deliberately NOT mirrored into src/ (D-02)",
        "lambda appended as a 129th input row, so d_in = 129 rather than 128 (D-03)",
        "BoundedThetaTransform ported from src/amortized/architecture.jl into the spike trainer (F2 remedy)",
    )

    # =====================================================================================
    # EXECUTABLE self-checks (cheap; no inference, no simulation, no file writes)
    # =====================================================================================
    # D-01: the fresh DEV seed is not any reserved, burned or derived gate seed.
    @assert !(UInt64(P11_DEV_SEED) in _p11_forbidden())
    @assert P11_SALT ∉ (0x9E37_79B9_7F4A_7C15,   # HOLDOUT_SALT
                        0xD1B5_4A32_D192_ED03,   # FOLD_SALT / CBS_SPLIT_SALT
                        0xBF58_476D_1CE4_E5B9,   # VAL_SALT
                        PROD_SALT, AMEND_SALT)
    # A duplicate inside the forbidden set would silently shrink it; catch that.
    @assert length(unique(_p11_forbidden())) == length(_p11_forbidden())
    # The recomputation reproduces the frozen file's OWN derived seeds -- proof that the
    # single-draw `_derive` above matches the gate's redraw-loop rule, so nothing is missed.
    @assert PROD_SEED_V2 == GateV2.PROD_SEED_V2
    @assert PROD_SEED    == GateV2.PROD_SEED
    # Ladder integrity: sorted, and spanning exactly the declared λ range.
    @assert issorted(SC2_RUNGS) && first(SC2_RUNGS) == LAMBDA_MIN && last(SC2_RUNGS) == LAMBDA_MAX
    @assert issorted(SC3_SHIFT_RUNGS) && first(SC3_SHIFT_RUNGS) > LAMBDA_MAX
    @assert issorted(SC3_EPS_RUNGS)
    @assert LAMBDA_MIN > 0.0
    # The pre-registered N clears the band arithmetic it was derived from.
    @assert N_PER_RUNG >= N_MIN_DERIVED
    @assert 0.0 < SC2_TOST_DELTA < SC2_COVERAGE_NOMINAL
    # The F5 mixture read succeeded and is a proper mixture.
    @assert length(P11_IMSIZE_SET) == length(P11_IMSIZE_WEIGHTS)
    @assert isapprox(sum(P11_IMSIZE_WEIGHTS), 1.0; atol = 1e-12)
    @assert !(P11_PROBE_COMPARABILITY_IMSIZE in P11_IMSIZE_SET)  # probe-only arm, not a train arm
    # The read-time λ choices are on, or inside, the trained ladder.
    @assert SC3_READ_LAMBDA_PRIMARY in SC2_RUNGS && SC3_READ_LAMBDA_SECONDARY in SC2_RUNGS
    # The declared-deviation list is a fixed-length enumeration, not an open-ended note.
    @assert length(P11_RESEARCH_NET_DEVIATIONS) == 4
    @assert P11_ITERATION_ALLOWANCE == 1
    # Fixtures must never ride a counter a reported number rides on.
    @assert P11_FIXTURE_COUNTER ∉ (P11_PROBE_COUNTER, P11_DATAGEN_COUNTER, P11_LADDER_COUNTER,
                                   P11_BREAKDOWN_COUNTER, P11_ATTENUATION_COUNTER,
                                   P11_REALIMAGE_COUNTER)
end

# =========================================================================================
# TIER 2 -- PROBE-DERIVED CONSTANTS. APPENDED, NEVER EDITED.
# =========================================================================================
#
# PROVENANCE DISCLOSURE (the gate_consts_8_v2.jl:15-26 voice). Everything in this block is a
# MEASUREMENT, taken from exactly one artifact:
#
#     spike/validation/p11_probe_report.jld2
#     generated   = 2026-07-25T21:26:36.787Z
#     elapsed_min = 7.108699997266133
#     produced by spike/validation/run_p11_probe.jl on the RESERVED stream
#                 p11_rng(P11_PROBE_COUNTER = 1) off P11_DEV_SEED, in ONE reported run
#
# THE TIER-1 BLOCK ABOVE IS BYTE-UNCHANGED since its own commit (d336699) -- this file grew,
# it did not move. That is checkable mechanically (`git diff` on the appending commit shows
# additions only) and it is re-asserted as literals in spike/test/test_p11_consts.jl, so an
# accidental Tier-1 edit during an append breaks the suite loudly rather than quietly rewriting
# the pre-registration.
#
# NOTHING HAS BEEN TRAINED AGAINST THIS FILE. No net exists yet in Phase 11; the probe is
# simulator-only and consumed no trained model. These constants are therefore measurements of
# the FORWARD MODEL's sensitivity, frozen before the first training run, exactly as D-06
# requires ("the probe supplies the SC2 threshold; the threshold is not guessed").
#
# WHAT IS JUDGMENT AND WHAT IS MEASUREMENT, KEPT VISIBLY SEPARATE. Every attenuation allowance
# below is a TIER-1 constant fixed before the probe ran; only the multiplicand is Tier 2. That
# is why SC2_SPEARMAN_FLOOR is written as a PRODUCT and never as a collapsed decimal.
#
# THE ABORT CRITERION IS NOT RE-STATED HERE. `P11_PROBE_S_FLOOR` / `P11_PROBE_SPAN_FLOOR` are
# Tier 1 and stay Tier 1; this block records only the two raw numbers they are evaluated against,
# so the verdict is reproducible from the file rather than from a narrative.
#
# Guarded as a SECOND block keyed on :SC2_SPEARMAN_FLOOR -- a distinct sentinel from Tier 1's
# :P11_DEV_SEED, so the two tiers can never collide on re-include.

if !isdefined(@__MODULE__, :SC2_SPEARMAN_FLOOR)

    # --- The two raw numbers the Tier-1 abort criterion is evaluated against ---------------
    # report field "S_probe": corspearman(SC2_RUNGS, mean dRho_eq per rung) on the F5 mixture
    # arm. 1.0 means dRho_eq was strictly increasing across all seven SC2 rungs.
    const P11_PROBE_S_MEASURED    = 1.0
    # report field "ladder_span": maximum(dRho_eq) - minimum(dRho_eq) over SC2_RUNGS, F5 arm.
    # Measured 0.0469 against the Tier-1 floor of 0.02 -- above it, by a factor of ~2.3.
    const P11_PROBE_SPAN_MEASURED = 0.046884564903982365

    # --- The SC2(a) monotonicity floor ------------------------------------------------------
    # WRITTEN AS A PRODUCT ON PURPOSE. SC2_SPEARMAN_ATTENUATION (Tier 1, 0.5) is the judgment
    # allowance for the fact that the width-vs-lambda response is a COMPOSITION of (summary
    # displacement <- lambda), which is what the probe measured, with (posterior width <-
    # summary), which is a learned and lossy stage the probe cannot see. P11_PROBE_S_MEASURED
    # is the measurement. Collapsing the two into a single decimal would hide which half is
    # which, which is precisely the confusion D-04 exists to prevent.
    const SC2_SPEARMAN_FLOOR = SC2_SPEARMAN_ATTENUATION * P11_PROBE_S_MEASURED

    # --- The dRho_eq calibration axis, frozen -----------------------------------------------
    # report fields "drho_eq_slope_f5" / "drho_eq_intercept_f5" / "drho_eq_r2_f5": ordinary
    # least squares of the mean paired 2-norm on the four pre-registered dRho rungs,
    #     ||ds||_2 = P11_DRHO_EQ_SLOPE * dRho + P11_DRHO_EQ_INTERCEPT
    # inverted by the ladder as dRho_eq(x) = (x - intercept) / slope. REPORTING ONLY -- no
    # pass/fail reads these -- but frozen here so the ladder's dRho_eq axis cannot drift.
    # NOTE for anyone reading a ladder table: because the fit has a positive intercept, a ZERO
    # displacement maps to a slightly NEGATIVE dRho_eq (-intercept/slope = -0.0028). That is an
    # artifact of the affine inversion, not a negative effect.
    const P11_DRHO_EQ_SLOPE     = 6.458505483772493
    const P11_DRHO_EQ_INTERCEPT = 0.018126542853519556
    const P11_DRHO_EQ_R2        = 0.9999313512418424

    # report fields "drho_eq_slope_256" / "..._intercept_256" / "..._r2_256": the same fit on the
    # 256-squared comparability arm, kept ONLY so the pre-registered probe is comparable with the
    # indicative research run. It is a probe arm, never a training or ladder arm.
    const P11_DRHO_EQ_SLOPE_256     = 7.318546621771334
    const P11_DRHO_EQ_INTERCEPT_256 = 0.10002486401477406
    const P11_DRHO_EQ_R2_256        = 0.9997658248160909

    # --- Which image-size arm the ladder runs on --------------------------------------------
    # THIS CONSTANT RECORDS THAT THE PROBE DID NOT OVERTURN A TIER-1 CHOICE; it does not make
    # one. The F5 mixture was already locked in Tier 1 (P11_IMSIZE_SET / P11_IMSIZE_WEIGHTS,
    # read from the frozen gate file), and the sec-C10(c) branch-2 escape hatch -- move the
    # ladder to different image sizes -- would only have been taken had the abort criterion
    # fired. It did not: S_probe = 1.0 >= 0.9 and span = 0.0469 >= 0.02 on the mixture itself.
    const P11_LADDER_IMSIZE_ARM = :f5_mixture

    # --- The SC1g lambda-ablation tripwire factor -------------------------------------------
    # report field "lambda_ratio": dRho_eq(LAMBDA_MAX) / dRho_eq(LAMBDA_MIN) across SC2_RUNGS on
    # the F5 arm. The arithmetic, recorded so it can be re-derived from the artifact:
    #     raw measured ratio                            = 5.003943268109952
    #     attenuated by SC2_SPEARMAN_ATTENUATION (0.5)   = 2.501971634054976
    #     floored at 1.05                               = 2.501971634054976   (floor not binding)
    # The 1.05 floor exists so the tripwire can NEVER become vacuous: a small measured ratio,
    # halved, could otherwise land below 1.0 and turn the material inequality into a tautology.
    # NOTE ON STRICTNESS, stated in advance: this is a MATERIAL bar derived from the forward
    # model's own sensitivity, not a "conditioning is alive at all" bar. It supersedes the 1.15
    # placeholder in spike/test/test_lambda_ablation.jl, and it is deliberately harder to clear.
    # The measured ratio is if anything CONSERVATIVE: dRho_eq(LAMBDA_MIN = 0.25) sits on the
    # interpolation-onset plateau, which inflates the denominator and shrinks the ratio.
    const P11_LAMBDA_ABLATION_FACTOR = 2.501971634054976

    # --- EXECUTABLE self-checks (cheap; no inference, no simulation, no file writes) ---------
    @assert 0.0 <= SC2_SPEARMAN_FLOOR <= 1.0
    @assert SC2_SPEARMAN_FLOOR == SC2_SPEARMAN_ATTENUATION * P11_PROBE_S_MEASURED
    @assert -1.0 <= P11_PROBE_S_MEASURED <= 1.0        # it is a rank correlation
    @assert P11_PROBE_SPAN_MEASURED > 0.0
    @assert P11_DRHO_EQ_SLOPE > 0.0 && P11_DRHO_EQ_SLOPE_256 > 0.0
    @assert 0.0 <= P11_DRHO_EQ_R2 <= 1.0 && 0.0 <= P11_DRHO_EQ_R2_256 <= 1.0
    @assert P11_LAMBDA_ABLATION_FACTOR >= 1.05          # the tripwire can never be vacuous
    @assert P11_LADDER_IMSIZE_ARM === :f5_mixture
end
