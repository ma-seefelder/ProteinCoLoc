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

# spike/validation/p12_consts.jl --- Phase-12 Tier-1 pre-registration (R-5).
#
# THE ANTI-SNOOPING CONTRACT (R-5). Every ladder, every N, every pass/fail threshold and
# every seed the spatial-colocalization-map phase consumes is LOCKED in this ONE file and
# committed BEFORE anything runs -- before the lattice kernels exist, before the stage-1
# simulator edit, before the D-12 Stage-1 ridge, before a single training pair is drawn.
# Changing a value below after a run would be "tune until it passes" data-snooping; these
# constants ARE the committed pre-registration.
#
# TWO-TIER STRUCTURE -- READ THIS BEFORE EDITING ANYTHING.
#   TIER 1 is the single guard block in THIS file, keyed on :P12_DEV_SEED, committed before
#     anything in Phase 12 has run.
#   TIER 2 IS NOT ONE BLOCK. It is a SEQUENCE of append-only blocks, each keyed on its OWN
#     sentinel const, because a second append that reused an earlier block's sentinel would
#     be silently SKIPPED once that block exists, and an append with no guard at all is a
#     hard `const` redefinition on re-include. THREE are reserved here, with the plan that
#     opens each:
#       :P12_CHOSEN_PRIOR  -- opened by 12-15 (P12_CHOSEN_PRIOR, P12_K_PROD, P12_D_PROD)
#       :P12_N_LOW         -- opened by 12-15 (P12_N_LOW). CONFIRMED 2026-07-29 as option
#                             (ii): the identified/vacuous boundary in 12-18's SBC row
#                             classes is a MEASUREMENT off 12-15's truncation curve, not a
#                             Tier-1 guess. It gets its OWN sentinel rather than riding in
#                             the :P12_CHOSEN_PRIOR block even though the same plan writes
#                             both, because 12-15's NONE-BEATS-ABLATION branch appends
#                             NEITHER, and two separately-guarded blocks make "which branch
#                             fired" readable from the file instead of inferred.
#       :P12_FISHERZ_NEFF  -- opened by 12-16 (P12_FISHERZ_NEFF)
#   NO CONSTANT IN EITHER TIER IS EVER *EDITED*. Constants are only ever APPENDED, so the
#   git history of this file is itself the pre-registration audit trail. A diff that
#   MODIFIES a line below is, by construction, a pre-registration breach.
#
# REPORTED IS NOT GATED, AND THE DISTINCTION IS MACHINE-CHECKABLE. Only the D-12 Stage-1 and
# Stage-2 thresholds gate. The atom-mass, epsilon-identifiability, vacuous-column-shrinkage
# and S-4 / D-06 guard bars are declared in advance so they are COMPARABLE, but they decide
# nothing. That split is not left to prose: it is two disjoint constant-name tuples
# (P12_GATING_CONSTANTS / P12_REPORTING_ONLY_CONSTANTS) with an in-file emptiness assertion
# on their intersection and an independent test in spike/test/test_p12_consts.jl.
#
# SEED DISCIPLINE (R-5). P12_DEV_SEED is a FRESH Random123 stream, EXECUTABLY proven disjoint
# from every reserved stream in the repository: the spike training stream (NPE_MASTER_SEED),
# the spike validation stream (VAL_MASTER_SEED), the fixture stream (VAL_FIX_SEED), the
# productionization datagen stream (DEFAULT_MASTER_SEED), the shipped ratio-pairing stream
# (RATIO_PAIR_SEED), the corpus stream (CORPUS_MASTER_SEED), Phase 11's and Phase 13's DEV and
# FIXTURE streams, every DEV key any prior spike / diagnostic / research probe has already
# OBSERVED, and all eight seeds of BOTH frozen gate families (PROD_SEED, PROD_SEED_V2) --
# which are DERIVED, not literal, and so are RECOMPUTED here rather than trusted from a
# comment. Phase 11's pre-registration predates Phase 13 and therefore does not forbid it;
# this file closes that gap in the other direction.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local constants; reaches src/ not at all and
# test/gate/ only READ-ONLY, through an isolated module. Adds NO package: this file imports
# only Random123 and Distributions, both already in spike/Project.toml.
#
# Guarded as ONE Tier-1 block keyed on :P12_DEV_SEED so a re-include under a test file that
# also loads another phase's constants is a silent no-op.

# The frozen amended pre-registration, loaded into an ISOLATED module: it defines the same
# const names as `gate_consts_8.jl`, and we want ONLY its imsize mixture and its two derived
# ship-gate seed families from it. Nothing is run against it.
#
# THIS SITS OUTSIDE THE GUARD BLOCK ON PURPOSE: Julia rejects a `module` expression that is
# not at top level ("syntax: \"module\" expression not at top level"), so it cannot live
# inside the guard block below. Idempotency is instead provided by the guarded-include idiom
# every caller uses:
#     isdefined(@__MODULE__, :P12_DEV_SEED) || include(".../p12_consts.jl")
#
# (Same precedent as `module GateV2` in p11_consts.jl:75-77 and `module _GC` in
# spike/p13/consts.jl:96-98. NOTE THE PATH DEPTH: from spike/validation/ the repo root is two
# levels up.)
module GateV2P12
    include(joinpath(@__DIR__, "..", "..", "test", "gate", "gate_consts_8_v2.jl"))
end

if !isdefined(@__MODULE__, :P12_DEV_SEED)
    import Random123: Philox4x        # counter-based RNG; the fresh disjoint DEV stream
    import Distributions: Uniform     # the r1 prior is a bound distribution OBJECT, see §4

    # =====================================================================================
    # 1. FORBIDDEN seeds (R-5) -- every stream a Phase-12 run must NOT consume
    # =====================================================================================
    # Literals below were read from their defining sites, not from any summary document.
    const NPE_MASTER_SEED     = 0x0000_0000_00C0_FFEE  # spike TRAINING stream (npe/train_npe.jl:65)
    const VAL_MASTER_SEED     = 0x0000_0000_5BC0_FFEE  # spike VALIDATION stream (validation/consts.jl:76)
    const VAL_FIX_SEED        = 0x0000_0000_00F1_F7ED  # spike FIXTURE stream (validation/consts.jl:79 = 0xF1F7ED)
    const DEFAULT_MASTER_SEED = 0x0000_0000_0000_0001  # productionization DATAGEN stream
    const RATIO_PAIR_SEED     = 0x0000_0000_004A_7107  # shipped ratio pairing (src/amortized/train_ratio.jl:60)
    const CORPUS_MASTER_SEED  = 0x0000_0000_00C0_5EED  # corpus stream (corpus/manifest.csv, D-17)
    const P11_DEV_SEED        = 0x0000_0000_0B11_DE71  # spike/validation/p11_consts.jl:144
    const P13_DEV_SEED        = 0x0000_0000_0B13_DE71  # spike/p13/consts.jl:188
    const P13_FIX_SEED        = 0x0000_0000_0B13_F1F7  # spike/p13/consts.jl:189

    # DEV seeds burned by the F2 bounded-theta remedy diagnostic (07-CALIBRATION-FINDINGS),
    # copied from p11_consts.jl:93.
    const F2_DEV_SEEDS    = (0x0000_0000_0DE7_C0DE, 0x0000_0000_DE7C_0DE2)
    # DEV seeds burned by earlier spike work (de-risk / smoke / ablation streams),
    # copied from p11_consts.jl:95-97.
    const SPIKE_DEV_SEEDS = (0x0000_0000_00CA_11B0, 0x0000_0000_00CA_11B1,
                             0x0000_0000_00CA_11B2, 0x0000_0000_00CA_11B3,
                             0x0000_0000_00DE_C0DE, 0x0000_0000_0DE7_C0D3)
    # Keys consumed by Phase 11's INDICATIVE research probe (2026-07-25), copied from
    # p11_consts.jl:101-102. Pre-observed streams can never govern a reported number.
    const P11_RESEARCH_BURNED = (0x0000_0000_0000_BEEF, 0x0000_0000_0000_CAFE,
                                 0x0000_0000_0000_FEED, 0x0000_0000_0000_DEAD)

    # Every DEV key any prior spike, diagnostic or research probe has already OBSERVED, copied
    # VERBATIM (all 25 entries, in declaration order) from spike/p13/consts.jl:120-146.
    const P13_BURNED_DEV_SEEDS = (
        0x0000_0000_0DE7_C0DE,   # F2 bounded-theta remedy diagnostic (07-CALIBRATION-FINDINGS)
        0x0000_0000_DE7C_0DE2,   # F2 remedy, second key
        0x0000_0000_0DE7_C0D3,   # F2 remedy, third key
        0x0000_0000_00B7_A770,   # spike 006/008 BF-attrition keys (widened from 0x00B7A77x)
        0x0000_0000_00B7_A771,
        0x0000_0000_00B7_A772,
        0x0000_0000_00B7_A773,
        0x0000_0000_00CA_11B0,   # spike 010 calibration-vs-imsize keys
        0x0000_0000_00CA_11B1,
        0x0000_0000_00CA_11B2,
        0x0000_0000_00CA_11B3,
        0x0000_0000_00DE_C0DE,   # spike 013 mu-prior-truncation training key
        0x0000_0000_00BF_5014,   # spike 014 BF simulation-validation key
        0x0000_0000_00A1_CE01,   # misc burned spike / gate-script keys
        0x0000_0000_00A1_CE02,
        0x0000_0000_C0FF_EE5B,
        0x0000_00C0_FFEE_5BC0,
        0x0000_0000_00A7_0115,
        0x0000_0000_00A7_0213,
        0x0000_0000_0134_D8F3,
        0x0000_0000_0000_BEEF,   # the four 11-RESEARCH indicative-probe burns
        0x0000_0000_0000_CAFE,
        0x0000_0000_0000_FEED,
        0x0000_0000_0000_DEAD,
        0x0000_0000_0000_0007,   # the bare global-RNG seed the same probe used
    )

    # The v1 and v2 ship-gate seeds are DERIVED, not literal. RECOMPUTE them from the frozen
    # mixing constants so they can be forbidden -- a comment naming a derived seed is not
    # evidence. Nothing here seeds anything; these values exist only to be excluded.
    const PROD_SALT    = 0x94D0_49BB_1331_11EB  # v1 gate salt   (gate_consts_8_v2.jl:266)
    const PROD_MASTER  = 0x0000_0000_09E3_779B  # v1 gate master (gate_consts_8_v2.jl:267)
    const AMEND_SALT   = 0xC4CE_B9FE_1A85_EC53  # v2 gate salt   (gate_consts_8_v2.jl:295)
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
        _p12_forbidden_list() -> Vector{UInt64}

    The forbidden-seed list in declaration order, so an accidental duplicate stays visible.

    P13_BURNED_DEV_SEEDS IS THE BURNED-KEY SUPERSET AND IS LISTED ONCE. All twelve entries of
    `F2_DEV_SEEDS`, `SPIKE_DEV_SEEDS` and `P11_RESEARCH_BURNED` also appear inside it (they are
    the same historical burns, re-collected by Phase 13 into one tuple). Concatenating all four
    tuples would therefore put twelve keys in this list TWICE, which is exactly what the
    duplicate assertion at the foot of this block exists to catch -- the assertion would fire on
    a redundancy rather than on a real collision, and the whole check would have to be dropped
    to make the file load. So the three Phase-11 tuples are declared above (read from their
    defining sites, in the pre-registration record, and independently asserted in
    spike/test/test_p12_consts.jl) and their SUBSUMPTION by P13_BURNED_DEV_SEEDS is PROVED
    below rather than assumed. The forbidden SET is identical either way; only the double-count
    is removed.
    """
    _p12_forbidden_list() = vcat(
        UInt64[UInt64(0),
               UInt64(NPE_MASTER_SEED),     UInt64(VAL_MASTER_SEED),
               UInt64(VAL_FIX_SEED),        UInt64(DEFAULT_MASTER_SEED),
               UInt64(RATIO_PAIR_SEED),     UInt64(CORPUS_MASTER_SEED),
               UInt64(P11_DEV_SEED),        UInt64(P13_DEV_SEED),
               UInt64(P13_FIX_SEED)],
        collect(UInt64, P13_BURNED_DEV_SEEDS),
        # DERIVED, NEVER TRUSTED FROM A COMMENT: both frozen gate families are RECOMPUTED
        # above and forbidden from the recomputed values.
        collect(UInt64, values(PROD_SEED)),
        collect(UInt64, values(PROD_SEED_V2)),
    )

    """
        _p12_forbidden() -> Set{UInt64}

    Every seed a Phase-12 stream must not be: zero, the six named reserved streams, Phase 11's
    and Phase 13's DEV / FIXTURE streams, every prior burned spike key, and all eight
    RECOMPUTED ship-gate seeds. `P12_DEV_SEED` and `P12_FIX_SEED` are asserted absent from this
    set at the foot of this block AND, independently, in `spike/test/test_p12_consts.jl`.
    """
    _p12_forbidden() = Set{UInt64}(_p12_forbidden_list())

    # =====================================================================================
    # 2. The salt inventory, and the FRESH Phase-12 streams (R-5)
    # =====================================================================================
    # The salt inventory of the whole repository. A fresh salt must not be any of these,
    # because a repeated (seed, salt) pair is a repeated Philox key, hence a repeated stream.
    # The first eight entries are Phase 13's P13_REPO_SALTS (spike/p13/consts.jl:150-159),
    # copied unchanged.
    #
    # THE NINTH ENTRY IS A GAP PHASE 13 LEFT OPEN, AND IT IS A COLLISION WORTH RECORDING.
    # `0x2545_F491_4F6C_DD1D` is BOTH `P13_SALT` (spike/p13/consts.jl:192) AND
    # `P11_DATAGEN_SALT` (spike/data/p11_generate.jl:90) -- one literal, two phases, two
    # independent streams keyed off it. Phase 13's own inventory omits it (a fresh salt is
    # asserted against the inventory, so the inventory never lists the salt being minted), so
    # nothing in the repository would have caught a third phase reaching for the same value.
    # This is precisely the kind of shared key the inventory exists to surface, so it is
    # recorded here rather than left to be rediscovered.
    #
    # `P11_SALT` (0xA24B_AED4_663E_E121) needs no separate append: Phase 13 already carries it
    # as the sixth entry.
    const P12_REPO_SALTS = (
        0x9E37_79B9_7F4A_7C15,   # HOLDOUT_SALT      (spike/data/seeding.jl:45)
        0xD1B5_4A32_D192_ED03,   # FOLD / CBS_SPLIT  (spike/data/seeding.jl:46)
        0xBF58_476D_1CE4_E5B9,   # VAL_SALT
        0x94D0_49BB_1331_11EB,   # PROD_SALT         (gate_consts_8_v2.jl:266)
        0xC4CE_B9FE_1A85_EC53,   # AMEND_SALT        (gate_consts_8_v2.jl:295)
        0xA24B_AED4_663E_E121,   # P11_SALT          (p11_consts.jl:143)
        0xA5A5_A5A5_A5A5_A5A5,   # COSTES
        0xA5A5_A5A5_DEAD_BEEF,   # COSTES, second
        0x2545_F491_4F6C_DD1D,   # P13_SALT (p13/consts.jl:192) == P11_DATAGEN_SALT
                                 # (data/p11_generate.jl:90) -- ONE literal, TWO phases
    )

    # --- The FRESH Phase-12 streams -------------------------------------------------------
    # A distinct (seed, salt) pair gives a different Philox key, hence a disjoint stream; the
    # assertions at the foot of this block make the disjointness EXPLICIT rather than merely
    # overwhelmingly probable.
    const P12_DEV_SEED     = 0x0000_0000_0B12_DE71  # "0B12" = phase 12, "DE71" = DEV-1 (R-5)
    const P12_FIX_SEED     = 0x0000_0000_0B12_F1F7  # FIXTURE stream, mirroring the VAL_FIX_SEED
                                                    # discipline: fixtures must NEVER consume
                                                    # (and so never pre-observe) a reported stream.
    const P12_SALT         = 0xFF51_AFD7_ED55_8CCD  # Murmur3 finalizer constant; not any repo salt
    const P12_DATAGEN_SALT = 0xC2B2_AE3D_27D4_EB4F  # Murmur3 finalizer constant, second word; the
                                                    # TRAINING-POOL draw stream is carved by its
                                                    # SALT (mirroring P11_DATAGEN_SALT) because the
                                                    # counter word must carry the per-sample GLOBAL
                                                    # INDEX -- that is what makes the pool
                                                    # byte-identical for any thread count.

    "The Phase-12 REPORTED RNG: `Philox4x` keyed by `(P12_DEV_SEED xor P12_SALT, counter)`."
    p12_rng(counter::Integer = 0) =
        Philox4x(UInt64, (UInt64(P12_DEV_SEED) ⊻ P12_SALT, UInt64(counter)))

    "The Phase-12 FIXTURE RNG: same construction off `P12_FIX_SEED`, a provably disjoint key."
    p12_fix_rng(counter::Integer = 0) =
        Philox4x(UInt64, (UInt64(P12_FIX_SEED) ⊻ P12_SALT, UInt64(counter)))

    # RESERVED counters, ONE PER REPORTED RUNNER, so two activities can never silently share a
    # sub-stream. Mirrors the VAL_FIX_SEED discipline of spike/validation/consts.jl:77-79 and
    # p11_consts.jl:146-159: fixtures ride their own SEED and their own counter, and NEVER
    # consume a counter a reported number rides on.
    const P12_SIM02_COUNTER      = 1   # run_p12_sim02.jl        -- per-region SIM-02 consistency
    const P12_STAGE1_COUNTER     = 2   # run_p12_stage1_ridge.jl -- the D-12 Stage-1 GATE
    const P12_ELL_RIDGE_COUNTER  = 3   # run_p12_ell_ridge.jl    -- correlation-length recovery
    const P12_EPS_RIDGE_COUNTER  = 4   # run_p12_eps_ridge.jl    -- chromatic-eps identifiability
    const P12_MINISPIKE_COUNTER  = 5   # the D-12 Stage-2 mini-spike training run
    const P12_DATAGEN_COUNTER    = 6   # p12_generate.jl         -- training-pair generation
    const P12_SBC_COUNTER        = 7   # run_p12_sbc.jl          -- rank/uniformity tables
    const P12_COVERAGE_COUNTER   = 8   # run_p12_coverage.jl     -- coverage + scoring rule
    const P12_GUARDS_COUNTER     = 9   # run_p12_guards.jl       -- the S-4 / D-06 guards
    const P12_FIXTURE_COUNTER    = 99  # FIXTURES ONLY -- never a reported counter

    # =====================================================================================
    # 3. Lattice / geometry (D-01, D-04, D-08)
    # =====================================================================================
    const P12_G      = 8    # the frozen 8x8 patch grid; 64 regions per image
    const P12_K_DEV  = 63   # D-04 HIGH-RANK arm: every non-constant DCT-II mode of an 8x8
                            # lattice (64 modes minus the constant mode, which is the global
                            # term c0). "High rank" here means "no truncation at all", so the
                            # spike arm cannot be accused of having chosen a flattering rank.
    # theta row budget: the global term, the 63 deviation coefficients, the 7 non-rho nuisances
    # (spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise, chromatic_eps)
    # and the correlation-length row, which is APPENDED LAST for the reason prior.jl:87-92 gives:
    # downstream code indexes theta POSITIONALLY, so appending keeps every existing index valid.
    const P12_D_MINISPIKE = 1 + P12_K_DEV + 7 + 1     # = 72

    # THE ARMS ARE PARAMETRIZED BY INDUCED LAG-1 CORRELATION, NOT BY alpha OR ell (R-8).
    # Measured induced lag-1 correlation runs 0.136 at alpha = 0.5 to 0.947 at alpha = 0.999, so
    # a uniform prior on alpha would put roughly 90 % of its mass on "no spatial structure" and
    # would make D-08's correlation length unidentifiable BY PARAMETRIZATION RATHER THAN BY
    # PHYSICS -- the S-1 prior-echo trap, which this milestone has already met once in Phase 11.
    # Parametrizing both arms by induced r1 is what makes the CAR-vs-GP comparison a question
    # about the kernels rather than about their coordinates.
    const P12_R1_MIN    = 0.05
    const P12_R1_MAX    = 0.95
    const P12_R1_LADDER = (0.05, 0.25, 0.50, 0.75, 0.95)   # the five reported rungs
    # D-10 / R-4: the matched ablation RETRAINS at zero induced spatial correlation. This is a
    # TRAINING arm, not a read-time pin -- pinning r1 = 0 at read time on a spatially-trained net
    # would compare a net against itself under a distribution it was never fitted to.
    const P12_ABLATION_R1 = 0.0

    const P12_GP_JITTER          = 1e-8            # diagonal jitter for the GP Cholesky
    const P12_R1_SOLVE_TOL       = 1e-6            # bisection tolerance, induced r1 -> (alpha | ell)
    const P12_CAR_ALPHA_BRACKET  = (0.0, 0.99999)  # CAR arm: proper-CAR bracket, 1.0 is the
                                                   # improper (intrinsic) limit and is excluded
    const P12_GP_ELL_BRACKET     = (0.05, 20.0)    # GP arm: correlation length in lattice units,
                                                   # from far below one cell to far above the grid

    # =====================================================================================
    # 4. THE r1 PRIOR ITSELF -- pi(theta), not a convenience alias
    # =====================================================================================
    # This is the spatial field's prior. It is Tier 1 by the same logic that puts MU_PRIOR in
    # spike/simulator/prior.jl:40: 12-07 draws from it, and 12-15's whole CAR-vs-GP comparison
    # rests on "both arms drawn from the SAME P12_R1_PRIOR" at matched induced lag-1 correlation.
    # Binding the distribution OBJECT once, here, is what stops a plan constructing its own and
    # quietly comparing two different priors -- the R-8 trap in a different guise.
    #
    # CONFIRMED 2026-07-29: uniform on r1, not a Beta and not a log-uniform. A uniform on the
    # induced quantity is the honest default given no prior belief about biological smoothness,
    # which is also D-08's stated reason for INFERRING the correlation length rather than
    # conditioning on it.
    const P12_R1_PRIOR = Uniform(P12_R1_MIN, P12_R1_MAX)

    # =====================================================================================
    # 5. SIM-02, per region (REPORTED)
    # =====================================================================================
    # IDENTICAL to prior.jl:46's SIM02_W1_TOL, fixed BEFORE the field simulator runs, and
    # deliberately NOT reverse-engineered from any measured value. The assertion is on the MAX
    # over the 64 regions, never the mean: a mean would let one badly-mismatched region hide
    # behind 63 good ones, and the per-region claim is exactly what the spatial map is selling.
    const P12_SIM02_W1_TOL_PERREGION = 0.10
    const P12_SIM02_M                = 20_000   # draws behind each per-region Wasserstein-1

    # =====================================================================================
    # 6. D-12 STAGE 1 -- THE CHEAP GATE THAT PROTECTS THE MINI-SPIKE (**GATING**)
    # =====================================================================================
    const P12_STAGE1_N               = 25_000
    const P12_STAGE1_IMSIZE          = (512, 512)
    const P12_STAGE1_TEST_FRACTION   = 0.25
    const P12_STAGE1_RATIO_CEILING   = 0.95
    const P12_STAGE1_CONTROL_CEILING = 0.5
    const P12_STAGE1_MIN_RUNGS       = 1
    const P12_STAGE1_N_TEST_PER_RUNG_DERIVED = 1250

    # DERIVATION OF P12_STAGE1_RATIO_CEILING (recorded verbatim, the way p11_consts.jl:213-241
    # records the N_min derivation). With P12_STAGE1_TEST_FRACTION of P12_STAGE1_N split across
    # the five P12_R1_LADDER rungs:
    #     n_test = 0.25 * 25 000 / 5 = 1250 DATASETS per rung
    #     sampling sd of an RMSE ratio at that n ~= 1/sqrt(2 * n_test) = 1/sqrt(2500) = 0.0200
    #     require 2.5 sd of separation from 1.0:  1 - 2.5 * 0.0200 = 0.950
    #
    # THE UNIT THE BAR IS APPLIED TO IS NOT THE UNIT IT WAS DERIVED FOR, AND THAT IS RECORDED
    # HERE BECAUSE THE DERIVATION ABOVE WILL BE QUOTED VERBATIM IN 12-STAGE1-VERDICT.md. The
    # arithmetic is the sampling sd of a SINGLE-REGION RMSE ratio, while 12-11 gates on
    # `ratio_mean`, the mean over 64 regions, whose sd is no larger. The gate is therefore
    # CONSERVATIVE, not mis-scaled -- but a reader must not take the derivation as describing the
    # gated statistic. 12-11 additionally reports `ratio_max` against the same bar so the
    # per-region reading exists.
    #
    # WHY THE COUNT IS DATASETS AND NOT REGION-DRAWS. Each test dataset contributes 64
    # region-draws, but the 64 regions of one image share every nuisance and the global term, so
    # counting region-draws as independent would overstate the power of this gate by up to 8x in
    # sd terms. The derivation deliberately counts the INDEPENDENT unit.
    #
    # THE COMPUTE CONSEQUENCE, IN THE SAME BREATH. The research sized the Stage-1 ridge at ~6 min
    # for 5 000 samples at 512-squared, so 25 000 is ~30 min -- still an order of magnitude
    # cheaper than the ~1.8 h mini-spike this gate exists to protect. That ratio IS the whole
    # point of D-12 Stage 1.

    # =====================================================================================
    # 7. D-12 STAGE 2 -- CALIBRATION *AND* A PROPER SCORING RULE (**GATING**)
    # =====================================================================================
    # UNCHANGED FROM spike/validation/consts.jl:44-46 ON PURPOSE. The settled fix for the
    # Phase-7 "M = 2000 is over-powered" finding is the TEST (a nuisance-appropriate equivalence
    # test instead of a point null on rows that are not identified), not the power.
    const P12_SBC_M    = 2000
    const P12_SBC_L    = 999
    const P12_SBC_BINS = 50

    # TARGETS ONLY -- the global term c0 and the derived Delta-rho column -- get the strict
    # point null. Those are the quantities the map is sold on; a miss there is a real defect.
    const P12_SBC_TARGET_KS_ALPHA = 0.01
    # The 7 nuisances, the correlation-length row and the high-index deviation coefficients get
    # a nuisance-appropriate EQUIVALENCE test instead. TOST inverts the over-power incentive:
    # larger M makes equivalence EASIER to establish, because the interval shrinks toward the
    # point estimate. Pitfall 8 is explicit that a point null on a non-identified row produces a
    # loud, uninformative failure; this project has already paid for that twice.
    const P12_SBC_TOST_DELTA = 0.10
    const P12_SBC_TOST_ALPHA = 0.05

    const P12_COVERAGE_NOMINAL           = 0.90
    const P12_STAGE2_COVERAGE_TOST_DELTA = 0.03
    const P12_STAGE2_N_MIN               = 271
    const Z_TWO_SIDED_90                 = 1.644854   # Phi^-1(0.95); the 90 % two-sided limb

    # DERIVATION OF P12_STAGE2_N_MIN, RECORDED VERBATIM -- AND LABELLED PRECISELY, BECAUSE TWO
    # DIFFERENT INTERVAL CONSTRUCTIONS ARE IN PLAY AND AN EARLIER DRAFT COLLAPSED THEM INTO ONE
    # WORD. The formula below is the WALD (normal-approximation) HALF-WIDTH. Require it, at
    # p-hat = P12_COVERAGE_NOMINAL, to fit inside a +/- P12_STAGE2_COVERAGE_TOST_DELTA band:
    #     z * sqrt(p(1-p)/N) = 1.644854 * sqrt(0.09)/sqrt(N) = 0.49346/sqrt(N) <= 0.03
    #     => sqrt(N) >= 16.449 => N >= 270.6 => N_min = 271
    #
    # N_MIN IS *SIZED* WITH THE WALD HALF-WIDTH AS A STANDARD APPROXIMATION; THE INTERVAL
    # REPORTED AND TESTED IS *WILSON* (spike/validation/p11_stats.jl:49-64), whose own docstring
    # records that Wald is deliberately NOT used because near p = 0.90 with n in the low hundreds
    # it covers poorly and can extend past 1. Do not call this sizing "Wilson-derived" -- Wilson
    # is a different construction and yields a different N -- and do not call the phase
    # "Wald-based" either, which would contradict p11_stats.jl. The number 271 is correct for what
    # the formula computes. Tier 1 is append-only, so the LABEL is fixed here, before it freezes.
    #
    # NOTE ON THE NAME `Z_TWO_SIDED_90`: spike/validation/p11_consts.jl:209 binds the SAME
    # unprefixed name to 1.644853627 (more digits of the same quantile). The two files are never
    # loaded into one module by the test suite, and under Julia 1.12 a same-name const rebinding
    # is legal rather than an error, so this is recorded as a fact rather than repaired -- both
    # files are frozen pre-registrations and neither may be edited.
    #
    # THE PRE-REGISTERED PASS RULE, IN WORDS: a coverage-only win with NO scoring-rule
    # improvement is NOT a pass, and a miscalibrated arm is DISQUALIFIED rather than "better".
    # P12_STAGE2_LOGSCORE_MIN is in nats per held-out region: the minimum mean improvement in
    # log predictive density of the spatial arm over the MATCHED ablation.
    const P12_STAGE2_LOGSCORE_MIN = 0.02

    # =====================================================================================
    # 8. Training augmentation (Pitfall 4)
    # =====================================================================================
    # MCAR random-region masking: k regions are masked per training sample, k drawn uniformly
    # from this set (0 to 8 of the 64 regions, i.e. 0 - 12.5 %). Its purpose is that D-09's
    # held-out-region configuration is IN-DISTRIBUTION at read time rather than a covariate
    # shift the net has never seen.
    const P12_MASK_K_SET = 0:8

    # =====================================================================================
    # 9. REPORTING-ONLY bars -- declared in advance so they are COMPARABLE, but NOT gates
    # =====================================================================================
    # S-4: a field whose radial component carries at least half its energy cannot be linearly
    # separated from a purely radial nuisance (the stage-6 chromatic term is exactly such a
    # nuisance). Reported so the confound is visible; it decides nothing.
    const P12_RADIAL_ENERGY_CEILING = 0.50
    # D-06 half-cell offset arm: relative degradation in per-region RMSE and in coverage when the
    # lattice is offset by half a cell against the patch grid.
    const P12_OFFSET_GRID_TOL = 0.10
    # The F3 boolean vacuity flag. THE COMMENT MUST NAME BOTH ESTIMATORS IT IS APPLIED TO,
    # BECAUSE THEY ARE NOT THE SAME QUANTITY:
    #   12-18's        `shrinkage`                = post_sd     / prior_sd  -- the width of a
    #                                                neural flow's posterior;
    #   12-06's/12-13's `ridge_residual_shrinkage` = residual_sd / prior_sd  -- the out-of-sample
    #                                                RMSE of a linear point predictor.
    # ONE FLOOR IS APPLIED TO BOTH DELIBERATELY: it is a bar on WHAT FRACTION OF THE PRIOR SPREAD
    # SURVIVES, which both estimate. But they coincide only under a well-specified
    # linear-Gaussian model, so the two are reported under DIFFERENT key names and read as a
    # pair, NEVER averaged. Documenting this floor against one estimator while applying it to the
    # other is this milestone's recorded unit-error class; naming both here is what prevents it.
    # The flag gates nothing.
    const P12_VACUOUS_SHRINKAGE_FLOOR = 0.90

    # =====================================================================================
    # 10. Budget
    # =====================================================================================
    const P12_N_PAIRS         = 50_000    # matches P11_N_PAIRS, so the comparison against the
                                          # Phase-11 lane is CAPACITY-CONTROLLED
    const P12_MINISPIKE_N     = 10_000    # the D-12 Stage-2 mini-spike training pool
    const P12_MINISPIKE_IMSIZE = (512, 512)
    const P12_DATAGEN_WALLCLOCK_CEILING_MIN   = 150  # EXCEEDING THIS IS A BLOCKER, RECORDED AS
    const P12_MINISPIKE_WALLCLOCK_CEILING_MIN = 150  # SUCH -- never licence to downgrade the arm
    # BINDING INVARIANT (F5): rank/coverage claims hold only under the joint the estimator was
    # TRAINED on, so train-joint must equal eval-joint. Both are READ from the frozen amended
    # pre-registration through the isolated GateV2P12 module above and are NEVER retyped here --
    # a retyped mixture is a mixture that can silently drift.
    const P12_IMSIZE_SET     = GateV2P12.SBC_IMSIZE_SET
    const P12_IMSIZE_WEIGHTS = GateV2P12.SBC_IMSIZE_WEIGHTS

    # =====================================================================================
    # 11. The iteration ledger, declared in advance
    # =====================================================================================
    # Phase 7 was amended TWICE after seeing results and Phase 11's gate was amended again; the
    # credibility cost of that is the reason this constant exists. ONE documented iteration is
    # authorised, and it is authorised HERE, before any result exists. A SECOND iteration is NOT
    # authorised by this file and cannot be authorised by amending this file. The single
    # allowance is cheapest spent BEFORE the mini-spike -- on the D-12 Stage-1 arm -- never on
    # relaxing a threshold after training.
    const P12_ITERATION_ALLOWANCE = 1

    # =====================================================================================
    # 12. TWO repo roots. CONFLATING THEM DEFEATS 12-17's DURABILITY RULE
    # =====================================================================================

    """
        p12_repo_root() -> String

    The ENCLOSING checkout: walk up from `@__DIR__` until a `.git` entry is found, testing
    `isdir(...) || isfile(...)` because IN A GIT WORKTREE `.git` IS A *FILE* carrying a gitdir
    pointer, not a directory (`spike/test/test_p13_result.jl:198-201` already handles both and is
    the precedent this mirrors).
    """
    function p12_repo_root()
        d = normpath(@__DIR__)
        while true
            g = joinpath(d, ".git")
            (isdir(g) || isfile(g)) && return d
            p = dirname(d)
            p == d && break
            d = p
        end
        return normpath(joinpath(@__DIR__, "..", ".."))
    end

    const P12_REPO_ROOT = p12_repo_root()

    """
        p12_primary_checkout() -> String

    The PRIMARY checkout -- a *different* directory from [`p12_repo_root`](@ref) whenever the
    caller runs inside a linked worktree.

    THE WALK CANNOT ANSWER THIS QUESTION, AND THE REASON IS PRECISE: IT *ACCEPTS* A WORKTREE'S
    `.git` FILE RATHER THAN *FOLLOWING* IT. `isdir(...) || isfile(...)` establishes only that a
    `.git` entry exists; the `test_p13_result.jl:198-201` precedent it comes from is answering
    "is this a checkout at all, so may I shell out to git" -- an AVAILABILITY question -- not
    "where is the durable checkout". Inside a worktree the first `.git` the walk meets is the
    worktree's pointer file, so `p12_repo_root()` returns the worktree root, and a bundle written
    there is exactly what `git worktree remove --force` deletes.

    Resolved with git rather than by parsing the pointer: `--git-common-dir` is the PRIMARY
    `.git` even when called from a linked worktree, so its parent is the primary checkout.

    There is deliberately NO `occursin("worktree", path)` substring test anywhere in this file.
    It is redundant once the root is git-derived, and it is actively harmful: inside a worktree
    `P12_REPO_ROOT` itself contains the substring, so a criterion demanding both "under the root"
    and "no `worktree` component" is unsatisfiable and would make 12-17 throw unconditionally --
    potentially after hours of training. THE GIT-DERIVED ROOT *IS* THE CHECK.
    """
    function p12_primary_checkout()
        if Sys.which("git") === nothing
            @warn "p12_primary_checkout: git is unavailable, so primary-checkout resolution is " *
                  "DEGRADED to the enclosing checkout. A bundle written here will NOT survive a " *
                  "worktree cleanup (12-17 durability rule)." enclosing = P12_REPO_ROOT
            return P12_REPO_ROOT
        end
        common = readchomp(`git -C $(P12_REPO_ROOT) rev-parse --path-format=absolute --git-common-dir`)
        return normpath(dirname(common))
    end

    const P12_PRIMARY_CHECKOUT = p12_primary_checkout()

    # WHICH CONSUMER USES WHICH, STATED HERE SO NO PLAN HAS TO GUESS:
    #   12-17 and 12-14 Task 3 use P12_PRIMARY_CHECKOUT. Trained bundles and pools are GITIGNORED
    #     and must survive a worktree cleanup, so they must land in the durable tree. Their
    #     assertion is "under P12_PRIMARY_CHECKOUT and containing no `worktree` path component".
    #   12-19 uses P12_REPO_ROOT. `test/test_images/` and `corpus/` are TRACKED, so they exist in
    #     whatever checkout is executing, and the sealed-holdout refusal should confine paths to
    #     the checkout the run is actually reading from. Using the primary root here would be
    #     wrong in the opposite direction.
    # Both are Tier 1 rather than test-local because runners depend on them for CORRECTNESS, and
    # one use is this phase's only irreversible risk: 12-19's executable refusal is what keeps the
    # Phase-16 sealed holdout sealed. 12-05 READS these bindings rather than re-deriving its own.

    # =====================================================================================
    # 13. Shared helper: the equal-tailed interval
    # =====================================================================================
    """
        _p12_interval(v, nominal) -> (lo, hi)

    The `nominal`-level equal-tailed interval of a draw vector. PORTED (three lines) from
    `spike/validation/run_p11_coverage.jl:93-101`.

    IT LIVES HERE, IN THE WAVE-1 FILE EVERY PLAN ALREADY INCLUDES, FOR A REASON WORTH RECORDING.
    The Phase-11 original cannot be reused by including its file: `run_p11_coverage.jl` has no
    `*_LOAD_ONLY` guard and ends in a bare `main()`, so an include to reach the helper would
    execute the entire Phase-11 coverage run. And the two Phase-12 consumers sit two waves apart
    -- 12-15 needs it at wave 8, while 12-16 does not create `p12_coverage.jl` until wave 10 --
    so putting it in the coverage file would leave 12-15 with a helper that does not exist yet,
    and putting a copy in each would invite exactly the drift a shared helper prevents.

    It is a pure function of its arguments, it gates nothing, and it is NOT a constant: the
    append-only Tier-2 discipline does not apply to it.
    """
    function _p12_interval(v::AbstractVector, nominal::Real)
        a = (1 - nominal) / 2
        s = sort(collect(v))
        n = length(s)
        lo = s[max(1, floor(Int, a * n))]
        hi = s[min(n, ceil(Int, (1 - a) * n))]
        return lo, hi
    end

    # =====================================================================================
    # 14. DESCOPE ROUTING, MADE EXECUTABLE
    # =====================================================================================
    # The routing existing only as prose in two documents is what makes a wrong dispatch
    # possible: six plans would each have to re-read the same paragraph and agree. ONE function
    # with four call sites is the fix.

    "The one anchored line `12-STAGE1-VERDICT.md` must carry, and the only thing read from it."
    const P12_STAGE1_VERDICT_PATTERN = r"^VERDICT: (PROCEED|DESCOPE)$"m
    const P12_STAGE1_VERDICT_FILE    = "12-STAGE1-VERDICT.md"

    """
        p12_stage1_verdict(; dir) -> :proceed | :descope | :absent

    Parse the unique anchored line `^VERDICT: (PROCEED|DESCOPE)\$` out of
    `12-STAGE1-VERDICT.md` under `dir`. Returns `:absent` when the file does not exist, and
    THROWS when the file exists but does not carry exactly one such line -- it NEVER defaults to
    either outcome, because a silently-defaulted routing decision is the failure this function
    exists to remove.

    THE `dir` KEYWORD IS REQUIRED AT DECLARATION TIME and is deliberately not deferred "until
    something needs it": 12-14's testset 7 has to exercise the absent / DESCOPE / PROCEED
    branches without writing into the real planning directory, so it must point the lookup at a
    temp dir -- and Tier 1 is APPEND-ONLY, so a signature declared without the keyword here could
    never be widened later without a pre-registration breach.

    THE DEFAULT IS `P12_REPO_ROOT`, NOT `P12_PRIMARY_CHECKOUT`, AND THAT IS THE §12 RULE APPLIED
    TO ITS OWN CASE. `12-STAGE1-VERDICT.md` is a TRACKED planning document that 12-11 writes and
    commits -- not a gitignored artefact that must survive a worktree cleanup. §12 is explicit
    that tracked paths use `P12_REPO_ROOT` "so they exist in whatever checkout is executing" and
    that the primary root "would be wrong in the opposite direction". An earlier draft of this
    signature defaulted to `P12_PRIMARY_CHECKOUT`, and inside a worktree the failure is silent
    and total: 12-11 writes and commits the verdict in the worktree, the primary checkout does
    not have it until merge, this function returns `:absent`, and `p12_require_proceed` -- the
    first statement of 12-17, 12-18, 12-19 and 12-20 -- throws on all four after a Stage-1 gate
    that legitimately PASSED, with 12-14's spend guard throwing too. `P12_REPO_ROOT` always finds
    the file the current run just wrote.
    """
    function p12_stage1_verdict(; dir = joinpath(P12_REPO_ROOT, ".planning", "phases",
                                                 "12-spatial-colocalization-map"))
        path = joinpath(dir, P12_STAGE1_VERDICT_FILE)
        isfile(path) || return :absent
        # CRLF is normalised before matching: an anchored `$` will not match across a stray \r,
        # and this repository is developed on Windows.
        src = replace(read(path, String), "\r\n" => "\n")
        ms = collect(eachmatch(P12_STAGE1_VERDICT_PATTERN, src))
        length(ms) == 1 || error(
            "p12_stage1_verdict: $path exists but carries $(length(ms)) lines matching " *
            "$(P12_STAGE1_VERDICT_PATTERN); exactly one is required. Refusing to default to " *
            "either outcome -- a routing decision must be read, never guessed.")
        return ms[1].captures[1] == "PROCEED" ? :proceed : :descope
    end

    """
        p12_require_proceed(plan_id; dir)

    Throw unless `p12_stage1_verdict(; dir)` is `:proceed`. THE FOUR DESCOPE-EXCLUDED RUNNERS --
    12-17, 12-18, 12-19 and 12-20 -- call this as their FIRST statement. 12-14 does NOT: it runs
    on both paths and applies its own narrower `arm = :none` rule instead.
    """
    function p12_require_proceed(plan_id::AbstractString;
                                 dir = joinpath(P12_REPO_ROOT, ".planning", "phases",
                                                "12-spatial-colocalization-map"))
        v = p12_stage1_verdict(; dir = dir)
        v === :proceed && return :proceed
        error("$plan_id is excluded on the DESCOPE path. 12-11's routing sentence, verbatim: " *
              "\"If the D-12 Stage-1 gate does not clear its pre-registered bars the phase " *
              "DESCOPES: the mini-spike is not trained, and 12-17, 12-18, 12-19 and 12-20 do " *
              "not run.\" " *
              "$(P12_STAGE1_VERDICT_FILE) under $dir reports :$(v) " *
              "(:absent means the Stage-1 gate has not been adjudicated yet).")
    end

    # =====================================================================================
    # 15. Declared research-net deviations
    # =====================================================================================
    # Naming the model changes IN ADVANCE is the R-5 discipline: the Phase-7 precedent is that
    # UNNAMED model changes are how credibility gets spent. The report must enumerate exactly
    # these five and no others; anything else that changes is an undeclared deviation.
    const P12_RESEARCH_NET_DEVIATIONS = (
        "a CNN over the (G,G,2,K) reshaped summary replaces the shipped MLP over 128 rows (D-02)",
        "theta widens to 1 + K + 7 + 1 rows, with DCT coefficients replacing the scalar rho_true (D-01/D-08)",
        "the lattice-prior field replaces scalar rho in forward-model stage 1 (D-05/D-06)",
        "MCAR random region masking added as a TRAINING-TIME augmentation (Pitfall 4)",
        "per-row z-scoring retained as the DEFAULT arm; shared-scalar scaling is declared as the " *
        "alternative, taken only if the CNN underperforms the MLP control, and whichever is used " *
        "MUST be recorded (R-3)",
    )

    # =====================================================================================
    # 16. GATED vs REPORTED, made machine-checkable
    # =====================================================================================
    # Two disjoint tuples of constant NAMES. A number that decides something is in the first; a
    # number that is merely comparable is in the second; nothing is in both, and the assertion at
    # the foot of this block is what makes that a property of the file rather than of a paragraph.
    const P12_GATING_CONSTANTS = (:P12_STAGE1_RATIO_CEILING, :P12_STAGE1_CONTROL_CEILING,
                                  :P12_STAGE1_MIN_RUNGS, :P12_SIM02_W1_TOL_PERREGION,
                                  :P12_SBC_TARGET_KS_ALPHA, :P12_SBC_TOST_DELTA,
                                  :P12_STAGE2_COVERAGE_TOST_DELTA, :P12_STAGE2_LOGSCORE_MIN,
                                  :P12_STAGE2_N_MIN)
    const P12_REPORTING_ONLY_CONSTANTS = (:P12_RADIAL_ENERGY_CEILING, :P12_OFFSET_GRID_TOL,
                                          :P12_VACUOUS_SHRINKAGE_FLOOR)

    # =====================================================================================
    # 17. Stage-1 golden regression (D-06) -- THREE constants, and the difference between a
    #     CRITERION and an EXCUSE is written into the comment because that is what it is
    # =====================================================================================
    # A theta carrying NO `rho_field` must reproduce the pre-edit simulator BYTE-FOR-BYTE.
    # `==` is the contract, in the voice of p11_consts.jl's P11_STAGE6_EXACT. Do NOT soften a
    # failing assertion in the consumer to an `isapprox`.
    const P12_STAGE1_EXACT = true
    # A CRITERION, NOT A FALLBACK. A constant-VALUED field still goes through the separable
    # upsample, whose row weights `1-t` and `t` sum to 1 only up to IEEE rounding, so a constant
    # field agrees with the scalar path to a few ULP and *cannot* be bit-identical. This is
    # recorded here because 12-VALIDATION.md's row reads as though it could be.
    const P12_STAGE1_CONSTFIELD_TOL = 1e-12
    # THE ESCAPE HATCH, FOR THE **EXACT** ARM ONLY, in the voice of p11_consts.jl:314-317: ANY
    # use of it is a RECORDED DEVIATION, because it means an upstream numerics change -- not a pass.
    const P12_STAGE1_GOLDEN_TOLERANCE_FALLBACK = 1e-12

    # =====================================================================================
    # EXECUTABLE self-checks (cheap; no inference, no simulation, no file writes)
    # =====================================================================================
    # R-5: the fresh streams are not any reserved, burned or derived gate seed.
    @assert !(UInt64(P12_DEV_SEED) in _p12_forbidden()) "P12_DEV_SEED collides with a forbidden seed"
    @assert !(UInt64(P12_FIX_SEED) in _p12_forbidden()) "P12_FIX_SEED collides with a forbidden seed"
    @assert UInt64(P12_DEV_SEED) != UInt64(P12_FIX_SEED) "the fixture stream must differ from the reported stream"
    # A duplicate inside the forbidden list would silently shrink the set; catch that.
    @assert length(unique(_p12_forbidden_list())) == length(_p12_forbidden_list()) "duplicate entry in the forbidden-seed list"
    # ... and the three Phase-11 burned-key tuples really ARE subsumed by Phase 13's superset,
    # which is what licenses listing the superset once instead of concatenating all four.
    @assert issubset(Set(UInt64.(F2_DEV_SEEDS)),        Set(UInt64.(P13_BURNED_DEV_SEEDS)))
    @assert issubset(Set(UInt64.(SPIKE_DEV_SEEDS)),     Set(UInt64.(P13_BURNED_DEV_SEEDS)))
    @assert issubset(Set(UInt64.(P11_RESEARCH_BURNED)), Set(UInt64.(P13_BURNED_DEV_SEEDS)))
    @assert length(P13_BURNED_DEV_SEEDS) == 25
    # Both fresh salts are absent from the repository-wide inventory.
    @assert P12_SALT ∉ P12_REPO_SALTS "P12_SALT reuses an existing repository salt"
    @assert P12_DATAGEN_SALT ∉ P12_REPO_SALTS "P12_DATAGEN_SALT reuses an existing repository salt"
    @assert P12_SALT != P12_DATAGEN_SALT
    # The recomputation reproduces the frozen file's OWN derived seeds -- proof that the
    # single-draw `_derive` above matches the gate's redraw-loop rule, so nothing is missed.
    @assert PROD_SEED_V2 == GateV2P12.PROD_SEED_V2
    @assert PROD_SEED    == GateV2P12.PROD_SEED
    # Fixtures must never ride a counter a reported number rides on.
    @assert P12_FIXTURE_COUNTER ∉ (P12_SIM02_COUNTER, P12_STAGE1_COUNTER, P12_ELL_RIDGE_COUNTER,
                                   P12_EPS_RIDGE_COUNTER, P12_MINISPIKE_COUNTER,
                                   P12_DATAGEN_COUNTER, P12_SBC_COUNTER, P12_COVERAGE_COUNTER,
                                   P12_GUARDS_COUNTER)
    # Lattice integrity, and the theta row budget the whole phase indexes positionally.
    @assert P12_D_MINISPIKE == 1 + P12_K_DEV + 7 + 1 == 72
    @assert P12_K_DEV == P12_G^2 - 1
    # Ladder integrity: sorted, and spanning exactly the declared r1 bracket.
    @assert issorted(P12_R1_LADDER)
    @assert first(P12_R1_LADDER) == P12_R1_MIN && last(P12_R1_LADDER) == P12_R1_MAX
    @assert 0.0 < P12_R1_MIN < P12_R1_MAX < 1.0
    @assert P12_ABLATION_R1 < P12_R1_MIN          # the ablation is OUTSIDE the prior, by design
    @assert first(P12_CAR_ALPHA_BRACKET) < last(P12_CAR_ALPHA_BRACKET) < 1.0
    @assert 0.0 < first(P12_GP_ELL_BRACKET) < last(P12_GP_ELL_BRACKET)
    # The prior is the bound distribution OBJECT, over exactly the declared bracket.
    @assert P12_R1_PRIOR isa Uniform
    @assert minimum(P12_R1_PRIOR) == P12_R1_MIN && maximum(P12_R1_PRIOR) == P12_R1_MAX
    # The Stage-1 arithmetic is re-derivable from the file, not only from the comment.
    @assert P12_STAGE1_N_TEST_PER_RUNG_DERIVED ==
            round(Int, P12_STAGE1_TEST_FRACTION * P12_STAGE1_N / length(P12_R1_LADDER))
    @assert round(1 - 2.5 / sqrt(2 * P12_STAGE1_N_TEST_PER_RUNG_DERIVED); digits = 3) ==
            P12_STAGE1_RATIO_CEILING
    @assert 0.0 < P12_STAGE1_CONTROL_CEILING < P12_STAGE1_RATIO_CEILING < 1.0
    @assert 1 <= P12_STAGE1_MIN_RUNGS <= length(P12_R1_LADDER)
    # The recorded N_min derivation reproduces: 0.49346/sqrt(N) <= delta => N >= 270.6 => 271.
    @assert ceil(Int, (Z_TWO_SIDED_90 * sqrt(0.09) / P12_STAGE2_COVERAGE_TOST_DELTA)^2) ==
            P12_STAGE2_N_MIN
    @assert 0.0 < P12_STAGE2_COVERAGE_TOST_DELTA < P12_COVERAGE_NOMINAL
    @assert (P12_SBC_L + 1) % P12_SBC_BINS == 0   # rank-bin evenness (Pitfall 4)
    @assert 0.0 < P12_SBC_TARGET_KS_ALPHA < P12_SBC_TOST_ALPHA
    @assert P12_STAGE2_LOGSCORE_MIN > 0.0
    # Masking is a proper subset of the 64 regions, and 0 (no masking) is in-distribution.
    @assert first(P12_MASK_K_SET) == 0 && last(P12_MASK_K_SET) < P12_G^2
    # The F5 mixture read succeeded and is a proper mixture.
    @assert length(P12_IMSIZE_SET) == length(P12_IMSIZE_WEIGHTS)
    @assert isapprox(sum(P12_IMSIZE_WEIGHTS), 1.0; atol = 1e-12)
    # The declared-deviation list is a fixed-length enumeration, not an open-ended note.
    @assert length(P12_RESEARCH_NET_DEVIATIONS) == 5
    @assert all(d -> d isa AbstractString && !isempty(d), P12_RESEARCH_NET_DEVIATIONS)
    @assert P12_ITERATION_ALLOWANCE == 1
    # REPORTED IS NOT GATED, asserted rather than asserted-in-prose.
    @assert isempty(intersect(P12_GATING_CONSTANTS, P12_REPORTING_ONLY_CONSTANTS))
    # Both roots resolved to real directories, and the primary one is where git says it is.
    @assert isdir(P12_REPO_ROOT)
    @assert isdir(P12_PRIMARY_CHECKOUT)
    @assert P12_STAGE1_EXACT === true

    # --- CLOSING NEGATIVE ASSERTIONS: none of the THREE Tier-2 blocks exists yet -----------
    # THREE SENTINELS MEANS THREE ASSERTIONS, one per BLOCK and not one per PLAN. Each is
    # removed in the SAME commit that appends its own block -- 12-15 opens :P12_CHOSEN_PRIOR and
    # :P12_N_LOW, 12-16 opens :P12_FISHERZ_NEFF -- and only then. A shared or missing assertion
    # would leave an append with a test-side negative assertion to retire but no file-side
    # counterpart, which reads as an omission to a later auditor. The mirror of these three lives
    # in spike/test/test_p12_consts.jl testset 8, so the file side and the test side stay
    # symmetric.
    @assert !isdefined(@__MODULE__, :P12_CHOSEN_PRIOR)
    @assert !isdefined(@__MODULE__, :P12_N_LOW)
    # :P12_FISHERZ_NEFF's negative assertion was REMOVED 2026-07-31 by 12-16, in the SAME COMMIT
    # that appended its Tier-2 block (block 5, at the foot of this file) and its mirror in
    # spike/test/test_p12_consts.jl testset "Tier 2 does not exist yet". Dropping only one of the
    # three -- file-side assertion, test-side assertion, block -- is the documented way this goes
    # wrong: it leaves an append with a retired counterpart on one side only, which reads as an
    # omission to a later auditor. The other two sentinels stay closed: 12-15's NONE-BEATS-ABLATION
    # branch appended NEITHER :P12_CHOSEN_PRIOR nor :P12_N_LOW, and this plan opens neither.
end

# =========================================================================================
# TIER-2 APPEND (block 4 of 4) -- the D-12 STAGE-1 CONTROL ADJUDICATION
# Sentinel: :P12_STAGE1_CONTROL_ADJUDICATION.  Opened 2026-07-31 BY USER RULING.
# Provenance artifact: spike/validation/p12_stage1_report.jld2 (generated 2026-07-29T19:14:32.274Z,
# runner committed at e88b97f BEFORE any pool existed).  Full reasoning:
# .planning/phases/12-spatial-colocalization-map/12-STAGE1-VERDICT.md.
# =========================================================================================
#
# WHY A FOURTH BLOCK EXISTS WHEN THE HEADER AT :33-47 RESERVES THREE. Three were reserved at
# FREEZE TIME, when the only foreseen Tier-2 appends were measurements two later plans would
# need. This one was not foreseen: it records the adjudication of a Tier-1 GATE COMPONENT that
# the run showed to be mis-scaled. The header block above is deliberately left BYTE-UNCHANGED
# rather than amended to say "four" -- editing the pre-registration's own description of itself
# to accommodate a later event is precisely the move this file exists to make impossible. The
# header is therefore correct AS A RECORD OF WHAT WAS FORESEEN, and this comment is the record
# of what was not.
#
# WHAT WAS ADJUDICATED, in one paragraph, so this file is readable without the verdict document.
# The D-12 Stage-1 gate is an AND of two components. `borrowing_ok` -- the half the phase is
# about -- passed at 5 of 5 rungs (`ratio_mean` 0.926/0.847/0.719/0.562/0.336 against
# P12_STAGE1_RATIO_CEILING = 0.95). `control_live` FAILED (max control ratio 0.732 against
# P12_STAGE1_CONTROL_CEILING = 0.5), and the frozen rule reads a failed control as "the harness
# is not live". THE RULING IS THAT THAT PREMISE IS FALSE HERE, on two measured facts:
#   (a) LIVENESS IS PROVEN, by a certificate the gate never asked for. The own-row ridge attains
#       the analytic information limit sqrt(1 - corr^2) to within 5e-4 at EVERY rung -- a bound
#       computed with no ridge, no split and no standardizer, i.e. without any part of the
#       machinery the control exists to audit, so no defect in that machinery could manufacture
#       agreement with it. An optimal estimator is not a dead one.
#   (b) THE CONTROL, AS SPECIFIED, IS NOT A LIVENESS TEST. The own-row limit is ABOVE 0.5 at four
#       of the five rungs, so nothing could reach 0.5 from a region's own row there; the only
#       thing that carries the control under 0.5 at all is BORROWING from the other 63 regions --
#       and borrowing is the very effect the gate exists to measure. The control therefore passes
#       exactly where `borrowing_ok` is strongest and fails where it is weakest. It is a second,
#       harder borrowing test wearing a positive control's name, and it cannot certify the
#       instrument independently of the effect under test.
# The ceiling is mis-scaled for the unit it is applied to; the harness is provably live by (a);
# the gate reduces to `borrowing_ok`; the verdict is PROCEED.
#
# THE CEILING KEEPS ITS VALUE. `P12_STAGE1_CONTROL_CEILING = 0.5` at :354 is NOT edited and NOT
# superseded by a friendlier number -- it stands as the historical record of what was
# pre-registered, exactly as `SPEEDUP_GATE = 100.0` and `P11_LAMBDA_ABLATION_FACTOR = 2.502` kept
# theirs after their own bars were found mis-specified. Nothing below is a REPLACEMENT BAR.
#
# READ THIS BEFORE USING ANY CONSTANT BELOW: EVERY ONE IS A MEASUREMENT OR A RECORD. NONE IS A
# THRESHOLD. None enters `P12_GATING_CONSTANTS` (which is Tier-1 and append-only, so it could not
# be extended even if that were wanted) and none may be applied as a pass/fail bar by any later
# plan. The Stage-1 gate was adjudicated ONCE, by a human, on the numbers below; re-applying a
# post-hoc criterion to the same data would be the data-snooping this whole file forbids.
#
# WHY THE 0.5 CEILING WAS REACHABLE-LOOKING AND THE 0.95 ONE WAS NOT MIS-SCALED: the ratio
# ceiling carries its derivation verbatim at :358 AND is re-derived by an executable assertion at
# :763-764, so a wrong value could not have loaded. The control ceiling has NO derivation block
# anywhere, and its only guard at :765 is the ordering constraint 0 < CONTROL < RATIO < 1, which
# 0.4, 0.5 and 0.6 satisfy identically. That asymmetry is the defect, and it is why the two
# ceilings were adjudicated differently.
if !isdefined(@__MODULE__, :P12_STAGE1_CONTROL_ADJUDICATION)

    # --- The ruling itself, as a readable symbol ------------------------------------------
    const P12_STAGE1_CONTROL_ADJUDICATION = :ceiling_mis_scaled_control_at_information_limit
    const P12_STAGE1_CONTROL_ADJUDICATION_DATE     = "2026-07-31"
    const P12_STAGE1_CONTROL_ADJUDICATION_BY       = :user
    const P12_STAGE1_CONTROL_ADJUDICATION_ARTIFACT = "spike/validation/p12_stage1_report.jld2"
    const P12_STAGE1_CONTROL_ADJUDICATION_DOC      =
        ".planning/phases/12-spatial-colocalization-map/12-STAGE1-VERDICT.md"

    # --- The measurements the ruling rests on, one entry per P12_R1_LADDER rung ------------
    # All four vectors are copied at FULL Float64 precision from the artifact keys of the same
    # name, read back out of the .jld2 -- not retyped from a console log or a summary document.

    # `own_row_bound_mean`: sqrt(1 - corr(row_r, z_field[r])^2), the RMSE-over-prior ratio of the
    # BEST POSSIBLE linear predictor of region r's lattice value from region r's OWN summary row.
    # A property of the DATA. No ridge, no penalty, no split, no standardizer -- which is exactly
    # why it can audit the machinery: no defect in that machinery could manufacture agreement
    # with a bound that does not use it.
    const P12_STAGE1_OWNROW_LIMIT_MEASURED = (0.765640750897169, 0.7120078394184097,
                                              0.6366720262496517, 0.5555131021382556,
                                              0.4718299883347116)
    # `ownrow_ratio_mean`: the SAME quantity as estimated by the ridge under audit.
    const P12_STAGE1_OWNROW_RIDGE_MEASURED = (0.7658810882318206, 0.7123565992857308,
                                              0.6371620688064296, 0.5558956237524629,
                                              0.4719855253907083)
    # `control_ratio_mean`: the full-128-row positive control -- region r's own row INCLUDED.
    # This is the vector `control_live` was computed from, and its maximum (0.73188 at r1 = 0.05)
    # is the number that failed the 0.5 ceiling.
    const P12_STAGE1_CONTROL_RATIO_MEASURED = (0.7318805770440442, 0.6774010699701734,
                                               0.5982727231625359, 0.49658948181593837,
                                               0.3236859998843231)
    # `global_control_ratio`: the Phase-11-COMPARABLE benchmark -- recovery of the IMAGE-LEVEL
    # (not per-region) field mean. 12-CONTEXT S-2 records Phase 11 recovering global rho_true at
    # 0.157 on a harness known to work; the comparable rung here is r1 = 0.95, where the global
    # level of a near-constant field is the analogue of a constant rho, and it measures 0.19687
    # at the noisiest single image size in the F5 set. THIS IS THE ONLY LEGITIMATE COMPARISON TO
    # PHASE 11'S 0.157 -- S-2 warns explicitly that 0.157 must NOT be read as a per-region
    # expectation, and reading it as one is how the 0.5 ceiling most plausibly acquired its value.
    const P12_STAGE1_GLOBAL_CONTROL_MEASURED = (0.7717256500340796, 0.6140011119059539,
                                                0.38141373438720366, 0.23233993788672802,
                                                0.1968675170141333)

    # The rungs at which P12_STAGE1_CONTROL_CEILING lies BELOW the OWN-ROW information limit --
    # i.e. where no estimator of any kind, correct or broken or neural, could reach 0.5 FROM A
    # REGION'S OWN ROW. FOUR of five, not three. Both DERIVED from the constants above by the
    # assertion block below, not counted by hand.
    #
    # THE QUALIFIER "FROM A REGION'S OWN ROW" IS LOAD-BEARING AND MUST NOT BE DROPPED. The
    # own-row bound is NOT an upper bound on the full-128 control, which sees all 128 rows and
    # legitimately BEATS it at every rung by borrowing from the other 63 regions. A condensed
    # retelling of this ruling stated "the limit is above the bar at the three shortest rungs",
    # which is true but incomplete -- the limit (0.5555) is above the bar at r1 = 0.75 as well,
    # where the control nonetheless CLEARED the bar (0.4966) precisely because it borrows. Both
    # tuples are recorded so the two facts can never again be collapsed into one sentence.
    const P12_STAGE1_OWNROW_LIMIT_ABOVE_CEILING_RUNGS = (0.05, 0.25, 0.50, 0.75)
    const P12_STAGE1_CONTROL_CLEARS_CEILING_RUNGS     = (0.75, 0.95)

    # --- EXECUTABLE self-checks: the ruling's premises, re-derived from this file -----------
    # These do not gate anything. They make the ARGUMENT falsifiable at include time: if a later
    # edit ever made one of these false, the reasoning in 12-STAGE1-VERDICT.md would no longer
    # follow from the numbers, and the file would refuse to load rather than quietly disagree
    # with its own verdict.
    @assert length(P12_STAGE1_OWNROW_LIMIT_MEASURED)   == length(P12_R1_LADDER)
    @assert length(P12_STAGE1_OWNROW_RIDGE_MEASURED)   == length(P12_R1_LADDER)
    @assert length(P12_STAGE1_CONTROL_RATIO_MEASURED)  == length(P12_R1_LADDER)
    @assert length(P12_STAGE1_GLOBAL_CONTROL_MEASURED) == length(P12_R1_LADDER)
    # PREMISE 1 -- the estimator ATTAINS the information limit at every rung (optimal, not
    # broken). The published claim is "to within 5e-4"; assert exactly that.
    @assert all(P12_STAGE1_OWNROW_RIDGE_MEASURED .>= P12_STAGE1_OWNROW_LIMIT_MEASURED)
    @assert maximum(P12_STAGE1_OWNROW_RIDGE_MEASURED .- P12_STAGE1_OWNROW_LIMIT_MEASURED) < 5e-4
    # PREMISE 2 -- the ceiling sits BELOW the own-row information limit at four of the five
    # rungs, so at those four no own-row estimator could have cleared it...
    @assert Tuple(P12_R1_LADDER[collect(P12_STAGE1_OWNROW_LIMIT_MEASURED) .>
                                P12_STAGE1_CONTROL_CEILING]) ==
            P12_STAGE1_OWNROW_LIMIT_ABOVE_CEILING_RUNGS
    @assert length(P12_STAGE1_OWNROW_LIMIT_ABOVE_CEILING_RUNGS) == 4
    # ... and the full-128 control cleared it at exactly the two LONGEST rungs, which is the
    # substance of the defect: the only thing that carries the control under 0.5 is BORROWING,
    # and borrowing is the effect the gate exists to measure. A positive control that can only
    # pass when the effect under test is present is not an independent liveness certificate.
    @assert Tuple(P12_R1_LADDER[collect(P12_STAGE1_CONTROL_RATIO_MEASURED) .<=
                                P12_STAGE1_CONTROL_CEILING]) ==
            P12_STAGE1_CONTROL_CLEARS_CEILING_RUNGS
    # PREMISE 3 -- borrowing genuinely helps: the full-128 control beats the own-row limit at
    # EVERY rung, which is what carries it under the ceiling at the two rungs where that is
    # physically possible at all.
    @assert all(P12_STAGE1_CONTROL_RATIO_MEASURED .< P12_STAGE1_OWNROW_LIMIT_MEASURED)
    # PREMISE 4 -- the number that failed really is the maximum of the control vector, so the
    # gate arithmetic recorded in the artifact is reproduced here rather than asserted in prose.
    @assert maximum(P12_STAGE1_CONTROL_RATIO_MEASURED) > P12_STAGE1_CONTROL_CEILING
    # PREMISE 5 -- the Phase-11 benchmark reproduces where the two are comparable (r1 = 0.95),
    # and is monotone in r1 for the reason the physics predicts.
    @assert last(P12_STAGE1_GLOBAL_CONTROL_MEASURED) < 0.25
    @assert issorted(P12_STAGE1_GLOBAL_CONTROL_MEASURED; rev = true)
    # THE TIER-1 RECORD IS INTACT: the ceiling keeps its pre-registered value, and the ruling
    # spent no iteration. Both are asserted HERE, in the block that adjudicated them, so that a
    # later edit to either would fail at the point where its meaning is written down.
    @assert P12_STAGE1_CONTROL_CEILING == 0.5
    @assert P12_ITERATION_ALLOWANCE    == 1
    # The adjudication is a RECORD, never a bar: none of its names may join the gating tuple.
    @assert :P12_STAGE1_CONTROL_ADJUDICATION ∉ P12_GATING_CONSTANTS
    @assert :P12_STAGE1_OWNROW_LIMIT_MEASURED ∉ P12_GATING_CONSTANTS
    @assert :P12_STAGE1_CONTROL_RATIO_MEASURED ∉ P12_GATING_CONSTANTS
    # The verdict this adjudication produced, recorded as a SYMBOL rather than checked by
    # calling `p12_stage1_verdict()` here. That call was written first and then deliberately
    # removed: it would make this constants file -- which every Phase-12 runner and test
    # includes -- fail to LOAD whenever the planning document is absent, e.g. in a checkout with
    # a partial .planning tree. The dependency must run consts -> (nothing), and routing ->
    # document; inverting it would let a docs-only change break the spike. The file-side and
    # document-side records are instead cross-checked by spike/test/test_p12_consts.jl, which is
    # the right place for a check that needs the working tree.
    const P12_STAGE1_VERDICT_ADJUDICATED = :proceed
    @assert P12_STAGE1_VERDICT_ADJUDICATED === :proceed
end

# =========================================================================================
# TIER-2 APPEND (block 5) -- THE FISHER-z OBSERVATION-NOISE MODEL
# Sentinel: :P12_FISHERZ_NEFF.  Opened by 12-16, which is the plan the header at :47 reserved
# it for -- unlike block 4, this one WAS foreseen at freeze time.
# Provenance artifact: spike/validation/p12_coverage_sim_report.jld2
#   generated     = 2026-07-31T21:46:02.339Z
#   elapsed_min   = 19.988
#   produced by   spike/validation/run_p12_coverage.jl on the RESERVED stream
#                 p12_rng(P12_COVERAGE_COUNTER = 8) off P12_DEV_SEED, in ONE reported run
#   scored bundle spike/artifacts/p12/p12_ablation_none.jld2
#                 sha256 4bdf5b030ae0a02eb7052479a53dc20bedbac5fb9d432986f26fae9112b6f017
# =========================================================================================
#
# THIS IS A DECLARED MODELLING ASSUMPTION, NOT A THRESHOLD, AND NOTHING BELOW GATES ANYTHING.
# `n_eff` parametrizes the per-region observation-noise model whose JOINT with the posterior D-09
# validates. No branch anywhere compares a measured quantity against it, none of these names joins
# `P12_GATING_CONSTANTS` (which is Tier 1 and append-only, so it could not be extended even if that
# were wanted), and a later plan may NOT apply any of them as a pass/fail bar.
#
# WHAT WAS MEASURED. Under the Fisher-z model the per-region correlation entry satisfies
# `Var(atanh(r)) = 1/(n_eff - 3)`. `calibrate_neff` draws `n_theta = 5` prior parameters at each
# image size and observes EACH of them `n_obs = 20` times through independent channel-pair
# simulations, so within a block the parameter is FIXED and only the observation varies -- the
# spread is observation noise and nothing else. A design that redrew the field each time would
# measure the prior convolved with the noise and would fit an `n_eff` far too small.
#
# THE PER-SIZE VALUES ARE THE ONES THE RUN APPLIED. `n_eff` moves by more than 12x across
# `P12_IMSIZE_SET` -- 36.8 at 512^2 against 454.9 at 2048^2 -- exactly the image-size dependence
# Pitfall 2 measured for the per-region noise sd, and the reason `12-16-PLAN.md:168` asks for a
# per-size calibration rather than one number. THE IMAGE SIZE IS OBSERVED AT READ TIME, NOT LATENT,
# so conditioning on it is free and is also available on a real image. Applying the pooled constant
# instead would make every 512^2 interval too narrow and every 2048^2 interval too wide, and the
# pooled coverage could then sit at nominal as the AVERAGE OF TWO OPPOSITE MISCALIBRATIONS while
# neither size was calibrated. `P12_FISHERZ_NEFF.applied` records which was used.
#
# THE POOLED SCALAR IS RECORDED BUT WAS NOT APPLIED, and `P12_FISHERZ_NEFF_VARZ_FIT_RESIDUAL` is
# the evidence for why: it is the largest gap between a size's REALIZED Var(z) and the variance the
# pooled constant implies, and at 0.0181 it is LARGER THAN THE POOLED VARIANCE ITSELF (0.01146).
# A single constant does not fit these four sizes, and that fact is recorded as a measurement
# rather than left for a reader to rediscover.
if !isdefined(@__MODULE__, :P12_FISHERZ_NEFF)

    # ONE BINDING CARRYING BOTH, SO THE POOLED VALUE CANNOT BE REACHED BY ACCIDENT. A bare
    # `P12_FISHERZ_NEFF = 90.2` beside a separate per-size Dict would let a later caller pick the
    # scalar because it is shorter to type, which is precisely the choice this measurement says is
    # wrong. Reading `.pooled` is then a deliberate act, and `.applied` says what the run did.
    const P12_FISHERZ_NEFF = (
        by_imsize = Dict{String,Float64}(
            "(512, 512)"    => 36.84437519087758,
            "(1024, 1024)"  => 128.307574931308,
            "(1376, 1028)"  => 166.42790882779167,
            "(2048, 2048)"  => 454.9270535543171,
        ),
        pooled  = 90.22383775087685,
        applied = :per_imsize,
    )

    # The realized Var(z) each entry above was inverted from, at FULL Float64 precision, read back
    # out of the artifact rather than retyped from a console log.
    const P12_FISHERZ_VARZ_BY_IMSIZE = Dict{String,Float64}(
        "(512, 512)"    => 0.02954700727551148,
        "(1024, 1024)"  => 0.0079803635219035,
        "(1376, 1028)"  => 0.006118905927222789,
        "(2048, 2048)"  => 0.0022127464867066427,
    )
    # REPORTED, GATES NOTHING. See the header: larger than the pooled variance itself.
    const P12_FISHERZ_NEFF_VARZ_FIT_RESIDUAL = 0.01808225147267538
    # The calibration design, recorded so the conditioning is readable from the file.
    const P12_FISHERZ_NEFF_N_THETA = 5
    const P12_FISHERZ_NEFF_N_OBS   = 20
    const P12_FISHERZ_NEFF_ARTIFACT  = "spike/validation/p12_coverage_sim_report.jld2"
    const P12_FISHERZ_NEFF_GENERATED = "2026-07-31T21:46:02.339Z"

    # --- EXECUTABLE self-checks: the block re-derives its own arithmetic ---------------------
    # These gate nothing. They make the block FALSIFIABLE at include time, the same discipline
    # block 4 applies to its own premises.

    # (1) Every entry admits a positive Fisher-z variance. `n_eff <= 3` would make 1/(n_eff-3)
    # non-positive and every predictive interval undefined.
    @assert all(v -> v > 3.0, values(P12_FISHERZ_NEFF.by_imsize))
    @assert P12_FISHERZ_NEFF.pooled > 3.0
    # (2) The keys are EXACTLY the frozen F5 mixture -- train-joint == eval-joint. A missing size
    # would silently fall back at read time; an extra one would be a size nothing was trained on.
    @assert Set(keys(P12_FISHERZ_NEFF.by_imsize)) == Set(string.(P12_IMSIZE_SET))
    @assert Set(keys(P12_FISHERZ_VARZ_BY_IMSIZE)) == Set(string.(P12_IMSIZE_SET))
    # (3) Each n_eff really is the inverse of its recorded variance, so the two records cannot
    # drift apart into a pair that looks consistent and is not.
    for k in keys(P12_FISHERZ_VARZ_BY_IMSIZE)
        @assert isapprox(P12_FISHERZ_NEFF.by_imsize[k],
                         3 + 1 / P12_FISHERZ_VARZ_BY_IMSIZE[k]; rtol = 1e-9)
    end
    # (4) THE POOLED VALUE IS RE-DERIVED FROM THE PER-SIZE ONES, not trusted. It is fitted on the
    # MEAN VARIANCE, never averaged over the per-size n_eff -- n_eff is a nonlinear function of the
    # variance, so a mean of n_eff values is not the n_eff of the mean variance, and recording the
    # wrong one would misstate what a reader reaching for `.pooled` would get.
    @assert isapprox(P12_FISHERZ_NEFF.pooled,
                     3 + 1 / (sum(values(P12_FISHERZ_VARZ_BY_IMSIZE)) /
                              length(P12_FISHERZ_VARZ_BY_IMSIZE)); rtol = 1e-9)
    # (5) THE SIZE DEPENDENCE IS REAL, AND THAT IS WHY THE PER-SIZE CONSTANT WAS APPLIED. The
    # largest n_eff exceeds the smallest by more than 10x, and the pooled misfit exceeds the pooled
    # variance itself -- both asserted rather than described.
    @assert maximum(values(P12_FISHERZ_NEFF.by_imsize)) >
            10 * minimum(values(P12_FISHERZ_NEFF.by_imsize))
    @assert P12_FISHERZ_NEFF_VARZ_FIT_RESIDUAL > 1 / (P12_FISHERZ_NEFF.pooled - 3)
    @assert P12_FISHERZ_NEFF.applied === :per_imsize
    # (6) A MEASUREMENT, NEVER A BAR. None of these names may join the gating tuple.
    @assert :P12_FISHERZ_NEFF ∉ P12_GATING_CONSTANTS
    @assert :P12_FISHERZ_VARZ_BY_IMSIZE ∉ P12_GATING_CONSTANTS
    @assert :P12_FISHERZ_NEFF_VARZ_FIT_RESIDUAL ∉ P12_GATING_CONSTANTS
    # (7) The Tier-1 record is intact and the ruling spent no allowance, asserted HERE in the block
    # that appended, so a later edit fails where its meaning is written down.
    @assert P12_COVERAGE_NOMINAL == 0.90
    @assert P12_STAGE2_N_MIN == 271
    @assert P12_ITERATION_ALLOWANCE == 1
end
