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

# spike/p13/consts.jl --- Phase-13 Tier-1 pre-registration (D-04).
#
# THE ANTI-SNOOPING CONTRACT (D-04). EVERY THRESHOLD, LADDER, SEED, STATISTIC AND
# DESIGN CHOICE THAT ANY REPORTED PHASE-13 NUMBER WILL BE SCORED AGAINST IS LOCKED
# IN THIS ONE FILE AND COMMITTED BEFORE ANY PHASE-13 CODE RUNS AND BEFORE A SINGLE
# LABELLED DATUM EXISTS. Phase 13 promotes "mutually exclusive" to a first-class
# hypothesis on the back of a network that does not yet exist; the only thing that
# makes its eventual numbers citable is that the bars were set FIRST. Changing a
# value below after a run would be "tune until it passes" data-snooping. These
# constants ARE the committed pre-registration.
#
# TWO-TIER STRUCTURE -- READ THIS BEFORE EDITING ANYTHING.
#   TIER 1 is the single guard block in THIS file, committed before anything runs.
#   TIER 2 is `P13_TAU` -- the ONE number D-06 says must be MEASURED, not chosen. It
#     is DELIBERATELY ABSENT from this file. Plan 13-10 runs the tau probe against
#     the POST-Phase-11 simulator and APPENDS the measured value in a SECOND guard
#     block keyed on a DIFFERENT sentinel const, carrying a one-line provenance
#     comment naming the probe artifact it came from.
#   NO CONSTANT IN EITHER BLOCK IS EVER *EDITED*. Constants are only ever APPENDED,
#   so the git history of this file is itself the pre-registration audit trail. A
#   diff that MODIFIES a line below is, by construction, a pre-registration breach.
#   The Tier-1 self-check at the foot of this block asserts `P13_TAU` is ABSENT; the
#   executor of plan 13-10 must RELOCATE that assertion (never delete it) so the
#   Tier-2 block lands below it.
#
# PROVENANCE DISCLOSURE (the gate_consts_8_v2.jl:15-26 voice). Indicative, throw-away
# measurements were taken while researching this phase -- class masses under pi, the
# truncated-Cauchy atom probabilities, the real-fixture ingestion numbers and the
# measured alpha ladders. Some magnitudes below were therefore chosen with an
# approximate sense of scale already in hand. That is precisely why every value here
# is argued ONLY from properties of the DESIGN -- AUC standard-error arithmetic, the
# empty-bin behaviour of the shared binning routine, the resolution the fixed summary
# has, conventional traffic-light bands, the boost/complement algebra -- and NEVER from
# what any net would pass. The keys those runs consumed are listed as BURNED below and
# are forbidden here. NOTHING REPORTED HAS BEEN RUN AGAINST THIS FILE. It is a
# specification, frozen before use.
#
# SEED DISCIPLINE (D-01). P13_DEV_SEED and P13_FIX_SEED are FRESH Random123 streams,
# EXECUTABLY proven disjoint from every reserved stream in the repository: the spike
# training stream (NPE_MASTER_SEED), the spike validation stream (VAL_MASTER_SEED), the
# spike fixture stream (VAL_FIX_SEED), the productionization datagen stream
# (DEFAULT_MASTER_SEED), the shipped ratio-pairing stream (RATIO_PAIR_SEED), the
# provenance-manifest stream (CORPUS_MASTER_SEED), the Phase-11 DEV stream
# (P11_DEV_SEED -- forbidden EVEN THOUGH PHASE 11 HAS NOT RUN YET), every prior burned
# spike DEV seed, and all eight seeds of BOTH frozen gate families (PROD_SEED,
# PROD_SEED_V2) -- which are DERIVED, not literal, and are therefore RECOMPUTED here
# rather than trusted from a comment.
#
# EXECUTION BLOCKS ON PHASE 11 (D-02). Phase 13 trains on the registration-aware
# `zt` of Phase 11's research NPE and conditions on Phase 11's lambda input (D-03).
# Phase 13 must NEVER fall back to the shipped grid-8 basis: that is a different input
# surface and would invalidate D-02, D-03 and the D-12 continuity framing. The ONE
# sanctioned read of the shipped grid-8 bundle is as a frozen OOD COMPARISON REFERENCE
# (P13_REAL_OOD_COMPARISON), never as an evidence basis.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local constants; touches no src/,
# adds no dependency, and reaches test/gate/ only READ-ONLY through an isolated module.
#
# Guarded as ONE Tier-1 block keyed on :P13_DEV_SEED so a re-include under a test file
# or a runner is a silent no-op -- a redefinition to the same value would otherwise
# warn on a `const`.

# The frozen AMENDED grid-8 ship-gate pre-registration, loaded into an ISOLATED module
# purely to READ its derived seed families (PROD_SEED / PROD_SEED_V2) and its F5 image-size
# mixture. NO GATE IS RUN. THE FILE IS NOT MODIFIED. (spike-014 precedent,
# .planning/spikes/014-bf-sim-validation/run_sim_bf.jl:40-45; house twin at
# spike/validation/p11_consts.jl:66-77.)
#
# THIS SITS OUTSIDE THE GUARD BLOCK ON PURPOSE: Julia rejects a `module` expression that
# is not at top level ("syntax: \"module\" expression not at top level"), so it cannot
# live inside the guard block below. Idempotency is instead provided by the guarded-include
# idiom every caller uses:
#     isdefined(@__MODULE__, :P13_DEV_SEED) || include(".../p13/consts.jl")
#
# NOTE THE PATH DEPTH: from spike/p13/ the repo root is TWO levels up, not three.
module _GC
    include(joinpath(@__DIR__, "..", "..", "test", "gate", "gate_consts_8_v2.jl"))
end

if !isdefined(@__MODULE__, :P13_DEV_SEED)
    import Random123: Philox4x     # counter-based RNG; the fresh disjoint Phase-13 streams

    # =====================================================================================
    # A. FORBIDDEN seeds and the FRESH Phase-13 streams (D-01)
    # =====================================================================================
    # Literals below were read from their defining sites, not from any summary document.
    const DEFAULT_MASTER_SEED = 0x0000_0000_0000_0001  # productionization DATAGEN (datagen.jl:52)
    const NPE_MASTER_SEED     = 0x0000_0000_00C0_FFEE  # spike TRAINING stream (npe/train_npe.jl:65)
    const VAL_MASTER_SEED     = 0x0000_0000_5BC0_FFEE  # spike VALIDATION stream (validation/consts.jl:76)
    const VAL_FIX_SEED        = 0x0000_0000_00F1_F7ED  # spike FIXTURE stream (validation/consts.jl:79 = 0xF1F7ED)
    const RATIO_PAIR_SEED     = 0x0000_0000_004A_7107  # shipped ratio pairing (src/amortized/train_ratio.jl:60)
    const CORPUS_MASTER_SEED  = 0x0000_0000_00C0_5EED  # provenance-manifest stream (manifest.csv header)
    # FORBIDDEN EVEN THOUGH PHASE 11 HAS NOT RUN YET (13-RESEARCH K1): Phase 11's Tier-1
    # pre-registration already claims this stream, and Phase 13 executes AFTER it (D-02).
    # Reserving it now is what stops two research lanes silently sharing a key.
    const P11_DEV_SEED        = 0x0000_0000_0B11_DE71  # spike/validation/p11_consts.jl:144

    # Every DEV key any prior spike, diagnostic or research probe has already OBSERVED. A
    # pre-observed stream can never govern a reported number (D-01), so each is burned.
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

    # The salt inventory of the whole repository. A fresh salt must not be any of these,
    # because a repeated (seed, salt) pair is a repeated Philox key, hence a repeated stream.
    const P13_REPO_SALTS = (
        0x9E37_79B9_7F4A_7C15,   # HOLDOUT_SALT      (spike/data/seeding.jl:45)
        0xD1B5_4A32_D192_ED03,   # FOLD / CBS_SPLIT  (spike/data/seeding.jl:46)
        0xBF58_476D_1CE4_E5B9,   # VAL_SALT
        0x94D0_49BB_1331_11EB,   # PROD_SALT         (gate_consts_8_v2.jl:266)
        0xC4CE_B9FE_1A85_EC53,   # AMEND_SALT        (gate_consts_8_v2.jl:295)
        0xA24B_AED4_663E_E121,   # P11_SALT          (p11_consts.jl:143)
        0xA5A5_A5A5_A5A5_A5A5,   # COSTES
        0xA5A5_A5A5_DEAD_BEEF,   # COSTES, second
    )

    "The forbidden-seed list in declaration order, so an accidental duplicate stays visible."
    _p13_forbidden_list() = vcat(
        UInt64[UInt64(0),
               UInt64(DEFAULT_MASTER_SEED), UInt64(NPE_MASTER_SEED),
               UInt64(VAL_MASTER_SEED),     UInt64(VAL_FIX_SEED),
               UInt64(RATIO_PAIR_SEED),     UInt64(CORPUS_MASTER_SEED),
               UInt64(P11_DEV_SEED)],
        collect(UInt64, P13_BURNED_DEV_SEEDS),
        # DERIVED, NEVER TRUSTED FROM A COMMENT: both frozen gate families are RECOMPUTED by
        # the isolated read-only `_GC` module above and forbidden from the recomputed values.
        collect(UInt64, values(_GC.PROD_SEED)),
        collect(UInt64, values(_GC.PROD_SEED_V2)),
    )

    """
        _p13_forbidden() -> Set{UInt64}

    Every seed a Phase-13 stream must not be: zero, the seven named reserved streams, the
    Phase-11 DEV stream, every prior burned spike key, and all eight RECOMPUTED ship-gate
    seeds. `P13_DEV_SEED` and `P13_FIX_SEED` are asserted absent from this set at the foot
    of this block AND, independently, in `spike/test/test_p13_consts.jl`.
    """
    _p13_forbidden() = Set{UInt64}(_p13_forbidden_list())

    # --- The FRESH Phase-13 streams -------------------------------------------------------
    # A distinct (seed, salt) pair gives a different Philox key, hence a disjoint stream; the
    # assertions below make the disjointness EXPLICIT rather than merely overwhelmingly likely.
    const P13_DEV_SEED = 0x0000_0000_0B13_DE71  # "0B13" = phase 13, "DE71" = DEV-1 (mirrors P11)
    const P13_FIX_SEED = 0x0000_0000_0B13_F1F7  # FIXTURE stream, mirroring the VAL_FIX_SEED
                                                # discipline: fixtures must NEVER consume (and so
                                                # never pre-observe) a reported stream.
    const P13_SALT     = 0x2545_F491_4F6C_DD1D  # xorshift64* multiplier; not any repo salt

    "The Phase-13 REPORTED RNG: `Philox4x` keyed by `(P13_DEV_SEED xor P13_SALT, counter)`."
    p13_rng(counter::Integer = 0) =
        Philox4x(UInt64, (UInt64(P13_DEV_SEED) ⊻ P13_SALT, UInt64(counter)))

    "The Phase-13 FIXTURE RNG: same construction off `P13_FIX_SEED`, a provably disjoint key."
    p13_fix_rng(counter::Integer = 0) =
        Philox4x(UInt64, (UInt64(P13_FIX_SEED) ⊻ P13_SALT, UInt64(counter)))

    # RESERVED counters, one per activity, so two activities can never silently share a
    # sub-stream. The fixture counter rides the FIXTURE seed and is never a reported counter.
    const P13_TAU_COUNTER        = 1   # D-06 tau probe (simulator-only)
    const P13_DATAGEN_COUNTER    = 2   # labelled three-way training-pair generation
    const P13_GATE_COUNTER       = 3   # D-12 / D-13 reported gate evaluation
    const P13_ALPHA_COUNTER      = 4   # D-16 alpha-graded series (simulated substrate draws)
    const P13_CONTINUITY_COUNTER = 5   # D-12 binary-NRE continuity check (reported, not gated)
    const P13_FIXTURE_COUNTER    = 99  # FIXTURES ONLY -- never a reported counter

    """
        P13_GLOBAL_RNG_DISCIPLINE :: String

    The two-layer RNG rule this phase must follow, recorded as a constant so it is part of
    the pre-registration rather than folklore.
    """
    const P13_GLOBAL_RNG_DISCIPLINE = """
    TWO LAYERS, ALWAYS BOTH.
    (1) Counter-based per-index draws come from p13_rng(counter) / p13_fix_rng(counter).
    (2) Random.seed! derived from P13_DEV_SEED MUST be called immediately before EVERY
        `train` call and before EVERY posterior-draw loop, because Flux's
        DataLoader(shuffle = true) draws from the GLOBAL RNG and would otherwise make the
        run irreproducible no matter how carefully layer (1) is threaded
        (13-RESEARCH Pitfall 6). Never Random.seed!(i) inside a loop; never a bare
        global-RNG draw in generation code.
    """

    # =====================================================================================
    # B. The F5 image-size mixture -- READ from the frozen gate file, NEVER retyped
    # =====================================================================================
    # BINDING INVARIANT (F5): rank/coverage/calibration claims hold only under the joint the
    # estimator was TRAINED on, so train-joint MUST equal eval-joint. Both are therefore READ
    # through the isolated `_GC` module above and are NEVER retyped here -- a retyped mixture
    # is a mixture that can silently drift out of agreement with the net it scores.
    const P13_IMSIZE_SET     = _GC.SBC_IMSIZE_SET
    const P13_IMSIZE_WEIGHTS = _GC.SBC_IMSIZE_WEIGHTS

    # =====================================================================================
    # C. The D-05 hypothesis cut
    # =====================================================================================
    # VARIANT (a), "tau-contrast", in full:
    #     COLOC      iff  rho_s >  tau  AND  rho_s - rho_c >  tau
    #     EXCLUSION  iff  rho_s < -tau  AND  rho_s - rho_c < -tau
    #     RANDOM     otherwise
    # One tau governs BOTH factors, so both are stated at the summary's MEASURED resolution --
    # which is D-06's entire point. The rejected variant (b) tests the contrast with a strict
    # inequality and no dead zone, making a 0.001 contrast count as "stands apart" and adding
    # label noise exactly where the heads must be calibrated.
    #
    # THE COUNTER-EXAMPLE THAT MAKES THIS DECISION LOAD-BEARING. A sample at rho = 0.8 under a
    # control at rho = 0.9 has delta_rho < 0. A naive three-way extension of delta_rho
    # (src/amortized/infer.jl:117-121) would publish that as "mutually exclusive" while the
    # sample is strongly colocalized. "Less colocalized than the control" is NOT segregation.
    # The two-factor cut on the rho_sample LEVEL is what prevents this.
    const P13_CUT_VARIANT = :tau_contrast

    """
        p13_tau() -> Float64

    The Tier-2 tripwire. `P13_TAU` is DELIBERATELY ABSENT from the Tier-1 block: D-06 requires
    it to be MEASURED, not chosen. Calling this before the probe has run fails loudly.
    """
    p13_tau() = isdefined(@__MODULE__, :P13_TAU) ? P13_TAU : error("""
        P13_TAU is Tier-2 and does not exist yet.

        D-06 requires tau to be MEASURED from the fixed patch-correlation summary's own
        resolution, not chosen. Plan 13-10 must run spike/p13/run_tau_probe.jl against the
        POST-Phase-11 simulator (P13_TAU_REQUIRE_POST_P11_SIMULATOR) and APPEND the measured
        value in a SECOND guard block below the Tier-1 block in spike/p13/consts.jl.

        Generating labelled three-way data before that point violates D-06: the class cut
        would be drawn at a tau that had not yet been measured.
        """)

    # =====================================================================================
    # D. The D-06 tau probe specification (locked BEFORE the probe consumes a stream)
    # =====================================================================================
    # The probe measures A(delta) = AUC( mbar | rho = 0  vs  mbar | rho = -delta ) and sets
    # tau to the smallest delta on this grid clearing P13_TAU_AUC. Report the WHOLE curve.
    const P13_TAU_DELTA_GRID = (0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20)
    const P13_TAU_AUC        = 0.90   # the "reliably distinguish" bar, pre-registered as a
                                      # NUMBER so "reliably" cannot be argued after the fact
    const P13_TAU_R          = 400    # draws per arm; AUC standard error ~0.03 at A ~ 0.9, so
                                      # the 0.90 bar is resolved to +/- 0.03
    # The statistic is FIXED HERE so it cannot be selected to make a threshold convenient:
    # mbar = the mean of the PRESENT continuous summary rows 1:G^2, weighted by the mask rows
    # G^2+1:2G^2 (so absent patches carry no weight). It is FIT-FREE -- hence honestly
    # pre-registerable -- and it is the same quantity the frozen `ghat` calibration was built
    # on (spike/simulator/ghat.jl:34, E[mu | rho_true]). A LEARNED readout (an LDA on the 64
    # rows, or the NPE's own rho-hat) would be more powerful and would report a SMALLER tau,
    # but it would make tau a property of a trained MODEL rather than of the SUMMARY, which
    # contradicts D-06's wording. Report a learned number as a secondary column if wanted.
    const P13_TAU_STATISTIC  = :mbar_mask_weighted
    # THE SINGLE MOST IMPORTANT DIFFERENCE FROM THE PHASE-11 D-06 PROBE. Phase 11's probe is
    # deliberately PAIRED (same latent field, only geometry changes) because it measures
    # DISPLACEMENT. Tau measures DISCRIMINABILITY FOR A SINGLE UNSEEN IMAGE, where the latent
    # field and the nuisances are unknown. A paired design would report the resolution of a
    # counterfactual the deployed tool never has, and would UNDERSTATE tau by a large factor.
    const P13_TAU_DESIGN     = :unpaired
    # Measure A(0 vs -delta) AND A(0 vs +delta) and take the MAX. The prior atoms are
    # asymmetric by ~2.5x (0.048527 at rho = -0.99 vs 0.019108 at +0.99), so the resolution
    # may be too; a single symmetric tau is what D-05 assumes, and the max is the conservative
    # reading.
    const P13_TAU_BOTH_DIRECTIONS = true
    # Report mean AND median AUC over bootstrap resamples: one degenerate all-missing grid can
    # move mbar by O(1) through the encode_d01 missing-to-zero path.
    const P13_TAU_BOOTSTRAP_B = 1000
    # A RULE, NOT A NUMBER. Phase 11 has not landed, so what is locked is the RULE plus an
    # EXPECTED value; the realized number is asserted against Phase-11's LAMBDA_MAX at probe
    # time (plan 13-10). Resolution degrades with registration uncertainty, so a single tau
    # needs a stated lambda, and the widest rung is the most conservative tau.
    const P13_TAU_REFERENCE_LAMBDA_RULE     = :widest_rung
    const P13_TAU_REFERENCE_LAMBDA_EXPECTED = 3.0
    # THE ABORT CRITERION MUST BE IMMUNE TO ITS OWN RESULT. If A(delta) < P13_TAU_AUC for
    # EVERY delta in the grid then tau is not measurable at this summary's resolution: STOP
    # AND REPORT. Extending the grid afterwards would make "below resolution" unfalsifiable.
    const P13_TAU_ABORT_EXTEND_GRID = false
    # Freezing a tau measured against the PRE-Phase-11 simulator would be a pre-registration
    # error: the deployment distribution is wider (widened SHIFT_PRIOR, chromatic eps), so such
    # a tau is OPTIMISTIC. The probe blocks on Phase 11's SIMULATOR merge, not on its training
    # run (13-RESEARCH E4).
    const P13_TAU_REQUIRE_POST_P11_SIMULATOR = true

    # =====================================================================================
    # E. The D-07 correction design and its verification bars
    # =====================================================================================
    # PRIMARY (D-07-i): stratify CLASS FREQUENCIES ONLY, by whole-class accept/reject, so
    # pi(theta | class) is preserved exactly. THE THEOREM THAT FORCES THIS CHOICE: the scalar
    # per-head correction log(q_pos/q_neg) is exact IFF stratification changed only class
    # FREQUENCIES and not the WITHIN-class theta shape. D-07's literal within-class
    # rho-stratification is therefore PROVABLY NOT repairable by a scalar -- it needs the
    # importance-weighted fallback below. Section F5's negative control demonstrates this
    # empirically rather than asserting it.
    const P13_STRATIFICATION = :class_frequency
    # Target frequencies. The correction is MEASURED from the realized labels regardless and is
    # NEVER assumed to be 0: at tau = 0.10 the pi-level per-head log-odds are log(q_C/q_R) =
    # -0.5333 and log(q_E/q_R) = -0.1654, against the shipped binary net's -0.0102, i.e.
    # 15-50x larger. It cannot be waved away.
    const P13_TARGET_CLASS_FREQ = (exclusion = 1//3, random = 1//3, coloc = 1//3)
    # PRE-DECLARED FALLBACK (D-07-ii): within-range rho-stratification with importance-weighted
    # BCE. Needs the closed-form pi-density of rho plus explicit ATOM handling at rho = -0.99,
    # where a density ratio is undefined and the atom must be its own stratum.
    const P13_STRATIFICATION_FALLBACK = :importance_weighted
    # ONE documented iteration is authorised, and it is authorised HERE, before any result
    # exists. A SECOND iteration is NOT authorised by this file and cannot be authorised by
    # amending this file. Phase 7 was amended TWICE after seeing results and the credibility
    # cost of that is the reason this constant exists.
    const P13_ITERATION_ALLOWANCE = 1
    const P13_ITERATION_TRIGGER = """
    THE ONE PRE-DECLARED CONDITION for spending P13_ITERATION_ALLOWANCE: the D-12 evaluation
    shows the exclusion-class per-class AUC below P13_AUC_FLOOR_EXCLUSION SPECIFICALLY AT HIGH
    |rho| (i.e. the deep-exclusion tail near the -0.99 atom), while the mid-range is resolved.
    The allowance is then spent on switching to P13_STRATIFICATION_FALLBACK (D-07-ii) and
    re-running ONCE. It is NEVER spent on relaxing a threshold after training. A second
    iteration is not authorised by this file.
    """
    # --- F5 verification pass bars (locked, else the verification is unfalsifiable) --------
    # The verification trains the two-head trunk on a CLOSED-FORM conjugate-Gaussian toy whose
    # three-way log Bayes factors are analytic, then compares the CORRECTED head logits against
    # them. No quadrature, no new dependency.
    const P13_F5_CORR_MIN     = 0.99   # cor(corrected logit, analytic log BF) floor
    const P13_F5_MAXABS_TOL   = 0.25   # max |Delta log BF| over the central band
    const P13_F5_CENTRAL_FRAC = 0.90   # the central fraction of the Z grid the bars apply to
    # DELIBERATELY SKEWED toy class frequencies, so a missing correction is impossible to miss.
    const P13_F5_SKEW_FREQ    = (0.60, 0.25, 0.15)
    # Without the negative control the test cannot distinguish "the correction works" from "the
    # correction is unnecessary here": the UNCORRECTED logit must FAIL the same bars by roughly
    # |log(q_A/q_B)|.
    const P13_F5_NEGCTRL_REQUIRED = true

    # =====================================================================================
    # F. The D-09/D-10 architecture and the COPIED training recipe
    # =====================================================================================
    # COPIED VERBATIM from src/amortized/train_ratio.jl:53-54 and :155-163. Copying rather than
    # choosing is what D-10's attribution argument rests on: any difference in the 2-way overlap
    # is then attributable to the HEAD alone, not to a re-tuned recipe.
    const P13_SUMMARY_WIDTH   = 256      # RATIO_SUMMARY_WIDTH  (train_ratio.jl:54)
    const P13_NUM_SUMMARIES   = 64       # RATIO_NUM_SUMMARIES  (train_ratio.jl:53)
    const P13_TRAIN_N         = 48_000   # n                    (train_ratio.jl:156)
    const P13_EPOCHS          = 300      # epochs               (train_ratio.jl:157)
    const P13_BATCHSIZE       = 128      # batchsize            (train_ratio.jl:157)
    const P13_LR              = 2.5e-4   # learning_rate        (train_ratio.jl:158)
    const P13_WEIGHT_DECAY    = 1e-4     # weight_decay         (train_ratio.jl:158)
    const P13_VAL_FRAC        = 0.15     # val_frac             (train_ratio.jl:159)
    const P13_STOPPING_EPOCHS = 40       # stopping_epochs      (train_ratio.jl:159)
    const P13_USE_GPU         = false    # CPU-only is the reproducible baseline (CLAUDE.md)
    # The A5 route-R1 smoke gate that must pass before ANY real training run: 2 epochs on a
    # 200-sample toy, proving NeuralEstimators' `train` accepts the custom estimator subtype
    # AND the custom loss.
    const P13_SMOKE_N      = 200
    const P13_SMOKE_EPOCHS = 2
    # LAMBDA PLACEMENT. The correct form is
    #     Z_pair = vcat(pair_encode(Zs_128, Zc_128), lambda)
    # i.e. lambda appended AFTER the pair encoding. The two wrong forms are
    #   (i)  pair_encode(Zs_129, Zc_129) -- throws on pair_encode's `iseven` guard, because 129
    #        is odd. That is the GOOD failure and it is why the guard is worth relying on.
    #   (ii) vcat(pair_encode(Zs_128, Zc_128), lambda, lambda) -- lambda counted twice, silent.
    const P13_LAMBDA_PLACEMENT = :append_after_pair_encode
    # NEVER WRITE 321 AS A LITERAL. Phase 11 may deliver a conditioning vector of length != 1,
    # so the input width is DERIVED as ratio_input_dim(G) + n_cond, where
    # n_cond = length(encode_lambda(lambda)). A hard-coded 321 would silently be wrong.
    const P13_INPUT_WIDTH_RULE = :ratio_input_dim_plus_ncond

    # =====================================================================================
    # G. The D-12 amended gate (SC2 is amended; declare the amendment BEFORE results)
    # =====================================================================================
    # Labelled evaluation pairs. At the measured class masses under variant (a) this puts
    # >= 2500 samples in each head's RESTRICTED evaluation set, clearing P13_MIN_EVAL_PER_HEAD
    # with margin (each head is scored only on its own class pair, per D-11).
    const P13_GATE_M = 4000
    # Threshold-free one-vs-random discrimination is THE gate. The shipped binary analogue
    # reached AUC 0.994 in spike 014, so 0.90 is a FLOOR and not a target.
    const P13_AUC_FLOOR_COLOC     = 0.90
    const P13_AUC_FLOOR_EXCLUSION = 0.90
    # The 3x3 confusion matrix at argmax(logBF_C, 0, logBF_E) is REPORTED as DESCRIPTIVE and is
    # explicitly NOT a decision rule. Decisions and abstention are Phase 14's scope; a
    # descriptive argmax quoted as a classifier would pre-empt that phase's design.
    const P13_CONFUSION_RULE = :argmax_descriptive
    # The shipped binary-NRE continuity check at a fixed lambda is REPORTED, NEVER GATED,
    # because D-02 gives the two nets different input surfaces -- a disagreement is then a
    # statement about the basis change, not a failure of either net.
    const P13_CONTINUITY_GATED = false

    # =====================================================================================
    # H. The D-13 calibration gate
    # =====================================================================================
    # RE-DECLARED LOCALLY, never imported from spike/validation/consts.jl. Every
    # test/gate/gate_consts_*.jl re-declares these; each pre-registration is self-contained, so
    # it can be read and audited without chasing another file's current value.
    const P13_ECE_GREEN  = 0.05   # ECE <= GREEN  -> :green
    const P13_ECE_YELLOW = 0.10   # ECE <= YELLOW -> :yellow, else :red
    # Matches _bin_calibration's default so the number stays comparable with spike 014's 0.0188.
    const P13_ECE_NBINS  = 10
    # ECE IS THE GATE STATISTIC. MCE IS EXPLICITLY NOT.
    # THE EMPTY-BIN TRAP, IN FULL: _bin_calibration (spike/validation/sbc.jl:82-130) scores an
    # EMPTY bin as predicted_rate = midpoint, observed_rate = 0.0. That contributes ZERO ECE
    # weight (the ECE term is weighted by count_i/total) but FULL MCE weight (MCE is an
    # unweighted max over bins). On a well-separated three-way problem many of the 10 bins will
    # be empty, so MCE would be DOMINATED BY EMPTY BINS and is not a usable gate statistic
    # here. MCE is REPORTED beside its empty-bin count and is NEVER gated. _bin_calibration is
    # shared, gate-lineage code and must NOT be hand-modified to work around this.
    const P13_GATE_STATISTIC = :ece
    # ECE is biased downward at small n and when bins are sparse, so a green ECE on a small
    # evaluation set is not evidence. Pre-register the minimum per-head evaluation count.
    const P13_MIN_EVAL_PER_HEAD = 1000
    # A green ECE on a head whose AUC is below this is a VACUOUS PASS -- a head that says
    # "0.5, always" is perfectly calibrated and perfectly useless. Such a result must be
    # LABELLED as a vacuous pass in the report, not quoted as calibration evidence.
    const P13_VACUOUS_AUC_FLOOR = 0.60
    # D-14: the traffic light is a BAND, not a p-value, so there is no over-power pathology of
    # the kind the M = 2000 SBC point-null had, and no equivalence machinery is needed.
    const P13_TOST_REQUIRED = false

    # =====================================================================================
    # I. The D-16 semi-synthetic alpha ladder (mask-based disjoint reassignment)
    # =====================================================================================
    const P13_ALPHA_LADDER = (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0)
    # THE STRICTLY-POSITIVE BACKGROUND FLOOR. b = quantile(vec(y), P13_ALPHA_BG_QUANTILE),
    # ASSERTED > 0 at runtime. Measured 0.00624 on the positive fixture and 0.00591 on the
    # negative one, both far below the Otsu threshold. An image with more than 5% exact-zero
    # pixels would silently give b = 0 and re-open the _exclude_zero trap, which is exactly why
    # the assertion exists rather than a comment.
    const P13_ALPHA_BG_QUANTILE = 0.05
    # The FROZEN, ALREADY-SHIPPED Otsu rule: M = _calculate_mask(mci)[1]
    # (src/LoadImages.jl:235-241), computed ONCE from the UNMODIFIED ch1 and held FIXED across
    # the whole ladder. Computing it per alpha would make the rungs incomparable, the same
    # reason local_map.jl:208 rejects tile-local Otsu. Using the shipped helper is what makes
    # D-16's segmentation rule a pre-registered choice at zero cost.
    const P13_ALPHA_MASK_RULE   = :calculate_mask_ch1_once
    const P13_ALPHA_N_IMAGES    = 64      # the SIMULATED arm's image count only
    const P13_ALPHA_INVARIANT_RTOL = 1e-10  # sum(y_out) ~ sum(y): total intensity is preserved
    const P13_ALPHA_GATED       = false   # D-15/D-16 are supporting evidence, never a gate

    # TWO-VALUED AND PRE-REGISTERED. The amended D-15 makes the six committed microscopy TIFFs
    # under test/test_images/ the source images D-16's mask-based series is ALSO constructed
    # from, so the ladder has two substrates and the transform is ONE code path typed on
    # AbstractMatrix{Float64} that knows nothing about where its input came from.
    #
    # THE TWO ARMS ARE REPORTED SEPARATELY AND ARE NEVER AVERAGED, because their alpha = 0
    # semantics differ: on simulated substrate the ladder runs random -> exclusion BY
    # CONSTRUCTION, while on real substrate it runs MODERATELY COLOCALIZED -> near-exclusion
    # (measured mbar = +0.3292 positive / +0.2481 negative at alpha = 0). Both are informative;
    # they are not the same experiment.
    #
    # EXECUTING RUNNER PER ARM: :simulated is plan 13-13's run_alpha_series.jl; :real is plan
    # 13-16's run_p13_realimage.jl.
    const P13_ALPHA_SUBSTRATE = (:simulated, :real)

    # THE CORRECTED ZERO INVARIANT. The naive all-pixels-strictly-positive form is FALSE on real
    # data and MUST NEVER BE WRITTEN: the source files already contain exact-zero pixels before
    # any alpha (2 px in positive_c2, 3 px in negative_c2, of 1 414 528 -- 13-RESEARCH J5.3
    # finding 4, Pitfall 2b). Writing it would fail Wave 0 on a CORRECT algorithm at alpha = 0,
    # where the output IS the input. The assertion is
    #     count(iszero, y_out) == count(iszero, y)
    # -- "alpha introduces no NEW zero" -- and the SOURCE zero count is recorded in the
    # artifact. This matters because the zero SET, not the zero count, drives _exclude_zero
    # (src/colocalization.jl:154-169), which DROPS any pixel where EITHER channel is zero.
    const P13_ALPHA_ZERO_INVARIANT = :count_iszero_preserved

    # THE PRE-REGISTERED MASK-FRACTION GUARD. Outside these bounds the runner fires a NAMED
    # error rather than producing an uninterpretable ladder. MEASURED BASIS: on a
    # permissive-threshold sweep a 50% mask gives boost 2.678, while the real Otsu masks at
    # 13.49% (positive) and 22.92% (negative) give 1.509 and 1.598 -- comfortable, with a
    # complement large enough to absorb the redistributed mass.
    const P13_ALPHA_MASK_FRACTION_BOUNDS = (0.01, 0.40)

    # maximum(y_out) <= this bound is ASSERTED AND RECORDED, NEVER SILENTLY CLAMPED. The
    # summary is a per-patch CORRELATION and is therefore invariant to a global affine rescale;
    # the boost bites only through the mask/complement contrast WITHIN a patch, which is the
    # intended mechanism. The measured post-boost maximum stayed at or below 0.75 in every
    # swept arm.
    const P13_ALPHA_MAX_VALUE_BOUND = 1.0

    # THE TRAP THAT DOMINATES THIS DESIGN. _exclude_zero (src/colocalization.jl:154-169) drops
    # any pixel where EITHER channel is zero, so ZEROING ch2 inside the mask DELETES those
    # pixels instead of creating anti-correlation: the measured correlation would then be
    # computed over the complement only (~unchanged, i.e. random), object-dense patches could
    # fall below the 15-survivor floor and go `missing`, and the ladder would read FLAT for a
    # reason that has nothing to do with segregation. Hence the strictly-positive floor b.
    # HONEST CAVEAT ON THE ENDPOINT: at alpha = 1 the mask region is not EMPTY, it is at b. The
    # endpoint is "ch2 reduced to background wherever ch1 has objects", not "ch2 absent" -- do
    # not write "fully disjoint" without that qualifier.

    # =====================================================================================
    # I2. The D-15 real-image QUALITATIVE arm (D-15 AMENDED 2026-07-25)
    # =====================================================================================
    # EVERY KNOB THE REAL-IMAGE ARM READS, SELECTS OR HEADLINES ON IS LOCKED HERE, BEFORE THE
    # SIX COMMITTED TIFFS ARE EVER OPENED -- so no lambda and no ladder can be chosen after
    # seeing a real-image verdict.
    #
    # THE PHYSICAL ANCHORS OF THE PROVENANCE MANIFEST ARE NOT USED. Both physical-primary rows
    # are split = sealed_holdout, sha256 = "PENDING-FETCH", bytes = 0, reserved for the Phase-16
    # blind evaluation. The SANCTIONED SUBSTITUTE is the six committed microscopy TIFFs under
    # test/test_images/, which are READ-ONLY INPUT to this phase.

    # The three sentinels that make the arm's posture MACHINE-CHECKED rather than merely
    # written. The six TIFFs carry NO colocalization ground-truth label, so the arm can show the
    # three-way verdicts behave SENSIBLY on real microscopy but CANNOT show they are CORRECT.
    # No real-image quantity has a pass/fail threshold anywhere in this file -- and if no
    # threshold exists, the arm cannot accidentally become a gate. Labelled real segregation
    # validation is deferred to Phase 16's blind evaluation.
    const P13_REAL_IS_GATED          = false
    const P13_REAL_QUALITATIVE_ONLY  = true
    const P13_REAL_READ_ONLY         = true

    # The pre-registered arms, mirroring 11-10-PLAN.md's shape. The primary channel pair is the
    # EXACT pair the frozen `ghat` calibration anchor was measured on; the second arm is
    # labelled REDUNDANCY and is never a separate claim.
    const P13_REAL_CONDITIONS      = ("positive", "negative")
    const P13_REAL_SAMPLE          = "positive"
    const P13_REAL_CONTROL         = "negative"
    const P13_REAL_CHANNEL_PAIR    = (1, 2)
    const P13_REAL_REDUNDANCY_PAIR = (1, 3)

    # GRID TRUNCATION, STATED EXPLICITLY BECAUSE A READER WHO COMPUTES 1028/8 WILL EXPECT
    # 128.5-px patches. patch() TRUNCATES to the largest exact multiple
    # (src/colocalization.jl:37-43). At G = 8 the column axis divides exactly, the bottom 4 rows
    # of the row axis (0.39%) are not used, every patch holds 22 016 px -- three orders of
    # magnitude above the 15-survivor floor -- and ZERO patches went `missing` on either
    # fixture. THEREFORE NO CROP, RESIZE, PAD OR num_patches CHANGE EXISTS OR MAY BE ADDED:
    # adding one would be a new unregistered preprocessing choice and would diverge from the
    # frozen _real_anchor() path (Pitfall 2c).
    #
    # TRANSPOSITION FOOTNOTE. IMSIZE_SET (spike/data/seeding.jl:52-58) carries this frame size
    # WITH THE AXES SWAPPED, at weight 0.03 -- so the training distribution deliberately
    # includes it. Per-patch pixel count and truncation are identical; only the column-major
    # ordering of the 8x8 grid corresponds to a transposed spatial layout. The A14 falsifier
    # checks that in one line.
    const P13_REAL_IMSIZE               = (1028, 1376)
    const P13_REAL_GRID_TRUNCATION_ROWS = 4

    # THE FROZEN LINEAGE ANCHORS, recorded at spike/simulator/ghat.jl:39 and reproduced exactly
    # during research. The ingestion REGRESSES against these, which is what makes the real
    # substrate a lineage-carrying reference point in this project's own simulator calibration
    # rather than an arbitrary set of TIFFs.
    const P13_REAL_ANCHOR_MBAR = (positive = 0.3292, negative = 0.2481)
    # An ingestion REGRESSION tolerance -- explicitly NOT a Phase-13 pass/fail bar. Documented
    # exemption 1 of 2 in the not-a-gate source assertion.
    const P13_REAL_ANCHOR_TOL  = 1e-3

    # The real arm deliberately reuses the simulated ladder so the two curves are comparable
    # RUNG FOR RUNG. The measured summary-level sign crossing sits at alpha ~ 0.49 (positive)
    # and ~ 0.46 (negative) -- strictly INTERIOR to this grid rather than pinned at an endpoint,
    # which is what makes the crossing point informative.
    const P13_REAL_ALPHA_GRID = P13_ALPHA_LADDER

    # A RULE PLUS AN EXPECTED VALUE, exactly as P13_TAU_REFERENCE_LAMBDA_RULE is handled,
    # because Phase 11 has not landed; the realized LAMBDA_MAX is asserted against Phase-11's
    # constant at read time (plan 13-16). This closes 13-RESEARCH Open Question 6.
    #
    # THE REASONING. A real image has NO KNOWN lambda: its registration error is unmeasured and
    # unmeasurable from the file, so there is no correct value to plug in. The arm therefore
    # reports the WHOLE sweep and headlines the widest, most conservative rung. A conservative
    # lambda WEAKENS the evidence, so a verdict that survives it is the credible one.
    #
    # [ASSUMED] (A11): widefield chromatic aberration plus filter-cube/stage repeatability on
    # these particular microscopes is at least 1 px, which is why LAMBDA_MIN = 0.25 is not a
    # credible operating point for them. This is a DOMAIN JUDGEMENT about typical widefield
    # systems, NOT a measurement of this instrument, and it GATES NOTHING.
    #
    # EXPLICITLY FORBIDDEN: estimating lambda from the images (e.g. a cross-correlation peak
    # offset). That is a new estimator with its own validation burden, inside a phase that is
    # not about registration.
    const P13_REAL_LAMBDA_READS_RULE       = :full_phase11_ladder
    const P13_REAL_LAMBDA_HEADLINE_RULE    = :widest_rung
    const P13_REAL_LAMBDA_HEADLINE_EXPECTED = 3.0

    # 13-RESEARCH Open Question 7 answered YES. The Phase-13 net's OOD verdict on these
    # fixtures is measured and reported BESIDE the shipped net's: agreement corroborates the
    # named limit, divergence is a finding about what the registration-aware basis buys on real
    # data. The shipped artifacts/grid_8/ood_nulls_8.jld2 is read as a frozen COMPARISON
    # REFERENCE ONLY and is never an evidence basis -- D-02's no-grid-8-fallback rule governs
    # the EVIDENCE read, not the quoting of a comparison.
    const P13_REAL_OOD_COMPARISON = true
    # The measured shipped-net numbers, frozen so the comparison has a stable reference. BOTH
    # fixtures are flagged out-of-distribution at 2.4x over threshold in BOTH read directions.
    # THIS IS NOT A BUG AND NOT A REASON TO SKIP THE CHECK: it is the OOD channel doing its job,
    # it is arguably the most honest single number the real-data arm produces, it must NEVER be
    # suppressed or "fixed" (Pitfall 13), and every real-image result is printed next to its OOD
    # verdict and its lambda.
    const P13_REAL_OOD_SHIPPED_DENSITY = 433.69
    # A RECORDED reference value belonging to the SHIPPED net, not a Phase-13 bar. Documented
    # exemption 2 of 2 in the not-a-gate source assertion.
    const P13_REAL_OOD_SHIPPED_THRESHOLD = 179.14

    # 13-RESEARCH Open Question 8 answered YES: carry the honesty note verbatim.
    const P13_REAL_NAMING_CORRECTION = """
    NAMING CORRECTION (must appear in the report). The folder names positive/ and negative/ are
    the original package's biological test conditions, NOT colocalization labels. The
    "negative" pair measures mean patch correlation +0.2481 (rho_true ~ +0.215) and is
    therefore a POSITIVELY CORRELATED pair, not an anti-correlated one. Nothing in this phase
    may treat negative/ as an exclusion example, and a reader who assumes otherwise will
    misread every real-image figure in the phase.
    """

    const P13_REAL_SUBSTITUTION_RECORD = """
    TARGET-SUBSTITUTION RECORD (11-10-PLAN.md house shape). D-15 originally named the physical
    mitochondria anchor of the frozen provenance manifest as the real-data check. Verification
    showed both physical-primary rows are sha256 = "PENDING-FETCH", bytes = 0 and
    split = sealed_holdout -- unfetched, and reserved behind the anti-snooping control for the
    Phase 16 blind evaluation. Consuming one here would irreversibly burn Phase 16 on the very
    hypothesis Phase 16 exists to evaluate. The sanctioned substitute is the six committed
    microscopy TIFFs under test/test_images/. The sealed-holdout accessor is NOT called -- not
    in a script, not in a test, not behind a flag -- and no seal-break escape hatch exists.
    """

    # =====================================================================================
    # J. Decoupling, deferrals and DECLARED deviations
    # =====================================================================================
    const P13_SRC_UNTOUCHED = true
    # The one real segregated anchor is TRIPLY unavailable: sha256 = "PENDING-FETCH", the data
    # directory is empty, and split = sealed_holdout behind the anti-snooping seal. Phase 13
    # does NOT open the seal.
    const P13_PHYSICAL_ANCHOR_DEFERRED_TO = "Phase 16"
    # D-08 is right that `log_bf_simplex` (src/results.jl:172) is a misnomer -- the three values
    # are (log BF(C:R), 0, log BF(E:R)), two free numbers plus a structural zero, not a point on
    # a 2-simplex. The rename is RECORDED here and DEFERRED to a future productionization phase
    # (13-RESEARCH G3 Open Question 1). src/ stays byte-unchanged in this phase.
    const P13_RESULTS_RENAME_DEFERRED = true
    # Naming the model changes IN ADVANCE is the D-04 discipline. The report must enumerate
    # exactly these four and no others; anything else that changes is an UNDECLARED deviation.
    const P13_DECLARED_DEVIATIONS = (
        "a custom ThreeWayEvidenceNet <: NeuralEstimator subtype instead of RatioEstimator, because RatioEstimator's loss is hard-coded and a passed loss is silently ignored in v0.2.1",
        "a single Dense(num_summaries, 2) two-logit head instead of a one-logit model-index head",
        "a masked two-term BCE loss with per-head participation weights, so each head trains only on its own class pair (D-11)",
        "a per-head evidence-scale correction measured with measure_head_log_odds, which is NOT the binary `- log_prior_odds` term and must not be copied from it",
    )

    # =====================================================================================
    # EXECUTABLE self-checks (cheap; no inference, no simulation, no file writes)
    # =====================================================================================
    # D-01: both fresh streams are disjoint from every reserved, burned and RECOMPUTED seed.
    @assert !(UInt64(P13_DEV_SEED) in _p13_forbidden()) "P13_DEV_SEED collides with a forbidden seed"
    @assert !(UInt64(P13_FIX_SEED) in _p13_forbidden()) "P13_FIX_SEED collides with a forbidden seed"
    @assert UInt64(P13_DEV_SEED) != UInt64(P13_FIX_SEED) "the fixture stream must differ from the reported stream"
    @assert !(P13_SALT in P13_REPO_SALTS) "P13_SALT reuses an existing repository salt"
    # A duplicate inside the forbidden list would silently shrink the set; catch that.
    @assert length(_p13_forbidden()) == length(_p13_forbidden_list()) "duplicate entry in the forbidden-seed list"
    # The recompute path actually ran (a comment naming a derived seed is not evidence).
    @assert any(v -> UInt64(v) in _p13_forbidden(), values(_GC.PROD_SEED))
    @assert any(v -> UInt64(v) in _p13_forbidden(), values(_GC.PROD_SEED_V2))
    # Fixtures must never ride a counter a reported number rides on.
    @assert P13_FIXTURE_COUNTER ∉ (P13_TAU_COUNTER, P13_DATAGEN_COUNTER, P13_GATE_COUNTER,
                                   P13_ALPHA_COUNTER, P13_CONTINUITY_COUNTER)
    # Ladder and grid integrity.
    @assert issorted(P13_TAU_DELTA_GRID) && first(P13_TAU_DELTA_GRID) > 0.0
    @assert issorted(P13_ALPHA_LADDER) && first(P13_ALPHA_LADDER) == 0.0 && last(P13_ALPHA_LADDER) == 1.0
    @assert sum(values(P13_TARGET_CLASS_FREQ)) == 1
    @assert sum(P13_IMSIZE_WEIGHTS) ≈ 1.0
    @assert length(P13_IMSIZE_SET) == length(P13_IMSIZE_WEIGHTS)
    # D-13: the band is ordered and ECE -- not MCE -- is the gate statistic.
    @assert P13_ECE_GREEN < P13_ECE_YELLOW
    @assert P13_GATE_STATISTIC === :ece
    # D-16: two substrates, the CORRECTED zero invariant, and a non-degenerate mask guard.
    @assert P13_ALPHA_SUBSTRATE == (:simulated, :real)
    @assert P13_ALPHA_ZERO_INVARIANT === :count_iszero_preserved
    @assert 0.0 < first(P13_ALPHA_MASK_FRACTION_BOUNDS) < last(P13_ALPHA_MASK_FRACTION_BOUNDS) < 1.0
    @assert 0.0 < P13_ALPHA_BG_QUANTILE < 1.0
    # D-15: the not-a-gate posture, the shared ladder, and the arm assignment.
    @assert P13_REAL_IS_GATED == false && P13_REAL_QUALITATIVE_ONLY == true && P13_REAL_READ_ONLY == true
    @assert P13_REAL_ALPHA_GRID == P13_ALPHA_LADDER
    @assert P13_REAL_SAMPLE in P13_REAL_CONDITIONS && P13_REAL_CONTROL in P13_REAL_CONDITIONS &&
            P13_REAL_SAMPLE != P13_REAL_CONTROL
    @assert P13_REAL_LAMBDA_HEADLINE_RULE === :widest_rung &&
            P13_REAL_LAMBDA_READS_RULE === :full_phase11_ladder
    @assert P13_REAL_OOD_COMPARISON == true &&
            P13_REAL_OOD_SHIPPED_DENSITY > P13_REAL_OOD_SHIPPED_THRESHOLD
    @assert occursin("biological", P13_REAL_NAMING_CORRECTION) &&
            occursin("0.2481", P13_REAL_NAMING_CORRECTION)
    @assert occursin("Phase 16", P13_REAL_SUBSTITUTION_RECORD)
    # Deferrals and the fixed-length declared-deviation enumeration.
    @assert P13_SRC_UNTOUCHED == true && P13_RESULTS_RENAME_DEFERRED == true
    @assert P13_PHYSICAL_ANCHOR_DEFERRED_TO == "Phase 16"
    @assert length(P13_DECLARED_DEVIATIONS) == 4
    @assert P13_ITERATION_ALLOWANCE == 1
    # TIER 2 MUST NOT BE PRESENT IN THE TIER-1 COMMIT. Plan 13-10 does not DELETE this
    # assertion -- it RELOCATES it above the appended Tier-2 block, so the Tier-1 guard still
    # proves that tau did not exist when the bars were set.
    # RELOCATED BY PLAN 13-10, NOT DELETED. At the Tier-1 commit
    # c42cc8e4a3d33b6ad95e3b5d5af58d84066ab4ad ("feat(13-01): lock the Phase-13 Tier-1
    # pre-registration") this block ended with the statement
    #     @assert !isdefined(@__MODULE__, :P13_TAU) "P13_TAU is Tier-2 and must be absent from the Tier-1 commit"
    # and it HELD: `git show c42cc8e:spike/p13/consts.jl` carries that line, and no P13_TAU
    # existed anywhere in the repository at that commit. THAT is the pre-registration
    # guarantee -- the bars were set while tau did not exist -- and it is now a matter of
    # git history rather than a live runtime check, because the Tier-2 block appended below
    # legitimately supplies the value. The EQUIVALENT, Tier-2-AWARE live assertion sits
    # immediately above that block: tau may exist only in the company of its probe
    # provenance, so a tau smuggled in from anywhere other than the appended block still
    # fails loudly.
end

# =========================================================================================
# THE RELOCATED TIER-1 SELF-CHECK (plan 13-10)
# =========================================================================================
# The Tier-1 block asserted that tau did not exist. That statement is preserved verbatim in
# the comment at the foot of that block, together with the commit at which it held. What
# survives as EXECUTABLE code is its purpose, restated so it is still falsifiable now that a
# measured tau legitimately exists: TAU MAY ONLY EVER ARRIVE TOGETHER WITH ITS PROBE
# PROVENANCE. A P13_TAU defined without P13_TAU_PROBE_ARTIFACT is a tau that came from
# somewhere other than the appended Tier-2 block below -- a hand-edit, a REPL, a caller
# setting it -- which is precisely the failure the original assertion existed to catch.
# It sits ABOVE the Tier-2 block (outside the Tier-1 guard, so it also runs on a re-include)
# exactly as the Tier-1 header requires.
@assert (!isdefined(@__MODULE__, :P13_TAU) ||
         isdefined(@__MODULE__, :P13_TAU_PROBE_ARTIFACT)) "P13_TAU exists without its probe provenance: tau was set outside the appended Tier-2 block"

# =========================================================================================
# TIER 2 -- PROBE-DERIVED, APPENDED AFTER THE MEASUREMENT
# =========================================================================================
# THIS IS THE ONE NUMBER IN PHASE 13 THAT WAS MEASURED RATHER THAN CHOSEN, AND THE ORDER OF
# OPERATIONS IS WHAT MAKES IT LEGITIMATE.
#
# The probe ran on 2026-07-27 (UTC), reported by `spike/p13/run_tau_probe.jl` against the
# POST-Phase-11 simulator, and persisted its whole A(delta) curve to
# `spike/p13/tau_probe_report.jld2`. The provenance it recorded:
#
#   repo HEAD at the probe run  cecb69c58770a5b45bec8151c657d9024b999f80
#   simulator commit            78dc37f517ad1b8e1e71afa963691f7ddefe5f63
#   simulator guard             post_p11 = true; theta arity 8 with chromatic_eps present,
#                               shift-prior half-width 3.0 == expected 3.0
#   stream                      P13_DEV_SEED at counter P13_TAU_COUNTER, R = 400 per arm,
#                               unpaired arms on two provably disjoint Philox keys
#
# THE MEASUREMENT, on the frozen grid, as magnitudes max(A, 1-A), both directions:
#
#   delta   A_neg      A_pos      A = max    bootstrap median
#   0.02    0.529025   0.592806   0.592806   0.592806
#   0.03    0.566481   0.625869   0.625869   0.625184
#   0.05    0.639375   0.689925   0.689925   0.688956
#   0.075   0.723788   0.760544   0.760544   0.761053
#   0.10    0.800338   0.819356   0.819356   0.821203
#   0.15    0.916356   0.905794   0.916356   0.917831   <-- FIRST delta to clear the bar
#   0.20    0.978556   0.953612   0.978556   0.978541
#
# No arm produced a degenerate all-absent patch grid, so every AUC rests on the full 400
# draws per arm. The curve is monotone increasing in delta, which is the shape the design
# predicts and the absence of the "statistic is not tracking rho at all" failure mode.
#
# tau = 0.15: the SMALLEST delta on the frozen P13_TAU_DELTA_GRID whose A(delta) reached
# the pre-registered bar P13_TAU_AUC. The grid was not extended, the bar was not moved, R
# was not raised and no seed was rerolled -- and none of that was needed, because the bar
# was cleared on the grid as committed.
#
# THIS BLOCK WAS COMMITTED SEPARATELY, AFTER THE ARTIFACT COMMIT
# 581602a73a3f180f490f1c421c43028aa7cf9d47 that carries the runner and the curve. That
# separation is the auditable evidence that tau was MEASURED BEFORE IT WAS LOCKED: the
# measurement exists in git history at a commit that contains no Tier-2 value, so no reader
# has to take the ordering on trust.
#
# Guarded on :P13_TAU, a DIFFERENT sentinel from the Tier-1 block's :P13_DEV_SEED, so this
# block is independently idempotent under a re-include.
if !isdefined(@__MODULE__, :P13_TAU)
    const P13_TAU = 0.15                    # first grid delta clearing the bar; A(tau) = 0.916356
    const P13_TAU_MEASURED_AUC = 0.91635625 # the realized A(tau) = max(A_neg, A_pos) at tau
    const P13_TAU_PROBE_SHA = "cecb69c58770a5b45bec8151c657d9024b999f80"
    const P13_TAU_SIMULATOR_SHA = "78dc37f517ad1b8e1e71afa963691f7ddefe5f63"
    const P13_TAU_PROBE_ARTIFACT = "spike/p13/tau_probe_report.jld2"
    # TAU IS REALLY TAU-OF-LAMBDA, so the reference registration uncertainty must travel with
    # the number or the number is not interpretable. P13_TAU_REFERENCE_LAMBDA_RULE
    # (= :widest_rung) was RESOLVED at probe time against Phase-11's own frozen ladder: the
    # widest SC2 rung is LAMBDA_MAX = 3.0, which is also the half-width of the live
    # SHIFT_PRIOR the probe marginalized over. Because the probe integrates over the FULL
    # shift prior rather than conditioning on a rung, the measured tau is already the tau at
    # the widest registration uncertainty -- the conservative reading the rule asks for.
    const P13_TAU_REFERENCE_LAMBDA = 3.0

    # Executable, not decorative: tau must be a grid point (never an interpolated or rounded
    # value) and must actually have cleared the frozen bar.
    @assert P13_TAU in P13_TAU_DELTA_GRID "P13_TAU is not a point on the frozen pre-registered grid"
    @assert P13_TAU_MEASURED_AUC >= P13_TAU_AUC "P13_TAU_MEASURED_AUC does not clear the frozen P13_TAU_AUC bar"
end
