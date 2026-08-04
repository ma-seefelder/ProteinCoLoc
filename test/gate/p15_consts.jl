#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/p15_consts.jl --- the PHASE-15 pre-registration (calibration operating envelope).
#
# This file is the executable form of
# `.planning/phases/15-calibration-operating-envelope-and-ci-gate/15-CONTEXT.md`, decisions
# D-01..D-14 together with its appended amendments D-02a, D-04a, D-08a and D-10a. Read that
# document before reading this file.
#
# ============================ WHY THIS IS A FORK, NOT AN EXTENSION ============================
# It MIRRORS `gate_consts_8_v2.jl` and deliberately does NOT extend it. The amended v2
# pre-registration is bound to exactly ONE run on the fresh `PROD_SEED_V2[8]` (`run_gate.jl:47-54`,
# amendment protocol §6.3), and THAT RUN IS SPENT — it produced
# `artifacts/amended_v2/grid_8/gate_report_8.jld2` with `seed = 1135605683775656488`. Extending
# the v2 file would entangle this phase with a closed amendment and would consume a stream that is
# no longer fresh.
#
# The amended RULE SET is therefore reproduced here BY BYTE COPY rather than restated, so that
# "Phase 15 scores its SBC arm under the same amended rules the shipped report was scored under"
# is a per-constant EQUALITY CHECK (`test/gate/test_p15_consts.jl`, testset 4) rather than a claim.
# `GATE_CONSTS_VERSION` and `GATE_AMENDMENT_DOC` are the ONLY shared-name values allowed to differ.
#
# LOADING: this file defines the SAME const names as `gate_consts_8.jl` / `gate_consts_8_v2.jl`
# and is guarded on :SBC_M, so no two of them may EVER be loaded into the same module. It is
# selected only by the explicit `run_gate.jl --consts p15_consts.jl` override, never by the
# default `gate_consts_<G>.jl` path (D-01) — consuming a pre-registration must stay an auditable
# invocation and never a side effect.
#
# =============================== THE LICENSING STANDARD FOR A BAR =============================
# `.planning/STATE.md` records that this milestone has concluded "the bar was wrong, not the
# model" FOUR times, and that the count is now high enough that "we pre-registered it" no longer
# settles an argument on its own. The standard each of those amendments had to meet — shown, FROM
# EVIDENCE INDEPENDENT OF THE MACHINERY UNDER AUDIT, to measure something other than what it named
# — applies to every bar below. So every bar below states its derivation NEXT TO IT. A bar with no
# derivation is precisely the liability the Stage-1-ceiling adjudication already cost this project.
#
# ================================ TWO-TIER STRUCTURE — READ FIRST =============================
#   TIER 1 is the single guard block in THIS file, keyed on :SBC_M, committed before any Phase-15
#     rung exists: the ECE break criterion (D-04a), the OOD fire-rate margin (D-07), the ladder
#     endpoints for all seven axes (D-03), the CI golden tolerance (D-09), the D-13 declarations.
#   TIER 2 IS NOT ONE BLOCK. It is a SEQUENCE of append-only blocks, each keyed on its OWN
#     sentinel const, because a second append that reused an earlier block's sentinel would be
#     silently SKIPPED once that block exists, and an append with no guard at all is a hard `const`
#     redefinition on re-include. THREE sentinels are RESERVED here, with the plan that opens each:
#       :P15_GOLDEN_RANK_DIGEST   -- opened by 15-06 (the committed D-09 fast-gate golden)
#       :P15_ECE_ANCHOR_MEASURED  -- opened by 15-07 (the rung-0 ECE anchor, measured at
#                                    P15_SWEEP_M on the P15 stream; the D-04a break anchor)
#       :P15_ITERATION_SPENT      -- opened ONLY if D-13's single iteration allowance is spent
#   NO CONSTANT IN EITHER TIER IS EVER *EDITED*. Constants are only ever APPENDED, so the git
#   history of this file is itself the pre-registration audit trail. A diff that MODIFIES a line
#   below is, by construction, a pre-registration breach (CONVENTIONS C-03 declaration-line diff:
#   `git diff -U0 -- test/gate/p15_consts.jl | grep -E '^[+-]\s*const '`).
#
# SEED DISCIPLINE (D-11): a FRESH `PROD_SEED_P15`, derived by the SAME salted-Philox redraw rule
# from a DISTINCT (P15_MASTER, P15_SALT) pair, re-drawn until disjoint from every burned stream in
# the repository — including BOTH frozen gate families, which are DERIVED and are therefore
# RECOMPUTED below rather than trusted from a comment. `prod_seed`/`prod_rng` alias onto the P15
# stream so `gate_global_seed(G)` (`harness.jl:82`), the GLOBAL half of reproducibility, rides it
# too. Asserted in-file AND independently in `test/gate/test_p15_consts.jl`, whose negative fixture
# SHOWS the guard firing on a deliberate collision (CONVENTIONS C-04 rule 3).
#
# NOTHING HAS BEEN RUN AGAINST THIS FILE AT THE TIME OF COMMIT. It is a specification, frozen
# before use.
#
# DEPENDENCIES: Random123 only, exactly as every other gate pre-registration. The frozen-file
# SHA-256 pins are LITERALS here (recomputed and compared in the testset) so this file performs no
# file I/O at load time and a pin can never be tautologically self-satisfying.

if !isdefined(@__MODULE__, :SBC_M)
    import Random123: Philox4x     # counter-based RNG for the fresh disjoint per-grid stream

    # The ONLY two shared-name values allowed to differ from `gate_consts_8_v2.jl`. Every other
    # shared constant is asserted per-name equal to the v2 file in `test_p15_consts.jl` testset 4.
    const GATE_CONSTS_VERSION = 15
    const GATE_CONSTS_GRID    = 8
    const GATE_AMENDMENT_DOC  =
        ".planning/phases/15-calibration-operating-envelope-and-ci-gate/15-CONTEXT.md"

    # =====================================================================================
    # SBC — reported amended 8×8 run
    # =====================================================================================
    # UNCHANGED from v1 (deliberately: changing M would change power in a direction the author
    # knew he needed, which is exactly the move a post-hoc amendment must not make).
    const SBC_M          = 2000        # pre-registered SBC draws (θ*~π → simulate → CPU-infer)
    const SBC_L          = 999         # posterior draws per SBC draw; L+1 = 1000
    const SBC_BINS       = 50          # ECE/MCE RELIABILITY bins ONLY (sbc_calibration) — v1 value
    const SBC_ECE_GREEN  = 0.05        # ECE ≤ GREEN  → :green traffic-light verdict  (v1 value)
    const SBC_ECE_YELLOW = 0.10        # ECE ≤ YELLOW → :yellow, else :red            (v1 value)

    # --- A1 : multiplicity control (FWER, not per-test α) --------------------------------
    #
    # DEFECT. `sbc_gate` required all 8 per-parameter KS tests to clear α = 0.05 independently.
    # For a PERFECTLY CALIBRATED net that conjunction fails with probability
    #     FWER = 1 − (1 − 0.05)^8 = 1 − 0.66342 = 0.33658 ≈ 33.7 %
    # — one run in three. This is arithmetic on the number of tests (8, fixed by
    # SBC_PARAM_LABELS) and the frozen per-test α (0.05); it needs no result to derive.
    #
    # RULE. Holm–Bonferroni step-down at FWER = 0.05, applied SEPARATELY to the 8 KS p-values
    # and to the 8 χ² p-values. Sorted p₍₁₎ ≤ … ≤ p₍₈₎; adjusted p̃₍ᵢ₎ = min(1, running-max of
    # (8 − j + 1)·p₍ⱼ₎ for j ≤ i); the arm PASSES iff all(p̃ > 0.05), i.e. Holm rejects nothing.
    # Step-down critical values: 0.00625 0.007143 0.008333 0.01 0.0125 0.016667 0.025 0.05.
    #
    # WHY HOLM. Valid under ARBITRARY dependence — required, because the 8 columns share the same
    # M simulated datasets and the same posterior draws, so they are dependent by construction
    # with an unknown structure. Uniformly at least as powerful as plain Bonferroni (step-down
    # dominance). Šidák is rejected (assumes independence). Benjamini–Hochberg is rejected: the
    # gate's claim is a CONJUNCTION ("all 8 calibrated"), for which family-wise — not
    # false-discovery — control is the right error concept, and BH would be MORE permissive here,
    # the wrong direction for an amendment written after a failure.
    const SBC_MULTIPLICITY = :holm     # Holm–Bonferroni step-down (see holm_adjusted / holm_pass)
    const SBC_N_TESTS      = 8         # 7 θ + the paired-draw Δρ column (D-01)
    const SBC_KS_FWER      = 0.05      # family-wise error rate target, KS family
    const SBC_CHI2_FWER    = 0.05      # family-wise error rate target, χ² family
    # NOTE: v1's SBC_KS_ALPHA / SBC_CHI2_ALPHA are RETIRED as gate inputs — deliberately NOT
    # redefined here, so any code still reading a per-test α against this file errors loudly
    # instead of silently applying an uncorrected conjunction.

    # --- A1b : χ² rank-bin count, decoupled from the ECE binning -------------------------
    #
    # DEFECT. v1 used ONE constant (SBC_BINS = 50) for two unrelated statistics: the χ² rank bins
    # and the ECE reliability bins. Neither could then be set on its own merits.
    #
    # RULE. SBC_CHI2_BINS = 20 for χ² ONLY; SBC_BINS = 50 retained for ECE ONLY.
    #
    # WHY 20. (i) DF ALLOCATION: the alternatives SBC exists to detect are LOW-frequency — a
    # location shift is a monotone ramp, over/under-confidence is a U / inverted-U. A 50-bin χ²
    # spends 49 df, of which only the first few carry calibration meaning; the rest dilute power
    # against those alternatives while adding power against high-frequency structure that has no
    # calibration interpretation. 20 bins (19 df) reallocates the same M toward the alternatives
    # of interest — a power REALLOCATION, not a uniform reduction, and its direction is fixed by
    # which alternatives matter, a question answerable before any data. (ii) ASYMPTOTIC ACCURACY:
    # expected count 2000/20 = 100 per bin vs 40 at 50 bins. (iii) EVENNESS: L+1 = 1000 is
    # divisible by 20 (exactly 50 ranks per bin), preserving v1's exact-evenness requirement.
    # HONEST LIMIT: the decoupling is forced; the value 20 is convention-informed judgement, not
    # uniquely derived (see amendment §5).
    #
    # NOT FIXED HERE (deliberately): the point null of exact uniformity is known false a priori
    # for any finite-capacity flow. A practical-equivalence band would address that, and is
    # exactly the escape hatch a post-hoc amendment must NOT introduce. M stays 2000.
    const SBC_CHI2_BINS = 20

    # --- F5 : the gate's simulate distribution MUST equal the training joint --------------
    #
    # DEFECT. SBC rank uniformity holds only under the joint p(θ,y) = π(θ)p(y|θ) the estimator was
    # TRAINED on. v1's scalar SBC_IMSIZE = (256,256) was inherited from a documented
    # COMPUTE-BUDGET deviation in run 07-05, not from a scientific decision. D-08 anchors the
    # design on the real paper regime (1376×1028); D-09 makes image dimensions a configurable
    # parameter to be validated across a realistic range. Training moves to a realistic mixture;
    # if the gate does not move with it the SBC arm proves calibration on a distribution the net
    # will never see. This is a VALIDITY CONDITION OF THE METHOD, derivable with no data at all.
    #
    # RULE. A pre-registered categorical mixture, drawn PER SBC DRAW from the gate rng via the
    # SHARED ProteinCoLoc.sample_imsize (never a hand-rolled second sampler — the gate and
    # training joints must not be able to drift apart in implementation).
    #
    #   size          w      role                                    measured s/pair (32 thr)
    #   512×512      0.40    cheap lower bracket, bounds mixture cost   0.0274
    #   1024×1024    0.25    intermediate bracket below the anchor      0.0766
    #   1376×1028    0.25    THE D-08 REAL-DATA ANCHOR                  0.0974
    #   2048×2048    0.10    upper bracket (cost ≈ 2.965× the anchor)   ≈0.29 (area-extrapolated)
    #
    # 256² is EXCLUDED: below the realistic acquisition range, and its presence in v1 was the
    # budget artifact being corrected.
    # E[cost] = .40·.0274 + .25·.0766 + .25·.0974 + .10·.29 = 0.08346 s/pair
    #         ⇒ 50 000 pairs ≈ 4173 s ≈ 1.16 h  (CHEAPER than 50k at the anchor alone, ≈1.35 h).
    # HONEST LIMIT: the weights are reasoned (anchor mass, bracketing, cost) but not uniquely
    # determined. What is NOT judgement is that gate mixture == training mixture.
    const SBC_IMSIZE_SET     = ((512, 512), (1024, 1024), (1376, 1028), (2048, 2048))
    const SBC_IMSIZE_WEIGHTS = (0.40, 0.25, 0.25, 0.10)          # Σ = 1.0 (asserted below)

    # SENTINEL, not a size. v1 threaded SBC_IMSIZE as a SCALAR tuple into
    # `sim.simulate_pair(rng, θ; imsize = ...)`. Any code path that still does so must FAIL LOUDLY
    # rather than silently revert to a single image size and quietly invalidate the SBC proof.
    const SBC_IMSIZE = :mixture

    # BINDING INVARIANT (mechanically enforced, not a comment): SBC_IMSIZE_SET /
    # SBC_IMSIZE_WEIGHTS must equal the `imsize_set` / `imsize_weights` the gated net was TRAINED
    # with, read back via `ProteinCoLoc.training_imsize_provenance(...)` (persisted per F6). On
    # any mismatch — or on `recorded = false`, which includes ALL three currently frozen bundles —
    # the gate must refuse to run with status `:provenance_mismatch`. A net whose training
    # image-size distribution is `:unknown` CANNOT be SBC-gated under this amendment.
    const SBC_REQUIRE_IMSIZE_PROVENANCE = true

    # --- SBC : fast fixture (per-task smoke; NEVER the reported numbers) ------------------
    const SBC_FIX_M    = 8             # tiny SBC-draw count for the quick harness smoke
    const SBC_FIX_L    = 15            # tiny posterior-draw count (fixture; L+1 = 16)
    const SBC_FIX_BINS = 4             # tiny rank-bin count (fixture; 16 divisible by 4)

    # =====================================================================================
    # A2 — Amortized Bayes factor (BF gate; non-clamped baseline)
    # =====================================================================================
    #
    # DEFECT 1 (precision). v1 ran n = BF_SWEEP_N = 25 paired draws (recorded n 15–20 after
    # dropping non-finite pairs) and compared a POINT estimate to 0.95. The Fisher-z 95 % CI for a
    # correlation near 0.95 is ±0.07–0.10 wide at those n: the gate cannot discriminate 0.94 from
    # 0.95. Derivable from n and Var(atanh r̂) ≈ 1/(n−3) alone — no result required.
    #
    # DEFECT 2 (extreme-value instability). `max|Δ logBF|` is a MAXIMUM: E[max] is non-decreasing
    # in n for any non-degenerate distribution. A fixed tolerance on a maximum therefore silently
    # TIGHTENS whenever n is raised — which this amendment must do — for reasons unrelated to
    # agreement quality. A property of `maximum` as an estimator; again no result required.
    #
    # PRE-REGISTERED PRECISION TARGET: the two-sided 95 % Fisher-z CI for ρ, AT ρ = BF_CORR_MIN,
    # must have half-width ≤ 0.03 on the correlation scale — 3× finer than the [0.90,1.00] region
    # the decision lives in. DERIVATION (the LOWER limb binds, atanh compresses upward):
    #     atanh(0.95) = ½ln(1.95/0.05) = ½ln 39 = 1.831781
    #     atanh(0.92) = ½ln(1.92/0.08) = ½ln 24 = 1.589027
    #     Δz = 0.242754 ;  √(n−3) ≥ 1.959964/0.242754 = 8.07386
    #     n − 3 ≥ 65.19 ⇒ n ≥ 68.19 ⇒ n_min = 69          (= bf_required_n(0.95, 0.03))
    # Set BF_GATE_N = 100 (n_min rounded up; margin without another doubling of cost). CHECK:
    #     1.959964/√97 = 0.199001 ; z ∈ [1.632780, 2.030782] ; r ∈ [0.92646, 0.96613]
    #     half-widths 0.02354 / 0.01613 — both ≤ 0.03 ✓
    # The target is stated AT THE THRESHOLD, where discrimination is required; it is not claimed
    # to hold uniformly (at ρ = 0.90 the lower half-width at n = 100 is ≈0.045).
    const BF_CORR_MIN    = 0.95        # UNCHANGED from v1 — the threshold value is NOT re-tuned
    const BF_GATE_N      = 100         # paired draws ATTEMPTED for the correlation estimate
    const BF_GATE_N_MIN  = 69          # derived minimum; fewer FINITE pairs ⇒ arm is :invalid,
                                       # never :pass (closes the silent-attrition failure mode)
    const BF_CI_HALFWIDTH_TARGET = 0.03    # the pre-registered precision target n was derived from
    const BF_CI_LEVEL            = 0.95    # one-sided lower-confidence-bound level

    # DECISION RULE (CI-based, one-sided, THREE-WAY):
    #     LB = tanh(atanh(r̂) − 1.644854/√(n−3)) ;  UB = tanh(atanh(r̂) + 1.644854/√(n−3))
    #     :pass         if LB ≥ BF_CORR_MIN
    #     :fail         if UB <  BF_CORR_MIN
    #     :inconclusive otherwise            (corr_pass is TRUE only for :pass)
    # WHY. A ship gate must require EVIDENCE FOR the claim "ρ ≥ 0.95", not merely failure to
    # refute it; the burden of proof belongs on the thing being shipped, and a lower confidence
    # bound is the standard formalisation. NOTE THE DIRECTION: at n = 100 this demands
    #     r̂ ≥ tanh(1.831781 + 1.644854/√97) = tanh(1.998797) = 0.96394
    # versus r̂ ≥ 0.95 in v1 — STRICTLY HARDER. An amendment written to manufacture a pass would
    # not have done this. `:inconclusive` is an honesty gain, not a leniency gain: `passed` is
    # false for it, exactly as for `:fail`. REJECTED ALTERNATIVE: a point rule at larger n would
    # be MORE permissive in expectation while inheriting v1's silence about estimation error.
    const BF_DECISION_RULE = :lower_confidence_bound

    # `max|Δ logBF|`: KEPT, restated on an n-stable estimand, maximum still REPORTED.
    #   GATING   : quantile(|Δ logBF|, BF_LOGBF_Q) ≤ BF_LOGBF_TOL
    #   REPORTED : max|Δ logBF| (plus q50), so the tail stays inspectable and auditable
    # The TOLERANCE VALUE is NOT touched (0.5 nats ≈ one Jeffreys evidence category); only the
    # estimand changes, from an n-dependent extreme to an n-stable order statistic. Second,
    # independent reason for the quantile: the gate already scores against the NON-CLAMPED
    # baseline `kde_log_bf_unclamped` (introduced precisely to remove the 1e-8 clamp-tail
    # artifact), but a KDE log-density evaluated in the far tail of a finite sample is
    # numerically unstable REGARDLESS of clamping — the BASELINE's own tail variance is large
    # there, so one tail draw can dominate a maximum without indicating any disagreement in the
    # regime that matters. That is an argument about the baseline estimator's definition, not
    # about any observed value. The maximum is reported so the artifact stays visible rather than
    # being defined away.
    const BF_LOGBF_TOL = 0.5           # UNCHANGED from v1 — NOT re-tuned
    const BF_LOGBF_Q   = 0.95          # gated quantile of |Δ logBF| (95th of 100 order stats)

    # Δρ SWEEP (a separate REPORTING product, carried over verbatim). v1 overloaded BF_SWEEP_N as
    # the gate's n (`bf_gate(...; n = BF_SWEEP_N)`) — that overloading is the mechanism by which
    # the sample size was never chosen on statistical grounds. BF_GATE_N now governs the gate.
    const BF_SWEEP_LO = -0.6
    const BF_SWEEP_HI = 0.8
    const BF_SWEEP_N  = 25

    # =====================================================================================
    # OOD / misspecification flag — CARRIED OVER VERBATIM (this arm PASSED and is NOT amended)
    # =====================================================================================
    # Re-run only because the F5 image-size mixture changes its input distribution.
    const OOD_ID_QUANTILE = 0.95       # ID score quantile = operating point (~5% ID FPR)
    const OOD_KS_EPS      = 0.05       # negative-control KS-invariance statistic bound
    const OOD_AUC_MIN     = 0.80       # positive-control separability floor (ROC AUC)
    const OOD_GRID_LEVELS = 4          # misspec magnitudes per family
    const OOD_PP_REPS     = 50         # posterior-predictive mismatch replicate count

    # =====================================================================================
    # FORBIDDEN seeds (anti-snooping, Pitfall 3) — the amended gate must reuse NONE of them
    # =====================================================================================
    const NPE_MASTER_SEED     = 0x0000_0000_00C0_FFEE   # spike TRAINING stream    (forbidden)
    const VAL_MASTER_SEED     = 0x0000_0000_5BC0_FFEE   # spike VALIDATION stream  (forbidden)
    const DEFAULT_MASTER_SEED = 0x0000_0000_0000_0001   # productionization DATAGEN (forbidden)
    # DEV seeds burned by the F2 bounded-θ remedy DIAGNOSTIC (07-CALIBRATION-FINDINGS):
    const DEV_SEEDS = (0x0000_0000_0DE7_C0DE, 0x0000_0000_DE7C_0DE2)

    # --- burned streams the v2 forbidden set predates, added for Phase 15 (D-11) --------------
    # Each literal below was VERIFIED against its defining site in this checkout before being
    # written here; none is copied from a summary document. A comment naming a seed is not
    # evidence, which is why the two DERIVED gate families further down are recomputed instead.
    const VAL_FIX_SEED       = 0x0000_0000_00F1_F7ED  # spike FIXTURE stream (spike/validation/consts.jl:79)
    const RATIO_PAIR_SEED    = 0x0000_0000_004A_7107  # shipped ratio pairing (src/amortized/train_ratio.jl:60)
    const CORPUS_MASTER_SEED = 0x0000_0000_00C0_5EED  # provenance-manifest stream (corpus/manifest.csv)
    const P11_DEV_SEED       = 0x0000_0000_0B11_DE71  # spike/validation/p11_consts.jl
    const P12_DEV_SEED       = 0x0000_0000_0B12_DE71  # spike/validation/p12_consts.jl:251
    const P12_FIX_SEED       = 0x0000_0000_0B12_F1F7  # spike/validation/p12_consts.jl:252
    const P13_DEV_SEED       = 0x0000_0000_0B13_DE71  # spike/p13/consts.jl:197
    const P13_FIX_SEED       = 0x0000_0000_0B13_F1F7  # spike/p13/consts.jl:198

    # v1 mixing constants, replicated verbatim so the v1 PROD_SEED[G] values can be RECOMPUTED
    # here and added to the forbidden set (the amended seed must be disjoint from every reported
    # v1 gate stream, not merely from the spike streams).
    const PROD_SALT   = 0x94D0_49BB_1331_11EB
    const PROD_MASTER = 0x0000_0000_09E3_779B

    """
        _derive_prod_seed(G) -> UInt64

    The v1 per-grid gate seed rule, replicated VERBATIM from `gate_consts_template.jl` so the v1
    seeds can be recomputed and treated as forbidden. Not used to seed anything here.
    """
    function _derive_prod_seed(G::Integer)
        rng = Philox4x(UInt64, (UInt64(PROD_MASTER) ⊻ PROD_SALT, UInt64(G)))
        s = rand(rng, UInt64)
        while s in (VAL_MASTER_SEED, NPE_MASTER_SEED, UInt64(0))
            s = rand(rng, UInt64)
        end
        return s
    end

    "The v1 (original pre-registration) gate seeds — recomputed here purely to forbid them."
    const PROD_SEED = Dict{Int,UInt64}(G => _derive_prod_seed(G) for G in (4, 8, 16, 32))

    # =====================================================================================
    # FRESH amended gate stream (same salted-Philox construction, DISTINCT mixing constants)
    # =====================================================================================
    # AMEND_SALT is a Murmur3 finalizer constant, distinct from PROD_SALT (0x94D049BB133111EB),
    # the spike VAL_SALT (0xBF58476D1CE4E5B9), HOLDOUT_SALT (0x9E3779B97F4A7C15) and FOLD_SALT
    # (0xD1B54A32D192ED03). A distinct (master, salt) pair gives a different Philox key, hence a
    # different stream; the redraw loop below then makes disjointness from every named seed
    # explicit rather than merely overwhelmingly probable.
    const AMEND_SALT   = 0xC4CE_B9FE_1A85_EC53
    const AMEND_MASTER = 0x0000_0000_C0DE_2026

    """
        _forbidden_seeds() -> Tuple

    Every seed the amended gate stream must not be: 0, the spike training/validation streams, the
    productionization datagen stream, the two F2-remedy DEV seeds, and ALL FOUR v1 `PROD_SEED[G]`.
    """
    _forbidden_seeds() = (UInt64(0), NPE_MASTER_SEED, VAL_MASTER_SEED, DEFAULT_MASTER_SEED,
                          DEV_SEEDS..., values(PROD_SEED)...)

    """
        _derive_prod_seed_v2(G) -> UInt64

    The FRESH amended per-grid gate seed: a full-width Philox draw keyed on
    `(AMEND_MASTER ⊻ AMEND_SALT, G)`, re-drawn until it is not in `_forbidden_seeds()`. Same
    construction as the template's `_derive_prod_seed`, distinct mixing constants, strictly larger
    forbidden set.
    """
    function _derive_prod_seed_v2(G::Integer)
        rng = Philox4x(UInt64, (UInt64(AMEND_MASTER) ⊻ AMEND_SALT, UInt64(G)))
        s = rand(rng, UInt64)
        while s in _forbidden_seeds()
            s = rand(rng, UInt64)
        end
        return s
    end

    """
        PROD_SEED_V2 :: Dict{Int,UInt64}

    The amended gate seeds. `PROD_SEED_V2[8]` was THE seed for the single confirmation run of the
    amended grid-8 gate (protocol §6 of the amendment: ONE run, no iteration). The sibling grids
    are carried for reference; grids 4 and 16 were NOT re-gated by that amendment.

    **SPENT.** That one run has happened: it produced
    `artifacts/amended_v2/grid_8/gate_report_8.jld2`, whose `report.seed` is `PROD_SEED_V2[8]`.
    Recomputed here for exactly one purpose — to be FORBIDDEN in `_p15_forbidden()` below.
    """
    const PROD_SEED_V2 = Dict{Int,UInt64}(G => _derive_prod_seed_v2(G) for G in (4, 8, 16, 32))

    "The amended gate seed for grid `G` (derived on the fly for a novel grid)."
    prod_seed_v2(G::Integer) = get(() -> _derive_prod_seed_v2(G), PROD_SEED_V2, Int(G))

    "The amended per-grid ship-gate RNG: `Philox4x` keyed by `(prod_seed_v2(G), 0)`."
    prod_rng_v2(G::Integer) = Philox4x(UInt64, (prod_seed_v2(G), UInt64(0)))

    # =====================================================================================
    # FRESH PHASE-15 stream (D-11) — same construction again, a THIRD distinct (master, salt)
    # =====================================================================================
    # WHY A NEW KEY AND NOT A NEW INDEX RANGE. CONVENTIONS C-04: an index-keyed pool drawn at the
    # same Philox key is a SUPERSET of what was already drawn there, not a held-out set — a
    # "fresh range" on `PROD_SEED_V2[8]` would silently re-derive numbers the spent confirmation
    # run already observed. Only a DISTINCT (master, salt) pair changes the Philox KEY, and only a
    # different key gives structural (rather than probabilistic) immunity. The redraw loop below
    # then makes disjointness from every NAMED seed explicit rather than merely overwhelmingly
    # probable.

    """
        P15_REPO_SALTS :: NTuple{11,UInt64}

    Every salt already in use anywhere in this repository. `P15_SALT` must be none of them: a
    repeated (master ⊻ salt) pair is a repeated Philox key, hence a repeated stream, no matter how
    fresh the seed derived from it looks. Inventory verified at `spike/validation/p12_consts.jl:234-245`.
    """
    const P15_REPO_SALTS = (
        0x9E37_79B9_7F4A_7C15,   # HOLDOUT_SALT           (spike/data/seeding.jl)
        0xD1B5_4A32_D192_ED03,   # FOLD / CBS_SPLIT_SALT  (spike/data/seeding.jl)
        0xBF58_476D_1CE4_E5B9,   # VAL_SALT               (spike validation stream)
        0x94D0_49BB_1331_11EB,   # PROD_SALT              (v1 ship gate, above)
        0xC4CE_B9FE_1A85_EC53,   # AMEND_SALT             (v2 amended ship gate, above)
        0xA24B_AED4_663E_E121,   # P11_SALT               (spike/validation/p11_consts.jl)
        0xA5A5_A5A5_A5A5_A5A5,   # COSTES                 (randomization salt)
        0xA5A5_A5A5_DEAD_BEEF,   # COSTES, second word
        0x2545_F491_4F6C_DD1D,   # P13_SALT == P11_DATAGEN_SALT — ONE literal, TWO phases
        0xFF51_AFD7_ED55_8CCD,   # P12_SALT               (spike/validation/p12_consts.jl:255)
        0xC2B2_AE3D_27D4_EB4F,   # P12_DATAGEN_SALT       (spike/validation/p12_consts.jl:256)
    )

    # A degski64 finalizer constant. Chosen for the same reason PROD_SALT (SplitMix64) and
    # AMEND_SALT (Murmur3) were: a published high-avalanche mixing word, from a DIFFERENT hash
    # family than any entry of P15_REPO_SALTS, so it is not a near-neighbour of an existing key.
    # Asserted below and independently in the testset to be absent from P15_REPO_SALTS.
    const P15_SALT   = 0x8EBC_6AF0_9C88_C6E3
    const P15_MASTER = 0x0000_0000_0B15_2026   # "0B15" = phase 15, "2026" = the milestone year

    """
        _p15_forbidden() -> Tuple

    Every seed a Phase-15 stream must NOT be: zero, the six named reserved streams, the two
    F2-remedy DEV seeds, Phase 11/12/13's DEV and FIXTURE streams, and ALL EIGHT ship-gate seeds of
    BOTH frozen families — which are DERIVED above, never literal, and so enter this set from their
    RECOMPUTED values. `PROD_SEED_V2[8]` in particular is the SPENT amended-gate seed the Phase-15
    domain constraint names by hand.
    """
    _p15_forbidden() = (UInt64(0), NPE_MASTER_SEED, VAL_MASTER_SEED, VAL_FIX_SEED,
                        DEFAULT_MASTER_SEED, RATIO_PAIR_SEED, CORPUS_MASTER_SEED, DEV_SEEDS...,
                        P11_DEV_SEED, P12_DEV_SEED, P12_FIX_SEED, P13_DEV_SEED, P13_FIX_SEED,
                        values(PROD_SEED)..., values(PROD_SEED_V2)...)

    """
        _derive_prod_seed_p15(G) -> UInt64

    The FRESH Phase-15 per-grid gate seed: a full-width Philox draw keyed on
    `(P15_MASTER ⊻ P15_SALT, G)`, re-drawn until it is not in `_p15_forbidden()`. Same construction
    as `_derive_prod_seed` / `_derive_prod_seed_v2` above (deliberately — the rule is the frozen
    part), distinct mixing constants, strictly larger forbidden set.
    """
    function _derive_prod_seed_p15(G::Integer)
        rng = Philox4x(UInt64, (UInt64(P15_MASTER) ⊻ P15_SALT, UInt64(G)))
        s = rand(rng, UInt64)
        while s in _p15_forbidden()
            s = rand(rng, UInt64)
        end
        return s
    end

    """
        PROD_SEED_P15 :: Dict{Int,UInt64}

    The Phase-15 gate seeds. `PROD_SEED_P15[8]` is THE seed every reported Phase-15 number rides:
    the envelope sweep, the OOD-crossing arm and the rung-0 ECE anchor. The sibling grids are
    carried for symmetry with the two frozen families; D-10a drops the `grid_16` contrast, so no
    grid other than 8 is swept by this phase.
    """
    const PROD_SEED_P15 = Dict{Int,UInt64}(G => _derive_prod_seed_p15(G) for G in (4, 8, 16, 32))

    "The Phase-15 gate seed for grid `G` (derived on the fly for a novel grid)."
    prod_seed_p15(G::Integer) = get(() -> _derive_prod_seed_p15(G), PROD_SEED_P15, Int(G))

    "The Phase-15 per-grid gate RNG: `Philox4x` keyed by `(prod_seed_p15(G), 0)`."
    prod_rng_p15(G::Integer) = Philox4x(UInt64, (prod_seed_p15(G), UInt64(0)))

    # Phase 15 consumes the P15 stream. These aliases (mirroring `gate_consts_8_v2.jl:342-343`)
    # exist so the EXISTING gate machinery — `sbc_gate`, `bf_gate`, `ood_gate`, all of which
    # default to `rng = prod_rng(G)` — rides the FRESH stream when this file is the loaded
    # pre-registration, with no change to any of them.
    #
    # THEY ALSO CARRY THE GLOBAL HALF OF REPRODUCIBILITY. `harness.jl:82` defines
    # `gate_global_seed(G) = rand(Philox4x(UInt64, (prod_seed(G), UInt64(1))), UInt64)` and
    # `seed_gate_global!` pins `Random.default_rng()` to it, because `sampleposterior` threads no
    # `rng` and draws its flow base samples GLOBALLY (`harness.jl:49-60`). Aliasing `prod_seed`
    # here is therefore what makes BOTH halves of every Phase-15 arm ride P15 rather than the spent
    # v2 stream — and the testset asserts the two global seeds differ, so this is checked, not
    # assumed.
    prod_seed(G::Integer) = prod_seed_p15(G)
    prod_rng(G::Integer)  = prod_rng_p15(G)

    # =====================================================================================
    # Executable form of the amended decision rules (pure functions; nothing is run here)
    # =====================================================================================

    # Standard normal quantiles, hard-coded so this pre-registration file carries no statistical
    # dependency beyond Random123.
    const Z_TWO_SIDED_95 = 1.959963985   # Φ⁻¹(0.975)
    const Z_ONE_SIDED_95 = 1.644853627   # Φ⁻¹(0.950)

    """
        holm_adjusted(p) -> Vector{Float64}

    Holm–Bonferroni adjusted p-values for a family of `m = length(p)` tests, returned in the INPUT
    order. Sort ascending, multiply p₍ⱼ₎ by `(m − j + 1)`, take the running maximum (monotonicity
    enforcement) and cap at 1. Valid under arbitrary dependence between the tests.
    """
    function holm_adjusted(p::AbstractVector{<:Real})
        m   = length(p)
        ord = sortperm(collect(float.(p)))
        adj = Vector{Float64}(undef, m)
        running = 0.0
        for (j, i) in enumerate(ord)
            running = max(running, (m - j + 1) * float(p[i]))
            adj[i]  = min(1.0, running)
        end
        return adj
    end

    """
        holm_pass(p; fwer = 0.05) -> Bool

    `true` iff the Holm–Bonferroni procedure at family-wise error rate `fwer` rejects NO null in
    the family — i.e. every adjusted p-value exceeds `fwer`. This is the amended SBC uniformity
    verdict (A1), applied separately to the 8 KS and the 8 χ² p-values.
    """
    holm_pass(p::AbstractVector{<:Real}; fwer::Real = SBC_KS_FWER) =
        all(>(fwer), holm_adjusted(p))

    """
        fisher_z_ci(r, n; level = 0.95, tail = :two) -> (lo, hi)

    Fisher-z confidence interval for a Pearson correlation: `tanh(atanh(r) ± z·/√(n−3))`, with
    `z· = Φ⁻¹(0.975)` for `tail = :two` and `Φ⁻¹(0.95)` for `tail = :one`. Returns `(-1.0, 1.0)`
    (i.e. no information) for `n ≤ 3` or a non-finite `r`.
    """
    function fisher_z_ci(r::Real, n::Integer; level::Real = BF_CI_LEVEL, tail::Symbol = :two)
        (isfinite(r) && n > 3) || return (lo = -1.0, hi = 1.0)
        level == BF_CI_LEVEL || error("fisher_z_ci: only the pre-registered level $BF_CI_LEVEL " *
                                      "is admissible in a frozen gate")
        zc = tail === :two ? Z_TWO_SIDED_95 : Z_ONE_SIDED_95
        z  = atanh(clamp(float(r), -0.999999, 0.999999))
        d  = zc / sqrt(n - 3)
        return (lo = tanh(z - d), hi = tanh(z + d))
    end

    """
        bf_required_n(r0 = BF_CORR_MIN, halfwidth = BF_CI_HALFWIDTH_TARGET) -> Int

    The sample size at which the two-sided 95 % Fisher-z CI for a correlation of `r0` has a
    half-width ≤ `halfwidth` on the correlation scale. The LOWER limb binds (atanh compresses
    upward), so `n = ceil(3 + (Φ⁻¹(0.975) / (atanh(r0) − atanh(r0 − halfwidth)))²)`. At the
    pre-registered target this returns **69**, the derivation behind `BF_GATE_N_MIN`.
    """
    function bf_required_n(r0::Real = BF_CORR_MIN, halfwidth::Real = BF_CI_HALFWIDTH_TARGET)
        d = atanh(float(r0)) - atanh(float(r0) - float(halfwidth))
        return ceil(Int, 3 + (Z_TWO_SIDED_95 / d)^2)
    end

    """
        bf_corr_verdict(r, n; thr = BF_CORR_MIN) -> Symbol

    The amended three-way BF correlation verdict (A2): `:pass` iff the ONE-SIDED 95 % lower
    confidence bound on ρ is ≥ `thr`; `:fail` iff the corresponding upper bound is < `thr`;
    `:inconclusive` when the interval straddles `thr`. `:invalid` if fewer than `BF_GATE_N_MIN`
    finite pairs survived — attrition must never be able to return the gate to the underpowered
    regime. Only `:pass` sets `corr_pass`.
    """
    function bf_corr_verdict(r::Real, n::Integer; thr::Real = BF_CORR_MIN)
        n >= BF_GATE_N_MIN || return :invalid
        isfinite(r) || return :invalid
        ci = fisher_z_ci(r, n; tail = :one)
        ci.lo >= thr && return :pass
        ci.hi <  thr && return :fail
        return :inconclusive
    end

    # #####################################################################################
    # #                             PHASE-15 TIER-1 BARS                                  #
    # #                                                                                   #
    # # Everything from here to the self-check block is frozen BEFORE any Phase-15 rung    #
    # # exists. EVERY bar carries its derivation in the comment adjacent to it. A bar with #
    # # no derivation is the precise liability the Stage-1-ceiling adjudication cost this  #
    # # project, and the licensing standard in this file's header applies to all of them.  #
    # #####################################################################################

    # =====================================================================================
    # Sweep design (D-03, and the research compute budget)
    # =====================================================================================
    #
    # M IS A PRECISION PROPERTY HERE, NOT A POWER LEVER. `gate_consts_8_v2.jl:66-68` froze
    # SBC_M = 2000 explicitly because "changing M would change power in a direction the author knew
    # he needed" — a concern about a SIGNIFICANCE gate. D-04 deliberately moved the break statistic
    # off significance and onto an EFFECT SIZE (ECE), which retires exactly that objection: ECE's
    # dependence on M is resolution, and it is disclosed rather than exploited.
    #
    # DERIVATION of M = 1000, from the measured ECE null (P15_ECE_NULL_* below): the smallest true
    # ECE degradation the ladder can resolve is ≈ the margin — 0.00893 at M = 1000 versus 0.0067 at
    # M = 2000 for TWICE the compute. Against the shipped ρ_true ECE of 0.03713 a resolution of
    # 0.00893 detects a ~24 % relative degradation, which is the granularity a domain-of-
    # applicability map is read at. Below M = 500 the resolution (0.013+) approaches the shipped
    # anchor's own magnitude and the ladder stops discriminating at all.
    const P15_SWEEP_M = 1000            # SBC draws per rung (resolution, not power — see above)
    const P15_RUNGS   = 5               # rungs per axis; endpoints FROZEN per D-13, spacing is
                                        # Claude's discretion per the CONTEXT ruling
    const P15_RUNG0   = 0               # the shared IN-PRIOR anchor rung: no misspecification, no
                                        # extrapolation. Measured ONCE per net at P15_SWEEP_M on
                                        # the P15 stream and appended under :P15_ECE_ANCHOR_MEASURED

    # DERIVATION: `07-GATE-AMENDMENT.md:140-149` measures E[cost] = 0.08346 s/pair over the F5
    # mixture AT 32 THREADS, and the research cost model puts the single-threaded arm at ~13×
    # slower. A single-threaded sweep is therefore a budget ERROR, not a style choice, and
    # JULIA_NUM_THREADS must be pinned in every runner and CI job that touches this phase.
    const P15_SWEEP_MIN_THREADS = 16    # below this the sweep cannot finish inside its own ceiling

    # DERIVATION: the reconciled cost model (which reproduces the ONE observed gate wall-clock to
    # within 5 % — predicted 29.4 min vs 30.9 min measured) predicts ~8.5 h at 32 threads for
    # 1 rung-0 + SEVEN axes × 5 rungs at M = 1000. The ceiling carries ~1.41× headroom, which is
    # there for research assumption A2: misspecification generators are taken to average ~1.3×
    # `simulate_pair`, and that factor is ESTIMATED, not measured.
    #
    # THE BASIS OF THIS NUMBER, STATED EXACTLY. 8.5 h is the SEVEN-axis figure. `15-RESEARCH.md`'s
    # 7.3 h is a SIX-axis estimate computed before D-02a added the seventh (pure-offset
    # autofluorescence) axis; quoting 7.3 h / 1.6× here would be a bar whose stated justification
    # does not match its actual basis, which is the exact failure this project has already had to
    # correct four times.
    const P15_SWEEP_WALLCLOCK_CEILING_HOURS = 12.0   # ~1.41× over the 7-axis ~8.5 h prediction

    # =====================================================================================
    # The break criterion (D-04, D-04a) — ECE, anchored at rung 0, margin from the null
    # =====================================================================================
    #
    # WHAT ECE IS HERE, EXACTLY (derived from `sbc.jl:302-315` + `_bin_calibration`, and verified
    # numerically against the shipped report to the last float digit):
    #     ECE = (1/19) · Σ_{α ∈ {0.05, 0.10, …, 0.95}} |α − ĉ(α)|,   ĉ(α) = empirical central-
    # interval coverage at nominal level α from u = (rank + 0.5)/(L + 1). Each α lands in its own
    # reliability bin, so the statistic does not depend on SBC_BINS at all.
    #
    # ECE HAS A POSITIVE EXPECTATION EVEN UNDER PERFECT CALIBRATION — it is a mean of 19 absolute
    # deviations of binomial proportions. That null is computable with NO net, NO simulation and NO
    # data: draw M uniform ranks, evaluate the identity above, repeat. Hence the design below.
    const P15_ECE_NULL_B    = 20_000     # Monte-Carlo replicates behind the two bars below
    const P15_ECE_NULL_SEED = 20260804   # Philox4x(UInt64, (seed, 0)); the whole measurement is
                                         # reproducible from these two numbers plus M and L alone

    # MEASURED under exactly that design (B = P15_ECE_NULL_B, seed = P15_ECE_NULL_SEED,
    # M = P15_SWEEP_M, L = SBC_L, uniform ranks on 0:L), rounded to 5 dp:
    const P15_ECE_NULL_MEAN = 0.01033    # E[ECE₀]; closed form 0.336/√M = 0.010625 agrees to 2.8 %
    const P15_ECE_NULL_Q95  = 0.01926    # q95(ECE₀) at the same M

    # INDEPENDENT CROSS-CHECK, carried so the pinned values are auditable against a second
    # measurement rather than only against themselves: `15-RESEARCH.md` §"The ECE Null
    # Distribution" measured the same quantities with B = 2000 (a different, unrecorded seed) and
    # obtained 0.01038 / 0.01884. The mean agrees to 5e-5; the q95 differs by 0.0004, which is
    # ~0.6 SE of a 95th percentile estimated from only 2000 replicates. The 20 000-replicate values
    # are pinned because they are the ones the RECORDED design reproduces — a bar must be
    # regenerable from what the file says produced it.
    const P15_ECE_NULL_MEAN_RESEARCH = 0.01038   # 15-RESEARCH.md, B = 2000
    const P15_ECE_NULL_Q95_RESEARCH  = 0.01884   # 15-RESEARCH.md, B = 2000

    # THE MARGIN. Derivation: `q95(ECE₀) − E[ECE₀]` at the sweep's own M. Both terms are arithmetic
    # on the TEST DESIGN; neither is an outcome. It is an EFFECT SIZE, so unlike a p-value it does
    # not grow teeth as M rises.
    #
    # OPERATING CHARACTERISTIC, STATED UP FRONT. q95 was chosen over q99 KNOWING the cost: a
    # genuinely-still-calibrated rung crosses it ~5 % of the time, i.e. ≈1.75 expected FALSE breaks
    # across the 35 axis × rung tests (0.05 × 7 axes × 5 rungs). THIS IS WHY THE DOMAIN MAP MUST BE
    # READ PER AXIS AND NEVER PER RUNG, and that instruction has to appear wherever the map is
    # reported. (1.75 / 35 are the SEVEN-axis figures; ~1.5 / ~30 in earlier drafts are the stale
    # six-axis numbers from before D-02a.)
    const P15_ECE_MARGIN = 0.00893       # = P15_ECE_NULL_Q95 − P15_ECE_NULL_MEAN (asserted below)

    """
        p15_break_threshold(ece_rung0) -> Float64

    The Phase-15 break threshold for one column: the rung-0 ECE plus `P15_ECE_MARGIN`.

    THE ANCHOR IS A TIER-2 MEASUREMENT, not a shipped number (D-04a). It is `ECE(rung 0)` measured
    at `P15_SWEEP_M` on the P15 stream for the net under test, appended under the
    `:P15_ECE_ANCHOR_MEASURED` sentinel by 15-07.

    WHY THE SHIPPED `gate_report_8.jld2` ECEs (ρ_true 0.03713, Δρ 0.00429) ARE NOT THE ANCHOR:
    Δρ's shipped value sits BELOW its own null mean of 0.00751 at M = 2000. Anchoring there would
    declare a break on ≈90 % of PERFECTLY CALIBRATED rungs — D-04's own failure mode ("a criterion
    whose pass/fail tracks something other than what it was meant to measure") reappearing inside
    the anchor. Anchoring at rung 0 at the sweep's own M makes the null's noise floor COMMON-MODE
    between anchor and rung, so it cancels instead of accumulating.
    """
    p15_break_threshold(ece_rung0::Real) = float(ece_rung0) + P15_ECE_MARGIN

    # THE ENVELOPE IS DEFINED BY TWO COLUMNS ONLY (D-05). `SBC_PARAM_LABELS` (`sbc.jl:123-124`)
    # orders the columns 1 = ρ_true, 2 = spillover, 3 = autofluorescence, 4 = label_efficiency,
    # 5 = shift_dx, 6 = shift_dy, 7 = noise, 8 = Δρ. Columns 1 and 8 are the SBC-calibrated targets
    # (randomized-rank atom handling) and the quantities a user acts on.
    const P15_BREAK_COLUMNS = (1, 8)                   # GATING
    # The six nuisance marginals are MEASURED AND REPORTED at every rung and gate NOTHING. Their
    # ~0.08 SD drift is an ALREADY-DISCLOSED named v2.0 limit (`docs/amortized.md`); gating on it
    # would re-fail the tool for a reason already on the record and would collapse the envelope to
    # near-zero on every axis, producing no map at all.
    const P15_NUISANCE_COLUMNS = (2, 3, 4, 5, 6, 7)    # REPORTED, never gating
    const P15_REPORTED_COLUMNS = (1, 2, 3, 4, 5, 6, 7, 8)   # all eight are reported every rung

    # The 19 nominal levels the ECE identity averages over — the same grid `sbc_coverage` defaults
    # to (`sbc.jl`), stated here so the statistic is fully specified by this file.
    const P15_COVERAGE_LEVELS = 0.05:0.05:0.95

    # MCE IS FORBIDDEN IN THIS PHASE (D-04a). Derivation: `_bin_calibration` walks all 50
    # reliability bins including the 31 EMPTY ones, where `predicted_rate` is the bin midpoint and
    # `observed_rate` is 0. The maximum of |predicted − observed| is therefore attained on an empty
    # bin and is pinned at 0.99 on ALL EIGHT columns of the shipped report regardless of the net.
    # It is an artifact of the binning, not a measurement of anything. No Phase-15 code may read
    # it; 15-04 turns this constant into an executable source guard.
    const P15_MCE_FORBIDDEN = true      # reading `.mce` anywhere in Phase 15 is a defect

    # =====================================================================================
    # The OOD crossing arm (D-06, D-07)
    # =====================================================================================
    #
    # DERIVATION of the sample sizes: `fire_rate` at a rung is a binomial proportion with n = n_pos.
    # At the shipped `n_pos = 40` its standard error at p ≈ 0.05 is √(0.05·0.95/40) = 0.034, so ANY
    # margin below ~0.10 sits inside one SE of the baseline and is not a detection at all. The OOD
    # arm runs NO posterior draws (`gate_ood_roc` computes only `maha_score`/`noise_score`), so
    # raising n to 200 costs ≈20 min for the entire experiment — the cheapest possible fix.
    const P15_OOD_N_ID  = 200           # in-distribution images for the MEASURED baseline
    const P15_OOD_N_FIT = 200           # images the TRAIN-only null is fitted on
    const P15_OOD_N_POS = 200           # positives per axis × rung (the binomial n above)

    """
        _p15_binom_quantile(n, p, q) -> Int

    The smallest `k` with `P(Binomial(n, p) ≤ k) ≥ q`, computed by accumulating the exact binomial
    PMF through the recurrence `P(k) = P(k-1)·((n-k+1)/k)·(p/(1-p))` from `P(0) = (1-p)^n`.

    Hand-rolled ON PURPOSE: this pre-registration carries no statistical dependency beyond Random123
    (the same reason `Z_TWO_SIDED_95` is a literal above), and the point of the helper is that
    `P15_OOD_FIRE_MARGIN` is DERIVED at load time rather than typed as a number someone chose.
    """
    function _p15_binom_quantile(n::Integer, p::Real, q::Real)
        pf, qf = float(p), float(q)
        (0.0 < pf < 1.0) || error("_p15_binom_quantile: p must be in (0, 1), got $p")
        (0.0 < qf < 1.0) || error("_p15_binom_quantile: q must be in (0, 1), got $q")
        pk  = (1 - pf)^n                 # P(X = 0)
        cdf = pk
        k   = 0
        while cdf < qf && k < n
            k  += 1
            pk *= (n - k + 1) / k * pf / (1 - pf)
            cdf += pk
        end
        return k
    end

    # DERIVATION: the one-sided binomial 95th percentile of the ID fire count at n = P15_OOD_N_POS
    # under the design ID rate `1 − OOD_ID_QUANTILE = 0.05`, expressed as a rate and taken relative
    # to that same rate:  15/200 − 0.05 = 0.025.  Asserted below via `_p15_binom_quantile`, so the
    # bar is reproduced from the design rather than trusted as a literal.
    #
    # THE RULE COMPARES AGAINST THE MEASURED BASELINE, NEVER THE DESIGN VALUE (D-07):
    #     L_ood(axis) = first rung r with fire_rate[axis][r] > id_fire_rate + P15_OOD_FIRE_MARGIN
    # where `id_fire_rate` is the value `gate_ood_roc` MEASURES on this net on this stream. Using
    # `1 − OOD_ID_QUANTILE` instead would score the detector against its nominal operating point
    # rather than its realised one — the same class of error as anchoring on a shipped number.
    const P15_OOD_FIRE_MARGIN = 0.025   # = q95(Binom(200, 0.05))/200 − 0.05 (asserted below)

    # =====================================================================================
    # Ladders and mechanism typing (D-02, D-02a, D-03)
    # =====================================================================================
    #
    # THE MECHANISM TYPING IS THE FINDING, NOT BOOKKEEPING. An IN-PRIOR axis that breaks coverage
    # is a TRAINING failure — the net saw those θ and still miscalibrates. An OUT-OF-MODEL axis
    # that breaks is a SCOPE limit — expected, honest, and exactly what the OOD flag exists for.
    # One number covering both would produce a domain map nobody can act on.

    # The simulator prior extrema, in the units each axis's ladder moves, cited per entry from
    # `src/amortized/simulator.jl:92-95`. These are LITERALS so this file stays dependency-free;
    # `test_p15_consts.jl` asserts each equals `maximum(ProteinCoLoc.<PRIOR>)`, so the literal
    # cannot drift away from the simulator without a test failing.
    const P15_PRIOR_BOUNDARY = (
        spillover        = 0.2,    # maximum(SPILLOVER_PRIOR)        = Uniform(0.0, 0.2)
        registration     = 1.0,    # maximum(SHIFT_PRIOR)            = Uniform(-1.0, 1.0), px
        autofluorescence = 0.1,    # maximum(AUTOFLUORESCENCE_PRIOR) = Uniform(0.0, 0.1)
    )

    # RUNG 2 IS THE PRIOR BOUNDARY on every in-prior axis (D-03): a break INSIDE prior support is a
    # materially different — and more serious — result than a break outside it, and a ladder that
    # starts at the edge cannot tell them apart. Rung 1 therefore sits strictly inside.
    #
    # ENDPOINT DERIVATION, spillover: rung 5 = 1.0 is the SIMULATOR'S OWN guard
    # (`simulator.jl:232`, `0 ≤ spillover ≤ 1`). The endpoint is the model's limit, not a choice —
    # the ladder cannot be extended without the forward model throwing.
    const P15_SPILLOVER_LADDER = (0.1, 0.2, 0.35, 0.6, 1.0)      # θ.spillover override

    # ENDPOINT DERIVATION, registration: rung 5 = 8 px is 8× the prior boundary and 2.7× the widest
    # shift Phase 11 researched (3 px). Beyond it the vacated border fill (`BG_FLOOR`,
    # `simulator.jl:183`) occupies enough of the frame that the image is dominated by the fill
    # rather than by the registration error, so a larger rung would measure the fill.
    const P15_SHIFT_LADDER = (0.5, 1.0, 2.0, 4.0, 8.0)           # |(shift_dx, shift_dy)| in px

    # ENDPOINT DERIVATION, autofluorescence: rung 5 = 0.8 is 8× the prior boundary (matching the
    # registration ladder's reach in boundary units) and is the magnitude the shipped compound
    # `background` family already reaches at its level 2 (`bleed = 0.1 + 0.3·level`), so the two
    # axes stay comparable in offset units at their far ends. `simulator.jl:234` guards only
    # `autofluorescence ≥ 0`, so there is no model limit to stop at.
    const P15_AF_LADDER = (0.05, 0.1, 0.2, 0.4, 0.8)             # pure additive offset

    # The four SHIPPED families take their own internal `level` formulas (`misspec.jl:104-184`) and
    # are swept at their native integer levels; there is no prior-support quantity to re-express.
    # HONESTY NOTE FOR THE REPORT: `misspec_texture`'s smoothing saturates at σ = 0.8 for
    # level ≥ 4.4 (`σ = max(0.8, 3.0 − 0.5·level)`), so its rung 5 differs from rung 4 in SPOT
    # COUNT only, not in sharpness. And `misspec_optics` moves σx as well as σy, so that axis is
    # "aberrated/out-of-focus", not pure anisotropy.
    const P15_NATIVE_LEVELS = 1:P15_RUNGS

    # The seven axes in REPORT ORDER: the four shipped families in their fixed Phase-5 order first,
    # then the three in-prior axes.
    const P15_AXES = (:texture, :noise, :optics, :background,
                      :spillover, :registration, :autofluorescence)

    # D-02a's reclassification, recorded WITH ITS EVIDENCE. `background` is OUT-OF-MODEL ONLY: its
    # level-1 bleed is already 0.4 against a prior max of 0.1 (4× the boundary before the ladder
    # starts), AND it composes a multiplicative illumination gradient and a radial vignette the
    # simulator has NO parameter for — so it cannot express "at the prior edge" in any units and
    # D-03 was unsatisfiable on it as written. Adding the seventh, pure-offset `autofluorescence`
    # axis is what makes D-03 hold on every in-prior axis WITHOUT EXCEPTION; that route was chosen
    # over amending D-03 to "where markable" precisely so no exception has to be carried.
    const P15_AXIS_MECHANISM = (
        texture          = :out_of_model,   # puncta mixture; not expressible at any θ
        noise            = :out_of_model,   # spikes + heavy tails beyond Poisson+Gaussian
        optics           = :out_of_model,   # σ_psf = 1.3 is FIXED, not a prior (simulator.jl:179)
        background       = :out_of_model,   # D-02a: gradient + vignette have no θ at all
        spillover        = :in_prior,       # Uniform(0.0, 0.2)
        registration     = :in_prior,       # Uniform(-1.0, 1.0) px, per axis
        autofluorescence = :in_prior,       # Uniform(0.0, 0.1)
    )

    # `nothing` means "THIS AXIS HAS NO PRIOR SUPPORT TO MARK", never "unmeasured". Key order is
    # kept identical to `P15_AXIS_MECHANISM` so the two are comparable by `keys` (asserted below).
    const P15_PRIOR_BOUNDARY_RUNG = (
        texture          = nothing,
        noise            = nothing,
        optics           = nothing,
        background       = nothing,   # D-02a: out-of-model, so no boundary rung is carried
        spillover        = 2,
        registration     = 2,
        autofluorescence = 2,
    )

    # =====================================================================================
    # D-13 declarations — written down BEFORE any rung can exist
    # =====================================================================================
    #
    # WHY THE FREEZE. "Extend until it fails" is searching for the failure you want. The Phase-12
    # ruling refused the symmetric move in as many words: "that would be tuning a model until it
    # passes its own calibration gate — a gate that has then stopped measuring anything."
    const P15_ITERATION_ALLOWANCE = 1   # EXACTLY ONE, and CONSTRAINED: spendable on ladder
                                        # EXTENSION ONLY, and only in the "never breaks within the
                                        # ladder" case, recorded as spent by appending the Tier-2
                                        # :P15_ITERATION_SPENT block. Any other extension is an
                                        # AMENDMENT and must be argued as one, against the same
                                        # licensing standard as this file's header.

    "The pre-declared reading of the EMPTY-envelope degenerate outcome (D-13)."
    const P15_DEGENERATE_EMPTY_ENVELOPE_READING =
        "Breaks at rung 1 on every axis means the tool has no usable operating range at the swept " *
        "resolution. That is a REPORTABLE NEGATIVE about the shipped net's domain of " *
        "applicability, not a phase failure, and it is written up as one — the Phase-11 and " *
        "Phase-12 negative-but-useful closures are the precedent. It does NOT license widening " *
        "the break margin, lowering M, or re-anchoring the threshold."

    "The pre-declared reading of the UNBOUNDED-envelope degenerate outcome (D-13)."
    const P15_DEGENERATE_UNBOUNDED_ENVELOPE_READING =
        "Never breaking within the ladder means the operating envelope EXCEEDS the swept range. " *
        "The honest report is then the BOUND — 'calibration holds at least out to rung 5 on this " *
        "axis' — not a wider search. This is the ONLY case in which P15_ITERATION_ALLOWANCE may " *
        "be spent, and spending it extends the ladder ONCE and is recorded as spent."

    # SC2 IS THREE-VALUED AND DECLARED PER AXIS (D-06), because "the flag never fired" has two
    # OPPOSITE meanings — working correctly where nothing was wrong, versus failing to warn where
    # something was — and a two-valued verdict collapses them. Collapsing them is what would force
    # a fifth after-the-fact amendment. Note the asymmetry that motivates it: the three in-prior
    # axes are inside the distribution the OOD null was FITTED on, so the detector may correctly
    # never fire on them.
    const P15_SC2_STATES = (:protective, :late, :silent_but_safe)

    "The pre-registered SC2 verdict rule (D-06)."
    const P15_SC2_PASS_RULE =
        "Per axis: PROTECTIVE if L_ood <= L_break (the flag fires at or before the rung where " *
        "coverage breaks); LATE if L_ood > L_break (coverage breaks while the flag is still " *
        "quiet — the failure SC2 exists to catch); SILENT-BUT-SAFE if coverage never breaks " *
        "within the swept ladder, so no warning was owed. SC2 PASSES IFF NO AXIS IS LATE."

    "The pre-registered response to a LATE axis (D-14). Written before any axis can be late."
    const P15_LATE_AXIS_RESPONSE =
        "A LATE axis is REPORTED AS A NAMED LIMIT — 'the OOD flag does not protect against X' — " *
        "and the phase CLOSES on it. That is a real, publishable finding on the Phase-11/12 " *
        "negative-but-useful precedent. EXPLICITLY FORBIDDEN in response: lowering " *
        "OOD_ID_QUANTILE until the flag fires early enough (tuning a detector until it passes " *
        "its own gate), and building a new detector channel for the failing mechanism (new " *
        "capability, its own phase). No retuning, no retraining."

    # =====================================================================================
    # D-10a — the dropped grid_16 contrast, carried as a NAMED LIMIT
    # =====================================================================================
    const P15_GRID16_CONTRAST = :dropped

    "The evidence for dropping the grid_16 contrast, and the limit that is carried instead (D-10a)."
    const P15_GRID16_DROP_REASON =
        "artifacts/grid_16/npe_16.jld2's `meta` is (:grid, :n_pairs, :use_gpu) — " *
        "`training_imsize_provenance` is ABSENT ENTIRELY, not merely unrecorded — so " *
        "SBC_REQUIRE_IMSIZE_PROVENANCE = true makes run_gate return :provenance_mismatch and run " *
        "NO arm. It is moreover a v1-era net (256^2, unbounded theta), so bypassing the guard " *
        "would confound grid x imsize-regime x theta-space in a single comparison and the " *
        "resulting number could not answer the grid-dependence question the contrast exists to " *
        "answer. CARRIED LIMIT: grid-dependence of the operating envelope is UNTESTED; the " *
        "contrast was ATTEMPTED and BLOCKED by missing provenance on the only available " *
        "higher-dimension net. This is a named limit, not a silent scope cut."

    # The net that carries the REPORTED envelope (D-10's surviving half), and the artifacts tree it
    # lives in. NOTE: `run_gate`'s DEFAULT root is the v1 `artifacts/` tree, so every Phase-15
    # invocation must pass `--artifacts-root artifacts/amended_v2` explicitly or it will silently
    # gate the wrong net.
    const P15_REPORTED_NET   = "artifacts/amended_v2/grid_8"
    const P15_ARTIFACTS_ROOT = "artifacts/amended_v2"

    # =====================================================================================
    # Frozen-file pins — the executable form of this phase's do-not-touch constraints
    # =====================================================================================
    #
    # NORMALIZATION (load-bearing): the pinned digest is SHA-256 of the file's bytes with CRLF
    # newlines normalized to LF. `.gitattributes` sets `* text=auto`, so the WORKING-TREE bytes of
    # every text file here are CRLF on Windows and LF on Linux. A raw-bytes pin would therefore
    # fail on the ubuntu-latest CI this phase exists to build, for a reason that has nothing to do
    # with the file's content — a gate that fails for the wrong reason gets disabled, which is the
    # outcome D-09 explicitly warns about. Normalizing makes the pin content-addressed in the same
    # sense git's own blob identity is.
    const P15_FROZEN_HASH_NORMALIZATION = :lf   # CRLF -> LF before hashing; see above

    """
        P15_FROZEN_FILE_SHA256 :: NamedTuple

    Repo-relative path + LF-normalized SHA-256 for every file Phase 15 must NOT touch.
    `test/gate/test_p15_consts.jl` recomputes each on every test run. Why each is pinned:

      • `misspec.jl` — appending to `OOD_FAMILIES` would silently REWRITE the shipped, PASSED OOD
        arm's `combined_auc` and its pass conjunction (D-02a). New families must be passed through
        `gate_ood_roc`'s existing `families =` keyword as a local merge, never appended.
      • `sbc.jl`, `harness.jl` — the statistics and the seeded draw→simulate→infer chain the
        reported numbers of three phases already ride.
      • `gate_consts_8_v2.jl` — the frozen amended pre-registration this file forks; if it moved,
        the per-constant inheritance check in testset 4 would be comparing against something else.
      • `Project.toml` / `Manifest.toml` — `test/runtests.jl:32-49` asserts the EXACT strings
        "0.2.1" and "0.16.10" while `[compat]` permits a wider range, so ANY re-resolve is a coin
        flip against the co-resolution hard gate. D-08a's ruling — CI installs via
        `Pkg.instantiate()` from the committed manifest and never resolves — is what these two pins
        make enforceable.
      • `spike/Project.toml` / `spike/Manifest.toml` — `spike/test/test_p12_decoupling.jl:176-178`
        already asserts the spike environment byte-unchanged; pinned here so a Phase-15 change
        cannot break a Phase-12 guard from the other side.
    """
    const P15_FROZEN_FILE_SHA256 = (
        misspec = (path   = "test/gate/misspec.jl",
                   sha256 = "f908dcadba63113438d0a216b0ee9559f82d833ccef90b286e6409fdfd1ecfbe"),
        sbc     = (path   = "test/gate/sbc.jl",
                   sha256 = "ed6dc886b2486ceb97495aa0e43fafcd571c70782351853650761988053dd1aa"),
        harness = (path   = "test/gate/harness.jl",
                   sha256 = "b10bca3c616918cd29b39dfa6bc06126d81bacaecc91fe917d0f192b94b41523"),
        gate_consts_8_v2 = (path   = "test/gate/gate_consts_8_v2.jl",
                            sha256 = "fb8617be363a3d20179c94f2b782cbb1bb4cf9a5dd2f25b01a13bf79a22a00e2"),
        project  = (path   = "Project.toml",
                    sha256 = "abed47c830052669306bdad290c51b94eccbf69fbb6bcee067add51d38304fc0"),
        manifest = (path   = "Manifest.toml",
                    sha256 = "3a75e078ed39b191fdceff8b1b5f81f5a39ebda33c3e6dfa6923aa5bb7500af8"),
        spike_project  = (path   = "spike/Project.toml",
                          sha256 = "06f483b2ec3ef9b4c96133435aaebf0dbab0e927ca69db5b03d6b011be578a8d"),
        spike_manifest = (path   = "spike/Manifest.toml",
                          sha256 = "0b32fa2827439e1ca4e0953360e99d1dd3f0268cdc80b41988ec39e42fad4678"),
    )

    # =====================================================================================
    # Tier split — "reported is not gated" made MACHINE-CHECKABLE (D-04a, D-05)
    # =====================================================================================
    # This split is not left to prose. It is two disjoint constant-name tuples with an emptiness
    # assertion on their intersection, mirroring `p12_consts.jl`'s
    # P12_GATING_CONSTANTS / P12_REPORTING_ONLY_CONSTANTS.
    const P15_KS_REPORTED   = true      # KS p-value computed and reported at EVERY rung, gates nothing
    const P15_CHI2_REPORTED = true      # χ² p-value likewise — D-04 keeps both auditable rather than
                                        # hidden, precisely because they are not the break statistic

    "The constants that DECIDE something. Nothing outside this tuple may change a Phase-15 verdict."
    const P15_GATING_CONSTANTS = (:P15_ECE_MARGIN, :P15_BREAK_COLUMNS, :P15_OOD_FIRE_MARGIN,
                                  :P15_GOLDEN_ECE_TOL)

    "Declared in advance so they are COMPARABLE across rungs, but they decide NOTHING."
    const P15_REPORTING_ONLY_CONSTANTS = (:P15_NUISANCE_COLUMNS, :P15_COVERAGE_LEVELS,
                                          :P15_KS_REPORTED, :P15_CHI2_REPORTED,
                                          :P15_AXIS_MECHANISM, :P15_PRIOR_BOUNDARY_RUNG)

    # =====================================================================================
    # The D-09 fast-gate golden bars (consumed by 15-06; frozen here with everything else)
    # =====================================================================================
    # The fast tier reuses the existing SMOKE FIXTURE scale, so it costs seconds and cannot become
    # a second reported run by accident.
    const P15_GOLDEN_M    = SBC_FIX_M   # 8
    const P15_GOLDEN_L    = SBC_FIX_L   # 15  (L + 1 = 16)
    const P15_GOLDEN_BINS = SBC_FIX_BINS # 4  (16 divisible by 4)

    # A CONCRETE size, not the `:mixture` sentinel: `gate_imsize` (`harness.jl:158-171`) returns a
    # concrete tuple outright and consumes ZERO numbers from the rng, so the golden's draw sequence
    # cannot shift if the mixture weights ever change. 256/8 = 32 ≥ 4, which clears
    # `_guard_misspec_imsize`'s patch-survivor floor.
    const P15_GOLDEN_IMSIZE = (256, 256)

    # DERIVATION: under a PINNED Julia (1.12.6) and the committed manifest, the only admissible
    # run-to-run drift in this statistic is float-summation reassociation, which is of order 1e-12
    # at these magnitudes. 1e-8 leaves four orders of headroom above that while still failing on any
    # real change to the net, the harness or the statistic. A JULIA VERSION BUMP IS A LEGITIMATE
    # RE-BLESS REASON and must be NAMED AS SUCH in the commit that re-blesses; a golden that can be
    # silently regenerated is a rubber stamp (D-09).
    const P15_GOLDEN_ECE_TOL = 1e-8     # golden-vs-recomputed tolerance for the fast CI tier

    "What the fast tier does and does NOT prove (D-09). Stated so SC3 is not oversold."
    const P15_GOLDEN_HONESTY =
        "At M = P15_GOLDEN_M = 8 the ECE null mean is ~0.336/sqrt(8) ~= 0.12 — an order of " *
        "magnitude above any miscalibration worth detecting. The fast gate therefore detects " *
        "CHANGE (a bit-level regression in the net, the harness or the statistic), NOT " *
        "miscalibration. Miscalibration is detected by the SLOW tier at reported scale. Any " *
        "claim that the per-push CI proves calibration oversells SC3."

    # =====================================================================================
    # Self-checks on the frozen constants (cheap, no inference, no simulation)
    # =====================================================================================
    @assert length(SBC_IMSIZE_SET) == length(SBC_IMSIZE_WEIGHTS)
    @assert isapprox(sum(SBC_IMSIZE_WEIGHTS), 1.0; atol = 1e-12)
    @assert all(>(0), SBC_IMSIZE_WEIGHTS)
    @assert (1376, 1028) in SBC_IMSIZE_SET            # the D-08 real-data anchor is present
    @assert !((256, 256) in SBC_IMSIZE_SET)           # the 07-05 budget artifact is excluded
    @assert (SBC_L + 1) % SBC_CHI2_BINS == 0          # exact rank-bin evenness (χ²)
    @assert (SBC_L + 1) % SBC_BINS      == 0          # exact rank-bin evenness (v1 ECE binning)
    @assert SBC_N_TESTS == 8
    @assert bf_required_n() == BF_GATE_N_MIN          # the derivation reproduces the frozen n_min
    @assert BF_GATE_N >= BF_GATE_N_MIN
    let ci = fisher_z_ci(BF_CORR_MIN, BF_GATE_N; tail = :two)
        @assert (BF_CORR_MIN - ci.lo) <= BF_CI_HALFWIDTH_TARGET   # precision target met at n
        @assert (ci.hi - BF_CORR_MIN) <= BF_CI_HALFWIDTH_TARGET
    end
    # The amended correlation rule is STRICTLY HARDER than v1's point rule at the gate's n.
    @assert bf_corr_verdict(BF_CORR_MIN, BF_GATE_N) != :pass
    # Seed disjointness (the runtests testset asserts the same set independently).
    @assert !(PROD_SEED_V2[8] in _forbidden_seeds())
    @assert length(unique(values(PROD_SEED_V2))) == 4
    @assert isempty(intersect(Set(values(PROD_SEED_V2)), Set(values(PROD_SEED))))

    # --- Phase-15 seed discipline (D-11) --------------------------------------------------
    # `test/gate/test_p15_consts.jl` asserts the same set independently AND shows the guard
    # firing on a deliberately-colliding derivation (CONVENTIONS C-04 rule 3).
    @assert !(PROD_SEED_P15[8] in _p15_forbidden())
    @assert length(unique(values(PROD_SEED_P15))) == 4
    @assert isempty(intersect(Set(values(PROD_SEED_P15)), Set(values(PROD_SEED))))
    @assert isempty(intersect(Set(values(PROD_SEED_P15)), Set(values(PROD_SEED_V2))))
    @assert PROD_SEED_P15[8] != PROD_SEED_V2[8]        # the SPENT seed, named by D-11
    @assert !(P15_SALT in P15_REPO_SALTS)              # a DISTINCT KEY, not a distinct index range
    @assert (UInt64(P15_MASTER) ⊻ P15_SALT) != (UInt64(AMEND_MASTER) ⊻ AMEND_SALT)
    @assert (UInt64(P15_MASTER) ⊻ P15_SALT) != (UInt64(PROD_MASTER) ⊻ PROD_SALT)
    # The aliases really do point at P15 — this is what carries `gate_global_seed` (harness.jl:82).
    @assert prod_seed(8) == PROD_SEED_P15[8]
    @assert rand(Philox4x(UInt64, (prod_seed(8), UInt64(1))), UInt64) !=
            rand(Philox4x(UInt64, (prod_seed_v2(8), UInt64(1))), UInt64)

    # --- Phase-15 Tier-1 bars: the derivations, re-executed at load time -------------------
    # The break margin IS the arithmetic in its own comment, not a number resembling it.
    @assert isapprox(P15_ECE_MARGIN, P15_ECE_NULL_Q95 - P15_ECE_NULL_MEAN; atol = 1e-12)
    @assert P15_ECE_NULL_Q95 > P15_ECE_NULL_MEAN > 0
    # The closed form 0.336/√M reproduces the measured null mean (the two are independent).
    @assert abs(P15_ECE_NULL_MEAN - 0.336 / sqrt(P15_SWEEP_M)) / P15_ECE_NULL_MEAN < 0.05
    # …and so does the independent B = 2000 research measurement, to within its own MC error.
    @assert abs(P15_ECE_NULL_MEAN - P15_ECE_NULL_MEAN_RESEARCH) < 0.0015
    @assert abs(P15_ECE_NULL_Q95  - P15_ECE_NULL_Q95_RESEARCH)  < 0.0030
    # The break threshold is the anchor PLUS the margin, never the margin alone.
    @assert p15_break_threshold(0.0) == P15_ECE_MARGIN

    # The OOD margin is DERIVED from the binomial null at the pre-registered n, not typed.
    @assert _p15_binom_quantile(200, 0.05, 0.95) == 15
    @assert isapprox(P15_OOD_FIRE_MARGIN,
                     _p15_binom_quantile(P15_OOD_N_POS, 1 - OOD_ID_QUANTILE, 0.95) /
                     P15_OOD_N_POS - (1 - OOD_ID_QUANTILE); atol = 1e-12)
    @assert P15_OOD_N_POS >= 200 && P15_OOD_N_ID >= 200 && P15_OOD_N_FIT >= 200

    # Ladders: strictly increasing, right length, rung 2 IS the prior boundary (D-03), and the
    # spillover endpoint respects the simulator's own guard (simulator.jl:232).
    for lad in (P15_SPILLOVER_LADDER, P15_SHIFT_LADDER, P15_AF_LADDER)
        @assert length(lad) == P15_RUNGS
        @assert issorted(collect(lad)) && allunique(lad)
        @assert all(>(0), lad)
    end
    @assert P15_SPILLOVER_LADDER[2] == P15_PRIOR_BOUNDARY.spillover
    @assert P15_SHIFT_LADDER[2]     == P15_PRIOR_BOUNDARY.registration
    @assert P15_AF_LADDER[2]        == P15_PRIOR_BOUNDARY.autofluorescence
    @assert P15_SPILLOVER_LADDER[1] <  P15_PRIOR_BOUNDARY.spillover        # rung 1 is INSIDE
    @assert P15_SHIFT_LADDER[1]     <  P15_PRIOR_BOUNDARY.registration
    @assert P15_AF_LADDER[1]        <  P15_PRIOR_BOUNDARY.autofluorescence
    @assert P15_SPILLOVER_LADDER[end] <= 1.0        # `0 ≤ spillover ≤ 1` (simulator.jl:232)

    # Axes and mechanism typing.
    @assert length(P15_AXES) == 7
    @assert allunique(P15_AXES)
    @assert keys(P15_AXIS_MECHANISM) == keys(P15_PRIOR_BOUNDARY_RUNG)
    @assert Tuple(keys(P15_AXIS_MECHANISM)) == P15_AXES
    @assert all(m -> m in (:in_prior, :out_of_model), values(P15_AXIS_MECHANISM))
    # Exactly the in-prior axes carry a boundary rung, and it is rung 2 on each (D-02a/D-03).
    for ax in P15_AXES
        if P15_AXIS_MECHANISM[ax] === :in_prior
            @assert P15_PRIOR_BOUNDARY_RUNG[ax] == 2
            @assert haskey(P15_PRIOR_BOUNDARY, ax)
        else
            @assert P15_PRIOR_BOUNDARY_RUNG[ax] === nothing
        end
    end
    @assert P15_AXIS_MECHANISM.background === :out_of_model     # D-02a's reclassification
    @assert length(P15_NATIVE_LEVELS) == P15_RUNGS

    # D-13 declarations exist and are non-empty BEFORE any rung does.
    @assert P15_ITERATION_ALLOWANCE == 1
    @assert length(P15_SC2_STATES) == 3 && :late in P15_SC2_STATES
    @assert all(!isempty, (P15_DEGENERATE_EMPTY_ENVELOPE_READING,
                           P15_DEGENERATE_UNBOUNDED_ENVELOPE_READING,
                           P15_SC2_PASS_RULE, P15_LATE_AXIS_RESPONSE, P15_GRID16_DROP_REASON,
                           P15_GOLDEN_HONESTY))
    @assert P15_GRID16_CONTRAST === :dropped
    @assert P15_MCE_FORBIDDEN

    # TIER SPLIT: gating and reporting-only names are disjoint, and every named constant exists.
    @assert isempty(intersect(P15_GATING_CONSTANTS, P15_REPORTING_ONLY_CONSTANTS))
    @assert allunique(P15_GATING_CONSTANTS) && allunique(P15_REPORTING_ONLY_CONSTANTS)
    # COLUMN SPLIT: every gating target is also reported, and no nuisance marginal gates (D-05).
    @assert isempty(intersect(P15_BREAK_COLUMNS, P15_NUISANCE_COLUMNS))
    @assert sort(vcat(collect(P15_BREAK_COLUMNS), collect(P15_NUISANCE_COLUMNS))) ==
            collect(P15_REPORTED_COLUMNS)
    @assert Tuple(sort(collect(intersect(P15_BREAK_COLUMNS, P15_REPORTED_COLUMNS)))) == (1, 8)
    @assert length(P15_COVERAGE_LEVELS) == 19

    # Golden bars ride the smoke fixture, and the fixture's rank binning stays exact.
    @assert (P15_GOLDEN_L + 1) % P15_GOLDEN_BINS == 0
    @assert P15_GOLDEN_IMSIZE isa Tuple{Int,Int}                # concrete: consumes no rng
    @assert P15_GOLDEN_IMSIZE[1] ÷ GATE_CONSTS_GRID >= 4        # clears _guard_misspec_imsize
    @assert 0 < P15_GOLDEN_ECE_TOL < 1e-6

    # Frozen-file pins are well-formed (their VALUES are checked against disk in the testset).
    @assert length(P15_FROZEN_FILE_SHA256) == 8
    @assert all(e -> length(e.sha256) == 64 && !isempty(e.path), P15_FROZEN_FILE_SHA256)
    @assert allunique(map(e -> e.path, values(P15_FROZEN_FILE_SHA256)))
    @assert P15_FROZEN_HASH_NORMALIZATION === :lf

    # Budget bars are self-consistent.
    @assert P15_SWEEP_M > 0 && P15_RUNGS > 0 && P15_RUNG0 == 0
    @assert P15_SWEEP_MIN_THREADS >= 16
    @assert P15_SWEEP_WALLCLOCK_CEILING_HOURS > 8.5     # the 7-axis prediction it must cover
end
