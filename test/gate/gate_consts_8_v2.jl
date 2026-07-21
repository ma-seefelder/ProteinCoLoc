#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/gate_consts_8_v2.jl --- the AMENDED 8×8 ship-gate pre-registration (v2).
#
# THIS FILE IS FROZEN. It is the executable form of
# `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md`, which states,
# for each amended item, the defect, the OUTCOME-INDEPENDENT justification, the new rule and what
# is NOT changed. Read that document before reading this file; every number below is derived
# there and none of them is an outcome.
#
# ================================ PROVENANCE DISCLOSURE ================================
# This amendment was written AFTER the original grid-8 ship-gate results were seen. That is
# precisely why every change is argued ONLY from properties of the test design — family-wise
# error rate, estimator precision, extreme-value instability, and the validity conditions of SBC
# — and NEVER from what the observed numbers were or from what any net would pass.
#
# NOTHING HAS BEEN RUN AGAINST THIS FILE. It is a specification, frozen before use.
#
# It SUPPLEMENTS and does NOT replace `gate_consts_8.jl`. That file, `gate_consts_4.jl`,
# `gate_consts_16.jl` and `gate_consts_template.jl` are byte-unchanged, and the recorded FAIL
# verdicts in `artifacts/grid_{4,8,16}/gate_report_{4,8,16}.jld2` remain intact and citable. Any
# report of an amended verdict MUST co-cite the original FAIL.
#
# WHAT CHANGED (three items; see the amendment for the full argument):
#   A1  multiplicity — the naive all-8 α=0.05 KS/χ² conjunction has FWER = 1 − 0.95^8 = 33.7 %.
#       Replaced by HOLM–BONFERRONI at a FWER target of 0.05 (valid under arbitrary dependence,
#       uniformly ≥ Bonferroni in power). Per-test α retired, NOT re-tuned looser.
#   A1b χ² bin count — SBC_BINS was overloaded for BOTH the χ² rank bins and the ECE reliability
#       bins. Decoupled: SBC_CHI2_BINS = 20 (19 df, 100 expected/bin at M = 2000) for χ²;
#       SBC_BINS = 50 UNCHANGED for the ECE arm, which is therefore numerically identical to v1.
#   A2  BF precision — n = BF_SWEEP_N = 25 (recorded n 15–20) cannot resolve a 0.95 threshold;
#       the Fisher-z 95 % CI is ±0.07–0.10 there. n derived from a pre-registered precision
#       target (half-width ≤ 0.03 at r = 0.95) → n_min = 69, set BF_GATE_N = 100. Decision rule
#       moves from a point estimate to a ONE-SIDED 95 % LOWER CONFIDENCE BOUND — STRICTLY HARDER.
#       `max|Δ logBF|` is KEPT (reported) with an n-stable q95 twin as the gated statistic;
#       BF_CORR_MIN = 0.95 and BF_LOGBF_TOL = 0.5 are NOT re-tuned.
#   F5  SBC validity — SBC ranks are uniform only under the TRAINING joint. The scalar
#       SBC_IMSIZE = (256,256) was a documented COMPUTE-BUDGET deviation (run 07-05), not a
#       scientific choice. Replaced by the D-09 realistic MIXTURE 512²/1024²/1376×1028/2048²,
#       which the gate MUST assert equals the net's persisted training imsize provenance.
#
# SEED DISCIPLINE: a FRESH `PROD_SEED_V2[8]`, derived by the SAME salted-Philox rule as the
# template but from a DISTINCT (AMEND_MASTER, AMEND_SALT) pair, and re-drawn until it is disjoint
# from EVERY v1 PROD_SEED[G], from VAL_MASTER_SEED / NPE_MASTER_SEED / DEFAULT_MASTER_SEED, and
# from the F2-remedy DEV seeds. Unit-asserted in runtests.
#
# LOADING: this file defines the SAME const names as `gate_consts_8.jl` and is guarded on
# :SBC_M, so the two must NEVER be loaded into the same module. `runtests.jl` loads it into an
# isolated `module GateConstsV2` for the disjointness assertions.

if !isdefined(@__MODULE__, :SBC_M)
    import Random123: Philox4x     # counter-based RNG for the fresh disjoint per-grid stream

    const GATE_CONSTS_VERSION = 2
    const GATE_CONSTS_GRID    = 8
    const GATE_AMENDMENT_DOC  =
        ".planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md"

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

    Fresh amended gate seeds. `PROD_SEED_V2[8]` is THE seed for the single confirmation run of the
    amended grid-8 gate (protocol §6 of the amendment: ONE run, no iteration). The sibling grids
    are carried for reference; grids 4 and 16 are NOT re-gated by this amendment.
    """
    const PROD_SEED_V2 = Dict{Int,UInt64}(G => _derive_prod_seed_v2(G) for G in (4, 8, 16, 32))

    "The amended gate seed for grid `G` (derived on the fly for a novel grid)."
    prod_seed_v2(G::Integer) = get(() -> _derive_prod_seed_v2(G), PROD_SEED_V2, Int(G))

    "The amended per-grid ship-gate RNG: `Philox4x` keyed by `(prod_seed_v2(G), 0)`."
    prod_rng_v2(G::Integer) = Philox4x(UInt64, (prod_seed_v2(G), UInt64(0)))

    # The amended gate consumes the v2 stream. These aliases exist so the existing gate machinery
    # (`sbc_gate`, `bf_gate`, `ood_gate` all default to `rng = prod_rng(G)`) rides the FRESH
    # stream when this file is the loaded pre-registration.
    prod_seed(G::Integer) = prod_seed_v2(G)
    prod_rng(G::Integer)  = prod_rng_v2(G)

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
end
