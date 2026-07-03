#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/gate_consts_template.jl --- FRESH-per-grid pre-registration template (D-05).
#
# THE ANTI-SNOOPING CONTRACT (Pitfall 3 / T-7-08). Every M/L/bin-count and every pass/fail
# threshold a per-grid ship-gate scores against is LOCKED here and committed BEFORE any reported
# run. Each shipped grid gets its OWN `gate_consts_<G>.jl` — a COPY of this template re-expressed
# FRESH (its own M/L/imsize + its own disjoint `PROD_SEED[G]`) BEFORE that grid's gate runs.
# Changing any value below after a reported run would be "tune until calibrated" data-snooping.
#
# SEED DISCIPLINE (D-05, the headline delta from the spike): the per-grid gate stream is a FRESH
# Random123 stream `prod_rng(G)`, PROVABLY DISJOINT from — and NEVER reusing —
#   • `NPE_MASTER_SEED = 0xC0FFEE`   (the spike TRAINING stream), and
#   • `VAL_MASTER_SEED = 0x5BC0FFEE` (the spike VALIDATION stream).
# The net was iterated against BOTH, so reusing either for a productionization gate is snooping.
# `PROD_SEED[G]` is a full-width Philox draw keyed on `(PROD_MASTER ⊻ PROD_SALT, G)`, asserted
# distinct from the two forbidden seeds (see the runtests disjointness test).
#
# BF / OOD ride-alongs (07-RESEARCH §Validation Architecture, Memo §5): the BF gate scores the
# amortized log-BF against the NON-CLAMPED KDE baseline (`kde_log_bf_unclamped`, bf.jl) so
# `max|Δ logBF|` is free of the 1e-8 clamp-tail artifact; the OOD gate runs with `with_pp = true`
# (the re-enabled posterior-predictive channel).
#
# Guarded as ONE block keyed on :SBC_M so re-inclusion (harness.jl / sbc.jl / run_gate.jl each
# guard this include) is a silent no-op — a redefinition to the same value would otherwise warn on
# a `const`.

if !isdefined(@__MODULE__, :SBC_M)
    import Random123: Philox4x     # counter-based RNG for the fresh disjoint per-grid stream

    # --- SBC (reported run) : re-express FRESH per grid -------------------------------------
    const SBC_M          = 2000        # pre-registered SBC draws (θ*~π → simulate → CPU-infer)
    const SBC_L          = 999         # posterior draws per SBC draw; L+1 = 1000 …
    const SBC_BINS       = 50          # … divisible by SBC_BINS (rank-bin evenness)
    const SBC_KS_ALPHA   = 0.05        # KS rank-uniformity pass threshold
    const SBC_CHI2_ALPHA = 0.05        # χ² rank-uniformity pass threshold
    const SBC_ECE_GREEN  = 0.05        # ECE ≤ GREEN  → :green traffic-light verdict
    const SBC_ECE_YELLOW = 0.10        # ECE ≤ YELLOW → :yellow, else :red
    const SBC_IMSIZE     = (256, 256)  # reported simulate_pair image size (raise per grid so
                                       # fine grids clear the ≥15-survivor patch floor)

    # --- SBC : fast fixture (per-task smoke; NEVER the reported numbers) --------------------
    const SBC_FIX_M      = 8           # tiny SBC-draw count for the quick harness smoke
    const SBC_FIX_L      = 15          # tiny posterior-draw count (fixture; L+1 = 16)
    const SBC_FIX_BINS   = 4           # tiny rank-bin count (fixture; 16 divisible by 4)

    # --- Amortized Bayes factor (BF gate; non-clamped baseline) -----------------------------
    const BF_CORR_MIN    = 0.95        # cor(amortized logBF, non-clamped KDE baseline) ≥
    const BF_LOGBF_TOL   = 0.5         # max |Δ logBF| over the Δρ sweep ≤
    const BF_SWEEP_LO    = -0.6        # held-out Δρ sweep lower bound
    const BF_SWEEP_HI    = 0.8         # held-out Δρ sweep upper bound
    const BF_SWEEP_N     = 25          # Δρ sweep grid points

    # --- OOD / misspecification flag (OOD gate; with_pp = true) -----------------------------
    const OOD_ID_QUANTILE = 0.95       # ID score quantile = operating point (~5% ID FPR)
    const OOD_KS_EPS      = 0.05       # negative-control KS-invariance statistic bound
    const OOD_AUC_MIN     = 0.80       # positive-control separability floor (ROC AUC)
    const OOD_GRID_LEVELS = 4          # misspec magnitudes per family
    const OOD_PP_REPS     = 50         # posterior-predictive mismatch replicate count

    # --- FORBIDDEN seeds (anti-snooping, Pitfall 3) ----------------------------------------
    # The spike training + validation streams. A productionization gate must NEVER reuse either.
    const NPE_MASTER_SEED = 0x0000_0000_00C0_FFEE   # spike training stream   (forbidden)
    const VAL_MASTER_SEED = 0x0000_0000_5BC0_FFEE   # spike validation stream (forbidden)

    # --- Fresh disjoint productionization gate stream (D-05) --------------------------------
    # PROD_SALT is a SplitMix64 mixing constant, distinct from the spike VAL_SALT
    # (0xBF58476D1CE4E5B9), FOLD_SALT and HOLDOUT_SALT. PROD_SEED[G] is a full-width Philox draw
    # keyed on (PROD_MASTER ⊻ PROD_SALT, G), guaranteed distinct from the two forbidden seeds so
    # the reported per-grid numbers never reuse a stream the net was trained/validated against.
    const PROD_SALT   = 0x94D0_49BB_1331_11EB
    const PROD_MASTER = 0x0000_0000_09E3_779B

    """
        _derive_prod_seed(G) -> UInt64

    The fresh, disjoint per-grid gate seed: a full-width Philox draw keyed on
    `(PROD_MASTER ⊻ PROD_SALT, G)`, re-drawn on the astronomically unlikely event it collides with
    a forbidden seed (or 0). Distinct for each grid, and distinct from `VAL_MASTER_SEED` /
    `NPE_MASTER_SEED` (unit-asserted).
    """
    function _derive_prod_seed(G::Integer)
        rng = Philox4x(UInt64, (UInt64(PROD_MASTER) ⊻ PROD_SALT, UInt64(G)))
        s = rand(rng, UInt64)
        while s in (VAL_MASTER_SEED, NPE_MASTER_SEED, UInt64(0))
            s = rand(rng, UInt64)
        end
        return s
    end

    """
        PROD_SEED :: Dict{Int,UInt64}

    Fresh disjoint gate seeds for the shipped grid family (D-04): `PROD_SEED[G]` for
    `G ∈ (4, 8, 16, 32)`. Each is bracket-indexed by grid and used to build the reported gate RNG
    `prod_rng(G)`. A novel grid derives its seed on the fly via `prod_seed`.
    """
    const PROD_SEED = Dict{Int,UInt64}(G => _derive_prod_seed(G) for G in (4, 8, 16, 32))

    """
        prod_seed(G) -> UInt64

    The fresh disjoint gate seed for grid `G`: the pre-registered `PROD_SEED[G]` for a shipped
    grid, else derived on the fly (same rule) for a novel grid.
    """
    prod_seed(G::Integer) = get(() -> _derive_prod_seed(G), PROD_SEED, Int(G))

    """
        prod_rng(G) -> Philox4x

    The per-grid ship-gate RNG: a Random123 `Philox4x` keyed by `(prod_seed(G), 0)`. Its draw
    stream is PROVABLY DISJOINT from the spike training/validation streams (fresh full-width seed);
    two constructions with the same `G` are bit-identical, so the gate is reproducible.
    """
    prod_rng(G::Integer) = Philox4x(UInt64, (prod_seed(G), UInt64(0)))
end
