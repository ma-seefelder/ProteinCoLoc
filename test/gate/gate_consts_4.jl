#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/gate_consts_4.jl --- the FRESH 4×4 ship-gate pre-registration (D-05, PROD-02).
#
# This is the committed, locked pre-registration for the 4×4 coarse-but-robust grid — instantiated
# from gate_consts_template.jl and committed BEFORE the reported 4×4 gate runs (anti-snooping,
# Pitfall 3 / T-7-08). Every M/L/bin count and every pass/fail threshold below is FROZEN here;
# changing any value after the reported run would be "tune until calibrated" data-snooping.
#
# 4×4 is the coarse-but-robust grid (summary_dim(4)=32, cont_rows(4)=16). Each 64×64 patch at
# SBC_IMSIZE=(256,256) carries ≥4096 px — far above the ≥15-survivor floor — so the 4×4 grid is
# robust at any image size and the reported 256² gate needs no per-grid image-size increase (unlike
# the fine grids). The reported SBC/BF/OOD thresholds are the same Phase-5-class values the memo §6
# pass conditions are stated against; they are re-expressed FRESH here for the 4×4 grid.
#
# SEED DISCIPLINE (the headline D-05 delta): the reported 4×4 gate rides the FRESH disjoint stream
# prod_rng(4) = Philox4x(PROD_SEED[4], 0), keyed on (PROD_MASTER ⊻ PROD_SALT, 4). PROD_SEED[4] is a
# full-width Philox draw, PROVABLY DISJOINT from — and never reusing —
#   • the spike TRAINING stream NPE_MASTER_SEED = 0xC0FFEE,
#   • the spike VALIDATION stream VAL_MASTER_SEED = 0x5BC0FFEE, and
#   • the productionization DATAGEN stream DEFAULT_MASTER_SEED = 0x1 (this grid's training pool).
# It is a distinct draw from PROD_SEED[8] (the grid index G=4 keys a different Philox counter) and
# is asserted distinct from the two forbidden spike seeds (runtests disjointness test).
#
# Guarded as ONE block keyed on :SBC_M so re-inclusion (run_gate.jl / harness.jl / sbc.jl each guard
# their include) is a silent no-op.

if !isdefined(@__MODULE__, :SBC_M)
    import Random123: Philox4x     # counter-based RNG for the fresh disjoint per-grid stream

    # --- SBC (reported 4×4 run) ------------------------------------------------------------
    const SBC_M          = 2000        # pre-registered SBC draws (θ*~π → simulate → CPU-infer)
    const SBC_L          = 999         # posterior draws per SBC draw; L+1 = 1000 …
    const SBC_BINS       = 50          # … divisible by SBC_BINS (rank-bin evenness)
    const SBC_KS_ALPHA   = 0.05        # KS rank-uniformity pass threshold
    const SBC_CHI2_ALPHA = 0.05        # χ² rank-uniformity pass threshold
    const SBC_ECE_GREEN  = 0.05        # ECE ≤ GREEN  → :green traffic-light verdict
    const SBC_ECE_YELLOW = 0.10        # ECE ≤ YELLOW → :yellow, else :red
    const SBC_IMSIZE     = (256, 256)  # reported simulate_pair image size (robust for 4×4:
                                       # 64×64 px/patch ≥4096 px, far above the ≥15-survivor floor)

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
    const NPE_MASTER_SEED = 0x0000_0000_00C0_FFEE   # spike training stream   (forbidden)
    const VAL_MASTER_SEED = 0x0000_0000_5BC0_FFEE   # spike validation stream (forbidden)

    # --- Fresh disjoint productionization gate stream (D-05) --------------------------------
    # PROD_SALT / PROD_MASTER are the SAME committed mixing constants as the template, so
    # PROD_SEED[4] is the reproducible, pre-registered fresh 4×4 gate seed. The grid index G=4 keys
    # a different Philox counter than G=8, so PROD_SEED[4] is a distinct draw from PROD_SEED[8].
    const PROD_SALT   = 0x94D0_49BB_1331_11EB
    const PROD_MASTER = 0x0000_0000_09E3_779B

    """
        _derive_prod_seed(G) -> UInt64

    The fresh, disjoint per-grid gate seed: a full-width Philox draw keyed on
    `(PROD_MASTER ⊻ PROD_SALT, G)`, re-drawn on the astronomically unlikely event it collides with
    a forbidden seed (or 0). Distinct from `VAL_MASTER_SEED` / `NPE_MASTER_SEED` (unit-asserted).
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

    Fresh disjoint gate seeds for the shipped grid family (D-04). `PROD_SEED[4]` is the reported
    4×4 gate seed used to build `prod_rng(4)`; the sibling grids are carried for reference.
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

    The per-grid ship-gate RNG: a Random123 `Philox4x` keyed by `(prod_seed(G), 0)`, PROVABLY
    DISJOINT from the spike training/validation streams. Reproducible for a fixed `G`.
    """
    prod_rng(G::Integer) = Philox4x(UInt64, (prod_seed(G), UInt64(0)))
end
