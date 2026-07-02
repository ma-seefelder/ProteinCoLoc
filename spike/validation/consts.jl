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

# spike/validation/consts.jl --- Phase-5 pre-registered constants (SBC/BF/OOD).
#
# THE ANTI-SNOOPING CONTRACT (SBC-03 / D-02). Every M/L/bin-count and every
# pass/fail threshold the Phase-5 validation bundle scores against is LOCKED in
# this ONE file and committed BEFORE any reported run. Plans 05-02 (BF) and 05-03
# (OOD) READ from here; the reported Wave-3 runs (05-04) consume VAL_MASTER_SEED.
# Changing any value below after a reported run would be "tune until calibrated"
# data-snooping (STATE Blockers, Phase 5) -- these are the committed pre-registration.
#
# SEED DISCIPLINE (D-02): VAL_MASTER_SEED is a FRESH Random123 stream, PROVABLY
# DISJOINT from the training/holdout NPE_MASTER_SEED = 0xC0FFEE (train_npe.jl:65).
# The reported numbers must never reuse the stream the net was trained on. A
# separate VAL_FIX_SEED drives the fast per-task fixtures so they never consume
# (and thus never leak / pre-observe) the reserved reported stream.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local constants; touches no src/.
#
# Guarded as ONE block keyed on :SBC_M so re-inclusion under runtests.jl (each of
# harness.jl / sbc.jl / test_sbc.jl guards this include) is a silent no-op --
# a redefinition to the same value would otherwise warn on a `const`.

if !isdefined(@__MODULE__, :SBC_M)
    # --- SBC (SBC-01/02/03) : reported run --------------------------------------
    const SBC_M          = 2000     # pre-registered SBC draws (θ*~π→simulate→infer)
    const SBC_L          = 999      # posterior draws per SBC draw; L+1 = 1000 …
    const SBC_BINS       = 50       # … divisible by SBC_BINS (rank-bin evenness, Pitfall 4)
    const SBC_KS_ALPHA   = 0.05     # KS rank-uniformity pass threshold
    const SBC_CHI2_ALPHA = 0.05     # χ² rank-uniformity pass threshold
    const SBC_ECE_GREEN  = 0.05     # ECE ≤ GREEN  → :green traffic-light verdict
    const SBC_ECE_YELLOW = 0.10     # ECE ≤ YELLOW → :yellow, else :red
    const SBC_IMSIZE     = (256, 256)  # reported simulate_pair image size (dims ≥ 64)

    # --- SBC : fast fixture (per-task gate; NEVER the reported numbers) ----------
    const SBC_FIX_M      = 50       # tiny SBC-draw count for the quick gate
    const SBC_FIX_L      = 99       # tiny posterior-draw count (fixture)
    const SBC_FIX_BINS   = 10       # tiny rank-bin count (fixture; L+1 = 100 divisible)

    # --- Amortized Bayes factor (BF-01/02, D-08) --------------------------------
    const BF_CORR_MIN    = 0.95     # D-08a: cor(amortized logBF, KDE-baseline BF) ≥
    const BF_LOGBF_TOL   = 0.5      # D-08b: max |Δ logBF| over the Δρ sweep ≤
    const BF_SWEEP_LO    = -0.6     # D-08 held-out Δρ sweep lower bound
    const BF_SWEEP_HI    = 0.8      # D-08 held-out Δρ sweep upper bound
    const BF_SWEEP_N     = 25       # D-08 Δρ sweep grid points

    # --- OOD / misspecification flag (OOD-01..04, D-03..06) ---------------------
    const OOD_ID_QUANTILE = 0.95    # ID score quantile = operating point (~5% ID FPR, D-06)
    const OOD_KS_EPS      = 0.05    # D-04 negative-control KS-invariance statistic bound
    const OOD_AUC_MIN     = 0.80    # positive-control separability floor (ROC AUC)
    const OOD_GRID_LEVELS = 4       # D-03 misspec magnitudes per family
    const OOD_PP_REPS     = 50      # posterior-predictive mismatch replicate count

    # --- Pre-registered seeds (D-02) --------------------------------------------
    # RESERVED reported stream, consumed ONLY by the Wave-3 run (05-04). DISJOINT
    # from NPE_MASTER_SEED = 0xC0FFEE (train_npe.jl:65) so the reported SBC/BF/OOD
    # numbers never reuse the stream the net was trained/held-out on (D-02).
    const VAL_MASTER_SEED = 0x5BC0FFEE   # ≠ 0xC0FFEE (NPE_MASTER_SEED): fresh disjoint stream
    # Fixture-only stream for the fast per-task gates. Distinct from VAL_MASTER_SEED
    # so the quick gates never consume (never pre-observe) the reported stream.
    const VAL_FIX_SEED    = 0xF1F7ED
end
