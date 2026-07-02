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

# spike/test/test_ood.jl --- Phase-5 OOD / misspecification-flag gate (OOD-01/02, D-03..06).
#
# The FAST per-task gate. It exercises the OOD machinery on a SMALL fixture (64×64, drawn
# from the FIXTURE stream val_rng(VAL_FIX_SEED) so the quick gate NEVER consumes the
# reserved reported VAL_MASTER_SEED stream, D-02). The FULL reported ROC grid at the
# LOCKED constants (OOD_GRID_LEVELS levels × 4 families at SBC_IMSIZE) is owned by Wave-3
# (05-04 run_ood.jl), which calls `ood_roc_over_grid`/`plot_ood_roc` on the frozen nulls.
#
# TRAIN-ONLY FREEZE (SC5/D-06): the Mahalanobis null is fit on the TRAIN pool ONLY; the
# ID-quantile OPERATING POINT is committed on a HELD-OUT ID pool (disjoint from the null
# fit) so it reflects true out-of-sample ID behaviour (an in-sample threshold is
# optimistically low — a fresh ID draw then exceeds it). Misspecified images are FULLY
# EXTERNAL and only ever SCORED.
#
# HONEST SEPARABILITY (reported in 05-03-SUMMARY.md, NOT tuned to pass): the fixed 8×8
# patch-Pearson summary is a CORRELATION statistic, so the density channel robustly
# detects families that ALTER correlation structure (optics/PSF — AUC≈1.0 at all levels)
# but is a MEASURED blind spot for detector-NOISE mismatch (AUC≈0 — the corruption
# decorrelates patches toward a central summary), a second named blind spot orthogonal to
# the D-04 correlation-preserving one. The fixture uses OPTICS as the robust positive; the
# full four-family characterization is the Wave-3 reported run's job.
#
# Mirrors test_sbc.jl / test_bf.jl: header → using Test → guarded include of ood.jl →
# fixture computed ONCE → one outer @testset. NO figure call (figures are script artifacts).

using Test
using Random
using Statistics

# Unit under test: the OOD surface (pulls ood.jl → harness.jl → consts.jl transitively).
isdefined(@__MODULE__, :fit_ood_nulls) || include(joinpath(@__DIR__, "..", "validation", "ood.jl"))

# Fixture-scale separability floor for the fast grid. The REPORTED gate (05-04) asserts
# the ROC AUC ≥ OOD_AUC_MIN (0.80) at the locked constants; here we assert the robust
# positive family clears a fixture floor.
const OOD_FIX_AUC_FLOOR = 0.9

# --- Frozen model + OOD fixture built ONCE (fixture stream; never the reported one) -----
const OOD_FIX_M = load_frozen_model()

function _build_ood_fixture(m)
    imsize = (64, 64); reps = 12; ppN = 150
    rng = val_rng(VAL_FIX_SEED)

    # (1) TRAIN pool → fit the Mahalanobis null (SC5). 150 > 64 continuous rows ⇒ the
    #     64×64 covariance is well-conditioned (not ridge-dominated).
    Ztrain = reduce(hcat, [id_summary(m, rng; imsize = imsize) for _ in 1:150])
    nulls  = fit_ood_nulls(Ztrain)

    # (2) HELD-OUT ID pool → the pre-registered operating point (out-of-sample, honest).
    mahaval  = [maha_score(nulls, id_summary(m, rng; imsize = imsize)) for _ in 1:60]
    maha_thr = id_threshold(mahaval)                 # q = OOD_ID_QUANTILE
    Zpp      = [id_summary(m, rng; imsize = imsize) for _ in 1:12]
    idpp     = [pp_mismatch_score(m, z; reps = reps, N = ppN, imsize = imsize,
                                  rng = val_rng(VAL_FIX_SEED)) for z in Zpp]
    pp_thr   = id_threshold(idpp)

    # (3) A comfortably-CENTRAL ID base pair (min maha of 30 fresh ID pairs) for the
    #     negative-control transforms — a representative ID image, robustly below threshold.
    idpairs  = [(θ = sample_prior(rng); simulate_pair(rng, θ; imsize = imsize)) for _ in 1:30]
    basepair = idpairs[argmin([maha_score(nulls, frozen_summary(m, p)) for p in idpairs])]
    Zbase    = frozen_summary(m, basepair)

    # (4) A strong POSITIVE control — optics/PSF mismatch (robustly separable, SC5-external).
    θo   = sample_prior(rng)
    Zopt = frozen_summary(m, misspec_optics(rng, θo; imsize = imsize, level = OOD_GRID_LEVELS))

    # (5) Fixture ROC over the grid (maha-only for speed; the reported run adds PP + PNG).
    roc = ood_roc_over_grid(m, nulls; levels = OOD_GRID_LEVELS, n_id = 20, n_pos = 6,
                            rng = val_rng(VAL_FIX_SEED), imsize = imsize, with_pp = false)

    return (imsize = imsize, reps = reps, ppN = ppN, Ztrain = Ztrain, nulls = nulls,
            maha_thr = maha_thr, pp_thr = pp_thr, id_fpr = mean(mahaval .> maha_thr),
            basepair = basepair, Zbase = Zbase, Zopt = Zopt, roc = roc)
end

const OOD_FIX = _build_ood_fixture(OOD_FIX_M)

# The OR-fused flag at the fixture's pre-registered thresholds. `reps`/`N`/`imsize` MUST
# match the values `pp_thr` was computed with (deterministic PP stream), so the PP channel
# is scored on the same footing as its threshold.
_ood_fix_flag(Z) = ood_flag(OOD_FIX.nulls, OOD_FIX_M, Z;
                            maha_thr = OOD_FIX.maha_thr, pp_thr = OOD_FIX.pp_thr,
                            rng = val_rng(VAL_FIX_SEED), reps = OOD_FIX.reps,
                            imsize = OOD_FIX.imsize, N = OOD_FIX.ppN)

@testset "Phase 5 — OOD / Misspecification Flag (OOD-01/02, D-03..06)" verbose = true begin

    @testset "SC4 (OOD-01) OR-flag fires on misspec, quiet in-distribution" begin
        # (a) the summary-density channel fires on the strong external positive control …
        @test maha_score(OOD_FIX.nulls, OOD_FIX.Zopt) > OOD_FIX.maha_thr
        # … and so the OR-combined flag (D-05) fires on it.
        @test _ood_fix_flag(OOD_FIX.Zopt) == true

        # (b) a central in-distribution image stays BELOW the operating point on the
        #     density channel, and the OR-flag stays QUIET on it.
        @test maha_score(OOD_FIX.nulls, OOD_FIX.Zbase) <= OOD_FIX.maha_thr
        @test _ood_fix_flag(OOD_FIX.Zbase) == false

        # (c) the operating point is an HONEST out-of-sample ID quantile: the held-out ID
        #     false-positive rate sits at ~1−OOD_ID_QUANTILE (not the optimistic in-sample 0).
        @test OOD_FIX.id_fpr <= 0.15
    end

    @testset "SC5 (D-04) summary-orthogonal negative control — the named blind spot" begin
        # Each correlation-preserving transform leaves the 8×8 patch-Pearson VALUE
        # distribution KS-invariant (< OOD_KS_EPS) — the MEASURED, NAMED blind spot — and
        # the flag stays QUIET on a central ID image (by construction).
        #
        # Nuance (honest, documented in 05-03-SUMMARY.md): affine leaves the summary VECTOR
        # bit-identical (Pearson is exactly scale/shift invariant), so the flag is provably
        # unchanged. Rotation (a grid isometry) and block-permute preserve only the value
        # MULTISET and permute the vector; because the density null encodes spatial neighbour
        # covariance, the position-aware Mahalanobis has RESIDUAL sensitivity to spatial
        # rearrangement (largest for a random block-permute). On a central ID base all three
        # stay quiet; a strongly OOD-leaning ID base plus block-permute could tip over — the
        # honest residual, not a hidden failure.
        negctrls = (("affine", (p -> negctrl_affine(p; a = 2.0, b = 0.5))),
                    ("rotate", (p -> negctrl_rotate_flip(p; k = 1))),
                    ("block",  (p -> negctrl_block_permute(val_rng(VAL_FIX_SEED), p; blocks = 8))))
        for (name, tf) in negctrls
            ks = verify_summary_invariance(OOD_FIX.basepair, tf)
            @test ks < OOD_KS_EPS                                      # summary KS-invariant (D-04)
            @test _ood_fix_flag(frozen_summary(OOD_FIX_M, tf(OOD_FIX.basepair))) == false  # flag quiet
        end

        # affine is the EXACT blind spot: the standardized summary is bit-identical on the
        # continuous rows, so the flag decision is provably invariant (strongest guarantee).
        Zaff = frozen_summary(OOD_FIX_M, negctrl_affine(OOD_FIX.basepair; a = 2.0, b = 0.5))
        @test Zaff[1:64] == OOD_FIX.Zbase[1:64]
    end

    @testset "SC5/D-06 nulls fit on TRAIN only + pre-registered operating point" begin
        # The Mahalanobis null was fit on the TRAIN continuous rows ONLY (T-05-07): its mean
        # equals the train continuous-row mean, and it never touched the mask rows (Pitfall 2)
        # nor any misspec/eval image.
        @test OOD_FIX.nulls.cont == collect(1:64)
        @test OOD_FIX.nulls.μS ≈ vec(mean(Float64.(OOD_FIX.Ztrain[1:64, :]); dims = 2))

        # the pass/fail thresholds the REPORTED run (05-04) scores against are the committed
        # pre-registered values — asserted-as-locked here, NEVER tuned to pass.
        @test OOD_ID_QUANTILE == 0.95     # D-06 operating-point quantile (~5% ID FPR)
        @test OOD_KS_EPS      == 0.05     # D-04 negative-control KS-invariance bound
        @test OOD_AUC_MIN     == 0.80     # D-03 positive-control separability floor (reported gate)
        @test OOD_GRID_LEVELS == 4        # D-03 misspec magnitudes per family
        # fixtures never consume the reserved reported stream (D-02).
        @test VAL_FIX_SEED != VAL_MASTER_SEED
    end

    @testset "OOD-02 controlled ROC over the misspec grid + post-hoc Youden reference" begin
        roc = OOD_FIX.roc
        # (a) the ROC AUC of the robust positive family clears the fixture floor at every
        #     magnitude; the REPORTED gate (05-04) scores this against OOD_AUC_MIN (0.80).
        @test length(roc.maha_auc[:optics]) == OOD_GRID_LEVELS
        @test all(a -> a >= OOD_FIX_AUC_FLOOR, roc.maha_auc[:optics])
        @test 0.0 <= roc.combined_auc <= 1.0

        # (b) the operating point reported is the PRE-REGISTERED ID quantile (D-06) …
        @test roc.id_quantile == OOD_ID_QUANTILE
        @test roc.auc_min == OOD_AUC_MIN
        # … and Youden-J is present as a POST-HOC REFERENCE ONLY — never used as a gate
        #    (no pass/fail assertion compares against roc.youden; D-06).
        @test isfinite(roc.youden.j)
        @test 0.0 <= roc.youden.tpr <= 1.0 && 0.0 <= roc.youden.fpr <= 1.0
    end

    @testset "OOD ran CPU-only (D-10)" begin
        @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
    end

end
