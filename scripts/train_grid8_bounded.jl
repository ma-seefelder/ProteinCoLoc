#########################################################################################
# scripts/train_grid8_bounded.jl --- the F2 bounded-θ CONTROLLED EXPERIMENT (diagnostic).
#
# Retrains grid 8 with the SAME configuration as the frozen run 07-05
# (`imsize_set = ((256,256),)`, `n_pairs = 50_000`, default architecture/optimizer recipe,
# `master_seed = DEFAULT_MASTER_SEED`, `use_gpu = false`) so the ONLY intended difference is the
# bounded (logit-of-prior-box) θ representation introduced for 07-CALIBRATION-FINDINGS F2.
#
# The bundle is written to a SEPARATE artifact root; `artifacts/grid_8/` (the frozen, gated
# bundle) is NEVER touched. This is a DEVELOPMENT run: it is evaluated on a DEV seed only and
# MUST NOT be run through the pre-registered ship-gate.
#
# Run:  julia --project -t auto scripts/train_grid8_bounded.jl
#########################################################################################
using ProteinCoLoc

const ROOT = normpath(joinpath(@__DIR__, "..", "artifacts", "bounded_theta"))

@info "bounded-θ grid-8 retrain (single-variable experiment vs run 07-05)" ROOT nthreads=Threads.nthreads()
t0 = time()
b = ProteinCoLoc._train_grid_pipeline(8;
        imsize_set     = ((256, 256),),     # identical to 07-05 (matches SBC_IMSIZE)
        n_pairs        = 50_000,            # identical to 07-05
        use_gpu        = false,             # identical to 07-05 (persisted meta says use_gpu=false)
        artifacts_root = ROOT,              # SEPARATE store — artifacts/grid_8/ untouched
        skip_if_done   = true,
        verbose        = true)
@info "done" minutes=(time() - t0) / 60 grid=b.grid
