#########################################################################################
# scripts/train_grid8_amended_v2.jl --- the AMENDED-REGIME grid-8 retrain (07-GATE-AMENDMENT §6.1).
#
# Protocol §6 step 1 ("Retrain first"): grid 8 is retrained with
#   • the BOUNDED (logit-of-prior-box) θ representation now implemented in
#     `src/amortized/architecture.jl` (F2 remedy, commit 31a4733), and
#   • the F5 REALISTIC IMAGE-SIZE MIXTURE — `SBC_IMSIZE_SET` / `SBC_IMSIZE_WEIGHTS` READ DIRECTLY
#     FROM the frozen `test/gate/gate_consts_8_v2.jl` (never retyped here: the amendment's binding
#     invariant is that the TRAINING joint equals the GATE joint, and the only way to guarantee
#     that mechanically is for both to read the same frozen constants),
# everything else identical to run 07-05 / the bounded-θ diagnostic: `n_pairs = 50_000`,
# `use_gpu = false`, default architecture + optimizer recipe, `master_seed = DEFAULT_MASTER_SEED`.
#
# The bundle is written to a NEW artifact root, `artifacts/amended_v2/grid_8/`. The frozen stores
# `artifacts/grid_{4,8,16}/` and `artifacts/bounded_theta/grid_8/` are NEVER touched.
#
# TRAINING ONLY. This script does NOT run any gate arm. The amended ship-gate is bound by
# protocol §6.3 to exactly ONE run on the fresh `PROD_SEED_V2[8]`, which is a separate, deliberate
# step — invoking it here would consume the single confirmatory opportunity.
#
# Run:  julia --project -t auto scripts/train_grid8_amended_v2.jl
#########################################################################################
using ProteinCoLoc

# The frozen amended pre-registration, loaded into an ISOLATED module: it defines the same const
# names as `gate_consts_8.jl`, and we want ONLY its imsize mixture here. Nothing is run against it.
module GateV2
    include(joinpath(@__DIR__, "..", "test", "gate", "gate_consts_8_v2.jl"))
end

const ROOT = normpath(joinpath(@__DIR__, "..", "artifacts", "amended_v2"))

# The F5 mixture, taken verbatim from the frozen file (NOT retyped).
const AMENDED_IMSIZE_SET     = GateV2.SBC_IMSIZE_SET
const AMENDED_IMSIZE_WEIGHTS = GateV2.SBC_IMSIZE_WEIGHTS

# Guard: if the frozen file is ever edited, this run must not silently train a different joint.
@assert length(AMENDED_IMSIZE_SET) == length(AMENDED_IMSIZE_WEIGHTS)
@assert isapprox(sum(AMENDED_IMSIZE_WEIGHTS), 1.0; atol = 1e-12)
@assert (1376, 1028) in AMENDED_IMSIZE_SET          # the D-08 real-data anchor
@assert !((256, 256) in AMENDED_IMSIZE_SET)         # the 07-05 compute-budget artifact, excluded

@info "amended-regime grid-8 retrain (bounded θ + F5 imsize mixture)" ROOT AMENDED_IMSIZE_SET AMENDED_IMSIZE_WEIGHTS nthreads = Threads.nthreads()
t0 = time()
b = ProteinCoLoc._train_grid_pipeline(8;
        imsize_set     = AMENDED_IMSIZE_SET,        # F5 mixture, read from the frozen v2 consts
        imsize_weights = AMENDED_IMSIZE_WEIGHTS,
        n_pairs        = 50_000,            # identical to 07-05
        use_gpu        = false,             # identical to 07-05
        artifacts_root = ROOT,              # NEW store — every prior artifact untouched
        skip_if_done   = true,
        verbose        = true)
elapsed = (time() - t0) / 60
@info "training done" minutes = elapsed grid = b.grid

# --- provenance READBACK (the amendment's SBC_REQUIRE_IMSIZE_PROVENANCE invariant) -------------
npe_path = joinpath(ROOT, "grid_8", "npe_8.jld2")
prov = ProteinCoLoc.training_imsize_provenance(ProteinCoLoc.load_estimator(npe_path))
@info "training_imsize_provenance readback" prov
@assert prov.recorded
@assert Tuple(prov.imsize_set) == AMENDED_IMSIZE_SET
@assert Tuple(prov.imsize_weights) == AMENDED_IMSIZE_WEIGHTS
@assert prov.imsize_source === :generate_samples
@info "provenance matches the frozen SBC_IMSIZE_SET / SBC_IMSIZE_WEIGHTS exactly"
