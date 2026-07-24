#########################################################################################
# Spike 012 --- the CAPACITY LEVER: grid-8 retrain with a HIGHER-CAPACITY FLOW.
#
# QUESTION. After the bounded-θ fix, `autofluorescence` and `label_efficiency` still show a
# small F2 marginal LOCATION drift (~0.05–0.07 posterior-SD) on the amended_v2 net. These are
# NON-identified nuisances (shrinkage ≈ 1.0). Does raising the FLOW's capacity reduce that
# marginal drift below the 0.10-SD equivalence margin (07-NUISANCE-SBC-SPEC §3) — while keeping
# the coloc targets ρ_true / Δρ clean?
#
# THE LEVER (flow only; the summary net governs identifiability, not the marginal — left alone):
#     num_coupling_layers  10 → 16   (NPE_COUPLING)
#     flow_width          128 → 256   (NPE_FLOW_WIDTH)
#     flow_depth            2 →  3    (NPE_FLOW_DEPTH)
# passed as `npe_kwargs` into `_train_grid_pipeline` → `train_npe` → `build_estimator`. src/ is
# NOT edited (the knobs are already keyword-plumbed through the whole pipeline).
#
# EVERYTHING ELSE IDENTICAL TO amended_v2 (scripts/train_grid8_amended_v2.jl):
#   • the F5 image-size MIXTURE, READ DIRECTLY from the frozen test/gate/gate_consts_8_v2.jl
#     (never retyped — the training joint must equal the gate joint by construction),
#   • n_pairs = 50_000, use_gpu = false, master_seed = DEFAULT_MASTER_SEED (⇒ the 50k datagen
#     pool is BYTE-IDENTICAL to amended_v2; only the flow capacity differs),
#   • default optimizer recipe / epochs / early stopping.
#
# NEW artifact root: artifacts/spike012_highcap/grid_8/. amended_v2 and every frozen bundle are
# NEVER touched (sha256 verified before/after in the spike README).
#
# TRAINING ONLY. No gate arm is run; no PROD seed is consumed.
#
# Run:  julia --project -t auto .planning/spikes/012-capacity-lever/train_highcap.jl
#########################################################################################
using ProteinCoLoc

# The frozen amended pre-registration, loaded into an ISOLATED module: we want ONLY its imsize
# mixture. Nothing is run against it.
module GateV2
    include(joinpath(@__DIR__, "..", "..", "..", "test", "gate", "gate_consts_8_v2.jl"))
end

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "artifacts", "spike012_highcap"))

# The F5 mixture, taken verbatim from the frozen file (NOT retyped).
const IMSIZE_SET     = GateV2.SBC_IMSIZE_SET
const IMSIZE_WEIGHTS = GateV2.SBC_IMSIZE_WEIGHTS

# Guard: if the frozen file is ever edited, this run must not silently train a different joint.
@assert length(IMSIZE_SET) == length(IMSIZE_WEIGHTS)
@assert isapprox(sum(IMSIZE_WEIGHTS), 1.0; atol = 1e-12)
@assert (1376, 1028) in IMSIZE_SET          # the D-08 real-data anchor
@assert !((256, 256) in IMSIZE_SET)         # the 07-05 compute-budget artifact, excluded

# The capacity lever (flow only). These are the amended_v2 defaults + the increments.
const HIGHCAP_NPE_KWARGS = (; num_coupling_layers = 16, flow_depth = 3, flow_width = 256)
@assert HIGHCAP_NPE_KWARGS.num_coupling_layers > ProteinCoLoc.NPE_COUPLING     # 16 > 10
@assert HIGHCAP_NPE_KWARGS.flow_width          > ProteinCoLoc.NPE_FLOW_WIDTH   # 256 > 128
@assert HIGHCAP_NPE_KWARGS.flow_depth          > ProteinCoLoc.NPE_FLOW_DEPTH   # 3 > 2

@info "spike-012 high-capacity grid-8 retrain (flow lever, everything else = amended_v2)" ROOT IMSIZE_SET IMSIZE_WEIGHTS HIGHCAP_NPE_KWARGS baseline_flow = (ProteinCoLoc.NPE_COUPLING, ProteinCoLoc.NPE_FLOW_DEPTH, ProteinCoLoc.NPE_FLOW_WIDTH) nthreads = Threads.nthreads()

t0 = time()
b = ProteinCoLoc._train_grid_pipeline(8;
        imsize_set     = IMSIZE_SET,            # F5 mixture, read from the frozen v2 consts
        imsize_weights = IMSIZE_WEIGHTS,
        n_pairs        = 50_000,                # identical to amended_v2
        use_gpu        = false,                 # identical to amended_v2
        artifacts_root = ROOT,                  # NEW store — every prior artifact untouched
        npe_kwargs     = HIGHCAP_NPE_KWARGS,    # THE CAPACITY LEVER (flow only)
        skip_if_done   = true,
        verbose        = true)                  # per-epoch train/val risk → the training curve
elapsed = (time() - t0) / 60
@info "training done" minutes = elapsed grid = b.grid

# --- provenance READBACK (mirrors the amended_v2 script's invariant check) ----------------------
npe_path = joinpath(ROOT, "grid_8", "npe_8.jld2")
loaded   = ProteinCoLoc.load_estimator(npe_path)
prov     = ProteinCoLoc.training_imsize_provenance(loaded)
@info "training_imsize_provenance readback" prov
@assert prov.recorded
@assert Tuple(prov.imsize_set) == IMSIZE_SET
@assert Tuple(prov.imsize_weights) == IMSIZE_WEIGHTS
@assert prov.imsize_source === :generate_samples

# --- confirm the persisted architecture actually carries the raised flow capacity ---------------
@info "persisted architecture" arch = loaded.arch
@assert loaded.arch.num_coupling_layers == 16
@assert loaded.arch.flow_depth == 3
@assert loaded.arch.flow_width == 256
# the summary net (identifiability knobs) must be UNCHANGED from amended_v2 defaults
@assert loaded.arch.dstar == ProteinCoLoc.NPE_DSTAR
@assert loaded.arch.depth == ProteinCoLoc.NPE_DEPTH
@assert loaded.arch.width == ProteinCoLoc.NPE_WIDTH
@info "high-capacity flow persisted; summary net unchanged — spike-012 training complete"
