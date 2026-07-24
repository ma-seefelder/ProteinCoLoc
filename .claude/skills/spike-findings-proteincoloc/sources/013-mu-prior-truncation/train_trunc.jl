#########################################################################################
# Spike 013 --- the SIMULATOR-side μ-prior TRUNCATION lever: grid-8 retrain.
#
# QUESTION. The ρ_true ±0.99 prior atoms (~6.5% of draws; Spike 010 addendum) come from the
# Truncated-Cauchy μ-prior placing mass on μ the forward model cannot realize; `ghat` clamps those
# to ±0.99. TRUNCATING the μ-prior to the ghat-achievable support [GHAT_MU_MIN, GHAT_MU_MAX] =
# [-0.67976, 0.847149] makes ρ_true CONTINUOUS over [-0.99, 0.99] (no clamp atoms) — so ρ_true SBC
# is valid WITHOUT randomized ranks. Does it work at the root, and what does it COST near ±0.9?
#
# THE LEVER (SIMULATOR-side; injected datagen, src/ NOT edited):
#     μ* ~ Truncated(Cauchy(0,0.3), -1, 1)                         [baseline: ProteinCoLoc.sample_prior]
#   → μ* ~ Truncated(Cauchy(0,0.3), GHAT_MU_MIN, GHAT_MU_MAX)      [this run: sample_prior_trunc]
# The 6 nuisances reuse the SAME src/ prior constants unchanged (see trunc_datagen.jl).
#
# EVERYTHING ELSE IDENTICAL TO amended_v2 (scripts/train_grid8_amended_v2.jl):
#   • the F5 image-size MIXTURE, READ DIRECTLY from the frozen test/gate/gate_consts_8_v2.jl,
#   • n_pairs = 50_000, use_gpu = false, default architecture + optimizer recipe + epochs,
#   • the ONLY change is the μ-prior (⇒ ρ_true support), injected via the `datagen` seam.
#
# NEW artifact root: artifacts/spike013_trunc/grid_8/. amended_v2, spike012_highcap and every frozen
# bundle are NEVER touched (sha256 verified before/after in the spike README).
#
# SEED. master_seed = a FRESH DEV seed, asserted disjoint from every PROD/VAL/NPE/DEFAULT seed and
# every dev seed burned so far. PROD_SEED_V2 is NOT consumed. TRAINING ONLY — no gate arm is run.
#
# PROVENANCE. Injecting `datagen` records the training imsize provenance as `:unknown` (the pipeline
# cannot see inside the closure). That is expected and handled: the eval uses the SAME known mixture
# explicitly (read from the frozen consts) rather than reading it back.
#
# Run:  julia --project -t auto .planning/spikes/013-mu-prior-truncation/train_trunc.jl
#########################################################################################
using ProteinCoLoc
include(joinpath(@__DIR__, "trunc_datagen.jl"))

# The frozen amended pre-registration, in an ISOLATED module: we want ONLY its imsize mixture and
# the forbidden-seed sets. Nothing is run against it.
module GateV2
    include(joinpath(@__DIR__, "..", "..", "..", "test", "gate", "gate_consts_8_v2.jl"))
end

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "artifacts", "spike013_trunc"))

# The F5 mixture, taken verbatim from the frozen file (NOT retyped).
const IMSIZE_SET     = GateV2.SBC_IMSIZE_SET
const IMSIZE_WEIGHTS = GateV2.SBC_IMSIZE_WEIGHTS
@assert length(IMSIZE_SET) == length(IMSIZE_WEIGHTS)
@assert isapprox(sum(IMSIZE_WEIGHTS), 1.0; atol = 1e-12)
@assert (1376, 1028) in IMSIZE_SET          # the D-08 real-data anchor
@assert !((256, 256) in IMSIZE_SET)         # the 07-05 compute-budget artifact, excluded

# ---- FRESH DEV training seed, asserted disjoint from EVERYTHING pre-registered / already burned ---
const TRAIN_SEED = UInt64(0x00DEC0DE)       # fresh; NOT DEFAULT_MASTER_SEED (0x1), NOT any dev seed used
const NAMED_FORBIDDEN = UInt64[
    0x5BC0FFEE,      # VAL_MASTER_SEED
    0xC0FFEE,        # NPE_MASTER_SEED
    0x1,             # DEFAULT_MASTER_SEED
    0xDE7C0DE, 0xDE7C0DE2, 0x0DE7C0D3,                          # burned dev seeds
    0x00B7A771, 0x00B7A772, 0x00CA11B0, 0x00CA11B1,             # burned dev seeds
    0x00A1CE01, 0x00CA11B2, 0x00CA11B3,                         # burned dev seeds
]
const FORBIDDEN = Set{UInt64}(vcat(NAMED_FORBIDDEN,
                                   collect(values(GateV2.PROD_SEED_V2)),
                                   collect(values(GateV2.PROD_SEED))))
@assert !(TRAIN_SEED in FORBIDDEN) "TRAIN_SEED $(repr(TRAIN_SEED)) collides with a forbidden seed"
@assert !(TRAIN_SEED in values(GateV2.PROD_SEED_V2)) "TRAIN_SEED collides with a PROD_SEED_V2 value"

@info "spike-013 μ-prior-truncation grid-8 retrain (SIMULATOR-side lever, everything else = amended_v2)" ROOT IMSIZE_SET IMSIZE_WEIGHTS mu_prior = MU_PRIOR_TRUNC TRAIN_SEED = repr(TRAIN_SEED) nthreads = Threads.nthreads()

# The injected truncated-prior datagen closure (the ONLY delta vs amended_v2's real datagen).
trunc_datagen = () -> generate_samples_trunc(50_000; grid = 8, master_seed = TRAIN_SEED,
                                             imsize_set = IMSIZE_SET, imsize_weights = IMSIZE_WEIGHTS)

t0 = time()
b = ProteinCoLoc._train_grid_pipeline(8;
        n_pairs        = 50_000,                # informational (datagen is injected)
        use_gpu        = false,                 # identical to amended_v2
        artifacts_root = ROOT,                  # NEW store — every prior artifact untouched
        datagen        = trunc_datagen,         # THE SIMULATOR-SIDE LEVER (truncated μ-prior)
        skip_if_done   = true,
        verbose        = true)                  # per-epoch train/val risk → the training curve
elapsed = (time() - t0) / 60
@info "training done" minutes = elapsed grid = b.grid

# --- provenance READBACK: injected datagen ⇒ :unknown (expected; handled explicitly in the eval) --
npe_path = joinpath(ROOT, "grid_8", "npe_8.jld2")
loaded   = ProteinCoLoc.load_estimator(npe_path)
prov     = ProteinCoLoc.training_imsize_provenance(loaded)
@info "training_imsize_provenance readback (injected datagen ⇒ :unknown by design)" prov
@assert prov.imsize_source === :injected_datagen
@assert prov.recorded == false

# --- confirm the persisted architecture is the amended_v2 DEFAULT (only the prior changed) --------
@info "persisted architecture (must equal amended_v2 defaults; only the μ-prior differs)" arch = loaded.arch
@assert loaded.arch.num_coupling_layers == ProteinCoLoc.NPE_COUPLING
@assert loaded.arch.flow_depth == ProteinCoLoc.NPE_FLOW_DEPTH
@assert loaded.arch.flow_width == ProteinCoLoc.NPE_FLOW_WIDTH
@assert loaded.arch.dstar == ProteinCoLoc.NPE_DSTAR
@assert loaded.arch.depth == ProteinCoLoc.NPE_DEPTH
@assert loaded.arch.width == ProteinCoLoc.NPE_WIDTH
@info "spike-013 training complete — default architecture, truncated μ-prior pool"
