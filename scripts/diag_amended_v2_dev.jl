#########################################################################################
# scripts/diag_amended_v2_dev.jl --- DIAGNOSTIC ONLY. NOT a gate run.
#
# A development-time smoke of the amended-regime grid-8 bundle
# (`artifacts/amended_v2/grid_8/`): the net loads CPU-resident, its persisted training image-size
# provenance reads back as the frozen `SBC_IMSIZE_SET`/`SBC_IMSIZE_WEIGHTS`, and one
# draw→simulate→CPU-infer step at EACH of the four mixture image sizes produces finite posterior
# draws. Nothing here is scored against any pre-registered threshold.
#
# SEED DISCIPLINE (07-GATE-AMENDMENT protocol §6.3/§6.6). This runs on a DEV seed that is asserted
# disjoint from EVERY seed the project has burned or reserved: `PROD_SEED_V2[G]` (the single
# confirmatory amended-gate stream — untouched by this script), every v1 `PROD_SEED[G]`,
# `VAL_MASTER_SEED`, `NPE_MASTER_SEED`, `DEFAULT_MASTER_SEED`, and the two F2-remedy DEV seeds.
#
# THE AMENDED SHIP-GATE IS NOT INVOKED HERE. It is bound to exactly ONE run on `PROD_SEED_V2[8]`.
#
# Run:  julia --project -t auto scripts/diag_amended_v2_dev.jl
#########################################################################################
using ProteinCoLoc
using Random
import Statistics: median, std
import Random123: Philox4x

const GATEDIR = normpath(joinpath(@__DIR__, "..", "test", "gate"))

# Load the FROZEN amended pre-registration (constants only) then the gate harness, so the
# diagnostic exercises exactly the code path the confirmatory run will take.
include(joinpath(GATEDIR, "gate_consts_8_v2.jl"))
include(joinpath(GATEDIR, "harness.jl"))

# --- the DEV seed, and its disjointness from every reserved stream ---------------------------
const DIAG_DEV_SEED = 0x0000_0000_0DE7_C0D3      # ≠ 0x0DE7C0DE / 0xDE7C0DE2 (F2 remedy dev seeds)
const RESERVED = (UInt64(0), NPE_MASTER_SEED, VAL_MASTER_SEED, DEFAULT_MASTER_SEED,
                  DEV_SEEDS..., values(PROD_SEED)..., values(PROD_SEED_V2)...)
@assert !(DIAG_DEV_SEED in RESERVED) "DIAG_DEV_SEED collides with a reserved stream"
@assert DIAG_DEV_SEED != PROD_SEED_V2[8]
@info "DIAGNOSTIC (DEV seed, NOT the gate)" dev_seed = repr(DIAG_DEV_SEED) reserved_n = length(RESERVED)

const ROOT = normpath(joinpath(@__DIR__, "..", "artifacts", "amended_v2"))
m = ProteinCoLoc.load_estimator(joinpath(ROOT, "grid_8", "npe_8.jld2"))

prov = ProteinCoLoc.training_imsize_provenance(m)
@info "training_imsize_provenance" prov
@assert prov.recorded
@assert Tuple(prov.imsize_set) == SBC_IMSIZE_SET
@assert Tuple(prov.imsize_weights) == SBC_IMSIZE_WEIGHTS

# The binding invariant the amended gate will enforce, checked here read-only.
chk = assert_imsize_provenance(m)
@info "assert_imsize_provenance (the gate's refusal check)" ok = chk.ok status = chk.status
@assert chk.ok && chk.status === :ok

# --- one draw→simulate→CPU-infer step at EACH mixture size (DEV stream) -----------------------
Random.seed!(DIAG_DEV_SEED)                      # posterior draws come from the GLOBAL stream
rng = Philox4x(UInt64, (DIAG_DEV_SEED, UInt64(0)))
for isz in SBC_IMSIZE_SET
    t = draw_simulate_infer(m, rng; G = 8, imsize = isz, N = 64)
    ρ = @view t.draws[1, :]
    @info "DIAGNOSTIC draw" imsize = t.imsize realised = size(t.Z) all_finite = all(isfinite, t.draws) ρ_true = round(t.θ.ρ_true; digits = 4) ρ_post_median = round(median(ρ); digits = 4) ρ_post_sd = round(std(ρ); digits = 4)
    @assert all(isfinite, t.draws)
    @assert all(-1 .<= ρ .<= 1)                  # the bounded (logit-of-prior-box) θ-space holds
end

# A per-draw MIXTURE sample, to confirm the shared sampler runs on the real set.
let szs = [gate_imsize(rng) for _ in 1:200]
    @info "DIAGNOSTIC realised mixture over 200 draws" counts = realised_imsize_counts(szs)
    @assert all(s -> s in SBC_IMSIZE_SET, szs)
end

@info "DIAGNOSTIC complete — NO gate arm was run, PROD_SEED_V2[8] is untouched"
