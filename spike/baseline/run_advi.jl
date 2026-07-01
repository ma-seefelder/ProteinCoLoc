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

# spike/baseline/run_advi.jl --- 04-04: the modern AdvancedVI ADVI ground-truth
# baseline (NPE-02/NPE-03). Runs ENTIRELY inside the isolated `spike/baseline/`
# env (D-01); Turing must never enter the pinned spike env and src/ is never
# edited (the @model is lifted read-only in model.jl).
#
# WHAT THIS PORTS (src/bayes.jl L308-329, READ-ONLY port target):
#   REMOVED  q = vi(m, ADVI(num_latent, iter))         # ADVI(n,iter) + 2-arg vi both gone
#   REMOVED  DynamicPPL.syms(VarInfo(m))               # gone; params read by @varname now
#   SKIPPED  sample(m, Prior(), posterior_samples)     # 100k prior chain = dead cost (INFER-3,
#                                                       #   D-08: the ADVI clock is vi() ONLY)
#   MODERN   result = vi(m, q_meanfield_gaussian, ITER; adtype=AutoForwardDiff())
#            draws  = rand(result, N)                  # constrained VarNamedTuple draws
#
# VERSIONS RESOLVED IN THIS ENV (verified, NOT the RESEARCH-anticipated ~0.7):
#   Turing 0.45.0, AdvancedVI 0.6.2, DynamicPPL 0.41.8, Bijectors 0.15.24, ForwardDiff 1.4.1.
#
# CONSTRAINED SPACE (Pitfall 3 / T-04-CONSTR): in Turing 0.45 `rand(result::VIResult, N)`
# returns `DynamicPPL.VarNamedTuple`s of RAW (constrained) parameter values -- it applies
# the inverse link/bijector internally via `result.ldf` before handing back marginals. So
# the μ draws arrive ALREADY in the model's `Truncated(...,-1,1)` support; no manual
# `Bijectors` inverse transform is needed on this API. `_assert_constrained!` HARD-CHECKS
# that invariant per pair -- a violation would mean a genuine transform/model bug, not a
# missing unconstrain step, and must surface rather than be silently ghat-clamped.
#
# DECOUPLING: `include(model.jl)` pulls the read-only src/ @model lift, the 04-01
# resimulate hook (now active -- its simulator deps were added to THIS isolated env), and
# the frozen ghat / holdout generator. All new state lives under spike/baseline/.

using Turing        # vi, q_meanfield_gaussian, q_fullrank_gaussian, AutoForwardDiff, @varname
using JLD2          # jldsave / jldopen (atomic artifact writer, cache.jl idiom)
using Statistics    # quantile (ρ-space 90% intervals), extrema
using Dates         # timestamp in the artifact meta

# model.jl: build_coloc_model + the lifted coloc_model + (now-active) resimulate hook,
# which transitively brings resimulate_holdout / load_holdout / the simulator / seeding.
include(joinpath(@__DIR__, "model.jl"))

# ghat (μ→ρ, frozen SIM-02 map) and _write_holdout (holdout generator) are needed for the
# ρ-space precompute and to materialize the reserved holdout the baseline re-simulates
# from. Guard on isdefined so a re-include (or a prior transitive include) is a no-op.
isdefined(Main, :ghat)           || include(joinpath(@__DIR__, "..", "simulator", "ghat.jl"))
isdefined(Main, :_write_holdout) || include(joinpath(@__DIR__, "..", "data", "generate.jl"))

# ============================ Pre-registered baseline defaults ====================
# A5 (RESEARCH): the ADVI baseline configuration, fixed BEFORE any reported run so the
# NPE-02 RMSE / NPE-03 speedup comparison is not tuned post-hoc.
const CHANNELS          = [1, 2]                 # 2-channel selector (fixed spike scope)
const NUM_PATCHES       = 8                      # fixed 8×8 patch grid (CLAUDE.md constraint)
const ITER              = 1000                   # vi() max iterations (A5)
const N                 = 10_000                 # posterior draws per pair (A5)
const ADTYPE            = AutoForwardDiff()      # safe default AD for this small model (A5)
const FAMILY            = q_meanfield_gaussian   # mean-field VI; q_fullrank_gaussian is the
                                                 #   documented A5 sensitivity option
const NPE_MASTER_SEED   = 0xC0FFEE               # pre-registered seed (04-01 constant)
const HOLDOUT_H         = 20                     # reserved holdout stacks (NPE-02 ≥20 → 10 pairs)
const DEFAULT_HOLDOUT_DIR = joinpath(@__DIR__, "holdout")   # baseline-local holdout fixture
const ARTIFACT_PATH     = joinpath(@__DIR__, "advi_artifact.jld2")  # the D-02 cross-env hand-off
const CONSTR_TOL        = 1e-6                    # numerical slack on the [-1,1] support check

# VarName handles for the two per-condition global means read out of each draw. The nested
# per-image/per-patch parameters are marginalized out; only μ_sample / μ_control map to ρ.
const _VN_MU_SAMPLE  = @varname(μ_sample)
const _VN_MU_CONTROL = @varname(μ_control)

# git ref of the tree the read-only @model was lifted from (artifact provenance, T-04-ARTIFACT).
const BAYES_GIT_REF = try
    strip(read(`git -C $(@__DIR__) rev-parse HEAD`, String))
catch
    "unknown"
end

# ============================ constrained-support guard ===========================
"""
    _assert_constrained!(v, name) -> v

Hard-assert that every draw in `v` lies in the model's constrained support `[-1,1]`
(within `CONSTR_TOL`). Pitfall 3 / T-04-CONSTR: `rand(VIResult, N)` on Turing 0.45 already
returns constrained draws, so this is an INVARIANT CHECK — a failure signals a real
transform/model bug, not a missing unconstrain step, and must halt rather than be masked by
`ghat`'s endpoint clamping.
"""
function _assert_constrained!(v::AbstractVector{<:Real}, name::AbstractString)
    lo, hi = extrema(v)
    @assert (lo >= -1 - CONSTR_TOL && hi <= 1 + CONSTR_TOL) (
        "ADVI $name draws fell outside the constrained support [-1,1]: extrema=($lo,$hi). " *
        "Turing 0.45 rand(VIResult) must yield constrained draws; a violation indicates a " *
        "bijector/model bug (Pitfall 3 / T-04-CONSTR).")
    return v
end

# ============================ single-pair ADVI ====================================
"""
    advi_pair(mci_sample, mci_control; channels=CHANNELS, num_patches=NUM_PATCHES,
              iter=ITER, posterior_samples=N, adtype=ADTYPE, family=FAMILY)
        -> (mu_sample, mu_control, wall_clock)

Run the modern-API ADVI on ONE (sample, control) pair of RAW `MultiChannelImage`s.

Builds the read-only lifted `coloc_model` via `build_coloc_model` (the exact object
`colocalization()` runs ADVI on), fits it with `vi(m, family, iter; adtype=...)`, and draws
`posterior_samples` constrained samples. Returns the full per-condition μ SAMPLE VECTORS
(length `posterior_samples`, so the spike can recompute Δρ by MC differencing) and the
`vi()`-only wall-clock (D-08: the prior chain is skipped; only `vi()` is timed).

The μ draws are asserted to lie in `[-1,1]` before return (T-04-CONSTR).
"""
function advi_pair(mci_sample, mci_control;
                   channels = CHANNELS, num_patches = NUM_PATCHES,
                   iter = ITER, posterior_samples = N,
                   adtype = ADTYPE, family = FAMILY)
    m = build_coloc_model(mci_sample, mci_control, channels, num_patches)

    # Time ONLY the vi() call (D-08). The removed 100k `sample(m, Prior(), …)` prior chain
    # (INFER-3 dead cost) is NOT run at all.
    t0  = time_ns()
    res = vi(m, family, iter; adtype = adtype, show_progress = false)
    wall_clock = (time_ns() - t0) / 1e9

    draws      = rand(res, posterior_samples)            # Vector{VarNamedTuple}, constrained
    mu_sample  = Float64[d[_VN_MU_SAMPLE]  for d in draws]
    mu_control = Float64[d[_VN_MU_CONTROL] for d in draws]

    _assert_constrained!(mu_sample,  "μ_sample")
    _assert_constrained!(mu_control, "μ_control")

    return (mu_sample = mu_sample, mu_control = mu_control, wall_clock = wall_clock)
end

# ============================ holdout materialization =============================
"""
    ensure_holdout(dir=DEFAULT_HOLDOUT_DIR; master_seed=NPE_MASTER_SEED, H=HOLDOUT_H) -> dir

Ensure a reserved-holdout `holdout.jld2` exists at `dir`, generating `H` stacks from the
DISJOINT `holdout_rng(master_seed, ·)` namespace (D-10) via the UNCHANGED simulator if
absent. Idempotent: an existing holdout is reused untouched so the baseline stays
reproducible and `resimulate_holdout` replays the identical byte-for-byte stacks. Returns
the holdout dir (the argument `run_all_pairs` / the artifact writer consume).
"""
function ensure_holdout(dir = DEFAULT_HOLDOUT_DIR;
                        master_seed = NPE_MASTER_SEED, H = HOLDOUT_H)
    isdir(dir) || mkpath(dir)
    isfile(joinpath(dir, "holdout.jld2")) || _write_holdout(dir, master_seed, H)
    return dir
end
