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

# spike/comparator/run_comparator.jl --- CMP-08: seeded entry point (D-12)
#
# THE callable deliverable of Phase 9. A single seeded entry point that composes the
# whole cross-method comparator harness behind one function:
#
#   build_shared_inputs  (seeded, regime-labelled shared MultiChannelImages, D-02/D-03)
#     -> build_table      (classical estimator battery + divergence/traffic-light, D-09/D-10)
#     -> tapqir_anchor    (OPTIONAL Tapqir sanity anchor, NEVER gating the run, D-07/D-08)
#     -> write_table      (atomic content-addressed CSV+JLD2 artifact, D-05)
#     -> comparator_audit (BayesInteractomics-style band-count + KS-uniformity report, D-11)
#
# REPRODUCIBILITY (D-12): everything downstream of the Philox `MASTER_SEED` is a pure
# function of `(master_seed, index)` via the `sample_rng`/`costes_rng` streams, so the
# emitted table + content-hash artifact dir are BIT-REPRODUCIBLE across runs and thread
# counts. Two same-seed `run_comparator` calls resolve to the IDENTICAL artifact dir and
# the IDENTICAL `Costes_p` column (the 09-06 D-13 determinism gate).
#
# GRACEFUL DEGRADATION (D-08 / T-09-16): the optional Tapqir anchor NEVER fails the run.
# `tapqir_anchor` self-catches every error and returns `status=:skipped`; this wrapper
# additionally defends with its own try/catch, so a Python/Conda hiccup can only ever
# yield a `:skipped` status stored in the run meta — the classical table always ships.
#
# DECOUPLING (hard constraint, CLAUDE.md): the FROZEN `src/` math (patch/correlation/
# MultiChannelImage) is reached ONLY transitively through the read-only `include()` in
# `spike/contract.jl`. This file contains NO `include(...src/...)` of its own — a second
# src/ include would double-define the frozen functions and break the decoupling proof.
# All code lives under `spike/`; no `src/`, no root Project.toml is touched.

# --- ORDER MATTERS: contract.jl FIRST (brings build_mci + the read-only src/ coupling
#     into scope), then the comparator modules in dependency order. Mirrors
#     spike/data/generate.jl:42-51's guarded-include idiom exactly. Each guard skips the
#     include when a sibling already pulled the symbol into scope, so this file loads both
#     standalone and after runtests.jl already loaded a comparator module.
isdefined(@__MODULE__, :build_mci)            || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :MASTER_SEED)          || include(joinpath(@__DIR__, "config.jl"))
isdefined(@__MODULE__, :manders)              || include(joinpath(@__DIR__, "classical.jl"))
isdefined(@__MODULE__, :build_shared_inputs)  || include(joinpath(@__DIR__, "inputs.jl"))
isdefined(@__MODULE__, :build_table)          || include(joinpath(@__DIR__, "table.jl"))
isdefined(@__MODULE__, :comparator_audit)     || include(joinpath(@__DIR__, "audit.jl"))
isdefined(@__MODULE__, :tapqir_anchor)        || include(joinpath(@__DIR__, "tapqir_bridge.jl"))

"""
    run_comparator(; master_seed=MASTER_SEED, n_per_regime=8, imsize=(128,128),
                   inputs=nothing, outdir=mktempdir(), run_tapqir=true) -> NamedTuple

The single seeded Phase-9 entry point (D-12). Builds the shared inputs, runs the full
classical estimator battery + divergence/traffic-light assembly, optionally runs the
Tapqir sanity anchor, writes the atomic content-addressed artifact, and renders the
audit. Returns

    (table::DataFrame,        # the per-input comparison table (build_table)
     artifact_dir::String,    # the content-hash artifact dir (write_table)
     audit::String,           # the BayesInteractomics-style Markdown audit (comparator_audit)
     tapqir::NamedTuple)      # (status, value, reason) from the optional anchor

Inputs:
  * when `inputs === nothing`, the shared set is GENERATED seeded via
    `build_shared_inputs(; master_seed, n_per_regime, imsize)` — `n_per_regime` inputs
    per `(:coloc, :random, :exclusion)` regime, each a pure function of
    `(master_seed, index)` (D-02);
  * when `inputs` is a `Vector{<:MultiChannelImage}`, the D-03 passthrough is used
    (regime `:unknown`, no fabricated ground truth) and the estimator battery is still
    seeded off the `master_seed` kwarg so the Costes null stays reproducible.

`run_tapqir=false` skips the (out-of-process) Tapqir probe entirely and records
`status=:skipped` — used by the determinism gate to keep the run subprocess-free and
fully seeded. A `:skipped` Tapqir status is a fully green outcome (D-08 / T-09-16).

Two calls with the same `master_seed`, `n_per_regime`, `imsize`, and `run_tapqir` produce
a BYTE-IDENTICAL table and the IDENTICAL `artifact_dir` (content-hash), independent of
thread count (D-12).
"""
function run_comparator(; master_seed::Integer = MASTER_SEED,
                        n_per_regime::Int = 8,
                        imsize::Tuple{Int,Int} = (128, 128),
                        inputs = nothing,
                        outdir::AbstractString = mktempdir(),
                        run_tapqir::Bool = true)
    # --- 1. Shared inputs (seeded generation, or the D-03 passthrough) ------------
    shared = inputs === nothing ?
        build_shared_inputs(; master_seed = master_seed,
                            n_per_regime = n_per_regime, imsize = imsize) :
        build_shared_inputs(inputs)

    # --- 2. Classical estimator battery + divergence/traffic-light table ----------
    # Seed the Costes null off the run's master_seed (concrete even for the passthrough,
    # whose shared.master_seed is `missing`), so Costes_p stays reproducible either way.
    df = build_table(shared; master_seed = master_seed)

    # --- 3. Optional Tapqir sanity anchor (NEVER gates the run, D-07/D-08) ---------
    # tapqir_anchor already self-catches; the extra try/catch is defence-in-depth so a
    # skip can only ever record a status, never raise out of run_comparator (T-09-16).
    tapqir = if run_tapqir
        try
            tapqir_anchor(; tol = TAPQIR_TOL)
        catch e
            (status = :skipped, value = missing, reason = sprint(showerror, e))
        end
    else
        (status = :skipped, value = missing, reason = "run_tapqir=false")
    end

    # --- 4. Atomic content-addressed artifact (CSV + JLD2) ------------------------
    # The hashed config carries the master seed, the input count, and the Tapqir status
    # (all deterministic), so two same-seed runs resolve to the identical artifact dir.
    config = (master_seed = master_seed,
              n_inputs = nrow(df),
              tapqir_status = tapqir.status)
    artifact_dir = write_table(df, config; outdir = outdir)

    # --- 5. Audit report (band counts + Costes-p KS-uniformity) -------------------
    audit = comparator_audit(df)
    # Persist the audit alongside the artifact for the human-readable deliverable. Its
    # timestamp is NOT part of the content hash (hash is over CMP_SRC_FILES + config), so
    # this write never perturbs the determinism of the artifact dir.
    write(joinpath(artifact_dir, "audit.md"), audit)

    return (table = df, artifact_dir = artifact_dir, audit = audit, tapqir = tapqir)
end
