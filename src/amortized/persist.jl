#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder

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

#############################################################################################
# src/amortized/persist.jl --- CPU-resident Flux.state persistence (PROD-01, Pitfall 4, T-7-07).
#
# The atomic `.tmp` → reopen-integrity-`@assert` → `mv(...; force=true)` wrapper is carried
# VERBATIM from the spike (spike/npe/train_npe.jl:160-180 save_npe; spike/data/cache.jl write_*),
# but the PAYLOAD is migrated (Pitfall 4 / T-7-07):
#
#   whole-object  `estimator = result.estimator`   (spike)
#        ↓
#   CPU-resident  `model_state = Flux.state(cpu(est))` + architecture metadata   (shipped)
#
# WHY (T-7-07 reproducibility trust boundary): a GPU-trained estimator carries GPU-resident
# parameter arrays and device-specific object graphs. `Flux.state(cpu(est))` extracts a
# device-independent nested-array snapshot; `load_*` rebuilds the architecture with
# `build_estimator`/`build_ratio_estimator` and applies `Flux.loadmodel!`, so the reloaded net is
# CPU-resident and the pre-registered CPU ship-gate + shipped inference are byte-reproducible
# regardless of the training hardware. This also narrows the deserialization surface (T-7-01):
# parameter arrays, not an arbitrary object graph.
#
# The frozen OOD null fits (density + noise + PP references) are plain arrays/NamedTuples (with a
# LinearAlgebra `Cholesky`), NOT a Flux model, so they persist DIRECTLY through the SAME atomic
# wrapper. Every artifact carries a `schema_version` key.
#
# INTERRUPTION SAFETY: a crash before the `mv` leaves only a discardable `.tmp`, never a
# half-written artifact — the SKIP-IF-DONE integrity check (`_estimator_ok`/`_ood_nulls_ok`)
# treats a torn file as not-done.
#############################################################################################

import Flux
import JLD2

# On-disk schema tags (bump on any breaking layout change).
const ESTIMATOR_SCHEMA = 1
const OOD_NULLS_SCHEMA = 1

# ============================ NPE estimator (CPU-resident Flux.state) ========================

"""
    save_estimator(path, est, θzt, zt, arch; meta = (;)) -> String

Atomically persist a trained NPE estimator CPU-RESIDENT (Pitfall 4 / T-7-07): write
`Flux.state(cpu(est))` (a device-independent nested-array snapshot — NOT the whole
`PosteriorEstimator` object), the architecture metadata `arch` (`d_in`, `D`, `dstar`, `depth`,
`width`, coupling/flow params) needed to rebuild the net, and the FROZEN `θzt`/`zt` transforms to
`path.tmp`, reopen to integrity-`@assert` the `model_state` key, then `mv(...; force=true)`
(filesystem-atomic). Returns `path`.
"""
function save_estimator(path, est, θzt, zt, arch; meta = (;))
    mkpath(dirname(path))
    tmp   = path * ".tmp"
    state = Flux.state(Flux.cpu(est))                # CPU-resident nested-array snapshot
    JLD2.jldsave(tmp;
        schema_version  = ESTIMATOR_SCHEMA,
        model_state     = state,
        arch            = arch,
        theta_transform = θzt,
        zt              = zt,
        meta            = meta)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "model_state") "save_estimator: integrity check failed, $tmp missing model_state"
    end
    mv(tmp, path; force = true)                      # atomic commit
    return path
end

"""
    load_estimator(path) -> NamedTuple

Reload an NPE estimator persisted by `save_estimator`, CPU-RESIDENT: rebuild the architecture via
`build_estimator(arch.d_in, arch.D; ...)` (fresh init) then `Flux.loadmodel!(est, model_state)` to
copy the frozen parameters into it. Returns `(estimator, θzt, zt, arch, meta)`; the estimator is
on CPU, so `posterior_for(estimator, Z; use_gpu=false)` is device-independent.
"""
function load_estimator(path)
    isfile(path) || error("load_estimator: no file at $path")
    d = JLD2.load(path)
    haskey(d, "model_state") || error("load_estimator: $path missing model_state key")
    a   = d["arch"]
    est = build_estimator(a.d_in, a.D; dstar = a.dstar, depth = a.depth, width = a.width,
                          num_coupling_layers = a.num_coupling_layers,
                          flow_depth = a.flow_depth, flow_width = a.flow_width)
    Flux.loadmodel!(est, d["model_state"])
    return (estimator = est, θzt = d["theta_transform"], zt = d["zt"],
            arch = a, meta = get(d, "meta", nothing))
end

"""
    _estimator_ok(path) -> Bool

SKIP-IF-DONE integrity predicate for a persisted estimator: `true` iff `path` exists and reopens
with the `model_state` + `arch` keys present. A torn `.tmp`-turned-final (any open/read error) is
reported not-done so the pipeline recomputes rather than loading a corrupt artifact.
"""
function _estimator_ok(path)
    isfile(path) || return false
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "model_state") && haskey(f, "arch")
        end
    catch e
        e isa InterruptException && rethrow()
        return false
    end
end

# ============================ NRE ratio estimator (CPU-resident) =============================

"""
    save_ratio(path, ratio; meta = (;)) -> String

Atomically persist the trained ratio HANDLE `ratio` (as returned by `train_ratio`) CPU-RESIDENT:
`Flux.state(cpu(ratio.estimator))` + the topology metadata (`input_dim`, `num_summaries`,
`summary_width`) needed to rebuild it + the MEASURED `log_prior_odds` (Pitfall 5) + the D-07
`split_threshold`, through the atomic `.tmp`→integrity-`@assert`→`mv` wrapper. Returns `path`.
"""
function save_ratio(path, ratio; meta = (;))
    mkpath(dirname(path))
    tmp   = path * ".tmp"
    state = Flux.state(Flux.cpu(ratio.estimator))
    JLD2.jldsave(tmp;
        schema_version  = ESTIMATOR_SCHEMA,
        model_state     = state,
        input_dim       = ratio.input_dim,
        num_summaries   = ratio.num_summaries,
        summary_width   = ratio.summary_width,
        log_prior_odds  = ratio.log_prior_odds,
        split_threshold = RATIO_SPLIT_THRESHOLD,
        meta            = meta)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "model_state") "save_ratio: integrity check failed, $tmp missing model_state"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    load_ratio(path) -> NamedTuple

Reload a ratio handle persisted by `save_ratio`, CPU-RESIDENT: rebuild via
`build_ratio_estimator(input_dim; num_summaries, width)` then `Flux.loadmodel!`. Returns
`(estimator, log_prior_odds, num_summaries, summary_width, input_dim)` — the same handle shape
`train_ratio` returns, so `amortized_log_bf(handle.estimator, …, handle.log_prior_odds)` runs.
"""
function load_ratio(path)
    isfile(path) || error("load_ratio: no file at $path")
    d = JLD2.load(path)
    haskey(d, "model_state") || error("load_ratio: $path missing model_state key")
    est = build_ratio_estimator(d["input_dim"];
                                num_summaries = d["num_summaries"], width = d["summary_width"])
    Flux.loadmodel!(est, d["model_state"])
    return (estimator = est, log_prior_odds = d["log_prior_odds"],
            num_summaries = d["num_summaries"], summary_width = d["summary_width"],
            input_dim = d["input_dim"])
end

function _ratio_ok(path)
    isfile(path) || return false
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "model_state") && haskey(f, "input_dim")
        end
    catch e
        e isa InterruptException && rethrow()
        return false
    end
end

# ============================ OOD null fits (plain data, atomic) =============================

"""
    save_ood_nulls(path, ood_nulls; meta = (;)) -> String

Atomically persist the frozen OOD null fits `ood_nulls` (a NamedTuple carrying at least
`:density = (cont, μS, C)`; optionally `:noise`, `:zref`, `:thr`) through the SAME atomic
`.tmp`→integrity-`@assert`→`mv` wrapper. The nulls are plain arrays / `NamedTuple`s (plus a
LinearAlgebra `Cholesky`), NOT a Flux model, so they persist DIRECTLY — no `Flux.state` needed.
Returns `path`.
"""
function save_ood_nulls(path, ood_nulls; meta = (;))
    mkpath(dirname(path))
    tmp = path * ".tmp"
    JLD2.jldsave(tmp;
        schema_version = OOD_NULLS_SCHEMA,
        ood_nulls      = ood_nulls,
        meta           = meta)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "ood_nulls") "save_ood_nulls: integrity check failed, $tmp missing ood_nulls"
    end
    mv(tmp, path; force = true)
    return path
end

"""
    load_ood_nulls(path) -> NamedTuple

Reload the frozen OOD null fits persisted by `save_ood_nulls`. Integrity-checks the `ood_nulls`
key on open; returns the stored NamedTuple unchanged (the `maha_score`/`noise_score`/`ood_verdict`
scorers consume it directly).
"""
function load_ood_nulls(path)
    isfile(path) || error("load_ood_nulls: no file at $path")
    d = JLD2.load(path)
    haskey(d, "ood_nulls") || error("load_ood_nulls: $path missing ood_nulls key")
    return d["ood_nulls"]
end

function _ood_nulls_ok(path)
    isfile(path) || return false
    try
        return JLD2.jldopen(path, "r") do f
            haskey(f, "ood_nulls")
        end
    catch e
        e isa InterruptException && rethrow()
        return false
    end
end
