#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/run_slow_calibration.jl --- the SLOW-TIER reported-scale calibration run.
#
#   julia --project=. --threads=auto test/gate/run_slow_calibration.jl [--out DIR]
#   P15_SLOW_M=1000 julia --project=. --threads=auto test/gate/run_slow_calibration.jl --out slow-report
#
# ============================ WHAT THIS TIER DETECTS, AND WHY IT EXISTS =======================
# THE FAST TIER (`ci_golden.jl`, `.github/workflows/ci.yml`) DETECTS **CHANGE**. THIS TIER DETECTS
# **MISCALIBRATION**. That split is not a nicety: at the fast tier's `P15_GOLDEN_M = 8` draws the
# ECE null mean is ~0.12, an order of magnitude above any miscalibration worth detecting, so the
# fast gate's numbers carry no calibration information at all — they are a bit-level fingerprint.
# This script runs the SBC/coverage arm at REPORTED scale (`P15_SWEEP_M` by default), which is the
# scale at which an ECE means something. `P15_GOLDEN_HONESTY` in `p15_consts.jl` is the frozen
# wording of the same statement and is printed verbatim by every run of this script.
#
# ================================ WHAT IT DELIBERATELY DOES NOT DO ============================
# It does NOT run the 7-axis envelope sweep. See the header of
# `.github/workflows/calibration-slow.yml` for the full reasoning (4-core hosted runners, a 6 h job
# cap, a ~8.5 h sweep at 32 threads). The sweep is a one-off local run whose artifact is a committed
# report.
#
# It is DELIBERATELY NOT ROUTED THROUGH `run_gate.jl`'s CLI, for two independent reasons:
#   (1) `run_gate`'s default `artifacts_root` is the v1 `artifacts/` tree, and on a CI runner there
#       is NO in-repo `artifacts/` tree at all — the shipped net has to be resolved through
#       `ProteinCoLoc.estimator_for(8)`'s three integrity-checked paths instead of a path;
#   (2) under `SBC_REQUIRE_IMSIZE_PROVENANCE = true` `run_gate` REFUSES to run with status
#       `:provenance_mismatch` on any net whose training image-size distribution was not persisted.
#       That check is still COMPUTED and REPORTED here (see `slow_provenance`) — it is simply not
#       allowed to silently turn the calibration tier into a no-op.
# `write_gate_report` is reused from `run_gate.jl` (included for that function alone) so the report
# this tier uploads is written through the same atomic `.tmp` -> reopen-integrity -> `mv` wrapper
# every other gate report uses.
#
# ================================ THE ONE EXTERNALLY-INFLUENCED VALUE =========================
# `P15_SLOW_M` is the `workflow_dispatch` input, and it is the ONLY value in this phase that anyone
# outside the repository can influence. It is read from the ENVIRONMENT — never interpolated into a
# shell command line by the workflow — parsed with `tryparse` and HARD-BOUNDED to
# `P15_SLOW_M_MIN:P15_SLOW_M_MAX`. The check runs at LOAD time, before the package is even loaded,
# so an out-of-range value cannot burn runner minutes before being rejected (T-15-07, T-15-05).
#
# NO STDLIB BEYOND WHAT `[deps]` ALREADY CARRIES. `Printf`, `Dates` and `SHA` are stdlibs that are
# NOT `[deps]` of this package, and importing one aborts the whole suite under `Pkg.test()`'s
# sandbox (15-01 deviation 2, 15-06 deviations 2-3). Number formatting below therefore uses `round`
# + `rpad`, and the timestamp uses `Base.Libc.strftime`.

# =============================================================================================
# THE BOUNDS CHECK, FIRST — before `using`, before any include, before any runner minute is spent
# =============================================================================================

"Lower bound on the dispatch-supplied SBC draw count. Below this the ECE is pure null noise."
const P15_SLOW_M_MIN = 100

"""
Upper bound on the dispatch-supplied SBC draw count. This is a DENIAL-OF-SERVICE bound, not a
statistical one: at the reconciled cost model's ~0.656 s/draw at 32 threads, and roughly 4-6x that
on a 4-core hosted runner, M = 5000 already sits at the edge of the 6 h hosted job cap. A dispatch
input must not be able to ask for a job that cannot finish.
"""
const P15_SLOW_M_MAX = 5000

"""
    p15_parse_slow_m(raw) -> Union{Int,Nothing}

Parse and BOUNDS-CHECK the `P15_SLOW_M` environment value. `nothing` means "not supplied" (the
caller then falls back to the pre-registered `P15_SWEEP_M`); anything supplied must be an integer
inside `P15_SLOW_M_MIN:P15_SLOW_M_MAX` or this throws, naming the bound it violated.

Deliberately a pure function of its argument so the guard can be exercised without an environment.
"""
function p15_parse_slow_m(raw)
    raw === nothing && return nothing
    s = strip(String(raw))
    isempty(s) && return nothing
    m = tryparse(Int, s)
    m === nothing && error(
        "run_slow_calibration: P15_SLOW_M = $(repr(s)) is not an integer. The slow tier's SBC " *
        "draw count must be an integer in $(P15_SLOW_M_MIN):$(P15_SLOW_M_MAX).")
    (P15_SLOW_M_MIN <= m <= P15_SLOW_M_MAX) || error(
        "run_slow_calibration: P15_SLOW_M = $m is outside the admissible range " *
        "$(P15_SLOW_M_MIN):$(P15_SLOW_M_MAX). Below $(P15_SLOW_M_MIN) the ECE is indistinguishable " *
        "from its own null; above $(P15_SLOW_M_MAX) the job cannot finish inside the 6 h " *
        "GitHub-hosted wall-clock cap. This value reaches Julia through the environment only and " *
        "is never interpolated into a shell command.")
    return m
end

"The dispatch-supplied M, or `nothing` when the workflow input was left empty."
const P15_SLOW_M_REQUESTED = p15_parse_slow_m(get(ENV, "P15_SLOW_M", nothing))

# =============================================================================================
# Only now: the package and the gate machinery
# =============================================================================================

using ProteinCoLoc
import Artifacts

# Ordered, guarded includes: pre-registration -> harness -> sbc -> `write_gate_report`.
# Guarded so this file is both standalone-runnable and includable from a module that has already
# loaded the Phase-15 pre-registration (`p15_consts.jl` defines the same const names as every other
# gate pre-registration, so no two of them may share a module).
if !isdefined(@__MODULE__, :SBC_M)
    include(joinpath(@__DIR__, "p15_consts.jl"))
end
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "harness.jl"))
isdefined(@__MODULE__, :sbc_gate)            || include(joinpath(@__DIR__, "sbc.jl"))
isdefined(@__MODULE__, :write_gate_report)   || include(joinpath(@__DIR__, "run_gate.jl"))

# The Phase-15 bars exist ONLY in the Phase-15 pre-registration. If some OTHER consts file was
# loaded into this module first, the guard above would have SKIPPED the include and every bar below
# would silently be the wrong one — fail loudly instead (the same guard `ci_golden.jl` carries).
if !(isdefined(@__MODULE__, :GATE_CONSTS_VERSION) && GATE_CONSTS_VERSION == 15 &&
     isdefined(@__MODULE__, :P15_SWEEP_M))
    error("run_slow_calibration.jl: the loaded pre-registration is not the Phase-15 one " *
          "(GATE_CONSTS_VERSION = " *
          (isdefined(@__MODULE__, :GATE_CONSTS_VERSION) ? string(GATE_CONSTS_VERSION) : "undefined") *
          "). Load `p15_consts.jl` — into an isolated module if another consts file is already " *
          "present — before including this file.")
end

"The SBC draw count this run will use: the dispatch input when supplied, else the pre-registered M."
const SLOW_M = P15_SLOW_M_REQUESTED === nothing ? P15_SWEEP_M : P15_SLOW_M_REQUESTED

"The grid this tier gates. v2.0 ships the 8x8 reference grid only (D-10)."
const SLOW_G = 8

"The file this tier writes into `--out`. Uploaded by the workflow as the `slow-report` artifact."
const SLOW_REPORT_NAME = "slow_calibration_report.jld2"

"The document that governs re-blessing. A red slow tier is NEVER a re-bless reason."
const SLOW_REBLESS_DOC =
    ".planning/phases/15-calibration-operating-envelope-and-ci-gate/15-REBLESS.md"

# =============================================================================================
# OUTPUT DIRECTORY — anywhere but `artifacts/`
# =============================================================================================

"""
    slow_out_dir(args) -> String

The `--out DIR` output directory (default: the working directory), REFUSING any path that lies
inside an `artifacts` directory. The shipped bundles and the reported gate report live there; a CI
job must never be able to overwrite a reported artifact with a run whose M came from a dispatch box.
"""
function slow_out_dir(args = ARGS)
    dir = "."
    for i in 1:(length(args) - 1)
        args[i] == "--out" && (dir = args[i + 1])
    end
    abs = abspath(dir)
    any(p -> lowercase(p) == "artifacts", splitpath(abs)) && error(
        "run_slow_calibration: --out $(repr(dir)) resolves inside an `artifacts` directory " *
        "($abs). This tier writes a RUN REPORT, never a shipped artifact, and must not be able to " *
        "overwrite one.")
    return abs
end

# =============================================================================================
# THE NET
# =============================================================================================

"""
    slow_net() -> NamedTuple

Resolve the shipped 8x8 net and its content identity, returning
`(m, bundle, dir, resolved_via, artifact_hash, net_content_digest, meta_source)`.

`ProteinCoLoc.estimator_for(8)` is used rather than a path because it resolves in three
INTEGRITY-CHECKED steps (`src/registry.jl:126-161`): the content-addressed artifact store, then the
in-repo `artifacts/amended_v2/grid_8` dev dir WITH a tree-sha1 verification against `Artifacts.toml`,
then a lazy download of the release asset. Only the third path exists on a CI runner.

THE FIELD MISMATCH THE SHIM EXISTS FOR — do NOT "simplify" it away. `estimator_for` returns an
`EstimatorBundle` STRUCT with fields `(grid, npe, ratio, ood_nulls, zt, θzt, calibration)`
(`src/registry.jl:59-67`), while the gate harness reads the NamedTuple `(estimator, θzt, zt, ...)`
that `load_estimator` returns (`src/amortized/persist.jl:101-102`). The posterior net is `b.npe`;
there is no `b.estimator` field. This mirrors `p15_golden_net` in `ci_golden.jl` deliberately, so
the fast and slow tiers gate the SAME resolved object.

The content digest is a fold of `Base.hash` over sorted filenames and file bytes: mode-independent
(unlike `Pkg.GitTools.tree_hash`, which disagrees with itself across resolution paths) and
import-free. Non-cryptographic by construction — it catches accident and regression, not an
adversary; adversarial substitution is what the store's content addressing and `_verify_tree_sha1`
cover.

ONE FIELD DIFFERS FROM `ci_golden.jl`'S SHIM, AND ON PURPOSE: `meta`. `EstimatorBundle` has no
`meta` field (`src/registry.jl:59-67`), so a shim that hard-codes `meta = nothing` makes
`training_imsize_provenance` return `recorded = false` for EVERY net — a fact about the shim, not
about the net, and it would be printed here as if it were a finding. The persisted meta is
therefore re-read from the resolved `npe_G.jld2` and threaded in, so the F5 binding invariant
(gate mixture == training mixture) is actually CHECKABLE on this tier. The fast tier does not need
it: it asserts a digest, not a distributional theorem. When the meta cannot be read the field falls
back to `nothing` and `meta_source` says so, rather than the check silently reporting a mismatch.
"""
function slow_net()
    b = ProteinCoLoc.estimator_for(SLOW_G)

    toml     = ProteinCoLoc._artifacts_toml()
    expected = Artifacts.artifact_hash("grid_$(SLOW_G)", toml)
    expected === nothing &&
        error("run_slow_calibration: Artifacts.toml has no `grid_$(SLOW_G)` entry at $toml.")

    devdir = ProteinCoLoc._dev_bundle_dir(SLOW_G)
    dir, via = if Artifacts.artifact_exists(expected)
        (Artifacts.artifact_path(expected), :artifact_store)
    elseif isdir(devdir)
        (devdir, :dev_dir)
    else
        (nothing, :unresolved)
    end

    digest = if dir === nothing
        ""
    else
        h = zero(UInt)
        for f in sort(readdir(dir))
            h = hash(f, h)
            h = hash(read(joinpath(dir, f)), h)
        end
        string(h)
    end

    # The persisted training metadata, re-read from the SAME resolved directory (see the docstring).
    npe_path = dir === nothing ? nothing : joinpath(dir, "npe_$(SLOW_G).jld2")
    meta, meta_source = if npe_path !== nothing && isfile(npe_path)
        try
            (ProteinCoLoc.load_estimator(npe_path).meta, :npe_jld2)
        catch err
            @warn "run_slow_calibration: could not re-read the persisted training metadata; the " *
                  "F5 image-size provenance check will report `:unknown` rather than a mismatch" *
                  " (this is a MISSING check, not a passing one)" exception = err
            (nothing, :unreadable)
        end
    else
        (nothing, :no_resolved_dir)
    end

    return (m = (estimator = b.npe, θzt = b.θzt, zt = b.zt, arch = nothing, meta = meta),
            bundle = b, dir = dir, resolved_via = via, meta_source = meta_source,
            artifact_hash = string(expected), net_content_digest = digest)
end

# =============================================================================================
# THE ANCHOR (Tier-2 block `:P15_ECE_ANCHOR_MEASURED`, opened by 15-07)
# =============================================================================================

"""
    slow_anchor_ece(col) -> Union{Float64,Nothing}

The rung-0 ECE anchor for column `col`, or `nothing` when the `:P15_ECE_ANCHOR_MEASURED` Tier-2
block has not been appended to `p15_consts.jl` yet.

THE ANCHOR IS A TIER-2 MEASUREMENT, NOT A SHIPPED NUMBER (D-04a): it is `ECE(rung 0)` measured at
`P15_SWEEP_M` on the P15 stream for the net under test. Anchoring on the shipped
`gate_report_8.jld2` ECEs instead would declare a break on ~90 % of perfectly calibrated rungs,
because the shipped Δρ ECE sits BELOW its own null mean.

Several plausible shapes are accepted because 15-07 owns that block's layout and this file must not
dictate it: a `Vector`/`Tuple` indexed by column, an `AbstractDict` keyed by column index or by the
`SBC_PARAM_LABELS` label (as `String` or `Symbol`), or a `NamedTuple` keyed by the label. An anchor
that is present but in NONE of those shapes is a defect, not a pass — `slow_main` exits 4 on it
rather than quietly skipping the comparison.
"""
function slow_anchor_ece(col::Integer)
    isdefined(@__MODULE__, :P15_ECE_ANCHOR_MEASURED) || return nothing
    a = getfield(@__MODULE__, :P15_ECE_ANCHOR_MEASURED)
    lab = SBC_PARAM_LABELS[col]
    if a isa AbstractDict
        for k in (Int(col), lab, Symbol(lab), string(lab))
            haskey(a, k) && return float(a[k])
        end
        return :unreadable
    elseif a isa NamedTuple
        s = Symbol(lab)
        hasproperty(a, s) && return float(getproperty(a, s))
        return :unreadable
    elseif a isa Union{AbstractVector,Tuple}
        length(a) >= col || return :unreadable
        v = a[col]
        return v isa Real ? float(v) : :unreadable
    end
    return :unreadable
end

# =============================================================================================
# PRINTING
# =============================================================================================

"Fixed-width float, without `Printf` (a stdlib that is not a `[dep]`; see this file's header)."
_f(x, d = 5) = rpad(string(round(float(x); digits = d)), d + 4)

"`YYYY-MM-DDTHH:MM:SSZ`, without `Dates` (same reason)."
_slow_now() = Base.Libc.strftime("%Y-%m-%dT%H:%M:%SZ", time())

"""
    slow_print_header(io, out, threads)

The framing every run of this tier must carry: the fast/slow split stated VERBATIM from the frozen
`P15_GOLDEN_HONESTY`, plus the run's own parameters.
"""
function slow_print_header(io::IO, out::AbstractString, threads::Integer)
    println(io, "=" ^ 92)
    println(io, "SLOW-TIER CALIBRATION RUN — reported-scale SBC/coverage on the shipped grid-",
                SLOW_G, " net")
    println(io, "=" ^ 92)
    println(io)
    println(io, "WHAT THE TWO TIERS DO (frozen wording, P15_GOLDEN_HONESTY):")
    println(io, "  ", P15_GOLDEN_HONESTY)
    println(io)
    println(io, "  So, stated for this run: the FAST tier detects CHANGE. THIS tier detects")
    println(io, "  MISCALIBRATION. A green fast tier says nothing about calibration; this is the")
    println(io, "  run that does.")
    println(io)
    println(io, "  started      = ", _slow_now())
    println(io, "  M            = ", SLOW_M,
                P15_SLOW_M_REQUESTED === nothing ?
                    "  (pre-registered P15_SWEEP_M; no dispatch input supplied)" :
                    "  (from the P15_SLOW_M dispatch input, bounds-checked to " *
                    "$(P15_SLOW_M_MIN):$(P15_SLOW_M_MAX))")
    println(io, "  L            = ", SBC_L, "   bins = ", SBC_BINS,
                "   chi2_bins = ", _gate_chi2_bins(SBC_BINS))
    println(io, "  imsize       = ", repr(SBC_IMSIZE), " over ", SBC_IMSIZE_SET,
                " at ", SBC_IMSIZE_WEIGHTS)
    println(io, "  threads      = ", threads,
                threads < P15_SWEEP_MIN_THREADS ?
                    "  (BELOW P15_SWEEP_MIN_THREADS = $(P15_SWEEP_MIN_THREADS) — expected on a " *
                    "4-core hosted runner; this costs WALL-CLOCK ONLY, never a number)" : "")
    println(io, "  julia        = ", VERSION)
    println(io, "  gate seed    = ", prod_seed(SLOW_G),
                "   global seed = ", gate_global_seed(SLOW_G))
    println(io, "  out          = ", out)
    println(io)
    return nothing
end

"""
    slow_print_columns(io, g)

The eight per-column lines: label, ECE, traffic light, KS p, chi2 p, shrinkage. All eight are
REPORTED (`P15_REPORTED_COLUMNS`); only columns `P15_BREAK_COLUMNS` gate anything (D-05), and the
line marks which is which so the table cannot be read as eight gates.
"""
function slow_print_columns(io::IO, g)
    println(io, "PER-COLUMN SBC RESULTS  (M = ", g.M, ", L = ", g.L, ")")
    println(io, "  ", rpad("#", 3), rpad("column", 20), rpad("ECE", 9), rpad("light", 8),
                rpad("ks_p", 9), rpad("chi2_p", 9), rpad("shrinkage", 10), "role")
    println(io, "  ", "-" ^ 86)
    for (i, p) in enumerate(g.per_param)
        role = i in P15_BREAK_COLUMNS ? "GATING" : "reported only"
        println(io, "  ", rpad(i, 3), rpad(string(p.label), 20), _f(p.ece), rpad(string(p.light), 8),
                    _f(p.ks_p), _f(p.chi2_p), rpad(_f(p.shrinkage), 10), role)
    end
    println(io)
    println(io, "  The six nuisance marginals gate NOTHING (D-05). Their drift is an already-")
    println(io, "  disclosed named v2.0 limit; gating on it would re-fail the tool for a reason")
    println(io, "  already on the record.")
    if !isempty(g.vacuous_params)
        println(io, "  VACUOUS (posterior ~= prior, shrinkage > ", g.vacuous_cutoff, "): ",
                    join(g.vacuous_params, ", "),
                    " — such a column passes rank uniformity BY CONSTRUCTION while having learned")
        println(io, "  nothing. Reported, never gating.")
    end
    println(io)
    return nothing
end

"""
    slow_print_coverage(io, g)

The nominal-vs-empirical coverage curve for the two GATING columns. This is the diagnostic an ECE
summarises: a single number cannot say whether a column is over- or under-confident, and the curve
can.
"""
function slow_print_coverage(io::IO, g)
    println(io, "COVERAGE CURVES for the gating columns ", P15_BREAK_COLUMNS,
                "  (nominal -> empirical central-interval coverage)")
    for c in P15_BREAK_COLUMNS
        cov = sbc_coverage(collect(view(g.ranks, :, c));
                           levels = collect(P15_COVERAGE_LEVELS), L = g.L)
        println(io, "  column ", c, " (", SBC_PARAM_LABELS[c], ")")
        println(io, "    nominal   ", join((rpad(round(x; digits = 2), 6) for x in cov.nominal)))
        println(io, "    empirical ", join((rpad(round(x; digits = 3), 6) for x in cov.empirical)))
        println(io, "    residual  ",
                    join((rpad(round(e - n; digits = 3), 6)
                          for (n, e) in zip(cov.nominal, cov.empirical))))
    end
    println(io)
    return nothing
end

# =============================================================================================
# ENTRY POINT
# =============================================================================================

"""
    slow_main(args = ARGS) -> Int

Run the reported-scale SBC arm once, print everything, write the report, return an exit code:

  0  ran; every GATING column is within `p15_break_threshold(anchor)` of its anchor — OR no anchor
     has been measured yet, in which case the run is REPORTED and the absence is stated. An honest
     status, not a silent pass.
  1  a GATING column's ECE exceeds its break threshold — the miscalibration this tier exists for.
  4  the `:P15_ECE_ANCHOR_MEASURED` block is present but in a shape this script cannot read. A
     defect, deliberately not degraded into a pass.

A red run of this tier is NEVER a re-bless reason ($(SLOW_REBLESS_DOC)); the fast tier's golden is
not even consulted here.
"""
function slow_main(args = ARGS)
    out = slow_out_dir(args)
    mkpath(out)
    threads = Threads.nthreads()
    slow_print_header(stdout, out, threads)

    net = slow_net()
    println("NET RESOLUTION")
    println("  resolved_via       = ", net.resolved_via,
            "  (recorded, never asserted — it differs between a dev checkout and a runner)")
    println("  artifact_hash      = ", net.artifact_hash)
    println("  net_content_digest = ", net.net_content_digest)

    # The F5 binding invariant, COMPUTED AND REPORTED rather than enforced here. `run_gate` refuses
    # outright on `:provenance_mismatch`; this tier reports it loudly instead, because refusing
    # would turn the one job that can detect miscalibration into a no-op on exactly the nets whose
    # training provenance predates F6. Whoever reads this report needs the flag, not a blank page.
    prov = assert_imsize_provenance(net.m)
    println("  imsize provenance  = ", prov.status, "  (recorded = ",
            hasproperty(prov.training, :recorded) ? prov.training.recorded : "unknown",
            ", meta read from ", net.meta_source, ")")
    if net.meta_source !== :npe_jld2
        println("    NOTE: the persisted training metadata could not be read on this run, so the")
        println("    provenance status above is the ABSENCE of a check, not a failed one.")
    end
    if !prov.ok
        println("    WARNING: the gate's simulate distribution is NOT verified equal to the net's")
        println("    TRAINING image-size distribution. SBC rank uniformity is a theorem about the")
        println("    joint the estimator was trained under, so the numbers below are conditional on")
        println("    that mixture being right. Reported, not silently dropped.")
    end
    println()

    t0 = time()
    g  = sbc_gate(net.m; G = SLOW_G,
                  M = SLOW_M, L = SBC_L, bins = SBC_BINS,
                  chi2_bins = _gate_chi2_bins(SBC_BINS),
                  imsize = SBC_IMSIZE,                    # the pre-registered F5 mixture
                  imsize_set = SBC_IMSIZE_SET, imsize_weights = SBC_IMSIZE_WEIGHTS,
                  sim = default_simulator(), rng = prod_rng(SLOW_G))
    elapsed = time() - t0
    println("SBC arm complete in ", round(elapsed / 60; digits = 2), " min (",
            round(elapsed / SLOW_M; digits = 3), " s per draw at ", threads, " threads)")
    println("  realised image sizes: ", g.realised_imsize_counts)
    println()

    slow_print_columns(stdout, g)
    slow_print_coverage(stdout, g)

    # --- the break comparison ---------------------------------------------------------------
    anchors    = Any[slow_anchor_ece(c) for c in P15_BREAK_COLUMNS]
    unreadable = any(a -> a === :unreadable, anchors)
    have       = !unreadable && all(a -> a !== nothing, anchors)
    broken     = String[]
    rows       = NamedTuple[]

    println("BREAK CRITERION (D-04a): ECE(column) vs anchor + P15_ECE_MARGIN = ", P15_ECE_MARGIN)
    if unreadable
        println("  The :P15_ECE_ANCHOR_MEASURED block is PRESENT but this script cannot read its")
        println("  shape. That is a defect in the wiring between 15-07 and this file, not a pass.")
    elseif !have
        println("  NO ANCHOR MEASURED YET. The :P15_ECE_ANCHOR_MEASURED Tier-2 block has not been")
        println("  appended to test/gate/p15_consts.jl, so there is nothing to compare against.")
        println("  This run is a REPORT, not a verdict — stated rather than passed silently.")
        for (k, c) in enumerate(P15_BREAK_COLUMNS)
            push!(rows, (column = c, label = String(SBC_PARAM_LABELS[c]),
                         ece = g.per_param[c].ece, anchor = nothing,
                         threshold = nothing, delta = nothing, broke = false))
        end
    else
        for (k, c) in enumerate(P15_BREAK_COLUMNS)
            a   = anchors[k]::Float64
            thr = p15_break_threshold(a)
            e   = g.per_param[c].ece
            b   = e > thr
            b && push!(broken, "$(SBC_PARAM_LABELS[c]) (column $c): ECE $(e) > threshold $(thr)")
            println("  column ", c, " (", rpad(SBC_PARAM_LABELS[c], 8), ")  anchor ", _f(a),
                    "  threshold ", _f(thr), "  measured ", _f(e),
                    "  delta ", _f(e - a), b ? "  BREAK" : "  within")
            push!(rows, (column = c, label = String(SBC_PARAM_LABELS[c]), ece = e, anchor = a,
                         threshold = thr, delta = e - a, broke = b))
        end
    end
    println()

    status = unreadable ? :anchor_unreadable :
             !have      ? :reported_no_anchor :
             isempty(broken) ? :within_envelope : :broke

    report = (kind = :slow_calibration, schema = 1, status = status,
              grid = SLOW_G, M = SLOW_M, L = SBC_L,
              m_source = P15_SLOW_M_REQUESTED === nothing ? :pre_registered : :dispatch_input,
              seed = prod_seed(SLOW_G), global_seed = gate_global_seed(SLOW_G),
              julia_version = string(VERSION), threads = threads,
              started_utc = _slow_now(), elapsed_seconds = elapsed,
              artifact_hash = net.artifact_hash, net_content_digest = net.net_content_digest,
              resolved_via = net.resolved_via, meta_source = net.meta_source,
              imsize_provenance_status = prov.status,
              imsize_provenance_training = prov.training,
              break_columns = P15_BREAK_COLUMNS, ece_margin = P15_ECE_MARGIN,
              break_rows = rows, broke = broken,
              sbc = g,
              honesty = P15_GOLDEN_HONESTY)

    path = write_gate_report(joinpath(out, SLOW_REPORT_NAME), report)
    println("REPORT WRITTEN: ", path)
    println("  This path is the ALIVENESS EVIDENCE for the slow tier — see the header of")
    println("  .github/workflows/calibration-slow.yml. A dispatch-only workflow nobody dispatches")
    println("  is silently dead.")
    println()

    if unreadable
        println("EXIT 4 — anchor present but unreadable. Reconcile this script's `slow_anchor_ece`")
        println("with the shape 15-07 appended under :P15_ECE_ANCHOR_MEASURED.")
        return 4
    elseif !have
        println("EXIT 0 — REPORTED, NOT GATED. The slow tier is running before the anchor was")
        println("measured. Nothing here says the net is calibrated; nothing here says it is not.")
        return 0
    elseif isempty(broken)
        println("EXIT 0 — every gating column is within its break threshold at M = ", SLOW_M, ".")
        return 0
    else
        println("EXIT 1 — MISCALIBRATION DETECTED on ", length(broken), " gating column(s):")
        for b in broken
            println("  * ", b)
        end
        println()
        println("This is the tier that can say that, and this is it saying it. Do NOT re-bless the")
        println("fast-tier golden in response — it is not involved. Read ", SLOW_REBLESS_DOC, ".")
        return 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(slow_main(ARGS))
end
