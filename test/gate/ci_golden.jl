#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/ci_golden.jl --- the D-09 FAST-TIER calibration golden.
#
#   julia --project=. --threads=auto test/gate/ci_golden.jl
#   julia --project=. --threads=auto test/gate/ci_golden.jl --rebless --reason "..."
#
# ============================ WHAT THIS GATE DETECTS, AND WHAT IT DOES NOT ====================
# IT DETECTS **CHANGE**, NOT MISCALIBRATION. It runs the SBC fixture at `P15_GOLDEN_M = 8` draws.
# At M = 8 the ECE null mean is ~0.336/sqrt(8) ~= 0.12 — an order of magnitude ABOVE any
# miscalibration worth detecting — so the ECE number produced here carries NO calibration
# information whatsoever. What it does carry is a bit-level fingerprint of the whole inference
# path: the resolved net, the summary chain, the standardizers, the flow forward pass and the
# draw order. A regression in any of them moves the rank digest.
#
# MISCALIBRATION IS DETECTED BY THE **SLOW TIER**, at reported scale (`P15_SWEEP_M`). Any claim
# that this per-push job proves calibration oversells SC3. `P15_GOLDEN_HONESTY` in
# `p15_consts.jl` is the frozen wording of the same statement; keep the two in agreement.
#
# ================================ WHY IT IS DETERMINISTIC AT ALL ==============================
# BOTH halves of the gate's randomness are pinned, and both are load-bearing:
#   (1) the counter-based Philox stream `prod_rng(8)` drives the prior draw and the forward
#       simulation — the property that let Phase 11's lost 54 MB pool be rebuilt byte-identically
#       at zero iteration cost;
#   (2) `seed_gate_global!(8)` (called inside `sbc_ranks_and_spread`) pins Julia's GLOBAL RNG,
#       which is where `sampleposterior` draws the flow's base samples from — it threads no `rng`
#       (see the GLOBAL-RNG SEEDING block in `harness.jl`).
# `P15_GOLDEN_IMSIZE` is a CONCRETE `(256, 256)`, not the `:mixture` sentinel, so `gate_imsize`
# takes its scalar branch and consumes ZERO numbers from the rng; the draw sequence therefore
# cannot shift if the mixture weights ever change.
#
# THE HONEST CAVEAT: `seed_gate_global!` ultimately calls `Random.seed!` on Julia's default
# task-local RNG, whose stream semantics are NOT a cross-version contract and have changed across
# Julia minors. A bit-exact golden is therefore a **JULIA-PATCH** assertion as much as a code
# assertion. The workflow pins `version: '1.12.6'`; a Julia patch bump is a LEGITIMATE re-bless
# reason and must be named as such. Float-summation reassociation is absorbed by
# `P15_GOLDEN_ECE_TOL`; the integer rank digest has no tolerance by design.
#
# ==================================== RE-BLESSING ============================================
# Governed by `.planning/phases/15-calibration-operating-envelope-and-ci-gate/15-REBLESS.md`.
# `--rebless` REFUSES without a non-empty reason, APPENDS a new sentinel-keyed Tier-2 block to
# `p15_consts.jl` (never edits or deletes an earlier one), and prints old and new side by side for
# the commit message. A golden that can be silently regenerated is a rubber stamp.
#
# DELIBERATELY NOT ROUTED THROUGH `run_gate.jl`: under `SBC_REQUIRE_IMSIZE_PROVENANCE = true` it
# asserts the gate mixture against the net's training provenance, and this golden uses one
# concrete size on purpose. It is also not `run_envelope`: the golden is not a reported run.

using ProteinCoLoc
import Artifacts

# --- Ordered, guarded includes: pre-registration -> harness -> sbc ---------------------------
# Guarded so this file is both standalone-runnable and includable from a module that has already
# loaded the Phase-15 pre-registration (`test_p15_golden.jl` does exactly that, in an ISOLATED
# module — `p15_consts.jl` defines the same const names as every other gate pre-registration).
if !isdefined(@__MODULE__, :SBC_M)
    include(joinpath(@__DIR__, "p15_consts.jl"))
end
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "harness.jl"))
isdefined(@__MODULE__, :sbc_gate)            || include(joinpath(@__DIR__, "sbc.jl"))

# The golden's bars (`P15_GOLDEN_*`) exist ONLY in the Phase-15 pre-registration. If some OTHER
# consts file was loaded into this module first, the guard above would have SKIPPED the include
# and every bar below would silently be the wrong one — fail loudly instead.
if !(isdefined(@__MODULE__, :GATE_CONSTS_VERSION) && GATE_CONSTS_VERSION == 15 &&
     isdefined(@__MODULE__, :P15_GOLDEN_M))
    error("ci_golden.jl: the loaded pre-registration is not the Phase-15 one " *
          "(GATE_CONSTS_VERSION = " *
          (isdefined(@__MODULE__, :GATE_CONSTS_VERSION) ? string(GATE_CONSTS_VERSION) : "undefined") *
          "). Load `p15_consts.jl` — into an isolated module if another consts file is already " *
          "present — before including this file.")
end

"The pre-registration file the golden blocks are appended to. Resolved from THIS file's directory."
const P15_GOLDEN_CONSTS_PATH = joinpath(@__DIR__, "p15_consts.jl")

"The re-bless procedure. `--rebless` prints this path when it refuses."
const P15_REBLESS_DOC =
    ".planning/phases/15-calibration-operating-envelope-and-ci-gate/15-REBLESS.md"

"""
The ONLY reasons that legitimise a re-bless (D-09). Printed verbatim when `--rebless` refuses,
so the refusal teaches the rule rather than merely enforcing it.
"""
const P15_REBLESS_LEGITIMATE_REASONS = (
    "a Julia PATCH-VERSION bump — name the old and the new version in the reason",
    "a deliberate, reviewed change to the inference path (net, summary chain, standardizers, flow)",
    "a re-shipped net (a new grid_8 artifact hash)",
)

"The one reason that is NEVER legitimate. Quoted by the refusal message and by 15-REBLESS.md."
const P15_REBLESS_ILLEGITIMATE_REASON =
    "a red CI run nobody has explained. The gate going red IS the gate working. Re-blessing to " *
    "make it green is the \"tune until it passes\" move this project has already refused four times."

"Upper bound on the Tier-2 golden block search. A sequence longer than this is a smell, not a need."
const P15_GOLDEN_MAX_BLOCKS = 99

"Block 1 owns the bare sentinel `:P15_GOLDEN_RANK_DIGEST`; block n > 1 owns `..._Vn`."
_p15_golden_suffix(n::Integer) = n == 1 ? "" : "_V$(Int(n))"

# =============================================================================================
# MEASUREMENT
# =============================================================================================

"""
    p15_golden_net() -> NamedTuple

Resolve the shipped 8×8 net and its content identity.

`ProteinCoLoc.estimator_for(8)` is used rather than a bare path because it resolves in three
INTEGRITY-CHECKED steps (`src/registry.jl:126-161`): the content-addressed artifact store, then
the in-repo `artifacts/amended_v2/grid_8` dev dir WITH an explicit tree-sha1 verification against
`Artifacts.toml`, then a lazy download of the release asset. So it works on this machine, on a dev
checkout and on a runner with no in-repo `artifacts/` tree, and it never loads an unverified model.

THE FIELD MISMATCH THE SHIM EXISTS FOR — do NOT "simplify" it away. `estimator_for` returns an
`EstimatorBundle` STRUCT with fields `(grid, npe, ratio, ood_nulls, zt, θzt, calibration)`
(`src/registry.jl:59-67`). The gate harness's `m` is the NamedTuple
`(estimator, θzt, zt, arch, meta)` that `ProteinCoLoc.load_estimator` returns
(`src/amortized/persist.jl:101-102`), and the harness reads `m.estimator` / `m.zt` / `m.θzt`
(`harness.jl:197-242`). The posterior net is `b.npe`, NOT `b.estimator` — there is no such field.
"""
function p15_golden_net()
    b = ProteinCoLoc.estimator_for(8)          # three-path, hash-verified resolution
    m = (estimator = b.npe, θzt = b.θzt, zt = b.zt, arch = nothing, meta = nothing)

    toml     = ProteinCoLoc._artifacts_toml()
    expected = Artifacts.artifact_hash("grid_8", toml)
    expected === nothing && error("ci_golden: Artifacts.toml has no `grid_8` entry at $toml.")

    # Which of the three paths actually served the bytes. RECORDED but NEVER asserted: it
    # legitimately differs between this machine (store / dev dir) and a runner (download).
    devdir = ProteinCoLoc._dev_bundle_dir(8)
    dir, via = if Artifacts.artifact_exists(expected)
        (Artifacts.artifact_path(expected), :artifact_store)
    elseif isdir(devdir)
        (devdir, :dev_dir)
    else
        (nothing, :unresolved)
    end

    # An INDEPENDENT digest of the bundle's BYTES, so a silent net substitution fails HERE rather
    # than surfacing downstream as unexplained numeric drift. `artifact_hash` alone cannot do that
    # job: it is read out of `Artifacts.toml`, so it detects a changed PIN, not changed bytes.
    #
    # WHY NOT A GIT TREE HASH (measured, not assumed). `Pkg.GitTools.tree_hash` is FILE-MODE
    # SENSITIVE, and the two resolution paths hand back different modes for byte-identical files:
    # the artifact store's copy is read-only (0o444) and the in-repo dev copy is writable (0o666),
    # which on this checkout yields tree hashes `7d72b47e…` and `90e6b63a…` for four files that
    # compare byte-EQUAL. Asserting a tree hash would therefore fail on the resolution path alone
    # — precisely the spurious red D-09 warns gets a gate disabled.
    #
    # WHY `hash` AND NOT SHA-256. The `SHA` stdlib is not a `[dep]`, and importing it aborts the
    # whole suite under `Pkg.test()`'s sandbox (15-01, deviation 2). `Base.hash` needs no import,
    # is content-based and mode-independent, and is the SAME primitive the rank digest already
    # stakes determinism on under the pinned Julia. STATED PLAINLY: this is a non-cryptographic
    # digest. It catches accident and regression, not an adversary — adversarial substitution is
    # what the store's content addressing and `_verify_tree_sha1` cover.
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

    return (m = m, bundle = b, dir = dir, resolved_via = via,
            artifact_hash = string(expected), net_content_digest = digest)
end

"""
    p15_golden_measure() -> NamedTuple

Run the fixture-scale SBC ONCE and return everything the golden asserts plus everything a failure
report needs: `(artifact_hash, net_content_digest, resolved_via, rank_digest, ece, labels, ranks,
global_seed, julia_version)`.

The MCE field of `sbc_gate`'s per-column results is read NOWHERE in this file: it is pinned at
exactly 0.99 on all eight columns by empty reliability bins, so it is an artifact, not a
measurement (`P15_MCE_FORBIDDEN`).
"""
function p15_golden_measure()
    net = p15_golden_net()
    g = sbc_gate(net.m; G = 8,
                 M = P15_GOLDEN_M, L = P15_GOLDEN_L, bins = P15_GOLDEN_BINS,
                 chi2_bins = _gate_chi2_bins(P15_GOLDEN_BINS),
                 imsize = P15_GOLDEN_IMSIZE,        # concrete: the scalar branch, zero rng draws
                 imsize_set = nothing, imsize_weights = nothing,
                 sim = default_simulator(), rng = prod_rng(8))
    return (artifact_hash      = net.artifact_hash,
            net_content_digest = net.net_content_digest,
            resolved_via  = net.resolved_via,
            rank_digest   = string(hash(g.ranks)),
            ece           = Float64[p.ece for p in g.per_param],
            labels        = String[string(p.label) for p in g.per_param],
            ranks         = g.ranks,
            global_seed   = gate_global_seed(8),
            julia_version = string(VERSION))
end

# =============================================================================================
# THE COMMITTED GOLDEN (Tier-2 blocks in p15_consts.jl)
# =============================================================================================

"""
    p15_golden_blessed(mod = @__MODULE__) -> Union{NamedTuple, Nothing}

The HIGHEST-numbered Tier-2 golden block currently defined, or `nothing` if none is. Reading the
highest is what makes `--rebless` an APPEND: a new block supersedes without any line above it
being rewritten, and the git history of `p15_consts.jl` stays the audit trail.
"""
function p15_golden_blessed(mod::Module = @__MODULE__)
    for n in P15_GOLDEN_MAX_BLOCKS:-1:1
        s   = _p15_golden_suffix(n)
        sym = Symbol("P15_GOLDEN_RANK_DIGEST", s)
        isdefined(mod, sym) || continue
        fld(k) = getfield(mod, Symbol("P15_GOLDEN_", k, s))
        has(k) = isdefined(mod, Symbol("P15_GOLDEN_", k, s))
        return (block         = n,
                sentinel      = String(sym),
                rank_digest   = fld("RANK_DIGEST"),
                artifact_hash = fld("ARTIFACT_HASH"),
                net_content_digest = has("NET_CONTENT_DIGEST") ? fld("NET_CONTENT_DIGEST") : nothing,
                ece           = fld("ECE"),
                julia_version = fld("JULIA_VERSION"),
                blessed_on    = fld("BLESSED_ON"),
                reason        = has("REASON") ? fld("REASON") : "initial blessing, plan 15-06")
    end
    return nothing
end

"The block number a re-bless would open next (1 when nothing is blessed yet)."
function p15_golden_next_block(mod::Module = @__MODULE__)
    g = p15_golden_blessed(mod)
    return g === nothing ? 1 : g.block + 1
end

# =============================================================================================
# THE THREE ASSERTIONS — cheapest first, ALL failures collected
# =============================================================================================

"""
    p15_golden_check(measured, golden) -> NamedTuple

Compare a `p15_golden_measure()` result against a committed golden block. Returns
`(ok, failures)`, where `failures` is a vector of human-readable strings — ALL of them, never just
the first, because "the net changed AND the ranks moved" and "only the ranks moved" are different
diagnoses and stopping early hides which one you have.

Order is cheapest-and-most-diagnostic first:
  1. `artifact_hash` — a NET SWAP. Exact match. Not a tolerance question and never a re-bless
     without a stated reason.
  2. `net_content_digest` — the SAME identity recomputed FROM THE BYTES, so a substitution that
     leaves `Artifacts.toml` untouched still fails here and not as downstream numeric drift.
  3. `rank_digest`  — exact, zero tolerance. Catches any change in draw order, the summary chain,
     the standardizers or the flow forward pass.
  4. `ece[1:8]`     — within `P15_GOLDEN_ECE_TOL` (float reassociation headroom only).
"""
function p15_golden_check(measured, golden)
    failures = String[]

    if measured.artifact_hash != golden.artifact_hash
        push!(failures,
            "ARTIFACT HASH CHANGED — the resolved grid_8 bundle is NOT the one this golden was " *
            "blessed against.\n    golden   = $(golden.artifact_hash)\n    measured = " *
            "$(measured.artifact_hash)\n    A net swap is a RE-BLESS-WITH-REASON event, not a " *
            "tolerance question. See $(P15_REBLESS_DOC).")
    end

    if golden.net_content_digest !== nothing && !isempty(golden.net_content_digest) &&
       measured.net_content_digest != golden.net_content_digest
        push!(failures,
            "NET CONTENT DIGEST CHANGED — the resolved bundle's BYTES differ from the blessed " *
            "ones even though Artifacts.toml may be untouched.\n" *
            "    golden   = $(golden.net_content_digest)\n" *
            "    measured = $(measured.net_content_digest)\n" *
            "    A net swap is a RE-BLESS-WITH-REASON event. See $(P15_REBLESS_DOC).")
    end

    if measured.rank_digest != golden.rank_digest
        push!(failures,
            "RANK DIGEST CHANGED — the SBC rank table is not bit-identical.\n" *
            "    golden   = $(golden.rank_digest)\n    measured = $(measured.rank_digest)\n" *
            "    This has NO tolerance by design. Either the inference path changed, or the Julia " *
            "patch version did (`Random.seed!` stream semantics are not a cross-version contract; " *
            "golden was blessed on Julia $(golden.julia_version), this run is " *
            "$(measured.julia_version)).")
    end

    n = min(length(measured.ece), length(golden.ece))
    for i in 1:n
        d = abs(measured.ece[i] - golden.ece[i])
        d <= P15_GOLDEN_ECE_TOL && continue
        lab = i <= length(measured.labels) ? measured.labels[i] : "column $i"
        push!(failures,
            "ECE[$i] ($lab) moved by $(d), above P15_GOLDEN_ECE_TOL = $(P15_GOLDEN_ECE_TOL).\n" *
            "    golden   = $(golden.ece[i])\n    measured = $(measured.ece[i])")
    end
    length(measured.ece) == length(golden.ece) ||
        push!(failures, "ECE vector length changed: golden $(length(golden.ece)), " *
                        "measured $(length(measured.ece)).")

    return (ok = isempty(failures), failures = failures)
end

"""
    p15_golden_table(measured, golden; io = stdout)

Print the OLD and NEW values side by side. Used on any failure and by `--rebless`, in a form that
pastes into a commit message.
"""
function p15_golden_table(measured, golden; io::IO = stdout)
    println(io, "  field                     golden                          measured")
    println(io, "  ------------------------  ------------------------------  ------------------------------")
    row(name, a, b) = println(io, "  ", rpad(name, 24), "  ", rpad(string(a), 30), "  ", string(b))
    row("artifact_hash", golden === nothing ? "-" : golden.artifact_hash, measured.artifact_hash)
    row("net_content_digest", golden === nothing ? "-" :
        (golden.net_content_digest === nothing ? "(not blessed)" : golden.net_content_digest),
        measured.net_content_digest)
    row("rank_digest", golden === nothing ? "-" : golden.rank_digest, measured.rank_digest)
    row("julia_version", golden === nothing ? "-" : golden.julia_version, measured.julia_version)
    for i in eachindex(measured.ece)
        lab = i <= length(measured.labels) ? measured.labels[i] : "col$i"
        gv  = golden === nothing || i > length(golden.ece) ? "-" : golden.ece[i]
        row("ece[$i] $lab", gv, measured.ece[i])
    end
    println(io, "  resolved_via = ", measured.resolved_via,
                "  (recorded, NOT asserted — it differs by machine)")
    println(io, "  global_seed  = ", measured.global_seed)
    return nothing
end

# =============================================================================================
# BLESSING — always an APPEND, never an edit
# =============================================================================================

"Today's date as `YYYY-MM-DD`, without the `Dates` stdlib (not a [dep]; see `p15_golden_net`)."
_p15_today() = Base.Libc.strftime("%Y-%m-%d", time())

"""
    p15_golden_block_text(measured; block, reason) -> String

The paste-ready Tier-2 block for `measured`, guarded on its OWN sentinel. Block 1 owns the bare
`:P15_GOLDEN_RANK_DIGEST` sentinel reserved by `p15_consts.jl`; every later block owns `..._Vn`,
because a second append reusing an earlier block's sentinel would be SILENTLY SKIPPED once that
block exists, and an append with no guard at all is a hard `const` redefinition on re-include.
"""
function p15_golden_block_text(measured; block::Integer = 1, reason::AbstractString = "")
    s   = _p15_golden_suffix(block)
    ece = join(("        " * repr(x) for x in measured.ece), ",\n")
    r   = replace(strip(String(reason)), "\n" => " ")
    hdr = block == 1 ?
        "# TIER 2 / block 1 --- the D-09 fast-gate golden. INITIAL BLESSING, plan 15-06." :
        "# TIER 2 / block $(block) --- the D-09 fast-gate golden, RE-BLESSED."
    why = block == 1 ?
        "# REASON: initial blessing, plan 15-06 (no golden existed; the fast tier had nothing to\n" *
        "#         assert against)." :
        "# REASON (required by `--rebless`, D-09): $(r)"
    return string(
        "\n",
        "# =========================================================================================\n",
        hdr, "\n",
        "# =========================================================================================\n",
        block == 1 ?
            "# APPENDED at the end of the file. Nothing above it was touched (D-09).\n" :
            "# APPENDED, never edited: block $(block - 1) above is left BYTE-IDENTICAL on purpose,\n" *
            "# because the git history of this file IS the pre-registration audit trail (D-09).\n",
        why, "\n",
        "# Measured by `julia --project=. --threads=auto test/gate/ci_golden.jl", block == 1 ? "" : " --rebless", "`\n",
        "# on Julia $(measured.julia_version), against the grid_8 bundle resolved via ",
        string(measured.resolved_via), ".\n",
        "# WHAT IT ASSERTS: change, NOT calibration — see P15_GOLDEN_HONESTY and the header of\n",
        "# `ci_golden.jl`. Re-blessing is governed by $(P15_REBLESS_DOC).\n",
        "if !isdefined(@__MODULE__, :P15_GOLDEN_RANK_DIGEST", s, ")\n",
        "    # `string(hash(ranks))` over the M×8 integer rank table. NO tolerance, by design.\n",
        "    const P15_GOLDEN_RANK_DIGEST", s, "   = \"", measured.rank_digest, "\"\n",
        "    # The `grid_8` git-tree-sha1 pinned in Artifacts.toml — a net swap fails here FIRST.\n",
        "    const P15_GOLDEN_ARTIFACT_HASH", s, " = \"", measured.artifact_hash, "\"\n",
        "    # The same identity recomputed FROM THE BYTES — mode-independent, so it does not\n",
        "    # depend on which of the three resolution paths served them (see `p15_golden_net`).\n",
        "    const P15_GOLDEN_NET_CONTENT_DIGEST", s, " = \"", measured.net_content_digest, "\"\n",
        "    # The eight per-column ECEs, asserted within P15_GOLDEN_ECE_TOL. At M = 8 these carry\n",
        "    # NO calibration information; they are a numeric fingerprint. Order = SBC_PARAM_LABELS.\n",
        "    const P15_GOLDEN_ECE", s, "           = [\n", ece, "]\n",
        "    # A bit-exact golden is a JULIA-PATCH assertion too: `Random.seed!` stream semantics\n",
        "    # are not a cross-version contract. A patch bump is a LEGITIMATE re-bless reason.\n",
        "    const P15_GOLDEN_JULIA_VERSION", s, " = \"", measured.julia_version, "\"\n",
        "    const P15_GOLDEN_BLESSED_ON", s, "    = \"", _p15_today(), "\"\n",
        block == 1 ? "" :
            "    const P15_GOLDEN_REASON$(s)         = \"$(replace(r, "\"" => "'"))\"\n",
        "end\n")
end

"""
    p15_golden_append!(path, text) -> String

APPEND `text` to `path` and return it. Nothing above the appended bytes is read, rewritten or
truncated — the file is opened in append mode, which is what makes "appends, never edits" a
property of the mechanism rather than a promise about the caller.

Line endings are matched to whatever the file already uses, so a Windows working tree
(`.gitattributes` `* text=auto` + `core.autocrlf=true`) does not end up with a mixed file.
"""
function p15_golden_append!(path::AbstractString, text::AbstractString)
    crlf = occursin("\r\n", read(path, String))
    out  = crlf ? replace(text, "\r\n" => "\n", "\n" => "\r\n") : text
    open(path, "a") do io
        write(io, out)
    end
    return path
end

# =============================================================================================
# ENTRY POINT
# =============================================================================================

"""
    p15_rebless_reason(args) -> Union{String, Nothing}

The re-bless reason from `--reason "..."` or the `P15_REBLESS_REASON` environment variable, or
`nothing` when neither carries a non-empty string. The env var exists so the reason can be
supplied by a script without shell-quoting games; it is NOT a way to skip stating one.
"""
function p15_rebless_reason(args = ARGS)
    for i in 1:(length(args) - 1)
        args[i] == "--reason" && !isempty(strip(args[i + 1])) && return String(strip(args[i + 1]))
    end
    env = get(ENV, "P15_REBLESS_REASON", "")
    return isempty(strip(env)) ? nothing : String(strip(env))
end

"Print the refusal that teaches the rule (D-09). Touches no file."
function p15_rebless_refusal(io::IO = stdout)
    println(io, "REFUSED: --rebless requires a non-empty reason.")
    println(io, "  Supply it as `--reason \"...\"` or in the P15_REBLESS_REASON environment variable.")
    println(io, "  NO FILE WAS TOUCHED.")
    println(io)
    println(io, "  A re-bless is legitimate ONLY for:")
    for r in P15_REBLESS_LEGITIMATE_REASONS
        println(io, "    - ", r)
    end
    println(io)
    println(io, "  It is NEVER legitimate for:")
    println(io, "    - ", P15_REBLESS_ILLEGITIMATE_REASON)
    println(io)
    println(io, "  Procedure: ", P15_REBLESS_DOC)
    return nothing
end

"""
    p15_golden_main(args = ARGS) -> Int

The script entry point. Exit codes, all of them meaningful:

  0  the measured values match the highest committed golden block
  1  a mismatch (the side-by-side table is printed first)
  2  NO golden block exists yet — the measured values are printed as a paste-ready Tier-2 block
     and nothing is written. Blessing is an explicit commit, never a side effect of a run.
  3  `--rebless` was passed without a reason. No file was touched.
"""
function p15_golden_main(args = ARGS)
    rebless = "--rebless" in args

    if rebless
        reason = p15_rebless_reason(args)
        if reason === nothing
            p15_rebless_refusal()
            return 3
        end
        measured = p15_golden_measure()
        old      = p15_golden_blessed()
        block    = p15_golden_next_block()
        println("RE-BLESSING the D-09 fast-gate golden — block $(block).")
        println("Reason: ", reason)
        println()
        p15_golden_table(measured, old)
        println()
        p15_golden_append!(P15_GOLDEN_CONSTS_PATH,
                           p15_golden_block_text(measured; block = block, reason = reason))
        println("Appended block $(block) to ", P15_GOLDEN_CONSTS_PATH)
        println("Nothing above it was modified. Commit the append ALONE, with the reason and the")
        println("table above in the commit message (", P15_REBLESS_DOC, ").")
        return 0
    end

    measured = p15_golden_measure()
    golden   = p15_golden_blessed()

    if golden === nothing
        println("NO GOLDEN BLESSED YET. Nothing was written — a golden is blessed by an explicit")
        println("commit, never as a side effect of running this script.")
        println()
        println("Paste the block below at the END of test/gate/p15_consts.jl and commit it:")
        println(p15_golden_block_text(measured; block = 1))
        return 2
    end

    r = p15_golden_check(measured, golden)
    if r.ok
        println("D-09 fast-gate golden OK (block $(golden.block), blessed $(golden.blessed_on) on ",
                "Julia $(golden.julia_version)).")
        println("  rank_digest   = ", measured.rank_digest)
        println("  artifact_hash = ", measured.artifact_hash)
        println("  NOTE: this gate detects CHANGE, not miscalibration (M = $(P15_GOLDEN_M)).")
        return 0
    end

    println("D-09 FAST-GATE GOLDEN FAILED — $(length(r.failures)) mismatch(es) against block ",
            golden.block, ".")
    println()
    p15_golden_table(measured, golden)
    println()
    for f in r.failures
        println("  * ", f)
    end
    println()
    println("A red golden is the gate WORKING. Do NOT re-bless to make it green.")
    println("Read ", P15_REBLESS_DOC, " before doing anything else.")
    return 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(p15_golden_main(ARGS))
end
