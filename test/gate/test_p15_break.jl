#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/test_p15_break.jl --- fixtures for every Phase-15 DECISION RULE, plus two source
# guards, written BEFORE any rung exists to argue with them.
#
# `test/gate/p15_envelope.jl` is the executable form of the break criterion (D-04/D-04a/D-05) and
# the SC2 verdict (D-06/D-07/D-13/D-14). Those are the two most easily-corrupted decisions in this
# phase: the anchor is arithmetic that looks like a preference, and "the flag never fired" has two
# opposite meanings. This file pins each of them to a fixture now, while no result exists that a
# reader could want them to reach.
#
# NO NET, NO FORWARD MODEL, NO INFERENCE. Every number below is either arithmetic on the test
# design or a hand-built record. The file runs in seconds.
#
# ISOLATED MODULE, NOT STYLE. `p15_envelope.jl` loads `p15_consts.jl`, which defines the same
# const names (`SBC_M`, `prod_seed`, ...) that `run_gate.jl` has already loaded into the test
# module through the template, and every gate consts file is guarded on `:SBC_M`. A direct include
# would be silently skipped and every assertion below would read the WRONG pre-registration.
#
# RUNS BOTH WAYS: standalone (`julia --project=. --threads=auto test/gate/test_p15_break.jl`) and
# `include`d from `test/runtests.jl`. Every path is built from `@__DIR__`, never from the process
# working directory.

using Test
using ProteinCoLoc
import Random

module P15BreakEnv
    include(joinpath(@__DIR__, "p15_envelope.jl"))
end

const P15BRK_GATE_DIR = @__DIR__

# =============================================================================================
# The two source-guard predicates, and why their tokens are ASSEMBLED rather than written
# =============================================================================================
#
# A guard that greps for a token cannot contain that token, or it flags itself and every run
# after it is a lie in one direction or the other. Both banned tokens are therefore built from
# fragments at load time, so the literal never appears anywhere in this file — including in the
# failure messages and in the deliberately-violating fixture of the last testset.

const _P15BRK_FIELD_TOKEN     = "m" * "ce"
const _P15BRK_QUANTILE_TOKEN  = "OOD_ID_" * "QUANTILE"
const _P15BRK_THRESHOLD_TOKEN = "id_" * "threshold"

const _P15BRK_FIELD_RX  = Regex("\\b" * _P15BRK_FIELD_TOKEN * "\\b")
const _P15BRK_RETUNE_RX = Regex("(" * _P15BRK_QUANTILE_TOKEN * "|" *
                                      _P15BRK_THRESHOLD_TOKEN * ")\\s*=")

"""
    _p15brk_code_of(line) -> String

The CODE part of one source line: everything up to the first `#` that is not inside a double-
quoted string.

CONVENTIONS C-05: a STRUCTURAL BAN is checked on COMMENT-STRIPPED source, because "a ban a
comment can satisfy is not a ban" — and, symmetrically, a ban a comment can VIOLATE is not a ban
either. Every file scanned below legitimately DISCUSSES what it forbids; that discussion is the
point, so the check must not see it. Note this must strip TRAILING comments and not only
whole-line ones: `p15_consts.jl:698` declares the prohibition itself with the field name in a
trailing comment on a `const` line, and a whole-line-only filter would flag the very file that
states the rule.

Limitation, stated rather than glossed: the string tracking is per line, so a `#` inside a
multi-line docstring cuts the line. That errs toward stripping prose, which is the correct
direction for a structural ban.
"""
function _p15brk_code_of(line::AbstractString)
    out = IOBuffer()
    in_string = false
    escaped   = false
    for ch in line
        if in_string
            if escaped
                escaped = false
            elseif ch == '\\'
                escaped = true
            elseif ch == '"'
                in_string = false
            end
            print(out, ch)
        else
            ch == '#' && break
            ch == '"' && (in_string = true)
            print(out, ch)
        end
    end
    return String(take!(out))
end

"Every (line number, raw line) in `src` whose CODE part matches `rx`."
function _p15brk_scan(src::AbstractString, rx::Regex)
    hits = Tuple{Int,String}[]
    for (i, raw) in enumerate(split(src, '\n'))
        occursin(rx, _p15brk_code_of(raw)) && push!(hits, (i, strip(String(raw))))
    end
    return hits
end

_p15brk_scan_file(path::AbstractString, rx::Regex) = _p15brk_scan(read(path, String), rx)

"Read a constant whose NAME is assembled at run time, so the literal never appears in this file."
_p15brk_quantile_of(m::Module) = getproperty(m, Symbol(_P15BRK_QUANTILE_TOKEN))

_p15brk_path(name::AbstractString) = joinpath(P15BRK_GATE_DIR, name)

"Every `test/gate/test_p15_*.jl` that exists right now — sibling plans add theirs without edits."
const P15BRK_TEST_FILES = sort(filter(f -> startswith(f, "test_p15_") && endswith(f, ".jl"),
                                      readdir(P15BRK_GATE_DIR)))

"The files the field ban is checked over (D-04a). Absent siblings are simply not in the list."
const P15BRK_FIELD_SCAN = filter(isfile, vcat(
    [_p15brk_path("p15_consts.jl"), _p15brk_path("p15_misspec.jl"),
     _p15brk_path("p15_envelope.jl")],
    [_p15brk_path(f) for f in P15BRK_TEST_FILES]))

"""
The files the detector-retuning ban is checked over as a ZERO-HIT ban (D-14).

`p15_consts.jl` is deliberately NOT in this list, and is handled separately and MORE STRICTLY in
the same testset. It byte-forks the amended-v2 rule set, which RECORDS the shipped detector
operating point as part of the pre-registration; a zero-hit ban there would either be unsatisfiable
or would force editing a frozen file. Value equality against the shipped constant is the stronger
check anyway: absence of the line would not catch a CHANGED number, and equality does.
"""
const P15BRK_RETUNE_SCAN = filter(isfile, vcat(
    [_p15brk_path("p15_misspec.jl"), _p15brk_path("p15_envelope.jl"),
     _p15brk_path("ci_golden.jl")],
    [_p15brk_path(f) for f in P15BRK_TEST_FILES]))

# =============================================================================================

@testset "Phase-15 break criterion and SC2 verdict" begin

    # -----------------------------------------------------------------------------------------
    @testset "1 — ECE identity reproduces the production statistic (D-04)" begin
        rng = Random.Xoshiro(20260804)
        L   = P15BreakEnv.SBC_L

        uniform_draw = rand(rng, 0:L, 600)
        # OVER-DISPERSED: mass pushed to both rank extremes, i.e. posteriors far too narrow.
        overdispersed = [rand(rng, Bool) ? rand(rng, 0:49) : rand(rng, (L - 49):L)
                         for _ in 1:600]
        # DEGENERATE: every rank at the floor. Coverage is 0 at every level, so ECE = mean(alpha).
        degenerate = zeros(Int, 600)

        for (name, ranks) in (("uniform", uniform_draw),
                              ("over-dispersed", overdispersed),
                              ("degenerate", degenerate))
            ours = P15BreakEnv.p15_ece(ranks, L)
            prod = P15BreakEnv.sbc_calibration(collect(ranks); L = L,
                                               n_bins = P15BreakEnv.SBC_BINS).ece
            @test isapprox(ours, prod; atol = 1e-12)
            name == "degenerate" && @test isapprox(ours, 0.5; atol = 1e-12)
        end

        # THE PROPERTY THAT MAKES THE STATISTIC WELL DEFINED: with 19 nominal levels spaced 0.05,
        # each level lands in its own reliability bin for any bin count fine enough, the empty bins
        # carry weight 0, and the ECE therefore does not depend on `SBC_BINS` at all. (At 10 bins
        # two levels would share a bin and the invariance genuinely breaks — the claim is about
        # resolution, not about magic.)
        for ranks in (uniform_draw, overdispersed, degenerate)
            eces = [P15BreakEnv.sbc_calibration(collect(ranks); L = L, n_bins = b).ece
                    for b in (20, 50, 100)]
            @test isapprox(eces[1], eces[2]; atol = 1e-12)
            @test isapprox(eces[2], eces[3]; atol = 1e-12)
            @test isapprox(P15BreakEnv.p15_ece(ranks, L), eces[2]; atol = 1e-12)
        end
    end

    # -----------------------------------------------------------------------------------------
    @testset "2 — the ECE null is design arithmetic, not an outcome (D-04a)" begin
        M, L = P15BreakEnv.P15_SWEEP_M, P15BreakEnv.SBC_L
        @test M == 1000

        # B = 2000 here rather than the pinned P15_ECE_NULL_B = 20_000, so the testset stays in
        # seconds. The tolerances below are Monte-Carlo error at THIS B, derived, not chosen:
        # sd(ECE_0) ~ 0.0045 at M = 1000, so SE(mean) ~ 1.0e-4 and SE(q95) ~ 2.1e-4. The bars used
        # are >= 6 SE.
        null = P15BreakEnv.p15_ece_null(M, L; B = 2000)
        @test isapprox(null.mean, P15BreakEnv.P15_ECE_NULL_MEAN; atol = 0.0015)
        @test isapprox(null.q,    P15BreakEnv.P15_ECE_NULL_Q95;  atol = 0.0030)

        # THE MARGIN IS THE NULL'S OWN SPREAD, so it regenerates from the design alone.
        @test isapprox(null.q - null.mean, P15BreakEnv.P15_ECE_MARGIN; atol = 0.0015)
        @test isapprox(P15BreakEnv.P15_ECE_NULL_Q95 - P15BreakEnv.P15_ECE_NULL_MEAN,
                       P15BreakEnv.P15_ECE_MARGIN; atol = 1e-9)

        # The closed form E[ECE_0] ~= 0.336/sqrt(M) that the derivation comment in p15_consts.jl
        # cites, checked at three M rather than only at the swept one.
        for m in (250, 1000, 2000)
            mean_m = P15BreakEnv.p15_ece_null(m, L; B = 2000).mean
            closed = 0.336 / sqrt(m)
            @test abs(mean_m - closed) / closed < 0.05
        end

        # ------------------------------------------------------------------------------------
        # D-04a MADE VISIBLE. 0.0042894736842104715 is the SHIPPED Delta-rho ECE from
        # artifacts/amended_v2/grid_8/gate_report_8.jld2, measured at M = 2000. 0.00751 is the
        # null MEAN of a PERFECTLY CALIBRATED estimator at that same M. The shipped value sits
        # BELOW its own null mean, so a break threshold anchored there would fire on roughly 90 %
        # of perfectly calibrated rungs. THAT is why the anchor is rung 0 at the sweep's own M.
        @test 0.0042894736842104715 < 0.00751

        # …and the same statement measured rather than quoted, so it does not rest on two literals.
        @test 0.0042894736842104715 < P15BreakEnv.p15_ece_null(2000, L; B = 2000).mean

        # The break threshold is the anchor PLUS the margin, and the margin alone is never it.
        @test P15BreakEnv.p15_break_threshold_for(0.0) == P15BreakEnv.P15_ECE_MARGIN
        @test P15BreakEnv.p15_break_threshold_for(0.02) == 0.02 + P15BreakEnv.P15_ECE_MARGIN
    end

    # -----------------------------------------------------------------------------------------
    @testset "3 — the shipped-anchor reader is honest (D-04a)" begin
        @test P15BreakEnv.p15_read_shipped_anchor(joinpath(P15BRK_GATE_DIR,
                  "definitely_not_a_report_$(rand(Random.Xoshiro(1), 1:10^9)).jld2")) === nothing

        report = joinpath(P15BreakEnv.P15_ENVELOPE_REPO_ROOT, P15BreakEnv.P15_ARTIFACTS_ROOT,
                          "grid_8", "gate_report_8.jld2")
        if isfile(report)
            anchor = P15BreakEnv.p15_read_shipped_anchor(report)
            @test anchor !== nothing
            @test isapprox(anchor.ece["ρ_true"], 0.03713157894736839; atol = 1e-12)
            @test isapprox(anchor.ece["Δρ"],     0.0042894736842104715; atol = 1e-12)
            @test anchor.M == 2000
            # The shipped M is NOT the sweep's M — which is the whole reason it is not the anchor.
            @test anchor.M != P15BreakEnv.P15_SWEEP_M
        else
            @info "shipped gate report absent (artifacts/ is git-ignored) — the value assertions " *
                  "in testset 3 were SKIPPED, not failed; the absent-path branch still ran" report
        end
    end

    # -----------------------------------------------------------------------------------------
    @testset "4 — the break reduction reads ONLY columns 1 and 8 (D-05)" begin
        thr = 0.05
        # Columns 2..7 are CATASTROPHICALLY miscalibrated at every rung (0.5, i.e. maximal), while
        # the two gating columns stay clean at 0.01. If the nuisance marginals gated anything at
        # all, this would break at rung 1.
        mkrow(rung; c1 = 0.01, c8 = 0.01) =
            (rung = rung, ece = [c1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, c8])

        clean = [mkrow(k) for k in 1:5]
        b0 = P15BreakEnv.p15_axis_break(clean, thr)
        @test b0.L_break === nothing
        @test b0.per_column[1] === nothing
        @test b0.per_column[8] === nothing

        # Delta-rho (column 8) degrades at rung 3.
        c8only = [k == 3 ? mkrow(k; c8 = 0.2) : mkrow(k) for k in 1:5]
        b1 = P15BreakEnv.p15_axis_break(c8only, thr)
        @test b1.L_break == 3
        @test b1.per_column[8] == 3
        @test b1.per_column[1] === nothing

        # rho_true (column 1) also degrades, and EARLIER: the envelope ends where EITHER target
        # stops being calibrated, so the earliest crossing across the two gating columns wins.
        both = [k == 3 ? mkrow(k; c8 = 0.2) :
                k == 2 ? mkrow(k; c1 = 0.2) : mkrow(k) for k in 1:5]
        b2 = P15BreakEnv.p15_axis_break(both, thr)
        @test b2.L_break == 2
        @test b2.per_column[1] == 2
        @test b2.per_column[8] == 3

        # The gating tuple is asserted at call time, so a widened gate cannot pass silently.
        @test P15BreakEnv.P15_BREAK_COLUMNS == (1, 8)
        @test length(P15BreakEnv.P15_NUISANCE_COLUMNS) == 6
        @test isempty(intersect(P15BreakEnv.P15_BREAK_COLUMNS,
                                P15BreakEnv.P15_NUISANCE_COLUMNS))

        # Rung 0 is the shared ANCHOR, not a ladder rung, and feeding it in is an error.
        @test_throws AssertionError P15BreakEnv.p15_axis_break(
            [mkrow(P15BreakEnv.P15_RUNG0), mkrow(1)], thr)

        # The reporting record carries the six nuisance ECEs and never the maximum-error field.
        rec = P15BreakEnv.p15_report_only((axis = :spillover, rung = 2,
                                           ece = collect(0.01:0.01:0.08),
                                           ks_p = ones(8), chi2_p = ones(8)))
        @test length(rec.nuisance_ece) == 6
        @test !hasproperty(rec, Symbol(_P15BRK_FIELD_TOKEN))
    end

    # -----------------------------------------------------------------------------------------
    @testset "5 — SC2 is three-valued, per axis, with both LATE sub-cases (D-06, D-13)" begin
        @test P15BreakEnv.P15_SC2_STATES == (:protective, :late, :silent_but_safe)

        # PROTECTIVE — the flag fires before, and exactly at, the breaking rung.
        @test P15BreakEnv.p15_axis_verdict(3, 2) === :protective
        @test P15BreakEnv.p15_axis_verdict(3, 3) === :protective
        # LATE, sub-case (a) — the flag fires, but only AFTER coverage broke.
        @test P15BreakEnv.p15_axis_verdict(2, 4) === :late
        # LATE, sub-case (b) — the flag NEVER fires while coverage did break. This is the case a
        # two-valued verdict would silently score as a pass.
        @test P15BreakEnv.p15_axis_verdict(2, nothing) === :late
        # SILENT-BUT-SAFE — coverage never broke, so no warning was owed, fired or not.
        @test P15BreakEnv.p15_axis_verdict(nothing, nothing) === :silent_but_safe
        @test P15BreakEnv.p15_axis_verdict(nothing, 3) === :silent_but_safe

        for (lb, lo) in ((3, 2), (3, 3), (2, 4), (2, nothing), (nothing, nothing), (nothing, 1))
            @test P15BreakEnv.p15_axis_verdict(lb, lo) in P15BreakEnv.P15_SC2_STATES
        end

        @test P15BreakEnv.p15_sc2_pass(Dict(:a => :protective, :b => :silent_but_safe))
        @test !P15BreakEnv.p15_sc2_pass(Dict(:a => :protective, :b => :late))
        @test !P15BreakEnv.p15_sc2_pass(Dict(:a => :late))

        axes = P15BreakEnv.P15_AXES
        @test length(axes) == 7

        # DEGENERATE READING 1 (D-13): every axis breaks at rung 1 — an empty envelope.
        empty_env = P15BreakEnv.p15_axis_verdicts(
            Dict{Symbol,Any}(a => 1 for a in axes),
            Dict{Symbol,Any}(a => nothing for a in axes))
        @test empty_env.degenerate === :empty_envelope
        @test !empty_env.sc2_pass                      # every axis broke unwarned, so all are LATE

        # DEGENERATE READING 2 (D-13): no axis breaks anywhere — an unbounded envelope.
        unbounded = P15BreakEnv.p15_axis_verdicts(
            Dict{Symbol,Any}(a => nothing for a in axes),
            Dict{Symbol,Any}(a => nothing for a in axes))
        @test unbounded.degenerate === :unbounded_envelope
        @test unbounded.sc2_pass
        @test all(a -> unbounded.verdicts[a] === :silent_but_safe, axes)

        # NEITHER degenerate reading: one axis breaks, and is protected.
        mixed = P15BreakEnv.p15_axis_verdicts(
            Dict{Symbol,Any}(a => (a === :noise ? 3 : nothing) for a in axes),
            Dict{Symbol,Any}(a => (a === :noise ? 2 : nothing) for a in axes))
        @test mixed.degenerate === nothing
        @test mixed.verdicts[:noise] === :protective
        @test mixed.sc2_pass

        # EVERY axis record carries its mechanism label and its prior boundary rung, so the map
        # reads as two overlaid stories: where TRAINING is insufficient (an in-prior axis that
        # breaks) and where the MODEL runs out (an out-of-model axis that breaks).
        for a in axes
            rec = mixed.axes[a]
            @test rec.axis === a
            @test rec.mechanism === P15BreakEnv.P15_AXIS_MECHANISM[a]
            @test rec.mechanism in (:in_prior, :out_of_model)
            @test rec.prior_boundary_rung === P15BreakEnv.P15_PRIOR_BOUNDARY_RUNG[a]
        end
        @test P15BreakEnv.P15_AXIS_MECHANISM[:spillover] === :in_prior
        @test P15BreakEnv.P15_PRIOR_BOUNDARY_RUNG[:spillover] == 2
        @test P15BreakEnv.P15_AXIS_MECHANISM[:background] === :out_of_model
        @test P15BreakEnv.P15_PRIOR_BOUNDARY_RUNG[:background] === nothing

        # An axis missing from either input is an ERROR, never a silent "no break": that is how an
        # incomplete sweep would otherwise be read as an unbounded envelope.
        @test_throws ErrorException P15BreakEnv.p15_axis_verdicts(
            Dict{Symbol,Any}(:noise => 1), Dict{Symbol,Any}(a => nothing for a in axes))
    end

    # -----------------------------------------------------------------------------------------
    @testset "6 — the fire-rate rule uses the MEASURED baseline (D-07)" begin
        # DELIBERATELY NOT the design value: the detector's nominal in-distribution fire rate is
        # 1 - 0.95 = 0.05, and this fixture's MEASURED baseline is more than twice that. A rule
        # written against the nominal value would cross at rung 1 here; the pre-registered rule
        # crosses at rung 3.
        measured_baseline = 0.11
        @test measured_baseline != 1 - _p15brk_quantile_of(ProteinCoLoc)

        rates = [0.10, 0.12, 0.30, 0.60, 0.90]
        @test P15BreakEnv.p15_ood_rung(rates, measured_baseline) == 3
        # …and what the rejected nominal-baseline rule would have said, shown rather than asserted
        # about: it fires two rungs early.
        @test P15BreakEnv.p15_ood_rung(rates, 1 - _p15brk_quantile_of(ProteinCoLoc)) == 1

        # A FLAT fixture exactly at the measured baseline is NOT a firing. An "any image fires"
        # rule would have returned rung 1 on every axis and proved nothing.
        @test P15BreakEnv.p15_ood_rung(fill(0.05, 5), 0.05) === nothing
        @test P15BreakEnv.p15_ood_rung(fill(measured_baseline, 5), measured_baseline) === nothing

        # The margin is binomial design arithmetic, not a chosen number.
        @test P15BreakEnv.P15_OOD_FIRE_MARGIN == 0.025
        @test isapprox(P15BreakEnv.P15_OOD_FIRE_MARGIN,
                       P15BreakEnv._p15_binom_quantile(P15BreakEnv.P15_OOD_N_POS, 0.05, 0.95) /
                       P15BreakEnv.P15_OOD_N_POS - 0.05; atol = 1e-12)

        # C-06: the branch fixtures sit clearly INSIDE their branches, never on the frontier —
        # half the margin below and above, so no assertion here is decided by float representation.
        half = P15BreakEnv.P15_OOD_FIRE_MARGIN / 2
        @test P15BreakEnv.p15_ood_rung([0.05 + half], 0.05) === nothing
        @test P15BreakEnv.p15_ood_rung([0.05 + P15BreakEnv.P15_OOD_FIRE_MARGIN + half], 0.05) == 1
    end

    # -----------------------------------------------------------------------------------------
    @testset "7 — source guard: the maximum-calibration-error field is never read (D-04a)" begin
        @test P15BreakEnv.P15_MCE_FORBIDDEN === true
        @test !isempty(P15BRK_FIELD_SCAN)
        for path in P15BRK_FIELD_SCAN
            hits = _p15brk_scan_file(path, _P15BRK_FIELD_RX)
            @test isempty(hits)
            isempty(hits) || @error(
                "Phase-15 code reads the `$(_P15BRK_FIELD_TOKEN)` field. That statistic is " *
                "pinned at exactly 0.99 on ALL EIGHT columns of the shipped report by the 31 " *
                "EMPTY reliability bins, where the predicted rate is the bin midpoint and the " *
                "observed rate is 0. It is an ARTIFACT OF THE BINNING, not a measurement of the " *
                "net, and quoting it as evidence would be reporting a number that means nothing.",
                file = path, lines = hits)
        end
        # The scan really did look at the file that DECLARES the prohibition, in a trailing
        # comment on a `const` line — the case a whole-line-only comment filter would flag.
        @test _p15brk_path("p15_consts.jl") in P15BRK_FIELD_SCAN
    end

    # -----------------------------------------------------------------------------------------
    @testset "8 — source guard: the detector is never retuned (D-14)" begin
        @test !isempty(P15BRK_RETUNE_SCAN)
        for path in P15BRK_RETUNE_SCAN
            hits = _p15brk_scan_file(path, _P15BRK_RETUNE_RX)
            @test isempty(hits)
            isempty(hits) || @error(
                "Phase-15 code assigns a detector parameter. Lowering the in-distribution " *
                "operating quantile until a LATE axis passes is TUNING A DETECTOR UNTIL IT " *
                "PASSES ITS OWN GATE, and building a new detector channel for the failing " *
                "mechanism is new capability with its own phase. Both are forbidden by D-14: a " *
                "LATE axis is reported as a NAMED LIMIT and the phase closes on it.",
                file = path, lines = hits)
        end

        # `p15_consts.jl` is held to the STRONGER check instead of the zero-hit ban: it byte-forks
        # the amended-v2 rule set, so it RECORDS the shipped operating point on purpose. Absence
        # of that line would not catch a CHANGED number; equality with the shipped constant does.
        @test _p15brk_quantile_of(P15BreakEnv) == _p15brk_quantile_of(ProteinCoLoc)
        @test _p15brk_quantile_of(P15BreakEnv) == 0.95
        # …and it may record ONLY that constant: no threshold assignment of any kind.
        @test isempty(_p15brk_scan_file(_p15brk_path("p15_consts.jl"),
                                        Regex(_P15BRK_THRESHOLD_TOKEN * "\\s*=")))
        for (_, line) in _p15brk_scan_file(_p15brk_path("p15_consts.jl"), _P15BRK_RETUNE_RX)
            @test startswith(line, "const ")
        end
    end

    # -----------------------------------------------------------------------------------------
    @testset "9 — both guards are SHOWN TO FIRE (CONVENTIONS C-04 rule 3)" begin
        # "An assurance never shown non-empty is not a test." The two testsets above are green on
        # a compliant tree, which on its own is equally consistent with a predicate that can never
        # match anything. So run the SAME predicate functions against a file built to violate both.
        mktempdir() do dir
            violating = joinpath(dir, "violating_fixture.jl")
            write(violating, join([
                "x = something." * _P15BRK_FIELD_TOKEN,          # MUST fire: reads the field
                _P15BRK_QUANTILE_TOKEN * " = 0.5",               # MUST fire: retunes the detector
                "# " * _P15BRK_QUANTILE_TOKEN * " = 0.9",        # whole-line comment: must NOT
                "z = 1  # " * _P15BRK_FIELD_TOKEN * " here",     # trailing comment: must NOT
                "w = 2",                                         # clean
            ], "\n"))

            field_hits  = _p15brk_scan_file(violating, _P15BRK_FIELD_RX)
            retune_hits = _p15brk_scan_file(violating, _P15BRK_RETUNE_RX)
            @test length(field_hits)  == 1
            @test length(retune_hits) == 1
            @test first(field_hits)[1]  == 1     # the code line, not the comment lines
            @test first(retune_hits)[1] == 2

            # And a clean file yields nothing, so the predicates are not trivially everything.
            clean = joinpath(dir, "clean_fixture.jl")
            write(clean, "a = 1\nb = a + 1\n")
            @test isempty(_p15brk_scan_file(clean, _P15BRK_FIELD_RX))
            @test isempty(_p15brk_scan_file(clean, _P15BRK_RETUNE_RX))
        end
    end
end
