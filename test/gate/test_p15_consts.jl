#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/test_p15_consts.jl --- the INDEPENDENT assertions on the Phase-15 pre-registration.
#
# `test/gate/p15_consts.jl` already self-checks at load time. That is not enough on its own: a
# self-check lives in the same file it checks, so a single wrong edit can move both halves
# together. This file re-states the load-bearing properties from OUTSIDE the file, against
# `gate_consts_8_v2.jl` loaded separately, and against the simulator and the working tree.
#
# THE GUARD IS SHOWN TO FIRE. CONVENTIONS C-04 rule 3: "an `isempty(intersect(...))` never shown
# non-empty is an assurance, not a test." Testset 2 constructs a derivation that DELIBERATELY
# collides with the spent amended-gate stream and asserts the forbidden set catches it, so the
# disjointness machinery is demonstrated rather than merely exercised on a passing case.
#
# ISOLATED MODULES ARE MANDATORY, NOT STYLE. `p15_consts.jl` and `gate_consts_8_v2.jl` define the
# SAME const names (`SBC_M`, `prod_seed`, …) as the pre-registration `run_gate.jl` has already
# loaded into the test module, and each is guarded on `:SBC_M`, so loading either directly would
# be silently skipped and every comparison below would be reading the wrong file. The module names
# here (`P15Consts`, `P15RefV2`) deliberately avoid `runtests.jl`'s existing `GateConstsV2`.
#
# RUNS BOTH WAYS: standalone (`julia --project=. --threads=auto test/gate/test_p15_consts.jl`) and
# `include`d from `test/runtests.jl`. Every path is built from `@__DIR__`, never from the process
# working directory.

using Test
using ProteinCoLoc
import Random123

module P15Consts
    include(joinpath(@__DIR__, "p15_consts.jl"))
end

module P15RefV2
    include(joinpath(@__DIR__, "gate_consts_8_v2.jl"))
end

"The repository root, resolved from this file rather than from the process working directory."
const P15_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# =============================================================================================
# A SELF-CONTAINED SHA-256, and why this file does not `import SHA`
# =============================================================================================
#
# `SHA` is a Julia stdlib and loads fine under `julia --project=.`, but `Pkg.test()` runs the suite
# in a SANDBOX environment built from the package's `[deps]` plus its `[targets] test` extras —
# where `SHA` is not a direct dependency and `import SHA` fails outright. The obvious fix, adding
# `SHA` to `[extras]`, would EDIT `Project.toml`, which is one of the eight files this phase pins
# as frozen (D-08a: the co-resolution hard gate asserts exact version strings that a re-resolve can
# break). Pinning a file and then editing it in the same commit is not a pre-registration.
#
# So the digest is computed here, in ~50 lines of pure Julia with no dependency at all, and is
# VERIFIED AGAINST THE PUBLISHED NIST TEST VECTORS in the testset below — including the two-block
# and three-block padding paths, which are where a hand-rolled implementation actually goes wrong.
# When `SHA` happens to be loadable (any non-sandbox run) the testset additionally cross-checks
# this implementation against it on the real pinned files.

const _P15_SHA256_K = UInt32[
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2]

_p15_rotr(x::UInt32, n::Int) = (x >> n) | (x << (32 - n))

"""
    _p15_sha256(data) -> String

SHA-256 of `data`, as a 64-character lowercase hex string (FIPS 180-4). Dependency-free; see the
comment block above for why this is not `SHA.sha256`. Verified against the NIST vectors in the
"Phase-15 frozen files are byte-unchanged" testset before any pin is compared against it.
"""
function _p15_sha256(data::AbstractVector{UInt8})
    msg    = Vector{UInt8}(data)
    bitlen = UInt64(length(data)) * 8
    push!(msg, 0x80)
    while length(msg) % 64 != 56
        push!(msg, 0x00)
    end
    for i in 7:-1:0
        push!(msg, UInt8((bitlen >> (8i)) & 0xff))
    end

    h = UInt32[0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
               0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]
    w = Vector{UInt32}(undef, 64)
    for base in 1:64:length(msg)
        @inbounds for i in 0:15
            o = base + 4i
            w[i + 1] = (UInt32(msg[o]) << 24) | (UInt32(msg[o + 1]) << 16) |
                       (UInt32(msg[o + 2]) << 8) | UInt32(msg[o + 3])
        end
        @inbounds for i in 17:64
            s0 = _p15_rotr(w[i - 15], 7)  ⊻ _p15_rotr(w[i - 15], 18) ⊻ (w[i - 15] >> 3)
            s1 = _p15_rotr(w[i - 2], 17)  ⊻ _p15_rotr(w[i - 2], 19)  ⊻ (w[i - 2] >> 10)
            w[i] = w[i - 16] + s0 + w[i - 7] + s1
        end
        a, b, c, d, e, f, g, hh = h[1], h[2], h[3], h[4], h[5], h[6], h[7], h[8]
        @inbounds for i in 1:64
            S1  = _p15_rotr(e, 6) ⊻ _p15_rotr(e, 11) ⊻ _p15_rotr(e, 25)
            ch  = (e & f) ⊻ (~e & g)
            t1  = hh + S1 + ch + _P15_SHA256_K[i] + w[i]
            S0  = _p15_rotr(a, 2) ⊻ _p15_rotr(a, 13) ⊻ _p15_rotr(a, 22)
            maj = (a & b) ⊻ (a & c) ⊻ (b & c)
            t2  = S0 + maj
            hh, g, f, e = g, f, e, d + t1
            d, c, b, a  = c, b, a, t1 + t2
        end
        h[1] += a; h[2] += b; h[3] += c; h[4] += d
        h[5] += e; h[6] += f; h[7] += g; h[8] += hh
    end
    return join(string(x; base = 16, pad = 8) for x in h)
end

_p15_sha256(s::AbstractString) = _p15_sha256(Vector{UInt8}(codeunits(s)))

"""
    _p15_lf_sha256(path) -> String

SHA-256 of a text file with CRLF newlines normalized to LF.

WHY NORMALIZE. `.gitattributes` sets `* text=auto`, so the working-tree bytes of every text file
in this repository are CRLF on Windows and LF on Linux. A raw-bytes pin would therefore fail on
the ubuntu-latest CI this phase exists to build, for a reason that has nothing to do with the
file's contents — and a gate that fails for the wrong reason gets disabled, which is precisely the
outcome D-09 warns about. Normalizing makes the pin content-addressed in the same sense git's own
blob identity is, while still catching any real change to the file.
"""
function _p15_lf_sha256(path::AbstractString)
    return _p15_sha256(Vector{UInt8}(codeunits(replace(String(read(path)), "\r\n" => "\n"))))
end

"""
    _p15_sha_stdlib() -> Union{Module,Nothing}

The `SHA` stdlib if this environment can load it, else `nothing`. `Pkg.test()`'s sandbox cannot
(see the comment block above), so the cross-check that uses this is OPTIONAL and its absence is
reported rather than silently passed over. A plain top-level `import SHA` would abort the whole
suite in the sandbox, which is why this goes through `Base.require` inside a `try`.
"""
function _p15_sha_stdlib()
    try
        return Base.require(@__MODULE__, :SHA)
    catch
        return nothing
    end
end

@testset "Phase-15 seed disjointness (D-11)" begin
    P     = P15Consts
    fresh = P.PROD_SEED_P15[8]

    # --- the six named reserved streams, plus zero ---------------------------------------
    @test fresh != 0
    @test fresh != P.NPE_MASTER_SEED                   # 0xC0FFEE    spike TRAINING
    @test fresh != P.VAL_MASTER_SEED                   # 0x5BC0FFEE  spike VALIDATION
    @test fresh != P.VAL_FIX_SEED                      # 0xF1F7ED    spike FIXTURE
    @test fresh != P.DEFAULT_MASTER_SEED               # 0x1         productionization DATAGEN
    @test fresh != P.RATIO_PAIR_SEED                   # 0x4A7107    shipped ratio pairing
    @test fresh != P.CORPUS_MASTER_SEED                # 0xC05EED    provenance manifest

    # --- burned DEV / FIXTURE streams -----------------------------------------------------
    @test !(fresh in P.DEV_SEEDS)                      # the two F2-remedy diagnostic keys
    for s in (P.P11_DEV_SEED, P.P12_DEV_SEED, P.P12_FIX_SEED, P.P13_DEV_SEED, P.P13_FIX_SEED)
        @test fresh != s
    end
    # The literals are the ones their defining sites carry (a comment naming a seed is not
    # evidence; these were read from disk when the pre-registration was written).
    @test P.VAL_FIX_SEED        === 0x0000_0000_00F1_F7ED
    @test P.RATIO_PAIR_SEED     === 0x0000_0000_004A_7107
    @test P.CORPUS_MASTER_SEED  === 0x0000_0000_00C0_5EED
    @test P.P11_DEV_SEED        === 0x0000_0000_0B11_DE71
    @test P.P12_DEV_SEED        === 0x0000_0000_0B12_DE71
    @test P.P12_FIX_SEED        === 0x0000_0000_0B12_F1F7
    @test P.P13_DEV_SEED        === 0x0000_0000_0B13_DE71
    @test P.P13_FIX_SEED        === 0x0000_0000_0B13_F1F7

    # --- BOTH frozen ship-gate families, all eight seeds ----------------------------------
    # THE SPENT SEED, NAMED EXPLICITLY. `run_gate.jl:47-54` binds the amended pre-registration to
    # exactly ONE run on PROD_SEED_V2[8]; that run produced
    # artifacts/amended_v2/grid_8/gate_report_8.jld2 and the seed is consumed.
    @test fresh != P.PROD_SEED_V2[8]
    for G in (4, 8, 16, 32)
        @test fresh != P.PROD_SEED[G]
        @test fresh != P.PROD_SEED_V2[G]
        @test P.PROD_SEED_P15[G] != P.PROD_SEED[G]
        @test P.PROD_SEED_P15[G] != P.PROD_SEED_V2[G]
    end
    @test isempty(intersect(Set(values(P.PROD_SEED_P15)), Set(values(P.PROD_SEED))))
    @test isempty(intersect(Set(values(P.PROD_SEED_P15)), Set(values(P.PROD_SEED_V2))))
    @test length(unique(values(P.PROD_SEED_P15))) == 4          # distinct per grid
    @test !(fresh in P._p15_forbidden())

    # The recomputed v2 family really is the frozen one (so forbidding it means something).
    @test P.PROD_SEED_V2 == P15RefV2.PROD_SEED_V2
    @test P.PROD_SEED    == P15RefV2.PROD_SEED

    # --- C-04 immunity is a KEY difference, not an index-range difference ------------------
    @test !(P.P15_SALT in P.P15_REPO_SALTS)
    @test length(P.P15_REPO_SALTS) == 11 && allunique(P.P15_REPO_SALTS)
    @test (UInt64(P.P15_MASTER) ⊻ P.P15_SALT) !=
          (UInt64(P15RefV2.AMEND_MASTER) ⊻ P15RefV2.AMEND_SALT)
    @test (UInt64(P.P15_MASTER) ⊻ P.P15_SALT) !=
          (UInt64(P15RefV2.PROD_MASTER) ⊻ P15RefV2.PROD_SALT)
    # The two frozen salts are in the inventory the fresh salt is checked against.
    @test P15RefV2.AMEND_SALT in P.P15_REPO_SALTS
    @test P15RefV2.PROD_SALT  in P.P15_REPO_SALTS
end

@testset "Phase-15 seed guard FIRES (C-04 rule 3)" begin
    # THE NEGATIVE FIXTURE. Build the SAME derivation the pre-registration uses, but keyed on the
    # SPENT amended-gate pair and with NO redraw loop, and show its first draw lands inside
    # `_p15_forbidden()`. Without this the disjointness testset above only ever demonstrates that
    # a passing case passes — an assurance, not a test.
    _colliding(master, salt, G) =
        rand(Random123.Philox4x(UInt64, (UInt64(master) ⊻ UInt64(salt), UInt64(G))), UInt64)

    collided = _colliding(P15RefV2.AMEND_MASTER, P15RefV2.AMEND_SALT, 8)
    @test collided == P15RefV2.PROD_SEED_V2[8]        # it really is the spent seed
    @test collided in P15Consts._p15_forbidden()      # …and the guard CATCHES it
    @test collided != P15Consts.PROD_SEED_P15[8]      # the redraw loop moved P15 off it

    # The same demonstration on the v1 family, so the guard is shown to catch both.
    v1 = _colliding(P15RefV2.PROD_MASTER, P15RefV2.PROD_SALT, 8)
    @test v1 == P15RefV2.PROD_SEED[8]
    @test v1 in P15Consts._p15_forbidden()

    # And a seed the guard must NOT catch, so the set is not trivially everything.
    @test !(P15Consts.PROD_SEED_P15[8] in P15Consts._p15_forbidden())
end

@testset "Phase-15 stream and its global half (D-11, Pitfall 5)" begin
    P = P15Consts
    @test P.prod_seed(8) == P.PROD_SEED_P15[8]
    @test P.prod_rng(8) isa Random123.Philox4x
    @test rand(P.prod_rng(8), UInt64, 4) == rand(P.prod_rng(8), UInt64, 4)   # reproducible
    @test rand(P.prod_rng(8), UInt64, 4) != rand(P.prod_rng(4), UInt64, 4)   # per-grid disjoint
    @test P.prod_seed(8) != P15RefV2.prod_seed(8)      # the aliases do NOT ride the spent stream

    # THE GLOBAL HALF (harness.jl:49-60, :82). `sampleposterior` threads no rng and draws its flow
    # base samples from the GLOBAL stream, so `gate_global_seed(G)` is the second, easily-dropped
    # half of reproducibility. It is derived from `prod_seed(G)`, which means aliasing `prod_seed`
    # onto P15 is what carries it. Recompute the harness rule in BOTH modules and assert they part.
    _global_half(M, G) = rand(Random123.Philox4x(UInt64, (M.prod_seed(G), UInt64(1))), UInt64)
    @test _global_half(P, 8) != _global_half(P15RefV2, 8)
    @test _global_half(P, 8) == _global_half(P, 8)                     # deterministic
    @test _global_half(P, 8) != _global_half(P, 4)                     # per-grid
    @test _global_half(P, 8) != P.prod_seed(8)                         # counter-1 sibling, not the seed
end

@testset "Phase-15 inherits the amended v2 rule set exactly (D-01)" begin
    # WHAT THIS MAKES CHECKABLE: "Phase 15 scores its SBC arm under the same amended rules the
    # shipped report was scored under" — otherwise a claim, here a per-constant equality.
    skip = (:GATE_CONSTS_VERSION, :GATE_AMENDMENT_DOC)
    refnames = filter(names(P15RefV2; all = true)) do n
        s = String(n)
        occursin(r"^[A-Z][A-Z0-9_]*$", s) && isdefined(P15RefV2, n) && isconst(P15RefV2, n)
    end
    # Guard against a filter that silently matches nothing and turns this into a vacuous pass.
    @test length(refnames) >= 40
    @test :SBC_M in refnames && :PROD_SEED_V2 in refnames && :BF_CORR_MIN in refnames

    for n in refnames
        n in skip && continue
        @test isdefined(P15Consts, n)
        @test getfield(P15Consts, n) == getfield(P15RefV2, n)
    end

    # The two deliberate differences, asserted as differences so they cannot drift back silently.
    @test P15Consts.GATE_CONSTS_VERSION == 15 && P15RefV2.GATE_CONSTS_VERSION == 2
    @test P15Consts.GATE_AMENDMENT_DOC != P15RefV2.GATE_AMENDMENT_DOC
    @test P15Consts.GATE_AMENDMENT_DOC ==
          ".planning/phases/15-calibration-operating-envelope-and-ci-gate/15-CONTEXT.md"
    @test isfile(joinpath(P15_REPO_ROOT, P15Consts.GATE_AMENDMENT_DOC))
    @test P15Consts.GATE_CONSTS_GRID == 8

    # The consequences that make the inheritance load-bearing for the SBC arm.
    @test P15Consts.SBC_IMSIZE === :mixture
    @test P15Consts.SBC_REQUIRE_IMSIZE_PROVENANCE
    @test P15Consts.SBC_MULTIPLICITY === :holm
    @test P15Consts.SBC_L == 999 && P15Consts.SBC_BINS == 50
end

@testset "Phase-15 ladders derive from the simulator prior (D-03)" begin
    P = P15Consts
    # The literals in the pre-registration ARE the simulator's prior extrema — so the file can stay
    # dependency-free without letting its numbers drift away from `src/amortized/simulator.jl`.
    @test P.P15_PRIOR_BOUNDARY.spillover        == maximum(ProteinCoLoc.SPILLOVER_PRIOR)
    @test P.P15_PRIOR_BOUNDARY.autofluorescence == maximum(ProteinCoLoc.AUTOFLUORESCENCE_PRIOR)
    @test P.P15_PRIOR_BOUNDARY.registration     == maximum(ProteinCoLoc.SHIFT_PRIOR)

    for lad in (P.P15_SPILLOVER_LADDER, P.P15_SHIFT_LADDER, P.P15_AF_LADDER)
        @test length(lad) == P.P15_RUNGS
        @test issorted(collect(lad)) && allunique(lad)     # strictly increasing
    end
    # RUNG 2 IS THE BOUNDARY on every in-prior axis, and rung 1 is strictly inside it: a ladder
    # that starts at the edge cannot distinguish "breaks where we trained" from "breaks past it".
    @test P.P15_SPILLOVER_LADDER[2] == P.P15_PRIOR_BOUNDARY.spillover
    @test P.P15_SHIFT_LADDER[2]     == P.P15_PRIOR_BOUNDARY.registration
    @test P.P15_AF_LADDER[2]        == P.P15_PRIOR_BOUNDARY.autofluorescence
    @test P.P15_SPILLOVER_LADDER[1] < P.P15_PRIOR_BOUNDARY.spillover
    @test P.P15_SHIFT_LADDER[1]     < P.P15_PRIOR_BOUNDARY.registration
    @test P.P15_AF_LADDER[1]        < P.P15_PRIOR_BOUNDARY.autofluorescence

    # The spillover endpoint is the SIMULATOR'S OWN guard (simulator.jl:232), not a choice: the
    # ladder cannot be extended past it without `simulate_pair` throwing. Shown by construction.
    @test P.P15_SPILLOVER_LADDER[end] <= 1.0
    let θ = (ρ_true = 0.5, spillover = 1.0 + 1e-6, autofluorescence = 0.05,
             label_efficiency = 0.9, shift_dx = 0.0, shift_dy = 0.0, noise = 0.5,
             chromatic_eps = 0.0)
        @test_throws ArgumentError ProteinCoLoc.simulate_pair(
            Random123.Philox4x(UInt64, (UInt64(1), UInt64(0))), θ; imsize = (64, 64))
    end

    # Mechanism typing (D-02, D-02a) — the seven axes, and background reclassified out-of-model.
    @test length(P.P15_AXES) == 7
    @test Tuple(keys(P.P15_AXIS_MECHANISM)) == P.P15_AXES
    @test keys(P.P15_AXIS_MECHANISM) == keys(P.P15_PRIOR_BOUNDARY_RUNG)
    @test P.P15_AXIS_MECHANISM.background       === :out_of_model
    @test P.P15_AXIS_MECHANISM.autofluorescence === :in_prior
    @test P.P15_AXIS_MECHANISM.optics           === :out_of_model   # σ_psf is FIXED, not a prior
    @test P.P15_PRIOR_BOUNDARY_RUNG.optics      === nothing         # nothing = nothing to MARK
    @test P.P15_PRIOR_BOUNDARY_RUNG.background  === nothing
    @test count(==(:in_prior), values(P.P15_AXIS_MECHANISM)) == 3
    for ax in P.P15_AXES
        if P.P15_AXIS_MECHANISM[ax] === :in_prior
            @test P.P15_PRIOR_BOUNDARY_RUNG[ax] == 2
        else
            @test P.P15_PRIOR_BOUNDARY_RUNG[ax] === nothing
        end
    end
    # The four shipped families are exactly the out-of-model axes (misspec.jl:181-184).
    @test Tuple(k for k in P.P15_AXES if P.P15_AXIS_MECHANISM[k] === :out_of_model) ==
          (:texture, :noise, :optics, :background)
end

@testset "Phase-15 tier split and D-13 declarations (D-12, D-13)" begin
    P = P15Consts
    # "Reported is not gated" is machine-checkable, not prose.
    @test isempty(intersect(P.P15_GATING_CONSTANTS, P.P15_REPORTING_ONLY_CONSTANTS))
    for n in (P.P15_GATING_CONSTANTS..., P.P15_REPORTING_ONLY_CONSTANTS...)
        @test isdefined(P, n)
    end
    @test :P15_ECE_MARGIN in P.P15_GATING_CONSTANTS
    @test :P15_NUISANCE_COLUMNS in P.P15_REPORTING_ONLY_CONSTANTS
    @test P.P15_KS_REPORTED && P.P15_CHI2_REPORTED          # computed, reported, NOT gating
    @test P.P15_MCE_FORBIDDEN                               # empty-bin artifact, pinned at 0.99

    # D-05 columns: two gate, six are reported only, together they are all eight.
    @test P.P15_BREAK_COLUMNS == (1, 8)                     # ρ_true and Δρ (sbc.jl:123-124)
    @test isempty(intersect(P.P15_BREAK_COLUMNS, P.P15_NUISANCE_COLUMNS))
    @test sort(vcat(collect(P.P15_BREAK_COLUMNS), collect(P.P15_NUISANCE_COLUMNS))) == collect(1:8)

    # D-13, written down before any rung exists.
    @test P.P15_ITERATION_ALLOWANCE == 1
    @test P.P15_SC2_STATES == (:protective, :late, :silent_but_safe)
    @test length(P.P15_SC2_STATES) == 3 && allunique(P.P15_SC2_STATES)
    for s in (P.P15_DEGENERATE_EMPTY_ENVELOPE_READING, P.P15_DEGENERATE_UNBOUNDED_ENVELOPE_READING,
              P.P15_SC2_PASS_RULE, P.P15_LATE_AXIS_RESPONSE, P.P15_GRID16_DROP_REASON,
              P.P15_GOLDEN_HONESTY)
        @test s isa String && !isempty(s)
    end
    # The SC2 rule and the LATE response say the things D-06/D-14 require them to say.
    @test occursin("LATE", P.P15_SC2_PASS_RULE)
    @test occursin("SILENT-BUT-SAFE", P.P15_SC2_PASS_RULE)
    @test occursin("OOD_ID_QUANTILE", P.P15_LATE_AXIS_RESPONSE)   # naming what is FORBIDDEN

    # D-10a: dropped with its evidence, carried as a named limit rather than silently omitted.
    @test P.P15_GRID16_CONTRAST === :dropped
    @test occursin("training_imsize_provenance", P.P15_GRID16_DROP_REASON)
    @test occursin("UNTESTED", P.P15_GRID16_DROP_REASON)
    @test P.P15_REPORTED_NET == "artifacts/amended_v2/grid_8"
    @test P.P15_ARTIFACTS_ROOT == "artifacts/amended_v2"

    # The D-09 golden bars ride the smoke fixture and are honest about what they prove.
    @test P.P15_GOLDEN_M == P.SBC_FIX_M && P.P15_GOLDEN_L == P.SBC_FIX_L
    @test P.P15_GOLDEN_BINS == P.SBC_FIX_BINS
    @test (P.P15_GOLDEN_L + 1) % P.P15_GOLDEN_BINS == 0
    @test P.P15_GOLDEN_IMSIZE == (256, 256)
    @test P.P15_GOLDEN_ECE_TOL == 1e-8
    @test occursin("CHANGE", P.P15_GOLDEN_HONESTY)
end

@testset "Phase-15 frozen files are byte-unchanged (D-02a, D-08a)" begin
    P = P15Consts

    # --- the digest itself is verified BEFORE any pin is compared against it ---------------
    # Published FIPS 180-4 / NIST test vectors. The three lengths are chosen to exercise the three
    # padding paths a hand-rolled SHA-256 gets wrong: empty (pad-only block), 3 bytes (single
    # block), 56 bytes (length spills into a SECOND block) and 112 bytes (three blocks).
    @test _p15_sha256("") ==
          "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    @test _p15_sha256("abc") ==
          "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    @test _p15_sha256("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
          "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
    @test _p15_sha256("abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmnhijklmno" *
                      "ijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu") ==
          "cf5b16a778af8380036ce59e7b0492370b249b11e8f07a51afac45037afee9d1"
    @test length(_p15_sha256("x")) == 64
    @test _p15_sha256("x") != _p15_sha256("y")            # it is not a constant function

    # Optional cross-check against the SHA stdlib. Available under `julia --project=.`, NOT inside
    # `Pkg.test()`'s sandbox — reported honestly rather than silently skipped.
    let sha = _p15_sha_stdlib()
        if sha === nothing
            @info "Phase-15: SHA stdlib not loadable in this environment (Pkg.test sandbox) — " *
                  "the hand-rolled digest is verified by the NIST vectors above only"
        else
            for entry in P.P15_FROZEN_FILE_SHA256
                bytes = Vector{UInt8}(codeunits(replace(
                    String(read(joinpath(P15_REPO_ROOT, entry.path))), "\r\n" => "\n")))
                @test _p15_sha256(bytes) == bytes2hex(sha.sha256(bytes))
            end
        end
    end

    @test P.P15_FROZEN_HASH_NORMALIZATION === :lf
    @test length(P.P15_FROZEN_FILE_SHA256) == 8
    for entry in P.P15_FROZEN_FILE_SHA256
        path = joinpath(P15_REPO_ROOT, entry.path)
        @test isfile(path)
        # WHY THIS MATTERS, stated in the failure message rather than left to the reader:
        #   misspec.jl   — appending to OOD_FAMILIES silently REWRITES the shipped, PASSED OOD
        #                  arm's combined_auc and its pass conjunction (D-02a).
        #   Manifest.toml — a re-resolved manifest can satisfy [compat] and still fail the exact
        #                  co-resolution strings runtests.jl:32-49 asserts (D-08a).
        @test _p15_lf_sha256(path) == entry.sha256
    end
    # The pins really do cover the four files this phase must not touch, plus both environments.
    pinned = Set(e.path for e in P.P15_FROZEN_FILE_SHA256)
    for p in ("test/gate/misspec.jl", "test/gate/sbc.jl", "test/gate/harness.jl",
              "test/gate/gate_consts_8_v2.jl", "Project.toml", "Manifest.toml",
              "spike/Project.toml", "spike/Manifest.toml")
        @test p in pinned
    end
    # The normalization is what makes the pin platform-independent: on a CRLF checkout the RAW
    # digest differs from the pinned one, and that difference must not be mistaken for a change.
    @test _p15_lf_sha256(joinpath(P15_REPO_ROOT, "test/gate/misspec.jl")) ==
          P.P15_FROZEN_FILE_SHA256.misspec.sha256
end

@testset "Phase-15 ECE margin is design arithmetic (D-04a)" begin
    P = P15Consts
    # The margin IS q95 − mean of the measured null; nothing here is an outcome of any net.
    @test P.P15_ECE_MARGIN ≈ P.P15_ECE_NULL_Q95 - P.P15_ECE_NULL_MEAN atol = 1e-12
    @test P.P15_ECE_NULL_Q95 > P.P15_ECE_NULL_MEAN > 0
    @test P.P15_ECE_NULL_B >= 2000 && P.P15_ECE_NULL_SEED > 0
    # The closed form E[ECE₀] ≈ 0.336/√M reproduces the measurement — two independent routes.
    @test abs(P.P15_ECE_NULL_MEAN - 0.336 / sqrt(P.P15_SWEEP_M)) / P.P15_ECE_NULL_MEAN < 0.05
    # …as does the independent B = 2000 research probe, within its own Monte-Carlo error.
    @test abs(P.P15_ECE_NULL_MEAN - P.P15_ECE_NULL_MEAN_RESEARCH) < 0.0015
    @test abs(P.P15_ECE_NULL_Q95  - P.P15_ECE_NULL_Q95_RESEARCH)  < 0.0030

    # The break threshold is a FORM over a rung-0 anchor, never a stored number.
    @test P.p15_break_threshold(0.0) == P.P15_ECE_MARGIN
    @test P.p15_break_threshold(0.03) ≈ 0.03 + P.P15_ECE_MARGIN atol = 1e-12
    # D-04a's arithmetic, made visible: the shipped Δρ ECE sits BELOW its own M=2000 null mean, so
    # anchoring there would declare a break on ~90 % of perfectly calibrated rungs.
    @test 0.0042894736842104715 < 0.00751

    # The OOD fire margin is the binomial null at the pre-registered n, recomputed here.
    @test P._p15_binom_quantile(200, 0.05, 0.95) == 15
    @test P._p15_binom_quantile(P.P15_OOD_N_POS, 1 - P.OOD_ID_QUANTILE, 0.95) /
          P.P15_OOD_N_POS - (1 - P.OOD_ID_QUANTILE) ≈ P.P15_OOD_FIRE_MARGIN atol = 1e-12
    # Sanity anchors for the hand-rolled quantile (median and both tails of a symmetric case).
    @test P._p15_binom_quantile(10, 0.5, 0.5) == 5
    @test P._p15_binom_quantile(100, 0.5, 0.999) > 50
    @test P._p15_binom_quantile(200, 0.05, 0.05) < 10
    # …and the sizing argument behind n_pos = 200: at the shipped n_pos = 40 one binomial SE is
    # 0.034, which swallows the margin entirely; at 200 it is 0.015, which does not.
    @test sqrt(0.05 * 0.95 / 40) > P.P15_OOD_FIRE_MARGIN
    @test sqrt(0.05 * 0.95 / P.P15_OOD_N_POS) < P.P15_OOD_FIRE_MARGIN

    # Budget bars.
    @test P.P15_SWEEP_M == 1000 && P.P15_RUNGS == 5 && P.P15_RUNG0 == 0
    @test P.P15_SWEEP_MIN_THREADS >= 16
    @test P.P15_SWEEP_WALLCLOCK_CEILING_HOURS > 8.5
end
