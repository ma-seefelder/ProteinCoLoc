# Phase 14: Decision and Abstention Layer — Pattern Map

**Mapped:** 2026-08-03
**Files analyzed:** 16 new files (7 module + 4 runner + 5 test) — list reconciled against
`14-RESEARCH.md` §D and `14-VALIDATION.md` Wave 0 (see *File-list reconciliation* below)
**Analogs found:** 16 / 16 (13 exact/role-match, 3 pattern-only)
**Language:** Julia. `spike/p14/` does not exist yet (`ls spike/p14` → no such directory).

---

## 0. THE HARD CONSTRAINT, FIRST

**`spike/Project.toml` and `spike/Manifest.toml` are asserted byte-unchanged by a RUNNING test.**
Every pattern below was checked against the frozen sixteen. **No pattern in this document introduces
a package.** The only packages any recommended analog uses are `Test`, `Pkg`, `Statistics`, `Dates`,
`SHA`, `LinearAlgebra` (all stdlib, no `Project.toml` entry needed — they ride the default
`LOAD_PATH` `@stdlib` entry, per `test_p12_decoupling.jl:171-174`) plus `JLD2`, `StatsBase`,
`Random123`, `Flux`, `NeuralEstimators`, `Distributions`, `CairoMakie`, `Images` — all eight of
which ARE among the frozen sixteen.

Verbatim, from `spike/test/test_p12_decoupling.jl:170-183` — **the environment freeze**:

```julia
@testset "the spike environment is byte-frozen" begin
    # 12-PATTERNS §0.1: the spatial map's stdlib reach (LinearAlgebra, SparseArrays) rides
    # the default LOAD_PATH `@stdlib` entry and needs no Project.toml entry, so ANY movement
    # in these two files means a `Pkg.add`/`Pkg.resolve` happened -- which is the thing the
    # frozen manifest exists to prevent.
    if _HAS_GIT
        @test success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                          dir = P12_REPO_ROOT))
        @test isempty(readchomp(Cmd(`git status --porcelain -- spike/Project.toml spike/Manifest.toml`;
                                    dir = P12_REPO_ROOT)))
    else
        @info "git unavailable — skipping the executable spike-environment freeze assertions"
    end
end
```

Verbatim, `:185-208` — **the exactly-sixteen dependency set** (note the banner: *"each phase states
its own claim in its own file"*, which is the licence for `test_p14_decoupling.jl` to restate it):

```julia
@testset "the dependency name set is exactly the frozen sixteen (resolve-risk clause k)" begin
    # (k) RESOLVE-RISK GATE (Phase 12). Re-asserted FROM THIS FILE rather than by editing
    # `runtests.jl:114-127`: the Phase-11 block there is a committed record of what Phase 11
    # froze, and each phase states its own claim in its own file.
    @test Set(keys(Pkg.project().dependencies)) == Set([
        "BenchmarkTools", "CSV", "CairoMakie", "CoordinateTransformations", "DataFrames",
        "Distributions", "Flux", "HypothesisTests", "ImageFiltering", "ImageTransformations",
        "Images", "Interpolations", "JLD2", "NeuralEstimators", "Random123", "StatsBase"])
    @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"

    # ... Phase-12-specific half: each banned package asserted BY NAME, one test each, so the
    # failure says WHICH one appeared.
    @test !haskey(Pkg.project().dependencies, "AbstractGPs")
    @test !haskey(Pkg.project().dependencies, "KernelFunctions")
    @test !haskey(Pkg.project().dependencies, "GaussianRandomFields")
end
```

> **Phase-14 analogue of the by-name half:** `@test !haskey(Pkg.project().dependencies,
> "ConformalPrediction")` — plus `MLJ` / `ROCAnalysis` / `MultipleTesting`, which
> `14-RESEARCH.md` §"Explicitly NOT to be added" names. The by-name form is the pattern: one
> `@test` each so the failure names the offender.

Verbatim, `:291-296` — **the corpus tree byte-unchanged** (D-03's bound), inside the sealed-holdout
testset:

```julia
    if _HAS_GIT
        @test success(Cmd(`git diff --quiet HEAD -- corpus`; dir = P12_REPO_ROOT))
        @test isempty(readchomp(Cmd(`git status --porcelain -- corpus`; dir = P12_REPO_ROOT)))
    else
        @info "git unavailable — skipping the executable corpus/ byte-equality assertion"
    end
```

Verbatim, `:351-353` — **CPU-only**:

```julia
@testset "decoupling checks ran CPU-only" begin
    @test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))
end
```

And `:152-168`, **`src/` byte-unchanged, including the untracked case** — the assertion Phase 14 must
also restate because D-01 keeps it in the research lane:

```julia
@testset "src/ is byte-unchanged (CLAUDE.md hard constraint)" begin
    if _HAS_GIT
        @test success(Cmd(`git diff --quiet HEAD -- src`; dir = P12_REPO_ROOT))
        @test success(Cmd(`git diff --quiet HEAD -- src/results.jl`; dir = P12_REPO_ROOT))
        @test success(Cmd(`git diff --quiet HEAD -- src/amortized/summary.jl`; dir = P12_REPO_ROOT))
        @test success(Cmd(`git diff --quiet HEAD -- src/amortized/architecture.jl`; dir = P12_REPO_ROOT))
        @test success(Cmd(`git diff --quiet HEAD -- src/amortized/local_map.jl`; dir = P12_REPO_ROOT))
        # THE UNTRACKED CASE, WHICH `git diff HEAD` CANNOT SEE. A NEW file dropped under
        # `src/` is invisible to a diff against HEAD and would pass every assertion above
        # while being exactly the productionization-by-stealth this constraint forbids.
        @test isempty(readchomp(Cmd(`git status --porcelain -- src`; dir = P12_REPO_ROOT)))
    else
        @info "git unavailable — skipping the executable src/ byte-equality assertions"
    end
end
```

**Do not break these. `test_p12_decoupling.jl` scans `("spike/p12", r"^.*\.jl$")` and
`("spike/test", r"^test_p12_.*\.jl$")` — a `spike/p14/` tree and `test_p14_*.jl` files are NOT in
its scanned surface (`:78-88`), so Phase 14 files cannot make the Phase-12 corpus-token scan red.
The three assertions above are repo-global, however, and Phase 14 CAN break all three.**

---

## 1. File-list reconciliation (RESEARCH vs VALIDATION vs the task brief)

The task brief's file list is **correct and complete except for two additions**, both of which come
from `14-VALIDATION.md`:

| File | Source | Note |
|---|---|---|
| `spike/p14/consts.jl` | VALIDATION Wave 0 (names all 15 constants) | ✓ as briefed |
| `spike/p14/posterior.jl` | RESEARCH §Code Examples `p14_class_posterior` | ✓ |
| `spike/p14/fdr.jl` | RESEARCH §Code Examples `p14_bayes_fdr` | ✓ |
| `spike/p14/conformal.jl` | RESEARCH §Code Examples `p14_conformal_quantile` / `_set` / `p14_lac_score` | ✓ |
| `spike/p14/fuse.jl` | RESEARCH §Code Examples `p14_fuse` + `p14_ood_state` (§D.4) | ✓ |
| `spike/p14/result.jl` | RESEARCH §D.6 (`P14Result` + `P14BatchDecision`) | ✓ |
| `spike/p14/decide.jl` | RESEARCH §D.6 (`decide_coloc`, batch-level) | ✓ |
| `run_p14_conformal.jl` | VALIDATION SC1-d | ✓ |
| `run_p14_riskcoverage.jl` | VALIDATION SC3-a/b/c | ✓ |
| `run_p14_ood_arm.jl` | VALIDATION SC3-d | ✓ |
| `run_p14_fdr_check.jl` | VALIDATION SC1-b/c | ✓ |
| `test_p14_consts.jl`, `_decoupling.jl`, `_provenance.jl`, `_fdr.jl`, `_conformal.jl`, `_fuse.jl`, `_posterior.jl` | VALIDATION | ✓ |
| **`spike/test/fixtures/` additions** | **VALIDATION Wave 0** — "the unequal-posterior class-order fixture and the FDR discriminating fixture" | **MISSING from the brief's table.** Existing fixtures are `.jld2` goldens (`p11_stage6_golden.jld2`, `p12_stage1_golden.jld2`); the Phase-14 fixtures described are *hand-built literals*, so the analog is the in-file `const RES_FIX_*` idiom of `test_p13_result.jl:68-80`, not a `.jld2`. |
| **A `p14_ood_state` home** | RESEARCH §D.4 places it with the fusion. Recommend `fuse.jl`; SC2-c greps for `\.flag == false` **outside `p14_ood_state`**, so the two must be greppable together. |

Nothing in RESEARCH or VALIDATION contradicts the brief. `run_p14_fdr_check.jl` is listed in
VALIDATION's SC1-b/c rows but **not** in its "Before `/bm:verify-work`" runner list (which names only
conformal / riskcoverage / ood_arm) — flagged for the planner, not resolved here.

---

## 2. File Classification

| New file | Role | Data flow | Closest analog (absolute path) | Match |
|---|---|---|---|---|
| `spike/p14/consts.jl` | config / pre-registration | frozen constants | `C:\...\ProteinCoLoc\spike\p13\consts.jl` | **exact** |
| `spike/p14/posterior.jl` | utility (numeric) | transform | `C:\...\ProteinCoLoc\spike\p13\net.jl` (`three_way_probs`, `three_way_log_bf`) + `spike\p13\labels.jl` (`three_way_label`) | role-match |
| `spike/p14/fdr.jl` | utility (numeric) | batch transform | `C:\...\ProteinCoLoc\spike\validation\ood.jl:319` (`roc_auc`) — the house *hand-rolled, no-dependency, sort+cumsum numeric helper* | role-match |
| `spike/p14/conformal.jl` | utility (numeric) | batch transform | `C:\...\ProteinCoLoc\spike\validation\ood.jl:344` (`id_threshold`) + `spike\simulator\ghat.jl:110` (`ghat`, frozen clamped lookup) | role-match |
| `spike/p14/fuse.jl` | service (decision rule) | transform | `C:\...\ProteinCoLoc\src\amortized\ood.jl:363` (`ood_verdict` OR-fusion) + `C:\...\ProteinCoLoc\src\amortized\local_map.jl:93` (`_recorded_ood_threshold`) + `spike\p13\labels.jl:121` (`three_way_label`, the branch-per-decision shape) | role-match |
| `spike/p14/result.jl` | model (result type) | container | `C:\...\ProteinCoLoc\spike\p13\result.jl` | **exact** |
| `spike/p14/decide.jl` | controller (entry point) | request-response | `C:\...\ProteinCoLoc\src\amortized\api.jl:107` (`colocalization_amortized`, the src-shaped signature) + `spike\p13\run_p13_realimage.jl` (spike-lane composition) | role-match |
| `spike/p14/run_p14_*.jl` (×4) | runner script | batch → `.jld2` | `C:\...\ProteinCoLoc\spike\p13\run_three_way_gate.jl` | **exact** |
| `spike/test/test_p14_consts.jl` | test | assertion | `C:\...\ProteinCoLoc\spike\test\test_p13_consts.jl` | **exact** |
| `spike/test/test_p14_decoupling.jl` | test | source scan + git | `C:\...\ProteinCoLoc\spike\test\test_p12_decoupling.jl` | **exact** |
| `spike/test/test_p14_provenance.jl` | test | artifact load + sha | `C:\...\ProteinCoLoc\spike\test\test_p13_consts.jl` §"tau is measured and locked (Tier 2)" (`:121-144`) + `spike\p13\run_three_way_gate.jl:361-370` (the pinned-both-sides sha assertions) | role-match |
| `spike/test/test_p14_fdr.jl` / `_conformal.jl` / `_posterior.jl` | test | fixture assertion | `C:\...\ProteinCoLoc\spike\test\test_p13_result.jl` (fixtures-once + one outer `@testset`) | **exact** |
| `spike/test/test_p14_fuse.jl` | test | truth-table | `C:\...\ProteinCoLoc\spike\test\test_p13_labels.jl` (**NOT VERIFIED** — not opened this session; chosen by name/role. `test_p13_result.jl` is the verified fallback.) | pattern-only |

---

## 3. Pattern Assignments

### 3.1 `spike/p14/consts.jl` — analog `spike/p13/consts.jl`

**(a) The two-tier banner + append-never-overwrite rule** (`spike/p13/consts.jl:21-44`), verbatim —
this is the structure Wave 0 asks for:

```julia
# spike/p13/consts.jl --- Phase-13 Tier-1 pre-registration (D-04).
#
# THE ANTI-SNOOPING CONTRACT (D-04). EVERY THRESHOLD, LADDER, SEED, STATISTIC AND
# DESIGN CHOICE THAT ANY REPORTED PHASE-13 NUMBER WILL BE SCORED AGAINST IS LOCKED
# IN THIS ONE FILE AND COMMITTED BEFORE ANY PHASE-13 CODE RUNS ...
#
# TWO-TIER STRUCTURE -- READ THIS BEFORE EDITING ANYTHING.
#   TIER 1 is the single guard block in THIS file, committed before anything runs.
#   TIER 2 is `P13_TAU` -- the ONE number D-06 says must be MEASURED, not chosen. It
#     is DELIBERATELY ABSENT from this file. Plan 13-10 runs the tau probe ... and
#     APPENDS the measured value in a SECOND guard block keyed on a DIFFERENT sentinel
#     const, carrying a one-line provenance comment naming the probe artifact ...
#   NO CONSTANT IN EITHER BLOCK IS EVER *EDITED*. Constants are only ever APPENDED,
#   so the git history of this file is itself the pre-registration audit trail. A
#   diff that MODIFIES a line below is, by construction, a pre-registration breach.
```

**(b) The guard-block sentinel rule** (`:102-109`) — **this is a recorded failure; do not repeat
it.** The sentinel must be a name **this file alone declares**, never a seed (later phases mirror
seeds for disjointness, and guarding on a mirrored name made a full-suite run skip ~104 constants):

```julia
# THE BODY SENTINEL IS `P13_DECLARED_DEVIATIONS` (declared unconditionally at the foot of this
# block), NOT `P13_DEV_SEED`. A seed is the WORST possible sentinel in this project: R-5 obliges
# every Phase-N pre-registration to re-declare prior phases' seeds in order to assert stream
# disjointness ... Guarding on the mirrored name made a full-suite run skip this entire body ...
if !isdefined(@__MODULE__, :P13_DECLARED_DEVIATIONS)
    import Random123: Philox4x     # counter-based RNG; the fresh disjoint Phase-13 streams
```

For Phase 14 the sentinel must therefore be a P14-only name (e.g. `P14_DECLARED_DEVIATIONS`), and
**every caller uses the guarded include** — `isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) ||
include(joinpath(@__DIR__, "consts.jl"))`.

**(c) The forbidden-seed inventory + fresh streams** (`:116-218`), condensed — copy the *shape*:

```julia
    const DEFAULT_MASTER_SEED = 0x0000_0000_0000_0001  # productionization DATAGEN (datagen.jl:52)
    const NPE_MASTER_SEED     = 0x0000_0000_00C0_FFEE  # spike TRAINING stream (npe/train_npe.jl:65)
    const VAL_MASTER_SEED     = 0x0000_0000_5BC0_FFEE  # spike VALIDATION stream (validation/consts.jl:76)
    const VAL_FIX_SEED        = 0x0000_0000_00F1_F7ED  # spike FIXTURE stream (validation/consts.jl:79)
    const RATIO_PAIR_SEED     = 0x0000_0000_004A_7107  # shipped ratio pairing (src/amortized/train_ratio.jl:60)
    const CORPUS_MASTER_SEED  = 0x0000_0000_00C0_5EED  # provenance-manifest stream (manifest.csv header)
    const P11_DEV_SEED        = 0x0000_0000_0B11_DE71  # spike/validation/p11_consts.jl:144
    const P13_BURNED_DEV_SEEDS = ( ... 25 burned keys ... )
    const P13_REPO_SALTS = ( ... 8 salts ... )   # a repeated (seed,salt) pair is a repeated Philox key

    "The forbidden-seed list in declaration order, so an accidental duplicate stays visible."
    _p13_forbidden_list() = vcat(
        UInt64[UInt64(0), UInt64(DEFAULT_MASTER_SEED), ... UInt64(P11_DEV_SEED)],
        collect(UInt64, P13_BURNED_DEV_SEEDS),
        # DERIVED, NEVER TRUSTED FROM A COMMENT: both frozen gate families are RECOMPUTED by
        # the isolated read-only `_GC` module above and forbidden from the recomputed values.
        collect(UInt64, values(_GC.PROD_SEED)),
        collect(UInt64, values(_GC.PROD_SEED_V2)),
    )
    _p13_forbidden() = Set{UInt64}(_p13_forbidden_list())

    const P13_DEV_SEED = 0x0000_0000_0B13_DE71  # "0B13" = phase 13, "DE71" = DEV-1 (mirrors P11)
    const P13_FIX_SEED = 0x0000_0000_0B13_F1F7  # FIXTURE stream ... fixtures must NEVER consume
                                                # (and so never pre-observe) a reported stream.
    const P13_SALT     = 0x2545_F491_4F6C_DD1D  # xorshift64* multiplier; not any repo salt

    "The Phase-13 REPORTED RNG: `Philox4x` keyed by `(P13_DEV_SEED xor P13_SALT, counter)`."
    p13_rng(counter::Integer = 0) =
        Philox4x(UInt64, (UInt64(P13_DEV_SEED) ⊻ P13_SALT, UInt64(counter)))
    p13_fix_rng(counter::Integer = 0) =
        Philox4x(UInt64, (UInt64(P13_FIX_SEED) ⊻ P13_SALT, UInt64(counter)))

    const P13_TAU_COUNTER        = 1   # ... one per activity ...
    const P13_FIXTURE_COUNTER    = 99  # FIXTURES ONLY -- never a reported counter
```

The isolated-module read of the frozen gate constants sits **outside** the guard block, and the file
says why (`:91-100`): *"Julia rejects a `module` expression that is not at top level."* Note also the
path-depth reminder repeated in every `spike/p13/` file: **from `spike/p13/` the repo root is TWO
levels up, not three.**

**(d) The seed-disjointness assertion idiom — executable self-checks at the foot of the block**
(`:729-745`), verbatim:

```julia
    # D-01: both fresh streams are disjoint from every reserved, burned and RECOMPUTED seed.
    @assert !(UInt64(P13_DEV_SEED) in _p13_forbidden()) "P13_DEV_SEED collides with a forbidden seed"
    @assert !(UInt64(P13_FIX_SEED) in _p13_forbidden()) "P13_FIX_SEED collides with a forbidden seed"
    @assert UInt64(P13_DEV_SEED) != UInt64(P13_FIX_SEED) "the fixture stream must differ from the reported stream"
    @assert !(P13_SALT in P13_REPO_SALTS) "P13_SALT reuses an existing repository salt"
    # A duplicate inside the forbidden list would silently shrink the set; catch that.
    @assert length(_p13_forbidden()) == length(_p13_forbidden_list()) "duplicate entry in the forbidden-seed list"
    # The recompute path actually ran (a comment naming a derived seed is not evidence).
    @assert any(v -> UInt64(v) in _p13_forbidden(), values(_GC.PROD_SEED))
    @assert any(v -> UInt64(v) in _p13_forbidden(), values(_GC.PROD_SEED_V2))
    # Fixtures must never ride a counter a reported number rides on.
    @assert P13_FIXTURE_COUNTER ∉ (P13_TAU_COUNTER, P13_DATAGEN_COUNTER, P13_GATE_COUNTER,
                                   P13_ALPHA_COUNTER, P13_CONTINUITY_COUNTER)
```

The **key-word cross-product** form VALIDATION's "Seeds" row demands already exists at
`spike/p13/run_three_way_gate.jl:159-165` — copy it:

```julia
# The key words that must differ, asserted rather than assumed. "Obviously distinct" is how two
# lanes end up sharing a Philox key ...
@assert (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_GATE_COUNTER)) !=
        (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_DATAGEN_COUNTER)) "the gate key word collides with the training-pool key word"
@assert (UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(P13_GATE_COUNTER)) !=
        (UInt64(P13_FIX_SEED) ⊻ P13_SALT ⊻ UInt64(P13_FIXTURE_COUNTER)) "the gate key word collides with the fixture family"
```

**(e) The iteration-allowance pattern** (`:353-365`), verbatim — Wave 0 lists
`P14_ITERATION_ALLOWANCE`, and this is its exact shape (an allowance **plus** a named,
pre-declared trigger, plus the sentence that forbids self-amendment):

```julia
    # ONE documented iteration is authorised, and it is authorised HERE, before any result
    # exists. A SECOND iteration is NOT authorised by this file and cannot be authorised by
    # amending this file. Phase 7 was amended TWICE after seeing results and the credibility
    # cost of that is the reason this constant exists.
    const P13_ITERATION_ALLOWANCE = 1
    const P13_ITERATION_TRIGGER = """
    THE ONE PRE-DECLARED CONDITION for spending P13_ITERATION_ALLOWANCE: the D-12 evaluation
    shows the exclusion-class per-class AUC below P13_AUC_FLOOR_EXCLUSION SPECIFICALLY AT HIGH
    |rho| ... The allowance is then spent on switching to P13_STRATIFICATION_FALLBACK (D-07-ii)
    and re-running ONCE. It is NEVER spent on relaxing a threshold after training. A second
    iteration is not authorised by this file.
    """
```

`P13_ITERATION_TRIGGER` is then **written into the artifact** by the runner
(`run_three_way_gate.jl:728-729`: `iteration_trigger_fired = trigger_fired,
iteration_trigger_text = P13_ITERATION_TRIGGER`) and asserted non-trivial by the test
(`test_p13_consts.jl:151-154`). Copy all three sites.

**(f) The Tier-2 append + τ's freeze ordering** (`:790-872`). This is the exact structure D-07 is
modelled on. The relocated (never deleted) Tier-1 self-check, verbatim `:800-803`:

```julia
# It sits ABOVE the Tier-2 block (outside the Tier-1 guard, so it also runs on a re-include)
# exactly as the Tier-1 header requires.
@assert (!isdefined(@__MODULE__, :P13_TAU) ||
         isdefined(@__MODULE__, :P13_TAU_PROBE_ARTIFACT)) "P13_TAU exists without its probe provenance: tau was set outside the appended Tier-2 block"
```

The freeze-ordering record, verbatim `:842-846` (**the lines the brief asked for**):

```julia
# THIS BLOCK WAS COMMITTED SEPARATELY, AFTER THE ARTIFACT COMMIT
# 581602a73a3f180f490f1c421c43028aa7cf9d47 that carries the runner and the curve. That
# separation is the auditable evidence that tau was MEASURED BEFORE IT WAS LOCKED: the
# measurement exists in git history at a commit that contains no Tier-2 value, so no reader
# has to take the ordering on trust.
```

and the Tier-2 block itself, `:853-872`:

```julia
if !isdefined(@__MODULE__, :P13_TAU)
    const P13_TAU = 0.15                    # first grid delta clearing the bar; A(tau) = 0.916356
    const P13_TAU_MEASURED_AUC = 0.91635625 # the realized A(tau) = max(A_neg, A_pos) at tau
    const P13_TAU_PROBE_SHA = "cecb69c58770a5b45bec8151c657d9024b999f80"
    const P13_TAU_SIMULATOR_SHA = "78dc37f517ad1b8e1e71afa963691f7ddefe5f63"
    const P13_TAU_PROBE_ARTIFACT = "spike/p13/tau_probe_report.jld2"
    const P13_TAU_REFERENCE_LAMBDA = 3.0

    # Executable, not decorative: tau must be a grid point (never an interpolated or rounded
    # value) and must actually have cleared the frozen bar.
    @assert P13_TAU in P13_TAU_DELTA_GRID "P13_TAU is not a point on the frozen pre-registered grid"
    @assert P13_TAU_MEASURED_AUC >= P13_TAU_AUC "P13_TAU_MEASURED_AUC does not clear the frozen P13_TAU_AUC bar"
end
```

**(g) The fail-loudly Tier-2 accessor** (`:266-282`) — D-07 wants τ *loaded*, so Phase 14's analogue
is a loader that errors rather than defaulting. The idiom:

```julia
    """
        p13_tau() -> Float64

    The Tier-2 tripwire. `P13_TAU` is DELIBERATELY ABSENT from the Tier-1 block: D-06 requires
    it to be MEASURED, not chosen. Calling this before the probe has run fails loudly.
    """
    p13_tau() = isdefined(@__MODULE__, :P13_TAU) ? P13_TAU : error("""
        P13_TAU is Tier-2 and does not exist yet. ...
        """)
```

**Phase-14 note:** `P13_TAU` is already in scope wherever `p13/consts.jl` is included, so a Phase-14
τ loader can cross-check the loaded artifact value against it for free (RESEARCH §D.3 recommends
treating the **net artifact** as source of truth and cross-asserting the other three).

---

### 3.2 `spike/p14/result.jl` — analog `spike/p13/result.jl`

**(a) The load order + read-only reach into `src/`** (`:63-84`), verbatim:

```julia
# ORDER MATTERS (spike/contract.jl:34-39, the documented load order): StatsBase + Statistics must
# be in scope before src/colocalization.jl is reached transitively, and Images before
# src/LoadImages.jl (its convenience constructor calls Images.otsu_threshold). ...
using StatsBase      # corspearman/corkendall for the transitively-reached correlation() Dict
using Statistics     # mean, quantile
using Images         # otsu_threshold (transitive src/LoadImages.jl requirement)

# --- Guarded includes, in dependency order --------------------------------------------------
# READ-ONLY reach into src/ for the supertype, `_iface_error`, `OODVerdict` and `CalibrationMeta`.
# NOTE THE PATH DEPTH: from spike/p13/ the repo root is TWO levels up, not three.
isdefined(@__MODULE__, :AbstractColocResult) ||
    include(joinpath(@__DIR__, "..", "..", "src", "results.jl"))
isdefined(@__MODULE__, :P13_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :_bin_calibration) ||
    include(joinpath(@__DIR__, "..", "validation", "sbc.jl"))
```

**(b) The guard-block split — a Julia-1.12 docstring trap, already paid for once** (`:86-93`),
verbatim. **A Phase-14 result type must copy this split or its docstrings silently vanish:**

```julia
# THE GUARD BLOCK COVERS ONLY THE `const`s AND THE `struct`, NOT THE FUNCTIONS, AND THAT SPLIT IS
# DELIBERATE. Julia 1.12 DROPS every docstring written inside an `if ... end` block -- the parser
# emits the `Core.@doc` call but the docsystem never registers it, so `@doc f` returns `nothing`
# for a function documented inside a guard (verified on 1.12.6). ... Method redefinition under a
# re-include is silent and harmless in Julia; only `const` and `struct` redefinition would warn or
# throw, and those are what the guard protects.
if !isdefined(@__MODULE__, :ThreeHypothesisColocResult)
```
…and the matching close at `:166`: `end # if !isdefined(...) -- consts + struct only`.

**(c) The inner-constructor validation idiom** (`:142-164`), verbatim — the brief's item 3:

```julia
struct ThreeHypothesisColocResult <: AbstractColocResult
    grid             :: Int
    posterior        :: Matrix{Float64}
    log_bf_vs_random :: ThreeWayLogBF
    ood              :: OODVerdict
    calibration      :: CalibrationMeta
    meta             :: NamedTuple

    function ThreeHypothesisColocResult(grid::Integer, posterior::AbstractMatrix{<:Real},
                                        log_bf_vs_random::NamedTuple, ood::OODVerdict,
                                        calibration::CalibrationMeta, meta::NamedTuple)
        lbf = _p13_as_three_way_logbf(log_bf_vs_random)
        # STRUCTURAL, NOT MERELY CONVENTIONAL. D-08 emits evidence AGAINST the random reference,
        # so the reference entry is the reference scored against itself: identically zero, for
        # every input, by construction. A non-zero value here means the caller built the triple
        # against some other reference (or fabricated a third measurement), and every downstream
        # comparison would silently be on a different scale.
        lbf.random === 0.0 || throw(ArgumentError(
            "ThreeHypothesisColocResult (D-08): log_bf_vs_random.random must be exactly 0.0 " *
            "-- the random class is the REFERENCE, scored against itself. Got $(lbf.random)."))
        return new(Int(grid), Matrix{Float64}(posterior), lbf, ood, calibration, meta)
    end
end
```

**(d) The key-order normalizer that makes composition safe** (`:179-184`) — RESEARCH §D.1 flags the
key-order trap (`three_way_log_bf` returns `(coloc, exclusion, random)`; `ThreeWayLogBF` stores
`(coloc, random, exclusion)`). **Always go through this; never build the NamedTuple positionally:**

```julia
function _p13_as_three_way_logbf(nt::NamedTuple)
    Set(keys(nt)) == Set(P13_LOG_BF_KEYS) || throw(ArgumentError(
        "ThreeHypothesisColocResult (D-08): the evidence triple must carry exactly the keys " *
        "$(P13_LOG_BF_KEYS); got $(keys(nt))."))
    return ThreeWayLogBF((Float64(nt.coloc), Float64(nt.random), Float64(nt.exclusion)))
end
```

with `const P13_LOG_BF_KEYS = (:coloc, :random, :exclusion)` (`:103`) and
`const ThreeWayLogBF = NamedTuple{P13_LOG_BF_KEYS, Tuple{Float64,Float64,Float64}}` (`:115`).
RESEARCH §D.2 **recommends composition (option a)**: carry a `ThreeHypothesisColocResult` as a
*field*, so the `random === 0.0` invariant is enforced once. That is compatible with everything
above and requires no re-statement.

**(e) The lossy-accessor documentation precedent** (`:204-220`) — this is the model for
`is_ood(r) = r.ood_state === :fired`:

```julia
"""
    bayes_factor(r::ThreeHypothesisColocResult) -> Float64

Return `r.log_bf_vs_random.coloc`, i.e. `log BF(coloc : random)`.

THIS CHOICE IS DOCUMENTED RATHER THAN SILENT, because a single-number accessor on a three-way
result is inherently lossy and the reader deserves to know which number they got. ...
A caller who wants THE THREE-WAY EVIDENCE must use [`log_bf_vs_random`](@ref) ...
"""
bayes_factor(r::ThreeHypothesisColocResult) = r.log_bf_vs_random.coloc
```

**(f) The do-not-fabricate accessor** (`:240-243`) — Phase 14's `delta_rho` must do exactly this:

```julia
function delta_rho(r::ThreeHypothesisColocResult)
    haskey(r.meta, :delta_rho_draws) || _iface_error(r, :delta_rho)
    return r.meta.delta_rho_draws
end
```
(docstring: *"DO NOT FABRICATE A DELTA RHO … Letting the documented interface error fire is the
honest answer"*, with D-05's ρ_s = 0.8 / ρ_c = 0.9 counter-example.)

**(g) `p13_calibration_meta`** (`spike/p13/result.jl:337-352`) — the required-keys-merged-over-caller
idiom, if the Phase-14 result carries calibration:

```julia
function p13_calibration_meta(cal::CalibrationResult; grid::Integer, auc::Real,
                              gate::NamedTuple = NamedTuple())
    empty_bins = count(iszero, cal.bin_counts)
    required = (statistic              = P13_GATE_STATISTIC,
                ece_green              = Float64(P13_ECE_GREEN),
                ece_yellow             = Float64(P13_ECE_YELLOW),
                n_bins                 = Int(P13_ECE_NBINS),
                mce_reported_not_gated = true,
                empty_bins             = empty_bins,
                auc                    = Float64(auc),
                vacuous                = vacuous_pass(cal.ece, auc),
                verdict                = p13_traffic_light(cal.ece))
    return CalibrationMeta(cal.bin_midpoints, cal.predicted_rate, cal.observed_rate,
                           cal.bin_counts, cal.ece, cal.mce, Int(grid),
                           merge(gate, required))
end
```
Two contracts to carry: **`auc` is a REQUIRED keyword, not an option** (the vacuous-pass guard is
only real if it cannot be omitted), and **caller `gate` is merged UNDER `required`** so provenance is
preserved but required keys can never be overwritten. Note the **empty-bin trap** documented at
`:321-329`: `_bin_calibration` scores an empty bin as `predicted_rate = midpoint, observed_rate =
0.0`, which is zero ECE weight but **full MCE weight** — so ECE gates and MCE is reported.

---

### 3.3 The `AbstractColocResult` accessor contract — `src/results.jl:35-79` (verbatim)

```julia
abstract type AbstractColocResult end     # :43

_iface_error(r::AbstractColocResult, f::Symbol) = error(                       # :47-49
    "`$(f)` is not implemented for $(typeof(r)). Every `AbstractColocResult` subtype must " *
    "implement the accessor interface: delta_rho, bayes_factor, is_ood, posterior_draws.")

delta_rho(r::AbstractColocResult)       = _iface_error(r, :delta_rho)        # :57
bayes_factor(r::AbstractColocResult)    = _iface_error(r, :bayes_factor)     # :65
is_ood(r::AbstractColocResult)          = _iface_error(r, :is_ood)           # :72
posterior_draws(r::AbstractColocResult) = _iface_error(r, :posterior_draws)  # :79
```

**Exactly four accessors are required** — `delta_rho`, `bayes_factor`, `is_ood`, `posterior_draws`.
Supporting value types in the same file: `OODVerdict(score::Float64, flag::Bool,
per_channel::NamedTuple)` at `:93-97`, and the eight-field `CalibrationMeta` at `:117-126`
(`bin_midpoints, predicted_rate, observed_rate, bin_counts, ece, mce, grid, gate`).

New Phase-14 accessors layer on top and are **never bolted-on fields** (Phase-7 D-02):
`decision`, `abstain_reason`, `ood_state`, `conformal_set`, `null_posterior`, `cross_method`,
`local_map` (RESEARCH §D.6).

**The subtype-relation assertions to copy** (`spike/test/test_p13_result.jl:84-94`):

```julia
@testset "is a new AbstractColocResult subtype (Phase-7 D-02)" begin
    @test ThreeHypothesisColocResult <: AbstractColocResult
    @test isabstracttype(AbstractColocResult)
    @test isconcretetype(ThreeHypothesisColocResult)
    @test _res() isa AbstractColocResult
    # D-02's other half: the new variant does NOT bolt fields onto the shipped struct.
    @test fieldnames(AmortizedColocResult) == (:grid, :posterior, :delta_rho_draws,
                                               :log_bayes_factor, :ood, :calibration, :meta)
end
```

---

### 3.4 Reusable helpers — **no dependency may be added, so these MUST be reused**

| Helper | Exact path : line | Signature | Notes |
|---|---|---|---|
| `roc_auc` | `C:\...\ProteinCoLoc\spike\validation\ood.jl:319` | `roc_auc(neg_scores::AbstractVector, pos_scores::AbstractVector) -> (fpr, tpr, auc)` | Tie-aware Mann–Whitney U. Returns `(Float64[0,1], Float64[0,1], NaN)` when either arm is empty. **This is the SC3-c "in-repo `roc_auc`".** |
| `roc_auc` (shipped twin) | `C:\...\ProteinCoLoc\src\amortized\ood.jl:277` | same | Identical implementation. The **spike-lane** copy is the one the p13 runner includes. |
| AUC-only wrapper | `C:\...\ProteinCoLoc\spike\p13\run_three_way_gate.jl:215` | `p13_gate_auc(neg_scores, pos_scores) = roc_auc(neg_scores, pos_scores)[3]` | Copy this one-liner rather than re-indexing everywhere. |
| `id_threshold` | `C:\...\ProteinCoLoc\spike\validation\ood.jl:344` | `id_threshold(id_scores::AbstractVector; q::Real = OOD_ID_QUANTILE) -> Float64` | `= quantile(id_scores, q)`. `OOD_ID_QUANTILE = 0.95`. **The pre-registered operating point — NOT Youden-J.** |
| Spearman | `StatsBase.corspearman` | `corspearman(x, y) -> Float64` | Used at `spike\p13\run_three_way_gate.jl:601`; `spike\validation\run_p11_probe.jl:59` comments *"corspearman (never hand-rolled, PATTERNS S6)"*. **StatsBase is in the frozen sixteen. Do not hand-roll SC3-b's ρ.** |
| `_bin_calibration` | `C:\...\ProteinCoLoc\spike\validation\sbc.jl:82` | `_bin_calibration(posterior_probs::AbstractVector{<:Real}, empirical_positive::AbstractVector{Bool}; n_bins::Int = 10) -> CalibrationResult` | Equal-width bins over [0,1], right endpoint included in the last bin. **Shared gate-lineage code — NEVER hand-modified.** |
| `CalibrationResult` | `C:\...\ProteinCoLoc\spike\validation\sbc.jl:65-72` | 6 fields: `bin_midpoints, predicted_rate, observed_rate, bin_counts, ece, mce` | Same order as `CalibrationMeta`'s first six. |
| `sbc_traffic_light` | `C:\...\ProteinCoLoc\spike\validation\sbc.jl:139` | `sbc_traffic_light(ece::Real) -> Symbol` | Phase 13 **re-declared its own** band locally (`p13_traffic_light`, `result.jl:276`) rather than importing this — *"every pre-registration in this project is SELF-CONTAINED"*. Follow the re-declare pattern if Phase 14 needs a band. |
| `ghat` | `C:\...\ProteinCoLoc\spike\simulator\ghat.jl:110` | `ghat(μ::Real) -> Float64` | Frozen monotone clamped piecewise-linear μ → ρ_true; knots at `:93`/`:95`, clamp bounds `GHAT_MU_MIN=-0.67976` / `GHAT_MU_MAX=0.847149` at `:99-100`. **SC2-d requires every `patch_correlation(` compared to τ to be wrapped in `ghat(`.** |
| `three_way_label` | `C:\...\ProteinCoLoc\spike\p13\labels.jl:121` | `three_way_label(rho_sample::Real, rho_control::Real; tau::Real = p13_tau()) -> ThreeWayClass` | Throws `ArgumentError` on `tau <= 0`. **Loaded, never reimplemented.** |
| `class_masses` | `C:\...\ProteinCoLoc\spike\p13\labels.jl:301` | `class_masses(classes::AbstractVector{ThreeWayClass})` | The re-derivation surface for VALIDATION's "Prior re-derivation" row. |
| `three_way_log_bf` / `three_way_probs` / `encode_conditioned_pair` | `spike\p13\net.jl:349` / `:372` / `:429` | see RESEARCH §D.1 | `encode_conditioned_pair` self-asserts its built length against `three_way_input_dim(G, length(c))` (`net.jl:434-436`) — **never write 321 as a literal.** |
| `p13_real_id_reference` | `C:\...\ProteinCoLoc\spike\p13\run_p13_realimage.jl:412` | `p13_real_id_reference(pool_dir) -> NamedTuple` | The OOD-null construction for the spike lane: reads `shard_*.jld2` from `h.meta.pool_dir`, stacks **both** `p.Zs` and `p.Zc`, fits the density null, takes the copied shipped ID quantile. **Reuse verbatim** (RESEARCH §D.4). Note its `bad(note)` early-return idiom — every failure path returns `(available = false, …, note = "why")` instead of throwing. |

---

### 3.5 `spike/p14/fuse.jl` — analogs `src/amortized/ood.jl:363` + `src/amortized/local_map.jl:93`

**(a) The shipped OR-fusion tail, verbatim** (`src/amortized/ood.jl:385-397`) — this is where
`flag == false` becomes two-valued, and it is the mechanic D-06 exists to handle:

```julia
    zref  = haskey(ood_nulls, :zref) ? ood_nulls.zref : nothing
    fused = -Inf
    per   = NamedTuple()
    for (name, s) in channels
        z     = zref === nothing ? float(s) : _robust_z(zref, name, s)
        fused = max(fused, z)
        per   = merge(per, NamedTuple{(name,)}((s,)))
    end

    thr  = haskey(ood_nulls, :thr) ? ood_nulls.thr : nothing
    flag = thr === nothing ? false : fused > thr
    return OODVerdict(fused, flag, per)
```

**(b) The `nothing`-returning three-valued resolver to mirror** (`src/amortized/local_map.jl:93-108`),
verbatim. Its acceptance predicate is the one a Phase-14 resolver must match exactly:

```julia
function _recorded_ood_threshold(artifact_dir, grid)
    path = joinpath(artifact_dir, "gate_report_$(Int(grid)).jld2")
    isfile(path) || return nothing
    rep = try
        d = JLD2.load(path)
        get(d, "report", nothing)
    catch e
        e isa InterruptException && rethrow()
        nothing
    end
    (rep isa NamedTuple && haskey(rep, :ood)) || return nothing
    o = rep.ood
    (o isa NamedTuple && haskey(o, :id_threshold)) || return nothing
    t = o.id_threshold
    return (t isa Real && isfinite(t)) ? float(t) : nothing
end
```
Note `e isa InterruptException && rethrow()` inside the `catch` — the house idiom; a bare `catch`
that swallows Ctrl-C is the bug it prevents. Docstring (`:89-91`): *"Returns `nothing` — never a
substitute value … A missing operating point must leave the flag HONESTLY INERT rather than
silently default to something arbitrary."*

**(c) The NOT-CHECKED semantics, quoted** (`src/amortized/local_map.jl:68-71`):

> `meta::NamedTuple`: provenance (`N`, `ood_score`, `n_sentinel`, `tile_size`, `channels`,
> `sentinel`, `ood_thr`). `ood_thr` is the recorded pre-registered operating point the flags were
> decided against; **`nothing` means none was available, so a `false` flag on a non-sentinel tile
> means "NOT CHECKED", not "in distribution".**

and `:116-119`: *"Without one, `ood_verdict` short-circuits to `flag = false` and the per-tile flag
is VACUOUS: it then fires only structurally (sentinel tiles), never as a detector, which silently
understates risk."*

**(d) The branch-per-decision, cite-the-decision shape** — `three_way_label`
(`spike\p13\labels.jl:121-131`) is the in-repo model for `p14_fuse`: validate the input first, then
early-return one branch per rule, with the counter-example that makes the rule load-bearing recorded
in the docstring (*"`three_way_label(0.8, 0.9; tau = 0.10) === RANDOM  # NEVER EXCLUSION`"*).
RESEARCH §Code Examples already gives `p14_fuse` in exactly this shape.

**(e) `LocalColocMap` has no uncertainty field** (`src/amortized/local_map.jl:73-79`:
`grid, tiles, delta_rho, ood_flag, meta`) — the direct evidence for D-04. Carry it as a labelled
uncontrolled-display field; it never enters the FDR arithmetic.

---

### 3.6 `spike/p14/fdr.jl`, `conformal.jl`, `posterior.jl` — no direct logic analog

There is **no existing FDR or conformal code in this repo**. RESEARCH §Code Examples supplies all
three implementations already; the *pattern* to copy is the house numeric-helper shape, of which
`roc_auc` (`spike\validation\ood.jl:309-335`) is the canonical instance:

- a docstring that states the returned tuple, the convention (*"HIGHER score ⇒ more OOD"*), and
  **why no package was added** (*"No dependency (~15 lines)"* — mirrors `14-RESEARCH`'s
  *"The rule is `sort` + `cumsum` + `findlast`. Three lines."*);
- an explicit degenerate-input branch **before** the arithmetic (`(P == 0 || N == 0) && return …`);
- pure `Base`/`Statistics` operations, no allocation of a package-level abstraction.

`ghat` (`spike\simulator\ghat.jl:110-121`) is the second instance and the closer analog for
`p14_conformal_quantile`: a frozen, clamped, order-statistic-style lookup whose docstring names the
clamping behaviour explicitly (*"μ beyond [GHAT_MU_MIN, GHAT_MU_MAX] clamps ρ_true to the nearest
swept endpoint"*). RESEARCH's `p14_conformal_quantile` correspondingly **throws** rather than clamps
when `alpha < 1/(n+1)` — and SC1-e asserts that throw at n = 30, α = 0.03.

For `posterior.jl`, the class-key discipline comes straight from `spike\p13\result.jl:98-103`:

```julia
# The KEY ORDER of the stored NamedTuple is (:coloc, :random, :exclusion). The sketch's positional
# `NTuple{3,Float64}` was documented as `{null, coloc, anti-coloc}` -- REFERENCE CLASS FIRST --
# which no reader would guess ... A NamedTuple removes the hazard outright: every read site is by
# NAME, so a reordering cannot change a number.
const P13_LOG_BF_KEYS = (:coloc, :random, :exclusion)
```
This is why VALIDATION's "Class order" row demands an **unequal** (0.7 / 0.2 / 0.1) fixture.

---

### 3.7 `spike/p14/decide.jl` — analogs `src/amortized/api.jl:107` + the spike composition lane

The **src-shaped signature** to imitate (`src/amortized/api.jl:107-110`):

```julia
function colocalization_amortized(img::MultiChannelImage, control::MultiChannelImage,
                                  channels::AbstractVector{<:Integer};
                                  num_patches::Integer = 8, N::Integer = 2000,
                                  use_gpu::Bool = false)
```

Copy four things from it:

1. **Validation-first (V5)**, `:111-120` — arity, channel range, positivity, each with a message
   naming the function: `throw(ArgumentError("colocalization_amortized: \`channels\` must select
   exactly 2 channels, got $(length(channels)).")`
2. **A `# Named limits (non-negotiable honesty — carried into the manuscript)` docstring section**
   (`:96-105`). Phase 14's equivalent must carry D-03's simulator-derived conformal guarantee, D-04's
   uncontrolled tiles, and D-01's "NOT shipped in v2.0".
3. **The result is constructed once, at the end, with a provenance `meta` NamedTuple** (`:145-147`).
4. **`# Returns` naming the accessors** the caller should read it through.

Batch shape: `decide_coloc` is batch-level (FDR is a batch quantity) → `P14BatchDecision` with
`fdr_scope = :decided_subset_only` (RESEARCH §D.6).

The net-loading preamble (RESEARCH §D.1) is the guarded-include chain:

```julia
isdefined(@__MODULE__, :ThreeWayEvidenceNet) || include(joinpath(@__DIR__, "..", "p13", "net.jl"))
isdefined(@__MODULE__, :three_way_label)     || include(joinpath(@__DIR__, "..", "p13", "labels.jl"))
isdefined(@__MODULE__, :ThreeHypothesisColocResult) ||
    include(joinpath(@__DIR__, "..", "p13", "result.jl"))
isdefined(@__MODULE__, :load_p13_basis)      || include(joinpath(@__DIR__, "..", "p13", "preconditions.jl"))
h = load_three_way(joinpath(@__DIR__, "..", "p13", "three_way_net.jld2"))
```

---

### 3.8 `spike/p14/run_p14_*.jl` — analog `spike/p13/run_three_way_gate.jl`

**(a) The anti-snooping banner + the order-of-operations rule** (`:21-62`), the load-bearing part
verbatim:

```julia
# EVERY THRESHOLD THIS RUNNER SCORES AGAINST WAS COMMITTED TO `spike/p13/consts.jl` BEFORE THE
# NET EXISTED. ... Not one of them is a literal here, and the test-side acceptance criteria
# assert that, so a value cannot be "temporarily" inlined and left behind.
#
# AN HONEST SHORTFALL IS A PHASE-13 FINDING, NEVER A LICENCE TO ACT. ...
#
# THE ORDER OF OPERATIONS IS THE POINT: read the thresholds, compute the statistics, PERSIST THE
# ARTIFACT, print the headline, then assert. The artifact is written before anything can throw,
# so a FAILING gate still leaves a complete, self-describing report on disk.
#
# THIS FILE IS A RUNNER, NOT A LIBRARY AND NOT A TEST. It composes surfaces that already exist
# and are already tested. Including it does nothing; it must be run deliberately:
#
#     julia --project=spike -t auto spike/p13/run_three_way_gate.jl
```
The figure is likewise rendered **before** the gate (`:768-772`) for the same reason: *"the evidence
of an honest failure must survive the failure."*

**(b) The atomic `.jld2` save + integrity reopen idiom** (`:171-190`), verbatim:

```julia
function _p13_gate_save_report(path; kwargs...)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f
        @assert haskey(f, "auc_coloc") "_p13_gate_save_report: integrity check failed ($tmp)"
        @assert haskey(f, "confusion") "_p13_gate_save_report: integrity check failed ($tmp)"
    end
    mv(tmp, path; force = true)
    return path
end
```
The same idiom, for models, is `save_three_way` / `load_three_way` (`spike\p13\net.jl:459-514`) with
its **narrow deserialization surface** rule: store `Flux.state(net)` + architecture integers, never
the estimator object; `load_three_way` rebuilds the topology and errors on a schema mismatch.

**(c) The git blob sha helper** (`:200-206`) — points the artifact *into* the audit trail:

```julia
function p13_blob_sha(path)
    bytes = read(path)
    ctx   = SHA.SHA1_CTX()
    SHA.update!(ctx, Vector{UInt8}("blob $(length(bytes))\0"))
    SHA.update!(ctx, bytes)
    return bytes2hex(SHA.digest!(ctx))
end
```
plus `p13_consts_sha() = bytes2hex(SHA.sha256(read(joinpath(@__DIR__, "consts.jl"))))`
(`spike\p13\net.jl:524`).

**(d) The provenance-metadata block written into every artifact** (`:738-765`), verbatim structure —
**this is the shape a Phase-14 `.jld2` must copy**:

```julia
        # --- EVERY THRESHOLD THIS RUN WAS SCORED AGAINST, written INTO the artifact ---
        P13_GATE_M = P13_GATE_M,
        P13_AUC_FLOOR_COLOC = P13_AUC_FLOOR_COLOC,
        ...
        P13_ITERATION_ALLOWANCE = P13_ITERATION_ALLOWANCE,
        P13_TAU = P13_TAU, P13_TAU_REFERENCE_LAMBDA = P13_TAU_REFERENCE_LAMBDA,
        P13_CUT_VARIANT = P13_CUT_VARIANT,
        # --- provenance: the stream, the net, the frozen pre-registration ---
        master_seed = UInt64(P13_DEV_SEED), salt = UInt64(P13_SALT),
        gate_counter = Int(P13_GATE_COUNTER), datagen_counter = Int(P13_DATAGEN_COUNTER),
        net_path = net_path, net_head_log_odds = hlo,
        net_consts_sha = h.consts_sha, net_meta = h.meta,
        consts_sha = p13_consts_sha(), consts_git_blob_sha = p13_blob_sha(consts_path),
        phase11_net = P13_PHASE11_NET, phase11_sha = p13_phase11_sha(),
        grid = G, reported = reported,
        amendment = "13-SC2-AMENDMENT.md -- any report citing this result must cite the " *
                    "original ROADMAP SC1/SC2/SC3 alongside it",
        generated = string(Dates.now(Dates.UTC)) * "Z")
    verbose && println("      persisted (before any verdict) -> $report_path")
```
Note the `amendment =` key — **Phase 14 amends SC1 (D-02) and SC2 (D-05) and must carry the same
field.** Note also the per-item score columns (`:730-733`), written *"so a later plan need not
re-run the gate"*.

**(e) The run-time decoupling proof, step 0 of `main`** (`:270-280`), verbatim — every Phase-14
runner should open with it:

```julia
    # Asserted while the run happens, not only in review: a long reported run that started on a
    # clean tree and finished on a dirty one would be a repudiation hole (T-13-39).
    src_clean = success(Cmd(`git diff --quiet HEAD -- src`; dir = P13_GATE_REPO_ROOT))
    env_clean = success(Cmd(`git diff --quiet HEAD -- spike/Project.toml spike/Manifest.toml`;
                            dir = P13_GATE_REPO_ROOT))
    @assert src_clean "D-01 decoupling breach: `git diff --quiet HEAD -- src` FAILED ..."
    @assert env_clean "D-01 decoupling breach: spike/Project.toml or spike/Manifest.toml is modified ..."
```

**(f) The smoke-mode redirect** (`:264-268`) — a non-reported run never overwrites the reported
artifact:

```julia
    reported = (m == P13_GATE_M) && (min_per_head == P13_MIN_EVAL_PER_HEAD)
    if !reported
        report_path = replace(report_path, ".jld2" => "_smoke.jld2")
        fig_path    = replace(fig_path, ".png" => "_smoke.png")
    end
```

**(g) The locked-threshold banner** (`:296-355`) — the runner prints every frozen constant, its
source, and which are GATED vs REPORTED-NOT-GATED, before drawing anything. SC1-c ("REPORTED, NOT
GATED — no bar") and SC3's four judgement-call bars both need this treatment.

**(h) The guarded script entry point** (`net.jl:566` / `run_three_way_gate.jl:874`):
`if abspath(PROGRAM_FILE) == @__FILE__` — including a runner must do nothing.

---

### 3.9 `spike/test/test_p14_consts.jl` — analog `spike/test/test_p13_consts.jl`

**(a) The file layout, quoted from its own banner** (`:34-36`): *"Mirrors test_sbc.jl / test_bf.jl:
license header -> using Test -> guarded include of the unit under test -> source read ONCE at module
level -> one outer @testset. No training, no simulation, no figure call -- this gate is
milliseconds."*

**(b) The comment-stripped source read at module level** (`:44-50`), verbatim:

```julia
const P13_CONSTS_SRC = read(joinpath(@__DIR__, "..", "p13", "consts.jl"), String)
const P13_CONSTS_CODE = join(
    filter(l -> !startswith(strip(l), "#"), split(P13_CONSTS_SRC, '\n')), '\n')
```

**(c) The literal-value + disjointness testsets** (`:54-94`), verbatim:

```julia
    @testset "seeds are fresh and disjoint" begin
        # The locked literals. A changed seed is a changed experiment.
        @test P13_DEV_SEED == 0x0000_0000_0B13_DE71
        @test P13_FIX_SEED == 0x0000_0000_0B13_F1F7
        @test P13_SALT     == 0x2545_F491_4F6C_DD1D
        @test !(UInt64(P13_DEV_SEED) in _p13_forbidden())
        @test !(UInt64(P13_FIX_SEED) in _p13_forbidden())
        # Pairwise, by name, so a failure says WHICH stream was collided with.
        for reserved in (VAL_MASTER_SEED, NPE_MASTER_SEED, VAL_FIX_SEED, RATIO_PAIR_SEED,
                         CORPUS_MASTER_SEED, DEFAULT_MASTER_SEED, P11_DEV_SEED)
            @test UInt64(P13_DEV_SEED) != UInt64(reserved)
            @test UInt64(P13_FIX_SEED) != UInt64(reserved)
        end
        @test UInt64(P13_DEV_SEED) != UInt64(P13_FIX_SEED)
        @test !(P13_SALT in P13_REPO_SALTS)
        @test any(v -> UInt64(v) in _p13_forbidden(), values(_GC.PROD_SEED))
        @test any(v -> UInt64(v) in _p13_forbidden(), values(_GC.PROD_SEED_V2))
        @test length(_p13_forbidden()) == length(_p13_forbidden_list())
    end

    @testset "streams are counter-separated" begin
        # Distinct reserved counters must not share a sub-stream.
        @test rand(p13_rng(P13_TAU_COUNTER), UInt64) != rand(p13_rng(P13_DATAGEN_COUNTER), UInt64)
        # The fixture stream never overlaps the reported one, even at the same counter.
        @test rand(p13_fix_rng(P13_FIXTURE_COUNTER), UInt64) !=
              rand(p13_rng(P13_FIXTURE_COUNTER), UInt64)
        # Reproducibility: same construction, same draw.
        @test rand(p13_rng(P13_TAU_COUNTER), UInt64) == rand(p13_rng(P13_TAU_COUNTER), UInt64)
        @test P13_FIXTURE_COUNTER ∉ (P13_TAU_COUNTER, P13_DATAGEN_COUNTER, P13_GATE_COUNTER,
                                     P13_ALPHA_COUNTER, P13_CONTINUITY_COUNTER)
    end
```
The **draw-a-number** form (`rand(rng, UInt64) != rand(rng2, UInt64)`) is the executable proof of
stream separation; VALIDATION's "Seeds" row extends it to the **full Philox key-word cross-product**
against the P13 family — the key words are `P13_DEV_SEED ⊻ P13_SALT ⊻ c` for
`c ∈ {1,2,3,4,5,99}` and the bare `P13_DEV_SEED ⊻ P13_SALT`.

**(d) The Tier-2 provenance testset** (`:121-144`) — the direct analog for `test_p14_provenance.jl`:

```julia
    @testset "tau is measured and locked (Tier 2)" begin
        @test isdefined(@__MODULE__, :P13_TAU)
        @test P13_TAU in P13_TAU_DELTA_GRID
        @test P13_TAU_MEASURED_AUC >= P13_TAU_AUC
        # Both provenance shas are full 40-character hex object names ...
        @test occursin(r"^[0-9a-f]{40}$", P13_TAU_PROBE_SHA)
        @test occursin(r"^[0-9a-f]{40}$", P13_TAU_SIMULATOR_SHA)
        @test P13_TAU_REFERENCE_LAMBDA == P13_TAU_REFERENCE_LAMBDA_EXPECTED
        @test P13_TAU_PROBE_ARTIFACT == "spike/p13/tau_probe_report.jld2"
        @test p13_tau() === P13_TAU
    end
```

**(e) The iteration-allowance assertions** (`:150-154`):

```julia
        @test P13_ITERATION_ALLOWANCE == 1
        @test P13_ITERATION_TRIGGER isa AbstractString
        @test !isempty(strip(P13_ITERATION_TRIGGER))
        @test occursin("high", lowercase(P13_ITERATION_TRIGGER))
```

---

### 3.10 `spike/test/test_p14_provenance.jl` — the D-07 sha table

**The pinned-both-sides assertion, verbatim** (`spike\p13\run_three_way_gate.jl:361-370`) — this is
the exact form D-07's "assert the divergence" requires, **and its rejected widening is recorded so
an executor cannot re-introduce it**:

```julia
    # THE TRAIN-TIME PROVENANCE GUARD ... The two sides are now pinned to NAMED, DATED literals
    # instead of merely to each other, which is STRICTLY STRONGER: the old single equality would
    # have passed silently after a retrain plus an undisclosed edit. They are EXPECTED to differ
    # here ... and that divergence is asserted rather than tolerated. ... The widened `||` form is
    # REJECTED: it would accept any future drift silently.
    @assert h.consts_sha     == P13_CONSTS_SHA.pre_amendment  "the net was NOT trained under the frozen pre-registration: ..."
    @assert p13_consts_sha() == P13_CONSTS_SHA.post_amendment "spike/p13/consts.jl is at NEITHER sha256 this amendment authorises ...: a third value is drift, not an amendment"
```

with the two literals at `spike\p13\net.jl:556-559`:

```julia
const P13_CONSTS_SHA = (
    pre_amendment  = "100e97a37bb470bb7cb4fbdbd8e33f219d83adcf61fe398a90cf3ad3391f3fa6",
    post_amendment = "d6a6e63bf19bb13ad5c61eeaaf9fc18fab27c5aa7ba879bad5c493a5c71ae247",
)
```

and the rejection recorded in that docstring (`net.jl:549-554`): *"REJECTED, and recorded so a later
reader knows it was considered: the widened form `h.consts_sha == p13_consts_sha() || h.consts_sha
== PRE_AMENDMENT_SHA`. Once that disjunction exists it accepts ANY future drift silently."*

**⚠ THE TRAP.** `tau_probe_report.jld2` records a **third** `consts_sha256`
(`640ec90e…aca4cf`), matching neither literal — correct, because the probe ran before the Tier-2
append. VALIDATION's D-07 row therefore demands `probe["consts_sha256"] ∉ values(P13_CONSTS_SHA)` —
**assert the divergence**; widening it to a disjunction is itself the failure.

The `@testset` skeleton for the git half is `test_p13_result.jl:190-212`, verbatim:

```julia
        _has_git = Sys.which("git") !== nothing &&
                   (isdir(joinpath(P13_REPO_ROOT, ".git")) ||
                    isfile(joinpath(P13_REPO_ROOT, ".git")))
        if _has_git
            @test success(Cmd(`git diff --quiet HEAD -- src/results.jl`; dir = P13_REPO_ROOT))
            @test success(Cmd(`git diff --quiet HEAD -- src`; dir = P13_REPO_ROOT))
        else
            @info "git unavailable — skipping the executable src/ byte-equality assertion"
        end
```
(`isdir(...) || isfile(...)` because **inside a linked worktree `.git` is a FILE** —
`test_p12_decoupling.jl:64`.)

---

### 3.11 `spike/test/test_p14_decoupling.jl` — analog `spike/test/test_p12_decoupling.jl`

Beyond §0's four quoted testsets, copy these five mechanisms:

**(a) The self-scanning problem and its two answers** (`:34-45`), verbatim — the file matches its own
scan pattern, so a naive token ban is red by construction:

```julia
# THE SOURCE SCAN SCANS *ITSELF*, AND THAT IS WHY THREE NEEDLES ARE BUILT BY CONCATENATION.
#   * the ALLOWLIST (`_P12_CORPUS_METADATA_ALLOWED`) exempts the four files that legitimately
#     mention corpus METADATA ...;
#   * the three needles that would still match this file even under the allowlist
#     (`open_` + `sealed_holdout`, and the two spellings of the sealed IMAGE directory) are
#     assembled at run time from fragments, so the contiguous literal never appears in this
#     source at all.
const _SEALED_OPENER            = "open_" * "sealed_holdout"
const _SEALED_IMAGE_DIR_LITERAL = "cor" * "pus/data"
const _SEALED_IMAGE_DIR_SYMBOL  = "CORPUS_" * "DATA_DIR"
```

**(b) The comment-stripper and the grow-with-the-phase source enumeration** (`:71-106`):

```julia
_strip_comment_lines(src::AbstractString) =
    join(filter(l -> !startswith(strip(l), "#"), split(src, '\n')), '\n')

const _P12_SOURCE_PATTERNS = (
    ("spike/validation", r"^p12_.*\.jl$"), ..., ("spike/p12", r"^.*\.jl$"),
    ("spike/test",       r"^test_p12_.*\.jl$"), ...)

"Repo-relative, forward-slashed paths of every Phase-12 source that EXISTS right now."
function _p12_sources()
    out = String[]
    for (rel, pat) in _P12_SOURCE_PATTERNS
        d = joinpath(P12_REPO_ROOT, rel)
        isdir(d) || continue
        for f in readdir(d)
            occursin(pat, f) || continue
            p = rel * "/" * f
            p in out || push!(out, p)
        end
    end
    return sort(out)
end
```
*"Files that do not exist yet are simply absent from the enumeration: … the scan grows with the phase
instead of asserting against a frozen list that goes stale."* Phase-14 patterns:
`("spike/p14", r"^.*\.jl$")` and `("spike/test", r"^test_p14_.*\.jl$")`.

**(c) The allowlist declared as data, with the "adding an entry is a DECISION" note** (`:109-118`):

```julia
const _P12_CORPUS_METADATA_ALLOWED = ("spike/validation/p12_consts.jl", ...)
const _P12_CORPUS_TOKENS = ("corpus", "sealed_holdout", "manifest.csv", "CORPUS_MASTER_SEED")
const _IMAGE_LOADERS = ("load(", "FileIO.load", "build_mci")
```

**(d) The four sealed-holdout checks** (`:229-296`), in order: (a) no image load resolves onto the
corpus in ANY source, allowlisted or not; (b) non-allowlisted sources may not mention it at all
(comment-stripped); (c) allowlisted files get the narrower property (metadata yes, sealed image
directory never); (d) the **positive** `test/test_images` assertion so the three negatives cannot be
satisfied vacuously. The self-consistency check at `:233-234` is worth copying literally:

```julia
        # The allowlist must describe files that really are in the scanned surface, or an
        # exemption could silently protect nothing (a typo, or a renamed file).
        @test all(a -> !isfile(joinpath(P12_REPO_ROOT, a)) || a in sources,
                  _P12_CORPUS_METADATA_ALLOWED)
```
and the failures collect **offending lines, not booleans** (`:238-244`: `push!(load_hits, "$p: $ln")`
then `@test load_hits == String[]`) — *"so a failure shows the code, not just a boolean."*

**(e) The cache-root separation testset** (`:321-349`) — the Runtime State Inventory says Phase-14
caches go under a new `cache/p14/` root:

```julia
        p11_root = normpath(joinpath(P12_REPO_ROOT, "spike", "data", "cache", "p11"))
        p12_root = normpath(joinpath(P12_REPO_ROOT, "spike", "data", "cache", "p12"))
        @test p11_root != p12_root
        @test !startswith(p11_root, p12_root * Base.Filesystem.path_separator)
        ...
        _p12_gen = joinpath(P12_REPO_ROOT, "spike", "data", "p12_generate.jl")
        if isfile(_p12_gen)
            _p12_gen_body = _strip_comment_lines(read(_p12_gen, String))
            @test occursin("joinpath(@__DIR__, \"cache\", \"p12\")", _p12_gen_body)
            @test !occursin("\"p11\"", _p12_gen_body)
        else
            @info "spike/data/p12_generate.jl does not exist yet (12-09) — ..."
        end
```
The `if isfile(...) … else @info …` **arms-itself-later** idiom is the house answer to "the file
lands in a later wave"; use it for every Phase-14 source that does not exist at Wave 0.

**(f) The source-grep technique for SC2-c / SC2-d** — `spike\test\test_p13_tau.jl:102-104, 275-315`:

```julia
const TAU_SRC   = read(joinpath(@__DIR__, "..", "p13", "tau_probe.jl"), String)
const TAU_LINES = filter(l -> !startswith(strip(l), "#"), split(TAU_SRC, '\n'))
const TAU_CODE  = join(TAU_LINES, '\n')
...
@test !any(l -> occursin("P13_TAU_DELTA_GRID", l) && occursin(op, l), TAU_LINES)   # per-LINE form
@test occursin("roc_auc(", TAU_CODE)          # the helper IS used ...
@test !occursin("function roc_auc", TAU_CODE) # ... and NOT re-implemented
@test !occursin("corpus/", TAU_CODE)
@test !occursin("corpus\\", TAU_CODE)         # both path separators
```
The **per-line** variant (`TAU_LINES`) is what SC2-d needs (a `patch_correlation(` call compared to τ
on the same line); the **positive-plus-negative pair** (`occursin("roc_auc(")` *and*
`!occursin("function roc_auc")`) is exactly the no-reimplementation assertion Phase 14 should carry
for `roc_auc` / `corspearman` / `_bin_calibration`.

---

### 3.12 `spike/test/test_p14_fdr.jl` / `_conformal.jl` / `_fuse.jl` / `_posterior.jl`

Analog `spike\test\test_p13_result.jl` — **fixtures built ONCE at module level, one outer
`@testset ... verbose = true`, no simulation** (`:68-82`):

```julia
# --- Fixtures built ONCE --------------------------------------------------------------------
const RES_FIX_CAL  = _bin_calibration(fill(0.5, 100), vcat(trues(50), falses(50));
                                      n_bins = P13_ECE_NBINS)
const RES_FIX_META = p13_calibration_meta(RES_FIX_CAL; grid = 8, auc = 0.97)
const RES_FIX_OOD  = OODVerdict(1.23, false, (; density = 1.23, noise = 0.4))

"Build a valid three-hypothesis result; `meta` is where control draws live, if any."
_res(; coloc = 2.5, exclusion = -3.25, random = 0.0, meta = NamedTuple()) =
    ThreeHypothesisColocResult(8, zeros(7, 5),
                               (coloc = coloc, random = random, exclusion = exclusion),
                               RES_FIX_OOD, RES_FIX_META, meta)

@testset "P13 three-hypothesis result type (D-08)" verbose = true begin
```
The `_res(; kwargs...)` **keyword builder** is the pattern for SC2-a's 36-cell truth table: one
builder, 36 `@test` calls over the product, so a failure names the cell.

Also copy the `===`-not-`==` note (`:96-99`): *"`===` and not `==`: `-0.0 == 0.0` is true but they
are different bit patterns"* — directly relevant to `logbf.random === 0.0` and to
`p14_conformal_quantile`'s `q̂ != Statistics.quantile(s, 0.9)` assertion.

Every test file ends with the CPU-only testset (`test_p13_result.jl:214-216`, identical to
`test_p12_decoupling.jl:351-353`).

**`spike/test/fixtures/`:** existing contents are `p11_stage6_golden.jld2` and
`p12_stage1_golden.jld2` — **`.jld2` goldens captured by `capture_p1X_golden.jl` scripts**, not
hand-written literals. Phase 14's two fixtures are hand-built literals, so the in-file `const`
pattern above is the right analog; a `.jld2` under `fixtures/` would need a capture script. **NOT
VERIFIED:** I did not open `capture_p12_golden.jl`, so the capture-script pattern is unmapped.

---

## 4. Shared Patterns (apply to every Phase-14 file)

### 4.1 The license header
Every `spike/` and `src/` file opens with the 19-line AGPL block
(`spike\test\test_p12_decoupling.jl:1-19`, byte-identical across files). Copy it verbatim.

### 4.2 The guarded include
`isdefined(@__MODULE__, :SENTINEL) || include(joinpath(@__DIR__, ...))` — used at every include site
in the lane. The sentinel is a name the included file **alone** declares (§3.1b). The spike lane is a
**flat top-level namespace**: when the whole suite runs in one process, all files land in the same
module, so re-`const`-ing a name is a real collision (`test_p13_result.jl:58-66` guards
`P13_REPO_ROOT` for exactly this reason).

### 4.3 The path-depth reminder
From `spike/p13/` (hence `spike/p14/`) **the repo root is TWO levels up, not three**:
`joinpath(@__DIR__, "..", "..", "src", "results.jl")`. Stated in three separate files.

### 4.4 Argument validation
`throw(ArgumentError("<function_name>: <what was required>; got $(actual)"))` — function name first,
requirement, then the observed value. See `three_way_label` (`labels.jl:122-124`),
`ThreeHypothesisColocResult` (`result.jl:159-161`), `encode_conditioned_pair` (`net.jl:434-436`),
`colocalization_amortized` (`api.jl:112-120`), and RESEARCH's `p14_bayes_fdr` / `p14_fuse` /
`p14_conformal_quantile` examples — all identical.

### 4.5 Assertion messages name the breach and its decision
`@assert <cond> "D-01 decoupling breach: …"`, `"P13_DEV_SEED collides with a forbidden seed"`,
`"the gate key word collides with the training-pool key word"`. Never a bare `@assert`.

### 4.6 Comments record what was REJECTED and why
`P13_CONSTS_SHA`'s "REJECTED, and recorded so a later reader knows it was considered"; the Tier-1
self-check "RELOCATED BY PLAN 13-10, NOT DELETED"; `test_p12_decoupling.jl:279-283` correcting a
stale comment in place with *"true when written, false since"*. **Corrections are appended, never
overwritten** (CONTEXT §Specific Ideas). Phase 14 must carry the SC1/SC2 amendment citation in code,
artifact and report.

### 4.7 CPU-only, always
`@test !any(id -> occursin("CUDA", id.name), keys(Base.loaded_modules))` closes every test file; no
runner imports CUDA.

### 4.8 `runtests.jl` is NOT a gate
It exits 1 at the Phase-4 `SPEEDUP_GATE`, masking every later include block. Per-file runs
(`julia --project=spike spike/test/test_p14_<name>.jl`) are the only reliable signal
(`14-VALIDATION.md` Test Infrastructure).

### 4.9 Execute on the main working tree
No git worktrees — the OOD null fit needs the ~50 MB gitignored `spike/data/cache/p13/7c65a1a9…`
pool, which does not exist in a fresh worktree and was destroyed once already this milestone.

---

## 5. No Analog Found

| File | Role | Data flow | Reason |
|---|---|---|---|
| `spike/p14/fdr.jl` (the FDR *rule*) | utility | batch transform | **No FDR code exists in this repo.** Use `14-RESEARCH.md` §Code Examples `p14_bayes_fdr` verbatim; take only the *file shape* from `spike\validation\ood.jl:309-335`. |
| `spike/p14/conformal.jl` (the conformal *construction*) | utility | batch transform | **No conformal code exists; `ConformalPrediction.jl` appears nowhere but in the ROADMAP sentence.** Use RESEARCH §Code Examples verbatim; file shape from `roc_auc` / `ghat`. |
| `spike/test/fixtures/` Phase-14 fixtures | fixture | — | Existing fixtures are `.jld2` goldens with capture scripts; the Phase-14 fixtures are hand-built literals. Use the in-file `const RES_FIX_*` pattern instead. Capture-script pattern **NOT VERIFIED**. |

---

## 6. Conventions

**Convention derivation skipped** — `node bin/gsd-tools.cjs verify conventions --derive` returned
`{"mode":"derive","skipped":true,"reason":"no-readable-files","axes":[]}` both repo-wide and at
`--scope spike`. This is a Julia repository; the deriver's file-type set does not include `.jl`, so
no 4-axis table (file-name casing / identifier casing / export style / import style) can be produced
from the shared module. Recorded rather than fabricated.

**Observed conventions, derived by hand from the files read this session** (not a majority-vote
computation — treat as descriptive, and prefer the local directory's style in any conflict):

| Axis | Dominant | Evidence | Status |
|---|---|---|---|
| File-name casing | `snake_case.jl`, phase-prefixed (`p13_*.jl`, `run_p13_*.jl`, `test_p13_*.jl`) | every file in `spike/p13/`, `spike/test/`, `spike/validation/` | named contract |
| Identifier casing — functions | `snake_case`, phase-prefixed for phase-local (`p13_calibration_meta`, `p13_gate_rng`), leading `_` for file-private (`_p13_as_three_way_logbf`, `_bin_calibration`, `_iface_error`, `_recorded_ood_threshold`) | uniform across `src/` and `spike/` | named contract |
| Identifier casing — constants | `SCREAMING_SNAKE`, phase-prefixed (`P13_TAU`, `P14_ALPHA_CONFORMAL`) | `consts.jl` files, `gate_consts_*.jl` | named contract |
| Identifier casing — types | `UpperCamelCase` (`ThreeHypothesisColocResult`, `OODVerdict`, `CalibrationMeta`, `LocalColocMap`) | `src/results.jl`, `spike/p13/result.jl` | named contract |
| Export / module style | **No modules, no `export`, in the spike lane.** Flat top-level definitions loaded by `include`; the only `module` in `spike/p13/` is `_GC`, an *isolation* device, and its file says why it must sit outside the guard block. `src/` is a package with real exports. | `spike/p13/consts.jl:98-100`; `test_p13_result.jl:58-66` (the flat-namespace collision note) | **contested by directory, intentionally** |
| Import style | `using X` for the umbrella (`using Test`, `using JLD2`, `using Statistics`), `import X: sym` when only one symbol is wanted (`import Random123: Philox4x`), `import StatsBase` + qualified calls when the name would collide (`StatsBase.corspearman`) | `run_three_way_gate.jl:63-69`; `consts.jl:110`; `result.jl:68-70` | named contract |

**Contested hotspots (author's choice).** The `src/` ↔ `spike/` split is this repository's prototype
*intentional* contested axis, and it is the exact analogue of the CJS↔SDK dual resolver: `src/**` is
a **Julia package** — modules, `export`, a `Project.toml`, public API in `src/amortized/api.jl`;
`spike/**` is a **flat script lane** — no modules, no exports, guarded `include`s, a separate frozen
`Project.toml`. Each half is internally consistent per-directory and contested only repo-wide.
**Phase 14 lives in `spike/p14/`, so it takes the spike half: no `module`, no `export`, guarded
includes, phase-prefixed names.** Reviewers and planners should match the directory's local style
rather than a repo-wide dominant. Note that `src/` is additionally **byte-frozen** here by D-01 and by
a running test, so the question of which style to write in `src/` does not arise this phase.

---

## Metadata

**Analog search scope:** `spike/p13/`, `spike/test/`, `spike/validation/`, `spike/simulator/`,
`spike/comparator/` (via RESEARCH), `src/`, `src/amortized/`
**Files read this session:** 18 (`14-CONTEXT.md`, `14-VALIDATION.md`, `14-RESEARCH.md` §§Arch-map /
Stack / D / Runtime-inventory / Code-Examples, `spike/test/test_p12_decoupling.jl` (full),
`spike/p13/consts.jl` (§§A-D, F, J, self-checks, Tier-2), `spike/p13/result.jl`,
`spike/p13/net.jl:420-583`, `spike/p13/run_three_way_gate.jl` (banner, helpers, main §§0-4, save
block), `spike/p13/labels.jl:100-135`, `spike/test/test_p13_consts.jl:20-169`,
`spike/test/test_p13_result.jl:20-100, 180-219`, `spike/test/test_p13_tau.jl` (grep lines),
`spike/validation/sbc.jl:55-144`, `spike/validation/ood.jl:300-349`, `spike/simulator/ghat.jl:93-121`,
`src/results.jl:30-129`, `src/amortized/ood.jl:318-398`, `src/amortized/local_map.jl:55-125`,
`src/amortized/api.jl:74-149`, `spike/p13/run_p13_realimage.jl:412-445`)
**Not opened (gaps, marked NOT VERIFIED above):** `spike/test/test_p13_labels.jl`,
`spike/test/capture_p12_golden.jl`, `spike/p13/preconditions.jl`, `spike/comparator/classical.jl`
(signatures taken from RESEARCH §D.5, which records them as VERIFIED on disk this session),
`spike/p12/result.jl`
**Pattern extraction date:** 2026-08-03
