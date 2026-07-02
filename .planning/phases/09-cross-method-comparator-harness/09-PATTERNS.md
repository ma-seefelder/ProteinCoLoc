# Phase 9: Cross-Method Comparator Harness - Pattern Map

**Mapped:** 2026-07-02
**Files analyzed:** 8 (6 new, 2 modified) — spike-local only, NO `src/` edits
**Analogs found:** 8 / 8 (all files have a strong in-repo or sibling-repo analog)

> **Read-only reuse discipline (CLAUDE.md hard constraint):** every analog below is
> consumed read-only. The comparator reaches `src/colocalization.jl` (`patch`,
> `correlation`) and `MultiChannelImage` ONLY transitively through `spike/contract.jl`'s
> existing `include()` — never a second `include` of `src/`, never an edit.
>
> **CACHE-INVALIDATION TRAP (RESEARCH Pitfall 1 / hashguard.jl:48):** `spike/data/encode.jl`
> is entry 5 of `HASH_SRC_FILES`. Any byte change to it (even a "factor Manders into a
> callable" refactor) flips `cache_hash` and auto-invalidates the completed ≥50k Phase-3
> cache + the trained Phase-4 NPE. **Do NOT plan any edit to `encode.jl`.** D-04's intent is
> satisfied by reproducing the formula in the NEW module and pinning it with an equality test.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `spike/comparator/classical.jl` (new) | utility (estimator wrappers: manders, pearson_whole, spearman/patch_correlation, costes_p) | transform | `spike/data/encode.jl:encode_aug` (+ `spike/data/seeding.jl` for Costes RNG salt) | exact (formula reuse) |
| `spike/comparator/config.jl` (new) | config (pre-declared consts / thresholds) | (declarative) | `spike/test/test_npe.jl:70-96` pre-registered consts + `spike/simulator/prior.jl:40-57` module consts | exact |
| `spike/comparator/table.jl` (new) | service (DataFrame assembly + divergence/traffic-light + CSV/JLD2 writer, content-hash dir) | transform + file-I/O | `spike/data/cache.jl` (atomic hash-named JLD2) + BayesInteractomics `copula_diagnostics.jl:88` (traffic-light) | role-match (composite) |
| `spike/comparator/tapqir_bridge.jl` (new) | provider (optional external bridge, graceful skip-with-flag) | request-response (out-of-process) | `spike/baseline/Project.toml` (isolated sub-env precedent) + runtests.jl CUDA graceful-degradation gate | role-match |
| `spike/comparator/run_comparator.jl` (new) | controller / entry point (seeded driver) | batch | `spike/data/generate.jl` (guarded-include driver + `generate_sample` seeded loop) | exact |
| `spike/test/test_comparator.jl` (new) | test | request-response | `spike/test/test_npe.jl` (scaffold + pre-registered consts + fixture-oracle idiom) | exact |
| `spike/test/runtests.jl` (modify) | test (harness aggregator) | request-response | `spike/test/runtests.jl` itself (append `include` + extend resolve-risk block d/e/f) | exact (self) |
| `spike/Project.toml` (modify) | config (promote DataFrames + CSV transitive→direct) | (declarative) | `spike/Project.toml` itself + Phase-3 JLD2/Random123 promotion precedent | exact (self) |

**Note on file layout (Claude's discretion, D-05/D-09):** the task file list folds `costes_p`
into `classical.jl` and the divergence/audit into `table.jl`. RESEARCH §Recommended Structure
suggested finer splits (`costes.jl`, `inputs.jl`, `audit.jl`). Either is acceptable; the
excerpts below are grouped by the task's file list and each notes where a split is natural.

## Pattern Assignments

### `spike/comparator/classical.jl` (utility, transform)

**Analog:** `spike/data/encode.jl:encode_aug` (lines 96-137) — the Manders M1/M2 +
whole-image Pearson formulas and the `_safe` NaN-not-throw reduction idiom. **Reproduce
byte-for-byte in the new module; do NOT edit `encode.jl`** (hash-guarded).

**`_safe` degenerate-guard idiom to copy** (`encode.jl:80`):
```julia
_safe(f, v::AbstractVector) = isempty(v) ? NaN : f(v)
```

**Manders + whole-image Pearson formula to reproduce exactly** (`encode.jl:99-109`) —
same `mci.otsu_threshold` source is load-bearing (Specific Idea: comparator and
`encode_aug` must never silently diverge on M1/M2):
```julia
t1  = mci.otsu_threshold[1]
t2  = mci.otsu_threshold[2]
s1 = sum(ch1)
s2 = sum(ch2)
M1 = s1 > 0 ? sum(ch1 .* (ch2 .> t2)) / s1 : NaN   # frac. ch1 signal where ch2 present
M2 = s2 > 0 ? sum(ch2 .* (ch1 .> t1)) / s2 : NaN   # frac. ch2 signal where ch1 present
pearson_whole = cor(vec(ch1), vec(ch2))
```

**Pearson/Spearman patch-grid reuse** — reach the FROZEN `src/` math through
`contract.jl`, NOT a new `include`. Mirror `spike/contract.jl:patch_summary` (lines 85-90)
which already calls `patch.(...,8)` and `correlation(...; method=...)`:
```julia
# src/colocalization.jl:221 correlation() builds Dict(:pearson=>cor, :spearman=>corspearman,
#   :kendall=>corkendall) at CALL time — so StatsBase must be in scope (contract.jl already
#   `using StatsBase`). Patches with ≤15 surviving pixels become `missing` (src:235);
#   _exclude_zero (src:154/187) drops 0.0/NaN/missing first. Whole-image scalar per D-04:
spearman_whole(mci) = corspearman(vec(mci.data[1]), vec(mci.data[2]))
# patch-grid variant (mean over non-missing patches), if reported:
#   ρ = correlation(patch.([x,y],8)...; method=:spearman); mean(collect(skipmissing(ρ)))
```

**Costes-p (NEW, D-05) — seeded RNG salt idiom** copied from `spike/data/seeding.jl:45-46,86-87`
(disjoint XOR-salted Philox stream so Costes scrambles are bit-reproducible and independent
of the input-generation stream):
```julia
using Random123
const COSTES_SALT = 0x...   # nonzero, distinct from HOLDOUT_SALT/FOLD_SALT
costes_rng(master_seed, idx) = Philox4x(UInt64, (UInt64(master_seed) ⊻ COSTES_SALT, UInt64(idx)))
# Costes et al. 2004: block-scramble ch2 (block side ≈ resolution element, ~7px for σ_psf=1.3),
# recompute cor(vec(ch1),vec(scr)); p = (#{r_scr ≥ r_obs} + 1)/(n + 1)  (Davison–Hinkley +1/+1).
# BLOCK scramble, NOT pixel scramble (RESEARCH Pitfall 2 — pixel scramble biases p→0).
```

**Consistency test hook (D-04 intent, no encode.jl edit):** `test_comparator.jl` must assert
`manders(mci)` equals the M1/M2 that `encode_aug` produces for the same `mci` (slice the two
known moment positions from `encode_aug`'s output — see AUG layout at `encode.jl:128-134`, M1/M2
are the first two moments, i.e. indices 129,130 of the 142-vec).

---

### `spike/comparator/config.jl` (config, declarative)

**Analog:** `spike/test/test_npe.jl:70-96` (pre-registered `const` block committed BEFORE any
run — the anti-data-snooping discipline D-14 requires) and `spike/simulator/prior.jl:40-57`
(module-level `const` priors/tolerances with documented rationale).

**Pre-registered-consts idiom to copy** (`test_npe.jl:70-78`):
```julia
# --- Pre-registered constants (declared BEFORE any reported run; RESEARCH A2/A3) --
const NPE_RMSE_TOLERANCE   = 1.2       # D-09 ...
const SPEEDUP_GATE         = 100.0     # NPE-03 ...
const NPE_MASTER_SEED      = 0xC0FFEE  # Random123 reproducibility for the Phase-4 runs
```

**Phase-9 consts to declare (D-14, documented in artifact header):**
`COSTES_N_SCRAMBLE=200`, `COSTES_BLOCK_PX≈7`, `COSTES_SALT`, `DIVERGENCE_WARN`, `DIVERGENCE_FAIL`,
`MASTER_SEED`, `TAPQIR_TOL`, `TAPQIR_PUBLISHED_VALUE`. Mirror `prior.jl:42-46`'s "fixed BEFORE
the run, NOT tuned to pass" comment style.

---

### `spike/comparator/table.jl` (service, transform + file-I/O)

**Analog A — content-addressed atomic JLD2 artifact:** `spike/data/cache.jl` + `spike/data/hashguard.jl`.

**Content-hash naming to mirror** (`hashguard.jl:76-84`, `cache.jl:147-164`) — the artifact dir
is named by a SHA-256 over the config + sources so it is bit-reproducible and auto-separating:
```julia
function cache_hash(config::NamedTuple; src_files = HASH_SRC_FILES)
    ctx = SHA.SHA256_CTX()
    for f in src_files; SHA.update!(ctx, read(f)); SHA.update!(ctx, UInt8[0x00]); end
    SHA.update!(ctx, Vector{UInt8}(canonical(config)))
    return bytes2hex(SHA.digest!(ctx))
end
```
> **Trap:** do NOT add `encode.jl` (or any comparator source) to the EXISTING `HASH_SRC_FILES`
> in `hashguard.jl` (that would re-hash the Phase-3 cache). The comparator table needs its OWN
> hash over its OWN config + comparator sources — a separate `cache_hash` call with a distinct
> `src_files` list, or a fresh helper in `table.jl`.

**Atomic write idiom to copy** (`cache.jl:97-108`): `jldsave` to `.tmp` → reopen integrity-check →
`mv(tmp, path; force=true)`. CSV writer is analogous (write `.tmp`, `mv`).

**Analog B — divergence / traffic-light column:** BayesInteractomics `copula_diagnostics.jl:88`
status idiom (the exact `:pass/:warn/:fail` ternary; render as `:green/:amber/:red`):
```julia
status = ks_t < ks_warn ? :pass : ks_t < ks_fail ? :warn : :fail
# Phase-9 form (thresholds from config.jl, pre-declared — D-10/D-14):
traffic_light(d) = d < DIVERGENCE_WARN ? :green : d < DIVERGENCE_FAIL ? :amber : :red
```

**Analog C — tidy DataFrame table type:** BayesInteractomics `_bin_calibration`
(`calibration.jl:22-72`) builds parallel `Float64[]`/`Int[]` column vectors then packs a struct;
the DataFrame equivalent is one row per shared input with columns
`{input_id, regime, Costes_p, M1, M2, Pearson, Spearman, NPE_ρ̂, NPE_OOD_flag, divergence, light}`
(D-09). NPE columns are OPTIONAL (`missing` when Phase-4 artifacts absent — RESEARCH OQ3).

---

### `spike/comparator/tapqir_bridge.jl` (provider, out-of-process request-response)

**Analog A — isolated sub-env precedent:** `spike/baseline/Project.toml` — the spike ALREADY
isolates a heavy, resolve-risky stack (Turing/AdvancedVI/ForwardDiff) in its own
`Project.toml`+`Manifest.toml` under `spike/baseline/`. Replicate exactly for PythonCall/CondaPkg
under `spike/comparator/tapqir_env/`. **Never add PythonCall/CondaPkg to the main
`spike/Project.toml`** (RESEARCH Pitfall 4/5 — NeuralEstimators v0.2.1 co-resolve landmine).

**Analog B — graceful skip-with-flag idiom:** the runtests.jl CUDA graceful-degradation gate
(`runtests.jl:52-58`) proves the "optional heavy stack must never hard-fail the core path"
discipline. The bridge returns a status NamedTuple; on any failure it skips, classical table
stays green (D-08):
```julia
function tapqir_anchor(; tol = TAPQIR_TOL)
    avail = try _tapqir_env_available() catch; false end
    avail || return (status = :skipped, reason = "Tapqir/Python env unavailable", value = missing)
    try
        recovered = _run_tapqir_tutorial()          # isolated sub-env / subprocess (RESEARCH OQ2)
        return (status = abs(recovered - TAPQIR_PUBLISHED_VALUE) ≤ tol ? :passed : :failed,
                value = recovered)
    catch e
        return (status = :skipped, reason = sprint(showerror, e), value = missing)
    end
end
```
> RESEARCH OQ1/A3: exact tutorial dataset + published anchor value/tolerance are UNCONFIRMED —
> run the tutorial once in the isolated env during Wave-1 to pin `TAPQIR_PUBLISHED_VALUE`/`TAPQIR_TOL`
> in `config.jl`. Until then the D-13 test accepts "anchor within tol OR skipped-with-flag."

---

### `spike/comparator/run_comparator.jl` (controller, batch)

**Analog:** `spike/data/generate.jl` — the seeded driver idiom: guarded `isdefined`-includes to
compose the frozen chain, then a keyed per-sample loop.

**Guarded-include composition to copy** (`generate.jl` head — contract.jl FIRST):
```julia
isdefined(@__MODULE__, :build_mci)     || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :sample_prior)  || include(joinpath(@__DIR__, "..", "simulator", "prior.jl"))
isdefined(@__MODULE__, :simulate_pair) || include(joinpath(@__DIR__, "..", "simulator", "forward.jl"))
isdefined(@__MODULE__, :HOLDOUT_SALT)  || include(joinpath(@__DIR__, "..", "data", "seeding.jl"))
```

**Seeded shared-input generation to copy** (`generate.jl:generate_sample`) — the θ→image→mci path
with a Philox stream keyed by `(master_seed, idx)`; each row carries ground-truth θ (D-02):
```julia
rng = sample_rng(master_seed, idx)          # spike/data/seeding.jl:75 (Philox4x, thread-independent)
θ   = sample_prior(rng)                       # simulator/prior.jl:71 — ρ_true carries the regime label
mci = build_mci(simulate_pair(rng, θ; imsize=isz))   # contract.jl:60 + forward.jl:110 (UNCHANGED)
```
D-03: also accept a caller-supplied `Vector{MultiChannelImage}` (skip generation) so the Phase-8
corpus can be fed later without redesign.

---

### `spike/test/test_comparator.jl` (test)

**Analog:** `spike/test/test_npe.jl` — the whole file is the template.

- **Header + guarded-include + consts-before-gate scaffold:** `test_npe.jl:45-96`.
- **Wave-0 skipped-scaffold idiom** (one `@test_skip true` per SC, tagged with the CMP-* req):
  mirror how `test_npe.jl` enumerates SC1..SC5 as named child testsets under one
  `@testset ... verbose = true`.
- **Fixture-oracle idiom for the divergence/traffic-light rule** (`test_npe.jl:274-280, 303-313`):
  test `traffic_light`/divergence against HAND-BUILT inputs with KNOWN expected colors (an
  independent oracle), not by re-deriving the function's own formula — catches a `<`/`>` flip
  copied into both code and test.
- **Determinism gate (D-13a):** two seeded `run_comparator` runs → byte-identical JLD2 + identical
  content-hash dir. Mirror the BIT-EXACT `==` (not `≈`) rule at `test_npe.jl:408-427`.
- **Finiteness gate (D-13b):** Costes_p, M1, M2, Pearson, Spearman all `isfinite` on a
  non-degenerate fixture.
- **Manders-equality gate (D-04):** `manders(mci)` == the M1/M2 slice of `encode_aug(mci, M)`.
- **Tapqir gate (D-13c):** `status ∈ (:passed,:skipped)` (accept skip when env absent).

---

### `spike/test/runtests.jl` (modify — test aggregator)

**Analog:** the file itself. Two additive edits, both matching existing patterns:

1. **Append the include** after `test_npe.jl` (mirror lines 92-103):
   ```julia
   include(joinpath(@__DIR__, "test_comparator.jl"))
   ```
2. **Extend the resolve-risk block (d/e/f → add g)** at `runtests.jl:59-86` — re-assert the
   NeuralEstimators v0.2.1 pin AFTER DataFrames/CSV promotion, and assert PythonCall/CondaPkg
   are ABSENT from the main env (they live only in `tapqir_env/`):
   ```julia
   # (g) RESOLVE-RISK GATE (Phase 9): DataFrames/CSV promoted transitive→direct must NOT
   #     downgrade the pin; PythonCall/CondaPkg must NEVER enter the main spike env.
   @test haskey(_deps_by_name, "DataFrames")
   @test haskey(_deps_by_name, "CSV")
   @test !haskey(Pkg.project().dependencies, "PythonCall")
   @test !haskey(Pkg.project().dependencies, "CondaPkg")
   @test Pkg.dependencies()[Base.UUID("38f6df31-6b4a-4144-b2af-7ace2da57606")].version == v"0.2.1"
   ```
   (UUID/pattern identical to the existing (d)/(e)/(f) assertions — copy their shape verbatim.)

---

### `spike/Project.toml` (modify — config)

**Analog:** the file itself + the Phase-3 JLD2/Random123 promotion precedent (already present as
direct deps, lines 12,14). Add two lines to `[deps]` (UUIDs from the resolved Manifest — they are
already transitive, so promotion should NOT change resolved versions, RESEARCH A2):
```
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
```
> The DataFrames UUID is confirmed in `spike/baseline/Project.toml:4`. Re-freeze
> `spike/Manifest.toml` via `julia --project=spike -e 'using Pkg; Pkg.add(["DataFrames","CSV"])'`.
> This is a Wave-0 env task mirroring Phase-3 03-01. **Do NOT touch the root `Project.toml`.**

## Shared Patterns

### NaN-safe degenerate-input reduction (D-06)
**Source:** `spike/data/encode.jl:80` (`_safe`) and `spike/contract.jl:100-107` (`induced_mu`
returns NaN on empty `skipmissing`).
**Apply to:** every estimator in `classical.jl` — return a NaN-safe scalar, never throw, on a
degenerate (all-zero / all-missing) input.

### Reproducible seeding: Philox keyed by (seed, idx) with a disjoint salt (D-12)
**Source:** `spike/data/seeding.jl:45-46` (salt consts), `:75-87` (`sample_rng`/`holdout_rng`).
**Apply to:** the θ-grid generation (`run_comparator.jl`, reuse `sample_rng`) and the Costes
scrambles (`classical.jl`, NEW `COSTES_SALT` stream). Thread- and order-independent by construction.

### Read-only `src/` reach via contract.jl only
**Source:** `spike/contract.jl:45-47` (the sole `include` of `src/`), `:85-90` (`patch_summary`
using `patch`/`correlation`).
**Apply to:** all Pearson/Spearman/patch math. Never add a second `include("../src/...")`
(double-defines the frozen functions — RESEARCH Anti-Pattern).

### Content-addressed atomic artifact (D-09/D-12)
**Source:** `spike/data/cache.jl:97-134` (atomic `.tmp`→integrity-check→`mv`) + `hashguard.jl:76-84`.
**Apply to:** the table writer in `table.jl` — but with its OWN hash config/src list, never by
mutating the Phase-3 `HASH_SRC_FILES`.

### Traffic-light status ternary (D-10/D-11)
**Source:** BayesInteractomics `copula_diagnostics.jl:88`; `types.jl:391` (`status::Symbol`);
`types.jl:688` (`:pass→PASS` render).
**Apply to:** the divergence column in `table.jl` and the audit summary.

### Audit report as an IOBuffer Markdown table (D-11 / SC3)
**Source:** BayesInteractomics `predictive_checks.jl:886-924` (`generate_diagnostics_report`:
`IOBuffer`, `println(io, "| Metric | Value |")` …) and `_ks_test_uniform` (`:729-736`).
**Apply to:** the audit summary (band counts per traffic-light color + optional KS-uniformity on
Costes p-values), written alongside the CSV/JLD2 artifact.

### Isolated sub-env for a resolve-risky external stack
**Source:** `spike/baseline/Project.toml` (Turing/AdvancedVI isolated from the main spike env).
**Apply to:** `spike/comparator/tapqir_env/` (PythonCall/CondaPkg). `.gitignore` the CondaPkg
materialized env dir; commit only `tapqir_env/Project.toml`+`Manifest.toml`.

### Pre-registered constants (anti-data-snooping, D-14)
**Source:** `spike/test/test_npe.jl:70-78` and `spike/simulator/prior.jl:42-46`.
**Apply to:** `config.jl` — all thresholds committed before the first table, documented in the
artifact header.

### Resolve-risk env gate (NeuralEstimators v0.2.1 pin)
**Source:** `spike/test/runtests.jl:59-86` (gates d/e/f, keyed by NeuralEstimators UUID
`38f6df31-6b4a-4144-b2af-7ace2da57606`).
**Apply to:** the new gate (g) after DataFrames/CSV promotion; assert PythonCall/CondaPkg absent.

## No Analog Found

| File | Role | Data Flow | Reason |
|------|------|-----------|--------|
| — | — | — | Every Phase-9 file has a strong in-repo or sibling-repo analog. |

The single genuinely NEW algorithm (Costes block-scramble p-value) has no in-repo analog for the
*statistic itself*, but its seeding, NaN-safety, and pre-registration all reuse established
spike patterns (seeding.jl salts, encode.jl `_safe`, test_npe.jl consts). Cite Costes et al.
2004 (Biophys J 86:3993) for the block-randomization scheme; planner should use RESEARCH
§Pattern 2 for the concrete recipe (block side ≈7px, N=200, Davison–Hinkley +1/+1).

## Conventions

> Deterministic derivation SKIPPED (`gsd-tools verify conventions` → `no-readable-files`:
> the tool does not parse `.jl` sources). Table below is derived manually from the 8 spike/
> and src/ files read this session — treat as observed majority patterns, not a tool contract.

| Axis | Dominant | Share | Entropy | Status |
|------|----------|-------|---------|--------|
| File-name casing | lower_snake / lowercase (`encode.jl`, `run_comparator.jl`, `test_npe.jl`) | ~100% | low | named contract |
| Identifier casing | snake_case functions & consts, `UPPER_SNAKE` module consts (`cache_hash`, `sample_rng`, `HASH_SRC_FILES`, `COSTES_SALT`) | ~95% | low | named contract |
| Export style | no `module`/`export` — plain top-level `function` defs composed via guarded `include()` (`isdefined(@__MODULE__, :x) || include(...)`) | ~100% in spike/ | low | named contract |
| Import style | `using` at file top (`using JLD2`, `using Random123`), commented with the exact symbols used; frozen `src/` reached ONLY via `contract.jl` include | ~100% | low | named contract |

**Additional named contracts observed in every spike source (match them):**
- **GNU AGPL license header** (the 19-line `#=...=#` block) opens every `.jl` file — new comparator
  files must carry it verbatim (see any analog, lines 1-19).
- **Module-doc comment** after the header: `# spike/<path> --- <REQ-ID>: <one-line purpose>` then a
  DECOUPLING note restating the hard constraint (e.g. `encode.jl:40-42`, `cache.jl:36-37`).
- **Requirement/decision inline tags** (`D-04`, `CR-02`, `T-03-06`, `RESEARCH Pitfall N`) annotate
  non-obvious choices — the planner should mint `CMP-*` tags (RESEARCH §Phase Requirements) and use
  them the same way.

**Contested hotspots (author's choice):** none within `spike/`. The repo-wide CJS↔SDK dual-resolver
prototype (`bin/lib/**` CommonJS vs `sdk/src/**` ESM) does NOT exist here — this is a single-language
Julia project, so there is no dual-module split to match. The one legitimate internal boundary is
**main-env vs isolated-sub-env** (`spike/Project.toml` vs `spike/baseline/Project.toml` /
`spike/comparator/tapqir_env/Project.toml`): each is internally consistent; the Tapqir bridge MUST
live in its own sub-env, never the main one.

## Metadata

**Analog search scope:** `spike/` (contract, data, simulator, npe, test, baseline), `src/colocalization.jl`,
and sibling `~/Documents/GitHub/BayesInteractomics/src/diagnostics/`.
**Files scanned:** 14 (encode.jl, contract.jl, hashguard.jl, cache.jl, seeding.jl, forward.jl,
prior.jl, generate.jl, runtests.jl, test_npe.jl, Project.toml, baseline/Project.toml,
src/colocalization.jl; BI calibration.jl/types.jl/copula_diagnostics.jl/predictive_checks.jl).
**Pattern extraction date:** 2026-07-02
