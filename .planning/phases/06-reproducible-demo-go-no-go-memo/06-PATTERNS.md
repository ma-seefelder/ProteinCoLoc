# Phase 6: Reproducible Demo + Go/No-Go Memo - Pattern Map

**Mapped:** 2026-07-03
**Files analyzed:** 2 (1 Julia entry script, 1 markdown memo)
**Analogs found:** 2 / 2 (both exact/strong role matches; the memo composites several analogs)

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `spike/demo.jl` | entry-script (orchestrator) | batch / request-response (two-tier: fast load vs `--full` recompute) | `spike/02_simulator_demo.jl` (top-level idiom) + `spike/validation/run_sbc.jl` (gate/report structure) | exact (idiom) + role-match (gate) |
| `.planning/phases/06-.../06-GO-NO-GO-MEMO.md` | documentation (decision memo) | transform (report on-disk numbers → prose verdict) | `.planning/phases/05-.../05-04-SUMMARY.md` (Phase SUMMARY prose) + `spike/NOTES.md` (decision-record tables) | role-match |

**Key insight (carried from RESEARCH):** `demo.jl` is ~90% orchestration + ~10% new code — it *composes* frozen read surfaces (`harness.jl`, `run_*.jl`, `*_report.jld2`) and re-implements nothing. Every pattern below already exists as a tested surface; copy it, do not invent.

---

## Pattern Assignments

### `spike/demo.jl` (entry-script, two-tier batch)

**Primary analog:** `spike/02_simulator_demo.jl` (top-level seeded CPU-only script idiom)
**Secondary analog:** `spike/validation/run_sbc.jl` (report-load / gate / `_save` structure)
**Composition surfaces:** `spike/validation/harness.jl`, `spike/data/seeding.jl`, `spike/validation/consts.jl`, `spike/npe/train_npe.jl` (`load_npe`), `spike/validation/train_ratio.jl` (`load_ratio`)

#### 1. AGPL header (copy verbatim, lines 1-19 of every spike file)
Every `spike/*.jl` opens with the identical 19-line `#= … =#` AGPL block. Copy it byte-for-byte to the top of `demo.jl`.
```julia
# Source: spike/02_simulator_demo.jl:1-19 (identical in 00_smoke.jl, harness.jl, run_sbc.jl, seeding.jl, consts.jl)
#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder
...
You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=#
```
Follow the header with a `# spike/demo.jl --- …` block-comment banner describing purpose + the PowerShell-safe verify command, exactly as `02_simulator_demo.jl:21-40` and `run_sbc.jl:21-41` do.

#### 2. `using` block + headless CairoMakie activate (only if plotting) (lines 42-48)
```julia
# Source: spike/02_simulator_demo.jl:42-48
using CairoMakie
using StatsBase
using Statistics
using Random
using Distributions
CairoMakie.activate!()                       # headless raster/vector, no OpenGL/display
```
`demo.jl` additionally needs `using JLD2, Test, Dates` (from `run_sbc.jl:44-47`). If `demo.jl` renders no new composite figure, `CairoMakie` may be omitted (reuse existing `figures/*.png` per RESEARCH "Don't Hand-Roll").

#### 3. Guarded, order-dependent `include()` of composition surfaces (idempotent)
The harness establishes the canonical guarded-include idiom — `isdefined(@__MODULE__, :symbol) || include(...)` — so `demo.jl` loads cleanly standalone AND under `runtests.jl`.
```julia
# Source: spike/validation/harness.jl:51-58 (ORDER MATTERS: consts → trainer/infer → simulator → contract → encode → seeding)
isdefined(@__MODULE__, :VAL_MASTER_SEED) || include(joinpath(@__DIR__, "validation", "consts.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "validation", "harness.jl"))
# harness.jl transitively pulls load_npe (train_npe.jl), posterior_for (infer.jl),
# sample_prior/simulate_pair, build_mci (contract.jl → read-only src/), encode_d01, sample_rng.
```
Note the `@__DIR__` + `joinpath` path convention (never bare relative paths) used everywhere: `02_simulator_demo.jl:51-53`, `harness.jl:51-58`, `run_sbc.jl:50-51`.

#### 4. Bare-ARGS `--full` flag (NO ArgParse) + script-vs-include guard
```julia
# Source: spike/npe/run_thread_sweep.jl:215 (bare-ARGS + PROGRAM_FILE guard)
if abspath(PROGRAM_FILE) == (@__FILE__) && !isempty(ARGS) && ARGS[1] == "--worker"
```
Adapt to:
```julia
const FULL = ("--full" in ARGS)
# optional: wrap the driver in `if abspath(PROGRAM_FILE) == @__FILE__ … end`
# (train_ratio.jl:369 uses this exact guard so an include never triggers the run)
```

#### 5. Step-0 git-status decoupling assertion (SC2, DEMO-02) — argv form, no shell string
```julia
# Source: RESEARCH Pattern 4 (verified `git status --porcelain -- src/ …` currently empty this session).
# Julia backticks pass argv directly (no shell) — keep args as literal tokens, never interpolate.
out = read(`git status --porcelain -- src/ src/bayes.jl src/colocalization.jl`, String)
@assert isempty(strip(out)) "DECOUPLING VIOLATED: src/ is dirty:\n$out"
```
Run this FIRST, before any compute (fail fast). Scope the hard assertion to the three `src/` paths — a whole-tree `git status` false-positives on `.planning/` churn (RESEARCH Pitfall 4). Any whole-tree read is informational only, never a gate.

#### 6. Fast-tier chain proof on the FIXTURE seed (never the reserved stream)
Mirror `test_sbc.jl` exactly: fixture constants + `val_rng(VAL_FIX_SEED)`, load the frozen model ONCE, one `draw_simulate_infer` pass.
```julia
# Source: spike/test/test_sbc.jl:40-42 + spike/validation/harness.jl:78-113
m = load_frozen_model()                                   # loads trained_npe.jld2 ONCE (frozen zt/θzt)
r = draw_simulate_infer(m, val_rng(VAL_FIX_SEED); imsize = (64,64), N = SBC_FIX_L)
# r.θ (true 7-tuple), r.draws (7×N physical-θ posterior). CPU-only; NO fit(ZScoreTransform).
```
**Reproducibility assert (RESEARCH Open Q3):** draw the fixture chain twice on the same `VAL_FIX_SEED` and `@assert` the two posterior draws are `==` (Philox4x is deterministic) — the cheapest honest SC1 proof.
Anti-pattern (RESEARCH): never touch `val_rng(VAL_MASTER_SEED)` in the fast tier; never call `fit(ZScoreTransform,...)` (transforms are frozen in `m.zt`/`m.θzt`).

#### 7. Fast-tier headline-number load from `*_report.jld2` (VERIFIED keys) — guard for fresh clone
```julia
# Source: RESEARCH "Code Examples" (report keys read from disk this session); reports are GITIGNORED.
const SBC_REPORT = joinpath(@__DIR__, "validation", "sbc_report.jld2")
if isfile(SBC_REPORT)
    sbc = jldopen(SBC_REPORT, "r") do f
        (labels=f["labels"], ks_p=f["ks_p"], chi2_p=f["chi2_p"], ece=f["ece"],
         verdict=f["verdict"], M=f["SBC_M"], L=f["SBC_L"])
    end
else
    @warn "sbc_report.jld2 absent (gitignored); run --full to regenerate (post-iteration numbers)"
end
```
Verified keys — `sbc_report.jld2`: `rank_table,labels,ks_p,chi2_p,ece,mce,verdict,cov_nominal,cov_empirical,SBC_M,SBC_L,SBC_BINS,…,caption,elapsed_min,generated`; `bf_report.jld2`: `targets,logbf_amortized,logbf_kde,corr,max_abs_err,log_prior_odds,…`; `ood_report.jld2`: `families,levels,maha_auc,fam_best_auc,combined_auc,id_fire_rate,neg_ks,pp_channel_note,…`.
**Assumption A1 (RESEARCH):** print `bf_report["max_abs_err"]` so the memo can confirm the `10.95` figure before quoting it. **Pitfall 1:** these are POST-ITERATION (set 2) numbers — the demo table must annotate them as such, never as pre-registered.

#### 8. `--full` dispatch to the reported gates (subprocess, independent exit codes)
```julia
# Source: RESEARCH Open Q2 recommendation (each run_*.jl calls main() at file scope + @testset hard gate)
if FULL
    for s in ("run_sbc.jl", "run_bf.jl", "run_ood.jl")
        run(`julia --project=spike $(joinpath(@__DIR__, "validation", s))`)  # each exits nonzero on gate fail
    end
end
```
Prefer subprocess spawn over in-process `include` so one gate's `@testset` throw does not abort the others (RESEARCH Open Q2). The `run(\`julia --project=… script\`)` subprocess form is already used in `run_bf.jl:99` for the isolated baseline.

#### 9. Success-criteria table + self-assert close (PowerShell-safe, no POSIX test)
```julia
# Source: spike/validation/run_sbc.jl:141-160 (rpad column table) + 02_simulator_demo.jl:148-152 (@assert self-check)
println("-"^78)
println(rpad("criterion", 12), rpad("check", 40), "status")
# … one row per SC1/SC2/SC3 + DEMO-01/02/03 …
@assert isfile(MEMO_PATH) "demo FAILED: missing memo $MEMO_PATH"   # DEMO-03 presence gate
println("demo OK: …")
```
The script self-asserts (`@assert`/`@test`) so a missing artifact fails the run — the verify command is simply `julia --project=spike spike/demo.jl` (no POSIX `test -f`; CLAUDE.md Windows-tauglich).

#### 10. Atomic JLD2 report save (ONLY relevant if `--full` writes a demo-level summary; otherwise skip)
```julia
# Source: spike/validation/run_sbc.jl:56-65 (tmp → integrity-check → mv; identical in run_bf.jl:60-69, run_ood.jl:78-87)
function _save_report(path; kwargs...)
    mkpath(dirname(path)); tmp = path * ".tmp"
    jldsave(tmp; kwargs...)
    JLD2.jldopen(tmp, "r") do f; @assert haskey(f, "<key>") "integrity check failed ($tmp)"; end
    mv(tmp, path; force = true); return path
end
```
`demo.jl` primarily *reads* reports; only reach for this if it persists its own composite artifact.

---

### `.planning/phases/06-.../06-GO-NO-GO-MEMO.md` (documentation, transform)

**Primary analog:** `.planning/phases/05-.../05-04-SUMMARY.md` (Phase SUMMARY prose — the source of the verbatim Set 1 numbers).
**Secondary analog:** `spike/NOTES.md` (decision-record with labeled tables + "Explicitly excluded" rationale rows — the closest in-repo precedent for a decision memo with tables and honest tradeoffs).

The memo is prose, not code, so there is no line-excerpt to copy structurally — but two concrete content sources are load-bearing and must be transcribed, not recomputed:

**Set 1 (PRIMARY, pre-registered, all FAIL) — transcribe verbatim from `05-04-SUMMARY.md`** (RESEARCH "Go/No-Go Memo — Verbatim Numbers"). NOT reproducible from current artifacts (iterations overwrote the nets) — cite from the frozen summary text, never from `*_report.jld2`.
- SBC: ρ_true KS p=2.81e-40, χ² p=5.27e-146, ECE 0.189 (red); Δρ ECE 0.192 (red); 6/8 params green.
- BF: corr 0.8861, max|Δ logBF| 3.598 (both gates FAIL); log_prior_odds −0.0294.
- OOD: pooled AUC 0.730 (FAIL); noise family AUC 0.0 (named blind spot); ID fire-rate 0.05.

**Set 2 (POST-HOC, confirmatory) — from on-disk `*_report.jld2` + STATE.md `[Phase 5-iter1/iter2]`.** Frame explicitly as post-hoc (D-03/D-04): SBC ρ_true ECE 0.0164 green; BF corr 0.9358 (clamped-KDE tail artifact); OOD combined AUC 1.0. Name the three mitigations: `consts.jl` byte-unchanged vs `e9c91d3`, DEV seed `0xDE7C0DE`, single confirmatory VAL run each.

**Memo required content (DEMO-03, from CONTEXT D-01..D-08):** Clean Go verdict; both number sets in order (set1 primary → set2 post-hoc); explicit falsification condition (D-02); data-snooping exposure named + mitigations (D-04); independent confirmation ship-gate (D-05); Phase 8–16 build-out DAG with Phase 11 ∥ 13 lead, Phase 12 follows 11 (D-07/D-08); SBC "calibrated under the simulator, paired with OOD" framing (carried-forward hard requirement). Length 2–3 pages.

**Location:** `.planning/phases/06-reproducible-demo-go-no-go-memo/06-GO-NO-GO-MEMO.md` (RESEARCH Open Q1 recommendation — keeps the `commit_docs` git flow consistent and `src/` untouched). `spike/GO_NO_GO.md` is an equally-decoupled alternative; planner's call.

---

## Shared Patterns

### Reproducible seeding (Philox4x, never `Random.seed!`)
**Source:** `spike/data/seeding.jl:75-76` (`sample_rng`) + `spike/validation/harness.jl:93-94` (`val_rng`)
**Apply to:** the `demo.jl` fast-tier chain (use `val_rng(VAL_FIX_SEED)`).
```julia
# seeding.jl:75-76
sample_rng(master_seed::Integer, idx::Integer) = Philox4x(UInt64, (UInt64(master_seed), UInt64(idx)))
# harness.jl:93-94 — the validation stream (XOR-salted, disjoint from training/holdout/fold)
val_rng(master_seed = VAL_MASTER_SEED) = Philox4x(UInt64, (UInt64(master_seed) ⊻ VAL_SALT, UInt64(0)))
```
Discipline (`consts.jl:76-79`): `VAL_MASTER_SEED=0x5BC0FFEE` is the RESERVED reported stream; `VAL_FIX_SEED=0xF1F7ED` is the fixture stream. Fast tier MUST use the fixture seed so it never pre-observes the reported stream.

### Frozen-model load (single source, transforms never re-fit)
**Source:** `spike/validation/harness.jl:78-80` → `spike/npe/train_npe.jl:195-201` (`load_npe`) + `spike/validation/train_ratio.jl:356-364` (`load_ratio`)
**Apply to:** both demo tiers.
```julia
# harness.jl:78-80
load_frozen_model(path = joinpath(@__DIR__, "..", "npe", "trained_npe.jld2")) = load_npe(path)
# load_npe returns (estimator, θzt, zt, variant, d_in, meta); load_ratio returns (estimator, zt, log_prior_odds, num_summaries, meta)
```
Both loaders integrity-check the `"estimator"` key on open and `error()` on a missing file — `demo.jl` inherits that guard for the frozen nets (which ARE committed / present; blocking if deleted).

### CPU-only (`use_gpu = false` everywhere)
**Source:** `spike/00_smoke.jl:66-67`, `spike/validation/harness.jl:110` (`posterior_for(...; use_gpu=false)`)
**Apply to:** every NeuralEstimators call reached by `demo.jl` (all go through the harness, which already passes `use_gpu=false`). Never load CUDA; `test_sbc.jl:114-116` asserts no CUDA module is loaded.

### Anti-snooping guard assertions
**Source:** `spike/validation/run_sbc.jl:74-75` (also `run_bf.jl:78-79`)
**Apply to:** any `--full` path in `demo.jl` (inherited when it shells out to `run_*.jl`).
```julia
@assert VAL_MASTER_SEED != NPE_MASTER_SEED "reported stream must be disjoint from training stream (D-02)"
@assert VAL_MASTER_SEED != VAL_FIX_SEED    "reported stream must be disjoint from fixture stream (D-02)"
```

---

## No Analog Found

None. Both new files have strong existing analogs (the entry-script idiom is exact; the memo composites `05-04-SUMMARY.md` + `NOTES.md`). No new libraries, no unprecedented data flow.

---

## Conventions

Repo-wide automated derivation **skipped** (`gsd-tools verify conventions --derive --scope spike` returned `{skipped: true, reason: "no-readable-files"}` — the deterministic module targets JS/TS/Python identifier axes and does not parse Julia `.jl` sources). Conventions below are derived by direct majority-vote observation across the 8 spike files read this session (`02_simulator_demo.jl`, `00_smoke.jl`, `seeding.jl`, `harness.jl`, `consts.jl`, `run_sbc.jl`, `test_sbc.jl`, `run_thread_sweep.jl`).

| Axis | Dominant | Share | Entropy | Status |
|------|----------|-------|---------|--------|
| File-name casing | `snake_case.jl`, numeric-prefixed entry scripts (`00_smoke.jl`, `02_simulator_demo.jl`) | ~100% | low | named contract |
| Identifier casing | `snake_case` functions/locals; `SCREAMING_SNAKE` module-level `const`; Greek/Unicode where domain-meaningful (`ρ_true`, `θzt`, `Δρ`) | ~95% | low | named contract |
| Export / module style | Flat top-level functions, NO `module` wrapper; guarded `isdefined(@__MODULE__, :sym) \|\| include(...)` for idempotent loading | ~100% (spike/) | low | named contract |
| Import / include style | `using Pkg` at top; intra-spike deps via `include(joinpath(@__DIR__, ...))` with explicit path segments; `src/` reached ONLY through read-only `contract.jl` | ~100% | low | named contract |

Additional named contracts (deviation is a regression, not author's choice):
- **AGPL header** — identical 19-line `#= … =#` block opens every `.jl` file (verbatim copy).
- **Purpose banner** — a `# spike/<file>.jl --- <one-line purpose>` comment block with decoupling + PowerShell-safe-verify notes follows the header.
- **Self-assert close** — scripts end with `@assert isfile(...)`/`@test` + a `println("… OK: …")` line; the verify command is the bare `julia --project=spike <script>` (no POSIX `test -f`).
- **Pre-registration lock** — all thresholds/seeds live in `consts.jl`, guarded as one `if !isdefined(@__MODULE__, :SBC_M)` block; `demo.jl` reads them, never redefines.

**Contested hotspots (author's choice):** None material within `spike/`. The spike is a single internally-consistent CJS-analog subtree — flat-function, no-module, guarded-include style throughout. (The plugin's own `bin/lib/**` CJS vs `sdk/src/**` ESM dual-resolver split is the prototype intentional-contested case; no equivalent split exists here — `spike/` is uniformly one style, and `demo.jl` should match it exactly.) The one hard boundary is directional, not stylistic: `spike/` code may `include` `src/` read-only via `contract.jl` but must never write it (the decoupling proof enforces this).

---

## Metadata

**Analog search scope:** `spike/` (entry scripts, `validation/`, `data/`, `npe/`, `test/`), `.planning/phases/*/` (SUMMARY precedent).
**Files scanned/read:** `02_simulator_demo.jl`, `00_smoke.jl`, `data/seeding.jl`, `validation/harness.jl`, `validation/consts.jl`, `validation/run_sbc.jl`, `validation/run_bf.jl` (head), `validation/run_ood.jl` (head), `test/test_sbc.jl`, `npe/run_thread_sweep.jl` (tail), `npe/train_npe.jl` (`load_npe`), `validation/train_ratio.jl` (`load_ratio`), `spike/NOTES.md` (head).
**Verified this session:** `git status --porcelain -- src/ …` empty (SC2 satisfiable); report keys; header/idiom line numbers.
**Pattern extraction date:** 2026-07-03
