# Phase 6: Reproducible Demo + Go/No-Go Memo - Research

**Researched:** 2026-07-03
**Domain:** Reproducible Julia entry-script authoring (spike/demo.jl) + scientific decision-memo authoring; no new libraries, no net retraining
**Confidence:** HIGH (all findings verified against on-disk artifacts and existing spike code this session)

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **D-01: Memo verdict = Clean Go** — recommend full productionization (Phase 7). The memo must be explicit that the literal pre-registered gates did **not** pass as written; the Clean Go is justified by reading each residual failure to a characterized, non-method cause (χ²-over-power at M=2000, residual data-scale gap, clamped-KDE baseline artifact), not by claiming a pass. Backing evidence to cite: core thesis proven (Phase 4 amortized >100× faster than ADVI at ρ recovery corr 0.983); OOD PASS (pooled AUC 1.0 post-iteration); SBC ECE green on all 8 params (ρ_true 0.016, was 0.19 red); BF near-miss corr 0.936 with diagnosed clamped-baseline cause.
- **D-02: State the falsification condition explicitly** — what result would have forced a No-Go (SBC ECE staying red after capacity+data iteration, OR BF mid-range disagreement, OR an OOD family undetectable by any summary-orthogonal channel).
- **D-03: Report BOTH number sets, clearly labeled, in this order:** (1) the **original pre-registered gate run** at locked consts on the fresh `VAL_MASTER_SEED` (all three gates FAIL) as the *primary pre-registered result*; then (2) the **post-iteration confirmatory numbers** from the two bounded retrains, explicitly framed as *post-hoc*.
- **D-04: Name the data-snooping exposure plainly** (frozen net retrained twice *after* the pre-registered result was seen) + the three mitigations: `consts.jl` byte-unchanged (verified vs commit `e9c91d3`), model selection on a disjoint DEV seed (`0xDE7C0DE`), single confirmatory VAL run per iteration (no retry-to-pass). Do NOT present post-iteration numbers as if they were the pre-registered result.
- **D-05: Independent confirmation ship-gate** — the Clean Go commits Phase 7 to a fresh-seed, re-pre-registered SBC/BF/OOD confirmation run (new locked consts, new disjoint seed, single run) that must reproduce the calibration/OOD/BF story before amortized inference ships into `src/`. This is a **hard ship-gate**, not a nice-to-have.
- **D-06: Two-tier demo.jl** — fast default (`julia --project=spike spike/demo.jl`): loads frozen artifacts (`spike/npe/trained_npe.jld2`, `spike/validation/trained_ratio.jld2`) + reported `*_report.jld2` for headline numbers, chains every layer end-to-end at fixture scale (small M/L, short Δρ sweep, OOD subset) to prove the pipeline reproduces, tabulates the success-criteria table. A `--full` flag re-runs the reported-scale pipeline (M=2000 SBC, full BF sweep, full OOD grid). The decoupling proof (SC2) runs inside `demo.jl` as a `git status --porcelain -- src/ src/bayes.jl src/colocalization.jl` assertion.
- **D-07: Memo commits to the existing Phase 8–16 roadmap DAG** as full-build-out direction; feature work is not blocked behind a calibration close-out — the one hard precondition is the D-05 ship-gate (inside Phase 7).
- **D-08: Prioritize Phase 11 (registration/chromatic uncertainty) ∥ Phase 13 (three-hypothesis amortized BF) as lead axes; Phase 12 (spatial coloc map) follows 11.** Hierarchy/3D/multi-channel are the longer-horizon build-out SC3 asks for, sequenced after the Wave-B differentiators.

### Claude's Discretion
- Exact layout and section order of the memo (within the 2–3 page bound), figure selection from the reported `figures/` set, and the demo's success-criteria table formatting.
- Whether the decoupling proof also greps the two manuscript pipeline dirs by name or relies on the `src/`-clean assertion plus a whole-tree `git status` — planner/researcher's call. **(Research recommendation below: rely on the `src/`-clean + whole-tree assertion; no separate manuscript-pipeline directories exist to grep — see Finding G1.)**

### Deferred Ideas (OUT OF SCOPE)
- The independent confirmation ship-gate is *named* by the memo (D-05) but **executed in Phase 7** — Phase 6 reports and decides, it does not retrain or re-run gates.
- Re-enabling the OOD posterior-predictive channel in the reported OR-fusion (`run_ood.jl` still runs `with_pp=false`) — a Phase-7 hardening item.
- RxInfer independent cross-check (BACK-01) — deferred post-Go; the memo may cite it as a "nice-to-have for the paper."
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| DEMO-01 | `spike/demo.jl` chains the full pipeline reproducibly from a fixed seed and tabulates success criteria | Entry-script idiom fixed by `spike/02_simulator_demo.jl` + `spike/00_smoke.jl`; two-tier fixture/reported split fixed by Phase-5 `test_*.jl` vs `run_*.jl`; seeding via `sample_rng`/`val_rng` (Philox4x). All frozen artifacts + reports verified present on disk (Findings A1–A2). CLI `--full` pattern established (Finding C1). |
| DEMO-02 | Main repo and both manuscript pipelines are demonstrably untouched (decoupling proof) | `git status --porcelain` shell-out from Julia; `src/` currently clean (verified). No separate manuscript-pipeline dirs exist — decoupling reduces to `src/`-clean + whole-tree assertion (Finding G1). |
| DEMO-03 | 2–3-page Go/No-Go memo reports metrics + a concrete full-build-out decision | Original pre-registered FAIL numbers verbatim in `05-04-SUMMARY.md` (Finding B1); post-iteration confirmatory numbers in STATE.md + reproducible from on-disk `*_report.jld2` (Finding B2). Figure inventory in `spike/validation/figures/` verified (Finding B3). Build-out DAG in ROADMAP §Progress (D-07/D-08). |
</phase_requirements>

## Summary

Phase 6 is a **documentation + single-script** phase. It ships exactly three deliverables — `spike/demo.jl` (a two-tier seeded end-to-end runner), an in-script decoupling assertion, and a 2–3-page Go/No-Go memo — and it deliberately does **not** retrain nets or re-run gates (Phase 5 is closed; the memo *reports* its results). No new Julia packages are introduced; every dependency (`JLD2`, `CairoMakie`, `Random123`, `Test`) is already in the pinned `spike/Manifest.toml`. Because nothing is installed, there is no package-legitimacy audit surface.

The single most consequential finding for planning is a **provenance split in the on-disk artifacts**: the `*_report.jld2` files under `spike/validation/` currently hold the **post-iteration confirmatory numbers** (verified this session: SBC ρ_true ECE=0.0164 green, BF corr=0.9358, OOD combined AUC=1.0), **not** the original pre-registered FAIL run. The original pre-registered numbers (SBC ρ_true KS p=2.81e-40 / ECE 0.189 red, BF corr 0.886 / max|Δ|=3.60, OOD pooled AUC 0.730) survive **only as text** in `05-04-SUMMARY.md` because the iterations overwrote `trained_npe.jld2` (committed `200e971`) and `trained_ratio.jld2`. This means D-03's *primary* number set (set 1) is **not reproducible from current artifacts** — `demo.jl --full` would regenerate set 2, not set 1. The memo must therefore cite set 1 from the frozen `05-04-SUMMARY.md` record, and the demo's success-criteria table must be explicit that its loaded headline numbers are the post-iteration (set 2) values.

**Primary recommendation:** Build `demo.jl` as a top-level AGPL-headered CPU-only script mirroring `02_simulator_demo.jl`, gated on `"--full" in ARGS`; in fast mode chain prior→simulator→summary→NPE/NRE→SBC/BF/OOD at fixture scale (using the `test_*.jl` `VAL_FIX_SEED` streams, not the reserved `VAL_MASTER_SEED`), load the three `*_report.jld2` for headline numbers, print a success-criteria table, and run the `git status` decoupling assertion. Author the memo to report set 1 (from `05-04-SUMMARY.md`) as primary and set 2 (from the on-disk reports / STATE.md) as post-hoc, with the D-04 mitigations and D-05 ship-gate.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| End-to-end pipeline chaining (prior→sim→summary→NPE/NRE→SBC/BF/OOD) | Spike script (`spike/demo.jl`) | Existing spike modules (harness/sbc/bf/ood) | demo.jl *composes* the frozen read surfaces; it re-implements nothing (mirrors how `run_*.jl` drive the machinery) |
| Reproducible seeding | `spike/data/seeding.jl` + `harness.jl` | — | `sample_rng`/`val_rng` (Philox4x) are the canonical primitives; demo.jl threads a fixed seed through these, never `Random.seed!` |
| Frozen-artifact / report loading | `JLD2` (fast tier) | — | Headline numbers loaded from `*_report.jld2`; nets from `.jld2` via `load_npe`/`load_ratio` |
| Decoupling proof (git status) | OS shell-out from Julia | — | Read-only certification; no repo tier is mutated |
| Go/No-Go memo | Documentation (`.planning/` or `spike/`) | — | A prose artifact; it reports, it does not compute |
| Figure selection | Reuse `spike/validation/figures/*.png` | CairoMakie (only if a new composite panel is wanted) | Reported-scale figures already exist; regeneration is optional |

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| JLD2 | pinned in `spike/Manifest.toml` | Load frozen nets + `*_report.jld2` headline numbers | Already the spike serialization choice; `jldopen(f,"r")` read pattern used by every `run_*.jl` |
| Random123 | v1.7.1 (pinned) | `Philox4x` counter-based seeding via `sample_rng`/`val_rng` | The reproducibility primitive (D-11); thread/order-independent |
| Test | stdlib | Hard-gate `@testset` idiom for the fixture-scale checks + `--full` reported gates | Every `run_*.jl` / `test_*.jl` uses it; exits nonzero on fail |
| CairoMakie | pinned | Optional: render a composite demo figure headless | `02_simulator_demo.jl` idiom (`CairoMakie.activate!()`); only if a new panel is desired |

**No packages are installed in this phase.** `spike/Project.toml` + `spike/Manifest.toml` stay byte-unchanged (verify with `git diff`). The whole point of the pinned Manifest is CPU-only reproducibility (SC1).

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Bare `"--full" in ARGS` flag parse | `ArgParse.jl` | ArgParse is **not** in the pinned Manifest — adding it violates the no-install constraint. Use the bare-ARGS idiom already present in `spike/npe/run_thread_sweep.jl` (`ARGS[1] == "--worker"`). |
| Reuse existing `figures/*.png` | Regenerate via CairoMakie | Regeneration costs minutes and risks nondeterministic font/raster diffs; the reported PNGs already exist and are the memo's evidence. Prefer reuse. |

**Installation:** none. Confirm environment is intact:
```bash
julia --project=spike -e 'using Pkg; Pkg.status()'   # must show NeuralEstimators 0.2.1, Flux 0.16.10, no CUDA
```

## Package Legitimacy Audit

**Not applicable — this phase installs no external packages.** `spike/Project.toml`/`spike/Manifest.toml` must remain byte-unchanged (hard decoupling constraint). Any planned task that adds a dependency is out of scope and should be rejected in plan-check.

## Architecture Patterns

### System Data-Flow Diagram (demo.jl fast tier)

```
                    julia --project=spike spike/demo.jl   [--full]
                                   │
                    ┌──────────────┴───────────────┐
                    │  parse ARGS → FULL = "--full"  │
                    └──────────────┬───────────────┘
                                   │
              ┌────────────────────┼─────────────────────────────┐
              │ (0) DECOUPLING ASSERTION                          │
              │  read(`git status --porcelain -- src/ …`) → @assert empty │
              └────────────────────┬─────────────────────────────┘
                                   │
        ┌──────────────────────────┼──────────────────────────┐
        │ (1) CHAIN PROOF  (fixed seed, fixture scale)          │
        │   sample_prior(rng) → simulate_pair → build_mci        │
        │     → patch_summary → encode_d01 → standardize(m.zt)   │
        │     → posterior_for (NPE)  ─┐                          │
        │     → amortized_log_bf (NRE)┤ one forward pass each    │
        │     → ood score channels    ┘                          │
        │   rng = val_rng(VAL_FIX_SEED)  (NOT the reserved stream)│
        └──────────────────────────┬──────────────────────────┘
                                   │
        ┌──────────────────────────┼──────────────────────────┐
        │ (2) HEADLINE NUMBERS                                   │
        │   fast:  jldopen(sbc_report/bf_report/ood_report.jld2) │
        │   --full: run_sbc.jl / run_bf.jl / run_ood.jl (reported)│
        └──────────────────────────┬──────────────────────────┘
                                   │
        ┌──────────────────────────┼──────────────────────────┐
        │ (3) SUCCESS-CRITERIA TABLE → stdout                    │
        │   SC1 reproduces?  SC2 decoupled?  SC3 memo exists?     │
        │   + DEMO-01/02/03 pass/fail  → self-assert, exit code   │
        └───────────────────────────────────────────────────────┘
```

### Pattern 1: Top-level seeded CPU-only entry script
**What:** AGPL header → `using` block → `CairoMakie.activate!()` (if plotting) → guarded `include()` of the composition surfaces → fixed-seed config consts → compute → `save`/print → `@assert isfile(...)` / `@test` self-assert at the end. PowerShell-safe verification: the script self-asserts (no POSIX `test -f`); the verify command is just `julia --project=spike spike/demo.jl`.
**When to use:** the whole of `demo.jl`.
**Example (idiom from `spike/02_simulator_demo.jl`):**
```julia
# Source: spike/02_simulator_demo.jl (verified this session)
using CairoMakie, StatsBase, Statistics, Random, Distributions
CairoMakie.activate!()                       # headless, no display
include(joinpath(@__DIR__, "contract.jl"))               # build_mci, patch_summary
include(joinpath(@__DIR__, "simulator", "forward.jl"))   # simulate_pair
include(joinpath(@__DIR__, "simulator", "prior.jl"))     # sample_prior
# ... compute ...
@assert isfile(PNG_PATH) "demo FAILED: missing $PNG_PATH"
println("demo OK: ...")
```

### Pattern 2: Compose, don't re-implement (drive the frozen read surfaces)
**What:** demo.jl should `include` and call `harness.jl` (`load_frozen_model`, `draw_simulate_infer`, `val_rng`), `sbc.jl`, `bf.jl`, `ood.jl` — exactly as `run_*.jl` do — rather than re-deriving inference. For `--full`, `include` the three `run_*.jl` (each `main()` is called at file scope) or shell out to them as subprocesses.
**When to use:** both tiers.
**Example (composition surface from `harness.jl`, verified):**
```julia
# Source: spike/validation/harness.jl (verified this session)
m = load_frozen_model()                       # loads spike/npe/trained_npe.jld2 ONCE
r = draw_simulate_infer(m, val_rng(VAL_FIX_SEED); imsize=(64,64), N=99)
# r.θ (true), r.draws (7×N physical-θ posterior)  — CPU-only, frozen m.zt/m.θzt
```

### Pattern 3: Bare-ARGS flag (no ArgParse dependency)
**What:** `const FULL = "--full" in ARGS`. Optionally also gate script-vs-include with `if abspath(PROGRAM_FILE) == @__FILE__`.
**When to use:** top of `demo.jl`.
**Example (idiom from `spike/npe/run_thread_sweep.jl`, verified):**
```julia
# Source: spike/npe/run_thread_sweep.jl:215 (verified this session)
if abspath(PROGRAM_FILE) == (@__FILE__) && !isempty(ARGS) && ARGS[1] == "--worker"
```
Adapt to: `const FULL = ("--full" in ARGS)`.

### Pattern 4: git-status decoupling assertion (SC2)
**What:** shell out to git with **explicit argv (no shell string interpolation)**, capture stdout, assert empty.
**When to use:** step 0 of demo.jl (fail fast before any compute).
**Example:**
```julia
# Recommended pattern — argv form avoids any shell-injection surface
out = read(`git status --porcelain -- src/ src/bayes.jl src/colocalization.jl`, String)
@assert isempty(strip(out)) "DECOUPLING VIOLATED: src/ is dirty:\n$out"
# Optional whole-tree read-only sanity (excluding .planning/ churn) — see Finding G1
```
Julia backtick command literals pass argv directly to the process (no shell), so there is no command-injection risk here — keep the arguments as separate literal tokens, never build the command from an interpolated string.

### Anti-Patterns to Avoid
- **Consuming the reserved `VAL_MASTER_SEED` in the fast tier.** The fast fixture chain must use `VAL_FIX_SEED` (via `val_rng(VAL_FIX_SEED)`) so it never pre-observes the reported stream — exactly the discipline `test_sbc.jl` follows. Only `--full` reported gates touch `VAL_MASTER_SEED`.
- **Re-fitting standardization or retraining.** demo.jl loads `m.zt`/`m.θzt` frozen via `load_frozen_model`; it never calls `fit(ZScoreTransform,...)`. No net is trained (Phase 5 closed).
- **Editing `consts.jl`, `Project.toml`, `Manifest.toml`, or anything under `src/`.** All are hard-locked. demo.jl reaches `src/` only through read-only `include()` via `contract.jl`.
- **Adding `ArgParse` or any package** to satisfy `--full` parsing. Use bare `ARGS`.
- **Presenting the on-disk report numbers as the pre-registered result.** They are the post-iteration (set 2) numbers (Finding B2). The memo's primary set 1 comes from `05-04-SUMMARY.md`.

### Recommended Project Structure (delta only)
```
spike/
├── demo.jl              # NEW — the Phase-6 deliverable (top-level, mirrors 02_simulator_demo.jl)
├── npe/trained_npe.jld2         # frozen NPE (committed 200e971 = iter1 retrain)  [fast-tier load]
├── validation/
│   ├── trained_ratio.jld2       # frozen NRE (gitignored, present, iter2)          [fast-tier load]
│   ├── sbc_report.jld2          # gitignored, present — POST-ITERATION numbers      [fast-tier load]
│   ├── bf_report.jld2           # gitignored, present — POST-ITERATION numbers      [fast-tier load]
│   ├── ood_report.jld2          # gitignored, present — POST-ITERATION numbers      [fast-tier load]
│   ├── run_sbc.jl / run_bf.jl / run_ood.jl   # reported gates driven by --full
│   └── figures/*.png            # reported-scale figures for the memo
.planning/phases/06-.../
└── 06-GO-NO-GO-MEMO.md   # NEW — or spike/GO_NO_GO.md; the 2–3-page memo (planner's call on location)
```

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Reproducible seeding | A new `Random.seed!` scheme | `sample_rng` / `val_rng` (Philox4x) from `seeding.jl`/`harness.jl` | Thread/order-independent, already the project contract; reinventing risks a stream that overlaps the training seed |
| θ→sim→summary→infer step | New inference glue in demo.jl | `draw_simulate_infer` / `draw_simulate_infer_paired` (`harness.jl`) | Freezes `m.zt`/`m.θzt` correctly (Pitfall 5); one CPU-only forward pass |
| Reported SBC/BF/OOD run | Re-derive M=2000 loop | `run_sbc.jl` / `run_bf.jl` / `run_ood.jl` under `--full` | Each already owns its artifact, figures, and a real pre-registered gate |
| CLI flag parsing | ArgParse dependency | `"--full" in ARGS` | No install allowed; bare-ARGS idiom already in the repo |
| Headline number retrieval | Recompute from scratch in fast mode | `jldopen(*_report.jld2, "r")` | Reports already persisted; keys verified below |
| Decoupling check | Diff files by hand | `read(\`git status --porcelain -- src/ …\`, String)` + `@assert isempty` | One line, authoritative, matches D-06 |

**Key insight:** Every capability demo.jl needs already exists as a frozen, tested surface. Phase 6 is 90% orchestration + prose, ~10% new code.

## Runtime State Inventory

> demo.jl loads pre-existing runtime state; this inventory maps what the fast tier depends on and its provenance.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data (frozen nets) | `spike/npe/trained_npe.jld2` (4.9 MB, **git-committed** `200e971` = iter1 retrain); `spike/validation/trained_ratio.jld2` (1.1 MB, **gitignored but present**, iter2 retrain) | demo.jl loads read-only via `load_npe`/`load_ratio`. **Note:** committed net = iter1, NOT the original Phase-4 net. |
| Stored data (headline reports) | `sbc_report.jld2` (140 KB), `bf_report.jld2` (2.7 KB), `ood_report.jld2` (21 KB) — all **gitignored but present**; all hold **POST-ITERATION confirmatory** numbers (VERIFIED, Finding B2) | Fast tier loads these for the success-criteria table. Planner must ensure the memo does NOT mislabel these as pre-registered. |
| Stored data (BF sweep draws) | `bf_sweep_draws.jld2` (569 KB, gitignored, present) | Consumed by the BF baseline reproduction under `--full`; not needed by fast tier. |
| Live service config | None — fully offline, CPU-only, local artifacts | None. |
| OS-registered state | None. | None. |
| Secrets/env vars | None. | None. |
| Build artifacts | `spike/Manifest.toml` pins the env; must stay byte-unchanged. | Verify `git diff spike/Manifest.toml` empty at phase end. |

**Reproducibility caveat (VERIFIED):** the original pre-registered FAIL run is **not** reproducible from current on-disk artifacts — `trained_npe.jld2` (committed `200e971`) and `trained_ratio.jld2` were overwritten by the two bounded iterations. `demo.jl --full` regenerates the **post-iteration** (set 2) numbers, not set 1. Set 1 lives only as verbatim text in `05-04-SUMMARY.md`. This is the primary planning constraint for D-03.

## Common Pitfalls

### Pitfall 1: Treating the on-disk reports as the pre-registered result
**What goes wrong:** demo.jl / the memo loads `*_report.jld2` and labels those numbers as the pre-registered gate outcome. They are the **post-iteration** numbers (SBC ρ_true ECE=0.016 green, BF corr=0.936, OOD AUC=1.0).
**Why it happens:** the iterations overwrote the artifacts in place; filenames are unchanged.
**How to avoid:** memo set 1 = verbatim from `05-04-SUMMARY.md` (Finding B1); set 2 = on-disk reports / STATE.md, labeled post-hoc. demo's table must annotate which set its loaded numbers belong to.
**Warning signs:** any SBC ECE shown as green in a section labeled "pre-registered"; BF corr shown as 0.936 (not 0.886) under "primary result".

### Pitfall 2: Fast tier consuming the reserved stream
**What goes wrong:** demo.jl fixture chain uses `val_rng(VAL_MASTER_SEED)`, pre-observing the reported stream and eroding the anti-snooping contract.
**How to avoid:** use `VAL_FIX_SEED` in fast mode (mirror `test_sbc.jl`).
**Warning signs:** `VAL_MASTER_SEED` appearing outside a `--full`/reported path.

### Pitfall 3: Non-finite posterior θ̂ crashing the OOD/PP path
**What goes wrong:** on strongly-OOD inputs the frozen NPE can return non-finite θ̂; the PP re-simulation path throws ("all θ fields must be finite"). `run_ood.jl` runs `with_pp=false` for this reason; the iter1 finite-guard exists in `ood.jl` but the reported fusion still excludes PP.
**How to avoid:** demo.jl's OOD fixture path should mirror `run_ood.jl` (`with_pp=false`) — do not enable PP fusion. (Re-enabling PP is an explicit Phase-7 deferred item.)
**Warning signs:** a demo run aborting inside `simulate_pair` during the OOD leg.

### Pitfall 4: `git status` false-positive from `.planning/` churn
**What goes wrong:** a whole-tree `git status` assertion trips on unrelated `.planning/` edits (e.g. the current `.planning/HANDOFF.json` deletion + untracked `05-PATTERNS.md`) and reports a decoupling "violation" that is not a `src/` mutation.
**How to avoid:** scope the hard assertion to `-- src/ src/bayes.jl src/colocalization.jl` (this path-limited status is currently clean — VERIFIED). If a whole-tree sanity print is added, make it informational, not a gate, or path-exclude `.planning/`, `spike/`, `manuscript/`.
**Warning signs:** demo failing SC2 while `git status -- src/` is actually empty.

### Pitfall 5: Long fast-tier runtime defeating the "fast default" intent
**What goes wrong:** fixture-scale still uses 256×256 imsize or M=2000, making the "fast" default take minutes.
**How to avoid:** use fixture constants — `SBC_FIX_M`/`SBC_FIX_L` at 64×64 (per `test_sbc.jl`), a short Δρ sweep, an OOD subset. Keep the default in the seconds-to-low-minutes range; reserve reported scale for `--full`.

## Code Examples

### Loading headline numbers (fast tier) — VERIFIED report keys
```julia
# Source: verified by reading the on-disk artifacts this session (2026-07-03)
using JLD2
sbc = jldopen("spike/validation/sbc_report.jld2","r") do f
    (labels=f["labels"], ks_p=f["ks_p"], chi2_p=f["chi2_p"], ece=f["ece"],
     verdict=f["verdict"], caption=f["caption"], M=f["SBC_M"], L=f["SBC_L"])
end
bf  = jldopen("spike/validation/bf_report.jld2","r") do f
    (corr=f["corr"], max_abs_err=f["max_abs_err"],
     logbf_amortized=f["logbf_amortized"], logbf_kde=f["logbf_kde"], targets=f["targets"])
end
ood = jldopen("spike/validation/ood_report.jld2","r") do f
    (maha_auc=f["maha_auc"], combined_auc=f["combined_auc"], fam_best_auc=f["fam_best_auc"],
     id_fire_rate=f["id_fire_rate"], neg_ks=f["neg_ks"], pp_channel_note=f["pp_channel_note"])
end
```
**Verified report keys:**
- `sbc_report.jld2`: `rank_table, labels, ks_p, chi2_p, ece, mce, verdict, cov_nominal, cov_empirical, SBC_M, SBC_L, SBC_BINS, SBC_KS_ALPHA, SBC_CHI2_ALPHA, SBC_ECE_GREEN, SBC_ECE_YELLOW, SBC_IMSIZE, VAL_MASTER_SEED, caption, elapsed_min, generated`
- `bf_report.jld2`: `targets, logbf_amortized, logbf_kde, corr, max_abs_err, log_prior_odds, BF_CORR_MIN, BF_LOGBF_TOL, BF_SWEEP_LO, BF_SWEEP_HI, BF_SWEEP_N, VAL_MASTER_SEED, generated`
- `ood_report.jld2`: `families, levels, maha_auc, fam_best_auc, fire_rate, combined_auc, roc_fpr, roc_tpr, maha_thr, id_fire_rate, youden_j, youden_fpr, youden_tpr, pp_channel_viable, pp_channel_note, neg_ks, neg_ks_max, neg_fire_rate, OOD_AUC_MIN, OOD_KS_EPS, OOD_ID_QUANTILE, OOD_GRID_LEVELS, OOD_PP_REPS, OOD_IMSIZE, n_id, n_pos, VAL_MASTER_SEED, generated`

### `--full` dispatch to reported gates
```julia
const FULL = ("--full" in ARGS)
if FULL
    # each run_*.jl calls main() at file scope and is a hard gate (exits nonzero on fail)
    include(joinpath(@__DIR__, "validation", "run_sbc.jl"))
    include(joinpath(@__DIR__, "validation", "run_bf.jl"))
    include(joinpath(@__DIR__, "validation", "run_ood.jl"))
else
    # fast: load *_report.jld2 (above) + run a fixture-scale chain proof on VAL_FIX_SEED
end
```

## Go/No-Go Memo — Verbatim Numbers (for DEMO-03)

### Set 1 — Original pre-registered gate run (PRIMARY, all three FAIL)
Source: `05-04-SUMMARY.md` (`.planning/phases/05-.../05-04-SUMMARY.md`). These are **not** reproducible from current artifacts — cite from the frozen summary. `[VERIFIED: 05-04-SUMMARY.md read this session]`

**SBC (M=2000, L=999, bins=50, imsize 256×256, VAL_MASTER_SEED) — FAIL every parameter:**
| parameter | KS p | χ² p | ECE | verdict |
|-----------|------|------|-----|---------|
| ρ_true | 2.81e-40 | 5.27e-146 | 0.189 | red |
| spillover | 3.30e-20 | 1.54e-50 | 0.0459 | green |
| autofluorescence | 2.42e-15 | 3.51e-38 | 0.0469 | green |
| label_efficiency | 1.35e-6 | 7.54e-36 | 0.0302 | green |
| shift_dx | 1.03e-9 | 2.77e-33 | 0.0499 | green |
| shift_dy | 3.80e-12 | 1.89e-49 | 0.068 | yellow |
| noise | 2.81e-4 | ~0 (2.96e-370) | 0.0408 | green |
| Δρ | 9.49e-40 | 1.45e-147 | 0.192 | red |

**BF:** reported-scale ratio net (n=3000, ep=100): corr(amortized, KDE) = **0.8861**, max|Δ logBF| = **3.598** → both D-08a (≥0.95) and D-08b (≤0.5) FAIL. `log_prior_odds = −0.0294`.

**OOD:** density (Mahalanobis) best-per-family AUC — texture 1.0 (PASS), optics 0.983 (PASS), background 0.96 (PASS), **noise 0.0 (FAIL — named blind spot)**; **pooled AUC = 0.730 < 0.80 → FAIL**; held-out ID fire-rate = 0.05 (honest ~5% FPR). Negative controls: affine provably blind & quiet (correct); rotate/block density-flag fires ~93% (position-aware residual). PP channel reported **non-viable** (non-finite θ̂ out-of-distribution).

### Set 2 — Post-iteration confirmatory numbers (POST-HOC, D-03/D-04)
Sources: STATE.md `[Phase 5-iter1]`/`[Phase 5-iter2]` + on-disk `*_report.jld2` (VERIFIED this session, Finding B2). Frame explicitly as post-hoc; consts.jl byte-unchanged vs `e9c91d3`; DEV seed `0xDE7C0DE`; single confirmatory VAL run each.

- **SBC (iter1 retrain — higher-capacity flow, same 50k cache):** ECE **green on all 8 params** (ρ_true ECE=0.0164; was 0.189 red). KS now passes 3/8 (spillover p=0.068, autofluorescence p=0.196, noise p=0.395); ρ_true KS p=0.00323, Δρ KS p=0.00518 still reject. Overall strict all-8 KS∧χ² conjunction still FAILS (χ² hyper-sensitive at M=2000; label_efficiency/shift genuinely non-uniform). `[VERIFIED: sbc_report.jld2]`
- **BF (iter2 retrain — NRE on same 50k cache + A7 difference encoding):** corr **0.9358** (↑ from 0.886; DEV 0.862→0.903), max|Δ logBF| **10.95** (STATE) — the on-disk `bf_report.jld2` confirms corr=0.9358. Large max|Δ| diagnosed as a **clamped-KDE-baseline tail artifact** (baseline hits its `_clampp=1e-8` floor → logBF ±18.41 at |Δρ|≳0.4 while the bounded NRE stays ~±8; mid-range agrees within ~1–2). `[VERIFIED: bf_report.jld2 corr; CITED: STATE.md iter2 for max|Δ|]`
- **OOD (iter2 — added image-noise Channel 3, OR-fused with density):** all 4 families best-AUC **1.0**, **combined pooled AUC = 1.0**, ID fire-rate 0.05, negative controls KS-invariant + density-quiet. Remaining blind spot narrowed to transforms orthogonal to BOTH the 8×8 correlation summary AND the image-noise features (e.g. pure positive-affine rescale). `[VERIFIED: ood_report.jld2 combined_auc=1.0, maha_auc]`

### Speed / accuracy thesis (restate per NPE-03)
Amortized inference >100× faster than per-dataset ADVI (median 325× on the forward-pass clock; ~16× on the full posterior-sample workload) at comparable RMSE; realized ADVI baseline ~0.5 s/pair. ρ recovery corr 0.983 (Phase 4). `[CITED: REQUIREMENTS.md NPE-03]`

### Figures available for the memo (`spike/validation/figures/`, VERIFIED present)
- Per-parameter SBC: `sbc_rank_hist_<param>.png`, `sbc_coverage_<param>.png`, `sbc_reliability_<param>.png` for all 8 (ρ_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise, Δρ) — these are the **iter1** (post-iteration) figures (mtime 15:51).
- BF: `bf_agreement.png` (iter2, mtime 16:51).
- OOD: `ood_roc_combined.png` (iter2 fused detector, mtime 17:00) and the older `ood_roc.png` (12:11, density-only original).
- Older single-panel SBC figures (`sbc_rank_hist.png`, `sbc_coverage.png`, `sbc_reliability.png`, mtime 12:11) predate the iterations.
- Also available: `spike/figures/plausibility.png` (SIM-04 simulator behavior).

### Full-build-out decision (D-07/D-08) — commit to the ROADMAP DAG
Source: ROADMAP §Progress "Execution Order" + STATE Roadmap Evolution. `[CITED: ROADMAP.md]`
- Waves downstream of the Phase-6 Go: **A {8, 9, 10}** (parallelizable now; 9 & 10 already complete), **B {11 ∥ 13, then 12}**, **C {14, 15}**, **D {16}**.
- Lead axes (D-08): **Phase 11** (registration + chromatic uncertainty as latent) ∥ **Phase 13** (three-hypothesis amortized BF) — the axes that most differentiate v2.0 from Tapqir/Costes/Manders per Phase-10 positioning. **Phase 12** (spatial coloc map, GP/CAR) follows 11 (shared PosteriorEstimator/training code — serialized to avoid a merge collision; descope-to-v2.1 candidate).
- Longer horizon (SC3 "hierarchy/3D/multi-channel"): sequenced after the Wave-B differentiators; explicitly Out-of-Scope for v2.0 per REQUIREMENTS §Out of Scope.

## State of the Art

Not applicable — no external technology moved. This phase uses the frozen, pinned spike stack (NeuralEstimators 0.2.1, Flux 0.16.10) established in Phase 1 and unchanged since. `[CITED: spike/Project.toml + STATE.md]`

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | BF confirmatory max\|Δ logBF\| = 10.95 (only corr=0.9358 was read directly from `bf_report.jld2`; 10.95 is from STATE.md prose — `max_abs_err` key exists but its value was not printed this session) | Memo Set 2 | Low — a headline-figure typo; the planner should have demo.jl print `bf_report["max_abs_err"]` to confirm before the memo quotes it |
| A2 | Memo file location (`.planning/phases/06-.../` vs `spike/`) is left to planner discretion; no locked location exists | Structure | Low — cosmetic; both keep `src/` untouched |
| A3 | The "two manuscript pipelines" of DEMO-02 have no distinct repo directories and equate to the `src/` package + (protected) new `manuscript/` tree | Finding G1 | Medium — if a hidden/legacy pipeline dir exists outside the tracked tree, the `src/`-only assertion would under-cover it; mitigated by the optional whole-tree sanity read |

**Note:** Set 1 numbers, Set 2 SBC/OOD numbers, all report keys, artifact presence, and `git status -- src/` cleanliness were **VERIFIED** directly this session — they are not assumptions.

## Open Questions

1. **Where should the memo live?**
   - What we know: it must not touch `src/`; `.planning/phases/06-.../06-GO-NO-GO-MEMO.md` and `spike/GO_NO_GO.md` both satisfy decoupling.
   - Recommendation: put it in `.planning/phases/06-reproducible-demo-go-no-go-memo/` alongside the phase artifacts (keeps `commit_docs` git flow consistent); optionally symlink/copy a rendered PDF later. Planner's call.

2. **Should `--full` `include` the `run_*.jl` in-process or spawn subprocesses?**
   - What we know: each `run_*.jl` calls `main()` at file scope and exits nonzero on gate fail via `@testset`. `include`-ing all three in one process shares loaded modules (faster) but a failing `@testset` throws rather than sets an exit code cleanly.
   - Recommendation: spawn each as a subprocess (`run(\`julia --project=spike spike/validation/run_sbc.jl\`)`) so each hard gate's exit code is captured independently and one failure doesn't abort the others' reporting. Fast tier stays in-process.

3. **Does the fast-tier chain need to assert numerical reproduction, or just "runs clean"?**
   - What we know: SC1 says "chains the full pipeline reproducibly." The cheapest honest proof is bit-identical re-draw on a fixed `VAL_FIX_SEED` (Philox4x is deterministic).
   - Recommendation: have demo.jl run the fixture chain twice on the same seed and `@assert` the two posterior draws are `==` (or `isapprox`), proving reproducibility without needing to match reported-scale numbers.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Julia | running demo.jl | ✓ (pinned) | 1.12.6 via juliaup override + `spike/.julia-version` | — |
| `spike/` env (NeuralEstimators/Flux/JLD2/Random123/CairoMakie/Test) | all of demo.jl | ✓ | pinned in `spike/Manifest.toml` (NeuralEstimators 0.2.1, Flux 0.16.10, no CUDA) | — |
| `git` (for SC2 decoupling assertion) | DEMO-02 | ✓ (repo is a git repo) | system git | If absent: fall back to a tracked-hash comparison of `src/` files, but git is present |
| Frozen nets `trained_npe.jld2` / `trained_ratio.jld2` | fast + full tiers | ✓ present | iter1 / iter2 | none — **blocking if deleted** |
| `*_report.jld2` (sbc/bf/ood) | fast-tier headline numbers | ✓ present (gitignored) | post-iteration | `--full` regenerates (post-iteration numbers only) |

**Missing dependencies with no fallback:** none — all fast-tier inputs verified present on disk this session.
**Missing dependencies with fallback:** none material.

**Caveat:** because the reports + `trained_ratio.jld2` are **gitignored**, a fresh clone will NOT have them. `demo.jl` should detect missing reports and, if absent, either instruct the user to run `--full` (regenerates post-iteration numbers) or degrade the fast tier gracefully. Plan a guarded `isfile(...)` check with a clear message rather than a raw load error.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Julia stdlib `Test` (`@testset`/`@test`) + `@assert` self-checks (`02_simulator_demo.jl` idiom) |
| Config file | none — `spike/test/runtests.jl` is the aggregate gate; demo.jl self-asserts |
| Quick run command | `julia --project=spike spike/demo.jl` (fast tier — must run clean, seconds-to-minutes) |
| Full suite command | `julia --project=spike spike/demo.jl --full` (re-runs reported-scale gates) |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| DEMO-01 | demo.jl chains prior→sim→summary→NPE/NRE→SBC/BF/OOD from a fixed seed and prints a success-criteria table | smoke + self-assert | `julia --project=spike spike/demo.jl` (exit 0) | ❌ Wave 0 — demo.jl to be created |
| DEMO-01 | fast-tier chain is bit-reproducible on `VAL_FIX_SEED` | unit (in-script `@assert`) | same run; twin-draw equality assert | ❌ Wave 0 |
| DEMO-02 | `git status --porcelain -- src/ src/bayes.jl src/colocalization.jl` is empty | in-script `@assert` | same run (step 0) | ✅ currently clean (VERIFIED) — assertion logic ❌ Wave 0 |
| DEMO-03 | memo exists, is 2–3 pages, reports BOTH number sets + verdict + falsification condition + build-out DAG | doc review + optional in-demo `@assert isfile(memo)` + keyword grep | manual review + `julia … spike/demo.jl` asserting the memo file is present | ❌ Wave 0 — memo to be authored |
| DEMO-01 (full) | reported gates re-run at locked consts | integration (hard gates, exit nonzero on fail) | `julia --project=spike spike/demo.jl --full` | ✅ run_*.jl exist; `--full` wiring ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** `julia --project=spike spike/demo.jl` (fast tier green + decoupling assert passes).
- **Per wave merge:** `julia --project=spike spike/test/runtests.jl` (existing aggregate gate stays green — demo.jl must not break it) + `git diff --exit-code spike/Project.toml spike/Manifest.toml spike/validation/consts.jl` (locked files unchanged).
- **Phase gate:** `julia --project=spike spike/demo.jl` clean AND `git status -- src/` empty AND memo present with both number sets, before `/gsd:verify-work`.

### Wave 0 Gaps
- [ ] `spike/demo.jl` — the deliverable (DEMO-01/02); does not yet exist.
- [ ] `.planning/phases/06-.../06-GO-NO-GO-MEMO.md` (or `spike/GO_NO_GO.md`) — the memo (DEMO-03); does not yet exist.
- [ ] A guarded `isfile()` check in demo.jl for the gitignored reports (fresh-clone robustness).
- [ ] (Optional) an in-demo memo-content assertion (grep for both number sets, "falsification", "Clean Go").
- Framework install: none — `Test` is stdlib, everything else pinned.

## Security Domain

`security_enforcement` is not set in `.planning/config.json` (treated as enabled), but this phase's surface is an **offline, CPU-only, local documentation + script** task with no auth, no network, no untrusted input, no persistence of secrets. Most ASVS categories are N/A.

### Applicable ASVS Categories
| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | — (no auth surface) |
| V3 Session Management | no | — |
| V4 Access Control | no | — |
| V5 Input Validation | minimal | The only external input is `ARGS` (`--full` flag) — exact-match compare, no interpolation |
| V6 Cryptography | no | — (no crypto; seeds are reproducibility keys, not secrets) |
| V12 Files/Resources | minimal | demo.jl writes only under `spike/` (figures/stdout); loads read-only artifacts |

### Known Threat Patterns for this task
| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Command injection via the git shell-out | Tampering/Elevation | Use Julia backtick argv form (`\`git status --porcelain -- src/ …\``) — no shell string interpolation; arguments are literal tokens (Pattern 4). Do NOT build the command from an interpolated string. |
| Accidental `src/` mutation defeating the decoupling proof | Tampering | The decoupling assertion IS the control; keep demo.jl's `src/` access read-only via `include()` through `contract.jl`. |
| Dependency tampering (adding a package) | Supply chain | No installs; `git diff` on `Project.toml`/`Manifest.toml` must be empty at phase gate. |

## Project Constraints (from CLAUDE.md)

- **Hard decoupling:** spike lives entirely under `spike/` with its own `Project.toml`/`Manifest.toml`; `src/` and both manuscript pipelines must remain provably untouched. demo.jl reaches `src/` only via read-only `include()`. — *This is DEMO-02's entire purpose.*
- **CPU-only baseline:** `use_gpu=false` everywhere; the spike must run without CUDA. demo.jl inherits this from the harness/infer surfaces.
- **Reproducibility:** everything reproducible from a fixed Random123 seed; use `sample_rng`/`val_rng` (Philox4x), never `Random.seed!`.
- **Windows-tauglich:** PowerShell-safe verification — the script self-asserts (`@assert`/`@test`), no POSIX `test -f`; verify command is just `julia --project=spike spike/demo.jl`.
- **Summary dimension fixed** (8×8) during the spike — demo.jl does not parameterize it.
- **GSD workflow enforcement:** file-changing tools must run through a GSD command (this phase is planned work — proceed via plan-phase/execute-phase).
- **No new dependencies** in the spike env for this phase (reinforced by the pinned-Manifest reproducibility constraint).

## Sources

### Primary (HIGH confidence — verified this session)
- `spike/02_simulator_demo.jl`, `spike/00_smoke.jl` — entry-script idiom (AGPL header, CairoMakie.activate!, guarded includes, self-assert).
- `spike/data/seeding.jl`, `spike/validation/harness.jl` — `sample_rng`/`val_rng` Philox4x seeding, `load_frozen_model`, `draw_simulate_infer(_paired)`.
- `spike/validation/consts.jl` — locked pre-registration (VAL_MASTER_SEED=0x5BC0FFEE, VAL_FIX_SEED=0xF1F7ED, all thresholds).
- `spike/validation/run_sbc.jl` (read in full), `run_bf.jl` (head) — reported-gate structure, atomic report save, `@testset` hard gate.
- `spike/test/test_sbc.jl` — fixture-scale pattern (SBC_FIX_M/L at 64×64 on VAL_FIX_SEED).
- `spike/npe/run_thread_sweep.jl:215` — bare-ARGS flag idiom; `abspath(PROGRAM_FILE)==@__FILE__` guard.
- On-disk `sbc_report.jld2` / `bf_report.jld2` / `ood_report.jld2` — keys + values read directly (Finding B2): SBC ρ_true ECE=0.0164, BF corr=0.9358, OOD combined_auc=1.0 (post-iteration).
- `.planning/phases/05-.../05-04-SUMMARY.md` — Set 1 verbatim pre-registered FAIL numbers (Finding B1).
- `git status --porcelain -- src/ …` → empty (SC2 currently satisfiable); `git ls-files` → no manuscript-pipeline dirs (Finding G1); `git log -- trained_npe.jld2` → committed at iter1 `200e971`.
- `.planning/config.json` — `nyquist_validation: true`, no `security_enforcement` key.

### Secondary (CITED)
- `.planning/STATE.md` `[Phase 5-iter1]`/`[Phase 5-iter2]` — Set 2 post-iteration confirmatory numbers + mitigations (consts.jl vs `e9c91d3`, DEV seed 0xDE7C0DE).
- `.planning/ROADMAP.md` §Progress / Execution Order — Phase 8–16 DAG, Wave A/B/C/D (D-07/D-08).
- `.planning/REQUIREMENTS.md` — DEMO-01/02/03, NPE-03 speed restatement, Out-of-Scope.

## Metadata

**Confidence breakdown:**
- Standard stack / no-install: HIGH — verified against pinned Project/Manifest and existing script imports.
- Architecture / demo.jl patterns: HIGH — direct idioms copied from `02_simulator_demo.jl`, `run_sbc.jl`, `harness.jl`, `run_thread_sweep.jl`.
- Number-set provenance (the critical finding): HIGH — report values read directly from disk and cross-checked against 05-04-SUMMARY.md and STATE.md.
- Manuscript-pipeline decoupling scope (G1/A3): MEDIUM — no tracked pipeline dirs found, but a legacy/untracked pipeline outside the tree cannot be fully excluded.

**Research date:** 2026-07-03
**Valid until:** 2026-08-02 (stable — no fast-moving external deps; only risk is further net iterations overwriting the on-disk reports, which would change Set 2)
