# Phase 9: Cross-Method Comparator Harness - Research

**Researched:** 2026-07-02
**Domain:** Classical fluorescence-colocalization estimators (Costes-p, Manders M1/M2, Pearson, Spearman), a Tapqir/CoSMoS sanity bridge, and a seeded/reproducible comparison-table harness — all spike-local, reusing the BayesInteractomics comparator/audit pattern.
**Confidence:** HIGH (in-repo reuse surfaces, seeding, cache, table pattern), MEDIUM (Costes scheme parameters, Tapqir tutorial specifics)

## Summary

Phase 9 is almost entirely an **integration/assembly** phase over code that already exists in this repo and in the sibling BayesInteractomics repo. Three of the four classical estimators (Manders M1/M2, whole-image Pearson) are already implemented inline in `spike/data/encode.jl:encode_aug`; Pearson/Spearman over the patch grid are already in the frozen `src/colocalization.jl:correlation`; the shared-input generator (`simulate_pair` + `build_mci`) with ground-truth θ, the Random123/Philox seeding, and the content-addressed JLD2 cache are all shipped and green from Phases 2–4. The genuinely **new** code is small and well-scoped: a Costes randomization p-value, a divergence/traffic-light column, a tidy DataFrame→CSV/JLD2 table writer, an optional Tapqir bridge, and a seeded entry point — all under a new `spike/comparator/` module.

Two findings materially shape the plan. **First**, `DataFrames` and `CSV` are present in `spike/Manifest.toml` only as *transitive* deps (pulled by NeuralEstimators); they are **not** in `spike/Project.toml [deps]`. Using them directly requires a Wave-0 env task that promotes them to direct deps, re-freezes the Manifest, and extends the resolve-risk gate — exactly the pattern Phase 3 used to promote JLD2. **Second, and critically:** `spike/data/encode.jl` is a member of the Phase-3 cache content-hash list (`spike/data/hashguard.jl:HASH_SRC_FILES`, line 48). Editing `encode.jl` at all — even a byte-neutral refactor to "factor Manders into a shared callable" — flips `cache_hash` and **auto-invalidates the entire ≥50k-sample Phase-3 training cache and the trained Phase-4 NPE**. The literal reading of D-04 ("factor into a callable rather than copy-paste") therefore collides with a completed, expensive artifact. The safe resolution (below) keeps `encode.jl` byte-identical and enforces D-04's *intent* (consistency with `encode_aug`) via a shared new-module callable pinned by an equality test.

**Primary recommendation:** Build a new self-contained `spike/comparator/` module (estimators + Costes-p + table + divergence + seeded entry point) that consumes `Vector{MultiChannelImage}` read-only; do **not** touch `encode.jl` or `src/`; isolate the Tapqir bridge in its own sub-environment with skip-with-flag graceful degradation; and pin all "knows when classics are wrong" thresholds as pre-declared constants (D-14) mirroring the Phase-4/`test_npe.jl` pre-registration idiom.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

- **D-01:** Shared input container is the existing **`MultiChannelImage`** (2-channel), built via `spike/contract.jl:build_mci`, so every estimator AND the NPE consume byte-identical inputs (SC1).
- **D-02:** Primary shared input set is **simulator-generated** via `spike/simulator/forward.jl:simulate_pair` over a **seeded θ grid** spanning colocalization regimes (coloc / random / exclusion across the ρ_true range). The simulator's θ gives each row a **ground-truth regime label**.
- **D-03:** Keep the harness input-source-agnostic: accept any `Vector{MultiChannelImage}` so the Phase-8 external corpus / CBS can be fed later without redesign, but take **no Phase-8 dependency now**.
- **D-04:** **Reuse, do not reimplement**: Manders M1/M2 (Otsu from `mci.otsu_threshold`) from `encode_aug`, factored into a callable rather than copy-paste so comparator and augmented encoder stay consistent; Pearson/Spearman from `src/colocalization.jl:correlation` (and/or whole-image `cor`); report the **whole-image** scalar per method matching `encode_aug`'s `pearson_whole`.
- **D-05:** **Implement new: Costes-p** — the Costes randomization significance p-value (block/pixel-scramble null → fraction of scrambles with correlation ≥ observed). New code lives **only in the spike-local comparator module**; seed with the spike RNG.
- **D-06:** All estimator wrappers live in a new spike module, e.g. `spike/comparator/classical.jl`, each returning a NaN-safe scalar with the same degenerate-input guarding style as `encode_aug` (`_safe`, guarded reductions).
- **D-07:** The Tapqir bridge is a **minimal sanity anchor only** (SC2): reproduce a **published Tapqir example/tutorial dataset** and check the recovered quantity against the published value. **NOT** a comparator column on our simulator images (Tapqir models CoSMoS single-molecule spot data, a different regime).
- **D-08:** Integrate via **PythonCall/CondaPkg**, isolated in the spike. Must **degrade gracefully**: if the Python/Tapqir env is unavailable the harness **skips the Tapqir anchor with a clear flag** and the classical battery still runs green. Tapqir is not a hard gate on the classical table.
- **D-09:** Emit a **tidy per-input `DataFrame`**: one row per shared input, columns = `{input_id, regime/θ ground-truth, Costes_p, M1, M2, Pearson, Spearman, NPE_ρ̂ (or Δρ), NPE_OOD_flag}`. Persist as a versioned artifact (CSV for humans + JLD2 for exact reload), seeded and content-addressed like the Phase-3 cache.
- **D-10:** Add a **disagreement/divergence indicator** — a derived column/flag marking rows where a classical verdict diverges from the simulator ground-truth (and, where available, from the NPE). Present with the **BayesInteractomics traffic-light pattern** (green/amber/red), not a bare number.
- **D-11:** Reuse the **BayesInteractomics comparator/audit table + diagnostics pattern** (`_bin_calibration`/`CalibrationResult`/traffic-light thresholds; `_ks_test_uniform`/`model_diagnostics` reporting style) for table assembly and audit summary (SC3).
- **D-12:** Single seeded entry point `spike/comparator/run_comparator.jl` that builds shared inputs, runs all estimators + optional Tapqir anchor, writes the table. Seed with **Random123 (Philox)** exactly as the rest of the spike; table is bit-reproducible across runs and thread counts.
- **D-13:** A spike test (wired into `spike/test`) asserts: (a) determinism of the table across two seeded runs; (b) Costes-p, M1, M2, Pearson, Spearman all emit **finite** values on a non-degenerate fixture (SC1); (c) the Tapqir bridge either reproduces the published anchor within tolerance OR is cleanly skipped-with-flag when env is absent (SC2).
- **D-14:** The "knows when classics are wrong" thresholds (divergence cutoff, traffic-light bands) must be **fixed/pre-declared** in the harness config, not tuned post-hoc. Document chosen thresholds + rationale in the artifact header.

### Claude's Discretion

- Exact Costes randomization scheme (whole-image pixel scramble vs block/Van Steensel shift) and number of scrambles — pick the standard, cite it, keep it seeded.
- Table column ordering, CSV vs additional Arrow output, figure rendering (if any).
- Module/file layout under `spike/comparator/` and how estimators are registered.

### Deferred Ideas (OUT OF SCOPE)

- Run classical estimators + NPE on the Phase-8 physical corpus / CBS for the final head-to-head — **Phase 16**.
- Costes automatic threshold (Costes regression) determination beyond Otsu — keep Otsu for M1/M2 now.
- Full Tapqir integration as a live comparator column — explicitly rejected as invalid for our data regime (D-07).
- Manuscript prose framing of the "when classics are wrong" delta — **Phase 10**.
</user_constraints>

<phase_requirements>
## Phase Requirements

Phase 9 requirement IDs are TBD in ROADMAP/REQUIREMENTS. Proposed **CMP-*** IDs mapped to the 3 success criteria and the 14 locked decisions:

| ID | Description | Maps to | Research Support |
|----|-------------|---------|------------------|
| **CMP-01** | Classical battery (Costes-p, Manders M1/M2, Pearson, Spearman) runs on one shared `MultiChannelImage` input set and emits a per-method comparison table | SC1; D-01,D-04,D-06,D-09 | Reuse surfaces in `encode_aug`/`correlation`; table via DataFrame (§Standard Stack, §Code Examples) |
| **CMP-02** | Costes randomization significance p-value implemented new, spike-local, seeded (block-scramble null; p = fraction of scrambles with r ≥ r_obs) | SC1; D-05 | Costes et al. 2004 scheme (§Pattern 2, §Code Examples) |
| **CMP-03** | Shared inputs are simulator-generated over a seeded θ grid with ground-truth regime labels; harness signature accepts any `Vector{MultiChannelImage}` | SC1; D-02,D-03 | `simulate_pair`/`sample_prior`/`build_mci` (§Reusable Assets) |
| **CMP-04** | Divergence / traffic-light "knows when classics are wrong" column with **pre-declared** thresholds | SC1(positioning); D-10,D-14 | BayesInteractomics status idiom + Phase-4 pre-registration (§Pattern 3, §Pitfall 3) |
| **CMP-05** | Table persisted as tidy `DataFrame` + CSV + JLD2, content-addressed / seeded like the Phase-3 cache | SC1; D-09,D-12 | Phase-3 `cache.jl`/`hashguard.jl` (§Reusable Assets, §Pitfall 1) |
| **CMP-06** | Tapqir bridge reproduces a published Tapqir tutorial example within tolerance **OR** skips-with-flag when env absent | SC2; D-07,D-08 | Tapqir eLife 2022 + isolated env / PythonCall (§Pattern 4, §Environment Availability) |
| **CMP-07** | Reuse BayesInteractomics comparator/audit pattern (`CalibrationResult`/`_bin_calibration`/traffic-light; `_ks_test_uniform`/report style) for table assembly + audit summary | SC3; D-11 | BayesInteractomics diagnostics (§Reusable Assets, §Code Examples) |
| **CMP-08** | Single seeded entry point `run_comparator.jl`; Random123 Philox; bit-reproducible across runs and thread counts | SC3; D-12 | `seeding.jl` (§Reusable Assets) |
| **CMP-09** | Spike test gate asserts determinism, finiteness, and Tapqir anchor-or-skip | SC1,SC2; D-13 | `test_npe.jl`/`runtests.jl` scaffold idiom (§Validation Architecture) |

*The planner should propose these CMP-* IDs to append to REQUIREMENTS.md §Phase 9 and the traceability table.*
</phase_requirements>

## Project Constraints (from CLAUDE.md)

Actionable directives extracted from the project CLAUDE.md — treat with the same authority as locked decisions:

- **Decoupling (hard):** All new code under `spike/`, with the spike's own `Project.toml`/`Manifest.toml`. **No `src/` edits, no root `Project.toml`/`Manifest.toml` edits.** `src/` is consumed read-only via the existing `include()` boundary (`spike/contract.jl`).
- **Platform:** Windows-tauglich; optional external stacks (CUDA today, Python/Tapqir here) must **degrade gracefully to a working core path** — never a hard failure.
- **Reproducibility:** Everything reproducible from a fixed seed; the spike standard is **Random123 (Philox4x)**, keyed by `(master_seed, index)`, thread- and order-independent.
- **GSD workflow enforcement:** Do not make direct repo edits outside a GSD workflow.
- **No global-env installs of spike deps** (violates decoupling) — use `Pkg.activate("spike")` and its own manifest.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Shared-input generation (θ grid → images) | Simulator (`spike/simulator/`) | Seeding (`spike/data/seeding.jl`) | Already the spike's shared-input source with ground-truth θ (D-02) |
| Input container / summary contract | Contract (`spike/contract.jl`) | frozen `src/` (read-only) | `build_mci`/`patch_summary` are the frozen shared-input contract (D-01) |
| Classical estimator math | **New `spike/comparator/`** | reuse `encode_aug` formula + `correlation` | New callables; D-04 reuse-not-reimplement; must not edit hashed `encode.jl` |
| Costes randomization p-value | **New `spike/comparator/`** | `seeding.jl` (Philox) | Missing everywhere; spike-local + seeded (D-05) |
| Table assembly / divergence / traffic-light | **New `spike/comparator/`** | BayesInteractomics pattern (template) | D-09/D-10/D-11 |
| Persistence (CSV + JLD2, content hash) | **New `spike/comparator/`** | reuse Phase-3 `cache.jl` idiom | D-09/D-12 |
| Tapqir sanity anchor | **New isolated sub-env** (Python) | PythonCall/CondaPkg | Different data regime; optional, graceful-skip (D-07/D-08) |
| Reproducible seeding | `spike/data/seeding.jl` | Random123 | Reuse `sample_rng`/Philox (D-12) |
| Test gate | `spike/test/` | stdlib `Test` | D-13; plug into existing `runtests.jl` |

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| DataFrames.jl | 1.7.x (already resolved transitively) | Tidy per-input comparison table (D-09) | De-facto Julia tabular container; already the BayesInteractomics diagnostics table type; **must be promoted to a direct dep** (currently transitive via NeuralEstimators) `[CITED: spike/Manifest.toml]` |
| CSV.jl | 0.10.x (already resolved transitively) | Human-readable CSV export (D-09) | Standard DataFrame writer; **must be promoted to a direct dep** `[CITED: spike/Manifest.toml]` |
| JLD2.jl | already a direct dep | Exact-reload artifact + content-addressed dir (D-09/D-12) | Already the Phase-3 cache backend; Windows-clean, pure-Julia `[VERIFIED: spike/Project.toml]` |
| Random123.jl | already a direct dep | Philox4x seeded, thread-independent RNG (D-12) | Already the spike seeding primitive `[VERIFIED: spike/Project.toml]` |
| Statistics / StatsBase | already direct deps | `cor`, `corspearman`, median/IQR, NaN-safe reductions | Already used by `encode_aug`/`correlation` `[VERIFIED: spike/Project.toml]` |
| SHA (stdlib) | stdlib | Content hash for the artifact dir (mirror `hashguard.jl`) | Already used by Phase-3 cache; no `Pkg.add` `[VERIFIED: spike/data/hashguard.jl]` |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| PythonCall.jl + CondaPkg.jl | latest | Bridge to the Python Tapqir package (D-08) | Tapqir anchor only; **isolate in its own sub-environment** (see Pitfall 5) |
| CairoMakie.jl | already a direct dep | Optional divergence/traffic-light figure (Claude's discretion) | Only if a figure is wanted; headless spike-local pattern already established |
| HypothesisTests.jl | already a direct dep | Optional KS-uniformity in the audit summary (mirrors `_ks_test_uniform`) | If the audit reports a uniformity statistic |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| DataFrames.jl | Plain `NamedTuple` vectors + manual CSV | Loses the BayesInteractomics table/`leftjoin` idiom D-11 wants; DataFrames already resolved |
| PythonCall/CondaPkg for Tapqir | Subprocess call to a pre-existing Tapqir CLI env | Subprocess avoids any co-resolve risk entirely; viable graceful-skip alternative (see Open Questions) |
| Costes block scramble | Whole-image pixel scramble | Pixel scramble destroys spatial autocorrelation → over-optimistic (too-small) p-values; block scramble is the canonical Costes null (§Pattern 2) |
| Content-hash artifact dir | Timestamped dir | Loses bit-reproducibility / auto-invalidation; Phase-3 already proves the hash pattern |

**Installation (Wave-0 env task — mirror Phase-3 03-01):**

```bash
# In spike/ env only — never touch root Project.toml
julia --project=spike -e 'using Pkg; Pkg.add(["DataFrames","CSV"])'   # promote transitive→direct
# Tapqir bridge deps go in an ISOLATED sub-env (see Pitfall 5), NOT the main spike env:
#   spike/comparator/tapqir_env/  with its own Project.toml pinning PythonCall + CondaPkg
```

**Version verification:** DataFrames and CSV are already pinned in `spike/Manifest.toml` (resolved transitively via NeuralEstimators). Promoting them to direct deps should **not** change the resolved versions; the Wave-0 resolve-risk gate must re-assert NeuralEstimators stays at v0.2.1 (the recurring GLMakie-class downgrade landmine — see Pitfall 4).

## Package Legitimacy Audit

> This phase installs no packages from npm/PyPI on the core path. Julia deps (DataFrames, CSV, JLD2, Random123) are established General-registry packages already resolved in the pinned `spike/Manifest.toml`. slopcheck targets npm/PyPI/crates and is not the relevant gate for the Julia General registry; the relevant gate is the **resolve-risk test** (NeuralEstimators pin preserved) already in `runtests.jl`.

| Package | Registry | Age | Source Repo | Verdict | Disposition |
|---------|----------|-----|-------------|---------|-------------|
| DataFrames.jl | Julia General | ~9 yrs | github.com/JuliaData/DataFrames.jl | established | Approved (promote transitive→direct) |
| CSV.jl | Julia General | ~8 yrs | github.com/JuliaData/CSV.jl | established | Approved (promote transitive→direct) |
| PythonCall.jl / CondaPkg.jl | Julia General | ~4 yrs | github.com/JuliaPy/PythonCall.jl | established | Approved — **isolated sub-env only** |
| tapqir | PyPI | eLife-published 2022, Apache-2.0, v1.1.19 (2023-07-07) | github.com/gelles-brandeis/tapqir | established, peer-reviewed | Approved — optional, isolated, graceful-skip `[CITED: github.com/gelles-brandeis/tapqir]` |

**Packages removed due to slop verdict:** none.
**Packages flagged suspicious:** none. (slopcheck not run — Julia General registry ecosystem; the resolve-risk gate is the applicable integrity check.)

## Architecture Patterns

### System Architecture Diagram

```
                         master_seed (Random123 Philox)
                                    │
              ┌─────────────────────┼──────────────────────────┐
              ▼                                                 ▼
   sample_prior(rng)  ──►  θ grid (seeded, spans              (independent stream for
   (spike/simulator)       coloc/random/exclusion)             Costes scrambles per input)
              │            + ground-truth regime label
              ▼
   simulate_pair(rng,θ) ──► [ch1,ch2] ──► build_mci ──► Vector{MultiChannelImage}   ◄── (D-03: OR any external
   (spike/simulator)                      (contract.jl)     = SHARED INPUTS               Vector{MultiChannelImage})
                                                                  │
                    ┌─────────────────────────┬─────────────────┼───────────────────┬───────────────────┐
                    ▼                          ▼                 ▼                   ▼                   ▼
             manders(mci)              pearson_whole(mci)   spearman(mci)      costes_p(rng,mci)   NPE ρ̂/Δρ + OOD
             M1,M2 (Otsu from          cor(vec,vec)         patch-grid          NEW: block-           (Phase-4 nets,
             mci.otsu_threshold)       (= encode_aug)       correlation()       scramble null         where available)
                    └──────────────────────────┴─────────────────┴───────────────────┴───────────────────┘
                                                                  │
                                                                  ▼
                                              per-input row  (input_id, regime, Costes_p, M1, M2,
                                                              Pearson, Spearman, NPE_ρ̂, NPE_OOD_flag)
                                                                  │
                                                                  ▼
                             divergence + traffic-light column  (pre-declared thresholds, D-14)
                             "classical verdict vs ground-truth (and NPE)"  → green/amber/red
                                                                  │
                              ┌───────────────────────────────────┼───────────────────────────────┐
                              ▼                                     ▼                               ▼
                    tidy DataFrame                        CSV (humans)                    JLD2 (exact reload)
                                                        content-addressed dir (SHA-256 over sources+config)
                                                                  │
                    ┌─────────────────────────────────────────────┘
                    ▼
          audit summary (BayesInteractomics report style: counts per traffic-light band, KS-uniformity)

   ── SEPARATE, OPTIONAL BRANCH (never gates the table) ─────────────────────────────────────────
   tapqir_bridge()  ──►  isolated Python env available?  ──yes──► run published tutorial example ──► check recovered
   (skip-with-flag)                    │                          (CoSMoS spot data, NOT our images)   value ± tol
                                       └──no──► skip_flag=true, anchor_status=:skipped (classical table still green)
```

### Recommended Project Structure

```
spike/comparator/
├── config.jl            # PRE-DECLARED thresholds + constants (D-14): COSTES_N_SCRAMBLE,
│                        #   COSTES_BLOCK_PX, DIVERGENCE_CUTOFF, traffic-light bands, MASTER_SEED
├── classical.jl         # NaN-safe estimator callables: manders, pearson_whole, spearman_whole,
│                        #   patch_correlation; each takes a MultiChannelImage (D-04/D-06)
├── costes.jl            # NEW seeded Costes randomization p-value (D-05)
├── inputs.jl            # seeded θ-grid → Vector{MultiChannelImage} builder; regime labels (D-02);
│                        #   also accepts a caller-supplied Vector{MultiChannelImage} (D-03)
├── table.jl             # row assembly → DataFrame; divergence + traffic-light column (D-09/D-10);
│                        #   CSV + JLD2 writers; content-hash dir (mirror cache.jl/hashguard.jl)
├── audit.jl             # BayesInteractomics-style report: band counts, KS-uniformity (D-11)
├── tapqir_bridge.jl     # optional isolated Python bridge; skip-with-flag (D-07/D-08)
├── run_comparator.jl    # single seeded entry point (D-12)
└── tapqir_env/          # ISOLATED sub-env: own Project.toml/Manifest (PythonCall+CondaPkg)
spike/test/
└── test_comparator.jl   # D-13 gate, included from runtests.jl (mirror test_npe.jl scaffold)
```

### Pattern 1: Read-only reuse via the existing include() boundary

**What:** The comparator reaches `src/colocalization.jl:correlation`, `patch`, and `MultiChannelImage` only through `spike/contract.jl`'s already-established read-only `include()` — never a new `include` of `src/`, never an edit.
**When to use:** All Pearson/Spearman/patch-grid math and the `MultiChannelImage` type.
**Example:**
```julia
# spike/comparator/classical.jl
# contract.jl is ALREADY included by the harness entry point; it brings
# MultiChannelImage + patch + correlation (frozen src/) into scope read-only.
# Do NOT add a second include of src/colocalization.jl.
function patch_correlation(mci; method::Symbol = :spearman)
    x, y = mci.data[1], mci.data[2]
    xp, yp = patch.([x, y], 8)                       # src/colocalization.jl:37 (frozen)
    ρ = correlation(xp, yp; method = method)         # src/colocalization.jl:221 (frozen)
    vals = collect(skipmissing(ρ))
    return isempty(vals) ? NaN : Statistics.mean(vals)
end
```

### Pattern 2: Costes randomization significance p-value (NEW, D-05)

**What:** The canonical Costes et al. (2004) null: scramble one channel in **blocks the size of the resolution element (≈ PSF/spot size)** to preserve within-block autocorrelation, recompute the correlation coefficient, repeat N times, and report `p = (#{r_scramble ≥ r_obs} + 1) / (N + 1)`. Block scrambling (not pixel scrambling) is essential — scrambling individual pixels destroys spatial autocorrelation and yields falsely tiny p-values.
**When to use:** The Costes_p column of the table.
**Recommended concrete scheme (Claude's discretion, D-05):**
- **Statistic:** whole-image Pearson `cor(vec(ch1), vec(ch2))` (matches `encode_aug`'s `pearson_whole` so Costes-p and Pearson are the same statistic under H0 vs observed).
- **Null:** partition ch2 into square blocks of side `COSTES_BLOCK_PX`, randomly permute the block positions (a grid permutation), recompute r. Block side ≈ the resolution element; for this simulator the PSF is `σ_psf = 1.3 px`, so a **block side of ~7 px (≈ 4–6·σ_psf, one FWHM-ish resolution cell)** is a defensible, cite-anchored choice. Pre-declare it.
- **N:** `COSTES_N_SCRAMBLE = 200` — the number Costes et al. originally proposed for a stable null; larger N is more reliable but 200 is the standard `[CITED: Costes et al. 2004, Biophys J 86:3993]`.
- **Seeding:** draw block permutations from a Philox stream keyed by `(master_seed, input_index)` in a **disjoint salt namespace** (mirror `HOLDOUT_SALT`/`FOLD_SALT` in `seeding.jl`) so Costes-p is bit-reproducible and independent of the input-generation stream.
- **p-value:** the `+1/+1` (Davison–Hinkley) correction avoids a p=0 artifact and is standard for Monte-Carlo permutation tests.

### Pattern 3: Divergence + traffic-light with pre-declared thresholds (D-10/D-14)

**What:** A derived column marking rows where a classical verdict contradicts the ground-truth regime (and, where available, the NPE). Rendered as `:green`/`:amber`/`:red` via the exact BayesInteractomics status idiom: `status = x < warn ? :green : x < fail ? :amber : :red` with `warn`/`fail` **fixed in `config.jl` before any table is produced**.
**When to use:** The positioning payload — the whole reason for the phase.
**Example (mirrors `copula_diagnostics.jl:88`):**
```julia
# spike/comparator/config.jl  — pre-declared BEFORE any run (D-14), like test_npe.jl consts
const DIVERGENCE_WARN = 0.30   # |classical_verdict − ground_truth_regime_score| amber above this
const DIVERGENCE_FAIL = 0.60   # red above this
# spike/comparator/table.jl
traffic_light(d) = d < DIVERGENCE_WARN ? :green : d < DIVERGENCE_FAIL ? :amber : :red
```

### Pattern 4: Optional external bridge with skip-with-flag (D-07/D-08)

**What:** The Tapqir anchor is a separate branch that never gates the classical table. The bridge probes availability first; on any failure (no Python, no isolated env, no Tapqir, import error) it returns `(status = :skipped, reason = ...)` and the harness continues green — identical in spirit to the CPU/CUDA graceful-degradation the spike already enforces.
**Example:**
```julia
# spike/comparator/tapqir_bridge.jl
function tapqir_anchor(; tol = TAPQIR_TOL)
    avail = try _tapqir_env_available() catch; false end
    avail || return (status = :skipped, reason = "Tapqir/Python env unavailable", value = missing)
    try
        recovered = _run_tapqir_tutorial()          # isolated env / PythonCall
        ok = abs(recovered - TAPQIR_PUBLISHED_VALUE) ≤ tol
        return (status = ok ? :passed : :failed, value = recovered)
    catch e
        return (status = :skipped, reason = sprint(showerror, e), value = missing)
    end
end
```

### Anti-Patterns to Avoid

- **Editing `encode.jl` to "factor out" Manders.** `encode.jl ∈ HASH_SRC_FILES` → invalidates the entire Phase-3 cache + Phase-4 NPE. Put the callable in the new module; enforce consistency by an equality test (Pitfall 1).
- **Whole-image pixel scramble for Costes.** Destroys autocorrelation → false significance. Use block scramble (Pattern 2).
- **Adding a second `include("../src/colocalization.jl")`.** Double-defines the frozen functions; reach them through `contract.jl` only.
- **Tapqir as a comparator column on our images.** Explicitly rejected (D-07): wrong data regime (CoSMoS single-molecule spots vs diffuse 2-channel fields).
- **Tuning divergence/traffic-light thresholds against the produced table.** Data-snooping (D-14); pre-declare in `config.jl`.
- **Adding PythonCall/CondaPkg to the main `spike/Project.toml`.** Co-resolve risk against the NeuralEstimators v0.2.1 pin (Pitfall 4/5); isolate.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Pearson/Spearman/patch grid | A new correlation loop | `src/colocalization.jl:correlation` via `contract.jl` | Frozen, tested, `_exclude_zero` + ≥15-px floor already handled |
| Manders M1/M2 | A fresh threshold+ratio impl | The exact `encode_aug` formula, re-expressed as a callable in the new module | D-04 consistency; same `mci.otsu_threshold` source (Specific Idea) |
| Content-addressed reproducible artifact | Timestamp/uuid dirs | Mirror `spike/data/cache.jl` + `hashguard.jl` (SHA-256 over sources+config) | Bit-reproducibility + auto-invalidation already proven in Phase 3 |
| Seeded, thread-independent RNG | `Random.seed!` | `Random123` Philox keyed by `(seed, idx)` with a disjoint salt | Order/thread independence (D-12); the spike standard |
| Table + traffic-light + audit report | Bespoke printing | BayesInteractomics `CalibrationResult`/`_bin_calibration`/`_ks_test_uniform`/report idiom | D-11; a proven template in the sibling repo |
| NaN-safe reductions on degenerate inputs | ad-hoc `isempty` checks scattered | The `_safe(f, v)` idiom from `encode.jl` | D-06; one guarded reduction style, matches `encode_aug` |
| Colocalization significance p-value | Simple shuffle/quadgk | Costes block-scramble (Pattern 2) | The field-standard null; the whole point of the Costes column |

**Key insight:** This phase is ~80% wiring existing, tested code and ~20% new (Costes-p, divergence column, table writer, Tapqir bridge). The risk is not algorithmic difficulty — it is (a) accidentally invalidating the completed Phase-3/4 artifacts by touching hashed source, and (b) letting the optional Python stack destabilize the pinned spike environment.

## Runtime State Inventory

> Phase 9 is additive (new module), but it *reads* content-hashed artifacts and adds an optional external stack, so the runtime-state check applies.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | Phase-3 JLD2 training cache is content-addressed by a SHA-256 that **includes `spike/data/encode.jl` and the two frozen `src/` files** (`hashguard.jl:HASH_SRC_FILES`). Any edit to those files silently invalidates the cache dir name → next generation run rebuilds ≥50k samples and the NPE must be retrained. | **Do not edit `encode.jl` or `src/`.** New callables go in `spike/comparator/`. Consistency with `encode_aug` enforced by an equality test, not a refactor. |
| Live service config | None — no external service holds Phase-9 state. | None (verified: no daemon/DB in the spike). |
| OS-registered state | None — no scheduled tasks/services. | None. |
| Secrets / env vars | Tapqir bridge may read a `CONDA`/`JULIA_CONDAPKG_*` env or a path to a pre-built Tapqir env. No secrets. | Document any env var the bridge honors; default to graceful skip if unset. |
| Build artifacts / installed packages | `spike/Manifest.toml` will change when DataFrames/CSV are promoted to direct deps (re-freeze). A new isolated `spike/comparator/tapqir_env/Manifest.toml` will be created. CondaPkg materializes a Conda env on disk (large, machine-local) — must be git-ignored. | Re-freeze `spike/Manifest.toml`; `.gitignore` the CondaPkg env dir; commit the isolated `tapqir_env/Project.toml` + `Manifest.toml` as the reproducibility artifact. |

**Nothing found in categories 2 & 3:** None — verified by inspection of the spike (no daemons, DBs, or OS registrations; the spike is a self-contained Julia project run via `julia --project=spike`).

## Common Pitfalls

### Pitfall 1: Refactoring `encode.jl` invalidates the Phase-3/4 cache (HIGH severity)

**What goes wrong:** D-04 says "factor [Manders] into a callable rather than copy-paste." Taken literally (edit `encode_aug` to call a shared `manders()`), this changes `encode.jl`'s bytes.
**Why it happens:** `spike/data/encode.jl` is entry 5 of `HASH_SRC_FILES` (`hashguard.jl:48`). `cache_hash` folds the raw bytes of every listed file; any change flips the digest that *names the cache directory*, so the completed ≥50k-sample cache and the trained NPE become "stale/absent" and would have to be regenerated (expensive, and Phases 3–4 are already `Complete`).
**How to avoid:** Keep `encode.jl` byte-identical. Put the shared callable (`manders(mci) -> (M1, M2)`) in the **new** `spike/comparator/classical.jl`, reproducing the exact `encode_aug` formula (same `mci.otsu_threshold` source). Enforce D-04's *intent* with a test asserting `manders(mci)` equals the M1/M2 that `encode_aug` computes for the same `mci` (extract via the known moment positions or a direct recompute) — consistency by test, not by shared source.
**Warning signs:** A plan task that says "edit `encode.jl`" or "refactor `encode_aug`"; a suddenly-different cache dir hash; a forced 50k regeneration in the plan.

### Pitfall 2: Costes pixel-scramble instead of block-scramble (MEDIUM)

**What goes wrong:** Scrambling individual pixels gives a null with no spatial autocorrelation, so almost any real image looks "significant" → p-values biased toward 0.
**Why it happens:** Pixel scramble is simpler to code.
**How to avoid:** Block scramble with block side ≈ the resolution element (Pattern 2); pre-declare `COSTES_BLOCK_PX`. Cite Costes et al. 2004.
**Warning signs:** Costes_p ≈ 0 for random-regime inputs where Pearson ≈ 0.

### Pitfall 3: Post-hoc threshold tuning (data-snooping) (MEDIUM)

**What goes wrong:** Choosing the divergence cutoff / traffic-light bands after seeing which rows you "want" to flag red.
**Why it happens:** The positioning goal ("classic is wrong here") tempts fitting thresholds to the story.
**How to avoid:** D-14 — declare `DIVERGENCE_WARN`/`FAIL` and Costes/α thresholds as `const`s in `config.jl` before the first table, documented in the artifact header. Mirror the `test_npe.jl` "pre-registered constants declared in Wave-0" idiom (STATE Blockers: the same hazard flagged for Phase 5).
**Warning signs:** Threshold constants edited in the same commit that reports table results.

### Pitfall 4: Resolve-risk — a new dep downgrades the NeuralEstimators v0.2.1 pin (MEDIUM-HIGH)

**What goes wrong:** Adding DataFrames/CSV (or, worse, PythonCall/CondaPkg) to `spike/Project.toml` triggers a co-resolve that caps NeuralEstimators below the pinned v0.2.1 — the exact GLMakie/Turing landmine documented in Phases 1–4 (`runtests.jl` checks (d)/(e)/(f)).
**Why it happens:** The Julia resolver silently downgrades to satisfy compat.
**How to avoid:** DataFrames/CSV are already resolved transitively at compatible versions, so promotion should be inert — but the Wave-0 resolve-risk gate must re-assert `Pkg.dependencies()[NeuralEstimators UUID].version == v"0.2.1"` after promotion (extend the existing test block). Keep PythonCall/CondaPkg **out of the main env** entirely (Pitfall 5).
**Warning signs:** The resolve-risk testset fails after `Pkg.add`; NeuralEstimators resolves to 0.1.x.

### Pitfall 5: PythonCall/CondaPkg destabilizing the spike env or hard-failing on Windows (MEDIUM)

**What goes wrong:** Putting the Python stack in the main spike env risks co-resolve downgrades and makes the *classical* battery depend on a heavy, Windows-fragile Conda materialization.
**Why it happens:** Convenience of a single environment.
**How to avoid:** Isolate exactly as `spike/baseline/` isolates Turing/AdvancedVI (own `Project.toml`/`Manifest.toml`, D-01 precedent). Invoke the Tapqir bridge either by activating the isolated sub-env or via a subprocess; probe availability and skip-with-flag on any failure (D-08). The classical table must run green with zero Python present.
**Warning signs:** `test_comparator.jl` cannot pass without a Conda env; classical results blocked by a Python import error.

### Pitfall 6: Column-major / vec() ordering mismatch with the NPE (LOW)

**What goes wrong:** The optional NPE_ρ̂ column is computed from a summary vector whose element order differs from the trained net's expected `encode_d01` layout.
**Why it happens:** Re-deriving the summary instead of reusing `patch_summary`/`encode_d01`.
**How to avoid:** If including the NPE column, build its input via the *same* `patch_summary`→`encode_d01` path (rows 1:64 vals col-major, 65:128 mask) the net was trained on; reuse `spike/npe/infer.jl:posterior_for`.
**Warning signs:** NPE_ρ̂ nonsensical vs regime label on inputs where classicals agree.

## Code Examples

### Manders M1/M2 as a callable (reproduces `encode_aug`, NEW module, does NOT edit encode.jl)

```julia
# spike/comparator/classical.jl  — byte-for-byte the encode_aug formula (encode.jl:102-106)
_safe(f, v::AbstractVector) = isempty(v) ? NaN : f(v)   # mirror encode.jl:80

"manders(mci) -> (M1, M2) using the SAME Otsu thresholds as encode_aug (mci.otsu_threshold)."
function manders(mci)
    ch1, ch2 = mci.data[1], mci.data[2]
    t1, t2   = mci.otsu_threshold[1], mci.otsu_threshold[2]
    s1, s2   = sum(ch1), sum(ch2)
    M1 = s1 > 0 ? sum(ch1 .* (ch2 .> t2)) / s1 : NaN
    M2 = s2 > 0 ? sum(ch2 .* (ch1 .> t1)) / s2 : NaN
    return (M1, M2)
end

"whole-image Pearson (= encode_aug's pearson_whole, encode.jl:109)."
pearson_whole(mci) = cor(vec(mci.data[1]), vec(mci.data[2]))
```

### Costes block-scramble p-value (NEW, seeded)

```julia
# spike/comparator/costes.jl
using Random123
const COSTES_SALT = 0xA5A5A5A5DEADBEEF   # disjoint stream (mirror seeding.jl salts)
costes_rng(master_seed, idx) = Philox4x(UInt64, (UInt64(master_seed) ⊻ COSTES_SALT, UInt64(idx)))

"""
    costes_p(master_seed, idx, mci; n=COSTES_N_SCRAMBLE, block=COSTES_BLOCK_PX) -> Float64

Costes randomization significance: block-scramble ch2, recompute whole-image Pearson,
p = (#{r_scr ≥ r_obs} + 1)/(n + 1). Block side ≈ resolution element preserves autocorrelation.
"""
function costes_p(master_seed, idx, mci; n = COSTES_N_SCRAMBLE, block = COSTES_BLOCK_PX)
    ch1, ch2 = mci.data[1], mci.data[2]
    r_obs = cor(vec(ch1), vec(ch2))
    isfinite(r_obs) || return NaN
    rng = costes_rng(master_seed, idx)
    ge = 0
    for _ in 1:n
        scr = _block_permute(rng, ch2, block)      # grid-block position permutation
        r = cor(vec(ch1), vec(scr))
        (isfinite(r) && r ≥ r_obs) && (ge += 1)
    end
    return (ge + 1) / (n + 1)                       # Davison–Hinkley +1/+1
end
```

### Traffic-light + audit counts (BayesInteractomics idiom)

```julia
# spike/comparator/table.jl  (status idiom from copula_diagnostics.jl:88)
traffic_light(d) = d < DIVERGENCE_WARN ? :green : d < DIVERGENCE_FAIL ? :amber : :red
# audit.jl: band counts + KS-uniformity on p-values, printed as a Markdown table
# mirroring generate_diagnostics_report (predictive_checks.jl:886) and _ks_test_uniform (…:729).
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Costes automatic threshold (regression on 2D histogram) for M1/M2 | Otsu per channel (project choice) | D-04 / Deferred | Keep Otsu now; Costes-regression thresholds deferred |
| Van Steensel CCF block-shift significance | Costes block-randomization p-value | Costes et al. 2004 superseded/complemented van Steensel | We implement Costes; block scramble is the standard null `[CITED: imagej.net colocalization-analysis]` |
| Binary "colocalized?" decisions | Probabilistic per-spot output (Tapqir/cosmos) | Ordabayev et al. eLife 2022 | Tapqir is the modern SBI-adjacent baseline we anchor against `[CITED: elifesciences.org/articles/73860]` |

**Deprecated/outdated:** Pixel-level (non-block) randomization for Costes significance — known to over-report significance; not used.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | `COSTES_N_SCRAMBLE = 200` and block side ≈ 7 px (≈ resolution element for σ_psf=1.3) are appropriate for this simulator | Pattern 2 | Too-few scrambles → coarse p; wrong block size → biased p. Mitigated by pre-declaring and citing; tune only via a documented pre-registered choice. |
| A2 | DataFrames/CSV promotion to direct deps does not change resolved versions or downgrade NeuralEstimators | Standard Stack | If it downgrades, the resolve-risk gate fails — caught in Wave-0, not silently. |
| A3 | The Tapqir bridge should run against the **eLife 2022 tutorial "sample data set"** (Gelles-lab CoSMoS example distributed with Tapqir Part II) and check a recovered scalar (e.g. average target-specific spot probability / a `cosmos` model parameter) within a tolerance | Pattern 4, Environment | The exact dataset name + published numeric value + tolerance are **unconfirmed**; must be pinned by actually running the tutorial once. Until then the bridge should assert "reproduces documented tutorial output" generically. `[ASSUMED]` |
| A4 | The optional NPE_ρ̂ / OOD columns (D-09) can reuse the Phase-4 `infer.jl` surface directly on comparator inputs | Pitfall 6 | If the trained net expects a specific standardization frozen from training folds, the column needs the loader's frozen stats; may require the Phase-3 standardization artifact. |
| A5 | Isolating PythonCall/CondaPkg in a sub-env (like `spike/baseline/`) is sufficient for graceful degradation on Windows | Pitfall 5 | If PythonCall must load in the main process to marshal data, a subprocess-based bridge may be needed instead. |

## Open Questions (RESOLVED)

1. **Exact Tapqir tutorial dataset + published anchor value (SC2).**
   - What we know: Tapqir (Ordabayev et al., eLife 2022; `gelles-brandeis/tapqir`, Apache-2.0, v1.1.19) ships a Part II tutorial that runs the `cosmos` model on a "sample data set"; the model outputs per-frame/per-location target-specific spot probabilities and global parameters.
   - What's unclear: the precise sample-dataset identifier, which scalar to check, and its published value + acceptable tolerance.
   - Recommendation: In the first planning wave, run the tutorial once in the isolated env to capture the recovered scalar; pin it as `TAPQIR_PUBLISHED_VALUE`/`TAPQIR_TOL` in `config.jl` (pre-declared). Until captured, the D-13 test accepts "anchor passed within tol OR skipped-with-flag."
   - RESOLVED: Plan 09-04 Task 2 attempts the tutorial run to capture-and-pin `TAPQIR_PUBLISHED_VALUE`/`TAPQIR_TOL` in `config.jl`, or leaves a documented `NaN` clean-skip when Conda is unavailable; the 09-06 D-13 gate accepts `status ∈ (:passed,:skipped)`.

2. **Bridge mechanism: in-process PythonCall vs subprocess.**
   - What we know: BayesInteractomics has a CondaPkg precedent (project stack notes); PythonCall works in-process.
   - What's unclear: whether in-process PythonCall in the main harness risks env contamination on Windows.
   - Recommendation: Prefer a **subprocess** invocation of a script run under the isolated `tapqir_env` (or a pre-built Conda env), returning the recovered scalar via stdout/JSON — fully decoupled and trivially skip-with-flag. Keep in-process PythonCall as fallback.
   - RESOLVED: Plan 09-04 Task 1 mandates subprocess-only invocation under the isolated `tapqir_env` (no in-process PythonCall in the main harness); the acceptance criteria assert subprocess isolation.

3. **Does the divergence column require the NPE, or is ground-truth-only sufficient for SC1?**
   - What we know: D-10 says diverge "from the simulator ground-truth (and, where available, from the NPE)."
   - Recommendation: Make the NPE columns **optional** (present when the Phase-4 net + frozen standardization are available; `missing` otherwise). SC1 and the core positioning are satisfiable from ground-truth alone; the NPE comparison is an enhancement, keeping Phase 9 independent of Phase-4 artifact availability.
   - RESOLVED: Plan 09-05 Task 1 makes the NPE columns optional (`missing` when no NPE artifact is supplied); divergence is computed from ground-truth regime alone, keeping Phase 9 independent of Phase-4 artifact availability.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Julia + `spike/` env | Whole harness | ✓ | Julia 1.12.6 pinned | — |
| DataFrames.jl / CSV.jl | Table (D-09) | ✓ (transitive; promote to direct) | resolved in Manifest | Plain NamedTuple + manual CSV |
| JLD2 / Random123 / StatsBase / Statistics / SHA | Table, seeding, estimators, hash | ✓ (direct/stdlib) | pinned | — |
| Frozen `src/colocalization.jl`, `LoadImages.jl` | Pearson/Spearman/patch, MultiChannelImage | ✓ (read-only via contract.jl) | baseline f581d95 | — |
| Phase-4 trained NPE + frozen standardization | Optional NPE_ρ̂/OOD columns | ? (present if Phase-4 artifacts on disk) | — | Emit `missing` for NPE columns |
| Python + Conda + Tapqir | Tapqir anchor (D-08) | ✗ (not installed; optional) | tapqir 1.1.19 (target) | **Skip-with-flag** (D-08) — core path stays green |

**Missing dependencies with no fallback:** none — the classical battery + table have no hard external dependency.
**Missing dependencies with fallback:** Python/Tapqir → skip-with-flag; Phase-4 NPE artifacts → `missing` NPE columns.

## Validation Architecture

> `workflow.nyquist_validation` config key not located; treated as **enabled** (absent = enabled).

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Julia stdlib `Test` (`@testset`/`@test`) |
| Config file | none — single entry `spike/test/runtests.jl` |
| Quick run command | `julia --project=spike spike/test/runtests.jl` |
| Full suite command | `julia --project=spike spike/test/runtests.jl` (all phase testsets included) |

New testset `spike/test/test_comparator.jl` is `include`d from `runtests.jl` (append after `test_npe.jl`), mirroring the Phase-2/3/4 scaffold-then-fill idiom. Wave-0 lands named `@test_skip` placeholders for every SC; later waves replace them with real gates.

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| CMP-01 | Costes_p, M1, M2, Pearson, Spearman all **finite** on a non-degenerate fixture (D-13b) | unit | `julia --project=spike spike/test/runtests.jl` (test_comparator) | ❌ Wave 0 |
| CMP-02 | Costes-p uses block scramble; p bit-reproducible under fixed seed; p∈[1/(n+1),1] | unit | same | ❌ Wave 0 |
| CMP-02 | `manders(mci)` **equals** the M1/M2 `encode_aug` computes for the same mci (D-04 consistency, no encode.jl edit) | unit/consistency | same | ❌ Wave 0 |
| CMP-03 | Seeded θ grid → `Vector{MultiChannelImage}` reproducible; regime labels attached; harness accepts external vector | unit | same | ❌ Wave 0 |
| CMP-04 | Traffic-light bands use pre-declared consts; `traffic_light` monotone; red only above `DIVERGENCE_FAIL` | unit | same | ❌ Wave 0 |
| CMP-05/08 | **Determinism:** two seeded `run_comparator` runs produce byte-identical table (JLD2) + identical content-hash dir (D-13a); thread-count independence | integration | same | ❌ Wave 0 |
| CMP-06 | Tapqir anchor: reproduces published value ± tol **OR** `status=:skipped` cleanly when env absent (D-13c) | integration | same | ❌ Wave 0 |
| CMP-07 | Audit summary emits band counts + KS-uniformity in the BayesInteractomics report style | unit | same | ❌ Wave 0 |
| — | **Resolve-risk gate:** NeuralEstimators stays v0.2.1 after DataFrames/CSV promotion; Turing/PythonCall absent from main env | env | same (extend existing block d/e/f) | ⚠️ extend `runtests.jl` |

### Sampling Rate
- **Per task commit:** `julia --project=spike spike/test/runtests.jl` (fast: small fixture, small N for Costes in tests).
- **Per wave merge:** full suite green.
- **Phase gate:** full suite green + a real `run_comparator` produces the CSV/JLD2 artifact with a stable content hash before `/gsd:verify-work`.

### Wave 0 Gaps
- [ ] `spike/test/test_comparator.jl` — SC1/SC2/SC3 skipped scaffolds + pre-registered `const`s (mirror `test_npe.jl` Wave-0).
- [ ] `spike/comparator/config.jl` — pre-declared thresholds (`COSTES_N_SCRAMBLE`, `COSTES_BLOCK_PX`, `DIVERGENCE_WARN/FAIL`, `MASTER_SEED`, `TAPQIR_TOL`) committed BEFORE any reported run (D-14).
- [ ] Extend `runtests.jl` resolve-risk block to re-assert NeuralEstimators v0.2.1 after DataFrames/CSV promotion and to assert PythonCall/CondaPkg absent from the main env.
- [ ] Env: promote DataFrames + CSV to direct deps; re-freeze `spike/Manifest.toml`.
- [ ] `.gitignore` the CondaPkg materialized env dir; commit isolated `tapqir_env/Project.toml`+`Manifest.toml`.

## Security Domain

> `security_enforcement` config key not located; this is a local, offline scientific spike with no network service, no auth, no untrusted user input. ASVS categories are largely N/A. The one relevant surface:

| Concern | Applies | Standard Control |
|---------|---------|------------------|
| Untrusted parameter/input validation | yes (mild) | `simulate_pair` already validates θ ranges + imsize; comparator estimators guard degenerate inputs (`_safe`, NaN-not-throw) — continue that discipline |
| Supply-chain (new deps) | yes | Resolve-risk gate + pinned Manifest; Tapqir is Apache-2.0 peer-reviewed; isolate the Python env |
| Code execution via Python bridge | yes (mild) | Run Tapqir only in the isolated env / subprocess; never `eval` untrusted content; graceful-skip on failure |

No V2/V3/V4 (auth/session/access-control) surface — offline single-user research tool.

## Sources

### Primary (HIGH confidence)
- In-repo (read this session): `spike/data/encode.jl` (Manders/Pearson formula + `_safe`), `spike/contract.jl` (`build_mci`/`patch_summary`), `spike/simulator/forward.jl` (`simulate_pair`) + `prior.jl` (`sample_prior`), `spike/data/seeding.jl` (Philox salts), `spike/data/cache.jl` + `hashguard.jl` (content-hash cache; **`encode.jl ∈ HASH_SRC_FILES`**), `src/colocalization.jl:correlation` (frozen), `spike/test/runtests.jl` + `test_npe.jl` (scaffold + resolve-risk gate + pre-registered consts), `spike/Project.toml` + `Manifest.toml` (DataFrames/CSV transitive).
- Sibling repo (read this session): BayesInteractomics `src/diagnostics/{calibration.jl,types.jl,predictive_checks.jl,copula_diagnostics.jl}` — `_bin_calibration`/`CalibrationResult`, `_ks_test_uniform`, `generate_diagnostics_report`, and the `:pass/:warn/:fail` traffic-light status idiom.

### Secondary (MEDIUM confidence)
- Costes et al. 2004, "Automatic and Quantitative Measurement of Protein-Protein Colocalization in Live Cells," Biophys J 86:3993–4003 — block randomization, ~200 rounds, P>95% down to ~3% coloc `[CITED: sciencedirect.com/science/article/pii/S0006349504744392]`.
- ImageJ colocalization analysis + JACoP docs — Costes block scramble vs van Steensel block shift `[CITED: imagej.net/imaging/colocalization-analysis]`.
- Ordabayev et al. 2022, eLife 11:e73860 — Tapqir/cosmos CoSMoS model, probabilistic per-spot output `[CITED: elifesciences.org/articles/73860]`; `github.com/gelles-brandeis/tapqir` (Apache-2.0, v1.1.19, 2023-07-07).

### Tertiary (LOW confidence)
- Tapqir tutorial specifics (exact sample dataset name, published anchor value, tolerance) — not confirmed; captured as A3/Open Question 1 (`tapqir.readthedocs.io` Part II referenced but dataset details not extracted).

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all Julia deps confirmed in the pinned Manifest; DataFrames/CSV promotion path proven by the Phase-3 JLD2 precedent.
- Architecture / reuse surfaces: HIGH — exact signatures and the cache-hash constraint read directly from source this session.
- Costes scheme parameters: MEDIUM — method is canonical and cited; block size / N are defensible pre-declared choices (A1).
- Tapqir anchor specifics: LOW-MEDIUM — package + paper confirmed; exact tutorial dataset/value unconfirmed (A3, OQ1).

**Research date:** 2026-07-02
**Valid until:** 2026-08-01 (stable in-repo surfaces; Tapqir tutorial specifics should be confirmed by running it once during planning)
