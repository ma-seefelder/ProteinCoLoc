# ProteinCoLoc `src/` Refactoring Proposal (consolidated)

**Status:** analysis complete — **proposals only, nothing implemented** (per user directive).
**Evidence:** 5-dimension analysis, every finding adversarially verified (46 agents); hot-path
claims micro-benchmarked. Full detail per finding in `00N-*/README.md`; index in `MANIFEST.md`.

> **Sequencing constraint:** CLAUDE.md forbids editing `src/` during the active v2.0 (AmortizedColoc)
> spike. Treat this as a **separate `src/`-hardening effort** — run it *before* resuming the spike or
> during the v2.0 **productionization** phase. Do not interleave with the decoupled spike. Several
> items (lift `@model`, document prior ranges, seed threading) directly *enable* v2.0, so doing them
> at the productionization boundary is natural.

---

## Tier 0 — Correctness bugs (fix first, independent of any refactor)

| # | Bug | Location | Fix | Effort |
|---|-----|----------|-----|--------|
| B1 | **Δρ mixes posterior μ_sample with prior μ_control** → wrong headline effect size in `result.txt`; inconsistent with the BF in the same row. Found by **4/5** dimensions. | `utils.jl:368` | `posterior_samples.:μ_sample .- posterior_samples.:μ_control` (or `delta_rho(posterior)`) | S |
| B2 | **Posterior path uses removed Turing API** (`ADVI(n,iter)`, `vi(m,alg)`, `DynamicPPL.syms`) → likely MethodError/UndefVarError on the installed stack (0.43–0.45). | `bayes.jl:323-329`, `bayes.jl:44` | Resolve Manifest mismatch; port to `vi(m, q0, iter; adtype=…)` and `res.q` (VIResult on 0.44/0.45). | M |
| B3 | **`convert_posterior_samples` magic `1:10`** mislabels 2 per-image local params; mislabeled cols escape to the posterior CSV (10 cols vs prior's 8). | `bayes.jl:42-50` | Select the 8 globals explicitly (`1:8` + named vector). Note CSV col count 10→8. | S |
| B4 | **`minmax_norm!` mutates shared image data** → inference silently depends on plot flags/order. | `plot.jl:149` | `minmax_norm!(copy(img_data[i]))` (or a non-mutating variant). | S |
| B5 | **Fisher-z docstring is false** (dead `#tanh`); trap if uncommented (would corrupt ν/σ/τ). *Not a live bug.* | `bayes.jl:40,52-53` | Delete dead line; correct docstring to "raw [-1,1] space, no transform." | S |

B1 is the single most important finding of the spike (1-line fix, wrong scientific number).
B2 is the biggest *runnability* risk — confirm the posterior path runs at all before other bayes work.

---

## Tier 1 — Performance (evidence-backed; behavior-preserving)

| # | Change | Location | Measured / expected | Effort |
|---|--------|----------|---------------------|--------|
| P1 | **Rewrite `correlation()`** — drop per-call `Dict` dynamic dispatch + `collect(view)` + Union buffers; single-pass filter into two reused `Float64` buffers via a function barrier `_corr(f::F,…)`. (= TYPE-6.) | `colocalization.jl:221-239` | **3.4–6.0× faster, 135–1394× less memory** (benchmarked) | M |
| P2 | **Hoist `InterpKDE`** out of the `quadgk` integrand in `compute_BayesFactor` and `bayes_rangeplot` (currently rebuilds a 2048-pt spline per node). | `bayes.jl:113-124`, `plot.jl:682-699` | Large constant-factor speedup of BF/range-plot | S |
| P3 | **`patch()` dense eltype** `Array{T,4}` instead of forced Union (complements P1). Update `test/runtests.jl:131`. | `colocalization.jl:46,86` | ~11% patch mem; unlocks P1's full win | S |
| P4 | **Dedupe correlation** computed for the patched-correlation plot vs inference (compute once, share). Do **after** B4. | `utils.jl:200-238`, `bayes.jl:164-173` | Removes one full patch+corr pass/channel-pair (bounded by ADVI) | M |
| P5 | **`_prepare_data` cleanup** — collect per-image vectors directly into `Vector{Vector{Float64}}`; drop the Union 3D array + aliasing `fill` idiom. | `bayes.jl:164-201` | Modest (runs 2×/run) | S |

*Honest framing (verifier corrections):* `patch.([x,y])` (PERF-3) and the abstract `get_images`
container (PERF-7) are **cleanups, not speedups**; threading (PERF-8) is bounded by `n_images` and
dwarfed by ADVI — and is only safe if P1's buffers stay thread-local.

---

## Tier 2 — Extensibility & structure hierarchy

| # | Change | Location | Why | Effort |
|---|--------|----------|-----|--------|
| E1 | **Lift the `@model` to module scope** (`@model coloc_model`); split `fit_prior`/`fit_posterior`; move the orphan docstring onto it. Verified no outer-scope capture → behavior-preserving. | `bayes.jl:234-339` | Unit-testable; **reusable prior ranges for v2.0 SBI**; breaks the monolith | M |
| E2 | **Real image interface** — erroring fallbacks on the abstract type + concrete impls on `MultiChannelImage`; keep `ArgumentError`. | `LoadImages.jl:31-85` | A new image backend becomes implementable | S |
| E3 | **Route all field access through accessors** (`image_data`, `image_name`, …); add `image_name!` setter. | `bayes.jl`,`plot.jl`,`utils.jl`,`main.jl`,`LoadImages.jl` | Makes E2 actually substitutable | M |
| E4 | **`MultiChannelImageStack <: AbstractVector{T}`** (define `size`+`getindex`, delete hand-rolled methods); add `AbstractImageStack`; relax `::MultiChannelImageStack` signatures. Drop `mutable` (TYPE-8). | `LoadImages.jl:148-196`, `bayes.jl:234-236` | Completes the container protocol; enables synthetic stacks (v2.0) | M |
| E5 | **Parametrize `CoLocResult{St,A}` + `AbstractCoLocResult`** — annotate `advi_result::A`, **without** renaming fields. Update `typeof(x)==CoLocResult` tests. | `bayes.jl:75-82` | Type-clean, dispatchable result hierarchy | M |
| E6 | **Concrete container eltype** in `get_images` (`{Float64,String,Float64,Int}`) and `main.jl` (`eltype(stack.img)`). | `utils.jl:63`, `main.jl:140` | Robustness vs future struct-param changes | S |
| E7 | *(deferred, L)* **Joint N-channel model.** Pairwise N-channel already works via fan-out; joint modeling needs the model + `.μ_sample`/`.μ_control` consumers reworked. | `bayes.jl:267-303`, consumers | True multi-channel inference | L |

---

## Tier 3 — Maintenance, DRY & data visualization (Goal 3a)

| # | Change | Location | Effort |
|---|--------|----------|--------|
| M1 | **`delta_rho(::CoLocResult)` helper** — replaces 7 inline Δρ copies; structurally prevents B1-class drift. | `bayes.jl` + 7 sites | S |
| M2 | **Centralize BF math** — `_tail_prob`/`_bayes_factor` (or `bf_at`) shared by `compute_BayesFactor` + `bayes_rangeplot`; keep the `ε>1e-5` warning + NaN guard. (= DRY-1 + VIZ-4; combine with P2.) | `bayes.jl:109-137`, `plot.jl:680-699` | M |
| M3 | **Shared plotting theme + helpers** — `cm_to_pt`, `coloc_figure`, scoped `with_theme(COLOC_THEME)`; delete dead `cm_to_px`/`calculate_font_size`. **Scoped, not global `set_theme!`.** | `plot.jl` | M |
| M4 | **Real error logging** — `catch e; @warn … exception=(e, catch_backtrace())` at the 5 plot sites. | `utils.jl` (5 sites) | S |
| M5 | **`plot_images` dispatch** via `_plotter(::Val{…})`; keep the shared per-image `try/catch`; fix the `$suffix_` interp bug. | `utils.jl:125-159` | S |
| M6 | **`min_pairs::Int=15` keyword** for the patch-inclusion threshold (default preserves behavior). | `colocalization.jl:235` | S |
| M7 | **Draw the 8 prior globals directly** instead of a full Chains then discard (INFER-3); fold `num_latent` cleanup (INFER-2) into B2. | `bayes.jl:308-323` | M |

---

## Tier 4 — Visualization backend & provenance (Goal 3b)

| # | Change | Location | Effort |
|---|--------|----------|--------|
| V1 | **CairoMakie for file output** — headless/CI/GUI-batch safe; restores the advertised `.svg`/`.pdf` export; change `endswith||@error` guards to `throw(ArgumentError)`. | `plot.jl`, `ProteinCoLoc.jl:28`, `Project.toml` | M |
| R1 | **Thread a seed** — `seed=nothing` kwarg on `start_analysis` → `rng=Xoshiro(seed)` into `sample`/`rand`/shuffles; `Random.seed!` before `vi`. (Depends on B2.) | `main.jl`, `bayes.jl`, `LoadImages.jl` | M |
| R2 | **Run-manifest** — `write_manifest()` emitting `run_manifest.toml`: Julia + pkg versions, seed, timestamp, git commit, all args, per-input SHA-256. | `main.jl` | M |
| R3 | **Embed provenance in `CoLocResult`** — `provenance::NamedTuple` (`cor_method`, `n_iter`, `n_posterior_samples`) via a defaulted constructor. (Couples with E5.) | `bayes.jl` | M |
| O1 | **Structured result writer** — one named/typed `DataFrame` row + `CSV.write(...; append=isfile)`; keeps ints as ints, full precision, no header/value drift. | `utils.jl:353-401` | M |
| O2 | **Fix output guard** — guard on `result.txt`/the manifest (not the GUI-only `log.txt`); open `result.txt` `"w"` once per run, append rows within the run. | `main.jl:119`, `utils.jl` | S |

---

## Suggested execution order (dependencies)

1. **B2** (port Turing API; resolve Manifest) — *gate: does the posterior even run?*
2. **B1, B4, B3, B5** — correctness bugs (B1 is 1 line; pair B1 with M1 `delta_rho`).
3. **P1 (+P3, +TYPE-6 dispatch), P2** — the measured performance wins.
4. **E1** (lift `@model`) — unlocks testing + v2.0 prior reuse; then **E2/E3, E4, E5/E6**.
5. **M1–M7** DRY/maintenance (M2 with P2; M1 with B1).
6. **V1, R1→R2→R3, O1, O2** — viz backend + provenance (R1 before R2 before R3).
7. *(deferred)* **E7** joint N-channel model.

## Cross-cutting notes

- **Same fix, multiple findings:** P1 ≡ TYPE-6 (correlation dispatch); B1 ≡ POST-3 ≡ VIZ-2 ≡ BUG-1;
  M2 ≡ DRY-1 ≡ VIZ-4. Implement once.
- **Tests:** P3 and E5 require updating assertions in `test/runtests.jl` (lines 131 and 203–298).
- **Not byte-identical:** V1 (CairoMakie rasterizer) and the bug fixes (B1, B3) intentionally change
  outputs — regenerate any cached result files / golden images.
