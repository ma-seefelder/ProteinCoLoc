---
spike: 001
name: hotpath-perf-memory
type: analysis
validates: "Given the patch→correlation→_prepare_data pipeline + KDE/Bayes-factor path, when profiled & read, then the allocation/type-stability hotspots are identified and an evidence-backed redesign is proposed"
verdict: VALIDATED
related: [002, 003, 004]
tags: [performance, memory, type-stability]
---

# Spike 001: Hot-path performance & memory

## What This Validates

Given the package's hottest code path (`patch` → `correlation` → `_prepare_data`, plus the
KDE/`quadgk` Bayes-factor path), when micro-benchmarked and read for type stability and
allocations, then the dominant hotspots are identified and a behavior-preserving redesign is
proposed with measured before/after numbers.

## Research / Method

- Read `colocalization.jl` (patch/correlation/_exclude_zero), `bayes.jl` (_prepare_data,
  compute_BayesFactor), `LoadImages.jl`, `utils.jl`, `main.jl`.
- Ran a standalone micro-benchmark (`bench_hotpath.jl`, results in `bench_results.txt`) using only
  stdlib `Statistics.cor` + built-in `@allocated`/`@elapsed` — no real images, no heavy deps.
- Every finding was adversarially verified against the actual code (see verdicts).

## How to Run

```bash
cd .planning/spikes/001-hotpath-perf-memory
julia --startup-file=no bench_hotpath.jl
```

## Results — measured (Julia 1.12.6, synthetic 50%-zero images)

The headline is **`correlation()`**: dropping the per-call `Dict` dynamic dispatch + the
`vec(collect(view(...)))` Union materializations + Union working buffers, fused into a single-pass
filter over two reused `Float64` buffers:

| Image | patches | `correlation` current | rewritten | speedup / mem |
|-------|---------|-----------------------|-----------|---------------|
| 1024² | 8×8 | 24.7 ms / 35.7 MB | 4.1 ms / 0.26 MB | **6.0× / 136×** |
| 2048² | 8×8 | 96.2 ms / 142.6 MB | 28.3 ms / 1.05 MB | **3.4× / 136×** |
| 1024² | 32×32 | 23.6 ms / 36.2 MB | 10.1 ms / 0.03 MB | **2.35× / 1394×** |

`patch()` Union→dense is only ~1.1× on its own (mem 1.12×); its real cost is realized downstream.

## Findings (all verified)

| ID | Finding | Verdict | Sev | Effort |
|----|---------|---------|-----|--------|
| PERF-1 | `correlation()`: per-call `Dict` dynamic dispatch + `collect(view)` + Union buffers | ✓ confirmed | High | M |
| PERF-5 | `pdf(::UnivariateKDE,x)` rebuilds the InterpKDE spline on **every** `quadgk` node | ✓ confirmed | High | S |
| PERF-2 | `patch()` always allocates `Array{Union{Float64,Missing},4}` even for dense input | ⚠ needs-nuance | Med | S |
| PERF-4 | `_prepare_data` builds an intermediate Union 3D array + per-image `skipmissing/collect/filter` | ✓ confirmed | Med | S |
| PERF-6 | Patched-correlation plot and inference recompute the **same** correlation | ✓ confirmed | Med | M |
| PERF-3 | `patch.([x,y], n)` broadcast allocates throwaway containers | ⚠ needs-nuance | Low | S |
| PERF-7 | `get_images` builds a Vector with an unbound 4th type param (abstract eltype) | ⚠ overstated | Low | S |
| PERF-8 | Per-image correlation is embarrassingly parallel but serial | ⚠ needs-nuance | Low | M |

### PERF-1 — `correlation()` rewrite (the win) — ✓ confirmed (high)
`colocalization.jl:221-239`. `cor_dict = Dict(:pearson=>cor,…)` is rebuilt every call and
`cor_func` is inferred `::Function` (abstract) → the inner `cor_func(a,b)` is a **dynamic
dispatch**; each patch pair does `vec(collect(view(x,i,j,:,:)))` (two Union-matrix allocations) +
`_exclude_zero` (two more). Fix: resolve the function once and call through a function barrier
`_corr(f::F, x, y) where F` that filters zero/NaN/missing in one pass into two reused `Float64`
buffers and slices with `view(buf,1:k)`. Numerically identical (Pearson is order-invariant; same
`≤15` cutoff). `_exclude_zero` becomes dead code on this path.
*Verifier note:* the "1394×" upper bound comes from high-`np` rows; the two headline rows are ~135×.

### PERF-5 — InterpKDE rebuilt per integrand evaluation — ✓ confirmed (high)
`bayes.jl:113-124`, `plot.jl:682-699`. KernelDensity defines `pdf(k::UnivariateKDE,x) =
pdf(InterpKDE(k),x)`, so `x->pdf(kde_dist,x)` reconstructs a 2048-point quadratic B-spline on
**every** `quadgk` node. `bayes_rangeplot`'s comment says "cached KDEs" but only the `UnivariateKDE`
is cached — the spline is still rebuilt thousands of times across its threshold loop. Fix: hoist
`ik = InterpKDE(kde(v))` once and integrate `pdf(ik,x)`. Numerically identical.
*Verifier note:* default range is 21 thresholds (not 160), still hundreds–thousands of rebuilds.

### PERF-2 — `patch()` Union eltype — ⚠ needs-nuance (med)
`colocalization.jl:46,86`. Both methods hard-code the Union 4D array despite the `where T` param;
production data is dense `Float64`, so widening is needless and forces Union reads in `correlation`.
Fix: `Array{T,4}(undef,…)`. **Caveats:** memory drop is ~11% (1.12×), **not** "half"; and
`test/runtests.jl:131` asserts the Union type and must be updated. Best realized together with PERF-1.

### PERF-4 — `_prepare_data` intermediate Union array — ✓ confirmed (med)
`bayes.jl:164-201`. Builds `sample_image::Array{Union,3}`, copies each correlation matrix into a
slice, reshapes, then `filter(!isnan, collect(skipmissing(row)))` per image; also uses the
aliasing-prone `fill(Vector{Float64}(), n)` idiom. Fix: collect each image's valid correlations
directly into `Vector{Vector{Float64}}(undef,n)`. **Caveats:** `reshape` is free (negligible);
`_prepare_data` runs only twice per run, so savings are modest, not a hot-path win.

### PERF-6 — duplicate correlation (plot vs inference) — ✓ confirmed (med)
`utils.jl:200-238`, `plot.jl:193-195`, `bayes.jl:164-173`. With the (default) patched-correlation
plot on, `plot()` computes `correlation(patch(x),patch(y))` per image, then `colocalization →
_prepare_data` computes it **again** at the same `num_patches`. Fix: compute once per
(stack, channel-pair, num_patches) and share. **Caveats:** total wall-clock gain is bounded by the
dominant ADVI cost; and a correct shared path must reconcile the `minmax_norm!` mutation (see CORR-2).

### PERF-3 / PERF-7 / PERF-8 — smaller / corrected
- **PERF-3** (`patch.([x,y],n)`): real but the allocations are tiny pointer-vectors (no data copy)
  and the destructure is already type-stable → reframe as a **readability** cleanup, not perf.
- **PERF-7** (abstract `get_images` container): only the unbound `I`/`pixel_size` is affected, and
  `pixel_dimensions` is never called; `MultiChannelImage` is a `mutable struct` (boxes regardless).
  Real value is robustness/inference-clarity, **not** measurable speed. Use `Int`, not `Int64`
  (32-bit portability). Couples with TYPE-5.
- **PERF-8** (threading): the loop is **already** thread-safe; the genuine risk is the *opposite* —
  if PERF-1 introduces a shared scratch buffer, parallelizing creates a data race. Keep buffers
  thread-local. ROI is modest (bounded by `n_images`, dwarfed by ADVI). Minimal form:
  `Threads.@threads` on the existing `bayes.jl:168-173` loop (writes are to disjoint slices).

## Correctness bugs noticed here (detail in Spike 005)
- **CORR-2** (`plot.jl:149`): `minmax_norm!` mutates `image.data` in place; because plotting runs
  before inference, enabling plots silently changes inference inputs. Fix: `minmax_norm!(copy(...))`.

## Signal for the Build
Do **PERF-1 + PERF-5 first** (largest, lowest-risk, behavior-preserving wins). PERF-2 rides along
with PERF-1 (and TYPE-6 is the same dispatch fix). PERF-6 only after CORR-2. PERF-7/PERF-8 are
cleanups, not speedups — frame them honestly.
