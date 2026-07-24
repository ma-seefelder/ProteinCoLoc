---
spike: 004
name: visualization-reuse
type: analysis
validates: "Given plot.jl + utils.jl plotting (1.1k lines), when audited, then duplication/coupling is identified and a reusable, themed, headless-safe plotting layer is proposed"
verdict: VALIDATED
related: [001, 003, 005]
tags: [visualization, dry, maintenance]
---

# Spike 004: Visualization layer — reuse & maintainability

## What This Validates

Given the plotting layer (`plot.jl` 736 lines + the `plot_images`/`generate_plots` orchestration in
`utils.jl`), when audited for duplication, coupling and maintainability, then a reusable, themed,
headless-safe plotting layer is proposed.

## Research / Method

Read all of `plot.jl` (7 plot functions) and the `utils.jl` orchestration. Each finding verified
against the code (grep-confirmed counts for the duplication claims).

## Findings (all verified)

| ID | Finding | Verdict | Sev | Bug? | Effort |
|----|---------|---------|-----|------|--------|
| VIZ-1 | `Δρ = μ_sample .- μ_control` hand-coded in **7** places | ✓ confirmed | Med | no | S |
| VIZ-2 | `generate_txt` Δρ uses posterior μ_sample minus **prior** μ_control | ✓ confirmed | High | **yes** | S |
| VIZ-6 | Bare `try/catch` swallows every plot exception into a generic `@warn` | ✓ confirmed | Med | no | S |
| VIZ-4 | `bayes_rangeplot` recomputes KDEs + re-implements the BF odds | ✓ confirmed | Med | no | M |
| VIZ-7 | Hard GLMakie dep: no headless support, advertised `.svg` export is dead | ✓ confirmed | Med | no | M |
| VIZ-3 | Figure-sizing/theming boilerplate duplicated; 2 dead helpers; no shared Theme | ⚠ needs-nuance | Med | no | M |
| VIZ-5 | `plot_images` dispatches plot type on a `Symbol` via if/elseif | ⚠ needs-nuance | Low | no | S |

### VIZ-1 — extract a `delta_rho` helper — ✓ confirmed (med)
`bayes.jl:110-111`, `plot.jl:527,597-598,680-681`, `utils.jl:368`. Grep confirms **8** copies of the
`μ_sample .- μ_control` pattern; 7 are the correct within-distribution form, the 8th
(`utils.jl:368`) is the buggy one (VIZ-2). A single `delta_rho(r::CoLocResult) = r.posterior.μ_sample
.- r.posterior.μ_control` removes the duplication **and structurally prevents the drift that caused
the bug**. (Finder's "shadowing" caveat is spurious — `delta_rho` and `Δρ` are distinct identifiers.)

### VIZ-2 — Δρ wrong distribution — ✓ confirmed (high, BUG)
`utils.jl:368,358-359`. Same bug as POST-3/CORR-1/BUG-1 (found by 4 dimensions). Subtracting
independent prior & posterior draws inflates variance and biases the mean toward the prior control
location (`Cauchy(0,0.3)`, centered ~0). Fix: `posterior_samples.:μ_control` (or `delta_rho`).

### VIZ-6 — error-swallowing `try/catch` — ✓ confirmed (med)
`utils.jl:138-145,147-154,248-252,257-261,266-271`. Five plot calls use a bare `catch` (no exception
var) emitting only "The X plot could not be generated." The object + backtrace are discarded, so a
real defect (renamed column, save failure, the `.svg` issue, OOM) is indistinguishable from an
expected skip. The expensive inference (`compute_BayesFactor`, `colocalization`) runs **outside** the
try/catch, so these blocks only mask rendering errors. Fix:
`catch e; @warn "…failed" exception=(e, catch_backtrace())` at all five sites.

### VIZ-4 — `bayes_rangeplot` re-does inference work — ✓ confirmed (med)
`plot.jl:682-699`, `bayes.jl:113-136`, `utils.jl:245,258`. The rangeplot refits `kde(Δρ_post/prior)`
and re-derives the odds ratio that `compute_BayesFactor` already computed; in the standard run both
run back-to-back per channel pair (~100k-sample KDEs fit twice). Implementations have already
diverged (rangeplot drops the `ε` warning, adds a NaN guard). Fix: a shared `bf_at(prior_kde,
post_kde, ρ0)`; keep `compute_BayesFactor`'s `(bf,p_post,p_prior)` return + `ε>1e-5` warning and the
rangeplot's `a>0 ? … : NaN` guard. Combine with PERF-5 (hoist `InterpKDE`). Overlaps DRY-1.

### VIZ-7 — GLMakie / headless / dead `.svg` — ✓ confirmed (med)
`plot.jl:123-212,134,404,175`, `ProteinCoLoc.jl:28`, `Project.toml:9`. Every figure is built/saved
with GLMakie (OpenGL/GLFW — needs a live graphics context), unsuitable for headless servers/CI/GUI
batch runs that only need files. `plot()` validates `.svg` (`plot.jl:134`) but GLMakie can only
rasterize → `GLMakie.save("x.svg",fig)` throws: the advertised `.svg` path is **dead on arrival**.
Fix: render file output with **CairoMakie** (PNG+SVG+PDF, no GL context; already used in `spike/`),
and change the `endswith(...) || @error` guards to `throw(ArgumentError(...))` (`@error` only logs,
doesn't halt). **Caveats:** not byte-identical (different rasterizer); the single `activate!()`
coupling is overstated (GLMakie is the sole backend, auto-activated).

### VIZ-3 — theming boilerplate + dead helpers — ⚠ needs-nuance (med)
`plot.jl:463-465,600-605,714-720` (cm→pt + identical `Figure` kwargs ×3), `colormap=(:viridis,0.3)`
×10, `xticklabelrotation=deg2rad(60)` ×7, two near-identical `Legend` blocks (`519-524`,`632-637`);
`cm_to_px` (`57-61`) and `calculate_font_size` (`81-83`) are **dead** (zero call sites). Fix: a
`cm_to_pt` helper, a `coloc_figure(w,h)` helper, a shared `COLOC_THEME`, delete the dead functions.
**Critical caveat:** use **scoped `with_theme(COLOC_THEME) do … end`** around only the three Bayesian
plots — a global `set_theme!` would alter the image plots (`plot`/`plot_mask`/`_local_correlation_plot`
use a black background, `fontsize=12`, unrotated axes). Also verify ax5's deliberate grey density
(`plot.jl:538`) still wins over a themed `Density` colormap.

### VIZ-5 — `plot_images` Symbol dispatch — ⚠ needs-nuance (low)
`utils.jl:125-159`. The `if plot_type==:local_correlation … elseif … else error` branches differ
only in which function is called; bodies (loop/filename/try-catch) are duplicated. A `Val`-based
`_plotter(::Val{:patched_correlation})=plot` is reasonable. **Caveats:** the finder's sketch (a)
drops the per-image try/catch (would abort the whole loop on one failure) and (b) has a
string-interp bug (`$suffix_`). Keep one shared try/catch inside the loop; `Val(symbol)` on a
runtime Symbol is dynamic (fine — plotting isn't perf-critical). Symbols are internal-only.

## Correctness bug noticed here (detail in Spike 005)
- **VIZ-2 = the Δρ bug** (also POST-3/CORR-1/BUG-1).

## Signal for the Build
VIZ-1 (`delta_rho`) + VIZ-2 (the bug fix) go together and kill the bug class. VIZ-6 (real error
logging) is a one-liner ×5 with big debuggability payoff. VIZ-3/-4/-5 are DRY wins (scoped theme!).
VIZ-7 (CairoMakie) is the strategic move for headless/CI/manuscript `.svg`+`.pdf` output.
