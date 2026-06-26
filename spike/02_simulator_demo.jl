#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder
E-Mail: manuel.seefelder@uni-ulm.de
Postal address: Department of Gene Therapy, University of Ulm, Helmholzstr. 8/1, 89081 Ulm, Germany

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=#

# spike/02_simulator_demo.jl --- SIM-04 plausibility evidence (D-13): the seeded,
# headless demo that regenerates the static CairoMakie figures under spike/figures/.
#
# Mirrors the spike/00_smoke.jl script shape (AGPL header, fixed-seed top-level
# CPU-only script, correctness self-assert at the end). It exercises the FULL
# sample_prior -> simulate_pair -> build_mci -> summary chain (Waves 1-3) and renders:
#   (row 1) image heatmaps for a high-ρ and a low-ρ channel pair,
#   (row 2) the ρ_true -> mean-patch-correlation monotonicity sweep,
#           the induced-μ histogram overlaid on the Turing μ-prior density,
#           the spillover-effect curve, and the sub-pixel-shift-effect curve.
#
# DECOUPLING (CLAUDE.md / D-13): CairoMakie is spike-local; this script writes ONLY
# under spike/figures/ and includes src/ read-only via contract.jl. The baseline
# f581d95 (src/ + root manifests) stays byte-identical.
#
# PowerShell-safe verification (02-04-PLAN critical note 3): the SCRIPT itself
# self-asserts `@assert isfile(...)` on its PNG+PDF outputs so a missing figure
# fails the run -- the automated verify command is simply
#     julia --project=spike spike/02_simulator_demo.jl
# (no POSIX `test -f`).

using CairoMakie
using StatsBase
using Statistics
using Random
using Distributions

CairoMakie.activate!()                       # headless raster/vector, no OpenGL/display

# --- sample -> simulate -> summary chain (Waves 1-3), read-only src/ via contract --
include(joinpath(@__DIR__, "contract.jl"))               # build_mci, patch_summary, induced_mu
include(joinpath(@__DIR__, "simulator", "forward.jl"))   # simulate_pair
include(joinpath(@__DIR__, "simulator", "prior.jl"))     # sample_prior, MU_PRIOR, ghat, GHAT_MU_*

# --- Demo configuration (deterministic, CPU-only) -------------------------------
const DEMO_RNG    = Random.Xoshiro(2026)     # D-14: explicit, threaded rng
const DEMO_IMSIZE = (256, 256)               # modest size -> finishes in minutes on CPU
const DEMO_NREP   = 3                         # Monte-Carlo replicates averaged per sweep point

# θ with the 6 nuisances pinned at their prior MEDIANS, isolating ρ_true's effect
# (matches prior.jl's nuisance ranges: spillover U(0,.2), autofluor U(0,.1),
#  label_eff U(.6,1), shift U(-1,1), noise U(0,1)).
_θ_med(ρ) = (ρ_true = ρ, spillover = 0.1, autofluorescence = 0.05,
             label_efficiency = 0.8, shift_dx = 0.0, shift_dy = 0.0, noise = 0.5)

# Mean per-patch correlation through the REAL summary, averaged over DEMO_NREP sims.
function mean_patch_corr(rng::AbstractRNG, θ; imsize = DEMO_IMSIZE, n = DEMO_NREP)
    vals = Float64[]
    for _ in 1:n
        m = induced_mu(build_mci(simulate_pair(rng, θ; imsize = imsize)))
        isfinite(m) && push!(vals, m)
    end
    return mean(vals)
end

# --- (a) ρ_true monotonicity sweep (≥15 points) ---------------------------------
ρ_grid    = collect(range(-0.9, 0.9; length = 15))
mean_corr = [mean_patch_corr(DEMO_RNG, _θ_med(ρ)) for ρ in ρ_grid]

# --- (b) spillover-effect sweep at fixed ρ_true ---------------------------------
const ρ_FIX = 0.3
spill_grid = collect(range(0.0, 0.2; length = 8))
corr_spill = [mean_patch_corr(DEMO_RNG, merge(_θ_med(ρ_FIX), (spillover = s,))) for s in spill_grid]

# --- (b) sub-pixel |shift|-effect sweep at fixed ρ_true -------------------------
shift_steps = collect(range(0.0, 1.0; length = 8))           # per-axis shift
shift_mag   = sqrt.(2.0 .* shift_steps .^ 2)                  # Euclidean |shift|
corr_shift  = [mean_patch_corr(DEMO_RNG, merge(_θ_med(ρ_FIX),
                  (shift_dx = s, shift_dy = s))) for s in shift_steps]

# --- (c) induced-μ samples under the FULL prior (sample_prior) -------------------
μ_samples = Float64[]
for _ in 1:150
    θ = sample_prior(DEMO_RNG)
    m = induced_mu(build_mci(simulate_pair(DEMO_RNG, θ; imsize = DEMO_IMSIZE)))
    isfinite(m) && push!(μ_samples, m)
end

# --- Image panels: a high-ρ and a low-ρ channel pair (shared seed) ---------------
pair_hi = simulate_pair(Random.Xoshiro(7), _θ_med(0.7);  imsize = DEMO_IMSIZE)
pair_lo = simulate_pair(Random.Xoshiro(7), _θ_med(-0.7); imsize = DEMO_IMSIZE)

# --- Assemble the figure ---------------------------------------------------------
fig = Figure(size = (1400, 760))

heatmap(fig[1, 1], pair_hi[1]; colormap = :viridis,
        axis = (title = "high ρ_true=+0.7  ch1", aspect = DataAspect()))
heatmap(fig[1, 2], pair_hi[2]; colormap = :viridis,
        axis = (title = "high ρ_true=+0.7  ch2", aspect = DataAspect()))
heatmap(fig[1, 3], pair_lo[1]; colormap = :viridis,
        axis = (title = "low ρ_true=-0.7  ch1", aspect = DataAspect()))
heatmap(fig[1, 4], pair_lo[2]; colormap = :viridis,
        axis = (title = "low ρ_true=-0.7  ch2", aspect = DataAspect()))

ax_sweep = CairoMakie.Axis(fig[2, 1]; title = "(a) ρ_true → mean patch-corr",
                xlabel = "ρ_true", ylabel = "mean patch-corr")
scatterlines!(ax_sweep, ρ_grid, mean_corr)

ax_mu = CairoMakie.Axis(fig[2, 2]; title = "(c) induced μ vs Turing μ-prior",
             xlabel = "induced μ", ylabel = "density")
hist!(ax_mu, μ_samples; normalization = :pdf, bins = 20, color = (:steelblue, 0.6))
let target = Truncated(Cauchy(0.0, 0.3), GHAT_MU_MIN, GHAT_MU_MAX),
    xs = range(GHAT_MU_MIN, GHAT_MU_MAX; length = 200)
    lines!(ax_mu, xs, pdf.(target, xs); color = :firebrick, linewidth = 2)
end

ax_spill = CairoMakie.Axis(fig[2, 3]; title = "(b) spillover ↑ → corr ↑",
                xlabel = "spillover", ylabel = "mean patch-corr")
scatterlines!(ax_spill, spill_grid, corr_spill)

ax_shift = CairoMakie.Axis(fig[2, 4]; title = "(b) |shift| ↑ → corr ↓",
                xlabel = "|sub-pixel shift|", ylabel = "mean patch-corr")
scatterlines!(ax_shift, shift_mag, corr_shift)

# --- Save (raster + vector); create spike/figures/ if absent ---------------------
const FIG_DIR = joinpath(@__DIR__, "figures")
isdir(FIG_DIR) || mkpath(FIG_DIR)
const PNG_PATH = joinpath(FIG_DIR, "plausibility.png")
const PDF_PATH = joinpath(FIG_DIR, "plausibility.pdf")
save(PNG_PATH, fig)                          # raster; format inferred from extension
save(PDF_PATH, fig)                          # vector

# --- Correctness self-assert (00_smoke.jl idiom; PowerShell-safe, no POSIX test) -
@assert isfile(PNG_PATH) "demo FAILED: missing $PNG_PATH"
@assert isfile(PDF_PATH) "demo FAILED: missing $PDF_PATH"

println("demo OK: wrote ", PNG_PATH, " and ", PDF_PATH)
println("  monotonicity Spearman(ρ_true, mean-corr) = ", corspearman(ρ_grid, mean_corr))
println("  spillover effect  Δcorr = ", corr_spill[end] - corr_spill[1])
println("  |shift| effect    Δcorr = ", corr_shift[1] - corr_shift[end])
