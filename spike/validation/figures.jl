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

# spike/validation/figures.jl --- Phase-5 shared CairoMakie figure surface.
#
# The ONE figure module for the whole validation bundle, OWNED by this plan (05-01)
# so the Wave-2 plans (05-02 BF, 05-03 OOD) only CALL these functions and never edit
# figures.jl. Follows the Phase-2 `02_simulator_demo.jl` headless-save idiom
# (CairoMakie.activate!() → build Figure → save PNG). Every figure is written under
# spike/validation/figures/ (created on demand); PNGs are gitignored regenerable
# artifacts, not committed deliverables.
#
# HARD RULE: figure generation NEVER runs inside the fast test gate (Phase-4 pattern:
# figures are script artifacts, not gate assertions). test_sbc.jl calls NONE of these.
#
# DECOUPLING (CLAUDE.md): spike-local; touches no src/. Guarded include of sbc.jl
# brings the CalibrationResult type into scope for the reliability figure.

using CairoMakie
using Statistics

CairoMakie.activate!()   # headless raster/vector backend (no OpenGL/display)

isdefined(@__MODULE__, :CalibrationResult) || include(joinpath(@__DIR__, "sbc.jl"))

# --- shared output directory (created on demand) --------------------------------
const VAL_FIG_DIR = joinpath(@__DIR__, "figures")
_val_fig_path(name) = (isdir(VAL_FIG_DIR) || mkpath(VAL_FIG_DIR); joinpath(VAL_FIG_DIR, name))

"""
    plot_rank_histogram(ranks; bins = SBC_BINS, L = SBC_L, label = "θ",
                        filename = "sbc_rank_hist.png") -> String

SBC rank histogram for one parameter with the expected-uniform band. A calibrated
net gives a flat histogram; ∪/∩/slope shapes flag under-/over-/biased dispersion.
Saves a PNG under spike/validation/figures/ and returns its path. Captioned with
`SBC_CAPTION` (SBC-04).
"""
function plot_rank_histogram(ranks::AbstractVector{<:Integer};
                             bins::Integer = SBC_BINS, L::Integer = SBC_L,
                             label::AbstractString = "θ",
                             filename::AbstractString = "sbc_rank_hist.png")
    fig = Figure(size = (640, 460))
    ax  = CairoMakie.Axis(fig[1, 1]; title = "SBC rank histogram — $label",
                          subtitle = SBC_CAPTION,
                          xlabel = "rank (0:$L)", ylabel = "count")
    hist!(ax, ranks; bins = bins, color = (:steelblue, 0.7))
    expected = length(ranks) / bins
    hlines!(ax, [expected]; color = :firebrick, linestyle = :dash)
    path = _val_fig_path(filename)
    save(path, fig)
    return path
end

"""
    plot_coverage_curve(nominal, empirical; filename = "sbc_coverage.png") -> String

Nominal-vs-empirical coverage curve with the y=x calibration diagonal. Saves a PNG
and returns its path.
"""
function plot_coverage_curve(nominal::AbstractVector, empirical::AbstractVector;
                             filename::AbstractString = "sbc_coverage.png")
    fig = Figure(size = (520, 500))
    ax  = CairoMakie.Axis(fig[1, 1]; title = "SBC coverage", subtitle = SBC_CAPTION,
                          xlabel = "nominal coverage", ylabel = "empirical coverage")
    lines!(ax, [0, 1], [0, 1]; color = :gray, linestyle = :dash)   # calibration diagonal
    scatterlines!(ax, nominal, empirical; color = :steelblue)
    path = _val_fig_path(filename)
    save(path, fig)
    return path
end

"""
    plot_reliability_ece(cal::CalibrationResult; filename = "sbc_reliability.png") -> String

Reliability diagram (predicted vs observed rate) annotated with ECE/MCE and the
traffic-light verdict color (green/yellow/red on the pre-registered SBC_ECE_* cutoffs).
Saves a PNG and returns its path.
"""
function plot_reliability_ece(cal::CalibrationResult;
                              filename::AbstractString = "sbc_reliability.png")
    verdict = sbc_traffic_light(cal.ece)
    vcolor  = verdict === :green ? :seagreen : (verdict === :yellow ? :orange : :firebrick)
    fig = Figure(size = (520, 500))
    ax  = CairoMakie.Axis(fig[1, 1];
        title = "SBC reliability — ECE=$(round(cal.ece; digits=3)) " *
                "MCE=$(round(cal.mce; digits=3)) [$verdict]",
        subtitle = SBC_CAPTION,
        xlabel = "predicted rate", ylabel = "observed rate")
    lines!(ax, [0, 1], [0, 1]; color = :gray, linestyle = :dash)
    scatterlines!(ax, cal.predicted_rate, cal.observed_rate; color = vcolor)
    path = _val_fig_path(filename)
    save(path, fig)
    return path
end

"""
    plot_roc(fpr, tpr, auc; filename = "ood_roc.png") -> String

ROC curve for the OOD detector (Plan 05-03 CALLS this; owned here). Annotated with the
AUC. Saves a PNG and returns its path.
"""
function plot_roc(fpr::AbstractVector, tpr::AbstractVector, auc::Real;
                  filename::AbstractString = "ood_roc.png")
    fig = Figure(size = (520, 500))
    ax  = CairoMakie.Axis(fig[1, 1]; title = "OOD ROC — AUC=$(round(auc; digits=3))",
                          xlabel = "false positive rate", ylabel = "true positive rate")
    lines!(ax, [0, 1], [0, 1]; color = :gray, linestyle = :dash)   # chance line
    lines!(ax, fpr, tpr; color = :steelblue)
    path = _val_fig_path(filename)
    save(path, fig)
    return path
end

"""
    plot_bf_agreement(logbf_amortized, logbf_kde; filename = "bf_agreement.png") -> String

Amortized-vs-KDE-baseline log Bayes-factor agreement scatter with the y=x line (Plan
05-02 CALLS this; owned here). Saves a PNG and returns its path.
"""
function plot_bf_agreement(logbf_amortized::AbstractVector, logbf_kde::AbstractVector;
                           filename::AbstractString = "bf_agreement.png")
    fig = Figure(size = (520, 500))
    ax  = CairoMakie.Axis(fig[1, 1]; title = "amortized vs KDE-baseline log BF",
                          xlabel = "KDE-baseline log BF", ylabel = "amortized log BF")
    lo = min(minimum(logbf_amortized), minimum(logbf_kde))
    hi = max(maximum(logbf_amortized), maximum(logbf_kde))
    lines!(ax, [lo, hi], [lo, hi]; color = :gray, linestyle = :dash)
    scatter!(ax, logbf_kde, logbf_amortized; color = :steelblue)
    path = _val_fig_path(filename)
    save(path, fig)
    return path
end
