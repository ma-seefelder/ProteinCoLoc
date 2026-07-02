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

# spike/comparator/audit.jl --- CMP-07: audit summary (D-11)
#
# The BayesInteractomics-style audit summary over a built comparison table (D-11):
# a Markdown report counting the rows in each green/amber/red/gray traffic-light band
# (how many rows the classics get "wrong" vs the simulator ground truth) plus a
# KS-uniformity statistic on the finite Costes-p column (the well-specified null
# expects Uniform(0,1) p-values). Renders via an `IOBuffer` Markdown table in the
# `predictive_checks.jl:generate_diagnostics_report` style.
#
# The `_ks_uniform` helper MIRRORS BayesInteractomics `_ks_test_uniform`
# (predictive_checks.jl:729-736): the max absolute deviation between the empirical
# CDF and the Uniform(0,1) CDF over the finite p-values — reimplemented here rather
# than imported (decoupling; the sibling package is a read-only template).
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local, read-only over the DataFrame;
# no src/ reach, no Phase-3/Phase-5 files touched.

using DataFrames   # nrow / column access over the comparison table (D-09)
using Dates        # now() timestamp in the report header (stdlib)

"""
    _ks_uniform(pvals) -> Float64

The Kolmogorov-Smirnov statistic for `Uniform(0,1)`: the maximum absolute deviation
between the empirical CDF of the finite entries of `pvals` and the Uniform(0,1) CDF.
Small values indicate p-values consistent with uniformity (a well-specified Costes
null). Missing/NaN/non-finite entries are filtered; an empty finite set returns `NaN`.
Mirrors BayesInteractomics `_ks_test_uniform` (predictive_checks.jl:729-736).
"""
function _ks_uniform(pvals)
    valid = Float64[float(x) for x in pvals if !ismissing(x) && isfinite(x)]
    isempty(valid) && return NaN
    sorted    = sort(valid)
    n         = length(sorted)
    ecdf_vals = (1:n) ./ n
    return maximum(abs.(ecdf_vals .- sorted))
end

"""
    comparator_audit(df) -> String

The BayesInteractomics-style audit summary (D-11) for a comparison table `df` built
by `build_table`. Returns a Markdown report string with

  * a **traffic-light band-count** table (green/amber/red/gray rows + fractions) —
    the red/amber rows are the concrete "classic is wrong here" cases (D-10); and
  * a **KS-uniformity** statistic on the finite `Costes_p` column — the well-specified
    Costes null expects Uniform(0,1) p-values, so a large KS deviation flags
    misspecification.

Rendered via an `IOBuffer` in the `generate_diagnostics_report` Markdown-table style.
"""
function comparator_audit(df)
    io = IOBuffer()

    # Header
    println(io, "# Cross-Method Comparator Audit")
    println(io, "Generated: $(Dates.now()) | Package: ProteinCoLoc v2.0 (spike/comparator)")
    println(io)

    total = nrow(df)

    # --- Traffic-light band counts (D-10) ---
    counts = Dict{Symbol,Int}(:green => 0, :amber => 0, :red => 0, :gray => 0)
    for l in df.light
        counts[l] = get(counts, l, 0) + 1
    end

    println(io, "## Traffic-Light Band Counts")
    println(io)
    println(io, "Red/amber rows are the cases where the classical verdict diverges from the")
    println(io, "simulator ground-truth regime (the \"classic is wrong here\" rows, D-10).")
    println(io)
    println(io, "| Band | Rows | Fraction |")
    println(io, "|------|------|----------|")
    for b in (:green, :amber, :red, :gray)
        c    = get(counts, b, 0)
        frac = total > 0 ? round(c / total; digits = 3) : 0.0
        println(io, "| $b | $c | $frac |")
    end
    println(io, "| **total** | $total | 1.0 |")
    println(io)

    # --- Costes-p KS-uniformity ---
    ks      = _ks_uniform(df.Costes_p)
    n_fin   = count(x -> !ismissing(x) && isfinite(x), df.Costes_p)
    ks_str  = isnan(ks) ? "NaN (no finite Costes-p)" : string(round(ks; digits = 4))

    println(io, "## Costes-p KS-Uniformity")
    println(io)
    println(io, "A well-specified Costes randomization null yields Uniform(0,1) p-values;")
    println(io, "a large KS deviation from uniformity flags misspecification.")
    println(io)
    println(io, "| Metric | Value |")
    println(io, "|--------|-------|")
    println(io, "| KS deviation from Uniform(0,1) | $ks_str |")
    println(io, "| Finite Costes-p count | $n_fin / $total |")
    println(io)

    return String(take!(io))
end
