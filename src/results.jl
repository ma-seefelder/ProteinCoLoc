#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder

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

#############################################################################################
# D-02 result-type hierarchy (v2.0 "AmortizedColoc" breaking release)
#
# `AbstractColocResult` is the single governing supertype for every colocalization result.
# Downstream/plotting code dispatches on the SHARED ACCESSOR INTERFACE below
# (`delta_rho`, `bayes_factor`, `is_ood`, `posterior_draws`) — NOT on concrete fields — so
# new result variants slot in as new subtypes without touching consumers.
#
# Shipped subtype:   `AmortizedColocResult` (the public amortized NPE/NRE path).
# Internal subtype:  `AdviColocResult`      (the Turing/ADVI reference path; defined in
#                    ext/ProteinCoLocTuringExt.jl, weakdep-gated, NOT exported — D-01/D-03).
#
# Future extension points (Phase 12 spatial map, Phase 13 three-hypothesis BF) are sketched
# as comments at the bottom of this file — deliberately NOT implemented here.
#############################################################################################

"""
    abstract type AbstractColocResult

Supertype for every colocalization-analysis result. Concrete subtypes must implement the
shared accessor interface (`delta_rho`, `bayes_factor`, `is_ood`, `posterior_draws`) so that
plotting, reporting and decision code can stay type-agnostic and dispatch on the interface
rather than on backend-specific fields.
"""
abstract type AbstractColocResult end

# --- Shared accessor interface (defined on the supertype; concrete subtypes override) ------

_iface_error(r::AbstractColocResult, f::Symbol) = error(
    "`$(f)` is not implemented for $(typeof(r)). Every `AbstractColocResult` subtype must " *
    "implement the accessor interface: delta_rho, bayes_factor, is_ood, posterior_draws.")

"""
    delta_rho(r::AbstractColocResult)

Return the Δρ (sample − control colocalization contrast) carried by the result. For draw-based
results this is the vector of posterior Δρ draws.
"""
delta_rho(r::AbstractColocResult) = _iface_error(r, :delta_rho)

"""
    bayes_factor(r::AbstractColocResult)

Return the colocalization Bayes factor evidence carried by the result (amortized results
return the log Bayes factor).
"""
bayes_factor(r::AbstractColocResult) = _iface_error(r, :bayes_factor)

"""
    is_ood(r::AbstractColocResult)

Return whether the input was flagged out-of-distribution / misspecified by the OOD detector.
"""
is_ood(r::AbstractColocResult) = _iface_error(r, :is_ood)

"""
    posterior_draws(r::AbstractColocResult)

Return the posterior draws backing the result (parameter × draw layout).
"""
posterior_draws(r::AbstractColocResult) = _iface_error(r, :posterior_draws)

# --- Supporting value types ---------------------------------------------------------------

"""
    OODVerdict(score, flag, per_channel)

Out-of-distribution / misspecification verdict attached to an `AmortizedColocResult`.

# Fields
- `score::Float64`: fused OOD detector statistic (robust-z / distance).
- `flag::Bool`: `true` if the input crossed the pre-registered ID operating point.
- `per_channel::NamedTuple`: per-detector-channel breakdown (e.g. density, noise).
"""
struct OODVerdict
    score::Float64
    flag::Bool
    per_channel::NamedTuple
end

"""
    CalibrationMeta(bin_midpoints, predicted_rate, observed_rate, bin_counts, ece, mce, grid, gate)

The grid's PASSED D-05 ship-gate calibration report. The first six fields are the binned
reliability-diagram result ported (as a pattern, not a dependency) from the spike's
`CalibrationResult` (`spike/validation/sbc.jl`); `grid` and `gate` add provenance so a result
can prove which pre-registered gate it was calibrated under.

# Fields
- `bin_midpoints::Vector{Float64}`: reliability-bin midpoints.
- `predicted_rate::Vector{Float64}`: mean predicted probability per bin.
- `observed_rate::Vector{Float64}`: observed positive rate per bin.
- `bin_counts::Vector{Int}`: samples per bin.
- `ece::Float64`: Expected Calibration Error.
- `mce::Float64`: Maximum Calibration Error.
- `grid::Int`: patch grid `G` this calibration belongs to.
- `gate::NamedTuple`: gate provenance (e.g. `(; seed, consts_hash, passed)`).
"""
struct CalibrationMeta
    bin_midpoints::Vector{Float64}
    predicted_rate::Vector{Float64}
    observed_rate::Vector{Float64}
    bin_counts::Vector{Int}
    ece::Float64
    mce::Float64
    grid::Int
    gate::NamedTuple
end

# --- Shipped concrete subtype: the amortized (NPE/NRE) result ------------------------------

"""
    AmortizedColocResult

Result of the public amortized colocalization path: a single forward pass through the
grid-specific NPE (posterior) and NRE (Bayes factor), with an OOD flag and the grid's passed
calibration report. This is the v2.0 primary public result type.

# Fields
- `grid::Int`: patch grid `G` (the registry key) the estimator was trained for.
- `posterior::Matrix{Float64}`: 7×N physical-θ posterior draws (row 1 = `ρ_true`).
- `delta_rho_draws::Vector{Float64}`: Δρ (sample − control) posterior draws.
- `log_bayes_factor::Float64`: amortized log Bayes factor (coloc vs null) from the NRE.
- `ood::OODVerdict`: out-of-distribution / misspecification verdict.
- `calibration::CalibrationMeta`: the grid's PASSED D-05 ship-gate calibration report.
- `meta::NamedTuple`: run metadata (e.g. `(; N)`, channels, seed).
"""
struct AmortizedColocResult <: AbstractColocResult
    grid            :: Int
    posterior       :: Matrix{Float64}
    delta_rho_draws :: Vector{Float64}
    log_bayes_factor:: Float64
    ood             :: OODVerdict
    calibration     :: CalibrationMeta
    meta            :: NamedTuple
end

delta_rho(r::AmortizedColocResult)       = r.delta_rho_draws
bayes_factor(r::AmortizedColocResult)    = r.log_bayes_factor
is_ood(r::AmortizedColocResult)          = r.ood.flag
posterior_draws(r::AmortizedColocResult) = r.posterior

#############################################################################################
# Extension points for future roadmap result variants — SKETCH ONLY, do NOT implement here.
#
# Each future variant slots in as a NEW `AbstractColocResult` subtype plus (optionally) NEW
# accessor functions layered on top of the shared interface — never as bolted-on fields on the
# existing structs (D-02).
#
# Phase 13 — Three-Hypothesis Amortized Bayes Factor:
#   struct ThreeHypothesisColocResult <: AbstractColocResult
#       grid            :: Int
#       posterior       :: Matrix{Float64}
#       log_bf_simplex  :: NTuple{3, Float64}   # {null, coloc, anti-coloc}
#       ood             :: OODVerdict
#       calibration     :: CalibrationMeta
#       meta            :: NamedTuple
#   end
#   # new interface extension (dispatches on the supertype for the generic case):
#   bayes_factor_simplex(r::ThreeHypothesisColocResult) = r.log_bf_simplex
#
# Phase 12 — Spatial per-region Δρ map (GP/CAR lattice):
#   struct SpatialColocResult <: AbstractColocResult
#       grid            :: Int
#       region_delta_rho:: Matrix{Float64}      # per-region Δρ lattice
#       region_sd       :: Matrix{Float64}      # per-region posterior sd
#       ood             :: OODVerdict
#       calibration     :: CalibrationMeta
#       meta            :: NamedTuple
#   end
#   # new interface extensions:
#   delta_rho_map(r::SpatialColocResult)   = r.region_delta_rho
#   uncertainty_map(r::SpatialColocResult) = r.region_sd
#############################################################################################
