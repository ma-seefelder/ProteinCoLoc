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
# MOVED (v2.0 breaking release, D-01 / D-03).
#
# The Turing/ADVI reference colocalization path formerly defined here — the `@model`, ADVI
# `colocalization()`, `compute_BayesFactor()`, `convert_posterior_samples()`, `_prepare_data()`
# and the `CoLocResult` struct — now lives in `ext/ProteinCoLocTuringExt.jl` as the internal,
# weakdep-gated reference implementation (`CoLocResult` renamed to `AdviColocResult`, a subtype
# of `AbstractColocResult`). It loads ONLY when `Turing` is present and is NOT exported.
#
# This file is intentionally no longer `include`d by `src/ProteinCoLoc.jl`. The public v2.0
# surface is the amortized NPE/NRE path; see `src/results.jl` for the result-type hierarchy.
#############################################################################################
