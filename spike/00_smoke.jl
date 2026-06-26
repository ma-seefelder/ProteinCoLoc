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

# spike/00_smoke.jl --- Phase 1 environment + smoke gate (ENV-02).
#
# Highest-value de-risking step of the whole spike: prove the pre-1.0
# NeuralEstimators v0.2.1 + Flux v0.16 stack trains a normalising-flow
# posterior estimator and recovers a KNOWN parameter on this Windows
# machine, CPU-only (D-04). The gate is "green + correctness" (D-03):
# it does not merely run -- it asserts the recovered posterior mean is
# within tolerance of the true value, which catches a silently-wrong API
# or backend. Synthetic 1-parameter Gaussian only; needs nothing from src/.
#
# API written against the INSTALLED v0.2.1 docstrings (pre-1.0 drift, P4):
#   * PosteriorEstimator(summary_network, q)  -- q is a constructed
#     NormalisingFlow INSTANCE passed POSITIONALLY (the `q=` keyword form
#     of the convenience constructor expects a TYPE, not an instance).
#   * NormalisingFlow(d; num_summaries = dstar)
#   * sampleposterior(est, Z; N = ...)        -- N is a KEYWORD.
#   * train(...; use_gpu = false)             -- use_gpu DEFAULTS to true.

using NeuralEstimators, Flux, Distributions, Random

Random.seed!(2026)   # determinism: stdlib Random (Random123 is a later phase)

# --- Toy model: theta ~ Normal(0,1);  Z | theta ~ Normal(theta, sigma^2),
#     observed as m i.i.d. replicates per data set. -----------------------
const d             = 1        # one unknown parameter (the mean theta)
const m             = 200      # replicates per data set (informative data)
const sigma         = 1.0f0    # known noise scale
const num_summaries = 8        # summary-statistic dimension (>= d)

# prior sampler: returns a d x K matrix of theta draws (theta ~ N(0,1))
prior_sampler(K) = Float32.(randn(1, K))

# simulator: given a d x K parameter matrix, return an m x K matrix whose
# k-th column is one data set of m sorted replicates Z | theta_k.
function gaussian_simulator(theta::AbstractMatrix)
    K = size(theta, 2)
    Float32.(reduce(hcat, [theta[1, k] .+ sigma .* sort(randn(m)) for k in 1:K]))
end

# summary network: maps an m-vector data set to num_summaries summaries.
network = Chain(Dense(m, 64, gelu), Dense(64, 64, gelu), Dense(64, num_summaries))
q       = NormalisingFlow(d; num_summaries = num_summaries)
est     = PosteriorEstimator(network, q)

# CPU-only training (D-04: use_gpu = false is mandatory).
est = train(est, prior_sampler, gaussian_simulator;
            K = 8000, epochs = 60, use_gpu = false, verbose = false)

# --- Inference on data generated from a KNOWN theta -----------------------
theta_true = 0.7f0
Z_obs      = Float32.(reshape(theta_true .+ sigma .* sort(randn(m)), m, 1))
draws      = sampleposterior(est, Z_obs; N = 1000)   # d x N matrix
mu_hat     = posteriormean(draws)[1]                 # mean over the N draws

# Correctness gate (D-03): generous tolerance -- its job is to separate a
# correct CPU backend from a silently-broken one, not to certify precision.
tol = 0.3
@assert abs(mu_hat - theta_true) < tol "smoke FAILED: recovered $mu_hat vs true $theta_true (tol $tol)"

println("smoke OK: recovered mu_hat = ", mu_hat, " vs theta_true = ", theta_true,
        "  (|delta| = ", abs(mu_hat - theta_true), " < tol = ", tol, ")")
