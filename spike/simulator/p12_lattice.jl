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

# spike/simulator/p12_lattice.jl --- the spatial prior's ARITHMETIC ONLY: the two covariance
# arms (CAR and GP), the per-cell rescale that makes D-05 true, the induced-lag-1
# reparametrization that makes D-08 a fair question, and the exact orthonormal DCT-II basis
# that lets D-01 and D-07 hold at once.
#
# THIS FILE DRAWS NO theta AND KNOWS NOTHING ABOUT rho. It produces 64x64 covariance matrices,
# a closure that draws marginally-standard 64-vectors from one, and an orthonormal rotation of
# the 8x8 lattice. Turning a drawn z-field into a rho-field -- the Normal CDF, the MU_PRIOR
# quantile, ghat -- is `p12_prior.jl`'s job (D-05), not this file's. Keeping the two apart is
# what lets the copula be tested against a field whose marginals are already known to be right.
#
# THE COLUMN-MAJOR INDEXING CONTRACT IS INHERITED, NOT CHOSEN, AND MUST NEVER BE TRANSPOSED.
# `src/amortized/summary.jl:72-76` (`encode_d01`) builds `vec(coalesce.(M, 0.0))` COLUMN-MAJOR
# from the GxG matrix `patch_summary` returns, so summary row `(j-1)*G + i` IS patch `(i, j)`
# and `reshape(rows[1:G^2], G, G)` recovers that matrix with NO transpose. Every adjacency,
# radius and basis index below therefore goes through `p12_idx`. A transposed variant would
# rotate the lattice relative to the image and would pass every symmetric test in the suite --
# which is precisely why the contract is asserted (test_p12_lattice.jl testset 8) rather than
# merely written down here.
#
# THE sqrt(diag) RESCALE IS NOT AN OPTIMISATION; IT IS THE STEP THAT MAKES D-05's
# MARGINAL-PRESERVATION CLAIM TRUE. A CAR covariance does NOT have a unit diagonal: at
# alpha = 0.95 the measured per-cell marginal sd runs 0.658 at an edge against 0.992 in the
# interior (12-RESEARCH.md Pattern 1, measured 2026-07-27). Feeding an unrescaled draw through
# the copula gives each region a DIFFERENT effective prior on rho -- a position-dependent
# prior, i.e. exactly the radial/edge artifact class this phase exists to detect. The rescale
# lives in `field_sampler` (in the DRAW), not in Sigma, so Sigma legitimately keeps its
# non-unit diagonal and only a draw can see the difference.
#
# ZERO NEW DEPENDENCIES (T-12-09). `LinearAlgebra` and `SparseArrays` are stdlibs, reachable
# through Julia's default LOAD_PATH `@stdlib` entry exactly as `run_p11_recovery.jl:50-53`
# reaches `LinearAlgebra`. Neither is in `spike/Project.toml` and neither may be added --
# `spike/test/runtests.jl:123-126` asserts the exact 16-name dependency set, and a resolve
# could dislodge the NeuralEstimators 0.2.1 pin. Nothing below installs a package, mutates
# the load path, or reaches for a third-party Gaussian-process / random-field library: at
# n = 64 the dense stdlib path is the whole job, and every such name is deliberately absent
# from this file so its absence stays greppable.
#
# Run (self-check):
#     julia --project=spike -e 'include("spike/simulator/p12_lattice.jl")'

using LinearAlgebra
using SparseArrays

# ORDER MATTERS: the Tier-1 pre-registration (P12_G, the r1 ladder and brackets, the GP
# jitter, the bisection tolerance) before any arithmetic that consumes it. Guarded for
# idempotency, the `run_p11_recovery.jl:55-56` idiom.
isdefined(@__MODULE__, :P12_G) || include(joinpath(@__DIR__, "..", "validation", "p12_consts.jl"))

# =====================================================================================
# 1. The indexing contract
# =====================================================================================

"""
    p12_idx(i, j; G = P12_G) -> Int

Flat index of lattice cell `(i, j)` in the COLUMN-MAJOR ordering `vec()` produces.

This is not a free choice. `src/amortized/summary.jl:72-76` encodes the `G×G` patch-correlation
matrix as `vec(coalesce.(M, 0.0))`, so summary row `(j-1)*G + i` is patch `(i, j)` and
`reshape(rows[1:G^2], G, G) == M` with no transpose. Every lattice construction in this file
routes through this one function so the contract is stated exactly once.
"""
p12_idx(i, j; G = P12_G) = (j - 1) * G + i

# =====================================================================================
# 2. The lattice and the two covariance arms
# =====================================================================================

"""
    lattice_adjacency(G = P12_G) -> SparseMatrixCSC{Float64}

The `G²×G²` 4-neighbour adjacency of the lattice, in `p12_idx` (column-major) order.
Symmetric, zero diagonal; each undirected edge appears as two stored entries.
"""
function lattice_adjacency(G::Integer = P12_G)
    G >= 1 || throw(ArgumentError("lattice_adjacency: G must be ≥ 1, got $G"))
    W = spzeros(Float64, G * G, G * G)
    for j in 1:G, i in 1:G, (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
        i2, j2 = i + di, j + dj
        (1 <= i2 <= G && 1 <= j2 <= G) || continue
        W[p12_idx(i, j; G = G), p12_idx(i2, j2; G = G)] = 1.0
    end
    return W
end

"""
    car_sigma(G, α) -> Symmetric{Float64}

The CAR arm: `Σ = (D − αW)⁻¹` with `W` the 4-neighbour adjacency and `D` its degree matrix.
`Q = D − αW` is SPD for every `0 ≤ α < 1` (strict diagonal dominance); `α = 1` is the improper
intrinsic limit and is excluded by the guard, matching `P12_CAR_ALPHA_BRACKET`.

TWO CONSEQUENCES OF `D` NOT BEING A MULTIPLE OF `I`, both load-bearing later:

1. Boundary cells have degree 2 (corners) or 3 (edges) against 4 in the interior, so the
   MARGINAL VARIANCES OF `Σ` ARE NON-UNIFORM — the reason `field_sampler` must rescale, and
   the reason the ablation cannot simply reuse `Σ` with α = 0 rounded off.
2. The 2-D DCT-II basis diagonalizes `Q` EXACTLY ONLY IN THE α → 1 LIMIT, where `Q` degenerates
   to the free-boundary graph Laplacian. For α < 1 it is an excellent but approximate
   diagonalization. This costs nothing here: `p12_dct_matrix` is used because it is an exact
   orthonormal, smoothness-ordered basis, which it is for every α — never because it is claimed
   to be the eigenbasis of the realized `Q` (12-RESEARCH.md Pattern 5, caveat stated there).
"""
function car_sigma(G::Integer, α::Real)
    (0 <= α < 1) || throw(ArgumentError(
        "car_sigma: α must satisfy 0 ≤ α < 1 (proper CAR; α = 1 is the improper intrinsic " *
        "limit and is excluded by P12_CAR_ALPHA_BRACKET), got $α"))
    W = lattice_adjacency(G)
    D = Diagonal(vec(sum(W; dims = 2)))
    Σ = inv(Symmetric(Matrix(D - α * W)))
    return Symmetric(Σ)
end

"""
    gp_sigma(G, ℓ; jitter = P12_GP_JITTER, kernel = :matern32) -> Symmetric{Float64}

The GP arm: a stationary correlation kernel evaluated on the lattice cell centres, in the SAME
column-major order as `p12_idx`, plus the pre-registered diagonal jitter.

THE DEFAULT KERNEL IS MATÉRN ν = 3/2, NOT SQUARED-EXPONENTIAL. Measured on this lattice
(12-RESEARCH.md Pattern 1, 2026-07-27) the SE kernel reaches `cond(K) = 5.5e9` at ℓ = 8 —
within a few orders of magnitude of losing the Cholesky to rounding at the top of
`P12_GP_ELL_BRACKET`. Matérn-3/2 is far better conditioned at the same effective range because
its spectral density decays polynomially rather than as a Gaussian. `:se` stays reachable by
keyword so the mini-spike can report both arms, but it is not the default.

`jitter` IS TIER-1 PRE-REGISTERED (`P12_GP_JITTER`) AND MUST NOT BE TUNED AFTER A FAILURE.
Raising jitter until a Cholesky succeeds is a pre-registration breach dressed as a numerical
fix, so a failed factorization THROWS with the failing `(ℓ, kernel, jitter)` named in the
message instead (T-12-16).
"""
function gp_sigma(G::Integer, ℓ::Real; jitter::Real = P12_GP_JITTER, kernel::Symbol = :matern32)
    G >= 1 || throw(ArgumentError("gp_sigma: G must be ≥ 1, got $G"))
    ℓ > 0   || throw(ArgumentError("gp_sigma: ℓ must be > 0, got $ℓ"))
    kernel in (:matern32, :se) || throw(ArgumentError(
        "gp_sigma: unknown kernel $kernel (expected :matern32, the default, or :se)"))
    P = [(float(i), float(j)) for j in 1:G for i in 1:G]   # column-major, matches p12_idx
    n = length(P)
    K = Matrix{Float64}(undef, n, n)
    @inbounds for b in 1:n, a in 1:n
        d2 = (P[a][1] - P[b][1])^2 + (P[a][2] - P[b][2])^2
        if kernel === :se
            K[a, b] = exp(-d2 / (2 * ℓ^2))
        else
            s = sqrt(3 * d2) / ℓ                            # Matérn 3/2: (1+s)exp(−s)
            K[a, b] = (1 + s) * exp(-s)
        end
    end
    Σ = Symmetric(K + jitter * I)
    F = cholesky(Σ; check = false)
    issuccess(F) || throw(ArgumentError(
        "gp_sigma: Cholesky failed at (ℓ = $ℓ, kernel = $kernel, jitter = $jitter). " *
        "P12_GP_JITTER is Tier-1 pre-registered: do NOT raise it to make this pass."))
    return Σ
end

# =====================================================================================
# 3. The draw, and the rescale that makes D-05 correct
# =====================================================================================

"""
    field_sampler(Σ) -> (rng -> Vector{Float64})

Return a closure drawing one lattice field from `Σ`, RESCALED to per-cell marginal sd 1:

    rng -> (L * randn(rng, n)) ./ sqrt.(diag(Σ))

THE `./ sd` IS THE POINT OF THIS FUNCTION. `Σ` does not have a unit diagonal — at CAR α = 0.95
the marginal sd measured 0.658 at an edge cell against 0.992 in the interior. Without the
rescale the copula in `p12_prior.jl` pushes a non-standard normal through `Φ`, so the per-region
`u` is not `Uniform(0,1)`, so the induced per-region prior on ρ is NOT `MU_PRIOR ∘ ĝ` and varies
with POSITION in the image: SIM-02 then fails at the boundary and nowhere else, and the phase
acquires the exact radial/edge artifact it exists to guard against.

The rescale is deliberately in the SAMPLER and not folded into `Σ`, because `induced_lag1`,
`car_alpha_for_r1` and `gp_ell_for_r1` all need the raw covariance and rescale it themselves;
a `Σ` that had already been normalized would make the "is the rescale applied?" question
untestable — the correlation matrix `Σ ./ (sd*sd')` has a unit diagonal for ANY SPD `Σ`,
whether or not any sampler ever divides. Only a DRAW can see the difference.

The `1e-12I` below is a numerical-symmetry epsilon on the factorization of an already-SPD
matrix (12-RESEARCH.md Pattern 1, verbatim). It is NOT the GP model jitter — that is
`P12_GP_JITTER`, applied inside `gp_sigma`, and it is pre-registered.
"""
function field_sampler(Σ::AbstractMatrix)
    sd = sqrt.(diag(Σ))
    all(>(0), sd) || throw(ArgumentError("field_sampler: Σ has a non-positive diagonal entry"))
    L = cholesky(Σ + 1e-12I).L
    n = size(Σ, 1)
    return rng -> (L * randn(rng, n)) ./ sd
end

# =====================================================================================
# 4. The one knob both arms are addressed by (R-8)
# =====================================================================================

"""
    induced_lag1(Σ; G = isqrt(size(Σ,1))) -> Float64

Mean Pearson correlation over all 4-neighbour pairs of the RESCALED covariance
`R = Σ ./ (sd * sd')`, i.e. the lag-1 spatial correlation a draw from `field_sampler(Σ)`
actually exhibits.

THIS IS THE SINGLE QUANTITY BOTH ARMS ARE PARAMETRIZED BY (R-8). CAR's α is a broken knob at
G = 8: measured induced r₁ is 0.136 at α = 0.5, 0.438 at α = 0.95, 0.692 at α = 0.99 and 0.947
at α = 0.999, so a uniform prior on α would put roughly 90 % of its mass on "no spatial
structure" and would make D-08's correlation length unidentifiable BY PARAMETRIZATION RATHER
THAN BY PHYSICS. Solving for α (or ℓ) at a requested r₁ makes the CAR-vs-GP comparison a
question about kernel SHAPE rather than about two incomparable coordinate systems.

Works on an EMPIRICAL covariance too, which is how the tests check a sampler rather than a
formula.
"""
function induced_lag1(Σ::AbstractMatrix; G::Integer = isqrt(size(Σ, 1)))
    G * G == size(Σ, 1) || throw(ArgumentError(
        "induced_lag1: Σ is $(size(Σ,1))×$(size(Σ,2)); expected G²×G² with G = $G"))
    sd = sqrt.(diag(Σ))
    R  = Σ ./ (sd * sd')
    W  = lattice_adjacency(G)
    rows, cols, _ = findnz(W)
    s = 0.0
    @inbounds for k in eachindex(rows)
        s += R[rows[k], cols[k]]
    end
    return s / length(rows)
end

"""
    _p12_bisect_r1(f, lo, hi, target, who) -> Float64

Bisect the monotone map `f : parameter → induced r₁` for `f(x) = target`, to
`P12_R1_SOLVE_TOL`. Asserts monotonicity ON THE BRACKET ENDPOINTS on entry (a silently
decreasing arm would otherwise return a plausible wrong number), and throws a message naming
the achievable range when `target` is out of reach rather than returning an endpoint.
"""
function _p12_bisect_r1(f, lo::Real, hi::Real, target::Real, who::AbstractString)
    a, b = float(lo), float(hi)
    fa, fb = f(a), f(b)
    fa < fb || throw(ArgumentError(
        "$who: induced_lag1 must be strictly increasing across the bracket ($lo, $hi), " *
        "but f(lo) = $fa and f(hi) = $fb"))
    (fa <= target <= fb) || throw(ArgumentError(
        "$who: requested r1 = $target is outside the achievable range [$fa, $fb] on the " *
        "pre-registered bracket ($lo, $hi)"))
    for _ in 1:200
        m  = 0.5 * (a + b)
        fm = f(m)
        abs(fm - target) <= P12_R1_SOLVE_TOL && return m
        fm < target ? (a = m) : (b = m)
        (b - a) <= 4 * eps(max(abs(a), abs(b))) && break   # bracket collapsed to float precision
    end
    return 0.5 * (a + b)
end

"""
    car_alpha_for_r1(r1; G = P12_G) -> Float64

The CAR α whose induced lag-1 correlation is `r1`, by bisection over
`P12_CAR_ALPHA_BRACKET` to `P12_R1_SOLVE_TOL`.
"""
car_alpha_for_r1(r1::Real; G::Integer = P12_G) =
    _p12_bisect_r1(α -> induced_lag1(car_sigma(G, α); G = G),
                   P12_CAR_ALPHA_BRACKET[1], P12_CAR_ALPHA_BRACKET[2], r1, "car_alpha_for_r1")

"""
    gp_ell_for_r1(r1; G = P12_G, kernel = :matern32) -> Float64

The GP correlation length ℓ whose induced lag-1 correlation is `r1`, by bisection over
`P12_GP_ELL_BRACKET` to `P12_R1_SOLVE_TOL`.
"""
gp_ell_for_r1(r1::Real; G::Integer = P12_G, kernel::Symbol = :matern32) =
    _p12_bisect_r1(ℓ -> induced_lag1(gp_sigma(G, ℓ; kernel = kernel); G = G),
                   P12_GP_ELL_BRACKET[1], P12_GP_ELL_BRACKET[2], r1, "gp_ell_for_r1")

"""
    lattice_sigma(arm::Symbol, r1; G = P12_G, kernel = :matern32) -> Symmetric{Float64}

The covariance of the requested arm AT the requested induced lag-1 correlation:
`:car`, `:gp`, or `:none`.

`:none` — and equally any arm at `r1 == P12_ABLATION_R1` — returns `Symmetric(Matrix(1.0I))`,
the D-10 matched-ablation covariance. The ablation is therefore a GENUINE DRAW FROM AN
INDEPENDENT FIELD taken through the same sampler, copula and datagen path as the spatial arms,
not a special case bolted on at read time: `field_sampler` of the identity draws i.i.d.
standard normals, so every downstream step is byte-for-byte the spatial code with r₁ = 0.
"""
function lattice_sigma(arm::Symbol, r1::Real; G::Integer = P12_G, kernel::Symbol = :matern32)
    arm in (:car, :gp, :none) || throw(ArgumentError(
        "lattice_sigma: unknown arm $arm (expected :car, :gp or :none)"))
    (arm === :none || r1 == P12_ABLATION_R1) && return Symmetric(Matrix(1.0I, G * G, G * G))
    arm === :car && return car_sigma(G, car_alpha_for_r1(r1; G = G))
    return gp_sigma(G, gp_ell_for_r1(r1; G = G, kernel = kernel); kernel = kernel)
end

# =====================================================================================
# 5. The deviation basis: 2-D DCT-II (D-01 and D-07 at the same time)
# =====================================================================================
# WHY A BASIS CHANGE IS NOT A D-07 VIOLATION. D-07 forbids scoring a DETERMINISTIC LOSSY
# transform of the sampled parameter -- the `ghat` clamp (which manufactures atoms) or a cell
# mean (a projection that discards information). The DCT-II matrix below is ORTHONORMAL and
# EXACTLY INVERTIBLE, so RANKING A DCT COEFFICIENT IS RANKING A SAMPLED PARAMETER: an
# orthogonal rotation discards nothing, and `p12_idct(p12_dct(F)) == F` to float tolerance.
# The field is sampled, expressed in this basis, and the coefficients ARE the theta rows.
#
# D-04's "high rank" therefore has a concrete meaning rather than a judgement call:
# K = G^2 - 1 = 63 (`P12_K_DEV`) -- every non-constant mode, i.e. NO TRUNCATION AT ALL, so the
# spike arm cannot be accused of having picked a flattering rank. Truncation is a PRODUCTION
# choice made later (12-15), against a measurable reconstruction error.

"""
    p12_dct_matrix(G = P12_G) -> Matrix{Float64}

The orthonormal 1-D DCT-II matrix, `C[p+1, i] = c_p · cos(π p (i − ½) / G)` with
`c_0 = √(1/G)` and `c_p = √(2/G)` for `p ≥ 1`, so `C'C == I` to float tolerance.
"""
function p12_dct_matrix(G::Integer = P12_G)
    G >= 1 || throw(ArgumentError("p12_dct_matrix: G must be ≥ 1, got $G"))
    C = Matrix{Float64}(undef, G, G)
    @inbounds for p in 0:(G - 1), i in 1:G
        c = p == 0 ? sqrt(1 / G) : sqrt(2 / G)
        C[p + 1, i] = c * cos(π * p * (i - 0.5) / G)
    end
    return C
end

"""
    p12_dct(F) -> Matrix{Float64}

Forward 2-D DCT-II of the `G×G` field `F`: `Ĉ = C F C'`. Row `p+1`, column `q+1` of the result
is the coefficient of the separable mode `(p, q)`. Parseval holds exactly (up to rounding):
`sum(abs2, p12_dct(F)) == sum(abs2, F)`.

The DC coefficient satisfies the identity `p12_dct(F)[1,1] == G · mean(F)` (because
`c_0² · Σᵢⱼ F = (1/G)·G²·mean(F)`), which is what makes the global θ row a clean rescaling of
the field mean rather than a separately-defined quantity.
"""
function p12_dct(F::AbstractMatrix)
    size(F, 1) == size(F, 2) || throw(ArgumentError(
        "p12_dct: expected a square G×G field, got $(size(F,1))×$(size(F,2))"))
    C = p12_dct_matrix(size(F, 1))
    return C * F * C'
end

"""
    p12_idct(Ĉ) -> Matrix{Float64}

Inverse 2-D DCT-II: `F = C' Ĉ C`. Exact inverse of `p12_dct` to float tolerance.
"""
function p12_idct(Ĉ::AbstractMatrix)
    size(Ĉ, 1) == size(Ĉ, 2) || throw(ArgumentError(
        "p12_idct: expected a square G×G coefficient array, got $(size(Ĉ,1))×$(size(Ĉ,2))"))
    C = p12_dct_matrix(size(Ĉ, 1))
    return C' * Ĉ * C
end

"""
    p12_dct_vec(z; G = isqrt(length(z))) -> Vector{Float64}

`p12_dct` on the COLUMN-MAJOR flat field vector, returning a column-major flat coefficient
vector, in **flat mode-index order**.

!!! warning "THIS IS NOT THE ORDER θ ROWS ARE STORED IN — ONE MORE STEP FOLLOWS"
    An earlier version of this docstring called the return value *"the form θ rows are actually
    stored in"*. **THAT WAS FALSE, and it is the sentence that licensed the wrong-basis defect**
    repaired in `f039729`. `p12_theta_column` (`spike/data/p12_generate.jl:193`) applies
    [`p12_dct_order`](@ref) to this output before storing it:

        c = p12_dct_vec(vec(draw.z_field); G = G)[p12_dct_order(G)]

    so **θ rows are in SMOOTHNESS order, this vector is in FLAT order, and the two differ** —
    `p12_dct_order(8)` begins `[1, 2, 9, 10, 3, 17, …]`, not `[1, 2, 3, …]`.

    To go from θ rows back to a field, use [`p12_region_field`](@ref) (`spike/p12/result.jl`),
    which inverts the permutation first. See the warning on [`p12_idct_vec`](@ref).
"""
function p12_dct_vec(z::AbstractVector; G::Integer = isqrt(length(z)))
    G * G == length(z) || throw(ArgumentError(
        "p12_dct_vec: length(z) = $(length(z)) is not G² for G = $G"))
    return vec(p12_dct(reshape(collect(float.(z)), G, G)))
end

"""
    p12_idct_vec(c; G = isqrt(length(c))) -> Vector{Float64}

Inverse of [`p12_dct_vec`](@ref), on the same column-major **flat mode-index** layout.

!!! danger "DO NOT PASS A θ COLUMN TO THIS FUNCTION — USE `p12_region_field`"
    **θ rows are NOT in this layout.** They are `p12_dct_vec(...)[p12_dct_order(G)]` —
    permuted into SMOOTHNESS order by `p12_theta_column` — and calling this function on them
    reconstructs a **DIFFERENT FIELD**, silently.

    **A θ column is indistinguishable from a flat coefficient vector by every property a
    caller can check**: same length `G²`, same element type, same plausible magnitudes. Nothing
    throws, nothing is `NaN`, and because a permutation is ORTHOGONAL the error is not a small
    perturbation — the reconstruction is wrong by O(1). Measured on a real draw: exact to
    `3.1e-15` through [`p12_region_field`](@ref), wrong by `3.17` through this function.

    Use **[`p12_region_field`](@ref)** (`spike/p12/result.jl`), which scatters each coefficient
    back through [`p12_dct_order`](@ref) *before* the inverse transform. This function is
    correct **only** for a vector that came out of `p12_dct_vec` and was never permuted — e.g.
    the truncation curve in `run_p12_minispike.jl`, which works in flat index space throughout.

    **WHY THIS WARNING EXISTS.** The sentence above it was *locally true and licensed a false
    inference*, which is worse than a missing check: a careful reader who verified the docstring
    was **confirmed in the error**. `p12_theta_column` documents that it permutes and
    `p12_region_field` documents that it inverts, but before `f039729` **neither named the
    other**, so the composition was undocumented while all three functions were individually
    correct. That gap produced the phase's second wrong-space defect
    (`12-STAGE1-VERDICT.md` §7.4).
"""
function p12_idct_vec(c::AbstractVector; G::Integer = isqrt(length(c)))
    G * G == length(c) || throw(ArgumentError(
        "p12_idct_vec: length(c) = $(length(c)) is not G² for G = $G"))
    return vec(p12_idct(reshape(collect(float.(c)), G, G)))
end

"""
    p12_dct_order(G = P12_G) -> Vector{Int}

The permutation of the `G²` flat coefficient indices that sorts the modes by SMOOTHNESS:
ascending total degree `p + q`, ties broken by `p² + q²` (the isotropic frequency of the
separable mode), then by flat index.

FROZEN ORDERING, SO "THE FIRST K MODES" IS A DEFINED OBJECT. The DC coefficient `(p,q) = (0,0)`
has total degree 0 and is therefore always index 1, which is what lets the global term be
`c₀` and the 63 deviation rows be `p12_dct_order(G)[2:end]` without any further convention.
"""
function p12_dct_order(G::Integer = P12_G)
    G >= 1 || throw(ArgumentError("p12_dct_order: G must be ≥ 1, got $G"))
    keyed = [(p + q, p^2 + q^2, p12_idx(p + 1, q + 1; G = G))
             for q in 0:(G - 1) for p in 0:(G - 1)]
    return Int[k[3] for k in sort(keyed)]
end

# =====================================================================================
# 6. Lattice geometry for the S-4 radial-energy guard
# =====================================================================================

"""
    radial_basis(G = P12_G) -> Matrix{Float64}

The `G²×2` design matrix `[1  radius]` over lattice cell centres, radius measured from the grid
centre `((G+1)/2, (G+1)/2)`, in `p12_idx` (column-major) order.

Lives here rather than in the guard runner so the S-4 radial-energy guard does not re-derive
the lattice geometry — a second derivation is a second chance to transpose it.
"""
function radial_basis(G::Integer = P12_G)
    G >= 1 || throw(ArgumentError("radial_basis: G must be ≥ 1, got $G"))
    c = (G + 1) / 2
    X = Matrix{Float64}(undef, G * G, 2)
    @inbounds for j in 1:G, i in 1:G
        r = hypot(i - c, j - c)
        X[p12_idx(i, j; G = G), 1] = 1.0
        X[p12_idx(i, j; G = G), 2] = r
    end
    return X
end
