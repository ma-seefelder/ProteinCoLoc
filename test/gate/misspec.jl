#########################################################################################
# This file is part of ProteinCoLoc.jl, licensed under the MIT License (MIT).
# See LICENSE.md in the project root for license information.
# Author: Dr. rer. nat. Manuel Seefelder
#########################################################################################
#
# test/gate/misspec.jl --- the per-grid OOD CONTROL SIMULATORS (D-03/D-04) for the ship-gate.
#
# `src/amortized/ood.jl` is deliberately the OOD READ SURFACE only; its header states that the
# "controlled misspecification-grid ROC EXPERIMENT (misspec_* families, negative controls,
# `ood_roc_over_grid`) is the per-grid D-05 SHIP-GATE's job". THIS is that file: it supplies the
# positive-control (misspecified) forward models and the summary-orthogonal negative controls the
# `--ood` arm of `run_gate.jl` needs in order to actually SCORE a separability AUC.
#
# PROMOTED (grid-generalized, physics/level semantics BYTE-UNCHANGED) from the proven Phase-5 spike:
#   spike/validation/ood.jl:392-410  _guard_misspec_imsize / _puncta_field
#   spike/validation/ood.jl:412-502  misspec_texture/noise/optics/background + OOD_FAMILIES
#   spike/validation/ood.jl:513-581  negctrl_affine / negctrl_rotate_flip / negctrl_block_permute
#                                    + verify_summary_invariance
#   spike/validation/ood.jl:612-699  ood_roc_over_grid  -> `gate_ood_roc` (density∨noise fusion)
#
# NOTHING IS RE-TUNED. The four families, their four magnitude levels, the level→strength maps and
# the negative-control transforms are the ESTABLISHED Phase-5 design, copied across. The ONLY
# deltas are the D-05 grid generalization:
#   (1) the patch grid `G` is a parameter (the spike hardcoded 8) — it drives the ≥64px imsize
#       guard, the frozen summary read chain and the negative-control block size,
#   (2) the frozen summary is read through the harness chain
#       `patch_summary(build_mci(pair), G) → encode_d01 → standardize_summary(·, m.zt, :min)`
#       (the loaded `m` bundle has no `:grid` field, so `ProteinCoLoc.frozen_summary` is not used),
#   (3) every inference is CPU-only.
#
# REPORTED DETECTOR = the OR-FUSED density∨noise channel (the Phase-5 iter2 design: the image-noise
# channel closes the detector-noise blind spot of a correlation-only summary). The POSTERIOR-
# PREDICTIVE channel is EXCLUDED from the reported fusion, exactly as the reported spike run did
# (`with_pp = false`): it costs OOD_PP_REPS re-simulations per scored sample, which is
# ~40k extra 256² simulations per grid — a resolution/cost decision, not a threshold change.
#
# NAMED BLIND SPOT (D-04, honesty capstone): the fixed G×G patch-Pearson summary is a CORRELATION
# statistic, so correlation-preserving transforms (affine exactly; rotate-flip / block-permute as a
# value multiset) are PROVABLY invisible to it. The negative controls MEASURE that blind spot; a
# near-0.5 separability on them is EXPECTED and reported, never hidden.

using ProteinCoLoc
import Statistics: mean, median, var, std
import Random: randperm
import ImageFiltering: imfilter, Kernel
import HypothesisTests: ApproximateTwoSampleKSTest

isdefined(@__MODULE__, :SBC_IMSIZE)      || include(joinpath(@__DIR__, "gate_consts_template.jl"))
isdefined(@__MODULE__, :draw_simulate_infer) || include(joinpath(@__DIR__, "harness.jl"))

# Median absolute deviation (robust scale) — hand-rolled, mirrors ood.jl's internal `_mad`.
_gate_mad(v) = (mv = median(v); median(abs.(v .- mv)))

"""
    gate_summary(m, pair, G) -> Vector

The FROZEN read chain for an arbitrary (possibly misspecified, fully external) 2-channel image
`pair`, grid-parametrized: `build_mci → patch_summary(·, G) → encode_d01 →
standardize_summary(·, m.zt, :min)`. Byte-identical to what `draw_simulate_infer` does to an
in-distribution image, so a misspecified image is scored with NO re-fitting (frozen-stats
discipline, Pitfall 5).
"""
gate_summary(m, pair, G::Integer) = ProteinCoLoc.standardize_summary(
    ProteinCoLoc.encode_d01(ProteinCoLoc.patch_summary(ProteinCoLoc.build_mci(pair), G)),
    m.zt, :min)

# ============================================================================
# Positive controls — the four misspecification families (D-03)
# ============================================================================
#
# All four produce FULLY EXTERNAL image pairs (SC5) the shared-latent smooth-Gaussian-field
# simulator CANNOT produce, each swept over OOD_GRID_LEVELS magnitudes for a dose-response ROC.

# imsize guard, grid-generalized: each patch must keep >15 px so patches are non-missing.
function _guard_misspec_imsize(imsize, G::Integer)
    (imsize[1] ≥ 64 && imsize[2] ≥ 64) ||
        throw(ArgumentError("misspec imsize dims must be ≥ 64 for the $(G)×$(G) patch grid, " *
                            "got $imsize"))
    (imsize[1] ÷ G ≥ 4 && imsize[2] ÷ G ≥ 4) ||
        throw(ArgumentError("misspec imsize $imsize gives <16 px per $(G)×$(G) patch"))
    return nothing
end

# Sparse Gaussian-puncta field: `nspots` bright impulses smoothed at scale σ — the building block
# of the texture family (sharp, sparse spots the SMOOTH-field simulator structurally cannot make).
function _puncta_field(rng, imsize, nspots::Integer, σ::Real)
    f = zeros(Float64, imsize...)
    H, W = imsize
    for _ in 1:nspots
        f[rand(rng, 1:H), rand(rng, 1:W)] += 0.5 + rand(rng)
    end
    return imfilter(f, Kernel.gaussian(σ))
end

"""
    misspec_texture(rng, θ; imsize = SBC_IMSIZE, level = OOD_GRID_LEVELS, G = 8)

Family 1 — TEXTURE-model mismatch. A NEW puncta/granular generator (shared + private Gaussian-spot
mixture, NOT `simulate_pair`): higher `level` ⇒ more, sharper spots ⇒ further from the smooth-field
training manifold. `θ.ρ_true` drives the shared spot weight (with `sign(ρ)` on channel 2) so the OOD
image still carries a colocalization structure. The canonical real-world OOD positive control.
"""
function misspec_texture(rng, θ; imsize = SBC_IMSIZE, level::Integer = OOD_GRID_LEVELS,
                         G::Integer = 8)
    _guard_misspec_imsize(imsize, G)
    ρ = clamp(θ.ρ_true, -1.0, 1.0)
    σ = max(0.8, 3.0 - 0.5 * level)              # sharper (more un-smooth) at higher level
    a = sqrt(abs(ρ)); b = sqrt(1.0 - abs(ρ))
    shared = _puncta_field(rng, imsize, 30 * level, σ)
    p1     = _puncta_field(rng, imsize, 20 * level, σ)
    p2     = _puncta_field(rng, imsize, 20 * level, σ)
    ch1 = ProteinCoLoc.BG_FLOOR .+ a .* shared .+ b .* p1
    ch2 = ProteinCoLoc.BG_FLOOR .+ sign(ρ) .* a .* shared .+ b .* p2
    return [Matrix{Float64}(max.(ch1, 0.0)), Matrix{Float64}(max.(ch2, 0.0))]
end

"""
    misspec_noise(rng, θ; imsize = SBC_IMSIZE, level = OOD_GRID_LEVELS, G = 8)

Family 2 — NOISE-model mismatch. Corrupts a `simulate_pair` base with salt-and-pepper spikes plus
heavy-tailed (cubic) additive noise BEYOND the simulator's Poisson+Gaussian model; `level` scales
both the spike fraction and the amplitude.
"""
function misspec_noise(rng, θ; imsize = SBC_IMSIZE, level::Integer = OOD_GRID_LEVELS,
                       G::Integer = 8)
    _guard_misspec_imsize(imsize, G)
    base  = ProteinCoLoc.simulate_pair(rng, θ; imsize = imsize)
    frac  = 0.02 * level
    scale = 2.0 * level
    out = map(base) do ch
        c  = copy(ch)
        mx = maximum(c)
        @inbounds for idx in eachindex(c)
            rand(rng) < frac && (c[idx] = rand(rng) < 0.5 ? 0.0 : mx * (1 + scale))
        end
        c .+= scale .* mx .* (rand(rng, size(c)...) .- 0.5) .^ 3   # heavy-tailed, non-Gaussian
        Matrix{Float64}(max.(c, 0.0))
    end
    return collect(out)
end

"""
    misspec_optics(rng, θ; imsize = SBC_IMSIZE, level = OOD_GRID_LEVELS, G = 8)

Family 3 — OPTICS/PSF mismatch. Convolves a `simulate_pair` base with a strongly ANISOTROPIC /
aberrated Gaussian PSF (σy ≫ σx) far beyond the fixed isotropic `σ_psf`; `level` scales the
anisotropy the simulator cannot produce.
"""
function misspec_optics(rng, θ; imsize = SBC_IMSIZE, level::Integer = OOD_GRID_LEVELS,
                        G::Integer = 8)
    _guard_misspec_imsize(imsize, G)
    base = ProteinCoLoc.simulate_pair(rng, θ; imsize = imsize)
    σy   = ProteinCoLoc.σ_psf + 2.0 * level       # ≫ σ_psf — aberrated/out-of-focus
    σx   = ProteinCoLoc.σ_psf + 0.2 * level
    return [Matrix{Float64}(max.(imfilter(ch, Kernel.gaussian((σy, σx))), 0.0)) for ch in base]
end

"""
    misspec_background(rng, θ; imsize = SBC_IMSIZE, level = OOD_GRID_LEVELS, G = 8)

Family 4 — BACKGROUND/illumination mismatch. Applies a multiplicative illumination gradient +
radial vignette and an autofluorescence bleed BEYOND the simulator's `Uniform(0, 0.1)` offset;
`level` scales the gradient, vignette and bleed.
"""
function misspec_background(rng, θ; imsize = SBC_IMSIZE, level::Integer = OOD_GRID_LEVELS,
                            G::Integer = 8)
    _guard_misspec_imsize(imsize, G)
    base = ProteinCoLoc.simulate_pair(rng, θ; imsize = imsize)
    H, W = imsize
    gy   = reshape(range(-1.0, 1.0; length = H), H, 1)   # H×1
    gx   = reshape(range(-1.0, 1.0; length = W), 1, W)   # 1×W
    grad  = 1.0 .+ (0.5 * level) .* gy .+ (0.3 * level) .* gx           # H×W
    vign  = 1.0 .- (0.2 * level) .* (gy .^ 2 .+ gx .^ 2)               # H×W
    bleed = 0.1 + 0.3 * level                                          # beyond U(0, 0.1)
    return [Matrix{Float64}(max.(ch .* max.(grad, 0.0) .* max.(vign, 0.0) .+ bleed, 0.0))
            for ch in base]
end

# The four positive-control families, in the FIXED Phase-5 order (the reported ROC grid, D-03).
const OOD_FAMILIES = (texture    = misspec_texture,
                      noise      = misspec_noise,
                      optics     = misspec_optics,
                      background = misspec_background)

"""
    misspec_simulator(gen; level, G) -> NamedTuple

Wrap a positive-control generator into the `(sample_prior, simulate_pair, build_mci)` SIMULATOR
CONTRACT `draw_simulate_infer` / `run_gate(...; pos_sim = …)` consume, at a fixed misspecification
`level`. The prior stays π(θ) (the misspecification is in the FORWARD MODEL, not the prior — that
is what "model misspecification" means).
"""
misspec_simulator(gen; level::Integer = OOD_GRID_LEVELS, G::Integer = 8) = (
    sample_prior  = ProteinCoLoc.sample_prior,
    simulate_pair = (rng, θ; imsize = SBC_IMSIZE) -> gen(rng, θ; imsize = imsize,
                                                         level = level, G = G),
    build_mci     = ProteinCoLoc.build_mci,
)

# ============================================================================
# Summary-orthogonal NEGATIVE controls (D-04) — the named blind spot
# ============================================================================

"""
    negctrl_affine(pair; a = 2.0, b = 0.5)

Per-channel positive affine `a·x + b` (a > 0) on the SIGNAL pixels (x > 0), leaving exact-zero
pixels at zero. Pearson is exactly scale/shift invariant, so the summary is UNCHANGED — the
strongest theoretical guarantee. Preserving the zero set is load-bearing (a blanket `+b` would turn
read-noise zeros into signal and change which pixels `_exclude_zero` drops per patch).
"""
function negctrl_affine(pair; a::Real = 2.0, b::Real = 0.5)
    a > 0 || throw(ArgumentError("negctrl_affine: slope a must be > 0, got $a"))
    return [Matrix{Float64}(map(v -> v > 0 ? a * v + b : v, ch)) for ch in pair]
end

"""
    negctrl_rotate_flip(pair; k = 1)

Global 90°·k rotation (no interpolation). On a square image this permutes WHICH patch is where but
preserves the SET of per-patch correlation values, so the summary DISTRIBUTION is unchanged.
"""
function negctrl_rotate_flip(pair; k::Integer = 1)
    rot(x) = (kk = mod(k, 4); kk == 0 ? x : kk == 1 ? rotl90(x) : kk == 2 ? rot180(x) : rotr90(x))
    return [Matrix{Float64}(rot(ch)) for ch in pair]
end

"""
    negctrl_block_permute(rng, pair; blocks = 8)

Distribution-preserving spatial rearrangement: partition into `blocks × blocks` tiles and apply the
SAME random tile permutation to BOTH channels. Each patch's paired pixels move together, so the
per-patch correlation values are preserved as a SET. The gate calls it with `blocks = G` so the
tiles align with the grid's patches (the spike's `blocks = 8` on its 8×8 grid, generalized).
"""
function negctrl_block_permute(rng, pair; blocks::Integer = 8)
    H, W = size(pair[1])
    # NON-DIVISIBLE SIZES (07-GATE-AMENDMENT §4, F5 mixture). v1 ran this control at the single
    # size (256,256), which is divisible by every gate grid, so an exact-divisibility `throw` was
    # sufficient. The amended pre-registration draws from a REALISTIC mixture that includes the
    # D-08 anchor 1376×1028, and 1028 is NOT divisible by 8 — the control would simply throw on a
    # quarter of the draws. The transform is therefore generalized: the largest ALIGNED
    # `blocks × blocks` region of equal tiles is permuted and the ≤ blocks−1 leftover rows/columns
    # are left in place. This is the minimal generalization that keeps the transform a genuine
    # rearrangement of whole patches (the property that makes it summary-orthogonal); the residual
    # strip is measured, not assumed away, because `verify_summary_invariance` still scores the
    # realised KS statistic of the whole transformed image.
    bh = H ÷ blocks; bw = W ÷ blocks
    (bh >= 1 && bw >= 1) ||
        throw(ArgumentError("negctrl_block_permute: $((H, W)) smaller than blocks=$blocks"))
    perm = randperm(rng, blocks * blocks)             # SAME permutation for both channels
    function permute(x)
        out = copy(x)                                  # leftover strip stays exactly as it was
        for (dst, src) in enumerate(perm)
            di, dj = fldmod1(dst, blocks); si, sj = fldmod1(src, blocks)
            out[(di-1)*bh+1:di*bh, (dj-1)*bw+1:dj*bw] = x[(si-1)*bh+1:si*bh, (sj-1)*bw+1:sj*bw]
        end
        return out
    end
    return [Matrix{Float64}(permute(ch)) for ch in pair]
end

"""
    verify_summary_invariance(pair, transform, G) -> Float64

Empirically verify (D-04) that `transform` leaves the `G×G` patch-Pearson summary statistically
UNCHANGED: the two-sample KS statistic (`ApproximateTwoSampleKSTest.δ`) between the multisets of
non-missing correlation values before/after. A correlation-preserving transform gives δ ≈ 0; call
sites assert `δ < OOD_KS_EPS`. Returns `NaN` if either summary is fully missing.
"""
function verify_summary_invariance(pair, transform, G::Integer)
    _vals(p) = Float64.(collect(skipmissing(
        ProteinCoLoc.patch_summary(ProteinCoLoc.build_mci(p), G))))
    v0 = _vals(pair)
    v1 = _vals(transform(pair))
    (isempty(v0) || isempty(v1)) && return NaN
    return ApproximateTwoSampleKSTest(v0, v1).δ
end

"""
    gate_negctrls(G; rng_factory) -> NamedTuple

The three summary-orthogonal negative-control transforms for grid `G` (block size = `G`, so the
tiles align with the patch grid). `rng_factory()` supplies the RNG the block permutation draws from.
"""
gate_negctrls(G::Integer; rng_factory) = (
    affine = (p -> negctrl_affine(p; a = 2.0, b = 0.5)),
    rotate = (p -> negctrl_rotate_flip(p; k = 1)),
    block  = (p -> negctrl_block_permute(rng_factory(), p; blocks = G)),
)

# ============================================================================
# The controlled ROC over the misspecification grid (D-03/D-04/D-06)
# ============================================================================

"""
    gate_ood_roc(m, density; G, rng, families = OOD_FAMILIES, levels = OOD_GRID_LEVELS,
                 n_id = 200, n_pos = 40, n_fit = n_id, imsize = SBC_IMSIZE) -> NamedTuple

The per-grid controlled ROC/AUC experiment — the grid-generalized port of the spike's
`ood_roc_over_grid` (`with_pp = false`, the reported Phase-5 configuration).

Draws, SEQUENTIALLY off the fresh disjoint per-grid stream `rng`:
  (A) `n_fit` TRAIN-ONLY ID pairs → the image-noise Mahalanobis null + the per-channel robust-z
      reference (median/MAD) that puts density and noise on a comparable footing for OR-fusion.
      Fit on ID (training-distribution) draws ONLY — misspec data is external and only ever SCORED.
  (B) `n_id` fresh ID negatives → the density ID operating point `id_threshold` (`OOD_ID_QUANTILE`,
      the PRE-REGISTERED gate operating point) and the fused operating point.
  (C) for each family × level, `n_pos` fully-external positives, scored through the FROZEN chain.

`fused_auc[fam][lvl]` is the REPORTED (gated) AUC of the OR-fused density∨noise detector;
`density_auc` / `noise_auc` are the standalone channels, returned for transparency. `combined_auc`
is the pooled strongest-level fused AUC (the spike's headline number). CPU-only.
"""
function gate_ood_roc(m, density; G::Integer, rng, families = OOD_FAMILIES,
                      levels::Integer = OOD_GRID_LEVELS, n_id::Integer = 200,
                      n_pos::Integer = 40, n_fit::Integer = n_id, imsize = SBC_IMSIZE,
                      imsize_set = GATE_IMSIZE_SET, imsize_weights = GATE_IMSIZE_WEIGHTS)
    sp(r) = ProteinCoLoc.sample_prior(r)
    # F5 (07-GATE-AMENDMENT §4): the OOD arm's ID pool must come from the SAME joint the net was
    # trained on, so the image size is drawn PER PAIR off this rng through the shared sampler. The
    # realised sizes are recorded and returned, never reconstructed.
    imsizes = Tuple{Int,Int}[]
    function idpair(r)
        θ   = sp(r)
        isz = gate_imsize(r; imsize = imsize, imsize_set = imsize_set,
                          imsize_weights = imsize_weights)
        push!(imsizes, isz)
        return ProteinCoLoc.simulate_pair(r, θ; imsize = isz)
    end

    # --- (A) TRAIN-ONLY noise null + robust-z fusion reference -------------------------------
    fit_pairs = [idpair(rng) for _ in 1:n_fit]
    nnull     = ProteinCoLoc.fit_noise_null(reduce(hcat, [ProteinCoLoc.noise_features(p)
                                                          for p in fit_pairs]))
    fit_maha  = [ProteinCoLoc.maha_score(density, gate_summary(m, p, G)) for p in fit_pairs]
    fit_noise = [ProteinCoLoc.noise_score(nnull, p) for p in fit_pairs]
    med_d = median(fit_maha);  scl_d = _gate_mad(fit_maha)  + 1e-12
    med_n = median(fit_noise); scl_n = _gate_mad(fit_noise) + 1e-12
    # OR-fusion (D-05): per-sample max of the two channels' robust-z scores.
    zfuse(sd, sn) = max((sd - med_d) / scl_d, (sn - med_n) / scl_n)

    # --- (B) ID negatives + the PRE-REGISTERED operating points ------------------------------
    id_pairs = [idpair(rng) for _ in 1:n_id]
    id_maha  = [ProteinCoLoc.maha_score(density, gate_summary(m, p, G)) for p in id_pairs]
    id_noise = [ProteinCoLoc.noise_score(nnull, p) for p in id_pairs]
    id_fused = [zfuse(id_maha[k], id_noise[k]) for k in 1:n_id]
    maha_thr  = ProteinCoLoc.id_threshold(id_maha;  q = OOD_ID_QUANTILE)
    fused_thr = ProteinCoLoc.id_threshold(id_fused; q = OOD_ID_QUANTILE)

    # --- (C) the family × level positive-control grid ----------------------------------------
    fam_names   = keys(families)
    fused_auc   = Dict{Symbol,Vector{Float64}}()
    density_auc = Dict{Symbol,Vector{Float64}}()
    noise_auc   = Dict{Symbol,Vector{Float64}}()
    fire_rate   = Dict{Symbol,Vector{Float64}}()
    pooled_pos  = Float64[]                       # strongest level, all families (headline ROC)

    for fam in fam_names
        gen = families[fam]
        fused_auc[fam] = Float64[]; density_auc[fam] = Float64[]
        noise_auc[fam] = Float64[]; fire_rate[fam]   = Float64[]
        for lvl in 1:levels
            pm = Float64[]; pn = Float64[]; pf = Float64[]; fires = 0
            for _ in 1:n_pos
                θ   = sp(rng)
                isz = gate_imsize(rng; imsize = imsize, imsize_set = imsize_set,
                                  imsize_weights = imsize_weights)
                push!(imsizes, isz)
                img = gen(rng, θ; imsize = isz, level = lvl, G = G)
                ms  = ProteinCoLoc.maha_score(density, gate_summary(m, img, G))
                ns  = ProteinCoLoc.noise_score(nnull, img)
                fs  = zfuse(ms, ns)
                push!(pm, ms); push!(pn, ns); push!(pf, fs)
                fires += (fs > fused_thr) ? 1 : 0
                lvl == levels && push!(pooled_pos, fs)
            end
            push!(fused_auc[fam],   ProteinCoLoc.roc_auc(id_fused, pf)[3])   # the GATED AUC
            push!(density_auc[fam], ProteinCoLoc.roc_auc(id_maha,  pm)[3])
            push!(noise_auc[fam],   ProteinCoLoc.roc_auc(id_noise, pn)[3])
            push!(fire_rate[fam],   fires / n_pos)
        end
    end

    fpr, tpr, combined = ProteinCoLoc.roc_auc(id_fused, pooled_pos)
    return (levels = collect(1:levels), families = collect(fam_names),
            fused_auc = fused_auc, density_auc = density_auc, noise_auc = noise_auc,
            fire_rate = fire_rate, combined_auc = combined,
            id_threshold = maha_thr, fused_threshold = fused_thr,
            id_fire_rate = mean(id_fused .> fused_thr),
            youden = ProteinCoLoc.youden_j(fpr, tpr),          # POST-HOC reference, NEVER the gate
            noise_null = nnull, imsizes = imsizes,
            zref = (med_d = med_d, scl_d = scl_d, med_n = med_n, scl_n = scl_n),
            n_id = n_id, n_pos = n_pos, n_fit = n_fit)
end
