---
phase: 02-forward-simulator-summary-contract
verified: 2026-06-26T00:00:00Z
status: passed
score: 4/4 must-haves verified
has_blocking_gaps: false
overrides_applied: 0
---

# Phase 2: Forward Simulator + Summary Contract Verification Report

**Phase Goal:** A physics forward-simulator emits real MultiChannelImage pairs from θ that
the *unmodified* existing summary functions ingest, with a prior provably consistent with the
Turing model so every downstream comparison stays valid.
**Verified:** 2026-06-26
**Status:** passed
**Re-verification:** No — initial verification

---

## Hard Gate Result

`julia --project=spike spike/test/runtests.jl` — **exit 0**

| Testset | Pass | Total | Time |
|---------|------|-------|------|
| NeuralEstimators CPU smoke | 5 | 5 | 3m38.9s |
| SIM-03 summary contract | 14 | 14 | 2.2s |
| SIM-01 forward pipeline | 12 | 12 | 2.0s |
| SIM-03 on simulator output | 10 | 10 | 0.0s |
| D-15 monotone ρ_true → induced-correlation | 5 | 5 | 0.4s |
| SIM-02 prior consistency | 10 | 10 | 37.8s |
| SIM-04 plausibility | 21 | 21 | 4.3s |
| **TOTAL** | **77** | **77** | **~4m26s** |

The SIM-02 testset runs N=300 pairs at 512² (the metric_N calibration budget), driving the
~37.6 s subtimer. Suite did not cancel.

Decoupling invariant: `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` → **exit 0**
(`DECOUPLING: root untouched`)

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | SIM-01: simulate_pair(rng,θ;imsize) runs all 7 D-06 stages returning Vector{Matrix{Float64}} length 2 from an explicit rng::AbstractRNG | VERIFIED | spike/simulator/forward.jl — 187 lines, all 7 stages in order (shared-latent→thinning→PSF→spillover→autofluor→warp-shift→Poisson+Gaussian), returns `[Matrix{Float64}(max.(ch1,0.0)), Matrix{Float64}(max.(ch2,0.0))]`; "SIM-01 forward pipeline" testset 12/12 |
| 2 | SIM-01: θ is the 7-tuple (ρ_true, spillover, autofluorescence, label_efficiency, shift_dx, shift_dy, noise); PSF is a fixed Gaussian nuisance NOT a θ component | VERIFIED | forward.jl defines `θ_BASE = (ρ_true=..., spillover=..., autofluorescence=..., label_efficiency=..., shift_dx=..., shift_dy=..., noise=...)`; PSF width is module constant `σ_psf=1.3` not read from θ |
| 3 | SIM-01: Stage 1 is spike-validated shared-latent correlated smooth Gaussian-field generator with sign(ρ) flip enabling anti-correlation | VERIFIED | forward.jl lines 146–147: `ch1=_softplus.(a.*L .+ b.*ε1)`, `ch2=_softplus.(sign(ρ).*a.*L .+ b.*ε2)`; softplus/log1p present; "D-15 monotone" testset asserts μ_ind[1]<0 at ρ=-0.7; 5/5 pass |
| 4 | SIM-02: sample_prior's induced μ matches Turing μ-prior Truncated(Cauchy(0,0.3),-1,1) within pre-declared KS/Wasserstein tolerance; calibration is a real measured sweep+fit | VERIFIED | ghat.jl header: W1=0.052<0.10 (tol), KS=0.101<0.20 (tol), sweep ρ_true∈[-0.99,0.99] 25pts n_per=100 512²; Spearman(ρ_grid,E[μ])=1.0; test "induced μ matches..." 1/1 pass at N=300 512² |
| 5 | SIM-02: ĝ spans past ±0.9 because generator induces only μ≈±0.76 at knob ±0.9 (D-15 rationale documented) | VERIFIED | ghat.jl: GHAT_MU_KNOTS spans [-0.680,+0.847]; GHAT_RHO_KNOTS endpoints are ±0.99; NOTES.md §3 documents D-15 rationale explicitly |
| 6 | SIM-02: Real anchor extracted faithfully via LoadImages.jl load_tiff (NOT RGB→luminance); negative tail documented PRIOR-ONLY | VERIFIED | calibration.jl _real_anchor() calls `load_tiff(joinpath(_REPO_ROOT,"test","test_images",cond,"$(cond)_c1.tif"))`; ghat.jl: positive μ=0.3292, negative μ=0.2481 (both positive → "negative tail: physically reachable = false ⇒ PRIOR-ONLY"); NOTES.md §3 §"Real anchor (D-16)" |
| 7 | SIM-02: σ/τ/ν consistency checked against Turing scale/tail hyperpriors | VERIFIED | ghat.jl: pooled SD=0.317 (inside Truncated(Cauchy(0.1,0.3),1e-4,1)), excess kurtosis=-0.052 (near-Gaussian, consistent with large-ν Exponential); NOTES.md §3 §"σ/τ/ν consistency" |
| 8 | SIM-03: simulator output builds valid MultiChannelImage flowing through UNCHANGED patch()/correlation() via read-only include() | VERIFIED | contract.jl includes src/LoadImages.jl and src/colocalization.jl read-only; no function named `patch` or `correlation` in spike/ (grep returns zero matches); "SIM-03 on simulator output" 10/10 |
| 9 | SIM-04: ρ_true sweep shows monotone Spearman ≥ threshold; spillover raises and shift lowers correlation as paired effects; output is valid MultiChannelImage | VERIFIED | "SIM-04 plausibility" 21/21: Spearman 1.0≥0.95 over 15-pt sweep; paired spillover mean(Δ)>0.02 p<0.05; paired shift mean(Δ)>0.01 p<0.05; perturbation runs valid MCI ≤2 missing |
| 10 | Decoupling: root Project.toml/Manifest.toml/src/ byte-identical to baseline f581d95 | VERIFIED | `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` exits 0 |

**Score:** 4/4 SIM requirements verified (10 underlying truths, all VERIFIED)

---

### Deferred Items

None.

---

## Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `spike/contract.jl` | include() coupling boundary + build_mci/patch_summary/induced_mu helpers | VERIFIED | 108 lines; includes src/LoadImages.jl and src/colocalization.jl; defines build_mci, patch_summary (not Base.summary — CR-01 fix), induced_mu with NaN guard (CR-02 fix) |
| `spike/simulator/forward.jl` | simulate_pair(rng,θ;imsize) — 7-stage D-06 pipeline | VERIFIED | 187 lines; all 7 stages; shared-latent correlated fields (stage 1); PSF as fixed constant; explicit rng threading; θ/imsize validation (≥64 dims, finite, in-range nuisances) |
| `spike/simulator/calibration.jl` | offline ρ_true sweep + PAVA isotonic fit + ĝ inverse + evidence + faithful anchor | VERIFIED | 347 lines; PAVA isotonic regression implemented; _real_anchor() uses load_tiff; emit_ghat() emits frozen artifact; corspearman used for monotonicity metric |
| `spike/simulator/ghat.jl` | Frozen ĝ knots/coefficients + ghat(μ) evaluator + GHAT_MU_MIN/MAX | VERIFIED | Auto-generated; 25 knot pairs GHAT_MU_KNOTS[-0.680,+0.847]/GHAT_RHO_KNOTS[-0.99,+0.99]; ghat(μ) evaluator; calibration evidence header |
| `spike/simulator/prior.jl` | MU_PRIOR const + sample_prior(rng) using frozen ĝ | VERIFIED | MU_PRIOR=Truncated(Cauchy(0.0,0.3),-1.0,1.0); sample_prior draws μ*~MU_PRIOR then ρ_true=ghat(μ*); 6 nuisances from Uniform priors; SIM02_W1_TOL=0.10 |
| `spike/test/test_simulator.jl` | SIM-01/02/03/04 testsets | VERIFIED | 304 lines; 5 testsets (SIM-03, SIM-01, SIM-03-on-output, D-15, SIM-02, SIM-04); all WR fixes applied |
| `spike/02_simulator_demo.jl` | Seeded headless demo regenerating plausibility figures | VERIFIED | 156 lines; CairoMakie.activate!(); 4 independent per-sweep RNGs (WR-03 fix); @assert isfile self-check on png+pdf |
| `spike/figures/plausibility.png` | Rendered SIM-04 evidence (raster) | VERIFIED | File exists at spike/figures/plausibility.png |
| `spike/figures/plausibility.pdf` | Rendered SIM-04 evidence (vector) | VERIFIED | File exists at spike/figures/plausibility.pdf |
| `spike/NOTES.md §3` | SIM-02 calibration documentation | VERIFIED | §3 documents ĝ method, span-past-±0.9 rationale (D-15), W1/KS evidence, size-invariance (max|Δ|=0.032), real anchor + negative-tail PRIOR-ONLY verdict (D-16), nuisance ranges, σ/τ/ν check |

---

## Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| spike/contract.jl | src/colocalization.jl | include() read-only | VERIFIED | Line 47: `include(joinpath(@__DIR__, "..", "src", "colocalization.jl"))` |
| spike/contract.jl | src/LoadImages.jl | include() read-only | VERIFIED | Line 46: `include(joinpath(@__DIR__, "..", "src", "LoadImages.jl"))` |
| spike/test/runtests.jl | spike/test/test_simulator.jl | include() | VERIFIED | Line 71: `include(joinpath(@__DIR__, "test_simulator.jl"))` |
| spike/simulator/forward.jl | spike/contract.jl | output Vector{Matrix{Float64}} feeds build_mci/patch_summary | VERIFIED | simulate_pair tested via build_mci in "SIM-03 on simulator output" testset |
| spike/simulator/forward.jl | ImageFiltering.imfilter / Kernel.gaussian | fixed Gaussian PSF + shared-latent field smoothing | VERIFIED | Lines 85, 158–159: imfilter + Kernel.gaussian in _smooth_field and PSF stage |
| spike/simulator/prior.jl | spike/simulator/ghat.jl | include() frozen calibration map | VERIFIED | Line 37: `include(joinpath(@__DIR__, "ghat.jl"))` |
| spike/simulator/calibration.jl | spike/contract.jl | induced_mu measured through the real summary | VERIFIED | calibration.jl line 56 includes contract.jl; calls induced_mu() |
| spike/test/test_simulator.jl | spike/simulator/forward.jl + prior.jl + contract.jl | SIM-04 quantitative gate via corspearman | VERIFIED | corspearman appears at lines 191, 260 for ĝ-monotonicity and SIM-04(a) gates |
| spike/02_simulator_demo.jl | spike/figures/ | CairoMakie save() | VERIFIED | Lines 145–146: save(PNG_PATH, fig), save(PDF_PATH, fig); @assert isfile verifies both |

---

## Data-Flow Trace (Level 4)

The core simulator chain is: `sample_prior(rng) → θ → simulate_pair(rng,θ) → Vector{Matrix{Float64}} → build_mci → MultiChannelImage → patch_summary → 8×8 ρ → induced_mu → Float64`.

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|--------------------|--------|
| spike/contract.jl `patch_summary` | 8×8 ρ matrix | src/colocalization.jl `patch()`+`correlation()` called on mci.data[1]/[2] | Yes — real patch computation on channel arrays | FLOWING |
| spike/contract.jl `induced_mu` | Float64 mean | `mean(skipmissing(patch_summary(mci)))` | Yes — mean of real per-patch correlations | FLOWING |
| spike/simulator/forward.jl `simulate_pair` | ch1, ch2 matrices | shared-latent smooth fields + 7 physics stages | Yes — all entries finite ≥ 0, small-positive background | FLOWING |
| spike/simulator/prior.jl `sample_prior` | θ NamedTuple | ghat(rand(rng, MU_PRIOR)) + Uniform nuisance draws | Yes — ρ_true from frozen ĝ, all 6 nuisances sampled | FLOWING |
| spike/simulator/ghat.jl `ghat` | ρ_true Float64 | GHAT_MU_KNOTS/GHAT_RHO_KNOTS literal arrays from offline calibration | Yes — 25-point measured sweep, not asserted range copy | FLOWING |

---

## Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| Full test suite exits 0 | `julia --project=spike spike/test/runtests.jl` | 77/77 passing, exit 0 | PASS |
| Decoupling invariant holds | `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` | exit 0 | PASS |
| NeuralEstimators stays v0.2.1 | UUID check via Pkg.dependencies() | v0.2.1 confirmed | PASS |
| No CUDA installed | filter on CUDA name in Pkg.dependencies() | 0 CUDA deps found | PASS |
| No patch/correlation reimplementation | grep `function (patch|correlation)` in spike/ | 0 matches | PASS |
| No forbidden packages | grep NormalizingFlows/InvertibleNetworks/RxInfer/BayesFlow in spike/ | 0 matches | PASS |

---

## Probe Execution

No phase-specific probes declared. The hard gate is `julia --project=spike spike/test/runtests.jl` (run above, exit 0).

---

## Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| SIM-01 | 02-01, 02-02 | simulate_pair(θ) → MultiChannelImage with 7-stage pipeline, explicit rng, 7-tuple θ, shared-latent stage-1 generator (D-15) | SATISFIED | spike/simulator/forward.jl verified; "SIM-01 forward pipeline" 12/12 |
| SIM-02 | 02-03 | Prior π(θ) consistent with Turing @model μ-prior via measured induced-μ calibration; documented in NOTES.md §3 | SATISFIED | ghat.jl: W1=0.052<0.10, KS=0.101<0.20; "SIM-02 prior consistency" 10/10; NOTES.md §3 complete |
| SIM-03 | 02-01, 02-02 | Output verifies as valid MultiChannelImage so unchanged correlation()/patch() apply | SATISFIED | contract.jl read-only include(); no patch/correlation reimplementation; "SIM-03 on simulator output" 10/10 |
| SIM-04 | 02-04 | Plausibility plots confirm ρ_true↑→patch-corr↑; spillover and shift visibly affect the pair | SATISFIED | "SIM-04 plausibility" 21/21; Spearman 1.0≥0.95; paired effects significant; figures at spike/figures/ |

All four phase-2 requirements are SATISFIED. No orphaned requirements.

---

## Anti-Patterns Found

No TBD, FIXME, XXX, TODO, HACK, or PLACEHOLDER markers in any spike/ .jl file (grep returns zero matches).

Code review issues CR-01/CR-02 and WR-01..WR-07 identified in 02-REVIEW.md were all resolved in the committed code:

| Issue | Severity | Resolution confirmed |
|-------|----------|----------------------|
| CR-01: Base.summary extension wrong return type | Critical | Fixed — contract.jl uses `patch_summary(mci)` not `Base.summary` |
| CR-02: induced_mu crash on all-missing summary | Critical | Fixed — forward.jl validates imsize≥64; contract.jl `induced_mu` returns NaN on empty collection |
| WR-01: calibrate() mean/std on empty μs | Warning | Fixed — calibration.jl has `if isempty(μs)` guard with `@warn` |
| WR-02: D-15 monotone test one simulation per point | Warning | Fixed — test averages 3 replicates per ρ via `for i in 0:2` |
| WR-03: Shared mutable DEMO_RNG across sweeps | Warning | Fixed — 02_simulator_demo.jl uses RHO_RNG/SPILL_RNG/SHIFT_RNG/MU_RNG |
| WR-04: SIM-02 test N=120 near tolerance | Warning | Fixed — test uses N=300 (comment WR-04 in test_simulator.jl) |
| WR-05: _w1_test no empty guard | Warning | Fixed — test_simulator.jl `_w1_test` has `n < 5 && return NaN` |
| WR-06: Translation axis convention reversed | Warning | Fixed — forward.jl line 172: `Translation(θ.shift_dy, θ.shift_dx)` |
| WR-07: Nuisance range validation missing | Warning | Fixed — forward.jl validates spillover∈[0,1], autofluorescence≥0, label_efficiency∈[0,1], noise≥0 |
| IN-01: Duplicate "(b)" comment label | Info | Fixed — demo.jl uses (b)/(c)/(d) labels correctly |
| IN-02: ghat.jl prerequisite undocumented | Info | Fixed — test_simulator.jl line 175 has the "IN-02: prior.jl transitively includes frozen ghat.jl" comment |
| IN-03: eps() denominator floor too small | Info | Fixed — forward.jl uses `max(std(f), 1e-10)` |

---

## Documented Caveats (Not Failures)

**D-03 near-Gaussian induced summary (known minor caveat):** The per-patch-correlation distribution
under the prior has excess kurtosis ≈ −0.05 — approximately Gaussian, not heavy-tailed. This does
not contradict the Turing ν ~ Exponential() prior (which permits near-Gaussian at large ν). ν is
"induced and checked, honestly near-Gaussian here, not set." Documented in NOTES.md §3 §"σ/τ/ν
consistency" and in ghat.jl evidence header. The SIM-02 gate is on μ (W1/KS), which passes.

**NOTES.md minor text inconsistency:** NOTES.md §3 line "SIM-02 prior consistency (W1 < tol at 512²,
N=120)" describes an earlier test plan; the committed test uses N=300 (more rigorous). The code and
gate are correct; the NOTES.md text is slightly stale.

**Negative μ-prior tail is PRIOR-ONLY (by design):** Real fluorescence images (positive and
negative controls) both yield positive induced μ (0.329 / 0.248). The simulator can generate
anti-correlated outputs (D-15 sign-flip verified), but no real data supports the negative tail.
Documented in NOTES.md §3 and ghat.jl evidence header per D-16. The SIM-02 consistency claim is
correctly scoped to the physically-realized μ range.

---

## Human Verification Required

None. All phase-2 criteria are quantitative and verified programmatically. The plausibility figure
aesthetics are supplementary to the passing SIM-04 quantitative gate.

---

## Gaps Summary

No gaps. All four SIM requirements satisfied, hard gate exits 0 (77/77), decoupling invariant
holds, no forbidden packages, no debt markers, all code-review issues resolved.

---

## Commit Verification

All SUMMARY-declared commits confirmed in `git log`:

| Commit | Plan | Description |
|--------|------|-------------|
| caf33e1 | 02-01 Task 1 | feat: extend spike env + resolve-risk/CUDA/pin guards |
| f68cc4d | 02-01 Task 2 RED | test: failing SIM-03 summary-contract testset |
| 329d93e | 02-01 Task 2 GREEN | feat: include() contract boundary, SIM-03 proven |
| caa8017 | 02-02 Task 1 | feat: implement simulate_pair 7-stage forward simulator |
| e11f8b2 | 02-02 Task 2 | test: SIM-01/SIM-03/D-15 on real simulator output |
| ef51011 | 02-03 Task 1 | feat: induced-μ calibration sweep + frozen monotone ĝ |
| edf3df1 | 02-03 Task 2 RED | test: failing SIM-02 prior-consistency testset |
| 0959b49 | 02-03 Task 2 GREEN | feat: sample_prior via frozen ĝ (SIM-02 GREEN) |
| 3e45819 | 02-03 Task 2 docs | docs: complete NOTES §3 — SIM-02 calibration evidence |
| 6eac6ad | 02-04 Task 1 | feat: seeded plausibility demo + CairoMakie figures |
| 55e330d | 02-04 Task 2 | test: SIM-04 quantitative plausibility gate |

---

_Verified: 2026-06-26_
_Verifier: Claude (gsd-verifier)_
