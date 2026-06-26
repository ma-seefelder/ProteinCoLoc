---
phase: 01-environment-smoke-gate
verified: 2026-06-26T17:45:00Z
status: passed
score: 8/8 must-haves verified
has_blocking_gaps: false
re_verification: false
---

# Phase 1: Environment + Smoke Gate Verification Report

**Phase Goal:** A reproducible, isolated spike environment exists and the pre-1.0 NeuralEstimators + Flux stack is proven to work on this machine — the highest-risk unknown is retired before any other investment.
**Verified:** 2026-06-26
**Status:** PASSED
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Req | Status | Evidence |
|---|-------|-----|--------|----------|
| T1 | Baseline ref `f581d95` recorded before spike work; root files byte-identical to it | ENV-01 | VERIFIED | `01-BASELINE.md` records `commit-as-baseline` decision, hash `f581d95dbc3596c5d9a064d6325ed70ca43a79c0`; `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` exits 0 |
| T2 | Parent package accessible read-only; include() fallback documented with verbatim conflict evidence | ENV-01 | VERIFIED | `spike/NOTES.md` §4 records the Pkg.develop downgrade (0.2.1→0.1.4) + MethodError evidence + include() fallback decision; pre-authorized by D-01 in Plan 01-04 must_haves (explicit OR condition) |
| T3 | Isolated spike env has exactly 3 deps (NeuralEstimators, Flux, Distributions); no top-level CUDA installed | ENV-02 | VERIFIED | `spike/Project.toml` contains exactly those 3 deps; `grep -c "^\[\[deps\.CUDA\]\]" spike/Manifest.toml` returns 0 |
| T4 | `spike/00_smoke.jl` trains PosteriorEstimator (NormalisingFlow) CPU-only and passes correctness assertion | ENV-02 | VERIFIED | Smoke run: `mu_hat=0.79911894` vs `theta_true=0.7`, `|delta|=0.099 < tol=0.3`; `use_gpu = false` on line 67 |
| T5 | `spike/test/runtests.jl` wraps smoke in @testset, asserts correctness + CUDA absence; exits 0 | ENV-02 | VERIFIED | `julia --project=spike spike/test/runtests.jl` exits 0; 4/4 tests pass in 3m32.3s |
| T6 | `spike/Manifest.toml` committed with exact NeuralEstimators 0.2.1, Flux 0.16.10, julia_version 1.12.6; no top-level CUDA | ENV-03 | VERIFIED | `git ls-files --error-unmatch spike/Manifest.toml` succeeds; versions confirmed via grep |
| T7 | Julia 1.12.6 enforced via juliaup directory override + documented in `spike/.julia-version` | ENV-03 | VERIFIED | `juliaup override status` shows `ProteinCoLoc\spike → 1.12.6`; `cat spike/.julia-version` returns `1.12.6` |
| T8 | `spike/NOTES.md` records NeuralEstimators.jl as default SBI engine and BayesFlow (PythonCall) as fallback-only, with excluded packages named | ENV-04 | VERIFIED | §1 states "Default SBI engine: NeuralEstimators.jl (v0.2.1)" and "Fallback only: BayesFlow (Python, via PythonCall/CondaPkg)"; excluded packages table present |

**Score:** 8/8 truths verified

**Note on T2 (SC-1 wording vs actual):** The ROADMAP SC-1 text says "available read-only via `Pkg.develop`". The actual coupling is via the pre-authorized D-01 include() fallback, not Pkg.develop. This is NOT a failure: Plan 01-04's `must_haves.truths` explicitly states the OR condition — "via `Pkg.develop(path="..")`, **OR** the documented include() fallback is in place when the parent dep tree does not co-resolve with Flux 0.16" — and the fallback IS in place with verbatim conflict evidence. The PLAN is the authoritative truth source and the fallback is a complete, documented outcome.

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `.planning/phases/01-environment-smoke-gate/01-BASELINE.md` | Baseline ref + verify command | VERIFIED | Contains `f581d95`, reconciliation decision, and `git diff --quiet f581d95 -- ...` command |
| `spike/Project.toml` | Minimal isolated deps: NeuralEstimators, Flux, Distributions | VERIFIED | Exactly 3 deps, no CUDA, no heavy root deps |
| `spike/Manifest.toml` | Committed reproducibility artifact, exact versions, no top-level CUDA | VERIFIED | git-tracked; NeuralEstimators 0.2.1, Flux 0.16.10, julia_version 1.12.6; `[[deps.CUDA]]` count = 0 |
| `spike/00_smoke.jl` | CPU-only NPE smoke with `use_gpu=false` + correctness assertion | VERIFIED | 82 lines; AGPL header; `use_gpu = false` line 67; `@assert abs(mu_hat - theta_true) < tol` line 78 |
| `spike/test/runtests.jl` | @testset wrapping smoke + CUDA absence guard | VERIFIED | `@testset "NeuralEstimators CPU smoke" verbose = true`; includes 00_smoke.jl; 3 CUDA-absence @tests |
| `spike/.julia-version` | Contains literal `1.12.6` | VERIFIED | Content: `1.12.6` |
| `spike/NOTES.md` | Stack decision + Julia pin + coupling outcome + decoupling proof | VERIFIED | §1 (stack), §2 (Julia pin), §4 (coupling), §5 (decoupling proof); all git-tracked |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `spike/test/runtests.jl` | `spike/00_smoke.jl` | `include(joinpath(@__DIR__, "..", "00_smoke.jl"))` | WIRED | Line 46 — include confirmed present |
| `spike/00_smoke.jl` | `NeuralEstimators.train` | `use_gpu = false` at the train call | WIRED | Line 67: `train(est, prior_sampler, gaussian_simulator; K = 8000, epochs = 60, use_gpu = false, verbose = false)` |
| `spike/NOTES.md` | `spike/Manifest.toml` + juliaup override | Records `juliaup`, `BayesFlow`, `fallback` | WIRED | §2 records juliaup mechanism and package pin; §1 records fallback |
| `spike/NOTES.md` | Parent package (`src/`) | `include()` fallback documented | WIRED | §4 states `include("../src/colocalization.jl")` pattern for Phase 2+ |
| `01-BASELINE.md` | git ref `f581d95` | Recorded hash + verify command | WIRED | Hash `f581d95dbc3596c5d9a064d6325ed70ca43a79c0` present; `git diff --quiet f581d95 -- ...` command recorded |

---

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| ENV-02: NPE smoke trains + recovers parameter | `julia --project=spike spike/test/runtests.jl` | `4/4 pass, mu_hat=0.79911894, |delta|=0.099 < tol=0.3` in 3m32.3s | PASS |
| ENV-01: Root decoupling proof | `git diff --quiet f581d95 -- Project.toml Manifest.toml src/` | exit 0, `root untouched` | PASS |
| ENV-03: Exact versions in Manifest | `grep version spike/Manifest.toml` for NeuralEstimators + Flux | `0.2.1` and `0.16.10` | PASS |
| ENV-03: Julia pin enforced | `juliaup override status` | `ProteinCoLoc\spike → 1.12.6` | PASS |
| ENV-02: No top-level CUDA | `grep -c "^\[\[deps\.CUDA\]\]" spike/Manifest.toml` | `0` | PASS |
| ENV-01: ProteinCoLoc not in Manifest (develop reverted) | `grep -i ProteinCoLoc spike/Manifest.toml` | no match | PASS |

---

### Requirements Coverage

| Requirement | Source Plans | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| ENV-01 | Plans 01-01, 01-04 | Isolated spike/, main package accessible read-only, root untouched | SATISFIED | Baseline recorded; include() fallback documented; `git diff f581d95` exits 0; `src/` clean |
| ENV-02 | Plan 01-02 | <30-line smoke trains PosteriorEstimator, runs sampleposterior, green | SATISFIED | 4/4 tests pass, exit 0, mu_hat within tolerance; CUDA absent |
| ENV-03 | Plan 01-03 | Manifest pinned+committed; smoke is hard gate | SATISFIED | NeuralEstimators 0.2.1, Flux 0.16.10, julia_version 1.12.6 committed; juliaup override 1.12.6 |
| ENV-04 | Plan 01-03 | Stack-decision note: NeuralEstimators default, BayesFlow fallback-only | SATISFIED | NOTES.md §1 with default/fallback statement + excluded packages table |

**Orphaned requirements:** None. All four phase requirements (ENV-01–04) are claimed in plans and evidenced in code.

---

### Anti-Patterns Found

| File | Pattern | Severity | Assessment |
|------|---------|----------|------------|
| (none) | — | — | No TBD/FIXME/XXX/HACK/PLACEHOLDER markers found in any phase-modified file |

---

### Human Verification Required

None. All phase-1 acceptance criteria are verifiable programmatically (file presence + content inspection + `git diff` + test run). The smoke test was executed in this session and produced a deterministic, seeded result.

---

### Gaps Summary

No gaps. All 8 must-have truths verified, all artifacts substantive and wired, all key links confirmed, ENV-01 through ENV-04 satisfied, hard gate (smoke) passed with exit 0.

---

## Commit Audit

All commits referenced in summaries verified present in git history:

| Commit | Summary reference | Verified |
|--------|------------------|---------|
| `f581d95` | Baseline freeze | EXISTS |
| `56f85c3` | Minimal spike env | EXISTS |
| `6124c62` | CPU-only NPE smoke + Test gate | EXISTS |
| `9a20327` | Julia 1.12.6 pin + .julia-version | EXISTS |
| `0283a46` | NOTES.md stack decision | EXISTS |
| `c2717d9` | Coupling outcome (include fallback) | EXISTS |
| `b4ae100` | Decoupling proof | EXISTS |

---

_Verified: 2026-06-26T17:45:00Z_
_Verifier: Claude (gsd-verifier)_
