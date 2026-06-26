# Phase 1: Environment + Smoke Gate - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-26
**Phase:** 1-environment-smoke-gate
**Areas discussed:** Parent coupling, Smoke strength, Reproducibility scope, Gate durability

---

## Parent-package coupling

| Option | Description | Selected |
|--------|-------------|----------|
| Pkg.develop, include() fallback | Try `Pkg.develop(path="..")` first (clean boundary; Phase 4 ADVI needs the Turing model); fall back to `include()` of specific src files if GLMakie 0.10.5 + Turing + Flux 0.16 fail to co-resolve | ✓ |
| Pkg.develop only | Always Pkg.develop the full package; hard dependency on the whole tree resolving | |
| include() src files | Lightweight include() only; avoids dep-tree resolve but bypasses the package boundary | |

**User's choice:** Pkg.develop with include() fallback
**Notes:** Preserves the package boundary that Phase 4's ADVI baseline needs, while keeping Phase 1 unblocked if the parent dep tree (heavy GLMakie/Turing) won't co-resolve with Flux 0.16.

---

## Smoke strength

| Option | Description | Selected |
|--------|-------------|----------|
| Green + correctness assert | train + sampleposterior must run AND posterior recovers the known 1-param Gaussian mean within tolerance | ✓ |
| Runs green only | train + sampleposterior complete without error; no numerical check | |

**User's choice:** Green + correctness assertion
**Notes:** A stronger gate that catches a silently-wrong API/backend, not just a crash.

---

## Reproducibility scope

| Option | Description | Selected |
|--------|-------------|----------|
| Packages + Julia version | Commit spike/Manifest.toml AND record exact Julia version (.julia-version / juliaup channel) | ✓ |
| Packages only | Commit spike/Manifest.toml only; rely on ambient Julia | |

**User's choice:** Packages + Julia version
**Notes:** Main repo declares no Julia compat bound; NeuralEstimators/RxInfer are version-sensitive; reproducibility is a DoD item.

---

## Gate durability

| Option | Description | Selected |
|--------|-------------|----------|
| Re-runnable harness | Wrap smoke in spike/test/runtests.jl (stdlib Test) so the hard gate is re-checkable before each later phase; no external CI | ✓ |
| One-off script | spike/00_smoke.jl run manually; moment-in-time check | |

**User's choice:** Re-runnable harness
**Notes:** Makes "green smoke is the hard gate for all later phases" an actually-repeatable guard.

---

## Claude's Discretion

- Exact tolerance value for the correctness assertion.
- File/dir layout within `spike/` beyond the named files.
- Exact juliaup / `.julia-version` mechanism.
- `NOTES.md` vs a dedicated `STACK-DECISION.md` for the ENV-04 stack note.

## Deferred Ideas

- CUDA/GPU acceleration — only if Phase 4 CPU training is slow.
- External CI for the smoke gate — out of scope; local re-runnable harness instead.
- RxInfer baseline (BACK-01) — v2/deferred.
