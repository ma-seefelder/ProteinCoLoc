# Phase 6 — Deferred / Out-of-Scope Discoveries

Items discovered during Phase-6 execution that are NOT caused by this phase's changes and
are out of scope per the executor SCOPE BOUNDARY. Logged, not fixed.

## D6-DEFER-01 — Flaky aggregate test: test_ood.jl SC5 block-permute negative control

**Discovered during:** Plan 06-01 Task 3 verification (`julia --project=spike spike/test/runtests.jl`).

**Symptom:** `test_ood.jl` → "SC5 (D-04) summary-orthogonal negative control — the named
blind spot" intermittently reports 1 failure (block-permute OR-flag expected `false`) when
run inside the aggregate `runtests.jl`, but passes 7/7 when run standalone. Re-running the
aggregate suite alternates pass/fail.

**Root cause (pre-existing, NOT demo.jl):** `runtests.jl` does not load `spike/demo.jl`
(verified — no reference), so demo.jl cannot affect it. The OOD fixture's OR-flag path
(`_ood_fix_flag` → `ood_flag` → `pp_mismatch_score`) calls `posterior_for`/`sampleposterior`,
which draws from Julia's **global default RNG**. The aggregate suite never calls
`Random.seed!`, and Julia seeds `Random.default_rng()` non-deterministically at each process
start, so the PP channel's threshold/score — and thus the block-permute OR-flag on a
borderline base image — varies run-to-run. The test_ood.jl source itself flags this as an
honest residual: "a strongly OOD-leaning ID base plus block-permute could tip over — the
honest residual, not a hidden failure."

**Why out of scope:** Not caused by Phase-6 changes (demo.jl is standalone, never included
by runtests.jl). The flakiness lives entirely in the pre-existing `test_ood.jl` fixture.

**Suggested Phase-7 fix (do not apply in Phase 6):** seed the global RNG deterministically
at the top of `test_ood.jl`'s `_build_ood_fixture` (e.g. `Random.seed!(VAL_FIX_SEED)`) so the
PP channel is reproducible in the aggregate suite, matching the discipline `spike/demo.jl`
adopts for its NPE twin-run proof. Note: `spike/demo.jl` runs the OOD fixture with
`with_pp = false`, so it is unaffected by this global-RNG dependence and its OOD twin-run
proof is deterministic.
