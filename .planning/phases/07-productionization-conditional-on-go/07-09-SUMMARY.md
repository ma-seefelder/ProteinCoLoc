# 07-09 Summary — 32×32 CAP decision (training correctly skipped)

**Plan:** 07-09 (wave 7, `autonomous: false`) — conditional 32×32 pipeline + ship-vs-cap decision.
**Outcome:** **CAP the family.** 32×32 is NOT trained and NOT shipped. Recorded in
`gate-32x32.md`.

## What happened (deviation from the plan's train-then-decide flow — authorized)

The 07-09-PLAN was written before the calibration investigation and assumed 32×32 would be trained
and then gated, with the user choosing ship-with-caveat vs cap based on the gate. Reality changed
that ordering: the **Go/No-Go decision (2026-07-24, Option A — GO with named limits,
`07-GO-NO-GO-UPDATE.md`) determined the disposition up front.** The GO rests on the 8×8 grid via
post-hoc re-analysis (randomized-rank SBC) + simulation-based BF (spike 014); it does not depend on
any fine grid. The user's binding cap-the-family decision therefore made the conditional 32×32
training run unnecessary, and it was **correctly skipped**.

This is the plan's `cap-the-family` branch (Task 3), which the plan itself names as an expected
outcome rather than a failure. The plan is complete via that branch.

## Tasks

- **Task 1 (train 32×32 via `_train_grid_pipeline(32)`): SKIPPED.** No `julia -t auto` training run;
  no `artifacts/grid_32/*.jld2`. Rationale: a 32×32 gate would run the KDE-BF / non-atom-corrected
  apparatus that spikes 006–014 proved defective, so it would fail like 16×16 and contribute nothing
  — at large CPU cost.
- **Task 2 (fresh 32×32 pre-registration + gate): SKIPPED.** No `test/gate/gate_consts_32.jl`, no
  `PROD_SEED[32]` drawn, no `gate_report_32.jld2`. No pre-registered seed consumed.
- **Task 3 (user ship-vs-cap decision): DONE — CAP.** Recorded in `gate-32x32.md` with rationale.

## Deliverable

- `.planning/phases/07-productionization-conditional-on-go/gate-32x32.md` — CAP-DECISION note
  (3-point rationale: GO rests on 8×8; a 32×32 gate would run the defective apparatus; user
  cap-the-family decision). Cites `07-GO-NO-GO-UPDATE.md`.

## Consequences for 07-10

- `_SHIPPED_GRIDS = (8,)` — only 8×8 ships. 4/16/32 excluded, 64 dropped.
- 32 ∉ `_SHIPPED_GRIDS`; no `grid_32` entry in `Artifacts.toml`.

## Constraints honored

- No training. No gate run. No pre-registered seed consumed. `spike/` byte-untouched. Existing
  frozen artifacts (`amended_v2`, `grid_4/8/16`, `spike01x`) byte-identical. 64×64 stays dropped.
