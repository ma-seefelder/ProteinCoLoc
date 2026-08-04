# Phase 15: Calibration Operating Envelope and CI Gate - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in `15-CONTEXT.md` — this log preserves how they were reached.

**Date:** 2026-08-04
**Phase:** 15-calibration-operating-envelope-and-ci-gate
**Mode:** discuss `--analyze --all` (all gray areas auto-selected; trade-off table before each
question group)
**Areas discussed:** Sweep axes, Break criterion, OOD-before-break evidence, CI substrate,
Net + pool lane, Pre-registration and failure discipline

---

## Area 1 — Sweep axes

**Trade-off presented:** SC1's literal four (spillover/PSF/autofluorescence/registration) vs the four
shipped `OOD_FAMILIES` (texture/noise/optics/background) vs a mechanism-typed union.
Recommendation: the union, on the grounds that the in-prior/out-of-model typing is the finding.

| Question | Options presented | Selected |
|---|---|---|
| Which nuisance axes does the sweep cover? | Union typed by mechanism (Rec.) / SC1's literal four / shipped four only / SC1 four + noise | **Union typed by mechanism** → D-02 |
| How far do in-prior axes push? | Within AND beyond prior (Rec.) / beyond only / within only | **Within AND beyond**, prior boundary a marked rung → D-03 |

**Scouting inputs that shaped the options:** `misspec.jl:181-183` (spillover and registration have no
generator); `simulator.jl:92,138-142,179` (three of SC1's four axes are inferred nuisances inside the
prior; only `σ_psf` is fixed).

---

## Area 2 — Break criterion

**Trade-off presented:** KS/χ² p-value vs coverage band vs ECE vs ECE-primary-with-both-reported.
Recommendation: ECE primary, anchored to the shipped net's own in-prior ECE, with band and p-value
reported but non-gating.

| Question | Options presented | Selected |
|---|---|---|
| What statistic declares a break? | ECE primary + band/p reported (Rec.) / band only / KS-χ² p / any-of-three union | **ECE primary, band + p non-gating** → D-04 |
| Which target defines the envelope? | Δρ and ρ_true only (Rec.) / all eight equally / Δρ only | **Δρ and ρ_true only** → D-05 |

**Prior-decision inputs:** the Phase-7 over-powered-χ² diagnosis at M = 2000; the Phase-12 ruling
that its coverage band tracked training duration rather than calibration; the ~0.08 SD nuisance
drift already disclosed as a named v2.0 limit.

---

## Area 3 — OOD-before-break (SC2)

**Trade-off presented:** per-axis strict crossing vs pooled margin vs per-axis three-valued vs
report-only. Recommendation: per-axis three-valued with all states pre-declared, because "the flag
never fired" has two opposite meanings.

| Question | Options presented | Selected |
|---|---|---|
| How is SC2 declared? | Per-axis three-valued (Rec.) / per-axis strict / pooled margin / report-only | **Per-axis three-valued**, PROTECTIVE / LATE / SILENT-BUT-SAFE → D-06 |
| What counts as "firing"? | Fire-rate above ID baseline (Rec.) / fixed rate ≥50% / any single image / AUC-based | **Fire-rate above measured ID baseline** → D-07 |

**Constraint raised during the area:** `OOD_ID_QUANTILE = 0.95` means ~5% of in-distribution images
fire by construction, so an any-image rule would pass SC2 trivially on every axis.

---

## Area 4 — CI substrate (SC3)

**Trade-off presented:** golden-value regression vs derived tolerance band vs two-tier vs
assertion-only-in-runtests. Recommendation: two-tier with a golden fast tier, noting that
counter-based Philox seeding makes a deterministic golden viable here where it usually would not be.

| Question | Options presented | Selected |
|---|---|---|
| What is the CI substrate? | GH Actions two-tier (Rec.) / GH Actions fast only / local hook + runtests / two-tier on Windows + Linux | **GitHub Actions, two-tier** → D-08 |
| How is drift detected? | Golden-value tight tolerance (Rec.) / derived tolerance band / artifact hash only / golden + hash | **Golden-value, tight tolerance** → D-09 |

**Scouting input:** no `.github/`, no active hooks, no Makefile exist — SC3's CI had to be built, not
configured. `SBC_M = 2000 × SBC_L = 999` over a 512²–2048² mixture does not fit any push-triggered
job; `SBC_FIX_M = 8` already exists as the right-shaped fixture.

Windows-matrix CI was offered as an option and not selected; carried to Deferred rather than dropped.

---

## Area 5 — Net and draw provenance

**Constraint surfaced before the questions:** `run_gate.jl:47-52` binds `gate_consts_8_v2.jl` to
exactly one run on `PROD_SEED_V2[8]`, and that run is spent (it produced
`artifacts/amended_v2/grid_8/gate_report_8.jld2`). Phase 15 cannot reuse that seed.

**Trade-off presented:** shipped net only vs all three grids vs shipped + one contrast.
Recommendation: shipped `grid_8` reported, `grid_16` as one secondary contrast, because higher
summary dimension is the informative direction.

| Question | Options presented | Selected |
|---|---|---|
| Which net(s)? | grid_8 + grid_16 contrast (Rec.) / grid_8 only / all three grids | **grid_8 reported + grid_16 contrast** → D-10 |
| How are draws generated (C-04)? | Fresh Phase-15 gate seed with asserted disjointness (Rec.) / carve from pool with dual index assertion / validation stream | **Fresh `PROD_SEED_P15`, disjointness asserted** → D-11 |

**Rejected-with-reason:** pool carving inherits val-block model-selection contamination while passing
every index check; the `VAL_MASTER_SEED` stream drove early stopping for the shipped net.

---

## Area 6 — Pre-registration and failure discipline

**Trade-off presented:** full two-tier `p15_consts.jl` vs bars-only vs extending
`gate_consts_8_v2.jl`. Recommendation: a full new two-tier file, given the four prior
"the bar was wrong" amendments and the stated licensing standard for each.

| Question | Options presented | Selected |
|---|---|---|
| How is Phase 15 pre-registered? | Full two-tier `p15_consts.jl` (Rec.) / bars only / extend `gate_consts_8_v2.jl` | **Full two-tier `p15_consts.jl`** → D-12 |
| What is authorized on a degenerate envelope? | Pre-declare both readings + frozen endpoints + 1 iteration (Rec.) / same with zero iterations / unlimited ladder extension | **Both readings pre-declared, endpoints frozen, `P15_ITERATION_ALLOWANCE = 1`** (ladder extension only, on "never breaks") → D-13 |
| What happens if an axis is LATE? | Report as named limit and close (Rec.) / retune OOD threshold / add a detector channel | **Report as named limit and close** → D-14 |

Both non-selected responses to the LATE question were carried to Deferred with their reasons, so the
refusal is on the record rather than implicit.

---

## Deferred Ideas Captured

- New OOD detector channel for a LATE axis (v2.1; forbidden as a Phase-15 response)
- Retuning `OOD_ID_QUANTILE`
- `grid_4` sweep / full three-grid envelope trend
- Per-region envelope mapping (blocked on Phase 12)
- `windows-latest` in the fast-tier CI matrix
- Scheduled (nightly/weekly) slow tier instead of release-triggered

## Claude's Discretion Items

Generator functional forms, ladder resolution and spacing, the concrete golden statistic, report
schema/naming, and whether the `grid_16` contrast sweeps all six axes or a stated subset.

## Scope Creep

None raised — the discussion stayed inside the ROADMAP boundary. Adjacent capabilities that surfaced
(new detector channels, per-region mapping) were captured as Deferred rather than folded in.
