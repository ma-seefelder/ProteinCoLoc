# Phase 14: Decision and Abstention Layer - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in `14-CONTEXT.md` — this log preserves how they were reached.

**Date:** 2026-08-03
**Phase:** 14-decision-and-abstention-layer
**Mode:** discuss `--analyze --all` (all gray areas auto-selected; trade-off table before each question)
**Areas discussed:** Evidence lane, Conformal machinery, Calibration set, FDR unit, Abstention fusion,
Not-checked handling, Threshold provenance

---

## Scouting findings presented before questioning

Five findings were surfaced from disk before any question was asked. Three of them constrained what
could be built:

| # | Finding | Evidence | Effect |
|---|---|---|---|
| F1 | The three abstention inputs sit in different lanes | `src/amortized/ood.jl` (shipped) vs `spike/comparator/classical.jl` vs `spike/p13/` (DEV seed) | Framed the D-01 lane question |
| F2 | `ConformalPrediction.jl` exists nowhere in the repo | One grep hit, in the ROADMAP sentence itself | Made SC1's named library a live decision, not an inheritance |
| F3 | SC2 contradicts Phase 9's purpose | Phase 9 = "knows when the classics are wrong"; SC2 makes disagreement trigger silence | Forced the D-05 amendment |
| F4 | Phase 13 already owns τ | `spike/p13/tau_probe.jl` — "the whole three-way cut hangs off ONE number", measured not chosen | Framed D-07 |
| F5 | Two dependency phases closed NEGATIVE | `12-VERIFICATION.md` (no calibrated per-region uncertainty); `11-CLOSURE.md` (no reship) | Removed per-region FDR from scope (D-04) |

Two constraints were verified as **running tests**, not prose:
- `spike/test/test_p12_decoupling.jl:176-178` — `spike/Project.toml` / `Manifest.toml` byte-unchanged.
- `spike/test/test_p12_decoupling.jl:210` — the Phase-16 sealed holdout is not consumed.

Corpus counted directly: **32 data rows, 2 `sealed_holdout` (anchors-only), 30 open** → split-conformal
floor of α ≥ 1/31 ≈ 0.032.

---

## Questions and selections

### Area 1 — Evidence lane

**Options presented:** Bridge (spike/, src-shaped API) *[recommended]* / Spike lane, no bridge /
Ship in src/, 2-way only / Ship in src/, promote the 3-way net first

**Selected:** **Bridge — spike/, src-shaped API**

Recommendation reasoning given: Phase 13 chose the research lane specifically to keep the Phase-7 GO
untouched (its D-01), and nothing since changed that calculus — but a decision layer that cannot reach
`src/` makes Phase 16 harder. A src-shaped signature costs little now and preserves both options.
→ **D-01**

### Area 2 — Conformal machinery

**Options presented:** Hybrid (Bayesian FDR + hand-rolled conformal) *[recommended]* / Hand-rolled
conformal only / Bayesian FDR only / Use `ConformalPrediction.jl` as SC1 says

**Selected:** **Hybrid — Bayesian FDR primary + hand-rolled split conformal**

Recommendation reasoning given: the named library is the one option with a concrete verifiable cost
(breaks a running guard) and the least payoff — it wraps MLJ models, not a NeuralEstimators posterior,
for what amounts to a sorted-score quantile. Bayesian FDR is the honest primary because a calibrated
posterior is v2.0's central claim; conformal earns its place only as the misspecification hedge.
→ **D-02, and SC1 amended** (library name dropped; conformal met in substance)

### Area 3 — Conformal calibration set

**Options presented:** Simulator primary + corpus as reported check *[recommended]* / Re-partition the
corpus now / Corpus only (α ≥ 0.032) / Simulator only

**Selected:** **Simulator primary + corpus as reported check**

Tension surfaced: conformal's value is coverage that survives the model being wrong, so simulator
calibration is circular exactly where it matters — but the non-circular alternative consumes Phase 16's
blind evaluation rows. Reporting both numbers shows the sim-vs-real gap rather than hiding it.
→ **D-03**, with the circularity recorded as a named manuscript limit

### Area 4 — FDR unit

**Options presented:** Per image pair, tiles shown but not FDR-controlled *[recommended]* / Per pair
only, no tile output / Per experiment or condition group / Per tile or region

**Selected:** **Per image pair; tiles displayed but explicitly NOT FDR-controlled**

Evidence cited: `LocalColocMap` carries `delta_rho` and `ood_flag` but **no uncertainty field**.
Claiming per-region FDR after Phase 12 returned NO would repeat the closed-vs-succeeded conflation
corrected in ROADMAP at `e29c846`.
→ **D-04**

### Area 5 — Abstention fusion

**Options presented:** Asymmetric — disagreement abstains only with OOD *[recommended]* / Four-state
output / Literal OR as SC2 is written / Drop disagreement from abstention entirely

**Selected:** **Asymmetric fusion**

The F3 tension drove this: SC2's literal OR silences the tool exactly where Phase 9 says it should
speak, and three OR'd triggers compound toward an abstention rate that swallows the batch.
→ **D-05, and SC2 amended**

*Considered and not chosen:* the four-state output (`DECIDED-CONTRA-CLASSICAL`) — more transparent,
but a non-standard state needing extra defence. Recorded in Deferred as the fallback if reviewers
push back on the asymmetry.

### Area 5b — Not-checked OOD handling

**Options presented:** Abstain by default with explicit opt-out *[recommended]* / Hard error on missing
operating point / Decide but flag prominently

**Selected:** **Abstain by default, explicit opt-out**

Evidence cited: `local_map.jl:70-72` states outright that a `false` flag on a non-sentinel tile means
"NOT CHECKED, not in distribution"; `_recorded_ood_threshold` returns `nothing` when the gate report is
absent, unreadable, or non-finite. Reasoning given: the prevented failure is silent and scientific, the
imposed cost is loud and operational.
→ **D-06**

### Area 6 — Threshold provenance

**Options presented:** Inherit τ from Phase 13's artifact with provenance asserted *[recommended]* /
Inherit τ and pre-register α_FDR too / τ as user parameter with Phase 13 default / Re-derive τ in
Phase 14

**Selected:** **Inherit τ by loading the artifact, assert provenance**

Three numbers separated explicitly: τ (inherited, measured, never tuned here), α_FDR (user parameter
per SC1), α_conformal (pre-registered separately, not to be conflated with α_FDR).
→ **D-07**

---

## Scope creep redirected

None raised by the user. Four items were moved to Deferred by analysis rather than by request:
per-region FDR (blocked by Phase 12), corpus re-partitioning (pre-registration matter), promoting the
three-way net to `src/` (reopens the Phase-7 GO), and the four-state output (considered, not chosen).

## Claude's discretion recorded

Conformal score function; the concrete Bayesian-FDR estimator form; result type shape and field names;
whether the SC3 risk-coverage curve runs on simulator draws, corpus rows, or both.

## Amendments this discussion makes to ROADMAP Phase 14

| Criterion | Original | Amended by |
|---|---|---|
| SC1 | "using conformal sets (**ConformalPrediction.jl**) + decision-risk" | **D-02** — hand-rolled; library name dropped, conformal met in substance |
| SC2 | "OOD **∨** cross-method disagreement **∨** ambiguous conformal set" | **D-05** — asymmetric; disagreement alone reports rather than silences |

Both amendments must be cited alongside any Phase-14 result.
