# Phase 14: Decision and Abstention Layer - Context

**Gathered:** 2026-08-03 (discuss --analyze --all)
**Status:** Ready for planning

<domain>
## Phase Boundary

Turn calibrated posteriors + the three-way Bayes factor into an actionable batch decision
**{coloc / not / ABSTAIN}** at a controlled Bayesian FDR, abstaining exactly when the tool should be
silent.

**Five scouting findings reshape this phase. Read them before planning — three of them change what
can be built at all.**

1. **The three abstention inputs are not on the same footing.** SC2 wants
   `OOD ∨ cross-method disagreement ∨ ambiguous conformal set`. Verified on disk:
   **OOD is shipped and exported** (`src/amortized/ood.jl` — `ood_flag`, `ood_verdict`, `OODVerdict`,
   `is_ood`); **cross-method exists but only in the spike** (`spike/comparator/classical.jl` —
   `manders`, `costes_p`, `patch_correlation`); the **three-way BF is real and verified** (Phase 13,
   5.5/6 must-haves, goal achieved) but lives in `spike/p13/` on a **DEV seed in the research lane**.
   The *shipped* `bayes_factor` in `src/` is still **two-way**.

2. **`ConformalPrediction.jl` is named in SC1 and exists nowhere in the repo.** A single grep hit, in
   the ROADMAP sentence itself. It is in neither `Project.toml`. SC1's named mechanism is
   aspirational, not inherited — and adding it would break a running guard (see Constraints).

3. **SC2 as written contradicts Phase 9's entire purpose.** Phase 9 exists so v2.0 can be positioned
   as *"knows when the classics are wrong, not merely agrees with them."* SC2 makes disagreement with
   the classics a trigger for **silence** — so the tool would abstain precisely in the cases that
   demonstrate its headline differentiator. **Resolved by D-05, which amends SC2.**

4. **Phase 13 already owns one threshold.** `spike/p13/tau_probe.jl`: *"the whole three-way cut hangs
   off ONE number"*, and τ is **measured rather than chosen**. Phase 14 must inherit it, not
   re-derive it (D-07).

5. **Two dependency phases closed NEGATIVE, and that removes a success criterion's foundation.**
   Phase 12 delivered **no calibrated per-region uncertainty** (`12-VERIFICATION.md`: goal NOT
   achieved; the shipped deliverable is the neutralized ablation; pooled LRO coverage 0.9707 outside
   its pre-registered [0.87, 0.93]). Phase 11 never reshipped a registration-aware net. Consequently
   **a per-region FDR has nothing calibrated underneath it** — see D-04.

**Out of scope:** per-region FDR control (blocked on Phase 12's NO), promoting the three-way net into
`src/` (reopens the Phase-7 GO), corpus re-partitioning (a pre-registration matter), and the blind
external evaluation itself (Phase 16).

</domain>

<decisions>
## Implementation Decisions

### Evidence lane and integration surface

- **D-01: Bridge — build in `spike/`, against a `src/`-shaped API.** `decide_coloc(...)` is built in
  the spike research lane consuming Phase 13's three-way net (`spike/p13/three_way_net.jld2`,
  `spike/p13/net.jl`), but given the signature it would have in `src/` — taking a bundle + images and
  returning an `AbstractColocResult` subtype. This keeps SC1's three-way claim honest, keeps the
  Phase-7 GO untouched (matching Phase 13's D-01 research-lane posture), and makes later promotion
  mechanical rather than a rewrite.
  **Consequence for planning:** the decision layer is NOT shipped in v2.0. Say so plainly in any
  report; do not let "src-shaped" be read as "in src".

### Decision machinery

- **D-02: Hybrid — Bayesian FDR primary, hand-rolled split conformal as the hedge.**
  Bayesian FDR on the calibrated posterior is the primary rule in-distribution (this is the natural
  fit: a calibrated `P(H₀ | data)` gives direct expected-FDR control by sorting and thresholding, and
  SBC already backs the calibration). A **small hand-rolled split-conformal set** is the
  distribution-free hedge under misspecification. **`ConformalPrediction.jl` is NOT added.**
  **SC1 is AMENDED:** "conformal sets (ConformalPrediction.jl)" → conformal sets met *in substance*,
  hand-rolled; the library name is dropped. Rationale: the named library is the one option with a
  concrete verifiable cost (breaks the frozen-env guard) and the least payoff — split conformal is a
  sorted-score quantile, and the library wraps MLJ models rather than a NeuralEstimators posterior.

- **D-03: Conformal calibrates on held-out SIMULATOR draws; the corpus is a REPORTED CHECK only.**
  Fit conformal scores on held-out simulator draws (tight α, costs no corpus data), then **report
  observed coverage on the 30 open corpus rows beside it** as a named check.
  **This is deliberately not the strongest available claim.** Corpus-only calibration would be
  genuinely non-circular but caps α ≥ 1/31 ≈ **0.032** (n=30) and would **consume Phase 16's blind
  evaluation data**. Reporting both numbers shows the sim-vs-real gap instead of hiding it behind
  either source alone.
  **Named limit to carry into the manuscript:** the headline conformal guarantee is simulator-derived
  and therefore inherits the simulator's misspecification; the corpus check is what bounds that.

### Decision unit

- **D-04: FDR is controlled PER IMAGE PAIR. Tiles are displayed but explicitly NOT FDR-controlled.**
  Per-pair is the only unit where the calibration provably holds — the posterior, the three-way BF,
  and the SBC proof are all per-pair. The per-tile map from `local_coloc_map` is still surfaced, but
  **labelled in the result as uncontrolled display, not a controlled call.**
  **Rationale, stated so it is not re-litigated:** claiming per-region FDR after Phase 12 returned NO
  would be exactly the "counting a closed negative as a success" failure corrected in ROADMAP at
  `e29c846`. `LocalColocMap` carries `delta_rho` and `ood_flag` but **no uncertainty field** — there
  is nothing to control against.

### Abstention rule

- **D-05: ASYMMETRIC fusion. SC2 is AMENDED.**
  - `OOD` fires → **ABSTAIN**
  - ambiguous conformal set → **ABSTAIN**
  - cross-method disagreement **alone** → **DECIDE**, and record the disagreement as a named field on
    the result. Disagreeing with Costes/Manders is the product thesis, not a failure mode.
  - cross-method disagreement **∧** OOD → **ABSTAIN** (independent reason to distrust ourselves)

  **SC2's literal `∨` is replaced by this rule.** Rationale: the literal OR silences the tool exactly
  where Phase 9 says it should speak, and three OR'd triggers compound — a tool that is usually
  silent is not useful. This preserves SC2's protective intent without gutting Phase 9.
  *Considered and not chosen:* a four-state output adding `DECIDED-CONTRA-CLASSICAL` — see Deferred.

- **D-06: A `false` OOD flag that means NOT CHECKED must never read as "in distribution."**
  The OOD input is **three-valued**: fired / clear / not-checked. `_recorded_ood_threshold`
  (`src/amortized/local_map.jl:93`) returns `nothing` when the gate report is absent, unreadable, or
  records a non-finite threshold, and `local_map.jl:70-72` states outright that *"a `false` flag on a
  non-sentinel tile means NOT CHECKED, not in distribution."*
  **Not-checked ABSTAINS by default**, with an explicit opt-out keyword the caller must set
  deliberately. The reason for silence is recorded in the result so it is auditable.
  **Rationale:** the failure this prevents is silent and scientific; the cost it imposes is loud and
  operational.

### Threshold provenance

- **D-07: τ is INHERITED from Phase 13's artifact by loading it, with provenance asserted.**
  Load τ from the Phase-13 report rather than copying the literal into a Phase-14 constants file;
  **fail loudly** if the artifact is missing or its hash does not match. Copying the value is how two
  numbers silently diverge; re-deriving guarantees a second, different cut on the same hypothesis
  space.
  - **τ** (three-way cut) — inherited, measured by Phase 13, never tuned here.
  - **α_FDR** — a **user parameter**, per SC1's "user-set Bayesian FDR".
  - **α_conformal** — the miscoverage level for the hedge, **pre-registered separately**. Distinct
    from α_FDR and must not be conflated in code or prose.

### Claude's Discretion

- Exact conformal score function (posterior-probability-based vs BF-based) — pick during planning and
  justify against the calibration evidence.
- The concrete Bayesian-FDR estimator form (direct posterior-probability sorting vs a decision-risk
  formulation) — both satisfy D-02; pick the one that composes better with the three-way output.
- Result type shape and field names, provided it subtypes `AbstractColocResult` per the established
  hierarchy.
- Whether the risk-coverage curve (SC3) is computed on simulator draws, corpus rows, or both — but
  D-03's circularity caveat applies to whichever is chosen and must be stated.

</decisions>

<specifics>
## Specific Ideas

- The tool disagreeing with Costes/Manders is a **feature to display**, not an error to suppress —
  this is the Phase 9 positioning and D-05 exists to protect it.
- Abstention should be **explainable**: every ABSTAIN carries which trigger fired (OOD / ambiguous /
  not-checked / disagreement∧OOD). A silent abstention is nearly as bad as a wrong call.
- Follow the project's established discipline on numbers: **measured rather than chosen**,
  pre-registered before results exist, and **corrections appended rather than overwritten**.

</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Direct inputs (what this phase consumes)
- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-VERIFICATION.md` — Phase 13 goal
  achieved (5.5/6); confirms the three-way net exists and what it emits.
- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-CONTEXT.md` — D-01 research-lane
  posture, D-02 `zt` basis, D-08 structural-zero `random`, D-12/D-15 amendments.
- `.planning/phases/13-three-hypothesis-amortized-bayes-factor/13-REPORT.md` — what was compared and
  what was declined (note: §14 point 5 undercounts missing `status:` fields — nine files, not two).
- `spike/p13/net.jl`, `spike/p13/result.jl`, `spike/p13/three_way_net.jld2` — the evidence net,
  `ThreeHypothesisColocResult <: AbstractColocResult`.
- `spike/p13/tau_probe.jl` — τ provenance; the number D-07 inherits.
- `src/amortized/ood.jl` — `ood_flag`, `ood_verdict`, `OODVerdict`; the shipped OOD channel.
- `spike/comparator/classical.jl` — `manders`, `costes_p`, `patch_correlation`; the cross-method arm.

### Constraints that bind this phase (violating these breaks running tests)
- `spike/test/test_p12_decoupling.jl:176-178` — **`spike/Project.toml` and `spike/Manifest.toml` must
  be byte-unchanged vs HEAD.** This is why D-02 forbids adding `ConformalPrediction.jl`.
- `spike/test/test_p12_decoupling.jl:210` — **the Phase-16 sealed holdout must not be consumed.**
  Bounds D-03.
- `src/amortized/local_map.jl:70-72, :93` — the NOT-CHECKED semantics that D-06 exists to handle.
- `corpus/manifest.csv`, `corpus/anchor_rows.jl` — 32 data rows, 2 `sealed_holdout` (anchors-only,
  D-09), 30 open. The n=30 that caps α in D-03.

### Why two success criteria are amended, and why per-region is out of scope
- `.planning/phases/12-spatial-colocalization-map/12-VERIFICATION.md` — Phase 12 goal NOT achieved; no
  calibrated per-region uncertainty exists. Grounds D-04.
- `.planning/ROADMAP.md` Phase 12 entry + Closure block (as amended at `e29c846`) — the
  closed-vs-succeeded distinction this phase must not repeat.
- `.planning/phases/11-registration-and-chromatic-uncertainty-as-latent/11-CLOSURE.md` — Phase 11
  closed negative-but-useful; no registration-aware reship.
- `.planning/ROADMAP.md` Phase 9 entry — the "knows when the classics are wrong" positioning that
  D-05 protects.
- `.planning/ROADMAP.md` Phase 14 entry — SC1 and SC2 as originally written; **both are amended here
  (D-02, D-05) and the amendments must be cited alongside any Phase-14 result.**

### Project-level discipline
- `./.claude/skills/spike-findings-proteincoloc/SKILL.md` — seed discipline (DEV seeds asserted
  disjoint from `PROD_SEED_V2`, `VAL_MASTER_SEED`, `NPE_MASTER_SEED`, `CORPUS_MASTER_SEED`), and the
  rule that gate changes are pre-registration matters rather than in-phase decisions.
- `.planning/phases/07-productionization-conditional-on-go/07-GATE-AMENDMENT.md` — §6.4 blocks
  iterating the current amended gate.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `src/amortized/ood.jl` — `ood_flag` (already an OR-fusion of Mahalanobis + posterior-predictive
  channels), `ood_verdict`, `OODVerdict`, exported `is_ood`. Consume directly; do not reimplement.
- `spike/p13/net.jl` + `three_way_net.jld2` — `ThreeWayEvidenceNet <: NeuralEstimator`, shared trunk +
  `Dense(64,2)` head, emits `log BF(coloc:random)` and `log BF(exclusion:random)` in one pass.
- `spike/p13/result.jl` — `ThreeHypothesisColocResult <: AbstractColocResult` with an inner
  constructor enforcing `random === 0.0` (structural zero, D-08). The result-type pattern to follow.
- `spike/comparator/classical.jl` — `manders`, `costes_p`, `patch_correlation`; the cross-method
  signals D-05 consumes.
- `src/amortized/local_map.jl` — `LocalColocMap` (`grid`, `tiles`, `delta_rho`, `ood_flag`, `meta`).
  Note: **no uncertainty field** — the direct evidence for D-04.

### Established Patterns
- **`AbstractColocResult` subtype hierarchy** — every result type subtypes it; accessors
  (`delta_rho`, `bayes_factor`, `is_ood`, `posterior_draws`) are the public surface.
- **Research-lane / DEV-seed discipline** — spike work runs on seeds asserted disjoint from all
  pre-registered seeds; never consume a pre-registered seed, never write a gate report from one.
- **Pre-registered constants files** — `p12_consts.jl`, `spike/p13/consts.jl`; numbers are frozen
  before results exist and corrections are **appended, never overwritten**.
- **Executable constraint guards** — this project encodes its hard constraints as running tests
  (`test_p12_decoupling.jl`), not as prose. Phase 14 should do the same for D-03's corpus bound and
  D-06's not-checked rule.

### Integration Points
- `decide_coloc(...)` — new, in `spike/`, with a `src/`-shaped signature (D-01).
- Reads: the Phase-13 net + τ artifact, the shipped OOD channel, the comparator outputs.
- Returns: an `AbstractColocResult` subtype carrying the call, the abstention trigger (when silent),
  the cross-method disagreement field (D-05), and the uncontrolled per-tile display (D-04).

</code_context>

<deferred>
## Deferred Ideas

- **Per-region / per-tile FDR control** — blocked on Phase 12's NO (no calibrated per-region
  uncertainty). Revisit in v2.1 alongside SPAT-07, which was itself deferred there.
- **Re-partitioning the corpus for a non-circular conformal guarantee** — genuinely better science
  (removes D-03's circularity caveat) but shrinks Phase 16's blind evaluation set and needs a
  pre-registered split ratio. A pre-registration matter, not a Phase-14 decision.
- **Promoting the three-way net into `src/`** — would make the decision layer shippable, but reopens
  the Phase-7 GO and pulls a DEV-seed research artifact onto the production path. Deliberately not
  done here (D-01).
- **Four-state output adding `DECIDED-CONTRA-CLASSICAL`** — considered during discussion as an
  alternative to D-05's asymmetric rule. Rejected for v2.0 as a non-standard output state needing
  extra defence, but it is the more transparent design if reviewers push back on the asymmetry.
- **Adding `ConformalPrediction.jl`** — would require amending or exempting the frozen-env guard and
  accepting the MLJ dependency tree. Revisit only if the hand-rolled conformal proves insufficient.

</deferred>

---

*Phase: 14-decision-and-abstention-layer*
*Context gathered: 2026-08-03*
