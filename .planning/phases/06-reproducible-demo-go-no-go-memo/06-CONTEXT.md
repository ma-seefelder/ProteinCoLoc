# Phase 6: Reproducible Demo + Go/No-Go Memo - Context

**Gathered:** 2026-07-03
**Status:** Ready for planning

> **Provenance note.** The user launched `/gsd:discuss-phase 6` interactively but stepped
> away before answering. To avoid stalling the spike's closing phase, the decisions below
> were selected by Claude in the spirit of `--auto` — each is the *recommended* option,
> grounded in the Phase-5 artifacts and post-iteration STATE. They are **overridable**:
> re-run `/gsd:discuss-phase 6` (choose "Update it") or edit this file before planning.
> The three flagged as **[USER-OWNED]** — the memo verdict, the honesty framing, and the
> full-build-out call — are the ones most worth a human confirmation, because they are
> scientific/strategic judgments, not implementation mechanics.

<domain>
## Phase Boundary

Close the spike with a **falsifiable, reproducible verdict**. Three deliverables, no more:

1. **`spike/demo.jl`** — a single seeded script that chains every pipeline layer
   (prior → simulator → summary → NPE/NRE → SBC/BF/OOD) reproducibly from a fixed
   Random123 seed, CPU-only under the pinned `spike/Manifest.toml`, and tabulates the
   success criteria (DEMO-01, DEMO-02).
2. **Decoupling proof** — `git status` on `src/`, `src/bayes.jl`, `src/colocalization.jl`
   is demonstrably clean; the spike provably never edited the main package or either
   manuscript pipeline (DEMO-02, SC2).
3. **A 2–3-page Go/No-Go memo** reporting the metrics, framing SBC as "calibrated under
   the simulator" **paired with the OOD result**, and stating a **concrete full-build-out
   decision** (DEMO-03, SC3).

**Out of scope:** any `src/` edits or productionization (Phase 7 — starts only on a Go);
retraining/re-validating the nets (Phase 5 is closed — the memo *reports* its results, it
does not re-open them); the external corpus (Phase 8) and comparator harness (Phase 9,
already complete). **No `src/` edits** — the demo reaches `src/` only via read-only
`include()`, and the decoupling proof exists precisely to certify that.

Requirements: DEMO-01, DEMO-02, DEMO-03.
</domain>

<decisions>
## Implementation Decisions

### Memo verdict — the central call **[USER-OWNED]**
- **D-01:** The memo renders a **Conditional Go**, not a clean Go and not a No-Go.
  Rationale from the evidence as it actually stands after two bounded iterations
  (consts.jl byte-unchanged throughout):
  - The **core thesis is proven** — Phase 4 demonstrated amortized inference >100× faster
    than per-dataset ADVI at comparable point accuracy (ρ recovery corr 0.983).
  - **OOD: PASS** — after adding the noise-sensitive channel, pooled AUC 1.0, all 4
    families 1.0, ID fire-rate 0.05, negative controls behave.
  - **SBC: calibrated but the literal pre-registered gate rejects** — ECE is green on all
    8 parameters (ρ_true 0.016; was 0.19 red), but the M=2000 KS∧χ² conjunction still
    fails. Diagnosed cause: an over-powered χ² sub-criterion at M=2000 plus two genuinely
    non-uniform nuisance params (shift, label_efficiency) that need more training data.
  - **BF: near-miss, failure diagnosed as a baseline artifact** — corr 0.936 (bar 0.95);
    max|Δ logBF| large only because the clamped-KDE baseline hits its 1e-8 floor at the
    sweep tails (|Δρ|≳0.4), while the bounded NRE stays correct. Mid-range agrees within 1–2.
  - **The Conditional Go is gated on a bounded, freshly-pre-registered calibration
    close-out** (see D-07) that must clear before any `src/` integration ships in Phase 7.
  This verdict is honest: it does **not** claim the pre-registered gates passed as written;
  it reads the evidence as "method validated; two pre-registration design flaws + one
  data-scale gap; none are method-killers."
- **D-02:** The memo states the **falsification condition explicitly** — what result would
  have forced a No-Go (e.g. SBC ECE staying red after capacity+data iteration, or BF
  mid-range disagreement, or an OOD family undetectable by any summary-orthogonal channel)
  — so the Go is falsifiable, not rationalized.

### Pre-registration honesty in reporting **[USER-OWNED]**
- **D-03:** The memo reports **BOTH number sets, clearly labeled and in this order:**
  (1) the **original pre-registered gate run** at locked consts on the fresh
  `VAL_MASTER_SEED` (all three gates FAIL) as the *primary pre-registered result*; then
  (2) the **post-iteration confirmatory numbers** from the two bounded retrains, explicitly
  framed as *post-hoc*.
- **D-04:** The memo **names the data-snooping exposure plainly**: the frozen net was
  retrained twice *after* the pre-registered result was seen. It then states the three
  mitigations that keep the iterations defensible — `consts.jl` byte-unchanged (verified
  vs commit `e9c91d3`), model selection on a **disjoint DEV seed** (`0xDE7C0DE`), and a
  **single confirmatory VAL run** per iteration (no retry-to-pass). It does **not** present
  the post-iteration numbers as if they were the pre-registered result.
- **D-05:** Because the post-iteration numbers are post-hoc, the memo's forward gate (D-07)
  requires the close-out to be **re-pre-registered fresh** (new locked consts, new disjoint
  seed) so the Phase-7 ship-gate is clean of this spike's snooping exposure.

### demo.jl reproducibility scope
- **D-06:** **Two-tier demo.** Default (`julia --project=spike spike/demo.jl`) is **fast**:
  it loads the frozen artifacts (`spike/npe/trained_npe.jld2`,
  `spike/validation/trained_ratio.jld2`) and chains **every layer** end-to-end from the
  fixed seed at fixture scale (small M/L, short Δρ sweep, OOD subset) to *prove the pipeline
  reproduces*, tabulates the success-criteria table, and **loads the reported
  `*_report.jld2` artifacts for the headline numbers** (so the memo's figures are the
  reported-scale ones without a multi-minute recompute). A `--full` flag re-runs the
  reported-scale pipeline (M=2000 SBC, full BF sweep, full OOD grid) for a from-scratch
  faithful reproduction. Rationale: SC1 demands CPU-only reproducibility under the pinned
  Manifest; the fast default keeps `demo.jl` runnable in the "milliseconds-to-minutes"
  spirit of the whole thesis, while `--full` gives the auditable faithful path. This mirrors
  the fixture-vs-reported split already established in Phase 5
  (`test_*.jl` fixtures vs `run_*.jl` reported scripts).
  - The decoupling proof (SC2) runs inside `demo.jl` as an assertion: shell out to
    `git status --porcelain -- src/ src/bayes.jl src/colocalization.jl` and assert empty.

### Full-build-out decision **[USER-OWNED]**
- **D-07:** The memo commits to the **existing Phase 8–16 roadmap DAG as the full-build-out
  direction**, gated by a **calibration close-out** that must pass *before* Wave-B feature
  work begins: retrain at larger cache scale (toward 200k) to clear the residual SBC
  non-uniformity on shift/label_efficiency, re-express the BF gate against a non-clamped
  baseline, and **re-pre-register** SBC/BF thresholds for that run.
- **D-08:** Within the DAG, the memo **prioritizes the axes that most differentiate v2.0
  from Tapqir/Costes/Manders** (per the Phase-10 manuscript positioning): **registration-
  and-chromatic-uncertainty-as-latent (Phase 11)** and the **three-hypothesis amortized BF
  (Phase 13)** as the lead feature axes, with the **spatial coloc map (Phase 12)** following
  11. Hierarchy/3D/multi-channel are named as the *longer-horizon* build-out the memo's SC3
  asks for, sequenced after the calibration close-out and the Wave-B differentiators.

### Claude's Discretion
- Exact layout and section order of the memo (within the 2–3 page bound), figure selection
  from the reported `figures/` set, and the demo's success-criteria table formatting.
- Whether the decoupling proof also greps the two manuscript pipeline dirs by name or relies
  on the `src/`-clean assertion plus a whole-tree `git status` — planner/researcher's call.
</decisions>

<specifics>
## Specific Ideas

- The memo's SBC framing is a **carried-forward hard requirement** (Phase 5 D-02): SBC is
  reported "calibrated under the simulator" and **must be paired with the OOD result** in
  the same breath — SBC alone is never presented as a real-data guarantee.
- The verdict should read like a **pre-registration success story even though the literal
  gates failed**: the discipline (locked consts, disjoint seeds, honest negatives) is what
  lets the memo make a *credible* Conditional Go instead of a hand-wave.
</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & success criteria
- `.planning/ROADMAP.md` §"Phase 6: Reproducible Demo + Go/No-Go Memo" — SC1–3, DEMO-01/02/03.
- `.planning/ROADMAP.md` §"Phase 7" — the Go-conditional productionization contract the
  memo's verdict feeds; §"Progress"/"Execution Order" — the Phase 8–16 DAG the full-build-out
  decision (D-07/D-08) must reference.

### Phase-5 results the memo reports (do NOT re-open)
- `.planning/phases/05-validation-bundle-sbc-amortized-bf-ood/05-04-SUMMARY.md` — the
  **original pre-registered trifecta run** with verbatim numbers (all three gates FAIL);
  the primary pre-registered result for D-03.
- `.planning/STATE.md` §"Accumulated Context" → Decisions (the `[Phase 5-iter1]` and
  `[Phase 5-iter2]` entries) — the **post-iteration confirmatory numbers** (OOD→PASS,
  SBC ECE-green, BF corr 0.936 + clamped-baseline diagnosis) for D-03/D-04.
- `spike/validation/consts.jl` — the **byte-unchanged pre-registration** (locked M/L/bins,
  all SBC/BF/OOD thresholds, `VAL_MASTER_SEED=0x5BC0FFEE`); the honesty anchor for D-04.
- `.planning/phases/05-validation-bundle-sbc-amortized-bf-ood/05-CONTEXT.md` §decisions D-02 —
  the "SBC calibrated-under-the-simulator, paired with OOD" reporting discipline.

### Project-level
- `.planning/PROJECT.md` — core value statement + Key Decisions table (verdict must align
  with the stated >100× / SBC-proof / honest-OOD thesis).
- `spike/NOTES.md` — recorded NeuralEstimators-default / BayesFlow-fallback decisions (ENV-03/04).
</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **No `spike/demo.jl` exists yet** — it is the Phase-6 deliverable. The seeded-entry-script
  pattern is set by `spike/00_smoke.jl` and `spike/02_simulator_demo.jl` (AGPL header, guarded
  includes, CPU-only, headless figure save) — `demo.jl` follows that idiom at top level.
- **Reported gates to compose:** `spike/validation/run_sbc.jl`, `run_bf.jl`, `run_ood.jl` —
  each already owns its artifact + figures + a real hard gate (exits nonzero on fail). `demo.jl`
  drives these (or their fixture-scale equivalents) rather than re-implementing.
- **Shared harness:** `spike/validation/harness.jl` (`draw_simulate_infer`, paired Δρ path).
- **Frozen artifacts the fast tier loads:** `spike/npe/trained_npe.jld2`,
  `spike/validation/trained_ratio.jld2`, and the gitignored reported `*_report.jld2` +
  `spike/validation/figures/*.png`.

### Established Patterns
- **Reproducibility primitive:** `spike/data/seeding.jl` — `sample_rng(master_seed, idx)` =
  `Philox4x` keyed by `(seed, global_index)`, thread/order-independent. `demo.jl`'s fixed seed
  threads through this, not `Random.seed!`.
- **Pinned env:** `spike/Project.toml` + `spike/Manifest.toml` (NeuralEstimators 0.2.1, Flux
  0.16.10, no CUDA) — `demo.jl` runs `--project=spike`, CPU-only (`use_gpu=false` everywhere).
- **Fixture-vs-reported split:** Phase 5's `test_*.jl` (tiny M/L fixtures) vs `run_*.jl`
  (reported scale) is the exact precedent for D-06's two-tier demo.

### Integration Points
- **Read-only `src/` coupling:** `spike/contract.jl` reaches `src/` via `include()` only; the
  decoupling proof (SC2) asserts `git status` on `src/`, `src/bayes.jl`, `src/colocalization.jl`
  stays clean — the certification that this coupling never mutated.
</code_context>

<deferred>
## Deferred Ideas

- **The calibration close-out itself** (200k-scale retrain, non-clamped BF baseline,
  re-pre-registered thresholds) is *named* by the memo (D-07) but **executed in Phase 7**,
  not here — Phase 6 reports and decides, it does not retrain.
- **Re-enabling the OOD posterior-predictive channel** in the reported OR-fusion (the θ̂
  finite-guard from iter1 made it viable but `run_ood.jl` still runs `with_pp=false`) —
  a Phase-7 hardening item, not a Phase-6 deliverable.
- **RxInfer independent cross-check** of the ADVI posterior — already deferred post-Go
  (STATE Deferred Items, BACK-01); the memo may cite it as a "nice-to-have for the paper."

### Reviewed Todos (not folded)
None — `todo.match-phase 6` returned zero matches.
</deferred>

---

*Phase: 06-reproducible-demo-go-no-go-memo*
*Context gathered: 2026-07-03 (Claude-selected recommended defaults; user stepped away mid-discuss)*
