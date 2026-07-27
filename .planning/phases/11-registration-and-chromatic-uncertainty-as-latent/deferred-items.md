# Phase 11 — deferred items

Out-of-scope discoveries logged during execution. Not fixed here.

## From plan 11-02 (chromatic ε as an 8th θ column)

- **`scripts/diag_bounded_theta_dev.jl:62-65` still hardcodes a θ arity of 7**
  (`Matrix{Int}(undef, M, 7)`, `Matrix{Float64}(undef, M, 7)`, `zeros(Float64, 7)`,
  `Matrix{Float64}(undef, M, 7)`).
  It is a developer diagnostic script, not part of the package or of `Pkg.test()`, so nothing
  currently breaks — but it would under-count columns if re-run against the extended prior.
  Not fixed here: `scripts/` is outside plan 11-02's file scope, and touching it would widen the
  D-11/D-12 provenance commit for no correctness gain.

- **`test/gate/sbc.jl` `SBC_PARAM_LABELS` is still 8 labels (7 θ + Δρ)** and
  `sbc_ranks_and_spread` still builds an `M×8` rank table over `for p in 1:7`.
  This is deliberate and matches `11-RESEARCH.md` §A4 row 25: the src-side ship gate is not
  re-run in Phase 11 (D-01), the frozen pre-registration must not be edited, and the loop is
  arity-safe against a wider posterior (it reads the first 7 rows of an 8-row draw matrix).
  Revisit only in a phase that actually re-runs the src gate.

## From plan 11-04 (research-net scaffold)

- **`spike/test/test_simulator.jl:199-200` still pins the PRE-11-02 seven-field θ**, so
  `julia --project=spike spike/test/runtests.jl` currently exits 1:

  ```
  sample_prior: 7-field θ, in-range ρ_true, deterministic (D-14): Test Failed at
    spike/test/test_simulator.jl:200
    Expression: keys(θ) == (:ρ_true, :spillover, :autofluorescence, :label_efficiency,
                            :shift_dx, :shift_dy, :noise)
    Evaluated:  (:ρ_true, …, :noise, :chromatic_eps) == (:ρ_true, …, :noise)
  ```

  **Pre-existing, and not caused by plan 11-04.** `spike/simulator/prior.jl` gained
  `chromatic_eps` in `ca02b0e` (plan 11-02); `spike/test/test_simulator.jl` was last touched in
  `02c3bee` (Phase 2) and was never updated to the post-edit arity the way its sibling
  `spike/test/test_p11_forced_theta.jl` was (11-02 deviation 3). The failure is present at plan
  11-04's worktree base `2aae21f` with no 11-04 file loaded.

  **Blast radius is larger than one assertion:** `test_simulator.jl` THROWS, so
  `runtests.jl` aborts there and the five later includes (`test_data_pipeline.jl`,
  `test_npe.jl`, `test_sbc.jl`, `test_bf.jl`, `test_ood.jl`, `test_comparator.jl`) do not run
  at all.

  Not fixed here: `spike/test/test_simulator.jl` is outside plan 11-04's declared file set, and
  plan 11-03 was executing in parallel against the same phase, so an undeclared edit risked a
  merge collision on a file neither plan owns. The fix itself is mechanical and mirrors what
  11-02 already did for `test_p11_forced_theta.jl` — extend the expected tuple with
  `:chromatic_eps`, keep the assertion a tripwire rather than a restatement, and update the
  testset name, which still says "7-field θ".

  The Phase-11 clause (i) dependency gate that plan 11-04 added to `runtests.jl` is
  **unaffected and green** (`D-04 CPU-only (no CUDA dependency, none loaded)`: 23/23 pass,
  up from 21 assertions).

  **RESOLVED in `b298c00`** (orchestrator post-merge integration gate, wave 3). The fix was
  larger than the one assertion above: five further sites hardcoded the arity and each was
  masked by the previous abort. See that commit for the full list. All arity sites are now
  derived from `sample_prior`. Phase 3 is 70/70; Phase 5 (SBC/BF/OOD) and Phase 9 (comparator)
  are green. One gate remains red — the next item.

## From the wave-3 post-merge integration gate (orchestrator, commit `b298c00`)

- **BLOCKER — `SPEEDUP_GATE` (NPE-03) now FAILS: median speedup 92.5× and 84.0× on two
  consecutive runs, against a gate of `> 100.0×`** (`spike/test/test_npe.jl:230`,
  gate defined at `:74`).

  This is the **only** remaining failure in `spike/test/runtests.jl`. It was invisible before
  `b298c00` because `test_simulator.jl` aborted the suite four includes earlier.

  **Reproducible, not noise:** 92.50× then 83.97× on two separate serial runs with no other
  agents active. Both are well below the gate.

  **Most likely cause, stated as a hypothesis and NOT yet confirmed by measurement:** the spike
  NPE now trains an **8-marginal** flow instead of a 7-marginal one, because `train_fold`
  derives the flow width from the (now 8-column) θ. A wider flow is a slower forward pass, and
  the gate measures `t_advi / t_npe`. Confirming this needs an A/B of the same benchmark at
  `D_flow = 7` vs `8`, which has not been run.

  **Scope — this does NOT invalidate the shipped >100× claim.** The shipped `amended_v2/grid_8`
  bundle still carries a **7-marginal** flow (`NPE_D = 7`, unmoved per D-16), and its
  `Artifacts.toml` pin is byte-unchanged. The failing number comes from a spike-side net that,
  post-D-09, is no longer the shipped one. The claim in `CLAUDE.md` (">100× faster than the
  existing per-dataset ADVI") is about the shipped path and is untouched by this.

  **Deliberately NOT fixed, and `SPEEDUP_GATE` deliberately NOT lowered.** Lowering a
  pre-registered threshold after seeing the result it failed is exactly the "amended twice,
  credibility spent" pattern Phase 7 was burned by (see `07-GATE-AMENDMENT.md` §6.4 and user
  memory `phase7-amended-gate-diagnosis`). The honest options are (a) confirm the 7-vs-8 flow
  hypothesis and record it as a named limit scoped to the research net, or (b) accept that the
  legacy Phase-4 speedup gate is not meaningful against an 8-column prior and retire it
  explicitly rather than relax it. Both are decisions for the user, not for execution.

  **Carry into the Phase-11 report (plan 11-11).** The report must state that adding
  `chromatic_eps` cost measurable inference speed on the spike net, and must be scrupulous
  that the shipped bundle's speed claim is unaffected — the same "which model did this number
  come from" discipline D-16 already requires.

---

## Structural: `git worktree remove --force` destroys gitignored bulk artifacts

**Found:** 2026-07-27, between plan 11-07 and the SC1g diagnosis.

**What happened.** Plan 11-07's executor generated the 50,000-pair training pool
(~54 MB, 56 min of CPU) into `spike/data/cache/p11/`, which `.gitignore:407`
(`spike/data/cache/*`) excludes. The executor ran in a git worktree. Standard wave cleanup runs
`git worktree remove --force`, which deleted the entire worktree directory **including the
untracked pool**. The pre-removal rescue step in `execute-phase.md` only rescues `*SUMMARY.md`,
so nothing caught it. The loss was silent: 11-07-SUMMARY.md's "available for reuse, at no further
compute cost" line was true when written and false an hour later.

**Why it recurs.** The hazard is structural, not a mistake anyone made. Any expensive artifact that
is (a) deliberately gitignored because it is bulk and regenerable, and (b) produced inside a
worktree, is destroyed by cleanup by construction. Phase 13 has the same shape.

**Proposed fixes, in order of preference.**
1. **Write expensive caches outside the worktree.** Resolve cache roots against the primary
   worktree (`git worktree list --porcelain | head -1`) rather than the cwd, so a cache written by
   a worktree agent lands in the main checkout and survives cleanup. Fixes it at the source and
   needs no workflow change.
2. **Extend the rescue step beyond `*SUMMARY.md`** to a declared artifact list — e.g. a plan
   frontmatter key naming paths to preserve before removal.
3. **At minimum, warn.** Have cleanup refuse (or loudly report) when a worktree contains untracked
   files above a size threshold, instead of deleting them silently.

**Mitigating property, specific to this pool.** Regeneration is byte-identical: the generating
config is a pure function of byte-unchanged constants and sampling is counter-based Philox keyed
per sample index. Restoration cost 52 min and was verified exact (all four realized image-size
counts reproduced). So this instance cost only compute — but an artifact with any nondeterminism
would have been unrecoverable.

---

## Structural: STATE.md `## Current Position` is single-slot

**Found:** 2026-07-27, during concurrent execution of phases 11 and 13.

**What happened.** Phases 11 and 13 executed in parallel. The block holds one position, so
phase 13's executor overwrote phase 11's line and STATE.md showed no trace of an active, blocked
phase 11. Restored additively in commit `40dbb4b` (10 insertions, 0 deletions, phase 13's block
byte-unchanged) with a note that both entries are authoritative — but that is a convention, not a
mechanism, and the next writer can still clobber it.

**Proposed fix.** Make the block **per-phase** rather than single-slot: one subsection per active
phase, keyed by phase number, with writers required to update only their own subsection.

```markdown
## Current Position

### Phase 11 — registration-and-chromatic-uncertainty-as-latent
Status: BLOCKED | Plan: 7 of 11

### Phase 13 — three-hypothesis-amortized-bayes-factor
Status: EXECUTING | Plan: 6 of 16
```

This makes concurrent writes non-overlapping by construction instead of relying on each agent's
discipline. **Not restructured mid-flight** — doing so while two phases are actively writing
STATE.md would cause exactly the contention it aims to prevent. It should be done between
milestones, when no phase is executing.
