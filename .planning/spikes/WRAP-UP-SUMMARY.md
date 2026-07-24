# Spike Wrap-Up Summary

**Date:** 2026-07-24
**Spikes packaged:** 13 (001–013)
**Output skill:** `.claude/skills/spike-findings-proteincoloc/` (SKILL.md + 3 references + sources/)

The 13 spikes were synthesized into a persistent project skill so future build conversations inherit
the verified conclusions instead of re-deriving them. Findings are grouped into **three areas**.

## The three areas

### 1. `src/` refactoring roadmap — spikes 001–005 (`references/src-refactoring-analysis.md`)
A *separate, earlier* session: an analysis-only refactoring roadmap for the v1.0 `src/` core (5
dimensions, 46 verifier agents, hot-path micro-benchmarked). **Nothing implemented**; gated by the
"don't touch `src/` during the v2.0 spike" rule. Consolidated roadmap in `PROPOSAL.md`.
- **Δρ correctness bug** (`utils.jl:368`, found by 4/5 dimensions) — the headline effect size mixes
  posterior μ_sample with **prior** μ_control. 1-line fix, wrong scientific number.
- **Hot-path rewrite** — `correlation()` via a function barrier: **3.4–6.0× faster, 135–1394× less
  memory** (measured). Hoist `InterpKDE` out of the `quadgk` integrand.
- Real type hierarchy (`AbstractVector` stack, parametrized `CoLocResult` without renames); lift the
  Turing `@model` to module scope (enables v2.0 prior reuse); CairoMakie headless viz; seed +
  run-manifest provenance. **Posterior path uses removed Turing API — confirm it runs at all first.**

### 2. Bayes-factor ship-gate — spikes 006, 007, 008 (`references/bayes-factor-gate.md`)
Why the amended grid-8 gate's BF arm returned 58/100 finite pairs (< pre-registered 69).
- **006:** attrition is **100% baseline-side**, boundary-triggered (KDE `p_post` → {0,1} → ±Inf on
  one-sided posteriors), concentrated at large |Δρ|. The **NRE never** goes non-finite.
- **007:** survivors-only filter biases r̂ down ~0.03–0.05 and computes correlation on the ambiguous
  middle; the clamped all-pairs number is a two-valued saturated column, also unsound.
- **008 (PARTIAL):** the log-space remedy proved the **KDE baseline is the weak link, not the NRE** —
  the clamp *masked* a genuine order-of-magnitude NRE-vs-KDE tail disagreement. A KDE BF from L=999
  draws is **not a valid reference past |log-BF| ≈ log(L) ≈ 6.9**. Memo §5 "clamp artifact" framing
  is INCOMPLETE. Build guidance (restrict regime / raise L / compare sign-rank / reconsider the tol)
  is **unapplied** — a pre-registration matter under `07-GATE-AMENDMENT.md` §6.

### 3. SBC calibration of the NPE — spikes 009, 010, 011, 012, 013 (`references/sbc-calibration.md`)
Why the amended gate's SBC arm rejected every parameter.
- **009:** shrinkage measures **width, not location** — a vacuous parameter can still reject on a
  location shift. Two coexisting failure modes: location drift vs shape.
- **010 (+ADDENDUM):** ρ_true "overconfidence" **RETRACTED** — it is a **prior-atom artifact**
  (ghat clamps ~6.5% of Truncated-Cauchy μ draws onto exactly ±0.99; predicted 6.76% vs observed
  6.45%). Atoms excluded / randomized ranks → ρ_true **calibrated** (KS 0.34–0.77). **Coloc targets
  ρ_true/Δρ are calibrated.** Image size is a covariate via out-of-support generalization, not drift.
- **011:** residual nuisance failures = **0.06-SD marginal location drift on NON-identified
  nuisances**, exactly at the M=2000 over-power edge (M=2000 resolves 0.059 SD @ 50% power; M=500
  needs 0.118).
- **012 (PARTIAL):** raising flow capacity **REDISTRIBUTES** the drift (0.02–0.08 SD band unchanged),
  does not reduce it; targets stay calibrated. Not a reliable lever.
- **013:** μ-prior truncation removes the atoms (ρ_true SBC natively clean, KS 0.342) **but** doubles
  |ρ|≥0.95 bias and breaks ADVI comparability; nuisance drift unchanged. **Randomized ranks achieve
  the same SBC validity at neither cost.**
- **Recommended approach:** randomized ranks for atom-carrying params (ρ_true only), a
  nuisance-appropriate equivalence test for non-identified nuisances (spec draft
  `07-NUISANCE-SBC-SPEC-DRAFT.md`), strict point-null SBC at M=2000 for the coloc TARGETS only. §6.4
  blocks iterating the current amended_v2 gate; any change bundles into one final amendment on a
  retrained model.

## Honesty ledger

- **PARTIAL verdicts:** 008 (BF remedy failed — the failure is the finding), 012 (capacity
  redistributes rather than reduces).
- **Retraction:** 010's parent README "overconfident posterior" conclusion is superseded by its
  ADDENDUM (prior-atom artifact). Read the addendum as settled.
- **All gate remedies are UNAPPLIED** — pre-registration matters under `07-GATE-AMENDMENT.md` §6.
- **Single DEV seed** for each of 006–013 (asserted disjoint from all pre-registered seeds); the
  ~±1.5-z seed variance applies to per-parameter pass/fail. Robust results are the drift-magnitude
  bands and target calibration, not individual verdicts.
- **`src/` untouched** throughout 006–013 (custom datagen injected via the `datagen` seam; baseline
  artifacts verified byte-identical, npe_8 sha256 `198bb078…`).

## Artifacts produced by this wrap-up

- `.claude/skills/spike-findings-proteincoloc/SKILL.md`
- `.claude/skills/spike-findings-proteincoloc/references/{src-refactoring-analysis,bayes-factor-gate,sbc-calibration}.md`
- `.claude/skills/spike-findings-proteincoloc/sources/001-…/ … /013-…/` (READMEs + `.jl` + small
  figures/logs; `.jld2` data binaries excluded)
- `.planning/spikes/WRAP-UP-SUMMARY.md` (this file), updated `.planning/spikes/CONVENTIONS.md`
- CLAUDE.md routing line → `Skill("spike-findings-proteincoloc")`
