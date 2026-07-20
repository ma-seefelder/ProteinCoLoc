---
phase: 8
slug: external-physical-ground-truth-corpus
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-07-20
---

# Phase 8 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Derived from `08-RESEARCH.md` → "## Validation Architecture".

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (`@testset`/`@test`) — the repo's existing pattern (`spike/test/runtests.jl`, `src`/`test` runners) |
| **Config file** | none — Julia `Test`; new corpus tests wired into a `corpus/test/runtests.jl` (Wave 0) |
| **Quick run command** | `julia --project=. corpus/test/runtests.jl` (offline subset — no network) |
| **Full suite command** | same, with network available (adds one live small-CBS-fetch smoke) |
| **Estimated runtime** | ~30 seconds offline; +~15s for the online CBS-fetch smoke |

---

## Sampling Rate

- **After every task commit:** Run the offline unit subset (manifest schema, tier guard, holdout guard, hashing, skip-with-flag, fixture conversion) — no network, always green.
- **After every plan wave:** Run the full offline suite.
- **Before `/gsd:verify-work`:** Offline suite green; live CBS-fetch smoke green when online (skip-with-flag otherwise).
- **Max feedback latency:** ~30 seconds

---

## Per-Task Verification Map

| Success Criterion / Decision | Behavior | Test Type | Automated Command | File Exists | Status |
|---------|----------|-----------|-------------------|-------------|--------|
| SC3 | Manifest schema valid (required columns; tier ∈ {physical-primary, simulated-secondary}; split ∈ {sealed_holdout, dev, eval}) | unit | `runtests.jl` manifest testset | ❌ W0 | ⬜ pending |
| SC2 / D-06 | Tier guard: default accessor never blends tiers; mixing requires explicit `mixed_tier` opt-in | unit | tier-separation testset | ❌ W0 | ⬜ pending |
| D-09 | Split guard: default corpus iterator never yields a `sealed_holdout` row; anchors only via `open_sealed_holdout` | unit | holdout-guard testset | ❌ W0 | ⬜ pending |
| D-05 | Hash MISMATCH = hard error (flip a byte on a fixture → expect `error`) | unit | integrity testset | ❌ W0 | ⬜ pending |
| D-05 | Unreachable = skip-with-flag (bad URL → `(status=:skipped)`, green) | unit | skip-with-flag testset | ❌ W0 | ⬜ pending |
| SC1 / D-04 | Content hash deterministic + names artifact (same bytes ⇒ same sha256) | unit | hashing testset | ❌ W0 | ⬜ pending |
| SC1 / D-03 | Fixture two-channel TIFF → `MultiChannelImage` → `patch(·,8)`/`correlation` apply unchanged, finite result | unit | conversion testset (synthetic fixture, no network) | ❌ W0 | ⬜ pending |
| D-04 | No image bytes staged in git (corpus data dir gitignored) | unit | gitignore/staging assertion | ❌ W0 | ⬜ pending |
| SC2 / D-07 | CBS live fetch smoke (one small set) verifies + tier-tags — SKIPS cleanly offline | integration | full-suite only | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `corpus/test/runtests.jl` — the testsets above
- [ ] A committed **synthetic two-channel TIFF fixture** (generated in-test or tiny tracked file) so conversion + hashing + integrity tests run with ZERO network
- [ ] `.gitignore` scoped so the corpus **data** dir is excluded but the manifest / fetch-script / tests are tracked (mirror the Phase-2 scoped-negation idiom)
- [ ] Confirm the corpus module's environment (reuse spike/root deps — no new adds expected; `SHA`/`Downloads` are stdlib)

---

## Manual-Only Verifications

| Behavior | Criterion | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Physical-anchor accession selection (positive tandem-FP + matched segregated) meets D-03 predicate (biological coloc/segregation label + open license) | SC1 / D-01 / D-02 / D-03 | Requires human judgement on biological ground-truth label + license terms; no automated oracle. Gate with `checkpoint:human-verify`. | Confirm each pinned accession/DOI resolves, carries an unambiguous biological coloc/segregation label in its source publication, and has an open/redistributable license before recording its sha256 in the manifest. |

---

## Validation Sign-Off

- [ ] All tasks have automated verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (synthetic fixture, runtests.jl, scoped .gitignore)
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
