---
phase: 14
slug: decision-and-abstention-layer
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-08-03
---

# Phase 14 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Derived from `14-RESEARCH.md` § *Validation Architecture*. The SC → validation map below is the
> authoritative contract; the per-task map is populated by the planner and updated during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Julia stdlib `Test` (`@testset` / `@test`) — no install needed |
| **Config file** | none — plain scripts under `spike/test/` |
| **Quick run command** | `julia --project=spike spike/test/test_p14_<name>.jl` |
| **Full suite command** | ⚠ `julia --project=spike spike/test/runtests.jl` **exits 1 at the Phase-4 `SPEEDUP_GATE`, masking every later include block** (`deferred-items.md` D-13-A, `13-REPORT.md` §11 item 6). **Per-file runs are the only reliable signal.** Do not use `runtests.jl` as a gate. |
| **Per-wave command** | Loop every `spike/test/test_p14_*.jl` file individually; all must exit 0 |
| **Estimated runtime** | unit files < 10 s each (no simulation); integration runners are minutes (simulator draws) |
| **Fixture stream** | `p14_fix_rng(P14_FIXTURE_COUNTER)`, mirroring the `p13_fix_rng` pattern — Phase 14 needs its own `P14_FIX_SEED` so a quick gate never pre-observes a reported stream |

---

## Sampling Rate

- **After every task commit:** run the relevant `test_p14_*.jl` unit file(s) — each < 10 s, no simulation.
- **After every plan wave:** run **all** `test_p14_*.jl` unit files, per-file (never via `runtests.jl`).
- **Before `/bm:verify-work`:** all unit files green **plus** the **four** integration runners
  (`run_p14_fdr_check.jl`, `run_p14_conformal.jl`, `run_p14_riskcoverage.jl`, `run_p14_ood_arm.jl`)
  producing persisted `.jld2` artifacts with recorded seeds, counters and sha provenance.
  *(Reconciliation: `run_p14_fdr_check.jl` backs SC1-b/SC1-c and was missing from an earlier
  three-runner list. Four is the correct count.)*
- **Max feedback latency:** < 10 s for unit; integration runners are wave-boundary only.

---

## SC → Validation Map (authoritative)

SC1 and SC2 are **AMENDED** by D-02 and D-05 (`14-CONTEXT.md`, `14-DISCUSSION-LOG.md`). The rows below
validate the amended criteria, never the originals.

| ID | Behaviour | Measured on | Pre-registered number | Falsified if | Type | Command |
|---|---|---|---|---|---|---|
| **SC1-a** | `p14_bayes_fdr` implements the running-**mean** prefix rule | hand-built fixtures, no simulation | `v = [0.01, 0.02, 0.15]`, α = 0.10 ⇒ `k* = 3`, `fdp = 0.06`; `issorted(run)` on 10⁴ random ascending vectors | a discriminating fixture returns the running-**max** answer, or `run` is unsorted | unit | `test_p14_fdr.jl` |
| **SC1-b** | Realized FDP among DECIDED-and-accepted items tracks α_FDR | fresh unstratified draws, `n_eval = 2000` at `P14_EVAL_COUNTER` | over α_FDR ∈ {0.01, 0.05, 0.10, 0.20}: realized FDP within the binomial 95 % band of the rule's *predicted* `fdp` | realized FDP exceeds the upper band at any grid α | integration | `run_p14_fdr_check.jl` |
| **SC1-c** | Prior sensitivity of the FDR estimate is quantified, not hidden | same batch, π_C swept over {0.05, 0.10, 0.217, 0.40, 0.70} | **REPORTED, NOT GATED** — no bar | failing to report it is the failure | integration | same runner |
| **SC1-d** | Split-conformal marginal coverage ≥ 1 − α_conf | calibrate `n_cal = 2000` @ `P14_CAL_COUNTER`; measure `n_eval = 2000` @ `P14_EVAL_COUNTER`; both **unstratified** | coverage ≥ 1 − α_conf − 3·√(α(1−α)/n_eval); at α_conf = 0.10, n = 2000 ⇒ **≥ 0.880** | realized coverage falls below the band | integration | `run_p14_conformal.jl` |
| **SC1-e** | `p14_conformal_quantile` returns the **order statistic**, not an interpolated quantile | fixture | n = 9, α = 0.10 ⇒ `q̂ == maximum(s)` and `q̂ != Statistics.quantile(s, 0.9)` for non-degenerate `s`; at n = 30 the α = 0.03 call **throws** | either equality flips | unit | `test_p14_conformal.jl` |
| **SC1-f** ⟵ *rewritten per D-03a* | The real-data check runs on the **six committed microscopy TIFFs**, and the corpus is recorded as unavailable-by-design rather than silently skipped | `test/test_images/{positive,negative}/*_c{1,2,3}.tif`, read-only; `corpus/manifest.csv` **metadata only** (never `corpus/data/`) | the six TIFFs produce a recorded decision + abstention reason each; the corpus row records `corpus_images_available == 0` as **UNFETCHED BY DESIGN** (`.gitignore:451-453`) with the tier breakdown (30 `simulated-secondary`, 2 `physical-primary` = sealed holdout) | a coverage number is emitted from corpus data; the sealed holdout is read; the corpus is reported as "gone"/"empty" rather than unfetched; or the six-TIFF result is presented as a coverage claim | integration | `run_p14_real_images.jl` + `test_p14_decoupling.jl` |
| **SC2-a** | `p14_fuse` implements the D-05 asymmetry exactly | fixture — all 3 × 3 × 2 × 2 = 36 combinations | full truth table cell by cell; in particular `(:clear, :singleton, true, false) ⇒ (:decide, :decided_contra_classical)` | any cell differs | unit | `test_p14_fuse.jl` |
| **SC2-b** | D-06: missing / non-finite OOD threshold ⇒ `:not_checked` ⇒ abstain by default | fixture `ood_nulls` with `:thr` absent / `Inf` / `NaN` / finite | all three non-finite/absent cases ⇒ `:not_checked` ⇒ `:abstain`, `reason = :ood_not_checked`; opt-out flips to `:decide` | any case reads as `:clear` | unit | `test_p14_fuse.jl` |
| **SC2-c** | D-06 anti-regression: no source uses `verdict.flag == false` as an in-distribution proxy | grep of comment-stripped Phase-14 sources | zero hits for `!.*\.flag` / `\.flag == false` outside `p14_ood_state` | any hit | unit (source grep) | `test_p14_decoupling.jl` |
| **SC2-d** | Every classical-vs-amortized comparison passes through `ghat` | grep of comment-stripped Phase-14 sources | every `patch_correlation(` call site compared to τ is wrapped in `ghat(` | any bare comparison | unit (source grep) | `test_p14_decoupling.jl` |
| **SC3-a** | Confidence ordering beats random (selective skill) | `n_eval = 2000` unstratified | `skill = (AURC_rand − AURC)/(AURC_rand − AURC_oracle) ≥ P14_SKILL_FLOOR = 0.60` | skill below floor | integration | `run_p14_riskcoverage.jl` |
| **SC3-b** | Selective risk rises with coverage (monotone trend) | same | `spearman(coverage, selective_risk) ≥ 0.95` over coverage ∈ **[0.20, 1.00]**, floor frozen before the run | ρ < 0.95, **or** the coverage floor is changed after seeing the curve | integration | same runner |
| **SC3-c** | Abstention concentrates on hard cases | same | `AUC_hard ≥ 0.80` via the in-repo `roc_auc` | AUC below floor | integration | same runner |
| **SC3-d** | Abstention concentrates on OOD cases | `OOD_FAMILIES` at the strongest level vs in-distribution, matched n | `abstain_rate(misspec) − abstain_rate(ID) ≥ 0.50` | difference below floor | integration | `run_p14_ood_arm.jl` |
| **D-07** | τ is loaded, four-way agreed, sha-pinned | the four Phase-13 artifacts | `net.tau == probe.tau == gate.P13_TAU == P13_TAU == 0.15`; `net.consts_sha == P13_CONSTS_SHA.pre_amendment`; `p13_consts_sha() == .post_amendment`; **`probe.consts_sha256 ∉ values(P13_CONSTS_SHA)`** (the divergence is asserted, not papered over) | any equality fails, or the divergence assertion is widened to a disjunction | unit | `test_p14_provenance.jl` |
| **D-01/env** | `src/`, `spike/Project.toml`, `spike/Manifest.toml`, `corpus/`, `artifacts/` untouched; 16 deps exactly; no CUDA | git + `Pkg` | all assertions from `test_p12_decoupling.jl:152-208, 291-296, 351-353`, restated for Phase 14 | any diff | unit | `test_p14_decoupling.jl` |
| **Seeds** | P14 streams provably disjoint from every P13 and forbidden stream | constants | P14 seeds ∉ `_p13_forbidden()`; P14 counters pairwise distinct; **full cross-product of Philox key words pairwise distinct** vs the P13 family | any collision | unit | `test_p14_consts.jl` |
| **Class order** | Every named field lands on its intended class | fixture with deliberately **UNEQUAL** posteriors (0.7 / 0.2 / 0.1) | `p.coloc`, `p.random`, `p.exclusion` each equal their intended value | any permutation | unit | `test_p14_posterior.jl` |
| **Prior re-derivation** | `pi_class_masses` is re-derivable, not taken on trust | `sample_prior` + `three_way_label` on a fixture stream | re-derived masses match `h.meta.pi_class_masses` to the fixture's MC error | mismatch beyond MC error | unit | `test_p14_posterior.jl` |

**Why the class-order row exists:** five different class orderings coexist in this codebase
(`three_way_probs`, `three_way_log_bf`, `ThreeWayLogBF`, `confusion_class_order`, `pi_class_masses`).
ECE, AUC and FDR are all invariant to a *consistent* relabelling, so a positional-index bug is
invisible to every headline metric — this is the exact shape of the Phase-12 permutation artifact.
Equal-probability fixtures cannot catch it; the fixture must be unequal.

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| *(populated by gsd-planner — every task must map to at least one SC row above)* | | | | | | | | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

**Planner contract:** every task in every `14-*-PLAN.md` must name the SC row(s) it is verified by, and
no three consecutive tasks may run without an automated verify.

---

## Wave 0 Requirements

- [ ] `spike/p14/consts.jl` — the Tier-1 pre-registration: `P14_DEV_SEED`, `P14_FIX_SEED`,
      `P14_CAL_COUNTER`, `P14_EVAL_COUNTER`, `P14_OOD_COUNTER`, `P14_FIXTURE_COUNTER`,
      `P14_ALPHA_CONFORMAL`, `P14_N_CAL`, `P14_N_EVAL`, `P14_SKILL_FLOOR`, `P14_SPEARMAN_FLOOR`,
      `P14_COVERAGE_FLOOR`, `P14_AUC_HARD_FLOOR`, `P14_OOD_MARGIN_FLOOR`, `P14_ITERATION_ALLOWANCE`.
      **Frozen and committed in a commit that precedes every result commit**, using the Tier-1/Tier-2
      guard-block structure of `spike/p13/consts.jl` and the append-never-overwrite rule.
- [ ] `spike/test/test_p14_consts.jl` — literal-value + disjointness assertions (mirrors
      `test_p13_consts.jl`).
- [ ] `spike/test/test_p14_decoupling.jl` — environment / corpus / `src` guards (mirrors
      `test_p12_decoupling.jl`, with the Phase-14 allowlist and the two new source-grep checks
      SC2-c / SC2-d).
- [ ] `spike/test/test_p14_provenance.jl` — the D-07 four-way τ + sha table.
- [ ] `spike/test/fixtures/` — the unequal-posterior class-order fixture and the FDR discriminating
      fixture.
- [ ] No framework install needed (`Test` is stdlib).

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| The four SC3 bars (0.60 / 0.95 / 0.80 / 0.50) and the 0.20 coverage floor are **judgement calls with no derivation** | SC3-a..d | No principled derivation exists; this project has recorded four cases of an underived bar measuring the wrong thing | The numbers are frozen in `spike/p14/consts.jl` **before any result exists** and are labelled in the report as judgement calls. A user may overrule any of them, but only *before* the freeze commit — never after a result. |
| Assumption A2 — Bayesian posterior-expected FDP is immune to data-dependent pre-selection | SC1-b, D-05 | It is the entire justification for abstain-then-sort; standard conditioning, but load-bearing | The report must state the derivation explicitly so a referee can check it rather than take it on trust. |
| Real-data check substrate (**D-03a**, supersedes D-03's corpus half) | SC1-f | The corpus holds **no unsealed physical ground truth** — 30 of 32 rows are `tier: simulated-secondary` (the CBS benchmark, i.e. another simulator) and the only 2 `physical-primary` rows are Phase 16's `sealed_holdout`. Its bytes are also unfetched by design. So a corpus check would be neither physical nor runnable. | Use the **six committed microscopy TIFFs** under `test/test_images/`, read-only — the same substitution Phase 13 sanctioned at `spike/p13/real_images.jl:71,145` *"precisely so that Phase 16's blind corpus stays blind."* Report it as an **illustration, not a coverage claim**: six TIFFs in two conditions bounds nothing tightly, and the report must say so. The withdrawn `α ≥ 1/31 ≈ 0.032` figure must not appear anywhere. |

---

## Execution Environment Note

**Execute on the main working tree — no git worktrees.** The OOD null fit needs the ~50 MB gitignored
`cache/p13/` pool, which does not exist inside a fresh worktree and was destroyed once already this
milestone (`worktree-cleanup-destroys-gitignored-bulk-data`).

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] `runtests.jl` is NOT used as a gate anywhere (it exits 1 at the Phase-4 `SPEEDUP_GATE`)
- [ ] Feedback latency < 10s for unit files
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
