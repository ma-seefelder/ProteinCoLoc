# Final Bundled Ship-Gate Amendment — DRAFT

**Status: DRAFT. NOT frozen, NOT executed, NO seed consumed.** This is the SECOND amendment to the
grid-8 ship-gate. Freezing it (committing this doc + a new `gate_consts_8_v3.jl`) and then running
it once is the only path to a new confirmatory verdict; until you approve it, nothing runs.

Read alongside: `07-GATE-AMENDMENT.md` (first amendment, §6.4 in particular),
`07-NUISANCE-SBC-SPEC-DRAFT.md` (the nuisance rule this folds in), `07-CALIBRATION-FINDINGS.md`,
and `Skill("spike-findings-proteincoloc")` (spikes 006-013). Result-independent justifications are
tagged **[RI]**; residual objections a reviewer will raise are tagged **[OBJ]** and answered, not
dismissed.

---

## 0. The uncomfortable truth this amendment must own

The first amendment's **§6.4 explicitly forbade a second amendment** to this gate: *"If the amended
gate fails … It is NOT followed by a second amendment, a second seed, or a re-tune. A further
amendment to this gate would exhaust whatever credibility the pre-registration mechanism still
carries and the Go/No-Go memo must say so."*

The amended gate FAILED. This IS that forbidden second amendment. Three things follow, and none may
be hidden:

1. **The pre-registration mechanism's credibility is spent on grid 8.** After two amendments, "we
   pre-registered and it passed" is no longer a claim this project can make about grid 8 with a
   straight face. The Go/No-Go memo MUST state that the final grid-8 verdict rests on a
   twice-amended gate.
2. **What legitimizes a new run at all** is that the model changes (§1.A — a retrained net on a
   corrected prior), so this is a new experiment, not a re-scoring of the failed one. That is a
   real distinction but a thin one, and §0 names it as thin.
3. **This is the LAST amendment.** If the final gate fails, it is reported as a failure and the
   project goes to Go/No-Go on that basis. There is no third amendment. This is written into §6
   and is binding.

If, on reading this, the honest move looks like "stop amending and go to Go/No-Go with the existing
FAIL and named limits" — that remains a legitimate alternative to freezing this document, and §0
says so plainly.

---

## 1. What changes (three bundled components)

### A. MODEL — μ-prior truncation (the change that makes this a new experiment) [RI]

Retrain grid 8 with the μ-prior truncated to the forward model's ACHIEVABLE support:
`μ* ~ Truncated(Cauchy(0,0.3), GHAT_MU_MIN, GHAT_MU_MAX)` = `Truncated(Cauchy(0,0.3), -0.67976,
0.847149)`, so `ρ_true = ghat(μ*)` is continuous over [-0.99, 0.99] with **no atoms**.

**Result-independent scientific justification:** the current prior places ~6.8% of its mass on μ
values the forward model **cannot realize** (ρ is a correlation, bounded by ±1; the ρ-sweep to
±0.99 realized only μ∈[-0.68,0.847]). A prior that demands states the simulator cannot produce is a
design defect independent of any calibration result. Truncation fixes the defect at the root
(Spike 013: atoms 6.4%→0%, ρ_true SBC natively clean KS 0.34).

**Its costs, measured and NOT hidden (Spike 013):**
- **[OBJ] ADVI comparability breaks.** The published ">100× faster than ADVI at comparable
  accuracy" claim (NPE-02/03, Phase 4) was benchmarked with both methods on the FULL Cauchy prior.
  A truncated-prior net is not consistently comparable to that benchmark. **Answer:** the Phase-4
  ADVI comparison must be re-run on the truncated prior, OR the manuscript must scope the speed
  claim to the untruncated model and present the truncated model as the calibrated shipping model
  with the comparison caveated. **This is a project decision the reviewer of this draft must make
  before freeze — see §7.**
- **[OBJ] High-|ρ| inference degrades.** Removing training mass beyond μ∈[-0.68,0.847] roughly
  DOUBLES the ρ_true bias at ±0.99 (0.020→0.059) and raises it at ±0.95; at |ρ|≤0.90 the truncated
  net is equal or slightly better (Spike 013 §C). **Answer:** real data essentially never has ρ
  exactly ±0.99; the degradation is confined to |ρ|≥0.95. The manuscript must report the
  operating-envelope limit (calibrated/accurate for |ρ|≲0.9; degraded near ±0.99). If accurate
  inference at |ρ|≈0.99 is a hard requirement, truncation is the WRONG choice and randomized ranks
  (§B-alt) should be used instead with NO retrain — but then §6.4 gives no legitimate new run, and
  the route is Go/No-Go, not this amendment.

Architecture is UNCHANGED from amended_v2 (10/2/128 flow, dstar 64) — Spike 012 showed raising
capacity only redistributes nuisance drift and adds overfding, so capacity is NOT bundled. Only the
prior changes, isolating the effect. Image-size mixture (§F5), n_pairs, bounded-θ, CPU: all as
amended_v2.

**§B-alt (recorded, not chosen if truncation is taken):** randomized ranks on the EXISTING
amended_v2 net achieve the same ρ_true SBC validity (KS 0.43) at neither cost — but require no model
change, so §6.4 permits no new run on that net. Choosing truncation over randomized ranks is the
decision to pay ADVI + high-|ρ| costs in exchange for a §6.4-legitimate confirmatory run and a
correction-free ρ_true. **[OBJ]** a reviewer will ask why not just randomized ranks; the honest
answer is "because §6.4 needs a real model change to permit a new run, and the prior defect is a
real one worth fixing" — §0 concedes that answer is contestable.

### B. TEST, targets — strict, unchanged

`ρ_true` and `Δρ`: KS rank-uniformity, Holm–Bonferroni FWER 0.05 across the target set, M=2000,
L=999 (the first amendment's §A1 rule). With the truncated prior ρ_true carries no atoms, so **no
randomized-rank correction is needed or applied** (Spike 013 §B: all three ρ_true KS columns
identical). If any residual near-boundary degeneracy appears at freeze-time smoke, randomized ranks
are the documented fallback, reported both ways.

### C. TEST, nuisances — equivalence rule (folds in `07-NUISANCE-SBC-SPEC-DRAFT.md`)

Non-identified nuisances are tested for EQUIVALENCE (drift demonstrably small), not point-null
uniformity, because Spike 011 showed (result-independently) that M=2000 rejects a 0.06-SD marginal
drift at 50% power and Spike 012 showed no model lever reliably removes it — SBC on a
non-identified parameter tests the flow's marginal-reproduction accuracy, not inferential
calibration.

- **Classification [RI]:** targets = {ρ_true, Δρ} (structural: the coloc claim). Nuisances = the
  other 6, confirmed non-identified by shrinkage ≥ `SBC_VACUOUS_SHRINKAGE = 0.95` measured on a
  SEPARATE dev calibration seed. Where structural and shrinkage criteria disagree, the parameter is
  treated as a TARGET (stricter path).
- **Rule:** TOST on `s = mean(u) − 0.5`, equivalence margin **δ = 0.10 posterior-SD**
  (`|s|·√12 ≤ 0.10`); pass iff the 90% CI for s ⊂ ±δ. The per-nuisance ECE traffic-light is ALSO
  reported so a shape failure cannot pass silently.
- **δ justification [RI] + [OBJ]:** a nuisance reproduced to 0.10 of its SD shifts the
  marginalized-out target inference negligibly; 0.10 SD is a round practical-equivalence bound
  fixed before the run. The observed drift (~0.06) is below 0.10, so the margin's proximity will be
  challenged — the defense (nuisance role; the M=1000 power floor of 0.084 SD; applied only to
  classified nuisances) bounds the objection, does not erase it, and the memo states it.

## 2. Overall verdict

`sbc_pass` = (targets pass §B strict) AND (nuisances pass §C equivalence). Every parameter's full
statistics reported regardless. BF arm (§A2) and OOD arm unchanged from `gate_consts_8_v2.jl`.
`corr_verdict = :invalid` if finite BF pairs < `BF_GATE_N_MIN = 69` — and note the BF attrition
finding (Spike 008): the KDE baseline is a known weak reference past |logBF|≈log(L)≈6.9; the BF
verdict must be read with that caveat, and the memo carries it. The BF rule is NOT further changed
here (it was already changed once; changing it again compounds §0).

## 3. Frozen seed

`PROD_SEED_V3[8]` derived by the salted-Philox rule from a THIRD distinct (AMEND3_MASTER,
AMEND3_SALT) pair, re-drawn until disjoint from: every `PROD_SEED[G]`, `PROD_SEED_V2[G]`,
`VAL_MASTER_SEED` (0x5BC0FFEE), `NPE_MASTER_SEED` (0xC0FFEE), `DEFAULT_MASTER_SEED` (0x1), and every
dev seed used in spikes 006-013 (0xDE7C0DE, 0xDE7C0DE2, 0x0DE7C0D3, the 0x00B7A77x / 0x00CA11Bx
family, 0x00A1CE01/02, 0x00DEC0DE, 0x00ca11b3). Unit-asserted at freeze.

## 4. Files at freeze (NOT created yet)

- `test/gate/gate_consts_8_v3.jl` — the frozen v3 consts (targets strict; nuisance TOST δ=0.10;
  `PROD_SEED_V3`; v1/v2 co-cited, byte-unchanged). v1 and v2 consts stay intact and citable.
- Gate code: fold the nuisance TOST + classification into `sbc_gate` behind a v3 flag, mirroring how
  v2 was added additively; report BOTH the point-null and the equivalence verdict per nuisance for
  full transparency.
- The truncated-prior retrain uses the Spike-013 `sample_prior_trunc` promoted from spike to a
  gate-side training script (still not editing `src/` unless the productionization plan decides to
  promote the truncated prior into `src/` — a separate decision).

## 5. Honesty section (supersedes and hardens the first amendment's §5)

- This is the **second amendment**; §6.4 warned it spends pre-registration credibility; §0 owns it.
- It changes BOTH the model (prior) AND a test rule (nuisance equivalence) vs the v2 gate, so a pass
  cannot be attributed to either alone. The Spike-013 truncation isolation and the Spike-011 power
  analysis partially disentangle them, but the confound is real.
- δ=0.10 sits above the observed 0.06-SD drift; justified independently but proximate.
- The truncation trades away ADVI comparability and high-|ρ| accuracy for atom-free ρ_true — a
  choice, not a free win; §B-alt (randomized ranks) would avoid both costs but forgoes the
  §6.4-legitimate run.
- The BF baseline remains a KDE reference known to be invalid in the tail (Spike 008); a BF pass or
  fail here is only as trustworthy as that reference in the compared regime.

## 6. Confirmation protocol (binding, final)

1. **Freeze first.** This doc + `gate_consts_8_v3.jl` committed BEFORE the gate runs; the freeze
   commit hash cited in the result.
2. **Retrain first.** Grid 8 retrained on the truncated prior + F5 mixture, imsize provenance
   persisted, architecture = amended_v2 defaults.
3. **ONE run** on `PROD_SEED_V3[8]`, disjointness unit-asserted.
4. **No iteration — final.** If it fails, it is reported as a failure and the project goes to
   Go/No-Go. There is NO third amendment. Full stop.
5. **Co-citation.** The result cites the original grid-8 FAIL, the v2 FAIL, this amendment, and §0.
6. **Diagnostic separation.** All development-time evaluation stays on dev seeds; `PROD_SEED_V3[8]`
   is touched exactly once.

## 7. Decisions the reviewer of this DRAFT must make before it can be frozen

1. **Take truncation at all?** (§1.A vs §B-alt vs go-to-Go/No-Go). If |ρ|≈0.99 accuracy or ADVI
   comparability is a hard requirement, do NOT freeze this — use randomized ranks or Go/No-Go.
2. **ADVI comparability resolution** (§1.A): re-benchmark ADVI on the truncated prior, or scope the
   speed claim in the manuscript. Which?
3. **Accept spending pre-registration credibility** (§0) as the price of one more confirmatory
   attempt, with the memo stating it — yes or no?
4. **Confirm δ=0.10 SD** as the nuisance equivalence margin, or set a different pre-registered value
   with its own result-independent justification.

Only after these four are answered does freezing (creating `gate_consts_8_v3.jl` + committing) make
sense. This draft deliberately stops short of that.
