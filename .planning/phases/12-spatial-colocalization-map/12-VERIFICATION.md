---
phase: 12-spatial-colocalization-map
verified: 2026-08-03T11:07:19Z
status: gaps_found
score: 8/9 must-haves verified (6 verified, 2 accepted by recorded user ruling as deferred-to-v2.1, 1 failed)
has_blocking_gaps: true
overrides_applied: 2
overrides:
  - must_have: "SPAT-07: Per-region calibration is reported as atom-free Gaussian-space SBC plus randomized-rank ρ-space SBC, with nuisance-appropriate equivalence testing on non-identified rows."
    reason: "User ruling, 2026-07-31 (12-D13-AUTHORISATION.md §9): the authorised 12-16 scope is SPAT-05 + SPAT-06 only; SPAT-07 is deferred to v2.1. Recorded on the record as `spat07_scope = :deferred` in the 12-16 artifact, not as a silent omission."
    accepted_by: "user (recorded in 12-D13-AUTHORISATION.md §9, ruling dated 2026-07-31)"
    accepted_at: "2026-07-31T00:00:00Z"
  - must_have: "SPAT-08: The chromatic-radial confound, the offset-grid artifact and the ε = 0 ablation are measured and reported as named guards."
    reason: "User ruling, 2026-07-31 (12-D13-AUTHORISATION.md §7): SPAT-08 is deferred to v2.1 on the Stage-2 descope path, extending the deferral 12-11-PLAN.md:255-259 pre-authorised for a Stage-1 descope. 12-06 delivered one precursor measurement (chromatic-ε linear identifiability, in scope, not deferred) but the full three-guard suite (radial-energy, offset-grid, ε=0 evaluation ablation) requires 12-20, which the user did not authorise."
    accepted_by: "user (recorded in 12-D13-AUTHORISATION.md §7, ruling dated 2026-07-31)"
    accepted_at: "2026-07-31T00:00:00Z"
gaps:
  - truth: "SPAT-06: Leave-region-out predictive coverage on real microscopy images is measured for the spatial model against a matched retrained ablation, on both calibration and a proper scoring rule."
    status: failed
    severity: blocking
    reason: "Two of the requirement's three named elements did not happen: (1) 'real microscopy images' — 12-19, the real-image arm, never ran; it is hard-blocked on a missing artifact (p12_train_full_report.jld2, confirmed absent on disk) that only the never-executed 12-17 would have produced. (2) 'the spatial model' — no spatial model exists; the mini-spike returned NONE-BEATS-ABLATION in all three runs (12-15, plus two reruns), so the only trained deliverable is the neutralized ablation. What WAS delivered (12-16, in scope, not deferred) is a leave-region-out predictive-coverage measurement on SIMULATED data comparing the ablation against a prior-only floor (no spatial arm to compare). That reduced measurement's own pre-registered criterion also failed: pooled coverage 0.9707 lies outside the band [0.87, 0.93] by a between-dataset-cluster interval of [0.9675, 0.9738], clear of the band by roughly 0.037. 12-16-SUMMARY.md states this in its own words: 'SPAT-06 IS NOT MET.' This is a genuinely reported negative result, not a stub or a hidden gap, and per the user's own recorded discipline (12-D13-AUTHORISATION.md §10.1, §10.5) it is NOT to be re-tuned by adjusting epochs, N, thresholds or the noise model."
    artifacts:
      - path: "spike/validation/p12_coverage.jl"
        issue: "The predictive-coverage machinery (p12_coloc_map, the D-09 leave-region-out construction) exists and runs correctly, but it was only ever exercised on the ablation arm against a simulated joint, never against a spatial arm and never on the six committed real TIFFs."
      - path: "spike/validation/p12_coverage_sim_report.jld2"
        issue: "The artifact this claim rests on (generated 2026-07-31T23:46, per 12-16) records coverage=0.9707, outside [0.87, 0.93] — the pooled, gated quantity fails its own pre-registered band."
    missing:
      - "A trained spatial-arm bundle to compare against the ablation (none exists at any scale — the mini-spike never selected one)."
      - "Execution of 12-19 against the six real committed TIFFs (test/test_images/{positive,negative}/*_c{1,2,3}.tif)."
    note: "This gap is already fully adjudicated and closed on the record — it is NOT presented here as new work to schedule. 12-D13-AUTHORISATION.md §10.5 records a user ruling ('Document only') explicitly declining further training or threshold changes to chase a passing coverage band. This verification report does not recommend re-running, retuning, or executing 12-17/12-18/12-19/12-20 to close this gap; it reports the true state so it is not silently read as satisfied."
deferred: []
human_verification: []
---

# Phase 12: Spatial Colocalization Map — Verification Report

**Phase Goal:** "Replace exchangeable patch pooling with a spatial lattice prior over the correlation
grid, producing an amortized per-region Δρ map with calibrated per-region uncertainty — the feature
that makes v2.0 'spatial' and differentiates it from Tapqir"

**Verified:** 2026-08-03T11:07:19Z
**Status:** gaps_found
**Re-verification:** No — initial verification

---

## Headline answer, stated first because a reader should not have to find it

**NO — the phase did not achieve its stated goal.** A per-region Δρ map mechanism was built and
works (`p12_coloc_map`, verified below), but two of the three things the goal sentence actually
promises are false on the record:

1. **"Spatial"** — no spatial prior (CAR or GP) ever beat the neutralized, non-spatial ablation.
   `NONE-BEATS-ABLATION` fired in three independent mini-spike runs (`12-MINISPIKE-VERDICT.md`,
   `12-15-RERUN-VERDICT.md`). The shipped deliverable is the ablation model — the *same* network
   trained with spatial borrowing switched off — wearing the `SpatialColocResult` type. There is no
   trained spatial arm anywhere in this phase's artifacts, at any scale.
2. **"Calibrated per-region uncertainty"** — `12-D13-AUTHORISATION.md` §9.1 states this explicitly:
   *"THE PHASE MAKES NO PER-REGION CALIBRATION CLAIM BEYOND PREDICTIVE COVERAGE."* SBC-based
   calibration (SPAT-07) is deferred to v2.1 by recorded user ruling and never ran. The one
   calibration-adjacent claim that *was* in scope for this phase — leave-region-out predictive
   coverage (SPAT-06) — was measured and **failed its own pre-registered band** (0.9707 against
   [0.87, 0.93]). `12-16-SUMMARY.md` records this in its own words as "SPAT-06 IS NOT MET."

**This is not read here as an execution failure.** Every one of these outcomes was reached through a
pre-registered, executed measurement, and every deviation from the original 20-plan schedule
(12-17..12-20 never running) is the direct, explicit consequence of recorded user rulings dated
2026-07-31 and 2026-08-03 (`12-D13-AUTHORISATION.md`, `12-STAGE1-VERDICT.md`,
`12-MINISPIKE-VERDICT.md`, `12-15-RERUN-VERDICT.md`). The project's own record already states the
honest conclusion in stronger language than this report needs to add: *"no arm was calibrated at this
scale"* and *"we cannot exclude that longer training would have brought the spatial arms into the
coverage band"* — the latter fixed as the exact sentence to carry into the manuscript
(`12-D13-AUTHORISATION.md` §10.5). No re-run, retune, or further phase-12 execution is recommended
by this report.

---

## Per-requirement verification (SPAT-01 .. SPAT-09)

### SPAT-01 — lattice prior + marginal-preserving copula (D-03, D-05, R-8)

**DELIVERED.** Both prior arms are implemented as dense linear algebra, not `AbstractGPs.jl` (per
`12-SC3-AMENDMENT.md` §2, a pre-registered dependency-resolution decision):
`car_sigma(G, α)` at `spike/simulator/p12_lattice.jl:125` and `gp_sigma(G, ℓ; kernel=:matern32)` at
`:153`. The marginal-preserving copula (`ghat(quantile(MU_PRIOR, Φ(z)))`) is exercised per-region and
validated against the pre-registered `P12_SIM02_W1_TOL_PERREGION = 0.10` bar at M = 20,000 across all
11 arm×rung configurations: max per-region W1 = **0.00675** (`spike/validation/p12_sim02_report.jld2`,
`12-10-SUMMARY.md`) — roughly 15× headroom, so SIM-02 holds per region as required.

The "chosen by the pre-registered mini-spike" clause resolved to **neither** arm: the mini-spike ran
three times (unseeded, seeded 18-epoch control, seeded 100-epoch treatment) and returned
`NONE-BEATS-ABLATION` all three times (`12-MINISPIKE-VERDICT.md`, `12-15-RERUN-VERDICT.md`). That is a
valid branch of the pre-registered `select_prior` rule, not an unbuilt apparatus, so the requirement
(build the arms, run the pre-registered comparison) is DELIVERED even though the comparison's answer
was negative.

### SPAT-02 — spatially-varying forward simulator, bit-for-bit at constant ρ (D-06)

**DELIVERED.** `spike/simulator/p12_lattice.jl`'s `_p12_upsample_local` (bilinear/separable) and
`spike/simulator/forward.jl`'s stage-1 optional `rho_field` binding implement the spatial rendering.
Bit-for-bit reproduction is proven against a pre-edit golden fixture captured one commit earlier
(`spike/test/fixtures/p12_stage1_golden.jld2`, `spike/test/capture_p12_golden.jl`): a θ with no
`rho_field` matches the golden bytes **12/12** exactly across six fixture keys, and a constant-valued
field matches the golden **exactly** (`maximum(abs, Δ) = 0.0`) against a tolerance of `1e-12`
(`12-08-SUMMARY.md`). Note `forward.jl` here is `spike/simulator/forward.jl` — a spike-local file, not
`src/`.

### SPAT-03 — CNN summary net reading the `G×G×2` lattice (D-02, R-7)

**DELIVERED.** `build_p12_summary_net` at `spike/npe/p12_architecture.jl:253` is the lean CNN specified
by `12-SC3-AMENDMENT.md` §1: `Conv((3,3),2=>16,gelu)` → `Conv((3,3),16=>32,gelu)` →
`Conv((1,1),32=>8,gelu)` channel bottleneck. Summary content is unchanged (128-dim: 64 continuous +
64 mask rows via `reshape_summary`), matching D-02's requirement that only the network's *view*
changes, not the summary's information content.

### SPAT-04 — global term + DCT-II per-region deviation field, correlation length in θ (D-01, D-04, D-07, D-08, R-1, R-2)

**DELIVERED, with a documented and repaired downstream defect noted for the record.** `P12_THETA_ROWS`
/ `p12_theta_index` (`spike/npe/p12_architecture.jl:343,361`) assemble the global term, the per-region
deviation field in the DCT-II basis (`p12_dct_vec` / `p12_idct_vec`, exact-invertibility asserted:
`p12_idct∘p12_dct == id`), and the correlation length r₁ as part of θ. A wrong-basis defect **was**
found in a downstream *scorer* (`run_p12_minispike.jl:326` applied `p12_idct_vec` directly to a
permuted θ column, never inverting the DCT-order permutation) — this affected reported skill/coverage
numbers in the mini-spike verdicts but did **not** touch the trained bundles or the head architecture
itself (training consumes θ as stored and never reconstructs a field). Repaired at `f039729`,
documented as a named family instance (`12-STAGE1-VERDICT.md` §7.4(a), §10.4).

### SPAT-05 — spike-local `coloc_map` returning per-region Δρ + per-region uncertainty (SC2 as scoped)

**DELIVERED.** `p12_coloc_map` at `spike/validation/p12_coverage.jl:990` returns a
`SpatialColocResultSpike` (`spike/p12/result.jl:231`) carrying `region_delta_rho` (primary
deliverable), `region_sd` (the per-region uncertainty on the difference), plus the exposed
`region_rho_sample` / `region_rho_control` maps and their own uncertainties. Each of the sample and
control stacks is one amortized forward pass — verified structurally: exactly two call sites of a
single `_p12cov_single_stack` helper (never two inlined copies) inside `p12_coloc_map`, with exactly
one `sampleposterior(` call inside that helper (`12-16-SUMMARY.md` §9, Testset 8). Scoped to `spike/`
per `12-SC3-AMENDMENT.md` §3 (a pre-registered lane decision, not a shortcut) — `src/results.jl`
stays byte-unchanged as its comment-only sketch.

**Named caveat that travels with this artifact:** the model behind `p12_coloc_map` in every artifact
this phase produced is the non-spatial D-13 ablation, not a trained spatial arm — the mechanism is
verified; the underlying model is not "spatial" in the sense the phase goal means.

### SPAT-06 — leave-region-out predictive coverage on real images, spatial vs. matched ablation (D-09, D-10, D-11, R-4)

**FAILED, as specified.** See the `gaps` entry in the frontmatter for full detail. Summary: no real
microscopy images were scored (12-19 never ran — hard-blocked on a missing artifact only the
never-run 12-17 would produce); no spatial model exists to compare against the ablation
(`NONE-BEATS-ABLATION`); the reduced substitute that *did* run (12-16, in scope, simulated data,
ablation vs. a prior-only floor) failed its own pre-registered coverage band (0.9707 outside
[0.87, 0.93]). `12-16-SUMMARY.md` states plainly: "SPAT-06 IS NOT MET."

### SPAT-07 — atom-free Gaussian-space SBC + randomized-rank ρ-space SBC + nuisance equivalence testing (R-2, Pitfalls 7-8)

**DEFERRED-TO-v2.1 — PASSED (override).** User ruling 2026-07-31 (`12-D13-AUTHORISATION.md` §9).
Confirmed on disk: `spike/test/test_p12_sbc.jl` (39 lines) is a deliberate, self-declaring empty
scaffold — `# P12_PENDING_SCAFFOLD, implemented by plan 12-18`, `@info "P12_PENDING_SCAFFOLD: ...
plan 12-18 implements the row-partition, randomized-rank and uniformity tests."` No `p12_sbc.jl`
runner exists anywhere in `spike/`. 12-18 never ran. The deferral is explicit, dated, and recorded as
`spat07_scope = :deferred` in the 12-16 artifact, not a silent gap.

### SPAT-08 — chromatic-radial confound, offset-grid artifact, ε = 0 ablation as named guards (S-4, D-06)

**DEFERRED-TO-v2.1 — PASSED (override), with one precursor measurement genuinely delivered.** User
ruling 2026-07-31 (`12-D13-AUTHORISATION.md` §7). 12-06 (in scope, not deferred) measured chromatic-ε
linear identifiability on the existing Phase-11 pool: `eps_ratio = 0.96307` (barely identifiable, 3.7%
better than the empirical prior-mean baseline) against the pre-registered vacuity floor
`P12_VACUOUS_SHRINKAGE_FLOOR = 0.90` — `ridge_residual_shrinkage = 0.96870 > 0.90` ⇒ `vacuous = true`
(`spike/validation/p12_eps_ridge_report.jld2`, `12-06-SUMMARY.md`). That is a real, in-scope
measurement, but it is only one narrow piece of S-4; the three named guards SPAT-08 actually requires
(radial-energy report, offset-grid arm, ε = 0 evaluation ablation) all live in 12-20, which never ran.
No `p12` file under `spike/` implements a radial-energy or offset-grid guard.

### SPAT-09 — two-stage pre-registered descope trigger governing the phase; on descope, ablation ships (D-12, D-13, R-5)

**DELIVERED, with a documented gate-wiring defect the user had to route around manually.** The net
outcome is real: the phase *was* governed by a two-stage descope logic, and the ablation model *did*
ship as the deliverable — `spike/validation/p12_ablation_report.jld2` and
`spike/artifacts/p12/p12_ablation_none.jld2` exist on disk and are the artifacts `p12_coloc_map`
consumes. **But the automated code path that was supposed to carry this never fired.** `12-CONTEXT.md`
defines "D-12 Stage 1" as the CAR-vs-GP mini-spike itself; the *implemented* Stage-1 gate
(`p12_stage1_verdict()`, read from `12-STAGE1-VERDICT.md`) is a different measurement (12-11's ridge
probe) that returned `:proceed`. When the mini-spike then returned `NONE-BEATS-ABLATION` — precisely
the condition `12-CONTEXT.md` names as the descope trigger — no plan step's routing logic recognized
it: `12-14` Task 3 (the only coded producer of the descope deliverable) is keyed on
`p12_stage1_verdict() === :descope` and is a no-op under `:proceed`; `12-17` is impossible as written
(keys on `P12_CHOSEN_PRIOR`, which is asserted absent); `12-16` is specified to throw with neither
bundle present. This divergence is documented at length as a named recurring defect family
(`12-STAGE1-VERDICT.md` §7.2, §7.6) and was resolved **not by fixing the code**, but by a standalone,
explicitly-dated user authorization (`12-D13-AUTHORISATION.md`) that created the D-13 ablation bundle
as a new record without touching the gate files. The requirement's net effect ("on descope, the
ablation ships") is achieved; the mechanism that was supposed to achieve it automatically is broken
and undocumented as fixed — it remains broken on the record by deliberate choice (`12-D13-AUTHORISATION.md`
§1: "not touched").

---

## Hard project constraints — verified independently

| Constraint | Check | Result |
|---|---|---|
| `src/` untouched by the spike | `git log --oneline -- src/ \| head` | Last commit touching `src/` is `ca02b0e` (2026-07-25, Phase 11-02) — **before** Phase 12 began (2026-07-29). No Phase-12 commit touches `src/`. **VERIFIED.** |
| `spike/Project.toml` / `spike/Manifest.toml` unchanged by Phase 12 | `git log --oneline -5 -- spike/Project.toml spike/Manifest.toml` | Last touch `7b97369` (Phase 9-01). **VERIFIED byte-unchanged since before Phase 12.** |
| `P12_CHOSEN_PRIOR`, `P12_N_LOW`, `P12_K_PROD` still absent/asserted-absent | `grep -n "P12_CHOSEN_PRIOR\|P12_N_LOW\|P12_K_PROD" spike/validation/p12_consts.jl` | All three appear only in comments and in `@assert !isdefined(...)` guards (`p12_consts.jl:798-799`), never as `const` definitions. **VERIFIED absent, exactly as the D-13 route requires.** |
| Phase-11 pool untouched | `12-11-SUMMARY.md` self-check | Byte-identical sizes and mtimes across all six shard files, before/after. |

---

## Required Artifacts

| Artifact | Expected | Status | Details |
|---|---|---|---|
| `spike/simulator/p12_lattice.jl` | CAR + GP kernels, copula | ✓ VERIFIED | `car_sigma`, `gp_sigma`, `_p12_upsample_local` all present and substantive (479 lines) |
| `spike/simulator/forward.jl` (spike-local) | spatial ρ rendering, bit-identical fallback | ✓ VERIFIED | Golden-fixture-proven 12/12 exact match |
| `spike/npe/p12_architecture.jl` | lean CNN summary net, θ row assembly | ✓ VERIFIED | `build_p12_summary_net`, `P12_THETA_ROWS` present (399 lines) |
| `spike/npe/train_p12_npe.jl` | training pipeline | ✓ VERIFIED | 626 lines; exercised by 12-14/12-15/D-13 runs |
| `spike/validation/p12_coverage.jl` | `p12_coloc_map`, D-09 LRO construction | ✓ VERIFIED | 1168 lines; `p12_coloc_map` at :990 |
| `spike/p12/result.jl` | `SpatialColocResultSpike` | ✓ VERIFIED | struct present with all three maps + uncertainties |
| `spike/validation/p12_stage1_report.jld2` | Stage-1 gate measurement | ✓ VERIFIED (present, adjudicated by user ruling) | |
| `spike/validation/p12_minispike_report*.jld2` (+ 4 rerun/repair variants) | mini-spike measurements | ✓ VERIFIED (present, 5 runs total) | |
| `spike/validation/p12_ablation_report.jld2`, `spike/artifacts/p12/p12_ablation_none.jld2` | D-13 descope deliverable | ✓ VERIFIED | present, consumed by 12-16 |
| `spike/validation/p12_coverage_sim_report.jld2` | SPAT-06 LRO measurement | ✓ VERIFIED present, ✗ FAILS its own band | coverage 0.9707 outside [0.87, 0.93] |
| `spike/test/test_p12_sbc.jl` | SPAT-07 SBC tests | ✗ SCAFFOLD ONLY | self-declared `P12_PENDING_SCAFFOLD`, 12-18 never ran |
| `spike/validation/run_p12_*` for radial-energy/offset-grid/ε=0 guards | SPAT-08 full guard suite | ✗ MISSING | only the ε-identifiability precursor (12-06) exists; 12-20 never ran |
| `p12_train_full_report.jld2` | production-scale (50k) spatial bundle | ✗ MISSING (confirmed by search) | required by 12-19; hard-blocks the real-image arm |

---

## Data-Flow Trace (Level 4) — `p12_coloc_map`

| Artifact | Data source | Produces real data | Status |
|---|---|---|---|
| `SpatialColocResultSpike.region_delta_rho` | `_p12cov_single_stack(bundle, mci)` → `sampleposterior` on a trained NPE bundle | Yes — bundle is the real (non-stub) D-13 ablation net, verified sha256-checked before scoring (`12-16-SUMMARY.md` §10) | ✓ FLOWING (mechanically); model behind it is the non-spatial ablation, not a spatial arm — see SPAT-05 caveat |
| `p12_coverage_sim_report.jld2` coverage figures | fresh draws through `p12_rng(P12_COVERAGE_COUNTER)` on the validation salt (structurally disjoint from the training/datagen salt, asserted at load time) | Yes — real simulated draws, not cached training data | ✓ FLOWING, but the number produced fails the pre-registered band |

---

## Requirements Coverage

| Requirement | Source Plan(s) | Status | Evidence |
|---|---|---|---|
| SPAT-01 | 12-03, 12-07, 12-09, 12-10, 12-15 | ✓ SATISFIED | see per-requirement section above |
| SPAT-02 | 12-07, 12-08, 12-09 | ✓ SATISFIED | golden-fixture bit-for-bit proof |
| SPAT-03 | 12-04, 12-09, 12-14, (12-17 not run) | ✓ SATISFIED | CNN present, summary content unchanged |
| SPAT-04 | 12-03, 12-04, 12-13, 12-14, 12-15, (12-17 not run) | ✓ SATISFIED (downstream scorer defect found + repaired, head unaffected) | DCT-II basis, exact invertibility asserted |
| SPAT-05 | 12-12, 12-16 | ✓ SATISFIED | `p12_coloc_map` verified structurally and by data-flow trace |
| SPAT-06 | 12-02, 12-16, (12-19 not run) | ✗ BLOCKED | real-image arm never ran; no spatial model exists; delivered substitute fails its own band |
| SPAT-07 | 12-16, (12-18 not run) | PASSED (override) | deferred to v2.1, user ruling 2026-07-31 |
| SPAT-08 | 12-06, (12-20 not run) | PASSED (override) | deferred to v2.1, user ruling 2026-07-31; ε-identifiability precursor delivered |
| SPAT-09 | 12-01, 12-02, 12-05, 12-11, 12-14, 12-15 | ✓ SATISFIED (net effect achieved; automated routing mechanism documented-broken) | D-13 ablation shipped via standalone user authorization |

No orphaned requirements: all nine SPAT IDs are claimed by at least one plan, and REQUIREMENTS.md's
Phase-12 mapping (`SPAT-01..SPAT-09`) matches exactly.

---

## Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|---|---|---|---|---|
| `spike/test/test_p12_sbc.jl` | 21 | `P12_PENDING_SCAFFOLD` marker, no failing assertion | ℹ️ Info | Deliberate, self-documenting, deferred-by-ruling scaffold — not a defect, already accounted for as SPAT-07's override |

No `TBD`, `FIXME`, or `XXX` markers found in any of the phase's core implementation files
(`p12_lattice.jl`, `p12_prior.jl`, `p12_architecture.jl`, `train_p12_npe.jl`, `p12_generate.jl`,
`p12_coverage.jl`, `result.jl`, `run_p12_minispike.jl`, `run_p12_ablation.jl`) — swept explicitly, zero
matches.

---

## Behavioral Spot-Checks

Not run as live commands — this is a Julia research spike whose "runnable" outputs are the committed
`.jld2` artifacts, and the project's own multiply-repeated, independently re-derived (never
copy-pasted-from-stdout) numbers across `12-STAGE1-VERDICT.md`, `12-MINISPIKE-VERDICT.md`,
`12-15-RERUN-VERDICT.md`, and `12-16-SUMMARY.md` already constitute exactly this kind of check, cross-
verified here against the artifact set and source files. Re-running the Julia suite was judged
unnecessary and out of scope for this verification (no code changes to validate; artifacts and their
sha256/mtime provenance already checked above).

## Probe Execution

SKIPPED — no `scripts/*/tests/probe-*.sh` convention in this repository; the phase's own runners
(`run_p12_*.jl`) are the analogous mechanism and are addressed under Required Artifacts / Behavioral
Spot-Checks above.

## Human Verification Required

None. Every must-have in this phase is a numeric/artifact claim resolvable from committed `.jld2`
files and source code, not a visual or UX judgment.

---

## Gaps Summary

**One blocking gap: SPAT-06.** Two of its three named elements (real microscopy images, a trained
spatial model) never happened, and the reduced substitute that did run failed its own pre-registered
coverage band. This is the direct, load-bearing reason the phase's stated goal — "calibrated
per-region uncertainty" — is not achieved: SPAT-06 was the only in-scope (non-deferred) calibration-
adjacent claim, and it failed.

**Two items resolved by recorded override, not by gap-closure work: SPAT-07, SPAT-08.** Both are
explicitly deferred to v2.1 by dated user rulings on the record (`12-D13-AUTHORISATION.md` §7, §9).
Treating these as ordinary "must fix" gaps would misrepresent a closed, human-adjudicated scope
decision as an open defect.

**This report does not recommend executing 12-17, 12-18, 12-19, or 12-20, and does not recommend
re-running, retuning, or extending any existing measurement.** The project's own record has already
closed these questions through explicit, dated user rulings — including a ruling (2026-08-03,
`12-D13-AUTHORISATION.md` §10.5) that declined two offered follow-up options ("one pre-declared
longer training run" and "a no-compute examination of the early-stopping traces") in favor of
documenting the limit as-is. The purpose of recording SPAT-06 as a blocking gap here is to state the
true state of the goal accurately, not to propose new work — this phase's own governance has already
decided that no new work is authorized.

---

_Verified: 2026-08-03T11:07:19Z_
_Verifier: Claude (gsd-verifier)_
