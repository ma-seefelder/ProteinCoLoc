# Phase 11 — D-12 Provenance, Decoupling and Shipped-Bundle Regression

**Run:** 2026-07-25 19:30–19:50 UTC, CPU-only, Windows 11 / Julia 1.12.6, single worktree
**Worktree HEAD at run time:** `8d7e45c63798ef035b7728c62148fc9e23a121d9` (plan 11-03, task 1 landed)
**Shas under audit:**

```
PHASE11_BASE_SHA / PRE_EPS_SHA (fixture-embedded pin)  17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348
PRE_EPS_SHA (immediate first parent of the ε commit)   a2ced26d93c7409a465aa50320e188868b3ab462
EPS_COMMIT_SHA (the D-11/D-12 commit)                  ca02b0e8a696269a2b8710b11ca3073edc5ae6c7
Artifacts.toml grid_8 pin (git-tree-sha1)              90e6b63a8a234d067b407fefd7914f2ae4845448
```

**Protocol:** every command below was executed and its **literal output** transcribed. Nothing here
is inferred from code reading. Where a check could not be run in the form the plan specified, the
form actually run is stated and the difference is called out rather than papered over. Empty output
is written as `(empty)`.

---

## Overall verdict: **PASS, with one pre-existing finding recorded (not caused by Phase 11)**

| Group | Check | Verdict |
|---|---|---|
| 1 | D-12 pinned-sha provenance | **PASS** — the pin is the ε commit's first parent, and both quoted shas predate `chromatic_eps` |
| 2 | The five §F20c decoupling commands | **PASS** — all four diffs empty, root suite green |
| 3 | Shipped-bundle regression | **PASS** — `estimator_for(8)` loads, persisted `D == 7`, `θzt` 7-element, `colocalization_amortized(...)` returns a finite Δρ on real images, `_theta_tuple(v7)` still feeds `simulate_pair` |
| 3-F | `Artifacts.toml` tree-sha1 vs the *installed store copy* | **FINDING (pre-existing, out of scope)** — the four bundle blobs are byte-identical to the pinned tree; the tree-sha1 differs by the git **file-mode bit only**, a Windows artifact of how the pin was generated. See §3.6 |
| 4 | The narrower Q8 decoupling claim | **ADOPTED**, with one correction to the plan's enumeration (a fourth file is in the diff) |

---

## 1. D-12 sha verification

### 1.1 The pin is the ε commit's first parent

```console
$ git rev-parse "ca02b0e8a696269a2b8710b11ca3073edc5ae6c7^"
a2ced26d93c7409a465aa50320e188868b3ab462
```

**Expected:** `PRE_EPS_SHA`. **Observed:** `a2ced26d93c7409a465aa50320e188868b3ab462` — the
first-parent sha recorded in `11-02-SUMMARY.md`. **PASS.**

Two shas are legitimately in play here and both are quoted in named limit #8, because 11-02 was
given two authorities that disagreed (the fixture-embedded `PHASE11_BASE_SHA = 17ebd1e…` and the
plan's "first parent of the ε commit" criterion = `a2ced26…`). They are reconciled by byte-identity,
proven directly:

```console
$ git diff --stat 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348 a2ced26d93c7409a465aa50320e188868b3ab462 -- src/amortized/simulator.jl
(empty)
```

`src/amortized/simulator.jl` is byte-identical at both shas, so
`git show <SHA>:src/amortized/simulator.jl` reconstructs the same training-distribution source
against either. The provenance-by-reference claim is therefore true as written under either reading.

### 1.2 Both quoted shas predate the parameter

```console
$ git show "17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348:src/amortized/simulator.jl" | grep -c chromatic_eps
0

$ git show "a2ced26d93c7409a465aa50320e188868b3ab462:src/amortized/simulator.jl" | grep -c chromatic_eps
0
```

**Expected:** `0`. **Observed:** `0` for both. **PASS.** The pinned source genuinely predates
`chromatic_eps`; the named limit does not point at a commit that already contained the parameter.

### 1.3 The shas are actually quoted in the named limit

```console
$ grep -c "17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348" docs/amortized.md
2

$ grep -c "a2ced26d93c7409a465aa50320e188868b3ab462" docs/amortized.md
1
```

**Expected:** at least 1 for the pin. **Observed:** 2 occurrences of the pin (the statement and the
`git show` reconstruction command) and 1 of the parent. **PASS.**

**Group 1 verdict: PASS.** Threat **T-11-09** (repudiation — a wrong sha silently converting named
limit #8 into a false statement) is closed by measurement, not by argument.

---

## 2. The five §F20c decoupling commands

All run against `PHASE11_BASE_SHA = 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348`.

### 2.1 The v1.0 manuscript-pipeline sources are untouched

```console
$ git diff --stat 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348..HEAD -- src/bayes.jl src/colocalization.jl src/LoadImages.jl src/plot.jl
(empty)
```

**Expected:** empty. **Observed:** empty. **PASS.**

### 2.2 Nothing outside `src/amortized/` changed inside `src/`

```console
$ git diff --name-only 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348..HEAD -- src/ | grep -v '^src/amortized/'
(empty)
```

**Expected:** empty. **Observed:** empty (`grep` exit status 1, i.e. no lines matched). **PASS.**

For completeness, the unfiltered list — this is the exact blast radius inside `src/`:

```console
$ git diff --name-only 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348..HEAD -- src/
src/amortized/datagen.jl
src/amortized/ood.jl
src/amortized/simulator.jl
src/amortized/train_npe.jl
```

Four files, all under `src/amortized/`. Note `train_npe.jl` — see §4.2; the plan's Q8 wording
enumerates only three and is corrected there.

### 2.3 The real image bytes are untouched

```console
$ git diff --stat 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348..HEAD -- test/test_images/
(empty)

$ git ls-files -s test/test_images/ | sha256sum
eaee22f9185460fd910788026e7a33511f6de373f313459b96bb43f3fc1ef276 *-
```

**Expected:** empty diff. **Observed:** empty. **PASS.** The digest above is the recorded fingerprint
of the read-only image inputs (index entries: mode, blob sha, stage, path for all six TIFFs).

### 2.4 The shipped artifact pin is unchanged

```console
$ git diff --stat 17ebd1edf6f7f7af51e4ae7f9f0c3faa3fc9d348..HEAD -- Artifacts.toml
(empty)
```

**Expected:** empty. **Observed:** empty. **PASS.** Threat **T-11-11** (tampering with shipped
artifact integrity) has its first line of defence intact: the pin nobody may move has not moved.

### 2.5 The v1 manuscript tests still pass against the edited tree

```console
$ julia --project -t auto -e 'using Pkg; Pkg.test()'
...
Test Summary:                  | Pass  Total  Time
co-resolution gate (Finding 1) |    4      4  2.0s
Test Summary:         | Pass  Total  Time
D-02 result hierarchy |   10     10  0.3s
Test Summary:                                | Pass  Total  Time
LoadImages                                   |   32     32  7.8s
  load_tiff                                  |    3      3  3.8s
  MultiChannelImage constructor              |    6      6  2.6s
  MultiChannelImage constructor (5 channels) |    5      5  0.6s
  calculate and apply mask                   |    7      7  0.4s
  MultiChannelImageStack                     |   11     11  0.4s
Test Summary: | Pass  Total  Time
patch         |    3      3  0.3s
Test Summary:   | Pass  Total  Time
Colocalization  |   29     29  1.0s
...
Test Summary:      | Pass  Total  Time
chromatic ε (D-09) |   25     25  3.3s
...
Test Summary:                                                    | Pass  Total  Time
amended gate machinery (07-GATE-AMENDMENT required code changes) |  254    254  6.8s
     Testing ProteinCoLoc tests passed
EXIT=0
```

**Expected:** green. **Observed:** green, exit 0, zero failures across every testset. **PASS.**

This is the substantive command of the five: `test/runtests.jl` exercises the real-image
`MultiChannelImage` read path (`LoadImages` — 32/32, `load_tiff` — 3/3), so a green suite is direct
evidence the v1 read path is intact after the `src/` edit.

**Group 2 verdict: PASS.** All four diffs empty, suite green.

---

## 3. Shipped-bundle regression (the highest-consequence check)

Run read-only in a single `julia --project` invocation. RESEARCH §A5 argues **from code** that the
8th θ column cannot break the shipped bundle load (`persist.jl:97` reads `a.D` off disk; `θzt` is
frozen data). An argument is not a regression test; this section is the test.

### 3.1 Which resolution branch was taken — recorded, not assumed

```console
[G3.0] Artifacts.toml            : <worktree>/Artifacts.toml
[G3.0] pinned git-tree-sha1      : 90e6b63a8a234d067b407fefd7914f2ae4845448
[G3.0] artifact_exists(expected) : true
[G3.0] isdir(dev bundle dir)     : false  (<worktree>/artifacts/amended_v2/grid_8)
[G3.0] RESOLUTION BRANCH TAKEN   : installed artifact store (Artifacts.artifact_path, hash-verified at install)
```

**The branch taken was branch 1 of `_lazy_load_from_artifact!`'s three (`registry.jl:137-147`): the
already-installed content-addressed artifact store.** Two consequences are stated plainly rather
than glossed:

- **No download happened.** No network access was made and no `ensure_artifact_installed` call ran.
  Nothing here is evidence about the GitHub-Release download path.
- **`_verify_tree_sha1` did NOT run on this branch.** It is called only on branch 2 (the in-repo dev
  fallback). The plan's acceptance criterion says exercising `estimator_for(8)` exercises
  `_verify_tree_sha1`; on the branch actually taken it does not. Branch 1 relies on Pkg's
  install-time verification of the content-addressed store. Section 3.6 supplies the independent
  integrity evidence that this section would otherwise be asserting without proof.

### 3.2 `estimator_for(8)` resolves

```console
[G3.1] estimator_for(8)          : loaded  (16.8 s)
[G3.1] typeof(bundle)            : ProteinCoLoc.EstimatorBundle
[G3.1] bundle.grid               : 8
```

**PASS.** The shipped read surface still resolves after the `src/` edit.

### 3.3 The persisted architecture still reports `D == 7`

```console
[G3.2] npe artifact              : C:\Users\Manuel\.julia\artifacts\90e6b63a…\npe_8.jld2
[G3.2] persisted arch.D          : 7   (expect 7)
[G3.2] persisted arch.d_in       : 128
[G3.2] theta_prior_bounds() len  : 8   (now 8)
[G3.2] ProteinCoLoc.NPE_D        : 7
```

**PASS.** This is the concrete form of RESEARCH §A5's argument: the prior's θ arity is now **8**,
and the shipped net's flow is still a **7**-marginal flow, because `load_estimator` rebuilds from
`arch.D` read off disk (`persist.jl:97`) and never from `theta_prior_bounds()`. `NPE_D` and
`build_estimator`'s default `D` are also still 7, so nothing about the shipped default moved.

Three quantities that merely coincided at 7 are now distinguishable, and the report must keep them
distinct: the **prior's** θ arity (8), the **shipped bundle's** frozen flow marginal count (7, read
off disk), and `NPE_D` / `build_estimator`'s default (7).

### 3.4 The frozen θ transform is 7-element persisted data

```console
[G3.3] typeof(bundle.θzt)        : ProteinCoLoc.BoundedThetaTransform{Float64, StatsBase.ZScoreTransform{Float64, Vector{Float64}}}
[G3.3] length(bundle.θzt.lo)     : 7   (expect 7)
[G3.3] length(bundle.θzt.hi)     : 7   (expect 7)
[G3.3] bundle.θzt.lo             : [-0.99, 0.0, 0.0, 0.6, -1.0, -1.0, 0.0]
[G3.3] bundle.θzt.hi             : [0.99, 0.2, 0.1, 1.0, 1.0, 1.0, 1.0]
```

**PASS.** The `lo`/`hi` vectors are the **frozen** prior box baked in at training time, not a
recomputation — note `lo[5..6] = -1.0` / `hi[5..6] = 1.0`, i.e. the **un-widened** `SHIFT_PRIOR`
(D-02's widening is spike-only, ruling Q2), and no 8th entry despite `theta_prior_bounds()` now
returning 8 rows.

### 3.5 `colocalization_amortized(...)` end to end on the committed real images

Sample `test/test_images/positive/positive_c{1,2,3}.tif`, control
`test/test_images/negative/negative_c{1,2,3}.tif`, both constructed through the frozen
`MultiChannelImage(name, paths, channels)` convenience constructor (`src/LoadImages.jl:432-446`),
channels `[1, 2]`, `num_patches = 8`.

```console
[G3.4] sample  : positive      (1028, 1376)
[G3.4] control : negative  (1028, 1376)
[G3.4] elapsed                   : 8.951 s
[G3.4] typeof(result)            : AmortizedColocResult
[G3.4] result isa AmortizedColocResult : true
[G3.4] n Δρ draws                : 2000
[G3.4] all(isfinite, Δρ draws)   : true
[G3.4] Δρ point estimate (mean)  : 0.11021147163771093
[G3.4] Δρ median                 : 0.10957340896129608
[G3.4] Δρ 90% CI                 : [-0.07727977633476257, 0.30126190185546875]
[G3.4] log Bayes factor          : 1.987931489944458
[G3.4] OOD flag (is_ood)         : true
[G3.4] OOD verdict               : OODVerdict(433.68840347024536, true, (density = 433.68840347024536,))
```

**PASS** on the mechanical criteria: the return type is `AmortizedColocResult` and all 2000 Δρ draws
are finite.

> **THIS IS A SMOKE REGRESSION, NOT A SCIENTIFIC CLAIM.** The point estimate `Δρ ≈ 0.110` and the
> log Bayes factor `≈ 1.99` are recorded **only** so a future run can detect that the shipped read
> path started answering differently. They are **not** evidence of colocalization in these images
> and must not be quoted as such. In particular the OOD flag is **true** (density Mahalanobis
> statistic 433.7, far above the recorded ID operating point): the shipped estimator is telling us
> these real 1028×1376 images sit outside its training summary distribution, which is exactly the
> named limit the v2.0 GO already carries. A flagged result is not a result. What this section
> proves is that the pipeline **runs and returns finite numbers after the `src/` edit** — nothing
> more.

### 3.6 `Artifacts.toml` integrity — one pre-existing finding, fully characterised

Because §3.1 established that `_verify_tree_sha1` did not run on the branch taken, the integrity
property was verified independently. Doing so surfaced a **pre-existing** discrepancy that is
recorded here in full rather than omitted.

```console
$ julia --project -e '... Pkg.GitTools.tree_hash ...'
dev tree_hash   = 90e6b63a8a234d067b407fefd7914f2ae4845448  (== Artifacts.toml pin)
store tree_hash = 7d72b47e8e862f4f328a4e0d9bfaf328b10d94db  (mode-bit only difference)
```

The installed store copy hashes to a **different** tree-sha1 than the `Artifacts.toml` pin. The
bundle bytes are nevertheless identical — proven per file with SHA-256:

| File | store SHA-256 | in-repo dev SHA-256 | identical |
|---|---|---|---|
| `gate_report_8.jld2` | `af3acf7c86db97f6d93e5c66624111c612fdf17fb3293f0f2ff82c00dda2a70f` | same | **yes** |
| `npe_8.jld2` | `198bb0783d35023fc4be6c04740589f6acb7d72fa047f3867fb91358b9042fed` | same | **yes** |
| `ood_nulls_8.jld2` | `5793e7c5f1117b0c728fed354b63a445aad7326542fecc2f541d803c1ed48a8d` | same | **yes** |
| `ratio_8.jld2` | `0d0108e7f6d6adb73342e2c98561c8554d8953a9a2ba079560535a6389438e56` | same | **yes** |

**Root cause, measured.** A git tree hash covers the file **mode** as well as the content, and on
Windows `Pkg.GitTools.gitmode` derives that mode from `Sys.isexecutable`, which tracks the
write bit:

```console
--- C:/Users/Manuel/.julia/artifacts/90e6b63a…            (installed store, read-only)
  gate_report_8.jld2  gitmode=100644  isexec=false  filemode=100444
  npe_8.jld2          gitmode=100644  isexec=false  filemode=100444
  ood_nulls_8.jld2    gitmode=100644  isexec=false  filemode=100444
  ratio_8.jld2        gitmode=100644  isexec=false  filemode=100444
--- .../artifacts/amended_v2/grid_8                       (in-repo dev tree, writable)
  gate_report_8.jld2  gitmode=100755  isexec=true   filemode=100666
  npe_8.jld2          gitmode=100755  isexec=true   filemode=100666
  ood_nulls_8.jld2    gitmode=100755  isexec=true   filemode=100666
  ratio_8.jld2        gitmode=100755  isexec=true   filemode=100666
```

The pin `90e6b63a…` therefore encodes mode `100755`, because it was generated on Windows against a
**writable** working-tree copy. Julia's artifact installation makes store files read-only, flipping
the derived mode to `100644` and hence the tree hash. Running the real guard both ways confirms it:

```console
_verify_tree_sha1(devdir) -> C:/Users/.../artifacts/amended_v2/grid_8        # PASSES
_verify_tree_sha1(store)  -> ERRORED: Artifact integrity check FAILED for the 8×8 estimator bundle:
                             expected git-tree-sha1 90e6b63a8a234d067b407fefd7914f2ae4845448
                             (Artifacts.toml), got 7d72b47e8e862f4f328a4e0d9bfaf328b10d94db …
```

**Scope and disposition.**

- **Not caused by Phase 11.** `Artifacts.toml` is byte-unchanged since `PHASE11_BASE_SHA` (§2.4) and
  the artifact store lives outside the repository entirely. This condition predates every Phase-11
  commit.
- **Not a byte-integrity failure.** All four blobs are byte-identical to the pinned tree. What
  differs is a POSIX permission bit that Windows synthesises.
- **Latent, not currently firing.** On this machine branch 1 is taken and `_verify_tree_sha1` is
  never called, so the mismatch is invisible. It would fire only if the store copy were removed
  *and* the in-repo dev dir were present — and in that configuration the dev dir hashes to the pin,
  so it would **pass**. The failure mode that would actually bite is a non-Windows machine
  reproducing the pin from a writable tree.
- **Deliberately not fixed here.** Any repair moves `Artifacts.toml`, which would break the D-16
  "no artifact change" boundary and the §2.4 byte-unchanged claim this very document rests on.
  Logged to `deferred-items.md` for a future release-engineering pass.

### 3.7 `_theta_tuple(v7)` still feeds `simulate_pair` (risk R3, threat T-11-08)

```console
[G3.5] _theta_tuple(v7) keys     : (:ρ_true, :spillover, :autofluorescence, :label_efficiency, :shift_dx, :shift_dy, :noise, :chromatic_eps)
[G3.5] length(θ7)                : 8
[G3.5] hasproperty chromatic_eps : true
[G3.5] simulate_pair(θ7) accepted: true 2×(128, 128)  all-finite=true
[G3.5] _theta_tuple(v8) keys     : (:ρ_true, :spillover, :autofluorescence, :label_efficiency, :shift_dx, :shift_dy, :noise, :chromatic_eps)
[G3.5] out7 == out8 (ε = 0)      : true
```

**PASS.** A **7**-long posterior-mean vector — exactly what the shipped 7-marginal flow produces —
still produces a θ that `simulate_pair` accepts. `_theta_tuple` emits `chromatic_eps` at its
`length(v) >= 8 ? … : 0.0` guard, so the short vector yields the explicit zero, and the output is
**bit-identical** (`==`, no tolerance) to feeding the 8-long vector with `chromatic_eps = 0.0`. The
posterior-predictive OOD channel and every other legacy `simulate_pair` call site are therefore
unaffected.

**Group 3 verdict: PASS**, with the §3.6 finding recorded and deferred.

---

## 4. The decoupling claim (orchestrator ruling Q8)

### 4.1 The literal claim, quoted

`.planning/ROADMAP.md:126`, verbatim:

> *"2. The main repo and both manuscript pipelines are demonstrably untouched — `git status` on
> `src/`, `bayes.jl`, `colocalization.jl` is clean (decoupling proof)"*

**That claim is no longer available to Phase 11, and this document says so plainly.** D-11
deliberately edits `src/amortized/simulator.jl` in order to keep the spike and productionized
simulators from diverging, and D-12 records the provenance cost of doing so. `git status` on `src/`
is therefore not clean across the phase, by design. Pretending otherwise would be the exact kind of
claim the named-limits mechanism exists to prevent.

### 4.2 The narrower claim actually adopted

> **The v1.0 manuscript pipeline files (`src/bayes.jl`, `src/colocalization.jl`,
> `src/LoadImages.jl`, `src/plot.jl`) and the real image bytes (`test/test_images/`) are
> byte-unchanged; only `src/amortized/` is edited, and within it only `simulator.jl` (the D-11 `ε`
> extension and the D-10 composed warp), plus a defensive read in `ood.jl`, the θ-arity constants in
> `datagen.jl`, and a derived-arity repair in `train_npe.jl`. `Artifacts.toml` and the shipped
> artifact are unchanged.**

Every clause of that sentence is backed by a recorded command in §2. It is honesty-by-naming, in the
same voice as the existing named limits — a narrower claim that is *checkable*, in place of a wider
claim that is *false*.

**One correction to the plan's wording.** The plan (and RESEARCH §F20c) enumerates three edited files
inside `src/amortized/`: `simulator.jl`, `ood.jl`, `datagen.jl`. The measured diff (§2.2) contains a
**fourth**, `src/amortized/train_npe.jl`. That file is in the diff because of plan 11-02's recorded
deviation 1: after the prior gained an 8th row, `fit_theta_transform` demanded 8 θ rows while
`build_estimator` was still constructing a 7-marginal flow from `NPE_D`, so the src-side training
path accepted **no** θ arity at all. The repair derives the flow's marginal count from
`size(θtr, 1)` instead of `NPE_D`. `NPE_D`, `build_estimator`'s default `D` and
`src/amortized/architecture.jl` are untouched (§3.3 measures all three), so the shipped read surface
is unaffected — but the enumeration must name four files, not three, or the claim is inaccurate.

### 4.3 Expected work is not pipeline breakage

The root suite went **red during plan 11-02** on θ-arity literals in `test/runtests.jl`,
`test/gpu_smoke.jl` and `test/test_local_map.jl`. RESEARCH §A4 row 24 predicted exactly this: *"the
suite **will** go red on θ-arity literals until those are updated — that is expected work, not
evidence of pipeline breakage, and the report must distinguish the two."*

The two are distinguished here as follows. 11-02 took a **pre-edit baseline** suite run (green)
before touching anything, updated the arity literals to
`length(ProteinCoLoc.theta_prior_bounds())` rather than to a fresh literal `8`, and reached green
again. §2.5 above re-confirms green at this commit. A green → red → green transition with the red
phase attributable to enumerated test-side literals is *not* the same event as a manuscript pipeline
breaking, and this document does not let the two be conflated.

### 4.4 The `skip_if_done` caveat (RESEARCH §B7)

`_train_grid_pipeline`'s skip branch (`src/amortized/pipeline.jl:163-166`) checks only **artifact
file integrity** (`_estimator_ok` / `_ratio_ok` / `_ood_nulls_ok`), **not** the datagen digest — so a
post-edit `train_and_register(8)` against the existing `artifacts/amended_v2/` root would **load the
old bundle and not retrain**. That is the desired behaviour under D-01 (no reship, no gate
reopened), but it means the digest change is **invisible at that call site**, and nobody should
mistake that silence for "nothing changed". The training-data cache directory has in fact moved for
two independent reasons — the hashed `simulator.jl` source bytes and the hashed
`generating_config.theta_dim` field going 7 → 8 — and the pool the shipped net was trained from is
now orphaned. Orphaned caches are **retained, not pruned**: they are the byte-level fallback if the
git-pinned reconstruction of §1 is ever disputed.

---

## 5. What this document establishes

1. Named limit #8's pinned sha `17ebd1e…` is provably the ε commit's first-parent lineage
   (byte-identical to `a2ced26…` for `simulator.jl`) and provably predates `chromatic_eps`. The
   provenance-by-reference repair is a true statement, not an assumed one.
2. The v1.0 manuscript pipeline sources, the real image bytes and `Artifacts.toml` are byte-unchanged
   since `PHASE11_BASE_SHA`, and the root suite is green against the edited tree.
3. The shipped `amended_v2/grid_8` bundle still loads, still reports a persisted `D == 7` and a
   7-element frozen `BoundedThetaTransform`, and `colocalization_amortized(...)` still returns a
   finite Δρ on the committed real images — so the 8th θ column did not reach the shipped read
   surface. RESEARCH §A5's code argument is now a measurement.
4. A 7-row posterior mean still round-trips through `_theta_tuple` into a `simulate_pair`-valid θ,
   bit-identically to the explicit `chromatic_eps = 0` call.
5. The `Artifacts.toml` tree-sha1 pin does not match the *installed store copy* on this Windows
   machine, for a file-mode reason with byte-identical contents. Pre-existing, latent, characterised
   in §3.6, deferred.
6. The ROADMAP's literal decoupling claim is unavailable to Phase 11 and has been replaced by a
   narrower checkable one, with the plan's three-file enumeration corrected to four.

## 6. What was deliberately NOT done

- **No download path was exercised.** Branch 1 (installed store) was taken; nothing here is evidence
  about the GitHub-Release fetch.
- **`_verify_tree_sha1` was not exercised by `estimator_for(8)`** on the branch taken; it was
  exercised separately (§3.6), and the result is reported as measured rather than as hoped.
- **No scientific claim is made from §3.5.** The Δρ and log-BF numbers are smoke-regression
  fingerprints on an OOD-flagged image pair, not a colocalization result.
- **Nothing was fixed.** The §3.6 finding is recorded and deferred; repairing it would move
  `Artifacts.toml` and break the D-16 boundary this document is auditing.
- **No package was installed** and no dependency was added to either environment (T-11-SC).

---

*Phase: 11-registration-and-chromatic-uncertainty-as-latent*
*Plan: 11-03, task 2*
*Recorded: 2026-07-25*
