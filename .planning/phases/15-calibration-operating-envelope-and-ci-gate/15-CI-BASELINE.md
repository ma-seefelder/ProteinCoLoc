# Phase 15 — CI Baseline (measured)

**Produced by:** 15-02, Task 2 (`checkpoint:human-verify`, executed autonomously in an unattended run)
**Date:** 2026-08-04
**Workflow:** `.github/workflows/headless-smoke.yml` ("headless load smoke")
**Branch used to trigger:** `gsd/p15-ci-smoke` (throwaway CI branch; never `main`, never `gsd/v2.0-milestone`)
**Repository:** `github.com/ma-seefelder/ProteinCoLoc`

Everything below is copied from real run logs. No number, URL or version string in this file was
inferred, remembered or reconstructed.

---

## 1. Headline verdicts

| Assumption | Statement | VERDICT |
|---|---|---|
| **A4** | `Pkg.instantiate()` on the committed manifest installs the pinned versions without re-resolving | **CONFIRMED** |
| **A5** | GLMakie precompiles and loads on headless `ubuntu-latest` given xvfb + the listed apt packages | **CONFIRMED** |
| *(not an assumption — discovered)* | `estimator_for(8)` resolves the shipped `grid_8` bundle on a runner with no in-repo `artifacts/` tree | **REFUTED — see §5. This is the finding of this plan.** |

**Net effect on the phase:** the headless CI tier is UNBLOCKED (A4 + A5 both hold, and they were the
named highest-risk unknowns). Any CI job that calls `ProteinCoLoc.estimator_for(8)` — which includes
the D-09 fast-tier golden as designed in 15-06 and the slow tier in 15-08 — is **BLOCKED** until a
`v2.0.0` release carrying `grid_8.tar.gz` is published. See §5 for the exact required act.

---

## 2. Runs observed

| # | Cache state | Run URL | Job | Job wall-clock | Overall result |
|---|---|---|---|---|---|
| 1 | **cold** | https://github.com/ma-seefelder/ProteinCoLoc/actions/runs/30935050140 | 92078969512 | **15 m 21 s** (17:41:15Z → 17:56:36Z) | red — failed at the artifact step only |
| 2 | **warm** | https://github.com/ma-seefelder/ProteinCoLoc/actions/runs/30936970798 | 92085434127 | **3 m 07 s** (18:05:30Z → 18:08:37Z) | red — failed at the artifact step only |
| 3 | warm (xvfb control) | https://github.com/ma-seefelder/ProteinCoLoc/actions/runs/30937534659 | 92087338957 | **1 m 41 s** (18:12:47Z → 18:14:28Z) | red — failed at the artifact step only |

Run 3 carried a temporary extra control step (§4) and was reverted afterwards; its timings are
listed because they are the cleanest warm measurement (run 2's `checkout` step took an anomalous
1 m 31 s against 4 s in runs 1 and 3, so run 2's total is inflated by ~1 m 27 s of runner-side
network variance that is not attributable to this workflow).

**Runner:** GitHub-hosted `ubuntu-latest` = **Ubuntu 24.04.4 LTS**, runner version `2.336.0`,
Azure region westus2, 4-core. `julia-actions/setup-julia@v3` reported **`Julia Version 1.12.6`**.

### Per-step timings

| Step | Cold (run 1) | Warm (run 3) |
|---|---|---|
| Set up job | 2 s | 1 s |
| `actions/checkout@v4` | 4 s | 4 s |
| apt install (10 packages) | 20 s | 22 s |
| `julia-actions/setup-julia@v3` | 9 s | 9 s |
| `julia-actions/cache@v3` (restore) | 1 s | 8 s |
| `Pkg.instantiate()` | **11 m 41 s** | **3 s** |
| load under `xvfb-run` | **2 m 28 s** | **14 s** |
| version assertion | 1 s | 1 s |
| `estimator_for(8)` | 20 s (fail) | 22 s (fail) |
| cache save (post) | 11 s | 9 s |

**Cache:** `julia-actions/cache@v3` saved and restored **~537 MB (563,274,657 B)** at ~75 MB/s, under
key `julia-cache;workflow=headless load smoke;job=load;os=Linux;run_id=…;run_attempt=1`.

### Minute budget the two later workflows can plan against

- **Fast tier (D-08/D-09) fixed overhead: ≈ 1 m on a warm cache, ≈ 15 m on a cold one.**
  That overhead is checkout + apt + setup-julia + cache restore + instantiate + one package load.
  Whatever the golden script itself costs is *added* to that.
- The dominant cold cost is `Pkg.instantiate()` (11 m 41 s), reduced to 3 s warm — so the cache is
  worth ~230× on that step and must not be disabled.
- Package load is 2 m 28 s cold / 14 s warm; `cache-compiled` (default `true`) is what buys this.
- Cold cost recurs whenever `Project.toml`/`Manifest.toml` change, and whenever the 7-day Actions
  cache eviction window lapses on a quiet branch. Budget the cold number, not the warm one, for a
  worst-case slow-tier estimate.
- The slow tier's `JULIA_NUM_THREADS` must be set explicitly; `ubuntu-latest` is 4-core, which is
  **8×** below the 32 threads the ~8.5 h sweep budget assumes. A full sweep does not fit a
  GitHub-hosted runner inside the 6 h job cap.

---

## 3. The apt package list that worked

Installed **before any Julia step**, exactly as committed in
`.github/workflows/headless-smoke.yml`:

```
xorg-dev  mesa-utils  xvfb  libgl1  freeglut3-dev
libxrandr-dev  libxinerama-dev  libxcursor-dev  libxi-dev  libxext-dev
```

**Iterations required: ZERO.** This list worked on the first attempt. `apt-get` reported
`2 upgraded, 87 newly installed, 0 to remove and 81 not upgraded.` with no unmet-dependency error
and no `Unable to locate package` on Ubuntu 24.04.4. There is no second attempt to record.

**Sufficiency vs necessity — stated precisely.** This list is proven **sufficient**. It is
**NOT proven minimal**: no run removed a package to see what breaks, so any claim that all ten are
required would be unmeasured. Keeping all ten is nonetheless the right call — it is the list
Makie's own headless CI uses, it costs ~20 s, and narrowing it on a single observation would trade
a known-good recipe for an untested one.

---

## 4. Is `xvfb-run` required? — measured, not assumed

**YES, for every step that loads `ProteinCoLoc`. NO, for steps that do not.**

This was tested with an explicit control rather than inferred. Run 3 (30937534659) carried a
temporary `continue-on-error` step running the identical load *without* the prefix. It failed:

```
ERROR: InitError: Exception[GLFW.GLFWError(65550, "X11: The DISPLAY environment variable is missing"),
                            ErrorException("glfwInit failed")]
  [1] __init__()
    @ GLFW ~/.julia/packages/GLFW/wA4ue/src/GLFW.jl:69
during initialization of module GLFW
```

So the requirement is real and it bites at GLFW's `__init__`, i.e. at *module initialisation*, which
is reached by any `using ProteinCoLoc` because `src/ProteinCoLoc.jl:28` imports GLMakie
unconditionally. The control step was removed once measured; the finding is recorded as a comment
beside the steps it constrains.

Per-step requirement, as observed:

| Step | Prefixed? | Outcome |
|---|---|---|
| `Pkg.instantiate()` | no | success — Pkg does not load the package |
| `using ProteinCoLoc` | **yes** | success; **fails without the prefix** (control above) |
| version assertion via `Pkg.dependencies()` | no | success — does not load the package |
| `estimator_for(8)` | **yes** | reached the artifact fetch, so the prefix was not the failure cause |

**Trap worth recording for whoever writes the next workflow:** a `continue-on-error: true` step
reports `conclusion: success` in the jobs API even when it actually failed. Read the step's *log*,
not its conclusion, or the control silently inverts its own answer.

---

## 5. THE FINDING — the shipped `grid_8` artifact is not fetchable

Both runs failed at the same step, and the cause is not a workflow defect:

```
Downloading artifact: grid_8
     Failure artifact: grid_8
ERROR: Unable to automatically download/install artifact 'grid_8' from sources listed in
       '/home/runner/work/ProteinCoLoc/ProteinCoLoc/Artifacts.toml'.
Sources attempted:
- https://pkg.julialang.org/artifact/90e6b63a8a234d067b407fefd7914f2ae4845448
    Error: RequestError: HTTP/2 404
- https://github.com/ma-seefelder/ProteinCoLoc/releases/download/v2.0.0/grid_8.tar.gz
    Error: RequestError: HTTP/2 404
```

Confirmed independently of the runner:

- `curl -L` on the release asset URL → **HTTP 404**
- `curl -L` on the package-server mirror → **HTTP 404**
- `gh release list` returns only `v1.0.1` (2024-07-10) and `v1.0.0.compiled` (2024-01-25).
  **There is no `v2.0.0` release at all.**

`Artifacts.toml` declares `grid_8` with `git-tree-sha1 = 90e6b63a8a234d067b407fefd7914f2ae4845448`,
`sha256 = 17162904c69b30e2a99a6dd934434553d08ebaadbfd55dea664b36efcf92b6d5`, `lazy = true`, pointing
at that non-existent asset. The files themselves exist only in the author's working checkout at
`artifacts/amended_v2/grid_8/` (`npe_8.jld2`, `ratio_8.jld2`, `ood_nulls_8.jld2`,
`gate_report_8.jld2`), which is git-ignored.

**What this means, in plain terms.** Nobody who clones this repository can obtain the shipped
network. `estimator_for(8)` works on the author's machine only because the git-ignored dev directory
happens to be there. The v2.0 headline claim — amortized inference from the shipped net — is not
reproducible by a third party today, and the CI runner is simply the first machine honest enough to
say so. The 15-RESEARCH note that "the release tarball contains the amended_v2 net" was verified
against a *local* tarball, not against the published release.

**Why the step was left failing rather than softened.** Making it non-blocking, or deleting it,
would turn a true statement about a broken shipping path into a green badge. The red run is the
accurate report.

**Required next step (NOT taken here — out of this plan's scope, and a public, irreversible act
that needs the author's decision):** publish a `v2.0.0` release with `grid_8.tar.gz` built from
`artifacts/amended_v2/grid_8/`, and verify the uploaded tarball's `sha256` equals
`17162904c69b30e2a99a6dd934434553d08ebaadbfd55dea664b36efcf92b6d5` and its unpacked tree-sha1 equals
`90e6b63a8a234d067b407fefd7914f2ae4845448` — otherwise `ensure_artifact_installed` will reject it and
the step will still be red, for a different reason.

**Downstream gating:** 15-06 (fast-tier golden) and 15-08 (slow tier) both depend on
`estimator_for(8)`. Neither can go green until the release exists. They can still be *written*
against this proven environment recipe; they cannot be *proven* against it.

---

## 6. Hardening items deferred to pre-v2.0-release

| Item | Status | Note |
|---|---|---|
| **SHA-pin the third-party actions** | **DEFERRED — do before v2.0 release** | `actions/checkout@v4`, `julia-actions/setup-julia@v3` and `julia-actions/cache@v3` are pinned to *major tags*, which are mutable refs: the owner can repoint them at new code that runs with this workflow's token. Replace each with a full 40-char commit SHA plus a `# vX.Y.Z` comment. Recorded here rather than silently skipped (threat T-15-06). |
| `actions/checkout@v4` Node 20 deprecation | observed, no action needed yet | The runner annotates: "Node.js 20 is deprecated … forced to run on Node.js 24: actions/checkout@v4". Resolve naturally when SHA-pinning to a v5 release. |
| Windows fast-tier runner | deferred (already in 15-CONTEXT Deferred) | The Windows-tauglich constraint is still asserted in prose only. |

---

## 7. Exactness check

The apt list recorded in §3 is the list in `.github/workflows/headless-smoke.yml` at this branch's
HEAD, character-for-character. The workflow was edited three times across the runs above:

1. `c5ff0f3` — original file → run 1 (cold).
2. `02689e6` — header comment recording the artifact finding → run 2 (warm). No step changed.
3. `738684b` — added the temporary xvfb control step → run 3.
4. (after run 3) — control step removed, its measured result kept as a comment. **No step that
   executed in any run above was altered by this commit**, so every timing and every verdict here
   still describes the file at HEAD.
