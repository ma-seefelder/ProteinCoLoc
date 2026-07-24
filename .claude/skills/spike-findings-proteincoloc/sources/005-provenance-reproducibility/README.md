---
spike: 005
name: provenance-reproducibility
type: analysis
validates: "Given the result struct + result.txt/CSV outputs, when audited, then the absence of provenance/reproducibility metadata (and the output-handling hazards) is identified and a run-manifest + provenance-aware result type is proposed"
verdict: VALIDATED
related: [003, 004]
tags: [provenance, reproducibility, output]
---

# Spike 005: Provenance tracking & reproducibility

## What This Validates

Given the analysis outputs (`result.txt`, per-channel CSVs) and the `CoLocResult` struct, when
audited, then the (near-total) absence of reproducibility metadata and the output-handling hazards
are identified and a provenance mechanism is proposed.

## Research / Method

Read `main.jl` (`start_analysis` outputs + guards), `utils.jl` (`generate_txt`, `compute_stats`),
`bayes.jl` (`CoLocResult`, sampling/`vi`/`rand`), `LoadImages.jl` (shuffle RNG). Verified; finder
citations to `colocalization.jl` for the sampling code were **corrected to `bayes.jl`**.

## Findings (all verified)

| ID | Finding | Verdict | Sev | Bug? | Effort |
|----|---------|---------|-----|------|--------|
| BUG-1 | Reported Δρ mixes posterior μ_sample with **prior** μ_control | ✓ confirmed | High | **yes** | S |
| REPRO-1 | No RNG/seed threaded through sampling/shuffling — runs non-reproducible | ⚠ needs-nuance | High | no | M |
| PROV-1 | No run-manifest: zero reproducibility metadata in any output | ✓ confirmed | High | no | M |
| OUT-1 | Hand-rolled `\|`-delimited writer with drift-prone parallel header/value lists | ✓ confirmed | Med | no | M |
| OUT-2 | Ineffective `log.txt` guard + append-mode `result.txt` → silent cross-run contamination | ✓ confirmed | Med | no | S |
| PROV-2 | `CoLocResult` does not carry the parameters that produced it | ⚠ needs-nuance | Med | no | M |

### BUG-1 — Δρ wrong headline effect size — ✓ confirmed (high, BUG)
`utils.jl:368` vs `bayes.jl:110-111`. The Δρ mean/median/95%-CI written to `result.txt` subtract a
**prior** draw of μ_control from a **posterior** draw of μ_sample (two independent distributions),
so the central scientific number is wrong and inconsistent with the Bayes factor in the same row.
Fix: `posterior_samples.:μ_sample .- posterior_samples.:μ_control`. **This is the most important
finding of the whole spike** — independently confirmed by all of perf/bayes/viz/prov. Previously
generated result files are invalid and should be regenerated after the fix.

### REPRO-1 — no seed → non-reproducible — ⚠ needs-nuance (high)
`bayes.jl:308,325,328` (sampling — **not** colocalization.jl), `LoadImages.jl:291,326` (shuffle),
`main.jl:83-106` (no seed param). Every stochastic step (`sample(m,Prior(),N)`, `vi(...)`,
`rand(q,N)`, `shuffle!`, `randperm`) uses the global RNG unseeded, so the posterior, BF, CI, and the
shuffled null control all change run-to-run — a published Bayes factor cannot be reproduced. Fix:
add `seed::Union{Nothing,Integer}=nothing` to `start_analysis`, build one `rng = Xoshiro(seed)` in
`colocalization` and pass it (`sample(rng,…)`, `rand(rng,…)`); for `vi` (no rng arg) call
`Random.seed!(seed)` immediately before; thread `rng` into `shuffle_pixels`/`shuffle_blocks!`.
`seed=nothing` preserves current behavior. Also import `Random.Xoshiro/seed!` (only `shuffle!,
randperm` are currently imported). Pairs with PROV-1/PROV-2. (Note: depends on INFER-1 being fixed
first, since the `vi` call is currently broken.)

### PROV-1 — no run-manifest — ✓ confirmed (high)
`main.jl:108-216`, `utils.jl:353-401`. Nothing records Julia/package versions, seed, timestamp, git
commit, input paths/hashes, or the result-shaping params (`number_patches`, `number_iterations`,
`number_posterior_samples`, `cor_method`, `shuffle`/`shuffle_method`). `MultiChannelImage.path`
exists but is never serialized; the result.txt header is only `channel_1|channel_2|ρ_threshold|bf`
+ stats — a result file is uninterpretable after the fact. Fix: a `write_manifest(output_folder;…)`
called once at the start of `start_analysis`, emitting `run_manifest.toml`/`.json` with
`string(VERSION)`, `Pkg.dependencies()` (or Manifest), seed, `Dates.now()`, git commit if available,
every `start_analysis` arg, and per-input `bytes2hex(SHA.sha256(read(path)))`. Uses stdlib
TOML/SHA/Dates/Pkg. This also gives OUT-2's guard a file the package actually produces. (The seed
field presupposes REPRO-1.)

### OUT-1 — fragile hand-rolled writer — ✓ confirmed (med)
`utils.jl:371-385,388-399,356`. `generate_txt` `join()`s a 24-element header list and a separate
24-element value list with a `|` delimiter via raw `open/write` — the two are maintained
independently and silently drift on any column change. `round.([...];digits=6)` over a heterogeneous
vector promotes integer channel ids to `"1.0"`/`"2.0"` and truncates `bf` precision. The signature
says `file="result.csv"` but callers pass `result.txt` with pipe content. Fix: build one
`DataFrame` row with named/typed columns and `CSV.write(file, row; append=isfile(file))` (CSV/
DataFrames already imported). **Caveat:** this changes the delimiter `|`→`,` (and the GUI may read
`result.txt`); to stay drop-in, keep `.txt` + `delim='|'` while still using the single-schema writer.

### OUT-2 — guard + append contamination — ✓ confirmed (med)
`main.jl:119,171-179`, `utils.jl:371,388`. `start_analysis` errors if `log.txt` exists, but
`log.txt` is written **only by the out-of-scope GUI** (`gui.jl:723`) — for a direct caller the guard
never fires. Meanwhile `result.txt` is opened in append mode (header only when absent), so re-running
into the same folder **silently appends** another run's rows (no run id/separator), while the
per-channel CSVs are **overwritten** — the two artifacts become mutually inconsistent. Fix: guard on
`result.txt` (or on the PROV-1 manifest); open `result.txt` `"w"` once before the channel loop and
append rows **within** the run (append is legitimately needed across channel combinations). Optionally
run-scoped/timestamped filenames tied to the manifest run id.

### PROV-2 — `CoLocResult` lacks its own parameters — ⚠ needs-nuance (med)
`bayes.jl:75-82` (struct), `bayes.jl:310-318,331-335` (construction — **not** colocalization.jl). The
struct stores only `img,control,channels,num_patches,posterior,advi_result`; it omits `cor_method`,
`iter`, `posterior_samples`, so reporting code must re-thread these separately (drift risk). Fix:
add a `provenance::NamedTuple` (or fields) populated at construction from args in scope, with a
defaulted outer constructor so existing positional calls keep working. **Corrections:** there is no
`seed` arg to capture yet (add via REPRO-1 first); `ρ_threshold` belongs to the BF step, not the
run, so don't store it as run provenance. Realistic fields: `cor_method`, `n_iter`,
`n_posterior_samples`. Couples with TYPE-1 (the struct is being parametrized anyway).

## Signal for the Build
**BUG-1 is the top priority of the entire spike** (1-line fix, wrong headline number). Then the
provenance trio in order: REPRO-1 (seed) → PROV-1 (manifest) → PROV-2 (embed in result), since each
builds on the previous. OUT-1/OUT-2 harden the output layer. Provenance + seed directly serve the
project's stated reproducibility goal and the manuscript's defensibility.
