---
spike: 006
name: bf-attrition-mechanism
type: standard
validates: "Given the amended BF gate's pair loop, when a pair is dropped, then the non-finite value originates in the non-clamped KDE baseline (p_post at a {0,1} boundary), not in the amortized NRE"
verdict: VALIDATED
related: [007, 008]
tags: [bayes-factor, ship-gate, kde, numerics, attrition]
---

# Spike 006: BF-Arm Attrition Mechanism

## What This Validates

**Given** the amended grid-8 ship-gate's BF pair loop, **when** a pair is discarded by the
`isfinite(a) && isfinite(b)` filter (`test/gate/run_gate.jl:~22`), **then** the non-finite value
comes from the non-clamped KDE baseline reaching an exact probability boundary — not from the
amortized NRE.

Motivation: the confirmatory amended gate obtained only **58 of 100** finite pairs, below the
pre-registered `BF_GATE_N_MIN = 69`, returning `corr_verdict = :invalid`
(`.planning/phases/07-productionization-conditional-on-go/gate-8x8-amended.md` §2).

## Research

No external research: pure numerics over existing code, no new dependencies.

The mechanism is documented in the source itself. `src/amortized/bf.jl:113-124`
(`_p_gt_threshold_unclamped`) returns `clamp(1 - p_le, 0.0, 1.0)` with **no** `1e-8` floor, and the
docstring of `kde_log_bf_unclamped` (`src/amortized/bf.jl:128-138`) states plainly:

> "May return `±Inf` for a perfectly one-sided posterior; that is the honest, un-floored value the
> gate is meant to see."

So the ±Inf were **intended** (T-7-06 / Memo §5 removed the clamp so `max|Δ logBF|` would be free
of the saturation artifact). What was not anticipated is that the gate then *filters those pairs
out*, which starves `n` and — see Spike 007 — biases the survivors.

## How to Run

```bash
SPIKE006_N=100 julia --project=. -t auto .planning/spikes/006-bf-attrition-mechanism/measure_attrition.jl
julia --project=.planning/spikes/figenv .planning/spikes/006-bf-attrition-mechanism/figures.jl
```

Figures use an isolated environment (`.planning/spikes/figenv`) because `CairoMakie` is not a root
dependency and the root co-resolution is fragile (the GLMakie/Makie 0.21-vs-0.24 conflict cost a
Phase-7 plan). The root `Project.toml`/`Manifest.toml` are **not** touched.

## What to Expect

Attrition around 40%, entirely baseline-side, with every dropped pair's `p_post` sitting exactly
on 0.0 or 1.0.

## Investigation Trail

1. **Read the filter.** `bf_gate` keeps a pair only if both log-BF values are finite.
2. **Read the baseline.** `_p_gt_threshold_unclamped` clamps only for numerical-domain safety
   (QuadGK overshoot), which maps roundoff onto an *exact* boundary → `log(p/(1-p))` = ±Inf.
3. **Instrumented replication** (`measure_attrition.jl`) on a DEV seed (`0x00B7A771`, asserted
   disjoint from `PROD_SEED_V2` and all reported/dev seeds), same net, same amended consts.
4. **n=20 smoke** reproduced 40% attrition, 8/8 baseline-side. Scaled to n=100.
5. **Followed a surprise:** 31 of the 60 *kept* pairs also have fully one-sided draws
   (`frac_pos ∈ {0,1}`). One-sidedness is therefore **necessary but not sufficient** — the drop
   additionally requires the QuadGK integral of the KDE to underflow to exactly 1.0 (or 0.0).
   The boundary is a floating-point event, not a clean statistical one.

## Results

**VERDICT: VALIDATED.** n = 100 attempted, DEV seed `0x00B7A771`.

| quantity | value |
|---|---|
| dropped | **40 / 100 (40.0%)** — the confirmatory gate saw 42/100 on a different seed |
| amortized (NRE) non-finite | **0** |
| baseline (KDE) non-finite | **40** |
| both non-finite | 0 |
| dropped with `p_post == 1.0` exactly | 37 |
| dropped with `p_post == 0.0` exactly | 3 |
| dropped with `p_post` strictly inside (0,1) | **0** |

**Every single drop is a boundary event, and every one is baseline-side.** The amortized estimator
never produced a non-finite log-BF.

### Surprise 1 — the drops are not random, they are the high-signal pairs

| group | median \|Δρ_true\| | mean | n |
|---|---|---|---|
| dropped | **0.691** | 0.709 | 40 |
| kept | **0.255** | 0.358 | 60 |

The gate discards precisely the pairs where the two channels differ most — where a Bayes factor is
easiest and most interesting — and computes its correlation on the ambiguous middle. Quantified in
Spike 007.

### Surprise 2 — attrition scales with image size, linking this to finding F5

| image size | n | dropped |
|---|---|---|
| 512² | 42 | 14 (33.3%) |
| 1024² | 30 | 15 (50.0%) |
| 1376×1028 | 21 | 6 (28.6%) |
| 2048² | 7 | **5 (71.4%)** |

Larger images → more pixels per patch → sharper posteriors → more fully one-sided draw sets → more
boundary hits. The move to the realistic image-size mixture (finding F5) therefore *aggravated* an
attrition problem that the 256²-only regime partly masked. The 1376×1028 cell breaks monotonicity
(28.6% on n=21) — sample sizes per cell are small, so read the trend, not the individual cells.

### Figure

`attrition.png` — (A) source of non-finiteness, (B) |Δρ_true| dropped vs kept, (C) `p_post`
histogram showing the dropped mass at exactly 1.0, (D) attrition rate rising with |Δρ_true|.

### What this does NOT establish

- It does not show the *gate's* 42% and this 40% are the same draws — different seeds, same
  mechanism and magnitude.
- It does not by itself prove the survivors give a biased r̂ (that is Spike 007).
- It says nothing about which baseline should replace the current one (Spike 008).
