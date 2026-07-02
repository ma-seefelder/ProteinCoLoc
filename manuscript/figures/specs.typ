// figures/specs.typ — figure specification entries (D-10). Each spec reserves a figure's slot
// with its ID, working title, panel list, the phase that produces its data, the claim row(s) it
// supports, and a status placeholder — so later phases swap in the rendered figure without
// renumbering. Every spec names its future figures/<name>.yaml + <name>.typ target so later
// phases clone the F1 pattern (D-14/D-16).
//
// This file is #included by main.typ under the `= Figure Specifications` heading, so its
// top-level content (the placeholder #figure blocks) is emitted here. In Plan 04 Task 1 all
// seven figures (F1–F7) are enumerated as PENDING specs; Task 2 then renders the ONE worked
// exemplar F1 (speedup-vs-RMSE) end-to-end from figures/f1_speedup.{yaml,typ}. F2–F7 stay specs,
// each naming its per-figure YAML target for the producing phase to fill.
//
// Pattern: RESEARCH Pattern 3 — a captioned dashed placeholder carrying the spec text, with a
// unique <fig:...> label. Labels are DEFINED here but not @-referenced (Pitfall 2: an @ref to a
// missing target fails the build; a defined-but-unreferenced label is fine).

// ---------------------------------------------------------------------------
// F1 — Speedup vs RMSE  (Phase 4; figures/f1_speedup.{yaml,typ})
// F1 carries the Phase-4 HONEST caveat (median ~325× forward-pass, ~16× full-workload, ADVI
// ~0.5 s/pair, RMSE parity, wider NPE intervals) — never the bare >100× headline (RESEARCH
// §Honest Headline Numbers). Task 2 renders this spec end-to-end.
// ---------------------------------------------------------------------------
#figure(
  rect(width: 100%, stroke: (dash: "dashed"), inset: 8pt)[
    #set align(left)
    *SPEC F1 — Speedup vs RMSE (amortized NPE vs per-dataset ADVI)* \
    Panels: (a) speedup factor (log) vs ρ-space RMSE, NPE regimes against the ADVI 1× baseline;
    (b) paired 90% interval width, NPE vs ADVI. \
    Data: Phase 4. Supports: C1, C2. \
    Honest caveat: median ~325× forward-pass (min pair ~125×), ~16× on the full posterior-sample
    (N=2000) workload, against a realized ADVI baseline ~0.5 s/pair — reported PAIRED with RMSE
    parity (0.1271 vs 0.1270) and the wider NPE 90% intervals (0.884 vs 0.605). NOT a bare >100×. \
    YAML: figures/f1_speedup.yaml + figures/f1_speedup.typ. Status: PENDING.
  ],
  caption: [Placeholder — the producing phase renders the panel via its per-figure YAML.],
) <fig:f1_speedup>

// ---------------------------------------------------------------------------
// F2 — SBC rank histograms + coverage  (Phase 5; figures/f2_sbc.{yaml,typ})
// ---------------------------------------------------------------------------
#figure(
  rect(width: 100%, stroke: (dash: "dashed"), inset: 8pt)[
    #set align(left)
    *SPEC F2 — SBC rank histograms + coverage* \
    Panels: (a) per-parameter SBC rank histograms with the uniform band; (b) empirical coverage
    vs nominal credible level; (c) KS / χ² rank-uniformity p-values. \
    Data: Phase 5. Supports: C3. \
    Caveat: calibration is UNDER THE SIMULATOR — always reported paired with the OOD result (F3);
    state as "calibrated under the simulator", never unqualified. \
    YAML: figures/f2_sbc.yaml + figures/f2_sbc.typ. Status: PENDING.
  ],
  caption: [Placeholder — Phase 5 renders the SBC panel via its per-figure YAML.],
) <fig:f2_sbc>

// ---------------------------------------------------------------------------
// F3 — OOD ROC over the misspecification grid  (Phase 5; figures/f3_ood.{yaml,typ})
// ---------------------------------------------------------------------------
#figure(
  rect(width: 100%, stroke: (dash: "dashed"), inset: 8pt)[
    #set align(left)
    *SPEC F3 — OOD / misspecification ROC* \
    Panels: (a) ROC of the OOD flag over the misspecification grid; (b) detection rate by
    misspecification type, with the summary-orthogonal BLIND SPOT named explicitly. \
    Data: Phase 5. Supports: C4. \
    Caveat: the flag is structurally bounded by the fixed patch-correlation summary —
    summary-orthogonal misspecifications are provably undetectable (blind spot NAMED, not hidden). \
    YAML: figures/f3_ood.yaml + figures/f3_ood.typ. Status: PENDING.
  ],
  caption: [Placeholder — Phase 5 renders the OOD ROC panel via its per-figure YAML.],
) <fig:f3_ood>

// ---------------------------------------------------------------------------
// F4 — Cross-method comparison  (Phase 9 + Phase-8 corpus; figures/f4_comparator.{yaml,typ})
// ---------------------------------------------------------------------------
#figure(
  rect(width: 100%, stroke: (dash: "dashed"), inset: 8pt)[
    #set align(left)
    *SPEC F4 — Cross-method comparison (vs Costes / Manders / Pearson)* \
    Panels: (a) calibrated disagreement with classical coefficients where they mislead;
    (b) agreement in the well-specified regime; (c) external-corpus validation summary. \
    Data: Phase 9 comparator harness + Phase 8 physical corpus. Supports: C8, C9. \
    Caveat: framed as regime-appropriate, not blanket superiority; corpus is NOT simulator-generated
    (guards the circularity of validating only under the simulator). \
    YAML: figures/f4_comparator.yaml + figures/f4_comparator.typ. Status: PENDING.
  ],
  caption: [Placeholder — Phases 8/9 render the comparator panel via its per-figure YAML.],
) <fig:f4_comparator>

// ---------------------------------------------------------------------------
// F5 — Registration-UQ posterior widening  (Phase 11; figures/f5_registration.{yaml,typ})
// ---------------------------------------------------------------------------
#figure(
  rect(width: 100%, stroke: (dash: "dashed"), inset: 8pt)[
    #set align(left)
    *SPEC F5 — Registration-UQ posterior widening* \
    Panels: (a) posterior width with fixed vs inferred sub-pixel registration; (b) coverage
    restored once registration is promoted to an inferred latent. \
    Data: Phase 11. Supports: C6. \
    Caveat: the posterior widens HONESTLY under registration uncertainty rather than assuming
    fixed registration. \
    YAML: figures/f5_registration.yaml + figures/f5_registration.typ. Status: PENDING.
  ],
  caption: [Placeholder — Phase 11 renders the registration-UQ panel via its per-figure YAML.],
) <fig:f5_registration>

// ---------------------------------------------------------------------------
// F6 — Spatial per-region Δρ map + uncertainty  (Phase 12; figures/f6_spatial.{yaml,typ})
// ---------------------------------------------------------------------------
#figure(
  rect(width: 100%, stroke: (dash: "dashed"), inset: 8pt)[
    #set align(left)
    *SPEC F6 — Spatial per-region Δρ map + uncertainty* \
    Panels: (a) per-region Δρ map over the correlation lattice; (b) matched per-region uncertainty
    map. \
    Data: Phase 12 (descope-to-v2.1 candidate). Supports: C7. \
    Caveat: framed as complementary to Tapqir's sparse single-molecule regime, not "better". \
    YAML: figures/f6_spatial.yaml + figures/f6_spatial.typ. Status: PENDING.
  ],
  caption: [Placeholder — Phase 12 renders the spatial-map panel via its per-figure YAML.],
) <fig:f6_spatial>

// ---------------------------------------------------------------------------
// F7 — Three-hypothesis Bayes factor simplex  (Phase 13; figures/f7_threebf.{yaml,typ})
// ---------------------------------------------------------------------------
#figure(
  rect(width: 100%, stroke: (dash: "dashed"), inset: 8pt)[
    #set align(left)
    *SPEC F7 — Three-hypothesis Bayes factor simplex* \
    Panels: (a) amortized 3-way Bayes-factor simplex (coloc / anti-coloc / null); (b) validation
    against the KDE Bayes-factor baseline in the well-specified regime. \
    Data: Phase 13 (extends the 2-way NRE BF of Phase 5). Supports: C5. \
    Caveat: the NRE log-ratio is the amortized log Bayes factor; validated against the KDE BF
    baseline where well-specified. \
    YAML: figures/f7_threebf.yaml + figures/f7_threebf.typ. Status: PENDING.
  ],
  caption: [Placeholder — Phase 13 renders the three-hypothesis BF panel via its per-figure YAML.],
) <fig:f7_threebf>
