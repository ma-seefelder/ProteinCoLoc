// claims.typ — the claim spine as DATA + a render function (D-06, RESEARCH Pattern 1).
//
// `#let`-only (no bare top-level content) so `#import "claims.typ": claims, claim_table`
// pulls definitions and emits NOTHING at the import site (Pitfall 4). The rendered table is
// produced by main.typ via `#figure(claim_table(claims), ...) <tab:claims>`.
//
// The `claims` array is itself the diffable plain-source artifact; `claim_table` renders it.
// Schema per row (all string-valued for a clean table + easy diff):
//   id     — claim identifier (e.g. "C1")
//   claim  — one-line claim statement (differentiator headline)
//   axis   — differentiator axis (Amortization | Calibration/SBC | Speed | ...)
//   phase  — supporting / owning phase(s) (e.g. "4,7") — D-07 artifact linkage (owning phase)
//   fig    — figure/experiment spec ID that substantiates it (e.g. "F1") — D-07 artifact linkage
//   status — "supported (measured)" | "supported (caveat)" | "reserved slot" (honesty posture)
//   notes  — the honest caveat carried VERBATIM (Phase-4 numbers; "under the simulator"; OOD
//            blind spot) so the honesty slots exist by construction and cannot be quietly
//            dropped by later phases (RESEARCH §"Honest Headline Numbers"; Pitfall 5).
//
// This is the FULL seeded headline claim set (Plan 10-02) mapped from CONTEXT D-06 and the
// RESEARCH §"Phase Requirements → claim→phase→status" table. Each row names the artifact
// (fig spec ID + owning/supporting phase) that will substantiate it (D-07). The status column
// separates supported-measured from reserved-slot so the skeleton never overclaims (Pitfall 5).

#let claims = (
  (
    id: "C1",
    claim: "Amortized ms-scale colocalization inference in a single forward pass",
    axis: "Amortization",
    phase: "4,7",
    fig: "F1",
    status: "supported (measured)",
    notes: "NPE trained once; new datasets inferred in one amortized forward pass (ms). "
      + "The cleanest, least-contestable delta vs per-dataset ADVI/SVI.",
  ),
  (
    id: "C2",
    claim: ">100x speedup at comparable RMSE (caveated — not a bare >100x vs minutes)",
    axis: "Speed",
    phase: "4",
    fig: "F1",
    status: "supported (caveat)",
    notes: "Median ~325x forward-pass (min pair ~125x); ~16x on the full posterior-sample "
      + "(N=2000) workload; against a realized ADVI baseline ~0.5 s/pair (seconds, not minutes). "
      + "Only meaningful PAIRED with RMSE parity (NPE 0.1271 vs ADVI 0.1270) and reported with "
      + "the wider NPE 90% intervals (0.884 vs ADVI 0.605 — mean-field VI under-disperses).",
  ),
  (
    id: "C3",
    claim: "Formal SBC rank-uniformity + coverage calibration proof",
    axis: "Calibration/SBC",
    phase: "5",
    fig: "F2",
    status: "reserved slot",
    notes: "Calibration is UNDER THE SIMULATOR (SBC pending Phase 5); always reported PAIRED "
      + "with the OOD result. State as 'calibrated under the simulator', never 'calibrated' "
      + "unqualified.",
  ),
  (
    id: "C4",
    claim: "Honest OOD / misspecification flag",
    axis: "OOD",
    phase: "5",
    fig: "F3",
    status: "reserved slot; blind-spot named",
    notes: "The flag is structurally bounded by the fixed patch-correlation summary: "
      + "summary-orthogonal misspecifications are provably undetectable — this OOD blind spot "
      + "is NAMED, not hidden. Never claim 'detects all misspecification'.",
  ),
  (
    id: "C5",
    claim: "Amortized Bayes factor for model comparison (2-way to 3-way)",
    axis: "Model comparison",
    phase: "5,13",
    fig: "F7",
    status: "reserved slot",
    notes: "NRE log-ratio as the amortized log Bayes factor: 2-way coloc-vs-null (Phase 5), "
      + "extended to 3-way (Phase 13). Validated against the KDE BF baseline in the "
      + "well-specified regime.",
  ),
  (
    id: "C6",
    claim: "Registration uncertainty promoted to an inferred latent",
    axis: "Registration-UQ",
    phase: "11",
    fig: "F5",
    status: "reserved slot",
    notes: "Sub-pixel registration inferred as a latent so the posterior widens honestly under "
      + "registration uncertainty, rather than assuming fixed registration.",
  ),
  (
    id: "C7",
    claim: "Amortized spatial per-region colocalization map with per-region UQ",
    axis: "Spatial map",
    phase: "12",
    fig: "F6",
    status: "reserved slot",
    notes: "Per-region Delta-rho map over a dense correlation lattice (Phase 12; descope-to-v2.1 "
      + "candidate). Framed as complementary to Tapqir's sparse single-molecule regime, not "
      + "'better'.",
  ),
  (
    id: "C8",
    claim: "Cross-method: knows when the classics (Costes/Manders/Pearson) are wrong",
    axis: "Comparator",
    phase: "9",
    fig: "F4",
    status: "reserved slot",
    notes: "Comparator harness (Phase 9) demonstrates calibrated disagreement with classical "
      + "coefficients where they mislead. Framed as regime-appropriate, not blanket superiority.",
  ),
  (
    id: "C9",
    claim: "Non-circular external validation on a physical corpus",
    axis: "External validation",
    phase: "8",
    fig: "F4",
    status: "reserved slot",
    notes: "Validation on real images (Phase 8) that are NOT generated by the training simulator, "
      + "guarding against the circularity of only validating under the simulator (pairs with C3).",
  ),
)

// Render the claim array as a Typst table. `fr` columns distribute available width so the long
// claim/notes text wraps instead of overflowing the page. Column count (7) matches the header.
#let claim_table(claims) = table(
  columns: (auto, 2fr, auto, auto, auto, auto, 3fr),
  align: (left + top),
  table.header([ID], [Claim], [Axis], [Phase], [Fig], [Status], [Notes / honest caveat]),
  ..claims
    .map(c => (c.id, c.claim, c.axis, c.phase, c.fig, c.status, c.notes))
    .flatten(),
)
