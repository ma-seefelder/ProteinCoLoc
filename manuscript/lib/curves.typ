// curves.typ — manuscript-owned discrimination + reliability curve module.
//
// Promotion of two validated spike prototypes into one manuscript lib module:
//   * Spike 004.1 `curve_templates.typ` — ROC / Precision–Recall family.
//   * Spike 004.2 `calibration_curve.typ` — reliability diagram + D-11 honesty
//     contract (rescaled comparator scores †-tagged).
//
// This is the single source of truth for every curve figure in the manuscript
// (Phase-4 synthetic head-to-head 04-07, ablation 04-08, Phase-5 cross-study).
// It is re-exported through `lib/template.typ`, so every section and figure
// reaches `roc-curve` / `pr-curve` / `roc-pr-duo` / `reliability-curve` via the
// one canonical `#import "../lib/template.typ": *`. The lilaq 0.6.0 toolchain is
// locked here (D-13, MET-04); no figure may fork to an external raster plotter.
//
// MODULE CONTRACT (root-independent, caller passes the palette in):
//   The module never calls yaml() itself. The caller loads /colours.yaml once
//   and passes `palette` into every plotting function, so the module is
//   root-independent and unit-testable.
//
//   #import "../lib/template.typ": *
//   #let palette = yaml("/colours.yaml")
//
//   // ROC / PR (model variants → palette.model_series, or per-series color:):
//   #let series = series-from-csv(
//     palette,
//     roc-csv: csv("benchmark_synthetic_data/roc.csv"),
//     pr-csv:  csv("benchmark_synthetic_data/pr.csv"),
//     metrics-csv: csv("benchmark_synthetic_data/metrics.csv"),
//     models: (
//       (key: "saintexpress",       label: [SAINTexpress],       color: rgb(palette.tool_series.saintexpress)),
//       (key: "bayesinteractomics", label: [BayesInteractomics], color: rgb(palette.tool_series.bayesinteractomics)),
//     ),
//   )
//   #roc-pr-duo(series, palette: palette, pos-rate: 0.3464)
//   #curve-legend(series, palette: palette)
//
//   // Reliability (tools → palette.tool_series):
//   #let calib = series-from-calib(
//     palette,
//     bins-csv:    csv("benchmark_synthetic_data/calibration_bins.csv"),
//     metrics-csv: csv("benchmark_synthetic_data/calibration_metrics.csv"),
//     tools: (
//       (key: "saintexpress",       label: [SAINTexpress]),
//       (key: "bayesinteractomics", label: [BayesInteractomics]),
//     ),
//   )
//   #reliability-panel(calib, palette: palette)
//
// Compile callers with `typst compile --root . …` so `/colours.yaml` resolves.

#import "@preview/lilaq:0.6.0" as lq

// ===========================================================================
// SHARED — number formatting. Pad to exactly two decimals so legends read
// evenly. Identical in both spike modules; kept as a single definition here.
// ===========================================================================

#let fmt2(x) = {
  let s = str(calc.round(x, digits: 2))
  if s.contains(".") {
    let p = s.split(".")
    if p.at(1).len() < 2 { s = p.at(0) + "." + p.at(1) + "0" }
  } else { s = s + ".00" }
  s
}

// ===========================================================================
// SHARED — legend swatch. Identical in both spike modules; kept once.
// ===========================================================================

#let _swatch(color, thickness) = box(
  width: 2.2em, height: thickness, fill: color, baseline: -0.3em,
)

// ###########################################################################
// ROC / PRECISION–RECALL  (promoted from Spike 004.1 curve_templates.typ)
// ###########################################################################

// ===========================================================================
// LOCKED STYLE — the ROC/PR visual contract. 7×6.5cm single-panel.
// ===========================================================================

#let curve-style = (
  panel:        (width: 7cm, height: 6.5cm),  // canonical single-panel size
  stroke-main:  1.2pt,                         // calibrated / model-of-interest
  stroke-muted: 0.9pt,                         // jaggy / uncalibrated baseline
  stroke-ref:   0.7pt,                         // chance reference line
  ref-dash:     "dashed",
  // Per-model-series stroke weight. dnn_only is downweighted because the
  // uncalibrated prior is jaggy at low recall and otherwise dominates the eye
  // (finding carried over from Spike 001.1 v2).
  stroke-for: (
    dnn_only:         0.9pt,
    meta_learner:     1.2pt,
    meta_learner_unc: 1.2pt,
  ),
)

// Map a model-series key → its semantic colour from /colours.yaml.
// Falls back to axis_chrome for any key not in the model_series block.
// NOTE: the four benchmark TOOLS live under palette.tool_series, not
// model_series (RESEARCH Pitfall 6). Callers plotting tools pass an explicit
// per-series `color: rgb(palette.tool_series.<key>)` via series-from-csv,
// which honours `m.at("color", default: ...)`; this default lookup is only the
// fallback for BayesInteractomics' OWN internal model variants (model_series).
#let series-color(key, palette) = {
  let ms = palette.model_series
  if key in ms { rgb(ms.at(key)) } else { rgb(palette.reference.axis_chrome) }
}

#let series-stroke(key) = {
  let sf = curve-style.stroke-for
  if key in sf { sf.at(key) } else { curve-style.stroke-main }
}

// ===========================================================================
// Data loading — pull one model's (x, y) out of a long-format curve CSV in the
// canonical schema shared by Spike 001.1 (real) and 004.1 (synthetic):
//   ROC: row(model), fpr, tpr, threshold
//   PR : row(model), precision, recall, threshold   (plot recall on x)
// Boundary rows carry Inf in the threshold column but valid corner coords; we
// keep them so curves anchor at (0,0)/(1,1).
// ===========================================================================

#let _xy(raw, name, x-idx, y-idx) = {
  let r = raw.slice(1).filter(row => row.at(0) == name)
  (r.map(row => float(row.at(x-idx))), r.map(row => float(row.at(y-idx))))
}

#let _metric(raw, name, col) = {
  let hit = raw.slice(1).filter(row => row.at(0) == name)
  if hit.len() == 0 { none } else { float(hit.at(0).at(col)) }
}

// Build the `series` array consumed by every plotting function from the three
// CSVs + a model spec list. Each entry is self-describing (colour, stroke,
// curves, metrics) so the plotting functions stay dumb and uniform.
#let series-from-csv(
  palette,
  roc-csv: none, pr-csv: none, metrics-csv: none,
  models: (),
) = {
  models.map(m => {
    let key = m.key
    let (fpr, tpr) = if roc-csv == none { ((), ()) } else { _xy(roc-csv, key, 1, 2) }
    let (rec, prec) = if pr-csv == none { ((), ()) } else { _xy(pr-csv, key, 2, 1) }
    (
      key:     key,
      label:   m.label,
      color:   m.at("color", default: series-color(key, palette)),
      stroke:  m.at("stroke", default: series-stroke(key)),
      roc:     (fpr, tpr),
      pr:      (rec, prec),
      auc-roc: if metrics-csv == none { none } else { _metric(metrics-csv, key, 3) },
      auc-pr:  if metrics-csv == none { none } else { _metric(metrics-csv, key, 4) },
    )
  })
}

// ===========================================================================
// ROC panel — diagonal chance line + one stroke per series.
// ===========================================================================

#let roc-curve(
  series,
  palette: none,
  width: curve-style.panel.width,
  height: curve-style.panel.height,
  title: [ROC],
) = {
  let c-ref = rgb(palette.reference.chance_line)
  let plots = series.map(s => lq.plot(
    s.roc.at(0), s.roc.at(1),
    stroke: (paint: s.color, thickness: s.stroke), mark: none, label: none,
  ))
  lq.diagram(
    width: width, height: height,
    xlim: (0, 1), ylim: (0, 1),
    xlabel: [False positive rate], ylabel: [True positive rate],
    title: title, legend: none,
    lq.plot((0, 1), (0, 1),
      stroke: (paint: c-ref, thickness: curve-style.stroke-ref, dash: curve-style.ref-dash),
      mark: none, label: none),
    ..plots,
  )
}

// ===========================================================================
// PR panel — horizontal chance line at the class prior + one stroke per series.
// ===========================================================================

#let pr-curve(
  series,
  palette: none,
  pos-rate: 0.5,
  width: curve-style.panel.width,
  height: curve-style.panel.height,
  title: [Precision–Recall],
) = {
  let c-ref = rgb(palette.reference.chance_line)
  let plots = series.map(s => lq.plot(
    s.pr.at(0), s.pr.at(1),
    stroke: (paint: s.color, thickness: s.stroke), mark: none, label: none,
  ))
  lq.diagram(
    width: width, height: height,
    xlim: (0, 1), ylim: (0, 1),
    xlabel: [Recall], ylabel: [Precision],
    title: title, legend: none,
    lq.plot((0, 1), (pos-rate, pos-rate),
      stroke: (paint: c-ref, thickness: curve-style.stroke-ref, dash: curve-style.ref-dash),
      mark: none, label: none),
    ..plots,
  )
}

// ===========================================================================
// Shared ROC/PR legend strip — swatch | name | AUC-ROC · AUC-PR. Sits between
// the panel grid and the caption (chance references are explained in the
// caption, never in the legend — Spike 001.1 v2 convention).
// ===========================================================================

#let curve-legend(series, palette: none, show-metrics: true) = {
  let chrome = rgb(palette.reference.axis_chrome)
  let cells = series.map(s => {
    let metric = if show-metrics and s.auc-roc != none and s.auc-pr != none {
      text(fill: chrome, size: 0.92em)[AUC-ROC #fmt2(s.auc-roc) · AUC-PR #fmt2(s.auc-pr)]
    } else { [] }
    (_swatch(s.color, s.stroke), s.label, metric)
  }).flatten()
  align(center)[
    #grid(
      columns: (auto, auto, auto),
      column-gutter: 1.2em, row-gutter: 0.55em,
      align: (center + horizon, left + horizon, left + horizon),
      ..cells,
    )
  ]
}

// ===========================================================================
// Convenience: the standard two-panel ROC + PR block used by main figures.
// One call renders the locked-layout duo (ROC left, PR right).
// ===========================================================================

#let roc-pr-duo(series, palette: none, pos-rate: 0.5, gutter: 8mm) = align(center)[
  #grid(
    columns: (auto, auto), column-gutter: gutter,
    roc-curve(series, palette: palette),
    pr-curve(series, palette: palette, pos-rate: pos-rate),
  )
]

// ###########################################################################
// CALIBRATION / RELIABILITY  (promoted from Spike 004.2 calibration_curve.typ)
//
// HONESTY CONTRACT (D-11, locked with the user 2026-05-29):
//   posterior_prob / SaintScore / MiST score are kept native on [0,1]; only
//   CompPASS scoreWD is min-max rescaled (flagged `rescaled` per series). The
//   caller's caption must state that a rescaled comparator x-axis is a score,
//   not a probability — only BayesInteractomics' diagonal is literally a
//   calibration statement. The legend tags rescaled tools with a † marker.
// ###########################################################################

// ===========================================================================
// LOCKED STYLE — the reliability visual contract. SQUARE 7×7cm ([0,1]²).
// Distinct from curve-style (ROC/PR is 7×6.5cm) — different visual contracts.
// ===========================================================================

#let calib-style = (
  panel:        (width: 7cm, height: 7cm),  // SQUARE — reliability is [0,1]²
  stroke-main:  1.1pt,                        // comparator series
  stroke-hero:  1.5pt,                        // BayesInteractomics (the posterior)
  stroke-ref:   0.7pt,                        // diagonal perfect-calibration line
  ref-dash:     "dashed",
  mark:         "o",
  mark-size:    4pt,
  // Per-tool stroke weight — the hero posterior is drawn boldest; comparators
  // recede. Any key not listed falls back to stroke-main.
  stroke-for: (
    bayesinteractomics: 1.5pt,
    saintexpress:       1.1pt,
    comppass:           1.1pt,
    mist:               1.1pt,
  ),
)

// Map a tool key → its semantic colour from /colours.yaml `tool_series`.
// Falls back to axis_chrome for any key outside the tool_series block.
// NAMESPACED (calib-) so the tool_series lookup never collides with the
// ROC/PR series-color which reads model_series.
#let calib-series-color(key, palette) = {
  let ts = palette.tool_series
  if key in ts { rgb(ts.at(key)) } else { rgb(palette.reference.axis_chrome) }
}

// NAMESPACED (calib-) so the calibration stroke table never collides with the
// ROC/PR series-stroke which reads curve-style.stroke-for.
#let calib-series-stroke(key) = {
  let sf = calib-style.stroke-for
  if key in sf { sf.at(key) } else { calib-style.stroke-main }
}

// ===========================================================================
// Data loading — pull one tool's reliability points + summary metrics out of
// the long-format Spike-004.2 CSVs:
//   bins   : row(tool), p_pred, p_obs, bin_count, bin_index   (plot p_pred vs p_obs)
//   metrics: row(tool), n, n_pos, pos_rate, ece, mce, brier, rescaled
// Bins arrive already sorted by bin_index from calibration.jl; we resort by
// p_pred defensively so the connecting line is always left-to-right.
// ===========================================================================

#let _calib-bins-for(raw, name) = {
  let r = raw.slice(1).filter(row => row.at(0) == name)
  let trip = r.map(row => (float(row.at(1)), float(row.at(2)), int(row.at(3))))
  let srt = trip.sorted(key: t => t.at(0))
  (srt.map(t => t.at(0)), srt.map(t => t.at(1)), srt.map(t => t.at(2)))
}

// NAMESPACED (calib-) — returns the RAW cell (string), unlike the ROC/PR
// _metric which casts to float. series-from-calib casts selectively because
// column 7 (rescaled) is a boolean-as-string, not a number.
#let _calib-metric(raw, name, col) = {
  let hit = raw.slice(1).filter(row => row.at(0) == name)
  if hit.len() == 0 { none } else { hit.at(0).at(col) }
}

// Build the `series` array consumed by every plotting function. Each entry is
// self-describing (colour, stroke, points, metrics, rescaled flag) so the
// plotting functions stay dumb and uniform — same idiom as the ROC/PR side.
#let series-from-calib(
  palette,
  bins-csv: none, metrics-csv: none,
  tools: (),
) = {
  tools.map(t => {
    let key = t.key
    let (xs, ys, counts) = if bins-csv == none { ((), (), ()) } else { _calib-bins-for(bins-csv, key) }
    let m(c) = if metrics-csv == none { none } else {
      let v = _calib-metric(metrics-csv, key, c)
      if v == none { none } else { float(v) }
    }
    let rescaled = if metrics-csv == none { false } else {
      _calib-metric(metrics-csv, key, 7) == "true"
    }
    (
      key:      key,
      label:    t.label,
      color:    t.at("color", default: calib-series-color(key, palette)),
      stroke:   t.at("stroke", default: calib-series-stroke(key)),
      pts:      (xs, ys, counts),
      ece:      m(4),
      mce:      m(5),
      brier:    m(6),
      rescaled: rescaled,
    )
  })
}

// ===========================================================================
// Reliability panel — dashed perfect-calibration diagonal + one line-with-
// markers per tool. Points above the diagonal = under-confident, below =
// over-confident.
// ===========================================================================

#let reliability-curve(
  series,
  palette: none,
  width: calib-style.panel.width,
  height: calib-style.panel.height,
  title: [Reliability],
  xlabel: [Predicted probability (rescaled score†)],
  ylabel: [Observed interactor fraction],
) = {
  let c-ref = rgb(palette.reference.chance_line)
  let plots = series.map(s => lq.plot(
    s.pts.at(0), s.pts.at(1),
    stroke: (paint: s.color, thickness: s.stroke),
    mark: calib-style.mark, mark-size: calib-style.mark-size, color: s.color,
    label: none,
  ))
  lq.diagram(
    width: width, height: height,
    xlim: (0, 1), ylim: (0, 1),
    xlabel: xlabel, ylabel: ylabel,
    title: title, legend: none,
    lq.plot((0, 1), (0, 1),
      stroke: (paint: c-ref, thickness: calib-style.stroke-ref, dash: calib-style.ref-dash),
      mark: none, label: none),
    ..plots,
  )
}

// ===========================================================================
// Shared reliability legend strip — swatch | name(†) | ECE · MCE · Brier. The
// † marks a tool whose x-axis is a rescaled score, not a probability (its
// metrics are indicative only). Perfect-calibration diagonal is explained in
// the caption, never the legend — same convention as the ROC/PR side.
// ===========================================================================

#let calib-legend(series, palette: none, show-metrics: true) = {
  let chrome = rgb(palette.reference.axis_chrome)
  let cells = series.map(s => {
    let name = if s.rescaled { [#s.label#super[†]] } else { s.label }
    let metric = if show-metrics and s.ece != none and s.mce != none and s.brier != none {
      text(fill: chrome, size: 0.92em)[ECE #fmt2(s.ece) · MCE #fmt2(s.mce) · Brier #fmt2(s.brier)]
    } else { [] }
    (_swatch(s.color, s.stroke), name, metric)
  }).flatten()
  align(center)[
    #grid(
      columns: (auto, auto, auto),
      column-gutter: 1.2em, row-gutter: 0.55em,
      align: (center + horizon, left + horizon, left + horizon),
      ..cells,
    )
  ]
}

// ===========================================================================
// Convenience: the standard single-panel reliability block + legend in one
// call (panel centred above the legend strip).
// ===========================================================================

#let reliability-panel(series, palette: none, title: [Reliability]) = align(center)[
  #reliability-curve(series, palette: palette, title: title)
  #v(0.4em)
  #calib-legend(series, palette: palette)
]
