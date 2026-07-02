// f1_speedup.typ — Figure F1 (speedup-vs-RMSE), the ONE worked reference figure (D-16).
//
// FULLY YAML-DRIVEN (D-14): every graphical value — panel geometry, axis limits/ticks/labels,
// series colours (as palette ROLES), stroke thickness, marks, legend layout — is read from
// figures/f1_speedup.yaml through lib/figstyle.typ. Nothing graphical is hard-coded here; numeric
// data lives in f1_speedup_data/speedup.csv. Restyle F1 by editing the YAML, never this file.
//
// Standalone preview:  typst compile --root manuscript manuscript/figures/f1_speedup.typ out.pdf
// Integrated:          #import "f1_speedup.typ": f1-speedup-figure  (in figures/specs.typ)
//
// NOTE (units): the YAML holds panel geometry as plain cm numbers; we apply the unit here as
// `n * 10mm` (= n cm). The DIMENSION VALUE still comes entirely from the YAML — the `10mm` is a
// unit-conversion constant, not a hard-coded panel size (D-14).

#import "@preview/lilaq:0.6.0" as lq
#import "../lib/figstyle.typ": figstyle

#let fs   = figstyle("f1_speedup")
#let data = csv("f1_speedup_data/speedup.csv")

// --- helpers ---------------------------------------------------------------
#let _cm = n => n * 10mm                                   // YAML cm number -> length
#let _dash = d => if d == "solid" { none } else { d }      // YAML dash -> lilaq dash
#let _mark = m => if m == "none" { none } else { m }       // YAML mark -> lilaq mark

// series lookup by `key` — all styling flows from the YAML
#let _ser = (:)
#for s in fs.series { _ser.insert(s.key, s) }

// (rmse, speedup) for a CSV regime (cols: method,regime,rmse,interval_width,wall_clock_s,speedup)
#let _row = regime => {
  let r = data.slice(1).filter(row => row.at(1) == regime)
  (float(r.at(0).at(2)), float(r.at(0).at(5)))
}
#let (advi-x, advi-y) = _row("per-dataset")
#let (full-x, full-y) = _row("full-workload")
#let (fwd-x,  fwd-y)  = _row("forward-pass")

// a single styled point pulled from the named YAML series
#let _point(key, x, y) = {
  let s = _ser.at(key)
  lq.plot((x,), (y,), stroke: none, mark: _mark(s.mark), color: (fs.color)(s.color), label: s.label)
}

// --- panel (every parameter from fs) ---------------------------------------
#let _base = _ser.at("advi_baseline")

#let _panel = lq.diagram(
  width:  _cm(fs.panel.width),
  height: _cm(fs.panel.height),
  xlim: fs.axes.xlim,
  ylim: fs.axes.ylim,
  yscale: fs.axes.yscale,
  xlabel: fs.axes.xlabel,
  ylabel: fs.axes.ylabel,
  title: fs.axes.title,
  xaxis: (ticks: fs.axes.xticks),
  yaxis: (ticks: fs.axes.yticks),
  legend: none,
  // ADVI 1x reference line across the x-range (styling from the YAML series)
  lq.plot(fs.axes.xlim, (advi-y, advi-y), mark: none,
    stroke: (paint: (fs.color)(_base.color), thickness: _base.thickness * 1pt, dash: _dash(_base.dash))),
  // data points (RMSE parity, three speedup regimes)
  _point("advi", advi-x, advi-y),
  _point("npe_full", full-x, full-y),
  _point("npe_forward", fwd-x, fwd-y),
)

// --- manual legend (labels/colours/cols all from the YAML) ------------------
#let _swatch(s) = box(width: 0.9em, height: 0.9em, radius: 1pt, fill: (fs.color)(s.color), baseline: 0.12em)
#let _leg-entry(s) = box[#_swatch(s) #h(3pt)#s.label]
#let _legend = {
  set text(size: fs.text.size * 1pt)
  align(center, grid(
    columns: (auto,) * fs.legend.cols,
    column-gutter: 1.6em, row-gutter: 0.5em,
    align: left + horizon,
    ..fs.series.map(_leg-entry),
  ))
}

#let _body = {
  if fs.legend.position == "top" { _legend; v(2mm); _panel }
  else if fs.legend.position == "bottom" { _panel; v(2mm); _legend }
  else { _panel }
}

// Exported figure body — global scale from the YAML.
#let f1-speedup-figure = scale(
  x: fs.panel.scale * 100%, y: fs.panel.scale * 100%, reflow: true, _body,
)

// --- standalone preview (discarded when this file is #imported as a module) --
#set page(width: auto, height: auto, margin: 4pt)
#set text(font: ("Libertinus Serif", "Linux Libertine", "DejaVu Serif", "Times New Roman"),
  size: fs.text.size * 1pt, lang: "en", region: "GB")
#f1-speedup-figure
