// ============================================================
// figstyle.typ — the per-figure YAML convention loader (D-14 / D-15).
//
// GOAL (D-14): every figure gets a sibling `manuscript/figures/<name>.yaml` declaring ALL
// graphical parameters, so the user changes any visual aspect of a figure by editing YAML —
// nothing graphical is hard-coded in the figure's `.typ`. Colours are named by palette ROLE
// (resolved against `/colours.yaml`, the single source of hex) with an optional per-element
// "#rrggbb" override.
//
// USAGE (in figures/<name>.typ, compiled with --root manuscript):
//   #import "../lib/figstyle.typ": figstyle
//   #let fs = figstyle("f1_speedup")          // loads figures/f1_speedup.yaml
//   #let w  = fs.panel.width * 1cm            // geometry values are plain numbers (units below)
//   #let c  = (fs.color)("tool_series.bayesinteractomics")  // role -> rgb; or (fs.color)("#5b21b6")
//   ... fs.axes.xlim, fs.axes.xlabel, fs.series, fs.legend.position ...
//
// NOTE (Typst idiom): `color` is a FUNCTION stored in the returned dict. Typst forbids calling a
// dict-stored function with method syntax (`fs.color(x)` errors), so wrap the field access in
// parentheses: `(fs.color)(role_or_hex)`. All other keys are plain values: `fs.panel.width`, etc.
//
// `#let`-only — `#import` emits NOTHING (Pitfall 4).
//
// ------------------------------------------------------------
// YAML SCHEMA (every key + type + default). Any key you omit falls back to the default here;
// later phases override only what they need. Units are documented per key (YAML holds plain
// numbers/strings — the figure `.typ` applies the unit, e.g. `fs.panel.width * 1cm`).
//
//   panel:                      # figure geometry
//     width:   <number>  cm     (default 8)    overall panel width  in cm
//     height:  <number>  cm     (default 6)    overall panel height in cm
//     gutter:  <number>  pt      (default 10)   inter-panel grid gutter in pt
//     scale:   <number>  unitless(default 1.0) global scale factor
//
//   axes:
//     xlim:   <array[2] number | null>  (default null)  x-axis limits (lo, hi)
//     ylim:   <array[2] number | null>  (default null)  y-axis limits (lo, hi)
//     xticks: <array number | null>     (default null)  explicit x tick positions
//     yticks: <array number | null>     (default null)  explicit y tick positions
//     xlabel: <string | null>           (default null)  x-axis label
//     ylabel: <string | null>           (default null)  y-axis label
//     title:  <string | null>           (default null)  panel/axes title
//
//   series:  <array of dicts>  (default ())  per-series styling; each entry may set:
//     color:     <string>  palette role ("tool_series.hero") OR "#rrggbb" override
//     thickness: <number>  stroke thickness in pt        (suggested default 1.5)
//     dash:      <string>  "solid" | "dashed" | "dotted" (suggested default "solid")
//     mark:      <string>  marker glyph / name           (suggested default "none")
//     bar_width: <number>  bar width (data units or fraction)
//     offset:    <number>  bar/group offset
//     label:     <string>  legend label for the series
//   (series entries are merged wholesale — the YAML list REPLACES the default empty list.)
//
//   legend:
//     position: <string>  "top" | "bottom" | "left" | "right" | "none"  (default "top")
//     cols:     <number>  legend columns                                 (default 1)
//
//   text:
//     size: <number>  body text size in pt  (default 9)
//
// (fs.color)(role_or_hex) -> rgb:   (call with parentheses — see NOTE above)
//   * "#rrggbb"            -> rgb("#rrggbb")            (explicit per-figure override, D-14)
//   * "group.role"         -> rgb(palette.group.role)  (dotted path into /colours.yaml)
//   * an already-resolved colour value -> returned unchanged
// ============================================================

// Documented defaults covering ALL graphical groups (D-14). Geometry as plain numbers so YAML
// round-trips; the figure `.typ` applies units.
#let _DEFAULTS = (
  panel: (width: 8, height: 6, gutter: 10, scale: 1.0),
  axes: (
    xlim: none, ylim: none, xticks: none, yticks: none,
    xlabel: none, ylabel: none, title: none,
  ),
  series: (),
  legend: (position: "top", cols: 1),
  text: (size: 9),
)

// Deep-merge `over` onto `base`: nested dicts merge recursively; scalars/arrays in `over`
// replace those in `base` wholesale (so a YAML `series:` list REPLACES the default list).
#let _merge(base, over) = {
  if type(base) == dictionary and type(over) == dictionary {
    let out = base
    for (k, v) in over {
      if k in out {
        out.insert(k, _merge(out.at(k), v))
      } else {
        out.insert(k, v)
      }
    }
    out
  } else {
    over
  }
}

// Traverse a dotted path ("tool_series.bayesinteractomics") into the palette dict.
#let _resolve(palette, path) = {
  let node = palette
  for part in path.split(".") {
    node = node.at(part)
  }
  node
}

// figstyle(name): load /colours.yaml + /figures/<name>.yaml, merge over defaults, and return a
// dict exposing the merged panel/axes/series/legend/text, a `color(role_or_hex)` helper, the raw
// `palette`, and the raw `spec`.
#let figstyle(name) = {
  let palette = yaml("/colours.yaml")
  let spec = yaml("/figures/" + name + ".yaml")
  let merged = _merge(_DEFAULTS, spec)
  // Named `_color` (not `color`) so the `type(role) == color` guard still sees Typst's
  // built-in `color` type rather than shadowing it with this closure.
  let _color = role => {
    if type(role) == color { role }
    else if type(role) == str and role.starts-with("#") { rgb(role) }
    else if type(role) == str { rgb(_resolve(palette, role)) }
    else { role }
  }
  (
    panel: merged.panel,
    axes: merged.axes,
    series: merged.series,
    legend: merged.legend,
    text: merged.text,
    color: _color,
    palette: palette,
    spec: spec,
  )
}
