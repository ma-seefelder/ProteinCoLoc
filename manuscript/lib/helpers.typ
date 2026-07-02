// Shared layout helpers used by both main.typ and supplement.typ.

// ---------- monospace gene symbols ----------
#let g(x) = text(font: ("DejaVu Sans Mono", "Consolas"), size: 0.85em, x)

// ---------- inline figure + panel cross-reference ----------
// Renders "Figure 1A" (panel letter attached) rather than "Figure 1, panel A".
// Usage: #figp(<fig-overview>, [C])  or  #figp(<fig-hap40>, [A,B])
#let figp(lbl, p) = [#ref(lbl)#p]

// ---------- panel composition (bold A/B identifiers, no overlap) ----------
#let plab(s) = text(weight: "bold", size: 11pt, s)
#let panel2(la, a, lb, b, ratio: (1fr, 1fr)) = grid(
  columns: ratio, column-gutter: 10pt, row-gutter: 2pt,
  align: left + bottom,
  plab(la), plab(lb),
  align(center + bottom, a), align(center + bottom, b),
)
#let panel3(la, a, lb, b, lc, c) = grid(
  columns: (1fr, 1fr, 1fr), column-gutter: 10pt, row-gutter: 2pt,
  align: left + bottom,
  plab(la), plab(lb), plab(lc),
  align(center + bottom, a), align(center + bottom, b), align(center + bottom, c),
)

// ---------- full-width (page-spanning) figure helper ----------
// kind: image is forced so a float whose body contains a table (e.g. the HAP40
// engagement heatmap) is still numbered/captioned as a Figure, not auto-detected
// as a Table.
#let widefig(body, caption: none, lbl: none) = [
  #place(top, scope: "parent", float: true)[
    #figure(body, caption: caption, kind: image)
    #if lbl != none { label(lbl) }
  ]
]

// ---------- breakable supplementary table helper ----------
// kind: "supplement_table" gets a top caption + breakable block from the template;
// cross-references render as "Supplementary Table N".
#let supptable(body, caption: none, lbl: none) = [
  #figure(body, caption: caption, kind: "supplement_table", supplement: [Supplementary Table])
  #if lbl != none { label(lbl) }
]

// ---------- large full-page LANDSCAPE table helper ----------
// For main-text tables too wide even for the two-column parent scope (e.g. the claim
// spine, the related-work capability matrix): emit a DEDICATED landscape page (single
// column) so the table gets the full page height as its width. kind: table -> numbered
// "Table N" with the template's top caption. This is the standard placement for oversized
// tables — prefer it over cramming a wide table into the two-column body.
#let landscapetable(body, caption: none, lbl: none) = [
  #page(flipped: true, columns: 1)[
    #figure(body, caption: caption, kind: table)
    #if lbl != none { label(lbl) }
  ]
]

// ---------- supplementary figure helper (own counter, spans both columns) ----------
// kind: "supplement_figure" gets an independent "Supplementary Figure N" counter;
// cross-references render accordingly. Floating variant for the two-column main body.
#let suppfig(body, caption: none, lbl: none) = [
  #place(top, scope: "parent", float: true)[
    #figure(body, caption: caption, kind: "supplement_figure", supplement: [Supplementary Figure])
    #if lbl != none { label(lbl) }
  ]
]

// Non-floating variant: places the figure inline in document order. Used in the
// single-column supplement so figures appear where written (no empty heading page).
#let suppfigblock(body, caption: none, lbl: none) = [
  #figure(body, caption: caption, kind: "supplement_figure", supplement: [Supplementary Figure])
  #if lbl != none { label(lbl) }
]
