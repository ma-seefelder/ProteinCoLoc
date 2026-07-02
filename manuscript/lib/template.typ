#import "@preview/abbr:0.2.3": *
#import "@preview/zero:0.5.0": *
#import "@preview/smartaref:0.1.0": *
#import "helpers.typ": *
#import "curves.typ": *

#import "@preview/wordometer:0.1.5": total-words, word-count

// Placeholder text for drafting — blue lorem to highlight unfinished sections.
#let bluelorem(n) = text(fill: blue, lorem(n))

// Customize list of figures/tables to show short caption in outline
#let in-outline = state("in-outline", false)
#show outline: it => {
  in-outline.update(true)
  it
  in-outline.update(false)
}

#let flex-caption(long, short) = context if in-outline.get() { short } else { long }

#let refmissing(content) = text(fill: red, content)

// The project function defines how the document looks.
#let run-in-heading(title) = text(weight: "bold", title + ".") + h(0.3em)

#let project(
  title: "",
  abstract: [],
  authors: (),
  aff: (),
  keywords: [],
  corresponding: [],
  header: [],
  body,
) = {
  // Set the document's basic properties.
  set document(author: authors.map(a => a.name), title: title)

  // Page setup
  set page(
    numbering: "1",
    number-align: end,
    header: [
      #grid(
        columns: (5fr, 1fr, 2fr),
        gutter: 4pt,
        align: (left + top, center + horizon, right + top),
        "", "", text(smallcaps(header), size: 9pt),
      )
      #line(length: 100%, stroke: 0.25pt)
    ],
    header-ascent: 20%,
    paper: "a4",
    margin: (x: 2cm, y: 2.5cm), // Explicit margins are usually good
  )

  // Styling of inline fraction
  show math.equation.where(block: false): set math.frac(style: "horizontal")

  // Font definitions
  let body-font = ("Libertinus Serif", "Linux Libertine", "DejaVu Serif", "Times New Roman")
  let sans-font = ("Libertinus Sans", "Linux Biolinum", "DejaVu Sans", "Arial")

  // Set body font family.
  set text(font: body-font, lang: "en", region: "GB", size: 10pt)

  // Match math typesetting to the Libertinus family.
  show math.equation: set text(font: "Libertinus Math")

  // Style the abbreviations
  config(style: it => text(fill: black, it), space-char: " ")

  // Set paragraph spacing and heading styles
  set block(above: 1.2em, below: 1.2em)
  show heading: set text(font: sans-font)
  set par(leading: 0.75em, justify: true)

  // Title row.
  align(center)[
    #block(text(font: sans-font, weight: 900, 1.5em, title))
  ]

  // Author information.
  let author-block(authors) = {
    v(0.5em)
    align(left)[
      #text(size: 8pt)[
        #(
          authors
            .map(author => [
              #author.name
              #super(author.affiliation)
            ])
            .join(", ")
        )
      ]
    ]
    v(0.25em)
  }

  let affiliation-table(affiliations) = {
    set text(size: 8pt)
    table(
      columns: (auto, 1fr),
      align: (right, left),
      stroke: none,
      inset: (x: 0pt, y: 2pt),
      ..for affiliation in affiliations {
        ([#super(affiliation.id)], [#affiliation.dept])
      }
    )
    v(0.5em)
  }

  author-block(authors)
  affiliation-table(aff)

  // Table styling
  show table.cell: set text(size: 9pt)
  show table.cell.where(y: 0): set text(weight: "bold")
  show table.cell.where(x: 0): set text(weight: "bold")
  set table(
    inset: 5pt,
    stroke: (x, y) => if y == 0 {
      (bottom: 0.7pt + black)
    } else {
      none
    },
  )

  // Figure caption styling
  show figure.caption: set text(size: 9pt, weight: "bold")
  show figure.caption: set align(left)

  show figure.where(kind: table): set figure.caption(position: top)

  // Supplement table specific
  show figure.where(kind: "supplement_table"): set figure.caption(position: top)
  show figure.where(kind: "supplement_table"): set block(breakable: true)

  // Math equations
  set math.equation(numbering: "(1)")
  show math.equation.where(block: true): set text(size: 9pt)


  // Abstract
  if abstract != [] {
    block(
      inset: (x: 2.5%),
      [
        #set text(size: 0.9em, weight: 600)
        #line(length: 100%, stroke: 0.5pt)
        #align(center)[#text(size: 1.1em)[Abstract]]
        #abstract
        #line(length: 100%, stroke: 0.5pt)
      ],
    )
  }

  // Keywords
  if keywords != [] {
    v(0.5em)
    text(weight: 900, size: 9pt)[Keywords: ]
    text(size: 9pt)[#keywords]
    v(0.5em)
  }

  // Corresponding authors
  if corresponding != [] {
    text(weight: 900, size: 9pt)[Corresponding author: ]
    text(size: 9pt)[#corresponding]
  }

  // Divider
  v(1em)

  // Main body setup
  set par(justify: true)
  set text(weight: "regular")
  set page(columns: 2)
  set par.line(numbering: "1")

  body
}
