// ProteinCoLoc v2.0 — AmortizedColoc manuscript entry.
//
// Founded on the BayesInteractomics_Method_paper Typst template (D-13): lib/, styles/,
// fonts/, colours.yaml and build.ps1 are COPIED in from that read-only source; main.typ is
// re-authored here (not copied) using the template's wiring.
//
// Venue-neutral (D-02) — candidate venues, NOT locked: Bioinformatics | Nature Methods |
// eLife-style methods paper. Journal template lock-in is deferred to Phase 16.
//
// Typst: the source template pinned 0.14; this machine runs 0.15.0 (build.ps1 defaults to the
// installed `typst`, keeping the $env:TYPST override — D-17). Template packages (abbr, zero,
// smartaref, wordometer, lilaq, alexandria) auto-download from the Typst registry.
//
// Build (exit-code-correct, D-17): `pwsh -File manuscript/build.ps1` (primary) or
// `bash manuscript/build.sh` (fallback), each -> manuscript/main.pdf.

#import "lib/template.typ": *
#import "lib/tiered-bib.typ": *
#import "lib/helpers.typ": landscapetable
#import "claims.typ": claims, claim_table

#show: project.with(
  title: [ProteinCoLoc v2.0 — AmortizedColoc: amortized, calibrated colocalization inference],
  authors: (
    (name: "Manuel Seefelder", email: "seefelder.manuel@gmail.com", affiliation: "1"),
  ),
  aff: (
    (id: "1", dept: "Department of Gene Therapy, Ulm University"),
  ),
  abstract: [
    // TODO(phase 16): author the final abstract. Compile-safe placeholder below.
    _Draft slot._ ProteinCoLoc v2.0 trains a neural posterior estimator and a neural ratio
    estimator on the package's existing patch-correlation summary statistic, delivering
    calibrated colocalization posteriors and Bayes factors in milliseconds per dataset — with
    a simulation-based-calibration coverage proof (under the simulator) and an honest
    out-of-distribution flag.
  ],
  keywords: [Colocalization, Amortized inference, Simulation-based calibration, Neural posterior estimation, Bayes factor, Fluorescence microscopy],
  corresponding: [Manuel Seefelder, #link("mailto:seefelder.manuel@gmail.com"), Department of Gene Therapy, Ulm University],
  header: "Seefelder 2026",
)

// Tiered bibliography (Nature-Methods two-tier: main "References" + "Methods references").
// Each reference appears once, in the highest tier it is cited in; numbering is global.
// Related Work already cites 4 keys (H tier), so the main "References" list renders NON-empty;
// the Methods list stays empty until Phase 5/7 wire @cites under the M tier. Each key cited must
// exist in refs.bib (a dangling @cite fails the build).
#show: tiered-bib(enabled: true, bib: "refs.bib", style: "styles/nature.csl", read: path => read(path))

// ---- Body ----
#include "sections/introduction.typ"

= Related Work

#include "sections/related_work.typ"

= Results

#include "sections/results.typ"

= Discussion

#include "sections/discussion.typ"

// Main-text references (H tier).
#[
  #set par.line(numbering: none)  // reference lists are not line-numbered
  #print-references(tiers: ("H",), read: path => read(path))
]

= Methods
#set-tier("M")  // citations below default to the Methods tier

#include "sections/methods.typ"

// Methods references (M tier).
#[
  #set par.line(numbering: none)
  #print-references(tiers: ("M",), read: path => read(path))
]

= Claims

// Wide 6-column spine — placed on its own landscape page (D-user: large tables may be landscape
// on a separate page) so every column reads without cramping.
#landscapetable(
  claim_table(claims),
  caption: [Claim spine — the contract mapping each scientific claim to its differentiator
    axis, supporting phase(s), figure/experiment ID, and status. Later plans populate the full
    claim set; a status column distinguishes _supported (measured)_ from _reserved slot_.],
  lbl: "tab:claims",
)

= Figure Specifications

#include "figures/specs.typ"

= Data and Code Availability

#include "sections/availability.typ"
