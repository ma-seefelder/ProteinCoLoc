// ============================================================
// tiered-bib.typ — toggleable publisher-tiered bibliography
//   H (main) -> M (methods) -> S (supplement)
// A reference appears ONCE, in the highest tier it is cited in.
// Numbering is global & continuous across the three lists.
// Graduated from .planning/spikes/009-tiered-bibliography (VALIDATED).
//
// Usage in main.typ (after the `project` show rule):
//   #import "lib/tiered-bib.typ": *
//   #show: tiered-bib(enabled: true, bib: "literature.bib",
//                     style: "styles/nature.csl", read: path => read(path))
//   ... @key citations as usual ...
//   = Methods
//   #set-tier("M")                                  // start of Methods
//   #print-references(tiers: ("H", "M"), read: path => read(path))
//   #include "supplement.typ"
// and in supplement.typ:
//   #set-tier("S")                                  // top of the SI
//   #print-references(tiers: ("S",), read: path => read(path))
//
// IMPORTANT: pass `read: path => read(path)` from the CALLING file. Typst
// resolves `read()`/bibliography paths relative to the file the call is written
// in; this module lives in lib/, but the .bib lives in the project root, so the
// caller's read closure (resolving relative to main.typ / supplement.typ) is
// threaded through for both key-parsing and Alexandria.
//
// Toggle the whole feature with `enabled:`. When false the document uses a
// single native #bibliography() (emitted by the H call) and #set-tier is an
// inert no-op — body text is byte-identical either way.
//
// NOTE: apply with `#show: tiered-bib(...)`, NOT `.with(...)` — tiered-bib
// already returns the `body => {...}` show-rule function.
// ============================================================
#import "@preview/alexandria:0.2.2" as alex

#let _RANK  = ("H": 0, "M": 1, "S": 2)
#let _NAMES = ("H", "M", "S")
#let _P = "r-"   // internal Alexandria prefix; authors never type it

#let _tier = state("tb-tier", "H")
#let _cfg  = state("tb-cfg", (enabled: false, bib: none, style: "ieee"))

// Drop at the start of a section. No-op (state only) when disabled.
#let set-tier(t) = _tier.update(t)

// Group consecutive citations into ONE collapsed citation ("4,5" / "4–6")
// instead of concatenated bare numbers ("45"). Alexandria does NOT auto-group
// adjacent @refs (native Typst does), so multi-source citations must be wrapped:
//   #cites[@Foo2020@Bar2021]            -> ^{4,5}
//   #cites[@Foo2020@Bar2021@Baz2022]    -> ^{4–6}   (CSL collapses the range)
// Single citations need no wrapper — keep writing @key as usual.
// Mode-aware: routes through Alexandria's citegroup when enabled; falls back to
// native cites (which Typst auto-groups) when disabled, so output is consistent.
#let cites(body) = context {
  let kids = if body.func() == [].func() { body.children } else { (body,) }
  let keys = kids
    .filter(x => x.func() in (ref, cite))
    .map(x => if x.func() == ref { str(x.target) } else { str(x.key) })
  if _cfg.get().enabled {
    alex.citegroup(prefix: _P, keys.map(k => ref(label(_P + k))).join())
  } else {
    keys.map(k => cite(label(k))).join()
  }
}

// Parse citation keys from one or more .bib paths (read via the caller-supplied
// `rd`) so citation routing and the tier scan ignore figure / equation /
// section cross-references.
#let _bibkeys(paths, rd) = {
  let keys = ()
  for p in paths {
    for m in rd(p).matches(regex("@\w+\s*\{\s*([^,\s]+)\s*,")) {
      keys.push(m.captures.at(0))
    }
  }
  keys
}

#let tiered-bib(
  enabled: false,
  bib: none,
  style: "ieee",
  read: path => read(path),
  titles: ([References], [Methods references], [Supplementary references]),
) = body => {
  let rd = read
  _cfg.update(c => (enabled: enabled, bib: bib, style: style, titles: titles))
  if not enabled {
    // OFF: plain Typst. Native citations; one native bibliography via print-references().
    body
  } else {
    // ON: register a NON-empty Alexandria prefix. An EMPTY prefix would make
    // Alexandria intercept EVERY @ref — including @fig:/@sec: — and crash with
    // "key fig:demo does not exist in the bibliography". So we register "r-"
    // and route only genuine bib keys to it; cross-refs are left untouched.
    let keyset = _bibkeys(if type(bib) == array { bib } else { (bib,) }, rd)
    show: alex.alexandria(prefix: _P, read: rd)              // inner rule (runs 2nd)
    show ref: it => {                                        // our rule (runs 1st)
      let k = str(it.target)
      if k.starts-with(_P) { return it }   // already routed -> Alexandria renders it
      if k not in keyset { return it }      // @fig:/@sec:/@eq: cross-ref -> leave alone
      // emit a ref (NOT cite) to avoid re-triggering this rule (show-rule depth);
      // Alexandria's inner rule then turns it into the numbered citation.
      ref(label(_P + k))
    }
    body
  }
}

// Emits one or more reference lists.
//   tiers: which tiers to render here, in order (default all three).
//   read:  caller's `path => read(path)` (see header note on path resolution).
//   size:  body text size for the entries (paper uses 8pt).
// OFF mode emits a single native bibliography, but only when "H" is among the
// requested tiers, so a separate `tiers: ("S",)` call in the SI stays empty.
#let print-references(tiers: _NAMES, read: path => read(path), size: 8pt) = context {
  let rd = read
  let cfg = _cfg.get()
  let _title(i) = heading(numbering: none, outlined: false, level: 1, cfg.titles.at(i))
  if not cfg.enabled {
    if "H" in tiers {
      _title(0)
      // Pass bytes (read via the caller's `rd`) so the native bibliography +
      // CSL resolve relative to the CALLER, not this lib/ module. A built-in
      // style name (no ".csl") is passed through unchanged.
      let bibarg = if type(cfg.bib) == array { cfg.bib.map(p => bytes(rd(p))) } else { bytes(rd(cfg.bib)) }
      let stylearg = if type(cfg.style) == str and cfg.style.ends-with(".csl") { bytes(rd(cfg.style)) } else { cfg.style }
      text(size: size, bibliography(bibarg, title: none, full: false, style: stylearg))
    }
  } else {
    let paths = if type(cfg.bib) == array { cfg.bib } else { (cfg.bib,) }
    let keyset = _bibkeys(paths, rd)
    // Highest tier (min rank) each key is cited in. We scan the ORIGINAL @key refs
    // (target in keyset); the routed `r-` refs and cross-refs are ignored.
    let rank = (:)
    for r in query(selector(ref)) {
      let k = str(r.target)
      if k in keyset {
        let rr = _RANK.at(_tier.at(r.location()))
        if k not in rank or rr < rank.at(k) { rank.insert(k, rr) }
      }
    }
    let tier-of(key) = _NAMES.at(rank.at(key, default: 0))
    alex.load-bibliography(cfg.bib, prefix: _P, style: cfg.style)
    context {
      let b = alex.get-bibliography(_P)
      // b.references are in global first-appearance (= numbering) order, each with a
      // precomputed first-field ("[1]"...). Filtering preserves the global numbers.
      let pick(t) = (..b, references: b.references.filter(r => tier-of(str(r.key)) == t))
      for (i, t) in _NAMES.enumerate() {
        if t not in tiers { continue }
        let sub = pick(t)
        if sub.references.len() > 0 {
          _title(i)
          text(size: size, alex.render-bibliography(sub, title: none))
        }
      }
    }
  }
}
