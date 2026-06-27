// Label + colour helpers (GR-2.1) shared by skeleton.typ / decorations.typ.
// Pure formatting: condensed-formula subscripts, isotope/charge/oxidation
// superscripts, Jmol CPK element colours.

#import "@preview/cetz:0.5.2"

// ── Jmol CPK element colors ─────────────────────────────────────────────────
#let cpk-colors = (
  H: rgb("#000000"), C: rgb("#000000"), N: rgb("#3050F8"), O: rgb("#FF0D0D"),
  F: rgb("#90E050"), Cl: rgb("#1FF01F"), Br: rgb("#A62929"), I: rgb("#940094"),
  S: rgb("#FFD123"), P: rgb("#FF8000"), B: rgb("#FFB5B5"),
)
#let element-color(elem, color, overrides) = {
  if not color { return black }
  if elem in overrides { return overrides.at(elem) }
  cpk-colors.at(elem, default: black)
}

// ── label text: subscript digit runs, charge superscript ────────────────────
#let render-formula(text, script-size: 0.85em) = {
  let parts = ()
  let buf = ""
  for ch in text.clusters() {
    if ch >= "0" and ch <= "9" {
      if buf != "" { parts.push(buf); buf = "" }
      parts.push(sub(size: script-size, ch))
    } else {
      buf += ch
    }
  }
  if buf != "" { parts.push(buf) }
  parts.join()
}

// Split a simple condensed formula into element groups (sym + subscript count).
// Returns none for labels with parens/brackets (not simply reversible).
#let element-groups(text) = {
  if "(" in text or "[" in text or "]" in text { return none }
  let groups = ()
  let chars = text.clusters()
  let i = 0
  while i < chars.len() {
    let c = chars.at(i)
    if upper(c) == c and lower(c) != c {
      // uppercase letter starts a group
      let sym = c
      i += 1
      while i < chars.len() and lower(chars.at(i)) == chars.at(i) and upper(chars.at(i)) != chars.at(i) {
        sym += chars.at(i)
        i += 1
      }
      let d = ""
      while i < chars.len() and chars.at(i) >= "0" and chars.at(i) <= "9" {
        d += chars.at(i)
        i += 1
      }
      groups.push((sym: sym, sub: d))
    } else {
      i += 1
    }
  }
  groups
}

// GR-2.1.4: oxidation number as a superscript Roman numeral.
#let roman(n) = {
  let neg = n < 0
  let n = calc.abs(n)
  if n == 0 { return "0" }
  let table = ((1000, "M"), (900, "CM"), (500, "D"), (400, "CD"), (100, "C"),
    (90, "XC"), (50, "L"), (40, "XL"), (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I"))
  let s = ""
  for (v, sym) in table {
    while n >= v { s += sym; n -= v }
  }
  if neg { "−" + s } else { s }
}

// GR-2.1.6: reversed labels read the element groups in reverse (CH3 -> H3C),
// keeping each subscript with its element. `reversed` only applies to simple
// element-only labels.
#let format-label(atom, reversed: false, oxidation: none, script-size: 0.85em) = {
  let groups = element-groups(atom.text)
  let body = if reversed and groups != none {
    groups.rev().map(g => g.sym + if g.sub != "" { sub(size: script-size, g.sub) }).join()
  } else {
    render-formula(atom.text, script-size: script-size)
  }
  // GR-2.1.3: isotope as a left superscript
  if atom.isotope != none {
    body = super(size: script-size, str(atom.isotope)) + body
  }
  // GR-2.1.4: oxidation number as a superscript Roman numeral
  if oxidation != none {
    body = body + super(size: script-size, roman(oxidation))
  }
  if atom.charge != 0 {
    let n = calc.abs(atom.charge)
    let sign = if atom.charge > 0 { "+" } else { "−" }
    body = body + super(size: script-size, (if n > 1 { str(n) } else { "" }) + sign)
  }
  body
}
