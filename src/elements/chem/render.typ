// Coordinate renderer: LayoutOut (from the WASM core) -> CeTZ drawing.
//
// Reuses alchemist's drawing idiom (bonds + labels) but driven by absolute
// coordinates from the layout engine, not turtle angles. Each atom is exposed
// as a named CeTZ anchor so hand-drawn annotations / cross-links can attach.

#import "@preview/cetz:0.5.2"
#import cetz.draw: *

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

// ── one bond between two points ─────────────────────────────────────────────
#let draw-bond(a, b, bond, config) = {
  // GR-6: in circle mode, aromatic ring bonds are plain single lines.
  if config.at("aromatic", default: none) == "circle" and bond.at("aromatic", default: false) {
    line(a, b, stroke: config.stroke)
    return
  }
  let gap = config.double-gap
  let sh = config.double-shorten
  let dx = b.at(0) - a.at(0)
  let dy = b.at(1) - a.at(1)
  let len = calc.sqrt(dx * dx + dy * dy)
  if len == 0 { return }
  let (ux, uy) = (dx / len, dy / len)
  let (nx, ny) = (-uy, ux) // perpendicular
  let stroke = config.stroke

  let kind = bond.kind
  if kind == "single" {
    line(a, b, stroke: stroke)
  } else if kind == "double" {
    if bond.inner.x == 0 and bond.inner.y == 0 {
      // symmetric
      line(
        (a.at(0) + nx * gap / 2, a.at(1) + ny * gap / 2),
        (b.at(0) + nx * gap / 2, b.at(1) + ny * gap / 2),
        stroke: stroke,
      )
      line(
        (a.at(0) - nx * gap / 2, a.at(1) - ny * gap / 2),
        (b.at(0) - nx * gap / 2, b.at(1) - ny * gap / 2),
        stroke: stroke,
      )
    } else {
      // main line + inner offset line toward ring centroid
      line(a, b, stroke: stroke)
      let ox = bond.inner.x * gap
      let oy = bond.inner.y * gap
      line(
        (a.at(0) + ox + ux * sh, a.at(1) + oy + uy * sh),
        (b.at(0) + ox - ux * sh, b.at(1) + oy - uy * sh),
        stroke: stroke,
      )
    }
  } else if kind == "triple" {
    line(a, b, stroke: stroke)
    line(
      (a.at(0) + nx * gap, a.at(1) + ny * gap),
      (b.at(0) + nx * gap, b.at(1) + ny * gap),
      stroke: stroke,
    )
    line(
      (a.at(0) - nx * gap, a.at(1) - ny * gap),
      (b.at(0) - nx * gap, b.at(1) - ny * gap),
      stroke: stroke,
    )
  } else if kind == "dative" {
    // GR-1.7 coordination/dative bond: an arrow from donor to acceptor
    line(a, b, stroke: stroke, mark: (end: ">"))
  } else if kind.starts-with("cram-filled") {
    // solid wedge: narrow at a, wide at b
    let w = config.wedge-width
    line(
      a,
      (b.at(0) + nx * w / 2, b.at(1) + ny * w / 2),
      (b.at(0) - nx * w / 2, b.at(1) - ny * w / 2),
      close: true,
      fill: black,
      stroke: stroke,
    )
  } else if kind.starts-with("cram-dashed") {
    // hashed wedge
    let n = config.hash-count
    for i in range(1, n + 1) {
      let t = i / (n + 1)
      let cx = a.at(0) + dx * t
      let cy = a.at(1) + dy * t
      let hw = config.wedge-width * t / 2
      line(
        (cx + nx * hw, cy + ny * hw),
        (cx - nx * hw, cy - ny * hw),
        stroke: stroke,
      )
    }
  } else {
    line(a, b, stroke: stroke)
  }
}

// ── the full molecule from a LayoutOut dict ─────────────────────────────────
#let draw-chem(layout, name: "mol", config: (:)) = {
  let cfg = (
    // GR (Brecher 2008, p.299): bond thickness ≈ the atom-label font stroke
    // (the leg of a capital "H"). Using an em-relative width keeps that match at
    // any text size; avoid < ¼× or > 4× the label stroke.
    stroke: 0.06em + black,
    double-gap: 0.13,
    double-shorten: 0.16,
    wedge-width: 0.28,
    hash-count: 6,
    label-pad: 1pt,
    label-gap: 0.16, // bond clearance (coord units) at a labeled atom
    color: false,
    atom-colors: (:),
    lone-pairs: none, // none | "dots" | "lines"
    lp-offset: 0.34,
    lp-dot-r: 0.028,
    lp-gap: 0.07,
    lp-line: 0.12,
    aromatic: none, // none (Kekulé) | "circle"
    aromatic-radius: 0.62,
    show-all-h: false,
    oxidation: (:), // GR-2.1.4: atom-id (int or str) -> oxidation number
    ionic-bonds: false, // GR-7.1: dotted lines between disconnected components
    script-size: 0.85em, // size of sub/superscripts (charge, count, oxidation)
    brackets: none, // GR-5.7: enclose in [ ]; value = charge shown outside (e.g. -2)
    delocalize: (), // GR-5.4: list of (from, to) atom-id pairs drawn delocalized
    variable-attach: (), // GR-9.4: list of (from-atom-id, (ring-atom-ids...))
    vertical: (), // GR-2.1.7: atom-ids whose label stacks vertically (multi-line)
    bent: (), // GR-1.5: list of (from, to) atom-id pairs drawn with a kink
    bent-offset: 0.3, // perpendicular kink size for bent bonds
    multi-centre: (), // GR-1.9: list of (from-atom-id, (atom-ids...)) hapto bonds
    pseudo: (), // GR-12: list of (from, to) atom-id pairs drawn as a wavy connector
    partial-charge: (:), // GR-5.6: atom-id (str) -> "+" or "-" drawn as δ+/δ-
    arrows: (), // electron-pushing curly arrows: list of (from-id, to-id) or (from, to, side)
    arrow-paint: rgb("#c00"),
  ) + config

  let pos(i) = {
    let a = layout.atoms.at(i)
    (a.pos.x, a.pos.y)
  }

  // an atom is drawn as a text label (not a bare skeletal vertex)
  let is-labeled(i) = {
    let a = layout.atoms.at(i)
    not a.skeletal or (cfg.show-all-h and a.at("implicit_h", default: 0) > 0)
  }

  // pull a bond's endpoints back to leave a gap at labeled atoms
  let trim-ends(fi, ti) = {
    let a = pos(fi)
    let b = pos(ti)
    let dx = b.at(0) - a.at(0)
    let dy = b.at(1) - a.at(1)
    let len = calc.sqrt(dx * dx + dy * dy)
    if len == 0 { return (a, b) }
    let g = cfg.label-gap
    let (ux, uy) = (dx / len, dy / len)
    let na = if is-labeled(fi) { (a.at(0) + ux * g, a.at(1) + uy * g) } else { a }
    let nb = if is-labeled(ti) { (b.at(0) - ux * g, b.at(1) - uy * g) } else { b }
    (na, nb)
  }

  // bonds (under labels)
  let is-deloc(f, t) = cfg.delocalize.any(p => (p.at(0) == f and p.at(1) == t) or (p.at(0) == t and p.at(1) == f))
  let is-bent(f, t) = cfg.bent.any(p => (p.at(0) == f and p.at(1) == t) or (p.at(0) == t and p.at(1) == f))
  let is-pseudo(f, t) = cfg.pseudo.any(p => (p.at(0) == f and p.at(1) == t) or (p.at(0) == t and p.at(1) == f))
  for bond in layout.bonds {
    let (a, b) = trim-ends(bond.from, bond.to)
    if is-pseudo(bond.from, bond.to) {
      // GR-12 pseudobond: a wavy connector (squiggle) between two atoms
      let dx = b.at(0) - a.at(0)
      let dy = b.at(1) - a.at(1)
      let len = calc.sqrt(dx * dx + dy * dy)
      let (ux, uy) = if len == 0 { (0, 0) } else { (dx / len, dy / len) }
      let (nx, ny) = (-uy, ux)
      let waves = int(calc.max(3, calc.round(len / 0.18)))
      let amp = 0.09
      let steps = waves * 6
      let pts = range(steps + 1).map(k => {
        let t = k / steps
        let s = calc.sin(t * waves * calc.pi) * amp
        (a.at(0) + ux * len * t + nx * s, a.at(1) + uy * len * t + ny * s)
      })
      line(..pts, stroke: cfg.stroke)
    } else if is-bent(bond.from, bond.to) {
      // GR-1.5 bent bond: a polyline with a perpendicular kink at the midpoint
      let dx = b.at(0) - a.at(0)
      let dy = b.at(1) - a.at(1)
      let len = calc.sqrt(dx * dx + dy * dy)
      let (nx, ny) = if len == 0 { (0, 0) } else { (-dy / len, dx / len) }
      let mid = ((a.at(0) + b.at(0)) / 2 + nx * cfg.bent-offset, (a.at(1) + b.at(1)) / 2 + ny * cfg.bent-offset)
      line(a, mid, b, stroke: cfg.stroke)
    } else if is-deloc(bond.from, bond.to) {
      // GR-5.4 delocalized bond: solid line + parallel dashed line (½ order)
      let dx = b.at(0) - a.at(0)
      let dy = b.at(1) - a.at(1)
      let len = calc.sqrt(dx * dx + dy * dy)
      let (nx, ny) = if len == 0 { (0, 0) } else { (-dy / len, dx / len) }
      let g = cfg.double-gap
      // offset the dashed line toward the ring centroid when known, else +normal
      let s = if bond.inner.x != 0 or bond.inner.y != 0 {
        if bond.inner.x * nx + bond.inner.y * ny < 0 { -1 } else { 1 }
      } else { 1 }
      line(a, b, stroke: cfg.stroke)
      let sh = cfg.double-shorten
      line(
        (a.at(0) + nx * g * s + dx / len * sh, a.at(1) + ny * g * s + dy / len * sh),
        (b.at(0) + nx * g * s - dx / len * sh, b.at(1) + ny * g * s - dy / len * sh),
        stroke: (paint: black, thickness: cfg.stroke.thickness, dash: "dashed"),
      )
    } else {
      draw-bond(a, b, bond, cfg)
    }
  }

  // GR-7.1: dotted ionic bonds between disconnected components (salts)
  if cfg.ionic-bonds {
    let n = layout.atoms.len()
    // union-find over covalent bonds to label connected components
    let parent = range(n)
    let find(x) = {
      let r = x
      while parent.at(r) != r { r = parent.at(r) }
      r
    }
    for bond in layout.bonds {
      let ra = find(bond.from)
      let rb = find(bond.to)
      if ra != rb { parent.at(ra) = rb }
    }
    let comp = range(n).map(find)
    // for each unordered pair of distinct components, connect their closest atoms
    let roots = comp.dedup()
    for i in range(roots.len()) {
      for j in range(i + 1, roots.len()) {
        let (ra, rb) = (roots.at(i), roots.at(j))
        let best = none
        for a in range(n) {
          if comp.at(a) != ra { continue }
          for b in range(n) {
            if comp.at(b) != rb { continue }
            let pa = pos(a)
            let pb = pos(b)
            let d = (pa.at(0) - pb.at(0)) * (pa.at(0) - pb.at(0)) + (pa.at(1) - pb.at(1)) * (pa.at(1) - pb.at(1))
            if best == none or d < best.at(0) { best = (d, a, b) }
          }
        }
        if best != none {
          line(pos(best.at(1)), pos(best.at(2)), stroke: (paint: gray, thickness: cfg.stroke.thickness, dash: "dotted"))
        }
      }
    }
  }

  // GR-9.4: variable attachment — a bond from a substituent into a ring centre,
  // indicating it attaches at an unspecified ring position.
  for (from, ring-ids) in cfg.variable-attach {
    let n = ring-ids.len()
    let cx = ring-ids.map(i => pos(i).at(0)).sum() / n
    let cy = ring-ids.map(i => pos(i).at(1)).sum() / n
    let a = pos(from)
    // stop a bit short of the centre so the line ends inside the ring
    let dx = cx - a.at(0)
    let dy = cy - a.at(1)
    line(a, (a.at(0) + dx * 0.7, a.at(1) + dy * 0.7), stroke: cfg.stroke)
  }

  // GR-1.9: multi-centre (η/hapto) bond — a full bond from an atom to the
  // centroid of a set of atoms (e.g. metallocene ring → metal).
  for (from, ids) in cfg.multi-centre {
    let n = ids.len()
    let cx = ids.map(i => pos(i).at(0)).sum() / n
    let cy = ids.map(i => pos(i).at(1)).sum() / n
    line(pos(from), (cx, cy), stroke: cfg.stroke)
  }

  // GR-6: aromatic delocalization circles
  if cfg.aromatic == "circle" {
    for ring in layout.at("aromatic_rings", default: ()) {
      circle(
        (ring.center.x, ring.center.y),
        radius: ring.radius * cfg.aromatic-radius,
        stroke: cfg.stroke,
      )
    }
  }

  // lone pairs (under labels, drawn from atom center outward)
  if cfg.lone-pairs != none {
    for atom in layout.atoms {
      let clr = element-color(atom.element, cfg.color, cfg.atom-colors)
      let p = (atom.pos.x, atom.pos.y)
      for dir in atom.lone_pair_dirs {
        let (dx, dy) = (dir.x, dir.y)
        let (px, py) = (p.at(0) + dx * cfg.lp-offset, p.at(1) + dy * cfg.lp-offset)
        let (nx, ny) = (-dy, dx) // perpendicular
        if cfg.lone-pairs == "dots" {
          circle((px + nx * cfg.lp-gap, py + ny * cfg.lp-gap), radius: cfg.lp-dot-r, fill: clr, stroke: none)
          circle((px - nx * cfg.lp-gap, py - ny * cfg.lp-gap), radius: cfg.lp-dot-r, fill: clr, stroke: none)
        } else {
          line(
            (px + nx * cfg.lp-line, py + ny * cfg.lp-line),
            (px - nx * cfg.lp-line, py - ny * cfg.lp-line),
            stroke: cfg.stroke,
          )
        }
      }
    }
  }

  // radicals (GR-5.3): single dots, always shown
  for atom in layout.atoms {
    let clr = element-color(atom.element, cfg.color, cfg.atom-colors)
    let p = (atom.pos.x, atom.pos.y)
    for dir in atom.at("radical_dirs", default: ()) {
      circle(
        (p.at(0) + dir.x * cfg.lp-offset, p.at(1) + dir.y * cfg.lp-offset),
        radius: cfg.lp-dot-r,
        fill: clr,
        stroke: none,
      )
    }
  }

  // GR-5.6: partial charges, drawn as a δ+/δ- glyph above the atom
  for (id, sign) in cfg.partial-charge {
    let a = layout.atoms.at(int(id))
    let s = if sign == "+" { "+" } else { "−" }
    content(
      (a.pos.x, a.pos.y + 0.42),
      text(size: cfg.script-size * 1.1, [δ] + super(size: cfg.script-size, s)),
    )
  }

  // atom labels + named anchors
  for atom in layout.atoms {
    let p = (atom.pos.x, atom.pos.y)
    let ih = atom.at("implicit_h", default: 0)
    let as-label = not atom.skeletal or (cfg.show-all-h and ih > 0)
    let anchorname = name + "-a" + str(atom.id)
    if not as-label {
      hide(circle(p, radius: 0.001, name: anchorname))
    } else {
      let clr = element-color(atom.element, cfg.color, cfg.atom-colors)
      let dir = atom.at("label_dir", default: "center")
      let ox = cfg.oxidation.at(str(atom.id), default: none)
      let vert = cfg.vertical.contains(atom.id) and not atom.skeletal
      let ss = cfg.script-size
      // Element groups: skeletal show-all-h carbons synthesize C(+Hn); others
      // parse the condensed label (CH3 -> [(C,""),(H,3)]).
      let groups = if atom.skeletal {
        ((sym: atom.element, sub: ""),) + (
          if ih == 1 { ((sym: "H", sub: ""),) } else if ih > 1 { ((sym: "H", sub: str(ih)),) } else { () }
        )
      } else { element-groups(atom.text) }
      // Plain element-only labels: anchor the *bonded* element symbol on the
      // atom point and let the H-count / further groups trail off to the side,
      // so the bond vertex sits on the symbol (not the whole-string box centre)
      // and subscripts don't drag the text up or sideways.
      let plain = groups != none and atom.charge == 0 and atom.isotope == none and ox == none and not vert
      if plain {
        let reversed = dir == "left"
        let ordered = if reversed { groups.rev() } else { groups }
        let ci = if reversed { ordered.len() - 1 } else { 0 }
        let conn = ordered.at(ci)
        let lefts = ordered.slice(0, ci)
        let rights = ordered.slice(ci + 1)
        let gc(g) = g.sym + if g.sub != "" { sub(size: ss, g.sub) }
        let pad = cfg.label-pad
        // white mask; pad only the outer edges so adjacent pieces stay tight.
        let wht(b, lp, rp) = box(
          fill: white,
          inset: (left: if lp { pad } else { 0pt }, right: if rp { pad } else { 0pt }, top: pad, bottom: pad),
          text(fill: clr, b),
        )
        let has-right = conn.sub != "" or rights.len() > 0
        let has-left = lefts.len() > 0
        // GR-2.1.5: the bond points to the CENTER of the connecting element
        // symbol (positioned where the bond would have ended if unlabeled).
        content(p, wht(conn.sym, not has-left, not has-right), anchor: "center", name: anchorname)
        if has-right {
          let rb = {
            if conn.sub != "" { sub(size: ss, conn.sub) }
            rights.map(gc).join()
          }
          content(anchorname + ".base-east", wht(rb, false, true), anchor: "base-west")
        }
        if has-left {
          content(anchorname + ".base-west", wht(lefts.map(gc).join(), true, false), anchor: "base-east")
        }
      } else {
        // Fallback for charged / isotopic / oxidation / vertically-stacked labels.
        let reversed = dir == "left"
        let anchor = if dir == "left" { "east" } else if dir == "right" { "west" } else { "center" }
        let vgroups = if vert { element-groups(atom.text) } else { none }
        let body = if vert and vgroups != none {
          text(fill: clr, stack(
            dir: ttb,
            spacing: 0.15em,
            ..vgroups.map(g => g.sym + if g.sub != "" { sub(size: ss, g.sub) }),
          ))
        } else {
          text(fill: clr, format-label(atom, reversed: reversed, oxidation: ox, script-size: ss))
        }
        content(p, box(fill: white, inset: cfg.label-pad, body), anchor: if vert { "center" } else { anchor }, name: anchorname)
      }
    }
  }

  // Electron-pushing curly arrows (drawn on top). The arc is offset perpendicular
  // to the bond so it stays clear of the bond lines, and curves in to point at the
  // target atom — the same geometry as `curly-arrow`, but addressed by atom id so
  // the user never touches Cetz anchors.
  for arr in cfg.arrows {
    let a = pos(arr.at(0))
    let b = pos(arr.at(1))
    let side = if arr.len() > 2 { arr.at(2) } else { 1 }
    let dx = b.at(0) - a.at(0)
    let dy = b.at(1) - a.at(1)
    let len = calc.max(calc.sqrt(dx * dx + dy * dy), 0.0001)
    let (ux, uy) = (dx / len, dy / len)
    let (px, py) = (-uy * side, ux * side)
    let (off, bend, pad) = (0.3, 0.55, 0.16)
    let p0 = (a.at(0) + ux * pad + px * off, a.at(1) + uy * pad + py * off)
    let p3 = (b.at(0) - ux * pad * 0.2 + px * off * 0.55, b.at(1) - uy * pad * 0.2 + py * off * 0.55)
    let c1 = (p0.at(0) + ux * len * 0.2 + px * bend, p0.at(1) + uy * len * 0.2 + py * bend)
    let c2 = (p3.at(0) + px * bend, p3.at(1) + py * bend)
    bezier(p0, p3, c1, c2, stroke: 0.7pt + cfg.arrow-paint, mark: (end: ">", scale: 0.85))
  }

  // GR-5.7 / GR-5.4: enclose a polyatomic ion in square brackets with the
  // (delocalized) charge placed outside, to the upper right.
  if cfg.brackets != none {
    let xs = layout.atoms.map(a => a.pos.x)
    let ys = layout.atoms.map(a => a.pos.y)
    let pad = 0.45
    let (x0, x1) = (calc.min(..xs) - pad, calc.max(..xs) + pad)
    let (y0, y1) = (calc.min(..ys) - pad, calc.max(..ys) + pad)
    let tick = 0.22
    // left bracket
    line((x0 + tick, y1), (x0, y1), (x0, y0), (x0 + tick, y0), stroke: cfg.stroke)
    // right bracket
    line((x1 - tick, y1), (x1, y1), (x1, y0), (x1 - tick, y0), stroke: cfg.stroke)
    // charge label outside the upper-right corner
    let q = cfg.brackets
    if q != 0 {
      let m = calc.abs(q)
      let sign = if q > 0 { "+" } else { "−" }
      content(
        (x1 + 0.12, y1),
        text(super(size: cfg.script-size, (if m > 1 { str(m) } else { "" }) + sign)),
        anchor: "south-west",
      )
    }
  }
}
