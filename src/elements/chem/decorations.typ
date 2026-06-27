// Non-bond GR decorations drawn as a thin overlay on top of the shared
// skeleton (skeleton.typ). Everything here is a pure coordinate overlay that has
// no counterpart in alchemist's link/fragment drawing: ionic dotted bonds,
// brackets, aromatic circles, electron-pushing arrows, partial/oxidation charges,
// hapto / variable-attachment, and the special bond *styles* (delocalized, bent,
// wavy) that replace a plain bond.

#import "@preview/cetz:0.5.2"
#import cetz.draw: *
#import "../../utils/utils.typ": convert-length
#import "labels.typ": element-color

// bond styles handled here instead of by the skeleton's normal link drawing.
#let is-special(cfg, f, t) = {
  let has(list) = list.any(p => (p.at(0) == f and p.at(1) == t) or (p.at(0) == t and p.at(1) == f))
  has(cfg.at("delocalize", default: ())) or has(cfg.at("bent", default: ())) or has(cfg.at("pseudo", default: ()))
}

#let draw-decorations(layout, name: "mol", config: (:), scale: 1.6) = {
  let cfg = (
    stroke: 0.06em + black,
    double-gap: 0.13,
    double-shorten: 0.16,
    color: false,
    atom-colors: (:),
    aromatic: none,
    aromatic-radius: 0.62,
    oxidation: (:),
    ionic-bonds: false,
    script-size: 0.85em,
    brackets: none,
    delocalize: (),
    variable-attach: (),
    bent: (),
    bent-offset: 0.3,
    multi-centre: (),
    pseudo: (),
    partial-charge: (:),
    arrows: (),
    arrow-paint: rgb("#c00"),
  ) + config
  get-ctx(ctx => {
  // resolve a font-relative bond length (e.g. atom-sep = 3em) to canvas-unit
  // floats so the manual vector geometry below stays float-valued.
  let s = if type(scale) == length { convert-length(ctx, scale) } else { scale }
  let stroke = cfg.stroke
  let P(i) = (layout.atoms.at(i).pos.x * s, layout.atoms.at(i).pos.y * s)
  let g = cfg.double-gap * s
  let sh = cfg.double-shorten * s

  // ── special bond styles (replace the plain bond) ──────────────────────────
  for b in layout.bonds {
    let key(list) = list.any(p => (p.at(0) == b.from and p.at(1) == b.to) or (p.at(0) == b.to and p.at(1) == b.from))
    let a = P(b.from)
    let c = P(b.to)
    let dx = c.at(0) - a.at(0)
    let dy = c.at(1) - a.at(1)
    let len = calc.max(calc.sqrt(dx * dx + dy * dy), 1e-9)
    let (ux, uy) = (dx / len, dy / len)
    let (nx, ny) = (-uy, ux)
    if key(cfg.pseudo) {
      let waves = int(calc.max(3, calc.round(len / (0.18 * s))))
      let amp = 0.09 * s
      let steps = waves * 6
      let pts = range(steps + 1).map(k => {
        let t = k / steps
        let q = calc.sin(t * waves * calc.pi) * amp
        (a.at(0) + ux * len * t + nx * q, a.at(1) + uy * len * t + ny * q)
      })
      line(..pts, stroke: stroke)
    } else if key(cfg.bent) {
      let mid = ((a.at(0) + c.at(0)) / 2 + nx * cfg.bent-offset * s, (a.at(1) + c.at(1)) / 2 + ny * cfg.bent-offset * s)
      line(a, mid, c, stroke: stroke)
    } else if key(cfg.delocalize) {
      let sgn = if b.inner.x != 0 or b.inner.y != 0 {
        if b.inner.x * nx + b.inner.y * ny < 0 { -1 } else { 1 }
      } else { 1 }
      line(a, c, stroke: stroke)
      line(
        (a.at(0) + nx * g * sgn + ux * sh, a.at(1) + ny * g * sgn + uy * sh),
        (c.at(0) + nx * g * sgn - ux * sh, c.at(1) + ny * g * sgn - uy * sh),
        stroke: (paint: black, thickness: stroke.thickness, dash: "dashed"),
      )
    }
  }

  // ── GR-7.1 ionic dotted bonds between disconnected components ──────────────
  if cfg.ionic-bonds {
    let n = layout.atoms.len()
    let parent = range(n)
    let find(x) = {
      let r = x
      while parent.at(r) != r { r = parent.at(r) }
      r
    }
    for b in layout.bonds {
      let ra = find(b.from)
      let rb = find(b.to)
      if ra != rb { parent.at(ra) = rb }
    }
    let comp = range(n).map(find)
    let roots = comp.dedup()
    for i in range(roots.len()) {
      for j in range(i + 1, roots.len()) {
        let (ra, rb) = (roots.at(i), roots.at(j))
        let best = none
        for x in range(n) {
          if comp.at(x) != ra { continue }
          for y in range(n) {
            if comp.at(y) != rb { continue }
            let px = P(x)
            let py = P(y)
            let d = (px.at(0) - py.at(0)) * (px.at(0) - py.at(0)) + (px.at(1) - py.at(1)) * (px.at(1) - py.at(1))
            if best == none or d < best.at(0) { best = (d, x, y) }
          }
        }
        if best != none {
          line(P(best.at(1)), P(best.at(2)), stroke: (paint: gray, thickness: stroke.thickness, dash: "dotted"))
        }
      }
    }
  }

  // ── GR-9.4 variable attachment (bond into a ring centre) ───────────────────
  for (from, ring-ids) in cfg.variable-attach {
    let nn = ring-ids.len()
    let cx = ring-ids.map(i => P(i).at(0)).sum() / nn
    let cy = ring-ids.map(i => P(i).at(1)).sum() / nn
    let a = P(from)
    line(a, (a.at(0) + (cx - a.at(0)) * 0.7, a.at(1) + (cy - a.at(1)) * 0.7), stroke: stroke)
  }

  // ── GR-1.9 multi-centre (hapto / η) bond to a centroid ─────────────────────
  for (from, ids) in cfg.multi-centre {
    let nn = ids.len()
    let cx = ids.map(i => P(i).at(0)).sum() / nn
    let cy = ids.map(i => P(i).at(1)).sum() / nn
    line(P(from), (cx, cy), stroke: stroke)
  }

  // ── GR-6 aromatic delocalization circles ──────────────────────────────────
  if cfg.aromatic == "circle" {
    for ring in layout.at("aromatic_rings", default: ()) {
      circle((ring.center.x * s, ring.center.y * s), radius: ring.radius * s * cfg.aromatic-radius, stroke: stroke)
    }
  }

  // ── GR-5.6 partial charges (δ+/δ-) ─────────────────────────────────────────
  for (id, sign) in cfg.partial-charge {
    let a = P(int(id))
    let sg = if sign == "+" { "+" } else { "−" }
    content((a.at(0), a.at(1) + 0.42 * s), text(size: cfg.script-size * 1.1, [δ] + super(size: cfg.script-size, sg)))
  }

  // ── electron-pushing curly arrows (addressed by atom id) ───────────────────
  for arr in cfg.arrows {
    let a = P(arr.at(0))
    let b = P(arr.at(1))
    let side = if arr.len() > 2 { arr.at(2) } else { 1 }
    let dx = b.at(0) - a.at(0)
    let dy = b.at(1) - a.at(1)
    let len = calc.max(calc.sqrt(dx * dx + dy * dy), 1e-9)
    let (ux, uy) = (dx / len, dy / len)
    let (px, py) = (-uy * side, ux * side)
    let (off, bend, pad) = (0.3 * s, 0.55 * s, 0.16 * s)
    let p0 = (a.at(0) + ux * pad + px * off, a.at(1) + uy * pad + py * off)
    let p3 = (b.at(0) - ux * pad * 0.2 + px * off * 0.55, b.at(1) - uy * pad * 0.2 + py * off * 0.55)
    let c1 = (p0.at(0) + ux * len * 0.2 + px * bend, p0.at(1) + uy * len * 0.2 + py * bend)
    let c2 = (p3.at(0) + px * bend, p3.at(1) + py * bend)
    bezier(p0, p3, c1, c2, stroke: 0.7pt + cfg.arrow-paint, mark: (end: ">", scale: 0.85))
  }

  // ── GR-5.7 enclose a polyatomic ion in brackets ───────────────────────────
  if cfg.brackets != none {
    let xs = layout.atoms.map(a => a.pos.x * s)
    let ys = layout.atoms.map(a => a.pos.y * s)
    let pad = 0.45 * s
    let (x0, x1) = (calc.min(..xs) - pad, calc.max(..xs) + pad)
    let (y0, y1) = (calc.min(..ys) - pad, calc.max(..ys) + pad)
    let tick = 0.22 * s
    line((x0 + tick, y1), (x0, y1), (x0, y0), (x0 + tick, y0), stroke: stroke)
    line((x1 - tick, y1), (x1, y1), (x1, y0), (x1 - tick, y0), stroke: stroke)
    let q = cfg.brackets
    if q != 0 {
      let m = calc.abs(q)
      let sg = if q > 0 { "+" } else { "−" }
      content((x1 + 0.12 * s, y1), text(super(size: cfg.script-size, (if m > 1 { str(m) } else { "" }) + sg)), anchor: "south-west")
    }
  }
  })
}
