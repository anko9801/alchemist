// Coordinate renderer that reuses alchemist's own drawing primitives.
//
// LayoutOut (from the WASM engine) is rendered by placing each atom as a native
// alchemist fragment at its absolute coordinate and drawing the bonds with
// alchemist's own link draw functions, centre-to-centre with a uniform gap at
// labelled atoms — no turtle geometry. Reuses alchemist's fragment/link/lewis
// drawing and config; the engine only supplies coordinates.

#import "@preview/cetz:0.5.2"
#import "../../drawer.typ": default-ctx
#import "../../drawer/fragment.typ" as fragment
#import "../../default.typ": default
#import "../../utils/utils.typ": convert-length
#import "../../elements/links.typ" as links
#import "../../elements/lewis.typ": lewis-double, lewis-line, lewis-single
#import "labels.typ": format-label, element-color, element-groups, group-math
#import "style.typ": chem-defaults
#import "geometry.typ": atom-anchor, atom-pos, box-trim, label-box, resolve-scale

// membership test for a bond {from,to} in an unordered atom-pair list
#let pair-in(list, f, t) = list.any(p => (p.at(0) == f and p.at(1) == t) or (p.at(0) == t and p.at(1) == f))

// A dative/coordination bond: a single line with an arrowhead (GR-1.7).
#let dative = links.build-link((length, ctx, cetz-ctx, args) => {
  import cetz.draw: *
  line((0, 0), (length, 0), stroke: args.at("stroke", default: ctx.config.single.stroke), mark: (end: ">"))
})

// A double bond maps onto alchemist's own `double` link via its `offset` option:
// "center" for a symmetric (acyclic) bond, "left"/"right" to put the shortened
// inner line on the ring-centroid side. Returns the offset for a given engine
// `inner` vector expressed in the bond-local frame (z = cross(bond, inner)).
#let double-offset(z, zero) = if zero { "center" } else if z > 0 { "left" } else { "right" }

// engine bond `kind` -> alchemist link function. Wedges are swapped right<->left:
// the engine names a wedge narrow-at-`from` (the stereocentre), while alchemist's
// `cram-*-right` draws wide-at-`from`, so the opposite-handed variant gives the
// chemically correct narrow-at-stereocentre triangle.
#let link-fn = (
  "single": links.single,
  "double": links.double,
  "triple": links.triple,
  "dative": dative,
  "cram-filled-right": links.cram-filled-left,
  "cram-filled-left": links.cram-filled-right,
  "cram-dashed-right": links.cram-dashed-left,
  "cram-dashed-left": links.cram-dashed-right,
  "cram-hollow-right": links.cram-hollow-left,
  "cram-hollow-left": links.cram-hollow-right,
)

// is this atom rendered as a glyph (label) rather than a bare skeletal vertex?
#let is-labeled(atom, cfg) = {
  not atom.skeletal or (cfg.show-all-h and atom.at("implicit_h", default: 0) > 0)
}

// an alchemist fragment dict; `count`/`empty` follow from the sub-atom list
#let make-fragment(name, atoms, colors: none, lewis: (), vertical: false, empty: false) = (
  type: "fragment", name: name, atoms: atoms, colors: colors,
  links: (:), lewis: lewis, vertical: vertical, count: atoms.len(), empty: empty,
)

// Build a native alchemist fragment dict for one engine atom.
#let atom-fragment(atom, name, cfg) = {
  if not is-labeled(atom, cfg) {
    return make-fragment(name, ((none, true),), empty: true)
  }
  // lone pairs / radicals as alchemist lewis elements, from engine directions
  let lewis = ()
  if cfg.lone-pairs != none {
    let lp = if cfg.lone-pairs == "lines" { lewis-line } else { lewis-double }
    for d in atom.lone_pair_dirs { lewis.push(lp(angle: calc.atan2(d.x, d.y))) }
  }
  for d in atom.at("radical_dirs", default: ()) {
    lewis.push(lewis-single(angle: calc.atan2(d.x, d.y), offset: "center"))
  }
  let clr = element-color(atom.element, cfg.color, cfg.atom-colors)
  let colors = if clr != black { clr } else { none }
  let ox = cfg.oxidation.at(str(atom.id), default: none)
  let reversed = atom.at("label_dir", default: "center") == "left"
  let groups = element-groups(atom.text)
  let atoms(order) = order.map(gp => (group-math(gp.sym, gp.sub), true))

  if groups != none {
    let heavy = groups.filter(x => x.sym == atom.element)
    let rest = groups.filter(x => x.sym != atom.element)
    // GR-2.1.7 vertical: element groups stacked, heavy atom first (the connector)
    if cfg.vertical.contains(atom.id) and not atom.skeletal {
      return make-fragment(name, atoms(heavy + rest), colors: colors, lewis: lewis, vertical: true)
    }
    // GR-2.1.5 plain: groups left-to-right, heavy atom at the connecting end (first
    // for a right label "CH2", last for a reversed left label "H2C"), so the bond
    // meets the *element symbol* on the coordinate — not the whole-"CH2"-box centre.
    if atom.charge == 0 and atom.isotope == none and ox == none {
      return make-fragment(name, atoms(if reversed { rest + heavy } else { heavy + rest }), colors: colors, lewis: lewis)
    }
  }
  // charged / isotopic / oxidation labels: one GR-formatted blob
  make-fragment(name, ((format-label(atom, reversed: reversed, oxidation: ox), true),), colors: colors, lewis: lewis)
}

// Render a LayoutOut to CeTZ drawables (compose inside any canvas).
#let draw-skeleton-core(layout, name: "mol", config: (:)) = {
  let cfg = chem-defaults + config
  // one shared alchemist config so every alchemist link (single/double/triple/
  // cram/dative) strokes at the same weight as the manual special-bond styles
  let alch = default
  alch.fragment-color = none
  alch.single.stroke = cfg.stroke
  alch.double.stroke = cfg.stroke
  alch.triple.stroke = cfg.stroke

  let aname(i) = atom-anchor(name, i)
  let atom(i) = layout.atoms.at(i)
  let frags = layout.atoms.map(a => atom-fragment(a, aname(a.id), cfg))
  // sub-atom index of the heavy atom every bond should connect to (the element
  // symbol, not a trailing H): last group for reversed labels, else the first.
  let conn-idx(i) = if atom(i).at("label_dir", default: "center") == "left" {
    frags.at(i).count - 1
  } else { 0 }

  {
    // place every atom at its absolute coordinate as a fragment, with its heavy
    // atom (the connecting sub-atom) centred on the coordinate. Drawn first so the
    // bonds below can measure each label's box and trim themselves to clear it.
    cetz.draw.get-ctx(cetz-ctx => {
      import cetz.draw: *
      let s = resolve-scale(cetz-ctx, cfg.scale)
      let ctx = (..default-ctx, config: alch)
      for (i, a) in layout.atoms.enumerate() {
        let frag = frags.at(i)
        let coord = atom-pos(a, s)
        let conn = conn-idx(a.id)
        group(name: aname(a.id), anchor: "conn", {
          anchor("default", coord)
          fragment.draw-fragment-text(ctx, frag, coord)
          anchor("conn", (name: str(conn), anchor: "mid"))
        })
        fragment.draw-fragment-lewis(ctx, aname(a.id), frag.count, frag.lewis)
        // semantic anchor for DSL `:label` (mechanism arrows / cross-links)
        if a.at("label", default: none) != none {
          anchor(name + "-" + a.label, coord)
        }
      }
    })
    // draw the bonds, shortening the *bond length* itself at each labelled end so
    // it stops a uniform gap short of the label box (ChemDraw / old-renderer way).
    // Because the length is trimmed before the link is drawn, both lines of a
    // double bond come out equal, every end is a clean cap perpendicular to the
    // bond, and no white box or halo is painted over anything.
    cetz.draw.get-ctx(cetz-ctx => {
      import cetz.draw: *
      let s = resolve-scale(cetz-ctx, cfg.scale)
      let margin = convert-length(cetz-ctx, cfg.label-clearance)
      let lctx = (..default-ctx, config: alch)
      // label box of each atom relative to its coordinate (none for bare vertices)
      let boxes = layout.atoms.map(a => if is-labeled(a, cfg) {
        let (x, y) = atom-pos(a, s)
        label-box(cetz-ctx, aname(a.id), x, y)
      } else { none })
      let trim(bx, ux, uy) = box-trim(bx, ux, uy, margin)
      let stroke = cfg.stroke
      let dashed = (paint: black, thickness: stroke.thickness, dash: "dashed")
      // delocalized inner line reuses the *same* offset as a real double bond
      let g = convert-length(cetz-ctx, alch.double.gap)
      let bent-off = cfg.bent-kink * s
      for b in layout.bonds {
        let a = atom-pos(atom(b.from), s)
        let c = atom-pos(atom(b.to), s)
        let dx = c.at(0) - a.at(0)
        let dy = c.at(1) - a.at(1)
        let len = calc.max(calc.sqrt(dx * dx + dy * dy), 1e-9)
        let (ux, uy) = (dx / len, dy / len)
        let ta = trim(boxes.at(b.from), ux, uy)
        let tc = trim(boxes.at(b.to), -ux, -uy)
        let len2 = len - ta - tc
        if len2 <= 0 { continue }
        let a2 = (a.at(0) + ux * ta, a.at(1) + uy * ta)
        // z = cross(bond, inner): the side the offset/inner line goes to
        let z = dx * b.inner.y - dy * b.inner.x
        let side = if z > 0 { 1 } else { -1 }
        // everything is drawn in the bond-local frame: origin at the trimmed start,
        // +x along the bond over `len2`, +y perpendicular (left).
        scope({
          set-origin(a2)
          rotate(calc.atan2(dx, dy))
          if cfg.aromatic == "circle" and b.at("aromatic", default: false) {
            // GR-6 circle mode: aromatic ring bonds are plain single lines
            line((0, 0), (len2, 0), stroke: stroke)
          } else if pair-in(cfg.delocalize, b.from, b.to) {
            // GR-6 partial (delocalized) bond: solid line + dashed inner line
            line((0, 0), (len2, 0), stroke: stroke)
            line((g, side * g), (len2 - g, side * g), stroke: dashed)
          } else if pair-in(cfg.bent, b.from, b.to) {
            // GR-1.5 bent bond: a perpendicular kink at the midpoint
            line((0, 0), (len2 / 2, bent-off), (len2, 0), stroke: stroke)
          } else if pair-in(cfg.pseudo, b.from, b.to) {
            // GR-12 pseudobond: a wavy connector
            let waves = int(calc.max(3, calc.round(len2 / (cfg.wave-period * s))))
            let amp = cfg.wave-amplitude * s
            let steps = waves * 6
            let pts = range(steps + 1).map(k => {
              let t = k / steps
              (t * len2, calc.sin(t * waves * calc.pi) * amp)
            })
            line(..pts, stroke: stroke)
          } else {
            // ordinary bond via alchemist's own link functions
            let (link, override) = if b.kind == "double" {
              let zero = b.inner.x == 0 and b.inner.y == 0
              (links.double, (offset: double-offset(z, zero)))
            } else {
              (link-fn.at(b.kind, default: links.single), (:))
            }
            (link().first().draw)(len2, lctx, cetz-ctx, override: override)
          }
        })
      }
    })
  }
}
