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
#import "labels.typ": format-label, element-color, element-groups
#import "decorations.typ": is-special

// A dative/coordination bond: a single line with an arrowhead (GR-1.7).
#let dative = links.build-link((length, ctx, cetz-ctx, args) => {
  import cetz.draw: *
  line((0, 0), (length, 0), stroke: args.at("stroke", default: ctx.config.single.stroke), mark: (end: ">"))
})

// GR double bond. `side` (+1/-1) puts the second line toward the ring centroid as
// a shortened inner line (skeletal convention); `side: 0` draws it symmetric (for
// non-ring double bonds, e.g. C=O). Matches the previous renderer's clean style
// rather than alchemist's full-length offset, which looks cramped in small rings.
#let gr-double(side) = links.build-link((length, ctx, cetz-ctx, args) => {
  import cetz.draw: *
  let stroke = args.at("stroke", default: ctx.config.double.stroke)
  let g = convert-length(cetz-ctx, args.at("gap", default: ctx.config.double.gap))
  if side == 0 {
    line((0, g / 2), (length, g / 2), stroke: stroke)
    line((0, -g / 2), (length, -g / 2), stroke: stroke)
  } else {
    let sh = g * 1.4 // shorten the inner line at both ends
    line((0, 0), (length, 0), stroke: stroke)
    line((sh, side * g), (length - sh, side * g), stroke: stroke)
  }
})

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

// Uniform gap kept between an atom label and the bonds touching it: each bond is
// trimmed to the edge of the measured label box plus this margin, so the same
// breathing room surrounds every atom (glyph-aware, ChemDraw ~1.6pt @ 10pt).
#let label-clearance = 0.18em

// Build a native alchemist fragment dict for one engine atom.
#let atom-fragment(atom, name, cfg) = {
  if not is-labeled(atom, cfg) {
    return (
      type: "fragment", name: name, atoms: ((none, true),), colors: none,
      links: (:), lewis: (), vertical: false, count: 1, empty: true,
    )
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
  let vertical = cfg.vertical.contains(atom.id) and not atom.skeletal
  let ox = cfg.oxidation.at(str(atom.id), default: none)
  let reversed = atom.at("label_dir", default: "center") == "left"

  // A vertical label (GR-2.1.7) is split into element groups so they stack
  // top-to-bottom; the heavy atom is first (= the connecting sub-atom).
  if vertical {
    let g = element-groups(atom.text)
    if g != none {
      let heavy = g.filter(x => x.sym == atom.element)
      let groups = heavy + g.filter(x => x.sym != atom.element)
      let atoms = groups.map(gp => ({ text(gp.sym); if gp.sub != "" { sub(gp.sub) } }, true))
      return (
        type: "fragment", name: name, atoms: atoms, colors: colors,
        links: (:), lewis: lewis, vertical: true, count: groups.len(), empty: false,
      )
    }
  }
  // A plain label (GR-2.1.5) is split into element groups placed left-to-right with
  // the heavy atom at the connecting end (first for a right label "CH2", last for a
  // reversed left label "H2C"), so the bond meets the *element symbol* on the atom
  // coordinate — not the centre of the whole "CH2" box.
  let g = element-groups(atom.text)
  if g != none and atom.charge == 0 and atom.isotope == none and ox == none {
    let heavy = g.filter(x => x.sym == atom.element)
    let rest = g.filter(x => x.sym != atom.element)
    let groups = if reversed { rest + heavy } else { heavy + rest }
    let atoms = groups.map(gp => ({ text(gp.sym); if gp.sub != "" { sub(gp.sub) } }, true))
    return (
      type: "fragment", name: name, atoms: atoms, colors: colors,
      links: (:), lewis: lewis, vertical: false, count: groups.len(), empty: false,
    )
  }
  (
    type: "fragment", name: name,
    atoms: ((format-label(atom, reversed: reversed, oxidation: ox), true),),
    colors: colors, links: (:), lewis: lewis, vertical: false, count: 1, empty: false,
  )
}

// Render a LayoutOut to CeTZ drawables (compose inside any canvas).
#let draw-skeleton-core(layout, name: "mol", config: (:)) = {
  let cfg = (
    // one engine bond-length unit maps to alchemist's `atom-sep`, so generated
    // molecules sit at the same scale as hand-drawn ones. Resolved to canvas-unit
    // floats at draw time (a plain number is used as-is). Styling (strokes,
    // margins, lewis, colours) is shared via alchemist's `default` config.
    scale: default.atom-sep,
    label-margin: 0.16em, // clearance beyond the label bbox (ChemDraw ~1.6pt @ 10pt)
    color: false,
    atom-colors: (:),
    lone-pairs: none,
    show-all-h: false,
    vertical: (),
    oxidation: (:),
    delocalize: (),
    bent: (),
    pseudo: (),
  ) + config
  let alch = default
  alch.fragment-color = none

  let aname(i) = name + "-a" + str(i)
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
      let s = if type(cfg.scale) == length { convert-length(cetz-ctx, cfg.scale) } else { cfg.scale }
      let ctx = (..default-ctx, config: alch)
      for (i, a) in layout.atoms.enumerate() {
        let frag = frags.at(i)
        let coord = (a.pos.x * s, a.pos.y * s)
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
      let s = if type(cfg.scale) == length { convert-length(cetz-ctx, cfg.scale) } else { cfg.scale }
      let margin = convert-length(cetz-ctx, label-clearance)
      let lctx = (..default-ctx, config: alch)
      // label box of atom `i` expressed relative to its coordinate (canvas units)
      let box-of(i) = {
        if not is-labeled(atom(i), cfg) { return none }
        let (_, nw) = cetz.coordinate.resolve(cetz-ctx, (name: aname(i), anchor: "north-west"))
        let (_, se) = cetz.coordinate.resolve(cetz-ctx, (name: aname(i), anchor: "south-east"))
        let cx = atom(i).pos.x * s
        let cy = atom(i).pos.y * s
        (x0: nw.at(0) - cx, x1: se.at(0) - cx, y0: se.at(1) - cy, y1: nw.at(1) - cy)
      }
      // distance from the coordinate to where the ray (ux,uy) leaves the box, + margin
      let trim(bx, ux, uy) = {
        if bx == none { return 0 }
        let big = 1e6
        let tx = if ux > 1e-6 { bx.x1 / ux } else if ux < -1e-6 { bx.x0 / ux } else { big }
        let ty = if uy > 1e-6 { bx.y1 / uy } else if uy < -1e-6 { bx.y0 / uy } else { big }
        calc.min(tx, ty) + margin
      }
      let boxes = layout.atoms.map(a => box-of(a.id))
      for b in layout.bonds {
        if is-special(cfg, b.from, b.to) { continue }
        let a = (atom(b.from).pos.x * s, atom(b.from).pos.y * s)
        let c = (atom(b.to).pos.x * s, atom(b.to).pos.y * s)
        let dx = c.at(0) - a.at(0)
        let dy = c.at(1) - a.at(1)
        let len = calc.max(calc.sqrt(dx * dx + dy * dy), 1e-9)
        let (ux, uy) = (dx / len, dy / len)
        let ta = trim(boxes.at(b.from), ux, uy)
        let tc = trim(boxes.at(b.to), -ux, -uy)
        let len2 = len - ta - tc
        if len2 <= 0 { continue }
        let a2 = (a.at(0) + ux * ta, a.at(1) + uy * ta)
        let drawer = if b.kind == "double" {
          let z = dx * b.inner.y - dy * b.inner.x
          let side = if b.inner.x == 0 and b.inner.y == 0 { 0 } else if z > 0 { 1 } else { -1 }
          gr-double(side)
        } else {
          link-fn.at(b.kind, default: links.single)
        }
        scope({
          set-origin(a2)
          rotate(calc.atan2(dx, dy))
          (drawer().first().draw)(len2, lctx, cetz-ctx)
        })
      }
    })
  }
}
