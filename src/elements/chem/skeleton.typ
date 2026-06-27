// Coordinate renderer that reuses alchemist's own drawing primitives.
//
// LayoutOut (from the WASM engine) is rendered by *placing each atom as a native
// alchemist fragment at its absolute coordinate* and then driving the existing
// `draw-link-decoration` for the bonds — no turtle geometry. This keeps a single
// drawing path (and therefore identical CeTZ anchors) shared with hand-drawn
// `skeletize` molecules; the engine only supplies coordinates.

#import "@preview/cetz:0.5.2"
#import "../../drawer.typ": default-ctx, draw-link-decoration
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

  // Plain (uncharged, no isotope/oxidation, non-vertical) condensed labels are
  // split into one sub-atom per element group so a bond connects to the right
  // element; everything else stays a single GR-formatted blob.
  let groups = if atom.charge == 0 and atom.isotope == none and ox == none and not vertical and not atom.skeletal {
    element-groups(atom.text)
  } else { none }
  if groups != none and groups.len() > 1 {
    let ordered = if reversed { groups.rev() } else { groups }
    let atoms = ordered.map(gp => ({ text(gp.sym); if gp.sub != "" { sub(gp.sub) } }, true))
    return (
      type: "fragment", name: name, atoms: atoms, colors: colors,
      links: (:), lewis: lewis, vertical: vertical, count: ordered.len(), empty: false,
    )
  }
  (
    type: "fragment", name: name, atoms: ((format-label(atom, reversed: reversed, oxidation: ox), true),),
    colors: colors, links: (:), lewis: lewis, vertical: vertical, count: 1, empty: false,
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

  {
    // place every atom at its absolute coordinate as a fragment
    cetz.draw.get-ctx(cetz-ctx => {
      // resolve a font-relative bond length (atom-sep = 3em) to canvas-unit floats
      let s = if type(cfg.scale) == length { convert-length(cetz-ctx, cfg.scale) } else { cfg.scale }
      let ctx = default-ctx
      ctx.config = alch
      for (i, a) in layout.atoms.enumerate() {
        ctx.last-anchor = (type: "coord", anchor: (a.pos.x * s, a.pos.y * s))
        let (c2, drawing) = fragment.draw-fragment-elements(frags.at(i), ctx)
        ctx = c2
        drawing
        // semantic anchor for DSL `:label` (mechanism arrows / cross-links)
        if a.at("label", default: none) != none {
          cetz.draw.anchor(name + "-" + a.label, (a.pos.x * s, a.pos.y * s))
        }
      }
    })
    // draw the plain bonds turtle-free via the shared link decorator
    {
      let ctx = default-ctx
      ctx.config = alch
      ctx.links = layout.bonds
        .filter(b => not is-special(cfg, b.from, b.to))
        .map(b => {
          let args = (:)
          if b.kind == "double" and (b.inner.x != 0 or b.inner.y != 0) {
            let f = atom(b.from).pos
            let t = atom(b.to).pos
            let z = (t.x - f.x) * b.inner.y - (t.y - f.y) * b.inner.x
            args.offset = if z > 0 { "left" } else { "right" }
          }
          let draw = link-fn.at(b.kind, default: links.single)
          let fsk = not is-labeled(atom(b.from), cfg)
          let tsk = not is-labeled(atom(b.to), cfg)
          (
            type: "link", hide: false,
            name: name + "-l" + str(b.from) + "-" + str(b.to),
            from-name: aname(b.from), to-name: aname(b.to),
            from: none, to: none,
            from-pos: if fsk { (name: aname(b.from)) } else { (name: aname(b.from), anchor: "mid") },
            to-pos: none,
            draw: draw(..args).first().draw,
            over: none, override: (:),
            ignore-from-margins: fsk, ignore-to-margins: tsk,
          )
        })
      // register placed fragments as hooks so the decorator can resolve them
      for (i, a) in layout.atoms.enumerate() {
        ctx.hooks.insert(aname(a.id), (
          type: "fragment", name: aname(a.id),
          count: frags.at(i).count, vertical: frags.at(i).vertical,
          empty: not is-labeled(a, cfg),
        ))
      }
      draw-link-decoration(ctx).at(1)
    }
  }
}
