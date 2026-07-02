// Shared drawing geometry for the chem renderer.
//
// A line that touches an atom — whether it is a bond (skeleton.typ) or an
// overlay such as a hapto / variable-attachment line (decorations.typ) — must be
// pulled back so it clears that atom's label. That trimming logic lives here so
// both callers use exactly the same rule and no line ever runs into a glyph.

#import "@preview/cetz:0.5.2"

// The label box of the atom drawn as group `name`, expressed relative to its
// coordinate (cx, cy) in canvas units. Works for any placed atom (a bare vertex
// resolves to a near-zero box, so it is simply not trimmed beyond the margin).
#let label-box(cetz-ctx, name, cx, cy) = {
  let (_, nw) = cetz.coordinate.resolve(cetz-ctx, (name: name, anchor: "north-west"))
  let (_, se) = cetz.coordinate.resolve(cetz-ctx, (name: name, anchor: "south-east"))
  (x0: nw.at(0) - cx, x1: se.at(0) - cx, y0: se.at(1) - cy, y1: nw.at(1) - cy)
}

// Distance from the coordinate to where the ray (ux, uy) leaves `box`, plus a
// uniform margin. `box == none` (an atom with no drawn label) means no trim.
#let box-trim(box, ux, uy, margin) = {
  if box == none { return 0 }
  let big = 1e6
  let tx = if ux > 1e-6 { box.x1 / ux } else if ux < -1e-6 { box.x0 / ux } else { big }
  let ty = if uy > 1e-6 { box.y1 / uy } else if uy < -1e-6 { box.y0 / uy } else { big }
  calc.min(tx, ty) + margin
}
