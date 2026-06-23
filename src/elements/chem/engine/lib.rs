//! Alchemist molecule engine — Rust/WASM core.
//!
//! Pipeline: source string (DSL or SMILES) -> Graph IR -> layout -> LayoutOut JSON.
//! The Typst side maps the coordinates to native alchemist elements.

mod front;
mod graph;
mod layout;
mod output;

use serde::Deserialize;
use wasm_minimal_protocol::*;

initiate_protocol!();

#[derive(Deserialize)]
struct Input {
    #[serde(default = "default_format")]
    format: String,
    source: String,
    #[serde(default)]
    options: Options,
}

fn default_format() -> String {
    "dsl".to_string()
}

#[derive(Deserialize, Default)]
struct Options {
    /// "iupac" (default) | "as-written" | "none"
    #[serde(default)]
    orientation: Option<String>,
    /// extra rotation in degrees applied after orientation
    #[serde(default)]
    rotation: f64,
}

/// Stage 1 of the coordgen pipeline: parse the source and emit the molecular
/// graph as the whitespace-separated integer format the coordgen plugin reads:
///   nAtoms  Z...  nBonds  (a b order)...
#[wasm_func]
pub fn cg_input(input: &[u8]) -> Result<Vec<u8>, String> {
    use std::fmt::Write;
    let input: Input =
        serde_json::from_slice(input).map_err(|e| format!("invalid input json: {e}"))?;
    let g = front::parse(&input.format, &input.source)?;
    let mut s = String::new();
    let _ = write!(s, "{}\n", g.n());
    for node in &g.nodes {
        let _ = write!(s, "{} ", graph::atomic_number(&node.element));
    }
    let _ = write!(s, "\n{}\n", g.bonds.len());
    for b in &g.bonds {
        let _ = write!(s, "{} {} {}\n", b.a, b.b, b.kind.order());
    }
    Ok(s.into_bytes())
}

/// Stage 2 of the coordgen pipeline: given the original source (to re-derive the
/// graph deterministically) and the coordinates produced by the coordgen plugin
/// ("x y" per atom line), apply IUPAC orientation + stereo and build LayoutOut.
#[wasm_func]
pub fn finish(meta: &[u8], coords: &[u8]) -> Result<Vec<u8>, String> {
    let input: Input =
        serde_json::from_slice(meta).map_err(|e| format!("invalid input json: {e}"))?;
    let mut graph = front::parse(&input.format, &input.source)?;
    let text = std::str::from_utf8(coords).map_err(|e| format!("coords not utf8: {e}"))?;
    let mut pts: Vec<(f64, f64)> = Vec::with_capacity(graph.n());
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut it = line.split_whitespace();
        let x: f64 = it.next().and_then(|v| v.parse().ok()).ok_or("bad coord x")?;
        let y: f64 = it.next().and_then(|v| v.parse().ok()).ok_or("bad coord y")?;
        pts.push((x, y));
    }
    if pts.len() != graph.n() {
        return Err(format!("coord count {} != atom count {}", pts.len(), graph.n()));
    }
    let iupac = input.options.orientation.as_deref().unwrap_or("iupac") == "iupac";
    layout::from_coords(&mut graph, &pts, iupac, input.options.rotation.to_radians());
    Ok(serde_json::to_vec(&output::build(&graph)).map_err(|e| format!("serialize failed: {e}"))?)
}

