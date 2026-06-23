//! Hand-written OpenSMILES parser -> Graph IR.
//!
//! Covers the organic subset + bracket atoms (isotope, charge, H-count,
//! chirality), bonds (`- = # :` and directional `/ \`), branches, ring closures
//! (single digit and `%nn`), the disconnection dot, and lowercase aromatic atoms
//! with Kekulization. Each SMILES atom becomes one graph vertex; carbons render
//! as bare skeletal vertices, heteroatoms as element + hydrogen labels.

use crate::graph::{standard_valence, BondKind, Graph, Node};

#[derive(Clone, Copy, PartialEq)]
enum BondTok {
    Default,
    Single,
    Double,
    Triple,
    Aromatic,
    Up,
    Down,
}

impl BondTok {
    fn kind(self) -> BondKind {
        match self {
            BondTok::Double => BondKind::Double,
            BondTok::Triple => BondKind::Triple,
            _ => BondKind::Single,
        }
    }
}

struct Builder {
    g: Graph,
    aromatic_atom: Vec<bool>,
    explicit_h: Vec<Option<u8>>,
    // temp bonds: (a, b, tok)
    bonds: Vec<(usize, usize, BondTok)>,
}

impl Builder {
    fn new() -> Self {
        Builder {
            g: Graph::default(),
            aromatic_atom: Vec::new(),
            explicit_h: Vec::new(),
            bonds: Vec::new(),
        }
    }

    fn add_atom(
        &mut self,
        element: String,
        aromatic: bool,
        charge: i8,
        isotope: Option<u16>,
        explicit_h: Option<u8>,
        chirality: i8,
    ) -> usize {
        let id = self.g.add_node(Node {
            text: None,
            element,
            charge,
            isotope,
            hcount: 0,
            label: None,
            aromatic,
            chirality,
            preceding: false,
            h_explicit: true,
            pos: None,
            ring_ids: Vec::new(),
        });
        self.aromatic_atom.push(aromatic);
        self.explicit_h.push(explicit_h);
        id
    }
}

struct P {
    c: Vec<char>,
    i: usize,
}

impl P {
    fn peek(&self) -> Option<char> {
        self.c.get(self.i).copied()
    }
    fn bump(&mut self) -> Option<char> {
        let c = self.peek();
        if c.is_some() {
            self.i += 1;
        }
        c
    }
}

const ORGANIC2: &[&str] = &["Cl", "Br"];
const ORGANIC1: &[char] = &['B', 'C', 'N', 'O', 'P', 'S', 'F', 'I'];
const AROMATIC1: &[char] = &['b', 'c', 'n', 'o', 'p', 's'];

fn parse_uint(p: &mut P) -> Option<u32> {
    let start = p.i;
    while matches!(p.peek(), Some(c) if c.is_ascii_digit()) {
        p.i += 1;
    }
    if p.i == start {
        None
    } else {
        p.c[start..p.i].iter().collect::<String>().parse().ok()
    }
}

fn parse_bracket(p: &mut P, b: &mut Builder) -> Result<usize, String> {
    // already consumed '['
    let isotope = parse_uint(p).map(|n| n as u16);
    // symbol
    let (element, aromatic) = match p.peek() {
        Some(c) if c.is_ascii_uppercase() => {
            let mut s = c.to_string();
            p.i += 1;
            if matches!(p.peek(), Some(l) if l.is_ascii_lowercase()) {
                s.push(p.bump().unwrap());
            }
            (s, false)
        }
        Some(c) if c.is_ascii_lowercase() => {
            p.i += 1;
            (c.to_ascii_uppercase().to_string(), true)
        }
        Some('*') => {
            p.i += 1;
            ("*".to_string(), false)
        }
        _ => return Err("bad bracket atom symbol".into()),
    };
    // chirality @ (anti) / @@ (clockwise)
    let mut chirality = 0i8;
    if p.peek() == Some('@') {
        p.i += 1;
        if p.peek() == Some('@') {
            p.i += 1;
            chirality = 1;
        } else {
            chirality = -1;
        }
        // skip a stereo-class tail like @TH1, but not the H-count
        while matches!(p.peek(), Some(c) if c.is_ascii_alphanumeric()) {
            if p.peek() == Some('H') {
                break;
            }
            p.i += 1;
        }
    }
    // hydrogens
    let mut hcount = 0u8;
    if p.peek() == Some('H') {
        p.i += 1;
        hcount = parse_uint(p).unwrap_or(1) as u8;
    }
    // charge
    let mut charge = 0i8;
    match p.peek() {
        Some('+') => {
            p.i += 1;
            if let Some(n) = parse_uint(p) {
                charge = n as i8;
            } else {
                charge = 1;
                while p.peek() == Some('+') {
                    charge += 1;
                    p.i += 1;
                }
            }
        }
        Some('-') => {
            p.i += 1;
            if let Some(n) = parse_uint(p) {
                charge = -(n as i8);
            } else {
                charge = -1;
                while p.peek() == Some('-') {
                    charge -= 1;
                    p.i += 1;
                }
            }
        }
        _ => {}
    }
    if p.peek() != Some(']') {
        return Err("unterminated bracket atom".into());
    }
    p.i += 1;
    Ok(b.add_atom(element, aromatic, charge, isotope, Some(hcount), chirality))
}

fn parse_organic(p: &mut P, b: &mut Builder) -> Option<usize> {
    // two-letter first
    if p.i + 1 < p.c.len() {
        let two: String = p.c[p.i..p.i + 2].iter().collect();
        if ORGANIC2.contains(&two.as_str()) {
            p.i += 2;
            return Some(b.add_atom(two, false, 0, None, None, 0));
        }
    }
    match p.peek() {
        Some(c) if ORGANIC1.contains(&c) => {
            p.i += 1;
            Some(b.add_atom(c.to_string(), false, 0, None, None, 0))
        }
        Some(c) if AROMATIC1.contains(&c) => {
            p.i += 1;
            Some(b.add_atom(c.to_ascii_uppercase().to_string(), true, 0, None, None, 0))
        }
        Some('*') => {
            p.i += 1;
            Some(b.add_atom("*".to_string(), false, 0, None, None, 0))
        }
        _ => None,
    }
}

pub fn parse(source: &str) -> Result<Graph, String> {
    let mut p = P {
        c: source.trim().chars().collect(),
        i: 0,
    };
    let mut b = Builder::new();
    let mut prev: Option<usize> = None;
    let mut pending = BondTok::Default;
    let mut stack: Vec<Option<usize>> = Vec::new();
    // ring closure digit -> (atom, bond tok)
    let mut rings: std::collections::HashMap<u32, (usize, BondTok)> = std::collections::HashMap::new();

    while let Some(ch) = p.peek() {
        match ch {
            '(' => {
                p.i += 1;
                stack.push(prev);
            }
            ')' => {
                p.i += 1;
                prev = stack.pop().ok_or("unbalanced ')' in SMILES")?;
            }
            '-' => {
                p.i += 1;
                pending = BondTok::Single;
            }
            '=' => {
                p.i += 1;
                pending = BondTok::Double;
            }
            '#' => {
                p.i += 1;
                pending = BondTok::Triple;
            }
            ':' => {
                p.i += 1;
                pending = BondTok::Aromatic;
            }
            '/' => {
                p.i += 1;
                pending = BondTok::Up;
            }
            '\\' => {
                p.i += 1;
                pending = BondTok::Down;
            }
            '.' => {
                p.i += 1;
                prev = None;
                pending = BondTok::Default;
            }
            '%' => {
                p.i += 1;
                let d1 = p.bump().and_then(|c| c.to_digit(10));
                let d2 = p.bump().and_then(|c| c.to_digit(10));
                let (Some(d1), Some(d2)) = (d1, d2) else {
                    return Err("bad %nn ring bond".into());
                };
                close_ring(&mut b, &mut rings, prev, pending, d1 * 10 + d2)?;
                pending = BondTok::Default;
            }
            c if c.is_ascii_digit() => {
                p.i += 1;
                close_ring(&mut b, &mut rings, prev, pending, c.to_digit(10).unwrap())?;
                pending = BondTok::Default;
            }
            '[' => {
                p.i += 1;
                let id = parse_bracket(&mut p, &mut b)?;
                link(&mut b, &mut prev, pending, id);
                pending = BondTok::Default;
            }
            c if c.is_alphabetic() || c == '*' => {
                let Some(id) = parse_organic(&mut p, &mut b) else {
                    return Err(format!("unexpected character {:?}", c));
                };
                link(&mut b, &mut prev, pending, id);
                pending = BondTok::Default;
            }
            c if c.is_whitespace() => {
                p.i += 1;
            }
            c => return Err(format!("unexpected character {:?}", c)),
        }
    }

    if !stack.is_empty() {
        return Err("unclosed '(' in SMILES".into());
    }
    if pending != BondTok::Default {
        return Err("SMILES ends with a dangling bond".into());
    }
    if !rings.is_empty() {
        return Err("unclosed ring bond in SMILES".into());
    }

    finalize(b)
}

fn link(b: &mut Builder, prev: &mut Option<usize>, pending: BondTok, id: usize) {
    if let Some(p) = *prev {
        b.g.nodes[id].preceding = true;
        let tok = if pending == BondTok::Default {
            if b.aromatic_atom[p] && b.aromatic_atom[id] {
                BondTok::Aromatic
            } else {
                BondTok::Single
            }
        } else {
            pending
        };
        b.bonds.push((p, id, tok));
    }
    *prev = Some(id);
}

fn close_ring(
    b: &mut Builder,
    rings: &mut std::collections::HashMap<u32, (usize, BondTok)>,
    prev: Option<usize>,
    pending: BondTok,
    n: u32,
) -> Result<(), String> {
    let cur = prev.ok_or("ring bond before any atom")?;
    if let Some((open, open_tok)) = rings.remove(&n) {
        let tok = if pending != BondTok::Default {
            pending
        } else if open_tok != BondTok::Default {
            open_tok
        } else if b.aromatic_atom[open] && b.aromatic_atom[cur] {
            BondTok::Aromatic
        } else {
            BondTok::Single
        };
        b.bonds.push((open, cur, tok));
    } else {
        rings.insert(n, (cur, pending));
    }
    Ok(())
}

fn finalize(mut b: Builder) -> Result<Graph, String> {
    let n = b.g.n();
    // Kekulize aromatic bonds: assign alternating doubles.
    let mut order = vec![BondTok::Single; b.bonds.len()];
    for (i, &(_, _, tok)) in b.bonds.iter().enumerate() {
        order[i] = tok;
    }
    kekulize(&b, &mut order);

    // Create real graph bonds, preserving / \ direction for cis/trans.
    let dirs: Vec<i8> = b
        .bonds
        .iter()
        .map(|&(_, _, t)| match t {
            BondTok::Up => 1,
            BondTok::Down => -1,
            _ => 0,
        })
        .collect();
    for (i, &(a, c, _)) in b.bonds.iter().enumerate() {
        let idx = b.g.add_bond(a, c, order[i].kind());
        b.g.bonds[idx].direction = dirs[i];
    }

    // Implicit H + label text per atom.
    for i in 0..n {
        let elem = b.g.nodes[i].element.clone();
        let charge = b.g.nodes[i].charge;
        let bond_sum: i16 = b.g.adj[i]
            .iter()
            .map(|&(_, bi)| b.g.bonds[bi].kind.order() as i16)
            .sum();
        let h = match b.explicit_h[i] {
            Some(h) => h,
            None => match standard_valence(&elem) {
                Some(v) => (v + charge_adjust(&elem, charge) - bond_sum).max(0) as u8,
                None => 0,
            },
        };
        b.g.nodes[i].hcount = h;
        // label: carbons stay skeletal unless charged/isotope; heteroatoms labelled
        let skeletal = elem == "C" && charge == 0 && b.g.nodes[i].isotope.is_none();
        if !skeletal {
            b.g.nodes[i].text = Some(make_label(&elem, h));
        }
    }

    Ok(b.g)
}

/// Valence adjustment for charged heteroatoms (e.g. N+ has valence 4, O- valence 1).
fn charge_adjust(elem: &str, charge: i8) -> i16 {
    match elem {
        "N" | "P" => charge as i16, // N+ -> 4 bonds, so +1 to available
        "O" | "S" => charge as i16,
        "C" => -(charge.abs() as i16),
        _ => 0,
    }
}

fn make_label(elem: &str, h: u8) -> String {
    if h == 0 {
        elem.to_string()
    } else if h == 1 {
        format!("{elem}H")
    } else {
        format!("{elem}H{h}")
    }
}

/// Greedy + augmenting Kekulization: aromatic atoms that need a double bond are
/// matched along aromatic edges.
fn kekulize(b: &Builder, order: &mut [BondTok]) {
    let n = b.g.n();
    // which atoms need a double bond in the pi system
    let mut needs = vec![false; n];
    for i in 0..n {
        if !b.aromatic_atom[i] {
            continue;
        }
        // an atom already carrying a double/triple bond (e.g. exocyclic c(=O))
        // does not need a ring double bond.
        let has_multiple = b
            .bonds
            .iter()
            .any(|&(a, c, t)| (a == i || c == i) && matches!(t, BondTok::Double | BondTok::Triple));
        if has_multiple {
            continue;
        }
        let elem = &b.g.nodes[i].element;
        let charge = b.g.nodes[i].charge;
        // total connectivity (any bond + implicit H) distinguishes pyridine-type
        // from pyrrole-type aromatic atoms.
        let total_deg = b.bonds.iter().filter(|&&(a, c, _)| a == i || c == i).count()
            + b.g.nodes[i].hcount as usize;
        needs[i] = match elem.as_str() {
            "C" => charge == 0,
            // Pyridine-type N/P carries only its two ring bonds (no H/substituent)
            // and takes a ring double bond; pyrrole-type N (3-connected via an H or
            // substituent) donates its lone pair to the ring and stays all-single.
            "N" | "P" => charge > 0 || total_deg <= 2,
            "B" => true,
            _ => false, // O, S donate a lone pair
        };
    }
    // aromatic edges
    let arom_edges: Vec<usize> = b
        .bonds
        .iter()
        .enumerate()
        .filter(|(_, &(_, _, t))| t == BondTok::Aromatic)
        .map(|(i, _)| i)
        .collect();

    let mut matched = vec![false; n];
    // greedy by ascending available degree
    let mut atoms: Vec<usize> = (0..n).filter(|&i| needs[i]).collect();
    atoms.sort_by_key(|&i| {
        arom_edges
            .iter()
            .filter(|&&e| b.bonds[e].0 == i || b.bonds[e].1 == i)
            .count()
    });
    for &i in &atoms {
        if matched[i] {
            continue;
        }
        // find an aromatic edge to an unmatched needing neighbor
        for &e in &arom_edges {
            let (a, c, _) = b.bonds[e];
            let j = if a == i {
                c
            } else if c == i {
                a
            } else {
                continue;
            };
            if needs[j] && !matched[j] {
                order[e] = BondTok::Double;
                matched[i] = true;
                matched[j] = true;
                break;
            }
        }
    }
    // remaining aromatic edges stay single
    for &e in &arom_edges {
        if order[e] == BondTok::Aromatic {
            order[e] = BondTok::Single;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(s: &str) -> Graph {
        parse(s).unwrap()
    }

    #[test]
    fn ethanol() {
        let m = g("CCO");
        assert_eq!(m.nodes.len(), 3);
        assert_eq!(m.bonds.len(), 2);
        assert_eq!(m.nodes[2].element, "O");
        assert_eq!(m.nodes[2].text.as_deref(), Some("OH")); // O + 1 implicit H
    }

    #[test]
    fn benzene_kekulized() {
        let m = g("c1ccccc1");
        assert_eq!(m.nodes.len(), 6);
        assert_eq!(m.bonds.len(), 6);
        assert_eq!(
            m.bonds.iter().filter(|b| b.kind == BondKind::Double).count(),
            3,
            "benzene kekulizes to 3 double bonds"
        );
    }

    #[test]
    fn bracket_charge() {
        let m = g("[NH4+]");
        assert_eq!(m.nodes[0].charge, 1);
        assert_eq!(m.nodes[0].text.as_deref(), Some("NH4"));
    }

    #[test]
    fn ring_closure() {
        let m = g("C1CCCCC1"); // cyclohexane
        assert_eq!(m.nodes.len(), 6);
        assert_eq!(m.bonds.len(), 6);
    }

    #[test]
    fn branches_and_double() {
        let m = g("CC(=O)O"); // acetic acid
        assert_eq!(m.nodes.len(), 4);
        assert!(m.bonds.iter().any(|b| b.kind == BondKind::Double));
    }
}
