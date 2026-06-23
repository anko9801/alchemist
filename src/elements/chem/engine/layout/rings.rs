//! Ring perception (SSSR) and ring-system placement.

use std::collections::{HashSet, VecDeque};
use std::f64::consts::PI;

use super::Pt;
use crate::graph::Graph;

// ── SSSR ─────────────────────────────────────────────────────────────────────

pub fn find_rings(g: &Graph) -> Vec<Vec<usize>> {
    let target = cycle_rank(g);
    if target == 0 {
        return Vec::new();
    }
    let words = (g.bonds.len() + 63) / 64;
    let mut seen: HashSet<Vec<u64>> = HashSet::new();
    let mut cand: Vec<(Vec<usize>, Vec<u64>)> = Vec::new();

    for (bi, b) in g.bonds.iter().enumerate() {
        if let Some(ring) = shortest_path_excluding(g, b.a, b.b, bi) {
            if ring.len() < 3 {
                continue;
            }
            let bits = ring_bits(g, &ring, words);
            if seen.insert(bits.clone()) {
                cand.push((ring, bits));
            }
        }
    }
    cand.sort_by(|a, b| a.0.len().cmp(&b.0.len()).then_with(|| a.1.cmp(&b.1)));

    let mut basis: Vec<Vec<u64>> = Vec::new();
    let mut rings = Vec::new();
    for (ring, bits) in cand {
        if add_independent(&mut basis, bits) {
            rings.push(ring);
            if rings.len() == target {
                break;
            }
        }
    }
    rings
}

fn cycle_rank(g: &Graph) -> usize {
    let n = g.n();
    if n == 0 {
        return 0;
    }
    let mut seen = vec![false; n];
    let mut comps = 0;
    for s in 0..n {
        if seen[s] {
            continue;
        }
        comps += 1;
        let mut st = vec![s];
        seen[s] = true;
        while let Some(u) = st.pop() {
            for &(v, _) in &g.adj[u] {
                if !seen[v] {
                    seen[v] = true;
                    st.push(v);
                }
            }
        }
    }
    g.bonds.len() + comps - n
}

fn shortest_path_excluding(g: &Graph, from: usize, to: usize, excl: usize) -> Option<Vec<usize>> {
    let n = g.n();
    let mut parent: Vec<Option<usize>> = vec![None; n];
    parent[from] = Some(from);
    let mut q = VecDeque::new();
    q.push_back(from);
    while let Some(u) = q.pop_front() {
        for &(v, bi) in &g.adj[u] {
            if bi == excl || parent[v].is_some() {
                continue;
            }
            parent[v] = Some(u);
            if v == to {
                let mut path = vec![to];
                let mut cur = to;
                while cur != from {
                    cur = parent[cur]?;
                    path.push(cur);
                }
                path.reverse();
                return Some(path);
            }
            q.push_back(v);
        }
    }
    None
}

fn ring_bits(g: &Graph, ring: &[usize], words: usize) -> Vec<u64> {
    let mut bits = vec![0u64; words];
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        if let Some(bi) = g.bond_between(a, b) {
            bits[bi / 64] |= 1 << (bi % 64);
        }
    }
    bits
}

fn add_independent(basis: &mut Vec<Vec<u64>>, bits: Vec<u64>) -> bool {
    let mut cand = bits;
    for e in basis.iter() {
        if let Some(p) = pivot(e) {
            if bit(&cand, p) {
                xor(&mut cand, e);
            }
        }
    }
    let Some(p) = pivot(&cand) else { return false };
    for e in basis.iter_mut() {
        if bit(e, p) {
            xor(e, &cand);
        }
    }
    basis.push(cand);
    basis.sort_by_key(|b| pivot(b).unwrap_or(usize::MAX));
    true
}

fn pivot(bits: &[u64]) -> Option<usize> {
    bits.iter()
        .enumerate()
        .find(|(_, w)| **w != 0)
        .map(|(i, w)| i * 64 + w.trailing_zeros() as usize)
}
fn bit(bits: &[u64], b: usize) -> bool {
    bits[b / 64] & (1 << (b % 64)) != 0
}
fn xor(a: &mut [u64], b: &[u64]) {
    for (x, y) in a.iter_mut().zip(b) {
        *x ^= y;
    }
}

pub fn tag_ring_membership(g: &mut Graph) {
    let rings = g.rings.clone();
    for (ri, ring) in rings.iter().enumerate() {
        for &a in ring {
            if !g.nodes[a].ring_ids.contains(&ri) {
                g.nodes[a].ring_ids.push(ri);
            }
        }
        for i in 0..ring.len() {
            let a = ring[i];
            let b = ring[(i + 1) % ring.len()];
            if let Some(bi) = g.bond_between(a, b) {
                if !g.bonds[bi].ring_ids.contains(&ri) {
                    g.bonds[bi].ring_ids.push(ri);
                }
            }
        }
    }
}

// ── ring geometry ────────────────────────────────────────────────────────────

fn ring_has_edge(ring: &[usize], a: usize, b: usize) -> bool {
    (0..ring.len()).any(|i| {
        let x = ring[i];
        let y = ring[(i + 1) % ring.len()];
        (x == a && y == b) || (x == b && y == a)
    })
}
pub fn largest_gap(sorted: &[f64]) -> (f64, f64) {
    let mut best_start = sorted[0];
    let mut best = 0.0;
    for i in 0..sorted.len() {
        let s = sorted[i];
        let e = if i + 1 < sorted.len() {
            sorted[i + 1]
        } else {
            sorted[0] + 2.0 * PI
        };
        if e - s > best {
            best = e - s;
            best_start = s;
        }
    }
    (best_start, best)
}

// ── ring-inner directions (double-bond side, GR-1.10) ────────────────────────

pub fn ring_inner_dirs(g: &Graph, coords: &[Pt]) -> Vec<(f64, f64)> {
    let mut dirs = vec![(0.0, 0.0); g.bonds.len()];
    for (bi, b) in g.bonds.iter().enumerate() {
        if b.kind.order() != 2 {
            continue;
        }
        // GR-1.10: ring double bonds offset toward the ring centroid;
        // asymmetric acyclic double bonds offset toward the more-substituted side;
        // symmetric / terminal double bonds stay centered (inner = 0).
        if let Some(ring) = best_ring_for_bond(g, b.a, b.b) {
            let n = ring.len() as f64;
            let cx = ring.iter().map(|&a| coords[a].0).sum::<f64>() / n;
            let cy = ring.iter().map(|&a| coords[a].1).sum::<f64>() / n;
            let mx = (coords[b.a].0 + coords[b.b].0) / 2.0;
            let my = (coords[b.a].1 + coords[b.b].1) / 2.0;
            let (vx, vy) = (cx - mx, cy - my);
            let l = (vx * vx + vy * vy).sqrt();
            if l > 1e-6 {
                dirs[bi] = (vx / l, vy / l);
            }
        } else {
            dirs[bi] = asymmetric_offset(g, coords, b.a, b.b);
        }
    }
    dirs
}

/// GR-1.10.1: pick the perpendicular side of an acyclic double bond that carries
/// more substituents. Returns (0,0) for symmetric or terminal double bonds.
fn asymmetric_offset(g: &Graph, coords: &[Pt], a: usize, b: usize) -> (f64, f64) {
    // cumulated double bonds (allene/cumulene): an sp centre carries two double
    // bonds; draw the bond centred (symmetric) rather than offset to one side.
    let doubles_at = |x: usize| {
        g.adj[x]
            .iter()
            .filter(|&&(_, bi)| g.bonds[bi].kind.order() == 2)
            .count()
    };
    if doubles_at(a) >= 2 || doubles_at(b) >= 2 {
        return (0.0, 0.0);
    }
    let pa = coords[a];
    let pb = coords[b];
    let (dx, dy) = (pb.0 - pa.0, pb.1 - pa.1);
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-6 {
        return (0.0, 0.0);
    }
    let (nx, ny) = (-dy / len, dx / len); // perpendicular
    let mx = (pa.0 + pb.0) / 2.0;
    let my = (pa.1 + pb.1) / 2.0;
    let mut net = 0.0;
    for &c in &[a, b] {
        for &(v, _) in &g.adj[c] {
            if v == a || v == b {
                continue;
            }
            let s = (coords[v].0 - mx) * nx + (coords[v].1 - my) * ny;
            net += s;
        }
    }
    if net > 1e-6 {
        (nx, ny)
    } else if net < -1e-6 {
        (-nx, -ny)
    } else {
        (0.0, 0.0)
    }
}

fn best_ring_for_bond<'a>(g: &'a Graph, a: usize, b: usize) -> Option<&'a Vec<usize>> {
    g.rings
        .iter()
        .filter(|r| ring_has_edge(r, a, b))
        .max_by(|x, y| {
            unsat(g, x)
                .cmp(&unsat(g, y))
                .then_with(|| y.len().cmp(&x.len()))
        })
}

fn unsat(g: &Graph, ring: &[usize]) -> usize {
    (0..ring.len())
        .filter_map(|i| g.bond_between(ring[i], ring[(i + 1) % ring.len()]))
        .filter(|&bi| g.bonds[bi].kind.order() == 2)
        .count()
}
