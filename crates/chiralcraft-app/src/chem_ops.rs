//! Advanced ChiralCraft interactive canvas & chemistry operations:
//! - Ring-interior double bond inward normal detection (ChiralCraft 1:1 visual bond rendering)
//! - Freehand polygon Lasso hit-testing
//! - Mass Fragmentation bond-cleavage exact mass & formula computation
//! - Structure Perspective 3D tilt projection
//! - Clean Up Reaction & Automatic Reaction Atom Mapping
//! - Interactive directory browser helper

use std::collections::{HashMap, VecDeque};
use std::fs;
use std::path::Path;

use chiralcraft_chem::{ChemAtom, ChemBond, ChemMolecule};
use chiralcraft_doc::DocumentScene;

use crate::icons::format_chemical_formula;

/// For each bond index belonging to a 3..=8-membered ring, returns a unit vector `(nx, ny)`
/// pointing inward toward the smallest ring's centroid so double bonds render with
/// ChiralCraft's authentic shortened inner parallel line.
pub fn ring_bond_inward_normals(mol: &ChemMolecule) -> HashMap<usize, (f64, f64)> {
    let n = mol.atoms.len();
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.begin_atom < n && b.end_atom < n {
            adj[b.begin_atom].push((b.end_atom, bi));
            adj[b.end_atom].push((b.begin_atom, bi));
        }
    }

    let mut normals = HashMap::new();
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.begin_atom >= n || b.end_atom >= n {
            continue;
        }
        // BFS shortest path from begin_atom to end_atom excluding bond `bi` (max depth 7 => ring <= 8)
        let mut parent: Vec<Option<usize>> = vec![None; n];
        let mut depth: Vec<usize> = vec![usize::MAX; n];
        let mut q = VecDeque::new();
        depth[b.begin_atom] = 0;
        q.push_back(b.begin_atom);

        let mut found = false;
        while let Some(u) = q.pop_front() {
            if depth[u] >= 7 {
                continue;
            }
            for &(v, ebi) in &adj[u] {
                if ebi == bi {
                    continue;
                }
                if depth[v] == usize::MAX {
                    depth[v] = depth[u] + 1;
                    parent[v] = Some(u);
                    if v == b.end_atom {
                        found = true;
                        break;
                    }
                    q.push_back(v);
                }
            }
            if found {
                break;
            }
        }
        if !found {
            continue;
        }
        let mut cycle = vec![b.end_atom];
        let mut cur = b.end_atom;
        while let Some(p) = parent[cur] {
            cycle.push(p);
            if p == b.begin_atom {
                break;
            }
            cur = p;
        }
        if cycle.len() < 3 || cycle.len() > 8 {
            continue;
        }
        let cx = cycle.iter().map(|&i| mol.atoms[i].x).sum::<f64>() / (cycle.len() as f64);
        let cy = cycle.iter().map(|&i| mol.atoms[i].y).sum::<f64>() / (cycle.len() as f64);
        let a1 = &mol.atoms[b.begin_atom];
        let a2 = &mol.atoms[b.end_atom];
        let dx = a2.x - a1.x;
        let dy = a2.y - a1.y;
        let len = dx.hypot(dy).max(1e-6);
        let mut nx = -dy / len;
        let mut ny = dx / len;
        let mx = 0.5 * (a1.x + a2.x);
        let my = 0.5 * (a1.y + a2.y);
        if (cx - mx) * nx + (cy - my) * ny < 0.0 {
            nx = -nx;
            ny = -ny;
        }
        normals.insert(bi, (nx, ny));
    }
    normals
}

/// Ray-casting point-in-polygon check for genuine freehand Lasso selection.
pub fn point_in_polygon(px: f64, py: f64, poly: &[(f64, f64)]) -> bool {
    if poly.len() < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = poly.len() - 1;
    for i in 0..poly.len() {
        let (xi, yi) = poly[i];
        let (xj, yj) = poly[j];
        let intersects = ((yi > py) != (yj > py))
            && (px < (xj - xi) * (py - yi) / ((yj - yi).abs().max(1e-12) * (yj - yi).signum()) + xi);
        if intersects {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn segments_intersect(a: (f64, f64), b: (f64, f64), c: (f64, f64), d: (f64, f64)) -> bool {
    fn orient(p: (f64, f64), q: (f64, f64), r: (f64, f64)) -> f64 {
        (q.0 - p.0) * (r.1 - p.1) - (q.1 - p.1) * (r.0 - p.0)
    }
    let o1 = orient(a, b, c);
    let o2 = orient(a, b, d);
    let o3 = orient(c, d, a);
    let o4 = orient(c, d, b);
    (o1 * o2 < 0.0) && (o3 * o4 < 0.0)
}

/// ChiralCraft `Mass Fragmentation` tool:
/// Cuts all bonds crossed by `(p1, p2)` and computes exact formula & monoisotopic exact mass
/// for each resulting radical fragment (matching ChiralCraft's homolytic bond cleavage report).
pub fn fragment_molecule_by_line(
    mol: &ChemMolecule,
    p1: (f64, f64),
    p2: (f64, f64),
) -> Vec<(String, f64, (f64, f64))> {
    let mut cut_bonds = Vec::new();
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.begin_atom < mol.atoms.len() && b.end_atom < mol.atoms.len() {
            let a1 = (mol.atoms[b.begin_atom].x, mol.atoms[b.begin_atom].y);
            let a2 = (mol.atoms[b.end_atom].x, mol.atoms[b.end_atom].y);
            if segments_intersect(p1, p2, a1, a2) {
                cut_bonds.push(bi);
            }
        }
    }
    if cut_bonds.is_empty() {
        return Vec::new();
    }

    let n = mol.atoms.len();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        if cut_bonds.contains(&bi) || b.begin_atom >= n || b.end_atom >= n {
            continue;
        }
        adj[b.begin_atom].push(b.end_atom);
        adj[b.end_atom].push(b.begin_atom);
    }

    let mut visited = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut comp = Vec::new();
        let mut q = VecDeque::new();
        visited[start] = true;
        q.push_back(start);
        while let Some(u) = q.pop_front() {
            comp.push(u);
            for &v in &adj[u] {
                if !visited[v] {
                    visited[v] = true;
                    q.push_back(v);
                }
            }
        }
        // Build sub-molecule with radical on cleaved atoms (homolytic cleavage, no extra H gained!)
        let mut idx_map = HashMap::new();
        let mut sub = ChemMolecule::default();
        for (new_i, &old_i) in comp.iter().enumerate() {
            idx_map.insert(old_i, new_i);
            let mut a: ChemAtom = mol.atoms[old_i].clone();
            let cleaved_orders: i32 = cut_bonds
                .iter()
                .filter(|&&bi| mol.bonds[bi].begin_atom == old_i || mol.bonds[bi].end_atom == old_i)
                .map(|&bi| mol.bonds[bi].order.clamp(1, 3))
                .sum();
            if cleaved_orders > 0 {
                a.radical = a.radical.saturating_add(cleaved_orders).min(2);
            }
            sub.atoms.push(a);
        }
        for (bi, b) in mol.bonds.iter().enumerate() {
            if cut_bonds.contains(&bi) {
                continue;
            }
            if let (Some(&ni), Some(&nj)) = (idx_map.get(&b.begin_atom), idx_map.get(&b.end_atom)) {
                let mut nb: ChemBond = b.clone();
                nb.begin_atom = ni;
                nb.end_atom = nj;
                sub.bonds.push(nb);
            }
        }
        let mut fp = sub.compute_formula_properties();
        let formula = format_chemical_formula(&fp.get_formula());
        let em = fp.get_exact_mass().val;
        let cx = sub.atoms.iter().map(|a| a.x).sum::<f64>() / (sub.atoms.len().max(1) as f64);
        let cy = sub.atoms.iter().map(|a| a.y).sum::<f64>() / (sub.atoms.len().max(1) as f64);
        out.push((formula, em, (cx, cy)));
    }
    out
}

/// ChiralCraft `Structure Perspective` tool: applies 3D yaw/pitch foreshortening around centroid.
pub fn apply_perspective_tilt(mol: &mut ChemMolecule, dx: f64, dy: f64) {
    if mol.atoms.is_empty() {
        return;
    }
    let cx = mol.atoms.iter().map(|a| a.x).sum::<f64>() / (mol.atoms.len() as f64);
    let cy = mol.atoms.iter().map(|a| a.y).sum::<f64>() / (mol.atoms.len() as f64);
    let yaw = (dx * 0.015).clamp(-1.1, 1.1);
    let pitch = (dy * 0.015).clamp(-1.1, 1.1);
    let (sy, cyw) = yaw.sin_cos();
    let (sp, cp) = pitch.sin_cos();
    let focal = 260.0;
    for a in &mut mol.atoms {
        let x0 = a.x - cx;
        let y0 = a.y - cy;
        let z0 = a.z_coord;
        // Rotate around Y (yaw)
        let x1 = x0 * cyw + z0 * sy;
        let z1 = -x0 * sy + z0 * cyw;
        // Rotate around X (pitch)
        let y2 = y0 * cp - z1 * sp;
        let z2 = y0 * sp + z1 * cp;
        let scale = focal / (focal - z2).max(40.0);
        a.x = cx + x1 * scale;
        a.y = cy + y2 * scale;
        a.z_coord = z2;
    }
}

/// ChiralCraft `Structure > Clean Up Reaction` (`P13` parity):
/// Aligns reactants to the left of the reaction arrow and products to the right along a common midline
/// with exact `15.10 pt` margins and `30.20 pt` minimum arrow length.
pub fn clean_up_reaction(scene: &mut DocumentScene) -> bool {
    if scene.molecules.len() < 2 {
        return false;
    }
    let (arrow_mid_x, arrow_y) = if let Some(arr) = scene.arrows.first_mut() {
        let mid_x = 0.5 * (arr.tail.0 + arr.head.0);
        let mid_y = 0.5 * (arr.tail.1 + arr.head.1);
        let half_len = ((arr.head.0 - arr.tail.0).hypot(arr.head.1 - arr.tail.1) * 0.5).max(15.10);
        arr.tail = (mid_x - half_len, mid_y);
        arr.head = (mid_x + half_len, mid_y);
        (mid_x, mid_y)
    } else {
        // Auto-place reaction arrow between first and second half of molecules
        let non_empty: Vec<&ChemMolecule> = scene
            .molecules
            .iter()
            .filter(|m| !m.atoms.is_empty())
            .collect();
        let denom = non_empty.len().max(1) as f64;
        let all_cx: f64 = non_empty
            .iter()
            .map(|m| m.atoms.iter().map(|a| a.x).sum::<f64>() / (m.atoms.len() as f64))
            .sum::<f64>()
            / denom;
        let all_cy: f64 = non_empty
            .iter()
            .map(|m| m.atoms.iter().map(|a| a.y).sum::<f64>() / (m.atoms.len() as f64))
            .sum::<f64>()
            / denom;
        (all_cx, all_cy)
    };

    let (arrow_left, arrow_right) = if let Some(arr) = scene.arrows.first() {
        (arr.tail.0.min(arr.head.0), arr.tail.0.max(arr.head.0))
    } else {
        (arrow_mid_x - 15.10, arrow_mid_x + 15.10)
    };

    let mut left_mols = Vec::new();
    let mut right_mols = Vec::new();
    for (mi, m) in scene.molecules.iter().enumerate() {
        if m.atoms.is_empty() {
            continue;
        }
        let cx = m.atoms.iter().map(|a| a.x).sum::<f64>() / (m.atoms.len() as f64);
        if cx <= arrow_mid_x {
            left_mols.push((mi, cx));
        } else {
            right_mols.push((mi, cx));
        }
    }
    left_mols.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    right_mols.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let margin = 15.10;
    let mut cursor_x = arrow_left - margin;
    for (mi, _) in left_mols {
        let m = &mut scene.molecules[mi];
        let min_x = m.atoms.iter().map(|a| a.x).fold(f64::INFINITY, f64::min);
        let max_x = m.atoms.iter().map(|a| a.x).fold(f64::NEG_INFINITY, f64::max);
        let min_y = m.atoms.iter().map(|a| a.y).fold(f64::INFINITY, f64::min);
        let max_y = m.atoms.iter().map(|a| a.y).fold(f64::NEG_INFINITY, f64::max);
        let dx = cursor_x - max_x;
        let dy = arrow_y - 0.5 * (min_y + max_y);
        for a in &mut m.atoms {
            a.x += dx;
            a.y += dy;
        }
        cursor_x -= (max_x - min_x) + margin * 1.4;
    }

    let mut cursor_x = arrow_right + margin;
    for (mi, _) in right_mols {
        let m = &mut scene.molecules[mi];
        let min_x = m.atoms.iter().map(|a| a.x).fold(f64::INFINITY, f64::min);
        let max_x = m.atoms.iter().map(|a| a.x).fold(f64::NEG_INFINITY, f64::max);
        let min_y = m.atoms.iter().map(|a| a.y).fold(f64::INFINITY, f64::min);
        let max_y = m.atoms.iter().map(|a| a.y).fold(f64::NEG_INFINITY, f64::max);
        let dx = cursor_x - min_x;
        let dy = arrow_y - 0.5 * (min_y + max_y);
        for a in &mut m.atoms {
            a.x += dx;
            a.y += dy;
        }
        cursor_x += (max_x - min_x) + margin * 1.4;
    }
    true
}

/// ChiralCraft `Structure > Map Reaction Atoms`:
/// Automatically pairs atoms on the reactant side (left of arrow) with matching-element atoms
/// on the product side (right of arrow), assigning sequential map IDs `1..=k`.
pub fn auto_map_reaction_atoms(scene: &DocumentScene) -> HashMap<(usize, usize), u32> {
    let mut map = HashMap::new();
    if scene.molecules.len() < 2 {
        return map;
    }
    let split_x = if let Some(arr) = scene.arrows.first() {
        0.5 * (arr.tail.0 + arr.head.0)
    } else {
        let xs: Vec<f64> = scene
            .molecules
            .iter()
            .filter(|m| !m.atoms.is_empty())
            .map(|m| m.atoms.iter().map(|a| a.x).sum::<f64>() / (m.atoms.len() as f64))
            .collect();
        xs.iter().sum::<f64>() / (xs.len().max(1) as f64)
    };

    let mut lhs_atoms: Vec<(usize, usize, i32)> = Vec::new();
    let mut rhs_atoms: Vec<(usize, usize, i32)> = Vec::new();
    for (mi, m) in scene.molecules.iter().enumerate() {
        if m.atoms.is_empty() {
            continue;
        }
        let cx = m.atoms.iter().map(|a| a.x).sum::<f64>() / (m.atoms.len() as f64);
        for (ai, a) in m.atoms.iter().enumerate() {
            if cx <= split_x {
                lhs_atoms.push((mi, ai, a.z));
            } else {
                rhs_atoms.push((mi, ai, a.z));
            }
        }
    }

    let mut rhs_used = vec![false; rhs_atoms.len()];
    let mut next_num = 1_u32;
    for &(lmi, lai, lz) in &lhs_atoms {
        if let Some(pos) = rhs_atoms
            .iter()
            .enumerate()
            .position(|(idx, &(_, _, rz))| !rhs_used[idx] && rz == lz)
        {
            rhs_used[pos] = true;
            let (rmi, rai, _) = rhs_atoms[pos];
            map.insert((lmi, lai), next_num);
            map.insert((rmi, rai), next_num);
            next_num += 1;
        }
    }
    map
}

/// Lists subdirectories and supported chemical files in `dir` for the Open/Save browser.
pub fn list_directory_entries(dir: &Path) -> Vec<(String, bool)> {
    let mut entries = Vec::new();
    if let Ok(rd) = fs::read_dir(dir) {
        for e in rd.flatten() {
            let name = e.file_name().to_string_lossy().to_string();
            if name.starts_with('.') {
                continue;
            }
            let is_dir = e.file_type().map(|t| t.is_dir()).unwrap_or(false);
            if is_dir || is_supported_chem_file(&name) {
                entries.push((name, is_dir));
            }
        }
    }
    entries.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    entries
}

fn is_supported_chem_file(name: &str) -> bool {
    let low = name.to_lowercase();
    if low == "poscar" || low == "contcar" || low.ends_with(".vasp") || low.ends_with(".poscar") {
        return true;
    }
    [".cdxml", ".cdx", ".cds", ".mol", ".sdf", ".rxn", ".ct", ".cml", ".tgf", ".skc", ".xyz", ".cif", ".pdb", ".out", ".log", ".inp", ".gjf", ".jdx", ".svg", ".png", ".pdf", ".eps", ".emf", ".sln", ".inchi"]
        .iter()
        .any(|ext| low.ends_with(ext))
}

/// Snaps a dragged bond endpoint `(ex, ey)` to ChiralCraft's `Fixed Lengths` and `Fixed Angles` (`15°` steps).
pub fn snap_bond_endpoint(
    sx: f64,
    sy: f64,
    ex: f64,
    ey: f64,
    fixed_lengths: bool,
    fixed_angles: bool,
    bond_length_pt: f64,
) -> (f64, f64) {
    let dx = ex - sx;
    let dy = ey - sy;
    let raw_len = dx.hypot(dy);
    if raw_len < 1e-4 {
        return (ex, ey);
    }
    let mut ang = dy.atan2(dx);
    if fixed_angles {
        let step = 15.0_f64.to_radians();
        ang = (ang / step).round() * step;
    }
    let len = if fixed_lengths { bond_length_pt } else { raw_len };
    (sx + len * ang.cos(), sy + len * ang.sin())
}

/// Standard ChiralCraft CPK element display color (`Colors > By Element`).
pub fn cpk_element_color(z: i32) -> egui::Color32 {
    use egui::Color32;
    match z {
        7 => Color32::from_rgb(0, 0, 220),       // N = blue
        8 => Color32::from_rgb(220, 0, 0),       // O = red
        9 => Color32::from_rgb(0, 160, 110),     // F = teal-green
        15 => Color32::from_rgb(210, 90, 0),     // P = orange
        16 => Color32::from_rgb(185, 150, 0),    // S = dark yellow
        17 => Color32::from_rgb(0, 155, 0),      // Cl = green
        35 => Color32::from_rgb(155, 50, 30),    // Br = brown
        53 => Color32::from_rgb(140, 0, 180),    // I = purple
        5 => Color32::from_rgb(190, 90, 120),    // B = rose
        14 => Color32::from_rgb(110, 110, 140),  // Si = slate
        _ => Color32::BLACK,
    }
}

/// Scales selected molecules around `(anchor_x, anchor_y)` by `(sx, sy)` (interactive 8-handle resize).
pub fn scale_selected_molecules(
    mols: &mut [ChemMolecule],
    selected: &[usize],
    anchor_x: f64,
    anchor_y: f64,
    sx: f64,
    sy: f64,
) {
    for &mi in selected {
        if let Some(m) = mols.get_mut(mi) {
            for a in &mut m.atoms {
                a.x = anchor_x + (a.x - anchor_x) * sx;
                a.y = anchor_y + (a.y - anchor_y) * sy;
            }
        }
    }
}

/// Rotates selected molecules around `(cx, cy)` by `angle_rad` (interactive top rotation handle).
pub fn rotate_selected_molecules(
    mols: &mut [ChemMolecule],
    selected: &[usize],
    cx: f64,
    cy: f64,
    angle_rad: f64,
) {
    let (s, c) = angle_rad.sin_cos();
    for &mi in selected {
        if let Some(m) = mols.get_mut(mi) {
            for a in &mut m.atoms {
                let dx = a.x - cx;
                let dy = a.y - cy;
                a.x = cx + dx * c - dy * s;
                a.y = cy + dx * s + dy * c;
            }
        }
    }
}

/// Parses a standard InChI string (`InChI=1S/<formula>/c<conn>...`) into a 2D `ChemMolecule`
/// for `Edit > Paste Special > InChI...`.
pub fn inchi_to_molecule(inchi: &str) -> Result<ChemMolecule, String> {
    let s = inchi.trim();
    if !s.starts_with("InChI=") {
        return Err("InChI string must start with 'InChI='".to_string());
    }
    let layers: Vec<&str> = s.split('/').collect();
    if layers.len() < 2 {
        return Err("Missing formula layer in InChI".to_string());
    }
    let formula = layers[1];
    // Parse heavy atoms from Hill formula layer (supporting multi-component `.`-separated fragments)
    const MAX_INCHI_ATOMS: usize = 10_000;
    let mut heavy_zs: Vec<i32> = Vec::new();
    let mut comp_offsets: Vec<usize> = Vec::new();
    for raw_comp in formula.split('.') {
        let mut comp_mult = 1usize;
        let mut rest = raw_comp;
        let lead_digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
        if !lead_digits.is_empty() {
            comp_mult = lead_digits.parse::<usize>().unwrap_or(1).clamp(1, 200);
            rest = &rest[lead_digits.len()..];
        }
        let mut comp_zs: Vec<i32> = Vec::new();
        let chars: Vec<char> = rest.chars().collect();
        let mut i = 0usize;
        while i < chars.len() {
            if chars[i].is_ascii_uppercase() {
                let mut sym = String::new();
                sym.push(chars[i]);
                i += 1;
                while i < chars.len() && chars[i].is_ascii_lowercase() {
                    sym.push(chars[i]);
                    i += 1;
                }
                let mut num_s = String::new();
                while i < chars.len() && chars[i].is_ascii_digit() {
                    num_s.push(chars[i]);
                    i += 1;
                }
                let cnt: usize = if num_s.is_empty() {
                    1
                } else {
                    num_s.parse::<usize>().map_err(|_| format!("Invalid atom count '{num_s}' in InChI formula"))?
                };
                if cnt > MAX_INCHI_ATOMS {
                    return Err(format!("InChI atom count exceeds maximum ({MAX_INCHI_ATOMS})"));
                }
                if sym != "H" {
                    let z = chiralcraft_chem::get_element_by_symbol(&sym)
                        .map(|e| e.z)
                        .ok_or_else(|| format!("Unknown element '{sym}' in InChI formula"))?;
                    if comp_zs.len().saturating_add(cnt) > MAX_INCHI_ATOMS {
                        return Err(format!("InChI atom count exceeds maximum ({MAX_INCHI_ATOMS})"));
                    }
                    for _ in 0..cnt {
                        comp_zs.push(z);
                    }
                }
            } else {
                i += 1;
            }
        }
        for _ in 0..comp_mult {
            if heavy_zs.len().saturating_add(comp_zs.len()) > MAX_INCHI_ATOMS {
                return Err(format!("InChI atom count exceeds maximum ({MAX_INCHI_ATOMS})"));
            }
            comp_offsets.push(heavy_zs.len());
            heavy_zs.extend_from_slice(&comp_zs);
        }
    }
    if heavy_zs.is_empty() {
        return Err("No heavy atoms found in InChI formula".to_string());
    }
    let mut mol = ChemMolecule::default();
    for (idx, &z) in heavy_zs.iter().enumerate() {
        let mut a = ChemAtom::new((idx + 1) as u32, z);
        a.x = (idx as f64) * 14.4;
        a.y = ((idx % 2) as f64) * 10.0;
        mol.atoms.push(a);
    }
    // Parse `/c` connectivity layer (e.g. `c1-2-3(4)5-1;1-2`)
    if let Some(c_layer) = layers.iter().find(|l| l.starts_with('c')) {
        let mut next_bid = (mol.atoms.len() + 1) as u32;
        let mut add_bond = |m: &mut ChemMolecule, u: usize, v: usize| {
            if u < m.atoms.len() && v < m.atoms.len() && u != v
                && !m.bonds.iter().any(|b| (b.begin_atom == u && b.end_atom == v) || (b.begin_atom == v && b.end_atom == u))
            {
                m.bonds.push(ChemBond {
                    id: next_bid,
                    begin_atom: u,
                    end_atom: v,
                    order: 1,
                    display: chiralcraft_chem::BondDisplay::Solid,
                    visible: true,
                    ignore_warnings: false,
                });
                next_bid += 1;
            }
        };
        // Expand `N*` multipliers (e.g. `/2C2H6O/c2*1-2-3`) so every repeated component gets its bonds
        let mut segs: Vec<&str> = Vec::new();
        for seg in c_layer[1..].split(';') {
            match seg.split_once('*') {
                Some((n, body)) if !n.is_empty() && n.chars().all(|c| c.is_ascii_digit()) => {
                    let reps = n.parse::<usize>().unwrap_or(1).clamp(1, comp_offsets.len().max(1));
                    segs.extend(std::iter::repeat(body).take(reps));
                }
                _ => segs.push(seg),
            }
        }
        for (comp_i, seg) in segs.into_iter().enumerate() {
            let offset = comp_offsets.get(comp_i).copied().unwrap_or(0);
            let bchars: Vec<char> = seg.chars().collect();
            let mut pos = 0usize;
            let mut stack: Vec<usize> = Vec::new();
            let mut prev: Option<usize> = None;
            while pos < bchars.len() {
                match bchars[pos] {
                    '0'..='9' => {
                        let mut ns = String::new();
                        while pos < bchars.len() && bchars[pos].is_ascii_digit() {
                            ns.push(bchars[pos]);
                            pos += 1;
                        }
                        if let Ok(n1) = ns.parse::<usize>() {
                            if n1 > 0 {
                                let idx = offset + (n1 - 1);
                                if let Some(p) = prev {
                                    add_bond(&mut mol, p, idx);
                                }
                                prev = Some(idx);
                            }
                        }
                    }
                    '-' => { pos += 1; }
                    '(' => {
                        if let Some(p) = prev {
                            stack.push(p);
                        }
                        pos += 1;
                    }
                    ',' => {
                        prev = stack.last().copied();
                        pos += 1;
                    }
                    ')' => {
                        prev = stack.pop();
                        pos += 1;
                    }
                    _ => { pos += 1; }
                }
            }
        }
    }
    Ok(mol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inchi_component_multiplier_bonds_every_copy() {
        // Two ethanol molecules: 6 heavy atoms, 2 bonds per copy
        let mol = inchi_to_molecule("InChI=1S/2C2H6O/c2*1-2-3/h2*3H,2H2,1H3").unwrap();
        assert_eq!(mol.atoms.len(), 6);
        assert_eq!(mol.bonds.len(), 4);
        assert!(mol.bonds.iter().any(|b| (b.begin_atom, b.end_atom) == (3, 4)));
        assert!(mol.bonds.iter().any(|b| (b.begin_atom, b.end_atom) == (4, 5)));
    }
}
