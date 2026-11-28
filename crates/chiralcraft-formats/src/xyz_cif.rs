//! Pure-Rust 3D `.xyz` & `.cif` (Crystallographic Information File) Importers
//!
//! - Parses Cartesian `.xyz` coordinates and crystallographic fractional `.cif` coordinates
//!   (full triclinic/monoclinic/orthorhombic/hexagonal/cubic fractional→Cartesian transform)
//! - Infers covalent bonds automatically via standard Pyykkö/Cordero covalent radii
//! - Preserves full `z_coord` (`has_3d_coords = true`) for 3D perspective tilt & Rear Opacity depth cueing

use chiralcraft_chem::{get_element_by_symbol, BondDisplay, ChemAtom, ChemMolecule};

/// Crystallographic unit cell in the same scaled coordinate frame as the imported `ChemMolecule`.
#[derive(Debug, Clone, PartialEq)]
pub struct CrystalUnitCell {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
    pub origin: (f64, f64, f64),
    pub va: (f64, f64, f64),
    pub vb: (f64, f64, f64),
    pub vc: (f64, f64, f64),
}

impl CrystalUnitCell {
    /// Returns the 8 corner vertices of the parallelepiped unit cell.
    pub fn corners(&self) -> [(f64, f64, f64); 8] {
        let (ox, oy, oz) = self.origin;
        let add = |u: (f64, f64, f64), v: (f64, f64, f64)| (u.0 + v.0, u.1 + v.1, u.2 + v.2);
        let c000 = (ox, oy, oz);
        let c100 = add(c000, self.va);
        let c010 = add(c000, self.vb);
        let c001 = add(c000, self.vc);
        let c110 = add(c100, self.vb);
        let c101 = add(c100, self.vc);
        let c011 = add(c010, self.vc);
        let c111 = add(c110, self.vc);
        [c000, c100, c010, c110, c001, c101, c011, c111]
    }

    /// Standard 12 edges indexing into `corners()`.
    pub fn edge_indices() -> [(usize, usize); 12] {
        [
            (0, 1), (0, 2), (0, 4),
            (1, 3), (1, 5),
            (2, 3), (2, 6),
            (4, 5), (4, 6),
            (3, 7), (5, 7), (6, 7),
        ]
    }
}

/// Approximate single-bond covalent radius in Ångströms for automatic 3D bond perception.
fn covalent_radius_angstrom(z: i32) -> f64 {
    match z {
        1 => 0.31,
        5 => 0.84,
        6 => 0.76,
        7 => 0.71,
        8 => 0.66,
        9 => 0.57,
        14 => 1.11,
        15 => 1.07,
        16 => 1.05,
        17 => 1.02,
        26 => 1.32,
        27 => 1.26,
        28 => 1.24,
        29 => 1.32,
        30 => 1.22,
        35 => 1.20,
        44 => 1.46,
        45 => 1.42,
        46 => 1.39,
        53 => 1.39,
        77 => 1.41,
        78 => 1.36,
        79 => 1.36,
        _ if z <= 10 => 0.75,
        _ if z <= 18 => 1.10,
        _ => 1.35,
    }
}

/// Infer bonds from 3D Cartesian coordinates (in Å) using covalent radii, then scale `(x, y, z)` to `bond_length_pt`.
/// Returns `(scale, (cx, cy, cz))` so crystallographic unit cell vectors can be transformed into the same coordinate frame.
fn perceive_3d_bonds_and_scale_with_transform(
    mol: &mut ChemMolecule,
    coords_ang: &[(f64, f64, f64)],
    target_blen: f64,
    perceive_bonds: bool,
) -> (f64, (f64, f64, f64)) {
    let n = mol.atoms.len();
    let mut next_bid = (mol.bonds.len() as u32) + (n as u32) + 1;
    let mut bond_lens_ang = Vec::new();

    // Deuterium/tritium pseudo-atoms (225/226) bond like hydrogen.
    let eff_z = |z: i32| if z == 225 || z == 226 { 1 } else { z };
    if perceive_bonds {
        for i in 0..n {
            let zi = eff_z(mol.atoms[i].z);
            let ri = covalent_radius_angstrom(zi);
            let (xi, yi, zci) = coords_ang[i];
            for j in (i + 1)..n {
                let zj = eff_z(mol.atoms[j].z);
                if zi == 1 && zj == 1 && n > 2 {
                    continue;
                }
                let rj = covalent_radius_angstrom(zj);
                let (xj, yj, zcj) = coords_ang[j];
                let d = ((xj - xi).powi(2) + (yj - yi).powi(2) + (zcj - zci).powi(2)).sqrt();
                let max_d = (ri + rj) * 1.22;
                if d > 0.42 && d <= max_d {
                    let ord = if zi == 6 && zj == 6 && d < 1.24 {
                        3
                    } else if zi == 6 && zj == 6 && d < 1.37 {
                        2
                    } else if ((zi == 6 && zj == 8) || (zi == 8 && zj == 6)) && d < 1.27 {
                        2
                    } else if ((zi == 6 && zj == 7) || (zi == 7 && zj == 6)) && d < 1.19 {
                        3
                    } else if ((zi == 6 && zj == 7) || (zi == 7 && zj == 6)) && d < 1.32 {
                        2
                    } else {
                        1
                    };
                    mol.add_bond(next_bid, i, j, ord, BondDisplay::Solid);
                    next_bid += 1;
                    if zi > 1 && zj > 1 {
                        bond_lens_ang.push(d);
                    }
                }
            }
        }
    } else {
        for b in &mol.bonds {
            if b.begin_atom < n && b.end_atom < n {
                let (xi, yi, zci) = coords_ang[b.begin_atom];
                let (xj, yj, zcj) = coords_ang[b.end_atom];
                let d = ((xj - xi).powi(2) + (yj - yi).powi(2) + (zcj - zci).powi(2)).sqrt();
                if d > 0.4 {
                    bond_lens_ang.push(d);
                }
            }
        }
    }

    let median_ang = if !bond_lens_ang.is_empty() {
        bond_lens_ang.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        bond_lens_ang[bond_lens_ang.len() / 2]
    } else {
        1.50
    };
    let scale = target_blen / median_ang.max(0.5);

    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for &(x, y, z) in coords_ang {
        cx += x;
        cy += y;
        cz += z;
    }
    let nf = (n.max(1)) as f64;
    cx /= nf;
    cy /= nf;
    cz /= nf;

    for (i, &(x, y, z)) in coords_ang.iter().enumerate() {
        mol.atoms[i].x = 220.0 + (x - cx) * scale;
        mol.atoms[i].y = 170.0 - (y - cy) * scale;
        mol.atoms[i].z_coord = (z - cz) * scale;
    }
    mol.has_3d_coords = true;
    (scale, (cx, cy, cz))
}

fn perceive_3d_bonds_and_scale(mol: &mut ChemMolecule, coords_ang: &[(f64, f64, f64)], target_blen: f64) {
    let _ = perceive_3d_bonds_and_scale_with_transform(mol, coords_ang, target_blen, true);
}

/// Parses an `.xyz` 3D molecular coordinate file into a `ChemMolecule` with perceived bonds and `has_3d_coords = true`.
pub fn xyz_to_molecule(xyz: &str) -> Result<ChemMolecule, String> {
    let lines: Vec<&str> = xyz.lines().collect();
    if lines.len() < 3 {
        return Err("XYZ file too short (expected atom count, comment, and coordinates)".to_string());
    }
    let expected_n: usize = lines[0]
        .trim()
        .parse()
        .map_err(|_| "Invalid atom count on line 1 of XYZ file".to_string())?;
    let mut mol = ChemMolecule::new();
    let mut coords = Vec::new();
    for (idx, line) in lines.iter().skip(2).enumerate() {
        if coords.len() >= expected_n {
            break;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 4 {
            continue;
        }
        // Atomic numbers ("6") are a common XYZ element form; otherwise strip label suffixes ("C12", "O2-").
        let z = match toks[0].parse::<i32>() {
            Ok(z) if (1..=118).contains(&z) => z,
            _ => {
                let sym = toks[0].trim_end_matches(|c: char| c.is_ascii_digit() || c == '+' || c == '-');
                match sym {
                    "D" => 225,
                    "T" => 226,
                    _ => get_element_by_symbol(sym).map(|e| e.z).unwrap_or(6),
                }
            }
        };
        let parse_coord = |t: &str, axis: &str| -> Result<f64, String> {
            t.parse::<f64>()
                .ok()
                .filter(|v| v.is_finite())
                .ok_or_else(|| format!("Bad {axis} coord at XYZ line {}", idx + 3))
        };
        let x = parse_coord(toks[1], "X")?;
        let y = parse_coord(toks[2], "Y")?;
        let zc = parse_coord(toks[3], "Z")?;
        mol.add_atom(ChemAtom::new((coords.len() + 1) as u32, z));
        coords.push((x, y, zc));
    }
    if mol.atoms.is_empty() {
        return Err("No valid atom coordinates found in XYZ file".to_string());
    }
    perceive_3d_bonds_and_scale(&mut mol, &coords, 14.40);
    Ok(mol)
}

fn parse_cif_float(s: &str) -> Option<f64> {
    // Strip crystallographic standard uncertainty in parentheses, e.g. "12.345(6)" -> 12.345
    let clean = s.split('(').next()?.trim();
    clean.parse::<f64>().ok().filter(|v| v.is_finite())
}

/// CIF data names are case-insensitive, and DDL2/mmCIF writes `_atom_site.fract_x` for `_atom_site_fract_x`.
fn normalize_cif_tag(tag: &str) -> String {
    tag.to_ascii_lowercase().replace('.', "_")
}

/// Splits one CIF line into value tokens: honours `'...'`/`"..."` quoting (a quote only closes when
/// followed by whitespace) and stops at a `#` comment outside quotes.
fn cif_line_tokens(line: &str) -> Vec<&str> {
    let b = line.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        while i < b.len() && b[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= b.len() || b[i] == b'#' {
            break;
        }
        let q = b[i];
        if q == b'\'' || q == b'"' {
            let start = i + 1;
            let mut j = start;
            while j < b.len() && !(b[j] == q && (j + 1 == b.len() || b[j + 1].is_ascii_whitespace())) {
                j += 1;
            }
            out.push(&line[start..j.min(b.len())]);
            i = j + 1;
        } else {
            let start = i;
            while i < b.len() && !b[i].is_ascii_whitespace() {
                i += 1;
            }
            out.push(&line[start..i]);
        }
    }
    out
}

fn extract_cif_element(raw: &str) -> i32 {
    let trimmed = raw.trim();
    if let Ok(z) = trimmed.parse::<i32>() {
        if (1..=118).contains(&z) {
            return z;
        }
    }
    let letters: String = trimmed
        .chars()
        .skip_while(|c| !c.is_ascii_alphabetic())
        .take_while(|c| c.is_ascii_alphabetic())
        .collect();
    if letters.eq_ignore_ascii_case("D") {
        return 225;
    }
    if letters.eq_ignore_ascii_case("T") {
        return 226;
    }
    if letters.len() >= 2 {
        let two = format!(
            "{}{}",
            letters.chars().next().unwrap().to_ascii_uppercase(),
            letters.chars().nth(1).unwrap().to_ascii_lowercase()
        );
        if let Some(el) = get_element_by_symbol(&two) {
            return el.z;
        }
    }
    if let Some(ch) = letters.chars().next() {
        let one = ch.to_ascii_uppercase().to_string();
        if let Some(el) = get_element_by_symbol(&one) {
            return el.z;
        }
    }
    6
}

/// Extracts atomic number from PDB columns 76..78 (explicit element) or 12..16 (atom name),
/// avoiding false matches of standard organic/protein atom names (`CA`, `CD`, `NE`, `HG`, `1HB`) as metals/noble gases!
fn extract_pdb_element(line: &str) -> i32 {
    if let Some(elem_col) = line.get(76..78).map(|s| s.trim()).filter(|s| !s.is_empty()) {
        return extract_cif_element(elem_col);
    }
    let name_col = line.get(12..16).unwrap_or(" C  ");
    // In PDB format, columns 12-13 (first two chars of 12..16) hold the right-justified element symbol
    // unless col 12 is a digit (e.g. "1HB ", "2HG1" -> hydrogen!)
    let c0 = name_col.as_bytes().first().copied().unwrap_or(b' ');
    if c0.is_ascii_digit() {
        let second = name_col.as_bytes().get(1).copied().unwrap_or(b'H') as char;
        if second.eq_ignore_ascii_case(&'H') {
            return 1;
        }
        if second.eq_ignore_ascii_case(&'D') {
            return 225;
        }
    }
    let trimmed = name_col.trim();
    // 4-character hydrogen names start in column 13 ("HG21", "HE21", "HO5'", "HD11"); a real
    // two-letter element there (Hg, He, Ho) is at most 2 characters long ("HG  ").
    if (c0 == b'H' || c0 == b'h') && trimmed.len() > 2 {
        return 1;
    }
    // Standard protein/nucleic organic atom names: CA/CB/CG/CD/CE/CZ, NA/NB/ND/NE/NZ, OA/OD/OE/OG/OH/OP, HA/HB/HD/HE/HG/HH/HZ, SD/SG
    if c0 == b' ' {
        if let Some(first_alpha) = trimmed.chars().find(|c| c.is_ascii_alphabetic()) {
            let one = first_alpha.to_ascii_uppercase().to_string();
            if matches!(one.as_str(), "C" | "N" | "O" | "H" | "S" | "P" | "F" | "I" | "B") {
                if let Some(el) = get_element_by_symbol(&one) {
                    return el.z;
                }
            }
        }
    }
    extract_cif_element(trimmed)
}

/// Parses a `.cif` (Crystallographic Information File) asymmetric unit into a 3D `ChemMolecule`
/// AND its scaled `CrystalUnitCell` using the exact crystallographic fractional-to-Cartesian matrix.
pub fn cif_to_molecule_with_cell(cif: &str) -> Result<(ChemMolecule, CrystalUnitCell), String> {
    let mut a = 10.0_f64;
    let mut b = 10.0_f64;
    let mut c = 10.0_f64;
    let mut alpha = 90.0_f64;
    let mut beta = 90.0_f64;
    let mut gamma = 90.0_f64;

    let lines: Vec<&str> = cif.lines().collect();
    for line in &lines {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut parts = t.split_whitespace();
        let Some(key) = parts.next() else { continue };
        let Some(val) = parts.next() else { continue };
        match normalize_cif_tag(key).as_str() {
            "_cell_length_a" => { if let Some(v) = parse_cif_float(val) { a = v; } }
            "_cell_length_b" => { if let Some(v) = parse_cif_float(val) { b = v; } }
            "_cell_length_c" => { if let Some(v) = parse_cif_float(val) { c = v; } }
            "_cell_angle_alpha" => { if let Some(v) = parse_cif_float(val) { alpha = v; } }
            "_cell_angle_beta" => { if let Some(v) = parse_cif_float(val) { beta = v; } }
            "_cell_angle_gamma" => { if let Some(v) = parse_cif_float(val) { gamma = v; } }
            _ => {}
        }
    }
    if [a, b, c].iter().any(|&l| l <= 0.0) || [alpha, beta, gamma].iter().any(|&g| g <= 0.0 || g >= 180.0) {
        return Err(format!("Invalid CIF unit cell: a={a} b={b} c={c} alpha={alpha} beta={beta} gamma={gamma}"));
    }

    // Build fractional -> Cartesian orthogonalization matrix (standard crystallographic convention)
    let ar = alpha.to_radians();
    let br = beta.to_radians();
    let gr = gamma.to_radians();
    let cos_a = ar.cos();
    let cos_b = br.cos();
    let cos_g = gr.cos();
    let sin_g = gr.sin().max(1e-6);
    let vol_factor = (1.0 - cos_a * cos_a - cos_b * cos_b - cos_g * cos_g + 2.0 * cos_a * cos_b * cos_g)
        .max(1e-9)
        .sqrt();
    let m11 = a;
    let m12 = b * cos_g;
    let m13 = c * cos_b;
    let m22 = b * sin_g;
    let m23 = c * (cos_a - cos_b * cos_g) / sin_g;
    let m33 = c * vol_factor / sin_g;

    let mut mol = ChemMolecule::new();
    let mut coords = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        if lines[i].trim().eq_ignore_ascii_case("loop_") {
            i += 1;
            let mut cols = Vec::new();
            while i < lines.len() && lines[i].trim().starts_with('_') {
                // A tag line holds just the data name (anything after it is a comment or stray text).
                cols.push(normalize_cif_tag(lines[i].split_whitespace().next().unwrap_or("")));
                i += 1;
            }
            let ix = cols.iter().position(|k| k == "_atom_site_fract_x");
            let iy = cols.iter().position(|k| k == "_atom_site_fract_y");
            let iz = cols.iter().position(|k| k == "_atom_site_fract_z");
            let isym = cols
                .iter()
                .position(|k| k == "_atom_site_type_symbol")
                .or_else(|| cols.iter().position(|k| k == "_atom_site_label"));
            if let (Some(ix), Some(iy), Some(iz), Some(isym)) = (ix, iy, iz, isym) {
                // Loop values form one token stream (a row may wrap lines; blank/comment lines are allowed);
                // the loop ends at the next data name or reserved word.
                let mut toks: Vec<&str> = Vec::new();
                while i < lines.len() {
                    let row = lines[i].trim();
                    let lower = row.get(..7).unwrap_or(row).to_ascii_lowercase();
                    if row.starts_with('_')
                        || ["loop_", "data_", "stop_", "save_", "global_"].iter().any(|w| lower.starts_with(w))
                    {
                        break;
                    }
                    toks.extend(cif_line_tokens(lines[i]));
                    i += 1;
                }
                for rec in toks.chunks_exact(cols.len()) {
                    if let (Some(fx), Some(fy), Some(fz)) = (
                        parse_cif_float(rec[ix]),
                        parse_cif_float(rec[iy]),
                        parse_cif_float(rec[iz]),
                    ) {
                        let z = extract_cif_element(rec[isym]);
                        let cx = m11 * fx + m12 * fy + m13 * fz;
                        let cy = m22 * fy + m23 * fz;
                        let cz = m33 * fz;
                        mol.add_atom(ChemAtom::new((coords.len() + 1) as u32, z));
                        coords.push((cx, cy, cz));
                    }
                }
                if !mol.atoms.is_empty() {
                    break;
                }
            }
            continue;
        }
        i += 1;
    }

    if mol.atoms.is_empty() {
        return Err("No _atom_site_fract_x/y/z records found in CIF file".to_string());
    }
    let (scale, (cx, cy, cz)) = perceive_3d_bonds_and_scale_with_transform(&mut mol, &coords, 14.40, true);
    let cell = CrystalUnitCell {
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
        origin: (220.0 - cx * scale, 170.0 + cy * scale, -cz * scale),
        va: (m11 * scale, 0.0, 0.0),
        vb: (m12 * scale, -m22 * scale, 0.0),
        vc: (m13 * scale, -m23 * scale, m33 * scale),
    };
    Ok((mol, cell))
}

/// Parses a `.cif` (Crystallographic Information File) asymmetric unit into a 3D `ChemMolecule`.
pub fn cif_to_molecule(cif: &str) -> Result<ChemMolecule, String> {
    cif_to_molecule_with_cell(cif).map(|(m, _)| m)
}

/// Parses a Protein Data Bank (`.pdb`) structure (`ATOM` / `HETATM` + optional `CONECT` records)
/// into a 3D `ChemMolecule` (`has_3d_coords = true`).
pub fn pdb_to_molecule(pdb: &str) -> Result<ChemMolecule, String> {
    let mut mol = ChemMolecule::new();
    let mut coords = Vec::new();
    let mut serial_to_idx: std::collections::HashMap<i32, usize> = std::collections::HashMap::new();
    let mut conect_pairs: Vec<(i32, i32)> = Vec::new();
    // Alternate locations (col 17): keep blank plus the first altLoc seen, not every conformer.
    let mut first_alt: Option<u8> = None;

    let mut in_first_model = true;
    for line in pdb.lines() {
        if line.starts_with("ENDMDL") && !mol.atoms.is_empty() {
            in_first_model = false;
            continue;
        }
        if !in_first_model && (line.starts_with("ATOM  ") || line.starts_with("HETATM")) {
            continue;
        }
        if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
            let serial: i32 = line
                .get(6..11)
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or((coords.len() + 1) as i32);
            let x = line.get(30..38).and_then(|s| s.trim().parse::<f64>().ok());
            let y = line.get(38..46).and_then(|s| s.trim().parse::<f64>().ok());
            let zc = line.get(46..54).and_then(|s| s.trim().parse::<f64>().ok());
            let (Some(x), Some(y), Some(zc)) = (x, y, zc) else {
                continue;
            };
            if !x.is_finite() || !y.is_finite() || !zc.is_finite() {
                continue;
            }
            let alt = line.as_bytes().get(16).copied().unwrap_or(b' ');
            if alt != b' ' {
                if *first_alt.get_or_insert(alt) != alt {
                    continue;
                }
            }
            let z = extract_pdb_element(line);
            let idx = mol.atoms.len();
            serial_to_idx.insert(serial, idx);
            mol.add_atom(ChemAtom::new((idx + 1) as u32, z));
            coords.push((x, y, zc));
        } else if line.starts_with("CONECT") {
            // Fixed-width 5-column PDB CONECT fields (handles >= 10,000 atoms where fields touch).
            // Only cols 7-31 (atom + 4 bonded atoms) are covalent bonds; legacy cols 32-61 are H-bonds/salt bridges.
            let mut toks: Vec<i32> = Vec::new();
            let mut pos = 6usize;
            let lim = line.len().min(31);
            while pos < lim {
                let end = (pos + 5).min(lim);
                if let Some(chunk) = line.get(pos..end) {
                    if let Ok(v) = chunk.trim().parse::<i32>() {
                        toks.push(v);
                    }
                }
                pos = end;
            }
            if toks.len() >= 2 {
                let a = toks[0];
                // Keep both directions: many writers list a bond only on one atom's record.
                for &b in &toks[1..] {
                    conect_pairs.push((a, b));
                }
            }
        }
    }
    if mol.atoms.is_empty() {
        return Err("No ATOM/HETATM coordinate records found in PDB file".to_string());
    }
    // Repeated entries in one atom's CONECT records encode bond order (count per direction, max of both).
    let mut directed: std::collections::HashMap<(usize, usize), i32> = std::collections::HashMap::new();
    for (sa, sb) in conect_pairs {
        if let (Some(&ia), Some(&ib)) = (serial_to_idx.get(&sa), serial_to_idx.get(&sb)) {
            if ia != ib {
                *directed.entry((ia, ib)).or_insert(0) += 1;
            }
        }
    }
    let mut conect_bonds: Vec<(usize, usize, i32)> = Vec::new();
    let mut seen: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    let mut keys: Vec<(usize, usize)> = directed.keys().copied().collect();
    keys.sort_unstable();
    for (ia, ib) in keys {
        let key = (ia.min(ib), ia.max(ib));
        if seen.insert(key) {
            let n1 = directed.get(&(ia, ib)).copied().unwrap_or(0);
            let n2 = directed.get(&(ib, ia)).copied().unwrap_or(0);
            conect_bonds.push((key.0, key.1, n1.max(n2).clamp(1, 3)));
        }
    }
    // CONECT usually covers only HETATM/non-standard links; atoms it never mentions (e.g. the protein)
    // still need distance-based bonds, but perceived bonds touching CONECT atoms are discarded.
    let mut covered = vec![false; mol.atoms.len()];
    for &(a, b, _) in &conect_bonds {
        covered[a] = true;
        covered[b] = true;
    }
    let _ = perceive_3d_bonds_and_scale_with_transform(&mut mol, &coords, 14.40, true);
    if !conect_bonds.is_empty() {
        mol.bonds.retain(|b| !covered[b.begin_atom] && !covered[b.end_atom]);
        let mut next_bid = mol.bonds.iter().map(|b| b.id).max().unwrap_or(0).max(mol.atoms.len() as u32) + 1;
        for (a, b, ord) in conect_bonds {
            mol.add_bond(next_bid, a, b, ord, BondDisplay::Solid);
            next_bid += 1;
        }
    }
    Ok(mol)
}

/// Parses a VASP `POSCAR` / `CONTCAR` crystal file into a 3D `ChemMolecule` + `CrystalUnitCell`.
pub fn poscar_to_molecule_with_cell(poscar: &str) -> Result<(ChemMolecule, CrystalUnitCell), String> {
    let lines: Vec<&str> = poscar.lines().map(|l| l.trim()).collect();
    if lines.len() < 8 {
        return Err("POSCAR/CONTCAR file too short".to_string());
    }
    let raw_scale: f64 = lines[1]
        .split_whitespace()
        .next()
        .and_then(|s| s.parse().ok())
        .filter(|v: &f64| v.is_finite() && *v != 0.0)
        .ok_or_else(|| "Invalid scale factor on line 2 of POSCAR".to_string())?;
    let parse_raw_vec3 = |line: &str| -> Result<(f64, f64, f64), String> {
        let vals: Vec<f64> = line
            .split_whitespace()
            .take(3)
            .filter_map(|t| t.parse::<f64>().ok())
            .filter(|v| v.is_finite())
            .collect();
        if vals.len() < 3 {
            return Err(format!("Expected 3 lattice vector components in POSCAR line: {line}"));
        }
        Ok((vals[0], vals[1], vals[2]))
    };
    let va_raw = parse_raw_vec3(lines[2])?;
    let vb_raw = parse_raw_vec3(lines[3])?;
    let vc_raw = parse_raw_vec3(lines[4])?;
    // Per VASP spec: negative scale_factor is target unit-cell volume (Å³), not linear scale!
    let scale_factor = if raw_scale < 0.0 {
        let det = (va_raw.0 * (vb_raw.1 * vc_raw.2 - vb_raw.2 * vc_raw.1)
            - va_raw.1 * (vb_raw.0 * vc_raw.2 - vb_raw.2 * vc_raw.0)
            + va_raw.2 * (vb_raw.0 * vc_raw.1 - vb_raw.1 * vc_raw.0))
            .abs();
        if det < 1e-12 {
            return Err("Degenerate POSCAR lattice vectors for volume scaling".to_string());
        }
        ((-raw_scale) / det).cbrt()
    } else {
        raw_scale
    };
    let va_ang = (va_raw.0 * scale_factor, va_raw.1 * scale_factor, va_raw.2 * scale_factor);
    let vb_ang = (vb_raw.0 * scale_factor, vb_raw.1 * scale_factor, vb_raw.2 * scale_factor);
    let vc_ang = (vc_raw.0 * scale_factor, vc_raw.1 * scale_factor, vc_raw.2 * scale_factor);

    // Line 6 may be element symbols (VASP 5+) or counts (VASP 4)
    let tok5: Vec<&str> = lines[5].split_whitespace().collect();
    let (species, counts_line_idx): (Vec<i32>, usize) = if tok5.iter().all(|t| t.parse::<usize>().is_ok()) {
        (vec![6; tok5.len()], 5)
    } else {
        let zs = tok5.iter().map(|s| extract_cif_element(s)).collect();
        (zs, 6)
    };
    let counts: Vec<usize> = lines
        .get(counts_line_idx)
        .ok_or_else(|| "Missing atom counts line in POSCAR".to_string())?
        .split_whitespace()
        .filter_map(|t| t.parse::<usize>().ok())
        .collect();
    if counts.is_empty() {
        return Err("No atom counts found in POSCAR".to_string());
    }
    let mut mode_idx = counts_line_idx + 1;
    if lines
        .get(mode_idx)
        .map(|l| l.to_ascii_lowercase().starts_with('s'))
        .unwrap_or(false)
    {
        mode_idx += 1; // Skip "Selective dynamics"
    }
    let mode_line = lines
        .get(mode_idx)
        .ok_or_else(|| "Missing Direct/Cartesian line in POSCAR".to_string())?;
    let is_direct = !mode_line.to_ascii_lowercase().starts_with('c')
        && !mode_line.to_ascii_lowercase().starts_with('k');

    let avail_lines = lines.len().saturating_sub(mode_idx + 1);
    let total_atoms = counts
        .iter()
        .try_fold(0usize, |acc, &c| acc.checked_add(c))
        .ok_or_else(|| "POSCAR atom count overflow".to_string())?;
    if total_atoms > avail_lines || total_atoms > 200_000 {
        return Err(format!(
            "POSCAR atom count ({total_atoms}) exceeds available coordinate lines ({avail_lines})"
        ));
    }
    let mut atom_zs = Vec::with_capacity(total_atoms);
    for (si, &cnt) in counts.iter().enumerate() {
        let z = species.get(si).copied().unwrap_or(6);
        for _ in 0..cnt {
            atom_zs.push(z);
        }
    }
    let mut mol = ChemMolecule::new();
    let mut coords = Vec::new();
    for (k, &z) in atom_zs.iter().enumerate() {
        let Some(line) = lines.get(mode_idx + 1 + k) else {
            break;
        };
        let vals: Vec<f64> = line
            .split_whitespace()
            .take(3)
            .filter_map(|t| t.parse::<f64>().ok())
            .collect();
        if vals.len() < 3 {
            continue;
        }
        let (x, y, zc) = if is_direct {
            (
                vals[0] * va_ang.0 + vals[1] * vb_ang.0 + vals[2] * vc_ang.0,
                vals[0] * va_ang.1 + vals[1] * vb_ang.1 + vals[2] * vc_ang.1,
                vals[0] * va_ang.2 + vals[1] * vb_ang.2 + vals[2] * vc_ang.2,
            )
        } else {
            (vals[0] * scale_factor, vals[1] * scale_factor, vals[2] * scale_factor)
        };
        mol.add_atom(ChemAtom::new((coords.len() + 1) as u32, z));
        coords.push((x, y, zc));
    }
    if mol.atoms.is_empty() {
        return Err("No valid atom positions found in POSCAR".to_string());
    }
    let (scale, (cx, cy, cz)) = perceive_3d_bonds_and_scale_with_transform(&mut mol, &coords, 14.40, true);
    let norm = |v: (f64, f64, f64)| (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt();
    let dot = |u: (f64, f64, f64), v: (f64, f64, f64)| u.0 * v.0 + u.1 * v.1 + u.2 * v.2;
    let a = norm(va_ang).max(1e-6);
    let b = norm(vb_ang).max(1e-6);
    let c = norm(vc_ang).max(1e-6);
    let alpha = (dot(vb_ang, vc_ang) / (b * c)).clamp(-1.0, 1.0).acos().to_degrees();
    let beta = (dot(va_ang, vc_ang) / (a * c)).clamp(-1.0, 1.0).acos().to_degrees();
    let gamma = (dot(va_ang, vb_ang) / (a * b)).clamp(-1.0, 1.0).acos().to_degrees();
    let cell = CrystalUnitCell {
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
        origin: (220.0 - cx * scale, 170.0 + cy * scale, -cz * scale),
        va: (va_ang.0 * scale, -va_ang.1 * scale, va_ang.2 * scale),
        vb: (vb_ang.0 * scale, -vb_ang.1 * scale, vb_ang.2 * scale),
        vc: (vc_ang.0 * scale, -vc_ang.1 * scale, vc_ang.2 * scale),
    };
    Ok((mol, cell))
}

/// Parses a VASP `POSCAR` / `CONTCAR` file into a 3D `ChemMolecule`.
pub fn poscar_to_molecule(poscar: &str) -> Result<ChemMolecule, String> {
    poscar_to_molecule_with_cell(poscar).map(|(m, _)| m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xyz_and_cif_3d_import_and_bonding() {
        let water_xyz = "\
3
Water molecule
O   0.0000   0.0000   0.0000
H   0.7570   0.5860   0.0000
H  -0.7570   0.5860   0.0000
";
        let w = xyz_to_molecule(water_xyz).expect("xyz water");
        assert_eq!(w.atoms.len(), 3);
        assert_eq!(w.bonds.len(), 2);
        assert!(w.has_3d_coords);

        let sample_cif = "\
data_ethanol
_cell_length_a   10.000(2)
_cell_length_b   10.000(2)
_cell_length_c   10.000(2)
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.000 0.000 0.000
C2 C 0.154 0.000 0.000
O1 O 0.220 0.125 0.050
";
        let (eth, cell) = cif_to_molecule_with_cell(sample_cif).expect("cif ethanol");
        assert_eq!(eth.atoms.len(), 3);
        assert_eq!(eth.bonds.len(), 2);
        assert!(eth.has_3d_coords);
        assert!(eth.atoms[2].z_coord.abs() > 0.1);
        assert!((cell.a - 10.0).abs() < 1e-6);
        assert_eq!(cell.corners().len(), 8);

        let sample_pdb = "\
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C2  LIG A   1       1.520   0.000   0.000  1.00  0.00           C
ATOM      3  O1  LIG A   1       2.150   1.210   0.450  1.00  0.00           O
CONECT    1    2
CONECT    2    3
END
";
        let pdb_mol = pdb_to_molecule(sample_pdb).expect("pdb ethanol");
        assert_eq!(pdb_mol.atoms.len(), 3);
        assert_eq!(pdb_mol.bonds.len(), 2);
        assert!(pdb_mol.has_3d_coords);

        let sample_poscar = "\
Si2 diamond demo
1.0
5.43 0.00 0.00
0.00 5.43 0.00
0.00 0.00 5.43
Si
2
Direct
0.00 0.00 0.00
0.25 0.25 0.25
";
        let (si_mol, si_cell) = poscar_to_molecule_with_cell(sample_poscar).expect("poscar Si");
        assert_eq!(si_mol.atoms.len(), 2);
        assert_eq!(si_mol.atoms[0].z, 14);
        assert_eq!(si_mol.bonds.len(), 1);
        assert!((si_cell.a - 5.43).abs() < 1e-6);
    }

    #[test]
    fn test_xyz_cif_pdb_audit_regressions() {
        // XYZ: numeric atomic numbers, deuterium bonds like H, NaN rejected.
        let xyz = "3\nwater-d2\n8 0.0 0.0 0.0\nD 0.96 0.0 0.0\nD -0.24 0.93 0.0\n";
        let m = xyz_to_molecule(xyz).unwrap();
        assert_eq!((m.atoms[0].z, m.atoms[1].z), (8, 225));
        assert_eq!(m.bonds.len(), 2, "O-D bonds only (no D-D)");
        assert!(xyz_to_molecule("1\n\nC nan 0 0\n").is_err());
        // CIF: mmCIF dotted/uppercase tags, quoted labels, a row wrapped over two lines, blank line inside loop.
        let cif = "data_x\n_CELL.LENGTH_A 5.0\n_cell.length_b 5.0\n_cell.length_c 5.0\nloop_\n_atom_site.label\n_atom_site.type_symbol\n_atom_site.fract_x\n_atom_site.fract_y\n_atom_site.fract_z\n'C 1' C 0.0 0.0 0.0\n\nO1 O 0.25\n 0.0 0.0\n";
        let (cm, cell) = cif_to_molecule_with_cell(cif).unwrap();
        assert_eq!(cm.atoms.len(), 2);
        assert_eq!(cm.atoms[1].z, 8);
        assert!((cell.a - 5.0).abs() < 1e-9);
        assert!(cif_to_molecule("_cell_length_a 0\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nC1 0 0 0\n").is_err());
        // PDB: 4-char H names are hydrogen, altLoc B dropped, one-directional CONECT kept,
        // protein atoms without CONECT still perceived, H-bond columns ignored.
        let pdb = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00
ATOM      2  CA  ALA A   1       1.460   0.000   0.000  1.00  0.00
ATOM      3 HG21 ALA A   1       1.800   1.000   0.000  1.00  0.00
ATOM      4  CB AALA A   1       2.000  -1.400   0.000  0.50  0.00
ATOM      5  CB BALA A   1       2.100  -1.300   0.200  0.50  0.00
HETATM    6  C1  LIG B   1      10.000   0.000   0.000  1.00  0.00
HETATM    7  O1  LIG B   1      11.200   0.000   0.000  1.00  0.00
CONECT    6    7    7                                  1
END
";
        let p = pdb_to_molecule(pdb).unwrap();
        assert_eq!(p.atoms.len(), 6, "altLoc B skipped");
        assert_eq!(p.atoms[2].z, 1, "HG21 is hydrogen, not mercury");
        let has = |a: usize, b: usize| p.bonds.iter().find(|x| (x.begin_atom, x.end_atom) == (a.min(b), a.max(b)) || (x.begin_atom, x.end_atom) == (a.max(b), a.min(b)));
        assert!(has(0, 1).is_some(), "protein N-CA perceived even though CONECT exists");
        assert_eq!(has(4, 5).map(|b| b.order), Some(2), "repeated CONECT entry = double bond");
        assert!(has(4, 0).is_none(), "H-bond column (atom 1) ignored");
    }
}
