//! ChiralCraft Biopolymer (HELM / FASTA / Official Monomer JSON) & 3MF 3D Export!
//!
//! - Parses shipped ChiralCraft `natural-monomers.json` (and `pistoia-alliance-monomers.json`)
//! - Converts FASTA Peptide / DNA / RNA and HELM v1/v2 notation to linked 2D `ChemMolecule`
//! - Serializes `ChemMolecule` sequences back to HELM (`PEPTIDE1{A.C.D}$$$$V2.0`) and FASTA
//! - Exports valid 3MF 3D model XML (`Ball and Stick` and `Stick` modes)

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use chiralcraft_chem::{get_element, BondDisplay, ChemAtom, ChemBond, ChemMolecule};

use crate::smiles_to_molecule;

#[derive(Debug, Clone)]
pub struct MonomerEntry {
    pub symbol: String,
    pub name: String,
    pub polymer_type: String, // "PEPTIDE" | "RNA" | "CHEM"
    pub natural_analog: String,
    pub smiles: String,
}

#[derive(Debug, Clone, Default)]
pub struct MonomerLibrary {
    pub entries: Vec<MonomerEntry>,
    pub by_symbol: HashMap<(String, String), usize>, // (polymer_type, symbol) -> idx
}

impl MonomerLibrary {
    /// Loads official ChiralCraft monomer JSON (`natural-monomers.json`) using a fast zero-dep parser.
    pub fn load_from_json_path(path: &Path) -> Result<Self, String> {
        let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
        Ok(Self::parse_json(&text))
    }

    pub fn parse_json(text: &str) -> Self {
        let mut lib = Self::default();
        // Split top-level objects cleanly by `"id":` field boundaries inside the JSON array
        let mut pos = 0usize;
        let bytes = text.as_bytes();
        while let Some(rel) = text[pos..].find("\"polymerType\"") {
            let anchor = pos + rel;
            // Find enclosing `{` before anchor and matching `}` after anchor
            let obj_start = text[..anchor].rfind("\n  {").or_else(|| text[..anchor].rfind('{')).unwrap_or(0);
            let mut depth = 0i32;
            let mut in_str = false;
            let mut esc = false;
            let mut obj_end = text.len();
            for (k, &b) in bytes[obj_start..].iter().enumerate() {
                if in_str {
                    if esc {
                        esc = false;
                    } else if b == b'\\' {
                        esc = true;
                    } else if b == b'"' {
                        in_str = false;
                    }
                } else {
                    match b {
                        b'"' => in_str = true,
                        b'{' => depth += 1,
                        b'}' => {
                            depth -= 1;
                            if depth == 0 {
                                obj_end = obj_start + k + 1;
                                break;
                            }
                        }
                        _ => {}
                    }
                }
            }
            let chunk = &text[obj_start..obj_end];
            let symbol = extract_json_string(chunk, "symbol").unwrap_or_default();
            let name = extract_json_string(chunk, "name").unwrap_or_else(|| symbol.clone());
            let polymer_type = extract_json_string(chunk, "polymerType").unwrap_or_else(|| "PEPTIDE".to_string());
            let natural_analog = extract_json_string(chunk, "naturalAnalog").unwrap_or_else(|| symbol.clone());
            let smiles = extract_json_string(chunk, "smiles").unwrap_or_default();
            if !symbol.is_empty() {
                let key = (polymer_type.clone(), symbol.clone());
                if !lib.by_symbol.contains_key(&key) {
                    let idx = lib.entries.len();
                    lib.entries.push(MonomerEntry {
                        symbol,
                        name,
                        polymer_type,
                        natural_analog,
                        smiles,
                    });
                    lib.by_symbol.insert(key, idx);
                }
            }
            pos = obj_end.max(anchor + 13);
            if pos >= text.len() {
                break;
            }
        }
        if lib.entries.is_empty() {
            lib.populate_builtin_defaults();
        }
        lib
    }

    pub fn populate_builtin_defaults(&mut self) {
        for &(sym, name, smi) in STANDARD_AMINO_ACIDS {
            let idx = self.entries.len();
            self.entries.push(MonomerEntry {
                symbol: sym.to_string(),
                name: name.to_string(),
                polymer_type: "PEPTIDE".to_string(),
                natural_analog: sym.to_string(),
                smiles: smi.to_string(),
            });
            self.by_symbol.insert(("PEPTIDE".to_string(), sym.to_string()), idx);
        }
        for &(sym, name, smi) in STANDARD_NUCLEOBASES {
            let idx = self.entries.len();
            self.entries.push(MonomerEntry {
                symbol: sym.to_string(),
                name: name.to_string(),
                polymer_type: "RNA".to_string(),
                natural_analog: sym.to_string(),
                smiles: smi.to_string(),
            });
            self.by_symbol.insert(("RNA".to_string(), sym.to_string()), idx);
        }
    }
}

fn extract_json_string(obj: &str, key: &str) -> Option<String> {
    let pat = format!("\"{key}\"");
    let kpos = obj.find(&pat)?;
    let rest = &obj[kpos + pat.len()..];
    let cpos = rest.find(':')?;
    let val = rest[cpos + 1..].trim_start();
    if !val.starts_with('"') {
        return None;
    }
    let mut out = String::new();
    let mut esc = false;
    for ch in val[1..].chars() {
        if esc {
            match ch {
                'n' => out.push('\n'),
                't' => out.push('\t'),
                '"' => out.push('"'),
                '\\' => out.push('\\'),
                other => out.push(other),
            }
            esc = false;
        } else if ch == '\\' {
            esc = true;
        } else if ch == '"' {
            return Some(out);
        } else {
            out.push(ch);
        }
    }
    None
}

const STANDARD_AMINO_ACIDS: &[(&str, &str, &str)] = &[
    ("A", "Alanine", "CC(N)C(=O)O"),
    ("R", "Arginine", "NC(=N)NCCCC(N)C(=O)O"),
    ("N", "Asparagine", "NC(=O)CC(N)C(=O)O"),
    ("D", "Aspartic Acid", "O=C(O)CC(N)C(=O)O"),
    ("C", "Cysteine", "SCC(N)C(=O)O"),
    ("E", "Glutamic Acid", "O=C(O)CCC(N)C(=O)O"),
    ("Q", "Glutamine", "NC(=O)CCC(N)C(=O)O"),
    ("G", "Glycine", "NCC(=O)O"),
    ("H", "Histidine", "c1c[nH]cn1CC(N)C(=O)O"),
    ("I", "Isoleucine", "CCC(C)C(N)C(=O)O"),
    ("L", "Leucine", "CC(C)CC(N)C(=O)O"),
    ("K", "Lysine", "NCCCCC(N)C(=O)O"),
    ("M", "Methionine", "CSCCC(N)C(=O)O"),
    ("F", "Phenylalanine", "c1ccccc1CC(N)C(=O)O"),
    ("P", "Proline", "O=C(O)C1CCCN1"),
    ("S", "Serine", "OCC(N)C(=O)O"),
    ("T", "Threonine", "CC(O)C(N)C(=O)O"),
    ("W", "Tryptophan", "c1ccc2[nH]cc(CC(N)C(=O)O)c2c1"),
    ("Y", "Tyrosine", "Oc1ccc(CC(N)C(=O)O)cc1"),
    ("V", "Valine", "CC(C)C(N)C(=O)O"),
];

const STANDARD_NUCLEOBASES: &[(&str, &str, &str)] = &[
    ("A", "Adenine", "Nc1ncnc2[nH]cnc12"),
    ("C", "Cytosine", "Nc1cc[nH]c(=O)n1"),
    ("G", "Guanine", "Nc1nc2[nH]cnc2c(=O)[nH]1"),
    ("T", "Thymine", "Cc1c[nH]c(=O)[nH]c1=O"),
    ("U", "Uracil", "O=c1cc[nH]c(=O)[nH]1"),
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PolymerKind {
    Peptide,
    Dna,
    Rna,
}

/// Converts a FASTA / 1-letter sequence (`>header\nACDEF` or `ACDEF`) into a linked biopolymer
/// `ChemMolecule` with monomer residue labels (`A`, `C`, `D`, ...) laid out cleanly left-to-right!
pub fn sequence_to_biopolymer_molecule(
    seq_input: &str,
    kind: PolymerKind,
    bond_length_pt: f64,
) -> Result<ChemMolecule, String> {
    let letters: Vec<String> = seq_input
        .lines()
        .filter(|l| !l.trim().starts_with('>'))
        .flat_map(|l| l.chars())
        .filter(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase().to_string())
        .collect();
    if letters.is_empty() {
        return Err("No valid monomer letters found in sequence".to_string());
    }
    build_monomer_chain(&letters, kind, bond_length_pt)
}

/// Parses a HELM v1/v2 string (e.g. `PEPTIDE1{A.C.D.E.F}$$$$V2.0` or `RNA1{R(A)P.R(U)P.R(G)}$$$$V2.0`)
/// into a linked biopolymer `ChemMolecule` (including cyclic/disulfide cross-links if specified).
pub fn helm_to_molecule(helm: &str, bond_length_pt: f64) -> Result<ChemMolecule, String> {
    let trimmed = helm.trim();
    let sections: Vec<&str> = trimmed.split('$').collect();
    let poly_part = sections.first().copied().unwrap_or(trimmed);
    let first_poly = poly_part.split('|').next().unwrap_or(poly_part);
    let brace_start = first_poly
        .find('{')
        .ok_or_else(|| "HELM string must contain `{monomers}`".to_string())?;
    let brace_end = first_poly
        .rfind('}')
        .ok_or_else(|| "HELM string missing closing `}`".to_string())?;
    let header = first_poly[..brace_start].to_ascii_uppercase();
    let kind = if header.starts_with("RNA") {
        PolymerKind::Rna
    } else if header.starts_with("DNA") {
        PolymerKind::Dna
    } else {
        PolymerKind::Peptide
    };
    let body = &first_poly[brace_start + 1..brace_end];
    let monomers: Vec<String> = body
        .split('.')
        .map(|m| m.trim().trim_matches(['[', ']']).to_string())
        .filter(|m| !m.is_empty())
        .collect();
    if monomers.is_empty() {
        return Err("No monomers found inside HELM `{...}`".to_string());
    }
    let mut mol = build_monomer_chain(&monomers, kind, bond_length_pt)?;
    // Parse optional HELM connection section (e.g. `PEPTIDE1,PEPTIDE1,1:R3-5:R3`)
    if let Some(conn_sec) = sections.get(1) {
        let mut next_bid = (mol.atoms.len() + mol.bonds.len() + 10) as u32;
        for conn in conn_sec.split('|').filter(|c| !c.trim().is_empty()) {
            let cols: Vec<&str> = conn.split(',').collect();
            if cols.len() >= 3 {
                if let Some((lhs, rhs)) = cols[2].split_once('-') {
                    let i1 = lhs.split(':').next().and_then(|s| s.parse::<usize>().ok());
                    let i2 = rhs.split(':').next().and_then(|s| s.parse::<usize>().ok());
                    if let (Some(a), Some(b)) = (i1, i2) {
                        if a >= 1 && b >= 1 && a <= mol.atoms.len() && b <= mol.atoms.len() && a != b {
                            mol.bonds.push(ChemBond {
                                id: next_bid,
                                begin_atom: a - 1,
                                end_atom: b - 1,
                                order: 1,
                                display: BondDisplay::Solid,
                                visible: true,
                                ignore_warnings: false,
                            });
                            next_bid += 1;
                        }
                    }
                }
            }
        }
    }
    Ok(mol)
}

fn build_monomer_chain(
    monomers: &[String],
    kind: PolymerKind,
    bond_length_pt: f64,
) -> Result<ChemMolecule, String> {
    let step = (bond_length_pt * 1.65).max(20.0);
    let mut mol = ChemMolecule::default();
    let per_row = 12usize;
    for (i, sym) in monomers.iter().enumerate() {
        let row = i / per_row;
        let col = i % per_row;
        let eff_col = if row % 2 == 0 { col } else { per_row - 1 - col };
        let mut a = ChemAtom::new((i + 1) as u32, 6);
        a.x = 120.0 + (eff_col as f64) * step;
        a.y = 140.0 + (row as f64) * (step * 1.4);
        a.label = Some(match kind {
            PolymerKind::Peptide => {
                if sym.len() == 1 {
                    one_to_three_aa(sym).unwrap_or(sym).to_string()
                } else {
                    sym.clone()
                }
            }
            PolymerKind::Dna => {
                if sym.len() == 1 { format!("d{sym}") } else { sym.clone() }
            }
            PolymerKind::Rna => sym.clone(),
        });
        mol.atoms.push(a);
        if i > 0 {
            mol.bonds.push(ChemBond {
                id: (monomers.len() + i + 1) as u32,
                begin_atom: i - 1,
                end_atom: i,
                order: 1,
                display: BondDisplay::Solid,
                visible: true,
                ignore_warnings: false,
            });
        }
    }
    Ok(mol)
}

/// Expands a single monomer symbol into a full all-atom `ChemMolecule` from `MonomerLibrary`.
pub fn monomer_to_all_atom_molecule(lib: &MonomerLibrary, polymer_type: &str, symbol: &str) -> Option<ChemMolecule> {
    let idx = *lib.by_symbol.get(&(polymer_type.to_string(), symbol.to_string()))?;
    let entry = &lib.entries[idx];
    // Strip atom-map annotations (`[H:1]`, `[OH:2]`) for clean SMILES parsing
    let cleaned = entry
        .smiles
        .replace("[H:1]", "")
        .replace("[OH:2]", "O")
        .replace("[H:3]", "")
        .replace("[*:1]", "")
        .replace("[*:2]", "");
    smiles_to_molecule(&cleaned).ok()
}

/// Serializes a biopolymer or peptide molecule into **HELM v2.0** notation (`PEPTIDE1{A.C.D}$$$$V2.0`).
pub fn molecule_to_helm(mol: &ChemMolecule) -> String {
    let mut tokens: Vec<String> = Vec::new();
    let mut is_nucleic = false;
    for a in &mol.atoms {
        if let Some(lbl) = &a.label {
            if let Some(one) = three_to_one_aa(lbl) {
                tokens.push(one.to_string());
            } else if lbl.starts_with('d') && lbl.len() == 2 && matches!(&lbl[1..], "A" | "C" | "G" | "T") {
                is_nucleic = true;
                tokens.push(format!("[dR]({})P", &lbl[1..]));
            } else if matches!(lbl.as_str(), "A" | "C" | "G" | "U" | "T") {
                is_nucleic = true;
                tokens.push(format!("R({lbl})P"));
            } else {
                tokens.push(if lbl.len() > 1 { format!("[{lbl}]") } else { lbl.clone() });
            }
        }
    }
    if tokens.is_empty() {
        let smi = crate::molecule_to_smiles(mol);
        return format!("CHEM1{{[{smi}]}}$$$$V2.0");
    }
    if is_nucleic {
        format!("RNA1{{{}}}$$$$V2.0", tokens.join("."))
    } else {
        format!("PEPTIDE1{{{}}}$$$$V2.0", tokens.join("."))
    }
}

/// Serializes a biopolymer molecule into **FASTA** format (`>ChiralCraft_Sequence\nACDEF...`).
pub fn molecule_to_fasta(mol: &ChemMolecule, title: &str) -> String {
    let mut seq = String::new();
    for a in &mol.atoms {
        if let Some(lbl) = &a.label {
            if let Some(one) = three_to_one_aa(lbl) {
                seq.push_str(one);
            } else if lbl.starts_with('d') && lbl.len() == 2 {
                seq.push_str(&lbl[1..]);
            } else if lbl.len() == 1 && lbl.chars().all(|c| c.is_ascii_alphabetic()) {
                seq.push_str(&lbl.to_ascii_uppercase());
            } else if let Some(inner) = lbl.strip_prefix("R(").and_then(|s| s.split(')').next()) {
                seq.push_str(inner);
            }
        }
    }
    if seq.is_empty() {
        seq.push('X');
    }
    format!(">{title}\n{seq}\n")
}

/// Cleans up a biopolymer sequence into a neat IUPAC/ChiralCraft serpentine grid (`Structure > Clean Up Biopolymer`).
pub fn clean_up_biopolymer(mol: &mut ChemMolecule, bond_length_pt: f64) {
    let step = (bond_length_pt * 1.65).max(20.0);
    let per_row = 10usize;
    let ox = mol.atoms.first().map(|a| a.x).unwrap_or(120.0);
    let oy = mol.atoms.first().map(|a| a.y).unwrap_or(140.0);
    for (i, a) in mol.atoms.iter_mut().enumerate() {
        let row = i / per_row;
        let col = i % per_row;
        let eff_col = if row % 2 == 0 { col } else { per_row - 1 - col };
        a.x = ox + (eff_col as f64) * step;
        a.y = oy + (row as f64) * (step * 1.45);
    }
}

fn one_to_three_aa(s: &str) -> Option<&'static str> {
    match s {
        "A" => Some("Ala"), "R" => Some("Arg"), "N" => Some("Asn"), "D" => Some("Asp"),
        "C" => Some("Cys"), "E" => Some("Glu"), "Q" => Some("Gln"), "G" => Some("Gly"),
        "H" => Some("His"), "I" => Some("Ile"), "L" => Some("Leu"), "K" => Some("Lys"),
        "M" => Some("Met"), "F" => Some("Phe"), "P" => Some("Pro"), "S" => Some("Ser"),
        "T" => Some("Thr"), "W" => Some("Trp"), "Y" => Some("Tyr"), "V" => Some("Val"),
        _ => None,
    }
}

fn three_to_one_aa(s: &str) -> Option<&'static str> {
    match s {
        "Ala" => Some("A"), "Arg" => Some("R"), "Asn" => Some("N"), "Asp" => Some("D"),
        "Cys" => Some("C"), "Glu" => Some("E"), "Gln" => Some("Q"), "Gly" => Some("G"),
        "His" => Some("H"), "Ile" => Some("I"), "Leu" => Some("L"), "Lys" => Some("K"),
        "Met" => Some("M"), "Phe" => Some("F"), "Pro" => Some("P"), "Ser" => Some("S"),
        "Thr" => Some("T"), "Trp" => Some("W"), "Tyr" => Some("Y"), "Val" => Some("V"),
        _ => None,
    }
}

/// Generates a valid **3MF 3D Model XML** (`Edit > Copy As > 3MF (Ball and Stick)` / `3MF (Stick)`).
pub fn molecule_to_3mf_xml(mol: &ChemMolecule, ball_and_stick: bool) -> String {
    let mut verts: Vec<(f64, f64, f64)> = Vec::new();
    let mut tris: Vec<(usize, usize, usize)> = Vec::new();

    // Scale 14.4 pt -> ~10 mm for 3D printing
    let scale = 10.0 / 14.4;
    let cx = if mol.atoms.is_empty() { 0.0 } else { mol.atoms.iter().map(|a| a.x).sum::<f64>() / (mol.atoms.len() as f64) };
    let cy = if mol.atoms.is_empty() { 0.0 } else { mol.atoms.iter().map(|a| a.y).sum::<f64>() / (mol.atoms.len() as f64) };

    let mut add_octahedron = |x: f64, y: f64, z: f64, r: f64| {
        let b = verts.len();
        verts.push((x + r, y, z));
        verts.push((x - r, y, z));
        verts.push((x, y + r, z));
        verts.push((x, y - r, z));
        verts.push((x, y, z + r));
        verts.push((x, y, z - r));
        for &(v1, v2, v3) in &[
            (0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
            (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5),
        ] {
            tris.push((b + v1, b + v2, b + v3));
        }
    };

    for a in &mol.atoms {
        let ax = (a.x - cx) * scale;
        let ay = -(a.y - cy) * scale;
        let az = a.z_coord * scale;
        let r = if ball_and_stick {
            let cov = get_element(a.z).map(|e| if e.z == 1 { 1.4 } else { 2.4 }).unwrap_or(2.2);
            cov
        } else {
            1.15
        };
        add_octahedron(ax, ay, az, r);
    }

    // Add 6-sided bond prisms
    for b in &mol.bonds {
        if b.begin_atom >= mol.atoms.len() || b.end_atom >= mol.atoms.len() {
            continue;
        }
        let a1 = &mol.atoms[b.begin_atom];
        let a2 = &mol.atoms[b.end_atom];
        let p1 = ((a1.x - cx) * scale, -(a1.y - cy) * scale, a1.z_coord * scale);
        let p2 = ((a2.x - cx) * scale, -(a2.y - cy) * scale, a2.z_coord * scale);
        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let dz = p2.2 - p1.2;
        let len = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-6);
        let ux = dx / len;
        let uy = dy / len;
        let uz = dz / len;
        let (rx, ry, rz) = if uz.abs() < 0.9 { (-uy, ux, 0.0) } else { (1.0, 0.0, 0.0) };
        let rlen = (rx * rx + ry * ry + rz * rz).sqrt().max(1e-6);
        let nx = rx / rlen;
        let ny = ry / rlen;
        let nz = rz / rlen;
        let bx = uy * nz - uz * ny;
        let by = uz * nx - ux * nz;
        let bz = ux * ny - uy * nx;

        let rad = if ball_and_stick { 0.85 } else { 1.10 };
        let base = verts.len();
        for k in 0..6 {
            let ang = (k as f64) * std::f64::consts::TAU / 6.0;
            let (sa, ca) = ang.sin_cos();
            let ox = (nx * ca + bx * sa) * rad;
            let oy = (ny * ca + by * sa) * rad;
            let oz = (nz * ca + bz * sa) * rad;
            verts.push((p1.0 + ox, p1.1 + oy, p1.2 + oz));
            verts.push((p2.0 + ox, p2.1 + oy, p2.2 + oz));
        }
        for k in 0..6 {
            let kn = (k + 1) % 6;
            let v0 = base + k * 2;
            let v1 = base + k * 2 + 1;
            let v2 = base + kn * 2;
            let v3 = base + kn * 2 + 1;
            tris.push((v0, v1, v3));
            tris.push((v0, v3, v2));
        }
    }

    let mut xml = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <model unit=\"millimeter\" xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\">\n\
         <resources><object id=\"1\" type=\"model\"><mesh><vertices>\n",
    );
    for (x, y, z) in &verts {
        xml.push_str(&format!("<vertex x=\"{x:.3}\" y=\"{y:.3}\" z=\"{z:.3}\"/>\n"));
    }
    xml.push_str("</vertices><triangles>\n");
    for (v1, v2, v3) in &tris {
        xml.push_str(&format!("<triangle v1=\"{v1}\" v2=\"{v2}\" v3=\"{v3}\"/>\n"));
    }
    xml.push_str("</triangles></mesh></object></resources>\n<build><item objectid=\"1\"/></build>\n</model>\n");
    xml
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_helm_fasta_monomer_json_and_3mf() {
        let pep = sequence_to_biopolymer_molecule(">TestPeptide\nACDEF", PolymerKind::Peptide, 14.4).unwrap();
        assert_eq!(pep.atoms.len(), 5);
        assert_eq!(pep.atoms[0].label.as_deref(), Some("Ala"));
        assert_eq!(molecule_to_helm(&pep), "PEPTIDE1{A.C.D.E.F}$$$$V2.0");
        assert!(molecule_to_fasta(&pep, "Pep").contains("ACDEF"));

        let cyclic = helm_to_molecule("PEPTIDE1{C.A.G.K.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0", 14.4).unwrap();
        assert_eq!(cyclic.atoms.len(), 5);
        assert_eq!(cyclic.bonds.len(), 5); // 4 backbone + 1 disulfide crosslink!

        let m3 = molecule_to_3mf_xml(&pep, true);
        assert!(m3.contains("<model") && m3.contains("<triangle"));

        let json_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../research/application/natural-monomers.json");
        if json_path.exists() {
            let lib = MonomerLibrary::load_from_json_path(&json_path).unwrap();
            assert!(lib.entries.len() >= 20, "natural-monomers.json should load >=20 monomers, got {}", lib.entries.len());
        }
    }
}
