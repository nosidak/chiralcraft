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
        // One string-aware scan: a monomer is any object that directly owns a `"polymerType":` key.
        // (Indentation/minification independent; nested `rgroups` objects are never mistaken for it.)
        let bytes = text.as_bytes();
        let mut stack: Vec<(usize, bool)> = Vec::new();
        let mut objects: Vec<(usize, usize)> = Vec::new();
        let mut i = 0usize;
        while i < bytes.len() {
            match bytes[i] {
                b'"' => {
                    let Some((s, used)) = decode_json_string(&text[i + 1..]) else { break };
                    let end = i + 1 + used;
                    if s == "polymerType" && text[end..].trim_start().starts_with(':') {
                        if let Some(top) = stack.last_mut() {
                            top.1 = true;
                        }
                    }
                    i = end;
                    continue;
                }
                b'{' => stack.push((i, false)),
                b'}' => {
                    if let Some((start, is_monomer)) = stack.pop() {
                        if is_monomer {
                            objects.push((start, i + 1));
                        }
                    }
                }
                _ => {}
            }
            i += 1;
        }
        for (obj_start, obj_end) in objects {
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

/// Decodes a JSON string body (`s` starts just after the opening quote). Returns the decoded text and
/// the byte length consumed including the closing quote, or `None` if unterminated.
fn decode_json_string(s: &str) -> Option<(String, usize)> {
    fn hex4(it: &mut std::str::CharIndices<'_>) -> Option<u32> {
        let mut v = 0u32;
        for _ in 0..4 {
            v = v * 16 + it.next()?.1.to_digit(16)?;
        }
        Some(v)
    }
    let mut out = String::new();
    let mut it = s.char_indices();
    while let Some((i, ch)) = it.next() {
        match ch {
            '"' => return Some((out, i + 1)),
            '\\' => match it.next()?.1 {
                'n' => out.push('\n'),
                't' => out.push('\t'),
                'r' => out.push('\r'),
                'b' => out.push('\u{8}'),
                'f' => out.push('\u{c}'),
                'u' => {
                    let hi = hex4(&mut it)?;
                    let mut c = char::from_u32(hi);
                    if (0xD800..0xDC00).contains(&hi) {
                        // UTF-16 surrogate pair `\uD83D\uDE00`
                        let mut look = it.clone();
                        c = None;
                        if look.next().map(|x| x.1) == Some('\\') && look.next().map(|x| x.1) == Some('u') {
                            if let Some(lo) = hex4(&mut look).filter(|lo| (0xDC00..0xE000).contains(lo)) {
                                c = char::from_u32(0x10000 + ((hi - 0xD800) << 10) + (lo - 0xDC00));
                                it = look;
                            }
                        }
                    }
                    out.push(c.unwrap_or('\u{FFFD}'));
                }
                other => out.push(other), // `"`, `\`, `/`
            },
            c => out.push(c),
        }
    }
    None
}

/// Returns the string value of a key owned directly by this JSON object (not by nested objects/arrays).
fn extract_json_string(obj: &str, key: &str) -> Option<String> {
    let bytes = obj.as_bytes();
    let mut depth = 0i32;
    let mut i = 0usize;
    while i < bytes.len() {
        match bytes[i] {
            b'"' => {
                let (s, used) = decode_json_string(&obj[i + 1..])?;
                let end = i + 1 + used;
                if depth == 1 && s == key {
                    if let Some(v) = obj[end..].trim_start().strip_prefix(':') {
                        return v.trim_start().strip_prefix('"').and_then(decode_json_string).map(|x| x.0);
                    }
                }
                i = end;
                continue;
            }
            b'{' | b'[' => depth += 1,
            b'}' | b']' => depth -= 1,
            _ => {}
        }
        i += 1;
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
    // First FASTA record only (a second `>` header starts another sequence); `;` lines are comments.
    let mut body_lines: Vec<&str> = Vec::new();
    for l in seq_input.lines().map(str::trim) {
        if l.starts_with('>') {
            if body_lines.iter().any(|b| !b.is_empty()) {
                break;
            }
            continue;
        }
        if !l.starts_with(';') {
            body_lines.push(l);
        }
    }
    let body = body_lines.join(" ");
    let trimmed = body.trim();
    let tokens: Vec<&str> = trimmed
        .split(|c: char| c == '-' || c == '.' || c.is_whitespace())
        .filter(|s| s.chars().any(|c| c.is_ascii_alphanumeric()))
        .collect();
    // Named-monomer mode ("Ala-Cys", "dA-dC", "ALA.GLY") only when every token is a monomer name;
    // otherwise `-`/`.` are FASTA alignment gaps ("AC-DE" = A,C,D,E).
    let is_name = |s: &str| {
        s.chars().count() == 1
            || (kind == PolymerKind::Peptide && canonical_three(s).is_some())
            || (s.chars().any(|c| c.is_ascii_lowercase()) && s.chars().any(|c| c.is_ascii_uppercase()))
    };
    let named_mode = (trimmed.contains('-') || trimmed.contains('.'))
        && tokens.iter().any(|s| s.chars().count() > 1)
        && tokens.iter().all(|s| is_name(s));
    let letters: Vec<String> = if named_mode {
        tokens
            .iter()
            .map(|s| {
                if s.len() == 1 {
                    s.to_ascii_uppercase()
                } else if let Some(c) = canonical_three(s).filter(|_| kind == PolymerKind::Peptide) {
                    c.to_string()
                } else {
                    s.to_string()
                }
            })
            .collect()
    } else {
        trimmed
            .chars()
            .filter(|c| c.is_ascii_alphabetic())
            .map(|c| c.to_ascii_uppercase().to_string())
            .collect()
    };
    if letters.is_empty() {
        return Err("No valid monomer letters found in sequence".to_string());
    }
    build_monomer_chain(&letters, kind, bond_length_pt)
}

/// Parses a HELM v1/v2 string (e.g. `PEPTIDE1{A.C.D.E.F}$$$$V2.0` or `RNA1{R(A)P.R(U)P.R(G)}$$$$V2.0`)
/// into a linked biopolymer `ChemMolecule` (including cyclic/disulfide cross-links if specified).
pub fn helm_to_molecule(helm: &str, bond_length_pt: f64) -> Result<ChemMolecule, String> {
    let trimmed = helm.trim();
    let sections: Vec<&str> = split_helm_top_level(trimmed, '$');
    let poly_part = sections.first().copied().unwrap_or(trimmed);

    let mut mol = ChemMolecule::default();
    let mut poly_atom_map: HashMap<(String, usize), usize> = HashMap::new();
    let mut next_aid = 1u32;
    let mut next_bid = 1000u32;
    let step = (bond_length_pt * 1.65).max(20.0);
    let per_row = 12usize;
    let mut chain_y_offset = 0.0_f64;

    for poly_str in split_helm_top_level(poly_part, '|').into_iter().filter(|s| !s.trim().is_empty()) {
        let brace_start = poly_str
            .find('{')
            .ok_or_else(|| "HELM string must contain `{monomers}`".to_string())?;
        let brace_end = poly_str
            .rfind('}')
            .ok_or_else(|| "HELM string missing closing `}`".to_string())?;
        if brace_end <= brace_start {
            return Err("Malformed HELM braces".to_string());
        }
        let poly_id = poly_str[..brace_start].trim().to_ascii_uppercase();
        let kind = if poly_id.starts_with("RNA") {
            PolymerKind::Rna
        } else if poly_id.starts_with("DNA") {
            PolymerKind::Dna
        } else {
            PolymerKind::Peptide
        };
        let body = &poly_str[brace_start + 1..brace_end];
        let monomers: Vec<String> = split_helm_top_level(body, '.')
            .into_iter()
            .map(|m| strip_helm_brackets(m.trim()).to_string())
            .filter(|m| !m.is_empty())
            .collect();
        if monomers.is_empty() {
            return Err("No monomers found inside HELM `{...}`".to_string());
        }
        let mut prev_idx: Option<usize> = None;
        for (i, sym) in monomers.iter().enumerate() {
            let row = i / per_row;
            let col = i % per_row;
            let eff_col = if row % 2 == 0 { col } else { per_row - 1 - col };
            let mut a = ChemAtom::new(next_aid, 6);
            next_aid = next_aid.saturating_add(1);
            a.x = 120.0 + (eff_col as f64) * step;
            a.y = 140.0 + chain_y_offset + (row as f64) * (step * 1.4);
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
            let idx = mol.add_atom(a);
            poly_atom_map.insert((poly_id.clone(), i + 1), idx);
            if let Some(p) = prev_idx {
                mol.bonds.push(ChemBond {
                    id: next_bid,
                    begin_atom: p,
                    end_atom: idx,
                    order: 1,
                    display: BondDisplay::Solid,
                    visible: true,
                    ignore_warnings: false,
                });
                next_bid = next_bid.saturating_add(1);
            }
            prev_idx = Some(idx);
        }
        let rows = monomers.len().div_ceil(per_row).max(1);
        chain_y_offset += (rows as f64 + 0.8) * (step * 1.4);
    }
    if mol.atoms.is_empty() {
        return Err("No polymers found in HELM string".to_string());
    }
    // Parse optional HELM connection section (e.g. `PEPTIDE1,PEPTIDE1,1:R3-5:R3`)
    if let Some(conn_sec) = sections.get(1) {
        for conn in conn_sec.split('|').filter(|c| !c.trim().is_empty()) {
            let cols: Vec<&str> = conn.split(',').collect();
            if cols.len() >= 3 {
                let p1 = cols[0].trim().to_ascii_uppercase();
                let p2 = cols[1].trim().to_ascii_uppercase();
                if let Some((lhs, rhs)) = cols[2].split_once('-') {
                    let i1 = lhs.split(':').next().and_then(|s| s.trim().parse::<usize>().ok());
                    let i2 = rhs.split(':').next().and_then(|s| s.trim().parse::<usize>().ok());
                    if let (Some(a1), Some(b1)) = (i1, i2) {
                        let idx_a = poly_atom_map.get(&(p1, a1)).copied();
                        let idx_b = poly_atom_map.get(&(p2, b1)).copied();
                        if let (Some(ia), Some(ib)) = (idx_a, idx_b) {
                            let dup = mol.bonds.iter().any(|b| {
                                (b.begin_atom, b.end_atom) == (ia, ib) || (b.begin_atom, b.end_atom) == (ib, ia)
                            });
                            if ia != ib && !dup {
                                mol.bonds.push(ChemBond {
                                    id: next_bid,
                                    begin_atom: ia,
                                    end_atom: ib,
                                    order: 1,
                                    display: BondDisplay::Solid,
                                    visible: true,
                                    ignore_warnings: false,
                                });
                                next_bid = next_bid.saturating_add(1);
                            }
                        }
                    }
                }
            }
        }
    }
    // Bond ids follow atom ids (the fixed 1000 base collided with atom ids for >999 monomers).
    let n_atoms = mol.atoms.len() as u32;
    for (k, b) in mol.bonds.iter_mut().enumerate() {
        b.id = n_atoms + k as u32 + 1;
    }
    Ok(mol)
}

/// Splits on `sep` outside `[...]` (inline SMILES such as `[C.O]` or `[C$C]` stays one monomer).
fn split_helm_top_level(s: &str, sep: char) -> Vec<&str> {
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut start = 0usize;
    for (i, c) in s.char_indices() {
        match c {
            '[' => depth += 1,
            ']' => depth = (depth - 1).max(0),
            _ if c == sep && depth == 0 => {
                out.push(&s[start..i]);
                start = i + c.len_utf8();
            }
            _ => {}
        }
    }
    out.push(&s[start..]);
    out
}

/// `[Orn]` → `Orn`, but `[dR](A)P` stays intact (its first bracket closes before the end).
fn strip_helm_brackets(m: &str) -> &str {
    if !(m.starts_with('[') && m.ends_with(']')) {
        return m;
    }
    let mut depth = 0i32;
    for (i, c) in m.char_indices() {
        match c {
            '[' => depth += 1,
            ']' => {
                depth -= 1;
                if depth == 0 {
                    return if i + 1 == m.len() { &m[1..i] } else { m };
                }
            }
            _ => {}
        }
    }
    m
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
    smiles_to_molecule(&strip_rgroup_atoms(&entry.smiles)).ok()
}

/// Removes HELM R-group attachment atoms (`[*:n]`, `[H:n]` for any n; `[OH:n]` → `O`) and the empty
/// branches `()` they leave behind (`C([*:2])=O` → `C=O`), so the monomer SMILES parses cleanly.
fn strip_rgroup_atoms(smi: &str) -> String {
    let mut out = String::with_capacity(smi.len());
    let mut rest = smi;
    while let Some(p) = rest.find('[') {
        out.push_str(&rest[..p]);
        let after = &rest[p..];
        let Some(q) = after.find(']') else {
            out.push_str(after);
            rest = "";
            break;
        };
        let inner = &after[1..q];
        let (sym, map) = inner.split_once(':').unwrap_or((inner, ""));
        let is_map = !map.is_empty() && map.bytes().all(|b| b.is_ascii_digit());
        match sym {
            "*" | "H" if is_map => {}
            "OH" if is_map => out.push('O'),
            _ => out.push_str(&after[..=q]),
        }
        rest = &after[q + 1..];
    }
    out.push_str(rest);
    while out.contains("()") {
        out = out.replace("()", "");
    }
    out
}

/// Serializes a biopolymer or peptide molecule into **HELM v2.0** notation (`PEPTIDE1{A.C.D}$$$$V2.0`).
pub fn molecule_to_helm(mol: &ChemMolecule) -> String {
    let mut tokens: Vec<String> = Vec::new();
    let mut atom_to_monomer: HashMap<usize, usize> = HashMap::new();
    let mut is_nucleic = false;
    for (ai, a) in mol.atoms.iter().enumerate() {
        if let Some(lbl) = &a.label {
            let tok = if let Some(one) = three_to_one_aa(lbl) {
                one.to_string()
            } else if lbl.starts_with("R(") || lbl.starts_with("[dR](") || lbl.starts_with("dR(") {
                is_nucleic = true;
                lbl.clone()
            } else if lbl.starts_with('d') && lbl.len() == 2 && matches!(&lbl[1..], "A" | "C" | "G" | "T") {
                is_nucleic = true;
                format!("[dR]({})P", &lbl[1..])
            } else if matches!(lbl.as_str(), "A" | "C" | "G" | "U" | "T") {
                is_nucleic = true;
                format!("R({lbl})P")
            } else if lbl.len() > 1 && strip_helm_brackets(lbl).len() == lbl.len() {
                format!("[{lbl}]")
            } else {
                lbl.clone()
            };
            tokens.push(tok);
            atom_to_monomer.insert(ai, tokens.len());
        }
    }
    if tokens.is_empty() {
        let smi = crate::molecule_to_smiles(mol);
        return format!("CHEM1{{[{smi}]}}$$$$V2.0");
    }
    let poly_id = if is_nucleic { "RNA1" } else { "PEPTIDE1" };
    // Preserve non-backbone cross-links (cyclic / disulfide bonds)
    let mut conns: Vec<String> = Vec::new();
    for b in &mol.bonds {
        if let (Some(&m1), Some(&m2)) = (atom_to_monomer.get(&b.begin_atom), atom_to_monomer.get(&b.end_atom)) {
            if m1.abs_diff(m2) > 1 {
                let (lo, hi) = if m1 < m2 { (m1, m2) } else { (m2, m1) };
                conns.push(format!("{poly_id},{poly_id},{lo}:R3-{hi}:R3"));
            }
        }
    }
    format!("{poly_id}{{{}}}${}$$$V2.0", tokens.join("."), conns.join("|"))
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
            } else if let Some(inner) = ["R(", "dR(", "[dR]("]
                .iter()
                .find_map(|p| lbl.strip_prefix(p))
                .and_then(|s| s.split(')').next())
            {
                seq.push_str(inner);
            }
        }
    }
    if seq.is_empty() {
        seq.push('X');
    }
    // A newline in the title would turn the rest of it into sequence data.
    let title = title.replace(['\r', '\n'], " ");
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

/// Case-insensitive 3-letter amino-acid code → canonical "Ala" spelling.
fn canonical_three(s: &str) -> Option<&'static str> {
    if s.len() != 3 {
        return None;
    }
    "ARNDCEQGHILKMFPSTWYV"
        .chars()
        .filter_map(|c| one_to_three_aa(c.encode_utf8(&mut [0u8; 4])))
        .find(|t| t.eq_ignore_ascii_case(s))
}

fn three_to_one_aa(s: &str) -> Option<&'static str> {
    let s = canonical_three(s)?;
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

    #[test]
    fn test_biopolymer_audit_regressions() {
        // FASTA alignment gaps are not monomer separators; only the first record; `;` comments skipped.
        let g = sequence_to_biopolymer_molecule(">s1\n;comment\nAC-DE\n>s2\nWWW", PolymerKind::Peptide, 14.4).unwrap();
        assert_eq!(g.atoms.len(), 4);
        assert_eq!(molecule_to_fasta(&g, "a\nb"), ">a b\nACDE\n");
        // Named mode, case-insensitive 3-letter codes.
        let n = sequence_to_biopolymer_molecule("ALA-gly-Cys", PolymerKind::Peptide, 14.4).unwrap();
        assert_eq!(molecule_to_helm(&n), "PEPTIDE1{A.G.C}$$$$V2.0");
        // DNA round trip: writer emits `[dR](A)P`, reader must not mangle it.
        let d = sequence_to_biopolymer_molecule("ACG", PolymerKind::Dna, 14.4).unwrap();
        let helm = molecule_to_helm(&d);
        assert_eq!(helm, "RNA1{[dR](A)P.[dR](C)P.[dR](G)P}$$$$V2.0");
        let back = helm_to_molecule(&helm, 14.4).unwrap();
        assert_eq!(back.atoms[0].label.as_deref(), Some("[dR](A)P"));
        assert_eq!(molecule_to_helm(&back), helm);
        assert!(molecule_to_fasta(&back, "x").contains("ACG"));
        // Bracketed names stripped once; inline SMILES with '.'/'$' kept whole; duplicate link ignored.
        let h = helm_to_molecule("PEPTIDE1{[Orn].[C.O].[C$C]}$PEPTIDE1,PEPTIDE1,1:R2-2:R1$$$V2.0", 14.4).unwrap();
        assert_eq!(h.atoms.len(), 3);
        assert_eq!(h.atoms[0].label.as_deref(), Some("Orn"));
        assert_eq!(h.atoms[1].label.as_deref(), Some("C.O"));
        assert_eq!(h.bonds.len(), 2);
        // Bond ids never collide with atom ids.
        let long = format!("PEPTIDE1{{{}}}$$$$V2.0", vec!["A"; 1200].join("."));
        let lm = helm_to_molecule(&long, 14.4).unwrap();
        let ids: std::collections::HashSet<u32> = lm.atoms.iter().map(|a| a.id).chain(lm.bonds.iter().map(|b| b.id)).collect();
        assert_eq!(ids.len(), lm.atoms.len() + lm.bonds.len());
        // JSON: escapes, nested objects, minified layout.
        let lib = MonomerLibrary::parse_json(
            r#"[{"rgroups":[{"symbol":"R1"}],"symbol":"X\u00e9\/1","polymerType":"PEPTIDE","smiles":"C([*:3])[*:1]","name":"a\"b\ud83d\ude00"}]"#,
        );
        assert_eq!(lib.entries.len(), 1);
        assert_eq!(lib.entries[0].symbol, "Xé/1");
        assert_eq!(lib.entries[0].name, "a\"b😀");
        assert_eq!(strip_rgroup_atoms(&lib.entries[0].smiles), "C");
        assert_eq!(strip_rgroup_atoms("[*:1]N[C@@H](C)C([*:2])=O"), "N[C@@H](C)C=O");
        assert_eq!(strip_rgroup_atoms("C(=O)[OH:2]"), "C(=O)O");
    }
}
