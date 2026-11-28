//! Live PubChem (PUG-REST) & ChEMBL (EBI REST) Chemical Database Search & Structure Import
//!
//! Zero external HTTP/JSON crate dependencies:
//! - Minimal recursive-descent JSON parser with full `\uXXXX` and escape support
//! - Non-blocking background queries via system `curl` (`std::thread` + `std::sync::mpsc`)
//! - Supports PubChem & ChEMBL search by Name / CAS / Formula / CID / ChEMBL ID, Exact SMILES,
//!   2D Tanimoto Similarity (configurable threshold), and Substructure
//! - Curated 2D MOL/SDF import with IUPAC skeletal hydrogen pruning (preserving chiral H & isotopes)
//!   and automatic bond-length scaling to the active document stylesheet (`14.4 pt` ACS, `10.8 pt` Nature, etc.)

use std::collections::{BTreeMap, HashMap, HashSet};
use std::process::Command;
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;

use chiralcraft_chem::{get_element, BondDisplay, CcVal, ChemMolecule};
use chiralcraft_doc::clean_up_structure_preserving_stereo;
use chiralcraft_formats::{molfile_to_molecule, molfile_v3000_to_molecule, smiles_to_molecule};
use egui::{Color32, FontId, Pos2, Rect, RichText, Sense, Stroke, Vec2};

use crate::chem_ops::cpk_element_color;
use crate::palettes::CD_WINDOW_BG;

// ============================================================================
// 1. Minimal Zero-Dependency JSON Parser
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub enum JsonValue {
    Null,
    Bool(bool),
    Number(f64),
    String(String),
    Array(Vec<JsonValue>),
    Object(BTreeMap<String, JsonValue>),
}

impl JsonValue {
    pub fn parse(input: &str) -> Result<Self, String> {
        let mut p = JsonParser {
            bytes: input.as_bytes(),
            pos: 0,
            depth: 0,
        };
        let v = p.parse_value()?;
        p.skip_ws();
        if p.pos < p.bytes.len() {
            return Err(format!("Trailing characters in JSON at byte {}", p.pos));
        }
        Ok(v)
    }

    pub fn get(&self, key: &str) -> Option<&JsonValue> {
        match self {
            JsonValue::Object(map) => map.get(key),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            JsonValue::String(s) => Some(s.as_str()),
            _ => None,
        }
    }

    pub fn as_string_lossy(&self) -> Option<String> {
        match self {
            JsonValue::String(s) => Some(s.clone()),
            JsonValue::Number(n) => {
                if (n.fract()).abs() < 1e-9 && n.abs() < 1e15 {
                    Some(format!("{}", *n as i64))
                } else {
                    Some(format!("{n:.2}"))
                }
            }
            JsonValue::Bool(b) => Some(b.to_string()),
            _ => None,
        }
    }

    pub fn as_u64(&self) -> Option<u64> {
        match self {
            JsonValue::Number(n) if *n >= 0.0 => Some(*n as u64),
            JsonValue::String(s) => s.trim().parse::<u64>().ok(),
            _ => None,
        }
    }

    pub fn as_i64(&self) -> Option<i64> {
        match self {
            JsonValue::Number(n) => Some(*n as i64),
            JsonValue::String(s) => s.trim().parse::<f64>().ok().map(|f| f as i64),
            _ => None,
        }
    }

    pub fn as_array(&self) -> Option<&[JsonValue]> {
        match self {
            JsonValue::Array(v) => Some(v.as_slice()),
            _ => None,
        }
    }
}

struct JsonParser<'a> {
    bytes: &'a [u8],
    pos: usize,
    depth: usize,
}

impl<'a> JsonParser<'a> {
    const MAX_DEPTH: usize = 128;

    fn skip_ws(&mut self) {
        while self.pos < self.bytes.len() {
            match self.bytes[self.pos] {
                b' ' | b'\t' | b'\r' | b'\n' => self.pos += 1,
                _ => break,
            }
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn parse_value(&mut self) -> Result<JsonValue, String> {
        self.skip_ws();
        match self.peek() {
            Some(b'{') => self.parse_object(),
            Some(b'[') => self.parse_array(),
            Some(b'"') => self.parse_string().map(JsonValue::String),
            Some(b't') => self.expect_lit(b"true", JsonValue::Bool(true)),
            Some(b'f') => self.expect_lit(b"false", JsonValue::Bool(false)),
            Some(b'n') => self.expect_lit(b"null", JsonValue::Null),
            Some(c) if c == b'-' || c.is_ascii_digit() => self.parse_number(),
            Some(c) => Err(format!("Unexpected JSON char '{}' at {}", c as char, self.pos)),
            None => Err("Unexpected end of JSON".to_string()),
        }
    }

    fn expect_lit(&mut self, lit: &[u8], val: JsonValue) -> Result<JsonValue, String> {
        if self.bytes.get(self.pos..self.pos + lit.len()) == Some(lit) {
            self.pos += lit.len();
            Ok(val)
        } else {
            Err(format!("Invalid JSON literal at {}", self.pos))
        }
    }

    fn parse_object(&mut self) -> Result<JsonValue, String> {
        if self.depth >= Self::MAX_DEPTH {
            return Err("JSON nesting depth exceeded".to_string());
        }
        self.depth += 1;
        self.pos += 1; // '{'
        let mut map = BTreeMap::new();
        self.skip_ws();
        if self.peek() == Some(b'}') {
            self.pos += 1;
            self.depth -= 1;
            return Ok(JsonValue::Object(map));
        }
        loop {
            self.skip_ws();
            if self.peek() != Some(b'"') {
                self.depth -= 1;
                return Err(format!("Expected JSON string key at {}", self.pos));
            }
            let key = self.parse_string()?;
            self.skip_ws();
            if self.peek() != Some(b':') {
                self.depth -= 1;
                return Err(format!("Expected ':' at {}", self.pos));
            }
            self.pos += 1;
            let val = self.parse_value()?;
            map.insert(key, val);
            self.skip_ws();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b'}') => {
                    self.pos += 1;
                    self.depth -= 1;
                    return Ok(JsonValue::Object(map));
                }
                _ => {
                    self.depth -= 1;
                    return Err(format!("Expected ',' or '}}' at {}", self.pos));
                }
            }
        }
    }

    fn parse_array(&mut self) -> Result<JsonValue, String> {
        if self.depth >= Self::MAX_DEPTH {
            return Err("JSON nesting depth exceeded".to_string());
        }
        self.depth += 1;
        self.pos += 1; // '['
        let mut arr = Vec::new();
        self.skip_ws();
        if self.peek() == Some(b']') {
            self.pos += 1;
            self.depth -= 1;
            return Ok(JsonValue::Array(arr));
        }
        loop {
            let v = self.parse_value()?;
            arr.push(v);
            self.skip_ws();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b']') => {
                    self.pos += 1;
                    self.depth -= 1;
                    return Ok(JsonValue::Array(arr));
                }
                _ => {
                    self.depth -= 1;
                    return Err(format!("Expected ',' or ']' at {}", self.pos));
                }
            }
        }
    }

    fn parse_string(&mut self) -> Result<String, String> {
        self.pos += 1; // '"'
        let mut out = String::new();
        while self.pos < self.bytes.len() {
            let b = self.bytes[self.pos];
            self.pos += 1;
            match b {
                b'"' => return Ok(out),
                b'\\' => {
                    if self.pos >= self.bytes.len() {
                        return Err("Unterminated escape in JSON string".to_string());
                    }
                    let esc = self.bytes[self.pos];
                    self.pos += 1;
                    match esc {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'b' => out.push('\u{0008}'),
                        b'f' => out.push('\u{000C}'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'u' => {
                            let hex = self
                                .bytes
                                .get(self.pos..self.pos + 4)
                                .ok_or_else(|| "Truncated \\u escape".to_string())?;
                            let hex_s = std::str::from_utf8(hex).map_err(|_| "Bad \\u UTF-8")?;
                            // from_str_radix alone would accept a leading '+' (e.g. "\u+041")
                            if !hex.iter().all(|c| c.is_ascii_hexdigit()) {
                                return Err(format!("Bad \\u hex {hex_s}"));
                            }
                            let code = u32::from_str_radix(hex_s, 16)
                                .map_err(|_| format!("Bad \\u hex {hex_s}"))?;
                            self.pos += 4;
                            // Handle UTF-16 surrogate pairs
                            if (0xD800..=0xDBFF).contains(&code)
                                && self.bytes.get(self.pos..self.pos + 2) == Some(b"\\u")
                            {
                                if let Some(lo_hex) = self.bytes.get(self.pos + 2..self.pos + 6) {
                                    if let Ok(lo_s) = std::str::from_utf8(lo_hex) {
                                        if let Ok(lo) = u32::from_str_radix(lo_s, 16) {
                                            if (0xDC00..=0xDFFF).contains(&lo) {
                                                self.pos += 6;
                                                let cp = 0x10000
                                                    + (((code - 0xD800) << 10) | (lo - 0xDC00));
                                                if let Some(ch) = char::from_u32(cp) {
                                                    out.push(ch);
                                                    continue;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            // Unpaired surrogates are not scalar values: substitute U+FFFD instead of dropping them
                            out.push(char::from_u32(code).unwrap_or('\u{FFFD}'));
                        }
                        other => out.push(other as char),
                    }
                }
                _ if b < 0x80 => out.push(b as char),
                _ => {
                    // Multi-byte UTF-8 sequence starting at self.pos - 1 (O(1) decode of 1 codepoint)
                    let start = self.pos - 1;
                    let expected_len = if (b & 0xE0) == 0xC0 {
                        2
                    } else if (b & 0xF0) == 0xE0 {
                        3
                    } else if (b & 0xF8) == 0xF0 {
                        4
                    } else {
                        return Err("Invalid UTF-8 lead byte in JSON string".to_string());
                    };
                    let slice = self
                        .bytes
                        .get(start..start + expected_len)
                        .ok_or_else(|| "Truncated UTF-8 in JSON string".to_string())?;
                    let s = std::str::from_utf8(slice)
                        .map_err(|_| "Invalid UTF-8 in JSON string".to_string())?;
                    let ch = s
                        .chars()
                        .next()
                        .ok_or_else(|| "Invalid UTF-8 in JSON string".to_string())?;
                    out.push(ch);
                    self.pos = start + expected_len;
                }
            }
        }
        Err("Unterminated JSON string".to_string())
    }

    fn parse_number(&mut self) -> Result<JsonValue, String> {
        let start = self.pos;
        if self.peek() == Some(b'-') {
            self.pos += 1;
        }
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() || matches!(c, b'.' | b'e' | b'E' | b'+' | b'-') {
                self.pos += 1;
            } else {
                break;
            }
        }
        let s = std::str::from_utf8(&self.bytes[start..self.pos])
            .map_err(|_| "Invalid number UTF-8".to_string())?;
        let n = s
            .parse::<f64>()
            .map_err(|_| format!("Invalid JSON number '{s}'"))?;
        if !n.is_finite() {
            return Err(format!("Non-finite JSON number '{s}'"));
        }
        Ok(JsonValue::Number(n))
    }
}

// ============================================================================
// 2. Database Types, URL Encoding, & Parsers
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DbSourceFilter {
    Both,
    PubChem,
    ChEMBL,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DbSource {
    PubChem,
    ChEMBL,
}

impl DbSource {
    pub fn label(self) -> &'static str {
        match self {
            DbSource::PubChem => "PubChem",
            DbSource::ChEMBL => "ChEMBL",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DbSearchMode {
    NameOrId,
    ExactSmiles,
    Similarity2D,
    Substructure,
}

impl DbSearchMode {
    pub fn label(self) -> &'static str {
        match self {
            DbSearchMode::NameOrId => "Name / CAS / ID",
            DbSearchMode::ExactSmiles => "Exact SMILES",
            DbSearchMode::Similarity2D => "2D Similarity",
            DbSearchMode::Substructure => "Substructure",
        }
    }
}

#[derive(Debug, Clone)]
pub struct DbHit {
    pub source: DbSource,
    pub id: String,
    pub cid: Option<u64>,
    pub name: String,
    pub iupac_name: Option<String>,
    pub formula: String,
    pub molecular_weight: String,
    pub smiles: String,
    pub inchikey: Option<String>,
    pub max_phase: Option<String>,
    pub first_approval: Option<i64>,
    pub alogp: Option<String>,
    pub psa: Option<String>,
    pub ro5_violations: Option<i64>,
    pub similarity_pct: Option<String>,
    pub molfile: Option<String>,
    pub preview_mol: Option<ChemMolecule>,
    pub sdf_loading: bool,
}

impl DbHit {
    pub fn format_caption(&self) -> String {
        let mut lines = Vec::new();
        let display_name = if !self.name.is_empty() {
            self.name.clone()
        } else {
            self.iupac_name.clone().unwrap_or_else(|| self.id.clone())
        };
        lines.push(format!("{} ({})", display_name, self.id));
        if !self.formula.is_empty() || !self.molecular_weight.is_empty() {
            lines.push(format!(
                "Formula: {}   MW: {}",
                if self.formula.is_empty() { "—" } else { &self.formula },
                if self.molecular_weight.is_empty() { "—" } else { &self.molecular_weight }
            ));
        }
        if let Some(ref ph) = self.max_phase {
            if let Some(yr) = self.first_approval {
                lines.push(format!("Clinical Status: {ph} (First Approved: {yr})"));
            } else {
                lines.push(format!("Clinical Status: {ph}"));
            }
        }
        lines.join("\n")
    }
}

pub fn url_encode(s: &str) -> String {
    let mut out = String::with_capacity(s.len() * 2);
    for &b in s.as_bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~') {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{:02X}", b));
        }
    }
    out
}

/// Remove non-essential explicit hydrogens from curated 2D MOLfiles (like PubChem SDF)
/// so imported structures adhere to clean IUPAC skeletal conventions.
/// Preserves explicit `H` when:
/// - attached via a stereobond (`WedgeBegin`, `WedgedHashBegin`, `Wavy`, etc.)
/// - non-standard isotope (`isotope != 0`, e.g. D/T)
/// - nonzero formal charge or radical
/// - unbonded or bonded solely to other H atoms (e.g., `H₂`)
pub fn strip_nonessential_explicit_hydrogens(mol: &mut ChemMolecule) {
    let n = mol.atoms.len();
    if n <= 1 {
        return;
    }
    let mut keep = vec![true; n];
    for (ai, atom) in mol.atoms.iter().enumerate() {
        if atom.z != 1 || atom.charge != 0 || atom.isotope != 0 || atom.radical != 0 {
            continue;
        }
        let mut incident = Vec::new();
        for b in &mol.bonds {
            if b.begin_atom == ai {
                incident.push((b.end_atom, b.order, b.display));
            } else if b.end_atom == ai {
                incident.push((b.begin_atom, b.order, b.display));
            }
        }
        // Only prune single-bonded terminal H attached to a heavy atom via a plain Solid bond
        if incident.len() == 1 {
            let (nbr, order, disp) = incident[0];
            let nbr_is_heavy = mol.atoms.get(nbr).map(|a| a.z > 1).unwrap_or(false);
            let is_plain_single = order == 1 && matches!(disp, BondDisplay::Solid);
            if nbr_is_heavy && is_plain_single {
                keep[ai] = false;
            }
        }
    }
    if keep.iter().all(|&k| k) {
        return;
    }
    let mut remap: Vec<Option<usize>> = vec![None; n];
    let mut new_atoms = Vec::with_capacity(n);
    for (old_i, atom) in mol.atoms.iter().cloned().enumerate() {
        if keep[old_i] {
            remap[old_i] = Some(new_atoms.len());
            new_atoms.push(atom);
        }
    }
    let mut new_bonds = Vec::with_capacity(mol.bonds.len());
    for b in &mol.bonds {
        if let (Some(na), Some(nb)) = (
            remap.get(b.begin_atom).copied().flatten(),
            remap.get(b.end_atom).copied().flatten(),
        ) {
            let mut nbond = b.clone();
            nbond.begin_atom = na;
            nbond.end_atom = nb;
            new_bonds.push(nbond);
        }
    }
    mol.atoms = new_atoms;
    mol.bonds = new_bonds;
}

/// Normalize a molecule's 2D coordinates so its median bond length matches `target_bond_length_pt`
/// (e.g., `14.4 pt` for ACS Document 1996, `10.8 pt` for Nature Research) and center it around `(cx, cy)`.
pub fn normalize_molecule_to_stylesheet(
    mol: &mut ChemMolecule,
    target_bond_length_pt: f64,
    center_x: f64,
    center_y: f64,
) {
    if mol.atoms.is_empty() {
        return;
    }
    let mut lens: Vec<f64> = mol
        .bonds
        .iter()
        .filter_map(|b| {
            let a1 = mol.atoms.get(b.begin_atom)?;
            let a2 = mol.atoms.get(b.end_atom)?;
            let d = (a2.x - a1.x).hypot(a2.y - a1.y);
            if d > 1e-4 { Some(d) } else { None }
        })
        .collect();
    if !lens.is_empty() {
        lens.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median = lens[lens.len() / 2];
        if median > 1e-4 {
            let scale = target_bond_length_pt / median;
            for a in &mut mol.atoms {
                a.x *= scale;
                a.y *= scale;
            }
        }
    }
    let (mut x1, mut y1, mut x2, mut y2) = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for a in &mol.atoms {
        x1 = x1.min(a.x);
        y1 = y1.min(a.y);
        x2 = x2.max(a.x);
        y2 = y2.max(a.y);
    }
    if x1.is_finite() {
        let cur_cx = 0.5 * (x1 + x2);
        let cur_cy = 0.5 * (y1 + y2);
        let dx = center_x - cur_cx;
        let dy = center_y - cur_cy;
        for a in &mut mol.atoms {
            a.x += dx;
            a.y += dy;
        }
    }
}

/// Build a clean, stylesheet-scaled `ChemMolecule` from a `DbHit` (preferring curated 2D MOLfile,
/// falling back to SMILES + 2D coordinate cleanup).
pub fn hit_to_molecule(
    hit: &DbHit,
    bond_length_pt: f64,
    center_x: f64,
    center_y: f64,
) -> Result<ChemMolecule, String> {
    let mut mol = if let Some(ref mf) = hit.molfile {
        let parsed = if mf.contains("V3000") {
            molfile_v3000_to_molecule(mf)
        } else {
            molfile_to_molecule(mf)
        };
        match parsed {
            Ok(mut m) if !m.atoms.is_empty() => {
                strip_nonessential_explicit_hydrogens(&mut m);
                normalize_molecule_to_stylesheet(&mut m, bond_length_pt, center_x, center_y);
                m
            }
            _ => {
                if hit.smiles.is_empty() {
                    return Err("MOLfile parse failed and no SMILES available".to_string());
                }
                let mut m = smiles_to_molecule(&hit.smiles)?;
                clean_up_structure_preserving_stereo(&mut m, &[], bond_length_pt);
                normalize_molecule_to_stylesheet(&mut m, bond_length_pt, center_x, center_y);
                m
            }
        }
    } else if !hit.smiles.is_empty() {
        let mut m = smiles_to_molecule(&hit.smiles)?;
        clean_up_structure_preserving_stereo(&mut m, &[], bond_length_pt);
        normalize_molecule_to_stylesheet(&mut m, bond_length_pt, center_x, center_y);
        m
    } else {
        return Err("No MOLfile or SMILES available for this compound".to_string());
    };
    // Ensure all atom/bond visibility flags are active
    for a in &mut mol.atoms {
        a.visible = true;
    }
    for b in &mut mol.bonds {
        b.visible = true;
    }
    Ok(mol)
}

fn format_chembl_phase(raw: Option<&JsonValue>) -> Option<String> {
    let s = raw?.as_string_lossy()?;
    let v: f64 = s.trim().parse().ok()?;
    if v >= 3.9 {
        Some("Approved (Phase IV)".to_string())
    } else if v >= 2.9 {
        Some("Phase III".to_string())
    } else if v >= 1.9 {
        Some("Phase II".to_string())
    } else if v >= 0.9 {
        Some("Phase I".to_string())
    } else if v > 0.0 {
        Some("Early Phase I".to_string())
    } else if v < 0.0 {
        // ChEMBL uses -1 for "clinical phase unknown"
        Some("Unknown Phase".to_string())
    } else {
        Some("Preclinical".to_string())
    }
}

/// Parse PubChem PUG-REST `PropertyTable` JSON into `Vec<DbHit>`.
pub fn parse_pubchem_properties_json(body: &str) -> Result<Vec<DbHit>, String> {
    let root = JsonValue::parse(body)?;
    if let Some(fault) = root.get("Fault") {
        let msg = fault
            .get("Message")
            .and_then(|v| v.as_str())
            .unwrap_or("PubChem query failed");
        let details = fault
            .get("Details")
            .and_then(|v| v.as_array())
            .and_then(|a| a.first())
            .and_then(|v| v.as_str());
        return Err(match details {
            Some(d) => format!("PubChem: {msg} ({d})"),
            None => format!("PubChem: {msg}"),
        });
    }
    let props = root
        .get("PropertyTable")
        .and_then(|t| t.get("Properties"))
        .and_then(|p| p.as_array())
        .ok_or_else(|| "PubChem: No properties returned".to_string())?;

    let mut hits = Vec::new();
    for item in props {
        let cid = item.get("CID").and_then(|v| v.as_u64());
        let id = cid
            .map(|c| format!("CID {c}"))
            .unwrap_or_else(|| "CID ?".to_string());
        let title = item
            .get("Title")
            .and_then(|v| v.as_str())
            .unwrap_or("")
            .to_string();
        let iupac = item
            .get("IUPACName")
            .and_then(|v| v.as_str())
            .map(|s| s.to_string());
        let formula = item
            .get("MolecularFormula")
            .and_then(|v| v.as_str())
            .unwrap_or("")
            .to_string();
        let mw = item
            .get("MolecularWeight")
            .and_then(|v| v.as_string_lossy())
            .unwrap_or_default();
        let smiles = item
            .get("SMILES")
            .or_else(|| item.get("IsomericSMILES"))
            .or_else(|| item.get("CanonicalSMILES"))
            .or_else(|| item.get("ConnectivitySMILES"))
            .and_then(|v| v.as_str())
            .unwrap_or("")
            .to_string();
        let inchikey = item
            .get("InChIKey")
            .and_then(|v| v.as_str())
            .map(|s| s.to_string());
        let name = if !title.is_empty() {
            title
        } else {
            iupac.clone().unwrap_or_else(|| id.clone())
        };
        let mut hit = DbHit {
            source: DbSource::PubChem,
            id,
            cid,
            name,
            iupac_name: iupac,
            formula,
            molecular_weight: mw,
            smiles,
            inchikey,
            max_phase: None,
            first_approval: None,
            alogp: None,
            psa: None,
            ro5_violations: None,
            similarity_pct: None,
            molfile: None,
            preview_mol: None,
            sdf_loading: false,
        };
        // Pre-build a fast SMILES preview until curated 2D SDF arrives
        if !hit.smiles.is_empty() {
            if let Ok(m) = hit_to_molecule(&hit, 14.4, 0.0, 0.0) {
                hit.preview_mol = Some(m);
            }
        }
        hits.push(hit);
    }
    Ok(hits)
}

/// Parse PubChem autocomplete JSON (`dictionary_terms.compound`).
pub fn parse_pubchem_autocomplete_json(body: &str) -> Vec<String> {
    let Ok(root) = JsonValue::parse(body) else {
        return Vec::new();
    };
    root.get("dictionary_terms")
        .and_then(|d| d.get("compound"))
        .and_then(|a| a.as_array())
        .map(|arr| {
            arr.iter()
                .filter_map(|v| v.as_str().map(|s| s.to_string()))
                .collect()
        })
        .unwrap_or_default()
}

/// Parse a single ChEMBL molecule JSON object into `DbHit`.
fn parse_chembl_molecule_obj(item: &JsonValue) -> Option<DbHit> {
    let chembl_id = item
        .get("molecule_chembl_id")
        .and_then(|v| v.as_str())?
        .to_string();
    let pref_name = item
        .get("pref_name")
        .and_then(|v| v.as_str())
        .unwrap_or("")
        .to_string();
    let structs = item.get("molecule_structures");
    let smiles = structs
        .and_then(|s| s.get("canonical_smiles"))
        .and_then(|v| v.as_str())
        .unwrap_or("")
        .to_string();
    let molfile = structs
        .and_then(|s| s.get("molfile"))
        .and_then(|v| v.as_str())
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.to_string());
    let inchikey = structs
        .and_then(|s| s.get("standard_inchi_key"))
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());
    let props = item.get("molecule_properties");
    let formula = props
        .and_then(|p| p.get("full_molformula"))
        .and_then(|v| v.as_str())
        .unwrap_or("")
        .to_string();
    let mw = props
        .and_then(|p| p.get("full_mwt").or_else(|| p.get("mw_freebase")))
        .and_then(|v| v.as_string_lossy())
        .unwrap_or_default();
    let alogp = props
        .and_then(|p| p.get("alogp"))
        .and_then(|v| v.as_string_lossy());
    let psa = props
        .and_then(|p| p.get("psa"))
        .and_then(|v| v.as_string_lossy());
    let ro5 = props
        .and_then(|p| p.get("num_ro5_violations"))
        .and_then(|v| v.as_i64());
    let max_phase = format_chembl_phase(item.get("max_phase"));
    let first_approval = item.get("first_approval").and_then(|v| v.as_i64());
    let similarity_pct = item
        .get("similarity")
        .and_then(|v| v.as_string_lossy())
        .map(|s| format!("{s}%"));

    let display_name = if !pref_name.is_empty() {
        pref_name
    } else {
        chembl_id.clone()
    };

    let mut hit = DbHit {
        source: DbSource::ChEMBL,
        id: chembl_id,
        cid: None,
        name: display_name,
        iupac_name: None,
        formula,
        molecular_weight: mw,
        smiles,
        inchikey,
        max_phase,
        first_approval,
        alogp,
        psa,
        ro5_violations: ro5,
        similarity_pct,
        molfile,
        preview_mol: None,
        sdf_loading: false,
    };
    if let Ok(m) = hit_to_molecule(&hit, 14.4, 0.0, 0.0) {
        hit.preview_mol = Some(m);
    }
    Some(hit)
}

/// Parse ChEMBL search/similarity/substructure/single-molecule JSON response into `Vec<DbHit>`.
pub fn parse_chembl_molecules_json(body: &str) -> Result<Vec<DbHit>, String> {
    let root = JsonValue::parse(body)?;
    if let Some(arr) = root.get("molecules").and_then(|v| v.as_array()) {
        let mut hits = Vec::new();
        for item in arr {
            if let Some(h) = parse_chembl_molecule_obj(item) {
                hits.push(h);
            }
        }
        return Ok(hits);
    }
    // Single-molecule endpoint (`/molecule/CHEMBL25.json`)
    if let Some(h) = parse_chembl_molecule_obj(&root) {
        return Ok(vec![h]);
    }
    Err("ChEMBL: Unexpected response format".to_string())
}

// ============================================================================
// 3. Non-Blocking HTTP Client (via system `curl`)
// ============================================================================

fn http_get(url: &str, accept: &str, timeout_secs: u64) -> Result<String, String> {
    let output = Command::new("curl")
        .arg("-sS")
        .arg("-L")
        .arg("--max-time")
        .arg(timeout_secs.to_string())
        .arg("-H")
        .arg(format!("Accept: {accept}"))
        .arg("-H")
        .arg("User-Agent: ChiralCraft/1.0")
        .arg("-w")
        .arg("\n__HTTP_STATUS__:%{http_code}")
        .arg(url)
        .output()
        .map_err(|e| format!("Failed to invoke curl: {e}"))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("Network error: {}", stderr.trim()));
    }
    let raw = String::from_utf8_lossy(&output.stdout);
    let (body, code) = if let Some(idx) = raw.rfind("\n__HTTP_STATUS__:") {
        let code_str = raw[idx + "\n__HTTP_STATUS__:".len()..].trim();
        let code = code_str.parse::<u16>().unwrap_or(0);
        (&raw[..idx], code)
    } else {
        (raw.as_ref(), 200)
    };
    if code == 404 {
        return Err("Not found (HTTP 404)".to_string());
    }
    if code >= 400 {
        return Err(format!("HTTP {code}: {}", body.chars().take(160).collect::<String>()));
    }
    Ok(body.to_string())
}

const PUBCHEM_PROPS: &str =
    "Title,IUPACName,MolecularFormula,MolecularWeight,SMILES,ConnectivitySMILES,InChIKey";

fn is_inchikey(q: &str) -> bool {
    let b = q.as_bytes();
    b.len() == 27
        && b[14] == b'-'
        && b[25] == b'-'
        && b[0..14].iter().all(|c| c.is_ascii_uppercase())
        && b[15..25].iter().all(|c| c.is_ascii_uppercase())
        && b[26].is_ascii_uppercase()
}

fn is_formula_like(q: &str) -> bool {
    if q.len() < 2 || q.len() > 24 || q.contains(' ') {
        return false;
    }
    let has_upper = q.chars().any(|c| c.is_ascii_uppercase());
    let has_digit = q.chars().any(|c| c.is_ascii_digit());
    let all_ok = q.chars().all(|c| c.is_ascii_alphanumeric());
    // Avoid mistaking simple names with digits like "VitaminB12" or "CHEMBL25"
    has_upper
        && has_digit
        && all_ok
        && !q.to_ascii_uppercase().starts_with("CHEMBL")
        && !q.to_ascii_uppercase().starts_with("CID")
        && q.chars().filter(|c| c.is_ascii_lowercase()).count() <= 3
}

fn fetch_pubchem_sdf(cid: u64) -> Result<String, String> {
    let url = format!("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF");
    http_get(&url, "chemical/x-mdl-sdfile", 12)
}

pub fn query_pubchem_public(query: &str, mode: DbSearchMode, sim_threshold: u32) -> Result<Vec<DbHit>, String> {
    query_pubchem(query, mode, sim_threshold)
}

pub fn query_chembl_public(query: &str, mode: DbSearchMode, sim_threshold: u32) -> Result<Vec<DbHit>, String> {
    query_chembl(query, mode, sim_threshold)
}

fn query_pubchem(query: &str, mode: DbSearchMode, sim_threshold: u32) -> Result<Vec<DbHit>, String> {
    let q = query.trim();
    if q.is_empty() {
        return Ok(Vec::new());
    }
    let enc = url_encode(q);
    let mut hits = match mode {
        DbSearchMode::NameOrId => {
            let stripped_cid = q
                .strip_prefix("CID:")
                .or_else(|| q.strip_prefix("cid:"))
                .or_else(|| q.strip_prefix("CID "))
                .or_else(|| q.strip_prefix("cid "))
                .unwrap_or(q)
                .trim();
            if !stripped_cid.is_empty() && stripped_cid.chars().all(|c| c.is_ascii_digit()) {
                let url = format!(
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{stripped_cid}/property/{PUBCHEM_PROPS}/JSON"
                );
                let body = http_get(&url, "application/json", 12)?;
                parse_pubchem_properties_json(&body)?
            } else if is_inchikey(q) {
                let url = format!(
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{enc}/property/{PUBCHEM_PROPS}/JSON"
                );
                let body = http_get(&url, "application/json", 12)?;
                parse_pubchem_properties_json(&body)?
            } else if is_formula_like(q) {
                let url = format!(
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{enc}/property/{PUBCHEM_PROPS}/JSON?MaxRecords=12"
                );
                let body = http_get(&url, "application/json", 14)?;
                parse_pubchem_properties_json(&body)?
            } else {
                // Exact name / synonym / CAS match first (query param: PubChem rejects %2F in path segments)
                let url = format!(
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/property/{PUBCHEM_PROPS}/JSON?name={enc}"
                );
                let mut res = match http_get(&url, "application/json", 12) {
                    Ok(body) => parse_pubchem_properties_json(&body).unwrap_or_default(),
                    Err(_) => Vec::new(),
                };
                // Also try PubChem autocomplete if exact match returned empty (partial names like "ibuprof")
                if res.is_empty() {
                    let ac_url = format!(
                        "https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{enc}/json?limit=6"
                    );
                    if let Ok(ac_body) = http_get(&ac_url, "application/json", 8) {
                        let terms = parse_pubchem_autocomplete_json(&ac_body);
                        let mut seen_cids = HashSet::new();
                        for term in terms.into_iter().take(5) {
                            let t_url = format!(
                                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/property/{PUBCHEM_PROPS}/JSON?name={}",
                                url_encode(&term)
                            );
                            if let Ok(tb) = http_get(&t_url, "application/json", 8) {
                                if let Ok(mut th) = parse_pubchem_properties_json(&tb) {
                                    for h in th.drain(..) {
                                        if h.cid.map(|c| seen_cids.insert(c)).unwrap_or(true) {
                                            res.push(h);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                res
            }
        }
        DbSearchMode::ExactSmiles => {
            let url = format!(
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/property/{PUBCHEM_PROPS}/JSON?smiles={enc}"
            );
            let body = http_get(&url, "application/json", 12)?;
            parse_pubchem_properties_json(&body)?
        }
        DbSearchMode::Similarity2D => {
            let thr = sim_threshold.clamp(60, 99);
            let url = format!(
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/property/{PUBCHEM_PROPS}/JSON?smiles={enc}&Threshold={thr}&MaxRecords=15"
            );
            let body = http_get(&url, "application/json", 15)?;
            parse_pubchem_properties_json(&body)?
        }
        DbSearchMode::Substructure => {
            let url = format!(
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smiles/property/{PUBCHEM_PROPS}/JSON?smiles={enc}&MaxRecords=15"
            );
            let body = http_get(&url, "application/json", 15)?;
            parse_pubchem_properties_json(&body)?
        }
    };

    // Eagerly fetch curated 2D SDF for the top 3 PubChem hits so their previews & stereo are crisp immediately!
    for h in hits.iter_mut().take(3) {
        if let Some(cid) = h.cid {
            if let Ok(sdf) = fetch_pubchem_sdf(cid) {
                h.molfile = Some(sdf);
                if let Ok(m) = hit_to_molecule(h, 14.4, 0.0, 0.0) {
                    h.preview_mol = Some(m);
                }
            }
        }
    }
    Ok(hits)
}

fn query_chembl(query: &str, mode: DbSearchMode, sim_threshold: u32) -> Result<Vec<DbHit>, String> {
    let q = query.trim();
    if q.is_empty() {
        return Ok(Vec::new());
    }
    let enc = url_encode(q);
    let url = match mode {
        DbSearchMode::NameOrId => {
            let up = q.to_ascii_uppercase();
            if up.starts_with("CHEMBL") && up[6..].chars().all(|c| c.is_ascii_digit()) && up.len() > 6 {
                format!("https://www.ebi.ac.uk/chembl/api/data/molecule/{up}.json")
            } else if is_inchikey(q) {
                format!(
                    "https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__standard_inchi_key={enc}&limit=10"
                )
            } else {
                format!(
                    "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q={enc}&limit=12"
                )
            }
        }
        DbSearchMode::ExactSmiles => {
            format!(
                "https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__canonical_smiles__flexmatch={enc}&limit=10"
            )
        }
        DbSearchMode::Similarity2D => {
            let thr = sim_threshold.clamp(40, 100);
            format!(
                "https://www.ebi.ac.uk/chembl/api/data/similarity/{enc}/{thr}.json?limit=15"
            )
        }
        DbSearchMode::Substructure => {
            format!(
                "https://www.ebi.ac.uk/chembl/api/data/substructure/{enc}.json?limit=15"
            )
        }
    };
    let body = http_get(&url, "application/json", 15)?;
    let mut hits = parse_chembl_molecules_json(&body)?;
    if mode == DbSearchMode::NameOrId {
        rank_chembl_hits(&mut hits, q);
    }
    Ok(hits)
}

/// Rank ChEMBL text-search hits so exact preferred-name matches and approved single-component drugs
/// appear above unnamed analogs and multi-component salts.
pub fn rank_chembl_hits(hits: &mut [DbHit], query: &str) {
    let ql = query.trim().to_ascii_lowercase();
    let score = |h: &DbHit| -> (i32, i32, i32, i32) {
        let nl = h.name.trim().to_ascii_lowercase();
        let exact = if nl == ql { 2 } else if nl.starts_with(&ql) { 1 } else { 0 };
        let has_pref = if h.name != h.id { 1 } else { 0 };
        let phase = match h.max_phase.as_deref() {
            Some("Approved (Phase IV)") => 4,
            Some("Phase III") => 3,
            Some("Phase II") => 2,
            Some("Phase I") => 1,
            _ => 0,
        };
        let non_salt = if h.smiles.contains('.') { 0 } else { 1 };
        (exact, phase, has_pref, non_salt)
    };
    hits.sort_by(|a, b| score(b).cmp(&score(a)));
}

// ============================================================================
// 4. Async State & Background Worker Messages
// ============================================================================

#[derive(Debug)]
enum DbWorkerMsg {
    SearchFinished {
        req_id: u64,
        hits: Vec<DbHit>,
        warnings: Vec<String>,
    },
    PubChemSdfReady {
        cid: u64,
        molfile: Result<String, String>,
    },
}

pub struct ChemDbState {
    pub query: String,
    pub source: DbSourceFilter,
    pub mode: DbSearchMode,
    pub sim_threshold: u32,
    pub hits: Vec<DbHit>,
    pub selected_idx: Option<usize>,
    pub searching: bool,
    pub status: String,
    req_counter: u64,
    tx: Sender<DbWorkerMsg>,
    rx: Receiver<DbWorkerMsg>,
}

impl Default for ChemDbState {
    fn default() -> Self {
        let (tx, rx) = mpsc::channel();
        Self {
            query: String::from("aspirin"),
            source: DbSourceFilter::Both,
            mode: DbSearchMode::NameOrId,
            sim_threshold: 80,
            hits: Vec::new(),
            selected_idx: None,
            searching: false,
            status: "Enter a compound name, CAS number, CID, ChEMBL ID, formula, or SMILES and press Search.".to_string(),
            req_counter: 0,
            tx,
            rx,
        }
    }
}

impl ChemDbState {
    /// Launch a non-blocking search on a background thread.
    pub fn start_search(&mut self, ctx: &egui::Context) {
        let q = self.query.trim().to_string();
        if q.is_empty() {
            self.status = "Enter a search query first.".to_string();
            return;
        }
        self.req_counter += 1;
        let req_id = self.req_counter;
        self.searching = true;
        self.status = format!(
            "Searching {} ({}) for \"{}\"...",
            match self.source {
                DbSourceFilter::Both => "PubChem + ChEMBL",
                DbSourceFilter::PubChem => "PubChem",
                DbSourceFilter::ChEMBL => "ChEMBL",
            },
            self.mode.label(),
            q
        );

        let source = self.source;
        let mode = self.mode;
        let sim = self.sim_threshold;
        let tx = self.tx.clone();
        let repaint_ctx = ctx.clone();

        thread::spawn(move || {
            let mut hits = Vec::new();
            let mut warnings = Vec::new();

            match source {
                DbSourceFilter::PubChem => match query_pubchem(&q, mode, sim) {
                    Ok(mut h) => hits.append(&mut h),
                    Err(e) => warnings.push(e),
                },
                DbSourceFilter::ChEMBL => match query_chembl(&q, mode, sim) {
                    Ok(mut h) => hits.append(&mut h),
                    Err(e) => warnings.push(e),
                },
                DbSourceFilter::Both => {
                    // Run PubChem and ChEMBL in parallel for snappy combined results!
                    let q1 = q.clone();
                    let q2 = q.clone();
                    let h1 = thread::spawn(move || query_pubchem(&q1, mode, sim));
                    let h2 = thread::spawn(move || query_chembl(&q2, mode, sim));
                    match h1.join() {
                        Ok(Ok(mut ph)) => hits.append(&mut ph),
                        Ok(Err(e)) => warnings.push(format!("PubChem: {e}")),
                        Err(_) => warnings.push("PubChem worker panicked".to_string()),
                    }
                    match h2.join() {
                        Ok(Ok(mut ch)) => hits.append(&mut ch),
                        Ok(Err(e)) => warnings.push(format!("ChEMBL: {e}")),
                        Err(_) => warnings.push("ChEMBL worker panicked".to_string()),
                    }
                }
            }

            let _ = tx.send(DbWorkerMsg::SearchFinished {
                req_id,
                hits,
                warnings,
            });
            repaint_ctx.request_repaint();
        });
    }

    /// Fetch curated 2D SDF on demand for a selected PubChem hit that didn't get eagerly loaded.
    pub fn ensure_pubchem_sdf(&mut self, hit_idx: usize, ctx: &egui::Context) {
        let Some(hit) = self.hits.get_mut(hit_idx) else {
            return;
        };
        if hit.source != DbSource::PubChem || hit.molfile.is_some() || hit.sdf_loading {
            return;
        }
        let Some(cid) = hit.cid else {
            return;
        };
        hit.sdf_loading = true;
        let tx = self.tx.clone();
        let repaint_ctx = ctx.clone();
        thread::spawn(move || {
            let res = fetch_pubchem_sdf(cid);
            let _ = tx.send(DbWorkerMsg::PubChemSdfReady { cid, molfile: res });
            repaint_ctx.request_repaint();
        });
    }

    /// Poll background worker responses (call once per UI frame).
    pub fn poll(&mut self) {
        while let Ok(msg) = self.rx.try_recv() {
            match msg {
                DbWorkerMsg::SearchFinished {
                    req_id,
                    hits,
                    warnings,
                } => {
                    if req_id != self.req_counter {
                        continue; // stale search superseded by a newer one
                    }
                    self.searching = false;
                    let n = hits.len();
                    self.hits = hits;
                    self.selected_idx = if n > 0 { Some(0) } else { None };
                    if n > 0 {
                        let pc = self.hits.iter().filter(|h| h.source == DbSource::PubChem).count();
                        let ch = self.hits.iter().filter(|h| h.source == DbSource::ChEMBL).count();
                        self.status = if warnings.is_empty() {
                            format!("Found {n} compound(s) (PubChem: {pc}, ChEMBL: {ch})")
                        } else {
                            format!("Found {n} compound(s) • Note: {}", warnings.join("; "))
                        };
                    } else if !warnings.is_empty() {
                        self.status = format!("No matches found ({})", warnings.join("; "));
                    } else {
                        self.status = "No matching compounds found.".to_string();
                    }
                }
                DbWorkerMsg::PubChemSdfReady { cid, molfile } => {
                    for h in &mut self.hits {
                        if h.cid == Some(cid) {
                            h.sdf_loading = false;
                            if let Ok(ref sdf) = molfile {
                                h.molfile = Some(sdf.clone());
                                if let Ok(m) = hit_to_molecule(h, 14.4, 0.0, 0.0) {
                                    h.preview_mol = Some(m);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// ============================================================================
// 5. 2D Vector Structure Preview Renderer
// ============================================================================

fn paint_molecule_preview(painter: &egui::Painter, rect: Rect, mol: &ChemMolecule) {
    painter.rect_filled(rect, 2.0, Color32::WHITE);
    painter.rect_stroke(
        rect,
        2.0,
        Stroke::new(1.0_f32, Color32::from_gray(150)),
        egui::StrokeKind::Inside,
    );
    if mol.atoms.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "No 2D coordinates",
            FontId::proportional(11.0),
            Color32::GRAY,
        );
        return;
    }

    let (mut x1, mut y1, mut x2, mut y2) = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for a in &mol.atoms {
        x1 = x1.min(a.x);
        y1 = y1.min(a.y);
        x2 = x2.max(a.x);
        y2 = y2.max(a.y);
    }
    let w = (x2 - x1).max(12.0);
    let h = (y2 - y1).max(12.0);
    let pad = 20.0_f32;
    let avail_w = (rect.width() - pad * 2.0).max(20.0) as f64;
    let avail_h = (rect.height() - pad * 2.0).max(20.0) as f64;
    let scale = (avail_w / w).min(avail_h / h).min(2.6);
    let cx = 0.5 * (x1 + x2);
    let cy = 0.5 * (y1 + y2);
    let to_pos = |ax: f64, ay: f64| -> Pos2 {
        Pos2::new(
            rect.center().x + ((ax - cx) * scale) as f32,
            rect.center().y + ((ay - cy) * scale) as f32,
        )
    };

    let stroke_w = (1.35_f32).max((scale * 0.07) as f32).min(2.0);

    // Draw bonds
    for b in &mol.bonds {
        let (Some(a1), Some(a2)) = (mol.atoms.get(b.begin_atom), mol.atoms.get(b.end_atom)) else {
            continue;
        };
        let p1 = to_pos(a1.x, a1.y);
        let p2 = to_pos(a2.x, a2.y);
        let dv = p2 - p1;
        let len = dv.length();
        if len < 1e-3 {
            continue;
        }
        let dir = dv / len;
        let perp = Vec2::new(-dir.y, dir.x);
        let col = Color32::from_gray(25);

        match b.display {
            BondDisplay::WedgeBegin | BondDisplay::WedgeEnd => {
                let (tip, base) = if b.display == BondDisplay::WedgeBegin {
                    (p1, p2)
                } else {
                    (p2, p1)
                };
                let hw = (len * 0.14).clamp(2.2, 4.2);
                painter.add(egui::Shape::convex_polygon(
                    vec![tip, base + perp * hw, base - perp * hw],
                    col,
                    Stroke::NONE,
                ));
            }
            BondDisplay::WedgedHashBegin | BondDisplay::WedgedHashEnd => {
                let (tip, base) = if b.display == BondDisplay::WedgedHashBegin {
                    (p1, p2)
                } else {
                    (p2, p1)
                };
                let hw = (len * 0.14).clamp(2.2, 4.2);
                let steps = 5;
                for k in 1..=steps {
                    let t = (k as f32) / (steps as f32);
                    let mid = tip + (base - tip) * t;
                    let w_k = hw * t;
                    painter.line_segment(
                        [mid - perp * w_k, mid + perp * w_k],
                        Stroke::new(stroke_w * 0.85, col),
                    );
                }
            }
            _ => match b.order {
                2 => {
                    let off = (len * 0.09).clamp(1.6, 3.0);
                    painter.line_segment([p1 + perp * off, p2 + perp * off], Stroke::new(stroke_w, col));
                    painter.line_segment([p1 - perp * off, p2 - perp * off], Stroke::new(stroke_w, col));
                }
                3 => {
                    let off = (len * 0.11).clamp(2.0, 3.5);
                    painter.line_segment([p1, p2], Stroke::new(stroke_w, col));
                    painter.line_segment([p1 + perp * off, p2 + perp * off], Stroke::new(stroke_w, col));
                    painter.line_segment([p1 - perp * off, p2 - perp * off], Stroke::new(stroke_w, col));
                }
                6 => {
                    let off = (len * 0.09).clamp(1.6, 2.8);
                    painter.line_segment([p1, p2], Stroke::new(stroke_w, col));
                    painter.line_segment(
                        [p1 + dir * (len * 0.15) + perp * off, p2 - dir * (len * 0.15) + perp * off],
                        Stroke::new(stroke_w * 0.8, Color32::from_gray(90)),
                    );
                }
                _ => {
                    painter.line_segment([p1, p2], Stroke::new(stroke_w, col));
                }
            },
        }
    }

    // Draw heteroatom / charged atom labels with clean white knockout
    let font_sz = (scale * 0.65).clamp(8.5, 13.0) as f32;
    let font = FontId::proportional(font_sz);
    let val = CcVal::compute(mol);

    for (ai, a) in mol.atoms.iter().enumerate() {
        let deg = mol
            .bonds
            .iter()
            .filter(|b| b.begin_atom == ai || b.end_atom == ai)
            .count();
        let show = a.z != 6 || deg == 0 || a.charge != 0 || a.isotope != 0 || a.label.is_some();
        if !show {
            continue;
        }
        let sym = a
            .label
            .clone()
            .unwrap_or_else(|| get_element(a.z).map(|e| e.symbol.to_string()).unwrap_or_else(|| "?".to_string()));
        let h_cnt = val.implicit_h(ai);
        let mut lbl = sym;
        if a.z != 6 && h_cnt > 0 {
            lbl.push('H');
            if h_cnt > 1 {
                for ch in h_cnt.to_string().chars() {
                    lbl.push(match ch {
                        '0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
                        '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉',
                        _ => ch,
                    });
                }
            }
        }
        if a.charge != 0 {
            let abs_c = a.charge.unsigned_abs();
            if abs_c > 1 {
                for ch in abs_c.to_string().chars() {
                    lbl.push(match ch {
                        '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
                        '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
                        _ => ch,
                    });
                }
            }
            lbl.push(if a.charge > 0 { '⁺' } else { '⁻' });
        }
        let pos = to_pos(a.x, a.y);
        let bg_w = (lbl.chars().count() as f32) * (font_sz * 0.52) + 3.0;
        let bg_h = font_sz + 1.5;
        painter.rect_filled(
            Rect::from_center_size(pos, Vec2::new(bg_w, bg_h)),
            2.0,
            Color32::WHITE,
        );
        painter.text(
            pos,
            egui::Align2::CENTER_CENTER,
            lbl,
            font.clone(),
            cpk_element_color(a.z),
        );
    }
}

// ============================================================================
// 6. Floating Interactive Search Panel (`Search > Search PubChem / ChEMBL...`)
// ============================================================================

#[derive(Debug, Clone)]
pub enum ChemDbAction {
    None,
    UseCanvasSelection,
    ImportHit { hit: DbHit, with_caption: bool },
    ImportBatch(Vec<DbHit>),
    CopyText { label: String, text: String },
}

pub fn show_chem_db_window(
    ctx: &egui::Context,
    open: &mut bool,
    state: &mut ChemDbState,
    has_canvas_mol: bool,
) -> ChemDbAction {
    let mut act = ChemDbAction::None;
    let mut trigger_search = false;
    let mut need_sdf_for: Option<usize> = None;

    egui::Window::new("Search PubChem / ChEMBL")
        .open(open)
        .default_width(740.0)
        .default_height(470.0)
        .resizable(true)
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            // Row 1: Source + Mode + Similarity Slider
            ui.horizontal_wrapped(|ui| {
                ui.strong("Database:");
                ui.selectable_value(&mut state.source, DbSourceFilter::Both, "Both");
                ui.selectable_value(&mut state.source, DbSourceFilter::PubChem, "PubChem");
                ui.selectable_value(&mut state.source, DbSourceFilter::ChEMBL, "ChEMBL");
                ui.separator();
                ui.strong("Mode:");
                for m in [
                    DbSearchMode::NameOrId,
                    DbSearchMode::ExactSmiles,
                    DbSearchMode::Similarity2D,
                    DbSearchMode::Substructure,
                ] {
                    ui.selectable_value(&mut state.mode, m, m.label());
                }
            });

            // Row 2: Search input + buttons
            ui.horizontal(|ui| {
                ui.label("Query:");
                let resp = ui.add_sized(
                    [310.0, 22.0],
                    egui::TextEdit::singleline(&mut state.query)
                        .hint_text("e.g. aspirin, ibuprofen, 50-78-2, CID 2244, CHEMBL25, C9H8O4"),
                );
                if resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                    trigger_search = true;
                }
                if ui
                    .add_enabled(!state.searching, egui::Button::new(if state.searching { "Searching..." } else { "Search" }))
                    .clicked()
                {
                    trigger_search = true;
                }
                if ui
                    .add_enabled(has_canvas_mol && !state.searching, egui::Button::new("From Canvas"))
                    .on_hover_text("Load SMILES from the active/selected structure on the canvas and search")
                    .clicked()
                {
                    act = ChemDbAction::UseCanvasSelection;
                }
                if state.mode == DbSearchMode::Similarity2D {
                    ui.add(
                        egui::Slider::new(&mut state.sim_threshold, 60..=99)
                            .suffix("%")
                            .text("Tanimoto"),
                    );
                }
            });

            // Status bar
            ui.horizontal(|ui| {
                if state.searching {
                    ui.spinner();
                }
                ui.label(RichText::new(&state.status).small().color(Color32::from_gray(60)));
            });
            ui.separator();

            // Main Split Area: Left Results Table | Right 2D Preview + Metadata
            ui.columns(2, |cols| {
                // LEFT: Results List
                cols[0].strong(format!("Results ({})", state.hits.len()));
                egui::ScrollArea::vertical()
                    .id_salt("chem_db_hits_scroll")
                    .max_height(290.0)
                    .show(&mut cols[0], |ui| {
                        egui::Grid::new("chem_db_hits_grid")
                            .striped(true)
                            .num_columns(5)
                            .spacing([8.0, 3.0])
                            .show(ui, |ui| {
                                ui.strong("DB");
                                ui.strong("ID");
                                ui.strong("Name");
                                ui.strong("Formula");
                                ui.strong("MW");
                                ui.end_row();

                                for (idx, h) in state.hits.iter().enumerate() {
                                    let is_sel = state.selected_idx == Some(idx);
                                    let badge_col = match h.source {
                                        DbSource::PubChem => Color32::from_rgb(20, 90, 165),
                                        DbSource::ChEMBL => Color32::from_rgb(0, 120, 105),
                                    };
                                    ui.label(RichText::new(h.source.label()).small().color(badge_col));
                                    if ui.selectable_label(is_sel, &h.id).clicked() {
                                        state.selected_idx = Some(idx);
                                        if h.source == DbSource::PubChem && h.molfile.is_none() {
                                            need_sdf_for = Some(idx);
                                        }
                                    }
                                    let short_name: String = h.name.chars().take(24).collect();
                                    let name_resp = ui
                                        .selectable_label(is_sel, short_name)
                                        .on_hover_text(&h.name);
                                    if name_resp.clicked() {
                                        state.selected_idx = Some(idx);
                                        if h.source == DbSource::PubChem && h.molfile.is_none() {
                                            need_sdf_for = Some(idx);
                                        }
                                    }
                                    if name_resp.double_clicked() {
                                        act = ChemDbAction::ImportHit {
                                            hit: h.clone(),
                                            with_caption: false,
                                        };
                                    }
                                    ui.label(&h.formula);
                                    ui.label(&h.molecular_weight);
                                    ui.end_row();
                                }
                            });
                    });

                // RIGHT: 2D Vector Structure Preview + Metadata Card
                let ui = &mut cols[1];
                if let Some(sel) = state.selected_idx.and_then(|i| state.hits.get(i)) {
                    ui.horizontal(|ui| {
                        ui.strong(&sel.name);
                        ui.label(RichText::new(format!("({})", sel.id)).small().color(Color32::DARK_GRAY));
                    });
                    let (rect, _) = ui.allocate_exact_size(Vec2::new(330.0, 168.0), Sense::hover());
                    if let Some(ref pmol) = sel.preview_mol {
                        paint_molecule_preview(ui.painter(), rect, pmol);
                    } else {
                        ui.painter().rect_filled(rect, 2.0, Color32::WHITE);
                        ui.painter().rect_stroke(
                            rect,
                            2.0,
                            Stroke::new(1.0_f32, Color32::from_gray(150)),
                            egui::StrokeKind::Inside,
                        );
                        ui.painter().text(
                            rect.center(),
                            egui::Align2::CENTER_CENTER,
                            if sel.sdf_loading { "Fetching 2D SDF..." } else { "No structure preview" },
                            FontId::proportional(11.0),
                            Color32::GRAY,
                        );
                    }

                    ui.add_space(4.0);
                    egui::ScrollArea::vertical()
                        .id_salt("chem_db_meta_scroll")
                        .max_height(110.0)
                        .show(ui, |ui| {
                            ui.label(format!(
                                "Formula: {}   |   MW: {}   |   2D MOL: {}",
                                if sel.formula.is_empty() { "—" } else { &sel.formula },
                                if sel.molecular_weight.is_empty() { "—" } else { &sel.molecular_weight },
                                if sel.molfile.is_some() {
                                    "Curated 2D"
                                } else if sel.sdf_loading {
                                    "Loading..."
                                } else {
                                    "SMILES layout"
                                }
                            ));
                            if let Some(ref ph) = sel.max_phase {
                                let yr = sel
                                    .first_approval
                                    .map(|y| format!(" (First Approved: {y})"))
                                    .unwrap_or_default();
                                ui.colored_label(
                                    Color32::from_rgb(0, 110, 80),
                                    format!("Clinical Phase: {ph}{yr}"),
                                );
                            }
                            if sel.alogp.is_some() || sel.psa.is_some() || sel.ro5_violations.is_some() {
                                ui.label(format!(
                                    "ALogP: {}   |   TPSA: {} Å²   |   Ro5 Violations: {}",
                                    sel.alogp.as_deref().unwrap_or("—"),
                                    sel.psa.as_deref().unwrap_or("—"),
                                    sel.ro5_violations
                                        .map(|v| v.to_string())
                                        .unwrap_or_else(|| "—".to_string())
                                ));
                            }
                            if let Some(ref sim) = sel.similarity_pct {
                                ui.label(format!("Similarity: {sim}"));
                            }
                            if let Some(ref iupac) = sel.iupac_name {
                                ui.label(RichText::new(format!("IUPAC: {iupac}")).small());
                            }
                            if let Some(ref ik) = sel.inchikey {
                                ui.label(RichText::new(format!("InChIKey: {ik}")).small().monospace());
                            }
                            if !sel.smiles.is_empty() {
                                ui.label(RichText::new(format!("SMILES: {}", sel.smiles)).small().monospace());
                            }
                        });
                } else {
                    ui.centered_and_justified(|ui| {
                        ui.label(
                            RichText::new("Select a search hit to preview its 2D structure and properties.\nDouble-click any row to import directly onto the canvas.")
                                .color(Color32::DARK_GRAY),
                        );
                    });
                }
            });

            ui.separator();
            // Bottom Action Bar
            let sel_hit = state.selected_idx.and_then(|i| state.hits.get(i)).cloned();
            ui.horizontal_wrapped(|ui| {
                if ui
                    .add_enabled(sel_hit.is_some(), egui::Button::new("Import Structure"))
                    .clicked()
                {
                    if let Some(h) = sel_hit.clone() {
                        act = ChemDbAction::ImportHit {
                            hit: h,
                            with_caption: false,
                        };
                    }
                }
                if ui
                    .add_enabled(sel_hit.is_some(), egui::Button::new("Import + Label"))
                    .on_hover_text("Import 2D structure and place Name, ID, Formula, MW caption underneath")
                    .clicked()
                {
                    if let Some(h) = sel_hit.clone() {
                        act = ChemDbAction::ImportHit {
                            hit: h,
                            with_caption: true,
                        };
                    }
                }
                ui.separator();
                if ui
                    .add_enabled(
                        sel_hit.as_ref().map(|h| !h.smiles.is_empty()).unwrap_or(false),
                        egui::Button::new("Copy SMILES"),
                    )
                    .clicked()
                {
                    if let Some(h) = sel_hit.as_ref() {
                        act = ChemDbAction::CopyText {
                            label: "SMILES".to_string(),
                            text: h.smiles.clone(),
                        };
                    }
                }
                if ui
                    .add_enabled(
                        sel_hit.as_ref().and_then(|h| h.inchikey.as_ref()).is_some(),
                        egui::Button::new("Copy InChIKey"),
                    )
                    .clicked()
                {
                    if let Some(ik) = sel_hit.as_ref().and_then(|h| h.inchikey.clone()) {
                        act = ChemDbAction::CopyText {
                            label: "InChIKey".to_string(),
                            text: ik,
                        };
                    }
                }
                ui.separator();
                let batch_ok = state.hits.len() > 1;
                if ui
                    .add_enabled(batch_ok, egui::Button::new(format!("Import Grid (Top {})", state.hits.len().min(6))))
                    .on_hover_text("Import up to 6 search hits arranged neatly in a grid with labels")
                    .clicked()
                {
                    act = ChemDbAction::ImportBatch(state.hits.iter().take(6).cloned().collect());
                }
            });
        });

    if trigger_search {
        state.start_search(ctx);
    }
    if let Some(idx) = need_sdf_for {
        state.ensure_pubchem_sdf(idx, ctx);
    }
    act
}

// ============================================================================
// 6b. Non-Blocking Live PubChem IUPAC + Common Name Resolver (Structure → Name)
// ============================================================================

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ResolvedName {
    pub pubchem_title: Option<String>,
    pub iupac_name: Option<String>,
    pub cid: Option<u64>,
    pub pending: bool,
    pub failed: bool,
}

impl ResolvedName {
    /// Formats an authoritative display name combining IUPAC name & common/drug name when both exist.
    pub fn format_display(&self) -> Option<String> {
        let iupac = self
            .iupac_name
            .as_deref()
            .map(str::trim)
            .filter(|s| !s.is_empty());
        let title = self
            .pubchem_title
            .as_deref()
            .map(str::trim)
            .filter(|s| !s.is_empty());
        match (iupac, title) {
            (Some(iu), Some(ti)) if !iu.eq_ignore_ascii_case(ti) => Some(format!("{iu} ({ti})")),
            (Some(iu), _) => Some(iu.to_string()),
            (None, Some(ti)) => Some(ti.to_string()),
            (None, None) => None,
        }
    }
}

#[derive(Debug)]
struct NameWorkerReply {
    smiles: String,
    caption_id: Option<u32>,
    resolved: ResolvedName,
}

pub struct LiveNameResolver {
    cache: HashMap<String, ResolvedName>,
    tx: Sender<NameWorkerReply>,
    rx: Receiver<NameWorkerReply>,
}

impl Default for LiveNameResolver {
    fn default() -> Self {
        let (tx, rx) = mpsc::channel();
        Self {
            cache: HashMap::new(),
            tx,
            rx,
        }
    }
}

impl LiveNameResolver {
    /// Returns cached PubChem name (if resolved) for a given SMILES string.
    pub fn get_cached(&self, smiles: &str) -> Option<&ResolvedName> {
        self.cache.get(smiles.trim())
    }

    pub fn is_pending(&self, smiles: &str) -> bool {
        self.cache
            .get(smiles.trim())
            .map(|r| r.pending)
            .unwrap_or(false)
    }

    pub fn has_any_pending(&self) -> bool {
        self.cache.values().any(|r| r.pending)
    }

    /// Triggers a non-blocking PubChem `Title + IUPACName` lookup by SMILES (deduplicated via cache).
    /// Returns `Some(display_name)` immediately if already cached!
    pub fn resolve(
        &mut self,
        smiles: &str,
        caption_id: Option<u32>,
        ctx: &egui::Context,
    ) -> Option<String> {
        let key = smiles.trim().to_string();
        if key.is_empty() {
            return None;
        }
        if let Some(existing) = self.cache.get(&key) {
            if let Some(disp) = existing.format_display() {
                return Some(disp);
            }
            if existing.pending {
                return None;
            }
        }
        self.cache.insert(
            key.clone(),
            ResolvedName {
                pending: true,
                ..Default::default()
            },
        );
        let tx = self.tx.clone();
        let repaint_ctx = ctx.clone();
        thread::spawn(move || {
            let enc = url_encode(&key);
            let url = format!(
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/property/Title,IUPACName/JSON?smiles={enc}"
            );
            let resolved = match http_get(&url, "application/json", 8)
                .and_then(|body| parse_pubchem_properties_json(&body))
            {
                Ok(hits) if !hits.is_empty() => {
                    let h = &hits[0];
                    ResolvedName {
                        pubchem_title: if !h.name.is_empty() && !h.name.starts_with("CID ") {
                            Some(h.name.clone())
                        } else {
                            None
                        },
                        iupac_name: h.iupac_name.clone(),
                        cid: h.cid,
                        pending: false,
                        failed: false,
                    }
                }
                _ => ResolvedName {
                    pending: false,
                    failed: true,
                    ..Default::default()
                },
            };
            let _ = tx.send(NameWorkerReply {
                smiles: key,
                caption_id,
                resolved,
            });
            repaint_ctx.request_repaint();
        });
        None
    }

    /// Poll finished lookups, update the SMILES cache, and return `(caption_id, smiles, resolved)` events.
    pub fn poll(&mut self) -> Vec<(Option<u32>, String, ResolvedName)> {
        let mut out = Vec::new();
        while let Ok(rep) = self.rx.try_recv() {
            self.cache.insert(rep.smiles.clone(), rep.resolved.clone());
            out.push((rep.caption_id, rep.smiles, rep.resolved));
        }
        out
    }
}

// ============================================================================
// 7. Unit Tests (JSON, PubChem/ChEMBL Parsing, H-Pruning, Bond Normalization)
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use chiralcraft_chem::ChemAtom;

    #[test]
    fn test_json_parser_escapes_and_nested() {
        let s = r#"{"a": 123, "b": "line1\nline2\t\"q\"\u00B5", "c": [true, false, null, -4.5]}"#;
        let v = JsonValue::parse(s).expect("json parse");
        assert_eq!(v.get("a").and_then(|x| x.as_u64()), Some(123));
        assert_eq!(v.get("b").and_then(|x| x.as_str()), Some("line1\nline2\t\"q\"µ"));
        let arr = v.get("c").and_then(|x| x.as_array()).unwrap();
        assert_eq!(arr.len(), 4);
        assert_eq!(arr[0], JsonValue::Bool(true));
        assert_eq!(arr[2], JsonValue::Null);
    }

    #[test]
    fn test_parse_pubchem_properties_json() {
        let sample = r#"{
          "PropertyTable": {
            "Properties": [
              {
                "CID": 2244,
                "MolecularFormula": "C9H8O4",
                "MolecularWeight": "180.16",
                "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                "IUPACName": "2-acetyloxybenzoic acid",
                "Title": "Aspirin"
              }
            ]
          }
        }"#;
        let hits = parse_pubchem_properties_json(sample).expect("parse pubchem");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].source, DbSource::PubChem);
        assert_eq!(hits[0].id, "CID 2244");
        assert_eq!(hits[0].name, "Aspirin");
        assert_eq!(hits[0].formula, "C9H8O4");
        assert_eq!(hits[0].molecular_weight, "180.16");
        assert!(hits[0].preview_mol.is_some());
    }

    #[test]
    fn test_parse_chembl_molecules_json() {
        let sample = r#"{
          "molecules": [
            {
              "molecule_chembl_id": "CHEMBL25",
              "pref_name": "ASPIRIN",
              "max_phase": "4.0",
              "first_approval": 1950,
              "molecule_properties": {
                "full_molformula": "C9H8O4",
                "full_mwt": "180.16",
                "alogp": "1.31",
                "psa": "63.60",
                "num_ro5_violations": 0
              },
              "molecule_structures": {
                "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "standard_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                "molfile": "\n  RDKit  2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n"
              }
            }
          ]
        }"#;
        let hits = parse_chembl_molecules_json(sample).expect("parse chembl");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].source, DbSource::ChEMBL);
        assert_eq!(hits[0].id, "CHEMBL25");
        assert_eq!(hits[0].name, "ASPIRIN");
        assert_eq!(hits[0].max_phase.as_deref(), Some("Approved (Phase IV)"));
        assert_eq!(hits[0].first_approval, Some(1950));
        assert_eq!(hits[0].ro5_violations, Some(0));
        assert!(hits[0].molfile.is_some());
        assert_eq!(hits[0].preview_mol.as_ref().unwrap().atoms.len(), 2);
    }

    #[test]
    fn test_strip_nonessential_hydrogens_preserves_stereo_and_isotopes() {
        let mut mol = ChemMolecule::new();
        let mut c = ChemAtom::new(1, 6);
        c.x = 0.0;
        c.y = 0.0;
        mol.add_atom(c);
        // Normal explicit H (should be stripped)
        let mut h_plain = ChemAtom::new(2, 1);
        h_plain.x = 1.0;
        mol.add_atom(h_plain);
        // Wedged chiral H (MUST be preserved)
        let mut h_stereo = ChemAtom::new(3, 1);
        h_stereo.x = -1.0;
        mol.add_atom(h_stereo);
        // Deuterium (isotope 2, MUST be preserved)
        let mut h_deut = ChemAtom::new(4, 1);
        h_deut.x = 0.0;
        h_deut.y = 1.0;
        h_deut.isotope = 2;
        mol.add_atom(h_deut);

        mol.add_bond(10, 0, 1, 1, BondDisplay::Solid);
        mol.add_bond(11, 0, 2, 1, BondDisplay::WedgeBegin);
        mol.add_bond(12, 0, 3, 1, BondDisplay::Solid);

        strip_nonessential_explicit_hydrogens(&mut mol);
        assert_eq!(mol.atoms.len(), 3, "Only plain H should be stripped");
        assert_eq!(mol.bonds.len(), 2);
        assert!(mol.bonds.iter().any(|b| b.display == BondDisplay::WedgeBegin));
        assert!(mol.atoms.iter().any(|a| a.isotope == 2));
    }

    #[test]
    fn test_normalize_molecule_to_stylesheet_scales_bond_length() {
        let mut mol = ChemMolecule::new();
        let mut a1 = ChemAtom::new(1, 6);
        a1.x = 0.0;
        a1.y = 0.0;
        let mut a2 = ChemAtom::new(2, 6);
        a2.x = 0.82; // RDKit 0.82 Å bond
        a2.y = 0.0;
        mol.add_atom(a1);
        mol.add_atom(a2);
        mol.add_bond(1, 0, 1, 1, BondDisplay::Solid);

        normalize_molecule_to_stylesheet(&mut mol, 14.4, 200.0, 150.0);
        let d = (mol.atoms[1].x - mol.atoms[0].x).hypot(mol.atoms[1].y - mol.atoms[0].y);
        assert!((d - 14.4).abs() < 1e-6, "Expected 14.4 pt bond, got {d}");
        let cx = 0.5 * (mol.atoms[0].x + mol.atoms[1].x);
        assert!((cx - 200.0).abs() < 1e-6);
    }

    #[test]
    fn test_rank_chembl_hits_prefers_exact_approved_drug_over_unnamed_salt() {
        let salt = DbHit {
            source: DbSource::ChEMBL,
            id: "CHEMBL5665374".to_string(),
            cid: None,
            name: "CHEMBL5665374".to_string(),
            iupac_name: None,
            formula: "C17H20ClNO3".to_string(),
            molecular_weight: "321.80".to_string(),
            smiles: "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5.Cl".to_string(),
            inchikey: None,
            max_phase: None,
            first_approval: None,
            alogp: None,
            psa: None,
            ro5_violations: None,
            similarity_pct: None,
            molfile: None,
            preview_mol: None,
            sdf_loading: false,
        };
        let morphine = DbHit {
            source: DbSource::ChEMBL,
            id: "CHEMBL70".to_string(),
            cid: None,
            name: "MORPHINE".to_string(),
            iupac_name: None,
            formula: "C17H19NO3".to_string(),
            molecular_weight: "285.34".to_string(),
            smiles: "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5".to_string(),
            inchikey: None,
            max_phase: Some("Approved (Phase IV)".to_string()),
            first_approval: Some(1984),
            alogp: None,
            psa: None,
            ro5_violations: None,
            similarity_pct: None,
            molfile: None,
            preview_mol: None,
            sdf_loading: false,
        };
        let mut hits = vec![salt, morphine];
        rank_chembl_hits(&mut hits, "morphine");
        assert_eq!(hits[0].id, "CHEMBL70");
        assert_eq!(hits[0].name, "MORPHINE");
    }

    #[test]
    fn test_pubchem_explicit_h_sdf_import_to_clean_acs_skeletal_molecule() {
        // Real PubChem 2D SDF format with explicit hydrogens (ethanol: 3 heavy + 6 H = 9 atoms)
        let sdf = "\
702
  -OEChem-11262600572D

  9  8  0     0  0  0  0  0  0999 V2000
    2.5369   -0.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2690   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.4030    0.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.8030    0.7200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.8059    0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2690   -0.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.8690   -0.7200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  4  1  0  0  0  0
  2  3  1  0  0  0  0
  2  5  1  0  0  0  0
  2  6  1  0  0  0  0
  3  7  1  0  0  0  0
  3  8  1  0  0  0  0
  3  9  1  0  0  0  0
M  END
$$$$
";
        let hit = DbHit {
            source: DbSource::PubChem,
            id: "CID 702".to_string(),
            cid: Some(702),
            name: "Ethanol".to_string(),
            iupac_name: Some("ethanol".to_string()),
            formula: "C2H6O".to_string(),
            molecular_weight: "46.07".to_string(),
            smiles: "CCO".to_string(),
            inchikey: None,
            max_phase: None,
            first_approval: None,
            alogp: None,
            psa: None,
            ro5_violations: None,
            similarity_pct: None,
            molfile: Some(sdf.to_string()),
            preview_mol: None,
            sdf_loading: false,
        };
        let mol = hit_to_molecule(&hit, 14.4, 200.0, 160.0).expect("import ethanol");
        assert_eq!(mol.atoms.len(), 3, "6 plain explicit H atoms must be pruned");
        assert_eq!(mol.bonds.len(), 2);
        let bl = (mol.atoms[1].x - mol.atoms[0].x).hypot(mol.atoms[1].y - mol.atoms[0].y);
        assert!((bl - 14.4).abs() < 0.05, "Bond length should scale to 14.4 pt, got {bl}");
        let val = CcVal::compute(&mol);
        // O has 1 implicit H, CH2 has 2, CH3 has 3 -> total 6 implicit H
        let total_h: u32 = (0..mol.atoms.len()).map(|i| val.implicit_h(i) as u32).sum();
        assert_eq!(total_h, 6);
    }

    #[test]
    fn test_resolved_name_formatting() {
        let r = ResolvedName {
            pubchem_title: Some("Aspirin".to_string()),
            iupac_name: Some("2-acetyloxybenzoic acid".to_string()),
            cid: Some(2244),
            pending: false,
            failed: false,
        };
        assert_eq!(
            r.format_display().as_deref(),
            Some("2-acetyloxybenzoic acid (Aspirin)")
        );
        let same = ResolvedName {
            pubchem_title: Some("Ethanol".to_string()),
            iupac_name: Some("ethanol".to_string()),
            cid: Some(702),
            pending: false,
            failed: false,
        };
        assert_eq!(same.format_display().as_deref(), Some("ethanol"));
    }

    #[test]
    fn test_json_unicode_escape_edge_cases() {
        // Valid surrogate pair -> U+1F600
        let v = JsonValue::parse(r#""\uD83D\uDE00""#).unwrap();
        assert_eq!(v.as_str(), Some("\u{1F600}"));
        // Lone surrogates become U+FFFD rather than vanishing
        let v = JsonValue::parse(r#""a\uD83Db""#).unwrap();
        assert_eq!(v.as_str(), Some("a\u{FFFD}b"));
        let v = JsonValue::parse(r#""\uDE00""#).unwrap();
        assert_eq!(v.as_str(), Some("\u{FFFD}"));
        // Non-hex digits (incl. '+') rejected
        assert!(JsonValue::parse(r#""\u+041""#).is_err());
    }

    #[test]
    fn test_chembl_unknown_phase() {
        assert_eq!(
            format_chembl_phase(Some(&JsonValue::Number(-1.0))).as_deref(),
            Some("Unknown Phase")
        );
        assert_eq!(
            format_chembl_phase(Some(&JsonValue::String("0.5".into()))).as_deref(),
            Some("Early Phase I")
        );
    }
}
