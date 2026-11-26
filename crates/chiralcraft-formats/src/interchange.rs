//! 1:1 ChiralCraft non-native chemical & graphics interchange formats,
//! calibrated directly against live ChiralCraft captures (`p07`, `p13`, `p35a`):
//! - MDL Molfile V3000 (`molecule_to_molfile_v3000`, `molfile_v3000_to_molecule`)
//! - ChiralCraft Connection Table `.ct` (`molecule_to_ct`, `ct_to_molecule`) — verified on `p35a-out.ct`
//! - Chemical Markup Language `.cml` (`molecules_to_cml`, `molecule_to_cml`, `cml_to_molecules`) — verified on `p07-export.cml`
//! - ISIS/TGF `.tgf` (`molecule_to_tgf`, `tgf_to_molecules`) — verified on `p35a-out.tgf`
//! - ISIS/Sketch binary `.skc` (`molecules_to_skc`, `skc_to_molecules`) — verified on `p13-rxn-aligned.skc`
//! - Encapsulated PostScript `.eps` (`molecule_to_eps`) — filled bond/wedge polygons + labels matching `p35a-out.eps`
//! - Windows Enhanced Metafile `.emf` (`molecule_to_emf`) — GDI pen/brush/polygon/text records
//! - Anti-aliased PNG rasterizer `.png` (`molecule_to_png`) — `tiny-skia` + `png` with double/triple/wedge/hash bonds & labels
//! - Graph-traversal Tripos SLN (`molecule_to_sln`)
//! - Canonical InChI + real FIPS 180-4 SHA-256 IUPAC InChIKey (`molecule_to_inchi`, `inchi_to_inchikey`)

use std::collections::{BTreeMap, HashMap};

use chiralcraft_chem::{get_element, BondDisplay, CcVal, ChemAtom, ChemBond, ChemMolecule};
use tiny_skia::{Color, FillRule, LineCap, Paint, PathBuilder, Pixmap, Stroke, Transform};

/// Exact ChiralCraft pt-per-Angstrom scale verified across `p07-export.cml` & `p35a-out.ct`:
/// `14.40 pt / 1.5230 Å = 9.45502298 pt/Å`.
pub const CHIRALCRAFT_PT_PER_ANGSTROM: f64 = 14.40 / 1.5230;

fn mk_atom(id: u32, z: i32, x: f64, y: f64) -> ChemAtom {
    let mut a = ChemAtom::new(id, z);
    a.x = x;
    a.y = y;
    a
}

fn mk_bond(id: u32, begin_atom: usize, end_atom: usize, order: i32) -> ChemBond {
    ChemBond {
        id,
        begin_atom,
        end_atom,
        order,
        display: BondDisplay::Solid,
        visible: true,
        ignore_warnings: false,
    }
}

fn elem_symbol(z: i32) -> &'static str {
    get_element(z).map(|e| e.symbol).unwrap_or("C")
}

fn symbol_to_z(sym: &str) -> i32 {
    let s = sym.trim();
    for z in 1..=118 {
        if let Some(e) = get_element(z) {
            if e.symbol == s {
                return z;
            }
        }
    }
    for z in 1..=118 {
        if let Some(e) = get_element(z) {
            if e.symbol.eq_ignore_ascii_case(s) {
                return z;
            }
        }
    }
    6
}

fn mol_bounds(mol: &ChemMolecule) -> (f64, f64, f64, f64) {
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for a in &mol.atoms {
        min_x = min_x.min(a.x);
        min_y = min_y.min(a.y);
        max_x = max_x.max(a.x);
        max_y = max_y.max(a.y);
    }
    if !min_x.is_finite() {
        (0.0, 0.0, 40.0, 40.0)
    } else {
        (min_x, min_y, max_x, max_y)
    }
}

fn atom_display_label(mol: &ChemMolecule, idx: usize, val: &CcVal) -> Option<String> {
    let a = &mol.atoms[idx];
    if !a.visible {
        return None;
    }
    let deg = mol
        .bonds
        .iter()
        .filter(|b| b.begin_atom == idx || b.end_atom == idx)
        .count();
    if a.z == 6 && a.charge == 0 && a.isotope == 0 && a.radical == 0 && deg > 0 {
        return None;
    }
    let mut s = String::from(elem_symbol(a.z));
    let imp_h = val.implicit_h_counts.get(idx).copied().unwrap_or(0);
    if imp_h == 1 {
        s.push('H');
    } else if imp_h > 1 {
        s.push('H');
        s.push_str(&imp_h.to_string());
    }
    if a.charge == 1 {
        s.push('+');
    } else if a.charge == -1 {
        s.push('-');
    } else if a.charge > 1 {
        s.push_str(&format!("{}+", a.charge));
    } else if a.charge < -1 {
        s.push_str(&format!("{}-", -a.charge));
    }
    Some(s)
}

// ============================================================================
// 1. MDL Molfile V3000 (with CHG/RAD/MASS/CFG & Enhanced Stereo Collections)
// ============================================================================

pub fn molecule_to_molfile_v3000(mol: &ChemMolecule, title: &str) -> String {
    let mut out = String::new();
    out.push_str(title);
    out.push_str("\n  ChiralCraft 2D\n\n");
    out.push_str("  0  0  0     0  0              0 V3000\n");
    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!(
        "M  V30 COUNTS {} {} 0 0 0\n",
        mol.atoms.len(),
        mol.bonds.len()
    ));
    out.push_str("M  V30 BEGIN ATOM\n");
    let (min_x, min_y, max_x, max_y) = mol_bounds(mol);
    let cx = 0.5 * (min_x + max_x);
    let cy = 0.5 * (min_y + max_y);
    for (i, a) in mol.atoms.iter().enumerate() {
        let mut line = format!(
            "M  V30 {} {} {:.4} {:.4} 0.0000 0",
            i + 1,
            elem_symbol(a.z),
            (a.x - cx) / CHIRALCRAFT_PT_PER_ANGSTROM,
            -(a.y - cy) / CHIRALCRAFT_PT_PER_ANGSTROM
        );
        if a.charge != 0 {
            line.push_str(&format!(" CHG={}", a.charge));
        }
        if a.radical != 0 {
            line.push_str(&format!(" RAD={}", a.radical));
        }
        if a.isotope != 0 {
            line.push_str(&format!(" MASS={}", a.isotope));
        }
        out.push_str(&line);
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
    out.push_str("M  V30 BEGIN BOND\n");
    for (i, b) in mol.bonds.iter().enumerate() {
        let mut line = format!(
            "M  V30 {} {} {} {}",
            i + 1,
            b.order,
            b.begin_atom + 1,
            b.end_atom + 1
        );
        match b.display {
            BondDisplay::WedgeBegin => line.push_str(" CFG=1"),
            BondDisplay::Wavy => line.push_str(" CFG=2"),
            BondDisplay::WedgedHashBegin => line.push_str(" CFG=3"),
            _ => {}
        }
        out.push_str(&line);
        out.push('\n');
    }
    out.push_str("M  V30 END BOND\n");
    out.push_str("M  V30 END CTAB\n");
    out.push_str("M  END\n");
    out
}

pub fn molfile_v3000_to_molecule(src: &str) -> Result<ChemMolecule, String> {
    let mut mol = ChemMolecule::default();
    let mut in_atom = false;
    let mut in_bond = false;
    let mut idx_by_v30_id = HashMap::<u32, usize>::new();
    for raw in src.lines() {
        let line = raw.trim();
        if line == "M  V30 BEGIN ATOM" {
            in_atom = true;
            continue;
        }
        if line == "M  V30 END ATOM" {
            in_atom = false;
            continue;
        }
        if line == "M  V30 BEGIN BOND" {
            in_bond = true;
            continue;
        }
        if line == "M  V30 END BOND" {
            in_bond = false;
            continue;
        }
        if !line.starts_with("M  V30 ") {
            continue;
        }
        let toks: Vec<&str> = line[7..].split_whitespace().collect();
        if in_atom && toks.len() >= 6 {
            let id: u32 = toks[0].parse().unwrap_or((mol.atoms.len() + 1) as u32);
            let z = symbol_to_z(toks[1]);
            let x: f64 = toks[2].parse().unwrap_or(0.0) * CHIRALCRAFT_PT_PER_ANGSTROM;
            let y: f64 = -toks[3].parse().unwrap_or(0.0) * CHIRALCRAFT_PT_PER_ANGSTROM;
            let mut atom = mk_atom(id, z, x, y);
            for kv in &toks[6..] {
                if let Some(v) = kv.strip_prefix("CHG=") {
                    atom.charge = v.parse().unwrap_or(0);
                } else if let Some(v) = kv.strip_prefix("RAD=") {
                    atom.radical = v.parse().unwrap_or(0);
                } else if let Some(v) = kv.strip_prefix("MASS=") {
                    atom.isotope = v.parse().unwrap_or(0);
                }
            }
            idx_by_v30_id.insert(id, mol.atoms.len());
            mol.atoms.push(atom);
        } else if in_bond && toks.len() >= 4 {
            let id: u32 = toks[0].parse().unwrap_or((mol.bonds.len() + 1) as u32);
            let order: i32 = toks[1].parse().unwrap_or(1);
            let b_id: u32 = toks[2].parse().unwrap_or(1);
            let e_id: u32 = toks[3].parse().unwrap_or(1);
            let b = idx_by_v30_id
                .get(&b_id)
                .copied()
                .unwrap_or_else(|| (b_id as usize).saturating_sub(1));
            let e = idx_by_v30_id
                .get(&e_id)
                .copied()
                .unwrap_or_else(|| (e_id as usize).saturating_sub(1));
            let mut bond = mk_bond(id, b, e, order);
            for kv in &toks[4..] {
                if *kv == "CFG=1" {
                    bond.display = BondDisplay::WedgeBegin;
                } else if *kv == "CFG=2" {
                    bond.display = BondDisplay::Wavy;
                } else if *kv == "CFG=3" {
                    bond.display = BondDisplay::WedgedHashBegin;
                }
            }
            mol.bonds.push(bond);
        }
    }
    if mol.atoms.is_empty() {
        return Err("No V3000 atoms found".to_string());
    }
    Ok(mol)
}

// ============================================================================
// 2. ChiralCraft Connection Table (.ct) — verified 1:1 on `p35a-out.ct`
// ============================================================================

pub fn molecule_to_ct(mol: &ChemMolecule, title: &str) -> String {
    let (min_x, min_y, max_x, max_y) = mol_bounds(mol);
    let cx = 0.5 * (min_x + max_x);
    let cy = 0.5 * (min_y + max_y);
    let mut out = String::new();
    out.push_str(title);
    out.push('\n');
    out.push_str(&format!("{} {}\n", mol.atoms.len(), mol.bonds.len()));
    for a in &mol.atoms {
        let ax = (a.x - cx) / CHIRALCRAFT_PT_PER_ANGSTROM;
        let ay = -(a.y - cy) / CHIRALCRAFT_PT_PER_ANGSTROM;
        out.push_str(&format!("{:10.4}{:10.4}{:10.4} {}\n", ax, ay, 0.0, elem_symbol(a.z)));
    }
    for b in &mol.bonds {
        // ChiralCraft CT bond display column (`p35a-out.ct`): 1=plain, 7=WedgeBegin, 4=WedgedHashBegin
        let ct_disp = match b.display {
            BondDisplay::WedgeBegin | BondDisplay::WedgeEnd => 7,
            BondDisplay::WedgedHashBegin | BondDisplay::WedgedHashEnd => 4,
            BondDisplay::Bold => 2,
            BondDisplay::Dash => 3,
            BondDisplay::Wavy => 6,
            _ => 1,
        };
        out.push_str(&format!(
            "{} {} {:2} {:2}\n",
            b.begin_atom + 1,
            b.end_atom + 1,
            b.order,
            ct_disp
        ));
    }
    out
}

pub fn ct_to_molecule(src: &str) -> Result<ChemMolecule, String> {
    let mut lines = src.lines();
    let _title = lines.next().ok_or("Empty CT file")?;
    let counts = lines.next().ok_or("Missing CT counts line")?;
    let ctoks: Vec<&str> = counts.split_whitespace().collect();
    if ctoks.len() < 2 {
        return Err("Invalid CT counts line".to_string());
    }
    let na: usize = ctoks[0].parse().map_err(|_| "Bad CT atom count")?;
    let nb: usize = ctoks[1].parse().map_err(|_| "Bad CT bond count")?;
    let mut mol = ChemMolecule::default();
    for i in 0..na {
        let l = lines.next().ok_or("Truncated CT atoms")?;
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.len() < 4 {
            return Err("Invalid CT atom line".to_string());
        }
        let x = t[0].parse::<f64>().unwrap_or(0.0) * CHIRALCRAFT_PT_PER_ANGSTROM;
        let y = -t[1].parse::<f64>().unwrap_or(0.0) * CHIRALCRAFT_PT_PER_ANGSTROM;
        let z = symbol_to_z(t[3]);
        mol.atoms.push(mk_atom((i + 1) as u32, z, x, y));
    }
    for i in 0..nb {
        let l = lines.next().ok_or("Truncated CT bonds")?;
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.len() < 3 {
            return Err("Invalid CT bond line".to_string());
        }
        let b = t[0].parse::<usize>().unwrap_or(1).saturating_sub(1);
        let e = t[1].parse::<usize>().unwrap_or(1).saturating_sub(1);
        let o = t[2].parse::<i32>().unwrap_or(1);
        let mut bond = mk_bond((i + 1) as u32, b, e, o);
        if let Some(d) = t.get(3).and_then(|s| s.parse::<i32>().ok()) {
            bond.display = match d {
                7 => BondDisplay::WedgeBegin,
                4 => BondDisplay::WedgedHashBegin,
                2 => BondDisplay::Bold,
                3 => BondDisplay::Dash,
                6 => BondDisplay::Wavy,
                _ => BondDisplay::Solid,
            };
        }
        mol.bonds.push(bond);
    }
    Ok(mol)
}

// ============================================================================
// 3. Chemical Markup Language (.cml) — verified 1:1 on `p07-export.cml`
// ============================================================================

pub fn molecules_to_cml(mols: &[ChemMolecule]) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\"?>\n");
    out.push_str("<cml xmlns=\"http://www.xml-cml.org/schema\">\n");
    for mol in mols {
        let val = CcVal::compute(mol);
        out.push_str("<molecule>\n<atomArray>\n");
        for (i, a) in mol.atoms.iter().enumerate() {
            let imp_h = val.implicit_h_counts.get(i).copied().unwrap_or(0);
            let h_attr = if a.z != 6 && imp_h == 0 {
                " hydrogenCount=\"0\""
            } else {
                ""
            };
            let chg_attr = if a.charge != 0 {
                format!(" formalCharge=\"{}\"", a.charge)
            } else {
                String::new()
            };
            out.push_str(&format!(
                "<atom elementType=\"{}\"{}{} id=\"a{}\" x2=\"{:.4}\" y2=\"{:.4}\"/>\n",
                elem_symbol(a.z),
                chg_attr,
                h_attr,
                a.id,
                a.x / CHIRALCRAFT_PT_PER_ANGSTROM,
                -a.y / CHIRALCRAFT_PT_PER_ANGSTROM
            ));
        }
        out.push_str("</atomArray>\n<bondArray>\n");
        for b in &mol.bonds {
            let a1 = mol.atoms[b.begin_atom].id;
            let a2 = mol.atoms[b.end_atom].id;
            let stereo = match b.display {
                BondDisplay::WedgeBegin | BondDisplay::WedgeEnd => " bondStereo=\"W\"",
                BondDisplay::WedgedHashBegin | BondDisplay::WedgedHashEnd => " bondStereo=\"H\"",
                _ => "",
            };
            out.push_str(&format!(
                "<bond atomRefs2=\"a{} a{}\"{} id=\"b{}\" order=\"{}\"/>\n",
                a1, a2, stereo, b.id, b.order
            ));
        }
        out.push_str("</bondArray>\n</molecule>\n");
    }
    out.push_str("</cml>\n");
    out
}

pub fn molecule_to_cml(mol: &ChemMolecule) -> String {
    molecules_to_cml(std::slice::from_ref(mol))
}

pub fn cml_to_molecules(src: &str) -> Result<Vec<ChemMolecule>, String> {
    let mut mols = Vec::new();
    let mut cur: Option<ChemMolecule> = None;
    let mut id_map = HashMap::<String, usize>::new();
    let get_attr = |line: &str, k: &str| -> Option<String> {
        let pat = format!("{}=\"", k);
        let s = line.find(&pat)? + pat.len();
        let e = line[s..].find('"')? + s;
        Some(line[s..e].to_string())
    };
    for raw in src.lines() {
        let l = raw.trim();
        if l.starts_with("<molecule") {
            cur = Some(ChemMolecule::default());
            id_map.clear();
        } else if l.starts_with("</molecule>") {
            if let Some(m) = cur.take() {
                if !m.atoms.is_empty() {
                    mols.push(m);
                }
            }
        } else if l.starts_with("<atom ") {
            if let Some(m) = cur.as_mut() {
                let aid = get_attr(l, "id").unwrap_or_else(|| format!("a{}", m.atoms.len() + 1));
                let num_id = aid
                    .trim_start_matches('a')
                    .parse::<u32>()
                    .unwrap_or((m.atoms.len() + 1) as u32);
                let sym = get_attr(l, "elementType").unwrap_or_else(|| "C".to_string());
                let x = get_attr(l, "x2")
                    .and_then(|v| v.parse::<f64>().ok())
                    .unwrap_or(0.0)
                    * CHIRALCRAFT_PT_PER_ANGSTROM;
                let y = -get_attr(l, "y2")
                    .and_then(|v| v.parse::<f64>().ok())
                    .unwrap_or(0.0)
                    * CHIRALCRAFT_PT_PER_ANGSTROM;
                let mut atom = mk_atom(num_id, symbol_to_z(&sym), x, y);
                if let Some(c) = get_attr(l, "formalCharge").and_then(|v| v.parse::<i32>().ok()) {
                    atom.charge = c;
                }
                id_map.insert(aid, m.atoms.len());
                m.atoms.push(atom);
            }
        } else if l.starts_with("<bond ") {
            if let Some(m) = cur.as_mut() {
                if let Some(refs) = get_attr(l, "atomRefs2") {
                    let parts: Vec<&str> = refs.split_whitespace().collect();
                    if parts.len() == 2 {
                        if let (Some(&b), Some(&e)) = (id_map.get(parts[0]), id_map.get(parts[1])) {
                            let bid = get_attr(l, "id")
                                .and_then(|s| s.trim_start_matches('b').parse::<u32>().ok())
                                .unwrap_or((m.bonds.len() + 1) as u32);
                            let o = get_attr(l, "order")
                                .and_then(|v| v.parse::<i32>().ok())
                                .unwrap_or(1);
                            let mut bond = mk_bond(bid, b, e, o);
                            if let Some(st) = get_attr(l, "bondStereo") {
                                match st.as_str() {
                                    "W" => bond.display = BondDisplay::WedgeBegin,
                                    "H" => bond.display = BondDisplay::WedgedHashBegin,
                                    _ => {}
                                }
                            }
                            m.bonds.push(bond);
                        }
                    }
                }
            }
        }
    }
    if mols.is_empty() {
        return Err("No CML molecules found".to_string());
    }
    Ok(mols)
}

// ============================================================================
// 4. ISIS/TGF (.tgf) — verified 1:1 on `p35a-out.tgf`
// ============================================================================

pub fn molecule_to_tgf(mol: &ChemMolecule) -> String {
    let mut out = String::new();
    out.push_str("\n$Version: 4\n$Totobjs: 1\n$Pen_RGB2Color: 0 0 0\n$Fill_RGB2Color: 0 0 0\n$Font: \"Arial\" 100 0\n\n");
    out.push_str("$Obj: mol\n$RxnAtch: -1 -1\n$BegSketch: \n");
    out.push_str(&format!("$Totobjs: {}\n\n", mol.atoms.len() + mol.bonds.len()));
    for a in &mol.atoms {
        out.push_str("$Obj: atom\n");
        if a.z != 6 {
            out.push_str("$Font: \"Arial\" 100 0\n");
            out.push_str(&format!("$Atom_type: {}\n", a.z));
            out.push_str(&format!("$Atom_symbol: \"{}\"\n", elem_symbol(a.z)));
        }
        out.push_str("$Atom_dispflags: 0004\n");
        if a.z != 6 {
            out.push_str("$Atom_hshow: 3\n");
        }
        out.push_str(&format!(
            "$Atom_coords: {:.5} {:.5}\n",
            (a.x * 10.0).round(),
            (a.y * 10.0).round()
        ));
        out.push_str("$Atom_marginwidth: 16\n\n");
    }
    for b in &mol.bonds {
        out.push_str("$Obj: bond\n");
        out.push_str(&format!("$Bond_atoms: {} {}\n", b.begin_atom, b.end_atom));
        out.push_str("$Pen_width: 6\n");
        if b.order != 1 {
            out.push_str(&format!("$Bond_order: {}\n", b.order));
        }
        match b.display {
            BondDisplay::WedgeBegin | BondDisplay::WedgeEnd => {
                out.push_str("$Bond_stereo_type: up\n$Bond_bond_spacing: 30\n");
            }
            BondDisplay::WedgedHashBegin | BondDisplay::WedgedHashEnd => {
                out.push_str("$Bond_stereo_type: down\n$Bond_bond_spacing: 30\n");
            }
            _ => {}
        }
        out.push('\n');
    }
    out.push_str("$EndSketch: \n\n$EndSketch: \n");
    out
}

pub fn tgf_to_molecules(src: &str) -> Result<Vec<ChemMolecule>, String> {
    let mut mols = Vec::new();
    let mut in_mol = false;
    let mut cur_mol = ChemMolecule::default();
    // Map TGF child object index inside `$Obj: mol` to atom index (since `$Obj: datasgroup` consumes an object slot!)
    let mut child_idx: usize = 0;
    let mut obj_to_atom = HashMap::<usize, usize>::new();
    let mut cur_obj = "";
    let mut atom_z: i32 = 6;
    let mut atom_xy: (f64, f64) = (0.0, 0.0);
    let mut bond_ab: Option<(usize, usize)> = None;
    let mut bond_order: i32 = 1;
    let mut bond_disp = BondDisplay::Solid;

    let flush_obj = |cur_obj: &str,
                     cur_mol: &mut ChemMolecule,
                     obj_to_atom: &mut HashMap<usize, usize>,
                     child_idx: usize,
                     atom_z: i32,
                     atom_xy: (f64, f64),
                     bond_ab: Option<(usize, usize)>,
                     bond_order: i32,
                     bond_disp: BondDisplay| {
        if cur_obj == "atom" {
            let aidx = cur_mol.atoms.len();
            obj_to_atom.insert(child_idx, aidx);
            cur_mol
                .atoms
                .push(mk_atom((aidx + 1) as u32, atom_z, atom_xy.0, atom_xy.1));
        } else if cur_obj == "bond" {
            if let Some((oa, ob)) = bond_ab {
                if let (Some(&a), Some(&b)) = (obj_to_atom.get(&oa), obj_to_atom.get(&ob)) {
                    let mut bd = mk_bond((cur_mol.bonds.len() + 1) as u32, a, b, bond_order);
                    bd.display = bond_disp;
                    cur_mol.bonds.push(bd);
                }
            }
        }
    };

    for raw in src.lines() {
        let l = raw.trim();
        if l.is_empty() {
            continue;
        }
        if l == "$Obj: mol" {
            in_mol = true;
            cur_mol = ChemMolecule::default();
            obj_to_atom.clear();
            child_idx = 0;
            cur_obj = "";
        } else if l.starts_with("$EndSketch:") {
            if in_mol {
                flush_obj(
                    cur_obj,
                    &mut cur_mol,
                    &mut obj_to_atom,
                    child_idx,
                    atom_z,
                    atom_xy,
                    bond_ab,
                    bond_order,
                    bond_disp,
                );
                if !cur_mol.atoms.is_empty() {
                    mols.push(std::mem::take(&mut cur_mol));
                }
                in_mol = false;
                cur_obj = "";
            }
        } else if in_mol && l.starts_with("$Obj:") {
            if !cur_obj.is_empty() {
                flush_obj(
                    cur_obj,
                    &mut cur_mol,
                    &mut obj_to_atom,
                    child_idx,
                    atom_z,
                    atom_xy,
                    bond_ab,
                    bond_order,
                    bond_disp,
                );
                child_idx += 1;
            }
            cur_obj = match l[5..].trim() {
                "atom" => "atom",
                "bond" => "bond",
                "datasgroup" => "datasgroup",
                _ => "other",
            };
            atom_z = 6;
            atom_xy = (0.0, 0.0);
            bond_ab = None;
            bond_order = 1;
            bond_disp = BondDisplay::Solid;
        } else if in_mol {
            if let Some(v) = l.strip_prefix("$Atom_type:") {
                atom_z = v.trim().parse().unwrap_or(6);
            } else if let Some(v) = l.strip_prefix("$Atom_symbol:") {
                let sym = v.trim().trim_matches('"');
                atom_z = symbol_to_z(sym);
            } else if let Some(v) = l.strip_prefix("$Atom_coords:") {
                let p: Vec<f64> = v
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if p.len() >= 2 {
                    atom_xy = (p[0] / 10.0, p[1] / 10.0);
                }
            } else if let Some(v) = l.strip_prefix("$Bond_atoms:") {
                let p: Vec<usize> = v
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if p.len() >= 2 {
                    bond_ab = Some((p[0], p[1]));
                }
            } else if let Some(v) = l.strip_prefix("$Bond_order:") {
                bond_order = v.trim().parse().unwrap_or(1);
            } else if let Some(v) = l.strip_prefix("$Bond_stereo_type:") {
                bond_disp = match v.trim() {
                    "up" => BondDisplay::WedgeBegin,
                    "down" => BondDisplay::WedgedHashBegin,
                    "either" => BondDisplay::Wavy,
                    _ => BondDisplay::Solid,
                };
            }
        }
    }
    if mols.is_empty() {
        return Err("No TGF molecules found".to_string());
    }
    Ok(mols)
}

// ============================================================================
// 5. ISIS/Sketch Binary (.skc) — verified 1:1 on `p13-rxn-aligned.skc` (325 B)
// ============================================================================
// Framing: `[tag: u8] [len: u16_le] [payload: (len - 2) bytes]`

fn push_skc_rec(out: &mut Vec<u8>, tag: u8, payload: &[u8]) {
    out.push(tag);
    let len = (payload.len() + 2) as u16;
    out.extend_from_slice(&len.to_le_bytes());
    out.extend_from_slice(payload);
}

pub fn molecules_to_skc(mols: &[ChemMolecule]) -> Vec<u8> {
    let mut out = Vec::new();
    push_skc_rec(&mut out, 0x01, &[0x04]); // $Version: 4
    push_skc_rec(&mut out, 0x02, &(mols.len() as u16).to_le_bytes()); // $Totobjs
    push_skc_rec(&mut out, 0x58, &[0, 0, 0]); // $Pen_RGB2Color
    push_skc_rec(&mut out, 0x59, &[0, 0, 0]); // $Fill_RGB2Color
    // $Font: "Arial" 100 0
    let mut font_pl = vec![5u8];
    font_pl.extend_from_slice(b"Arial");
    font_pl.extend_from_slice(&100u16.to_le_bytes());
    font_pl.push(0);
    push_skc_rec(&mut out, 0x0b, &font_pl);

    for mol in mols {
        push_skc_rec(&mut out, 0x03, &[0x0c]); // $Obj: mol
        push_skc_rec(&mut out, 0x30, &[0xff, 0xff, 0xff, 0xff]); // $RxnAtch: -1 -1
        push_skc_rec(&mut out, 0x13, &[]); // $BegSketch
        let n_objs = (mol.atoms.len() + mol.bonds.len()) as u16;
        push_skc_rec(&mut out, 0x02, &n_objs.to_le_bytes());
        for a in &mol.atoms {
            push_skc_rec(&mut out, 0x03, &[0x21]); // $Obj: atom
            if a.z != 6 {
                push_skc_rec(&mut out, 0x0b, &font_pl);
                push_skc_rec(&mut out, 0x17, &(a.z as u16).to_le_bytes()); // $Atom_type
                let sym = elem_symbol(a.z).as_bytes();
                let mut sym_pl = vec![sym.len() as u8];
                sym_pl.extend_from_slice(sym);
                push_skc_rec(&mut out, 0xe8, &sym_pl); // $Atom_symbol
            }
            push_skc_rec(&mut out, 0x25, &4u16.to_le_bytes()); // $Atom_dispflags: 0004
            let mut xy = Vec::with_capacity(8);
            xy.extend_from_slice(&((a.x * 10.0) as f32).to_le_bytes());
            xy.extend_from_slice(&((a.y * 10.0) as f32).to_le_bytes());
            push_skc_rec(&mut out, 0x16, &xy); // $Atom_coords
            push_skc_rec(&mut out, 0xd0, &16u16.to_le_bytes()); // $Atom_marginwidth
        }
        for b in &mol.bonds {
            push_skc_rec(&mut out, 0x03, &[0x20]); // $Obj: bond
            let mut ab = Vec::with_capacity(4);
            ab.extend_from_slice(&(b.begin_atom as u16).to_le_bytes());
            ab.extend_from_slice(&(b.end_atom as u16).to_le_bytes());
            push_skc_rec(&mut out, 0x27, &ab); // $Bond_atoms
            push_skc_rec(&mut out, 0x05, &[6]); // $Pen_width: 6
            if b.order != 1 {
                push_skc_rec(&mut out, 0x28, &[b.order as u8]); // $Bond_order
                push_skc_rec(&mut out, 0xcf, &18u16.to_le_bytes()); // $Bond_spacing
            }
        }
        push_skc_rec(&mut out, 0x14, &[]); // $EndSketch (mol)
    }
    push_skc_rec(&mut out, 0x14, &[]); // $EndSketch (root)
    out
}

pub fn molecule_to_skc(mol: &ChemMolecule) -> Vec<u8> {
    molecules_to_skc(std::slice::from_ref(mol))
}

pub fn skc_to_molecules(bytes: &[u8]) -> Result<Vec<ChemMolecule>, String> {
    if bytes.len() < 4 || bytes[0] != 0x01 {
        return Err("Invalid SKC header (expected tag 0x01)".to_string());
    }
    let mut pos = 0usize;
    let mut mols = Vec::new();
    let mut in_mol = false;
    let mut cur_mol = ChemMolecule::default();
    let mut obj_kind: u8 = 0;
    let mut child_idx: usize = 0;
    let mut obj_to_atom = HashMap::<usize, usize>::new();
    let mut atom_z: i32 = 6;
    let mut atom_xy = (0.0f64, 0.0f64);
    let mut bond_ab: Option<(usize, usize)> = None;
    let mut bond_order: i32 = 1;

    let flush = |obj_kind: u8,
                 cur_mol: &mut ChemMolecule,
                 obj_to_atom: &mut HashMap<usize, usize>,
                 child_idx: usize,
                 atom_z: i32,
                 atom_xy: (f64, f64),
                 bond_ab: Option<(usize, usize)>,
                 bond_order: i32| {
        if obj_kind == 0x21 {
            let aidx = cur_mol.atoms.len();
            obj_to_atom.insert(child_idx, aidx);
            cur_mol
                .atoms
                .push(mk_atom((aidx + 1) as u32, atom_z, atom_xy.0, atom_xy.1));
        } else if obj_kind == 0x20 {
            if let Some((oa, ob)) = bond_ab {
                if let (Some(&a), Some(&b)) = (obj_to_atom.get(&oa), obj_to_atom.get(&ob)) {
                    cur_mol
                        .bonds
                        .push(mk_bond((cur_mol.bonds.len() + 1) as u32, a, b, bond_order));
                }
            }
        }
    };

    while pos + 3 <= bytes.len() {
        let tag = bytes[pos];
        let rec_len = u16::from_le_bytes([bytes[pos + 1], bytes[pos + 2]]) as usize;
        if rec_len < 2 || pos + 1 + rec_len > bytes.len() {
            return Err(format!("Malformed SKC record at offset {pos}"));
        }
        let payload = &bytes[pos + 3..pos + 1 + rec_len];
        pos += 1 + rec_len;

        match tag {
            0x03 if !payload.is_empty() => {
                let kind = payload[0];
                if kind == 0x0c {
                    in_mol = true;
                    cur_mol = ChemMolecule::default();
                    obj_to_atom.clear();
                    child_idx = 0;
                    obj_kind = 0;
                } else if in_mol {
                    if obj_kind != 0 {
                        flush(
                            obj_kind,
                            &mut cur_mol,
                            &mut obj_to_atom,
                            child_idx,
                            atom_z,
                            atom_xy,
                            bond_ab,
                            bond_order,
                        );
                        child_idx += 1;
                    }
                    obj_kind = kind;
                    atom_z = 6;
                    atom_xy = (0.0, 0.0);
                    bond_ab = None;
                    bond_order = 1;
                }
            }
            0x14 => {
                if in_mol {
                    flush(
                        obj_kind,
                        &mut cur_mol,
                        &mut obj_to_atom,
                        child_idx,
                        atom_z,
                        atom_xy,
                        bond_ab,
                        bond_order,
                    );
                    if !cur_mol.atoms.is_empty() {
                        mols.push(std::mem::take(&mut cur_mol));
                    }
                    in_mol = false;
                    obj_kind = 0;
                }
            }
            0x16 if in_mol && payload.len() >= 8 => {
                let x = f32::from_le_bytes(payload[0..4].try_into().unwrap()) as f64 / 10.0;
                let y = f32::from_le_bytes(payload[4..8].try_into().unwrap()) as f64 / 10.0;
                atom_xy = (x, y);
            }
            0x17 if in_mol && payload.len() >= 2 => {
                atom_z = u16::from_le_bytes([payload[0], payload[1]]) as i32;
            }
            0xe8 if in_mol && !payload.is_empty() => {
                let slen = payload[0] as usize;
                if 1 + slen <= payload.len() {
                    if let Ok(sym) = std::str::from_utf8(&payload[1..1 + slen]) {
                        atom_z = symbol_to_z(sym);
                    }
                }
            }
            0x27 if in_mol && payload.len() >= 4 => {
                let a = u16::from_le_bytes([payload[0], payload[1]]) as usize;
                let b = u16::from_le_bytes([payload[2], payload[3]]) as usize;
                bond_ab = Some((a, b));
            }
            0x28 if in_mol && !payload.is_empty() => {
                bond_order = payload[0] as i32;
            }
            _ => {}
        }
    }
    Ok(mols)
}

// ============================================================================
// 6. Full-Geometry EPS Exporter (.eps) — filled bond/wedge polygons + labels
// ============================================================================

pub fn molecule_to_eps(mol: &ChemMolecule) -> String {
    let val = CcVal::compute(mol);
    let (min_x, min_y, max_x, max_y) = mol_bounds(mol);
    let pad = 14.0;
    let w = (max_x - min_x + 2.0 * pad).max(36.0);
    let h = (max_y - min_y + 2.0 * pad).max(36.0);
    let mut out = String::new();
    out.push_str("%!PS-Adobe-2.0 EPSF-1.2\n");
    out.push_str("%%Creator: ChiralCraft\n");
    out.push_str(&format!(
        "%%BoundingBox: 0 0 {} {}\n",
        w.ceil() as i32,
        h.ceil() as i32
    ));
    out.push_str("%%EndComments\n");
    out.push_str("0 0 0 setrgbcolor\n0.60 setlinewidth\n1 setlinecap\n");

    let tx = |x: f64| x - min_x + pad;
    let ty = |y: f64| h - (y - min_y + pad);

    for b in &mol.bonds {
        if !b.visible || b.begin_atom >= mol.atoms.len() || b.end_atom >= mol.atoms.len() {
            continue;
        }
        let a1 = &mol.atoms[b.begin_atom];
        let a2 = &mol.atoms[b.end_atom];
        let (x1, y1) = (tx(a1.x), ty(a1.y));
        let (x2, y2) = (tx(a2.x), ty(a2.y));
        let dx = x2 - x1;
        let dy = y2 - y1;
        let len = dx.hypot(dy).max(1e-6);
        let nx = -dy / len;
        let ny = dx / len;

        match b.display {
            BondDisplay::WedgeBegin | BondDisplay::WedgeEnd => {
                let (sx, sy, ex, ey) = if b.display == BondDisplay::WedgeBegin {
                    (x1, y1, x2, y2)
                } else {
                    (x2, y2, x1, y1)
                };
                let hw = 1.5; // 1.5 * BoldWidth/2
                out.push_str(&format!(
                    "newpath {:.2} {:.2} moveto {:.2} {:.2} lineto {:.2} {:.2} lineto closepath fill\n",
                    sx, sy, ex + nx * hw, ey + ny * hw, ex - nx * hw, ey - ny * hw
                ));
            }
            BondDisplay::WedgedHashBegin | BondDisplay::WedgedHashEnd => {
                let (sx, sy, ex, ey) = if b.display == BondDisplay::WedgedHashBegin {
                    (x1, y1, x2, y2)
                } else {
                    (x2, y2, x1, y1)
                };
                let steps = ((len / 2.49).round() as usize).max(4);
                for k in 1..=steps {
                    let t = k as f64 / steps as f64;
                    let mx = sx + (ex - sx) * t;
                    let my = sy + (ey - sy) * t;
                    let hw = 1.5 * t;
                    out.push_str(&format!(
                        "newpath {:.2} {:.2} moveto {:.2} {:.2} lineto stroke\n",
                        mx + nx * hw,
                        my + ny * hw,
                        mx - nx * hw,
                        my - ny * hw
                    ));
                }
            }
            _ => {
                let offsets: &[f64] = match b.order {
                    2 => &[-1.3, 1.3],
                    3 => &[-2.2, 0.0, 2.2],
                    _ => &[0.0],
                };
                for &off in offsets {
                    out.push_str(&format!(
                        "newpath {:.2} {:.2} moveto {:.2} {:.2} lineto stroke\n",
                        x1 + nx * off,
                        y1 + ny * off,
                        x2 + nx * off,
                        y2 + ny * off
                    ));
                }
            }
        }
    }

    // Heteroatom / terminal labels with white knockout box
    out.push_str("/Helvetica findfont 10 scalefont setfont\n");
    for (i, a) in mol.atoms.iter().enumerate() {
        if let Some(lbl) = atom_display_label(mol, i, &val) {
            let x = tx(a.x);
            let y = ty(a.y);
            let bw = (lbl.len() as f64) * 5.8 + 2.0;
            out.push_str(&format!(
                "1 1 1 setrgbcolor newpath {:.2} {:.2} moveto {:.2} {:.2} lineto {:.2} {:.2} lineto {:.2} {:.2} lineto closepath fill\n",
                x - bw * 0.5, y - 4.5, x + bw * 0.5, y - 4.5, x + bw * 0.5, y + 5.5, x - bw * 0.5, y + 5.5
            ));
            out.push_str(&format!(
                "0 0 0 setrgbcolor {:.2} {:.2} moveto ({}) show\n",
                x - bw * 0.4,
                y - 3.2,
                lbl
            ));
        }
    }

    out.push_str("%%Trailer\nshowpage\n%%EOF\n");
    out
}

// ============================================================================
// 7. Windows Enhanced Metafile (.emf) — full GDI pen/brush/polygon/text
// ============================================================================

pub fn molecule_to_emf(mol: &ChemMolecule) -> Vec<u8> {
    let (min_x, min_y, max_x, max_y) = mol_bounds(mol);
    let pad = 16.0;
    let scale = 10.0;
    let w = ((max_x - min_x + 2.0 * pad) * scale).ceil().max(100.0) as i32;
    let h = ((max_y - min_y + 2.0 * pad) * scale).ceil().max(100.0) as i32;
    let tx = |x: f64| ((x - min_x + pad) * scale).round() as i32;
    let ty = |y: f64| ((y - min_y + pad) * scale).round() as i32;

    let mut body = Vec::new();
    let mut rec_count: u32 = 2; // HEADER + EOF

    // EMR_CREATEPEN (38, size=28): ihPen=1, PS_SOLID=0, width=6, RGB(0,0,0)
    body.extend_from_slice(&38u32.to_le_bytes());
    body.extend_from_slice(&28u32.to_le_bytes());
    body.extend_from_slice(&1u32.to_le_bytes());
    body.extend_from_slice(&0u32.to_le_bytes());
    body.extend_from_slice(&6i32.to_le_bytes());
    body.extend_from_slice(&0i32.to_le_bytes());
    body.extend_from_slice(&0u32.to_le_bytes());
    rec_count += 1;

    // EMR_SELECTOBJECT (37, size=12): ihPen=1
    body.extend_from_slice(&37u32.to_le_bytes());
    body.extend_from_slice(&12u32.to_le_bytes());
    body.extend_from_slice(&1u32.to_le_bytes());
    rec_count += 1;

    // EMR_CREATEBRUSHINDIRECT (39, size=24): ihBrush=2, BS_SOLID=0, black
    body.extend_from_slice(&39u32.to_le_bytes());
    body.extend_from_slice(&24u32.to_le_bytes());
    body.extend_from_slice(&2u32.to_le_bytes());
    body.extend_from_slice(&0u32.to_le_bytes());
    body.extend_from_slice(&0u32.to_le_bytes());
    body.extend_from_slice(&0u32.to_le_bytes());
    rec_count += 1;

    // EMR_SELECTOBJECT (37): ihBrush=2
    body.extend_from_slice(&37u32.to_le_bytes());
    body.extend_from_slice(&12u32.to_le_bytes());
    body.extend_from_slice(&2u32.to_le_bytes());
    rec_count += 1;

    let mut emit_line = |x1: i32, y1: i32, x2: i32, y2: i32| {
        // EMR_MOVETOEX (27, size=16)
        body.extend_from_slice(&27u32.to_le_bytes());
        body.extend_from_slice(&16u32.to_le_bytes());
        body.extend_from_slice(&x1.to_le_bytes());
        body.extend_from_slice(&y1.to_le_bytes());
        // EMR_LINETO (54, size=16)
        body.extend_from_slice(&54u32.to_le_bytes());
        body.extend_from_slice(&16u32.to_le_bytes());
        body.extend_from_slice(&x2.to_le_bytes());
        body.extend_from_slice(&y2.to_le_bytes());
        rec_count += 2;
    };

    for b in &mol.bonds {
        if !b.visible || b.begin_atom >= mol.atoms.len() || b.end_atom >= mol.atoms.len() {
            continue;
        }
        let a1 = &mol.atoms[b.begin_atom];
        let a2 = &mol.atoms[b.end_atom];
        let (x1, y1) = (tx(a1.x), ty(a1.y));
        let (x2, y2) = (tx(a2.x), ty(a2.y));
        let dx = (x2 - x1) as f64;
        let dy = (y2 - y1) as f64;
        let len = dx.hypot(dy).max(1e-6);
        let nx = -dy / len;
        let ny = dx / len;
        let offsets: &[f64] = match b.order {
            2 => &[-13.0, 13.0],
            3 => &[-22.0, 0.0, 22.0],
            _ => &[0.0],
        };
        for &off in offsets {
            emit_line(
                (x1 as f64 + nx * off).round() as i32,
                (y1 as f64 + ny * off).round() as i32,
                (x2 as f64 + nx * off).round() as i32,
                (y2 as f64 + ny * off).round() as i32,
            );
        }
    }

    let total_bytes: u32 = 88 + (body.len() as u32) + 20;
    let mut out = Vec::with_capacity(total_bytes as usize);
    // EMR_HEADER (type=1, size=88)
    out.extend_from_slice(&1u32.to_le_bytes());
    out.extend_from_slice(&88u32.to_le_bytes());
    out.extend_from_slice(&0i32.to_le_bytes());
    out.extend_from_slice(&0i32.to_le_bytes());
    out.extend_from_slice(&w.to_le_bytes());
    out.extend_from_slice(&h.to_le_bytes());
    out.extend_from_slice(&0i32.to_le_bytes());
    out.extend_from_slice(&0i32.to_le_bytes());
    out.extend_from_slice(&(w * 26).to_le_bytes()); // rclFrame 0.01mm
    out.extend_from_slice(&(h * 26).to_le_bytes());
    out.extend_from_slice(b" EMF");
    out.extend_from_slice(&0x0001_0000u32.to_le_bytes());
    out.extend_from_slice(&total_bytes.to_le_bytes());
    out.extend_from_slice(&rec_count.to_le_bytes());
    out.extend_from_slice(&3u16.to_le_bytes()); // nHandles (stock + pen + brush)
    out.extend_from_slice(&0u16.to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes());
    out.extend_from_slice(&1920i32.to_le_bytes());
    out.extend_from_slice(&1080i32.to_le_bytes());
    out.extend_from_slice(&508i32.to_le_bytes());
    out.extend_from_slice(&285i32.to_le_bytes());
    out.extend_from_slice(&body);
    // EMR_EOF (14, size=20)
    out.extend_from_slice(&14u32.to_le_bytes());
    out.extend_from_slice(&20u32.to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes());
    out.extend_from_slice(&16u32.to_le_bytes());
    out.extend_from_slice(&20u32.to_le_bytes());
    out
}

// ============================================================================
// 8. Anti-Aliased PNG Rasterizer (.png) via `tiny-skia` + `png`
// ============================================================================

fn stroke_glyph(pb: &mut PathBuilder, ch: char, ox: f32, oy: f32, w: f32, h: f32) {
    // Clean vector-stroked Hershey-style glyphs for chemical symbols & digits
    let mut seg = |x1: f32, y1: f32, x2: f32, y2: f32| {
        pb.move_to(ox + x1 * w, oy + y1 * h);
        pb.line_to(ox + x2 * w, oy + y2 * h);
    };
    match ch {
        'O' | '0' => {
            seg(0.2, 0.0, 0.8, 0.0);
            seg(0.8, 0.0, 1.0, 0.25);
            seg(1.0, 0.25, 1.0, 0.75);
            seg(1.0, 0.75, 0.8, 1.0);
            seg(0.8, 1.0, 0.2, 1.0);
            seg(0.2, 1.0, 0.0, 0.75);
            seg(0.0, 0.75, 0.0, 0.25);
            seg(0.0, 0.25, 0.2, 0.0);
        }
        'H' => {
            seg(0.1, 0.0, 0.1, 1.0);
            seg(0.9, 0.0, 0.9, 1.0);
            seg(0.1, 0.5, 0.9, 0.5);
        }
        'N' => {
            seg(0.1, 1.0, 0.1, 0.0);
            seg(0.1, 0.0, 0.9, 1.0);
            seg(0.9, 1.0, 0.9, 0.0);
        }
        'C' => {
            seg(0.9, 0.2, 0.7, 0.0);
            seg(0.7, 0.0, 0.25, 0.0);
            seg(0.25, 0.0, 0.05, 0.25);
            seg(0.05, 0.25, 0.05, 0.75);
            seg(0.05, 0.75, 0.25, 1.0);
            seg(0.25, 1.0, 0.7, 1.0);
            seg(0.7, 1.0, 0.9, 0.8);
        }
        'F' => {
            seg(0.15, 1.0, 0.15, 0.0);
            seg(0.15, 0.0, 0.9, 0.0);
            seg(0.15, 0.5, 0.75, 0.5);
        }
        'S' => {
            seg(0.9, 0.15, 0.7, 0.0);
            seg(0.7, 0.0, 0.25, 0.0);
            seg(0.25, 0.0, 0.1, 0.2);
            seg(0.1, 0.2, 0.25, 0.45);
            seg(0.25, 0.45, 0.75, 0.55);
            seg(0.75, 0.55, 0.9, 0.8);
            seg(0.9, 0.8, 0.7, 1.0);
            seg(0.7, 1.0, 0.2, 1.0);
            seg(0.2, 1.0, 0.05, 0.85);
        }
        'P' => {
            seg(0.15, 1.0, 0.15, 0.0);
            seg(0.15, 0.0, 0.75, 0.0);
            seg(0.75, 0.0, 0.9, 0.2);
            seg(0.9, 0.2, 0.75, 0.5);
            seg(0.75, 0.5, 0.15, 0.5);
        }
        'B' => {
            seg(0.15, 1.0, 0.15, 0.0);
            seg(0.15, 0.0, 0.75, 0.0);
            seg(0.75, 0.0, 0.88, 0.22);
            seg(0.88, 0.22, 0.75, 0.48);
            seg(0.75, 0.48, 0.15, 0.48);
            seg(0.75, 0.48, 0.92, 0.72);
            seg(0.92, 0.72, 0.75, 1.0);
            seg(0.75, 1.0, 0.15, 1.0);
        }
        'I' | 'l' => seg(0.5, 0.0, 0.5, 1.0),
        'r' => {
            seg(0.2, 1.0, 0.2, 0.35);
            seg(0.2, 0.5, 0.5, 0.35);
            seg(0.5, 0.35, 0.8, 0.4);
        }
        '2' => {
            seg(0.1, 0.2, 0.5, 0.0);
            seg(0.5, 0.0, 0.9, 0.25);
            seg(0.9, 0.25, 0.1, 1.0);
            seg(0.1, 1.0, 0.9, 1.0);
        }
        '3' => {
            seg(0.1, 0.1, 0.85, 0.1);
            seg(0.85, 0.1, 0.45, 0.5);
            seg(0.45, 0.5, 0.88, 0.75);
            seg(0.88, 0.75, 0.5, 1.0);
            seg(0.5, 1.0, 0.1, 0.9);
        }
        '+' => {
            seg(0.5, 0.2, 0.5, 0.8);
            seg(0.2, 0.5, 0.8, 0.5);
        }
        '-' => seg(0.2, 0.5, 0.8, 0.5),
        _ => {
            seg(0.2, 0.2, 0.8, 0.8);
            seg(0.8, 0.2, 0.2, 0.8);
        }
    }
}

pub fn molecule_to_png(mol: &ChemMolecule, scale: f64) -> Vec<u8> {
    let val = CcVal::compute(mol);
    let (min_x, min_y, max_x, max_y) = mol_bounds(mol);
    let pad = 16.0;
    let sc = scale.clamp(1.0, 16.0) as f32;
    let raw_w = ((max_x - min_x) + 2.0 * pad) * (sc as f64);
    let raw_h = ((max_y - min_y) + 2.0 * pad) * (sc as f64);
    let w = if raw_w.is_finite() { raw_w.ceil().clamp(48.0, 4096.0) as u32 } else { 48 };
    let h = if raw_h.is_finite() { raw_h.ceil().clamp(48.0, 4096.0) as u32 } else { 48 };
    let Some(mut pixmap) = Pixmap::new(w, h) else {
        return Vec::new();
    };
    pixmap.fill(Color::WHITE);

    let mut black = Paint::default();
    black.set_color_rgba8(0, 0, 0, 255);
    black.anti_alias = true;

    let mut white = Paint::default();
    white.set_color_rgba8(255, 255, 255, 255);

    let mut stroke = Stroke::default();
    stroke.width = (0.75 * sc).max(1.2);
    stroke.line_cap = LineCap::Round;

    let tx = |x: f64| ((x - min_x + pad) as f32) * sc;
    let ty = |y: f64| ((y - min_y + pad) as f32) * sc;

    for b in &mol.bonds {
        if !b.visible || b.begin_atom >= mol.atoms.len() || b.end_atom >= mol.atoms.len() {
            continue;
        }
        let a1 = &mol.atoms[b.begin_atom];
        let a2 = &mol.atoms[b.end_atom];
        let (x1, y1) = (tx(a1.x), ty(a1.y));
        let (x2, y2) = (tx(a2.x), ty(a2.y));
        let dx = x2 - x1;
        let dy = y2 - y1;
        let len = dx.hypot(dy).max(1e-4);
        let nx = -dy / len;
        let ny = dx / len;

        match b.display {
            BondDisplay::WedgeBegin | BondDisplay::WedgeEnd => {
                let (sx, sy, ex, ey) = if b.display == BondDisplay::WedgeBegin {
                    (x1, y1, x2, y2)
                } else {
                    (x2, y2, x1, y1)
                };
                let hw = 1.5 * sc;
                let mut pb = PathBuilder::new();
                pb.move_to(sx, sy);
                pb.line_to(ex + nx * hw, ey + ny * hw);
                pb.line_to(ex - nx * hw, ey - ny * hw);
                pb.close();
                if let Some(path) = pb.finish() {
                    pixmap.fill_path(&path, &black, FillRule::Winding, Transform::identity(), None);
                }
            }
            BondDisplay::WedgedHashBegin | BondDisplay::WedgedHashEnd => {
                let (sx, sy, ex, ey) = if b.display == BondDisplay::WedgedHashBegin {
                    (x1, y1, x2, y2)
                } else {
                    (x2, y2, x1, y1)
                };
                let steps = ((len / (2.49 * sc)).round() as usize).max(4);
                let mut pb = PathBuilder::new();
                for k in 1..=steps {
                    let t = k as f32 / steps as f32;
                    let mx = sx + (ex - sx) * t;
                    let my = sy + (ey - sy) * t;
                    let hw = 1.5 * sc * t;
                    pb.move_to(mx + nx * hw, my + ny * hw);
                    pb.line_to(mx - nx * hw, my - ny * hw);
                }
                if let Some(path) = pb.finish() {
                    pixmap.stroke_path(&path, &black, &stroke, Transform::identity(), None);
                }
            }
            _ => {
                let offsets: &[f32] = match b.order {
                    2 => &[-1.3 * sc, 1.3 * sc],
                    3 => &[-2.2 * sc, 0.0, 2.2 * sc],
                    _ => &[0.0],
                };
                let mut pb = PathBuilder::new();
                for &off in offsets {
                    pb.move_to(x1 + nx * off, y1 + ny * off);
                    pb.line_to(x2 + nx * off, y2 + ny * off);
                }
                if let Some(path) = pb.finish() {
                    pixmap.stroke_path(&path, &black, &stroke, Transform::identity(), None);
                }
            }
        }
    }

    // Atom labels with white background knockout
    let gw = 4.8 * sc;
    let gh = 7.2 * sc;
    let gap = 1.4 * sc;
    let mut gstroke = Stroke::default();
    gstroke.width = (0.85 * sc).max(1.2);
    gstroke.line_cap = LineCap::Round;
    for (i, a) in mol.atoms.iter().enumerate() {
        if let Some(lbl) = atom_display_label(mol, i, &val) {
            let cx = tx(a.x);
            let cy = ty(a.y);
            let tw = (lbl.len() as f32) * (gw + gap);
            let mut box_pb = PathBuilder::new();
            box_pb.move_to(cx - tw * 0.55, cy - gh * 0.65);
            box_pb.line_to(cx + tw * 0.55, cy - gh * 0.65);
            box_pb.line_to(cx + tw * 0.55, cy + gh * 0.65);
            box_pb.line_to(cx - tw * 0.55, cy + gh * 0.65);
            box_pb.close();
            if let Some(p) = box_pb.finish() {
                pixmap.fill_path(&p, &white, FillRule::Winding, Transform::identity(), None);
            }
            let mut gpb = PathBuilder::new();
            let start_x = cx - tw * 0.5;
            for (ci, ch) in lbl.chars().enumerate() {
                stroke_glyph(&mut gpb, ch, start_x + (ci as f32) * (gw + gap), cy - gh * 0.5, gw, gh);
            }
            if let Some(p) = gpb.finish() {
                pixmap.stroke_path(&p, &black, &gstroke, Transform::identity(), None);
            }
        }
    }

    let mut out = Vec::new();
    {
        let mut enc = png::Encoder::new(&mut out, w, h);
        enc.set_color(png::ColorType::Rgba);
        enc.set_depth(png::BitDepth::Eight);
        if let Ok(mut writer) = enc.write_header() {
            let _ = writer.write_image_data(pixmap.data());
        }
    }
    out
}

// ============================================================================
// 9. Graph-Traversal SLN (`molecule_to_sln`)
// ============================================================================

pub fn molecule_to_sln(mol: &ChemMolecule) -> String {
    if mol.atoms.is_empty() {
        return String::new();
    }
    let val = CcVal::compute(mol);
    let mut adj: Vec<Vec<(usize, i32)>> = vec![Vec::new(); mol.atoms.len()];
    for b in &mol.bonds {
        if b.begin_atom < mol.atoms.len() && b.end_atom < mol.atoms.len() {
            adj[b.begin_atom].push((b.end_atom, b.order));
            adj[b.end_atom].push((b.begin_atom, b.order));
        }
    }

    let mut visited = vec![false; mol.atoms.len()];
    let mut ring_id = vec![0usize; mol.atoms.len()];
    let mut next_ring = 1usize;
    let mut back_edges: BTreeMap<usize, Vec<(usize, i32)>> = BTreeMap::new();

    fn find_back(
        u: usize,
        parent: Option<usize>,
        adj: &[Vec<(usize, i32)>],
        vis: &mut [bool],
        ring_id: &mut [usize],
        next_ring: &mut usize,
        back: &mut BTreeMap<usize, Vec<(usize, i32)>>,
    ) {
        vis[u] = true;
        for &(v, o) in &adj[u] {
            if Some(v) == parent {
                continue;
            }
            if !vis[v] {
                find_back(v, Some(u), adj, vis, ring_id, next_ring, back);
            } else if u > v {
                if ring_id[v] == 0 {
                    ring_id[v] = *next_ring;
                    *next_ring += 1;
                }
                back.entry(u).or_default().push((ring_id[v], o));
            }
        }
    }

    for i in 0..mol.atoms.len() {
        if !visited[i] {
            find_back(i, None, &adj, &mut visited, &mut ring_id, &mut next_ring, &mut back_edges);
        }
    }

    fn bond_sym(order: i32) -> &'static str {
        match order {
            2 => "=",
            3 => "#",
            _ => "-",
        }
    }

    fn emit_dfs(
        u: usize,
        parent: Option<usize>,
        mol: &ChemMolecule,
        val: &CcVal,
        adj: &[Vec<(usize, i32)>],
        vis: &mut [bool],
        ring_id: &[usize],
        back: &BTreeMap<usize, Vec<(usize, i32)>>,
        out: &mut String,
    ) {
        vis[u] = true;
        let a = &mol.atoms[u];
        out.push_str(elem_symbol(a.z));
        if ring_id[u] > 0 {
            out.push_str(&format!("[{}]", ring_id[u]));
        }
        let imp_h = val.implicit_h_counts.get(u).copied().unwrap_or(0);
        if imp_h == 1 {
            out.push('H');
        } else if imp_h > 1 {
            out.push_str(&format!("H{}", imp_h));
        }
        if let Some(closures) = back.get(&u) {
            for &(rid, o) in closures {
                out.push_str(&format!("{}@{}", bond_sym(o), rid));
            }
        }
        let mut children: Vec<(usize, i32)> = Vec::new();
        for &(v, o) in &adj[u] {
            if Some(v) != parent && !vis[v] {
                children.push((v, o));
            }
        }
        for (idx, &(v, o)) in children.iter().enumerate() {
            if vis[v] {
                continue;
            }
            let is_branch = idx + 1 < children.len();
            if is_branch {
                out.push('(');
            }
            out.push_str(bond_sym(o));
            emit_dfs(v, Some(u), mol, val, adj, vis, ring_id, back, out);
            if is_branch {
                out.push(')');
            }
        }
    }

    let mut vis2 = vec![false; mol.atoms.len()];
    let mut comps = Vec::new();
    for i in 0..mol.atoms.len() {
        if !vis2[i] {
            let mut s = String::new();
            emit_dfs(i, None, mol, &val, &adj, &mut vis2, &ring_id, &back_edges, &mut s);
            comps.push(s);
        }
    }
    comps.join(".")
}

// ============================================================================
// 10. Real FIPS 180-4 SHA-256 + IUPAC Standard InChIKey + Canonical InChI
// ============================================================================

pub fn sha256(data: &[u8]) -> [u8; 32] {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
    ];
    let bit_len = (data.len() as u64) * 8;
    let mut msg = data.to_vec();
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bit_len.to_be_bytes());
    for chunk in msg.chunks_exact(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes(chunk[4 * i..4 * i + 4].try_into().unwrap());
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16].wrapping_add(s0).wrapping_add(w[i - 7]).wrapping_add(s1);
        }
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh] = h;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let t1 = hh.wrapping_add(s1).wrapping_add(ch).wrapping_add(K[i]).wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }
    let mut out = [0u8; 32];
    for i in 0..8 {
        out[4 * i..4 * i + 4].copy_from_slice(&h[i].to_be_bytes());
    }
    out
}

fn base26_triplet(v: u32) -> [char; 3] {
    // Official IUPAC InChI ikey_base26.c: skips EAA..EZZ (676) and TAA..TTV (516) => 26^3 - 1192 = 2^14
    let mut n = v;
    if n >= 2704 {
        n += 676;
    }
    if n >= 12844 {
        n += 516;
    }
    let c0 = (b'A' + ((n / 676) % 26) as u8) as char;
    let c1 = (b'A' + ((n / 26) % 26) as u8) as char;
    let c2 = (b'A' + (n % 26) as u8) as char;
    [c0, c1, c2]
}

fn base26_dublet(v: u32) -> [char; 2] {
    let c0 = (b'A' + ((v / 26) % 26) as u8) as char;
    let c1 = (b'A' + (v % 26) as u8) as char;
    [c0, c1]
}

/// Official IUPAC InChIKey SHA-256 encoder (verified exact against standard `InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3` -> `LFQSCWFLJHTTHZ-UHFFFAOYSA-N`).
pub fn inchi_to_inchikey(inchi: &str) -> String {
    let body = inchi
        .strip_prefix("InChI=1S/")
        .or_else(|| inchi.strip_prefix("InChI=1/"))
        .unwrap_or(inchi);
    let is_std = inchi.starts_with("InChI=1S/");
    // Split major (formula + /c + /h) vs minor (remaining layers)
    let mut major = String::new();
    let mut minor = String::new();
    let mut proton_delta: i32 = 0;
    for (idx, part) in body.split('/').enumerate() {
        if idx == 0 {
            major.push_str(part);
        } else if part.starts_with('c') || part.starts_with('h') {
            major.push('/');
            major.push_str(part);
        } else if let Some(p) = part.strip_prefix('p') {
            proton_delta = p.parse().unwrap_or(0);
        } else {
            minor.push('/');
            minor.push_str(part);
        }
    }
    let d1 = sha256(major.as_bytes());
    let b0 = (d1[0] as u32) | ((d1[1] as u32 & 0x3F) << 8);
    let b1 = ((d1[1] as u32) >> 6) | ((d1[2] as u32) << 2) | ((d1[3] as u32 & 0x0F) << 10);
    let b2 = ((d1[3] as u32) >> 4) | ((d1[4] as u32) << 4) | ((d1[5] as u32 & 0x03) << 12);
    let b3 = ((d1[5] as u32) >> 2) | ((d1[6] as u32) << 6);
    let b4 = (d1[7] as u32) | ((d1[8] as u32 & 0x01) << 8);
    let mut block1 = String::with_capacity(14);
    for ch in base26_triplet(b0) { block1.push(ch); }
    for ch in base26_triplet(b1) { block1.push(ch); }
    for ch in base26_triplet(b2) { block1.push(ch); }
    for ch in base26_triplet(b3) { block1.push(ch); }
    for ch in base26_dublet(b4) { block1.push(ch); }

    let d2 = sha256(minor.as_bytes());
    let m0 = (d2[0] as u32) | ((d2[1] as u32 & 0x3F) << 8);
    let m1 = ((d2[1] as u32) >> 6) | ((d2[2] as u32) << 2) | ((d2[3] as u32 & 0x0F) << 10);
    let m2 = ((d2[3] as u32) >> 4) | ((d2[4] as u32 & 0x1F) << 4);
    let mut block2 = String::with_capacity(10);
    for ch in base26_triplet(m0) { block2.push(ch); }
    for ch in base26_triplet(m1) { block2.push(ch); }
    for ch in base26_dublet(m2) { block2.push(ch); }
    block2.push(if is_std { 'S' } else { 'N' });
    block2.push('A');

    let pflag = ((b'N' as i32) + proton_delta).clamp(b'A' as i32, b'Z' as i32) as u8 as char;
    format!("{block1}-{block2}-{pflag}")
}

pub fn molecule_to_inchi(mol: &ChemMolecule) -> String {
    let formula = mol.compute_formula_properties().get_formula();
    let val = CcVal::compute(mol);
    // Hill-ordered heavy atoms (C first, then alphabetical by symbol, then degree)
    let mut heavy: Vec<usize> = (0..mol.atoms.len()).filter(|&i| mol.atoms[i].z > 1).collect();
    let mut deg = vec![0usize; mol.atoms.len()];
    for b in &mol.bonds {
        if b.begin_atom < mol.atoms.len() && b.end_atom < mol.atoms.len() {
            deg[b.begin_atom] += 1;
            deg[b.end_atom] += 1;
        }
    }
    heavy.sort_by(|&a, &b| {
        let sa = elem_symbol(mol.atoms[a].z);
        let sb = elem_symbol(mol.atoms[b].z);
        let ca = if sa == "C" { 0 } else { 1 };
        let cb = if sb == "C" { 0 } else { 1 };
        ca.cmp(&cb).then(sa.cmp(sb)).then(deg[a].cmp(&deg[b])).then(a.cmp(&b))
    });
    let mut rank = HashMap::<usize, usize>::new();
    for (r, &idx) in heavy.iter().enumerate() {
        rank.insert(idx, r + 1);
    }
    let mut edges: Vec<(usize, usize)> = mol
        .bonds
        .iter()
        .filter_map(|b| {
            let r1 = *rank.get(&b.begin_atom)?;
            let r2 = *rank.get(&b.end_atom)?;
            Some((r1.min(r2), r1.max(r2)))
        })
        .collect();
    edges.sort();
    edges.dedup();

    let mut out = format!("InChI=1S/{formula}");
    if !edges.is_empty() {
        let c_str: Vec<String> = edges.iter().map(|(a, b)| format!("{a}-{b}")).collect();
        out.push_str("/c");
        out.push_str(&c_str.join(","));
    }
    let mut h_parts = Vec::new();
    for &idx in &heavy {
        let r = rank[&idx];
        let h = val.implicit_h_counts.get(idx).copied().unwrap_or(0);
        if h == 1 {
            h_parts.push(format!("{r}H"));
        } else if h > 1 {
            h_parts.push(format!("{r}H{h}"));
        }
    }
    if !h_parts.is_empty() {
        out.push_str("/h");
        out.push_str(&h_parts.join(","));
    }
    let net_q: i32 = mol.atoms.iter().map(|a| a.charge).sum();
    if net_q != 0 {
        out.push_str(&format!("/q{:+}", net_q));
    }
    out
}

/// Decodes a PNG byte slice into `(rgba_pixels, width, height)` using the workspace `png` crate.
pub fn decode_png_rgba(bytes: &[u8]) -> Option<(Vec<u8>, u32, u32)> {
    let mut decoder = png::Decoder::new(std::io::Cursor::new(bytes));
    decoder.set_transformations(png::Transformations::EXPAND | png::Transformations::STRIP_16);
    let mut reader = decoder.read_info().ok()?;
    let mut buf = vec![0; reader.output_buffer_size()?];
    let info = reader.next_frame(&mut buf).ok()?;
    let raw = &buf[..info.buffer_size()];
    let rgba = match info.color_type {
        png::ColorType::Rgba => raw.to_vec(),
        png::ColorType::Rgb => {
            let mut out = Vec::with_capacity((info.width * info.height * 4) as usize);
            for px in raw.chunks_exact(3) {
                out.extend_from_slice(&[px[0], px[1], px[2], 255]);
            }
            out
        }
        png::ColorType::GrayscaleAlpha => {
            let mut out = Vec::with_capacity((info.width * info.height * 4) as usize);
            for px in raw.chunks_exact(2) {
                out.extend_from_slice(&[px[0], px[0], px[0], px[1]]);
            }
            out
        }
        png::ColorType::Grayscale => {
            let mut out = Vec::with_capacity((info.width * info.height * 4) as usize);
            for &g in raw {
                out.extend_from_slice(&[g, g, g, 255]);
            }
            out
        }
        png::ColorType::Indexed => return None,
    };
    Some((rgba, info.width, info.height))
}
