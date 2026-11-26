//! Markush R-Group Library Enumerator & SAR Table Generator,
//! Offline Retrosynthetic Disconnection Advisor,
//! and Publication TOC Graphic / Reaction Scheme Composer.

use std::collections::{BTreeMap, VecDeque};

use chiralcraft_chem::{
    compute_analysis_report, compute_medchem_report, BondDisplay, CcAromaticityInfo, ChemAtom,
    ChemMolecule, MedChemReport,
};
use chiralcraft_formats::{molecule_to_smiles, molecules_to_sdf, smiles_to_molecule, SdfRecord};

use crate::{
    clean_up_structure_preserving_stereo, CaptionNode, DocumentScene, LabelJustification,
    SceneArrow, SceneGraphic,
};

// ============================================================================
// 1. Markush / R-Group Library Enumerator & SAR Grid
// ============================================================================

#[derive(Debug, Clone)]
pub struct EnumeratedAnalog {
    pub id: String,
    pub r1_label: String,
    pub r2_label: String,
    pub smiles: String,
    pub formula: String,
    pub medchem: MedChemReport,
    pub molecule: ChemMolecule,
}

fn substituent_to_smiles(tok: &str) -> Option<String> {
    let t = tok.trim();
    if t.is_empty() {
        return None;
    }
    let mapped = match t {
        "H" => return Some(String::new()),
        "Me" | "CH3" => "C",
        "Et" | "CH2CH3" => "CC",
        "Pr" | "nPr" => "CCC",
        "iPr" => "C(C)C",
        "Bu" | "nBu" => "CCCC",
        "tBu" => "C(C)(C)C",
        "F" => "F",
        "Cl" => "Cl",
        "Br" => "Br",
        "I" => "I",
        "OH" => "O",
        "OMe" | "MeO" => "OC",
        "OEt" | "EtO" => "OCC",
        "CF3" | "F3C" => "C(F)(F)F",
        "OCF3" => "OC(F)(F)F",
        "CN" | "NC" => "C#N",
        "NO2" | "O2N" => "[N+](=O)[O-]",
        "NH2" | "H2N" => "N",
        "NMe2" => "N(C)C",
        "NHAc" => "NC(=O)C",
        "COOH" | "CO2H" => "C(=O)O",
        "CO2Me" | "COOMe" => "C(=O)OC",
        "CONH2" => "C(=O)N",
        "Ph" => "c1ccccc1",
        "Bn" => "Cc1ccccc1",
        "Py" | "Pyridyl" => "c1ccncc1",
        "Cy" => "C1CCCCC1",
        "Morph" | "Morpholino" => "N1CCOCC1",
        "Pip" | "Piperidino" => "N1CCCCC1",
        other => other,
    };
    Some(mapped.to_string())
}

fn replace_r_site(mol: &mut ChemMolecule, r_idx: usize, sub_smiles: &str) {
    if r_idx >= mol.atoms.len() {
        return;
    }
    if sub_smiles.is_empty() {
        // Replace R with implicit H: remove R atom and its bonds
        mol.bonds.retain(|b| b.begin_atom != r_idx && b.end_atom != r_idx);
        mol.atoms.remove(r_idx);
        for b in &mut mol.bonds {
            if b.begin_atom > r_idx { b.begin_atom -= 1; }
            if b.end_atom > r_idx { b.end_atom -= 1; }
        }
        return;
    }
    let Ok(sub) = smiles_to_molecule(sub_smiles) else {
        return;
    };
    if sub.atoms.is_empty() {
        return;
    }
    // Replace R atom with atom 0 of substituent, then append remaining substituent atoms & bonds
    mol.atoms[r_idx].z = sub.atoms[0].z;
    mol.atoms[r_idx].charge = sub.atoms[0].charge;
    mol.atoms[r_idx].label = None;

    let mut map = vec![0usize; sub.atoms.len()];
    map[0] = r_idx;
    let mut next_aid = mol.atoms.iter().map(|a| a.id).max().unwrap_or(1) + 1;
    let mut next_bid = mol.bonds.iter().map(|b| b.id).max().unwrap_or(100) + 1;
    for si in 1..sub.atoms.len() {
        let mut a = sub.atoms[si].clone();
        a.id = next_aid;
        next_aid += 1;
        a.x = mol.atoms[r_idx].x + (si as f64) * 8.0;
        a.y = mol.atoms[r_idx].y + ((si % 2) as f64) * 6.0;
        map[si] = mol.atoms.len();
        mol.add_atom(a);
    }
    for b in &sub.bonds {
        if b.begin_atom < sub.atoms.len() && b.end_atom < sub.atoms.len() {
            mol.add_bond(next_bid, map[b.begin_atom], map[b.end_atom], b.order, b.display);
            next_bid += 1;
        }
    }
}

fn find_r_atoms(mol: &ChemMolecule, tag: &str) -> Vec<usize> {
    let tag_u = tag.to_ascii_uppercase();
    mol.atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| {
            a.label
                .as_deref()
                .map(|l| {
                    let lu = l.trim().to_ascii_uppercase();
                    lu == tag_u || (tag_u == "R1" && lu == "R")
                })
                .unwrap_or(false)
        })
        .map(|(i, _)| i)
        .collect()
}

/// Enumerates a Markush / R-group library from `scaffold` containing `R1` (or `R`) and optional `R2` labels.
/// If no `R1` atom exists in `scaffold`, attaches `R1` to the last atom automatically so users can enumerate any drawn core!
pub fn enumerate_rgroup_library(
    scaffold: &ChemMolecule,
    r1_csv: &str,
    r2_csv: &str,
    bond_len_pt: f64,
) -> Vec<EnumeratedAnalog> {
    let r1_list: Vec<String> = r1_csv
        .split(|c: char| c == ',' || c == ';' || c == '\n')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    let r2_list_raw: Vec<String> = r2_csv
        .split(|c: char| c == ',' || c == ';' || c == '\n')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if scaffold.atoms.is_empty() || r1_list.is_empty() {
        return Vec::new();
    }

    let mut base = scaffold.clone();
    let has_r1 = !find_r_atoms(&base, "R1").is_empty();
    let has_r2 = !find_r_atoms(&base, "R2").is_empty();
    if !has_r1 {
        // Add an explicit R1 attachment stub on atom 0 so enumeration works even without pre-labeling
        let aid = base.atoms.iter().map(|a| a.id).max().unwrap_or(1) + 1;
        let bid = base.bonds.iter().map(|b| b.id).max().unwrap_or(100) + 1;
        let mut r_atom = ChemAtom::new(aid, 6);
        r_atom.x = base.atoms[0].x + bond_len_pt;
        r_atom.y = base.atoms[0].y;
        r_atom.label = Some("R1".to_string());
        let ri = base.atoms.len();
        base.add_atom(r_atom);
        base.add_bond(bid, 0, ri, 1, BondDisplay::Solid);
    }

    let r2_iter: Vec<String> = if has_r2 && !r2_list_raw.is_empty() {
        r2_list_raw
    } else {
        vec!["—".to_string()]
    };

    let mut out = Vec::new();
    let mut seq = 1usize;
    for r1 in &r1_list {
        let Some(r1_smi) = substituent_to_smiles(r1) else { continue };
        for r2 in &r2_iter {
            let mut m = base.clone();
            // Replace R2 sites in descending index order first, then R1 sites
            if has_r2 && r2 != "—" {
                if let Some(r2_smi) = substituent_to_smiles(r2) {
                    let mut idxs = find_r_atoms(&m, "R2");
                    idxs.sort_unstable_by(|a, b| b.cmp(a));
                    for idx in idxs {
                        replace_r_site(&mut m, idx, &r2_smi);
                    }
                }
            }
            let mut r1_idxs = find_r_atoms(&m, "R1");
            r1_idxs.sort_unstable_by(|a, b| b.cmp(a));
            for idx in r1_idxs {
                replace_r_site(&mut m, idx, &r1_smi);
            }
            clean_up_structure_preserving_stereo(&mut m, &[], bond_len_pt);
            let rep = compute_analysis_report(&m);
            let medchem = compute_medchem_report(&m);
            let smiles = molecule_to_smiles(&m);
            out.push(EnumeratedAnalog {
                id: format!("AN-{seq:03}"),
                r1_label: r1.clone(),
                r2_label: r2.clone(),
                smiles,
                formula: rep.formula,
                medchem,
                molecule: m,
            });
            seq += 1;
            if out.len() >= 64 {
                return out;
            }
        }
    }
    out
}

pub fn export_sar_csv(analogs: &[EnumeratedAnalog]) -> String {
    let mut s = String::from("ID,R1,R2,Formula,SMILES,MW,LogP,tPSA,HBD,HBA,Fsp3,CNS_MPO,QED,SA,Lipinski,PAINS_Alerts\n");
    for a in analogs {
        s.push_str(&format!(
            "{},{},{},{},\"{}\",{:.2},{:.2},{:.2},{},{},{:.2},{:.2},{:.2},{:.2},{},{}\n",
            a.id,
            a.r1_label,
            a.r2_label,
            a.formula,
            a.smiles,
            a.medchem.mw,
            a.medchem.log_p,
            a.medchem.tpsa,
            a.medchem.hbd,
            a.medchem.hba,
            a.medchem.fsp3,
            a.medchem.cns_mpo,
            a.medchem.qed,
            a.medchem.sa_score,
            if a.medchem.lipinski_pass { "PASS" } else { "FAIL" },
            a.medchem.alerts.len()
        ));
    }
    s
}

pub fn export_sar_sdf(analogs: &[EnumeratedAnalog]) -> String {
    let recs: Vec<SdfRecord> = analogs
        .iter()
        .map(|a| {
            let mut props = BTreeMap::new();
            props.insert("R1".into(), a.r1_label.clone());
            props.insert("R2".into(), a.r2_label.clone());
            props.insert("MW".into(), format!("{:.2}", a.medchem.mw));
            props.insert("LogP".into(), format!("{:.2}", a.medchem.log_p));
            props.insert("tPSA".into(), format!("{:.2}", a.medchem.tpsa));
            props.insert("CNS_MPO".into(), format!("{:.2}", a.medchem.cns_mpo));
            props.insert("QED".into(), format!("{:.2}", a.medchem.qed));
            props.insert("PAINS_Alerts".into(), a.medchem.alerts.len().to_string());
            SdfRecord {
                title: a.id.clone(),
                molecule: a.molecule.clone(),
                properties: props,
            }
        })
        .collect();
    molecules_to_sdf(&recs)
}

// ============================================================================
// 2. Offline Retrosynthetic Disconnection Advisor
// ============================================================================

#[derive(Debug, Clone)]
pub struct RetrosyntheticHit {
    pub reaction_name: &'static str,
    pub bond_idx: usize,
    pub atom_u: usize,
    pub atom_v: usize,
    pub reagents: &'static str,
    pub score: u8,
    pub synthon_a: ChemMolecule,
    pub synthon_b: ChemMolecule,
    pub synthons_summary: String,
}

fn split_molecule_at_bond(
    mol: &ChemMolecule,
    cut_bi: usize,
    cap_u_z: Option<i32>,
    cap_v_z: Option<i32>,
    bond_len: f64,
) -> Option<(ChemMolecule, ChemMolecule)> {
    let n = mol.atoms.len();
    let bcut = mol.bonds.get(cut_bi)?;
    let (u0, v0) = (bcut.begin_atom, bcut.end_atom);

    // BFS from u0 without traversing cut_bi; if v0 is reachable, bond is in a ring (skip simple 2-piece split)
    let mut adj = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        if bi != cut_bi && b.begin_atom < n && b.end_atom < n {
            adj[b.begin_atom].push(b.end_atom);
            adj[b.end_atom].push(b.begin_atom);
        }
    }
    let mut side_a = vec![false; n];
    let mut q = VecDeque::new();
    side_a[u0] = true;
    q.push_back(u0);
    while let Some(cur) = q.pop_front() {
        for &nb in &adj[cur] {
            if !side_a[nb] {
                side_a[nb] = true;
                q.push_back(nb);
            }
        }
    }
    if side_a[v0] {
        return None; // Ring bond; keep acyclic strategic disconnections crisp
    }

    let build_half = |want_a: bool, attach_orig: usize, cap_z: Option<i32>| -> ChemMolecule {
        let mut out = ChemMolecule::new();
        let mut map = vec![usize::MAX; n];
        for (i, a) in mol.atoms.iter().enumerate() {
            if side_a[i] == want_a {
                map[i] = out.atoms.len();
                out.add_atom(a.clone());
            }
        }
        let mut next_bid = 500u32;
        for (bi, b) in mol.bonds.iter().enumerate() {
            if bi == cut_bi || b.begin_atom >= n || b.end_atom >= n {
                continue;
            }
            if side_a[b.begin_atom] == want_a && side_a[b.end_atom] == want_a {
                out.add_bond(next_bid, map[b.begin_atom], map[b.end_atom], b.order, b.display);
                next_bid += 1;
            }
        }
        if let Some(cz) = cap_z {
            let att = map[attach_orig];
            if att < out.atoms.len() {
                let nid = out.atoms.iter().map(|a| a.id).max().unwrap_or(1) + 1;
                let mut ca = ChemAtom::new(nid, cz);
                ca.x = out.atoms[att].x + bond_len;
                ca.y = out.atoms[att].y;
                if cz == 5 {
                    ca.label = Some("B(OH)₂".to_string());
                }
                let ci = out.atoms.len();
                out.add_atom(ca);
                out.add_bond(next_bid, att, ci, 1, BondDisplay::Solid);
            }
        }
        clean_up_structure_preserving_stereo(&mut out, &[], bond_len);
        out
    };

    let ma = build_half(true, u0, cap_u_z);
    let mb = build_half(false, v0, cap_v_z);
    if ma.atoms.is_empty() || mb.atoms.is_empty() {
        return None;
    }
    Some((ma, mb))
}

/// Suggests classic strategic retrosynthetic disconnections for `mol` (Amide, Ester, Suzuki,
/// Buchwald-Hartwig, Sulfonamide, Williamson Ether, Reductive Amination, Friedel-Crafts, Wittig).
pub fn suggest_disconnections(mol: &ChemMolecule, bond_len: f64) -> Vec<RetrosyntheticHit> {
    let n = mol.atoms.len();
    let arom = CcAromaticityInfo::perceive(mol);
    let mut adj = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.begin_atom < n && b.end_atom < n {
            adj[b.begin_atom].push((b.end_atom, bi, b.order));
            adj[b.end_atom].push((b.begin_atom, bi, b.order));
        }
    }
    let has_carbonyl = |idx: usize| -> bool {
        mol.atoms[idx].z == 6 && adj[idx].iter().any(|&(k, _, o)| mol.atoms[k].z == 8 && o == 2)
    };
    let has_sulfonyl = |idx: usize| -> bool {
        mol.atoms[idx].z == 16
            && adj[idx].iter().filter(|&&(k, _, o)| mol.atoms[k].z == 8 && o == 2).count() >= 2
    };

    let mut hits = Vec::new();
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.begin_atom >= n || b.end_atom >= n {
            continue;
        }
        let (u, v) = (b.begin_atom, b.end_atom);
        let (zu, zv) = (mol.atoms[u].z, mol.atoms[v].z);

        let rule: Option<(&'static str, &'static str, u8, Option<i32>, Option<i32>)> = if b.order == 1 {
            if has_carbonyl(u) && zv == 7 {
                Some(("Amide Coupling (HATU / EDC)", "HATU, DIPEA, DMF, rt", 98, Some(8), None))
            } else if has_carbonyl(v) && zu == 7 {
                Some(("Amide Coupling (HATU / EDC)", "HATU, DIPEA, DMF, rt", 98, None, Some(8)))
            } else if has_carbonyl(u) && zv == 8 {
                Some(("Esterification (EDC / DMAP)", "EDC·HCl, DMAP, DCM, rt", 94, Some(8), None))
            } else if has_carbonyl(v) && zu == 8 {
                Some(("Esterification (EDC / DMAP)", "EDC·HCl, DMAP, DCM, rt", 94, None, Some(8)))
            } else if has_sulfonyl(u) && zv == 7 {
                Some(("Sulfonamide Coupling", "Et3N or Pyridine, DCM, 0 °C → rt", 95, Some(17), None))
            } else if has_sulfonyl(v) && zu == 7 {
                Some(("Sulfonamide Coupling", "Et3N or Pyridine, DCM, 0 °C → rt", 95, None, Some(17)))
            } else if zu == 6 && zv == 6 && arom.is_atom_aromatic(u) && arom.is_atom_aromatic(v) {
                Some(("Suzuki–Miyaura Biaryl Cross-Coupling", "Pd(PPh3)4 (5 mol%), K2CO3, THF/H2O, 80 °C", 96, Some(35), Some(5)))
            } else if zu == 6 && arom.is_atom_aromatic(u) && zv == 7 && !has_carbonyl(u) {
                Some(("Buchwald–Hartwig C–N Amination", "Pd2(dba)3, BrettPhos, NaOtBu, dioxane, 100 °C", 90, Some(35), None))
            } else if zv == 6 && arom.is_atom_aromatic(v) && zu == 7 && !has_carbonyl(v) {
                Some(("Buchwald–Hartwig C–N Amination", "Pd2(dba)3, BrettPhos, NaOtBu, dioxane, 100 °C", 90, None, Some(35)))
            } else if zu == 8 && zv == 6 && !has_carbonyl(v) && !arom.is_atom_aromatic(v) {
                Some(("Williamson Ether Synthesis", "NaH or K2CO3, DMF, 60 °C", 85, None, Some(35)))
            } else if zv == 8 && zu == 6 && !has_carbonyl(u) && !arom.is_atom_aromatic(u) {
                Some(("Williamson Ether Synthesis", "NaH or K2CO3, DMF, 60 °C", 85, Some(35), None))
            } else if zu == 7 && zv == 6 && !has_carbonyl(v) && !arom.is_atom_aromatic(v) {
                Some(("Reductive Amination", "NaBH(OAc)3, AcOH, DCE, rt", 82, None, Some(8)))
            } else if zv == 7 && zu == 6 && !has_carbonyl(u) && !arom.is_atom_aromatic(u) {
                Some(("Reductive Amination", "NaBH(OAc)3, AcOH, DCE, rt", 82, Some(8), None))
            } else if has_carbonyl(u) && zv == 6 && arom.is_atom_aromatic(v) {
                Some(("Friedel–Crafts Acylation", "AlCl3, DCM, 0 °C → rt", 80, Some(17), None))
            } else if has_carbonyl(v) && zu == 6 && arom.is_atom_aromatic(u) {
                Some(("Friedel–Crafts Acylation", "AlCl3, DCM, 0 °C → rt", 80, None, Some(17)))
            } else {
                None
            }
        } else if b.order == 2 && zu == 6 && zv == 6 && !arom.is_atom_aromatic(u) && !arom.is_atom_aromatic(v) {
            Some(("Wittig Olefination / Cross-Metathesis", "Grubbs II (3 mol%), DCM, 40 °C  or  Ph3P=CHR", 78, Some(8), Some(8)))
        } else {
            None
        };

        if let Some((name, reagents, score, cap_u, cap_v)) = rule {
            if let Some((sa, sb)) = split_molecule_at_bond(mol, bi, cap_u, cap_v, bond_len) {
                let fa = compute_analysis_report(&sa).formula;
                let fb = compute_analysis_report(&sb).formula;
                hits.push(RetrosyntheticHit {
                    reaction_name: name,
                    bond_idx: bi,
                    atom_u: u,
                    atom_v: v,
                    reagents,
                    score,
                    synthon_a: sa,
                    synthon_b: sb,
                    synthons_summary: format!("{fa}  +  {fb}"),
                });
            }
        }
    }
    hits.sort_by(|a, b| b.score.cmp(&a.score));
    hits
}

/// Stamps `Target  ⇒  Synthon A  +  Synthon B` with suggested forward conditions onto `scene`.
pub fn stamp_retrosynthesis_scheme(
    scene: &mut DocumentScene,
    next_id: &mut u32,
    target: &ChemMolecule,
    hit: &RetrosyntheticHit,
    origin_x: f64,
    origin_y: f64,
) {
    let shift_center = |m: &mut ChemMolecule, cx: f64, cy: f64| {
        if m.atoms.is_empty() { return; }
        let (mut x1, mut y1, mut x2, mut y2) = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for a in &m.atoms {
            x1 = x1.min(a.x); y1 = y1.min(a.y); x2 = x2.max(a.x); y2 = y2.max(a.y);
        }
        let dx = cx - 0.5 * (x1 + x2);
        let dy = cy - 0.5 * (y1 + y2);
        for a in &mut m.atoms {
            a.x += dx;
            a.y += dy;
        }
    };

    let mut t = target.clone();
    let mut sa = hit.synthon_a.clone();
    let mut sb = hit.synthon_b.clone();
    shift_center(&mut t, origin_x, origin_y);
    shift_center(&mut sa, origin_x + 150.0, origin_y);
    shift_center(&mut sb, origin_x + 255.0, origin_y);

    scene.molecules.push(t);
    scene.molecules.push(sa);
    scene.molecules.push(sb);

    let aid = *next_id; *next_id += 1;
    scene.arrows.push(SceneArrow {
        id: aid,
        tail: (origin_x + 52.0, origin_y),
        head: (origin_x + 102.0, origin_y),
        arrowhead_head: "Retrosynthetic".to_string(),
        arrowhead_tail: None,
        no_go: None,
        dipole: false,
    });
    let cid_plus = *next_id; *next_id += 1;
    scene.captions.push(CaptionNode {
        id: cid_plus,
        text: "+".to_string(),
        x: origin_x + 202.0,
        y: origin_y,
        element: None,
        generic_nickname: None,
        warning: None,
        justification: LabelJustification::Center,
    });
    let cid_top = *next_id; *next_id += 1;
    scene.captions.push(CaptionNode {
        id: cid_top,
        text: hit.reaction_name.to_string(),
        x: origin_x + 77.0,
        y: origin_y - 16.0,
        element: None,
        generic_nickname: None,
        warning: None,
        justification: LabelJustification::Center,
    });
    let cid_bot = *next_id; *next_id += 1;
    scene.captions.push(CaptionNode {
        id: cid_bot,
        text: hit.reagents.to_string(),
        x: origin_x + 77.0,
        y: origin_y + 18.0,
        element: None,
        generic_nickname: None,
        warning: None,
        justification: LabelJustification::Center,
    });
}

// ============================================================================
// 3. Publication TOC Graphic Presets, Arrow Condition Stamper & Serpentine Wrapper
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TocPreset {
    pub name: &'static str,
    pub width_cm: f64,
    pub height_cm: f64,
}

impl TocPreset {
    pub fn presets() -> &'static [TocPreset] {
        &[
            TocPreset { name: "ACS / JACS TOC (8.25 × 4.45 cm)", width_cm: 8.25, height_cm: 4.45 },
            TocPreset { name: "Nature Research (9.0 × 9.0 cm)", width_cm: 9.0, height_cm: 9.0 },
            TocPreset { name: "RSC TOC (8.0 × 4.0 cm)", width_cm: 8.0, height_cm: 4.0 },
            TocPreset { name: "Angewandte Chemie (11.0 × 4.0 cm)", width_cm: 11.0, height_cm: 4.0 },
            TocPreset { name: "Widescreen Slide 16:9 (16.0 × 9.0 cm)", width_cm: 16.0, height_cm: 9.0 },
        ]
    }

    pub fn size_pt(&self) -> (f64, f64) {
        (self.width_cm * 28.3465, self.height_cm * 28.3465)
    }
}

/// Inserts an exact journal-dimensioned TOC Graphic frame onto `scene` and returns its `(x1, y1, x2, y2)` bbox.
pub fn insert_toc_frame(
    scene: &mut DocumentScene,
    next_id: &mut u32,
    preset: TocPreset,
    origin_x: f64,
    origin_y: f64,
) -> (f64, f64, f64, f64) {
    let (w, h) = preset.size_pt();
    let bbox = (origin_x, origin_y, origin_x + w, origin_y + h);
    let gid = *next_id; *next_id += 1;
    scene.graphics.push(SceneGraphic {
        id: gid,
        graphic_type: "Rectangle".to_string(),
        sub_type: Some("Plain".to_string()),
        bbox,
    });
    let cid = *next_id; *next_id += 1;
    scene.captions.push(CaptionNode {
        id: cid,
        text: preset.name.to_string(),
        x: origin_x + 4.0,
        y: origin_y - 10.0,
        element: None,
        generic_nickname: None,
        warning: None,
        justification: LabelJustification::Left,
    });
    bbox
}

/// Returns how many molecules have atoms extending outside `toc_bbox`.
pub fn count_molecules_overflowing_frame(scene: &DocumentScene, bbox: (f64, f64, f64, f64)) -> usize {
    let (x1, y1, x2, y2) = bbox;
    scene
        .molecules
        .iter()
        .filter(|m| {
            m.atoms
                .iter()
                .any(|a| a.x < x1 || a.x > x2 || a.y < y1 || a.y > y2)
        })
        .count()
}

/// Stamps clean, centered above-arrow (reagents) and below-arrow (solvent, temp, time, yield) captions.
pub fn stamp_arrow_conditions(
    scene: &mut DocumentScene,
    next_id: &mut u32,
    arrow_idx: usize,
    above_text: &str,
    below_text: &str,
) -> bool {
    let Some(arr) = scene.arrows.get(arrow_idx).cloned() else {
        return false;
    };
    let mx = 0.5 * (arr.tail.0 + arr.head.0);
    let my = 0.5 * (arr.tail.1 + arr.head.1);
    if !above_text.trim().is_empty() {
        let cid = *next_id; *next_id += 1;
        scene.captions.push(CaptionNode {
            id: cid,
            text: above_text.trim().to_string(),
            x: mx,
            y: my - 13.0,
            element: None,
            generic_nickname: None,
            warning: None,
            justification: LabelJustification::Center,
        });
    }
    if !below_text.trim().is_empty() {
        let cid = *next_id; *next_id += 1;
        scene.captions.push(CaptionNode {
            id: cid,
            text: below_text.trim().to_string(),
            x: mx,
            y: my + 14.0,
            element: None,
            generic_nickname: None,
            warning: None,
            justification: LabelJustification::Center,
        });
    }
    true
}

/// Wraps all molecules in `scene` into a tidy multi-row publication reaction scheme (`per_row` steps per row)
/// with auto-generated connecting reaction arrows!
pub fn wrap_serpentine_reaction_scheme(
    scene: &mut DocumentScene,
    next_id: &mut u32,
    per_row: usize,
) -> usize {
    let n = scene.molecules.len();
    if n < 2 {
        return 0;
    }
    let cols = per_row.max(2);
    let col_step = 148.0;
    let row_step = 115.0;
    let start_x = 95.0;
    let start_y = 110.0;

    let mut centers = Vec::with_capacity(n);
    for (i, m) in scene.molecules.iter_mut().enumerate() {
        let row = i / cols;
        let col_in_row = i % cols;
        // Serpentine: even rows left-to-right, odd rows right-to-left!
        let col = if row % 2 == 0 { col_in_row } else { (cols - 1).saturating_sub(col_in_row) };
        let tx = start_x + (col as f64) * col_step;
        let ty = start_y + (row as f64) * row_step;
        if !m.atoms.is_empty() {
            let (mut x1, mut y1, mut x2, mut y2) = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
            for a in &m.atoms {
                x1 = x1.min(a.x); y1 = y1.min(a.y); x2 = x2.max(a.x); y2 = y2.max(a.y);
            }
            let dx = tx - 0.5 * (x1 + x2);
            let dy = ty - 0.5 * (y1 + y2);
            for a in &mut m.atoms {
                a.x += dx;
                a.y += dy;
            }
            let r = 0.5 * (x2 - x1).hypot(y2 - y1).max(24.0);
            centers.push((tx, ty, r));
        } else {
            centers.push((tx, ty, 20.0));
        }
    }

    scene.arrows.clear();
    for k in 0..(centers.len() - 1) {
        let (x1, y1, r1) = centers[k];
        let (x2, y2, r2) = centers[k + 1];
        let dvx = x2 - x1;
        let dvy = y2 - y1;
        let dist = dvx.hypot(dvy).max(1.0);
        let ux = dvx / dist;
        let uy = dvy / dist;
        let p1 = (r1 + 8.0).min(dist * 0.35);
        let p2 = (r2 + 8.0).min(dist * 0.35);
        let aid = *next_id; *next_id += 1;
        scene.arrows.push(SceneArrow {
            id: aid,
            tail: (x1 + ux * p1, y1 + uy * p1),
            head: (x2 - ux * p2, y2 - uy * p2),
            arrowhead_head: "Full".to_string(),
            arrowhead_tail: None,
            no_go: None,
            dipole: false,
        });
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_markush_retrosynthesis_and_toc_tools() {
        // 1. Markush R-group enumeration on benzoic scaffold
        let mut core = smiles_to_molecule("c1ccccc1C").unwrap();
        core.atoms[6].label = Some("R1".to_string());
        let analogs = enumerate_rgroup_library(&core, "H, Me, F, OMe, CF3", "", 14.40);
        assert_eq!(analogs.len(), 5);
        let csv = export_sar_csv(&analogs);
        assert!(csv.contains("CNS_MPO") && csv.contains("AN-005"));
        let sdf = export_sar_sdf(&analogs);
        assert!(sdf.contains("$$$$"));

        // 2. Retrosynthesis advisor detects Amide Coupling on acetanilide (CC(=O)Nc1ccccc1)
        let acetanilide = smiles_to_molecule("CC(=O)Nc1ccccc1").unwrap();
        let hits = suggest_disconnections(&acetanilide, 14.40);
        assert!(
            hits.iter().any(|h| h.reaction_name.contains("Amide")),
            "Expected Amide disconnection for acetanilide, got: {:?}",
            hits.iter().map(|h| h.reaction_name).collect::<Vec<_>>()
        );

        // 3. TOC Frame + Serpentine scheme wrapper
        let mut scene = DocumentScene::default();
        let mut nid = 100u32;
        let bbox = insert_toc_frame(&mut scene, &mut nid, TocPreset::presets()[0], 20.0, 20.0);
        assert_eq!(scene.graphics.len(), 1);
        assert!((bbox.2 - bbox.0 - 8.25 * 28.3465).abs() < 0.1);
        for _ in 0..4 {
            scene.molecules.push(smiles_to_molecule("c1ccccc1").unwrap());
        }
        assert_eq!(wrap_serpentine_reaction_scheme(&mut scene, &mut nid, 2), 4);
        assert_eq!(scene.arrows.len(), 3);
    }
}
