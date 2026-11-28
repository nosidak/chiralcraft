//! Complete 98-binding ChiralCraft Context-Sensitive Hotkey Engine
//! implementing atom, bond, and generic canvas hotkeys:
//! - Atom Context (57 bindings): 33 `LABELTEXT`, 17 `SPROUT`, `CHARGE +/-`, `ATOMNUMBER`, `ATTACHMENTPOINT`, `SELECT`
//! - Bond Context (27 bindings): `BONDORDER`, `BONDDISPLAY`, `BONDDOUBLEPOSITION`, `FUSERING`, `SELECT`
//! - Generic Context (14 bindings): `TOOLMODE` (11), `DIALOG` (3)

use chiralcraft_chem::{get_element_by_symbol, BondDisplay, ChemAtom, ChemBond, ChemMolecule};
use chiralcraft_doc::{stamp_template, BuiltinTemplate, InteractiveDocument};

use crate::icons::Tool;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DoubleBondAlign {
    Auto,
    Left,
    Right,
    Center,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HotkeyOutcome {
    Handled(String),
    SwitchTool(Tool),
    OpenChooseNickname,
    OpenPropertiesDialog,
    SelectMolecule(usize),
    SetDoubleAlign((usize, usize), DoubleBondAlign),
    ToggleAtomNumbers,
    NotMatched,
}

/// Computes the cleanest outward sprout angle from atom `ai` in `mol`.
pub fn outward_sprout_angle(mol: &ChemMolecule, ai: usize) -> f64 {
    let Some(a) = mol.atoms.get(ai) else {
        return -std::f64::consts::FRAC_PI_6;
    };
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut deg = 0usize;
    for b in &mol.bonds {
        let nb = if b.begin_atom == ai {
            Some(b.end_atom)
        } else if b.end_atom == ai {
            Some(b.begin_atom)
        } else {
            None
        };
        if let Some(ni) = nb.and_then(|i| mol.atoms.get(i)) {
            sx += ni.x - a.x;
            sy += ni.y - a.y;
            deg += 1;
        }
    }
    if deg == 0 {
        -std::f64::consts::FRAC_PI_6
    } else if deg == 1 {
        sy.atan2(sx) + 2.0 * std::f64::consts::FRAC_PI_3
    } else {
        (-sy).atan2(-sx)
    }
}

fn sprout_atom_bond(
    mol: &mut ChemMolecule,
    ai: usize,
    ang: f64,
    blen: f64,
    z: i32,
    order: i32,
    display: BondDisplay,
    next_id: &mut u32,
) -> usize {
    let (ax, ay) = (mol.atoms[ai].x, mol.atoms[ai].y);
    let nid = *next_id;
    *next_id += 1;
    let bid = *next_id;
    *next_id += 1;
    let mut na = ChemAtom::new(nid, z);
    na.x = ax + blen * ang.cos();
    na.y = ay + blen * ang.sin();
    let new_idx = mol.atoms.len();
    mol.atoms.push(na);
    mol.bonds.push(ChemBond {
        id: bid,
        begin_atom: ai,
        end_atom: new_idx,
        order,
        display,
        visible: true,
        ignore_warnings: false,
    });
    new_idx
}

fn attach_ring_to_atom(
    mol: &mut ChemMolecule,
    ai: usize,
    n_sides: usize,
    aromatic: bool,
    cpd: bool,
    pendant: bool,
    blen: f64,
    next_id: &mut u32,
) {
    let ang0 = outward_sprout_angle(mol, ai);
    let root = if pendant {
        sprout_atom_bond(mol, ai, ang0, blen, 6, 1, BondDisplay::Solid, next_id)
    } else {
        ai
    };
    let (ax, ay) = (mol.atoms[root].x, mol.atoms[root].y);
    let r = blen / (2.0 * (std::f64::consts::PI / (n_sides as f64)).sin());
    let rcx = ax + r * ang0.cos();
    let rcy = ay + r * ang0.sin();
    let start_ang = (ay - rcy).atan2(ax - rcx);
    let step = std::f64::consts::TAU / (n_sides as f64);
    let mut prev = root;
    for k in 1..n_sides {
        let ang = start_ang + (k as f64) * step;
        let nid = *next_id;
        *next_id += 1;
        let mut na = ChemAtom::new(nid, 6);
        na.x = rcx + r * ang.cos();
        na.y = rcy + r * ang.sin();
        let ni = mol.atoms.len();
        mol.atoms.push(na);
        let bid = *next_id;
        *next_id += 1;
        let ord = if aromatic && k % 2 == 1 {
            2
        } else if cpd && (k == 1 || k == 3) {
            2
        } else {
            1
        };
        mol.bonds.push(ChemBond {
            id: bid,
            begin_atom: prev,
            end_atom: ni,
            order: ord,
            display: BondDisplay::Solid,
            visible: true,
            ignore_warnings: false,
        });
        prev = ni;
    }
    let bid = *next_id;
    *next_id += 1;
    mol.bonds.push(ChemBond {
        id: bid,
        begin_atom: prev,
        end_atom: root,
        order: 1,
        display: BondDisplay::Solid,
        visible: true,
        ignore_warnings: false,
    });
}

fn fuse_ring_to_bond(
    mol: &mut ChemMolecule,
    bi: usize,
    n_sides: usize,
    aromatic: bool,
    cpd: bool,
    blen: f64,
    next_id: &mut u32,
) {
    let Some(b) = mol.bonds.get(bi).cloned() else {
        return;
    };
    if b.begin_atom >= mol.atoms.len() || b.end_atom >= mol.atoms.len() {
        return;
    }
    let p1 = (mol.atoms[b.begin_atom].x, mol.atoms[b.begin_atom].y);
    let p2 = (mol.atoms[b.end_atom].x, mol.atoms[b.end_atom].y);
    let mx = 0.5 * (p1.0 + p2.0);
    let my = 0.5 * (p1.1 + p2.1);
    let raw_len = (p2.0 - p1.0).hypot(p2.1 - p1.1);
    let edge_len = raw_len.max(blen);
    let cx = mol.atoms.iter().map(|a| a.x).sum::<f64>() / (mol.atoms.len().max(1) as f64);
    let cy = mol.atoms.iter().map(|a| a.y).sum::<f64>() / (mol.atoms.len().max(1) as f64);
    let (ex, ey) = if raw_len > 1e-6 {
        ((p2.0 - p1.0) / raw_len, (p2.1 - p1.1) / raw_len)
    } else {
        (1.0, 0.0)
    };
    let mut nx = -ey;
    let mut ny = ex;
    if (mx + nx - cx).hypot(my + ny - cy) < (mx - nx - cx).hypot(my - ny - cy) {
        nx = -nx;
        ny = -ny;
    }
    let apothem = edge_len / (2.0 * (std::f64::consts::PI / (n_sides as f64)).tan());
    let r = edge_len / (2.0 * (std::f64::consts::PI / (n_sides as f64)).sin());
    let rcx = mx + nx * apothem;
    let rcy = my + ny * apothem;
    let a_start = (p2.1 - rcy).atan2(p2.0 - rcx);
    let mut step = std::f64::consts::TAU / (n_sides as f64);
    let test_ang = a_start + step;
    if (rcx + r * test_ang.cos() - p1.0).hypot(rcy + r * test_ang.sin() - p1.1) < 1.0 {
        step = -step;
    }
    let shared_or_adj_double = b.order >= 2
        || mol.bonds.iter().any(|ob| {
            ob.order >= 2
                && (ob.begin_atom == b.begin_atom
                    || ob.end_atom == b.begin_atom
                    || ob.begin_atom == b.end_atom
                    || ob.end_atom == b.end_atom)
        });
    let mut prev_idx = b.end_atom;
    for k in 1..=(n_sides - 2) {
        let ang = a_start + (k as f64) * step;
        let nid = *next_id;
        *next_id += 1;
        let mut na = ChemAtom::new(nid, 6);
        na.x = rcx + r * ang.cos();
        na.y = rcy + r * ang.sin();
        let new_idx = mol.atoms.len();
        mol.atoms.push(na);
        let bid = *next_id;
        *next_id += 1;
        let ord = if aromatic {
            if shared_or_adj_double {
                if k % 2 == 0 { 2 } else { 1 }
            } else if k % 2 == 1 {
                2
            } else {
                1
            }
        } else if cpd && (if shared_or_adj_double { k == 2 } else { k == 1 || k == 3 }) {
            // Cyclopentadiene needs two non-adjacent C=C; if the fused atoms already carry one, only k2 is safe
            2
        } else {
            1
        };
        mol.bonds.push(ChemBond {
            id: bid,
            begin_atom: prev_idx,
            end_atom: new_idx,
            order: ord,
            display: BondDisplay::Solid,
            visible: true,
            ignore_warnings: false,
        });
        prev_idx = new_idx;
    }
    let bid = *next_id;
    *next_id += 1;
    let closing_ord = if aromatic && !shared_or_adj_double && n_sides % 2 == 0 {
        2
    } else {
        1
    };
    mol.bonds.push(ChemBond {
        id: bid,
        begin_atom: prev_idx,
        end_atom: b.begin_atom,
        order: closing_ord,
        display: BondDisplay::Solid,
        visible: true,
        ignore_warnings: false,
    });
}

/// Element of the atom through which a hotkey nickname attaches; `None` for generic labels (R, X, ?).
fn nickname_attachment_z(lbl: &str) -> Option<i32> {
    match lbl {
        "Et" | "Me" | "Ph" | "CF3" | "CO2Me" | "Cbz" | "Boc" | "Fmoc" => Some(6),
        "NO2" | "N3" => Some(7),
        "OMe" => Some(8),
        "MgBr" => Some(12),
        _ => None,
    }
}

/// Applies an exact case-sensitive hotkey character from `hotkeys.xml` for the Atom context.
pub fn apply_atom_hotkey(doc: &mut InteractiveDocument, mi: usize, ai: usize, ch: char) -> HotkeyOutcome {
    let blen = doc.settings.bond_length_pt;
    // 1. Exact 33 `LABELTEXT` bindings from `hotkeys.xml`
    let label_opt = match ch {
        'A' => Some("Ac"),
        'b' => Some("Br"),
        'B' => Some("B"),
        'c' => Some("C"),
        'C' => Some("Cl"),
        'd' => Some("D"),
        'e' => Some("Et"),
        'E' => Some("CO2Me"),
        'f' => Some("F"),
        'F' => Some("CF3"),
        'h' => Some("H"),
        'H' => Some("Cbz"),
        'i' => Some("I"),
        'l' => Some("Cl"),
        'L' => Some("Li"),
        'm' => Some("Me"),
        'M' => Some("MgBr"),
        'n' | 'w' => Some("N"),
        'N' => Some("NO2"),
        'o' | 'q' => Some("O"),
        'O' => Some("OMe"),
        'Q' => Some("Fmoc"),
        'p' => Some("P"),
        'P' => Some("Ph"),
        'r' => Some("R"),
        's' => Some("S"),
        'S' => Some("Si"),
        'x' => Some("X"),
        'y' => Some("Boc"),
        'Z' => Some("N3"),
        '!' => Some("?"),
        _ => None,
    };
    if let Some(lbl) = label_opt {
        doc.apply_scene_tx("Atom Hotkey", |s, _| {
            if let Some(a) = s.molecules.get_mut(mi).and_then(|m| m.atoms.get_mut(ai)) {
                if lbl == "D" {
                    a.z = 1;
                    a.isotope = 2;
                    a.label = Some("D".to_string());
                } else if lbl == "Ac" {
                    // ChemDraw 'A' hotkey inserts Acetyl ("Ac") nickname, NOT Actinium (Z=89)
                    a.z = 6;
                    a.isotope = 0;
                    a.label = Some("Ac".to_string());
                } else if let Some(el) = get_element_by_symbol(lbl) {
                    a.z = el.z;
                    a.isotope = 0;
                    a.label = None;
                } else {
                    // Attachment element of the group, so formula/valence don't inherit the previous atom's element
                    if let Some(z) = nickname_attachment_z(lbl) {
                        a.z = z;
                    }
                    a.isotope = 0;
                    a.label = Some(lbl.to_string());
                }
            }
        });
        return HotkeyOutcome::Handled(format!("Atom label → {lbl}"));
    }

    // 2. Charge, Atom Number, Attachment Point, Select
    match ch {
        '+' => {
            doc.apply_scene_tx("Charge +", |s, _| {
                if let Some(a) = s.molecules.get_mut(mi).and_then(|m| m.atoms.get_mut(ai)) {
                    a.charge += 1;
                }
            });
            return HotkeyOutcome::Handled("Atom charge +1".to_string());
        }
        '-' => {
            doc.apply_scene_tx("Charge -", |s, _| {
                if let Some(a) = s.molecules.get_mut(mi).and_then(|m| m.atoms.get_mut(ai)) {
                    a.charge -= 1;
                }
            });
            return HotkeyOutcome::Handled("Atom charge -1".to_string());
        }
        '\'' => return HotkeyOutcome::ToggleAtomNumbers,
        '.' => {
            doc.apply_scene_tx("Attachment Point", |s, _| {
                if let Some(a) = s.molecules.get_mut(mi).and_then(|m| m.atoms.get_mut(ai)) {
                    a.label = Some("*".to_string());
                }
            });
            return HotkeyOutcome::Handled("Added attachment point (*)".to_string());
        }
        'g' => return HotkeyOutcome::SelectMolecule(mi),
        _ => {}
    }

    // 3. Exact `SPROUT` bindings from `hotkeys.xml`
    match ch {
        '0' | '1' => {
            doc.apply_scene_tx("Sprout Bond", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        let ang = outward_sprout_angle(m, ai);
                        sprout_atom_bond(m, ai, ang, blen, 6, 1, BondDisplay::Solid, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted C-C bond".to_string())
        }
        '2' => {
            doc.apply_scene_tx("Sprout Carbonyl", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        let ang = outward_sprout_angle(m, ai);
                        sprout_atom_bond(m, ai, ang, blen, 8, 2, BondDisplay::Solid, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted carbonyl (=O)".to_string())
        }
        '3' | 'a' | 'J' => {
            let pendant = ch == 'J';
            doc.apply_scene_tx("Sprout Benzene", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        attach_ring_to_atom(m, ai, 6, true, false, pendant, blen, nid);
                    }
                }
            });
            HotkeyOutcome::Handled(if pendant { "Sprouted phenyl" } else { "Sprouted benzene" }.to_string())
        }
        '4' => {
            doc.apply_scene_tx("Sprout Wedge", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        let ang = outward_sprout_angle(m, ai);
                        sprout_atom_bond(m, ai, ang, blen, 6, 1, BondDisplay::WedgeBegin, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted wedged bond".to_string())
        }
        '5' => {
            doc.apply_scene_tx("Sprout Hashed Wedge", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        let ang = outward_sprout_angle(m, ai);
                        sprout_atom_bond(m, ai, ang, blen, 6, 1, BondDisplay::WedgedHashBegin, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted hashed-wedged bond".to_string())
        }
        '6' => {
            doc.apply_scene_tx("Sprout Cyclohexane", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        attach_ring_to_atom(m, ai, 6, false, false, false, blen, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted cyclohexane".to_string())
        }
        '7' => {
            doc.apply_scene_tx("Sprout Cyclopentane", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        attach_ring_to_atom(m, ai, 5, false, false, false, blen, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted cyclopentane".to_string())
        }
        '8' => {
            doc.apply_scene_tx("Sprout Methylidene", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        let ang = outward_sprout_angle(m, ai);
                        sprout_atom_bond(m, ai, ang, blen, 6, 2, BondDisplay::Solid, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted methylidene (=CH2)".to_string())
        }
        '9' => {
            doc.apply_scene_tx("Sprout Dimethyl", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        let ang = outward_sprout_angle(m, ai);
                        sprout_atom_bond(m, ai, ang - 0.52, blen, 6, 1, BondDisplay::Solid, nid);
                        sprout_atom_bond(m, ai, ang + 0.52, blen, 6, 1, BondDisplay::Solid, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted gem-dimethyl".to_string())
        }
        'z' => {
            doc.apply_scene_tx("Sprout Alkyne", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        let ang = outward_sprout_angle(m, ai);
                        let c1 = sprout_atom_bond(m, ai, ang, blen, 6, 1, BondDisplay::Solid, nid);
                        sprout_atom_bond(m, c1, ang, blen, 6, 3, BondDisplay::Solid, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted alkyne".to_string())
        }
        'v' => {
            doc.apply_scene_tx("Sprout Cyclopropyl", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        attach_ring_to_atom(m, ai, 3, false, false, true, blen, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted cyclopropyl".to_string())
        }
        'u' => {
            doc.apply_scene_tx("Sprout Cyclobutyl", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        attach_ring_to_atom(m, ai, 4, false, false, true, blen, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted cyclobutyl".to_string())
        }
        'K' => {
            doc.apply_scene_tx("Sprout tert-Butyl", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        let ang = outward_sprout_angle(m, ai);
                        let c1 = sprout_atom_bond(m, ai, ang, blen, 6, 1, BondDisplay::Solid, nid);
                        sprout_atom_bond(m, c1, ang, blen, 6, 1, BondDisplay::Solid, nid);
                        sprout_atom_bond(m, c1, ang - std::f64::consts::FRAC_PI_2, blen, 6, 1, BondDisplay::Solid, nid);
                        sprout_atom_bond(m, c1, ang + std::f64::consts::FRAC_PI_2, blen, 6, 1, BondDisplay::Solid, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted tert-butyl".to_string())
        }
        'k' => {
            doc.apply_scene_tx("Sprout Sulfonyl", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        let ang = outward_sprout_angle(m, ai);
                        let s_idx = sprout_atom_bond(m, ai, ang, blen, 16, 1, BondDisplay::Solid, nid);
                        sprout_atom_bond(m, s_idx, ang - 1.05, blen, 8, 2, BondDisplay::Solid, nid);
                        sprout_atom_bond(m, s_idx, ang + 1.05, blen, 8, 2, BondDisplay::Solid, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted sulfonyl (-SO2-)".to_string())
        }
        'j' => {
            doc.apply_scene_tx("Sprout Cyclopentadienyl", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    if ai < m.atoms.len() {
                        attach_ring_to_atom(m, ai, 5, false, true, true, blen, nid);
                    }
                }
            });
            HotkeyOutcome::Handled("Sprouted cyclopentadienyl".to_string())
        }
        _ => HotkeyOutcome::NotMatched,
    }
}

/// Applies an exact case-sensitive hotkey character from `hotkeys.xml` for the Bond context.
pub fn apply_bond_hotkey(doc: &mut InteractiveDocument, mi: usize, bi: usize, ch: char) -> HotkeyOutcome {
    let blen = doc.settings.bond_length_pt;
    match ch {
        '1' | '2' | '3' => {
            let ord = (ch as u8 - b'0') as i32;
            doc.apply_scene_tx("Bond Order", |s, _| {
                if let Some(b) = s.molecules.get_mut(mi).and_then(|m| m.bonds.get_mut(bi)) {
                    b.order = ord;
                    b.display = BondDisplay::Solid;
                }
            });
            HotkeyOutcome::Handled(format!("Set bond order to {ord}"))
        }
        'l' => HotkeyOutcome::SetDoubleAlign((mi, bi), DoubleBondAlign::Left),
        'r' => HotkeyOutcome::SetDoubleAlign((mi, bi), DoubleBondAlign::Right),
        'c' => HotkeyOutcome::SetDoubleAlign((mi, bi), DoubleBondAlign::Center),
        'g' => HotkeyOutcome::SelectMolecule(mi),
        'd' | 'D' | 'H' | 'h' | 'W' | 'b' | 'B' | 'w' | 'y' => {
            let (ord, disp) = match ch {
                'd' => (1, BondDisplay::Dash),
                'D' => (2, BondDisplay::Dash),
                'H' => (1, BondDisplay::Hash),
                'h' | 'W' => (1, BondDisplay::WedgedHashBegin),
                'w' => (1, BondDisplay::WedgeBegin),
                'b' => (1, BondDisplay::Bold),
                'B' => (2, BondDisplay::Bold),
                'y' => (1, BondDisplay::Wavy),
                _ => (1, BondDisplay::Solid),
            };
            doc.apply_scene_tx("Bond Display", |s, _| {
                if let Some(b) = s.molecules.get_mut(mi).and_then(|m| m.bonds.get_mut(bi)) {
                    b.order = ord;
                    b.display = disp;
                }
            });
            HotkeyOutcome::Handled(format!("Bond style updated ({ch})"))
        }
        'v' | '4' | '5' | '6' | '7' | '8' | 'a' | 'z' => {
            let (sides, arom, cpd) = match ch {
                'v' => (3, false, false),
                '4' => (4, false, false),
                '5' => (5, false, false),
                '6' => (6, false, false),
                '7' => (7, false, false),
                '8' => (8, false, false),
                'a' => (6, true, false),
                'z' => (5, false, true),
                _ => (6, false, false),
            };
            doc.apply_scene_tx("Fuse Ring Hotkey", |s, nid| {
                if let Some(m) = s.molecules.get_mut(mi) {
                    fuse_ring_to_bond(m, bi, sides, arom, cpd, blen, nid);
                }
            });
            HotkeyOutcome::Handled(format!("Fused {sides}-membered ring onto bond"))
        }
        '9' | '0' => {
            let pos = doc
                .scene
                .molecules
                .get(mi)
                .and_then(|m| m.bonds.get(bi).and_then(|b| m.atoms.get(b.begin_atom)))
                .map(|a| (a.x + blen * 1.5, a.y))
                .unwrap_or((180.0, 150.0));
            let chair = stamp_template(BuiltinTemplate::CyclohexaneChair, pos.0, pos.1, blen);
            doc.stamp_molecule_tx(chair, "Chair Hotkey");
            HotkeyOutcome::Handled("Added Cyclohexane Chair".to_string())
        }
        _ => HotkeyOutcome::NotMatched,
    }
}

/// Applies an exact case-sensitive hotkey character from `hotkeys.xml` for the Generic context.
pub fn apply_generic_hotkey(ch: char) -> HotkeyOutcome {
    match ch {
        'j' => HotkeyOutcome::SwitchTool(Tool::Benzene),
        'J' => HotkeyOutcome::SwitchTool(Tool::Cyclopentadiene),
        'x' => HotkeyOutcome::SwitchTool(Tool::SolidBond),
        'X' => HotkeyOutcome::SwitchTool(Tool::AcyclicChain),
        'e' => HotkeyOutcome::SwitchTool(Tool::Arrow),
        'E' => HotkeyOutcome::SwitchTool(Tool::ChemSymbol),
        't' => HotkeyOutcome::SwitchTool(Tool::TextLabel),
        'T' => HotkeyOutcome::SwitchTool(Tool::Bracket),
        'g' => HotkeyOutcome::SwitchTool(Tool::TlcPlate),
        'G' => HotkeyOutcome::SwitchTool(Tool::Orbital),
        ' ' => HotkeyOutcome::SwitchTool(Tool::Marquee),
        '3' => HotkeyOutcome::SwitchTool(Tool::Cyclopropane),
        '4' => HotkeyOutcome::SwitchTool(Tool::Cyclobutane),
        '5' => HotkeyOutcome::SwitchTool(Tool::Cyclopentane),
        '6' => HotkeyOutcome::SwitchTool(Tool::Cyclohexane),
        '7' => HotkeyOutcome::SwitchTool(Tool::Cycloheptane),
        '8' => HotkeyOutcome::SwitchTool(Tool::Cyclooctane),
        '=' => HotkeyOutcome::OpenChooseNickname,
        '/' | '?' => HotkeyOutcome::OpenPropertiesDialog,
        _ => HotkeyOutcome::NotMatched,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nickname_hotkeys_set_attachment_element() {
        assert_eq!(nickname_attachment_z("Et"), Some(6));
        assert_eq!(nickname_attachment_z("NO2"), Some(7));
        assert_eq!(nickname_attachment_z("OMe"), Some(8));
        assert_eq!(nickname_attachment_z("R"), None);
        // Every hotkey label that is neither an element nor generic must have an attachment element
        for ch in ('a'..='z').chain('A'..='Z') {
            let mut doc = InteractiveDocument::new_saved_seed("t", chiralcraft_doc::DocumentSettings::acs_document_1996());
            let mut mol = ChemMolecule::default();
            let mut a = ChemAtom::new(1, 16);
            a.x = 0.0;
            mol.atoms.push(a);
            doc.scene.molecules.push(mol);
            apply_atom_hotkey(&mut doc, 0, 0, ch);
            let atom = &doc.scene.molecules[0].atoms[0];
            if let Some(lbl) = atom.label.as_deref() {
                if !matches!(lbl, "R" | "X" | "?" | "D") {
                    assert_ne!(atom.z, 16, "hotkey {ch:?} label {lbl} kept previous element S");
                }
            }
        }
    }

    #[test]
    fn test_fused_cyclopentadiene_has_two_nonadjacent_double_bonds() {
        let mut mol = ChemMolecule::default();
        for (i, x) in [0.0, 14.4].iter().enumerate() {
            let mut a = ChemAtom::new(i as u32 + 1, 6);
            a.x = *x;
            mol.atoms.push(a);
        }
        mol.add_bond(3, 0, 1, 1, BondDisplay::Solid);
        let mut nid = 10;
        fuse_ring_to_bond(&mut mol, 0, 5, false, true, 14.4, &mut nid);
        assert_eq!(mol.atoms.len(), 5);
        let doubles: Vec<&ChemBond> = mol.bonds.iter().filter(|b| b.order == 2).collect();
        assert_eq!(doubles.len(), 2);
        let (a, b) = (doubles[0], doubles[1]);
        assert!(![a.begin_atom, a.end_atom].iter().any(|x| *x == b.begin_atom || *x == b.end_atom));
    }
}
