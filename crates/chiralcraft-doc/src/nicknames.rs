//! 1:1 Mutable Single-Table Nickname Engine (`537` official ChiralCraft shipped factory nicknames)
//! plus real graph-splicing `Expand Label` / `Contract Label` (`P17`, `P22`, `P31`, `P33d`, `P35c`).
//!
//! Verified runtime semantics (`P35c`):
//! - Built-in and user-defined nicknames live in a single mutable table loaded from `SHIPPED_NICKNAMES`.
//! - Redefining a name that matches an existing nickname AND an element symbol (e.g. `"Ac"` = acetyl + Actinium)
//!   requires two sequential confirmations: `ExistingNickname` (`[Replace]`) and `ElementSymbol` (`[Override]`).
//! - Deleting a nickname via `List Nicknames > Delete` removes it from the table completely so `"Ac"` falls
//!   through to periodic table element `Ac` (`Z=89`, Actinium).

use std::collections::{BTreeMap, HashSet};

use chiralcraft_chem::{
    get_element, BondDisplay, ChemAtom, ChemBond, ChemMolecule, SHIPPED_NICKNAMES,
};

#[derive(Debug, Clone, PartialEq)]
pub struct NicknameDef {
    pub name: String,
    pub fragment: ChemMolecule,
    /// Index of the attachment atom in `fragment` (bonds to external parent).
    pub attach_atom_idx: usize,
    /// Attachment bond vector pointing from `attach_atom_idx` toward the ExternalConnectionPoint.
    pub attach_dx: f64,
    pub attach_dy: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NicknameCollision {
    None,
    ExistingNickname,
    ElementSymbol(i32),
    BothNicknameAndElement(i32),
}

#[derive(Debug, Clone, PartialEq)]
pub enum LabelResolution {
    Nickname(NicknameDef),
    Element(i32),
    Unresolved,
}

/// Metadata tracking an expanded nickname residue (ChiralCraft `<chemicalproperty Name="Residue">`),
/// enabling symmetric `Contract Label` (`P22`, `P31`).
#[derive(Debug, Clone, PartialEq)]
pub struct ExpandedResidue {
    pub name: String,
    pub parent_atom_id: u32,
    pub expanded_atom_ids: Vec<u32>,
    pub anchor_x: f64,
    pub anchor_y: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NicknameTable {
    entries: BTreeMap<String, NicknameDef>,
}

impl Default for NicknameTable {
    fn default() -> Self {
        Self::with_factory_defaults()
    }
}

impl NicknameTable {
    /// Loads all `537` official ChiralCraft shipped factory nicknames verbatim.
    pub fn with_factory_defaults() -> Self {
        let mut entries = BTreeMap::new();
        for sn in SHIPPED_NICKNAMES {
            let mut frag = ChemMolecule::default();
            for (i, a) in sn.atoms.iter().enumerate() {
                let mut ca = ChemAtom::new((i + 1) as u32, a.z);
                ca.charge = a.charge;
                ca.isotope = a.isotope;
                ca.x = a.x;
                ca.y = a.y;
                frag.atoms.push(ca);
            }
            for (i, b) in sn.bonds.iter().enumerate() {
                let disp = match b.disp {
                    1 => BondDisplay::WedgeBegin,
                    2 => BondDisplay::WedgedHashBegin,
                    3 => BondDisplay::WedgeEnd,
                    4 => BondDisplay::WedgedHashEnd,
                    _ => BondDisplay::Solid,
                };
                frag.bonds.push(ChemBond {
                    id: (i + 1) as u32,
                    begin_atom: b.b as usize,
                    end_atom: b.e as usize,
                    order: b.order,
                    display: disp,
                    visible: true,
                    ignore_warnings: false,
                });
            }
            entries.insert(
                sn.name.to_string(),
                NicknameDef {
                    name: sn.name.to_string(),
                    fragment: frag,
                    attach_atom_idx: sn.attach_idx as usize,
                    attach_dx: sn.attach_dx,
                    attach_dy: sn.attach_dy,
                },
            );
        }
        Self { entries }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn lookup(&self, name: &str) -> Option<&NicknameDef> {
        self.entries.get(name)
    }

    pub fn list_names(&self) -> Vec<String> {
        self.entries.keys().cloned().collect()
    }

    /// Check ChiralCraft's two-stage collision prompts (`P33d`, `P35c`).
    pub fn check_collision(&self, name: &str) -> NicknameCollision {
        let has_nick = self.entries.contains_key(name);
        let elem_z = (1..=118).find(|&z| {
            get_element(z)
                .map(|e| e.symbol == name)
                .unwrap_or(false)
        });
        match (has_nick, elem_z) {
            (true, Some(z)) => NicknameCollision::BothNicknameAndElement(z),
            (true, None) => NicknameCollision::ExistingNickname,
            (false, Some(z)) => NicknameCollision::ElementSymbol(z),
            (false, None) => NicknameCollision::None,
        }
    }

    /// Define or override a nickname in the single mutable table.
    pub fn define_nickname(
        &mut self,
        name: &str,
        fragment: ChemMolecule,
        attach_atom_idx: usize,
        confirm_replace_nickname: bool,
        confirm_override_element: bool,
    ) -> Result<(), NicknameCollision> {
        let col = self.check_collision(name);
        match col {
            NicknameCollision::None => {}
            NicknameCollision::ExistingNickname if confirm_replace_nickname => {}
            NicknameCollision::ElementSymbol(_) if confirm_override_element => {}
            NicknameCollision::BothNicknameAndElement(_)
                if confirm_replace_nickname && confirm_override_element => {}
            other => return Err(other),
        }
        self.entries.insert(
            name.to_string(),
            NicknameDef {
                name: name.to_string(),
                fragment,
                attach_atom_idx,
                attach_dx: -14.4,
                attach_dy: 0.0,
            },
        );
        Ok(())
    }

    /// Delete a nickname from the single mutable table (`List Nicknames > Delete`, `P35c`).
    /// Does NOT revert to a hidden factory layer (`"Ac"` falls through to Actinium `Z=89`).
    pub fn delete_nickname(&mut self, name: &str) -> bool {
        self.entries.remove(name).is_some()
    }

    /// Resolve a label string: nickname table first, then periodic table element symbol.
    pub fn resolve_label(&self, label: &str) -> LabelResolution {
        if let Some(def) = self.entries.get(label) {
            LabelResolution::Nickname(def.clone())
        } else if let Some(z) = (1..=118).find(|&z| {
            get_element(z)
                .map(|e| e.symbol == label)
                .unwrap_or(false)
        }) {
            LabelResolution::Element(z)
        } else {
            LabelResolution::Unresolved
        }
    }

    /// Real graph-splicing `Structure > Expand Label` (`P17`, `P31`, `P35c`):
    /// Replaces `atom_idx` in `mol` with the nickname's full chemical graph, oriented along
    /// the incident parent bond vector and scaled to `target_bond_len` (default `14.40 pt`).
    /// Returns an `ExpandedResidue` handle for symmetric `contract_residue_group`.
    pub fn expand_label_at_atom(
        &self,
        mol: &mut ChemMolecule,
        atom_idx: usize,
        label: &str,
        target_bond_len: f64,
        next_id: &mut u32,
    ) -> Option<ExpandedResidue> {
        let normalized: String = label
            .trim()
            .chars()
            .map(|c| match c {
                '₀' => '0', '₁' => '1', '₂' => '2', '₃' => '3', '₄' => '4',
                '₅' => '5', '₆' => '6', '₇' => '7', '₈' => '8', '₉' => '9',
                other => other,
            })
            .collect();
        let canonical = match normalized.as_str() {
            "O2N" => "NO2",
            "MeO" => "OMe",
            "EtO" => "OEt",
            "iPrO" => "OiPr",
            "tBuO" => "OtBu",
            "AcO" => "OAc",
            "PhO" => "OPh",
            "BnO" => "OBn",
            "TBSO" => "OTBS",
            "TMSO" => "OTMS",
            "TsO" => "OTs",
            "MsO" => "OMs",
            "TfO" => "OTf",
            "F3C" => "CF3",
            "Cl3C" => "CCl3",
            "NC" => "CN",
            "HO2C" => "CO2H",
            "MeO2C" => "CO2Me",
            "EtO2C" => "CO2Et",
            "BocHN" => "NHBoc",
            "CbzHN" => "NHCbz",
            "FmocHN" => "NHFmoc",
            "AcHN" => "NHAc",
            "MeHN" => "NHMe",
            "Me2N" => "NMe2",
            "Et2N" => "NEt2",
            "O2S" => "SO2",
            other => other,
        };
        let owned_fallback;
        let def: &NicknameDef = if let Some(d) = self
            .entries
            .get(label)
            .or_else(|| self.entries.get(normalized.as_str()))
            .or_else(|| self.entries.get(canonical))
        {
            d
        } else if let Some(fb) = builtin_functional_group_def(canonical) {
            owned_fallback = fb;
            &owned_fallback
        } else {
            return None;
        };
        if atom_idx >= mol.atoms.len() || def.fragment.atoms.is_empty() {
            return None;
        }
        let target_x = mol.atoms[atom_idx].x;
        let target_y = mol.atoms[atom_idx].y;

        // Find parent neighbor (if attached by a single bond) to orient the nickname cleanly
        let parent_nb = mol
            .bonds
            .iter()
            .find_map(|b| {
                if b.begin_atom == atom_idx {
                    Some(b.end_atom)
                } else if b.end_atom == atom_idx {
                    Some(b.begin_atom)
                } else {
                    None
                }
            });
        let parent_atom_id = parent_nb.map(|i| mol.atoms[i].id).unwrap_or(0);

        // Compute rotation angle & scale from nickname's ExternalConnectionPoint vector to parent vector
        let def_len = def.attach_dx.hypot(def.attach_dy);
        let scale = if def_len > 1e-3 {
            target_bond_len / def_len
        } else {
            1.0
        };
        let rot = if let Some(p_idx) = parent_nb {
            let pdx = mol.atoms[p_idx].x - target_x;
            let pdy = mol.atoms[p_idx].y - target_y;
            let want_ang = pdy.atan2(pdx);
            let def_ang = def.attach_dy.atan2(def.attach_dx);
            want_ang - def_ang
        } else {
            0.0
        };
        let (sin_r, cos_r) = rot.sin_cos();

        let attach_src = def.attach_atom_idx.min(def.fragment.atoms.len() - 1);
        let mut idx_map = vec![0usize; def.fragment.atoms.len()];
        let mut expanded_atom_ids = Vec::with_capacity(def.fragment.atoms.len());

        // Replace `atom_idx` in-place with the attachment atom of the nickname
        {
            let src = &def.fragment.atoms[attach_src];
            let dst = &mut mol.atoms[atom_idx];
            dst.z = src.z;
            dst.charge = src.charge;
            dst.isotope = src.isotope;
            dst.radical = src.radical;
            // Stale H-count/valence from the label atom would corrupt perception of the new attach atom.
            dst.valence_override = src.valence_override;
            dst.explicit_num_hydrogens = src.explicit_num_hydrogens;
            dst.label = None;
            dst.is_query = false;
            dst.is_unspecified_label = false;
            idx_map[attach_src] = atom_idx;
            expanded_atom_ids.push(dst.id);
        }

        // Append remaining nickname atoms (subtract attachment atom origin so user-defined nicknames never drift)
        let attach_x = def.fragment.atoms.get(attach_src).map(|a| a.x).unwrap_or(0.0);
        let attach_y = def.fragment.atoms.get(attach_src).map(|a| a.y).unwrap_or(0.0);
        for (i, src) in def.fragment.atoms.iter().enumerate() {
            if i == attach_src {
                continue;
            }
            let id = *next_id;
            *next_id = next_id.saturating_add(1);
            let dx = src.x - attach_x;
            let dy = src.y - attach_y;
            let rx = (dx * cos_r - dy * sin_r) * scale;
            let ry = (dx * sin_r + dy * cos_r) * scale;
            let mut a = ChemAtom::new(id, src.z);
            a.charge = src.charge;
            a.isotope = src.isotope;
            a.radical = src.radical;
            a.x = target_x + rx;
            a.y = target_y + ry;
            idx_map[i] = mol.atoms.len();
            expanded_atom_ids.push(id);
            mol.atoms.push(a);
        }

        // Append nickname internal bonds
        for b in &def.fragment.bonds {
            if b.begin_atom >= idx_map.len() || b.end_atom >= idx_map.len() {
                continue;
            }
            let bid = *next_id;
            *next_id = next_id.saturating_add(1);
            mol.bonds.push(ChemBond {
                id: bid,
                begin_atom: idx_map[b.begin_atom],
                end_atom: idx_map[b.end_atom],
                order: b.order,
                display: b.display,
                visible: true,
                ignore_warnings: false,
            });
        }

        Some(ExpandedResidue {
            name: def.name.clone(),
            parent_atom_id,
            expanded_atom_ids,
            anchor_x: target_x,
            anchor_y: target_y,
        })
    }

    /// Real symmetric `Structure > Contract Label` (`P22`, `P31`):
    /// Collapses an `ExpandedResidue` back into a single labeled attachment atom.
    pub fn contract_residue_group(mol: &mut ChemMolecule, res: &ExpandedResidue) -> bool {
        if res.expanded_atom_ids.is_empty() {
            return false;
        }
        let anchor_id = res.expanded_atom_ids[0];
        let drop_ids: HashSet<u32> = res.expanded_atom_ids[1..].iter().copied().collect();
        let Some(anchor_old_idx) = mol.atoms.iter().position(|a| a.id == anchor_id) else {
            return false;
        };

        // Restore anchor atom as the contracted nickname placeholder
        {
            let a = &mut mol.atoms[anchor_old_idx];
            a.z = 6;
            a.charge = 0;
            a.isotope = 0;
            a.radical = 0;
            a.valence_override = None;
            a.explicit_num_hydrogens = None;
            a.x = res.anchor_x;
            a.y = res.anchor_y;
            a.label = Some(res.name.clone());
        }

        // Build old->new atom index map while removing dropped atoms
        let mut old_to_new: Vec<Option<usize>> = vec![None; mol.atoms.len()];
        let mut new_atoms = Vec::with_capacity(mol.atoms.len().saturating_sub(drop_ids.len()));
        for (old_i, a) in mol.atoms.drain(..).enumerate() {
            if !drop_ids.contains(&a.id) {
                old_to_new[old_i] = Some(new_atoms.len());
                new_atoms.push(a);
            }
        }
        mol.atoms = new_atoms;

        // Retain only bonds whose endpoints survived and remap indices
        let mut new_bonds = Vec::with_capacity(mol.bonds.len());
        for mut b in mol.bonds.drain(..) {
            if let (Some(nb), Some(ne)) = (
                old_to_new.get(b.begin_atom).copied().flatten(),
                old_to_new.get(b.end_atom).copied().flatten(),
            ) {
                b.begin_atom = nb;
                b.end_atom = ne;
                new_bonds.push(b);
            }
        }
        mol.bonds = new_bonds;
        true
    }
}

/// On-the-fly functional-group graph definitions for common abbreviations (`NO2`, `CF3`, `CO2H`,
/// `OMe`, `OEt`, `OAc`, `CN`, `CO2Me`, `CO2Et`, `NHAc`, `NMe2`) so `Structure > Expand Label`
/// works seamlessly on both forward and left-facing reversed labels without altering the 537 table count.
fn builtin_functional_group_def(name: &str) -> Option<NicknameDef> {
    let l = 14.4_f64;
    let c30 = 30.0_f64.to_radians().cos() * l;
    let s30 = 30.0_f64.to_radians().sin() * l;
    let mut m = ChemMolecule::default();
    let add_a = |m: &mut ChemMolecule, z: i32, charge: i32, x: f64, y: f64| -> usize {
        let mut a = ChemAtom::new((m.atoms.len() + 1) as u32, z);
        a.charge = charge;
        a.x = x;
        a.y = y;
        m.atoms.push(a);
        m.atoms.len() - 1
    };
    let add_b = |m: &mut ChemMolecule, b: usize, e: usize, order: i32| {
        m.bonds.push(ChemBond {
            id: (m.bonds.len() + 1) as u32,
            begin_atom: b,
            end_atom: e,
            order,
            display: BondDisplay::Solid,
            visible: true,
            ignore_warnings: false,
        });
    };
    match name {
        "NO2" => {
            let n = add_a(&mut m, 7, 1, 0.0, 0.0);
            let o1 = add_a(&mut m, 8, 0, c30, -s30);
            let o2 = add_a(&mut m, 8, -1, c30, s30);
            add_b(&mut m, n, o1, 2);
            add_b(&mut m, n, o2, 1);
        }
        "CF3" | "CCl3" => {
            let hz = if name == "CF3" { 9 } else { 17 };
            let c = add_a(&mut m, 6, 0, 0.0, 0.0);
            let h1 = add_a(&mut m, hz, 0, 0.0, -l);
            let h2 = add_a(&mut m, hz, 0, l, 0.0);
            let h3 = add_a(&mut m, hz, 0, 0.0, l);
            add_b(&mut m, c, h1, 1);
            add_b(&mut m, c, h2, 1);
            add_b(&mut m, c, h3, 1);
        }
        "CO2H" | "COOH" => {
            let c = add_a(&mut m, 6, 0, 0.0, 0.0);
            let o1 = add_a(&mut m, 8, 0, c30, -s30);
            let o2 = add_a(&mut m, 8, 0, c30, s30);
            add_b(&mut m, c, o1, 2);
            add_b(&mut m, c, o2, 1);
        }
        "CO2Me" => {
            let c = add_a(&mut m, 6, 0, 0.0, 0.0);
            let o1 = add_a(&mut m, 8, 0, c30, -s30);
            let o2 = add_a(&mut m, 8, 0, c30, s30);
            let me = add_a(&mut m, 6, 0, 2.0 * c30, 0.0);
            add_b(&mut m, c, o1, 2);
            add_b(&mut m, c, o2, 1);
            add_b(&mut m, o2, me, 1);
        }
        "CO2Et" => {
            let c = add_a(&mut m, 6, 0, 0.0, 0.0);
            let o1 = add_a(&mut m, 8, 0, c30, -s30);
            let o2 = add_a(&mut m, 8, 0, c30, s30);
            let c1 = add_a(&mut m, 6, 0, 2.0 * c30, 0.0);
            let c2 = add_a(&mut m, 6, 0, 3.0 * c30, s30);
            add_b(&mut m, c, o1, 2);
            add_b(&mut m, c, o2, 1);
            add_b(&mut m, o2, c1, 1);
            add_b(&mut m, c1, c2, 1);
        }
        "OMe" => {
            let o = add_a(&mut m, 8, 0, 0.0, 0.0);
            let c = add_a(&mut m, 6, 0, c30, -s30);
            add_b(&mut m, o, c, 1);
        }
        "OEt" => {
            let o = add_a(&mut m, 8, 0, 0.0, 0.0);
            let c1 = add_a(&mut m, 6, 0, c30, -s30);
            let c2 = add_a(&mut m, 6, 0, 2.0 * c30, 0.0);
            add_b(&mut m, o, c1, 1);
            add_b(&mut m, c1, c2, 1);
        }
        "OAc" => {
            let o = add_a(&mut m, 8, 0, 0.0, 0.0);
            let c = add_a(&mut m, 6, 0, c30, -s30);
            let o2 = add_a(&mut m, 8, 0, c30, -s30 - l);
            let me = add_a(&mut m, 6, 0, 2.0 * c30, 0.0);
            add_b(&mut m, o, c, 1);
            add_b(&mut m, c, o2, 2);
            add_b(&mut m, c, me, 1);
        }
        "CN" => {
            let c = add_a(&mut m, 6, 0, 0.0, 0.0);
            let n = add_a(&mut m, 7, 0, l, 0.0);
            add_b(&mut m, c, n, 3);
        }
        "NHAc" => {
            let n = add_a(&mut m, 7, 0, 0.0, 0.0);
            let c = add_a(&mut m, 6, 0, c30, -s30);
            let o = add_a(&mut m, 8, 0, c30, -s30 - l);
            let me = add_a(&mut m, 6, 0, 2.0 * c30, 0.0);
            add_b(&mut m, n, c, 1);
            add_b(&mut m, c, o, 2);
            add_b(&mut m, c, me, 1);
        }
        "NMe2" => {
            let n = add_a(&mut m, 7, 0, 0.0, 0.0);
            let m1 = add_a(&mut m, 6, 0, c30, -s30);
            let m2 = add_a(&mut m, 6, 0, c30, s30);
            add_b(&mut m, n, m1, 1);
            add_b(&mut m, n, m2, 1);
        }
        _ => return None,
    }
    Some(NicknameDef {
        name: name.to_string(),
        fragment: m,
        attach_atom_idx: 0,
        attach_dx: -l,
        attach_dy: 0.0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand_contract_reset_stale_h_count_and_valence() {
        let mut mol = ChemMolecule::new();
        let mut c = ChemAtom::new(1, 6);
        c.x = 0.0;
        mol.add_atom(c);
        let mut lab = ChemAtom::new(2, 6);
        lab.x = 14.4;
        lab.label = Some("OMe".to_string());
        lab.explicit_num_hydrogens = Some(0);
        lab.valence_override = Some(1);
        mol.add_atom(lab);
        mol.add_bond(10, 0, 1, 1, chiralcraft_chem::BondDisplay::Solid);
        let table = NicknameTable::with_factory_defaults();
        let mut nid = 100u32;
        let res = table.expand_label_at_atom(&mut mol, 1, "OMe", 14.4, &mut nid).expect("OMe expands");
        assert_eq!(mol.atoms[1].z, 8);
        assert_eq!(mol.atoms[1].explicit_num_hydrogens, None);
        assert_eq!(mol.atoms[1].valence_override, None);
        mol.atoms[1].explicit_num_hydrogens = Some(1);
        assert!(NicknameTable::contract_residue_group(&mut mol, &res));
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.atoms[1].label.as_deref(), Some("OMe"));
        assert_eq!(mol.atoms[1].explicit_num_hydrogens, None);
    }
}
