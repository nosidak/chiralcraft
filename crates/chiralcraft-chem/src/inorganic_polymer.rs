//! Organometallic & Coordination Chemistry Engine (18-Electron Rule / Green CBC)
//! and Polymer / Materials Repeat-Unit Calculator + BigSMILES Generator.

use crate::{compute_analysis_report, get_element, BondDisplay, ChemMolecule};

// ============================================================================
// 1. Organometallic & Coordination Chemistry (18-Electron Rule & Green CBC)
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub struct MetalComplexReport {
    pub metal_atom_idx: usize,
    pub metal_symbol: String,
    pub group_number: i32,
    pub oxidation_state: i32,
    pub d_electrons: i32,
    pub ligand_electrons: i32,
    pub total_valence_electrons: i32,
    pub coordination_number: usize,
    pub x_ligands: usize,
    pub l_ligands: usize,
    pub predicted_geometry: &'static str,
    pub rule_status: &'static str,
    pub ligand_details: Vec<String>,
}

/// Returns IUPAC d-block group number (`3..=12`) if `z` is a transition metal.
pub fn transition_metal_group(z: i32) -> Option<i32> {
    match z {
        21..=30 => Some(z - 18), // Sc(3) .. Zn(12)
        39..=48 => Some(z - 36), // Y(3) .. Cd(12)
        57 => Some(3),           // La
        72..=80 => Some(z - 68), // Hf(4) .. Hg(12)
        _ => None,
    }
}

/// Analyzes every transition metal center in `mol` using Green's Covalent Bond Classification (L/X)
/// to report formal oxidation state, d^n count, 18e⁻ / 16e⁻ electron count, and coordination geometry.
pub fn analyze_organometallic_complexes(mol: &ChemMolecule) -> Vec<MetalComplexReport> {
    let n = mol.atoms.len();
    let mut adj: Vec<Vec<(usize, i32, BondDisplay)>> = vec![Vec::new(); n];
    for b in &mol.bonds {
        if b.begin_atom < n && b.end_atom < n {
            adj[b.begin_atom].push((b.end_atom, b.order, b.display));
            adj[b.end_atom].push((b.begin_atom, b.order, b.display));
        }
    }

    let mut reports = Vec::new();
    for (mi, atom) in mol.atoms.iter().enumerate() {
        let Some(group) = transition_metal_group(atom.z) else {
            continue;
        };
        let sym = get_element(atom.z).map(|e| e.symbol).unwrap_or("M").to_string();
        let mut x_cnt = 0usize;
        let mut l_cnt = 0usize;
        let mut extra_e = 0i32;
        let mut details = Vec::new();

        for &(lig, order, disp) in &adj[mi] {
            let la = &mol.atoms[lig];
            let lsym = la
                .label
                .as_deref()
                .unwrap_or_else(|| get_element(la.z).map(|e| e.symbol).unwrap_or("L"));
            let l_deg = adj[lig].len();
            let lbl_low = lsym.to_ascii_lowercase();

            // Check multi-hapto label hints (Cp, Cp*, COD, bpy, phen)
            if lbl_low == "cp" || lbl_low == "cp*" || lbl_low.contains("η5") {
                // η5-Cp is L2X (6e- donor, +1 oxidation state)
                l_cnt += 2;
                x_cnt += 1;
                details.push(format!("{lsym} (η⁵, L₂X, 6e⁻)"));
                continue;
            }
            if lbl_low == "cod" || lbl_low == "nbd" || lbl_low == "bpy" || lbl_low == "dppe" || lbl_low == "dppf" {
                l_cnt += 2;
                details.push(format!("{lsym} (κ²/η⁴ bidentate, L₂, 4e⁻)"));
                continue;
            }

            if disp == BondDisplay::Dative || order == 5 {
                l_cnt += 1;
                details.push(format!("{lsym} (dative L-type, 2e⁻)"));
            } else if order == 2 {
                // Fischer/Schrock carbene or oxo/imido: treat M=O / M=N as X2 (4e-), M=C as L (2e- Schrock X2)
                if matches!(la.z, 8 | 7 | 16) && l_deg == 1 {
                    x_cnt += 2;
                    details.push(format!("={lsym} (X₂ oxo/imido, 4e⁻)"));
                } else {
                    l_cnt += 1;
                    details.push(format!("={lsym} (carbene L-type, 2e⁻)"));
                }
            } else if order >= 3 {
                x_cnt += 3;
                details.push(format!("≡{lsym} (X₃ carbyne/nitrido, 6e⁻)"));
            } else {
                // Single bond classification:
                // Halogens (F, Cl, Br, I), H, C (alkyl/aryl unless CO!), alkoxide O (deg<=2), amide N (deg==2 single) -> X-type
                // Phosphines P (deg>=3 or deg==1 PH3), carbonyl C=O, amines/ammines N (deg>=3, deg==1 neutral NH3, or aromatic/nitrile), ethers/aqua O (deg>=3) -> L-type
                let is_co = la.z == 6 && adj[lig].iter().any(|&(k, o, _)| mol.atoms[k].z == 8 && o >= 2);
                let is_nitrile_or_imine_n = la.z == 7
                    && adj[lig].iter().any(|&(k, o, _)| k != mi && (o >= 2 || o == 4));
                let is_l = is_co
                    || is_nitrile_or_imine_n
                    || (la.z == 15 && (l_deg >= 3 || (l_deg == 1 && la.charge == 0)))
                    || (la.z == 7 && (l_deg >= 3 || (l_deg == 1 && la.charge == 0)))
                    || (la.z == 8 && l_deg >= 3)
                    || (la.z == 16 && l_deg >= 3)
                    || lbl_low == "co"
                    || lbl_low.starts_with("pph")
                    || lbl_low.starts_with("pme")
                    || lbl_low.starts_with("pcy")
                    || lbl_low == "nh3"
                    || lbl_low == "h2o"
                    || lbl_low == "py"
                    || lbl_low == "thf"
                    || lbl_low == "mecn";
                if is_l {
                    l_cnt += 1;
                    details.push(format!("{lsym} (L-type neutral, 2e⁻)"));
                } else if transition_metal_group(la.z).is_some() {
                    // Metal-metal single bond contributes 1e- to each metal
                    extra_e += 1;
                    details.push(format!("M–{lsym} (metal-metal bond, 1e⁻)"));
                } else {
                    x_cnt += 1;
                    details.push(format!("{lsym} (X-type anionic, 2e⁻)"));
                }
            }
        }

        let ox = atom.charge + (x_cnt as i32);
        let d_elec = (group - ox).max(0);
        // Ionic method: total valence e- = d_elec + 2*(x_cnt + l_cnt) + extra_e
        // Equivalent to neutral method: group - charge + x_cnt + 2*l_cnt + extra_e!
        let lig_e = 2 * ((x_cnt + l_cnt) as i32) + extra_e;
        let total_e = d_elec + lig_e;
        let cn = adj[mi].len();

        let geometry = match cn {
            0 => "Unbound Metal Atom",
            1 => "Monocoordinate",
            2 => "Linear (d¹⁰ / 14e⁻ typical)",
            3 => "Trigonal Planar",
            4 => {
                if d_elec == 8 && matches!(atom.z, 45 | 46 | 77 | 78 | 79 | 28) {
                    "Square Planar (d⁸ 16e⁻ classic)"
                } else {
                    "Tetrahedral"
                }
            }
            5 => "Trigonal Bipyramidal / Square Pyramidal",
            6 => "Octahedral",
            _ => "High-Coordination Polyhedron",
        };

        let status = if total_e == 18 {
            "18e⁻ Closed-Shell (Saturated)"
        } else if total_e == 16 && d_elec == 8 {
            "16e⁻ Square-Planar Stable (d⁸)"
        } else if total_e < 18 {
            "Coordinatively Unsaturated (<18e⁻)"
        } else {
            "Hypervalent / >18e⁻ Complex"
        };

        reports.push(MetalComplexReport {
            metal_atom_idx: mi,
            metal_symbol: sym,
            group_number: group,
            oxidation_state: ox,
            d_electrons: d_elec,
            ligand_electrons: lig_e,
            total_valence_electrons: total_e,
            coordination_number: cn,
            x_ligands: x_cnt,
            l_ligands: l_cnt,
            predicted_geometry: geometry,
            rule_status: status,
            ligand_details: details,
        });
    }
    reports
}

// ============================================================================
// 2. Polymer & Materials Repeat-Unit Calculator + BigSMILES
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CopolymerTopology {
    Homopolymer,
    Random,
    Block,
    Alternating,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolymerReport {
    pub repeat_formula: String,
    pub repeat_mw: f64,
    pub effective_unit_mw: f64,
    pub end_group_mw: f64,
    pub dp: f64,
    pub mn: f64,
    pub mw: f64,
    pub pdi: f64,
    pub weight_frac_a: f64,
    pub big_smiles: String,
}

/// Computes polymer number-average ($M_n$), weight-average ($M_w$), copolymer weight fraction,
/// and standard **BigSMILES** notation from the drawn repeat unit.
pub fn compute_polymer_metrics(
    repeat_unit: &ChemMolecule,
    repeat_smiles: &str,
    dp: f64,
    pdi: f64,
    end_group_mw: f64,
    comonomer_mw: f64,
    mol_frac_a: f64,
    topology: CopolymerTopology,
) -> PolymerReport {
    let rep = compute_analysis_report(repeat_unit);
    // Subtract 2 H masses (2.016 Da) when converting drawn closed monomer into divalent repeat unit
    let ru_mw = (rep.mol_weight - 2.01588).max(12.0);
    let fa = mol_frac_a.clamp(0.01, 1.0);
    let fb = 1.0 - fa;
    let eff_mw = match topology {
        CopolymerTopology::Homopolymer => ru_mw,
        CopolymerTopology::Alternating => 0.5 * ru_mw + 0.5 * comonomer_mw.max(12.0),
        CopolymerTopology::Random | CopolymerTopology::Block => fa * ru_mw + fb * comonomer_mw.max(12.0),
    };
    let w_a = match topology {
        CopolymerTopology::Homopolymer => 1.0,
        CopolymerTopology::Alternating => ru_mw / (ru_mw + comonomer_mw.max(12.0)),
        CopolymerTopology::Random | CopolymerTopology::Block => (fa * ru_mw) / eff_mw.max(1e-6),
    };
    let dp_c = dp.max(1.0);
    let pdi_c = pdi.max(1.0);
    let mn = dp_c * eff_mw + end_group_mw.max(0.0);
    let mw = mn * pdi_c;

    let clean_smi = if repeat_smiles.trim().is_empty() { "CC" } else { repeat_smiles.trim() };
    let big_smiles = match topology {
        CopolymerTopology::Homopolymer => format!("{{[$]{clean_smi}[$]}}"),
        CopolymerTopology::Random => format!("{{[$]{clean_smi}[$],[$]C(C)([$])}}"),
        CopolymerTopology::Block => format!("{{[$]{clean_smi}[$]}}{{[$]C(C)[$]}}"),
        CopolymerTopology::Alternating => format!("{{[<]{clean_smi}[>],[>]C(C)[<]}}"),
    };

    PolymerReport {
        repeat_formula: rep.formula,
        repeat_mw: ru_mw,
        effective_unit_mw: eff_mw,
        end_group_mw: end_group_mw.max(0.0),
        dp: dp_c,
        mn,
        mw,
        pdi: pdi_c,
        weight_frac_a: w_a,
        big_smiles,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ChemAtom;

    #[test]
    fn test_cisplatin_16e_square_planar_and_polymer_metrics() {
        // Cisplatin: cis-[Pt(NH3)2Cl2] -> Pt(II), d8, 16e-, Square Planar!
        let mut cis = ChemMolecule::new();
        cis.add_atom(ChemAtom::new(1, 78)); // Pt (group 10)
        cis.add_atom(ChemAtom::new(2, 17)); // Cl (X)
        cis.add_atom(ChemAtom::new(3, 17)); // Cl (X)
        let mut n1 = ChemAtom::new(4, 7); n1.label = Some("NH3".into());
        let mut n2 = ChemAtom::new(5, 7); n2.label = Some("NH3".into());
        cis.add_atom(n1);
        cis.add_atom(n2);
        cis.add_bond(10, 0, 1, 1, BondDisplay::Solid);
        cis.add_bond(11, 0, 2, 1, BondDisplay::Solid);
        cis.add_bond(12, 0, 3, 1, BondDisplay::Dative);
        cis.add_bond(13, 0, 4, 1, BondDisplay::Dative);

        let reps = analyze_organometallic_complexes(&cis);
        assert_eq!(reps.len(), 1);
        assert_eq!(reps[0].oxidation_state, 2);
        assert_eq!(reps[0].d_electrons, 8);
        assert_eq!(reps[0].total_valence_electrons, 16);
        assert!(reps[0].predicted_geometry.contains("Square Planar"));

        // Polyethylene repeat unit (-CH2-CH2-)
        let mut pe = ChemMolecule::new();
        pe.add_atom(ChemAtom::new(1, 6));
        pe.add_atom(ChemAtom::new(2, 6));
        pe.add_bond(10, 0, 1, 1, BondDisplay::Solid);
        let poly = compute_polymer_metrics(&pe, "CC", 100.0, 1.15, 2.02, 104.15, 1.0, CopolymerTopology::Homopolymer);
        assert!((poly.repeat_mw - 28.05).abs() < 0.2);
        assert!((poly.mn - 2807.0).abs() < 20.0);
        assert!(poly.big_smiles.contains("{[$]CC[$]}"));
    }
}
