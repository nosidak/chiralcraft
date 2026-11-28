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
/// to report formal oxidation state, d^n count, 18e- / 16e- electron count, and coordination geometry.
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

        // Multi-hapto ligands drawn explicitly (e.g. η5-Cp as five M–C bonds, η2-alkene as two):
        // mutually bonded unlabeled carbons on the metal form one η^k ligand (L_{k/2} or L_{(k-1)/2}X).
        let hapto_cand: Vec<usize> = {
            let mut v: Vec<usize> = adj[mi]
                .iter()
                .filter(|&&(k, o, d)| {
                    mol.atoms[k].z == 6 && mol.atoms[k].label.is_none() && o == 1 && d != BondDisplay::Dative
                })
                .map(|&(k, _, _)| k)
                .collect();
            v.sort_unstable();
            v.dedup();
            v
        };
        let mut in_hapto = vec![false; n];
        let mut hapto_ligands = 0usize;
        let mut seen = vec![false; n];
        for &s in &hapto_cand {
            if seen[s] {
                continue;
            }
            let mut comp = vec![s];
            seen[s] = true;
            let mut q = 0;
            while q < comp.len() {
                let u = comp[q];
                q += 1;
                for &(w, _, _) in &adj[u] {
                    if !seen[w] && hapto_cand.binary_search(&w).is_ok() {
                        seen[w] = true;
                        comp.push(w);
                    }
                }
            }
            let k = comp.len();
            if k >= 2 {
                hapto_ligands += 1;
                for &c in &comp {
                    in_hapto[c] = true;
                }
                l_cnt += k / 2;
                x_cnt += k % 2;
                let kind = if k % 2 == 1 { format!("L{}X", k / 2) } else { format!("L{}", k / 2) };
                details.push(format!("η{k}-C{k} ({kind}, {}e-)", 2 * (k / 2) + 2 * (k % 2)));
            }
        }

        for &(lig, order, disp) in &adj[mi] {
            if in_hapto[lig] {
                continue;
            }
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
                details.push(format!("{lsym} (η5, L2X, 6e-)"));
                continue;
            }
            if lbl_low == "cod" || lbl_low == "nbd" || lbl_low == "bpy" || lbl_low == "dppe" || lbl_low == "dppf" {
                l_cnt += 2;
                details.push(format!("{lsym} (κ2/η4 bidentate, L2, 4e-)"));
                continue;
            }

            if disp == BondDisplay::Dative || order == 5 {
                l_cnt += 1;
                details.push(format!("{lsym} (dative L-type, 2e-)"));
            } else if transition_metal_group(la.z).is_some() {
                // Metal-metal bond (single/double/triple/quadruple) contributes `order` e- to each metal
                let mm = order.clamp(1, 4);
                extra_e += mm;
                details.push(format!("M–{lsym} (metal-metal bond order {mm}, {mm}e-)"));
            } else if order == 2 {
                // Fischer/Schrock carbene or oxo/imido: M=O / M=S (deg 1) or M=NR (deg <= 2) -> X2 (4e-), M=CR2 -> L (2e-)
                if (matches!(la.z, 8 | 16) && l_deg == 1) || (la.z == 7 && l_deg <= 2) {
                    x_cnt += 2;
                    details.push(format!("={lsym} (X2 oxo/imido, 4e-)"));
                } else {
                    l_cnt += 1;
                    details.push(format!("={lsym} (carbene L-type, 2e-)"));
                }
            } else if order >= 3 {
                x_cnt += 3;
                details.push(format!("≡{lsym} (X3 carbyne/nitrido, 6e-)"));
            } else {
                // Single bond classification:
                // Halogens (F, Cl, Br, I), H, C (alkyl/aryl unless CO!), alkoxide O (deg<=2), amide N (deg==2 single) -> X-type
                // Phosphines P (deg>=3 or deg==1 PH3), carbonyl C=O (deg==2), amines/ammines N (deg>=3, deg==1 neutral NH3, or aromatic/nitrile), ethers/aqua O (deg>=3) -> L-type
                let is_co = la.z == 6 && l_deg == 2 && adj[lig].iter().any(|&(k, o, _)| mol.atoms[k].z == 8 && o >= 2);
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
                    details.push(format!("{lsym} (L-type neutral, 2e-)"));
                } else {
                    x_cnt += 1;
                    details.push(format!("{lsym} (X-type anionic, 2e-)"));
                }
            }
        }

        let ox = atom.charge + (x_cnt as i32);
        let d_elec = (group - ox).max(0);
        // Neutral method: total valence e- = group - charge + x_cnt + 2*l_cnt + extra_e (never fabricates e- if ox > group)
        let lig_e = 2 * ((x_cnt + l_cnt) as i32) + extra_e;
        let total_e = (group - atom.charge + (x_cnt as i32) + 2 * (l_cnt as i32) + extra_e).max(0);
        let cn = adj[mi].iter().filter(|&&(k, _, _)| !in_hapto[k]).count() + hapto_ligands;

        let geometry = match cn {
            0 => "Unbound Metal Atom",
            1 => "Monocoordinate",
            2 => "Linear (d¹⁰ / 14e- typical)",
            3 => "Trigonal Planar",
            4 => {
                if d_elec == 8 && matches!(atom.z, 45 | 46 | 77 | 78 | 79 | 28) {
                    "Square Planar (d⁸ 16e- classic)"
                } else {
                    "Tetrahedral"
                }
            }
            5 => "Trigonal Bipyramidal / Square Pyramidal",
            6 => "Octahedral",
            _ => "High-Coordination Polyhedron",
        };

        let status = if total_e == 18 {
            "18e- Closed-Shell (Saturated)"
        } else if total_e == 16 && d_elec == 8 {
            "16e- Square-Planar Stable (d⁸)"
        } else if total_e < 18 {
            "Coordinatively Unsaturated (<18e-)"
        } else {
            "Hypervalent / >18e- Complex"
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
    // f64::clamp propagates NaN; treat a NaN fraction as pure A
    let fa = if mol_frac_a.is_nan() { 1.0 } else { mol_frac_a.clamp(0.01, 1.0) };
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
        CopolymerTopology::Alternating => format!("{{[<]{clean_smi}[<],[>]C(C)[>]}}"),
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
        let nan = compute_polymer_metrics(&pe, "CC", 100.0, 1.0, 0.0, 50.0, f64::NAN, CopolymerTopology::Random);
        assert!(nan.mn.is_finite());
    }

    #[test]
    fn test_ferrocene_explicit_eta5_bonds() {
        let mut m = ChemMolecule::new();
        m.add_atom(ChemAtom::new(1, 26)); // Fe
        for ring in 0..2 {
            let base = m.atoms.len();
            for j in 0..5 {
                m.add_atom(ChemAtom::new((base + j + 1) as u32, 6));
            }
            for j in 0..5 {
                let o = if j % 2 == 0 && j < 4 { 2 } else { 1 };
                m.add_bond((100 + ring * 10 + j) as u32, base + j, base + (j + 1) % 5, o, BondDisplay::Solid);
                m.add_bond((200 + ring * 10 + j) as u32, 0, base + j, 1, BondDisplay::Solid);
            }
        }
        let r = &analyze_organometallic_complexes(&m)[0];
        assert_eq!(r.oxidation_state, 2);
        assert_eq!(r.d_electrons, 6);
        assert_eq!(r.total_valence_electrons, 18);
        assert_eq!(r.coordination_number, 2);
    }
}
