//! Physiological pH Microspecies & Isoelectric Point (pI), Tautomer Enumerator,
//! Stereoisomer Enumerator, Resonance Structure Generator, and Curved-Arrow Electron Pusher.

use std::collections::BTreeSet;

use crate::{
    perceive_stereochemistry, BondDisplay, CcAromaticityInfo, CcVal, ChemMolecule,
};

fn build_adj(mol: &ChemMolecule) -> Vec<Vec<(usize, usize, i32)>> {
    let n = mol.atoms.len();
    let mut adj = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.begin_atom < n && b.end_atom < n {
            adj[b.begin_atom].push((b.end_atom, bi, b.order));
            adj[b.end_atom].push((b.begin_atom, bi, b.order));
        }
    }
    adj
}

// ============================================================================
// 1. Physiological pH Microspecies & Isoelectric Point (pI)
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub struct IonizableSite {
    pub atom_idx: usize,
    pub label: &'static str,
    pub pka: f64,
    pub is_acid: bool, // true: 0 -> -1 above pKa; false (base): +1 -> 0 above pKa
}

#[derive(Debug, Clone, PartialEq)]
pub struct MicrospeciesResult {
    pub ph: f64,
    pub net_charge: f64,
    pub formal_charge_sum: i32,
    pub isoelectric_point: Option<f64>,
    pub sites: Vec<IonizableSite>,
    pub molecule: ChemMolecule,
}

/// Identifies common organic acid/base ionizable sites with empirical macro-pKa values.
pub fn detect_ionizable_sites(mol: &ChemMolecule) -> Vec<IonizableSite> {
    let n = mol.atoms.len();
    let adj = build_adj(mol);
    let val = CcVal::compute(mol);
    let arom = CcAromaticityInfo::perceive(mol);
    let mut sites = Vec::new();

    let mut seen_guanidine_c = BTreeSet::new();
    for i in 0..n {
        let a = &mol.atoms[i];
        let z = a.z;
        let heavy_nbrs: Vec<(usize, i32)> = adj[i]
            .iter()
            .filter(|&&(nb, _, _)| mol.atoms[nb].z > 1)
            .map(|&(nb, _, o)| (nb, o))
            .collect();
        let imp_h = val.implicit_h_counts.get(i).copied().unwrap_or(0) as usize;
        let exp_h = adj[i].iter().filter(|&&(nb, _, _)| mol.atoms[nb].z == 1).count();
        let total_h = imp_h + exp_h;

        if z == 8 && heavy_nbrs.len() == 1 && heavy_nbrs[0].1 == 1 {
            let (nb, _) = heavy_nbrs[0];
            if mol.atoms[nb].z == 6 {
                let has_co = adj[nb].iter().any(|&(k, _, o)| k != i && mol.atoms[k].z == 8 && o == 2);
                if has_co {
                    sites.push(IonizableSite { atom_idx: i, label: "Carboxylic Acid", pka: 4.2, is_acid: true });
                } else if arom.is_atom_aromatic(nb) {
                    sites.push(IonizableSite { atom_idx: i, label: "Phenol", pka: 9.9, is_acid: true });
                }
            } else if matches!(mol.atoms[nb].z, 15 | 16) {
                sites.push(IonizableSite { atom_idx: i, label: "Sulfonic / Phosphoric Acid", pka: 1.9, is_acid: true });
            }
        } else if z == 16 && heavy_nbrs.len() == 1 && total_h >= 1 {
            sites.push(IonizableSite { atom_idx: i, label: "Thiol", pka: 8.5, is_acid: true });
        } else if z == 7 {
            let is_amide = heavy_nbrs.iter().any(|&(nb, _)| {
                mol.atoms[nb].z == 6 && adj[nb].iter().any(|&(k, _, o)| mol.atoms[k].z == 8 && o == 2)
            });
            let is_nitro = heavy_nbrs.iter().filter(|&&(nb, _)| mol.atoms[nb].z == 8).count() >= 2;
            let is_nitrile = heavy_nbrs.iter().any(|&(_, o)| o == 3);
            if is_amide || is_nitro || is_nitrile {
                continue;
            }
            let is_sulfonamide = heavy_nbrs.iter().any(|&(nb, _)| {
                mol.atoms[nb].z == 16 && adj[nb].iter().any(|&(k, _, o)| mol.atoms[k].z == 8 && o == 2)
            });
            if is_sulfonamide && total_h >= 1 {
                sites.push(IonizableSite { atom_idx: i, label: "Sulfonamide", pka: 10.1, is_acid: true });
                continue;
            }
            if arom.is_atom_aromatic(i) {
                if heavy_nbrs.len() == 2 && total_h == 0 {
                    sites.push(IonizableSite { atom_idx: i, label: "Pyridine / Imidazole-like N", pka: 5.3, is_acid: false });
                }
                continue;
            }
            let attached_arom = heavy_nbrs.iter().any(|&(nb, _)| arom.is_atom_aromatic(nb));
            let guanidine_c = heavy_nbrs.iter().find_map(|&(nb, _)| {
                (mol.atoms[nb].z == 6 && adj[nb].iter().filter(|&&(k, _, _)| mol.atoms[k].z == 7).count() >= 3)
                    .then_some(nb)
            });
            if let Some(gc) = guanidine_c {
                if seen_guanidine_c.insert(gc) {
                    sites.push(IonizableSite { atom_idx: i, label: "Guanidine", pka: 12.5, is_acid: false });
                }
            } else if attached_arom {
                sites.push(IonizableSite { atom_idx: i, label: "Aniline", pka: 4.6, is_acid: false });
            } else if heavy_nbrs.len() <= 3 && a.charge <= 0 {
                sites.push(IonizableSite { atom_idx: i, label: "Aliphatic Amine", pka: 9.8, is_acid: false });
            }
        }
    }
    sites
}

pub fn net_charge_at_ph(sites: &[IonizableSite], ph: f64) -> f64 {
    sites
        .iter()
        .map(|s| {
            if s.is_acid {
                -1.0 / (1.0 + 10.0_f64.powf(s.pka - ph))
            } else {
                1.0 / (1.0 + 10.0_f64.powf(ph - s.pka))
            }
        })
        .sum()
}

/// Computes the dominant protonation microspecies at `ph` (`1.0..=14.0`), fractional net charge, and isoelectric point `pI`.
pub fn compute_microspecies_at_ph(mol: &ChemMolecule, ph: f64) -> MicrospeciesResult {
    let sites = detect_ionizable_sites(mol);
    let mut out = mol.clone();
    let site_atoms: BTreeSet<usize> = sites.iter().map(|s| s.atom_idx).collect();
    let perm_charge: f64 = out
        .atoms
        .iter()
        .enumerate()
        .filter(|(i, _)| !site_atoms.contains(i))
        .map(|(_, a)| a.charge as f64)
        .sum();
    for s in &sites {
        out.atoms[s.atom_idx].charge = if s.is_acid {
            if ph > s.pka { -1 } else { 0 }
        } else if ph < s.pka {
            1
        } else {
            0
        };
    }
    let net_charge = perm_charge + net_charge_at_ph(&sites, ph);
    let formal_charge_sum = out.atoms.iter().map(|a| a.charge).sum();

    // Solve pI (pH where net_charge == 0) via bisection if charge crosses zero across 0..=14
    let q_lo = perm_charge + net_charge_at_ph(&sites, 0.0);
    let q_hi = perm_charge + net_charge_at_ph(&sites, 14.0);
    let isoelectric_point = if q_lo >= 0.0 && q_hi <= 0.0 && !sites.is_empty() {
        let mut lo = 0.0_f64;
        let mut hi = 14.0_f64;
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if perm_charge + net_charge_at_ph(&sites, mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(0.5 * (lo + hi))
    } else {
        None
    };

    MicrospeciesResult {
        ph,
        net_charge,
        formal_charge_sum,
        isoelectric_point,
        sites,
        molecule: out,
    }
}

// ============================================================================
// 2. Tautomer Enumerator (1,3-Prototropic Shifts: Keto-Enol, Amide-Imidic, Imine-Enamine)
// ============================================================================

fn bond_signature(mol: &ChemMolecule) -> String {
    let mut parts: Vec<String> = mol
        .bonds
        .iter()
        .map(|b| {
            let (u, v) = if b.begin_atom < b.end_atom {
                (b.begin_atom, b.end_atom)
            } else {
                (b.end_atom, b.begin_atom)
            };
            format!("{u}-{v}:{}:{:?}", b.order, b.display)
        })
        .collect();
    parts.sort();
    let chgs: Vec<String> = mol.atoms.iter().map(|a| format!("{}:{}", a.z, a.charge)).collect();
    format!("{}|{}", chgs.join(","), parts.join(";"))
}

/// Enumerates valid prototropic tautomers (keto–enol, imine–enamine, amide/lactam–imidic acid,
/// and 1,3-heterocyclic H-shifts).
pub fn enumerate_tautomers(mol: &ChemMolecule, max_count: usize) -> Vec<(String, ChemMolecule)> {
    let n = mol.atoms.len();
    let adj = build_adj(mol);
    let val = CcVal::compute(mol);
    let mut results = vec![("Input tautomer".to_string(), mol.clone())];
    let mut seen = BTreeSet::new();
    seen.insert(bond_signature(mol));

    // Match 1,3-triad: X(=)Y(-)Z-H  -->  H-X(-)Y(=)Z
    for y in 0..n {
        if !matches!(mol.atoms[y].z, 6 | 7) {
            continue;
        }
        for &(x, bx, ox) in &adj[y] {
            if ox != 2 || !matches!(mol.atoms[x].z, 6 | 7 | 8 | 16) {
                continue;
            }
            for &(z, bz, oz) in &adj[y] {
                if z == x || oz != 1 || !matches!(mol.atoms[z].z, 6 | 7 | 8 | 16) {
                    continue;
                }
                let zh = val.implicit_h_counts.get(z).copied().unwrap_or(0);
                if zh < 1 {
                    continue;
                }
                // Avoid breaking a benzene-like aromatic C-C unless heteroatom involved
                if mol.atoms[x].z == 6 && mol.atoms[z].z == 6 {
                    continue;
                }
                let mut cand = mol.clone();
                cand.bonds[bx].order = 1;
                cand.bonds[bz].order = 2;
                let sig = bond_signature(&cand);
                if seen.insert(sig) {
                    let desc = match (mol.atoms[x].z, mol.atoms[z].z) {
                        (8, 6) => "Enol tautomer",
                        (6, 8) => "Keto tautomer",
                        (8, 7) => "Imidic acid / Lactim tautomer",
                        (7, 8) => "Amide / Lactam tautomer",
                        (7, 6) => "Enamine tautomer",
                        (6, 7) => "Imine tautomer",
                        _ => "1,3-Prototropic tautomer",
                    };
                    results.push((desc.to_string(), cand));
                    if results.len() >= max_count.max(1) {
                        return results;
                    }
                }
            }
        }
    }
    results
}

// ============================================================================
// 3. Stereoisomer Enumerator (2^n Tetrahedral Centers)
// ============================================================================

/// Enumerates tetrahedral stereoisomers by toggling wedged/hashed stereo bonds at perceived
/// stereocenters (capped at `max_count`, typically 8–16).
pub fn enumerate_stereoisomers(mol: &ChemMolecule, max_count: usize) -> Vec<(String, ChemMolecule)> {
    let stereo = perceive_stereochemistry(mol);
    let n = mol.atoms.len();
    let adj = build_adj(mol);

    // Collect one representative exocyclic/terminal bond per tetrahedral center
    let mut center_bonds: Vec<(usize, usize)> = Vec::new();
    for (&ai, &code) in &stereo.atom_stereo {
        if code >= 2 && ai < n {
            // Prefer an existing stereo bond, otherwise first single bond from center
            let chosen = adj[ai]
                .iter()
                .find(|&&(_, bi, _)| !matches!(mol.bonds[bi].display, BondDisplay::Solid))
                .or_else(|| adj[ai].iter().find(|&&(_, _, o)| o == 1))
                .map(|&(_, bi, _)| bi);
            if let Some(bi) = chosen {
                if !center_bonds.iter().any(|&(_, b)| b == bi) {
                    center_bonds.push((ai, bi));
                }
            }
        }
    }
    // Fallback: also include any explicitly wedged/hashed bond not yet in `center_bonds`
    for (bi, b) in mol.bonds.iter().enumerate() {
        if matches!(
            b.display,
            BondDisplay::WedgeBegin | BondDisplay::WedgeEnd | BondDisplay::WedgedHashBegin | BondDisplay::WedgedHashEnd
        ) && !center_bonds.iter().any(|&(_, existing)| existing == bi)
        {
            center_bonds.push((b.begin_atom, bi));
        }
    }

    if center_bonds.is_empty() {
        return vec![("Achiral / No stereocenters detected".to_string(), mol.clone())];
    }

    let k = center_bonds.len().min(4); // up to 16 stereoisomers
    let total = (1usize << k).min(max_count.max(1));
    let mut out = Vec::with_capacity(total);
    for mask in 0..total {
        let mut iso = mol.clone();
        let mut tags = Vec::new();
        for bit in 0..k {
            let (ai, bi) = center_bonds[bit];
            let is_wedge = ((mask >> bit) & 1) == 0;
            let b = &mut iso.bonds[bi];
            if b.begin_atom == ai {
                b.display = if is_wedge { BondDisplay::WedgeBegin } else { BondDisplay::WedgedHashBegin };
            } else {
                b.display = if is_wedge { BondDisplay::WedgeEnd } else { BondDisplay::WedgedHashEnd };
            }
            tags.push(format!("C{}:{}", ai + 1, if is_wedge { "▲" } else { "┄" }));
        }
        let st = perceive_stereochemistry(&iso);
        let rs_tags: Vec<String> = center_bonds
            .iter()
            .take(k)
            .filter_map(|&(ai, _)| match st.atom_stereo.get(&ai).copied() {
                Some(2) => Some(format!("{}(R)", ai + 1)),
                Some(3) => Some(format!("{}(S)", ai + 1)),
                _ => None,
            })
            .collect();
        let label = if !rs_tags.is_empty() {
            format!("Isomer {}: {}", mask + 1, rs_tags.join(", "))
        } else {
            format!("Isomer {}: {}", mask + 1, tags.join(" "))
        };
        out.push((label, iso));
    }
    out
}

// ============================================================================
// 4. Resonance Structure Generator & Curved-Arrow Electron Pusher
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub struct ResonanceContributor {
    pub description: String,
    pub weight_score: i32,
    pub molecule: ChemMolecule,
}

/// Enumerates valid resonance contributors (polarized carbonyl/imine/nitro forms, lone-pair donation
/// into adjacent π-bonds, carboxylate delocalization, and alternate Kekulé aromatic rings).
pub fn enumerate_resonance_structures(mol: &ChemMolecule, max_count: usize) -> Vec<ResonanceContributor> {
    let n = mol.atoms.len();
    let adj = build_adj(mol);
    let mut seen = BTreeSet::new();
    seen.insert(bond_signature(mol));

    let score_mol = |m: &ChemMolecule| -> i32 {
        let chg_pen: i32 = m.atoms.iter().map(|a| a.charge.abs() * 15).sum();
        let neg_on_en: i32 = m
            .atoms
            .iter()
            .map(|a| if a.charge < 0 && matches!(a.z, 7 | 8 | 9 | 16) { 6 } else { 0 })
            .sum();
        100 - chg_pen + neg_on_en
    };

    let mut out = vec![ResonanceContributor {
        description: "Major / Input contributor".to_string(),
        weight_score: score_mol(mol),
        molecule: mol.clone(),
    }];

    // 1. Polarized heteronuclear double bond: X=Y -> X(+)-Y(-) (e.g., C=O -> C(+)-O(-))
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.order == 2 && b.begin_atom < n && b.end_atom < n {
            let (u, v) = (b.begin_atom, b.end_atom);
            let pairs = [(u, v), (v, u)];
            for (donor, acc) in pairs {
                if mol.atoms[donor].z == 6 && matches!(mol.atoms[acc].z, 7 | 8 | 16) {
                    let mut cand = mol.clone();
                    cand.bonds[bi].order = 1;
                    cand.atoms[donor].charge += 1;
                    cand.atoms[acc].charge -= 1;
                    if seen.insert(bond_signature(&cand)) {
                        out.push(ResonanceContributor {
                            description: format!("Polarized π-bond (atom #{}⁺ — #{}⁻)", donor + 1, acc + 1),
                            weight_score: score_mol(&cand),
                            molecule: cand,
                        });
                    }
                }
            }
        }
    }

    // 2. Allylic / Amide / Carboxylate 3-atom delocalization: D(-/lp)-Y=A -> D(+/0)=Y-A(-)
    for y in 0..n {
        for &(acc, ba, oa) in &adj[y] {
            if oa != 2 {
                continue;
            }
            for &(don, bd, od) in &adj[y] {
                if don == acc || od != 1 {
                    continue;
                }
                let don_no_mult = adj[don].iter().all(|&(_, _, o)| o == 1);
                let can_donate = mol.atoms[don].charge < 0
                    || (matches!(mol.atoms[don].z, 7 | 8 | 16) && mol.atoms[don].charge == 0 && don_no_mult);
                if can_donate {
                    let mut cand = mol.clone();
                    cand.bonds[bd].order = 2;
                    cand.bonds[ba].order = 1;
                    cand.atoms[don].charge += 1;
                    cand.atoms[acc].charge -= 1;
                    if seen.insert(bond_signature(&cand)) {
                        let desc = if mol.atoms[don].charge < 0 {
                            format!("Anionic delocalization (#{} → #{})", don + 1, acc + 1)
                        } else {
                            format!("Lone-pair conjugation (#{} → #{})", don + 1, acc + 1)
                        };
                        out.push(ResonanceContributor {
                            description: desc,
                            weight_score: score_mol(&cand),
                            molecule: cand,
                        });
                    }
                }
            }
        }
    }

    // 3. Alternate Kekulé aromatic ring resonance form (constrained valid re-kekulization so fused rings & heterocycles stay valence-valid)
    let arom = CcAromaticityInfo::perceive(mol);
    let arom_dbl_bonds: Vec<usize> = mol
        .bonds
        .iter()
        .enumerate()
        .filter(|&(bi, b)| arom.is_bond_aromatic(bi) && b.order == 2)
        .map(|(bi, _)| bi)
        .collect();
    if arom_dbl_bonds.len() >= 3 {
        let mut base_dek = mol.clone();
        crate::dekekulize_aromatic_bonds(&mut base_dek);
        for &ban_bi in &arom_dbl_bonds {
            let mut kek = base_dek.clone();
            // Force one previously-double aromatic bond to single, then find a valid Kekulé matching
            kek.bonds[ban_bi].order = 1;
            if crate::kekulize_aromatic_bonds(&mut kek) && seen.insert(bond_signature(&kek)) {
                out.push(ResonanceContributor {
                    description: "Alternate Kekulé aromatic contributor".to_string(),
                    weight_score: score_mol(&kek),
                    molecule: kek,
                });
                break;
            }
        }
    }

    out.sort_by(|a, b| b.weight_score.cmp(&a.weight_score));
    out.truncate(max_count.max(1));
    out
}

/// Applies a curved electron-pushing arrow (`from_atom -> to_atom`) and returns the resulting
/// mechanistic intermediate / resonance form with updated bond orders and formal charges!
pub fn apply_electron_push_arrow(
    mol: &ChemMolecule,
    from_atom: usize,
    to_atom: usize,
) -> Result<(ChemMolecule, String), String> {
    let n = mol.atoms.len();
    if from_atom >= n || to_atom >= n || from_atom == to_atom {
        return Err("Select two distinct atoms for an electron-pushing arrow".to_string());
    }
    let mut out = mol.clone();
    // Case A: Direct bond between from_atom and to_atom
    if let Some(bi) = out.bonds.iter().position(|b| {
        (b.begin_atom == from_atom && b.end_atom == to_atom)
            || (b.begin_atom == to_atom && b.end_atom == from_atom)
    }) {
        if out.bonds[bi].order >= 2 {
            // Heterolytic π-bond shift onto to_atom: A=B -> A(+)-B(-)
            out.bonds[bi].order -= 1;
            out.atoms[from_atom].charge += 1;
            out.atoms[to_atom].charge -= 1;
            return Ok((
                out,
                format!("Pushed π-electrons from #{} onto #{} (bond order {}, formal charges updated)", from_atom + 1, to_atom + 1, mol.bonds[bi].order - 1),
            ));
        } else {
            // Lone pair from `from_atom` forms a π-bond with `to_atom`: A-B -> A(+)=B(-)
            // If `to_atom` has an adjacent π-bond, push that π-bond out to maintain octet!
            out.bonds[bi].order += 1;
            out.atoms[from_atom].charge += 1;
            let adj_pi = out.bonds.iter().enumerate().find(|&(k, b)| {
                k != bi && b.order >= 2 && (b.begin_atom == to_atom || b.end_atom == to_atom)
            }).map(|(k, b)| (k, if b.begin_atom == to_atom { b.end_atom } else { b.begin_atom }));
            if let Some((pi_idx, sink)) = adj_pi {
                out.bonds[pi_idx].order -= 1;
                out.atoms[sink].charge -= 1;
                return Ok((
                    out,
                    format!("Conjugated push: lone pair #{} → π(#{}-#{}) → atom #{}⁻", from_atom + 1, from_atom + 1, to_atom + 1, sink + 1),
                ));
            }
            out.atoms[to_atom].charge -= 1;
            return Ok((
                out,
                format!("Formed π-bond between #{} and #{} from lone pair", from_atom + 1, to_atom + 1),
            ));
        }
    }

    // Case B: Nucleophilic attack between two non-bonded atoms (`from_atom` attacks `to_atom`)
    let next_id = out.bonds.iter().map(|b| b.id).max().unwrap_or(100) + 1;
    out.add_bond(next_id, from_atom, to_atom, 1, BondDisplay::Solid);
    out.atoms[from_atom].charge += 1;
    // If electrophile `to_atom` has a π-bond (e.g. C=O), break it onto its partner!
    let adj_pi = out.bonds.iter().enumerate().find(|&(k, b)| {
        k + 1 < out.bonds.len() && b.order >= 2 && (b.begin_atom == to_atom || b.end_atom == to_atom)
    }).map(|(k, b)| (k, if b.begin_atom == to_atom { b.end_atom } else { b.begin_atom }));
    if let Some((pi_idx, sink)) = adj_pi {
        out.bonds[pi_idx].order -= 1;
        out.atoms[sink].charge -= 1;
        return Ok((
            out,
            format!("Nucleophilic addition: #{} attacked #{}, π-electrons shifted to #{}⁻", from_atom + 1, to_atom + 1, sink + 1),
        ));
    }
    out.atoms[to_atom].charge -= 1;
    Ok((
        out,
        format!("Formed σ-bond via electron push #{} → #{}", from_atom + 1, to_atom + 1),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ChemAtom;

    #[test]
    fn test_glycine_zwitterion_ph_and_pi() {
        // Glycine: H2N-CH2-C(=O)-OH
        let mut gly = ChemMolecule::new();
        gly.add_atom(ChemAtom::new(1, 7)); // N
        gly.add_atom(ChemAtom::new(2, 6)); // CH2
        gly.add_atom(ChemAtom::new(3, 6)); // C=O
        gly.add_atom(ChemAtom::new(4, 8)); // =O
        gly.add_atom(ChemAtom::new(5, 8)); // -OH
        gly.add_bond(10, 0, 1, 1, BondDisplay::Solid);
        gly.add_bond(11, 1, 2, 1, BondDisplay::Solid);
        gly.add_bond(12, 2, 3, 2, BondDisplay::Solid);
        gly.add_bond(13, 2, 4, 1, BondDisplay::Solid);

        let phys = compute_microspecies_at_ph(&gly, 7.4);
        assert_eq!(phys.molecule.atoms[0].charge, 1, "Amine protonated at pH 7.4");
        assert_eq!(phys.molecule.atoms[4].charge, -1, "Carboxylate deprotonated at pH 7.4");
        assert_eq!(phys.formal_charge_sum, 0);
        let pi = phys.isoelectric_point.expect("amphoteric pI");
        assert!((pi - 7.0).abs() < 1.2, "Expected Glycine pI ~6-7, got {pi:.2}");
    }

    #[test]
    fn test_tautomer_resonance_and_electron_push() {
        // Acetone: CH3-C(=O)-CH3
        let mut ace = ChemMolecule::new();
        ace.add_atom(ChemAtom::new(1, 6));
        ace.add_atom(ChemAtom::new(2, 6));
        ace.add_atom(ChemAtom::new(3, 8));
        ace.add_atom(ChemAtom::new(4, 6));
        ace.add_bond(10, 0, 1, 1, BondDisplay::Solid);
        ace.add_bond(11, 1, 2, 2, BondDisplay::Solid);
        ace.add_bond(12, 1, 3, 1, BondDisplay::Solid);

        let tauts = enumerate_tautomers(&ace, 6);
        assert!(tauts.len() >= 2, "Acetone must enumerate its enol tautomer");
        let res = enumerate_resonance_structures(&ace, 6);
        assert!(res.len() >= 2, "Acetone must enumerate C(+)-O(-) resonance contributor");
        let (pushed, msg) = apply_electron_push_arrow(&ace, 1, 2).expect("push C=O");
        assert_eq!(pushed.bonds[1].order, 1);
        assert_eq!(pushed.atoms[1].charge, 1);
        assert_eq!(pushed.atoms[2].charge, -1);
        assert!(msg.contains("π-electrons"));
    }
}
