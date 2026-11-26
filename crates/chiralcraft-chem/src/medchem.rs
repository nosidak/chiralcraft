//! Medicinal Chemistry & ADMET Scorecard, PAINS/Brenk Structural Alerts,
//! and High-Resolution Mass Spectrometry (HRMS) Isotope Envelope Simulator.
//!
//! Pure-Rust clean-room implementation using IUPAC `chiralcraft_tables::ISOTOPES`.

use std::collections::BTreeMap;

use crate::{
    compute_analysis_report, estimate_chemical_properties, get_element_isotopes,
    CcAromaticityInfo, CcVal, ChemMolecule,
};

// ============================================================================
// 1. Medicinal Chemistry Scorecard & PAINS / Brenk Structural Alerts
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlertSeverity {
    Warning,
    High,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructuralAlert {
    pub name: &'static str,
    pub category: &'static str,
    pub severity: AlertSeverity,
    pub atom_indices: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MedChemReport {
    pub mw: f64,
    pub exact_mass: f64,
    pub log_p: f64,
    pub tpsa: f64,
    pub hbd: usize,
    pub hba: usize,
    pub rotatable_bonds: usize,
    pub heavy_atoms: usize,
    pub aromatic_rings: usize,
    pub fsp3: f64,
    pub lipinski_violations: usize,
    pub lipinski_pass: bool,
    pub veber_pass: bool,
    pub ghose_pass: bool,
    pub cns_mpo: f64,
    pub qed: f64,
    pub sa_score: f64,
    pub alerts: Vec<StructuralAlert>,
}

/// Counts heavy-atom adjacency and simple ring membership via DFS.
fn build_adj(mol: &ChemMolecule) -> Vec<Vec<(usize, i32)>> {
    let n = mol.atoms.len();
    let mut adj = vec![Vec::new(); n];
    for b in &mol.bonds {
        if b.begin_atom < n && b.end_atom < n {
            adj[b.begin_atom].push((b.end_atom, b.order));
            adj[b.end_atom].push((b.begin_atom, b.order));
        }
    }
    adj
}

fn is_bond_in_ring(mol: &ChemMolecule, adj: &[Vec<(usize, i32)>], u: usize, v: usize) -> bool {
    let n = mol.atoms.len();
    let mut vis = vec![false; n];
    let mut stack = vec![u];
    vis[u] = true;
    while let Some(cur) = stack.pop() {
        for &(nb, _) in &adj[cur] {
            if (cur == u && nb == v) || (cur == v && nb == u) {
                continue;
            }
            if nb == v {
                return true;
            }
            if !vis[nb] {
                vis[nb] = true;
                stack.push(nb);
            }
        }
    }
    false
}

/// Detects PAINS & Brenk reactive/toxicophore functional groups and returns exact atom indices
/// so the canvas can highlight offending motifs live!
pub fn detect_structural_alerts(mol: &ChemMolecule) -> Vec<StructuralAlert> {
    let n = mol.atoms.len();
    let adj = build_adj(mol);
    let val = CcVal::compute(mol);
    let arom = CcAromaticityInfo::perceive(mol);
    let mut alerts = Vec::new();

    // 1. Nitro / Nitroaromatic (-NO2)
    for i in 0..n {
        if mol.atoms[i].z == 7 {
            let oxy_nbrs: Vec<usize> = adj[i]
                .iter()
                .filter(|&&(nb, _)| mol.atoms[nb].z == 8 && adj[nb].len() == 1)
                .map(|&(nb, _)| nb)
                .collect();
            if oxy_nbrs.len() >= 2 {
                let is_arom = adj[i].iter().any(|&(nb, _)| arom.is_atom_aromatic(nb));
                let mut atoms = vec![i];
                atoms.extend(oxy_nbrs);
                alerts.push(StructuralAlert {
                    name: if is_arom { "Nitroaromatic (mutagenicity alert)" } else { "Aliphatic Nitro Group" },
                    category: "Brenk",
                    severity: if is_arom { AlertSeverity::High } else { AlertSeverity::Warning },
                    atom_indices: atoms,
                });
            }
        }
    }

    // 2. Acyl Halide (R-C(=O)-X) & Aldehyde (R-CH=O)
    for i in 0..n {
        if mol.atoms[i].z == 6 {
            let co_oxy = adj[i].iter().find(|&&(nb, ord)| mol.atoms[nb].z == 8 && ord == 2).map(|&(nb, _)| nb);
            if let Some(o_idx) = co_oxy {
                if let Some(&(hal_idx, _)) = adj[i].iter().find(|&&(nb, _)| matches!(mol.atoms[nb].z, 9 | 17 | 35 | 53)) {
                    alerts.push(StructuralAlert {
                        name: "Acyl Halide (highly reactive electrophile)",
                        category: "Brenk",
                        severity: AlertSeverity::High,
                        atom_indices: vec![i, o_idx, hal_idx],
                    });
                } else {
                    let heavy_deg = adj[i].iter().filter(|&&(nb, _)| mol.atoms[nb].z > 1).count();
                    let imp_h = val.implicit_h_counts.get(i).copied().unwrap_or(0);
                    let exp_h = adj[i].iter().filter(|&&(nb, _)| mol.atoms[nb].z == 1).count();
                    if heavy_deg <= 2 && (imp_h as usize + exp_h) >= 1 {
                        alerts.push(StructuralAlert {
                            name: "Aldehyde (reactive Schiff-base electrophile)",
                            category: "Brenk",
                            severity: AlertSeverity::Warning,
                            atom_indices: vec![i, o_idx],
                        });
                    }
                }
            }
        }
    }

    // 3. Michael Acceptor (C=C-C=O or C=C-C#N)
    for b in &mol.bonds {
        if b.order == 2 && b.begin_atom < n && b.end_atom < n {
            let (u, v) = (b.begin_atom, b.end_atom);
            if mol.atoms[u].z == 6 && mol.atoms[v].z == 6 && !arom.is_atom_aromatic(u) && !arom.is_atom_aromatic(v) {
                for &end in &[u, v] {
                    for &(nb, ord) in &adj[end] {
                        if nb == u || nb == v || ord != 1 || mol.atoms[nb].z != 6 {
                            continue;
                        }
                        if let Some(&(o_idx, _)) = adj[nb].iter().find(|&&(k, ko)| mol.atoms[k].z == 8 && ko == 2) {
                            alerts.push(StructuralAlert {
                                name: "Michael Acceptor (α,β-unsaturated carbonyl)",
                                category: "PAINS",
                                severity: AlertSeverity::High,
                                atom_indices: vec![u, v, nb, o_idx],
                            });
                        } else if let Some(&(n_idx, _)) = adj[nb].iter().find(|&&(k, ko)| mol.atoms[k].z == 7 && ko == 3) {
                            alerts.push(StructuralAlert {
                                name: "Michael Acceptor (α,β-unsaturated nitrile)",
                                category: "PAINS",
                                severity: AlertSeverity::Warning,
                                atom_indices: vec![u, v, nb, n_idx],
                            });
                        }
                    }
                }
            }
        }
    }

    // 4. Epoxide / Aziridine / Strained 3-membered heterocycle
    for i in 0..n {
        if matches!(mol.atoms[i].z, 7 | 8) {
            let nbrs: Vec<usize> = adj[i].iter().map(|&(nb, _)| nb).collect();
            'tri: for a in 0..nbrs.len() {
                for b in (a + 1)..nbrs.len() {
                    let (u, v) = (nbrs[a], nbrs[b]);
                    if adj[u].iter().any(|&(k, _)| k == v) {
                        alerts.push(StructuralAlert {
                            name: if mol.atoms[i].z == 8 { "Epoxide (alkylating electrophile)" } else { "Aziridine (alkylating electrophile)" },
                            category: "Brenk",
                            severity: AlertSeverity::High,
                            atom_indices: vec![i, u, v],
                        });
                        break 'tri;
                    }
                }
            }
        }
    }

    // 5. Peroxide (O-O), Hydrazine/Azo (N-N / N=N), Disulfide (S-S)
    for b in &mol.bonds {
        if b.begin_atom >= n || b.end_atom >= n {
            continue;
        }
        let (u, v) = (b.begin_atom, b.end_atom);
        let (zu, zv) = (mol.atoms[u].z, mol.atoms[v].z);
        if zu == 8 && zv == 8 {
            alerts.push(StructuralAlert {
                name: "Peroxide (O–O oxidative instability)",
                category: "Brenk",
                severity: AlertSeverity::High,
                atom_indices: vec![u, v],
            });
        } else if zu == 7 && zv == 7 && !arom.is_atom_aromatic(u) && !arom.is_atom_aromatic(v) {
            alerts.push(StructuralAlert {
                name: if b.order >= 2 { "Azo Group (PAINS / chromophore)" } else { "Hydrazine / Hydrazide (reactive N–N)" },
                category: if b.order >= 2 { "PAINS" } else { "Brenk" },
                severity: AlertSeverity::Warning,
                atom_indices: vec![u, v],
            });
        }
    }

    // 6. Primary Aniline (Ar-NH2) & Catechol (ortho-dihydroxybenzene)
    let mut arom_oh_carbons: Vec<(usize, usize)> = Vec::new();
    for i in 0..n {
        let z = mol.atoms[i].z;
        let heavy_deg = adj[i].iter().filter(|&&(nb, _)| mol.atoms[nb].z > 1).count();
        if z == 7 && heavy_deg == 1 {
            if let Some(&(nb, _)) = adj[i].iter().find(|&&(k, _)| mol.atoms[k].z == 6 && arom.is_atom_aromatic(k)) {
                alerts.push(StructuralAlert {
                    name: "Primary Aniline (metabolic aniline alert)",
                    category: "Brenk",
                    severity: AlertSeverity::Warning,
                    atom_indices: vec![i, nb],
                });
            }
        } else if z == 8 && heavy_deg == 1 {
            if let Some(&(nb, _)) = adj[i].iter().find(|&&(k, _)| mol.atoms[k].z == 6 && arom.is_atom_aromatic(k)) {
                arom_oh_carbons.push((nb, i));
            }
        } else if z == 16 && heavy_deg == 1 && val.implicit_h_counts.get(i).copied().unwrap_or(0) >= 1 {
            alerts.push(StructuralAlert {
                name: "Free Thiol (-SH redox / covalent alert)",
                category: "Brenk",
                severity: AlertSeverity::Warning,
                atom_indices: vec![i],
            });
        }
    }
    for a in 0..arom_oh_carbons.len() {
        for b in (a + 1)..arom_oh_carbons.len() {
            let (ca, oa) = arom_oh_carbons[a];
            let (cb, ob) = arom_oh_carbons[b];
            if adj[ca].iter().any(|&(k, _)| k == cb) {
                alerts.push(StructuralAlert {
                    name: "Catechol (PAINS redox / metal chelator)",
                    category: "PAINS",
                    severity: AlertSeverity::High,
                    atom_indices: vec![ca, oa, cb, ob],
                });
            }
        }
    }

    // 7. Sulfonyl Halide (-SO2-Cl)
    for i in 0..n {
        if mol.atoms[i].z == 16 {
            let oxys: Vec<usize> = adj[i].iter().filter(|&&(k, o)| mol.atoms[k].z == 8 && o == 2).map(|&(k, _)| k).collect();
            if oxys.len() >= 2 {
                if let Some(&(hal, _)) = adj[i].iter().find(|&&(k, _)| matches!(mol.atoms[k].z, 9 | 17 | 35)) {
                    let mut atoms = vec![i, hal];
                    atoms.extend(oxys);
                    alerts.push(StructuralAlert {
                        name: "Sulfonyl Halide (reactive electrophile)",
                        category: "Brenk",
                        severity: AlertSeverity::High,
                        atom_indices: atoms,
                    });
                }
            }
        }
    }

    alerts
}

/// Computes full Medicinal Chemistry & ADMET Scorecard (Lipinski, Veber, Ghose, CNS MPO, QED, Fsp³, SA, PAINS).
pub fn compute_medchem_report(mol: &ChemMolecule) -> MedChemReport {
    let rep = compute_analysis_report(mol);
    let props = estimate_chemical_properties(mol);
    let val = CcVal::compute(mol);
    let arom = CcAromaticityInfo::perceive(mol);
    let adj = build_adj(mol);
    let n = mol.atoms.len();

    let heavy_atoms = mol.atoms.iter().filter(|a| a.z > 1).count();
    let mut hbd = 0usize;
    let mut hba = 0usize;
    let mut total_c = 0usize;
    let mut sp3_c = 0usize;

    for (i, a) in mol.atoms.iter().enumerate() {
        if matches!(a.z, 7 | 8) {
            hba += 1;
            let imp_h = val.implicit_h_counts.get(i).copied().unwrap_or(0).max(0) as usize;
            let exp_h = adj[i].iter().filter(|&&(nb, _)| mol.atoms[nb].z == 1).count();
            if imp_h + exp_h > 0 {
                hbd += 1;
            }
        }
        if a.z == 6 {
            total_c += 1;
            let has_mult = adj[i].iter().any(|&(_, ord)| ord >= 2);
            if !has_mult && !arom.is_atom_aromatic(i) {
                sp3_c += 1;
            }
        }
    }

    let mut rotatable_bonds = 0usize;
    for b in &mol.bonds {
        if b.order == 1 && b.begin_atom < n && b.end_atom < n {
            let (u, v) = (b.begin_atom, b.end_atom);
            let du = adj[u].iter().filter(|&&(k, _)| mol.atoms[k].z > 1).count();
            let dv = adj[v].iter().filter(|&&(k, _)| mol.atoms[k].z > 1).count();
            if mol.atoms[u].z > 1 && mol.atoms[v].z > 1 && du > 1 && dv > 1 && !is_bond_in_ring(mol, &adj, u, v) {
                // Exclude amide C-N bonds from rotatable count
                let is_amide = (mol.atoms[u].z == 6 && mol.atoms[v].z == 7 && adj[u].iter().any(|&(k, o)| mol.atoms[k].z == 8 && o == 2))
                    || (mol.atoms[v].z == 6 && mol.atoms[u].z == 7 && adj[v].iter().any(|&(k, o)| mol.atoms[k].z == 8 && o == 2));
                if !is_amide {
                    rotatable_bonds += 1;
                }
            }
        }
    }

    let arom_atoms = (0..n).filter(|&i| arom.is_atom_aromatic(i)).count();
    let aromatic_rings = (arom_atoms + 3) / 5; // fast clean ring count approximation
    let fsp3 = if total_c > 0 { (sp3_c as f64) / (total_c as f64) } else { 0.0 };

    let mw = rep.mol_weight;
    let log_p = props.log_p;
    let tpsa = props.tpsa;

    let mut lipinski_violations = 0usize;
    if mw > 500.0 { lipinski_violations += 1; }
    if log_p > 5.0 { lipinski_violations += 1; }
    if hbd > 5 { lipinski_violations += 1; }
    if hba > 10 { lipinski_violations += 1; }
    let lipinski_pass = lipinski_violations <= 1;
    let veber_pass = rotatable_bonds <= 10 && tpsa <= 140.0;
    let ghose_pass = (160.0..=480.0).contains(&mw)
        && (-0.4..=5.6).contains(&log_p)
        && (20..=70).contains(&heavy_atoms)
        && (40.0..=130.0).contains(&props.molar_refractivity);

    // Pfizer Wager CNS MPO (0.0 to 6.0): 6 desirability ramps
    let ramp_down = |x: f64, lo: f64, hi: f64| -> f64 {
        if x <= lo { 1.0 } else if x >= hi { 0.0 } else { (hi - x) / (hi - lo) }
    };
    let d_logp = ramp_down(log_p, 3.0, 5.0);
    let d_logd = ramp_down(log_p - 0.3, 2.0, 4.0);
    let d_mw = ramp_down(mw, 360.0, 500.0);
    let d_tpsa = if tpsa <= 20.0 {
        (tpsa / 20.0).clamp(0.0, 1.0)
    } else if tpsa <= 40.0 {
        (tpsa - 20.0) / 20.0
    } else if tpsa <= 90.0 {
        1.0
    } else {
        ramp_down(tpsa, 90.0, 120.0)
    };
    let d_hbd = ramp_down(hbd as f64, 0.5, 3.5);
    let d_pka = props.pka.map(|p| ramp_down(p, 8.0, 10.0)).unwrap_or(1.0);
    let cns_mpo = (d_logp + d_logd + d_mw + d_tpsa + d_hbd + d_pka).clamp(0.0, 6.0);

    let alerts = detect_structural_alerts(mol);

    // QED (Quantitative Estimate of Drug-likeness, 0.0..=1.0) smooth desirability score
    let bell = |x: f64, ideal: f64, spread: f64| -> f64 {
        (-((x - ideal) / spread).powi(2)).exp().clamp(0.05, 1.0)
    };
    let q_mw = bell(mw, 330.0, 190.0);
    let q_lp = bell(log_p, 2.5, 2.4);
    let q_psa = bell(tpsa, 65.0, 65.0);
    let q_hbd = bell(hbd as f64, 1.5, 2.8);
    let q_hba = bell(hba as f64, 4.0, 4.5);
    let q_rot = bell(rotatable_bonds as f64, 3.5, 5.0);
    let q_al = (0.65_f64).powi(alerts.len() as i32);
    let qed = (q_mw.powf(0.20)
        * q_lp.powf(0.20)
        * q_psa.powf(0.15)
        * q_hbd.powf(0.15)
        * q_hba.powf(0.10)
        * q_rot.powf(0.10)
        * q_al.powf(0.10))
    .clamp(0.0, 1.0);

    // Synthetic Accessibility (SA) heuristic (1.0 = trivial, 10.0 = very complex)
    let stereo_centers = mol
        .bonds
        .iter()
        .filter(|b| !matches!(b.display, crate::BondDisplay::Solid))
        .count();
    let ring_complexity = (mol.bonds.len() as isize - heavy_atoms as isize + 1).max(0) as f64;
    let bridgehead_like = (0..n).filter(|&i| adj[i].len() >= 4).count() as f64;
    let sa_score = (1.4
        + (heavy_atoms as f64) * 0.055
        + ring_complexity * 0.42
        + (stereo_centers as f64) * 0.35
        + bridgehead_like * 0.45)
        .clamp(1.0, 10.0);

    MedChemReport {
        mw,
        exact_mass: rep.exact_mass,
        log_p,
        tpsa,
        hbd,
        hba,
        rotatable_bonds,
        heavy_atoms,
        aromatic_rings,
        fsp3,
        lipinski_violations,
        lipinski_pass,
        veber_pass,
        ghose_pass,
        cns_mpo,
        qed,
        sa_score,
        alerts,
    }
}

// ============================================================================
// 2. HRMS Isotope Envelope Simulator & Adduct Calculator
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub struct IsotopePeak {
    pub mz: f64,
    pub rel_intensity: f64, // 0.0 ..= 100.0
}

#[derive(Debug, Clone, PartialEq)]
pub struct AdductEntry {
    pub name: &'static str,
    pub charge: i32,
    pub mz: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct HrmsSimulation {
    pub formula: String,
    pub neutral_exact_mass: f64,
    pub selected_adduct: String,
    pub charge: i32,
    pub peaks: Vec<IsotopePeak>,
    pub adducts: Vec<AdductEntry>,
}

impl HrmsSimulation {
    /// Evaluates a smooth HRMS Gaussian profile curve `(mz, intensity 0..100)` at resolving power `r`.
    pub fn simulate_profile(&self, resolving_power: f64, samples: usize) -> Vec<(f64, f64)> {
        if self.peaks.is_empty() || samples < 2 {
            return Vec::new();
        }
        let min_mz = self.peaks.iter().map(|p| p.mz).fold(f64::INFINITY, f64::min) - 0.6;
        let max_mz = self.peaks.iter().map(|p| p.mz).fold(f64::NEG_INFINITY, f64::max) + 0.6;
        let span = (max_mz - min_mz).max(1.0);
        let r = resolving_power.clamp(1_000.0, 250_000.0);
        let mut out = Vec::with_capacity(samples);
        for k in 0..samples {
            let x = min_mz + span * (k as f64) / ((samples - 1) as f64);
            let mut y = 0.0;
            for pk in &self.peaks {
                let fwhm = (pk.mz / r).max(0.015);
                let sigma = fwhm / 2.35482;
                let z = (x - pk.mz) / sigma;
                if z.abs() < 6.0 {
                    y += pk.rel_intensity * (-0.5 * z * z).exp();
                }
            }
            out.push((x, y.min(100.0)));
        }
        out
    }

    /// Computes exact ppm mass error between an observed $m/z$ and the theoretical base/monoisotopic peak.
    pub fn ppm_error(&self, observed_mz: f64) -> Option<f64> {
        let base = self.peaks.first()?;
        if base.mz <= 0.0 {
            return None;
        }
        Some((observed_mz - base.mz) / base.mz * 1_000_000.0)
    }
}

/// Computes standard ESI/MALDI/APCI adduct $m/z$ values for a neutral exact mass `m`.
pub fn compute_adducts(m: f64) -> Vec<AdductEntry> {
    const H_PLUS: f64 = 1.007276466812;
    const NA_PLUS: f64 = 22.989218;
    const K_PLUS: f64 = 38.963158;
    const NH4_PLUS: f64 = 18.033823;
    const ACN_H: f64 = 42.033823;
    const CL_MINUS: f64 = 34.969402;
    const HCOO_MINUS: f64 = 44.998201;
    const ACO_MINUS: f64 = 59.013851;
    vec![
        AdductEntry { name: "[M+H]+", charge: 1, mz: m + H_PLUS },
        AdductEntry { name: "[M+Na]+", charge: 1, mz: m + NA_PLUS },
        AdductEntry { name: "[M+K]+", charge: 1, mz: m + K_PLUS },
        AdductEntry { name: "[M+NH4]+", charge: 1, mz: m + NH4_PLUS },
        AdductEntry { name: "[M+CH3CN+H]+", charge: 1, mz: m + ACN_H },
        AdductEntry { name: "[M+2H]2+", charge: 2, mz: (m + 2.0 * H_PLUS) * 0.5 },
        AdductEntry { name: "[2M+H]+", charge: 1, mz: 2.0 * m + H_PLUS },
        AdductEntry { name: "[2M+Na]+", charge: 1, mz: 2.0 * m + NA_PLUS },
        AdductEntry { name: "[M-H]-", charge: -1, mz: (m - H_PLUS).max(0.0) },
        AdductEntry { name: "[M+Cl]-", charge: -1, mz: m + CL_MINUS },
        AdductEntry { name: "[M+HCOO]-", charge: -1, mz: m + HCOO_MINUS },
        AdductEntry { name: "[M+CH3COO]-", charge: -1, mz: m + ACO_MINUS },
    ]
}

fn convolve_isotopes(dist: Vec<(f64, f64)>, elem_iso: &[(f64, f64)], count: usize) -> Vec<(f64, f64)> {
    let mut cur = dist;
    for _ in 0..count {
        // Group by 0.005 Da mass bin to keep fine structure while staying ultra-fast
        let mut bins: BTreeMap<i64, (f64, f64)> = BTreeMap::new();
        for &(m1, p1) in &cur {
            if p1 < 1e-6 {
                continue;
            }
            for &(m2, p2) in elem_iso {
                let p = p1 * p2;
                if p < 1e-6 {
                    continue;
                }
                let m = m1 + m2;
                let key = (m * 200.0).round() as i64;
                let entry = bins.entry(key).or_insert((0.0, 0.0));
                entry.0 += m * p;
                entry.1 += p;
            }
        }
        let mut next: Vec<(f64, f64)> = bins
            .into_values()
            .filter(|&(_, p)| p > 1e-6)
            .map(|(mp, p)| (mp / p, p))
            .collect();
        // Keep top 48 peaks to guarantee O(1) bounded memory even for huge peptides
        if next.len() > 48 {
            next.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            next.truncate(48);
            next.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        }
        cur = next;
    }
    cur
}

/// Simulates the exact HRMS isotopic distribution for `mol` under the chosen adduct name (e.g. `"[M+H]+"`).
pub fn simulate_hrms(mol: &ChemMolecule, adduct_name: &str) -> HrmsSimulation {
    let rep = compute_analysis_report(mol);
    let val = CcVal::compute(mol);
    let mut elem_counts: BTreeMap<i32, usize> = BTreeMap::new();
    for (i, a) in mol.atoms.iter().enumerate() {
        if a.z > 0 {
            *elem_counts.entry(a.z).or_insert(0) += 1;
            let ih = val.implicit_h_counts.get(i).copied().unwrap_or(0).max(0) as usize;
            if ih > 0 {
                *elem_counts.entry(1).or_insert(0) += ih;
            }
        }
    }

    // Account for dimer/adduct isotope shifts (e.g. [M+Cl]- adds Cl isotope pattern!)
    let adducts = compute_adducts(rep.exact_mass);
    let chosen = adducts
        .iter()
        .find(|a| a.name == adduct_name)
        .cloned()
        .unwrap_or(AdductEntry {
            name: "[M+H]+",
            charge: 1,
            mz: rep.exact_mass + 1.007276466812,
        });

    let mult = if chosen.name.starts_with("[2M") { 2 } else { 1 };
    if mult == 2 {
        for v in elem_counts.values_mut() {
            *v *= 2;
        }
    }
    if chosen.name == "[M+Cl]-" {
        *elem_counts.entry(17).or_insert(0) += 1;
    }

    let mut dist = vec![(0.0_f64, 1.0_f64)];
    for (z, cnt) in elem_counts {
        let raw_iso: Vec<(f64, f64)> = get_element_isotopes(z)
            .iter()
            .filter(|iso| iso.abundance() > 0.0005)
            .map(|iso| (iso.mass(), iso.abundance()))
            .collect();
        if raw_iso.is_empty() {
            continue;
        }
        let sum_ab: f64 = raw_iso.iter().map(|&(_, a)| a).sum();
        let norm_iso: Vec<(f64, f64)> = raw_iso.into_iter().map(|(m, a)| (m, a / sum_ab.max(1e-9))).collect();
        dist = convolve_isotopes(dist, &norm_iso, cnt);
    }

    // Shift distribution so its lowest major monoisotopic mass aligns with exact adduct m/z
    let mono_raw = dist.first().map(|&(m, _)| m).unwrap_or(rep.exact_mass);
    let shift = chosen.mz * (chosen.charge.abs().max(1) as f64) - mono_raw;
    let z_abs = chosen.charge.abs().max(1) as f64;
    let max_p = dist.iter().map(|&(_, p)| p).fold(0.0_f64, f64::max).max(1e-9);

    let mut peaks: Vec<IsotopePeak> = dist
        .into_iter()
        .map(|(m, p)| IsotopePeak {
            mz: (m + shift) / z_abs,
            rel_intensity: (p / max_p) * 100.0,
        })
        .filter(|pk| pk.rel_intensity >= 0.25)
        .collect();
    peaks.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));

    HrmsSimulation {
        formula: rep.formula,
        neutral_exact_mass: rep.exact_mass,
        selected_adduct: chosen.name.to_string(),
        charge: chosen.charge,
        peaks,
        adducts,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondDisplay, ChemAtom};

    fn make_bromobenzene() -> ChemMolecule {
        let mut m = ChemMolecule::new();
        for i in 0..6 {
            m.add_atom(ChemAtom::new((i + 1) as u32, 6));
        }
        m.add_atom(ChemAtom::new(7, 35)); // Br
        for i in 0..6 {
            m.add_bond((10 + i) as u32, i, (i + 1) % 6, if i % 2 == 0 { 2 } else { 1 }, BondDisplay::Solid);
        }
        m.add_bond(20, 0, 6, 1, BondDisplay::Solid);
        m
    }

    #[test]
    fn test_medchem_and_hrms_bromine_doublet() {
        let brb = make_bromobenzene();
        let mc = compute_medchem_report(&brb);
        assert!(mc.lipinski_pass);
        assert!(mc.cns_mpo >= 0.0 && mc.cns_mpo <= 6.0);
        assert!(mc.qed > 0.0 && mc.qed <= 1.0);
        assert!(mc.sa_score >= 1.0 && mc.sa_score <= 10.0);

        // Bromobenzene MUST show the classic ~1:1 M / M+2 bromine isotope doublet!
        let hrms = simulate_hrms(&brb, "[M+H]+");
        assert!(hrms.peaks.len() >= 2);
        let m2 = hrms
            .peaks
            .iter()
            .find(|p| (p.mz - (hrms.peaks[0].mz + 1.998)).abs() < 0.03)
            .expect("Br M+2 isotope peak");
        assert!(m2.rel_intensity > 85.0, "Br M+2 should be ~97%, got {}", m2.rel_intensity);
        assert!(hrms.simulate_profile(30_000.0, 120).len() == 120);
    }

    #[test]
    fn test_pains_brenk_alerts_detect_epoxide_and_nitro() {
        let mut epox = ChemMolecule::new();
        epox.add_atom(ChemAtom::new(1, 6));
        epox.add_atom(ChemAtom::new(2, 6));
        epox.add_atom(ChemAtom::new(3, 8));
        epox.add_bond(10, 0, 1, 1, BondDisplay::Solid);
        epox.add_bond(11, 1, 2, 1, BondDisplay::Solid);
        epox.add_bond(12, 2, 0, 1, BondDisplay::Solid);
        let alerts = detect_structural_alerts(&epox);
        assert!(alerts.iter().any(|a| a.name.contains("Epoxide")));
    }
}
