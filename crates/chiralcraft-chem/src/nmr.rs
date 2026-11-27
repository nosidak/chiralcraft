//! Pure-Rust ¹H-NMR and ¹³C-NMR Chemical Shift & Multiplet Prediction Engine
//! for ChiralCraft parity (`Structure > Predict 1H-NMR / 13C-NMR Shifts`).
//!
//! Computes chemically realistic δ (ppm) chemical shifts, vicinal ³J splitting multiplicities,
//! per-atom shift assignments, simulated Lorentzian spectrum curves, and publication-ready
//! NMR multiplet reports directly from `ChemMolecule` + `CcVal` + `CcAromaticityInfo`.

use std::collections::BTreeMap;

use crate::{CcAromaticityInfo, CcVal, ChemMolecule};

#[derive(Debug, Clone, PartialEq)]
pub struct NmrPeak {
    pub atom_idx: usize,
    pub shift_ppm: f64,
    pub protons: i32,
    pub multiplicity: String,
    pub j_hz: Vec<f64>,
    pub environment: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NmrPrediction {
    pub nucleus: &'static str,
    pub freq_mhz: f64,
    pub peaks: Vec<NmrPeak>,
    pub report_text: String,
}

impl NmrPrediction {
    /// Evaluates a Lorentzian sum spectrum over `n_points` spanning `[min_ppm, max_ppm]`.
    /// Returns `(ppm, intensity)` normalized to `0.0..=1.0`.
    pub fn simulate_spectrum(&self, min_ppm: f64, max_ppm: f64, n_points: usize) -> Vec<(f64, f64)> {
        let n = n_points.max(32);
        let gamma = if self.nucleus == "1H" { 0.035 } else { 0.45 };
        let mut pts = Vec::with_capacity(n);
        let mut max_i = 1e-6_f64;

        for k in 0..n {
            let ppm = min_ppm + (max_ppm - min_ppm) * (k as f64) / ((n - 1) as f64);
            let mut y = 0.0;
            for pk in &self.peaks {
                let weight = (pk.protons.max(1) as f64).min(9.0);
                // Render multiplet sub-lines for 1H-NMR
                let sublines: Vec<(f64, f64)> = if self.nucleus == "1H" && !pk.j_hz.is_empty() {
                    let dj = pk.j_hz[0] / self.freq_mhz;
                    match pk.multiplicity.as_str() {
                        "d" => vec![(-0.5 * dj, 0.5), (0.5 * dj, 0.5)],
                        "t" => vec![(-dj, 0.25), (0.0, 0.5), (dj, 0.25)],
                        "q" => vec![
                            (-1.5 * dj, 0.125),
                            (-0.5 * dj, 0.375),
                            (0.5 * dj, 0.375),
                            (1.5 * dj, 0.125),
                        ],
                        "dd" => {
                            let dj2 = pk.j_hz.get(1).copied().unwrap_or(2.2) / self.freq_mhz;
                            vec![
                                (-0.5 * dj - 0.5 * dj2, 0.25),
                                (-0.5 * dj + 0.5 * dj2, 0.25),
                                (0.5 * dj - 0.5 * dj2, 0.25),
                                (0.5 * dj + 0.5 * dj2, 0.25),
                            ]
                        }
                        "m" => vec![
                            (-1.2 * dj, 0.15),
                            (-0.6 * dj, 0.25),
                            (0.0, 0.30),
                            (0.6 * dj, 0.25),
                            (1.2 * dj, 0.15),
                        ],
                        _ => vec![(0.0, 1.0)],
                    }
                } else {
                    vec![(0.0, 1.0)]
                };
                for (off, rel) in sublines {
                    let dx = ppm - (pk.shift_ppm + off);
                    y += weight * rel * (gamma * gamma) / (dx * dx + gamma * gamma);
                }
            }
            if y > max_i {
                max_i = y;
            }
            pts.push((ppm, y));
        }
        for p in &mut pts {
            p.1 /= max_i;
        }
        pts
    }
}

fn bond_order_between(mol: &ChemMolecule, u: usize, v: usize) -> i32 {
    mol.bonds
        .iter()
        .find(|b| (b.begin_atom == u && b.end_atom == v) || (b.begin_atom == v && b.end_atom == u))
        .map(|b| b.order)
        .unwrap_or(1)
}

fn has_carbonyl_o(mol: &ChemMolecule, c_idx: usize) -> bool {
    mol.neighbors(c_idx)
        .iter()
        .any(|&(nb, bi)| mol.atoms[nb].z == 8 && mol.bonds[bi].order == 2)
}

fn is_carboxyl_carbon(mol: &ChemMolecule, c_idx: usize) -> bool {
    if mol.atoms[c_idx].z != 6 {
        return false;
    }
    let mut dbl_o = false;
    let mut sgl_o = false;
    for &(nb, bi) in &mol.neighbors(c_idx) {
        if mol.atoms[nb].z == 8 {
            if mol.bonds[bi].order == 2 {
                dbl_o = true;
            } else if mol.bonds[bi].order == 1 {
                sgl_o = true;
            }
        }
    }
    dbl_o && sgl_o
}

fn is_amide_carbon(mol: &ChemMolecule, c_idx: usize) -> bool {
    mol.atoms[c_idx].z == 6
        && has_carbonyl_o(mol, c_idx)
        && mol.neighbors(c_idx).iter().any(|&(nb, _)| mol.atoms[nb].z == 7)
}

fn is_nitro_n(mol: &ChemMolecule, n_idx: usize) -> bool {
    if mol.atoms[n_idx].z != 7 {
        return false;
    }
    let o_cnt = mol
        .neighbors(n_idx)
        .iter()
        .filter(|&&(nb, _)| mol.atoms[nb].z == 8)
        .count();
    o_cnt >= 2
}

/// Computes shortest graph distance (in bonds) between two atoms up to `max_d`.
fn graph_distance(mol: &ChemMolecule, start: usize, goal: usize, max_d: usize) -> Option<usize> {
    if start == goal {
        return Some(0);
    }
    let mut visited = vec![false; mol.atoms.len()];
    let mut q = std::collections::VecDeque::new();
    visited[start] = true;
    q.push_back((start, 0usize));
    while let Some((u, d)) = q.pop_front() {
        if d >= max_d {
            continue;
        }
        for &(v, _) in &mol.neighbors(u) {
            if v == goal {
                return Some(d + 1);
            }
            if !visited[v] {
                visited[v] = true;
                q.push_back((v, d + 1));
            }
        }
    }
    None
}

/// Predicts ¹H-NMR chemical shifts (400 MHz) and splitting patterns for every H-bearing atom.
pub fn predict_1h_nmr(mol: &ChemMolecule) -> NmrPrediction {
    let val = CcVal::compute(mol);
    let arom = CcAromaticityInfo::perceive(mol);
    let mut peaks = Vec::new();

    // Precompute aromatic substituent increments (ortho, meta, para) relative to benzene (7.27 ppm)
    let mut arom_substituents: Vec<(usize, (f64, f64, f64))> = Vec::new();
    for (i, a) in mol.atoms.iter().enumerate() {
        if arom.is_atom_aromatic(i) && a.z == 6 {
            for &(nb, _) in &mol.neighbors(i) {
                if !arom.is_atom_aromatic(nb) {
                    let inc = match mol.atoms[nb].z {
                        8 => {
                            if mol.neighbors(nb).iter().any(|&(nn, _)| nn != i && has_carbonyl_o(mol, nn)) {
                                (-0.25, 0.03, -0.13) // -OAc
                            } else {
                                (-0.56, -0.12, -0.45) // -OH / -OR
                            }
                        }
                        7 => {
                            if is_nitro_n(mol, nb) {
                                (0.95, 0.26, 0.38) // -NO2
                            } else if mol.neighbors(nb).iter().any(|&(nn, _)| nn != i && has_carbonyl_o(mol, nn)) {
                                (0.12, -0.07, -0.28) // -NHAc
                            } else {
                                (-0.75, -0.25, -0.65) // -NH2 / -NR2
                            }
                        }
                        6 => {
                            if is_carboxyl_carbon(mol, nb) {
                                (0.85, 0.18, 0.28) // -COOH / -COOR
                            } else if has_carbonyl_o(mol, nb) {
                                (0.62, 0.14, 0.21) // -CHO / -COR
                            } else if mol.neighbors(nb).iter().any(|&(nn, bi)| mol.atoms[nn].z == 7 && mol.bonds[bi].order == 3) {
                                (0.36, 0.18, 0.28) // -CN
                            } else {
                                (-0.18, -0.09, -0.18) // -Alkyl
                            }
                        }
                        9 => (-0.29, -0.02, -0.23),
                        17 => (0.03, -0.02, -0.09),
                        35 => (0.18, -0.08, -0.04),
                        53 => (0.39, -0.21, 0.00),
                        16 => (-0.08, -0.16, -0.22),
                        _ => (0.0, 0.0, 0.0),
                    };
                    arom_substituents.push((i, inc));
                }
            }
        }
    }

    for (i, a) in mol.atoms.iter().enumerate() {
        let impl_h = val.implicit_h_counts.get(i).copied().unwrap_or(0) as i32;
        let expl_h = mol
            .neighbors(i)
            .iter()
            .filter(|&&(nb, _)| mol.atoms[nb].z == 1)
            .count() as i32;
        let total_h = impl_h + expl_h;
        if total_h <= 0 || a.z == 1 {
            continue;
        }

        let nbrs = mol.neighbors(i);
        let (shift, env, is_exchangeable) = if a.z == 8 {
            // -OH: carboxylic acid vs phenol vs alcohol
            let attached_c = nbrs.iter().find(|&&(nb, _)| mol.atoms[nb].z == 6).map(|&(nb, _)| nb);
            if let Some(ci) = attached_c {
                if is_carboxyl_carbon(mol, ci) {
                    (11.80, "Carboxylic acid -COOH".to_string(), true)
                } else if arom.is_atom_aromatic(ci) {
                    (5.35, "Phenolic -OH".to_string(), true)
                } else {
                    (3.65, "Alcohol -OH".to_string(), true)
                }
            } else {
                (4.79, "Hydroxyl -OH".to_string(), true)
            }
        } else if a.z == 7 {
            let attached_c = nbrs.iter().find(|&&(nb, _)| mol.atoms[nb].z == 6).map(|&(nb, _)| nb);
            if let Some(ci) = attached_c {
                if has_carbonyl_o(mol, ci) {
                    (8.05, "Amide -NH-".to_string(), true)
                } else if arom.is_atom_aromatic(ci) {
                    (4.00, "Aniline -NH2".to_string(), true)
                } else {
                    (1.85, "Aliphatic amine -NH".to_string(), true)
                }
            } else {
                (2.20, "Amine -NH".to_string(), true)
            }
        } else if a.z == 16 {
            (1.65, "Thiol -SH".to_string(), true)
        } else if a.z == 6 {
            // Carbon-bound protons
            let max_order = nbrs.iter().map(|&(_, bi)| mol.bonds[bi].order).max().unwrap_or(1);
            if has_carbonyl_o(mol, i) {
                let mut s = 9.75;
                if nbrs.iter().any(|&(nb, _)| arom.is_atom_aromatic(nb)) {
                    s = 10.02; // Benzaldehyde
                }
                (s, "Aldehyde -CHO".to_string(), false)
            } else if arom.is_atom_aromatic(i) {
                let mut s = 7.27;
                // Heteroaromatic ring boost (pyridine-like imine N only, excluding pyrrole/indole -NH-)
                for (j, aj) in mol.atoms.iter().enumerate() {
                    if arom.is_atom_aromatic(j) && aj.z == 7 {
                        let j_nbrs = mol.neighbors(j);
                        let j_heavy_deg = j_nbrs.iter().filter(|&&(nb, _)| mol.atoms[nb].z > 1).count();
                        let j_h = val.implicit_h(j) as usize + j_nbrs.iter().filter(|&&(nb, _)| mol.atoms[nb].z == 1).count();
                        if j_heavy_deg == 2 && j_h == 0 {
                            match graph_distance(mol, i, j, 3) {
                                Some(1) => s += 1.32,
                                Some(2) => s += 0.02,
                                Some(3) => s += 0.42,
                                _ => {}
                            }
                        }
                    }
                }
                for &(sub_c, (do_, dm, dp)) in &arom_substituents {
                    if let Some(d) = graph_distance(mol, i, sub_c, 3) {
                        match d {
                            1 => s += do_,
                            2 => s += dm,
                            3 => s += dp,
                            _ => {}
                        }
                    }
                }
                (s.clamp(5.8, 9.4), "Aromatic CH".to_string(), false)
            } else if max_order == 3 {
                (2.55, "Terminal alkyne ≡CH".to_string(), false)
            } else if max_order == 2 {
                // Alkene proton (Pascual–Meier–Simon baseline 5.25 ppm)
                let mut s: f64 = 5.25;
                for &(nb, bi) in &nbrs {
                    if mol.bonds[bi].order == 1 {
                        s += match mol.atoms[nb].z {
                            8 => 1.22,
                            7 => 0.80,
                            17 => 1.08,
                            35 => 1.07,
                            6 if arom.is_atom_aromatic(nb) => 1.38,
                            6 if has_carbonyl_o(mol, nb) => 1.00,
                            6 => 0.45,
                            _ => 0.0,
                        };
                    }
                }
                (s.clamp(4.4, 8.2), "Alkene =CH".to_string(), false)
            } else {
                // sp3 CH4 / CH3 / CH2 / CH (Shoolery's additive increments!)
                let mut s: f64 = match total_h {
                    4 => 0.23,
                    3 => 0.88,
                    2 => 1.25,
                    _ => 1.55,
                };
                for &(nb, _) in &nbrs {
                    let z = mol.atoms[nb].z;
                    if z == 1 {
                        continue;
                    }
                    s += match z {
                        8 => {
                            // Ester/acyloxy O vs ether/alcohol O
                            if mol.neighbors(nb).iter().any(|&(nn, _)| nn != i && has_carbonyl_o(mol, nn)) {
                                2.95 // R-COO-CH2- (~4.15)
                            } else if mol.neighbors(nb).iter().any(|&(nn, _)| nn != i && arom.is_atom_aromatic(nn)) {
                                2.88 // Ar-O-CH3 (~3.76)
                            } else {
                                2.38 // HO-CH2- / R-O-CH2- (~3.63)
                            }
                        }
                        7 => {
                            if is_nitro_n(mol, nb) {
                                3.40
                            } else if mol.neighbors(nb).iter().any(|&(nn, _)| nn != i && has_carbonyl_o(mol, nn)) {
                                2.05
                            } else {
                                1.50
                            }
                        }
                        9 => 3.45,
                        17 => 2.53,
                        35 => 2.33,
                        53 => 2.15,
                        16 => 1.45,
                        6 => {
                            if is_carboxyl_carbon(mol, nb) {
                                1.22 // CH3-COOH (~2.10)
                            } else if has_carbonyl_o(mol, nb) {
                                1.24 // CH3-C=O (~2.12)
                            } else if arom.is_atom_aromatic(nb) {
                                1.44 // Ph-CH3 (~2.32)
                            } else if mol.neighbors(nb).iter().any(|&(nn, bi)| mol.atoms[nn].z == 7 && mol.bonds[bi].order == 3) {
                                1.12 // CH3-CN (~2.00)
                            } else if mol.neighbors(nb).iter().any(|&(_, bi)| mol.bonds[bi].order == 2) {
                                0.82 // Allylic (~1.70)
                            } else if mol.neighbors(nb).iter().any(|&(_, bi)| mol.bonds[bi].order == 3) {
                                0.95 // Propargylic
                            } else {
                                // Beta-inductive effect from heteroatoms on neighbor carbon
                                let beta_het = mol.neighbors(nb).iter().any(|&(nn, _)| nn != i && matches!(mol.atoms[nn].z, 7 | 8 | 9 | 17 | 35));
                                if beta_het { 0.30 } else { 0.05 }
                            }
                        }
                        _ => 0.0,
                    };
                }
                let kind = match total_h {
                    4 => "Methane CH4",
                    3 => "Aliphatic -CH3",
                    2 => "Aliphatic -CH2-",
                    _ => "Aliphatic >CH-",
                };
                (s.clamp(0.15, 7.8), kind.to_string(), false)
            }
        } else {
            (2.0, "Heteroatom H".to_string(), true)
        };

        // Determine vicinal 3J splitting from neighboring carbon-bound protons
        let (mult, js) = if is_exchangeable {
            ("s".to_string(), Vec::new())
        } else {
            let mut vicinal_groups: Vec<i32> = Vec::new();
            for &(nb, _) in &nbrs {
                if mol.atoms[nb].z == 6 {
                    let nb_h = val.implicit_h_counts.get(nb).copied().unwrap_or(0) as i32
                        + mol
                            .neighbors(nb)
                            .iter()
                            .filter(|&&(nn, _)| mol.atoms[nn].z == 1)
                            .count() as i32;
                    if nb_h > 0 {
                        vicinal_groups.push(nb_h);
                    }
                }
            }
            let total_vic: i32 = vicinal_groups.iter().sum();
            let j_typ = if arom.is_atom_aromatic(i) { 7.8 } else { 7.1 };
            match (vicinal_groups.len(), total_vic) {
                (0, _) | (_, 0) => ("s".to_string(), Vec::new()),
                (1, 1) => ("d".to_string(), vec![j_typ]),
                (1, 2) => ("t".to_string(), vec![j_typ]),
                (1, 3) => ("q".to_string(), vec![j_typ]),
                (2, 2) if arom.is_atom_aromatic(i) => ("t".to_string(), vec![7.8]),
                (_, 2) => ("t".to_string(), vec![j_typ]),
                _ => ("m".to_string(), vec![j_typ]),
            }
        };

        peaks.push(NmrPeak {
            atom_idx: i,
            shift_ppm: (shift * 100.0).round() / 100.0,
            protons: total_h,
            multiplicity: mult,
            j_hz: js,
            environment: env,
        });
    }

    // Detect chemically equivalent vicinal pairs (identical predicted shift & environment, e.g. benzene, ethane, ethylene) -> singlet
    let shift_by_atom: std::collections::HashMap<usize, (f64, String)> = peaks
        .iter()
        .map(|p| (p.atom_idx, (p.shift_ppm, p.environment.clone())))
        .collect();
    for p in &mut peaks {
        if p.multiplicity == "s" {
            continue;
        }
        let my_key = shift_by_atom.get(&p.atom_idx);
        let mut non_equiv_vic: Vec<i32> = Vec::new();
        for &(nb, _) in &mol.neighbors(p.atom_idx) {
            if mol.atoms[nb].z == 6 {
                let nb_key = shift_by_atom.get(&nb);
                if nb_key == my_key {
                    continue;
                }
                let nb_h = val.implicit_h_counts.get(nb).copied().unwrap_or(0) as i32
                    + mol.neighbors(nb).iter().filter(|&&(nn, _)| mol.atoms[nn].z == 1).count() as i32;
                if nb_h > 0 {
                    non_equiv_vic.push(nb_h);
                }
            }
        }
        let total_vic: i32 = non_equiv_vic.iter().sum();
        let j_typ = if arom.is_atom_aromatic(p.atom_idx) { 7.8 } else { 7.1 };
        let (m, js) = match (non_equiv_vic.len(), total_vic) {
            (0, _) | (_, 0) => ("s".to_string(), Vec::new()),
            (1, 1) => ("d".to_string(), vec![j_typ]),
            (1, 2) => ("t".to_string(), vec![j_typ]),
            (1, 3) => ("q".to_string(), vec![j_typ]),
            (2, 2) if arom.is_atom_aromatic(p.atom_idx) => ("t".to_string(), vec![7.8]),
            (_, 2) => ("t".to_string(), vec![j_typ]),
            _ => ("m".to_string(), vec![j_typ]),
        };
        p.multiplicity = m;
        p.j_hz = js;
    }

    peaks.sort_by(|a, b| b.shift_ppm.partial_cmp(&a.shift_ppm).unwrap_or(std::cmp::Ordering::Equal));

    // Group equivalent peaks within 0.03 ppm ONLY when multiplicity & environment match (never merge e.g. CH2 q + OH s)
    let mut grouped: BTreeMap<(i32, String, String), (f64, i32, String, Vec<f64>)> = BTreeMap::new();
    for p in &peaks {
        let bucket = (p.shift_ppm * 30.0).round() as i32;
        grouped
            .entry((bucket, p.multiplicity.clone(), p.environment.clone()))
            .and_modify(|e| e.1 += p.protons)
            .or_insert_with(|| (p.shift_ppm, p.protons, p.multiplicity.clone(), p.j_hz.clone()));
    }
    let parts: Vec<String> = grouped
        .values()
        .rev()
        .map(|(s, h, m, js)| {
            if js.is_empty() {
                format!("{s:.2} ({m}, {h}H)")
            } else {
                let jstr: Vec<String> = js.iter().map(|j| format!("{j:.1}")).collect();
                format!("{s:.2} ({m}, J = {} Hz, {h}H)", jstr.join(", "))
            }
        })
        .collect();
    let report_text = if parts.is_empty() {
        "1H NMR (400 MHz): no protons".to_string()
    } else {
        format!("1H NMR (400 MHz): δ {}.", parts.join(", "))
    };

    NmrPrediction {
        nucleus: "1H",
        freq_mhz: 400.0,
        peaks,
        report_text,
    }
}

/// Predicts ¹³C-NMR chemical shifts (100 MHz) for every carbon atom in `mol`.
pub fn predict_13c_nmr(mol: &ChemMolecule) -> NmrPrediction {
    let val = CcVal::compute(mol);
    let arom = CcAromaticityInfo::perceive(mol);
    let mut peaks = Vec::new();

    // Aromatic carbon substituent increments (ipso, ortho, meta, para) from benzene (128.5 ppm)
    let mut arom_subs: Vec<(usize, (f64, f64, f64, f64))> = Vec::new();
    for (i, a) in mol.atoms.iter().enumerate() {
        if arom.is_atom_aromatic(i) && a.z == 6 {
            for &(nb, _) in &mol.neighbors(i) {
                if !arom.is_atom_aromatic(nb) {
                    let inc = match mol.atoms[nb].z {
                        8 => {
                            if mol.neighbors(nb).iter().any(|&(nn, _)| nn != i && has_carbonyl_o(mol, nn)) {
                                (22.4, -7.1, 1.1, -3.2) // -OAc
                            } else {
                                (26.9, -12.7, 1.4, -7.3) // -OH / -OR
                            }
                        }
                        7 => {
                            if is_nitro_n(mol, nb) {
                                (19.9, -4.9, 0.9, 6.1) // -NO2
                            } else {
                                (18.2, -13.4, 0.8, -10.0) // -NH2
                            }
                        }
                        6 => {
                            if is_carboxyl_carbon(mol, nb) {
                                (2.1, 1.5, 0.0, 5.1) // -COOH / -COOR
                            } else if has_carbonyl_o(mol, nb) {
                                (8.6, 1.3, 0.6, 5.5) // -CHO / -COR
                            } else if mol.neighbors(nb).iter().any(|&(nn, bi)| mol.atoms[nn].z == 7 && mol.bonds[bi].order == 3) {
                                (-15.4, 3.6, 0.6, 3.9) // -CN
                            } else {
                                (9.3, 0.8, 0.0, -2.9) // -Alkyl
                            }
                        }
                        9 => (34.8, -12.9, 1.4, -4.5),
                        17 => (6.2, 0.4, 1.3, -1.9),
                        35 => (-5.5, 3.4, 1.7, -1.6),
                        53 => (-34.1, 8.7, 1.4, -1.4),
                        16 => (10.2, -1.5, 0.4, -3.6),
                        _ => (0.0, 0.0, 0.0, 0.0),
                    };
                    arom_subs.push((i, inc));
                }
            }
        }
    }

    for (i, a) in mol.atoms.iter().enumerate() {
        if a.z != 6 {
            continue;
        }
        let nbrs = mol.neighbors(i);
        let max_order = nbrs.iter().map(|&(_, bi)| mol.bonds[bi].order).max().unwrap_or(1);
        let impl_h = val.implicit_h_counts.get(i).copied().unwrap_or(0) as i32;
        let expl_h = nbrs.iter().filter(|&&(nb, _)| mol.atoms[nb].z == 1).count() as i32;
        let total_c_h = impl_h + expl_h;

        let (shift, env) = if has_carbonyl_o(mol, i) {
            if is_carboxyl_carbon(mol, i) {
                let is_ester = nbrs.iter().any(|&(nb, bi)| {
                    mol.atoms[nb].z == 8
                        && mol.bonds[bi].order == 1
                        && mol.neighbors(nb).iter().any(|&(nn, _)| nn != i && mol.atoms[nn].z == 6)
                });
                if is_ester {
                    (170.5, "Ester C=O".to_string())
                } else {
                    (176.2, "Carboxylic acid C=O".to_string())
                }
            } else if is_amide_carbon(mol, i) {
                (170.8, "Amide C=O".to_string())
            } else if total_c_h >= 1 {
                let s = if nbrs.iter().any(|&(nb, _)| arom.is_atom_aromatic(nb)) { 192.3 } else { 200.5 };
                (s, "Aldehyde C=O".to_string())
            } else {
                let s = if nbrs.iter().any(|&(nb, _)| arom.is_atom_aromatic(nb)) { 198.1 } else { 207.2 };
                (s, "Ketone C=O".to_string())
            }
        } else if nbrs.iter().any(|&(nb, bi)| mol.atoms[nb].z == 7 && mol.bonds[bi].order == 3) {
            (118.2, "Nitrile C≡N".to_string())
        } else if arom.is_atom_aromatic(i) {
            let mut s = 128.5;
            for (j, aj) in mol.atoms.iter().enumerate() {
                if arom.is_atom_aromatic(j) && aj.z == 7 {
                    let j_nbrs = mol.neighbors(j);
                    let j_heavy_deg = j_nbrs.iter().filter(|&&(nb, _)| mol.atoms[nb].z > 1).count();
                    let j_h = val.implicit_h(j) as usize + j_nbrs.iter().filter(|&&(nb, _)| mol.atoms[nb].z == 1).count();
                    if j_heavy_deg == 2 && j_h == 0 {
                        match graph_distance(mol, i, j, 3) {
                            Some(1) => s += 21.4,
                            Some(2) => s -= 4.7,
                            Some(3) => s += 7.4,
                            _ => {}
                        }
                    }
                }
            }
            for &(sub_c, (di, do_, dm, dp)) in &arom_subs {
                match graph_distance(mol, i, sub_c, 3) {
                    Some(0) => s += di,
                    Some(1) => s += do_,
                    Some(2) => s += dm,
                    Some(3) => s += dp,
                    _ => {}
                }
            }
            (s.clamp(95.0, 172.0), "Aromatic C".to_string())
        } else if max_order == 3 {
            let s = if total_c_h >= 1 { 68.5 } else { 81.2 };
            (s, "Alkyne C≡".to_string())
        } else if max_order == 2 {
            let mut s: f64 = 123.3;
            for &(nb, _) in &nbrs {
                if bond_order_between(mol, i, nb) == 1 {
                    s += match mol.atoms[nb].z {
                        8 => 28.0,
                        7 => 18.0,
                        17 => 3.0,
                        6 => 10.6,
                        _ => 0.0,
                    };
                }
            }
            (s.clamp(98.0, 162.0), "Alkene C=".to_string())
        } else {
            // Grant–Paul sp3 carbon increments (base -2.3 ppm)
            let mut s: f64 = -2.3;
            let mut alpha_c = 0usize;
            for &(nb, _) in &nbrs {
                match mol.atoms[nb].z {
                    1 => {}
                    6 => {
                        alpha_c += 1;
                        if has_carbonyl_o(mol, nb) {
                            s += 22.5;
                        } else if arom.is_atom_aromatic(nb) {
                            s += 23.0;
                        } else {
                            s += 9.1;
                            // Beta atoms only for aliphatic alpha-carbons (carbonyl/aromatic already include beta deshielding)
                            for &(bb, _) in &mol.neighbors(nb) {
                                if bb != i {
                                    match mol.atoms[bb].z {
                                        6 => s += 9.4,
                                        8 => s += 10.0,
                                        7 => s += 11.0,
                                        17 => s += 10.5,
                                        _ => {}
                                    }
                                }
                            }
                        }
                    }
                    8 => {
                        let is_ester_o = mol.neighbors(nb).iter().any(|&(nn, _)| nn != i && has_carbonyl_o(mol, nn));
                        s += if is_ester_o { 51.0 } else { 49.0 };
                    }
                    7 => s += 29.0,
                    9 => s += 68.0,
                    17 => s += 31.0,
                    35 => s += 20.0,
                    53 => s -= 6.0,
                    16 => s += 11.0,
                    _ => {}
                }
            }
            if alpha_c >= 3 {
                s -= 2.5; // steric branching correction
            }
            (s.clamp(-3.0, 105.0), "Aliphatic C".to_string())
        };

        peaks.push(NmrPeak {
            atom_idx: i,
            shift_ppm: (shift * 10.0).round() / 10.0,
            protons: 1,
            multiplicity: "s".to_string(),
            j_hz: Vec::new(),
            environment: env,
        });
    }

    peaks.sort_by(|a, b| b.shift_ppm.partial_cmp(&a.shift_ppm).unwrap_or(std::cmp::Ordering::Equal));

    let mut uniq_shifts: Vec<f64> = Vec::new();
    for p in &peaks {
        if uniq_shifts.last().map(|&l| (l - p.shift_ppm).abs() > 0.15).unwrap_or(true) {
            uniq_shifts.push(p.shift_ppm);
        }
    }
    let parts: Vec<String> = uniq_shifts.iter().map(|s| format!("{s:.1}")).collect();
    let report_text = if parts.is_empty() {
        "13C NMR (100 MHz): no carbons".to_string()
    } else {
        format!("13C NMR (100 MHz): δ {}.", parts.join(", "))
    };

    NmrPrediction {
        nucleus: "13C",
        freq_mhz: 100.0,
        peaks,
        report_text,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondDisplay, ChemAtom, ChemBond};

    fn make_ethanol() -> ChemMolecule {
        let mut m = ChemMolecule::default();
        m.atoms.push(ChemAtom::new(1, 6)); // CH3
        m.atoms.push(ChemAtom::new(2, 6)); // CH2
        m.atoms.push(ChemAtom::new(3, 8)); // OH
        m.bonds.push(ChemBond { id: 10, begin_atom: 0, end_atom: 1, order: 1, display: BondDisplay::Solid, visible: true, ignore_warnings: false });
        m.bonds.push(ChemBond { id: 11, begin_atom: 1, end_atom: 2, order: 1, display: BondDisplay::Solid, visible: true, ignore_warnings: false });
        m
    }

    fn make_benzene() -> ChemMolecule {
        let mut m = ChemMolecule::default();
        for i in 0..6 {
            m.atoms.push(ChemAtom::new((i + 1) as u32, 6));
        }
        for i in 0..6 {
            m.bonds.push(ChemBond {
                id: (10 + i) as u32,
                begin_atom: i,
                end_atom: (i + 1) % 6,
                order: if i % 2 == 0 { 2 } else { 1 },
                display: BondDisplay::Solid,
                visible: true,
                ignore_warnings: false,
            });
        }
        m
    }

    #[test]
    fn test_nmr_ethanol_and_benzene_realistic_shifts() {
        let eth = make_ethanol();
        let h = predict_1h_nmr(&eth);
        let ch3 = h.peaks.iter().find(|p| p.atom_idx == 0).unwrap();
        let ch2 = h.peaks.iter().find(|p| p.atom_idx == 1).unwrap();
        assert!((ch3.shift_ppm - 1.18).abs() < 0.25 && ch3.multiplicity == "t");
        assert!((ch2.shift_ppm - 3.68).abs() < 0.30 && ch2.multiplicity == "q");

        let c = predict_13c_nmr(&eth);
        let c_ch2 = c.peaks.iter().find(|p| p.atom_idx == 1).unwrap();
        let c_ch3 = c.peaks.iter().find(|p| p.atom_idx == 0).unwrap();
        assert!(c_ch2.shift_ppm > 50.0 && c_ch2.shift_ppm < 66.0);
        assert!(c_ch3.shift_ppm > 12.0 && c_ch3.shift_ppm < 24.0);

        let benz = make_benzene();
        let hb = predict_1h_nmr(&benz);
        assert!((hb.peaks[0].shift_ppm - 7.27).abs() < 0.05);
        let cb = predict_13c_nmr(&benz);
        assert!((cb.peaks[0].shift_ppm - 128.5).abs() < 0.2);
        assert!(!hb.simulate_spectrum(0.0, 12.0, 128).is_empty());
    }
}
