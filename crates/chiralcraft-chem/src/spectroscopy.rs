//! Multi-Spectroscopy Prediction Engine:
//! - Simulated FT-IR (`4000–400 cm⁻¹`, `%T` transmittance curve + clickable bond/atom band assignments)
//! - UV–Vis Chromophore & λmax Estimator (Woodward–Fieser + conjugated π / aromatic / push-pull rules + perceived color)
//! - 2D NMR Correlation Maps (`COSY ¹H–¹H`, `HSQC ¹H–¹³C`, `HMBC 2–3 bond ¹H–¹³C`)

use std::collections::BTreeSet;

use crate::{
    predict_13c_nmr, predict_1h_nmr, CcAromaticityInfo, CcVal, ChemMolecule,
};

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

// ============================================================================
// 1. Simulated FT-IR Spectrum (4000 - 400 cm^-1)
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub struct IrBand {
    pub wavenumber_cm1: f64,
    pub intensity: f64, // 0.0..=1.0 absorption depth
    pub width_cm1: f64,
    pub strength: &'static str, // "s", "m", "w", "br"
    pub assignment: &'static str,
    pub atom_indices: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IrPrediction {
    pub bands: Vec<IrBand>,
    pub report_text: String,
}

impl IrPrediction {
    /// Evaluates a smooth FT-IR %Transmittance curve `(cm^-1, %T 0..100)` from `4000` down to `400 cm^-1`.
    pub fn simulate_transmittance(&self, samples: usize) -> Vec<(f64, f64)> {
        if samples < 2 {
            return Vec::new();
        }
        let mut out = Vec::with_capacity(samples);
        for k in 0..samples {
            let wn = 4000.0 - 3600.0 * (k as f64) / ((samples - 1) as f64);
            let mut abs = 0.0_f64;
            for b in &self.bands {
                let dx = (wn - b.wavenumber_cm1) / b.width_cm1.max(8.0);
                // Lorentzian band shape
                abs += b.intensity / (1.0 + dx * dx);
            }
            let transmittance = (100.0 * (-abs * 2.1).exp()).clamp(4.0, 100.0);
            out.push((wn, transmittance));
        }
        out
    }
}

/// Predicts diagnostic FT-IR vibrational bands and exact atom indices for `mol`.
pub fn predict_ir_spectrum(mol: &ChemMolecule) -> IrPrediction {
    let n = mol.atoms.len();
    let adj = build_adj(mol);
    let val = CcVal::compute(mol);
    let arom = CcAromaticityInfo::perceive(mol);
    let mut bands = Vec::new();

    let mut has_sp3_ch = false;
    let mut has_arom_ch = false;
    let mut has_arom_ring = false;
    let mut sp3_c_atoms = Vec::new();
    let mut arom_atoms = Vec::new();

    for i in 0..n {
        let z = mol.atoms[i].z;
        let imp_h = val.implicit_h_counts.get(i).copied().unwrap_or(0);
        let exp_h = adj[i].iter().filter(|&&(k, _)| mol.atoms[k].z == 1).count() as u8;
        let tot_h = imp_h.saturating_add(exp_h);

        if z == 8 && tot_h >= 1 {
            // Check if carboxylic acid -C(=O)-OH vs phenol vs alcohol
            let is_acid = adj[i].iter().any(|&(c, _)| {
                mol.atoms[c].z == 6 && adj[c].iter().any(|&(o2, o)| o2 != i && mol.atoms[o2].z == 8 && o == 2)
            });
            let is_phenol = adj[i].iter().any(|&(c, _)| arom.is_atom_aromatic(c));
            if is_acid {
                bands.push(IrBand {
                    wavenumber_cm1: 3000.0,
                    intensity: 0.78,
                    width_cm1: 210.0,
                    strength: "br",
                    assignment: "Carboxylic Acid O–H stretch (very broad)",
                    atom_indices: vec![i],
                });
            } else {
                bands.push(IrBand {
                    wavenumber_cm1: if is_phenol { 3380.0 } else { 3340.0 },
                    intensity: 0.72,
                    width_cm1: 110.0,
                    strength: "br",
                    assignment: if is_phenol { "Phenol O–H stretch" } else { "Alcohol O–H stretch (H-bonded)" },
                    atom_indices: vec![i],
                });
            }
        } else if z == 7 && tot_h >= 1 {
            let is_amide = adj[i].iter().any(|&(c, _)| {
                mol.atoms[c].z == 6 && adj[c].iter().any(|&(o2, o)| mol.atoms[o2].z == 8 && o == 2)
            });
            bands.push(IrBand {
                wavenumber_cm1: if is_amide { 3310.0 } else { 3370.0 },
                intensity: 0.52,
                width_cm1: 42.0,
                strength: "m",
                assignment: if is_amide { "Amide N–H stretch" } else { "Amine N–H stretch" },
                atom_indices: vec![i],
            });
        } else if z == 6 {
            if arom.is_atom_aromatic(i) {
                has_arom_ring = true;
                arom_atoms.push(i);
                if tot_h >= 1 {
                    has_arom_ch = true;
                }
            } else if !adj[i].iter().any(|&(_, o)| o >= 2) && tot_h >= 1 {
                has_sp3_ch = true;
                if sp3_c_atoms.len() < 4 {
                    sp3_c_atoms.push(i);
                }
            }
            // Carbonyl C=O diagnostic frequencies!
            if let Some(&(o_idx, _)) = adj[i].iter().find(|&&(k, o)| mol.atoms[k].z == 8 && o == 2) {
                let has_hal = adj[i].iter().any(|&(k, _)| matches!(mol.atoms[k].z, 9 | 17 | 35));
                let has_n = adj[i].iter().any(|&(k, o)| mol.atoms[k].z == 7 && o == 1);
                let has_oh = adj[i].iter().any(|&(k, o)| {
                    k != o_idx && mol.atoms[k].z == 8 && o == 1 && val.implicit_h_counts.get(k).copied().unwrap_or(0) >= 1
                });
                let has_or = adj[i].iter().any(|&(k, o)| {
                    k != o_idx && mol.atoms[k].z == 8 && o == 1 && adj[k].len() >= 2
                });
                let conj = adj[i].iter().any(|&(k, _)| k != o_idx && (arom.is_atom_aromatic(k) || adj[k].iter().any(|&(_, o2)| o2 == 2)));
                let shift = if conj { -20.0 } else { 0.0 };
                let (wn, desc) = if has_hal {
                    (1800.0 + shift, "Acyl Halide C=O stretch")
                } else if has_or {
                    (1742.0 + shift, "Ester C=O stretch")
                } else if has_oh {
                    (1712.0 + shift, "Carboxylic Acid C=O stretch")
                } else if has_n {
                    (1665.0 + shift, "Amide I C=O stretch")
                } else if tot_h >= 1 {
                    (1728.0 + shift, "Aldehyde C=O stretch")
                } else {
                    (1715.0 + shift, "Ketone C=O stretch")
                };
                bands.push(IrBand {
                    wavenumber_cm1: wn,
                    intensity: 0.94,
                    width_cm1: 22.0,
                    strength: "s",
                    assignment: desc,
                    atom_indices: vec![i, o_idx],
                });
                if has_n {
                    bands.push(IrBand {
                        wavenumber_cm1: 1545.0,
                        intensity: 0.58,
                        width_cm1: 24.0,
                        strength: "m",
                        assignment: "Amide II N–H bend / C–N stretch",
                        atom_indices: vec![i, o_idx],
                    });
                }
            }
        } else if z == 7 {
            let oxys: Vec<usize> = adj[i].iter().filter(|&&(k, _)| mol.atoms[k].z == 8 && adj[k].len() == 1).map(|&(k, _)| k).collect();
            if oxys.len() >= 2 {
                let mut at = vec![i];
                at.extend(oxys);
                bands.push(IrBand { wavenumber_cm1: 1530.0, intensity: 0.85, width_cm1: 24.0, strength: "s", assignment: "Nitro N=O asym stretch", atom_indices: at.clone() });
                bands.push(IrBand { wavenumber_cm1: 1350.0, intensity: 0.76, width_cm1: 22.0, strength: "s", assignment: "Nitro N=O sym stretch", atom_indices: at });
            }
        } else if z == 16 {
            let oxys: Vec<usize> = adj[i].iter().filter(|&&(k, o)| mol.atoms[k].z == 8 && o == 2).map(|&(k, _)| k).collect();
            if oxys.len() >= 2 {
                bands.push(IrBand { wavenumber_cm1: 1325.0, intensity: 0.82, width_cm1: 22.0, strength: "s", assignment: "Sulfone / Sulfonamide S=O asym", atom_indices: vec![i, oxys[0], oxys[1]] });
                bands.push(IrBand { wavenumber_cm1: 1155.0, intensity: 0.86, width_cm1: 22.0, strength: "s", assignment: "Sulfone / Sulfonamide S=O sym", atom_indices: vec![i, oxys[0], oxys[1]] });
            } else if oxys.len() == 1 {
                bands.push(IrBand { wavenumber_cm1: 1050.0, intensity: 0.78, width_cm1: 22.0, strength: "s", assignment: "Sulfoxide S=O stretch", atom_indices: vec![i, oxys[0]] });
            }
        }
    }

    // Bond-based bands (C#N, C#C, C=C, C-O, C-F, C-Cl, C-Br)
    for b in &mol.bonds {
        if b.begin_atom >= n || b.end_atom >= n { continue; }
        let (u, v) = (b.begin_atom, b.end_atom);
        let (zu, zv) = (mol.atoms[u].z, mol.atoms[v].z);
        if b.order == 3 && ((zu == 6 && zv == 7) || (zu == 7 && zv == 6)) {
            bands.push(IrBand { wavenumber_cm1: 2235.0, intensity: 0.65, width_cm1: 16.0, strength: "m", assignment: "Nitrile C≡N stretch", atom_indices: vec![u, v] });
        } else if b.order == 3 && zu == 6 && zv == 6 {
            bands.push(IrBand { wavenumber_cm1: 2140.0, intensity: 0.45, width_cm1: 18.0, strength: "w", assignment: "Alkyne C≡C stretch", atom_indices: vec![u, v] });
        } else if b.order == 2 && zu == 6 && zv == 6 && !arom.is_atom_aromatic(u) && !arom.is_atom_aromatic(v) {
            bands.push(IrBand { wavenumber_cm1: 1648.0, intensity: 0.52, width_cm1: 20.0, strength: "m", assignment: "Alkene C=C stretch", atom_indices: vec![u, v] });
        } else if b.order == 1 && ((zu == 6 && zv == 8) || (zu == 8 && zv == 6)) {
            bands.push(IrBand { wavenumber_cm1: 1185.0, intensity: 0.64, width_cm1: 28.0, strength: "s", assignment: "C–O stretch (ether/ester/alcohol)", atom_indices: vec![u, v] });
        } else if b.order == 1 && ((zu == 6 && zv == 9) || (zu == 9 && zv == 6)) {
            bands.push(IrBand { wavenumber_cm1: 1140.0, intensity: 0.74, width_cm1: 26.0, strength: "s", assignment: "C–F stretch", atom_indices: vec![u, v] });
        } else if b.order == 1 && ((zu == 6 && zv == 17) || (zu == 17 && zv == 6)) {
            bands.push(IrBand { wavenumber_cm1: 745.0, intensity: 0.60, width_cm1: 22.0, strength: "m", assignment: "C–Cl stretch", atom_indices: vec![u, v] });
        } else if b.order == 1 && ((zu == 6 && zv == 35) || (zu == 35 && zv == 6)) {
            bands.push(IrBand { wavenumber_cm1: 620.0, intensity: 0.55, width_cm1: 22.0, strength: "m", assignment: "C–Br stretch", atom_indices: vec![u, v] });
        }
    }

    if has_arom_ch {
        bands.push(IrBand { wavenumber_cm1: 3055.0, intensity: 0.36, width_cm1: 24.0, strength: "w", assignment: "Aromatic C–H stretch", atom_indices: arom_atoms.iter().take(4).copied().collect() });
    }
    if has_sp3_ch {
        bands.push(IrBand { wavenumber_cm1: 2935.0, intensity: 0.55, width_cm1: 35.0, strength: "m", assignment: "Aliphatic sp³ C–H stretch", atom_indices: sp3_c_atoms });
    }
    if has_arom_ring {
        bands.push(IrBand { wavenumber_cm1: 1595.0, intensity: 0.52, width_cm1: 20.0, strength: "m", assignment: "Aromatic ring C=C stretch", atom_indices: arom_atoms.clone() });
        bands.push(IrBand { wavenumber_cm1: 1490.0, intensity: 0.50, width_cm1: 20.0, strength: "m", assignment: "Aromatic ring breathing", atom_indices: arom_atoms.clone() });
        bands.push(IrBand { wavenumber_cm1: 755.0, intensity: 0.62, width_cm1: 24.0, strength: "s", assignment: "Aromatic C–H out-of-plane bend", atom_indices: arom_atoms });
    }

    // Deduplicate closely spaced identical assignments and sort high -> low cm^-1
    bands.sort_by(|a, b| b.wavenumber_cm1.partial_cmp(&a.wavenumber_cm1).unwrap_or(std::cmp::Ordering::Equal));
    let mut dedup: Vec<IrBand> = Vec::new();
    for b in bands {
        if !dedup.iter().any(|d| d.assignment == b.assignment && (d.wavenumber_cm1 - b.wavenumber_cm1).abs() < 12.0) {
            dedup.push(b);
        }
    }
    let summary: Vec<String> = dedup
        .iter()
        .take(8)
        .map(|b| format!("{:.0} ({})", b.wavenumber_cm1, b.strength))
        .collect();
    let report_text = format!("IR (neat, cm⁻¹): νmax {}.", summary.join(", "));
    IrPrediction { bands: dedup, report_text }
}

// ============================================================================
// 2. UV-Vis Chromophore & λmax Estimator
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub struct UvVisPrediction {
    pub lambda_max_nm: f64,
    pub epsilon_est: f64,
    pub transition: &'static str,
    pub perceived_color_name: &'static str,
    pub perceived_rgb: (u8, u8, u8),
    pub chromophore: String,
}

/// Estimates primary UV–Vis absorption λmax (nm), molar absorptivity ε, and perceived solution color.
pub fn predict_uv_vis(mol: &ChemMolecule) -> UvVisPrediction {
    let n = mol.atoms.len();
    let adj = build_adj(mol);
    let arom = CcAromaticityInfo::perceive(mol);

    // Find largest connected conjugated π-system (aromatic + sp2/sp atoms)
    let is_pi: Vec<bool> = (0..n)
        .map(|i| arom.is_atom_aromatic(i) || adj[i].iter().any(|&(_, o)| o >= 2))
        .collect();
    let mut vis = vec![false; n];
    let mut max_pi_comp = 0usize;
    let mut max_arom_comp = 0usize;
    let mut has_azo = false;
    let mut has_nitro = false;
    let mut has_donor = false;
    let mut has_carbonyl = false;

    for b in &mol.bonds {
        if b.begin_atom < n && b.end_atom < n {
            let (zu, zv) = (mol.atoms[b.begin_atom].z, mol.atoms[b.end_atom].z);
            if b.order == 2 && zu == 7 && zv == 7 {
                has_azo = true;
            }
            if b.order == 2 && ((zu == 6 && zv == 8) || (zu == 8 && zv == 6)) {
                has_carbonyl = true;
            }
        }
    }
    for i in 0..n {
        if mol.atoms[i].z == 7 && adj[i].iter().filter(|&&(k, _)| mol.atoms[k].z == 8).count() >= 2 {
            has_nitro = true;
        }
        if matches!(mol.atoms[i].z, 7 | 8) && adj[i].iter().any(|&(k, o)| o == 1 && is_pi[k]) {
            has_donor = true;
        }
    }

    for i in 0..n {
        if !is_pi[i] || vis[i] {
            continue;
        }
        let mut stack = vec![i];
        vis[i] = true;
        let mut cnt = 0usize;
        let mut acnt = 0usize;
        while let Some(u) = stack.pop() {
            cnt += 1;
            if arom.is_atom_aromatic(u) { acnt += 1; }
            for &(v, _) in &adj[u] {
                if is_pi[v] && !vis[v] {
                    vis[v] = true;
                    stack.push(v);
                }
            }
        }
        if cnt > max_pi_comp {
            max_pi_comp = cnt;
            max_arom_comp = acnt;
        }
    }

    let (mut lmax, mut eps, trans, chromo): (f64, f64, &'static str, String) = if max_pi_comp >= 20 && mol.atoms.iter().filter(|a| a.z == 7).count() >= 4 {
        (418.0, 220_000.0, "π → π* (Soret band)", "Porphyrin / Macrocyclic tetrapyrrole".to_string())
    } else if has_azo && max_arom_comp >= 6 {
        (415.0 + if has_donor { 45.0 } else { 0.0 }, 24_000.0, "π → π* / n → π*", "Azobenzene / Azo chromophore".to_string())
    } else if max_arom_comp >= 14 {
        (378.0, 18_000.0, "π → π*", format!("Extended polyaromatic ({max_arom_comp} ring atoms)"))
    } else if max_arom_comp >= 10 {
        (312.0, 11_000.0, "π → π*", "Naphthalene / Bicyclic aromatic".to_string())
    } else if max_arom_comp >= 5 {
        let base = 254.0 + ((max_pi_comp.saturating_sub(6)) as f64) * 14.0;
        let push_pull = if has_nitro && has_donor { 115.0 } else if has_nitro { 45.0 } else if has_donor && has_carbonyl { 35.0 } else if has_donor { 22.0 } else { 0.0 };
        (base + push_pull, 14_000.0, "π → π*", if has_nitro && has_donor { "Push–pull donor/acceptor aromatic".into() } else { "Substituted aromatic chromophore".into() })
    } else if max_pi_comp >= 4 {
        // Conjugated diene / polyene / enone Woodward–Fieser approximation
        let extra_db = max_pi_comp.saturating_sub(4) / 2;
        (217.0 + (extra_db as f64) * 30.0 + if has_carbonyl { 12.0 } else { 0.0 }, 18_500.0, "π → π*", format!("Conjugated π-system ({max_pi_comp} sp² atoms)"))
    } else if has_carbonyl {
        (278.0, 25.0, "n → π* (weak)", "Isolated carbonyl".to_string())
    } else {
        (190.0, 100.0, "σ → σ* (far UV)", "Saturated / No UV–Vis chromophore".to_string())
    };
    lmax = lmax.clamp(185.0, 750.0);
    if lmax > 400.0 {
        eps = eps.max(15_000.0);
    }

    // Complementary perceived solution color from absorbed λmax!
    let (cname, rgb) = if lmax < 385.0 {
        ("Colorless (UV absorber)", (242, 246, 250))
    } else if lmax < 435.0 {
        ("Pale Yellow / Yellow", (248, 224, 92))
    } else if lmax < 480.0 {
        ("Yellow–Orange", (245, 162, 54))
    } else if lmax < 515.0 {
        ("Orange–Red", (228, 86, 56))
    } else if lmax < 570.0 {
        ("Crimson / Violet-Red", (195, 48, 128))
    } else if lmax < 620.0 {
        ("Indigo / Deep Blue", (56, 98, 212))
    } else {
        ("Teal / Emerald Green", (42, 168, 128))
    };

    UvVisPrediction {
        lambda_max_nm: lmax,
        epsilon_est: eps,
        transition: trans,
        perceived_color_name: cname,
        perceived_rgb: rgb,
        chromophore: chromo,
    }
}

// ============================================================================
// 3. 2D NMR Correlation Maps (COSY, HSQC, HMBC)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Nmr2dKind {
    Cosy,
    Hsqc,
    Hmbc,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Nmr2dCrossPeak {
    pub x_ppm: f64,
    pub y_ppm: f64,
    pub atom_x: usize,
    pub atom_y: usize,
    pub label: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Nmr2dPrediction {
    pub kind: Nmr2dKind,
    pub x_range: (f64, f64),
    pub y_range: (f64, f64),
    pub peaks: Vec<Nmr2dCrossPeak>,
}

/// Predicts 2D NMR cross-peaks for COSY (vicinal ¹H–¹H), HSQC (1-bond ¹H–¹³C), or HMBC (2–3 bond ¹H–¹³C).
pub fn predict_2d_nmr(mol: &ChemMolecule, kind: Nmr2dKind) -> Nmr2dPrediction {
    let n = mol.atoms.len();
    let adj = build_adj(mol);
    let h1 = predict_1h_nmr(mol);
    let c13 = predict_13c_nmr(mol);

    let h_by_atom: std::collections::HashMap<usize, f64> = h1
        .peaks
        .iter()
        .map(|p| (p.atom_idx, p.shift_ppm))
        .collect();
    let c_by_atom: std::collections::HashMap<usize, f64> = c13
        .peaks
        .iter()
        .map(|p| (p.atom_idx, p.shift_ppm))
        .collect();

    let mut peaks = Vec::new();
    match kind {
        Nmr2dKind::Cosy => {
            // Diagonal peaks + vicinal (bonded C-C) 3J H-H off-diagonal cross-peaks
            for (&ai, &sh) in &h_by_atom {
                peaks.push(Nmr2dCrossPeak {
                    x_ppm: sh,
                    y_ppm: sh,
                    atom_x: ai,
                    atom_y: ai,
                    label: format!("Diag #{} ({:.2} ppm)", ai + 1, sh),
                });
            }
            for b in &mol.bonds {
                if b.begin_atom < n && b.end_atom < n {
                    if let (Some(&s1), Some(&s2)) = (h_by_atom.get(&b.begin_atom), h_by_atom.get(&b.end_atom)) {
                        peaks.push(Nmr2dCrossPeak {
                            x_ppm: s1,
                            y_ppm: s2,
                            atom_x: b.begin_atom,
                            atom_y: b.end_atom,
                            label: format!("³J H#{} ↔ H#{} ({:.2}, {:.2})", b.begin_atom + 1, b.end_atom + 1, s1, s2),
                        });
                        peaks.push(Nmr2dCrossPeak {
                            x_ppm: s2,
                            y_ppm: s1,
                            atom_x: b.end_atom,
                            atom_y: b.begin_atom,
                            label: format!("³J H#{} ↔ H#{} ({:.2}, {:.2})", b.end_atom + 1, b.begin_atom + 1, s2, s1),
                        });
                    }
                }
            }
            Nmr2dPrediction { kind, x_range: (0.0, 10.5), y_range: (0.0, 10.5), peaks }
        }
        Nmr2dKind::Hsqc => {
            // Direct 1-bond ¹H–¹³C correlations on protonated carbons
            for (&ai, &sh) in &h_by_atom {
                if let Some(&sc) = c_by_atom.get(&ai) {
                    peaks.push(Nmr2dCrossPeak {
                        x_ppm: sh,
                        y_ppm: sc,
                        atom_x: ai,
                        atom_y: ai,
                        label: format!("¹J C#{}–H ({:.2}, {:.1})", ai + 1, sh, sc),
                    });
                }
            }
            Nmr2dPrediction { kind, x_range: (0.0, 10.5), y_range: (0.0, 210.0), peaks }
        }
        Nmr2dKind::Hmbc => {
            // 2-bond and 3-bond ¹H–¹³C correlations (graph distance 1 or 2 between protonated host atom & target C)
            let mut seen = BTreeSet::new();
            for (&ha, &sh) in &h_by_atom {
                // Distance 1 (2J_CH) and distance 2 (3J_CH)
                for &(n1, _) in &adj[ha] {
                    if let Some(&sc) = c_by_atom.get(&n1) {
                        if seen.insert((ha, n1)) {
                            peaks.push(Nmr2dCrossPeak {
                                x_ppm: sh,
                                y_ppm: sc,
                                atom_x: ha,
                                atom_y: n1,
                                label: format!("²J H#{} → C#{} ({:.2}, {:.1})", ha + 1, n1 + 1, sh, sc),
                            });
                        }
                    }
                    for &(n2, _) in &adj[n1] {
                        if n2 != ha {
                            if let Some(&sc) = c_by_atom.get(&n2) {
                                if seen.insert((ha, n2)) {
                                    peaks.push(Nmr2dCrossPeak {
                                        x_ppm: sh,
                                        y_ppm: sc,
                                        atom_x: ha,
                                        atom_y: n2,
                                        label: format!("³J H#{} → C#{} ({:.2}, {:.1})", ha + 1, n2 + 1, sh, sc),
                                    });
                                }
                            }
                        }
                    }
                }
            }
            Nmr2dPrediction { kind, x_range: (0.0, 10.5), y_range: (0.0, 210.0), peaks }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondDisplay, ChemAtom};

    #[test]
    fn test_ir_uv_and_2d_nmr_on_acetophenone_like() {
        // Ethanol: CH3-CH2-OH
        let mut eth = ChemMolecule::new();
        eth.add_atom(ChemAtom::new(1, 6));
        eth.add_atom(ChemAtom::new(2, 6));
        eth.add_atom(ChemAtom::new(3, 8));
        eth.add_bond(10, 0, 1, 1, BondDisplay::Solid);
        eth.add_bond(11, 1, 2, 1, BondDisplay::Solid);

        let ir = predict_ir_spectrum(&eth);
        assert!(ir.bands.iter().any(|b| b.assignment.contains("O–H")));
        assert_eq!(ir.simulate_transmittance(100).len(), 100);

        let cosy = predict_2d_nmr(&eth, Nmr2dKind::Cosy);
        assert!(cosy.peaks.iter().any(|p| p.atom_x == 0 && p.atom_y == 1));
        let hsqc = predict_2d_nmr(&eth, Nmr2dKind::Hsqc);
        assert!(hsqc.peaks.len() >= 2);
        let hmbc = predict_2d_nmr(&eth, Nmr2dKind::Hmbc);
        assert!(!hmbc.peaks.is_empty());

        let uv = predict_uv_vis(&eth);
        assert!(uv.lambda_max_nm < 220.0);
    }
}
