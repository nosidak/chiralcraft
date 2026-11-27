//! All-Discipline Flagship Engines in `chiralcraft-doc`:
//! 1. Reaction Energy Profile (`ΔG‡`) Data Model, Eyring Kinetics, Boltzmann Selectivity & SVG Renderer
//! 2. ELN Lab Reagent / Solvent Safety Handbook (GHS Hazards, Density, BP, Liquid Volume) & ACS/JOC SI Generator
//! 3. Native LaTeX (`\chemfig` / `\schemestart`) & TikZ Code Exporter
//! 4. Multi-Compound SDF / SMILES / CSV Virtual Screening Library Loader & Filter Exporter
//! 5. Built-In Clean-Room Organometallic Catalysts Constructors

use std::collections::BTreeMap;

use chiralcraft_chem::{
    compute_analysis_report, compute_medchem_report, get_element, predict_13c_nmr, predict_1h_nmr,
    predict_ir_spectrum, simulate_hrms, BondDisplay, CcVal, ChemAtom, ChemMolecule, MedChemReport,
};
use chiralcraft_formats::{
    molecule_to_smiles, molecules_to_sdf, sdf_to_molecules, smiles_to_molecule, SdfRecord,
};

use crate::{clean_up_structure_preserving_stereo, InteractiveDocument, StoichiometryRow};

// ============================================================================
// 1. Reaction Energy Profile (ΔG‡) & Eyring / Boltzmann Kinetics
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub struct EnergyState {
    pub label: String,
    pub energy: f64, // kcal/mol (or kJ/mol if profile.unit_kj)
    pub is_ts: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EnergyProfile {
    pub id: u32,
    pub bbox: (f64, f64, f64, f64),
    pub unit_kj: bool,
    pub temp_k: f64,
    pub pathway_a: Vec<EnergyState>,
    pub pathway_b: Vec<EnergyState>,
}

impl EnergyProfile {
    pub fn starter(id: u32, x: f64, y: f64) -> Self {
        Self {
            id,
            bbox: (x, y, x + 230.0, y + 135.0),
            unit_kj: false,
            temp_k: 298.15,
            pathway_a: vec![
                EnergyState { label: "SM".into(), energy: 0.0, is_ts: false },
                EnergyState { label: "TS1".into(), energy: 18.4, is_ts: true },
                EnergyState { label: "Int".into(), energy: 6.2, is_ts: false },
                EnergyState { label: "TS2".into(), energy: 15.1, is_ts: true },
                EnergyState { label: "Prod".into(), energy: -8.5, is_ts: false },
            ],
            pathway_b: vec![
                EnergyState { label: "SM".into(), energy: 0.0, is_ts: false },
                EnergyState { label: "TS1'".into(), energy: 20.6, is_ts: true },
                EnergyState { label: "Int'".into(), energy: 8.0, is_ts: false },
                EnergyState { label: "TS2'".into(), energy: 16.8, is_ts: true },
                EnergyState { label: "Prod'".into(), energy: -6.2, is_ts: false },
            ],
        }
    }

    pub fn unit_str(&self) -> &'static str {
        if self.unit_kj { "kJ/mol" } else { "kcal/mol" }
    }

    /// Maximum activation barrier ΔG‡ (converted to kcal/mol) along `pathway_a`.
    pub fn max_barrier_kcal(&self) -> f64 {
        let start = self.pathway_a.first().map(|s| s.energy).unwrap_or(0.0);
        let peak = self
            .pathway_a
            .iter()
            .map(|s| s.energy)
            .fold(start, f64::max);
        let raw = (peak - start).max(0.0);
        if self.unit_kj { raw / 4.184 } else { raw }
    }

    /// ΔΔG‡ (in kcal/mol) between Pathway B peak and Pathway A peak (if Pathway B is non-empty).
    pub fn ddg_kcal(&self) -> Option<f64> {
        if self.pathway_b.is_empty() || self.pathway_a.is_empty() {
            return None;
        }
        let peak_a = self.pathway_a.iter().map(|s| s.energy).fold(f64::NEG_INFINITY, f64::max);
        let peak_b = self.pathway_b.iter().map(|s| s.energy).fold(f64::NEG_INFINITY, f64::max);
        let diff = (peak_b - peak_a).abs();
        Some(if self.unit_kj { diff / 4.184 } else { diff })
    }
}

/// Eyring equation: returns `(rate_constant_s_inv, half_life_seconds)` for barrier `dg_kcal` at `temp_k`.
pub fn eyring_rate_and_halflife(dg_kcal: f64, temp_k: f64) -> (f64, f64) {
    const KB_OVER_H: f64 = 2.083_661_912e10; // s^-1 K^-1
    const R_KCAL: f64 = 1.987_204_258e-3; // kcal mol^-1 K^-1
    let t = temp_k.max(10.0);
    let k = KB_OVER_H * t * (-dg_kcal / (R_KCAL * t)).exp();
    let t_half = if k > 1e-30 { std::f64::consts::LN_2 / k } else { f64::INFINITY };
    (k, t_half)
}

/// Boltzmann selectivity from `ddg_kcal` at `temp_k`: returns `(major_pct, ee_pct)`.
pub fn boltzmann_selectivity(ddg_kcal: f64, temp_k: f64) -> (f64, f64) {
    const R_KCAL: f64 = 1.987_204_258e-3;
    let t = temp_k.max(10.0);
    // Numerically stable formulation: inv = exp(-|ΔΔG‡| / RT) ∈ [0, 1], never overflows to inf/inf = NaN!
    let inv = (-ddg_kcal.abs() / (R_KCAL * t)).exp();
    let major_pct = (1.0 / (1.0 + inv)) * 100.0;
    let ee_pct = ((1.0 - inv) / (1.0 + inv)) * 100.0;
    (major_pct, ee_pct)
}

pub fn format_halflife(sec: f64) -> String {
    if !sec.is_finite() || sec > 3.15e12 {
        "> 100,000 yr".to_string()
    } else if sec > 3.15e7 {
        format!("{:.1} yr", sec / 3.15e7)
    } else if sec > 86400.0 {
        format!("{:.1} d", sec / 86400.0)
    } else if sec > 3600.0 {
        format!("{:.1} h", sec / 3600.0)
    } else if sec > 60.0 {
        format!("{:.1} min", sec / 60.0)
    } else if sec >= 1e-3 {
        format!("{:.3} s", sec)
    } else {
        format!("{:.2e} s", sec)
    }
}

/// Renders all `EnergyProfile` diagrams in a document scene to clean SVG paths & labels.
pub fn render_energy_profiles_svg(profiles: &[EnergyProfile], ink: &str) -> String {
    let mut svg = String::new();
    for ep in profiles {
        let (x1, y1, x2, y2) = ep.bbox;
        let rx = x1.min(x2);
        let ry = y1.min(y2);
        let w = (x2 - x1).abs().max(120.0);
        let h = (y2 - y1).abs().max(80.0);
        let mut min_e = f64::INFINITY;
        let mut max_e = f64::NEG_INFINITY;
        for st in ep.pathway_a.iter().chain(ep.pathway_b.iter()) {
            min_e = min_e.min(st.energy);
            max_e = max_e.max(st.energy);
        }
        if !min_e.is_finite() {
            continue;
        }
        let pad_e = ((max_e - min_e) * 0.18).max(3.0);
        min_e -= pad_e;
        max_e += pad_e;
        let e_span = (max_e - min_e).max(1.0);
        let plot_l = rx + 26.0;
        let plot_r = rx + w - 10.0;
        let plot_t = ry + 14.0;
        let plot_b = ry + h - 18.0;

        // Axes
        svg.push_str(&format!(
            "<line x1=\"{plot_l:.2}\" y1=\"{plot_t:.2}\" x2=\"{plot_l:.2}\" y2=\"{plot_b:.2}\" stroke=\"{ink}\" stroke-width=\"0.9\"/>\n"
        ));
        svg.push_str(&format!(
            "<line x1=\"{plot_l:.2}\" y1=\"{plot_b:.2}\" x2=\"{plot_r:.2}\" y2=\"{plot_b:.2}\" stroke=\"{ink}\" stroke-width=\"0.9\"/>\n"
        ));
        svg.push_str(&format!(
            "<text x=\"{:.2}\" y=\"{:.2}\" font-size=\"6.5\" fill=\"{ink}\">ΔG ({})</text>\n",
            rx + 2.0, plot_t - 3.0, ep.unit_str()
        ));

        let draw_path = |svg: &mut String, states: &[EnergyState], other: &[EnergyState], color: &str, dashed: bool, is_secondary: bool| {
            let n = states.len();
            if n == 0 {
                return;
            }
            let step = (plot_r - plot_l) / (n as f64);
            let bar_hw = (step * 0.24).clamp(5.0, 14.0);
            let mut pts = Vec::with_capacity(n);
            for (i, st) in states.iter().enumerate() {
                let cx = plot_l + (i as f64 + 0.5) * step;
                let cy = plot_b - ((st.energy - min_e) / e_span) * (plot_b - plot_t);
                pts.push((cx, cy, st));
            }
            let dash_attr = if dashed { " stroke-dasharray=\"3,2\"" } else { "" };
            for w_pair in pts.windows(2) {
                let (x_a, y_a, _) = w_pair[0];
                let (x_b, y_b, _) = w_pair[1];
                let sx = x_a + bar_hw;
                let ex = x_b - bar_hw;
                let mx = 0.5 * (sx + ex);
                svg.push_str(&format!(
                    "<path d=\"M {sx:.2} {y_a:.2} C {mx:.2} {y_a:.2}, {mx:.2} {y_b:.2}, {ex:.2} {y_b:.2}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"0.9\"{dash_attr}/>\n"
                ));
            }
            for (i, &(cx, cy, st)) in pts.iter().enumerate() {
                svg.push_str(&format!(
                    "<line x1=\"{:.2}\" y1=\"{cy:.2}\" x2=\"{:.2}\" y2=\"{cy:.2}\" stroke=\"{color}\" stroke-width=\"2.1\"/>\n",
                    cx - bar_hw, cx + bar_hw
                ));
                let mut ly = cy - 4.0;
                if let Some(ost) = other.get(i) {
                    if (st.energy - ost.energy).abs() < 0.05 && st.label == ost.label {
                        if is_secondary {
                            continue;
                        }
                    } else if st.energy < ost.energy {
                        ly = cy + 8.5;
                    }
                }
                let esc_lbl = st.label.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
                svg.push_str(&format!(
                    "<text x=\"{cx:.2}\" y=\"{ly:.2}\" font-size=\"6.2\" text-anchor=\"middle\" fill=\"{color}\">{esc_lbl} ({:.1})</text>\n",
                    st.energy
                ));
            }
        };

        draw_path(&mut svg, &ep.pathway_a, &ep.pathway_b, "#1657B8", false, false);
        if !ep.pathway_b.is_empty() {
            draw_path(&mut svg, &ep.pathway_b, &ep.pathway_a, "#C83232", true, true);
        }
    }
    svg
}

// ============================================================================
// 2. ELN Lab Solvent & Reagent Safety Handbook + ACS/JOC SI Procedure Generator
// ============================================================================

#[derive(Debug, Clone, Copy)]
pub struct LabReagentInfo {
    pub name: &'static str,
    pub abbrev: &'static str,
    pub formula: &'static str,
    pub mw: f64,
    pub density_g_ml: Option<f64>, // None for solids
    pub bp_c: Option<f64>,
    pub ghs: &'static str,
}

pub const LAB_REAGENTS: &[LabReagentInfo] = &[
    LabReagentInfo { name: "Tetrahydrofuran", abbrev: "THF", formula: "C4H8O", mw: 72.11, density_g_ml: Some(0.889), bp_c: Some(66.0), ghs: "Flammable | Irritant | Peroxide-Former | CMR(Carc2)" },
    LabReagentInfo { name: "Dichloromethane", abbrev: "DCM", formula: "CH2Cl2", mw: 84.93, density_g_ml: Some(1.326), bp_c: Some(40.0), ghs: "CMR(Carc2) | Irritant | Toxic" },
    LabReagentInfo { name: "N,N-Dimethylformamide", abbrev: "DMF", formula: "C3H7NO", mw: 73.09, density_g_ml: Some(0.944), bp_c: Some(153.0), ghs: "Flammable | CMR(Repr1B) | Acutely Toxic" },
    LabReagentInfo { name: "Dimethyl Sulfoxide", abbrev: "DMSO", formula: "C2H6OS", mw: 78.13, density_g_ml: Some(1.100), bp_c: Some(189.0), ghs: "Skin-Penetrating Carrier | Combustible" },
    LabReagentInfo { name: "Acetonitrile", abbrev: "MeCN", formula: "C2H3N", mw: 41.05, density_g_ml: Some(0.786), bp_c: Some(82.0), ghs: "Flammable | Acutely Toxic | Irritant" },
    LabReagentInfo { name: "Ethyl Acetate", abbrev: "EtOAc", formula: "C4H8O2", mw: 88.11, density_g_ml: Some(0.902), bp_c: Some(77.0), ghs: "Flammable | Eye Irritant" },
    LabReagentInfo { name: "Toluene", abbrev: "PhMe", formula: "C7H8", mw: 92.14, density_g_ml: Some(0.867), bp_c: Some(111.0), ghs: "Flammable | CMR(Repr2) | Aspiration Hazard" },
    LabReagentInfo { name: "Diethyl Ether", abbrev: "Et2O", formula: "C4H10O", mw: 74.12, density_g_ml: Some(0.713), bp_c: Some(35.0), ghs: "Extremely Flammable | Peroxide-Former" },
    LabReagentInfo { name: "1,4-Dioxane", abbrev: "Dioxane", formula: "C4H8O2", mw: 88.11, density_g_ml: Some(1.033), bp_c: Some(101.0), ghs: "Flammable | CMR(Carc1B) | Peroxide-Former" },
    LabReagentInfo { name: "Methanol", abbrev: "MeOH", formula: "CH4O", mw: 32.04, density_g_ml: Some(0.792), bp_c: Some(65.0), ghs: "Flammable | Acutely Toxic (Optic Nerve)" },
    LabReagentInfo { name: "Ethanol", abbrev: "EtOH", formula: "C2H6O", mw: 46.07, density_g_ml: Some(0.789), bp_c: Some(78.0), ghs: "Flammable | Eye Irritant" },
    LabReagentInfo { name: "n-Hexane", abbrev: "Hex", formula: "C6H14", mw: 86.18, density_g_ml: Some(0.659), bp_c: Some(69.0), ghs: "Flammable | Neurotoxin (STOT RE2) | CMR" },
    LabReagentInfo { name: "Chloroform", abbrev: "CHCl3", formula: "CHCl3", mw: 119.38, density_g_ml: Some(1.489), bp_c: Some(61.0), ghs: "Acutely Toxic | CMR(Carc2)" },
    LabReagentInfo { name: "Triethylamine", abbrev: "Et3N", formula: "C6H15N", mw: 101.19, density_g_ml: Some(0.726), bp_c: Some(89.0), ghs: "Flammable | Corrosive | Acutely Toxic" },
    LabReagentInfo { name: "N,N-Diisopropylethylamine", abbrev: "DIPEA", formula: "C8H19N", mw: 129.24, density_g_ml: Some(0.742), bp_c: Some(127.0), ghs: "Flammable | Corrosive | Irritant" },
    LabReagentInfo { name: "Pyridine", abbrev: "Py", formula: "C5H5N", mw: 79.10, density_g_ml: Some(0.982), bp_c: Some(115.0), ghs: "Flammable | Acutely Toxic | Irritant" },
    LabReagentInfo { name: "Trifluoroacetic Acid", abbrev: "TFA", formula: "C2HF3O2", mw: 114.02, density_g_ml: Some(1.489), bp_c: Some(72.0), ghs: "Corrosive | Acutely Toxic (Inhalation)" },
    LabReagentInfo { name: "Thionyl Chloride", abbrev: "SOCl2", formula: "SOCl2", mw: 118.97, density_g_ml: Some(1.631), bp_c: Some(75.0), ghs: "Water-Reactive | Corrosive | Acutely Toxic" },
    LabReagentInfo { name: "Oxalyl Chloride", abbrev: "(COCl)2", formula: "C2Cl2O2", mw: 126.93, density_g_ml: Some(1.478), bp_c: Some(63.0), ghs: "Water-Reactive | Corrosive | Acutely Toxic (CO)" },
    LabReagentInfo { name: "n-Butyllithium (2.5 M Hex)", abbrev: "n-BuLi", formula: "C4H9Li", mw: 64.06, density_g_ml: Some(0.693), bp_c: None, ghs: "PYROPHORIC | Water-Reactive | Corrosive | Flammable" },
    LabReagentInfo { name: "Sodium Hydride (60% disp)", abbrev: "NaH", formula: "NaH", mw: 24.00, density_g_ml: None, bp_c: None, ghs: "Water-Reactive (H2 Gas) | Flammable Solid" },
    LabReagentInfo { name: "Lithium Aluminum Hydride", abbrev: "LiAlH4", formula: "LiAlH4", mw: 37.95, density_g_ml: None, bp_c: None, ghs: "Water-Reactive | PYROPHORIC | Corrosive" },
    LabReagentInfo { name: "Sodium Borohydride", abbrev: "NaBH4", formula: "NaBH4", mw: 37.83, density_g_ml: None, bp_c: None, ghs: "Water-Reactive | Acutely Toxic | Corrosive" },
    LabReagentInfo { name: "Potassium Carbonate", abbrev: "K2CO3", formula: "K2CO3", mw: 138.21, density_g_ml: None, bp_c: None, ghs: "Irritant (Skin/Eye)" },
    LabReagentInfo { name: "HATU", abbrev: "HATU", formula: "C10H15F6N6OP", mw: 380.23, density_g_ml: None, bp_c: None, ghs: "Sensitizer | Explosive Dust Caution" },
    LabReagentInfo { name: "EDC Hydrochloride", abbrev: "EDC·HCl", formula: "C8H18ClN3", mw: 191.70, density_g_ml: None, bp_c: None, ghs: "Irritant | Skin Sensitizer" },
    LabReagentInfo { name: "4-Dimethylaminopyridine", abbrev: "DMAP", formula: "C7H10N2", mw: 122.17, density_g_ml: None, bp_c: None, ghs: "Acutely Toxic (Dermal/Oral)" },
    LabReagentInfo { name: "Tetrakis(triphenylphosphine)Pd(0)", abbrev: "Pd(PPh3)4", formula: "C72H60P4Pd", mw: 1155.56, density_g_ml: None, bp_c: None, ghs: "Air/Light-Sensitive | Irritant" },
    LabReagentInfo { name: "m-CPBA (77%)", abbrev: "mCPBA", formula: "C7H5ClO3", mw: 172.57, density_g_ml: None, bp_c: None, ghs: "Organic Peroxide | Oxidizer | Corrosive" },
    LabReagentInfo { name: "Acetic Anhydride", abbrev: "Ac2O", formula: "C4H6O3", mw: 102.09, density_g_ml: Some(1.082), bp_c: Some(140.0), ghs: "Flammable | Corrosive | Acutely Toxic" },
];

pub fn lookup_lab_reagent(query: &str) -> Vec<&'static LabReagentInfo> {
    let q = query.trim().to_ascii_lowercase();
    LAB_REAGENTS
        .iter()
        .filter(|r| {
            q.is_empty()
                || r.name.to_ascii_lowercase().contains(&q)
                || r.abbrev.to_ascii_lowercase().contains(&q)
                || r.formula.to_ascii_lowercase().contains(&q)
                || r.ghs.to_ascii_lowercase().contains(&q)
        })
        .collect()
}

/// Generates a publication-ready ACS/JOC Supporting Information experimental procedure
/// + predicted ¹H NMR, ¹³C NMR, IR, and HRMS characterization paragraph!
pub fn generate_si_experimental_block(
    product: &ChemMolecule,
    product_name: &str,
    stoich: &[StoichiometryRow],
    solvent: &str,
    solvent_ml: f64,
    temp_str: &str,
    time_str: &str,
    yield_pct: f64,
) -> String {
    let rep = compute_analysis_report(product);
    let h1 = predict_1h_nmr(product);
    let c13 = predict_13c_nmr(product);
    let ir = predict_ir_spectrum(product);
    let hrms = simulate_hrms(product, "[M+H]+");
    let mh_mz = hrms.peaks.first().map(|p| p.mz).unwrap_or(rep.exact_mass + 1.0073);

    let reactants_desc: Vec<String> = stoich
        .iter()
        .filter(|r| r.role != "Product")
        .map(|r| format!("{} ({:.1} mg, {:.2} mmol, {:.2} equiv)", r.label, r.mass_mg, r.mmol, r.equivalents))
        .collect();
    let charge_str = if reactants_desc.is_empty() {
        "the starting substrate (1.00 mmol, 1.00 equiv) and reagents".to_string()
    } else {
        reactants_desc.join(", ")
    };

    let base_mmol = stoich.first().map(|r| r.mmol).unwrap_or(1.0);
    let prod_mmol = base_mmol * (yield_pct / 100.0);
    let prod_mg = prod_mmol * rep.mol_weight;

    format!(
        "General Procedure for {product_name} ({}):\n\
         To an oven-dried Schlenk flask equipped with a magnetic stir bar under N₂ was charged {charge_str} in anhydrous {solvent} ({solvent_ml:.1} mL). \
         The resulting mixture was stirred at {temp_str} for {time_str} until TLC analysis indicated complete consumption of starting material. \
         Upon cooling to room temperature, the reaction mixture was diluted with EtOAc and washed sequentially with saturated aqueous NH₄Cl and brine. \
         The organic layer was dried over anhydrous Na₂SO₄, filtered, and concentrated under reduced pressure. \
         Purification by flash column chromatography on silica gel afforded {product_name} ({prod_mg:.1} mg, {prod_mmol:.2} mmol, {yield_pct:.0}% yield).\n\n\
         Characterization Data:\n\
         {}\n\
         {}\n\
         {}\n\
         HRMS (ESI+): calcd for {} [M+H]⁺ {:.4}; found {:.4}.",
        rep.formula, h1.report_text, c13.report_text, ir.report_text, rep.formula, mh_mz, mh_mz
    )
}

// ============================================================================
// 3. Native LaTeX (`\chemfig` / `\schemestart`) & TikZ Exporter
// ============================================================================

/// Converts a `ChemMolecule` into compilable LaTeX `\chemfig{...}` syntax!
pub fn molecule_to_chemfig(mol: &ChemMolecule) -> String {
    let n = mol.atoms.len();
    if n == 0 {
        return "\\chemfig{}".to_string();
    }
    let val = CcVal::compute(mol);
    let mut adj: Vec<Vec<(usize, i32, BondDisplay)>> = vec![Vec::new(); n];
    for b in &mol.bonds {
        if b.begin_atom < n && b.end_atom < n {
            let rev_disp = match b.display {
                BondDisplay::WedgeBegin => BondDisplay::WedgeEnd,
                BondDisplay::WedgeEnd => BondDisplay::WedgeBegin,
                BondDisplay::WedgedHashBegin => BondDisplay::WedgedHashEnd,
                BondDisplay::WedgedHashEnd => BondDisplay::WedgedHashBegin,
                other => other,
            };
            adj[b.begin_atom].push((b.end_atom, b.order, b.display));
            adj[b.end_atom].push((b.begin_atom, b.order, rev_disp));
        }
    }

    let atom_tex = |idx: usize| -> String {
        let a = &mol.atoms[idx];
        if let Some(lbl) = &a.label {
            return lbl.clone();
        }
        if a.z == 6 && a.charge == 0 && !adj[idx].is_empty() {
            return String::new();
        }
        let sym = get_element(a.z).map(|e| e.symbol).unwrap_or("?");
        let h = val.implicit_h_counts.get(idx).copied().unwrap_or(0);
        let mut s = sym.to_string();
        if h == 1 {
            s.push('H');
        } else if h > 1 {
            s.push_str(&format!("H_{{{h}}}"));
        }
        if a.charge > 0 {
            if a.charge == 1 {
                s.push_str("^{+}");
            } else {
                s.push_str(&format!("^{{{}+}}", a.charge));
            }
        } else if a.charge < 0 {
            if a.charge == -1 {
                s.push_str("^{-}");
            } else {
                s.push_str(&format!("^{{{}-}}", -a.charge));
            }
        }
        s
    };

    let bond_tex = |u: usize, v: usize, order: i32, disp: BondDisplay| -> String {
        let dx = mol.atoms[v].x - mol.atoms[u].x;
        let dy = -(mol.atoms[v].y - mol.atoms[u].y); // LaTeX Y is up!
        let ang = (dy.atan2(dx).to_degrees().rem_euclid(360.0) / 15.0).round() as i32 * 15;
        let bsym = match disp {
            BondDisplay::WedgeBegin | BondDisplay::WedgeEnd => "<",
            BondDisplay::WedgedHashBegin | BondDisplay::WedgedHashEnd | BondDisplay::Hash => "<:",
            _ => match order {
                2 => "=",
                3 => "~",
                _ => "-",
            },
        };
        format!("{bsym}[:{ang}]")
    };

    let mut visited = vec![false; n];
    let mut ring_id = 0usize;
    let mut ring_tags: BTreeMap<usize, Vec<String>> = BTreeMap::new();
    let mut tree_edges: Vec<Vec<(usize, i32, BondDisplay)>> = vec![Vec::new(); n];
    let mut back_edges: Vec<(usize, usize, i32)> = Vec::new();

    // First pass: build DFS spanning tree & discover ring closure hooks
    fn dfs_tree(
        u: usize,
        parent: Option<usize>,
        adj: &[Vec<(usize, i32, BondDisplay)>],
        vis: &mut [bool],
        tree: &mut [Vec<(usize, i32, BondDisplay)>],
        backs: &mut Vec<(usize, usize, i32)>,
    ) {
        vis[u] = true;
        for &(v, ord, disp) in &adj[u] {
            if Some(v) == parent {
                continue;
            }
            if !vis[v] {
                tree[u].push((v, ord, disp));
                dfs_tree(v, Some(u), adj, vis, tree, backs);
            } else if u < v {
                backs.push((u, v, ord));
            }
        }
    }
    let mut comp_roots = Vec::new();
    for start in 0..n {
        if !visited[start] {
            comp_roots.push(start);
            dfs_tree(start, None, &adj, &mut visited, &mut tree_edges, &mut back_edges);
        }
    }
    for (u, v, ord) in back_edges {
        ring_id += 1;
        let bkind = if ord == 2 { ",=" } else { "" };
        ring_tags.entry(u).or_default().push(format!("?[r{ring_id}]"));
        ring_tags.entry(v).or_default().push(format!("?[r{ring_id}{bkind}]"));
    }

    fn emit(
        u: usize,
        tree: &[Vec<(usize, i32, BondDisplay)>],
        ring_tags: &BTreeMap<usize, Vec<String>>,
        atom_tex: &dyn Fn(usize) -> String,
        bond_tex: &dyn Fn(usize, usize, i32, BondDisplay) -> String,
        out: &mut String,
    ) {
        out.push_str(&atom_tex(u));
        if let Some(tags) = ring_tags.get(&u) {
            for t in tags {
                out.push_str(t);
            }
        }
        let children = &tree[u];
        for (idx, &(v, ord, disp)) in children.iter().enumerate() {
            let is_last = idx + 1 == children.len();
            if !is_last {
                out.push('(');
            }
            out.push_str(&bond_tex(u, v, ord, disp));
            emit(v, tree, ring_tags, atom_tex, bond_tex, out);
            if !is_last {
                out.push(')');
            }
        }
    }

    let mut parts = Vec::with_capacity(comp_roots.len());
    for root in comp_roots {
        let mut body = String::new();
        emit(root, &tree_edges, &ring_tags, &atom_tex, &bond_tex, &mut body);
        parts.push(format!("\\chemfig{{{body}}}"));
    }
    parts.join(" \\+ ")
}

/// Exports the full document scene as a compilable LaTeX `\schemestart ... \schemestop` snippet.
pub fn document_to_chemfig_scheme(doc: &InteractiveDocument) -> String {
    let mut out = String::from("% Generated by ChiralCraft — requires \\usepackage{chemfig}\n\\schemestart\n");
    for (i, m) in doc.scene.molecules.iter().enumerate() {
        if i > 0 {
            if i == 1 && !doc.scene.arrows.is_empty() {
                out.push_str("  \\arrow{->[Reagents][Conditions]}\n");
            } else {
                out.push_str("  \\+ \n");
            }
        }
        out.push_str("  ");
        out.push_str(&molecule_to_chemfig(m));
        out.push('\n');
    }
    out.push_str("\\schemestop\n");
    out
}

/// Exports the active document scene as standalone LaTeX `tikzpicture` vector graphics!
pub fn document_to_tikz(doc: &InteractiveDocument) -> String {
    let mut out = String::from("% Generated by ChiralCraft — requires \\usepackage{tikz}\n\\begin{tikzpicture}[x=0.035cm, y=-0.035cm, line cap=round]\n");
    for mol in &doc.scene.molecules {
        let val = CcVal::compute(mol);
        for b in &mol.bonds {
            if b.begin_atom >= mol.atoms.len() || b.end_atom >= mol.atoms.len() { continue; }
            let a1 = &mol.atoms[b.begin_atom];
            let a2 = &mol.atoms[b.end_atom];
            let style = match b.order {
                2 => "double, double distance=1.6pt",
                3 => "line width=1.4pt",
                _ => "line width=0.65pt",
            };
            out.push_str(&format!(
                "  \\draw[{style}] ({:.2},{:.2}) -- ({:.2},{:.2});\n",
                a1.x, a1.y, a2.x, a2.y
            ));
        }
        for (ai, a) in mol.atoms.iter().enumerate() {
            let lbl = a.label.clone().or_else(|| {
                if a.z != 6 || a.charge != 0 {
                    let sym = get_element(a.z).map(|e| e.symbol).unwrap_or("?");
                    let h = val.implicit_h_counts.get(ai).copied().unwrap_or(0);
                    Some(if h > 1 { format!("{sym}H$_{h}$") } else if h == 1 { format!("{sym}H") } else { sym.to_string() })
                } else {
                    None
                }
            });
            if let Some(txt) = lbl {
                out.push_str(&format!(
                    "  \\node[fill=white, inner sep=1pt, font=\\sffamily\\scriptsize] at ({:.2},{:.2}) {{{txt}}};\n",
                    a.x, a.y
                ));
            }
        }
    }
    for arr in &doc.scene.arrows {
        out.push_str(&format!(
            "  \\draw[->, >=stealth, line width=0.8pt] ({:.2},{:.2}) -- ({:.2},{:.2});\n",
            arr.tail.0, arr.tail.1, arr.head.0, arr.head.1
        ));
    }
    out.push_str("\\end{tikzpicture}\n");
    out
}

// ============================================================================
// 4. Multi-Compound SDF / SMILES / CSV Virtual Screening Library Engine
// ============================================================================

#[derive(Debug, Clone)]
pub struct ScreeningRecord {
    pub index: usize,
    pub title: String,
    pub smiles: String,
    pub formula: String,
    pub medchem: MedChemReport,
    pub molecule: ChemMolecule,
}

/// Loads a multi-compound `.sdf`, `.smi`, or `.csv` file and computes live MedChem & PAINS descriptors for every entry!
pub fn load_screening_library(path: &str, bond_len: f64) -> Result<Vec<ScreeningRecord>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let low = path.to_ascii_lowercase();
    let mut raw: Vec<(String, ChemMolecule)> = Vec::new();
    if low.ends_with(".sdf") || low.ends_with(".sd") {
        for (i, rec) in sdf_to_molecules(&text)?.into_iter().enumerate() {
            let title = if rec.title.trim().is_empty() { format!("CMPD-{:03}", i + 1) } else { rec.title };
            raw.push((title, rec.molecule));
        }
    } else {
        for (i, line) in text.lines().enumerate() {
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') || t.to_ascii_lowercase().starts_with("smiles") {
                continue;
            }
            let (smi, title) = if t.contains(',') {
                let mut parts = t.splitn(2, ',');
                let a = parts.next().unwrap_or("").trim().trim_matches('"');
                let b = parts.next().unwrap_or("").trim().trim_matches('"');
                // Support either "SMILES,Name" or "ID,SMILES"
                if smiles_to_molecule(a).is_ok() {
                    (a.to_string(), if b.is_empty() { format!("CMPD-{:03}", i + 1) } else { b.to_string() })
                } else {
                    (b.split(',').next().unwrap_or("").trim().trim_matches('"').to_string(), a.to_string())
                }
            } else {
                let mut parts = t.split_whitespace();
                let a = parts.next().unwrap_or("");
                let b = parts.next().unwrap_or("");
                (a.to_string(), if b.is_empty() { format!("CMPD-{:03}", i + 1) } else { b.to_string() })
            };
            if let Ok(mut m) = smiles_to_molecule(&smi) {
                clean_up_structure_preserving_stereo(&mut m, &[], bond_len);
                raw.push((title, m));
            }
        }
    }
    if raw.is_empty() {
        return Err("No valid chemical structures found in library file".to_string());
    }
    Ok(raw
        .into_iter()
        .enumerate()
        .map(|(idx, (title, molecule))| {
            let rep = compute_analysis_report(&molecule);
            let medchem = compute_medchem_report(&molecule);
            let smiles = molecule_to_smiles(&molecule);
            ScreeningRecord {
                index: idx + 1,
                title,
                smiles,
                formula: rep.formula,
                medchem,
                molecule,
            }
        })
        .collect())
}

/// Built-in starter screening library of 10 classic approved drugs for instant virtual-screening demos.
pub fn builtin_demo_screening_library(bond_len: f64) -> Vec<ScreeningRecord> {
    let items = [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
        ("Acetaminophen", "CC(=O)Nc1ccc(O)cc1"),
        ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
        ("Naproxen", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
        ("Lidocaine", "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
        ("Diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
        ("4-Nitroaniline (Alert)", "Nc1ccc([N+](=O)[O-])cc1"),
        ("Catechol (PAINS)", "Oc1ccccc1O"),
        ("Metformin", "CN(C)C(=N)NC(=N)N"),
    ];
    items
        .iter()
        .enumerate()
        .filter_map(|(i, &(title, smi))| {
            let mut m = smiles_to_molecule(smi).ok()?;
            clean_up_structure_preserving_stereo(&mut m, &[], bond_len);
            let rep = compute_analysis_report(&m);
            let medchem = compute_medchem_report(&m);
            Some(ScreeningRecord {
                index: i + 1,
                title: title.to_string(),
                smiles: smi.to_string(),
                formula: rep.formula,
                medchem,
                molecule: m,
            })
        })
        .collect()
}

pub fn export_screening_sdf(records: &[&ScreeningRecord]) -> String {
    let sdf_recs: Vec<SdfRecord> = records
        .iter()
        .map(|r| {
            let mut props = BTreeMap::new();
            props.insert("Formula".into(), r.formula.clone());
            props.insert("MW".into(), format!("{:.2}", r.medchem.mw));
            props.insert("LogP".into(), format!("{:.2}", r.medchem.log_p));
            props.insert("CNS_MPO".into(), format!("{:.2}", r.medchem.cns_mpo));
            props.insert("QED".into(), format!("{:.2}", r.medchem.qed));
            props.insert("PAINS_Alerts".into(), r.medchem.alerts.len().to_string());
            SdfRecord {
                title: r.title.clone(),
                molecule: r.molecule.clone(),
                properties: props,
            }
        })
        .collect();
    molecules_to_sdf(&sdf_recs)
}

// ============================================================================
// 5. Built-In Clean-Room Organometallic Catalysts Constructors
// ============================================================================

pub fn build_cisplatin(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let mut pt = ChemAtom::new(1, 78); pt.x = 0.0; pt.y = 0.0;
    let mut cl1 = ChemAtom::new(2, 17); cl1.x = -blen; cl1.y = -blen * 0.85;
    let mut cl2 = ChemAtom::new(3, 17); cl2.x = -blen; cl2.y = blen * 0.85;
    let mut n1 = ChemAtom::new(4, 7); n1.x = blen; n1.y = -blen * 0.85; n1.label = Some("NH3".into());
    let mut n2 = ChemAtom::new(5, 7); n2.x = blen; n2.y = blen * 0.85; n2.label = Some("NH3".into());
    for a in [pt, cl1, cl2, n1, n2] { m.add_atom(a); }
    m.add_bond(10, 0, 1, 1, BondDisplay::Solid);
    m.add_bond(11, 0, 2, 1, BondDisplay::Solid);
    m.add_bond(12, 3, 0, 5, BondDisplay::Dative);
    m.add_bond(13, 4, 0, 5, BondDisplay::Dative);
    m
}

pub fn build_pd_pph3_4(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let mut pd = ChemAtom::new(1, 46); pd.x = 0.0; pd.y = 0.0;
    m.add_atom(pd);
    let dirs = [(-1.0, -0.8), (1.0, -0.8), (1.0, 0.8), (-1.0, 0.8)];
    for (k, &(dx, dy)) in dirs.iter().enumerate() {
        let mut p = ChemAtom::new((k + 2) as u32, 15);
        p.x = dx * blen;
        p.y = dy * blen;
        p.label = Some("PPh3".into());
        m.add_atom(p);
        m.add_bond((10 + k) as u32, k + 1, 0, 5, BondDisplay::Dative);
    }
    m
}

pub fn build_wilkinson_catalyst(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let mut rh = ChemAtom::new(1, 45); rh.x = 0.0; rh.y = 0.0;
    let mut cl = ChemAtom::new(2, 17); cl.x = -blen * 1.1; cl.y = 0.0;
    let mut p1 = ChemAtom::new(3, 15); p1.x = 0.0; p1.y = -blen * 1.1; p1.label = Some("PPh3".into());
    let mut p2 = ChemAtom::new(4, 15); p2.x = blen * 1.1; p2.y = 0.0; p2.label = Some("PPh3".into());
    let mut p3 = ChemAtom::new(5, 15); p3.x = 0.0; p3.y = blen * 1.1; p3.label = Some("PPh3".into());
    for a in [rh, cl, p1, p2, p3] { m.add_atom(a); }
    m.add_bond(10, 0, 1, 1, BondDisplay::Solid);
    m.add_bond(11, 2, 0, 5, BondDisplay::Dative);
    m.add_bond(12, 3, 0, 5, BondDisplay::Dative);
    m.add_bond(13, 4, 0, 5, BondDisplay::Dative);
    m
}

pub fn build_vaska_complex(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let mut ir = ChemAtom::new(1, 77); ir.x = 0.0; ir.y = 0.0;
    let mut cl = ChemAtom::new(2, 17); cl.x = -blen * 1.1; cl.y = 0.0;
    let mut co = ChemAtom::new(3, 6); co.x = blen * 1.1; co.y = 0.0; co.label = Some("CO".into());
    let mut p1 = ChemAtom::new(4, 15); p1.x = 0.0; p1.y = -blen * 1.1; p1.label = Some("PPh3".into());
    let mut p2 = ChemAtom::new(5, 15); p2.x = 0.0; p2.y = blen * 1.1; p2.label = Some("PPh3".into());
    for a in [ir, cl, co, p1, p2] { m.add_atom(a); }
    m.add_bond(10, 0, 1, 1, BondDisplay::Solid);
    m.add_bond(11, 2, 0, 5, BondDisplay::Dative);
    m.add_bond(12, 3, 0, 5, BondDisplay::Dative);
    m.add_bond(13, 4, 0, 5, BondDisplay::Dative);
    m
}

pub fn build_ferrocene_sandwich(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let mut fe = ChemAtom::new(1, 26); fe.x = 0.0; fe.y = 0.0;
    let mut cp1 = ChemAtom::new(2, 6); cp1.x = 0.0; cp1.y = -blen * 1.15; cp1.label = Some("Cp".into());
    let mut cp2 = ChemAtom::new(3, 6); cp2.x = 0.0; cp2.y = blen * 1.15; cp2.label = Some("Cp".into());
    for a in [fe, cp1, cp2] { m.add_atom(a); }
    m.add_bond(10, 0, 1, 1, BondDisplay::Solid);
    m.add_bond(11, 0, 2, 1, BondDisplay::Solid);
    m
}

pub fn build_grubbs_core(blen: f64) -> ChemMolecule {
    let mut m = ChemMolecule::new();
    let mut ru = ChemAtom::new(1, 44); ru.x = 0.0; ru.y = 0.0;
    let mut cl1 = ChemAtom::new(2, 17); cl1.x = -blen; cl1.y = -blen * 0.7;
    let mut cl2 = ChemAtom::new(3, 17); cl2.x = -blen; cl2.y = blen * 0.7;
    let mut p1 = ChemAtom::new(4, 15); p1.x = 0.0; p1.y = -blen * 1.1; p1.label = Some("PCy3".into());
    let mut p2 = ChemAtom::new(5, 15); p2.x = 0.0; p2.y = blen * 1.1; p2.label = Some("PCy3".into());
    let mut carb = ChemAtom::new(6, 6); carb.x = blen * 1.1; carb.y = 0.0; carb.label = Some("CHPh".into());
    for a in [ru, cl1, cl2, p1, p2, carb] { m.add_atom(a); }
    m.add_bond(10, 0, 1, 1, BondDisplay::Solid);
    m.add_bond(11, 0, 2, 1, BondDisplay::Solid);
    m.add_bond(12, 3, 0, 5, BondDisplay::Dative);
    m.add_bond(13, 4, 0, 5, BondDisplay::Dative);
    m.add_bond(14, 0, 5, 2, BondDisplay::Solid);
    m
}

// ============================================================================
// 6. Lossless `<objecttag Name="chiralcraft_extensions">` Persistence
//    (Preserves Energy Profiles, Table Cells, and TLC Metadata in CDXML/SVG/PNG/PDF!)
// ============================================================================

fn esc_field(s: &str) -> String {
    s.replace('%', "%25")
        .replace('|', "%7C")
        .replace('\r', "%0D")
        .replace('\n', "%0A")
        .replace(';', "%3B")
        .replace(',', "%2C")
}

fn unesc_field(s: &str) -> String {
    let bytes = s.as_bytes();
    let mut out = String::with_capacity(s.len());
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            match &s[i..i + 3] {
                "%25" => { out.push('%'); i += 3; continue; }
                "%7C" => { out.push('|'); i += 3; continue; }
                "%0D" => { out.push('\r'); i += 3; continue; }
                "%0A" => { out.push('\n'); i += 3; continue; }
                "%3B" => { out.push(';'); i += 3; continue; }
                "%2C" => { out.push(','); i += 3; continue; }
                _ => {}
            }
        }
        if let Some(ch) = s[i..].chars().next() {
            out.push(ch);
            i += ch.len_utf8();
        } else {
            break;
        }
    }
    out
}

/// Serializes `energy_profiles`, `tables` (with cell text), and `tlc_plates` (with labels, colors, solvent)
/// into a compact base64-encoded CDXML `<objecttag Name="chiralcraft_extensions" Value="..."/>` element.
pub fn serialize_chiralcraft_extensions_cdxml(
    energy_profiles: &[EnergyProfile],
    tables: &[crate::SceneTable],
    tlc_plates: &[crate::SceneTlcPlate],
) -> Option<String> {
    if energy_profiles.is_empty() && tables.is_empty() && tlc_plates.is_empty() {
        return None;
    }
    let mut raw = String::new();
    for ep in energy_profiles {
        let fmt_states = |list: &[EnergyState]| -> String {
            list.iter()
                .map(|s| format!("{},{:.4},{}", esc_field(&s.label), s.energy, if s.is_ts { 1 } else { 0 }))
                .collect::<Vec<_>>()
                .join(";")
        };
        raw.push_str(&format!(
            "EP|{}|{:.2},{:.2},{:.2},{:.2}|{}|{:.2}|{}|{}\n",
            ep.id,
            ep.bbox.0,
            ep.bbox.1,
            ep.bbox.2,
            ep.bbox.3,
            if ep.unit_kj { 1 } else { 0 },
            ep.temp_k,
            fmt_states(&ep.pathway_a),
            fmt_states(&ep.pathway_b),
        ));
    }
    for tbl in tables {
        let rows_enc: Vec<String> = tbl
            .cells
            .iter()
            .map(|r| r.iter().map(|c| esc_field(c)).collect::<Vec<_>>().join(","))
            .collect();
        raw.push_str(&format!(
            "TBL|{}|{}|{}|{:.2},{:.2},{:.2},{:.2}|{}\n",
            tbl.id,
            tbl.rows,
            tbl.cols,
            tbl.bbox.0,
            tbl.bbox.1,
            tbl.bbox.2,
            tbl.bbox.3,
            rows_enc.join(";")
        ));
    }
    for tlc in tlc_plates {
        let rfs: Vec<String> = tlc.spot_rfs.iter().map(|v| format!("{v:.4}")).collect();
        let lbls: Vec<String> = tlc.spot_labels.iter().map(|l| esc_field(l)).collect();
        let cols: Vec<String> = tlc.spot_colors.iter().map(|&(r, g, b)| format!("{r}:{g}:{b}")).collect();
        raw.push_str(&format!(
            "TLC|{}|{:.2},{:.2},{:.2},{:.2}|{:.4}|{}|{}|{}|{}\n",
            tlc.id,
            tlc.bbox.0,
            tlc.bbox.1,
            tlc.bbox.2,
            tlc.bbox.3,
            tlc.origin_fraction,
            rfs.join(","),
            lbls.join(","),
            cols.join(","),
            esc_field(&tlc.solvent_label)
        ));
    }
    let b64 = chiralcraft_formats::base64_encode(raw.as_bytes());
    Some(format!(
        "<objecttag id=\"990001\" Name=\"chiralcraft_extensions\" Value=\"{b64}\"/>"
    ))
}

/// Restores extended ChiralCraft objects from `<objecttag Name="chiralcraft_extensions">` inside a parsed CDXML.
pub fn restore_chiralcraft_extensions_from_cdxml(
    cdxml: &chiralcraft_formats::CdxmlDocument,
    energy_profiles: &mut Vec<EnergyProfile>,
    tables: &mut Vec<crate::SceneTable>,
    tlc_plates: &mut Vec<crate::SceneTlcPlate>,
) {
    let Some(b64) = cdxml
        .elements
        .iter()
        .find(|e| e.tag == "objecttag" && e.attr("Name") == Some("chiralcraft_extensions"))
        .and_then(|e| e.attr("Value"))
    else {
        return;
    };
    let Ok(bytes) = chiralcraft_formats::base64_decode(b64) else {
        return;
    };
    let Ok(raw) = String::from_utf8(bytes) else {
        return;
    };

    let mut parsed_tables = Vec::new();
    let mut parsed_tlcs = Vec::new();
    let parse_bbox = |s: &str| -> Option<(f64, f64, f64, f64)> {
        let v: Vec<f64> = s.split(',').filter_map(|t| t.parse().ok()).collect();
        if v.len() == 4 { Some((v[0], v[1], v[2], v[3])) } else { None }
    };
    let parse_states = |s: &str| -> Vec<EnergyState> {
        if s.is_empty() { return Vec::new(); }
        s.split(';')
            .filter_map(|item| {
                let p: Vec<&str> = item.split(',').collect();
                if p.len() >= 3 {
                    Some(EnergyState {
                        label: unesc_field(p[0]),
                        energy: p[1].parse().ok()?,
                        is_ts: p[2] == "1",
                    })
                } else {
                    None
                }
            })
            .collect()
    };

    for line in raw.lines() {
        let parts: Vec<&str> = line.split('|').collect();
        match parts.first().copied() {
            Some("EP") if parts.len() >= 7 => {
                if let (Ok(id), Some(bbox), Ok(temp_k)) = (
                    parts[1].parse::<u32>(),
                    parse_bbox(parts[2]),
                    parts[4].parse::<f64>(),
                ) {
                    energy_profiles.push(EnergyProfile {
                        id,
                        bbox,
                        unit_kj: parts[3] == "1",
                        temp_k,
                        pathway_a: parse_states(parts[5]),
                        pathway_b: parse_states(parts[6]),
                    });
                }
            }
            Some("TBL") if parts.len() >= 6 => {
                if let (Ok(id), Ok(rows), Ok(cols), Some(bbox)) = (
                    parts[1].parse::<u32>(),
                    parts[2].parse::<usize>(),
                    parts[3].parse::<usize>(),
                    parse_bbox(parts[4]),
                ) {
                    if rows > 256 || cols > 256 || rows.saturating_mul(cols) > 16_384 {
                        continue;
                    }
                    let cells: Vec<Vec<String>> = if parts[5].is_empty() {
                        vec![vec![String::new(); cols]; rows]
                    } else {
                        parts[5]
                            .split(';')
                            .take(rows)
                            .map(|r| r.split(',').take(cols).map(unesc_field).collect())
                            .collect()
                    };
                    let mut tbl = crate::SceneTable { id, rows, cols, bbox, cells };
                    tbl.ensure_cells();
                    parsed_tables.push(tbl);
                }
            }
            Some("TLC") if parts.len() >= 8 => {
                if let (Ok(id), Some(bbox), Ok(orig)) = (
                    parts[1].parse::<u32>(),
                    parse_bbox(parts[2]),
                    parts[3].parse::<f64>(),
                ) {
                    let spot_rfs: Vec<f64> = parts[4].split(',').filter_map(|v| v.parse().ok()).collect();
                    let spot_labels: Vec<String> = if parts[5].is_empty() {
                        Vec::new()
                    } else {
                        parts[5].split(',').map(unesc_field).collect()
                    };
                    let spot_colors: Vec<(u8, u8, u8)> = parts[6]
                        .split(',')
                        .filter_map(|rgb| {
                            let c: Vec<u8> = rgb.split(':').filter_map(|v| v.parse().ok()).collect();
                            if c.len() == 3 { Some((c[0], c[1], c[2])) } else { None }
                        })
                        .collect();
                    let mut tlc = crate::SceneTlcPlate {
                        id,
                        bbox,
                        origin_fraction: orig,
                        spot_rfs,
                        spot_labels,
                        spot_colors,
                        solvent_label: unesc_field(parts[7]),
                    };
                    tlc.ensure_metadata();
                    parsed_tlcs.push(tlc);
                }
            }
            _ => {}
        }
    }
    if !parsed_tables.is_empty() {
        *tables = parsed_tables;
    }
    if !parsed_tlcs.is_empty() {
        *tlc_plates = parsed_tlcs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chiralcraft_chem::analyze_organometallic_complexes;

    #[test]
    fn test_energy_profile_eyring_latex_eln_and_organometallics() {
        // 1. Eyring kinetics & Boltzmann selectivity
        let (k, thalf) = eyring_rate_and_halflife(18.0, 298.15);
        assert!(k > 0.1 && k < 10.0, "18 kcal/mol at 298 K should be ~0.4 s^-1, got {k}");
        assert!(thalf > 0.05 && thalf < 10.0);
        let (maj, ee) = boltzmann_selectivity(1.364, 298.15);
        assert!((maj - 90.9).abs() < 1.0 && (ee - 81.8).abs() < 1.5);

        // 2. Organometallic catalysts satisfy 18e- / 16e- rules!
        let pd = build_pd_pph3_4(14.4);
        assert_eq!(analyze_organometallic_complexes(&pd)[0].total_valence_electrons, 18);
        let fe = build_ferrocene_sandwich(14.4);
        assert_eq!(analyze_organometallic_complexes(&fe)[0].total_valence_electrons, 18);
        let rh = build_wilkinson_catalyst(14.4);
        assert_eq!(analyze_organometallic_complexes(&rh)[0].total_valence_electrons, 16);

        // 3. LaTeX \chemfig & TikZ export
        let benz = smiles_to_molecule("c1ccccc1O").unwrap();
        let cf = molecule_to_chemfig(&benz);
        assert!(cf.starts_with("\\chemfig{") && cf.contains("OH"));
        let mut doc = InteractiveDocument::new_saved_seed("Demo", crate::DocumentSettings::acs_document_1996());
        doc.scene.molecules.push(benz.clone());
        let tikz = document_to_tikz(&doc);
        assert!(tikz.contains("\\begin{tikzpicture}"));

        // 4. ELN GHS lookup & SI procedure generation
        assert!(!lookup_lab_reagent("THF").is_empty());
        let si = generate_si_experimental_block(&benz, "Phenol", &[], "THF", 5.0, "rt", "2 h", 92.0);
        assert!(si.contains("HRMS (ESI+)") && si.contains("IR (neat"));

        // 5. Demo screening library
        let lib = builtin_demo_screening_library(14.4);
        assert_eq!(lib.len(), 10);

        // 6. Full Lossless Extension Persistence (Energy Profile + Table Cells + Colored TLC)
        //    across CDXML -> Editable SVG -> Editable PNG -> Editable PDF!
        let mut ext_doc = InteractiveDocument::new_saved_seed("FullExt", crate::DocumentSettings::acs_document_1996());
        ext_doc.scene.molecules.push(crate::build_adamantane_2d(14.4));
        ext_doc.scene.energy_profiles.push(EnergyProfile::starter(501, 40.0, 120.0));
        ext_doc.scene.tables.push(crate::SceneTable {
            id: 502,
            rows: 2,
            cols: 2,
            bbox: (20.0, 20.0, 120.0, 70.0),
            cells: vec![vec!["Ligand".into(), "ee (%)".into()], vec!["BINAP".into(), "98%".into()]],
        });
        ext_doc.scene.tlc_plates.push(crate::SceneTlcPlate {
            id: 503,
            bbox: (160.0, 20.0, 205.0, 100.0),
            origin_fraction: 0.16,
            spot_rfs: vec![0.22, 0.74],
            spot_labels: vec!["SM".into(), "Prod".into()],
            spot_colors: vec![(20, 20, 20), (205, 45, 45)],
            solvent_label: "EtOAc/Hex 1:3".into(),
        });

        let svg = ext_doc.to_svg_transparent();
        assert!(!svg.contains("NaN") && svg.contains("BINAP") && svg.contains("TS1") && svg.contains("EtOAc/Hex 1:3"));
        let from_svg = InteractiveDocument::from_svg(&svg).expect("reopen svg with extensions");
        assert_eq!(from_svg.scene.energy_profiles.len(), 1);
        assert_eq!(from_svg.scene.energy_profiles[0].pathway_a.len(), 5);
        assert_eq!(from_svg.scene.tables[0].cells[1][0], "BINAP");
        assert_eq!(from_svg.scene.tlc_plates[0].solvent_label, "EtOAc/Hex 1:3");
        assert_eq!(from_svg.scene.tlc_plates[0].spot_colors[1], (205, 45, 45));

        let pdf = from_svg.to_editable_pdf(false, 20.0);
        let from_pdf = InteractiveDocument::from_pdf(&pdf).expect("reopen pdf with extensions");
        assert_eq!(from_pdf.scene.energy_profiles.len(), 1);
        assert_eq!(from_pdf.scene.tables[0].cells[1][1], "98%");

        // 7. Visual SVG Benchmark Gallery Audit (zero NaN, valid XML structure across complex showcase figures)
        let benchmarks: Vec<(&str, ChemMolecule)> = vec![
            ("Adamantane", crate::build_adamantane_2d(14.4)),
            ("Haworth", crate::build_haworth_pyranose(14.4)),
            ("Porphine", crate::build_porphine_core(14.4)),
            ("Crown18", crate::build_crown_ether_18_6(14.4)),
            ("Ferrocene", build_ferrocene_sandwich(14.4)),
            ("PdCatalyst", build_pd_pph3_4(14.4)),
        ];
        for (name, mol) in benchmarks {
            let mut bdoc = InteractiveDocument::new_saved_seed(name, crate::DocumentSettings::acs_document_1996());
            bdoc.scene.molecules.push(mol);
            let s_light = bdoc.to_svg_styled(true, false, 0.45, 3.6);
            let s_dark = bdoc.to_svg_styled(true, true, 0.45, 3.6);
            assert!(s_light.starts_with("<svg ") && s_light.contains("</svg>") && !s_light.contains("NaN"), "{name} light SVG valid");
            assert!(s_dark.contains("#F2F5F9") && !s_dark.contains("NaN"), "{name} dark SVG valid");
        }
    }
}
