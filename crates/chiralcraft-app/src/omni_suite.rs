//! All-Discipline Omni-Suite Studio Windows:
//! 1. Multi-Spectroscopy Studio (FT-IR %T, UV-Vis λmax & Color Swatch, Interactive 2D NMR: COSY / HSQC / HMBC)
//! 2. Organometallic 18-Electron (Green CBC) & Polymer / BigSMILES Studio
//! 3. Reaction Energy Profile (ΔG‡) & Eyring / Boltzmann Kinetics Builder
//! 4. ELN Procedure Generator & Lab Solvent/Reagent Safety (GHS + Liquid Volume mL) Handbook
//! 5. Native LaTeX (`\chemfig` / `\schemestart`) & TikZ Exporter
//! 6. Multi-Compound SDF / SMILES Virtual Screening & Library Browser

use chiralcraft_chem::{
    analyze_organometallic_complexes, compute_polymer_metrics, predict_2d_nmr, predict_ir_spectrum,
    predict_uv_vis, ChemMolecule, CopolymerTopology, Nmr2dKind,
};
use chiralcraft_doc::{
    boltzmann_selectivity, build_cisplatin, build_ferrocene_sandwich, build_grubbs_core,
    build_pd_pph3_4, build_vaska_complex, build_wilkinson_catalyst, builtin_demo_screening_library,
    document_to_chemfig_scheme, document_to_tikz, export_screening_sdf, eyring_rate_and_halflife,
    format_halflife, generate_si_experimental_block, load_screening_library, lookup_lab_reagent,
    molecule_to_chemfig, EnergyProfile, EnergyState, InteractiveDocument, ScreeningRecord,
    StoichiometryRow,
};
use chiralcraft_formats::molecule_to_smiles;
use egui::{Color32, FontId, Pos2, RichText, Sense, Stroke, Vec2};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpecTab {
    FtIr,
    UvVis,
    Nmr2d,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InoPolyTab {
    Organometallic,
    Polymer,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LatexMode {
    ChemfigSingle,
    ChemfigScheme,
    TikzVector,
}

pub struct OmniSuiteState {
    // 1. Multi-Spectroscopy
    pub show_spectroscopy: bool,
    pub spec_tab: SpecTab,
    pub nmr2d_kind: Nmr2dKind,
    pub highlighted_spec_atoms: Vec<usize>,

    // 2. Organometallic & Polymer
    pub show_ino_poly: bool,
    pub ino_poly_tab: InoPolyTab,
    pub poly_dp: f64,
    pub poly_pdi: f64,
    pub poly_end_mw: f64,
    pub poly_comon_mw: f64,
    pub poly_mol_frac: f64,
    pub poly_topology: CopolymerTopology,

    // 3. Reaction Energy Profile & Eyring/Boltzmann
    pub show_energy_profile: bool,
    pub draft_profile: EnergyProfile,
    pub show_pathway_b: bool,

    // 4. ELN Procedure & Lab Safety / GHS Handbook
    pub show_eln_safety: bool,
    pub ghs_query: String,
    pub liq_calc_mmol: f64,
    pub eln_solvent: String,
    pub eln_solvent_ml: f64,
    pub eln_temp: String,
    pub eln_time: String,
    pub eln_yield: f64,
    pub eln_generated_si: String,

    // 5. LaTeX & TikZ Exporter
    pub show_latex_export: bool,
    pub latex_mode: LatexMode,

    // 6. Virtual Screening / SDF Browser
    pub show_screening: bool,
    pub screening_path: String,
    pub screening_records: Vec<ScreeningRecord>,
    pub screen_query: String,
    pub screen_max_mw: f64,
    pub screen_max_logp: f64,
    pub screen_min_mpo: f64,
    pub screen_lipinski_only: bool,
    pub screen_no_pains_only: bool,
}

impl Default for OmniSuiteState {
    fn default() -> Self {
        Self {
            show_spectroscopy: false,
            spec_tab: SpecTab::FtIr,
            nmr2d_kind: Nmr2dKind::Hsqc,
            highlighted_spec_atoms: Vec::new(),

            show_ino_poly: false,
            ino_poly_tab: InoPolyTab::Organometallic,
            poly_dp: 100.0,
            poly_pdi: 1.15,
            poly_end_mw: 2.02,
            poly_comon_mw: 104.15,
            poly_mol_frac: 0.50,
            poly_topology: CopolymerTopology::Homopolymer,

            show_energy_profile: false,
            draft_profile: EnergyProfile::starter(1, 140.0, 220.0),
            show_pathway_b: true,

            show_eln_safety: false,
            ghs_query: String::new(),
            liq_calc_mmol: 1.0,
            eln_solvent: "THF".to_string(),
            eln_solvent_ml: 5.0,
            eln_temp: "80 °C".to_string(),
            eln_time: "12 h".to_string(),
            eln_yield: 91.0,
            eln_generated_si: String::new(),

            show_latex_export: false,
            latex_mode: LatexMode::ChemfigSingle,

            show_screening: false,
            screening_path: String::new(),
            screening_records: Vec::new(),
            screen_query: String::new(),
            screen_max_mw: 600.0,
            screen_max_logp: 6.0,
            screen_min_mpo: 0.0,
            screen_lipinski_only: false,
            screen_no_pains_only: false,
        }
    }
}

pub enum OmniSuiteAction {
    None,
    StampMolecule(ChemMolecule, &'static str),
    StampEnergyProfile(EnergyProfile),
    PasteCaption(String),
    CopyText(String, String),
    ExportFile { path: String, content: String },
    StatusMsg(String),
}

// 1. Multi-Spectroscopy Studio (FT-IR, UV-Vis, 2D NMR)
pub fn show_spectroscopy_window(
    ctx: &egui::Context,
    state: &mut OmniSuiteState,
    active_mol: Option<&ChemMolecule>,
) -> OmniSuiteAction {
    if !state.show_spectroscopy {
        state.highlighted_spec_atoms.clear();
        return OmniSuiteAction::None;
    }
    let mut act = OmniSuiteAction::None;
    let mut open = state.show_spectroscopy;

    egui::Window::new("Multi-Spectroscopy Studio (FT-IR / UV-Vis / 2D NMR)")
        .open(&mut open)
        .default_pos([260.0, 80.0])
        .default_width(560.0)
        .resizable(true)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                for (lbl, t) in [
                    ("FT-IR (4000-400 cm^-1)", SpecTab::FtIr),
                    ("UV-Vis & Color", SpecTab::UvVis),
                    ("2D NMR (COSY / HSQC / HMBC)", SpecTab::Nmr2d),
                ] {
                    if ui.selectable_label(state.spec_tab == t, lbl).clicked() {
                        state.spec_tab = t;
                        state.highlighted_spec_atoms.clear();
                    }
                }
            });
            ui.separator();

            let Some(mol) = active_mol.filter(|m| !m.atoms.is_empty()) else {
                state.highlighted_spec_atoms.clear();
                ui.weak("Draw or select a molecule to simulate IR, UV-Vis, and 2D NMR spectra.");
                return;
            };

            match state.spec_tab {
                SpecTab::FtIr => {
                    let ir = predict_ir_spectrum(mol);
                    let (rect, _) = ui.allocate_exact_size(Vec2::new(535.0, 155.0), Sense::hover());
                    let p = ui.painter();
                    p.rect_filled(rect, 2.0, Color32::WHITE);
                    p.rect_stroke(rect, 2.0, Stroke::new(1.0_f32, Color32::from_gray(145)), egui::StrokeKind::Inside);
                    let plot = rect.shrink2(Vec2::new(18.0, 18.0));
                    let wn_to_x = |wn: f64| plot.min.x + (((4000.0 - wn) / 3600.0).clamp(0.0, 1.0) as f32) * plot.width();
                    let t_to_y = |tr: f64| plot.max.y - ((tr / 100.0).clamp(0.0, 1.0) as f32) * plot.height();
                    p.line_segment([Pos2::new(plot.min.x, plot.max.y), Pos2::new(plot.max.x, plot.max.y)], Stroke::new(1.0_f32, Color32::BLACK));
                    for wn in [4000.0, 3500.0, 3000.0, 2500.0, 2000.0, 1500.0, 1000.0, 500.0] {
                        let x = wn_to_x(wn);
                        p.line_segment([Pos2::new(x, plot.max.y), Pos2::new(x, plot.max.y + 4.0)], Stroke::new(1.0_f32, Color32::BLACK));
                        p.text(Pos2::new(x, plot.max.y + 5.0), egui::Align2::CENTER_TOP, format!("{wn:.0}"), FontId::proportional(8.5), Color32::BLACK);
                    }
                    let curve = ir.simulate_transmittance(320);
                    let pts: Vec<Pos2> = curve.iter().map(|&(wn, tr)| Pos2::new(wn_to_x(wn), t_to_y(tr))).collect();
                    for w in pts.windows(2) {
                        p.line_segment([w[0], w[1]], Stroke::new(1.25_f32, Color32::from_rgb(18, 68, 150)));
                    }
                    ui.small("Click any IR band to highlight its vibrating functional group on canvas:");
                    egui::ScrollArea::vertical().max_height(130.0).show(ui, |ui| {
                        for b in &ir.bands {
                            let is_sel = !b.atom_indices.is_empty() && state.highlighted_spec_atoms == b.atom_indices;
                            if ui.selectable_label(is_sel, format!("{:.0} cm^-1 ({}) — {}", b.wavenumber_cm1, b.strength, b.assignment)).clicked() {
                                state.highlighted_spec_atoms = if is_sel { Vec::new() } else { b.atom_indices.clone() };
                            }
                        }
                    });
                    if ui.button("Paste IR Line on Canvas").clicked() {
                        act = OmniSuiteAction::PasteCaption(ir.report_text);
                    }
                }
                SpecTab::UvVis => {
                    let uv = predict_uv_vis(mol);
                    ui.horizontal(|ui| {
                        let (sw, _) = ui.allocate_exact_size(Vec2::new(42.0, 42.0), Sense::hover());
                        ui.painter().rect_filled(sw, 6.0, Color32::from_rgb(uv.perceived_rgb.0, uv.perceived_rgb.1, uv.perceived_rgb.2));
                        ui.painter().rect_stroke(sw, 6.0, Stroke::new(1.0_f32, Color32::GRAY), egui::StrokeKind::Inside);
                        ui.vertical(|ui| {
                            ui.strong(format!("λmax ≈ {:.0} nm   (ε ≈ {:.0} M⁻¹cm^-1)", uv.lambda_max_nm, uv.epsilon_est));
                            ui.label(format!("Transition: {}  |  Chromophore: {}", uv.transition, uv.chromophore));
                            ui.label(format!("Perceived Solution Color: {}", uv.perceived_color_name));
                        });
                    });
                    // Mini UV-Vis absorption band plot (190 - 700 nm)
                    let (rect, _) = ui.allocate_exact_size(Vec2::new(535.0, 120.0), Sense::hover());
                    let p = ui.painter();
                    p.rect_filled(rect, 2.0, Color32::WHITE);
                    p.rect_stroke(rect, 2.0, Stroke::new(1.0_f32, Color32::from_gray(145)), egui::StrokeKind::Inside);
                    let plot = rect.shrink2(Vec2::new(18.0, 18.0));
                    let to_x = |nm: f64| plot.min.x + (((nm - 180.0) / 520.0).clamp(0.0, 1.0) as f32) * plot.width();
                    // Visible spectrum rainbow bar at bottom (380-700 nm)
                    for step in 0..32 {
                        let nm = 380.0 + (step as f64) * 10.0;
                        let col = match step {
                            0..=4 => Color32::from_rgb(130, 60, 210),
                            5..=10 => Color32::from_rgb(50, 120, 240),
                            11..=16 => Color32::from_rgb(45, 195, 110),
                            17..=21 => Color32::from_rgb(245, 220, 45),
                            22..=26 => Color32::from_rgb(245, 130, 35),
                            _ => Color32::from_rgb(215, 45, 45),
                        };
                        let r = egui::Rect::from_min_max(Pos2::new(to_x(nm), plot.max.y - 4.0), Pos2::new(to_x(nm + 10.0), plot.max.y));
                        p.rect_filled(r, 0.0, col);
                    }
                    let mut pts = Vec::new();
                    for k in 0..=180 {
                        let nm = 180.0 + (k as f64) * (520.0 / 180.0);
                        let z = (nm - uv.lambda_max_nm) / 28.0;
                        let y = (-0.5 * z * z).exp();
                        pts.push(Pos2::new(to_x(nm), plot.max.y - 6.0 - (y as f32) * (plot.height() - 10.0)));
                    }
                    for w in pts.windows(2) {
                        p.line_segment([w[0], w[1]], Stroke::new(1.5_f32, Color32::from_rgb(160, 30, 110)));
                    }
                    for nm in [200.0, 300.0, 400.0, 500.0, 600.0, 700.0] {
                        p.text(Pos2::new(to_x(nm), plot.max.y + 3.0), egui::Align2::CENTER_TOP, format!("{nm:.0}"), FontId::proportional(8.5), Color32::BLACK);
                    }
                }
                SpecTab::Nmr2d => {
                    ui.horizontal(|ui| {
                        for (lbl, k) in [
                            ("COSY (¹H–¹H)", Nmr2dKind::Cosy),
                            ("HSQC (¹J ¹H–¹³C)", Nmr2dKind::Hsqc),
                            ("HMBC (²/³J ¹H–¹³C)", Nmr2dKind::Hmbc),
                        ] {
                            if ui.selectable_label(state.nmr2d_kind == k, lbl).clicked() {
                                state.nmr2d_kind = k;
                                state.highlighted_spec_atoms.clear();
                            }
                        }
                    });
                    let pred2d = predict_2d_nmr(mol, state.nmr2d_kind);
                    let (rect, resp) = ui.allocate_exact_size(Vec2::new(320.0, 220.0), Sense::click());
                    let p = ui.painter();
                    p.rect_filled(rect, 2.0, Color32::WHITE);
                    p.rect_stroke(rect, 2.0, Stroke::new(1.0_f32, Color32::from_gray(140)), egui::StrokeKind::Inside);
                    let plot = rect.shrink2(Vec2::new(26.0, 18.0));
                    let x_to_px = |ppm: f64| plot.min.x + (((pred2d.x_range.1 - ppm) / (pred2d.x_range.1 - pred2d.x_range.0)).clamp(0.0, 1.0) as f32) * plot.width();
                    let y_to_px = |ppm: f64| plot.min.y + (((ppm - pred2d.y_range.0) / (pred2d.y_range.1 - pred2d.y_range.0)).clamp(0.0, 1.0) as f32) * plot.height();
                    // Reversed NMR axes: high ppm top-left / bottom-left standard
                    p.rect_stroke(plot, 0.0, Stroke::new(0.9_f32, Color32::BLACK), egui::StrokeKind::Inside);
                    if state.nmr2d_kind == Nmr2dKind::Cosy {
                        p.line_segment([plot.min, plot.max], Stroke::new(0.6_f32, Color32::from_gray(200)));
                    }
                    let click_pos = if resp.clicked() { resp.interact_pointer_pos() } else { None };
                    for pk in &pred2d.peaks {
                        let pt = Pos2::new(x_to_px(pk.x_ppm), y_to_px(pk.y_ppm));
                        let sel = state.highlighted_spec_atoms.contains(&pk.atom_x) && state.highlighted_spec_atoms.contains(&pk.atom_y);
                        let col = if sel { Color32::from_rgb(220, 35, 35) } else { Color32::from_rgb(20, 95, 195) };
                        p.circle_stroke(pt, 5.0, Stroke::new(1.0_f32, col));
                        p.circle_filled(pt, 2.5, col);
                        if let Some(cp) = click_pos {
                            if cp.distance(pt) < 9.0 {
                                state.highlighted_spec_atoms = vec![pk.atom_x, pk.atom_y];
                            }
                        }
                    }
                    ui.small("Click any contour cross-peak or row below to highlight BOTH coupled atoms on canvas:");
                    egui::ScrollArea::vertical().max_height(110.0).show(ui, |ui| {
                        for pk in &pred2d.peaks {
                            let sel = state.highlighted_spec_atoms == vec![pk.atom_x, pk.atom_y];
                            if ui.selectable_label(sel, &pk.label).clicked() {
                                state.highlighted_spec_atoms = vec![pk.atom_x, pk.atom_y];
                            }
                        }
                    });
                }
            }
        });
    state.show_spectroscopy = open;
    act
}

// 2. Organometallic 18-Electron & Polymer / BigSMILES Studio
pub fn show_ino_poly_window(
    ctx: &egui::Context,
    state: &mut OmniSuiteState,
    active_mol: Option<&ChemMolecule>,
    bond_len: f64,
) -> OmniSuiteAction {
    if !state.show_ino_poly {
        return OmniSuiteAction::None;
    }
    let mut act = OmniSuiteAction::None;
    let mut open = state.show_ino_poly;

    egui::Window::new("Organometallic 18e- & Polymer / BigSMILES Studio")
        .open(&mut open)
        .default_pos([300.0, 110.0])
        .default_width(480.0)
        .resizable(true)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut state.ino_poly_tab, InoPolyTab::Organometallic, "Organometallic 18e- / Green CBC");
                ui.selectable_value(&mut state.ino_poly_tab, InoPolyTab::Polymer, "Polymer Mn/Mw/D & BigSMILES");
            });
            ui.separator();

            match state.ino_poly_tab {
                InoPolyTab::Organometallic => {
                    ui.horizontal_wrapped(|ui| {
                        ui.small("Quick catalysts:");
                        if ui.small_button("Pd(PPh₃)₄").clicked() { act = OmniSuiteAction::StampMolecule(build_pd_pph3_4(bond_len), "Pd(PPh3)4"); }
                        if ui.small_button("Cisplatin").clicked() { act = OmniSuiteAction::StampMolecule(build_cisplatin(bond_len), "Cisplatin"); }
                        if ui.small_button("Wilkinson's").clicked() { act = OmniSuiteAction::StampMolecule(build_wilkinson_catalyst(bond_len), "Wilkinson Catalyst"); }
                        if ui.small_button("Vaska's").clicked() { act = OmniSuiteAction::StampMolecule(build_vaska_complex(bond_len), "Vaska Complex"); }
                        if ui.small_button("Ferrocene").clicked() { act = OmniSuiteAction::StampMolecule(build_ferrocene_sandwich(bond_len), "Ferrocene"); }
                        if ui.small_button("Grubbs Ru").clicked() { act = OmniSuiteAction::StampMolecule(build_grubbs_core(bond_len), "Grubbs Catalyst"); }
                    });
                    ui.separator();
                    let Some(mol) = active_mol.filter(|m| !m.atoms.is_empty()) else {
                        ui.weak("Draw or stamp a transition-metal complex to analyze 18e- count & oxidation state.");
                        return;
                    };
                    let reports = analyze_organometallic_complexes(mol);
                    if reports.is_empty() {
                        ui.weak("No d-block transition metal centers found in the active molecule.");
                    } else {
                        for r in &reports {
                            ui.group(|ui| {
                                let col = if r.total_valence_electrons == 18 || (r.total_valence_electrons == 16 && r.d_electrons == 8) {
                                    Color32::from_rgb(20, 145, 70)
                                } else {
                                    Color32::from_rgb(205, 120, 15)
                                };
                                ui.label(RichText::new(format!(
                                    "{} (Group {}) — {}e-  |  Oxidation State: {:+}  |  d^{}",
                                    r.metal_symbol, r.group_number, r.total_valence_electrons, r.oxidation_state, r.d_electrons
                                )).strong().color(col));
                                ui.label(format!("Status: {}  |  Geometry: {}", r.rule_status, r.predicted_geometry));
                                ui.label(format!("Green CBC: ML{}X{}  (CN = {})", r.l_ligands, r.x_ligands, r.coordination_number));
                                for ld in &r.ligand_details {
                                    ui.small(format!("  • {ld}"));
                                }
                            });
                        }
                    }
                }
                InoPolyTab::Polymer => {
                    let Some(mol) = active_mol.filter(|m| !m.atoms.is_empty()) else {
                        ui.weak("Draw a monomer / repeat unit on canvas to compute Mn, Mw, and BigSMILES.");
                        return;
                    };
                    ui.horizontal(|ui| {
                        ui.label("Topology:");
                        for (lbl, top) in [
                            ("Homopolymer", CopolymerTopology::Homopolymer),
                            ("Random", CopolymerTopology::Random),
                            ("Block", CopolymerTopology::Block),
                            ("Alternating", CopolymerTopology::Alternating),
                        ] {
                            ui.selectable_value(&mut state.poly_topology, top, lbl);
                        }
                    });
                    egui::Grid::new("poly_inputs").num_columns(4).spacing([8.0, 4.0]).show(ui, |ui| {
                        ui.label("DP (n):"); ui.add(egui::DragValue::new(&mut state.poly_dp).range(1.0..=100_000.0));
                        ui.label("Dispersity (D):"); ui.add(egui::DragValue::new(&mut state.poly_pdi).speed(0.02).range(1.0..=10.0));
                        ui.end_row();
                        ui.label("End-Group MW:"); ui.add(egui::DragValue::new(&mut state.poly_end_mw).range(0.0..=5000.0));
                        ui.label("Comonomer MW:"); ui.add(egui::DragValue::new(&mut state.poly_comon_mw).range(12.0..=2000.0));
                        ui.end_row();
                    });
                    if matches!(state.poly_topology, CopolymerTopology::Random | CopolymerTopology::Block) {
                        ui.add(egui::Slider::new(&mut state.poly_mol_frac, 0.05..=0.95).text("Mol frac A"));
                    }
                    let smi = molecule_to_smiles(mol);
                    let rep = compute_polymer_metrics(
                        mol,
                        &smi,
                        state.poly_dp,
                        state.poly_pdi,
                        state.poly_end_mw,
                        state.poly_comon_mw,
                        state.poly_mol_frac,
                        state.poly_topology,
                    );
                    ui.separator();
                    ui.strong(format!(
                        "Repeat Unit MW: {:.2} Da  |  Mn: {:.1} Da ({:.2} kDa)  |  Mw: {:.1} Da ({:.2} kDa)",
                        rep.repeat_mw, rep.mn, rep.mn / 1000.0, rep.mw, rep.mw / 1000.0
                    ));
                    ui.label(format!("Weight Fraction A: {:.1}%  |  D = {:.2}", rep.weight_frac_a * 100.0, rep.pdi));
                    ui.horizontal(|ui| {
                        ui.label("BigSMILES:");
                        ui.label(RichText::new(&rep.big_smiles).monospace());
                        if ui.small_button("Copy").clicked() {
                            act = OmniSuiteAction::CopyText(rep.big_smiles.clone(), "Copied BigSMILES".into());
                        }
                    });
                }
            }
        });
    state.show_ino_poly = open;
    act
}

// 3. Reaction Energy Profile (ΔG‡) & Eyring / Boltzmann Kinetics Builder
pub fn show_energy_profile_window(
    ctx: &egui::Context,
    state: &mut OmniSuiteState,
) -> OmniSuiteAction {
    if !state.show_energy_profile {
        return OmniSuiteAction::None;
    }
    let mut act = OmniSuiteAction::None;
    let mut open = state.show_energy_profile;

    egui::Window::new("Reaction Energy Profile (ΔG‡) & Kinetics Builder")
        .open(&mut open)
        .default_pos([260.0, 110.0])
        .default_width(500.0)
        .resizable(true)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.checkbox(&mut state.draft_profile.unit_kj, "Use kJ/mol");
                ui.checkbox(&mut state.show_pathway_b, "Competing Pathway B");
                ui.label("T (K):");
                ui.add(egui::DragValue::new(&mut state.draft_profile.temp_k).speed(1.0).range(77.0..=1200.0));
            });

            // Live Eyring & Boltzmann Readout!
            let dg = state.draft_profile.max_barrier_kcal();
            let (k_rate, t_half) = eyring_rate_and_halflife(dg, state.draft_profile.temp_k);
            ui.label(RichText::new(format!(
                "Eyring Kinetics (ΔG‡ = {dg:.2} kcal/mol @ {:.1} K):  k = {k_rate:.2e} s^-1  |  t1/2 = {}",
                state.draft_profile.temp_k,
                format_halflife(t_half)
            )).strong().color(Color32::from_rgb(20, 95, 185)));
            if state.show_pathway_b {
                if let Some(ddg) = state.draft_profile.ddg_kcal() {
                    let (maj, ee) = boltzmann_selectivity(ddg, state.draft_profile.temp_k);
                    ui.label(RichText::new(format!(
                        "Boltzmann Selectivity (ΔΔG‡ = {ddg:.2} kcal/mol):  {maj:.1} : {:.1}  ({ee:.1}% ee)",
                        100.0 - maj
                    )).color(Color32::from_rgb(175, 45, 45)));
                }
            }
            ui.separator();

            // Pathway A state editor
            let unit_s = state.draft_profile.unit_str();
            ui.horizontal(|ui| {
                ui.strong("Pathway A (Blue Solid)");
                if ui.small_button("+ State").clicked() {
                    let idx = state.draft_profile.pathway_a.len() + 1;
                    state.draft_profile.pathway_a.push(EnergyState { label: format!("S{idx}"), energy: 0.0, is_ts: false });
                }
            });
            let mut rem_a = None;
            let len_a = state.draft_profile.pathway_a.len();
            for (i, st) in state.draft_profile.pathway_a.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.add(egui::TextEdit::singleline(&mut st.label).desired_width(85.0));
                    ui.add(egui::DragValue::new(&mut st.energy).speed(0.2).fixed_decimals(1));
                    ui.label(unit_s);
                    ui.checkbox(&mut st.is_ts, "TS");
                    if len_a > 2 && ui.small_button("×").clicked() {
                        rem_a = Some(i);
                    }
                });
            }
            if let Some(ri) = rem_a { state.draft_profile.pathway_a.remove(ri); }

            if state.show_pathway_b {
                ui.separator();
                ui.horizontal(|ui| {
                    ui.strong("Pathway B (Red Dashed)");
                    if ui.small_button("+ State").clicked() {
                        let idx = state.draft_profile.pathway_b.len() + 1;
                        state.draft_profile.pathway_b.push(EnergyState { label: format!("S{idx}'"), energy: 2.0, is_ts: false });
                    }
                });
                let mut rem_b = None;
                let len_b = state.draft_profile.pathway_b.len();
                for (i, st) in state.draft_profile.pathway_b.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        ui.add(egui::TextEdit::singleline(&mut st.label).desired_width(85.0));
                        ui.add(egui::DragValue::new(&mut st.energy).speed(0.2).fixed_decimals(1));
                        ui.label(unit_s);
                        if len_b > 2 && ui.small_button("×").clicked() {
                            rem_b = Some(i);
                        }
                    });
                }
                if let Some(ri) = rem_b { state.draft_profile.pathway_b.remove(ri); }
            }

            ui.separator();
            if ui.button("Stamp Energy Profile on Canvas (SVG/PDF Ready)").clicked() {
                let mut prof = state.draft_profile.clone();
                if !state.show_pathway_b {
                    prof.pathway_b.clear();
                }
                act = OmniSuiteAction::StampEnergyProfile(prof);
            }
        });
    state.show_energy_profile = open;
    act
}

// 4. ELN Procedure Generator & Lab Safety / GHS Handbook
pub fn show_eln_safety_window(
    ctx: &egui::Context,
    state: &mut OmniSuiteState,
    active_mol: Option<&ChemMolecule>,
    active_name: &str,
    stoich_rows: &[StoichiometryRow],
) -> OmniSuiteAction {
    if !state.show_eln_safety {
        return OmniSuiteAction::None;
    }
    let mut act = OmniSuiteAction::None;
    let mut open = state.show_eln_safety;

    egui::Window::new("ELN Supporting Info & Lab Safety / GHS Handbook")
        .open(&mut open)
        .default_pos([210.0, 90.0])
        .default_width(600.0)
        .resizable(true)
        .show(ctx, |ui| {
            egui::CollapsingHeader::new("Solvent & Reagent Safety Handbook + Liquid Volume (mL) Calculator")
                .default_open(true)
                .show(ui, |ui| {
                    ui.horizontal(|ui| {
                        ui.label("Search:");
                        ui.add(egui::TextEdit::singleline(&mut state.ghs_query).desired_width(190.0).hint_text("THF, DCM, n-BuLi, Pyrophoric..."));
                        ui.label("Scale (mmol):");
                        ui.add(egui::DragValue::new(&mut state.liq_calc_mmol).speed(0.1).range(0.01..=500.0));
                    });
                    let hits = lookup_lab_reagent(&state.ghs_query);
                    egui::ScrollArea::vertical().max_height(155.0).show(ui, |ui| {
                        egui::Grid::new("ghs_tbl").striped(true).num_columns(6).spacing([8.0, 3.0]).show(ui, |ui| {
                            for h in ["Reagent", "MW", "Density", "BP °C", "Mass / Vol", "GHS Hazards"] {
                                ui.strong(h);
                            }
                            ui.end_row();
                            for r in hits {
                                let mg = state.liq_calc_mmol * r.mw;
                                let amt = match r.density_g_ml {
                                    Some(d) => format!("{:.1} mg ({:.3} mL)", mg, (mg / 1000.0) / d),
                                    None => format!("{:.1} mg (solid)", mg),
                                };
                                let hcol = if r.ghs.contains("PYROPHORIC") || r.ghs.contains("CMR") {
                                    Color32::from_rgb(205, 40, 40)
                                } else if r.ghs.contains("Corrosive") || r.ghs.contains("Toxic") {
                                    Color32::from_rgb(200, 115, 15)
                                } else {
                                    Color32::from_gray(95)
                                };
                                ui.label(format!("{} ({})", r.abbrev, r.name));
                                ui.label(format!("{:.2}", r.mw));
                                ui.label(r.density_g_ml.map(|d| format!("{d:.3}")).unwrap_or_else(|| "—".into()));
                                ui.label(r.bp_c.map(|b| format!("{b:.0}")).unwrap_or_else(|| "—".into()));
                                ui.label(amt);
                                ui.label(RichText::new(r.ghs).small().color(hcol));
                                ui.end_row();
                            }
                        });
                    });
                });

            ui.separator();
            ui.strong("ACS / JOC Supporting Information Experimental Generator");
            ui.horizontal_wrapped(|ui| {
                ui.label("Solvent:"); ui.add(egui::TextEdit::singleline(&mut state.eln_solvent).desired_width(60.0));
                ui.label("mL:"); ui.add(egui::DragValue::new(&mut state.eln_solvent_ml).range(0.1..=500.0));
                ui.label("Temp:"); ui.add(egui::TextEdit::singleline(&mut state.eln_temp).desired_width(55.0));
                ui.label("Time:"); ui.add(egui::TextEdit::singleline(&mut state.eln_time).desired_width(45.0));
                ui.label("Yield %:"); ui.add(egui::DragValue::new(&mut state.eln_yield).range(1.0..=100.0));
            });
            ui.horizontal(|ui| {
                if ui.button("Generate SI Procedure + Characterization").clicked() {
                    if let Some(m) = active_mol.filter(|m| !m.atoms.is_empty()) {
                        state.eln_generated_si = generate_si_experimental_block(
                            m,
                            active_name,
                            stoich_rows,
                            &state.eln_solvent,
                            state.eln_solvent_ml,
                            &state.eln_temp,
                            &state.eln_time,
                            state.eln_yield,
                        );
                    } else {
                        act = OmniSuiteAction::StatusMsg("Draw or select the product molecule first".into());
                    }
                }
                if ui.add_enabled(!state.eln_generated_si.is_empty(), egui::Button::new("Copy to Clipboard")).clicked() {
                    act = OmniSuiteAction::CopyText(state.eln_generated_si.clone(), "Copied ACS/JOC Supporting Info paragraph".into());
                }
            });
            if !state.eln_generated_si.is_empty() {
                egui::ScrollArea::vertical().max_height(145.0).show(ui, |ui| {
                    ui.add(egui::TextEdit::multiline(&mut state.eln_generated_si).desired_width(f32::INFINITY).desired_rows(7));
                });
            }
        });
    state.show_eln_safety = open;
    act
}

// 5. Native LaTeX (\chemfig / \schemestart) & TikZ Exporter
pub fn show_latex_export_window(
    ctx: &egui::Context,
    state: &mut OmniSuiteState,
    doc: &InteractiveDocument,
    active_mol: Option<&ChemMolecule>,
) -> OmniSuiteAction {
    if !state.show_latex_export {
        return OmniSuiteAction::None;
    }
    let mut act = OmniSuiteAction::None;
    let mut open = state.show_latex_export;

    egui::Window::new("LaTeX (\\chemfig) & TikZ Code Exporter")
        .open(&mut open)
        .default_pos([320.0, 120.0])
        .default_width(520.0)
        .resizable(true)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut state.latex_mode, LatexMode::ChemfigSingle, "\\chemfig (Active)");
                ui.selectable_value(&mut state.latex_mode, LatexMode::ChemfigScheme, "\\schemestart (Full Scheme)");
                ui.selectable_value(&mut state.latex_mode, LatexMode::TikzVector, "TikZ Vector");
            });
            let mut code = match state.latex_mode {
                LatexMode::ChemfigSingle => active_mol.map(molecule_to_chemfig).unwrap_or_else(|| "% Draw a molecule first".into()),
                LatexMode::ChemfigScheme => document_to_chemfig_scheme(doc),
                LatexMode::TikzVector => document_to_tikz(doc),
            };
            egui::ScrollArea::vertical().max_height(210.0).show(ui, |ui| {
                ui.add(egui::TextEdit::multiline(&mut code).font(egui::TextStyle::Monospace).desired_width(f32::INFINITY).desired_rows(10));
            });
            ui.horizontal(|ui| {
                if ui.button("Copy LaTeX Code").clicked() {
                    act = OmniSuiteAction::CopyText(code.clone(), "Copied LaTeX code to clipboard".into());
                }
                if ui.button("Save to /tmp/chiralcraft-export.tex").clicked() {
                    act = OmniSuiteAction::ExportFile {
                        path: "/tmp/chiralcraft-export.tex".into(),
                        content: code,
                    };
                }
            });
        });
    state.show_latex_export = open;
    act
}

// 6. Multi-Compound SDF / SMILES Virtual Screening & Library Browser
pub fn show_screening_window(
    ctx: &egui::Context,
    state: &mut OmniSuiteState,
    bond_len: f64,
) -> OmniSuiteAction {
    if !state.show_screening {
        return OmniSuiteAction::None;
    }
    let mut act = OmniSuiteAction::None;
    let mut open = state.show_screening;

    egui::Window::new("Multi-Compound SDF / SMILES Virtual Screening Browser")
        .open(&mut open)
        .default_pos([180.0, 85.0])
        .default_width(680.0)
        .resizable(true)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("File (.sdf/.smi/.csv):");
                ui.add(egui::TextEdit::singleline(&mut state.screening_path).desired_width(240.0));
                if ui.button("Load File").clicked() && !state.screening_path.trim().is_empty() {
                    match load_screening_library(state.screening_path.trim(), bond_len) {
                        Ok(recs) => {
                            let cnt = recs.len();
                            state.screening_records = recs;
                            act = OmniSuiteAction::StatusMsg(format!("Loaded {cnt} compounds for virtual screening"));
                        }
                        Err(e) => act = OmniSuiteAction::StatusMsg(format!("Library load error: {e}")),
                    }
                }
                if ui.button("Load Demo Drug Set (10)").clicked() {
                    state.screening_records = builtin_demo_screening_library(bond_len);
                }
            });
            ui.horizontal_wrapped(|ui| {
                ui.label("Filter:");
                ui.add(egui::TextEdit::singleline(&mut state.screen_query).desired_width(110.0).hint_text("name / SMILES"));
                ui.label("Max MW:"); ui.add(egui::DragValue::new(&mut state.screen_max_mw).range(50.0..=2000.0));
                ui.label("Max LogP:"); ui.add(egui::DragValue::new(&mut state.screen_max_logp).speed(0.2).range(-5.0..=12.0));
                ui.label("Min MPO:"); ui.add(egui::DragValue::new(&mut state.screen_min_mpo).speed(0.2).range(0.0..=6.0));
                ui.checkbox(&mut state.screen_lipinski_only, "Lipinski PASS");
                ui.checkbox(&mut state.screen_no_pains_only, "PAINS = 0");
            });
            ui.separator();

            let q = state.screen_query.trim().to_ascii_lowercase();
            let filtered: Vec<&ScreeningRecord> = state
                .screening_records
                .iter()
                .filter(|r| {
                    (q.is_empty() || r.title.to_ascii_lowercase().contains(&q) || r.smiles.to_ascii_lowercase().contains(&q))
                        && r.medchem.mw <= state.screen_max_mw
                        && r.medchem.log_p <= state.screen_max_logp
                        && r.medchem.cns_mpo >= state.screen_min_mpo
                        && (!state.screen_lipinski_only || r.medchem.lipinski_pass)
                        && (!state.screen_no_pains_only || r.medchem.alerts.is_empty())
                })
                .collect();

            ui.label(format!("Showing {} of {} compound(s)", filtered.len(), state.screening_records.len()));
            egui::ScrollArea::vertical().max_height(225.0).show(ui, |ui| {
                egui::Grid::new("screen_tbl").striped(true).num_columns(9).spacing([8.0, 3.0]).show(ui, |ui| {
                    for h in ["#", "Name", "Formula", "MW", "LogP", "MPO", "QED", "Alerts", "Canvas"] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for r in &filtered {
                        ui.label(format!("{}", r.index));
                        ui.label(&r.title);
                        ui.label(&r.formula);
                        ui.label(format!("{:.1}", r.medchem.mw));
                        ui.label(format!("{:.2}", r.medchem.log_p));
                        ui.label(format!("{:.2}", r.medchem.cns_mpo));
                        ui.label(format!("{:.2}", r.medchem.qed));
                        let al_cnt = r.medchem.alerts.len();
                        ui.label(RichText::new(format!("{al_cnt}")).color(if al_cnt == 0 { Color32::from_rgb(20, 145, 70) } else { Color32::from_rgb(205, 45, 45) }));
                        if ui.small_button("Open").clicked() {
                            act = OmniSuiteAction::StampMolecule(r.molecule.clone(), "Open Screened Compound");
                        }
                        ui.end_row();
                    }
                });
            });
            ui.separator();
            if ui.add_enabled(!filtered.is_empty(), egui::Button::new("Export Filtered Hits to /tmp/chiralcraft-hits.sdf")).clicked() {
                act = OmniSuiteAction::ExportFile {
                    path: "/tmp/chiralcraft-hits.sdf".into(),
                    content: export_screening_sdf(&filtered),
                };
            }
        });
    state.show_screening = open;
    act
}
