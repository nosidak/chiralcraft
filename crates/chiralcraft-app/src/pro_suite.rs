//! Flagship Medicinal Chemistry, HRMS, Isomer/Resonance, Markush SAR,
//! Retrosynthesis, and Publication Figure Composer Floating Studio Windows.

use chiralcraft_chem::{
    apply_electron_push_arrow, compute_medchem_report, compute_microspecies_at_ph,
    enumerate_resonance_structures, enumerate_stereoisomers, enumerate_tautomers, simulate_hrms,
    AlertSeverity, ChemMolecule,
};
use chiralcraft_doc::{
    count_molecules_overflowing_frame, enumerate_rgroup_library, export_sar_csv, export_sar_sdf,
    suggest_disconnections, DocumentScene, EnumeratedAnalog, RetrosyntheticHit, TocPreset,
};
use egui::{Color32, FontId, Pos2, RichText, Sense, Stroke, Vec2};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsomerTab {
    PhMicrospecies,
    Tautomers,
    Stereoisomers,
    Resonance,
    ElectronPush,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SarSortKey {
    Id,
    Mw,
    LogP,
    CnsMpo,
    Qed,
    Alerts,
}

pub struct ProSuiteState {
    // 1. MedChem / PAINS
    pub show_medchem: bool,
    pub highlight_pains_on_canvas: bool,
    pub active_pains_atoms: Vec<usize>,

    // 2. pH Microspecies / Tautomers / Stereoisomers / Resonance / Electron Push
    pub show_isomers: bool,
    pub isomer_tab: IsomerTab,
    pub target_ph: f64,
    pub push_from_atom: usize,
    pub push_to_atom: usize,
    pub push_feedback: String,

    // 3. HRMS Isotope Envelope & Adducts
    pub show_hrms: bool,
    pub hrms_adduct: String,
    pub hrms_resolution: f64,
    pub hrms_observed_mz: String,

    // 4. Markush R-Group Enumerator & SAR Grid
    pub show_markush: bool,
    pub r1_input: String,
    pub r2_input: String,
    pub sar_analogs: Vec<EnumeratedAnalog>,
    pub sar_sort: SarSortKey,

    // 5. Retrosynthesis Advisor
    pub show_retro: bool,
    pub retro_hits: Vec<RetrosyntheticHit>,
    pub highlighted_retro_bond: Option<usize>,

    // 6. Publication TOC & Scheme Composer
    pub show_pub_composer: bool,
    pub selected_toc_idx: usize,
    pub last_toc_bbox: Option<(f64, f64, f64, f64)>,
    pub cond_above: String,
    pub cond_below: String,
    pub serpentine_per_row: usize,

    // 7. JCAMP-DX Experimental NMR Overlay path input
    pub jcamp_path_buf: String,
}

impl Default for ProSuiteState {
    fn default() -> Self {
        Self {
            show_medchem: false,
            highlight_pains_on_canvas: true,
            active_pains_atoms: Vec::new(),

            show_isomers: false,
            isomer_tab: IsomerTab::PhMicrospecies,
            target_ph: 7.4,
            push_from_atom: 1,
            push_to_atom: 2,
            push_feedback: String::new(),

            show_hrms: false,
            hrms_adduct: "[M+H]+".to_string(),
            hrms_resolution: 30_000.0,
            hrms_observed_mz: String::new(),

            show_markush: false,
            r1_input: "H, Me, F, Cl, OMe, CF3, CN".to_string(),
            r2_input: "H, Me, Ph".to_string(),
            sar_analogs: Vec::new(),
            sar_sort: SarSortKey::Id,

            show_retro: false,
            retro_hits: Vec::new(),
            highlighted_retro_bond: None,

            show_pub_composer: false,
            selected_toc_idx: 0,
            last_toc_bbox: None,
            cond_above: "Pd(PPh3)4 (5 mol%), K2CO3".to_string(),
            cond_below: "THF/H2O, 80 °C, 12 h, 92%".to_string(),
            serpentine_per_row: 3,

            jcamp_path_buf: String::new(),
        }
    }
}

pub enum ProSuiteAction {
    None,
    ReplaceActiveMolecule(ChemMolecule, &'static str),
    StampMolecule(ChemMolecule, &'static str),
    StampMoleculeGrid(Vec<(String, ChemMolecule)>, &'static str),
    StampRetrosynthesis(RetrosyntheticHit),
    InsertTocFrame(TocPreset),
    StampArrowConditions { above: String, below: String },
    WrapSerpentine(usize),
    ExportTextFile { suggested_name: &'static str, content: String },
    StatusMsg(String),
}

// 1. MedChem / ADMET & PAINS Scorecard Window
pub fn show_medchem_window(
    ctx: &egui::Context,
    state: &mut ProSuiteState,
    active_mol: Option<&ChemMolecule>,
) {
    if !state.show_medchem {
        state.active_pains_atoms.clear();
        return;
    }
    let mut open = state.show_medchem;
    egui::Window::new("Medicinal Chemistry, ADMET & PAINS Scorecard")
        .open(&mut open)
        .default_pos([840.0, 70.0])
        .default_width(375.0)
        .resizable(true)
        .show(ctx, |ui| {
            let Some(mol) = active_mol.filter(|m| !m.atoms.is_empty()) else {
                state.active_pains_atoms.clear();
                ui.weak("Draw or select a molecule to inspect drug-likeness & PAINS alerts.");
                return;
            };
            let rep = compute_medchem_report(mol);
            // Collect all alert atoms for live canvas glow
            let mut alert_atoms: Vec<usize> = rep.alerts.iter().flat_map(|a| a.atom_indices.iter().copied()).collect();
            alert_atoms.sort_unstable();
            alert_atoms.dedup();
            state.active_pains_atoms = if state.highlight_pains_on_canvas { alert_atoms } else { Vec::new() };

            // Rule badges
            ui.horizontal(|ui| {
                let badge = |ui: &mut egui::Ui, name: &str, ok: bool| {
                    let col = if ok { Color32::from_rgb(20, 140, 70) } else { Color32::from_rgb(200, 45, 45) };
                    ui.label(RichText::new(format!("{} {}", if ok { "✓" } else { "✗" }, name)).strong().color(col));
                };
                badge(ui, &format!("Lipinski ({} viol)", rep.lipinski_violations), rep.lipinski_pass);
                ui.separator();
                badge(ui, "Veber", rep.veber_pass);
                ui.separator();
                badge(ui, "Ghose", rep.ghose_pass);
            });
            ui.separator();

            // Key composite scores
            ui.horizontal(|ui| {
                ui.label(format!("CNS MPO: {:.2} / 6.0", rep.cns_mpo));
                ui.add(egui::ProgressBar::new((rep.cns_mpo / 6.0) as f32).desired_width(140.0));
            });
            ui.horizontal(|ui| {
                ui.label(format!("QED Drug-Likeness: {:.2}", rep.qed));
                ui.add(egui::ProgressBar::new(rep.qed as f32).desired_width(140.0));
            });
            ui.add_space(4.0);

            egui::Grid::new("medchem_metrics").num_columns(4).spacing([14.0, 3.0]).show(ui, |ui| {
                ui.strong("MW:"); ui.label(format!("{:.2}", rep.mw));
                ui.strong("LogP:"); ui.label(format!("{:.2}", rep.log_p));
                ui.end_row();
                ui.strong("tPSA:"); ui.label(format!("{:.1} Å²", rep.tpsa));
                ui.strong("HBD / HBA:"); ui.label(format!("{} / {}", rep.hbd, rep.hba));
                ui.end_row();
                ui.strong("RotBonds:"); ui.label(format!("{}", rep.rotatable_bonds));
                ui.strong("Fsp³:"); ui.label(format!("{:.2}", rep.fsp3));
                ui.end_row();
                ui.strong("Heavy Atoms:"); ui.label(format!("{}", rep.heavy_atoms));
                ui.strong("SA Score:"); ui.label(format!("{:.2} / 10", rep.sa_score));
                ui.end_row();
            });

            ui.separator();
            ui.horizontal(|ui| {
                ui.strong(format!("PAINS / Brenk Alerts ({})", rep.alerts.len()));
                ui.checkbox(&mut state.highlight_pains_on_canvas, "Highlight on Canvas");
            });
            if rep.alerts.is_empty() {
                ui.label(RichText::new("✓ Zero PAINS or reactive toxicophore alerts detected").color(Color32::from_rgb(20, 140, 70)));
            } else {
                egui::ScrollArea::vertical().max_height(130.0).show(ui, |ui| {
                    for al in &rep.alerts {
                        let col = match al.severity {
                            AlertSeverity::High => Color32::from_rgb(205, 40, 40),
                            AlertSeverity::Warning => Color32::from_rgb(200, 125, 15),
                        };
                        ui.label(RichText::new(format!("[{}] {}", al.category, al.name)).color(col).strong());
                    }
                });
            }
        });
    state.show_medchem = open;
}

// 2. pH Microspecies, Tautomers, Stereoisomers, Resonance & Electron Push Window
pub fn show_isomers_window(
    ctx: &egui::Context,
    state: &mut ProSuiteState,
    active_mol: Option<&ChemMolecule>,
) -> ProSuiteAction {
    if !state.show_isomers {
        return ProSuiteAction::None;
    }
    let mut act = ProSuiteAction::None;
    let mut open = state.show_isomers;

    egui::Window::new("pH Microspecies, Tautomers, Stereoisomers & Resonance")
        .open(&mut open)
        .default_pos([280.0, 95.0])
        .default_width(470.0)
        .resizable(true)
        .show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                for (lbl, tab) in [
                    ("pH / pI", IsomerTab::PhMicrospecies),
                    ("Tautomers", IsomerTab::Tautomers),
                    ("Stereoisomers", IsomerTab::Stereoisomers),
                    ("Resonance", IsomerTab::Resonance),
                    ("Electron Push", IsomerTab::ElectronPush),
                ] {
                    if ui.selectable_label(state.isomer_tab == tab, lbl).clicked() {
                        state.isomer_tab = tab;
                    }
                }
            });
            ui.separator();

            let Some(mol) = active_mol.filter(|m| !m.atoms.is_empty()) else {
                ui.weak("Draw or select a molecule on canvas first.");
                return;
            };

            match state.isomer_tab {
                IsomerTab::PhMicrospecies => {
                    ui.horizontal(|ui| {
                        ui.label("Target pH:");
                        ui.add(egui::Slider::new(&mut state.target_ph, 1.0..=14.0).fixed_decimals(1));
                        if ui.small_button("7.4").clicked() { state.target_ph = 7.4; }
                        if ui.small_button("2.0").clicked() { state.target_ph = 2.0; }
                        if ui.small_button("11.0").clicked() { state.target_ph = 11.0; }
                    });
                    let micro = compute_microspecies_at_ph(mol, state.target_ph);
                    ui.label(format!(
                        "Net Charge at pH {:.1}: {:+.2} (formal: {:+})  |  Isoelectric Point (pI): {}",
                        micro.ph,
                        micro.net_charge,
                        micro.formal_charge_sum,
                        micro.isoelectric_point.map(|v| format!("{v:.2}")).unwrap_or_else(|| "N/A".into())
                    ));
                    if micro.sites.is_empty() {
                        ui.weak("No ionizable acid/base functional groups detected.");
                    } else {
                        for s in &micro.sites {
                            let state_str = if s.is_acid {
                                if state.target_ph > s.pka { "Deprotonated (-)" } else { "Neutral (H)" }
                            } else if state.target_ph < s.pka {
                                "Protonated (+)"
                            } else {
                                "Neutral"
                            };
                            ui.label(format!("• Atom #{} ({}, pKa ~{:.1}): {state_str}", s.atom_idx + 1, s.label, s.pka));
                        }
                    }
                    ui.add_space(4.0);
                    ui.horizontal(|ui| {
                        if ui.button("Apply Protonation State to Active Molecule").clicked() {
                            act = ProSuiteAction::ReplaceActiveMolecule(micro.molecule.clone(), "Apply pH Microspecies");
                        }
                        if ui.button("Stamp Microspecies Copy").clicked() {
                            act = ProSuiteAction::StampMolecule(micro.molecule, "Stamp pH Microspecies");
                        }
                    });
                }
                IsomerTab::Tautomers => {
                    let tauts = enumerate_tautomers(mol, 12);
                    ui.label(format!("Found {} prototropic tautomer(s):", tauts.len()));
                    egui::ScrollArea::vertical().max_height(180.0).show(ui, |ui| {
                        for (i, (desc, tm)) in tauts.iter().enumerate() {
                            ui.horizontal(|ui| {
                                ui.label(format!("#{}: {desc}", i + 1));
                                if ui.small_button("Stamp").clicked() {
                                    act = ProSuiteAction::StampMolecule(tm.clone(), "Stamp Tautomer");
                                }
                            });
                        }
                    });
                    if tauts.len() > 1 && ui.button("Stamp All Tautomers as Grid").clicked() {
                        act = ProSuiteAction::StampMoleculeGrid(tauts, "Stamp Tautomer Grid");
                    }
                }
                IsomerTab::Stereoisomers => {
                    let isos = enumerate_stereoisomers(mol, 16);
                    ui.label(format!("Enumerated {} stereoisomer(s):", isos.len()));
                    egui::ScrollArea::vertical().max_height(180.0).show(ui, |ui| {
                        for (lbl, im) in &isos {
                            ui.horizontal(|ui| {
                                ui.label(lbl);
                                if ui.small_button("Stamp").clicked() {
                                    act = ProSuiteAction::StampMolecule(im.clone(), "Stamp Stereoisomer");
                                }
                            });
                        }
                    });
                    if isos.len() > 1 && ui.button("Stamp All Stereoisomers as Grid").clicked() {
                        act = ProSuiteAction::StampMoleculeGrid(isos, "Stamp Stereoisomer Grid");
                    }
                }
                IsomerTab::Resonance => {
                    let res = enumerate_resonance_structures(mol, 10);
                    ui.label(format!("Ranked {} resonance contributor(s):", res.len()));
                    egui::ScrollArea::vertical().max_height(180.0).show(ui, |ui| {
                        for (i, r) in res.iter().enumerate() {
                            ui.horizontal(|ui| {
                                ui.label(format!("#{} (score {}): {}", i + 1, r.weight_score, r.description));
                                if ui.small_button("Stamp").clicked() {
                                    act = ProSuiteAction::StampMolecule(r.molecule.clone(), "Stamp Resonance Contributor");
                                }
                            });
                        }
                    });
                    if res.len() > 1 && ui.button("Stamp All Contributors").clicked() {
                        let pairs = res.into_iter().map(|r| (r.description, r.molecule)).collect();
                        act = ProSuiteAction::StampMoleculeGrid(pairs, "Stamp Resonance Grid");
                    }
                }
                IsomerTab::ElectronPush => {
                    ui.label("Curved-Arrow Mechanistic Stepper (turn on View > Atom Numbers to see indices):");
                    let max_a = mol.atoms.len().max(1);
                    state.push_from_atom = state.push_from_atom.clamp(1, max_a);
                    state.push_to_atom = state.push_to_atom.clamp(1, max_a);
                    ui.horizontal(|ui| {
                        ui.label("From Atom #:");
                        ui.add(egui::DragValue::new(&mut state.push_from_atom).range(1..=max_a));
                        ui.label("-> To Atom #:");
                        ui.add(egui::DragValue::new(&mut state.push_to_atom).range(1..=max_a));
                        if ui.button("Push Electrons -> Stamp Product").clicked() {
                            match apply_electron_push_arrow(mol, state.push_from_atom - 1, state.push_to_atom - 1) {
                                Ok((next_mol, msg)) => {
                                    state.push_feedback = msg;
                                    act = ProSuiteAction::StampMolecule(next_mol, "Electron-Pushing Product");
                                }
                                Err(e) => state.push_feedback = e,
                            }
                        }
                    });
                    if !state.push_feedback.is_empty() {
                        ui.label(RichText::new(&state.push_feedback).small());
                    }
                }
            }
        });
    state.show_isomers = open;
    act
}

// 3. HRMS Isotope Envelope & Adduct Simulator Window
pub fn show_hrms_window(
    ctx: &egui::Context,
    state: &mut ProSuiteState,
    active_mol: Option<&ChemMolecule>,
) {
    if !state.show_hrms {
        return;
    }
    let mut open = state.show_hrms;
    egui::Window::new("HRMS Isotope Envelope & Adduct Simulator")
        .open(&mut open)
        .default_pos([340.0, 120.0])
        .default_width(540.0)
        .resizable(true)
        .show(ctx, |ui| {
            let Some(mol) = active_mol.filter(|m| !m.atoms.is_empty()) else {
                ui.weak("Draw or select a molecule to simulate its HRMS isotopic pattern.");
                return;
            };
            let sim = simulate_hrms(mol, &state.hrms_adduct);
            ui.horizontal(|ui| {
                ui.strong(format!("Formula: {}  |  Neutral Exact Mass: {:.4}", sim.formula, sim.neutral_exact_mass));
                ui.separator();
                egui::ComboBox::from_id_salt("hrms_adduct_sel")
                    .selected_text(&state.hrms_adduct)
                    .show_ui(ui, |ui| {
                        for ad in &sim.adducts {
                            ui.selectable_value(&mut state.hrms_adduct, ad.name.to_string(), format!("{} ({:.4})", ad.name, ad.mz));
                        }
                    });
            });
            ui.horizontal(|ui| {
                ui.label("Resolution R:");
                ui.add(egui::Slider::new(&mut state.hrms_resolution, 5_000.0..=100_000.0).fixed_decimals(0));
                ui.label("Observed m/z:");
                ui.add(egui::TextEdit::singleline(&mut state.hrms_observed_mz).desired_width(82.0).hint_text("e.g. 181.0501"));
                if let Ok(obs) = state.hrms_observed_mz.trim().parse::<f64>() {
                    if let Some(err) = sim.ppm_error(obs) {
                        let col = if err.abs() <= 5.0 { Color32::from_rgb(20, 145, 70) } else { Color32::from_rgb(205, 45, 45) };
                        ui.label(RichText::new(format!("Δ = {err:+.2} ppm")).strong().color(col));
                    }
                }
            });

            // Simulated HRMS profile + centroid stick plot
            let (rect, _) = ui.allocate_exact_size(Vec2::new(515.0, 150.0), Sense::hover());
            let p = ui.painter();
            p.rect_filled(rect, 2.0, Color32::WHITE);
            p.rect_stroke(rect, 2.0, Stroke::new(1.0_f32, Color32::from_gray(145)), egui::StrokeKind::Inside);
            let plot = rect.shrink2(Vec2::new(18.0, 18.0));
            let base_y = plot.max.y - 10.0;
            p.line_segment([Pos2::new(plot.min.x, base_y), Pos2::new(plot.max.x, base_y)], Stroke::new(1.0_f32, Color32::BLACK));

            let prof = sim.simulate_profile(state.hrms_resolution, 280);
            if let (Some(first), Some(last)) = (prof.first(), prof.last()) {
                let min_x = first.0;
                let max_x = last.0.max(min_x + 0.5);
                let to_x = |mz: f64| plot.min.x + (((mz - min_x) / (max_x - min_x)).clamp(0.0, 1.0) as f32) * plot.width();
                let h_scale = (base_y - plot.min.y - 8.0).max(20.0);
                // Profile curve
                let pts: Vec<Pos2> = prof.iter().map(|&(mz, y)| Pos2::new(to_x(mz), base_y - (y as f32 / 100.0) * h_scale)).collect();
                for w in pts.windows(2) {
                    p.line_segment([w[0], w[1]], Stroke::new(1.2_f32, Color32::from_rgb(25, 95, 195)));
                }
                // Centroid sticks + labels for top peaks
                for pk in sim.peaks.iter().take(10) {
                    let x = to_x(pk.mz);
                    let y = base_y - (pk.rel_intensity as f32 / 100.0) * h_scale;
                    p.line_segment([Pos2::new(x, base_y), Pos2::new(x, y)], Stroke::new(1.4_f32, Color32::from_rgb(190, 35, 35)));
                    if pk.rel_intensity >= 10.0 {
                        p.text(Pos2::new(x, y - 3.0), egui::Align2::CENTER_BOTTOM, format!("{:.4}", pk.mz), FontId::proportional(8.5), Color32::BLACK);
                    }
                }
            }
            ui.add_space(4.0);
            egui::ScrollArea::vertical().max_height(115.0).show(ui, |ui| {
                egui::Grid::new("hrms_adduct_tbl").striped(true).num_columns(4).spacing([14.0, 2.0]).show(ui, |ui| {
                    ui.strong("Adduct"); ui.strong("m/z"); ui.strong("Adduct"); ui.strong("m/z");
                    ui.end_row();
                    for pair in sim.adducts.chunks(2) {
                        ui.label(pair[0].name);
                        ui.label(format!("{:.4}", pair[0].mz));
                        if pair.len() > 1 {
                            ui.label(pair[1].name);
                            ui.label(format!("{:.4}", pair[1].mz));
                        } else {
                            ui.label(""); ui.label("");
                        }
                        ui.end_row();
                    }
                });
            });
        });
    state.show_hrms = open;
}

// 4. Markush R-Group Library Enumerator & SAR Grid Window
pub fn show_markush_window(
    ctx: &egui::Context,
    state: &mut ProSuiteState,
    active_mol: Option<&ChemMolecule>,
    bond_len: f64,
) -> ProSuiteAction {
    if !state.show_markush {
        return ProSuiteAction::None;
    }
    let mut act = ProSuiteAction::None;
    let mut open = state.show_markush;

    egui::Window::new("Markush R-Group Library Enumerator & SAR Grid")
        .open(&mut open)
        .default_pos([200.0, 100.0])
        .default_width(620.0)
        .resizable(true)
        .show(ctx, |ui| {
            ui.small("Label atoms on your scaffold as R1 (or R) and R2 — or enumerate directly on any drawn core:");
            ui.horizontal(|ui| {
                ui.label("R1 list:");
                ui.add(egui::TextEdit::singleline(&mut state.r1_input).desired_width(210.0));
                ui.label("R2 list:");
                ui.add(egui::TextEdit::singleline(&mut state.r2_input).desired_width(140.0));
                if ui.button("Enumerate Library").clicked() {
                    if let Some(m) = active_mol.filter(|m| !m.atoms.is_empty()) {
                        state.sar_analogs = enumerate_rgroup_library(m, &state.r1_input, &state.r2_input, bond_len);
                    } else {
                        act = ProSuiteAction::StatusMsg("Draw a scaffold molecule first to enumerate R-groups".into());
                    }
                }
            });
            ui.horizontal(|ui| {
                ui.label("Sort by:");
                for (lbl, k) in [
                    ("ID", SarSortKey::Id),
                    ("MW", SarSortKey::Mw),
                    ("LogP", SarSortKey::LogP),
                    ("CNS MPO ↓", SarSortKey::CnsMpo),
                    ("QED ↓", SarSortKey::Qed),
                    ("Alerts", SarSortKey::Alerts),
                ] {
                    if ui.selectable_label(state.sar_sort == k, lbl).clicked() {
                        state.sar_sort = k;
                        match k {
                            SarSortKey::Id => state.sar_analogs.sort_by(|a, b| a.id.cmp(&b.id)),
                            SarSortKey::Mw => state.sar_analogs.sort_by(|a, b| a.medchem.mw.partial_cmp(&b.medchem.mw).unwrap_or(std::cmp::Ordering::Equal)),
                            SarSortKey::LogP => state.sar_analogs.sort_by(|a, b| a.medchem.log_p.partial_cmp(&b.medchem.log_p).unwrap_or(std::cmp::Ordering::Equal)),
                            SarSortKey::CnsMpo => state.sar_analogs.sort_by(|a, b| b.medchem.cns_mpo.partial_cmp(&a.medchem.cns_mpo).unwrap_or(std::cmp::Ordering::Equal)),
                            SarSortKey::Qed => state.sar_analogs.sort_by(|a, b| b.medchem.qed.partial_cmp(&a.medchem.qed).unwrap_or(std::cmp::Ordering::Equal)),
                            SarSortKey::Alerts => state.sar_analogs.sort_by(|a, b| a.medchem.alerts.len().cmp(&b.medchem.alerts.len())),
                        }
                    }
                }
            });
            ui.separator();
            egui::ScrollArea::vertical().max_height(220.0).show(ui, |ui| {
                egui::Grid::new("sar_grid").striped(true).num_columns(9).spacing([8.0, 3.0]).show(ui, |ui| {
                    for h in ["ID", "R1", "R2", "Formula", "MW", "LogP", "MPO", "QED", "Action"] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for a in &state.sar_analogs {
                        ui.label(&a.id);
                        ui.label(&a.r1_label);
                        ui.label(&a.r2_label);
                        ui.label(&a.formula);
                        ui.label(format!("{:.1}", a.medchem.mw));
                        ui.label(format!("{:.2}", a.medchem.log_p));
                        ui.label(format!("{:.2}", a.medchem.cns_mpo));
                        ui.label(format!("{:.2}", a.medchem.qed));
                        if ui.small_button("Stamp").clicked() {
                            act = ProSuiteAction::StampMolecule(a.molecule.clone(), "Stamp Analog");
                        }
                        ui.end_row();
                    }
                });
            });
            ui.separator();
            ui.horizontal(|ui| {
                ui.label(format!("{} analog(s)", state.sar_analogs.len()));
                if ui.add_enabled(!state.sar_analogs.is_empty(), egui::Button::new("Export SAR CSV")).clicked() {
                    act = ProSuiteAction::ExportTextFile {
                        suggested_name: "chiralcraft-sar-library.csv",
                        content: export_sar_csv(&state.sar_analogs),
                    };
                }
                if ui.add_enabled(!state.sar_analogs.is_empty(), egui::Button::new("Export Library SDF")).clicked() {
                    act = ProSuiteAction::ExportTextFile {
                        suggested_name: "chiralcraft-sar-library.sdf",
                        content: export_sar_sdf(&state.sar_analogs),
                    };
                }
            });
        });
    state.show_markush = open;
    act
}

// 5. Offline Retrosynthetic Disconnection Advisor Window
pub fn show_retrosynthesis_window(
    ctx: &egui::Context,
    state: &mut ProSuiteState,
    active_mol: Option<&ChemMolecule>,
    bond_len: f64,
) -> ProSuiteAction {
    if !state.show_retro {
        state.highlighted_retro_bond = None;
        return ProSuiteAction::None;
    }
    let mut act = ProSuiteAction::None;
    let mut open = state.show_retro;

    egui::Window::new("Retrosynthetic Disconnection Advisor")
        .open(&mut open)
        .default_pos([240.0, 140.0])
        .default_width(530.0)
        .resizable(true)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                if ui.button("Analyze Target Structure").clicked() {
                    if let Some(m) = active_mol.filter(|m| !m.atoms.is_empty()) {
                        state.retro_hits = suggest_disconnections(m, bond_len);
                        state.highlighted_retro_bond = state.retro_hits.first().map(|h| h.bond_idx);
                    }
                }
                ui.weak("Click a disconnection to highlight the strategic bond on canvas.");
            });
            ui.separator();
            if state.retro_hits.is_empty() {
                ui.weak("No acyclic strategic disconnections analyzed yet (click Analyze Target Structure).");
            } else {
                egui::ScrollArea::vertical().max_height(230.0).show(ui, |ui| {
                    for hit in &state.retro_hits {
                        let is_hl = state.highlighted_retro_bond == Some(hit.bond_idx);
                        ui.group(|ui| {
                            ui.horizontal(|ui| {
                                if ui.selectable_label(is_hl, RichText::new(hit.reaction_name).strong()).clicked() {
                                    state.highlighted_retro_bond = Some(hit.bond_idx);
                                }
                                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                    if ui.small_button("Stamp Retrosynthesis (=>)").clicked() {
                                        act = ProSuiteAction::StampRetrosynthesis(hit.clone());
                                    }
                                    ui.label(RichText::new(format!("Score {}", hit.score)).small().weak());
                                });
                            });
                            ui.label(RichText::new(format!("Synthons: {}  |  Forward: {}", hit.synthons_summary, hit.reagents)).small());
                        });
                    }
                });
            }
        });
    state.show_retro = open;
    act
}

// 6. Publication TOC Graphic & Multi-Step Scheme Composer Window
pub fn show_publication_composer_window(
    ctx: &egui::Context,
    state: &mut ProSuiteState,
    scene: &DocumentScene,
) -> ProSuiteAction {
    if !state.show_pub_composer {
        return ProSuiteAction::None;
    }
    let mut act = ProSuiteAction::None;
    let mut open = state.show_pub_composer;
    let presets = TocPreset::presets();
    state.selected_toc_idx = state.selected_toc_idx.min(presets.len() - 1);

    egui::Window::new("Publication TOC & Scheme Composer")
        .open(&mut open)
        .default_pos([780.0, 260.0])
        .default_width(390.0)
        .resizable(true)
        .show(ctx, |ui| {
            ui.strong("Journal TOC Graphic Frame");
            egui::ComboBox::from_id_salt("toc_preset_sel")
                .selected_text(presets[state.selected_toc_idx].name)
                .show_ui(ui, |ui| {
                    for (i, pr) in presets.iter().enumerate() {
                        ui.selectable_value(&mut state.selected_toc_idx, i, pr.name);
                    }
                });
            if ui.button("Insert TOC Boundary Frame").clicked() {
                act = ProSuiteAction::InsertTocFrame(presets[state.selected_toc_idx]);
            }
            if let Some(bbox) = state.last_toc_bbox {
                let over = count_molecules_overflowing_frame(scene, bbox);
                if over == 0 {
                    ui.label(RichText::new("✓ All molecules fit cleanly inside TOC frame").color(Color32::from_rgb(20, 145, 70)));
                } else {
                    ui.label(RichText::new(format!("⚠ {over} molecule(s) extend outside TOC boundary")).color(Color32::from_rgb(205, 115, 15)));
                }
            }

            ui.separator();
            ui.strong("Smart Reaction Arrow Conditions");
            ui.horizontal(|ui| {
                ui.label("Above:");
                ui.text_edit_singleline(&mut state.cond_above);
            });
            ui.horizontal(|ui| {
                ui.label("Below:");
                ui.text_edit_singleline(&mut state.cond_below);
            });
            if ui.button("Stamp Conditions on Active Arrow").clicked() {
                act = ProSuiteAction::StampArrowConditions {
                    above: state.cond_above.clone(),
                    below: state.cond_below.clone(),
                };
            }

            ui.separator();
            ui.strong("Serpentine Multi-Step Scheme Layout");
            ui.horizontal(|ui| {
                ui.label("Steps per row:");
                ui.add(egui::DragValue::new(&mut state.serpentine_per_row).range(2..=6));
                if ui.button("Auto-Wrap Scheme").clicked() {
                    act = ProSuiteAction::WrapSerpentine(state.serpentine_per_row);
                }
            });
        });
    state.show_pub_composer = open;
    act
}
