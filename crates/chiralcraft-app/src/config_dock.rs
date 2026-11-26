//! Persistent User Configuration (`~/.config/chiralcraft/config.json`),
//! Unified Right-Side Tabbed Studio Dock, and Computational Chemistry (ORCA / Gaussian / xTB) Studio UI.

use std::fs;
use std::path::PathBuf;

use chiralcraft_chem::{
    analyze_organometallic_complexes, compute_analysis_report, compute_medchem_report,
    estimate_chemical_properties, predict_ir_spectrum, predict_uv_vis, simulate_hrms,
    AlertSeverity, ChemMolecule,
};
use chiralcraft_doc::{CopyBackground, ExportOptions, SceneTable, SceneTlcPlate};
use chiralcraft_formats::{
    export_quantum_input, parse_quantum_output, QuantumJobSpec, QuantumJobType, QuantumOutputReport,
    QuantumProgram,
};
use egui::{Color32, RichText};

use crate::chem_db::JsonValue;

// ============================================================================
// 1. Persistent User Configuration (~/.config/chiralcraft/config.json)
// ============================================================================

#[derive(Debug, Clone)]
pub struct UserConfig {
    pub theme: String, // "light", "dark", "slide"
    pub png_dpi: u32,
    pub svg_padding_pt: f64,
    pub copy_bg: String, // "transparent", "white", "slide"
    pub rear_opacity: f32,
    pub crossing_clearance_pt: f64,
    pub show_right_dock: bool,
    pub recent_files: Vec<String>,
}

impl Default for UserConfig {
    fn default() -> Self {
        Self {
            theme: "light".to_string(),
            png_dpi: 300,
            svg_padding_pt: 20.0,
            copy_bg: "transparent".to_string(),
            rear_opacity: 0.45,
            crossing_clearance_pt: 3.6,
            show_right_dock: true,
            recent_files: Vec::new(),
        }
    }
}

fn config_file_path() -> Option<PathBuf> {
    if let Ok(xdg) = std::env::var("XDG_CONFIG_HOME") {
        if !xdg.trim().is_empty() {
            return Some(PathBuf::from(xdg).join("chiralcraft").join("config.json"));
        }
    }
    std::env::var("HOME")
        .ok()
        .map(|h| PathBuf::from(h).join(".config").join("chiralcraft").join("config.json"))
}

fn json_esc(s: &str) -> String {
    let mut o = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            c => o.push(c),
        }
    }
    o
}

impl UserConfig {
    pub fn load() -> Self {
        let mut cfg = Self::default();
        let Some(path) = config_file_path() else {
            return cfg;
        };
        let Ok(txt) = fs::read_to_string(&path) else {
            return cfg;
        };
        let Ok(root) = JsonValue::parse(&txt) else {
            return cfg;
        };
        if let Some(s) = root.get("theme").and_then(|v| v.as_str()) {
            cfg.theme = s.to_string();
        }
        if let Some(JsonValue::Number(n)) = root.get("png_dpi") {
            cfg.png_dpi = (*n as u32).clamp(72, 1200);
        }
        if let Some(JsonValue::Number(n)) = root.get("svg_padding_pt") {
            cfg.svg_padding_pt = n.clamp(4.0, 80.0);
        }
        if let Some(s) = root.get("copy_bg").and_then(|v| v.as_str()) {
            cfg.copy_bg = s.to_string();
        }
        if let Some(JsonValue::Number(n)) = root.get("rear_opacity") {
            cfg.rear_opacity = (*n as f32).clamp(0.15, 1.0);
        }
        if let Some(JsonValue::Number(n)) = root.get("crossing_clearance_pt") {
            cfg.crossing_clearance_pt = n.clamp(0.0, 8.0);
        }
        if let Some(JsonValue::Bool(b)) = root.get("show_right_dock") {
            cfg.show_right_dock = *b;
        }
        if let Some(JsonValue::Array(arr)) = root.get("recent_files") {
            cfg.recent_files = arr
                .iter()
                .filter_map(|v| v.as_str().map(|s| s.to_string()))
                .take(10)
                .collect();
        }
        cfg
    }

    pub fn save(&self) {
        let Some(path) = config_file_path() else {
            return;
        };
        if let Some(parent) = path.parent() {
            let _ = fs::create_dir_all(parent);
        }
        let recents: Vec<String> = self
            .recent_files
            .iter()
            .take(10)
            .map(|r| format!("\"{}\"", json_esc(r)))
            .collect();
        let json = format!(
            "{{\n  \"theme\": \"{}\",\n  \"png_dpi\": {},\n  \"svg_padding_pt\": {:.1},\n  \"copy_bg\": \"{}\",\n  \"rear_opacity\": {:.2},\n  \"crossing_clearance_pt\": {:.1},\n  \"show_right_dock\": {},\n  \"recent_files\": [{}]\n}}\n",
            json_esc(&self.theme),
            self.png_dpi,
            self.svg_padding_pt,
            json_esc(&self.copy_bg),
            self.rear_opacity,
            self.crossing_clearance_pt,
            self.show_right_dock,
            recents.join(", ")
        );
        let _ = fs::write(path, json);
    }

    pub fn push_recent(&mut self, path: &str) {
        if path.trim().is_empty() {
            return;
        }
        self.recent_files.retain(|p| p != path);
        self.recent_files.insert(0, path.to_string());
        self.recent_files.truncate(10);
        self.save();
    }

    pub fn to_export_options(&self) -> ExportOptions {
        let copy_bg = match self.copy_bg.as_str() {
            "white" => CopyBackground::White,
            "slide" => CopyBackground::DarkSlide,
            _ => CopyBackground::Transparent,
        };
        ExportOptions {
            png_dpi: self.png_dpi,
            svg_padding_pt: self.svg_padding_pt,
            copy_bg,
        }
    }
}

// ============================================================================
// 2. Quantum Chemistry Studio State (ORCA / Gaussian / xTB)
// ============================================================================

pub struct QuantumStudioState {
    pub show_window: bool,
    pub spec: QuantumJobSpec,
    pub import_path: String,
    pub parsed_report: Option<QuantumOutputReport>,
    pub parse_error: String,
}

impl Default for QuantumStudioState {
    fn default() -> Self {
        Self {
            show_window: false,
            spec: QuantumJobSpec::default(),
            import_path: String::new(),
            parsed_report: None,
            parse_error: String::new(),
        }
    }
}

pub enum DockQuantumAction {
    None,
    CopyText(String, String),
    ExportFile { path: String, content: String },
    StampOptimized3d(ChemMolecule),
    OpenFloatingStudio(&'static str),
    ConfigChanged,
}

pub fn render_quantum_panel_contents(
    ui: &mut egui::Ui,
    qstate: &mut QuantumStudioState,
    active_mol: Option<&ChemMolecule>,
) -> DockQuantumAction {
    let mut act = DockQuantumAction::None;
    ui.strong("Export DFT / Semi-Empirical Input");
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut qstate.spec.program, QuantumProgram::Orca, "ORCA (.inp)");
        ui.selectable_value(&mut qstate.spec.program, QuantumProgram::Gaussian, "Gaussian (.gjf)");
        ui.selectable_value(&mut qstate.spec.program, QuantumProgram::Xtb, "xTB (GFN2)");
    });
    egui::Grid::new("q_grid").num_columns(2).spacing([6.0, 3.0]).show(ui, |ui| {
        ui.label("Method:");
        ui.add(egui::TextEdit::singleline(&mut qstate.spec.method).desired_width(110.0));
        ui.end_row();
        ui.label("Basis:");
        ui.add(egui::TextEdit::singleline(&mut qstate.spec.basis).desired_width(110.0));
        ui.end_row();
        ui.label("Job:");
        egui::ComboBox::from_id_salt("q_job_type")
            .selected_text(format!("{:?}", qstate.spec.job_type))
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut qstate.spec.job_type, QuantumJobType::OptFreq, "Opt + Freq");
                ui.selectable_value(&mut qstate.spec.job_type, QuantumJobType::SinglePoint, "Single Point");
                ui.selectable_value(&mut qstate.spec.job_type, QuantumJobType::TdDft, "TD-DFT (UV-Vis)");
                ui.selectable_value(&mut qstate.spec.job_type, QuantumJobType::Nmr, "GIAO NMR");
            });
        ui.end_row();
        ui.label("Solvent:");
        egui::ComboBox::from_id_salt("q_solv")
            .selected_text(&qstate.spec.solvent)
            .show_ui(ui, |ui| {
                for s in ["None", "THF", "DCM", "Water", "DMSO", "MeCN", "Toluene"] {
                    ui.selectable_value(&mut qstate.spec.solvent, s.to_string(), s);
                }
            });
        ui.end_row();
        ui.label("Charge / Mult:");
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut qstate.spec.charge).range(-6..=6));
            ui.add(egui::DragValue::new(&mut qstate.spec.multiplicity).range(1..=9));
        });
        ui.end_row();
    });

    let mut preview = active_mol
        .filter(|m| !m.atoms.is_empty())
        .map(|m| export_quantum_input(m, &qstate.spec))
        .unwrap_or_else(|| "# Draw a molecule to generate an ORCA / Gaussian / xTB input deck".into());
    egui::ScrollArea::vertical().id_salt("q_prev_scroll").max_height(125.0).show(ui, |ui| {
        ui.add(egui::TextEdit::multiline(&mut preview).font(egui::TextStyle::Monospace).desired_width(f32::INFINITY).desired_rows(6));
    });
    ui.horizontal(|ui| {
        if ui.small_button("Copy Deck").clicked() {
            act = DockQuantumAction::CopyText(preview.clone(), "Copied quantum input deck".into());
        }
        let ext = match qstate.spec.program {
            QuantumProgram::Orca => "inp",
            QuantumProgram::Gaussian => "gjf",
            QuantumProgram::Xtb => "xyz",
        };
        if ui.small_button(format!("Save /tmp/chiralcraft.{ext}")).clicked() {
            act = DockQuantumAction::ExportFile {
                path: format!("/tmp/chiralcraft.{ext}"),
                content: preview,
            };
        }
    });

    ui.separator();
    ui.strong("Import ORCA (.out) / Gaussian (.log)");
    ui.horizontal(|ui| {
        ui.add(egui::TextEdit::singleline(&mut qstate.import_path).desired_width(145.0).hint_text("path/to/job.out"));
        if ui.small_button("Parse").clicked() && !qstate.import_path.trim().is_empty() {
            match fs::read_to_string(qstate.import_path.trim())
                .map_err(|e| e.to_string())
                .and_then(|s| parse_quantum_output(&s))
            {
                Ok(rep) => {
                    qstate.parsed_report = Some(rep);
                    qstate.parse_error.clear();
                }
                Err(e) => {
                    qstate.parse_error = e;
                }
            }
        }
    });
    if !qstate.parse_error.is_empty() {
        ui.label(RichText::new(&qstate.parse_error).small().color(Color32::from_rgb(205, 45, 45)));
    }
    if let Some(rep) = &qstate.parsed_report {
        ui.group(|ui| {
            ui.label(RichText::new(format!("{} — {}", rep.program, if rep.converged { "Converged ✓" } else { "Incomplete" })).strong());
            if let Some(e) = rep.energy_hartree {
                ui.small(format!("SCF Energy: {e:.6} Eh"));
            }
            if let (Some(h), Some(l), Some(g)) = (rep.homo_ev, rep.lumo_ev, rep.homo_lumo_gap_ev) {
                ui.small(format!("HOMO: {h:.2} eV | LUMO: {l:.2} eV | Gap: {g:.2} eV"));
            }
            if let Some(d) = rep.dipole_debye {
                ui.small(format!("Dipole: {d:.3} Debye"));
            }
            if !rep.ir_frequencies_cm1.is_empty() {
                ui.small(format!("DFT IR Modes: {} (max {:.0} cm⁻¹)", rep.ir_frequencies_cm1.len(), rep.ir_frequencies_cm1.iter().copied().fold(0.0_f64, f64::max)));
            }
            if let Some(opt) = &rep.optimized_molecule {
                if ui.small_button("Stamp DFT 3D Geometry on Canvas").clicked() {
                    act = DockQuantumAction::StampOptimized3d(opt.clone());
                }
            }
        });
    }
    act
}

pub fn show_quantum_window(
    ctx: &egui::Context,
    qstate: &mut QuantumStudioState,
    active_mol: Option<&ChemMolecule>,
) -> DockQuantumAction {
    if !qstate.show_window {
        return DockQuantumAction::None;
    }
    let mut act = DockQuantumAction::None;
    let mut open = qstate.show_window;
    egui::Window::new("Computational Chemistry Bridge (ORCA / Gaussian / xTB)")
        .open(&mut open)
        .default_pos([310.0, 100.0])
        .default_width(430.0)
        .resizable(true)
        .show(ctx, |ui| {
            act = render_quantum_panel_contents(ui, qstate, active_mol);
        });
    qstate.show_window = open;
    act
}

// ============================================================================
// 3. Unified Right-Side Tabbed Studio Dock
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DockTab {
    Analysis,
    MedChem,
    Spectra,
    Quantum,
    Inspector,
}

pub fn show_right_studio_dock(
    ctx: &egui::Context,
    open: &mut bool,
    tab: &mut DockTab,
    active_mol: Option<&ChemMolecule>,
    live_name: Option<&str>,
    qstate: &mut QuantumStudioState,
    tables: &mut [SceneTable],
    tlc_plates: &mut [SceneTlcPlate],
    rear_opacity: &mut f32,
    crossing_clearance_pt: &mut f64,
    export_opts: &mut ExportOptions,
) -> DockQuantumAction {
    if !*open {
        return DockQuantumAction::None;
    }
    let mut act = DockQuantumAction::None;
    egui::SidePanel::right("chiralcraft_studio_dock")
        .default_width(295.0)
        .min_width(240.0)
        .max_width(440.0)
        .resizable(true)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.strong("ChiralCraft Studio");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.small_button("×").on_hover_text("Hide Studio Dock").clicked() {
                        *open = false;
                    }
                });
            });
            ui.horizontal_wrapped(|ui| {
                for (lbl, t) in [
                    ("Props", DockTab::Analysis),
                    ("MedChem", DockTab::MedChem),
                    ("Spectra", DockTab::Spectra),
                    ("DFT", DockTab::Quantum),
                    ("Inspect", DockTab::Inspector),
                ] {
                    ui.selectable_value(tab, t, lbl);
                }
            });
            ui.separator();

            egui::ScrollArea::vertical().show(ui, |ui| match *tab {
                DockTab::Analysis => {
                    let Some(mol) = active_mol.filter(|m| !m.atoms.is_empty()) else {
                        ui.weak("Draw or import a structure to see live chemistry analysis.");
                        return;
                    };
                    let rep = compute_analysis_report(mol);
                    let props = estimate_chemical_properties(mol);
                    if let Some(nm) = live_name {
                        ui.label(RichText::new(nm).strong());
                    }
                    ui.label(format!("Formula: {}", rep.formula));
                    ui.label(format!("Exact Mass: {:.4}", rep.exact_mass));
                    ui.label(format!("Mol. Wt.: {:.2}", rep.mol_weight));
                    ui.small(format!("Elem. Anal.: {}", rep.elem_analysis));
                    ui.separator();
                    ui.label(format!("LogP: {:.2}  |  tPSA: {:.1} Å²", props.log_p, props.tpsa));
                    ui.label(format!("LogS: {:.2}  |  MR: {:.1}", props.log_s, props.molar_refractivity));
                    if let Some(pka) = props.pka {
                        ui.label(format!("Est. pKa: {pka:.2}"));
                    }
                    let metals = analyze_organometallic_complexes(mol);
                    if !metals.is_empty() {
                        ui.separator();
                        ui.strong("Organometallic 18e⁻");
                        for r in metals {
                            ui.small(format!("{}: {}e⁻ (ox {:+}, d^{}) — {}", r.metal_symbol, r.total_valence_electrons, r.oxidation_state, r.d_electrons, r.predicted_geometry));
                        }
                    }
                    ui.separator();
                    ui.horizontal_wrapped(|ui| {
                        if ui.small_button("3D Viewer").clicked() { act = DockQuantumAction::OpenFloatingStudio("3d"); }
                        if ui.small_button("pH / Isomers").clicked() { act = DockQuantumAction::OpenFloatingStudio("isomers"); }
                        if ui.small_button("Retrosynthesis").clicked() { act = DockQuantumAction::OpenFloatingStudio("retro"); }
                        if ui.small_button("LaTeX").clicked() { act = DockQuantumAction::OpenFloatingStudio("latex"); }
                    });
                }
                DockTab::MedChem => {
                    let Some(mol) = active_mol.filter(|m| !m.atoms.is_empty()) else {
                        ui.weak("Draw a structure for live Lipinski, CNS MPO, QED, and PAINS alerts.");
                        return;
                    };
                    let mc = compute_medchem_report(mol);
                    let ok_col = |ok: bool| if ok { Color32::from_rgb(20, 145, 70) } else { Color32::from_rgb(205, 45, 45) };
                    ui.label(RichText::new(format!("Lipinski: {} ({} viol)", if mc.lipinski_pass { "PASS" } else { "FAIL" }, mc.lipinski_violations)).strong().color(ok_col(mc.lipinski_pass)));
                    ui.label(format!("Veber: {}  |  Ghose: {}", if mc.veber_pass { "PASS" } else { "FAIL" }, if mc.ghose_pass { "PASS" } else { "FAIL" }));
                    ui.label(format!("CNS MPO: {:.2} / 6.0", mc.cns_mpo));
                    ui.add(egui::ProgressBar::new((mc.cns_mpo / 6.0) as f32));
                    ui.label(format!("QED: {:.2}  |  Fsp³: {:.2}  |  SA: {:.1}", mc.qed, mc.fsp3, mc.sa_score));
                    ui.separator();
                    ui.strong(format!("PAINS / Brenk Alerts ({})", mc.alerts.len()));
                    if mc.alerts.is_empty() {
                        ui.label(RichText::new("✓ Clean (0 alerts)").color(Color32::from_rgb(20, 145, 70)));
                    } else {
                        for al in &mc.alerts {
                            let c = if al.severity == AlertSeverity::High { Color32::from_rgb(205, 40, 40) } else { Color32::from_rgb(200, 125, 15) };
                            ui.label(RichText::new(format!("[{}] {}", al.category, al.name)).small().color(c));
                        }
                    }
                    ui.separator();
                    ui.horizontal_wrapped(|ui| {
                        if ui.small_button("Full Scorecard").clicked() { act = DockQuantumAction::OpenFloatingStudio("medchem"); }
                        if ui.small_button("Markush SAR").clicked() { act = DockQuantumAction::OpenFloatingStudio("markush"); }
                        if ui.small_button("Screen SDF").clicked() { act = DockQuantumAction::OpenFloatingStudio("screen"); }
                    });
                }
                DockTab::Spectra => {
                    let Some(mol) = active_mol.filter(|m| !m.atoms.is_empty()) else {
                        ui.weak("Draw a structure for live FT-IR, UV–Vis, HRMS, and NMR previews.");
                        return;
                    };
                    let ir = predict_ir_spectrum(mol);
                    let uv = predict_uv_vis(mol);
                    let hrms = simulate_hrms(mol, "[M+H]+");
                    ui.strong("HRMS [M+H]⁺");
                    if let Some(pk) = hrms.peaks.first() {
                        ui.label(format!("Base peak: {:.4} ({})", pk.mz, hrms.formula));
                    }
                    ui.separator();
                    ui.strong("UV–Vis");
                    ui.label(format!("λmax ≈ {:.0} nm — {}", uv.lambda_max_nm, uv.perceived_color_name));
                    ui.separator();
                    ui.strong("Diagnostic FT-IR");
                    for b in ir.bands.iter().take(5) {
                        ui.small(format!("{:.0} cm⁻¹ ({}) {}", b.wavenumber_cm1, b.strength, b.assignment));
                    }
                    ui.separator();
                    ui.horizontal_wrapped(|ui| {
                        if ui.small_button("IR / UV / 2D NMR").clicked() { act = DockQuantumAction::OpenFloatingStudio("spec"); }
                        if ui.small_button("HRMS Plot").clicked() { act = DockQuantumAction::OpenFloatingStudio("hrms"); }
                        if ui.small_button("1H NMR").clicked() { act = DockQuantumAction::OpenFloatingStudio("nmr1h"); }
                        if ui.small_button("13C NMR").clicked() { act = DockQuantumAction::OpenFloatingStudio("nmr13c"); }
                    });
                }
                DockTab::Quantum => {
                    act = render_quantum_panel_contents(ui, qstate, active_mol);
                }
                DockTab::Inspector => {
                    let mut cfg_changed = false;
                    ui.strong("3D & Crossing Rendering");
                    ui.horizontal(|ui| {
                        ui.label("Rear Opacity:");
                        cfg_changed |= ui.add(egui::Slider::new(rear_opacity, 0.15..=1.0).fixed_decimals(2)).changed();
                    });
                    ui.horizontal(|ui| {
                        ui.label("Bridge Gap:");
                        cfg_changed |= ui.add(egui::Slider::new(crossing_clearance_pt, 0.0..=8.0).fixed_decimals(1)).changed();
                    });
                    ui.separator();
                    ui.strong("Export & Copy");
                    ui.horizontal(|ui| {
                        ui.label("DPI:");
                        for dpi in [150u32, 300, 600] {
                            if ui.selectable_label(export_opts.png_dpi == dpi, format!("{dpi}")).clicked() {
                                export_opts.png_dpi = dpi;
                                cfg_changed = true;
                            }
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("Copy BG:");
                        for (lbl, bg) in [
                            ("Transp.", CopyBackground::Transparent),
                            ("White", CopyBackground::White),
                            ("Slide", CopyBackground::DarkSlide),
                        ] {
                            if ui.selectable_label(export_opts.copy_bg == bg, lbl).clicked() {
                                export_opts.copy_bg = bg;
                                cfg_changed = true;
                            }
                        }
                    });
                    if cfg_changed {
                        act = DockQuantumAction::ConfigChanged;
                    }
                    ui.separator();
                    ui.small(format!("Canvas Objects: {} table(s), {} TLC plate(s)", tables.len(), tlc_plates.len()));
                    ui.horizontal_wrapped(|ui| {
                        if ui.small_button("TLC / Table Editor").clicked() { act = DockQuantumAction::OpenFloatingStudio("inspector"); }
                        if ui.small_button("Energy Profile").clicked() { act = DockQuantumAction::OpenFloatingStudio("energy"); }
                        if ui.small_button("TOC Composer").clicked() { act = DockQuantumAction::OpenFloatingStudio("toc"); }
                        if ui.small_button("ELN & GHS").clicked() { act = DockQuantumAction::OpenFloatingStudio("eln"); }
                    });
                }
            });
        });
    act
}
