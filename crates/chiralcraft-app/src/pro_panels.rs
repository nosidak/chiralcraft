//! ChiralCraft Floating Panels & Toolbars:
//! - Interactive ¹H-NMR & ¹³C-NMR Spectrum Viewer (reversed ppm axis, Lorentzian plot, clickable peak↔atom assignment)
//! - HELM / Biopolymer Monomer Toolbar (loaded from official `natural-monomers.json`!)
//! - Structure Toolbar (`Check`, `Clean Up`, `S→N`, `N→S`, `¹H-NMR`, `¹³C-NMR`, `Stoich`, `AutoNum`, `Abs/OR/AND`)
//! - BioGraphics Toolbar (BioGraphics templates + Gel Electrophoresis Plate)
//! - Ring Fill & Highlight Color Picker Toolbar
//! - Live Reaction Stoichiometry Grid Window

use chiralcraft_chem::NmrPrediction;
use chiralcraft_doc::{
    format_stoichiometry_text, CopyBackground, DiffStatus, DocumentDiff, ExportOptions, SceneTable,
    SceneTlcPlate, StoichiometryRow,
};
use chiralcraft_formats::{JcampSpectrum, MonomerLibrary, PolymerKind};
use egui::{Color32, FontId, Pos2, RichText, Sense, Stroke, Vec2};

use crate::palettes::CD_WINDOW_BG;

#[derive(Debug, Clone)]
pub enum ProPanelAction {
    None,
    PasteText(String),
    HighlightAtom(Option<usize>),
    ClearNmrAssignment,
    StampNmrShifts(bool),
    ClearStampedNmrShifts,
    LoadJcampOverlay(String),
    ClearJcampOverlay,
    InsertBiopolymer { seq: String, kind: PolymerKind },
    InsertMonomerAtom { polymer_type: String, symbol: String },
    CleanUpBiopolymer,
    CopyHelm,
    CopyFasta,
    StructureCmd(&'static str),
    OpenBioGraphicsLibrary(&'static str),
    InsertGelPlate,
    SetRingFill(Option<Color32>),
    SetHighlightColor(Option<Color32>),
    ExportStoichiometryCsv(String),
}

/// Interactive ChiralCraft ¹H-NMR / ¹³C-NMR Spectrum Viewer Window
/// with optional JCAMP-DX (`.jdx`) experimental spectrum overlay.
pub fn show_nmr_window(
    ctx: &egui::Context,
    open: &mut bool,
    pred: &NmrPrediction,
    selected_atom: &mut Option<usize>,
    jcamp_overlay: Option<&JcampSpectrum>,
    jcamp_path_buf: &mut String,
) -> ProPanelAction {
    let mut act = ProPanelAction::None;
    let title = format!(
        "NMR {} Spectrum Estimation ({:.0} MHz)",
        if pred.nucleus == "1H" { "¹H" } else { "¹³C" },
        pred.freq_mhz
    );
    let (min_ppm, max_ppm) = if pred.nucleus == "1H" { (0.0, 12.5) } else { (0.0, 220.0) };

    egui::Window::new(title)
        .open(open)
        .default_pos([260.0, 74.0])
        .default_width(540.0)
        .default_height(360.0)
        .constrain_to(ctx.available_rect())
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            // 1. Reversed-axis IUPAC NMR spectrum plot (max_ppm on LEFT -> 0 ppm on RIGHT!)
            let (rect, _) = ui.allocate_exact_size(Vec2::new(515.0, 155.0), Sense::hover());
            let p = ui.painter();
            p.rect_filled(rect, 2.0, Color32::WHITE);
            p.rect_stroke(rect, 2.0, Stroke::new(1.0_f32, Color32::from_gray(140)), egui::StrokeKind::Inside);

            let plot = rect.shrink2(Vec2::new(16.0, 18.0));
            let base_y = plot.max.y - 10.0;
            let ppm_to_x = |ppm: f64| -> f32 {
                let frac = ((max_ppm - ppm) / (max_ppm - min_ppm)).clamp(0.0, 1.0) as f32;
                plot.min.x + frac * plot.width()
            };

            // Baseline + ticks
            p.line_segment([Pos2::new(plot.min.x, base_y), Pos2::new(plot.max.x, base_y)], Stroke::new(1.0_f32, Color32::BLACK));
            let tick_step = if pred.nucleus == "1H" { 1.0 } else { 20.0 };
            let mut tppm = 0.0;
            while tppm <= max_ppm + 0.01 {
                let tx = ppm_to_x(tppm);
                p.line_segment([Pos2::new(tx, base_y), Pos2::new(tx, base_y + 4.0)], Stroke::new(1.0_f32, Color32::BLACK));
                p.text(Pos2::new(tx, base_y + 6.0), egui::Align2::CENTER_TOP, format!("{tppm:.0}"), FontId::proportional(9.0), Color32::BLACK);
                tppm += tick_step;
            }
            p.text(Pos2::new(plot.center().x, rect.max.y - 2.0), egui::Align2::CENTER_BOTTOM, "Chemical Shift (ppm)", FontId::proportional(9.5), Color32::DARK_GRAY);

            // Optional JCAMP-DX Experimental NMR Overlay (emerald green trace behind simulated blue trace!)
            let h_scale = (base_y - plot.min.y - 10.0).max(20.0);
            if let Some(exp) = jcamp_overlay {
                let epts: Vec<Pos2> = exp
                    .points
                    .iter()
                    .filter(|&&(ppm, y)| ppm.is_finite() && y.is_finite() && ppm >= min_ppm && ppm <= max_ppm)
                    .map(|&(ppm, y)| Pos2::new(ppm_to_x(ppm), base_y - (y.clamp(0.0, 1.0) as f32) * h_scale * 0.92))
                    .collect();
                for w in epts.windows(2) {
                    p.line_segment([w[0], w[1]], Stroke::new(1.1_f32, Color32::from_rgba_unmultiplied(20, 155, 85, 190)));
                }
                p.text(
                    Pos2::new(plot.min.x + 4.0, plot.min.y + 2.0),
                    egui::Align2::LEFT_TOP,
                    format!("Exp: {}", exp.title),
                    FontId::proportional(8.5),
                    Color32::from_rgb(15, 135, 70),
                );
            }

            // Simulated Lorentzian spectrum trace
            let curve = pred.simulate_spectrum(min_ppm, max_ppm, 360);
            let spts: Vec<Pos2> = curve
                .iter()
                .map(|&(ppm, y)| Pos2::new(ppm_to_x(ppm), base_y - (y as f32) * h_scale))
                .collect();
            for w in spts.windows(2) {
                p.line_segment([w[0], w[1]], Stroke::new(1.2_f32, Color32::from_rgb(15, 50, 150)));
            }

            // Highlight selected peak in red
            if let Some(sel_ai) = *selected_atom {
                for pk in pred.peaks.iter().filter(|pk| pk.atom_idx == sel_ai) {
                    let x = ppm_to_x(pk.shift_ppm);
                    p.line_segment([Pos2::new(x, plot.min.y), Pos2::new(x, base_y)], Stroke::new(1.5_f32, Color32::RED));
                    p.text(Pos2::new(x, plot.min.y + 2.0), egui::Align2::CENTER_TOP, format!("{:.2}", pk.shift_ppm), FontId::proportional(9.5), Color32::RED);
                }
            }

            ui.add_space(4.0);
            ui.label(RichText::new("Click any peak row to highlight its assigned atom on the canvas:").small());
            egui::ScrollArea::vertical().max_height(125.0).show(ui, |ui| {
                egui::Grid::new("nmr_peak_table").striped(true).num_columns(5).spacing([12.0, 2.0]).show(ui, |ui| {
                    ui.strong("Atom #");
                    ui.strong("Shift (ppm)");
                    ui.strong("Mult.");
                    ui.strong("H / C");
                    ui.strong("Environment");
                    ui.end_row();
                    for pk in &pred.peaks {
                        let is_sel = *selected_atom == Some(pk.atom_idx);
                        if ui.selectable_label(is_sel, format!("#{}", pk.atom_idx + 1)).clicked() {
                            *selected_atom = if is_sel { None } else { Some(pk.atom_idx) };
                            act = ProPanelAction::HighlightAtom(*selected_atom);
                        }
                        ui.label(if pred.nucleus == "1H" { format!("{:.2}", pk.shift_ppm) } else { format!("{:.1}", pk.shift_ppm) });
                        ui.label(&pk.multiplicity);
                        ui.label(format!("{}", pk.protons));
                        ui.label(&pk.environment);
                        ui.end_row();
                    }
                });
            });
            ui.separator();
            ui.horizontal(|ui| {
                if ui.button("Paste Protocol").clicked() {
                    act = ProPanelAction::PasteText(pred.report_text.clone());
                }
                if ui.button("Stamp δ Labels on Canvas").clicked() {
                    act = ProPanelAction::StampNmrShifts(pred.nucleus == "13C");
                }
                if ui.button("Clear δ Labels").clicked() {
                    act = ProPanelAction::ClearStampedNmrShifts;
                }
                if ui.button("Clear Assignment").clicked() {
                    *selected_atom = None;
                    act = ProPanelAction::ClearNmrAssignment;
                }
            });
            ui.horizontal(|ui| {
                ui.label("JCAMP-DX (.jdx):");
                ui.add(egui::TextEdit::singleline(jcamp_path_buf).desired_width(175.0).hint_text("experimental.jdx"));
                if ui.button("Browse...").clicked() {
                    if let Some(p) = crate::native_dialogs::pick_open_file(
                        "Open Experimental NMR JCAMP-DX",
                        &["jdx", "dx", "jcamp"],
                    ) {
                        *jcamp_path_buf = p.clone();
                        act = ProPanelAction::LoadJcampOverlay(p);
                    }
                }
                if ui.button("Overlay").clicked() && !jcamp_path_buf.trim().is_empty() {
                    act = ProPanelAction::LoadJcampOverlay(jcamp_path_buf.trim().to_string());
                }
                if jcamp_overlay.is_some() && ui.button("Clear Overlay").clicked() {
                    act = ProPanelAction::ClearJcampOverlay;
                }
            });
        });
    act
}

/// Floating `HELM Monomers` Biopolymer Toolbar (`View > Show HELM Monomer Toolbar`).
pub fn show_helm_monomer_window(
    ctx: &egui::Context,
    open: &mut bool,
    lib: &MonomerLibrary,
    seq_buf: &mut String,
    kind: &mut PolymerKind,
) -> ProPanelAction {
    let mut act = ProPanelAction::None;
    egui::Window::new("HELM Monomers")
        .open(open)
        .default_width(305.0)
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(kind, PolymerKind::Peptide, "PEPTIDE");
                ui.selectable_value(kind, PolymerKind::Rna, "RNA");
                ui.selectable_value(kind, PolymerKind::Dna, "DNA");
            });
            ui.horizontal(|ui| {
                ui.label("Seq:");
                ui.text_edit_singleline(seq_buf);
                if ui.button("Build").clicked() && !seq_buf.trim().is_empty() {
                    act = ProPanelAction::InsertBiopolymer {
                        seq: seq_buf.clone(),
                        kind: *kind,
                    };
                }
            });
            ui.separator();
            let target_type = match *kind {
                PolymerKind::Peptide => "PEPTIDE",
                PolymerKind::Rna | PolymerKind::Dna => "RNA",
            };
            ui.label(RichText::new("Click to append to sequence • Shift-click to stamp all-atom monomer:").small());
            egui::ScrollArea::vertical().max_height(155.0).show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    let shift = ui.input(|i| i.modifiers.shift);
                    for entry in lib.entries.iter().filter(|e| e.polymer_type == target_type).take(48) {
                        let btn = ui.add_sized([34.0, 22.0], egui::Button::new(&entry.symbol))
                            .on_hover_text(format!("{} ({})", entry.name, entry.symbol));
                        if btn.clicked() {
                            if shift {
                                act = ProPanelAction::InsertMonomerAtom {
                                    polymer_type: entry.polymer_type.clone(),
                                    symbol: entry.symbol.clone(),
                                };
                            } else {
                                if entry.symbol.len() > 1 || seq_buf.contains('-') {
                                    if !seq_buf.is_empty() && !seq_buf.ends_with('-') {
                                        if !seq_buf.contains('-') && seq_buf.chars().all(|c| c.is_ascii_alphabetic()) {
                                            let expanded: String = seq_buf
                                                .chars()
                                                .map(|c| c.to_string())
                                                .collect::<Vec<_>>()
                                                .join("-");
                                            *seq_buf = expanded;
                                        }
                                        seq_buf.push('-');
                                    }
                                }
                                seq_buf.push_str(&entry.symbol);
                            }
                        }
                    }
                });
            });
            ui.separator();
            ui.horizontal(|ui| {
                if ui.small_button("Clean Sequence").clicked() { act = ProPanelAction::CleanUpBiopolymer; }
                if ui.small_button("Copy HELM").clicked() { act = ProPanelAction::CopyHelm; }
                if ui.small_button("Copy FASTA").clicked() { act = ProPanelAction::CopyFasta; }
            });
        });
    act
}

/// Top docked `Structure Toolbar` (`View > Show Structure Toolbar`).
pub fn show_structure_toolbar(ctx: &egui::Context) -> ProPanelAction {
    let mut act = ProPanelAction::None;
    egui::TopBottomPanel::top("structure_toolbar")
        .frame(egui::Frame::side_top_panel(&ctx.style()).fill(CD_WINDOW_BG).inner_margin(egui::Margin::symmetric(6, 2)))
        .show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                for (lbl, cmd, tip) in [
                    ("Check", "check", "Check Structure"),
                    ("Clean Up", "cleanup", "Clean Up Structure (Shift+Ctrl+K)"),
                    ("Clean Rxn", "clean_rxn", "Clean Up Reaction (Shift+Ctrl+X)"),
                    ("S→N", "struct_to_name", "Convert Structure to Name (Alt+Ctrl+N)"),
                    ("N→S", "name_to_struct", "Convert Name to Structure (Shift+Ctrl+N)"),
                    ("¹H-NMR", "nmr_1h", "Predict ¹H-NMR Shifts"),
                    ("¹³C-NMR", "nmr_13c", "Predict ¹³C-NMR Shifts"),
                    ("Map Rxn", "map_rxn", "Map Reaction Atoms"),
                    ("AutoNum", "autonum", "Autonumber Reaction"),
                    ("Stoich Grid", "stoich", "Open Reaction Stoichiometry Grid"),
                    ("abs", "stereo_abs", "Enhanced Stereochemistry: Absolute"),
                    ("or1", "stereo_or", "Enhanced Stereochemistry: Create OR Group"),
                    ("&1", "stereo_and", "Enhanced Stereochemistry: Create AND Group"),
                ] {
                    if ui.small_button(lbl).on_hover_text(tip).clicked() {
                        act = ProPanelAction::StructureCmd(cmd);
                    }
                }
            });
        });
    act
}

/// Floating `BioGraphics Toolbar` (`View > Show BioGraphics Toolbar`).
pub fn show_biographics_toolbar(ctx: &egui::Context, open: &mut bool) -> ProPanelAction {
    let mut act = ProPanelAction::None;
    egui::Window::new("BioGraphics Toolbar")
        .open(open)
        .resizable(false)
        .collapsible(false)
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                for (lbl, lib_name) in [
                    ("BioGraphics & Lab", "BioGraphics & Lab Shapes"),
                    ("Lab Glassware", "Lab Glassware & Apparatus"),
                    ("Haworth & Macro", "Haworth & Macrocycles"),
                    ("Catalysts", "Organometallic Catalysts"),
                    ("Amino Acids", "Amino Acids"),
                    ("Carbohydrates", "Carbohydrates"),
                    ("Nucleobases", "Nucleobases"),
                    ("Aromatics", "Aromatics & Heterocycles"),
                    ("Pharmaceuticals", "Pharmaceuticals"),
                ] {
                    if ui.small_button(lbl).clicked() {
                        act = ProPanelAction::OpenBioGraphicsLibrary(lib_name);
                    }
                }
                if ui.small_button("Gel Plate (GEP)").clicked() {
                    act = ProPanelAction::InsertGelPlate;
                }
            });
        });
    act
}

/// Floating `Ring Fill & Highlight Colors` Toolbar (`View > Show Ring Fill Colors Toolbar`).
pub fn show_ring_fill_toolbar(ctx: &egui::Context, open: &mut bool) -> ProPanelAction {
    let mut act = ProPanelAction::None;
    let swatches = [
        ("Pale Blue", Color32::from_rgba_unmultiplied(180, 215, 255, 160)),
        ("Pale Green", Color32::from_rgba_unmultiplied(185, 240, 195, 160)),
        ("Pale Yellow", Color32::from_rgba_unmultiplied(255, 245, 170, 175)),
        ("Pale Pink", Color32::from_rgba_unmultiplied(255, 200, 210, 165)),
        ("Pale Orange", Color32::from_rgba_unmultiplied(255, 220, 175, 165)),
        ("Pale Lavender", Color32::from_rgba_unmultiplied(220, 200, 255, 165)),
        ("Pale Gray", Color32::from_rgba_unmultiplied(215, 215, 220, 165)),
    ];
    egui::Window::new("Ring Fill & Highlight Colors")
        .open(open)
        .resizable(false)
        .collapsible(false)
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            ui.label(RichText::new("Ring Fill Color:").small().strong());
            ui.horizontal(|ui| {
                if ui.small_button("None").clicked() {
                    act = ProPanelAction::SetRingFill(None);
                }
                for (name, col) in swatches {
                    let (r, resp) = ui.allocate_exact_size(Vec2::splat(18.0), Sense::click());
                    ui.painter().rect_filled(r, 1.0, col);
                    ui.painter().rect_stroke(r, 1.0, Stroke::new(1.0_f32, Color32::GRAY), egui::StrokeKind::Inside);
                    if resp.on_hover_text(format!("Ring Fill: {name}")).clicked() {
                        act = ProPanelAction::SetRingFill(Some(col));
                    }
                }
            });
            ui.add_space(4.0);
            ui.label(RichText::new("Highlight Color:").small().strong());
            ui.horizontal(|ui| {
                if ui.small_button("None").clicked() {
                    act = ProPanelAction::SetHighlightColor(None);
                }
                for (name, col) in swatches {
                    let (r, resp) = ui.allocate_exact_size(Vec2::splat(18.0), Sense::click());
                    ui.painter().rect_filled(r, 1.0, col);
                    ui.painter().rect_stroke(r, 1.0, Stroke::new(1.0_f32, Color32::GRAY), egui::StrokeKind::Inside);
                    if resp.on_hover_text(format!("Highlight: {name}")).clicked() {
                        act = ProPanelAction::SetHighlightColor(Some(col));
                    }
                }
            });
        });
    act
}

/// Live Reaction `Stoichiometry Grid` Window.
pub fn show_stoichiometry_window(
    ctx: &egui::Context,
    open: &mut bool,
    base_mmol: &mut f64,
    rows: &mut Vec<StoichiometryRow>,
) -> ProPanelAction {
    let mut act = ProPanelAction::None;
    egui::Window::new("Stoichiometry Grid")
        .open(open)
        .default_width(590.0)
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Limiting Scale (mmol):");
                if ui.add(egui::DragValue::new(base_mmol).speed(0.05).range(0.001..=10000.0)).changed() {
                    for r in rows.iter_mut() {
                        r.recalculate(*base_mmol);
                    }
                }
            });
            ui.separator();
            if rows.is_empty() {
                ui.label("Draw reactants and products (separated by a reaction arrow) to populate the stoichiometry table.");
            } else {
                egui::Grid::new("stoich_grid").striped(true).num_columns(8).spacing([8.0, 3.0]).show(ui, |ui| {
                    for h in ["Compound", "Role", "Formula", "MW", "Eq.", "mmol", "Mass (mg)", "Yield (%)"] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for r in rows.iter_mut() {
                        ui.label(&r.label);
                        ui.label(&r.role);
                        ui.label(&r.formula);
                        ui.label(format!("{:.2}", r.mw));
                        let mut changed = ui.add(egui::DragValue::new(&mut r.equivalents).speed(0.05).range(0.01..=100.0)).changed();
                        ui.label(format!("{:.3}", r.mmol));
                        ui.label(format!("{:.1}", r.mass_mg));
                        if r.role == "Product" {
                            changed |= ui.add(egui::DragValue::new(&mut r.yield_pct).speed(1.0).range(0.0..=100.0)).changed();
                        } else {
                            ui.label("—");
                        }
                        if changed {
                            r.recalculate(*base_mmol);
                        }
                        ui.end_row();
                    }
                });
            }
            ui.separator();
            ui.horizontal(|ui| {
                if ui.button("Paste Table on Canvas").clicked() && !rows.is_empty() {
                    act = ProPanelAction::PasteText(format_stoichiometry_text(rows, false));
                }
                if ui.button("Export CSV").clicked() && !rows.is_empty() {
                    act = ProPanelAction::ExportStoichiometryCsv(format_stoichiometry_text(rows, true));
                }
            });
        });
    act
}

/// Interactive Object Inspector for TLC Plates (editable spots/Rf/colors/solvent),
/// Tables (editable cell grid), 3D Depth Cueing ("Rear Opacity"), and Crossing-Bond Clearance.
pub fn show_object_inspector_window(
    ctx: &egui::Context,
    open: &mut bool,
    tables: &mut [SceneTable],
    tlc_plates: &mut [SceneTlcPlate],
    rear_opacity: &mut f32,
    crossing_clearance_pt: &mut f64,
    export_opts: &mut ExportOptions,
) -> bool {
    let mut modified = false;
    egui::Window::new("Object Inspector (TLC / Tables / 3D / Export)")
        .open(open)
        .default_pos([420.0, 120.0])
        .default_width(370.0)
        .constrain_to(ctx.available_rect())
        .resizable(true)
        .show(ctx, |ui| {
            egui::CollapsingHeader::new("Rendering & 3D Depth")
                .default_open(true)
                .show(ui, |ui| {
                    ui.horizontal(|ui| {
                        ui.label("3D Rear Opacity:");
                        if ui.add(egui::Slider::new(rear_opacity, 0.15..=1.0).fixed_decimals(2)).changed() {
                            modified = true;
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("Crossing Bond Gap (pt):");
                        if ui.add(egui::Slider::new(crossing_clearance_pt, 0.0..=8.0).fixed_decimals(1)).changed() {
                            modified = true;
                        }
                    });
                });

            ui.separator();
            egui::CollapsingHeader::new("Export & Copy Preferences")
                .default_open(true)
                .show(ui, |ui| {
                    ui.horizontal(|ui| {
                        ui.label("PNG DPI:");
                        for dpi in [150u32, 300, 600] {
                            if ui.selectable_label(export_opts.png_dpi == dpi, format!("{dpi}")).clicked() {
                                export_opts.png_dpi = dpi;
                                modified = true;
                            }
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("Padding (pt):");
                        if ui.add(egui::Slider::new(&mut export_opts.svg_padding_pt, 8.0..=60.0).fixed_decimals(0)).changed() {
                            modified = true;
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("Copy Background:");
                        for (lbl, bg) in [
                            ("Transparent", CopyBackground::Transparent),
                            ("White", CopyBackground::White),
                            ("Dark Slide", CopyBackground::DarkSlide),
                        ] {
                            if ui.selectable_label(export_opts.copy_bg == bg, lbl).clicked() {
                                export_opts.copy_bg = bg;
                                modified = true;
                            }
                        }
                    });
                });

            ui.separator();
            egui::CollapsingHeader::new(format!("TLC Plates ({})", tlc_plates.len()))
                .default_open(true)
                .show(ui, |ui| {
                    if tlc_plates.is_empty() {
                        ui.weak("No TLC plates on canvas. Use the TLC tool or Edit > Insert Object.");
                    }
                    let palette: [(&str, (u8, u8, u8)); 6] = [
                        ("Ink", (20, 20, 20)),
                        ("Blue", (25, 95, 205)),
                        ("Red", (205, 45, 45)),
                        ("Green", (25, 145, 70)),
                        ("Violet", (135, 50, 175)),
                        ("Amber", (215, 130, 20)),
                    ];
                    for (pi, tlc) in tlc_plates.iter_mut().enumerate() {
                        tlc.ensure_metadata();
                        ui.group(|ui| {
                            ui.horizontal(|ui| {
                                ui.strong(format!("TLC Plate #{}", pi + 1));
                                if ui.small_button("+ Spot").clicked() {
                                    let idx = tlc.spot_rfs.len() + 1;
                                    tlc.spot_rfs.push(0.50);
                                    tlc.spot_labels.push(format!("L{idx}"));
                                    tlc.spot_colors.push((20, 20, 20));
                                    modified = true;
                                }
                            });
                            ui.horizontal(|ui| {
                                ui.label("Eluent:");
                                if ui.text_edit_singleline(&mut tlc.solvent_label).changed() {
                                    modified = true;
                                }
                            });
                            ui.horizontal(|ui| {
                                ui.label("Baseline:");
                                if ui.add(egui::Slider::new(&mut tlc.origin_fraction, 0.05..=0.35).fixed_decimals(2)).changed() {
                                    modified = true;
                                }
                            });
                            let mut remove_idx = None;
                            for si in 0..tlc.spot_rfs.len() {
                                ui.horizontal(|ui| {
                                    if let Some(lbl) = tlc.spot_labels.get_mut(si) {
                                        if ui.add(egui::TextEdit::singleline(lbl).desired_width(48.0)).changed() {
                                            modified = true;
                                        }
                                    }
                                    ui.label("Rf:");
                                    if ui.add(egui::DragValue::new(&mut tlc.spot_rfs[si]).speed(0.01).range(0.02..=0.98).fixed_decimals(2)).changed() {
                                        modified = true;
                                    }
                                    if let Some(col) = tlc.spot_colors.get_mut(si) {
                                        for (cname, rgb) in palette {
                                            let active = *col == rgb;
                                            let c32 = egui::Color32::from_rgb(rgb.0, rgb.1, rgb.2);
                                            let btn = egui::Button::new(" ")
                                                .fill(c32)
                                                .stroke(if active {
                                                    egui::Stroke::new(2.0_f32, egui::Color32::YELLOW)
                                                } else {
                                                    egui::Stroke::new(0.5_f32, egui::Color32::GRAY)
                                                })
                                                .min_size(egui::vec2(14.0, 14.0));
                                            if ui.add(btn).on_hover_text(cname).clicked() {
                                                *col = rgb;
                                                modified = true;
                                            }
                                        }
                                    }
                                    if tlc.spot_rfs.len() > 1 && ui.small_button("×").clicked() {
                                        remove_idx = Some(si);
                                    }
                                });
                            }
                            if let Some(ri) = remove_idx {
                                tlc.spot_rfs.remove(ri);
                                if ri < tlc.spot_labels.len() { tlc.spot_labels.remove(ri); }
                                if ri < tlc.spot_colors.len() { tlc.spot_colors.remove(ri); }
                                modified = true;
                            }
                        });
                    }
                });

            ui.separator();
            egui::CollapsingHeader::new(format!("Tables ({})", tables.len()))
                .default_open(true)
                .show(ui, |ui| {
                    if tables.is_empty() {
                        ui.weak("No tables on canvas. Use the Table tool or Edit > Insert Object.");
                    }
                    for (ti, tbl) in tables.iter_mut().enumerate() {
                        tbl.ensure_cells();
                        ui.group(|ui| {
                            ui.horizontal(|ui| {
                                ui.strong(format!("Table #{} ({}×{})", ti + 1, tbl.rows, tbl.cols));
                                if ui.small_button("+Row").clicked() && tbl.rows < 12 {
                                    tbl.rows += 1;
                                    tbl.ensure_cells();
                                    modified = true;
                                }
                                if ui.small_button("-Row").clicked() && tbl.rows > 1 {
                                    tbl.rows -= 1;
                                    tbl.ensure_cells();
                                    modified = true;
                                }
                                if ui.small_button("+Col").clicked() && tbl.cols < 8 {
                                    tbl.cols += 1;
                                    tbl.ensure_cells();
                                    modified = true;
                                }
                                if ui.small_button("-Col").clicked() && tbl.cols > 1 {
                                    tbl.cols -= 1;
                                    tbl.ensure_cells();
                                    modified = true;
                                }
                            });
                            egui::Grid::new(format!("tbl_edit_{ti}")).spacing([4.0, 3.0]).show(ui, |ui| {
                                for r in 0..tbl.rows {
                                    for c in 0..tbl.cols {
                                        if let Some(cell) = tbl.cells.get_mut(r).and_then(|row| row.get_mut(c)) {
                                            if ui.add(egui::TextEdit::singleline(cell).desired_width(68.0)).changed() {
                                                modified = true;
                                            }
                                        }
                                    }
                                    ui.end_row();
                                }
                            });
                        });
                    }
                });
        });
    modified
}

pub enum DiffCompareAction {
    None,
    CompareTabs(usize, usize),
    CompareActiveWithFile(String),
}

/// Chemistry-Aware Document Diff / Compare Window (compares open tabs or current tab vs a file on disk).
pub fn show_diff_compare_window(
    ctx: &egui::Context,
    open: &mut bool,
    tab_names: &[String],
    left_tab: &mut usize,
    right_tab: &mut usize,
    file_path_buf: &mut String,
    last_diff: Option<&DocumentDiff>,
) -> DiffCompareAction {
    let mut action = DiffCompareAction::None;
    egui::Window::new("Compare Documents (Chemistry-Aware Diff)")
        .open(open)
        .default_pos([240.0, 95.0])
        .default_width(540.0)
        .constrain_to(ctx.available_rect())
        .resizable(true)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Left Tab:");
                if !tab_names.is_empty() {
                    *left_tab = (*left_tab).min(tab_names.len() - 1);
                    *right_tab = (*right_tab).min(tab_names.len() - 1);
                    egui::ComboBox::from_id_salt("diff_left_tab")
                        .selected_text(&tab_names[*left_tab])
                        .show_ui(ui, |ui| {
                            for (i, n) in tab_names.iter().enumerate() {
                                ui.selectable_value(left_tab, i, n);
                            }
                        });
                    ui.label("Right Tab:");
                    egui::ComboBox::from_id_salt("diff_right_tab")
                        .selected_text(&tab_names[*right_tab])
                        .show_ui(ui, |ui| {
                            for (i, n) in tab_names.iter().enumerate() {
                                ui.selectable_value(right_tab, i, n);
                            }
                        });
                    if ui.button("Compare Tabs").clicked() {
                        action = DiffCompareAction::CompareTabs(*left_tab, *right_tab);
                    }
                }
            });
            ui.horizontal(|ui| {
                ui.label("Or compare active tab with file:");
                ui.add(egui::TextEdit::singleline(file_path_buf).desired_width(220.0).hint_text("path/to/other.cdxml or .svg"));
                if ui.button("Compare File").clicked() && !file_path_buf.trim().is_empty() {
                    action = DiffCompareAction::CompareActiveWithFile(file_path_buf.trim().to_string());
                }
            });
            ui.separator();
            if let Some(diff) = last_diff {
                ui.horizontal(|ui| {
                    ui.strong(format!("{}  vs  {}", diff.left_name, diff.right_name));
                    if !diff.has_changes() {
                        ui.label(RichText::new("✓ Chemically identical").color(Color32::from_rgb(20, 145, 70)));
                    }
                });
                egui::ScrollArea::vertical().max_height(260.0).show(ui, |ui| {
                    for e in &diff.entries {
                        let (tag, col) = match e.status {
                            DiffStatus::Unchanged => ("=", Color32::GRAY),
                            DiffStatus::Added => ("+ ADDED", Color32::from_rgb(20, 145, 70)),
                            DiffStatus::Removed => ("- REMOVED", Color32::from_rgb(200, 45, 45)),
                            DiffStatus::Modified => ("~ MODIFIED", Color32::from_rgb(205, 130, 15)),
                        };
                        ui.horizontal(|ui| {
                            ui.label(RichText::new(tag).strong().color(col));
                            ui.label(&e.summary);
                        });
                    }
                    for c in &diff.added_captions {
                        ui.label(RichText::new(format!("+ Caption: {c}")).color(Color32::from_rgb(20, 145, 70)));
                    }
                    for c in &diff.removed_captions {
                        ui.label(RichText::new(format!("- Caption: {c}")).color(Color32::from_rgb(200, 45, 45)));
                    }
                    if diff.arrow_delta != 0 || diff.table_delta != 0 || diff.tlc_delta != 0 {
                        ui.label(format!(
                            "Objects delta — Arrows: {:+}, Tables: {:+}, TLC Plates: {:+}",
                            diff.arrow_delta, diff.table_delta, diff.tlc_delta
                        ));
                    }
                });
            } else {
                ui.weak("Select two open tabs or enter a file path and click Compare.");
            }
        });
    action
}
