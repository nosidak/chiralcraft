//! Pixel-faithful ChiralCraft floating palettes matching live captures:
//! - `Tools` (`tools-20261118.png`): 18x2 iconic grid with tear-off flyout sub-palettes & vector `.ctp` previews
//! - `Periodic Table` (`p23-periodic-table.png`): full 18-col + f-block color-coded table
//! - `Analysis` (`p08-analysis-0.png`): checkboxes, Decimals spinner, mini m/z plot, `[Paste]`
//! - `Chemical Properties` (`p08-chem-properties.png`): all 13 checkboxes, `[Paste]`, `[Report]`
//! - `Info` (`p25-info-window.png`): Pointer X/Y/dX/dY, 3x3 Selection anchor grid, Other Angle/Dist/%
//! - `Character Map`: Greek & scientific chemistry symbols (`View > Show Character Map Window`)

use chiralcraft_chem::{compute_analysis_report, estimate_chemical_properties, get_element, ChemMolecule};
use chiralcraft_doc::{SceneCurve, TemplateLibrary};
use egui::{Color32, FontId, Pos2, Rect, RichText, Sense, Stroke, Vec2};

use crate::icons::{draw_steel_blue_gradient, format_chemical_formula, tool_icon_button, Tool};

pub const CD_WINDOW_BG: Color32 = Color32::from_rgb(240, 240, 240);
pub const CD_DESK_GRAY: Color32 = Color32::from_rgb(160, 160, 160);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActiveFlyout {
    None,
    Accessories,
    Arrows,
    Brackets,
    Chains,
    ChemSymbols,
    Chromatography,
    Colors,
    Curves,
    DottedDashed,
    Shapes,
    Fragmentation,
    MultipleBonds,
    Orbitals,
    PenTools,
    Windows,
    Templates,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PaletteCommand {
    None,
    StampTemplate(usize, usize),
    SetColor(Color32),
    ToggleColorByElement,
    OpenColorPicker,
    ToggleCurve(&'static str),
    WindowAction(&'static str),
}

pub struct AnalysisPanelState {
    pub chk_name: bool,
    pub chk_formula: bool,
    pub chk_exact: bool,
    pub chk_mw: bool,
    pub chk_mz: bool,
    pub chk_elem: bool,
    pub decimals: usize,
}

impl Default for AnalysisPanelState {
    fn default() -> Self {
        Self {
            chk_name: true,
            chk_formula: true,
            chk_exact: true,
            chk_mw: true,
            chk_mz: true,
            chk_elem: true,
            decimals: 2,
        }
    }
}

#[derive(Debug, Clone)]
pub enum AnalysisAction {
    None,
    Paste(String),
    PasteNameCentered,
    ResolvePubChemName,
    LinkCaption(chiralcraft_doc::LinkedCaptionKind),
}

pub struct PropsPanelState {
    pub checks: [bool; 13],
}

impl Default for PropsPanelState {
    fn default() -> Self {
        Self { checks: [true; 13] }
    }
}

/// Returns `(fill_color, text_color)` matching `p23-periodic-table.png`.
fn ptable_cell_colors(z: i32) -> (Color32, Color32) {
    match z {
        6 => (Color32::BLACK, Color32::WHITE),
        7 => (Color32::from_rgb(0, 0, 255), Color32::WHITE),
        8 => (Color32::from_rgb(255, 0, 0), Color32::WHITE),
        9 => (Color32::from_rgb(64, 154, 114), Color32::WHITE),
        11 => (Color32::from_rgb(184, 52, 184), Color32::BLACK),
        15 => (Color32::from_rgb(220, 38, 108), Color32::WHITE),
        16 => (Color32::from_rgb(226, 218, 64), Color32::BLACK),
        17 => (Color32::from_rgb(0, 148, 0), Color32::WHITE),
        26 => (Color32::from_rgb(0, 82, 0), Color32::WHITE),
        28 => (Color32::from_rgb(96, 96, 96), Color32::WHITE),
        29 => (Color32::from_rgb(206, 134, 54), Color32::BLACK),
        35 => (Color32::from_rgb(154, 72, 72), Color32::WHITE),
        47 => (Color32::from_rgb(142, 152, 174), Color32::BLACK),
        79 => (Color32::from_rgb(255, 210, 0), Color32::WHITE),
        _ => (Color32::WHITE, Color32::BLACK),
    }
}

/// Full 18-column ChiralCraft Periodic Table window (`p23-periodic-table.png`).
pub fn show_periodic_table_window(
    ctx: &egui::Context,
    open: &mut bool,
    active_element_z: &mut Option<i32>,
) -> Option<i32> {
    let mut picked = None;
    let main_rows: [[i32; 18]; 7] = [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
        [3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 6, 7, 8, 9, 10],
        [11, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13, 14, 15, 16, 17, 18],
        [19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
        [37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54],
        [55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86],
        [87, 88, 89, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118],
    ];
    let f_rows: [[i32; 14]; 2] = [
        [58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71],
        [90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103],
    ];

    egui::Window::new("Periodic Table")
        .open(open)
        .resizable(false)
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            ui.spacing_mut().item_spacing = Vec2::splat(0.0);
            let mut draw_cell = |ui: &mut egui::Ui, z: i32| {
                let (rect, resp) = ui.allocate_exact_size(Vec2::new(22.0, 22.0), Sense::click());
                if z == 0 {
                    return;
                }
                let (fill, fg) = ptable_cell_colors(z);
                let p = ui.painter();
                p.rect_filled(rect, 0.0, fill);
                let border = if *active_element_z == Some(z) {
                    Stroke::new(2.0_f32, Color32::from_rgb(0, 120, 215))
                } else {
                    Stroke::new(1.0_f32, Color32::BLACK)
                };
                p.rect_stroke(rect, 0.0, border, egui::StrokeKind::Inside);
                let sym = get_element(z).map(|e| e.symbol).unwrap_or("?");
                p.text(
                    rect.center(),
                    egui::Align2::CENTER_CENTER,
                    sym,
                    FontId::proportional(10.5),
                    fg,
                );
                if resp.clicked() {
                    *active_element_z = Some(z);
                    picked = Some(z);
                }
                let tip = get_element(z)
                    .map(|e| format!("{} (Z={}, AW {:.3})", e.name, e.z, e.atomic_weight()))
                    .unwrap_or_default();
                resp.on_hover_text(tip);
            };

            for row in &main_rows {
                ui.horizontal(|ui| {
                    for &z in row {
                        draw_cell(ui, z);
                    }
                });
            }
            ui.add_space(6.0);
            for row in &f_rows {
                ui.horizontal(|ui| {
                    ui.add_space(66.0);
                    for &z in row {
                        draw_cell(ui, z);
                    }
                });
            }
        });
    picked
}

/// Floating `Character Map` window (`View > Show Character Map Window`).
/// Returns `Some(sym)` if user clicks a symbol to insert as a caption.
pub fn show_character_map_window(ctx: &egui::Context, open: &mut bool) -> Option<String> {
    let mut picked = None;
    let syms = [
        "α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "λ", "μ", "ν", "π", "ρ", "σ", "τ", "φ", "χ", "ψ", "ω",
        "Δ", "Σ", "Ω", "±", "°", "Å", "·", "•", "†", "‡", "≈", "≠", "≤", "≥", "∞", "→", "←", "↔", "⇌",
        "↑", "↓", "‰", "™", "®", "©", "µ", "×", "÷", "∂",
    ];
    egui::Window::new("Character Map")
        .open(open)
        .resizable(false)
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            egui::Grid::new("charmap_grid").num_columns(8).spacing([2.0, 2.0]).show(ui, |ui| {
                for (i, s) in syms.iter().enumerate() {
                    if ui.add_sized([24.0, 22.0], egui::Button::new(*s)).clicked() {
                        picked = Some((*s).to_string());
                    }
                    if (i + 1) % 8 == 0 {
                        ui.end_row();
                    }
                }
            });
        });
    picked
}

/// Floating `Analysis` window matching `p08-analysis-0.png`, enhanced with live IUPAC / PubChem naming.
pub fn show_analysis_window(
    ctx: &egui::Context,
    open: &mut bool,
    state: &mut AnalysisPanelState,
    mol: Option<&ChemMolecule>,
    live_name: Option<&str>,
    pubchem_pending: bool,
) -> AnalysisAction {
    let mut action = AnalysisAction::None;
    egui::Window::new("Analysis")
        .open(open)
        .resizable(false)
        .default_width(295.0)
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            let Some(m) = mol else {
                ui.label("Select or draw a structure.");
                return;
            };
            let rep = compute_analysis_report(m);
            let d = state.decimals;
            let pretty_formula = format_chemical_formula(&rep.formula);
            if let Some(nm) = live_name {
                ui.horizontal_wrapped(|ui| {
                    ui.checkbox(&mut state.chk_name, "Name:");
                    ui.strong(nm);
                    if pubchem_pending {
                        ui.spinner();
                    }
                });
                ui.separator();
            }
            ui.horizontal(|ui| {
                ui.checkbox(&mut state.chk_formula, "Formula:");
                ui.strong(&pretty_formula);
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut state.chk_exact, "Exact Mass:");
                ui.label(format!("{:.*}", d, rep.exact_mass));
                ui.add_space(16.0);
                ui.label("Decimals:");
                ui.add(egui::DragValue::new(&mut state.decimals).range(0..=6));
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut state.chk_mw, "Mol. Wt.:");
                ui.label(format!("{:.*}", d, rep.mol_weight));
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut state.chk_mz, "m/z:");
                ui.vertical(|ui| {
                    for &(mass, pct) in &rep.mz_peaks {
                        ui.label(format!("{:.*} ({:.1}%)", d, mass, pct));
                    }
                });
                let (srect, _) = ui.allocate_exact_size(Vec2::new(62.0, 34.0), Sense::hover());
                let p = ui.painter();
                p.line_segment(
                    [Pos2::new(srect.min.x, srect.max.y - 2.0), Pos2::new(srect.max.x, srect.max.y - 2.0)],
                    Stroke::new(1.0_f32, Color32::BLACK),
                );
                for (k, &(_, pct)) in rep.mz_peaks.iter().enumerate() {
                    let x = srect.min.x + 10.0 + (k as f32) * 8.0;
                    let h = ((pct / 100.0) as f32) * (srect.height() - 6.0);
                    p.line_segment(
                        [Pos2::new(x, srect.max.y - 2.0), Pos2::new(x, srect.max.y - 2.0 - h)],
                        Stroke::new(1.5_f32, Color32::from_gray(70)),
                    );
                }
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut state.chk_elem, "Elem. Anal.:");
                ui.label(&rep.elem_analysis);
            });
            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                if ui.button("  Paste  ").clicked() {
                    let mut lines = Vec::new();
                    if state.chk_name {
                        if let Some(nm) = live_name {
                            lines.push(format!("Name: {nm}"));
                        }
                    }
                    if state.chk_formula { lines.push(format!("Chemical Formula: {}", pretty_formula)); }
                    if state.chk_exact { lines.push(format!("Exact Mass: {:.*}", d, rep.exact_mass)); }
                    if state.chk_mw { lines.push(format!("Molecular Weight: {:.*}", d, rep.mol_weight)); }
                    if state.chk_elem { lines.push(format!("Elemental Analysis: {}", rep.elem_analysis)); }
                    action = AnalysisAction::Paste(lines.join("\n"));
                }
                if ui.button("Paste Name").on_hover_text("Stamp centered IUPAC name below structure (Alt+Ctrl+N)").clicked() {
                    action = AnalysisAction::PasteNameCentered;
                }
                if ui
                    .add_enabled(!pubchem_pending, egui::Button::new("PubChem Name"))
                    .on_hover_text("Fetch authoritative IUPAC & common drug name from PubChem")
                    .clicked()
                {
                    action = AnalysisAction::ResolvePubChemName;
                }
            });
            ui.add_space(3.0);
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Live Linked:").small().strong());
                if ui.small_button("Link Name").on_hover_text("Attach auto-updating IUPAC name caption below structure").clicked() {
                    action = AnalysisAction::LinkCaption(chiralcraft_doc::LinkedCaptionKind::IupacName);
                }
                if ui.small_button("Link Formula").on_hover_text("Attach auto-updating formula caption").clicked() {
                    action = AnalysisAction::LinkCaption(chiralcraft_doc::LinkedCaptionKind::Formula);
                }
                if ui.small_button("Link MW").on_hover_text("Attach auto-updating MW caption").clicked() {
                    action = AnalysisAction::LinkCaption(chiralcraft_doc::LinkedCaptionKind::MolecularWeight);
                }
                if ui.small_button("Link Summary").on_hover_text("Attach live Name | Formula | MW caption").clicked() {
                    action = AnalysisAction::LinkCaption(chiralcraft_doc::LinkedCaptionKind::FullSummary);
                }
            });
        });
    action
}

/// Floating `Chemical Properties` window matching `p08-chem-properties.png` (all 13 checkboxes + `[Paste]` / `[Report]`).
pub fn show_properties_window(
    ctx: &egui::Context,
    open: &mut bool,
    state: &mut PropsPanelState,
    mol: Option<&ChemMolecule>,
) -> Option<String> {
    let mut paste_text = None;
    egui::Window::new("Chemical Properties")
        .open(open)
        .resizable(false)
        .default_width(245.0)
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            let Some(m) = mol else {
                ui.label("Select or draw a structure.");
                return;
            };
            let pr = estimate_chemical_properties(m);
            let rows: [(&str, String); 13] = [
                ("Boiling Point:", format!("{:.2} [K]", pr.boiling_point_k)),
                ("Melting Point:", format!("{:.2} [K]", pr.melting_point_k)),
                ("Critical Temp:", format!("{:.2} [K]", pr.critical_temp_k)),
                ("Critical Pres:", format!("{:.2} [Bar]", pr.critical_pres_bar)),
                ("Critical Vol:", format!("{:.1} [cm³/mol]", pr.critical_vol)),
                ("Gibbs Energy:", format!("{:.1} [kJ/mol]", pr.gibbs_energy)),
                ("Log P:", format!("{:.2}", pr.log_p)),
                ("MR:", format!("{:.2} [cm³/mol]", pr.molar_refractivity)),
                ("Henry's Law:", format!("{:.2}", pr.henry_law)),
                ("Heat of Form:", format!("{:.2} [kJ/mol]", pr.heat_of_form)),
                ("tPSA:", format!("{:.2}", pr.tpsa)),
                ("LogS:", format!("{:.3}", pr.log_s)),
                ("pKa:", pr.pka.map(|v| format!("{v:.3}")).unwrap_or_else(|| "N/A".to_string())),
            ];
            egui::Grid::new("chem_props_grid").num_columns(2).spacing([8.0, 2.0]).show(ui, |ui| {
                for (i, (lbl, val)) in rows.iter().enumerate() {
                    ui.checkbox(&mut state.checks[i], *lbl);
                    ui.label(val);
                    ui.end_row();
                }
            });
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                if ui.button("   Paste   ").clicked() {
                    let mut lines = Vec::new();
                    for (i, (lbl, val)) in rows.iter().enumerate() {
                        if state.checks[i] {
                            lines.push(format!("{lbl} {val}"));
                        }
                    }
                    paste_text = Some(lines.join("\n"));
                }
                if ui.button("   Report   ").clicked() {
                    let mut lines = Vec::new();
                    for (lbl, val) in &rows {
                        lines.push(format!("{lbl} {val}"));
                    }
                    paste_text = Some(lines.join("\n"));
                }
            });
        });
    paste_text
}

/// Floating `Info` window matching `p25-info-window.png`.
pub fn show_info_window(
    ctx: &egui::Context,
    open: &mut bool,
    cursor_doc: Option<(f64, f64)>,
    drag_delta: Option<(f64, f64)>,
    sel_bbox: Option<(f64, f64, f64, f64)>,
    zoom: f32,
) {
    egui::Window::new("Info")
        .open(open)
        .resizable(false)
        .default_width(135.0)
        .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
        .show(ctx, |ui| {
            let pt_to_cm = 1.0 / 28.3465;
            ui.label(RichText::new("Pointer:").strong());
            let (px, py) = cursor_doc.unwrap_or((0.0, 0.0));
            let (dx, dy) = drag_delta.unwrap_or((0.0, 0.0));
            ui.label(format!("   X: {:.3} cm", px * pt_to_cm));
            ui.label(format!("   Y: {:.3} cm", py * pt_to_cm));
            ui.label(format!("  dX: {:.3} cm", dx * pt_to_cm));
            ui.label(format!("  dY: {:.3} cm", dy * pt_to_cm));
            ui.add_space(6.0);
            ui.label(RichText::new("Selection:").strong());
            ui.horizontal(|ui| {
                let (grect, _) = ui.allocate_exact_size(Vec2::splat(30.0), Sense::hover());
                let p = ui.painter();
                for r in 0..3 {
                    for c in 0..3 {
                        let cell = Rect::from_min_size(
                            grect.min + Vec2::new((c as f32) * 10.0, (r as f32) * 10.0),
                            Vec2::splat(7.0),
                        );
                        if r == 0 && c == 0 {
                            p.rect_filled(cell, 0.0, Color32::BLACK);
                        } else {
                            p.rect_filled(cell, 0.0, Color32::WHITE);
                            p.rect_stroke(cell, 0.0, Stroke::new(1.0_f32, Color32::BLACK), egui::StrokeKind::Inside);
                        }
                    }
                }
                ui.vertical(|ui| {
                    if let Some((x1, y1, x2, y2)) = sel_bbox {
                        ui.label(format!("X: {:.3} cm", x1 * pt_to_cm));
                        ui.label(format!("Y: {:.3} cm", y1 * pt_to_cm));
                        ui.label(format!("W: {:.3} cm", (x2 - x1) * pt_to_cm));
                        ui.label(format!("H: {:.3} cm", (y2 - y1) * pt_to_cm));
                    } else {
                        ui.label("X: —");
                        ui.label("Y: —");
                    }
                });
            });
            ui.add_space(6.0);
            ui.label(RichText::new("Other:").strong());
            let dist = dx.hypot(dy);
            let ang = dy.atan2(dx).to_degrees();
            ui.label(format!("  Angle: {:.1}°", if dist > 0.1 { ang } else { 0.0 }));
            ui.label(format!("  Dist.: {:.2} cm", dist * pt_to_cm));
            ui.label(format!("      %: {:.0}%", zoom * 100.0));
        });
}

/// Draws a miniature vector thumbnail of a `.ctp` template molecule inside a tile!
fn draw_molecule_thumbnail(ui: &mut egui::Ui, mol: &ChemMolecule, idx: usize) -> egui::Response {
    let (rect, resp) = ui.allocate_exact_size(Vec2::splat(46.0), Sense::click());
    let p = ui.painter();
    p.rect_filled(
        rect,
        1.0,
        if resp.hovered() {
            Color32::from_rgb(222, 238, 252)
        } else {
            Color32::WHITE
        },
    );
    p.rect_stroke(rect, 1.0, Stroke::new(1.0_f32, Color32::from_gray(180)), egui::StrokeKind::Inside);

    if !mol.atoms.is_empty() {
        let mut x1 = f64::INFINITY;
        let mut y1 = f64::INFINITY;
        let mut x2 = f64::NEG_INFINITY;
        let mut y2 = f64::NEG_INFINITY;
        for a in &mol.atoms {
            x1 = x1.min(a.x);
            y1 = y1.min(a.y);
            x2 = x2.max(a.x);
            y2 = y2.max(a.y);
        }
        let inner = rect.shrink(6.0);
        let w = (x2 - x1).max(1.0);
        let h = (y2 - y1).max(1.0);
        let scale = ((inner.width() as f64) / w).min((inner.height() as f64) / h);
        let ox = inner.center().x as f64 - 0.5 * (x1 + x2) * scale;
        let oy = inner.center().y as f64 - 0.5 * (y1 + y2) * scale;
        let map = |x: f64, y: f64| Pos2::new((ox + x * scale) as f32, (oy + y * scale) as f32);

        for b in &mol.bonds {
            if b.begin_atom < mol.atoms.len() && b.end_atom < mol.atoms.len() {
                let p1 = map(mol.atoms[b.begin_atom].x, mol.atoms[b.begin_atom].y);
                let p2 = map(mol.atoms[b.end_atom].x, mol.atoms[b.end_atom].y);
                let d = p2 - p1;
                let len = d.length().max(1e-3);
                let n = Vec2::new(-d.y / len, d.x / len);
                if b.order == 2 {
                    p.line_segment([p1 + n * 1.1, p2 + n * 1.1], Stroke::new(1.0_f32, Color32::BLACK));
                    p.line_segment([p1 - n * 1.1, p2 - n * 1.1], Stroke::new(1.0_f32, Color32::BLACK));
                } else {
                    p.line_segment([p1, p2], Stroke::new(1.0_f32, Color32::BLACK));
                }
            }
        }
    }
    let f = format_chemical_formula(&mol.compute_formula_properties().get_formula());
    resp.on_hover_text(format!("#{} ({})", idx + 1, f))
}

/// Draws a miniature vector thumbnail of a graphic/curve `.ctp` template tile!
fn draw_graphic_thumbnail(ui: &mut egui::Ui, curves: &[SceneCurve], idx: usize) -> egui::Response {
    let (rect, resp) = ui.allocate_exact_size(Vec2::splat(46.0), Sense::click());
    let p = ui.painter();
    p.rect_filled(
        rect,
        1.0,
        if resp.hovered() { Color32::from_rgb(222, 238, 252) } else { Color32::WHITE },
    );
    p.rect_stroke(rect, 1.0, Stroke::new(1.0_f32, Color32::from_gray(180)), egui::StrokeKind::Inside);

    let mut x1 = f64::INFINITY;
    let mut y1 = f64::INFINITY;
    let mut x2 = f64::NEG_INFINITY;
    let mut y2 = f64::NEG_INFINITY;
    for c in curves {
        for &(x, y) in &c.points {
            x1 = x1.min(x);
            y1 = y1.min(y);
            x2 = x2.max(x);
            y2 = y2.max(y);
        }
    }
    if x1.is_finite() {
        let inner = rect.shrink(6.0);
        let w = (x2 - x1).max(1.0);
        let h = (y2 - y1).max(1.0);
        let scale = ((inner.width() as f64) / w).min((inner.height() as f64) / h);
        let ox = inner.center().x as f64 - 0.5 * (x1 + x2) * scale;
        let oy = inner.center().y as f64 - 0.5 * (y1 + y2) * scale;
        let map = |x: f64, y: f64| Pos2::new((ox + x * scale) as f32, (oy + y * scale) as f32);
        for c in curves {
            for w in c.points.windows(2) {
                p.line_segment([map(w[0].0, w[0].1), map(w[1].0, w[1].1)], Stroke::new(1.0_f32, Color32::BLACK));
            }
            if c.closed {
                if let (Some(&f), Some(&l)) = (c.points.first(), c.points.last()) {
                    p.line_segment([map(l.0, l.1), map(f.0, f.1)], Stroke::new(1.0_f32, Color32::BLACK));
                }
            }
        }
    }
    resp.on_hover_text(format!("Graphic Template #{}", idx + 1))
}

/// Iconic 18x2 `Tools` palette matching `tools-20261118.png` + all 15 `Other Toolbars` tear-offs + all 26 `.ctp` libraries.
pub fn show_tools_palette(
    ctx: &egui::Context,
    show_main_tools: bool,
    tool: &mut Tool,
    active_element_z: &mut Option<i32>,
    active_flyout: &mut ActiveFlyout,
    ctp_libraries: &[TemplateLibrary],
    selected_ctp_idx: &mut usize,
) -> PaletteCommand {
    let mut cmd = PaletteCommand::None;

    if show_main_tools {
        egui::SidePanel::left("cd_tools_palette")
            .exact_width(64.0)
            .resizable(false)
            .frame(
                egui::Frame::side_top_panel(&ctx.style())
                    .fill(CD_WINDOW_BG)
                    .inner_margin(egui::Margin::same(3)),
            )
            .show(ctx, |ui| {
                // Exact steel-blue gradient "Tools" mini title bar + minimize box (`tools-20261118.png`)
                let (trect, _) = ui.allocate_exact_size(Vec2::new(58.0, 18.0), Sense::hover());
                draw_steel_blue_gradient(ui.painter(), trect);
                ui.painter().text(
                    trect.left_center() + Vec2::new(3.0, 0.0),
                    egui::Align2::LEFT_CENTER,
                    "Tools",
                    FontId::proportional(11.0),
                    Color32::BLACK,
                );
                let min_box = Rect::from_min_size(Pos2::new(trect.max.x - 15.0, trect.min.y + 2.0), Vec2::splat(14.0));
                ui.painter().rect_filled(min_box, 0.0, Color32::from_gray(236));
                ui.painter().rect_stroke(min_box, 0.0, Stroke::new(1.0_f32, Color32::from_gray(120)), egui::StrokeKind::Inside);
                ui.painter().line_segment(
                    [Pos2::new(min_box.min.x + 3.0, min_box.max.y - 3.0), Pos2::new(min_box.max.x - 3.0, min_box.max.y - 3.0)],
                    Stroke::new(1.8_f32, Color32::BLACK),
                );
                ui.add_space(2.0);
                ui.spacing_mut().item_spacing = Vec2::splat(2.0);

                let rows: [((Tool, ActiveFlyout, &str), (Tool, ActiveFlyout, &str)); 18] = [
                    ((Tool::Lasso, ActiveFlyout::None, "Lasso (Freehand Selection)"), (Tool::Marquee, ActiveFlyout::None, "Marquee [Space]")),
                    ((Tool::Perspective, ActiveFlyout::None, "Structure Perspective (3D Tilt Drag)"), (Tool::Fragmentation, ActiveFlyout::Fragmentation, "Mass Fragmentation")),
                    ((Tool::SolidBond, ActiveFlyout::None, "Solid Bond"), (Tool::Eraser, ActiveFlyout::None, "Eraser")),
                    ((Tool::MultipleBond, ActiveFlyout::MultipleBonds, "Multiple Bonds"), (Tool::TextLabel, ActiveFlyout::None, "Text [t]")),
                    ((Tool::DashedBond, ActiveFlyout::DottedDashed, "Dotted/Dashed Bonds"), (Tool::PenCurve, ActiveFlyout::PenTools, "Pen Tools")),
                    ((Tool::HashedBond, ActiveFlyout::None, "Hashed Bond"), (Tool::Arrow, ActiveFlyout::Arrows, "Arrows [e]")),
                    ((Tool::HashWedgeBond, ActiveFlyout::None, "Hashed Wedged Bond"), (Tool::Orbital, ActiveFlyout::Orbitals, "Orbitals [G]")),
                    ((Tool::BoldBond, ActiveFlyout::None, "Bold Bond"), (Tool::ShapeRect, ActiveFlyout::Shapes, "Drawing Elements")),
                    ((Tool::WedgeBond, ActiveFlyout::None, "Wedged Bond"), (Tool::Bracket, ActiveFlyout::Brackets, "Brackets [T]")),
                    ((Tool::HollowWedgeBond, ActiveFlyout::None, "Hollow Wedged Bond"), (Tool::ChemSymbol, ActiveFlyout::ChemSymbols, "Chemical Symbols [E]")),
                    ((Tool::WavyBond, ActiveFlyout::None, "Wavy Bond"), (Tool::ReactionMap, ActiveFlyout::None, "Reaction Atom Map")),
                    ((Tool::Table, ActiveFlyout::None, "Table (2x2)"), (Tool::TlcPlate, ActiveFlyout::Chromatography, "Chromatography Tools [g]")),
                    ((Tool::AcyclicChain, ActiveFlyout::Chains, "Chains [X]"), (Tool::Templates, ActiveFlyout::Templates, "Templates (.ctp Libraries)")),
                    ((Tool::Cyclopropane, ActiveFlyout::None, "Cyclopropane [3]"), (Tool::Cyclobutane, ActiveFlyout::None, "Cyclobutane [4]")),
                    ((Tool::Cyclopentane, ActiveFlyout::None, "Cyclopentane [5]"), (Tool::Cyclohexane, ActiveFlyout::None, "Cyclohexane [6]")),
                    ((Tool::Cycloheptane, ActiveFlyout::None, "Cycloheptane [7]"), (Tool::Cyclooctane, ActiveFlyout::None, "Cyclooctane [8]")),
                    ((Tool::Chair1, ActiveFlyout::None, "Cyclohexane Chair 1"), (Tool::Chair2, ActiveFlyout::None, "Cyclohexane Chair 2")),
                    ((Tool::Cyclopentadiene, ActiveFlyout::None, "Cyclopentadiene [J]"), (Tool::Benzene, ActiveFlyout::None, "Benzene [j]")),
                ];

                egui::ScrollArea::vertical().show(ui, |ui| {
                    for ((lt, lfly, ltip), (rt, rfly, rtip)) in rows {
                        ui.horizontal(|ui| {
                            let l_disp = if lt == Tool::MultipleBond && matches!(*tool, Tool::DoubleBond | Tool::TripleBond | Tool::DativeBond | Tool::QuadrupleBond) {
                                *tool
                            } else {
                                lt
                            };
                            let l_active = *tool == lt
                                || (lt == Tool::MultipleBond && matches!(*tool, Tool::DoubleBond | Tool::TripleBond | Tool::DativeBond | Tool::QuadrupleBond));

                            let r_disp = match rt {
                                Tool::Arrow if matches!(*tool, Tool::Arrow | Tool::EquilibriumArrow | Tool::NoGoArrow | Tool::RetrosyntheticArrow | Tool::CurvedArrow) => *tool,
                                Tool::Orbital if matches!(*tool, Tool::Orbital | Tool::OrbitalS | Tool::OrbitalD) => *tool,
                                Tool::ShapeRect if matches!(*tool, Tool::ShapeRect | Tool::ShapeOval | Tool::ShapeLine) => *tool,
                                Tool::Bracket if matches!(*tool, Tool::Bracket | Tool::RoundBracket | Tool::CurlyBracket) => *tool,
                                Tool::ChemSymbol if matches!(*tool, Tool::ChemSymbol | Tool::SymbolMinus | Tool::SymbolRadical | Tool::SymbolLonePair) => *tool,
                                _ => rt,
                            };
                            let r_active = *tool == r_disp;

                            let lr = tool_icon_button(ui, l_disp, l_active, lfly != ActiveFlyout::None, ltip);
                            if lr.secondary_clicked() && lfly != ActiveFlyout::None {
                                *active_flyout = if *active_flyout == lfly { ActiveFlyout::None } else { lfly };
                            } else if lr.clicked() {
                                *active_element_z = None;
                                if lfly != ActiveFlyout::None {
                                    *active_flyout = if *active_flyout == lfly { ActiveFlyout::None } else { lfly };
                                }
                                if lt == Tool::MultipleBond && !matches!(*tool, Tool::DoubleBond | Tool::TripleBond | Tool::DativeBond | Tool::QuadrupleBond) {
                                    *tool = Tool::DoubleBond;
                                } else if lt != Tool::MultipleBond {
                                    *tool = lt;
                                }
                            }

                            let rr = tool_icon_button(ui, r_disp, r_active, rfly != ActiveFlyout::None, rtip);
                            if rr.secondary_clicked() && rfly != ActiveFlyout::None {
                                *active_flyout = if *active_flyout == rfly { ActiveFlyout::None } else { rfly };
                            } else if rr.clicked() {
                                *active_element_z = None;
                                if rfly != ActiveFlyout::None {
                                    *active_flyout = if *active_flyout == rfly { ActiveFlyout::None } else { rfly };
                                }
                                if rt != Tool::Templates {
                                    *tool = r_disp;
                                }
                            }
                        });
                    }
                });
            });
    }

    // ---- Tear-Off Floating Sub-Palettes (ALL 15 Prime `Other Toolbars` + All 26 `.ctp` Templates) ----
    let show_mini_flyout = |title: &str, y_pos: f32, items: &[(Tool, &str)], tool: &mut Tool, active_flyout: &mut ActiveFlyout| {
        let mut keep = true;
        egui::Window::new(title)
            .open(&mut keep)
            .resizable(false)
            .collapsible(false)
            .default_pos(Pos2::new(72.0, y_pos))
            .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
            .show(ctx, |ui| {
                ui.horizontal_wrapped(|ui| {
                    for &(t, tip) in items {
                        if tool_icon_button(ui, t, *tool == t, false, tip).clicked() {
                            *tool = t;
                        }
                    }
                });
            });
        if !keep {
            *active_flyout = ActiveFlyout::None;
        }
    };

    match *active_flyout {
        ActiveFlyout::Accessories => show_mini_flyout(
            "Accessories",
            120.0,
            &[
                (Tool::Table, "Table (2x2)"),
                (Tool::TlcPlate, "TLC Plate"),
                (Tool::ReactionMap, "Reaction Atom Map"),
                (Tool::AcyclicChain, "Acyclic Chain"),
                (Tool::Perspective, "Structure Perspective"),
                (Tool::Fragmentation, "Mass Fragmentation"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::MultipleBonds => show_mini_flyout(
            "Multiple Bonds",
            140.0,
            &[
                (Tool::DoubleBond, "Double Bond"),
                (Tool::TripleBond, "Triple Bond"),
                (Tool::DativeBond, "Dative Bond"),
                (Tool::QuadrupleBond, "Quadruple Bond"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::DottedDashed => show_mini_flyout(
            "Dotted/Dashed Bond Tools",
            165.0,
            &[
                (Tool::DashedBond, "Dashed Bond"),
                (Tool::HashedBond, "Hashed Bond"),
                (Tool::HashWedgeBond, "Hashed Wedged Bond"),
                (Tool::HollowWedgeBond, "Hollow Wedged Bond"),
                (Tool::WavyBond, "Wavy Bond"),
                (Tool::BoldBond, "Bold Bond"),
                (Tool::WedgeBond, "Wedged Bond"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::Arrows => show_mini_flyout(
            "Arrows",
            195.0,
            &[
                (Tool::Arrow, "Reaction Arrow"),
                (Tool::EquilibriumArrow, "Equilibrium Arrow"),
                (Tool::RetrosyntheticArrow, "Retrosynthetic Arrow"),
                (Tool::NoGoArrow, "No-Go (Cross) Arrow"),
                (Tool::CurvedArrow, "Curved Mechanism Arrow"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::Orbitals => show_mini_flyout(
            "Orbitals",
            225.0,
            &[
                (Tool::Orbital, "p Orbital"),
                (Tool::OrbitalS, "s Orbital"),
                (Tool::OrbitalD, "d Orbital"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::Shapes => show_mini_flyout(
            "Drawing Elements",
            250.0,
            &[
                (Tool::ShapeRect, "Shadow Rectangle"),
                (Tool::ShapeOval, "Oval / Circle"),
                (Tool::ShapeLine, "Straight Line"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::Brackets => show_mini_flyout(
            "Brackets",
            280.0,
            &[
                (Tool::Bracket, "Square Brackets [ ]"),
                (Tool::RoundBracket, "Parentheses ( )"),
                (Tool::CurlyBracket, "Curly Braces { }"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::ChemSymbols => show_mini_flyout(
            "Chemical Symbols",
            305.0,
            &[
                (Tool::ChemSymbol, "Circle Plus (+)"),
                (Tool::SymbolMinus, "Circle Minus (-)"),
                (Tool::SymbolRadical, "Radical Dot (•)"),
                (Tool::SymbolLonePair, "Lone Pair (:)"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::Chains => show_mini_flyout(
            "Chains",
            330.0,
            &[
                (Tool::AcyclicChain, "Acyclic Chain"),
                (Tool::SolidBond, "Solid Bond"),
                (Tool::Cyclopropane, "Cyclopropane"),
                (Tool::Cyclobutane, "Cyclobutane"),
                (Tool::Cyclopentane, "Cyclopentane"),
                (Tool::Cyclohexane, "Cyclohexane"),
                (Tool::Benzene, "Benzene"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::Chromatography => show_mini_flyout(
            "Chromatography Tools",
            335.0,
            &[
                (Tool::TlcPlate, "TLC Plate"),
                (Tool::Table, "Data Table"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::Fragmentation => show_mini_flyout(
            "Fragmentation Tools",
            120.0,
            &[
                (Tool::Fragmentation, "Mass Fragmentation (Homolytic Cleavage)"),
                (Tool::Eraser, "Eraser"),
                (Tool::Lasso, "Lasso"),
                (Tool::Marquee, "Marquee"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::PenTools => show_mini_flyout(
            "Pen Tools",
            175.0,
            &[
                (Tool::PenCurve, "Pen / Curve Tool"),
                (Tool::ShapeLine, "Straight Line"),
                (Tool::CurvedArrow, "Curved Mechanism Arrow"),
            ],
            tool,
            active_flyout,
        ),
        ActiveFlyout::Colors => {
            let mut keep = true;
            egui::Window::new("Colors")
                .open(&mut keep)
                .resizable(false)
                .collapsible(false)
                .default_pos(Pos2::new(72.0, 210.0))
                .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
                .show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        for (tip, col) in [
                            ("Foreground (Black)", Color32::BLACK),
                            ("Other # 1 (Red)", Color32::from_rgb(255, 0, 0)),
                            ("Other # 2 (Yellow)", Color32::from_rgb(240, 220, 0)),
                            ("Other # 3 (Green)", Color32::from_rgb(0, 210, 0)),
                            ("Other # 4 (Cyan)", Color32::from_rgb(0, 220, 235)),
                            ("Other # 5 (Blue)", Color32::from_rgb(0, 0, 255)),
                            ("Other # 6 (Magenta)", Color32::from_rgb(255, 0, 255)),
                        ] {
                            let (r, resp) = ui.allocate_exact_size(Vec2::splat(18.0), Sense::click());
                            ui.painter().rect_filled(r, 0.0, col);
                            ui.painter().rect_stroke(r, 0.0, Stroke::new(1.0_f32, Color32::GRAY), egui::StrokeKind::Inside);
                            if resp.on_hover_text(tip).clicked() {
                                cmd = PaletteCommand::SetColor(col);
                            }
                        }
                    });
                    ui.horizontal(|ui| {
                        if ui.small_button("By Element").clicked() { cmd = PaletteCommand::ToggleColorByElement; }
                        if ui.small_button("Other...").clicked() { cmd = PaletteCommand::OpenColorPicker; }
                    });
                });
            if !keep { *active_flyout = ActiveFlyout::None; }
        }
        ActiveFlyout::Curves => {
            let mut keep = true;
            egui::Window::new("Curves")
                .open(&mut keep)
                .resizable(false)
                .collapsible(false)
                .default_pos(Pos2::new(72.0, 240.0))
                .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
                .show(ctx, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        for (lbl, key) in [
                            ("Plain", "plain"), ("Dashed", "dashed"), ("Bold", "bold"), ("Doubled", "doubled"),
                            ("Arrow End", "arrow_end"), ("Arrow Start", "arrow_start"),
                            ("Closed", "closed"), ("Filled", "filled"), ("Shaded", "shaded"), ("Faded", "faded"),
                        ] {
                            if ui.small_button(lbl).clicked() { cmd = PaletteCommand::ToggleCurve(key); }
                        }
                    });
                });
            if !keep { *active_flyout = ActiveFlyout::None; }
        }
        ActiveFlyout::Windows => {
            let mut keep = true;
            egui::Window::new("Windows")
                .open(&mut keep)
                .resizable(false)
                .collapsible(false)
                .default_pos(Pos2::new(72.0, 270.0))
                .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
                .show(ctx, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        for (lbl, key) in [
                            ("Cascade", "cascade"), ("Tile V", "tile_v"), ("Tile H", "tile_h"),
                            ("Analysis", "analysis"), ("Properties", "properties"),
                            ("Periodic Table", "ptable"), ("Info", "info"), ("Char Map", "charmap"),
                        ] {
                            if ui.small_button(lbl).clicked() { cmd = PaletteCommand::WindowAction(key); }
                        }
                    });
                });
            if !keep { *active_flyout = ActiveFlyout::None; }
        }
        ActiveFlyout::Templates => {
            if !ctp_libraries.is_empty() {
                if *selected_ctp_idx >= ctp_libraries.len() {
                    *selected_ctp_idx = 0;
                }
                let mut open = true;
                egui::Window::new("Templates")
                    .open(&mut open)
                    .default_pos(Pos2::new(72.0, 220.0))
                    .frame(egui::Frame::window(&ctx.style()).fill(CD_WINDOW_BG))
                    .show(ctx, |ui| {
                        let cur_name = ctp_libraries[*selected_ctp_idx].name.clone();
                        egui::ComboBox::from_id_salt("ctp_flyout_combo")
                            .width(220.0)
                            .selected_text(cur_name)
                            .show_ui(ui, |ui| {
                                for (i, lib) in ctp_libraries.iter().enumerate() {
                                    ui.selectable_value(
                                        selected_ctp_idx,
                                        i,
                                        format!("{} ({})", lib.name, lib.total_count()),
                                    );
                                }
                            });
                        ui.separator();
                        let lib = &ctp_libraries[*selected_ctp_idx];
                        egui::ScrollArea::vertical().max_height(260.0).show(ui, |ui| {
                            egui::Grid::new("tpl_grid").num_columns(4).spacing([4.0, 4.0]).show(ui, |ui| {
                                if !lib.templates.is_empty() {
                                    for (ti, mol) in lib.templates.iter().enumerate() {
                                        if draw_molecule_thumbnail(ui, mol, ti).clicked() {
                                            cmd = PaletteCommand::StampTemplate(*selected_ctp_idx, ti);
                                        }
                                        if (ti + 1) % 4 == 0 { ui.end_row(); }
                                    }
                                } else if !lib.graphic_templates.is_empty() {
                                    for (ti, curves) in lib.graphic_templates.iter().enumerate() {
                                        if draw_graphic_thumbnail(ui, curves, ti).clicked() {
                                            cmd = PaletteCommand::StampTemplate(*selected_ctp_idx, ti);
                                        }
                                        if (ti + 1) % 4 == 0 { ui.end_row(); }
                                    }
                                } else {
                                    ui.label("Empty library — draw a structure and open New Templates to populate.");
                                }
                            });
                        });
                    });
                if !open {
                    *active_flyout = ActiveFlyout::None;
                }
            }
        }
        ActiveFlyout::None => {}
    }

    cmd
}
