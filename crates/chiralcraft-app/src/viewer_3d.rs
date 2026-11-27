//! Interactive Pure-Rust 3D Molecular & Crystal Viewer (`egui` Software Renderer)
//!
//! - Drag-to-orbit 3D rotation (`yaw`, `pitch`), zoom, and `a`/`b`/`c` axis snap
//! - Render modes: **Ball & Stick**, **Stick**, **Space-Filling (CPK)**, **Wireframe**
//! - CPK crystallographic element colors with radial specular shading & depth-sorted Painter's algorithm
//! - Optional crystallographic **Unit Cell wireframe** for `.cif` and VASP `POSCAR`/`CONTCAR`
//! - **"Snapshot to 2D Canvas"** stamps the exact rotated 3D projection onto the 2D drawing sheet
//!   preserving `z_coord` for Rear Opacity depth cueing and crossing-bond bridge gaps!

use chiralcraft_chem::{
    get_element, measure_angle_deg, measure_dihedral_deg, measure_distance_angstrom,
    optimize_geometry_3d, ChemMolecule,
};
use chiralcraft_doc::lift_molecule_2d_to_3d;
use chiralcraft_formats::CrystalUnitCell;
use egui::{Color32, FontId, Pos2, Sense, Stroke, Vec2};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RenderMode3d {
    BallAndStick,
    Stick,
    SpaceFilling,
    Wireframe,
}

#[derive(Debug, Clone)]
pub struct VibrationalAnimState {
    pub base_mol_3d: ChemMolecule,
    pub wavenumber_cm1: f64,
    pub label: String,
    pub mode_vectors: Vec<[f64; 3]>,
    pub amplitude_pt: f64,
    pub speed_hz: f64,
    pub playing: bool,
    pub show_vectors: bool,
}

pub struct Viewer3dState {
    pub open: bool,
    pub yaw: f64,
    pub pitch: f64,
    pub zoom: f32,
    pub mode: RenderMode3d,
    pub show_unit_cell: bool,
    pub show_labels: bool,
    pub dark_bg: bool,
    pub unit_cell: Option<CrystalUnitCell>,
    pub measure_atoms: Vec<usize>,
    pub last_opt_summary: String,
    pub vib_anim: Option<VibrationalAnimState>,
}

impl Default for Viewer3dState {
    fn default() -> Self {
        Self {
            open: false,
            yaw: 0.45,
            pitch: 0.30,
            zoom: 1.0,
            mode: RenderMode3d::BallAndStick,
            show_unit_cell: true,
            show_labels: false,
            dark_bg: true,
            unit_cell: None,
            measure_atoms: Vec::new(),
            last_opt_summary: String::new(),
            vib_anim: None,
        }
    }
}

pub enum Viewer3dAction {
    None,
    SnapshotTo2d(ChemMolecule),
    LiftActiveTo3d,
    ApplyOptimized3d(ChemMolecule, String),
}

/// CPK crystallographic element RGB color for 3D spheres & sticks.
pub fn vesta_element_color(z: i32) -> Color32 {
    match z {
        1 => Color32::from_rgb(238, 238, 238),  // H
        5 => Color32::from_rgb(255, 181, 181),  // B
        6 => Color32::from_rgb(90, 96, 106),    // C
        7 => Color32::from_rgb(48, 80, 248),    // N
        8 => Color32::from_rgb(235, 38, 38),    // O
        9 => Color32::from_rgb(144, 224, 80),   // F
        14 => Color32::from_rgb(218, 175, 128), // Si
        15 => Color32::from_rgb(255, 128, 0),   // P
        16 => Color32::from_rgb(235, 215, 45),  // S
        17 => Color32::from_rgb(31, 220, 31),   // Cl
        26 => Color32::from_rgb(224, 102, 51),  // Fe
        28 => Color32::from_rgb(80, 208, 80),   // Ni
        29 => Color32::from_rgb(200, 128, 51),  // Cu
        30 => Color32::from_rgb(125, 128, 176), // Zn
        35 => Color32::from_rgb(166, 41, 41),   // Br
        46 => Color32::from_rgb(0, 105, 133),   // Pd
        53 => Color32::from_rgb(148, 0, 148),   // I
        78 => Color32::from_rgb(185, 185, 205), // Pt
        79 => Color32::from_rgb(245, 205, 35),  // Au
        _ => Color32::from_rgb(155, 130, 195),
    }
}

fn vdw_radius_factor(z: i32) -> f32 {
    match z {
        1 => 0.65,
        6 => 1.00,
        7 => 0.92,
        8 => 0.88,
        9 => 0.84,
        14 => 1.22,
        15 => 1.12,
        16 => 1.10,
        17 => 1.05,
        35 => 1.15,
        53 => 1.28,
        _ => 1.18,
    }
}

/// Rotates a centered 3D point `(x, y, z)` by `(yaw, pitch)`.
pub fn rotate_point_3d(x: f64, y: f64, z: f64, yaw: f64, pitch: f64) -> (f64, f64, f64) {
    let (sy, cy) = yaw.sin_cos();
    let (sp, cp) = pitch.sin_cos();
    let x1 = x * cy + z * sy;
    let z1 = -x * sy + z * cy;
    let y2 = y * cp - z1 * sp;
    let z2 = y * sp + z1 * cp;
    (x1, y2, z2)
}

/// Produces a rotated 3D `ChemMolecule` snapshot centered at `(220, 170)` ready to stamp onto the 2D canvas.
pub fn snapshot_rotated_molecule(mol: &ChemMolecule, yaw: f64, pitch: f64) -> ChemMolecule {
    let lifted = lift_molecule_2d_to_3d(mol);
    if lifted.atoms.is_empty() {
        return lifted;
    }
    let nf = lifted.atoms.len() as f64;
    let cx = lifted.atoms.iter().map(|a| a.x).sum::<f64>() / nf;
    let cy = lifted.atoms.iter().map(|a| a.y).sum::<f64>() / nf;
    let cz = lifted.atoms.iter().map(|a| a.z_coord).sum::<f64>() / nf;
    let mut out = lifted;
    for a in &mut out.atoms {
        let (rx, ry, rz) = rotate_point_3d(a.x - cx, a.y - cy, a.z_coord - cz, yaw, pitch);
        a.x = 220.0 + rx;
        a.y = 170.0 + ry;
        a.z_coord = rz;
    }
    out.has_3d_coords = true;
    out
}

pub fn show_3d_viewer_window(
    ctx: &egui::Context,
    state: &mut Viewer3dState,
    active_mol: Option<&ChemMolecule>,
) -> Viewer3dAction {
    if !state.open {
        return Viewer3dAction::None;
    }
    let mut action = Viewer3dAction::None;
    let mut open = state.open;

    egui::Window::new("3D Molecular & Crystal Viewer")
        .open(&mut open)
        .default_pos([340.0, 74.0])
        .default_size([450.0, 415.0])
        .constrain_to(ctx.available_rect())
        .resizable(true)
        .show(ctx, |ui| {
            // Top controls
            ui.horizontal_wrapped(|ui| {
                for (lbl, m) in [
                    ("Ball & Stick", RenderMode3d::BallAndStick),
                    ("Stick", RenderMode3d::Stick),
                    ("Space-Fill", RenderMode3d::SpaceFilling),
                    ("Wire", RenderMode3d::Wireframe),
                ] {
                    if ui.selectable_label(state.mode == m, lbl).clicked() {
                        state.mode = m;
                    }
                }
                ui.separator();
                if ui.small_button("a").on_hover_text("Align along a-axis").clicked() {
                    state.yaw = std::f64::consts::FRAC_PI_2;
                    state.pitch = 0.0;
                }
                if ui.small_button("b").on_hover_text("Align along b-axis").clicked() {
                    state.yaw = 0.0;
                    state.pitch = std::f64::consts::FRAC_PI_2;
                }
                if ui.small_button("c").on_hover_text("Align along c-axis").clicked() {
                    state.yaw = 0.0;
                    state.pitch = 0.0;
                }
                if ui.small_button("Iso").on_hover_text("Reset isometric view").clicked() {
                    state.yaw = 0.45;
                    state.pitch = 0.30;
                    state.zoom = 1.0;
                }
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut state.show_labels, "Atom Symbols");
                ui.checkbox(&mut state.dark_bg, "Dark Stage");
                if state.unit_cell.is_some() {
                    ui.checkbox(&mut state.show_unit_cell, "Unit Cell");
                }
                ui.add(egui::Slider::new(&mut state.zoom, 0.35..=3.5).text("Zoom"));
            });

            // Optional Vibrational Normal-Mode Animation Controls (from xTB / ORCA / Gaussian / Molden IR)
            let mut stop_vib = false;
            if let Some(anim) = &mut state.vib_anim {
                ui.horizontal_wrapped(|ui| {
                    let lbl = if anim.playing { "Pause" } else { "Play" };
                    if ui.small_button(lbl).clicked() {
                        anim.playing = !anim.playing;
                    }
                    ui.checkbox(&mut anim.show_vectors, "Vectors");
                    ui.add(
                        egui::Slider::new(&mut anim.amplitude_pt, 1.5..=14.0)
                            .fixed_decimals(1)
                            .text("Amp"),
                    );
                    ui.colored_label(
                        Color32::from_rgb(50, 195, 255),
                        format!("v = {:.1} cm^-1 ({})", anim.wavenumber_cm1, anim.label),
                    );
                    if ui.small_button("× Stop").clicked() {
                        stop_vib = true;
                    }
                });
                if anim.playing {
                    ctx.request_repaint();
                }
            }
            if stop_vib {
                state.vib_anim = None;
            }

            let avail = Vec2::new(ui.available_width().max(260.0), (ui.available_height() - 36.0).max(200.0));
            let (rect, resp) = ui.allocate_exact_size(avail, Sense::click_and_drag());
            if resp.dragged() {
                let d = resp.drag_delta();
                state.yaw += (d.x as f64) * 0.012;
                state.pitch = (state.pitch + (d.y as f64) * 0.012).clamp(-1.55, 1.55);
            }
            if resp.hovered() {
                let scroll = ui.input(|i| i.raw_scroll_delta.y);
                if scroll.abs() > 0.1 {
                    state.zoom = (state.zoom * (1.0 + scroll * 0.002)).clamp(0.35, 3.5);
                }
            }

            let p = ui.painter_at(rect);
            let bg = if state.dark_bg {
                Color32::from_rgb(20, 24, 32)
            } else {
                Color32::from_rgb(246, 248, 252)
            };
            p.rect_filled(rect, 4.0, bg);
            p.rect_stroke(rect, 4.0, Stroke::new(1.0_f32, Color32::from_gray(90)), egui::StrokeKind::Inside);

            let anim_phase = ui.input(|i| i.time);
            let display_mol_owned: Option<ChemMolecule> = if let Some(anim) = &state.vib_anim {
                let mut m = anim.base_mol_3d.clone();
                let s = if anim.playing {
                    (anim_phase * std::f64::consts::TAU * anim.speed_hz).sin()
                } else {
                    0.85
                };
                for (idx, a) in m.atoms.iter_mut().enumerate() {
                    if let Some(v) = anim.mode_vectors.get(idx) {
                        a.x += s * anim.amplitude_pt * v[0];
                        a.y -= s * anim.amplitude_pt * v[1];
                        a.z_coord += s * anim.amplitude_pt * v[2];
                    }
                }
                Some(m)
            } else {
                active_mol.filter(|m| !m.atoms.is_empty()).map(lift_molecule_2d_to_3d)
            };

            if let Some(mol) = display_mol_owned {
                if !mol.atoms.is_empty() {
                    let nf = mol.atoms.len() as f64;
                    let cx = mol.atoms.iter().map(|a| a.x).sum::<f64>() / nf;
                    let cy = mol.atoms.iter().map(|a| a.y).sum::<f64>() / nf;
                    let cz = mol.atoms.iter().map(|a| a.z_coord).sum::<f64>() / nf;

                    let mut max_r = 12.0_f64;
                    let mut rot_pts = Vec::with_capacity(mol.atoms.len());
                    for a in &mol.atoms {
                        let (rx, ry, rz) = rotate_point_3d(a.x - cx, a.y - cy, a.z_coord - cz, state.yaw, state.pitch);
                        max_r = max_r.max((rx * rx + ry * ry + rz * rz).sqrt());
                        rot_pts.push((rx, ry, rz));
                    }
                    let base_scale = (rect.width().min(rect.height()) * 0.36 / (max_r as f32)) * state.zoom;
                    let center = rect.center();
                    let project = |(rx, ry, rz): (f64, f64, f64)| -> (Pos2, f32) {
                        let persp = 1.0 + (rz as f32) / ((max_r as f32) * 5.0 + 60.0);
                        let sx = center.x + (rx as f32) * base_scale * persp;
                        let sy = center.y + (ry as f32) * base_scale * persp;
                        (Pos2::new(sx, sy), persp)
                    };

                    // Optional crystallographic unit cell wireframe
                    if state.show_unit_cell {
                        if let Some(cell) = &state.unit_cell {
                            let corners: Vec<Pos2> = cell
                                .corners()
                                .iter()
                                .map(|&(x, y, z)| {
                                    let r = rotate_point_3d(x - cx, y - cy, z - cz, state.yaw, state.pitch);
                                    project(r).0
                                })
                                .collect();
                            let cell_col = if state.dark_bg {
                                Color32::from_rgba_unmultiplied(120, 195, 255, 150)
                            } else {
                                Color32::from_rgba_unmultiplied(35, 95, 175, 160)
                            };
                            for (u, v) in CrystalUnitCell::edge_indices() {
                                p.line_segment([corners[u], corners[v]], Stroke::new(1.1_f32, cell_col));
                            }
                        }
                    }

                    // Depth-sorted draw list (bonds + atoms) via Painter's algorithm
                    enum DrawItem {
                        Bond { z_avg: f64, a: usize, b: usize, order: i32 },
                        Atom { z: f64, idx: usize },
                    }
                    let mut items = Vec::with_capacity(mol.atoms.len() + mol.bonds.len());
                    for b in &mol.bonds {
                        if b.begin_atom < rot_pts.len() && b.end_atom < rot_pts.len() {
                            let z_avg = 0.5 * (rot_pts[b.begin_atom].2 + rot_pts[b.end_atom].2);
                            items.push(DrawItem::Bond { z_avg, a: b.begin_atom, b: b.end_atom, order: b.order });
                        }
                    }
                    for (idx, pt) in rot_pts.iter().enumerate() {
                        items.push(DrawItem::Atom { z: pt.2, idx });
                    }
                    items.sort_by(|x, y| {
                        let zx = match x { DrawItem::Bond { z_avg, .. } => *z_avg, DrawItem::Atom { z, .. } => *z };
                        let zy = match y { DrawItem::Bond { z_avg, .. } => *z_avg, DrawItem::Atom { z, .. } => *z };
                        zx.partial_cmp(&zy).unwrap_or(std::cmp::Ordering::Equal)
                    });

                    for item in items {
                        match item {
                            DrawItem::Bond { a, b, order, .. } => {
                                if state.mode == RenderMode3d::SpaceFilling {
                                    continue;
                                }
                                let (p1, s1) = project(rot_pts[a]);
                                let (p2, s2) = project(rot_pts[b]);
                                let mid = p1 + (p2 - p1) * 0.5;
                                let c1 = vesta_element_color(mol.atoms[a].z);
                                let c2 = vesta_element_color(mol.atoms[b].z);
                                let base_w = match state.mode {
                                    RenderMode3d::Wireframe => 1.6_f32,
                                    RenderMode3d::Stick => 4.2_f32 * ((s1 + s2) * 0.5) * state.zoom.sqrt(),
                                    RenderMode3d::BallAndStick => 3.0_f32 * ((s1 + s2) * 0.5) * state.zoom.sqrt(),
                                    RenderMode3d::SpaceFilling => 0.0,
                                };
                                let d = p2 - p1;
                                let len = d.length().max(1e-3);
                                let n = Vec2::new(-d.y / len, d.x / len);
                                let offsets: &[f32] = if order == 2 {
                                    &[-base_w * 0.55, base_w * 0.55]
                                } else if order >= 3 {
                                    &[-base_w * 0.85, 0.0, base_w * 0.85]
                                } else {
                                    &[0.0]
                                };
                                for &off in offsets {
                                    let shift = n * off;
                                    p.line_segment([p1 + shift, mid + shift], Stroke::new(base_w.max(1.2), c1));
                                    p.line_segment([mid + shift, p2 + shift], Stroke::new(base_w.max(1.2), c2));
                                }
                            }
                            DrawItem::Atom { idx, .. } => {
                                if state.mode == RenderMode3d::Wireframe {
                                    continue;
                                }
                                let atom = &mol.atoms[idx];
                                let (pos, persp) = project(rot_pts[idx]);
                                let base_col = vesta_element_color(atom.z);
                                let rad = match state.mode {
                                    RenderMode3d::BallAndStick => 7.5 * vdw_radius_factor(atom.z) * persp * state.zoom.sqrt(),
                                    RenderMode3d::Stick => 2.5 * persp * state.zoom.sqrt(),
                                    RenderMode3d::SpaceFilling => 19.0 * vdw_radius_factor(atom.z) * persp * state.zoom.sqrt(),
                                    RenderMode3d::Wireframe => 0.0,
                                };
                                // Shaded sphere: dark rim + base fill + bright specular highlight
                                let rim = Color32::from_rgb(
                                    (base_col.r() as f32 * 0.45) as u8,
                                    (base_col.g() as f32 * 0.45) as u8,
                                    (base_col.b() as f32 * 0.45) as u8,
                                );
                                p.circle_filled(pos, rad, base_col);
                                p.circle_stroke(pos, rad, Stroke::new(0.8_f32, rim));
                                if rad > 3.5 {
                                    let hi = Color32::from_rgba_unmultiplied(255, 255, 255, 125);
                                    p.circle_filled(pos + Vec2::new(-rad * 0.30, -rad * 0.30), rad * 0.32, hi);
                                }
                                if state.show_labels && rad > 4.0 {
                                    let sym = get_element(atom.z).map(|e| e.symbol).unwrap_or("?");
                                    let lcol = if atom.z == 1 || atom.z == 16 || atom.z == 79 {
                                        Color32::BLACK
                                    } else {
                                        Color32::WHITE
                                    };
                                    p.text(pos, egui::Align2::CENTER_CENTER, sym, FontId::proportional((rad * 0.95).clamp(7.5, 13.0)), lcol);
                                }
                            }
                        }
                    }

                    // Draw vibrational displacement vectors when animating a quantum normal mode
                    if let Some(anim) = &state.vib_anim {
                        if anim.show_vectors {
                            let vec_col = Color32::from_rgb(255, 205, 45);
                            for (idx, pt) in rot_pts.iter().enumerate() {
                                if let Some(v) = anim.mode_vectors.get(idx) {
                                    let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                                    if mag < 0.04 {
                                        continue;
                                    }
                                    let tip_raw = rotate_point_3d(
                                        (mol.atoms[idx].x - cx) + v[0] * anim.amplitude_pt * 1.8,
                                        (mol.atoms[idx].y - cy) - v[1] * anim.amplitude_pt * 1.8,
                                        (mol.atoms[idx].z_coord - cz) + v[2] * anim.amplitude_pt * 1.8,
                                        state.yaw,
                                        state.pitch,
                                    );
                                    let p0 = project(*pt).0;
                                    let p1 = project(tip_raw).0;
                                    p.line_segment([p0, p1], Stroke::new(1.8_f32, vec_col));
                                    p.circle_filled(p1, 2.6, vec_col);
                                }
                            }
                        }
                    }
                    // Click-to-measure: pick nearest projected atom on primary click!
                    state.measure_atoms.retain(|&i| i < mol.atoms.len());
                    if resp.clicked() {
                        if let Some(cpos) = resp.interact_pointer_pos() {
                            let mut best: Option<(usize, f32)> = None;
                            for (idx, pt) in rot_pts.iter().enumerate() {
                                let sp = project(*pt).0;
                                let d = sp.distance(cpos);
                                if d < 16.0 && best.map(|(_, bd)| d < bd).unwrap_or(true) {
                                    best = Some((idx, d));
                                }
                            }
                            if let Some((hit_i, _)) = best {
                                if let Some(pos_i) = state.measure_atoms.iter().position(|&x| x == hit_i) {
                                    state.measure_atoms.remove(pos_i);
                                } else {
                                    if state.measure_atoms.len() >= 4 {
                                        state.measure_atoms.clear();
                                    }
                                    state.measure_atoms.push(hit_i);
                                }
                            }
                        }
                    }

                    // Highlight measured atoms in bright cyan rings
                    for (ord_k, &mai) in state.measure_atoms.iter().enumerate() {
                        if let Some(&pt) = rot_pts.get(mai) {
                            let sp = project(pt).0;
                            p.circle_stroke(sp, 11.0, Stroke::new(2.0_f32, Color32::from_rgb(50, 230, 255)));
                            p.text(
                                sp + Vec2::new(11.0, -11.0),
                                egui::Align2::LEFT_BOTTOM,
                                format!("M{}", ord_k + 1),
                                FontId::proportional(9.5),
                                Color32::from_rgb(50, 230, 255),
                            );
                        }
                    }
                    // Live 3D measurement readout (Distance / Angle / Dihedral)
                    let meas_txt = match state.measure_atoms.as_slice() {
                        [a, b] => measure_distance_angstrom(&mol, *a, *b)
                            .map(|d| format!("Distance #{}-#{}: {d:.3} Å", a + 1, b + 1)),
                        [a, b, c] => measure_angle_deg(&mol, *a, *b, *c)
                            .map(|ang| format!("Angle #{}-#{}-#{}: {ang:.1}°", a + 1, b + 1, c + 1)),
                        [a, b, c, d] => measure_dihedral_deg(&mol, *a, *b, *c, *d)
                            .map(|dih| format!("Dihedral #{}-#{}-#{}-#{}: {dih:.1}°", a + 1, b + 1, c + 1, d + 1)),
                        _ => None,
                    };
                    if let Some(mt) = meas_txt {
                        p.text(
                            Pos2::new(rect.min.x + 8.0, rect.min.y + 8.0),
                            egui::Align2::LEFT_TOP,
                            mt,
                            FontId::proportional(11.0),
                            Color32::from_rgb(80, 235, 255),
                        );
                    }

                    // Subtle footer hint inside viewport
                    let hint = if let Some(cell) = &state.unit_cell {
                        format!("Drag: orbit • Click 2–4 atoms: measure • Cell: a={:.2} b={:.2} c={:.2} Å", cell.a, cell.b, cell.c)
                    } else {
                        format!("Drag: orbit • Click 2–4 atoms: measure • {} atoms, {} bonds", mol.atoms.len(), mol.bonds.len())
                    };
                    p.text(
                        Pos2::new(rect.min.x + 8.0, rect.max.y - 6.0),
                        egui::Align2::LEFT_BOTTOM,
                        hint,
                        FontId::proportional(10.0),
                        if state.dark_bg { Color32::from_gray(175) } else { Color32::from_gray(90) },
                    );
                } else {
                    p.text(rect.center(), egui::Align2::CENTER_CENTER, "Draw or open a molecule to inspect in 3D", FontId::proportional(12.0), Color32::GRAY);
                }
            } else {
                p.text(rect.center(), egui::Align2::CENTER_CENTER, "Draw or open a molecule (.mol, .xyz, .cif, .pdb, POSCAR) to inspect in 3D", FontId::proportional(12.0), Color32::GRAY);
            }

            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                let has_mol = active_mol.map(|m| !m.atoms.is_empty()).unwrap_or(false);
                if ui.add_enabled(has_mol, egui::Button::new("Optimize 3D (UFF)")).clicked() {
                    if let Some(m) = active_mol {
                        let lifted = lift_molecule_2d_to_3d(m);
                        let (opt, rep) = optimize_geometry_3d(&lifted, 180);
                        let summary = format!(
                            "UFF 3D Opt ({} iters): E {:.1} -> {:.2} kcal/mol",
                            rep.iterations, rep.initial_energy_kcal, rep.final_energy_kcal
                        );
                        state.last_opt_summary = summary.clone();
                        action = Viewer3dAction::ApplyOptimized3d(opt, summary);
                    }
                }
                if ui.add_enabled(has_mol, egui::Button::new("Snapshot to 2D Canvas")).clicked() {
                    if let Some(m) = active_mol {
                        action = Viewer3dAction::SnapshotTo2d(snapshot_rotated_molecule(m, state.yaw, state.pitch));
                    }
                }
                if ui.add_enabled(has_mol, egui::Button::new("2D -> 3D Lift")).clicked() {
                    action = Viewer3dAction::LiftActiveTo3d;
                }
                if !state.measure_atoms.is_empty() && ui.small_button("Clear Measure").clicked() {
                    state.measure_atoms.clear();
                }
            });
            if !state.last_opt_summary.is_empty() {
                ui.small(&state.last_opt_summary);
            }
        });

    state.open = open;
    action
}

#[cfg(test)]
mod tests {
    use super::*;
    use chiralcraft_formats::smiles_to_molecule;

    #[test]
    fn test_3d_viewer_snapshot_preserves_atoms_and_3d_coords() {
        let m = smiles_to_molecule("C1CCCCC1").expect("cyclohexane");
        let snap = snapshot_rotated_molecule(&m, 0.6, 0.4);
        assert_eq!(snap.atoms.len(), m.atoms.len());
        assert!(snap.has_3d_coords);
        assert!(snap.atoms.iter().any(|a| a.z_coord.abs() > 0.5));
    }
}
