//! ChiralCraft Scene Models:
//! - Enhanced Stereochemistry (`Abs`, `Or(n)`, `And(n)`)
//! - Live Reaction Stoichiometry Grid (`MW`, `Eq.`, `mmol`, `Mass`, `Density`, `Volume`, `% Yield`)
//! - Gel Electrophoresis Plate (`GEP`) with MW markers & migration bands
//! - Ring Fill & Atom Highlight colors + Reaction AutoNumber labels

use std::collections::HashMap;

use chiralcraft_chem::{compute_analysis_report, BondDisplay, ChemAtom, ChemBond, ChemMolecule};

use crate::{DocumentScene, SceneArrow};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnhancedStereoKind {
    Absolute,
    Or(u32),
    And(u32),
}

impl EnhancedStereoKind {
    pub fn display_label(self) -> String {
        match self {
            Self::Absolute => "abs".to_string(),
            Self::Or(n) => format!("or{n}"),
            Self::And(n) => format!("&{n}"),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StoichiometryRow {
    pub mol_idx: usize,
    pub role: String, // "Reactant" | "Product"
    pub label: String,
    pub formula: String,
    pub mw: f64,
    pub equivalents: f64,
    pub mmol: f64,
    pub mass_mg: f64,
    pub density_g_ml: f64,
    pub volume_ul: f64,
    pub yield_pct: f64,
}

impl StoichiometryRow {
    /// Recomputes `mmol`, `mass_mg`, and `volume_ul` from `base_mmol` and `equivalents` / `yield_pct`.
    pub fn recalculate(&mut self, base_mmol: f64) {
        if self.role == "Reactant" {
            self.mmol = base_mmol * self.equivalents;
            self.mass_mg = self.mmol * self.mw;
        } else {
            self.mmol = base_mmol * self.equivalents * (self.yield_pct / 100.0);
            self.mass_mg = self.mmol * self.mw;
        }
        self.volume_ul = if self.density_g_ml > 0.01 {
            self.mass_mg / self.density_g_ml
        } else {
            0.0
        };
    }
}

/// Builds a live Reaction Stoichiometry Grid from the molecules & reaction arrow in `scene`.
pub fn build_stoichiometry_rows(scene: &DocumentScene, base_mmol: f64) -> Vec<StoichiometryRow> {
    let split_x = scene
        .arrows
        .first()
        .map(|a| 0.5 * (a.tail.0 + a.head.0))
        .unwrap_or(f64::INFINITY);
    let mut rows = Vec::new();
    let mut r_idx = 1usize;
    let mut p_idx = 1usize;
    for (mi, mol) in scene.molecules.iter().enumerate() {
        if mol.atoms.is_empty() {
            continue;
        }
        let cx = mol.atoms.iter().map(|a| a.x).sum::<f64>() / (mol.atoms.len() as f64);
        let is_product = split_x.is_finite() && cx > split_x;
        let rep = compute_analysis_report(mol);
        let (role, label) = if is_product {
            let l = format!("Product {p_idx}");
            p_idx += 1;
            ("Product".to_string(), l)
        } else {
            let l = format!("Reactant {r_idx}");
            r_idx += 1;
            ("Reactant".to_string(), l)
        };
        let mut row = StoichiometryRow {
            mol_idx: mi,
            role,
            label,
            formula: rep.formula,
            mw: rep.mol_weight,
            equivalents: 1.0,
            mmol: base_mmol,
            mass_mg: base_mmol * rep.mol_weight,
            density_g_ml: 1.0,
            volume_ul: base_mmol * rep.mol_weight,
            yield_pct: 100.0,
        };
        row.recalculate(base_mmol);
        rows.push(row);
    }
    rows
}

/// Formats the Stoichiometry Grid as a clean text table for pasting onto the canvas or CSV export.
pub fn format_stoichiometry_text(rows: &[StoichiometryRow], csv: bool) -> String {
    let sep = if csv { "," } else { " | " };
    let mut out = format!(
        "Compound{sep}Role{sep}Formula{sep}MW{sep}Eq.{sep}mmol{sep}Mass (mg){sep}Yield (%)\n"
    );
    for r in rows {
        out.push_str(&format!(
            "{}{sep}{}{sep}{}{sep}{:.2}{sep}{:.2}{sep}{:.3}{sep}{:.1}{sep}{:.1}\n",
            r.label, r.role, r.formula, r.mw, r.equivalents, r.mmol, r.mass_mg, r.yield_pct
        ));
    }
    out
}

#[derive(Debug, Clone, PartialEq)]
pub struct GelLane {
    pub name: String,
    /// Relative migration positions `0.05..=0.95` from top wells downward
    pub bands: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SceneGelPlate {
    pub id: u32,
    pub bbox: (f64, f64, f64, f64),
    pub lanes: Vec<GelLane>,
}

impl SceneGelPlate {
    pub fn default_at(id: u32, x: f64, y: f64) -> Self {
        Self {
            id,
            bbox: (x, y, x + 96.0, y + 84.0),
            lanes: vec![
                GelLane { name: "M".to_string(), bands: vec![0.18, 0.32, 0.48, 0.64, 0.80] },
                GelLane { name: "1".to_string(), bands: vec![0.34, 0.72] },
                GelLane { name: "2".to_string(), bands: vec![0.46, 0.68] },
                GelLane { name: "3".to_string(), bands: vec![0.28, 0.52, 0.78] },
            ],
        }
    }
}

/// Finds all 3..=8-membered ring vertex polygons in `mol` for `Colors > Ring Fill Color`.
pub fn extract_ring_polygons(mol: &ChemMolecule) -> Vec<Vec<(f64, f64)>> {
    let rings = crate::find_chordless_small_rings(mol, 8);
    let mut polys = Vec::new();
    for (verts, _) in rings {
        if verts.len() >= 3 {
            let pts: Vec<(f64, f64)> = verts
                .iter()
                .filter_map(|&vi| mol.atoms.get(vi).map(|a| (a.x, a.y)))
                .collect();
            if pts.len() >= 3 {
                polys.push(pts);
            }
        }
    }
    polys
}

// ============================================================================
// Crossing-Bond Clearance (Bridge Gaps) & 3D Rear Opacity Helpers
// ============================================================================

/// Detects all strict segment-segment intersections between non-adjacent bonds in `mol`.
/// Returns a map from the **rear** bond index to the intersection coordinates `(ix, iy)`
/// where a front bond crosses over it, so renderers can cleanly split the rear bond with a bridge gap!
pub fn find_bond_crossings(mol: &ChemMolecule) -> HashMap<usize, Vec<(f64, f64)>> {
    let mut rear_gaps: HashMap<usize, Vec<(f64, f64)>> = HashMap::new();
    let nb = mol.bonds.len();
    for i in 0..nb {
        let bi = &mol.bonds[i];
        if !bi.visible || bi.begin_atom >= mol.atoms.len() || bi.end_atom >= mol.atoms.len() {
            continue;
        }
        let a1 = &mol.atoms[bi.begin_atom];
        let a2 = &mol.atoms[bi.end_atom];
        let (x1, y1, x2, y2) = (a1.x, a1.y, a2.x, a2.y);
        for j in (i + 1)..nb {
            let bj = &mol.bonds[j];
            if !bj.visible || bj.begin_atom >= mol.atoms.len() || bj.end_atom >= mol.atoms.len() {
                continue;
            }
            // Skip adjacent bonds sharing a vertex
            if bi.begin_atom == bj.begin_atom
                || bi.begin_atom == bj.end_atom
                || bi.end_atom == bj.begin_atom
                || bi.end_atom == bj.end_atom
            {
                continue;
            }
            let b1 = &mol.atoms[bj.begin_atom];
            let b2 = &mol.atoms[bj.end_atom];
            let (x3, y3, x4, y4) = (b1.x, b1.y, b2.x, b2.y);
            let denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
            if denom.abs() < 1e-6 {
                continue;
            }
            let t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
            let u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom;
            if t > 0.14 && t < 0.86 && u > 0.14 && u < 0.86 {
                let ix = x1 + t * (x2 - x1);
                let iy = y1 + t * (y2 - y1);
                // Choose which bond is behind:
                // 1) If 3D coords exist, lower average Z is rear;
                // 2) Wedged/Bold bonds stay in front of plain/hashed bonds;
                // 3) Otherwise higher index is rear (deterministic).
                let zi = 0.5 * (a1.z_coord + a2.z_coord);
                let zj = 0.5 * (b1.z_coord + b2.z_coord);
                let is_front_style = |d: BondDisplay| {
                    matches!(d, BondDisplay::WedgeBegin | BondDisplay::WedgeEnd | BondDisplay::Bold)
                };
                let rear_idx = if mol.has_3d_coords && (zi - zj).abs() > 0.25 {
                    if zi < zj { i } else { j }
                } else if is_front_style(bi.display) && !is_front_style(bj.display) {
                    j
                } else if is_front_style(bj.display) && !is_front_style(bi.display) {
                    i
                } else {
                    j
                };
                rear_gaps.entry(rear_idx).or_default().push((ix, iy));
            }
        }
    }
    rear_gaps
}

/// Splits a segment `(x1, y1) -> (x2, y2)` into visible sub-segments leaving a clean `gap_radius`
/// clearance around every crossing point in `crossings`.
pub fn split_segment_for_crossings(
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
    crossings: &[(f64, f64)],
    gap_radius: f64,
) -> Vec<((f64, f64), (f64, f64))> {
    if crossings.is_empty() || gap_radius <= 0.1 {
        return vec![((x1, y1), (x2, y2))];
    }
    let dx = x2 - x1;
    let dy = y2 - y1;
    let len = dx.hypot(dy);
    if len < 1e-4 {
        return vec![((x1, y1), (x2, y2))];
    }
    let dt = (gap_radius / len).clamp(0.02, 0.35);
    let mut intervals: Vec<(f64, f64)> = crossings
        .iter()
        .map(|&(ix, iy)| {
            let t = (((ix - x1) * dx + (iy - y1) * dy) / (len * len)).clamp(0.0, 1.0);
            ((t - dt).max(0.0), (t + dt).min(1.0))
        })
        .collect();
    intervals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut out = Vec::new();
    let mut cur = 0.0_f64;
    for (g0, g1) in intervals {
        if g0 > cur + 0.02 {
            out.push(((x1 + dx * cur, y1 + dy * cur), (x1 + dx * g0, y1 + dy * g0)));
        }
        cur = cur.max(g1);
    }
    if cur < 0.98 {
        out.push(((x1 + dx * cur, y1 + dy * cur), (x2, y2)));
    }
    out
}

/// Computes normalized Z depth `[0.0 = furthest back, 1.0 = frontmost]` for every atom in `mol`,
/// mapping to opacity in `[rear_opacity, 1.0]` when `mol.has_3d_coords` is true.
pub fn compute_3d_atom_opacities(mol: &ChemMolecule, rear_opacity: f32) -> Vec<f32> {
    let n = mol.atoms.len();
    let ro = rear_opacity.clamp(0.0, 1.0);
    if !mol.has_3d_coords || ro >= 0.995 || n == 0 {
        return vec![1.0; n];
    }
    let mut min_z = f64::INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    for a in &mol.atoms {
        min_z = min_z.min(a.z_coord);
        max_z = max_z.max(a.z_coord);
    }
    let span = max_z - min_z;
    if span < 1e-4 {
        return vec![1.0; n];
    }
    mol.atoms
        .iter()
        .map(|a| {
            let t = ((a.z_coord - min_z) / span).clamp(0.0, 1.0) as f32;
            ro + t * (1.0 - ro)
        })
        .collect()
}

// ============================================================================
// Catalytic Cycle Builder & Clean-Up
// ============================================================================

/// Arranges `mol_indices` (3..=8 intermediates) neatly around an ellipse `(cx, cy, rx, ry)`
/// and adds curved reaction arrows connecting consecutive intermediates in the cycle!
pub fn arrange_catalytic_cycle(
    scene: &mut DocumentScene,
    next_id: &mut u32,
    mol_indices: &[usize],
    cx: f64,
    cy: f64,
    rx: f64,
    ry: f64,
    clockwise: bool,
) -> usize {
    let valid: Vec<usize> = mol_indices
        .iter()
        .copied()
        .filter(|&mi| mi < scene.molecules.len() && !scene.molecules[mi].atoms.is_empty())
        .collect();
    let n = valid.len();
    if n < 2 {
        return 0;
    }

    // Remove previous cycle arrows so Clean Up Cycle is idempotent
    scene.arrows.retain(|a| a.arrowhead_head != "Curved");

    let mut centers = Vec::with_capacity(n);
    for (k, &mi) in valid.iter().enumerate() {
        let dir = if clockwise { 1.0 } else { -1.0 };
        let ang = -std::f64::consts::FRAC_PI_2 + dir * (k as f64) * std::f64::consts::TAU / (n as f64);
        let tx = cx + rx * ang.cos();
        let ty = cy + ry * ang.sin();
        if let Some(m) = scene.molecules.get_mut(mi) {
            let (mut x1, mut y1, mut x2, mut y2) = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
            for a in &m.atoms {
                x1 = x1.min(a.x);
                y1 = y1.min(a.y);
                x2 = x2.max(a.x);
                y2 = y2.max(a.y);
            }
            if x1.is_finite() {
                let dx = tx - 0.5 * (x1 + x2);
                let dy = ty - 0.5 * (y1 + y2);
                for a in &mut m.atoms {
                    a.x += dx;
                    a.y += dy;
                }
            }
            let rad = 0.5 * (x2 - x1).hypot(y2 - y1).max(20.0);
            centers.push((tx, ty, rad));
        }
    }

    // Connect consecutive intermediates along the cycle with clean curved arrows trimmed to bounding radii
    for k in 0..centers.len() {
        let (x1, y1, r1) = centers[k];
        let (x2, y2, r2) = centers[(k + 1) % centers.len()];
        let dvx = x2 - x1;
        let dvy = y2 - y1;
        let dist = dvx.hypot(dvy).max(1.0);
        let ux = dvx / dist;
        let uy = dvy / dist;
        let pad1 = (r1 + 10.0).min(dist * 0.35);
        let pad2 = (r2 + 10.0).min(dist * 0.35);
        let tail = (x1 + ux * pad1, y1 + uy * pad1);
        let head = (x2 - ux * pad2, y2 - uy * pad2);
        let aid = *next_id;
        *next_id += 1;
        scene.arrows.push(SceneArrow {
            id: aid,
            tail,
            head,
            arrowhead_head: "Curved".to_string(),
            arrowhead_tail: None,
            no_go: None,
            dipole: false,
        });
    }
    n
}

/// Inserts a starter 4-step Pd-catalyzed cross-coupling catalytic cycle ready to edit!
pub fn insert_starter_catalytic_cycle(scene: &mut DocumentScene, next_id: &mut u32, cx: f64, cy: f64) -> Vec<usize> {
    let labels = ["L₂Pd(0)", "Ar–Pd(II)–X", "Ar–Pd(II)–R", "L₂Pd(0) + Ar–R"];
    let mut idxs = Vec::new();
    for lbl in labels {
        let mut m = ChemMolecule::new();
        let aid = *next_id;
        *next_id += 1;
        let mut a = ChemAtom::new(aid, 46); // Pd
        a.label = Some(lbl.to_string());
        a.x = cx;
        a.y = cy;
        m.add_atom(a);
        // Add tiny invisible structural anchor or single ligand stub so bbox is comfortable
        let bid = *next_id;
        *next_id += 1;
        let _ = (bid, ChemBond { id: 0, begin_atom: 0, end_atom: 0, order: 1, display: BondDisplay::Solid, visible: true, ignore_warnings: false });
        idxs.push(scene.molecules.len());
        scene.molecules.push(m);
    }
    arrange_catalytic_cycle(scene, next_id, &idxs, cx, cy, 105.0, 88.0, true);
    idxs
}

// ============================================================================
// Live Linked Property & Name Captions
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinkedCaptionKind {
    Formula,
    MolecularWeight,
    ExactMass,
    IupacName,
    FullSummary,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinkedCaption {
    pub caption_id: u32,
    pub molecule_idx: usize,
    pub kind: LinkedCaptionKind,
    pub decimals: usize,
}

pub fn format_linked_caption_text(
    mol: &ChemMolecule,
    kind: LinkedCaptionKind,
    decimals: usize,
    live_name: &str,
) -> String {
    let rep = compute_analysis_report(mol);
    match kind {
        LinkedCaptionKind::Formula => format!("Formula: {}", rep.formula),
        LinkedCaptionKind::MolecularWeight => format!("MW: {:.*}", decimals, rep.mol_weight),
        LinkedCaptionKind::ExactMass => format!("Exact Mass: {:.*}", decimals, rep.exact_mass),
        LinkedCaptionKind::IupacName => live_name.to_string(),
        LinkedCaptionKind::FullSummary => format!(
            "{}  |  {}  |  MW: {:.*}",
            live_name, rep.formula, decimals, rep.mol_weight
        ),
    }
}

// ============================================================================
// Export Preferences (DPI, Margin, Copy Background)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CopyBackground {
    Transparent,
    White,
    DarkSlide,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExportOptions {
    pub png_dpi: u32,
    pub svg_padding_pt: f64,
    pub copy_bg: CopyBackground,
}

impl Default for ExportOptions {
    fn default() -> Self {
        Self {
            png_dpi: 300,
            svg_padding_pt: 20.0,
            copy_bg: CopyBackground::Transparent,
        }
    }
}

impl ExportOptions {
    /// Converts target DPI (`150`, `300`, `600`) to raster pixel scale factor relative to 72 pt/in.
    pub fn png_scale(&self) -> f64 {
        (self.png_dpi.clamp(72, 1200) as f64) / 72.0
    }
}

// ============================================================================
// 2D -> 3D Conformer Lift (Pure-Rust Stereochemistry & Hybridization Preserving)
// ============================================================================

/// Lifts a 2D `ChemMolecule` into a clean 3D conformer (`has_3d_coords = true`) preserving:
/// - Planar sp²/aromatic rings & double bonds
/// - Tetrahedral Z-elevation on wedged (`+z`) and hashed (`-z`) stereobonds
/// - Gentle alternating chair-like puckering on saturated 6-membered rings and sp³ chains
pub fn lift_molecule_2d_to_3d(mol: &ChemMolecule) -> ChemMolecule {
    let mut out = mol.clone();
    if out.atoms.is_empty() {
        return out;
    }
    if out.has_3d_coords && out.atoms.iter().any(|a| a.z_coord.abs() > 0.25) {
        return out;
    }
    let n = out.atoms.len();
    let mut is_unsat = vec![false; n];
    let mut adj: Vec<Vec<(usize, BondDisplay)>> = vec![Vec::new(); n];
    for b in &out.bonds {
        if b.begin_atom < n && b.end_atom < n {
            if b.order >= 2 {
                is_unsat[b.begin_atom] = true;
                is_unsat[b.end_atom] = true;
            }
            adj[b.begin_atom].push((b.end_atom, b.display));
            adj[b.end_atom].push((b.begin_atom, match b.display {
                BondDisplay::WedgeBegin => BondDisplay::WedgeEnd,
                BondDisplay::WedgeEnd => BondDisplay::WedgeBegin,
                BondDisplay::WedgedHashBegin => BondDisplay::WedgedHashEnd,
                BondDisplay::WedgedHashEnd => BondDisplay::WedgedHashBegin,
                other => other,
            }));
        }
    }

    // 1. Honor explicit stereo wedges / hashes first
    let mut z = vec![0.0_f64; n];
    let mut pinned = vec![false; n];
    for b in &out.bonds {
        if b.begin_atom >= n || b.end_atom >= n {
            continue;
        }
        let dlen = (out.atoms[b.end_atom].x - out.atoms[b.begin_atom].x)
            .hypot(out.atoms[b.end_atom].y - out.atoms[b.begin_atom].y)
            .max(10.0);
        let lift = dlen * 0.62;
        match b.display {
            BondDisplay::WedgeBegin | BondDisplay::HollowWedgeBegin => {
                z[b.end_atom] = lift;
                pinned[b.end_atom] = true;
            }
            BondDisplay::WedgeEnd | BondDisplay::HollowWedgeEnd => {
                z[b.begin_atom] = lift;
                pinned[b.begin_atom] = true;
            }
            BondDisplay::WedgedHashBegin | BondDisplay::Hash => {
                z[b.end_atom] = -lift;
                pinned[b.end_atom] = true;
            }
            BondDisplay::WedgedHashEnd => {
                z[b.begin_atom] = -lift;
                pinned[b.begin_atom] = true;
            }
            _ => {}
        }
    }

    // 2. Apply alternating sp3 tetrahedral / chair-like Z puckering on non-unsaturated atoms
    let mut visited = vec![false; n];
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut queue = std::collections::VecDeque::new();
        visited[start] = true;
        queue.push_back((start, 0usize));
        while let Some((u, depth)) = queue.pop_front() {
            if !pinned[u] && !is_unsat[u] && adj[u].len() >= 2 {
                let sign = if depth % 2 == 0 { 1.0 } else { -1.0 };
                let amp = if adj[u].len() >= 3 { 4.8 } else { 3.2 };
                z[u] = sign * amp;
            }
            for &(v, _) in &adj[u] {
                if !visited[v] {
                    visited[v] = true;
                    queue.push_back((v, depth + 1));
                }
            }
        }
    }

    // Center mean Z to 0
    let mean_z = z.iter().sum::<f64>() / (n as f64);
    for (i, a) in out.atoms.iter_mut().enumerate() {
        a.z_coord = z[i] - mean_z;
    }
    out.has_3d_coords = true;
    out
}

// ============================================================================
// Smart Reverse-Order Group Typing (Left-Facing Bond End Abbreviations)
// ============================================================================

/// Recognizes reversed left-facing functional group abbreviations (e.g. `O2N`, `MeO`, `AcO`,
/// `F3C`, `HO`, `H2N`, `BocHN`, `HO2C`) and returns `(canonical_nickname, formatted_display)`.
pub fn resolve_reverse_abbreviation(raw: &str) -> Option<(&'static str, &'static str)> {
    let clean: String = raw.trim().replace(['-', '–', '—'], "");
    match clean.as_str() {
        "HO" => Some(("OH", "HO")),
        "HS" => Some(("SH", "HS")),
        "H2N" | "H₂N" => Some(("NH2", "H₂N")),
        "H3N+" | "H₃N⁺" => Some(("NH3+", "H₃N⁺")),
        "O2N" | "O₂N" => Some(("NO2", "O₂N")),
        "ON" => Some(("NO", "ON")),
        "NC" => Some(("CN", "NC")),
        "MeO" => Some(("OMe", "MeO")),
        "EtO" => Some(("OEt", "EtO")),
        "PrO" => Some(("OPr", "PrO")),
        "iPrO" => Some(("OiPr", "iPrO")),
        "BuO" => Some(("OBu", "BuO")),
        "tBuO" => Some(("OtBu", "tBuO")),
        "PhO" => Some(("OPh", "PhO")),
        "BnO" => Some(("OBn", "BnO")),
        "AcO" => Some(("OAc", "AcO")),
        "BzO" => Some(("OBz", "BzO")),
        "MsO" => Some(("OMs", "MsO")),
        "TsO" => Some(("OTs", "TsO")),
        "TfO" => Some(("OTf", "TfO")),
        "F3C" | "F₃C" => Some(("CF3", "F₃C")),
        "Cl3C" | "Cl₃C" => Some(("CCl3", "Cl₃C")),
        "F3CO" | "F₃CO" => Some(("OCF3", "F₃CO")),
        "MeS" => Some(("SMe", "MeS")),
        "PhS" => Some(("SPh", "PhS")),
        "Me2N" | "Me₂N" => Some(("NMe2", "Me₂N")),
        "Et2N" | "Et₂N" => Some(("NEt2", "Et₂N")),
        "AcHN" => Some(("NHAc", "AcHN")),
        "BocHN" => Some(("NHBoc", "BocHN")),
        "CbzHN" => Some(("NHCbz", "CbzHN")),
        "FmocHN" => Some(("NHFmoc", "FmocHN")),
        "TsHN" => Some(("NHTs", "TsHN")),
        "HO2C" | "HO₂C" | "HOOC" => Some(("COOH", "HO₂C")),
        "MeO2C" | "MeO₂C" => Some(("CO2Me", "MeO₂C")),
        "EtO2C" | "EtO₂C" => Some(("CO2Et", "EtO₂C")),
        "tBuO2C" | "tBuO₂C" => Some(("CO2tBu", "tBuO₂C")),
        "H2NOC" | "H₂NOC" => Some(("CONH2", "H₂NOC")),
        "OHC" => Some(("CHO", "OHC")),
        "Me3Si" | "Me₃Si" => Some(("TMS", "Me₃Si")),
        "TBDMSO" => Some(("OTBDMS", "TBDMSO")),
        "TBSO" => Some(("OTBS", "TBSO")),
        "TIPSO" => Some(("OTIPS", "TIPSO")),
        "THPO" => Some(("OTHP", "THPO")),
        "MOMO" => Some(("OMOM", "MOMO")),
        _ => None,
    }
}

// ============================================================================
// Chemistry-Aware Document Diff / Compare
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiffStatus {
    Unchanged,
    Added,
    Removed,
    Modified,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MolDiffEntry {
    pub status: DiffStatus,
    pub left_index: Option<usize>,
    pub right_index: Option<usize>,
    pub left_formula: String,
    pub right_formula: String,
    pub left_smiles: String,
    pub right_smiles: String,
    pub left_mw: f64,
    pub right_mw: f64,
    pub summary: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DocumentDiff {
    pub left_name: String,
    pub right_name: String,
    pub entries: Vec<MolDiffEntry>,
    pub added_captions: Vec<String>,
    pub removed_captions: Vec<String>,
    pub arrow_delta: isize,
    pub table_delta: isize,
    pub tlc_delta: isize,
}

impl DocumentDiff {
    pub fn has_changes(&self) -> bool {
        self.entries.iter().any(|e| e.status != DiffStatus::Unchanged)
            || !self.added_captions.is_empty()
            || !self.removed_captions.is_empty()
            || self.arrow_delta != 0
            || self.table_delta != 0
            || self.tlc_delta != 0
    }

    pub fn format_text_report(&self) -> String {
        let mut out = format!("Chemical Diff: {}  vs  {}\n", self.left_name, self.right_name);
        let unchanged = self.entries.iter().filter(|e| e.status == DiffStatus::Unchanged).count();
        let modified = self.entries.iter().filter(|e| e.status == DiffStatus::Modified).count();
        let added = self.entries.iter().filter(|e| e.status == DiffStatus::Added).count();
        let removed = self.entries.iter().filter(|e| e.status == DiffStatus::Removed).count();
        out.push_str(&format!(
            "Molecules: {unchanged} unchanged, {modified} modified, {added} added, {removed} removed\n"
        ));
        for e in &self.entries {
            let tag = match e.status {
                DiffStatus::Unchanged => "  [=]",
                DiffStatus::Modified => "  [~]",
                DiffStatus::Added => "  [+]",
                DiffStatus::Removed => "  [-]",
            };
            out.push_str(&format!("{tag} {}\n", e.summary));
        }
        for c in &self.added_captions {
            out.push_str(&format!("  [+caption] {c}\n"));
        }
        for c in &self.removed_captions {
            out.push_str(&format!("  [-caption] {c}\n"));
        }
        if self.arrow_delta != 0 || self.table_delta != 0 || self.tlc_delta != 0 {
            out.push_str(&format!(
                "Objects delta: arrows {:+}, tables {:+}, TLC plates {:+}\n",
                self.arrow_delta, self.table_delta, self.tlc_delta
            ));
        }
        out
    }
}

/// Compares two `InteractiveDocument` instances chemically by SMILES, molecular formula, MW,
/// stereochemistry, captions, and reaction objects rather than raw XML lines.
pub fn compare_documents(left: &crate::InteractiveDocument, right: &crate::InteractiveDocument) -> DocumentDiff {
    let summarize = |m: &ChemMolecule| -> (String, String, f64, usize) {
        let rep = compute_analysis_report(m);
        let smi = chiralcraft_formats::molecule_to_smiles(m);
        let stereo_cnt = m
            .bonds
            .iter()
            .filter(|b| !matches!(b.display, BondDisplay::Solid))
            .count();
        (rep.formula, smi, rep.mol_weight, stereo_cnt)
    };

    let l_info: Vec<_> = left.scene.molecules.iter().map(summarize).collect();
    let r_info: Vec<_> = right.scene.molecules.iter().map(summarize).collect();
    let mut l_matched = vec![false; l_info.len()];
    let mut r_matched = vec![false; r_info.len()];
    let mut entries = Vec::new();

    // Pass 1: exact chemical match (same SMILES + same stereo bond count)
    for (li, (lf, ls, lmw, lst)) in l_info.iter().enumerate() {
        if let Some(ri) = r_info.iter().enumerate().position(|(ri, (rf, rs, _, rst))| {
            !r_matched[ri] && lf == rf && ls == rs && lst == rst
        }) {
            l_matched[li] = true;
            r_matched[ri] = true;
            entries.push(MolDiffEntry {
                status: DiffStatus::Unchanged,
                left_index: Some(li),
                right_index: Some(ri),
                left_formula: lf.clone(),
                right_formula: r_info[ri].0.clone(),
                left_smiles: ls.clone(),
                right_smiles: r_info[ri].1.clone(),
                left_mw: *lmw,
                right_mw: r_info[ri].2,
                summary: format!("Mol #{}: {} ({}) unchanged", li + 1, lf, ls),
            });
        }
    }

    // Pass 2: pair remaining by position/index as Modified
    let rem_l: Vec<usize> = (0..l_info.len()).filter(|&i| !l_matched[i]).collect();
    let rem_r: Vec<usize> = (0..r_info.len()).filter(|&i| !r_matched[i]).collect();
    let pair_n = rem_l.len().min(rem_r.len());
    for k in 0..pair_n {
        let li = rem_l[k];
        let ri = rem_r[k];
        let (lf, ls, lmw, lst) = &l_info[li];
        let (rf, rs, rmw, rst) = &r_info[ri];
        let dmw = rmw - lmw;
        let reason = if ls == rs && lst != rst {
            format!("stereo changed ({lst} -> {rst} stereobonds)")
        } else {
            format!("{lf} ({ls}) -> {rf} ({rs}), ΔMW {dmw:+.2}")
        };
        entries.push(MolDiffEntry {
            status: DiffStatus::Modified,
            left_index: Some(li),
            right_index: Some(ri),
            left_formula: lf.clone(),
            right_formula: rf.clone(),
            left_smiles: ls.clone(),
            right_smiles: rs.clone(),
            left_mw: *lmw,
            right_mw: *rmw,
            summary: format!("Mol #{} -> #{}: {reason}", li + 1, ri + 1),
        });
    }
    for &li in rem_l.iter().skip(pair_n) {
        let (lf, ls, lmw, _) = &l_info[li];
        entries.push(MolDiffEntry {
            status: DiffStatus::Removed,
            left_index: Some(li),
            right_index: None,
            left_formula: lf.clone(),
            right_formula: String::new(),
            left_smiles: ls.clone(),
            right_smiles: String::new(),
            left_mw: *lmw,
            right_mw: 0.0,
            summary: format!("Removed Mol #{}: {} ({}, MW {:.2})", li + 1, lf, ls, lmw),
        });
    }
    for &ri in rem_r.iter().skip(pair_n) {
        let (rf, rs, rmw, _) = &r_info[ri];
        entries.push(MolDiffEntry {
            status: DiffStatus::Added,
            left_index: None,
            right_index: Some(ri),
            left_formula: String::new(),
            right_formula: rf.clone(),
            left_smiles: String::new(),
            right_smiles: rs.clone(),
            left_mw: 0.0,
            right_mw: *rmw,
            summary: format!("Added Mol #{}: {} ({}, MW {:.2})", ri + 1, rf, rs, rmw),
        });
    }

    let l_caps: Vec<String> = left.scene.captions.iter().map(|c| c.text.clone()).collect();
    let r_caps: Vec<String> = right.scene.captions.iter().map(|c| c.text.clone()).collect();
    let added_captions = r_caps.iter().filter(|c| !l_caps.contains(c)).cloned().collect();
    let removed_captions = l_caps.iter().filter(|c| !r_caps.contains(c)).cloned().collect();

    DocumentDiff {
        left_name: left.name.clone(),
        right_name: right.name.clone(),
        entries,
        added_captions,
        removed_captions,
        arrow_delta: (right.scene.arrows.len() as isize) - (left.scene.arrows.len() as isize),
        table_delta: (right.scene.tables.len() as isize) - (left.scene.tables.len() as isize),
        tlc_delta: (right.scene.tlc_plates.len() as isize) - (left.scene.tlc_plates.len() as isize),
    }
}
