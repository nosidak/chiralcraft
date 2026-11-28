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
// 1. Multi-Backend FT-IR Spectroscopy Engine (4000 - 400 cm^-1)
//    Supports Instant Rule-Based, GFN2-xTB / GFN-FF, ORCA/Gaussian DFT, Molden,
//    and optional feature-gated ML backends via a unified IrPrediction model.
// ============================================================================

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IrSpectrumSource {
    RuleBased,
    XtbGfn2,
    XtbGfnFf,
    OrcaDft,
    GaussianDft,
    Molden,
    MlModel(String),
}

impl IrSpectrumSource {
    pub fn label(&self) -> String {
        match self {
            IrSpectrumSource::RuleBased => "Rule-Based (Instant)".to_string(),
            IrSpectrumSource::XtbGfn2 => "GFN2-xTB (Harmonic)".to_string(),
            IrSpectrumSource::XtbGfnFf => "GFN-FF (Harmonic)".to_string(),
            IrSpectrumSource::OrcaDft => "ORCA DFT (Harmonic)".to_string(),
            IrSpectrumSource::GaussianDft => "Gaussian DFT (Harmonic)".to_string(),
            IrSpectrumSource::Molden => "Molden (Harmonic)".to_string(),
            IrSpectrumSource::MlModel(name) => format!("ML ({name})"),
        }
    }

    pub fn is_quantum(&self) -> bool {
        matches!(
            self,
            IrSpectrumSource::XtbGfn2
                | IrSpectrumSource::XtbGfnFf
                | IrSpectrumSource::OrcaDft
                | IrSpectrumSource::GaussianDft
                | IrSpectrumSource::Molden
        )
    }

    pub fn recommended_scale_factor(&self) -> f64 {
        match self {
            IrSpectrumSource::RuleBased | IrSpectrumSource::MlModel(_) => 1.0,
            IrSpectrumSource::XtbGfn2 | IrSpectrumSource::XtbGfnFf => 0.970,
            IrSpectrumSource::OrcaDft | IrSpectrumSource::GaussianDft | IrSpectrumSource::Molden => 0.967,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IrLineShape {
    Lorentzian,
    Gaussian,
    PseudoVoigt,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IrYAxis {
    TransmittancePct,
    Absorbance,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IrRenderParams {
    pub scale_factor: f64,
    /// Global FWHM (cm⁻¹) used for quantum/stick spectra (or when `use_per_band_widths == false`).
    pub fwhm_cm1: f64,
    /// When true, uses each band's own `width_cm1` HWHM (preserves broad O–H vs sharp C=O in rule mode).
    pub use_per_band_widths: bool,
    pub line_shape: IrLineShape,
    pub y_axis: IrYAxis,
    pub x_min_cm1: f64,
    pub x_max_cm1: f64,
}

impl Default for IrRenderParams {
    fn default() -> Self {
        Self {
            scale_factor: 1.0,
            fwhm_cm1: 18.0,
            use_per_band_widths: true,
            line_shape: IrLineShape::Lorentzian,
            y_axis: IrYAxis::TransmittancePct,
            x_min_cm1: 400.0,
            x_max_cm1: 4000.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct IrBand {
    /// Raw unscaled wavenumber (cm⁻¹). Negative values represent imaginary modes.
    pub wavenumber_cm1: f64,
    /// Normalized absorption depth `0.0..=1.0` used for plotting and relative strength.
    pub intensity: f64,
    /// Physical IR intensity in `km/mol` when computed by xTB / ORCA / Gaussian.
    pub intensity_km_mol: Option<f64>,
    /// Band HWHM (cm⁻¹) for rule-based broadening.
    pub width_cm1: f64,
    /// Spectroscopic strength descriptor (`"vs"`, `"s"`, `"m"`, `"w"`, `"br"`, `"imag"`).
    pub strength: String,
    pub assignment: String,
    pub symmetry: Option<String>,
    pub is_imaginary: bool,
    pub mode_index: Option<usize>,
    /// Canvas heavy-atom indices highlighted when this band is selected.
    pub atom_indices: Vec<usize>,
    /// Optional per-atom kinetic-energy / displacement weights `(atom_idx, 0.0..=1.0)`.
    pub atom_weights: Vec<(usize, f64)>,
    /// Optional 3D Cartesian normal-mode displacements per atom (in explicit-H 3D geometry order).
    pub normal_mode_vectors: Option<Vec<[f64; 3]>>,
}

impl IrBand {
    pub fn rule(
        wavenumber_cm1: f64,
        intensity: f64,
        width_cm1: f64,
        strength: &str,
        assignment: impl Into<String>,
        atom_indices: Vec<usize>,
    ) -> Self {
        Self {
            wavenumber_cm1,
            intensity: intensity.clamp(0.0, 1.0),
            intensity_km_mol: None,
            width_cm1,
            strength: strength.to_string(),
            assignment: assignment.into(),
            symmetry: None,
            is_imaginary: wavenumber_cm1 < 0.0,
            mode_index: None,
            atom_indices,
            atom_weights: Vec::new(),
            normal_mode_vectors: None,
        }
    }

    pub fn scaled_wavenumber(&self, scale_factor: f64) -> f64 {
        if self.is_imaginary {
            self.wavenumber_cm1
        } else {
            self.wavenumber_cm1 * scale_factor
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct IrPrediction {
    pub source: IrSpectrumSource,
    pub bands: Vec<IrBand>,
    pub imaginary_count: usize,
    pub default_scale_factor: f64,
    pub report_text: String,
    /// Explicit-H 3D equilibrium geometry matching `IrBand::normal_mode_vectors` (for 3D vibration animation).
    pub mode_geometry_3d: Option<ChemMolecule>,
}

/// Converts `%Transmittance` (`0..=100`) to `Absorbance` (`-log10(T)`).
#[inline]
pub fn transmittance_pct_to_absorbance(t_pct: f64) -> f64 {
    let t = (t_pct / 100.0).clamp(1e-4, 1.0);
    -t.log10()
}

/// Converts `Absorbance` (`>= 0`) to `%Transmittance` (`0..=100`).
#[inline]
pub fn absorbance_to_transmittance_pct(abs: f64) -> f64 {
    (100.0 * 10.0_f64.powf(-abs.max(0.0))).clamp(0.5, 100.0)
}

fn eval_lineshape(wn: f64, center: f64, hwhm: f64, shape: IrLineShape) -> f64 {
    let g = hwhm.max(2.0);
    let dx = (wn - center) / g;
    match shape {
        IrLineShape::Lorentzian => 1.0 / (1.0 + dx * dx),
        IrLineShape::Gaussian => (-std::f64::consts::LN_2 * dx * dx).exp(),
        IrLineShape::PseudoVoigt => {
            let lor = 1.0 / (1.0 + dx * dx);
            let gau = (-std::f64::consts::LN_2 * dx * dx).exp();
            0.5 * lor + 0.5 * gau
        }
    }
}

impl IrPrediction {
    /// Backwards-compatible default `%Transmittance` curve `(cm^-1, %T 0..100)` from `4000` down to `400 cm^-1`.
    pub fn simulate_transmittance(&self, samples: usize) -> Vec<(f64, f64)> {
        let params = IrRenderParams {
            scale_factor: self.default_scale_factor,
            use_per_band_widths: !self.source.is_quantum(),
            y_axis: IrYAxis::TransmittancePct,
            ..IrRenderParams::default()
        };
        self.simulate_curve(&params, samples)
    }

    /// Evaluates a smooth Absorbance curve `(cm^-1, A >= 0)` from `x_max_cm1` down to `x_min_cm1`.
    pub fn simulate_absorbance(&self, params: &IrRenderParams, samples: usize) -> Vec<(f64, f64)> {
        if samples < 2 {
            return Vec::new();
        }
        // Accept a swapped range (min > max) instead of collapsing it to a 10 cm⁻¹ window.
        let (a, b) = (params.x_min_cm1.min(params.x_max_cm1), params.x_min_cm1.max(params.x_max_cm1));
        let hi = b.max(a + 10.0);
        let lo = a.min(hi - 10.0);
        let span = hi - lo;
        let global_hwhm = (params.fwhm_cm1 * 0.5).max(2.0);
        let mut out = Vec::with_capacity(samples);
        for k in 0..samples {
            let wn = hi - span * (k as f64) / ((samples - 1) as f64);
            let mut abs = 0.0_f64;
            for b in &self.bands {
                if b.is_imaginary || b.wavenumber_cm1 <= 0.0 || b.intensity <= 1e-6 {
                    continue;
                }
                let center = b.scaled_wavenumber(params.scale_factor);
                let hwhm = if params.use_per_band_widths && !self.source.is_quantum() {
                    b.width_cm1.max(4.0)
                } else {
                    global_hwhm
                };
                if (wn - center).abs() > hwhm * 18.0 {
                    continue;
                }
                abs += b.intensity * eval_lineshape(wn, center, hwhm, params.line_shape);
            }
            out.push((wn, abs));
        }
        out
    }

    /// Evaluates either `%Transmittance` or `Absorbance` curve `(cm^-1, y)` from `x_max_cm1` down to `x_min_cm1`.
    pub fn simulate_curve(&self, params: &IrRenderParams, samples: usize) -> Vec<(f64, f64)> {
        let abs_curve = self.simulate_absorbance(params, samples);
        match params.y_axis {
            IrYAxis::Absorbance => abs_curve,
            IrYAxis::TransmittancePct => abs_curve
                .into_iter()
                .map(|(wn, a)| (wn, absorbance_to_transmittance_pct(a * 0.92)))
                .collect(),
        }
    }

    /// Formats an ACS/JOC-style `IR (cm^-1): νmax ...` experimental summary line.
    pub fn format_report(&self, scale_factor: f64) -> String {
        let mut real_bands: Vec<&IrBand> = self
            .bands
            .iter()
            .filter(|b| !b.is_imaginary && b.wavenumber_cm1 >= 400.0 && b.intensity >= 0.15)
            .collect();
        if self.source.is_quantum() {
            real_bands.sort_by(|a, b| {
                b.intensity
                    .partial_cmp(&a.intensity)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            real_bands.truncate(10);
            real_bands.sort_by(|a, b| {
                b.wavenumber_cm1
                    .partial_cmp(&a.wavenumber_cm1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        } else {
            real_bands.truncate(8);
        }
        let summary: Vec<String> = real_bands
            .iter()
            .map(|b| format!("{:.0} ({})", b.scaled_wavenumber(scale_factor), b.strength))
            .collect();
        let src_tag = match &self.source {
            IrSpectrumSource::RuleBased => "neat".to_string(),
            other => format!("calc, {}", other.label()),
        };
        format!("IR ({src_tag}, cm⁻¹): νmax {}.", summary.join(", "))
    }

    /// Exports the broadened spectrum curve as a 2-column CSV (`wavenumber_cm1,value`).
    pub fn to_csv(&self, params: &IrRenderParams, samples: usize) -> String {
        let header = match params.y_axis {
            IrYAxis::TransmittancePct => "wavenumber_cm1,transmittance_pct\n",
            IrYAxis::Absorbance => "wavenumber_cm1,absorbance\n",
        };
        let mut out = String::from(header);
        for (wn, y) in self.simulate_curve(params, samples) {
            out.push_str(&format!("{wn:.2},{y:.5}\n"));
        }
        out
    }

    /// Exports the broadened spectrum as a standard IUPAC JCAMP-DX 5.00 IR `.jdx` file.
    pub fn to_jcamp_dx(&self, title: &str, params: &IrRenderParams, samples: usize) -> String {
        let curve = self.simulate_curve(params, samples.max(16));
        let first_x = curve.first().map(|p| p.0).unwrap_or(params.x_max_cm1);
        let last_x = curve.last().map(|p| p.0).unwrap_or(params.x_min_cm1);
        let yunits = match params.y_axis {
            IrYAxis::TransmittancePct => "TRANSMITTANCE",
            IrYAxis::Absorbance => "ABSORBANCE",
        };
        let mut out = format!(
            "##TITLE={title}\n##JCAMP-DX=5.00\n##DATA TYPE=INFRARED SPECTRUM\n##ORIGIN=ChiralCraft ({})\n##XUNITS=1/CM\n##YUNITS={yunits}\n##FIRSTX={first_x:.4}\n##LASTX={last_x:.4}\n##XFACTOR=1.0\n##YFACTOR=1.0\n##NPOINTS={}\n##XYDATA=(X++(Y..Y))\n",
            self.source.label(),
            curve.len()
        );
        for chunk in curve.chunks(8) {
            if let Some(first) = chunk.first() {
                out.push_str(&format!("{:.3}", first.0));
                for &(_, y) in chunk {
                    let y_val = if params.y_axis == IrYAxis::TransmittancePct {
                        y / 100.0
                    } else {
                        y
                    };
                    out.push_str(&format!(" {y_val:.5}"));
                }
                out.push('\n');
            }
        }
        out.push_str("##END=\n");
        out
    }
}

/// Pluggable backend trait for FT-IR prediction engines (pure-Rust rules, xTB, DFT, or optional ML).
pub trait IrBackend {
    fn name(&self) -> &str;
    fn predict(&self, mol: &ChemMolecule) -> Result<IrPrediction, String>;
}

/// Instant pure-Rust functional-group correlation backend.
pub struct RuleBasedIrBackend;

impl IrBackend for RuleBasedIrBackend {
    fn name(&self) -> &str {
        "Rule-Based (Pure Rust)"
    }

    fn predict(&self, mol: &ChemMolecule) -> Result<IrPrediction, String> {
        Ok(predict_ir_spectrum(mol))
    }
}

/// Feature-gated extension seam (`--features ml-ir`) for plugging in an ONNX / Chemprop-IR model
/// without adding ML runtime dependencies to the default pure-Rust ChiralCraft binary.
#[cfg(feature = "ml-ir")]
pub trait MlIrModelHook: IrBackend + Send + Sync {
    fn model_path(&self) -> &str;
}

// ============================================================================
// 1b. Experimental Peak Picking & Quantitative Spectral Comparison (HQI / r)
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub struct IrPeak {
    pub wavenumber_cm1: f64,
    pub y_value: f64,
    pub prominence: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IrComparisonMetrics {
    /// Hit Quality Index (`0.0..=100.0`, squared cosine similarity in absorbance space).
    pub hqi_pct: f64,
    /// Mean-centered Pearson correlation (`-1.0..=1.0`).
    pub pearson_r: f64,
    /// Uncentered cosine similarity (`0.0..=1.0`).
    pub cosine_sim: f64,
    pub overlap_points: usize,
}

/// Detects prominent absorbance maxima from an IR curve `(wavenumber_cm1, absorbance)`.
pub fn pick_ir_peaks(curve_wn_abs: &[(f64, f64)], min_prominence: f64, max_peaks: usize) -> Vec<IrPeak> {
    if curve_wn_abs.len() < 5 || max_peaks == 0 {
        return Vec::new();
    }
    let mut pts = curve_wn_abs.to_vec();
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let max_y = pts.iter().map(|p| p.1).fold(0.0_f64, f64::max);
    if max_y <= 1e-9 {
        return Vec::new();
    }
    let thresh = min_prominence.max(0.02) * max_y;
    let n = pts.len();
    let win = (n / 40).clamp(2, 18);
    let mut found = Vec::new();
    for i in 1..(n - 1) {
        let y = pts[i].1;
        if y >= pts[i - 1].1 && y > pts[i + 1].1 && y >= thresh {
            let lo = i.saturating_sub(win);
            let hi = (i + win).min(n - 1);
            let left_min = pts[lo..i].iter().map(|p| p.1).fold(y, f64::min);
            let right_min = pts[(i + 1)..=hi].iter().map(|p| p.1).fold(y, f64::min);
            let prom = y - left_min.max(right_min);
            if prom >= thresh * 0.45 {
                found.push(IrPeak {
                    wavenumber_cm1: pts[i].0,
                    y_value: y,
                    prominence: prom,
                });
            }
        }
    }
    found.sort_by(|a, b| b.prominence.partial_cmp(&a.prominence).unwrap_or(std::cmp::Ordering::Equal));
    let mut kept: Vec<IrPeak> = Vec::new();
    for pk in found {
        if !kept.iter().any(|k| (k.wavenumber_cm1 - pk.wavenumber_cm1).abs() < 18.0) {
            kept.push(pk);
            if kept.len() >= max_peaks {
                break;
            }
        }
    }
    kept.sort_by(|a, b| b.wavenumber_cm1.partial_cmp(&a.wavenumber_cm1).unwrap_or(std::cmp::Ordering::Equal));
    kept
}

fn sample_sorted_linear(sorted_pts: &[(f64, f64)], x: f64) -> Option<f64> {
    if sorted_pts.len() < 2 {
        return None;
    }
    if x < sorted_pts.first()?.0 || x > sorted_pts.last()?.0 {
        return None;
    }
    let idx = sorted_pts.partition_point(|&(px, _)| px < x);
    if idx == 0 {
        return Some(sorted_pts[0].1);
    }
    if idx >= sorted_pts.len() {
        return Some(sorted_pts[sorted_pts.len() - 1].1);
    }
    let (x0, y0) = sorted_pts[idx - 1];
    let (x1, y1) = sorted_pts[idx];
    if (x1 - x0).abs() < 1e-12 {
        Some(y0)
    } else {
        let t = (x - x0) / (x1 - x0);
        Some(y0 + t * (y1 - y0))
    }
}

/// Computes Hit Quality Index (HQI %), Pearson correlation `r`, and cosine similarity between
/// two IR absorbance spectra `(wavenumber_cm1, absorbance)` over `[wn_min, wn_max]`.
pub fn compare_ir_spectra(
    calc_wn_abs: &[(f64, f64)],
    exp_wn_abs: &[(f64, f64)],
    wn_min: f64,
    wn_max: f64,
) -> IrComparisonMetrics {
    let mut a: Vec<(f64, f64)> = calc_wn_abs.iter().copied().filter(|(x, y)| x.is_finite() && y.is_finite()).collect();
    let mut b: Vec<(f64, f64)> = exp_wn_abs.iter().copied().filter(|(x, y)| x.is_finite() && y.is_finite()).collect();
    a.sort_by(|u, v| u.0.partial_cmp(&v.0).unwrap_or(std::cmp::Ordering::Equal));
    b.sort_by(|u, v| u.0.partial_cmp(&v.0).unwrap_or(std::cmp::Ordering::Equal));
    if a.len() < 2 || b.len() < 2 {
        return IrComparisonMetrics {
            hqi_pct: 0.0,
            pearson_r: 0.0,
            cosine_sim: 0.0,
            overlap_points: 0,
        };
    }
    let lo = wn_min.max(a.first().unwrap().0).max(b.first().unwrap().0);
    let hi = wn_max.min(a.last().unwrap().0).min(b.last().unwrap().0);
    if hi <= lo + 4.0 {
        return IrComparisonMetrics {
            hqi_pct: 0.0,
            pearson_r: 0.0,
            cosine_sim: 0.0,
            overlap_points: 0,
        };
    }

    let step = 2.0_f64;
    let count = (((hi - lo) / step).floor() as usize).clamp(10, 2000);
    let mut va = Vec::with_capacity(count);
    let mut vb = Vec::with_capacity(count);
    for k in 0..=count {
        let wn = lo + (hi - lo) * (k as f64) / (count as f64);
        if let (Some(ya), Some(yb)) = (sample_sorted_linear(&a, wn), sample_sorted_linear(&b, wn)) {
            va.push(ya.max(0.0));
            vb.push(yb.max(0.0));
        }
    }
    let m = va.len();
    if m < 4 {
        return IrComparisonMetrics {
            hqi_pct: 0.0,
            pearson_r: 0.0,
            cosine_sim: 0.0,
            overlap_points: m,
        };
    }

    let dot: f64 = va.iter().zip(&vb).map(|(x, y)| x * y).sum();
    let na: f64 = va.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb: f64 = vb.iter().map(|y| y * y).sum::<f64>().sqrt();
    let cos = if na > 1e-12 && nb > 1e-12 {
        (dot / (na * nb)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let hqi_pct = (cos * cos * 100.0).clamp(0.0, 100.0);

    let ma = va.iter().sum::<f64>() / (m as f64);
    let mb = vb.iter().sum::<f64>() / (m as f64);
    let mut cov = 0.0;
    let mut sda = 0.0;
    let mut sdb = 0.0;
    for (x, y) in va.iter().zip(&vb) {
        let dx = x - ma;
        let dy = y - mb;
        cov += dx * dy;
        sda += dx * dx;
        sdb += dy * dy;
    }
    let pearson_r = if sda > 1e-12 && sdb > 1e-12 {
        (cov / (sda.sqrt() * sdb.sqrt())).clamp(-1.0, 1.0)
    } else {
        0.0
    };

    IrComparisonMetrics {
        hqi_pct,
        pearson_r,
        cosine_sim: cos,
        overlap_points: m,
    }
}

/// Shortest cycle length containing `start` (up to `max_len`), used for carbonyl ring-strain rules.
fn smallest_ring_size(adj: &[Vec<(usize, i32)>], start: usize, max_len: usize) -> Option<usize> {
    let mut best: Option<usize> = None;
    for &(nbr, _) in &adj[start] {
        // BFS from `nbr` to `start` without using the direct edge `(start, nbr)`
        let mut dist = vec![usize::MAX; adj.len()];
        let mut q = std::collections::VecDeque::new();
        dist[nbr] = 1;
        q.push_back(nbr);
        while let Some(u) = q.pop_front() {
            let d = dist[u];
            if d >= max_len || best.map(|b| d + 1 >= b).unwrap_or(false) {
                continue;
            }
            for &(v, _) in &adj[u] {
                if u == nbr && v == start {
                    continue;
                }
                if v == start {
                    let cyc = d + 1;
                    best = Some(best.map(|b| b.min(cyc)).unwrap_or(cyc));
                    break;
                }
                if dist[v] == usize::MAX {
                    dist[v] = d + 1;
                    q.push_back(v);
                }
            }
        }
    }
    best
}

/// Predicts diagnostic FT-IR vibrational bands and exact atom indices for `mol` (Tier 1 Instant Rule Engine).
pub fn predict_ir_spectrum(mol: &ChemMolecule) -> IrPrediction {
    let n = mol.atoms.len();
    let adj = build_adj(mol);
    let val = CcVal::compute(mol);
    let arom = CcAromaticityInfo::perceive(mol);
    let mut bands: Vec<IrBand> = Vec::new();

    let mut has_sp3_ch = false;
    let mut has_ch2 = false;
    let mut has_ch3 = false;
    let mut has_arom_ch = false;
    let mut has_arom_ring = false;
    let mut sp3_c_atoms = Vec::new();
    let mut ch2_atoms = Vec::new();
    let mut ch3_atoms = Vec::new();
    let mut arom_atoms = Vec::new();

    let tot_h_of = |idx: usize| -> u8 {
        let imp_h = val.implicit_h_counts.get(idx).copied().unwrap_or(0);
        let exp_h = adj[idx].iter().filter(|&&(k, _)| mol.atoms[k].z == 1).count() as u8;
        imp_h.saturating_add(exp_h)
    };

    let mut handled_anhydride_o: BTreeSet<usize> = BTreeSet::new();

    for i in 0..n {
        let z = mol.atoms[i].z;
        let tot_h = tot_h_of(i);

        // 1. O-H stretches (carboxylic acid, phenol, alcohol)
        if z == 8 && tot_h >= 1 {
            let is_acid = adj[i].iter().any(|&(c, _)| {
                mol.atoms[c].z == 6 && adj[c].iter().any(|&(o2, o)| o2 != i && mol.atoms[o2].z == 8 && o == 2)
            });
            let is_phenol = adj[i].iter().any(|&(c, _)| arom.is_atom_aromatic(c));
            if is_acid {
                bands.push(IrBand::rule(
                    3000.0,
                    0.80,
                    210.0,
                    "br",
                    "Carboxylic Acid O–H stretch (very broad)",
                    vec![i],
                ));
            } else {
                bands.push(IrBand::rule(
                    if is_phenol { 3380.0 } else { 3340.0 },
                    0.74,
                    110.0,
                    "br",
                    if is_phenol { "Phenol O–H stretch" } else { "Alcohol O–H stretch (H-bonded)" },
                    vec![i],
                ));
            }
        // 2. N-H stretches (primary NH2 doublet vs secondary N-H, amine vs amide)
        } else if z == 7 && tot_h >= 1 {
            let is_amide = adj[i].iter().any(|&(c, _)| {
                mol.atoms[c].z == 6 && adj[c].iter().any(|&(o2, o)| mol.atoms[o2].z == 8 && o == 2)
            });
            if tot_h >= 2 {
                if is_amide {
                    bands.push(IrBand::rule(3355.0, 0.54, 34.0, "m", "Primary Amide N–H asym stretch", vec![i]));
                    bands.push(IrBand::rule(3185.0, 0.50, 38.0, "m", "Primary Amide N–H sym stretch", vec![i]));
                } else {
                    bands.push(IrBand::rule(3400.0, 0.48, 36.0, "m", "Primary Amine N–H asym stretch", vec![i]));
                    bands.push(IrBand::rule(3320.0, 0.46, 36.0, "m", "Primary Amine N–H sym stretch", vec![i]));
                }
                bands.push(IrBand::rule(1615.0, 0.52, 24.0, "m", "Primary NH₂ scissoring bend", vec![i]));
            } else {
                bands.push(IrBand::rule(
                    if is_amide { 3305.0 } else { 3345.0 },
                    0.52,
                    40.0,
                    "m",
                    if is_amide { "Secondary Amide N–H stretch" } else { "Secondary Amine N–H stretch" },
                    vec![i],
                ));
            }
        // 3. Thiol S-H stretch
        } else if z == 16 && tot_h >= 1 {
            bands.push(IrBand::rule(2565.0, 0.28, 22.0, "w", "Thiol S–H stretch", vec![i]));
        // 4. Carbon centers (aromatic, sp3, carbonyls, anhydrides, aldehydes, terminal alkynes, cumulenes)
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
                if tot_h == 2 {
                    has_ch2 = true;
                    if ch2_atoms.len() < 4 {
                        ch2_atoms.push(i);
                    }
                } else if tot_h >= 3 {
                    has_ch3 = true;
                    if ch3_atoms.len() < 4 {
                        ch3_atoms.push(i);
                    }
                }
            }

            // Terminal alkyne ≡C–H stretch
            if tot_h >= 1 && adj[i].iter().any(|&(k, o)| o == 3 && mol.atoms[k].z == 6) {
                bands.push(IrBand::rule(3300.0, 0.72, 18.0, "s", "Terminal Alkyne ≡C–H stretch", vec![i]));
            }

            // Cumulene / Allene C=C=C or Isocyanate -N=C=O
            let dbl_nbrs: Vec<usize> = adj[i].iter().filter(|&&(_, o)| o == 2).map(|&(k, _)| k).collect();
            if dbl_nbrs.len() == 2 {
                let (za, zb) = (mol.atoms[dbl_nbrs[0]].z, mol.atoms[dbl_nbrs[1]].z);
                let mut at = vec![i, dbl_nbrs[0], dbl_nbrs[1]];
                if (za == 7 && zb == 8) || (za == 8 && zb == 7) {
                    bands.push(IrBand::rule(2270.0, 0.92, 20.0, "vs", "Isocyanate –N=C=O asym stretch", at));
                    continue;
                } else if za == 6 && zb == 6 {
                    bands.push(IrBand::rule(1955.0, 0.68, 20.0, "m", "Allene C=C=C asym stretch", at));
                } else if za == 8 && zb == 8 {
                    at.sort_unstable();
                    bands.push(IrBand::rule(2349.0, 0.90, 18.0, "vs", "CO₂ asym stretch", at));
                    continue;
                }
            }

            // Metal carbonyl M–C≡O (triple or double bonded to O and bonded to a d-block transition metal)
            if let Some(&(o_idx, ord)) = adj[i].iter().find(|&&(k, o)| mol.atoms[k].z == 8 && o >= 2) {
                let bonded_metal = adj[i].iter().find(|&&(k, _)| {
                    let kz = mol.atoms[k].z;
                    (21..=30).contains(&kz) || (39..=48).contains(&kz) || (72..=80).contains(&kz)
                });
                if let Some(&(m_idx, _)) = bonded_metal {
                    if ord == 3 || adj[i].len() == 2 {
                        bands.push(IrBand::rule(
                            1985.0,
                            0.94,
                            18.0,
                            "vs",
                            "Metal Carbonyl M–C≡O stretch",
                            vec![m_idx, i, o_idx],
                        ));
                        continue;
                    }
                }
            }

            // Carbonyl C=O diagnostic frequencies (with ring-strain, anhydride, and conjugation shifts)
            if let Some(&(o_idx, _)) = adj[i].iter().find(|&&(k, o)| mol.atoms[k].z == 8 && o == 2) {
                // Check for anhydride: C(=O)-O-C(=O)
                let anhy_o = adj[i].iter().find(|&&(k, o)| {
                    k != o_idx
                        && mol.atoms[k].z == 8
                        && o == 1
                        && adj[k].iter().any(|&(c2, _)| {
                            c2 != i
                                && mol.atoms[c2].z == 6
                                && adj[c2].iter().any(|&(o3, o3ord)| o3 != k && mol.atoms[o3].z == 8 && o3ord == 2)
                        })
                });
                if let Some(&(bridge_o, _)) = anhy_o {
                    if handled_anhydride_o.insert(bridge_o) {
                        let ring_sz = smallest_ring_size(&adj, i, 7);
                        let (hi_wn, lo_wn, desc) = if matches!(ring_sz, Some(5)) {
                            (1855.0, 1785.0, "Cyclic 5-Ring Anhydride C=O")
                        } else {
                            (1820.0, 1755.0, "Anhydride C=O")
                        };
                        bands.push(IrBand::rule(hi_wn, 0.88, 20.0, "s", format!("{desc} sym stretch"), vec![i, o_idx, bridge_o]));
                        bands.push(IrBand::rule(lo_wn, 0.86, 20.0, "s", format!("{desc} asym stretch"), vec![i, o_idx, bridge_o]));
                    }
                    continue;
                }

                let has_hal = adj[i].iter().any(|&(k, _)| matches!(mol.atoms[k].z, 9 | 17 | 35));
                let has_n = adj[i].iter().any(|&(k, o)| mol.atoms[k].z == 7 && o == 1);
                let has_oh = adj[i].iter().any(|&(k, o)| {
                    k != o_idx && mol.atoms[k].z == 8 && o == 1 && tot_h_of(k) >= 1
                });
                let has_or = adj[i].iter().any(|&(k, o)| {
                    k != o_idx
                        && mol.atoms[k].z == 8
                        && o == 1
                        && adj[k].iter().filter(|&&(m, _)| mol.atoms[m].z > 1).count() >= 2
                });
                let conj = adj[i].iter().any(|&(k, _)| {
                    k != o_idx && (arom.is_atom_aromatic(k) || adj[k].iter().any(|&(_, o2)| o2 == 2))
                });
                let ring_sz = smallest_ring_size(&adj, i, 7);
                let ring_shift = match ring_sz {
                    Some(3) => 100.0,
                    Some(4) => 65.0,
                    Some(5) => 30.0,
                    _ => 0.0,
                };
                let conj_shift = if conj { -22.0 } else { 0.0 };

                let (wn, desc) = if has_hal {
                    (1800.0 + conj_shift, "Acyl Halide C=O stretch")
                } else if has_or {
                    if matches!(ring_sz, Some(4)) {
                        (1820.0 + conj_shift, "β-Lactone C=O stretch (ring-strained)")
                    } else if matches!(ring_sz, Some(5)) {
                        (1772.0 + conj_shift, "γ-Lactone C=O stretch")
                    } else {
                        (1742.0 + conj_shift, "Ester C=O stretch")
                    }
                } else if has_oh {
                    (1712.0 + conj_shift, "Carboxylic Acid C=O stretch")
                } else if has_n {
                    if matches!(ring_sz, Some(4)) {
                        (1755.0 + conj_shift, "β-Lactam C=O stretch (ring-strained)")
                    } else if matches!(ring_sz, Some(5)) {
                        (1700.0 + conj_shift, "γ-Lactam C=O stretch")
                    } else {
                        (1665.0 + conj_shift, "Amide I C=O stretch")
                    }
                } else if tot_h >= 1 {
                    (1728.0 + conj_shift, "Aldehyde C=O stretch")
                } else if matches!(ring_sz, Some(4)) {
                    (1780.0 + conj_shift, "Cyclobutanone C=O stretch (ring-strained)")
                } else if matches!(ring_sz, Some(5)) {
                    (1745.0 + conj_shift, "Cyclopentanone C=O stretch")
                } else {
                    (1715.0 + conj_shift + ring_shift, "Ketone C=O stretch")
                };

                bands.push(IrBand::rule(wn, 0.94, 20.0, "vs", desc, vec![i, o_idx]));

                if tot_h >= 1 && !has_hal && !has_n && !has_oh && !has_or {
                    // Aldehyde Fermi doublet (C-H stretch)
                    bands.push(IrBand::rule(2825.0, 0.44, 22.0, "m", "Aldehyde C–H Fermi doublet (high)", vec![i]));
                    bands.push(IrBand::rule(2720.0, 0.42, 22.0, "m", "Aldehyde C–H Fermi doublet (low)", vec![i]));
                }
                let has_amide_nh = adj[i].iter().any(|&(k, o)| mol.atoms[k].z == 7 && o == 1 && tot_h_of(k) >= 1);
                if has_amide_nh && !matches!(ring_sz, Some(4)) {
                    bands.push(IrBand::rule(1545.0, 0.58, 24.0, "m", "Amide II N–H bend / C–N stretch", vec![i, o_idx]));
                }
            }
        // 5. Nitro & Azide nitrogen centers
        } else if z == 7 {
            let oxys: Vec<usize> = adj[i]
                .iter()
                .filter(|&&(k, _)| mol.atoms[k].z == 8 && adj[k].len() == 1)
                .map(|&(k, _)| k)
                .collect();
            if oxys.len() >= 2 {
                let mut at = vec![i];
                at.extend(oxys);
                bands.push(IrBand::rule(1530.0, 0.85, 24.0, "s", "Nitro N=O asym stretch", at.clone()));
                bands.push(IrBand::rule(1350.0, 0.76, 22.0, "s", "Nitro N=O sym stretch", at));
            }
            // Azide -N=N+=N- (non-aromatic linear central N with two N neighbors, at least one multiple bond)
            let n_nbrs: Vec<(usize, i32)> = adj[i].iter().copied().filter(|&(k, _)| mol.atoms[k].z == 7).collect();
            if !arom.is_atom_aromatic(i) && n_nbrs.len() == 2 && n_nbrs.iter().any(|&(_, o)| o >= 2) {
                bands.push(IrBand::rule(
                    2135.0,
                    0.90,
                    18.0,
                    "vs",
                    "Azide –N₃ asym stretch",
                    vec![n_nbrs[0].0, i, n_nbrs[1].0],
                ));
            }
        // 6. Sulfur (sulfone/sulfonamide/sulfoxide)
        } else if z == 16 {
            let oxys: Vec<usize> = adj[i]
                .iter()
                .filter(|&&(k, o)| mol.atoms[k].z == 8 && o == 2)
                .map(|&(k, _)| k)
                .collect();
            if oxys.len() >= 2 {
                bands.push(IrBand::rule(1325.0, 0.82, 22.0, "s", "Sulfone / Sulfonamide S=O asym", vec![i, oxys[0], oxys[1]]));
                bands.push(IrBand::rule(1155.0, 0.86, 22.0, "s", "Sulfone / Sulfonamide S=O sym", vec![i, oxys[0], oxys[1]]));
            } else if oxys.len() == 1 {
                bands.push(IrBand::rule(1050.0, 0.78, 22.0, "s", "Sulfoxide S=O stretch", vec![i, oxys[0]]));
            }
        // 7. Phosphorus (phosphine oxide / phosphate P=O)
        } else if z == 15 {
            if let Some(&(o_idx, _)) = adj[i].iter().find(|&&(k, o)| mol.atoms[k].z == 8 && o == 2) {
                bands.push(IrBand::rule(1250.0, 0.84, 24.0, "s", "Phosphoryl P=O stretch", vec![i, o_idx]));
            }
        }
    }

    // Bond-based bands (C#N, C#C, C=C + alkene OOP wag, C-O, Si-O, C-F, C-Cl, C-Br)
    for b in &mol.bonds {
        if b.begin_atom >= n || b.end_atom >= n {
            continue;
        }
        let (u, v) = (b.begin_atom, b.end_atom);
        let (zu, zv) = (mol.atoms[u].z, mol.atoms[v].z);
        if b.order == 3 && ((zu == 6 && zv == 7) || (zu == 7 && zv == 6)) {
            bands.push(IrBand::rule(2235.0, 0.65, 16.0, "m", "Nitrile C≡N stretch", vec![u, v]));
        } else if b.order == 3 && zu == 6 && zv == 6 {
            bands.push(IrBand::rule(2140.0, 0.45, 18.0, "w", "Alkyne C≡C stretch", vec![u, v]));
        } else if b.order == 2 && zu == 6 && zv == 6 && !arom.is_atom_aromatic(u) && !arom.is_atom_aromatic(v) {
            bands.push(IrBand::rule(1648.0, 0.52, 20.0, "m", "Alkene C=C stretch", vec![u, v]));
            let hu = tot_h_of(u);
            let hv = tot_h_of(v);
            if hu >= 1 || hv >= 1 {
                bands.push(IrBand::rule(3080.0, 0.34, 22.0, "w", "Alkene =C–H stretch", vec![u, v]));
            }
            if (hu == 2 && hv == 1) || (hu == 1 && hv == 2) {
                bands.push(IrBand::rule(990.0, 0.58, 18.0, "m", "Vinyl =C–H OOP bend (990)", vec![u, v]));
                bands.push(IrBand::rule(910.0, 0.64, 18.0, "s", "Vinyl =CH₂ OOP wag (910)", vec![u, v]));
            } else if (hu == 2 && hv == 0) || (hu == 0 && hv == 2) {
                bands.push(IrBand::rule(890.0, 0.64, 18.0, "s", "Geminal =CH₂ OOP wag", vec![u, v]));
            } else if hu == 1 && hv == 1 {
                bands.push(IrBand::rule(968.0, 0.60, 18.0, "s", "Disubstituted Alkene =C–H OOP bend", vec![u, v]));
            }
        } else if b.order == 1 && ((zu == 6 && zv == 8) || (zu == 8 && zv == 6)) {
            bands.push(IrBand::rule(1185.0, 0.64, 28.0, "s", "C–O stretch (ether/ester/alcohol)", vec![u, v]));
        } else if b.order == 1 && ((zu == 14 && zv == 8) || (zu == 8 && zv == 14)) {
            bands.push(IrBand::rule(1080.0, 0.86, 30.0, "vs", "Siloxane / Silane Si–O stretch", vec![u, v]));
        } else if b.order == 1 && ((zu == 6 && zv == 9) || (zu == 9 && zv == 6)) {
            bands.push(IrBand::rule(1140.0, 0.74, 26.0, "s", "C–F stretch", vec![u, v]));
        } else if b.order == 1 && ((zu == 6 && zv == 17) || (zu == 17 && zv == 6)) {
            bands.push(IrBand::rule(745.0, 0.60, 22.0, "m", "C–Cl stretch", vec![u, v]));
        } else if b.order == 1 && ((zu == 6 && zv == 35) || (zu == 35 && zv == 6)) {
            bands.push(IrBand::rule(620.0, 0.55, 22.0, "m", "C–Br stretch", vec![u, v]));
        }
    }

    if has_arom_ch {
        bands.push(IrBand::rule(
            3055.0,
            0.36,
            24.0,
            "w",
            "Aromatic C–H stretch",
            arom_atoms.iter().take(4).copied().collect(),
        ));
    }
    if has_sp3_ch {
        bands.push(IrBand::rule(2935.0, 0.55, 35.0, "m", "Aliphatic sp³ C–H stretch", sp3_c_atoms));
    }
    if has_ch2 {
        bands.push(IrBand::rule(1465.0, 0.44, 22.0, "m", "Aliphatic CH₂ scissoring bend", ch2_atoms));
    }
    if has_ch3 {
        bands.push(IrBand::rule(1375.0, 0.40, 20.0, "m", "Aliphatic CH₃ umbrella bend", ch3_atoms));
    }
    if has_arom_ring {
        bands.push(IrBand::rule(1595.0, 0.52, 20.0, "m", "Aromatic ring C=C stretch", arom_atoms.clone()));
        bands.push(IrBand::rule(1490.0, 0.50, 20.0, "m", "Aromatic ring breathing", arom_atoms.clone()));

        // Refine benzene 6-ring OOP substitution pattern when a single 6-carbon aromatic ring is present
        let sub_positions: Vec<usize> = arom_atoms
            .iter()
            .copied()
            .filter(|&a| {
                adj[a]
                    .iter()
                    .any(|&(nb, _)| mol.atoms[nb].z != 1 && !arom.is_atom_aromatic(nb))
            })
            .collect();
        if arom_atoms.len() == 6 && sub_positions.len() == 1 {
            bands.push(IrBand::rule(750.0, 0.68, 18.0, "s", "Monosubstituted Benzene C–H OOP bend", arom_atoms.clone()));
            bands.push(IrBand::rule(695.0, 0.64, 18.0, "s", "Monosubstituted Benzene ring puckering", arom_atoms));
        } else if arom_atoms.len() == 6 && sub_positions.len() == 2 {
            // Determine ortho (dist 1), meta (dist 2), or para (dist 3) within the ring
            let (s0, s1) = (sub_positions[0], sub_positions[1]);
            let mut dist = vec![usize::MAX; n];
            let mut q = std::collections::VecDeque::new();
            dist[s0] = 0;
            q.push_back(s0);
            while let Some(u) = q.pop_front() {
                for &(v, _) in &adj[u] {
                    if arom.is_atom_aromatic(v) && dist[v] == usize::MAX {
                        dist[v] = dist[u] + 1;
                        q.push_back(v);
                    }
                }
            }
            match dist[s1] {
                1 => bands.push(IrBand::rule(750.0, 0.66, 20.0, "s", "ortho-Disubstituted Benzene C–H OOP bend", arom_atoms)),
                2 => {
                    bands.push(IrBand::rule(780.0, 0.64, 20.0, "s", "meta-Disubstituted Benzene C–H OOP bend", arom_atoms.clone()));
                    bands.push(IrBand::rule(690.0, 0.58, 18.0, "m", "meta-Disubstituted Benzene ring bend", arom_atoms));
                }
                3 => bands.push(IrBand::rule(820.0, 0.68, 20.0, "s", "para-Disubstituted Benzene C–H OOP bend", arom_atoms)),
                _ => bands.push(IrBand::rule(755.0, 0.62, 24.0, "s", "Aromatic C–H out-of-plane bend", arom_atoms)),
            }
        } else if has_arom_ch {
            bands.push(IrBand::rule(755.0, 0.62, 24.0, "s", "Aromatic C–H out-of-plane bend", arom_atoms));
        }
    }

    // Deduplicate closely spaced identical assignments and sort high -> low cm^-1
    bands.sort_by(|a, b| b.wavenumber_cm1.partial_cmp(&a.wavenumber_cm1).unwrap_or(std::cmp::Ordering::Equal));
    let mut dedup: Vec<IrBand> = Vec::new();
    for b in bands {
        if !dedup.iter().any(|d| d.assignment == b.assignment && (d.wavenumber_cm1 - b.wavenumber_cm1).abs() < 12.0) {
            dedup.push(b);
        }
    }
    let mut pred = IrPrediction {
        source: IrSpectrumSource::RuleBased,
        bands: dedup,
        imaginary_count: 0,
        default_scale_factor: 1.0,
        report_text: String::new(),
        mode_geometry_3d: None,
    };
    pred.report_text = pred.format_report(1.0);
    pred
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

    // Find largest connected conjugated π-system (aromatic + sp2/sp atoms, excluding tetrahedral sulfone/phosphoryl S/P)
    let is_pi: Vec<bool> = (0..n)
        .map(|i| {
            if arom.is_atom_aromatic(i) {
                true
            } else if matches!(mol.atoms[i].z, 15 | 16) && adj[i].len() >= 3 {
                false
            } else {
                adj[i].iter().any(|&(_, o)| matches!(o, 2 | 3 | 6))
            }
        })
        .collect();
    let mut vis = vec![false; n];
    let mut max_pi_comp = 0usize;
    let mut max_arom_comp = 0usize;
    let mut best_comp: Vec<usize> = Vec::new();
    let mut has_any_carbonyl = false;

    for b in &mol.bonds {
        if b.begin_atom < n && b.end_atom < n {
            let (zu, zv) = (mol.atoms[b.begin_atom].z, mol.atoms[b.end_atom].z);
            if b.order == 2 && ((zu == 6 && zv == 8) || (zu == 8 && zv == 6)) {
                has_any_carbonyl = true;
            }
        }
    }

    for i in 0..n {
        if !is_pi[i] || vis[i] {
            continue;
        }
        let mut stack = vec![i];
        vis[i] = true;
        let mut comp = Vec::new();
        let mut acnt = 0usize;
        while let Some(u) = stack.pop() {
            comp.push(u);
            if arom.is_atom_aromatic(u) { acnt += 1; }
            for &(v, _) in &adj[u] {
                if is_pi[v] && !vis[v] {
                    vis[v] = true;
                    stack.push(v);
                }
            }
        }
        if comp.len() > max_pi_comp || (comp.len() == max_pi_comp && acnt > max_arom_comp) {
            max_pi_comp = comp.len();
            max_arom_comp = acnt;
            best_comp = comp;
        }
    }

    let in_best: BTreeSet<usize> = best_comp.iter().copied().collect();
    let mut has_azo = false;
    let mut has_nitro = false;
    let mut has_donor = false;
    let mut has_carbonyl = has_any_carbonyl;
    if !in_best.is_empty() {
        has_carbonyl = false;
        for b in &mol.bonds {
            if b.begin_atom < n && b.end_atom < n && in_best.contains(&b.begin_atom) && in_best.contains(&b.end_atom) {
                let (zu, zv) = (mol.atoms[b.begin_atom].z, mol.atoms[b.end_atom].z);
                if b.order == 2 && zu == 7 && zv == 7 {
                    has_azo = true;
                }
                if b.order == 2 && ((zu == 6 && zv == 8) || (zu == 8 && zv == 6)) {
                    has_carbonyl = true;
                }
            }
        }
        for &u in &best_comp {
            if mol.atoms[u].z == 7 && adj[u].iter().filter(|&&(k, _)| mol.atoms[k].z == 8).count() >= 2 {
                has_nitro = true;
            }
        }
        for i in 0..n {
            // Lone-pair donor (NR2/OR) only: π-atoms (nitro N/O, pyridine N, azo N, C=O) are not donors.
            if matches!(mol.atoms[i].z, 7 | 8) && !is_pi[i] {
                // Donor must be attached directly to an aromatic/alkene carbon of best_comp (NOT a carboxylic -OH on C=O)
                let attached_to_arom_or_alkene = adj[i].iter().any(|&(k, o)| {
                    if o != 1 || !in_best.contains(&k) || mol.atoms[k].z != 6 {
                        return false;
                    }
                    let is_carbonyl_c = mol.atoms[k].z == 6 && adj[k].iter().any(|&(m, mo)| mol.atoms[m].z == 8 && mo == 2);
                    !is_carbonyl_c
                });
                if attached_to_arom_or_alkene {
                    has_donor = true;
                }
            }
        }
        if !has_carbonyl {
            has_carbonyl = has_any_carbonyl;
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

    let max_h_ppm = h1.peaks.iter().map(|p| p.shift_ppm).fold(10.5_f64, f64::max).ceil().max(10.5);
    let max_c_ppm = c13.peaks.iter().map(|p| p.shift_ppm).fold(210.0_f64, f64::max).ceil().max(210.0);

    let mut peaks = Vec::new();
    match kind {
        Nmr2dKind::Cosy => {
            // Diagonal peaks + vicinal (bonded C-C) 3J H-H off-diagonal cross-peaks (exclude exchangeable OH/NH/SH)
            for (&ai, &sh) in &h_by_atom {
                if mol.atoms.get(ai).map(|a| a.z == 6).unwrap_or(false) {
                    peaks.push(Nmr2dCrossPeak {
                        x_ppm: sh,
                        y_ppm: sh,
                        atom_x: ai,
                        atom_y: ai,
                        label: format!("Diag #{} ({:.2} ppm)", ai + 1, sh),
                    });
                }
            }
            for b in &mol.bonds {
                if b.begin_atom < n && b.end_atom < n && b.begin_atom != b.end_atom
                    && mol.atoms[b.begin_atom].z == 6 && mol.atoms[b.end_atom].z == 6
                {
                    if let (Some(&s1), Some(&s2)) = (h_by_atom.get(&b.begin_atom), h_by_atom.get(&b.end_atom)) {
                        if (s1 - s2).abs() > 1e-3 {
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
            }
            Nmr2dPrediction { kind, x_range: (0.0, max_h_ppm), y_range: (0.0, max_h_ppm), peaks }
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
            Nmr2dPrediction { kind, x_range: (0.0, max_h_ppm), y_range: (0.0, max_c_ppm), peaks }
        }
        Nmr2dKind::Hmbc => {
            // 2-bond (all distance 1 first) then 3-bond (distance 2) ¹H–¹³C correlations on C-H protons
            let mut seen = BTreeSet::new();
            for (&ha, &sh) in &h_by_atom {
                if !mol.atoms.get(ha).map(|a| a.z == 6).unwrap_or(false) {
                    continue;
                }
                // Pass 1: Distance 1 (2J_CH)
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
                }
                // Pass 2: Distance 2 (3J_CH)
                for &(n1, _) in &adj[ha] {
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
            Nmr2dPrediction { kind, x_range: (0.0, max_h_ppm), y_range: (0.0, max_c_ppm), peaks }
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

    #[test]
    fn test_expanded_rule_ir_ring_strain_and_hqi_comparison() {
        // Cyclobutanone: 4-membered ring ketone should shift C=O up to ~1780 cm^-1
        let mut cb = ChemMolecule::new();
        for id in 1..=4 {
            cb.add_atom(ChemAtom::new(id, 6));
        }
        cb.add_atom(ChemAtom::new(5, 8));
        cb.add_bond(10, 0, 1, 1, BondDisplay::Solid);
        cb.add_bond(11, 1, 2, 1, BondDisplay::Solid);
        cb.add_bond(12, 2, 3, 1, BondDisplay::Solid);
        cb.add_bond(13, 3, 0, 1, BondDisplay::Solid);
        cb.add_bond(14, 0, 4, 2, BondDisplay::Solid);

        let ir = predict_ir_spectrum(&cb);
        let co = ir.bands.iter().find(|b| b.assignment.contains("C=O")).expect("C=O band");
        assert!(co.wavenumber_cm1 >= 1770.0, "Cyclobutanone C=O expected >=1770, got {}", co.wavenumber_cm1);

        let params = IrRenderParams::default();
        let abs_curve = ir.simulate_absorbance(&params, 400);
        let self_cmp = compare_ir_spectra(&abs_curve, &abs_curve, 400.0, 4000.0);
        assert!(self_cmp.hqi_pct > 99.9, "Self HQI must be ~100%, got {}", self_cmp.hqi_pct);
        assert!(self_cmp.pearson_r > 0.999);

        let peaks = pick_ir_peaks(&abs_curve, 0.1, 10);
        assert!(peaks.iter().any(|p| (p.wavenumber_cm1 - 1780.0).abs() < 25.0));

        let jdx = ir.to_jcamp_dx("Cyclobutanone", &params, 120);
        assert!(jdx.contains("##DATA TYPE=INFRARED SPECTRUM") && jdx.contains("##XYDATA=(X++(Y..Y))"));

        let swapped = IrRenderParams { x_min_cm1: params.x_max_cm1, x_max_cm1: params.x_min_cm1, ..params.clone() };
        let sw = ir.simulate_absorbance(&swapped, 400);
        assert!((sw[0].0 - sw[sw.len() - 1].0).abs() > 1000.0, "swapped range must not collapse");
    }

    fn benzene_with(subs: &[(usize, u8)]) -> ChemMolecule {
        // Kekulé benzene; each (ring_idx, z) adds a substituent atom; z=7 with nitro flag via extra Os below.
        let mut m = ChemMolecule::new();
        for id in 0..6 {
            m.add_atom(ChemAtom::new(id + 1, 6));
        }
        for i in 0..6 {
            m.add_bond(20 + i as u32, i, (i + 1) % 6, if i % 2 == 0 { 2 } else { 1 }, BondDisplay::Solid);
        }
        let mut id = 100u32;
        for &(r, kind) in subs {
            let nidx = m.add_atom(ChemAtom::new(id, 7));
            m.add_bond(id + 500, r, nidx, 1, BondDisplay::Solid);
            id += 1;
            if kind == 1 {
                // nitro: N(=O)O
                for order in [2, 1] {
                    let o = m.add_atom(ChemAtom::new(id, 8));
                    m.add_bond(id + 500, nidx, o, order, BondDisplay::Solid);
                    id += 1;
                }
            }
        }
        m
    }

    #[test]
    fn test_uv_nitro_group_is_acceptor_not_donor() {
        let nitro = predict_uv_vis(&benzene_with(&[(0, 1)]));
        let pna = predict_uv_vis(&benzene_with(&[(0, 1), (3, 0)]));
        assert!(!nitro.chromophore.contains("Push"), "nitrobenzene is not push–pull: {}", nitro.chromophore);
        assert!(pna.chromophore.contains("Push"), "p-nitroaniline is push–pull: {}", pna.chromophore);
        assert!(nitro.lambda_max_nm < pna.lambda_max_nm);
    }
}
