//! Pure-Rust JCAMP-DX (`.jdx` / `.dx`) and Lab CSV/DPT Experimental Spectrum Parser.
//!
//! Supports:
//! - Standard IUPAC JCAMP-DX 4.24 / 5.00 / 6.00 `##XYDATA=(X++(Y..Y))`, `##PEAK TABLE`, `##XYPOINTS`
//! - Full **ASDF** decompression: AFFN, PAC (packed signs), SQZ (`@A-Ia-i`), DIF (`%J-Rj-r`), DUP (`S-Zs`)
//! - **Infrared (FT-IR)** spectra (`1/CM`, `MICROMETERS` → `cm⁻¹`, `%TRANSMITTANCE` / `TRANSMITTANCE` ↔ `ABSORBANCE`)
//! - **NMR** spectra (`PPM`, `HZ` → `ppm` via `##OBSERVE FREQUENCY`)
//! - 2-column instrument exports (`.csv`, `.tsv`, `.dpt`, `.txt` from Bruker OPUS, Thermo OMNIC, PerkinElmer)

use chiralcraft_chem::{absorbance_to_transmittance_pct, transmittance_pct_to_absorbance};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JcampDataType {
    Infrared,
    Nmr,
    UvVis,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JcampXUnit {
    Cm1,
    Micrometers,
    Ppm,
    Hz,
    Nanometers,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JcampYUnit {
    Transmittance,    // 0..1 (or 0..100 auto-detected)
    PctTransmittance, // 0..100
    Absorbance,       // 0..~3+
    Arbitrary,
}

#[derive(Debug, Clone, PartialEq)]
pub struct JcampSpectrum {
    pub title: String,
    pub nucleus: String,
    pub observe_freq_mhz: Option<f64>,
    pub data_type: JcampDataType,
    pub x_units: JcampXUnit,
    pub y_units: JcampYUnit,
    /// Normalized `(x, 0.0..=1.0)` sorted by ascending `x` (used by NMR overlay).
    pub points: Vec<(f64, f64)>,
    /// Physical `(x, y)` after `XFACTOR`/`YFACTOR` and `MICROMETERS -> cm^-1` / `HZ -> ppm` conversion, sorted ascending by `x`.
    pub raw_points: Vec<(f64, f64)>,
}

impl JcampSpectrum {
    pub fn is_infrared(&self) -> bool {
        if self.data_type == JcampDataType::Infrared {
            return true;
        }
        if matches!(self.x_units, JcampXUnit::Cm1 | JcampXUnit::Micrometers) {
            return true;
        }
        // Heuristic fallback when headers omit DATA TYPE/XUNITS: x spans IR wavenumber range
        if let (Some(first), Some(last)) = (self.raw_points.first(), self.raw_points.last()) {
            return last.0 > 400.0 && last.0 <= 6000.0 && first.0 >= 100.0 && (last.0 - first.0) > 250.0;
        }
        false
    }

    /// Returns `(wavenumber_cm1, transmittance_pct 0..100)` sorted by ascending wavenumber.
    pub fn ir_transmittance_curve(&self) -> Vec<(f64, f64)> {
        if self.raw_points.is_empty() {
            return Vec::new();
        }
        let max_y = self
            .raw_points
            .iter()
            .map(|&(_, y)| y)
            .fold(f64::NEG_INFINITY, f64::max);
        let min_y = self
            .raw_points
            .iter()
            .map(|&(_, y)| y)
            .fold(f64::INFINITY, f64::min);

        match self.y_units {
            JcampYUnit::PctTransmittance => self
                .raw_points
                .iter()
                .map(|&(x, y)| (x, y.clamp(0.0, 100.0)))
                .collect(),
            JcampYUnit::Transmittance => {
                let scale = if max_y <= 1.5 { 100.0 } else { 1.0 };
                self.raw_points
                    .iter()
                    .map(|&(x, y)| (x, (y * scale).clamp(0.0, 100.0)))
                    .collect()
            }
            JcampYUnit::Absorbance => self
                .raw_points
                .iter()
                .map(|&(x, y)| (x, absorbance_to_transmittance_pct((y - min_y.min(0.0)).max(0.0))))
                .collect(),
            JcampYUnit::Arbitrary => {
                // Infer from distribution: if median > 50% of max and dips go down, treat as %T; else Absorbance
                if max_y > 1.5 && max_y <= 105.0 && min_y >= 0.0 {
                    self.raw_points
                        .iter()
                        .map(|&(x, y)| (x, y.clamp(0.0, 100.0)))
                        .collect()
                } else if max_y <= 1.5 && min_y >= 0.0 {
                    let mut ys: Vec<f64> = self.raw_points.iter().map(|&(_, y)| y).collect();
                    ys.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                    let med = ys[ys.len() / 2];
                    if med > 0.55 {
                        self.raw_points
                            .iter()
                            .map(|&(x, y)| (x, (y * 100.0).clamp(0.0, 100.0)))
                            .collect()
                    } else {
                        self.raw_points
                            .iter()
                            .map(|&(x, y)| (x, absorbance_to_transmittance_pct(y.max(0.0))))
                            .collect()
                    }
                } else {
                    let span = (max_y - min_y).max(1e-9);
                    self.raw_points
                        .iter()
                        .map(|&(x, y)| (x, absorbance_to_transmittance_pct(((y - min_y) / span) * 1.5)))
                        .collect()
                }
            }
        }
    }

    /// Returns `(wavenumber_cm1, absorbance >= 0)` sorted by ascending wavenumber.
    pub fn ir_absorbance_curve(&self) -> Vec<(f64, f64)> {
        if self.raw_points.is_empty() {
            return Vec::new();
        }
        match self.y_units {
            JcampYUnit::Absorbance => {
                let min_y = self
                    .raw_points
                    .iter()
                    .map(|&(_, y)| y)
                    .fold(f64::INFINITY, f64::min);
                let shift = if min_y < 0.0 { -min_y } else { 0.0 };
                self.raw_points
                    .iter()
                    .map(|&(x, y)| (x, (y + shift).max(0.0)))
                    .collect()
            }
            _ => self
                .ir_transmittance_curve()
                .into_iter()
                .map(|(x, t_pct)| (x, transmittance_pct_to_absorbance(t_pct)))
                .collect(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum AsdfToken {
    Abs(f64),
    Dif(f64),
    Dup(usize),
}

fn sqz_digit(ch: char) -> Option<(f64, u8)> {
    match ch {
        '@' => Some((1.0, 0)),
        'A'..='I' => Some((1.0, (ch as u8) - b'A' + 1)),
        'a'..='i' => Some((-1.0, (ch as u8) - b'a' + 1)),
        _ => None,
    }
}

fn dif_digit(ch: char) -> Option<(f64, u8)> {
    match ch {
        '%' => Some((1.0, 0)),
        'J'..='R' => Some((1.0, (ch as u8) - b'J' + 1)),
        'j'..='r' => Some((-1.0, (ch as u8) - b'j' + 1)),
        _ => None,
    }
}

fn dup_digit(ch: char) -> Option<usize> {
    match ch {
        'S'..='Z' => Some(((ch as u8) - b'S' + 1) as usize),
        's' => Some(9),
        _ => None,
    }
}

/// Decodes one JCAMP-DX `##XYDATA=(X++(Y..Y))` line supporting AFFN, PAC, SQZ, DIF, and DUP encodings.
/// Returns `(x_start_raw, y_values_raw, used_dif_mode)`.
fn decode_asdf_line(line: &str) -> Option<(f64, Vec<f64>, bool)> {
    let body = line.split("$$").next().unwrap_or(line).trim();
    if body.is_empty() {
        return None;
    }

    let chars: Vec<char> = body.chars().collect();
    let len = chars.len();
    let mut i = 0usize;
    let mut tokens: Vec<AsdfToken> = Vec::new();
    let mut used_dif = false;

    while i < len {
        let c = chars[i];
        if c.is_whitespace() || c == ',' || c == ';' {
            i += 1;
            continue;
        }

        if let Some((sign, first_d)) = sqz_digit(c) {
            i += 1;
            let mut buf = String::new();
            if sign < 0.0 {
                buf.push('-');
            }
            buf.push((b'0' + first_d) as char);
            while i < len && (chars[i].is_ascii_digit() || chars[i] == '.') {
                buf.push(chars[i]);
                i += 1;
            }
            if let Ok(v) = buf.parse::<f64>() {
                tokens.push(AsdfToken::Abs(v));
            }
            continue;
        }

        if let Some((sign, first_d)) = dif_digit(c) {
            used_dif = true;
            i += 1;
            let mut buf = String::new();
            if sign < 0.0 {
                buf.push('-');
            }
            buf.push((b'0' + first_d) as char);
            while i < len && (chars[i].is_ascii_digit() || chars[i] == '.') {
                buf.push(chars[i]);
                i += 1;
            }
            if let Ok(v) = buf.parse::<f64>() {
                tokens.push(AsdfToken::Dif(v));
            }
            continue;
        }

        if let Some(first_d) = dup_digit(c) {
            i += 1;
            let mut count = first_d;
            while i < len && chars[i].is_ascii_digit() {
                count = count.saturating_mul(10).saturating_add((chars[i] as u8 - b'0') as usize);
                i += 1;
            }
            if count >= 2 {
                tokens.push(AsdfToken::Dup(count));
            }
            continue;
        }

        // Standard AFFN / PAC number (starts with digit, '+', '-', or '.')
        if c.is_ascii_digit() || c == '+' || c == '-' || c == '.' {
            let start = i;
            let mut has_dot = c == '.';
            let mut has_exp = false;
            i += 1;
            while i < len {
                let ch = chars[i];
                if ch.is_ascii_digit() {
                    i += 1;
                } else if ch == '.' && !has_dot && !has_exp {
                    has_dot = true;
                    i += 1;
                } else if (ch == 'E' || ch == 'e')
                    && !has_exp
                    && i + 1 < len
                    && (chars[i + 1] == '+' || chars[i + 1] == '-')
                {
                    // Scientific notation only when explicit sign (E+/E-);
                    // otherwise 'E'/'e' followed by a digit (e.g. "1700.0E0K0") is a JCAMP SQZ pseudo-digit (+5 / -5)!
                    has_exp = true;
                    i += 2;
                } else {
                    break;
                }
            }
            let s: String = chars[start..i].iter().collect();
            if let Ok(v) = s.parse::<f64>() {
                if v.is_finite() {
                    tokens.push(AsdfToken::Abs(v));
                }
            }
            continue;
        }

        // Skip unrecognized separator/noise character
        i += 1;
    }

    if tokens.len() < 2 {
        return None;
    }

    let x_start = match tokens[0] {
        AsdfToken::Abs(v) => v,
        _ => return None,
    };

    let mut ys: Vec<f64> = Vec::new();
    let mut last_dif: Option<f64> = None;
    for &tok in &tokens[1..] {
        match tok {
            AsdfToken::Abs(v) => {
                ys.push(v);
                last_dif = None;
            }
            AsdfToken::Dif(d) => {
                let prev = ys.last().copied().unwrap_or(0.0);
                ys.push(prev + d);
                last_dif = Some(d);
            }
            AsdfToken::Dup(count) => {
                let repeats = count.saturating_sub(1).min(4096);
                if let Some(d) = last_dif {
                    for _ in 0..repeats {
                        let prev = ys.last().copied().unwrap_or(0.0);
                        ys.push(prev + d);
                    }
                } else if let Some(&val) = ys.last() {
                    for _ in 0..repeats {
                        ys.push(val);
                    }
                }
            }
        }
    }

    if ys.is_empty() {
        None
    } else {
        Some((x_start, ys, used_dif))
    }
}

/// Parses a JCAMP-DX (`.jdx` / `.dx`) text file (IR, NMR, or UV-Vis, including ASDF/SQZ/DIF/DUP compression).
pub fn parse_jcamp_dx(text: &str) -> Result<JcampSpectrum, String> {
    let mut title = String::from("Experimental JCAMP-DX");
    let mut nucleus = String::from("1H");
    let mut obs_freq: Option<f64> = None;
    let mut data_type = JcampDataType::Unknown;
    let mut x_units = JcampXUnit::Unknown;
    let mut y_units = JcampYUnit::Arbitrary;
    let mut first_x: Option<f64> = None;
    let mut last_x: Option<f64> = None;
    let mut delta_x: Option<f64> = None;
    let mut x_factor = 1.0_f64;
    let mut y_factor = 1.0_f64;
    let mut npoints: Option<usize> = None;

    let mut in_xydata = false;
    let mut in_peak_table = false;
    let mut raw_xy: Vec<(f64, f64)> = Vec::new();
    let mut row_ys: Vec<(f64, Vec<f64>, bool)> = Vec::new();

    for line in text.lines() {
        let trim = line.trim();
        if trim.is_empty() || trim.starts_with("$$") {
            continue;
        }
        if trim.starts_with("##") {
            in_xydata = false;
            in_peak_table = false;
            let upper = trim.to_ascii_uppercase();
            let val = trim.split_once('=').map(|(_, v)| v.trim()).unwrap_or("");
            let val_up = val.to_ascii_uppercase();
            if upper.starts_with("##TITLE=") {
                if !val.is_empty() {
                    title = val.to_string();
                }
            } else if upper.starts_with("##DATA TYPE=") {
                if val_up.contains("INFRARED") || val_up.contains("IR") {
                    data_type = JcampDataType::Infrared;
                } else if val_up.contains("NMR") {
                    data_type = JcampDataType::Nmr;
                } else if val_up.contains("UV") {
                    data_type = JcampDataType::UvVis;
                }
            } else if upper.starts_with("##.OBSERVE NUCLEUS=") || upper.starts_with("##NUCLEUS=") {
                let clean = val.trim_start_matches('^').trim();
                if !clean.is_empty() {
                    nucleus = clean.to_string();
                }
            } else if upper.starts_with("##.OBSERVE FREQUENCY=") || upper.starts_with("##SF=") {
                obs_freq = val.split_whitespace().next().and_then(|s| s.parse().ok());
            } else if upper.starts_with("##XUNITS=") {
                if val_up.contains("1/CM") || val_up.contains("CM-1") || val_up.contains("WAVENUMBER") {
                    x_units = JcampXUnit::Cm1;
                } else if val_up.contains("MICROMETER") || val_up.contains("MICRON") || val_up == "UM" {
                    x_units = JcampXUnit::Micrometers;
                } else if val_up.contains("HZ") {
                    x_units = JcampXUnit::Hz;
                } else if val_up.contains("PPM") {
                    x_units = JcampXUnit::Ppm;
                } else if val_up.contains("NANOMETER") || val_up == "NM" {
                    x_units = JcampXUnit::Nanometers;
                }
            } else if upper.starts_with("##YUNITS=") {
                if val_up.contains("ABSORB") {
                    y_units = JcampYUnit::Absorbance;
                } else if val_up.contains("%") && val_up.contains("TRANS") {
                    y_units = JcampYUnit::PctTransmittance;
                } else if val_up.contains("TRANS") {
                    y_units = JcampYUnit::Transmittance;
                }
            } else if upper.starts_with("##FIRSTX=") {
                first_x = val.split_whitespace().next().and_then(|s| s.parse().ok());
            } else if upper.starts_with("##LASTX=") {
                last_x = val.split_whitespace().next().and_then(|s| s.parse().ok());
            } else if upper.starts_with("##DELTAX=") {
                delta_x = val.split_whitespace().next().and_then(|s| s.parse().ok());
            } else if upper.starts_with("##XFACTOR=") {
                if let Some(v) = val.split_whitespace().next().and_then(|s| s.parse().ok()) {
                    x_factor = v;
                }
            } else if upper.starts_with("##YFACTOR=") {
                if let Some(v) = val.split_whitespace().next().and_then(|s| s.parse().ok()) {
                    y_factor = v;
                }
            } else if upper.starts_with("##NPOINTS=") {
                npoints = val.split_whitespace().next().and_then(|s| s.parse().ok());
            } else if upper.starts_with("##END=") {
                if !row_ys.is_empty() || !raw_xy.is_empty() {
                    break;
                }
            } else if upper.starts_with("##XYDATA=") {
                if val_up.contains("XY..XY") {
                    in_peak_table = true;
                } else {
                    in_xydata = true;
                }
            } else if upper.starts_with("##PEAK TABLE=") || upper.starts_with("##XYPOINTS=") {
                if row_ys.is_empty() {
                    in_peak_table = true;
                }
            }
            continue;
        }

        if in_xydata {
            if let Some((rx, ys, used_dif)) = decode_asdf_line(trim) {
                row_ys.push((
                    rx * x_factor,
                    ys.into_iter().map(|y| y * y_factor).collect(),
                    used_dif,
                ));
            }
        } else if in_peak_table {
            let body = trim.split("$$").next().unwrap_or(trim);
            for chunk in body.split(';') {
                let parts: Vec<f64> = chunk
                    .split(|c: char| c == ',' || c.is_whitespace())
                    .filter_map(|t| t.trim().parse::<f64>().ok())
                    .collect();
                for pair in parts.chunks(2) {
                    if pair.len() == 2 {
                        raw_xy.push((pair[0] * x_factor, pair[1] * y_factor));
                    }
                }
            }
        }
    }

    if !row_ys.is_empty() {
        let fx = first_x.unwrap_or_else(|| row_ys.first().map(|r| r.0).unwrap_or(10.0));
        let lx = last_x.unwrap_or_else(|| row_ys.last().map(|r| r.0).unwrap_or(0.0));
        // In JCAMP-DX DIF mode, the last point of line `k` is a checkpoint repeated as the first point of line `k+1`.
        let any_dif = row_ys.iter().any(|r| r.2);
        let effective_total: usize = if any_dif {
            let mut tot = 0usize;
            for (idx, (_, ys, used_dif)) in row_ys.iter().enumerate() {
                let drop_cp = *used_dif && idx + 1 < row_ys.len() && ys.len() > 1;
                tot += ys.len() - usize::from(drop_cp);
            }
            tot
        } else {
            row_ys.iter().map(|(_, ys, _)| ys.len()).sum()
        };
        let denom = (npoints.unwrap_or(effective_total).max(2) - 1) as f64;
        let computed_dx = (lx - fx) / denom;
        let dx = match delta_x {
            Some(d) if d.abs() > 1e-12 => {
                // Honor sign of (LASTX - FIRSTX) in case file writes unsigned DELTAX for descending X!
                if computed_dx < 0.0 && d > 0.0 {
                    -d
                } else if computed_dx > 0.0 && d < 0.0 {
                    -d
                } else {
                    d
                }
            }
            _ => computed_dx,
        };
        for (row_idx, (row_x, ys, used_dif)) in row_ys.iter().enumerate() {
            let drop_checkpoint = *used_dif && row_idx + 1 < row_ys.len() && ys.len() > 1;
            let take_n = ys.len() - usize::from(drop_checkpoint);
            for (k, &y) in ys.iter().take(take_n).enumerate() {
                raw_xy.push((row_x + (k as f64) * dx, y));
            }
        }
    }

    if raw_xy.is_empty() {
        return Err("No XYDATA or PEAK TABLE numeric points found in JCAMP-DX file".to_string());
    }

    // Unit conversions
    if x_units == JcampXUnit::Hz {
        if let Some(sf) = obs_freq.filter(|&f| f > 1.0) {
            for (x, _) in &mut raw_xy {
                *x /= sf;
            }
            x_units = JcampXUnit::Ppm;
        }
    } else if x_units == JcampXUnit::Micrometers {
        raw_xy.retain(|&(x, y)| x.is_finite() && y.is_finite() && x > 0.1);
        for (x, _) in &mut raw_xy {
            *x = 10_000.0 / *x;
        }
        x_units = JcampXUnit::Cm1;
    }
    raw_xy.retain(|&(x, y)| x.is_finite() && y.is_finite());
    if raw_xy.is_empty() {
        return Err("No valid finite XYDATA points remain after unit conversion".to_string());
    }

    raw_xy.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let min_y = raw_xy.iter().map(|&(_, y)| y).fold(f64::INFINITY, f64::min);
    let max_y = raw_xy.iter().map(|&(_, y)| y).fold(f64::NEG_INFINITY, f64::max);
    let y_span = (max_y - min_y).max(1e-9);
    let points: Vec<(f64, f64)> = raw_xy
        .iter()
        .map(|&(x, y)| (x, ((y - min_y) / y_span).clamp(0.0, 1.0)))
        .collect();

    Ok(JcampSpectrum {
        title,
        nucleus,
        observe_freq_mhz: obs_freq,
        data_type,
        x_units,
        y_units,
        points,
        raw_points: raw_xy,
    })
}

/// Parses 2-column lab FT-IR `.csv`, `.tsv`, `.dpt`, or `.txt` files (Bruker OPUS, Thermo OMNIC, PerkinElmer).
pub fn parse_ir_csv_dpt(text: &str, title: &str) -> Result<JcampSpectrum, String> {
    let mut x_units = JcampXUnit::Cm1;
    let mut y_units = JcampYUnit::Arbitrary;
    let mut raw_xy: Vec<(f64, f64)> = Vec::new();

    for line in text.lines() {
        let trim = line.trim();
        if trim.is_empty() || trim.starts_with('#') || trim.starts_with("//") {
            let up = trim.to_ascii_uppercase();
            if up.contains("ABSORB") {
                y_units = JcampYUnit::Absorbance;
            } else if up.contains("%T") || up.contains("TRANS") {
                y_units = if up.contains('%') {
                    JcampYUnit::PctTransmittance
                } else {
                    JcampYUnit::Transmittance
                };
            }
            continue;
        }
        let toks: Vec<&str> = trim
            .split(|c: char| c == ',' || c == ';' || c == '\t' || c.is_whitespace())
            .filter(|s| !s.is_empty())
            .collect();
        if toks.len() >= 2 {
            if let (Ok(x), Ok(y)) = (toks[0].parse::<f64>(), toks[1].parse::<f64>()) {
                if x.is_finite() && y.is_finite() {
                    raw_xy.push((x, y));
                }
                continue;
            }
        }
        // Header row inspection
        let up = trim.to_ascii_uppercase();
        if up.contains("ABSORB") || up.contains(",ABS") || up.contains(" ABS") {
            y_units = JcampYUnit::Absorbance;
        } else if up.contains("%T") || up.contains("TRANSMITTANCE_PCT") {
            y_units = JcampYUnit::PctTransmittance;
        } else if up.contains("TRANS") {
            y_units = JcampYUnit::Transmittance;
        }
        if up.contains("MICROMETER") || up.contains("MICRON") || up.contains("(UM)") {
            x_units = JcampXUnit::Micrometers;
        }
    }

    if raw_xy.len() < 4 {
        return Err("No valid 2-column IR wavenumber/intensity pairs found in CSV/DPT file".to_string());
    }
    if x_units == JcampXUnit::Micrometers {
        raw_xy.retain(|&(x, _)| x > 0.1);
        for (x, _) in &mut raw_xy {
            *x = 10_000.0 / *x;
        }
        x_units = JcampXUnit::Cm1;
    }
    if raw_xy.is_empty() {
        return Err("No valid wavenumber/intensity points remain after micrometer conversion".to_string());
    }
    raw_xy.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let min_y = raw_xy.iter().map(|&(_, y)| y).fold(f64::INFINITY, f64::min);
    let max_y = raw_xy.iter().map(|&(_, y)| y).fold(f64::NEG_INFINITY, f64::max);
    let span = (max_y - min_y).max(1e-9);
    let points = raw_xy
        .iter()
        .map(|&(x, y)| (x, ((y - min_y) / span).clamp(0.0, 1.0)))
        .collect();

    Ok(JcampSpectrum {
        title: if title.trim().is_empty() {
            "Experimental IR CSV".to_string()
        } else {
            title.to_string()
        },
        nucleus: String::new(),
        observe_freq_mhz: None,
        data_type: JcampDataType::Infrared,
        x_units,
        y_units,
        points,
        raw_points: raw_xy,
    })
}

/// Auto-detects and parses either JCAMP-DX (`.jdx`/`.dx`) or 2-column lab IR (`.csv`/`.tsv`/`.dpt`/`.txt`).
pub fn parse_experimental_ir_file(path: &str, text: &str) -> Result<JcampSpectrum, String> {
    let low = path.to_ascii_lowercase();
    let trimmed = text.trim_start();
    if low.ends_with(".jdx") || low.ends_with(".dx") || low.ends_with(".jcamp") || trimmed.starts_with("##") {
        parse_jcamp_dx(text)
    } else {
        let title = std::path::Path::new(path)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("Experimental IR");
        parse_ir_csv_dpt(text, title)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jcamp_dx_xydata_and_peak_table_parsing() {
        let sample = "\
##TITLE=Ethanol 1H NMR
##JCAMP-DX=5.00
##.OBSERVE NUCLEUS=^1H
##.OBSERVE FREQUENCY=400.0
##XUNITS=PPM
##FIRSTX=4.0
##LASTX=1.0
##NPOINTS=7
##XYDATA=(X++(Y..Y))
4.0 0.05 0.90 0.10 0.02
2.0 0.03 0.15 1.00
##END=
";
        let spec = parse_jcamp_dx(sample).expect("parse jcamp");
        assert_eq!(spec.title, "Ethanol 1H NMR");
        assert_eq!(spec.nucleus, "1H");
        assert_eq!(spec.points.len(), 7);
        assert!(spec.points.iter().any(|&(_, y)| (y - 1.0).abs() < 1e-6));
    }

    #[test]
    fn test_jcamp_asdf_sqz_dif_dup_and_ir_units() {
        // Real-style NIST compressed IR block with SQZ ('I'=+9, 'E'=+5), DIF ('k'=-2, 'K'=+2), DUP ('T'=2x)
        // Line 1: 1700 I0(+90) k0(-20 -> 70) T(repeat -20 -> 50 checkpoint)
        // Line 2: 1708 E0(+50 checkpoint) K0(+20 -> 70)
        let nist_sample = "\
##TITLE=Acetone IR (NIST ASDF)
##JCAMP-DX=4.24
##DATA TYPE=INFRARED SPECTRUM
##XUNITS=1/CM
##YUNITS=TRANSMITTANCE
##XFACTOR=1.0
##YFACTOR=0.01
##FIRSTX=1700.0
##LASTX=1712.0
##DELTAX=4.0
##NPOINTS=4
##XYDATA=(X++(Y..Y))
1700I0k0T
1708E0K0
##END=
";
        let spec = parse_jcamp_dx(nist_sample).expect("parse ASDF IR");
        assert!(spec.is_infrared());
        assert_eq!(spec.data_type, JcampDataType::Infrared);
        assert_eq!(spec.y_units, JcampYUnit::Transmittance);
        let tr = spec.ir_transmittance_curve();
        assert_eq!(tr.len(), 4); // 1700 (90%), 1704 (70%), 1708 (50%), 1712 (70%)
        assert!((tr[0].0 - 1700.0).abs() < 1e-3 && (tr[0].1 - 90.0).abs() < 1e-3);
        assert!((tr[1].0 - 1704.0).abs() < 1e-3 && (tr[1].1 - 70.0).abs() < 1e-3);
        assert!((tr[2].0 - 1708.0).abs() < 1e-3 && (tr[2].1 - 50.0).abs() < 1e-3);
        assert!((tr[3].0 - 1712.0).abs() < 1e-3 && (tr[3].1 - 70.0).abs() < 1e-3);

        let abs = spec.ir_absorbance_curve();
        // At 50% transmittance, A = -log10(0.5) ≈ 0.30103
        assert!((abs[2].1 - 0.30103).abs() < 1e-3);
    }

    #[test]
    fn test_jcamp_micrometers_and_csv_dpt_parser() {
        let um_jdx = "\
##TITLE=Micron IR
##DATA TYPE=INFRARED SPECTRUM
##XUNITS=MICROMETERS
##YUNITS=ABSORBANCE
##PEAK TABLE=(XY..XY)
5.0,0.85; 10.0,0.40
##END=
";
        let sp = parse_jcamp_dx(um_jdx).expect("parse micron jdx");
        assert_eq!(sp.x_units, JcampXUnit::Cm1);
        assert!((sp.raw_points[0].0 - 1000.0).abs() < 1e-3); // 10000 / 10 um
        assert!((sp.raw_points[1].0 - 2000.0).abs() < 1e-3); // 10000 / 5 um

        let csv = "Wavenumber(cm-1),%Transmittance\n4000,98.5\n3000,72.0\n1715,12.4\n1000,85.0\n500,94.0\n";
        let sp2 = parse_experimental_ir_file("sample.csv", csv).expect("parse csv");
        assert!(sp2.is_infrared());
        assert_eq!(sp2.y_units, JcampYUnit::PctTransmittance);
        assert_eq!(sp2.raw_points.len(), 5);
    }

    #[test]
    fn test_jcamp_descending_unsigned_deltax_pac_scinotation_and_multiblock() {
        // Descending FIRSTX > LASTX with positive DELTAX, scientific notation, PAC signs, and trailing ##END= + second block
        let jdx = "\
##TITLE=Descending PAC IR
##DATA TYPE=INFRARED SPECTRUM
##XUNITS=1/CM
##YUNITS=ABSORBANCE
##FIRSTX=4000.0
##LASTX=3980.0
##DELTAX=10.0
##NPOINTS=3
##XYDATA=(X++(Y..Y))
4000.0 1.5E-02+0.45-0.02
##END=
##TITLE=Second Block Should Be Ignored
##PEAK TABLE=(XY..XY)
1000,9.99
##END=
";
        let sp = parse_jcamp_dx(jdx).expect("parse descending PAC");
        assert_eq!(sp.raw_points.len(), 3);
        // Sorted ascending: 3980 (-0.02), 3990 (0.45), 4000 (0.015)
        assert!((sp.raw_points[0].0 - 3980.0).abs() < 1e-3);
        assert!((sp.raw_points[1].0 - 3990.0).abs() < 1e-3 && (sp.raw_points[1].1 - 0.45).abs() < 1e-4);
        assert!((sp.raw_points[2].0 - 4000.0).abs() < 1e-3 && (sp.raw_points[2].1 - 0.015).abs() < 1e-4);
    }
}
