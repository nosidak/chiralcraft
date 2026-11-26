//! Pure-Rust JCAMP-DX (`.jdx` / `.dx`) Experimental NMR Spectrum Parser.
//!
//! Reads standard IUPAC JCAMP-DX `##XYDATA=(X++(Y..Y))` and `##PEAK TABLE=(XY..XY)` blocks,
//! converts Hz to ppm when `##XUNITS=HZ` + `##OBSERVE FREQUENCY` are present, and normalizes
//! intensities to `[0.0, 1.0]` for direct visual overlay behind ChiralCraft's predicted NMR curve.

#[derive(Debug, Clone, PartialEq)]
pub struct JcampSpectrum {
    pub title: String,
    pub nucleus: String,
    pub observe_freq_mhz: Option<f64>,
    /// `(ppm, normalized_intensity 0.0..=1.0)` sorted by ascending ppm.
    pub points: Vec<(f64, f64)>,
}

/// Parses a JCAMP-DX (`.jdx` / `.dx`) text file into a normalized `(ppm, intensity)` spectrum.
pub fn parse_jcamp_dx(text: &str) -> Result<JcampSpectrum, String> {
    let mut title = String::from("Experimental JCAMP-DX");
    let mut nucleus = String::from("1H");
    let mut obs_freq: Option<f64> = None;
    let mut x_units_hz = false;
    let mut first_x: Option<f64> = None;
    let mut last_x: Option<f64> = None;
    let mut x_factor = 1.0_f64;
    let mut y_factor = 1.0_f64;
    let mut npoints: Option<usize> = None;

    let mut in_xydata = false;
    let mut in_peak_table = false;
    let mut raw_xy: Vec<(f64, f64)> = Vec::new();
    let mut row_ys: Vec<(f64, Vec<f64>)> = Vec::new();

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
            if upper.starts_with("##TITLE=") {
                if !val.is_empty() { title = val.to_string(); }
            } else if upper.starts_with("##.OBSERVE NUCLEUS=") || upper.starts_with("##NUCLEUS=") {
                let clean = val.trim_start_matches('^').trim();
                if !clean.is_empty() { nucleus = clean.to_string(); }
            } else if upper.starts_with("##.OBSERVE FREQUENCY=") || upper.starts_with("##SF=") {
                obs_freq = val.split_whitespace().next().and_then(|s| s.parse().ok());
            } else if upper.starts_with("##XUNITS=") {
                x_units_hz = val.to_ascii_uppercase().contains("HZ");
            } else if upper.starts_with("##FIRSTX=") {
                first_x = val.split_whitespace().next().and_then(|s| s.parse().ok());
            } else if upper.starts_with("##LASTX=") {
                last_x = val.split_whitespace().next().and_then(|s| s.parse().ok());
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
            } else if upper.starts_with("##XYDATA=") {
                in_xydata = true;
            } else if upper.starts_with("##PEAK TABLE=") || upper.starts_with("##XYPOINTS=") {
                in_peak_table = true;
            }
            continue;
        }

        if in_xydata {
            // Strip inline JCAMP comments starting with $$
            let body = trim.split("$$").next().unwrap_or(trim);
            let nums: Vec<f64> = body
                .split(|c: char| c.is_whitespace() || c == ',' || c == ';')
                .filter_map(|t| t.parse::<f64>().ok())
                .collect();
            if nums.len() >= 2 {
                row_ys.push((nums[0] * x_factor, nums[1..].iter().map(|&y| y * y_factor).collect()));
            }
        } else if in_peak_table {
            let body = trim.split("$$").next().unwrap_or(trim);
            for pair in body.split(';').flat_map(|s| s.split_whitespace().collect::<Vec<_>>().chunks(2).map(|c| c.join(",")).collect::<Vec<_>>()) {
                let parts: Vec<f64> = pair
                    .split(',')
                    .filter_map(|t| t.trim().parse::<f64>().ok())
                    .collect();
                if parts.len() >= 2 {
                    raw_xy.push((parts[0] * x_factor, parts[1] * y_factor));
                }
            }
        }
    }

    if !row_ys.is_empty() {
        let total_ys: usize = row_ys.iter().map(|(_, ys)| ys.len()).sum();
        let fx = first_x.unwrap_or_else(|| row_ys.first().map(|r| r.0).unwrap_or(10.0));
        let lx = last_x.unwrap_or_else(|| row_ys.last().map(|r| r.0).unwrap_or(0.0));
        let denom = (npoints.unwrap_or(total_ys).max(2) - 1) as f64;
        let dx = (lx - fx) / denom;
        for (row_x, ys) in row_ys {
            for (k, y) in ys.into_iter().enumerate() {
                raw_xy.push((row_x + (k as f64) * dx, y));
            }
        }
    }

    if raw_xy.is_empty() {
        return Err("No XYDATA or PEAK TABLE numeric points found in JCAMP-DX file".to_string());
    }

    // Convert Hz -> ppm if needed
    if x_units_hz {
        if let Some(sf) = obs_freq.filter(|&f| f > 1.0) {
            for (x, _) in &mut raw_xy {
                *x /= sf;
            }
        }
    }

    let min_y = raw_xy.iter().map(|&(_, y)| y).fold(f64::INFINITY, f64::min);
    let max_y = raw_xy.iter().map(|&(_, y)| y).fold(f64::NEG_INFINITY, f64::max);
    let y_span = (max_y - min_y).max(1e-9);
    let mut points: Vec<(f64, f64)> = raw_xy
        .into_iter()
        .map(|(x, y)| (x, ((y - min_y) / y_span).clamp(0.0, 1.0)))
        .collect();
    points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    Ok(JcampSpectrum {
        title,
        nucleus,
        observe_freq_mhz: obs_freq,
        points,
    })
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
}
