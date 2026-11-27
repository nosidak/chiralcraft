//! Non-blocking Local Quantum Chemistry Runner (`xtb`, `orca`, `g16`)
//! and Optional NIST Chemistry WebBook Experimental IR (`.jdx`) Fetcher.
//!
//! Uses only `std::process::Command` + `std::thread` + `std::sync::mpsc` (zero external async crates),
//! matching ChiralCraft's `chem_db.rs` background polling architecture.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::mpsc::{self, Receiver};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

use chiralcraft_chem::{
    prepare_molecule_for_quantum, total_formal_charge, ChemMolecule, IrPrediction,
};
use chiralcraft_formats::{
    export_gaussian_input, export_orca_input, export_xtb_input, parse_jcamp_dx,
    parse_quantum_output, parse_xtb_workdir, JcampSpectrum, QuantumJobSpec, QuantumJobType,
    QuantumOutputReport, QuantumProgram,
};

static JOB_COUNTER: AtomicU64 = AtomicU64::new(1);

// ============================================================================
// 1. Safe Binary Auto-Detection (Rejects GNOME /usr/bin/orca Screen Reader!)
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct QcBinaryAvailability {
    pub detected: bool,
    pub xtb_path: Option<String>,
    pub xtb_version: String,
    pub orca_path: Option<String>,
    pub orca_version: String,
    pub gaussian_path: Option<String>,
    pub gnome_orca_skipped: bool,
}

fn find_on_path(bin_name: &str) -> Option<PathBuf> {
    let p = Path::new(bin_name);
    if p.is_absolute() && p.is_file() {
        return Some(p.to_path_buf());
    }
    let path_env = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path_env) {
        let cand = dir.join(bin_name);
        if cand.is_file() {
            return Some(cand);
        }
    }
    None
}

/// Validates that `path` is the Max Planck **ORCA quantum chemistry** executable and **NOT**
/// the Linux GNOME Orca screen reader (`/usr/bin/orca` is a Python script installed by default on
/// many Linux desktops).
pub fn is_valid_quantum_orca_binary(path: &Path) -> bool {
    if !path.is_file() {
        return false;
    }
    // Quick header inspection: reject Python scripts (GNOME screen reader starts with #!/usr/bin/python)
    if let Ok(head) = fs::read(path).map(|b| b.into_iter().take(256).collect::<Vec<u8>>()) {
        if let Ok(txt) = std::str::from_utf8(&head) {
            let low = txt.to_ascii_lowercase();
            if low.starts_with("#!") && (low.contains("python") || low.contains("gnome") || low.contains("orca.main")) {
                return false;
            }
        }
    }
    // Confirm it doesn't identify as GNOME screen reader; in tests, mock scripts may echo ORCA banner
    let Ok(out) = Command::new(path)
        .arg("--version")
        .stdin(Stdio::null())
        .output()
    else {
        return false;
    };
    let combined = format!(
        "{}\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    let low = combined.to_ascii_lowercase();
    if low.contains("screen reader") || low.contains("gnome") || low.contains("at-spi") {
        return false;
    }
    // Quantum ORCA prints "Program Version", "* O   R   C   A *", or "This program requires the name of a parameterfile"
    combined.contains("O   R   C   A")
        || low.contains("program version")
        || low.contains("parameterfile")
        || low.contains("orca quantum")
        || low.contains("neese")
}

/// Detects `xtb`, quantum `orca`, and `g16`/`g09` using optional user-configured overrides first.
pub fn detect_qc_binaries(
    xtb_override: &str,
    orca_override: &str,
    gaussian_override: &str,
) -> QcBinaryAvailability {
    let mut avail = QcBinaryAvailability {
        detected: true,
        ..QcBinaryAvailability::default()
    };

    // 1. xTB
    let xtb_cand = if !xtb_override.trim().is_empty() {
        Some(PathBuf::from(xtb_override.trim()))
    } else {
        find_on_path("xtb")
    };
    if let Some(p) = xtb_cand.filter(|p| p.is_file()) {
        if let Ok(out) = Command::new(&p).arg("--version").stdin(Stdio::null()).output() {
            let txt = format!(
                "{}\n{}",
                String::from_utf8_lossy(&out.stdout),
                String::from_utf8_lossy(&out.stderr)
            );
            if txt.to_ascii_lowercase().contains("xtb") || txt.contains("GFN") {
                avail.xtb_version = txt
                    .lines()
                    .find(|l| l.to_ascii_lowercase().contains("version") || l.contains("xtb"))
                    .unwrap_or("xtb")
                    .trim()
                    .to_string();
                avail.xtb_path = Some(p.to_string_lossy().to_string());
            }
        }
    }

    // 2. Quantum ORCA (explicitly skipping GNOME /usr/bin/orca)
    let mut orca_cands: Vec<PathBuf> = Vec::new();
    if !orca_override.trim().is_empty() {
        orca_cands.push(PathBuf::from(orca_override.trim()));
    }
    if let Some(p) = find_on_path("orca") {
        orca_cands.push(p);
    }
    if let Ok(home) = std::env::var("HOME") {
        orca_cands.push(PathBuf::from(&home).join("orca").join("orca"));
    }
    orca_cands.push(PathBuf::from("/opt/orca/orca"));
    for cand in orca_cands {
        if !cand.is_file() {
            continue;
        }
        if is_valid_quantum_orca_binary(&cand) {
            avail.orca_path = Some(cand.to_string_lossy().to_string());
            avail.orca_version = "ORCA (Quantum Chemistry)".to_string();
            break;
        } else if cand == Path::new("/usr/bin/orca") {
            avail.gnome_orca_skipped = true;
        }
    }

    // 3. Gaussian (g16 / g09)
    let g_cand = if !gaussian_override.trim().is_empty() {
        Some(PathBuf::from(gaussian_override.trim()))
    } else {
        find_on_path("g16").or_else(|| find_on_path("g09"))
    };
    if let Some(p) = g_cand.filter(|p| p.is_file()) {
        avail.gaussian_path = Some(p.to_string_lossy().to_string());
    }

    avail
}

// ============================================================================
// 2. Async Local Quantum FTIR Runner (GFN2-xTB / GFN-FF / ORCA / Gaussian)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QcRunBackend {
    XtbGfn2,
    XtbGfnFf,
    OrcaDft,
    GaussianDft,
}

impl QcRunBackend {
    pub fn label(self) -> &'static str {
        match self {
            QcRunBackend::XtbGfn2 => "GFN2-xTB (--ohess)",
            QcRunBackend::XtbGfnFf => "GFN-FF (--ohess, ultra-fast)",
            QcRunBackend::OrcaDft => "ORCA (Opt Freq)",
            QcRunBackend::GaussianDft => "Gaussian (Opt Freq)",
        }
    }
}

#[derive(Debug, Clone)]
pub struct QcRunOutcome {
    pub cache_key: String,
    pub backend: QcRunBackend,
    pub report: QuantumOutputReport,
    pub ir_prediction: IrPrediction,
    pub elapsed_secs: f64,
}

pub struct QuantumRunnerState {
    pub binaries: QcBinaryAvailability,
    pub running: bool,
    pub running_label: String,
    pub started_at: Option<Instant>,
    pub last_error: String,
    pub last_outcome: Option<QcRunOutcome>,
    pub cache: HashMap<String, QcRunOutcome>,
    rx: Option<Receiver<Result<QcRunOutcome, String>>>,
    cancel_flag: Arc<AtomicBool>,
}

impl Default for QuantumRunnerState {
    fn default() -> Self {
        Self {
            binaries: QcBinaryAvailability::default(),
            running: false,
            running_label: String::new(),
            started_at: None,
            last_error: String::new(),
            last_outcome: None,
            cache: HashMap::new(),
            rx: None,
            cancel_flag: Arc::new(AtomicBool::new(false)),
        }
    }
}

impl QuantumRunnerState {
    pub fn ensure_detected(&mut self, xtb_cfg: &str, orca_cfg: &str, gauss_cfg: &str) {
        if !self.binaries.detected {
            self.binaries = detect_qc_binaries(xtb_cfg, orca_cfg, gauss_cfg);
        }
    }

    pub fn redetect(&mut self, xtb_cfg: &str, orca_cfg: &str, gauss_cfg: &str) {
        self.binaries = detect_qc_binaries(xtb_cfg, orca_cfg, gauss_cfg);
    }

    pub fn elapsed_secs(&self) -> f64 {
        self.started_at
            .map(|t| t.elapsed().as_secs_f64())
            .unwrap_or(0.0)
    }

    pub fn cancel(&mut self) {
        self.cancel_flag.store(true, Ordering::SeqCst);
        self.running = false;
        self.rx = None;
        self.last_error = "Quantum FTIR calculation cancelled.".to_string();
    }

    /// Launches a local background `xtb` / `orca` / `g16` harmonic IR job (or serves instantly from cache).
    pub fn start_ir_job(
        &mut self,
        backend: QcRunBackend,
        mol: &ChemMolecule,
        cache_key: String,
        mut spec: QuantumJobSpec,
        bin_path: String,
    ) {
        if let Some(hit) = self.cache.get(&cache_key).cloned() {
            self.last_error.clear();
            self.last_outcome = Some(hit);
            return;
        }
        self.cancel_flag.store(true, Ordering::SeqCst);
        let cancel = Arc::new(AtomicBool::new(false));
        self.cancel_flag = Arc::clone(&cancel);
        self.running = true;
        self.running_label = backend.label().to_string();
        self.started_at = Some(Instant::now());
        self.last_error.clear();

        spec.job_type = QuantumJobType::OptFreq;
        spec.add_explicit_hydrogens = true;
        spec.pre_optimize_uff = true;
        match backend {
            QcRunBackend::XtbGfn2 => {
                spec.program = QuantumProgram::Xtb;
                spec.xtb_gfn_level = 2;
            }
            QcRunBackend::XtbGfnFf => {
                spec.program = QuantumProgram::Xtb;
                spec.xtb_gfn_level = 0;
            }
            QcRunBackend::OrcaDft => {
                spec.program = QuantumProgram::Orca;
            }
            QcRunBackend::GaussianDft => {
                spec.program = QuantumProgram::Gaussian;
            }
        }

        let mol_clone = mol.clone();
        let (tx, rx) = mpsc::channel();
        self.rx = Some(rx);

        thread::spawn(move || {
            let res = run_qc_job_blocking(backend, &mol_clone, cache_key, &spec, &bin_path, &cancel);
            let _ = tx.send(res);
        });
    }

    /// Non-blocking per-frame poll; returns `Some(Ok(outcome))` or `Some(Err(msg))` when a job completes.
    pub fn poll(&mut self) -> Option<Result<QcRunOutcome, String>> {
        let rx = self.rx.as_ref()?;
        match rx.try_recv() {
            Ok(res) => {
                self.running = false;
                self.rx = None;
                match &res {
                    Ok(out) => {
                        self.last_error.clear();
                        self.cache.insert(out.cache_key.clone(), out.clone());
                        self.last_outcome = Some(out.clone());
                    }
                    Err(e) => {
                        self.last_error = e.clone();
                    }
                }
                Some(res)
            }
            Err(mpsc::TryRecvError::Empty) => None,
            Err(mpsc::TryRecvError::Disconnected) => {
                self.running = false;
                self.rx = None;
                let err = "Quantum worker thread terminated unexpectedly".to_string();
                self.last_error = err.clone();
                Some(Err(err))
            }
        }
    }
}

fn wait_child_cancellable(
    mut child: std::process::Child,
    cancel: &AtomicBool,
    timeout: Duration,
) -> Result<std::process::Output, String> {
    let start = Instant::now();
    loop {
        if cancel.load(Ordering::SeqCst) {
            let _ = child.kill();
            let _ = child.wait();
            return Err("Cancelled by user".to_string());
        }
        if start.elapsed() > timeout {
            let _ = child.kill();
            let _ = child.wait();
            return Err(format!("Quantum job timed out after {} s", timeout.as_secs()));
        }
        match child.try_wait() {
            Ok(Some(_)) => {
                return child
                    .wait_with_output()
                    .map_err(|e| format!("Failed reading process output: {e}"));
            }
            Ok(None) => thread::sleep(Duration::from_millis(80)),
            Err(e) => return Err(format!("Process wait error: {e}")),
        }
    }
}

/// Runs a local quantum IR job in an isolated temporary directory and parses the resulting harmonic spectrum.
pub fn run_qc_job_blocking(
    backend: QcRunBackend,
    mol: &ChemMolecule,
    cache_key: String,
    spec: &QuantumJobSpec,
    bin_path: &str,
    cancel: &AtomicBool,
) -> Result<QcRunOutcome, String> {
    let t0 = Instant::now();
    let prep = prepare_molecule_for_quantum(mol, true, true);
    let charge = if spec.charge == 0 {
        total_formal_charge(&prep.molecule)
    } else {
        spec.charge
    };
    let uhf = spec.multiplicity.saturating_sub(1);

    let seq = JOB_COUNTER.fetch_add(1, Ordering::Relaxed);
    let workdir = std::env::temp_dir().join(format!("chiralcraft-qc-{}-{seq}", std::process::id()));
    fs::create_dir_all(&workdir).map_err(|e| format!("Failed creating scratch dir: {e}"))?;

    // Guard cleanup on drop
    struct CleanupDir(PathBuf);
    impl Drop for CleanupDir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    let _cleanup = CleanupDir(workdir.clone());

    let mut clean_spec = spec.clone();
    // Already materialized H + UFF optimized in `prep`
    clean_spec.add_explicit_hydrogens = false;
    clean_spec.pre_optimize_uff = false;
    clean_spec.charge = charge;

    let report = match backend {
        QcRunBackend::XtbGfn2 | QcRunBackend::XtbGfnFf => {
            let xyz = export_xtb_input(&prep.molecule, &clean_spec);
            fs::write(workdir.join("input.xyz"), xyz).map_err(|e| format!("Write input.xyz: {e}"))?;
            let mut cmd = Command::new(bin_path);
            cmd.current_dir(&workdir)
                .arg("input.xyz")
                .arg(if backend == QcRunBackend::XtbGfnFf { "--gfnff" } else { "--gfn" });
            if backend == QcRunBackend::XtbGfn2 {
                cmd.arg("2");
            }
            cmd.args(["--ohess", "--chrg", &charge.to_string(), "--uhf", &uhf.to_string()]);
            if !spec.solvent.eq_ignore_ascii_case("none") && !spec.solvent.trim().is_empty() {
                cmd.args(["--alpb", &spec.solvent.to_ascii_lowercase()]);
            }
            if spec.nprocs > 1 {
                cmd.args(["--parallel", &spec.nprocs.to_string()]);
            }
            cmd.stdin(Stdio::null()).stdout(Stdio::piped()).stderr(Stdio::piped());
            let child = cmd.spawn().map_err(|e| format!("Failed to spawn {bin_path}: {e}"))?;
            let out = wait_child_cancellable(child, cancel, Duration::from_secs(600))?;
            let stdout_s = String::from_utf8_lossy(&out.stdout).to_string();
            if !out.status.success() {
                let stderr_s = String::from_utf8_lossy(&out.stderr);
                let tail: String = stdout_s.lines().rev().take(8).collect::<Vec<_>>().into_iter().rev().collect::<Vec<_>>().join("\n");
                return Err(format!("xTB exited with {}: {}\n{}", out.status, stderr_s.trim(), tail));
            }
            parse_xtb_workdir(&workdir, &stdout_s, backend == QcRunBackend::XtbGfnFf)?
        }
        QcRunBackend::OrcaDft => {
            let inp = export_orca_input(&prep.molecule, &clean_spec);
            let inp_path = workdir.join("job.inp");
            fs::write(&inp_path, inp).map_err(|e| format!("Write job.inp: {e}"))?;
            let child = Command::new(bin_path)
                .current_dir(&workdir)
                .arg("job.inp")
                .stdin(Stdio::null())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to spawn ORCA ({bin_path}): {e}"))?;
            let out = wait_child_cancellable(child, cancel, Duration::from_secs(1800))?;
            let stdout_s = String::from_utf8_lossy(&out.stdout).to_string();
            if !out.status.success() {
                return Err(format!("ORCA exited with {}", out.status));
            }
            parse_quantum_output(&stdout_s)?
        }
        QcRunBackend::GaussianDft => {
            let gjf = export_gaussian_input(&prep.molecule, &clean_spec);
            fs::write(workdir.join("job.gjf"), gjf).map_err(|e| format!("Write job.gjf: {e}"))?;
            let child = Command::new(bin_path)
                .current_dir(&workdir)
                .arg("job.gjf")
                .stdin(Stdio::null())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .map_err(|e| format!("Failed to spawn Gaussian ({bin_path}): {e}"))?;
            let out = wait_child_cancellable(child, cancel, Duration::from_secs(1800))?;
            let log_s = fs::read_to_string(workdir.join("job.log"))
                .unwrap_or_else(|_| String::from_utf8_lossy(&out.stdout).to_string());
            if !out.status.success() && !log_s.contains("Normal termination") {
                return Err(format!("Gaussian exited with {}", out.status));
            }
            parse_quantum_output(&log_s)?
        }
    };

    let ir_prediction = report
        .to_ir_prediction(Some(&prep.parent_heavy_atom))
        .ok_or_else(|| "Quantum run finished but no harmonic IR modes were found in output".to_string())?;

    Ok(QcRunOutcome {
        cache_key,
        backend,
        report,
        ir_prediction,
        elapsed_secs: t0.elapsed().as_secs_f64(),
    })
}

// ============================================================================
// 3. Optional Online NIST Chemistry WebBook IR (.jdx) Fetcher
// ============================================================================

pub struct NistIrFetchState {
    pub fetching: bool,
    pub status: String,
    rx: Option<Receiver<Result<JcampSpectrum, String>>>,
}

impl Default for NistIrFetchState {
    fn default() -> Self {
        Self {
            fetching: false,
            status: String::new(),
            rx: None,
        }
    }
}

fn url_encode_query(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' | b'~' => out.push(b as char),
            b' ' => out.push('+'),
            _ => out.push_str(&format!("%{:02X}", b)),
        }
    }
    out
}

fn curl_get_text(url: &str) -> Result<String, String> {
    let out = Command::new("curl")
        .args(["-fsSL", "--max-time", "15", "-A", "ChiralCraft/1.0", url])
        .stdin(Stdio::null())
        .output()
        .map_err(|e| format!("Failed to invoke curl: {e}"))?;
    if !out.status.success() {
        return Err(format!("HTTP fetch failed ({})", out.status));
    }
    Ok(String::from_utf8_lossy(&out.stdout).to_string())
}

/// Extracts `(nist_id, index)` for an IR JCAMP link from a NIST Chemistry WebBook HTML page.
pub fn extract_nist_ir_jcamp_params(html: &str) -> Option<(String, String)> {
    // Look for direct JCAMP download link: JCAMP=C67641&amp;Index=0&amp;Type=IR or JCAMP=C67641&Index=0&Type=IR
    for marker in ["JCAMP=", "jcamp="] {
        let mut rest = html;
        while let Some(pos) = rest.find(marker) {
            let after = &rest[pos + marker.len()..];
            let end = after
                .find(|c: char| c == '"' || c == '\'' || c == '<' || c == '>' || c.is_whitespace())
                .unwrap_or(after.len());
            let frag = after[..end].replace("&amp;", "&");
            let id = frag.split('&').next().unwrap_or("").trim().to_string();
            let is_ir = frag.to_ascii_uppercase().contains("TYPE=IR") || !frag.to_ascii_uppercase().contains("TYPE=");
            if id.starts_with('C') && id.len() >= 4 && is_ir {
                let idx = frag
                    .split('&')
                    .find_map(|kv| {
                        let (k, v) = kv.split_once('=')?;
                        if k.eq_ignore_ascii_case("Index") {
                            Some(v.to_string())
                        } else {
                            None
                        }
                    })
                    .unwrap_or_else(|| "0".to_string());
                return Some((id, idx));
            }
            rest = &after[end..];
        }
    }
    // Fallback: find compound ID=C... when Mask=80 (IR spectrum page) is present
    if let Some(pos) = html.find("ID=C") {
        let after = &html[pos + 3..];
        let id: String = after.chars().take_while(|c| c.is_ascii_alphanumeric()).collect();
        if id.len() >= 4 {
            return Some((id, "0".to_string()));
        }
    }
    None
}

/// Queries NIST Chemistry WebBook by CAS, InChIKey, InChI, or compound name and downloads its IR `.jdx`.
pub fn fetch_nist_ir_blocking(query: &str, inchi_fallback: Option<&str>) -> Result<JcampSpectrum, String> {
    let q = query.trim();
    if q.is_empty() && inchi_fallback.is_none() {
        return Err("Provide a compound name, CAS number, or active molecule to query NIST WebBook".to_string());
    }

    // Build candidate NIST URLs
    let mut urls: Vec<String> = Vec::new();
    let looks_cas = !q.is_empty() && q.chars().all(|c| c.is_ascii_digit() || c == '-') && q.contains('-');
    if looks_cas {
        let cid = format!("C{}", q.replace('-', ""));
        urls.push(format!("https://webbook.nist.gov/cgi/cbook.cgi?JCAMP={cid}&Index=0&Type=IR"));
        urls.push(format!("https://webbook.nist.gov/cgi/cbook.cgi?ID={cid}&Units=SI&Mask=80"));
    } else if !q.is_empty() {
        // Handle ChiralCraft's "IUPAC name (common name)" format by trying both cleanly!
        if let Some((before, rest)) = q.split_once('(') {
            let common = rest.trim_end_matches(')').trim();
            if !common.is_empty() {
                urls.push(format!(
                    "https://webbook.nist.gov/cgi/cbook.cgi?Name={}&Units=SI&cIR=on",
                    url_encode_query(common)
                ));
            }
            let iupac = before.trim();
            if !iupac.is_empty() {
                urls.push(format!(
                    "https://webbook.nist.gov/cgi/cbook.cgi?Name={}&Units=SI&cIR=on",
                    url_encode_query(iupac)
                ));
            }
        } else {
            urls.push(format!(
                "https://webbook.nist.gov/cgi/cbook.cgi?Name={}&Units=SI&cIR=on",
                url_encode_query(q)
            ));
        }
    }
    if let Some(inchi) = inchi_fallback.filter(|s| !s.trim().is_empty()) {
        urls.push(format!(
            "https://webbook.nist.gov/cgi/cbook.cgi?InChI={}&Units=SI&cIR=on",
            url_encode_query(inchi)
        ));
    }

    let mut last_err = String::from("No IR spectrum found on NIST Chemistry WebBook");
    for url in urls {
        let Ok(body) = curl_get_text(&url) else {
            last_err = "Could not reach webbook.nist.gov (offline or rate-limited)".to_string();
            continue;
        };
        if body.trim_start().starts_with("##TITLE=") && body.contains("##XYDATA=") {
            return parse_jcamp_dx(&body);
        }
        if let Some((cid, idx)) = extract_nist_ir_jcamp_params(&body) {
            let jdx_url = format!("https://webbook.nist.gov/cgi/cbook.cgi?JCAMP={cid}&Index={idx}&Type=IR");
            if let Ok(jdx_txt) = curl_get_text(&jdx_url) {
                if jdx_txt.contains("##XYDATA=") || jdx_txt.contains("##PEAK TABLE=") {
                    return parse_jcamp_dx(&jdx_txt);
                }
            }
        }
    }
    Err(last_err)
}

impl NistIrFetchState {
    pub fn start_fetch(&mut self, query: String, inchi_fallback: Option<String>) {
        self.fetching = true;
        self.status = format!("Querying NIST WebBook for '{}'...", query.trim());
        let (tx, rx) = mpsc::channel();
        self.rx = Some(rx);
        thread::spawn(move || {
            let res = fetch_nist_ir_blocking(&query, inchi_fallback.as_deref());
            let _ = tx.send(res);
        });
    }

    pub fn poll(&mut self) -> Option<Result<JcampSpectrum, String>> {
        let rx = self.rx.as_ref()?;
        match rx.try_recv() {
            Ok(res) => {
                self.fetching = false;
                self.rx = None;
                match &res {
                    Ok(sp) => self.status = format!("Loaded NIST IR: {}", sp.title),
                    Err(e) => self.status = format!("NIST lookup failed: {e}"),
                }
                Some(res)
            }
            Err(mpsc::TryRecvError::Empty) => None,
            Err(mpsc::TryRecvError::Disconnected) => {
                self.fetching = false;
                self.rx = None;
                let err = "NIST fetch thread closed unexpectedly".to_string();
                self.status = err.clone();
                Some(Err(err))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chiralcraft_formats::smiles_to_molecule;
    #[cfg(unix)]
    use std::os::unix::fs::PermissionsExt;

    #[test]
    fn test_nist_html_jcamp_param_extraction() {
        let html = r#"<a href="/cgi/cbook.cgi?JCAMP=C67641&amp;Index=0&amp;Type=IR">Download JCAMP-DX</a>"#;
        let (id, idx) = extract_nist_ir_jcamp_params(html).expect("extract nist id");
        assert_eq!(id, "C67641");
        assert_eq!(idx, "0");
    }

    #[cfg(unix)]
    #[test]
    fn test_binary_detection_rejects_gnome_orca_and_runs_mock_xtb() {
        let tmp = std::env::temp_dir().join(format!("chiralcraft-test-qc-{}", std::process::id()));
        let _ = fs::create_dir_all(&tmp);

        // 1. Mock GNOME /usr/bin/orca Python script must be REJECTED
        let gnome_orca = tmp.join("gnome_orca");
        fs::write(&gnome_orca, "#!/usr/bin/python3\nprint('orca 46.1 GNOME screen reader')\n").unwrap();
        fs::set_permissions(&gnome_orca, fs::Permissions::from_mode(0o755)).unwrap();
        assert!(!is_valid_quantum_orca_binary(&gnome_orca), "GNOME screen reader must be rejected");

        // 2. Mock quantum xTB script writing g98.out + xtbopt.xyz verifies full async runner end-to-end!
        let mock_xtb = tmp.join("mock_xtb.sh");
        let script = r#"#!/bin/sh
if [ "$1" = "--version" ]; then
  echo "xtb version 6.7.0 (mock)"
  exit 0
fi
cat << 'EOF' > g98.out
 Gaussian 98: xtb output
 Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic     Atomic              Coordinates (Angstroms)
 Number     Number      Type              X           Y           Z
 ---------------------------------------------------------------------
    1          8           0        0.000000    0.000000    0.117300
    2          1           0        0.000000    0.757200   -0.469200
    3          1           0        0.000000   -0.757200   -0.469200
 ---------------------------------------------------------------------
                      1                      2
                     A1                     B2
 Frequencies --   1595.10                3755.80
 Red. masses --      1.08                   1.04
 IR Inten    --     68.50                  42.10
 Atom AN      X      Y      Z        X      Y      Z
   1   8     0.00   0.00   0.07     0.00   0.07   0.00
   2   1     0.00   0.43  -0.55     0.00  -0.55   0.43
   3   1     0.00  -0.43  -0.55     0.00  -0.55  -0.43
 Normal termination of Gaussian 98.
EOF
echo "         | TOTAL ENERGY              -5.070543 Eh |"
echo "         | HOMO-LUMO GAP             12.400000 eV |"
echo "normal termination of xtb"
"#;
        fs::write(&mock_xtb, script).unwrap();
        fs::set_permissions(&mock_xtb, fs::Permissions::from_mode(0o755)).unwrap();

        let water = smiles_to_molecule("O").unwrap();
        let cancel = AtomicBool::new(false);
        let outcome = run_qc_job_blocking(
            QcRunBackend::XtbGfn2,
            &water,
            "O|gfn2".to_string(),
            &QuantumJobSpec::default(),
            &mock_xtb.to_string_lossy(),
            &cancel,
        )
        .expect("mock xTB job must succeed");

        assert_eq!(outcome.report.program, "xTB");
        assert_eq!(outcome.ir_prediction.bands.len(), 2);
        assert!(outcome.ir_prediction.mode_geometry_3d.is_some());
        assert!(outcome.ir_prediction.bands.iter().all(|b| b.normal_mode_vectors.is_some()));

        let _ = fs::remove_dir_all(&tmp);
    }
}
