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
        #[cfg(windows)]
        {
            if !bin_name.to_ascii_lowercase().ends_with(".exe") {
                let cand_exe = dir.join(format!("{bin_name}.exe"));
                if cand_exe.is_file() {
                    return Some(cand_exe);
                }
            }
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
    // Quick 256-byte header inspection (without loading a huge binary into memory)
    if let Ok(mut f) = fs::File::open(path) {
        use std::io::Read;
        let mut head = [0u8; 256];
        if let Ok(n) = f.read(&mut head) {
            if let Ok(txt) = std::str::from_utf8(&head[..n]) {
                let low = txt.to_ascii_lowercase();
                if low.starts_with("#!")
                    && (low.contains("python") || low.contains("gnome") || low.contains("orca.main"))
                {
                    return false;
                }
            }
        }
    }
    // Confirm it doesn't identify as GNOME screen reader; in tests, mock scripts may echo ORCA banner
    let Some(out) = output_with_timeout(Command::new(path).arg("--version"), VERSION_PROBE_TIMEOUT) else {
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
        if let Some(out) = output_with_timeout(Command::new(&p).arg("--version"), VERSION_PROBE_TIMEOUT) {
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
            // Supersede any in-flight job, else its later result would overwrite this cached one.
            self.cancel_flag.store(true, Ordering::SeqCst);
            self.rx = None;
            self.running = false;
            self.started_at = None;
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

/// Spawns `cmd` as leader of its own process group (unix) so cancel/timeout can also reach
/// grandchildren such as ORCA's MPI workers.
fn spawn_in_own_group(cmd: &mut Command) -> std::io::Result<std::process::Child> {
    #[cfg(unix)]
    {
        use std::os::unix::process::CommandExt;
        cmd.process_group(0);
    }
    cmd.spawn()
}

/// Kills a child spawned by `spawn_in_own_group` and (unix) its whole process group:
/// SIGTERM, short grace, then SIGKILL; `child.kill()` remains the fallback (and the Windows path).
fn kill_process_tree(child: &mut std::process::Child) {
    #[cfg(unix)]
    {
        let pgid = format!("-{}", child.id());
        let signal_group = |sig: &str| {
            let _ = Command::new("kill")
                .args(["-s", sig, "--", &pgid])
                .stdin(Stdio::null())
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status();
        };
        signal_group("TERM");
        let grace = Instant::now();
        while grace.elapsed() < Duration::from_millis(500) {
            if matches!(child.try_wait(), Ok(Some(_))) {
                break;
            }
            thread::sleep(Duration::from_millis(25));
        }
        // Even if the leader exited on TERM, grandchildren may linger in the group.
        signal_group("KILL");
    }
    let _ = child.kill();
    let _ = child.wait();
}

/// Like `Command::output()` but bounded: the process (group) is killed after `timeout`.
fn output_with_timeout(cmd: &mut Command, timeout: Duration) -> Option<std::process::Output> {
    use std::io::Read;
    let mut child = spawn_in_own_group(cmd.stdin(Stdio::null()).stdout(Stdio::piped()).stderr(Stdio::piped())).ok()?;
    // Drain pipes on threads so a chatty process cannot block on a full pipe.
    let mut so = child.stdout.take()?;
    let mut se = child.stderr.take()?;
    let t_out = thread::spawn(move || {
        let mut v = Vec::new();
        let _ = so.read_to_end(&mut v);
        v
    });
    let t_err = thread::spawn(move || {
        let mut v = Vec::new();
        let _ = se.read_to_end(&mut v);
        v
    });
    let start = Instant::now();
    let status = loop {
        match child.try_wait() {
            Ok(Some(st)) => break st,
            Ok(None) if start.elapsed() < timeout => thread::sleep(Duration::from_millis(20)),
            _ => {
                kill_process_tree(&mut child);
                return None;
            }
        }
    };
    // A lingering grandchild could hold the pipes open; reap the group if they don't close in time.
    while !(t_out.is_finished() && t_err.is_finished()) && start.elapsed() < timeout {
        thread::sleep(Duration::from_millis(20));
    }
    if !(t_out.is_finished() && t_err.is_finished()) {
        kill_process_tree(&mut child);
    }
    Some(std::process::Output {
        status,
        stdout: t_out.join().unwrap_or_default(),
        stderr: t_err.join().unwrap_or_default(),
    })
}

const VERSION_PROBE_TIMEOUT: Duration = Duration::from_secs(3);

fn wait_child_status_cancellable(
    mut child: std::process::Child,
    cancel: &AtomicBool,
    timeout: Duration,
) -> Result<std::process::ExitStatus, String> {
    let start = Instant::now();
    loop {
        if cancel.load(Ordering::SeqCst) {
            kill_process_tree(&mut child);
            return Err("Cancelled by user".to_string());
        }
        if start.elapsed() > timeout {
            kill_process_tree(&mut child);
            return Err(format!("Quantum job timed out after {} s", timeout.as_secs()));
        }
        match child.try_wait() {
            Ok(Some(status)) => return Ok(status),
            Ok(None) => thread::sleep(Duration::from_millis(50)),
            Err(e) => return Err(format!("Process wait error: {e}")),
        }
    }
}

fn create_private_qc_workdir() -> Result<PathBuf, String> {
    let base = std::env::temp_dir();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    for attempt in 0..32u32 {
        let seq = JOB_COUNTER.fetch_add(1, Ordering::Relaxed);
        let cand = base.join(format!(
            "chiralcraft-qc-{}-{nanos}-{seq}-{attempt}",
            std::process::id()
        ));
        // Create with 0o700 atomically (no window where the scratch dir is world-accessible).
        let mut builder = fs::DirBuilder::new();
        #[cfg(unix)]
        {
            use std::os::unix::fs::DirBuilderExt;
            builder.mode(0o700);
        }
        match builder.create(&cand) {
            Ok(()) => return Ok(cand),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(format!("Failed creating scratch dir: {e}")),
        }
    }
    Err("Failed creating unique quantum scratch dir".to_string())
}

/// The child runs with `current_dir(workdir)`, so a relative path like `bin/xtb` would be resolved
/// against the scratch dir. Anchor such paths to our cwd; bare names (`xtb`) still use PATH lookup.
fn resolve_bin_for_workdir(bin_path: &str) -> PathBuf {
    let p = PathBuf::from(bin_path.trim());
    if p.is_relative() && p.components().count() > 1 {
        if let Ok(cwd) = std::env::current_dir() {
            return cwd.join(p);
        }
    }
    p
}

fn normalize_xtb_solvent(raw: &str) -> Result<Option<String>, String> {
    let t = raw.trim();
    if t.is_empty() || t.eq_ignore_ascii_case("none") {
        return Ok(None);
    }
    let low = t.to_ascii_lowercase();
    let mapped = match low.as_str() {
        "dcm" => "ch2cl2",
        "mecn" => "acetonitrile",
        other => other,
    };
    if mapped.starts_with('-')
        || !mapped
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
    {
        return Err(format!("Invalid xTB solvent name: '{raw}'"));
    }
    Ok(Some(mapped.to_string()))
}

fn run_command_logged(
    mut cmd: Command,
    workdir: &Path,
    cancel: &AtomicBool,
    timeout: Duration,
    label: &str,
) -> Result<(std::process::ExitStatus, String, String), String> {
    let out_path = workdir.join("stdout.log");
    let err_path = workdir.join("stderr.log");
    let out_file = fs::File::create(&out_path).map_err(|e| format!("Create stdout.log: {e}"))?;
    let err_file = fs::File::create(&err_path).map_err(|e| format!("Create stderr.log: {e}"))?;
    let child = spawn_in_own_group(
        cmd.stdin(Stdio::null())
            .stdout(Stdio::from(out_file))
            .stderr(Stdio::from(err_file)),
    )
    .map_err(|e| format!("Failed to spawn {label}: {e}"))?;
    let status = wait_child_status_cancellable(child, cancel, timeout)?;
    let stdout_s = fs::read_to_string(&out_path).unwrap_or_default();
    let stderr_s = fs::read_to_string(&err_path).unwrap_or_default();
    Ok((status, stdout_s, stderr_s))
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

    let workdir = create_private_qc_workdir()?;
    let resolved_bin = resolve_bin_for_workdir(bin_path);
    let resolved_bin = resolved_bin.to_string_lossy(); // input was &str, so lossless
    let bin_path: &str = &resolved_bin;

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
            let xtb_solvent = normalize_xtb_solvent(&spec.solvent)?;
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
            if let Some(solv) = xtb_solvent {
                cmd.args(["--alpb", &solv]);
            }
            if spec.nprocs > 1 {
                cmd.args(["--parallel", &spec.nprocs.to_string()]);
            }
            let (status, stdout_s, stderr_s) =
                run_command_logged(cmd, &workdir, cancel, Duration::from_secs(600), bin_path)?;
            if !status.success() {
                let tail: String = stdout_s.lines().rev().take(8).collect::<Vec<_>>().into_iter().rev().collect::<Vec<_>>().join("\n");
                return Err(format!("xTB exited with {status}: {}\n{}", stderr_s.trim(), tail));
            }
            parse_xtb_workdir(&workdir, &stdout_s, backend == QcRunBackend::XtbGfnFf)?
        }
        QcRunBackend::OrcaDft => {
            let inp = export_orca_input(&prep.molecule, &clean_spec);
            let inp_path = workdir.join("job.inp");
            fs::write(&inp_path, inp).map_err(|e| format!("Write job.inp: {e}"))?;
            let mut cmd = Command::new(bin_path);
            cmd.current_dir(&workdir).arg("job.inp");
            let (status, stdout_s, _) = run_command_logged(
                cmd,
                &workdir,
                cancel,
                Duration::from_secs(1800),
                &format!("ORCA ({bin_path})"),
            )?;
            if !status.success() {
                return Err(format!("ORCA exited with {status}"));
            }
            parse_quantum_output(&stdout_s)?
        }
        QcRunBackend::GaussianDft => {
            let gjf = export_gaussian_input(&prep.molecule, &clean_spec);
            fs::write(workdir.join("job.gjf"), gjf).map_err(|e| format!("Write job.gjf: {e}"))?;
            let mut cmd = Command::new(bin_path);
            cmd.current_dir(&workdir).arg("job.gjf");
            let (status, stdout_s, _) = run_command_logged(
                cmd,
                &workdir,
                cancel,
                Duration::from_secs(1800),
                &format!("Gaussian ({bin_path})"),
            )?;
            let log_s = fs::read_to_string(workdir.join("job.log")).unwrap_or(stdout_s);
            if !status.success() && !log_s.contains("Normal termination") {
                return Err(format!("Gaussian exited with {status}"));
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
        .args([
            "-fsSL",
            "--globoff",
            "--proto",
            "=https",
            "--proto-redir",
            "=https",
            "--max-time",
            "15",
            "-A",
            "ChiralCraft/1.0",
            "--",
            url,
        ])
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
    let valid_cid = |s: &str| {
        s.starts_with('C') && s.len() >= 4 && s[1..].bytes().all(|b| b.is_ascii_digit())
    };
    let valid_idx = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
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
            if valid_cid(&id) && is_ir {
                let idx = frag
                    .split('&')
                    .find_map(|kv| {
                        let (k, v) = kv.split_once('=')?;
                        if k.eq_ignore_ascii_case("Index") && valid_idx(v.trim()) {
                            Some(v.trim().to_string())
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
        if valid_cid(&id) {
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
        // Handle ChiralCraft's "IUPAC name (common name)" suffix cleanly without mangling "(R)-..." or "2-(acetoxy)..."
        if q.ends_with(')') {
            if let Some((before, rest)) = q.rsplit_once(" (") {
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
            }
        }
        urls.push(format!(
            "https://webbook.nist.gov/cgi/cbook.cgi?Name={}&Units=SI&cIR=on",
            url_encode_query(q)
        ));
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

    #[test]
    fn test_relative_bin_anchored_and_workdir_private() {
        let cwd = std::env::current_dir().unwrap();
        assert_eq!(resolve_bin_for_workdir("bin/xtb"), cwd.join("bin/xtb"));
        assert_eq!(resolve_bin_for_workdir("xtb"), PathBuf::from("xtb"));
        assert_eq!(resolve_bin_for_workdir("/opt/xtb"), PathBuf::from("/opt/xtb"));
        let wd = create_private_qc_workdir().unwrap();
        #[cfg(unix)]
        assert_eq!(fs::metadata(&wd).unwrap().permissions().mode() & 0o777, 0o700);
        let _ = fs::remove_dir(&wd);
    }

    #[cfg(unix)]
    #[test]
    fn test_version_probe_bounded_and_cancel_kills_grandchildren() {
        let tmp = std::env::temp_dir().join(format!("chiralcraft-test-pg-{}", std::process::id()));
        let _ = fs::create_dir_all(&tmp);
        // Hanging `--version` must not block detection forever.
        let hang = tmp.join("hang.sh");
        fs::write(&hang, "#!/bin/sh\nsleep 30\n").unwrap();
        fs::set_permissions(&hang, fs::Permissions::from_mode(0o755)).unwrap();
        let t0 = Instant::now();
        assert!(!is_valid_quantum_orca_binary(&hang));
        assert!(t0.elapsed() < Duration::from_secs(10));

        // Cancel must also kill a backgrounded grandchild (MPI-worker stand-in).
        let pidfile = tmp.join("gc.pid");
        let mut cmd = Command::new("sh");
        cmd.arg("-c").arg(format!("sleep 30 & echo $! > '{}'; wait", pidfile.display()));
        let child = spawn_in_own_group(&mut cmd).unwrap();
        while !fs::read_to_string(&pidfile).map(|s| s.ends_with('\n')).unwrap_or(false) {
            thread::sleep(Duration::from_millis(20));
        }
        let gc_pid = fs::read_to_string(&pidfile).unwrap().trim().to_string();
        let cancel = AtomicBool::new(true);
        assert!(wait_child_status_cancellable(child, &cancel, Duration::from_secs(60)).is_err());
        let alive = || Command::new("kill").args(["-0", &gc_pid]).stderr(Stdio::null()).status().map(|s| s.success()).unwrap_or(false);
        let t1 = Instant::now();
        while alive() && t1.elapsed() < Duration::from_secs(3) {
            thread::sleep(Duration::from_millis(50));
        }
        assert!(!alive(), "grandchild {gc_pid} survived cancel");
        let _ = fs::remove_dir_all(&tmp);
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

        // 3. A cache hit while another job is in flight must supersede (cancel) that job.
        let mut st = QuantumRunnerState::default();
        st.cache.insert("O|gfn2".to_string(), outcome);
        let slow = tmp.join("slow.sh");
        fs::write(&slow, "#!/bin/sh\nsleep 30\n").unwrap();
        fs::set_permissions(&slow, fs::Permissions::from_mode(0o755)).unwrap();
        st.start_ir_job(QcRunBackend::XtbGfn2, &water, "other".into(), QuantumJobSpec::default(), slow.to_string_lossy().into());
        assert!(st.running);
        st.start_ir_job(QcRunBackend::XtbGfn2, &water, "O|gfn2".into(), QuantumJobSpec::default(), slow.to_string_lossy().into());
        assert!(!st.running && st.poll().is_none());
        assert_eq!(st.last_outcome.as_ref().map(|o| o.cache_key.as_str()), Some("O|gfn2"));

        let _ = fs::remove_dir_all(&tmp);
    }
}
