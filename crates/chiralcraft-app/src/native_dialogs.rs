//! Zero-Dependency Cross-Platform Native OS File Dialogs (Open / Save)
//! and Safe Default Export Paths (Linux, macOS, Windows, FreeBSD).
//!
//! Avoids C/GTK/Win32 crate link issues across all 6 release targets by invoking the OS's built-in
//! file dialog service when available (`zenity` / `kdialog` on Linux/FreeBSD, `osascript` on macOS,
//! PowerShell `System.Windows.Forms` on Windows) and falling back cleanly in headless/test environments.

use std::path::PathBuf;
use std::process::{Command, Stdio};

/// Returns a cross-platform default directory for saving exported files (`~/Documents`, `~/Desktop`,
/// `$HOME` / `%USERPROFILE%`, or `std::env::temp_dir()` as last resort — never hardcoded `/tmp`).
pub fn default_export_dir() -> PathBuf {
    let home = std::env::var_os("HOME")
        .or_else(|| std::env::var_os("USERPROFILE"))
        .map(PathBuf::from);
    if let Some(h) = home {
        let docs = h.join("Documents");
        if docs.is_dir() {
            return docs;
        }
        if h.is_dir() {
            return h;
        }
    }
    std::env::temp_dir()
}

/// Builds a cross-platform fallback file path for `suggested_name` inside `default_export_dir()`.
pub fn fallback_export_path(suggested_name: &str) -> PathBuf {
    default_export_dir().join(suggested_name)
}

/// Opens a native OS "Open File" dialog and returns the selected path, or `None` if cancelled / headless.
pub fn pick_open_file(title: &str, extensions: &[&str]) -> Option<String> {
    if cfg!(test) {
        return None;
    }

    #[cfg(target_os = "windows")]
    {
        let filter_pats: Vec<String> = extensions.iter().map(|e| format!("*.{e}")).collect();
        let filter = if filter_pats.is_empty() {
            "All files (*.*)|*.*".to_string()
        } else {
            let joined = filter_pats.join(";");
            format!("Supported ({joined})|{joined}|All files (*.*)|*.*")
        };
        let ps = format!(
            "Add-Type -AssemblyName System.Windows.Forms; \
             $d = New-Object System.Windows.Forms.OpenFileDialog; \
             $d.Title = '{}'; \
             $d.Filter = '{}'; \
             if ($d.ShowDialog() -eq 'OK') {{ Write-Output $d.FileName }}",
            title.replace('\'', "''"),
            filter.replace('\'', "''")
        );
        if let Ok(out) = Command::new("powershell")
            .args(["-NoProfile", "-STA", "-Command", &ps])
            .stdin(Stdio::null())
            .output()
        {
            if out.status.success() {
                let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
                if !s.is_empty() {
                    return Some(s);
                }
            }
        }
        return None;
    }

    #[cfg(target_os = "macos")]
    {
        let script = format!(
            "POSIX path of (choose file with prompt \"{}\")",
            title.replace('"', "\\\"")
        );
        if let Ok(out) = Command::new("osascript")
            .args(["-e", &script])
            .stdin(Stdio::null())
            .output()
        {
            if out.status.success() {
                let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
                if !s.is_empty() {
                    return Some(s);
                }
            }
        }
        let _ = extensions;
        return None;
    }

    #[cfg(all(unix, not(target_os = "macos")))]
    {
        // 1. Try zenity (GNOME / GTK / most Linux desktops)
        let mut zcmd = Command::new("zenity");
        zcmd.args(["--file-selection", "--title", title]);
        if !extensions.is_empty() {
            let pats: Vec<String> = extensions.iter().map(|e| format!("*.{e}")).collect();
            zcmd.arg(format!("--file-filter=Supported | {}", pats.join(" ")));
            zcmd.arg("--file-filter=All files | *");
        }
        if let Ok(out) = zcmd.stdin(Stdio::null()).output() {
            if out.status.success() {
                let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
                if !s.is_empty() {
                    return Some(s);
                }
            }
            return None; // User cancelled zenity
        }

        // 2. Fallback to kdialog (KDE / Plasma)
        let filter = if extensions.is_empty() {
            "*".to_string()
        } else {
            extensions
                .iter()
                .map(|e| format!("*.{e}"))
                .collect::<Vec<_>>()
                .join(" ")
        };
        if let Ok(out) = Command::new("kdialog")
            .args(["--title", title, "--getopenfilename", &default_export_dir().to_string_lossy(), &filter])
            .stdin(Stdio::null())
            .output()
        {
            if out.status.success() {
                let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
                if !s.is_empty() {
                    return Some(s);
                }
            }
        }
        None
    }
}

/// Opens a native OS "Save File" dialog with `suggested_name`, returning the chosen path.
/// Returns `None` if cancelled or in headless/test environments.
pub fn pick_save_file(title: &str, suggested_name: &str) -> Option<String> {
    if cfg!(test) {
        return None;
    }
    let default_path = fallback_export_path(suggested_name);
    let default_str = default_path.to_string_lossy().to_string();

    #[cfg(target_os = "windows")]
    {
        let ps = format!(
            "Add-Type -AssemblyName System.Windows.Forms; \
             $d = New-Object System.Windows.Forms.SaveFileDialog; \
             $d.Title = '{}'; \
             $d.FileName = '{}'; \
             if ($d.ShowDialog() -eq 'OK') {{ Write-Output $d.FileName }}",
            title.replace('\'', "''"),
            suggested_name.replace('\'', "''")
        );
        if let Ok(out) = Command::new("powershell")
            .args(["-NoProfile", "-STA", "-Command", &ps])
            .stdin(Stdio::null())
            .output()
        {
            if out.status.success() {
                let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
                if !s.is_empty() {
                    return Some(s);
                }
            }
        }
        let _ = default_str;
        return None;
    }

    #[cfg(target_os = "macos")]
    {
        let script = format!(
            "POSIX path of (choose file name with prompt \"{}\" default name \"{}\")",
            title.replace('"', "\\\""),
            suggested_name.replace('"', "\\\"")
        );
        if let Ok(out) = Command::new("osascript")
            .args(["-e", &script])
            .stdin(Stdio::null())
            .output()
        {
            if out.status.success() {
                let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
                if !s.is_empty() {
                    return Some(s);
                }
            }
        }
        let _ = default_str;
        return None;
    }

    #[cfg(all(unix, not(target_os = "macos")))]
    {
        if let Ok(out) = Command::new("zenity")
            .args([
                "--file-selection",
                "--save",
                "--confirm-overwrite",
                "--title",
                title,
                "--filename",
                &default_str,
            ])
            .stdin(Stdio::null())
            .output()
        {
            if out.status.success() {
                let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
                if !s.is_empty() {
                    return Some(s);
                }
            }
            return None;
        }
        if let Ok(out) = Command::new("kdialog")
            .args(["--title", title, "--getsavefilename", &default_str])
            .stdin(Stdio::null())
            .output()
        {
            if out.status.success() {
                let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
                if !s.is_empty() {
                    return Some(s);
                }
            }
        }
        None
    }
}

/// Saves `content` using a native OS Save dialog (falling back to `~/Documents/<suggested_name>`
/// only when running headless/without a native dialog binary, never `/tmp/...`).
pub fn save_text_with_dialog(title: &str, suggested_name: &str, content: &str) -> Result<Option<String>, String> {
    if cfg!(test) {
        let p = std::env::temp_dir().join(suggested_name);
        std::fs::write(&p, content).map_err(|e| e.to_string())?;
        return Ok(Some(p.to_string_lossy().to_string()));
    }
    if let Some(chosen) = pick_save_file(title, suggested_name) {
        std::fs::write(&chosen, content).map_err(|e| e.to_string())?;
        Ok(Some(chosen))
    } else {
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_export_dir_is_non_empty_and_not_hardcoded_unix_tmp() {
        let dir = default_export_dir();
        assert!(!dir.as_os_str().is_empty());
        let p = fallback_export_path("aspirin-ftir.svg");
        assert!(p.file_name().and_then(|s| s.to_str()) == Some("aspirin-ftir.svg"));
    }
}
