# Architecture & Developer Guide

ChiralCraft is engineered as a clean, layered **5-crate pure-Rust workspace** where chemistry algorithms, file codecs, and document state are completely decoupled from the desktop GUI.

## Workspace Crate Dependency Graph

```text
                ┌──────────────────────────┐
                │     chiralcraft-app      │  eframe/egui GUI, 53 tools, 3D viewer,
                │  (GUI + CLI + MCP bin)   │  PubChem/ChEMBL, Studio Dock, CLI, MCP
                └────────────┬─────────────┘
                             │
                ┌────────────▼─────────────┐
                │     chiralcraft-doc      │  DocumentScene, Undo/Redo transactions,
                │   (Scene & Publishing)   │  16 stylesheets, 11 libraries, SVG, LaTeX,
                └──────┬────────────┬──────┘  Markush, Retrosynthesis, ΔG‡, ELN SI
                       │            │
       ┌───────────────▼──┐      ┌──▼───────────────────┐
       │chiralcraft-formats│      │   chiralcraft-chem   │  Valence, aromaticity, CIP,
       │ (22+ Codecs & IO) │      │ (Cheminformatics &   │  NMR/IR/UV/2D, MedChem/HRMS,
       └───────────────┬──┘      │  Force Field Engine) │  18e⁻ CBC, UFF-lite 3D, IUPAC
                       │         └──────────┬───────────┘
                       └──────────┬─────────┘
                        ┌─────────▼──────────┐
                        │ chiralcraft-tables │  118 elements, IUPAC isotopes,
                        │ (Reference Data)   │  valence rules, 537 nicknames
                        └────────────────────┘
```

## Key Design Principles

1. **Zero External Runtime Dependencies**: No Python, RDKit, Java, OpenBabel, or Electron. Everything compiles into a single stripped native binary (~16 MB on Linux x86_64) via `cargo build -p chiralcraft-app --release`.
2. **Transaction-Backed Document Model**: Every canvas mutation runs through `InteractiveDocument::apply_scene_tx`, guaranteeing consistent Undo/Redo history, `is_dirty()` tracking, and clean `Revert`.
3. **Single Source of Truth for GUI, CLI & MCP**: Analysis, NMR, FT-IR, UV–Vis, MedChem, HRMS, Retrosynthesis, 18e⁻ organometallics, Quantum I/O, LaTeX export, Chemical Diff, and Editable SVG/PNG/PDF generation all call the exact same library functions.
4. **Lossless Media Round-Trip**: `to_cdxml()` serializes core CDXML plus `<objecttag Name="chiralcraft_extensions">` (energy profiles, table cells, TLC metadata), which is then embedded verbatim into SVG (`<chiralcraft:document>`), PNG (`tEXt` chunk), and PDF 1.4 (`%CHIRALCRAFT_CDXML_B64:`).
5. **Clean-Room & Safety Hygiene**: Zero proprietary assets or code; no production `.unwrap()`/`panic!` in the desktop app; panic hook writes diagnostics to the user cache/state directory (`~/.cache/chiralcraft/crash.log` or `$XDG_STATE_HOME/chiralcraft/crash.log`) without corrupting saved documents (atomic write-and-rename for saves and SVG round-trips).

## Module Map (`crates/chiralcraft-app/src/`)

| Module | Purpose |
|---|---|
| `main.rs` | `ChiralCraftApp` state, MDI canvas, 53-tool interaction engine, 11 top menus, dialogs, CLI entrypoint |
| `io_bridge.rs` | Unified loader/saver across 22+ formats, SMILES fallback, atomic writes, Wayland/X11 SVG clipboard |
| `config_dock.rs` | Persistent `~/.config/chiralcraft/config.json`, unified Right-Side Studio Dock, Quantum Bridge UI |
| `command_palette.rs` | Fuzzy `:` / `F4` Command Palette + parameterized command parser |
| `viewer_3d.rs` | Interactive 3D orbital viewer, UFF optimizer controls, crystal unit-cell wireframe, distance/angle/dihedral measurement, vibrational mode animation |
| `pro_panels.rs` | ¹H/¹³C NMR viewer + JCAMP-DX overlay, HELM/FASTA builder, BioGraphics & Ring Fill toolbars, Stoichiometry, Diff, Object Inspector |
| `pro_suite.rs` | MedChem/ADMET/PAINS Scorecard, HRMS Isotope Simulator, Isomers/pH/Tautomers/Resonance, Markush SAR, Retrosynthesis, TOC Composer |
| `omni_suite.rs` | Multi-Spectroscopy (FT-IR/UV–Vis/2D NMR), Organometallic 18e⁻ & Polymer/BigSMILES, Energy Profile Builder, ELN/GHS Safety, LaTeX Export, Screening Browser |
| `palettes.rs` | Left Tools Palette, 15 flyout toolbars, Analysis, Properties, Periodic Table, Info & Character Map windows |
| `chem_db.rs` | Zero-dependency JSON parser & live PubChem / ChEMBL REST client + hybrid Structure ↔ Name resolver |
| `chem_ops.rs` | Ring-interior double bond normals, mass fragmentation cleavage, perspective tilt, reaction cleanup/mapping, bond snapping |
| `hotkeys.rs` | 98 context-sensitive atom, bond & generic hotkeys |
| `icons.rs` | Pure-vector toolbar and tool icons |
| `native_dialogs.rs` | Zero-dependency cross-platform native OS Open/Save file dialogs (`zenity`/`kdialog`, `osascript`, PowerShell) |
| `quantum_runner.rs` | Non-blocking local quantum chemistry runner (`xtb`, `orca`, `g16`) and optional NIST WebBook experimental IR fetcher |
| `mcp.rs` | Zero-dependency stdio JSON-RPC 2.0 MCP server (19 tools) |

## Testing & Verification

```bash
# Run the full workspace test suite (unit, round-trip, MCP NDJSON, and headless full-frame egui GUI tests)
cargo test --workspace --all-features

# Build stripped optimized release binary
cargo build -p chiralcraft-app --release

# Install Linux desktop binary, icons, MIME types, and .desktop launcher
./scripts/install-linux-desktop.sh

# Regenerate showcase GUI screenshots directly from GPU framebuffer
./target/release/chiralcraft --screenshot hero docs/images/chiralcraft-hero.png
```
