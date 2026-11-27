//! Command Palette (`:` / `Ctrl+Shift+P`) & Searchable Shortcut Finder
//!
//! Supports both instant fuzzy command search (every major drawing, analysis, theme,
//! 3D, diff, and export action + its shortcut) and parameterized one-liner commands:
//! - `rotate 30` / `rot -45`
//! - `scale 150`
//! - `smiles c1ccccc1O`
//! - `name ibuprofen`
//! - `style nature` / `style acs` / `style rsc` / `style wiley`
//! - `theme light` / `theme dark` / `theme slide`
//! - `cycle circle` / `cycle ellipse` / `cycle starter`
//! - `nmr 1h` / `nmr 13c` / `nmr stamp`
//! - `align left|center|right|top|bottom`
//! - `dpi 150|300|600`
//! - `3d` / `diff` / `inspector` / `db`

use egui::{Color32, Key, RichText};

#[derive(Debug, Clone, PartialEq)]
pub enum PaletteCommandAction {
    None,
    RotateDeg(f64),
    ScalePct(f64),
    InsertSmiles(String),
    InsertName(String),
    ApplyStyle(String),
    SetTheme(&'static str),
    CatalyticCycle(&'static str),
    PredictNmr { carbon13: bool },
    StampNmrLabels { carbon13: bool },
    ClearNmrLabels,
    Open3dViewer,
    OpenDiffCompare,
    OpenObjectInspector,
    OpenChemDb,
    CleanUpStructure,
    CleanUpReaction,
    CheckStructure,
    StructureToName,
    CopyEditableSvg,
    Align(&'static str),
    SetExportDpi(u32),
    InsertSpecialCore(&'static str),
    ToggleStereo,
    ToggleRulers,
    OpenSaveAs(&'static str),
    LinkCaption(&'static str),
    OpenMedChem,
    OpenIsomers(&'static str),
    OpenHrms,
    OpenMarkush,
    OpenRetro,
    OpenPubComposer,
    OpenSpectroscopy(&'static str),
    OpenInoPoly(&'static str),
    OpenEnergyProfile,
    OpenElnSafety,
    OpenLatexExport,
    OpenScreening,
    Optimize3dNow,
    OpenQuantumBridge,
    ToggleStudioDock,
}

pub struct PaletteEntry {
    pub label: &'static str,
    pub shortcut: &'static str,
    pub hint: &'static str,
    pub action: PaletteCommandAction,
}

pub fn builtin_palette_entries() -> Vec<PaletteEntry> {
    vec![
        PaletteEntry { label: "Clean Up Structure", shortcut: "Shift+Ctrl+K", hint: "Beautify active structure preserving stereo", action: PaletteCommandAction::CleanUpStructure },
        PaletteEntry { label: "Clean Up Reaction", shortcut: "Shift+Ctrl+X", hint: "Align reactants, arrow, and products", action: PaletteCommandAction::CleanUpReaction },
        PaletteEntry { label: "Check Structure (Valence & Geometry)", shortcut: "", hint: "Run chemistry validation warnings", action: PaletteCommandAction::CheckStructure },
        PaletteEntry { label: "Structure to Name (Hybrid Offline + PubChem)", shortcut: "Alt+Ctrl+N", hint: "Generate IUPAC/common name caption", action: PaletteCommandAction::StructureToName },
        PaletteEntry { label: "Copy Editable SVG (WordCraft / DeckCraft)", shortcut: "Ctrl+C", hint: "Copy lossless transparent SVG to clipboard", action: PaletteCommandAction::CopyEditableSvg },
        PaletteEntry { label: "Open 3D Molecular & Crystal Viewer", shortcut: "", hint: "Interactive rotatable Ball-and-Stick / CPK / Unit Cell viewer", action: PaletteCommandAction::Open3dViewer },
        PaletteEntry { label: "Compare Documents (Chemical Diff)", shortcut: "", hint: "Chemistry-aware diff between two tabs or files", action: PaletteCommandAction::OpenDiffCompare },
        PaletteEntry { label: "Show Object Inspector (TLC / Tables / 3D / Export)", shortcut: "", hint: "Edit TLC plates, tables, rear opacity, and DPI", action: PaletteCommandAction::OpenObjectInspector },
        PaletteEntry { label: "Search PubChem & ChEMBL Database", shortcut: "Shift+Ctrl+F", hint: "Live online compound search and 2D import", action: PaletteCommandAction::OpenChemDb },
        PaletteEntry { label: "MedChem, ADMET & PAINS Scorecard", shortcut: "", hint: "Lipinski, CNS MPO, QED, Fsp3, SA, PAINS alerts", action: PaletteCommandAction::OpenMedChem },
        PaletteEntry { label: "HRMS Isotope Envelope & Adduct Simulator", shortcut: "", hint: "True isotope distribution, ESI+/ESI- adducts & ppm error", action: PaletteCommandAction::OpenHrms },
        PaletteEntry { label: "pH Microspecies & Isoelectric Point (pI)", shortcut: "", hint: "Physiological pH protonation state slider", action: PaletteCommandAction::OpenIsomers("ph") },
        PaletteEntry { label: "Enumerate Tautomers / Stereoisomers / Resonance", shortcut: "", hint: "Prototropic tautomers, 2^n stereoisomers, resonance & electron push", action: PaletteCommandAction::OpenIsomers("tautomers") },
        PaletteEntry { label: "Markush R-Group Library & SAR Table", shortcut: "", hint: "Enumerate R1/R2 analogs with CSV/SDF export", action: PaletteCommandAction::OpenMarkush },
        PaletteEntry { label: "Retrosynthetic Disconnection Advisor", shortcut: "", hint: "Suggest strategic bonds & stamp ⇒ synthons", action: PaletteCommandAction::OpenRetro },
        PaletteEntry { label: "Publication TOC Graphic & Scheme Composer", shortcut: "", hint: "Journal TOC frames, arrow conditions & serpentine wrap", action: PaletteCommandAction::OpenPubComposer },
        PaletteEntry { label: "Multi-Spectroscopy Studio (FT-IR / UV–Vis / 2D NMR)", shortcut: "", hint: "Simulate IR %T, UV λmax, COSY, HSQC & HMBC", action: PaletteCommandAction::OpenSpectroscopy("ir") },
        PaletteEntry { label: "Organometallic 18-Electron & Green CBC Analyzer", shortcut: "", hint: "Oxidation state, d^n, 18e-/16e- count & catalysts", action: PaletteCommandAction::OpenInoPoly("metal") },
        PaletteEntry { label: "Polymer Mn / Mw / Dispersity & BigSMILES", shortcut: "", hint: "Repeat unit calculator & copolymer BigSMILES", action: PaletteCommandAction::OpenInoPoly("polymer") },
        PaletteEntry { label: "Reaction Energy Profile (ΔG‡) & Eyring Kinetics", shortcut: "", hint: "Build free-energy diagram + t1/2 & Boltzmann %ee", action: PaletteCommandAction::OpenEnergyProfile },
        PaletteEntry { label: "ELN Supporting Info & Lab Safety / GHS Handbook", shortcut: "", hint: "60+ solvent/reagent hazards, liquid mL & ACS SI generator", action: PaletteCommandAction::OpenElnSafety },
        PaletteEntry { label: "Export LaTeX (\\chemfig / \\schemestart / TikZ)", shortcut: "", hint: "Compilable LaTeX & TikZ vector code", action: PaletteCommandAction::OpenLatexExport },
        PaletteEntry { label: "Multi-Compound SDF / SMILES Screening Browser", shortcut: "", hint: "Filter compound libraries by MW, LogP, MPO & PAINS", action: PaletteCommandAction::OpenScreening },
        PaletteEntry { label: "Optimize 3D Geometry (UFF Force-Field)", shortcut: "", hint: "Run molecular mechanics energy minimization", action: PaletteCommandAction::Optimize3dNow },
        PaletteEntry { label: "Quantum Chemistry Bridge (ORCA / Gaussian / xTB)", shortcut: "", hint: "Export DFT decks & parse .out/.log HOMO-LUMO/IR", action: PaletteCommandAction::OpenQuantumBridge },
        PaletteEntry { label: "Toggle Right-Side Studio Dock", shortcut: "", hint: "Show or hide unified tabbed chemistry sidebar", action: PaletteCommandAction::ToggleStudioDock },
        PaletteEntry { label: "Predict 1H-NMR Spectrum", shortcut: "", hint: "400 MHz proton NMR viewer", action: PaletteCommandAction::PredictNmr { carbon13: false } },
        PaletteEntry { label: "Predict 13C-NMR Spectrum", shortcut: "", hint: "100 MHz carbon-13 NMR viewer", action: PaletteCommandAction::PredictNmr { carbon13: true } },
        PaletteEntry { label: "Stamp 1H-NMR Shift Labels on Canvas", shortcut: "", hint: "Place δ (ppm) badges beside atoms", action: PaletteCommandAction::StampNmrLabels { carbon13: false } },
        PaletteEntry { label: "Stamp 13C-NMR Shift Labels on Canvas", shortcut: "", hint: "Place 13C δ (ppm) badges beside carbons", action: PaletteCommandAction::StampNmrLabels { carbon13: true } },
        PaletteEntry { label: "Clear Stamped NMR Shift Labels", shortcut: "", hint: "Remove δ shift captions from canvas", action: PaletteCommandAction::ClearNmrLabels },
        PaletteEntry { label: "Theme: Classic Light (White Canvas)", shortcut: "", hint: "Classic light desktop chrome + white sheet", action: PaletteCommandAction::SetTheme("light") },
        PaletteEntry { label: "Theme: Dark UI + White Canvas", shortcut: "", hint: "Dark chrome with print-ready white sheet", action: PaletteCommandAction::SetTheme("dark") },
        PaletteEntry { label: "Theme: Dark UI + Dark Slide Canvas", shortcut: "", hint: "Dark slide canvas + white-ink SVG export", action: PaletteCommandAction::SetTheme("slide") },
        PaletteEntry { label: "Catalytic Cycle: Insert 4-Step Pd Starter", shortcut: "", hint: "Insert editable 4-intermediate catalytic cycle", action: PaletteCommandAction::CatalyticCycle("starter") },
        PaletteEntry { label: "Catalytic Cycle: Arrange Selected as Circle", shortcut: "", hint: "Arrange intermediates in a circle with curved arrows", action: PaletteCommandAction::CatalyticCycle("circle") },
        PaletteEntry { label: "Catalytic Cycle: Arrange Selected as Ellipse", shortcut: "", hint: "Arrange intermediates in an ellipse with curved arrows", action: PaletteCommandAction::CatalyticCycle("ellipse") },
        PaletteEntry { label: "Insert Core: Haworth β-D-Glucopyranose", shortcut: "", hint: "Pyranose Haworth projection", action: PaletteCommandAction::InsertSpecialCore("pyranose") },
        PaletteEntry { label: "Insert Core: Haworth β-D-Ribofuranose", shortcut: "", hint: "Furanose Haworth projection", action: PaletteCommandAction::InsertSpecialCore("furanose") },
        PaletteEntry { label: "Insert Core: Porphine Macrocycle", shortcut: "", hint: "24-atom porphine core", action: PaletteCommandAction::InsertSpecialCore("porphine") },
        PaletteEntry { label: "Insert Core: 18-Crown-6", shortcut: "", hint: "Crown ether macrocycle", action: PaletteCommandAction::InsertSpecialCore("crown") },
        PaletteEntry { label: "Insert Core: Adamantane Cage (3D)", shortcut: "", hint: "Tetrahedral diamondoid cage with bridge gaps", action: PaletteCommandAction::InsertSpecialCore("adamantane") },
        PaletteEntry { label: "Insert Core: Norbornane", shortcut: "", hint: "Bicyclo[2.2.1]heptane", action: PaletteCommandAction::InsertSpecialCore("norbornane") },
        PaletteEntry { label: "Attach Live Linked Formula Caption", shortcut: "", hint: "Auto-updates formula when molecule is edited", action: PaletteCommandAction::LinkCaption("formula") },
        PaletteEntry { label: "Attach Live Linked IUPAC Name Caption", shortcut: "", hint: "Auto-updates name when molecule is edited", action: PaletteCommandAction::LinkCaption("name") },
        PaletteEntry { label: "Export Editable Vector PDF...", shortcut: "", hint: "Save lossless editable PDF 1.4", action: PaletteCommandAction::OpenSaveAs("export.pdf") },
        PaletteEntry { label: "Export Editable PNG (300/600 DPI)...", shortcut: "", hint: "Save PNG with embedded ChiralCraft CDXML", action: PaletteCommandAction::OpenSaveAs("export.png") },
        PaletteEntry { label: "Export Editable SVG...", shortcut: "", hint: "Save WordCraft/DeckCraft editable SVG", action: PaletteCommandAction::OpenSaveAs("export.svg") },
        PaletteEntry { label: "Toggle Stereochemistry Labels (R/S/E/Z)", shortcut: "", hint: "Show or hide CIP labels", action: PaletteCommandAction::ToggleStereo },
        PaletteEntry { label: "Toggle Rulers", shortcut: "F11", hint: "Show or hide cm rulers", action: PaletteCommandAction::ToggleRulers },
        PaletteEntry { label: "Align Left Edges", shortcut: "", hint: "Align selected structures left", action: PaletteCommandAction::Align("left") },
        PaletteEntry { label: "Align Horizontal Centers", shortcut: "", hint: "Center selected structures horizontally", action: PaletteCommandAction::Align("center") },
        PaletteEntry { label: "Align Top Edges", shortcut: "", hint: "Align selected structures top", action: PaletteCommandAction::Align("top") },
    ]
}

/// Parses parameterized one-line commands like `rotate 30`, `scale 150`, `smiles CCO`,
/// `name aspirin`, `style nature`, `theme dark`, `dpi 600`, `nmr 13c`.
pub fn parse_direct_palette_command(input: &str) -> Option<PaletteCommandAction> {
    let trimmed = input.trim().trim_start_matches(':').trim();
    if trimmed.is_empty() {
        return None;
    }
    let mut parts = trimmed.splitn(2, char::is_whitespace);
    let cmd = parts.next()?.to_ascii_lowercase();
    let arg = parts.next().unwrap_or("").trim();
    match cmd.as_str() {
        "rotate" | "rot" => arg
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
            .map(PaletteCommandAction::RotateDeg),
        "scale" => arg
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite() && (1.0..=10_000.0).contains(v))
            .map(PaletteCommandAction::ScalePct),
        "smiles" | "smi" if !arg.is_empty() => Some(PaletteCommandAction::InsertSmiles(arg.to_string())),
        "name" | "iupac" if !arg.is_empty() => Some(PaletteCommandAction::InsertName(arg.to_string())),
        "style" | "stylesheet" if !arg.is_empty() => Some(PaletteCommandAction::ApplyStyle(arg.to_string())),
        "theme" => match arg.to_ascii_lowercase().as_str() {
            "light" | "classic" => Some(PaletteCommandAction::SetTheme("light")),
            "dark" | "chrome" => Some(PaletteCommandAction::SetTheme("dark")),
            "slide" | "darkslide" | "dark-slide" => Some(PaletteCommandAction::SetTheme("slide")),
            _ => None,
        },
        "cycle" => match arg.to_ascii_lowercase().as_str() {
            "circle" => Some(PaletteCommandAction::CatalyticCycle("circle")),
            "ellipse" => Some(PaletteCommandAction::CatalyticCycle("ellipse")),
            "" | "starter" | "4" => Some(PaletteCommandAction::CatalyticCycle("starter")),
            _ => None,
        },
        "nmr" => match arg.to_ascii_lowercase().as_str() {
            "1h" | "h" | "" => Some(PaletteCommandAction::PredictNmr { carbon13: false }),
            "13c" | "c" => Some(PaletteCommandAction::PredictNmr { carbon13: true }),
            "stamp" | "stamp1h" => Some(PaletteCommandAction::StampNmrLabels { carbon13: false }),
            "stamp13c" => Some(PaletteCommandAction::StampNmrLabels { carbon13: true }),
            "clear" => Some(PaletteCommandAction::ClearNmrLabels),
            _ => None,
        },
        "3d" | "view3d" => Some(PaletteCommandAction::Open3dViewer),
        "diff" | "compare" => Some(PaletteCommandAction::OpenDiffCompare),
        "inspector" | "tlc" => Some(PaletteCommandAction::OpenObjectInspector),
        "db" | "pubchem" | "chembl" => Some(PaletteCommandAction::OpenChemDb),
        "medchem" | "pains" | "admet" | "lipinski" => Some(PaletteCommandAction::OpenMedChem),
        "hrms" | "ms" | "isotope" | "adducts" => Some(PaletteCommandAction::OpenHrms),
        "ph" | "pi" | "microspecies" => Some(PaletteCommandAction::OpenIsomers("ph")),
        "tautomers" | "tautomer" => Some(PaletteCommandAction::OpenIsomers("tautomers")),
        "stereoisomers" | "stereo" => Some(PaletteCommandAction::OpenIsomers("stereo")),
        "resonance" => Some(PaletteCommandAction::OpenIsomers("resonance")),
        "push" | "mechanism" => Some(PaletteCommandAction::OpenIsomers("push")),
        "markush" | "rgroup" | "sar" => Some(PaletteCommandAction::OpenMarkush),
        "retro" | "retrosynthesis" => Some(PaletteCommandAction::OpenRetro),
        "toc" | "scheme" | "composer" => Some(PaletteCommandAction::OpenPubComposer),
        "ir" | "ftir" | "spec" | "spectra" => Some(PaletteCommandAction::OpenSpectroscopy("ir")),
        "uv" | "uvvis" => Some(PaletteCommandAction::OpenSpectroscopy("uv")),
        "cosy" | "hsqc" | "hmbc" | "2dnmr" => Some(PaletteCommandAction::OpenSpectroscopy("2d")),
        "18e" | "metal" | "organometallic" => Some(PaletteCommandAction::OpenInoPoly("metal")),
        "polymer" | "bigsmiles" => Some(PaletteCommandAction::OpenInoPoly("polymer")),
        "energy" | "eyring" | "dg" | "kinetics" => Some(PaletteCommandAction::OpenEnergyProfile),
        "eln" | "ghs" | "safety" | "si" => Some(PaletteCommandAction::OpenElnSafety),
        "latex" | "chemfig" | "tikz" | "tex" => Some(PaletteCommandAction::OpenLatexExport),
        "screen" | "sdf" | "library" => Some(PaletteCommandAction::OpenScreening),
        "opt3d" | "uff" | "minimize" => Some(PaletteCommandAction::Optimize3dNow),
        "orca" | "gaussian" | "xtb" | "dft" | "quantum" => Some(PaletteCommandAction::OpenQuantumBridge),
        "dock" | "sidebar" => Some(PaletteCommandAction::ToggleStudioDock),
        "dpi" => arg.parse::<u32>().ok().map(PaletteCommandAction::SetExportDpi),
        "align" => match arg.to_ascii_lowercase().as_str() {
            "left" => Some(PaletteCommandAction::Align("left")),
            "center" | "hcenter" => Some(PaletteCommandAction::Align("center")),
            "right" => Some(PaletteCommandAction::Align("right")),
            "top" => Some(PaletteCommandAction::Align("top")),
            "bottom" => Some(PaletteCommandAction::Align("bottom")),
            _ => None,
        },
        "export" => match arg.to_ascii_lowercase().as_str() {
            "pdf" => Some(PaletteCommandAction::OpenSaveAs("export.pdf")),
            "png" => Some(PaletteCommandAction::OpenSaveAs("export.png")),
            "svg" => Some(PaletteCommandAction::OpenSaveAs("export.svg")),
            "cdxml" => Some(PaletteCommandAction::OpenSaveAs("export.cdxml")),
            _ => None,
        },
        _ => None,
    }
}

pub struct CommandPaletteState {
    pub open: bool,
    pub query: String,
    pub selected_idx: usize,
    pub just_opened: bool,
}

impl Default for CommandPaletteState {
    fn default() -> Self {
        Self {
            open: false,
            query: String::new(),
            selected_idx: 0,
            just_opened: false,
        }
    }
}

impl CommandPaletteState {
    pub fn open_fresh(&mut self) {
        self.open = true;
        self.query.clear();
        self.selected_idx = 0;
        self.just_opened = true;
    }
}

pub fn show_command_palette(
    ctx: &egui::Context,
    state: &mut CommandPaletteState,
) -> PaletteCommandAction {
    if !state.open {
        return PaletteCommandAction::None;
    }
    let mut chosen = PaletteCommandAction::None;
    let mut close_now = false;

    if ctx.input(|i| i.key_pressed(Key::Escape)) {
        state.open = false;
        return PaletteCommandAction::None;
    }

    let entries = builtin_palette_entries();
    let q_low = state.query.trim().trim_start_matches(':').to_ascii_lowercase();
    let direct = parse_direct_palette_command(&state.query);
    let filtered: Vec<&PaletteEntry> = entries
        .iter()
        .filter(|e| {
            q_low.is_empty()
                || e.label.to_ascii_lowercase().contains(&q_low)
                || e.hint.to_ascii_lowercase().contains(&q_low)
                || e.shortcut.to_ascii_lowercase().contains(&q_low)
        })
        .collect();

    if !filtered.is_empty() {
        state.selected_idx = state.selected_idx.min(filtered.len() - 1);
    } else {
        state.selected_idx = 0;
    }

    if ctx.input(|i| i.key_pressed(Key::ArrowDown)) && !filtered.is_empty() {
        state.selected_idx = (state.selected_idx + 1) % filtered.len();
    }
    if ctx.input(|i| i.key_pressed(Key::ArrowUp)) && !filtered.is_empty() {
        state.selected_idx = if state.selected_idx == 0 {
            filtered.len() - 1
        } else {
            state.selected_idx - 1
        };
    }

    egui::Window::new("Command Palette (:  or  F4)")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_TOP, [0.0, 72.0])
        .default_width(520.0)
        .min_width(520.0)
        .show(ctx, |ui| {
            let te = ui.add(
                egui::TextEdit::singleline(&mut state.query)
                    .hint_text("Type a command (e.g. rotate 30, smiles CCO, name aspirin, style nature, 3d, diff)...")
                    .desired_width(f32::INFINITY),
            );
            if state.just_opened {
                te.request_focus();
                state.just_opened = false;
            }
            let enter = ctx.input(|i| i.key_pressed(Key::Enter));
            if let Some(ref direct_act) = direct {
                ui.add_space(4.0);
                if ui
                    .selectable_label(
                        true,
                        RichText::new(format!("Run command: {}", state.query.trim())).strong(),
                    )
                    .clicked()
                    || enter
                {
                    chosen = direct_act.clone();
                    close_now = true;
                }
            }
            ui.separator();
            egui::ScrollArea::vertical().max_height(290.0).show(ui, |ui| {
                for (i, entry) in filtered.iter().enumerate() {
                    let is_sel = direct.is_none() && i == state.selected_idx;
                    let resp = ui.horizontal(|ui| {
                        let btn = ui.selectable_label(is_sel, RichText::new(entry.label).strong());
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if !entry.shortcut.is_empty() {
                                ui.label(RichText::new(entry.shortcut).monospace().color(Color32::from_rgb(70, 120, 185)));
                            }
                            ui.label(RichText::new(entry.hint).small().weak());
                        });
                        btn
                    });
                    if resp.inner.clicked() || (enter && is_sel && chosen == PaletteCommandAction::None) {
                        chosen = entry.action.clone();
                        close_now = true;
                    }
                }
                if filtered.is_empty() && direct.is_none() {
                    ui.weak("No matching command. Try: rotate 45, scale 150, smiles c1ccccc1, name caffeine, style nature, 3d, diff");
                }
            });
        });

    if close_now {
        state.open = false;
    }
    chosen
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_command_palette_direct_parser() {
        assert_eq!(parse_direct_palette_command(":rotate 45"), Some(PaletteCommandAction::RotateDeg(45.0)));
        assert_eq!(parse_direct_palette_command("scale 150"), Some(PaletteCommandAction::ScalePct(150.0)));
        assert_eq!(parse_direct_palette_command("smiles c1ccccc1"), Some(PaletteCommandAction::InsertSmiles("c1ccccc1".into())));
        assert_eq!(parse_direct_palette_command("name aspirin"), Some(PaletteCommandAction::InsertName("aspirin".into())));
        assert_eq!(parse_direct_palette_command("style nature"), Some(PaletteCommandAction::ApplyStyle("nature".into())));
        assert_eq!(parse_direct_palette_command("theme slide"), Some(PaletteCommandAction::SetTheme("slide")));
        assert_eq!(parse_direct_palette_command("nmr stamp"), Some(PaletteCommandAction::StampNmrLabels { carbon13: false }));
        assert_eq!(parse_direct_palette_command("dpi 600"), Some(PaletteCommandAction::SetExportDpi(600)));
        assert_eq!(parse_direct_palette_command("3d"), Some(PaletteCommandAction::Open3dViewer));
        assert_eq!(parse_direct_palette_command("diff"), Some(PaletteCommandAction::OpenDiffCompare));
    }
}
