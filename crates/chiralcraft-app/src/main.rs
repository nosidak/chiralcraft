//! ChiralCraft — Native Desktop GUI (`eframe`/`egui`) & CLI Application
//!
//! Styled and wired for 1:1 ChiralCraft visual and functional parity:
//! - Exact Light Windows theme (`#F0F0F0` chrome, `#A0A0A0` MDI desk, steel-blue `#B9D1EA` frames)
//! - Exact 11-menu bar (`File`, `Edit`, `View`, `Object`, `Structure`, `Text`, `Curves`, `Colors`, `Search`, `Window`, `Help`)
//! - Exact unified General + Style top toolbar (`p08-opened.png`) with custom vector icons (zero emojis)
//! - Exact 18x2 iconic `Tools` palette (`tools-20261118.png`) with tear-off flyout sub-palettes & live vector `.ctp` previews
//! - Floating `Analysis`, `Chemical Properties`, 18-col `Periodic Table`, `Info`, and `Character Map` palettes
//! - Ring-interior shortened double bonds, inline on-canvas atom editing with Unicode subscripts,
//!   freehand polygon Lasso selection, 8-handle selection box, Mass Fragmentation, Structure Perspective,
//!   Reaction Atom Mapping, Clean Up Reaction, Multi-document MDI tabs, 537 official nicknames, and 12 formats

mod chem_db;
mod chem_ops;
mod command_palette;
mod config_dock;
mod hotkeys;
mod icons;
mod io_bridge;
mod mcp;
mod native_dialogs;
mod omni_suite;
mod palettes;
mod pro_panels;
mod pro_suite;
mod quantum_runner;
mod viewer_3d;

use std::collections::{BTreeMap, HashMap};
use std::env;
use std::fs;
use std::path::PathBuf;

use chiralcraft_chem::{
    analyze_organometallic_complexes, compute_medchem_report, estimate_chemical_properties,
    get_element, get_element_by_symbol, is_fallback_chemical_name, molecule_to_name,
    optimize_geometry_3d, perceive_stereochemistry, predict_13c_nmr, predict_1h_nmr,
    predict_ir_spectrum, predict_uv_vis, simulate_hrms, BondDisplay, CcAromaticityInfo, CcVal,
    ChemAtom, ChemBond, ChemMolecule, NameDictionary, NmrPrediction,
};
use chiralcraft_doc::{
    align_molecules, arrange_catalytic_cycle, build_adamantane_2d, build_crown_ether_18_6,
    build_haworth_furanose, build_haworth_pyranose, build_norbornane_2d, build_porphine_core,
    build_stoichiometry_rows, center_on_page, clean_up_structure_preserving_stereo,
    compare_documents, compute_3d_atom_opacities, distribute_molecules, document_to_chemfig_scheme,
    expand_generic_structure_from_cdxml, extract_ring_polygons, find_bond_crossings,
    flatten_molecule_to_2d, flip_atoms_preserving_stereo, format_linked_caption_text,
    insert_starter_catalytic_cycle, insert_toc_frame, join_atoms, join_molecules_at_atoms,
    lift_molecule_2d_to_3d, load_screening_library, resolve_reverse_abbreviation,
    split_segment_for_crossings, stamp_arrow_conditions, stamp_retrosynthesis_scheme,
    stamp_template, suggest_disconnections, svg_contains_chiralcraft_document,
    toggle_hydrogen_display, wrap_serpentine_reaction_scheme, AlignMode, AromaticCircle,
    BuiltinTemplate, CaptionNode, CopyBackground, DistributeAxis, DocumentDiff, DocumentSettings,
    EnhancedStereoKind, ExportOptions, FlipAxis, InteractiveDocument, LabelJustification,
    LinkedCaption, LinkedCaptionKind, NicknameCollision, SceneArrow, SceneCurve, SceneGelPlate,
    SceneGraphic, SceneTable, SceneTlcPlate, StoichiometryRow, TemplateLibrary,
};
use chiralcraft_formats::{
    clean_up_biopolymer, decode_png_rgba, export_quantum_input, helm_to_molecule,
    inchi_to_inchikey, molecule_to_3mf_xml, molecule_to_fasta, molecule_to_helm, molecule_to_inchi,
    molecule_to_molfile_v2000, molecule_to_molfile_v3000, molecule_to_sln,
    molecule_to_smiles, molecules_to_sdf, molfile_to_molecule, molfile_v3000_to_molecule,
    monomer_to_all_atom_molecule, parse_experimental_ir_file, parse_jcamp_dx, parse_quantum_output,
    sdf_to_molecules,
    sequence_to_biopolymer_molecule, smiles_to_molecule, CdxmlDocument, JcampSpectrum,
    MonomerLibrary, PolymerKind, QuantumJobSpec, QuantumProgram, SdfRecord,
};

const APP_ICON_PNG: &[u8] = include_bytes!("../../../assets/chiralcraft-icon-128.png");

fn load_app_icon_data() -> Option<egui::IconData> {
    let (rgba, width, height) = decode_png_rgba(APP_ICON_PNG)?;
    Some(egui::IconData { rgba, width, height })
}
use egui::{Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};

use crate::chem_db::{
    hit_to_molecule, show_chem_db_window, ChemDbAction, ChemDbState, DbHit, DbSearchMode,
    LiveNameResolver,
};
use crate::chem_ops::{
    apply_perspective_tilt, auto_map_reaction_atoms, clean_up_reaction, cpk_element_color,
    fragment_molecule_by_line, inchi_to_molecule, list_directory_entries, point_in_polygon,
    ring_bond_inward_normals, rotate_selected_molecules, scale_selected_molecules,
    snap_bond_endpoint,
};
use crate::hotkeys::{
    apply_atom_hotkey, apply_bond_hotkey, apply_generic_hotkey, DoubleBondAlign, HotkeyOutcome,
};
use crate::icons::{
    draw_steel_blue_gradient, format_chemical_formula, toolbar_icon_button, Tool, ToolbarIcon,
};
use crate::palettes::{
    show_analysis_window, show_character_map_window, show_info_window, show_periodic_table_window,
    show_properties_window, show_tools_palette, ActiveFlyout, AnalysisAction, AnalysisPanelState,
    PaletteCommand, PropsPanelState, CD_DESK_GRAY, CD_WINDOW_BG,
};
use crate::command_palette::{show_command_palette, CommandPaletteState, PaletteCommandAction};
use crate::config_dock::{
    show_quantum_window, show_right_studio_dock, DockQuantumAction, DockTab, QuantumStudioState,
    UserConfig,
};
use crate::io_bridge::{
    atomic_write_bytes, atomic_write_text, copy_svg_to_os_clipboard, load_document_any,
    load_document_any_with_cell, load_document_or_smiles, paste_svg_from_os_clipboard,
    save_document_any, save_document_any_styled,
};
use crate::omni_suite::{
    show_eln_safety_window, show_energy_profile_window, show_ino_poly_window,
    show_latex_export_window, show_screening_window, show_spectroscopy_window, InoPolyTab,
    IrEngineMode, OmniSuiteAction, OmniSuiteState, SpecTab,
};
use crate::pro_panels::{
    show_biographics_toolbar, show_diff_compare_window, show_helm_monomer_window, show_nmr_window,
    show_object_inspector_window, show_ring_fill_toolbar, show_stoichiometry_window,
    show_structure_toolbar, DiffCompareAction, ProPanelAction,
};
use crate::pro_suite::{
    show_hrms_window, show_isomers_window, show_markush_window, show_medchem_window,
    show_publication_composer_window, show_retrosynthesis_window, IsomerTab, ProSuiteAction,
    ProSuiteState,
};
use crate::quantum_runner::{NistIrFetchState, QcRunBackend, QuantumRunnerState};
use crate::viewer_3d::{show_3d_viewer_window, VibrationalAnimState, Viewer3dAction, Viewer3dState};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum UiTheme {
    ClassicLight,
    DarkChrome,
    DarkSlide,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ActiveDialog {
    None,
    OpenFile,
    InsertFile,
    InsertObject,
    SaveAs,
    ApplyStylesheetFile,
    PasteSmiles,
    PasteNameAsStructure,
    PasteInchi,
    PasteMolText,
    PasteCdxml,
    PasteEditableSvg,
    PasteBiopolymer(PolymerKind),
    PasteHelm,
    AddToNameDict,
    ScaleDialog,
    RotateDialog,
    PageSetup,
    DocumentSettings,
    DocumentProperties,
    ObjectSettings,
    Preferences,
    ListNicknames,
    DefineNickname,
    AtomProperties { mol_idx: usize, atom_idx: usize },
    BondProperties { mol_idx: usize, bond_idx: usize },
    BracketProperties,
    RLogicQuery,
    LineSpacing,
    FontSizeOther,
    CheckStructureDialog,
    ColorPicker,
    AttachData,
    Annotate,
    FindElement,
    EditCaption { cap_idx: usize },
    NewCaption { x_tenths: i32, y_tenths: i32 },
    HelpTopic { title: String, body: String },
    Shortcuts,
    Licensing,
    About,
    /// Active tab has unsaved changes and its close was requested.
    ConfirmCloseTab,
    /// The window close was requested while documents have unsaved changes.
    ConfirmExit,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TextAlign {
    Left,
    Center,
    Right,
    Justify,
    StackedAbove,
    Automatic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MdiDragMode {
    None,
    MoveTitle,
    ResizeE,
    ResizeS,
    ResizeSE,
    ResizeW,
    ResizeN,
    ScrollV,
    ScrollH,
}

static NEXT_TAB_ID: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);

fn crash_log_path() -> PathBuf {
    let base = std::env::var_os("XDG_STATE_HOME")
        .or_else(|| std::env::var_os("XDG_CACHE_HOME"))
        .or_else(|| std::env::var_os("LOCALAPPDATA"))
        .map(PathBuf::from)
        .or_else(|| {
            std::env::var_os("HOME")
                .or_else(|| std::env::var_os("USERPROFILE"))
                .map(|h| PathBuf::from(h).join(".cache"))
        })
        .unwrap_or_else(std::env::temp_dir);
    base.join("chiralcraft").join("crash.log")
}

struct OpenDocTab {
    tab_id: u64,
    doc: InteractiveDocument,
    file_path: Option<String>,
    meta_title: String,
    meta_author: String,
    meta_subject: String,
    meta_comments: String,
    paper_size: String,
    landscape: bool,
    page_width_pt: f64,
    page_height_pt: f64,
    selected_mols: Vec<usize>,
    groups: Vec<Vec<usize>>,
    reaction_atom_maps: HashMap<(usize, usize), u32>,
    reaction_auto_numbers: HashMap<usize, String>,
    enhanced_stereo: HashMap<(usize, usize), EnhancedStereoKind>,
    next_or_group: u32,
    next_and_group: u32,
    ring_fills: HashMap<usize, Color32>,
    mol_highlights: HashMap<usize, Color32>,
    gel_plates: Vec<SceneGelPlate>,
    double_bond_align: HashMap<(usize, usize), DoubleBondAlign>,
    show_atom_numbers: bool,
    bracket_usage: String,
    bracket_repeat_pattern: String,
    bracket_repeat_count: String,
    next_rxn_map_id: u32,
    annotations: BTreeMap<String, String>,
    attached_data: Vec<(String, String)>,
    linked_captions: Vec<LinkedCaption>,
    win_pos: Vec2,
    win_size: Vec2,
    scroll: Vec2,
    maximized: bool,
    minimized: bool,
    /// Side-table state observed at each undo-history position (`side_history[k]` = state at cursor `k`),
    /// so Undo/Redo/Revert restore fills, maps and links together with the scene they index into.
    side_history: Vec<Option<MolSideTables>>,
    side_cursor: usize,
}

/// Molecule/atom/bond-indexed tables kept beside the document scene (not part of `DocumentScene`).
#[derive(Clone)]
struct MolSideTables {
    groups: Vec<Vec<usize>>,
    reaction_atom_maps: HashMap<(usize, usize), u32>,
    reaction_auto_numbers: HashMap<usize, String>,
    enhanced_stereo: HashMap<(usize, usize), EnhancedStereoKind>,
    ring_fills: HashMap<usize, Color32>,
    mol_highlights: HashMap<usize, Color32>,
    double_bond_align: HashMap<(usize, usize), DoubleBondAlign>,
    linked_captions: Vec<LinkedCaption>,
    /// Selection is molecule-indexed too: restoring it keeps Undo from selecting a neighbour.
    selected_mols: Vec<usize>,
}

impl OpenDocTab {
    fn side_tables(&self) -> MolSideTables {
        MolSideTables {
            groups: self.groups.clone(),
            reaction_atom_maps: self.reaction_atom_maps.clone(),
            reaction_auto_numbers: self.reaction_auto_numbers.clone(),
            enhanced_stereo: self.enhanced_stereo.clone(),
            ring_fills: self.ring_fills.clone(),
            mol_highlights: self.mol_highlights.clone(),
            double_bond_align: self.double_bond_align.clone(),
            linked_captions: self.linked_captions.clone(),
            selected_mols: self.selected_mols.clone(),
        }
    }

    fn set_side_tables(&mut self, s: MolSideTables) {
        self.groups = s.groups;
        self.reaction_atom_maps = s.reaction_atom_maps;
        self.reaction_auto_numbers = s.reaction_auto_numbers;
        self.enhanced_stereo = s.enhanced_stereo;
        self.ring_fills = s.ring_fills;
        self.mol_highlights = s.mol_highlights;
        self.double_bond_align = s.double_bond_align;
        self.linked_captions = s.linked_captions;
        self.selected_mols = s.selected_mols;
    }

    /// Record the side tables for the current undo position. Run once per frame (before input is
    /// handled), so the entry at the pre-action cursor still holds the pre-action tables.
    fn sync_side_history(&mut self) {
        let c = self.doc.undo_cursor();
        if c < self.side_cursor {
            // History was reset/replaced outside Undo (e.g. a new document); start over.
            self.side_history.clear();
        } else if c > self.side_cursor {
            // New transaction(s) recorded: the old redo branch is gone.
            self.side_history.truncate(self.side_cursor + 1);
        }
        if self.side_history.len() <= c {
            self.side_history.resize(c + 1, None);
        }
        self.side_history[c] = Some(self.side_tables());
        self.side_cursor = c;
    }

    /// `Edit > Undo` / `Edit > Redo` including the side tables. Returns whether the document changed.
    fn undo_redo(&mut self, redo: bool) -> bool {
        self.sync_side_history();
        let changed = if redo { self.doc.redo() } else { self.doc.undo() };
        if changed {
            let c = self.doc.undo_cursor();
            if let Some(Some(s)) = self.side_history.get(c).cloned() {
                self.set_side_tables(s);
            }
            self.side_cursor = c;
        }
        changed
    }

    /// `File > Revert` including the side tables as they were at the last save.
    fn revert(&mut self) -> bool {
        self.sync_side_history();
        let saved = self.doc.saved_undo_cursor().and_then(|sc| self.side_history.get(sc).cloned().flatten());
        let changed = self.doc.revert();
        if changed {
            if let Some(s) = saved {
                self.set_side_tables(s);
            }
            self.side_history.clear();
            self.side_cursor = 0;
            self.sync_side_history();
        }
        changed
    }

    /// Re-key molecule-indexed side tables after `removed` molecule indices were deleted from the
    /// scene, so later molecules keep their fills/maps/linked captions instead of inheriting a neighbour's.
    fn remap_after_mol_removal(&mut self, removed: &[usize]) {
        if removed.is_empty() {
            return;
        }
        let mut rm = removed.to_vec();
        rm.sort_unstable();
        rm.dedup();
        self.remap_mol_indices(|mi| match rm.binary_search(&mi) {
            Ok(_) => None,
            Err(shift) => Some(mi - shift),
        });
    }

    /// Re-key side tables after molecule `from` was removed and re-inserted at index `to`.
    fn remap_after_mol_move(&mut self, from: usize, to: usize) {
        self.remap_mol_indices(|mi| {
            Some(if mi == from {
                to
            } else if from < to && mi > from && mi <= to {
                mi - 1
            } else if to < from && mi >= to && mi < from {
                mi + 1
            } else {
                mi
            })
        });
    }

    /// Apply a molecule-index mapping (`None` = molecule gone) to every molecule-indexed side table.
    fn remap_mol_indices(&mut self, map: impl Fn(usize) -> Option<usize> + Copy) {
        fn rekey<K: Copy + Eq + std::hash::Hash, V>(
            m: &mut HashMap<K, V>,
            f: impl Fn(K) -> Option<K>,
        ) {
            *m = std::mem::take(m).into_iter().filter_map(|(k, v)| f(k).map(|k| (k, v))).collect();
        }
        rekey(&mut self.reaction_atom_maps, |(mi, ai)| map(mi).map(|m| (m, ai)));
        rekey(&mut self.enhanced_stereo, |(mi, ai)| map(mi).map(|m| (m, ai)));
        rekey(&mut self.double_bond_align, |(mi, bi)| map(mi).map(|m| (m, bi)));
        rekey(&mut self.reaction_auto_numbers, map);
        rekey(&mut self.ring_fills, map);
        rekey(&mut self.mol_highlights, map);
        self.selected_mols = self.selected_mols.iter().filter_map(|&mi| map(mi)).collect();
        for grp in &mut self.groups {
            *grp = grp.iter().filter_map(|&mi| map(mi)).collect();
        }
        self.groups.retain(|g| g.len() >= 2);
        self.linked_captions.retain_mut(|lc| match map(lc.molecule_idx) {
            Some(m) => {
                lc.molecule_idx = m;
                true
            }
            None => false,
        });
    }

    /// Re-key side tables after molecule `drop` was appended into `keep` (its atoms/bonds offset by
    /// `atom_off`/`bond_off`) and then removed from the scene.
    fn remap_after_mol_merge(&mut self, keep: usize, drop: usize, atom_off: usize, bond_off: usize) {
        fn moveinto<V>(m: &mut HashMap<(usize, usize), V>, keep: usize, drop: usize, off: usize) {
            *m = std::mem::take(m)
                .into_iter()
                .map(|((mi, i), v)| if mi == drop { ((keep, i + off), v) } else { ((mi, i), v) })
                .collect();
        }
        moveinto(&mut self.reaction_atom_maps, keep, drop, atom_off);
        moveinto(&mut self.enhanced_stereo, keep, drop, atom_off);
        moveinto(&mut self.double_bond_align, keep, drop, bond_off);
        for lc in &mut self.linked_captions {
            if lc.molecule_idx == drop {
                lc.molecule_idx = keep;
            }
        }
        // The merged structure keeps its own colours, else adopts the absorbed one's.
        for m in [&mut self.ring_fills, &mut self.mol_highlights] {
            if let Some(c) = m.get(&drop).copied() {
                m.entry(keep).or_insert(c);
            }
        }
        self.remap_after_mol_removal(&[drop]);
    }

    /// Re-key atom-indexed side tables after atom `ai` of molecule `mi` was deleted.
    fn remap_after_atom_removal(&mut self, mi: usize, ai: usize) {
        fn shift<V>(m: &mut HashMap<(usize, usize), V>, mi: usize, ai: usize) {
            *m = std::mem::take(m)
                .into_iter()
                .filter_map(|((m0, a), v)| match (m0 == mi, a.cmp(&ai)) {
                    (true, std::cmp::Ordering::Equal) => None,
                    (true, std::cmp::Ordering::Greater) => Some(((m0, a - 1), v)),
                    _ => Some(((m0, a), v)),
                })
                .collect();
        }
        shift(&mut self.reaction_atom_maps, mi, ai);
        shift(&mut self.enhanced_stereo, mi, ai);
    }

    /// Re-key bond-indexed side tables after bond `bi` of molecule `mi` was deleted.
    fn remap_after_bond_removal(&mut self, mi: usize, bi: usize) {
        self.double_bond_align = std::mem::take(&mut self.double_bond_align)
            .into_iter()
            .filter_map(|((m0, b), v)| match (m0 == mi, b.cmp(&bi)) {
                (true, std::cmp::Ordering::Equal) => None,
                (true, std::cmp::Ordering::Greater) => Some(((m0, b - 1), v)),
                _ => Some(((m0, b), v)),
            })
            .collect();
    }

    fn new_at(doc: InteractiveDocument, offset_idx: usize) -> Self {
        let cascade = ((offset_idx % 8) as f32) * 24.0;
        let title = doc.name.clone();
        Self {
            tab_id: NEXT_TAB_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
            doc,
            file_path: None,
            meta_title: title,
            meta_author: String::new(),
            meta_subject: String::new(),
            meta_comments: String::new(),
            paper_size: "US Letter".to_string(),
            landscape: false,
            page_width_pt: 540.0,
            page_height_pt: 719.75,
            selected_mols: Vec::new(),
            groups: Vec::new(),
            reaction_atom_maps: HashMap::new(),
            reaction_auto_numbers: HashMap::new(),
            enhanced_stereo: HashMap::new(),
            next_or_group: 1,
            next_and_group: 1,
            ring_fills: HashMap::new(),
            mol_highlights: HashMap::new(),
            gel_plates: Vec::new(),
            double_bond_align: HashMap::new(),
            show_atom_numbers: false,
            bracket_usage: "SRU".to_string(),
            bracket_repeat_pattern: "Head-to-Tail".to_string(),
            bracket_repeat_count: "n".to_string(),
            next_rxn_map_id: 1,
            annotations: BTreeMap::new(),
            attached_data: Vec::new(),
            linked_captions: Vec::new(),
            win_pos: Vec2::new(24.0 + cascade, 18.0 + cascade),
            win_size: Vec2::new(760.0, 620.0),
            scroll: Vec2::ZERO,
            maximized: true,
            minimized: false,
            side_history: Vec::new(),
            side_cursor: 0,
        }
    }

    fn effective_page_size_pt(&self) -> (f64, f64) {
        if self.landscape {
            (self.page_height_pt, self.page_width_pt)
        } else {
            (self.page_width_pt, self.page_height_pt)
        }
    }

    fn new(doc: InteractiveDocument) -> Self {
        Self::new_at(doc, 0)
    }
}

struct ChiralCraftApp {
    tabs: Vec<OpenDocTab>,
    active_tab: usize,
    next_untitled_num: usize,
    tool: Tool,
    active_flyout: ActiveFlyout,
    zoom: f32,
    mdi_drag: MdiDragMode,
    sel_handle_drag: Option<usize>,
    sel_drag_prev: Option<(f64, f64)>,
    sel_drag_start_pt: Option<(f64, f64)>,
    sel_drag_start_bbox: Option<(f64, f64, f64, f64)>,
    sel_drag_start_mols: Option<Vec<ChemMolecule>>,
    last_desk_size: Vec2,
    pending_pubchem_caps: HashMap<u32, u64>,
    /// (dialog key, gesture) of the last coalescable dialog edit; see `coalesced_scene_tx`.
    undo_coalesce: Option<(egui::Id, egui::Id)>,
    roundtrip_saved_on_close: bool,
    /// The user chose to quit despite unsaved documents; let the next close request through.
    exit_confirmed: bool,
    clipboard_mols: Vec<ChemMolecule>,
    hovered_atom: Option<(usize, usize)>,
    hovered_bond: Option<(usize, usize)>,
    drag_start_canvas: Option<(f64, f64)>,
    drag_cur_canvas: Option<(f64, f64)>,
    drag_start_atom: Option<(usize, usize)>,
    dragging_single_atom: Option<(usize, usize)>,
    dragging_selected_move: bool,
    dragging_caption: Option<usize>,
    drag_ctrl_copy: bool,
    lasso_pts: Vec<(f64, f64)>,
    eraser_pts: Vec<(f64, f64)>,
    inline_atom_edit: Option<(usize, usize)>,
    inline_label_buf: String,
    inline_just_opened: bool,
    cursor_doc_pos: Option<(f64, f64)>,
    show_tools: bool,
    show_general_toolbar: bool,
    show_style_toolbar: bool,
    show_object_toolbar: bool,
    show_structure_toolbar: bool,
    show_helm_toolbar: bool,
    show_biographics_toolbar: bool,
    show_ring_fill_toolbar: bool,
    show_stoichiometry: bool,
    show_nmr_window: bool,
    nmr_prediction: Option<NmrPrediction>,
    nmr_mol_idx: usize,
    nmr_selected_atom: Option<usize>,
    stoich_base_mmol: f64,
    stoich_rows: Vec<StoichiometryRow>,
    monomer_lib: MonomerLibrary,
    helm_seq_buf: String,
    helm_polymer_kind: PolymerKind,
    name_dict: NameDictionary,
    name_input: String,
    helm_input: String,
    dict_name_input: String,
    show_reaction_interp: bool,
    show_analysis: bool,
    show_properties: bool,
    show_warnings: bool,
    show_stereo: bool,
    show_periodic_table: bool,
    show_info_window: bool,
    show_char_map: bool,
    show_rulers: bool,
    show_crosshair: bool,
    color_by_element: bool,
    analysis_state: AnalysisPanelState,
    props_state: PropsPanelState,
    active_element_z: Option<i32>,
    active_color: Color32,
    text_bold: bool,
    text_italic: bool,
    text_underline: bool,
    text_formula: bool,
    text_sub: bool,
    text_sup: bool,
    text_align: TextAlign,
    line_spacing_pct: f32,
    rlogic_group: String,
    rlogic_occurrence: String,
    rlogic_rest_h: bool,
    insert_obj_kind: usize,
    insert_obj_rows: usize,
    insert_obj_cols: usize,
    insert_obj_text: String,
    curve_filled: bool,
    curve_dashed: bool,
    curve_bold: bool,
    curve_doubled: bool,
    curve_shaded: bool,
    curve_faded: bool,
    curve_arrow_end: bool,
    curve_arrow_start: bool,
    curve_half_start: i8,
    curve_half_end: i8,
    last_repeat_action: Option<String>,
    check_warn_idx: usize,
    dialog: ActiveDialog,
    browser_dir: PathBuf,
    file_path_buf: String,
    smiles_input: String,
    inchi_input: String,
    multiline_input: String,
    scale_pct_input: String,
    rotate_deg_input: String,
    label_input: String,
    nick_search: String,
    nick_def_name: String,
    nick_pending_collision: Option<NicknameCollision>,
    nick_confirm_replace: bool,
    nick_confirm_override: bool,
    data_key_input: String,
    data_val_input: String,
    find_elem_input: String,
    status_msg: String,
    ctp_libraries: Vec<TemplateLibrary>,
    stylesheets: Vec<(String, DocumentSettings)>,
    selected_ctp_idx: usize,
    style_initialized: bool,
    show_chem_db: bool,
    chem_db_state: ChemDbState,
    live_name_resolver: LiveNameResolver,
    /// Active round-trip embedded SVG target path when launched with `--edit-svg <path.svg>`.
    roundtrip_svg_path: Option<String>,
    /// Last editable SVG pushed to OS clipboard so external SVG pastes from WordCraft/DeckCraft take priority.
    last_copied_svg: Option<String>,
    /// Cached About-dialog ChiralCraft emblem texture.
    app_icon_tex: Option<egui::TextureHandle>,
    /// Independent UI & canvas appearance theme (`ClassicLight`, `DarkChrome`, `DarkSlide`).
    ui_theme: UiTheme,
    /// 3D depth cueing opacity for rear atoms/bonds (`0.15..=1.0`, `1.0` = off).
    rear_opacity: f32,
    /// Bridge clearance gap radius in points where front bonds cross over rear bonds.
    crossing_clearance_pt: f64,
    /// Floating Object Inspector window for TLC plates, tables, and 3D depth controls.
    show_object_inspector: bool,
    /// Export & Copy resolution / padding / background preferences.
    export_opts: ExportOptions,
    /// Command Palette (`:` / `F4`) fuzzy action & shortcut launcher.
    cmd_palette: CommandPaletteState,
    /// Interactive 3D Molecular & Crystal Viewer state.
    viewer_3d: Viewer3dState,
    /// Chemistry-Aware Document Diff / Compare window state.
    show_diff_window: bool,
    diff_left_tab: usize,
    diff_right_tab: usize,
    diff_file_input: String,
    last_diff: Option<DocumentDiff>,
    /// Flagship MedChem, HRMS, Isomer/Resonance, Markush SAR, Retrosynthesis & TOC suite state.
    pro_suite: ProSuiteState,
    /// Optional loaded experimental JCAMP-DX spectrum overlaid behind the predicted NMR plot.
    jcamp_overlay: Option<JcampSpectrum>,
    /// All-Discipline Omni-Suite state (FT-IR/UV/2D-NMR, 18e-/Polymer, Energy Profiles, ELN/GHS, LaTeX, Screening).
    omni_suite: OmniSuiteState,
    /// Persistent user preferences (~/.config/chiralcraft/config.json) + Recent Files.
    user_config: UserConfig,
    /// Unified Right-Side Tabbed Studio Dock.
    show_right_dock: bool,
    dock_tab: DockTab,
    /// Computational Chemistry Bridge (ORCA / Gaussian / xTB) state.
    quantum: QuantumStudioState,
    /// Non-blocking local quantum runner (`xtb`, `orca`, `g16`) for real harmonic FT-IR.
    qc_runner: QuantumRunnerState,
    /// Optional online NIST Chemistry WebBook experimental IR (`.jdx`) fetcher.
    nist_ir_fetch: NistIrFetchState,
    /// Optional deterministic GPU screenshot capture path (`chiralcraft --screenshot <preset> <out.png>`).
    screenshot_out_path: Option<String>,
    screenshot_frames: u8,
}

impl ChiralCraftApp {
    fn new() -> Self {
        let doc = InteractiveDocument::new_saved_seed(
            "Untitled ACS Document 1996-1",
            DocumentSettings::acs_document_1996(),
        );
        let ctp_libraries = TemplateLibrary::builtin_libraries();
        let stylesheets = DocumentSettings::builtin_journal_stylesheets();
        let mut monomer_lib = MonomerLibrary::default();
        monomer_lib.populate_builtin_defaults();
        let cwd = env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
        let user_config = if cfg!(test) { UserConfig::default() } else { UserConfig::load() };
        let initial_theme = match user_config.theme.as_str() {
            "dark" => UiTheme::DarkChrome,
            "slide" => UiTheme::DarkSlide,
            _ => UiTheme::ClassicLight,
        };
        let initial_export_opts = user_config.to_export_options();
        let initial_rear_op = user_config.rear_opacity;
        let initial_clearance = user_config.crossing_clearance_pt;
        let initial_dock = user_config.show_right_dock;
        Self {
            tabs: vec![OpenDocTab::new(doc)],
            active_tab: 0,
            next_untitled_num: 2,
            tool: Tool::SolidBond,
            active_flyout: ActiveFlyout::None,
            zoom: 2.0,
            mdi_drag: MdiDragMode::None,
            sel_handle_drag: None,
            sel_drag_prev: None,
            sel_drag_start_pt: None,
            sel_drag_start_bbox: None,
            sel_drag_start_mols: None,
            last_desk_size: Vec2::new(1180.0, 680.0),
            pending_pubchem_caps: HashMap::new(),
            undo_coalesce: None,
            roundtrip_saved_on_close: false,
            exit_confirmed: false,
            clipboard_mols: Vec::new(),
            hovered_atom: None,
            hovered_bond: None,
            drag_start_canvas: None,
            drag_cur_canvas: None,
            drag_start_atom: None,
            dragging_single_atom: None,
            dragging_selected_move: false,
            dragging_caption: None,
            drag_ctrl_copy: false,
            lasso_pts: Vec::new(),
            eraser_pts: Vec::new(),
            inline_atom_edit: None,
            inline_label_buf: String::new(),
            inline_just_opened: false,
            cursor_doc_pos: None,
            show_tools: true,
            show_general_toolbar: true,
            show_style_toolbar: true,
            show_object_toolbar: false,
            show_structure_toolbar: false,
            show_helm_toolbar: false,
            show_biographics_toolbar: false,
            show_ring_fill_toolbar: false,
            show_stoichiometry: false,
            show_nmr_window: false,
            nmr_prediction: None,
            nmr_mol_idx: 0,
            nmr_selected_atom: None,
            stoich_base_mmol: 1.0,
            stoich_rows: Vec::new(),
            monomer_lib,
            helm_seq_buf: String::from("ACDEF"),
            helm_polymer_kind: PolymerKind::Peptide,
            name_dict: NameDictionary::with_builtins(),
            name_input: String::from("4-nitrobenzoic acid"),
            helm_input: String::from("PEPTIDE1{A.C.D.E.F}$$$$V2.0"),
            dict_name_input: String::new(),
            show_reaction_interp: false,
            show_analysis: false,
            show_properties: false,
            show_warnings: true,
            show_stereo: false,
            show_periodic_table: false,
            show_info_window: false,
            show_char_map: false,
            show_rulers: false,
            show_crosshair: false,
            color_by_element: false,
            analysis_state: AnalysisPanelState::default(),
            props_state: PropsPanelState::default(),
            active_element_z: None,
            active_color: Color32::BLACK,
            text_bold: false,
            text_italic: false,
            text_underline: false,
            text_formula: true,
            text_sub: false,
            text_sup: false,
            text_align: TextAlign::Left,
            line_spacing_pct: 100.0,
            rlogic_group: "R1".to_string(),
            rlogic_occurrence: "1".to_string(),
            rlogic_rest_h: false,
            insert_obj_kind: 0,
            insert_obj_rows: 2,
            insert_obj_cols: 2,
            insert_obj_text: "Reaction Conditions".to_string(),
            curve_filled: false,
            curve_dashed: false,
            curve_bold: false,
            curve_doubled: false,
            curve_shaded: false,
            curve_faded: false,
            curve_arrow_end: false,
            curve_arrow_start: false,
            curve_half_start: 0,
            curve_half_end: 0,
            last_repeat_action: None,
            check_warn_idx: 0,
            dialog: ActiveDialog::None,
            browser_dir: cwd,
            file_path_buf: String::from("export.cdxml"),
            smiles_input: String::from("CC(=O)Oc1ccccc1C(=O)O"),
            inchi_input: String::from("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"),
            multiline_input: String::new(),
            scale_pct_input: String::from("200"),
            rotate_deg_input: String::from("45"),
            label_input: String::new(),
            nick_search: String::new(),
            nick_def_name: String::new(),
            nick_pending_collision: None,
            nick_confirm_replace: false,
            nick_confirm_override: false,
            data_key_input: String::from("Tag"),
            data_val_input: String::from("104"),
            find_elem_input: String::from("O"),
            status_msg: "Ready".to_string(),
            ctp_libraries,
            stylesheets,
            selected_ctp_idx: 0,
            style_initialized: false,
            show_chem_db: false,
            chem_db_state: ChemDbState::default(),
            live_name_resolver: LiveNameResolver::default(),
            roundtrip_svg_path: None,
            last_copied_svg: None,
            app_icon_tex: None,
            ui_theme: initial_theme,
            rear_opacity: initial_rear_op,
            crossing_clearance_pt: initial_clearance,
            show_object_inspector: false,
            export_opts: initial_export_opts,
            cmd_palette: CommandPaletteState::default(),
            viewer_3d: Viewer3dState::default(),
            show_diff_window: false,
            diff_left_tab: 0,
            diff_right_tab: 0,
            diff_file_input: String::new(),
            last_diff: None,
            pro_suite: ProSuiteState::default(),
            jcamp_overlay: None,
            omni_suite: OmniSuiteState::default(),
            user_config,
            show_right_dock: initial_dock,
            dock_tab: DockTab::Analysis,
            quantum: QuantumStudioState::default(),
            qc_runner: QuantumRunnerState::default(),
            nist_ir_fetch: NistIrFetchState::default(),
            screenshot_out_path: None,
            screenshot_frames: 0,
        }
    }

    /// Launch ChiralCraft pre-loaded with a deterministic showcase scene and save a native GPU screenshot.
    fn new_showcase_preset(preset: &str, out_path: &str) -> Self {
        let mut app = Self::new();
        app.screenshot_out_path = Some(out_path.to_string());
        app.screenshot_frames = 0;
        app.ui_theme = UiTheme::ClassicLight;
        app.show_right_dock = true;
        app.zoom = 1.45;
        app.tab_mut().maximized = true;
        app.show_stereo = true;
        let blen = app.doc().settings.bond_length_pt;

        let place_smi = |smi: &str, cx: f64, cy: f64, blen: f64| -> ChemMolecule {
            let mut m = smiles_to_molecule(smi).unwrap_or_default();
            clean_up_structure_preserving_stereo(&mut m, &[], blen);
            let (x1, y1, x2, y2) = mol_bbox(&m);
            let dx = cx - 0.5 * (x1 + x2);
            let dy = cy - 0.5 * (y1 + y2);
            for a in &mut m.atoms {
                a.x += dx;
                a.y += dy;
            }
            m
        };

        match preset {
            "medchem" => {
                app.doc_mut().name = "Lead_Optimization_SAR.cdxml".to_string();
                app.show_warnings = false;
                let m = place_smi("CC(C)Cc1ccc(C(C)C(=O)O)cc1", 130.0, 125.0, blen);
                let m2 = place_smi("Brc1ccc(NC(=O)C)cc1", 130.0, 235.0, blen);
                app.doc_mut().stamp_molecule_tx(m, "Ibuprofen");
                app.doc_mut().stamp_molecule_tx(m2, "4-Bromoacetanilide");
                app.tab_mut().selected_mols = vec![1];
                app.dock_tab = DockTab::MedChem;
                app.pro_suite.show_hrms = true;
                app.status_msg = "MedChem ADMET Scorecard + HRMS Isotope Envelope (79Br/81Br doublet)".to_string();
            }
            "spectroscopy" => {
                app.doc_mut().name = "Aspirin_Spectroscopy.cdxml".to_string();
                app.show_warnings = false;
                let m = place_smi("CC(=O)Oc1ccccc1C(=O)O", 110.0, 175.0, blen);
                app.doc_mut().stamp_molecule_tx(m, "Aspirin");
                app.tab_mut().selected_mols = vec![0];
                app.stamp_nmr_shift_labels(false);
                app.dock_tab = DockTab::Spectra;
                app.omni_suite.show_spectroscopy = true;
                app.omni_suite.spec_tab = SpecTab::FtIr;
                app.status_msg = "1H-NMR On-Canvas Shift Badges + FT-IR / UV-Vis / 2D NMR Studio".to_string();
            }
            "3d-quantum" => {
                app.doc_mut().name = "Adamantane_DFT_3D.cdxml".to_string();
                app.show_warnings = false;
                let mut adam = build_adamantane_2d(blen);
                let (x1, y1, x2, y2) = mol_bbox(&adam);
                for a in &mut adam.atoms {
                    a.x += 130.0 - 0.5 * (x1 + x2);
                    a.y += 155.0 - 0.5 * (y1 + y2);
                }
                let lifted = lift_molecule_2d_to_3d(&adam);
                let (opt, rep) = optimize_geometry_3d(&lifted, 120);
                app.doc_mut().stamp_molecule_tx(opt, "Adamantane 3D");
                app.tab_mut().selected_mols = vec![0];
                app.viewer_3d.open = true;
                app.viewer_3d.last_opt_summary = format!(
                    "UFF 3D Opt ({} iters): E {:.1} -> {:.2} kcal/mol",
                    rep.iterations, rep.initial_energy_kcal, rep.final_energy_kcal
                );
                app.dock_tab = DockTab::Quantum;
                app.status_msg = "Interactive 3D Viewer + UFF-Lite Optimizer + ORCA/Gaussian/xTB Deck".to_string();
            }
            "synthesis" => {
                app.doc_mut().name = "Retrosynthesis_EnergyProfile.cdxml".to_string();
                app.show_warnings = false;
                let target = place_smi("CC(=O)Nc1ccccc1", 115.0, 88.0, blen);
                let hits = suggest_disconnections(&target, blen);
                app.doc_mut().stamp_molecule_tx(target.clone(), "Acetanilide");
                app.tab_mut().selected_mols = vec![0];
                if let Some(h0) = hits.first().cloned() {
                    app.pro_suite.retro_hits = hits;
                    app.pro_suite.highlighted_retro_bond = Some(h0.bond_idx);
                    app.doc_mut().apply_scene_tx("Retrosynthesis", |s, nid| {
                        stamp_retrosynthesis_scheme(s, nid, &target, &h0, 65.0, 175.0);
                        s.energy_profiles.push(chiralcraft_doc::EnergyProfile::starter(901, 45.0, 240.0));
                    });
                }
                app.pro_suite.show_retro = true;
                app.dock_tab = DockTab::Analysis;
                app.status_msg = "Retrosynthesis Disconnection Advisor + Reaction Energy Profile (ΔG‡)".to_string();
            }
            "organometallic" => {
                app.doc_mut().name = "Organometallic_Catalysis.cdxml".to_string();
                app.show_warnings = false;
                let cat_lib = app.ctp_libraries.iter().find(|l| l.name == "Organometallic Catalysts").cloned();
                if let Some(lib) = cat_lib {
                    if let Some(pd) = lib.stamp(0, 130.0, 120.0) {
                        app.doc_mut().stamp_molecule_tx(pd, "Pd(PPh3)4");
                    }
                    if let Some(pt) = lib.stamp(1, 130.0, 235.0) {
                        app.doc_mut().stamp_molecule_tx(pt, "Cisplatin");
                    }
                }
                app.tab_mut().selected_mols = vec![0];
                app.omni_suite.show_ino_poly = true;
                app.dock_tab = DockTab::Analysis;
                app.status_msg = "Organometallic 18-Electron Green CBC Analyzer & Catalyst Library".to_string();
            }
            "dark-slide" => {
                app.ui_theme = UiTheme::DarkSlide;
                app.show_warnings = false;
                app.doc_mut().name = "Keynote_Slide_Dark.svg".to_string();
                let mut haw = build_haworth_pyranose(blen);
                let (hx1, hy1, hx2, hy2) = mol_bbox(&haw);
                for a in &mut haw.atoms {
                    a.x += 105.0 - 0.5 * (hx1 + hx2);
                    a.y += 120.0 - 0.5 * (hy1 + hy2);
                }
                let mut adam = build_adamantane_2d(blen);
                let (ax1, ay1, ax2, ay2) = mol_bbox(&adam);
                for a in &mut adam.atoms {
                    a.x += 235.0 - 0.5 * (ax1 + ax2);
                    a.y += 120.0 - 0.5 * (ay1 + ay2);
                }
                app.doc_mut().stamp_molecule_tx(haw, "Haworth Pyranose");
                app.doc_mut().stamp_molecule_tx(adam, "Adamantane");
                app.doc_mut().apply_scene_tx("Dark Slide Cycle", |s, nid| {
                    insert_starter_catalytic_cycle(s, nid, 175.0, 265.0);
                });
                app.cmd_palette.open_fresh();
                app.dock_tab = DockTab::Inspector;
                app.status_msg = "Dark Slide Mode + Crossing-Bond Bridge Gaps + Command Palette (: / F4)".to_string();
            }
            _ => {
                // "hero": Clean publication reaction scheme + Aspirin + Ibuprofen + TLC plate + Energy Profile + Right Studio Dock
                app.doc_mut().name = "Synthesis_Scheme_ACS1996.cdxml".to_string();
                app.show_warnings = false;
                let r1 = place_smi("Oc1ccccc1C(=O)O", 92.0, 105.0, blen);
                let p1 = place_smi("CC(=O)Oc1ccccc1C(=O)O", 275.0, 105.0, blen);
                let ibu = place_smi("CC(C)Cc1ccc(C(C)C(=O)O)cc1", 110.0, 235.0, blen);
                app.doc_mut().stamp_molecule_tx(r1, "Salicylic Acid");
                app.doc_mut().stamp_molecule_tx(p1, "Aspirin");
                app.doc_mut().stamp_molecule_tx(ibu, "Ibuprofen");
                app.tab_mut().selected_mols = vec![1];
                app.doc_mut().apply_scene_tx("Hero Scheme Elements", |s, nid| {
                    let aid = *nid; *nid += 1;
                    s.arrows.push(SceneArrow {
                        id: aid,
                        tail: (145.0, 105.0),
                        head: (218.0, 105.0),
                        arrowhead_head: "Full".to_string(),
                        arrowhead_tail: None,
                        no_go: None,
                        dipole: false,
                    });
                    stamp_arrow_conditions(s, nid, 0, "Ac2O, H3PO4 (cat.)", "85 °C, 15 min, 92%");
                    let tid = *nid; *nid += 1;
                    s.tlc_plates.push(SceneTlcPlate {
                        id: tid,
                        bbox: (415.0, 62.0, 468.0, 148.0),
                        origin_fraction: 0.16,
                        spot_rfs: vec![0.28, 0.54, 0.76],
                        spot_labels: vec!["SM".into(), "Co".into(), "Prod".into()],
                        spot_colors: vec![(20, 20, 20), (25, 95, 205), (205, 45, 45)],
                        solvent_label: "EtOAc/Hex 1:3".into(),
                    });
                    s.energy_profiles.push(chiralcraft_doc::EnergyProfile::starter(905, 215.0, 185.0));
                });
                app.attach_linked_caption(LinkedCaptionKind::FullSummary);
                app.dock_tab = DockTab::Analysis;
                app.status_msg = "ChiralCraft 1.0 — ACS 1996 • Live Linked Captions • TLC • Energy Profile".to_string();
            }
        }
        app.doc_mut().save();
        app
    }

    /// Launch ChiralCraft pre-loaded with a file (or in `--edit-svg` round-trip mode for WordCraft/DeckCraft).
    fn new_with_file(path: &str, roundtrip_svg: bool) -> Self {
        let mut app = Self::new();
        let file_exists_nonempty = fs::metadata(path).map(|m| m.len() > 0).unwrap_or(false);
        let mut bind_roundtrip = false;
        if file_exists_nonempty {
            match load_document_any_with_cell(path) {
                Ok((doc, cell_opt)) => {
                    let mut tab = OpenDocTab::new(doc);
                    tab.file_path = Some(path.to_string());
                    tab.maximized = true;
                    app.tabs = vec![tab];
                    app.active_tab = 0;
                    if let Some(cell) = cell_opt {
                        app.viewer_3d.unit_cell = Some(cell);
                        app.viewer_3d.open = true;
                    }
                    bind_roundtrip = roundtrip_svg;
                    app.status_msg = if roundtrip_svg {
                        format!("Round-trip editing {path} — Ctrl+S or close window to apply edits")
                    } else {
                        format!("Opened {path}")
                    };
                }
                Err(e) => {
                    // Do NOT bind `file_path` or `roundtrip_svg_path` when loading an existing non-empty
                    // file failed — otherwise closing/saving would overwrite the user's file with an empty canvas!
                    app.tabs[0].maximized = true;
                    app.status_msg = format!("Failed to open {path}: {e} (opened blank canvas; original file left untouched)");
                }
            }
        } else {
            app.tabs[0].file_path = Some(path.to_string());
            app.tabs[0].maximized = true;
            bind_roundtrip = roundtrip_svg;
            app.status_msg = format!("New embedded diagram -> {path}");
        }
        if bind_roundtrip {
            app.roundtrip_svg_path = Some(path.to_string());
        }
        app
    }

    /// Atomically flush the round-trip tab's document back to `roundtrip_svg_path` (called on Ctrl+S & window exit).
    fn save_roundtrip_svg(&mut self) -> Result<(), String> {
        let Some(ref path) = self.roundtrip_svg_path.clone() else {
            return Ok(());
        };
        // Write the tab bound to the host SVG, not whichever tab happens to be active
        // (other tabs may have been opened since launch).
        let Some(ti) = self.tabs.iter().position(|t| t.file_path.as_deref() == Some(path.as_str())) else {
            return Ok(());
        };
        let dark_ink = self.ui_theme == UiTheme::DarkSlide;
        let svg = self.tabs[ti].doc.to_svg_styled(true, dark_ink, self.rear_opacity, self.crossing_clearance_pt);
        atomic_write_text(path, &svg)?;
        self.tabs[ti].doc.save();
        self.status_msg = format!("Applied editable SVG updates to {path}");
        Ok(())
    }

    #[inline]
    fn tab(&self) -> &OpenDocTab {
        &self.tabs[self.active_tab]
    }

    #[inline]
    fn tab_mut(&mut self) -> &mut OpenDocTab {
        &mut self.tabs[self.active_tab]
    }

    #[inline]
    fn doc(&self) -> &InteractiveDocument {
        &self.tabs[self.active_tab].doc
    }

    #[inline]
    fn doc_mut(&mut self) -> &mut InteractiveDocument {
        let t = &mut self.tabs[self.active_tab];
        // A step was recorded since the last sync (several edits in one frame): snapshot its side
        // tables now, before this access can record the next one.
        if t.doc.undo_cursor() != t.side_cursor {
            t.sync_side_history();
        }
        &mut t.doc
    }

    /// Indices of molecules targeted by an edit/transform: selected molecules if any, else all.
    fn target_mol_indices(&self) -> Vec<usize> {
        if !self.tab().selected_mols.is_empty() {
            self.tab().selected_mols.clone()
        } else {
            (0..self.doc().scene.molecules.len()).collect()
        }
    }

    /// Primary active molecule (first selected molecule, or first molecule in document).
    fn active_mol_index(&self) -> Option<usize> {
        self.tab()
            .selected_mols
            .first()
            .copied()
            .filter(|&mi| mi < self.doc().scene.molecules.len())
            .or_else(|| if self.doc().scene.molecules.is_empty() { None } else { Some(0) })
    }

    fn active_mol(&self) -> Option<&ChemMolecule> {
        self.active_mol_index().and_then(|mi| self.doc().scene.molecules.get(mi))
    }

    /// Selection-aware `Edit > Copy` — copies structures internally AND pushes Editable SVG (`image/svg+xml`)
    /// onto the OS clipboard for direct pasting into WordCraft / DeckCraft / LibreOffice!
    fn do_copy_action(&mut self) {
        let targets = self.target_mol_indices();
        let copied: Vec<ChemMolecule> = targets
            .iter()
            .filter_map(|&mi| self.doc().scene.molecules.get(mi).cloned())
            .collect();
        if !copied.is_empty() {
            let n = copied.len();
            self.clipboard_mols = copied;
            let (transparent, dark_ink) = match self.export_opts.copy_bg {
                CopyBackground::Transparent => (true, self.ui_theme == UiTheme::DarkSlide),
                CopyBackground::White => (false, false),
                CopyBackground::DarkSlide => (false, true),
            };
            let svg = self.doc().selection_to_svg_styled(
                &self.tab().selected_mols,
                transparent,
                dark_ink,
                self.rear_opacity,
                self.crossing_clearance_pt,
            );
            self.last_copied_svg = Some(svg.clone());
            match copy_svg_to_os_clipboard(&svg) {
                Ok(backend) => {
                    self.status_msg = format!(
                        "Copied {n} structure(s) as Editable SVG ({backend}) — ready to paste into WordCraft/DeckCraft"
                    );
                }
                Err(_) => {
                    self.status_msg = format!("Copied {n} structure(s)");
                }
            }
        }
    }

    /// Explicit `Edit > Copy As > Editable SVG (WordCraft / DeckCraft)`.
    fn do_copy_editable_svg(&mut self, ctx: &egui::Context) {
        let (transparent, dark_ink) = match self.export_opts.copy_bg {
            CopyBackground::Transparent => (true, self.ui_theme == UiTheme::DarkSlide),
            CopyBackground::White => (false, false),
            CopyBackground::DarkSlide => (false, true),
        };
        let svg = self.doc().selection_to_svg_styled(
            &self.tab().selected_mols,
            transparent,
            dark_ink,
            self.rear_opacity,
            self.crossing_clearance_pt,
        );
        self.last_copied_svg = Some(svg.clone());
        let mirror_path = std::env::temp_dir().join("chiralcraft-copy.svg");
        let _ = atomic_write_text(&mirror_path.to_string_lossy(), &svg);
        match copy_svg_to_os_clipboard(&svg) {
            Ok(backend) => {
                self.status_msg = format!(
                    "Copied Editable SVG ({} bytes via {backend}) — paste into WordCraft/DeckCraft",
                    svg.len()
                );
            }
            Err(_) => {
                ctx.copy_text(svg.clone());
                self.status_msg = format!(
                    "Copied Editable SVG XML ({} bytes) to clipboard",
                    svg.len()
                );
            }
        }
    }

    /// Unified `Edit > Paste`: prioritizes newer Editable ChiralCraft SVG from external apps (WordCraft/DeckCraft),
    /// otherwise pastes internal clipboard molecules!
    fn do_paste_action(&mut self) {
        if let Some(svg) = paste_svg_from_os_clipboard() {
            let is_external = self.last_copied_svg.as_deref() != Some(svg.as_str());
            if svg_contains_chiralcraft_document(&svg) && (is_external || self.clipboard_mols.is_empty()) {
                match InteractiveDocument::from_svg(&svg) {
                    Ok(imported) => {
                        let count = imported.scene.molecules.len();
                        let start = self.doc().scene.molecules.len();
                        stamp_molecules_tx(self.doc_mut(), imported.scene.molecules, "Paste Editable SVG");
                        self.tab_mut().selected_mols = (start..start + count).collect();
                        self.status_msg = format!("Pasted {count} editable structure(s) from SVG clipboard");
                        return;
                    }
                    Err(e) => {
                        self.status_msg = format!("SVG paste error: {e}");
                        return;
                    }
                }
            }
        }
        if !self.clipboard_mols.is_empty() {
            let mut clip = self.clipboard_mols.clone();
            let n = clip.len();
            for a in clip.iter_mut().flat_map(|m| m.atoms.iter_mut()) {
                a.x += 20.0;
                a.y += 20.0;
            }
            let start = self.doc().scene.molecules.len();
            stamp_molecules_tx(self.doc_mut(), clip, "Paste");
            self.tab_mut().selected_mols = (start..start + n).collect();
            self.status_msg = format!("Pasted {n} structure(s)");
        }
    }

    /// Selection-aware `Edit > Cut`.
    fn do_cut_action(&mut self) {
        let mut targets = self.target_mol_indices();
        if targets.is_empty() {
            return;
        }
        self.clipboard_mols = targets
            .iter()
            .filter_map(|&mi| self.doc().scene.molecules.get(mi).cloned())
            .collect();
        let n = self.clipboard_mols.len();
        targets.sort_unstable();
        targets.dedup();
        targets.reverse();
        let nmols = self.doc().scene.molecules.len();
        targets.retain(|&mi| mi < nmols);
        let removed = targets.clone();
        self.doc_mut().apply_scene_tx("Cut", |s, _| {
            for mi in targets {
                s.molecules.remove(mi);
            }
        });
        self.tab_mut().remap_after_mol_removal(&removed);
        self.tab_mut().selected_mols.clear();
        self.sanitize_state();
        self.status_msg = format!("Cut {n} structure(s)");
    }

    /// Selection-aware `Edit > Clear` (`Del`): removes selected molecules, or clears whole scene if none selected.
    fn do_clear_action(&mut self) {
        if !self.tab().selected_mols.is_empty() {
            let mut sel = self.tab().selected_mols.clone();
            sel.sort_unstable();
            sel.dedup();
            sel.reverse();
            let nmols = self.doc().scene.molecules.len();
            sel.retain(|&mi| mi < nmols);
            let n = sel.len();
            let removed = sel.clone();
            self.doc_mut().apply_scene_tx("Clear Selection", |s, _| {
                for mi in sel {
                    s.molecules.remove(mi);
                }
            });
            self.tab_mut().remap_after_mol_removal(&removed);
            self.tab_mut().selected_mols.clear();
            self.sanitize_state();
            self.status_msg = format!("Cleared {n} selected structure(s)");
        } else {
            self.doc_mut().apply_scene_tx("Clear", |s, _| {
                s.molecules.clear();
                s.arrows.clear();
                s.curves.clear();
                s.graphics.clear();
                s.aromatic_circles.clear();
                s.captions.clear();
                s.tables.clear();
                s.tlc_plates.clear();
                s.energy_profiles.clear();
            });
            self.tab_mut().gel_plates.clear();
            self.tab_mut().selected_mols.clear();
            self.tab_mut().reaction_atom_maps.clear();
            self.sanitize_state();
            self.status_msg = "Cleared document".to_string();
        }
    }

    /// Real disk-backed `File > Save` (`Ctrl+S`): writes to existing `file_path` (or round-trip SVG) or prompts `Save As`.
    fn do_save_action(&mut self) {
        if self.roundtrip_svg_path.is_some() && self.tab().file_path == self.roundtrip_svg_path {
            if let Err(e) = self.save_roundtrip_svg() {
                self.status_msg = format!("Round-trip save failed: {e}");
            }
            return;
        }
        let current = self.tab().file_path.clone();
        if let Some(path) = current.filter(|p| !save_as_format_is_lossy(self.doc(), p)) {
            let dark_ink = self.ui_theme == UiTheme::DarkSlide;
            match save_document_any_styled(
                self.doc(),
                &path,
                self.export_opts,
                dark_ink,
                self.rear_opacity,
                self.crossing_clearance_pt,
            ) {
                Ok(()) => {
                    self.doc_mut().save();
                    self.status_msg = format!("Saved {path}");
                }
                Err(e) => {
                    self.status_msg = format!("Save failed: {e}");
                }
            }
        } else if let Some(path) = self.tab().file_path.clone() {
            // Current file's format cannot hold the whole drawing: re-saving there and marking the
            // document clean would silently lose data, so ask for a lossless target instead.
            self.file_path_buf = std::path::Path::new(&path).with_extension("cdxml").to_string_lossy().into_owned();
            self.status_msg = format!("{path} cannot hold the whole drawing — choose a CDXML/CDX file to save");
            self.dialog = ActiveDialog::SaveAs;
        } else {
            self.file_path_buf = format!(
                "{}.cdxml",
                self.doc().name.replace(' ', "_")
            );
            self.dialog = ActiveDialog::SaveAs;
        }
    }

    /// Applies an alignment mode to selected molecules (if >=2 selected) or all molecules.
    fn do_align_action(&mut self, label: &str, mode: AlignMode) {
        let sel = if self.tab().selected_mols.len() >= 2 {
            self.tab().selected_mols.clone()
        } else {
            (0..self.doc().scene.molecules.len()).collect()
        };
        if sel.len() < 2 {
            self.status_msg = "Align requires at least 2 structures".to_string();
            return;
        }
        self.doc_mut().apply_scene_tx(label, |s, _| {
            let mut subset: Vec<ChemMolecule> = sel.iter().filter_map(|&mi| s.molecules.get(mi).cloned()).collect();
            align_molecules(&mut subset, mode);
            for (&mi, updated) in sel.iter().zip(subset.into_iter()) {
                if let Some(slot) = s.molecules.get_mut(mi) {
                    *slot = updated;
                }
            }
        });
        self.last_repeat_action = Some(format!("Align:{label}"));
        self.status_msg = format!("Aligned ({label})");
    }

    /// Applies horizontal/vertical distribution to selected molecules (if >=3 selected) or all molecules.
    fn do_distribute_action(&mut self, axis: DistributeAxis) {
        let sel = if self.tab().selected_mols.len() >= 3 {
            self.tab().selected_mols.clone()
        } else {
            (0..self.doc().scene.molecules.len()).collect()
        };
        let label = match axis {
            DistributeAxis::Horizontal => "Distribute Horizontally",
            DistributeAxis::Vertical => "Distribute Vertically",
        };
        if sel.len() < 3 {
            self.status_msg = "Distribute requires at least 3 structures".to_string();
            return;
        }
        self.doc_mut().apply_scene_tx(label, |s, _| {
            let mut subset: Vec<ChemMolecule> = sel.iter().filter_map(|&mi| s.molecules.get(mi).cloned()).collect();
            distribute_molecules(&mut subset, axis);
            for (&mi, updated) in sel.iter().zip(subset.into_iter()) {
                if let Some(slot) = s.molecules.get_mut(mi) {
                    *slot = updated;
                }
            }
        });
        self.last_repeat_action = Some(label.to_string());
        self.status_msg = label.to_string();
    }

    /// Selection-aware `Object > Center on Page` using the active document's real page size.
    fn do_center_on_page_action(&mut self) {
        let (pw, ph) = self.tab().effective_page_size_pt();
        let targets = self.target_mol_indices();
        if targets.is_empty() {
            return;
        }
        self.doc_mut().apply_scene_tx("Center on Page", |s, _| {
            let mut subset: Vec<ChemMolecule> = targets.iter().filter_map(|&mi| s.molecules.get(mi).cloned()).collect();
            center_on_page(&mut subset, pw, ph);
            for (&mi, updated) in targets.iter().zip(subset.into_iter()) {
                if let Some(slot) = s.molecules.get_mut(mi) {
                    *slot = updated;
                }
            }
        });
        self.last_repeat_action = Some("Center on Page".to_string());
        self.status_msg = "Centered on page".to_string();
    }

    /// Selection-aware Flip / Rotate 180°.
    fn do_flip_or_rot180(&mut self, label: &str, axis: FlipAxis, invert_stereo: bool) {
        let targets = self.target_mol_indices();
        if targets.is_empty() {
            return;
        }
        self.doc_mut().apply_scene_tx(label, |s, _| {
            for mi in targets {
                if let Some(m) = s.molecules.get_mut(mi) {
                    let all: Vec<usize> = (0..m.atoms.len()).collect();
                    flip_atoms_preserving_stereo(m, &all, axis, invert_stereo);
                }
            }
        });
        self.last_repeat_action = Some(label.to_string());
        self.status_msg = label.to_string();
    }

    /// Selection-aware `Structure > Clean Up Structure`.
    fn do_clean_up_structure(&mut self) {
        let blen = self.doc().settings.bond_length_pt;
        let targets = self.target_mol_indices();
        if targets.is_empty() {
            return;
        }
        self.doc_mut().apply_scene_tx("Structure Clean Up", |s, _| {
            for mi in targets {
                if let Some(m) = s.molecules.get_mut(mi) {
                    clean_up_structure_preserving_stereo(m, &[], blen);
                }
            }
        });
        self.last_repeat_action = Some("Clean Up Structure".to_string());
        self.status_msg = "Cleaned up structure".to_string();
    }

    /// Selection-aware `Object > Flatten`.
    fn do_flatten_action(&mut self) {
        let targets = self.target_mol_indices();
        if targets.is_empty() {
            return;
        }
        self.doc_mut().apply_scene_tx("Flatten", |s, _| {
            for mi in targets {
                if let Some(m) = s.molecules.get_mut(mi) {
                    flatten_molecule_to_2d(m);
                }
            }
        });
        self.last_repeat_action = Some("Flatten".to_string());
        self.status_msg = "Flattened to 2D".to_string();
    }

    /// Selection-aware `Structure > Toggle Hydrogen Display`.
    fn do_toggle_hydrogens(&mut self) {
        let targets = self.target_mol_indices();
        if targets.is_empty() {
            return;
        }
        self.doc_mut().apply_scene_tx("Toggle Hydrogen Display", |s, _| {
            for mi in targets {
                if let Some(m) = s.molecules.get_mut(mi) {
                    toggle_hydrogen_display(m);
                }
            }
        });
        self.last_repeat_action = Some("Toggle Hydrogen Display".to_string());
        self.status_msg = "Toggled explicit hydrogens".to_string();
    }

    /// Selection-aware `Object > Scale...`.
    fn do_scale_targets(&mut self, factor: f64) {
        let targets = self.target_mol_indices();
        if targets.is_empty() {
            return;
        }
        if let Some((x1, y1, x2, y2)) = self.selection_bbox() {
            let cx = 0.5 * (x1 + x2);
            let cy = 0.5 * (y1 + y2);
            self.doc_mut().apply_scene_tx("Scale", |s, _| {
                scale_selected_molecules(&mut s.molecules, &targets, cx, cy, factor, factor);
            });
            self.last_repeat_action = Some(format!("Scale:{:.4}", factor));
            self.status_msg = format!("Scaled by {:.0}%", factor * 100.0);
        }
    }

    /// Selection-aware `Object > Rotate...`.
    fn do_rotate_targets(&mut self, deg_cw: f64) {
        let targets = self.target_mol_indices();
        if targets.is_empty() {
            return;
        }
        if let Some((x1, y1, x2, y2)) = self.selection_bbox() {
            let cx = 0.5 * (x1 + x2);
            let cy = 0.5 * (y1 + y2);
            let rad = deg_cw.to_radians();
            self.doc_mut().apply_scene_tx("Rotate", |s, _| {
                rotate_selected_molecules(&mut s.molecules, &targets, cx, cy, rad);
            });
            self.last_repeat_action = Some(format!("Rotate:{:.4}", deg_cw));
            self.status_msg = format!("Rotated {deg_cw:.1}° CW");
        }
    }

    // ---- ChiralCraft Helpers (Structure↔Name, NMR, Stoich, Enhanced Stereo, Biopolymers) ----
    /// Compute the live chemical name of `mol`, preferring cached PubChem IUPAC+Common name if available,
    /// otherwise instant pure-Rust IUPAC systematic / trivial naming.
    fn compute_live_name(&self, mol: &ChemMolecule) -> (String, String, bool) {
        let smi = molecule_to_smiles(mol);
        if let Some(res) = self.live_name_resolver.get_cached(&smi) {
            if let Some(disp) = res.format_display() {
                return (disp, smi, false);
            }
        }
        let offline = molecule_to_name(mol, &smi, &self.name_dict);
        let fallback = is_fallback_chemical_name(&offline);
        (offline, smi, fallback)
    }

    /// Stamp a centered chemical name directly beneath the active structure (undoable).
    /// Returns `Some((caption_id, smiles, stamped_name))`.
    fn do_structure_to_name(&mut self) -> Option<(u32, String, String)> {
        let m = self.active_mol().cloned()?;
        let (name, smi, _) = self.compute_live_name(&m);
        let (x1, _, x2, y2) = mol_bbox(&m);
        let cx = 0.5 * (x1 + x2);
        let cy = y2 + 18.0;
        let mut cap_id = 0u32;
        self.doc_mut().apply_scene_tx("Convert Structure to Name", |scene, next_id| {
            cap_id = *next_id;
            *next_id += 1;
            scene.captions.push(CaptionNode {
                id: cap_id,
                text: name.clone(),
                x: cx,
                y: cy,
                element: None,
                generic_nickname: None,
                warning: None,
                justification: LabelJustification::Center,
            });
        });
        self.status_msg = format!("Structure to Name: {name}");
        Some((cap_id, smi, name))
    }

    /// Hybrid `Structure → Name`: stamps the offline IUPAC name immediately, and if `force_pubchem`
    /// is set or the structure is complex (`Compound (...)`), resolves the authoritative PubChem
    /// IUPAC + Common/Drug name in the background and upgrades the caption in place!
    fn do_structure_to_name_hybrid(&mut self, ctx: &egui::Context, force_pubchem: bool) {
        let Some((cap_id, smi, name)) = self.do_structure_to_name() else {
            return;
        };
        if force_pubchem || is_fallback_chemical_name(&name) {
            let tab_id = self.tab().tab_id;
            if let Some(cached) = self.live_name_resolver.resolve(&smi, Some(cap_id), ctx) {
                for cap in &mut self.doc_mut().scene.captions {
                    if cap.id == cap_id {
                        cap.text = cached.clone();
                    }
                }
                self.status_msg = format!("PubChem IUPAC Name: {cached}");
            } else {
                self.pending_pubchem_caps.insert(cap_id, tab_id);
                self.status_msg = format!("Structure to Name: {name} (resolving via PubChem...)");
            }
        }
    }

    fn do_name_to_structure(&mut self, name: &str) -> Result<(), String> {
        let smi = self.name_dict.name_to_smiles(name)?;
        let mut m = smiles_to_molecule(&smi)?;
        let blen = self.doc().settings.bond_length_pt;
        clean_up_structure_preserving_stereo(&mut m, &[], blen);
        for a in &mut m.atoms {
            a.x += 200.0;
            a.y += 170.0;
        }
        self.doc_mut().stamp_molecule_tx(m, "Convert Name to Structure");
        self.status_msg = format!("Converted \"{}\" → structure", name.trim());
        Ok(())
    }

    fn do_predict_nmr(&mut self, is_13c: bool) {
        let Some(mi) = self.active_mol_index() else {
            return;
        };
        let mol = self.doc().scene.molecules[mi].clone();
        let pred = if is_13c { predict_13c_nmr(&mol) } else { predict_1h_nmr(&mol) };
        self.status_msg = pred.report_text.clone();
        self.nmr_prediction = Some(pred);
        self.nmr_mol_idx = mi;
        self.nmr_selected_atom = None;
        self.show_nmr_window = true;
    }

    fn do_autonumber_reaction(&mut self) {
        let split_x = self
            .doc()
            .scene
            .arrows
            .first()
            .map(|a| 0.5 * (a.tail.0 + a.head.0))
            .unwrap_or(f64::INFINITY);
        let mut nums = HashMap::new();
        let mut idx = 1usize;
        let mut order: Vec<(usize, f64)> = self
            .doc()
            .scene
            .molecules
            .iter()
            .enumerate()
            .filter(|(_, m)| !m.atoms.is_empty())
            .map(|(mi, m)| {
                let cx = m.atoms.iter().map(|a| a.x).sum::<f64>() / (m.atoms.len() as f64);
                (mi, cx)
            })
            .collect();
        order.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        for (mi, cx) in order {
            let tag = if split_x.is_finite() && cx > split_x {
                format!("{idx} (Product)")
            } else {
                format!("{idx}")
            };
            nums.insert(mi, tag);
            idx += 1;
        }
        let n = nums.len();
        self.tab_mut().reaction_auto_numbers = nums;
        self.status_msg = format!("Autonumbered {n} compound(s) in reaction");
    }

    fn do_mark_enhanced_stereo(&mut self, kind: EnhancedStereoKind) {
        let targets = self.target_mol_indices();
        let mut marked = 0usize;
        for mi in targets {
            if let Some(mol) = self.doc().scene.molecules.get(mi).cloned() {
                let st = perceive_stereochemistry(&mol);
                for &ai in st.atom_stereo.keys() {
                    self.tab_mut().enhanced_stereo.insert((mi, ai), kind);
                    marked += 1;
                }
                if st.atom_stereo.is_empty() && !mol.atoms.is_empty() {
                    self.tab_mut().enhanced_stereo.insert((mi, 0), kind);
                    marked += 1;
                }
            }
        }
        self.show_stereo = true;
        self.status_msg = format!("Marked {marked} center(s) as {}", kind.display_label());
    }

    fn do_open_stoichiometry_grid(&mut self) {
        self.stoich_rows = build_stoichiometry_rows(&self.doc().scene, self.stoich_base_mmol);
        self.show_stoichiometry = true;
        self.status_msg = format!("Stoichiometry Grid updated ({} row(s))", self.stoich_rows.len());
    }

    /// Compute a clean, non-overlapping canvas center `(x, y)` for newly imported DB structures.
    fn next_import_center(&self) -> (f64, f64) {
        let mols = &self.doc().scene.molecules;
        if mols.is_empty() {
            return (200.0, 160.0);
        }
        let mut max_x = 0.0_f64;
        let mut max_y = 0.0_f64;
        let mut sum_y = 0.0_f64;
        let mut cnt = 0usize;
        for m in mols {
            if !m.atoms.is_empty() {
                let (_, y1, x2, y2) = mol_bbox(m);
                max_x = max_x.max(x2);
                max_y = max_y.max(y2);
                sum_y += 0.5 * (y1 + y2);
                cnt += 1;
            }
        }
        let (pw, _ph) = self.tab().effective_page_size_pt();
        let next_x = max_x + 95.0;
        if next_x < pw - 90.0 {
            (next_x, if cnt > 0 { sum_y / (cnt as f64) } else { 160.0 })
        } else {
            (180.0, max_y + 100.0)
        }
    }

    /// Import a single PubChem/ChEMBL hit onto the active canvas (scaled to active stylesheet).
    fn do_import_db_hit(&mut self, hit: &DbHit, with_caption: bool) {
        let blen = self.doc().settings.bond_length_pt;
        let (cx, cy) = self.next_import_center();
        match hit_to_molecule(hit, blen, cx, cy) {
            Ok(mol) => {
                let (x1, _, x2, y2) = mol_bbox(&mol);
                let label = format!("Import {} ({})", hit.name, hit.id);
                let caption_text = if with_caption { Some(hit.format_caption()) } else { None };
                let line_step = 12.0 * (self.line_spacing_pct as f64 / 100.0).clamp(0.6, 2.5);
                let mut new_mi = 0usize;
                self.doc_mut().apply_scene_tx(&label, |scene, next_id| {
                    let mut m = mol;
                    for a in &mut m.atoms {
                        a.id = *next_id;
                        *next_id += 1;
                    }
                    for b in &mut m.bonds {
                        b.id = *next_id;
                        *next_id += 1;
                    }
                    new_mi = scene.molecules.len();
                    scene.molecules.push(m);
                    if let Some(ref txt) = caption_text {
                        let cap_x = 0.5 * (x1 + x2);
                        let cap_y = y2 + 18.0;
                        for (i, line) in txt.lines().enumerate() {
                            let cid = *next_id;
                            *next_id += 1;
                            scene.captions.push(CaptionNode {
                                id: cid,
                                text: line.to_string(),
                                x: cap_x,
                                y: cap_y + (i as f64) * line_step,
                                element: None,
                                generic_nickname: None,
                                warning: None,
                                justification: LabelJustification::Center,
                            });
                        }
                    }
                });
                self.tab_mut().selected_mols = vec![new_mi];
                self.status_msg = format!(
                    "Imported {} ({}) from {}",
                    hit.name,
                    hit.id,
                    hit.source.label()
                );
            }
            Err(e) => {
                self.status_msg = format!("Import failed for {}: {e}", hit.id);
            }
        }
    }

    /// Import multiple PubChem/ChEMBL hits neatly arranged in a 3-column grid with captions.
    fn do_import_db_batch(&mut self, hits: &[DbHit]) {
        if hits.is_empty() {
            return;
        }
        let blen = self.doc().settings.bond_length_pt;
        let base_y = self
            .selection_bbox()
            .map(|(_, _, _, y2)| y2 + 95.0)
            .unwrap_or(130.0);
        let col_spacing = 175.0_f64;
        let row_spacing = 155.0_f64;
        let mut prepared: Vec<(ChemMolecule, String)> = Vec::new();
        for (idx, h) in hits.iter().enumerate() {
            let col = (idx % 3) as f64;
            let row = (idx / 3) as f64;
            let cx = 120.0 + col * col_spacing;
            let cy = base_y + row * row_spacing;
            if let Ok(m) = hit_to_molecule(h, blen, cx, cy) {
                let short: String = h.name.chars().take(28).collect();
                prepared.push((m, format!("{short}\n{}", h.id)));
            }
        }
        if prepared.is_empty() {
            self.status_msg = "No importable 2D structures in batch".to_string();
            return;
        }
        let count = prepared.len();
        let mut new_indices = Vec::with_capacity(count);
        self.doc_mut().apply_scene_tx("Import Chemical DB Grid", |scene, next_id| {
            for (mut m, cap) in prepared {
                let (x1, _, x2, y2) = mol_bbox(&m);
                for a in &mut m.atoms {
                    a.id = *next_id;
                    *next_id += 1;
                }
                for b in &mut m.bonds {
                    b.id = *next_id;
                    *next_id += 1;
                }
                new_indices.push(scene.molecules.len());
                scene.molecules.push(m);
                let cap_x = 0.5 * (x1 + x2);
                for (li, line) in cap.lines().enumerate() {
                    let cid = *next_id;
                    *next_id += 1;
                    scene.captions.push(CaptionNode {
                        id: cid,
                        text: line.to_string(),
                        x: cap_x,
                        y: y2 + 16.0 + (li as f64) * 11.5,
                        element: None,
                        generic_nickname: None,
                        warning: None,
                        justification: LabelJustification::Center,
                    });
                }
            }
        });
        self.tab_mut().selected_mols = new_indices;
        self.status_msg = format!("Imported {count} structure(s) from database grid");
    }

    /// Take the active/selected canvas structure's SMILES and search PubChem / ChEMBL immediately.
    fn do_search_canvas_in_db(&mut self, ctx: &egui::Context, mode: DbSearchMode) {
        if let Some(m) = self.active_mol() {
            let smi = molecule_to_smiles(m);
            if smi.trim().is_empty() {
                self.status_msg = "Could not derive SMILES from active structure".to_string();
                return;
            }
            self.chem_db_state.query = smi;
            self.chem_db_state.mode = mode;
            self.show_chem_db = true;
            self.chem_db_state.start_search(ctx);
        } else {
            self.status_msg = "Draw or select a structure on the canvas first".to_string();
        }
    }

    fn handle_pro_panel_action(&mut self, ctx: &egui::Context, act: ProPanelAction) {
        match act {
            ProPanelAction::None => {}
            ProPanelAction::PasteText(txt) => self.paste_caption_tx(txt),
            ProPanelAction::HighlightAtom(ai) => {
                self.nmr_selected_atom = ai;
            }
            ProPanelAction::ClearNmrAssignment => {
                self.nmr_selected_atom = None;
                self.tab_mut().annotations.remove("SpectrumAssignment");
                self.status_msg = "Cleared spectrum-structure assignment".to_string();
            }
            ProPanelAction::StampNmrShifts(carbon13) => {
                self.stamp_nmr_shift_labels(carbon13);
            }
            ProPanelAction::ClearStampedNmrShifts => {
                self.clear_nmr_shift_labels();
            }
            ProPanelAction::LoadJcampOverlay(path) => {
                match fs::read_to_string(&path).map_err(|e| e.to_string()).and_then(|s| parse_jcamp_dx(&s)) {
                    Ok(spec) => {
                        self.status_msg = format!("Loaded experimental JCAMP-DX overlay: {} ({} pts)", spec.title, spec.points.len());
                        self.jcamp_overlay = Some(spec);
                    }
                    Err(e) => {
                        self.status_msg = format!("JCAMP-DX load error: {e}");
                    }
                }
            }
            ProPanelAction::ClearJcampOverlay => {
                self.jcamp_overlay = None;
                self.status_msg = "Cleared experimental JCAMP-DX NMR overlay".to_string();
            }
            ProPanelAction::InsertBiopolymer { seq, kind } => {
                let blen = self.doc().settings.bond_length_pt;
                match sequence_to_biopolymer_molecule(&seq, kind, blen) {
                    Ok(m) => {
                        let cnt = m.atoms.len();
                        self.doc_mut().stamp_molecule_tx(m, "Insert Biopolymer");
                        self.status_msg = format!("Built {cnt}-monomer {:?} sequence", kind);
                    }
                    Err(e) => self.status_msg = e,
                }
            }
            ProPanelAction::InsertMonomerAtom { polymer_type, symbol } => {
                if let Some(mut m) = monomer_to_all_atom_molecule(&self.monomer_lib, &polymer_type, &symbol) {
                    let blen = self.doc().settings.bond_length_pt;
                    clean_up_structure_preserving_stereo(&mut m, &[], blen);
                    for a in &mut m.atoms { a.x += 190.0; a.y += 160.0; }
                    self.doc_mut().stamp_molecule_tx(m, "Stamp Monomer");
                    self.status_msg = format!("Stamped all-atom monomer {symbol}");
                }
            }
            ProPanelAction::CleanUpBiopolymer => {
                let blen = self.doc().settings.bond_length_pt;
                let targets = self.target_mol_indices();
                self.doc_mut().apply_scene_tx("Clean Up Biopolymer", |s, _| {
                    for mi in targets {
                        if let Some(m) = s.molecules.get_mut(mi) {
                            clean_up_biopolymer(m, blen);
                        }
                    }
                });
                self.status_msg = "Cleaned up biopolymer layout".to_string();
            }
            ProPanelAction::CopyHelm => {
                if let Some(m) = self.active_mol() {
                    let h = molecule_to_helm(m);
                    ctx.copy_text(h.clone());
                    self.status_msg = format!("Copied HELM: {h}");
                }
            }
            ProPanelAction::CopyFasta => {
                if let Some(m) = self.active_mol() {
                    let f = molecule_to_fasta(m, &self.doc().name);
                    ctx.copy_text(f.clone());
                    self.status_msg = "Copied FASTA sequence".to_string();
                }
            }
            ProPanelAction::StructureCmd(cmd) => match cmd {
                "check" => { self.check_warn_idx = 0; self.dialog = ActiveDialog::CheckStructureDialog; }
                "cleanup" => self.do_clean_up_structure(),
                "clean_rxn" => { self.doc_mut().apply_scene_tx("Clean Up Reaction", |s, _| { clean_up_reaction(s); }); }
                "struct_to_name" => self.do_structure_to_name_hybrid(ctx, false),
                "name_to_struct" => self.dialog = ActiveDialog::PasteNameAsStructure,
                "nmr_1h" => self.do_predict_nmr(false),
                "nmr_13c" => self.do_predict_nmr(true),
                "map_rxn" => {
                    let maps = auto_map_reaction_atoms(&self.doc().scene);
                    self.tab_mut().reaction_atom_maps = maps;
                }
                "autonum" => self.do_autonumber_reaction(),
                "stoich" => self.do_open_stoichiometry_grid(),
                "stereo_abs" => self.do_mark_enhanced_stereo(EnhancedStereoKind::Absolute),
                "stereo_or" => {
                    let n = self.tab().next_or_group;
                    self.tab_mut().next_or_group += 1;
                    self.do_mark_enhanced_stereo(EnhancedStereoKind::Or(n));
                }
                "stereo_and" => {
                    let n = self.tab().next_and_group;
                    self.tab_mut().next_and_group += 1;
                    self.do_mark_enhanced_stereo(EnhancedStereoKind::And(n));
                }
                _ => {}
            },
            ProPanelAction::OpenBioGraphicsLibrary(name) => {
                if let Some(idx) = self.ctp_libraries.iter().position(|l| l.name.eq_ignore_ascii_case(name)) {
                    self.selected_ctp_idx = idx;
                    self.active_flyout = ActiveFlyout::Templates;
                    self.status_msg = format!("Opened BioGraphics library: {name}");
                }
            }
            ProPanelAction::InsertGelPlate => {
                let id = (self.tab().gel_plates.len() + 1) as u32;
                self.tab_mut().gel_plates.push(SceneGelPlate::default_at(id, 140.0, 120.0));
                self.status_msg = "Inserted Gel Electrophoresis Plate (GEP)".to_string();
            }
            ProPanelAction::SetRingFill(col_opt) => {
                let targets = self.target_mol_indices();
                for mi in targets {
                    if let Some(c) = col_opt {
                        self.tab_mut().ring_fills.insert(mi, c);
                    } else {
                        self.tab_mut().ring_fills.remove(&mi);
                    }
                }
                self.status_msg = if col_opt.is_some() { "Applied Ring Fill Color".to_string() } else { "Cleared Ring Fill Color".to_string() };
            }
            ProPanelAction::SetHighlightColor(col_opt) => {
                let targets = self.target_mol_indices();
                for mi in targets {
                    if let Some(c) = col_opt {
                        self.tab_mut().mol_highlights.insert(mi, c);
                    } else {
                        self.tab_mut().mol_highlights.remove(&mi);
                    }
                }
                self.status_msg = if col_opt.is_some() { "Applied Highlight Color".to_string() } else { "Cleared Highlight Color".to_string() };
            }
            ProPanelAction::ExportStoichiometryCsv(csv) => {
                ctx.copy_text(csv.clone());
                match native_dialogs::save_text_with_dialog(
                    "Export Stoichiometry CSV",
                    "chiralcraft-stoichiometry.csv",
                    &csv,
                ) {
                    Ok(Some(path)) => {
                        self.status_msg = format!("Exported Stoichiometry CSV to {path} & clipboard");
                    }
                    Ok(None) => {
                        self.status_msg = "Copied Stoichiometry CSV to clipboard".to_string();
                    }
                    Err(e) => {
                        self.status_msg = format!("CSV export error: {e} (copied to clipboard)");
                    }
                }
            }
        }
    }

    fn reset_transient_interaction_state(&mut self) {
        self.mdi_drag = MdiDragMode::None;
        self.sel_handle_drag = None;
        self.sel_drag_prev = None;
        self.sel_drag_start_pt = None;
        self.sel_drag_start_bbox = None;
        self.sel_drag_start_mols = None;
        self.inline_atom_edit = None;
        self.hovered_atom = None;
        self.hovered_bond = None;
        self.drag_start_canvas = None;
        self.drag_cur_canvas = None;
        self.drag_start_atom = None;
        self.dragging_single_atom = None;
        self.dragging_selected_move = false;
        self.dragging_caption = None;
        self.drag_ctrl_copy = false;
        self.lasso_pts.clear();
        self.eraser_pts.clear();
        self.nmr_selected_atom = None;
        self.nmr_prediction = None;
        self.viewer_3d.vib_anim = None;
        self.viewer_3d.measure_atoms.clear();
        self.pro_suite.highlighted_retro_bond = None;
        self.omni_suite.highlighted_spec_atoms.clear();
    }

    fn new_document_tab(&mut self, settings: DocumentSettings) {
        let name = format!("Untitled {}-{}", settings.stylesheet_name, self.next_untitled_num);
        self.next_untitled_num += 1;
        let d = InteractiveDocument::new_saved_seed(&name, settings);
        let idx = self.tabs.len();
        self.tabs.push(OpenDocTab::new_at(d, idx));
        self.active_tab = self.tabs.len() - 1;
        self.reset_transient_interaction_state();
        self.status_msg = format!("Created {name}");
    }

    /// Close the active tab, asking first when it has unsaved changes.
    fn close_active_tab(&mut self) {
        if self.doc().is_dirty() {
            self.dialog = ActiveDialog::ConfirmCloseTab;
        } else {
            self.force_close_active_tab();
        }
    }

    /// Close the active tab, discarding any unsaved changes.
    fn force_close_active_tab(&mut self) {
        let closing_path = self.tabs.get(self.active_tab).and_then(|t| t.file_path.clone());
        if closing_path.is_some() && closing_path == self.roundtrip_svg_path {
            self.roundtrip_svg_path = None;
        }
        if self.tabs.len() > 1 {
            self.tabs.remove(self.active_tab);
            if self.active_tab >= self.tabs.len() {
                self.active_tab = self.tabs.len() - 1;
            }
        } else {
            let d = InteractiveDocument::new_saved_seed(
                "Untitled Document-1",
                DocumentSettings::acs_document_1996(),
            );
            self.tabs[0] = OpenDocTab::new_at(d, 0);
            self.roundtrip_svg_path = None;
        }
        self.dialog = ActiveDialog::None;
        self.reset_transient_interaction_state();
    }

    /// Guarantees all cached indices (`active_tab`, `hovered_atom`, `hovered_bond`,
    /// `inline_atom_edit`, `selected_mols`, `groups`, `reaction_atom_maps`, side tables) are strictly in-bounds.
    fn sanitize_state(&mut self) {
        if self.tabs.is_empty() {
            let d = InteractiveDocument::new_saved_seed(
                "Untitled Document-1",
                DocumentSettings::acs_document_1996(),
            );
            self.tabs.push(OpenDocTab::new_at(d, 0));
        }
        if self.active_tab >= self.tabs.len() {
            self.active_tab = self.tabs.len() - 1;
        }
        let mols = &self.tabs[self.active_tab].doc.scene.molecules;
        let valid_atom = |pair: (usize, usize)| {
            mols.get(pair.0).map(|m| pair.1 < m.atoms.len()).unwrap_or(false)
        };
        let valid_bond = |pair: (usize, usize)| {
            mols.get(pair.0).map(|m| pair.1 < m.bonds.len()).unwrap_or(false)
        };
        if let Some(p) = self.hovered_atom {
            if !valid_atom(p) { self.hovered_atom = None; }
        }
        if let Some(p) = self.hovered_bond {
            if !valid_bond(p) { self.hovered_bond = None; }
        }
        if let Some(p) = self.drag_start_atom {
            if !valid_atom(p) { self.drag_start_atom = None; }
        }
        if let Some(p) = self.dragging_single_atom {
            if !valid_atom(p) { self.dragging_single_atom = None; }
        }
        if let Some(p) = self.inline_atom_edit {
            if !valid_atom(p) { self.inline_atom_edit = None; }
        }
        let nmols = mols.len();
        let tab = &mut self.tabs[self.active_tab];
        tab.selected_mols.retain(|&mi| mi < nmols);
        for grp in &mut tab.groups {
            grp.retain(|&mi| mi < nmols);
        }
        tab.groups.retain(|g| g.len() >= 2);
        let atom_lens: Vec<usize> = tab.doc.scene.molecules.iter().map(|m| m.atoms.len()).collect();
        let bond_lens: Vec<usize> = tab.doc.scene.molecules.iter().map(|m| m.bonds.len()).collect();
        tab.reaction_atom_maps.retain(|&(mi, ai), _| atom_lens.get(mi).map(|&n| ai < n).unwrap_or(false));
        tab.enhanced_stereo.retain(|&(mi, ai), _| atom_lens.get(mi).map(|&n| ai < n).unwrap_or(false));
        tab.double_bond_align.retain(|&(mi, bi), _| bond_lens.get(mi).map(|&n| bi < n).unwrap_or(false));
        tab.reaction_auto_numbers.retain(|&mi, _| mi < nmols);
        tab.ring_fills.retain(|&mi, _| mi < nmols);
        tab.mol_highlights.retain(|&mi, _| mi < nmols);
        let cap_ids: std::collections::HashSet<u32> = tab.doc.scene.captions.iter().map(|c| c.id).collect();
        tab.linked_captions.retain(|lc| lc.molecule_idx < nmols && cap_ids.contains(&lc.caption_id));
    }

    fn apply_theme(ctx: &egui::Context, theme: UiTheme) {
        match theme {
            UiTheme::ClassicLight => {
                let mut visuals = egui::Visuals::light();
                visuals.panel_fill = CD_WINDOW_BG;
                visuals.window_fill = CD_WINDOW_BG;
                visuals.extreme_bg_color = Color32::WHITE;
                visuals.widgets.noninteractive.bg_fill = CD_WINDOW_BG;
                visuals.widgets.inactive.bg_fill = Color32::from_rgb(245, 245, 245);
                visuals.widgets.hovered.bg_fill = Color32::from_rgb(224, 238, 252);
                visuals.widgets.active.bg_fill = Color32::from_rgb(194, 224, 250);
                visuals.selection.bg_fill = Color32::from_rgb(0, 120, 215);
                visuals.window_stroke = Stroke::new(1.0_f32, Color32::from_rgb(120, 148, 180));
                ctx.set_visuals(visuals);
            }
            UiTheme::DarkChrome | UiTheme::DarkSlide => {
                let mut visuals = egui::Visuals::dark();
                visuals.panel_fill = Color32::from_rgb(34, 38, 46);
                visuals.window_fill = Color32::from_rgb(38, 42, 50);
                visuals.extreme_bg_color = Color32::from_rgb(24, 27, 33);
                visuals.widgets.inactive.bg_fill = Color32::from_rgb(50, 56, 66);
                visuals.widgets.hovered.bg_fill = Color32::from_rgb(66, 86, 114);
                visuals.widgets.active.bg_fill = Color32::from_rgb(55, 110, 175);
                visuals.selection.bg_fill = Color32::from_rgb(48, 128, 218);
                visuals.window_stroke = Stroke::new(1.0_f32, Color32::from_rgb(80, 96, 120));
                ctx.set_visuals(visuals);
            }
        }
    }

    fn effective_ink_color(&self) -> Color32 {
        if self.ui_theme == UiTheme::DarkSlide && self.active_color == Color32::BLACK {
            Color32::from_rgb(242, 245, 249)
        } else {
            self.active_color
        }
    }

    fn canvas_bg_color(&self) -> Color32 {
        if self.ui_theme == UiTheme::DarkSlide {
            Color32::from_rgb(30, 34, 42)
        } else {
            Color32::WHITE
        }
    }

    /// Attaches a live-updating property or IUPAC-name caption directly below the active molecule.
    fn attach_linked_caption(&mut self, kind: LinkedCaptionKind) {
        let Some(mi) = self.active_mol_index() else {
            self.status_msg = "Draw or select a molecule first to attach a live caption".to_string();
            return;
        };
        let Some(mol) = self.doc().scene.molecules.get(mi).cloned() else {
            return;
        };
        let (x1, _, x2, y2) = mol_bbox(&mol);
        let cx = 0.5 * (x1 + x2);
        let existing_for_mol = self
            .tab()
            .linked_captions
            .iter()
            .filter(|lc| lc.molecule_idx == mi)
            .count();
        let cy = y2 + 16.0 + (existing_for_mol as f64) * 14.0;
        let (live_name, _, _) = self.compute_live_name(&mol);
        let decimals = self.analysis_state.decimals;
        let text = format_linked_caption_text(&mol, kind, decimals, &live_name);
        let mut new_cap_id = 0u32;
        self.doc_mut().apply_scene_tx("Attach Live Linked Caption", |scene, next_id| {
            let cid = *next_id;
            *next_id += 1;
            new_cap_id = cid;
            scene.captions.push(CaptionNode {
                id: cid,
                text,
                x: cx,
                y: cy,
                element: None,
                generic_nickname: None,
                warning: None,
                justification: LabelJustification::Center,
            });
        });
        self.tab_mut().linked_captions.push(LinkedCaption {
            caption_id: new_cap_id,
            molecule_idx: mi,
            kind,
            decimals,
        });
        self.status_msg = "Attached live-linked caption (updates automatically on edit)".to_string();
    }

    /// Keeps every attached `LinkedCaption` synchronized live with its molecule's formula, MW, mass, or name.
    fn sync_linked_captions(&mut self) {
        if self.tab().linked_captions.is_empty() {
            return;
        }
        let nmols = self.doc().scene.molecules.len();
        self.tab_mut().linked_captions.retain(|lc| lc.molecule_idx < nmols);
        let links = self.tab().linked_captions.clone();
        let mut updates: Vec<(u32, String)> = Vec::new();
        for lc in &links {
            if let Some(mol) = self.doc().scene.molecules.get(lc.molecule_idx) {
                let (live_name, _, _) = self.compute_live_name(mol);
                let txt = format_linked_caption_text(mol, lc.kind, lc.decimals, &live_name);
                updates.push((lc.caption_id, txt));
            }
        }
        let caps = &mut self.doc_mut().scene.captions;
        for (cid, new_txt) in updates {
            if let Some(c) = caps.iter_mut().find(|c| c.id == cid) {
                if c.text != new_txt {
                    c.text = new_txt;
                }
            }
        }
    }

    fn insert_special_core(&mut self, label: &str, mut mol: ChemMolecule) {
        let (x1, y1, x2, y2) = mol_bbox(&mol);
        let dx = 220.0 - 0.5 * (x1 + x2);
        let dy = 185.0 - 0.5 * (y1 + y2);
        for a in &mut mol.atoms {
            a.x += dx;
            a.y += dy;
        }
        self.doc_mut().stamp_molecule_tx(mol, label);
        self.status_msg = format!("Inserted {label}");
    }

    fn arrange_selected_catalytic_cycle(&mut self, elliptical: bool) {
        let idxs = self.target_mol_indices();
        if idxs.len() < 2 {
            self.status_msg = "Need at least 2 molecules to arrange a catalytic cycle (or insert a starter cycle)".to_string();
            return;
        }
        let (rx, ry) = if elliptical { (130.0, 88.0) } else { (105.0, 105.0) };
        let mut arranged = 0usize;
        self.doc_mut().apply_scene_tx("Arrange Catalytic Cycle", |scene, next_id| {
            arranged = arrange_catalytic_cycle(scene, next_id, &idxs, 250.0, 210.0, rx, ry, true);
        });
        self.status_msg = format!("Arranged {arranged} intermediates into a catalytic cycle");
    }

    /// Stamps predicted ¹H or ¹³C NMR δ (ppm) shift captions directly beside assigned atoms on canvas,
    /// placing each badge radially outward from the ring/molecule centroid and resolving collisions.
    fn stamp_nmr_shift_labels(&mut self, carbon13: bool) {
        let Some(mi) = self.active_mol_index() else {
            self.status_msg = "Draw or select a molecule first to stamp NMR shifts".to_string();
            return;
        };
        let Some(mol) = self.doc().scene.molecules.get(mi).cloned() else {
            return;
        };
        let pred = if carbon13 { predict_13c_nmr(&mol) } else { predict_1h_nmr(&mol) };
        if pred.peaks.is_empty() {
            self.status_msg = format!("No {} NMR peaks to stamp", pred.nucleus);
            return;
        }
        let nf = mol.atoms.len().max(1) as f64;
        let mx = mol.atoms.iter().map(|a| a.x).sum::<f64>() / nf;
        let my = mol.atoms.iter().map(|a| a.y).sum::<f64>() / nf;

        struct PlacedBadge {
            text: String,
            cx: f64,
            cy: f64,
            w: f64,
            h: f64,
            just: LabelJustification,
        }
        let mut badges: Vec<PlacedBadge> = Vec::new();
        let mut seen_atoms = std::collections::HashSet::new();
        for pk in &pred.peaks {
            if !seen_atoms.insert(pk.atom_idx) {
                continue;
            }
            let Some(a) = mol.atoms.get(pk.atom_idx) else {
                continue;
            };
            // Compute outward unit vector away from bonded neighbors and molecule centroid
            let mut vx = 0.0_f64;
            let mut vy = 0.0_f64;
            for (nb, _) in mol.neighbors(pk.atom_idx) {
                if let Some(na) = mol.atoms.get(nb) {
                    let dx = a.x - na.x;
                    let dy = a.y - na.y;
                    let d = (dx * dx + dy * dy).sqrt().max(1e-3);
                    vx += dx / d;
                    vy += dy / d;
                }
            }
            // Blend in centroid vector so internal/ring atoms also point outward
            let cdx = a.x - mx;
            let cdy = a.y - my;
            let cd = (cdx * cdx + cdy * cdy).sqrt();
            if cd > 1e-3 {
                vx += 0.65 * (cdx / cd);
                vy += 0.65 * (cdy / cd);
            }
            let vlen = (vx * vx + vy * vy).sqrt();
            let (ux, uy) = if vlen > 1e-3 {
                (vx / vlen, vy / vlen)
            } else {
                (0.7, -0.7)
            };

            let txt = if carbon13 {
                format!("δ {:.1}", pk.shift_ppm)
            } else {
                format!("δ {:.2} ({})", pk.shift_ppm, pk.multiplicity)
            };
            let w = (txt.chars().count() as f64) * 5.2 + 8.0;
            let h = 12.0;
            let just = if ux > 0.28 {
                LabelJustification::Left
            } else if ux < -0.28 {
                LabelJustification::Right
            } else {
                LabelJustification::Center
            };
            let anchor_x = a.x + ux * 16.5;
            let anchor_y = a.y + uy * 15.0;
            let cx = match just {
                LabelJustification::Left => anchor_x + w * 0.5,
                LabelJustification::Right => anchor_x - w * 0.5,
                _ => anchor_x,
            };
            badges.push(PlacedBadge {
                text: txt,
                cx,
                cy: anchor_y,
                w,
                h,
                just,
            });
        }

        // Iteratively push badges clear of molecule atoms and other badges
        for _ in 0..20 {
            let mut moved = false;
            // 1. Repel from all atoms in the molecule so badges never sit on top of rings/bonds
            for b in &mut badges {
                for a in &mol.atoms {
                    let dx = b.cx - a.x;
                    let dy = b.cy - a.y;
                    let rx = b.w * 0.5 + 7.5;
                    let ry = b.h * 0.5 + 7.5;
                    if dx.abs() < rx && dy.abs() < ry {
                        let ox = rx - dx.abs();
                        let oy = ry - dy.abs();
                        if oy <= ox {
                            let sgn = if dy >= 0.0 { 1.0 } else { -1.0 };
                            b.cy += sgn * (oy * 0.65);
                        } else {
                            let sgn = if dx >= 0.0 { 1.0 } else { -1.0 };
                            b.cx += sgn * (ox * 0.65);
                        }
                        moved = true;
                    }
                }
            }
            // 2. Repel badge-vs-badge bounding boxes
            let n = badges.len();
            for i in 0..n {
                for j in (i + 1)..n {
                    let dx = badges[j].cx - badges[i].cx;
                    let dy = badges[j].cy - badges[i].cy;
                    let rx = 0.5 * (badges[i].w + badges[j].w) + 4.0;
                    let ry = 0.5 * (badges[i].h + badges[j].h) + 3.0;
                    if dx.abs() < rx && dy.abs() < ry {
                        let ox = rx - dx.abs();
                        let oy = ry - dy.abs();
                        if oy <= ox * 1.25 {
                            let sgn = if dy >= 0.0 { 1.0 } else { -1.0 };
                            let shift = oy * 0.55 + 0.5;
                            badges[i].cy -= sgn * shift;
                            badges[j].cy += sgn * shift;
                        } else {
                            let sgn = if dx >= 0.0 { 1.0 } else { -1.0 };
                            let shift = ox * 0.55 + 0.5;
                            badges[i].cx -= sgn * shift;
                            badges[j].cx += sgn * shift;
                        }
                        moved = true;
                    }
                }
            }
            if !moved {
                break;
            }
        }

        let mut stamped = 0usize;
        self.doc_mut().apply_scene_tx("Stamp NMR Shift Labels", |scene, next_id| {
            // Replace any previously stamped δ labels cleanly
            scene.captions.retain(|c| !c.text.starts_with('δ'));
            for b in &badges {
                let cid = *next_id;
                *next_id += 1;
                let x = match b.just {
                    LabelJustification::Left => b.cx - b.w * 0.5,
                    LabelJustification::Right => b.cx + b.w * 0.5,
                    _ => b.cx,
                };
                scene.captions.push(CaptionNode {
                    id: cid,
                    text: b.text.clone(),
                    x,
                    y: b.cy,
                    element: None,
                    generic_nickname: None,
                    warning: None,
                    justification: b.just,
                });
                stamped += 1;
            }
        });
        self.status_msg = format!("Stamped {stamped} {} NMR δ shift label(s) on canvas", pred.nucleus);
    }

    fn clear_nmr_shift_labels(&mut self) {
        self.doc_mut().apply_scene_tx("Clear NMR Shift Labels", |scene, _| {
            scene.captions.retain(|c| !c.text.starts_with('δ'));
        });
        self.status_msg = "Cleared stamped NMR δ labels".to_string();
    }

    fn execute_palette_command(&mut self, ctx: &egui::Context, act: PaletteCommandAction) {
        match act {
            PaletteCommandAction::None => {}
            PaletteCommandAction::RotateDeg(deg) => self.do_rotate_targets(deg),
            PaletteCommandAction::ScalePct(pct) => {
                if pct > 1.0 {
                    self.do_scale_targets(pct / 100.0);
                }
            }
            PaletteCommandAction::InsertSmiles(smi) => {
                match smiles_to_molecule(&smi) {
                    Ok(mut m) => {
                        let blen = self.doc().settings.bond_length_pt;
                        clean_up_structure_preserving_stereo(&mut m, &[], blen);
                        for a in &mut m.atoms { a.x += 180.0; a.y += 150.0; }
                        self.doc_mut().stamp_molecule_tx(m, "Insert SMILES");
                        self.status_msg = format!("Inserted SMILES: {smi}");
                    }
                    Err(e) => self.status_msg = format!("SMILES error: {e}"),
                }
            }
            PaletteCommandAction::InsertName(name) => {
                let blen = self.doc().settings.bond_length_pt;
                if let Ok(smi) = self.name_dict.name_to_smiles(&name) {
                    if let Ok(mut m) = smiles_to_molecule(&smi) {
                        clean_up_structure_preserving_stereo(&mut m, &[], blen);
                        for a in &mut m.atoms { a.x += 180.0; a.y += 150.0; }
                        self.doc_mut().stamp_molecule_tx(m, "Name to Structure");
                        self.status_msg = format!("Inserted \"{name}\" ({smi})");
                    }
                } else {
                    self.chem_db_state.query = name;
                    self.show_chem_db = true;
                    self.status_msg = "Name not in offline dictionary — opened PubChem/ChEMBL search".to_string();
                }
            }
            PaletteCommandAction::ApplyStyle(q) => {
                let qlow = q.to_ascii_lowercase();
                if let Some((lbl, st)) = self.stylesheets.iter().find(|(n, _)| n.to_ascii_lowercase().contains(&qlow)).cloned() {
                    replace_settings_tx(self.doc_mut(), st, "Apply Stylesheet");
                    self.do_clean_up_structure();
                    self.status_msg = format!("Applied stylesheet: {lbl}");
                } else {
                    self.status_msg = format!("No stylesheet matching \"{q}\"");
                }
            }
            PaletteCommandAction::SetTheme(t) => {
                let th = match t {
                    "dark" => UiTheme::DarkChrome,
                    "slide" => UiTheme::DarkSlide,
                    _ => UiTheme::ClassicLight,
                };
                self.ui_theme = th;
                Self::apply_theme(ctx, th);
                self.status_msg = format!("Switched theme to {t}");
            }
            PaletteCommandAction::CatalyticCycle(mode) => match mode {
                "circle" => self.arrange_selected_catalytic_cycle(false),
                "ellipse" => self.arrange_selected_catalytic_cycle(true),
                _ => {
                    self.doc_mut().apply_scene_tx("Insert Starter Catalytic Cycle", |s, nid| {
                        insert_starter_catalytic_cycle(s, nid, 250.0, 210.0);
                    });
                    self.status_msg = "Inserted 4-step catalytic cycle".to_string();
                }
            },
            PaletteCommandAction::PredictNmr { carbon13 } => self.do_predict_nmr(carbon13),
            PaletteCommandAction::StampNmrLabels { carbon13 } => self.stamp_nmr_shift_labels(carbon13),
            PaletteCommandAction::ClearNmrLabels => self.clear_nmr_shift_labels(),
            PaletteCommandAction::Open3dViewer => { self.viewer_3d.open = true; }
            PaletteCommandAction::OpenDiffCompare => {
                self.show_diff_window = true;
                if self.tabs.len() >= 2 {
                    self.diff_left_tab = 0;
                    self.diff_right_tab = 1;
                    self.last_diff = Some(compare_documents(&self.tabs[0].doc, &self.tabs[1].doc));
                }
            }
            PaletteCommandAction::OpenObjectInspector => { self.show_object_inspector = true; }
            PaletteCommandAction::OpenChemDb => { self.show_chem_db = true; }
            PaletteCommandAction::CleanUpStructure => self.do_clean_up_structure(),
            PaletteCommandAction::CleanUpReaction => {
                self.doc_mut().apply_scene_tx("Clean Up Reaction", |s, _| { clean_up_reaction(s); });
                self.status_msg = "Cleaned up reaction".to_string();
            }
            PaletteCommandAction::CheckStructure => {
                self.check_warn_idx = 0;
                self.dialog = ActiveDialog::CheckStructureDialog;
            }
            PaletteCommandAction::StructureToName => self.do_structure_to_name_hybrid(ctx, false),
            PaletteCommandAction::CopyEditableSvg => self.do_copy_editable_svg(ctx),
            PaletteCommandAction::Align(a) => match a {
                "left" => self.do_align_action("Left edges", AlignMode::Left),
                "center" => self.do_align_action("L/R centers", AlignMode::HorizontalCenter),
                "right" => self.do_align_action("Right edges", AlignMode::Right),
                "top" => self.do_align_action("Top edges", AlignMode::Top),
                _ => self.do_align_action("Bottom edges", AlignMode::Bottom),
            },
            PaletteCommandAction::SetExportDpi(dpi) => {
                self.export_opts.png_dpi = dpi.clamp(72, 1200);
                self.status_msg = format!("Set PNG export resolution to {} DPI", self.export_opts.png_dpi);
            }
            PaletteCommandAction::InsertSpecialCore(core) => {
                let blen = self.doc().settings.bond_length_pt;
                match core {
                    "pyranose" => self.insert_special_core("Haworth β-D-Glucopyranose", build_haworth_pyranose(blen)),
                    "furanose" => self.insert_special_core("Haworth β-D-Ribofuranose", build_haworth_furanose(blen)),
                    "porphine" => self.insert_special_core("Porphine", build_porphine_core(blen)),
                    "crown" => self.insert_special_core("18-Crown-6", build_crown_ether_18_6(blen)),
                    "adamantane" => self.insert_special_core("Adamantane", build_adamantane_2d(blen)),
                    _ => self.insert_special_core("Norbornane", build_norbornane_2d(blen)),
                }
            }
            PaletteCommandAction::ToggleStereo => { self.show_stereo = !self.show_stereo; }
            PaletteCommandAction::ToggleRulers => { self.show_rulers = !self.show_rulers; }
            PaletteCommandAction::OpenSaveAs(fname) => {
                self.file_path_buf = fname.to_string();
                self.dialog = ActiveDialog::SaveAs;
            }
            PaletteCommandAction::LinkCaption(kind) => {
                let k = if kind == "name" { LinkedCaptionKind::IupacName } else { LinkedCaptionKind::Formula };
                self.attach_linked_caption(k);
            }
            PaletteCommandAction::OpenMedChem => { self.pro_suite.show_medchem = true; }
            PaletteCommandAction::OpenHrms => { self.pro_suite.show_hrms = true; }
            PaletteCommandAction::OpenIsomers(tab) => {
                self.pro_suite.show_isomers = true;
                self.pro_suite.isomer_tab = match tab {
                    "tautomers" => IsomerTab::Tautomers,
                    "stereo" => IsomerTab::Stereoisomers,
                    "resonance" => IsomerTab::Resonance,
                    "push" => IsomerTab::ElectronPush,
                    _ => IsomerTab::PhMicrospecies,
                };
            }
            PaletteCommandAction::OpenMarkush => { self.pro_suite.show_markush = true; }
            PaletteCommandAction::OpenRetro => {
                self.pro_suite.show_retro = true;
                if let Some(m) = self.active_mol().cloned() {
                    let blen = self.doc().settings.bond_length_pt;
                    self.pro_suite.retro_hits = suggest_disconnections(&m, blen);
                    self.pro_suite.highlighted_retro_bond = self.pro_suite.retro_hits.first().map(|h| h.bond_idx);
                }
            }
            PaletteCommandAction::OpenPubComposer => { self.pro_suite.show_pub_composer = true; }
            PaletteCommandAction::OpenSpectroscopy(tab) => {
                self.omni_suite.show_spectroscopy = true;
                self.omni_suite.spec_tab = match tab {
                    "uv" => SpecTab::UvVis,
                    "2d" => SpecTab::Nmr2d,
                    _ => SpecTab::FtIr,
                };
            }
            PaletteCommandAction::OpenInoPoly(tab) => {
                self.omni_suite.show_ino_poly = true;
                self.omni_suite.ino_poly_tab = if tab == "polymer" {
                    InoPolyTab::Polymer
                } else {
                    InoPolyTab::Organometallic
                };
            }
            PaletteCommandAction::OpenEnergyProfile => { self.omni_suite.show_energy_profile = true; }
            PaletteCommandAction::OpenElnSafety => { self.omni_suite.show_eln_safety = true; }
            PaletteCommandAction::OpenLatexExport => { self.omni_suite.show_latex_export = true; }
            PaletteCommandAction::OpenScreening => { self.omni_suite.show_screening = true; }
            PaletteCommandAction::Optimize3dNow => {
                if let Some(mi) = self.active_mol_index() {
                    let mut summary = String::new();
                    self.doc_mut().apply_scene_tx("UFF 3D Optimize", |s, _| {
                        if let Some(m) = s.molecules.get_mut(mi) {
                            let lifted = lift_molecule_2d_to_3d(m);
                            let (opt, rep) = optimize_geometry_3d(&lifted, 180);
                            summary = format!("UFF 3D Opt ({} iters): E {:.1} → {:.2} kcal/mol", rep.iterations, rep.initial_energy_kcal, rep.final_energy_kcal);
                            *m = opt;
                        }
                    });
                    if !summary.is_empty() {
                        self.viewer_3d.last_opt_summary = summary.clone();
                        self.viewer_3d.open = true;
                        self.status_msg = summary;
                    }
                }
            }
            PaletteCommandAction::OpenQuantumBridge => { self.quantum.show_window = true; }
            PaletteCommandAction::ToggleStudioDock => {
                self.show_right_dock = !self.show_right_dock;
                self.persist_user_config();
            }
        }
    }

    fn persist_user_config(&mut self) {
        if cfg!(test) {
            return;
        }
        self.user_config.theme = match self.ui_theme {
            UiTheme::ClassicLight => "light".into(),
            UiTheme::DarkChrome => "dark".into(),
            UiTheme::DarkSlide => "slide".into(),
        };
        self.user_config.png_dpi = self.export_opts.png_dpi;
        self.user_config.svg_padding_pt = self.export_opts.svg_padding_pt;
        self.user_config.copy_bg = match self.export_opts.copy_bg {
            CopyBackground::Transparent => "transparent".into(),
            CopyBackground::White => "white".into(),
            CopyBackground::DarkSlide => "slide".into(),
        };
        self.user_config.rear_opacity = self.rear_opacity;
        self.user_config.crossing_clearance_pt = self.crossing_clearance_pt;
        self.user_config.show_right_dock = self.show_right_dock;
        self.user_config.save();
    }

    fn open_path_into_new_tab(&mut self, path: &str) {
        let low = path.to_ascii_lowercase();
        if low.ends_with(".jdx") || low.ends_with(".dx") {
            match fs::read_to_string(path).map_err(|e| e.to_string()).and_then(|s| parse_jcamp_dx(&s)) {
                Ok(spec) => {
                    self.status_msg = format!("Loaded experimental JCAMP-DX overlay: {}", spec.title);
                    self.jcamp_overlay = Some(spec);
                    self.show_nmr_window = true;
                }
                Err(e) => self.status_msg = format!("JCAMP-DX open error: {e}"),
            }
            return;
        }
        if low.ends_with(".out") || low.ends_with(".log") {
            if let Ok(txt) = fs::read_to_string(path) {
                if let Ok(qrep) = parse_quantum_output(&txt) {
                    self.quantum.parsed_report = Some(qrep);
                    self.quantum.import_path = path.to_string();
                    self.quantum.show_window = true;
                }
            }
        }
        match load_document_any_with_cell(path) {
            Ok((d, cell_opt)) => {
                let idx = self.tabs.len();
                let mut tab = OpenDocTab::new_at(d, idx);
                tab.file_path = Some(path.to_string());
                self.tabs.push(tab);
                self.active_tab = self.tabs.len() - 1;
                self.reset_transient_interaction_state();
                if let Some(cell) = cell_opt {
                    self.viewer_3d.unit_cell = Some(cell);
                    self.viewer_3d.open = true;
                }
                if !cfg!(test) {
                    self.user_config.push_recent(path);
                }
                self.status_msg = format!("Opened {path}");
            }
            Err(e) => self.status_msg = format!("Open failed ({path}): {e}"),
        }
    }

    fn handle_dock_quantum_action(&mut self, ctx: &egui::Context, act: DockQuantumAction) {
        match act {
            DockQuantumAction::None => {}
            DockQuantumAction::CopyText(txt, msg) => {
                ctx.copy_text(txt);
                self.status_msg = msg;
            }
            DockQuantumAction::SaveExportFile { title, suggested_name, content } => {
                match native_dialogs::save_text_with_dialog(title, &suggested_name, &content) {
                    Ok(Some(path)) => self.status_msg = format!("Saved {path}"),
                    Ok(None) => self.status_msg = "Save cancelled".into(),
                    Err(e) => self.status_msg = format!("Save error: {e}"),
                }
            }
            DockQuantumAction::StampOptimized3d(m) => {
                self.doc_mut().stamp_molecule_tx(m, "DFT Optimized 3D Geometry");
                self.viewer_3d.open = true;
                self.status_msg = "Stamped DFT-optimized 3D geometry onto canvas".into();
            }
            DockQuantumAction::PlotQuantumIr(ir_pred) => {
                self.omni_suite.ir_render.scale_factor = ir_pred.default_scale_factor;
                self.omni_suite.ir_render.use_per_band_widths = false;
                self.omni_suite.ir_engine = IrEngineMode::QuantumImported;
                self.omni_suite.ir_quantum_mol_sig = self
                    .active_mol()
                    .map(molecule_to_smiles)
                    .unwrap_or_default();
                self.omni_suite.ir_status_msg = format!(
                    "Plotted {} ({} modes)",
                    ir_pred.source.label(),
                    ir_pred.bands.len()
                );
                self.omni_suite.ir_quantum_pred = Some(ir_pred);
                self.omni_suite.spec_tab = SpecTab::FtIr;
                self.omni_suite.show_spectroscopy = true;
                self.status_msg = "Loaded quantum harmonic IR into Multi-Spectroscopy Studio".into();
            }
            DockQuantumAction::RunLocalQc(backend) => {
                self.launch_local_quantum_ir(backend);
            }
            DockQuantumAction::RedetectQcBinaries => {
                self.qc_runner.redetect(
                    &self.user_config.xtb_path,
                    &self.user_config.orca_path,
                    &self.user_config.gaussian_path,
                );
                self.persist_user_config();
                self.status_msg = format!(
                    "QC binaries: xTB={} | ORCA={}",
                    if self.qc_runner.binaries.xtb_path.is_some() { "OK" } else { "none" },
                    if self.qc_runner.binaries.orca_path.is_some() { "OK" } else { "none" }
                );
            }
            DockQuantumAction::OpenFloatingStudio(which) => match which {
                "3d" => { self.viewer_3d.open = true; }
                "isomers" => { self.pro_suite.show_isomers = true; }
                "retro" => {
                    self.pro_suite.show_retro = true;
                    if let Some(m) = self.active_mol().cloned() {
                        let blen = self.doc().settings.bond_length_pt;
                        self.pro_suite.retro_hits = suggest_disconnections(&m, blen);
                        self.pro_suite.highlighted_retro_bond = self.pro_suite.retro_hits.first().map(|h| h.bond_idx);
                    }
                }
                "latex" => { self.omni_suite.show_latex_export = true; }
                "medchem" => { self.pro_suite.show_medchem = true; }
                "markush" => { self.pro_suite.show_markush = true; }
                "screen" => { self.omni_suite.show_screening = true; }
                "spec" => { self.omni_suite.show_spectroscopy = true; }
                "hrms" => { self.pro_suite.show_hrms = true; }
                "nmr1h" => { self.do_predict_nmr(false); }
                "nmr13c" => { self.do_predict_nmr(true); }
                "inspector" => { self.show_object_inspector = true; }
                "energy" => { self.omni_suite.show_energy_profile = true; }
                "toc" => { self.pro_suite.show_pub_composer = true; }
                "eln" => { self.omni_suite.show_eln_safety = true; }
                _ => {}
            },
            DockQuantumAction::ConfigChanged => {
                self.persist_user_config();
            }
        }
    }

    fn handle_omni_suite_action(&mut self, ctx: &egui::Context, act: OmniSuiteAction) {
        match act {
            OmniSuiteAction::None => {}
            OmniSuiteAction::StampMolecule(mut m, label) => {
                let (x1, y1, x2, y2) = mol_bbox(&m);
                let dx = 220.0 - 0.5 * (x1 + x2);
                let dy = 180.0 - 0.5 * (y1 + y2);
                for a in &mut m.atoms {
                    a.x += dx;
                    a.y += dy;
                }
                self.doc_mut().stamp_molecule_tx(m, label);
                self.status_msg = format!("Stamped: {label}");
            }
            OmniSuiteAction::StampEnergyProfile(mut prof) => {
                self.doc_mut().apply_scene_tx("Stamp Energy Profile", |scene, next_id| {
                    prof.id = *next_id;
                    *next_id += 1;
                    scene.energy_profiles.push(prof);
                });
                self.status_msg = "Stamped Reaction Energy Profile onto canvas".to_string();
            }
            OmniSuiteAction::PasteCaption(txt) => self.paste_caption_tx(txt),
            OmniSuiteAction::CopyText(txt, msg) => {
                ctx.copy_text(txt);
                self.status_msg = msg;
            }
            OmniSuiteAction::SaveExportFile { title, suggested_name, content } => {
                match native_dialogs::save_text_with_dialog(title, &suggested_name, &content) {
                    Ok(Some(path)) => self.status_msg = format!("Exported {path}"),
                    Ok(None) => self.status_msg = "Export cancelled".into(),
                    Err(e) => self.status_msg = format!("Export error: {e}"),
                }
            }
            OmniSuiteAction::StatusMsg(m) => self.status_msg = m,
            OmniSuiteAction::RunQuantumIr(backend) => {
                self.launch_local_quantum_ir(backend);
            }
            OmniSuiteAction::CancelQuantumIr => {
                self.qc_runner.cancel();
                self.status_msg = "Cancelled local quantum FTIR job".into();
            }
            OmniSuiteAction::ImportQuantumIrFile(path) => {
                match fs::read_to_string(&path)
                    .map_err(|e| e.to_string())
                    .and_then(|s| parse_quantum_output(&s))
                {
                    Ok(rep) => {
                        if let Some(ir_pred) = rep.to_ir_prediction(None) {
                            let n_modes = ir_pred.bands.len();
                            let src = ir_pred.source.label();
                            self.omni_suite.ir_render.scale_factor = ir_pred.default_scale_factor;
                            self.omni_suite.ir_render.use_per_band_widths = false;
                            self.omni_suite.ir_engine = IrEngineMode::QuantumImported;
                            self.omni_suite.ir_quantum_mol_sig = self
                                .active_mol()
                                .map(molecule_to_smiles)
                                .unwrap_or_default();
                            self.omni_suite.ir_selected_band = None;
                            self.omni_suite.ir_status_msg = format!("Loaded {src}: {n_modes} harmonic modes");
                            self.omni_suite.ir_quantum_pred = Some(ir_pred);
                            self.quantum.parsed_report = Some(rep);
                            self.status_msg = format!("Imported {src} IR ({n_modes} modes) from {path}");
                        } else {
                            self.omni_suite.ir_status_msg = "No IR vibrational modes found in quantum output".into();
                            self.status_msg = self.omni_suite.ir_status_msg.clone();
                        }
                    }
                    Err(e) => {
                        self.omni_suite.ir_status_msg = format!("QC import error: {e}");
                        self.status_msg = self.omni_suite.ir_status_msg.clone();
                    }
                }
            }
            OmniSuiteAction::LoadExpIrFile(path) => {
                match fs::read_to_string(&path)
                    .map_err(|e| e.to_string())
                    .and_then(|s| parse_experimental_ir_file(&path, &s))
                {
                    Ok(spec) => {
                        self.status_msg = format!("Overlaid experimental IR: {} ({} pts)", spec.title, spec.raw_points.len());
                        self.omni_suite.ir_exp_overlay = Some(spec);
                    }
                    Err(e) => {
                        self.status_msg = format!("Experimental IR load failed: {e}");
                    }
                }
            }
            OmniSuiteAction::FetchNistIr(query) => {
                let inchi = self.active_mol().map(molecule_to_inchi);
                self.nist_ir_fetch.start_fetch(query, inchi);
                self.status_msg = "Querying NIST Chemistry WebBook for experimental IR...".into();
            }
            OmniSuiteAction::Animate3dVibMode {
                geom_3d,
                wavenumber_cm1,
                label,
                vectors,
            } => {
                self.viewer_3d.open = true;
                self.viewer_3d.vib_anim = Some(VibrationalAnimState {
                    base_mol_3d: geom_3d,
                    wavenumber_cm1,
                    label: label.clone(),
                    mode_vectors: vectors,
                    amplitude_pt: 6.5,
                    speed_hz: 1.6,
                    playing: true,
                    show_vectors: true,
                });
                self.status_msg = format!("Animating 3D normal mode: {wavenumber_cm1:.1} cm^-1 ({label})");
            }
        }
    }

    fn launch_local_quantum_ir(&mut self, backend: QcRunBackend) {
        self.qc_runner.ensure_detected(
            &self.user_config.xtb_path,
            &self.user_config.orca_path,
            &self.user_config.gaussian_path,
        );
        let Some(mol) = self.active_mol().cloned() else {
            self.status_msg = "Draw or select a molecule first".into();
            return;
        };
        let bin_opt = match backend {
            QcRunBackend::XtbGfn2 | QcRunBackend::XtbGfnFf => self.qc_runner.binaries.xtb_path.clone(),
            QcRunBackend::OrcaDft => self.qc_runner.binaries.orca_path.clone(),
            QcRunBackend::GaussianDft => self.qc_runner.binaries.gaussian_path.clone(),
        };
        let Some(bin_path) = bin_opt else {
            self.status_msg = format!("{} binary not found — configure path in DFT dock", backend.label());
            return;
        };
        let smi = molecule_to_smiles(&mol);
        let spec = self.quantum.spec.clone();
        let cache_key = format!(
            "{smi}|{backend:?}|{}|{}|{}|{}|{}",
            spec.method, spec.basis, spec.solvent, spec.charge, spec.multiplicity
        );
        let lifted = lift_molecule_2d_to_3d(&mol);
        self.qc_runner.start_ir_job(backend, &lifted, cache_key, spec, bin_path);
        if let Some(hit) = self.qc_runner.last_outcome.clone().filter(|_| !self.qc_runner.running) {
            self.omni_suite.ir_render.scale_factor = hit.ir_prediction.default_scale_factor;
            self.omni_suite.ir_render.use_per_band_widths = false;
            self.omni_suite.ir_engine = IrEngineMode::QuantumImported;
            self.omni_suite.ir_status_msg = format!(
                "Cached {} ({} modes)",
                hit.ir_prediction.source.label(),
                hit.ir_prediction.bands.len()
            );
            self.omni_suite.ir_quantum_mol_sig = smi;
            self.omni_suite.ir_quantum_pred = Some(hit.ir_prediction);
            self.quantum.parsed_report = Some(hit.report);
            self.status_msg = "Loaded cached quantum FTIR spectrum".into();
        } else {
            self.omni_suite.show_spectroscopy = true;
            self.omni_suite.spec_tab = SpecTab::FtIr;
            self.status_msg = format!("Running {} harmonic FTIR in background...", backend.label());
        }
    }

    fn handle_pro_suite_action(&mut self, act: ProSuiteAction) {
        match act {
            ProSuiteAction::None => {}
            ProSuiteAction::ReplaceActiveMolecule(m, label) => {
                if let Some(mi) = self.active_mol_index() {
                    self.doc_mut().apply_scene_tx(label, |s, _| {
                        if let Some(slot) = s.molecules.get_mut(mi) {
                            *slot = m;
                        }
                    });
                    self.status_msg = format!("Applied: {label}");
                }
            }
            ProSuiteAction::StampMolecule(mut m, label) => {
                let (x1, y1, x2, y2) = mol_bbox(&m);
                let dx = 280.0 - 0.5 * (x1 + x2);
                let dy = 210.0 - 0.5 * (y1 + y2);
                for a in &mut m.atoms {
                    a.x += dx;
                    a.y += dy;
                }
                self.doc_mut().stamp_molecule_tx(m, label);
                self.status_msg = format!("Stamped: {label}");
            }
            ProSuiteAction::StampMoleculeGrid(items, label) => {
                let cols = 3usize;
                self.doc_mut().apply_scene_tx(label, |scene, next_id| {
                    for (idx, (cap_txt, mut m)) in items.into_iter().enumerate() {
                        let gx = 105.0 + ((idx % cols) as f64) * 140.0;
                        let gy = 115.0 + ((idx / cols) as f64) * 115.0;
                        let (x1, y1, x2, y2) = mol_bbox(&m);
                        let dx = gx - 0.5 * (x1 + x2);
                        let dy = gy - 0.5 * (y1 + y2);
                        for a in &mut m.atoms {
                            a.x += dx;
                            a.y += dy;
                        }
                        let (_, _, _, my2) = mol_bbox(&m);
                        scene.molecules.push(m);
                        let cid = *next_id; *next_id += 1;
                        scene.captions.push(CaptionNode {
                            id: cid,
                            text: cap_txt,
                            x: gx,
                            y: my2 + 13.0,
                            element: None,
                            generic_nickname: None,
                            warning: None,
                            justification: LabelJustification::Center,
                        });
                    }
                });
                self.status_msg = format!("Stamped grid: {label}");
            }
            ProSuiteAction::StampRetrosynthesis(hit) => {
                if let Some(target) = self.active_mol().cloned() {
                    self.doc_mut().apply_scene_tx("Stamp Retrosynthesis", |scene, next_id| {
                        stamp_retrosynthesis_scheme(scene, next_id, &target, &hit, 110.0, 260.0);
                    });
                    self.status_msg = format!("Stamped retrosynthesis: {}", hit.reaction_name);
                }
            }
            ProSuiteAction::LoadExampleTarget(name, smi) => {
                let blen = self.doc().settings.bond_length_pt;
                match smiles_to_molecule(smi) {
                    Ok(mut m) => {
                        clean_up_structure_preserving_stereo(&mut m, &[], blen);
                        let (x1, y1, x2, y2) = mol_bbox(&m);
                        let dx = 220.0 - 0.5 * (x1 + x2);
                        let dy = 175.0 - 0.5 * (y1 + y2);
                        for a in &mut m.atoms {
                            a.x += dx;
                            a.y += dy;
                        }
                        self.doc_mut().stamp_molecule_tx(m.clone(), &format!("Load Example: {name}"));
                        let new_mi = self.doc().scene.molecules.len().saturating_sub(1);
                        self.tab_mut().selected_mols = vec![new_mi];
                        self.pro_suite.retro_hits = suggest_disconnections(&m, blen);
                        self.pro_suite.highlighted_retro_bond =
                            self.pro_suite.retro_hits.first().map(|h| h.bond_idx);
                        self.pro_suite.retro_mol_sig.clear();
                        self.status_msg = format!(
                            "Loaded example target {name} ({} disconnection(s) found)",
                            self.pro_suite.retro_hits.len()
                        );
                    }
                    Err(e) => {
                        self.status_msg = format!("Failed to load example {name}: {e}");
                    }
                }
            }
            ProSuiteAction::InsertTocFrame(preset) => {
                let mut bbox = (40.0, 40.0, 200.0, 160.0);
                self.doc_mut().apply_scene_tx("Insert TOC Graphic Frame", |scene, next_id| {
                    bbox = insert_toc_frame(scene, next_id, preset, 45.0, 55.0);
                });
                self.pro_suite.last_toc_bbox = Some(bbox);
                self.status_msg = format!("Inserted {}", preset.name);
            }
            ProSuiteAction::StampArrowConditions { above, below } => {
                let mut ok = false;
                self.doc_mut().apply_scene_tx("Stamp Arrow Conditions", |scene, next_id| {
                    if scene.arrows.is_empty() {
                        let aid = *next_id; *next_id += 1;
                        scene.arrows.push(SceneArrow {
                            id: aid,
                            tail: (160.0, 175.0),
                            head: (240.0, 175.0),
                            arrowhead_head: "Full".to_string(),
                            arrowhead_tail: None,
                            no_go: None,
                            dipole: false,
                        });
                    }
                    let last = scene.arrows.len() - 1;
                    ok = stamp_arrow_conditions(scene, next_id, last, &above, &below);
                });
                self.status_msg = if ok { "Stamped reaction conditions on arrow".into() } else { "No arrow available".into() };
            }
            ProSuiteAction::WrapSerpentine(per_row) => {
                let mut count = 0usize;
                self.doc_mut().apply_scene_tx("Serpentine Scheme Layout", |scene, next_id| {
                    count = wrap_serpentine_reaction_scheme(scene, next_id, per_row);
                });
                self.status_msg = format!("Wrapped {count} structure(s) into serpentine reaction scheme");
            }
            ProSuiteAction::ExportTextFile { suggested_name, content } => {
                match native_dialogs::save_text_with_dialog("Export File", &suggested_name, &content) {
                    Ok(Some(path)) => self.status_msg = format!("Exported {path}"),
                    Ok(None) => self.status_msg = "Export cancelled".into(),
                    Err(e) => self.status_msg = format!("Export error: {e}"),
                }
            }
            ProSuiteAction::StatusMsg(msg) => {
                self.status_msg = msg;
            }
        }
    }

    fn screen_to_doc(&self, page_origin: Pos2, pos: Pos2) -> (f64, f64) {
        let x = ((pos.x - page_origin.x) / self.zoom) as f64;
        let y = ((pos.y - page_origin.y) / self.zoom) as f64;
        (x, y)
    }

    fn doc_to_screen(&self, page_origin: Pos2, x: f64, y: f64) -> Pos2 {
        Pos2::new(
            page_origin.x + (x as f32) * self.zoom,
            page_origin.y + (y as f32) * self.zoom,
        )
    }

    /// Zoom-aware hit radius for atom vertices (~10 screen px, clamped to sensible doc-pt bounds).
    fn atom_hit_tol(&self) -> f64 {
        (10.0 / (self.zoom as f64).max(0.4)).clamp(4.5, 16.0)
    }

    /// Zoom-aware hit radius for bond segments (~8.5 screen px, clamped to sensible doc-pt bounds).
    fn bond_hit_tol(&self) -> f64 {
        (8.5 / (self.zoom as f64).max(0.4)).clamp(4.0, 14.0)
    }

    fn hit_atom(&self, dx: f64, dy: f64, tol: f64) -> Option<(usize, usize)> {
        let mut best: Option<((usize, usize), f64)> = None;
        for (mi, mol) in self.doc().scene.molecules.iter().enumerate() {
            for (ai, a) in mol.atoms.iter().enumerate() {
                if !a.visible {
                    continue;
                }
                let d = (a.x - dx).hypot(a.y - dy);
                if d <= tol && best.map(|(_, bd)| d < bd).unwrap_or(true) {
                    best = Some(((mi, ai), d));
                }
            }
        }
        best.map(|(idx, _)| idx)
    }

    fn hit_bond(&self, dx: f64, dy: f64, tol: f64) -> Option<(usize, usize)> {
        for (mi, mol) in self.doc().scene.molecules.iter().enumerate() {
            for (bi, b) in mol.bonds.iter().enumerate() {
                if !b.visible || b.begin_atom >= mol.atoms.len() || b.end_atom >= mol.atoms.len() {
                    continue;
                }
                let a1 = &mol.atoms[b.begin_atom];
                let a2 = &mol.atoms[b.end_atom];
                let vx = a2.x - a1.x;
                let vy = a2.y - a1.y;
                let l2 = vx * vx + vy * vy;
                if l2 < 1e-6 {
                    continue;
                }
                let t = (((dx - a1.x) * vx + (dy - a1.y) * vy) / l2).clamp(0.18, 0.82);
                let px = a1.x + t * vx;
                let py = a1.y + t * vy;
                if (px - dx).hypot(py - dy) <= tol {
                    return Some((mi, bi));
                }
            }
        }
        None
    }

    fn hit_molecule_bbox(&self, dx: f64, dy: f64) -> Option<usize> {
        for (mi, mol) in self.doc().scene.molecules.iter().enumerate().rev() {
            if mol.atoms.is_empty() {
                continue;
            }
            let (x1, y1, x2, y2) = mol_bbox(mol);
            if dx >= x1 - 6.0 && dx <= x2 + 6.0 && dy >= y1 - 6.0 && dy <= y2 + 6.0 {
                return Some(mi);
            }
        }
        None
    }

    fn hit_caption(&self, dx: f64, dy: f64) -> Option<usize> {
        for (ci, c) in self.doc().scene.captions.iter().enumerate().rev() {
            let w = (c.text.chars().count() as f64 * 5.6).max(24.0);
            let (x1, x2) = match c.justification {
                LabelJustification::Center | LabelJustification::Above | LabelJustification::Below => {
                    (c.x - 0.5 * w - 5.0, c.x + 0.5 * w + 5.0)
                }
                LabelJustification::Right => (c.x - w - 5.0, c.x + 5.0),
                LabelJustification::Left | LabelJustification::Auto => (c.x - 5.0, c.x + w + 5.0),
            };
            if dx >= x1 && dx <= x2 && (dy - c.y).abs() <= 9.0 {
                return Some(ci);
            }
        }
        None
    }

    fn expand_selection_with_groups(&mut self) {
        let mut expanded = self.tab().selected_mols.clone();
        let mut changed = true;
        while changed {
            changed = false;
            for grp in &self.tab().groups {
                if grp.iter().any(|mi| expanded.contains(mi)) {
                    for &mi in grp {
                        if !expanded.contains(&mi) {
                            expanded.push(mi);
                            changed = true;
                        }
                    }
                }
            }
        }
        self.tab_mut().selected_mols = expanded;
    }

    fn hit_selected_molecule(&self, dx: f64, dy: f64) -> bool {
        for &mi in &self.tab().selected_mols {
            if let Some(mol) = self.doc().scene.molecules.get(mi) {
                let (x1, y1, x2, y2) = mol_bbox(mol);
                if dx >= x1 - 6.0 && dx <= x2 + 6.0 && dy >= y1 - 6.0 && dy <= y2 + 6.0 {
                    return true;
                }
            }
        }
        false
    }

    fn selection_bbox(&self) -> Option<(f64, f64, f64, f64)> {
        let indices: Vec<usize> = if self.tab().selected_mols.is_empty() {
            (0..self.doc().scene.molecules.len()).collect()
        } else {
            self.tab().selected_mols.clone()
        };
        let mut x1 = f64::INFINITY;
        let mut y1 = f64::INFINITY;
        let mut x2 = f64::NEG_INFINITY;
        let mut y2 = f64::NEG_INFINITY;
        for mi in indices {
            if let Some(m) = self.doc().scene.molecules.get(mi) {
                for a in &m.atoms {
                    x1 = x1.min(a.x);
                    y1 = y1.min(a.y);
                    x2 = x2.max(a.x);
                    y2 = y2.max(a.y);
                }
            }
        }
        if x1.is_finite() { Some((x1, y1, x2, y2)) } else { None }
    }

    fn current_bond_spec(&self) -> (i32, BondDisplay) {
        match self.tool {
            Tool::MultipleBond | Tool::DoubleBond => (2, BondDisplay::Solid),
            Tool::TripleBond => (3, BondDisplay::Solid),
            Tool::QuadrupleBond => (4, BondDisplay::Solid),
            Tool::WedgeBond => (1, BondDisplay::WedgeBegin),
            Tool::HashWedgeBond => (1, BondDisplay::WedgedHashBegin),
            Tool::HollowWedgeBond => (1, BondDisplay::HollowWedgeBegin),
            Tool::BoldBond => (1, BondDisplay::Bold),
            Tool::DashedBond => (1, BondDisplay::Dash),
            Tool::HashedBond => (1, BondDisplay::Hash),
            Tool::WavyBond => (1, BondDisplay::Wavy),
            Tool::DativeBond => (1, BondDisplay::Dative),
            _ => (1, BondDisplay::Solid),
        }
    }

    fn is_bond_tool(&self) -> bool {
        matches!(
            self.tool,
            Tool::SolidBond
                | Tool::MultipleBond
                | Tool::DoubleBond
                | Tool::TripleBond
                | Tool::QuadrupleBond
                | Tool::DativeBond
                | Tool::DashedBond
                | Tool::HashedBond
                | Tool::HashWedgeBond
                | Tool::BoldBond
                | Tool::WedgeBond
                | Tool::HollowWedgeBond
                | Tool::WavyBond
        )
    }

    fn ring_spec(tool: Tool) -> Option<(usize, bool, bool)> {
        match tool {
            Tool::Cyclopropane => Some((3, false, false)),
            Tool::Cyclobutane => Some((4, false, false)),
            Tool::Cyclopentane => Some((5, false, false)),
            Tool::Cyclohexane => Some((6, false, false)),
            Tool::Cycloheptane => Some((7, false, false)),
            Tool::Cyclooctane => Some((8, false, false)),
            Tool::Cyclopentadiene => Some((5, false, true)),
            Tool::Benzene => Some((6, true, false)),
            _ => None,
        }
    }

    fn build_polygon_ring_angled(
        n_sides: usize,
        aromatic: bool,
        cpd: bool,
        cx: f64,
        cy: f64,
        blen: f64,
        rot_ang: f64,
    ) -> ChemMolecule {
        let mut p = ChemMolecule::default();
        let r = blen / (2.0 * (std::f64::consts::PI / (n_sides as f64)).sin());
        for i in 0..n_sides {
            let ang = (i as f64) * std::f64::consts::TAU / (n_sides as f64) + rot_ang;
            let mut a = ChemAtom::new((i + 1) as u32, 6);
            a.x = cx + r * ang.cos();
            a.y = cy + r * ang.sin();
            p.atoms.push(a);
        }
        for i in 0..n_sides {
            let ord = if aromatic && i % 2 == 0 {
                2
            } else if cpd && (i == 1 || i == 3) {
                2
            } else {
                1
            };
            p.bonds.push(ChemBond {
                id: (i + 1) as u32,
                begin_atom: i,
                end_atom: (i + 1) % n_sides,
                order: ord,
                display: BondDisplay::Solid,
                visible: true,
                ignore_warnings: false,
            });
        }
        p
    }

    #[cfg(test)]
    fn build_polygon_ring(n_sides: usize, aromatic: bool, cpd: bool, cx: f64, cy: f64, blen: f64) -> ChemMolecule {
        Self::build_polygon_ring_angled(n_sides, aromatic, cpd, cx, cy, blen, -std::f64::consts::FRAC_PI_2)
    }

    /// Computes the ideal next bond angle off atom `ai` in `mol`:
    /// - `deg == 0`: standard ACS -30 deg
    /// - `deg == 1`: linear 180 deg for alkynes/nitriles, otherwise proper alternating zig-zag
    /// - `deg >= 2`: exact midpoint of the largest angular gap between existing bonds (never overlaps!)
    fn best_new_bond_angle(mol: &ChemMolecule, ai: usize, ord: i32, chain_rad: f64) -> f64 {
        let ax = mol.atoms[ai].x;
        let ay = mol.atoms[ai].y;
        let mut nbrs: Vec<(usize, f64, i32)> = Vec::new();
        for b in &mol.bonds {
            let nb = if b.begin_atom == ai {
                Some(b.end_atom)
            } else if b.end_atom == ai {
                Some(b.begin_atom)
            } else {
                None
            };
            if let Some(n_i) = nb {
                if n_i < mol.atoms.len() {
                    let ang = (mol.atoms[n_i].y - ay).atan2(mol.atoms[n_i].x - ax);
                    nbrs.push((n_i, ang, b.order));
                }
            }
        }
        if nbrs.is_empty() {
            return -std::f64::consts::FRAC_PI_6;
        }
        if nbrs.len() == 1 {
            let (nb, ang_nb, b_ord) = nbrs[0];
            if b_ord == 3 || ord == 3 {
                return ang_nb + std::f64::consts::PI;
            }
            // Check neighbor's own neighbors to alternate zig-zag properly along chains!
            let mut nb_other: Option<(f64, f64)> = None;
            for b in &mol.bonds {
                let o = if b.begin_atom == nb && b.end_atom != ai {
                    Some(b.end_atom)
                } else if b.end_atom == nb && b.begin_atom != ai {
                    Some(b.begin_atom)
                } else {
                    None
                };
                if let Some(oi) = o {
                    if oi < mol.atoms.len() {
                        nb_other = Some((mol.atoms[oi].x, mol.atoms[oi].y));
                        break;
                    }
                }
            }
            let cand_plus = ang_nb + chain_rad;
            let cand_minus = ang_nb - chain_rad;
            if let Some((ox, oy)) = nb_other {
                // Vector from nb -> ai and vector from nb -> nb_other
                let ux = ax - mol.atoms[nb].x;
                let uy = ay - mol.atoms[nb].y;
                let vx = ox - mol.atoms[nb].x;
                let vy = oy - mol.atoms[nb].y;
                let cross_prev = ux * vy - uy * vx;
                let c1 = ux * cand_plus.sin() - uy * cand_plus.cos();
                // Want new bond on the SAME side as nb_other for trans zig-zag (opposite cross sign to incoming)
                // Wait: in trans zig-zag (A-B-C-D), B->A and C->D point to OPPOSITE sides of B->C!
                // Here ux,uy = B->C, vx,vy = B->A. Cross(B->C, B->A) and Cross(B->C, C->D) should have OPPOSITE signs!
                // Note: Cross(u, d) = ux * dy - uy * dx = ux * sin(cand) - uy * cos(cand).
                if cross_prev.abs() > 1e-4 {
                    return if (c1 * cross_prev) < 0.0 { cand_plus } else { cand_minus };
                }
            }
            // No grandparent (2-atom ethane): pick whichever candidate stays closer to horizontal (+0 / PI)
            if cand_plus.cos().abs() >= cand_minus.cos().abs() {
                cand_plus
            } else {
                cand_minus
            }
        } else {
            // 2+ neighbors: find largest angular gap in [0, 2*PI) and bisect it!
            let tau = std::f64::consts::TAU;
            let mut angles: Vec<f64> = nbrs
                .iter()
                .map(|&(_, a, _)| a.rem_euclid(tau))
                .collect();
            angles.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let mut best_gap = -1.0_f64;
            let mut best_mid = 0.0_f64;
            for i in 0..angles.len() {
                let a0 = angles[i];
                let a1 = if i + 1 < angles.len() {
                    angles[i + 1]
                } else {
                    angles[0] + tau
                };
                let gap = a1 - a0;
                if gap > best_gap {
                    best_gap = gap;
                    best_mid = (a0 + 0.5 * gap).rem_euclid(tau);
                }
            }
            best_mid
        }
    }

    /// Computes vertices for a spiro ring fused at atom `(mi, ai)` with optional directional target `(tx, ty)`.
    fn compute_spiro_ring_vertices(
        mol: &ChemMolecule,
        ai: usize,
        n_sides: usize,
        blen: f64,
        dir_target: Option<(f64, f64)>,
    ) -> Vec<(f64, f64)> {
        let ax = mol.atoms[ai].x;
        let ay = mol.atoms[ai].y;
        let dir_ang = if let Some((tx, ty)) = dir_target.filter(|&(tx, ty)| (tx - ax).hypot(ty - ay) > 3.0) {
            (ty - ay).atan2(tx - ax)
        } else {
            Self::best_new_bond_angle(mol, ai, 1, 120.0_f64.to_radians())
        };
        let r = blen / (2.0 * (std::f64::consts::PI / (n_sides as f64)).sin());
        let rcx = ax + r * dir_ang.cos();
        let rcy = ay + r * dir_ang.sin();
        let start_ang = (ay - rcy).atan2(ax - rcx);
        let step = std::f64::consts::TAU / (n_sides as f64);
        (0..n_sides)
            .map(|k| {
                let ang = start_ang + (k as f64) * step;
                (rcx + r * ang.cos(), rcy + r * ang.sin())
            })
            .collect()
    }

    /// Computes vertices for an edge-fused ring across bond `(mi, bi)` with optional side hint `(hx, hy)`.
    fn compute_edge_fused_ring_vertices(
        mol: &ChemMolecule,
        bi: usize,
        n_sides: usize,
        blen: f64,
        side_hint: Option<(f64, f64)>,
    ) -> Vec<(f64, f64)> {
        let b = &mol.bonds[bi];
        let p1 = (mol.atoms[b.begin_atom].x, mol.atoms[b.begin_atom].y);
        let p2 = (mol.atoms[b.end_atom].x, mol.atoms[b.end_atom].y);
        let mx = 0.5 * (p1.0 + p2.0);
        let my = 0.5 * (p1.1 + p2.1);
        let raw_len = (p2.0 - p1.0).hypot(p2.1 - p1.1);
        let edge_len = raw_len.max(blen);
        let (ex, ey) = if raw_len > 1e-6 {
            ((p2.0 - p1.0) / raw_len, (p2.1 - p1.1) / raw_len)
        } else {
            (1.0, 0.0)
        };
        let mut nx = -ey;
        let mut ny = ex;

        if let Some((hx, hy)) = side_hint {
            let dot = (hx - mx) * nx + (hy - my) * ny;
            if dot < -0.5 {
                nx = -nx;
                ny = -ny;
            }
        } else {
            // Default: place on the side away from the molecule centroid
            let cx = mol.atoms.iter().map(|a| a.x).sum::<f64>() / (mol.atoms.len().max(1) as f64);
            let cy = mol.atoms.iter().map(|a| a.y).sum::<f64>() / (mol.atoms.len().max(1) as f64);
            if (mx + nx - cx).hypot(my + ny - cy) < (mx - nx - cx).hypot(my - ny - cy) {
                nx = -nx;
                ny = -ny;
            }
        }
        let build_side = |sx_n: f64, sy_n: f64| -> Vec<(f64, f64)> {
            let apothem = edge_len / (2.0 * (std::f64::consts::PI / (n_sides as f64)).tan());
            let r = edge_len / (2.0 * (std::f64::consts::PI / (n_sides as f64)).sin());
            let rcx = mx + sx_n * apothem;
            let rcy = my + sy_n * apothem;
            let a_start = (p2.1 - rcy).atan2(p2.0 - rcx);
            let mut step = std::f64::consts::TAU / (n_sides as f64);
            let test_ang = a_start + step;
            if (rcx + r * test_ang.cos() - p1.0).hypot(rcy + r * test_ang.sin() - p1.1) < 1.0 {
                step = -step;
            }
            let mut verts = Vec::with_capacity(n_sides);
            verts.push(p2);
            for k in 1..=(n_sides - 2) {
                let ang = a_start + (k as f64) * step;
                verts.push((rcx + r * ang.cos(), rcy + r * ang.sin()));
            }
            verts.push(p1);
            verts
        };

        let overlap_tol = (blen * 0.35).max(4.0);
        let count_overlaps = |vs: &[(f64, f64)]| -> usize {
            vs.iter()
                .skip(1)
                .take(n_sides.saturating_sub(2))
                .filter(|&&(vx, vy)| mol.atoms.iter().any(|a| (a.x - vx).hypot(a.y - vy) <= overlap_tol))
                .count()
        };

        let cand = build_side(nx, ny);
        if count_overlaps(&cand) >= n_sides.saturating_sub(2) {
            let alt = build_side(-nx, -ny);
            if count_overlaps(&alt) < count_overlaps(&cand) {
                return alt;
            }
        }
        cand
    }

    /// Commits a ring placement (spiro on atom, edge-fused on bond, or free polygon on canvas)
    /// supporting both single-click and drag-to-orient / drag-to-pick-side!
    fn place_ring_interaction(
        &mut self,
        n_sides: usize,
        aromatic: bool,
        cpd: bool,
        sx: f64,
        sy: f64,
        drag_end: Option<(f64, f64)>,
    ) {
        let blen = self.doc().settings.bond_length_pt;
        let atol = self.atom_hit_tol();
        let btol = self.bond_hit_tol();

        // 1. Atom-hit -> Spiro Fusion (oriented toward drag_end if dragged)
        if let Some((mi, ai)) = self.hit_atom(sx, sy, atol) {
            self.doc_mut().apply_scene_tx("Spiro Ring", |scene, next_id| {
                if let Some(mol) = scene.molecules.get_mut(mi) {
                    let existing_deg = mol
                        .bonds
                        .iter()
                        .filter(|eb| eb.begin_atom == ai || eb.end_atom == ai)
                        .count();
                    let verts = Self::compute_spiro_ring_vertices(mol, ai, n_sides, blen, drag_end);
                    let mut prev_idx = ai;
                    for (k, &(vx, vy)) in verts.iter().enumerate().skip(1) {
                        let nid = *next_id;
                        *next_id += 1;
                        let mut na = ChemAtom::new(nid, 6);
                        na.x = vx;
                        na.y = vy;
                        let new_idx = mol.atoms.len();
                        mol.atoms.push(na);
                        let bid = *next_id;
                        *next_id += 1;
                        let ord = if aromatic {
                            // If spiro carbon already has >=2 bonds, avoid placing a double bond on `ai`
                            if existing_deg >= 2 {
                                if k % 2 == 0 { 2 } else { 1 }
                            } else if k % 2 == 1 {
                                2
                            } else {
                                1
                            }
                        } else if cpd && (k == 2 || k == 4) {
                            2
                        } else {
                            1
                        };
                        mol.bonds.push(ChemBond {
                            id: bid,
                            begin_atom: prev_idx,
                            end_atom: new_idx,
                            order: ord,
                            display: BondDisplay::Solid,
                            visible: true,
                            ignore_warnings: false,
                        });
                        prev_idx = new_idx;
                    }
                    let bid = *next_id;
                    *next_id += 1;
                    mol.bonds.push(ChemBond {
                        id: bid,
                        begin_atom: prev_idx,
                        end_atom: ai,
                        order: 1,
                        display: BondDisplay::Solid,
                        visible: true,
                        ignore_warnings: false,
                    });
                }
            });
            self.status_msg = format!("Fused spiro {n_sides}-membered ring");
            return;
        }

        // 2. Bond-hit -> Edge Ring Fusion (on the side of drag_end if dragged!)
        if let Some((mi, bi)) = self.hit_bond(sx, sy, btol) {
            self.doc_mut().apply_scene_tx("Fuse Ring", |scene, next_id| {
                if let Some(mol) = scene.molecules.get_mut(mi) {
                    let b = mol.bonds[bi].clone();
                    let shared_or_adj_double = b.order >= 2
                        || mol.bonds.iter().any(|ob| {
                            ob.order >= 2
                                && (ob.begin_atom == b.begin_atom
                                    || ob.end_atom == b.begin_atom
                                    || ob.begin_atom == b.end_atom
                                    || ob.end_atom == b.end_atom)
                        });
                    let verts = Self::compute_edge_fused_ring_vertices(mol, bi, n_sides, blen, drag_end);
                    let merge_tol = (blen * 0.32).max(4.0);
                    let mut prev_idx = b.end_atom;
                    for (k, &(vx, vy)) in verts.iter().enumerate().skip(1).take(n_sides - 2) {
                        // Auto-weld to an existing atom if fusing into a bay (e.g. phenanthrene/pyrene)!
                        let existing = mol.atoms.iter().enumerate().find_map(|(idx, a)| {
                            if idx != prev_idx && (a.x - vx).hypot(a.y - vy) <= merge_tol {
                                Some(idx)
                            } else {
                                None
                            }
                        });
                        let target_idx = if let Some(ei) = existing {
                            ei
                        } else {
                            let nid = *next_id;
                            *next_id += 1;
                            let mut na = ChemAtom::new(nid, 6);
                            na.x = vx;
                            na.y = vy;
                            let new_idx = mol.atoms.len();
                            mol.atoms.push(na);
                            new_idx
                        };
                        let already_bonded = mol.bonds.iter().any(|eb| {
                            (eb.begin_atom == prev_idx && eb.end_atom == target_idx)
                                || (eb.begin_atom == target_idx && eb.end_atom == prev_idx)
                        });
                        if !already_bonded {
                            let bid = *next_id;
                            *next_id += 1;
                            let ord = if aromatic {
                                if shared_or_adj_double {
                                    if k % 2 == 0 { 2 } else { 1 }
                                } else if k % 2 == 1 {
                                    2
                                } else {
                                    1
                                }
                            } else if cpd && (k == 1 || k == 3) {
                                2
                            } else {
                                1
                            };
                            mol.bonds.push(ChemBond {
                                id: bid,
                                begin_atom: prev_idx,
                                end_atom: target_idx,
                                order: ord,
                                display: BondDisplay::Solid,
                                visible: true,
                                ignore_warnings: false,
                            });
                        }
                        prev_idx = target_idx;
                    }
                    let already_closing = mol.bonds.iter().any(|eb| {
                        (eb.begin_atom == prev_idx && eb.end_atom == b.begin_atom)
                            || (eb.begin_atom == b.begin_atom && eb.end_atom == prev_idx)
                    });
                    if !already_closing && prev_idx != b.begin_atom {
                        let bid = *next_id;
                        *next_id += 1;
                        let closing_ord = if aromatic && !shared_or_adj_double && n_sides % 2 == 0 {
                            2
                        } else {
                            1
                        };
                        mol.bonds.push(ChemBond {
                            id: bid,
                            begin_atom: prev_idx,
                            end_atom: b.begin_atom,
                            order: closing_ord,
                            display: BondDisplay::Solid,
                            visible: true,
                            ignore_warnings: false,
                        });
                    }
                }
            });
            self.status_msg = format!("Fused {n_sides}-membered ring across bond");
            return;
        }

        // 3. Empty canvas -> Free ring (rotated to drag angle if dragged!)
        let rot_ang = drag_end
            .filter(|&(ex, ey)| (ex - sx).hypot(ey - sy) > 4.0)
            .map(|(ex, ey)| (ey - sy).atan2(ex - sx))
            .unwrap_or(-std::f64::consts::FRAC_PI_2);
        let ring = Self::build_polygon_ring_angled(n_sides, aromatic, cpd, sx, sy, blen, rot_ang);
        self.doc_mut().stamp_molecule_tx(ring, "Ring");
    }

    fn build_acyclic_chain(sx: f64, sy: f64, ex: f64, ey: f64, blen: f64) -> ChemMolecule {
        let dx = ex - sx;
        let dy = ey - sy;
        let dist = dx.hypot(dy);
        let step_x = blen * 0.8660254;
        let step_y = blen * 0.5;
        let n_bonds = ((dist / step_x).round() as usize).clamp(1, 30);
        let base_ang = dy.atan2(dx);
        let (sin_a, cos_a) = base_ang.sin_cos();

        let mut m = ChemMolecule::default();
        for i in 0..=n_bonds {
            let lx = (i as f64) * step_x;
            let ly = if i % 2 == 0 { 0.0 } else { -step_y };
            let mut a = ChemAtom::new((i + 1) as u32, 6);
            a.x = sx + lx * cos_a - ly * sin_a;
            a.y = sy + lx * sin_a + ly * cos_a;
            m.atoms.push(a);
            if i > 0 {
                m.bonds.push(ChemBond {
                    id: (n_bonds + 2 + i) as u32,
                    begin_atom: i - 1,
                    end_atom: i,
                    order: 1,
                    display: BondDisplay::Solid,
                    visible: true,
                    ignore_warnings: false,
                });
            }
        }
        m
    }

    fn open_inline_atom_edit(&mut self, mi: usize, ai: usize) {
        if let Some(a) = self.doc().scene.molecules.get(mi).and_then(|m| m.atoms.get(ai)) {
            self.inline_label_buf = a
                .label
                .clone()
                .unwrap_or_else(|| get_element(a.z).map(|e| e.symbol).unwrap_or("C").to_string());
            self.inline_atom_edit = Some((mi, ai));
            self.inline_just_opened = true;
        }
    }

    fn commit_inline_atom_edit(&mut self) {
        let Some((mi, ai)) = self.inline_atom_edit.take() else {
            return;
        };
        let txt = self.inline_label_buf.trim().to_string();
        if txt.is_empty() {
            return;
        }
        // Smart Reverse-Order Group Typing (e.g. HO, H2N, O2N, MeO, AcO, F3C, BocHN, HO2C)
        let rev_match = resolve_reverse_abbreviation(&txt);
        self.doc_mut().apply_scene_tx("Label Change", |scene, _| {
            if let Some(a) = scene.molecules.get_mut(mi).and_then(|m| m.atoms.get_mut(ai)) {
                if let Some((canon, disp)) = rev_match {
                    match canon {
                        "OH" => { a.z = 8; a.label = None; }
                        "SH" => { a.z = 16; a.label = None; }
                        "NH2" => { a.z = 7; a.label = None; }
                        _ => {
                            // Set true bonding heteroatom Z for accurate connectivity while preserving left-facing display!
                            if canon.starts_with('O') {
                                a.z = 8;
                            } else if canon.starts_with('N') {
                                a.z = 7;
                            } else if canon.starts_with('S') {
                                a.z = 16;
                            } else if canon == "TMS" {
                                a.z = 14;
                            } else {
                                a.z = 6;
                            }
                            a.label = Some(disp.to_string());
                        }
                    }
                } else if let Some(el) = get_element_by_symbol(&txt) {
                    a.z = el.z;
                    a.label = None;
                } else {
                    a.label = Some(txt.clone());
                }
            }
        });
        if let Some((canon, disp)) = rev_match {
            self.status_msg = format!("Smart reverse group \"{disp}\" -> {canon}");
        } else {
            self.status_msg = format!("Set atom label \"{txt}\"");
        }
    }

    fn toggle_aromatic_display_tx(&mut self) {
        self.doc_mut().apply_scene_tx("Toggle Aromatic Display", |scene, _| {
            if !scene.aromatic_circles.is_empty() {
                scene.aromatic_circles.clear();
                return;
            }
            for (mi, mol) in scene.molecules.iter().enumerate() {
                let arom = CcAromaticityInfo::perceive(mol);
                let arom_atoms: Vec<usize> = (0..mol.atoms.len())
                    .filter(|&i| arom.is_atom_aromatic(i))
                    .collect();
                if arom_atoms.len() >= 5 {
                    let cx = arom_atoms.iter().map(|&i| mol.atoms[i].x).sum::<f64>()
                        / (arom_atoms.len() as f64);
                    let cy = arom_atoms.iter().map(|&i| mol.atoms[i].y).sum::<f64>()
                        / (arom_atoms.len() as f64);
                    let avg_r = arom_atoms
                        .iter()
                        .map(|&i| (mol.atoms[i].x - cx).hypot(mol.atoms[i].y - cy))
                        .sum::<f64>()
                        / (arom_atoms.len() as f64);
                    scene.aromatic_circles.push(AromaticCircle {
                        molecule_idx: Some(mi),
                        center_x: cx,
                        center_y: cy,
                        radius: avg_r * 0.625,
                    });
                }
            }
        });
    }

    fn paste_caption_tx(&mut self, text: String) {
        let (cx, cy) = self
            .selection_bbox()
            .map(|(x1, _, _, y2)| (x1, y2 + 22.0))
            .unwrap_or((120.0, 260.0));
        let line_step = 13.0 * (self.line_spacing_pct as f64 / 100.0).clamp(0.5, 3.0);
        self.doc_mut().apply_scene_tx("Paste Caption", |scene, next_id| {
            for (i, line) in text.lines().enumerate() {
                let id = *next_id;
                *next_id += 1;
                scene.captions.push(CaptionNode {
                    id,
                    text: line.to_string(),
                    x: cx,
                    y: cy + (i as f64) * line_step,
                    element: None,
                    generic_nickname: None,
                    warning: None,
                    justification: LabelJustification::Left,
                });
            }
        });
        self.status_msg = "Pasted caption onto canvas".to_string();
    }

    fn handle_canvas_click(&mut self, dx: f64, dy: f64, shift: bool) {
        let blen = self.doc().settings.bond_length_pt;
        let atol = self.atom_hit_tol();
        let btol = self.bond_hit_tol();
        if let Some(el_z) = self.active_element_z {
            if let Some((mi, ai)) = self.hit_atom(dx, dy, atol) {
                self.doc_mut().apply_scene_tx("Change Element", |scene, _| {
                    if let Some(a) = scene.molecules.get_mut(mi).and_then(|m| m.atoms.get_mut(ai)) {
                        a.z = el_z;
                        a.label = None;
                    }
                });
                self.status_msg = format!(
                    "Set atom to {}",
                    get_element(el_z).map(|e| e.symbol).unwrap_or("?")
                );
                return;
            }
        }

        if matches!(self.tool, Tool::Marquee | Tool::Lasso) {
            if let Some(mi) = self.hit_molecule_bbox(dx, dy) {
                if shift {
                    if let Some(pos) = self.tab().selected_mols.iter().position(|&x| x == mi) {
                        self.tab_mut().selected_mols.remove(pos);
                    } else {
                        self.tab_mut().selected_mols.push(mi);
                    }
                } else {
                    self.tab_mut().selected_mols = vec![mi];
                }
                self.expand_selection_with_groups();
                self.status_msg = format!("Selected {} structure(s)", self.tab().selected_mols.len());
            } else if !shift {
                self.tab_mut().selected_mols.clear();
            }
            return;
        }

        if matches!(self.tool, Tool::Chair1 | Tool::Chair2) {
            let mut m = stamp_template(BuiltinTemplate::CyclohexaneChair, dx, dy, blen);
            if self.tool == Tool::Chair2 {
                let all: Vec<usize> = (0..m.atoms.len()).collect();
                flip_atoms_preserving_stereo(&mut m, &all, FlipAxis::Horizontal, false);
            }
            self.doc_mut().stamp_molecule_tx(m, "Cyclohexane Chair");
            self.status_msg = "Stamped Cyclohexane Chair".to_string();
            return;
        }

        if let Some((n_sides, aromatic, cpd)) = Self::ring_spec(self.tool) {
            self.place_ring_interaction(n_sides, aromatic, cpd, dx, dy, None);
            return;
        }

        match self.tool {
            Tool::TextLabel => {
                if let Some((mi, ai)) = self.hit_atom(dx, dy, atol) {
                    self.open_inline_atom_edit(mi, ai);
                } else if let Some(ci) = self.hit_caption(dx, dy) {
                    self.label_input = self.doc().scene.captions[ci].text.clone();
                    self.dialog = ActiveDialog::EditCaption { cap_idx: ci };
                } else {
                    self.label_input.clear();
                    self.dialog = ActiveDialog::NewCaption {
                        x_tenths: (dx * 10.0).round() as i32,
                        y_tenths: (dy * 10.0).round() as i32,
                    };
                }
            }
            Tool::Orbital | Tool::OrbitalS | Tool::OrbitalD => {
                let st = match self.tool {
                    Tool::OrbitalS => "s",
                    Tool::OrbitalD => "d",
                    _ => "p",
                };
                self.doc_mut().apply_scene_tx("Add Orbital", |scene, next_id| {
                    let id = *next_id;
                    *next_id += 1;
                    scene.graphics.push(SceneGraphic {
                        id,
                        graphic_type: "Orbital".to_string(),
                        sub_type: Some(st.to_string()),
                        bbox: (dx - 6.0, dy - 14.0, dx + 6.0, dy + 14.0),
                    });
                });
            }
            Tool::ChemSymbol | Tool::SymbolMinus | Tool::SymbolRadical | Tool::SymbolLonePair => {
                if let Some((mi, ai)) = self.hit_atom(dx, dy, 6.0) {
                    let t = self.tool;
                    self.doc_mut().apply_scene_tx("Atom Symbol", |scene, _| {
                        if let Some(a) = scene.molecules.get_mut(mi).and_then(|m| m.atoms.get_mut(ai)) {
                            match t {
                                Tool::ChemSymbol => a.charge += 1,
                                Tool::SymbolMinus => a.charge -= 1,
                                Tool::SymbolRadical => a.radical = if a.radical == 2 { 0 } else { 2 },
                                Tool::SymbolLonePair => a.radical = if a.radical == 1 { 0 } else { 1 },
                                _ => {}
                            }
                        }
                    });
                } else {
                    let st = match self.tool {
                        Tool::SymbolMinus => "CircleMinus",
                        Tool::SymbolRadical => "Radical",
                        Tool::SymbolLonePair => "LonePair",
                        _ => "CirclePlus",
                    };
                    self.doc_mut().apply_scene_tx("Add Symbol", |scene, next_id| {
                        let id = *next_id;
                        *next_id += 1;
                        scene.graphics.push(SceneGraphic {
                            id,
                            graphic_type: "Symbol".to_string(),
                            sub_type: Some(st.to_string()),
                            bbox: (dx - 5.0, dy - 5.0, dx + 5.0, dy + 5.0),
                        });
                    });
                }
            }
            Tool::Table => {
                self.doc_mut().apply_scene_tx("Insert Table", |scene, next_id| {
                    let id = *next_id;
                    *next_id += 1;
                    scene.tables.push(SceneTable {
                        id,
                        rows: 2,
                        cols: 2,
                        bbox: (dx, dy, dx + 96.0, dy + 62.0),
                        cells: vec![
                            vec!["Entry".to_string(), "Yield (%)".to_string()],
                            vec!["1".to_string(), "92".to_string()],
                        ],
                    });
                });
                self.show_object_inspector = true;
            }
            Tool::TlcPlate => {
                self.doc_mut().apply_scene_tx("Add TLC Plate", |scene, next_id| {
                    let id = *next_id;
                    *next_id += 1;
                    scene.tlc_plates.push(SceneTlcPlate {
                        id,
                        bbox: (dx, dy, dx + 42.0, dy + 78.0),
                        origin_fraction: 0.15,
                        spot_rfs: vec![0.25, 0.55, 0.78],
                        spot_labels: vec!["SM".to_string(), "Co".to_string(), "Rxn".to_string()],
                        spot_colors: vec![(20, 20, 20), (25, 95, 205), (205, 45, 45)],
                        solvent_label: "EtOAc/Hex 1:4".to_string(),
                    });
                });
                self.show_object_inspector = true;
            }
            Tool::Eraser => {
                if let Some((mi, ai)) = self.hit_atom(dx, dy, 6.0) {
                    let mut mol_removed = false;
                    self.doc_mut().apply_scene_tx("Erase Atom", |scene, _| {
                        if let Some(mol) = scene.molecules.get_mut(mi) {
                            mol.bonds.retain(|b| b.begin_atom != ai && b.end_atom != ai);
                            mol.atoms.remove(ai);
                            for b in &mut mol.bonds {
                                if b.begin_atom > ai { b.begin_atom -= 1; }
                                if b.end_atom > ai { b.end_atom -= 1; }
                            }
                            if mol.atoms.is_empty() {
                                scene.molecules.remove(mi);
                                mol_removed = true;
                            }
                        }
                    });
                    self.tab_mut().remap_after_atom_removal(mi, ai);
                    if mol_removed {
                        self.tab_mut().remap_after_mol_removal(&[mi]);
                    }
                    self.sanitize_state();
                } else if let Some((mi, bi)) = self.hit_bond(dx, dy, 5.0) {
                    self.doc_mut().apply_scene_tx("Erase Bond", |scene, _| {
                        if let Some(mol) = scene.molecules.get_mut(mi) {
                            if bi < mol.bonds.len() {
                                mol.bonds.remove(bi);
                            }
                        }
                    });
                    self.tab_mut().remap_after_bond_removal(mi, bi);
                } else if let Some(ci) = self.hit_caption(dx, dy) {
                    self.doc_mut().apply_scene_tx("Erase Caption", |scene, _| {
                        if ci < scene.captions.len() {
                            scene.captions.remove(ci);
                        }
                    });
                } else {
                    // Only push an Undo transaction when an actual non-molecule object is hit!
                    let hit_arrow = self.doc().scene.arrows.iter().position(|a| {
                        let (vx, vy) = (a.head.0 - a.tail.0, a.head.1 - a.tail.1);
                        let l2 = (vx * vx + vy * vy).max(1e-6);
                        let t = (((dx - a.tail.0) * vx + (dy - a.tail.1) * vy) / l2).clamp(0.0, 1.0);
                        (a.tail.0 + t * vx - dx).hypot(a.tail.1 + t * vy - dy) <= 7.0
                    });
                    let hit_curve = self.doc().scene.curves.iter().position(|c| {
                        c.points.iter().any(|&(px, py)| (px - dx).hypot(py - dy) <= 10.0)
                    });
                    let hit_graphic = self.doc().scene.graphics.iter().position(|g| {
                        let (x1, y1, x2, y2) = (g.bbox.0.min(g.bbox.2), g.bbox.1.min(g.bbox.3), g.bbox.0.max(g.bbox.2), g.bbox.1.max(g.bbox.3));
                        dx >= x1 - 5.0 && dx <= x2 + 5.0 && dy >= y1 - 5.0 && dy <= y2 + 5.0
                    });
                    let hit_table = self.doc().scene.tables.iter().position(|t| {
                        dx >= t.bbox.0 - 4.0 && dx <= t.bbox.2 + 4.0 && dy >= t.bbox.1 - 4.0 && dy <= t.bbox.3 + 4.0
                    });
                    let hit_tlc = self.doc().scene.tlc_plates.iter().position(|t| {
                        dx >= t.bbox.0 - 4.0 && dx <= t.bbox.2 + 4.0 && dy >= t.bbox.1 - 4.0 && dy <= t.bbox.3 + 4.0
                    });
                    let hit_ep = self.doc().scene.energy_profiles.iter().position(|ep| {
                        let (x1, y1, x2, y2) = (ep.bbox.0.min(ep.bbox.2), ep.bbox.1.min(ep.bbox.3), ep.bbox.0.max(ep.bbox.2), ep.bbox.1.max(ep.bbox.3));
                        dx >= x1 - 4.0 && dx <= x2 + 4.0 && dy >= y1 - 4.0 && dy <= y2 + 4.0
                    });
                    let hit_gel = self.tab().gel_plates.iter().position(|gp| {
                        let (x1, y1, x2, y2) = (gp.bbox.0.min(gp.bbox.2), gp.bbox.1.min(gp.bbox.3), gp.bbox.0.max(gp.bbox.2), gp.bbox.1.max(gp.bbox.3));
                        dx >= x1 - 4.0 && dx <= x2 + 4.0 && dy >= y1 - 4.0 && dy <= y2 + 4.0
                    });
                    if let Some(ai) = hit_arrow {
                        self.doc_mut().apply_scene_tx("Erase Arrow", |s, _| { if ai < s.arrows.len() { s.arrows.remove(ai); } });
                    } else if let Some(ci) = hit_curve {
                        self.doc_mut().apply_scene_tx("Erase Curve", |s, _| { if ci < s.curves.len() { s.curves.remove(ci); } });
                    } else if let Some(gi) = hit_graphic {
                        self.doc_mut().apply_scene_tx("Erase Graphic", |s, _| { if gi < s.graphics.len() { s.graphics.remove(gi); } });
                    } else if let Some(ti) = hit_table {
                        self.doc_mut().apply_scene_tx("Erase Table", |s, _| { if ti < s.tables.len() { s.tables.remove(ti); } });
                    } else if let Some(ti) = hit_tlc {
                        self.doc_mut().apply_scene_tx("Erase TLC Plate", |s, _| { if ti < s.tlc_plates.len() { s.tlc_plates.remove(ti); } });
                    } else if let Some(ei) = hit_ep {
                        self.doc_mut().apply_scene_tx("Erase Energy Profile", |s, _| { if ei < s.energy_profiles.len() { s.energy_profiles.remove(ei); } });
                    } else if let Some(gi) = hit_gel {
                        self.tab_mut().gel_plates.remove(gi);
                        self.status_msg = "Erased Gel Electrophoresis Plate".to_string();
                    }
                }
            }
            Tool::Templates => {
                if let Some(lib) = self.ctp_libraries.get(self.selected_ctp_idx).cloned() {
                    let lib_name = lib.name.clone();
                    if let Some(m) = lib.stamp(0, dx, dy) {
                        self.doc_mut().stamp_molecule_tx(m, "Template Stamp");
                        self.status_msg = format!("Stamped {lib_name} #1");
                    } else if let Some(curves) = lib.stamp_graphic(0, dx, dy) {
                        self.doc_mut().apply_scene_tx("Graphic Template Stamp", |s, nid| {
                            for mut c in curves {
                                c.id = *nid;
                                *nid += 1;
                                s.curves.push(c);
                            }
                        });
                        self.status_msg = format!("Stamped {lib_name} graphic #1");
                    } else {
                        self.active_flyout = ActiveFlyout::Templates;
                    }
                } else {
                    self.active_flyout = ActiveFlyout::Templates;
                }
            }
            Tool::AcyclicChain => {
                if let Some((mi, ai)) = self.hit_atom(dx, dy, atol) {
                    let (ax, ay) = (self.doc().scene.molecules[mi].atoms[ai].x, self.doc().scene.molecules[mi].atoms[ai].y);
                    let ang = Self::best_new_bond_angle(&self.doc().scene.molecules[mi], ai, 1, self.doc().settings.chain_angle_deg.to_radians());
                    let chain = Self::build_acyclic_chain(ax, ay, ax + blen * 3.2 * ang.cos(), ay + blen * 3.2 * ang.sin(), blen);
                    let added = chain.atoms.len().saturating_sub(1);
                    self.doc_mut().apply_scene_tx("Attach Chain", |scene, next_id| {
                        if let Some(mol) = scene.molecules.get_mut(mi) {
                            let mut prev_idx = ai;
                            for ca in chain.atoms.iter().skip(1) {
                                let nid = *next_id;
                                *next_id += 1;
                                let mut na = ChemAtom::new(nid, 6);
                                na.x = ca.x;
                                na.y = ca.y;
                                let new_idx = mol.atoms.len();
                                mol.atoms.push(na);
                                let bid = *next_id;
                                *next_id += 1;
                                mol.bonds.push(ChemBond {
                                    id: bid,
                                    begin_atom: prev_idx,
                                    end_atom: new_idx,
                                    order: 1,
                                    display: BondDisplay::Solid,
                                    visible: true,
                                    ignore_warnings: false,
                                });
                                prev_idx = new_idx;
                            }
                        }
                    });
                    self.status_msg = format!("Attached {added}-carbon chain to atom");
                } else {
                    let chain = Self::build_acyclic_chain(dx, dy, dx + blen * 3.2, dy, blen);
                    let nc = chain.atoms.len();
                    self.doc_mut().stamp_molecule_tx(chain, "Add Chain");
                    self.status_msg = format!("Added {nc}-carbon acyclic chain");
                }
            }
            Tool::PenCurve => {
                self.doc_mut().apply_scene_tx("Add Curve", |scene, next_id| {
                    let id = *next_id;
                    *next_id += 1;
                    scene.curves.push(SceneCurve {
                        id,
                        points: vec![
                            (dx - 24.0, dy),
                            (dx, dy - 22.0),
                            (dx, dy - 22.0),
                            (dx + 24.0, dy),
                        ],
                        closed: false,
                    });
                });
            }
            Tool::ShapeRect | Tool::ShapeOval | Tool::ShapeLine | Tool::Bracket | Tool::RoundBracket | Tool::CurlyBracket => {
                let (gtype, stype, bb) = match self.tool {
                    Tool::Bracket => ("Bracket", Some("Square".to_string()), (dx - 28.0, dy - 22.0, dx + 28.0, dy + 22.0)),
                    Tool::RoundBracket => ("Bracket", Some("Round".to_string()), (dx - 28.0, dy - 22.0, dx + 28.0, dy + 22.0)),
                    Tool::CurlyBracket => ("Bracket", Some("Curly".to_string()), (dx - 28.0, dy - 22.0, dx + 28.0, dy + 22.0)),
                    Tool::ShapeOval => ("Oval", Some("Circle".to_string()), (dx - 22.0, dy - 22.0, dx + 22.0, dy + 22.0)),
                    Tool::ShapeLine => ("Line", None, (dx - 26.0, dy, dx + 26.0, dy)),
                    _ => ("Rectangle", Some("Shadow".to_string()), (dx - 32.0, dy - 22.0, dx + 32.0, dy + 22.0)),
                };
                self.doc_mut().apply_scene_tx(gtype, |scene, next_id| {
                    let id = *next_id;
                    *next_id += 1;
                    scene.graphics.push(SceneGraphic {
                        id,
                        graphic_type: gtype.to_string(),
                        sub_type: stype,
                        bbox: bb,
                    });
                });
            }
            Tool::Arrow | Tool::NoGoArrow | Tool::EquilibriumArrow | Tool::RetrosyntheticArrow | Tool::CurvedArrow => {
                let (head, tail, nogo) = match self.tool {
                    Tool::NoGoArrow => ("Full".to_string(), None, Some("Cross".to_string())),
                    Tool::EquilibriumArrow => ("HalfLeft".to_string(), Some("HalfLeft".to_string()), None),
                    Tool::RetrosyntheticArrow => ("Retrosynthetic".to_string(), None, None),
                    Tool::CurvedArrow => ("Curved".to_string(), None, None),
                    _ => ("Full".to_string(), None, None),
                };
                self.doc_mut().apply_scene_tx("Add Arrow", |scene, next_id| {
                    let id = *next_id;
                    *next_id += 1;
                    scene.arrows.push(SceneArrow {
                        id,
                        tail: (dx - 26.0, dy),
                        head: (dx + 26.0, dy),
                        arrowhead_head: head,
                        arrowhead_tail: tail,
                        no_go: nogo,
                        dipole: false,
                    });
                });
            }
            t if self.is_bond_tool() => {
                let (ord, disp) = self.current_bond_spec();
                // ChiralCraft UX: atom vertices ALWAYS take hit priority over bond midpoints!
                if let Some((mi, ai)) = self.hit_atom(dx, dy, atol) {
                    let chain_rad = self.doc().settings.chain_angle_deg.to_radians();
                    self.doc_mut().apply_scene_tx("Bond", |scene, next_id| {
                        if let Some(mol) = scene.molecules.get_mut(mi) {
                            let ax = mol.atoms[ai].x;
                            let ay = mol.atoms[ai].y;
                            let ang = Self::best_new_bond_angle(mol, ai, ord, chain_rad);
                            let nx = ax + blen * ang.cos();
                            let ny = ay + blen * ang.sin();
                            let merge_tol = (blen * 0.35).max(4.5);
                            let existing_target = mol.atoms.iter().enumerate().find_map(|(idx, a)| {
                                if idx != ai && (a.x - nx).hypot(a.y - ny) <= merge_tol {
                                    Some(idx)
                                } else {
                                    None
                                }
                            });
                            let target_idx = if let Some(ei) = existing_target {
                                if mol.bonds.iter().any(|b| {
                                    (b.begin_atom == ai && b.end_atom == ei)
                                        || (b.begin_atom == ei && b.end_atom == ai)
                                }) {
                                    return;
                                }
                                ei
                            } else {
                                let nid = *next_id;
                                *next_id += 1;
                                let mut na = ChemAtom::new(nid, 6);
                                na.x = nx;
                                na.y = ny;
                                let new_idx = mol.atoms.len();
                                mol.atoms.push(na);
                                new_idx
                            };
                            let bid = *next_id;
                            *next_id += 1;
                            mol.bonds.push(ChemBond {
                                id: bid,
                                begin_atom: ai,
                                end_atom: target_idx,
                                order: ord,
                                display: disp,
                                visible: true,
                                ignore_warnings: false,
                            });
                        }
                    });
                } else if let Some((mi, bi)) = self.hit_bond(dx, dy, btol) {
                    let cur_order = self.doc().scene.molecules.get(mi).and_then(|m| m.bonds.get(bi)).map(|b| b.order).unwrap_or(1);
                    if matches!(t, Tool::DoubleBond | Tool::MultipleBond) && cur_order == 2 {
                        // ChemDraw UX: clicking an existing double bond with Double Bond cycles its alignment!
                        let cur = self.tab().double_bond_align.get(&(mi, bi)).copied().unwrap_or(DoubleBondAlign::Auto);
                        let next = match cur {
                            DoubleBondAlign::Auto => DoubleBondAlign::Center,
                            DoubleBondAlign::Center => DoubleBondAlign::Left,
                            DoubleBondAlign::Left => DoubleBondAlign::Right,
                            DoubleBondAlign::Right => DoubleBondAlign::Auto,
                        };
                        self.tab_mut().double_bond_align.insert((mi, bi), next);
                        self.status_msg = format!("Cycled double-bond alignment ({next:?})");
                    } else {
                        self.doc_mut().apply_scene_tx("Change Bond", |scene, _| {
                            if let Some(b) = scene.molecules.get_mut(mi).and_then(|m| m.bonds.get_mut(bi)) {
                                match t {
                                    Tool::SolidBond => {
                                        b.order = match b.order {
                                            1 => 2,
                                            2 => 3,
                                            _ => 1,
                                        };
                                        b.display = BondDisplay::Solid;
                                    }
                                    Tool::WedgeBond => {
                                        b.order = 1;
                                        b.display = if b.display == BondDisplay::WedgeBegin {
                                            BondDisplay::WedgeEnd
                                        } else {
                                            BondDisplay::WedgeBegin
                                        };
                                    }
                                    Tool::HashWedgeBond => {
                                        b.order = 1;
                                        b.display = if b.display == BondDisplay::WedgedHashBegin {
                                            BondDisplay::WedgedHashEnd
                                        } else {
                                            BondDisplay::WedgedHashBegin
                                        };
                                    }
                                    Tool::HollowWedgeBond => {
                                        b.order = 1;
                                        b.display = if b.display == BondDisplay::HollowWedgeBegin {
                                            BondDisplay::HollowWedgeEnd
                                        } else {
                                            BondDisplay::HollowWedgeBegin
                                        };
                                    }
                                    Tool::DativeBond => {
                                        if b.display == BondDisplay::Dative {
                                            std::mem::swap(&mut b.begin_atom, &mut b.end_atom);
                                        } else {
                                            b.order = 1;
                                            b.display = BondDisplay::Dative;
                                        }
                                    }
                                    _ => {
                                        b.order = ord;
                                        b.display = disp;
                                    }
                                }
                            }
                        });
                    }
                } else {
                    let ang = -std::f64::consts::FRAC_PI_6;
                    let mut m = ChemMolecule::default();
                    let mut a0 = ChemAtom::new(1, 6);
                    a0.x = dx;
                    a0.y = dy;
                    let mut a1 = ChemAtom::new(2, 6);
                    a1.x = dx + blen * ang.cos();
                    a1.y = dy + blen * ang.sin();
                    m.atoms.push(a0);
                    m.atoms.push(a1);
                    m.bonds.push(ChemBond {
                        id: 3,
                        begin_atom: 0,
                        end_atom: 1,
                        order: ord,
                        display: disp,
                        visible: true,
                        ignore_warnings: false,
                    });
                    self.doc_mut().stamp_molecule_tx(m, "Bond");
                }
            }
            _ => {}
        }
    }

    fn handle_canvas_drag_release(&mut self, sx: f64, sy: f64, ex: f64, ey: f64, shift: bool) {
        if (ex - sx).hypot(ey - sy) < 4.0 {
            self.handle_canvas_click(sx, sy, shift);
            return;
        }
        let blen = self.doc().settings.bond_length_pt;
        let atol = self.atom_hit_tol();
        let btol = self.bond_hit_tol();

        // 0. Ring & Chair tools support drag-to-orient / drag-to-pick-fusion-side (never silently drop!)
        if let Some((n_sides, aromatic, cpd)) = Self::ring_spec(self.tool) {
            self.place_ring_interaction(n_sides, aromatic, cpd, sx, sy, Some((ex, ey)));
            return;
        }
        if matches!(self.tool, Tool::Chair1 | Tool::Chair2) {
            let mut m = stamp_template(BuiltinTemplate::CyclohexaneChair, sx, sy, blen);
            if self.tool == Tool::Chair2 {
                let all: Vec<usize> = (0..m.atoms.len()).collect();
                flip_atoms_preserving_stereo(&mut m, &all, FlipAxis::Horizontal, false);
            }
            self.doc_mut().stamp_molecule_tx(m, "Cyclohexane Chair");
            self.status_msg = "Stamped Cyclohexane Chair".to_string();
            return;
        }

        match self.tool {
            Tool::Marquee | Tool::Lasso => {
                // Single-atom vertex drag (with auto-weld if dropped onto another atom in the same molecule!)
                if let Some((mi, ai)) = self.dragging_single_atom {
                    let drop_hit = self.hit_atom(ex, ey, atol);
                    let mut welded = false;
                    let mut dropped_bonds: Vec<usize> = Vec::new();
                    self.doc_mut().apply_scene_tx("Move Atom", |scene, _| {
                        if let Some(mol) = scene.molecules.get_mut(mi) {
                            if ai < mol.atoms.len() {
                                if let Some((tmi, tai)) = drop_hit {
                                    if tmi == mi && tai != ai && tai < mol.atoms.len() {
                                        // Weld atom `ai` into `tai`!
                                        for b in &mut mol.bonds {
                                            if b.begin_atom == ai { b.begin_atom = tai; }
                                            if b.end_atom == ai { b.end_atom = tai; }
                                        }
                                        // Remove self-loops and duplicate bonds
                                        let mut seen_pairs = std::collections::BTreeSet::new();
                                        for (bi, b) in mol.bonds.iter().enumerate() {
                                            let key = (b.begin_atom.min(b.end_atom), b.begin_atom.max(b.end_atom));
                                            if key.0 == key.1 || !seen_pairs.insert(key) {
                                                dropped_bonds.push(bi);
                                            }
                                        }
                                        for &bi in dropped_bonds.iter().rev() {
                                            mol.bonds.remove(bi);
                                        }
                                        mol.atoms.remove(ai);
                                        for b in &mut mol.bonds {
                                            if b.begin_atom > ai { b.begin_atom -= 1; }
                                            if b.end_atom > ai { b.end_atom -= 1; }
                                        }
                                        welded = true;
                                        return;
                                    }
                                }
                                mol.atoms[ai].x = ex;
                                mol.atoms[ai].y = ey;
                            }
                        }
                    });
                    if welded {
                        for &bi in dropped_bonds.iter().rev() {
                            self.tab_mut().remap_after_bond_removal(mi, bi);
                        }
                        self.tab_mut().remap_after_atom_removal(mi, ai);
                    }
                    self.sanitize_state();
                    self.status_msg = "Moved atom vertex".to_string();
                    return;
                }
                if let Some(ci) = self.dragging_caption {
                    let dx = ex - sx;
                    let dy = ey - sy;
                    self.doc_mut().apply_scene_tx("Move Caption", |scene, _| {
                        if let Some(cap) = scene.captions.get_mut(ci) {
                            cap.x += dx;
                            cap.y += dy;
                        }
                    });
                    return;
                }
                if self.dragging_selected_move && !self.tab().selected_mols.is_empty() {
                    let dx = ex - sx;
                    let dy = ey - sy;
                    let sel = self.tab().selected_mols.clone();
                    if self.drag_ctrl_copy {
                        self.doc_mut().apply_scene_tx("Duplicate", |scene, next_id| {
                            for mi in sel {
                                if let Some(mut clone) = scene.molecules.get(mi).cloned() {
                                    for a in &mut clone.atoms {
                                        a.id = *next_id;
                                        *next_id += 1;
                                        a.x += dx;
                                        a.y += dy;
                                    }
                                    for b in &mut clone.bonds {
                                        b.id = *next_id;
                                        *next_id += 1;
                                    }
                                    scene.molecules.push(clone);
                                }
                            }
                        });
                    } else {
                        self.doc_mut().apply_scene_tx("Move", |scene, _| {
                            for mi in sel {
                                if let Some(m) = scene.molecules.get_mut(mi) {
                                    for a in &mut m.atoms {
                                        a.x += dx;
                                        a.y += dy;
                                    }
                                }
                            }
                        });
                    }
                    return;
                }
                if !shift {
                    self.tab_mut().selected_mols.clear();
                }
                if self.tool == Tool::Lasso && self.lasso_pts.len() >= 3 {
                    let poly = self.lasso_pts.clone();
                    let mut picked = Vec::new();
                    for (mi, mol) in self.doc().scene.molecules.iter().enumerate() {
                        let inside_cnt = mol.atoms.iter().filter(|a| point_in_polygon(a.x, a.y, &poly)).count();
                        if !mol.atoms.is_empty() && (inside_cnt * 2 >= mol.atoms.len() || inside_cnt > 0 && inside_cnt * 2 > mol.atoms.len()) {
                            picked.push(mi);
                        }
                    }
                    for mi in picked {
                        if !self.tab().selected_mols.contains(&mi) {
                            self.tab_mut().selected_mols.push(mi);
                        }
                    }
                } else {
                    // ChiralCraft P05 >50% bounding-box overlap selection
                    let rx1 = sx.min(ex);
                    let ry1 = sy.min(ey);
                    let rx2 = sx.max(ex);
                    let ry2 = sy.max(ey);
                    let mut picked = Vec::new();
                    for (mi, mol) in self.doc().scene.molecules.iter().enumerate() {
                        if mol.atoms.is_empty() {
                            continue;
                        }
                        let (mx1, my1, mx2, my2) = mol_bbox(mol);
                        // Overlap fraction per axis; a degenerate (zero-span) axis — single atom,
                        // vertical/horizontal chain — counts as covered when its centre is inside.
                        let frac = |a1: f64, a2: f64, r1: f64, r2: f64| {
                            if a2 - a1 < 1.0 {
                                let c = (a1 + a2) * 0.5;
                                if c >= r1 && c <= r2 { 1.0 } else { 0.0 }
                            } else {
                                (r2.min(a2) - r1.max(a1)).max(0.0) / (a2 - a1)
                            }
                        };
                        if frac(mx1, mx2, rx1, rx2) > 0.5 && frac(my1, my2, ry1, ry2) > 0.5 {
                            picked.push(mi);
                        }
                    }
                    for mi in picked {
                        if !self.tab().selected_mols.contains(&mi) {
                            self.tab_mut().selected_mols.push(mi);
                        }
                    }
                }
                self.expand_selection_with_groups();
                self.status_msg = format!("Selected {} structure(s)", self.tab().selected_mols.len());
            }
            Tool::Perspective => {
                let targets: Vec<usize> = if self.tab().selected_mols.is_empty() {
                    (0..self.doc().scene.molecules.len()).collect()
                } else {
                    self.tab().selected_mols.clone()
                };
                self.doc_mut().apply_scene_tx("Structure Perspective", |scene, _| {
                    for mi in targets {
                        if let Some(m) = scene.molecules.get_mut(mi) {
                            apply_perspective_tilt(m, ex - sx, ey - sy);
                        }
                    }
                });
                self.status_msg = "Applied 3D Structure Perspective tilt".to_string();
            }
            Tool::Fragmentation => {
                let mut results = Vec::new();
                for mol in &self.doc().scene.molecules {
                    results.extend(fragment_molecule_by_line(mol, (sx, sy), (ex, ey)));
                }
                if results.is_empty() {
                    self.status_msg = "Mass Fragmentation: drag red line across one or more bonds".to_string();
                } else {
                    let count = results.len();
                    self.doc_mut().apply_scene_tx("Mass Fragmentation", |scene, next_id| {
                        let cid = *next_id;
                        *next_id += 1;
                        scene.graphics.push(SceneGraphic {
                            id: cid,
                            graphic_type: "CleavageLine".to_string(),
                            sub_type: None,
                            bbox: (sx, sy, ex, ey),
                        });
                        for (formula, em, (fx, fy)) in results {
                            let id = *next_id;
                            *next_id += 1;
                            scene.captions.push(CaptionNode {
                                id,
                                text: format!("Formula: {formula}  Exact Mass: {em:.2}"),
                                x: fx - 24.0,
                                y: fy - 14.0,
                                element: None,
                                generic_nickname: None,
                                warning: None,
                                justification: LabelJustification::Left,
                            });
                        }
                    });
                    self.status_msg = format!("Mass Fragmentation: cleaved into {count} fragment(s)");
                }
            }
            Tool::ReactionMap => {
                if let (Some(a1), Some(a2)) = (self.drag_start_atom, self.hit_atom(ex, ey, 7.0)) {
                    if a1 != a2 {
                        let id = self.tab().next_rxn_map_id;
                        self.tab_mut().next_rxn_map_id += 1;
                        self.tab_mut().reaction_atom_maps.insert(a1, id);
                        self.tab_mut().reaction_atom_maps.insert(a2, id);
                        self.status_msg = format!("Mapped reaction atom pair #{id}");
                    }
                }
            }
            Tool::Eraser => {
                // Drag-to-erase: erase all atoms & bonds swept along the eraser stroke in one Undo transaction!
                let mut stroke = self.eraser_pts.clone();
                if stroke.is_empty() {
                    stroke.push((sx, sy));
                    stroke.push((ex, ey));
                }
                let mut hit_atoms: Vec<(usize, usize)> = Vec::new();
                let mut hit_bonds: Vec<(usize, usize)> = Vec::new();
                for &(px, py) in &stroke {
                    if let Some(a) = self.hit_atom(px, py, atol) {
                        if !hit_atoms.contains(&a) {
                            hit_atoms.push(a);
                        }
                    } else if let Some(b) = self.hit_bond(px, py, btol) {
                        if !hit_bonds.contains(&b) {
                            hit_bonds.push(b);
                        }
                    }
                }
                if !hit_atoms.is_empty() || !hit_bonds.is_empty() {
                    let erased_cnt = hit_atoms.len() + hit_bonds.len();
                    hit_atoms.sort_by(|a, b| b.cmp(a));
                    let erased_atoms = hit_atoms.clone();
                    hit_bonds.sort_by(|a, b| b.cmp(a));
                    let erased_bonds = hit_bonds.clone();
                    let mut emptied: Vec<usize> = Vec::new();
                    self.doc_mut().apply_scene_tx("Erase Stroke", |scene, _| {
                        // Erase bonds first (reverse index order per molecule)
                        hit_bonds.sort_by(|a, b| b.cmp(a));
                        for (mi, bi) in hit_bonds {
                            if let Some(mol) = scene.molecules.get_mut(mi) {
                                if bi < mol.bonds.len() {
                                    mol.bonds.remove(bi);
                                }
                            }
                        }
                        // Erase atoms in reverse index order per molecule
                        hit_atoms.sort_by(|a, b| b.cmp(a));
                        for (mi, ai) in hit_atoms {
                            if let Some(mol) = scene.molecules.get_mut(mi) {
                                if ai < mol.atoms.len() {
                                    mol.bonds.retain(|b| b.begin_atom != ai && b.end_atom != ai);
                                    mol.atoms.remove(ai);
                                    for b in &mut mol.bonds {
                                        if b.begin_atom > ai { b.begin_atom -= 1; }
                                        if b.end_atom > ai { b.end_atom -= 1; }
                                    }
                                }
                            }
                        }
                        emptied = (0..scene.molecules.len()).filter(|&i| scene.molecules[i].atoms.is_empty()).collect();
                        scene.molecules.retain(|m| !m.atoms.is_empty());
                    });
                    for &(mi, bi) in &erased_bonds {
                        self.tab_mut().remap_after_bond_removal(mi, bi);
                    }
                    for &(mi, ai) in &erased_atoms {
                        self.tab_mut().remap_after_atom_removal(mi, ai);
                    }
                    self.tab_mut().remap_after_mol_removal(&emptied);
                    self.sanitize_state();
                    self.status_msg = format!("Erased {erased_cnt} atom(s)/bond(s)");
                }
            }
            Tool::AcyclicChain => {
                let start_atom = self.drag_start_atom.or_else(|| self.hit_atom(sx, sy, atol));
                if let Some((mi, ai)) = start_atom {
                    let (ax, ay) = self
                        .doc()
                        .scene
                        .molecules
                        .get(mi)
                        .and_then(|m| m.atoms.get(ai))
                        .map(|a| (a.x, a.y))
                        .unwrap_or((sx, sy));
                    let chain = Self::build_acyclic_chain(ax, ay, ex, ey, blen);
                    let added = chain.atoms.len().saturating_sub(1);
                    self.doc_mut().apply_scene_tx("Attach Chain", |scene, next_id| {
                        if let Some(mol) = scene.molecules.get_mut(mi) {
                            let mut prev_idx = ai;
                            for ca in chain.atoms.iter().skip(1) {
                                let nid = *next_id;
                                *next_id += 1;
                                let mut na = ChemAtom::new(nid, 6);
                                na.x = ca.x;
                                na.y = ca.y;
                                let new_idx = mol.atoms.len();
                                mol.atoms.push(na);
                                let bid = *next_id;
                                *next_id += 1;
                                mol.bonds.push(ChemBond {
                                    id: bid,
                                    begin_atom: prev_idx,
                                    end_atom: new_idx,
                                    order: 1,
                                    display: BondDisplay::Solid,
                                    visible: true,
                                    ignore_warnings: false,
                                });
                                prev_idx = new_idx;
                            }
                        }
                    });
                    self.status_msg = format!("Attached {added}-carbon chain");
                } else {
                    let chain = Self::build_acyclic_chain(sx, sy, ex, ey, blen);
                    let nc = chain.atoms.len();
                    self.doc_mut().stamp_molecule_tx(chain, "Add Chain");
                    self.status_msg = format!("Added {nc}-carbon acyclic chain");
                }
            }
            Tool::PenCurve => {
                self.doc_mut().apply_scene_tx("Add Curve", |scene, next_id| {
                    let id = *next_id;
                    *next_id += 1;
                    let mx = 0.5 * (sx + ex);
                    let my = 0.5 * (sy + ey) - 18.0;
                    scene.curves.push(SceneCurve {
                        id,
                        points: vec![(sx, sy), (mx, my), (mx, my), (ex, ey)],
                        closed: false,
                    });
                });
            }
            Tool::ShapeRect | Tool::ShapeOval | Tool::ShapeLine | Tool::Bracket | Tool::RoundBracket | Tool::CurlyBracket => {
                let (gtype, stype, bb) = match self.tool {
                    Tool::Bracket => ("Bracket", Some("Square".to_string()), (sx.min(ex), sy.min(ey), sx.max(ex), sy.max(ey))),
                    Tool::RoundBracket => ("Bracket", Some("Round".to_string()), (sx.min(ex), sy.min(ey), sx.max(ex), sy.max(ey))),
                    Tool::CurlyBracket => ("Bracket", Some("Curly".to_string()), (sx.min(ex), sy.min(ey), sx.max(ex), sy.max(ey))),
                    Tool::ShapeOval => ("Oval", Some("Circle".to_string()), (sx.min(ex), sy.min(ey), sx.max(ex), sy.max(ey))),
                    Tool::ShapeLine => ("Line", None, (sx, sy, ex, ey)),
                    _ => ("Rectangle", Some("Shadow".to_string()), (sx.min(ex), sy.min(ey), sx.max(ex), sy.max(ey))),
                };
                self.doc_mut().apply_scene_tx(gtype, |scene, next_id| {
                    let id = *next_id;
                    *next_id += 1;
                    scene.graphics.push(SceneGraphic {
                        id,
                        graphic_type: gtype.to_string(),
                        sub_type: stype,
                        bbox: bb,
                    });
                });
            }
            Tool::Arrow | Tool::NoGoArrow | Tool::EquilibriumArrow | Tool::RetrosyntheticArrow | Tool::CurvedArrow => {
                let (head, tail, nogo) = match self.tool {
                    Tool::NoGoArrow => ("Full".to_string(), None, Some("Cross".to_string())),
                    Tool::EquilibriumArrow => ("HalfLeft".to_string(), Some("HalfLeft".to_string()), None),
                    Tool::RetrosyntheticArrow => ("Retrosynthetic".to_string(), None, None),
                    Tool::CurvedArrow => ("Curved".to_string(), None, None),
                    _ => ("Full".to_string(), None, None),
                };
                self.doc_mut().apply_scene_tx("Add Arrow", |scene, next_id| {
                    let id = *next_id;
                    *next_id += 1;
                    scene.arrows.push(SceneArrow {
                        id,
                        tail: (sx, sy),
                        head: (ex, ey),
                        arrowhead_head: head,
                        arrowhead_tail: tail,
                        no_go: nogo,
                        dipole: false,
                    });
                });
            }
            _ if self.is_bond_tool() => {
                let (ord, disp) = self.current_bond_spec();
                let fixed_len = self.doc().settings.fixed_lengths;
                let fixed_ang = self.doc().settings.fixed_angles;
                let tool_now = self.tool;
                if let Some((mi, ai)) = self.drag_start_atom {
                    let mut end_hit = self.hit_atom(ex, ey, atol.max(6.5));
                    let (ax, ay) = self
                        .doc()
                        .scene
                        .molecules
                        .get(mi)
                        .and_then(|m| m.atoms.get(ai))
                        .map(|a| (a.x, a.y))
                        .unwrap_or((sx, sy));
                    let (nex, ney) = if end_hit.is_some() {
                        (ex, ey)
                    } else {
                        snap_bond_endpoint(ax, ay, ex, ey, fixed_len, fixed_ang, blen)
                    };
                    // Critical ChemDraw UX: if the SNAPPED endpoint lands on an existing atom (e.g. closing a ring!),
                    // automatically fuse/bond to that existing atom instead of creating a duplicate overlapping carbon!
                    if end_hit.is_none() {
                        let merge_tol = (blen * 0.38).max(5.5);
                        if let Some(snapped_hit) = self.hit_atom(nex, ney, merge_tol) {
                            if snapped_hit != (mi, ai) {
                                end_hit = Some(snapped_hit);
                            }
                        }
                    }
                    let mut merged_other = false;
                    // (keep, drop, atom_off, bond_off) recorded for side-table re-keying after a merge.
                    let merge_plan = end_hit.and_then(|(emi, _)| {
                        let mols = &self.doc().scene.molecules;
                        if emi == mi || mi >= mols.len() || emi >= mols.len() {
                            return None;
                        }
                        let (k, d) = if mi < emi { (mi, emi) } else { (emi, mi) };
                        Some((k, d, mols[k].atoms.len(), mols[k].bonds.len()))
                    });
                    self.doc_mut().apply_scene_tx("Bond", |scene, next_id| {
                        // Case A: Connecting atoms across TWO distinct molecules (e.g. drawn + imported!) -> merge!
                        if let Some((emi, eai)) = end_hit {
                            if emi != mi && mi < scene.molecules.len() && emi < scene.molecules.len() {
                                let (keep_i, drop_i, a_keep, a_drop) = if mi < emi {
                                    (mi, emi, ai, eai)
                                } else {
                                    (emi, mi, eai, ai)
                                };
                                let dropped = scene.molecules.remove(drop_i);
                                if let Some(mol) = scene.molecules.get_mut(keep_i) {
                                    let offset = mol.atoms.len();
                                    // Fresh ids: separately imported structures reuse 1..n, and
                                    // duplicate ids inside one molecule confuse id-keyed lookups.
                                    for mut a in dropped.atoms {
                                        a.id = *next_id;
                                        *next_id += 1;
                                        mol.atoms.push(a);
                                    }
                                    for mut b in dropped.bonds {
                                        b.id = *next_id;
                                        *next_id += 1;
                                        b.begin_atom += offset;
                                        b.end_atom += offset;
                                        mol.bonds.push(b);
                                    }
                                    let (b1, b2) = if mi < emi {
                                        (a_keep, offset + a_drop)
                                    } else {
                                        (offset + a_drop, a_keep)
                                    };
                                    let bid = *next_id;
                                    *next_id += 1;
                                    mol.bonds.push(ChemBond {
                                        id: bid,
                                        begin_atom: b1,
                                        end_atom: b2,
                                        order: ord,
                                        display: disp,
                                        visible: true,
                                        ignore_warnings: false,
                                    });
                                    merged_other = true;
                                }
                                return;
                            }
                        }
                        if let Some(mol) = scene.molecules.get_mut(mi) {
                            let target_idx = if let Some((emi, eai)) = end_hit {
                                if emi == mi && eai != ai && eai < mol.atoms.len() {
                                    // If a bond already connects ai <-> eai, upgrade/update it instead of duplicating!
                                    if let Some(existing) = mol.bonds.iter_mut().find(|b| {
                                        (b.begin_atom == ai && b.end_atom == eai)
                                            || (b.begin_atom == eai && b.end_atom == ai)
                                    }) {
                                        if tool_now == Tool::SolidBond {
                                            existing.order = match existing.order {
                                                1 => 2,
                                                2 => 3,
                                                _ => 1,
                                            };
                                            existing.display = BondDisplay::Solid;
                                        } else {
                                            existing.order = ord;
                                            existing.display = disp;
                                        }
                                        return;
                                    }
                                    eai
                                } else {
                                    return;
                                }
                            } else {
                                let nid = *next_id;
                                *next_id += 1;
                                let mut na = ChemAtom::new(nid, 6);
                                na.x = nex;
                                na.y = ney;
                                mol.atoms.push(na);
                                mol.atoms.len() - 1
                            };
                            let bid = *next_id;
                            *next_id += 1;
                            mol.bonds.push(ChemBond {
                                id: bid,
                                begin_atom: ai,
                                end_atom: target_idx,
                                order: ord,
                                display: disp,
                                visible: true,
                                ignore_warnings: false,
                            });
                        }
                    });
                    if merged_other {
                        if let Some((k, d, ao, bo)) = merge_plan {
                            self.tab_mut().remap_after_mol_merge(k, d, ao, bo);
                        }
                        self.sanitize_state();
                        self.status_msg = "Connected structures with bond".to_string();
                    }
                } else if let Some((emi, eai)) = self.hit_atom(ex, ey, 6.5) {
                    // Reverse drag from empty space onto an existing atom -> attach directly to that atom!
                    let (ax, ay) = self
                        .doc()
                        .scene
                        .molecules
                        .get(emi)
                        .and_then(|m| m.atoms.get(eai))
                        .map(|a| (a.x, a.y))
                        .unwrap_or((ex, ey));
                    let (nsx, nsy) = snap_bond_endpoint(ax, ay, sx, sy, fixed_len, fixed_ang, blen);
                    self.doc_mut().apply_scene_tx("Bond", |scene, next_id| {
                        if let Some(mol) = scene.molecules.get_mut(emi) {
                            let nid = *next_id;
                            *next_id += 1;
                            let mut na = ChemAtom::new(nid, 6);
                            na.x = nsx;
                            na.y = nsy;
                            let new_idx = mol.atoms.len();
                            mol.atoms.push(na);
                            let bid = *next_id;
                            *next_id += 1;
                            mol.bonds.push(ChemBond {
                                id: bid,
                                begin_atom: new_idx,
                                end_atom: eai,
                                order: ord,
                                display: disp,
                                visible: true,
                                ignore_warnings: false,
                            });
                        }
                    });
                } else {
                    let (nex, ney) = snap_bond_endpoint(sx, sy, ex, ey, fixed_len, fixed_ang, blen);
                    let mut m = ChemMolecule::default();
                    let mut a0 = ChemAtom::new(1, 6);
                    a0.x = sx;
                    a0.y = sy;
                    let mut a1 = ChemAtom::new(2, 6);
                    a1.x = nex;
                    a1.y = ney;
                    m.atoms.push(a0);
                    m.atoms.push(a1);
                    m.bonds.push(ChemBond {
                        id: 3,
                        begin_atom: 0,
                        end_atom: 1,
                        order: ord,
                        display: disp,
                        visible: true,
                        ignore_warnings: false,
                    });
                    self.doc_mut().stamp_molecule_tx(m, "Bond");
                }
            }
            _ => {}
        }
    }

    fn join_close_atoms_tx(&mut self) {
        let mut removed_mol: Option<usize> = None;
        let mut removed_atom: Option<(usize, usize)> = None;
        self.doc_mut().apply_scene_tx("Join", |scene, _| {
            if scene.molecules.len() >= 2 {
                let mut best: Option<(usize, usize, usize, usize, f64)> = None;
                for i in 0..scene.molecules.len() {
                    for j in (i + 1)..scene.molecules.len() {
                        for (ai, a) in scene.molecules[i].atoms.iter().enumerate() {
                            for (aj, b) in scene.molecules[j].atoms.iter().enumerate() {
                                let d = (a.x - b.x).hypot(a.y - b.y);
                                if d < 8.0 && best.map(|x| d < x.4).unwrap_or(true) {
                                    best = Some((i, ai, j, aj, d));
                                }
                            }
                        }
                    }
                }
                if let Some((m1, a1, m2, a2, _)) = best {
                    let mol_b = scene.molecules.remove(m2);
                    let _ = join_molecules_at_atoms(&mut scene.molecules[m1], a1, &mol_b, a2);
                    removed_mol = Some(m2);
                    return;
                }
            }
            for (mi, mol) in scene.molecules.iter_mut().enumerate() {
                let mut pair = None;
                for i in 0..mol.atoms.len() {
                    for j in (i + 1)..mol.atoms.len() {
                        if (mol.atoms[i].x - mol.atoms[j].x).hypot(mol.atoms[i].y - mol.atoms[j].y) < 6.0 {
                            pair = Some((i, j));
                            break;
                        }
                    }
                    if pair.is_some() { break; }
                }
                if let Some((i, j)) = pair {
                    if join_atoms(mol, i, j) {
                        removed_atom = Some((mi, j));
                    }
                    break;
                }
            }
        });
        if let Some(m2) = removed_mol {
            self.tab_mut().remap_after_mol_removal(&[m2]);
        }
        if let Some((mi, ai)) = removed_atom {
            self.tab_mut().remap_after_atom_removal(mi, ai);
        }
        self.sanitize_state();
    }

    fn run_repeat_action(&mut self) {
        let Some(act) = self.last_repeat_action.clone() else {
            return;
        };
        if let Some(rest) = act.strip_prefix("Scale:") {
            if let Ok(f) = rest.parse::<f64>() {
                self.do_scale_targets(f);
            }
            return;
        }
        if let Some(rest) = act.strip_prefix("Rotate:") {
            if let Ok(deg) = rest.parse::<f64>() {
                self.do_rotate_targets(deg);
            }
            return;
        }
        if let Some(rest) = act.strip_prefix("Align:") {
            let mode = match rest {
                "Left edges" => AlignMode::Left,
                "L/R centers" => AlignMode::HorizontalCenter,
                "Right edges" => AlignMode::Right,
                "Top edges" => AlignMode::Top,
                "T/B centers" => AlignMode::VerticalCenter,
                _ => AlignMode::Bottom,
            };
            self.do_align_action(rest, mode);
            return;
        }
        match act.as_str() {
            "Clean Up Structure" => self.do_clean_up_structure(),
            "Clean Up Reaction" => {
                self.doc_mut().apply_scene_tx("Clean Up Reaction", |s, _| { clean_up_reaction(s); });
            }
            "Flip Horizontal" => self.do_flip_or_rot180("Flip Horizontal", FlipAxis::Horizontal, false),
            "Flip Vertical" => self.do_flip_or_rot180("Flip Vertical", FlipAxis::Vertical, false),
            "Rotate 180° Horizontal" => self.do_flip_or_rot180("Rotate 180° Horizontal", FlipAxis::Horizontal, true),
            "Rotate 180° Vertical" => self.do_flip_or_rot180("Rotate 180° Vertical", FlipAxis::Vertical, true),
            "Center on Page" => self.do_center_on_page_action(),
            "Distribute Horizontally" => self.do_distribute_action(DistributeAxis::Horizontal),
            "Distribute Vertically" => self.do_distribute_action(DistributeAxis::Vertical),
            "Toggle Aromatic Display" => self.toggle_aromatic_display_tx(),
            "Toggle Hydrogen Display" => self.do_toggle_hydrogens(),
            "Flatten" => self.do_flatten_action(),
            "Join" => self.join_close_atoms_tx(),
            _ => {}
        }
        self.status_msg = format!("Repeated: {act}");
    }

    fn apply_shortcuts(&mut self, ctx: &egui::Context) {
        if self.dialog != ActiveDialog::None
            || self.inline_atom_edit.is_some()
            || self.cmd_palette.open
            || ctx.wants_keyboard_input()
        {
            return;
        }
        let open_palette = ctx.input(|i| {
            i.key_pressed(egui::Key::F4)
                || (!i.modifiers.ctrl
                    && !i.modifiers.alt
                    && i.events.iter().any(|ev| matches!(ev, egui::Event::Text(s) if s == ":")))
        });
        if open_palette {
            self.cmd_palette.open_fresh();
            return;
        }
        let mut do_undo = false;
        let mut do_redo = false;
        let mut do_repeat = false;
        let mut do_new = false;
        let mut do_name_to_struct = false;
        let mut do_struct_to_name = false;
        let mut do_open = false;
        let mut do_close = false;
        let mut do_save = false;
        let mut do_save_as = false;
        let mut do_page_setup = false;
        let mut do_cut = false;
        let mut do_copy = false;
        let mut do_paste = false;
        let mut do_sel_all = false;
        let mut do_invert_sel = false;
        let mut do_arom = false;
        let mut do_join = false;
        let mut do_cleanup = false;
        let mut do_cleanup_rxn = false;
        let mut do_copy_smiles = false;
        let mut do_copy_cdxml = false;
        let mut do_copy_mol = false;
        let mut do_copy_mol_v2000 = false;
        let mut do_paste_smiles_dlg = false;
        let mut do_paste_cdxml_dlg = false;
        let mut do_text_align: Option<TextAlign> = None;
        let mut do_toggle_bold = false;
        let mut do_toggle_italic = false;
        let mut do_toggle_underline = false;
        let mut do_toggle_sub = false;
        let mut do_toggle_sup = false;
        let mut do_toggle_formula = false;
        let mut do_group = false;
        let mut do_ungroup = false;
        let mut do_front = false;
        let mut do_back = false;
        let mut do_rot180_h = false;
        let mut do_rot180_v = false;
        let mut do_rot_dlg = false;
        let mut do_scale_dlg = false;
        let mut do_toggle_fixed_len = false;
        let mut do_toggle_fixed_ang = false;
        let mut do_toggle_crosshair = false;
        let mut do_toggle_rulers = false;
        let mut do_zoom_mode: Option<u8> = None;
        let mut do_help = false;
        let mut do_open_chem_db = false;
        let mut do_delete = false;
        let mut open_hover_label = false;
        let mut nudge: Option<(f64, f64)> = None;
        let mut typed_chars: Vec<char> = Vec::new();

        ctx.input(|i| {
            let ctrl = i.modifiers.ctrl || i.modifiers.command;
            let shift = i.modifiers.shift;
            let alt = i.modifiers.alt;
            // Exact ChiralCraft shortcuts verified from captured menus!
            if ctrl && shift && !alt && i.key_pressed(egui::Key::Z) { do_redo = true; }
            else if ctrl && !shift && !alt && i.key_pressed(egui::Key::Z) { do_undo = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::Y) { do_repeat = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::N) { do_name_to_struct = true; }
            else if ctrl && alt && !shift && i.key_pressed(egui::Key::N) { do_struct_to_name = true; }
            else if ctrl && !shift && !alt && i.key_pressed(egui::Key::N) { do_new = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::O) { do_open = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::W) { do_close = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::S) { do_save_as = true; }
            else if ctrl && !shift && !alt && i.key_pressed(egui::Key::S) { do_save = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::P) { do_page_setup = true; }
            else if ctrl && !shift && !alt && i.key_pressed(egui::Key::P) {
                do_save_as = true;
            }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::X) { do_cut = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::C) { do_copy = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::V) { do_paste = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::A) { do_sel_all = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::I) { do_invert_sel = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::J) { do_join = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::D) { do_copy_cdxml = true; }
            if ctrl && alt && !shift && i.key_pressed(egui::Key::C) { do_copy_smiles = true; }
            if ctrl && alt && shift && i.key_pressed(egui::Key::O) { do_copy_mol_v2000 = true; }
            else if ctrl && alt && !shift && i.key_pressed(egui::Key::O) { do_copy_mol = true; }
            if ctrl && alt && shift && i.key_pressed(egui::Key::P) { do_paste_cdxml_dlg = true; }
            else if ctrl && alt && !shift && i.key_pressed(egui::Key::P) { do_paste_smiles_dlg = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::L) { do_text_align = Some(TextAlign::Left); }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::C) { do_text_align = Some(TextAlign::Center); }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::R) { do_text_align = Some(TextAlign::Right); }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::J) { do_text_align = Some(TextAlign::Justify); }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::M) { do_text_align = Some(TextAlign::Automatic); }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::B) { do_toggle_bold = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::I) { do_toggle_italic = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::U) { do_toggle_underline = true; }
            if !ctrl && shift && !alt && i.key_pressed(egui::Key::F9) { do_toggle_formula = true; }
            else if !ctrl && !shift && !alt && i.key_pressed(egui::Key::F9) { do_toggle_sub = true; }
            if !ctrl && !shift && !alt && i.key_pressed(egui::Key::F10) { do_toggle_sup = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::L) { do_toggle_fixed_len = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::E) { do_toggle_fixed_ang = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::H) { do_toggle_crosshair = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::H) { do_rot180_h = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::V) { do_rot180_v = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::G) { do_ungroup = true; }
            else if ctrl && !shift && !alt && i.key_pressed(egui::Key::G) { do_group = true; }
            if ctrl && !shift && !alt && i.key_pressed(egui::Key::R) { do_rot_dlg = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::F) { do_open_chem_db = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::K) { do_cleanup = true; }
            else if ctrl && !shift && !alt && i.key_pressed(egui::Key::K) { do_scale_dlg = true; }
            if ctrl && shift && !alt && i.key_pressed(egui::Key::X) { do_cleanup_rxn = true; }
            if alt && !ctrl && i.key_pressed(egui::Key::K) { do_arom = true; }
            if i.key_pressed(egui::Key::F1) { do_help = true; }
            if i.key_pressed(egui::Key::F2) { do_front = true; }
            if i.key_pressed(egui::Key::F3) { do_back = true; }
            if i.key_pressed(egui::Key::F5) { do_zoom_mode = Some(5); }
            if i.key_pressed(egui::Key::F6) { do_zoom_mode = Some(6); }
            if i.key_pressed(egui::Key::F7) { do_zoom_mode = Some(7); }
            if i.key_pressed(egui::Key::F8) { do_zoom_mode = Some(8); }
            if i.key_pressed(egui::Key::F11) { do_toggle_rulers = true; }
            if !ctrl && !alt && i.key_pressed(egui::Key::Enter) { open_hover_label = true; }
            if !ctrl && !alt {
                for ev in &i.events {
                    if let egui::Event::Text(s) = ev {
                        for ch in s.chars() {
                            if !ch.is_control() {
                                typed_chars.push(ch);
                            }
                        }
                    }
                }
            }
            if !ctrl && !alt && i.key_pressed(egui::Key::Space) && !typed_chars.contains(&' ') {
                typed_chars.push(' ');
            }
            if i.key_pressed(egui::Key::Delete) || i.key_pressed(egui::Key::Backspace) {
                do_delete = true;
            }
            nudge = if i.key_pressed(egui::Key::ArrowRight) { Some((0.375, 0.0)) }
            else if i.key_pressed(egui::Key::ArrowLeft) { Some((-0.375, 0.0)) }
            else if i.key_pressed(egui::Key::ArrowDown) { Some((0.0, 0.375)) }
            else if i.key_pressed(egui::Key::ArrowUp) { Some((0.0, -0.375)) }
            else { None };
        });

        if do_undo { self.tab_mut().undo_redo(false); self.sanitize_state(); }
        if do_redo { self.tab_mut().undo_redo(true); self.sanitize_state(); }
        if do_repeat { self.run_repeat_action(); }
        if do_new { self.new_document_tab(DocumentSettings::acs_document_1996()); }
        if do_name_to_struct { self.dialog = ActiveDialog::PasteNameAsStructure; }
        if do_struct_to_name { self.do_structure_to_name_hybrid(ctx, false); }
        if do_open { self.dialog = ActiveDialog::OpenFile; }
        if do_close { self.close_active_tab(); }
        if do_save { self.do_save_action(); }
        if do_save_as {
            if self.file_path_buf.trim().is_empty() {
                self.file_path_buf = format!("{}.cdxml", self.doc().name.replace(' ', "_"));
            }
            self.dialog = ActiveDialog::SaveAs;
        }
        if do_page_setup { self.dialog = ActiveDialog::PageSetup; }
        if do_cut { self.do_cut_action(); }
        if do_copy { self.do_copy_action(); }
        if do_paste {
            self.do_paste_action();
        }
        if do_sel_all {
            let n = self.doc().scene.molecules.len();
            self.tab_mut().selected_mols = (0..n).collect();
        }
        if do_invert_sel {
            let n = self.doc().scene.molecules.len();
            let cur = self.tab().selected_mols.clone();
            self.tab_mut().selected_mols = (0..n).filter(|i| !cur.contains(i)).collect();
        }
        if do_toggle_fixed_len {
            let v = !self.doc().settings.fixed_lengths;
            self.doc_mut().settings.fixed_lengths = v;
        }
        if do_toggle_fixed_ang {
            let v = !self.doc().settings.fixed_angles;
            self.doc_mut().settings.fixed_angles = v;
        }
        if do_toggle_crosshair { self.show_crosshair = !self.show_crosshair; }
        if do_toggle_rulers { self.show_rulers = !self.show_rulers; }
        if let Some(zm) = do_zoom_mode {
            match zm {
                5 => self.zoom = 1.0,
                6 => {
                    let ws = if self.tab().maximized {
                        self.last_desk_size
                    } else {
                        self.tab().win_size
                    };
                    let (pw, ph) = self.tab().effective_page_size_pt();
                    self.zoom = ((ws.x - 30.0) / (pw as f32)).min((ws.y - 45.0) / (ph as f32)).clamp(0.4, 4.0);
                    self.tab_mut().scroll = Vec2::ZERO;
                }
                7 => self.zoom = (self.zoom * 1.25).min(6.0),
                _ => self.zoom = (self.zoom * 0.8).max(0.4),
            }
        }
        if do_group {
            let sel = self.tab().selected_mols.clone();
            if sel.len() >= 2 { self.tab_mut().groups.push(sel); }
        }
        if do_ungroup {
            let sel = self.tab().selected_mols.clone();
            self.tab_mut().groups.retain(|g| !g.iter().any(|mi| sel.contains(mi)));
        }
        if do_front {
            if let Some(&mi) = self.tab().selected_mols.first() {
                let mut moved_to = None;
                self.doc_mut().apply_scene_tx("Bring to Front", |s, _| {
                    if mi < s.molecules.len() {
                        let m = s.molecules.remove(mi);
                        s.molecules.push(m);
                        moved_to = Some(s.molecules.len() - 1);
                    }
                });
                if let Some(ni) = moved_to {
                    self.tab_mut().remap_after_mol_move(mi, ni);
                    self.tab_mut().selected_mols = vec![ni];
                }
            }
        }
        if do_back {
            if let Some(&mi) = self.tab().selected_mols.first() {
                let mut moved = false;
                self.doc_mut().apply_scene_tx("Send to Back", |s, _| {
                    if mi < s.molecules.len() {
                        let m = s.molecules.remove(mi);
                        s.molecules.insert(0, m);
                        moved = true;
                    }
                });
                if moved {
                    self.tab_mut().remap_after_mol_move(mi, 0);
                    self.tab_mut().selected_mols = vec![0];
                }
            }
        }
        if do_rot180_h {
            self.do_flip_or_rot180("Rotate 180° Horizontal", FlipAxis::Horizontal, true);
        }
        if do_rot180_v {
            self.do_flip_or_rot180("Rotate 180° Vertical", FlipAxis::Vertical, true);
        }
        if do_rot_dlg { self.dialog = ActiveDialog::RotateDialog; }
        if do_scale_dlg { self.dialog = ActiveDialog::ScaleDialog; }
        if do_help { self.dialog = ActiveDialog::Shortcuts; }
        if do_open_chem_db { self.show_chem_db = true; }
        if do_copy_smiles {
            if let Some(m) = self.active_mol() {
                let smi = molecule_to_smiles(m);
                ctx.copy_text(smi.clone());
                self.status_msg = format!("Copied SMILES: {smi}");
            }
        }
        if do_copy_cdxml {
            ctx.copy_text(self.doc().to_cdxml().to_xml_string());
            self.status_msg = "Copied CDXML Text (Ctrl+D)".to_string();
        }
        if do_copy_mol {
            if let Some(m) = self.active_mol() {
                ctx.copy_text(molecule_to_molfile_v3000(m, &self.doc().name));
                self.status_msg = "Copied MOL Text (Alt+Ctrl+O)".to_string();
            }
        }
        if do_copy_mol_v2000 {
            if let Some(m) = self.active_mol() {
                ctx.copy_text(molecule_to_molfile_v2000(m, &self.doc().name));
                self.status_msg = "Copied MOL V2000 Text (Alt+Shift+Ctrl+O)".to_string();
            }
        }
        if do_paste_smiles_dlg {
            self.dialog = ActiveDialog::PasteSmiles;
        }
        if do_paste_cdxml_dlg {
            self.multiline_input.clear();
            self.dialog = ActiveDialog::PasteCdxml;
        }
        if let Some(ta) = do_text_align {
            self.text_align = ta;
        }
        if do_toggle_bold { self.text_bold = !self.text_bold; }
        if do_toggle_italic { self.text_italic = !self.text_italic; }
        if do_toggle_underline { self.text_underline = !self.text_underline; }
        if do_toggle_sub { self.text_sub = !self.text_sub; if self.text_sub { self.text_sup = false; } }
        if do_toggle_sup { self.text_sup = !self.text_sup; if self.text_sup { self.text_sub = false; } }
        if do_toggle_formula { self.text_formula = !self.text_formula; }
        if do_delete && !self.tab().selected_mols.is_empty() {
            let mut sel = self.tab().selected_mols.clone();
            sel.sort_unstable();
            sel.dedup(); // a duplicate index would delete an extra, unselected molecule
            sel.reverse();
            let nmols = self.doc().scene.molecules.len();
            sel.retain(|&mi| mi < nmols);
            let removed = sel.clone();
            self.doc_mut().apply_scene_tx("Delete", |scene, _| {
                for mi in sel {
                    scene.molecules.remove(mi);
                }
            });
            self.tab_mut().remap_after_mol_removal(&removed);
            self.tab_mut().selected_mols.clear();
            self.sanitize_state();
        }
        if let Some((nx, ny)) = nudge {
            let targets: Vec<usize> = if self.tab().selected_mols.is_empty() {
                (0..self.doc().scene.molecules.len()).collect()
            } else {
                self.tab().selected_mols.clone()
            };
            self.doc_mut().apply_scene_tx("Nudge", |scene, _| {
                for mi in targets {
                    if let Some(m) = scene.molecules.get_mut(mi) {
                        for a in &mut m.atoms { a.x += nx; a.y += ny; }
                    }
                }
            });
        }
        if do_cleanup {
            self.do_clean_up_structure();
        }
        if do_cleanup_rxn {
            self.doc_mut().apply_scene_tx("Clean Up Reaction", |scene, _| { clean_up_reaction(scene); });
            self.last_repeat_action = Some("Clean Up Reaction".to_string());
        }
        if open_hover_label {
            if let Some((mi, ai)) = self.hovered_atom {
                self.open_inline_atom_edit(mi, ai);
            }
        }
        // Dispatch all 98 official `hotkeys.xml` context-sensitive bindings (falling back to generic hotkeys)!
        for ch in typed_chars {
            let outcome = if let Some((mi, ai)) = self.hovered_atom {
                match apply_atom_hotkey(self.doc_mut(), mi, ai, ch) {
                    HotkeyOutcome::NotMatched => apply_generic_hotkey(ch),
                    other => other,
                }
            } else if let Some((mi, bi)) = self.hovered_bond {
                match apply_bond_hotkey(self.doc_mut(), mi, bi, ch) {
                    HotkeyOutcome::NotMatched => apply_generic_hotkey(ch),
                    other => other,
                }
            } else {
                apply_generic_hotkey(ch)
            };
            match outcome {
                HotkeyOutcome::Handled(msg) => self.status_msg = msg,
                HotkeyOutcome::SwitchTool(t) => {
                    self.tool = t;
                    self.active_element_z = None;
                }
                HotkeyOutcome::OpenChooseNickname => {
                    self.dialog = ActiveDialog::ListNicknames;
                }
                HotkeyOutcome::OpenPropertiesDialog => {
                    if let Some((mi, ai)) = self.hovered_atom {
                        self.dialog = ActiveDialog::AtomProperties { mol_idx: mi, atom_idx: ai };
                    } else if !self.doc().scene.molecules.is_empty() && !self.doc().scene.molecules[0].atoms.is_empty() {
                        self.dialog = ActiveDialog::AtomProperties { mol_idx: 0, atom_idx: 0 };
                    }
                }
                HotkeyOutcome::SelectMolecule(mi) => {
                    self.tab_mut().selected_mols = vec![mi];
                }
                HotkeyOutcome::SetDoubleAlign(key, align) => {
                    self.tab_mut().double_bond_align.insert(key, align);
                    self.status_msg = format!("Double bond alignment → {align:?}");
                }
                HotkeyOutcome::ToggleAtomNumbers => {
                    let cur = self.tab().show_atom_numbers;
                    self.tab_mut().show_atom_numbers = !cur;
                    self.status_msg = format!("Atom numbers: {}", if !cur { "ON" } else { "OFF" });
                }
                HotkeyOutcome::NotMatched => {}
            }
        }
        if do_arom { self.toggle_aromatic_display_tx(); }
        if do_join { self.join_close_atoms_tx(); }
    }
}

fn mol_bbox(mol: &ChemMolecule) -> (f64, f64, f64, f64) {
    let mut x1 = f64::INFINITY;
    let mut y1 = f64::INFINITY;
    let mut x2 = f64::NEG_INFINITY;
    let mut y2 = f64::NEG_INFINITY;
    for a in &mol.atoms {
        x1 = x1.min(a.x);
        y1 = y1.min(a.y);
        x2 = x2.max(a.x);
        y2 = y2.max(a.y);
    }
    if x1.is_finite() { (x1, y1, x2, y2) } else { (0.0, 0.0, 0.0, 0.0) }
}

#[cfg(test)]
mod audit_a_tests {
    use super::*;

    fn tab_with_side_tables() -> OpenDocTab {
        let d = InteractiveDocument::new_saved_seed("t", DocumentSettings::acs_document_1996());
        let mut t = OpenDocTab::new_at(d, 0);
        for mi in 0..4 {
            t.ring_fills.insert(mi, Color32::from_rgb(mi as u8, 0, 0));
            t.reaction_atom_maps.insert((mi, 1), mi as u32 + 10);
            t.linked_captions.push(LinkedCaption {
                caption_id: mi as u32,
                molecule_idx: mi,
                kind: LinkedCaptionKind::Formula,
                decimals: 2,
            });
        }
        t.selected_mols = vec![1, 3];
        t.groups = vec![vec![0, 1], vec![2, 3]];
        t
    }

    #[test]
    fn mol_removal_rekeys_side_tables() {
        let mut t = tab_with_side_tables();
        t.remap_after_mol_removal(&[1]);
        assert_eq!(t.ring_fills.len(), 3);
        assert_eq!(t.ring_fills[&1], Color32::from_rgb(2, 0, 0));
        assert_eq!(t.ring_fills[&2], Color32::from_rgb(3, 0, 0));
        assert_eq!(t.reaction_atom_maps[&(2, 1)], 13);
        let lc: Vec<(u32, usize)> = t.linked_captions.iter().map(|l| (l.caption_id, l.molecule_idx)).collect();
        assert_eq!(lc, vec![(0, 0), (2, 1), (3, 2)]);
        assert_eq!(t.selected_mols, vec![2]);
        assert_eq!(t.groups, vec![vec![1, 2]]);
    }

    #[test]
    fn undo_redo_revert_restore_side_tables() {
        let mut t = tab_with_side_tables();
        t.doc.apply_scene_tx("Seed", |s, _| {
            for _ in 0..4 {
                s.molecules.push(ChemMolecule::default());
            }
        });
        t.doc.save();
        t.sync_side_history(); // frame boundary
        let before = t.ring_fills.clone();
        // Delete molecule 1 (as Edit > Clear does), then a frame passes.
        t.doc.apply_scene_tx("Clear Selection", |s, _| {
            s.molecules.remove(1);
        });
        t.remap_after_mol_removal(&[1]);
        t.sync_side_history();
        let after = t.ring_fills.clone();
        assert!(t.undo_redo(false));
        assert_eq!(t.ring_fills, before, "undo must re-key fills back onto the restored molecules");
        assert!(t.undo_redo(true));
        assert_eq!(t.ring_fills, after);
        t.ring_fills.insert(0, Color32::WHITE);
        assert!(t.revert());
        assert_eq!(t.ring_fills, before);
        assert_eq!(t.doc.scene.molecules.len(), 4);
    }

    #[test]
    fn paste_of_several_structures_is_one_undo_and_selects_them() {
        let mut app = ChiralCraftApp::new();
        let m = smiles_to_molecule("CC").unwrap();
        app.clipboard_mols = vec![m.clone(), m];
        app.do_paste_action();
        assert_eq!(app.doc().scene.molecules.len(), 2);
        assert_eq!(app.tab().selected_mols, vec![0, 1]);
        assert!(app.doc_mut().undo());
        assert!(app.doc().scene.molecules.is_empty());
    }

    #[test]
    fn closing_dirty_tab_asks_first() {
        let mut app = ChiralCraftApp::new();
        app.new_document_tab(DocumentSettings::acs_document_1996());
        app.doc_mut().apply_scene_tx("Bond", |s, _| s.molecules.push(ChemMolecule::default()));
        app.close_active_tab();
        assert_eq!(app.tabs.len(), 2, "dirty tab must not close silently");
        assert_eq!(app.dialog, ActiveDialog::ConfirmCloseTab);
        app.force_close_active_tab();
        assert_eq!(app.tabs.len(), 1);
    }

    #[test]
    fn two_steps_in_one_frame_keep_side_tables_for_each_undo() {
        let mut app = ChiralCraftApp::new();
        app.doc_mut().apply_scene_tx("Seed", |s, _| {
            for _ in 0..3 {
                s.molecules.push(ChemMolecule::default());
            }
        });
        app.tab_mut().sync_side_history(); // frame boundary
        app.tab_mut().ring_fills.insert(2, Color32::RED);
        app.tab_mut().selected_mols = vec![0];
        app.tab_mut().sync_side_history(); // frame boundary
        // F2 then F3 handled in the same frame (no frame-boundary sync between them).
        app.doc_mut().apply_scene_tx("Bring to Front", |s, _| {
            let m = s.molecules.remove(0);
            s.molecules.push(m);
        });
        app.tab_mut().remap_after_mol_move(0, 2);
        app.tab_mut().selected_mols = vec![2];
        app.doc_mut().apply_scene_tx("Send to Back", |s, _| {
            let m = s.molecules.remove(2);
            s.molecules.insert(0, m);
        });
        app.tab_mut().remap_after_mol_move(2, 0);
        app.tab_mut().selected_mols = vec![0];
        assert!(app.tab_mut().undo_redo(false));
        assert_eq!(app.tab().ring_fills.get(&1), Some(&Color32::RED), "fill must follow the molecule after Bring to Front");
        assert_eq!(app.tab().selected_mols, vec![2]);
        assert!(app.tab_mut().undo_redo(false));
        assert_eq!(app.tab().ring_fills.get(&2), Some(&Color32::RED));
        assert_eq!(app.tab().selected_mols, vec![0]);
    }

    #[test]
    fn mol_move_and_merge_rekey_side_tables() {
        let mut t = tab_with_side_tables();
        t.remap_after_mol_move(0, 3); // bring to front
        assert_eq!(t.ring_fills[&3], Color32::from_rgb(0, 0, 0));
        assert_eq!(t.ring_fills[&0], Color32::from_rgb(1, 0, 0));
        t.remap_after_mol_move(3, 0); // send back: identity again
        assert_eq!(t.ring_fills[&0], Color32::from_rgb(0, 0, 0));
        assert_eq!(t.ring_fills[&3], Color32::from_rgb(3, 0, 0));

        let mut t = tab_with_side_tables();
        t.remap_after_mol_merge(0, 2, 5, 4);
        assert_eq!(t.reaction_atom_maps[&(0, 6)], 12); // (2,1) -> (0, 1+5)
        assert_eq!(t.reaction_atom_maps[&(2, 1)], 13); // mol 3 shifted down
        assert!(t.linked_captions.iter().any(|l| l.caption_id == 2 && l.molecule_idx == 0));
        assert_eq!(t.ring_fills[&0], Color32::from_rgb(0, 0, 0), "keep's own fill wins");
        let mut u = tab_with_side_tables();
        u.ring_fills.remove(&0);
        u.remap_after_mol_merge(0, 2, 5, 4);
        assert_eq!(u.ring_fills[&0], Color32::from_rgb(2, 0, 0), "unfilled keep adopts the absorbed fill");

        t.remap_after_atom_removal(0, 1);
        assert!(!t.reaction_atom_maps.contains_key(&(0, 1)));
        assert_eq!(t.reaction_atom_maps[&(0, 5)], 12);
    }

    #[test]
    fn roundtrip_save_writes_bound_tab_not_active_tab() {
        let path = std::env::temp_dir().join(format!("cc_audit_rt_{}.svg", std::process::id()));
        let _ = fs::remove_file(&path);
        let ps = path.to_str().unwrap().to_string();
        let mut app = ChiralCraftApp::new_with_file(&ps, true);
        assert_eq!(app.roundtrip_svg_path.as_deref(), Some(ps.as_str()));
        app.new_document_tab(DocumentSettings::acs_document_1996());
        let m = smiles_to_molecule("CCO").unwrap();
        app.doc_mut().stamp_molecule_tx(m, "x");
        // Ctrl+S on the unrelated tab must not touch the host SVG.
        app.tabs[app.active_tab].file_path = None;
        app.do_save_action();
        assert!(!path.exists());
        app.save_roundtrip_svg().unwrap();
        let written = fs::read_to_string(&path).unwrap();
        let expected = app.tabs[0].doc.to_svg_styled(true, false, app.rear_opacity, app.crossing_clearance_pt);
        assert_eq!(written, expected);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn plain_save_to_lossy_format_falls_back_to_save_as() {
        let path = std::env::temp_dir().join(format!("cc_audit_lossy_{}.eps", std::process::id()));
        let _ = fs::remove_file(&path);
        let mut app = ChiralCraftApp::new();
        let m = smiles_to_molecule("CCO").unwrap();
        app.doc_mut().stamp_molecule_tx(m, "x");
        app.tab_mut().file_path = Some(path.to_string_lossy().into_owned());
        app.do_save_action();
        assert!(!path.exists(), "lossy re-save must not write and mark clean");
        assert!(matches!(app.dialog, ActiveDialog::SaveAs));
        assert!(app.file_path_buf.ends_with(".cdxml"));
    }
}

impl eframe::App for ChiralCraftApp {
    fn on_exit(&mut self, _gl: Option<&eframe::glow::Context>) {
        if self.roundtrip_svg_path.is_some() && !self.roundtrip_saved_on_close {
            self.roundtrip_saved_on_close = true;
            let _ = self.save_roundtrip_svg();
        }
    }

    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.render_ui_frame(ctx);
    }
}

impl ChiralCraftApp {
    fn render_ui_frame(&mut self, ctx: &egui::Context) {
        for tab in &mut self.tabs {
            tab.sync_side_history();
        }
        if !self.style_initialized {
            Self::apply_theme(ctx, self.ui_theme);
            self.style_initialized = true;
        }
        if let Some(out_path) = self.screenshot_out_path.clone() {
            let captured = ctx.input(|i| {
                i.events.iter().find_map(|ev| {
                    if let egui::Event::Screenshot { image, .. } = ev {
                        Some(image.clone())
                    } else {
                        None
                    }
                })
            });
            if let Some(img) = captured {
                if let Some(parent) = std::path::Path::new(&out_path).parent() {
                    if !parent.as_os_str().is_empty() {
                        let _ = fs::create_dir_all(parent);
                    }
                }
                let write_res = chiralcraft_formats::encode_rgba_png(
                    img.as_raw(),
                    img.width() as u32,
                    img.height() as u32,
                )
                .and_then(|png_bytes| atomic_write_bytes(&out_path, &png_bytes));
                if let Err(e) = write_res {
                    eprintln!("Screenshot export failed for {out_path}: {e}");
                    std::process::exit(1);
                }
                self.screenshot_out_path = None;
                self.exit_confirmed = true;
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            } else {
                self.screenshot_frames = self.screenshot_frames.saturating_add(1);
                if self.screenshot_frames == 6 {
                    ctx.send_viewport_cmd(egui::ViewportCommand::Screenshot(egui::UserData::default()));
                }
                ctx.request_repaint();
            }
        }
        if ctx.input(|i| i.viewport().close_requested())
            && self.roundtrip_svg_path.is_some()
            && !self.roundtrip_saved_on_close
        {
            match self.save_roundtrip_svg() {
                Ok(()) => self.roundtrip_saved_on_close = true,
                Err(e) => {
                    // Abort the close so the edits are not silently discarded.
                    ctx.send_viewport_cmd(egui::ViewportCommand::CancelClose);
                    self.status_msg = format!("Round-trip save failed: {e} — window kept open");
                }
            }
        }
        if ctx.input(|i| i.viewport().close_requested())
            && !self.exit_confirmed
            && self.tabs.iter().any(|t| t.doc.is_dirty())
        {
            ctx.send_viewport_cmd(egui::ViewportCommand::CancelClose);
            self.dialog = ActiveDialog::ConfirmExit;
            // Staying open: later round-trip edits must be written again on the next close.
            self.roundtrip_saved_on_close = false;
        }
        self.sanitize_state();
        self.sync_linked_captions();
        self.apply_shortcuts(ctx);

        // OS File Drag-and-Drop directly onto the ChiralCraft window!
        let dropped: Vec<PathBuf> = ctx.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .filter_map(|f| f.path.clone())
                .collect()
        });
        for p in dropped {
            if let Some(ps) = p.to_str() {
                self.open_path_into_new_tab(ps);
            }
        }

        // Round-trip embedded SVG editing banner (when launched by WordCraft / DeckCraft via `--edit-svg`)
        if let Some(rt_path) = self.roundtrip_svg_path.clone() {
            egui::TopBottomPanel::top("roundtrip_svg_banner")
                .frame(
                    egui::Frame::side_top_panel(&ctx.style())
                        .fill(Color32::from_rgb(224, 238, 255))
                        .inner_margin(egui::Margin::symmetric(8, 3)),
                )
                .show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        let fname = std::path::Path::new(&rt_path)
                            .file_name()
                            .and_then(|s| s.to_str())
                            .unwrap_or(&rt_path);
                        ui.strong(format!("Round-Trip Editable SVG: {fname}"));
                        ui.label("— Save (Ctrl+S) or close window to apply edits back to host document.");
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if ui.button("Save & Return").clicked() {
                                match self.save_roundtrip_svg() {
                                    Ok(()) => {
                                        self.roundtrip_saved_on_close = true;
                                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                                    }
                                    // Keep the window open so unsaved edits are not lost.
                                    Err(e) => self.status_msg = format!("Round-trip save failed: {e}"),
                                }
                            }
                            if ui.button("Save Now").clicked() {
                                if let Err(e) = self.save_roundtrip_svg() {
                                    self.status_msg = format!("Round-trip save failed: {e}");
                                }
                            }
                        });
                    });
                });
        }

        // ---- 1. Exact 11-Menu ChiralCraft Menu Bar (1:1 parity with all 11 captured menus) ----
        egui::TopBottomPanel::top("menu_bar")
            .frame(egui::Frame::side_top_panel(&ctx.style()).fill(Color32::WHITE).inner_margin(egui::Margin::symmetric(6, 2)))
            .show(ctx, |ui| {
                egui::menu::bar(ui, |ui| {
                    ui.menu_button("File", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        if ui.button("New ACS Document 1996\tCtrl+N").clicked() {
                            self.new_document_tab(DocumentSettings::acs_document_1996());
                            ui.close_menu();
                        }
                        if ui.button("Open...\tCtrl+O").clicked() { self.dialog = ActiveDialog::OpenFile; ui.close_menu(); }
                        ui.menu_button("Open Style Sheets", |ui| {
                            let styles = self.stylesheets.clone();
                            for (lbl, settings) in styles {
                                if ui.button(&lbl).clicked() {
                                    self.new_document_tab(settings);
                                    ui.close_menu();
                                }
                            }
                        });
                        ui.menu_button("Open Templates", |ui| {
                            let names: Vec<String> = self.ctp_libraries.iter().map(|l| l.name.clone()).collect();
                            egui::ScrollArea::vertical().max_height(420.0).show(ui, |ui| {
                                for (i, n) in names.iter().enumerate() {
                                    if ui.button(n).clicked() {
                                        self.selected_ctp_idx = i;
                                        self.active_flyout = ActiveFlyout::Templates;
                                        ui.close_menu();
                                    }
                                }
                                ui.separator();
                                if ui.button("New Templates").clicked() {
                                    let custom_idx = if let Some(idx) = self.ctp_libraries.iter().position(|l| l.name == "New Templates") {
                                        idx
                                    } else {
                                        self.ctp_libraries.push(TemplateLibrary {
                                            name: "New Templates".to_string(),
                                            templates: self.doc().scene.molecules.first().cloned().into_iter().collect(),
                                            graphic_templates: Vec::new(),
                                        });
                                        self.ctp_libraries.len() - 1
                                    };
                                    self.selected_ctp_idx = custom_idx;
                                    self.active_flyout = ActiveFlyout::Templates;
                                    self.status_msg = "Opened New Templates library".to_string();
                                    ui.close_menu();
                                }
                            });
                        });
                        ui.menu_button("Open Samples", |ui| {
                            for (lbl, smi) in [
                                ("Aspirin (C9H8O4)", "CC(=O)Oc1ccccc1C(=O)O"),
                                ("Caffeine (C8H10N4O2)", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
                                ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
                                ("Acetaminophen", "CC(=O)Nc1ccc(O)cc1"),
                                ("Naphthalene", "c1ccc2ccccc2c1"),
                                ("Testosterone Core", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"),
                                ("Penicillin G Core", "CC1(C)SC2C(N)C(=O)N2C1C(=O)O"),
                            ] {
                                if ui.button(lbl).clicked() {
                                    if let Ok(mut m) = smiles_to_molecule(smi) {
                                        let is_blank_untitled = self.doc().scene.molecules.is_empty()
                                            && self.doc().scene.graphics.is_empty()
                                            && self.doc().scene.captions.is_empty()
                                            && self.tab().file_path.is_none();
                                        if !is_blank_untitled {
                                            let d = InteractiveDocument::new_saved_seed(lbl, DocumentSettings::acs_document_1996());
                                            let idx = self.tabs.len();
                                            self.tabs.push(OpenDocTab::new_at(d, idx));
                                            self.active_tab = self.tabs.len() - 1;
                                            self.reset_transient_interaction_state();
                                        } else {
                                            self.doc_mut().name = lbl.to_string();
                                        }
                                        let blen = self.doc().settings.bond_length_pt;
                                        clean_up_structure_preserving_stereo(&mut m, &[], blen);
                                        for a in &mut m.atoms { a.x += 180.0; a.y += 150.0; }
                                        self.doc_mut().stamp_molecule_tx(m, lbl);
                                        self.doc_mut().save();
                                        self.status_msg = format!("Opened sample: {lbl}");
                                    }
                                    ui.close_menu();
                                }
                            }
                        });
                        ui.separator();
                        if ui.button("Close\tCtrl+W").clicked() { self.close_active_tab(); ui.close_menu(); }
                        let can_save = self.doc().is_dirty() || self.tab().file_path.is_none();
                        if ui.add_enabled(can_save, egui::Button::new("Save\tCtrl+S")).clicked() {
                            self.do_save_action();
                            ui.close_menu();
                        }
                        if ui.button("Save As\tShift+Ctrl+S").clicked() { self.dialog = ActiveDialog::SaveAs; ui.close_menu(); }
                        if ui.button("Compare Documents (Chemical Diff)...").clicked() {
                            self.show_diff_window = true;
                            if self.tabs.len() >= 2 {
                                self.diff_left_tab = 0;
                                self.diff_right_tab = 1;
                                self.last_diff = Some(compare_documents(&self.tabs[0].doc, &self.tabs[1].doc));
                            }
                            ui.close_menu();
                        }
                        if ui.button("Virtual Screening / SDF Library Browser...").clicked() {
                            self.omni_suite.show_screening = true;
                            ui.close_menu();
                        }
                        if ui.button("Export LaTeX (\\chemfig / TikZ)...").clicked() {
                            self.omni_suite.show_latex_export = true;
                            ui.close_menu();
                        }
                        if !self.user_config.recent_files.is_empty() {
                            let recents = self.user_config.recent_files.clone();
                            ui.menu_button("Open Recent", |ui| {
                                for rpath in recents {
                                    let lbl = std::path::Path::new(&rpath)
                                        .file_name()
                                        .and_then(|s| s.to_str())
                                        .unwrap_or(&rpath)
                                        .to_string();
                                    if ui.button(&lbl).on_hover_text(&rpath).clicked() {
                                        self.open_path_into_new_tab(&rpath);
                                        ui.close_menu();
                                    }
                                }
                            });
                        }
                        let can_revert = self.doc().is_dirty();
                        if ui.add_enabled(can_revert, egui::Button::new("Revert")).clicked() {
                            self.tab_mut().revert();
                            self.sanitize_state();
                            self.status_msg = "Reverted to last saved state".to_string();
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.button("Page Setup...\tShift+Ctrl+P").clicked() { self.dialog = ActiveDialog::PageSetup; ui.close_menu(); }
                        if ui.button("Print...\tCtrl+P").clicked() {
                            self.file_path_buf = "print_export.pdf".to_string();
                            self.dialog = ActiveDialog::SaveAs;
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.button("Document Settings...").clicked() { self.dialog = ActiveDialog::DocumentSettings; ui.close_menu(); }
                        if ui.button("Document Annotations...").clicked() { self.dialog = ActiveDialog::Annotate; ui.close_menu(); }
                        if ui.button("Document Properties...").clicked() { self.dialog = ActiveDialog::DocumentProperties; ui.close_menu(); }
                        ui.menu_button("Apply Document Settings from", |ui| {
                            if ui.button("Other...").clicked() {
                                self.dialog = ActiveDialog::ApplyStylesheetFile;
                                ui.close_menu();
                            }
                            ui.separator();
                            let styles = self.stylesheets.clone();
                            for (lbl, settings) in styles {
                                if ui.button(format!("{} ({:.2} pt)", lbl, settings.bond_length_pt)).clicked() {
                                    replace_settings_tx(self.doc_mut(), settings, "Apply Document Settings");
                                    self.status_msg = format!("Applied document settings from {lbl}");
                                    ui.close_menu();
                                }
                            }
                        });
                        if ui.button("Preferences...").clicked() { self.dialog = ActiveDialog::Preferences; ui.close_menu(); }
                        if ui.button("List Nicknames...").clicked() { self.dialog = ActiveDialog::ListNicknames; ui.close_menu(); }
                        ui.separator();
                        if ui.button("Exit ChiralCraft\tAlt+F4").clicked() {
                            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                            ui.close_menu();
                        }
                        }); }
                    });

                    ui.menu_button("Edit", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        let u_raw = self.doc().undo_menu_label();
                        let r_raw = self.doc().redo_menu_label();
                        let can_undo = !u_raw.contains("not available");
                        let can_redo = !r_raw.contains("not available");
                        let has_mols = !self.doc().scene.molecules.is_empty();
                        let has_any = has_mols
                            || !self.doc().scene.arrows.is_empty()
                            || !self.doc().scene.curves.is_empty()
                            || !self.doc().scene.graphics.is_empty()
                            || !self.doc().scene.captions.is_empty()
                            || !self.doc().scene.tables.is_empty()
                            || !self.doc().scene.tlc_plates.is_empty()
                            || !self.doc().scene.energy_profiles.is_empty()
                            || !self.tab().gel_plates.is_empty();
                        let has_clip = !self.clipboard_mols.is_empty();
                        let u_lbl = format!("{u_raw}\tCtrl+Z");
                        let r_lbl = format!("{r_raw}\tShift+Ctrl+Z");
                        if ui.add_enabled(can_undo, egui::Button::new(u_lbl)).clicked() { self.tab_mut().undo_redo(false); self.sanitize_state(); ui.close_menu(); }
                        if ui.add_enabled(can_redo, egui::Button::new(r_lbl)).clicked() { self.tab_mut().undo_redo(true); self.sanitize_state(); ui.close_menu(); }
                        ui.separator();
                        if ui.add_enabled(has_mols, egui::Button::new("Cut\tCtrl+X")).clicked() {
                            self.do_cut_action();
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Copy\tCtrl+C")).clicked() {
                            self.do_copy_action();
                            ui.close_menu();
                        }
                        let _ = has_clip;
                        if ui.button("Paste\tCtrl+V").clicked() {
                            self.do_paste_action();
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_any, egui::Button::new("Clear\tDel")).clicked() {
                            self.do_clear_action();
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Select All\tCtrl+A")).clicked() {
                            let n = self.doc().scene.molecules.len();
                            self.tab_mut().selected_mols = (0..n).collect();
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Invert Selection\tShift+Ctrl+I")).clicked() {
                            let n = self.doc().scene.molecules.len();
                            let cur = self.tab().selected_mols.clone();
                            self.tab_mut().selected_mols = (0..n).filter(|i| !cur.contains(i)).collect();
                            ui.close_menu();
                        }
                        let rep_lbl = match &self.last_repeat_action {
                            Some(a) => format!("Repeat {a}\tCtrl+Y"),
                            None => "Repeat not available\tCtrl+Y".to_string(),
                        };
                        if ui.add_enabled(self.last_repeat_action.is_some(), egui::Button::new(rep_lbl)).clicked() {
                            self.run_repeat_action();
                            ui.close_menu();
                        }
                        ui.separator();
                        ui.menu_button("Copy As", |ui| {
                            if ui.add_enabled(has_any, egui::Button::new("Editable SVG (WordCraft / DeckCraft)")).clicked() {
                                self.do_copy_editable_svg(ctx);
                                ui.close_menu();
                            }
                            ui.separator();
                            if ui.add_enabled(has_mols, egui::Button::new("SMILES\tAlt+Ctrl+C")).clicked() {
                                if let Some(m) = self.active_mol() {
                                    let smi = molecule_to_smiles(m);
                                    ctx.copy_text(smi.clone());
                                    self.status_msg = format!("Copied SMILES: {smi}");
                                }
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("SLN")).clicked() {
                                if let Some(m) = self.active_mol() {
                                    ctx.copy_text(molecule_to_sln(m));
                                    self.status_msg = "Copied SLN".to_string();
                                }
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("InChI")).clicked() {
                                if let Some(m) = self.active_mol() {
                                    let inc = molecule_to_inchi(m);
                                    ctx.copy_text(inc.clone());
                                    self.status_msg = format!("Copied InChI: {inc}");
                                }
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("InChI Key")).clicked() {
                                if let Some(m) = self.active_mol() {
                                    let key = inchi_to_inchikey(&molecule_to_inchi(m));
                                    ctx.copy_text(key.clone());
                                    self.status_msg = format!("Copied InChIKey: {key}");
                                }
                                ui.close_menu();
                            }
                            if ui.button("CDXML Text\tCtrl+D").clicked() {
                                ctx.copy_text(self.doc().to_cdxml().to_xml_string());
                                self.status_msg = "Copied CDXML".to_string();
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("MOL Text\tAlt+Ctrl+O")).clicked() {
                                if let Some(m) = self.active_mol() {
                                    ctx.copy_text(molecule_to_molfile_v3000(m, &self.doc().name));
                                    self.status_msg = "Copied MOL V3000".to_string();
                                }
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("MOL V2000 Text\tAlt+Shift+Ctrl+O")).clicked() {
                                if let Some(m) = self.active_mol() {
                                    ctx.copy_text(molecule_to_molfile_v2000(m, &self.doc().name));
                                    self.status_msg = "Copied MOL V2000".to_string();
                                }
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("HELM")).clicked() {
                                if let Some(m) = self.active_mol() {
                                    let h = molecule_to_helm(m);
                                    ctx.copy_text(h.clone());
                                    self.status_msg = format!("Copied HELM: {h}");
                                }
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("FASTA")).clicked() {
                                if let Some(m) = self.active_mol() {
                                    ctx.copy_text(molecule_to_fasta(m, &self.doc().name));
                                    self.status_msg = "Copied FASTA".to_string();
                                }
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("3MF (Ball and Stick)")).clicked() {
                                if let Some(m) = self.active_mol() {
                                    let xml = molecule_to_3mf_xml(m, true);
                                    ctx.copy_text(xml);
                                    self.status_msg = "Copied 3MF (Ball and Stick) 3D model to clipboard".to_string();
                                }
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("3MF (Stick)")).clicked() {
                                if let Some(m) = self.active_mol() {
                                    let xml = molecule_to_3mf_xml(m, false);
                                    ctx.copy_text(xml);
                                    self.status_msg = "Copied 3MF (Stick) 3D model to clipboard".to_string();
                                }
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("PNG (Editable)...")).clicked() {
                                match self.doc().to_editable_png(self.export_opts.png_dpi) {
                                    Ok(bytes) => {
                                        if let Some(chosen) = native_dialogs::pick_save_file("Save Editable PNG", "chiralcraft.png") {
                                            match atomic_write_bytes(&chosen, &bytes) {
                                                Ok(()) => self.status_msg = format!("Exported editable PNG ({} bytes, {} DPI) to {chosen}", bytes.len(), self.export_opts.png_dpi),
                                                Err(e) => self.status_msg = format!("PNG save failed: {e}"),
                                            }
                                        } else {
                                            self.status_msg = "PNG export cancelled".to_string();
                                        }
                                    }
                                    Err(e) => self.status_msg = format!("PNG export failed: {e}"),
                                }
                                ui.close_menu();
                            }
                        });
                        ui.menu_button("Paste Special", |ui| {
                            if ui.button("SMILES\tAlt+Ctrl+P").clicked() { self.dialog = ActiveDialog::PasteSmiles; ui.close_menu(); }
                            if ui.button("Name as Structure\tShift+Ctrl+N").clicked() { self.dialog = ActiveDialog::PasteNameAsStructure; ui.close_menu(); }
                            if ui.button("InChI").clicked() { self.dialog = ActiveDialog::PasteInchi; ui.close_menu(); }
                            if ui.button("Editable SVG (WordCraft / DeckCraft)...").clicked() {
                                self.multiline_input.clear();
                                self.dialog = ActiveDialog::PasteEditableSvg;
                                ui.close_menu();
                            }
                            if ui.button("CDXML Text\tAlt+Shift+Ctrl+P").clicked() {
                                self.multiline_input.clear();
                                self.dialog = ActiveDialog::PasteCdxml;
                                ui.close_menu();
                            }
                            if ui.button("MDL Text").clicked() {
                                self.multiline_input.clear();
                                self.dialog = ActiveDialog::PasteMolText;
                                ui.close_menu();
                            }
                            ui.separator();
                            if ui.button("Biopolymer").clicked() { self.dialog = ActiveDialog::PasteBiopolymer(PolymerKind::Peptide); ui.close_menu(); }
                            if ui.button("FASTA Peptide").clicked() { self.dialog = ActiveDialog::PasteBiopolymer(PolymerKind::Peptide); ui.close_menu(); }
                            if ui.button("FASTA DNA").clicked() { self.dialog = ActiveDialog::PasteBiopolymer(PolymerKind::Dna); ui.close_menu(); }
                            if ui.button("FASTA RNA").clicked() { self.dialog = ActiveDialog::PasteBiopolymer(PolymerKind::Rna); ui.close_menu(); }
                            if ui.button("HELM").clicked() { self.dialog = ActiveDialog::PasteHelm; ui.close_menu(); }
                        });
                        ui.separator();
                        if ui.button("Insert File...").clicked() { self.dialog = ActiveDialog::InsertFile; ui.close_menu(); }
                        if ui.button("Insert Object...").clicked() { self.dialog = ActiveDialog::InsertObject; ui.close_menu(); }
                        let has_embedded_obj = !self.doc().scene.graphics.is_empty()
                            || !self.doc().scene.tables.is_empty()
                            || !self.doc().scene.tlc_plates.is_empty()
                            || !self.doc().scene.energy_profiles.is_empty()
                            || !self.doc().scene.curves.is_empty()
                            || !self.doc().scene.captions.is_empty();
                        if has_embedded_obj {
                            ui.menu_button("Object", |ui| {
                                if ui.button("Edit Embedded Object...").clicked() {
                                    if !self.doc().scene.captions.is_empty() {
                                        let ci = self.doc().scene.captions.len() - 1;
                                        self.label_input = self.doc().scene.captions[ci].text.clone();
                                        self.dialog = ActiveDialog::EditCaption { cap_idx: ci };
                                    } else {
                                        self.dialog = ActiveDialog::InsertObject;
                                    }
                                    ui.close_menu();
                                }
                                if ui.button("Duplicate Last Object").clicked() {
                                    self.doc_mut().apply_scene_tx("Duplicate Object", |s, nid| {
                                        if let Some(mut g) = s.graphics.last().cloned() {
                                            g.id = *nid; *nid += 1;
                                            g.bbox = (g.bbox.0 + 14.0, g.bbox.1 + 14.0, g.bbox.2 + 14.0, g.bbox.3 + 14.0);
                                            s.graphics.push(g);
                                        } else if let Some(mut t) = s.tables.last().cloned() {
                                            t.id = *nid; *nid += 1;
                                            t.bbox = (t.bbox.0 + 14.0, t.bbox.1 + 14.0, t.bbox.2 + 14.0, t.bbox.3 + 14.0);
                                            s.tables.push(t);
                                        } else if let Some(mut p) = s.tlc_plates.last().cloned() {
                                            p.id = *nid; *nid += 1;
                                            p.bbox = (p.bbox.0 + 14.0, p.bbox.1 + 14.0, p.bbox.2 + 14.0, p.bbox.3 + 14.0);
                                            s.tlc_plates.push(p);
                                        } else if let Some(mut ep) = s.energy_profiles.last().cloned() {
                                            ep.id = *nid; *nid += 1;
                                            ep.bbox = (ep.bbox.0 + 14.0, ep.bbox.1 + 14.0, ep.bbox.2 + 14.0, ep.bbox.3 + 14.0);
                                            s.energy_profiles.push(ep);
                                        } else if let Some(mut c) = s.curves.last().cloned() {
                                            c.id = *nid; *nid += 1;
                                            for pt in &mut c.points { pt.0 += 14.0; pt.1 += 14.0; }
                                            s.curves.push(c);
                                        }
                                    });
                                    self.status_msg = "Duplicated last embedded object".to_string();
                                    ui.close_menu();
                                }
                                if ui.button("Delete Last Object").clicked() {
                                    self.doc_mut().apply_scene_tx("Delete Object", |s, _| {
                                        if s.graphics.pop().is_none()
                                            && s.tables.pop().is_none()
                                            && s.tlc_plates.pop().is_none()
                                            && s.energy_profiles.pop().is_none()
                                            && s.curves.pop().is_none()
                                        {
                                            let _ = s.captions.pop();
                                        }
                                    });
                                    self.status_msg = "Deleted last embedded object".to_string();
                                    ui.close_menu();
                                }
                            });
                        } else {
                            let _ = ui.add_enabled(false, egui::Button::new("Object"));
                        }
                        }); }
                    });

                    ui.menu_button("View", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        ui.checkbox(&mut self.show_crosshair, "Show Crosshair\tCtrl+H");
                        ui.checkbox(&mut self.show_rulers, "Show Rulers\tF11");
                        ui.separator();
                        ui.checkbox(&mut self.show_tools, "Show Main Toolbar");
                        ui.checkbox(&mut self.show_helm_toolbar, "Show HELM Monomer Toolbar");
                        ui.checkbox(&mut self.show_general_toolbar, "Show General Toolbar");
                        ui.checkbox(&mut self.show_style_toolbar, "Show Style Toolbar");
                        ui.checkbox(&mut self.show_object_toolbar, "Show Object Toolbar");
                        ui.checkbox(&mut self.show_structure_toolbar, "Show Structure Toolbar");
                        ui.checkbox(&mut self.show_ring_fill_toolbar, "Show Ring Fill Colors Toolbar");
                        ui.checkbox(&mut self.show_biographics_toolbar, "Show BioGraphics Toolbar");
                        ui.checkbox(&mut self.show_analysis, "Show Analysis Window");
                        ui.checkbox(&mut self.show_properties, "Show Chemical Properties Window");
                        ui.checkbox(&mut self.show_stoichiometry, "Show Stoichiometry Grid");
                        ui.checkbox(&mut self.show_object_inspector, "Show Object Inspector (TLC / Tables / 3D / Export)");
                        ui.checkbox(&mut self.viewer_3d.open, "Show 3D Molecular & Crystal Viewer");
                        if ui.checkbox(&mut self.show_right_dock, "Show Right-Side Studio Dock").changed() {
                            self.persist_user_config();
                        }
                        ui.checkbox(&mut self.pro_suite.show_pub_composer, "Show Publication TOC & Scheme Composer");
                        ui.checkbox(&mut self.omni_suite.show_spectroscopy, "Show Multi-Spectroscopy Studio (IR / UV / 2D NMR)");
                        ui.checkbox(&mut self.omni_suite.show_energy_profile, "Show Reaction Energy Profile Builder");
                        ui.checkbox(&mut self.omni_suite.show_eln_safety, "Show ELN & Lab Safety / GHS Handbook");
                        if ui.button("Command Palette...\t: / F4").clicked() {
                            self.cmd_palette.open_fresh();
                            ui.close_menu();
                        }
                        ui.checkbox(&mut self.show_info_window, "Show Info Window");
                        ui.checkbox(&mut self.show_periodic_table, "Show Periodic Table Window");
                        ui.checkbox(&mut self.show_char_map, "Show Character Map Window");
                        ui.menu_button("Appearance", |ui| {
                            for (lbl, th) in [
                                ("Classic Light (White Canvas)", UiTheme::ClassicLight),
                                ("Dark UI + White Canvas", UiTheme::DarkChrome),
                                ("Dark UI + Dark Slide Canvas", UiTheme::DarkSlide),
                            ] {
                                if ui.selectable_label(self.ui_theme == th, lbl).clicked() {
                                    self.ui_theme = th;
                                    Self::apply_theme(ctx, th);
                                    self.persist_user_config();
                                    self.status_msg = format!("Switched appearance to {lbl}");
                                    ui.close_menu();
                                }
                            }
                        });
                        ui.menu_button("Other Toolbars", |ui| {
                            for (lbl, fly) in [
                                ("Accessories", ActiveFlyout::Accessories),
                                ("Arrows", ActiveFlyout::Arrows),
                                ("Brackets", ActiveFlyout::Brackets),
                                ("Chains", ActiveFlyout::Chains),
                                ("Chemical Symbols", ActiveFlyout::ChemSymbols),
                                ("Chromatography Tools", ActiveFlyout::Chromatography),
                                ("Colors", ActiveFlyout::Colors),
                                ("Curves", ActiveFlyout::Curves),
                                ("Dotted/Dashed Bond Tools", ActiveFlyout::DottedDashed),
                                ("Drawing Elements", ActiveFlyout::Shapes),
                                ("Fragmentation Tools", ActiveFlyout::Fragmentation),
                                ("Multiple Bonds", ActiveFlyout::MultipleBonds),
                                ("Orbitals", ActiveFlyout::Orbitals),
                                ("Pen Tools", ActiveFlyout::PenTools),
                                ("Windows", ActiveFlyout::Windows),
                            ] {
                                let checked = self.active_flyout == fly;
                                if ui.selectable_label(checked, lbl).clicked() {
                                    self.active_flyout = if checked { ActiveFlyout::None } else { fly };
                                    ui.close_menu();
                                }
                            }
                        });
                        ui.menu_button("Templates", |ui| {
                            let names: Vec<String> = self.ctp_libraries.iter().map(|l| l.name.clone()).collect();
                            egui::ScrollArea::vertical().max_height(420.0).show(ui, |ui| {
                                for (i, n) in names.iter().enumerate() {
                                    let active = self.active_flyout == ActiveFlyout::Templates && self.selected_ctp_idx == i;
                                    if ui.selectable_label(active, n).clicked() {
                                        self.selected_ctp_idx = i;
                                        self.active_flyout = ActiveFlyout::Templates;
                                        ui.close_menu();
                                    }
                                }
                                ui.separator();
                                if ui.button("New Templates").clicked() {
                                    let custom_idx = if let Some(idx) = self.ctp_libraries.iter().position(|l| l.name == "New Templates") {
                                        idx
                                    } else {
                                        self.ctp_libraries.push(TemplateLibrary {
                                            name: "New Templates".to_string(),
                                            templates: self.doc().scene.molecules.first().cloned().into_iter().collect(),
                                            graphic_templates: Vec::new(),
                                        });
                                        self.ctp_libraries.len() - 1
                                    };
                                    self.selected_ctp_idx = custom_idx;
                                    self.active_flyout = ActiveFlyout::Templates;
                                    ui.close_menu();
                                }
                            });
                        });
                        ui.separator();
                        ui.checkbox(&mut self.show_warnings, "Show Chemical Warnings");
                        ui.checkbox(&mut self.show_reaction_interp, "Show Reaction Interpretation");
                        ui.checkbox(&mut self.tab_mut().show_atom_numbers, "Show Atom Numbers");
                        ui.separator();
                        if ui.add_enabled((self.zoom - 1.0).abs() > 0.01, egui::Button::new("Actual Size\tF5")).clicked() {
                            self.zoom = 1.0;
                            ui.close_menu();
                        }
                        if ui.button("Show Document\tF6").clicked() {
                            let ws = self.tab().win_size;
                            let (pw, ph) = self.tab().effective_page_size_pt();
                            self.zoom = ((ws.x - 30.0) / (pw as f32)).min((ws.y - 45.0) / (ph as f32)).clamp(0.4, 4.0);
                            self.tab_mut().scroll = Vec2::ZERO;
                            ui.close_menu();
                        }
                        if ui.button("Magnify\tF7").clicked() { self.zoom = (self.zoom * 1.25).min(6.0); ui.close_menu(); }
                        if ui.button("Reduce\tF8").clicked() { self.zoom = (self.zoom * 0.8).max(0.4); ui.close_menu(); }
                        }); }
                    });

                    ui.menu_button("Object", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        let has_mols = !self.doc().scene.molecules.is_empty();
                        let can_group = self.tab().selected_mols.len() >= 2;
                        let can_ungroup = self.tab().groups.iter().any(|g| g.iter().any(|mi| self.tab().selected_mols.contains(mi)));
                        if ui.button("Object Settings...").clicked() { self.dialog = ActiveDialog::ObjectSettings; ui.close_menu(); }
                        ui.menu_button("Apply Object Settings from", |ui| {
                            if ui.button("Other...").clicked() {
                                self.dialog = ActiveDialog::ApplyStylesheetFile;
                                ui.close_menu();
                            }
                            if ui.button("Current Document").clicked() {
                                self.do_clean_up_structure();
                                self.status_msg = "Applied Current Document object settings".to_string();
                                ui.close_menu();
                            }
                            ui.separator();
                            let styles = self.stylesheets.clone();
                            for (lbl, settings) in styles {
                                if ui.button(&lbl).clicked() {
                                    replace_settings_tx(self.doc_mut(), settings, "Apply Object Settings");
                                    self.do_clean_up_structure();
                                    self.status_msg = format!("Applied object settings from {lbl}");
                                    ui.close_menu();
                                }
                            }
                        });
                        ui.checkbox(&mut self.doc_mut().settings.fixed_lengths, "Fixed Lengths\tCtrl+L");
                        ui.checkbox(&mut self.doc_mut().settings.fixed_angles, "Fixed Angles\tCtrl+E");
                        ui.checkbox(&mut self.show_stereo, "Show Stereochemistry");
                        ui.separator();
                        if ui.button("Attach Data...").clicked() { self.dialog = ActiveDialog::AttachData; ui.close_menu(); }
                        if ui.button("Annotate...").clicked() { self.dialog = ActiveDialog::Annotate; ui.close_menu(); }
                        ui.separator();
                        if ui.add_enabled(has_mols, egui::Button::new("Center on Page")).clicked() {
                            self.do_center_on_page_action();
                            ui.close_menu();
                        }
                        ui.menu_button("Align", |ui| {
                            for (lbl, mode) in [
                                ("Left edges", AlignMode::Left),
                                ("L/R centers", AlignMode::HorizontalCenter),
                                ("Right edges", AlignMode::Right),
                            ] {
                                if ui.button(lbl).clicked() {
                                    self.do_align_action(lbl, mode);
                                    ui.close_menu();
                                }
                            }
                            ui.separator();
                            for (lbl, mode) in [
                                ("Top edges", AlignMode::Top),
                                ("T/B centers", AlignMode::VerticalCenter),
                                ("Bottom edges", AlignMode::Bottom),
                            ] {
                                if ui.button(lbl).clicked() {
                                    self.do_align_action(lbl, mode);
                                    ui.close_menu();
                                }
                            }
                        });
                        ui.menu_button("Distribute", |ui| {
                            if ui.button("Horizontally").clicked() {
                                self.do_distribute_action(DistributeAxis::Horizontal);
                                ui.close_menu();
                            }
                            if ui.button("Vertically").clicked() {
                                self.do_distribute_action(DistributeAxis::Vertical);
                                ui.close_menu();
                            }
                        });
                        ui.menu_button("Add Frame", |ui| {
                            for (lbl, gtype, stype) in [
                                ("Brackets", "Bracket", Some("Square")),
                                ("Parentheses", "Bracket", Some("Round")),
                                ("Braces", "Bracket", Some("Curly")),
                                ("Rectangle", "Rectangle", Some("Plain")),
                                ("Rounded Rectangle", "Rectangle", Some("Rounded")),
                                ("Rectangle with Shadow", "Rectangle", Some("Shadow")),
                                ("Rounded Rectangle with Shadow", "Rectangle", Some("RoundedShadow")),
                            ] {
                                if ui.button(lbl).clicked() {
                                    let (x1, y1, x2, y2) = self.selection_bbox().unwrap_or((120.0, 120.0, 260.0, 210.0));
                                    self.doc_mut().apply_scene_tx("Add Frame", |s, nid| {
                                        let id = *nid;
                                        *nid += 1;
                                        s.graphics.push(SceneGraphic {
                                            id,
                                            graphic_type: gtype.to_string(),
                                            sub_type: stype.map(|v| v.to_string()),
                                            bbox: (x1 - 10.0, y1 - 10.0, x2 + 10.0, y2 + 10.0),
                                        });
                                    });
                                    self.status_msg = format!("Added {lbl} frame");
                                    ui.close_menu();
                                }
                            }
                        });
                        ui.separator();
                        if ui.add_enabled(can_group, egui::Button::new("Group\tCtrl+G")).clicked() {
                            let sel = self.tab().selected_mols.clone();
                            self.tab_mut().groups.push(sel);
                            self.status_msg = "Grouped selected objects".to_string();
                            ui.close_menu();
                        }
                        if ui.add_enabled(can_ungroup, egui::Button::new("Ungroup\tShift+Ctrl+G")).clicked() {
                            let sel = self.tab().selected_mols.clone();
                            self.tab_mut().groups.retain(|g| !g.iter().any(|mi| sel.contains(mi)));
                            self.status_msg = "Ungrouped".to_string();
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Join\tCtrl+J")).clicked() {
                            self.join_close_atoms_tx();
                            self.last_repeat_action = Some("Join".to_string());
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.add_enabled(has_mols, egui::Button::new("Bring to Front\tF2")).clicked() {
                            if let Some(mi) = self.active_mol_index() {
                                let mut moved_to = None;
                                self.doc_mut().apply_scene_tx("Bring to Front", |s, _| {
                                    if mi < s.molecules.len() {
                                        let m = s.molecules.remove(mi);
                                        s.molecules.push(m);
                                        moved_to = Some(s.molecules.len() - 1);
                                    }
                                });
                                if let Some(ni) = moved_to {
                                    self.tab_mut().remap_after_mol_move(mi, ni);
                                    self.tab_mut().selected_mols = vec![ni];
                                }
                            }
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Send to Back\tF3")).clicked() {
                            if let Some(mi) = self.active_mol_index() {
                                let mut moved = false;
                                self.doc_mut().apply_scene_tx("Send to Back", |s, _| {
                                    if mi < s.molecules.len() {
                                        let m = s.molecules.remove(mi);
                                        s.molecules.insert(0, m);
                                        moved = true;
                                    }
                                });
                                if moved {
                                    self.tab_mut().remap_after_mol_move(mi, 0);
                                    self.tab_mut().selected_mols = vec![0];
                                }
                            }
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Flatten")).clicked() {
                            self.do_flatten_action();
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.add_enabled(has_mols, egui::Button::new("Flip Horizontal")).clicked() {
                            self.do_flip_or_rot180("Flip Horizontal", FlipAxis::Horizontal, false);
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Flip Vertical")).clicked() {
                            self.do_flip_or_rot180("Flip Vertical", FlipAxis::Vertical, false);
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Rotate 180° Horizontal\tShift+Ctrl+H")).clicked() {
                            self.do_flip_or_rot180("Rotate 180° Horizontal", FlipAxis::Horizontal, true);
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Rotate 180° Vertical\tShift+Ctrl+V")).clicked() {
                            self.do_flip_or_rot180("Rotate 180° Vertical", FlipAxis::Vertical, true);
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Rotate...\tCtrl+R")).clicked() { self.dialog = ActiveDialog::RotateDialog; ui.close_menu(); }
                        if ui.add_enabled(has_mols, egui::Button::new("Scale...\tCtrl+K")).clicked() { self.dialog = ActiveDialog::ScaleDialog; ui.close_menu(); }
                        }); }
                    });

                    ui.menu_button("Structure", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        let has_mols = !self.doc().scene.molecules.is_empty();
                        if ui.add_enabled(has_mols, egui::Button::new("Atom Properties...")).clicked() {
                            if let Some((mi, ai)) = self.hovered_atom {
                                self.dialog = ActiveDialog::AtomProperties { mol_idx: mi, atom_idx: ai };
                            } else if let Some(mi) = self.active_mol_index() {
                                if !self.doc().scene.molecules[mi].atoms.is_empty() {
                                    self.dialog = ActiveDialog::AtomProperties { mol_idx: mi, atom_idx: 0 };
                                }
                            }
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Bond Properties...")).clicked() {
                            if let Some((mi, bi)) = self.hovered_bond {
                                self.dialog = ActiveDialog::BondProperties { mol_idx: mi, bond_idx: bi };
                            } else if let Some(mi) = self.active_mol_index() {
                                if !self.doc().scene.molecules[mi].bonds.is_empty() {
                                    self.dialog = ActiveDialog::BondProperties { mol_idx: mi, bond_idx: 0 };
                                }
                            }
                            ui.close_menu();
                        }
                        if ui.button("Bracket Properties...").clicked() {
                            self.dialog = ActiveDialog::BracketProperties;
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.add_enabled(has_mols, egui::Button::new("Check Structure")).clicked() {
                            self.check_warn_idx = 0;
                            self.dialog = ActiveDialog::CheckStructureDialog;
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Clean Up Structure\tShift+Ctrl+K")).clicked() {
                            self.do_clean_up_structure();
                            ui.close_menu();
                        }
                        if ui.add_enabled(self.doc().scene.molecules.len() >= 2, egui::Button::new("Clean Up Reaction\tShift+Ctrl+X")).clicked() {
                            self.doc_mut().apply_scene_tx("Clean Up Reaction", |s, _| {
                                clean_up_reaction(s);
                            });
                            self.last_repeat_action = Some("Clean Up Reaction".to_string());
                            self.status_msg = "Cleaned up reaction alignment".to_string();
                            ui.close_menu();
                        }
                        ui.menu_button("Catalytic Cycle", |ui| {
                            if ui.add_enabled(self.doc().scene.molecules.len() >= 2, egui::Button::new("Arrange Selected as Circle")).clicked() {
                                self.arrange_selected_catalytic_cycle(false);
                                ui.close_menu();
                            }
                            if ui.add_enabled(self.doc().scene.molecules.len() >= 2, egui::Button::new("Arrange Selected as Ellipse")).clicked() {
                                self.arrange_selected_catalytic_cycle(true);
                                ui.close_menu();
                            }
                            ui.separator();
                            if ui.button("Insert 4-Step Pd Catalytic Cycle").clicked() {
                                self.doc_mut().apply_scene_tx("Insert Starter Catalytic Cycle", |s, nid| {
                                    insert_starter_catalytic_cycle(s, nid, 250.0, 210.0);
                                });
                                self.status_msg = "Inserted starter 4-step catalytic cycle".to_string();
                                ui.close_menu();
                            }
                        });
                        ui.menu_button("Insert Special Core", |ui| {
                            let blen = self.doc().settings.bond_length_pt;
                            if ui.button("Haworth β-D-Glucopyranose").clicked() {
                                self.insert_special_core("Haworth β-D-Glucopyranose", build_haworth_pyranose(blen));
                                ui.close_menu();
                            }
                            if ui.button("Haworth β-D-Ribofuranose").clicked() {
                                self.insert_special_core("Haworth β-D-Ribofuranose", build_haworth_furanose(blen));
                                ui.close_menu();
                            }
                            ui.separator();
                            if ui.button("Porphine Macrocycle").clicked() {
                                self.insert_special_core("Porphine Macrocycle", build_porphine_core(blen));
                                ui.close_menu();
                            }
                            if ui.button("18-Crown-6 Macrocycle").clicked() {
                                self.insert_special_core("18-Crown-6", build_crown_ether_18_6(blen));
                                ui.close_menu();
                            }
                            ui.separator();
                            if ui.button("Adamantane Bridged Cage").clicked() {
                                self.insert_special_core("Adamantane", build_adamantane_2d(blen));
                                ui.close_menu();
                            }
                            if ui.button("Norbornane Bridged Ring").clicked() {
                                self.insert_special_core("Norbornane", build_norbornane_2d(blen));
                                ui.close_menu();
                            }
                        });
                        if ui.add_enabled(has_mols, egui::Button::new("View Active Structure in 3D")).clicked() {
                            self.viewer_3d.open = true;
                            ui.close_menu();
                        }
                        if ui.button("Retrosynthetic Disconnection Advisor...").clicked() {
                            self.pro_suite.show_retro = true;
                            if let Some(m) = self.active_mol().cloned() {
                                let blen = self.doc().settings.bond_length_pt;
                                self.pro_suite.retro_hits = suggest_disconnections(&m, blen);
                                self.pro_suite.highlighted_retro_bond = self.pro_suite.retro_hits.first().map(|h| h.bond_idx);
                            }
                            ui.close_menu();
                        }
                        if ui.button("Multi-Spectroscopy Studio (FT-IR / UV–Vis / 2D NMR)...").clicked() {
                            self.omni_suite.show_spectroscopy = true;
                            ui.close_menu();
                        }
                        if ui.button("MedChem, ADMET & PAINS Scorecard...").clicked() {
                            self.pro_suite.show_medchem = true;
                            ui.close_menu();
                        }
                        if ui.button("Quantum Chemistry Bridge (ORCA / Gaussian / xTB)...").clicked() {
                            self.quantum.show_window = true;
                            ui.close_menu();
                        }
                        ui.menu_button("Advanced Chemistry Studios", |ui| {
                            if ui.button("HRMS Isotope Envelope & Adducts...").clicked() {
                                self.pro_suite.show_hrms = true;
                                ui.close_menu();
                            }
                            if ui.button("pH Microspecies, Tautomers, Stereoisomers & Resonance...").clicked() {
                                self.pro_suite.show_isomers = true;
                                ui.close_menu();
                            }
                            if ui.button("Markush R-Group Library & SAR Grid...").clicked() {
                                self.pro_suite.show_markush = true;
                                ui.close_menu();
                            }
                            ui.separator();
                            if ui.button("Organometallic 18e⁻ & Polymer / BigSMILES Studio...").clicked() {
                                self.omni_suite.show_ino_poly = true;
                                ui.close_menu();
                            }
                            if ui.button("Reaction Energy Profile (ΔG‡) & Eyring Kinetics...").clicked() {
                                self.omni_suite.show_energy_profile = true;
                                ui.close_menu();
                            }
                            if ui.button("ELN Supporting Info & GHS Safety Handbook...").clicked() {
                                self.omni_suite.show_eln_safety = true;
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("Optimize 3D Geometry (UFF)")).clicked() {
                                self.execute_palette_command(ctx, PaletteCommandAction::Optimize3dNow);
                                ui.close_menu();
                            }
                        });
                        ui.menu_button("NMR Shift Labels", |ui| {
                            if ui.add_enabled(has_mols, egui::Button::new("Stamp ¹H-NMR δ Labels on Canvas")).clicked() {
                                self.stamp_nmr_shift_labels(false);
                                ui.close_menu();
                            }
                            if ui.add_enabled(has_mols, egui::Button::new("Stamp ¹³C-NMR δ Labels on Canvas")).clicked() {
                                self.stamp_nmr_shift_labels(true);
                                ui.close_menu();
                            }
                            if ui.button("Clear Stamped δ Labels").clicked() {
                                self.clear_nmr_shift_labels();
                                ui.close_menu();
                            }
                        });
                        if ui.add_enabled(has_mols, egui::Button::new("Expand Label")).clicked() {
                            let targets = self.target_mol_indices();
                            let jobs: Vec<(usize, usize, String)> = targets
                                .into_iter()
                                .filter_map(|mi| self.doc().scene.molecules.get(mi).map(|m| (mi, m)))
                                .flat_map(|(mi, m)| {
                                    m.atoms.iter().enumerate().filter_map(move |(ai, a)| {
                                        a.label.as_ref().map(|l| (mi, ai, l.clone()))
                                    })
                                })
                                .collect();
                            let mut expanded_cnt = 0usize;
                            for (mi, ai, lbl) in jobs {
                                if self.doc_mut().expand_label_tx(mi, ai, &lbl) {
                                    expanded_cnt += 1;
                                }
                            }
                            self.status_msg = if expanded_cnt > 0 {
                                format!("Expanded {expanded_cnt} nickname/group label(s)")
                            } else {
                                "No expandable nickname or functional-group labels found".to_string()
                            };
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Contract Label")).clicked() {
                            if let Some(mi) = self.active_mol_index() {
                                self.doc_mut().contract_last_residue_tx(mi);
                            }
                            ui.close_menu();
                        }
                        if ui.button("Expand Generic Structure").clicked() {
                            let cdxml = self.doc().to_cdxml();
                            let expanded = expand_generic_structure_from_cdxml(&cdxml);
                            if expanded.is_empty() {
                                self.status_msg = "No generic ElementList sites found to expand".to_string();
                            } else {
                                let n = expanded.len();
                                let mut expanded = expanded;
                                for (k, m) in expanded.iter_mut().enumerate() {
                                    for a in &mut m.atoms {
                                        a.x += ((k % 3) as f64) * 90.0;
                                        a.y += 95.0 + ((k / 3) as f64) * 80.0;
                                    }
                                }
                                stamp_molecules_tx(self.doc_mut(), expanded, "Expand Generic");
                                self.status_msg = format!("Expanded generic structure into {n} concrete product(s)");
                            }
                            ui.close_menu();
                        }
                        if ui.button("Expand Generic Structure to SDF...").clicked() {
                            let cdxml = self.doc().to_cdxml();
                            let expanded = expand_generic_structure_from_cdxml(&cdxml);
                            let mols = if expanded.is_empty() {
                                self.doc().scene.molecules.clone()
                            } else {
                                expanded
                            };
                            let recs: Vec<SdfRecord> = mols
                                .into_iter()
                                .enumerate()
                                .map(|(i, molecule)| SdfRecord {
                                    title: format!("{}_expanded_{}", self.doc().name, i + 1),
                                    molecule,
                                    properties: BTreeMap::new(),
                                })
                                .collect();
                            let sdf = molecules_to_sdf(&recs);
                            match native_dialogs::save_text_with_dialog(
                                "Export Expanded Generic Structures (SDF)",
                                "expanded_generic.sdf",
                                &sdf,
                            ) {
                                Ok(Some(out)) => {
                                    self.status_msg = format!("Exported {} expanded structure(s) to {out}", recs.len());
                                }
                                Ok(None) => {
                                    self.status_msg = "Expanded SDF export cancelled".to_string();
                                }
                                Err(e) => {
                                    self.status_msg = format!("SDF export failed: {e}");
                                }
                            }
                            ui.close_menu();
                        }
                        if ui.button("Toggle Aromatic Display\tAlt+K").clicked() {
                            self.toggle_aromatic_display_tx();
                            self.last_repeat_action = Some("Toggle Aromatic Display".to_string());
                            ui.close_menu();
                        }
                        if ui.button("Toggle Hydrogen Display").clicked() {
                            self.do_toggle_hydrogens();
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.button("Add Multi-Center Attachment").clicked() {
                            if let Some((x1, y1, x2, y2)) = self.selection_bbox() {
                                let cx = 0.5 * (x1 + x2);
                                let cy = 0.5 * (y1 + y2);
                                self.doc_mut().apply_scene_tx("Multi-Center Attachment", |s, nid| {
                                    let id = *nid;
                                    *nid += 1;
                                    s.graphics.push(SceneGraphic {
                                        id,
                                        graphic_type: "MultiCenter".to_string(),
                                        sub_type: None,
                                        bbox: (cx - 3.0, cy - 3.0, cx + 3.0, cy + 3.0),
                                    });
                                });
                                self.status_msg = "Added Multi-Center Attachment at ring centroid".to_string();
                            }
                            ui.close_menu();
                        }
                        if ui.button("Add Variable Attachment").clicked() {
                            if let Some((x1, y1, x2, y2)) = self.selection_bbox() {
                                let cx = 0.5 * (x1 + x2);
                                let cy = 0.5 * (y1 + y2);
                                self.doc_mut().apply_scene_tx("Variable Attachment", |s, nid| {
                                    let id = *nid;
                                    *nid += 1;
                                    s.graphics.push(SceneGraphic {
                                        id,
                                        graphic_type: "VariableAttachment".to_string(),
                                        sub_type: None,
                                        bbox: (cx, cy, cx + 22.0, cy - 16.0),
                                    });
                                });
                                self.status_msg = "Added Variable Attachment marker".to_string();
                            }
                            ui.close_menu();
                        }
                        if ui.button("R-Logic Query...").clicked() {
                            self.dialog = ActiveDialog::RLogicQuery;
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Add 3D Property")).clicked() {
                            if let Some(m) = self.active_mol() {
                                let pr = estimate_chemical_properties(m);
                                self.paste_caption_tx(format!(
                                    "3D Molecular Volume: {:.1} Å³  |  tPSA: {:.2} Å²  |  MR: {:.2} cm³/mol",
                                    pr.molar_refractivity * 2.85,
                                    pr.tpsa,
                                    pr.molar_refractivity
                                ));
                                self.status_msg = "Computed and attached 3D Molecular Properties".to_string();
                            }
                            ui.close_menu();
                        }
                        ui.menu_button("Enhanced Stereochemistry", |ui| {
                            if ui.button("Clear All").clicked() {
                                self.tab_mut().enhanced_stereo.clear();
                                self.status_msg = "Cleared Enhanced Stereochemistry groups".to_string();
                                ui.close_menu();
                            }
                            if ui.button("Mark as Absolute (abs)").clicked() {
                                self.do_mark_enhanced_stereo(EnhancedStereoKind::Absolute);
                                ui.close_menu();
                            }
                            if ui.button("Create OR Group (or1)").clicked() {
                                let n = self.tab().next_or_group;
                                self.tab_mut().next_or_group += 1;
                                self.do_mark_enhanced_stereo(EnhancedStereoKind::Or(n));
                                ui.close_menu();
                            }
                            if ui.button("Create AND Group (&1)").clicked() {
                                let n = self.tab().next_and_group;
                                self.tab_mut().next_and_group += 1;
                                self.do_mark_enhanced_stereo(EnhancedStereoKind::And(n));
                                ui.close_menu();
                            }
                        });
                        if ui.add_enabled(has_mols, egui::Button::new("Clean Up Biopolymer")).clicked() {
                            self.handle_pro_panel_action(ctx, ProPanelAction::CleanUpBiopolymer);
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.button("Map Reaction Atoms").clicked() {
                            let maps = auto_map_reaction_atoms(&self.doc().scene);
                            let pairs = maps.len() / 2;
                            self.tab_mut().reaction_atom_maps = maps;
                            self.status_msg = format!("Mapped {pairs} reaction atom pair(s)");
                            ui.close_menu();
                        }
                        if ui.button("Clear Reaction Map").clicked() {
                            self.tab_mut().reaction_atom_maps.clear();
                            self.status_msg = "Cleared reaction atom map".to_string();
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Autonumber Reaction")).clicked() {
                            self.do_autonumber_reaction();
                            ui.close_menu();
                        }
                        if ui.add_enabled(!self.tab().reaction_auto_numbers.is_empty(), egui::Button::new("Clear Reaction Numbers")).clicked() {
                            self.tab_mut().reaction_auto_numbers.clear();
                            self.status_msg = "Cleared reaction compound numbers".to_string();
                            ui.close_menu();
                        }
                        if ui.button("Stoichiometry Grid...").clicked() {
                            self.do_open_stoichiometry_grid();
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.add_enabled(has_mols, egui::Button::new("Predict ¹H-NMR Shifts")).clicked() {
                            self.do_predict_nmr(false);
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Predict ¹³C-NMR Shifts")).clicked() {
                            self.do_predict_nmr(true);
                            ui.close_menu();
                        }
                        if ui.button("Make Spectrum-Structure Assignment").clicked() {
                            self.tab_mut().show_atom_numbers = true;
                            let n_atoms: usize = self.doc().scene.molecules.iter().map(|m| m.atoms.len()).sum();
                            self.tab_mut().annotations.insert(
                                "SpectrumAssignment".to_string(),
                                format!("{n_atoms} atom(s) linked"),
                            );
                            if self.nmr_prediction.is_none() && has_mols {
                                self.do_predict_nmr(false);
                            }
                            self.status_msg = format!("Linked {n_atoms} structure atom(s) for spectrum-structure assignment");
                            ui.close_menu();
                        }
                        if ui.button("Clear Spectrum-Structure Assignment").clicked() {
                            self.nmr_selected_atom = None;
                            self.tab_mut().annotations.remove("SpectrumAssignment");
                            self.status_msg = "Cleared spectrum-structure assignment".to_string();
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.add_enabled(has_mols, egui::Button::new("Add Structure to Dictionary...")).clicked() {
                            self.dict_name_input.clear();
                            self.dialog = ActiveDialog::AddToNameDict;
                            ui.close_menu();
                        }
                        if ui.button("Define Nickname...").clicked() {
                            self.nick_pending_collision = None;
                            self.nick_confirm_replace = false;
                            self.nick_confirm_override = false;
                            self.dialog = ActiveDialog::DefineNickname;
                            ui.close_menu();
                        }
                        if ui.button("Convert Name to Structure...\tShift+Ctrl+N").clicked() {
                            self.dialog = ActiveDialog::PasteNameAsStructure;
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Convert Structure to Name\tAlt+Ctrl+N")).clicked() {
                            self.do_structure_to_name_hybrid(ctx, false);
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mols, egui::Button::new("Resolve Name via PubChem")).clicked() {
                            self.do_structure_to_name_hybrid(ctx, true);
                            ui.close_menu();
                        }
                        }); }
                    });

                    ui.menu_button("Text", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        ui.menu_button("Font", |ui| {
                            for f in ["Arial", "Helvetica", "Symbol", "Times"] {
                                ui.selectable_value(&mut self.doc_mut().settings.label_font, f.to_string(), f);
                            }
                            ui.separator();
                            for f in ["Calibri", "Courier New", "Georgia", "Times New Roman", "Verdana"] {
                                ui.selectable_value(&mut self.doc_mut().settings.label_font, f.to_string(), f);
                            }
                        });
                        ui.menu_button("Style", |ui| {
                            if ui.button("Plain").clicked() {
                                self.text_bold = false;
                                self.text_italic = false;
                                self.text_underline = false;
                                self.text_sub = false;
                                self.text_sup = false;
                                self.text_formula = false;
                                ui.close_menu();
                            }
                            ui.checkbox(&mut self.text_bold, "Bold\tCtrl+B");
                            ui.checkbox(&mut self.text_italic, "Italic\tCtrl+I");
                            ui.checkbox(&mut self.text_underline, "Underline\tCtrl+U");
                            ui.separator();
                            if ui.checkbox(&mut self.text_sub, "Subscript\tF9").changed() && self.text_sub {
                                self.text_sup = false;
                            }
                            if ui.checkbox(&mut self.text_sup, "Superscript\tF10").changed() && self.text_sup {
                                self.text_sub = false;
                            }
                            ui.checkbox(&mut self.text_formula, "Formula\tShift+F9");
                        });
                        ui.menu_button("Size", |ui| {
                            for sz in [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 24.0] {
                                if ui.selectable_value(&mut self.doc_mut().settings.label_size_pt, sz, format!("{sz:.0}")).clicked() {
                                    self.doc_mut().settings.caption_size_pt = sz;
                                }
                            }
                            ui.separator();
                            if ui.button("Other...").clicked() {
                                self.dialog = ActiveDialog::FontSizeOther;
                                ui.close_menu();
                            }
                        });
                        ui.separator();
                        ui.selectable_value(&mut self.text_align, TextAlign::Left, "Flush Left\tShift+Ctrl+L");
                        ui.selectable_value(&mut self.text_align, TextAlign::Center, "Centered\tShift+Ctrl+C");
                        ui.selectable_value(&mut self.text_align, TextAlign::Right, "Flush Right\tShift+Ctrl+R");
                        ui.selectable_value(&mut self.text_align, TextAlign::Justify, "Justified\tShift+Ctrl+J");
                        ui.selectable_value(&mut self.text_align, TextAlign::StackedAbove, "Stacked Above");
                        ui.selectable_value(&mut self.text_align, TextAlign::Automatic, "Automatic\tShift+Ctrl+M");
                        ui.separator();
                        if ui.button("Line Spacing...").clicked() {
                            self.dialog = ActiveDialog::LineSpacing;
                            ui.close_menu();
                        }
                        }); }
                    });

                    ui.menu_button("Curves", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        if ui.button("Plain").clicked() {
                            self.curve_dashed = false;
                            self.curve_bold = false;
                            self.curve_doubled = false;
                            self.curve_filled = false;
                            self.curve_shaded = false;
                            self.curve_faded = false;
                            self.curve_arrow_start = false;
                            self.curve_arrow_end = false;
                            self.curve_half_start = 0;
                            self.curve_half_end = 0;
                            ui.close_menu();
                        }
                        ui.checkbox(&mut self.curve_dashed, "Dashed");
                        ui.checkbox(&mut self.curve_bold, "Bold");
                        ui.checkbox(&mut self.curve_doubled, "Doubled");
                        ui.checkbox(&mut self.curve_arrow_start, "Full Arrow at Start");
                        ui.checkbox(&mut self.curve_arrow_end, "Full Arrow at End");
                        let mut hsl = self.curve_half_start == -1;
                        if ui.checkbox(&mut hsl, "Half Arrow at Start Left").changed() { self.curve_half_start = if hsl { -1 } else { 0 }; }
                        let mut hsr = self.curve_half_start == 1;
                        if ui.checkbox(&mut hsr, "Half Arrow at Start Right").changed() { self.curve_half_start = if hsr { 1 } else { 0 }; }
                        let mut hel = self.curve_half_end == -1;
                        if ui.checkbox(&mut hel, "Half Arrow at End Left").changed() { self.curve_half_end = if hel { -1 } else { 0 }; }
                        let mut her = self.curve_half_end == 1;
                        if ui.checkbox(&mut her, "Half Arrow at End Right").changed() { self.curve_half_end = if her { 1 } else { 0 }; }
                        if ui.button("Closed").clicked() {
                            self.doc_mut().apply_scene_tx("Close Curve", |s, _| {
                                if let Some(c) = s.curves.last_mut() { c.closed = !c.closed; }
                            });
                            ui.close_menu();
                        }
                        ui.checkbox(&mut self.curve_filled, "Filled");
                        ui.checkbox(&mut self.curve_shaded, "Shaded");
                        ui.checkbox(&mut self.curve_faded, "Faded");
                        }); }
                    });

                    ui.menu_button("Colors", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        for (name, col) in [
                            ("Foreground", Color32::BLACK),
                            ("Other # 1", Color32::from_rgb(255, 0, 0)),
                            ("Other # 2", Color32::from_rgb(240, 220, 0)),
                            ("Other # 3", Color32::from_rgb(0, 210, 0)),
                            ("Other # 4", Color32::from_rgb(0, 220, 235)),
                            ("Other # 5", Color32::from_rgb(0, 0, 255)),
                            ("Other # 6", Color32::from_rgb(255, 0, 255)),
                        ] {
                            ui.horizontal(|ui| {
                                let (r, _) = ui.allocate_exact_size(Vec2::splat(13.0), Sense::hover());
                                ui.painter().rect_filled(r, 0.0, col);
                                ui.painter().rect_stroke(r, 0.0, Stroke::new(1.0_f32, Color32::GRAY), egui::StrokeKind::Inside);
                                if ui.button(name).clicked() {
                                    self.color_by_element = false;
                                    self.active_color = col;
                                    ui.close_menu();
                                }
                            });
                        }
                        if ui.button("Other...").clicked() {
                            self.dialog = ActiveDialog::ColorPicker;
                            ui.close_menu();
                        }
                        ui.separator();
                        ui.checkbox(&mut self.color_by_element, "By Element");
                        ui.separator();
                        ui.checkbox(&mut self.show_ring_fill_toolbar, "Show Ring Fill & Highlight Colors");
                        }); }
                    });

                    ui.menu_button("Search", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        let has_mol = self.active_mol().is_some();
                        if ui.button("Search PubChem / ChEMBL...   Shift+Ctrl+F").clicked() {
                            self.show_chem_db = true;
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mol, egui::Button::new("Lookup Selected Structure (Exact)")).clicked() {
                            self.do_search_canvas_in_db(ctx, DbSearchMode::ExactSmiles);
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mol, egui::Button::new("Similarity Search Selected Structure...")).clicked() {
                            self.do_search_canvas_in_db(ctx, DbSearchMode::Similarity2D);
                            ui.close_menu();
                        }
                        if ui.add_enabled(has_mol, egui::Button::new("Substructure Search Selected Structure...")).clicked() {
                            self.do_search_canvas_in_db(ctx, DbSearchMode::Substructure);
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.button("Copy PubChem Portal URL").clicked() {
                            ctx.copy_text("https://pubchem.ncbi.nlm.nih.gov".to_string());
                            self.status_msg = "Copied https://pubchem.ncbi.nlm.nih.gov to clipboard".to_string();
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.button("Find Nickname...").clicked() {
                            self.dialog = ActiveDialog::ListNicknames;
                            ui.close_menu();
                        }
                        if ui.button("Find Element on Canvas...").clicked() {
                            self.dialog = ActiveDialog::FindElement;
                            ui.close_menu();
                        }
                        }); }
                    });

                    ui.menu_button("Window", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        if ui.button("Cascade").clicked() {
                            for (i, t) in self.tabs.iter_mut().enumerate() {
                                let c = ((i % 8) as f32) * 26.0;
                                t.win_pos = Vec2::new(24.0 + c, 18.0 + c);
                                t.win_size = Vec2::new(740.0, 600.0);
                                t.maximized = false;
                                t.minimized = false;
                            }
                            ui.close_menu();
                        }
                        if ui.button("Tile Vertical").clicked() {
                            let n = self.tabs.len().max(1) as f32;
                            let avail_w = (self.last_desk_size.x - 16.0 - 6.0 * (n - 1.0)).max(300.0);
                            let avail_h = (self.last_desk_size.y - 16.0).max(200.0);
                            let w = (avail_w / n).max(260.0);
                            for (i, t) in self.tabs.iter_mut().enumerate() {
                                t.win_pos = Vec2::new(8.0 + (i as f32) * (w + 6.0), 8.0);
                                t.win_size = Vec2::new(w, avail_h);
                                t.maximized = false;
                                t.minimized = false;
                            }
                            ui.close_menu();
                        }
                        if ui.button("Tile Horizontal").clicked() {
                            let n = self.tabs.len().max(1) as f32;
                            let avail_w = (self.last_desk_size.x - 16.0).max(300.0);
                            let avail_h = (self.last_desk_size.y - 16.0 - 6.0 * (n - 1.0)).max(200.0);
                            let h = (avail_h / n).max(180.0);
                            for (i, t) in self.tabs.iter_mut().enumerate() {
                                t.win_pos = Vec2::new(8.0, 8.0 + (i as f32) * (h + 6.0));
                                t.win_size = Vec2::new(avail_w, h);
                                t.maximized = false;
                                t.minimized = false;
                            }
                            ui.close_menu();
                        }
                        let any_min = self.tabs.iter().any(|t| t.minimized);
                        if ui.add_enabled(any_min, egui::Button::new("Arrange Icons")).clicked() {
                            for (i, t) in self.tabs.iter_mut().enumerate() {
                                if t.minimized {
                                    t.win_pos = Vec2::new(8.0 + (i as f32) * 190.0, 640.0);
                                }
                            }
                            ui.close_menu();
                        }
                        ui.separator();
                        let tab_names: Vec<String> = self
                            .tabs
                            .iter()
                            .enumerate()
                            .map(|(idx, t)| format!("{} {}{}", idx + 1, t.doc.name, if t.doc.is_dirty() { " *" } else { "" }))
                            .collect();
                        for (idx, name) in tab_names.iter().enumerate() {
                            if ui.selectable_label(idx == self.active_tab, name).clicked() {
                                self.active_tab = idx;
                                self.tab_mut().minimized = false;
                                self.inline_atom_edit = None;
                                self.sanitize_state();
                                ui.close_menu();
                            }
                        }
                        }); }
                    });

                    ui.menu_button("Help", |ui| {
                        // Tall menus (Structure) overflow an 860 px window: keep every item reachable.
                        { let h = ui.ctx().screen_rect().height() - 48.0; egui::ScrollArea::vertical().max_height(h).min_scrolled_height(h).show(ui, |ui| {
                        if ui.button("Contents\tF1").clicked() {
                            self.dialog = ActiveDialog::Shortcuts;
                            ui.close_menu();
                        }
                        for (lbl, body) in [
                            ("Getting Started", "Draw bonds using the left Tools palette, hover atoms and press element hotkeys (o=O, n=N, f=F, C=Cl, b=Br, 2==O, A=Ac), press ':' or F4 for the Command Palette, or pick any of the 11 built-in Template libraries from View > Templates."),
                            ("Tutorials", "1. Click Benzene and click empty canvas to stamp.\n2. Hover an edge and click Benzene again to fuse Naphthalene.\n3. Hover an atom and press o, n, or Enter to type a nickname like Ph, Boc, O₂N, or CF₃.\n4. Use Structure > Clean Up Structure (Shift+Ctrl+K), Check Structure, or the Right-Side Studio Dock."),
                            ("ChiralCraft Basics", "ChiralCraft workspace: Drag MDI document title bars to move, drag edges/corners to resize, Ctrl+Wheel to zoom, right-click atoms/bonds for quick edits, and export to 22+ chemical, crystallographic, quantum, and lossless editable SVG/PNG/PDF formats."),
                            ("CLI & MCP Server", "Batch CLI & MCP server built directly into chiralcraft:\n  chiralcraft --help\n  chiralcraft analyze <file|SMILES>\n  chiralcraft convert <in> <out>\n  chiralcraft --mcp  (19-tool stdio JSON-RPC 2.0 MCP server)\nSupports CDX, CDXML, MOL, SDF, RXN, XYZ, CIF, PDB, POSCAR, JCAMP-DX, ORCA/Gaussian/xTB, and editable SVG/PNG/PDF."),
                            ("Recent Additions", "ChiralCraft highlights:\n• Live PubChem & ChEMBL search + offline/hybrid IUPAC Structure ↔ Name\n• ¹H/¹³C NMR + experimental JCAMP-DX overlay, FT-IR, UV–Vis & 2D NMR (COSY/HSQC/HMBC)\n• MedChem/ADMET/PAINS Scorecard, HRMS Isotopes, Isomers/pH, Markush SAR & Retrosynthesis\n• Interactive 3D Viewer, UFF-lite Optimizer, Crystal Unit Cells & ORCA/Gaussian/xTB Bridge\n• Organometallic 18e⁻ CBC, Polymers/BigSMILES, Reaction Energy Profiles (ΔG‡), LaTeX & ELN SI\n• Lossless Editable SVG/PNG/PDF round-trip + 19-tool MCP Server."),
                        ] {
                            if ui.button(lbl).clicked() {
                                self.dialog = ActiveDialog::HelpTopic {
                                    title: lbl.to_string(),
                                    body: body.to_string(),
                                };
                                ui.close_menu();
                            }
                        }
                        if ui.button("Hotkeys Cheat Sheet").clicked() {
                            self.dialog = ActiveDialog::Shortcuts;
                            ui.close_menu();
                        }
                        if ui.button("Shortcuts Cheat Sheet").clicked() {
                            self.dialog = ActiveDialog::Shortcuts;
                            ui.close_menu();
                        }
                        if ui.button("Diagnostics Info").clicked() {
                            let crash_log = crash_log_path();
                            self.dialog = ActiveDialog::HelpTopic {
                                title: "Diagnostics Info".to_string(),
                                body: format!(
                                    "ChiralCraft {} — Pure-Rust Chemical Structure Editor\nDiagnostics crash log: {}",
                                    env!("CARGO_PKG_VERSION"),
                                    crash_log.display()
                                ),
                            };
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.button("Licensing...").clicked() {
                            self.dialog = ActiveDialog::Licensing;
                            ui.close_menu();
                        }
                        if ui.button("About ChiralCraft...").clicked() {
                            self.dialog = ActiveDialog::About;
                            ui.close_menu();
                        }
                        }); }
                    });
                });
            });

        // ---- 2. Exact Unified General + Style Top Toolbar (`p08-opened.png`, 100% Vector Icons) ----
        if self.show_general_toolbar || self.show_style_toolbar {
            egui::TopBottomPanel::top("main_toolbar")
                .frame(egui::Frame::side_top_panel(&ctx.style()).fill(CD_WINDOW_BG).inner_margin(egui::Margin::symmetric(4, 2)))
                .show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        ui.spacing_mut().item_spacing.x = 2.0;
                        if self.show_general_toolbar {
                            if toolbar_icon_button(ui, ToolbarIcon::NewDoc, false, "New Document (Ctrl+N)").clicked() {
                                self.new_document_tab(DocumentSettings::acs_document_1996());
                            }
                            if toolbar_icon_button(ui, ToolbarIcon::Open, false, "Open...").clicked() {
                                self.dialog = ActiveDialog::OpenFile;
                            }
                            if toolbar_icon_button(ui, ToolbarIcon::Save, false, "Save (Ctrl+S)").clicked() {
                                self.do_save_action();
                            }
                            if toolbar_icon_button(ui, ToolbarIcon::Print, false, "Export Vector PDF / Print...").clicked() {
                                self.file_path_buf = "print_export.pdf".to_string();
                                self.dialog = ActiveDialog::SaveAs;
                            }
                            if toolbar_icon_button(ui, ToolbarIcon::Help, false, "Keyboard Shortcuts & Help").clicked() {
                                self.dialog = ActiveDialog::Shortcuts;
                            }
                            let u_lbl = self.doc().undo_menu_label();
                            let r_lbl = self.doc().redo_menu_label();
                            if toolbar_icon_button(ui, ToolbarIcon::Undo, false, &u_lbl).clicked() {
                                self.tab_mut().undo_redo(false);
                                self.sanitize_state();
                            }
                            if toolbar_icon_button(ui, ToolbarIcon::Redo, false, &r_lbl).clicked() {
                                self.tab_mut().undo_redo(true);
                                self.sanitize_state();
                            }
                            if toolbar_icon_button(ui, ToolbarIcon::Cut, false, "Cut").clicked() {
                                self.do_cut_action();
                            }
                            if toolbar_icon_button(ui, ToolbarIcon::Copy, false, "Copy").clicked() {
                                self.do_copy_action();
                            }
                            if toolbar_icon_button(ui, ToolbarIcon::Paste, false, "Paste (Ctrl+V)").clicked() {
                                self.do_paste_action();
                            }
                            ui.add_space(3.0);
                            egui::ComboBox::from_id_salt("zoom_combo")
                                .width(58.0)
                                .selected_text(format!("{:.0}%", self.zoom * 100.0))
                                .show_ui(ui, |ui| {
                                    for pct in [50, 75, 100, 150, 200, 300, 400] {
                                        ui.selectable_value(&mut self.zoom, (pct as f32) / 100.0, format!("{pct}%"));
                                    }
                                });
                            if toolbar_icon_button(ui, ToolbarIcon::ZoomIn, false, "Zoom In").clicked() {
                                self.zoom = (self.zoom * 1.25).min(6.0);
                            }
                            if toolbar_icon_button(ui, ToolbarIcon::ZoomOut, false, "Zoom Out").clicked() {
                                self.zoom = (self.zoom * 0.8).max(0.4);
                            }
                            if toolbar_icon_button(ui, ToolbarIcon::CleanUp, false, "Clean Up Structure (Shift+Ctrl+K)").clicked() {
                                self.do_clean_up_structure();
                            }
                        }
                        if self.show_style_toolbar {
                            ui.add_space(3.0);
                            let cur_font = self.doc().settings.label_font.clone();
                            egui::ComboBox::from_id_salt("font_combo")
                                .width(88.0)
                                .selected_text(&cur_font)
                                .show_ui(ui, |ui| {
                                    for f in ["Arial", "Helvetica", "Symbol", "Times", "Calibri", "Courier New", "Georgia", "Times New Roman", "Verdana"] {
                                        ui.selectable_value(&mut self.doc_mut().settings.label_font, f.to_string(), f);
                                    }
                                });
                            let cur_sz = self.doc().settings.label_size_pt;
                            egui::ComboBox::from_id_salt("size_combo")
                                .width(42.0)
                                .selected_text(format!("{cur_sz:.0}"))
                                .show_ui(ui, |ui| {
                                    for sz in [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 24.0] {
                                        if ui.selectable_value(&mut self.doc_mut().settings.label_size_pt, sz, format!("{sz:.0}")).clicked() {
                                            self.doc_mut().settings.caption_size_pt = sz;
                                        }
                                    }
                                });
                            ui.separator();
                            if toolbar_icon_button(ui, ToolbarIcon::AlignLeft, self.text_align == TextAlign::Left, "Flush Left").clicked() { self.text_align = TextAlign::Left; }
                            if toolbar_icon_button(ui, ToolbarIcon::AlignCenter, self.text_align == TextAlign::Center, "Centered").clicked() { self.text_align = TextAlign::Center; }
                            if toolbar_icon_button(ui, ToolbarIcon::AlignRight, self.text_align == TextAlign::Right, "Flush Right").clicked() { self.text_align = TextAlign::Right; }
                            if toolbar_icon_button(ui, ToolbarIcon::AlignJustify, self.text_align == TextAlign::Justify, "Justified").clicked() { self.text_align = TextAlign::Justify; }
                            ui.separator();
                            ui.toggle_value(&mut self.text_bold, egui::RichText::new("B").strong());
                            ui.toggle_value(&mut self.text_italic, egui::RichText::new("I").italics());
                            ui.toggle_value(&mut self.text_underline, egui::RichText::new("U").underline());
                            ui.separator();
                            ui.toggle_value(&mut self.text_formula, "CH₂").on_hover_text("Formula Subscripts");
                            if ui.toggle_value(&mut self.text_sub, "X₂").on_hover_text("Subscript").changed() && self.text_sub {
                                self.text_sup = false;
                            }
                            if ui.toggle_value(&mut self.text_sup, "X²").on_hover_text("Superscript").changed() && self.text_sup {
                                self.text_sub = false;
                            }
                            ui.separator();
                            let (crect, cresp) = ui.allocate_exact_size(Vec2::new(20.0, 18.0), Sense::click());
                            let inner = Rect::from_min_size(crect.min + Vec2::new(2.0, 2.0), Vec2::splat(13.0));
                            ui.painter().rect_filled(inner, 0.0, self.active_color);
                            ui.painter().rect_stroke(inner, 0.0, Stroke::new(1.0_f32, Color32::GRAY), egui::StrokeKind::Inside);
                            let br = crect.max - Vec2::new(1.0, 1.0);
                            ui.painter().add(egui::Shape::convex_polygon(
                                vec![br, br - Vec2::new(4.0, 0.0), br - Vec2::new(0.0, 4.0)],
                                Color32::BLACK,
                                Stroke::NONE,
                            ));
                            let cresp = cresp.on_hover_text("Click for custom color picker (Right-click for By Element CPK)");
                            if cresp.secondary_clicked() {
                                self.color_by_element = !self.color_by_element;
                            } else if cresp.clicked() {
                                self.dialog = ActiveDialog::ColorPicker;
                            }
                        }
                    });
                });
        }

        // ---- 2b. Optional Object Toolbar (`View > Show Object Toolbar`) ----
        if self.show_object_toolbar {
            egui::TopBottomPanel::top("object_toolbar")
                .frame(egui::Frame::side_top_panel(&ctx.style()).fill(CD_WINDOW_BG).inner_margin(egui::Margin::symmetric(6, 2)))
                .show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        if ui.small_button("Front").on_hover_text("Bring to Front (F2)").clicked() {
                            if let Some(mi) = self.active_mol_index() {
                                let mut new_idx = None;
                                self.doc_mut().apply_scene_tx("Bring to Front", |s, _| {
                                    if mi < s.molecules.len() {
                                        let m = s.molecules.remove(mi);
                                        s.molecules.push(m);
                                        new_idx = Some(s.molecules.len() - 1);
                                    }
                                });
                                if let Some(ni) = new_idx {
                                    self.tab_mut().selected_mols = vec![ni];
                                    self.status_msg = "Brought structure to front".to_string();
                                }
                            }
                        }
                        if ui.small_button("Back").on_hover_text("Send to Back (F3)").clicked() {
                            if let Some(mi) = self.active_mol_index() {
                                let mut moved = false;
                                self.doc_mut().apply_scene_tx("Send to Back", |s, _| {
                                    if mi < s.molecules.len() {
                                        let m = s.molecules.remove(mi);
                                        s.molecules.insert(0, m);
                                        moved = true;
                                    }
                                });
                                if moved {
                                    self.tab_mut().selected_mols = vec![0];
                                    self.status_msg = "Sent structure to back".to_string();
                                }
                            }
                        }
                        ui.separator();
                        if ui.small_button("Flip H").clicked() {
                            self.do_flip_or_rot180("Flip Horizontal", FlipAxis::Horizontal, false);
                        }
                        if ui.small_button("Flip V").clicked() {
                            self.do_flip_or_rot180("Flip Vertical", FlipAxis::Vertical, false);
                        }
                        if ui.small_button("Align L").clicked() { self.do_align_action("Left edges", AlignMode::Left); }
                        if ui.small_button("Align C").clicked() { self.do_align_action("L/R centers", AlignMode::HorizontalCenter); }
                        if ui.small_button("Align R").clicked() { self.do_align_action("Right edges", AlignMode::Right); }
                        if ui.small_button("Dist H").clicked() { self.do_distribute_action(DistributeAxis::Horizontal); }
                        if ui.small_button("Dist V").clicked() { self.do_distribute_action(DistributeAxis::Vertical); }
                        if ui.small_button("Center").clicked() { self.do_center_on_page_action(); }
                    });
                });
        }

        // ---- 2c. Optional ChiralCraft Structure Toolbar (`View > Show Structure Toolbar`) ----
        if self.show_structure_toolbar {
            let act = show_structure_toolbar(ctx);
            self.handle_pro_panel_action(ctx, act);
        }

        // ---- 3. Iconic 18x2 Left Tools Palette (`tools-20261118.png`) + ALL 15 `Other Toolbars` + 26 `.ctp` Templates ----
        match show_tools_palette(
            ctx,
            self.show_tools,
            &mut self.tool,
            &mut self.active_element_z,
            &mut self.active_flyout,
            &self.ctp_libraries,
            &mut self.selected_ctp_idx,
        ) {
            PaletteCommand::StampTemplate(lib_i, tpl_i) => {
                if let Some(lib) = self.ctp_libraries.get(lib_i) {
                    let lib_name = lib.name.clone();
                    if !lib.templates.is_empty() {
                        if let Some(m) = lib.stamp(tpl_i, 240.0, 200.0) {
                            self.doc_mut().stamp_molecule_tx(m, "Template Stamp");
                            self.status_msg = format!("Stamped {lib_name} #{}", tpl_i + 1);
                        }
                    } else if let Some(curves) = lib.stamp_graphic(tpl_i, 240.0, 200.0) {
                        self.doc_mut().apply_scene_tx("Graphic Template Stamp", |s, nid| {
                            for mut c in curves {
                                c.id = *nid;
                                *nid += 1;
                                s.curves.push(c);
                            }
                        });
                        self.status_msg = format!("Stamped {lib_name} graphic #{}", tpl_i + 1);
                    }
                }
            }
            PaletteCommand::SetColor(col) => {
                self.color_by_element = false;
                self.active_color = col;
                self.status_msg = "Updated active color".to_string();
            }
            PaletteCommand::ToggleColorByElement => {
                self.color_by_element = !self.color_by_element;
            }
            PaletteCommand::OpenColorPicker => {
                self.dialog = ActiveDialog::ColorPicker;
            }
            PaletteCommand::ToggleCurve(k) => match k {
                "plain" => {
                    self.curve_dashed = false;
                    self.curve_bold = false;
                    self.curve_doubled = false;
                    self.curve_filled = false;
                    self.curve_shaded = false;
                    self.curve_faded = false;
                    self.curve_arrow_start = false;
                    self.curve_arrow_end = false;
                    self.curve_half_start = 0;
                    self.curve_half_end = 0;
                }
                "dashed" => self.curve_dashed = !self.curve_dashed,
                "bold" => self.curve_bold = !self.curve_bold,
                "doubled" => self.curve_doubled = !self.curve_doubled,
                "arrow_end" => self.curve_arrow_end = !self.curve_arrow_end,
                "arrow_start" => self.curve_arrow_start = !self.curve_arrow_start,
                "closed" => {
                    self.doc_mut().apply_scene_tx("Close Curve", |s, _| {
                        if let Some(c) = s.curves.last_mut() { c.closed = !c.closed; }
                    });
                }
                "filled" => self.curve_filled = !self.curve_filled,
                "shaded" => self.curve_shaded = !self.curve_shaded,
                "faded" => self.curve_faded = !self.curve_faded,
                _ => {}
            },
            PaletteCommand::WindowAction(k) => match k {
                "cascade" => {
                    for (i, t) in self.tabs.iter_mut().enumerate() {
                        let c = ((i % 8) as f32) * 26.0;
                        t.win_pos = Vec2::new(24.0 + c, 18.0 + c);
                        t.win_size = Vec2::new(740.0, 600.0);
                        t.maximized = false;
                        t.minimized = false;
                    }
                }
                "tile_v" => {
                    let n = self.tabs.len().max(1) as f32;
                    let w = (1180.0 / n).max(300.0);
                    for (i, t) in self.tabs.iter_mut().enumerate() {
                        t.win_pos = Vec2::new(8.0 + (i as f32) * (w + 6.0), 8.0);
                        t.win_size = Vec2::new(w, 680.0);
                        t.maximized = false;
                        t.minimized = false;
                    }
                }
                "tile_h" => {
                    let n = self.tabs.len().max(1) as f32;
                    let h = (680.0 / n).max(200.0);
                    for (i, t) in self.tabs.iter_mut().enumerate() {
                        t.win_pos = Vec2::new(8.0, 8.0 + (i as f32) * (h + 6.0));
                        t.win_size = Vec2::new(1180.0, h);
                        t.maximized = false;
                        t.minimized = false;
                    }
                }
                "analysis" => self.show_analysis = !self.show_analysis,
                "properties" => self.show_properties = !self.show_properties,
                "ptable" => self.show_periodic_table = !self.show_periodic_table,
                "info" => self.show_info_window = !self.show_info_window,
                "charmap" => self.show_char_map = !self.show_char_map,
                _ => {}
            },
            PaletteCommand::None => {}
        }

        // ---- 4. Bottom Status Bar ----
        let status_live_name: Option<String> = self
            .active_mol()
            .cloned()
            .map(|m| self.compute_live_name(&m).0)
            .filter(|s| !s.is_empty());
        egui::TopBottomPanel::bottom("status_bar")
            .frame(egui::Frame::side_top_panel(&ctx.style()).fill(CD_WINDOW_BG).inner_margin(egui::Margin::symmetric(8, 2)))
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label(&self.status_msg);
                    if let Some(ref nm) = status_live_name {
                        ui.separator();
                        let short: String = if nm.chars().count() > 54 {
                            format!("{}…", nm.chars().take(53).collect::<String>())
                        } else {
                            nm.clone()
                        };
                        ui.label(egui::RichText::new(format!("Name: {short}")).strong().color(Color32::from_rgb(20, 75, 145)))
                            .on_hover_text(nm);
                    }
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        let warn_n = self.doc().check_structure_with_preferences().len();
                        ui.label(format!(
                            "Warnings: {warn_n}  |  Nicknames: {}  |  {} ({:.2} pt)  |  Docs: {}",
                            self.doc().nicknames.len(),
                            self.doc().settings.stylesheet_name,
                            self.doc().settings.bond_length_pt,
                            self.tabs.len()
                        ));
                    });
                });
            });

        // ---- 4b. Unified Right-Side Tabbed Studio Dock (MUST run before CentralPanel so canvas & floating windows never overlap it!) ----
        if self.show_right_dock && self.dock_tab == DockTab::Quantum {
            self.qc_runner.ensure_detected(
                &self.user_config.xtb_path,
                &self.user_config.orca_path,
                &self.user_config.gaussian_path,
            );
        }
        if self.show_right_dock {
            let active_mol_for_dock = self.active_mol().cloned();
            let live_nm_owned = active_mol_for_dock.as_ref().map(|m| self.compute_live_name(m).0);
            let mut dock_open = self.show_right_dock;
            let mut dock_tab = self.dock_tab;
            let mut rear_op = self.rear_opacity;
            let mut cross_gap = self.crossing_clearance_pt;
            let mut exp_opts = self.export_opts;
            let dock_act = {
                let doc = &mut self.tabs[self.active_tab].doc;
                show_right_studio_dock(
                    ctx,
                    &mut dock_open,
                    &mut dock_tab,
                    active_mol_for_dock.as_ref(),
                    live_nm_owned.as_deref(),
                    &mut self.quantum,
                    &mut self.user_config,
                    &self.qc_runner.binaries,
                    self.omni_suite.ir_quantum_pred.as_ref(),
                    &mut doc.scene.tables,
                    &mut doc.scene.tlc_plates,
                    &mut rear_op,
                    &mut cross_gap,
                    &mut exp_opts,
                )
            };
            let dock_toggled = dock_open != self.show_right_dock;
            self.show_right_dock = dock_open;
            self.dock_tab = dock_tab;
            self.rear_opacity = rear_op;
            self.crossing_clearance_pt = cross_gap;
            self.export_opts = exp_opts;
            if dock_toggled {
                self.persist_user_config();
            }
            self.handle_dock_quantum_action(ctx, dock_act);
        }

        // ---- 5. Central MDI Desk + Moveable/Resizeable Document Windows ----
        let desk_bg = match self.ui_theme {
            UiTheme::ClassicLight => CD_DESK_GRAY,
            UiTheme::DarkChrome | UiTheme::DarkSlide => Color32::from_rgb(24, 27, 33),
        };
        egui::CentralPanel::default()
            .frame(egui::Frame::central_panel(&ctx.style()).fill(desk_bg).inner_margin(egui::Margin::same(0)))
            .show(ctx, |ui| {
                let (desk_rect, response) =
                    ui.allocate_exact_size(ui.available_size(), Sense::click_and_drag());
                self.last_desk_size = desk_rect.size();
                let painter = ui.painter_at(desk_rect);
                painter.rect_filled(desk_rect, 0.0, desk_bg);

                let title_h = 22.0_f32;
                let sb_w = 13.0_f32;
                let border_w = 4.0_f32;

                // Draw inactive background MDI windows (click any to bring to front!)
                let mut bg_rects: Vec<(usize, Rect)> = Vec::new();
                let mut min_slot = 0usize;
                for (ti, t) in self.tabs.iter().enumerate() {
                    if ti == self.active_tab {
                        continue;
                    }
                    let fr = if t.minimized {
                        let r = Rect::from_min_size(
                            Pos2::new(desk_rect.min.x + 8.0 + (min_slot as f32) * 190.0, desk_rect.max.y - 28.0),
                            Vec2::new(182.0, 24.0),
                        );
                        min_slot += 1;
                        r
                    } else if t.maximized {
                        desk_rect.shrink(1.0)
                    } else {
                        Rect::from_min_size(desk_rect.min + t.win_pos, t.win_size)
                    };
                    painter.rect_filled(fr.translate(Vec2::splat(2.0)), 2.0, Color32::from_black_alpha(40));
                    painter.rect_filled(fr, 2.0, Color32::from_rgb(205, 215, 226));
                    painter.rect_stroke(fr, 2.0, Stroke::new(1.0_f32, Color32::from_gray(120)), egui::StrokeKind::Inside);
                    let tr = Rect::from_min_size(fr.min + Vec2::new(1.0, 1.0), Vec2::new((fr.width() - 2.0).max(10.0), title_h));
                    painter.rect_filled(tr, 1.0, Color32::from_rgb(212, 222, 234));
                    painter.text(
                        tr.left_center() + Vec2::new(8.0, 0.0),
                        egui::Align2::LEFT_CENTER,
                        &t.doc.name,
                        FontId::proportional(11.0),
                        Color32::from_gray(70),
                    );
                    if !t.minimized {
                        let vr = Rect::from_min_max(
                            Pos2::new(fr.min.x + border_w, fr.min.y + title_h + 2.0),
                            Pos2::new((fr.max.x - border_w).max(fr.min.x + 12.0), (fr.max.y - border_w).max(fr.min.y + title_h + 12.0)),
                        );
                        painter.rect_filled(vr, 0.0, Color32::from_gray(248));
                    }
                    bg_rects.push((ti, fr));
                }

                // Active MDI document geometry (moveable, resizeable, minimizable, maximizable!)
                let is_min = self.tab().minimized;
                let is_max = self.tab().maximized;
                if !is_max && !is_min {
                    let max_x = (desk_rect.width() - 120.0).max(0.0);
                    let max_y = (desk_rect.height() - 40.0).max(0.0);
                    let t = self.tab_mut();
                    t.win_pos.x = t.win_pos.x.clamp(0.0, max_x);
                    t.win_pos.y = t.win_pos.y.clamp(0.0, max_y);
                    t.win_size.x = t.win_size.x.clamp(260.0, desk_rect.width().max(260.0));
                    t.win_size.y = t.win_size.y.clamp(180.0, desk_rect.height().max(180.0));
                }

                let frame_rect = if is_min {
                    Rect::from_min_size(
                        Pos2::new(desk_rect.min.x + 8.0 + (min_slot as f32) * 190.0, desk_rect.max.y - 28.0),
                        Vec2::new(200.0, 24.0),
                    )
                } else if is_max {
                    desk_rect.shrink(1.0)
                } else {
                    Rect::from_min_size(desk_rect.min + self.tab().win_pos, self.tab().win_size)
                };

                let title_rect = Rect::from_min_size(
                    frame_rect.min + Vec2::new(1.0, 1.0),
                    Vec2::new((frame_rect.width() - 2.0).max(20.0), title_h),
                );
                let close_r = Rect::from_min_size(Pos2::new(title_rect.max.x - 24.0, title_rect.min.y + 3.0), Vec2::new(20.0, 15.0));
                let max_r = Rect::from_min_size(Pos2::new(title_rect.max.x - 46.0, title_rect.min.y + 3.0), Vec2::new(20.0, 15.0));
                let min_r = Rect::from_min_size(Pos2::new(title_rect.max.x - 68.0, title_rect.min.y + 3.0), Vec2::new(20.0, 15.0));

                // Draw active MDI window frame + steel-blue gradient title bar (`p08-opened.png`)
                painter.rect_filled(frame_rect.translate(Vec2::splat(3.0)), 2.0, Color32::from_black_alpha(65));
                painter.rect_filled(frame_rect, 2.0, Color32::from_rgb(166, 194, 224));
                painter.rect_stroke(frame_rect, 2.0, Stroke::new(1.0_f32, Color32::from_rgb(76, 108, 144)), egui::StrokeKind::Inside);
                draw_steel_blue_gradient(&painter, title_rect);
                let ico_r = Rect::from_min_size(title_rect.min + Vec2::new(5.0, 4.0), Vec2::new(11.0, 13.0));
                painter.rect_filled(ico_r, 0.0, Color32::WHITE);
                painter.rect_stroke(ico_r, 0.0, Stroke::new(1.0_f32, Color32::from_gray(90)), egui::StrokeKind::Inside);
                painter.text(
                    title_rect.left_center() + Vec2::new(20.0, 0.0),
                    egui::Align2::LEFT_CENTER,
                    format!("{}{}", self.doc().name, if self.doc().is_dirty() { " *" } else { "" }),
                    FontId::proportional(11.5),
                    Color32::BLACK,
                );
                for (r, bg, lbl, fg) in [
                    (min_r, Color32::from_rgb(214, 228, 244), "_", Color32::BLACK),
                    (max_r, Color32::from_rgb(214, 228, 244), if is_max { "❐" } else { "□" }, Color32::BLACK),
                    (close_r, Color32::from_rgb(198, 76, 66), "×", Color32::WHITE),
                ] {
                    painter.rect_filled(r, 1.5, bg);
                    painter.rect_stroke(r, 1.5, Stroke::new(1.0_f32, Color32::from_rgb(90, 114, 140)), egui::StrokeKind::Inside);
                    painter.text(r.center(), egui::Align2::CENTER_CENTER, lbl, FontId::proportional(10.5), fg);
                }

                // If minimized, allow clicking the minimized strip to restore
                if is_min {
                    if response.clicked() {
                        if let Some(pos) = response.interact_pointer_pos() {
                            if close_r.contains(pos) {
                                self.close_active_tab();
                            } else if frame_rect.contains(pos) {
                                self.tab_mut().minimized = false;
                            } else if let Some(&(ti, _)) = bg_rects.iter().rev().find(|(_, r)| r.contains(pos)) {
                                self.active_tab = ti;
                                self.tab_mut().minimized = false;
                                // Atom indices of the previous tab must not leak into the newly activated document.
                                self.inline_atom_edit = None;
                                self.sanitize_state();
                            }
                        }
                    }
                    return;
                }

                // Viewport inside the MDI frame + optional Rulers (`View > Show Rulers` / F11) + scrollbars + bottom-right resize grip
                let ruler_sz = if self.show_rulers { 18.0_f32 } else { 0.0_f32 };
                let view_rect = Rect::from_min_max(
                    Pos2::new(frame_rect.min.x + border_w + ruler_sz, frame_rect.min.y + title_h + 2.0 + ruler_sz),
                    Pos2::new((frame_rect.max.x - border_w - sb_w).max(frame_rect.min.x + 40.0 + ruler_sz), (frame_rect.max.y - border_w - sb_w).max(frame_rect.min.y + title_h + 40.0 + ruler_sz)),
                );
                let v_sb = Rect::from_min_size(Pos2::new(view_rect.max.x, view_rect.min.y), Vec2::new(sb_w, view_rect.height()));
                let h_sb = Rect::from_min_size(Pos2::new(view_rect.min.x, view_rect.max.y), Vec2::new(view_rect.width(), sb_w));
                let grip_rect = Rect::from_min_size(Pos2::new(view_rect.max.x, view_rect.max.y), Vec2::splat(sb_w + border_w));

                let (base_pw, base_ph) = self.tab().effective_page_size_pt();
                let wp = self.doc().scene.header_footer.width_pages.max(1) as f32;
                let hp = self.doc().scene.header_footer.height_pages.max(1) as f32;
                let page_w = (base_pw as f32) * wp * self.zoom;
                let page_h = (base_ph as f32) * hp * self.zoom;
                let max_scroll_x = (page_w - view_rect.width()).max(0.0);
                let max_scroll_y = (page_h - view_rect.height()).max(0.0);
                {
                    let t = self.tab_mut();
                    t.scroll.x = t.scroll.x.clamp(0.0, max_scroll_x);
                    t.scroll.y = t.scroll.y.clamp(0.0, max_scroll_y);
                }
                let scroll_px = self.tab().scroll;
                let page_origin = view_rect.min - scroll_px;
                let page_sheet = Rect::from_min_size(page_origin, Vec2::new(page_w, page_h));

                // Draw genuine horizontal & vertical cm rulers when `View > Show Rulers` is active!
                if self.show_rulers {
                    let h_ruler = Rect::from_min_max(
                        Pos2::new(view_rect.min.x, view_rect.min.y - ruler_sz),
                        Pos2::new(view_rect.max.x, view_rect.min.y),
                    );
                    let v_ruler = Rect::from_min_max(
                        Pos2::new(view_rect.min.x - ruler_sz, view_rect.min.y),
                        Pos2::new(view_rect.min.x, view_rect.max.y),
                    );
                    let corner_r = Rect::from_min_size(Pos2::new(v_ruler.min.x, h_ruler.min.y), Vec2::splat(ruler_sz));
                    painter.rect_filled(corner_r, 0.0, Color32::from_gray(236));
                    painter.rect_stroke(corner_r, 0.0, Stroke::new(1.0_f32, Color32::from_gray(150)), egui::StrokeKind::Inside);
                    painter.text(corner_r.center(), egui::Align2::CENTER_CENTER, "cm", FontId::proportional(7.5), Color32::DARK_GRAY);

                    let hp_r = ui.painter_at(h_ruler);
                    hp_r.rect_filled(h_ruler, 0.0, Color32::from_gray(246));
                    hp_r.rect_stroke(h_ruler, 0.0, Stroke::new(1.0_f32, Color32::from_gray(150)), egui::StrokeKind::Inside);
                    let cm_px = 28.3465_f32 * self.zoom;
                    let max_cm_x = ((base_pw as f32 * wp) / 28.3465).ceil() as i32;
                    for half in 0..=(max_cm_x * 2) {
                        let x = page_origin.x + (half as f32) * 0.5 * cm_px;
                        if x < h_ruler.min.x - 10.0 || x > h_ruler.max.x + 10.0 { continue; }
                        let tick_h = if half % 2 == 0 { 7.0 } else { 4.0 };
                        hp_r.line_segment([Pos2::new(x, h_ruler.max.y - tick_h), Pos2::new(x, h_ruler.max.y)], Stroke::new(1.0_f32, Color32::from_gray(80)));
                        if half % 2 == 0 {
                            hp_r.text(Pos2::new(x + 2.0, h_ruler.min.y + 1.0), egui::Align2::LEFT_TOP, format!("{}", half / 2), FontId::proportional(8.5), Color32::BLACK);
                        }
                    }

                    let vp_r = ui.painter_at(v_ruler);
                    vp_r.rect_filled(v_ruler, 0.0, Color32::from_gray(246));
                    vp_r.rect_stroke(v_ruler, 0.0, Stroke::new(1.0_f32, Color32::from_gray(150)), egui::StrokeKind::Inside);
                    let max_cm_y = ((base_ph as f32 * hp) / 28.3465).ceil() as i32;
                    for half in 0..=(max_cm_y * 2) {
                        let y = page_origin.y + (half as f32) * 0.5 * cm_px;
                        if y < v_ruler.min.y - 10.0 || y > v_ruler.max.y + 10.0 { continue; }
                        let tick_w = if half % 2 == 0 { 7.0 } else { 4.0 };
                        vp_r.line_segment([Pos2::new(v_ruler.max.x - tick_w, y), Pos2::new(v_ruler.max.x, y)], Stroke::new(1.0_f32, Color32::from_gray(80)));
                        if half % 2 == 0 {
                            vp_r.text(Pos2::new(v_ruler.min.x + 2.0, y + 1.0), egui::Align2::LEFT_TOP, format!("{}", half / 2), FontId::proportional(8.0), Color32::BLACK);
                        }
                    }
                    if let Some((cx, cy)) = self.cursor_doc_pos {
                        let sx = page_origin.x + (cx as f32) * self.zoom;
                        let sy = page_origin.y + (cy as f32) * self.zoom;
                        hp_r.line_segment([Pos2::new(sx, h_ruler.min.y), Pos2::new(sx, h_ruler.max.y)], Stroke::new(1.0_f32, Color32::from_rgb(0, 120, 215)));
                        vp_r.line_segment([Pos2::new(v_ruler.min.x, sy), Pos2::new(v_ruler.max.x, sy)], Stroke::new(1.0_f32, Color32::from_rgb(0, 120, 215)));
                    }
                }

                // Draw viewport background + themed page sheet clipped inside viewport
                let sheet_bg = self.canvas_bg_color();
                let vp_bg = if self.ui_theme == UiTheme::ClassicLight {
                    Color32::from_gray(225)
                } else {
                    Color32::from_rgb(20, 22, 27)
                };
                painter.rect_filled(view_rect, 0.0, vp_bg);
                let vp_painter = ui.painter_at(view_rect);
                vp_painter.rect_filled(page_sheet, 0.0, sheet_bg);
                vp_painter.rect_stroke(page_sheet, 0.0, Stroke::new(1.0_f32, Color32::from_gray(145)), egui::StrokeKind::Inside);
                // Multi-page tiled boundaries (`P17` WidthPages x HeightPages dashed lines)
                for px_i in 1..(self.doc().scene.header_footer.width_pages.max(1)) {
                    let x = page_origin.x + (px_i as f32) * (base_pw as f32) * self.zoom;
                    let mut y = page_sheet.min.y;
                    while y < page_sheet.max.y {
                        vp_painter.line_segment([Pos2::new(x, y), Pos2::new(x, (y + 6.0).min(page_sheet.max.y))], Stroke::new(1.0_f32, Color32::from_gray(180)));
                        y += 12.0;
                    }
                }
                for py_i in 1..(self.doc().scene.header_footer.height_pages.max(1)) {
                    let y = page_origin.y + (py_i as f32) * (base_ph as f32) * self.zoom;
                    let mut x = page_sheet.min.x;
                    while x < page_sheet.max.x {
                        vp_painter.line_segment([Pos2::new(x, y), Pos2::new((x + 6.0).min(page_sheet.max.x), y)], Stroke::new(1.0_f32, Color32::from_gray(180)));
                        x += 12.0;
                    }
                }
                painter.rect_stroke(view_rect, 0.0, Stroke::new(1.0_f32, Color32::from_gray(135)), egui::StrokeKind::Inside);

                // Proportional vertical & horizontal scrollbars + classic resize grip hatch
                painter.rect_filled(v_sb, 0.0, Color32::from_gray(240));
                painter.rect_filled(h_sb, 0.0, Color32::from_gray(240));
                let v_ratio = (view_rect.height() / page_h.max(1.0)).clamp(0.12, 1.0);
                let v_thumb_h = (v_sb.height() * v_ratio).max(22.0);
                let v_thumb_y = if max_scroll_y > 0.5 {
                    v_sb.min.y + (scroll_px.y / max_scroll_y) * (v_sb.height() - v_thumb_h)
                } else {
                    v_sb.min.y
                };
                let v_thumb = Rect::from_min_size(Pos2::new(v_sb.min.x + 2.0, v_thumb_y), Vec2::new(sb_w - 4.0, v_thumb_h));
                painter.rect_filled(v_thumb, 2.0, Color32::from_gray(178));

                let h_ratio = (view_rect.width() / page_w.max(1.0)).clamp(0.12, 1.0);
                let h_thumb_w = (h_sb.width() * h_ratio).max(22.0);
                let h_thumb_x = if max_scroll_x > 0.5 {
                    h_sb.min.x + (scroll_px.x / max_scroll_x) * (h_sb.width() - h_thumb_w)
                } else {
                    h_sb.min.x
                };
                let h_thumb = Rect::from_min_size(Pos2::new(h_thumb_x, h_sb.min.y + 2.0), Vec2::new(h_thumb_w, sb_w - 4.0));
                painter.rect_filled(h_thumb, 2.0, Color32::from_gray(178));

                // Corner resize grip
                painter.rect_filled(grip_rect, 0.0, Color32::from_gray(232));
                for k in 0..3 {
                    let off = 4.0 + (k as f32) * 3.5;
                    painter.line_segment(
                        [Pos2::new(grip_rect.max.x - off, grip_rect.max.y - 2.0), Pos2::new(grip_rect.max.x - 2.0, grip_rect.max.y - off)],
                        Stroke::new(1.0_f32, Color32::from_gray(120)),
                    );
                }

                // Helper to classify MDI border/title/scrollbar hit zones
                let classify_mdi_zone = |pos: Pos2| -> MdiDragMode {
                    if !is_max {
                        let edge = 7.0;
                        if grip_rect.expand(4.0).contains(pos) { return MdiDragMode::ResizeSE; }
                        let near_r = (pos.x - frame_rect.max.x).abs() <= edge && pos.y >= frame_rect.min.y && pos.y <= frame_rect.max.y;
                        let near_b = (pos.y - frame_rect.max.y).abs() <= edge && pos.x >= frame_rect.min.x && pos.x <= frame_rect.max.x;
                        let near_l = (pos.x - frame_rect.min.x).abs() <= edge && pos.y >= frame_rect.min.y && pos.y <= frame_rect.max.y;
                        let near_t = (pos.y - frame_rect.min.y).abs() <= 4.0 && pos.x >= frame_rect.min.x && pos.x <= frame_rect.max.x;
                        if near_r && near_b { return MdiDragMode::ResizeSE; }
                        if near_r { return MdiDragMode::ResizeE; }
                        if near_b { return MdiDragMode::ResizeS; }
                        if near_l { return MdiDragMode::ResizeW; }
                        if near_t { return MdiDragMode::ResizeN; }
                    }
                    if title_rect.contains(pos) && !close_r.contains(pos) && !max_r.contains(pos) && !min_r.contains(pos) {
                        return MdiDragMode::MoveTitle;
                    }
                    if v_sb.contains(pos) { return MdiDragMode::ScrollV; }
                    if h_sb.contains(pos) { return MdiDragMode::ScrollH; }
                    MdiDragMode::None
                };

                // Hover detection + cursor feedback + wheel scroll / Ctrl+wheel zoom
                if let Some(pos) = response.hover_pos() {
                    let zone = if self.mdi_drag != MdiDragMode::None { self.mdi_drag } else { classify_mdi_zone(pos) };
                    match zone {
                        MdiDragMode::ResizeSE => ctx.set_cursor_icon(egui::CursorIcon::ResizeNwSe),
                        MdiDragMode::ResizeE | MdiDragMode::ResizeW => ctx.set_cursor_icon(egui::CursorIcon::ResizeHorizontal),
                        MdiDragMode::ResizeS | MdiDragMode::ResizeN => ctx.set_cursor_icon(egui::CursorIcon::ResizeVertical),
                        MdiDragMode::MoveTitle => ctx.set_cursor_icon(if self.mdi_drag == MdiDragMode::MoveTitle { egui::CursorIcon::Grabbing } else { egui::CursorIcon::Grab }),
                        _ => {}
                    }
                    if view_rect.contains(pos) && self.mdi_drag == MdiDragMode::None {
                        let (dx, dy) = self.screen_to_doc(page_origin, pos);
                        let atol = self.atom_hit_tol();
                        let btol = self.bond_hit_tol();
                        self.cursor_doc_pos = Some((dx, dy));
                        self.hovered_atom = self.hit_atom(dx, dy, atol);
                        self.hovered_bond = if self.hovered_atom.is_none() {
                            self.hit_bond(dx, dy, btol)
                        } else {
                            None
                        };
                        let (scroll_delta, ctrl_down) = ctx.input(|i| (i.raw_scroll_delta, i.modifiers.ctrl || i.modifiers.command));
                        if scroll_delta.y.abs() > 0.1 || scroll_delta.x.abs() > 0.1 {
                            if ctrl_down {
                                // Smooth Zoom-to-Cursor: keep document point (dx, dy) anchored under pointer `pos`!
                                let old_z = self.zoom;
                                let new_z = (old_z * (1.0 + scroll_delta.y * 0.0025)).clamp(0.4, 6.0);
                                self.zoom = new_z;
                                let new_scroll_x = (dx as f32) * new_z - (pos.x - view_rect.min.x);
                                let new_scroll_y = (dy as f32) * new_z - (pos.y - view_rect.min.y);
                                let t = self.tab_mut();
                                t.scroll.x = new_scroll_x.max(0.0);
                                t.scroll.y = new_scroll_y.max(0.0);
                            } else {
                                let t = self.tab_mut();
                                t.scroll.y = (t.scroll.y - scroll_delta.y).clamp(0.0, max_scroll_y);
                                t.scroll.x = (t.scroll.x - scroll_delta.x).clamp(0.0, max_scroll_x);
                            }
                        }
                    } else {
                        self.hovered_atom = None;
                        self.hovered_bond = None;
                    }
                }
                if response.dragged_by(egui::PointerButton::Middle) {
                    let delta = response.drag_delta();
                    let t = self.tab_mut();
                    t.scroll.x = (t.scroll.x - delta.x).clamp(0.0, max_scroll_x);
                    t.scroll.y = (t.scroll.y - delta.y).clamp(0.0, max_scroll_y);
                }

                // Right-click context menu on atom / bond / canvas!
                // Capture the target at right-click time: hover keeps updating while the menu is open.
                let ctx_target_id = egui::Id::new("cc_canvas_context_target");
                if response.secondary_clicked() {
                    let tgt = (self.hovered_atom, self.hovered_bond);
                    ctx.data_mut(|d| d.insert_temp(ctx_target_id, tgt));
                }
                let (ctx_atom, ctx_bond): (Option<(usize, usize)>, Option<(usize, usize)>) = ctx
                    .data_mut(|d| d.get_temp(ctx_target_id))
                    .unwrap_or((self.hovered_atom, self.hovered_bond));
                response.context_menu(|ui| {
                    if let Some((mi, ai)) = ctx_atom {
                        ui.label(egui::RichText::new("Atom").strong());
                        if ui.button("Edit Label / Nickname... (Enter)").clicked() {
                            self.open_inline_atom_edit(mi, ai);
                            ui.close_menu();
                        }
                        if ui.button("Atom Properties...").clicked() {
                            self.dialog = ActiveDialog::AtomProperties { mol_idx: mi, atom_idx: ai };
                            ui.close_menu();
                        }
                        if ui.button("Charge +1").clicked() {
                            self.doc_mut().apply_scene_tx("Charge +", |s, _| {
                                if let Some(a) = s.molecules.get_mut(mi).and_then(|m| m.atoms.get_mut(ai)) { a.charge += 1; }
                            });
                            ui.close_menu();
                        }
                        if ui.button("Charge -1").clicked() {
                            self.doc_mut().apply_scene_tx("Charge -", |s, _| {
                                if let Some(a) = s.molecules.get_mut(mi).and_then(|m| m.atoms.get_mut(ai)) { a.charge -= 1; }
                            });
                            ui.close_menu();
                        }
                        if ui.button("Toggle Radical (•)").clicked() {
                            self.doc_mut().apply_scene_tx("Radical", |s, _| {
                                if let Some(a) = s.molecules.get_mut(mi).and_then(|m| m.atoms.get_mut(ai)) {
                                    a.radical = if a.radical == 2 { 0 } else { 2 };
                                }
                            });
                            ui.close_menu();
                        }
                        if let Some(lbl) = self
                            .doc()
                            .scene
                            .molecules
                            .get(mi)
                            .and_then(|m| m.atoms.get(ai))
                            .and_then(|a| a.label.clone())
                        {
                            if ui.button(format!("Expand Label \"{lbl}\"")).clicked() {
                                self.doc_mut().expand_label_tx(mi, ai, &lbl);
                                ui.close_menu();
                            }
                        }
                        ui.separator();
                    } else if let Some((mi, bi)) = ctx_bond {
                        ui.label(egui::RichText::new("Bond").strong());
                        for (lbl, ord, disp) in [
                            ("Single", 1, BondDisplay::Solid),
                            ("Double", 2, BondDisplay::Solid),
                            ("Triple", 3, BondDisplay::Solid),
                            ("Wedged", 1, BondDisplay::WedgeBegin),
                            ("Hashed Wedged", 1, BondDisplay::WedgedHashBegin),
                            ("Wavy (Either)", 1, BondDisplay::Wavy),
                            ("Dative", 1, BondDisplay::Dative),
                        ] {
                            if ui.button(lbl).clicked() {
                                self.doc_mut().apply_scene_tx("Bond Style", |s, _| {
                                    if let Some(b) = s.molecules.get_mut(mi).and_then(|m| m.bonds.get_mut(bi)) {
                                        b.order = ord;
                                        b.display = disp;
                                    }
                                });
                                ui.close_menu();
                            }
                        }
                        if ui.button("Bond Properties...").clicked() {
                            self.dialog = ActiveDialog::BondProperties { mol_idx: mi, bond_idx: bi };
                            ui.close_menu();
                        }
                        ui.separator();
                    }
                    if !self.doc().scene.molecules.is_empty() {
                        if ui.button("Cut\tCtrl+X").clicked() { self.do_cut_action(); ui.close_menu(); }
                        if ui.button("Copy\tCtrl+C").clicked() { self.do_copy_action(); ui.close_menu(); }
                    }
                    if !self.clipboard_mols.is_empty() {
                        if ui.button("Paste\tCtrl+V").clicked() {
                            let clip = self.clipboard_mols.clone();
                            for mut m in clip {
                                for a in &mut m.atoms { a.x += 20.0; a.y += 20.0; }
                                self.doc_mut().stamp_molecule_tx(m, "Paste");
                            }
                            ui.close_menu();
                        }
                    }
                    ui.separator();
                    if ui.button("Flip Horizontal").clicked() { self.do_flip_or_rot180("Flip Horizontal", FlipAxis::Horizontal, false); ui.close_menu(); }
                    if ui.button("Flip Vertical").clicked() { self.do_flip_or_rot180("Flip Vertical", FlipAxis::Vertical, false); ui.close_menu(); }
                    if ui.button("Clean Up Structure").clicked() {
                        self.do_clean_up_structure();
                        ui.close_menu();
                    }
                    if ui.button("Name Structure\tAlt+Ctrl+N").clicked() {
                        self.do_structure_to_name_hybrid(ctx, false);
                        ui.close_menu();
                    }
                    if ui.button("Resolve Name via PubChem").clicked() {
                        self.do_structure_to_name_hybrid(ctx, true);
                        ui.close_menu();
                    }
                    if ui.button("Check Structure").clicked() {
                        self.check_warn_idx = 0;
                        self.dialog = ActiveDialog::CheckStructureDialog;
                        ui.close_menu();
                    }
                });

                // Compute interactive 8 selection resize handles + top rotation knob when selection is active
                let sel_handles_screen: Option<(Rect, [Pos2; 8], Pos2)> = if !self.tab().selected_mols.is_empty() {
                    self.selection_bbox().map(|(x1, y1, x2, y2)| {
                        let s1 = self.doc_to_screen(page_origin, x1 - 5.0, y1 - 5.0);
                        let s2 = self.doc_to_screen(page_origin, x2 + 5.0, y2 + 5.0);
                        let sr = Rect::from_two_pos(s1, s2);
                        let hs = [
                            sr.left_top(), sr.center_top(), sr.right_top(),
                            sr.left_center(), sr.right_center(),
                            sr.left_bottom(), sr.center_bottom(), sr.right_bottom(),
                        ];
                        let rot_knob = sr.center_top() - Vec2::new(0.0, 14.0);
                        (sr, hs, rot_knob)
                    })
                } else {
                    None
                };

                let shift_held = ctx.input(|i| i.modifiers.shift);
                let alt_held = ctx.input(|i| i.modifiers.alt);
                let atol = self.atom_hit_tol();
                let mut double_click_consumed = false;
                if response.double_clicked() {
                    if let Some(pos) = response.interact_pointer_pos() {
                        if title_rect.contains(pos) && !close_r.contains(pos) && !max_r.contains(pos) && !min_r.contains(pos) {
                            self.tab_mut().maximized = !self.tab().maximized;
                            double_click_consumed = true;
                        } else if view_rect.contains(pos)
                            && matches!(self.tool, Tool::Marquee | Tool::Lasso | Tool::TextLabel)
                        {
                            let (dx, dy) = self.screen_to_doc(page_origin, pos);
                            if let Some((mi, ai)) = self.hit_atom(dx, dy, atol) {
                                self.open_inline_atom_edit(mi, ai);
                                double_click_consumed = true;
                            }
                        }
                    }
                }
                if response.drag_started_by(egui::PointerButton::Primary) {
                    // Drags start after egui's movement threshold; hit-test where the button went
                    // down, or a drag begun on an atom/handle starts beside it.
                    let press = ctx.input(|i| i.pointer.press_origin());
                    if let Some(pos) = press.or_else(|| response.interact_pointer_pos()) {
                        let zone = classify_mdi_zone(pos);
                        if zone != MdiDragMode::None {
                            self.mdi_drag = zone;
                        } else if let Some((_, hs, rot_knob)) = sel_handles_screen {
                            let mut hit_h = None;
                            if pos.distance(rot_knob) <= 7.0 {
                                hit_h = Some(8usize);
                            } else {
                                for (hi, &hp) in hs.iter().enumerate() {
                                    if pos.distance(hp) <= 6.5 {
                                        hit_h = Some(hi);
                                        break;
                                    }
                                }
                            }
                            if let Some(hi) = hit_h {
                                let start_pt = self.screen_to_doc(page_origin, pos);
                                self.sel_handle_drag = Some(hi);
                                self.sel_drag_prev = Some(start_pt);
                                self.sel_drag_start_pt = Some(start_pt);
                                self.sel_drag_start_bbox = self.selection_bbox();
                                self.sel_drag_start_mols = Some(self.doc().scene.molecules.clone());
                            } else if view_rect.contains(pos) {
                                let (dx, dy) = self.screen_to_doc(page_origin, pos);
                                self.drag_start_canvas = Some((dx, dy));
                                self.drag_cur_canvas = Some((dx, dy));
                                self.drag_start_atom = self.hit_atom(dx, dy, atol);
                                self.dragging_single_atom = None;
                                self.dragging_caption = None;
                                self.eraser_pts.clear();
                                if self.tool == Tool::Eraser {
                                    self.eraser_pts.push((dx, dy));
                                }
                                if matches!(self.tool, Tool::Marquee | Tool::Lasso) {
                                    if alt_held && self.drag_start_atom.is_some() {
                                        self.dragging_single_atom = self.drag_start_atom;
                                        self.dragging_selected_move = false;
                                    } else if self.hit_selected_molecule(dx, dy) {
                                        self.dragging_selected_move = true;
                                    } else if let Some((mi, ai)) = self.drag_start_atom {
                                        // Dragging an unselected atom vertex directly reshapes that vertex!
                                        self.dragging_single_atom = Some((mi, ai));
                                        self.dragging_selected_move = false;
                                    } else if let Some(mi) = self.hit_molecule_bbox(dx, dy) {
                                        if !shift_held {
                                            self.tab_mut().selected_mols = vec![mi];
                                        } else if !self.tab().selected_mols.contains(&mi) {
                                            self.tab_mut().selected_mols.push(mi);
                                        }
                                        self.expand_selection_with_groups();
                                        self.dragging_selected_move = true;
                                    } else if let Some(ci) = self.hit_caption(dx, dy) {
                                        self.dragging_selected_move = false;
                                        self.dragging_caption = Some(ci);
                                    } else {
                                        self.dragging_selected_move = false;
                                    }
                                } else {
                                    self.dragging_selected_move = false;
                                }
                                self.drag_ctrl_copy = ctx.input(|i| i.modifiers.ctrl || i.modifiers.command);
                                self.lasso_pts.clear();
                                if self.tool == Tool::Lasso && !self.dragging_selected_move && self.dragging_caption.is_none() && self.dragging_single_atom.is_none() {
                                    self.lasso_pts.push((dx, dy));
                                }
                            }
                        } else if view_rect.contains(pos) {
                            let (dx, dy) = self.screen_to_doc(page_origin, pos);
                            self.drag_start_canvas = Some((dx, dy));
                            self.drag_cur_canvas = Some((dx, dy));
                            self.drag_start_atom = self.hit_atom(dx, dy, atol);
                            self.dragging_single_atom = None;
                            self.dragging_caption = None;
                            self.eraser_pts.clear();
                            if self.tool == Tool::Eraser {
                                self.eraser_pts.push((dx, dy));
                            }
                            if matches!(self.tool, Tool::Marquee | Tool::Lasso) {
                                if let Some((mi, ai)) = self.drag_start_atom {
                                    self.dragging_single_atom = Some((mi, ai));
                                    self.dragging_selected_move = false;
                                } else if let Some(mi) = self.hit_molecule_bbox(dx, dy) {
                                    if !shift_held {
                                        self.tab_mut().selected_mols = vec![mi];
                                    } else if !self.tab().selected_mols.contains(&mi) {
                                        self.tab_mut().selected_mols.push(mi);
                                    }
                                    self.expand_selection_with_groups();
                                    self.dragging_selected_move = true;
                                } else if let Some(ci) = self.hit_caption(dx, dy) {
                                    self.dragging_selected_move = false;
                                    self.dragging_caption = Some(ci);
                                } else {
                                    self.dragging_selected_move = false;
                                }
                            } else {
                                self.dragging_selected_move = false;
                            }
                            self.drag_ctrl_copy = ctx.input(|i| i.modifiers.ctrl || i.modifiers.command);
                            self.lasso_pts.clear();
                            if self.tool == Tool::Lasso && !self.dragging_selected_move && self.dragging_caption.is_none() && self.dragging_single_atom.is_none() {
                                self.lasso_pts.push((dx, dy));
                            }
                        }
                    }
                }
                if response.dragged_by(egui::PointerButton::Primary) {
                    let delta = response.drag_delta();
                    if let Some(hi) = self.sel_handle_drag {
                        if let (
                            Some(pos),
                            Some((sx0, sy0)),
                            Some((bx1, by1, bx2, by2)),
                            Some(base_mols),
                        ) = (
                            response.interact_pointer_pos(),
                            self.sel_drag_start_pt,
                            self.sel_drag_start_bbox,
                            self.sel_drag_start_mols.clone(),
                        ) {
                            let (cx, cy) = self.screen_to_doc(page_origin, pos);
                            let mid_x = 0.5 * (bx1 + bx2);
                            let mid_y = 0.5 * (by1 + by2);
                            let sel = self.tab().selected_mols.clone();
                            let mut preview = base_mols;
                            if hi == 8 {
                                let a0 = (sy0 - mid_y).atan2(sx0 - mid_x);
                                let a1 = (cy - mid_y).atan2(cx - mid_x);
                                rotate_selected_molecules(&mut preview, &sel, mid_x, mid_y, a1 - a0);
                            } else {
                                let w = (bx2 - bx1).max(6.0);
                                let h = (by2 - by1).max(6.0);
                                let (ax, ay, sx, sy) = match hi {
                                    0 => (bx2, by2, ((bx2 - cx) / w).clamp(0.2, 5.0), ((by2 - cy) / h).clamp(0.2, 5.0)),
                                    1 => (mid_x, by2, 1.0, ((by2 - cy) / h).clamp(0.2, 5.0)),
                                    2 => (bx1, by2, ((cx - bx1) / w).clamp(0.2, 5.0), ((by2 - cy) / h).clamp(0.2, 5.0)),
                                    3 => (bx2, mid_y, ((bx2 - cx) / w).clamp(0.2, 5.0), 1.0),
                                    4 => (bx1, mid_y, ((cx - bx1) / w).clamp(0.2, 5.0), 1.0),
                                    5 => (bx2, by1, ((bx2 - cx) / w).clamp(0.2, 5.0), ((cy - by1) / h).clamp(0.2, 5.0)),
                                    6 => (mid_x, by1, 1.0, ((cy - by1) / h).clamp(0.2, 5.0)),
                                    _ => (bx1, by1, ((cx - bx1) / w).clamp(0.2, 5.0), ((cy - by1) / h).clamp(0.2, 5.0)),
                                };
                                scale_selected_molecules(&mut preview, &sel, ax, ay, sx, sy);
                            }
                            self.doc_mut().scene.molecules = preview;
                            self.sel_drag_prev = Some((cx, cy));
                        }
                    } else {
                        match self.mdi_drag {
                            MdiDragMode::MoveTitle => {
                                if self.tab().maximized && delta.length() > 3.0 {
                                    self.tab_mut().maximized = false;
                                }
                                self.tab_mut().win_pos += delta;
                            }
                            MdiDragMode::ResizeE => {
                                self.tab_mut().win_size.x = (self.tab().win_size.x + delta.x).max(260.0);
                            }
                            MdiDragMode::ResizeS => {
                                self.tab_mut().win_size.y = (self.tab().win_size.y + delta.y).max(180.0);
                            }
                            MdiDragMode::ResizeSE => {
                                self.tab_mut().win_size.x = (self.tab().win_size.x + delta.x).max(260.0);
                                self.tab_mut().win_size.y = (self.tab().win_size.y + delta.y).max(180.0);
                            }
                            MdiDragMode::ResizeW => {
                                let new_w = (self.tab().win_size.x - delta.x).max(260.0);
                                let dw = self.tab().win_size.x - new_w;
                                self.tab_mut().win_pos.x += dw;
                                self.tab_mut().win_size.x = new_w;
                            }
                            MdiDragMode::ResizeN => {
                                let new_h = (self.tab().win_size.y - delta.y).max(180.0);
                                let dh = self.tab().win_size.y - new_h;
                                self.tab_mut().win_pos.y += dh;
                                self.tab_mut().win_size.y = new_h;
                            }
                            MdiDragMode::ScrollV => {
                                let track = (v_sb.height() - v_thumb_h).max(1.0);
                                let t = self.tab_mut();
                                t.scroll.y = (t.scroll.y + delta.y * (max_scroll_y / track)).clamp(0.0, max_scroll_y);
                            }
                            MdiDragMode::ScrollH => {
                                let track = (h_sb.width() - h_thumb_w).max(1.0);
                                let t = self.tab_mut();
                                t.scroll.x = (t.scroll.x + delta.x * (max_scroll_x / track)).clamp(0.0, max_scroll_x);
                            }
                            MdiDragMode::None => {
                                if let Some(pos) = response.interact_pointer_pos() {
                                    let pt = self.screen_to_doc(page_origin, pos);
                                    self.drag_cur_canvas = Some(pt);
                                    if self.tool == Tool::Lasso && !self.dragging_selected_move && self.dragging_single_atom.is_none() {
                                        if self.lasso_pts.last().map(|&(lx, ly)| (lx - pt.0).hypot(ly - pt.1) > 2.0).unwrap_or(true) {
                                            self.lasso_pts.push(pt);
                                        }
                                    } else if self.tool == Tool::Eraser {
                                        if self.eraser_pts.last().map(|&(lx, ly)| (lx - pt.0).hypot(ly - pt.1) > 2.0).unwrap_or(true) {
                                            self.eraser_pts.push(pt);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if response.drag_stopped_by(egui::PointerButton::Primary) {
                    if let Some(hi) = self.sel_handle_drag.take() {
                        let final_mols = self.doc().scene.molecules.clone();
                        if let Some(base_mols) = self.sel_drag_start_mols.take() {
                            if final_mols != base_mols {
                                self.doc_mut().scene.molecules = base_mols;
                                let label = if hi == 8 { "Rotate Handle" } else { "Scale Handle" };
                                self.doc_mut().apply_scene_tx(label, |s, _| {
                                    s.molecules = final_mols;
                                });
                            }
                        }
                        self.sel_drag_prev = None;
                        self.sel_drag_start_pt = None;
                        self.sel_drag_start_bbox = None;
                    } else if self.mdi_drag != MdiDragMode::None {
                        self.mdi_drag = MdiDragMode::None;
                    } else if let (Some((sx, sy)), Some((ex, ey))) = (self.drag_start_canvas, self.drag_cur_canvas) {
                        self.handle_canvas_drag_release(sx, sy, ex, ey, shift_held);
                    }
                    self.drag_start_canvas = None;
                    self.drag_cur_canvas = None;
                    self.drag_start_atom = None;
                    self.dragging_single_atom = None;
                    self.eraser_pts.clear();
                    self.dragging_selected_move = false;
                    self.dragging_caption = None;
                    self.drag_ctrl_copy = false;
                    self.lasso_pts.clear();
                } else if response.clicked() && !double_click_consumed {
                    if let Some(pos) = response.interact_pointer_pos() {
                        if close_r.contains(pos) {
                            self.close_active_tab();
                        } else if max_r.contains(pos) {
                            self.tab_mut().maximized = !self.tab().maximized;
                        } else if min_r.contains(pos) {
                            self.tab_mut().minimized = true;
                        } else if view_rect.contains(pos) {
                            let (dx, dy) = self.screen_to_doc(page_origin, pos);
                            self.handle_canvas_click(dx, dy, shift_held);
                        } else if !frame_rect.contains(pos) {
                            if let Some(&(ti, _)) = bg_rects.iter().rev().find(|(_, r)| r.contains(pos)) {
                                self.active_tab = ti;
                                self.tab_mut().minimized = false;
                                self.inline_atom_edit = None;
                                self.sanitize_state();
                            }
                        }
                    }
                }

                // Clip drawing strictly inside the MDI viewport (`view_rect`)
                let pp = ui.painter_at(view_rect);
                let lw = (self.doc().settings.line_width_pt as f32 * self.zoom).max(1.2);
                let warns = if self.show_warnings {
                    self.doc().check_structure_with_preferences()
                } else {
                    Vec::new()
                };

                // Header / Footer
                let (eff_pw, eff_ph) = self.tab().effective_page_size_pt();
                if let Some(hdr) = &self.doc().scene.header_footer.header {
                    pp.text(self.doc_to_screen(page_origin, eff_pw * 0.5, 22.0), egui::Align2::CENTER_CENTER, hdr, FontId::proportional(10.0 * self.zoom * 0.6), Color32::DARK_GRAY);
                }
                if let Some(ftr) = &self.doc().scene.header_footer.footer {
                    pp.text(self.doc_to_screen(page_origin, eff_pw * 0.5, (eff_ph - 19.75).max(32.0)), egui::Align2::CENTER_CENTER, ftr, FontId::proportional(10.0 * self.zoom * 0.6), Color32::DARK_GRAY);
                }

                // Render molecules (with Ring Fills, Highlights, Crossing-Bond Clearance, 3D Depth Cueing, Ring-Interior Double Bonds, IUPAC Labels, Enhanced Stereo & NMR glow!)
                let default_ink = self.effective_ink_color();
                let canvas_bg = self.canvas_bg_color();
                let scale_alpha = |c: Color32, op: f32| -> Color32 {
                    if op >= 0.995 {
                        c
                    } else {
                        let a = ((c.a() as f32) * op.clamp(0.08, 1.0)).round() as u8;
                        Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), a)
                    }
                };
                for (mi, mol) in self.doc().scene.molecules.iter().enumerate() {
                    let is_sel = self.tab().selected_mols.contains(&mi);
                    // Professional: Molecule Highlight Color halo
                    if let Some(&hcol) = self.tab().mol_highlights.get(&mi) {
                        if !mol.atoms.is_empty() {
                            let (x1, y1, x2, y2) = mol_bbox(mol);
                            let hr = Rect::from_two_pos(
                                self.doc_to_screen(page_origin, x1 - 7.0, y1 - 7.0),
                                self.doc_to_screen(page_origin, x2 + 7.0, y2 + 7.0),
                            );
                            pp.rect_filled(hr, 6.0, hcol);
                        }
                    }
                    // Professional: Ring Fill Color polygons underneath bonds
                    if let Some(&rfcol) = self.tab().ring_fills.get(&mi) {
                        for poly in extract_ring_polygons(mol) {
                            let spts: Vec<Pos2> = poly.into_iter().map(|(x, y)| self.doc_to_screen(page_origin, x, y)).collect();
                            if spts.len() >= 3 {
                                pp.add(egui::Shape::convex_polygon(spts, rfcol, Stroke::NONE));
                            }
                        }
                    }
                    let val = CcVal::compute(mol);
                    let stereo = if self.show_stereo { Some(perceive_stereochemistry(mol)) } else { None };
                    let ring_normals = ring_bond_inward_normals(mol);
                    let crossings = find_bond_crossings(mol);
                    let atom_opacities = compute_3d_atom_opacities(mol, self.rear_opacity);

                    for (bi, b) in mol.bonds.iter().enumerate() {
                        if !b.visible || b.begin_atom >= mol.atoms.len() || b.end_atom >= mol.atoms.len() { continue; }
                        let a1 = &mol.atoms[b.begin_atom];
                        let a2 = &mol.atoms[b.end_atom];
                        let p1 = self.doc_to_screen(page_origin, a1.x, a1.y);
                        let p2 = self.doc_to_screen(page_origin, a2.x, a2.y);
                        let d = p2 - p1;
                        let len = d.length().max(1e-3);
                        let u = d / len;
                        let n = Vec2::new(-u.y, u.x);

                        // ChiralCraft blue bond hover glow + Retrosynthesis strategic bond highlight!
                        if self.hovered_bond == Some((mi, bi)) {
                            pp.line_segment([p1, p2], Stroke::new(5.5 * self.zoom, Color32::from_rgba_unmultiplied(120, 190, 255, 130)));
                        }
                        if self.pro_suite.show_retro
                            && self.active_mol_index() == Some(mi)
                            && self.pro_suite.highlighted_retro_bond == Some(bi)
                        {
                            pp.line_segment([p1, p2], Stroke::new(6.5 * self.zoom, Color32::from_rgba_unmultiplied(220, 45, 165, 145)));
                        }

                        let bond_op = 0.5 * (
                            atom_opacities.get(b.begin_atom).copied().unwrap_or(1.0)
                                + atom_opacities.get(b.end_atom).copied().unwrap_or(1.0)
                        );
                        let base_col = if is_sel {
                            Color32::from_rgb(0, 110, 235)
                        } else if self.color_by_element {
                            if a1.z != 6 { cpk_element_color(a1.z) } else if a2.z != 6 { cpk_element_color(a2.z) } else { default_ink }
                        } else {
                            default_ink
                        };
                        let col = scale_alpha(base_col, bond_op);
                        let b_cross = crossings.get(&bi).map(|v| v.as_slice()).unwrap_or(&[]);
                        let seg_pairs = split_segment_for_crossings(
                            a1.x,
                            a1.y,
                            a2.x,
                            a2.y,
                            b_cross,
                            self.crossing_clearance_pt,
                        );
                        let bw = (self.doc().settings.bold_width_pt as f32 * self.zoom).max(lw * 2.0);
                        let wedge_hw = 0.75 * bw; // IUPAC GR-1.3 / ChiralCraft F02: wedge base = 1.5 * BoldWidth
                        let hash_step_px = (self.doc().settings.hash_spacing_pt as f32 * self.zoom).max(2.0);
                        let bond_sep = (self.doc().settings.bond_length_pt * self.doc().settings.bond_spacing_pct / 100.0) as f32 * self.zoom;
                        match b.display {
                            BondDisplay::WedgeBegin | BondDisplay::WedgeEnd => {
                                let (s, e) = if b.display == BondDisplay::WedgeBegin { (p1, p2) } else { (p2, p1) };
                                pp.add(egui::Shape::convex_polygon(vec![s, e + n * wedge_hw, e - n * wedge_hw], col, Stroke::NONE));
                            }
                            BondDisplay::HollowWedgeBegin | BondDisplay::HollowWedgeEnd => {
                                let (s, e) = if b.display == BondDisplay::HollowWedgeBegin { (p1, p2) } else { (p2, p1) };
                                pp.add(egui::Shape::convex_polygon(vec![s, e + n * wedge_hw, e - n * wedge_hw], canvas_bg, Stroke::new(lw, col)));
                            }
                            BondDisplay::WedgedHashBegin | BondDisplay::WedgedHashEnd => {
                                let (s, e) = if b.display == BondDisplay::WedgedHashBegin { (p1, p2) } else { (p2, p1) };
                                let raw_steps = if len.is_finite() { (len / hash_step_px).round() } else { 4.0 };
                                let steps = (raw_steps.clamp(4.0, 256.0)) as usize;
                                for k in 1..=steps {
                                    let t = k as f32 / steps as f32;
                                    let m = s + (e - s) * t;
                                    let hw = wedge_hw * t;
                                    pp.line_segment([m + n * hw, m - n * hw], Stroke::new(lw, col));
                                }
                            }
                            BondDisplay::Bold => {
                                for ((sx, sy), (ex, ey)) in &seg_pairs {
                                    let sp1 = self.doc_to_screen(page_origin, *sx, *sy);
                                    let sp2 = self.doc_to_screen(page_origin, *ex, *ey);
                                    pp.line_segment([sp1, sp2], Stroke::new(bw, col));
                                }
                            }
                            BondDisplay::Wavy => {
                                let raw_steps = if len.is_finite() { (len / (1.8 * self.zoom)).round() } else { 8.0 };
                                let steps = (raw_steps.clamp(8.0, 256.0)) as usize;
                                let amp = (bw * 0.7).max(1.2 * self.zoom);
                                let mut prev = p1;
                                for k in 1..=steps {
                                    let t = (k as f32) / (steps as f32);
                                    let wave = (t * std::f32::consts::TAU * 3.5).sin() * amp;
                                    let cur = p1 + d * t + n * wave;
                                    pp.line_segment([prev, cur], Stroke::new(lw, col));
                                    prev = cur;
                                }
                            }
                            BondDisplay::Dative => {
                                for ((sx, sy), (ex, ey)) in &seg_pairs {
                                    let sp1 = self.doc_to_screen(page_origin, *sx, *sy);
                                    let sp2 = self.doc_to_screen(page_origin, *ex, *ey);
                                    pp.line_segment([sp1, sp2], Stroke::new(lw, col));
                                }
                                let hl = (bw * 2.0).max(3.5 * self.zoom);
                                pp.add(egui::Shape::convex_polygon(vec![p2, p2 - u * hl + n * wedge_hw, p2 - u * hl - n * wedge_hw], col, Stroke::NONE));
                            }
                            BondDisplay::Dash | BondDisplay::Hash => {
                                let raw_steps = if len.is_finite() { (len / hash_step_px).round() } else { 4.0 };
                                let steps = (raw_steps.clamp(4.0, 256.0)) as usize;
                                for k in 0..steps {
                                    let t0 = (k as f32) / (steps as f32);
                                    let t1 = ((k as f32) + 0.55) / (steps as f32);
                                    if b.display == BondDisplay::Hash {
                                        let m = p1 + d * ((k as f32 + 0.5) / steps as f32);
                                        pp.line_segment([m + n * wedge_hw, m - n * wedge_hw], Stroke::new(lw, col));
                                    } else {
                                        // IUPAC GR-1.3: dashed bond consists of collinear dash segments
                                        pp.line_segment([p1 + d * t0, p1 + d * t1.min(1.0)], Stroke::new(lw, col));
                                    }
                                }
                            }
                            _ => {
                                if b.order == 2 {
                                    let sep = bond_sep;
                                    let inset = len * 0.16;
                                    match self.tab().double_bond_align.get(&(mi, bi)).copied().unwrap_or(DoubleBondAlign::Auto) {
                                        DoubleBondAlign::Left => {
                                            pp.line_segment([p1, p2], Stroke::new(lw, col));
                                            pp.line_segment([p1 + u * inset + n * sep, p2 - u * inset + n * sep], Stroke::new(lw, col));
                                        }
                                        DoubleBondAlign::Right => {
                                            pp.line_segment([p1, p2], Stroke::new(lw, col));
                                            pp.line_segment([p1 + u * inset - n * sep, p2 - u * inset - n * sep], Stroke::new(lw, col));
                                        }
                                        DoubleBondAlign::Center => {
                                            let off = sep * 0.5;
                                            pp.line_segment([p1 + n * off, p2 + n * off], Stroke::new(lw, col));
                                            pp.line_segment([p1 - n * off, p2 - n * off], Stroke::new(lw, col));
                                        }
                                        DoubleBondAlign::Auto => {
                                            if let Some(&(rnx, rny)) = ring_normals.get(&bi) {
                                                let in_n = Vec2::new(rnx as f32, rny as f32);
                                                pp.line_segment([p1, p2], Stroke::new(lw, col));
                                                pp.line_segment([p1 + u * inset + in_n * sep, p2 - u * inset + in_n * sep], Stroke::new(lw, col));
                                            } else {
                                                let off = sep * 0.5;
                                                pp.line_segment([p1 + n * off, p2 + n * off], Stroke::new(lw, col));
                                                pp.line_segment([p1 - n * off, p2 - n * off], Stroke::new(lw, col));
                                            }
                                        }
                                    }
                                } else if b.order == 3 {
                                    // IUPAC GR-1.6.2: triple bond uses stylesheet BondSpacing with shortened flanking lines
                                    let inset = len * 0.12;
                                    pp.line_segment([p1, p2], Stroke::new(lw, col));
                                    pp.line_segment([p1 + u * inset + n * bond_sep, p2 - u * inset + n * bond_sep], Stroke::new(lw, col));
                                    pp.line_segment([p1 + u * inset - n * bond_sep, p2 - u * inset - n * bond_sep], Stroke::new(lw, col));
                                } else if b.order == 4 {
                                    for off in [-1.5 * bond_sep, -0.5 * bond_sep, 0.5 * bond_sep, 1.5 * bond_sep] {
                                        pp.line_segment([p1 + n * off, p2 + n * off], Stroke::new(lw, col));
                                    }
                                } else if b.order == 6 {
                                    // Delocalized aromatic (1.5) bond: solid outer bond + inset dashed inner bond
                                    let inset = len * 0.16;
                                    let in_n = if let Some(&(rnx, rny)) = ring_normals.get(&bi) {
                                        Vec2::new(rnx as f32, rny as f32)
                                    } else {
                                        n
                                    };
                                    pp.line_segment([p1, p2], Stroke::new(lw, col));
                                    let q1 = p1 + u * inset + in_n * bond_sep;
                                    let q2 = p2 - u * inset + in_n * bond_sep;
                                    let dashes = 4usize;
                                    for di in 0..dashes {
                                        let t0 = (di as f32) / (dashes as f32);
                                        let t1 = ((di as f32) + 0.58) / (dashes as f32);
                                        pp.line_segment(
                                            [q1 + (q2 - q1) * t0, q1 + (q2 - q1) * t1.min(1.0)],
                                            Stroke::new(lw, col),
                                        );
                                    }
                                } else {
                                    for ((sx, sy), (ex, ey)) in &seg_pairs {
                                        let sp1 = self.doc_to_screen(page_origin, *sx, *sy);
                                        let sp2 = self.doc_to_screen(page_origin, *ex, *ey);
                                        pp.line_segment([sp1, sp2], Stroke::new(lw, col));
                                    }
                                }
                            }
                        }

                        if let Some(st) = &stereo {
                            if let Some(&bcode) = st.bond_stereo.get(&bi) {
                                let btag = match bcode { 1 => Some("(E)"), 2 => Some("(Z)"), _ => None };
                                if let Some(t) = btag {
                                    pp.text(p1 + d * 0.5 + n * (7.0 * self.zoom), egui::Align2::CENTER_CENTER, t, FontId::proportional(7.5 * self.zoom), Color32::from_rgb(0, 100, 0));
                                }
                            }
                        }
                        if let Some(bw) = warns.iter().find(|w| w.object_id == b.id) {
                            let mid = p1 + d * 0.5;
                            let wr = Rect::from_center_size(mid, Vec2::splat(11.0 * self.zoom));
                            pp.rect_stroke(wr, 0.0, Stroke::new(1.5_f32, Color32::RED), egui::StrokeKind::Outside);
                            if self.hovered_bond == Some((mi, bi)) {
                                let tip_pos = mid + Vec2::new(10.0, -16.0);
                                let galley = pp.layout_no_wrap(bw.message.clone(), FontId::proportional(11.0), Color32::BLACK);
                                let bg_r = Rect::from_min_size(tip_pos, galley.size() + Vec2::new(10.0, 6.0));
                                pp.rect_filled(bg_r, 2.0, Color32::from_rgb(255, 255, 225));
                                pp.rect_stroke(bg_r, 2.0, Stroke::new(1.0_f32, Color32::from_rgb(180, 50, 50)), egui::StrokeKind::Inside);
                                pp.galley(tip_pos + Vec2::new(5.0, 3.0), galley, Color32::BLACK);
                            }
                        }
                    }

                    for (ai, a) in mol.atoms.iter().enumerate() {
                        if !a.visible { continue; }
                        let pos = self.doc_to_screen(page_origin, a.x, a.y);

                        // Signature ChiralCraft light-blue atom hover circle!
                        if self.hovered_atom == Some((mi, ai)) {
                            pp.circle_filled(pos, 5.2 * self.zoom, Color32::from_rgba_unmultiplied(120, 190, 255, 140));
                            pp.circle_stroke(pos, 5.2 * self.zoom, Stroke::new(1.0_f32, Color32::from_rgb(40, 120, 210)));
                        }

                        let nbr_idxs: Vec<usize> = mol
                            .bonds
                            .iter()
                            .filter_map(|b| {
                                if b.begin_atom == ai { Some(b.end_atom) }
                                else if b.end_atom == ai { Some(b.begin_atom) }
                                else { None }
                            })
                            .collect();
                        let deg = nbr_idxs.len();
                        let flip_left = deg == 1
                            && nbr_idxs[0] < mol.atoms.len()
                            && mol.atoms[nbr_idxs[0]].x > a.x + 1.0;
                        let to_sup = |num: i32| -> String {
                            num.to_string()
                                .chars()
                                .map(|ch| match ch {
                                    '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
                                    '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
                                    _ => ch,
                                })
                                .collect()
                        };
                        let (text_opt, anchor_mode) = if let Some(lbl) = &a.label {
                            (
                                Some(if self.text_formula { format_chemical_formula(lbl) } else { lbl.clone() }),
                                egui::Align2::CENTER_CENTER,
                            )
                        } else if a.z != 6 || a.charge != 0 || a.isotope != 0 || a.radical != 0 || deg == 0 {
                            let sym = get_element(a.z).map(|e| e.symbol).unwrap_or("?");
                            let h = val.implicit_h_counts.get(ai).copied().unwrap_or(0);
                            let mut s = String::new();
                            if flip_left && h > 0 && a.isotope == 0 && a.charge == 0 && a.radical == 0 {
                                // IUPAC GR-2.1.4: left-facing terminal heteroatom writes HO-, H2N-, HS-, H3C-
                                s.push('H');
                                if h > 1 { s.push_str(&h.to_string()); }
                                s.push_str(sym);
                                let disp = if self.text_formula { format_chemical_formula(&s) } else { s };
                                (Some(disp), egui::Align2::RIGHT_CENTER)
                            } else {
                                if a.isotope > 0 { s.push_str(&to_sup(a.isotope)); }
                                s.push_str(sym);
                                if h == 1 { s.push('H'); } else if h > 1 { s.push_str(&format!("H{h}")); }
                                let mut disp = if self.text_formula { format_chemical_formula(&s) } else { s };
                                // IUPAC GR-2.1.3: superscript formal charge & radical dot
                                match a.charge {
                                    1 => disp.push('⁺'),
                                    -1 => disp.push('⁻'),
                                    c if c > 1 => { disp.push_str(&to_sup(c)); disp.push('⁺'); }
                                    c if c < -1 => { disp.push_str(&to_sup(-c)); disp.push('⁻'); }
                                    _ => {}
                                }
                                if a.radical == 2 { disp.push('•'); } else if a.radical == 1 { disp.push(':'); }
                                let anchor = if deg == 1 && h > 0 && !flip_left {
                                    egui::Align2::LEFT_CENTER
                                } else {
                                    egui::Align2::CENTER_CENTER
                                };
                                (Some(disp), anchor)
                            }
                        } else {
                            (None, egui::Align2::CENTER_CENTER)
                        };
                        if let Some(lbl) = text_opt {
                            let fs = (self.doc().settings.label_size_pt as f32) * self.zoom * 0.92;
                            let margin_px = (self.doc().settings.margin_width_pt as f32) * self.zoom;
                            let half_char = fs * 0.34;
                            let text_pos = match anchor_mode {
                                egui::Align2::RIGHT_CENTER => pos + Vec2::new(half_char, 0.0),
                                egui::Align2::LEFT_CENTER => pos - Vec2::new(half_char, 0.0),
                                _ => pos,
                            };
                            let kw = (lbl.chars().count() as f32 * fs * 0.56).max(fs * 0.85) + margin_px * 1.1;
                            let kh = fs * 0.95 + margin_px * 0.8;
                            let kr = match anchor_mode {
                                egui::Align2::RIGHT_CENTER => Rect::from_min_size(Pos2::new(text_pos.x - kw + margin_px * 0.4, pos.y - kh * 0.5), Vec2::new(kw, kh)),
                                egui::Align2::LEFT_CENTER => Rect::from_min_size(Pos2::new(text_pos.x - margin_px * 0.4, pos.y - kh * 0.5), Vec2::new(kw, kh)),
                                _ => Rect::from_center_size(pos, Vec2::new(kw, kh)),
                            };
                            pp.rect_filled(kr, 2.0, canvas_bg);
                            let base_acol = if is_sel {
                                Color32::from_rgb(0, 110, 235)
                            } else if self.color_by_element {
                                cpk_element_color(a.z)
                            } else {
                                default_ink
                            };
                            let acol = scale_alpha(base_acol, atom_opacities.get(ai).copied().unwrap_or(1.0));
                            pp.text(text_pos, anchor_mode, &lbl, FontId::proportional(fs), acol);
                            if self.text_bold {
                                pp.text(text_pos + Vec2::new(0.6, 0.0), anchor_mode, &lbl, FontId::proportional(fs), acol);
                            }
                        }
                        if self.tab().show_atom_numbers {
                            pp.text(pos + Vec2::new(-7.0 * self.zoom, -6.5 * self.zoom), egui::Align2::RIGHT_BOTTOM, format!("{}", ai + 1), FontId::proportional(6.5 * self.zoom), Color32::from_gray(95));
                        }
                        if let Some(&map_num) = self.tab().reaction_atom_maps.get(&(mi, ai)) {
                            pp.text(pos + Vec2::new(6.0 * self.zoom, -6.0 * self.zoom), egui::Align2::LEFT_BOTTOM, format!("{map_num}"), FontId::proportional(6.5 * self.zoom), Color32::from_rgb(0, 80, 200));
                        }
                        if let Some(st) = &stereo {
                            if let Some(&code) = st.atom_stereo.get(&ai) {
                                let tag = match code { 2 => Some("(R)"), 3 => Some("(S)"), 4 => Some("(r)"), 5 => Some("(s)"), _ => None };
                                if let Some(t) = tag {
                                    pp.text(pos + Vec2::new(8.0 * self.zoom, -8.0 * self.zoom), egui::Align2::CENTER_CENTER, t, FontId::proportional(7.5 * self.zoom), Color32::from_rgb(0, 100, 0));
                                }
                            }
                        }
                        // Professional: Enhanced Stereochemistry (`abs`, `or1`, `&1`)
                        if let Some(&es) = self.tab().enhanced_stereo.get(&(mi, ai)) {
                            pp.text(
                                pos + Vec2::new(9.0 * self.zoom, 7.0 * self.zoom),
                                egui::Align2::LEFT_CENTER,
                                es.display_label(),
                                FontId::proportional(7.5 * self.zoom),
                                Color32::from_rgb(145, 25, 165),
                            );
                        }
                        // Live PAINS / Brenk structural alert atom highlight halo!
                        if self.pro_suite.show_medchem
                            && self.active_mol_index() == Some(mi)
                            && self.pro_suite.active_pains_atoms.contains(&ai)
                        {
                            pp.circle_filled(pos, 6.8 * self.zoom, Color32::from_rgba_unmultiplied(240, 75, 45, 95));
                            pp.circle_stroke(pos, 6.8 * self.zoom, Stroke::new(1.4_f32, Color32::from_rgb(210, 40, 30)));
                        }
                        // Live Multi-Spectroscopy (FT-IR / 2D NMR COSY/HSQC/HMBC) coupled atom highlight!
                        if self.omni_suite.show_spectroscopy
                            && self.active_mol_index() == Some(mi)
                            && self.omni_suite.highlighted_spec_atoms.contains(&ai)
                        {
                            pp.circle_filled(pos, 7.2 * self.zoom, Color32::from_rgba_unmultiplied(25, 185, 145, 105));
                            pp.circle_stroke(pos, 7.2 * self.zoom, Stroke::new(1.6_f32, Color32::from_rgb(15, 145, 110)));
                        }
                        // Professional: Live NMR Peak <-> Atom Highlight ring!
                        if self.show_nmr_window && self.nmr_mol_idx == mi && self.nmr_selected_atom == Some(ai) {
                            pp.circle_filled(pos, 7.0 * self.zoom, Color32::from_rgba_unmultiplied(255, 80, 40, 110));
                            pp.circle_stroke(pos, 7.0 * self.zoom, Stroke::new(1.8_f32, Color32::from_rgb(220, 20, 20)));
                        }
                        if let Some(aw) = warns.iter().find(|w| w.object_id == a.id) {
                            let wr = Rect::from_center_size(pos, Vec2::splat(11.0 * self.zoom));
                            pp.rect_stroke(wr, 0.0, Stroke::new(1.5_f32, Color32::RED), egui::StrokeKind::Outside);
                            if self.hovered_atom == Some((mi, ai)) {
                                let tip_pos = pos + Vec2::new(10.0, -18.0);
                                let galley = pp.layout_no_wrap(aw.message.clone(), FontId::proportional(11.0), Color32::BLACK);
                                let bg_r = Rect::from_min_size(tip_pos, galley.size() + Vec2::new(10.0, 6.0));
                                pp.rect_filled(bg_r, 2.0, Color32::from_rgb(255, 255, 225));
                                pp.rect_stroke(bg_r, 2.0, Stroke::new(1.0_f32, Color32::from_rgb(180, 50, 50)), egui::StrokeKind::Inside);
                                pp.galley(tip_pos + Vec2::new(5.0, 3.0), galley, Color32::BLACK);
                            }
                        }
                    }
                    // Professional: Reaction AutoNumber badge below molecule
                    if let Some(num_lbl) = self.tab().reaction_auto_numbers.get(&mi) {
                        if !mol.atoms.is_empty() {
                            let (x1, _, x2, y2) = mol_bbox(mol);
                            let npos = self.doc_to_screen(page_origin, 0.5 * (x1 + x2), y2 + 11.0);
                            pp.text(npos, egui::Align2::CENTER_TOP, num_lbl, FontId::proportional(9.5 * self.zoom), default_ink);
                        }
                    }
                }
                // Professional: Gel Electrophoresis Plates (GEP)
                for gep in &self.tab().gel_plates {
                    let r = Rect::from_two_pos(
                        self.doc_to_screen(page_origin, gep.bbox.0, gep.bbox.1),
                        self.doc_to_screen(page_origin, gep.bbox.2, gep.bbox.3),
                    );
                    let gel_bg = if self.ui_theme == UiTheme::DarkSlide {
                        Color32::from_rgb(36, 44, 58)
                    } else {
                        Color32::from_rgb(245, 248, 252)
                    };
                    pp.rect_filled(r, 1.0, gel_bg);
                    pp.rect_stroke(r, 1.0, Stroke::new(lw, default_ink), egui::StrokeKind::Outside);
                    let nlanes = gep.lanes.len().max(1) as f32;
                    let lane_w = r.width() / nlanes;
                    for (li, lane) in gep.lanes.iter().enumerate() {
                        let lx = r.min.x + (li as f32 + 0.5) * lane_w;
                        let well = Rect::from_center_size(Pos2::new(lx, r.min.y + 8.0 * self.zoom), Vec2::new(lane_w * 0.58, 4.0 * self.zoom));
                        pp.rect_filled(well, 0.5, Color32::from_gray(215));
                        pp.rect_stroke(well, 0.5, Stroke::new(0.9_f32, default_ink), egui::StrokeKind::Inside);
                        pp.text(Pos2::new(lx, r.min.y - 3.0 * self.zoom), egui::Align2::CENTER_BOTTOM, &lane.name, FontId::proportional(7.5 * self.zoom), default_ink);
                        let band_col = if self.ui_theme == UiTheme::DarkSlide {
                            Color32::from_rgb(130, 190, 255)
                        } else {
                            Color32::from_rgb(35, 55, 95)
                        };
                        for &pos_f in &lane.bands {
                            let by = r.min.y + (12.0 * self.zoom) + (pos_f as f32) * (r.height() - 16.0 * self.zoom);
                            let band_r = Rect::from_center_size(Pos2::new(lx, by), Vec2::new(lane_w * 0.62, 2.2 * self.zoom));
                            pp.rect_filled(band_r, 0.5, band_col);
                        }
                    }
                }

                // Signature ChiralCraft blue selection bounding box with 8 interactive resize handles + top rotation knob!
                if let Some((srect, hs, rot_knob)) = sel_handles_screen {
                    let blue = Color32::from_rgb(0, 120, 215);
                    pp.rect_stroke(srect, 0.0, Stroke::new(1.0_f32, blue), egui::StrokeKind::Outside);
                    for hp in hs {
                        let hr = Rect::from_center_size(hp, Vec2::splat(6.0));
                        pp.rect_filled(hr, 0.0, Color32::WHITE);
                        pp.rect_stroke(hr, 0.0, Stroke::new(1.0_f32, blue), egui::StrokeKind::Inside);
                    }
                    pp.line_segment([srect.center_top(), rot_knob], Stroke::new(1.0_f32, blue));
                    pp.circle_filled(rot_knob, 4.0, Color32::from_rgb(220, 240, 255));
                    pp.circle_stroke(rot_knob, 4.0, Stroke::new(1.2_f32, blue));
                }

                // Render Aromatic Circles, Styled Captions, Full-Featured Curves, Graphics, Tables, TLC Plates & Arrows
                for circ in &self.doc().scene.aromatic_circles {
                    pp.circle_stroke(self.doc_to_screen(page_origin, circ.center_x, circ.center_y), (circ.radius as f32) * self.zoom, Stroke::new(lw, default_ink));
                }
                for cap in &self.doc().scene.captions {
                    let cpos = self.doc_to_screen(page_origin, cap.x, cap.y);
                    let cfs = (self.doc().settings.caption_size_pt as f32) * self.zoom * 0.92;
                    let disp = if self.text_sub {
                        format_chemical_formula(&cap.text)
                    } else if self.text_sup {
                        cap.text.chars().map(|ch| match ch {
                            '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
                            '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
                            '+' => '⁺', '-' => '⁻', _ => ch,
                        }).collect()
                    } else if self.text_formula {
                        format_chemical_formula(&cap.text)
                    } else {
                        cap.text.clone()
                    };
                    let align = match cap.justification {
                        LabelJustification::Center | LabelJustification::Above | LabelJustification::Below => {
                            egui::Align2::CENTER_CENTER
                        }
                        LabelJustification::Right => egui::Align2::RIGHT_CENTER,
                        LabelJustification::Left | LabelJustification::Auto => match self.text_align {
                            TextAlign::Center | TextAlign::StackedAbove => egui::Align2::CENTER_CENTER,
                            TextAlign::Right => egui::Align2::RIGHT_CENTER,
                            TextAlign::Left | TextAlign::Justify | TextAlign::Automatic => egui::Align2::LEFT_CENTER,
                        },
                    };
                    let font = if self.text_italic { FontId::monospace(cfs) } else { FontId::proportional(cfs) };
                    let trect = pp.text(cpos, align, &disp, font.clone(), default_ink);
                    if self.text_bold {
                        pp.text(cpos + Vec2::new(0.7, 0.0), align, &disp, font, default_ink);
                    }
                    if self.text_underline {
                        pp.line_segment(
                            [Pos2::new(trect.min.x, trect.max.y), Pos2::new(trect.max.x, trect.max.y)],
                            Stroke::new(1.0_f32, default_ink),
                        );
                    }
                }
                for c in &self.doc().scene.curves {
                    let clw = if self.curve_bold { lw * 2.2 } else { lw };
                    let base_ccol = if self.active_color == Color32::BLACK {
                        default_ink
                    } else {
                        self.active_color
                    };
                    let ccol = if self.curve_faded {
                        Color32::from_rgba_unmultiplied(base_ccol.r(), base_ccol.g(), base_ccol.b(), 110)
                    } else {
                        base_ccol
                    };
                    let spts: Vec<Pos2> = c.points.iter().map(|&(x, y)| self.doc_to_screen(page_origin, x, y)).collect();
                    if (self.curve_filled || self.curve_shaded) && spts.len() >= 3 {
                        let fill = if self.curve_shaded {
                            Color32::from_rgba_unmultiplied(120, 160, 210, 150)
                        } else {
                            Color32::from_rgba_unmultiplied(185, 209, 234, 120)
                        };
                        pp.add(egui::Shape::convex_polygon(spts.clone(), fill, Stroke::new(clw, ccol)));
                    } else {
                        let draw_seg = |a: Pos2, b: Pos2| {
                            let d = b - a;
                            let slen = d.length().max(1e-3);
                            let sn = Vec2::new(-d.y / slen, d.x / slen);
                            let offsets: &[f32] = if self.curve_doubled { &[-1.3 * self.zoom, 1.3 * self.zoom] } else { &[0.0] };
                            for &off in offsets {
                                let pa = a + sn * off;
                                let pb = b + sn * off;
                                if self.curve_dashed {
                                    let raw_steps = if slen.is_finite() { (slen / (6.0 * self.zoom)).round() } else { 2.0 };
                                    let steps = (raw_steps.clamp(2.0, 256.0)) as usize;
                                    for k in 0..steps {
                                        let t0 = (k as f32) / (steps as f32);
                                        let t1 = ((k as f32) + 0.55) / (steps as f32);
                                        pp.line_segment([pa + d * t0, pa + d * t1.min(1.0)], Stroke::new(clw, ccol));
                                    }
                                } else {
                                    pp.line_segment([pa, pb], Stroke::new(clw, ccol));
                                }
                            }
                        };
                        for w in spts.windows(2) {
                            draw_seg(w[0], w[1]);
                        }
                        if c.closed {
                            if let (Some(&first), Some(&last)) = (spts.first(), spts.last()) {
                                draw_seg(last, first);
                            }
                        }
                    }
                    if (self.curve_arrow_end || self.curve_half_end != 0) && spts.len() >= 2 {
                        let p2 = spts[spts.len() - 1];
                        let p1 = spts[spts.len() - 2];
                        let u = (p2 - p1).normalized();
                        let n = Vec2::new(-u.y, u.x);
                        let hl = 5.5 * self.zoom;
                        let hw = 2.2 * self.zoom;
                        if self.curve_arrow_end {
                            pp.add(egui::Shape::convex_polygon(vec![p2, p2 - u * hl + n * hw, p2 - u * hl - n * hw], ccol, Stroke::NONE));
                        } else {
                            let side = self.curve_half_end as f32;
                            pp.line_segment([p2, p2 - u * hl + n * hw * side], Stroke::new(clw, ccol));
                        }
                    }
                    if (self.curve_arrow_start || self.curve_half_start != 0) && spts.len() >= 2 {
                        let p1 = spts[0];
                        let p2 = spts[1];
                        let u = (p1 - p2).normalized();
                        let n = Vec2::new(-u.y, u.x);
                        let hl = 5.5 * self.zoom;
                        let hw = 2.2 * self.zoom;
                        if self.curve_arrow_start {
                            pp.add(egui::Shape::convex_polygon(vec![p1, p1 - u * hl + n * hw, p1 - u * hl - n * hw], ccol, Stroke::NONE));
                        } else {
                            let side = self.curve_half_start as f32;
                            pp.line_segment([p1, p1 - u * hl + n * hw * side], Stroke::new(clw, ccol));
                        }
                    }
                }
                // View > Show Reaction Interpretation visual overlay (Reactants green, Products blue)
                if self.show_reaction_interp && !self.doc().scene.arrows.is_empty() {
                    let split_x = 0.5 * (self.doc().scene.arrows[0].tail.0 + self.doc().scene.arrows[0].head.0);
                    for m in &self.doc().scene.molecules {
                        if m.atoms.is_empty() { continue; }
                        let (x1, y1, x2, y2) = mol_bbox(m);
                        let cx = 0.5 * (x1 + x2);
                        let r = Rect::from_two_pos(
                            self.doc_to_screen(page_origin, x1 - 4.0, y1 - 4.0),
                            self.doc_to_screen(page_origin, x2 + 4.0, y2 + 4.0),
                        );
                        let (col, role) = if cx <= split_x {
                            (Color32::from_rgb(0, 150, 60), "Reactant")
                        } else {
                            (Color32::from_rgb(0, 95, 200), "Product")
                        };
                        pp.rect_stroke(r, 2.0, Stroke::new(1.2_f32, col), egui::StrokeKind::Outside);
                        pp.text(r.left_top() - Vec2::new(0.0, 4.0), egui::Align2::LEFT_BOTTOM, role, FontId::proportional(9.5), col);
                    }
                }
                for g in &self.doc().scene.graphics {
                    let p1 = self.doc_to_screen(page_origin, g.bbox.0, g.bbox.1);
                    let p2 = self.doc_to_screen(page_origin, g.bbox.2, g.bbox.3);
                    let r = Rect::from_two_pos(p1, p2);
                    match g.graphic_type.as_str() {
                        "CleavageLine" => {
                            pp.line_segment([p1, p2], Stroke::new(1.5_f32, Color32::RED));
                        }
                        "Line" => {
                            pp.line_segment([p1, p2], Stroke::new(lw, default_ink));
                        }
                        "Oval" => {
                            pp.circle_stroke(r.center(), r.width().min(r.height()) * 0.5, Stroke::new(lw, default_ink));
                        }
                        "MultiCenter" => {
                            pp.circle_filled(r.center(), 2.8 * self.zoom, Color32::from_rgb(0, 110, 210));
                        }
                        "VariableAttachment" => {
                            pp.line_segment([p1, p2], Stroke::new(lw, Color32::from_rgb(0, 110, 210)));
                            pp.circle_filled(p1, 2.2 * self.zoom, Color32::from_rgb(0, 110, 210));
                        }
                        "Orbital" => {
                            let c = r.center();
                            match g.sub_type.as_deref() {
                                Some("s") => { pp.circle_stroke(c, r.width() * 0.7, Stroke::new(lw, default_ink)); }
                                Some("d") => {
                                    for (ox, oy, f) in [(-0.4, -0.2, true), (0.4, 0.2, true), (-0.4, 0.2, false), (0.4, -0.2, false)] {
                                        let pc = Pos2::new(c.x + r.width() * ox, c.y + r.height() * oy);
                                        if f { pp.circle_filled(pc, r.width() * 0.35, default_ink); } else { pp.circle_stroke(pc, r.width() * 0.35, Stroke::new(lw, default_ink)); }
                                    }
                                }
                                _ => {
                                    pp.circle_filled(Pos2::new(c.x, c.y - r.height() * 0.25), r.width() * 0.45, default_ink);
                                    pp.circle_stroke(Pos2::new(c.x, c.y + r.height() * 0.25), r.width() * 0.45, Stroke::new(lw, default_ink));
                                }
                            }
                        }
                        "Symbol" => {
                            let c = r.center();
                            let rad = r.width().min(r.height()) * 0.5;
                            match g.sub_type.as_deref() {
                                Some("CircleMinus") => {
                                    pp.circle_stroke(c, rad, Stroke::new(lw, default_ink));
                                    pp.line_segment([Pos2::new(c.x - rad * 0.5, c.y), Pos2::new(c.x + rad * 0.5, c.y)], Stroke::new(lw, default_ink));
                                }
                                Some("Radical") => { pp.circle_filled(c, rad * 0.35, default_ink); }
                                Some("LonePair") => {
                                    pp.circle_filled(Pos2::new(c.x - rad * 0.35, c.y), rad * 0.25, default_ink);
                                    pp.circle_filled(Pos2::new(c.x + rad * 0.35, c.y), rad * 0.25, default_ink);
                                }
                                _ => {
                                    pp.circle_stroke(c, rad, Stroke::new(lw, default_ink));
                                    pp.line_segment([Pos2::new(c.x - rad * 0.5, c.y), Pos2::new(c.x + rad * 0.5, c.y)], Stroke::new(lw, default_ink));
                                    pp.line_segment([Pos2::new(c.x, c.y - rad * 0.5), Pos2::new(c.x, c.y + rad * 0.5)], Stroke::new(lw, default_ink));
                                }
                            }
                        }
                        "Bracket" => {
                            let arm = 6.0 * self.zoom;
                            match g.sub_type.as_deref() {
                                Some("Round") => {
                                    // Parentheses ( )
                                    for side in [-1.0_f32, 1.0_f32] {
                                        let bx = if side < 0.0 { r.min.x } else { r.max.x };
                                        let mut prev = Pos2::new(bx - side * arm, r.min.y);
                                        for k in 1..=10 {
                                            let t = (k as f32) / 10.0;
                                            let bulge = (t * std::f32::consts::PI).sin();
                                            let cur = Pos2::new(bx - side * arm * (1.0 - bulge), r.min.y + t * r.height());
                                            pp.line_segment([prev, cur], Stroke::new(lw, default_ink));
                                            prev = cur;
                                        }
                                    }
                                }
                                Some("Curly") => {
                                    // Braces { }
                                    for side in [-1.0_f32, 1.0_f32] {
                                        let bx = if side < 0.0 { r.min.x } else { r.max.x };
                                        let my = r.center().y;
                                        let pts = [
                                            Pos2::new(bx - side * arm, r.min.y),
                                            Pos2::new(bx - side * (arm * 0.35), r.min.y + r.height() * 0.12),
                                            Pos2::new(bx - side * (arm * 0.35), my - 4.0 * self.zoom),
                                            Pos2::new(bx + side * (arm * 0.25), my),
                                            Pos2::new(bx - side * (arm * 0.35), my + 4.0 * self.zoom),
                                            Pos2::new(bx - side * (arm * 0.35), r.max.y - r.height() * 0.12),
                                            Pos2::new(bx - side * arm, r.max.y),
                                        ];
                                        for w in pts.windows(2) {
                                            pp.line_segment([w[0], w[1]], Stroke::new(lw, default_ink));
                                        }
                                    }
                                }
                                _ => {
                                    pp.line_segment([r.min + Vec2::new(arm, 0.0), r.min], Stroke::new(lw, default_ink));
                                    pp.line_segment([r.min, Pos2::new(r.min.x, r.max.y)], Stroke::new(lw, default_ink));
                                    pp.line_segment([Pos2::new(r.min.x, r.max.y), Pos2::new(r.min.x + arm, r.max.y)], Stroke::new(lw, default_ink));
                                    pp.line_segment([Pos2::new(r.max.x - arm, r.min.y), Pos2::new(r.max.x, r.min.y)], Stroke::new(lw, default_ink));
                                    pp.line_segment([Pos2::new(r.max.x, r.min.y), r.max], Stroke::new(lw, default_ink));
                                    pp.line_segment([r.max, r.max - Vec2::new(arm, 0.0)], Stroke::new(lw, default_ink));
                                }
                            }
                        }
                        _ => {
                            let rad = match g.sub_type.as_deref() {
                                Some("Rounded") | Some("RoundedShadow") => 6.0 * self.zoom,
                                _ => 0.0,
                            };
                            let has_shadow = matches!(g.sub_type.as_deref(), Some("Shadow") | Some("RoundedShadow") | None);
                            if has_shadow {
                                let sh_col = if self.ui_theme == UiTheme::DarkSlide { Color32::from_gray(70) } else { Color32::from_gray(185) };
                                pp.rect_filled(r.translate(Vec2::splat(2.2 * self.zoom)), rad, sh_col);
                            }
                            pp.rect_stroke(r, rad, Stroke::new(lw, default_ink), egui::StrokeKind::Outside);
                        }
                    }
                }
                for tbl in &self.doc().scene.tables {
                    let r = Rect::from_two_pos(self.doc_to_screen(page_origin, tbl.bbox.0, tbl.bbox.1), self.doc_to_screen(page_origin, tbl.bbox.2, tbl.bbox.3));
                    pp.rect_stroke(r, 0.0, Stroke::new(lw, default_ink), egui::StrokeKind::Outside);
                    let cols = tbl.cols.max(1);
                    let rows = tbl.rows.max(1);
                    let cw = r.width() / (cols as f32);
                    let rh = r.height() / (rows as f32);
                    for c in 1..cols {
                        let x = r.min.x + cw * (c as f32);
                        pp.line_segment([Pos2::new(x, r.min.y), Pos2::new(x, r.max.y)], Stroke::new(lw, default_ink));
                    }
                    for row in 1..rows {
                        let y = r.min.y + rh * (row as f32);
                        pp.line_segment([Pos2::new(r.min.x, y), Pos2::new(r.max.x, y)], Stroke::new(lw, default_ink));
                    }
                    let cell_fs = (8.0 * self.zoom).clamp(7.0, 18.0);
                    for (ri, row_vals) in tbl.cells.iter().enumerate().take(rows) {
                        for (ci, txt) in row_vals.iter().enumerate().take(cols) {
                            if !txt.is_empty() {
                                let cx = r.min.x + (ci as f32 + 0.5) * cw;
                                let cy = r.min.y + (ri as f32 + 0.5) * rh;
                                pp.text(Pos2::new(cx, cy), egui::Align2::CENTER_CENTER, txt, FontId::proportional(cell_fs), default_ink);
                            }
                        }
                    }
                }
                for tlc in &self.doc().scene.tlc_plates {
                    let r = Rect::from_two_pos(self.doc_to_screen(page_origin, tlc.bbox.0, tlc.bbox.1), self.doc_to_screen(page_origin, tlc.bbox.2, tlc.bbox.3));
                    pp.rect_stroke(r, 0.0, Stroke::new(lw, default_ink), egui::StrokeKind::Outside);
                    let oy = r.max.y - r.height() * (tlc.origin_fraction as f32);
                    let fy = r.min.y + r.height() * 0.08;
                    pp.line_segment([Pos2::new(r.min.x, oy), Pos2::new(r.max.x, oy)], Stroke::new(1.0_f32, Color32::GRAY));
                    pp.line_segment([Pos2::new(r.min.x, fy), Pos2::new(r.max.x, fy)], Stroke::new(0.8_f32, Color32::from_gray(150)));
                    let lbl_fs = (6.8 * self.zoom).clamp(6.5, 14.0);
                    for (i, &rf) in tlc.spot_rfs.iter().enumerate() {
                        let sx = r.min.x + r.width() * ((i + 1) as f32) / ((tlc.spot_rfs.len() + 1) as f32);
                        let sy = oy - (rf as f32) * (oy - fy);
                        let scol = tlc.spot_colors.get(i).copied().map(|(cr, cg, cb)| {
                            if self.ui_theme == UiTheme::DarkSlide && cr < 40 && cg < 40 && cb < 40 {
                                default_ink
                            } else {
                                Color32::from_rgb(cr, cg, cb)
                            }
                        }).unwrap_or(default_ink);
                        pp.circle_filled(Pos2::new(sx, sy), 2.7 * self.zoom, scol);
                        if let Some(lbl) = tlc.spot_labels.get(i) {
                            if !lbl.is_empty() {
                                pp.text(Pos2::new(sx, r.max.y + 3.0 * self.zoom), egui::Align2::CENTER_TOP, lbl, FontId::proportional(lbl_fs), default_ink);
                            }
                        }
                    }
                    if !tlc.solvent_label.is_empty() {
                        pp.text(Pos2::new(r.center().x, r.min.y - 3.0 * self.zoom), egui::Align2::CENTER_BOTTOM, &tlc.solvent_label, FontId::proportional(lbl_fs), Color32::from_gray(120));
                    }
                }
                // Render Reaction Energy Profiles (ΔG‡ free-energy coordinate diagrams) live on canvas!
                for ep in &self.doc().scene.energy_profiles {
                    let r = Rect::from_two_pos(
                        self.doc_to_screen(page_origin, ep.bbox.0, ep.bbox.1),
                        self.doc_to_screen(page_origin, ep.bbox.2, ep.bbox.3),
                    );
                    let mut min_e = f64::INFINITY;
                    let mut max_e = f64::NEG_INFINITY;
                    for st in ep.pathway_a.iter().chain(ep.pathway_b.iter()) {
                        min_e = min_e.min(st.energy);
                        max_e = max_e.max(st.energy);
                    }
                    if !min_e.is_finite() { continue; }
                    let pad_e = ((max_e - min_e) * 0.18).max(3.0);
                    min_e -= pad_e;
                    max_e += pad_e;
                    let e_span = (max_e - min_e).max(1.0);
                    let pl = r.min.x + 24.0 * self.zoom;
                    let pr = r.max.x - 8.0 * self.zoom;
                    let pt = r.min.y + 14.0 * self.zoom;
                    let pb = r.max.y - 16.0 * self.zoom;
                    pp.line_segment([Pos2::new(pl, pt), Pos2::new(pl, pb)], Stroke::new(lw, default_ink));
                    pp.line_segment([Pos2::new(pl, pb), Pos2::new(pr, pb)], Stroke::new(lw, default_ink));
                    pp.text(
                        Pos2::new(r.min.x + 2.0 * self.zoom, pt - 3.0 * self.zoom),
                        egui::Align2::LEFT_BOTTOM,
                        format!("ΔG ({})", ep.unit_str()),
                        FontId::proportional((6.5 * self.zoom).clamp(7.0, 14.0)),
                        default_ink,
                    );
                    let draw_ep = |states: &[chiralcraft_doc::EnergyState], other: &[chiralcraft_doc::EnergyState], col: Color32, dashed: bool, is_secondary: bool| {
                        let n = states.len();
                        if n == 0 { return; }
                        let step = (pr - pl) / (n as f32);
                        let hw = (step * 0.24).clamp(4.0 * self.zoom, 14.0 * self.zoom);
                        let pts: Vec<(Pos2, &chiralcraft_doc::EnergyState)> = states
                            .iter()
                            .enumerate()
                            .map(|(i, st)| {
                                let cx = pl + (i as f32 + 0.5) * step;
                                let cy = pb - (((st.energy - min_e) / e_span) as f32) * (pb - pt);
                                (Pos2::new(cx, cy), st)
                            })
                            .collect();
                        for pair in pts.windows(2) {
                            let (p0, _) = pair[0];
                            let (p1, _) = pair[1];
                            let sx = p0.x + hw;
                            let ex = p1.x - hw;
                            let segs = 12usize;
                            let mut prev = Pos2::new(sx, p0.y);
                            for k in 1..=segs {
                                let t = (k as f32) / (segs as f32);
                                let smooth = t * t * (3.0 - 2.0 * t);
                                let cur = Pos2::new(sx + (ex - sx) * t, p0.y + (p1.y - p0.y) * smooth);
                                if !dashed || k % 2 == 1 {
                                    pp.line_segment([prev, cur], Stroke::new(lw * 1.1, col));
                                }
                                prev = cur;
                            }
                        }
                        let fs = (6.2 * self.zoom).clamp(6.5, 13.0);
                        for (i, (c, st)) in pts.into_iter().enumerate() {
                            pp.line_segment([Pos2::new(c.x - hw, c.y), Pos2::new(c.x + hw, c.y)], Stroke::new((lw * 2.4).max(2.2), col));
                            if let Some(ost) = other.get(i) {
                                if (st.energy - ost.energy).abs() < 0.05 && st.label == ost.label {
                                    if is_secondary {
                                        continue;
                                    }
                                    pp.text(Pos2::new(c.x, c.y - 3.0 * self.zoom), egui::Align2::CENTER_BOTTOM, format!("{} ({:.1})", st.label, st.energy), FontId::proportional(fs), col);
                                    continue;
                                }
                                if st.energy < ost.energy {
                                    pp.text(Pos2::new(c.x, c.y + 3.5 * self.zoom), egui::Align2::CENTER_TOP, format!("{} ({:.1})", st.label, st.energy), FontId::proportional(fs), col);
                                    continue;
                                }
                            }
                            pp.text(Pos2::new(c.x, c.y - 3.0 * self.zoom), egui::Align2::CENTER_BOTTOM, format!("{} ({:.1})", st.label, st.energy), FontId::proportional(fs), col);
                        }
                    };
                    draw_ep(&ep.pathway_a, &ep.pathway_b, Color32::from_rgb(22, 87, 184), false, false);
                    draw_ep(&ep.pathway_b, &ep.pathway_a, Color32::from_rgb(200, 50, 50), true, true);
                }
                for arr in &self.doc().scene.arrows {
                    let p1 = self.doc_to_screen(page_origin, arr.tail.0, arr.tail.1);
                    let p2 = self.doc_to_screen(page_origin, arr.head.0, arr.head.1);
                    let d = p2 - p1;
                    let len = d.length().max(1e-3);
                    let u = d / len;
                    let n = Vec2::new(-u.y, u.x);
                    let hl = 5.5 * self.zoom;
                    let hw = 2.2 * self.zoom;

                    if arr.arrowhead_tail.is_some() {
                        // Equilibrium double-shafted antiparallel harpoons
                        let sep = 1.8 * self.zoom;
                        pp.line_segment([p1 + n * sep, p2 + n * sep], Stroke::new(lw, default_ink));
                        pp.line_segment([p2 + n * sep, p2 + n * sep - u * hl + n * hw], Stroke::new(lw, default_ink));
                        pp.line_segment([p1 - n * sep, p2 - n * sep], Stroke::new(lw, default_ink));
                        pp.line_segment([p1 - n * sep, p1 - n * sep + u * hl - n * hw], Stroke::new(lw, default_ink));
                    } else if arr.arrowhead_head == "Retrosynthetic" {
                        let sep = 1.6 * self.zoom;
                        pp.line_segment([p1 + n * sep, p2 - u * 3.0 * self.zoom + n * sep], Stroke::new(lw, default_ink));
                        pp.line_segment([p1 - n * sep, p2 - u * 3.0 * self.zoom - n * sep], Stroke::new(lw, default_ink));
                        pp.line_segment([p2, p2 - u * hl + n * hw * 1.6], Stroke::new(lw, default_ink));
                        pp.line_segment([p2, p2 - u * hl - n * hw * 1.6], Stroke::new(lw, default_ink));
                    } else if arr.arrowhead_head == "Curved" {
                        let ctrl = p1 + d * 0.5 - n * (len * 0.35);
                        let mut prev = p1;
                        for k in 1..=14 {
                            let t = (k as f32) / 14.0;
                            let mt = 1.0 - t;
                            let cur = Pos2::new(
                                mt * mt * p1.x + 2.0 * mt * t * ctrl.x + t * t * p2.x,
                                mt * mt * p1.y + 2.0 * mt * t * ctrl.y + t * t * p2.y,
                            );
                            pp.line_segment([prev, cur], Stroke::new(lw, default_ink));
                            prev = cur;
                        }
                        let tu = (p2 - ctrl).normalized();
                        let tn = Vec2::new(-tu.y, tu.x);
                        pp.add(egui::Shape::convex_polygon(vec![p2, p2 - tu * hl + tn * hw, p2 - tu * hl - tn * hw], default_ink, Stroke::NONE));
                    } else {
                        pp.line_segment([p1, p2], Stroke::new(lw, default_ink));
                        pp.add(egui::Shape::convex_polygon(vec![p2, p2 - u * hl + n * hw, p2 - u * hl - n * hw], default_ink, Stroke::NONE));
                    }
                    if arr.no_go.as_deref() == Some("Cross") {
                        let m = p1 + d * 0.5;
                        let cr = 4.0 * self.zoom;
                        pp.line_segment([m - Vec2::new(cr, cr), m + Vec2::new(cr, cr)], Stroke::new(lw * 1.3, Color32::RED));
                        pp.line_segment([m - Vec2::new(cr, -cr), m + Vec2::new(cr, -cr)], Stroke::new(lw * 1.3, Color32::RED));
                    }
                }

                // View > Show Crosshair guide lines
                if self.show_crosshair {
                    if let Some(pos) = response.hover_pos() {
                        if view_rect.contains(pos) {
                            let cstroke = Stroke::new(0.8_f32, Color32::from_rgba_unmultiplied(100, 140, 200, 140));
                            pp.line_segment([Pos2::new(view_rect.min.x, pos.y), Pos2::new(view_rect.max.x, pos.y)], cstroke);
                            pp.line_segment([Pos2::new(pos.x, view_rect.min.y), Pos2::new(pos.x, view_rect.max.y)], cstroke);
                        }
                    }
                }

                // Live hover ghost preview (for Ring tools and Bond-tool atom hover before clicking!)
                if self.drag_start_canvas.is_none() && self.mdi_drag == MdiDragMode::None {
                    let blen = self.doc().settings.bond_length_pt;
                    let ghost_col = Color32::from_rgba_unmultiplied(0, 120, 225, 115);
                    if let Some((n_sides, _, _)) = Self::ring_spec(self.tool) {
                        if let Some((mi, ai)) = self.hovered_atom {
                            if let Some(mol) = self.doc().scene.molecules.get(mi) {
                                let verts = Self::compute_spiro_ring_vertices(mol, ai, n_sides, blen, None);
                                let spts: Vec<Pos2> = verts.iter().map(|&(x, y)| self.doc_to_screen(page_origin, x, y)).collect();
                                for i in 0..spts.len() {
                                    pp.line_segment([spts[i], spts[(i + 1) % spts.len()]], Stroke::new(1.3_f32, ghost_col));
                                }
                            }
                        } else if let Some((mi, bi)) = self.hovered_bond {
                            if let Some(mol) = self.doc().scene.molecules.get(mi) {
                                let verts = Self::compute_edge_fused_ring_vertices(mol, bi, n_sides, blen, None);
                                let spts: Vec<Pos2> = verts.iter().map(|&(x, y)| self.doc_to_screen(page_origin, x, y)).collect();
                                for i in 0..spts.len() {
                                    pp.line_segment([spts[i], spts[(i + 1) % spts.len()]], Stroke::new(1.3_f32, ghost_col));
                                }
                            }
                        }
                    } else if self.is_bond_tool() {
                        if let Some((mi, ai)) = self.hovered_atom {
                            if let Some(mol) = self.doc().scene.molecules.get(mi) {
                                if let Some(a) = mol.atoms.get(ai) {
                                    let (ord, _) = self.current_bond_spec();
                                    let ang = Self::best_new_bond_angle(mol, ai, ord, self.doc().settings.chain_angle_deg.to_radians());
                                    let sp1 = self.doc_to_screen(page_origin, a.x, a.y);
                                    let sp2 = self.doc_to_screen(page_origin, a.x + blen * ang.cos(), a.y + blen * ang.sin());
                                    pp.line_segment([sp1, sp2], Stroke::new(1.2_f32, Color32::from_rgba_unmultiplied(0, 120, 225, 95)));
                                }
                            }
                        }
                    }
                }

                // Live marquee / freehand lasso / single-atom / molecule & caption move ghost / ring drag / eraser / acyclic chain / snapped bond drag preview
                if let (Some((sx, sy)), Some((ex, ey))) = (self.drag_start_canvas, self.drag_cur_canvas) {
                    let p1 = self.doc_to_screen(page_origin, sx, sy);
                    let p2 = self.doc_to_screen(page_origin, ex, ey);
                    let blen = self.doc().settings.bond_length_pt;
                    let atol = self.atom_hit_tol();
                    let btol = self.bond_hit_tol();
                    if let Some((mi, ai)) = self.dragging_single_atom {
                        let ghost = Color32::from_rgba_unmultiplied(0, 120, 225, 195);
                        if let Some(mol) = self.doc().scene.molecules.get(mi) {
                            for b in &mol.bonds {
                                let other = if b.begin_atom == ai {
                                    Some(b.end_atom)
                                } else if b.end_atom == ai {
                                    Some(b.begin_atom)
                                } else {
                                    None
                                };
                                if let Some(oi) = other.and_then(|o| mol.atoms.get(o)) {
                                    pp.line_segment([self.doc_to_screen(page_origin, oi.x, oi.y), p2], Stroke::new(lw.max(1.5), ghost));
                                }
                            }
                            pp.circle_stroke(p2, 5.0 * self.zoom, Stroke::new(1.4_f32, ghost));
                        }
                    } else if self.dragging_selected_move && !self.tab().selected_mols.is_empty() {
                        let ddx = ex - sx;
                        let ddy = ey - sy;
                        let ghost = Color32::from_rgba_unmultiplied(0, 110, 215, 165);
                        for &mi in &self.tab().selected_mols {
                            if let Some(m) = self.doc().scene.molecules.get(mi) {
                                for b in &m.bonds {
                                    if let (Some(a1), Some(a2)) = (m.atoms.get(b.begin_atom), m.atoms.get(b.end_atom)) {
                                        pp.line_segment(
                                            [
                                                self.doc_to_screen(page_origin, a1.x + ddx, a1.y + ddy),
                                                self.doc_to_screen(page_origin, a2.x + ddx, a2.y + ddy),
                                            ],
                                            Stroke::new(lw, ghost),
                                        );
                                    }
                                }
                            }
                        }
                        if let Some((bx1, by1, bx2, by2)) = self.selection_bbox() {
                            let gr = Rect::from_two_pos(
                                self.doc_to_screen(page_origin, bx1 + ddx - 5.0, by1 + ddy - 5.0),
                                self.doc_to_screen(page_origin, bx2 + ddx + 5.0, by2 + ddy + 5.0),
                            );
                            pp.rect_stroke(gr, 0.0, Stroke::new(1.0_f32, ghost), egui::StrokeKind::Outside);
                        }
                    } else if let Some(ci) = self.dragging_caption {
                        if let Some(cap) = self.doc().scene.captions.get(ci) {
                            let gp = self.doc_to_screen(page_origin, cap.x + (ex - sx), cap.y + (ey - sy));
                            let galign = match cap.justification {
                                LabelJustification::Center | LabelJustification::Above | LabelJustification::Below => {
                                    egui::Align2::CENTER_CENTER
                                }
                                LabelJustification::Right => egui::Align2::RIGHT_CENTER,
                                LabelJustification::Left | LabelJustification::Auto => egui::Align2::LEFT_CENTER,
                            };
                            pp.text(
                                gp,
                                galign,
                                &cap.text,
                                FontId::proportional((self.doc().settings.caption_size_pt as f32) * self.zoom * 0.92),
                                Color32::from_rgba_unmultiplied(0, 110, 215, 175),
                            );
                        }
                    } else if self.tool == Tool::Lasso && !self.dragging_selected_move && self.lasso_pts.len() >= 2 {
                        let spts: Vec<Pos2> = self.lasso_pts.iter().map(|&(x, y)| self.doc_to_screen(page_origin, x, y)).collect();
                        for w in spts.windows(2) {
                            pp.line_segment([w[0], w[1]], Stroke::new(1.2_f32, Color32::from_rgb(0, 120, 215)));
                        }
                        if let (Some(&first), Some(&last)) = (spts.first(), spts.last()) {
                            pp.line_segment([last, first], Stroke::new(1.0_f32, Color32::from_rgba_unmultiplied(0, 120, 215, 130)));
                        }
                    } else if self.tool == Tool::Marquee && !self.dragging_selected_move {
                        pp.rect_stroke(Rect::from_two_pos(p1, p2), 0.0, Stroke::new(1.0_f32, Color32::from_rgb(0, 120, 215)), egui::StrokeKind::Outside);
                    } else if self.tool == Tool::Fragmentation {
                        pp.line_segment([p1, p2], Stroke::new(1.6_f32, Color32::RED));
                    } else if self.tool == Tool::Eraser && self.eraser_pts.len() >= 2 {
                        let spts: Vec<Pos2> = self.eraser_pts.iter().map(|&(x, y)| self.doc_to_screen(page_origin, x, y)).collect();
                        for w in spts.windows(2) {
                            pp.line_segment([w[0], w[1]], Stroke::new(2.0_f32, Color32::from_rgba_unmultiplied(220, 45, 45, 165)));
                        }
                    } else if let Some((n_sides, _, _)) = Self::ring_spec(self.tool) {
                        let ring_col = Color32::from_rgb(0, 110, 220);
                        let verts = if let Some((mi, ai)) = self.hit_atom(sx, sy, atol) {
                            self.doc().scene.molecules.get(mi).map(|m| Self::compute_spiro_ring_vertices(m, ai, n_sides, blen, Some((ex, ey))))
                        } else if let Some((mi, bi)) = self.hit_bond(sx, sy, btol) {
                            self.doc().scene.molecules.get(mi).map(|m| Self::compute_edge_fused_ring_vertices(m, bi, n_sides, blen, Some((ex, ey))))
                        } else {
                            let rot = if (ex - sx).hypot(ey - sy) > 4.0 { (ey - sy).atan2(ex - sx) } else { -std::f64::consts::FRAC_PI_2 };
                            let m = Self::build_polygon_ring_angled(n_sides, false, false, sx, sy, blen, rot);
                            Some(m.atoms.iter().map(|a| (a.x, a.y)).collect())
                        };
                        if let Some(vs) = verts {
                            let spts: Vec<Pos2> = vs.iter().map(|&(x, y)| self.doc_to_screen(page_origin, x, y)).collect();
                            for i in 0..spts.len() {
                                pp.line_segment([spts[i], spts[(i + 1) % spts.len()]], Stroke::new(1.5_f32, ring_col));
                            }
                        }
                    } else if self.tool == Tool::AcyclicChain {
                        let (ax, ay) = self
                            .drag_start_atom
                            .and_then(|(mi, ai)| self.doc().scene.molecules.get(mi).and_then(|m| m.atoms.get(ai)))
                            .map(|a| (a.x, a.y))
                            .unwrap_or((sx, sy));
                        let prev = Self::build_acyclic_chain(ax, ay, ex, ey, blen);
                        for b in &prev.bonds {
                            let a1 = &prev.atoms[b.begin_atom];
                            let a2 = &prev.atoms[b.end_atom];
                            pp.line_segment([self.doc_to_screen(page_origin, a1.x, a1.y), self.doc_to_screen(page_origin, a2.x, a2.y)], Stroke::new(1.2_f32, Color32::from_rgb(0, 120, 215)));
                        }
                        pp.text(p2 + Vec2::new(10.0, -10.0), egui::Align2::LEFT_CENTER, format!("{} C", prev.atoms.len()), FontId::proportional(12.0), Color32::from_rgb(0, 90, 210));
                    } else if self.is_bond_tool() {
                        let (ax, ay) = self
                            .drag_start_atom
                            .and_then(|(mi, ai)| self.doc().scene.molecules.get(mi).and_then(|m| m.atoms.get(ai)))
                            .map(|a| (a.x, a.y))
                            .unwrap_or((sx, sy));
                        let mut end_snap = self.hit_atom(ex, ey, atol.max(6.5));
                        let (mut nex, mut ney) = if let Some((emi, eai)) = end_snap {
                            self.doc()
                                .scene
                                .molecules
                                .get(emi)
                                .and_then(|m| m.atoms.get(eai))
                                .map(|a| (a.x, a.y))
                                .unwrap_or((ex, ey))
                        } else {
                            snap_bond_endpoint(
                                ax, ay, ex, ey,
                                self.doc().settings.fixed_lengths,
                                self.doc().settings.fixed_angles,
                                blen,
                            )
                        };
                        if end_snap.is_none() {
                            if let Some((emi, eai)) = self.hit_atom(nex, ney, (blen * 0.38).max(5.5)) {
                                if self.drag_start_atom != Some((emi, eai)) {
                                    end_snap = Some((emi, eai));
                                    if let Some(a) = self.doc().scene.molecules.get(emi).and_then(|m| m.atoms.get(eai)) {
                                        nex = a.x;
                                        ney = a.y;
                                    }
                                }
                            }
                        }
                        let sp1 = self.doc_to_screen(page_origin, ax, ay);
                        let sp2 = self.doc_to_screen(page_origin, nex, ney);
                        pp.line_segment([sp1, sp2], Stroke::new(1.5_f32, Color32::from_rgb(0, 100, 210)));
                        if end_snap.is_some() {
                            pp.circle_stroke(sp2, 5.5 * self.zoom, Stroke::new(1.3_f32, Color32::from_rgb(0, 120, 215)));
                        }
                    } else {
                        pp.line_segment([p1, p2], Stroke::new(1.2_f32, Color32::GRAY));
                    }
                }

                // Inline On-Canvas Atom Label Editing Box (signature ChiralCraft UX!)
                if let Some((mi, ai)) = self.inline_atom_edit {
                    if let Some(a) = self.doc().scene.molecules.get(mi).and_then(|m| m.atoms.get(ai)) {
                        let spos = self.doc_to_screen(page_origin, a.x, a.y);
                        let mut commit = false;
                        let mut cancel = false;
                        egui::Area::new(egui::Id::new("inline_atom_box"))
                            .fixed_pos(spos - Vec2::new(26.0, 11.0))
                            .order(egui::Order::Foreground)
                            .show(ctx, |ui| {
                                let te = ui.add(
                                    egui::TextEdit::singleline(&mut self.inline_label_buf)
                                        .desired_width(56.0)
                                        .font(FontId::proportional(12.0)),
                                );
                                if self.inline_just_opened {
                                    te.request_focus();
                                    self.inline_just_opened = false;
                                }
                                if ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                                    commit = true;
                                } else if ui.input(|i| i.key_pressed(egui::Key::Escape)) {
                                    cancel = true;
                                }
                            });
                        if commit {
                            self.commit_inline_atom_edit();
                        } else if cancel {
                            self.inline_atom_edit = None;
                        }
                    }
                }
            });

        // ---- 6. Authentic Floating ChiralCraft Palettes (`Analysis`, `Properties`, `Periodic Table`, `Info`, `Character Map`) ----
        let active_mol = self
            .tab()
            .selected_mols
            .first()
            .and_then(|&mi| self.doc().scene.molecules.get(mi))
            .or_else(|| self.doc().scene.molecules.first())
            .cloned();

        // Poll background PubChem Structure->Name lookups and upgrade any stamped canvas captions in place!
        let resolved_events = self.live_name_resolver.poll();
        if self.live_name_resolver.has_any_pending() {
            ctx.request_repaint_after(std::time::Duration::from_millis(80));
        }
        for (cap_id_opt, _smi, resolved) in resolved_events {
            if let Some(disp) = resolved.format_display() {
                if let Some(cid) = cap_id_opt {
                    let target_tab_id = self.pending_pubchem_caps.remove(&cid);
                    let target_tab = if let Some(tid) = target_tab_id {
                        self.tabs.iter_mut().find(|t| t.tab_id == tid)
                    } else {
                        self.tabs.get_mut(self.active_tab)
                    };
                    // Record as an undo step: mutating the scene in place left the tab clean and
                    // let Undo/Redo resurrect the stale offline name.
                    if let Some(tab) = target_tab {
                        if tab.doc.scene.captions.iter().any(|c| c.id == cid && c.text != disp) {
                            tab.doc.apply_scene_tx("PubChem Name", |s, _| {
                                for cap in &mut s.captions {
                                    if cap.id == cid {
                                        cap.text = disp.clone();
                                    }
                                }
                            });
                        }
                    }
                }
                self.status_msg = format!("PubChem IUPAC Name: {disp}");
            } else if let Some(cid) = cap_id_opt {
                self.pending_pubchem_caps.remove(&cid);
                self.status_msg = "PubChem lookup found no exact match (kept offline name)".to_string();
            }
        }

        if self.show_analysis {
            let live_info = active_mol.as_ref().map(|m| {
                let (nm, smi, _) = self.compute_live_name(m);
                let pend = self.live_name_resolver.is_pending(&smi);
                (nm, smi, pend)
            });
            let live_nm_ref = live_info.as_ref().map(|(nm, _, _)| nm.as_str());
            let is_pend = live_info.as_ref().map(|(_, _, p)| *p).unwrap_or(false);
            match show_analysis_window(
                ctx,
                &mut self.show_analysis,
                &mut self.analysis_state,
                active_mol.as_ref(),
                live_nm_ref,
                is_pend,
            ) {
                AnalysisAction::None => {}
                AnalysisAction::Paste(txt) => self.paste_caption_tx(txt),
                AnalysisAction::PasteNameCentered => self.do_structure_to_name_hybrid(ctx, false),
                AnalysisAction::ResolvePubChemName => {
                    if let Some((_, smi, _)) = live_info {
                        if let Some(cached) = self.live_name_resolver.resolve(&smi, None, ctx) {
                            self.status_msg = format!("PubChem IUPAC Name: {cached}");
                        } else {
                            self.status_msg = "Resolving IUPAC name via PubChem...".to_string();
                        }
                    }
                }
                AnalysisAction::LinkCaption(kind) => self.attach_linked_caption(kind),
            }
        }
        if self.show_properties {
            if let Some(txt) = show_properties_window(ctx, &mut self.show_properties, &mut self.props_state, active_mol.as_ref()) {
                self.paste_caption_tx(txt);
            }
        }
        if self.show_periodic_table {
            if let Some(z) = show_periodic_table_window(ctx, &mut self.show_periodic_table, &mut self.active_element_z) {
                self.status_msg = format!("Selected element {} (click any atom on canvas)", get_element(z).map(|e| e.symbol).unwrap_or("?"));
            }
        }
        if self.show_char_map {
            if let Some(sym) = show_character_map_window(ctx, &mut self.show_char_map) {
                self.paste_caption_tx(sym);
            }
        }
        if self.show_info_window {
            let dd = match (self.drag_start_canvas, self.drag_cur_canvas) {
                (Some((sx, sy)), Some((ex, ey))) => Some((ex - sx, ey - sy)),
                _ => None,
            };
            let sel_bb = self.selection_bbox();
            show_info_window(ctx, &mut self.show_info_window, self.cursor_doc_pos, dd, sel_bb, self.zoom);
        }

        // ---- 6b. ChiralCraft Floating Windows (`NMR`, `HELM Monomers`, `BioGraphics`, `Ring Fill`, `Stoichiometry`) ----
        let mut pro_actions = Vec::new();
        if self.show_nmr_window {
            let pred = self.nmr_prediction.clone().unwrap_or_else(|| {
                if let Some(ref m) = active_mol {
                    predict_1h_nmr(m)
                } else {
                    NmrPrediction {
                        nucleus: "1H",
                        freq_mhz: 400.0,
                        peaks: Vec::new(),
                        report_text: "No structure selected (showing experimental JCAMP-DX overlay if loaded).".to_string(),
                    }
                }
            });
            pro_actions.push(show_nmr_window(
                ctx,
                &mut self.show_nmr_window,
                &pred,
                &mut self.nmr_selected_atom,
                self.jcamp_overlay.as_ref(),
                &mut self.pro_suite.jcamp_path_buf,
            ));
        }
        if self.show_helm_toolbar {
            pro_actions.push(show_helm_monomer_window(
                ctx,
                &mut self.show_helm_toolbar,
                &self.monomer_lib,
                &mut self.helm_seq_buf,
                &mut self.helm_polymer_kind,
            ));
        }
        if self.show_biographics_toolbar {
            pro_actions.push(show_biographics_toolbar(ctx, &mut self.show_biographics_toolbar));
        }
        if self.show_ring_fill_toolbar {
            pro_actions.push(show_ring_fill_toolbar(ctx, &mut self.show_ring_fill_toolbar));
        }
        if self.show_stoichiometry {
            pro_actions.push(show_stoichiometry_window(
                ctx,
                &mut self.show_stoichiometry,
                &mut self.stoich_base_mmol,
                &mut self.stoich_rows,
            ));
        }
        for act in pro_actions {
            self.handle_pro_panel_action(ctx, act);
        }
        if self.show_object_inspector {
            let mut open = self.show_object_inspector;
            let mut rear_op = self.rear_opacity;
            let mut cross_gap = self.crossing_clearance_pt;
            let mut exp_opts = self.export_opts;
            let mut tables_edit = self.doc().scene.tables.clone();
            let mut tlc_edit = self.doc().scene.tlc_plates.clone();
            let changed = show_object_inspector_window(
                ctx,
                &mut open,
                &mut tables_edit,
                &mut tlc_edit,
                &mut rear_op,
                &mut cross_gap,
                &mut exp_opts,
            );
            self.show_object_inspector = open;
            self.rear_opacity = rear_op;
            self.crossing_clearance_pt = cross_gap;
            self.export_opts = exp_opts;
            if tables_edit != self.doc().scene.tables || tlc_edit != self.doc().scene.tlc_plates {
                self.doc_mut().apply_scene_tx("Object Inspector Edit", |s, _| {
                    s.tables = tables_edit;
                    s.tlc_plates = tlc_edit;
                });
            }
            if changed {
                self.persist_user_config();
                self.status_msg = "Updated Object Inspector settings".to_string();
            }
        }
        if self.viewer_3d.open {
            let active_clone = self.active_mol().cloned();
            match show_3d_viewer_window(ctx, &mut self.viewer_3d, active_clone.as_ref()) {
                Viewer3dAction::None => {}
                Viewer3dAction::SnapshotTo2d(m) => {
                    self.doc_mut().stamp_molecule_tx(m, "3D Snapshot to Canvas");
                    self.status_msg = "Stamped rotated 3D projection onto 2D canvas (with depth cueing)".to_string();
                }
                Viewer3dAction::LiftActiveTo3d => {
                    if let Some(mi) = self.active_mol_index() {
                        self.doc_mut().apply_scene_tx("Lift 2D to 3D", |s, _| {
                            if let Some(m) = s.molecules.get_mut(mi) {
                                *m = lift_molecule_2d_to_3d(m);
                            }
                        });
                        self.status_msg = "Lifted active structure to 3D conformer".to_string();
                    }
                }
                Viewer3dAction::ApplyOptimized3d(opt_mol, summary) => {
                    if let Some(mi) = self.active_mol_index() {
                        self.doc_mut().apply_scene_tx("UFF 3D Optimize", |s, _| {
                            if let Some(slot) = s.molecules.get_mut(mi) {
                                *slot = opt_mol;
                            }
                        });
                        self.status_msg = summary;
                    }
                }
            }
        }
        if self.show_diff_window {
            let tab_names: Vec<String> = self.tabs.iter().map(|t| t.doc.name.clone()).collect();
            let diff_act = show_diff_compare_window(
                ctx,
                &mut self.show_diff_window,
                &tab_names,
                &mut self.diff_left_tab,
                &mut self.diff_right_tab,
                &mut self.diff_file_input,
                self.last_diff.as_ref(),
            );
            match diff_act {
                DiffCompareAction::None => {}
                DiffCompareAction::CompareTabs(l, r) => {
                    if let (Some(lt), Some(rt)) = (self.tabs.get(l), self.tabs.get(r)) {
                        self.last_diff = Some(compare_documents(&lt.doc, &rt.doc));
                    }
                }
                DiffCompareAction::CompareActiveWithFile(path) => {
                    match load_document_any(&path) {
                        Ok(other) => {
                            self.last_diff = Some(compare_documents(self.doc(), &other));
                            self.status_msg = format!("Compared active document with {path}");
                        }
                        Err(e) => {
                            self.status_msg = format!("Diff load error: {e}");
                        }
                    }
                }
            }
        }
        if self.cmd_palette.open {
            let pact = show_command_palette(ctx, &mut self.cmd_palette);
            if pact != PaletteCommandAction::None {
                self.execute_palette_command(ctx, pact);
            }
        }

        // ---- 6d. Flagship MedChem, HRMS, Isomers/Resonance, Markush SAR, Retrosynthesis & TOC Suite ----
        let active_mol_clone = self.active_mol().cloned();
        let active_blen = self.doc().settings.bond_length_pt;
        show_medchem_window(ctx, &mut self.pro_suite, active_mol_clone.as_ref());
        show_hrms_window(ctx, &mut self.pro_suite, active_mol_clone.as_ref());
        let suite_acts = [
            show_isomers_window(ctx, &mut self.pro_suite, active_mol_clone.as_ref()),
            show_markush_window(ctx, &mut self.pro_suite, active_mol_clone.as_ref(), active_blen),
            show_retrosynthesis_window(ctx, &mut self.pro_suite, active_mol_clone.as_ref(), active_blen),
            show_publication_composer_window(ctx, &mut self.pro_suite, &self.tabs[self.active_tab].doc.scene),
        ];
        for sa in suite_acts {
            self.handle_pro_suite_action(sa);
        }

        // ---- 6e. All-Discipline Omni-Suite Windows (IR/UV/2D-NMR, 18e-/Polymer, Energy Profiles, ELN/GHS, LaTeX, Screening) ----
        if self.omni_suite.show_spectroscopy
            || self.quantum.show_window
            || (self.show_right_dock && self.dock_tab == DockTab::Quantum)
        {
            self.qc_runner.ensure_detected(
                &self.user_config.xtb_path,
                &self.user_config.orca_path,
                &self.user_config.gaussian_path,
            );
        }
        if let Some(res) = self.qc_runner.poll() {
            match res {
                Ok(out) => {
                    let n_modes = out.ir_prediction.bands.len();
                    let src = out.ir_prediction.source.label();
                    self.omni_suite.ir_render.scale_factor = out.ir_prediction.default_scale_factor;
                    self.omni_suite.ir_render.use_per_band_widths = false;
                    self.omni_suite.ir_engine = IrEngineMode::QuantumImported;
                    self.omni_suite.ir_selected_band = None;
                    self.omni_suite.ir_quantum_mol_sig = out
                        .cache_key
                        .split('|')
                        .next()
                        .unwrap_or("")
                        .to_string();
                    self.omni_suite.ir_status_msg = format!(
                        "{} ({src}) ready in {:.1}s ({n_modes} modes)",
                        out.backend.label(),
                        out.elapsed_secs
                    );
                    self.omni_suite.ir_quantum_pred = Some(out.ir_prediction);
                    self.quantum.parsed_report = Some(out.report);
                    self.status_msg = self.omni_suite.ir_status_msg.clone();
                }
                Err(e) => {
                    self.status_msg = format!("Quantum FTIR failed: {e}");
                }
            }
        }
        if let Some(res) = self.nist_ir_fetch.poll() {
            match res {
                Ok(sp) => {
                    self.status_msg = format!("Loaded NIST IR overlay: {}", sp.title);
                    self.omni_suite.ir_exp_overlay = Some(sp);
                }
                Err(e) => {
                    self.status_msg = format!("NIST IR fetch failed: {e}");
                }
            }
        }
        if self.qc_runner.running || self.nist_ir_fetch.fetching {
            ctx.request_repaint_after(std::time::Duration::from_millis(120));
        }

        let active_name_str = active_mol_clone
            .as_ref()
            .map(|m| self.compute_live_name(m).0)
            .unwrap_or_else(|| "Target Compound".to_string());
        let omni_acts = [
            show_spectroscopy_window(
                ctx,
                &mut self.omni_suite,
                active_mol_clone.as_ref(),
                &active_name_str,
                &self.qc_runner,
                &self.nist_ir_fetch,
            ),
            show_ino_poly_window(ctx, &mut self.omni_suite, active_mol_clone.as_ref(), active_blen),
            show_energy_profile_window(ctx, &mut self.omni_suite),
            show_eln_safety_window(ctx, &mut self.omni_suite, active_mol_clone.as_ref(), &active_name_str, &self.stoich_rows),
            show_latex_export_window(ctx, &mut self.omni_suite, &self.tabs[self.active_tab].doc, active_mol_clone.as_ref()),
            show_screening_window(ctx, &mut self.omni_suite, active_blen),
        ];
        for oa in omni_acts {
            self.handle_omni_suite_action(ctx, oa);
        }

        // ---- 6f. Quantum Chemistry Floating Window ----
        let q_win_act = show_quantum_window(
            ctx,
            &mut self.quantum,
            active_mol_clone.as_ref(),
            &mut self.user_config,
            &self.qc_runner.binaries,
        );
        self.handle_dock_quantum_action(ctx, q_win_act);

        // ---- 6c. Live PubChem / ChEMBL Chemical Database Search Window ----
        self.chem_db_state.poll();
        if self.chem_db_state.searching {
            ctx.request_repaint_after(std::time::Duration::from_millis(80));
        }
        if self.show_chem_db {
            let has_mol = self.active_mol().is_some();
            let db_act = show_chem_db_window(ctx, &mut self.show_chem_db, &mut self.chem_db_state, has_mol);
            match db_act {
                ChemDbAction::None => {}
                ChemDbAction::UseCanvasSelection => {
                    let mode = match self.chem_db_state.mode {
                        DbSearchMode::NameOrId => DbSearchMode::Similarity2D,
                        other => other,
                    };
                    self.do_search_canvas_in_db(ctx, mode);
                }
                ChemDbAction::ImportHit { hit, with_caption } => {
                    self.do_import_db_hit(&hit, with_caption);
                }
                ChemDbAction::ImportBatch(hits) => {
                    self.do_import_db_batch(&hits);
                }
                ChemDbAction::CopyText { label, text } => {
                    ctx.copy_text(text.clone());
                    self.status_msg = format!("Copied {label}: {text}");
                }
            }
        }

        // ---- 7. Modal Dialogs ----
        match self.dialog.clone() {
            ActiveDialog::ScaleDialog => {
                egui::Window::new("Scale Objects").collapsible(false).show(ctx, |ui| {
                    ui.horizontal(|ui| { ui.label("Scale (%):"); ui.text_edit_singleline(&mut self.scale_pct_input); });
                    ui.horizontal(|ui| {
                        if ui.button("Scale").clicked() {
                            if let Ok(pct) = self.scale_pct_input.trim().parse::<f64>() {
                                if pct > 5.0 && pct < 2000.0 { self.do_scale_targets(pct / 100.0); }
                            }
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::RotateDialog => {
                egui::Window::new("Rotate Objects").collapsible(false).show(ctx, |ui| {
                    ui.horizontal(|ui| { ui.label("Angle (degrees CW):"); ui.text_edit_singleline(&mut self.rotate_deg_input); });
                    ui.horizontal(|ui| {
                        if ui.button("Rotate").clicked() {
                            if let Ok(deg) = self.rotate_deg_input.trim().parse::<f64>() {
                                self.do_rotate_targets(deg);
                            }
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::NewCaption { x_tenths, y_tenths } => {
                egui::Window::new("Insert Text Caption").collapsible(false).show(ctx, |ui| {
                    ui.text_edit_singleline(&mut self.label_input);
                    ui.horizontal(|ui| {
                        if ui.button("OK").clicked() {
                            let t = self.label_input.clone();
                            self.doc_mut().apply_scene_tx("Caption", |s, nid| {
                                let id = *nid; *nid += 1;
                                s.captions.push(CaptionNode {
                                    id, text: t, x: (x_tenths as f64) / 10.0, y: (y_tenths as f64) / 10.0,
                                    element: None, generic_nickname: None, warning: None, justification: LabelJustification::Left,
                                });
                            });
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::PasteNameAsStructure => {
                egui::Window::new("Convert Name to Structure").collapsible(false).default_width(380.0).show(ctx, |ui| {
                    ui.label("Enter IUPAC or common chemical name (e.g. Aspirin, 4-nitrobenzoic acid, 2-chlorobutane, ethyl propanoate):");
                    let te = ui.text_edit_singleline(&mut self.name_input);
                    let enter = te.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
                    ui.horizontal(|ui| {
                        if ui.button("OK").clicked() || enter {
                            let nm = self.name_input.clone();
                            match self.do_name_to_structure(&nm) {
                                Ok(()) => self.dialog = ActiveDialog::None,
                                Err(e) => self.status_msg = e,
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::PasteBiopolymer(kind) => {
                let title = match kind {
                    PolymerKind::Peptide => "Paste Special — FASTA Peptide / Biopolymer",
                    PolymerKind::Dna => "Paste Special — FASTA DNA",
                    PolymerKind::Rna => "Paste Special — FASTA RNA",
                };
                egui::Window::new(title).collapsible(false).default_width(380.0).show(ctx, |ui| {
                    ui.label("Enter 1-letter sequence or FASTA block (e.g. >Seq1\\nACDEFGHIK):");
                    ui.text_edit_multiline(&mut self.helm_seq_buf);
                    ui.horizontal(|ui| {
                        if ui.button("Paste").clicked() {
                            let seq = self.helm_seq_buf.clone();
                            self.handle_pro_panel_action(ctx, ProPanelAction::InsertBiopolymer { seq, kind });
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::PasteHelm => {
                egui::Window::new("Paste Special — HELM").collapsible(false).default_width(400.0).show(ctx, |ui| {
                    ui.label("Enter HELM v1/v2 notation (e.g. PEPTIDE1{A.C.D.E.F}$$$$V2.0):");
                    ui.text_edit_singleline(&mut self.helm_input);
                    ui.horizontal(|ui| {
                        if ui.button("Paste").clicked() {
                            let blen = self.doc().settings.bond_length_pt;
                            match helm_to_molecule(&self.helm_input, blen) {
                                Ok(m) => {
                                    let cnt = m.atoms.len();
                                    self.doc_mut().stamp_molecule_tx(m, "Paste HELM");
                                    self.status_msg = format!("Pasted {cnt}-monomer HELM biopolymer");
                                    self.dialog = ActiveDialog::None;
                                }
                                Err(e) => self.status_msg = format!("HELM parse error: {e}"),
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::AddToNameDict => {
                egui::Window::new("Add Structure to Name Dictionary").collapsible(false).show(ctx, |ui| {
                    ui.horizontal(|ui| { ui.label("Name:"); ui.text_edit_singleline(&mut self.dict_name_input); });
                    ui.horizontal(|ui| {
                        if ui.button("Add").clicked() {
                            if let Some(m) = self.active_mol() {
                                let smi = molecule_to_smiles(m);
                                let nm = self.dict_name_input.trim().to_string();
                                if !nm.is_empty() {
                                    self.name_dict.add_entry(&nm, &smi);
                                    self.status_msg = format!("Added \"{nm}\" ({smi}) to Structure↔Name dictionary");
                                }
                            }
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::PasteSmiles => {
                egui::Window::new("Paste Special — SMILES").collapsible(false).show(ctx, |ui| {
                    ui.text_edit_singleline(&mut self.smiles_input);
                    ui.horizontal(|ui| {
                        if ui.button("Paste").clicked() {
                            match smiles_to_molecule(self.smiles_input.trim()) {
                                Ok(mut m) => {
                                    let blen = self.doc().settings.bond_length_pt;
                                    clean_up_structure_preserving_stereo(&mut m, &[], blen);
                                    for a in &mut m.atoms { a.x += 200.0; a.y += 180.0; }
                                    self.doc_mut().stamp_molecule_tx(m, "Paste SMILES");
                                    self.dialog = ActiveDialog::None;
                                }
                                Err(e) => self.status_msg = format!("SMILES error: {e}"),
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::PasteMolText => {
                egui::Window::new("Paste Special — MOL / SDF Text").collapsible(false).show(ctx, |ui| {
                    egui::ScrollArea::vertical().max_height(200.0).show(ui, |ui| {
                        ui.text_edit_multiline(&mut self.multiline_input);
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Paste").clicked() {
                            let txt = self.multiline_input.clone();
                            let mols_res = if txt.contains("$$$$") {
                                sdf_to_molecules(&txt).map(|v| v.into_iter().map(|r| r.molecule).collect())
                            } else if txt.contains("V3000") {
                                molfile_v3000_to_molecule(&txt).map(|m| vec![m])
                            } else {
                                molfile_to_molecule(&txt).map(|m| vec![m])
                            };
                            match mols_res {
                                Ok(mut mols) => {
                                    let bl = self.doc().settings.bond_length_pt;
                                    chiralcraft_doc::rescale_angstrom_molecules(&mut mols, bl);
                                    chiralcraft_doc::layout_imported_molecules(&mut mols, 540.0, 719.75);
                                    stamp_molecules_tx(self.doc_mut(), mols, "Paste MOL/SDF");
                                    self.dialog = ActiveDialog::None;
                                }
                                Err(e) => self.status_msg = format!("MOL/SDF parse error: {e}"),
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::PasteEditableSvg => {
                egui::Window::new("Paste Special — Editable ChiralCraft SVG").collapsible(false).show(ctx, |ui| {
                    egui::ScrollArea::vertical().max_height(200.0).show(ui, |ui| {
                        ui.text_edit_multiline(&mut self.multiline_input);
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Paste").clicked() {
                            match InteractiveDocument::from_svg(&self.multiline_input) {
                                Ok(d) => {
                                    let cnt = d.scene.molecules.len();
                                    stamp_molecules_tx(self.doc_mut(), d.scene.molecules, "Paste Editable SVG");
                                    self.status_msg = format!("Imported {cnt} editable structure(s) from SVG");
                                    self.dialog = ActiveDialog::None;
                                }
                                Err(e) => self.status_msg = format!("Editable SVG parse error: {e}"),
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::PasteCdxml => {
                egui::Window::new("Paste Special — CDXML").collapsible(false).show(ctx, |ui| {
                    egui::ScrollArea::vertical().max_height(200.0).show(ui, |ui| {
                        ui.text_edit_multiline(&mut self.multiline_input);
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Paste").clicked() {
                            match CdxmlDocument::parse(&self.multiline_input) {
                                Ok(cdxml) => {
                                    let d = InteractiveDocument::from_cdxml(&cdxml);
                                    stamp_molecules_tx(self.doc_mut(), d.scene.molecules, "Paste CDXML");
                                    self.dialog = ActiveDialog::None;
                                }
                                Err(e) => self.status_msg = format!("CDXML parse error: {e:?}"),
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::OpenFile => {
                egui::Window::new("Open").collapsible(false).default_width(420.0).show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        if ui.button("Up ..").clicked() {
                            if let Some(p) = self.browser_dir.parent() {
                                self.browser_dir = p.to_path_buf();
                            }
                        }
                        ui.label(self.browser_dir.display().to_string());
                    });
                    let entries = list_directory_entries(&self.browser_dir);
                    egui::ScrollArea::vertical().max_height(180.0).show(ui, |ui| {
                        for (name, is_dir) in entries {
                            let disp = if is_dir { format!("[DIR] {name}") } else { name.clone() };
                            if ui.selectable_label(false, disp).clicked() {
                                let next = self.browser_dir.join(&name);
                                if is_dir {
                                    self.browser_dir = next;
                                } else {
                                    self.file_path_buf = next.display().to_string();
                                }
                            }
                        }
                    });
                    ui.separator();
                    ui.label("File path (.cdxml, .cdx, .mol, .sdf, .smi, .xyz, .cif, .pdb, POSCAR, .jdx, .out/.log, .svg, .png, .pdf):");
                    ui.text_edit_singleline(&mut self.file_path_buf);
                    ui.horizontal(|ui| {
                        if ui.button("Browse OS...").clicked() {
                            if let Some(p) = native_dialogs::pick_open_file(
                                "Open Chemical Document",
                                &["cdxml", "cdx", "mol", "sdf", "rxn", "smi", "xyz", "cif", "pdb", "jdx", "out", "log", "molden", "svg", "png", "pdf"],
                            ) {
                                self.file_path_buf = p.clone();
                                self.open_path_into_new_tab(&p);
                                if !self.status_msg.starts_with("Open failed") && !self.status_msg.starts_with("JCAMP-DX open error") {
                                    self.dialog = ActiveDialog::None;
                                }
                            }
                        }
                        if ui.button("Open").clicked() {
                            let path = self.file_path_buf.trim().to_string();
                            if !path.is_empty() {
                                self.open_path_into_new_tab(&path);
                                if !self.status_msg.starts_with("Open failed") && !self.status_msg.starts_with("JCAMP-DX open error") {
                                    self.dialog = ActiveDialog::None;
                                }
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::EditCaption { cap_idx } => {
                egui::Window::new("Edit Text Caption").collapsible(false).show(ctx, |ui| {
                    ui.text_edit_singleline(&mut self.label_input);
                    ui.horizontal(|ui| {
                        if ui.button("OK").clicked() {
                            let t = self.label_input.clone();
                            self.doc_mut().apply_scene_tx("Edit Caption", |s, _| {
                                if let Some(c) = s.captions.get_mut(cap_idx) { c.text = t; }
                            });
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::PasteInchi => {
                egui::Window::new("Paste Special — InChI").collapsible(false).show(ctx, |ui| {
                    ui.text_edit_singleline(&mut self.inchi_input);
                    ui.horizontal(|ui| {
                        if ui.button("Paste").clicked() {
                            match inchi_to_molecule(self.inchi_input.trim()) {
                                Ok(mut m) => {
                                    let blen = self.doc().settings.bond_length_pt;
                                    clean_up_structure_preserving_stereo(&mut m, &[], blen);
                                    for a in &mut m.atoms { a.x += 200.0; a.y += 180.0; }
                                    self.doc_mut().stamp_molecule_tx(m, "Paste InChI");
                                    self.dialog = ActiveDialog::None;
                                }
                                Err(e) => self.status_msg = format!("InChI error: {e}"),
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::InsertFile => {
                egui::Window::new("Insert File into Active Document").collapsible(false).default_width(420.0).show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        if ui.button("Up ..").clicked() {
                            if let Some(p) = self.browser_dir.parent() {
                                self.browser_dir = p.to_path_buf();
                            }
                        }
                        ui.label(self.browser_dir.display().to_string());
                    });
                    let entries = list_directory_entries(&self.browser_dir);
                    egui::ScrollArea::vertical().max_height(160.0).show(ui, |ui| {
                        for (name, is_dir) in entries {
                            let disp = if is_dir { format!("[DIR] {name}") } else { name.clone() };
                            if ui.selectable_label(false, disp).clicked() {
                                let next = self.browser_dir.join(&name);
                                if is_dir { self.browser_dir = next; } else { self.file_path_buf = next.display().to_string(); }
                            }
                        }
                    });
                    ui.separator();
                    ui.label("File path (.cdxml, .cdx, .mol, .sdf, .rxn, .ct, .cml, .tgf, .skc):");
                    ui.text_edit_singleline(&mut self.file_path_buf);
                    ui.horizontal(|ui| {
                        if ui.button("Browse OS...").clicked() {
                            if let Some(p) = native_dialogs::pick_open_file(
                                "Insert Chemical File",
                                &["cdxml", "cdx", "mol", "sdf", "rxn", "ct", "cml", "tgf", "skc", "smi"],
                            ) {
                                self.file_path_buf = p;
                            }
                        }
                        if ui.button("Insert").clicked() {
                            match load_document_any(&self.file_path_buf) {
                                Ok(d) => {
                                    let cnt = d.scene.molecules.len();
                                    stamp_molecules_tx(self.doc_mut(), d.scene.molecules, "Insert File");
                                    self.status_msg = format!("Inserted {cnt} structure(s) from {}", self.file_path_buf);
                                    self.dialog = ActiveDialog::None;
                                }
                                Err(e) => self.status_msg = format!("Insert failed: {e}"),
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::AtomProperties { mol_idx, atom_idx } => {
                egui::Window::new("Atom Properties").collapsible(false).show(ctx, |ui| {
                    let orig = self
                        .doc()
                        .scene
                        .molecules
                        .get(mol_idx)
                        .and_then(|m| m.atoms.get(atom_idx))
                        .cloned();
                    if let Some(mut a) = orig.clone() {
                        ui.add(egui::DragValue::new(&mut a.z).range(1..=118).prefix("Element Z: "));
                        let sym = get_element(a.z).map(|e| format!("{} ({})", e.symbol, e.name)).unwrap_or_default();
                        ui.label(sym);
                        ui.add(egui::DragValue::new(&mut a.charge).range(-8..=8).prefix("Formal Charge: "));
                        ui.add(egui::DragValue::new(&mut a.isotope).range(0..=300).prefix("Isotope Mass (0=default): "));
                        egui::ComboBox::from_label("Radical")
                            .selected_text(match a.radical { 1 => "Singlet (:)", 2 => "Doublet (•)", 3 => "Triplet", _ => "None" })
                            .show_ui(ui, |ui| {
                                ui.selectable_value(&mut a.radical, 0, "None");
                                ui.selectable_value(&mut a.radical, 1, "Singlet (:)");
                                ui.selectable_value(&mut a.radical, 2, "Doublet (•)");
                                ui.selectable_value(&mut a.radical, 3, "Triplet");
                            });
                        ui.checkbox(&mut a.ignore_warnings, "Ignore Chemical Warnings");
                        ui.checkbox(&mut a.is_query, "Query / Generic Atom");
                        if Some(&a) != orig.as_ref() {
                            let key = egui::Id::new(("cc_atom_props", self.tab().tab_id, mol_idx, atom_idx));
                            let (gesture, mut last) = (undo_gesture_token(ctx), self.undo_coalesce);
                            coalesced_scene_tx(self.doc_mut(), gesture, &mut last, key, "Atom Properties", |s| {
                                if let Some(slot) = s.molecules.get_mut(mol_idx).and_then(|m| m.atoms.get_mut(atom_idx)) {
                                    *slot = a;
                                }
                            });
                            self.undo_coalesce = last;
                        }
                    } else {
                        ui.label("No atom selected.");
                    }
                    if ui.button("OK").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::BracketProperties => {
                egui::Window::new("Bracket Properties").collapsible(false).show(ctx, |ui| {
                    let mut usage = self.tab().bracket_usage.clone();
                    let mut pat = self.tab().bracket_repeat_pattern.clone();
                    let mut cnt = self.tab().bracket_repeat_count.clone();
                    egui::ComboBox::from_label("Bracket Usage").selected_text(&usage).show_ui(ui, |ui| {
                        for u in ["SRU", "Component", "Mixture", "Racemic", "Anomeric", "Multiple Group"] {
                            ui.selectable_value(&mut usage, u.to_string(), u);
                        }
                    });
                    egui::ComboBox::from_label("Repeat Pattern").selected_text(&pat).show_ui(ui, |ui| {
                        for p in ["Head-to-Tail", "Head-to-Head", "Either/Unknown"] {
                            ui.selectable_value(&mut pat, p.to_string(), p);
                        }
                    });
                    ui.horizontal(|ui| { ui.label("Repeat Count:"); ui.text_edit_singleline(&mut cnt); });
                    self.tab_mut().bracket_usage = usage;
                    self.tab_mut().bracket_repeat_pattern = pat;
                    self.tab_mut().bracket_repeat_count = cnt;
                    if ui.button("OK").clicked() {
                        self.status_msg = format!("Bracket: {} ({}, {})", self.tab().bracket_usage, self.tab().bracket_repeat_pattern, self.tab().bracket_repeat_count);
                        self.dialog = ActiveDialog::None;
                    }
                });
            }
            ActiveDialog::CheckStructureDialog => {
                let warns = self.doc().check_structure_with_preferences();
                egui::Window::new("Check Structure").collapsible(false).default_width(360.0).show(ctx, |ui| {
                    if warns.is_empty() {
                        ui.label("No chemical errors were found.");
                        if ui.button("OK").clicked() { self.dialog = ActiveDialog::None; }
                    } else {
                        let idx = self.check_warn_idx.min(warns.len() - 1);
                        let w = &warns[idx];
                        ui.strong(format!("Warning {} of {} (Object ID {}):", idx + 1, warns.len(), w.object_id));
                        ui.colored_label(Color32::from_rgb(190, 20, 20), &w.message);
                        ui.separator();
                        ui.horizontal(|ui| {
                            if ui.button("Ignore").clicked() {
                                let oid = w.object_id;
                                self.doc_mut().apply_scene_tx("Ignore Warning", |s, _| {
                                    for m in &mut s.molecules {
                                        for a in &mut m.atoms { if a.id == oid { a.ignore_warnings = true; } }
                                        for b in &mut m.bonds { if b.id == oid { b.ignore_warnings = true; } }
                                    }
                                });
                            }
                            if ui.button("Next").clicked() && idx + 1 < warns.len() {
                                self.check_warn_idx = idx + 1;
                            }
                            if ui.button("Stop").clicked() {
                                self.dialog = ActiveDialog::None;
                            }
                        });
                    }
                });
            }
            ActiveDialog::ColorPicker => {
                egui::Window::new("Colors — Custom Foreground").collapsible(false).show(ctx, |ui| {
                    ui.color_edit_button_srgba(&mut self.active_color);
                    ui.checkbox(&mut self.color_by_element, "Color By Element (CPK)");
                    if ui.button("OK").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::SaveAs => {
                egui::Window::new("Save As").collapsible(false).default_width(460.0).show(ctx, |ui| {
                    ui.label("Format (.cdxml, .cdx, .mol, .sdf, .rxn, .smi, .xyz, .svg, .png, .pdf, .eps, .emf, .cml, .ct, .tgf, .skc, .sln, .inchi):");
                    ui.text_edit_singleline(&mut self.file_path_buf);
                    ui.horizontal_wrapped(|ui| {
                        for ext in [".cdxml", ".cdx", ".mol", ".sdf", ".rxn", ".smi", ".xyz", ".svg", ".png", ".pdf", ".eps", ".emf", ".cml", ".ct", ".tgf", ".skc", ".sln", ".inchi", ".cds"] {
                            if ui.small_button(ext).clicked() {
                                let p = std::path::Path::new(&self.file_path_buf);
                                let ext_no_dot = ext.trim_start_matches('.');
                                let updated = if self.file_path_buf.trim().is_empty() {
                                    format!("export.{ext_no_dot}")
                                } else {
                                    p.with_extension(ext_no_dot).to_string_lossy().to_string()
                                };
                                self.file_path_buf = updated;
                            }
                        }
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Browse OS...").clicked() {
                            let suggested = if self.file_path_buf.trim().is_empty() {
                                "untitled.cdxml"
                            } else {
                                self.file_path_buf.trim()
                            };
                            if let Some(p) = native_dialogs::pick_save_file("Save Chemical Document As", suggested) {
                                self.file_path_buf = p;
                            }
                        }
                        // Typed paths bypass the native dialog's overwrite prompt: require a second click.
                        let overwrite_id = egui::Id::new("cc_save_as_overwrite_confirm");
                        let typed = self.file_path_buf.clone();
                        let confirmed: Option<String> = ctx.data_mut(|d| d.get_temp(overwrite_id));
                        let is_confirmed = confirmed.as_deref() == Some(typed.as_str());
                        let needs_confirm = !typed.trim().is_empty()
                            && std::path::Path::new(&typed).exists()
                            && self.tab().file_path.as_deref() != Some(typed.as_str())
                            && !is_confirmed;
                        let save_clicked = ui.button(if is_confirmed { "Overwrite" } else { "Save" }).clicked();
                        if save_clicked && needs_confirm {
                            ctx.data_mut(|d| d.insert_temp(overwrite_id, typed.clone()));
                            self.status_msg = format!("{typed} already exists — click Overwrite to replace it");
                        } else if save_clicked {
                            ctx.data_mut(|d| d.remove::<String>(overwrite_id));
                            let dark_ink = self.ui_theme == UiTheme::DarkSlide;
                            match save_document_any_styled(
                                self.doc(),
                                &self.file_path_buf,
                                self.export_opts,
                                dark_ink,
                                self.rear_opacity,
                                self.crossing_clearance_pt,
                            ) {
                                Ok(()) if save_as_format_is_lossy(self.doc(), &self.file_path_buf) => {
                                    // Export only: keep the tab's real file and dirty state (see save_as_format_is_lossy).
                                    self.status_msg = format!(
                                        "Exported {} (format cannot hold the whole drawing; document still unsaved)",
                                        self.file_path_buf
                                    );
                                    self.dialog = ActiveDialog::None;
                                }
                                Ok(()) => {
                                    let path = self.file_path_buf.clone();
                                    if self.roundtrip_svg_path.as_deref() != Some(path.as_str()) {
                                        self.roundtrip_svg_path = None;
                                    }
                                    self.tab_mut().file_path = Some(path.clone());
                                    if let Some(fname) = std::path::Path::new(&path).file_name().and_then(|s| s.to_str()) {
                                        self.doc_mut().name = fname.to_string();
                                    }
                                    self.doc_mut().save();
                                    if !cfg!(test) {
                                        self.user_config.push_recent(&path);
                                    }
                                    self.status_msg = format!("Saved {path}");
                                    self.dialog = ActiveDialog::None;
                                }
                                Err(e) => self.status_msg = format!("Save failed: {e}"),
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::PageSetup => {
                egui::Window::new("Page Setup").collapsible(false).default_width(340.0).show(ctx, |ui| {
                    let mut paper = self.tab().paper_size.clone();
                    egui::ComboBox::from_label("Paper Size").selected_text(&paper).show_ui(ui, |ui| {
                        for (name, w, h) in [
                            ("US Letter", 540.0_f64, 719.75_f64),
                            ("A4", 523.3_f64, 769.9_f64),
                            ("US Legal", 540.0_f64, 936.0_f64),
                            ("Tabloid", 720.0_f64, 1152.0_f64),
                        ] {
                            if ui.selectable_value(&mut paper, name.to_string(), name).clicked() {
                                self.tab_mut().page_width_pt = w;
                                self.tab_mut().page_height_pt = h;
                            }
                        }
                    });
                    self.tab_mut().paper_size = paper;
                    ui.horizontal(|ui| {
                        ui.label("Orientation:");
                        ui.radio_value(&mut self.tab_mut().landscape, false, "Portrait");
                        ui.radio_value(&mut self.tab_mut().landscape, true, "Landscape");
                    });
                    let mut wp = self.doc().scene.header_footer.width_pages;
                    let mut hp = self.doc().scene.header_footer.height_pages;
                    ui.horizontal(|ui| {
                        ui.add(egui::DragValue::new(&mut wp).range(1..=4).prefix("Pages Wide: "));
                        ui.add(egui::DragValue::new(&mut hp).range(1..=4).prefix("Pages Tall: "));
                    });
                    if wp != self.doc().scene.header_footer.width_pages || hp != self.doc().scene.header_footer.height_pages {
                        self.doc_mut().apply_scene_tx("Page Setup", |s, _| {
                            s.header_footer.width_pages = wp;
                            s.header_footer.height_pages = hp;
                        });
                    }
                    let (ew, eh) = self.tab().effective_page_size_pt();
                    ui.label(format!("Effective page: {:.1} × {:.1} pt", ew, eh));
                    if ui.button("OK").clicked() {
                        self.status_msg = format!(
                            "Page Setup: {} ({})",
                            self.tab().paper_size,
                            if self.tab().landscape { "Landscape" } else { "Portrait" }
                        );
                        self.dialog = ActiveDialog::None;
                    }
                });
            }
            ActiveDialog::DocumentProperties => {
                egui::Window::new("Document Properties").collapsible(false).default_width(380.0).show(ctx, |ui| {
                    let n_mols = self.doc().scene.molecules.len();
                    let n_atoms: usize = self.doc().scene.molecules.iter().map(|m| m.atoms.len()).sum();
                    let n_bonds: usize = self.doc().scene.molecules.iter().map(|m| m.bonds.len()).sum();
                    let n_objs = self.doc().scene.arrows.len()
                        + self.doc().scene.curves.len()
                        + self.doc().scene.graphics.len()
                        + self.doc().scene.captions.len()
                        + self.doc().scene.tables.len()
                        + self.doc().scene.tlc_plates.len();
                    let path_lbl = self.tab().file_path.clone().unwrap_or_else(|| "(Unsaved document)".to_string());
                    ui.label(format!("File: {path_lbl}"));
                    ui.label(format!("Stylesheet: {} ({:.2} pt)", self.doc().settings.stylesheet_name, self.doc().settings.bond_length_pt));
                    ui.label(format!("Contents: {n_mols} molecule(s), {n_atoms} atom(s), {n_bonds} bond(s), {n_objs} graphic/text object(s)"));
                    ui.separator();
                    ui.horizontal(|ui| { ui.label("Title:   "); ui.text_edit_singleline(&mut self.tab_mut().meta_title); });
                    ui.horizontal(|ui| { ui.label("Author:  "); ui.text_edit_singleline(&mut self.tab_mut().meta_author); });
                    ui.horizontal(|ui| { ui.label("Subject: "); ui.text_edit_singleline(&mut self.tab_mut().meta_subject); });
                    ui.label("Comments:");
                    ui.text_edit_multiline(&mut self.tab_mut().meta_comments);
                    if ui.button("OK").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::ObjectSettings => {
                egui::Window::new("Object Settings").collapsible(false).default_width(340.0).show(ctx, |ui| {
                    let settings_before = self.doc().settings.clone();
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.bond_length_pt).speed(0.2).range(4.0..=72.0).prefix("Fixed Length (pt): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.bond_spacing_pct).speed(0.5).range(5.0..=50.0).prefix("Bond Spacing (%): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.chain_angle_deg).speed(1.0).range(30.0..=180.0).prefix("Chain Angle (°): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.line_width_pt).speed(0.05).range(0.1..=6.0).prefix("Line Width (pt): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.bold_width_pt).speed(0.05).range(0.5..=10.0).prefix("Bold Width (pt): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.hash_spacing_pt).speed(0.05).range(0.5..=10.0).prefix("Hash Spacing (pt): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.margin_width_pt).speed(0.05).range(0.2..=10.0).prefix("Margin Width (pt): "));
                    if self.doc().settings != settings_before {
                        let key = egui::Id::new(("cc_object_settings", self.tab().tab_id));
                        let (gesture, mut last) = (undo_gesture_token(ctx), self.undo_coalesce);
                        coalesced_settings_changed(self.doc_mut(), gesture, &mut last, key, "Object Settings");
                        self.undo_coalesce = last;
                    }
                    if ui.button("OK").clicked() {
                        self.status_msg = format!("Updated Object Settings (Fixed Length {:.2} pt)", self.doc().settings.bond_length_pt);
                        self.dialog = ActiveDialog::None;
                    }
                });
            }
            ActiveDialog::DocumentSettings => {
                egui::Window::new("Document Settings").collapsible(false).default_width(360.0).show(ctx, |ui| {
                    ui.strong("Drawing");
                    let settings_before = self.doc().settings.clone();
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.bond_length_pt).speed(0.2).range(4.0..=72.0).prefix("Fixed Length (pt): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.bond_spacing_pct).speed(0.5).range(5.0..=50.0).prefix("Bond Spacing (%): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.chain_angle_deg).speed(1.0).range(30.0..=180.0).prefix("Chain Angle (°): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.line_width_pt).speed(0.05).range(0.1..=6.0).prefix("Line Width (pt): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.bold_width_pt).speed(0.05).range(0.5..=10.0).prefix("Bold Width (pt): "));
                    ui.add(egui::DragValue::new(&mut self.doc_mut().settings.hash_spacing_pt).speed(0.05).range(0.5..=10.0).prefix("Hash Spacing (pt): "));
                    let doc_settings_key = egui::Id::new(("cc_document_settings", self.tab().tab_id));
                    if self.doc().settings != settings_before {
                        let (gesture, mut last) = (undo_gesture_token(ctx), self.undo_coalesce);
                        coalesced_settings_changed(self.doc_mut(), gesture, &mut last, doc_settings_key, "Document Settings");
                        self.undo_coalesce = last;
                    }
                    ui.separator();
                    ui.strong("Header / Footer & Layout");
                    let mut wp = self.doc().scene.header_footer.width_pages;
                    let mut hp = self.doc().scene.header_footer.height_pages;
                    ui.horizontal(|ui| {
                        ui.add(egui::DragValue::new(&mut wp).range(1..=4).prefix("Pages Wide: "));
                        ui.add(egui::DragValue::new(&mut hp).range(1..=4).prefix("Pages Tall: "));
                    });
                    let mut hdr = self.doc().scene.header_footer.header.clone().unwrap_or_default();
                    let mut ftr = self.doc().scene.header_footer.footer.clone().unwrap_or_default();
                    ui.horizontal(|ui| { ui.label("Header:"); ui.text_edit_singleline(&mut hdr); });
                    ui.horizontal(|ui| { ui.label("Footer:"); ui.text_edit_singleline(&mut ftr); });
                    let new_hdr = if hdr.is_empty() { None } else { Some(hdr) };
                    let new_ftr = if ftr.is_empty() { None } else { Some(ftr) };
                    if wp != self.doc().scene.header_footer.width_pages
                        || hp != self.doc().scene.header_footer.height_pages
                        || new_hdr != self.doc().scene.header_footer.header
                        || new_ftr != self.doc().scene.header_footer.footer
                    {
                        let (gesture, mut last) = (undo_gesture_token(ctx), self.undo_coalesce);
                        coalesced_scene_tx(self.doc_mut(), gesture, &mut last, doc_settings_key, "Document Settings", |s| {
                            s.header_footer.width_pages = wp;
                            s.header_footer.height_pages = hp;
                            s.header_footer.header = new_hdr;
                            s.header_footer.footer = new_ftr;
                        });
                        self.undo_coalesce = last;
                    }
                    if ui.button("OK").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::Preferences => {
                egui::Window::new("Preferences — Warnings").collapsible(false).show(ctx, |ui| {
                    if ui.checkbox(&mut self.doc_mut().preferences.show_chemical_warnings, "Show Chemical Warnings").changed() {
                        self.show_warnings = self.doc().preferences.show_chemical_warnings;
                    }
                    let names = [
                        "Valence and Charge Errors", "Unbalanced Parentheses", "Invalid Isotopes",
                        "Isolated Bonds", "Atoms Near Other Bonds", "Unspecified Stereochemistry",
                        "Ambiguous Stereochemistry", "Stereobond Between Chiral Atoms", "Linear Atoms",
                        "Enhanced Relative Stereochemistry", "Uninterpretable Labels", "Crossing Bonds",
                        "Radical / Non-standard Valence",
                    ];
                    for (i, n) in names.iter().enumerate() {
                        ui.checkbox(&mut self.doc_mut().preferences.warnings_enabled[i], *n);
                    }
                    if ui.button("Done").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::ListNicknames => {
                let ncount = self.doc().nicknames.len();
                egui::Window::new(format!("List Nicknames ({ncount})")).show(ctx, |ui| {
                    ui.horizontal(|ui| { ui.label("Filter:"); ui.text_edit_singleline(&mut self.nick_search); });
                    let q = self.nick_search.to_lowercase();
                    let names: Vec<String> = self.doc().nicknames.list_names().into_iter().filter(|n| q.is_empty() || n.to_lowercase().contains(&q)).take(60).collect();
                    let mut to_delete = None;
                    egui::ScrollArea::vertical().max_height(260.0).show(ui, |ui| {
                        for n in names {
                            ui.horizontal(|ui| {
                                ui.label(&n);
                                if ui.small_button("Delete").clicked() { to_delete = Some(n); }
                            });
                        }
                    });
                    if let Some(d) = to_delete {
                        self.doc_mut().nicknames.delete_nickname(&d);
                        self.status_msg = format!("Deleted nickname \"{d}\"");
                    }
                    if ui.button("Done").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::DefineNickname => {
                egui::Window::new("Define Nickname").collapsible(false).show(ctx, |ui| {
                    ui.label("Nickname:");
                    ui.text_edit_singleline(&mut self.nick_def_name);
                    if let Some(col) = &self.nick_pending_collision {
                        match col {
                            NicknameCollision::ExistingNickname | NicknameCollision::BothNicknameAndElement(_) if !self.nick_confirm_replace => {
                                ui.colored_label(Color32::from_rgb(180, 90, 0), "That name is already in use by another nickname. Replace existing?");
                                ui.checkbox(&mut self.nick_confirm_replace, "Replace");
                            }
                            NicknameCollision::ElementSymbol(z) | NicknameCollision::BothNicknameAndElement(z) => {
                                ui.colored_label(Color32::RED, format!("Same as element symbol (Z={z}). Override element symbol?"));
                                ui.checkbox(&mut self.nick_confirm_override, "Override");
                            }
                            _ => {}
                        }
                    }
                    ui.horizontal(|ui| {
                        if ui.button("OK").clicked() {
                            if let Some(mol) = self.active_mol().cloned() {
                                let nm = self.nick_def_name.trim().to_string();
                                let rep = self.nick_confirm_replace;
                                let ovr = self.nick_confirm_override;
                                match self.doc_mut().nicknames.define_nickname(&nm, mol, 0, rep, ovr) {
                                    Ok(()) => { self.status_msg = format!("Defined \"{nm}\""); self.dialog = ActiveDialog::None; }
                                    Err(c) => self.nick_pending_collision = Some(c),
                                }
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::AttachData => {
                egui::Window::new("Attach Data").collapsible(false).show(ctx, |ui| {
                    ui.horizontal(|ui| { ui.label("Field Name:"); ui.text_edit_singleline(&mut self.data_key_input); });
                    ui.horizontal(|ui| { ui.label("Data Value:"); ui.text_edit_singleline(&mut self.data_val_input); });
                    ui.horizontal(|ui| {
                        if ui.button("Attach").clicked() {
                            let k = self.data_key_input.clone();
                            let v = self.data_val_input.clone();
                            self.tab_mut().attached_data.push((k.clone(), v.clone()));
                            self.paste_caption_tx(format!("{k}: {v}"));
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::Annotate => {
                egui::Window::new("Annotate").collapsible(false).show(ctx, |ui| {
                    ui.horizontal(|ui| { ui.label("Keyword:"); ui.text_edit_singleline(&mut self.data_key_input); });
                    ui.horizontal(|ui| { ui.label("Content:"); ui.text_edit_singleline(&mut self.data_val_input); });
                    ui.horizontal(|ui| {
                        if ui.button("OK").clicked() {
                            let k = self.data_key_input.clone();
                            let v = self.data_val_input.clone();
                            self.tab_mut().annotations.insert(k.clone(), v.clone());
                            self.status_msg = format!("Annotated {k}=\"{v}\"");
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::FindElement => {
                egui::Window::new("Find Element on Canvas").collapsible(false).show(ctx, |ui| {
                    ui.horizontal(|ui| { ui.label("Element Symbol:"); ui.text_edit_singleline(&mut self.find_elem_input); });
                    ui.horizontal(|ui| {
                        if ui.button("Find").clicked() {
                            let sym = self.find_elem_input.trim();
                            if let Some(el) = get_element_by_symbol(sym) {
                                let matches: Vec<usize> = self
                                    .doc()
                                    .scene
                                    .molecules
                                    .iter()
                                    .enumerate()
                                    .filter(|(_, m)| m.atoms.iter().any(|a| a.z == el.z))
                                    .map(|(mi, _)| mi)
                                    .collect();
                                let cnt = matches.len();
                                self.tab_mut().selected_mols = matches;
                                self.status_msg = format!("Found {} structure(s) containing {}", cnt, el.symbol);
                            }
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::Shortcuts => {
                egui::Window::new("ChiralCraft Shortcuts").collapsible(false).show(ctx, |ui| {
                    ui.label("• Hover Atom + O/N/C/F/S/P/I/H/B(Br)/L(Cl): Change element");
                    ui.label("• Hover Atom + Enter: Open inline label/nickname box");
                    ui.label("• Hover Bond + 1/2/3: Set Single / Double / Triple");
                    ui.label("• 3/4/5/6/7/8 / A: Select 3..8-ring / Benzene tool");
                    ui.label("• Space: Marquee tool | Shift+Ctrl+K: Clean Up Structure");
                    ui.label("• Alt+K: Toggle Aromatic Display | Ctrl+J: Join close atoms");
                    ui.label("• Right-click any flyout tool: Open tear-off sub-palette");
                    if ui.button("Close").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::InsertObject => {
                egui::Window::new("Insert Object").collapsible(false).default_width(340.0).show(ctx, |ui| {
                    ui.label("Object Type:");
                    for (idx, lbl) in [
                        (0, "ChiralCraft Data Table"),
                        (1, "TLC Chromatography Plate"),
                        (2, "Framed Graphic Box (Shadow Rectangle)"),
                        (3, "Text Caption Object"),
                        (4, "Chemical Document File (.cdxml / .cdx / .mol)..."),
                    ] {
                        ui.radio_value(&mut self.insert_obj_kind, idx, lbl);
                    }
                    ui.separator();
                    match self.insert_obj_kind {
                        0 => {
                            ui.horizontal(|ui| {
                                ui.add(egui::DragValue::new(&mut self.insert_obj_rows).range(1..=12).prefix("Rows: "));
                                ui.add(egui::DragValue::new(&mut self.insert_obj_cols).range(1..=12).prefix("Cols: "));
                            });
                        }
                        3 => {
                            ui.horizontal(|ui| {
                                ui.label("Text:");
                                ui.text_edit_singleline(&mut self.insert_obj_text);
                            });
                        }
                        _ => {}
                    }
                    ui.horizontal(|ui| {
                        if ui.button("OK").clicked() {
                            match self.insert_obj_kind {
                                0 => {
                                    let rows = self.insert_obj_rows.max(1);
                                    let cols = self.insert_obj_cols.max(1);
                                    self.doc_mut().apply_scene_tx("Insert Object (Table)", |s, nid| {
                                        let id = *nid; *nid += 1;
                                        s.tables.push(SceneTable {
                                            id,
                                            rows,
                                            cols,
                                            bbox: (120.0, 120.0, 120.0 + (cols as f64) * 48.0, 120.0 + (rows as f64) * 26.0),
                                            cells: vec![vec![String::new(); cols]; rows],
                                        });
                                    });
                                    self.show_object_inspector = true;
                                    self.status_msg = format!("Inserted {rows}x{cols} Table object");
                                    self.dialog = ActiveDialog::None;
                                }
                                1 => {
                                    self.doc_mut().apply_scene_tx("Insert Object (TLC Plate)", |s, nid| {
                                        let id = *nid; *nid += 1;
                                        s.tlc_plates.push(SceneTlcPlate {
                                            id,
                                            bbox: (150.0, 110.0, 196.0, 194.0),
                                            origin_fraction: 0.15,
                                            spot_rfs: vec![0.25, 0.55, 0.80],
                                            spot_labels: vec!["SM".to_string(), "Co".to_string(), "Rxn".to_string()],
                                            spot_colors: vec![(20, 20, 20), (25, 95, 205), (205, 45, 45)],
                                            solvent_label: "EtOAc/Hex 1:4".to_string(),
                                        });
                                    });
                                    self.show_object_inspector = true;
                                    self.status_msg = "Inserted TLC Plate object".to_string();
                                    self.dialog = ActiveDialog::None;
                                }
                                2 => {
                                    self.doc_mut().apply_scene_tx("Insert Object (Frame)", |s, nid| {
                                        let id = *nid; *nid += 1;
                                        s.graphics.push(SceneGraphic {
                                            id,
                                            graphic_type: "Rectangle".to_string(),
                                            sub_type: Some("Shadow".to_string()),
                                            bbox: (110.0, 110.0, 250.0, 190.0),
                                        });
                                    });
                                    self.status_msg = "Inserted Framed Box object".to_string();
                                    self.dialog = ActiveDialog::None;
                                }
                                3 => {
                                    let txt = self.insert_obj_text.clone();
                                    self.paste_caption_tx(txt);
                                    self.dialog = ActiveDialog::None;
                                }
                                _ => {
                                    self.dialog = ActiveDialog::InsertFile;
                                }
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::ApplyStylesheetFile => {
                egui::Window::new("Apply Settings From Stylesheet / Document").collapsible(false).default_width(420.0).show(ctx, |ui| {
                    ui.label("Path to .cds, .cdx, or .cdxml file:");
                    ui.text_edit_singleline(&mut self.file_path_buf);
                    ui.horizontal(|ui| {
                        if ui.button("Apply").clicked() {
                            match load_document_any(&self.file_path_buf) {
                                Ok(d) => {
                                    replace_settings_tx(self.doc_mut(), d.settings, "Apply Settings");
                                    self.status_msg = format!("Applied settings from {}", self.file_path_buf);
                                    self.dialog = ActiveDialog::None;
                                }
                                Err(e) => self.status_msg = format!("Failed to load settings: {e}"),
                            }
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::BondProperties { mol_idx, bond_idx } => {
                egui::Window::new("Bond Properties").collapsible(false).show(ctx, |ui| {
                    let orig = self
                        .doc()
                        .scene
                        .molecules
                        .get(mol_idx)
                        .and_then(|m| m.bonds.get(bond_idx))
                        .cloned();
                    if let Some(mut b) = orig.clone() {
                        egui::ComboBox::from_label("Bond Order")
                            .selected_text(format!("{}", b.order))
                            .show_ui(ui, |ui| {
                                for (ord, lbl) in [(1, "Single (1)"), (2, "Double (2)"), (3, "Triple (3)"), (4, "Quadruple (4)")] {
                                    ui.selectable_value(&mut b.order, ord, lbl);
                                }
                            });
                        egui::ComboBox::from_label("Display Style")
                            .selected_text(format!("{:?}", b.display))
                            .show_ui(ui, |ui| {
                                for (d, lbl) in [
                                    (BondDisplay::Solid, "Plain / Solid"),
                                    (BondDisplay::Bold, "Bold"),
                                    (BondDisplay::Dash, "Dashed"),
                                    (BondDisplay::Hash, "Hashed"),
                                    (BondDisplay::WedgeBegin, "Wedged"),
                                    (BondDisplay::WedgedHashBegin, "Hashed Wedged"),
                                    (BondDisplay::HollowWedgeBegin, "Hollow Wedged"),
                                    (BondDisplay::Wavy, "Wavy (Either)"),
                                    (BondDisplay::Dative, "Dative"),
                                ] {
                                    ui.selectable_value(&mut b.display, d, lbl);
                                }
                            });
                        ui.checkbox(&mut b.ignore_warnings, "Ignore Chemical Warnings");
                        if Some(&b) != orig.as_ref() {
                            let key = egui::Id::new(("cc_bond_props", self.tab().tab_id, mol_idx, bond_idx));
                            let (gesture, mut last) = (undo_gesture_token(ctx), self.undo_coalesce);
                            coalesced_scene_tx(self.doc_mut(), gesture, &mut last, key, "Bond Properties", |s| {
                                if let Some(slot) = s.molecules.get_mut(mol_idx).and_then(|m| m.bonds.get_mut(bond_idx)) {
                                    *slot = b;
                                }
                            });
                            self.undo_coalesce = last;
                        }
                    } else {
                        ui.label("No bond selected.");
                    }
                    if ui.button("OK").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::RLogicQuery => {
                egui::Window::new("R-Logic Query").collapsible(false).show(ctx, |ui| {
                    ui.horizontal(|ui| { ui.label("R-Group Label:"); ui.text_edit_singleline(&mut self.rlogic_group); });
                    ui.horizontal(|ui| { ui.label("Occurrence (e.g. 1, >0, 1-2):"); ui.text_edit_singleline(&mut self.rlogic_occurrence); });
                    ui.checkbox(&mut self.rlogic_rest_h, "Rest H (unoccupied sites are H)");
                    ui.horizontal(|ui| {
                        if ui.button("Apply").clicked() {
                            let tag = format!(
                                "R-Logic: {} occ={} restH={}",
                                self.rlogic_group,
                                self.rlogic_occurrence,
                                if self.rlogic_rest_h { "yes" } else { "no" }
                            );
                            self.tab_mut().annotations.insert("RLogic".to_string(), tag.clone());
                            self.paste_caption_tx(tag.clone());
                            self.status_msg = format!("Applied {tag}");
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::LineSpacing => {
                egui::Window::new("Line Spacing").collapsible(false).show(ctx, |ui| {
                    ui.add(egui::DragValue::new(&mut self.line_spacing_pct).range(50.0..=300.0).suffix("%").prefix("Line Spacing: "));
                    ui.horizontal(|ui| {
                        if ui.button("OK").clicked() {
                            self.status_msg = format!("Set caption line spacing to {:.0}%", self.line_spacing_pct);
                            self.dialog = ActiveDialog::None;
                        }
                        if ui.button("Cancel").clicked() { self.dialog = ActiveDialog::None; }
                    });
                });
            }
            ActiveDialog::FontSizeOther => {
                egui::Window::new("Other Font Size").collapsible(false).show(ctx, |ui| {
                    ui.add(
                        egui::DragValue::new(&mut self.doc_mut().settings.label_size_pt)
                            .range(4.0..=144.0)
                            .suffix(" pt")
                            .prefix("Size: "),
                    );
                    if ui.button("OK").clicked() {
                        let sz = self.doc().settings.label_size_pt;
                        self.doc_mut().settings.caption_size_pt = sz;
                        self.status_msg = format!("Set font size to {sz:.1} pt");
                        self.dialog = ActiveDialog::None;
                    }
                });
            }
            ActiveDialog::HelpTopic { title, body } => {
                egui::Window::new(&title).collapsible(false).default_width(420.0).show(ctx, |ui| {
                    ui.label(&body);
                    ui.add_space(6.0);
                    if ui.button("Close").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::Licensing => {
                egui::Window::new("ChiralCraft Licensing").collapsible(false).default_width(420.0).show(ctx, |ui| {
                    ui.strong(format!("ChiralCraft {} — GNU Affero General Public License v3.0", env!("CARGO_PKG_VERSION")));
                    ui.label("License: AGPL-3.0-or-later");
                    ui.label("Copyright (c) 2026 nosidak (https://github.com/nosidak/chiralcraft)");
                    ui.label("Full pure-Rust chemistry drawing, spectroscopy, MedChem, 3D/quantum, and 19-tool MCP studio.");
                    if ui.button("OK").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::About => {
                if self.app_icon_tex.is_none() {
                    if let Some((rgba, w, h)) = decode_png_rgba(APP_ICON_PNG) {
                        let img = egui::ColorImage::from_rgba_unmultiplied([w as usize, h as usize], &rgba);
                        self.app_icon_tex = Some(ctx.load_texture("chiralcraft_app_icon", img, egui::TextureOptions::LINEAR));
                    }
                }
                egui::Window::new("About ChiralCraft").collapsible(false).default_width(420.0).show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        if let Some(ref tex) = self.app_icon_tex {
                            ui.image((tex.id(), Vec2::splat(72.0)));
                            ui.add_space(8.0);
                        }
                        ui.vertical(|ui| {
                            ui.heading(format!("ChiralCraft {}", env!("CARGO_PKG_VERSION")));
                            ui.label("Pure-Rust scientific chemical structure & cheminformatics studio: IUPAC 2008 rendering, 16 journal stylesheets, 11 template libraries, live PubChem/ChEMBL search, Editable SVG/PNG/PDF round-trip, ¹H/¹³C/2D NMR, FT-IR, UV–Vis, MedChem/HRMS, Isomers, Markush SAR, Retrosynthesis, 3D UFF/Crystals, ORCA/Gaussian/xTB bridge, Organometallic 18e⁻, LaTeX, ELN SI & 19-tool MCP server.");
                            ui.small("AGPL-3.0-or-later • github.com/nosidak/chiralcraft");
                        });
                    });
                    ui.add_space(6.0);
                    if ui.button("OK").clicked() { self.dialog = ActiveDialog::None; }
                });
            }
            ActiveDialog::ConfirmCloseTab => {
                let title = self.doc().name.clone();
                egui::Window::new("Unsaved Changes").collapsible(false).resizable(false).show(ctx, |ui| {
                    ui.label(format!("Save changes to \"{title}\" before closing?"));
                    ui.horizontal(|ui| {
                        if ui.button("Save").clicked() {
                            self.dialog = ActiveDialog::None;
                            self.do_save_action();
                            // Closes only once the save really happened (a Save As prompt or a failed
                            // write leaves the tab open).
                            if !self.doc().is_dirty() {
                                self.force_close_active_tab();
                            }
                        }
                        if ui.button("Don't Save").clicked() {
                            self.force_close_active_tab();
                        }
                        if ui.button("Cancel").clicked() {
                            self.dialog = ActiveDialog::None;
                        }
                    });
                });
            }
            ActiveDialog::ConfirmExit => {
                let dirty: Vec<String> = self.tabs.iter().filter(|t| t.doc.is_dirty()).map(|t| t.doc.name.clone()).collect();
                egui::Window::new("Unsaved Changes").collapsible(false).resizable(false).show(ctx, |ui| {
                    ui.label(format!("{} document(s) have unsaved changes:", dirty.len()));
                    for n in &dirty {
                        ui.label(format!("• {n}"));
                    }
                    ui.horizontal(|ui| {
                        if ui.button("Discard Changes & Exit").clicked() {
                            self.exit_confirmed = true;
                            self.dialog = ActiveDialog::None;
                            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                        }
                        if ui.button("Cancel").clicked() {
                            self.dialog = ActiveDialog::None;
                        }
                    });
                });
            }
            ActiveDialog::None => {}
        }
    }
}

fn print_cli_help() {
    println!(
        "ChiralCraft {} — Pure-Rust Chemical Structure & Cheminformatics Studio\n\
         \n\
         USAGE:\n\
           chiralcraft [FILE]                       Launch desktop GUI (optionally opening FILE)\n\
           chiralcraft --edit-svg <FILE.svg>        Edit an embedded ChiralCraft SVG round-trip\n\
           chiralcraft --mcp                        Run stdio JSON-RPC 2.0 MCP server (19 tools)\n\
           chiralcraft --help | -h | --version | -V Show help or version\n\
         \n\
         CLI SUBCOMMANDS (accept any supported file OR raw SMILES):\n\
           convert <IN> <OUT>                       Convert between .cdxml/.cdx/.mol/.sdf/.rxn/.xyz/.cif/.pdb/POSCAR/.svg/.png/.pdf/.eps/.emf/.cml/.ct/.tgf/.skc/.sln/.inchi\n\
           analyze <IN>                             Formula, Exact Mass, MW, LogP, tPSA, LogS, pKa\n\
           diff <LEFT> <RIGHT>                      Chemistry-aware diff between two files\n\
           nmr <IN> [--13c]                         Predict 1H-NMR (default) or 13C-NMR shifts\n\
           ir <IN>                                  Predict FT-IR bands and UV-Vis lambda_max\n\
           medchem <IN>                             Lipinski, Veber, Ghose, CNS MPO, QED, SA, PAINS/Brenk\n\
           hrms <IN> [--adduct <ADDUCT>]            HRMS isotope envelope & ESI adducts (default [M+H]+)\n\
           retro <IN>                               Retrosynthesis disconnections & reagents\n\
           18e <IN>                                 Organometallic 18-electron & Green CBC analyzer\n\
           export-tex <IN>                          Export LaTeX \\chemfig / \\schemestart snippet\n\
           screen <LIB.smi|csv|sdf>                 Batch MedChem & PAINS screening table\n\
           quantum-in <IN> [--engine orca|gaussian|xtb]\n\
                                                    Generate ORCA, Gaussian, or GFN2-xTB input deck\n\
           quantum-parse <OUT|LOG>                  Parse ORCA .out or Gaussian .log energies & IR\n\
           inspect-svg | inspect-media <FILE>       Verify embedded editable ChiralCraft payload",
        env!("CARGO_PKG_VERSION")
    );
}

fn main() -> eframe::Result<()> {
    std::panic::set_hook(Box::new(|info| {
        let bt = std::backtrace::Backtrace::force_capture();
        let msg = format!("ChiralCraft panic: {info}\n{bt}\n");
        let crash_path = crash_log_path();
        if let Some(parent) = crash_path.parent() {
            let _ = fs::create_dir_all(parent);
        }
        let _ = fs::write(&crash_path, &msg);
        eprintln!("{msg}");
    }));
    let args: Vec<String> = env::args().collect();
    fn require_first_mol(doc: &InteractiveDocument) -> &ChemMolecule {
        doc.scene.molecules.first().unwrap_or_else(|| {
            eprintln!("Error: input document contains no molecules.");
            std::process::exit(1);
        })
    }
    match args.get(1).map(|s| s.as_str()) {
        Some("--help") | Some("-h") | Some("help") => {
            print_cli_help();
            Ok(())
        }
        Some("--version") | Some("-V") | Some("version") => {
            println!("chiralcraft {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        Some("convert") => {
            if args.len() < 4 {
                eprintln!("Usage: chiralcraft convert <IN> <OUT>");
                std::process::exit(2);
            }
            let doc = load_document_or_smiles(&args[2]).unwrap_or_else(|e| { eprintln!("Load error: {e}"); std::process::exit(1); });
            if let Err(e) = save_document_any(&doc, &args[3]) { eprintln!("Save error: {e}"); std::process::exit(1); }
            Ok(())
        }
        Some("analyze") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft analyze <FILE|SMILES>");
                std::process::exit(2);
            }
            let doc = load_document_or_smiles(&args[2]).unwrap_or_else(|e| { eprintln!("Load error: {e}"); std::process::exit(1); });
            let m = require_first_mol(&doc);
            let rep = chiralcraft_chem::compute_analysis_report(m);
            let p = estimate_chemical_properties(m);
            let pka_str = p.pka.map(|v| format!("{v:.2}")).unwrap_or_else(|| "N/A".to_string());
            println!("Formula: {} | Exact Mass: {:.4} | Mol. Wt.: {:.2} | Elem. Anal.: {}", rep.formula, rep.exact_mass, rep.mol_weight, rep.elem_analysis);
            println!("BP: {:.2} K | MP: {:.2} K | LogP: {:.2} | MR: {:.1} | tPSA: {:.2} | LogS: {:.3} | pKa: {pka_str}", p.boiling_point_k, p.melting_point_k, p.log_p, p.molar_refractivity, p.tpsa, p.log_s);
            Ok(())
        }
        Some("diff") => {
            if args.len() < 4 {
                eprintln!("Usage: chiralcraft diff <LEFT> <RIGHT>");
                std::process::exit(2);
            }
            let left = load_document_any(&args[2]).unwrap_or_else(|e| { eprintln!("Load error ({}): {e}", args[2]); std::process::exit(1); });
            let right = load_document_any(&args[3]).unwrap_or_else(|e| { eprintln!("Load error ({}): {e}", args[3]); std::process::exit(1); });
            let d = compare_documents(&left, &right);
            print!("{}", d.format_text_report());
            Ok(())
        }
        Some("nmr") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft nmr <FILE|SMILES> [--13c]");
                std::process::exit(2);
            }
            let carbon13 = args[2..].iter().any(|a| a == "--13c" || a == "-13c");
            let target = args[2..]
                .iter()
                .find(|a| *a != "--13c" && *a != "-13c")
                .map(|s| s.as_str())
                .unwrap_or_else(|| { eprintln!("Usage: chiralcraft nmr <FILE|SMILES> [--13c]"); std::process::exit(2); });
            let doc = load_document_or_smiles(target).unwrap_or_else(|e| { eprintln!("Load error: {e}"); std::process::exit(1); });
            let m = require_first_mol(&doc);
            let pred = if carbon13 { predict_13c_nmr(m) } else { predict_1h_nmr(m) };
            println!("{}", pred.report_text);
            for pk in &pred.peaks {
                println!("  atom #{}: δ {:.2} ppm ({}, n={}) — {}", pk.atom_idx + 1, pk.shift_ppm, pk.multiplicity, pk.protons, pk.environment);
            }
            Ok(())
        }
        Some("medchem") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft medchem <FILE|SMILES>");
                std::process::exit(2);
            }
            let doc = load_document_or_smiles(&args[2]).unwrap_or_else(|e| { eprintln!("Load error: {e}"); std::process::exit(1); });
            let m = require_first_mol(&doc);
            let mc = compute_medchem_report(m);
            println!("MW: {:.2} | LogP: {:.2} | tPSA: {:.1} | HBD: {} | HBA: {} | RotB: {} | Fsp3: {:.2}", mc.mw, mc.log_p, mc.tpsa, mc.hbd, mc.hba, mc.rotatable_bonds, mc.fsp3);
            println!("Lipinski: {} ({} viol) | Veber: {} | Ghose: {} | CNS MPO: {:.2}/6.0 | QED: {:.2} | SA: {:.2}/10", if mc.lipinski_pass { "PASS" } else { "FAIL" }, mc.lipinski_violations, if mc.veber_pass { "PASS" } else { "FAIL" }, if mc.ghose_pass { "PASS" } else { "FAIL" }, mc.cns_mpo, mc.qed, mc.sa_score);
            for al in &mc.alerts {
                println!("  ALERT [{}]: {}", al.category, al.name);
            }
            Ok(())
        }
        Some("hrms") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft hrms <FILE|SMILES> [--adduct <ADDUCT>]");
                std::process::exit(2);
            }
            let mut target: Option<&str> = None;
            let mut adduct: &str = "[M+H]+";
            let mut i = 2;
            while i < args.len() {
                if args[i] == "--adduct" || args[i] == "-a" {
                    if i + 1 >= args.len() {
                        eprintln!("Error: {} requires a value", args[i]);
                        std::process::exit(2);
                    }
                    adduct = &args[i + 1];
                    i += 2;
                } else if target.is_none() {
                    target = Some(&args[i]);
                    i += 1;
                } else {
                    adduct = &args[i];
                    i += 1;
                }
            }
            let target = target.unwrap_or_else(|| { eprintln!("Usage: chiralcraft hrms <FILE|SMILES> [--adduct <ADDUCT>]"); std::process::exit(2); });
            let doc = load_document_or_smiles(target).unwrap_or_else(|e| { eprintln!("Load error: {e}"); std::process::exit(1); });
            let m = require_first_mol(&doc);
            let sim = simulate_hrms(m, adduct);
            println!("Formula: {} | Neutral Exact Mass: {:.4} | Adduct: {}", sim.formula, sim.neutral_exact_mass, sim.selected_adduct);
            for pk in sim.peaks.iter().take(8) {
                println!("  m/z {:.4}  ({:.1}%)", pk.mz, pk.rel_intensity);
            }
            Ok(())
        }
        Some("retro") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft retro <FILE|SMILES>");
                std::process::exit(2);
            }
            let doc = load_document_or_smiles(&args[2]).unwrap_or_else(|e| { eprintln!("Load error: {e}"); std::process::exit(1); });
            let m = require_first_mol(&doc);
            let hits = suggest_disconnections(m, doc.settings.bond_length_pt);
            if hits.is_empty() {
                println!("No acyclic strategic disconnections suggested.");
            } else {
                for (i, h) in hits.iter().enumerate() {
                    println!("{}. {} (score {}) -> {} | Reagents: {}", i + 1, h.reaction_name, h.score, h.synthons_summary, h.reagents);
                }
            }
            Ok(())
        }
        Some("ir") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft ir <FILE|SMILES>");
                std::process::exit(2);
            }
            let doc = load_document_or_smiles(&args[2]).unwrap_or_else(|e| { eprintln!("Load error: {e}"); std::process::exit(1); });
            let m = require_first_mol(&doc);
            let ir = predict_ir_spectrum(m);
            let uv = predict_uv_vis(m);
            println!("{}", ir.report_text);
            for b in &ir.bands {
                println!("  {:.0} cm^-1 ({}) — {}", b.wavenumber_cm1, b.strength, b.assignment);
            }
            println!("UV-Vis: λmax ≈ {:.0} nm (ε ≈ {:.0}), {} ({})", uv.lambda_max_nm, uv.epsilon_est, uv.transition, uv.perceived_color_name);
            Ok(())
        }
        Some("18e") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft 18e <FILE|SMILES>");
                std::process::exit(2);
            }
            let doc = load_document_or_smiles(&args[2]).unwrap_or_else(|e| { eprintln!("Load error: {e}"); std::process::exit(1); });
            let m = require_first_mol(&doc);
            let reps = analyze_organometallic_complexes(m);
            if reps.is_empty() {
                println!("No d-block transition metal centers found.");
            } else {
                for r in reps {
                    println!("{} (Group {}): {}e- | Ox {:+} | d^{} | ML{}X{} | {} — {}", r.metal_symbol, r.group_number, r.total_valence_electrons, r.oxidation_state, r.d_electrons, r.l_ligands, r.x_ligands, r.predicted_geometry, r.rule_status);
                }
            }
            Ok(())
        }
        Some("export-tex") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft export-tex <FILE|SMILES>");
                std::process::exit(2);
            }
            let doc = load_document_or_smiles(&args[2]).unwrap_or_else(|e| { eprintln!("Load error: {e}"); std::process::exit(1); });
            let _ = require_first_mol(&doc);
            print!("{}", document_to_chemfig_scheme(&doc));
            Ok(())
        }
        Some("screen") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft screen <LIB.smi|csv|sdf>");
                std::process::exit(2);
            }
            let recs = load_screening_library(&args[2], 14.40).unwrap_or_else(|e| { eprintln!("Screen error: {e}"); std::process::exit(1); });
            println!("Screened {} compound(s):", recs.len());
            for r in recs {
                println!("  #{:03} {:<20} {} | MW {:6.1} | LogP {:5.2} | MPO {:.2} | QED {:.2} | Alerts {}", r.index, r.title, r.formula, r.medchem.mw, r.medchem.log_p, r.medchem.cns_mpo, r.medchem.qed, r.medchem.alerts.len());
            }
            Ok(())
        }
        Some("quantum-in") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft quantum-in <FILE|SMILES> [--engine orca|gaussian|xtb]");
                std::process::exit(2);
            }
            let mut target: Option<&str> = None;
            let mut engine: &str = "orca";
            let mut i = 2;
            while i < args.len() {
                if args[i] == "--engine" || args[i] == "--program" || args[i] == "-e" {
                    if i + 1 >= args.len() {
                        eprintln!("Error: {} requires a value (orca|gaussian|xtb)", args[i]);
                        std::process::exit(2);
                    }
                    engine = &args[i + 1];
                    i += 2;
                } else if target.is_none() {
                    target = Some(&args[i]);
                    i += 1;
                } else {
                    engine = &args[i];
                    i += 1;
                }
            }
            let target = target.unwrap_or_else(|| { eprintln!("Usage: chiralcraft quantum-in <FILE|SMILES> [--engine orca|gaussian|xtb]"); std::process::exit(2); });
            let mut spec = QuantumJobSpec::default();
            spec.program = match engine.to_ascii_lowercase().as_str() {
                "orca" => QuantumProgram::Orca,
                "gaussian" | "gjf" => QuantumProgram::Gaussian,
                "xtb" => QuantumProgram::Xtb,
                other => {
                    eprintln!("Error: unknown quantum engine '{other}' (expected orca, gaussian, or xtb)");
                    std::process::exit(2);
                }
            };
            let doc = load_document_or_smiles(target).unwrap_or_else(|e| { eprintln!("Load error: {e}"); std::process::exit(1); });
            let m = require_first_mol(&doc);
            let lifted = lift_molecule_2d_to_3d(m);
            print!("{}", export_quantum_input(&lifted, &spec));
            Ok(())
        }
        Some("quantum-parse") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft quantum-parse <OUT|LOG|MOLDEN>");
                std::process::exit(2);
            }
            let txt = fs::read_to_string(&args[2]).unwrap_or_else(|e| { eprintln!("Read error: {e}"); std::process::exit(1); });
            let rep = parse_quantum_output(&txt).unwrap_or_else(|e| { eprintln!("Parse error: {e}"); std::process::exit(1); });
            println!("Program: {} | Converged: {} | SCF (Eh): {:?} | HOMO-LUMO Gap (eV): {:?} | Dipole (D): {:?} | IR modes: {}", rep.program, rep.converged, rep.energy_hartree, rep.homo_lumo_gap_ev, rep.dipole_debye, rep.ir_modes.len().max(rep.ir_frequencies_cm1.len()));
            Ok(())
        }
        Some("inspect-svg") | Some("inspect-media") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft {} <FILE.svg|.png|.pdf>", args[1]);
                std::process::exit(2);
            }
            let low = args[2].to_ascii_lowercase();
            if !(low.ends_with(".svg") || low.ends_with(".png") || low.ends_with(".pdf")) {
                eprintln!("editable=false error=\"inspect-svg/inspect-media requires a .svg, .png, or .pdf file\"");
                std::process::exit(2);
            }
            let res = load_document_any(&args[2]);
            match res {
                Ok(doc) => {
                    println!(
                        "editable=true molecules={} captions={} stylesheet=\"{}\"",
                        doc.scene.molecules.len(),
                        doc.scene.captions.len(),
                        doc.settings.stylesheet_name
                    );
                    Ok(())
                }
                Err(e) => {
                    eprintln!("editable=false error=\"{e}\"");
                    std::process::exit(2);
                }
            }
        }
        Some("--mcp") => {
            mcp::run_mcp_stdio_server();
            Ok(())
        }
        Some("--screenshot") => {
            if args.len() < 4 {
                eprintln!("Usage: chiralcraft --screenshot <PRESET> <OUT.png>");
                std::process::exit(2);
            }
            let preset = args[2].clone();
            let out_path = args[3].clone();
            let mut vp = egui::ViewportBuilder::default()
                .with_inner_size([1360.0, 840.0])
                .with_title("ChiralCraft");
            if let Some(icon) = load_app_icon_data() {
                vp = vp.with_icon(std::sync::Arc::new(icon));
            }
            let opts = eframe::NativeOptions {
                viewport: vp,
                ..Default::default()
            };
            eframe::run_native(
                "ChiralCraft",
                opts,
                Box::new(move |cc| Ok(Box::new({ install_glyph_fallback(&cc.egui_ctx); ChiralCraftApp::new_showcase_preset(&preset, &out_path) }))),
            )
        }
        Some("--edit-svg") => {
            if args.len() < 3 {
                eprintln!("Usage: chiralcraft --edit-svg <FILE.svg>");
                std::process::exit(2);
            }
            let target_path = args[2].clone();
            let fname = std::path::Path::new(&target_path)
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("embedded.svg");
            let title = format!("ChiralCraft — Editing {fname}");
            let mut vp = egui::ViewportBuilder::default()
                .with_inner_size([1240.0, 820.0])
                .with_title(title.clone());
            if let Some(icon) = load_app_icon_data() {
                vp = vp.with_icon(std::sync::Arc::new(icon));
            }
            let opts = eframe::NativeOptions {
                viewport: vp,
                ..Default::default()
            };
            eframe::run_native(
                &title,
                opts,
                Box::new(move |cc| Ok(Box::new({ install_glyph_fallback(&cc.egui_ctx); ChiralCraftApp::new_with_file(&target_path, true) }))),
            )
        }
        Some(flag) if flag.starts_with('-') => {
            eprintln!("Unknown option: {flag}\nRun `chiralcraft --help` for usage.");
            std::process::exit(2);
        }
        Some(path) if !path.starts_with('-') => {
            let open_path = path.to_string();
            let mut vp = egui::ViewportBuilder::default()
                .with_inner_size([1360.0, 860.0])
                .with_title("ChiralCraft");
            if let Some(icon) = load_app_icon_data() {
                vp = vp.with_icon(std::sync::Arc::new(icon));
            }
            let opts = eframe::NativeOptions {
                viewport: vp,
                ..Default::default()
            };
            eframe::run_native(
                "ChiralCraft",
                opts,
                Box::new(move |cc| Ok(Box::new({ install_glyph_fallback(&cc.egui_ctx); ChiralCraftApp::new_with_file(&open_path, false) }))),
            )
        }
        _ => {
            let mut vp = egui::ViewportBuilder::default()
                .with_inner_size([1360.0, 860.0])
                .with_title("ChiralCraft");
            if let Some(icon) = load_app_icon_data() {
                vp = vp.with_icon(std::sync::Arc::new(icon));
            }
            let opts = eframe::NativeOptions {
                viewport: vp,
                ..Default::default()
            };
            eframe::run_native(
                "ChiralCraft",
                opts,
                Box::new(|cc| Ok(Box::new({ install_glyph_fallback(&cc.egui_ctx); ChiralCraftApp::new() }))),
            )
        }
    }
}

/// egui's proportional font (Ubuntu-Light) lacks ⁺ ⁻ → ⇌ ⊕ etc., which rendered as tofu boxes in
/// atom charges and labels; fall back to the bundled Hack font, which has them.
fn install_glyph_fallback(ctx: &egui::Context) {
    let mut fonts = egui::FontDefinitions::default();
    if let Some(prop) = fonts.families.get_mut(&egui::FontFamily::Proportional) {
        prop.push("Hack".to_owned());
    }
    ctx.set_fonts(fonts);
}

/// Stamps several molecules (multi-record paste/insert) as ONE undo transaction.
fn stamp_molecules_tx(doc: &mut InteractiveDocument, mols: Vec<ChemMolecule>, label: &str) {
    if mols.is_empty() {
        return;
    }
    doc.apply_scene_tx(label, |scene, next_id| {
        for mut mol in mols {
            for a in &mut mol.atoms {
                a.id = *next_id;
                *next_id += 1;
            }
            for b in &mut mol.bonds {
                b.id = *next_id;
                *next_id += 1;
            }
            scene.molecules.push(mol);
        }
    });
}

/// Identifies the widget gesture in progress: the pressed pointer (a drag) or else the
/// keyboard-focused widget (typing into a DragValue/TextEdit). `None` once the gesture ends.
fn undo_gesture_token(ctx: &egui::Context) -> Option<egui::Id> {
    ctx.input(|i| {
        if i.pointer.any_down() {
            i.pointer.press_origin().map(|p| egui::Id::new(("cc_press", p.x.to_bits(), p.y.to_bits())))
        } else {
            None
        }
    })
    .or_else(|| ctx.memory(|m| m.focused()))
}

/// Records a scene edit made from a live dialog widget. Edits continuing the same gesture
/// (same `key`, same drag/focus) are merged into the previous step, so one drag or one typed
/// value is ONE undo entry instead of one per frame/keystroke.
fn coalesced_scene_tx(
    doc: &mut InteractiveDocument,
    gesture: Option<egui::Id>,
    last: &mut Option<(egui::Id, egui::Id)>,
    key: egui::Id,
    label: &str,
    f: impl FnOnce(&mut chiralcraft_doc::DocumentScene),
) {
    let mut next = doc.scene.clone();
    f(&mut next);
    let merge = gesture.is_some()
        && *last == gesture.map(|g| (key, g))
        && doc.is_dirty()
        && doc.undo_menu_label() == format!("Undo {label}")
        && doc.redo_menu_label() == "Redo not available";
    if merge {
        // undo() restores the pre-gesture settings; keep the live widget values for this step.
        let live = doc.settings.clone();
        doc.undo();
        doc.settings = live;
    }
    doc.record_transaction(label, next);
    *last = gesture.map(|g| (key, g));
}

/// Records a settings-only step (scene unchanged); the doc's undo history snapshots settings,
/// so Undo restores them and the tab is marked unsaved.
fn coalesced_settings_changed(
    doc: &mut InteractiveDocument,
    gesture: Option<egui::Id>,
    last: &mut Option<(egui::Id, egui::Id)>,
    key: egui::Id,
    label: &str,
) {
    coalesced_scene_tx(doc, gesture, last, key, label, |_| {});
}

/// Replaces the document settings (stylesheet apply) as one undoable step.
fn replace_settings_tx(doc: &mut InteractiveDocument, settings: DocumentSettings, label: &str) {
    if doc.settings != settings {
        doc.settings = settings;
        doc.apply_scene_tx(label, |_, _| {});
    }
}

/// True when writing `doc` to `path` would drop content (format keeps only molecules, only the
/// first molecule, or cannot be re-opened). Such a "Save As" is an export: the tab must not be
/// retargeted to that file nor marked clean, or closing silently loses the rest of the drawing.
fn save_as_format_is_lossy(doc: &InteractiveDocument, path: &str) -> bool {
    let low = path.to_ascii_lowercase();
    let full = [".cdxml", ".cdx", ".cds", ".svg", ".png", ".pdf"];
    if full.iter().any(|e| low.ends_with(e)) {
        return false;
    }
    let write_only = [".eps", ".emf", ".sln", ".inchi"];
    if write_only.iter().any(|e| low.ends_with(e)) {
        return true;
    }
    let s = &doc.scene;
    let has_non_molecule = !s.arrows.is_empty()
        || !s.captions.is_empty()
        || !s.curves.is_empty()
        || !s.graphics.is_empty()
        || !s.tables.is_empty()
        || !s.tlc_plates.is_empty();
    let multi_mol_ok = [".sdf", ".sd", ".smi", ".smiles", ".rxn"].iter().any(|e| low.ends_with(e));
    has_non_molecule || (!multi_mol_ok && s.molecules.len() > 1)
}

#[cfg(test)]
mod audit_b_tests {
    use super::*;

    #[test]
    fn test_stamp_molecules_tx_is_single_undo_with_unique_ids() {
        let mut d = InteractiveDocument::new_saved_seed("t", DocumentSettings::acs_document_1996());
        let m1 = smiles_to_molecule("CCO").unwrap();
        let m2 = smiles_to_molecule("CC").unwrap();
        stamp_molecules_tx(&mut d, vec![m1, m2], "Paste MOL/SDF");
        assert_eq!(d.scene.molecules.len(), 2);
        let mut ids: Vec<u32> = d.scene.molecules.iter().flat_map(|m| m.atoms.iter().map(|a| a.id)).collect();
        let n = ids.len();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), n);
        d.undo();
        assert!(d.scene.molecules.is_empty(), "one Undo must remove the whole paste");
    }

    #[test]
    fn test_coalesced_scene_tx_one_undo_per_gesture() {
        let mut d = InteractiveDocument::new_saved_seed("t", DocumentSettings::acs_document_1996());
        d.scene.molecules.push(smiles_to_molecule("CCO").unwrap());
        d.save();
        let key = egui::Id::new("k");
        let (g1, g2) = (Some(egui::Id::new("drag1")), Some(egui::Id::new("drag2")));
        let mut last = None;
        for c in 1..=5 {
            coalesced_scene_tx(&mut d, g1, &mut last, key, "Atom Properties", |s| s.molecules[0].atoms[0].charge = c);
        }
        coalesced_scene_tx(&mut d, g2, &mut last, key, "Atom Properties", |s| s.molecules[0].atoms[0].isotope = 13);
        assert!(d.is_dirty());
        assert!(d.undo());
        assert_eq!(d.scene.molecules[0].atoms[0].charge, 5, "second gesture is its own step");
        assert!(d.undo());
        assert_eq!(d.scene.molecules[0].atoms[0].charge, 0, "whole first drag undone in one step");
        assert!(!d.undo());
        assert!(!d.is_dirty());
    }

    #[test]
    fn test_replace_settings_tx_marks_unsaved_only_on_change() {
        let mut d = InteractiveDocument::new_saved_seed("t", DocumentSettings::acs_document_1996());
        replace_settings_tx(&mut d, DocumentSettings::acs_document_1996(), "Apply");
        assert!(!d.is_dirty());
        replace_settings_tx(&mut d, DocumentSettings::standard_30pt_document(), "Apply");
        assert!(d.is_dirty());
    }

    #[test]
    fn test_save_as_lossy_format_detection() {
        let mut d = InteractiveDocument::new_saved_seed("t", DocumentSettings::acs_document_1996());
        d.scene.molecules.push(smiles_to_molecule("CCO").unwrap());
        assert!(!save_as_format_is_lossy(&d, "a.mol"));
        assert!(!save_as_format_is_lossy(&d, "a.CDXML"));
        assert!(save_as_format_is_lossy(&d, "a.eps"));
        d.scene.molecules.push(smiles_to_molecule("CC").unwrap());
        assert!(save_as_format_is_lossy(&d, "a.mol"));
        assert!(!save_as_format_is_lossy(&d, "a.sdf"));
        assert!(!save_as_format_is_lossy(&d, "a.cdx"));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_desktop_app_ring_fusion_spiro_chain_cdx_and_aromatic() {
        let mut app = ChiralCraftApp::new();
        assert_eq!(app.doc().nicknames.len(), 537);
        assert_eq!(app.ctp_libraries.len(), 11, "All 11 clean-room built-in template libraries must load");
        assert_eq!(app.stylesheets.len(), 16, "All 16 built-in journal stylesheets must load");
        // Verify clean-room BioGraphics & Lab Shapes graphic library loads & stamps curves!
        let biogfx = app.ctp_libraries.iter().find(|l| l.name == "BioGraphics & Lab Shapes").expect("BioGraphics & Lab Shapes must exist");
        assert!(biogfx.total_count() >= 8);
        let stamped_curves = biogfx.stamp_graphic(0, 200.0, 200.0).expect("BioGraphics curve must stamp");
        assert!(!stamped_curves.is_empty());
        // Starts blank like real ChiralCraft (`Untitled ACS Document 1996-1`)!
        assert_eq!(app.doc().scene.molecules.len(), 0);
        let blen = app.doc().settings.bond_length_pt;
        app.doc_mut().stamp_molecule_tx(stamp_template(BuiltinTemplate::Benzene, 180.0, 150.0, blen), "Benzene");
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), 6);

        // Ring-interior double bond detection on Benzene
        let rn = ring_bond_inward_normals(&app.doc().scene.molecules[0]);
        assert_eq!(rn.len(), 6);

        // Edge fusion -> Naphthalene (10 atoms)
        let b0 = app.doc().scene.molecules[0].bonds[0].clone();
        let mx = 0.5 * (app.doc().scene.molecules[0].atoms[b0.begin_atom].x + app.doc().scene.molecules[0].atoms[b0.end_atom].x);
        let my = 0.5 * (app.doc().scene.molecules[0].atoms[b0.begin_atom].y + app.doc().scene.molecules[0].atoms[b0.end_atom].y);
        app.tool = Tool::Benzene;
        app.handle_canvas_click(mx, my, false);
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), 10);

        // Undo + Spiro fusion (6 + 5 - 1 = 10 atoms)
        assert!(app.doc_mut().undo());
        let ax = app.doc().scene.molecules[0].atoms[0].x;
        let ay = app.doc().scene.molecules[0].atoms[0].y;
        app.tool = Tool::Cyclopentane;
        app.handle_canvas_click(ax, ay, false);
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), 10);

        // Acyclic chain
        app.tool = Tool::AcyclicChain;
        app.handle_canvas_drag_release(300.0, 100.0, 350.0, 100.0, false);
        assert_eq!(app.doc().scene.molecules[1].atoms.len(), 5);

        // Freehand Lasso polygon selection around chain (300..360, 80..120)
        app.tool = Tool::Lasso;
        app.lasso_pts = vec![(290.0, 75.0), (370.0, 75.0), (370.0, 125.0), (290.0, 125.0)];
        app.handle_canvas_drag_release(290.0, 75.0, 370.0, 125.0, false);
        assert_eq!(app.tab().selected_mols, vec![1]);

        // Mass Fragmentation cleavage across chain middle bond
        let frags = fragment_molecule_by_line(&app.doc().scene.molecules[1], (325.0, 70.0), (325.0, 130.0));
        assert_eq!(frags.len(), 2);

        // Reaction Atom Mapping & Clean Up Reaction
        app.tool = Tool::Arrow;
        app.handle_canvas_drag_release(230.0, 150.0, 270.0, 150.0, false);
        assert!(clean_up_reaction(&mut app.doc_mut().scene));
        let maps = auto_map_reaction_atoms(&app.doc().scene);
        assert!(maps.len() >= 10);

        // Multi-document MDI tabs, move/resize geometry, & stale-index crash protection
        app.new_document_tab(DocumentSettings::standard_30pt_document());
        assert_eq!(app.tabs.len(), 2);
        assert_eq!(app.doc().settings.bond_length_pt, 30.0);
        app.tab_mut().win_pos += Vec2::new(60.0, 40.0);
        app.tab_mut().win_size = Vec2::new(500.0, 400.0);
        app.tab_mut().maximized = true;
        app.hovered_atom = Some((99, 99));
        app.hovered_bond = Some((99, 99));
        app.inline_atom_edit = Some((99, 99));
        app.sanitize_state();
        assert_eq!(app.hovered_atom, None);
        assert_eq!(app.hovered_bond, None);
        assert_eq!(app.inline_atom_edit, None);
        app.close_active_tab();
        assert_eq!(app.tabs.len(), 1, "a clean tab closes without asking");

        // Aromatic toggle & binary CDX roundtrip
        app.toggle_aromatic_display_tx();
        assert!(!app.doc().scene.aromatic_circles.is_empty());

        let cdx_bytes = app.doc().to_cdxml().to_cdx().to_bytes();
        let reparsed = chiralcraft_formats::CdxDocument::parse(&cdx_bytes).unwrap().to_cdxml();
        let rt_doc = InteractiveDocument::from_cdxml(&reparsed);
        assert_eq!(rt_doc.scene.molecules.len(), 2);
    }

    #[test]
    fn test_hotkeys_xml_98_snapping_handles_inchi_and_rxn_sdf() {
        let mut app = ChiralCraftApp::new();
        let blen = app.doc().settings.bond_length_pt;
        app.doc_mut().stamp_molecule_tx(stamp_template(BuiltinTemplate::Benzene, 180.0, 150.0, blen), "Benzene");
        // Atom hotkeys: case-sensitive `C` -> Cl (Z=17), `F` -> CF3 label, `2` -> carbonyl sprout, `+` -> charge+1
        assert!(matches!(apply_atom_hotkey(app.doc_mut(), 0, 0, 'C'), HotkeyOutcome::Handled(_)));
        assert_eq!(app.doc().scene.molecules[0].atoms[0].z, 17);
        assert!(matches!(apply_atom_hotkey(app.doc_mut(), 0, 1, 'F'), HotkeyOutcome::Handled(_)));
        assert_eq!(app.doc().scene.molecules[0].atoms[1].label.as_deref(), Some("CF3"));
        let n_before = app.doc().scene.molecules[0].atoms.len();
        assert!(matches!(apply_atom_hotkey(app.doc_mut(), 0, 2, '2'), HotkeyOutcome::Handled(_)));
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), n_before + 1);
        assert_eq!(app.doc().scene.molecules[0].atoms.last().unwrap().z, 8);

        // Bond hotkeys: `w` -> WedgeBegin, `l` -> Left double align, `a` -> fuse benzene
        assert!(matches!(apply_bond_hotkey(app.doc_mut(), 0, 0, 'w'), HotkeyOutcome::Handled(_)));
        assert_eq!(app.doc().scene.molecules[0].bonds[0].display, BondDisplay::WedgeBegin);
        assert_eq!(apply_bond_hotkey(app.doc_mut(), 0, 0, 'l'), HotkeyOutcome::SetDoubleAlign((0, 0), DoubleBondAlign::Left));
        let n_pre_fuse = app.doc().scene.molecules[0].atoms.len();
        assert!(matches!(apply_bond_hotkey(app.doc_mut(), 0, 1, 'a'), HotkeyOutcome::Handled(_)));
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), n_pre_fuse + 4);

        // Generic hotkeys: `X` -> AcyclicChain, `=` -> ChooseNickname, `/` -> Properties
        assert_eq!(apply_generic_hotkey('X'), HotkeyOutcome::SwitchTool(Tool::AcyclicChain));
        assert_eq!(apply_generic_hotkey('='), HotkeyOutcome::OpenChooseNickname);
        assert_eq!(apply_generic_hotkey('/'), HotkeyOutcome::OpenPropertiesDialog);

        // Fixed Lengths & Fixed Angles snapping
        let (sx, sy) = snap_bond_endpoint(100.0, 100.0, 140.0, 104.0, true, true, 14.40);
        assert!(((sx - 100.0).hypot(sy - 100.0) - 14.40).abs() < 1e-6);
        assert!((sy - 100.0).abs() < 1e-6); // snapped to 0 deg!

        // Interactive scale & rotate handle math
        scale_selected_molecules(&mut app.doc_mut().scene.molecules, &[0], 180.0, 150.0, 2.0, 2.0);
        rotate_selected_molecules(&mut app.doc_mut().scene.molecules, &[0], 180.0, 150.0, std::f64::consts::FRAC_PI_2);

        // Paste Special > InChI parser (Ethanol: 2 C + 1 O = 3 heavy atoms, 2 bonds)
        let eth = inchi_to_molecule("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3").unwrap();
        assert_eq!(eth.atoms.len(), 3);
        assert_eq!(eth.bonds.len(), 2);
        assert_eq!(eth.atoms[2].z, 8);

        // Multi-molecule SDF & RXN serialization roundtrip via save/load_document_any
        let sdf_path = "/tmp/chiralcraft_test_parity.sdf";
        let rxn_path = "/tmp/chiralcraft_test_parity.rxn";
        app.doc_mut().stamp_molecule_tx(eth, "Ethanol");
        save_document_any(app.doc(), sdf_path).unwrap();
        let sdf_loaded = load_document_any(sdf_path).unwrap();
        assert_eq!(sdf_loaded.scene.molecules.len(), 2);
        save_document_any(app.doc(), rxn_path).unwrap();
        let rxn_loaded = load_document_any(rxn_path).unwrap();
        assert!(!rxn_loaded.scene.molecules.is_empty());
        let _ = fs::remove_file(sdf_path);
        let _ = fs::remove_file(rxn_path);

        // Selection-aware Cut/Copy/Clear/Flip/Scale/Rotate + Repeat + disk-backed Save + PageSetup!
        let mut app2 = ChiralCraftApp::new();
        let bl = app2.doc().settings.bond_length_pt;
        app2.doc_mut().stamp_molecule_tx(stamp_template(BuiltinTemplate::Benzene, 100.0, 100.0, bl), "M0");
        app2.doc_mut().stamp_molecule_tx(stamp_template(BuiltinTemplate::Benzene, 300.0, 100.0, bl), "M1");
        assert_eq!(app2.doc().scene.molecules.len(), 2);
        // Select ONLY molecule 1 and scale it — molecule 0 must remain untouched!
        let m0_before = app2.doc().scene.molecules[0].atoms[0].x;
        let m1_span_before = mol_bbox(&app2.doc().scene.molecules[1]).2 - mol_bbox(&app2.doc().scene.molecules[1]).0;
        app2.tab_mut().selected_mols = vec![1];
        app2.do_scale_targets(2.0);
        assert!((app2.doc().scene.molecules[0].atoms[0].x - m0_before).abs() < 1e-6);
        let m1_span_after = mol_bbox(&app2.doc().scene.molecules[1]).2 - mol_bbox(&app2.doc().scene.molecules[1]).0;
        assert!((m1_span_after - m1_span_before * 2.0).abs() < 1e-4);
        // Repeat scales molecule 1 by 2x again!
        app2.run_repeat_action();
        let m1_span_rep = mol_bbox(&app2.doc().scene.molecules[1]).2 - mol_bbox(&app2.doc().scene.molecules[1]).0;
        assert!((m1_span_rep - m1_span_before * 4.0).abs() < 1e-3);
        // Selection-aware Cut removes ONLY molecule 1!
        app2.do_cut_action();
        assert_eq!(app2.doc().scene.molecules.len(), 1);
        assert_eq!(app2.clipboard_mols.len(), 1);
        // Landscape Page Setup swaps effective width & height
        app2.tab_mut().landscape = true;
        let (ew, eh) = app2.tab().effective_page_size_pt();
        assert!(ew > eh);
        // Real disk-backed Save writes directly to file_path
        let save_p = "/tmp/chiralcraft_test_save.cdxml";
        app2.tab_mut().file_path = Some(save_p.to_string());
        app2.do_save_action();
        assert!(fs::read_to_string(save_p).unwrap().contains("<CDXML"));
        let _ = fs::remove_file(save_p);

        // Verify exact Nature Research & ACS 1996 binary .cds stylesheet parameter decoding!
        let (_, nr) = app2.stylesheets.iter().find(|(n, _)| n == "Nature Research Document").expect("Nature Research stylesheet");
        assert!((nr.bond_length_pt - 10.80).abs() < 0.05);
        assert!((nr.line_width_pt - 0.60).abs() < 0.05);
        assert!((nr.bold_width_pt - 1.56).abs() < 0.05);
        assert!((nr.label_size_pt - 6.0).abs() < 0.1);
        let (_, acs) = app2.stylesheets.iter().find(|(n, _)| n == "ACS Document 1996").expect("ACS stylesheet");
        assert!((acs.bond_length_pt - 14.40).abs() < 0.05);
        assert!((acs.line_width_pt - 0.60).abs() < 0.05);
        assert!((acs.bold_width_pt - 2.01).abs() < 0.05);
        assert!((acs.label_size_pt - 10.0).abs() < 0.1);

        // Verify IUPAC GR-1.6.2 180° linear alkyne sprouting!
        let mut app3 = ChiralCraftApp::new();
        app3.tool = Tool::TripleBond;
        app3.handle_canvas_drag_release(100.0, 100.0, 114.4, 100.0, false);
        app3.tool = Tool::SolidBond;
        app3.handle_canvas_click(114.4, 100.0, false);
        let a2 = &app3.doc().scene.molecules[0].atoms[2];
        assert!((a2.y - 100.0).abs() < 0.2 && a2.x > 125.0, "Alkyne sprout must be 180° linear: ({}, {})", a2.x, a2.y);

        // Verify ChiralCraft End-to-End Capabilities!
        let mut pro = ChiralCraftApp::new();
        assert!(pro.monomer_lib.entries.len() >= 20, "Official natural-monomers.json must load");
        // Name -> Structure -> Name + 1H/13C NMR
        pro.do_name_to_structure("Aspirin").unwrap();
        assert_eq!(pro.doc().scene.molecules.len(), 1);
        pro.do_structure_to_name();
        assert!(pro.doc().scene.captions.last().unwrap().text.contains("aspirin"));
        pro.do_predict_nmr(false);
        assert!(pro.show_nmr_window && pro.nmr_prediction.as_ref().unwrap().peaks.len() >= 5);
        pro.do_predict_nmr(true);
        assert_eq!(pro.nmr_prediction.as_ref().unwrap().nucleus, "13C");
        // Stoichiometry Grid, AutoNumber, Enhanced Stereo, Ring Fill, GEP
        pro.do_autonumber_reaction();
        assert_eq!(pro.tab().reaction_auto_numbers.len(), 1);
        pro.do_open_stoichiometry_grid();
        assert_eq!(pro.stoich_rows.len(), 1);
        pro.do_mark_enhanced_stereo(EnhancedStereoKind::Absolute);
        assert!(!pro.tab().enhanced_stereo.is_empty());
        assert!(!extract_ring_polygons(&pro.doc().scene.molecules[0]).is_empty());
    }

    #[test]
    fn test_drawn_and_imported_interactions_flawless() {
        use crate::chem_db::DbSource;
        let mut app = ChiralCraftApp::new();
        // 1. Import a curated PubChem-style molecule with label onto canvas
        let hit = DbHit {
            source: DbSource::PubChem,
            id: "CID 2244".to_string(),
            cid: Some(2244),
            name: "Aspirin".to_string(),
            iupac_name: Some("2-acetyloxybenzoic acid".to_string()),
            formula: "C9H8O4".to_string(),
            molecular_weight: "180.16".to_string(),
            smiles: "CC(=O)Oc1ccccc1C(=O)O".to_string(),
            inchikey: Some("BSYNRYMUTXBXSQ-UHFFFAOYSA-N".to_string()),
            max_phase: None,
            first_approval: None,
            alogp: None,
            psa: None,
            ro5_violations: None,
            similarity_pct: None,
            molfile: None,
            preview_mol: None,
            sdf_loading: false,
        };
        app.do_import_db_hit(&hit, true);
        assert_eq!(app.doc().scene.molecules.len(), 1);
        assert_eq!(app.tab().selected_mols, vec![0]);
        assert!(app.doc().scene.captions.len() >= 2, "Import + Label must create captions");

        // 2. Centered caption hit-test + Marquee caption drag!
        let cap0 = app.doc().scene.captions[0].clone();
        assert_eq!(cap0.justification, LabelJustification::Center);
        assert_eq!(app.hit_caption(cap0.x, cap0.y), Some(0));
        assert_eq!(app.hit_caption(cap0.x - 15.0, cap0.y), Some(0), "Left half of centered caption must hit");
        app.tool = Tool::Marquee;
        app.dragging_caption = Some(0);
        app.handle_canvas_drag_release(cap0.x, cap0.y, cap0.x + 25.0, cap0.y + 12.0, false);
        app.dragging_caption = None;
        assert!((app.doc().scene.captions[0].x - (cap0.x + 25.0)).abs() < 1e-6);

        // 3. Draw a second benzene ring next to imported Aspirin, then drag a bond ACROSS molecules to merge them!
        let bl = app.doc().settings.bond_length_pt;
        let a_imp = (app.doc().scene.molecules[0].atoms[0].x, app.doc().scene.molecules[0].atoms[0].y);
        app.doc_mut().stamp_molecule_tx(
            stamp_template(BuiltinTemplate::Benzene, a_imp.0 + 90.0, a_imp.1, bl),
            "Drawn Benzene",
        );
        assert_eq!(app.doc().scene.molecules.len(), 2);
        let total_atoms_before = app.doc().scene.molecules[0].atoms.len() + app.doc().scene.molecules[1].atoms.len();
        let b_drawn = (app.doc().scene.molecules[1].atoms[0].x, app.doc().scene.molecules[1].atoms[0].y);
        app.tool = Tool::SolidBond;
        app.drag_start_atom = app.hit_atom(a_imp.0, a_imp.1, 6.0);
        assert_eq!(app.drag_start_atom, Some((0, 0)));
        app.handle_canvas_drag_release(a_imp.0, a_imp.1, b_drawn.0, b_drawn.1, false);
        app.drag_start_atom = None;
        assert_eq!(app.doc().scene.molecules.len(), 1, "Cross-molecule bond drag must merge into 1 molecule");
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), total_atoms_before, "No duplicate floating atom created");

        // 4. Dragging a bond between two already-bonded atoms upgrades bond order instead of duplicating!
        let bonds_before = app.doc().scene.molecules[0].bonds.len();
        let b0 = app.doc().scene.molecules[0].bonds[0].clone();
        let p_a = (app.doc().scene.molecules[0].atoms[b0.begin_atom].x, app.doc().scene.molecules[0].atoms[b0.begin_atom].y);
        let p_b = (app.doc().scene.molecules[0].atoms[b0.end_atom].x, app.doc().scene.molecules[0].atoms[b0.end_atom].y);
        app.drag_start_atom = Some((0, b0.begin_atom));
        app.handle_canvas_drag_release(p_a.0, p_a.1, p_b.0, p_b.1, false);
        app.drag_start_atom = None;
        assert_eq!(app.doc().scene.molecules[0].bonds.len(), bonds_before, "Must NOT duplicate existing bond");
        assert_ne!(app.doc().scene.molecules[0].bonds[0].order, b0.order, "Existing bond order must cycle");

        // 5. Clicking empty space with Eraser must NOT push a no-op Undo step, and erasing caption works!
        app.tool = Tool::Eraser;
        let undo_before = app.doc().undo_menu_label();
        app.handle_canvas_click(5.0, 5.0, false);
        assert_eq!(app.doc().undo_menu_label(), undo_before, "Missed Eraser click must not pollute Undo stack");
        let caps_before = app.doc().scene.captions.len();
        let c_pos = (app.doc().scene.captions[0].x, app.doc().scene.captions[0].y);
        app.handle_canvas_click(c_pos.0, c_pos.1, false);
        assert_eq!(app.doc().scene.captions.len(), caps_before - 1);

        // 6. Undo + Redo roundtrip stays clean and sanitized!
        app.doc_mut().undo();
        app.sanitize_state();
        assert_eq!(app.doc().scene.captions.len(), caps_before);
        app.doc_mut().redo();
        app.sanitize_state();
        assert_eq!(app.doc().scene.captions.len(), caps_before - 1);
    }

    #[test]
    fn test_editable_svg_roundtrip_wordcraft_deckcraft_workflow() {
        // 1. Build a chiral structure with caption in ChiralCraft
        let mut app = ChiralCraftApp::new();
        app.do_name_to_structure("Aspirin").unwrap();
        app.do_structure_to_name();
        assert_eq!(app.doc().scene.molecules.len(), 1);
        assert_eq!(app.doc().scene.captions.len(), 1);
        let orig_atoms = app.doc().scene.molecules[0].atoms.len();
        let orig_bonds = app.doc().scene.molecules[0].bonds.len();

        // 2. Export Editable SVG (transparent for DeckCraft/WordCraft) & verify embedded metadata
        let svg_str = app.doc().selection_to_svg(&[], true);
        assert!(svg_str.contains("xmlns:chiralcraft=\"https://chiralcraft.dev/ns/1.0\""));
        assert!(svg_contains_chiralcraft_document(&svg_str));
        assert!(!svg_str.contains("fill=\"white\"/>\n<line"), "Transparent clipboard SVG omits outer white box");

        // 3. Save to temp .svg file (simulating WordCraft/DeckCraft double-click temp file)
        let svg_path = "/tmp/chiralcraft_roundtrip_test.svg";
        atomic_write_text(svg_path, &svg_str).unwrap();

        // 4. Launch ChiralCraft in `--edit-svg` round-trip mode, modify the structure, and flush on exit!
        let mut editor = ChiralCraftApp::new_with_file(svg_path, true);
        assert_eq!(editor.doc().scene.molecules.len(), 1);
        assert_eq!(editor.doc().scene.molecules[0].atoms.len(), orig_atoms);
        assert_eq!(editor.doc().scene.molecules[0].bonds.len(), orig_bonds);
        // User edits structure by adding a second molecule (Benzene)
        let bl = editor.doc().settings.bond_length_pt;
        editor.doc_mut().stamp_molecule_tx(stamp_template(BuiltinTemplate::Benzene, 340.0, 160.0, bl), "Add Benzene");
        // Closing ChiralCraft saves updated editable SVG atomically back to svg_path!
        editor.save_roundtrip_svg().unwrap();

        // 5. Host app (WordCraft / DeckCraft) reads updated SVG back -> contains BOTH molecules & valid SVG!
        let updated_svg = fs::read_to_string(svg_path).unwrap();
        assert!(updated_svg.starts_with("<svg ") && updated_svg.contains("</svg>"));
        let reloaded = load_document_any(svg_path).unwrap();
        assert_eq!(reloaded.scene.molecules.len(), 2, "Edited SVG must preserve added structure on round-trip");
        let _ = fs::remove_file(svg_path);
    }

    #[test]
    fn test_full_stress_templates_custom_labels_transparent_svg_and_empty_canvas_safety() {
        // 1. Every single clean-room built-in template across all 11 libraries must have valid non-zero 2D geometry!
        let libs = TemplateLibrary::builtin_libraries();
        assert_eq!(libs.len(), 11);
        for lib in &libs {
            assert!(lib.total_count() >= 5, "Library {} should have >=5 templates", lib.name);
            for (ti, mol) in lib.templates.iter().enumerate() {
                assert!(mol.atoms.len() >= 3, "{} #{} atoms", lib.name, ti);
                for (bi, b) in mol.bonds.iter().enumerate() {
                    let a1 = &mol.atoms[b.begin_atom];
                    let a2 = &mol.atoms[b.end_atom];
                    let d = (a2.x - a1.x).hypot(a2.y - a1.y);
                    assert!(
                        d > 1.0,
                        "Zero/collapsed bond in built-in template '{}' #{} bond {}: len={d:.3}",
                        lib.name, ti, bi
                    );
                }
            }
        }

        // 2. Full lossless Editable SVG round-trip: Custom atom label ("Boc") + Wedged stereo + Reaction Arrow + Centered Caption + Nature Research 10.8pt style + ZERO white sticker circles!
        let mut app = ChiralCraftApp::new();
        let (_, nr) = app.stylesheets.iter().find(|(n, _)| n == "Nature Research Document").unwrap().clone();
        app.doc_mut().settings = nr;
        let mut m = stamp_template(BuiltinTemplate::Benzene, 150.0, 150.0, 10.8);
        m.atoms[0].label = Some("Boc".to_string());
        m.bonds[0].display = BondDisplay::WedgeBegin;
        app.doc_mut().stamp_molecule_tx(m, "Labeled Chiral Ring");
        app.doc_mut().apply_scene_tx("Add Arrow & Caption", |s, nid| {
            let aid = *nid; *nid += 1;
            s.arrows.push(SceneArrow {
                id: aid, tail: (200.0, 150.0), head: (260.0, 150.0),
                arrowhead_head: "Full".to_string(), arrowhead_tail: None, no_go: None, dipole: false,
            });
            let cid = *nid; *nid += 1;
            s.captions.push(CaptionNode {
                id: cid, text: "Boc-Protected Core".to_string(), x: 150.0, y: 195.0,
                element: None, generic_nickname: None, warning: None, justification: LabelJustification::Center,
            });
        });
        let tsvg = app.doc().to_svg_transparent();
        assert!(!tsvg.contains("fill=\"white\""), "Transparent editable SVG must NOT contain white background/sticker shapes!");
        let roundtripped = InteractiveDocument::from_svg(&tsvg).expect("from_svg roundtrip");
        assert!((roundtripped.settings.bond_length_pt - 10.8).abs() < 0.05);
        assert_eq!(roundtripped.scene.molecules.len(), 1);
        assert_eq!(roundtripped.scene.molecules[0].atoms[0].label.as_deref(), Some("Boc"));
        assert_eq!(roundtripped.scene.molecules[0].bonds[0].display, BondDisplay::WedgeBegin);
        assert_eq!(roundtripped.scene.arrows.len(), 1);
        assert_eq!(roundtripped.scene.captions.len(), 1);
        assert_eq!(roundtripped.scene.captions[0].text, "Boc-Protected Core");

        // 3. Empty-canvas stress test: every action on an empty document must succeed gracefully with zero panics!
        let mut empty = ChiralCraftApp::new();
        empty.do_copy_action();
        empty.do_cut_action();
        empty.do_paste_action();
        empty.do_clear_action();
        empty.do_clean_up_structure();
        empty.do_flatten_action();
        empty.do_toggle_hydrogens();
        empty.do_scale_targets(1.5);
        empty.do_rotate_targets(45.0);
        empty.do_center_on_page_action();
        empty.do_structure_to_name();
        empty.do_predict_nmr(false);
        empty.do_predict_nmr(true);
        empty.do_autonumber_reaction();
        empty.do_open_stoichiometry_grid();
        assert_eq!(empty.doc().scene.molecules.len(), 0);

        // 4. Embedded 128x128 ChiralCraft app icon decodes cleanly to RGBA!
        let icon = load_app_icon_data().expect("embedded app icon must decode");
        assert_eq!((icon.width, icon.height), (128, 128));
        assert_eq!(icon.rgba.len(), 128 * 128 * 4);
    }

    #[test]
    fn test_reshiki_inspired_features_end_to_end() {
        // 1. Catalytic Cycle Builder (starter + elliptical arrange)
        let mut app = ChiralCraftApp::new();
        app.doc_mut().apply_scene_tx("Starter Cycle", |s, nid| {
            let idxs = insert_starter_catalytic_cycle(s, nid, 250.0, 210.0);
            assert_eq!(idxs.len(), 4);
        });
        assert_eq!(app.doc().scene.molecules.len(), 4);
        assert_eq!(app.doc().scene.arrows.len(), 4);
        app.arrange_selected_catalytic_cycle(true);
        assert_eq!(app.doc().scene.arrows.len(), 4, "Re-arranging cycle must be idempotent");

        // 2. Live Linked Property & Name Captions auto-update when molecule changes!
        let mut app2 = ChiralCraftApp::new();
        let blen = app2.doc().settings.bond_length_pt;
        app2.doc_mut().stamp_molecule_tx(stamp_template(BuiltinTemplate::Benzene, 180.0, 160.0, blen), "Benzene");
        app2.attach_linked_caption(LinkedCaptionKind::Formula);
        app2.attach_linked_caption(LinkedCaptionKind::IupacName);
        assert_eq!(app2.doc().scene.captions.len(), 2);
        assert!(app2.doc().scene.captions[0].text.contains("C6H6"));
        // Mutate an atom C -> N (pyridine) and verify sync_linked_captions updates automatically!
        app2.doc_mut().scene.molecules[0].atoms[0].z = 7;
        app2.sync_linked_captions();
        assert!(app2.doc().scene.captions[0].text.contains("C5H5N"), "Linked formula must auto-update: {}", app2.doc().scene.captions[0].text);
        assert!(app2.doc().scene.captions[1].text.to_lowercase().contains("pyridine"), "Linked name must auto-update: {}", app2.doc().scene.captions[1].text);

        // 3. Haworth sugars, Porphine, 18-Crown-6, Adamantane, Norbornane + Crossing-Bond Clearance
        let haworth = build_haworth_pyranose(blen);
        assert!(haworth.atoms.len() >= 10);
        assert!(haworth.bonds.iter().any(|b| b.display == BondDisplay::Bold));
        let porphine = build_porphine_core(blen);
        assert_eq!(porphine.atoms.len(), 24);
        let crown = build_crown_ether_18_6(blen);
        assert_eq!(crown.atoms.len(), 18);
        let adamantane = build_adamantane_2d(blen);
        assert_eq!(adamantane.atoms.len(), 10);
        let norbornane = build_norbornane_2d(blen);
        assert_eq!(norbornane.atoms.len(), 7);
        let crossings = find_bond_crossings(&adamantane);
        assert!(!crossings.is_empty(), "Adamantane 2D cage must detect crossing bonds for bridge clearance");

        // 4. Lab Glassware & Haworth libraries present in built-in 10 libraries
        assert!(app2.ctp_libraries.iter().any(|l| l.name == "Haworth & Macrocycles"));
        let glass = app2.ctp_libraries.iter().find(|l| l.name == "Lab Glassware & Apparatus").expect("Lab Glassware library");
        assert!(glass.graphic_templates.len() >= 6);

        // 5. 3D .xyz & .cif loading + 3D Depth Cueing + Dark-Slide white-ink SVG export
        let xyz_path = "/tmp/chiralcraft_test_3d.xyz";
        let _ = fs::write(
            xyz_path,
            "5\nMethane 3D\nC 0.0 0.0 0.0\nH 0.63 0.63 0.63\nH -0.63 -0.63 0.63\nH -0.63 0.63 -0.63\nH 0.63 -0.63 -0.63\n",
        );
        let doc3d = load_document_any(xyz_path).expect("load .xyz");
        assert_eq!(doc3d.scene.molecules.len(), 1);
        assert!(doc3d.scene.molecules[0].has_3d_coords);
        let ops = compute_3d_atom_opacities(&doc3d.scene.molecules[0], 0.35);
        assert!(ops.iter().any(|&o| o < 0.5) && ops.iter().any(|&o| o > 0.9));
        let dark_svg = doc3d.to_svg_styled(true, true, 0.35, 3.6);
        assert!(dark_svg.contains("#F2F5F9"), "Dark-slide SVG must render crisp light ink");
        assert!(dark_svg.contains("stroke-opacity="), "3D SVG must emit depth-cueing stroke-opacity");
        let _ = fs::remove_file(xyz_path);

        // 6. Editable TLC Plate & Table SVG rendering
        let mut tlc_doc = InteractiveDocument::new_saved_seed("TLC", DocumentSettings::acs_document_1996());
        tlc_doc.scene.tables.push(SceneTable {
            id: 1,
            rows: 2,
            cols: 2,
            bbox: (20.0, 20.0, 120.0, 70.0),
            cells: vec![vec!["Cat".into(), "Yield".into()], vec!["Pd(0)".into(), "96%".into()]],
        });
        tlc_doc.scene.tlc_plates.push(SceneTlcPlate {
            id: 2,
            bbox: (140.0, 20.0, 185.0, 100.0),
            origin_fraction: 0.15,
            spot_rfs: vec![0.3, 0.7],
            spot_labels: vec!["SM".into(), "Prod".into()],
            spot_colors: vec![(20, 20, 20), (205, 45, 45)],
            solvent_label: "EtOAc/Hex".into(),
        });
        let ext_svg = tlc_doc.to_svg();
        assert!(ext_svg.contains("Pd(0)") && ext_svg.contains("96%"));
        assert!(ext_svg.contains("EtOAc/Hex") && ext_svg.contains("rgb(205,45,45)"));
    }

    #[test]
    fn test_reshiki_wave2_all_7_features() {
        // 1. Editable PNG & Vector PDF lossless round-trip!
        let mut doc = InteractiveDocument::new_saved_seed("RoundTripMedia", DocumentSettings::acs_document_1996());
        let mut benz = stamp_template(BuiltinTemplate::Benzene, 160.0, 160.0, 14.40);
        benz.atoms[0].z = 7; // Pyridine
        doc.stamp_molecule_tx(benz, "Pyridine");

        let png_bytes = doc.to_editable_png(300).expect("editable png");
        let from_png = InteractiveDocument::from_png(&png_bytes).expect("reopen from png");
        assert_eq!(from_png.scene.molecules.len(), 1);
        assert_eq!(from_png.scene.molecules[0].atoms[0].z, 7);

        let pdf_bytes = doc.to_editable_pdf(false, 20.0);
        assert!(pdf_bytes.starts_with(b"%PDF-1.4"));
        let from_pdf = InteractiveDocument::from_pdf(&pdf_bytes).expect("reopen from pdf");
        assert_eq!(from_pdf.scene.molecules.len(), 1);
        assert_eq!(from_pdf.scene.molecules[0].atoms[0].z, 7);

        // 2. Chemistry-Aware Document Diff
        let mut doc_b = doc.clone();
        doc_b.scene.molecules[0].atoms[0].z = 6; // Back to benzene!
        let diff = compare_documents(&doc, &doc_b);
        assert!(diff.has_changes());
        let rep = diff.format_text_report();
        assert!(rep.contains("C5H5N") && rep.contains("C6H6"), "Diff report: {rep}");

        // 3. Smart Reverse-Order Group Typing (HO, H2N, O2N, MeO, AcO, F3C, BocHN)
        assert_eq!(resolve_reverse_abbreviation("O2N"), Some(("NO2", "O₂N")));
        assert_eq!(resolve_reverse_abbreviation("MeO"), Some(("OMe", "MeO")));
        assert_eq!(resolve_reverse_abbreviation("BocHN"), Some(("NHBoc", "BocHN")));
        let mut app = ChiralCraftApp::new();
        app.doc_mut().stamp_molecule_tx(stamp_template(BuiltinTemplate::Benzene, 180.0, 150.0, 14.40), "Benzene");
        app.inline_atom_edit = Some((0, 0));
        app.inline_label_buf = "HO".to_string();
        app.commit_inline_atom_edit();
        assert_eq!(app.doc().scene.molecules[0].atoms[0].z, 8);
        assert_eq!(app.doc().scene.molecules[0].atoms[0].label, None);
        app.inline_atom_edit = Some((0, 1));
        app.inline_label_buf = "O2N".to_string();
        app.commit_inline_atom_edit();
        assert_eq!(app.doc().scene.molecules[0].atoms[1].z, 7);
        assert_eq!(app.doc().scene.molecules[0].atoms[1].label.as_deref(), Some("O₂N"));

        // 4. On-Canvas NMR Shift Stamping & Clear
        let mut app_nmr = ChiralCraftApp::new();
        app_nmr.doc_mut().stamp_molecule_tx(smiles_to_molecule("CCO").unwrap(), "Ethanol");
        app_nmr.stamp_nmr_shift_labels(false);
        assert!(app_nmr.doc().scene.captions.iter().any(|c| c.text.starts_with('δ')));
        app_nmr.clear_nmr_shift_labels();
        assert!(!app_nmr.doc().scene.captions.iter().any(|c| c.text.starts_with('δ')));

        // 5. 2D -> 3D Conformer Lift + ExportOptions
        let lifted = lift_molecule_2d_to_3d(&smiles_to_molecule("C1CCCCC1").unwrap());
        assert!(lifted.has_3d_coords && lifted.atoms.iter().any(|a| a.z_coord.abs() > 1.0));
        assert!((ExportOptions { png_dpi: 600, ..Default::default() }.png_scale() - (600.0 / 72.0)).abs() < 1e-6);

        // 6. Command Palette execution (rotate + dpi + special core)
        app_nmr.execute_palette_command(&egui::Context::default(), PaletteCommandAction::SetExportDpi(600));
        assert_eq!(app_nmr.export_opts.png_dpi, 600);
        app_nmr.execute_palette_command(&egui::Context::default(), PaletteCommandAction::InsertSpecialCore("porphine"));
        assert_eq!(app_nmr.doc().scene.molecules.len(), 2);
    }

    #[test]
    fn test_paid_suite_all_7_flagship_features() {
        // 1. MedChem Scorecard + PAINS alert detection (4-nitroaniline)
        let nitro_an = smiles_to_molecule("Nc1ccc([N+](=O)[O-])cc1").unwrap();
        let mc = compute_medchem_report(&nitro_an);
        assert!(mc.lipinski_pass && mc.veber_pass);
        assert!(mc.alerts.iter().any(|a| a.name.contains("Nitroaromatic")));
        assert!(mc.alerts.iter().any(|a| a.name.contains("Aniline")));

        // 2. HRMS Isotope Envelope (Chlorobenzene -> 3:1 M / M+2 Cl doublet!)
        let clb = smiles_to_molecule("Clc1ccccc1").unwrap();
        let hrms = simulate_hrms(&clb, "[M+H]+");
        assert!(hrms.peaks.len() >= 2);
        let m2 = hrms.peaks.iter().find(|p| (p.mz - (hrms.peaks[0].mz + 1.997)).abs() < 0.03).expect("Cl M+2");
        assert!(m2.rel_intensity > 28.0 && m2.rel_intensity < 40.0);

        // 3. pH Microspecies, Tautomers, Stereoisomers, Resonance & Electron Pushing
        let gly = smiles_to_molecule("NCC(=O)O").unwrap();
        let phys = chiralcraft_chem::compute_microspecies_at_ph(&gly, 7.4);
        assert_eq!(phys.formal_charge_sum, 0);
        assert!(phys.isoelectric_point.is_some());

        // 4. Markush R-Group Enumeration + SAR Export
        let mut core = smiles_to_molecule("c1ccccc1C").unwrap();
        core.atoms[6].label = Some("R1".into());
        let analogs = chiralcraft_doc::enumerate_rgroup_library(&core, "H, F, OMe, CF3", "", 14.40);
        assert_eq!(analogs.len(), 4);

        // 5. Retrosynthesis Advisor + Stamping onto App Scene
        let mut app = ChiralCraftApp::new();
        let acet = smiles_to_molecule("CC(=O)Nc1ccccc1").unwrap();
        app.doc_mut().stamp_molecule_tx(acet.clone(), "Acetanilide");
        let hits = suggest_disconnections(&acet, 14.40);
        assert!(!hits.is_empty());
        app.handle_pro_suite_action(ProSuiteAction::StampRetrosynthesis(hits[0].clone()));
        assert!(app.doc().scene.arrows.iter().any(|a| a.arrowhead_head == "Retrosynthetic"));

        // 6. Publication TOC Frame + Arrow Conditions + Serpentine Scheme
        app.handle_pro_suite_action(ProSuiteAction::InsertTocFrame(chiralcraft_doc::TocPreset::presets()[0]));
        assert!(app.pro_suite.last_toc_bbox.is_some());
        app.handle_pro_suite_action(ProSuiteAction::StampArrowConditions {
            above: "HATU, DIPEA".into(),
            below: "DMF, rt, 94%".into(),
        });
        assert!(app.doc().scene.captions.iter().any(|c| c.text.contains("HATU")));
        app.handle_pro_suite_action(ProSuiteAction::WrapSerpentine(2));
        assert!(app.doc().scene.arrows.len() >= 3);

        // 7. JCAMP-DX Experimental NMR Overlay + MCP flagship tools
        let jdx = "##TITLE=Demo 1H\n##FIRSTX=8.0\n##LASTX=1.0\n##NPOINTS=4\n##XYDATA=(X++(Y..Y))\n8.0 0.1 0.9 0.2 0.0\n##END=\n";
        let spec = parse_jcamp_dx(jdx).expect("parse jcamp");
        assert_eq!(spec.points.len(), 4);
        let mcp_mc = mcp::handle_mcp_request(r#"{"jsonrpc":"2.0","id":10,"method":"tools/call","params":{"name":"medchem_scorecard","arguments":{"smiles":"CC(=O)Oc1ccccc1C(=O)O"}}}"#).unwrap();
        assert!(mcp_mc.contains("CNS MPO") && mcp_mc.contains("QED"));
        let mcp_hrms = mcp::handle_mcp_request(r#"{"jsonrpc":"2.0","id":11,"method":"tools/call","params":{"name":"simulate_hrms","arguments":{"smiles":"Brc1ccccc1"}}}"#).unwrap();
        assert!(mcp_hrms.contains("C6H5Br"));
        let mcp_retro = mcp::handle_mcp_request(r#"{"jsonrpc":"2.0","id":12,"method":"tools/call","params":{"name":"retrosynthesis_suggest","arguments":{"smiles":"CC(=O)Nc1ccccc1"}}}"#).unwrap();
        assert!(mcp_retro.contains("Amide Coupling"));
    }

    #[test]
    fn test_omni_suite_all_8_discipline_features() {
        // 1. Multi-Spectroscopy (FT-IR, UV-Vis, 2D NMR)
        let asp = smiles_to_molecule("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let ir = predict_ir_spectrum(&asp);
        assert!(ir.bands.iter().any(|b| b.assignment.contains("C=O")));
        let uv = predict_uv_vis(&asp);
        assert!(uv.lambda_max_nm > 240.0);
        assert!(!chiralcraft_chem::predict_2d_nmr(&asp, chiralcraft_chem::Nmr2dKind::Hsqc).peaks.is_empty());

        // 2. Organometallic 18e- & Built-in Catalyst Library (#11)
        let mut app = ChiralCraftApp::new();
        let cat_lib = app.ctp_libraries.iter().find(|l| l.name == "Organometallic Catalysts").expect("Catalysts lib");
        assert_eq!(cat_lib.templates.len(), 6);
        let pd_rep = analyze_organometallic_complexes(&cat_lib.templates[0]);
        assert_eq!(pd_rep[0].total_valence_electrons, 18);

        // 3. Reaction Energy Profile on Canvas & Exported SVG + Eyring/Boltzmann
        let prof = chiralcraft_doc::EnergyProfile::starter(1, 50.0, 50.0);
        app.handle_omni_suite_action(&egui::Context::default(), OmniSuiteAction::StampEnergyProfile(prof));
        assert_eq!(app.doc().scene.energy_profiles.len(), 1);
        let ep_svg = app.doc().to_svg();
        assert!(ep_svg.contains("TS1") && ep_svg.contains("kcal/mol"));

        // 4. 3D Force-Field UFF Optimization
        app.doc_mut().stamp_molecule_tx(smiles_to_molecule("CCC").unwrap(), "Propane");
        app.execute_palette_command(&egui::Context::default(), PaletteCommandAction::Optimize3dNow);
        assert!(app.doc().scene.molecules[0].has_3d_coords);

        // 5. LaTeX \chemfig & TikZ + ELN SI + Virtual Screening + MCP 15 tools
        let tex = document_to_chemfig_scheme(app.doc());
        assert!(tex.contains("\\schemestart") && tex.contains("\\chemfig{"));
        let demo_screen = chiralcraft_doc::builtin_demo_screening_library(14.40);
        assert_eq!(demo_screen.len(), 10);

        let mcp_spec = mcp::handle_mcp_request(r#"{"jsonrpc":"2.0","id":20,"method":"tools/call","params":{"name":"predict_spectroscopy","arguments":{"smiles":"CC(=O)Oc1ccccc1C(=O)O"}}}"#).unwrap();
        assert!(mcp_spec.contains("IR (neat") && mcp_spec.contains("HSQC"));
        let mcp_tex = mcp::handle_mcp_request(r#"{"jsonrpc":"2.0","id":21,"method":"tools/call","params":{"name":"export_latex","arguments":{"smiles":"c1ccccc1O"}}}"#).unwrap();
        assert!(mcp_tex.contains("chemfig") && mcp_tex.contains("tikzpicture"));
        let mcp_si = mcp::handle_mcp_request(r#"{"jsonrpc":"2.0","id":22,"method":"tools/call","params":{"name":"generate_si_procedure","arguments":{"smiles":"CCO","name":"Ethanol","solvent":"THF","yield_pct":95}}}"#).unwrap();
        assert!(mcp_si.contains("HRMS (ESI+)") && mcp_si.contains("95% yield"));
    }

    #[test]
    fn test_verification_audit_regressions() {
        // 1. Reversed label expansion (O₂N -> NO2, MeO -> OMe, F₃C -> CF3, BocHN -> NHBoc)
        let mut app = ChiralCraftApp::new();
        let blen = app.doc().settings.bond_length_pt;
        let mut ph = stamp_template(BuiltinTemplate::Benzene, 160.0, 150.0, blen);
        let mut sub = ChemAtom::new(99, 6);
        sub.x = ph.atoms[0].x - blen;
        sub.y = ph.atoms[0].y;
        sub.label = Some("O₂N".to_string());
        let sub_i = ph.add_atom(sub);
        ph.add_bond(99, 0, sub_i, 1, BondDisplay::Solid);
        app.doc_mut().stamp_molecule_tx(ph, "Nitrobenzene");
        assert!(app.doc_mut().expand_label_tx(0, sub_i, "O₂N"), "Reversed Unicode label O₂N must expand");
        assert!(app.doc().scene.molecules[0].atoms.len() > 7);
        assert!(app.doc_mut().undo(), "Expand Label must undo cleanly");
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), 7);

        // 2. Dark-slide SVG ink consistency (aromatic circle, arrow, table, TLC, energy profile) + XML escaping in header/footer
        app.toggle_aromatic_display_tx();
        app.doc_mut().apply_scene_tx("Add Dark-Slide Objects", |s, next_id| {
            let aid = *next_id; *next_id += 1;
            s.arrows.push(SceneArrow {
                id: aid,
                tail: (200.0, 150.0),
                head: (260.0, 150.0),
                arrowhead_head: "Full".to_string(),
                arrowhead_tail: None,
                no_go: None,
                dipole: false,
            });
            s.header_footer.header = Some("A & B < 50 °C > 0".to_string());
            s.energy_profiles.push(chiralcraft_doc::EnergyProfile::starter(500, 40.0, 220.0));
        });
        let dark_svg = app.doc().to_svg_styled(true, true, 0.45, 3.6);
        assert!(!dark_svg.contains("stroke=\"black\""), "Dark-slide SVG must not emit black strokes");
        assert!(!dark_svg.contains("fill=\"black\""), "Dark-slide SVG must not emit black fills");
        assert!(dark_svg.contains("A &amp; B &lt; 50 °C &gt; 0"), "Header must be XML-escaped");
        // Round-trip SVG must preserve energy profile & molecules
        let rt = InteractiveDocument::from_svg(&dark_svg).expect("Dark SVG must round-trip cleanly");
        assert_eq!(rt.scene.molecules.len(), 1);
        assert_eq!(rt.scene.energy_profiles.len(), 1);

        // 3. Disconnected salt + trivalent charge in LaTeX \chemfig
        let salt = smiles_to_molecule("[Fe+3].[Cl-].[Cl-].[Cl-]").unwrap();
        let tex = chiralcraft_doc::molecule_to_chemfig(&salt);
        assert!(tex.contains("Fe^{3+}") && tex.contains("\\+"), "Disconnected salt & +3 charge must serialize cleanly: {tex}");

        // 4. load_document_or_smiles accepts raw SMILES directly
        let smi_doc = load_document_or_smiles("CC(=O)Oc1ccccc1C(=O)O").expect("Raw SMILES must load");
        assert_eq!(smi_doc.scene.molecules.len(), 1);
        assert_eq!(chiralcraft_chem::compute_analysis_report(&smi_doc.scene.molecules[0]).formula, "C9H8O4");
    }

    #[test]
    fn test_exhaustive_gui_every_tool_window_dialog_and_action() {
        let mut app = ChiralCraftApp::new();
        let ctx = egui::Context::default();
        let run_frame = |app: &mut ChiralCraftApp| {
            let mut input = egui::RawInput::default();
            input.screen_rect = Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1440.0, 900.0)));
            let out = ctx.run(input, |c| app.render_ui_frame(c));
            assert!(!out.shapes.is_empty(), "Headless egui frame must produce rendered shapes");
        };

        // ---- 1. Render ALL floating windows & studios on an EMPTY canvas across all 3 themes & all 4 dock tabs ----
        app.show_analysis = true;
        app.show_properties = true;
        app.show_periodic_table = true;
        app.show_info_window = true;
        app.show_char_map = true;
        app.show_nmr_window = true; // Must NOT panic or vanish even with nmr_prediction == None!
        app.show_helm_toolbar = true;
        app.show_biographics_toolbar = true;
        app.show_ring_fill_toolbar = true;
        app.show_stoichiometry = true;
        app.show_object_inspector = true;
        app.show_object_toolbar = true;
        app.show_structure_toolbar = true;
        app.show_rulers = true;
        app.show_crosshair = true;
        app.show_reaction_interp = true;
        app.show_diff_window = true;
        app.show_chem_db = true;
        app.show_right_dock = true;
        app.viewer_3d.open = true;
        app.pro_suite.show_medchem = true;
        app.pro_suite.show_hrms = true;
        app.pro_suite.show_isomers = true;
        app.pro_suite.show_markush = true;
        app.pro_suite.show_retro = true;
        app.pro_suite.show_pub_composer = true;
        app.omni_suite.show_spectroscopy = true;
        app.omni_suite.show_ino_poly = true;
        app.omni_suite.show_energy_profile = true;
        app.omni_suite.show_eln_safety = true;
        app.omni_suite.show_latex_export = true;
        app.omni_suite.show_screening = true;
        app.quantum.show_window = true;
        app.cmd_palette.open_fresh();

        for theme in [UiTheme::ClassicLight, UiTheme::DarkChrome, UiTheme::DarkSlide] {
            app.ui_theme = theme;
            ChiralCraftApp::apply_theme(&ctx, theme);
            for dtab in [DockTab::Inspector, DockTab::MedChem, DockTab::Spectra, DockTab::Quantum] {
                app.dock_tab = dtab;
                run_frame(&mut app);
            }
        }
        app.cmd_palette.open = false;

        // ---- 2. Exercise EVERY SINGLE ONE of the 53 Tool::* variants (single click + drag + erase) ----
        let bond_tools = [
            Tool::SolidBond, Tool::DoubleBond, Tool::TripleBond, Tool::MultipleBond,
            Tool::DativeBond, Tool::QuadrupleBond, Tool::DashedBond, Tool::HashedBond,
            Tool::HashWedgeBond, Tool::BoldBond, Tool::WedgeBond, Tool::HollowWedgeBond, Tool::WavyBond,
        ];
        for (i, &t) in bond_tools.iter().enumerate() {
            let before = app.doc().scene.molecules.len();
            app.tool = t;
            app.handle_canvas_click(60.0 + (i as f64) * 28.0, 60.0, false);
            assert_eq!(app.doc().scene.molecules.len(), before + 1, "Bond tool {:?} single click must stamp bond", t);
        }
        // Drag a bond from an existing atom to empty space
        app.tool = Tool::WedgeBond;
        app.drag_start_atom = Some((0, 0));
        let a0 = (app.doc().scene.molecules[0].atoms[0].x, app.doc().scene.molecules[0].atoms[0].y);
        app.handle_canvas_drag_release(a0.0, a0.1, a0.0 + 20.0, a0.1 + 15.0, false);
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), 3);

        let ring_tools = [
            Tool::Cyclopropane, Tool::Cyclobutane, Tool::Cyclopentane, Tool::Cyclohexane,
            Tool::Cycloheptane, Tool::Cyclooctane, Tool::Chair1, Tool::Chair2,
            Tool::Cyclopentadiene, Tool::Benzene,
        ];
        for (i, &t) in ring_tools.iter().enumerate() {
            let before = app.doc().scene.molecules.len();
            app.tool = t;
            app.handle_canvas_click(70.0 + (i as f64) * 45.0, 140.0, false);
            assert_eq!(app.doc().scene.molecules.len(), before + 1, "Ring tool {:?} must stamp ring", t);
        }

        // Single-click default stamping for ALL arrow, shape, bracket, orbital, symbol, curve, chain, table, TLC, template tools!
        let arrow_tools = [Tool::Arrow, Tool::EquilibriumArrow, Tool::NoGoArrow, Tool::RetrosyntheticArrow, Tool::CurvedArrow];
        for (i, &t) in arrow_tools.iter().enumerate() {
            let before = app.doc().scene.arrows.len();
            app.tool = t;
            app.handle_canvas_click(80.0 + (i as f64) * 65.0, 210.0, false);
            assert_eq!(app.doc().scene.arrows.len(), before + 1, "Arrow tool {:?} single-click must stamp arrow", t);
        }

        let graphic_tools = [
            Tool::Orbital, Tool::OrbitalS, Tool::OrbitalD,
            Tool::ShapeRect, Tool::ShapeOval, Tool::ShapeLine,
            Tool::Bracket, Tool::RoundBracket, Tool::CurlyBracket,
            Tool::ChemSymbol, Tool::SymbolMinus, Tool::SymbolRadical, Tool::SymbolLonePair,
        ];
        for (i, &t) in graphic_tools.iter().enumerate() {
            let before = app.doc().scene.graphics.len();
            app.tool = t;
            app.handle_canvas_click(70.0 + (i as f64) * 35.0, 270.0, false);
            assert_eq!(app.doc().scene.graphics.len(), before + 1, "Graphic tool {:?} single-click must stamp object", t);
        }

        app.tool = Tool::PenCurve;
        let c_before = app.doc().scene.curves.len();
        app.handle_canvas_click(120.0, 320.0, false);
        assert_eq!(app.doc().scene.curves.len(), c_before + 1, "PenCurve single-click must add curve");

        app.tool = Tool::AcyclicChain;
        let m_before = app.doc().scene.molecules.len();
        app.handle_canvas_click(220.0, 320.0, false);
        assert_eq!(app.doc().scene.molecules.len(), m_before + 1, "AcyclicChain single-click must add chain");

        app.tool = Tool::Table;
        app.handle_canvas_click(320.0, 320.0, false);
        assert_eq!(app.doc().scene.tables.len(), 1);

        app.tool = Tool::TlcPlate;
        app.handle_canvas_click(430.0, 320.0, false);
        assert_eq!(app.doc().scene.tlc_plates.len(), 1);

        app.tool = Tool::TextLabel;
        app.handle_canvas_click(490.0, 320.0, false);
        assert!(matches!(app.dialog, ActiveDialog::NewCaption { .. }));
        app.dialog = ActiveDialog::None;

        // Tool::Templates across ALL 11 built-in template libraries!
        assert_eq!(app.ctp_libraries.len(), 11);
        for lib_i in 0..app.ctp_libraries.len() {
            app.selected_ctp_idx = lib_i;
            app.tool = Tool::Templates;
            let bm = app.doc().scene.molecules.len();
            let bc = app.doc().scene.curves.len();
            app.handle_canvas_click(150.0, 420.0, false);
            assert!(
                app.doc().scene.molecules.len() > bm || app.doc().scene.curves.len() > bc,
                "Template library #{} ({}) must stamp via Tool::Templates",
                lib_i, app.ctp_libraries[lib_i].name
            );
        }

        // Perspective, Fragmentation, ReactionMap, Marquee, Lasso
        app.tool = Tool::Marquee;
        app.handle_canvas_drag_release(0.0, 0.0, 500.0, 500.0, false);
        assert!(!app.tab().selected_mols.is_empty());
        app.tool = Tool::Perspective;
        app.handle_canvas_drag_release(200.0, 200.0, 225.0, 210.0, false);
        app.tool = Tool::ReactionMap;
        app.drag_start_atom = Some((0, 0));
        let b_pos = (app.doc().scene.molecules[1].atoms[0].x, app.doc().scene.molecules[1].atoms[0].y);
        app.handle_canvas_drag_release(a0.0, a0.1, b_pos.0, b_pos.1, false);
        assert_eq!(app.tab().reaction_atom_maps.len(), 2);

        // Energy Profile + Gel Plate + Eraser + Undo + Clear!
        app.handle_omni_suite_action(&ctx, OmniSuiteAction::StampEnergyProfile(chiralcraft_doc::EnergyProfile::starter(901, 50.0, 450.0)));
        app.handle_pro_panel_action(&ctx, ProPanelAction::InsertGelPlate);
        assert_eq!(app.doc().scene.energy_profiles.len(), 1);
        assert_eq!(app.tab().gel_plates.len(), 1);
        app.tool = Tool::Eraser;
        app.handle_canvas_click(90.0, 480.0, false); // Erases EnergyProfile
        assert_eq!(app.doc().scene.energy_profiles.len(), 0, "Eraser must erase EnergyProfile");
        assert!(app.doc_mut().undo(), "Erase EnergyProfile must undo cleanly");
        assert_eq!(app.doc().scene.energy_profiles.len(), 1);

        // ---- 3. Render POPULATED canvas with ALL 15 flyouts, ALL inner studio tabs, and ALL 40 ActiveDialog variants ----
        app.doc_mut().stamp_molecule_tx(smiles_to_molecule("CC(=O)Oc1ccccc1C(=O)O").unwrap(), "Aspirin");
        app.tab_mut().selected_mols = vec![app.doc().scene.molecules.len() - 1];
        app.do_predict_nmr(false);

        for fly in [
            ActiveFlyout::None, ActiveFlyout::Templates, ActiveFlyout::Arrows, ActiveFlyout::Brackets,
            ActiveFlyout::Chains, ActiveFlyout::ChemSymbols, ActiveFlyout::Chromatography,
            ActiveFlyout::Colors, ActiveFlyout::Curves, ActiveFlyout::DottedDashed,
            ActiveFlyout::Shapes, ActiveFlyout::Fragmentation, ActiveFlyout::MultipleBonds,
            ActiveFlyout::Orbitals, ActiveFlyout::PenTools, ActiveFlyout::Windows, ActiveFlyout::Accessories,
        ] {
            app.active_flyout = fly;
            run_frame(&mut app);
        }
        app.active_flyout = ActiveFlyout::None;

        for itab in [IsomerTab::PhMicrospecies, IsomerTab::Tautomers, IsomerTab::Stereoisomers, IsomerTab::Resonance, IsomerTab::ElectronPush] {
            app.pro_suite.isomer_tab = itab;
            run_frame(&mut app);
        }
        for stab in [SpecTab::FtIr, SpecTab::UvVis, SpecTab::Nmr2d] {
            app.omni_suite.spec_tab = stab;
            run_frame(&mut app);
        }
        for ip_tab in [InoPolyTab::Organometallic, InoPolyTab::Polymer] {
            app.omni_suite.ino_poly_tab = ip_tab;
            run_frame(&mut app);
        }

        let all_dialogs = [
            ActiveDialog::OpenFile, ActiveDialog::InsertFile, ActiveDialog::InsertObject,
            ActiveDialog::SaveAs, ActiveDialog::ApplyStylesheetFile, ActiveDialog::PasteSmiles,
            ActiveDialog::PasteNameAsStructure, ActiveDialog::PasteInchi, ActiveDialog::PasteMolText,
            ActiveDialog::PasteCdxml, ActiveDialog::PasteEditableSvg,
            ActiveDialog::PasteBiopolymer(PolymerKind::Peptide), ActiveDialog::PasteHelm,
            ActiveDialog::AddToNameDict, ActiveDialog::ScaleDialog, ActiveDialog::RotateDialog,
            ActiveDialog::PageSetup, ActiveDialog::DocumentSettings, ActiveDialog::DocumentProperties,
            ActiveDialog::ObjectSettings, ActiveDialog::Preferences, ActiveDialog::ListNicknames,
            ActiveDialog::DefineNickname, ActiveDialog::AtomProperties { mol_idx: 0, atom_idx: 0 },
            ActiveDialog::BondProperties { mol_idx: 0, bond_idx: 0 }, ActiveDialog::BracketProperties,
            ActiveDialog::RLogicQuery, ActiveDialog::LineSpacing, ActiveDialog::FontSizeOther,
            ActiveDialog::CheckStructureDialog, ActiveDialog::ColorPicker, ActiveDialog::AttachData,
            ActiveDialog::Annotate, ActiveDialog::FindElement, ActiveDialog::EditCaption { cap_idx: 0 },
            ActiveDialog::NewCaption { x_tenths: 1000, y_tenths: 1000 },
            ActiveDialog::HelpTopic { title: "Test".into(), body: "Body".into() },
            ActiveDialog::Shortcuts, ActiveDialog::Licensing, ActiveDialog::About,
            ActiveDialog::ConfirmExit,
        ];
        for dlg in all_dialogs {
            app.dialog = dlg;
            run_frame(&mut app);
        }
        app.dialog = ActiveDialog::None;

        // ---- 4. Verify Edit > Copy / Paste / Clear clears EVERYTHING including EnergyProfile & GelPlate ----
        app.do_copy_action();
        app.do_paste_action();
        app.tab_mut().selected_mols.clear();
        app.do_clear_action();
        assert!(app.doc().scene.molecules.is_empty());
        assert!(app.doc().scene.energy_profiles.is_empty(), "Clear must remove energy profiles");
        assert!(app.tab().gel_plates.is_empty(), "Clear must remove gel plates");
        assert!(app.doc_mut().undo(), "Clear must undo cleanly");
        assert!(!app.doc().scene.molecules.is_empty());
    }

    #[test]
    fn test_chemistry_accuracy_kekule_layout_formula_and_retrosynthesis() {
        // 1. SMILES aromatic rings must Kekulize to alternating 1/2 bonds (never order 6)
        //    and lay out with zero steric clashes / structure warnings (including ortho-substituted Aspirin!)
        for &(smi, label) in &[
            ("Oc1ccccc1C(=O)O", "Salicylic Acid"),
            ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
            ("CC(C)Cc1ccc(C(C)C(=O)O)cc1", "Ibuprofen"),
            ("CC(=O)Nc1ccccc1", "Acetanilide"),
            ("Brc1ccc(NC(=O)C)cc1", "4-Bromoacetanilide"),
        ] {
            let mut m = smiles_to_molecule(smi).expect(label);
            assert!(
                m.bonds.iter().all(|b| b.order == 1 || b.order == 2),
                "{label} must have Kekulized 1/2 bonds, found orders: {:?}",
                m.bonds.iter().map(|b| b.order).collect::<Vec<_>>()
            );
            assert!(
                m.bonds.iter().filter(|b| b.order == 2).count() >= 3,
                "{label} benzene ring must have 3 Kekulé double bonds"
            );
            clean_up_structure_preserving_stereo(&mut m, &[], 14.40);
            let warns = chiralcraft_doc::check_structure(&m, &[]);
            assert!(
                warns.is_empty(),
                "{label} should have 0 structure/overlap warnings after cleanup, got: {warns:?}"
            );
        }

        // 2. Formula formatter must subscript chemical counts only, never ordinary numbers/decimals/temperatures
        assert_eq!(format_chemical_formula("Ac2O, H3PO4 (cat.)"), "Ac₂O, H₃PO₄ (cat.)");
        assert_eq!(format_chemical_formula("85 °C, 15 min, 92%"), "85 °C, 15 min, 92%");
        assert_eq!(
            format_chemical_formula("2-acetoxybenzoic acid (aspirin) | C9H8O4 | MW: 180.16"),
            "2-acetoxybenzoic acid (aspirin) | C₉H₈O₄ | MW: 180.16"
        );
        assert_eq!(format_chemical_formula("δ 7.32 (dd)"), "δ 7.32 (dd)");
        assert_eq!(format_chemical_formula("C10H16"), "C₁₀H₁₆");
        assert_eq!(format_chemical_formula("Ca(OH)2"), "Ca(OH)₂");

        // 3. IUPAC naming distinguishes acetanilide derivatives from bare anilines
        let dict = NameDictionary::with_builtins();
        let acet = smiles_to_molecule("CC(=O)Nc1ccccc1").unwrap();
        assert_eq!(molecule_to_name(&acet, &molecule_to_smiles(&acet), &dict), "acetanilide");
        let br_acet = smiles_to_molecule("Brc1ccc(NC(=O)C)cc1").unwrap();
        assert_eq!(molecule_to_name(&br_acet, &molecule_to_smiles(&br_acet), &dict), "4-bromoacetanilide");

        // 4. Retrosynthetic Amide disconnection of Acetanilide yields non-overlapping Acetic Acid (CH3COOH) + Aniline
        let mut acet_2d = acet.clone();
        clean_up_structure_preserving_stereo(&mut acet_2d, &[], 14.40);
        let hits = suggest_disconnections(&acet_2d, 14.40);
        assert!(!hits.is_empty());
        let h0 = &hits[0];
        assert!(h0.reaction_name.contains("Amide"));
        let w_a = chiralcraft_doc::check_structure(&h0.synthon_a, &[]);
        let w_b = chiralcraft_doc::check_structure(&h0.synthon_b, &[]);
        assert!(w_a.is_empty(), "Synthon A must have zero overlap warnings: {w_a:?}");
        assert!(w_b.is_empty(), "Synthon B must have zero overlap warnings: {w_b:?}");
    }

    #[test]
    fn test_chemdraw_drawing_ux_all_8_fixes() {
        let mut app = ChiralCraftApp::new();
        let bl = app.doc().settings.bond_length_pt;

        // 1. Ring tool drag (>4pt) must NEVER silently drop the ring!
        app.tool = Tool::Benzene;
        app.handle_canvas_drag_release(120.0, 120.0, 145.0, 130.0, false);
        assert_eq!(app.doc().scene.molecules.len(), 1, "Dragged Benzene ring must be placed");
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), 6);

        // Edge-fuse a second ring via drag across bond 0 (choosing side)
        let b0 = app.doc().scene.molecules[0].bonds[0].clone();
        let p1 = (app.doc().scene.molecules[0].atoms[b0.begin_atom].x, app.doc().scene.molecules[0].atoms[b0.begin_atom].y);
        let p2 = (app.doc().scene.molecules[0].atoms[b0.end_atom].x, app.doc().scene.molecules[0].atoms[b0.end_atom].y);
        let mx = 0.5 * (p1.0 + p2.0);
        let my = 0.5 * (p1.1 + p2.1);
        app.handle_canvas_drag_release(mx, my, mx + 25.0, my - 20.0, false);
        assert_eq!(app.doc().scene.molecules.len(), 1, "Fused ring stays in same molecule");
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), 10, "Naphthalene has 10 carbons");

        // 2. SolidBond click-chain alternates zig-zag (does NOT curl into a triangle!) & quaternary carbon has zero overlap!
        app.doc_mut().apply_scene_tx("Clear", |s, _| { s.molecules.clear(); });
        app.tool = Tool::SolidBond;
        app.handle_canvas_click(100.0, 100.0, false); // 2 atoms
        for _ in 0..3 {
            let last = app.doc().scene.molecules[0].atoms.last().unwrap().clone();
            app.handle_canvas_click(last.x, last.y, false);
        }
        let m_chain = &app.doc().scene.molecules[0];
        assert_eq!(m_chain.atoms.len(), 5, "Pentane chain should have 5 carbons");
        // Verify chain extends horizontally (zig-zag), NOT curled back to atom 0!
        let end_dist = (m_chain.atoms[4].x - m_chain.atoms[0].x).hypot(m_chain.atoms[4].y - m_chain.atoms[0].y);
        assert!(end_dist > bl * 2.5, "Zig-zag pentane must extend >2.5 bond lengths, got {end_dist:.1}");

        // Click central atom twice more to make quaternary carbon -> no two neighbors overlap!
        let c2 = (m_chain.atoms[2].x, m_chain.atoms[2].y);
        app.handle_canvas_click(c2.0, c2.1, false);
        app.handle_canvas_click(c2.0, c2.1, false);
        let m_quat = &app.doc().scene.molecules[0];
        let mut min_pair = f64::INFINITY;
        for i in 0..m_quat.atoms.len() {
            for j in (i + 1)..m_quat.atoms.len() {
                min_pair = min_pair.min((m_quat.atoms[i].x - m_quat.atoms[j].x).hypot(m_quat.atoms[i].y - m_quat.atoms[j].y));
            }
        }
        assert!(min_pair > bl * 0.65, "Quaternary substituents must not overlap (min dist {min_pair:.2})");

        // 3. Closing a ring bond-by-bond via snapped endpoint auto-fuses into 6 atoms / 6 bonds (not 7 atoms!)
        app.doc_mut().apply_scene_tx("Clear", |s, _| { s.molecules.clear(); });
        // Build 5 edges of a regular hexagon, then drag the 6th bond from atom 5 toward atom 0 (slightly overshooting so raw cursor != atom 0, but snapped endpoint lands on atom 0!)
        let mut open_hex = ChiralCraftApp::build_polygon_ring(6, false, false, 200.0, 200.0, bl);
        open_hex.bonds.pop(); // remove closing bond (5 -> 0)
        app.doc_mut().stamp_molecule_tx(open_hex, "Open Hex");
        let a5 = (app.doc().scene.molecules[0].atoms[5].x, app.doc().scene.molecules[0].atoms[5].y);
        let a0 = (app.doc().scene.molecules[0].atoms[0].x, app.doc().scene.molecules[0].atoms[0].y);
        // Drag past a0 by 1.8x bond length along the a5->a0 direction so raw cursor misses a0,
        // while fixed-length snapped endpoint lands dead-on a0!
        let over_x = a5.0 + (a0.0 - a5.0) * 1.85;
        let over_y = a5.1 + (a0.1 - a5.1) * 1.85;
        app.tool = Tool::SolidBond;
        app.drag_start_atom = Some((0, 5));
        app.handle_canvas_drag_release(a5.0, a5.1, over_x, over_y, false);
        app.drag_start_atom = None;
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), 6, "Closing 6th bond must auto-fuse into 6 atoms (not 7!)");
        assert_eq!(app.doc().scene.molecules[0].bonds.len(), 6, "Ring must now have 6 bonds");

        // 4. AcyclicChain dragged off an existing atom attaches directly into that molecule!
        let before_atoms = app.doc().scene.molecules[0].atoms.len();
        app.tool = Tool::AcyclicChain;
        app.drag_start_atom = Some((0, 0));
        app.handle_canvas_drag_release(a0.0, a0.1, a0.0 + bl * 3.0, a0.1, false);
        app.drag_start_atom = None;
        assert_eq!(app.doc().scene.molecules.len(), 1, "Chain attached to atom must stay 1 molecule");
        assert!(app.doc().scene.molecules[0].atoms.len() > before_atoms + 1);

        // 5. Clicking a bond repeatedly with Wedge / Hash flips direction (Begin <-> End)!
        let wb = app.doc().scene.molecules[0].bonds[0].clone();
        let wm_x = 0.5 * (app.doc().scene.molecules[0].atoms[wb.begin_atom].x + app.doc().scene.molecules[0].atoms[wb.end_atom].x);
        let wm_y = 0.5 * (app.doc().scene.molecules[0].atoms[wb.begin_atom].y + app.doc().scene.molecules[0].atoms[wb.end_atom].y);
        app.tool = Tool::WedgeBond;
        app.handle_canvas_click(wm_x, wm_y, false);
        assert_eq!(app.doc().scene.molecules[0].bonds[0].display, BondDisplay::WedgeBegin);
        app.handle_canvas_click(wm_x, wm_y, false);
        assert_eq!(app.doc().scene.molecules[0].bonds[0].display, BondDisplay::WedgeEnd, "Second click must flip wedge direction");
        app.tool = Tool::HashWedgeBond;
        app.handle_canvas_click(wm_x, wm_y, false);
        assert_eq!(app.doc().scene.molecules[0].bonds[0].display, BondDisplay::WedgedHashBegin);
        app.handle_canvas_click(wm_x, wm_y, false);
        assert_eq!(app.doc().scene.molecules[0].bonds[0].display, BondDisplay::WedgedHashEnd, "Second click must flip hash direction");

        // 6. Single-atom drag moves ONLY that atom vertex!
        let old_a1 = (app.doc().scene.molecules[0].atoms[1].x, app.doc().scene.molecules[0].atoms[1].y);
        let old_a2 = (app.doc().scene.molecules[0].atoms[2].x, app.doc().scene.molecules[0].atoms[2].y);
        app.tool = Tool::Marquee;
        app.dragging_single_atom = Some((0, 1));
        app.handle_canvas_drag_release(old_a1.0, old_a1.1, old_a1.0 + 10.0, old_a1.1 + 8.0, false);
        app.dragging_single_atom = None;
        assert!((app.doc().scene.molecules[0].atoms[1].x - (old_a1.0 + 10.0)).abs() < 1e-6);
        assert!((app.doc().scene.molecules[0].atoms[2].x - old_a2.0).abs() < 1e-6, "Other atoms must remain fixed");

        // 7. Eraser drag-to-erase sweeps out atoms in one Undo step!
        let atoms_before_erase = app.doc().scene.molecules[0].atoms.len();
        let last_a = (app.doc().scene.molecules[0].atoms.last().unwrap().x, app.doc().scene.molecules[0].atoms.last().unwrap().y);
        app.tool = Tool::Eraser;
        app.eraser_pts = vec![(last_a.0 - 8.0, last_a.1), (last_a.0, last_a.1), (last_a.0 + 8.0, last_a.1)];
        app.handle_canvas_drag_release(last_a.0 - 8.0, last_a.1, last_a.0 + 8.0, last_a.1, false);
        app.eraser_pts.clear();
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), atoms_before_erase - 1);
        assert!(app.doc_mut().undo());
        assert_eq!(app.doc().scene.molecules[0].atoms.len(), atoms_before_erase);

        // 8. Pixel-scaled hit tolerance adapts cleanly to low vs high zoom!
        app.zoom = 0.75;
        let low_z_tol = app.atom_hit_tol();
        app.zoom = 2.0;
        let hi_z_tol = app.atom_hit_tol();
        assert!(low_z_tol > hi_z_tol * 2.0, "Low zoom must expand doc-space hit radius ({low_z_tol} vs {hi_z_tol})");
    }

    #[test]
    fn test_retrosynthesis_example_loader_and_non_overlapping_nmr_labels() {
        // 1. Default new document starts maximized (clean full-bleed white canvas)
        let mut app = ChiralCraftApp::new();
        assert!(app.tab().maximized, "New document tabs must start maximized");

        // 2. Retrosynthesis 1-click example loader stamps target & auto-populates disconnections
        app.handle_pro_suite_action(ProSuiteAction::LoadExampleTarget("Acetanilide", "CC(=O)Nc1ccccc1"));
        assert_eq!(app.doc().scene.molecules.len(), 1);
        assert!(!app.pro_suite.retro_hits.is_empty(), "Acetanilide must yield amide disconnection");
        assert!(app.pro_suite.highlighted_retro_bond.is_some());

        // 3. Stamped 1H-NMR δ badges on Aspirin must be non-overlapping and clear of atoms
        app.handle_pro_suite_action(ProSuiteAction::LoadExampleTarget("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"));
        app.stamp_nmr_shift_labels(false);
        let badges: Vec<&CaptionNode> = app
            .doc()
            .scene
            .captions
            .iter()
            .filter(|c| c.text.starts_with('δ'))
            .collect();
        assert!(badges.len() >= 4, "Aspirin must stamp multiple 1H-NMR δ labels");
        for i in 0..badges.len() {
            for j in (i + 1)..badges.len() {
                let dy = (badges[i].y - badges[j].y).abs();
                let dx = (badges[i].x - badges[j].x).abs();
                assert!(
                    dy > 9.0 || dx > 28.0,
                    "NMR badges '{}' and '{}' overlapped at ({:.1},{:.1}) vs ({:.1},{:.1})",
                    badges[i].text,
                    badges[j].text,
                    badges[i].x,
                    badges[i].y,
                    badges[j].x,
                    badges[j].y
                );
            }
        }
    }

    #[test]
    fn test_comprehensive_audit_app_regressions() {
        // 1. E/Z SMILES with slashes must load as SMILES, NOT be misclassified as a file path!
        let ez_doc = load_document_or_smiles("C/C=C/C").expect("E-2-butene SMILES must load");
        assert_eq!(ez_doc.scene.molecules.len(), 1);
        assert_eq!(ez_doc.scene.molecules[0].atoms.len(), 4);

        // 2. Atom hotkeys: 'a'/'3' sprouts Benzene with exactly 3 double bonds (not 4!); 'A' sets Acetyl ("Ac", Z=6), NOT Actinium (Z=89)!
        let mut app = ChiralCraftApp::new();
        let mut seed = ChemMolecule::default();
        seed.add_atom(ChemAtom::new(1, 6));
        app.doc_mut().stamp_molecule_tx(seed, "Seed");
        assert!(matches!(apply_atom_hotkey(app.doc_mut(), 0, 0, 'a'), HotkeyOutcome::Handled(_)));
        let m = &app.doc().scene.molecules[0];
        assert_eq!(m.atoms.len(), 6);
        assert_eq!(m.bonds.iter().filter(|b| b.order == 2).count(), 3, "Benzene hotkey must produce 3 double bonds, not 4");
        assert!(matches!(apply_atom_hotkey(app.doc_mut(), 0, 0, 'A'), HotkeyOutcome::Handled(_)));
        let a0 = &app.doc().scene.molecules[0].atoms[0];
        assert_eq!(a0.label.as_deref(), Some("Ac"));
        assert_eq!(a0.z, 6, "'A' hotkey must set Acetyl nickname (C), NOT Actinium (Z=89)");

        // 3. JSON parser depth limit, trailing garbage rejection, and multi-byte UTF-8 handling
        let deep = format!("{}{}", "[".repeat(200), "]".repeat(200));
        assert!(crate::chem_db::JsonValue::parse(&deep).is_err(), "Deeply nested JSON must fail gracefully");
        assert!(crate::chem_db::JsonValue::parse(r#"{"a":1} trailing"#).is_err(), "Trailing JSON tokens must be rejected");
        let utf8_v = crate::chem_db::JsonValue::parse(r#"{"name":"α-D-glucose (ΔG‡)"}"#).unwrap();
        assert_eq!(utf8_v.get("name").and_then(|v| v.as_str()), Some("α-D-glucose (ΔG‡)"));

        // 4. MCP JSON-RPC 2.0 compliance: notifications return None, non-object returns -32600, unknown tool returns -32602
        assert!(mcp::handle_mcp_request(r#"{"jsonrpc":"2.0","method":"notifications/progress","params":{}}"#).is_none());
        let bad_req = mcp::handle_mcp_request("[1,2,3]").unwrap();
        assert!(bad_req.contains("-32600"), "Non-object JSON-RPC must return -32600: {bad_req}");
        let bad_tool = mcp::handle_mcp_request(
            r#"{"jsonrpc":"2.0","id":7,"method":"tools/call","params":{"name":"nonexistent_tool","arguments":{}}}"#,
        ).unwrap();
        assert!(bad_tool.contains("-32602"), "Unknown tool must return -32602: {bad_tool}");
        let bad_nuc = mcp::handle_mcp_request(
            r#"{"jsonrpc":"2.0","id":8,"method":"tools/call","params":{"name":"predict_nmr","arguments":{"smiles":"CCO","nucleus":"19F"}}}"#,
        ).unwrap();
        assert!(bad_nuc.contains("\"isError\":true"), "Unsupported NMR nucleus must return isError:true");

        // 5. --edit-svg on an existing corrupt file must NOT bind roundtrip_svg_path or file_path (protects file from overwrite!)
        let corrupt = std::env::temp_dir().join(format!("chiralcraft-corrupt-{}.svg", std::process::id()));
        fs::write(&corrupt, "not a valid chiralcraft svg").unwrap();
        let app_bad = ChiralCraftApp::new_with_file(&corrupt.to_string_lossy(), true);
        assert!(app_bad.roundtrip_svg_path.is_none(), "Failed load on existing file must not bind roundtrip_svg_path");
        assert!(app_bad.tabs[0].file_path.is_none(), "Failed load on existing file must not bind file_path");
        let _ = fs::remove_file(&corrupt);

        // 6. Quantum runner: >128KB stdout never deadlocks pipe buffer & flag-injecting solvent is rejected!
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let tmp = std::env::temp_dir().join(format!("chiralcraft-bigout-{}", std::process::id()));
            let _ = fs::create_dir_all(&tmp);
            let mock = tmp.join("big_xtb.sh");
            let script = r#"#!/bin/sh
i=0
while [ $i -lt 3000 ]; do
  echo "SCF iteration line $i padding 01234567890123456789012345678901234567890123456789"
  i=$((i + 1))
done
cat << 'EOF' > g98.out
 Gaussian 98: xtb output
 Standard orientation:
 ---------------------------------------------------------------------
    1          8           0        0.000000    0.000000    0.117300
    2          1           0        0.000000    0.757200   -0.469200
    3          1           0        0.000000   -0.757200   -0.469200
 ---------------------------------------------------------------------
                      1
                     A1
 Frequencies --   1595.10
 Red. masses --      1.08
 IR Inten    --     68.50
 Atom AN      X      Y      Z
   1   8     0.00   0.00   0.07
   2   1     0.00   0.43  -0.55
   3   1     0.00  -0.43  -0.55
 Normal termination of Gaussian 98.
EOF
echo "normal termination of xtb"
"#;
            fs::write(&mock, script).unwrap();
            fs::set_permissions(&mock, fs::Permissions::from_mode(0o755)).unwrap();
            let water = smiles_to_molecule("O").unwrap();
            let cancel = std::sync::atomic::AtomicBool::new(false);
            let out = crate::quantum_runner::run_qc_job_blocking(
                QcRunBackend::XtbGfn2,
                &water,
                "O|big".into(),
                &QuantumJobSpec::default(),
                &mock.to_string_lossy(),
                &cancel,
            ).expect("Large (>128KB) stdout xTB job must complete without pipe deadlock");
            assert_eq!(out.ir_prediction.bands.len(), 1);

            let mut bad_spec = QuantumJobSpec::default();
            bad_spec.solvent = "--evil-flag".into();
            assert!(
                crate::quantum_runner::run_qc_job_blocking(
                    QcRunBackend::XtbGfn2,
                    &water,
                    "O|bad".into(),
                    &bad_spec,
                    &mock.to_string_lossy(),
                    &cancel,
                ).is_err(),
                "Flag-injecting solvent must be rejected"
            );
            let _ = fs::remove_dir_all(&tmp);
        }
    }
}
